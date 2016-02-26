/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010-2016                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
/*;  This program is free software; you can redistribute it and/or    */
/*;  modify it under the terms of the GNU General Public License as   */
/*;  published by the Free Software Foundation; either version 2 of   */
/*;  the License, or (at your option) any later version.              */
/*;                                                                   */
/*;  This program is distributed in the hope that it will be useful,  */
/*;  but WITHOUT ANY WARRANTY; without even the implied warranty of   */
/*;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    */
/*;  GNU General Public License for more details.                     */
/*;                                                                   */
/*;  You should have received a copy of the GNU General Public        */
/*;  License along with this program; if not, write to the Free       */
/*;  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,     */
/*;  MA 02139, USA.                                                   */
/*;                                                                   */
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitFArray.h"
#include "ObitFArrayUtil.h"
#include "ObitImageDesc.h"
#include "ObitImageSel.h"
#include "ObitImageMF.h"
#include "ObitFFT.h"
#include "ObitMem.h"
#include "ObitSystem.h"
#include "ObitBeamShape.h"
#include "ObitSpectrumFit.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageMF.c
 * ObitImageMF class function definitions.
 * This class is derived from the ObitData base class.
 *
 * This class contains an astronomical image and allows access.
 * An ObitImageMF is the front end to a persistent disk resident structure.
 *
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitImageMF";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitImageGetClass;

/** Flag to limit frequency division messages */
static gboolean doFreqDivMess = TRUE;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitImageMFClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitImageMFClassInfo myClassInfo = {FALSE};

/*---------------Private structures----------------*/
/* Spectrum fitting threaded function argument */
typedef struct {
  /* Number of spectral channels  */
  olong   nSpec;
  /* Array (nSpec) channel frequencies (hz) */
  odouble *Freq;
  /* Spectral index correction previousls applied to data */
  ofloat alpha;
  /* Array of sigma for each channel */
  ofloat *sigma;
  /** BeamShape object */
  ObitBeamShape *BeamShape;
  /* Array  (nSpec) of rows of pixel values in image in each of nSpec planes */
  ofloat **inData;
  /* Work array of flux for each channel */
  ofloat *workFlux;
  /* Work array of sigmarfor each channel */
  ofloat *workSigma;
  /* 0-rel row number */
  olong iy;
  /* Input image descriptor */
  ObitImageDesc *desc;
  /* Correct for primary Beam? */
  gboolean doPBCorr;
  /* Antenna diameter */
  ofloat antSize;
  /* Max. number of orders to fit  */
  olong  nOrder;
  /* (nOrder+1) Arrays of fitted values */
  ofloat **outData;
  /* Fitting argument */
  gpointer fitArg;
  /* Fitting result array */
  ofloat *fitResult;
  /* thread number , >0 -> no threading  */
  olong   ithread;
  /* Obit Thread object */
  ObitThread  *thread;
  /* Obit error stack object */
  ObitErr  *err;
  /* Fitting argument */
} FitSpecFuncArg;

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitImageMFInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitImageMFClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitImageMFClassInfoDefFn (gpointer inClass);

/** Private: Make arguments for Threaded CLEAN */
static olong MakeFitSpecArgs (ObitImageMF *in, olong maxThread,
			      ofloat antSize, olong nOrder,
			      FitSpecFuncArg ***args, 
			      ObitErr *err);

/** Private: Delete arguments for Threaded CLEAN */
static void KillFitSpecArgs (olong nargs, FitSpecFuncArg **args);

/** Private: Threaded Spectrum fitting   */
static gpointer ThreadFitSpec (gpointer arg);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitImageMF* newObitImageMF (gchar* name)
{
  ObitImageMF* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageMFClassInit();

  /* allocate/init structure */
  out = ObitMemAlloc0Name(sizeof(ObitImageMF), "ObitImageMF");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitImageMFInit((gpointer)out);

 return out;
} /* end newObitImageMF */

/**
 * Create a scratch file suitable for accepting the data to be read from in.
 * A scratch ImageMF is more or less the same as a normal ImageMF except that it is
 * automatically deleted on the final unreference.
 * The output will have the underlying files of the same type as in already 
 * allocated.
 * \param in  The object to copy (ObitImageMF), info may have 
 * \li ScrSize OBIT_int (?,1,1) Dimension of the desired scratch Image
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitImage* newObitImageMFScratch (ObitImage *in, ObitErr *err)
{
  /*const ObitClassInfo *ParentClass;*/
  ObitImage   *out=NULL;
  ObitImageMF *inMF=NULL, *outMF=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, size[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *outName;
  gchar *routine = "newObitImageMFScratch";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Make sure input ObitImageMF */
  Obit_retval_if_fail((ObitImageMFIsA(in)), err, NULL,
  		      "%s: Input image, %s NOT ImageMF", 
		      routine, in->name);

  /* Ensure in fully instantiated -assume OK if myIO exists */
  if (!in->myIO) ObitImageFullInstantiate (in, TRUE, err);
  if (err->error)Obit_traceback_val (err, routine, in->name, out);

  /* Create - derive object name */
  outName = g_strconcat ("Scratch Copy: ",in->name,NULL);
  out = (ObitImage*)newObitImageMF(outName);
  g_free(outName);

  /* Mark as scratch */
  out->isScratch = TRUE;

  /* deep copy any base class members - NO fails 
     ParentClass = myClassInfo.ParentClass;
     g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
     ParentClass->ObitCopy (in, out, err);*/

  /* Copy descriptor */
  out->myDesc = (gpointer)ObitImageDescCopy(in->myDesc, out->myDesc, err);

  /* Check if different size needed */
  if (ObitInfoListGetTest(in->info, "ScrSize", &type, dim, size)) {
    for (i=0; i<MIN (dim[0], IM_MAXDIM); i++) 
      if (size[i]>0) out->myDesc->inaxes[i] = size[i];
  }
 
  /* Force to float pixels */
  out->myDesc->bitpix=-32;

  /* Allocate underlying file */
  ObitSystemGetScratch (in->mySel->FileType, "MA", out->info, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  
  /* Register in the scratch file list */
  ObitSystemAddScratch ((Obit*)out, err);
  
  /* same size IO as input */
  dim[0] = 1;
  ObitInfoListPut (out->info, "IOBy", OBIT_long, dim, &in->myDesc->IOsize, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /* Copy MF stuff */
  inMF  = (ObitImageMF*)in;
  outMF = (ObitImageMF*)out;
  outMF->maxOrder = inMF->maxOrder;
  outMF->curOrder = inMF->curOrder;
  outMF->refFreq  = inMF->refFreq;
  outMF->fresh = inMF->fresh;
  outMF->nSpec = inMF->nSpec;
  if (outMF->BIFSpec)   g_free(outMF->BIFSpec);
  if (outMF->EIFSpec)   g_free(outMF->EIFSpec);
  if (outMF->BChanSpec) g_free(outMF->BChanSpec);
  if (outMF->EChanSpec) g_free(outMF->EChanSpec);
  if (outMF->specFreq)  g_free(outMF->specFreq);
  if (outMF->specFreqLo)  g_free(outMF->specFreqLo);
  if (outMF->specFreqHi)  g_free(outMF->specFreqHi);
  outMF->BIFSpec   = g_malloc0(inMF->nSpec*sizeof(olong));
  outMF->EIFSpec   = g_malloc0(inMF->nSpec*sizeof(olong));
  outMF->BChanSpec = g_malloc0(inMF->nSpec*sizeof(olong));
  outMF->EChanSpec = g_malloc0(inMF->nSpec*sizeof(olong));
  outMF->specFreq  = g_malloc0(inMF->nSpec*sizeof(odouble));
  outMF->specFreqLo= g_malloc0(inMF->nSpec*sizeof(odouble));
  outMF->specFreqHi= g_malloc0(inMF->nSpec*sizeof(odouble));
  for (i=0; i<inMF->nSpec; i++) {
    outMF->BIFSpec[i]   = inMF->BIFSpec[i];
    outMF->EIFSpec[i]   = inMF->EIFSpec[i];
    outMF->BChanSpec[i] = inMF->BChanSpec[i];
    outMF->EChanSpec[i] = inMF->EChanSpec[i];
    outMF->specFreq[i]  = inMF->specFreq[i];
    outMF->specFreqLo[i]= inMF->specFreqLo[i];
    outMF->specFreqHi[i]= inMF->specFreqHi[i];
  }

  /* Fully instantiate output */
  ObitImageFullInstantiate (out, FALSE, err);
  if (err->error)Obit_traceback_val (err, routine, out->name, out);
 
  return out;
} /* end newObitImageMFScratch */

/**
 * Constructor from ObitImage.
 * output will have same underlying file definition as in.
 * Adds planes as needed for spectral planes and 
 * changes frequency axis if norder >0.
 * Initializes class if needed on first call.
 * \param in     ObitImage to copy, MUST NOT have been fully instantiated
 * \param inData ObitUV from which data to be imaged.
 * \param norder Maximum order number of beam, 
 *        0=Normal, 1=spectral index, 2=also curvature
 * \param maxFBW Maximum fractional bandwidth at center of each IF
 * \param err    ObitErr for reporting errors.
 * \return the new object.
 */
ObitImageMF* ObitImageMFFromImage (ObitImage* in, ObitUV *inData, 
				   olong norder, ofloat maxFBW, 
				   ofloat alpha, odouble alphaRefF,
				   ObitErr *err)
{
  ObitImageMF* out = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "ObitImageMFFromImage";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitImageIsA(in));

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageMFClassInit();

  /* allocate/init structure */
  out = ObitMemAlloc0Name(sizeof(ObitImageMF), "ObitImageMF");

  /* initialize values */
  if (in->name!=NULL) out->name = g_strdup(in->name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitImageMFInit((gpointer)out);

  /* Copy info for same definition */
  out->info = ObitInfoListCopyData (in->info, out->info);

  /* Copy descriptor */
  out->myDesc = ObitImageDescCopy(in->myDesc, out->myDesc, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /* Actually adding spectral information? */
  out->myDesc->inaxes[out->myDesc->jlocf] = 1;
  if (norder>0) {
    /* Change frequency axis labeling */ 
    strncpy (out->myDesc->ctype[out->myDesc->jlocf], "SPECLNMF", IMLEN_KEYWORD-1);
    /* Save reference frequency */
    out->refFreq = out->myDesc->crval[out->myDesc->jlocf];  
    
    /* Add order to infolist */
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut(out->info, "NTERM", OBIT_long, dim, &norder);
  }  /* end Actually adding spectral information */

  /* Beams and things */
  out->myBeam = ObitRef(in->myBeam);
  ObitImageMFSetOrder (out, norder, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /* Set coarse channelization */
  ObitImageMFSetSpec (out, inData, maxFBW, alpha, alphaRefF, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  return out;
} /* end ObitImageMFFromImage */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitImageMFGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObiImagetGetWBClass */

/**
 * Test if two ObitImageMFs have the same underlying structures.
 * This test is done using values entered into the #ObitInfoList
 * in case the object has not yet been opened.
 * \param in1 First object to compare
 * \param in2 Second object to compare
 * \param err ObitErr for reporting errors.
 * \return TRUE if to objects have the same underlying structures
 * else FALSE
 */
gboolean ObitImageMFSame (ObitImageMF *in1, ObitImageMF *in2, ObitErr *err )
{
  /* Call ObitData function */
  return ObitDataSame ((ObitData*)in1, (ObitData*)in2, err);
} /* end ObitImageMFSame */

/**
 * Delete underlying files and the basic object.
 * \param inn  Pointer to object to be zapped.
 * \param err  ObitErr for reporting errors.
 * \return pointer for input object, NULL if deletion successful
 */
ObitImage* ObitImageMFZap (ObitImage *inn, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitImage *in = (ObitImage*)inn;
  gchar *routine = "ObitImageMFZap";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return inn;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* Any actual I/O? If not just delete object */
  if (in->mySel->FileType==OBIT_IO_MEM) return ObitImageMFUnref(in);

  /* Close if still active */
  if ((in->myStatus == OBIT_Active) || (in->myStatus == OBIT_Modified)){
   retCode = ObitIOClose(in->myIO, err);
   if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
     Obit_traceback_val (err, routine, in->name, inn);    
  }

  /* Ensure in fully instantiated */
  ObitErrLog(err); /* Show any pending messages as they may get lost */
  ObitImageFullInstantiate (inn, TRUE, err);
  /* If this fails, clear errors and assume it doesn't exist */
  if (err->error) { 
    ObitErrClearErr(err); 
    return ObitImageMFUnref(in); 
  }

  /* Delete Image and all tables  */
  ObitIOZap (in->myIO, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, inn);

  /* If it's scratch remove from list */
  if (in->isScratch) ObitSystemFreeScratch ((Obit*)in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, in);

  /* Delete object */
  in->isScratch = 0; /* Already deleted underlying structures */
  while (in) in = ObitImageMFUnref(in);
  
  return (ObitImage*)in;
} /* end ObitImageMFZap */

/**
 * Make a deep copy of input object.
 * Copies are made of complex members including disk files; these 
 * will be copied applying whatever selection is associated with the input.
 * Objects should be closed on input and will be closed on output.
 * In order for the disk file structures to be copied, the output file
 * must be sufficiently defined that it can be written.
 * The copy will be attempted but no errors will be logged until
 * both input and output have been successfully opened.
 * ObitInfoList and ObitThread members are only copied if the output object
 * didn't previously exist.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitImageMF* ObitImageMFCopy (ObitImageMF *in, ObitImageMF *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOCode iretCode, oretCode;
  gboolean oldExist;
  olong i;
  ObitHistory *inHist=NULL, *outHist=NULL;
  gchar *outName=NULL, *today=NULL;
  gchar *routine = "ObitImageMFCopy";

  /* error checks */
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Any actual I/O? Can't handle Memory only */
  if (in->mySel->FileType==OBIT_IO_MEM) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Image %s in memory only and I do not know how to copy it", 
		   routine, in->name);
    return NULL;
  }

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitImageMF(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes other additions only if out newly created */
  if (!oldExist) {
    /* copy */
    out->myDesc = ObitImageDescCopy(in->myDesc, out->myDesc, err);
    /* Don't copy selector */
    if (out->mySel) out->mySel = ObitUnref (out->mySel);
    out->mySel = newObitImageSel (out->name);
    /* Don't copy info */
    /*out->info = ObitInfoListUnref(out->info); */
    /*out->info = ObitInfoListRef(in->info); */
    /* Output will initially have no associated tables */
    out->tableList = ObitTableListUnref(out->tableList);
    out->tableList = newObitTableList(out->name);
    /* don't copy ObitThread  */
}

  /* If the output object was created this call it cannot be fully
     defined so we're done */
  if (!oldExist) return out;

  /* if input has file designated, copy data */
  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitImageOpen ((ObitImage*)in, OBIT_IO_ReadOnly, err);
  /* if it didn't work bail out */
  if ((iretCode!=OBIT_IO_OK) || (err->error)) {
    return out;
  }

  /* copy Descriptor - this time with full information */
  out->myDesc = ObitImageDescCopy(in->myDesc, out->myDesc, err);
  /* Float it */
  out->myDesc->bitpix = -32;

  /* Copy MF stuff */
  out->maxOrder = in->maxOrder;
  out->curOrder = in->curOrder;
  out->refFreq  = in->refFreq;
  out->fresh = in->fresh;
  out->nSpec = in->nSpec;
  if (out->BIFSpec)   g_free(out->BIFSpec);
  if (out->EIFSpec)   g_free(out->EIFSpec);
  if (out->BChanSpec) g_free(out->BChanSpec);
  if (out->EChanSpec) g_free(out->EChanSpec);
  if (out->specFreq)  g_free(out->specFreq);
  if (out->specFreqLo)  g_free(out->specFreqLo);
  if (out->specFreqHi)  g_free(out->specFreqHi);
  out->BIFSpec   = g_malloc0(in->nSpec*sizeof(olong));
  out->EIFSpec   = g_malloc0(in->nSpec*sizeof(olong));
  out->BChanSpec = g_malloc0(in->nSpec*sizeof(olong));
  out->EChanSpec = g_malloc0(in->nSpec*sizeof(olong));
  out->specFreq   = g_malloc0(in->nSpec*sizeof(odouble));
  out->specFreqLo = g_malloc0(in->nSpec*sizeof(odouble));
  out->specFreqHi = g_malloc0(in->nSpec*sizeof(odouble));
  for (i=0; i<in->nSpec; i++) {
    out->BIFSpec[i]   = in->BIFSpec[i];
    out->EIFSpec[i]   = in->EIFSpec[i];
    out->BChanSpec[i] = in->BChanSpec[i];
    out->EChanSpec[i] = in->EChanSpec[i];
    out->specFreq[i]  = in->specFreq[i];
    out->specFreqLo[i]= in->specFreqLo[i];
    out->specFreqHi[i]= in->specFreqHi[i];
  }

  /* Creation date today */
  today = ObitToday();
  strncpy (out->myDesc->date, today, IMLEN_VALUE);
  if (today) g_free(today);
 
  /* use same data buffer on input and output 
     so don't assign buffer for output */
  out->extBuffer = TRUE;

  /* test open output */
  ObitErrLog(err); /* Show any pending messages as they may get lost */
  oretCode = ObitImageOpen ((ObitImage*)out, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitImageOpen ((ObitImage*)out, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset external buffer */
    out->extBuffer = FALSE;
    return out;
  }

  /* Copy any history  unless Scratch */
  if (!in->isScratch && !out->isScratch) {
    inHist  = newObitDataHistory((ObitData*)in, OBIT_IO_ReadOnly, err);
    outHist = newObitDataHistory((ObitData*)out, OBIT_IO_WriteOnly, err);
    outHist = ObitHistoryCopy (inHist, outHist, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
    inHist  = ObitHistoryUnref(inHist);
    outHist = ObitHistoryUnref(outHist);
  }

  /* make sure the access sizes are the same */
  if (in->myDesc->IOsize != out->myDesc->IOsize) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Access sizes of two images differ", routine);
    Obit_log_error(err, OBIT_Error, 
		   "Objects %s %s", in->name, out->name);
    /* unset external buffer */
    out->extBuffer = FALSE;
    return out;
  }
  
  /* we're in business, copy */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    iretCode = ObitImageRead ((ObitImage*)in, in->image->array, err);
    if (iretCode!=OBIT_IO_OK) break;
    oretCode = ObitImageWrite ((ObitImage*)out, in->image->array, err);
  }
  
  /* unset external buffer */
  out->extBuffer = FALSE;
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, out);
  
  /* close files to be sure */
  iretCode = ObitImageClose ((ObitImage*)in, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, out);
  
  /* close files to be sure */
  oretCode = ObitImageClose ((ObitImage*)out, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, out->name, out);
  
  return out;
} /* end ObitImageMFCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an Image similar to the input one.
 * Output version set to floating pixels
 * \param in  The object to copy, info may have
 * \param out An existing object pointer for output, info may have
 * \li Size OBIT_int (?,1,1) Dimension of the desired  Image
 * \param err Error stack, returns if not empty.
 */
void ObitImageMFClone  (ObitImageMF *in, ObitImageMF *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOCode iretCode, oretCode;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, size[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ObitHistory *inHist=NULL, *outHist=NULL;
  gchar *today=NULL;
  gchar *exclude[]={"AIPS CC", "AIPS HI", "AIPS PL", "AIPS SL", NULL}; 
  gchar *routine = "ObitImageMFClone";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* Any actual I/O? Can't handle Memory only */
  if (in->mySel->FileType==OBIT_IO_MEM) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Image %s in memory only and I do not know how to clone it", 
		   routine, in->name);
    return;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes other additions */
  /* Don't copy selector */
  if (out->mySel) out->mySel = ObitUnref (out->mySel);
  out->mySel = newObitImageSel (out->name);
  /* Output will initially have no associated tables */
  out->tableList = ObitTableListUnref(out->tableList);
  out->tableList = newObitTableList(out->name);
  /* don't copy ObitThread  */

  /* Open to fully instantiate input and see if it's OK */
  iretCode = ObitImageOpen ((ObitImage*)in, OBIT_IO_ReadOnly, err);
  if ((iretCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, in->name);

  /* copy Descriptor */
  out->myDesc = ObitImageDescCopy(in->myDesc, out->myDesc, err);

  /* Check if different size needed */
  if (ObitInfoListGetTest(out->info, "Size", &type, dim, size)) {
    for (i=0; i<MIN (dim[0], IM_MAXDIM); i++) 
      if (size[i]>0) out->myDesc->inaxes[i] = size[i];
  }
 
  /* Copy MF stuff */
  out->maxOrder = in->maxOrder;
  out->curOrder = in->curOrder;
  out->refFreq  = in->refFreq;
  out->fresh = in->fresh;
  out->nSpec = in->nSpec;
  if (out->BIFSpec)   g_free(out->BIFSpec);
  if (out->EIFSpec)   g_free(out->EIFSpec);
  if (out->BChanSpec) g_free(out->BChanSpec);
  if (out->EChanSpec) g_free(out->EChanSpec);
  if (out->specFreq)  g_free(out->specFreq);
  if (out->specFreqLo)  g_free(out->specFreqLo);
  if (out->specFreqHi)  g_free(out->specFreqHi);
  out->BIFSpec   = g_malloc0(in->nSpec*sizeof(olong));
  out->EIFSpec   = g_malloc0(in->nSpec*sizeof(olong));
  out->BChanSpec = g_malloc0(in->nSpec*sizeof(olong));
  out->EChanSpec = g_malloc0(in->nSpec*sizeof(olong));
  out->specFreq   = g_malloc0(in->nSpec*sizeof(odouble));
  out->specFreqLo = g_malloc0(in->nSpec*sizeof(odouble));
  out->specFreqHi = g_malloc0(in->nSpec*sizeof(odouble));
  for (i=0; i<in->nSpec; i++) {
    out->BIFSpec[i]   = in->BIFSpec[i];
    out->EIFSpec[i]   = in->EIFSpec[i];
    out->BChanSpec[i] = in->BChanSpec[i];
    out->EChanSpec[i] = in->EChanSpec[i];
    out->specFreq[i]  = in->specFreq[i];
    out->specFreqLo[i]= in->specFreqLo[i];
    out->specFreqHi[i]= in->specFreqHi[i];
  }

  /* Creation date today */
  today = ObitToday();
  strncpy (out->myDesc->date, today, IMLEN_VALUE);
  if (today) g_free(today);
 
  /* Force to float pixels */
  out->myDesc->bitpix=-32;

  /* Open output */
  oretCode = ObitImageOpen ((ObitImage*)out, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitImageOpen ((ObitImage*)out, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset external buffer */
    out->extBuffer = FALSE;
    Obit_traceback_msg (err, routine, out->name);
  }

  /* Copy any history unless Scratch  */
  if (!in->isScratch && !out->isScratch) {
    inHist  = newObitDataHistory((ObitData*)in, OBIT_IO_ReadOnly, err);
    outHist = newObitDataHistory((ObitData*)out, OBIT_IO_WriteOnly, err);
    outHist = ObitHistoryCopy (inHist, outHist, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    inHist  = ObitHistoryUnref(inHist);
    outHist = ObitHistoryUnref(outHist);
  }

 /* Copy tables  */
  iretCode = ObitImageCopyTables ((ObitImage*)in, (ObitImage*)out, exclude, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Close files */
  ObitImageClose ((ObitImage*)in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  ObitImageClose ((ObitImage*)out, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
} /* end ObitImageMFClone */

/**
 * Ensures full instantiation of object - basically open to read/write header
 * and verify or create file.
 * If object has previously been opened, as demonstrated by the existance
 * of its myIO member, this operation is a no-op.
 * \param in     Pointer to object
 * \param exist  TRUE if object should previously exist, else FALSE
 * \param err    ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
void ObitImageMFFullInstantiate (ObitImageMF *in, gboolean exist, ObitErr *err)
{
  ObitIOAccess access;
  gchar *routine = "ObitImageMFFullInstantiate";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  if (in->myIO) return;  /* is this needed? */

  /* Open readonly if it should exist, else writeonly */
  if (exist) access = OBIT_IO_ReadOnly;
  else access = OBIT_IO_WriteOnly;
  in->extBuffer = TRUE;  /* Don't need to assign buffer here */

  /* Open and close */
  ObitImageOpen((ObitImage*)in, access, err);
  ObitImageClose((ObitImage*)in, err);
  if (err->error)Obit_traceback_msg (err, routine, in->name);
  in->extBuffer = FALSE;  /* May need buffer later */
} /* end ObitImageMFFullInstantiate */

/**
 * Sets maximum order and initializes work arrays
 * Must be called before object fully instantiated.
 * Increments the descriptor number of planes by order
 * \param in     Pointer to object, should be fully defined
 * \param order  Max. imaging order to use, 1=>alpha, min 0.
 * \param err    ObitErr for reporting errors.
 */
void ObitImageMFSetOrder (ObitImageMF *in, olong order, 
			  ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong nterm;
  /*gchar *routine = "ObitImageMFSetOrder";*/
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  
  in->maxOrder = MAX (0, order);
  in->curOrder = in->maxOrder;

  /* Increment number of planes by order */
  in->myDesc->inaxes[in->myDesc->jlocf] += in->maxOrder;

  /* Make note to self */
  nterm = order+1;
  ObitInfoListAlwaysPut (in->myDesc->info, "NTERM", OBIT_long, dim, &nterm);
} /* end  ObitImageMFSetOrder */

/**
 * Defines the IF and frequency ranges and number of coarse channels to image
 * Must be called before object fully instantiated.
 * Calculates the number of coarse channels and sets 
 * nSpec, specFreq, BIFSpec, EIFSpec, BChanSpec and EChanSpec.
 * Each IF divided into an integran number of coarse channels.
 * Increments the descriptor number of planes by nSpec
 * Frequency information added to image descriptor with keywords
 * NSPEC, FREQ001...
 * The low and high frequencies of each bin are FREL0001...,FREH0001...,
 * The value of the alpha, alphaRefF are similarly saved in the HEADER
 * \param in        Pointer to object, should be fully defined
 * \param uvdata    UV data to be imaged
 * \param maxFBW    Maximum fractional bandwidth at center of each IF
                    If < 0 then divide data by IF.
 * \param alpha     Prior spectral index
 * \param alphaRefF Reference frequency for alpha
 * \param err       ObitErr for reporting errors.
 */
void ObitImageMFSetSpec (ObitImageMF *in, ObitUV *inData, ofloat maxFBW,
			 ofloat alpha, odouble alphaRefF, ObitErr *err)
{
  olong iif, ichan, iclo, ichi, nSpec=0, nIF, nChan, incf, incif;
  olong fincf, fincif, i, count, count2=0, iCh;
  odouble sum=0.0, sum2=0.0, mxFreq, freqLo, freqHi;
  ObitUVDesc *uvdesc;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar keyword[12];
  gboolean done, lsb;
  olong IFBreak[1001], ChBreak[1001], nBreak=1000;  /* Frequency bin breaks */
  gchar *routine = "ObitImageMFSetSpec";
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Make sure inData fully instantiated */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* frequency tables if not defined */
  if (inData->myDesc->freqArr==NULL) ObitUVGetFreq (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  uvdesc = inData->myDesc;
  nChan = uvdesc->inaxes[uvdesc->jlocf];
  incf  = uvdesc->incf;
  nIF   = 1;
  incif = uvdesc->incs;
  if (uvdesc->jlocif>=0) {
    nIF = uvdesc->inaxes[uvdesc->jlocif];
    incif = uvdesc->incif;
  }
  lsb = uvdesc->cdelt[uvdesc->jlocf]<0.0;  /* Lower sideband? */

  /* determine coarse channels, use 0-rel values  */
  nSpec = 0;
  /* Only if in->curOrder>0 */
  if (in->curOrder>0) {
    nSpec = 1;   /* Number of break points including ends */
    IFBreak[0] = 0;
    ChBreak[0] = -1;
    freqLo = (uvdesc->freqIF[0] + (1.0-uvdesc->crpix[uvdesc->jlocf])*uvdesc->cdelt[uvdesc->jlocf]);
    if (lsb) mxFreq = (1.0-maxFBW) * freqLo;  /* Lowest frequency in first bin */
    else     mxFreq = (1.0+maxFBW) * freqLo;  /* Highest frequency in first bin */
    /* Loop over IFs */
    for (iif=0; iif<nIF; iif++) {
      /* Is this IF all in current bin? */
      freqLo = (uvdesc->freqIF[iif] + (1.0-uvdesc->crpix[uvdesc->jlocf])*uvdesc->cdelt[uvdesc->jlocf]);
      freqHi = (uvdesc->freqIF[iif] + (nChan-uvdesc->crpix[uvdesc->jlocf]+1)*uvdesc->cdelt[uvdesc->jlocf]);
      if (lsb) done = freqHi>mxFreq;   /* LSB */
      else     done = freqHi<mxFreq;   /* USB */
      /* By IF? */
      if (maxFBW<0.0) {
	IFBreak[nSpec]   = iif;
	ChBreak[nSpec++] = nChan;
      } else {  /* may split IFs */
	/* A break in this IF? */
	while (!done) {
	  /* Gap in frequency? */
	  if (((freqLo>mxFreq)&&(!lsb)) || ((freqLo<mxFreq)&&(lsb))) {
	    IFBreak[nSpec]   = iif-1;
	    ChBreak[nSpec++] = nChan-1;
	    if (lsb) mxFreq = (1.0-maxFBW) * freqLo;  /* Lowest frequency in next bin */
	    else     mxFreq = (1.0+maxFBW) * freqLo;  /* Highest frequency in next bin */
	  } else {  /* no gap */
	    iCh = (olong)(0.5 + (mxFreq-freqLo) / (uvdesc->cdelt[uvdesc->jlocf])) - 1;
	    IFBreak[nSpec]   = iif;
	    ChBreak[nSpec++] = MIN (nChan-1, iCh);
	    if (lsb) mxFreq -= maxFBW * freqLo;
	    else     mxFreq += maxFBW * freqLo;
	  }
	  /* This bin finished? */
	  if (lsb) done = freqHi>mxFreq;   /* LSB */
	  else     done = freqHi<mxFreq;   /* USB */
	}
      }
      /* Check blown array */
      Obit_return_if_fail((nSpec<=nBreak), err, 
			  "%s: Too many coarse spectral planes, >%d for %s", 
			  routine, nBreak-1,in->name);
    } /* end loop over IF */
    
   /* End */
    if (maxFBW<0.0) nSpec--;  /* One bin per IF? */
    /* all done? */
    if ((IFBreak[nSpec-1]<(nIF-1)) || (ChBreak[nSpec-1]<(nChan-1))) {
      /* No - add final */
      IFBreak[nSpec] = nIF-1;
      ChBreak[nSpec] = nChan-1;
    } else nSpec--;  /* No more needed */
    
    /* Make sure have enough spectra */
    Obit_return_if_fail((nSpec>in->curOrder), err,
                        "%s: TOO Few channels %d, for requested order, decrease maxFBW",
                        routine, nSpec);
    /* End if spectral */
  } else nSpec = 1;

  /* Create arrays */
  in->nSpec     = nSpec;
  in->BIFSpec   = g_malloc0(nSpec*sizeof(olong));
  in->EIFSpec   = g_malloc0(nSpec*sizeof(olong));
  in->BChanSpec = g_malloc0(nSpec*sizeof(olong));
  in->EChanSpec = g_malloc0(nSpec*sizeof(olong));
  in->specFreq  = g_malloc0(nSpec*sizeof(odouble));
  in->specFreqLo= g_malloc0(nSpec*sizeof(odouble));
  in->specFreqHi= g_malloc0(nSpec*sizeof(odouble));

  /* Channel and IF increments in frequency scaling array */
  fincf  = MAX (1, (uvdesc->incf  / 3) / uvdesc->inaxes[uvdesc->jlocs]);
  fincif = MAX (1, (uvdesc->incif / 3) / uvdesc->inaxes[uvdesc->jlocs]);

  /* Loop filling arrays - Only if in->curOrder>0 */
  if (in->curOrder>0) {
    for (i=0; i<nSpec; i++) {
      /* By IF? */
      if (maxFBW<0.0) {
	in->BIFSpec[i]   = i;
	in->EIFSpec[i]   = i;
	in->BChanSpec[i] = 0;
	in->EChanSpec[i] = nChan-1;
      } else {  /* Split across IFs */
	  in->BIFSpec[i]   = IFBreak[i];
	  in->BChanSpec[i] = ChBreak[i]+1;
	  if (in->BChanSpec[i]>=nChan) {  /* Into next IF? */
	    in->BChanSpec[i]  = 0;
	    in->BIFSpec[i]   += 1         ;
	  }
	  in->EIFSpec[i]     = IFBreak[i+1];
	  in->EChanSpec[i]   = ChBreak[i+1];
      }
    } /* end loop */
      /* end multi channel out */
  } else {
    in->BIFSpec[0]   = 0;
    in->EIFSpec[0]   = nIF-1;
    in->BChanSpec[0] = 0;
    in->EChanSpec[0] = nChan-1;
    in->refFreq = in->myDesc->crval[in->myDesc->jlocf];
  }

  /* average frequency */
  sum2 = 0.0; count2 = 0;
  for (i=0; i<nSpec; i++) {
    sum  = 0.0; count = 0;
    for (iif=in->BIFSpec[i]; iif<=in->EIFSpec[i]; iif++) {
      if (in->BIFSpec[i]==in->EIFSpec[i]) { /* both ends same IF */
	iclo = in->BChanSpec[i];
	ichi = in->EChanSpec[i];
      } else if (iif==in->BIFSpec[i]) {    /* First */
 	iclo = in->BChanSpec[i];
	ichi = nChan-1;
      } else if (iif==in->EIFSpec[i]) {    /* Last */
 	iclo = 0;
	ichi = in->EChanSpec[i];
      } else {                            /* Middle */
  	iclo = 0;
	ichi = nChan-1;
      }
      for (ichan=iclo; ichan<=ichi; ichan++) {
	count++;
	sum += uvdesc->freqArr[iif*fincif + ichan*fincf];
	count2++;
	sum2 += uvdesc->freqArr[iif*fincif + ichan*fincf];
      }
    } /* end IF loop */
    if (count>0) in->specFreq[i] = sum/count;
    else  in->specFreq[i] = -1.0;
    /* Low end of bin */
    ichan = in->BChanSpec[i];
    iif   = in->BIFSpec[i];
    in->specFreqLo[i] =  uvdesc->freqArr[iif*fincif + ichan*fincf];
    /* High end of bin */
    ichan = in->EChanSpec[i];
    iif   = in->EIFSpec[i];
    in->specFreqHi[i] =  uvdesc->freqArr[iif*fincif + ichan*fincf];
  }

  /* increment number of planes by nSpec (if>1) */
  if (nSpec>1) {
    in->myDesc->inaxes[in->myDesc->jlocf] += nSpec;
    
    /* Calculate & save reference frequency */
    in->refFreq = sum2/count2;
    in->myDesc->crval[in->myDesc->jlocf] = in->refFreq;  
    /* Approximage freq increment */
    in->myDesc->cdelt[in->myDesc->jlocf] = 
      (in->specFreq[nSpec-1] - in->specFreq[0]) / (nSpec-1);  
  }
  
  /* Tell channel division if err->prtLv>=2 */
  if ((err->prtLv>=2) && doFreqDivMess) {
    doFreqDivMess = FALSE; /* only once */
    Obit_log_error(err, OBIT_InfoErr, "Frequency division of %d IFs, %d chan",
		   nIF, nChan);
    for (i=0; i<nSpec; i++) {
      Obit_log_error(err, OBIT_InfoErr, 
		     "%3d Freq %8.4f GHz, IF/ch %3d/%4d to IF/ch %3d/%4d",
		     i+1, in->specFreq[i]*1.0e-9, 
		     in->BIFSpec[i]+1, in->BChanSpec[i]+1,
		     in->EIFSpec[i]+1, in->EChanSpec[i]+1);
    }
    ObitErrLog(err); 
  }

  /* Add frequency info to descriptor */
  ObitInfoListAlwaysPut (in->myDesc->info, "NSPEC", OBIT_long, dim, &nSpec);
  for (i=0; i<nSpec; i++) {
    sprintf (keyword, "FREQ%4.4d",i+1);
    ObitInfoListAlwaysPut (in->myDesc->info, keyword, OBIT_double, 
			   dim, &in->specFreq[i]);
    sprintf (keyword, "FREL%4.4d",i+1);
    ObitInfoListAlwaysPut (in->myDesc->info, keyword, OBIT_double, 
			   dim, &in->specFreqLo[i]);
    sprintf (keyword, "FREH%4.4d",i+1);
    ObitInfoListAlwaysPut (in->myDesc->info, keyword, OBIT_double, 
			   dim, &in->specFreqHi[i]);
 }

  /* Save Alpha */
  ObitInfoListAlwaysPut (in->myDesc->info, "ALPHA", OBIT_float, dim, &alpha);
  in->alpha = alpha;
 
  /* Save Alpha reference frequency */
  ObitInfoListAlwaysPut (in->myDesc->info, "RFALPHA", OBIT_double, dim, &alphaRefF);
  in->alphaRefF = alphaRefF;
 
} /* end  ObitImageMFSetSpec */

/**
 * Get info from the descriptor for coarse channels:
 * ALPHA, RFALPHA, NTERM NSPEC, FREQ001...
 * into alpha, maxOrder-1, nSpec, specFreq
 * \param in     Pointer to object, should be fully defined
 * \param alpha  Spectral index applied
 * \param err    ObitErr for reporting errors.
 */
void ObitImageMFGetSpec (ObitImageMF *in, ObitErr *err)
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  union ObitInfoListEquiv InfoReal; 
  odouble darr[10];
  gchar keyword[12];
  olong i, nSpec, nTerm;
  /*gchar *routine = "ObitImageMFSetSpec";*/
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Number of spectral channels */
  nSpec = 0;
  ObitInfoListGetTest (in->myDesc->info, "NSPEC", &type, dim, &nSpec);
  in->nSpec     = nSpec;
  if (nSpec<=0) return;

  /* Number of spectral terms */
  nTerm = 0;
  ObitInfoListGetTest (in->myDesc->info, "NTERM", &type, dim, &nTerm);
  in->maxOrder = nTerm - 1;

  /* Alpha */
  InfoReal.flt = 0.0;   type = OBIT_float;
  ObitInfoListGetTest(in->myDesc->info, "ALPHA", &type, dim, &InfoReal);
  if (type==OBIT_double) in->alpha = (ofloat)InfoReal.dbl;
  if (type==OBIT_float)  in->alpha = (ofloat)InfoReal.flt;

  /* Reference frequency if not set */
  if (in->refFreq<1.0) in->refFreq = in->myDesc->crval[in->myDesc->jlocf];

  /* Alpha reference frequency - default to reference frequency */
  darr[0] = in->refFreq;
  ObitInfoListGetTest (in->myDesc->info, "RFALPHA", &type, dim, darr);
  in->alphaRefF = darr[0];
  if (in->alphaRefF<=0.0) in->alphaRefF = in->refFreq;

  /* Create array */
  if (in->specFreq) g_free(in->specFreq);
  in->specFreq  = g_malloc0(nSpec*sizeof(odouble));
  in->specFreqLo= g_malloc0(nSpec*sizeof(odouble));
  in->specFreqHi= g_malloc0(nSpec*sizeof(odouble));

  /* Fetch frequencies */
  for (i=0; i<nSpec; i++) {
    in->specFreq[i] = 1.0;
    sprintf (keyword, "FREQ%4.4d",i+1);
    ObitInfoListGetTest (in->myDesc->info, keyword, &type, 
			 dim, &in->specFreq[i]);
    in->specFreqLo[i] = 1.0;
    sprintf (keyword, "FREL%4.4d",i+1);
    ObitInfoListGetTest (in->myDesc->info, keyword, &type, 
			 dim, &in->specFreqLo[i]);
    in->specFreqHi[i] = 1.0;
    sprintf (keyword, "FREH%4.4d",i+1);
    ObitInfoListGetTest (in->myDesc->info, keyword, &type, 
			 dim, &in->specFreqHi[i]);
  }
 
} /* end  ObitImageMFGetSpec */

/**
 * Zero blank combined image and higher order planes
 * Planes 1-1+maxOrder
 * \param in     Pointer to object, should be fully defined
 * \param err    ObitErr for reporting errors.
 */
void ObitImageMFBlank (ObitImageMF *in, ObitErr *err)
{
  ObitFArray *imPix=NULL;
  olong i, plane[5] = {1,1,1,1,1};
  gchar *routine = "ObitImageMFBlank";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageMFIsA(in));

  /* Create pixel array */
  imPix = ObitFArrayCreate ("Blank pixels", 2, in->myDesc->inaxes);
  ObitFArrayFill (imPix, 0.0);

  /* Loop writing */
  for (i=0; i<1+in->maxOrder; i++) {
    plane[0] = i+1;
    ObitImagePutPlane ((ObitImage*)in, imPix->array, plane, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end loop writing */
 
  imPix = ObitFArrayUnref(imPix);
} /* end ObitImageMFBlank */

/**
 * Make image or beam from combined spectral channels
 * Average using weighting by 1/sigma of plane
 * If the peak in the combined image is less than 10 times the 
 * RMS in the best channel image then a correction is made:
 * use average + 0.1 times most extreme channel value > 5 sigma.
 * Note: values of lambda>0.1 may lead to instabilities causing 
 * CLEAN not to converge.
 * The extrema array is then apodized using the FT of the CLEAN beam 
 * to remove out of band noise that can drive CLEAN into oscillations.
 * Extrema are the most discrepent points from the average.
 * The RMS is raw RMS about zero of residuals from the average channel 
 * flux density.
 * \param in     Pointer to object, should be fully defined
 * \param addExt Add contribution for most extreme points
 * \param err    ObitErr for reporting errors.
 */
void ObitImageMFCombine (ObitImageMF *in, gboolean addExt, ObitErr *err)
{
  ObitFArray *imPix=NULL, *extPix=NULL, *extConvl=NULL;
  olong i, plane[5] = {1,1,1,1,1}, pos[2];
  ofloat norm, sigma, sigClip=5.0, lambda=0.1, testSigma=10.0, beamArea, cells;
  ofloat wt, sumWt, p1, p2, minSigma, fblank = ObitMagicF();
  gboolean *bad=NULL, linPol=FALSE, doExt=FALSE;
  gchar *routine = "ObitImageMFCombine";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageMFIsA(in));

  /* Linear poln? - More aggressive extrema suppression */
  linPol = (in->myDesc->crval[ in->myDesc->jlocs]>=1.9) && 
           (in->myDesc->crval[ in->myDesc->jlocs]<=3.1);
  if (linPol) {
    testSigma = 5.0;  /* How many sigma residual is important */
    lambda    = 0.25; /* Reduction factor of extrama correction */
  }

  /* Loop accumulating average */
  sumWt    = 0.0;
  minSigma = 1.0e20;
  bad = g_malloc0(in->nSpec*sizeof(gboolean));
  for (i=1+in->maxOrder; i<1+in->maxOrder+in->nSpec; i++) {
    plane[0] = i+1;
    ObitImageGetPlane ((ObitImage*)in, NULL, plane, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Create pixel arrays for accumulation first pass */
    if (i==1+in->maxOrder) {
      imPix  = ObitFArrayCreate ("accum", 2, in->myDesc->inaxes);
      if (addExt) extPix = ObitFArrayCreate ("accum", 2, in->myDesc->inaxes);
    }
    /* Weight by 1/sigma */
    wt = ObitFArrayRMS(in->image);
    /* Values? */
    if ((wt==fblank) || (wt<=0.0)) {
      bad[i-1-in->maxOrder] = TRUE;
      continue;
    } else {
      bad[i-1-in->maxOrder] = FALSE;
    }
    minSigma = MIN (minSigma, wt);  /* Minimum valid sigma */
    wt = 1.0 / wt;
    sumWt += wt;
    ObitFArraySMul (in->image, wt);

    /* Accumulate */
    ObitFArrayAdd(imPix, in->image, imPix);

  } /* end loop accumulating */
 
  /* Normalize */
  norm = 1.0 / sumWt;
  ObitFArraySMul (imPix, norm);

  /* Add scaled extrema convolved with dirty beam 
   scaling by sigma keeps poor planes from overly affecting 
   the results */
  if (addExt && in->myBeam) {
    /* Loop accumulating extrema */
    for (i=1+in->maxOrder; i<1+in->maxOrder+in->nSpec; i++) {
      if (!bad[i-1-in->maxOrder]) {
	plane[0] = i+1;
	ObitImageGetPlane ((ObitImage*)in, NULL, plane, err);
	if (err->error) Obit_traceback_msg (err, routine, in->name);
	
	/* Subtract average */
	ObitFArraySub(in->image, imPix, in->image);
	
	/* Accumulate extrema - first clip below sigClip */
	sigma = ObitFArrayRMS0(in->image);
	ObitFArrayInClip (in->image, -sigClip*sigma, sigClip*sigma, 0.0);
	/* Weight by 1/sigma */
	wt = 1.0 / sigma;
	ObitFArraySMul (in->image, wt);   /* Scale by weight */
	/* ObitFArrayAbs (in->image);        Absolute value */
	ObitFArrayExtArr (extPix, in->image, extPix);
      } /* end if valid */
    } /* end loop accumulating */

    /* Normalize extrema image by minimum sigma */
    ObitFArraySMul (extPix, minSigma);

   /* Read dirty beam */
    plane[0] = 1;
    ObitImageGetPlane ((ObitImage*)in->myBeam, NULL, plane, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    /* Scale extrema - include effect of beam area */
    cells = fabs(in->myDesc->cdelt[0]);
    beamArea = 1.1331*(in->myDesc->beamMaj/cells)*(in->myDesc->beamMin/cells);
    /* If beam not available use 35 pixels */
    if (beamArea<0.1)  beamArea = 35.0;
    /* add if peak in imPix exceeds that in extPix - something occasionally goes wrong */
    p1 = ObitFArrayMaxAbs (imPix, pos);
    p2 = ObitFArrayMaxAbs (extPix, pos);
    doExt = (p1<p2/sqrt(in->nSpec)) && (p2>testSigma*minSigma);  /* Add extrema? */
    if (doExt) {
      ObitFArraySMul (extPix, lambda/beamArea); /* Scale residuals by beam area */
      /* Adjust signs to match imPix */
      ObitFArraySign (imPix, extPix);
      /* Convolve extrema with dirty beam */
      extConvl = ObitFArrayUtilConvolve (extPix, ((ObitImage*)in->myBeam)->image, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);

      /* Add extrema in a way to make CLEANing more likely */
       ObitFArrayAdd (imPix, extConvl, imPix);
    }
    if (err->prtLv>=4) {  /* Diagnostics */
      Obit_log_error(err, OBIT_InfoErr, "%s: Add extrema=%d pI=%f pRE=%f minSigma=%f",
		     routine, doExt, p1, p2, minSigma);
      ObitErrLog(err); 
    }
    /* Cleanup */
    ((ObitImage*)in->myBeam)->image  = ObitFArrayUnref(((ObitImage*)in->myBeam)->image);
    extConvl = ObitFArrayUnref(extConvl);
  }

  /* Write */
  plane[0] = 1;
  ObitImagePutPlane ((ObitImage*)in, imPix->array, plane, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  if (bad) g_free(bad);
  imPix     = ObitFArrayUnref(imPix);
  extPix    = ObitFArrayUnref(extPix);
  in->image = ObitFArrayUnref(in->image);
} /* end ObitImageMFCombine */

/**
 * Fit Spectra to coarse spectral planes in an ObitImageMF
 * Possibly uses threads, looping over rows.
 * \param in      Image to fit, info may have
 * \li nOrder OBIT_int (1,1,1) Maximum order to fit
 *            If absent use maxOrder
 * \param antSize If > 0 make primary beam corrections assuming antenna 
 *                diameter (m) antSize
 * \param err     Obit error stack object.
 */
void ObitImageMFFitSpec (ObitImageMF *in, ofloat antSize, ObitErr *err)
{
  olong ithread, maxThread, nThreads=1;
  FitSpecFuncArg **targs=NULL;
  ObitImage **inArray=NULL, **outArray=NULL; 
  ofloat sigma;
  olong i, iy, ny, nterm=1, nOrder, nTh, plane[5] = {1,1,1,1,1};
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOSize IOBy = OBIT_IO_byRow;
  gboolean OK;
  gchar *routine = "ObitImageMFFitSpec";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageMFIsA(in));

  /* Try reading spectral info if missing */
  if (in->nSpec<=0) ObitImageMFGetSpec (in, err);
  if (err->error) goto cleanup;

  /* Make sure have spectra */
  Obit_return_if_fail((in->nSpec>0), err, 
  		      "%s: NO spectra in %s", 
		      routine, in->name);

  /* Tell about it */
  Obit_log_error(err, OBIT_InfoErr, "Fitting pixel spectra");
  ObitErrLog(err); 

  nOrder = in->maxOrder;
  ObitInfoListGetTest (in->info, "nOrder", &type, dim, &nOrder);

  /* Create arrays of Images to handle input and output I/O */
  nterm = in->maxOrder+1;    /* Number of output spectral terms */
  inArray  = g_malloc0(in->nSpec*sizeof(ObitImage*));
  outArray = g_malloc0(nterm*sizeof(ObitImage*));
  trc[0] = in->myDesc->inaxes[0]; trc[1] = in->myDesc->inaxes[1]; 
  /* Input */
  for (i=0; i<in->nSpec; i++) {
    inArray[i] = newObitImage("Input");	
    /* Copy info */
    inArray[i]->info = ObitInfoListCopyData(in->info, inArray[i]->info);
    inArray[i]->extBuffer = TRUE;  /* external buffer */
    /* Set up for IO */
    dim[0] = 1;
    ObitInfoListAlwaysPut (inArray[i]->info, "IOBy", OBIT_long, dim, &IOBy); 
    dim[0] = 7;
    blc[2] = trc[2] = 2+in->maxOrder+i;
    ObitInfoListAlwaysPut (inArray[i]->info, "BLC", OBIT_long, dim, blc); 
    ObitInfoListAlwaysPut (inArray[i]->info, "TRC", OBIT_long, dim, trc);
    /* Open */
    ObitImageOpen (inArray[i], OBIT_IO_ReadOnly, err);
    if (err->error) goto cleanup;
  }
  /* output */
  for (i=0; i<nterm; i++) {
    outArray[i] = newObitImage("Output");
    outArray[i]->info = ObitInfoListCopyData(in->info, outArray[i]->info);
    outArray[i]->extBuffer = TRUE;
    /* Set up for IO */
    dim[0] = 1;
    ObitInfoListAlwaysPut (outArray[i]->info, "IOBy", OBIT_long, dim, &IOBy); 
    dim[0] = 7;
    blc[2] = trc[2] = i+1;
    ObitInfoListAlwaysPut (outArray[i]->info, "BLC", OBIT_long, dim, blc); 
    ObitInfoListAlwaysPut (outArray[i]->info, "TRC", OBIT_long, dim, trc);
    /* Open */
    ObitImageOpen (outArray[i], OBIT_IO_ReadWrite, err);
    if (err->error) goto cleanup;
  }

  /* Setup Threading */
  /* Only thread large cases */
  if (in->myDesc->inaxes[0]>200) maxThread = 1000;
  else maxThread = 1;
  nThreads = MakeFitSpecArgs (in, maxThread, antSize, nOrder, &targs, err);
  if (err->error) goto cleanup;

  /* Read through all input channels getting sigmas - save on thread arguments */
  for (i=0; i<in->nSpec; i++) {
    plane[0] = 2+in->maxOrder+i;
    ObitImageGetPlane ((ObitImage*)in, NULL, plane, err);
    sigma = ObitFArrayRMS(in->image);
    if (err->error) goto cleanup;
    /* Save on thread arguments */
    for (ithread=0; ithread<nThreads; ithread++)
      targs[ithread]->sigma[i] = sigma;
  } /* end loop getting sigmas */
  in->image = ObitFArrayUnref(in->image);  /* Free image buffer */

  /* Loop over rows - one row per thread */
  ny = in->myDesc->inaxes[1];
  for (iy=0; iy<ny; iy += nThreads) {

    nTh = MIN (nThreads, (ny-iy));
    /* Set up thread arguments */
    for (ithread=0; ithread<nTh; ithread++) {
      targs[ithread]->iy = iy+ithread;
    }

    /* Read data */
    for (i=0; i<in->nSpec; i++) {
      for (ithread=0; ithread<nTh; ithread++) {
	ObitImageRead(inArray[i], targs[ithread]->inData[i], err);
      }
      if (err->error) goto cleanup;
    }

    /* operation possibly in threads */
    OK = ObitThreadIterator (in->thread, nTh, 
			     (ObitThreadFunc)ThreadFitSpec, (gpointer**)targs);

    /* Check for problems */
    if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
    
    /* Write data */
    for (i=0; i<nterm; i++) {
      for (ithread=0; ithread<nTh; ithread++) {
	ObitImageWrite(outArray[i], targs[ithread]->outData[i], err);
      }
      if (err->error) goto cleanup;
   }

  } /* end loop over list */
  
  /* Cleanup */
 cleanup:
  KillFitSpecArgs (nThreads, targs);
  ObitThreadPoolFree (in->thread);  
  if (inArray) {
    for (i=0; i<in->nSpec; i++) {
      ObitImageClose (inArray[i], err);
      if (err->error) goto error;
      inArray[i] = ObitImageUnref(inArray[i]);
    }
    g_free(inArray);
  }
  if (outArray) {
    for (i=0; i<nterm; i++) {
      ObitImageClose (outArray[i], err);
      if (err->error) goto error;
      outArray[i] = ObitImageUnref(outArray[i]);
    }
    g_free(outArray);
  }
 error:
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitImageMFFitSpec */

/**
 * Fit Spectra to coarse spectral planes in an ObitImageMF, writes a separate output image
 * Possibly uses threads, looping over rows.
 * \param in      Image to fit, info may have
 * \li            nterm      int scalar Number of terms in spectral fit, 2=SI, 3=curvature [def 2]
 * \li            refFreq    double scalar Reference frequency for fit [def ref for in]
 * \li            maxChi2    float scalar Max. Chi Sq for accepting a partial spectrum [def 2.0]
 * \li            doError    boolean scalar If true do error analysis [def False]
 * \li            doBrokePow boolean scalar If true do broken power law (3 terms). [def False]
 * \li            calFract   float (?,1,1) Calibration error as fraction of flux
 *                           One per frequency or one for all, def 0.05
 * \li            doPBCor    boolean scalar If true do primary beam correction. [def False]
 * \li            PBmin      float (?,1,1) Minimum beam gain correction
 *                    One per frequency or one for all, def 0.05,
 *                    1.0 => no gain corrections
 * \li            antSize    float (?,1,1) Antenna diameter (m) for gain corr, 
 *                    One per frequency or one for all, def 25.0
 * \param out     Output Image into which to write parameters
 * \param err     Obit error stack object.
 */
void ObitImageMFFitSpec2 (ObitImageMF *in, ObitImageMF *out, ObitErr *err)
{
  ObitSpectrumFit* fitter=NULL;
  ofloat antsize, pbmin;
  olong i, plane[5] = {1,1,1,1,1};
  olong nOut, iplane, naxis[2], nterm=2; 
  odouble refFreq=-1.0;
  ofloat maxChi2=2.0, *PBmin=NULL, *antSize=NULL, *calFract=NULL; 
  gboolean doGain, doError=FALSE, doBrokePow=FALSE, doPBCor=FALSE;
  union ObitInfoListEquiv InfoReal; 
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1}, PBdim[MAXINFOELEMDIM], ASdim[MAXINFOELEMDIM];
  gchar *fitParms[] = {"refFreq","maxChi2","doError","doPBCor","doBrokePow","calFract","PBmin","antSize",
		       NULL};
  gchar  *today=NULL, *SPECLOGF = "SPECLOGF";
  gchar *routine = "ObitImageMFFitSpec2";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageMFIsA(in));

  /* Fully instantiate input */
  ObitImageMFFullInstantiate (in, TRUE, err);
  if (err->error) goto cleanup;

  /* Try reading spectral info if missing */
  if (in->nSpec<=0) ObitImageMFGetSpec (in, err);
  if (err->error) goto cleanup;

  /* Make sure have spectra */
  Obit_return_if_fail((in->nSpec>0), err, 
  		      "%s: NO spectra in %s", 
		      routine, in->name);

  /* Tell about it */
  Obit_log_error(err, OBIT_InfoErr, "Fitting pixel spectra");
  ObitErrLog(err); 

  /* Fitting parameters */
  ObitInfoListGetTest (in->info, "nterm",      &type, dim, &nterm);
  /* Get Reference frequency , default input ref. freq. */
  InfoReal.dbl = in->myDesc->crval[in->myDesc->jlocf]; 
  type = OBIT_double;  in->refFreq = InfoReal.dbl;
  ObitInfoListGetTest(in->info, "refFreq", &type, dim, &InfoReal);
  if (type==OBIT_float)       refFreq = InfoReal.flt;
  else if (type==OBIT_double) refFreq = (ofloat)InfoReal.dbl;
  /* Min Chi^2 for fit */
  InfoReal.flt = maxChi2; type = OBIT_float;
  ObitInfoListGetTest(in->info, "maxChi2", &type, dim, &InfoReal);
  if (type==OBIT_float)       maxChi2 = InfoReal.flt;
  else if (type==OBIT_double) maxChi2 = (ofloat)InfoReal.dbl;
  /* Want Error analysis? */
  InfoReal.itg = (olong)doError; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doError", &type, dim, &InfoReal);
  doError = InfoReal.itg;
   /* Want Broken power law ? */
  InfoReal.itg = (olong)doBrokePow; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doBrokePow", &type, dim, &InfoReal);
  doBrokePow = InfoReal.itg;
  /* Want primary beam correction? */
  InfoReal.itg = (olong)doPBCor; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doPBCor", &type, dim, &InfoReal);
  doPBCor = InfoReal.itg;
  ObitInfoListGetTest (in->info, "antSize",    &type, ASdim, (gpointer)&antSize);
  ObitInfoListGetP    (in->info, "PBmin",      &type, PBdim, (gpointer)&PBmin);
  
  /* Create Fitter */
  fitter = ObitSpectrumFitCreate("Fitter", nterm);
  fitter->refFreq    = refFreq;
  fitter->maxChi2    = maxChi2;
  fitter->doError    = doError;
  fitter->doBrokePow = doBrokePow;
  fitter->calFract   = calFract;

  /* Get fitting parameters, copy to fitter */
  ObitInfoListCopyList (in->info, fitter->info, fitParms);

  /* Determine number of frequency planes and initialize fitter */
  fitter->nfreq = in->nSpec;
  fitter->BeamShapes = g_malloc0(fitter->nfreq*sizeof(ObitBeamShape*));
  for (i=0; i<fitter->nfreq; i++) fitter->BeamShapes[i] = NULL;
  fitter->RMS        = g_malloc0(fitter->nfreq*sizeof(ofloat));
  fitter->calFract   = g_malloc0(fitter->nfreq*sizeof(ofloat));
  fitter->inFArrays  = g_malloc0(fitter->nfreq*sizeof(ObitFArray*));
  fitter->freqs      = g_malloc0(fitter->nfreq*sizeof(odouble));
  
  /* How many output planes? */
  if (fitter->doError) nOut = 1+fitter->nterm*2;
  else nOut = fitter->nterm;
  
  /* Define term arrays */
  fitter->outFArrays = g_malloc0((nOut)*sizeof(ObitFArray*));
  for (i=0; i<nOut; i++) fitter->outFArrays[i] = NULL;

  /* Image size */
  fitter->nx = in->myDesc->inaxes[0];
  fitter->ny = in->myDesc->inaxes[1];
  naxis[0] = (olong)fitter->nx;  naxis[1] = (olong)fitter->ny; 
  for (i=0; i<fitter->nfreq; i++) fitter->inFArrays[i]  = ObitFArrayCreate (NULL, 2, naxis);
  for (i=0; i<nOut; i++)          fitter->outFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);

  /* Calibration error */
  for (i=0; i<fitter->nfreq; i++) fitter->calFract[i] = 1.0e-5; /* default */
  if (ObitInfoListGetP(fitter->info, "calFract", &type, dim, (gpointer)&calFract)) {
    if (dim[0]>=fitter->nfreq) for (i=0; i<fitter->nfreq; i++) fitter->calFract[i] = calFract[i];
    else for (i=0; i<fitter->nfreq; i++) fitter->calFract[i] = calFract[0];
  }

  /* Output Image descriptor */
  out->myDesc = ObitImageDescCopy (in->myDesc, out->myDesc, err);
  if (err->error) goto cleanup;

  /* Change third axis to type "SPECLOGF" and leave the reference frequency
     as the "CRVAL" */
  out->myDesc->inaxes[out->myDesc->jlocf] =  nOut;
  out->myDesc->crval[out->myDesc->jlocf]  =  fitter->refFreq;
  out->myDesc->crpix[out->myDesc->jlocf]  =  1.0;
  out->myDesc->cdelt[out->myDesc->jlocf]  =  1.0;
  strncpy (out->myDesc->ctype[out->myDesc->jlocf], SPECLOGF, IMLEN_KEYWORD);
  out->myDesc->bitpix = -32;  /* Float it */

  /* Creation date today */
  today = ObitToday();
  strncpy (out->myDesc->date, today, IMLEN_VALUE);
  if (today) g_free(today);

  /* Fully instantiate output */
  ObitImageMFFullInstantiate (out, FALSE, err);
  if (err->error) goto cleanup;

  /* Loop reading planes */
  for (iplane=0; iplane<fitter->nfreq; iplane++) {
    plane[0] = iplane+in->maxOrder+2; 
    ObitImageGetPlane ((ObitImage*)in, fitter->inFArrays[iplane]->array, plane, err);
    if (err->error) goto cleanup;

    /* Get BeamShape */
    antsize = 25.0;
    if (antSize) antsize = antSize[0];
    if (antSize && (ASdim[0]>=fitter->nfreq)) antsize = antSize[iplane];
    pbmin   = 0.05;
    if (PBmin) pbmin = PBmin[0];
    if (PBmin && (PBdim[0]>=fitter->nfreq)) pbmin = PBmin[iplane];
    doGain = pbmin<0.999;
    fitter->BeamShapes[iplane] = ObitBeamShapeCreate ("BS", (ObitImage*)in, pbmin, antsize, doGain);

    /* Plane RMS */
    fitter->RMS[iplane] = ObitFArrayRMS(fitter->inFArrays[iplane]);
    /* Frequency */
    fitter->freqs[iplane] = in->specFreq[iplane];
    
  } /* end loop reading planes */

  /* Do actual fitting */
  ObitSpectrumFitter (fitter, err);
  if (err->error) goto cleanup;
  if (err->error) Obit_traceback_msg (err, routine, fitter->name);

  /* Update output header to reference Frequency */
  out->myDesc->crval[out->myDesc->jlocf] = fitter->refFreq;
  fitter->outDesc = ObitImageDescRef(in->myDesc); /* Desc to fitter */

  /* Write output */
  ObitSpectrumWriteOutput(fitter, (ObitImage*)out, err);
  if (err->error) goto cleanup;

  /* Cleanup */
 cleanup:
  fitter = ObitSpectrumFitUnref(fitter);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitImageMFFitSpec2 */

/**
 * Get pointer to wideband image beam
 * returns beam pointer and the plane in the image
 * \param inn    Image whose beam name is to be set 
 * \param beamNo Which order Beam [0-rel]
 *               0=dirty, 1=spectral index, 2=curvature...
 * \param plane  [out] Plane number for beam
 * \param err    Obit error structure
 * \return pointer to beam, NULL if not defined.
 */
ObitImage* ObitImageMFGetBeam (ObitImage *inn, olong beamNo, 
			       olong plane[5], ObitErr *err) 
{
  ObitImageMF *in = (ObitImageMF*)inn;
  ObitImage *theBeam;
  olong i;
  gchar *routine = "ObitImageMFGetBeam";

  /* error checks */
  g_assert (ObitImageMFIsA(inn));

  for (i=0; i<5; i++) plane[i] = 1;  /* Initialize plane */

  theBeam  = (ObitImage*)in->myBeam;

  plane[0] = MAX (1, beamNo+1);

  /* Make sure imaging this order */
  Obit_retval_if_fail((beamNo<=in->curOrder), err, NULL,
  		      "%s: Beam %d > %d not available for %s", 
		      routine, beamNo, in->curOrder, in->name);

  /* Make sure this plane exists */
  Obit_retval_if_fail(((theBeam) && (theBeam->myDesc) &&
		       (theBeam->myDesc->inaxes[2]>=plane[0])), err, NULL,
  		      "%s: Beam %d not available for %s", 
		      routine, beamNo, in->name);

  return theBeam;
} /* end ObitImageMFGetBeam */

/**
 * Get current highest order of image beam here 0
 * \param inn  Image whose beam order is sought
 * \return order number
 */
olong ObitImageMFGetBeamOrder (ObitImage *inn) 
{
  return 0;
} /* end ObitImageMFGetBeam */

/**
 * Get Frequency of current image plane, works for eitther Image or ImageMF
 * \param image  ImageMF in question
 * \return Freq Frequency (Hz), 0.0 if none in descriptor or list
 */
odouble ObitImageMFGetPlaneFreq (ObitImage *image) 
{
  odouble Freq=0.0;
  olong nterm =0,plane = image->myDesc->plane;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar keyword[20];

  /* Have ImageMF keywords? */
  if (!ObitInfoListGetTest(image->myDesc->info, "NTERM", &type, dim, &nterm)) {
    if (image->myDesc->jlocf>=0) 
      Freq = image->myDesc->crval[image->myDesc->jlocf];
    return Freq;
  }
  /* Get from header List or header if missing, term planes will get main header freq. */
  sprintf(keyword,"FREQ%4.4d", plane-nterm);
  if (!ObitInfoListGetTest(image->myDesc->info, keyword, &type, dim, &Freq)) {
    if (image->myDesc->jlocf>=0) 
      Freq = image->myDesc->crval[image->myDesc->jlocf];
  }
  /* DEBUG */
  return Freq;
} /* end ObitImageMFGetPlaneFreq */

/*-------Private functions called by ObitData class ------*/
/** Private: Zap */
static ObitData* ObitDataImageMFZap (ObitData *in, ObitErr *err)
{
  return (ObitData*)ObitImageMFZap ((ObitImage*)in, err);
} /* end ObitDataImageMFZap */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitImageMFClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();
  myClassInfo.hasScratch    = TRUE; /* Scratch files allowed */

  /* Set function pointers */
  ObitImageMFClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitImageMFClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitImageMFClassInfoDefFn ( gpointer inClass)
{
  ObitImageMFClassInfo *theClass = (ObitImageMFClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitImageMFClassInit;
  theClass->newObit       = (newObitFP)newObitImageMF;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitImageMFClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitImageMFGetClass;
  theClass->newObitImageScratch  = 
    (newObitImageScratchFP)newObitImageMFScratch;
  theClass->ObitCopy      = (ObitCopyFP)ObitImageMFCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitImageMFClear;
  theClass->ObitInit      = (ObitInitFP)ObitImageMFInit;
  theClass->ObitImageGetBeam = 
    (ObitImageGetBeamFP)ObitImageMFGetBeam;
  theClass->ObitImageGetBeamOrder = 
    (ObitImageGetBeamOrderFP)ObitImageMFGetBeamOrder;
  theClass->ObitImageGetPlaneFreq = 
    (ObitImageGetPlaneFreqFP)ObitImageMFGetPlaneFreq;

  /* Function pointers referenced from ObitData class */
  theClass->ObitDataZap     = (ObitDataZapFP)ObitDataImageMFZap;

} /* end ObitImageMFClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitImageMFInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageMF *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->nSpec     = 0;
  in->BIFSpec   = NULL;
  in->EIFSpec   = NULL;
  in->BChanSpec = NULL;
  in->EChanSpec = NULL;
  in->specFreq  = NULL;
  in->specFreqLo= NULL;
  in->specFreqHi= NULL;

} /* end ObitImageMFInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitImageMF* cast to an Obit*.
 */
void ObitImageMFClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageMF *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Delete underlying higher order beam files if isScratch */
  /* delete this class members */
  if (in->BIFSpec)   g_free(in->BIFSpec);
  if (in->EIFSpec)   g_free(in->EIFSpec);
  if (in->BChanSpec) g_free(in->BChanSpec);
  if (in->EChanSpec) g_free(in->EChanSpec);
  if (in->specFreq)  g_free(in->specFreq);
  if (in->specFreqLo)  g_free(in->specFreqLo);
  if (in->specFreqHi)  g_free(in->specFreqHi);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitImageMFClear */

/**
 * Make arguments for Threaded Spectrum fitting
 * \param in         MF Image to be fitted
 * \param maxThread  Maximum desirable no. threads
 * \param antSize    if >0 then make primary beam corrections
 * \param nOrder     Fitting order
 * \param args       [out] Created array of FitSpecFuncArg, 
 *                   delete with KillFitSpecArgs
 * \return number of elements in args.
 */
static olong MakeFitSpecArgs (ObitImageMF *image, olong maxThread,
			      ofloat antSize, olong nOrder,
			      FitSpecFuncArg ***args, 
			      ObitErr *err)
{
  olong i, j, nx, nThreads;
  ofloat pbmin = 0.01;
  gboolean doPBCor;
  gchar *routine = "MakeFitSpecArgs";
  
  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(image->thread));
  nThreads = MIN (nThreads, maxThread);
  nx = image->myDesc->inaxes[0];
  doPBCor = antSize>0.0;   /* Primary beam correction? */

  /* Initialize threadArg array */
  *args = g_malloc0(nThreads*sizeof(FitSpecFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*args)[i] = g_malloc0(sizeof(FitSpecFuncArg)); 
  
  for (i=0; i<nThreads; i++) {
    (*args)[i]->err       = err;
    (*args)[i]->desc      = ObitImageDescRef(image->myDesc);
    if (nThreads>1)
      (*args)[i]->ithread   = i;
    else
      (*args)[i]->ithread   = -1;
    (*args)[i]->thread    = image->thread;
    (*args)[i]->antSize   = antSize;
    (*args)[i]->doPBCorr  = doPBCor;
    (*args)[i]->BeamShape = ObitBeamShapeCreate ("BS", (ObitImage*)image, pbmin, antSize, doPBCor);
    (*args)[i]->nSpec     = image->nSpec;
    (*args)[i]->nOrder    = nOrder;
    (*args)[i]->alpha     = image->alpha;
    (*args)[i]->sigma     = g_malloc0(image->nSpec*sizeof(ofloat));
    (*args)[i]->workFlux  = g_malloc0(image->nSpec*sizeof(ofloat));
    (*args)[i]->workSigma = g_malloc0(image->nSpec*sizeof(ofloat));
    (*args)[i]->Freq      = g_malloc0(image->nSpec*sizeof(odouble));
    (*args)[i]->inData    = g_malloc0(image->nSpec*sizeof(ofloat*));
    for (j=0; j<image->nSpec; j++) {
      (*args)[i]->Freq[j]   = image->specFreq[j];
      (*args)[i]->inData[j] = g_malloc0(nx*sizeof(ofloat*));
    }
    (*args)[i]->outData    = g_malloc0((image->maxOrder+1)*sizeof(ofloat*));
    for (j=0; j<(image->maxOrder+1); j++) {
      (*args)[i]->outData[j] = g_malloc0(nx*sizeof(ofloat*));
    }
    (*args)[i]->fitArg    = ObitSpectrumFitMakeArg (image->nSpec, nOrder+1, 
						    image->refFreq, image->specFreq,
						    FALSE, 
						    &(*args)[i]->fitResult, err);
    if (err->error) Obit_traceback_val (err, routine, image->name, nThreads);
  }

  return nThreads;
} /*  end MakeFitSpecArgs */

/**
 * Delete arguments for Threaded FitSpec
 * \param nargs      number of elements in args.
 * \param args       Array of FitSpecFuncArg, type FitSpecFuncArg
 */
static void KillFitSpecArgs (olong nargs, FitSpecFuncArg **args)
{
  olong i, j;

  if (args==NULL) return;
  for (i=0; i<nargs; i++) {
    if (args[i]) {
      args[i]->desc      = ObitImageDescUnref(args[i]->desc);
      args[i]->BeamShape = ObitBeamShapeUnref(args[i]->BeamShape);
      if (args[i]->sigma)     g_free(args[i]->sigma);
      if (args[i]->workFlux)  g_free(args[i]->workFlux);
      if (args[i]->workSigma) g_free(args[i]->workSigma);
      if (args[i]->Freq)      g_free(args[i]->Freq);
      if (args[i]->outData) {
	for (j=0; j<(args[i]->nOrder+1); j++) 
	  if (args[i]->outData[j]) g_free(args[i]->outData[j]);
	g_free(args[i]->outData);
      }
      if (args[i]->inData) {
	for (j=0; j<args[i]->nSpec; j++) 
	  if (args[i]->inData[j]) g_free(args[i]->inData[j]);
	g_free(args[i]->inData);
      }
      if (args[i]->fitResult) g_free(args[i]->fitResult);
      ObitSpectrumFitKillArg(args[i]->fitArg);
      g_free(args[i]);
    }
  }
  g_free(args);
} /*  end KillFitSpecArgs */

/**
 * Threaded primary beam correction and spectrum fitting
 * Processes one row of an image
 * Arguments are given in the structure passed as arg
 * \param arg Pointer to CLEANFuncArg argument with elements:
 * \li alpha    Spectral index correction previously applied to data.
 * \li nSpec    Number of spectral channels
 * \li Freq     Array (nSpec) ofchannel frequencies (Hz).
 * \li sigma    Array of sigma for each channel
 * \li BeamShape BeamShape object
 * \li inData    Array  (nSpec) of rows of pixel values in image in each of nSpec planes
 * \li workFlux  Work array of flux for each channel
 * \li workSigma Work array of sigma for each channel
 * \li iy        0-rel row number
 * \li desc      Input image descriptor
 * \li doPBCorr  Correct for primary Beam?
 * \li antSize   Antenna diameter (m)
 * \li nOrder    Max. number of orders to fit 
 * \li outData   (nOrder+1) Arrays of fitted values 
 * \li fitArg    Fitting argument
 * \li fitResult Fitting result array
 * \li ithread   Thread number, <0 -> no threading 
 * \li thread    Obit Thread object 
 * \li err       Obit error stack object
 * \return NULL
 */
gpointer ThreadFitSpec (gpointer args)
{
  FitSpecFuncArg *largs  = (FitSpecFuncArg*)args;
  ofloat alpha             = largs->alpha;
  olong   nSpec            = largs->nSpec;
  odouble *Freq            = largs->Freq;
  ofloat *sigma            = largs->sigma;
  ObitBeamShape *BeamShape = largs->BeamShape;
  ofloat **inData          = largs->inData;
  ofloat *workFlux         = largs->workFlux;
  ofloat *workSigma        = largs->workSigma;
  olong  iy                = largs->iy;
  ObitImageDesc *desc      = largs->desc;
  gboolean doPBCorr        = largs->doPBCorr;
  /*ofloat antSize           = largs->antSize;*/
  olong  nOrder            = largs->nOrder;
  ofloat **outData         = largs->outData;
  gpointer fitArg          = largs->fitArg;
  ofloat *fitResult        = largs->fitResult;
  olong ithread            = largs->ithread;
  ObitThread *thread       = largs->thread;
  ObitErr  *err            = largs->err;

  /* Local */
  olong ix, i, nterm;
  ofloat PBCorr = 1.0, pixel[2];
  gboolean doJinc;
  ObitBeamShapeClassInfo *BSClass;
  odouble Angle=0.0, pos[2];
  ofloat fblank = ObitMagicF();
  gchar *routine = "ThreadFitSpec";

  nterm = nOrder+1;        /* Number of fitted spectral terms */
  doJinc = Freq[0]>1.0e9;  /* Use polynomial or Jinc for PB corr? */

  /* Loop over row */
  for (ix=0; ix<desc->inaxes[0]; ix++) {
    /* Primary beam correction? */
    if (doPBCorr) {
      /* Primary beam stuff - Distance from Center */
      pixel[0] = (ofloat)(ix+1.0); pixel[1] = (ofloat)(iy+1.0);
      ObitImageDescGetPos(desc, pixel, pos, err);
      if (err->error) {
	ObitThreadLock(thread);  /* Lock against other threads */
	Obit_log_error(err, OBIT_Error,"%s Error with position %s",
		       routine, routine);
	ObitThreadUnlock(thread); 
	return NULL;
      }
      BSClass = (ObitBeamShapeClassInfo*)(BeamShape->ClassInfo);
      Angle   = BSClass->ObitBeamShapeAngle(BeamShape, pos[0], pos[1], 0.0);
      
      /* Load arrays */
      for (i=0; i<nSpec; i++) {
	BeamShape->refFreq = Freq[i];  /* Set frequency */
	PBCorr  = BSClass->ObitBeamShapeGainSym(BeamShape, Angle);
	if (inData[i][ix]!=fblank) workFlux[i]  = inData[i][ix] / PBCorr;
	else workFlux[i] = fblank;
	if (sigma[i]>1.0e-9) workSigma[i] = sigma[i] / (PBCorr);
 	else {workSigma[i] = 1.0e10; workFlux[i] = fblank;}
      }
    } else { /* No PB correction */
      /* Load arrays */
      for (i=0; i<nSpec; i++) {
	workFlux[i]  = inData[i][ix];
	if (sigma[i]>1.0e-9) workSigma[i] = sigma[i];
  	else {workSigma[i] = 1.0e10; workFlux[i] = fblank;}
      }
    }
    
    /* Fit spectrum */
    ObitSpectrumFitSingleArg (fitArg, workFlux, workSigma, fitResult);

    /* Spectral index correction */
    if (nterm>=2) fitResult[1] += alpha;
    
    /* Save values */
    for (i=0; i<nterm; i++) {
      outData[i][ix] = fitResult[i];
      if (fitResult[0]==fblank)  outData[i][ix] = fblank;
    }
    
  }  /* end loop over row */

  if (err->error) goto finish;  /* Prior error? */


  /* Indicate completion */
 finish:
  if (largs->ithread>=0)
    ObitThreadPoolDone (thread, (gpointer)&ithread);
  return NULL;
} /* end ThreadFitSpec */

