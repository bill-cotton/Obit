/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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

#include "ObitImageDesc.h"
#include "ObitImageSel.h"
#include "ObitImageWB.h"
#include "ObitFArrayUtil.h"
#include "ObitFFT.h"
#include "ObitMem.h"
#include "ObitSystem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageWB.c
 * ObitImageWB class function definitions.
 * This class is derived from the ObitData base class.
 *
 * This class contains an astronomical image and allows access.
 * An ObitImageWB is the front end to a persistent disk resident structure.
 *
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitImageWB";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitImageGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitImageWBClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitImageWBClassInfo myClassInfo = {FALSE};

/*---------------Private structures----------------*/
/* Threaded function argument for image decomposition */
typedef struct {
  /* ObitThread to use */
  ObitThread *thread;
  /* First ObitFArray, residual * beam0 */
  ObitFArray *in1;
  /* Second ObitFArray residual * beam1 */
  ObitFArray *in2;
  /* Third ObitFArray residual * beam2 */
  ObitFArray *in3;
  /* ObitFArray for decomposition */
  ObitFArray *out;
  /* First element (1-rel) number */
  olong        first;
  /* Highest element (1-rel) number */
  olong        last;
  /* which spectral term estimated? */
  olong        iterm;
  /* which spectral order? */
  olong        norder;
  /* Scalar arguments */
  ofloat A00, A10, A01, A11, A20, A02, A21, A12, A22, iDelta;
  /* thread number  */
  olong        ithread;
} DecompFuncArg;

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitImageWBInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitImageWBClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitImageWBClassInfoDefFn (gpointer inClass);

/** Private: Correlate FArrays. */
static ObitFArray* ArrayCorrel (ObitFArray* in1, ObitFArray* in2, 
				ObitErr *err);

/** Private: Image decomposition */
static void doDecomp (ObitFArray* in1, ObitFArray* in2,ObitFArray* in3,
		      ObitFArray* out, olong iterm, olong norder,
		      ofloat iDelta, ofloat A00, 
		      ofloat A10, ofloat A01, ofloat A11, 
		      ofloat A20, ofloat A02, ofloat A21, ofloat A12, ofloat A22);

/** Private: Threaded Image decomposition */
static gpointer ThreadDecomp (gpointer arg);

/** Private: Make Threaded args */
static olong MakeDecompFuncArgs (ObitThread *thread, 
				 ObitFArray *in1,ObitFArray *in2,ObitFArray *in3,
				 ObitFArray *out, olong iterm, olong norder,
				 DecompFuncArg ***ThreadArgs);

/** Private: Delete Threaded args */
static void KillDecompFuncArgs (olong nargs, DecompFuncArg **ThreadArgs);
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitImageWB* newObitImageWB (gchar* name)
{
  ObitImageWB* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageWBClassInit();

  /* allocate/init structure */
  out = ObitMemAlloc0Name(sizeof(ObitImageWB), "ObitImageWB");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitImageWBInit((gpointer)out);

 return out;
} /* end newObitImageWB */

/**
 * Constructor from ObitImage.
 * output will have same underlying file definition as in.
 * Adds planes as needed for spectral planes and 
 * changes frequency axis if norder >0.
 * Initializes class if needed on first call.
 * \param in     ObitImage to copy, MUST NOT have been fully instantiated
 * \param norder Maximum order number of beam, 
 *        0=Normal, 1=spectral index, 2=also curvature
 * \param err    ObitErr for reporting errors.
 * \return the new object.
 */
ObitImageWB* ObitImageWBFromImage (ObitImage* in, olong norder, ObitErr *err)
{
  ObitImageWB* out = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "ObitImageWBFromImage";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitImageIsA(in));

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageWBClassInit();

  /* allocate/init structure */
  out = ObitMemAlloc0Name(sizeof(ObitImageWB), "ObitImageWB");

  /* initialize values */
  if (in->name!=NULL) out->name = g_strdup(in->name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitImageWBInit((gpointer)out);

  /* Copy info for same definition */
  out->info = ObitInfoListCopyData (in->info, out->info);

  /* Copy descriptor */
  out->myDesc = ObitImageDescCopy(in->myDesc, out->myDesc, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /* Actually adding spectral information? */
  if (norder>0) {
    /* Add norder planes */
    out->myDesc->inaxes[out->myDesc->jlocf] = 1 + norder;
    
    /* Change frequency axis labeling */ 
    strncpy (out->myDesc->ctype[out->myDesc->jlocf], "SPECLOGF", IMLEN_KEYWORD-1);
    /* Save reference frequency */
    out->refFreq = out->myDesc->crval[out->myDesc->jlocf];  
    
    /* Add order to infolist */
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut(out->info, "NTERM", OBIT_long, dim, &norder);
  }  /* end Actually adding spectral information */

  /* Beams and things */
  out->myBeam = ObitRef(in->myBeam);
  ObitImageWBSetOrder (out, norder, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  return out;
} /* end ObitImageWBFromImage */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitImageWBGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObiImagetGetWBClass */

/**
 * Test if two ObitImageWBs have the same underlying structures.
 * This test is done using values entered into the #ObitInfoList
 * in case the object has not yet been opened.
 * \param in1 First object to compare
 * \param in2 Second object to compare
 * \param err ObitErr for reporting errors.
 * \return TRUE if to objects have the same underlying structures
 * else FALSE
 */
gboolean ObitImageWBSame (ObitImageWB *in1, ObitImageWB *in2, ObitErr *err )
{
  /* Call ObitData function */
  return ObitDataSame ((ObitData*)in1, (ObitData*)in2, err);
} /* end ObitImageWBSame */

/**
 * Delete underlying files and the basic object.
 * \param inn  Pointer to object to be zapped.
 * \param err  ObitErr for reporting errors.
 * \return pointer for input object, NULL if deletion successful
 */
ObitImage* ObitImageWBZap (ObitImage *inn, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitImage *in = (ObitImage*)inn;
  gchar *routine = "ObitImageWBZap";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return inn;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* Any actual I/O? If not just delete object */
  if (in->mySel->FileType==OBIT_IO_MEM) return ObitImageWBUnref(in);

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
    return ObitImageWBUnref(in); 
  }

  /* Delete Image and all tables  */
  ObitIOZap (in->myIO, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, inn);

  /* If it's scratch remove from list */
  if (in->isScratch) ObitSystemFreeScratch ((Obit*)in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, in);

  /* Delete object */
  in->isScratch = 0; /* Already deleted underlying structures */
  while (in) in = ObitImageWBUnref(in);
  
  return (ObitImage*)in;
} /* end ObitImageWBZap */

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
ObitImageWB* ObitImageWBCopy (ObitImageWB *in, ObitImageWB *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOCode iretCode, oretCode;
  gboolean oldExist;
  ObitHistory *inHist=NULL, *outHist=NULL;
  gchar *outName=NULL, *today=NULL;
  gchar *routine = "ObitImageWBCopy";

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
    out = newObitImageWB(outName);
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
} /* end ObitImageWBCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an Image similar to the input one.
 * Output version set to floating pixels
 * \param in  The object to copy, info may have
 * \param out An existing object pointer for output, info may have
 * \li Size OBIT_int (?,1,1) Dimension of the desired  Image
 * \param err Error stack, returns if not empty.
 */
void ObitImageWBClone  (ObitImageWB *in, ObitImageWB *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOCode iretCode, oretCode;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, size[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ObitHistory *inHist=NULL, *outHist=NULL;
  gchar *today=NULL;
  gchar *exclude[]={"AIPS CC", "AIPS HI", "AIPS PL", "AIPS SL", NULL}; 
  gchar *routine = "ObitImageWBClone";

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
} /* end ObitImageWBClone */

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
void ObitImageWBFullInstantiate (ObitImageWB *in, gboolean exist, ObitErr *err)
{
  ObitIOAccess access;
  gchar *routine = "ObitImageWBFullInstantiate";

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
} /* end ObitImageWBFullInstantiate */

/**
 * Sets maximum order and initializes work arrays
 * \param in     Pointer to object, should be fully defined
 * \param order  Max. beam order to use, 0=> only dirty beam, min 0.
 * \param err    ObitErr for reporting errors.
 */
void ObitImageWBSetOrder (ObitImageWB *in, olong order, 
			  ObitErr *err)
{
  olong i, j, num, naxis[2]; /*ndim=2*/
  /*gchar *routine = "ObitImageWBSetOrder";*/
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  
  in->order = MAX (0, order);
  in->curOrder = in->order;
  num = MAX(1,(in->order+1));
  
  /* Define ResidArr  */
  naxis[0] = in->myDesc->inaxes[0];
  naxis[1] = in->myDesc->inaxes[1];
  in->ResidArr = g_malloc0(num*sizeof(ObitFArray*));
  for (i=0; i<num; i++) {
    /*in->ResidArr[i-1] = ObitFArrayCreate ("ResidArr", ndim, naxis);*/
  }

  /* Define BeamMatx  */
  in->BeamMatx = g_malloc0(num*sizeof(ObitFArray**));
  for (i=0; i<num; i++) {
    in->BeamMatx[i] = g_malloc0(num*sizeof(ObitFArray*));
      for (j=0; j<num; j++) {
	/*in->BeamMatx[i][j] =  ObitFArrayCreate ("BeamMatx", ndim, naxis);*/
      }
  }

} /* end  ObitImageWBSetOrder */

/**
 * Allocates and populates work arrays, 
 * replaces residual image with SW estimate.
 * Normal beam and residual for order = 0 (normal CLEAN)
 * \param in     Pointer to object, should be fully defined and "fresh"
 * \param err    ObitErr for reporting errors.
 */
void ObitImageWBMakeWork (ObitImageWB *in, ObitErr *err)
{
  ObitImage *myBeam = NULL;
  ObitFArray *imPix=NULL, *bmPix[5]={NULL,NULL,NULL,NULL,NULL}; 
  olong i, j, num, lbeam, order, plane[5] = {1,1,1,1,1}, pos[2];
  ofloat normFact;
  gchar *routine = "ObitImageWBMakeWork";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageWBIsA(in));
  myBeam = (ObitImage*)in->myBeam;
  Obit_return_if_fail((ObitImageIsA(myBeam)), err, 
		      "%s: Beam not defined for %s", 
		      routine, in->name);
  /* Make sure it's fresh */
  Obit_return_if_fail((in->fresh), err, 
		      "%s: Image %s is stale", 
		      routine, in->name);

  /* Deallocate any previous */
  ObitImageWBClearWork (in);

  /* Which order using? */
  order = MIN (in->order, in->curOrder);

  /* Read image and beams */
  num = MAX(1, (order+1));
  plane[0] = 1;
  ObitImageGetPlane ((ObitImage*)in, NULL, plane, err);
  imPix = ObitFArrayRef(in->image);  /* Pixel array */
  in->image = ObitFArrayUnref(in->image);
  for (i=0; i<num; i++) {
    plane[0] = i+1;
    ObitImageGetPlane (myBeam, NULL, plane, err);
    bmPix[i] = ObitFArrayRef(myBeam->image);  /* Pixel array */
    myBeam->image = ObitFArrayUnref(myBeam->image);
  }
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Basic output arrays */
  if (!in->ResidArr) in->ResidArr = g_malloc0(num*sizeof(ObitFArray**));
  if (!in->BeamMatx) {
    in->BeamMatx = g_malloc0(num*sizeof(ObitFArray**));
    for (i=0; i<num; i++)
      in->BeamMatx[i] = g_malloc0((num)*sizeof(ObitFArray*));
  }    

  /* By beam order */
  if (order<=0) {           /* Zero order */
    /* No convolutions only beam and image */
    in->BeamMatx[0][0] = ObitFArrayUnref(in->BeamMatx[0][0]);  /* in case */
    in->BeamMatx[0][0] = ObitFArrayCopy (bmPix[0], in->BeamMatx[0][0], err);
    in->ResidArr[0]    = ObitFArrayUnref(in->ResidArr[0]);  /* in case */
    in->ResidArr[0]    = ObitFArrayCopy (imPix, in->ResidArr[0], err);
  } else {                      /* Higher order - do convolutions */
    for (j=0; j<num; j++) {
      in->ResidArr[j] = ObitFArrayUnref(in->ResidArr[j]);  /* in case */
      in->ResidArr[j] = ArrayCorrel (imPix, bmPix[j], err);
    }
    for (j=0; j<num; j++) {
      for (i=0; i<num; i++) {
	in->BeamMatx[j][i] = ObitFArrayUnref(in->BeamMatx[j][i]);  /* in case */
	/* Some pairs of beams are identical - only compute one */
	if (j<=i)
	  in->BeamMatx[j][i] = ArrayCorrel (bmPix[j], bmPix[i], err);
	else
	  in->BeamMatx[j][i] = ObitFArrayRef(in->BeamMatx[i][j]);
      } /* end inner loop */
    } /* end outer loop */
  }

  /* Normalize everything by peak in Beam00 get residuals in approx Jy  */
  lbeam = in->BeamMatx[0][0]->naxis[0];
  pos[0] = lbeam/2; pos[1] = lbeam/2;
  normFact = *ObitFArrayIndex(in->BeamMatx[0][0], pos);
  if (normFact!=0.0) normFact = 1.0 / normFact;
  for (j=0; j<num; j++) {
    ObitFArraySMul (in->ResidArr[j], normFact);
    for (i=0; i<num; i++) {
      /* Some pairs of beams are identical - only compute one */
      if (j<=i)
	ObitFArraySMul (in->BeamMatx[j][i], normFact);
    } 
  }
  
  /* DEBUG
  ObitImageUtilArray2Image ("DbugCorrelResid0.fits", 0, in->ResidArr[0], err); 
  ObitImageUtilArray2Image ("DbugBeam0.fits",  0, bmPix[0], err);
  ObitImageUtilArray2Image ("DbugBeam00.fits", 0, in->BeamMatx[0][0], err);
  if (order>0) {
    ObitImageUtilArray2Image ("DbugCorrelResid1.fits", 0, in->ResidArr[1], err);
    ObitImageUtilArray2Image ("DbugBeam1.fits",  0, bmPix[1], err);
    ObitImageUtilArray2Image ("DbugBeam01.fits", 0, in->BeamMatx[0][1], err);
    ObitImageUtilArray2Image ("DbugBeam10.fits", 0, in->BeamMatx[1][0], err);
    ObitImageUtilArray2Image ("DbugBeam11.fits", 0, in->BeamMatx[1][1], err);
  }
  if (order>1) {
    ObitImageUtilArray2Image ("DbugCorrelResid2.fits", 0, in->ResidArr[2], err);
    ObitImageUtilArray2Image ("DbugBeam1.fits",  0, bmPix[1], err);
    ObitImageUtilArray2Image ("DbugBeam02.fits", 0, in->BeamMatx[0][2], err);
    ObitImageUtilArray2Image ("DbugBeam12.fits", 0, in->BeamMatx[1][2], err);
    ObitImageUtilArray2Image ("DbugBeam22.fits", 0, in->BeamMatx[2][2], err);
    ObitImageUtilArray2Image ("DbugBeam20.fits", 0, in->BeamMatx[2][0], err);
    ObitImageUtilArray2Image ("DbugBeam21.fits", 0, in->BeamMatx[2][1], err);
  }
  if (err->error) Obit_traceback_msg (err, routine, in->name);  */
  /* END DEBUG */

  /* Decompose image into SW spectral components */
  ObitImageWBDecomp (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  imPix = ObitFArrayUnref(imPix);
  for (i=0; i<num; i++) {
    bmPix[i] = ObitFArrayUnref(bmPix[i]);
  }
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitImageWBMakeWork */

/**
 * Decompose image to spectral components
 * Uses SW algorithm to decompose each pixel into flux, si, curvature
 * Note: this is not a deconvolution and replaces the "residual" pixels.
 * Does nothing if order = 0
 * \param in     Pointer to object, should be fully defined and "fresh"
 * \param err    ObitErr for reporting errors.
 */
void ObitImageWBDecomp (ObitImageWB *in, ObitErr *err)
{
  ofloat A00=0.0, A10=0.0, A01=0.0, A11=0.0;
  ofloat A20=0.0, A02=0.0, A21=0.0, A12=0.0, A22=0.0;
  ofloat iDelta;
  olong pos[2];
  ObitFArray *imPix=NULL;
  olong order, plane[5] = {1,1,1,1,1};
  gchar *routine = "ObitImageWBDecomp";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageWBIsA(in));
 
  /* Which order using? */
  order = MIN (in->order, in->curOrder);
  if (order<=0) return;  /* anything to do? */

  /* Are the work arrays populated? */
  if ((in->BeamMatx[0][0]==NULL) || (in->ResidArr[0]==NULL)) 
    ObitImageWBMakeWork (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Beam centers needed */
  pos[0] = in->BeamMatx[0][0]->naxis[0]/2; 
  pos[1] = in->BeamMatx[0][0]->naxis[1]/2; 
  A00 = *ObitFArrayIndex(in->BeamMatx[0][0], pos);
  if (order>0) {  /* First order */
    pos[0] = in->BeamMatx[0][1]->naxis[0]/2; 
    pos[1] = in->BeamMatx[0][1]->naxis[1]/2; 
    A01 = *ObitFArrayIndex(in->BeamMatx[0][1], pos);
    pos[0] = in->BeamMatx[1][0]->naxis[0]/2; 
    pos[1] = in->BeamMatx[1][0]->naxis[1]/2; 
    A10 = *ObitFArrayIndex(in->BeamMatx[1][0], pos);
    pos[0] = in->BeamMatx[1][1]->naxis[0]/2; 
    pos[1] = in->BeamMatx[1][1]->naxis[1]/2; 
    A11 = *ObitFArrayIndex(in->BeamMatx[1][1], pos);
  }
  if (order>1) {  /* second order */
    pos[0] = in->BeamMatx[2][0]->naxis[0]/2; 
    pos[1] = in->BeamMatx[2][0]->naxis[1]/2; 
    A20 = *ObitFArrayIndex(in->BeamMatx[2][0], pos);
    pos[0] = in->BeamMatx[0][2]->naxis[0]/2; 
    pos[1] = in->BeamMatx[0][2]->naxis[1]/2; 
    A02 = *ObitFArrayIndex(in->BeamMatx[0][2], pos);
    pos[0] = in->BeamMatx[2][1]->naxis[0]/2; 
    pos[1] = in->BeamMatx[2][1]->naxis[1]/2; 
    A21 = *ObitFArrayIndex(in->BeamMatx[2][1], pos);
    pos[0] = in->BeamMatx[1][2]->naxis[0]/2; 
    pos[1] = in->BeamMatx[1][2]->naxis[1]/2; 
    A12 = *ObitFArrayIndex(in->BeamMatx[1][2], pos);
    pos[0] = in->BeamMatx[2][2]->naxis[0]/2; 
    pos[1] = in->BeamMatx[2][2]->naxis[1]/2; 
    A22 = *ObitFArrayIndex(in->BeamMatx[2][2], pos);
  }

  /* Zero order = flux - read old */
  plane[0] = 1;
  ObitImageGetPlane ((ObitImage*)in, NULL, plane, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  imPix = ObitFArrayCopy(in->image, imPix, err);   /* Pixel array */
  in->extBuffer = TRUE;  /* Use separate buffer */

  if (order==1) { /* Zero order = flux  order = 1 */
    /* Check compatabilities */
    Obit_return_if_fail(ObitFArrayIsCompatable(imPix, in->ResidArr[0]), err,
			"%s: output and resid 0 FArrays  incompatible", 
			routine);
    Obit_return_if_fail(ObitFArrayIsCompatable(imPix, in->ResidArr[1]), err,
			"%s: output and resid 1 FArrays  incompatible", 
			routine);
    iDelta = 1.0 / (A00*A11 - A01*A10);  /* 1/Determinant */
    /* Decompose */
    doDecomp (in->ResidArr[0], in->ResidArr[1], in->ResidArr[2], imPix,
	      0, order, 
	      iDelta,  A00, A10, A01, A11, A20, A02, A21,  A12,  A22);
  } else {  /* Zero order = flux  norder = 2 */
  /* Check compatabilities */
    Obit_return_if_fail(ObitFArrayIsCompatable(imPix, in->ResidArr[0]), err,
			"%s: output and resid 0 FArrays  incompatible", 
			routine);
    Obit_return_if_fail(ObitFArrayIsCompatable(imPix, in->ResidArr[1]), err,
			"%s: output and resid 1 FArrays  incompatible", 
			routine);
    Obit_return_if_fail(ObitFArrayIsCompatable(imPix, in->ResidArr[2]), err,
			"%s: output and resid 2 FArrays  incompatible", 
			routine);
    /* 1/Determinant */
    iDelta = 1.0 / (A02*A02*A11 - 2.0*A01*A02*A12 + A01*A01*A22 + A00*(A12*A12-A11*A22));

    /* Decompose */
    doDecomp (in->ResidArr[0], in->ResidArr[1], in->ResidArr[2], imPix,
	      0, order, 
	      iDelta,  A00, A10, A01, A11, A20, A02, A21,  A12,  A22);
  } /* end flux */
  /* Write flux plane */
  plane[0] = 1;
  ObitImagePutPlane ((ObitImage*)in, imPix->array, plane, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* First order = spectral index */
  if (order==1) { /* First order =  spectral index  norder = 1 */
    iDelta = 1.0 / (A00*A11 - A01*A10);  /* 1/Determinant */
    /* Decompose */
    doDecomp (in->ResidArr[0], in->ResidArr[1], in->ResidArr[2], imPix,
	      1, order, 
	      iDelta,  A00, A10, A01, A11, A20, A02, A21,  A12,  A22);
  } else {  /* Zero order = si  norder = 2 */
    /* -1/Determinant */
    iDelta = 1.0 / (A02*A02*A11 - 2.0*A01*A02*A12 + A01*A01*A22 + A00*(A12*A12-A11*A22));

    /* Decompose */
    doDecomp (in->ResidArr[0], in->ResidArr[1], in->ResidArr[2], imPix,
	      1, order, 
	      iDelta,  A00, A10, A01, A11, A20, A02, A21,  A12,  A22);
  } /* end si */
  /* Write si plane */
  plane[0] = 2;
  ObitImagePutPlane ((ObitImage*)in, imPix->array, plane, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Second order = curvature */
  if (order>1) { /*     */
    /* -1/Determinant */
    iDelta = 1.0 / (A02*A02*A11 - 2.0*A01*A02*A12 + A01*A01*A22 + A00*(A12*A12-A11*A22)); /* OK */

    /* Decompose */
    doDecomp (in->ResidArr[0], in->ResidArr[1], in->ResidArr[2], imPix,
	      2, order, 
	      iDelta,  A00, A10, A01, A11, A20, A02, A21,  A12,  A22);
    /* Write curvature plane */
    plane[0] = 3;
    ObitImagePutPlane ((ObitImage*)in, imPix->array, plane, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end curvature */
  
  in->extBuffer = FALSE;  /* Reset separate buffer */

  /* It's now stale */
  in->fresh = FALSE;

  /* Cleanup */
  /* Ho No ObitImageWBClearWork (in);*/
  ObitFArrayUnref (imPix);
  
} /* end ObitImageWBDecomp */

/**
 * Deallocates work arrays
 * \param in     Pointer to object
 */
void ObitImageWBClearWork (ObitImageWB *in)
{
  olong i, j, num;

  g_assert (ObitIsA(in, &myClassInfo));

  /* If there is anything to delete, the image goes stale */
  
  num = MAX(1,(in->order+1));
  for (i=0; i<num; i++) {
    in->fresh = in->fresh && (in->ResidArr[i]==NULL); /* stale? */
    in->ResidArr[i] = ObitFArrayUnref(in->ResidArr[i]);
    for (j=0; j<num; j++) {
      in->fresh = in->fresh && (in->BeamMatx[i][j]==NULL);
      in->BeamMatx[i][j] = ObitFArrayUnref(in->BeamMatx[i][j]); /* stale? */
    }
  }

} /* end  ObitImageWBClearWork */


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
ObitImage* ObitImageWBGetBeam (ObitImage *inn, olong beamNo, 
			       olong plane[5], ObitErr *err) 
{
  ObitImageWB *in = (ObitImageWB*)inn;
  ObitImage *theBeam;
  olong i;
  gchar *routine = "ObitImageWBGetBeam";

  /* error checks */
  g_assert (ObitImageWBIsA(inn));

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
} /* end ObitImageWBGetBeam */

/**
 * Get current highest order of image beam
 * Used for Sault-Wieringa wideband imaging
 * \param inn  Image whose beam name is to be set 
 * \return order number
 */
olong ObitImageWBGetBeamOrder (ObitImage *inn) 
{
  ObitImageWB *in = (ObitImageWB*)inn;

  /* error checks */
  g_assert (ObitImageIsA(inn));

  return MIN (in->order, in->curOrder);
} /* end ObitImageWBGetBeam */

/*-------Private functions called by ObitData class ------*/
/** Private: Zap */
static ObitData* ObitDataImageWBZap (ObitData *in, ObitErr *err)
{
  return (ObitData*)ObitImageWBZap ((ObitImage*)in, err);
} /* end ObitDataImageWBZap */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitImageWBClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();
  myClassInfo.hasScratch    = TRUE; /* Scratch files allowed */

  /* Set function pointers */
  ObitImageWBClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitImageWBClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitImageWBClassInfoDefFn ( gpointer inClass)
{
  ObitImageWBClassInfo *theClass = (ObitImageWBClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitImageWBClassInit;
  theClass->newObit       = (newObitFP)newObitImageWB;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitImageWBClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitImageWBGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitImageWBCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitImageWBClear;
  theClass->ObitInit      = (ObitInitFP)ObitImageWBInit;
  theClass->ObitImageGetBeam = 
    (ObitImageGetBeamFP)ObitImageWBGetBeam;
  theClass->ObitImageGetBeamOrder = 
    (ObitImageGetBeamOrderFP)ObitImageWBGetBeamOrder;

  /* Function pointers referenced from ObitData class */
  theClass->ObitDataZap     = (ObitDataZapFP)ObitDataImageWBZap;

} /* end ObitImageWBClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitImageWBInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageWB *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->BeamMatx   = NULL;
  in->ResidArr   = NULL;
  in->fresh      = FALSE;
  in->order      = 0;
  in->curOrder   = 0;

} /* end ObitImageWBInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitImageWB* cast to an Obit*.
 */
void ObitImageWBClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageWB *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Delete underlying higher order beam files if isScratch */
  /* delete this class members */
  ObitImageWBClearWork (in);
  for (i=0; i<in->order; i++) {
    g_free(in->BeamMatx[i]);
  }
  g_free(in->BeamMatx);
  g_free(in->ResidArr);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitImageWBClear */

/**
 * Convolves two arrays
 * The two are at least doubles in size with zero padding with the 
 * smaller (if any) promoted to the size of the larger.
 * Output is trimmed to the size of the first input.
 * \param in1   First input 2D array
 * \param in2   Second input 2D array
 * \param err   ObitErr for reporting errors.
 * \return convolved 2D array, must be Unreffed when no longer needed.
 */
static ObitFArray* ArrayCorrel (ObitFArray* in1, ObitFArray* in2, 
			       ObitErr *err)
{
  ObitFArray *out = NULL;
  ObitFArray *pad1 = NULL, *pad2 = NULL, *convl = NULL;
  olong nx, ny, naxis[2], cen[2], blc[2], trc[2];
  ofloat factor;
  gchar *routine = "ObitImageWB:ArrayCorrel";

  /* error checks */
  if (err->error) return out;

  /* Get sizes of larger  */
  nx = MAX (in1->naxis[0], in2->naxis[0]);
  ny = MAX (in1->naxis[1], in2->naxis[1]);

  /* Now, what is needed for double size FFT? */
  naxis[0] = ObitFFTSuggestSize(2*nx);
  naxis[1] = ObitFFTSuggestSize(2*ny);

  /* Pad images */
  factor = 1.0;
  pad1 = ObitFArrayCreate ("Array1", 2, naxis);
  ObitFArrayPad (in1, pad1, factor);
  pad2 = ObitFArrayCreate ("Array2", 2, naxis);
  ObitFArrayPad (in2, pad2, factor);

  /* Correlate */
  convl = ObitFArrayUtilCorrel (pad1, pad2, err);
  if (err->error) goto cleanup;

  /* Extract output  */
  cen[0] = convl->naxis[0]/2; cen[1] = convl->naxis[1]/2;
  blc[0] = cen[0] - nx/2; blc[1] = cen[1] - ny/2;
  trc[0] = cen[0] + nx/2 - 1; trc[1] = cen[1] + ny/2 - 1;
  out = ObitFArraySubArr (convl, blc, trc, err);
  if (err->error) goto cleanup;

 cleanup:
  ObitFArrayUnref(pad1);
  ObitFArrayUnref(pad2);
  ObitFArrayUnref(convl);
  
  if (err->error) Obit_traceback_val (err, routine, in1->name, out);
  return out;
} /* end ArrayCorrel */

/**
 * Decompose FArrays into spectral terms possible using threads
 * \param in1     First ObitFArray, residual * beam0 
 * \param in2     Second ObitFArray residual * beam1
 * \param in3     Third ObitFArray residual * beam2
 * \param out     ObitFArray for decomposition 
 * \param iterm   which spectral term estimated?
 * \param norder  which spectral order?
 * \param iDelta  1/determanant
 * \param A00     Central beam00 value
 * \param A10     Central beam10 value
 * \param A01     Central beam01 value
 * \param A11     Central beam11 value
 * \param A20     Central beam20 value
 * \param A02     Central beam02 value
 * \param A21     Central beam21 value
 * \param A12     Central beam12 value
 * \param A22     Central beam22 value
 */
static void doDecomp (ObitFArray* in1, ObitFArray* in2, ObitFArray* in3,
		      ObitFArray* out, olong iterm, olong norder,
		      ofloat iDelta, ofloat A00, 
		      ofloat A10, ofloat A01, ofloat A11, 
		      ofloat A20, ofloat A02, ofloat A21, ofloat A12, ofloat A22)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  gboolean OK;
  DecompFuncArg **threadArgs;

  /* Initialize Threading */
  nThreads = MakeDecompFuncArgs (in1->thread, in1, in2, in3, out, iterm, norder,
				 &threadArgs);
  
  /* Divide up work */
  nElem = in1->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
  if (nElem<1000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem;
    threadArgs[i]->last    = hiElem;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    threadArgs[i]->iterm  = iterm;
    threadArgs[i]->norder = norder;
    threadArgs[i]->iDelta = iDelta;
    threadArgs[i]->A00    = A00;
    threadArgs[i]->A10    = A10;
    threadArgs[i]->A01    = A01;
    threadArgs[i]->A11    = A11;
    threadArgs[i]->A20    = A20;
    threadArgs[i]->A02    = A02;
    threadArgs[i]->A21    = A21;
    threadArgs[i]->A12    = A12;
    threadArgs[i]->A22    = A22;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  OK = ObitThreadIterator (in1->thread, nTh, 
			   (ObitThreadFunc)ThreadDecomp,
			   (gpointer**)threadArgs);

  /* Free local objects */
  KillDecompFuncArgs(nThreads, threadArgs);
} /* end doDecomp */


/**
 * Image decomposition function
 * Callable as thread
 * \param arg Pointer to DecompFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li in1      First ObitFArray, residual * beam0 
 * \li in2      Second ObitFArray residual * beam1
 * \li in3      Third ObitFArray residual * beam2
 * \li out      ObitFArray for decomposition 
 * \li iterm    which spectral term estimated?
 * \li norder   which spectral order?
 * \li iDelta   1/determanant
 * \li A00      Central beam value
 * \li A10      Central beam value
 * \li A01      Central beam value
 * \li A11      Central beam value
 * \li A20      Central beam value
 * \li A02      Central beam value
 * \li A21      Central beam value
 * \li A12      Central beam value
 * \li A22      Central beam value
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadDecomp (gpointer arg)
{
  /* Get arguments from structure */
  DecompFuncArg *largs = (DecompFuncArg*)arg;
  ObitFArray *in1       = largs->in1;
  ObitFArray *in2       = largs->in2;
  ObitFArray *in3       = largs->in3;
  ObitFArray *out       = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last-1;
  olong      iterm      = largs->iterm;
  olong      norder     = largs->norder;
  ofloat     iDelta     = largs->iDelta;
  ofloat     A00        = largs->A00;
  ofloat     A10        = largs->A10;
  ofloat     A01        = largs->A01;
  ofloat     A11        = largs->A11;
  ofloat     A02        = largs->A02;
  ofloat     A12        = largs->A12;
  ofloat     A22        = largs->A22;

  /* local */
  olong  i;
  ofloat R0, R1, R2;

  if (hiElem<loElem) goto finish;

  /* Only first order */
  if (norder==1) {
    if (iterm==0) {         /* flux */
      /* Loop over my elements */
      for (i=loElem; i<hiElem; i++) {
	R0 = in1->array[i];
	R1 = in2->array[i];
	out->array[i] = (A11*R0  - A01*R1) * iDelta;
      } /* end loop computing flux */
    } else if (iterm==1) {  /* Spectral index */
      /* Loop over my elements */
      for (i=loElem; i<hiElem; i++) {
	R0 = in1->array[i];
	R1 = in2->array[i];
	out->array[i] = (A00*R1  - A10*R0) * iDelta;
      } /* end loop computing spectral index */
    }
    
    /* Second order */
  } else if (norder==2) {
    if (iterm==0) {         /* flux */
      /* Loop over my elements */
      for (i=loElem; i<hiElem; i++) {
	R0 = in1->array[i];
	R1 = in2->array[i];
	R2 = in3->array[i];
	out->array[i] = (A12*A12*R0 + A01*A22*R1 - A12*( A02*R1 + A01*R2) + A11*(-A22*R0 + A02*R2)) * iDelta;
      } /* end loop computing flux */
    } else if (iterm==1) {  /* Spectral index */
      /* Loop over my elements */
      for (i=loElem; i<hiElem; i++) {
	R0 = in1->array[i];
	R1 = in2->array[i];
	R2 = in3->array[i];
	out->array[i] = (  A01*A22*R0 + A02*A02*R1 - A02*( A12*R0 + A01*R2) + A00*(-A22*R1 + A12*R2)) * iDelta;    
      } /* end loop computing spectral index */
    } else if (iterm==2) {  /* Curvature */
      /* Loop over my elements */
      for (i=loElem; i<hiElem; i++) {
	R0 = in1->array[i];
	R1 = in2->array[i];
	R2 = in3->array[i];
	out->array[i] =  (-A01*A12*R0 + A02*(A11*R0 - A01*R1) + A01*A01*R2 + A00*(A12*R1 - A11*R2)) * iDelta;
      } /* end loop computing spectral Curvature */
    }
  } /* End 2nd order */
  
  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /*  end ThreadDecomp */

/**
 * Make arguments for a Threaded ThreadDecompFunc?
 * \param thread     ObitThread object to be used
 * \param in         FA to be operated on
 * \param out        FA to the written
 * \param iterm      Which spectral term estimated?
 * \param norder     Which spectral order?
 * \param ThreadArgs[out] Created array of DecompFuncArg, 
 *                   delete with KillDecompFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeDecompFuncArgs (ObitThread *thread, 
				 ObitFArray *in1, ObitFArray *in2, ObitFArray *in3,
				 ObitFArray *out, olong iterm, olong norder,
				 DecompFuncArg ***ThreadArgs)

{
  olong i, nThreads;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(DecompFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(DecompFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->thread= ObitThreadRef(thread);
    (*ThreadArgs)[i]->in1   = ObitFArrayRef(in1);
    (*ThreadArgs)[i]->in2   = ObitFArrayRef(in2);
    if (in3) (*ThreadArgs)[i]->in3   = ObitFArrayRef(in3);
    else (*ThreadArgs)[i]->in3       = NULL;
    (*ThreadArgs)[i]->out   = ObitFArrayRef(out);
    (*ThreadArgs)[i]->first = 1;
    (*ThreadArgs)[i]->last  = in1->arraySize;
    (*ThreadArgs)[i]->iterm = iterm;
    (*ThreadArgs)[i]->norder= norder;
    (*ThreadArgs)[i]->ithread = i;
  }

  return nThreads;
} /*  end MakeInterpImageArgs */

/**
 * Delete arguments for ThreadDecompFunc
 * \param nargs      number of elements in ThreadArgs.
 * \param ThreadArgs Array of DecompFuncArg
 */
static void KillDecompFuncArgs (olong nargs, DecompFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      if (ThreadArgs[i]->thread) ObitThreadUnref(ThreadArgs[i]->thread);
      if (ThreadArgs[i]->in1)    ObitFArrayUnref(ThreadArgs[i]->in1);
      if (ThreadArgs[i]->in2)    ObitFArrayUnref(ThreadArgs[i]->in2);
      if (ThreadArgs[i]->in3)    ObitFArrayUnref(ThreadArgs[i]->in3);
      if (ThreadArgs[i]->out)    ObitFArrayUnref(ThreadArgs[i]->out);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillDecompFuncArgs */

