/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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

#include "ObitRMFit.h"
#include "ObitThread.h"
#ifdef HAVE_GSL
#include <gsl/gsl_fit.h>
#endif /* HAVE_GSL */ 
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif /* VELIGHT */
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitRMFit.c
 * ObitRMFit class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitRMFit";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitRMFitClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitRMFitClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitRMFitInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitRMFitClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitRMFitClassInfoDefFn (gpointer inClass);

/** Private: Do actual fitting */
static void Fitter (ObitRMFit* in, ObitErr *err);

/** Private: Write output image */
static void WriteOutput (ObitRMFit* in, ObitImage *outImage, 
			 ObitErr *err);


/** Private: Threaded fitting */
static gpointer ThreadNLRMFit (gpointer arg);

/*---------------Private structures----------------*/
/* FT threaded function argument */
typedef struct {
  /** ObitRMFit object */
  ObitRMFit *in;
  /* Obit error stack object */
  ObitErr      *err;
  /** First (1-rel) value in y to process this thread */
  olong        first;
  /** Highest (1-rel) value in y to process this thread  */
  olong        last;
  /** thread number, >0 -> no threading   */
  olong        ithread;
  /** max number of terms to fit  */
  olong        nterm;
  /** number of frequencies  */
  olong        nlamb2;
  /** Array of Lambda^2 (nlamb2) */
  ofloat *lamb2;
  /** Array of Q, U weights (1/RMS**2) per inFArrays (nlamb2) */
  ofloat *Qweight, *Uweight;
  /** Array of Q, U variance per inFArrays (nlamb2) */
  ofloat *Qvar, *Uvar;
  /** Array of Q, U pixel values being fitted (nlamb2) */
  ofloat *Qobs, *Uobs;
  /** min Q/U SNR */
  ofloat minQUSNR;
  /** Reference lambda^2 */
  ofloat refLamb2;
  /** Vector of guess/fitted coefficients, optional errors */
  ofloat *coef;
  /** Chi squared of fit */
  ofloat ChiSq;
  /** Do error analysis? */
  gboolean doError;
  /** work arrays */
  double *x, *y, *w;
} NLRMFitArg;

/** Private: Actual fitting */
static void NLRMFit (NLRMFitArg *arg);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitRMFit* newObitRMFit (gchar* name)
{
  ObitRMFit* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitRMFitClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitRMFit));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitRMFitInit((gpointer)out);

 return out;
} /* end newObitRMFit */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitRMFitGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitRMFitClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitRMFitGetClass */

/**
 * Make a deep copy of an ObitRMFit.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitRMFit* ObitRMFitCopy  (ObitRMFit *in, ObitRMFit *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, nOut;

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitRMFit(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nlamb2 = in->nlamb2;
  out->nterm = in->nterm;

  /* Arrays */
  if (out->QRMS) g_free(out->QRMS);
  out->QRMS = g_malloc0(out->nlamb2*sizeof(ofloat));
  if (in->QRMS)
    for (i=0; i<out->nlamb2; i++) out->QRMS[i] = in->QRMS[i];
  if (out->URMS) g_free(out->URMS);
  out->URMS = g_malloc0(out->nlamb2*sizeof(ofloat));
  if (in->URMS)
    for (i=0; i<out->nlamb2; i++) out->URMS[i] = in->URMS[i];
    
  /* reference this class members */
  if (out->outDesc) out->outDesc = ObitImageDescUnref(out->outDesc);
  if (in->outDesc)  out->outDesc = ObitImageDescRef(in->outDesc);
  if (out->inQFArrays) {
    for (i=0; i<in->nlamb2; i++) out->inQFArrays[i] = ObitFArrayUnref(out->inQFArrays[i]);
  }
  if (in->inQFArrays) {
    for (i=0; i<in->nlamb2; i++) out->inQFArrays[i] = ObitFArrayRef(in->inQFArrays[i]);
  }
  if (out->inUFArrays) {
    for (i=0; i<in->nlamb2; i++) out->inUFArrays[i] = ObitFArrayUnref(out->inUFArrays[i]);
  }
  if (in->inUFArrays) {
    for (i=0; i<in->nlamb2; i++) out->inUFArrays[i] = ObitFArrayRef(in->inUFArrays[i]);
  }


  /* How many output planes */
  if (in->doError) nOut = 1+in->nterm*2;
  else nOut = in->nterm;
  if (out->outFArrays) {
    for (i=0; i<nOut; i++) out->outFArrays[i] = ObitFArrayUnref(out->outFArrays[i]);
  }
  if (in->outFArrays) {
    for (i=0; i<nOut; i++) out->outFArrays[i] = ObitFArrayRef(in->outFArrays[i]);
  }

  return out;
} /* end ObitRMFitCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an RMFit similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitRMFitClone  (ObitRMFit *in, ObitRMFit *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i, nOut;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nlamb2 = in->nlamb2;
  out->nterm = in->nterm;
  
  /* Arrays */
  if (out->QRMS) g_free(out->QRMS);
  out->QRMS = g_malloc0(out->nlamb2*sizeof(ofloat));
  if (in->QRMS)
    for (i=0; i<out->nlamb2; i++) out->QRMS[i] = in->QRMS[i];
  if (out->URMS) g_free(out->URMS);
  out->URMS = g_malloc0(out->nlamb2*sizeof(ofloat));
  if (in->URMS)
    for (i=0; i<out->nlamb2; i++) out->URMS[i] = in->URMS[i];

  /* reference this class members */
  if (out->outDesc) out->outDesc = ObitImageDescUnref(out->outDesc);
  if (in->outDesc)  out->outDesc = ObitImageDescRef(in->outDesc);
  if (out->inQFArrays) {
    for (i=0; i<in->nlamb2; i++) out->inQFArrays[i] = ObitFArrayUnref(out->inQFArrays[i]);
  }
  if (in->inQFArrays) {
    for (i=0; i<in->nlamb2; i++) out->inQFArrays[i] = ObitFArrayRef(in->inQFArrays[i]);
  }
  if (out->inUFArrays) {
    for (i=0; i<in->nlamb2; i++) out->inUFArrays[i] = ObitFArrayUnref(out->inUFArrays[i]);
  }
  if (in->inUFArrays) {
    for (i=0; i<in->nlamb2; i++) out->inUFArrays[i] = ObitFArrayRef(in->inUFArrays[i]);
  }

  /* How many output planes */
  if (in->doError) nOut = 1+in->nterm*2;
  else nOut = in->nterm;
  if (out->outFArrays) {
    for (i=0; i<nOut; i++) out->outFArrays[i] = ObitFArrayUnref(out->outFArrays[i]);
  }
  if (in->outFArrays) {
    for (i=0; i<nOut; i++) out->outFArrays[i] = ObitFArrayRef(in->outFArrays[i]);
  }

} /* end ObitRMFitClone */

/**
 * Creates an ObitRMFit 
 * \param name   An optional name for the object.
 * \param nterm  Number of coefficients of powers of log(nu)
 * \return the new object.
 */
ObitRMFit* ObitRMFitCreate (gchar* name, olong nterm)
{
  ObitRMFit* out;

  /* Create basic structure */
  out = newObitRMFit (name);
  out->nterm = nterm;

  return out;
} /* end ObitRMFitCreate */

/**
 * Fit RM to an image cube pixels.
 * The third axis of the output image will be "SPECRM  " to indicate that the
 * planes are RM fit parameters.  The "CRVAL" on this axis will be the reference 
 * Frequency for the fitting.
 * Item "NTERM" is added to the output image descriptor
 * \param in       Spectral fitting object
 *                 Potential parameters on in->info:
 * \li "refLamb2" OBIT_double scalar Reference frequency for fit [def ref for inQImage]
 * \li "minQUSNR" OBIT_float  scalar min. SNR for Q and U pixels [def 3.0]
 * \li "doError" OBIT_boolean scalar If true do error analysis [def False]
 *
 * \param inQImage Q Image cube to be fitted
 * \param inUImage U Image cube to be fitted
 * \param outImage Image cube with fitted spectra.
 *                 Should be defined but not created.
 *                 Planes 1->nterm are coefficients per pixel
 *                 if doError:
 *                 Planes nterm+1->2*nterm are uncertainties in coefficients
 *                 Plane 2*nterm+1 = Chi squared of fit
 * \param err      Obit error stack object.
 */
void ObitRMFitCube (ObitRMFit* in, ObitImage *inQImage, ObitImage *inUImage, 
		    ObitImage *outImage, ObitErr *err)
{
  olong i, iplane, nOut, plane[5]={1,1,1,1,1}, noffset=0;
  olong naxis[2];
  ObitIOSize IOBy;
  ObitInfoType type;
  ObitIOCode retCode;
  union ObitInfoListEquiv InfoReal; 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  odouble freq;
  gchar *today=NULL, *SPECRM   = "SPECRM  ", keyword[9];
  gchar *routine = "ObitRMFitCube";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert(ObitImageIsA(inQImage));
  g_assert(ObitImageIsA(inUImage));
  g_assert(ObitImageIsA(outImage));

  /* Warn if no GSL implementation */
#ifndef HAVE_GSL
  Obit_log_error(err, OBIT_InfoWarn, "NO GSL available - results will be zeroes");
#endif

  /* Control parameters */
  /* Min Q/U pixel SNR for fit */
  InfoReal.flt = 3.0; type = OBIT_float;
  in->minQUSNR = 3.0;
  ObitInfoListGetTest(in->info, "minQUSNR", &type, dim, &InfoReal);
  if (type==OBIT_float) in->minQUSNR = InfoReal.flt;
  else if (type==OBIT_double) in->minQUSNR = (ofloat)InfoReal.dbl;

  /* Want Error analysis? */
  InfoReal.itg = (olong)FALSE; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doError", &type, dim, &InfoReal);
  in->doError = InfoReal.itg;

  /* Open input images to get info */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (inQImage->info, "IOBy", OBIT_long, dim, &IOBy);
  inQImage->extBuffer = TRUE;   /* Using inFArrays as I/O buffer */
  retCode = ObitImageOpen (inQImage, OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inQImage->name);
  ObitInfoListAlwaysPut (inUImage->info, "IOBy", OBIT_long, dim, &IOBy);
  inUImage->extBuffer = TRUE;   /* Using inFArrays as I/O buffer */
  retCode = ObitImageOpen (inUImage, OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inUImage->name);

  /* Check compatability */
  Obit_return_if_fail(((inQImage->myDesc->inaxes[0]==inUImage->myDesc->inaxes[0]) && 
		       (inQImage->myDesc->inaxes[1]==inUImage->myDesc->inaxes[1]) &&
		       (inQImage->myDesc->inaxes[2]==inUImage->myDesc->inaxes[2])), err,
		      "%s: Input images incompatable", routine);

  /* Get Reference frequency/lambda^2 , default from input Q ref. freq. */
  InfoReal.dbl = VELIGHT/inQImage->myDesc->crval[inQImage->myDesc->jlocf]; 
  InfoReal.dbl *= InfoReal.dbl;
  type = OBIT_double;
  in->refLamb2 = InfoReal.dbl;
  ObitInfoListGetTest(in->info, "refLamb2", &type, dim, &InfoReal);
  if (type==OBIT_float) in->refLamb2 = InfoReal.flt;
  else if (type==OBIT_double) in->refLamb2 = (ofloat)InfoReal.dbl;
  
  /* What plane does spectral data start on */
  if (!strncmp(inQImage->myDesc->ctype[inQImage->myDesc->jlocf], "SPECLNMF", 8)) {
    noffset = 2;  /* What plane does spectral data start on */
    ObitInfoListGetTest (inQImage->myDesc->info, "NTERM", &type, dim, &noffset);
  } else {   /* Normal spectral cube */
    noffset = 0;
  }

  /* Determine number of frequency planes and initialize in */
  in->nlamb2     = inQImage->myDesc->inaxes[inQImage->myDesc->jlocf]-noffset;
  in->QRMS       = g_malloc0(in->nlamb2*sizeof(ofloat));
  in->URMS       = g_malloc0(in->nlamb2*sizeof(ofloat));
  in->inQFArrays = g_malloc0(in->nlamb2*sizeof(ObitFArray*));
  in->inUFArrays = g_malloc0(in->nlamb2*sizeof(ObitFArray*));
  in->lamb2      = g_malloc0(in->nlamb2*sizeof(odouble));

  /* How many output planes? */
  if (in->doError) nOut = 1+in->nterm*2;
  else nOut = in->nterm;

   /* Define term arrays */
  in->outFArrays = g_malloc0((nOut)*sizeof(ObitFArray*));
  for (i=0; i<nOut; i++) in->outFArrays[i] = NULL;

 /* Image size */
  in->nx = inQImage->myDesc->inaxes[0];
  in->ny = inQImage->myDesc->inaxes[1];
  naxis[0] = (olong)in->nx;  naxis[1] = (olong)in->ny; 
  for (i=0; i<in->nlamb2; i++) in->inQFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);
  for (i=0; i<in->nlamb2; i++) in->inUFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);
  for (i=0; i<nOut; i++)       in->outFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);

  /* Output Image descriptor */
  outImage->myDesc = ObitImageDescCopy (inQImage->myDesc, in->outDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inQImage->name);

  /* Change third axis to type "SPECRM  " and leave the reference frequency
     as the "CRVAL" */
  outImage->myDesc->inaxes[outImage->myDesc->jlocf] =  nOut;
  outImage->myDesc->crval[outImage->myDesc->jlocf]  =  VELIGHT/sqrt(in->refLamb2);
  outImage->myDesc->crpix[outImage->myDesc->jlocf]  =  1.0;
  outImage->myDesc->cdelt[outImage->myDesc->jlocf]  =  1.0;
  strncpy (outImage->myDesc->ctype[outImage->myDesc->jlocf], SPECRM , IMLEN_KEYWORD);
  outImage->myDesc->bitpix = -32;  /* Float it */

  /* Creation date today */
  today = ObitToday();
  strncpy (outImage->myDesc->date, today, IMLEN_VALUE);
  if (today) g_free(today);

  /* Copy of output descriptor to in */
  in->outDesc = ObitImageDescCopy (outImage->myDesc, in->outDesc, err);

  /* Loop reading planes */
  for (iplane=0; iplane<in->nlamb2; iplane++) {
    /* Lambda^2 Check for MFImage outputs */
    if (!strncmp(inQImage->myDesc->ctype[inQImage->myDesc->jlocf], "SPECLNMF", 8)) {
	sprintf (keyword, "FREQ%4.4d",iplane+1);
	ObitInfoListGetTest (inQImage->myDesc->info, keyword, &type, dim, &freq);
       } else {   /* Normal spectral cube */
	freq = inQImage->myDesc->crval[inQImage->myDesc->jlocf] + 
	  inQImage->myDesc->cdelt[inQImage->myDesc->jlocf] * 
	  (inQImage->myDesc->plane - inQImage->myDesc->crpix[inQImage->myDesc->jlocf]);
      }
    in->lamb2[iplane] = (VELIGHT/freq)*(VELIGHT/freq);

    plane[0] = iplane+noffset;  /* Select correct plane */
    retCode = ObitImageGetPlane (inQImage, in->inQFArrays[iplane]->array, plane, err);
    retCode = ObitImageGetPlane (inUImage, in->inUFArrays[iplane]->array, plane, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) Obit_traceback_msg (err, routine, inQImage->name);

    /* Plane RMSes */
    in->QRMS[iplane] = ObitFArrayRMS(in->inQFArrays[iplane]);
    in->URMS[iplane] = ObitFArrayRMS(in->inUFArrays[iplane]);

    
  } /* end loop reading planes */

  /* Close inputs */
  retCode = ObitImageClose (inQImage, err);
  inQImage->extBuffer = FALSE;   /* May need I/O buffer later */
  retCode = ObitImageClose (inUImage, err);
  inUImage->extBuffer = FALSE;   /* May need I/O buffer later */
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inQImage->name);

  /* Do actual fitting */
  Fitter (in, err);
  if (err->error) Obit_traceback_msg (err, routine, inQImage->name);

  /* Update output header to reference Frequency */
  outImage->myDesc->crval[outImage->myDesc->jlocf] = VELIGHT/sqrt(in->refLamb2);
  in->outDesc->crval[in->outDesc->jlocf] = VELIGHT/sqrt(in->refLamb2);

  /* Write output */
  WriteOutput(in, outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inQImage->name);

} /* end ObitRMFitCube */

/**
 * Fit RM to an array of images
 * The third axis of the output image will be "SPECRM  " to indicate that then
 * planes are spectral fit parameters.  The "CRVAL" on this axis will be the reference 
 * Frequency for the fitting.
 * Item "NTERM" is added to the output image descriptor to give the maximum number 
 * of terms fitted
 * \param in       Spectral fitting object
 * \li "refLamb2" OBIT_double scalar Reference frequency for fit [def average of inputs]
 * \li "minQUSNR" OBIT_float scalar Max. min. SNR for Q and U pixels [def 3.0]
 * \li "doError" OBIT_boolean scalar If true do error analysis [def False]
 * \param nimage   Number of entries in imQArr
 * \param imQArr   Array of Q images to be fitted
 * \param imUArr   Array of U images to be fitted, same geometry as imQArr
 * \param outImage Image cube with fitted parameters.
 *                 Should be defined but not created.
 *                 Planes 1->nterm are coefficients per pixel
 *                 if doError:
 *                 Planes nterm+1->2*nterm are uncertainties in coefficients
 *                 Plane 2*nterm+1 = Chi squared of fit
 * \param err      Obit error stack object.
 */
void ObitRMFitImArr (ObitRMFit* in, olong nimage, 
		     ObitImage **imQArr, ObitImage **imUArr, 
		     ObitImage *outImage, ObitErr *err)
{
  olong i, iplane, nOut;
  olong naxis[2];
  ObitIOSize IOBy;
  ObitInfoType type;
  ObitIOCode retCode;
  union ObitInfoListEquiv InfoReal; 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat ipixel[2], opixel[2];
  gboolean bad;
  odouble avgLamb2, freq;
  gchar *today=NULL, *SPECRM   = "SPECRM  ";
  gchar *routine = "ObitRMFitCube";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert(ObitImageIsA(outImage));

  /* Warn if no GSL implementation */
#ifndef HAVE_GSL
  Obit_log_error(err, OBIT_InfoWarn, "NO GSL available - results will be zeroes");
#endif

  /* Control parameters */
  /* Min Q/U pixel SNR for fit */
  InfoReal.flt = 3.0; type = OBIT_float;
  in->minQUSNR = 3.0;
  ObitInfoListGetTest(in->info, "minQUSNR", &type, dim, &InfoReal);
  if (type==OBIT_float) in->minQUSNR = InfoReal.flt;
  else if (type==OBIT_double) in->minQUSNR = (ofloat)InfoReal.dbl;

  /* Want Error analysis? */
  InfoReal.itg = (olong)FALSE; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doError", &type, dim, &InfoReal);
  in->doError = InfoReal.itg;

  /* Determine number of lambda^2 planes and initialize in */
  in->nlamb2 = nimage;
  in->QRMS       = g_malloc0(in->nlamb2*sizeof(ofloat));
  in->URMS       = g_malloc0(in->nlamb2*sizeof(ofloat));
  in->inQFArrays = g_malloc0(in->nlamb2*sizeof(ObitFArray*));
  in->inUFArrays = g_malloc0(in->nlamb2*sizeof(ObitFArray*));
  in->lamb2      = g_malloc0(in->nlamb2*sizeof(odouble));

  /* How many output planes? */
  if (in->doError) nOut = 1+in->nterm*2;
  else nOut = in->nterm;
    
  /* Define term arrays */
  in->outFArrays = g_malloc0((nOut)*sizeof(ObitFArray*));
  for (i=0; i<nOut; i++) in->outFArrays[i] = NULL;

  /* Loop over images */
  for (iplane = 0; iplane<nimage; iplane++) {
    /* Open input image to get info */
    IOBy = OBIT_IO_byPlane;
    dim[0] = 1;
    ObitInfoListAlwaysPut (imQArr[iplane]->info, "IOBy", OBIT_long, dim, &IOBy);
    imQArr[iplane]->extBuffer = TRUE;   /* Using inFArrays as I/O buffer */
    retCode = ObitImageOpen (imQArr[iplane], OBIT_IO_ReadOnly, err);
    ObitInfoListAlwaysPut (imUArr[iplane]->info, "IOBy", OBIT_long, dim, &IOBy);
    imUArr[iplane]->extBuffer = TRUE;   /* Using inFArrays as I/O buffer */
    retCode = ObitImageOpen (imUArr[iplane], OBIT_IO_ReadOnly, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) 
      Obit_traceback_msg (err, routine, imQArr[iplane]->name);
    
    /* On first image initialize */
    if (iplane==0) {
      /* Image size */
      in->nx = imQArr[iplane]->myDesc->inaxes[0];
      in->ny = imQArr[iplane]->myDesc->inaxes[1];
      naxis[0] = (olong)in->nx;  naxis[1] = (olong)in->ny; 
      for (i=0; i<in->nlamb2; i++) in->inQFArrays[i]  = ObitFArrayCreate (NULL, 2, naxis);
      for (i=0; i<in->nlamb2; i++) in->inUFArrays[i]  = ObitFArrayCreate (NULL, 2, naxis);
      for (i=0; i<nOut; i++)       in->outFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);
      
      /* Output Image descriptor */
      outImage->myDesc = ObitImageDescCopy (imQArr[iplane]->myDesc, in->outDesc, err);
      if (err->error) Obit_traceback_msg (err, routine, imQArr[iplane]->name);

      /* Change third axis to type "SPECRM  " and leave the reference frequency
	 as the "CRVAL" */
      outImage->myDesc->inaxes[outImage->myDesc->jlocf] =  nOut;
      outImage->myDesc->crpix[outImage->myDesc->jlocf]  =  1.0;
      outImage->myDesc->cdelt[outImage->myDesc->jlocf]  =  1.0;
      strncpy (outImage->myDesc->ctype[outImage->myDesc->jlocf], SPECRM , IMLEN_KEYWORD);
      outImage->myDesc->bitpix = -32;  /* Float it */

      /* Creation date today */
      today = ObitToday();
      strncpy (outImage->myDesc->date, today, IMLEN_VALUE);
      if (today) g_free(today);

      /* Copy of output descriptor to in */
      in->outDesc = ObitImageDescCopy (outImage->myDesc, in->outDesc, err);

    } else { /* On subsequent images check for consistency */
      /* Check size of planes */
      Obit_return_if_fail(((imQArr[iplane]->myDesc->inaxes[0]==in->outDesc->inaxes[0]) && 
			   (imQArr[iplane]->myDesc->inaxes[1]==in->outDesc->inaxes[1])), err,
			  "%s: Image planes incompatible  %d!= %d or  %d!= %d", 
			  routine, imQArr[iplane]->myDesc->inaxes[0], in->outDesc->inaxes[0], 
			  imQArr[iplane]->myDesc->inaxes[1], in->outDesc->inaxes[1]) ;
      Obit_return_if_fail(((imUArr[iplane]->myDesc->inaxes[0]==in->outDesc->inaxes[0]) && 
			   (imUArr[iplane]->myDesc->inaxes[1]==in->outDesc->inaxes[1])), err,
			  "%s: Image planes incompatible  %d!= %d or  %d!= %d", 
			  routine, imUArr[iplane]->myDesc->inaxes[0], in->outDesc->inaxes[0], 
			  imUArr[iplane]->myDesc->inaxes[1], in->outDesc->inaxes[1]) ;

      /* Check alignment of pixels */
      ipixel[0] = 1.0; ipixel[1] = 1.0;
      bad = !ObitImageDescCvtPixel (imQArr[iplane]->myDesc, in->outDesc, ipixel, opixel, err);
      if (err->error) Obit_traceback_msg (err, routine, imQArr[iplane]->name);
      Obit_return_if_fail(!bad, err,
			  "%s: Image planes incompatible", routine);
      Obit_return_if_fail(((fabs(ipixel[0]-opixel[0])<0.01) &&
			   (fabs(ipixel[1]-opixel[1])<0.01)), err,
			  "%s: Image pixels not aligned %f!=%f or %f!=%f", 
			  routine, ipixel[0], opixel[0], ipixel[1], opixel[1]) ;
     } /* end consistency check */
    
    /* Read planes */
    retCode = ObitImageRead (imQArr[iplane], in->inQFArrays[iplane]->array, err);
    retCode = ObitImageRead (imUArr[iplane], in->inUFArrays[iplane]->array, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) Obit_traceback_msg (err, routine, imQArr[iplane]->name);
    
    /* Plane RMSes */
    in->QRMS[iplane] = ObitFArrayRMS(in->inQFArrays[iplane]);
    in->URMS[iplane] = ObitFArrayRMS(in->inUFArrays[iplane]);

    /* Lambda^2 */
    freq = imQArr[iplane]->myDesc->crval[imQArr[iplane]->myDesc->jlocf] + 
      imQArr[iplane]->myDesc->cdelt[imQArr[iplane]->myDesc->jlocf] * 
      (imQArr[iplane]->myDesc->plane - imQArr[iplane]->myDesc->crpix[imQArr[iplane]->myDesc->jlocf]);
    in->lamb2[iplane] = (VELIGHT/freq)*(VELIGHT/freq);

    /* Close inputs */
    retCode = ObitImageClose (imQArr[iplane], err);
    imQArr[iplane]->extBuffer = FALSE;   /* May need I/O buffer later */
    retCode = ObitImageClose (imUArr[iplane], err);
    imUArr[iplane]->extBuffer = FALSE;   /* May need I/O buffer later */
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) 
      Obit_traceback_msg (err, routine, imQArr[iplane]->name);
  } /* end loop over input images */

  /* Average Lambda^2 */
  avgLamb2 = 0.0;
  for (i=0; i<in->nlamb2; i++) avgLamb2 += in->lamb2[i];
  avgLamb2 /= in->nlamb2;

  /* Get Reference lambda^2 , default avg. input ref. freq. */
  InfoReal.dbl = avgLamb2; type = OBIT_double;
  in->refLamb2 = InfoReal.dbl;
  ObitInfoListGetTest(in->info, "refLamb2", &type, dim, &InfoReal);
  if (type==OBIT_float) in->refLamb2 = InfoReal.flt;
  else if (type==OBIT_double) in->refLamb2 = (ofloat)InfoReal.dbl;
  
  /* Update output header to reference Frequency */
  outImage->myDesc->crval[outImage->myDesc->jlocf] =VELIGHT/sqrt(in->refLamb2);
  in->outDesc->crval[in->outDesc->jlocf] = VELIGHT/sqrt(in->refLamb2);

  /* Do actual fitting */
  Fitter (in, err);
  if (err->error) Obit_traceback_msg (err, routine, imQArr[0]->name);

  /* Write output */
  WriteOutput(in, outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

} /* end ObitRMFitImArr */

/**
 * Fit single RM to Q, Uflux measurements
 * \param nlamb2    Number of entries in freq, flux, sigma
 * \param nterm     Number of coefficients of powers of log(nu) to fit
 * \param refLamb2  Reference lambda^2 (m^2)
 * \param lamb2     Array of lambda^2 (m^2)
 * \param qflux     Array of Q fluxes (Jy) same dim as lamb2
 * \param qsigma    Array of Q errors (Jy) same dim as lamb2
 * \param uflux     Array of U fluxes (Jy) same dim as lamb2
 * \param usigma    Array of U errors (Jy) same dim as lamb2
 * \param err      Obit error stack object.
 * \return  Array of fitter parameters, errors for each and Chi Squares of fit
 *          Initial terms are in Jy, other in log.
 */
ofloat* ObitRMFitSingle (olong nlamb2, olong nterm, odouble refLamb2, odouble *lamb2, 
			 ofloat *qflux, ofloat *qsigma, ofloat *uflux, ofloat *usigma, 
			 ObitErr *err)
{
  ofloat *out = NULL;
  olong i, j;
  NLRMFitArg *arg=NULL;
  gchar *routine = "ObitRMFitSingle";
  /* Warn if no GSL implementation */
#ifndef HAVE_GSL
  Obit_log_error(err, OBIT_InfoWarn, "NO GSL available - results will be zeroes");
#endif
  if (err->error) return out;


  /* Warn if ref. lambda^2 <= 0, set to 1.0 */
  if (refLamb2<=0.0) {
    refLamb2 = 1.0;
    Obit_log_error(err, OBIT_InfoWarn, 
		   "%s: Setting reference lambda^2 to 1", routine);
  }

  /* Create function argument */
  arg = g_malloc(sizeof(NLRMFitArg));
  arg->in             = NULL;     /* Not needed here */
  arg->nlamb2         = nlamb2;
  arg->nterm          = nterm;
  arg->doError        = TRUE;
  arg->minQUSNR       = 3.0;      /* min pixel SNR */
  arg->refLamb2       = refLamb2; /* Reference Frequency */
  arg->Qweight        = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Uweight        = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Qvar           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Uvar           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Qobs           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Uobs           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->lamb2          = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->x              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->y              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->w              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->coef           = g_malloc0(2*arg->nterm*sizeof(ofloat));
  for (i=0; i<nlamb2; i++) {
    arg->Qvar[i]    = qsigma[i]*qsigma[i];
    arg->Uvar[i]    = usigma[i]*usigma[i];
    arg->Qweight[i] = 1.0 / arg->Qvar[i];
    arg->Uweight[i] = 1.0 / arg->Uvar[i];
    arg->Qobs[i]    = qflux[i];
    arg->Uobs[i]    = uflux[i];
    arg->lamb2[i]   = lamb2[i];
  }

  /* Do fit */  
  NLRMFit(arg);
  
  /* get results parameters + errors */
  out = g_malloc0((2*arg->nterm+1)*sizeof(ofloat));
  for (j=0; j<nterm*2+1; j++) out[j] = arg->coef[j];
  
  /* Cleanup */
  if (arg->Qweight)   g_free(arg->Qweight);
  if (arg->Uweight)   g_free(arg->Uweight);
  if (arg->Qvar)      g_free(arg->Qvar);
  if (arg->Uvar)      g_free(arg->Uvar);
  if (arg->Qobs)      g_free(arg->Qobs);
  if (arg->Uobs)      g_free(arg->Uobs);
  if (arg->lamb2)     g_free(arg->lamb2);
  if (arg->x)         g_free(arg->x);
  if (arg->y)         g_free(arg->y);
  if (arg->w)         g_free(arg->w);
  if (arg->coef)      g_free(arg->coef);
  g_free(arg);

  return out;
} /* end ObitRMFitSingle */

/**
 * Make single spectrum fitting argument array
 * \param nlamb2   Number of entries in freq, flux, sigma
 * \param nterm    Number of coefficients of powers of log(nu) to fit
 * \param refLamb2 Reference lambda^2 (m^2)
 * \param lamb2    Array of lambda^2 (m^2)
 * \param out      Array for output results, should be g_freed when done
 * \param err      Obit error stack object.
 * \return  argument for single spectrum fitting, use ObitRMFitKillArg to dispose.
 */
gpointer ObitRMFitMakeArg (olong nlamb2, olong nterm, odouble refLamb2, 
			   odouble *lamb2, ofloat **out, ObitErr *err)
{
  olong i;
  NLRMFitArg *arg=NULL;
  gchar *routine = "ObitRMFitMakeArg";
  if (err->error) return out;

  /* Warn if too many terms asked for */
  if (nterm>2) {
      Obit_log_error(err, OBIT_InfoWarn, 
		    "%s: Asked for %d terms, will limit to 2", routine, nterm);
  }

  /* Warn if ref. lambda^2 <= 0, set to 1.0 */
  if (refLamb2<=0.0) {
    refLamb2 = 1.0;
    Obit_log_error(err, OBIT_InfoWarn, 
		   "%s: Setting reference Lambda^2 to 1", routine);
  }

  /* Create function argument */
  arg = g_malloc(sizeof(NLRMFitArg));
  arg->in             = NULL;     /* Not needed here */
  arg->nlamb2         = nlamb2;
  arg->nterm          = nterm;
  arg->doError        = TRUE;
  arg->minQUSNR       = 3.0;      /* min pixel SNR */
  arg->refLamb2       = refLamb2; /* Reference Frequency */
  arg->Qweight        = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Uweight        = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Qvar           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Uvar           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Qobs           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Uobs           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->lamb2          = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->x              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->y              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->w              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->coef           = g_malloc0(2*arg->nterm*sizeof(ofloat));
  for (i=0; i<nlamb2; i++) {
    arg->lamb2[i]        = lamb2[i];
  }

  /* output array */
  *out = (ofloat*)g_malloc0((2*arg->nterm+1)*sizeof(ofloat));

  return (gpointer)arg;
} /* end ObitRMFitMakeArg */

/**
 * Fit single spectrum to flux measurements using precomputed argument
 * \param aarg      pointer to argument for fitting
 * \param Qflux     Array of Q values to be fitted
 * \param Qsigma    Array of uncertainties of Qflux
 * \param Uflux     Array of U values to be fitted
 * \param Usigma    Array of uncertainties of Qflux
 * \param out       Result array at least 2*nterms+1 in size.
 *                  in order, fitted parameters, error estimates, chi sq of fit.
 */
void ObitRMFitSingleArg (gpointer aarg, 
			 ofloat *qflux, ofloat *qsigma, 
			 ofloat *uflux, ofloat *usigma,
			 ofloat *out)
{
  olong i, j;
  NLRMFitArg *arg=(NLRMFitArg*)aarg;
  ofloat fblank = ObitMagicF();
  gboolean allBad;

  /* Save flux array, sigma, Check if all data blanked  */
  allBad = TRUE;
  for (i=0; i<arg->nlamb2; i++) {
    if ((arg->Qobs[i]!=fblank) && 
	(arg->Uobs[i]!=fblank) && 
	((fabs(qflux[i])>arg->minQUSNR*qsigma[i]) ||
	 (fabs(uflux[i])>arg->minQUSNR*usigma[i]))) allBad = FALSE;
    arg->Qobs[i]    = qflux[i];
    arg->Uobs[i]    = uflux[i];
    arg->Qvar[i]    = qsigma[i];
    arg->Uvar[i]    = usigma[i];
    arg->Qweight[i] = 1.0 / arg->Qvar[i];
    arg->Uweight[i] = 1.0 / arg->Uvar[i];
  }
  
   /* Return fblanks/zeroes for no data */
  if (allBad) {
    out[0] = fblank;
    for (j=1; j<=arg->nterm*2; j++) out[j] = 0.0;
    return;
  }

  /* Fit */
  NLRMFit(arg);
  
  /* get results parameters + errors */
  for (j=0; j<arg->nterm*2; j++) out[j] = arg->coef[j];
  /* Chi squared */
  out[arg->nterm*2] = arg->ChiSq;
  
  return;
} /* end ObitRMFitSingleArg */

/**
 * Delete single fitting argument
 * \param aarg      pointer to argument to kill
 */
void ObitRMFitKillArg (gpointer aarg)
{
  NLRMFitArg *arg= (NLRMFitArg*)aarg;

  if (arg==NULL) return;
  if (arg->Qweight)   g_free(arg->Qweight);
  if (arg->Uweight)   g_free(arg->Uweight);
  if (arg->Qvar)      g_free(arg->Qvar);
  if (arg->Uvar)      g_free(arg->Uvar);
  if (arg->Qobs)      g_free(arg->Qobs);
  if (arg->Uobs)      g_free(arg->Uobs);
  if (arg->lamb2)     g_free(arg->lamb2);
  if (arg->x)         g_free(arg->x);
  if (arg->y)         g_free(arg->y);
  if (arg->w)         g_free(arg->w);
  if (arg->coef)      g_free(arg->coef);
  g_free(arg);

} /* end ObitRMFitKillArg */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitRMFitClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitRMFitClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitRMFitClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitRMFitClassInfoDefFn (gpointer inClass)
{
  ObitRMFitClassInfo *theClass = (ObitRMFitClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitRMFitClassInit;
  theClass->newObit       = (newObitFP)newObitRMFit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitRMFitClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitRMFitGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitRMFitCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitRMFitClear;
  theClass->ObitInit      = (ObitInitFP)ObitRMFitInit;
  theClass->ObitRMFitCreate = (ObitRMFitCreateFP)ObitRMFitCreate;
  theClass->ObitRMFitCube   = (ObitRMFitCubeFP)ObitRMFitCube;
  theClass->ObitRMFitImArr  = (ObitRMFitImArrFP)ObitRMFitImArr;
  theClass->ObitRMFitSingle = (ObitRMFitSingleFP)ObitRMFitSingle;
} /* end ObitRMFitClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitRMFitInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitRMFit *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread     = newObitThread();
  in->info       = newObitInfoList(); 
  in->nterm      = 2;
  in->nlamb2     = 0;
  in->minQUSNR   = 3.0;
  in->QRMS       = NULL;
  in->URMS       = NULL;
  in->outDesc    = NULL;
  in->inQFArrays = NULL;
  in->inUFArrays = NULL;
  in->outFArrays = NULL;
  in->lamb2      = NULL;

} /* end ObitRMFitInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitRMFit* cast to an Obit*.
 */
void ObitRMFitClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitRMFit *in = inn;
  olong i, nOut;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->QRMS) g_free(in->QRMS);
  if (in->URMS) g_free(in->URMS);
  in->thread = ObitThreadUnref(in->thread);
  in->info   = ObitInfoListUnref(in->info);
  if (in->outDesc) in->outDesc = ObitImageDescUnref(in->outDesc);
  if (in->inQFArrays) {
    for (i=0; i<in->nlamb2; i++) in->inQFArrays[i] = ObitFArrayUnref(in->inQFArrays[i]);
    g_free(in->inQFArrays);
  }
  if (in->inUFArrays) {
    for (i=0; i<in->nlamb2; i++) in->inUFArrays[i] = ObitFArrayUnref(in->inUFArrays[i]);
    g_free(in->inUFArrays);
  }

  /* How many output planes */
  if (in->doError) nOut = 1+in->nterm*2;
  else nOut = in->nterm;

  if (in->outFArrays) {
    for (i=0; i<nOut; i++) in->outFArrays[i] = ObitFArrayUnref(in->outFArrays[i]);
    g_free(in->outFArrays);
  }
  if (in->lamb2)     g_free(in->lamb2);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitRMFitClear */

/**
 * Does RMl fitting, input and output on input object
 * Work divided up amoung 1 or more threads
 * In each pixel an RM is fitted if the number of valid data points 
 * exceeds in->nterm, otherwise the pixel is blanked.
 * The results are stored in outFArrays:
 * \li first entry is the RM at the reference lambda^2
 * \li second is the EVLA (rad) at the reference lambda^2
 * \li entries nterm+1->nterm*2 RMS uncertainties on coefficients
 * \li entry 1+nterm*2 = the Chi squared of the fit
 *
 * \param  in  RMFit to fit
 * \param  err Obit error stack object.
 */
static void Fitter (ObitRMFit* in, ObitErr *err)
{
  olong i, loy, hiy, nyPerThread, nThreads;
  gboolean OK;
  NLRMFitArg **threadArgs;
  NLRMFitArg *args=NULL;
  ObitThreadFunc func=(ObitThreadFunc)ThreadNLRMFit;
  gchar *routine = "Fitter";

  /* error checks */
  if (err->error) return;

  /* Warn if too many terms asked for */
  if (in->nterm>2) {
      Obit_log_error(err, OBIT_InfoWarn, 
		    "%s: Asked for %d terms, will limit to 2", routine, in->nterm);
  }

  /* How many threads to use? */
  nThreads = MAX (1, ObitThreadNumProc(in->thread));

  /* Initialize threadArg array  */
  threadArgs = g_malloc0(nThreads*sizeof(NLRMFitArg*));
  for (i=0; i<nThreads; i++) 
    threadArgs[i] = g_malloc0(sizeof(NLRMFitArg)); 

  /* Set up thread arguments */
  for (i=0; i<nThreads; i++) {
    args = (NLRMFitArg*)threadArgs[i];
    args->in          = in;
    args->err         = err;
    args->nlamb2      = in->nlamb2;
    args->nterm       = in->nterm;
    args->doError     = in->doError;
    args->minQUSNR    = in->minQUSNR;     /*min pixel SNR */
    args->Qweight     = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->Uweight     = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->Qvar        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->Uvar        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->Qobs        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->Uobs        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->lamb2       = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->x           = g_malloc0(args->nlamb2*sizeof(double));
    args->y           = g_malloc0(args->nlamb2*sizeof(double));
    args->w           = g_malloc0(args->nlamb2*sizeof(double));
    if (args->doError)
      args->coef      = g_malloc0(2*args->nterm*sizeof(ofloat));
    else
      args->coef      = g_malloc0(args->nterm*sizeof(ofloat));
  }
  /* end initialize */
  
  /* Divide up work */
  nyPerThread = in->ny/nThreads;
  loy = 1;
  hiy = nyPerThread;
  hiy = MIN (hiy, in->ny);
  
  /* Set up thread arguments */
  for (i=0; i<nThreads; i++) {
    if (i==(nThreads-1)) hiy = in->ny;  /* Make sure do all */
    args = (NLRMFitArg*)threadArgs[i];
    args->first  = loy;
    args->last   = hiy;
    if (nThreads>1) args->ithread = i;
    else args->ithread = -1;
    /* Update which y */
    loy += nyPerThread;
    hiy += nyPerThread;
    hiy = MIN (hiy, in->ny);
  }
  
  /* Do operation possibly with threads */
  OK = ObitThreadIterator (in->thread, nThreads, func, (gpointer)threadArgs);
  
  /* Check for problems */
  if (!OK) {
    Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
  }

  /* Shut down any threading */
  ObitThreadPoolFree (in->thread);
  if (threadArgs) {
    for (i=0; i<nThreads; i++) {
      args = (NLRMFitArg*)threadArgs[i];
      if (args->Qweight)   g_free(args->Qweight);
      if (args->Uweight)   g_free(args->Uweight);
      if (args->Qvar)      g_free(args->Qvar);
      if (args->Uvar)      g_free(args->Uvar);
      if (args->Qobs)      g_free(args->Qobs);
      if (args->Uobs)      g_free(args->Uobs);
      if (args->lamb2)     g_free(args->lamb2);
      if (args->x)         g_free(args->x);
      if (args->y)         g_free(args->y);
      if (args->w)         g_free(args->w);
      if (args->coef)      g_free(args->coef);
     g_free(threadArgs[i]);
    }
    g_free(threadArgs);
  }

} /* end Fitter */

/**
 * Write contents on in to outImage
 * \param in       RM fitting object
 * \param outImage Image cube with fitted spectra.
 *                 Should be defined but not created.
 *                 Planes 1->nterm are coefficients per pixel
 *                 Planes nterm+1->2*nterm are uncertainties in coefficients
 * \param err      Obit error stack object.
 */
static void WriteOutput (ObitRMFit* in, ObitImage *outImage, 
			 ObitErr *err)
{
  olong iplane, nOut;
  ObitIOSize IOBy;
  olong  i, blc[IM_MAXDIM], trc[IM_MAXDIM];
  ObitIOCode retCode;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "WriteOutput";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert(ObitImageIsA(outImage));

  /* Open output- all, by plane */
  dim[0] = IM_MAXDIM;
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  ObitInfoListPut (outImage->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (outImage->info, "TRC", OBIT_long, dim, trc, err); 
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (outImage->info, "IOBy", OBIT_long, dim, &IOBy);
  outImage->extBuffer = TRUE;   /* Using outFArrays as I/O buffer */
  retCode = ObitImageOpen (outImage, OBIT_IO_ReadWrite, err);
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, outImage->name);

  /* save descriptive material */
  ObitImageDescCopyDesc (in->outDesc, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* How many output planes? */
  if (in->doError) nOut = 1+in->nterm*2;
  else nOut = in->nterm;

  /* Loop writing planes */
  for (iplane=0; iplane<nOut; iplane++) {
    retCode = ObitImageWrite (outImage, in->outFArrays[iplane]->array, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) Obit_traceback_msg (err, routine, outImage->name);

  } /* end loop writing planes */

  /* Save nterms on output descriptor */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (outImage->myDesc->info, "NTERM", OBIT_long, dim, &in->nterm);
  ObitInfoListAlwaysPut (((ObitImageDesc*)outImage->myIO->myDesc)->info, "NTERM", 
    OBIT_long, dim, &in->nterm);

  /* Close output */
  retCode = ObitImageClose (outImage, err);
  outImage->extBuffer = FALSE;   /* May need I/O buffer later */
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, outImage->name);

} /* end WriteOutput */

/**
 * Thread function to fit a portion of the image set
 * \param arg      NLRMFitArg structure
 */
static gpointer ThreadNLRMFit (gpointer arg)
{
  NLRMFitArg *larg      = (NLRMFitArg*)arg;
  ObitRMFit* in       = (ObitRMFit*)larg->in;
  olong lo            = larg->first-1;  /* First in y range */
  olong hi            = larg->last;     /* Highest in y range */
  gboolean doError    = larg->doError;  /* Error analysis? */
  ObitErr *err        = larg->err;

  olong ix, iy, indx, i, nOut;
  ofloat fblank = ObitMagicF();
  /*gchar *routine = "ThreadNLRMFit";*/

  /* error checks */
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Set up frequency info */
  larg->refLamb2 = in->refLamb2;
  for (i=0; i<in->nlamb2; i++) {
    larg->lamb2[i] = in->lamb2[i];
  }

  /* How many output planes */
  if (in->doError) nOut = 1+in->nterm*2;
  else nOut = in->nterm;

  /* Loop over pixels in Y */
  indx = lo*in->nx  -1;  /* Offset in pixel arrays */
  for (iy=lo; iy<hi; iy++) {
    /* Loop over pixels in X */
    for (ix=0; ix<in->nx; ix++) {
      indx ++;

      /* Collect values;  */
      for (i=0; i<in->nlamb2; i++) {
	larg->Qobs[i] = in->inQFArrays[i]->array[indx];
	larg->Uobs[i] = in->inUFArrays[i]->array[indx];
	/* Data valid? */
	if ((larg->Qobs[i]!=fblank) && 
	    (larg->Uobs[i]!=fblank) && 
	    ((fabs(larg->Qobs[i])>in->minQUSNR*in->QRMS[i]) ||
	     (fabs(larg->Uobs[i])>in->minQUSNR*in->URMS[i]))) {
	  /* Statistical weight */
	  larg->Qvar[i] = (in->QRMS[i]*in->QRMS[i]);
	  larg->Uvar[i] = (in->URMS[i]*in->URMS[i]);
	  larg->Qweight[i] = 1.0 / larg->Qvar[i];
	  larg->Uweight[i] = 1.0 / larg->Uvar[i];
	  /* End if datum valid */
	} else { /* invalid pixel */
	  larg->Qweight[i] = 0.0;
	  larg->Qvar[i] = 0.0;
	  larg->Uweight[i] = 0.0;
	  larg->Uvar[i] = 0.0;
	}
      } /* end loop over frequencies */
      
      /* Fit */
      NLRMFit(larg);
      
      /* Save to output */
      if (doError) {
	for (i=0; i<nOut; i++) 
	  in->outFArrays[i]->array[indx]  = larg->coef[i];
      } else { /* only values */
 	for (i=0; i<in->nterm; i++) 
	  in->outFArrays[i]->array[indx] = larg->coef[i];
      }

    } /* end x loop */
  } /* end y loop */

  /* Indicate completion */
  if (larg->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&larg->ithread);

  return NULL;
} /* end ThreadNLRMFit */

/**
 * Fit RM and EVPA at the reference lambda^2
 * Only fits for up to 5 terms
 * \param arg      NLRMFitArg structure
 *                 fitted parameters returned in arg->in->coef
 */
static void NLRMFit (NLRMFitArg *arg)
{
  olong i, nterm=arg->nterm, nvalid;
  ofloat fblank = ObitMagicF();
  odouble aarg, varaarg;
  double EVPA0, RM, varEVPA0, coVar, varRM, chi2;
  int status, numb;
 
  /* Initialize output */
  if (arg->doError) 
    for (i=0; i<2*nterm+1; i++) arg->coef[i] = 0.0;
  else
    for (i=0; i<nterm; i++) arg->coef[i]  = 0.0;
  /* Blank */
  arg->coef[0] =  arg->coef[1] = fblank;
  
  /* get EVLA values count valid data */
  nvalid = 0;
  numb = 0;
  for (i=0; i<arg->nlamb2; i++) {
    if ((arg->Qobs[i]!=fblank) && (arg->Qweight[i]>0.0) &&
	(arg->Uobs[i]!=fblank) && (arg->Uweight[i]>0.0)) {
      arg->x[numb] = arg->lamb2[i] - arg->refLamb2;
      arg->y[numb] = 0.5*atan2(arg->Uobs[i],arg->Qobs[i]);
      /* Weight */
      aarg = fabs(arg->Uobs[i]/arg->Qobs[i]);
      varaarg = aarg*aarg*(arg->Qvar[i]/(arg->Qobs[i]*arg->Qobs[i]) + 
			   arg->Uvar[i]/(arg->Uobs[i]*arg->Uobs[i]));
      arg->w[numb++] = (1.0+aarg*aarg)/varaarg;
    }
  }
  nvalid = numb;
  if (nvalid<=MAX(2,arg->nterm)) return;  /* enough good data? */

  /* try to unwrap EVPA */
  for (i=1; i<numb; i++) {
    if ((arg->y[i]-arg->y[i-1]) >  0.75*G_PI) arg->y[i] -= G_PI;
    if ((arg->y[i]-arg->y[i-1]) < -0.75*G_PI) arg->y[i] += G_PI;
  }

  /* Do fit */
#ifdef HAVE_GSL
  status = gsl_fit_wlinear (arg->x, 1, arg->w, 1, arg->y, 1, numb, 
			    &EVPA0, &RM, &varEVPA0, &coVar, &varRM, &chi2);
#else
  /* dummy values */
  EVPA0 = RM = varEVPA0 = coVar = varRM = chi2 = 0.0;
#endif /* HAVE_GSL */ 
  /* save values */
  if (nterm==1) arg->coef[0] = (ofloat)RM;
  else if (nterm==2) {
    arg->coef[0] = (ofloat)RM;
    arg->coef[1] = (ofloat)EVPA0;
    if (arg->doError) {
      arg->coef[2] = (ofloat)sqrt(varRM);
      arg->coef[3] = (ofloat)sqrt(varEVPA0);
      arg->coef[4] = (ofloat)chi2;
    }
  }
} /* end NLRMFit */

