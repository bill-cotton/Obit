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
#include "ObitSinCos.h"
#ifdef HAVE_GSL
#include <gsl/gsl_blas.h>
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


#ifdef HAVE_GSL
/** Private: Solver function evaluation */
static int RMFitFunc (const gsl_vector *x, void *params, 
		      gsl_vector *f);

/** Private: Solver Jacobian evaluation */
static int RMFitJac (const gsl_vector *x, void *params, 
		     gsl_matrix *J);

/** Private: Solver function + Jacobian evaluation */
static int RMFitFuncJac (const gsl_vector *x, void *params, 
			 gsl_vector *f, gsl_matrix *J);

#endif /* HAVE_GSL */ 

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
  /** number of valid data points in x, q, u, w  */
  olong        nvalid;
  /** maximum iteration  */
  olong        maxIter;
  /** acceptable iteration delta, rel and abs. */
  odouble minDelta;
  /** Array of Lambda^2 (nlamb2) */
  ofloat *lamb2;
  /** Array of Q, U weights (1/RMS) per inFArrays (nlamb2) */
  ofloat *Qweight, *Uweight;
  /** Array of Q, U variance per inFArrays (nlamb2) */
  ofloat *Qvar, *Uvar;
  /** Array of Q, U pixel values being fitted (nlamb2) */
  ofloat *Qobs, *Uobs;
  /** Array of polarized intensity pixel values being fitted (nlamb2) */
  ofloat *Pobs;
  /** min Q/U SNR */
  ofloat minQUSNR;
  /** min fraction of valid samples */
  ofloat minFrac;
  /** Reference lambda^2 */
  ofloat refLamb2;
  /** Vector of guess/fitted coefficients, optional errors */
  ofloat *coef;
  /** Chi squared of fit */
  ofloat ChiSq;
  /** Do error analysis? */
  gboolean doError;
  /** work arrays */
  double *x, *q, *u, *w;
  ofloat *wrk1, *wrk2, *wrk3;
#ifdef HAVE_GSL
  /** Fitting solver  */
  gsl_multifit_fdfsolver *solver;
  /** Fitting solver function structure */
  gsl_multifit_function_fdf *funcStruc;
  /** Fitting work vector */
  gsl_vector *work;
  /** Covariance matrix */
  gsl_matrix *covar;
#endif /* HAVE_GSL */ 
} NLRMFitArg;

/** Private: Actual fitting */
static void NLRMFit (NLRMFitArg *arg);

/** Private: coarse search */
static olong RMcoarse (NLRMFitArg *arg);

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
  ofloat s, c;

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
  ObitSinCosCalc(0.0, &s, &c);    /* Sine/cosine functions */

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
 * \li "minFrac"  OBIT_float  scalar min. fraction of planes included [def 0.5]
 * \li "doError"  OBIT_boolean scalar If true do error analysis [def False]
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
  Obit_log_error(err, OBIT_InfoWarn, "NO GSL available - results will be approximate");
#endif

  /* Control parameters */
  /* Min Q/U pixel SNR for fit */
  InfoReal.flt = 3.0; type = OBIT_float;
  in->minQUSNR = 3.0;
  ObitInfoListGetTest(in->info, "minQUSNR", &type, dim, &InfoReal);
  if (type==OBIT_float) in->minQUSNR = InfoReal.flt;
  else if (type==OBIT_double) in->minQUSNR = (ofloat)InfoReal.dbl;

  /* Min fraction of planes in the fit */
  InfoReal.flt = 0.5; type = OBIT_float;
  in->minFrac  = 0.5;
  ObitInfoListGetTest(in->info, "minFrac", &type, dim, &InfoReal);
  if (type==OBIT_float) in->minFrac = InfoReal.flt;
  else if (type==OBIT_double) in->minFrac = (ofloat)InfoReal.dbl;

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
	freq = inQImage->myDesc->crval[inQImage->myDesc->jlocf] + 
	  inQImage->myDesc->cdelt[inQImage->myDesc->jlocf] * 
	  (inQImage->myDesc->plane - inQImage->myDesc->crpix[inQImage->myDesc->jlocf]);
	ObitInfoListGetTest (inQImage->myDesc->info, keyword, &type, dim, &freq);
       } else {   /* Normal spectral cube */
	freq = inQImage->myDesc->crval[inQImage->myDesc->jlocf] + 
	  inQImage->myDesc->cdelt[inQImage->myDesc->jlocf] * 
	  (inQImage->myDesc->plane - inQImage->myDesc->crpix[inQImage->myDesc->jlocf]);
      }
    in->lamb2[iplane] = (VELIGHT/freq)*(VELIGHT/freq);

    plane[0] = iplane+noffset+1;  /* Select correct plane */
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
  Obit_log_error(err, OBIT_InfoWarn, "NO GSL available - results will be approximate");
#endif

  /* Control parameters */
  /* Min Q/U pixel SNR for fit */
  InfoReal.flt = 3.0; type = OBIT_float;
  in->minQUSNR = 3.0;
  ObitInfoListGetTest(in->info, "minQUSNR", &type, dim, &InfoReal);
  if (type==OBIT_float) in->minQUSNR = InfoReal.flt;
  else if (type==OBIT_double) in->minQUSNR = (ofloat)InfoReal.dbl;

  /* Min fraction of planes in the fit */
  InfoReal.flt = 0.5; type = OBIT_float;
  in->minFrac  = 0.5;
  ObitInfoListGetTest(in->info, "minFrac", &type, dim, &InfoReal);
  if (type==OBIT_float) in->minFrac = InfoReal.flt;
  else if (type==OBIT_double) in->minFrac = (ofloat)InfoReal.dbl;

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
      for (i=0; i<in->nlamb2; i++) in->inQFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);
      for (i=0; i<in->nlamb2; i++) in->inUFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);
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
  ofloat fblank = ObitMagicF();
  olong i, j;
  NLRMFitArg *arg=NULL;
  gchar *routine = "ObitRMFitSingle";
  /* GSL implementation */
#ifdef HAVE_GSL
  const gsl_multifit_fdfsolver_type *T=NULL;
#endif /* HAVE_GSL */ 
  /* Warn if no GSL implementation */
#ifndef HAVE_GSL
  Obit_log_error(err, OBIT_InfoWarn, "NO GSL available - results will be approximate");
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
  arg->maxIter        = 100;     /* max. number of iterations */
  arg->minDelta       = 1.0e-5;  /* Min step size */
  arg->nterm          = nterm;
  arg->doError        = TRUE;
  arg->minQUSNR       = 3.0;      /* min pixel SNR */
  arg->minFrac        = 0.5;      /* min fraction of samples */
  arg->refLamb2       = refLamb2; /* Reference Frequency */
  arg->Qweight        = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Uweight        = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Qvar           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Uvar           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Qobs           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Uobs           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Pobs           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->lamb2          = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->x              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->w              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->q              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->u              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->wrk1           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->wrk2           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->wrk3           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->coef           = g_malloc0(3*arg->nterm*sizeof(ofloat));
  for (i=0; i<nlamb2; i++) {
    arg->lamb2[i]   = lamb2[i];
    arg->Qvar[i]    = qsigma[i]*qsigma[i];
    arg->Uvar[i]    = usigma[i]*usigma[i];
    arg->Qobs[i]    = qflux[i];
    arg->Uobs[i]    = uflux[i];
    if ((arg->Qobs[i]!=fblank) && (arg->Uobs[i]!=fblank)) {
      arg->Pobs[i]    = sqrt (qflux[i]*qflux[i] + uflux[i]*uflux[i]);
      arg->Qweight[i] = 1.0 / qsigma[i];
      arg->Uweight[i] = 1.0 / usigma[i];
   } else{ 
      arg->Pobs[i]    = 0.0;
      arg->Qweight[i] = 0.0;
      arg->Uweight[i] = 0.0;
    }
  }
  /* GSL implementation */
#ifdef HAVE_GSL
  arg->solver = NULL;
  arg->covar  = NULL;
  arg->work = NULL;
  /* Setup solver */
  T = gsl_multifit_fdfsolver_lmder;
  arg->solver = gsl_multifit_fdfsolver_alloc(T, 2*arg->nlamb2, 2);

  /* Fitting function info */
  arg->funcStruc = g_malloc0(sizeof(gsl_multifit_function_fdf));
  arg->funcStruc->f      = &RMFitFunc;
  arg->funcStruc->df     = &RMFitJac;
  arg->funcStruc->fdf    = &RMFitFuncJac;
  arg->funcStruc->n      = 2*arg->nlamb2;
  arg->funcStruc->p      = 2;
  arg->funcStruc->params = arg;

  /* Set up work arrays */
  arg->covar = gsl_matrix_alloc(2, 2);
  arg->work  = gsl_vector_alloc(2);
#endif /* HAVE_GSL */ 

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
  if (arg->Pobs)      g_free(arg->Pobs);
  if (arg->lamb2)     g_free(arg->lamb2);
  if (arg->x)         g_free(arg->x);
  if (arg->w)         g_free(arg->w);
  if (arg->q)         g_free(arg->q);
  if (arg->u)         g_free(arg->u);
  if (arg->wrk1)      g_free(arg->wrk1);
  if (arg->wrk2)      g_free(arg->wrk2);
  if (arg->wrk3)      g_free(arg->wrk3);
  if (arg->coef)      g_free(arg->coef);
#ifdef HAVE_GSL
  if (arg->solver)   gsl_multifit_fdfsolver_free (arg->solver);
  if (arg->work)     gsl_vector_free(arg->work);
  if (arg->covar)    gsl_matrix_free(arg->covar);
  if (arg->funcStruc) g_free(arg->funcStruc);
#endif /* HAVE_GSL */
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
#ifdef HAVE_GSL
  const gsl_multifit_fdfsolver_type *T=NULL;
#endif /* HAVE_GSL */ 

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
  arg->maxIter        = 100;     /* max. number of iterations */
  arg->minDelta       = 1.0e-5;  /* Min step size */
  arg->nterm          = nterm;
  arg->doError        = TRUE;
  arg->minQUSNR       = 3.0;      /* min pixel SNR */
  arg->minFrac        = 0.5;      /* min fraction of samples */
  arg->refLamb2       = refLamb2; /* Reference Frequency */
  arg->Qweight        = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Uweight        = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Qvar           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Uvar           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Qobs           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Uobs           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->Pobs           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->lamb2          = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->x              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->w              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->q              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->u              = g_malloc0(arg->nlamb2*sizeof(double));
  arg->wrk1           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->wrk2           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->wrk3           = g_malloc0(arg->nlamb2*sizeof(ofloat));
  arg->coef           = g_malloc0(3*arg->nterm*sizeof(ofloat));
  for (i=0; i<nlamb2; i++) {
    arg->lamb2[i]        = lamb2[i];
  }
  /* GSL implementation */
#ifdef HAVE_GSL
  arg->solver = NULL; arg->covar = NULL; arg->work = NULL;
  /* Setup solver */
  T = gsl_multifit_fdfsolver_lmder;
  arg->solver = gsl_multifit_fdfsolver_alloc(T, 2*arg->nlamb2, 2);

  /* Fitting function info */
  arg->funcStruc = g_malloc0(sizeof(gsl_multifit_function_fdf));
  arg->funcStruc->f      = &RMFitFunc;
  arg->funcStruc->df     = &RMFitJac;
  arg->funcStruc->fdf    = &RMFitFuncJac;
  arg->funcStruc->n      = 2*arg->nlamb2;
  arg->funcStruc->p      = 2;
  arg->funcStruc->params = arg;

  /* Set up work arrays */
  arg->covar = gsl_matrix_alloc(2, 2);
  arg->work  = gsl_vector_alloc(2);
#endif /* HAVE_GSL */ 

  /* output array */
  *out = (ofloat*)g_malloc0((2*arg->nterm+1)*sizeof(ofloat));

  return (gpointer)arg;
} /* end ObitRMFitMakeArg */

/**
 * Fit single RM to measurements using precomputed argument
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
    if ((arg->Qobs[i]!=fblank) && (arg->Uobs[i]!=fblank)) {
      arg->Pobs[i]    = sqrt (qflux[i]*qflux[i] + uflux[i]*uflux[i]);
      /* arg->Qweight[i] = arg->Pobs[i] / arg->Qvar[i];
	 arg->Uweight[i] = arg->Pobs[i] / arg->Uvar[i]; */
      arg->Qweight[i] = 1.0 / qsigma[i];
      arg->Uweight[i] = 1.0 / usigma[i];
   } else {
      arg->Pobs[i]    = 0.0;
      arg->Qweight[i] = 0.0;
      arg->Uweight[i] = 0.0;
    }
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
  if (arg->Pobs)      g_free(arg->Pobs);
  if (arg->lamb2)     g_free(arg->lamb2);
  if (arg->x)         g_free(arg->x);
  if (arg->w)         g_free(arg->w);
  if (arg->q)         g_free(arg->q);
  if (arg->u)         g_free(arg->u);
  if (arg->wrk1)      g_free(arg->wrk1);
  if (arg->wrk2)      g_free(arg->wrk2);
  if (arg->wrk3)      g_free(arg->wrk3);
  if (arg->coef)      g_free(arg->coef);
#ifdef HAVE_GSL
  if (arg->solver)   gsl_multifit_fdfsolver_free (arg->solver);
  if (arg->work)     gsl_vector_free(arg->work);
  if (arg->covar)    gsl_matrix_free(arg->covar);
  if (arg->funcStruc) g_free(arg->funcStruc);
#endif /* HAVE_GSL */
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
  in->minFrac    = 0.5;
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
 * Does RM fitting, input and output on input object
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
  /* GSL implementation */
#ifdef HAVE_GSL
  const gsl_multifit_fdfsolver_type *T=NULL;
#endif /* HAVE_GSL */ 

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
    args->maxIter     = 100;     /* max. number of iterations */
    args->minDelta    = 1.0e-5;  /* Min step size */
    args->nterm       = in->nterm;
    args->doError     = in->doError;
    args->minQUSNR    = in->minQUSNR;     /* min pixel SNR */
    args->minFrac     = in->minFrac;      /* min fraction of samples  */
    args->Qweight     = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->Uweight     = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->Qvar        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->Uvar        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->Qobs        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->Uobs        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->Pobs        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->lamb2       = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->x           = g_malloc0(args->nlamb2*sizeof(double));
    args->w           = g_malloc0(args->nlamb2*sizeof(double));
    args->q           = g_malloc0(args->nlamb2*sizeof(double));
    args->u           = g_malloc0(args->nlamb2*sizeof(double));
    args->wrk1        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->wrk2        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->wrk3        = g_malloc0(args->nlamb2*sizeof(ofloat));
    if (args->doError)
      args->coef      = g_malloc0(3*args->nterm*sizeof(ofloat));
    else
      args->coef      = g_malloc0(args->nterm*sizeof(ofloat));
   /* GSL implementation */
#ifdef HAVE_GSL
    args->solver = NULL; args->covar = NULL; args->work = NULL;
  /* Setup solver */
    T = gsl_multifit_fdfsolver_lmder;
    args->solver = gsl_multifit_fdfsolver_alloc(T, 2*args->nlamb2, 2);
    
    /* Fitting function info */
    args->funcStruc = g_malloc0(sizeof(gsl_multifit_function_fdf));
    args->funcStruc->f      = &RMFitFunc;
    args->funcStruc->df     = &RMFitJac;
    args->funcStruc->fdf    = &RMFitFuncJac;
    args->funcStruc->n      = 2*args->nlamb2;
    args->funcStruc->p      = 2;
    args->funcStruc->params = args;
    
    /* Set up work arrays */
    args->covar = gsl_matrix_alloc(2, 2);
    args->work  = gsl_vector_alloc(2);
#endif /* HAVE_GSL */ 
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
      if (args->Pobs)      g_free(args->Pobs);
      if (args->lamb2)     g_free(args->lamb2);
      if (args->x)         g_free(args->x);
      if (args->w)         g_free(args->w);
      if (args->q)         g_free(args->q);
      if (args->u)         g_free(args->u);
      if (args->wrk1)      g_free(args->wrk1);
      if (args->wrk2)      g_free(args->wrk2);
      if (args->wrk3)      g_free(args->wrk3);
      if (args->coef)      g_free(args->coef);
#ifdef HAVE_GSL
      if (args->solver)    gsl_multifit_fdfsolver_free (args->solver);
      if (args->work)      gsl_vector_free(args->work);
      if (args->covar)     gsl_matrix_free(args->covar);
      if (args->funcStruc) g_free(args->funcStruc);
#endif /* HAVE_GSL */
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
	  larg->Pobs[i]    = sqrt (larg->Qobs[i]*larg->Qobs[i] + larg->Uobs[i]*larg->Uobs[i]);
	  /* larg->Qweight[i] = larg->Pobs[i] / larg->Qvar[i];
	     larg->Uweight[i] = larg->Pobs[i] / larg->Uvar[i];*/
	  larg->Qweight[i] = 1.0 / in->QRMS[i];
	  larg->Uweight[i] = 1.0 / in->URMS[i];
	  /* End if datum valid */
	} else { /* invalid pixel */
	  larg->Qweight[i] = 0.0;
	  larg->Qvar[i]    = 0.0;
	  larg->Uweight[i] = 0.0;
	  larg->Uvar[i]    = 0.0;
	  larg->Pobs[i]    = 0.0;
	}
	/* DEBUG
	if ((ix==669) && (iy==449)) { 
	  fprintf (stderr,"%3d q=%g u=%g qs=%g us=%g wt=%f\n",
		   i, larg->Qobs[i], larg->Uobs[i], in->QRMS[i], in->URMS[i], larg->Qweight[i]);
	} */
      } /* end loop over frequencies */
      
      /* Fit */
      NLRMFit(larg);
      /* DEBUG
      if ((ix==669) && (iy==449)) { 
	 fprintf (stderr,"ix=%4d iy=%4d\n",ix+1,iy+1);
	 fprintf (stderr,"RM=%f EVPA=%f chi2=%f\n",
		  larg->coef[0], larg->coef[1],larg->coef[4]);
      } */
      
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
 *                 RM, EVPA, sig RM, sig EVPA, chi2
 */
static void NLRMFit (NLRMFitArg *arg)
{
  olong iter=0, i, nterm=2, nvalid;
  ofloat sumwt, fblank = ObitMagicF();
  double chi2;
  int status, numb;
 
  /* Initialize output */
  if (arg->doError) 
    for (i=0; i<2*nterm+1; i++) arg->coef[i] = 0.0;
  else
    for (i=0; i<nterm; i++) arg->coef[i]  = 0.0;
  /* Blank */
  arg->coef[0] =  arg->coef[1] = fblank;
  
  /* try to unwrap EVPA, get data to be fitted, returns number of valid data 
     and crude fits */
  nvalid = RMcoarse (arg);
  numb   = nvalid;

  if (nvalid<=MAX(2,nterm)) return;  /* enough good data for fit? */
  /* High enough fraction of valid pixels? */
  if ((((ofloat)nvalid)/((ofloat)arg->nlamb2)) < arg->minFrac) return;

  /* Do fit */
#ifdef HAVE_GSL
    /* order EVPA, RM */
    gsl_vector_set(arg->work, 0, (double)arg->coef[1]);
    gsl_vector_set(arg->work, 1, (double)arg->coef[0]);
    arg->funcStruc->n      = arg->nlamb2*2;
    arg->funcStruc->p      = nterm;
    arg->funcStruc->params = arg;
    gsl_multifit_fdfsolver_set (arg->solver, arg->funcStruc, arg->work);
    iter = 0;
    /* iteration loop */
    do {
      iter++;
      status = gsl_multifit_fdfsolver_iterate(arg->solver);
      /*if (status) break;???*/

      status = gsl_multifit_test_delta (arg->solver->dx, arg->solver->x, 
					(double)arg->minDelta, 
					(double)arg->minDelta);
      /* DEBUG
      if (nvalid>nterm) {
	sumwt = (ofloat)gsl_blas_dnrm2(arg->solver->f);
	chi2 = (sumwt*sumwt)/(nvalid-nterm);
      } else chi2 = -1.0;
      for (i=0; i<nterm; i++) arg->coef[i] = (ofloat)gsl_vector_get(arg->solver->x, i);
      fprintf (stderr,"   iter=%d RM %f EVPA %f chi2 %g status %d\n",
	       iter, arg->coef[1], arg->coef[0], chi2, status); */
      /* end DEBUG */
    } while ((status==GSL_CONTINUE) && (iter<arg->maxIter));

    /* If it didn't work - bail */
    if ((status!=GSL_SUCCESS) && (status!=GSL_CONTINUE)) {
      fprintf (stderr, "Failed, status = %s\n", gsl_strerror(status));
      return;
    }
    
    /* normalized Chi squares */
    if (nvalid>nterm) {
      sumwt = (ofloat)gsl_blas_dnrm2(arg->solver->f);
      chi2 = (sumwt*sumwt)/(nvalid-nterm);
    } else chi2 = -1.0;

    /* Get fitted values - switch order to RM, EVPA*/
    arg->coef[0] = (ofloat)gsl_vector_get(arg->solver->x, 1);
    arg->coef[1] = (ofloat)gsl_vector_get(arg->solver->x, 0);
      
    /* Errors wanted? */
    if (arg->doError) {
      /* second argument removes degenerate col/row from Jacobean */
      gsl_multifit_covar (arg->solver->J, 1.0e-8, arg->covar);
      arg->coef[nterm+0] = sqrt(gsl_matrix_get(arg->covar, 1, 1));
      arg->coef[nterm+1] = sqrt(gsl_matrix_get(arg->covar, 0, 0));
      arg->coef[4] = chi2;
    } /* end of get errors */
      
#endif /* HAVE_GSL */ 
} /* end NLRMFit */

/**
 * Make crude estimate of RM, EVPA and populate argument arrays
 * \param arg      NLRMFitArg structure
 *                 Data to be fitted in x,q,u,w
 *                 RM, EVPA in coef
 * \return number of valid data
 */
static olong RMcoarse (NLRMFitArg *arg)
{
  olong i, j, ntest, nvalid, numb;
  ofloat minDL, maxDL, dRM, tRM, bestRM=0.0, fblank = ObitMagicF();
  ofloat bestQ, bestU, sumrQ, sumrU, penFact;
  double aarg, varaarg;
  odouble amb, res, best, test;
 

  /* get EVLA values count valid data*/
  nvalid = 0;
  numb   = 0;
  for (i=0; i<arg->nlamb2; i++) {
    if ((arg->Qobs[i]!=fblank) && (arg->Qweight[i]>0.0) &&
	(arg->Uobs[i]!=fblank) && (arg->Uweight[i]>0.0)) {
      arg->x[numb] = arg->lamb2[i] - arg->refLamb2;
      /* Weight */
      aarg    = fabs(arg->Uobs[i]/arg->Qobs[i]);
      varaarg = aarg*aarg*(arg->Qvar[i]/(arg->Qobs[i]*arg->Qobs[i]) + 
			   arg->Uvar[i]/(arg->Uobs[i]*arg->Uobs[i]));
      arg->w[numb] = (1.0+aarg*aarg)/varaarg;
      /* arg->q[numb] = arg->Qobs[i] * arg->w[numb];
	 arg->u[numb] = arg->Uobs[i] * arg->w[numb];*/
      arg->q[numb] = arg->Qobs[i];
      arg->u[numb] = arg->Uobs[i];
      /* DEBUG
      fprintf (stderr, "%3d x=%8.5f q=%8.5f u=%8.5f w=%8.5g \n ",
	       numb,arg->x[numb], arg->q[numb],arg->u[numb],arg->w[numb]);
      arg->Qweight[i] = arg->Uweight[i] = 1.0/20.0e-6;
      fprintf (stderr,"%3d l2=%g q=%g u=%g p=%f qwt=%g uwt=%g\n",
	       i,arg->lamb2[i]-arg->refLamb2, arg->Qobs[i], arg->Uobs[i],arg->Pobs[i], 
	       arg->Qweight[i], arg->Uweight[i]);*/
     numb++;
    }
  }
  nvalid = numb;
  if (nvalid<=MAX(2,arg->nterm)) return nvalid;  /* enough good data for fit? */
  /* High enough fraction of valid pixels? */
  if ((((ofloat)nvalid)/((ofloat)arg->nlamb2)) < arg->minFrac) return nvalid;
  arg->nvalid = nvalid;

  /* max and min delta lamb2 - assume lamb2 ordered */
  maxDL = fabs(arg->x[0]-arg->x[numb-1]);   /* range of actual delta lamb2 */
  minDL  = 1.0e20; 
  for (i=1; i<numb; i++) minDL = MIN (minDL, fabs(arg->x[i]-arg->x[i-1]));
  /* 
     ambiguity  = pi/min_dlamb2
     resolution = pi/max_dlamb2 = RM which gives 1 turn over  lamb2 range
   */
  amb = G_PI / fabs(minDL);
  res = G_PI / maxDL;
  dRM = 0.05 * res;   /* test RM interval */

  /* DEBUG
  fprintf (stderr,"amb = %f res= %f dRM=%f\n", amb, res, dRM); */
  /* end DEBUG */
  /* Test +/- half ambiguity every dRM */
  ntest = 1 + 0.5*amb/dRM;
  ntest =  MIN (ntest, 1001);   /* Some bounds */
  /* Penalty to downweight solutions away from zero, 
     weight 0.25 at edge of search */
  penFact = 0.25 / (0.5*ntest);
  best = -1.0e20;
  for (i=0; i<ntest; i++) {
    tRM = (i-ntest/2) * dRM;
    /* Loop over data samples - first phases to convert to ref Lamb2 */
    for (j=0; j<numb; j++) arg->wrk1[j]= -2*tRM * arg->x[j];
    /* sines/cosine */
    ObitSinCosVec (numb, arg->wrk1, arg->wrk3, arg->wrk2);
    /* DEBUG
    if ((fabs(tRM-133.0)<0.95*dRM) || (fabs(tRM-133.)<0.95*dRM)) {
      for (j=0; j<numb; j++) {
	test = (arg->q[j]*arg->wrk2[j])*(arg->q[j]*arg->wrk2[j]) + 
	  (arg->u[j]*arg->wrk3[j])*(arg->u[j]*arg->wrk3[j]);
	fprintf (stderr,"   i=%d j=%d t phase=%f obs phase %lf q=%f u=%lf test=%f\n",
		 i, j, arg->wrk1[j], 0.5*atan2(arg->u[j],arg->q[j]),
		 (arg->q[j]*arg->wrk2[j] - arg->u[j]*arg->wrk3[j])*1000, 
		 (arg->q[j]*arg->wrk3[j] + arg->u[j]*arg->wrk2[j])*1000, 
		 test);
      }
    } */
    /* end  DEBUG*/
    /* Find best sum of weighted amplitudes^2 converted to ref lambda^2 */
    test = 0.0; sumrQ = sumrU = 0.0;
    for (j=0; j<numb; j++) {
      sumrQ += arg->w[j]*(arg->q[j]*arg->wrk2[j] - arg->u[j]*arg->wrk3[j]); 
      sumrU += arg->w[j]*(arg->q[j]*arg->wrk3[j] + arg->u[j]*arg->wrk2[j]);
    } 
    test = sumrQ*sumrQ + sumrU*sumrU;
    /* Add penalty */
    test *= (1.0 - penFact*abs(i-ntest/2));
    /* DEBUG 
    fprintf (stderr," i=%d test=%g  tRM %f sum Q=%f sum U=%f pen %f\n",
	     i, test, tRM, sumrQ, sumrU,  penFact*abs(i-ntest/2));*/
  /* end DEBUG */
    if (test>best) {
      best = test;
      bestRM = tRM;
      bestQ = sumrQ;
      bestU = sumrU;
    }
  }

  /* Save */
  arg->coef[0] = bestRM;
  arg->coef[1] = 0.5 * atan2(bestU, bestQ);
  
  /* DEBUG 
  fprintf (stderr,"RM = %f EVPA= %f best=%f dRM=%f\n",
	   arg->coef[1],arg->coef[0], best, dRM); */
  /* end DEBUG */

  return nvalid;
} /* end RMcoarse */

#ifdef HAVE_GSL
/**
 * Function evaluator for spectral fitting solver
 * Evaluates (model-observed) / sigma
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of spectral_terms
 * \param param   Function parameter structure (NLRMFitArg)
 * \param f       Vector of (model-obs)/sigma for data points
 *                in order q, u each datum
 * \return completion code GSL_SUCCESS=OK
 */
static int RMFitFunc (const gsl_vector *x, void *params, 
		      gsl_vector *f)
{
  NLRMFitArg *args = (NLRMFitArg*)params;
  ofloat fblank = ObitMagicF();
  olong nlamb2  = args->nlamb2;
  double func;
  odouble RM, EVPA;
  size_t i, j;
  /* DEBUG
  odouble sum=0.0; */

  /* get model parameters */
  EVPA = gsl_vector_get(x, 0);
  RM   = gsl_vector_get(x, 1);
  /* DEBUG
  fprintf (stderr,"FitFunc RM=%f EVPA=%f\n", RM,EVPA); */

  /* First compute model phases */
  for (j=0; j<nlamb2; j++) 
    args->wrk1[j] = 2*(EVPA + RM * (args->lamb2[j]-args->refLamb2));
  /* then sine/cosine */
  ObitSinCosVec (nlamb2, args->wrk1, args->wrk3, args->wrk2);
  /* Loop over data - Residuals */
  for (i=0; i<nlamb2; i++) {
    /* Q model = args->Pobs[i] * args->wrk2[i]  */
    if ((args->Qobs[i]!=fblank) && (args->Qweight[i]>0.0)) {
      func = (args->Pobs[i] * args->wrk2[i] - args->Qobs[i]) * args->Qweight[i];
      /*  sum += func*func;   DEBUG
      fprintf (stderr,"FitFunc q i=%3ld func=%lg obs=%g mod=%g wt=%g p=%f\n", 
	       2*i, func, args->Uobs[i], args->Pobs[i]*args->wrk3[i],
	       args->Uweight[i],  args->Pobs[i]); */
      gsl_vector_set(f, 2*i, func);  /* Save function residual */
    } else {  /* Invalid data */
      func = 0.0;
      gsl_vector_set(f, 2*i, func);     /* Save function residual */
    }
    /* U model = args->Pobs[i] * args->wrk3[i]  */
    if ((args->Uobs[i]!=fblank) && (args->Uweight[i]>0.0)) {
      func = (args->Pobs[i] * args->wrk3[i] - args->Uobs[i]) * args->Uweight[i];
      /* sum += func*func;   DEBUG
      fprintf (stderr,"FitFunc u i=%3ld func=%lg obs=%g mod=%g wt=%g phs=%f\n", 
	       2*i+1, func, args->Uobs[i], args->Pobs[i]*args->wrk3[i],
	       args->Uweight[i], args->wrk1[i]*28.648); */
      gsl_vector_set(f, 2*i+1, func);  /* Save function residual */
    } else {  /* Invalid data */
      func = 0.0;
      gsl_vector_set(f, 2*i+1, func);     /* Save function residual */
    }
  } /* End loop over data */
  /* DEBUG
  fprintf (stderr,"FitFunc RMS=%g\n", sqrt(sum)); */
  
  return GSL_SUCCESS;
} /*  end RMFitFunc */

/**
 * Jacobian evaluator for spectral fitting solver
 * Evaluates partial derivatives of model wrt each parameter
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of spectral_terms
 * \param param   Function parameter structure (NLRMFitArg)
 * \param J       Jacobian matrix J[data_point, parameter]
 *                in order q, u each datum
 * \return completion code GSL_SUCCESS=OK
 */
static int RMFitJac (const gsl_vector *x, void *params, 
		     gsl_matrix *J)
{
  NLRMFitArg *args = (NLRMFitArg*)params;
  ofloat fblank = ObitMagicF();
  olong nlamb2  = args->nlamb2;
  odouble RM, EVPA;
  double jac;
  size_t i, j;

  /* get model parameters */
  EVPA = gsl_vector_get(x, 0);
  RM   = gsl_vector_get(x, 1);
  /* DEBUG 
  fprintf (stderr,"FitJac RM=%f EVPA=%f\n", RM,EVPA);*/

  /* First compute model phases */
  for (j=0; j<nlamb2; j++) 
    args->wrk1[j] = 2*(EVPA + RM * (args->lamb2[j]-args->refLamb2));
  /* then sine/cosine */
  ObitSinCosVec (nlamb2, args->wrk1, args->wrk3, args->wrk2);
  /* Loop over data - gradients */
  for (i=0; i<nlamb2; i++) {
    j = 0;    /* EVPA */
    /* d Q model/d EVPA = -2 args->Pobs[i] * args->wrk3[i]  */
    if ((args->Qobs[i]!=fblank) && (args->Qweight[i]>0.0)) {
      jac = -2*(args->Pobs[i] * args->wrk3[i]) * args->Qweight[i];
      gsl_matrix_set(J, 2*i, j, jac);  /* Save function gradient */
    } else {  /* Invalid data */
      jac = 0.0;
      gsl_matrix_set(J, 2*i, j, jac);     /* Save function gradient */
    }
    /* d U model/d EVPA = 2 args->Pobs[i] * args->wrk2[i]  */
    if ((args->Uobs[i]!=fblank) && (args->Uweight[i]>0.0)) {
      jac = 2*(args->Pobs[i] * args->wrk2[i]) * args->Uweight[i];
      gsl_matrix_set(J, 2*i+1, j, jac);  /* Save function gradient */
    } else {  /* Invalid data */
      jac = 0.0;
      gsl_matrix_set(J, 2*i+1, j, jac);     /* Save function gradient */
    }
    j = 1;   /* RM  */
    /* d Q model/d RM = -2 args->Pobs[i] * args->wrk3[i] * (args->lamb2[i]-args->refLamb2)  */
    if ((args->Qobs[i]!=fblank) && (args->Qweight[i]>0.0)) {
      jac = -2*(args->Pobs[i] * args->wrk3[i]) * 
	(args->lamb2[i]-args->refLamb2) * args->Qweight[i];
      gsl_matrix_set(J, 2*i, j, jac);  /* Save function gradient */
    } else {  /* Invalid data */
      jac = 0.0;
      gsl_matrix_set(J, 2*i, j, jac);     /* Save function gradient */
    }
    /* d U model/d RM = 2 args->Pobs[i] * args->wrk2[i] * (args->lamb2[i]-args->refLamb2) */
    if ((args->Uobs[i]!=fblank) && (args->Uweight[i]>0.0)) {
      jac = 2*(args->Pobs[i] * args->wrk2[i]) * 
	(args->lamb2[i]-args->refLamb2) * args->Uweight[i];
      gsl_matrix_set(J, 2*i+1, j, jac);  /* Save function gradient */
    } else {  /* Invalid data */
      jac = 0.0;
      gsl_matrix_set(J, 2*i+1, j, jac);     /* Save function gradient */
    }
  } /* End loop over data */
  
  return GSL_SUCCESS;
} /*  end RMFitJac */

/**
 * Function and Jacobian evaluator for spectral fitting solver
 * Function = (model-observed) / sigma
 * Jacobian =  partial derivatives of model wrt each parameter
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of spectral_terms
 * \param param   Function parameter structure (NLRMFitArg)
 * \param f       Vector of (model-obs)/sigma for data points
 *                in order q, u each datum
 * \param J       Jacobian matrix J[data_point, parameter]
 *                in order q, u each datum
 * \return completion code GSL_SUCCESS=OK
 */
static int RMFitFuncJac (const gsl_vector *x, void *params, 
			 gsl_vector *f, gsl_matrix *J)
{
  NLRMFitArg *args = (NLRMFitArg*)params;
  ofloat fblank = ObitMagicF();
  double func, jac;
  olong nlamb2  = args->nlamb2;
  odouble RM, EVPA;
  size_t i, j;
  /* DEBUG
  odouble sum=0.0; */

  /* get model parameters */
  EVPA = gsl_vector_get(x, 0);
  RM   = gsl_vector_get(x, 1);
  /* DEBUG
  fprintf (stderr,"FuncJac RM=%f EVPA=%f\n", RM,EVPA); */

  /* First compute model phases */
  for (j=0; j<nlamb2; j++) 
    args->wrk1[j] = 2*(EVPA + RM * (args->lamb2[j]-args->refLamb2));
  /* then sine/cosine */
  ObitSinCosVec (nlamb2, args->wrk1, args->wrk3, args->wrk2);
  /* Loop over data - gradients */
  for (i=0; i<nlamb2; i++) {
    if ((args->Qobs[i]!=fblank) && (args->Qweight[i]>0.0)) {
      /* Q model   = args->Pobs[i] * args->wrk2[i]  */
      func = (args->Pobs[i] * args->wrk2[i] - args->Qobs[i]) * args->Qweight[i];
      /* sum += func*func;   DEBUG */
      /* DEBUG
      fprintf (stderr,"FuncJac q i=%3ld func=%lg obs=%g mod=%g wt=%g p=%f\n", 
	       2*i, func, args->Qobs[i], args->Pobs[i] * args->wrk2[i],
	       args->Qweight[i], args->Pobs[i]); */
      gsl_vector_set(f, 2*i, func);  /* Save function residual */
      j = 0;    /* EVPA */
      /* d Q model/d EVPA = -2 args->Pobs[i] * args->wrk3[i]  */
      jac = -2*(args->Pobs[i] * args->wrk3[i]) * args->Qweight[i];
      /* DEBUG
      fprintf (stderr,"FuncJac q i=%3ld j=%3ld jac=%lf\n", 2*i, j, jac); */
      gsl_matrix_set(J, 2*i, j, jac);  /* Save function residual */
      j = 1;   /* RM  */
      /* d Q model/d RM = -2 args->Pobs[i] * args->wrk3[i] * (args->lamb2[i]-args->refLamb2)  */
      jac = -2*(args->Pobs[i] * args->wrk3[i]) * 
	(args->lamb2[i]-args->refLamb2) * args->Qweight[i];
      /* DEBUG 
      fprintf (stderr,"FuncJac q i=%3ld j=%3ld jac=%lf\n", 2*i+1, j, jac);*/
      gsl_matrix_set(J, 2*i, j, jac);  /* Save function gradient */
    } else {  /* Invalid data */
      func = 0.0;
      gsl_vector_set(f, 2*i, func);
      j = 0;    /* EVPA */
      jac = 0.0;
      gsl_matrix_set(J, 2*i, j, jac);
      j = 1;   /* RM  */
      gsl_matrix_set(J, 2*i, j, jac);
    }
    if ((args->Uobs[i]!=fblank) && (args->Uweight[i]>0.0)) {
      /* U model   = 2 args->Pobs[i] * args->wrk3[i]  */
      func = (args->Pobs[i] * args->wrk3[i] - args->Uobs[i]) * args->Uweight[i];
      /* sum += func*func;   DEBUG */
      /* DEBUG
      fprintf (stderr,"FuncJac u i=%3ld func=%lg obs=%g mod=%g wt=%g phs=%f\n", 
	       2*i+1, func, args->Uobs[i], args->Pobs[i]*args->wrk3[i],
	       args->Uweight[i], args->wrk1[i]*28.648); */
      gsl_vector_set(f, 2*i+1, func);  /* Save function residual */
      j = 0;    /* EVPA */
      /* d U model/d EVPA = 2 args->Pobs[i] * args->wrk2[i]  */
      jac = 2*(args->Pobs[i] * args->wrk2[i]) * args->Uweight[i];
      /* DEBUG
      fprintf (stderr,"FuncJac u i=%3ld j=%3ld jac=%lf\n", 2*i, j, jac); */
      gsl_matrix_set(J, 2*i+1, j, jac);  /* Save function gradient */
      j = 1;   /* RM  */
      /* d U model/d RM = args->Pobs[i] * args->wrk2[i] * (args->lamb2[i]-args->refLamb2) */
      jac = 2*(args->Pobs[i] * args->wrk2[i]) * 
	(args->lamb2[i]-args->refLamb2) * args->Uweight[i];
      /* DEBUG
      fprintf (stderr,"FuncJac u i=%3ld j=%3ld jac=%lf\n", 2*i+1, j, jac); */
      gsl_matrix_set(J, 2*i+1, j, jac);  /* Save function gradient */      
    } else {  /* Invalid data */
      func = 0.0;
      gsl_vector_set(f, 2*i+1, func);
      j = 0;    /* EVPA */
      jac = 0.0;
      gsl_matrix_set(J, 2*i+1, j, jac);
      j = 1;   /* RM  */
      gsl_matrix_set(J, 2*i+1, j, jac);
    }
  } /* End loop over data */
   /* DEBUG
  fprintf (stderr,"FuncJac RMS=%g\n", sqrt(sum)); */
 
  return GSL_SUCCESS;
} /*  end RMFitFuncJac */
#endif /* HAVE_GSL */ 
