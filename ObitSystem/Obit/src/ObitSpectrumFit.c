/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008-2010                                          */
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

#include "ObitSpectrumFit.h"
#include "ObitThread.h"
#ifdef HAVE_GSL
#include <gsl/gsl_blas.h>
#endif /* HAVE_GSL */ 
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSpectrumFit.c
 * ObitSpectrumFit class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitSpectrumFit";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitSpectrumFitClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSpectrumFitClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitSpectrumFitInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSpectrumFitClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSpectrumFitClassInfoDefFn (gpointer inClass);

/** Private: Do actual fitting */
static void Fitter (ObitSpectrumFit* in, ObitErr *err);

/** Private: Write output image */
static void WriteOutput (ObitSpectrumFit* in, ObitImage *outImage, 
			 ObitErr *err);

#ifdef HAVE_GSL
/** Private: Solver function evaluation */
static int SpecFitFunc (const gsl_vector *x, void *params, 
			gsl_vector *f);

/** Private: Solver Jacobian evaluation */
static int SpecFitJac (const gsl_vector *x, void *params, 
			gsl_matrix *J);

/** Private: Solver function + Jacobian evaluation */
static int SpecFitFuncJac (const gsl_vector *x, void *params, 
			gsl_vector *f, gsl_matrix *J);

/** Private: Solver function evaluation - broken power law */
static int SpecFitFuncBP (const gsl_vector *x, void *params, 
			  gsl_vector *f);

/** Private: Solver Jacobian evaluation - broken power law */
static int SpecFitJacBP (const gsl_vector *x, void *params, 
			 gsl_matrix *J);

/** Private: Solver function + Jacobian evaluation - broken power law */
static int SpecFitFuncJacBP (const gsl_vector *x, void *params, 
			     gsl_vector *f, gsl_matrix *J);
#endif /* HAVE_GSL */ 

/** Private: Threaded fitting */
static gpointer ThreadNLFit (gpointer arg);

/*---------------Private structures----------------*/
/* FT threaded function argument */
typedef struct {
  /** ObitSpectrumFit object */
  ObitSpectrumFit *in;
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
  /** number of terms being fitted */
  olong        fitTerm;
  /** number of frequencies  */
  olong        nfreq;
  /** Array of Nu per frequency point - broken power law */
  ofloat *nu;
  /** Array of log (Nu/Nu_0) per frequency point */
  ofloat *logNuOnu0;
  /** Array of weights (1/RMS**2) per inFArrays (nfreq) */
  ofloat *weight;
  /** Array of 1/RMS per inFArrays (nfreq) */
  ofloat *isigma;
  /** Array of pixel values being fitted per frequency point */
  ofloat *obs;
  /** maximum iteration  */
  olong        maxIter;
  /** acceptable iteration delta, rel and abs. */
  odouble minDelta;
  /** max acceptable normalized chi squares */
  ofloat maxChiSq;
  /** Reference frequency */
  ofloat refFreq;
  /** Vector of guess/fitted coefficients, optional errors */
  ofloat *coef;
  /** Chi squared of fit */
  ofloat ChiSq;
  /** Do error analysis? */
  gboolean doError;
  /** Do broken power law ? */
  gboolean doBrokePow;
  /** Do Primary beam correction? */
  gboolean doPBCorr;
#ifdef HAVE_GSL
  /** Fitting solver for two terms */
  gsl_multifit_fdfsolver *solver2;
  /** Fitting solver for three terms */
  gsl_multifit_fdfsolver *solver3;
  /** Fitting solver for four terms */
  gsl_multifit_fdfsolver *solver4;
  /** Fitting solver for five terms */
  gsl_multifit_fdfsolver *solver5;
  /** Fitting solver function structure */
  gsl_multifit_function_fdf *funcStruc;
  /** Fitting work vector of length 2 */
  gsl_vector *work2;
  /** Fitting work vector of length 3 */
  gsl_vector *work3;
  /** Fitting work vector of length 4 */
  gsl_vector *work4;
  /** Fitting work vector of length 5 */
  gsl_vector *work5;
  /** Covariance matrix for two terms */
  gsl_matrix *covar2;
  /** Covariance matrix for three terms */
  gsl_matrix *covar3;
  /** Covariance matrix for four terms */
  gsl_matrix *covar4;
  /** Covariance matrix for five terms */
  gsl_matrix *covar5;
#endif /* HAVE_GSL */ 
} NLFitArg;

/** Private: Actual fitting */
static void NLFit (NLFitArg *arg);

/** Private: Actual fitting broken power */
static void NLFitBP (NLFitArg *arg);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSpectrumFit* newObitSpectrumFit (gchar* name)
{
  ObitSpectrumFit* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSpectrumFitClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSpectrumFit));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSpectrumFitInit((gpointer)out);

 return out;
} /* end newObitSpectrumFit */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSpectrumFitGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSpectrumFitClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitSpectrumFitGetClass */

/**
 * Make a deep copy of an ObitSpectrumFit.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitSpectrumFit* ObitSpectrumFitCopy  (ObitSpectrumFit *in, ObitSpectrumFit *out, ObitErr *err)
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
    out = newObitSpectrumFit(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nfreq = in->nfreq;
  out->nterm = in->nterm;

  /* Arrays */
  if (out->RMS) g_free(out->RMS);
  out->RMS = g_malloc0(out->nfreq*sizeof(ofloat));
  if (in->RMS)
    for (i=0; i<out->nfreq; i++) out->RMS[i] = in->RMS[i];
  if (out->calFract) g_free(out->calFract);
  out->calFract = g_malloc0(out->nfreq*sizeof(ofloat));
  if (in->calFract)
    if (out->calFract) g_free(out->calFract);
    
  /* reference this class members */
  if (out->outDesc) out->outDesc = ObitImageDescUnref(out->outDesc);
  if (in->outDesc)  out->outDesc = ObitImageDescRef(in->outDesc);
  if (out->inFArrays) {
    for (i=0; i<in->nfreq; i++) out->inFArrays[i] = ObitFArrayUnref(out->inFArrays[i]);
  }
  if (in->inFArrays) {
    for (i=0; i<in->nfreq; i++) out->inFArrays[i] = ObitFArrayRef(in->inFArrays[i]);
  }

  if (out->BeamShapes) {
    for (i=0; i<in->nfreq; i++) out->BeamShapes[i] = ObitBeamShapeUnref(out->BeamShapes[i]);
  }
  if (in->BeamShapes) {
    for (i=0; i<in->nfreq; i++) out->BeamShapes[i] = ObitBeamShapeRef(in->BeamShapes[i]);
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
} /* end ObitSpectrumFitCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an SpectrumFit similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitSpectrumFitClone  (ObitSpectrumFit *in, ObitSpectrumFit *out, ObitErr *err)
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
  out->nfreq = in->nfreq;
  out->nterm = in->nterm;
  
  /* Arrays */
  if (out->RMS) g_free(out->RMS);
  out->RMS = g_malloc0(out->nfreq*sizeof(ofloat));
  if (in->RMS)
    for (i=0; i<out->nfreq; i++) out->RMS[i] = in->RMS[i];
  if (out->calFract) g_free(out->calFract);
  out->calFract = g_malloc0(out->nfreq*sizeof(ofloat));
  if (in->calFract)
    if (out->calFract) g_free(out->calFract);
    

  /* reference this class members */
  if (out->outDesc) out->outDesc = ObitImageDescUnref(out->outDesc);
  if (in->outDesc)  out->outDesc = ObitImageDescRef(in->outDesc);
  if (out->inFArrays) {
    for (i=0; i<in->nfreq; i++) out->inFArrays[i] = ObitFArrayUnref(out->inFArrays[i]);
  }
  if (in->inFArrays) {
    for (i=0; i<in->nfreq; i++) out->inFArrays[i] = ObitFArrayRef(in->inFArrays[i]);
  }

  if (out->BeamShapes) {
    for (i=0; i<in->nfreq; i++) out->BeamShapes[i] = ObitBeamShapeUnref(out->BeamShapes[i]);
  }
  if (in->BeamShapes) {
    for (i=0; i<in->nfreq; i++) out->BeamShapes[i] = ObitBeamShapeRef(in->BeamShapes[i]);
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

} /* end ObitSpectrumFitClone */

/**
 * Creates an ObitSpectrumFit 
 * \param name   An optional name for the object.
 * \param nterm  Number of coefficients of powers of log(nu)
 * \return the new object.
 */
ObitSpectrumFit* ObitSpectrumFitCreate (gchar* name, olong nterm)
{
  ObitSpectrumFit* out;

  /* Create basic structure */
  out = newObitSpectrumFit (name);
  out->nterm = nterm;

  return out;
} /* end ObitSpectrumFitCreate */

/**
 * Fit spectra to an image cube pixels.
 * The third axis of the output image will be "SPECLOGF" to indicate that then
 * planes are spectral fit parameters.  The "CRVAL" on this axis will be the reference 
 * Frequency for the fitting.
 * Item "NTERM" is added to the output image descriptor to give the maximum number 
 * Item ""BROKENPO" is added to the output image descriptor if fit is a broken power law
 * of terms fitted
 * \param in       Spectral fitting object
 *                 Potential parameters on in->info:
 * \li "refFreq" OBIT_double scalar Reference frequency for fit [def ref for inImage]
 * \li "maxChi2" OBIT_float scalar Max. Chi Sq for accepting a partial spectrum [def 1.5]
 * \li "doError" OBIT_boolean scalar If true do error analysis [def False]
 * \li "doPBCor" OBIT_boolean scalar If true do primary beam correction. [def False]
 * \li "doBrokePow" OBIT_boolean scalar If true do broken power law (3 terms). [def False]
 * \li "calFract" OBIT_float (?,1,1) Calibration error as fraction of flux
 *              One per frequency or one for all, def 1.0e-5
 * \li "PBmin"  OBIT_float (?,1,1) Minimum beam gain correction
 *              One per frequency or one for all, def 0.05, 1.0 => no gain corrections
 * \li "antSize" OBIT_float (?,1,1) Antenna diameter (m) for gain corr, 
 *              One per frequency or one for all, def 25.0
 *
 * \param inImage  Image cube to be fitted
 * \param outImage Image cube with fitted spectra.
 *                 Should be defined but not created.
 *                 Planes 1->nterm are coefficients per pixel
 *                 if doError:
 *                 Planes nterm+1->2*nterm are uncertainties in coefficients
 *                 Plane 2*nterm+1 = Chi squared of fit
 * \param err      Obit error stack object.
 */
void ObitSpectrumFitCube (ObitSpectrumFit* in, ObitImage *inImage, 
			  ObitImage *outImage, ObitErr *err)
{
  olong i, iplane, nOut;
  olong naxis[2];
  ObitIOSize IOBy;
  ObitInfoType type;
  ObitIOCode retCode;
  union ObitInfoListEquiv InfoReal; 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1}, PBdim[MAXINFOELEMDIM], ASdim[MAXINFOELEMDIM];
  ofloat *calFract, *PBmin, *antSize, pbmin, antsize;
  gboolean doGain;
  gchar *today=NULL, *SPECLOGF = "SPECLOGF";
  gchar *routine = "ObitSpectrumFitCube";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert(ObitImageIsA(inImage));
  g_assert(ObitImageIsA(outImage));

  /* Control parameters */
  /* Min Chi^2 for fit */
  InfoReal.flt = 1.5; type = OBIT_float;
  in->maxChi2 = 1.5;
  ObitInfoListGetTest(in->info, "maxChi2", &type, dim, &InfoReal);
  if (type==OBIT_float) in->maxChi2 = InfoReal.flt;
  else if (type==OBIT_double) in->maxChi2 = (ofloat)InfoReal.dbl;

  /* Want Error analysis? */
  InfoReal.itg = (olong)FALSE; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doError", &type, dim, &InfoReal);
  in->doError = InfoReal.itg;

  /* Want primary beam correction? */
  InfoReal.itg = (olong)FALSE; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doPBCor", &type, dim, &InfoReal);
  in->doPBCorr = InfoReal.itg;

  /* Want Broken power law ? */
  InfoReal.itg = (olong)FALSE; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doBrokePow", &type, dim, &InfoReal);
  in->doBrokePow = InfoReal.itg;

  /* Min PB gain */
  ObitInfoListGetP(in->info, "PBmin", &type, PBdim, (gpointer)&PBmin);
  /* Antenna diameter */
  ObitInfoListGetP(in->info, "antSize", &type, ASdim, (gpointer)&antSize);

  /* Open input image to get info */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (inImage->info, "IOBy", OBIT_long, dim, &IOBy);
  inImage->extBuffer = TRUE;   /* Using inFArrays as I/O buffer */
  retCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inImage->name);

  /* Get Reference frequency , default input ref. freq. */
  InfoReal.dbl = inImage->myDesc->crval[inImage->myDesc->jlocf]; 
  type = OBIT_double;
  in->refFreq = InfoReal.dbl;
  ObitInfoListGetTest(in->info, "refFreq", &type, dim, &InfoReal);
  if (type==OBIT_float) in->refFreq = InfoReal.flt;
  else if (type==OBIT_double) in->refFreq = (ofloat)InfoReal.dbl;
  
  /* Determine number of frequency planes and initialize in */
  in->nfreq = inImage->myDesc->inaxes[inImage->myDesc->jlocf];
  in->BeamShapes = g_malloc0(in->nfreq*sizeof(ObitBeamShape*));
  for (i=0; i<in->nfreq; i++) in->BeamShapes[i] = NULL;
  in->RMS        = g_malloc0(in->nfreq*sizeof(ofloat));
  in->calFract   = g_malloc0(in->nfreq*sizeof(ofloat));
  in->inFArrays  = g_malloc0(in->nfreq*sizeof(ObitFArray*));
  in->freqs      = g_malloc0(in->nfreq*sizeof(odouble));

  /* How many output planes? */
  if (in->doError) nOut = 1+in->nterm*2;
  else nOut = in->nterm;

   /* Define term arrays */
  in->outFArrays = g_malloc0((nOut)*sizeof(ObitFArray*));
  for (i=0; i<nOut; i++) in->outFArrays[i] = NULL;

 /* Image size */
  in->nx = inImage->myDesc->inaxes[0];
  in->ny = inImage->myDesc->inaxes[1];
  naxis[0] = (olong)in->nx;  naxis[1] = (olong)in->ny; 
  for (i=0; i<in->nfreq; i++) in->inFArrays[i]  = ObitFArrayCreate (NULL, 2, naxis);
  for (i=0; i<nOut; i++)      in->outFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);

  /* Calibration error */
  for (i=0; i<in->nfreq; i++) in->calFract[i] = 1.0e-5; /* default */
  if (ObitInfoListGetP(in->info, "calFract", &type, dim, (gpointer)&calFract)) {
    if (dim[0]>=in->nfreq) for (i=0; i<in->nfreq; i++) in->calFract[i] = calFract[i];
    else for (i=0; i<in->nfreq; i++) in->calFract[i] = calFract[0];
  }

  /* Output Image descriptor */
  outImage->myDesc = ObitImageDescCopy (inImage->myDesc, in->outDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Change third axis to type "SPECLOGF" and leave the reference frequency
     as the "CRVAL" */
  outImage->myDesc->inaxes[outImage->myDesc->jlocf] =  nOut;
  outImage->myDesc->crval[outImage->myDesc->jlocf]  =  in->refFreq;
  outImage->myDesc->crpix[outImage->myDesc->jlocf]  =  1.0;
  outImage->myDesc->cdelt[outImage->myDesc->jlocf]  =  1.0;
  strncpy (outImage->myDesc->ctype[outImage->myDesc->jlocf], SPECLOGF, IMLEN_KEYWORD);
  outImage->myDesc->bitpix = -32;  /* Float it */

  /* Creation date today */
  today = ObitToday();
  strncpy (outImage->myDesc->date, today, IMLEN_VALUE);
  if (today) g_free(today);

  /* Copy of output descriptor to in */
  in->outDesc = ObitImageDescCopy (outImage->myDesc, in->outDesc, err);

  /* Loop reading planes */
  for (iplane=0; iplane<in->nfreq; iplane++) {
    retCode = ObitImageRead (inImage, in->inFArrays[iplane]->array, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) Obit_traceback_msg (err, routine, inImage->name);

    /* Get BeamShape */
    antsize = 25.0;
    if (antSize) antsize = antSize[0];
    if (antSize && (ASdim[0]>=in->nfreq)) antsize = antSize[i];
    pbmin   = 0.05;
    if (PBmin) pbmin = PBmin[0];
    if (PBmin && (PBdim[0]>=in->nfreq)) pbmin = PBmin[i];
    doGain = pbmin<0.999;
    in->BeamShapes[iplane] = ObitBeamShapeCreate ("BS", inImage, pbmin, antsize, doGain);

    /* Plane RMS */
    in->RMS[iplane] = ObitFArrayRMS(in->inFArrays[iplane]);
    /* Frequency */
    in->freqs[iplane] = inImage->myDesc->crval[inImage->myDesc->jlocf] + 
      inImage->myDesc->cdelt[inImage->myDesc->jlocf] * 
      (inImage->myDesc->plane - inImage->myDesc->crpix[inImage->myDesc->jlocf]);
    
  } /* end loop reading planes */

  /* Close input */
  retCode = ObitImageClose (inImage, err);
  inImage->extBuffer = FALSE;   /* May need I/O buffer later */
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inImage->name);

  /* Average Frequency - no
  in->refFreq = 0.0;
  for (i=0; i<in->nfreq; i++) in->refFreq += in->freqs[i];
  in->refFreq /= in->nfreq;*/

  /* Do actual fitting */
  Fitter (in, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Update output header to reference Frequency */
  outImage->myDesc->crval[outImage->myDesc->jlocf] = in->refFreq;
  in->outDesc->crval[in->outDesc->jlocf] = in->refFreq;

  /* Write output */
  WriteOutput(in, outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

} /* end ObitSpectrumFitCube */

/**
 * Fit spectra to an array of images
 * The third axis of the output image will be "SPECLOGF" to indicate that then
 * planes are spectral fit parameters.  The "CRVAL" on this axis will be the reference 
 * Frequency for the fitting.
 * Item "NTERM" is added to the output image descriptor to give the maximum number 
 * of terms fitted
 * Item ""BROKENPO" is added to the output image descriptor if fit is a broken power law
 * \param in       Spectral fitting object
 * \li "refFreq" OBIT_double scalar Reference frequency for fit [def average of inputs]
 * \li "maxChi2" OBIT_float scalar Max. Chi Sq for accepting a partial spectrum [def 1.5]
 * \li "doError" OBIT_boolean scalar If true do error analysis [def False]
 * \li "doPBCor" OBIT_boolean scalar If true do primary beam correction.[def False]
 * \li "doBrokePow" OBIT_boolean scalar If true do broken power law (3 terms). [def False]
 * \li "calFract" OBIT_float (?,1,1) Calibration error as fraction of flux
 *              One per frequency or one for all, def 0.05
 * \li "PBmin"  OBIT_float (?,1,1) Minimum beam gain correction
 *              One per frequency or one for all, def 0.05, 1.0 => no gain corrections
 * \li "antSize" OBIT_float (?,1,1) Antenna diameter (m) for gain corr, 
 *              One per frequency or one for all, def 25.0
 * \param nimage   Number of entries in imArr
 * \param imArr    Array of images to be fitted
 * \param outImage Image cube with fitted spectra.
 *                 Should be defined but not created.
 *                 Planes 1->nterm are coefficients per pixel
 *                 if doError:
 *                 Planes nterm+1->2*nterm are uncertainties in coefficients
 *                 Plane 2*nterm+1 = Chi squared of fit
 * \param err      Obit error stack object.
 */
void ObitSpectrumFitImArr (ObitSpectrumFit* in, olong nimage, ObitImage **imArr, 
			   ObitImage *outImage, ObitErr *err)
{
  olong i, iplane, nOut;
  olong naxis[2];
  ObitIOSize IOBy;
  ObitInfoType type;
  ObitIOCode retCode;
  union ObitInfoListEquiv InfoReal; 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1}, PBdim[MAXINFOELEMDIM], ASdim[MAXINFOELEMDIM];
  ofloat *calFract=NULL, *PBmin=NULL, *antSize=NULL, pbmin, antsize, ipixel[2], opixel[2];
  gboolean doGain, bad;
  odouble avgFreq;
  gchar *today=NULL, *SPECLOGF = "SPECLOGF";
  gchar *routine = "ObitSpectrumFitCube";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert(ObitImageIsA(outImage));

  /* Control parameters */
  /* Max Chi^2 for fit */
   InfoReal.flt = 1.5; type = OBIT_float;
  in->maxChi2 = 1.5;
  ObitInfoListGetTest(in->info, "maxChi2", &type, dim, &InfoReal);
  if (type==OBIT_float) in->maxChi2 = InfoReal.flt;
  else if (type==OBIT_double) in->maxChi2 = (ofloat)InfoReal.dbl;

  /* Want Error analysis? */
  InfoReal.itg = (olong)FALSE; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doError", &type, dim, &InfoReal);
  in->doError = InfoReal.itg;

  /* Want primary beam correction? */
  InfoReal.itg = (olong)FALSE; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doPBCor", &type, dim, &InfoReal);
  in->doPBCorr = InfoReal.itg;

  /* Want Broken power law ? */
  InfoReal.itg = (olong)FALSE; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doBrokePow", &type, dim, &InfoReal);
  in->doBrokePow = InfoReal.itg;

  /* Min PB gain */
  ObitInfoListGetP(in->info, "PBmin", &type, PBdim, (gpointer)&PBmin);
  /* Antenna diameter */
  ObitInfoListGetP(in->info, "antSize", &type, ASdim, (gpointer)&antSize);

  /* Determine number of frequency planes and initialize in */
  in->nfreq = nimage;
  in->BeamShapes = g_malloc0(in->nfreq*sizeof(ObitBeamShape*));
  for (i=0; i<in->nfreq; i++) in->BeamShapes[i] = NULL;
  in->RMS        = g_malloc0(in->nfreq*sizeof(ofloat));
  in->calFract   = g_malloc0(in->nfreq*sizeof(ofloat));
  in->inFArrays  = g_malloc0(in->nfreq*sizeof(ObitFArray*));
  in->freqs      = g_malloc0(in->nfreq*sizeof(odouble));

  /* Calibration error */
  for (i=0; i<in->nfreq; i++) in->calFract[i] = 0.05; /* default */
  if (ObitInfoListGetP(in->info, "calFract", &type, dim, (gpointer)&calFract)) {
    if (dim[0]>=in->nfreq) for (i=0; i<in->nfreq; i++) in->calFract[i] = calFract[i];
    else for (i=0; i<in->nfreq; i++) in->calFract[i] = calFract[0];
  }
      
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
    ObitInfoListAlwaysPut (imArr[iplane]->info, "IOBy", OBIT_long, dim, &IOBy);
    imArr[iplane]->extBuffer = TRUE;   /* Using inFArrays as I/O buffer */
    retCode = ObitImageOpen (imArr[iplane], OBIT_IO_ReadOnly, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) 
      Obit_traceback_msg (err, routine, imArr[iplane]->name);
    
    /* On first image initialize */
    if (iplane==0) {
      /* Image size */
      in->nx = imArr[iplane]->myDesc->inaxes[0];
      in->ny = imArr[iplane]->myDesc->inaxes[1];
      naxis[0] = (olong)in->nx;  naxis[1] = (olong)in->ny; 
      for (i=0; i<in->nfreq; i++) in->inFArrays[i]  = ObitFArrayCreate (NULL, 2, naxis);
      for (i=0; i<nOut; i++)      in->outFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);
      
      /* Output Image descriptor */
      outImage->myDesc = ObitImageDescCopy (imArr[iplane]->myDesc, in->outDesc, err);
      if (err->error) Obit_traceback_msg (err, routine, imArr[iplane]->name);

      /* Change third axis to type "SPECLOGF" and leave the reference frequency
	 as the "CRVAL" */
      outImage->myDesc->inaxes[outImage->myDesc->jlocf] =  nOut;
      outImage->myDesc->crpix[outImage->myDesc->jlocf]  =  1.0;
      outImage->myDesc->cdelt[outImage->myDesc->jlocf]  =  1.0;
      strncpy (outImage->myDesc->ctype[outImage->myDesc->jlocf], SPECLOGF, IMLEN_KEYWORD);
      outImage->myDesc->bitpix = -32;  /* Float it */

      /* Creation date today */
      today = ObitToday();
      strncpy (outImage->myDesc->date, today, IMLEN_VALUE);
      if (today) g_free(today);

      /* Copy of output descriptor to in */
      in->outDesc = ObitImageDescCopy (outImage->myDesc, in->outDesc, err);

    } else { /* On subsequent images check for consistency */
      /* Check size of planes */
      Obit_return_if_fail(((imArr[iplane]->myDesc->inaxes[0]==in->outDesc->inaxes[0]) && 
			   (imArr[iplane]->myDesc->inaxes[1]==in->outDesc->inaxes[1])), err,
			  "%s: Image planes incompatible  %d!= %d or  %d!= %d", 
			  routine, imArr[iplane]->myDesc->inaxes[0], in->outDesc->inaxes[0], 
			  imArr[iplane]->myDesc->inaxes[1], in->outDesc->inaxes[1]) ;

      /* Check alignment of pixels */
      ipixel[0] = 1.0; ipixel[1] = 1.0;
      bad = !ObitImageDescCvtPixel (imArr[iplane]->myDesc, in->outDesc, ipixel, opixel, err);
      if (err->error) Obit_traceback_msg (err, routine, imArr[iplane]->name);
      Obit_return_if_fail(!bad, err,
			  "%s: Image planes incompatible", routine);
      Obit_return_if_fail(((fabs(ipixel[0]-opixel[0])<0.01) &&
			   (fabs(ipixel[1]-opixel[1])<0.01)), err,
			  "%s: Image pixels not aligned %f!=%f or %f!=%f", 
			  routine, ipixel[0], opixel[0], ipixel[1], opixel[1]) ;
     } /* end consistency check */
    
    /* Read plane */
    retCode = ObitImageRead (imArr[iplane], in->inFArrays[iplane]->array, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) Obit_traceback_msg (err, routine, imArr[iplane]->name);
    
    /* Get BeamShape */
    antsize = 25.0;
    if (antSize) antsize = antSize[0];
    if (antSize && (ASdim[0]>=in->nfreq)) antsize = antSize[iplane];
    pbmin   = 0.05;
    if (PBmin) pbmin = PBmin[0];
    if (PBmin && (PBdim[0]>=in->nfreq)) pbmin = PBmin[iplane];
    doGain = pbmin<0.999;
    in->BeamShapes[iplane] = ObitBeamShapeCreate ("BS", imArr[iplane], pbmin, antsize, doGain);
    
    /* Plane RMS */
    in->RMS[iplane] = ObitFArrayRMS(in->inFArrays[iplane]);

    /* Frequency */
    in->freqs[iplane] = imArr[iplane]->myDesc->crval[imArr[iplane]->myDesc->jlocf] + 
      imArr[iplane]->myDesc->cdelt[imArr[iplane]->myDesc->jlocf] * 
      (imArr[iplane]->myDesc->plane - imArr[iplane]->myDesc->crpix[imArr[iplane]->myDesc->jlocf]);
    
    /* Close input */
    retCode = ObitImageClose (imArr[iplane], err);
    imArr[iplane]->extBuffer = FALSE;   /* May need I/O buffer later */
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) 
      Obit_traceback_msg (err, routine, imArr[iplane]->name);
  } /* end loop over input images */

  /* Average Frequency */
  avgFreq = 0.0;
  for (i=0; i<in->nfreq; i++) avgFreq += in->freqs[i];
  avgFreq /= in->nfreq;

  /* Get Reference frequency , default avg. input ref. freq. */
  InfoReal.dbl = avgFreq; type = OBIT_double;
  in->refFreq = InfoReal.dbl;
  ObitInfoListGetTest(in->info, "refFreq", &type, dim, &InfoReal);
  if (type==OBIT_float) in->refFreq = InfoReal.flt;
  else if (type==OBIT_double) in->refFreq = (ofloat)InfoReal.dbl;
  
  /* Update output header to reference Frequency */
  outImage->myDesc->crval[outImage->myDesc->jlocf] = in->refFreq;
  in->outDesc->crval[in->outDesc->jlocf] = in->refFreq;

  /* Do actual fitting */
  Fitter (in, err);
  if (err->error) Obit_traceback_msg (err, routine, imArr[0]->name);

  /* Write output */
  WriteOutput(in, outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

} /* end ObitSpectrumFitImArr */

/**
 * Evaluate spectrum at a given frequency
 * \param in       Spectral fitting object
 * \param inImage  Spectral coefficient image
 *                 Planes 1->nterm are coefficients per pixel
 *                 Planes nterm+1->2*nterm are uncertainties in coefficients
 *                 Plane 2*nterm+1 = Chi squared of fit
 * \param outFreq  Output Frequency in Hz
 * \param outImage Image to write, must be defined but not yet exist.
 * \param err      Obit error stack object.
 */
void ObitSpectrumFitEval (ObitSpectrumFit* in, ObitImage *inImage, 
			  odouble outFreq, ObitImage *outImage, ObitErr *err)
{
  olong ix, iy, i, indx, nx, ny, nterm, jlocspec;
  olong naxis[2];
  ObitIOSize IOBy;
  ObitIOCode retCode;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitFArray **inArrays=NULL, *outArray=NULL;
  odouble refFreq, arg, aarg, sum;
  gchar *today=NULL, *SPECLOGF = "SPECLOGF";
  gchar *routine = "ObitSpectrumFitEval";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert(ObitImageIsA(inImage));
  g_assert(ObitImageIsA(outImage));

  /* Open input image to get info */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (inImage->info, "IOBy", OBIT_long, dim, &IOBy);
  inImage->extBuffer = TRUE;   /* Using inFArrays as I/O buffer */
  retCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inImage->name);

  /* Find "SPECLOGF" axis */
  jlocspec = -1;
  for (i=0; i<inImage->myDesc->naxis; i++) {
    if (!strncmp (inImage->myDesc->ctype[i], SPECLOGF, 8)) jlocspec = i;
  }
  Obit_return_if_fail((jlocspec>=0), err, 
		      "%s: No %s axis on %s", routine, SPECLOGF, inImage->name);
  

  /* Reference frequency is on inImage */
  refFreq = inImage->myDesc->crval[jlocspec];

  /* Output Image descriptor */
  outImage->myDesc = ObitImageDescCopy (inImage->myDesc,  outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  /* Frequency axis is actually fitted spectrum parameters and errors */
  outImage->myDesc->inaxes[outImage->myDesc->jlocf] = 1;
  outImage->myDesc->crval[outImage->myDesc->jlocf] = outFreq;

  /* Creation date today */
  today = ObitToday();
  strncpy (outImage->myDesc->date, today, IMLEN_VALUE);
  if (today) g_free(today);
  outImage->myDesc->bitpix = -32;  /* Float it */

  /* Determine number of frequency terms and initialize storage arrays 
   get nterm from input image descriptor with default number of pixels on
   SPECLOGF" axis */
  nterm = inImage->myDesc->inaxes[jlocspec];
  ObitInfoListGetTest (outImage->myDesc->info, "NTERM", &type, dim, &nterm);

  /* Image size */
  nx = inImage->myDesc->inaxes[0];
  ny = inImage->myDesc->inaxes[1];
  naxis[0] = (olong)nx;  naxis[1] = (olong)ny; 
  for (i=0; i<nterm; i++) inArrays[i]     = ObitFArrayCreate (NULL, 2, naxis);
  outArray = ObitFArrayCreate (NULL, 2, naxis);

  /* Loop reading planes */
  for (i=0; i<nterm; i++) {
    retCode = ObitImageRead (inImage, inArrays[i]->array, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;
  } /* end loop over planes */

  /* Loop over pixels evaluating spectrum */
  arg = log(outFreq/refFreq);
  /* Loop over pixels in Y */
  indx = -1;
  for (iy=0; iy<ny; iy++) {
    /* Loop over pixels in X */
    for (ix=0; ix<nx; ix++) {
      indx ++;
      sum = 0.0;
      aarg = arg;
      /* Sum polynomial in log(nv/nvRef) */
      for (i=1; i<nterm; i++) {
	sum += inArrays[i]->array[indx]*aarg;
	aarg *= arg;
      }
      /* Add value at reference frequency */
      sum = exp(sum) * inArrays[0]->array[indx];
      outArray->array[indx] = sum;
    } /* end loop in X */
  } /* end loop in y */

  /* Write output */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (outImage->info, "IOBy", OBIT_long, dim, &IOBy);
  inImage->extBuffer = TRUE;   /* Using outFArrays as I/O buffer */
  retCode = ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Write */
  retCode = ObitImageWrite (outImage, outArray->array, err);
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Close output */
  retCode = ObitImageClose (outImage, err);
  outImage->extBuffer = FALSE;   /* May need I/O buffer later */
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Cleanup */
 cleanup:
  if (inArrays) {
    for (i=0; i<nterm; i++) inArrays[i] = ObitFArrayUnref (inArrays[i]);
    g_free(inArrays);
  }
  if (outArray) outArray = ObitFArrayUnref (outArray);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
} /* end ObitSpectrumFitEval */

/**
 * Fit single spectrum to flux measurements
 * Without GSL, only the average intensity is determined (no spectrum)
 * \param nfreq    Number of entries in freq, flux, sigma
 * \param nterm    Number of coefficients of powers of log(nu) to fit
 * \param refFreq  Reference frequency (Hz)
 * \param freq     Array of Frequencies (Hz)
 * \param flux     Array of fluxes (Jy) same dim as freq
 * \param sigma    Array of errors (Jy) same dim as freq
 * \param doBrokePow  TRUE if a broken power law fit is desired (3 terms)
 * \param err      Obit error stack object.
 * \return  Array of fitter parameters, errors for each and Chi Squares of fit
 *          Initial terms are in Jy, other in log.
 */
ofloat* ObitSpectrumFitSingle (olong nfreq, olong nterm, odouble refFreq, 
			       odouble *freq, ofloat *flux, ofloat *sigma, 
			       gboolean doBrokePow, ObitErr *err)
{
  ofloat *out = NULL;
  olong i, j;
  NLFitArg *arg=NULL;
  gchar *routine = "ObitSpectrumFitSingle";
  /* GSL implementation */
#ifdef HAVE_GSL
  const gsl_multifit_fdfsolver_type *T=NULL;
#endif /* HAVE_GSL */ 

  if (err->error) return out;

  /* Warn if too many terms asked for */
  if (nterm>5) {
      Obit_log_error(err, OBIT_InfoWarn, 
		    "%s: Asked for %d terms, will limit to 5", routine, nterm);
  }

  /* Warn if incorrect number of terms */
  if (doBrokePow && (nterm!=3)) {
      Obit_log_error(err, OBIT_InfoWarn, 
		    "%s: Must have three terms for broken power law, not %d", 
		     routine, nterm);
  }

  /* Warn if ref. freq <= 0, set to 1.0e9 */
  if (refFreq<=0.0) {
    refFreq = 1.0e9;
    Obit_log_error(err, OBIT_InfoWarn, 
		   "%s: Setting reference Frequency to 1 GHz", routine);
  }

  /* Create function argument */
  arg = g_malloc(sizeof(NLFitArg));
  arg->in             = NULL;     /* Not needed here */
  arg->nfreq          = nfreq;
  arg->nterm          = nterm;
  arg->doError        = TRUE;
  arg->doBrokePow     = doBrokePow;
  arg->doPBCorr       = FALSE;
  arg->maxIter        = 100;
  arg->minDelta       = 1.0e-5;  /* Min step size */
  arg->maxChiSq       = 1.5;     /* max acceptable normalized chi squares */
  arg->refFreq        = refFreq; /* Reference Frequency */
  arg->weight         = g_malloc0(arg->nfreq*sizeof(ofloat));
  arg->isigma         = g_malloc0(arg->nfreq*sizeof(ofloat));
  arg->obs            = g_malloc0(arg->nfreq*sizeof(ofloat));
  arg->nu             = g_malloc0(arg->nfreq*sizeof(ofloat));
  arg->logNuOnu0      = g_malloc0(arg->nfreq*sizeof(ofloat));
  arg->coef           = g_malloc0(2*arg->nterm*sizeof(ofloat));
  for (i=0; i<nfreq; i++) {
    arg->isigma[i]    = 1.0 / sigma[i];
    arg->weight[i]    = arg->isigma[i]*arg->isigma[i];
    arg->obs[i]       = flux[i];
    arg->nu[i]        = freq[i];
    arg->logNuOnu0[i] = log(freq[i]/refFreq);
  }

  /* GSL implementation */
#ifdef HAVE_GSL
  arg->solver2 =  arg->solver3 =  arg->solver4 =  arg->solver5 = NULL;
  arg->covar2  =  arg->covar3  =  arg->covar4  =  arg->covar5  = NULL;
  arg->work2   =  arg->work3   =  arg->work4   =  arg->work5   = NULL;
  /* Setup solvers */
  T = gsl_multifit_fdfsolver_lmder;
  if (arg->nterm>=2)
    arg->solver2 = gsl_multifit_fdfsolver_alloc(T, arg->nfreq, 2);
  if (arg->nterm>=3)
    arg->solver3 = gsl_multifit_fdfsolver_alloc(T, arg->nfreq, 3);
  if (arg->nterm>=4)
    arg->solver4 = gsl_multifit_fdfsolver_alloc(T, arg->nfreq, 4);
  if (arg->nterm>=5)
    arg->solver5 = gsl_multifit_fdfsolver_alloc(T, arg->nfreq, 5);

  /* Fitting function info */
  arg->funcStruc = g_malloc0(sizeof(gsl_multifit_function_fdf));
  if (doBrokePow) {
    arg->funcStruc->f      = &SpecFitFuncBP;
    arg->funcStruc->df     = &SpecFitJacBP;
    arg->funcStruc->fdf    = &SpecFitFuncJacBP;
  } else {
    arg->funcStruc->f      = &SpecFitFunc;
    arg->funcStruc->df     = &SpecFitJac;
    arg->funcStruc->fdf    = &SpecFitFuncJac;
  }
  arg->funcStruc->n      = arg->nfreq;
  arg->funcStruc->p      = arg->nterm;
  arg->funcStruc->params = arg;

  /* Set up work arrays */
  if (arg->nterm>=2) {
    arg->covar2 = gsl_matrix_alloc(2, 2);
    arg->work2  = gsl_vector_alloc(2);
  }
  if (arg->nterm>=3) {
    arg->covar3 = gsl_matrix_alloc(3, 3);
    arg->work3  = gsl_vector_alloc(3);
  }
  if (arg->nterm>=4) {
    arg->covar4 = gsl_matrix_alloc(4, 4);
    arg->work4  = gsl_vector_alloc(4);
  }
  if (arg->nterm>=5) {
    arg->covar5 = gsl_matrix_alloc(5, 5);
    arg->work5  = gsl_vector_alloc(5);
  }
#endif /* HAVE_GSL */ 

  /* Fit - poly power or broken power */
  if (arg->doBrokePow) {
    /* Broken power law */
    arg->nterm = 2;  /* simple power law */
    arg->funcStruc->f      = &SpecFitFunc;
    arg->funcStruc->df     = &SpecFitJac;
    arg->funcStruc->fdf    = &SpecFitFuncJac;
    NLFit(arg);
    arg->nterm = 3;  /* broken power law */
    arg->funcStruc->f      = &SpecFitFuncBP;
    arg->funcStruc->df     = &SpecFitJacBP;
    arg->funcStruc->fdf    = &SpecFitFuncJacBP;
    NLFitBP(arg);
  } else {
    /* multi term power */
    NLFit(arg);
  }
  
  /* get results parameters + errors */
  out = g_malloc0((2*arg->nterm+1)*sizeof(ofloat));
  for (j=0; j<nterm*2; j++) out[j] = arg->coef[j];
  /* Chi squared */
  out[nterm*2] = arg->ChiSq;
  
  /* Cleanup */
  if (arg->weight)    g_free(arg->weight);
  if (arg->isigma)    g_free(arg->isigma);
  if (arg->obs)       g_free(arg->obs);
  if (arg->nu)        g_free(arg->nu);
  if (arg->logNuOnu0) g_free(arg->logNuOnu0);
  if (arg->coef)      g_free(arg->coef);
#ifdef HAVE_GSL
  if (arg->solver2)   gsl_multifit_fdfsolver_free (arg->solver2);
  if (arg->solver3)   gsl_multifit_fdfsolver_free (arg->solver3);
  if (arg->solver4)   gsl_multifit_fdfsolver_free (arg->solver4);
  if (arg->solver5)   gsl_multifit_fdfsolver_free (arg->solver5);
  if (arg->work2)     gsl_vector_free(arg->work2);
  if (arg->work3)     gsl_vector_free(arg->work3);
  if (arg->work4)     gsl_vector_free(arg->work4);
  if (arg->work5)     gsl_vector_free(arg->work5);
  if (arg->covar2)    gsl_matrix_free(arg->covar2);
  if (arg->covar3)    gsl_matrix_free(arg->covar3);
  if (arg->covar4)    gsl_matrix_free(arg->covar4);
  if (arg->covar5)    gsl_matrix_free(arg->covar5);
  if (arg->funcStruc) g_free(arg->funcStruc);
#endif /* HAVE_GSL */
  g_free(arg);

  return out;
} /* end ObitSpectrumFitSingle */

/**
 * Make single spectrum fitting argument array
 * Without GSL, only the average intensity is determined (no spectrum)
 * \param nfreq    Number of entries in freq, flux, sigma
 * \param nterm    Number of coefficients of powers of log(nu) to fit
 * \param refFreq  Reference frequency (Hz)
 * \param freq     Array of Frequencies (Hz)
 * \param flux     Array of fluxes (Jy) same dim as freq
 * \param sigma    Array of errors (Jy) same dim as freq
 * \param doBrokePow  TRUE if a broken power law fit is desired (3 terms)
 * \param out      Array for output results, should be g_freed when done
 * \param err      Obit error stack object.
 * \return  argument for single spectrum fitting, use ObitSpectrumFitKillArg to dispose.
 */
gpointer ObitSpectrumFitMakeArg (olong nfreq, olong nterm, odouble refFreq, 
				 odouble *freq, gboolean doBrokePow, 
				 ofloat **out, ObitErr *err)
{
  olong i;
  NLFitArg *arg=NULL;
  gchar *routine = "ObitSpectrumFitMakeArg";
  /* GSL implementation */
#ifdef HAVE_GSL
  const gsl_multifit_fdfsolver_type *T=NULL;
#endif /* HAVE_GSL */ 

  if (err->error) return out;

  /* Warn if too many terms asked for */
  if (nterm>5) {
      Obit_log_error(err, OBIT_InfoWarn, 
		    "%s: Asked for %d terms, will limit to 5", routine, nterm);
  }

  /* Warn if incorrect number of terms */
  if (doBrokePow && (nterm!=3)) {
      Obit_log_error(err, OBIT_InfoWarn, 
		    "%s: Must have three terms for broken power law, not %d", 
		     routine, nterm);
  }

  /* Warn if ref. freq <= 0, set to 1.0e9 */
  if (refFreq<=0.0) {
    refFreq = 1.0e9;
    Obit_log_error(err, OBIT_InfoWarn, 
		   "%s: Setting reference Frequency to 1 GHz", routine);
  }

  /* Create function argument */
  arg = g_malloc(sizeof(NLFitArg));
  arg->in             = NULL;     /* Not needed here */
  arg->nfreq          = nfreq;
  arg->nterm          = nterm;
  arg->doError        = TRUE;
  arg->doBrokePow     = doBrokePow;
  arg->doPBCorr       = FALSE;
  arg->maxIter        = 100;
  arg->minDelta       = 1.0e-5;  /* Min step size */
  arg->maxChiSq       = 3.0;     /* max acceptable normalized chi squares */
  arg->refFreq        = refFreq; /* Reference Frequency */
  arg->weight         = g_malloc0(arg->nfreq*sizeof(ofloat));
  arg->isigma         = g_malloc0(arg->nfreq*sizeof(ofloat));
  arg->obs            = g_malloc0(arg->nfreq*sizeof(ofloat));
  arg->nu             = g_malloc0(arg->nfreq*sizeof(ofloat));
  arg->logNuOnu0      = g_malloc0(arg->nfreq*sizeof(ofloat));
  arg->coef           = g_malloc0(2*arg->nterm*sizeof(ofloat));
  for (i=0; i<nfreq; i++) {
    arg->nu[i]        = freq[i];
    arg->logNuOnu0[i] = log(freq[i]/refFreq);
  }

  /* GSL implementation */
#ifdef HAVE_GSL
  arg->solver2 =  arg->solver3 =  arg->solver4 =  arg->solver5 = NULL;
  arg->covar2  =  arg->covar3  =  arg->covar4  =  arg->covar5  = NULL;
  arg->work2   =  arg->work3   =  arg->work4   =  arg->work5   = NULL;
  /* Setup solvers */
  T = gsl_multifit_fdfsolver_lmder;
  if (arg->nterm>=2)
    arg->solver2 = gsl_multifit_fdfsolver_alloc(T, arg->nfreq, 2);
  if (arg->nterm>=3)
    arg->solver3 = gsl_multifit_fdfsolver_alloc(T, arg->nfreq, 3);
  if (arg->nterm>=4)
    arg->solver4 = gsl_multifit_fdfsolver_alloc(T, arg->nfreq, 4);
  if (arg->nterm>=5)
    arg->solver5 = gsl_multifit_fdfsolver_alloc(T, arg->nfreq, 5);

  /* Fitting function info */
  arg->funcStruc = g_malloc0(sizeof(gsl_multifit_function_fdf));
  if (doBrokePow) {
    arg->funcStruc->f      = &SpecFitFuncBP;
    arg->funcStruc->df     = &SpecFitJacBP;
    arg->funcStruc->fdf    = &SpecFitFuncJacBP;
  } else {
    arg->funcStruc->f      = &SpecFitFunc;
    arg->funcStruc->df     = &SpecFitJac;
    arg->funcStruc->fdf    = &SpecFitFuncJac;
  }
  arg->funcStruc->n      = arg->nfreq;
  arg->funcStruc->p      = arg->nterm;
  arg->funcStruc->params = arg;

  /* Set up work arrays */
  if (arg->nterm>=2) {
    arg->covar2 = gsl_matrix_alloc(2, 2);
    arg->work2  = gsl_vector_alloc(2);
  }
  if (arg->nterm>=3) {
    arg->covar3 = gsl_matrix_alloc(3, 3);
    arg->work3  = gsl_vector_alloc(3);
  }
  if (arg->nterm>=4) {
    arg->covar4 = gsl_matrix_alloc(4, 4);
    arg->work4  = gsl_vector_alloc(4);
  }
  if (arg->nterm>=5) {
    arg->covar5 = gsl_matrix_alloc(5, 5);
    arg->work5  = gsl_vector_alloc(5);
  }
#endif /* HAVE_GSL */ 

  /* output array */
  *out = (ofloat*)g_malloc0((2*arg->nterm+1)*sizeof(ofloat));

  return (gpointer)arg;
} /* end ObitSpectrumFitMakeArg */

/**
 * Fit single spectrum to flux measurements using precomputed argument
 * \param aarg      pointer to argument for fitting
 * \param flux      Array of values to be fitted
 * \param sigma     Array of uncertainties of flux
 * \param out       Result array at least 2*nterms+1 in size.
 *                  in order, fitted parameters, error estimates, chi sq of fit.
 */
void ObitSpectrumFitSingleArg (gpointer aarg, ofloat *flux, ofloat *sigma,
			       ofloat *out)
{
  olong i, j;
  NLFitArg *arg=(NLFitArg*)aarg;
  ofloat fblank = ObitMagicF();
  gboolean allBad;

  /* Save flux array, sigma, Check if all data blanked  */
  allBad = TRUE;
   for (i=0; i<arg->nfreq; i++) arg->obs[i] = flux[i];
   for (i=0; i<arg->nfreq; i++) {
     arg->obs[i]    = flux[i];
     if (arg->obs[i]!=fblank) allBad = FALSE;
     arg->isigma[i] = 1.0 / sigma[i];
     arg->weight[i] = arg->isigma[i]*arg->isigma[i];
   }

   /* Return fblanks/zeroes for no data */
   if (allBad) {
     out[0] = fblank;
     for (j=1; j<=arg->nterm*2; j++) out[j] = 0.0;
    return;
  }

  /* Fit - poly power or broken power */
  if (arg->doBrokePow) {
    /* Broken power law */
    arg->nterm = 2;  /* simple power law */
    arg->funcStruc->f      = &SpecFitFunc;
    arg->funcStruc->df     = &SpecFitJac;
    arg->funcStruc->fdf    = &SpecFitFuncJac;
    NLFit(arg);
    arg->nterm = 3;  /* broken power law */
    arg->funcStruc->f      = &SpecFitFuncBP;
    arg->funcStruc->df     = &SpecFitJacBP;
    arg->funcStruc->fdf    = &SpecFitFuncJacBP;
    NLFitBP(arg);
  } else {
    /* multi term power */
    NLFit(arg);
  }
  
  /* get results parameters + errors */
  for (j=0; j<arg->nterm*2; j++) out[j] = arg->coef[j];
  /* Chi squared */
  out[arg->nterm*2] = arg->ChiSq;
  
  return;
} /* end ObitSpectrumFitSingleArg */

/**
 * Delete single fitting argument
 * \param aarg      pointer to argument to kill
 */
void ObitSpectrumFitKillArg (gpointer aarg)
{
  NLFitArg *arg= (NLFitArg*)aarg;

  if (arg==NULL) return;
  if (arg->weight)    g_free(arg->weight);
  if (arg->isigma)    g_free(arg->isigma);
  if (arg->obs)       g_free(arg->obs);
  if (arg->nu)        g_free(arg->nu);
  if (arg->logNuOnu0) g_free(arg->logNuOnu0);
  if (arg->coef)      g_free(arg->coef);
#ifdef HAVE_GSL
  if (arg->solver2)   gsl_multifit_fdfsolver_free (arg->solver2);
  if (arg->solver3)   gsl_multifit_fdfsolver_free (arg->solver3);
  if (arg->solver4)   gsl_multifit_fdfsolver_free (arg->solver4);
  if (arg->solver5)   gsl_multifit_fdfsolver_free (arg->solver5);
  if (arg->work2)     gsl_vector_free(arg->work2);
  if (arg->work3)     gsl_vector_free(arg->work3);
  if (arg->work4)     gsl_vector_free(arg->work4);
  if (arg->work5)     gsl_vector_free(arg->work5);
  if (arg->covar2)    gsl_matrix_free(arg->covar2);
  if (arg->covar3)    gsl_matrix_free(arg->covar3);
  if (arg->covar4)    gsl_matrix_free(arg->covar4);
  if (arg->covar5)    gsl_matrix_free(arg->covar5);
  if (arg->funcStruc) g_free(arg->funcStruc);
#endif /* HAVE_GSL */
  g_free(arg);

} /* end ObitSpectrumFitKillArg */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSpectrumFitClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSpectrumFitClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSpectrumFitClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSpectrumFitClassInfoDefFn (gpointer inClass)
{
  ObitSpectrumFitClassInfo *theClass = (ObitSpectrumFitClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSpectrumFitClassInit;
  theClass->newObit       = (newObitFP)newObitSpectrumFit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSpectrumFitClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSpectrumFitGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitSpectrumFitCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitSpectrumFitClear;
  theClass->ObitInit      = (ObitInitFP)ObitSpectrumFitInit;
  theClass->ObitSpectrumFitCreate = (ObitSpectrumFitCreateFP)ObitSpectrumFitCreate;
  theClass->ObitSpectrumFitCube   = (ObitSpectrumFitCubeFP)ObitSpectrumFitCube;
  theClass->ObitSpectrumFitImArr  = (ObitSpectrumFitImArrFP)ObitSpectrumFitImArr;
  theClass->ObitSpectrumFitEval   = (ObitSpectrumFitEvalFP)ObitSpectrumFitEval;
  theClass->ObitSpectrumFitSingle = (ObitSpectrumFitSingleFP)ObitSpectrumFitSingle;
} /* end ObitSpectrumFitClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitSpectrumFitInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSpectrumFit *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread     = newObitThread();
  in->info       = newObitInfoList(); 
  in->nterm      = 0;
  in->nfreq      = 0;
  in->doBrokePow = FALSE;
  in->maxChi2    = 1.5;
  in->RMS        = NULL;
  in->calFract   = NULL;
  in->outDesc    = NULL;
  in->inFArrays  = NULL;
  in->BeamShapes = NULL;
  in->outFArrays = NULL;
  in->freqs      = NULL;

} /* end ObitSpectrumFitInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitSpectrumFit* cast to an Obit*.
 */
void ObitSpectrumFitClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSpectrumFit *in = inn;
  olong i, nOut;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->RMS) g_free(in->RMS);
  if (in->calFract) g_free(in->calFract);
  in->thread = ObitThreadUnref(in->thread);
  in->info   = ObitInfoListUnref(in->info);
  if (in->outDesc) in->outDesc = ObitImageDescUnref(in->outDesc);
  if (in->inFArrays) {
    for (i=0; i<in->nfreq; i++) in->inFArrays[i] = ObitFArrayUnref(in->inFArrays[i]);
    g_free(in->inFArrays);
  }

  if (in->BeamShapes) {
    for (i=0; i<in->nfreq; i++) in->BeamShapes[i] = ObitBeamShapeUnref(in->BeamShapes[i]);
    g_free(in->BeamShapes);
  }

  /* How many output planes */
  if (in->doError) nOut = 1+in->nterm*2;
  else nOut = in->nterm;

  if (in->outFArrays) {
    for (i=0; i<nOut; i++) in->outFArrays[i] = ObitFArrayUnref(in->outFArrays[i]);
    g_free(in->outFArrays);
  }
  if (in->freqs)     g_free(in->freqs);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSpectrumFitClear */

/**
 * Does spectral fitting, input and output on input object
 * Work divided up amoung 1 or more threads
 * In each pixel determines weighted average and RMS and if the SNR
 * exceeds maxChi2 then the spectrum is fitted. 
 * The pixel value at the reference frequency is either the average pixel 
 * if the SNR is too low or 10^(first term) of fit.
 * The results are stored in outFArrays:
 * \li first entry is the pixel value (Jy/bm) at the reference freq
 * \li entries 1->nterm-1 are the coefficients of log(freq/refFreq)^n
 *     possibly magic value blanked for pixels with no signal
 * \li outFArrays[nterm] = RMS estimate for pixel values (Jy/bm) at refFreq
 * \li entries nterm+1->nterm*2 RMS uncertainties on higher order coefficients
 * \li entry 1+nterm*2 = the Chi squared of the fit
 *
 * \param  in  SpectrumFit to fit
 * \param  err Obit error stack object.
 */
static void Fitter (ObitSpectrumFit* in, ObitErr *err)
{
  olong i, loy, hiy, nyPerThread, nThreads;
  gboolean OK;
  NLFitArg **threadArgs;
  NLFitArg *args=NULL;
  ObitThreadFunc func=(ObitThreadFunc)ThreadNLFit;
  gchar *routine = "Fitter";
  /* GSL implementation */
#ifdef HAVE_GSL
  const gsl_multifit_fdfsolver_type *T=NULL;
#endif /* HAVE_GSL */ 

  /* error checks */
  if (err->error) return;

  /* Warn if too many terms asked for */
  if (in->nterm>5) {
      Obit_log_error(err, OBIT_InfoWarn, 
		    "%s: Asked for %d terms, will limit to 5", routine, in->nterm);
  }

  /* How many threads to use? */
  nThreads = MAX (1, ObitThreadNumProc(in->thread));

  /* Initialize threadArg array  */
  threadArgs = g_malloc0(nThreads*sizeof(NLFitArg*));
  for (i=0; i<nThreads; i++) 
    threadArgs[i] = g_malloc0(sizeof(NLFitArg)); 

  /* Set up thread arguments */
  for (i=0; i<nThreads; i++) {
    args = (NLFitArg*)threadArgs[i];
    args->in          = in;
    args->err         = err;
    args->nfreq       = in->nfreq;
    args->nterm       = in->nterm;
    args->doError     = in->doError;
    args->doBrokePow  = in->doBrokePow;
    args->doPBCorr    = in->doPBCorr;
    args->maxIter     = 100;
    args->minDelta    = 1.0e-2;          /* Min step size */
    args->maxChiSq    = in->maxChi2;     /* max acceptable normalized chi squares */
    args->weight      = g_malloc0(args->nfreq*sizeof(ofloat));
    args->isigma      = g_malloc0(args->nfreq*sizeof(ofloat));
    args->obs         = g_malloc0(args->nfreq*sizeof(ofloat));
    args->nu          = g_malloc0(args->nfreq*sizeof(ofloat));
    args->logNuOnu0   = g_malloc0(args->nfreq*sizeof(ofloat));
    if (args->doError)
      args->coef      = g_malloc0(2*args->nterm*sizeof(ofloat));
    else
      args->coef      = g_malloc0(args->nterm*sizeof(ofloat));
    /* GSL implementation */
#ifdef HAVE_GSL
    args->solver2 =  args->solver3 =  args->solver4 =  args->solver5 = NULL;
    args->covar2  =  args->covar3  =  args->covar4  =  args->covar5  = NULL;
    args->work2   =  args->work3   =  args->work4   =  args->work5   = NULL;
    /* Setup solvers */
    T = gsl_multifit_fdfsolver_lmder;
    if (args->nterm>=2)
      args->solver2 = gsl_multifit_fdfsolver_alloc(T, args->nfreq, 2);
    if (args->nterm>=3)
      args->solver3 = gsl_multifit_fdfsolver_alloc(T, args->nfreq, 3);
    if (args->nterm>=4)
      args->solver4 = gsl_multifit_fdfsolver_alloc(T, args->nfreq, 4);
    if (args->nterm>=5)
      args->solver5 = gsl_multifit_fdfsolver_alloc(T, args->nfreq, 5);
    
    /* Fitting function info */
    args->funcStruc = g_malloc0(sizeof(gsl_multifit_function_fdf));
    if (in->doBrokePow) {
      args->funcStruc->f      = &SpecFitFuncBP;
      args->funcStruc->df     = &SpecFitJacBP;
      args->funcStruc->fdf    = &SpecFitFuncJacBP;
    } else {
      args->funcStruc->f      = &SpecFitFunc;
      args->funcStruc->df     = &SpecFitJac;
      args->funcStruc->fdf    = &SpecFitFuncJac;
    }
    args->funcStruc->n      = args->nfreq;
    args->funcStruc->p      = args->nterm;
    args->funcStruc->params = args;
    
    /* Set up work arrays */
    if (args->nterm>=2) {
      args->covar2 = gsl_matrix_alloc(2, 2);
      args->work2  = gsl_vector_alloc(2);
    }
    if (args->nterm>=3) {
      args->covar3 = gsl_matrix_alloc(3, 3);
      args->work3  = gsl_vector_alloc(3);
    }
    if (args->nterm>=4) {
      args->covar4 = gsl_matrix_alloc(4, 4);
      args->work4  = gsl_vector_alloc(4);
    }
    if (args->nterm>=5) {
      args->covar5 = gsl_matrix_alloc(5, 5);
      args->work5  = gsl_vector_alloc(5);
    }
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
    args = (NLFitArg*)threadArgs[i];
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
      args = (NLFitArg*)threadArgs[i];
      if (args->weight)    g_free(args->weight);
      if (args->isigma)    g_free(args->isigma);
      if (args->obs)       g_free(args->obs);
      if (args->nu)        g_free(args->nu);
      if (args->logNuOnu0) g_free(args->logNuOnu0);
      if (args->coef)      g_free(args->coef);
#ifdef HAVE_GSL
      if (args->solver2)   gsl_multifit_fdfsolver_free (args->solver2);
      if (args->solver3)   gsl_multifit_fdfsolver_free (args->solver3);
      if (args->solver4)   gsl_multifit_fdfsolver_free (args->solver3);
      if (args->solver4)   gsl_multifit_fdfsolver_free (args->solver3);
      if (args->work2)     gsl_vector_free(args->work2);
      if (args->work3)     gsl_vector_free(args->work3);
      if (args->work4)     gsl_vector_free(args->work4);
      if (args->work5)     gsl_vector_free(args->work5);
      if (args->covar2)    gsl_matrix_free(args->covar2);
      if (args->covar3)    gsl_matrix_free(args->covar3);
      if (args->covar4)    gsl_matrix_free(args->covar4);
      if (args->covar5)    gsl_matrix_free(args->covar5);
      if (args->funcStruc) g_free(args->funcStruc);
#endif /* HAVE_GSL */
     g_free(threadArgs[i]);
    }
    g_free(threadArgs);
  }

} /* end Fitter */

/**
 * Write contents on in to outImage
 * \param in       Spectral fitting object
 * \param outImage Image cube with fitted spectra.
 *                 Should be defined but not created.
 *                 Planes 1->nterm are coefficients per pixel
 *                 Planes nterm+1->2*nterm are uncertainties in coefficients
 * \param err      Obit error stack object.
 */
static void WriteOutput (ObitSpectrumFit* in, ObitImage *outImage, 
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

  /* If doBrokePow set "BROKENPO" to TRUE on output descriptor */
  if (in->doBrokePow) {
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (outImage->myDesc->info, "BROKENPO", OBIT_bool, dim, &in->doBrokePow);
    ObitInfoListAlwaysPut (((ObitImageDesc*)outImage->myIO->myDesc)->info, "BROKENPO", 
			   OBIT_bool, dim, &in->doBrokePow);
  }

  /* Close output */
  retCode = ObitImageClose (outImage, err);
  outImage->extBuffer = FALSE;   /* May need I/O buffer later */
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, outImage->name);

} /* end WriteOutput */

/**
 * Thread function to fit a portion of the image set
 * \param arg      NLFitArg structure
 */
static gpointer ThreadNLFit (gpointer arg)
{
  NLFitArg *larg      = (NLFitArg*)arg;
  ObitSpectrumFit* in = (ObitSpectrumFit*)larg->in;
  olong lo            = larg->first-1;  /* First in y range */
  olong hi            = larg->last;     /* Highest in y range */
  gboolean doError    = larg->doError;  /* Error analysis? */
  gboolean doPBCorr   = larg->doPBCorr; /* Primary beam correction? */
  ObitErr *err        = larg->err;

  olong ix, iy, indx, i, nOut;
  odouble Angle=0.0, pos[2];
  ofloat pbfact, pixel[2];
  ofloat fblank = ObitMagicF();
  ObitBeamShapeClassInfo *BSClass;
  gchar *routine = "ThreadNLFit";

  /* error checks */
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Set up frequency info */
  larg->refFreq = in->refFreq;
  for (i=0; i<in->nfreq; i++) {
    larg->logNuOnu0[i] = log(in->freqs[i]/in->refFreq);
    larg->nu[i] = in->freqs[i];
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

      /* Primary beam correction? */
      if (doPBCorr) {
	/* Distance from Center */
	pixel[0] = (ofloat)(ix+1.0); pixel[1] = (ofloat)(iy+1.0);
	ObitImageDescGetPos(in->outDesc, pixel, pos, err);
	if (err->error) {
	  ObitThreadLock(in->thread);  /* Lock against other threads */
	  Obit_log_error(err, OBIT_Error,"%s Error with position %s",
			 routine, in->name);
	  ObitThreadUnlock(in->thread); 
	  return NULL;
	}
	BSClass = (ObitBeamShapeClassInfo*)(in->BeamShapes[0]->ClassInfo);
	Angle = BSClass->ObitBeamShapeAngle(in->BeamShapes[0], pos[0], pos[1], 0.0);
      }

      /* Collect values;  */
      for (i=0; i<in->nfreq; i++) {
	larg->obs[i] = in->inFArrays[i]->array[indx];
	  if (larg->obs[i]!=fblank) {
	    /* Statistical weight */
	    larg->isigma[i] = 1.0 / (in->RMS[i]*in->RMS[i] + 
				     in->calFract[i]*in->calFract[i]*
				     larg->obs[i]*larg->obs[i]);
	    larg->weight[i] = larg->isigma[i];
	    larg->isigma[i] = sqrt(larg->isigma[i]);
	    /* Primary beam correction */
	    if (doPBCorr) {
	      BSClass = (ObitBeamShapeClassInfo*)(in->BeamShapes[i]->ClassInfo);
	      pbfact  = BSClass->ObitBeamShapeGainSym(in->BeamShapes[i], Angle);
	      larg->obs[i] /= pbfact;
	      larg->weight[i]  *= pbfact*pbfact;
	    }
	    /* End if datum valid */
	  } else { /* invalid pixel */
	    larg->weight[i] = 0.0;
	    larg->isigma[i] = 0.0;
	  }
      } /* end loop over frequencies */
      
      /* initialize with values from last time */
      if (doError) {
	for (i=0; i<nOut-1; i++) 
	  larg->coef[i] = in->outFArrays[i]->array[indx];
      } else { /* only values */
 	for (i=0; i<in->nterm; i++) 
	  larg->coef[i] = in->outFArrays[i]->array[indx];
      }
      /* Fit - poly power or broken power */
      if (larg->doBrokePow) {
	/* Broken power law */
	in->nterm = 2;  /* simple power law */
	larg->funcStruc->f      = &SpecFitFunc;
	larg->funcStruc->df     = &SpecFitJac;
	larg->funcStruc->fdf    = &SpecFitFuncJac;
	NLFit(larg);
	in->nterm = 3;  /* broken power law */
	larg->funcStruc->f      = &SpecFitFuncBP;
	larg->funcStruc->df     = &SpecFitJacBP;
	larg->funcStruc->fdf    = &SpecFitFuncJacBP;
	NLFitBP(larg);
      } else {
	/* multi term power */
 	NLFit(larg);
      }
      
      /* Save to output */
      if (doError) {
	for (i=0; i<nOut-1; i++) 
	  in->outFArrays[i]->array[indx]  = larg->coef[i];
	in->outFArrays[nOut-1]->array[indx] = larg->ChiSq;
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
} /* end ThreadNLFit */

/**
 * Do non linear fit to a spectrum
 * Only fits for up to 5 terms
 * \param arg      NLFitArg structure
 *                 fitted parameters returned in arg->in->coef
 */
static void NLFit (NLFitArg *arg)
{
  olong iter=0, i, nterm=arg->nterm, nvalid, best;
  ofloat avg, delta, chi2Test, sigma, fblank = ObitMagicF();
  ofloat meanSNR, SNRperTerm=5.0;
  odouble sum, sumwt, sum2;
  gboolean isDone;
  int status;
#ifdef HAVE_GSL
  gsl_multifit_fdfsolver *solver=NULL;
  gsl_matrix *covar=NULL;
  gsl_vector *work=NULL;
#endif /* HAVE_GSL */ 
 
  /* Initialize output */
  if (arg->doError) 
    for (i=0; i<2*arg->nterm; i++) arg->coef[i] = 0.0;
  else
    for (i=0; i<arg->nterm; i++) arg->coef[i] = 0.0;
  
  /* determine weighted average, count valid data */
  sum = sumwt = 0.0;
  nvalid = 0;
  for (i=0; i<arg->nfreq; i++) {
    if ((arg->obs[i]!=fblank) && (arg->weight[i]>0.0)) {
      sum   += arg->weight[i] * arg->obs[i];
      sumwt += arg->weight[i];
      nvalid++;
    }
  }
  if (nvalid<=0) return;  /* any good data? */
  avg = sum/sumwt;

  /* Estimate of noise */
  sigma = 1.0 / sqrt(sumwt);

  /* Initial fit */
  arg->coef[0] = avg;
  best = 1;  /* best fit number of terms */
  
  /* determine chi squared, mean SNR */
  sum = sum2 = sumwt = 0.0;
  for (i=0; i<arg->nfreq; i++) {
    if ((arg->obs[i]!=fblank) && (arg->weight[i]>0.0)) {
      delta = (arg->obs[i]-avg);
      sumwt += arg->weight[i] * delta*delta;
      sum   += delta*delta;
      sum2  += fabs(arg->obs[i])*sqrt(arg->weight[i]);
    }
  }

  /* normalized Chi squares */
  if (nvalid>1) {
    arg->ChiSq = sumwt/(nvalid-1);
    meanSNR    = sum2/nvalid; /* mean SNR */
  } else {
    arg->ChiSq = -1.0;
    meanSNR    = 0.0; /* Only one value - don't need mean SNR */
  }

  /* Errors wanted? */
  if (arg->doError) arg->coef[arg->nterm] = sum/nvalid;

  /* Is this good enough? */
  isDone = (arg->ChiSq<0.0) || (arg->ChiSq<=arg->maxChiSq);
  if (meanSNR>(SNRperTerm*2.0)) isDone = FALSE;  /* Always try for high SNR */
  if (isDone) goto done;

  /* DEBUG 
  if (avg>2.0) {
    fprintf (stderr, "Found one %f\n", avg);
  } */

  /* Higher order terms do nonlinear least-squares fit */
  nterm = 2;
  while ((nterm<=MIN(5,arg->nterm)) && (nterm<=nvalid)) {
    arg->fitTerm = nterm;  /* How many actually being fitted */
#ifdef HAVE_GSL

    /* Set solver and covar */
    if (nterm==2) {
      solver = arg->solver2;
      covar  = arg->covar2;
      work   = arg->work2;
    } else if (nterm==3) {
      solver = arg->solver3;
      covar  = arg->covar3;
      work   = arg->work3;
    } else if (nterm==4) {
      solver = arg->solver4;
      covar  = arg->covar4;
      work   = arg->work4;
    } else if (nterm==5) {
      solver = arg->solver5;
      covar  = arg->covar5;
      work   = arg->work5;
    } 
    
    /* set initial guess - start with lower order fit */
    for (i=0; i<nterm; i++) gsl_vector_set(work, i, (double)arg->coef[i]);
    arg->funcStruc->n      = arg->nfreq;
    arg->funcStruc->p      = nterm;
    arg->funcStruc->params = arg;
    gsl_multifit_fdfsolver_set (solver, arg->funcStruc, work);
    iter = 0;
    
    /* iteration loop */
    do {
      iter++;
      status = gsl_multifit_fdfsolver_iterate(solver);
      /*if (status) break;???*/

      status = gsl_multifit_test_delta (solver->dx, solver->x, 
					(double)arg->minDelta, 
					(double)arg->minDelta);
    } while ((status==GSL_CONTINUE) && (iter<arg->maxIter));

    /* If it didn't work - bail */
    if ((status!=GSL_SUCCESS) && (status!=GSL_CONTINUE)) {
      fprintf (stderr, "Failed, status = %s\n", gsl_strerror(status));
      return;
    }
    
    /* normalized Chi squares */
    if (nvalid>nterm) {
      sumwt = (ofloat)gsl_blas_dnrm2(solver->f);
      chi2Test = (sumwt*sumwt)/(nvalid-nterm);
    } else chi2Test = -1.0;

    /* Did it improve over lower order? */
    if (chi2Test<arg->ChiSq) {
      best = nterm;
      arg->ChiSq = chi2Test;

      /* Get fitted values */
      for (i=0; i<nterm; i++) arg->coef[i] = (ofloat)gsl_vector_get(solver->x, i);
      
      /* Errors wanted? */
      if (arg->doError) {
	gsl_multifit_covar (solver->J, 0.0, covar);
	for (i=0; i<nterm; i++) {
	  arg->coef[arg->nterm+i] = sqrt(gsl_matrix_get(covar, i, i));
	  /* Clip to sanity range */
	  arg->coef[arg->nterm+i] = MAX (-1.0e5, MIN (1.0e5, arg->coef[arg->nterm+i]));
	}
      } /* end of get errors */
      
      /* Sanity check on spectral parameters [-3,3]
	 for (i=1; i<nterm; i++) arg->coef[i] = MAX(-3.0, MIN(3.0,arg->coef[i])); */
      
    }  /* End if better than lower order */

    /* Is this good enough? */
    isDone = (arg->ChiSq<0.0) || (arg->ChiSq<=arg->maxChiSq) || (chi2Test>arg->ChiSq);
    if ((meanSNR>(SNRperTerm*nterm)) && (nterm<arg->nterm)) 
      isDone = FALSE;  /* Always try for high SNR */
    if (isDone) goto done;
    /*    if ((arg->ChiSq>0.0) && (arg->ChiSq<=arg->maxChiSq)) goto done;*/

    nterm++;  /* next term */
  } /* end loop over adding terms */
#endif /* HAVE_GSL */ 
 done:
  /* sanity check, if flux < sigma, don't include higher order terms */
  if (fabs(arg->coef[0])<sigma)  
    for (i=1; i<arg->nterm; i++) arg->coef[i] = 0.0; 

  /*  Gonzo higher order fit */
  if ((fabs(arg->coef[1])>3.0) || ((nterm>2) && (fabs(arg->coef[2])>2.0))) {
    arg->coef[0] = avg;
    for (i=1; i<arg->nterm; i++) arg->coef[i] = 0.0; 
  }

} /* end NLFit */

/**
 * Do non linear fit to a broken power law spectrum
 * Only fits 3 terms
 * Broken power law version, above nv_b the spectrum steepens by half
 * Parameters:
 *   S_0   = Flux density at frequency nu_0
 *   alpha = spectral index
 *   nu_b  = Break frequency
 * When called the S_O and the alpha terms should have been estimated
 * (call to NLFit with 2 terms); the initial break frequency is set to 
 * the reference frequency. 
 * \param arg      NLFitArg structure
 *                 fitted parameters returned in arg->in->coef
 */
static void NLFitBP (NLFitArg *arg)
{
  olong iter=0, i, nterm, nvalid;
  ofloat avg, chi2Test, sigma, fblank = ObitMagicF();
  odouble sum, sumwt;
  int status;
#ifdef HAVE_GSL
  gsl_multifit_fdfsolver *solver;
  gsl_matrix *covar;
  gsl_vector *work;
#endif /* HAVE_GSL */ 
 
  /* determine weighted average, count valid data */
  sum = sumwt = 0.0;
  nvalid = 0;
  for (i=0; i<arg->nfreq; i++) {
    if ((arg->obs[i]!=fblank) && (arg->weight[i]>0.0)) {
      sum   += arg->weight[i] * arg->obs[i];
      sumwt += arg->weight[i];
      nvalid++;
    }
  }
  if (nvalid<=(arg->nterm)) return;  /* enough good data? */
  avg = sum/sumwt;
      
  /* Estimate of noise */
  sigma = 1.0 / sqrt(sumwt);

  /* Do nonlinear least-squares fit */
  nterm = 3;
  arg->fitTerm = nterm;  /* How many actually being fitted */
#ifdef HAVE_GSL

  /* Set solver and covar */
  solver = arg->solver3;
  covar  = arg->covar3;
  work   = arg->work3;
         
  /* set initial guess - start with lower order fit */
  /* Initial break frequency is reference freq */
  arg->coef[2] = arg->refFreq;
  for (i=0; i<nterm; i++) gsl_vector_set(work, i, (double)arg->coef[i]);
  arg->funcStruc->n      = arg->nfreq;
  arg->funcStruc->p      = nterm;
  arg->funcStruc->params = arg;
  gsl_multifit_fdfsolver_set (solver, arg->funcStruc, work);
  iter = 0;
      
  /* iteration loop */
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);
    /*if (status) break;???*/
    
    status = gsl_multifit_test_delta (solver->dx, solver->x, 
				      (double)arg->minDelta, 
				      (double)arg->minDelta);
  } while ((status==GSL_CONTINUE) && (iter<arg->maxIter));
  
  /* If it didn't work - bail */
  if ((status!=GSL_SUCCESS) && (status!=GSL_CONTINUE)) {
    fprintf (stderr, "Failed, status = %s\n", gsl_strerror(status));
    return;
  }
  
  /* normalized Chi squares */
  if (nvalid>nterm) {
    sumwt = (ofloat)gsl_blas_dnrm2(solver->f);
    chi2Test = (sumwt*sumwt)/(nvalid-nterm);
  } else chi2Test = -1.0;

  /* Is this better than a simple power law? */
  /*if (chi2Test<arg->ChiSq) {*/
  arg->ChiSq = chi2Test;
  
  /* Get fitted values */
  for (i=0; i<nterm; i++) arg->coef[i] = (ofloat)gsl_vector_get(solver->x, i);

  /* DEBUG 
  if (arg->coef[0]>50.0*sigma) {
    fprintf (stdout, "DEBUG %f %f %f \n",arg->coef[0], arg->coef[1], arg->coef[2]);
  }*/
  
  /* Errors wanted? */
  if (arg->doError) {
    gsl_multifit_covar (solver->J, 0.0, covar);
    for (i=0; i<nterm; i++) {
      arg->coef[arg->nterm+i] = sqrt(gsl_matrix_get(covar, i, i));
      /* Clip to sanity range
      arg->coef[arg->nterm+i] = MAX (-1.0e5, MIN (1.0e5, arg->coef[arg->nterm+i])); */
    }
  } /* end of get errors */
  
  
    /* Sanity check on spectral parameters [-3,3]
       for (i=1; i<nterm; i++) arg->coef[i] = MAX(-3.0, MIN(3.0,arg->coef[i])); */
    /* End if better than simple power law */
    /*} else {*/
    /* Worse, use simple power law Make break very high */
    /*    arg->coef[2] = 1.0e20;*/
    /*}*/
  
#endif /* HAVE_GSL */ 
  /* sanity check, if flux < sigma, don't include higher order terms */
  if (fabs(arg->coef[0])<sigma)  
    for (i=1; i<arg->nterm; i++) arg->coef[i] = 0.0; 

  /*  Gonzo higher order fit 
  if ((fabs(arg->coef[1])>3.0) || (arg->coef[2]>arg->nu[arg->nfreq-1])) {
    arg->coef[0] = avg;
    for (i=1; i<arg->nterm; i++) arg->coef[i] = 0.0; 
  }*/

} /* end NLFitBP */

#ifdef HAVE_GSL
/**
 * Function evaluator for spectral fitting solver
 * Evaluates (model-observed) / sigma
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of spectral_terms
 * \param param   Function parameter structure (NLFitArg)
 * \param f       Vector of (model-obs)/sigma for data points
 * \return completion code GSL_SUCCESS=OK
 */
static int SpecFitFunc (const gsl_vector *x, void *params, 
			gsl_vector *f)
{
  NLFitArg *args = (NLFitArg*)params;
  ofloat fblank = ObitMagicF();
  double flux, spec[10], func;
  odouble model, arg, argx, isigma, spect;
  size_t i, j;

  /* get model parameters */
  flux = gsl_vector_get(x, 0);
  for (j=1; j<args->fitTerm; j++) spec[j-1] = gsl_vector_get(x, j);

  /* Loop over spectral points */
  for (i=0; i<args->nfreq; i++) {
    if ((args->obs[i]!=fblank) && (args->isigma[i]>0.0)) {
      isigma = args->isigma[i];
      arg = 0.0;
      argx = args->logNuOnu0[i];
      for  (j=1; j<args->fitTerm; j++) {
	arg += spec[j-1]*argx;
	argx *= args->logNuOnu0[i];
      }
      spect = exp(arg);
      model = flux * spect;
      func =  (model - args->obs[i]) * isigma;
      gsl_vector_set(f, i, func);  /* Save function residual */
    } else {  /* Invalid data */
      func = 0.0;
      gsl_vector_set(f, i, func);     /* Save function residual */
    }
 } /* End loop over spectral points */

  return GSL_SUCCESS;
} /*  end SpecFitFunc */

/**
 * Jacobian evaluator for spectral fitting solver
 * Evaluates partial derivatives of model wrt each parameter
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of spectral_terms
 * \param param   Function parameter structure (NLFitArg)
 * \param J       Jacobian matrix J[data_point, parameter]
 * \return completion code GSL_SUCCESS=OK
 */
static int SpecFitJac (const gsl_vector *x, void *params, 
		       gsl_matrix *J)
{
  NLFitArg *args = (NLFitArg*)params;
  ofloat fblank = ObitMagicF();
  double flux, spec[10], jac;
  odouble model, arg, argx, isigma, spect;
  size_t i, j;

  /* get model parameters */
  flux = gsl_vector_get(x, 0);
  for (j=1; j<args->fitTerm; j++) spec[j-1] = gsl_vector_get(x, j);

  /* Loop over spectral points */
  for (i=0; i<args->nfreq; i++) {
    if ((args->obs[i]!=fblank) && (args->isigma[i]>0.0)) {
      isigma = args->isigma[i];
      arg = 0.0;
      argx = args->logNuOnu0[i];
      for  (j=1; j<args->fitTerm; j++) {
	arg += spec[j-1]*argx;
	argx *= args->logNuOnu0[i];
      }
      spect = exp(arg);
      model = flux * spect;
      /* Jacobian terms - first flux */
      jac = spect * isigma;
      gsl_matrix_set(J, i, 0, jac);   /* Save flux Jacobian */
      argx = args->logNuOnu0[i];
      for  (j=1; j<args->fitTerm; j++) {
	jac = model * isigma * argx;
	gsl_matrix_set(J, i, j, jac);  /* Save Jacobian */
 	argx *= args->logNuOnu0[i];
      }
    } else {  /* Invalid data */
      jac = 0.0;
      gsl_matrix_set(J, i, 0, jac);  /* Save flux Jacobian */
      for  (j=1; j<args->fitTerm; j++) {
	gsl_matrix_set(J, i, j, jac); /* Save Jacobian */
      }
    }
  } /* End loop over spectral points */

  return GSL_SUCCESS;
} /*  end SpecFitJac */

/**
 * Function and Jacobian evaluator for spectral fitting solver
 * Function = (model-observed) / sigma
 * Jacobian =  partial derivatives of model wrt each parameter
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of spectral_terms
 * \param param   Function parameter structure (NLFitArg)
 * \param f       Vector of (model-obs)/sigma for data points
 * \param J       Jacobian matrix J[data_point, parameter]
 * \return completion code GSL_SUCCESS=OK
 */
static int SpecFitFuncJac (const gsl_vector *x, void *params, 
			   gsl_vector *f, gsl_matrix *J)
{
  NLFitArg *args = (NLFitArg*)params;
  ofloat fblank = ObitMagicF();
  double flux, spec[10], func, jac;
  odouble model, arg, argx, isigma, spect;
  size_t i, j;

  /* get model parameters */
  flux = gsl_vector_get(x, 0);
  for (j=1; j<args->fitTerm; j++) spec[j-1] = gsl_vector_get(x, j);

  /* Loop over spectral points */
  for (i=0; i<args->nfreq; i++) {
    if ((args->obs[i]!=fblank) && (args->isigma[i]>0.0)) {
      isigma = args->isigma[i];
      arg = 0.0;
      argx = args->logNuOnu0[i];
      for  (j=1; j<args->fitTerm; j++) {
	arg += spec[j-1]*argx;
	argx *= args->logNuOnu0[i];
      }
      spect = exp(arg);
      model = flux * spect;
      func =  (model - args->obs[i]) * isigma;
      gsl_vector_set(f, i, func);  /* Save function residual */
      /* Jacobian terms - first flux */
      jac = spect * isigma;
      gsl_matrix_set(J, i, 0, jac);   /* Save flux Jacobian */
      argx = args->logNuOnu0[i];
      for  (j=1; j<args->fitTerm; j++) {
	jac = model * isigma * argx;
	gsl_matrix_set(J, i, j, jac);  /* Save Jacobian */
 	argx *= args->logNuOnu0[i];
     }
    } else {  /* Invalid data */
      func = 0.0;
      gsl_vector_set(f, i, func);     /* Save function residual */
      jac = 0.0;
      gsl_matrix_set(J, i, 0, jac);  /* Save flux Jacobian */
      for  (j=1; j<args->fitTerm; j++) {
	gsl_matrix_set(J, i, j, jac); /* Save Jacobian */
      }
    }
  } /* End loop over spectral points */

  return GSL_SUCCESS;
} /*  end SpecFitFuncJac */
/**
 * Function evaluator for spectral fitting solver
 * Broken power law version, above nv_b the spectrum steepens by half
 * Parameters:
 *   S_0   = Flux density at frequency nu_0
 *   alpha = spectral index
 *   nu_b  = Break frequency
 * Evaluates (model-observed) / sigma
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of spectral_terms
 * \param param   Function parameter structure (NLFitArg)
 * \param f       Vector of (model-obs)/sigma for data points
 * \return completion code GSL_SUCCESS=OK
 */
static int SpecFitFuncBP (const gsl_vector *x, void *params, 
			  gsl_vector *f)
{
  NLFitArg *args = (NLFitArg*)params;
  ofloat fblank = ObitMagicF();
  double flux, spec[10], func;
  odouble model, fluxb, arg, isigma, spect;
  size_t i, j;

  /* get model parameters */
  flux = gsl_vector_get(x, 0);
  /* Spectral terms  - alpha, nu_b */
  for (j=1; j<args->fitTerm; j++) spec[j-1] = gsl_vector_get(x, j);

  /* Loop over spectral points */
  for (i=0; i<args->nfreq; i++) {
    if ((args->obs[i]!=fblank) && (args->isigma[i]>0.0)) {
      /* Above or below the break? */
      if (args->nu[i]<spec[1]) { /* Below - simple power law */
	arg = spec[0] *  args->logNuOnu0[i];
	spect = exp(arg);
	model = flux * spect;
      } else { /* above */
      /* Get the flux density at the break */
	arg = spec[0] * log (spec[1] / args->refFreq);
	spect = exp(arg);
	fluxb = flux * spect;
	/* Flux density at frequency */
	arg = (spec[0]-0.5) * log (args->nu[i] / spec[1]);
	spect = exp(arg);
	model = fluxb * spect;
      } /* End of above or below the break */
      isigma = args->isigma[i];
      func =  (model - args->obs[i]) * isigma;
      gsl_vector_set(f, i, func);  /* Save function residual */
    } else {  /* Invalid data */
      func = 0.0;
      gsl_vector_set(f, i, func);     /* Save function residual */
    }
 } /* End loop over spectral points */

  return GSL_SUCCESS;
} /*  end SpecFitFuncBP */

/**
 * Jacobian evaluator for spectral fitting solver
 * Broken power law version, above nv_b the spectrum steepens by half
 * Parameters:
 *   S_0   = Flux density at frequency nu_0
 *   alpha = spectral index
 *   nu_b  = Break frequency
 * Evaluates partial derivatives of model wrt each parameter
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of spectral_terms
 * \param param   Function parameter structure (NLFitArg)
 * \param J       Jacobian matrix J[data_point, parameter]
 * \return completion code GSL_SUCCESS=OK
 */
static int SpecFitJacBP (const gsl_vector *x, void *params, 
			 gsl_matrix *J)
{
  NLFitArg *args = (NLFitArg*)params;
  ofloat fblank = ObitMagicF();
  double flux, spec[10], jac;
  odouble model, fluxb, arg, isigma, spect, spect2;
  odouble logNubONu0, logNuONub;
  size_t i, j;

  /* get model parameters */
  flux = gsl_vector_get(x, 0);
  /* Spectral terms - alpha, nu_b */
  for (j=1; j<args->fitTerm; j++) spec[j-1] = gsl_vector_get(x, j);

  /* Loop over spectral points */
  for (i=0; i<args->nfreq; i++) {
    if ((args->obs[i]!=fblank) && (args->isigma[i]>0.0)) {
      isigma = args->isigma[i];
      /* Above or below the break? */
      if (args->nu[i]<spec[1]) { /* Below - simple power law */
	arg = spec[0] *  args->logNuOnu0[i];
	spect = exp(arg);
	model = flux * spect;
	/* Jacobian terms - first flux */
	jac = spect * isigma;
	gsl_matrix_set(J, i, 0, jac);   /* Save flux Jacobian */
	/* Spectral index */
	jac = model * isigma * args->logNuOnu0[i];
	gsl_matrix_set(J, i, 1, jac);  /* Save spectral index Jacobian */
	/* Break frequency - no dependency */
	jac = 0.0;
	gsl_matrix_set(J, i, 2, jac);  /* Save spectral index Jacobian */
      } else { /* above */
	/* Get the flux density at the break */
	logNubONu0 = log (spec[1] / args->refFreq);
	arg = spec[0] * logNubONu0;
	spect = exp(arg);
	fluxb = flux * spect;
	/* Flux density at frequency */
	logNuONub = log (args->nu[i] / spec[1]);
	arg = (spec[0]-0.5) * logNuONub;
	spect2 = exp(arg);
	model  = fluxb * spect2;
	/* Jacobian terms - first flux */
	jac = spect * spect2 * isigma;
	gsl_matrix_set(J, i, 0, jac);   /* Save flux Jacobian */
	/* Spectral index */
	jac = (fluxb*logNubONu0 + model*logNuONub) * isigma;
	gsl_matrix_set(J, i, 1, jac);  /* Save spectral index Jacobian */
	/* Break frequency */
	jac = isigma * (fluxb*spec[0] - model*(spec[0]-0.5)) / spec[1];
	gsl_matrix_set(J, i, 2, jac);  /* Save break frequency Jacobian */
      } /* End of above or below the break */
    } else {  /* Invalid data */
      jac = 0.0;
      gsl_matrix_set(J, i, 0, jac);  /* Save flux Jacobian */
      for  (j=1; j<args->fitTerm; j++) {
	gsl_matrix_set(J, i, j, jac); /* Save Jacobian */
      }
    }
  } /* End loop over spectral points */

  return GSL_SUCCESS;
} /*  end SpecFitJacBP */

/**
 * Function and Jacobian evaluator for spectral fitting solver
 * Broken power law version, above nv_b the spectrum steepens by half
 * Parameters:
 *   S_0   = Flux density at frequency nu_0
 *   alpha = spectral index
 *   nu_b  = Break frequency
 * Function = (model-observed) / sigma
 * Jacobian =  partial derivatives of model wrt each parameter
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of spectral_terms
 * \param param   Function parameter structure (NLFitArg)
 * \param f       Vector of (model-obs)/sigma for data points
 * \param J       Jacobian matrix J[data_point, parameter]
 * \return completion code GSL_SUCCESS=OK
 */
static int SpecFitFuncJacBP (const gsl_vector *x, void *params, 
			   gsl_vector *f, gsl_matrix *J)
{
  NLFitArg *args = (NLFitArg*)params;
  ofloat fblank = ObitMagicF();
  double flux, spec[10], func, jac;
  odouble model, arg, fluxb, isigma, spect, spect2;
  odouble logNubONu0, logNuONub;
  size_t i, j;

  /* get model parameters */
  flux = gsl_vector_get(x, 0);
  /* Spectral terms - alpha, nu_b */
  for (j=1; j<args->fitTerm; j++) spec[j-1] = gsl_vector_get(x, j);

  /* Loop over spectral points */
  for (i=0; i<args->nfreq; i++) {
    if ((args->obs[i]!=fblank) && (args->isigma[i]>0.0)) {
      isigma = args->isigma[i];
      /* Above or below the break? */
      if (args->nu[i]<spec[1]) { /* Below - simple power law */
	arg = spec[0] *  args->logNuOnu0[i];
	spect = exp(arg);
	model = flux * spect;
	func =  (model - args->obs[i]) * isigma;
	gsl_vector_set(f, i, func);  /* Save function residual */
	/* Jacobian terms - first flux */
	jac = spect * isigma;
	gsl_matrix_set(J, i, 0, jac);   /* Save flux Jacobian */
	/* Spectral index */
	jac = model * isigma * args->logNuOnu0[i];
	gsl_matrix_set(J, i, 1, jac);  /* Save spectral index Jacobian */
	/* Break frequency - no dependency */
	jac = 0.0;
	gsl_matrix_set(J, i, 2, jac);  /* Save spectral index Jacobian */
      } else { /* above */
	/* Get the flux density at the break */
	logNubONu0 = log (spec[1] / args->refFreq);
	arg = spec[0] * logNubONu0;
	spect = exp(arg);
	fluxb = flux * spect;
	/* Flux density at frequency */
	logNuONub = log (args->nu[i] / spec[1]);
	arg = (spec[0]-0.5) * logNuONub;
	spect2 = exp(arg);
	model  = fluxb * spect2;
	func =  (model - args->obs[i]) * isigma;
	gsl_vector_set(f, i, func);  /* Save function residual */
	/* Jacobian terms - first flux */
	jac = spect * spect2 * isigma;
	gsl_matrix_set(J, i, 0, jac);   /* Save flux Jacobian */
	/* Spectral index */
	jac = (fluxb*logNubONu0 + model*logNuONub) * isigma;
	gsl_matrix_set(J, i, 1, jac);  /* Save spectral index Jacobian */
	/* Break frequency */
	jac = isigma * (fluxb*spec[0] - model*(spec[0]-0.5)) / spec[1];
	gsl_matrix_set(J, i, 2, jac);  /* Save break frequency Jacobian */
      } /* End of above or below the break */
    } else {  /* Invalid data */
      func = 0.0;
      gsl_vector_set(f, i, func);     /* Save function residual */
      jac = 0.0;
      gsl_matrix_set(J, i, 0, jac);  /* Save flux Jacobian */
      for  (j=1; j<args->fitTerm; j++) {
	gsl_matrix_set(J, i, j, jac); /* Save Jacobian */
      }
    }
  } /* End loop over spectral points */

  return GSL_SUCCESS;
} /*  end SpecFitFuncJacBP */
#endif /* HAVE_GSL */ 
