/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2025                                               */
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

#include "ObitFaraSyn.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFaraSyn.c
 * ObitFaraSyn class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitFaraSyn";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitFaraSynClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitFaraSynClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitFaraSynInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitFaraSynClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitFaraSynClassInfoDefFn (gpointer inClass);

/** Direct RMAna search */
static void DirectAnaSearch(ObitFaraSyn *in, ObitCArray *Accum, ObitFArray *workFA, 
			    ObitCArray *workCA, ObitErr *err);

/** Private: Write RMAna output */
static void WriteRMAna (ObitFaraSyn *in, ObitErr *err);

/*---------------Private structures----------------*/
/* FT threaded function argument */
typedef struct {
  /** ObitFaraSyn object */
  ObitFaraSyn *in;
  /* Obit error stack object */
  ObitErr      *err;
  /** First (1-rel) element to process this thread */
  olong        first;
  /** Highest (1-rel) element to process this thread  */
  olong        last;
  /** increment to process this thread  */
  olong        incr;
  /** current index in FArrays */
  olong        indx;
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
  ofloat *wrk1, *wrk2, *wrk3, *wrk4, *wrk5, *wrk6, *wrk7, *wrk8, *wrk9;
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
} FaraFitArg;

/** Private: Do actual fitting */
static void LSQFitter (ObitFaraSyn* in, ObitErr *err);

/** Private: Threaded fitting */
static gpointer ThreadFaraFit (gpointer arg);

/** Private: Actual least squares fitting */
static void FaraFit (FaraFitArg *arg);

#ifdef HAVE_GSL
/** Private: Solver function evaluation */
static int FaraFitFunc (const gsl_vector *x, void *params, 
			gsl_vector *f);

/** Private: Solver Jacobian evaluation */
static int FaraFitJac (const gsl_vector *x, void *params, 
		       gsl_matrix *J);

/** Private: Solver function + Jacobian evaluation */
static int FaraFitFuncJac (const gsl_vector *x, void *params, 
			   gsl_vector *f, gsl_matrix *J);

#endif /* HAVE_GSL */ 

/** Private: Determine amplitude for least squares fit */
static void FaraFitAmp (FaraFitArg *arg);


/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitFaraSyn* newObitFaraSyn (gchar* name)
{
  ObitFaraSyn* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFaraSynClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitFaraSyn));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitFaraSynInit((gpointer)out);

 return out;
} /* end newObitFaraSyn */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitFaraSynGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFaraSynClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitFaraSynGetClass */

/**
 * Make a deep copy of an ObitFaraSyn.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitFaraSyn* ObitFaraSynCopy  (ObitFaraSyn *in, ObitFaraSyn *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, oldnOut, oldnLamb2;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitFaraSyn(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  oldnLamb2     = out->nLamb2;
  oldnOut       = out->nOut;
  out->nLamb2   = in->nLamb2;
  out->refLamb2 = in->refLamb2;
  out->nOut     = in->nOut;
  out->qMedian  = in->qMedian;   out->qSigma   = in->qSigma; 
  out->uMedian  = in->uMedian;   out->uSigma   = in->uSigma; 
  out->minRMSyn = in->minRMSyn;  out->maxRMSyn = in->maxRMSyn;
  out->delRMSyn = in->delRMSyn;
  out->doWrite  = in->doWrite;
  out->minQUSNR = in->minQUSNR;
  out->minFrac  = in->minFrac;
  out->doError  = in->doError;
  out->doRMSyn  = in->doRMSyn;
  out->maxChi2  = in->maxChi2;
  /* Arrays of odouble or ofloats */
  if ((in->nLamb2>0) && (in->lamb2)) {
    if (out->lamb2) g_free(out->lamb2);
    out->lamb2 = g_malloc0(in->nLamb2*sizeof(odouble));
    for (i=0; i<in->nLamb2; i++) out->lamb2[i] = in->lamb2[i];
  }
  if ((in->nLamb2>0) && (in->QRMS)) {
    if (out->QRMS) g_free(out->QRMS);
    out->QRMS = g_malloc0(in->nLamb2*sizeof(ofloat));
    for (i=0; i<in->nLamb2; i++) out->QRMS[i] = in->QRMS[i];
  }
  if ((in->nLamb2>0) && (in->URMS)) {
    if (out->URMS) g_free(out->URMS);
    out->URMS = g_malloc0(in->nLamb2*sizeof(ofloat));
    for (i=0; i<in->nLamb2; i++) out->URMS[i] = in->URMS[i];
  }
  if ((in->nLamb2>0) && (in->plnWt)) {
    if (out->plnWt) g_free(out->plnWt);
    out->plnWt = g_malloc0(in->nLamb2*sizeof(ofloat));
    for (i=0; i<in->nLamb2; i++) out->plnWt[i] = in->plnWt[i];
  }
  /* Arrays of FArray or CArray - delete any old, reference any new */
  if (out->inQFArrays) {
    for (i=0; i<oldnLamb2; i++) 
      if (out->inQFArrays[i]) out->inQFArrays[i] = ObitFArrayUnref(out->inQFArrays[i]);
    g_free(out->inQFArrays);
  }
  if ((in->nLamb2>0) && (in->inQFArrays)) {
    out->inQFArrays = g_malloc0(in->nLamb2*sizeof(ObitFArray));
    for (i=0; i<in->nLamb2; i++) out->inQFArrays[i] = ObitFArrayRef(in->inQFArrays[i]);
  }
  if (out->inUFArrays) {
    for (i=0; i<oldnLamb2; i++)
      if (out->inUFArrays[i]) out->inUFArrays[i] = ObitFArrayUnref(out->inUFArrays[i]);
    g_free(out->inUFArrays);
  }
  if ((in->nLamb2>0) && (in->inUFArrays)) {
    out->inUFArrays = g_malloc0(in->nLamb2*sizeof(ObitFArray));
    for (i=0; i<in->nLamb2; i++) out->inUFArrays[i] = ObitFArrayRef(in->inUFArrays[i]);
  }
  if (out->outFArrays) {
    for (i=0; i<oldnOut; i++)
      if (out->outFArrays[i]) out->outFArrays[i] = ObitFArrayUnref(out->outFArrays[i]);
    g_free(out->outFArrays);
  }
  if ((in->nOut>0) && (in->outFArrays)) {
    out->outFArrays = g_malloc0(in->nOut*sizeof(ObitFArray));
    for (i=0; i<in->nOut; i++) out->outFArrays[i] = ObitFArrayRef(in->outFArrays[i]);
  }
  if (out->RMPlanes) {
    for (i=0; i<oldnOut; i++) 
      if (out->RMPlanes[i]) out->RMPlanes[i] = ObitCArrayUnref(out->RMPlanes[i]);
    g_free(out->RMPlanes);
  }
  if ((in->nOut>0) && (in->RMPlanes)) {
    out->RMPlanes = g_malloc0(in->nOut*sizeof(ObitCArray));
    for (i=0; i<in->nOut; i++) out->RMPlanes[i] = ObitCArrayRef(in->RMPlanes[i]);
  }
 /* Images */
 if (out->outAImage) out->outAImage = ObitImageUnref(out->outAImage);
 if (in->outAImage)  out->outAImage = ObitImageRef(in->outAImage);
 if (out->outPImage) out->outPImage = ObitImageUnref(out->outPImage);
 if (in->outPImage)  out->outPImage = ObitImageRef(in->outPImage);
 if (out->outFitRM)  out->outFitRM  = ObitImageUnref(out->outFitRM);
 if (in->outFitRM)   out->outFitRM  = ObitImageRef(in->outFitRM);

  return out;
} /* end ObitFaraSynCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an FaraSyn similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitFaraSynClone  (ObitFaraSyn *in, ObitFaraSyn *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i, oldnOut, oldnLamb2;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* Copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  oldnLamb2     = out->nLamb2;
  oldnOut       = out->nOut;
  out->nLamb2   = in->nLamb2;
  out->refLamb2 = in->refLamb2;
  out->nOut     = in->nOut;
  out->qMedian  = in->qMedian;   out->qSigma   = in->qSigma; 
  out->uMedian  = in->uMedian;   out->uSigma   = in->uSigma; 
  out->minRMSyn = in->minRMSyn;  out->maxRMSyn = in->maxRMSyn;
  out->delRMSyn = in->delRMSyn;
  out->doWrite  = in->doWrite;
  out->minQUSNR = in->minQUSNR;
  out->minFrac  = in->minFrac;
  out->doError  = in->doError;
  out->doRMSyn  = in->doRMSyn;
  out->maxChi2  = in->maxChi2;
  /* Arrays of odouble or ofloats */
  if ((in->nLamb2>0) && (in->lamb2)) {
    if (out->lamb2) g_free(out->lamb2);
    out->lamb2 = g_malloc0(in->nLamb2*sizeof(odouble));
    for (i=0; i<in->nLamb2; i++) out->lamb2[i] = in->lamb2[i];
  }
  if ((in->nLamb2>0) && (in->QRMS)) {
    if (out->QRMS) g_free(out->QRMS);
    out->QRMS = g_malloc0(in->nLamb2*sizeof(ofloat));
    for (i=0; i<in->nLamb2; i++) out->QRMS[i] = in->QRMS[i];
  }
  if ((in->nLamb2>0) && (in->URMS)) {
    if (out->URMS) g_free(out->URMS);
    out->URMS = g_malloc0(in->nLamb2*sizeof(ofloat));
    for (i=0; i<in->nLamb2; i++) out->URMS[i] = in->URMS[i];
  }
  if ((in->nLamb2>0) && (in->plnWt)) {
    if (out->plnWt) g_free(out->plnWt);
    out->plnWt = g_malloc0(in->nLamb2*sizeof(ofloat));
    for (i=0; i<in->nLamb2; i++) out->plnWt[i] = in->plnWt[i];
  }
  /* Arrays of FArray or CArray - delete any old, reference any new */
  if (out->inQFArrays) {
    for (i=0; i<oldnLamb2; i++) 
      if (out->inQFArrays[i]) out->inQFArrays[i] = ObitFArrayUnref(out->inQFArrays[i]);
    g_free(out->inQFArrays);
  }
  if ((in->nLamb2>0) && (in->inQFArrays)) {
    out->inQFArrays = g_malloc0(in->nLamb2*sizeof(ObitFArray));
    for (i=0; i<in->nLamb2; i++) out->inQFArrays[i] = ObitFArrayRef(in->inQFArrays[i]);
  }
  if (out->inUFArrays) {
    for (i=0; i<oldnLamb2; i++)
      if (out->inUFArrays[i]) out->inUFArrays[i] = ObitFArrayUnref(out->inUFArrays[i]);
    g_free(out->inUFArrays);
  }
  if ((in->nLamb2>0) && (in->inUFArrays)) {
    out->inUFArrays = g_malloc0(in->nLamb2*sizeof(ObitFArray));
    for (i=0; i<in->nLamb2; i++) out->inUFArrays[i] = ObitFArrayRef(in->inUFArrays[i]);
  }
  if (out->outFArrays) {
    for (i=0; i<oldnOut; i++)
      if (out->outFArrays[i]) out->outFArrays[i] = ObitFArrayUnref(out->outFArrays[i]);
    g_free(out->outFArrays);
  }
  if ((in->nOut>0) && (in->outFArrays)) {
    out->outFArrays = g_malloc0(in->nOut*sizeof(ObitFArray));
    for (i=0; i<in->nOut; i++) out->outFArrays[i] = ObitFArrayRef(in->outFArrays[i]);
  }
  if (out->RMPlanes) {
    for (i=0; i<oldnOut; i++) 
      if (out->RMPlanes[i]) out->RMPlanes[i] = ObitCArrayUnref(out->RMPlanes[i]);
    g_free(out->RMPlanes);
  }
  if ((in->nOut>0) && (in->RMPlanes)) {
    out->RMPlanes = g_malloc0(in->nOut*sizeof(ObitCArray));
    for (i=0; i<in->nOut; i++) out->RMPlanes[i] = ObitCArrayRef(in->RMPlanes[i]);
  }
 /* Images */
 if (out->outAImage) out->outAImage = ObitImageUnref(out->outAImage);
 if (in->outAImage)  out->outAImage = ObitImageRef(in->outAImage);
 if (out->outPImage) out->outPImage = ObitImageUnref(out->outPImage);
 if (in->outPImage)  out->outPImage = ObitImageRef(in->outPImage);
 if (out->outFitRM)  out->outFitRM  = ObitImageUnref(out->outFitRM);
 if (in->outFitRM)   out->outFitRM  = ObitImageRef(in->outFitRM);
} /* end ObitFaraSynClone */

/**
 * Creates an ObitFaraSyn 
 * Note: Not all parameters are used in all modes, 
 * give dummy values (NULL or 0) to those not needed.
 * \param   name       A name for the object.
 * \param   nLamb2     Number of Lambda^2 planes                                                   
 * \param   lamb2      Wavelength^2 (m^2) of planes                                                
 * \param   refLamb2   Reference Wavelength^2                                                      
 * \param   inQFArrays Array of nlamb2 FArrays with Q data                                         
 * \param   inUFArrays Array of nlamb2 FArrays with U data                                         
 * \param   plnWt      Weights for each Q/U plane                                                  
 * \param   QRMS       RMS of each Q plane                                                         
 * \param   URMS       RMS of each U plane                                                         
 * \param   nOut       Number of Faraday depth planes                                              
 * \param   minRMSyn   Min. Faraday depth (rad/m^2)                                                
 * \param   maxRMSyn   Max. Faraday depth (rad/m^2) - Faraday Analysis   
 * \param   delRMSyn   Increment in Faraday depth (rad/m^2)                                        
 * \param   doWrite    If true write Amplitudes to outAImage                                       
 * \param   RMPlanes   Array of nOut CArrays for Faraday depth                                     
 * \param   outAImage  Output amp cube if doWrite                                                  
 * \param   outPImage  Output phase cube if non NULL                                               
 * \param   outFitRM   Faraday analysis cube if non NULL                                               
 * \param   minFrac    Min. fraction of planes included - Faraday Analysis 
 * \param   minQUSNR   Min. SNR for Q and U pixels - Faraday Analysis 
 * \param   doError    If true do error analysis - Faraday Analysis 
 * \param   doRMSyn    If true do max RM synthesis - Faraday Analysis
 * \param   maxChi2    Max. Chi^2 for search - Faraday Analysis 
 * \return the new object.
 */
ObitFaraSyn* ObitFaraSynCreate 
  (gchar* name,	
   olong nLamb2, odouble *lamb2, odouble refLamb2,
   ObitFArray **inQFArrays, ObitFArray **inUFArrays,
   ofloat *plnWt, ofloat *QRMS, ofloat *URMS,
   olong nOut, ofloat minRMSyn, ofloat maxRMSyn, ofloat delRMSyn,
   ObitCArray **RMPlanes, 
   gboolean doWrite, ObitImage* outAImage, ObitImage* outPImage,
   ObitImage *outFitRM, ofloat minQUSNR, ofloat minFrac, 
   gboolean doError, gboolean doRMSyn, ofloat maxChi2)
{
  ObitFaraSyn* out;
  olong i;

  /* Create basic structure */
  out = newObitFaraSyn (name);
  /* Set members given */
  out->nLamb2   = nLamb2;
  out->refLamb2 = refLamb2;
  out->nOut     = nOut;
  out->minRMSyn = minRMSyn;
  out->maxRMSyn = maxRMSyn;
  out->delRMSyn = delRMSyn ;
  out->doWrite  = doWrite;
  out->minQUSNR = minQUSNR;
  out->minFrac  = minFrac;
  out->doError  = doError;
  out->doRMSyn  = doRMSyn;
  out->maxChi2  = maxChi2;
 /* Arrays of odouble or ofloats */
  if ((nLamb2>0) && (lamb2)) {
    out->lamb2 = g_malloc0(nLamb2*sizeof(odouble));
    for (i=0; i<nLamb2; i++) out->lamb2[i] = lamb2[i];
  }
  if ((nLamb2>0) && (QRMS)) {
    out->QRMS = g_malloc0(nLamb2*sizeof(ofloat));
    for (i=0; i<nLamb2; i++) out->QRMS[i] = QRMS[i];
  }
  if ((nLamb2>0) && (URMS)) {
    out->URMS = g_malloc0(nLamb2*sizeof(ofloat));
    for (i=0; i<nLamb2; i++) out->URMS[i] = URMS[i];
  }
   if ((nLamb2>0) && (plnWt)) {
    out->plnWt = g_malloc0(nLamb2*sizeof(ofloat));
    for (i=0; i<nLamb2; i++) out->plnWt[i] = plnWt[i];
  }
  /* Arrays of FArray  */
  if ((nLamb2>0) && (inQFArrays)) {
    out->inQFArrays = g_malloc0(nLamb2*sizeof(ObitFArray));
    for (i=0; i<nLamb2; i++) out->inQFArrays[i] = ObitFArrayRef(inQFArrays[i]);
  }
  if ((nLamb2>0) && (inUFArrays)) {
    out->inUFArrays = g_malloc0(nLamb2*sizeof(ObitFArray));
    for (i=0; i<nLamb2; i++) out->inUFArrays[i] = ObitFArrayRef(inUFArrays[i]);
  }
 if ((nOut>0) && (RMPlanes)) {
    out->RMPlanes = g_malloc0(nOut*sizeof(ObitCArray));
    for (i=0; i<nOut; i++) out->RMPlanes[i] = ObitCArrayRef(RMPlanes[i]);
  }
 /* Images */
 if (outAImage)  out->outAImage = ObitImageRef(outAImage);
 if (outPImage)  out->outPImage = ObitImageRef(outPImage);
 if (outFitRM)   out->outFitRM  = ObitImageRef(outFitRM);
 return out;
} /* end ObitFaraSynCreate */

/**
 * Set Q & U data 
 * \param   in         ObitFaraSyn object to update
 * \param   nLamb2     Number of Lambda^2 planes                                                   
 * \param   lamb2      Wavelength^2 (m^2) of planes                                                
 * \param   refLamb2   Reference Wavelength^2                                                      
 * \param   inQFArrays Array of nlamb2 FArrays with Q data                                         
 * \param   inUFArrays Array of nlamb2 FArrays with U data                                         
 * \param   plnWt      Weights for each Q/U plane                                                  
 * \param   QRMS       RMS of each Q plane                                                         
 * \param   URMS       RMS of each U plane     
 */                                                    
void ObitFaraSynSetQU (ObitFaraSyn* in, 
		       olong nLamb2, odouble *lamb2, odouble refLamb2,
		       ObitFArray **inQFArrays, ObitFArray **inUFArrays,
		       ofloat *plnWt, ofloat *QRMS, ofloat *URMS)
{
  olong i,oldnLamb2;

  /* Set members given */
  oldnLamb2    = in->nLamb2;
  in->nLamb2   = nLamb2;
  in->refLamb2 = refLamb2;
  /* Arrays of odouble or ofloats */
  if ((nLamb2>0) && (lamb2)) {
    if (in->lamb2) g_free(in->lamb2);
    in->lamb2 = g_malloc0(nLamb2*sizeof(odouble));
    for (i=0; i<nLamb2; i++) in->lamb2[i] = lamb2[i];
  }
  if ((nLamb2>0) && (QRMS)) {
    if (in->QRMS) g_free(in->QRMS);
    in->QRMS = g_malloc0(nLamb2*sizeof(ofloat));
    for (i=0; i<nLamb2; i++) in->QRMS[i] = QRMS[i];
  }
  if ((nLamb2>0) && (URMS)) {
    if (in->URMS) g_free(in->URMS);
    in->URMS = g_malloc0(nLamb2*sizeof(ofloat));
    for (i=0; i<nLamb2; i++) in->URMS[i] = URMS[i];
  }
   if ((nLamb2>0) && (plnWt)) {
    if (in->plnWt) g_free(in->plnWt);
    in->plnWt = g_malloc0(nLamb2*sizeof(ofloat));
    for (i=0; i<nLamb2; i++) in->plnWt[i] = plnWt[i];
  }
  /* Arrays of FArray or CArray - delete any old, reference any new */
  if (in->inQFArrays) {
    for (i=0; i<oldnLamb2; i++) 
      if (in->inQFArrays[i]) in->inQFArrays[i] = ObitFArrayUnref(in->inQFArrays[i]);
    g_free(in->inQFArrays);
  }
  if ((nLamb2>0) && (inQFArrays)) {
    in->inQFArrays = g_malloc0(nLamb2*sizeof(ObitFArray));
    for (i=0; i<nLamb2; i++) in->inQFArrays[i] = ObitFArrayRef(inQFArrays[i]);
  }
  if (in->inUFArrays) {
    for (i=0; i<oldnLamb2; i++)
      if (in->inUFArrays[i]) in->inUFArrays[i] = ObitFArrayUnref(in->inUFArrays[i]);
    g_free(in->inUFArrays);
  }
  if ((nLamb2>0) && (inUFArrays)) {
    in->inUFArrays = g_malloc0(nLamb2*sizeof(ObitFArray));
    for (i=0; i<nLamb2; i++) in->inUFArrays[i] = ObitFArrayRef(inUFArrays[i]);
  }
} /* end ObitFaraSynSetQU */

/**
 * Set RM Synthesis request/output 
 * Q/U arrays should have been set in ObitFaraSynCreate or ObitFaraSynSetQU
 * \param   in         ObitFaraSyn object to update
 * \param   nOut       Number of Faraday depth planes                                              
 * \param   minRMSyn   Min. Faraday depth (rad/m^2)                                                
 * \param   delRMSyn   Increment in Faraday depth (rad/m^2)                                        
 * \param   doWrite    If true write Amplitudes to outAImage                                       
 * \param   RMPlanes   Array of nOut CArrays for Faraday depth                                     
 * \param   outAImage  Output amp cube if doWrite                                                  
 * \param   outPImage  Output phase cube if non NULL          
 */                                     
 void ObitFaraSynSetRMSyn (ObitFaraSyn* in,
			   olong nOut, ofloat minRMSyn, ofloat delRMSyn,
			   ObitCArray **RMPlanes, 
			   gboolean doWrite, ObitImage* outAImage, 
			   ObitImage* outPImage)
{
  olong i, oldnOut;

  /* Set members given */
  oldnOut      = in->nOut;
  in->nOut     = nOut;
  in->minRMSyn = minRMSyn;
  in->delRMSyn = delRMSyn ;
  in->doWrite  = doWrite;
  /* Arrays of C/FArray  - delete any old, reference any new */
  if (in->RMPlanes) {
    for (i=0; i<oldnOut; i++) 
      if (in->RMPlanes[i]) in->RMPlanes[i] = ObitCArrayUnref(in->RMPlanes[i]);
    g_free(in->RMPlanes);
  }
  if ((nOut>0) && (RMPlanes)) {
    in->RMPlanes = g_malloc0(nOut*sizeof(ObitCArray));
    for (i=0; i<nOut; i++) in->RMPlanes[i] = ObitCArrayRef(RMPlanes[i]);
  }
  /* Images */
  if (in->outAImage) in->outAImage = ObitImageUnref(in->outAImage);
  if (outAImage)     in->outAImage = ObitImageRef(outAImage);
  if (in->outPImage) in->outPImage = ObitImageUnref(in->outPImage);
  if (outPImage)     in->outPImage = ObitImageRef(outPImage);
 } /* end ObitFaraSynSetRMSyn */

/**
 * Set RM Analysis request/output 
 * Q/U arrays should have been set in ObitFaraSynCreate or ObitFaraSynSetQU
 * \param   in         ObitFaraSyn object to update
 * \param   minRMSyn   Min. Faraday depth (rad/m^2)                                                
 * \param   maxRMSyn   Max. Faraday depth (rad/m^2)                                                
 * \param   delRMSyn   Increment in Faraday depth (rad/m^2)                                        
 * \param   outFitRM   Faraday analysis cube if non NULL                                               
 * \param   minQUSNR   Min. SNR for Q and U pixels - Faraday Analysis 
 * \param   minFrac    Min. fraction of planes included - Faraday Analysis
 * \param   doError    If true do error analysis - Faraday Analysis 
 * \param   doRMSyn    If true do max RM synthesis - Faraday Analysis
 * \param   maxChi2    Max. Chi^2 for search - Faraday Analysis 
 */                                     
 void ObitFaraSynSetRMAna (ObitFaraSyn* in,
			   ofloat minRMSyn, float maxRMSyn, ofloat delRMSyn,
			   ObitImage* outFitRM, ofloat minQUSNR, ofloat minFrac, 
			   gboolean doError, gboolean doRMSyn, ofloat maxChi2)
{
  /* Set members given */
  in->minRMSyn = minRMSyn;
  in->minRMSyn = maxRMSyn;
  in->delRMSyn = delRMSyn;
  in->minQUSNR = minQUSNR;
  in->minFrac  = minFrac;
  in->doError  = doError;
  in->doRMSyn  = doRMSyn;
  in->maxChi2  = maxChi2;
  /* Image */
  if (in->outFitRM) in->outAImage = ObitImageUnref(in->outFitRM);
  if (outFitRM)     in->outFitRM  = ObitImageRef(outFitRM);
 } /* end ObitFaraSynSetRMAna */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitFaraSynClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitFaraSynClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitFaraSynClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitFaraSynClassInfoDefFn (gpointer inClass)
{
  ObitFaraSynClassInfo *theClass = (ObitFaraSynClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitFaraSynClassInit;
  theClass->newObit       = (newObitFP)newObitFaraSyn;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitFaraSynClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitFaraSynGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitFaraSynCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitFaraSynClear;
  theClass->ObitInit      = (ObitInitFP)ObitFaraSynInit;
  theClass->ObitFaraSynCreate = (ObitFaraSynCreateFP)ObitFaraSynCreate;
  theClass->ObitFaraSynRMSyn = (ObitFaraSynRMSynFP)ObitFaraSynRMSyn;
  theClass->ObitFaraSynRMAna = (ObitFaraSynRMAnaFP)ObitFaraSynRMAna;

} /* end ObitFaraSynClassDefFn */

/**
 *  Core of Faraday synthesis - Faraday depth cube
 * \param   in         Input object using members:
 *        Input:
 * \li   nLamb2     Number of Lambda^2 planes
 * \li   lamb2      Wavelength^2 (m^2) of planes
 * \li   refLamb2   Reference Wavelength^2
 * \li   inQFArrays Array of nlamb2 FArrays with Q data
 * \li   inUFArrays Array of nlamb2 FArrays with U data
 * \li   plnWt      Weights for each Q/U plane
 * \li   QRMS       RMS of each Q plane
 * \li   URMS       RMS of each U plane
 * \li   nOut       Number of Faraday depth planes
 * \li   minRMSyn   Min. Faraday depth (rad/m^2)
 * \li   delRMSyn   Increment in Faraday depth (rad/m^2)
 * \li   doWrite    If true write Amplitudes to outAImage
 *                     Also phases to outPImage if non NULL
 *        Output:
 * \li   RMPlanes   Array of nOut CArrays for Faraday depth
 * \li   outAImage  Output amp cube if doWrite 
 * \li   outPImage  Output phase cube if non NULL
 * \param   err        Obit Error stack 
 */
void ObitFaraSynRMSyn (ObitFaraSyn *in, ObitErr *err)
{
  /* Local variables */
  olong nlamb2 = in->nLamb2;
  odouble *lamb2 = in->lamb2;
  odouble refLamb2 = in->refLamb2;
  olong nOut = in->nOut;
  ObitFArray **inQFArrays = in->inQFArrays, **inUFArrays = in->inUFArrays;
  ObitCArray **RMPlanes = in->RMPlanes;
  ObitImage *outAImage = in->outAImage, *outPImage = in->outPImage;
  ofloat *plnWt = in->plnWt, *QRMS = in->QRMS, *URMS = in->URMS;
  ofloat minRMSyn = in->minRMSyn;
  ofloat delRMSyn = in->delRMSyn;
  gboolean doWrite = in->doWrite;

  olong i, iRM, nThreads=1, naxis[2], nx, ny, plane[5]={1,1,1,1,1};
  ofloat RM, cmplx[2], sumWt, norm, *work=NULL, maxQ, maxU;
  gpointer threadArgs=NULL;
  ObitIOCode retCode;
  ObitFArray *workAmp=NULL;
  gchar *routine = "ObitFaraSynRMSyn";

  /* Statistics for filtering only median within 2 sigma */
  work    = g_malloc0(nlamb2*sizeof(ofloat));  /* work array */
  for (i=0; i<nlamb2; i++) work[i] = QRMS[i];
  in->qMedian = medianValue(work, 1, nlamb2);
  in->qSigma  = MedianSigma(nlamb2, work, in->qMedian);
  for (i=0; i<nlamb2; i++) work[i] = URMS[i];
  in->uMedian = medianValue(work, 1, nlamb2);
  in->uSigma  = MedianSigma(nlamb2, work, in->uMedian);
  if (work) g_free(work);
  /* Three sigma Upper limits on RMS */
  maxQ = in->qMedian+3*in->qSigma;
  maxU = in->uMedian+3*in->uSigma;

  /* Work array if needed */
  if (doWrite) {
    /* Image size */
    nx = outAImage->myDesc->inaxes[0];
    ny = outAImage->myDesc->inaxes[1];
    naxis[0] = (olong)nx;  naxis[1] = (olong)ny; 
    workAmp  = ObitFArrayCreate (NULL, 2, naxis);
  }
 
  /* Loop over RM planes */
  for (iRM=1; iRM<=nOut; iRM++) {
    /* Setup FT threading */
    if (threadArgs==NULL)
      nThreads = ObitCArrayMakeCAFuncArgs (in->thread, NULL, NULL, 
					   RMPlanes[0], inQFArrays[0], inUFArrays[0],
					   NULL, NULL, 2, 0, 0, 0, 0, 0, 0, &threadArgs);
    RM = minRMSyn + (iRM-1)*delRMSyn;
    /* Zero acccumulator */
    cmplx[0] = 0.0; cmplx[1] = 0.0;
    ObitCArrayFill(RMPlanes[iRM-1], cmplx);
    /* loop over input planes */
    sumWt = 0.0;
    for (i=0; i<nlamb2; i++) {
      /* Want this one ? */
      if ((plnWt[i]<=0.0) || !inQFArrays[i] || !inUFArrays[i] || 
	  (QRMS[i]>maxQ) || (URMS[i]>maxU)) continue;
      sincosf((ofloat)(-2.0*RM*(lamb2[i]-refLamb2)), &cmplx[1], &cmplx[0]); /* sin/cos factors */
      /* Rotate, accumulate this plane */
      ObitCArraySMulAccumTh (inQFArrays[i], inUFArrays[i], cmplx, RMPlanes[iRM-1],
			     nThreads, threadArgs);
      sumWt += plnWt[i];
     } /* Loop over lamb2 planes */
    /* Normalize - no threading conflict with ObitFArraySMul */
    norm = 1.0 / sumWt;
    ObitCArraySMul(RMPlanes[iRM-1], norm);

    /* Write if requested - these functions not threaded */
    if (doWrite) {
      ObitCArrayKillCAFuncArgs (nThreads, threadArgs);  /* These routines threaded */
      threadArgs = NULL;
      /* Get ampl - write to output if needed */
      ObitCArrayAmp(RMPlanes[iRM-1], workAmp);
      plane[0] = iRM;  /* Select correct plane */
      retCode = ObitImagePutPlane (outAImage, workAmp->array, plane, err);
      /* Also phase? */
      if (outPImage) {
	ObitCArrayPhase(RMPlanes[iRM-1], workAmp);
	retCode = ObitImagePutPlane (outPImage, workAmp->array, plane, err);
      } /* write phase */
      if ((retCode!=OBIT_IO_OK) || (err->error)) {
	Obit_traceback_msg (err, routine, outAImage->name);
	goto cleanup;
      } /* end write */
    }
   } /* end loop over RM */
  /* Shutdown FT threading if not already done */
  if (threadArgs) {
    ObitCArrayKillCAFuncArgs (nThreads, threadArgs);  /* cleanup */
    threadArgs = NULL;
  }
  cleanup:
  if (workAmp) ObitFArrayUnref(workAmp);
} /* end ObitFaraSynRMSyn */

/**
 *  Core of Faraday analysis -  Peak Faraday depth cube 
 * \param   in         Input object using members:
 *        Input:
 * \li   nLamb2     Number of Lambda^2 planes
 * \li   lamb2      Wavelength^2 (m^2) of planes
 * \li   refLamb2   Reference Wavelength^2
 * \li   inQFArrays Array of nlamb2 FArrays with Q data
 * \li   inUFArrays Array of nlamb2 FArrays with U data
 * \li   plnWt      Weights for each Q/U plane
 * \li   QRMS       RMS of each Q plane
 * \li   URMS       RMS of each U plane
 * \li   nOut       Number of planes in output.
 * \li   minRMSyn   Min. Faraday depth (rad/m^2)
 * \li   maxRMSyn   Max. Faraday depth (rad/m^2)
 * \li   delRMSyn   Increment in Faraday depth (rad/m^2)
 * \li   minQUSNR   Min SNR for Q and U pixels
 * \li   minFrac    Min. fraction of planes included
 * \li   doError    If true do error analysis
 * \li   doRMSyn    If true do max RM synthesis
 * \li   maxChi2    Max. Chi^2 for search
 *        Output:
 * \li   outFitRM   Output Faraday analysis cube
 * \param   err     Obit Error stack 
 */
void ObitFaraSynRMAna (ObitFaraSyn *in, ObitErr *err)
{
  /* Local variables */
  olong nlamb2 = in->nLamb2, nOut = in->nOut;
  ObitImage *outFitRM = in->outFitRM;
  ofloat *work, *QRMS = in->QRMS, *URMS = in->URMS;
  gboolean doRMSyn = in->doRMSyn;

  olong i, naxis[2];
  ObitCArray *Accum=NULL, *workCA=NULL;
  ObitFArray *workFA=NULL;
  /*ObitFArray **outFArrays;*/
  gchar *routine = "ObitFaraSynRMAna";

  /* Statistics for filtering  */
  work    = g_malloc0(nlamb2*sizeof(ofloat));  /* work array */
  for (i=0; i<nlamb2; i++) work[i] = QRMS[i];
  in->qMedian = medianValue(work, 1, nlamb2);
  in->qSigma  = MedianSigma(nlamb2, work, in->qMedian);
  for (i=0; i<nlamb2; i++) work[i] = URMS[i];
  in->uMedian = medianValue(work, 1, nlamb2);
  in->uSigma  = MedianSigma(nlamb2, work, in->uMedian);
  if (work) g_free(work);

  /* Define term arrays */
  naxis[0] = outFitRM->myDesc->inaxes[0];
  naxis[1] = outFitRM->myDesc->inaxes[1];
  in->outFArrays = g_malloc0((nOut)*sizeof(ObitFArray*));
  for (i=0; i<nOut; i++) in->outFArrays[i] = ObitFArrayCreate (NULL,2,naxis);

  /* Accumulation & work arrays */
  Accum   = ObitCArrayCreate (NULL, 2, naxis); /* Accumulator */
  workFA  = ObitFArrayCreate (NULL, 2, naxis); /* Peak amp^2 */
  workCA  = ObitCArrayCreate (NULL, 2, naxis); /* r/i of peak */

  /* Do initial direct search */
  DirectAnaSearch(in, Accum, workFA, workCA, err);
  if (err->error) {
    Obit_traceback_msg (err, routine, outFitRM->name);
    goto cleanup;
  }

  /* Least Squares fitting if needed */
  if (!doRMSyn) LSQFitter (in, err);
  if (err->error) {
    Obit_traceback_msg (err, routine, outFitRM->name);
    goto cleanup;
  }

  /* Save output */
  WriteRMAna (in, err);
  if (err->error) Obit_traceback_msg (err, routine, outFitRM->name);
  goto cleanup;  /* done*/

  cleanup:
  Accum   = ObitCArrayUnref(Accum);
  workFA  = ObitFArrayUnref(workFA);
  workCA  = ObitFArrayUnref(workCA);
} /* end ObitFaraSynRMAna */

/*---------------Private functions--------------------------*/
/**
 *  Direct search in Faraday depth for peak unwrapped PPol
 * \param   in         Input object using members:
 *        Input:
 * \li   nLamb2     Number of Lambda^2 planes
 * \li   lamb2      Wavelength^2 (m^2) of planes
 * \li   refLamb2   Reference Wavelength^2
 * \li   inQFArrays Array of nlamb2 FArrays with Q data
 * \li   inUFArrays Array of nlamb2 FArrays with U data
 * \li   plnWt      Weights for each Q/U plane
 * \li   QRMS       RMS of each Q plane
 * \li   URMS       RMS of each U plane
 * \li   nOut       Number of planes in output.
 * \li   minRMSyn   Min. Faraday depth (rad/m^2)
 * \li   maxRMSyn   Max. Faraday depth (rad/m^2)
 * \li   delRMSyn   Increment in Faraday depth (rad/m^2)
 * \li   minQUSNR   Min SNR for Q and U pixels
 * \li   minFrac    Min. fraction of planes included
 * \param Accum     Complex Accumulator work array
 * \param workFA    Floating work array
 * \param workCA    Complex work array
 *        Output:
 * \param   outFArrays Output arrays,planes:
 * \li      0     Peak Rotation measures (RM)
 * \li      2     PPol as amplitude squared, amplitude on return
 * \param   err     Obit Error stack 
 */
static void DirectAnaSearch(ObitFaraSyn *in, ObitCArray *Accum, ObitFArray *workFA, 
			    ObitCArray *workCA, ObitErr *err)
{
  /* Local variables */
  olong nlamb2 = in->nLamb2;
  odouble *lamb2 = in->lamb2, refLamb2 = in->refLamb2;
  ObitFArray **inQFArrays = in->inQFArrays, **inUFArrays = in->inUFArrays;
  ofloat *plnWt = in->plnWt, *QRMS = in->QRMS, *URMS = in->URMS;
  ofloat minRMSyn = in->minRMSyn;
  ofloat maxRMSyn = in->maxRMSyn;
  ofloat delRMSyn = in->delRMSyn;
  ObitFArray **outFArrays = in->outFArrays;

  olong i, iRM, nRM, nThreads=1;
  ofloat RM, cmplx[2], sumWt, norm, maxQ, maxU;
  ObitFArray *maxRM=NULL, *maxAmp2=NULL, *maxPh=NULL;
  gpointer threadArgs=NULL;
  /*gchar *routine = "DirectAnaSearch";*/

  /* 3 sigma Upper limits on RMS */
  maxQ = in->qMedian+3*in->qSigma;
  maxU = in->uMedian+3*in->uSigma;
  maxRM   = outFArrays[0]; /* alias for RM of peak amp2 */
  maxPh   = outFArrays[1]; /* alias for Peak phase */
  maxAmp2 = outFArrays[2]; /* alias for Peak amp^2 */

  /* Zero r/i of peak */
  cmplx[0] = cmplx[1] = 0.0;
  ObitCArrayFill(workCA, cmplx);

  /* Loop over RM range */
  nRM = (olong)(1.999+(maxRMSyn-minRMSyn)/delRMSyn);
  for (iRM=1; iRM<=nRM; iRM++) {
    /* Setup FT threading */
    if (threadArgs==NULL)
      nThreads = ObitCArrayMakeCAFuncArgs (in->thread, NULL, NULL, 
					   Accum, inQFArrays[0], inUFArrays[0],
					   NULL, NULL,2, 0, 0, 0, 0, 0, 0, &threadArgs);
    RM = minRMSyn + (iRM-1)*delRMSyn;
    /* Zero acccumulator */
    cmplx[0] = cmplx[1] = 0.0;
    ObitCArrayFill(Accum, cmplx);
    /* loop over input planes */
    sumWt = 0.0;
    for (i=0; i<nlamb2; i++) {
      /* Want this one ? */
      if ((plnWt[i]<=0.0) || !inQFArrays[i] || !inUFArrays[i] || 
	  (QRMS[i]>maxQ) || (URMS[i]>maxU)) continue;
      sincosf((ofloat)(-2.0*RM*(lamb2[i]-refLamb2)), &cmplx[1], &cmplx[0]); /* sin/cos factors */
      /* Rotate, accumulate this plane */
      ObitCArraySMulAccumTh (inQFArrays[i], inUFArrays[i], cmplx, Accum,
			     nThreads, threadArgs);
      sumWt += plnWt[i];
    } /* Loop over lamb2 planes */
    /* Normalize - no threading conflict with ObitFArraySMul */
    norm = 1.0 / sumWt;
    ObitCArraySMul(Accum, norm);  /* This needed? */
    
    /* Look for maximum Amp^2 */
    ObitCArrayAmp2 (Accum, workFA); /* Get amplitude squared */
    /* Save maxAmp2, RM, r/i */
    ObitCArrayMaxSetValues(workFA, maxAmp2, RM, maxRM, Accum, workCA); 
  } /* end loop over RM */
  /* Shutdown FT threading */
  if (threadArgs) {
    ObitCArrayKillCAFuncArgs (nThreads, threadArgs); 
    threadArgs = NULL;
  }
  /* Get phase */
  ObitCArrayPhase (workCA, maxPh);
  ObitFArraySMul(maxPh, 0.5);  /* times 0.5 to get EVPA */
} /* end DirectAnaSearch */

/**
 *  Write output for RMAna
 * \param   in         Input object using members:
 *        Input:
 * \li   nOut       Number of planes in output.
 * \li   doError    If true do error analysis
 *        Output:
 * \li   outFitRM   Output Faraday analysis cube
 * \li   outFArrays Output arrays,planes:
 * \li      0     Peak Rotation measures (RM)
 * \li      1     EVPA
 * \li      2     PPol as amplitude squared, amplitude on return
 * \li      3     if doError, the error estimates on RM
 * \li      4     if doError, the error estimates on PPol
 * \li      nOut  Chi^2 if not doRMSyn
 * \param   err        Obit Error stack 
 */
static void WriteRMAna (ObitFaraSyn *in, ObitErr *err)
{
  /* Local variables */
  olong nOut = in->nOut;
  gboolean doError = in->doError, doRMSyn = in->doRMSyn;
  ObitImage *outFitRM = in->outFitRM;
  ObitFArray **outFArrays = in->outFArrays;

  olong plane[5] = {1,1,1,1,1};
  ObitIOCode retCode;
  gchar *routine = "WriteRMAna";
  
  /* RM Plane */
  plane[0] = 1;
  retCode = ObitImagePutPlane (outFitRM, outFArrays[0]->array, plane, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_msg (err, routine, outFitRM->name);
    return;
  }
  /* EVPA plane */
  plane[0] = 2;
  retCode = ObitImagePutPlane (outFitRM, outFArrays[1]->array, plane, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_msg (err, routine, outFitRM->name);
    return;
  }
  /* PPol plane */
  if (doRMSyn) ObitFArraySqrt(outFArrays[2]); /* get amp plane for RMSyn */
  plane[0] = 3;
  retCode = ObitImagePutPlane (outFitRM, outFArrays[2]->array, plane, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_msg (err, routine, outFitRM->name);
    return;
  }
  /* Giving errors? */
  if (doError && (!doRMSyn)) {
    /* RM Error Plane */
    plane[0] = 4;
    retCode = ObitImagePutPlane (outFitRM, outFArrays[3]->array, plane, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) {
      Obit_traceback_msg (err, routine, outFitRM->name);
      return;
    }
    /* EVPA Error plane */
    plane[0] = 5;
    retCode = ObitImagePutPlane (outFitRM, outFArrays[4]->array, plane, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) {
      Obit_traceback_msg (err, routine, outFitRM->name);
      return;
    }
  } /* end doError */
  /* Chi^2 plane - always last */
  if (nOut>3) {
    plane[0] = nOut;
    retCode = ObitImagePutPlane (outFitRM, outFArrays[nOut-1]->array, plane, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) {
      Obit_traceback_msg (err, routine, outFitRM->name);
      return;
    }
  } /* end if need Chi^2 */
} /* end WriteRMAna  */

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitFaraSynInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFaraSyn *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread    = newObitThread();
  in->info      = newObitInfoList(); 
  in->nLamb2    = 0;
  in->lamb2     = NULL;
  in->refLamb2  = 1.0e-6;
  in->inQFArrays= NULL;  in->inUFArrays= NULL; in->outFArrays = NULL;
  in->plnWt     = NULL;
  in->QRMS      = NULL;  in->URMS      = NULL;
  in->nOut      = 0;
  in->qMedian   = -1.0;  in->qSigma = -1.0;
  in->uMedian   = -1.0;  in->uSigma = -1.0;
  in->minRMSyn  = 0.0;   in->maxRMSyn  = 0.0;  in->delRMSyn  = 0.0;
  in->RMPlanes  = NULL;
  in->outAImage = NULL;  in->outPImage = NULL;
  in->outFitRM  = NULL;
  in->doWrite   = FALSE;
  in->minQUSNR  = 3.0;
  in->minFrac   = 0.5;
  in->doError   = FALSE;
  in->doRMSyn   = FALSE;
  in->maxChi2   = 10.0;
} /* end ObitFaraSynInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitFaraSyn* cast to an Obit*.
 */
void ObitFaraSynClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFaraSyn *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread = ObitThreadUnref(in->thread);
  in->info   = ObitInfoListUnref(in->info);
  if (in->lamb2) {g_free(in->lamb2 ); in->lamb2  = NULL;}
  if (in->inQFArrays) {
    for (i=0; i<in->nLamb2; i++) 
      in->inQFArrays[i] = ObitFArrayUnref(in->inQFArrays[i]);
    g_free(in->inQFArrays); in->inQFArrays = NULL;
  }
  if (in->inUFArrays) {
    for (i=0; i<in->nLamb2; i++) 
      in->inUFArrays[i] = ObitFArrayUnref(in->inUFArrays[i]);
    g_free(in->inUFArrays); in->inUFArrays = NULL;
  }
   if (in->outFArrays) {
    for (i=0; i<in->nOut; i++) ObitFArrayUnref(in->outFArrays[i]);
    g_free(in->outFArrays);
  }
 if (in->plnWt) {g_free(in->plnWt ); in->plnWt = NULL;}
  if (in->QRMS)  {g_free(in->QRMS );  in->QRMS  = NULL;}
  if (in->URMS)  {g_free(in->URMS );  in->URMS  = NULL;}
  if (in->RMPlanes) {
    for (i=0; i<in->nOut; i++) 
      in->RMPlanes[i] = ObitCArrayUnref(in->RMPlanes[i]);
    g_free(in->RMPlanes); in->RMPlanes = NULL;
  }
  in->outAImage = ObitImageUnref(in->outAImage);
  in->outPImage = ObitImageUnref(in->outPImage);
  in->outFitRM  = ObitImageUnref(in->outFitRM);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
} /* end ObitFaraSynClear */

/**
 * Does least squares RM fitting, 
 * Work divided up amoung 1 or more threads using alternating pixels
 * In each pixel, an RM is fitted if the number of valid data points 
 * exceeds in->nterm, otherwise the pixel is blanked.
 * Initial guesses and fitted values are in outFArrays:
 * \li first entry is the RM (Jy) at the reference lambda^2
 * \li second is the EVLA (rad) at 0 lambda^2
 * \li third is the polarized amplitude at 0 lambda^2
 * \li highest is Chi^2 of the fit
 * if in->doError
 * \li 4th is least squares error for RM
 * \li 5th is least squares error for EVPA
 *
 * \param  in  ObitFaraFit to fit
 * \param  err Obit error stack object.
 */
static void LSQFitter (ObitFaraSyn* in, ObitErr *err)
{
  olong i, lo, hi, nelem, nPerThread, nThreads;
  gboolean OK, more=TRUE;
  FaraFitArg **threadArgs;
  FaraFitArg *args=NULL;
  ObitThreadFunc func=(ObitThreadFunc)ThreadFaraFit;
  gchar *routine = "LSQFitter";
  /* GSL implementation */
#ifdef HAVE_GSL
  const gsl_multifit_fdfsolver_type *T=NULL;
#endif /* HAVE_GSL */ 

  /* error checks */
  if (err->error) return;

  /* Message */
  Obit_log_error(err, OBIT_InfoErr, "Doing Least Squares fitting. ");

  /* How many threads to use? */
  nThreads = MAX (1, ObitThreadNumProc(in->thread));
  /* DEBUG */
  /* nThreads = 1;  DEBUG */

  /* Initialize threadArg array  */
  threadArgs = g_malloc0(nThreads*sizeof(FaraFitArg*));
  for (i=0; i<nThreads; i++) 
    threadArgs[i] = g_malloc0(sizeof(FaraFitArg)); 

  /* Set up thread arguments */
  for (i=0; i<nThreads; i++) {
    args = (FaraFitArg*)threadArgs[i];
    args->in          = in;
    args->err         = err;
    args->nlamb2      = in->nLamb2;
    args->maxIter     = 100;     /* max. number of iterations */
    args->minDelta    = 1.0e-5;  /* Min step size */
    args->nterm       = 2;
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
    args->wrk4        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->wrk5        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->wrk6        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->wrk7        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->wrk8        = g_malloc0(args->nlamb2*sizeof(ofloat));
    args->wrk9        = g_malloc0(args->nlamb2*sizeof(ofloat));
    if (args->doError)
      args->coef      = g_malloc0((1+3*args->nterm)*sizeof(ofloat));
    else
      args->coef      = g_malloc0((2+args->nterm)*sizeof(ofloat));
   /* GSL implementation */
#ifdef HAVE_GSL
    args->solver = NULL; args->covar = NULL; args->work = NULL;
  /* Setup solver */
    T = gsl_multifit_fdfsolver_lmder;
    args->solver = gsl_multifit_fdfsolver_alloc(T, 2*args->nlamb2, 2);
    
    /* Fitting function info */
    args->funcStruc = g_malloc0(sizeof(gsl_multifit_function_fdf));
    args->funcStruc->f      = &FaraFitFunc;
    args->funcStruc->df     = &FaraFitJac;
    args->funcStruc->fdf    = &FaraFitFuncJac;
    args->funcStruc->n      = 2*args->nlamb2;
    args->funcStruc->p      = 2;
    args->funcStruc->params = args;
    
    /* Set up work arrays */
    args->covar = gsl_matrix_alloc(2, 2);
    args->work  = gsl_vector_alloc(2);
#endif /* HAVE_GSL */ 
  }
  /* end initialize */
  
  /* Threading to use a sequence of pixels */
  nelem = in->outFitRM->myDesc->inaxes[0] * in->outFitRM->myDesc->inaxes[1];
  nPerThread = 20;  /* avoid timeout */
  lo = 1;
  hi = MIN(lo+nPerThread*nThreads-1, nelem);
  
  while (more) {
    /* Set up thread arguments */
    for (i=0; i<nThreads; i++) {
      /* if (i==(nThreads-1)) hiy = in->ny;  Make sure do all */
      args = (FaraFitArg*)threadArgs[i];
      args->first  = lo+i;
      args->last   = hi;
      args->incr   = nThreads;
      if (nThreads>1) args->ithread = i;
      else args->ithread = -1;
    } /* end loop setting up threads */
    /* Update range of elements */
    lo += nPerThread*nThreads;
    lo = MIN (lo, nelem);
    hi += nPerThread*nThreads;
    hi = MIN (hi, nelem);
  
    /* Do operation possibly with threads */
    func = (ObitThreadFunc)ThreadFaraFit;
    OK = ObitThreadIterator (in->thread, nThreads, func, (gpointer)threadArgs);
  
    /* Check for problems */
    if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
    more = (lo<nelem);  /* More passes? */
 } /* end outer (more) Loop */

  /* Shut down any threading */
  ObitThreadPoolFree (in->thread);
  if (threadArgs) {
    for (i=0; i<nThreads; i++) {
      args = (FaraFitArg*)threadArgs[i];
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
      if (args->wrk4)      g_free(args->wrk4);
      if (args->wrk5)      g_free(args->wrk5);
      if (args->wrk6)      g_free(args->wrk6);
      if (args->wrk7)      g_free(args->wrk7);
      if (args->wrk8)      g_free(args->wrk8);
      if (args->wrk9)      g_free(args->wrk9);
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
} /* end LSQFitter */

/**
 * Thread function to fit a portion of the image set
 * \param arg      FaraFitArg structure
 */
static gpointer ThreadFaraFit (gpointer arg)
{
  FaraFitArg *larg      = (FaraFitArg*)arg;
  ObitFaraSyn* in       = (ObitFaraSyn*)larg->in;
  olong lo            = larg->first-1;  /* First in element range */
  olong hi            = larg->last;     /* Highest in element range */
  olong inc           = larg->incr;     /* element increment */
  olong nOut          = larg->in->nOut; /* How many output planes */
  gboolean doError    = larg->doError;  /* Error analysis? */
  ObitErr *err        = larg->err;

  olong i, j, nterm=2, nValid=0;

  /* error checks */
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Anything to do? */
  if (hi<=lo) goto doneFit;

  /* Set up frequency info - subtract reference  */
  larg->refLamb2 = in->refLamb2;
  for (i=0; i<in->nLamb2; i++) 
    larg->lamb2[i] = in->lamb2[i] - in->refLamb2;

  /* Loop over pixels */
  for (j=lo; j<hi; j+=inc) {
    /* Collect values;  */
    nValid = 0;
    for (i=0; i<in->nLamb2; i++) {
      /* Planes may not exist */
      if (in->inQFArrays[i]) larg->Qobs[i] = in->inQFArrays[i]->array[j];
      else                   larg->Qobs[i] = 0.0;
      if (in->inUFArrays[i]) larg->Uobs[i] = in->inUFArrays[i]->array[j];
      else                   larg->Uobs[i] = 0.0;
     larg->indx = j;  /* save element index */
      /* Data valid? */
      if ((larg->Qobs[i]!=0.0) && (larg->Uobs[i]!=0.0) && 
	  ((fabs(larg->Qobs[i])>in->minQUSNR*in->QRMS[i]) ||
	   (fabs(larg->Uobs[i])>in->minQUSNR*in->URMS[i]))) {
	/* Statistical weight */
	nValid++;  /* How many good values */
	larg->Qvar[i] = (in->QRMS[i]*in->QRMS[i]);
	larg->Uvar[i] = (in->URMS[i]*in->URMS[i]);
	larg->Pobs[i]    = sqrt (larg->Qobs[i]*larg->Qobs[i] + larg->Uobs[i]*larg->Uobs[i]);
	larg->Qweight[i] = 1.0 / in->QRMS[i];
	larg->Uweight[i] = 1.0 / in->URMS[i];
	  /* End if datum valid */
      } else { /* invalid pixel */
	larg->Qobs[i]    = 0.0;
	larg->Uobs[i]    = 0.0;
	larg->Qweight[i] = 0.0;
	larg->Qvar[i]    = 0.0;
	larg->Uweight[i] = 0.0;
	larg->Uvar[i]    = 0.0;
	larg->Pobs[i]    = 0.0;
      }
    } /* end loop over frequencies */
    larg->nvalid = nValid;

    /* Fit */
    FaraFit(larg);
    /* Save to output */
    if (doError) {
      for (i=0; i<nOut; i++) 
	in->outFArrays[i]->array[larg->indx]  = larg->coef[i];
    } else { /* only values */
      for (i=0; i<(1+nterm); i++) 
	in->outFArrays[i]->array[larg->indx] = larg->coef[i];
    }
    
  } /* end pixel loop */

  /* Indicate completion */
 doneFit:
  if (larg->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&larg->ithread);

  return NULL;
} /* end ThreadFaraFit */

/**
 * Fit RM and EVPA at the reference lambda^2
 * \param arg      FaraFitArg structure
 *                 fitted parameters returned in arg->in->coef
 *                 RM, EVPA, amp, sig RM, sig EVPA, chi2
 *                 sig? if arg->doError
 */
static void FaraFit (FaraFitArg *arg)
{
  ObitFaraSyn *in = arg->in;
  olong nvalid    = arg->nvalid;
  olong indx      = arg->indx;
  olong iter=0, i, nterm=2;
  ofloat sumwt, chi2, fblank = ObitMagicF();
  int status;
 
  /* Initialize output */
  if (arg->doError) 
    for (i=0; i<2*nterm+2; i++) arg->coef[i] = 0.0;
  else
    for (i=0; i<=nterm; i++) arg->coef[i]  = 0.0;

  /* Blank initial values */
  arg->coef[0] = arg->coef[2] = arg->coef[3] = fblank;
  
  if (nvalid<=MAX(2,nterm)) return;  /* enough good data for fit? */
  /* High enough fraction of valid pixels? */
  if ((((ofloat)nvalid)/((ofloat)arg->nlamb2)) < arg->in->minFrac) return;
 
  /* Get initial values from RMSyn direct search */
  arg->coef[0] = in->outFArrays[0]->array[indx];
  arg->coef[1] = in->outFArrays[1]->array[indx];

  /* Do fit */
#ifdef HAVE_GSL
  gsl_matrix *J                  = NULL;
#if HAVE_GSL2==1  /* Newer GSL*/
  J = gsl_matrix_alloc (2*arg->nlamb2, 2);
#endif /* HAVE_GSL2 */
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
  } while ((status==GSL_CONTINUE) && (iter<arg->maxIter));
  
  /* If it didn't work - bail */
  if ((status!=GSL_SUCCESS) && (status!=GSL_CONTINUE)) {
    fprintf (stderr, "Failed, status = %s\n", gsl_strerror(status));
    arg->coef[0] = arg->coef[1] = fblank; return;
  }
  
  /* Get fitted values - switch order to RM, EVPA*/
  arg->coef[0] = (ofloat)gsl_vector_get(arg->solver->x, 1);
  arg->coef[1] = (ofloat)gsl_vector_get(arg->solver->x, 0);

  /* Amplitude to arg->coef[2] */
  FaraFitAmp (arg);
  
  /* Chi^2 */
  if (nvalid>nterm) {
    sumwt = (ofloat)gsl_blas_dnrm2(arg->solver->f);
    chi2 = (sumwt*sumwt)/(nvalid-nterm);
  } else chi2 = -1.0;
  arg->coef[arg->in->nOut-1] = chi2;

  /* Errors wanted? */
  if (arg->doError) {
#if HAVE_GSL2==1  /* Newer GSL*/
    gsl_multifit_fdfsolver_jac(arg->solver, J);
#else
    J = arg->solver->J;
#endif /* HAVE_GSL2 */ 
    /* second argument removes degenerate col/row from Jacobean */
    gsl_multifit_covar (J, 1.0e-8, arg->covar);
    arg->coef[nterm+1] = sqrt(gsl_matrix_get(arg->covar, 1, 1));
    arg->coef[nterm+2] = sqrt(gsl_matrix_get(arg->covar, 0, 0));
  } /* end of get errors */
  
  /* Cleanup */
#if HAVE_GSL2==1
  if (J) gsl_matrix_free (J);
#endif /* HAVE_GSL2 */ 
#endif /* HAVE_GSL */ 
} /* end FaraFit */

#ifdef HAVE_GSL
/**
 * Function evaluator for RM fitting solver
 * Evaluates (model-observed) / sigma
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of spectral_terms
 * \param param   Function parameter structure (FaraFitArg)
 * \param f       Vector of (model-obs)/sigma for data points
 *                in order q, u each datum
 * \return completion code GSL_SUCCESS=OK
 */
static int FaraFitFunc (const gsl_vector *x, void *params, 
		      gsl_vector *f)
{
  FaraFitArg *args = (FaraFitArg*)params;
  olong nlamb2  = args->nlamb2;
  ofloat *q_func = args->wrk4;  /* local alias */
  ofloat *u_func = args->wrk5;
  double func;
  odouble RM, EVPA;
  size_t i, j;

  /* get model parameters */
  EVPA = gsl_vector_get(x, 0);
  RM   = gsl_vector_get(x, 1);

  /* First compute model phases */
  for (j=0; j<nlamb2; j++) args->wrk1[j] = 2*(EVPA + RM * args->lamb2[j]);
  /* then sine/cosine */
  ObitSinCosVec (nlamb2, args->wrk1, args->wrk3, args->wrk2);
  /* Loop over data - calculating Residuals  */
  for (i=0; i<nlamb2; i++) {
    /* Q model = args->Pobs[i] * args->wrk2[i]  */
    q_func[i] = (args->Pobs[i] * args->wrk2[i] - args->Qobs[i]) * args->Qweight[i];
    /* U model = args->Pobs[i] * args->wrk3[i]  */
    u_func[i] = (args->Pobs[i] * args->wrk3[i] - args->Uobs[i]) * args->Uweight[i];
  } /* End loop over data */
  /* Loop over data - saving Residuals */
  for (i=0; i<nlamb2; i++) {
    func = q_func[i];   gsl_vector_set(f, 2*i, func);  /* Save function residual */
    /* U model = args->Pobs[i] * args->wrk3[i]  */
    func = u_func[i];  gsl_vector_set(f, 2*i+1, func);  /* Save function residual */
  } /* End loop over data */
  
  return GSL_SUCCESS;
} /*  end FaraFitFunc */

/**
 * Jacobian evaluator for spectral fitting solver
 * Evaluates partial derivatives of model wrt each parameter
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of spectral_terms
 * \param param   Function parameter structure (FaraFitArg)
 * \param J       Jacobian matrix J[data_point, parameter]
 *                in order q, u each datum
 * \return completion code GSL_SUCCESS=OK
 */
static int FaraFitJac (const gsl_vector *x, void *params, 
		     gsl_matrix *J)
{
  FaraFitArg *args = (FaraFitArg*)params;
  olong nlamb2       = args->nlamb2;
  ofloat *evpa_q_jac = args->wrk6;  /* local alias */
  ofloat *evpa_u_jac = args->wrk7;
  ofloat *rm_q_jac   = args->wrk8;  /* local alias */
  ofloat *rm_u_jac   = args->wrk8;
  odouble RM, EVPA;
  double jac;
  size_t i, j;

  /* get model parameters */
  EVPA = gsl_vector_get(x, 0);
  RM   = gsl_vector_get(x, 1);

  /* First compute model phases */
  for (j=0; j<nlamb2; j++) args->wrk1[j] = 2*(EVPA + RM * args->lamb2[j]);
  /* then sine/cosine */
  ObitSinCosVec (nlamb2, args->wrk1, args->wrk3, args->wrk2);
  /* Loop over data - calculating gradients - bad data will have zero weight */
  for (i=0; i<nlamb2; i++) {
    /* d Q model/d EVPA = -args->Pobs[i] * args->wrk3[i]  */
    evpa_q_jac[i] = -(args->Pobs[i] * args->wrk3[i]) * args->Qweight[i];
    /* d U model/d EVPA = args->Pobs[i] * args->wrk2[i]  */
    evpa_u_jac[i] = (args->Pobs[i] * args->wrk2[i]) * args->Uweight[i];
    /* d Q model/d RM = -2 args->Pobs[i] * args->wrk3[i] * (args->lamb2[i])  */
    rm_q_jac[i] = -2*(args->Pobs[i]*args->wrk3[i]) * (args->lamb2[i])*args->Qweight[i];
    /* d U model/d RM = 2 args->Pobs[i] * args->wrk2[i] * (args->lamb2[i]) */
    rm_u_jac[i] = 2*(args->Pobs[i] * args->wrk2[i]) * (args->lamb2[i]) * args->Uweight[i];
  } /* End loop over data */
  /* Loop over data - setting gradients */
  for (i=0; i<nlamb2; i++) {
    j = 0;    /* EVPA */
    /* d Q model/d EVPA = -args->Pobs[i] * args->wrk3[i]  */
    jac = evpa_q_jac[i];
    gsl_matrix_set(J, 2*i, j, jac);  /* Save function gradient */
    /* d U model/d EVPA = args->Pobs[i] * args->wrk2[i]  */
    jac = evpa_u_jac[i];
    gsl_matrix_set(J, 2*i+1, j, jac);  /* Save function gradient */
    j = 1;   /* RM  */
    /* d Q model/d RM = -2 args->Pobs[i] * args->wrk3[i] * (args->lamb2[i])  */
    jac = rm_q_jac[i];
    gsl_matrix_set(J, 2*i, j, jac);  /* Save function gradient */
    /* d U model/d RM = 2 args->Pobs[i] * args->wrk2[i] * (args->lamb2[i]) */
    jac = rm_u_jac[i];
    gsl_matrix_set(J, 2*i+1, j, jac);  /* Save function gradient */
  } /* End loop over data */
  
  return GSL_SUCCESS;
} /*  end FaraFitJac */

/**
 * Function and Jacobian evaluator for spectral fitting solver
 * Function = (model-observed) / sigma
 * Jacobian =  partial derivatives of model wrt each parameter
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of spectral_terms
 * \param param   Function parameter structure (FaraFitArg)
 * \param f       Vector of (model-obs)/sigma for data points
 *                in order q, u each datum
 * \param J       Jacobian matrix J[data_point, parameter]
 *                in order q, u each datum
 * \return completion code GSL_SUCCESS=OK
 */
static int FaraFitFuncJac (const gsl_vector *x, void *params, 
			 gsl_vector *f, gsl_matrix *J)
{
  FaraFitArg *args = (FaraFitArg*)params;
  double func, jac;
  olong nlamb2  = args->nlamb2;
  ofloat *q_func     = args->wrk4;  /* local alias */
  ofloat *u_func     = args->wrk5;
  ofloat *evpa_q_jac = args->wrk6;  /* local alias */
  ofloat *evpa_u_jac = args->wrk7;
  ofloat *rm_q_jac   = args->wrk8;  /* local alias */
  ofloat *rm_u_jac   = args->wrk8;
  odouble RM, EVPA;
  size_t i, j;

  /* get model parameters */
  EVPA = gsl_vector_get(x, 0);
  RM   = gsl_vector_get(x, 1);

  /* First compute model phases */
  for (j=0; j<nlamb2; j++) args->wrk1[j] = 2*(EVPA + RM * args->lamb2[j]);
  /* then sine/cosine */
  ObitSinCosVec (nlamb2, args->wrk1, args->wrk3, args->wrk2);
  /* Loop over data - calculating func, gradients - bad data will have zero weight */
  for (i=0; i<nlamb2; i++) {
    /* Q model   = args->Pobs[i] * args->wrk2[i]  */
    q_func[i] = (args->Pobs[i] * args->wrk2[i] - args->Qobs[i]) * args->Qweight[i];
    /* d Q model/d EVPA = -args->Pobs[i] * args->wrk3[i]  */
    evpa_q_jac[i] = -(args->Pobs[i] * args->wrk3[i]) * args->Qweight[i];
    /* d Q model/d RM = -2 args->Pobs[i] * args->wrk3[i] * (args->lamb2[i])  */
    rm_q_jac[i] = -2*(args->Pobs[i] * args->wrk3[i]) * args->lamb2[i] * args->Qweight[i];
    /* U model   = 2 args->Pobs[i] * args->wrk3[i]  */
    u_func[i] = (args->Pobs[i] * args->wrk3[i] - args->Uobs[i]) * args->Uweight[i];
    /* d U model/d EVPA = args->Pobs[i] * args->wrk2[i]  */
    evpa_u_jac[i] = (args->Pobs[i] * args->wrk2[i]) * args->Uweight[i];
    /* d U model/d RM = args->Pobs[i] * args->wrk2[i] * (args->lamb2[i]) */
    rm_u_jac[i] = 2*(args->Pobs[i] * args->wrk2[i]) * args->lamb2[i] * args->Uweight[i];
  } /* End loop over data */
 
  /* Loop over data - setting func, gradients */
  for (i=0; i<nlamb2; i++) {
    /* Q model   = args->Pobs[i] * args->wrk2[i]  */
    func = q_func[i];
    gsl_vector_set(f, 2*i, func);  /* Save function residual */
    j = 0;    /* EVPA */
    /* d Q model/d EVPA = -args->Pobs[i] * args->wrk3[i]  */
    jac = evpa_q_jac[i];
    gsl_matrix_set(J, 2*i, j, jac);  /* Save function residual */
    j = 1;   /* RM  */
    /* d Q model/d RM = -2 args->Pobs[i] * args->wrk3[i] * (args->lamb2[i])  */
    jac = rm_q_jac[i];
    gsl_matrix_set(J, 2*i, j, jac);  /* Save function gradient */
    /* U model   = 2 args->Pobs[i] * args->wrk3[i]  */
    func = u_func[i];
    gsl_vector_set(f, 2*i+1, func);  /* Save function residual */
    j = 0;    /* EVPA */
    /* d U model/d EVPA = args->Pobs[i] * args->wrk2[i]  */
    jac = evpa_u_jac[i];
    gsl_matrix_set(J, 2*i+1, j, jac);  /* Save function gradient */
    j = 1;   /* RM  */
    /* d U model/d RM = args->Pobs[i] * args->wrk2[i] * (args->lamb2[i]) */
    jac = rm_u_jac[i];
    gsl_matrix_set(J, 2*i+1, j, jac);  /* Save function gradient */      
  } /* End loop over data */
 
  return GSL_SUCCESS;
} /*  end FaraFitFuncJac */
#endif /* HAVE_GSL */ 
#include "ObitVecFunc.h"
/**
 * Determine the polarized amplitude corresponding to a solution
 * \param arg      FaraFitArg structure
 *                 fitted parameters in arg->in->coef
 * \return unwraped polarized amplitude in arg->coef[2]
 */
static void FaraFitAmp (FaraFitArg *arg)
{
  ofloat amp, fblank = ObitMagicF();
  olong nlamb2  = arg->nlamb2;
  odouble RM, EVPA, sumQW, sumQWt, sumUW, sumUWt;
  ofloat Qval, Uval;
  size_t i, j;

  /* get model parameters */
  EVPA = (double)arg->coef[1];
  RM   = (double)arg->coef[0];
  /* Valid? */
  arg->coef[2] = fblank; 
  if ((EVPA==fblank) || (RM==fblank)) return;
  sumQW = sumQWt = sumUW = sumUWt = 0.0;  /* Init sums */

  /* First compute model phases */
  for (j=0; j<nlamb2; j++)  arg->wrk1[j] = -2*(RM * (arg->lamb2[j]));
  /* then sine/cosine */
  ObitSinCosVec (nlamb2, arg->wrk1, arg->wrk3, arg->wrk2);

  /* Loop over data for amp  */
  for (i=0; i<nlamb2; i++) {
    /* Q model = arg->Pobs[i] * arg->wrk2[i]  */
    Qval = arg->Qobs[i]*arg->wrk2[i] - arg->Uobs[i]*arg->wrk3[i];
    sumQWt += arg->Qweight[i];
    sumQW  += Qval * arg->Qweight[i];
    /* U model = arg->Pobs[i] * arg->wrk3[i]  */
    Uval = arg->Qobs[i]*arg->wrk3[i] + arg->Uobs[i]*arg->wrk2[i];
    sumUWt += arg->Uweight[i];
    sumUW  += Uval * arg->Uweight[i];
  } /* End loop over data */

  /* Any valid data? */
  if ((sumQWt<=0.0) || (sumUWt<=0.0)) return;
  Qval = sumQW / sumQWt; Uval = sumUW / sumUWt; /* Average */
  amp = (ofloat)sqrt(Qval*Qval+Uval*Uval);

  /* Set output */
  arg->coef[2] = amp;
} /*  end FaraFitAmp */
