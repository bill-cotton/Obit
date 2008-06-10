/* $Id: ObitSpectrumFit.c,v 1.2 2008/05/06 14:03:02 bcotton Exp $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008                                               */
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
  olong i;

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

  if (out->outFArrays) {
    for (i=0; i<1+in->nterm*2; i++) out->outFArrays[i] = ObitFArrayUnref(out->outFArrays[i]);
  }
  if (in->outFArrays) {
    for (i=0; i<1+in->nterm*2; i++) out->outFArrays[i] = ObitFArrayRef(in->outFArrays[i]);
  }

  /* GSL implementation */
#ifdef HAVE_GSL
  out->X     = gsl_matrix_alloc((long)out->nfreq, (long)out->nterm);
  out->y     = gsl_vector_alloc((long)out->nfreq);
  out->w     = gsl_vector_alloc((long)out->nfreq);
  out->coef  = gsl_vector_alloc((long)out->nterm);
  out->covar = gsl_matrix_alloc((long)out->nterm, (long)out->nterm);
  out->work  = gsl_multifit_linear_alloc ((long)out->nfreq, (long)out->nterm);
#endif /* HAVE_GSL */ 
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
  olong i;

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

  if (out->outFArrays) {
    for (i=0; i<1+in->nterm*2; i++) out->outFArrays[i] = ObitFArrayUnref(out->outFArrays[i]);
  }
  if (in->outFArrays) {
    for (i=0; i<1+in->nterm*2; i++) out->outFArrays[i] = ObitFArrayRef(in->outFArrays[i]);
  }

  /* GSL implementation */
#ifdef HAVE_GSL
  out->X     = gsl_matrix_alloc((long)out->nfreq, (long)out->nterm);
  out->y     = gsl_vector_alloc((long)out->nfreq);
  out->w     = gsl_vector_alloc((long)out->nfreq);
  out->coef  = gsl_vector_alloc((long)out->nterm);
  out->covar = gsl_matrix_alloc((long)out->nterm, (long)out->nterm);
  out->work  = gsl_multifit_linear_alloc ((long)out->nfreq, (long)out->nterm);
#endif /* HAVE_GSL */ 

} /* end ObitSpectrumFitClone */

/**
 * Creates an ObitSpectrumFit 
 * \param name   An optional name for the object.
 * \param nterm  Number of coefficients of powers of log10(nu)
 * \return the new object.
 */
ObitSpectrumFit* ObitSpectrumFitCreate (gchar* name, olong nterm)
{
  ObitSpectrumFit* out;
  olong i;

  /* Create basic structure */
  out = newObitSpectrumFit (name);
  out->nterm = nterm;

  /* Define term arrays */
  out->outFArrays = g_malloc0((1+out->nterm*2)*sizeof(ObitFArray*));
  for (i=0; i<1+out->nterm*2; i++) out->outFArrays[i] = NULL;

  return out;
} /* end ObitSpectrumFitCreate */

/**
 * Fit spectra to an image cube pixels.
 * \param in       Spectral fitting object
 *                 Potential parameters on in->info:
 * \li "minSNR" OBIT_float scalar Minimum SNR for fitting spectrum [def 3.0]
 * \li "calFract" OBIT_float (?,1,1) Calibration error as fraction of flux
 *              One per frequency or one for all, def 0.05
 * \li "PBmin"  OBIT_float (?,1,1) Minimum beam gain correction
 *              One per frequency or one for all, def 0.05, 1.0 => no gain corrections
 * \li "antSize" OBIT_float (?,1,1) Antenna diameter (m) for gain corr, 
 *              One per frequency or one for all, def 25.0
 *
 * \param inImage  Image cube to be fitted
 * \param outImage Image cube with fitted spectra.
 *                 Should be defined but not created.
 *                 Planes 1->nterm are coefficients per pixel
 *                 Planes nterm+1->2*nterm are uncertainties in coefficients
 *                 Plane 2*nterm+1 = Chi squared of fit
 * \param err      Obit error stack object.
 */
void ObitSpectrumFitCube (ObitSpectrumFit* in, ObitImage *inImage, 
			  ObitImage *outImage, ObitErr *err)
{
  olong i, j, iplane;
  olong naxis[2];
  ObitIOSize IOBy;
  ObitInfoType type;
  ObitIOCode retCode;
  union ObitInfoListEquiv InfoReal; 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1}, PBdim[MAXINFOELEMDIM], ASdim[MAXINFOELEMDIM];
  ofloat *calFract, *PBmin, *antSize, pbmin, antsize;
  odouble refFreq, *freqs=NULL;
  double arg, xij;
  gboolean doGain;
  gchar *today=NULL;
  gchar *routine = "ObitSpectrumFitCube";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert(ObitImageIsA(inImage));
  g_assert(ObitImageIsA(outImage));

  /* Control parameters */
  /* Min SNR for fit */
  InfoReal.flt = 3.0; type = OBIT_float;
  in->minSNR = 3.0;
  ObitInfoListGetTest(in->info, "minSNR", &type, dim, &InfoReal);
  if (type==OBIT_float) in->minSNR = InfoReal.flt;
  else if (type==OBIT_double) in->minSNR = (ofloat)InfoReal.dbl;

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

  /* Determine number of frequency planes and initialize in */
  in->nfreq = inImage->myDesc->inaxes[inImage->myDesc->jlocf];
  in->BeamShapes = g_malloc0(in->nfreq*sizeof(ObitBeamShape*));
  for (i=0; i<in->nfreq; i++) in->BeamShapes[i] = NULL;
  in->RMS        = g_malloc0(in->nfreq*sizeof(ofloat));
  in->calFract   = g_malloc0(in->nfreq*sizeof(ofloat));
  in->inFArrays  = g_malloc0(in->nfreq*sizeof(ObitFArray*));
  freqs          = g_malloc0(in->nfreq*sizeof(odouble));
  /* GSL implementation */
#ifdef HAVE_GSL
  in->X     = gsl_matrix_alloc((long)in->nfreq, (long)in->nterm);
  in->y     = gsl_vector_alloc((long)in->nfreq);
  in->w     = gsl_vector_alloc((long)in->nfreq);
  in->coef  = gsl_vector_alloc((long)in->nterm);
  in->covar = gsl_matrix_alloc((long)in->nterm, (long)in->nterm);
  in->work  = gsl_multifit_linear_alloc ((long)in->nfreq, (long)in->nterm);
#endif /* HAVE_GSL */ 

  /* Image size */
  in->nx = inImage->myDesc->inaxes[0];
  in->ny = inImage->myDesc->inaxes[1];
  naxis[0] = (olong)in->nx;  naxis[1] = (olong)in->ny; 
  for (i=0; i<in->nfreq; i++) in->inFArrays[i]     = ObitFArrayCreate (NULL, 2, naxis);
  for (i=0; i<1+in->nterm*2; i++) in->outFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);

  /* Calibration error */
  for (i=0; i<in->nfreq; i++) in->calFract[i] = 0.05; /* default */
  if (ObitInfoListGetP(in->info, "calFract", &type, dim, (gpointer)&calFract)) {
    if (dim[0]>=in->nfreq) for (i=0; i<in->nfreq; i++) in->calFract[i] = calFract[i];
    else for (i=0; i<in->nfreq; i++) in->calFract[i] = calFract[0];
  }

  /* Output Image descriptor */
  in->outDesc = ObitImageDescCopy (inImage->myDesc, in->outDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  /* Frequency axis is actually fitted spectrum parameters and errors */
  in->outDesc->inaxes[in->outDesc->jlocf] =  1+in->nterm*2;
  in->outDesc->bitpix = -32;  /* Float it */
  /* Creation date today */
  today = ObitToday();
  strncpy (in->outDesc->date, today, IMLEN_VALUE);
  if (today) g_free(today);
  outImage->myDesc = ObitImageDescCopy (in->outDesc, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

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
    freqs[iplane] = inImage->myDesc->crval[inImage->myDesc->jlocf] + 
      inImage->myDesc->cdelt[inImage->myDesc->jlocf] * 
      (inImage->myDesc->crpix[inImage->myDesc->jlocf] - inImage->myDesc->plane);
    
  } /* end loop reading planes */

  /* Close input */
  retCode = ObitImageClose (inImage, err);
  inImage->extBuffer = FALSE;   /* May need I/O buffer later */
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inImage->name);

  /* Average Frequency */
  refFreq = 0.0;
  for (i=0; i<in->nfreq; i++) refFreq += freqs[i];
  refFreq /= in->nfreq;

  /* Set fitting matrix */
  /* GSL implementation */
#ifdef HAVE_GSL
  for (i=0; i<in->nfreq; i++) {
    xij = 1.0;
    arg = log10(freqs[i]/refFreq);
    for (j=0; j<in->nterm; j++) {
      gsl_matrix_set(in->X, i, j, xij);
      xij *= arg;
    }
  }
#endif /* HAVE_GSL */ 

  /* Do actual fitting */
  Fitter (in, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Update output header to reference Frequency */
  outImage->myDesc->crval[outImage->myDesc->jlocf] = refFreq;

  /* Write output */
  WriteOutput(in, outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  if (freqs) g_free(freqs);  /* Cleanup */

} /* end ObitSpectrumFitCube */

/**
 * Fit spectra to an array of images
 * \param in       Spectral fitting object
 * \param nimage   Number of entries in imArr
 * \param imArr    Array of images to be fitted
 * \param outImage Image cube with fitted spectra.
 *                 Should be defined but not created.
 *                 Planes 1->nterm are coefficients per pixel
 *                 Planes nterm+1->2*nterm are uncertainties in coefficients
 *                 Plane 2*nterm+1 = Chi squared of fit
 * \param err      Obit error stack object.
 */
void ObitSpectrumFitImArr (ObitSpectrumFit* in, olong nimage, ObitImage **imArr, 
			   ObitImage *outImage, ObitErr *err)
{
  olong i, j, iplane;
  olong naxis[2];
  ObitIOSize IOBy;
  ObitInfoType type;
  ObitIOCode retCode;
  union ObitInfoListEquiv InfoReal; 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1}, PBdim[MAXINFOELEMDIM], ASdim[MAXINFOELEMDIM];
  ofloat *calFract, *PBmin, *antSize, pbmin, antsize, ipixel[2], opixel[2];
  odouble refFreq, *freqs=NULL;
  double xij, arg;
  gboolean doGain, bad;
  gchar *today=NULL;
  gchar *routine = "ObitSpectrumFitCube";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert(ObitImageIsA(outImage));

  /* Control parameters */
  /* Min SNR for fit */
   InfoReal.flt = 3.0; type = OBIT_float;
  in->minSNR = 3.0;
  ObitInfoListGetTest(in->info, "minSNR", &type, dim, &InfoReal);
  if (type==OBIT_float) in->minSNR = InfoReal.flt;
  else if (type==OBIT_double) in->minSNR = (ofloat)InfoReal.dbl;

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
  freqs          = g_malloc0(in->nfreq*sizeof(odouble));
  /* GSL implementation */
#ifdef HAVE_GSL
  in->X     = gsl_matrix_alloc((long)in->nfreq, (long)in->nterm);
  in->y     = gsl_vector_alloc((long)in->nfreq);
  in->w     = gsl_vector_alloc((long)in->nfreq);
  in->coef  = gsl_vector_alloc((long)in->nterm);
  in->covar = gsl_matrix_alloc((long)in->nterm, (long)in->nterm);
  in->work  = gsl_multifit_linear_alloc ((long)in->nfreq, (long)in->nterm);
#endif /* HAVE_GSL */ 
    
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
      for (i=0; i<in->nfreq; i++) in->inFArrays[i]      = ObitFArrayCreate (NULL, 2, naxis);
      for (i=0; i<1+in->nterm*2; i++) in->outFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);
      
      /* Calibration error */
      for (i=0; i<in->nfreq; i++) in->calFract[i] = 0.05; /* default */
      if (ObitInfoListGetP(in->info, "calFract", &type, dim, (gpointer)&calFract)) {
	if (dim[0]>=in->nfreq) for (i=0; i<in->nfreq; i++) in->calFract[i] = calFract[i];
	else for (i=0; i<in->nfreq; i++) in->calFract[i] = calFract[0];
      }
      
      /* Output Image descriptor */
      in->outDesc = ObitImageDescCopy (imArr[iplane]->myDesc, in->outDesc, err);
      if (err->error) Obit_traceback_msg (err, routine, imArr[iplane]->name);
      /* Frequency axis is actually fitted spectrum parameters and errors */
      in->outDesc->inaxes[in->outDesc->jlocf] =  1+in->nterm*2;
      in->outDesc->bitpix = -32;  /* Float it */
      /* Creation date today */
      today = ObitToday();
      strncpy (in->outDesc->date, today, IMLEN_VALUE);
      if (today) g_free(today);
      outImage->myDesc = ObitImageDescCopy (in->outDesc, outImage->myDesc, err);
      if (err->error) Obit_traceback_msg (err, routine, imArr[iplane]->name);
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
    freqs[iplane] = imArr[iplane]->myDesc->crval[imArr[iplane]->myDesc->jlocf] + 
      imArr[iplane]->myDesc->cdelt[imArr[iplane]->myDesc->jlocf] * 
      (imArr[iplane]->myDesc->crpix[imArr[iplane]->myDesc->jlocf] - imArr[iplane]->myDesc->plane);
    
    /* Close input */
    retCode = ObitImageClose (imArr[iplane], err);
    imArr[iplane]->extBuffer = FALSE;   /* May need I/O buffer later */
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) 
      Obit_traceback_msg (err, routine, imArr[iplane]->name);
  } /* end loop over input images */

  /* Average Frequency */
  refFreq = 0.0;
  for (i=0; i<in->nfreq; i++) refFreq += freqs[i];
  refFreq /= in->nfreq;

  /* Set fitting matrix */
  /* GSL implementation */
#ifdef HAVE_GSL
  for (i=0; i<in->nfreq; i++) {
    xij = 1.0;
    arg = log10(freqs[i]/refFreq);
    for (j=0; j<in->nterm; j++) {
      gsl_matrix_set(in->X, i, j, xij);
      xij *= arg;
    }
  }
#endif /* HAVE_GSL */ 

  /* Do actual fitting */
  Fitter (in, err);
  if (err->error) Obit_traceback_msg (err, routine, imArr[0]->name);

  /* Update output header to reference Frequency */
  outImage->myDesc->crval[outImage->myDesc->jlocf] = refFreq;

  /* Write output */
  WriteOutput(in, outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  if (freqs) g_free(freqs);  /* Cleanup */

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
  olong ix, iy, i, indx, nx, ny, nterm;
  olong naxis[2];
  ObitIOSize IOBy;
  ObitIOCode retCode;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitFArray **inArrays=NULL, *outArray=NULL;
  odouble refFreq, arg, aarg, sum;
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

  /* Reference frequency is on inImage */
  refFreq = inImage->myDesc->crval[inImage->myDesc->jlocf];

  /* Output Image descriptor */
  outImage->myDesc = ObitImageDescCopy (inImage->myDesc,  outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  /* Frequency axis is actually fitted spectrum parameters and errors */
  outImage->myDesc->inaxes[outImage->myDesc->jlocf] = 1;
  outImage->myDesc->crval[outImage->myDesc->jlocf] = outFreq;

  /* Determine number of frequency terms and initialize storage arrays */
  nterm = (inImage->myDesc->inaxes[inImage->myDesc->jlocf]-1)/2;

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
  arg = log10(outFreq/refFreq);
  /* Loop over pixels in Y */
  indx = -1;
  for (iy=0; iy<ny; iy++) {
    /* Loop over pixels in X */
    for (ix=0; ix<nx; ix++) {
      indx ++;
      sum = 0.0;
      aarg = arg;
      /* Sum polynomial in log10(nv/nvRef) */
      for (i=1; i<nterm; i++) {
	sum += inArrays[i]->array[indx]*aarg;
	aarg *= arg;
      }
      /* Add value at reference frequency */
      sum = pow(10.0, sum) + inArrays[0]->array[indx];
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
 * Reference frequency of fit is 1 GHz.
 * \param nfreq    Number of entries in freq, flux, sigma
 * \param nterm    Number of coefficients of powers of log10(nu) to fit
 * \param freq     Frequency (Hz)
 * \param flux     Flux (Jy)
 * \param sigma    Errors (Jy)
 * \param err      Obit error stack object.
 * \return  Array of fitter parameters, errors for each and Chi Squares of fit
 *          Initial terms are in Jy, other in log10.
 */
gfloat* ObitSpectrumFitSingle (gint nfreq, olong nterm, odouble *freq, 
			    ofloat *flux, ofloat *sigma, ObitErr *err)
{
  ofloat *out = NULL;
  ofloat wt, fblank = ObitMagicF();
  double arg, xij, chisq;
  olong i, j;

  /* GSL implementation */
#ifdef HAVE_GSL
  gsl_matrix *X;
  gsl_vector *y;
  gsl_vector *w;
  gsl_vector *coef;
  gsl_matrix *covar;
  gsl_multifit_linear_workspace *work;
#else
  gchar *routine = "ObitSpectrumFitSingle";
#endif /* HAVE_GSL */

  if (err->error) return out;

  /* GSL implementation */
#ifdef HAVE_GSL
  X     = gsl_matrix_alloc((long)nfreq, (long)nterm);
  y     = gsl_vector_alloc((long)nfreq);
  w     = gsl_vector_alloc((long)nfreq);
  coef  = gsl_vector_alloc((long)nterm);
  covar = gsl_matrix_alloc((long)nterm, (long)nterm);
  work  = gsl_multifit_linear_alloc ((long)nfreq, (long)nterm);
  out   = g_malloc0((1+nterm*2)*sizeof(ofloat));

  /* Fill fitting arrays */
  for (i=0; i<nfreq; i++) {
    xij = 1.0;
    arg = log10(freq[i]*1.0e-9);
    for (j=0; j<nterm; j++) {
      gsl_matrix_set(X, i, j, xij);
      xij *= arg;
    }
  }

  /* Fill values for y, w */
  for (j=0; j<nfreq; j++) {
    if (flux[j]!=fblank) { /* Data valid */
      gsl_vector_set(y, j, (double)log10(flux[j]));
      wt = 1.0;
      if (sigma[j]>0.0) wt = 1.0 / (sigma[j]*sigma[j]);
      gsl_vector_set(w, j, (double)wt);
    } else {
      gsl_vector_set(y, j, 0.0);
      gsl_vector_set(w, j, 0.0);
    }
  }
  
  /* Fit */
  gsl_multifit_wlinear (X, w, y, coef, covar, &chisq, work);
  
  /* get results */
  for (j=0; j<nterm; j++) out[j] = (ofloat)gsl_vector_get(coef, j);
  /* Errors from diagonal terms of covarance */
  for (j=0; j<nterm; j++) out[nterm+j] = (ofloat)gsl_matrix_get(covar, j, j);
  /* First terms in Jy/bm rather than log */
  if (fabs(out[0]<20.0)) out[0] = (ofloat)pow(10.0, out[0]);
  else out[0] = fblank;
  if (fabs(out[nterm]<20.0)) out[nterm] = (ofloat)pow(10.0, out[nterm]);
  else out[nterm] = fblank;
  out[nterm*2]    = (ofloat)chisq;
  
  /* Cleanup */
  if (X)     gsl_matrix_free(X);
  if (y)     gsl_vector_free(y);
  if (w)     gsl_vector_free(w);
  if (coef)  gsl_vector_free(coef);
  if (covar) gsl_matrix_free(covar);
  if (work)  gsl_multifit_linear_free (work);

#else  /* No implementation */
  Obit_log_error(err, OBIT_Error, 
		 "%s: Not yet implemented", 
		 routine);
  return NULL;
#endif /* HAVE_GSL */

  return out;
} /* end ObitSpectrumFitSingle */

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
  in->minSNR     = 3.0;
  in->RMS        = NULL;
  in->calFract   = NULL;
  in->outDesc    = NULL;
  in->inFArrays  = NULL;
  in->BeamShapes = NULL;
  in->outFArrays = NULL;
#ifdef HAVE_GSL
  in->X     = NULL;
  in->y     = NULL;
  in->w     = NULL;
  in->coef  = NULL;
  in->covar = NULL;
  in->work  = NULL;
#endif /* HAVE_GSL */ 

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
  olong i;

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

  if (in->outFArrays) {
    for (i=0; i<1+in->nterm*2; i++) in->outFArrays[i] = ObitFArrayUnref(in->outFArrays[i]);
    g_free(in->outFArrays);
  }

#ifdef HAVE_GSL
  if (in->X)     gsl_matrix_free(in->X);
  if (in->y)     gsl_vector_free(in->y);
  if (in->w)     gsl_vector_free(in->w);
  if (in->coef)  gsl_vector_free(in->coef);
  if (in->covar) gsl_matrix_free(in->covar);
  if (in->work)  gsl_multifit_linear_free (in->work);
  in->X     = NULL;
  in->y     = NULL;
  in->w     = NULL;
  in->coef  = NULL;
  in->covar = NULL;
  in->work  = NULL;
#endif /* HAVE_GSL */ 

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSpectrumFitClear */

/**
 * Does spectral fitting, input and output on input object
 * In each pixel determines weighted average and RMS and if the SNR
 * exceeds minSNR then the spectrum is fitted. 
 * The pixel value at the reference frequency is either the average pixel 
 * if the SNR is too low or 10^(first term) of fit.
 * The results are stored in outFArrays:
 * \li first entry is the pixel value (Jy/bm) at the reference freq
 * \li entries 1->nterm-1 are the coefficients of log10(freq/refFreq)^n
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
  olong ix, iy, indx, i, j, count;
  odouble sum1, sum2, Angle, pos[2];
  double chisq;
  ofloat *val=NULL, *wt=NULL, *terms=NULL, pbfact, avg, sigma, wrms, pixel[2];
  ofloat fblank = ObitMagicF();
  ObitBeamShapeClassInfo *BSClass;
  gchar *routine = "Fitter ";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Work arrays */
  val   = g_malloc0(in->nfreq*sizeof(ofloat));
  wt    = g_malloc0(in->nfreq*sizeof(ofloat));
  terms = g_malloc0((1+in->nterm*2)*sizeof(ofloat));

  /* Loop over pixels in Y */
  indx = -1;
  for (iy=0; iy<in->ny; iy++) {
    /* Loop over pixels in X */
    for (ix=0; ix<in->nx; ix++) {
      indx ++;

      /* Distance from Center */
      pixel[0] = (ofloat)(ix+1.0); pixel[1] = (ofloat)(iy+1.0);
      ObitImageDescGetPos(in->outDesc, pixel, pos, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      BSClass = (ObitBeamShapeClassInfo*)(in->BeamShapes[0]->ClassInfo);
      Angle = BSClass->ObitBeamShapeAngle(in->BeamShapes[0], pos[0], pos[1], 0.0);
      /* Collect values; get weighted average */
      sum1 = sum2 = 0.0; count = 0;
      for (i=0; i<in->nfreq; i++) {
	val[i] = in->inFArrays[i]->array[indx];
	  if (val[i]!=fblank) {
	    /* Statistical weight */
	    wt[i] = 1.0 / (in->RMS[i]*in->RMS[i] + in->calFract[i]*in->calFract[i]*val[i]*val[i]);
	    /* Primary beam correction */
	    BSClass = (ObitBeamShapeClassInfo*)(in->BeamShapes[i]->ClassInfo);
	    pbfact  = BSClass->ObitBeamShapeGainSym(in->BeamShapes[i], Angle);
	    val[i] /= pbfact;
	    wt[i]  *= pbfact*pbfact;
	    count++;
	    sum1 += val[i]*wt[i];
	    sum2 += wt[i];
	    /* End if datum valid */
	  } else { /* invalid pixel */
	    wt[i] *= 0.0;
	  }
      } /* end loop over frequencies */
      
      /* Weighted average */
      if (count>0) {
	avg = sum1/sum2;
	sigma = 1.0 / sqrt(sum2);
      } else { /* No pixels */
	avg = fblank;
	sigma = 1.0e20;
      }
      
      /* Weighted RMS */
      if (avg!=fblank) {
	sum1 = sum2 = 0.0; count = 0;
	for (i=0; i<in->nfreq; i++) {
	  if (val[i]!=fblank) {
	    count++;
	    sum1 += (val[i]-avg)*(val[i]-avg)*wt[i];
	    sum2 += wt[i];
	  }
	} /* end loop over frequencies */
	if (count>1) wrms = sqrt (sum1/sum2);
	else wrms = fblank;
      } else {  /* No valid data */
	wrms = fblank;
      } /* end get pixel sigma */

      /* Anything worth fitting in this pixel? */
      if ((count>=in->nterm) && (avg!=fblank) && ((avg/sigma)>in->minSNR)) {
	
	/* GSL implementation */
#ifdef HAVE_GSL

	/* Fill values for y, w */
	for (j=0; j<in->nfreq; j++) {
	  if ((val[j]!=fblank) && (val[j]>sigma*0.1)) { /* Data valid */
	    gsl_vector_set(in->y, j, (double)log10(val[j]));
	    gsl_vector_set(in->w, j, (double)wt[j]);
	  } else {
	    gsl_vector_set(in->y, j, 0.0);
	    gsl_vector_set(in->w, j, 0.0);
	  }
	}

	/* Fit */
	gsl_multifit_wlinear (in->X, in->w, in->y, in->coef, in->covar, &chisq, in->work);
	
	/* get results */
	for (j=0; j<in->nterm; j++) terms[j] = (ofloat)gsl_vector_get(in->coef, j);
	/* Errors from diagonal terms of covarance */
	for (j=0; j<in->nterm; j++) terms[in->nterm+j] = (ofloat)gsl_matrix_get(in->covar, j, j);
	/* First terms in Jy/bm rather than log */
	if (fabs(terms[0]<20.0)) terms[0] = (ofloat)pow(10.0, terms[0]);
	else terms[0] = fblank;
	if (fabs(terms[in->nterm]<20.0)) terms[in->nterm] = (ofloat)pow(10.0, terms[in->nterm]);
	else terms[in->nterm] = fblank;
	terms[in->nterm*2] = (ofloat)chisq;
	
#endif /* HAVE_GSL */ 
      } /* end do fitting */
      else { /* no fitting higher terms zeroed */
	for (i=1; i<1+in->nterm*2; i++) terms[i] = 0.0;
	terms[0] = avg;
	terms[in->nterm] = wrms;
      }
      
      /* Save to output */
      for (i=0; i<1+in->nterm*2; i++) in->outFArrays[i]->array[indx] = terms[i];

    } /* end x loop */
  } /* end y loop */

  /* deallocate work arrays */
  if (val  ) g_free(val);
  if (wt)    g_free(wt);
  if (terms) g_free(terms);
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
  olong iplane;
  ObitIOSize IOBy;
  ObitIOCode retCode;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "WriteOutput";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert(ObitImageIsA(outImage));

  /* Open output */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (outImage->info, "IOBy", OBIT_long, dim, &IOBy);
  outImage->extBuffer = TRUE;   /* Using outFArrays as I/O buffer */
  retCode = ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, outImage->name);

  /* Loop writing planes */
  for (iplane=0; iplane<1+in->nterm*2; iplane++) {
    retCode = ObitImageWrite (outImage, in->outFArrays[iplane]->array, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) Obit_traceback_msg (err, routine, outImage->name);

  } /* end loop writing planes */

  /* Close output */
  retCode = ObitImageClose (outImage, err);
  outImage->extBuffer = FALSE;   /* May need I/O buffer later */
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, outImage->name);

} /* end WriteOutput */
