/* $Id$     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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

#include <math.h>
#include "ObitOTFGrid.h"
#include "ObitOTFSel.h"
#include "ObitFArrayUtil.h"
#include "ObitImageUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFGrid.c
 * ObitOTFGrid class function definitions.
 * GBT/OTF data class for gridding data into an image
 * This class is derived from the Obit base class.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitOTFGrid";

/**
 * ClassInfo structure ObitOTFGridClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitOTFGridClassInfo myClassInfo = {FALSE};

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/** Degrees to radians factor */
#ifndef DG2RAD  
#define DG2RAD G_PI / 180.0
#endif

/**  Radians to degrees factor */
#ifndef RAD2DG  
#define RAD2DG 180.0 / G_PI
#endif

/** Velocity of light */
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif
/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitOTFGridInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitOTFGridClear (gpointer in);

/** Private: Fill convolving function table */
static void ConvFunc (ObitOTFGrid* in, olong fnType, ofloat over, ofloat inParm[10]);

/** Private: Compute spherical wave functions */
static ofloat sphfn (olong ialf, olong im, olong iflag, ofloat eta);

/** Private: Set Class function pointers. */
static void ObitOTFGridClassInfoDefFn (gpointer inClass);

/** Private: Threaded prep/grid buffer */
static gpointer ThreadGridBuffer (gpointer arg);

/*---------------Private structures----------------*/
/* Gridding threaded function argument */
typedef struct {
  /* Gridding object */
  ObitOTFGrid *in;
  /* OTF data set to grid from current buffer */
  ObitOTF       *OTFin;
  /* First (1-rel) vis in otfdata buffer to process this thread */
  olong        first;
  /* Highest (1-rel) vis in otfdata buffer to process this thread  */
  olong        last;
  /* thread number , >0 -> no threading */
  olong        ithread;
  /* Temporary gridding array for thread */
  ObitFArray  *grid;
  /* Temporary weight gridding array for thread */
  ObitFArray  *wtgrid;
  /* Work arrays the size of ndetect */
  ofloat *xpos, *ypos;
} OTFGridFuncArg;

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitOTFGrid* newObitOTFGrid (gchar* name)
{
  ObitOTFGrid* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOTFGridClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitOTFGrid));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitOTFGridInit((gpointer)out);

 return out;
} /* end newObitOTFGrid */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitOTFGridGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOTFGridClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitOTFGridGetClass */

/**
 * Prepares for gridding OTF data of the type described by OTFin and
 * with derived image as described by imageOut.
 * Input data should be fully edited and calibrated.
 * The object OTFin will be opened during this call if it is not already open.
 * The output image should describe the center, size and grid spacing of the desired
 * image.
 * The gridding information should have been stored in the ObitInfoList on in:
 * \li "ConvType"  OBIT_long scalar = Convolving function type: [def=3]
 *                 0 = pillbox, 3 = Gaussian, 4 = Exp*Sinc, 5 = Spherodial wave
 * \li "ConvParm"  OBIT_float[10] = Convolving function parameters
 * \li "minWt"     OBIT_float (1,1,1) Minimum summed gridding convolution weight 
 *                 as a fraction of the maximum [def 0.01]
 * \li "Clip"      OBIT_float scalar = data values with abs. value larger are set zero weight
 *
 * Gridding convolution functions:
 * \li 0 = pillbox, 
 * \li 2 = Sinc, 
 *    Parm[0] = halfwidth in cells,
 *    Parm[1] = Expansion factor
 * \li 3 = Gaussian,
 *    Parm[0] = halfwidth in cells,[def 3.0]
 *    Parm[1] = Gaussian with as fraction or raw beam [def 1.0]
 * \li 4 = Exp*Sinc
 *    Parm[0] = halfwidth in cells, [def 2.0]
 *    Parm[1] = 1/sinc factor (cells) [def 1.55]
 *    Parm[2] = 1/exp factor (cells) [def 2.52]
 *    Parm[3] = exp power [def 2.0]
 * \li 5 = Spherodial wave
 *    Parm[0] = halfwidth in cells [def 3.0]
 *    Parm[1] = Alpha [def 5.0]
 *    Parm[2] = Expansion factor [not used]
 * \param in        Object to initialize
 * \param OTFin     Uv data object to be gridded.
 * \param imageDesc Descriptor for image to be derived.
 * \param err       ObitErr stack for reporting problems.
 */
void ObitOTFGridSetup (ObitOTFGrid *in, ObitOTF *OTFin, 
		      ObitImageDesc *imageDesc, ObitErr *err)
{
  ObitIOCode retCode;
  olong naxis[2];
  ofloat over, diam, lambda, cells_rad, temp, farr[10];
  olong i, convType;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean doCal, doCalSelect;
  ObitIOAccess access;
  gchar *routine = "ObitOTFGridSetup";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFGridIsA(in));
  g_assert (ObitOTFIsA(OTFin));
  g_assert (ObitImageDescIsA(imageDesc));
  if ((imageDesc->inaxes[0]<0) || (imageDesc->inaxes[1]<0)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: MUST fully define image descriptor %s",
		   routine,imageDesc->name);
    return;
  }

  /* get gridding information */
  /* minimum summed gridding weight */
  temp = 0.0001;
  ObitInfoListGetTest(in->info, "minWt", &type, dim, &temp);
  in->minWt = temp;

  /* Clipping level */
  temp = 1.0e20;
  ObitInfoListGetTest(in->info, "Clip", &type, dim, &temp);
  in->clip = temp;

  /* Calibration requested? */
  doCal = FALSE;
  ObitInfoListGetTest(OTFin->info, "doCalib", &type, dim, &doCal);
  doCalSelect = FALSE;
  ObitInfoListGetTest(OTFin->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCal || doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

 /* open OTF data to fully instantiate if not already open */
  if (in->myStatus==OBIT_Inactive) {
    retCode = ObitOTFOpen (OTFin, access, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Get imaging parameters from imageDesc */
  in->nxImage = imageDesc->inaxes[0];
  in->nyImage = imageDesc->inaxes[1];
  in->icenxImage = in->nxImage/2;
  in->icenyImage = in->nyImage/2;

  /* Center */
  in->raProj  = imageDesc->crval[imageDesc->jlocr];
  in->decProj = imageDesc->crval[imageDesc->jlocd];

  /* Beam size (as Gaussian sigma) in pixels - assume round */
  in->beamSize = OTFin->myDesc->beamSize / (2.355 * fabs (imageDesc->cdelt[imageDesc->jlocr]));

  /* Projection code */
  in->Proj = ObitOTFSkyModelProj (&imageDesc->ctype[imageDesc->jlocr][4]);

  /* create/resize grids as needed */
  naxis[0] = in->nxImage;
  naxis[1] = in->nyImage;

  if (in->grid==NULL) in->grid = ObitFArrayCreate ("OTF Grid", 2, naxis);
  /* reallocate if need be, zero in any case */
  else in->grid = ObitFArrayRealloc (in->grid, 2, naxis);

  if (in->gridWt==NULL) in->gridWt = ObitFArrayCreate ("OTF Wt Grid", 2, naxis);
  /* reallocate if need be, zero in any case */
  else in->gridWt = ObitFArrayRealloc (in->gridWt, 2, naxis);

   /* Scaling to cells */
  in->XScale = 1.0 / imageDesc->cdelt[0];
  in->YScale = 1.0 / imageDesc->cdelt[1];

  /* Allocate working arrays for detector positions */
  in->xpos = g_realloc(in->xpos, OTFin->geom->numberDetect*sizeof(ofloat));
  in->ypos = g_realloc(in->ypos, OTFin->geom->numberDetect*sizeof(ofloat));

  /* Oversampling factor in image plane relative to 1D Nyquist sampling. */
  diam = OTFin->myDesc->diameter;  /* telescope diameter */
  if (diam<1.0) diam = 100.0; /* default = 100 m */
  lambda = VELIGHT / OTFin->myDesc->crval[OTFin->myDesc->jlocf]; /* wavelength in meters */
  cells_rad = DG2RAD * fabs (imageDesc->cdelt[0]); /* cell spacing in image in radians */
  over = (lambda / (2.0 * diam)) / cells_rad;  /*??Over sampling factor */

  /* initialize convolving function table */
  /* pilbox (0) for testing (3 = Gaussian squared, 4=exp*sinc, 5=Spherodial wave) */

  /* User specified convoling parameters */
  for (i=0; i<10; i++) farr[i] = 0.0;  /* Default 0.0 */
  ObitInfoListGetTest(OTFin->info, "ConvParm", &type, dim, farr);
  for (i=0; i<10; i++)  in->convParm[i] = farr[i];
  convType = 3;  /* default - Gaussian */
  ObitInfoListGetTest(OTFin->info, "ConvType", &type, dim, &convType);
  in->convType = convType;

  ConvFunc(in, convType, over, in->convParm);
  /* fprintf (stderr," Oversampling image by %f %f %f\n",
	   over,206265.0*lambda / (2.0 * diam),  206265.0*cells_rad);DEBUG */

 }  /* end ObitOTFGridSetup */

/**
 * Read a OTF data object and accumulate to grid.
 * Buffering of data will use the buffers as defined on OTFin 
 * ("nRecPIO" in info member).
 * The OTFin object will be closed at the termination of this routine.
 * Requires setup by #ObitOTFGridCreate.
 * \param in      Gridding Object
 * \param OTFin   OTF data object to be gridded.
 *                Should be the same as passed to previous call to 
 *                #ObitOTFGridSetup for input in.
 * \param err     ObitErr stack for reporting problems.
 */
void ObitOTFGridReadOTF (ObitOTFGrid *in, ObitOTF *OTFin, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  ObitThreadFunc func=(ObitThreadFunc)ThreadGridBuffer;
  OTFGridFuncArg *args;
  olong i, nrec, lorec, hirec, nrecPerThread, nThreads, ndetect;
  gboolean done, OK;
  olong count;
  gchar *routine = "ObitOTFGridReadOTF";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFGridIsA(in));
  g_assert (ObitOTFIsA(OTFin));
  g_assert (ObitOTFDescIsA(OTFin->myDesc));

 /* OTFin should have been opened in  ObitOTFGridSetup */
  ndetect = OTFin->geom->numberDetect;    /* How many detectors */

  /* How many threads? */
  in->nThreads = MAX (1, ObitThreadNumProc(in->thread));

  /* Initialize threadArg array  */
  if (in->threadArgs==NULL) {
    in->threadArgs = g_malloc0(in->nThreads*sizeof(OTFGridFuncArg*));
    for (i=0; i<in->nThreads; i++) 
      in->threadArgs[i] = g_malloc0(sizeof(OTFGridFuncArg)); 
  } 
  
  /* Set up thread arguments */
  for (i=0; i<in->nThreads; i++) {
    args = (OTFGridFuncArg*)in->threadArgs[i];
    args->in    = in;
    args->OTFin = OTFin;
    args->xpos  = g_malloc0(ndetect*sizeof(ofloat));
    args->ypos  = g_malloc0(ndetect*sizeof(ofloat));
    if (i>0) {
      /* Need new zeroed arrays */
      args->grid = ObitFArrayCreate("Temp grid", in->grid->ndim,  in->grid->naxis);
      ObitFArrayFill (args->grid, 0.0);
      args->wtgrid = ObitFArrayCreate("Temp grid", in->gridWt->ndim,  in->gridWt->naxis);
      ObitFArrayFill (args->wtgrid, 0.0);
    } else {  /* Reference will do */
      args->grid   = ObitFArrayRef(in->grid);
      args->wtgrid = ObitFArrayRef(in->gridWt);
    }
  }
  /* end initialize */

  /* loop gridding data */
  done = (retCode != OBIT_IO_OK);
  count = 0;
  while (!done) {

    /* read buffer - applying calibration? */
    if ((OTFin->myIO->access == OBIT_IO_ReadOnly) || 
	(OTFin->myIO->access == OBIT_IO_ReadWrite))
      retCode = ObitOTFRead (OTFin, NULL, err);
    else
      retCode = ObitOTFReadSelect (OTFin, NULL, err);
    if (err->error) 
      Obit_traceback_msg (err, routine, in->name);
    done = (retCode == OBIT_IO_EOF); /* done? */
    count += OTFin->myDesc->numRecBuff;
    
    /* Divide up work */
    nrec = OTFin->myDesc->numRecBuff;
    if (nrec<1000) nThreads = 1;
    else nThreads = in->nThreads;
    nrecPerThread = nrec/nThreads;
    lorec = 1;
    hirec = nrecPerThread;
    hirec = MIN (hirec, nrec);

    /* Set up thread arguments */
    for (i=0; i<nThreads; i++) {
      if (i==(nThreads-1)) hirec = nrec;  /* Make sure do all */
      args = (OTFGridFuncArg*)in->threadArgs[i];
      args->first  = lorec;
      args->last   = hirec;
      if (nThreads>1) args->ithread = i;
      else args->ithread = -1;
      /* Update which rec */
      lorec += nrecPerThread;
      hirec += nrecPerThread;
      hirec = MIN (hirec, nrec);
    }

    /* Do  convolve and sum to grids operation on buffer possibly with threads */
    OK = ObitThreadIterator (in->thread, nThreads, func, in->threadArgs);
    
    /* Check for problems */
    if (!OK) {
      Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
      break;
    }
  } /* end loop reading/gridding data */

  /* Close data */
  retCode = ObitOTFClose (OTFin, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Accumulate thread grids if more than one */
  if (in->nThreads>1) {
    for (i=1; i<in->nThreads; i++) {
      args = (OTFGridFuncArg*)in->threadArgs[i];
      ObitFArrayAdd(in->grid, args->grid, in->grid);
      ObitFArrayAdd(in->gridWt, args->wtgrid, in->gridWt);
    }
  } /* end accumulating grids */

  /* Shut down any threading */
  ObitThreadPoolFree (in->thread);
  if (in->threadArgs) {
    for (i=0; i<in->nThreads; i++) {
      args = (OTFGridFuncArg*)in->threadArgs[i];
      if (args->grid) ObitFArrayUnref(args->grid);
      if (args->wtgrid) ObitFArrayUnref(args->wtgrid);
      if (args->xpos) g_free(args->xpos);
      if (args->ypos) g_free(args->ypos);
      g_free(in->threadArgs[i]);
    }
    g_free(in->threadArgs);
  }
  in->threadArgs = NULL;
  in->nThreads   = 0;

  /* Make sure some data processed */
  if (count<10) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: NO data selected for imaging %s",
		   routine,OTFin->myDesc->name);
    return;
  }
} /* end ObitOTFGridReadOTF  */

/**
 * Normalize the data grid by the weight grid.
 * Requires setup by #ObitOTFGridCreate and gridding by #ObitOTFGridReadOTF.
 * \param in        Object to initialize
 * \param array     Output image array.
 * \param imageDesc Descriptor for image to be derived.
 *                  Beam size corrected for effects of gridding convolution.
 * \param err       ObitErr stack for reporting problems.
 */
void ObitOTFGridNorm (ObitOTFGrid *in, ObitFArray *array,  ObitImageDesc *imageDesc, 
		      ObitErr *err)
{
  ofloat *grid=NULL, *weight=NULL, *image=NULL, minWt, fblank = ObitMagicF();
  ofloat BMout, BMcorr, BMnat, sumWt;
  olong i, pos[5] = {0,0,0,0,0};
  gboolean doScale;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar *routine = "ObitOTFGridNorm";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFGridIsA(in));
  g_assert (ObitFArrayIsA(array));
  /* Arrays compatable */
  g_assert (in->grid->ndim == array->ndim);
  g_assert (in->grid->naxis[0] == array->naxis[0]);
  g_assert (in->grid->naxis[1] == array->naxis[1]);

  /* Want beam normalized? */
  doScale = TRUE;
  ObitInfoListGetTest(in->info, "doScale", &type, dim, &doScale);

  /* data pointers */
  grid   = ObitFArrayIndex (in->grid, pos);
  weight = ObitFArrayIndex (in->gridWt, pos);
  image  = ObitFArrayIndex (array, pos);

  /* Additional scaling to image for change of units due to the increase
     of the beam size caused by the convolution to the imaging grid */
  /* Get beamsize from header in same units as in->beamSize*/
  /* Has beam been made? */
  BMnat = in->beamSize;  /* Natural (raw) beam size in pixels as FWHM */
  if ((in->fitBeamSize<=0.0) && (imageDesc->beamMaj>0.0)) /* Need beam size from image? */
    in->fitBeamSize = imageDesc->beamMaj;
  if (in->fitBeamSize>0.0)   /* beam been fitted? */
    BMout = in->fitBeamSize /  (2.355 * fabs (imageDesc->cdelt[imageDesc->jlocr]));
  else /* Use estimates of output size */
    BMout = sqrt (BMnat*BMnat + in->addBM*in->addBM);
  if (doScale) BMcorr = (BMout*BMout) / (BMnat*BMnat);
  else  BMcorr = 1.0;
  Obit_log_error(err, OBIT_InfoErr, 
		 "Correcting image by %f for new resolution",BMcorr);

  /* Loop over grid normalizing to output array */ 
  sumWt = 0.0;
  minWt = in->minWt * ObitFArrayMax(in->gridWt, pos);  /* Min wrt max */
  minWt = MAX (0.000001, minWt);
  for (i=0; i<array->arraySize; i++) {
    if (weight[i] > minWt) { 
      sumWt += weight[i];
      image[i] =  BMcorr * grid[i] / weight[i];
    } else { /* Low weight - blank */
      image[i] = fblank;
    }
  }

  /* Check that some data got imaged */
  if (sumWt<=0.0) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: NO DATA imaged in %s", routine,imageDesc->name);
    return;
  }
} /* end ObitOTFGridNorm */

/**
 * Calculates Point source response ("Beam")
 * If given a Beam it is assumed to be the actual instrumental response, 
 * else a default Gaussian is used.  The instrumental response including
 * the convolving function is calculated and left as the "myBeam" of image.
 * \param in      Gridding object with input beam size and Convolving fn. 
 *                Info includes:
 * \li "beamNx"   OBIT_int scalar Number of "x" pixels [def 32]
 * \li "beamNy"   OBIT_int scalar Number of "y" pixels [def 32]
 * \li "doScale"  OBIT_bool scalar If true, convolve/scale beam [def TRUE]
 *
 * \param image   Image to attach beam to,  Also fills in Beam size
 * \param Beam     If non NULL use as instrumental response beam 
 * \param err     ObitErr stack for reporting problems.
 */
void ObitOTFGridMakeBeam (ObitOTFGrid* in, ObitImage *image, 
			  ObitImage *Beam, ObitErr *err)
{
  ObitFArray *rawBeam=NULL, *convFunc=NULL, *beamArray=NULL;
  ObitInfoType type;
  ofloat FWHM=0.0, fitFWHM, peak, RMS, fcenter[2], *inArray, *outArray;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong naxis[2], center[2], nxCF, nyCF, ix, iy, ia, oa, ioff, ooff;
  olong ablc[2], atrc[2];
  gboolean odd, doScale;
  olong nx, ny;
  ObitIOSize IOBy;
  olong blc[IM_MAXDIM] = {1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0};
  ObitImage *theBeam=NULL;
  gchar *tname;
  gchar *routine = "ObitOTFGridMakeBeam";

  /* Get beam size */
  nx = 32;
  ObitInfoListGetTest(in->info, "beamNx", &type, dim, &nx);
  ny = 32;
  ObitInfoListGetTest(in->info, "beamNy", &type, dim, &ny);
  doScale = TRUE;
  ObitInfoListGetTest(in->info, "doScale", &type, dim, &doScale);
  
  /* Make beam at least as big as convolving func */
  nxCF = (in->convfn->naxis[0]-1) / in->convNperCell; 
  nyCF = (in->convfn->naxis[1]-1) / in->convNperCell; 
  nx = MAX (nx, nxCF);
  ny = MAX (ny, nyCF);
  nx = ObitFFTSuggestSize(nx);  /* FFT friendly size */
  ny = ObitFFTSuggestSize(ny);
  
  /* Make FArrays */
  naxis[0] = nx; naxis[1] = ny;
  convFunc = ObitFArrayCreate ("MakeBeamCFn", 2L, naxis);
    
  /* Actual instrumental response given? */
  if (Beam!=NULL) {  /* Read Beam */
    IOBy = OBIT_IO_byPlane;
    dim[0] = 1;
    ObitInfoListPut (Beam->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
    dim[0] = 7;
    ObitInfoListPut (Beam->info, "BLC", OBIT_long, dim, (gpointer)blc, err); 
    ObitInfoListPut (Beam->info, "TRC", OBIT_long, dim, (gpointer)trc, err);
    Beam->extBuffer = FALSE;
    ObitImageOpen  (Beam, OBIT_IO_ReadOnly, err); 
    ObitImageRead  (Beam, NULL, err);
    ObitImageClose (Beam, err); 
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    /* Extract beam, Window as 0-rel */
    odd = ((2*(nx/2)) != nx);
    if (odd) ablc[0] = (olong)(Beam->myDesc->crpix[0]+0.5) - nx/2-2;
    else     ablc[0] = (olong)(Beam->myDesc->crpix[0]+0.5) - nx/2-1;
    odd = ((2*(ny/2)) != ny);
    if (odd) ablc[1] = (olong)(Beam->myDesc->crpix[1]+0.5) - ny/2-2;
    else     ablc[1] = (olong)(Beam->myDesc->crpix[1]+0.5) - ny/2-1;
    atrc[0] = (olong)(Beam->myDesc->crpix[0]+0.5) + nx/2-2;
    atrc[1] = (olong)(Beam->myDesc->crpix[1]+0.5) + ny/2-2;
    rawBeam = ObitFArraySubArr(Beam->image, ablc, atrc, err);
    if (err->error) Obit_traceback_msg (err, routine, Beam->name);
    /* Free Image array? */
    Beam->image = ObitFArrayUnref(Beam->image);
  } else { /* Use default Gaussian */
    rawBeam  = ObitFArrayCreate ("MakeBeamRaw", 2L, naxis);
    /* Gaussian the size of raw beam */
    FWHM = in->beamSize * 2.355;
    center[0] = naxis[0]/2; center[1] = naxis[1]/2;
    ObitFArray2DCGauss (rawBeam, center, FWHM);
  } /* End Use default Gaussian */

  /* DEBUG point in center
     outArray = ObitFArrayIndex (convFunc, center);
     *outArray = 1.0; */

  /* Copy center of in->convfn to convFunc */
  naxis[0] = 0; naxis[1] = 0; 
  inArray  = ObitFArrayIndex (in->convfn, naxis); 
  outArray = ObitFArrayIndex (convFunc, naxis); 
  ioff = in->convNperCell/2 + in->convfn->naxis[0]*(in->convNperCell/2); 
  ooff = (ny/2 - nyCF/2) * ny + nx/2 - nxCF/2; 
  for (iy=0; iy<nyCF; iy++) { 
    ia = ioff + iy*(in->convfn->naxis[0]*in->convNperCell); 
    oa = ooff + iy*nx; 
    for (ix=0; ix<nxCF; ix++) { 
      outArray[oa] = inArray[ia]; 
      oa++; 
      ia += in->convNperCell; 
    }  /* end x loop */
  }  /* end y loop */


  /* Convolve if needed */
  if (doScale) {
    beamArray = ObitFArrayUtilConvolve (rawBeam, convFunc, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } else { /* Keep raw Beam */
    beamArray = ObitFArrayRef(rawBeam);
  }

  /* Normalize to peak */
  peak = ObitFArrayMax (beamArray, center);
  ObitFArraySMul (beamArray, 1.0 / (MAX (0.001, peak)));

  /* Fit beam if needed */
  /* in->beamSize = beam as Gaussian in pixels */
  peak = 1.0;
  fcenter[0] = (ofloat)center[0];
  fcenter[1] = (ofloat)center[1];
  fitFWHM = in->beamSize * 2.355;
  if (doScale) {
    RMS = ObitFArrayUtilFitCGauss(beamArray, &fitFWHM, fcenter, &peak, err);
  } else {
    RMS = 0.0;
  }
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  /* If we've done this before don't give message */
  if (!image->myBeam) {
    if (doScale) {
      Obit_log_error(err, OBIT_InfoErr, "Beam fit FWHM %f asec", 
		     fitFWHM* fabs(image->myDesc->cdelt[0])*3600.0);
    } else {
      Obit_log_error(err, OBIT_InfoErr, "Using beam FWHM %f asec", 
		     fitFWHM* fabs(image->myDesc->cdelt[0])*3600.0);
    }
  }

  /* fprintf (stderr,"DEBUG beam fit FWHM %f center %f %f peak %f RMS %f\n", */
  /*	   fitFWHM* fabs(image->myDesc->cdelt[0])*3600.0, fcenter[0], fcenter[1], peak, RMS); */
  /*peak =  ObitFArraySum(beamArray); */
  /*fprintf (stderr,"DEBUG sum of normalized PSF %f\n", peak); */
  
  /* Don't decrease raw resolution 
  fitFWHM = MAX (fitFWHM, FWHM);Not really */

  /* Save fitted beam */
  image->myDesc->beamMaj = fitFWHM * fabs(image->myDesc->cdelt[0]);
  image->myDesc->beamMin = fitFWHM * fabs(image->myDesc->cdelt[0]);
  image->myDesc->beamPA  = 0.0;
  in->fitBeamSize = fitFWHM * fabs(image->myDesc->cdelt[0]);
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(image->info, "fitBeamSize", OBIT_float, dim, &in->fitBeamSize);

  /* DEBUG - use raw beam rather than convolved one 
  beamArray = ObitFArrayRef(beamArray);
  beamArray = ObitFArrayRef(rawBeam);
  *//*  end DEBUG */

  /* Create/write "Beam" */
  if (!ObitImageIsA((ObitImage*)image->myBeam)) {
    tname = g_strconcat ("Beam for ",image->name, NULL);
    theBeam= newObitImage(tname);
    image->myBeam = (Obit*)theBeam;
    g_free(tname);
    theBeam->myDesc = 
      ObitImageDescCopy (image->myDesc, theBeam->myDesc, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
    theBeam->myDesc->inaxes[0] = beamArray->naxis[0];
    theBeam->myDesc->inaxes[1] = beamArray->naxis[1];
    odd = ((2*(nx/2)) != nx);
    if (odd) theBeam->myDesc->crpix[0]  = nx * 0.5 + 2.0;
    else     theBeam->myDesc->crpix[0]  = nx * 0.5 + 1.0;
    odd = ((2*(ny/2)) != ny);
    if (odd) theBeam->myDesc->crpix[1]  = ny * 0.5 + 2.0;
    else     theBeam->myDesc->crpix[1]  = ny * 0.5 + 1.0;
    theBeam->myDesc->crval[0] = 0.0;
    theBeam->myDesc->crval[1] = 0.0;
    /* external info */
    ObitImageSetBeamName (image, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
  } else {
    /* Old beam exists */
    theBeam = (ObitImage*)image->myBeam;
  }


  /* Save it to disk */
  theBeam->extBuffer = TRUE;  /* Don't need buffer */
  ObitImageOpen  (theBeam, OBIT_IO_WriteOnly, err);
  ObitImageWrite (theBeam, beamArray->array, err);
  ObitImageClose (theBeam, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  theBeam->extBuffer = FALSE;  /* Might need later */
  
  /* DEBUG
     ObitImageUtilArray2Image ("CF.fits", 1, beamArray, err);
     if (err->error) Obit_traceback_msg (err, routine, in->name); */

  /* Cleanup */
  rawBeam   = ObitFArrayUnref(rawBeam);
  convFunc  = ObitFArrayUnref(convFunc);
  beamArray = ObitFArrayUnref(beamArray);
  
} /* end ObitOTFGridMakeBeam */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitOTFGridClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitOTFGridClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitOTFGridClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitOTFGridClassInfoDefFn (gpointer inClass)
{
  ObitOTFGridClassInfo *theClass = (ObitOTFGridClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitOTFGridClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitOTFGridClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitOTFGridGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitOTFGridClear;
  theClass->ObitInit      = (ObitInitFP)ObitOTFGridInit;
  theClass->newObit       = (newObitFP)newObitOTFGrid;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitOTFGridSetup   = (ObitOTFGridSetupFP)ObitOTFGridSetup;
  theClass->ObitOTFGridReadOTF = (ObitOTFGridReadOTFFP)ObitOTFGridReadOTF;
  theClass->ObitOTFGridNorm    = (ObitOTFGridNormFP)ObitOTFGridNorm;

} /* end ObitOTFGridClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitOTFGridInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTFGrid *in = inn;
  olong i;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread       = newObitThread();
  in->info         = newObitInfoList(); 
  in->myStatus     = OBIT_Inactive;
  in->addBM        = 0.0;
  in->grid         = NULL;
  in->gridWt       = NULL;
  in->convfn       = NULL;
  in->xpos         = NULL;
  in->ypos         = NULL;
  in->minWt        = 0.01;
  in->nThreads     = 1;
  in->threadArgs   = NULL;
  for (i=0; i<10; i++)in->convParm[i] = 0.0;

} /* end ObitOTFGridInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitOTFGrid* cast to an Obit*.
 */
void ObitOTFGridClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitOTFGrid *in = inn;
  OTFGridFuncArg *args;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  in->grid      = ObitFArrayUnref(in->grid);  
  in->gridWt    = ObitFArrayUnref(in->grid);  
  in->convfn    = ObitFArrayUnref(in->convfn);
  if (in->xpos) g_free (in->xpos); in->xpos = NULL;
  if (in->ypos) g_free (in->ypos); in->ypos = NULL;
  if (in->threadArgs) {
    for (i=0; i<in->nThreads; i++) {
      args = (OTFGridFuncArg*)in->threadArgs[i];
      if (args->grid)   ObitFArrayUnref(args->grid);
      if (args->wtgrid) ObitFArrayUnref(args->wtgrid);
      if (args->xpos) g_free(args->xpos);
      if (args->ypos) g_free(args->ypos);
      g_free(in->threadArgs[i]);
    }
    g_free(in->threadArgs);
  }
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitOTFGridClear */


/**
 * Calculates 2D circular convolving function and attaches it to in.
 * Also sets addBM member of in to the additional Gaussian equivalent 
 * beam size (pixels) caused by the convolution function.
 * Note: addBM is determined heuristically from the observed broadening
 * of a simulated point source and must be redetermined if the convolution 
 * functions are modified.  
 * Algorithm lifted from AIPS.
 * \param in      Object with table to init.
 * \param fnType  Function type
 *                \li 0 = pillbox, 
 *                \li 2 = Sinc, 
 *                   inParm[0] = halfwidth in cells,[def 3.0]
 *	 	     inParm[1] = Expansion factor [def 1.55]
 *                \li 3 = Gaussian,
 *                   inParm[0] = halfwidth in cells,[def 3.0]
 *	 	     inParm[1] = Gaussian width as fraction of raw beam [def 1.0]
 *                \li 4 = Exp*Sinc
 *                   inParm[0] = halfwidth in cells, [def 2.0]
 *	 	     inParm[1] = 1/sinc factor (cells) [def 1.55]
 *	 	     inParm[2] = 1/exp factor (cells) [def 2.52]
 *	 	     inParm[3] = exp power [def 2.0]
 *                \li 5 = Spherodial wave
 *                   inParm[0] = halfwidth in cells [def 3.0]
 *	 	     inParm[1] = Alpha [def 5.0]
 *	 	     inParm[2] = Expansion factor [not used]
 * \param over    Oversampling factor in image plane relative
 *                to 1D Nyquist sampling.
 * \param inParms fnType dependent parameters
 */
static void ConvFunc (ObitOTFGrid* in, olong fnType, ofloat over, ofloat inParm[10])
{
  ofloat parm[4]; /* default parameters */
  ofloat xinc, eta, psi, p1, p2, u, umax, iumax, *convfnp;
  olong ialf, im, imeff, size, lim, bias, naxis[2];
  olong ix, iy, iaddr;
  ofloat ux, uy, width;

  /* error checks */
  g_assert (ObitOTFGridIsA(in));
  in->addBM = 0.0;   /* Beam broadening default - none */
  
  
  /*+++++++++++++++++ Pillbox ++++++++++++++++++++++++++++++++++++++++*/
  if (fnType==0) {
    in->addBM = 1.1864/over; /* heuristic value for beam broadening in pixels */
    in->addBM = 0.691; /* heuristic value for beam broadening in pixels */
    /* set parameters */
    imeff = 2;    /* width at critical sampling */
    imeff = MAX (2, MIN (8, imeff));  
    parm[0] = imeff / 2.0;
    im = imeff * over + 0.5;   /* Actual number of cells */
    im = 2 * (im / 2) + 1;     /* make odd */
    in->convWidth    = im; /* Width of convolving kernel in cells */
    in->convNperCell = 20;    /* Number of of tabulated points per cell in convfn */
    /*fprintf (stderr," DEBUG im %d parm[0] %f\n",im,parm[0]); *//* DEBUG */
    /*fprintf (stderr," DEBUG over %f\n",over );*/ /* DEBUG */
    
    /* allocate array*/
    lim = in->convWidth * in->convNperCell + 1;
    size = lim;
    naxis[0] = naxis[1] = size;
    in->convfn = ObitFArrayUnref(in->convfn);
    in->convfn = ObitFArrayCreate (in->name, 2L, naxis);
    
    /* get pointer to memory array */
    naxis[0] = naxis[1] = 0;
    convfnp = ObitFArrayIndex (in->convfn, naxis);
    
    /* fill function */
    xinc = 1.0 / ((over*0.5)*(ofloat)in->convNperCell); 
    umax = 1.0;
    bias = lim/2;
    for (iy=0; iy<lim; iy++) {
      uy = (iy-bias) * xinc;
      iaddr = iy * size;
      for (ix=0; ix<lim; ix++) {
	ux = (ix-bias) * xinc;
	u = sqrt (ux*ux + uy*uy);
	convfnp[iaddr] = 1.0;
	if (u == umax) convfnp[iaddr] = 0.5;
	else if (u > umax)  convfnp[iaddr] = 0.0;
	iaddr++;
      } /* end loop in x */
    } /* end loop in y */
    
  } else if (fnType==2) {
    /*+++++++++++++++++ Sinc ++++++++++++++++++++++++++++++++++++++++*/
    /* set parameters */
    if (inParm[0]<=0.0) 
      parm[0] = 3.0; /* Halfwidth of convolving function at critical sampling */
    else
      parm[0] = inParm[0];

    if (inParm[1]<=0.0) {
      parm[1] = 1.55; 
    } else
      parm[1] = inParm[1];

    in->addBM = 0.00;   /* heuristic value for beam broadening in pixels */

    /* set parameters */
    imeff = 1 + (olong)(2.0*parm[0]);  /* width at critical sampling */
    /* constrain range */
    imeff = MAX (4, MIN (8, imeff));  
    parm[0] = (imeff - 1.0) / 2.0;
    im = imeff  + 0.5;   /* Actual number of cells */
    im = 2 * (im / 2) + 1;     /* make odd */
    in->convWidth    = im;     /* Width of convolving kernel in cells */
    in->convNperCell = 20;     /* Number of of tabulated points per cell in convfn */
    p1 = G_PI / (over*parm[1]);
    
    /* allocate array */
    lim = in->convWidth * in->convNperCell + 1;
    size = lim;
    naxis[0] = naxis[1] = size;
    in->convfn = ObitFArrayUnref(in->convfn);
    in->convfn = ObitFArrayCreate (in->name, 2L, naxis);
    /* get pointer to memory array */
    naxis[0] = naxis[1] = 0;
    convfnp = ObitFArrayIndex (in->convfn, naxis);
    
    /* fill function */
    bias = lim/2;
    xinc = 1.0 / ((over*0.5)*(ofloat)in->convNperCell); 
    umax = parm[0]*over; 
    for (iy=0; iy<lim; iy++) {
      uy = (iy-bias) * xinc;
      iaddr = iy * size;
      for (ix=0; ix<lim; ix++) {
	ux = (ix-bias) * xinc;
	u = sqrt (ux*ux + uy*uy);
	convfnp[iaddr] = 0.0;
	
	/* trap center */
	if (u<xinc) convfnp[iaddr] = 1.0;
	else if (u <= umax) convfnp[iaddr] =  sin(u*p1) / (u*p1);
 	iaddr++;
      } /* end loop in x */
    } /* end loop in y */
    
    
  } else if (fnType==3) {
    /*+++++++++++++++++ Gaussian  +++++++++++++++++++++++++++++++++++++*/
    /* set parameters */
    if (inParm[0]<=0.0) 
       parm[0] = 3.0;   /* Effective range in terms of 1-D Nyquist sampling */
    else
      parm[0] = inParm[0];
    if (inParm[1]<=0.0) {
      width = 1.0;
      parm[1] = in->beamSize;
    } else {
      width = inParm[1];
      parm[1] = width * in->beamSize;
    }

    /* heuristic value for beam broadening in pixels */
    in->addBM = 2.355 * parm[1];  /* Function FWHM in pixels */

    /* set parameters */
    if (parm[1] <0.001) parm[1] = 1.0; /* sanity check */
    imeff = 1 + (olong)(2.0*parm[0]);  /* width at critical sampling */
    /* constrain range */
    imeff = MAX (4, MIN (8, imeff));  
    parm[0] = (imeff - 1.0) / 2.0;
    im = imeff * over * width + 0.5;/* Actual number of cells */
    im = 2 * (im / 2) + 1;  /* make odd */
    in->convWidth    = im;  /* make it odd */
    in->convNperCell = 20;  /* Number of of tabulated points per cell in convfn */
    p1 = -4.0 / (parm[1] * parm[1]);
    p1 = -1.0 * over * over / (2.0 * parm[1] * parm[1]);  /* Gaussian */
     /* fprintf (stderr," DEBUG beamsize %f %f %f\n",in->beamSize,in->addBM,width );  DEBUG */
     /*fprintf (stderr," DEBUG oversampling %f CF width  %d\n",over, im );   DEBUG */
    
    /* allocate array*/
    lim = in->convWidth * in->convNperCell + 1;
    size = lim;
    naxis[0] = naxis[1] = size;
    in->convfn = ObitFArrayUnref(in->convfn);
    in->convfn = ObitFArrayCreate (in->name, 2L, naxis);
    
    /* get pointer to memory array */
    naxis[0] = naxis[1] = 0;
    convfnp = ObitFArrayIndex (in->convfn, naxis);
    
    /* fill function */
    bias = lim/2;
    xinc = 1.0 / ((over*0.5)*(ofloat)in->convNperCell);
    umax = parm[0]*over; 
    for (iy=0; iy<lim; iy++) {
      uy = (iy-bias) * xinc;
      iaddr = iy * size;
      for (ix=0; ix<lim; ix++) {
	ux = (ix-bias) * xinc;
	u = sqrt (ux*ux + uy*uy);
	if (u > umax)  convfnp[iaddr] = 0.0;
	else convfnp[iaddr] = exp (p1 * u * u);
	iaddr++;
      } /* end loop in x */
    } /* end loop in y */
    
  } else if (fnType==4) {
    /*+++++++++++++++++ Exp Sinc ++++++++++++++++++++++++++++++++++++++++*/
    /* set parameters */
    if (inParm[0]<=0.0) 
       parm[0] = 2.0; /* Halfwidth of convolving function at critical sampling */
    else
      parm[0] = inParm[0];

    if (inParm[1]<=0.0) 
      parm[1] = 1.55;
    else
      parm[1] = inParm[1];

    if (inParm[2]<=0.0) 
      parm[2] = 2.52;
    else
      parm[2] = inParm[2];

    if (inParm[3]<=0.0) 
      parm[3] = 2.00;
    else
      parm[3] = inParm[3];

    /* heuristic value for beam broadening in pixels */
    in->addBM = 2.355 * in->beamSize * 0.2123; 

    imeff = 1 + (olong)(2.0*parm[0]);  /* width at critical sampling */
    /* constrain range */
    imeff = MAX (4, MIN (8, imeff));  
    parm[0] = (imeff - 1.0) / 2.0;
    im = imeff * over + 0.5;  /* Actual number of cells */
    im = 2 * (im / 2) + 1;    /* make odd */
    in->convWidth    = im;    /* Width of convolving kernel in cells */
    in->convNperCell = 20;    /* Number of of tabulated points per cell in convfn */
    p1 = G_PI  / (parm[1]);
    p2 = 1.0  / (parm[2]);
    /*   fprintf (stderr," DEBUG beamsize %f %f %f\n",in->beamSize,in->addBM,width ); DEBUG */
    /*fprintf (stderr," DEBUG oversampling %f CF width  %d\n",over, im );   DEBUG */
    
    /* allocate array */
    lim = in->convWidth * in->convNperCell + 1;
    size = lim;
    naxis[0] = naxis[1] = size;
    in->convfn = ObitFArrayUnref(in->convfn);
    in->convfn = ObitFArrayCreate (in->name, 2L, naxis);
    
    /* get pointer to memory array */
    naxis[0] = naxis[1] = 0;
    convfnp = ObitFArrayIndex (in->convfn, naxis);
    
    /* fill function */
    xinc = 1.0 / ((over*0.5)*(ofloat)in->convNperCell); 
    umax = parm[0]*over; 
    bias = lim/2;
    for (iy=0; iy<lim; iy++) {
      uy = (iy-bias) * xinc;
      iaddr = iy * size;
      for (ix=0; ix<lim; ix++) {
	ux = (ix-bias) * xinc;
	u = sqrt (ux*ux + uy*uy);
	convfnp[iaddr] = 0.0;
	/* trap center */
	if (u<xinc) convfnp[iaddr] = 1.0;
	else if (u <= umax) convfnp[iaddr] =  sin(u*p1) / (u*p1) *
			      exp (-pow ((u * p2), parm[3]));
	iaddr++;
      } /* end loop in x */
    } /* end loop in y */
    
    
  } else if (fnType==5) {
    
    /*+++++++++++++++++ Spherodial wave ++++++++++++++++++++++++++++++++*/
    /* set parameters */
    if (inParm[0]<=0.0) 
      parm[0] = 3.0; /* Halfwidth of convolving function at critical sampling */
    else
      parm[0] = inParm[0];

    if (inParm[1]<=0.0) 
      parm[1] = 1.00;
    else
      parm[1] = inParm[1];

    if (inParm[2]<=0.0) 
      parm[2] = 5.0;
    else
      parm[2] = inParm[2];

    /* heuristic value for beam broadening in pixels */
    in->addBM = 1.19;   /* heuristic value for beam broadening in pixels */

    /* set parameters */
    imeff = 1 + (olong)(2.0*parm[0]); /* width at critical sampling */
    /* constrain range */
    imeff = MAX (4, MIN (8, imeff));
    parm[0] = (imeff - 1.0) / 2.0;
    im = imeff * over + 0.5;  /* Actual number of cells */
    im = 2 * (im / 2) + 1;    /* make odd */
    in->convWidth    = im;    /* Width of convolving kernel in cells */
    in->convNperCell = 20;    /* Number of of tabulated points per cell in convfn */
    
    /* allocate array */
    lim = in->convWidth * in->convNperCell + 1;
    size = lim;
    naxis[0] = naxis[1] = size;
    in->convfn = ObitFArrayUnref(in->convfn);
    in->convfn = ObitFArrayCreate (in->name, 2L, naxis);
    
    /* get pointer to memory array */
    naxis[0] = naxis[1] = 0;
    convfnp = ObitFArrayIndex (in->convfn, naxis);
    
    xinc = 1.0 / ((over*0.5)*(ofloat)in->convNperCell);
    umax = parm[0]*over; 
    iumax = 1.0 / parm[0];
    bias = lim/2;
    ialf =  (olong)parm[1];
    for (iy=0; iy<lim; iy++) {
      uy = (iy-bias) * xinc;
      iaddr = iy * size;
      for (ix=0; ix<lim; ix++) {
	ux = (ix-bias) * xinc;
	u = sqrt (ux*ux + uy*uy);
	if (u > umax)  convfnp[iaddr] = 0.0;
	else {
	  eta = u * iumax;
	  psi = sphfn (ialf, imeff, 0, eta);
	  convfnp[iaddr] = psi;
	}
	iaddr++;
      } /* end loop in x */
    } /* end loop in y */
    
  } else { /* should never get here */
    g_error("Unknown convolving function type %d",fnType);
  }

  /* DEBUG - write as image 
#include "ObitImageUtil.h"
  ObitErr *err = newObitErr();
  ObitImageUtilArray2Image ("ConvolveFunc.fits",1,in->convfn, err); */
} /* end ConvFunc */

/**
 * Obsolete
 * Calculates convolving function and attaches it to in.
 * Compute Spherodial wave function convolving function table.
 * Algorithm lifted from AIPS.
 * \param in      Object with table to init.
 * \param fnType  Function type
 *                \li 0 = pillbox, 
 *                \li 3 = Gaussian squared, 
 *                \li 4 = Exp*Sinc
 *                \li 5 = Spherodial wave
 */
static void OldConvFunc (ObitOTFGrid* in, olong fnType)
{
  ofloat parm[4]; /* default parameters */
  ofloat xinc, eta, psi, p1, p2, u, absu, umax, *convfnp;
  olong ialf, im, nmax, i, size, lim, limit, bias, naxis[1];
  /*gfloat shit[701]; DEBUG */

  /* error checks */
  g_assert (ObitOTFGridIsA(in));


  /*+++++++++++++++++ Pillbox ++++++++++++++++++++++++++++++++++++++++*/
  if (fnType==0) {
   /* set parameters */
    parm[0] = 0.5;   /* AIPS defaults */
    in->convWidth      = 3; /* Width of convolving kernel in cells */
    in->convNperCell = 100; /* Number of of tabulated points per cell in convfn */

    /* allocate array*/
    lim = in->convWidth * in->convNperCell + 1;
    size = lim;
    naxis[0] = size;
    in->convfn = ObitFArrayUnref(in->convfn);
    in->convfn = ObitFArrayCreate (in->name, 1L, naxis);

    /* get pointer to memory array */
    naxis[0] = 0;
    convfnp = ObitFArrayIndex (in->convfn, naxis);

    /* fill function */
    xinc = 1.0 / (ofloat)in->convNperCell;
    umax = parm[0];
    bias = (in->convNperCell/2) * in->convWidth;
    for (i=0; i<lim; i++) {
      u = (i-bias) * xinc;
      absu = fabs (u);
      convfnp[i] = 1.0;
      if (absu == umax) convfnp[i] = 0.5;
      else if (absu > umax)  convfnp[i] = 0.0;
    }
    
   } else if (fnType==3) {
  /*+++++++++++++++++ Gaussian squared +++++++++++++++++++++++++++++++++++++*/
    /* set parameters */
     parm[0] = 3.0*in->beamSize; /* half width of convolving fn */
     parm[1] = in->beamSize;     /* Sigma of Gaussian in cells */
     parm[2] = 0.0;
     parm[3] = 0.0;
     if (parm[1] <0.001) parm[1] = 1.0; /* sanity check */
    in->convWidth    = 1.0 + 2*parm[0]; /* Width of convolving kernel in cells */
    in->convWidth    = 1 + (2*(in->convWidth / 2));  /* make it odd */
    in->convNperCell = 100; /* Number of of tabulated points per cell in convfn */
    p1 = -4.0 / (parm[1] * parm[1]);

    /* allocate array*/
    lim = in->convWidth * in->convNperCell + 1;
    size = lim;
    naxis[0] = size;
    in->convfn = ObitFArrayUnref(in->convfn);
    in->convfn = ObitFArrayCreate (in->name, 1L, naxis);

    /* get pointer to memory array */
    naxis[0] = 0;
    convfnp = ObitFArrayIndex (in->convfn, naxis);

    /* fill function */
    bias = (in->convNperCell/2) * in->convWidth;
    xinc = 1.0 / (ofloat)in->convNperCell;
    for (i=0; i<lim; i++) {
      u = (i - bias) * xinc;
      convfnp[i] = exp (p1 * u * u);
    }
    

   } else if (fnType==4) {
  /*+++++++++++++++++ Exp Sinc ++++++++++++++++++++++++++++++++++++++++*/
    /* set parameters */
     parm[0] = 3.0;   /* AIPS defaults */
     parm[1] = 1.55;
     parm[2] = 2.52;
     parm[3] = 2.00;
    in->convWidth    = 1.5 + 2*parm[0]; /* Width of convolving kernel in cells */
    in->convNperCell = 100; /* Number of of tabulated points per cell in convfn */
    p1 = G_PI / parm[1];
    p2 = 1.0 / parm[2];

    /* allocate array*/
    lim = in->convWidth * in->convNperCell + 1;
    size = lim;
    naxis[0] = size;
    in->convfn = ObitFArrayUnref(in->convfn);
    in->convfn = ObitFArrayCreate (in->name, 1L, naxis);

    /* get pointer to memory array */
    naxis[0] = 0;
    convfnp = ObitFArrayIndex (in->convfn, naxis);

    /* fill function */
    bias = (in->convNperCell/2) * in->convWidth;
    xinc = 1.0 / (ofloat)in->convNperCell;
    umax = parm[0];
    for (i=0; i<lim; i++) {
      u = (i - lim/2 - 1) * xinc;
      absu = fabs (u);
      convfnp[i] = 0.0;

      /* trap center */
      if (absu<xinc) convfnp[i] = 1.0;
      else if (absu <= umax) convfnp[i] =  sin(u*p1) / (u*p1) *
			       exp (-pow ((absu * p2), parm[3]));
    }
    

   } else if (fnType==5) {

    /*+++++++++++++++++ Spherodial wave ++++++++++++++++++++++++++++++++*/
    /* set parameters */
    in->convWidth    = 7;   /* Width of convolving kernel in cells */
    in->convNperCell = 100; /* Number of of tabulated points per cell in convfn */
    parm[0] = in->convWidth/ 2;
    parm[1] = 1.0;
    xinc = 1.0 / ((ofloat)in->convNperCell);
    
    /* allocate array*/
    lim = in->convWidth * in->convNperCell + 1;
    size = lim;
    naxis[0] = size;
    in->convfn = ObitFArrayUnref(in->convfn);
    in->convfn = ObitFArrayCreate (in->name, 1L, naxis);
    /* get pointer to memory array */
    naxis[0] = 0;
    convfnp = ObitFArrayIndex (in->convfn, naxis);
    
    nmax = parm[0]*in->convNperCell + 0.1;
    bias = (in->convNperCell/2) * in->convWidth;
    ialf = 2.0 * parm[1] + 1.1;
    im = 2.0 * parm[0] + 0.1;
    
    /* constrain range */
    im = MAX (4, MIN (8, im));
    
    /* compute half of the (symmetric) function */
    for (i=0; i<nmax; i++) {
      eta = (float)i / (float)(nmax - 1);
      psi = sphfn (ialf, im, 0, eta);
      /* DEBUG - use pillbox
	 if (i<in->convWidth/2) psi = 1.0;
	 else psi = 0.0;  */
      convfnp[bias+i] = psi;
    }
    
    /* Fill in other half */
    limit = bias-1;
    for (i=1; i<=limit; i++) convfnp[bias-i] = convfnp[bias+i];
    
  } /* end computing convolving fn */
  else { /* should never get here */
    g_error("Unknown convolving function type %d",fnType);
  }
} /* end OldConvFunc */


/**
 * Compute Spherodial wave function convolving function table.
 * Algorithm lifted from AIPS (author F. Schwab).
 * \param ialf  Selects the weighting exponent, alpha
 *              (IALF = 1, 2, 3, 4, and 5 correspond to
 *              alpha = 0, 1/2, 1, 3/2, and 2, resp.).
 * \param im    support width (4, 5, 6, 7, or 8)
 * \param iflag Chooses whether the spheroidal function itself, 
 *              or its Fourier transform, is to be approximated.  
 *              The latter is appropriate for gridding, and the former 
 *              for the u-v plane convolution.  
 *              The two differ by a factor (1-eta**2)**alpha.  
 *              iflag less than or equal to zero chooses the function
 *              appropriate for gridding, and iflag positive chooses 
 *              its Fourier transform.
 * \param eta   Eta, as the argument of the spheroidal function, 
 *              is a variable which ranges from 0 at the center of the 
 *              convoluting function to 1 at its edge (also from 0 at 
 *              the center of the gridding correction function to unity at
 *              the edge of the map).  range [0,1].
 * \return spherical wave function. -999.99 -> input error.
 */
static ofloat sphfn (olong ialf, olong im, olong iflag, ofloat eta)
{
  float psi=0.0, eta2, x;
  olong   j, ierr;
  static ofloat alpha[5] = {0.0, 0.5, 1.0, 1.5, 2.0};
  static ofloat p4[5][5] = {
    {1.584774e-2, -1.269612e-1,  2.333851e-1, -1.636744e-1, 5.014648e-2},
    {3.101855e-2, -1.641253e-1,  2.385500e-1, -1.417069e-1, 3.773226e-2},
    {5.007900e-2, -1.971357e-1,  2.363775e-1, -1.215569e-1, 2.853104e-2},
    {7.201260e-2, -2.251580e-1,  2.293715e-1, -1.038359e-1, 2.174211e-2},
    {9.585932e-2, -2.481381e-1,  2.194469e-1, -8.862132e-2, 1.672243e-2}};
  static ofloat q4[5][2] = {
    {4.845581e-1,  7.457381e-2},  {4.514531e-1,  6.458640e-2},
    {4.228767e-1,  5.655715e-2},  {3.978515e-1,  4.997164e-2},
    {3.756999e-1,  4.448800e-2}};
  static ofloat p5[5][7] = {
    {3.722238e-3, -4.991683e-2,  1.658905e-1, -2.387240e-1, 1.877469e-1, -8.159855e-2,  3.051959e-2},  
    {8.182649e-3, -7.325459e-2,  1.945697e-1, -2.396387e-1, 1.667832e-1, -6.620786e-2,  2.224041e-2},  
    {1.466325e-2, -9.858686e-2,  2.180684e-1, -2.347118e-1, 1.464354e-1, -5.350728e-2,  1.624782e-2},  
    {2.314317e-2, -1.246383e-1,  2.362036e-1, -2.257366e-1, 1.275895e-1, -4.317874e-2,  1.193168e-2},
    {3.346886e-2, -1.503778e-1,  2.492826e-1, -2.142055e-1, 1.106482e-1, -3.486024e-2,  8.821107e-3}};
  static ofloat q5[5] = 
    {2.418820e-1,  2.291233e-1,  2.177793e-1,  2.075784e-1, 1.983358e-1};
  static ofloat p6l[5][5] = {
    {5.613913e-2, -3.019847e-1,  6.256387e-1, -6.324887e-1, 3.303194e-1},  
    {6.843713e-2, -3.342119e-1,  6.302307e-1, -5.829747e-1, 2.765700e-1},  
    {8.203343e-2, -3.644705e-1,  6.278660e-1, -5.335581e-1, 2.312756e-1},  
    {9.675562e-2, -3.922489e-1,  6.197133e-1, -4.857470e-1, 1.934013e-1},
    {1.124069e-1, -4.172349e-1,  6.069622e-1, -4.405326e-1, 1.618978e-1}};
  static ofloat q6l[5][2] = {
    {9.077644e-1,  2.535284e-1},  {8.626056e-1,  2.291400e-1}, 
    {8.212018e-1,  2.078043e-1},  {7.831755e-1,  1.890848e-1},  
    {7.481828e-1, 1.726085e-1}};
  static ofloat p6u[5][5] = {
    {8.531865e-4, -1.616105e-2,  6.888533e-2, -1.109391e-1, 7.747182e-2},  
    {2.060760e-3, -2.558954e-2,  8.595213e-2, -1.170228e-1, 7.094106e-2},  
    {4.028559e-3, -3.697768e-2,  1.021332e-1, -1.201436e-1, 6.412774e-2},  
    {6.887946e-3, -4.994202e-2,  1.168451e-1, -1.207733e-1, 5.744210e-2},
    {1.071895e-2, -6.404749e-2,  1.297386e-1, -1.194208e-1, 5.112822e-2}};
  static ofloat q6u[5][2] = {
    {1.101270e+0,  3.858544e-1},  {1.025431e+0,  3.337648e-1},
    {9.599102e-1,  2.918724e-1},  {9.025276e-1,  2.575336e-1},
    {8.517470e-1,  2.289667e-1}};
  static ofloat p7l[5][5] = {
    {2.460495e-2, -1.640964e-1,  4.340110e-1, -5.705516e-1, 4.418614e-1},  
    {3.070261e-2, -1.879546e-1,  4.565902e-1, -5.544891e-1, 3.892790e-1},  
    {3.770526e-2, -2.121608e-1,  4.746423e-1, -5.338058e-1, 3.417026e-1},  
    {4.559398e-2, -2.362670e-1,  4.881998e-1, -5.098448e-1, 2.991635e-1},
    {5.432500e-2, -2.598752e-1,  4.974791e-1, -4.837861e-1, 2.614838e-1}};
  static ofloat q7l[5][2] = {
    {1.124957e+0,  3.784976e-1},  {1.075420e+0,  3.466086e-1},
    {1.029374e+0,  3.181219e-1},  {9.865496e-1,  2.926441e-1},
    {9.466891e-1,  2.698218e-1}};
  static  ofloat p7u[5][5] = {
    {1.924318e-4, -5.044864e-3,  2.979803e-2, -6.660688e-2, 6.792268e-2},  
    {5.030909e-4, -8.639332e-3,  4.018472e-2, -7.595456e-2, 6.696215e-2},  
    {1.059406e-3, -1.343605e-2,  5.135360e-2, -8.386588e-2, 6.484517e-2},  
    {1.941904e-3, -1.943727e-2,  6.288221e-2, -9.021607e-2, 6.193000e-2},
    {3.224785e-3, -2.657664e-2,  7.438627e-2, -9.500554e-2, 5.850884e-2}};
  static ofloat q7u[5][2] = {
    {1.450730e+0,  6.578685e-1},  {1.353872e+0,  5.724332e-1}, 
    {1.269924e+0,  5.032139e-1},  {1.196177e+0,  4.460948e-1},  
    {1.130719e+0,  3.982785e-1}};
  static ofloat p8l[5][6] = {
    {1.378030e-2, -1.097846e-1,  3.625283e-1, -6.522477e-1, 6.684458e-1, -4.703556e-1},  
    {1.721632e-2, -1.274981e-1,  3.917226e-1, -6.562264e-1, 6.305859e-1, -4.067119e-1},
    {2.121871e-2, -1.461891e-1,  4.185427e-1, -6.543539e-1, 5.904660e-1, -3.507098e-1},  
    {2.580565e-2, -1.656048e-1,  4.426283e-1, -6.473472e-1, 5.494752e-1, -3.018936e-1},
    {3.098251e-2, -1.854823e-1,  4.637398e-1, -6.359482e-1, 5.086794e-1, -2.595588e-1}};
  static ofloat q8l[5][2] = {
    {1.076975e+0,  3.394154e-1},  {1.036132e+0,  3.145673e-1},
    {9.978025e-1,  2.920529e-1},  {9.617584e-1,  2.715949e-1},
    {9.278774e-1,  2.530051e-1}};
  static ofloat p8u[5][6] = {
    {4.290460e-5, -1.508077e-3,  1.233763e-2, -4.091270e-2, 6.547454e-2, -5.664203e-2},  
    {1.201008e-4, -2.778372e-3,  1.797999e-2, -5.055048e-2, 7.125083e-2, -5.469912e-2},
    {2.698511e-4, -4.628815e-3,  2.470890e-2, -6.017759e-2, 7.566434e-2, -5.202678e-2},  
    {5.259595e-4, -7.144198e-3,  3.238633e-2, -6.946769e-2, 7.873067e-2, -4.889490e-2},
    {9.255826e-4, -1.038126e-2,  4.083176e-2, -7.815954e-2, 8.054087e-2, -4.552077e-2}};
 static ofloat q8u[5][2] = {
   {1.379457e+0,  5.786953e-1},  {1.300303e+0,  5.135748e-1},
   {1.230436e+0,  4.593779e-1},  {1.168075e+0,  4.135871e-1},
   {1.111893e+0,  3.744076e-1}};

   ierr = 0;
   /*  Check inputs. */
   if ((ialf<1) || (ialf>5)) ierr = 1;
   if ((im<4) || (im>8)) ierr = 2 + 10 * ierr;
   if (fabs(eta)>1.) ierr = 3 + 10 * ierr;
   if (ierr!=0) return 0.0;
   /*  So far, so good. */
   eta2 = eta*eta;
   j = ialf-1;

   /*  Branch on support width. */
   switch (im) {
 
   case 4:   /*  Support width = 4 cells. */
     x = eta2 - 1.0;
     psi = (p4[j][0] + x * (p4[j][1] + x * (p4[j][2] + x * (p4[j][3] + x * p4[j][4])))) / 
       (1.0 + x * (q4[j][0] + x * q4[j][1]));
     break;

   case 5:   /* Support width = 5 cells. */
     x = eta2 - 1.0;
     psi = (p5[j][0] + x * (p5[j][1] + x * (p5[j][2] + x * (p5[j][3] +  x * (p5[j][4] + x * (p5[j][5] + x * p5[j][6]))))))
       / (1.0 + x * q5[j]);
     break;

  case 6: /* Support width = 6 cells. */
       if (fabs(eta)<=0.75) {
	 x = eta2 - 0.5625;
	 psi = (p6l[j][0] + x * (p6l[j][1] + x * (p6l[j][2] + x * (p6l[j][3] + x * p6l[j][4])))) / 
	   (1.0 + x * (q6l[j][0] +  x * q6l[j][1]));
       } else {
	 x = eta2 - 1.0;
	 psi = (p6u[j][0] + x * (p6u[j][1] + x * (p6u[j][2] + x * (p6u[j][3] + x * p6u[j][4])))) / 
	   (1.0 + x * (q6u[j][0] + x * q6u[j][1]));
       }
     break;
     
   case 7: /* Support width = 7 cells. */
       if (fabs(eta)<=0.775) {
         x = eta2 - 0.600625;
         psi = (p7l[j][0] + x * (p7l[j][1] + x * (p7l[j][2] + x * (p7l[j][3] + x * p7l[j][4])))) / 
	   (1.0 + x * (q7l[j][0] + x * q7l[j][1]));
       } else {
	 x = eta2 - 1.0;
	 psi = (p7u[j][0] + x * (p7u[j][1] + x * (p7u[j][2] + x * (p7u[j][3] +  x * p7u[j][4])))) / 
	   (1.0 + x * (q7u[j][0] + x * q7u[j][1]));
       }
     break;

  case 8:      /* Support width = 8 cells. */
      if (abs(eta)<=0.775) {
	x = eta2 - .600625;
         psi = (p8l[j][0] + x * (p8l[j][1] + x * (p8l[j][2] + x * (p8l[j][3] + x * (p8l[j][4] + x * p8l[j][5]))))) /
           (1.0 + x * (q8l[j][0] + x * q8l[j][1]));
      } else {
	x = eta2 - 1.0;
      psi = (p8u[j][0] + x * (p8u[j][1] + x * (p8u[j][2] + x * (p8u[j][3] + x * (p8u[j][4] + x * p8u[j][5]))))) / 
	(1.0 + x * (q8u[j][0] + x * q8u[j][1]));
      }
      break;
      
   default: /* should never get here */
     g_assert_not_reached(); 
   }; /* end switch */

   /* Done? */
   if ((iflag>0) || (ialf==1) || (eta==0.0)) return psi;

   /* correction function */
   if (abs(eta) == 1.0) psi = 0.0;
   else psi = pow((1.0 - eta2), alpha[ialf-1]) * psi;

   return psi;
} /* end sphfn */

/** 
 * Convolves data in buffer on otfdata onto the grids
 * Arguments are given in the structure passed as arg
 * Can Run in Thread.
 * \param arg  Pointer to OTFGridFuncArg argument with elements
 * \li in     ObitOTFGrid object
 * \li OTFin   OTF data set to grid from current buffer 
 * \li first  First (1-rel) rec in OTFin buffer to process this thread
 * \li last   Highest (1-rel) rec in OTFin buffer to process this thread
 * \li ithread thread number, >0 -> no threading
 * \li grid   Data grid
 * \li wtgrid Weight grid
 * \li xpos, ypos, float arrays the size of ndetect
 */
static gpointer ThreadGridBuffer (gpointer arg)
{
  /* Get arguments from structure */
  OTFGridFuncArg *largs = (OTFGridFuncArg*)arg;
  ObitOTFGrid *in     = largs->in;
  ObitOTF *otfdata    = largs->OTFin;
  olong loRec         = largs->first-1;
  olong hiRec         = largs->last;
  ObitFArray *tgrid   = largs->grid;
  ObitFArray *tgridWt = largs->wtgrid;
  ofloat *xpos        = largs->xpos;
  ofloat *ypos        = largs->ypos;

  olong irec, nrec, ndet, idet, ncol, nrow, ix, iy, icx, icy;
  olong lGridRow, lGridCol, itemp, i;
  ofloat *grid, *gridWt, *gridStart, *gridTop, *gridWtStart, *gridWtTop;
  ofloat *convfnp, *x, *y, *rot, *rec, *cal, rtemp, clip, dataVal;
  ofloat xcell, ycell, wt, dataWt, fblank = ObitMagicF();
  olong  convx, convy;
  olong pos[] = {0,0,0,0,0};
  olong incdatawt;
  gboolean doDataWt;
  ObitOTFDesc *desc;
  ObitOTFArrayGeom *geom;

   /* how much data? */
  desc  = otfdata->myDesc;
  geom  = otfdata->geom;
  nrec  = hiRec - loRec + 1;;   /* number of data records */
  if (nrec<=0) goto finish; /* need something */

  /* number of detectors */
  ndet = 1;
  for (i=0; i<desc->naxis; i++) ndet *= MAX (1, desc->inaxes[i]);
  incdatawt = MAX (1, desc->incdatawt); /* increment in data-wt axis */
  ndet /= incdatawt;
  doDataWt = incdatawt>1;      /* Have Data-Wt axis? */
  clip = in->clip;             /* Clipping level */

  /* DEBUG
  doDataWt = FALSE; */

  /* initialize data pointers */
  x   = otfdata->buffer+desc->ilocra +loRec*desc->lrec;
  y   = otfdata->buffer+desc->ilocdec+loRec*desc->lrec;
  cal = otfdata->buffer+desc->iloccal+loRec*desc->lrec;
  rot = otfdata->buffer+desc->ilocrot+loRec*desc->lrec;
  rec = otfdata->buffer+desc->ilocdata+loRec*desc->lrec;
  
  lGridRow = tgrid->naxis[0];     /* length of grid row */
  lGridCol = tgrid->naxis[1];     /* length of grid column */

  /* beginning of the grids */
  pos[0] = 0;  pos[1] = 0;
  gridStart   = ObitFArrayIndex (tgrid, pos); 
  gridWtStart = ObitFArrayIndex (tgridWt, pos); 

  /* beginning of highest row */
  pos[1] = lGridCol-1;
  gridTop   = ObitFArrayIndex (tgrid, pos); 
  gridWtTop = ObitFArrayIndex (tgridWt, pos); 
  dataWt = 1.0;

  ncol = in->convWidth;
  nrow = in->convWidth;
  
  /* Loop over data records */
  for (irec=loRec; irec<hiRec; irec++) {

    /* Get detector positions projected onto image */
    ObitOTFArrayGeomProj (geom, *x, *y, *rot, in->raProj, in->decProj, in->Proj, 
			  xpos, ypos);
    /* loop over detectors */
    for (idet = 0; idet<ndet; idet++) {
      if ((rec[idet*incdatawt] == fblank) || (rec[idet*incdatawt+1]<=0.0)) continue;  /* flagged? */

      /* Get grid coordinates for this datum */
      xcell = xpos[idet] * in->XScale + in->icenxImage;
      ycell = ypos[idet] * in->YScale + in->icenyImage;

      /* get center cell */
      ix = (olong)(xcell + 0.5);
      iy = (olong)(ycell + 0.5);

      /* Bail out if too close to edge */
      if ((ix<in->convWidth) || (in->nxImage-ix<in->convWidth)) continue;
      if ((iy<in->convWidth) || (in->nyImage-iy<in->convWidth)) continue;

      /* back off half Kernel width */
      ix -= in->convWidth/2;
      iy -= in->convWidth/2;

      /* Starting convolution location, table has in->convNperCell points per cell */
      /* Determine fraction of the cell to get start location in convolving table. */
      if (xcell > 0.0) itemp = (olong)(xcell + 0.5);
      else itemp = ((olong)(xcell - 0.5));
      rtemp = in->convNperCell*(xcell - itemp + 0.5);
      if (rtemp > 0.0) rtemp += 0.5;
      else rtemp -= 0.5;
      convx = (olong)rtemp;
      
      /* now y convolving fn */
      if (ycell > 0.0) itemp = (olong)(ycell + 0.5);
      else itemp = ((olong)(ycell - 0.5));
      rtemp = in->convNperCell*(ycell - itemp + 0.5);
      if (rtemp > 0.0) rtemp += 0.5;
      else rtemp -= 0.5;
      convy = (olong)rtemp;
    
      /* Starting location in grids */
      pos[0] = ix;
      pos[1] = iy;
      grid   = ObitFArrayIndex (tgrid, pos); 
      gridWt = ObitFArrayIndex (tgridWt, pos);

      if (doDataWt) dataWt = rec[idet*incdatawt+1];
      else dataWt = 1.0;
      dataVal = rec[idet*incdatawt];
      if (fabs(dataVal)>clip) dataWt = 0.0;  /* Clip? */

      /* Do gridding */
      for (icy=0; icy<nrow; icy++) { /* Loop over rows */

	/* convolution fn pointer */
	pos[0] = convx; pos[1] = convy;
	convfnp = ObitFArrayIndex (in->convfn, pos);
	
	for (icx=0; icx<ncol; icx++) { /* Loop over columns */

	  wt = (*convfnp) * dataWt;
          if (dataWt<=0.0) wt = 0.0;
	  if (fabs(wt)>0.0) {
	    grid[icx]   += dataVal * wt;  /* sum data * wt */
	    /* gridWt[icx] += fabs(wt);       sum wt - not sure about fabs  */
	    gridWt[icx] += wt;         /* sum wt  */

	    /* DEBUG - look at values being put into a particular cell  
	    if (((ix+icx)==251) && ((iy+icy)==148)) {  
	      fprintf (stdout,"DEBUG %8d %2.2d %7.2f %7.2f d %6.2f w %8.2f dw %6.1f cf %7.3f g %9.2f gwt %9.2f f %6.2f\n", 
		       irec+(desc->firstRec), idet, 
		       xcell,  ycell, rec[idet*incdatawt], wt, dataWt, *convfnp, 
		       grid[icx], gridWt[icx], grid[icx]/gridWt[icx]);  
	    }  *//* end DEBUG */
	  }

	  convfnp  += in->convNperCell;
	} /* end loop over columns */

	grid   += lGridRow;
	gridWt += lGridRow;
	convy  += in->convNperCell;
      } /* end loop over rows */
	
    } /* end loop over detectors */

    /* update data pointers */
    x   += desc->lrec;
    y   += desc->lrec;
    rot += desc->lrec;
    cal += desc->lrec;
    rec += desc->lrec;
  } /* end loop over buffer */

  /* Indicate completion */
 finish:
  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadGridBuffer */

