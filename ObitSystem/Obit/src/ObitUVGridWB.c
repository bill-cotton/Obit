/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010-2014                                          */
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
#include "ObitUVGridWB.h"
#include "ObitThreadGrid.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVGridWB.c
 * ObitUVGridWB class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVGridWB";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitUVGridGetClass;

/** Degrees to radians factor */
#ifndef DG2RAD  
#define DG2RAD G_PI / 180.0
#endif

/**  Radians to degrees factor */
#ifndef RAD2DG  
#define RAD2DG 180.0 / G_PI
#endif

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitUVGridWBClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVGridWBClassInfo myClassInfo = {FALSE};

/*---------------Private structures----------------*/
/* FFT/gridding correction threaded function argument */
typedef struct {
  /* ObitThread with restart queue */
  ObitThread *thread;
  /* SkyModel with model components loaded (ObitSkyModelLoad) */
  ObitUVGrid *in;
  /* UV data set to model and subtract from current buffer */
  ObitFArray *array;
  /* If doGridCor=TRUE do gridding correction   */
  gboolean    doGridCor;
  /* thread number, >0 -> no threading   */
  olong       ithread;
} FFT2ImFuncArg;
/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVGridWBInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVGridWBClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVGridWBClassInfoDefFn (gpointer inClass);

/** Private: Prepare visibility data for gridding */
static void PrepBufferWB (ObitUVGrid* in, ObitUV *uvdata, 
			  olong loVis, olong hiVis);

/** Private: Threaded FFT/gridding correct */
static gpointer ThreadFFT2ImWB (gpointer arg);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVGridWB* newObitUVGridWB (gchar* name)
{
  ObitUVGridWB* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVGridWBClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVGridWB));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVGridWBInit((gpointer)out);

 return out;
} /* end newObitUVGridWB */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVGridWBGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVGridWBClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVGridWBGetClass */

/**
 * Prepares for gridding uv data of the type described by UVin and
 * with derived image as described by image.
 * Wideband (SW) version
 * Input data should be fully edited and calibrated, with any weighting applied 
 * and converted to the appropriate Stokes type.
 * The object UVin will be opened during this call if it is not already open.
 * image should describe the center, size and grid spacing of the desired
 * image.
 * The beams corresponding to each image should be made first using the
 * same ObitUVGridWB.
 * \param in       Object to initialize
 * \param UVin     Uv data object to be gridded.
 * \param imagee   Image (beam) to be gridded. (as Obit*)
 *                 Descriptor infoList entry "BeamTapr" gives any additional
 *                 tapering in degrees.
 * \param doBeam   TRUE is this is for a set of Beams.
 * \param err      ObitErr stack for reporting problems.
 */
void ObitUVGridWBSetup (ObitUVGrid *inn, ObitUV *UVin, Obit *imagee,
			gboolean doBeam, ObitErr *err)
{
  ObitIOCode retCode;
  ObitUVDesc *uvDesc;
  ObitUVGridWB *in = (ObitUVGridWB*)inn;
  ObitImageDesc *theDesc=NULL;
  ObitImage *image = (ObitImage*)imagee;
  ObitImage *myBeam;
  olong nx, ny, naxis[2];
  ofloat cellx, celly, dxyzc[3], xt, yt, zt, taper, BeamTaper=0.0;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean doCalSelect = FALSE;
  ObitIOAccess access;
  gchar *routine="ObitUVGridWBSetup";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVGridWBIsA(in));
  g_assert (ObitUVIsA(UVin));
  g_assert (ObitImageWBIsA(image));
  Obit_return_if_fail((image->myDesc->inaxes[0]>0) && 
		      (image->myDesc->inaxes[1]>0), err,
		      "%s: MUST fully define image descriptor %s",
		      routine, image->name);
  
  /* Need beam */
  myBeam = (ObitImage*)imagee;
  Obit_return_if_fail(ObitImageIsA(myBeam), err,
		      "%s: Beam for %s not defined", 
		      routine, image->name);
  
  /* Applying calibration or selection? */
  ObitInfoListGetTest(UVin->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;
  
  /* open uv data to fully instantiate if not already open */
  if (in->myStatus==OBIT_Inactive) {
    retCode = ObitUVOpen (UVin, access, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  uvDesc = UVin->myDesc;
  in->order   = ((ObitImageWB*)image)->order;    /* save spectral order */
  in->refFreq = ((ObitImageWB*)image)->refFreq;  /* reference frequency */

  /* Get source position if it's not already in header */
  if ((uvDesc->crval[uvDesc->jlocr]==0.0) && 
      (uvDesc->crval[uvDesc->jlocd]==0.0)) {
    ObitUVGetRADec (UVin, &uvDesc->crval[uvDesc->jlocr], 
			&uvDesc->crval[uvDesc->jlocd], err);
    if (err->error) Obit_traceback_msg (err, routine, UVin->name);
  }

  /* Beam, image dependent stuff */
  in->nxBeam = myBeam->myDesc->inaxes[0];
  in->nyBeam = myBeam->myDesc->inaxes[1];
  in->icenxBeam = in->nxBeam/2 + 1; 
  in->icenyBeam = in->nyBeam/2 + 1;
  in->nxImage = image->myDesc->inaxes[0];
  in->nyImage = image->myDesc->inaxes[1];
  in->icenxImage = in->nxImage/2 + 1;
  in->icenyImage = in->nyImage/2 + 1;
  
  /* Any additional tapering (deg) */
  ObitInfoListGetTest(image->myDesc->info, "BeamTapr", &type, dim, &BeamTaper);
  if (BeamTaper>0.0) {
    taper   = (1.0 / (((BeamTaper/2.35)/206265.))/(G_PI));
    in->BeamTaperUV = log(0.3)/(taper*taper);
  } else in->BeamTaperUV = 0.0;
  
  /* Get values by Beam/Image */
  in->doBeam = doBeam;
  if (doBeam) {
    theDesc = myBeam->myDesc;  /* Which descriptor in use */
    /* shift parameters */
    /* zeros for beam */
    in->dxc = 0.0;
    in->dyc = 0.0;
    in->dzc = 0.0;

  } else {
    /* shift parameters */
    theDesc = image->myDesc;  /* Which descriptor in use */
    ObitUVDescShiftPhase (uvDesc, image->myDesc, dxyzc, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    in->dxc = -dxyzc[0];
    in->dyc = -dxyzc[1];
    in->dzc = -dxyzc[2];
  }

  /* create/resize grid as needed */
  naxis[0] = 1 + theDesc->inaxes[0] / 2;
  naxis[1] = theDesc->inaxes[1];

  if (in->grid==NULL) {
    in->grid = ObitCArrayCreate ("UV Grid", 2, naxis);
    /* reallocate if need be, zero in any case */
  } else {
    in->grid = ObitCArrayRealloc (in->grid, 2, naxis);
  }
   /* Scaling to cells */
  nx = theDesc->inaxes[0];
  ny = theDesc->inaxes[1];
  cellx = (DG2RAD) * theDesc->cdelt[0]; /* x cells spacing in radians */
  celly = (DG2RAD) * theDesc->cdelt[1]; /* y cells spacing in radians */
  in->UScale =  nx * fabs(cellx);
  /* Flip sign on v to make maps come out upside down. */
  in->VScale = -ny * fabs(celly);
  in->WScale = 1.0;

  /* 3D rotation matrix */
  in->rotate = theDesc->crota[1] - uvDesc->crota[1]; /* rotation */
  in->do3Dmul = ObitUVDescShift3DMatrix (uvDesc, theDesc, in->URot3D, in->PRot3D);

  /* Rotate shift parameters if needed. */
  if (in->do3Dmul) {
    xt = (in->dxc)*in->PRot3D[0][0] + (in->dyc)*in->PRot3D[1][0] + (in->dzc)*in->PRot3D[2][0];
    yt = (in->dxc)*in->PRot3D[0][1] + (in->dyc)*in->PRot3D[1][1] + (in->dzc)*in->PRot3D[2][1];
    zt = (in->dxc)*in->PRot3D[0][2] + (in->dyc)*in->PRot3D[1][2] + (in->dzc)*in->PRot3D[2][2];
    /*fprintf (stderr,"scale %10.8f %10.8f %10.8f\n",in->UScale,in->VScale,in->WScale); */
    in->dxc = xt;
    in->dyc = yt;
    in->dzc = zt;
  }


  /* frequency tables if not defined */
  if ((uvDesc->freqArr==NULL) || (uvDesc->fscale==NULL)) {
    ObitUVGetFreq (UVin, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end setup frequency table */


 }  /* end ObitUVGridWBSetup */

/**
 * Perform half plane complex to real FFT, convert to center at the center order and
 * apply corrections for the convolution  function used in gridding
 * Requires setup by #ObitUVGridCreate and gridding by #ObitUVGridReadUV.
 * Writes image to disk.
 * \param in      Object to initialize
 *                info element "Channel" has plane number, def[1]
 * \param oout    Output image, should be open for call (as Obit*)
 * \param array   Output image array.
 * \param err     ObitErr stack for reporting problems.
 */
void ObitUVGridWBFFT2Im (ObitUVGrid *inn, Obit *oout, ObitErr *err)
{
  ObitUVGridWB *in = (ObitUVGridWB*)inn;
  ObitImage *out = (ObitImage*)oout;
  ofloat *ramp=NULL, *data=NULL, *imagep=NULL, *xCorrp=NULL, *yCorrp=NULL, fact;
  olong size, naxis[2], pos[5], pln, plane[5]={1,1,1,1,1};
  ObitFArray *xCorrTemp=NULL;
  ObitFArray *array = out->image;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong xdim[7];
  ObitUVGridClassInfo *gridClass = (ObitUVGridClassInfo*)in->ClassInfo; /* Gridder class */
  gchar *routine = "ObitUVGridWBFFT2Im";
  
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVGridIsA(in));
  g_assert (ObitFArrayIsA(array));
  /* Arrays compatable */
  g_assert (in->grid->ndim == array->ndim);
  g_assert (2*(in->grid->naxis[0]-1) == array->naxis[0]);
  g_assert (in->grid->naxis[1] == array->naxis[1]);
  
  /*ObitErrTimeLog(err, routine);  Add Timestamp */
  
  /* Beam or image? */
  if (in->doBeam) { 
    /* Making Beam */ 
    /* Create FFT object if not done before */
    if (in->FFTBeam==NULL) {
      xdim[0] = in->nxBeam; 
      xdim[1] = in->nyBeam; 
      in->FFTBeam = newObitFFT ("Beam FFT", OBIT_FFT_Reverse, OBIT_FFT_HalfComplex,
				2, xdim);
    }

    /* do FFT */
    ObitFFTC2R (in->FFTBeam, in->grid, array);
 
    /* reorder to center at center */
    ObitFArray2DCenter (array);

    /* Do gridding corrections */
    /* Create arrays / initialize if not done */
    if ((in->xCorrBeam==NULL) || (in->yCorrBeam==NULL)) {
      size = in->convWidth * in->convNperCell + 1;
      ramp = g_malloc0(2*size*sizeof(float));
      data = g_malloc0(2*size*sizeof(float));
      naxis[0] = in->nxBeam;
      in->xCorrBeam = ObitFArrayUnref(in->xCorrBeam); /* just in case */
      in->xCorrBeam = ObitFArrayCreate ("X Beam gridding correction", 1, naxis);
      naxis[0] = in->nyBeam;
      in->yCorrBeam = ObitFArrayUnref(in->yCorrBeam); /* just in case */
      in->yCorrBeam = ObitFArrayCreate ("Y Beam gridding correction", 1, naxis);
      
      /* X function */
      gridClass->GridCorrFn (inn, in->nxBeam, in->icenxBeam, data, ramp, in->xCorrBeam);
      
      /* If Y axis */
      gridClass->GridCorrFn (inn, in->nyBeam, in->icenyBeam, data, ramp, in->yCorrBeam);
    } /* end initialize correction functions */
    
    /* Normalization: use center value of beam - only order 0 */
    if (in->order==0) {
      pos[0] = in->icenxBeam-1; pos[1] = in->icenyBeam-1;
      imagep = ObitFArrayIndex(array, pos);
      xCorrp = ObitFArrayIndex(in->xCorrBeam, pos);
      pos[0] = in->icenyBeam-1;
      yCorrp = ObitFArrayIndex(in->yCorrBeam, pos);
      in->BeamNorm = (*imagep) * (*xCorrp) * (*yCorrp);
    }   
    /* MUST have beam peak for normalization */
    if (in->BeamNorm==0.0) {
      Obit_log_error(err, OBIT_Error, 
		     "ObitUVGridFFT2Im: MUST have made beam first: %s",
		     in->name);
      return;
    }
    
    /* Correct xCorr by normalization factor */
    fact = 1.0 / MAX (1.0e-20, in->BeamNorm);
    xCorrTemp = ObitFArrayCopy (in->xCorrBeam, NULL, err);
    ObitFArraySMul (xCorrTemp, fact);
    
    /* Do multiply */
    ObitFArrayMulColRow (array, xCorrTemp, in->yCorrBeam, array);
    
  } else { 
    /* Making Image */ 
    /* MUST have beam peak for normalization */
    if (in->BeamNorm==0.0) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: MUST have made beam first: %s",
		     routine, in->name);
      return;
    }
    
    /* Create FFT object if not done before */
    if (in->FFTImage==NULL) {
      xdim[0] = in->nxImage; 
      xdim[1] = in->nyImage; 
      in->FFTImage = newObitFFT ("Image FFT", OBIT_FFT_Reverse, OBIT_FFT_HalfComplex,
				 2, xdim);
    }
    
    /* do FFT */
    ObitFFTC2R (in->FFTImage, in->grid, array);
    
    /* reorder to cernter at center */
    ObitFArray2DCenter (array);
    
    /* Do gridding corrections */
    /* Create arrays / initialize if not done */
    if ((in->xCorrImage==NULL) || (in->yCorrImage==NULL)) {
      size = in->convWidth * in->convNperCell + 1;
      ramp = g_malloc0(2*size*sizeof(float));
      data = g_malloc0(2*size*sizeof(float));
      naxis[0] = in->nxImage;
      in->xCorrImage = ObitFArrayUnref(in->xCorrImage); /* just in case */
      in->xCorrImage = ObitFArrayCreate ("X Image gridding correction", 1, naxis);
      naxis[0] = in->nyImage;
      in->yCorrImage = ObitFArrayUnref(in->yCorrImage); /* just in case */
      in->yCorrImage = ObitFArrayCreate ("Y Image gridding correction", 1, naxis);
      
      /* X function */
      gridClass->GridCorrFn (inn, in->nxImage, in->icenxImage, data, ramp, in->xCorrImage);
      
      /* If Y axis */
      gridClass->GridCorrFn (inn, in->nyImage, in->icenyImage, data, ramp, in->yCorrImage);
    } /* end initialize correction functions */
    
    /* Normalization: use center value of beam */
    /* Correct xCorr by normalization factor */
    fact = 1.0 / MAX (1.0e-20, in->BeamNorm);
    xCorrTemp = ObitFArrayCopy (in->xCorrImage, NULL, err);
    ObitFArraySMul (xCorrTemp, fact);
    
    /* Do multiply   */
    ObitFArrayMulColRow (array, xCorrTemp, in->yCorrImage, array); /* DEBUG*/
    /* DEBUG try turning off gridding correction */

     /* DEBUG ObitFArraySMul (array, fact);  just normalize instead */
  } /* end make image */

  /* cleanup */
  xCorrTemp = ObitFArrayUnref(xCorrTemp);
  if (ramp) g_free (ramp); ramp = NULL;
  if (data) g_free (data); data = NULL;

  /* Write output */
  pln = 1;  /* Get channel/plane number */
  ObitInfoListGetTest(in->info, "Channel", &type, dim, &pln);
  plane[0] = pln;
  ObitImagePutPlane (out, array->array, plane, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
	
} /* end ObitUVGridWBFFT2Im */

 /**
 * Parallel perform half plane complex to real FFT, convert to center at the 
 * center order and apply corrections for the convolution  function used in gridding
 * Requires setup by #ObitUVGridWBCreate and gridding by #ObitUVGridWBReadUV.
 * If Beams are being made, there should be entries in in and array for both 
 * beam and image with the beam immediately prior to the associated image.
 * Apparently the threading in FFTW clashes with that in Obit so here the
 * FFTs are done sequentially 
 * Images written to disk
 * Wideband version:
 * \param nPar    Number of parallel griddings
 * \param in      Array of  objects to grid
 *                info element "Channel" has plane number, def[1]
 * \param oout    Array of output images,  pixel array elements must correspond 
 *                to those in in. (As Obit*)
 * \param err     ObitErr stack for reporting problems.
 */
void ObitUVGridWBFFT2ImPar (olong nPar, ObitUVGrid **inn, Obit **oout, ObitErr *err)
{
  ObitImage **out = (ObitImage**)oout;
  olong i, nTh, nnTh, off, nLeft, pos[5], xdim[7], pln, plane[5]={1,1,1,1,1};
  FFT2ImFuncArg *args=NULL;
  ObitFArray **array=NULL;
  ObitUVGridWB **in = (ObitUVGridWB**)inn;
  ObitThreadFunc func=(ObitThreadFunc)ThreadFFT2ImWB;
  gboolean OK;
  ofloat fact, *Corrp, BeamNorm=0.0;
  gchar *routine = "ObitUVGridWBFFT2ImPar";

  /* error checks */
  if (err->error) return;
  if (nPar<=0)    return;

  for (i=0; i<nPar; i++) {
    g_assert (ObitUVGridWBIsA(in[i]));
    g_assert (ObitImageIsA(out[i]));  
  }
  
  /* Create FArray array */
  array = g_malloc0(nPar*sizeof(ObitFArray*));

  /* FFTs */
  for (i=0; i<nPar; i++) {
    /* Create FFT object if not done before */
    if (in[i]->doBeam) { /* Beam? */
      if (in[i]->FFTBeam==NULL) {
	xdim[0] = in[i]->nxBeam; 
	xdim[1] = in[i]->nyBeam; 
	in[i]->FFTBeam = newObitFFT ("Beam FFT", OBIT_FFT_Reverse, OBIT_FFT_HalfComplex,
				     2, xdim);
	/* Reference or create beam array for higher order beam */
	if (in[i]->order==0)
	  array[i] = ObitFArrayRef(out[i]->image);
	else
	  array[i] = ObitFArrayCreate("Beam", 2, xdim);
      }
      
      /* do FFT */
      ObitFFTC2R (in[i]->FFTBeam, in[i]->grid, array[i]);
      
    } else { /* Image */
      /* Create FFT object if not done before */
      if (in[i]->FFTImage==NULL) {
	xdim[0] = in[i]->nxImage; 
	xdim[1] = in[i]->nyImage; 
	in[i]->FFTImage = newObitFFT ("Image FFT", OBIT_FFT_Reverse, OBIT_FFT_HalfComplex,
				      2, xdim);
	array[i] = ObitFArrayRef(out[i]->image);
      }
      
      /* do FFT */
      ObitFFTC2R (in[i]->FFTImage, in[i]->grid, array[i]);
    }
  } /* end loop doing FFTs */

  /* How many threads? */
  in[0]->nThreads = MAX (1, ObitThreadNumProc(in[0]->thread));
  in[0]->nThreads = MIN (nPar, in[0]->nThreads);

  /* Initialize threadArg array put all on in[0] */
  if (in[0]->threadArgs==NULL) {
    in[0]->threadArgs = g_malloc0(in[0]->nThreads*sizeof(FFT2ImFuncArg*));
    for (i=0; i<in[0]->nThreads; i++) 
      in[0]->threadArgs[i] = g_malloc0(sizeof(FFT2ImFuncArg)); 
  } 
  
  /* How many threads? */
  nTh = in[0]->nThreads;

  /* do jobs, doing nTh in parallel */
  off = 0;
  /* Set up thread arguments for first nTh girddings */
  for (i=0; i<nTh; i++) {
    args = (FFT2ImFuncArg*)in[0]->threadArgs[i];
    args->thread    = in[0]->thread;
    args->in        = inn[i+off];
    args->array     = array[i+off];
    /* args->doGridCor = !inn[i+off]->doBeam; No corrections for beams DEBUG */
    /* args->doGridCor = FALSE;   No corrections for anything DEBUG */
    args->doGridCor = TRUE;    /* Grid correct everything */
    if (nTh>1) args->ithread = i;
    else args->ithread = -1;
  }
  
  /* Do operation on buffer possibly with threads */
  OK = ObitThreadIterator (in[0]->thread, nTh, func, in[0]->threadArgs);
  
  /* Check for problems */
  if (!OK) {
    Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
    goto cleanup;
  }
  
  /* Loop over rest of images */
  nLeft = nPar - nTh;
  off   = nTh;
  while (nLeft>0) {
    nnTh = MIN (nTh, nLeft);  /* How many to do? */
    
    for (i=0; i<nnTh; i++) {
      args = (FFT2ImFuncArg*)in[0]->threadArgs[i];
      args->thread = in[0]->thread;
      args->in    = inn[i+off];
      args->array = array[i+off];
      /* args->doGridCor = !inn[i+off]->doBeam; No corrections for beams DEBUG */
      /* args->doGridCor = FALSE;   No corrections for anything DEBUG */
      args->doGridCor = TRUE;    /* Grid correct everything */
      if (nnTh>1) args->ithread = i;
      else args->ithread = -1;
    }
    
    /* Do operation on buffer possibly with threads */
    OK = ObitThreadIterator (in[0]->thread, nnTh, func, in[0]->threadArgs);
    
    /* Check for problems */
    if (!OK) {
      Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
      goto cleanup;
    }
    off   += nnTh;  /* update offset */
    nLeft -= nnTh;  /* update number left */
  } /* end loop over rest */
  
  /* Normalize - loop looking for an in entry with doBeam member set and order 0,
   the center peak is measured and used to normalize,  this peak is assumed
   to be the normalization for the subsequent image.
   if an image without corresponding beam is encountered, the BeamNorm
   member of it in[] is used to normalize */

  for (i=0; i<nPar; i++) {
    /* is this a 0 order beam? */
    if ((in[i]->doBeam) && (in[i]->order==0)) {
      pos[0] = in[i]->icenxBeam-1; pos[1] = in[i]->icenyBeam-1; pos[2] = 1;
      Corrp = ObitFArrayIndex(array[i], pos);
      BeamNorm = *Corrp;
      /* Check */
      if (BeamNorm==0.0) {
	Obit_log_error(err, OBIT_Error, "%s ERROR peak in beam is zero for: %s",
		       routine, in[i]->name);
	goto cleanup;
      }
      fact = 1.0 / MAX (1.0e-20, BeamNorm);
      ObitFArraySMul (array[i], fact);  /* Normalize beam */
      /* Save normalization on in[i,i+1] */
      if (!in[i+1]->doBeam) in[i+1]->BeamNorm = BeamNorm;
      in[i]->BeamNorm = BeamNorm;
      /* Write output */
      pln = 1+in[i]->order;  /* get channel/plane number */
      plane[0] = pln;
      ObitImagePutPlane (out[i], array[i]->array, plane, err);
      if (err->error) goto cleanup;
      i++;       /* Advance to next beam or image */
      in[i]->BeamNorm = BeamNorm;  /* Save */
    } /* end if beam */
    
   /*  Now image or higher order beam */
    if (in[i]->BeamNorm==0.0) in[i]->BeamNorm = BeamNorm;  /* Save */
    if (in[i]->BeamNorm==0.0) {
      Obit_log_error(err, OBIT_Error, "%s ERROR image normalization is zero for: %s",
		     routine, in[i]->name);
      goto cleanup;
    }
    fact = 1.0 / MAX (1.0e-20, in[i]->BeamNorm);
    ObitFArraySMul (array[i], fact);  /* Normalize image */
    
    /* Write output */
    if (!in[i]->doBeam) {  /* Image */
      pln = 1;  /* get channel/plane number */
      plane[0] = pln;
      ObitImagePutPlane (out[i], array[i]->array, plane, err);
    } else {  /* Beam */
      pln = 1+in[i]->order;  /* get channel/plane number */
      plane[0] = pln;
      ObitImagePutPlane (out[i], array[i]->array, plane, err);
    }
    if (err->error) goto cleanup;
  } /* end normalization loop */

  /*  cleanup */
 cleanup:
  if (array) {
    for (i=0; i<nPar; i++) {
      array[i] = ObitFArrayUnref( array[i]);
    }
    g_free(array);
  }
    if (err->error) Obit_traceback_msg (err, routine, out[0]->name);
} /* end ObitUVGridWBFFT2ImPar */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVGridWBClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVGridWBClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVGridWBClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVGridWBClassInfoDefFn (gpointer inClass)
{
  ObitUVGridWBClassInfo *theClass = (ObitUVGridWBClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit        = (ObitClassInitFP)ObitUVGridWBClassInit;
  theClass->ObitClassInfoDefFn   = (ObitClassInfoDefFnFP)ObitUVGridWBClassInfoDefFn;
  theClass->ObitGetClass         = (ObitGetClassFP)ObitUVGridWBGetClass;
  theClass->ObitUVGridSetup      = (ObitUVGridSetupFP)ObitUVGridWBSetup;
  theClass->ObitUVGridFFT2Im     = (ObitUVGridFFT2ImFP)ObitUVGridWBFFT2Im;
  theClass->ObitUVGridFFT2ImPar  = (ObitUVGridFFT2ImParFP)ObitUVGridWBFFT2ImPar;

} /* end ObitUVGridWBClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVGridWBInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVGridWB *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->order   = 0;
  in->refFreq = 0.0;
} /* end ObitUVGridWBInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVGridWB* cast to an Obit*.
 */
void ObitUVGridWBClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVGridWB *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVGridWBClear */

/** 
 * Reorders grid and do gridding correction.
 * NOTE: threading in FFTW apparently conflicts with Obit threads.
 * Arguments are given in the structure passed as arg
 * \param arg  Pointer to FFT2ImFuncArg argument with elements
 * \li thread    Thread with restart queue
 * \li in        Input ObitUVGrid object
 * \li array     Output ObitFArray Image array
 * \li doGridCor If TRUE do gridding correction
 * \li ithread   thread number, >0 -> no threading 
 */
static gpointer ThreadFFT2ImWB (gpointer arg)
{
  /* Get arguments from structure */
  FFT2ImFuncArg *largs = (FFT2ImFuncArg*)arg;
  ObitUVGrid *in     = largs->in;
  ObitFArray *array  = largs->array;
  gboolean doGridCor = largs->doGridCor;

  /* local */
  ofloat *ramp=NULL, *data=NULL;
  olong size, naxis[2];
  ObitUVGridClassInfo *gridClass = (ObitUVGridClassInfo*)in->ClassInfo; /* Gridder class */

  /* Beam or image? */
  if (in->doBeam) { 
    /* Making Beam */ 
    /* reorder to center at center */
    ObitFArray2DCenter (array);
    
    /* Do gridding corrections? */
    if (doGridCor) {
      /* Create arrays / initialize if not done */
      if ((in->xCorrBeam==NULL) || (in->yCorrBeam==NULL)) {
	size = in->convWidth * in->convNperCell + 1;
	ramp = g_malloc0(2*size*sizeof(float));
	data = g_malloc0(2*size*sizeof(float));
	naxis[0] = in->nxBeam;
	in->xCorrBeam = ObitFArrayUnref(in->xCorrBeam); /* just in case */
	in->xCorrBeam = ObitFArrayCreate ("X Beam gridding correction", 1, naxis);
	naxis[0] = in->nyBeam;
	in->yCorrBeam = ObitFArrayUnref(in->yCorrBeam); /* just in case */
	in->yCorrBeam = ObitFArrayCreate ("Y Beam gridding correction", 1, naxis);
	
	/* X function */
	gridClass->GridCorrFn (in, in->nxBeam, in->icenxBeam, data, ramp, in->xCorrBeam);
	
	/* If Y axis */
	gridClass->GridCorrFn (in, in->nyBeam, in->icenyBeam, data, ramp, in->yCorrBeam);
      } /* end initialize correction functions */
      
      /* Do multiply to make griding correction */
      ObitFArrayMulColRow (array, in->xCorrBeam, in->yCorrBeam, array);
    }
    
  } else { 
    /* Making Image */ 
    /* reorder to center at center */
    ObitFArray2DCenter (array);
    
    /* Do gridding corrections? */
    if (doGridCor) {
      /* Create arrays / initialize if not done */
      if ((in->xCorrImage==NULL) || (in->yCorrImage==NULL)) {
	size = in->convWidth * in->convNperCell + 1;
	ramp = g_malloc0(2*size*sizeof(float));
	data = g_malloc0(2*size*sizeof(float));
	naxis[0] = in->nxImage;
	in->xCorrImage = ObitFArrayUnref(in->xCorrImage); /* just in case */
	in->xCorrImage = ObitFArrayCreate ("X Image gridding correction", 1, naxis);
	naxis[0] = in->nyImage;
	in->yCorrImage = ObitFArrayUnref(in->yCorrImage); /* just in case */
	in->yCorrImage = ObitFArrayCreate ("Y Image gridding correction", 1, naxis);
	
	/* X function */
	gridClass->GridCorrFn (in, in->nxImage, in->icenxImage, data, ramp, in->xCorrImage);
	
	/* If Y axis */
	gridClass->GridCorrFn (in, in->nyImage, in->icenyImage, data, ramp, in->yCorrImage);
	
	/* Do multiply to make griding correction */
	ObitFArrayMulColRow (array, in->xCorrImage, in->yCorrImage, array);
      } /* end initialize correction functions */
    }
  } /* end make image */
  
  /* cleanup */
  if (ramp) g_free (ramp); ramp = NULL;
  if (data) g_free (data); data = NULL;

  /* Indicate completion */
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadFFT2ImWB */

