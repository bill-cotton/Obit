/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010-2011                                          */
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
#include <unistd.h>
#include "ObitUVGridMF.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVGridMF.c
 * ObitUVGridMF class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVGridMF";

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
 * ClassInfo structure ObitUVGridMFClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVGridMFClassInfo myClassInfo = {FALSE};

/*---------------Private structures----------------*/
/** Gridding threaded function argument */
typedef struct {
  /* ObitThread with restart queue */
  ObitThread *thread;
  /* SkyModel with model components loaded (ObitSkyModelLoad) */
  ObitUVGridMF *in;
  /* UV data set to model and subtract from current buffer */
  ObitUV       *UVin;
  /* Range of IF (0-rel) to process this thread  */
  olong        BIF, EIF;
  /* Range of channel (0-rel) to process this thread  */
  olong        BChan, EChan;
  /* thread number, >0 -> no threading   */
  olong        ithread;
  /* Gridding array for thread */
  ObitCArray  *grid;
  /* Array of u,v,w values */
  ObitFArray  *uvw;
  /* Number of floats in buffer   */
  olong        buffSize;
  /* I/O buffer to be coppied to UVin buffer */
  ofloat       *buffer;
} UVGridFuncArg;

/** gridding correction threaded function argument */
typedef struct {
  /* ObitThread with restart queue */
  ObitThread *thread;
  /* SkyModel with model components loaded (ObitSkyModelLoad) */
  ObitUVGrid *in;
  /* UV data set to model and subtract from current buffer */
  ObitFArray *array;
  /* thread number, >0 -> no threading   */
  olong       ithread;
} FFT2ImFuncArg;
/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVGridMFInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVGridMFClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVGridMFClassInfoDefFn (gpointer inClass);

/** Private: Grid a single image/Beam possibly with Threads */
static void GridOne (ObitUVGridMF* in, ObitUV *UVin, UVGridFuncArg **args, 
		     ObitThread *thread, ObitErr *err);

/** Private: Prepare visibility data for gridding */
static void PrepBufferMF (ObitUVGridMF* in, ObitUV *uvdata, olong BIF, olong EIF,
			  olong BChan, olong EChan, ObitFArray *uvw,
			  ObitCArray *accGrid);

/** Private: Grid visibility data */
static void GridBufferMF (ObitUVGridMF* in, ObitUV *uvdata, olong BIF, olong EIF,
			  olong BChan, olong EChan, ObitFArray *uvw, 
			  ObitCArray *accGrid);


/** Private: Threaded prep/grid buffer */
static gpointer ThreadUVGridMFBuffer (gpointer arg);

/** Private: Threaded FFT/gridding correct */
static gpointer ThreadFFT2ImMF (gpointer arg);

/** Private: Copy UVW portion of a data buffer */
static void copyUVW (UVGridFuncArg *arg);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVGridMF* newObitUVGridMF (gchar* name)
{
  ObitUVGridMF* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVGridMFClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVGridMF));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVGridMFInit((gpointer)out);

 return out;
} /* end newObitUVGridMF */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVGridMFGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVGridMFClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVGridMFGetClass */

/**
 * Prepares for gridding uv data of the type described by UVin and
 * with derived image as described by image.
 * Wideband version
 * Input data should be fully edited and calibrated, with any weighting applied 
 * and converted to the appropriate Stokes type.
 * The object UVin will be opened during this call if it is not already open.
 * image should describe the center, size and grid spacing of the desired
 * image.
 * The beams corresponding to each image should be made first using the
 * same ObitUVGridMF.
 * \param in       Object to initialize
 * \param UVin     Uv data object to be gridded.
 * \param imagee   Image (beam) to be gridded. (as Obit*)
 *                 Descriptor infoList entry "BeamTapr" gives any additional
 *                 tapering in degrees.
 * \param doBeam   TRUE is this is for a Beam.
 * \param err      ObitErr stack for reporting problems.
 */
void ObitUVGridMFSetup (ObitUVGrid *inn, ObitUV *UVin, Obit *imagee,
			gboolean doBeam, ObitErr *err)
{
  ObitIOCode retCode;
  ObitUVDesc *uvDesc;
  ObitUVGridMF *in = (ObitUVGridMF*)inn;
  ObitImageDesc *theDesc=NULL;
  ObitImageMF *image = (ObitImageMF*)imagee;
  ObitImage *myBeam;
  olong nx, ny, naxis[2], iSpec, size, nif, nfreq, nn;
  ofloat cellx, celly, dxyzc[3], xt, yt, zt, BeamTaper=0.0;
  ofloat *ramp=NULL, *data=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean doCalSelect = FALSE;
  ObitIOAccess access;
  ofloat Beam[3] = {0.0,0.0,0.0}, tarBeam[3], corrBeam[3];
  ofloat sigma2v, sigma2u, cpa, spa, taper;
  odouble tarFreq;
  gchar *routine="ObitUVGridMFSetup";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVGridMFIsA(in));
  g_assert (ObitUVIsA(UVin));
  g_assert (ObitImageMFIsA(image));
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
  in->nSpec    = ((ObitImageMF*)image)->nSpec;     /* save number of coarse channels */
  in->maxOrder = ((ObitImageMF*)image)->maxOrder;  /* save max imaging order */

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
    /* DEBUG taper   = (1.0 / (BeamTaper*DG2RAD/2.35)/(G_PI));
       in->BeamTaperUV = log(0.3)/(taper*taper);*/
    taper = BeamTaper*DG2RAD;
    in->BeamTaperUV = -taper*taper*2.15169;
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

  /* create/resize grids IF, channel ranges as needed */
  naxis[0] = 1 + theDesc->inaxes[0] / 2;
  naxis[1] = theDesc->inaxes[1];

  if (in->grids==NULL)     in->grids     = g_malloc0(in->nSpec*sizeof(ObitCArray*));
  if (in->BIFSpec==NULL)   in->BIFSpec   = g_malloc0(in->nSpec*sizeof(olong));
  if (in->EIFSpec==NULL)   in->EIFSpec   = g_malloc0(in->nSpec*sizeof(olong));
  if (in->BChanSpec==NULL) in->BChanSpec = g_malloc0(in->nSpec*sizeof(olong));
  if (in->EChanSpec==NULL) in->EChanSpec = g_malloc0(in->nSpec*sizeof(olong));
  if (in->BeamNorms==NULL) in->BeamNorms = g_malloc0(in->nSpec*sizeof(ofloat));
  for (iSpec=0; iSpec<in->nSpec; iSpec++) {
    if (in->grids[iSpec]==NULL) {
      in->grids[iSpec] = ObitCArrayCreate ("UV Grid", 2, naxis);
      /* reallocate if need be, zero in any case */
    } else {
      in->grids[iSpec] = ObitCArrayRealloc (in->grids[iSpec], 2, naxis);
    }
    in->BIFSpec[iSpec]   = image->BIFSpec[iSpec];
    in->EIFSpec[iSpec]   = image->EIFSpec[iSpec];
    in->BChanSpec[iSpec] = image->BChanSpec[iSpec];
    in->EChanSpec[iSpec] = image->EChanSpec[iSpec];
  } /* end loop over coarse channel */

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

  /* Gridding correction functions */
  if (doBeam) {
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
    GridCorrFn (inn, in->nxBeam, in->icenxBeam, data, ramp, in->xCorrBeam);
    
    /* If Y axis */
    GridCorrFn (inn, in->nyBeam, in->icenyBeam, data, ramp, in->yCorrBeam);
  } else { /* image */
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
    GridCorrFn (inn, in->nxImage, in->icenxImage, data, ramp, in->xCorrImage);
    
    /* If Y axis */
    GridCorrFn (inn, in->nyImage, in->icenyImage, data, ramp, in->yCorrImage);    
  }

  /* tapers for forcing beam size - include additional beam tapers here */
  if (uvDesc->jlocif>=0) nif = uvDesc->inaxes[uvDesc->jlocif];
  else  nif = 1;
  nfreq = uvDesc->inaxes[uvDesc->jlocf];
  nn = nif*nfreq;
  if (in->sigma1) g_free(in->sigma1);
  if (in->sigma2) g_free(in->sigma2);
  if (in->sigma3) g_free(in->sigma3);
  in->sigma1 = g_malloc0(nn*sizeof(ofloat));
  in->sigma2 = g_malloc0(nn*sizeof(ofloat));
  in->sigma3 = g_malloc0(nn*sizeof(ofloat));

  /* Get restoring beam */
  ObitInfoListGetTest(UVin->info, "Beam", &type, dim, Beam);
  tarFreq = UVin->myDesc->freqIF[0];  /* Frequency of lowest IF */

  for (iSpec=0; iSpec<nn; iSpec++) {
    /* Target beam size scaled to this frequency */
    tarBeam[0] = Beam[0] * tarFreq/uvDesc->freqArr[iSpec];
    tarBeam[1] = Beam[1] * tarFreq/uvDesc->freqArr[iSpec];
    tarBeam[2] = Beam[2];
    /* Correction beam including any additional */
    corrBeam[0] = sqrt (MAX(1.0e-10,Beam[0]*Beam[0]-tarBeam[0]*tarBeam[0]) + BeamTaper);
    corrBeam[1] = sqrt (MAX(1.0e-10,Beam[1]*Beam[1]-tarBeam[1]*tarBeam[1]) + BeamTaper);
    corrBeam[2] = Beam[2];
    /* 0.8 fudge factor 0.9 may be better */
    taper   = (0.8/ (((corrBeam[0]/2.35)/206265.))/(G_PI));
    sigma2u = log(0.3)/(taper*taper);
    taper   = (0.8/ (((corrBeam[1]/2.35)/206265.))/(G_PI));
    sigma2v = log(0.3)/(taper*taper);
    cpa     = cos(corrBeam[2]*DG2RAD);
    spa     = sin(corrBeam[2]*DG2RAD);
    in->sigma1[iSpec]  = (cpa*cpa*sigma2v + spa*spa*sigma2u);
    in->sigma2[iSpec]  = (spa*spa*sigma2v + cpa*cpa*sigma2u);
    in->sigma3[iSpec]  = 2.0*cpa*spa*(sigma2v - sigma2u);
  } /* End setting tapers */
 
  if (ramp) g_free (ramp); ramp = NULL;
  if (data) g_free (data); data = NULL;

}  /* end ObitUVGridMFSetup */

/**
 * Read a UV data object, applying any shift and accumulating to grid.
 * Buffering of data will use the buffers as defined on UVin 
 * ("nVisPIO" in info member).
 * The UVin object will be closed at the termination of this routine.
 * Requires setup by #ObitUVGridCreate.
 * The gridding information should have been stored in the ObitInfoList on in:
 * \li "Guardband" OBIT_float scalar = maximum fraction of U or v range allowed in grid.
 *             Default = 0.4.
 * \li "MaxBaseline" OBIT_float scalar = maximum baseline length in wavelengths.
 *             Default = 1.0e15.
 * \li "startChann" OBIT_long scalar = first channel (1-rel) in uv data to grid.
 *             Default = 1.
 * \li "numberChann" OBIT_long scalar = number of channels in uv data to grid.
 *             Default = all.
 * \param inn     Object to initialize
 * \param UVin    Uv data object to be gridded.
 *                Should be the same as passed to previous call to 
 *                #ObitUVGridSetup for input in.
 * \param err     ObitErr stack for reporting problems.
 */
void ObitUVGridMFReadUV (ObitUVGrid *inn, ObitUV *UVin, ObitErr *err)
{
  ObitUVGridMF *in = (ObitUVGridMF*)inn;
  ObitIOCode retCode = OBIT_IO_OK;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ofloat temp;
  olong   itemp, naxis[2];
  olong i;
  UVGridFuncArg *args=NULL;
  gboolean doCalSelect;
  gchar *routine="ObitUVGridMFReadUV";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVGridIsA(in));
  g_assert (ObitUVIsA(UVin));
  g_assert (ObitUVDescIsA(UVin->myDesc));
  g_assert (UVin->myDesc->fscale!=NULL); /* frequency scaling table */

   /* If more than one Stokes issue warning */
  if ((UVin->myDesc->jlocs>=0) && 
      (UVin->myDesc->inaxes[UVin->myDesc->jlocs]>1)) {
      Obit_log_error(err, OBIT_InfoWarn, 
		    "%s: More than one Stokes  ( %d) in data, ONLY USING FIRST", 
		     routine, UVin->myDesc->inaxes[UVin->myDesc->jlocs]);
  }

  /* get gridding information */
  /* guardband */
  temp = 0.4;
  /* temp = 0.1; debug */
  ObitInfoListGetTest(in->info, "Guardband", &type, dim, &temp);
  in->guardband = temp;
 
  /* baseline range */
  temp = 1.0e15;
  ObitInfoListGetTest(in->info, "MaxBaseline", &type, dim, &temp);
  in->blmax = temp;
  temp = 0.0;
  ObitInfoListGetTest(in->info, "MinBaseline", &type, dim, &temp);
  in->blmin = temp;

  /* Spectral channels to grid */
  itemp = 1;
  ObitInfoListGetTest(in->info, "startChann", &type, dim, &itemp);
  in->startChann = itemp;
  itemp = 0; /* all */
  ObitInfoListGetTest(in->info, "numberChann", &type, dim, &itemp);
  in->numberChann = itemp;

  /* Calibrating or selecting? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(UVin->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);

  /* UVin should have been opened in  ObitUVGridSetup */
  
  /* How many threads? threading over nSpec */
  in->nThreads = MAX (1, ObitThreadNumProc(in->thread));
  in->nThreads = MIN (in->nThreads, in->nSpec);

  /* Initialize threadArg array  */
  if (in->threadArgs==NULL) {
    in->threadArgs = g_malloc0(in->nThreads*sizeof(UVGridFuncArg*));
    for (i=0; i<in->nThreads; i++) 
      in->threadArgs[i] = g_malloc0(sizeof(UVGridFuncArg)); 
  } 
  
  /* Set up thread arguments */
  naxis[0] = 3; naxis[1] = UVin->mySel->nVisPIO;  /* Size of uvw buffer */
  for (i=0; i<in->nThreads; i++) {
    args = (UVGridFuncArg*)in->threadArgs[i];
    args->thread   = in->thread;
    args->in       = in;
    args->UVin     = UVin;
    args->grid     = in->grids[i];
    args->uvw      = ObitFArrayCreate("UVW", 2, naxis);
    args->buffSize = 0;
    args->buffer   = NULL;
  }
  /* end initialize */

  /* loop gridding data */
  while (retCode == OBIT_IO_OK) {
    
    /* read buffer */
    if (doCalSelect) retCode = ObitUVReadSelect (UVin, NULL, err);
    else retCode = ObitUVRead (UVin, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    GridOne(in, UVin, (UVGridFuncArg **)in->threadArgs, in->thread, err);
  } /* end loop reading/gridding data */

  /* Shut down any threading */
  ObitThreadPoolFree (in->thread);
  if (in->threadArgs) {
    for (i=0; i<in->nThreads; i++) {
      if (in->threadArgs[i]) {
	args = (UVGridFuncArg*)in->threadArgs[i];
	args->uvw  = ObitFArrayUnref(args->uvw);
	g_free(in->threadArgs[i]);
      }
    }
    g_free(in->threadArgs);
  }
  in->threadArgs = NULL;
  in->nThreads   = 0;

  /* Close data */
  retCode = ObitUVClose (UVin, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitUVGridMFReadUV  */

/**
 * Parallel read a UV data object, applying any shifts and accumulating to grids.
 * Buffering of data will use the buffers as defined on UVin 
 * ("nVisPIO" in info member).
 * The UVin object will be closed at the termination of this routine.
 * Requires setup by #ObitUVGridCreate.
 * The gridding information should have been stored in the ObitInfoList on in[0]:
 * \li "Guardband" OBIT_float scalar = maximum fraction of U or v range allowed in grid.
 *             Default = 0.4.
 * \li "MaxBaseline" OBIT_float scalar = maximum baseline length in wavelengths.
 *             Default = 1.0e15.
 * \li "startChann" OBIT_long scalar = first channel (1-rel) in uv data to grid.
 *             Default = 1.
 * \li "numberChann" OBIT_long scalar = number of channels in uv data to grid.
 *             Default = all.
 * \param nPar    Number of parallel griddings
 * \param inn     Array of  objects to grid
 *                Each should be initialized by ObitUVGridSetup
 *                To include beams, double nPar and set doBeam member 
 *                on one of each pair.
 * \param UVin    Array of UV data objects to be gridded.
 *                Should be the same as passed to previous call to 
 *                #ObitUVGridSetup for input in element.
 *                MUST all point to same data set with same selection
 *                but possible different calibration.
 *                All but [0] should be closed.
 * \param err     ObitErr stack for reporting problems.
 */
void ObitUVGridMFReadUVPar (olong nPar, ObitUVGrid **inn, ObitUV **UVin, ObitErr *err)
{
  ObitUVGridMF **in = (ObitUVGridMF**)inn;
  ObitIOCode retCode = OBIT_IO_OK;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ofloat temp;
  olong i, ip, itemp, naxis[2];
  olong nTh, nnTh, off, nCopy, nLeft, doCalib;
  UVGridFuncArg *args=NULL;
  gboolean doCalSelect;
  ObitUV **UVArr  = NULL;
  ofloat **buffers= NULL;
  gchar *routine="ObitUVGridMFReadUVPar";
  /* DEBUG 
  ObitFArray *dbgRArr=NULL, *dbgIArr=NULL;*/

  /* error checks */
  if (err->error) return;
  if (nPar<=0) return;
  for (ip=0; ip<nPar; ip++) {
    g_assert (ObitUVGridIsA(in[ip]));
    g_assert (ObitUVIsA(UVin[ip]));
    g_assert (ObitUVDescIsA(UVin[ip]->myDesc));
    g_assert (UVin[ip]->myDesc->fscale!=NULL); /* frequency scaling table */
  }

  /*  ObitErrTimeLog(err, routine);  Add Timestamp */

   /* If more than one Stokes issue warning */
  if ((UVin[0]->myDesc->jlocs>=0) && 
      (UVin[0]->myDesc->inaxes[UVin[0]->myDesc->jlocs]>1)) {
      Obit_log_error(err, OBIT_InfoWarn, 
		    "%s: More than one Stokes  ( %d) in data, ONLY USING FIRST", 
		     routine, UVin[0]->myDesc->inaxes[UVin[0]->myDesc->jlocs]);
  }

  /* get gridding information */
  /* guardband */
  temp = 0.4;
  /* temp = 0.1; debug */
  ObitInfoListGetTest(in[0]->info, "Guardband", &type, dim, &temp);
  for (ip=0; ip<nPar; ip++) in[ip]->guardband = temp;
 
  /* baseline range */
  temp = 1.0e15;
  ObitInfoListGetTest(in[0]->info, "MaxBaseline", &type, dim, &temp);
  for (ip=0; ip<nPar; ip++) in[ip]->blmax = temp;
  
  temp = 0.0;
  ObitInfoListGetTest(in[0]->info, "MinBaseline", &type, dim, &temp);
  for (ip=0; ip<nPar; ip++) in[ip]->blmin = temp;

  /* Spectral channels to grid */
  itemp = 1;
  ObitInfoListGetTest(in[0]->info, "startChann", &type, dim, &itemp);
  for (ip=0; ip<nPar; ip++) in[ip]->startChann = itemp;
  itemp = 0; /* all */
  ObitInfoListGetTest(in[0]->info, "numberChann", &type, dim, &itemp);
  for (ip=0; ip<nPar; ip++) in[ip]->numberChann = itemp;

  /* Calibrating and/or selecting? */
  doCalib = 0;
  ObitInfoListGetTest(UVin[0]->info, "doCalib", &type, dim, &doCalib);
  doCalSelect = FALSE;
  ObitInfoListGetTest(UVin[0]->info, "doCalSelect", &type, dim, &doCalSelect);

  /* UVin[0] should have been opened in  ObitUVGridSetup */
  
  /* How many threads? */
  in[0]->nThreads = MAX (1, ObitThreadNumProc(in[0]->thread));
  in[0]->nThreads = MIN (in[0]->nSpec, in[0]->nThreads);

  /* Initialize threadArg array put all on in[0] */
  if (in[0]->threadArgs==NULL) {
    in[0]->threadArgs = g_malloc0(in[0]->nThreads*sizeof(UVGridFuncArg*));
    for (i=0; i<in[0]->nThreads; i++) 
      in[0]->threadArgs[i] = g_malloc0(sizeof(UVGridFuncArg)); 
  } 
  
  /* Set up thread arguments */
  naxis[0] = 3; naxis[1] = UVin[0]->mySel->nVisPIO;  /* Size of uvw buffer */
  for (i=0; i<in[0]->nThreads; i++) {
    args = (UVGridFuncArg*)in[0]->threadArgs[i];
    args->thread   = in[0]->thread;
    args->uvw      = ObitFArrayCreate("UVW", 2, naxis);
    args->buffSize = 0;
    args->buffer   = NULL;
  }

  /* How many threads? */
  nTh = in[0]->nThreads;

  /* Array for UV data */
  UVArr = g_malloc0(2*sizeof(ObitUV*));

  /* Buffer array */
  buffers    = g_malloc0(2*sizeof(ofloat*));
  buffers[0] = g_malloc0(UVin[0]->bufferSize*sizeof(ofloat));
  buffers[1] = g_malloc0(UVin[0]->bufferSize*sizeof(ofloat));

  /* delete buffer on UVin[0] */
  g_free(UVin[0]->buffer); UVin[0]->buffer=NULL;
  UVin[0]->bufferSize = 0;

  ObitErrLog(err);

  /* loop gridding data */
  while (retCode == OBIT_IO_OK) {

    /* Initial buffer load */
    off = 0;
    UVArr[0] = UVin[0]; /* used for master buffer */
    UVArr[1] = UVin[1]; /* UVArr[1] is the one to actually be used */

    /* read buffer - first used as master - do copy in GridOne */
    nCopy = 1;
    if (doCalSelect) 
      retCode = ObitUVReadMultiSelect (nCopy, UVArr, buffers, err);
    else 
      retCode = ObitUVReadMulti (nCopy, UVArr, buffers, err);
    if (retCode==OBIT_IO_EOF) break;  /* Finished */
    if (err->error) goto cleanup;
    
    /* Set up thread arguments for first Beam/Image */
    for (i=0; i<nTh; i++) {
      args = (UVGridFuncArg*)in[0]->threadArgs[i];
      if (doCalib<=0) { /* Copy buffer info? */
	UVArr[1]->myDesc->numVisBuff =  UVArr[0]->myDesc->numVisBuff;
	args->buffer   = buffers[0];
	args->buffSize = UVArr[1]->myDesc->numVisBuff*UVArr[1]->myDesc->lrec;
      }
    }
    if (err->error) goto cleanup;

    /* Do operation on buffer possibly with threads to grid first image/beam */
    in[off]->nThreads = in[0]->nThreads;
    UVArr[1]->buffer  = buffers[1];
    GridOne(in[off], UVArr[1], (UVGridFuncArg **)in[0]->threadArgs, in[0]->thread, err);
    if (err->error) goto cleanup;

    /* reset buffers */
    UVin[0]->buffer  = buffers[0];     /* master buffer */
    UVArr[1]->buffer = NULL;

    /* Loop over rest of griddings */
    nLeft = nPar - 1;
    off   = 1;
    while (nLeft>0) {
      nnTh = MIN (1, nLeft);  /* How many to do? */
 
      /* Set up thread arguments for next griddings
	 UVArr[0] and buffers[0] not used for gridding */
      UVArr[1]   = UVin[off];

      /* reload buffers - first used as master - copy buffer in thread unless doCalib */
      nCopy = 1;
      if (doCalSelect) 
	retCode = ObitUVReReadMultiSelect (nCopy, UVArr, buffers, err);
      else 
	retCode = ObitUVReReadMulti (nCopy, UVArr, buffers, err);
      if (err->error) goto cleanup;
      
      /* Set up thread arguments for next Beam/Image */
      for (i=0; i<nTh; i++) {
	args = (UVGridFuncArg*)in[0]->threadArgs[i];
	if (doCalib<=0) { /* Copy buffer info? */
	  UVArr[1]->myDesc->numVisBuff =  UVArr[0]->myDesc->numVisBuff;
	  args->buffer = buffers[0];
	  args->buffSize = UVArr[1]->myDesc->numVisBuff*UVArr[1]->myDesc->lrec;
	}
      }
      if (err->error) goto cleanup;
      
      /* Do operation on buffer possibly with threads to grid next image/beam */
      in[off]->nThreads = in[0]->nThreads;
      UVArr[1]->buffer  = buffers[1];
      GridOne(in[off], UVArr[1], (UVGridFuncArg **)in[0]->threadArgs, in[0]->thread, err);
      if (err->error) goto cleanup;
      
      /* reset buffers */
      UVin[0]->buffer  = buffers[0]; /* master buffer */
      UVArr[1]->buffer = NULL;

      off   += 1;  /* update offset */
      nLeft -= 1;  /* update number left */
   } /* end loop over others */
  } /* end loop reading/gridding data */

  /*ObitErrTimeLog(err, "Stop Grid Loop");  DEBUG */
  ObitErrLog(err);

  /* DEBUG - look at first complex grid
  dbgRArr = ObitCArrayMakeF(in[0]->grids[0]);
  dbgIArr = ObitCArrayMakeF(in[0]->grids[0]);
  ObitCArrayReal (in[0]->grids[0], dbgRArr);
  ObitCArrayImag (in[0]->grids[0], dbgIArr);
  ObitImageUtilArray2Image ("DbugGridReal0.fits", 0, dbgRArr, err);  
  ObitImageUtilArray2Image ("DbugGridImag0.fits", 0, dbgIArr, err);  
  dbgRArr = ObitFArrayUnref(dbgRArr);
  dbgIArr = ObitFArrayUnref(dbgIArr); */
  /* end DEBUG */
		    
 /* Cleanup */
 cleanup:

  /* Shut down any threading */
  ObitThreadPoolFree (in[0]->thread);
  if (in[0]->threadArgs) {
    for (i=0; i<in[0]->nThreads; i++) {
      args = (UVGridFuncArg*)in[0]->threadArgs[i];
      args->uvw  = ObitFArrayUnref(args->uvw);
     if (in[0]->threadArgs[i]) g_free(in[0]->threadArgs[i]);
    }
    g_free(in[0]->threadArgs);
  }
  in[0]->threadArgs = NULL;
  in[0]->nThreads   = 0;

  /* Remove buffers from UV objects - they are deleted below */
  for (ip=0; ip<nPar; ip++) {
    UVin[ip]->buffer     = NULL;
    UVin[ip]->bufferSize = 0;
  }

  /* Close data */
  retCode = ObitUVClose (UVin[0], err);

 if (buffers) {
   g_free(buffers[0]);
   g_free(buffers[1]);
   g_free(buffers);
  }
  if (UVArr) g_free(UVArr);

  if (err->error) Obit_traceback_msg (err, routine, in[0]->name);

} /* end ObitUVGridMFReadUVPar  */

/**
 * Perform half plane complex to real FFT, convert to center at the center order and
 * apply corrections for the convolution  function used in gridding
 * Requires setup by #ObitUVGridCreate and gridding by #ObitUVGridReadUV.
 * Writes image to disk.
 * \param in      Object to initialize
 *                info element "Channel" has plane number, def[1]
 * \param oout    Output image as Obit*
 *                Uses info element "BeamNorms" on beam for channel normalization
 * \param err     ObitErr stack for reporting problems.
 */
void ObitUVGridMFFFT2Im (ObitUVGrid *inn, Obit *oout, ObitErr *err)
{
  ObitImage *out = (ObitImage*)oout;
  ofloat *imagep=NULL, *xCorrp=NULL, *yCorrp=NULL, fact, fblank = ObitMagicF();
  ObitUVGridMF *in = (ObitUVGridMF*)inn;
  olong j,pos[5], pln, plane[5]={1,1,1,1,1};
  ObitFArray *xCorrTemp=NULL;
  olong xdim[7];
  ObitFArray *array=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitImage *theBeam;
  ObitInfoType type;
  ObitImageClassInfo *imgClass;
  gchar *routine = "ObitUVGridMFFFT2Im";
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitUVGridIsA(in));
  g_assert (ObitImageIsA(out)); 
  
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

    /* Make image plane array */
    array = ObitFArrayCreate("Beam", 2, xdim);

    /* Loop over spectral planes */
    for (j=0; j<in->nSpec; j++) {
      /* do FFT */
      ObitFFTC2R (in->FFTBeam, in->grids[j], array);
      /* reorder to center at center */
      ObitFArray2DCenter (array);
 

      /* Do gridding corrections */
      /* Normalization: use center value of beam */
      pos[0] = in->icenxBeam-1; pos[1] = in->icenyBeam-1;
      imagep = ObitFArrayIndex(array, pos);
      xCorrp = ObitFArrayIndex(in->xCorrBeam, pos);
      pos[0] = in->icenyBeam-1;
      yCorrp = ObitFArrayIndex(in->yCorrBeam, pos);
      in->BeamNorms[j] = (*imagep) * (*xCorrp) * (*yCorrp);

      /* MUST have beam peak for normalization */
      if (in->BeamNorms[j]==0.0) {
	Obit_log_error(err, OBIT_Error, 
		       "ObitUVGridFFT2Im: MUST have made beam first: %s",
		       in->name);
	goto cleanup;
      }
      
      /* Bad plane? */
      if (in->BeamNorms[j]==fblank) {
	ObitFArrayFill (array, fblank); /* Blank fill */
      } else { /* OK */
	/* Correct xCorr by normalization factor */
	fact = 1.0 / MAX (1.0e-20, in->BeamNorms[j]);
	xCorrTemp = ObitFArrayCopy (in->xCorrBeam, xCorrTemp, err);
	ObitFArraySMul (xCorrTemp, fact);
	
	/* Do multiply */
	ObitFArrayMulColRow (array, xCorrTemp, in->yCorrBeam, array);
      }
      
      /* Write output */
      pln = 2+in->maxOrder+j;  /* Get channel/plane number */
      if (in->nSpec==1) pln = 1;       /* Normal imaging */
      plane[0] = pln;
      ObitImagePutPlane (out, array->array, plane, err);
      if (err->error) goto cleanup;
      
    } /* end loop over spectral planes */

    /* Determine combined beam */
    if (in->nSpec>1) ObitImageMFCombine ((ObitImageMF*)out, FALSE, err);
    if (err->error) goto cleanup;

    /* Save BeamNorms array on Beam */
    dim[0] = in->nSpec;  dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut(out->info, "BeamNorms", OBIT_float, dim, in->BeamNorms);
    
  } else {
    /* Making Image */ 
    /* Make image plane array */
    xdim[0] = in->nxImage; 
    xdim[1] = in->nyImage; 
    array = ObitFArrayCreate("Beam", 2, xdim);

    /* Fetch BeamNorms from beam */
    imgClass  = (ObitImageClassInfo*)out->ClassInfo;    /* Image class */
    theBeam   = imgClass->ObitImageGetBeam(out, 0, plane, err);
    ObitInfoListGetTest(theBeam->info, "BeamNorms", &type, dim, in->BeamNorms);
    if (err->error) goto cleanup;
 
 
    /* Loop over spectral planes */
    for (j=0; j<in->nSpec; j++) {
      /* MUST have beam peak for normalization */
      if (in->BeamNorms[j]==0.0) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: MUST have made beam first: %s",
		       routine, in->name);
	goto cleanup;
      }
      
      /* Bad plane? */
      if (in->BeamNorms[j]==fblank) {
	ObitFArrayFill (array, fblank); /* Blank fill */
      } else { /* OK */
	/* Create FFT object if not done before */
	if (in->FFTImage==NULL) {
	  xdim[0] = in->nxImage; 
	  xdim[1] = in->nyImage; 
	  in->FFTImage = newObitFFT ("Image FFT", OBIT_FFT_Reverse, OBIT_FFT_HalfComplex,
				     2, xdim);
	}
	
	ObitFFTC2R (in->FFTImage, in->grids[j], array);
	
	/* reorder to center at center */
	ObitFArray2DCenter (array);
	
	/* Do gridding corrections */
	/* Create arrays / initialize if not done */
	
	/* Normalization: use center value of beam */
	/* Correct xCorr by normalization factor */
	fact = 1.0 / MAX (1.0e-20, in->BeamNorms[j]);
	xCorrTemp = ObitFArrayCopy (in->xCorrImage, xCorrTemp, err);
	ObitFArraySMul (xCorrTemp, fact);
	
	/* Do multiply   */
	ObitFArrayMulColRow (array, xCorrTemp, in->yCorrImage, array); 
      } /* end if valid plane */

      /* Write output */
      pln = 2+in->maxOrder+j;  /* Get channel/plane number */
      if (in->nSpec==1) pln = 1;  /* Normal imaging */
      plane[0] = pln;
      ObitImagePutPlane (out, array->array, plane, err);
      if (err->error) goto cleanup;

    } /* end loop over spectral planes */

    /* Determine combined image */
    if (in->nSpec>1) ObitImageMFCombine ((ObitImageMF*)out, TRUE, err);
    if (err->error) goto cleanup;

  } /* end make image */
  
  /*  cleanup */
 cleanup:
  xCorrTemp = ObitFArrayUnref(xCorrTemp);
  array = ObitFArrayUnref(array);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
} /* end ObitUVGridMFFFT2Im */

 /**
 * Parallel perform half plane complex to real FFT, convert to center at the 
 * center order and apply corrections for the convolution  function used in gridding
 * Requires setup by #ObitUVGridMFCreate and gridding by #ObitUVGridMFReadUV.
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
 *                to those in in.  as Obit*
 *                Uses info element "BeamNorms" on beams for channel normalization
 * \param err     ObitErr stack for reporting problems.
 */
void ObitUVGridMFFFT2ImPar (olong nPar, ObitUVGrid **inn, Obit **oout, ObitErr *err)
{
  ObitImage **out = (ObitImage**)oout;
  olong i, ii, j, ip, nTh, nnTh, off, nLeft, pos[5], xdim[7], pln;
  olong narr, plane[5]={1,1,1,1,1};
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  FFT2ImFuncArg *args=NULL;
  ObitUVGridMF **in = (ObitUVGridMF**)inn;
  ObitThreadFunc func=(ObitThreadFunc)ThreadFFT2ImMF;
  ObitFArray **array=NULL;
  ObitImageClassInfo *imgClass;
  ObitImage *theBeam;
  gboolean OK;
  ofloat fact, *Corrp, fblank = ObitMagicF();
  gchar *routine = "ObitUVGridMFFFT2ImPar";
  /* DEBUG 
     ObitFArray *dbgRArr=NULL, *dbgIArr=NULL;*/

  /* error checks */
  if (err->error) return;
  if (nPar<=0)    return;

  for (i=0; i<nPar; i++) {
    g_assert (ObitUVGridMFIsA(in[i]));
    g_assert (ObitImageIsA(out[i]));  
  }

  /* Create FArray array */
  narr = nPar*in[0]->nSpec;
  array = g_malloc0(narr*sizeof(ObitFArray*));

  /* FFTs, image arrays  */
  ip = 0;
  for (i=0; i<nPar; i++) {
    /* Create FFT object if not done before */
    if (in[i]->doBeam) { /* Beam? */
      if (in[i]->FFTBeam==NULL) {
	xdim[0] = in[i]->nxBeam; 
	xdim[1] = in[i]->nyBeam; 
	in[i]->FFTBeam = newObitFFT ("Beam FFT", OBIT_FFT_Reverse, OBIT_FFT_HalfComplex,
				     2, xdim);
      }
      
      /* DEBUG - look at first complex grid 
      if (i==0) {
	dbgRArr = ObitCArrayMakeF(in[0]->grids[0]);
	dbgIArr = ObitCArrayMakeF(in[0]->grids[0]);
	ObitCArrayReal (in[0]->grids[0], dbgRArr);
	ObitCArrayImag (in[0]->grids[0], dbgIArr);
	ObitImageUtilArray2Image ("DbugGridReal1.fits", 0, dbgRArr, err);  
	ObitImageUtilArray2Image ("DbugGridImag1.fits", 0, dbgIArr, err);  
	dbgRArr = ObitFArrayUnref(dbgRArr);
	dbgIArr = ObitFArrayUnref(dbgIArr); 
      }   end DEBUG */
      /* Create array of pixel arrays then FFT */
      for (j=0; j<in[i]->nSpec; j++) {
	array[ip] = ObitFArrayCreate("Beam", 2, xdim);
	/* do FFT */
	ObitFFTC2R (in[i]->FFTBeam, in[i]->grids[j], array[ip++]);
      }
      /* DEBUG - look at first complex grid 
      if (i==0) {
	dbgRArr = ObitCArrayMakeF(in[0]->grids[0]);
	dbgIArr = ObitCArrayMakeF(in[0]->grids[0]);
	ObitCArrayReal (in[0]->grids[0], dbgRArr);
	ObitCArrayImag (in[0]->grids[0], dbgIArr);
	ObitImageUtilArray2Image ("DbugGridReal2.fits", 0, dbgRArr, err);  
	ObitImageUtilArray2Image ("DbugGridImag2.fits", 0, dbgIArr, err);  
	ObitImageUtilArray2Image ("DbugRawBeam.fits",  0, array[0], err);  
	dbgRArr = ObitFArrayUnref(dbgRArr);
	dbgIArr = ObitFArrayUnref(dbgIArr);
      }  end DEBUG */

      
    } else { /* Image */
      /* Create FFT object if not done before */
      if (in[i]->FFTImage==NULL) {
	xdim[0] = in[i]->nxImage; 
	xdim[1] = in[i]->nyImage; 
	in[i]->FFTImage = newObitFFT ("Image FFT", OBIT_FFT_Reverse, OBIT_FFT_HalfComplex,
				      2, xdim);
      }
      /* Create array of pixel arrays */
      for (j=0; j<in[i]->nSpec; j++) {
	array[ip] = ObitFArrayCreate("Image", 2, xdim);
	/* do FFT */
	ObitFFTC2R (in[i]->FFTImage, in[i]->grids[j], array[ip++]);
      }
    }
    
  } /* end loop doing FFTs */

  /* DEBUG 
  ObitImageUtilArray2Image ("DbugRawBeam0.fits",  0, array[0], err);  */
  /* END DEBUG */

  /*  Do gridding corrections threaded */
  /* How many threads? */
  in[0]->nThreads = MAX (1, ObitThreadNumProc(in[0]->thread));
  in[0]->nThreads = MIN (nPar*in[0]->nSpec, in[0]->nThreads);

  /* Initialize threadArg array put all on in[0] */
  if (in[0]->threadArgs==NULL) {
    in[0]->threadArgs = g_malloc0(in[0]->nThreads*sizeof(FFT2ImFuncArg*));
    for (i=0; i<in[0]->nThreads; i++) 
      in[0]->threadArgs[i] = g_malloc0(sizeof(FFT2ImFuncArg)); 
  } 
  
  /* How many threads? */
  nTh = in[0]->nThreads;

  /* Gridding corrections/normalization - do jobs, doing nTh in parallel */

  /* Loop over images - do nSpec planes in parallel using threading */
  for (i=0; i<nPar; i++) {
    nLeft = in[0]->nSpec;
    off   = i*nLeft;

    /* Loop over nSpec planes */
    while (nLeft>0) {
      nnTh = MIN (nTh, nLeft);  /* How many to do? */
      
      for (ii=0; ii<nnTh; ii++) {
	args = (FFT2ImFuncArg*)in[0]->threadArgs[ii];
	args->thread = in[0]->thread;
	args->in     = inn[i];
	args->array  = array[ii+off];
	if (nnTh>1) args->ithread = ii;
	else args->ithread = -1;
      }
      
      /* Do operation on buffer possibly with threads */
      OK = ObitThreadIterator (in[0]->thread, nnTh, func, in[0]->threadArgs);
      
      /* Check for problems */
      if (!OK) {
	Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
	goto cleanup;
      }
      nLeft -= nnTh;  /* update number left */
      off += nnTh;
    } /* end loop over rest */
  } /* End loop gridding correcting images */

  /* DEBUG
  ObitImageUtilArray2Image ("DbugRawBeam1.fits",  0, array[0], err);   */
  /* END DEBUG */

  /* Normalize - loop looking for an in entry with doBeam member set an,
   the center peak is measured and used to normalize,  this peak is assumed
   to be the normalization for the subsequent image.
   if an image without corresponding beam is encountered, the BeamNorm
   member of it in[] is used to normalize */


  for (i=0; i<nPar; i++) {
    out[i]->image = ObitFArrayUnref(out[i]->image);  /* Free buffer */
    /* is this a beam? */
    if (in[i]->doBeam) {
      pos[0] = in[i]->icenxBeam-1; pos[1] = in[i]->icenyBeam-1; pos[2] = 1;
      /* Loop over planes */
      for (j=0; j<in[i]->nSpec; j++) {
	Corrp = ObitFArrayIndex(array[i*in[i]->nSpec+j], pos);
	in[i]->BeamNorms[j] = *Corrp;
	/* Check */
	if (in[i]->BeamNorms[j]==0.0) {  /* No data? */
	  Obit_log_error(err, OBIT_InfoWarn, 
			 "%s Peak in beam  %d is zero for:%s, blank plane",
			 routine, j, in[i]->name);
	  in[i]->BeamNorms[j] = fblank;
	  /* DEBUG
	  ObitImageUtilArray2Image ("DbugRawBeamBad.fits",  0, array[i*in[i]->nSpec+j], err);   */
	  /* END DEBUG */
	  ObitFArrayFill (array[i*in[i]->nSpec+j], fblank); /* Blank fill */

	} else {  /* OK */
	  
	  /* Normalize */
	  fact = 1.0 / MAX (1.0e-20, in[i]->BeamNorms[j]);
	  ObitFArraySMul (array[i*in[i]->nSpec+j], fact);  /* Normalize beam */
	}
	
	/* Save normalization on in[i+1] */
	if (!in[i+1]->doBeam) in[i+1]->BeamNorms[j] = in[i]->BeamNorms[j];
	/* Save BeamNorms array on Beam */
	dim[0] = in[i]->nSpec;  dim[1] = dim[2] = dim[3] = dim[4] = 1;
	ObitInfoListAlwaysPut(out[i]->info, "BeamNorms", OBIT_float, dim, in[i]->BeamNorms);
	
	pln = 2+in[i]->maxOrder+j;      /* Get channel/plane number */
	if (in[i]->nSpec==1) pln = 1;   /* Normal imaging */
	plane[0] = pln;
	ObitImagePutPlane (out[i], array[i*in[i]->nSpec+j]->array, plane, err);
	if (err->error) goto cleanup;
	out[i]->image = ObitFArrayUnref(out[i]->image);  /* Free buffer */
      } /* end loop over planes */

      /* Determine combined beam */
      if (in[i]->nSpec>1) ObitImageMFCombine ((ObitImageMF*)out[i], FALSE, err);
      if (err->error) goto cleanup;

      i++;       /* Advance to image or next beam */
    } /* end if beam */
    /*  Now image */

    /* Fetch BeamNorms from beam */
    imgClass  = (ObitImageClassInfo*)out[i]->ClassInfo;    /* Image class */
    theBeam   = imgClass->ObitImageGetBeam(out[i], 0, plane, err);
    ObitInfoListGetTest(theBeam->info, "BeamNorms", &type, dim, in[i]->BeamNorms);
    if (err->error) goto cleanup;
    
    /* Loop over planes */
    for (j=0; j<in[i]->nSpec; j++) {
      if (in[i]->BeamNorms[j]==0.0) {
	Obit_log_error(err, OBIT_Error,
		       "%s ERROR image normalization is zero for: %s",
		       routine, in[i]->name);
	goto cleanup;
      } else if (in[i]->BeamNorms[j]==fblank) { /* Plane blanked? */

	ObitFArrayFill (array[i*in[i]->nSpec+j], fblank); /* Blank fill */
	
      } else {

	/* Normalize */
	fact = 1.0 / MAX (1.0e-20, in[i]->BeamNorms[j]);
	ObitFArraySMul (array[i*in[i]->nSpec+j], fact);  /* Normalize image */
      }
    
      /* Write output */
      pln = 2+j+in[i]->maxOrder;       /* Get channel/plane number */
      if (in[i]->nSpec==1) pln = 1;    /* Normal imaging */
      plane[0] = pln;
      ObitImagePutPlane (out[i], array[i*in[i]->nSpec+j]->array, plane, err);
      if (err->error) goto cleanup;
    } /* end loop over planes */

    /* Determine combined image */
    if (in[i]->nSpec>1) ObitImageMFCombine ((ObitImageMF*)out[i], TRUE, err);
    if (err->error) goto cleanup;
    
  } /* end normalization loop */
  
    /* DEBUG
    ObitImageUtilArray2Image ("DbugRawBeam1.fits",  0, array[0], err);   */
    /* END DEBUG */

  /*  cleanup */
 cleanup:
  if (array) {
    for (i=0; i<narr; i++) {
      array[i] = ObitFArrayUnref( array[i]);
    }
    g_free(array);
  }
    if (err->error) Obit_traceback_msg (err, routine, out[0]->name);
} /* end ObitUVGridMFFFT2ImPar */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVGridMFClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVGridMFClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVGridMFClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVGridMFClassInfoDefFn (gpointer inClass)
{
  ObitUVGridMFClassInfo *theClass = (ObitUVGridMFClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit        = (ObitClassInitFP)ObitUVGridMFClassInit;
  theClass->ObitClassInfoDefFn   = (ObitClassInfoDefFnFP)ObitUVGridMFClassInfoDefFn;
  theClass->ObitGetClass         = (ObitGetClassFP)ObitUVGridMFGetClass;
  theClass->ObitClear            = (ObitClearFP)ObitUVGridMFClear;
  theClass->ObitInit             = (ObitInitFP)ObitUVGridMFInit;
  theClass->ObitUVGridSetup      = (ObitUVGridSetupFP)ObitUVGridMFSetup;
  theClass->ObitUVGridReadUV     = (ObitUVGridReadUVFP)ObitUVGridMFReadUV;
  theClass->ObitUVGridReadUVPar  = (ObitUVGridReadUVParFP)ObitUVGridMFReadUVPar;
  theClass->ObitUVGridFFT2Im     = (ObitUVGridFFT2ImFP)ObitUVGridMFFFT2Im;
  theClass->ObitUVGridFFT2ImPar  = (ObitUVGridFFT2ImParFP)ObitUVGridMFFFT2ImPar;


} /* end ObitUVGridMFClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVGridMFInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVGridMF *in = inn;

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
  in->grids     = NULL;
  in->BeamNorms = NULL;
  in->sigma1    = NULL;
  in->sigma2    = NULL;
  in->sigma3    = NULL;
} /* end ObitUVGridMFInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVGridMF* cast to an Obit*.
 */
void ObitUVGridMFClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong iSpec;
  ObitUVGridMF *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->grids) {
    for (iSpec=0; iSpec<in->nSpec; iSpec++) {
      in->grids[iSpec] = ObitCArrayUnref(in->grids[iSpec]);
    }
    g_free(in->grids);
  }
  if (in->BIFSpec)   g_free(in->BIFSpec);
  if (in->EIFSpec)   g_free(in->EIFSpec);
  if (in->BChanSpec) g_free(in->BChanSpec);
  if (in->EChanSpec) g_free(in->EChanSpec);
  if (in->BeamNorms) g_free(in->BeamNorms);
  if (in->sigma1)    g_free(in->sigma1);
  if (in->sigma2)    g_free(in->sigma2);
  if (in->sigma3)    g_free(in->sigma3);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVGridMFClear */

 /**
 * Grid a buffer load of data into a single image
 * \param in      Gridding Object
 * \param UVin    UV data set to grid from current buffer
 * \param sargs   Array of arguments to use, Must be in->nThreads of these
 * \param thread  Thread object to use
 * \param err     ObitErr stack for reporting problems.
 */
static void GridOne (ObitUVGridMF* in, ObitUV *UVin, UVGridFuncArg **sargs, 
		     ObitThread *thread, ObitErr *err)
{
  olong i, iSpec, nvis, nTh, nThread, mTh;
  ObitThreadFunc func=(ObitThreadFunc)ThreadUVGridMFBuffer ;
  UVGridFuncArg *args;
  gboolean  OK;
  gsize size;
  gchar *routine="UVGridMF:GridOne";

  /* error checks */
  if (err->error) return;

  /* To thread or not to */
  nvis = UVin->myDesc->numVisBuff;
  if (nvis<100) nThread = 1;
  else nThread = in->nThreads;
  nThread = MIN (nThread, in->nSpec);

  /* Need to copy data? Only one copy needed per call. */
  args = sargs[0];
  if ((args->buffer!=NULL) && (args->buffSize>0)) {
    size = (UVin->myDesc->numVisBuff)*UVin->myDesc->lrec*sizeof(ofloat);
    memcpy (UVin->buffer, args->buffer, size);
  }

  /* Loop over coarse channel processing nThread at a time */
  for (iSpec=0; iSpec<in->nSpec; iSpec+=nThread) {
    nTh = MIN (nThread, in->nSpec-iSpec);
    mTh = MIN (nTh, (in->nSpec-iSpec));
    
    for (i=0; i<mTh; i++) {
      args = sargs[i];
      if (mTh>1) args->ithread = i;
      else args->ithread = -1;
      args->thread = thread;
      args->in     = in;
      args->UVin   = UVin;
      args->BIF    = in->BIFSpec[iSpec+i];
      args->EIF    = in->EIFSpec[iSpec+i];
      args->BChan  = in->BChanSpec[iSpec+i];
      args->EChan  = in->EChanSpec[iSpec+i];
      args->grid   = in->grids[iSpec+i];
    } /* end setting up args */
    
  /* Do operation on buffer possibly with threads */
  OK = ObitThreadIterator (thread, mTh, func, (gpointer)sargs);
    
  /* Check for problems */
  if (!OK) 
    Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
  } /* end Loop over coarse channel  */

} /* end GridOne  */

/**
 * Prepares a buffer load of visibility data for gridding:
 * \li rotate (u,v,w) if doing 3D imaging and a shift.
 * \li shift position if needed.
 * \li if doBeam then replace data with (1,0).
 * \li enforce guardband - no data near outer edges of grid 
 * \li Convert to cells at the reference frequency.
 * \li All data  converted to the positive V half plane.
 *
 * MF version, threading split by coarse channel
 * \param in       Object with grid to accumulate.
 * \param uvdata   Object with uvdata in buffer.
 * \param BIF      Low IF 0-rel to process
 * \param EIF      High IF 0-rel to process
 * \param BChan    Low channel 0-rel to process
 * \param EChan    High channel 0-rel to process
 * \param uvw      array of UVW values 
 */
void PrepBufferMF (ObitUVGridMF* in, ObitUV *uvdata, olong BIF, olong EIF,
		   olong BChan, olong EChan, ObitFArray *uvw, 
		   ObitCArray *accGrid)
{
  olong ivis, nvis, ifreq, nif, iif, nfreq, ifq;
  ofloat *u, *v, *w, *vis, *ifvis, *vvis;
  ofloat phase, cp, sp, vr, vi, uu, vv, ww, uf, vf, wf;
  ofloat bl2, blmax2, blmin2, wt, guardu, guardv;
  ObitUVDesc *desc;
  olong fincf, fincif, luvw;
  gboolean doShift, doFlag, flip;
  ofloat tape, taperWt;
  gboolean doTaper;

  /* error checks */
  g_assert (ObitUVGridIsA(in));
  g_assert (ObitUVIsA(uvdata));
  g_assert (uvdata->myDesc != NULL);
  g_assert (uvdata->buffer != NULL);

  /* how much data? */
  luvw  = uvw->naxis[0];  /* length of uvw entry */
  desc  = uvdata->myDesc;
  nvis  = desc->numVisBuff;
  if (nvis<=0) return; /* need something */
  nfreq = desc->inaxes[desc->jlocf];
  nif = 1;
  if (desc->jlocif>=0) nif = desc->inaxes[desc->jlocif];
  
  /* Channel and IF increments in frequency scaling array */
  fincf  = MAX (1, (desc->incf  / 3) / desc->inaxes[desc->jlocs]);
  fincif = MAX (1, (desc->incif / 3) / desc->inaxes[desc->jlocs]);

  /* initialize data pointers into buffer */
  u   = uvw->array;
  v   = uvw->array+desc->ilocv-desc->ilocu;
  w   = uvw->array+desc->ilocw-desc->ilocu;
  vis = uvdata->buffer + desc->nrparm + desc->incf*BChan + desc->incif*BIF;

  /* what needed */
  doShift = (in->dxc!=0.0) || (in->dyc!=0.0) || (in->dzc!=0.0);
  doShift = doShift && (!in->doBeam); /* no shift for beam */

  /* Baseline max, min values */
  blmax2 = in->blmax * in->blmax;
  blmin2 = in->blmin * in->blmin;

  /* guardband in wavelengths */
  guardu = ((1.0-in->guardband) * (ofloat)accGrid->naxis[0]) / fabs(in->UScale);
  guardv = ((1.0-in->guardband) * ((ofloat)accGrid->naxis[1])/2) / fabs(in->VScale);

  /* Tapering?  */
  doTaper = (fabs(in->sigma1[nfreq*nif/2])>1.0e-20) || (in->BeamTaperUV!=0.0);

  /* Loop over visibilities */
  for (ivis=0; ivis<nvis; ivis++) {

    /* check exterma */
    bl2 = (*u)*(*u) + (*v)*(*v);
    doFlag = ((bl2<blmin2) || (bl2>blmax2));

    /* rotate (u,v,w) if 3D */
    if (in->do3Dmul ) {
      uu = (*u)*in->URot3D[0][0] + (*v)*in->URot3D[0][1] + (*w)*in->URot3D[0][2];
      vv = (*u)*in->URot3D[1][0] + (*v)*in->URot3D[1][1] + (*w)*in->URot3D[1][2];
      ww = (*u)*in->URot3D[2][0] + (*v)*in->URot3D[2][1] + (*w)*in->URot3D[2][2];
      *u = uu;
      *v = vv;
      *w = ww;
    } /* end rotate u,v,w */
    
    /* in the correct half plane? */
    flip = (*u) <= 0.0;

    /* loop over IFs */
    ifvis = vis;
    for (iif=BIF; iif<=EIF; iif++) {

      /* loop over frequencies */
      vvis = ifvis;
      for (ifreq=BChan; ifreq<=EChan; ifreq++) {
	ifq = iif*fincif + ifreq*fincf;  /* index in IF/freq table */

	/* Scale coordinates to frequency */
	uf = *u * desc->fscale[ifq];
	vf = *v * desc->fscale[ifq];
	wf = *w * desc->fscale[ifq];

	/* Channel tapering weight */
	if (doTaper) {
	  /*Sfreq = Starget * tarFreq/desc->freqArr[ifq]; */
	  /* Scorr = sqrt(MAX (1.0e-10,Starget*Starget - Sfreq*Sfreq));  Correction size */
	  /*taper = 0.95 / (Scorr*tapeConst) ;*/
	  /* taper = (0.8/ ((Scorr/2.35)/206265.))/(G_PI); */
	  /* sigma2  = log(0.3)/(taper*taper); */
	  tape = (uf*uf*(in->sigma2[ifq]+in->BeamTaperUV) + 
		  vf*vf*(in->sigma1[ifq]+in->BeamTaperUV) + 
		  uf*vf*(in->sigma3[ifq]));
	  if (tape<-14.0) taperWt = 0.0; /* underflow */
	  else taperWt = exp(tape);
	  vvis[2] *= taperWt;  /* Apply to weight */
	}

	/* is this one wanted? */
	if (doFlag)  vvis[2] = 0.0;  /* baseline out of range? */
	wt = vvis[2];                /* data weight */
	if (wt <= 0.0) {vvis += desc->incf; continue;}
	
	/* shift position if needed */
	if (doShift) {
	  phase = (uf*in->dxc + vf*in->dyc + wf*in->dzc);
	  cp = cos(phase);
	  sp = sin(phase);
	  vr = vvis[0];
	  vi = vvis[1];
	  /* rotate phase of visibility */
	  vvis[0] = cp * vr - sp * vi;
	  vvis[1] = sp * vr + cp * vi;
	}
	
	/* Making a beam - if so replace data with (1,0) */
	if (in->doBeam) {
	  vvis[0] = 1.0;
	  vvis[1] = 0.0;
	}
	
	/* conjugate phase if needed */
	if (flip)  vvis[1] = - vvis[1];
	
	/* enforce guardband */
	if ((fabs(uf)>guardu) || (fabs(vf)>guardv)) vvis[2] = 0.0;
	
	vvis += desc->incf; /* visibility pointer */
      } /* end loop over frequencies */
      ifvis += desc->incif; /* visibility pointer */
    } /* Loop over IFs */

    /* Scale u,v,w to cells at reference frequency */
    if (flip) { /* put in other half plane */
      *u = -((*u) * in->UScale);
      *v = -((*v) * in->VScale);
      *w = -((*w) * in->WScale);
    } else { /* no flip */
      *u *= in->UScale;
      *v *= in->VScale;
      *w *= in->WScale;
    }

    /* update data pointers */
    u   += luvw;
    v   += luvw;
    w   += luvw;
    vis += desc->lrec;
  } /* end loop over visibilities */
} /* end PrepBufferMF */

/**
 * Convolves data in buffer on uvdata onto accGrid
 * Rows in the grid are in U and the data should have all been converted to the 
 * positive U half plane.
 * U, V, and W should be in cells and data not to be included on the grid should 
 * have zero weight.  Convolution functions must be created.
 * Details of data organization are set by FFTW, the zero v row is first and v=-1
 * row is last.
 *
 * MF version, threading split by coarse channel
 * \param inn      UVGrid Object 
 * \param uvdata   Object with uv data in buffer, prepared for gridding.
 * \param BIF      Low IF 0-rel to process
 * \param EIF      High IF 0-rel to process
 * \param BChan    Low channel 0-rel to process
 * \param EChan    High channel 0-rel to process
 * \param uvw      array of UVW values 
 * \param accGrid  Accumulation Grid
 */
void GridBufferMF (ObitUVGridMF* inn, ObitUV *uvdata, olong BIF, olong EIF,
		   olong BChan, olong EChan, ObitFArray *uvw, 
		   ObitCArray *accGrid)
{
  ObitUVGridMF *in = (ObitUVGridMF*)inn;
  olong ivis, nvis, ifreq, nfreq, ncol=0, iu, iv, iuu, ivv, icu, icv;
  olong iif, nif, ifq, lGridRow, lGridCol, itemp;
  ofloat *grid, *ggrid, *cconvu, *convu, *convv, *cconvv, *u, *v, *w, *vis, *vvis, *ifvis, *wt;
  ofloat *convfnp, *gridStart, *gridTop, visWtR, visWtI, visWtVR, visWtVI, rtemp, xtemp;
  ofloat uf, vf;
  olong fincf, fincif, luvw;
  olong pos[] = {0,0,0,0,0};
  ObitUVDesc *desc;

  /* error checks */
  g_assert (ObitUVGridMFIsA(in));
  g_assert (ObitUVIsA(uvdata));
  g_assert (uvdata->myDesc != NULL);
  desc  = uvdata->myDesc;
  nvis  = desc->numVisBuff;
  if (nvis<=0) return; /* need something */
  g_assert (uvdata->buffer != NULL);

  /* how much data? */
  luvw  = uvw->naxis[0];  /* length of uvw entry */
  nfreq = desc->inaxes[desc->jlocf];
  nif = 1;
  if (desc->jlocif>=0) nif = desc->inaxes[desc->jlocif];
 
  /* Channel and IF increments in frequency scaling array */
  fincf  = MAX (1, (desc->incf  / 3) / desc->inaxes[desc->jlocs]);
  fincif = MAX (1, (desc->incif / 3) / desc->inaxes[desc->jlocs]);

  /* initialize data pointers */
  u   = uvw->array;
  v   = uvw->array+desc->ilocv-desc->ilocu;
  w   = uvw->array+desc->ilocw-desc->ilocu;
  vis = uvdata->buffer + desc->nrparm + desc->incf*BChan + desc->incif*BIF;

  lGridRow = 2*accGrid->naxis[0]; /* length of row as floats */
  lGridCol = accGrid->naxis[1];   /* length of column */

  /* beginning of the grid */
  pos[0] = 0;  pos[1] = 0;
  gridStart = ObitCArrayIndex (accGrid, pos); 
  /* beginning of highest row */
  pos[1] = lGridCol-1;
  gridTop = ObitCArrayIndex (accGrid, pos); 

  /* convolution fn pointer */
  pos[0] = 0; pos[1] = 0;
  convfnp = ObitFArrayIndex (in->convfn, pos);

  /* Loop over visibilities */
  for (ivis=0; ivis<nvis; ivis++) {

    /* loop over IFs */
    ifvis = vis;
    for (iif=BIF; iif<=EIF; iif++) {

      /* loop over frequencies */
      vvis = ifvis;
      for (ifreq=BChan; ifreq<=EChan; ifreq++) {
	ifq = iif*fincif + ifreq*fincf;  /* index in IF/freq table */

	/* is this one wanted? */
	wt = vvis + 2; /* data weight */
	if (*wt <= 0.0) {vvis += desc->incf; continue;}

	/* data times weight */
	visWtR = vvis[0] * (*wt);
	visWtI = vvis[1] * (*wt);
	
	/* Scale u,v for frequency (w not used) */
	uf = *u * desc->fscale[ifq];
	vf = *v * desc->fscale[ifq];
	
	/* get center cell */
	if (vf > 0.0) iv = (olong)(vf + 0.5);
	else iv = (olong)(vf - 0.5);
	iu = (olong)(uf + 0.5);

	/* back off half Kernel width */
	iu -= in->convWidth/2;
	iv -= in->convWidth/2;
	
	/* Starting convolution location, table has in->convNperCell points per cell */
	/* Determine fraction of the cell to get start location in convolving table. */
	if (uf > 0.0) itemp = (olong)(uf + 0.5);
	else itemp = ((olong)(uf - 0.5));
	xtemp = in->convNperCell*(itemp - (uf) - 0.5);
	if (xtemp > 0.0) xtemp += 0.5;
	else xtemp -= 0.5;
	convu = convfnp + in->convNperCell + (olong)xtemp;
	
	/* now v convolving fn */
	if (vf > 0.0) itemp = (olong)(vf + 0.5);
	else itemp = ((olong)(vf - 0.5));
	rtemp = in->convNperCell*(itemp - (vf) - 0.5);
	if (rtemp > 0.0) rtemp += 0.5;
	else rtemp -= 0.5;
	convv = convfnp + in->convNperCell + (olong)rtemp;
	
	/* if too close to the center, have to break up and do conjugate halves */
	if (iu >= 0) { /* all in same half */
	  ncol = in->convWidth; /* Complex addressed as floats */
	  pos[0] = iu;
	  /* Do v center at the edges */
	  if (iv>=0) pos[1] = iv;
	  else pos[1] = iv + lGridCol;
	  
	} else { 
	  /* have to split - grid part in conjugate half */
	  /* FFTW uses only half of the first plane so split U */
	  iuu = -iu; /* hermitian */
	  ivv = -iv;
	  pos[0] = iuu;
	  /* Do v center at the edges */
	  if (ivv>=0) pos[1] = ivv;
	  else pos[1] = ivv + lGridCol;
	  /*grid = ObitCArrayIndex (accGrid, pos); pointer in grid */ 
	  if ((pos[0]<=accGrid->naxis[0]) && (pos[1]<=accGrid->naxis[1]) &&
	      (pos[0]>=0) && (pos[1]>=0) )
	    grid = accGrid->array+2*(pos[1]*accGrid->naxis[0]+pos[0]);
	  else
	    grid = NULL;

	  /* Ignore if outside grid */
	  if (grid!=NULL) {
	    ncol = iuu;
	    cconvv = convv;
	    for (icv=0; icv<in->convWidth; icv++) {
	      cconvu = convu;
	      visWtVR = visWtR * (*cconvv);
	      visWtVI = visWtI * (*cconvv);
	      /* Trickery with the v row to get data in the correct place for the FFT 
		 the following will only be triggered if the iv wraps */
	      if ((pos[1]-icv)==-1) {
		grid = gridTop+2*iuu; /* top of grid */
	      }
	      ggrid  = grid;
	      for (icu=0; icu<=ncol; icu++) {
		ggrid[0]   += visWtVR * (*cconvu);
		ggrid[1]   -= visWtVI * (*cconvu); /* conjugate */
		cconvu += in->convNperCell;  /* U Convolution kernel pointer */
		ggrid -= 2; /* gridding pointer - opposite of normal gridding */
	      } /* end inner u gridding loop */
	      cconvv += in->convNperCell;  /* V Convolution kernel pointer */
	      grid -= lGridRow; /* gridding pointer - reverse direction for conjugate */
	    } /* end outer v loop */
	    
	    /* set up for rest of grid */
	    ncol = (in->convWidth + iu); /* how many columns left? */
	    iu = 0;      /* by definition  start other half plane at iu=0 */
	    pos[0] = iu; 
	    /* Do v center at the edges */
	    if (iv>=0) pos[1] = iv;
	    else pos[1] = iv + lGridCol;
	    convu = convu + iuu * in->convNperCell; /* for other half in u */
	  } /* end if in grid */
	} /* End of dealing with conjugate portion */
	  
	/* main loop gridding - only if in grid */
	/*grid = ObitCArrayIndex (accGrid, pos);   pointer in grid */
	if ((pos[0]<=accGrid->naxis[0]) && (pos[1]<=accGrid->naxis[1]) &&
	    (pos[0]>=0) && (pos[1]>=0) )
	  grid = accGrid->array+2*(pos[1]*accGrid->naxis[0]+pos[0]);
	else
	  grid = NULL;

	if (grid!=NULL) {
	  for (icv=0; icv<in->convWidth; icv++) {
	    cconvu = convu;
	    visWtVR = visWtR * (*convv);
	    visWtVI = visWtI * (*convv);
	    /* Trickery with the v row to get data in the correct place for the FFT 
	       the following will only be triggered if the iv row goes non negative */
	    if ((iv<0) && ((iv+icv)==0)) grid = gridStart+2*iu; /* beginning of grid */
	    ggrid  = grid;
	    for (icu=0; icu<ncol; icu++) {
	      ggrid[0] += visWtVR * (*cconvu);  /* real */
	      ggrid[1] += visWtVI * (*cconvu) ; /* imag */

	      /* Hard core debug
	      if (ggrid-gridStart==72) {
		fprintf (stdout," reglr %10.5f %10.5f %3ld %10.5f %10.5f %15.5f %15.5f %5ld %5ld  %d\n",
			 uf, vf, ifq, visWtVR*(*cconvu), visWtVI*(*cconvu),ggrid[0],ggrid[1],icu,icv,ggrid-gridStart);
	      } */
	      cconvu += in->convNperCell;  /* Convolution kernel pointer */
	      ggrid += 2; /* gridding pointer */
	    } /* end inner gridding loop */
	    convv += in->convNperCell;  /* Convolution kernel pointer */
	    grid += lGridRow; /* gridding pointer */
	  } /* end outer gridding loop */
	} /* end if in grid */
	vvis += desc->incf; /* visibility pointer */
	
      } /* end loop over frequencies */
      ifvis += desc->incif; /* visibility pointer */
    } /* Loop over IFs */
    
    /* update data pointers */
    u   += luvw;
    v   += luvw;
    w   += luvw;
    vis += desc->lrec;
  } /* end loop over visibilities */
} /* end GridBufferMF */

/** 
 * Prepare and Grid a portion of the data buffer
 * Arguments are given in the structure passed as arg
 * Note the images and beams are not normalized.
 * \param arg  Pointer to UVGridFuncArg argument with elements
 * \li thread Thread with restart queue
 * \li in     ObitUVGrid object
 * \li UVin   UV data set to grid from current buffer
 * \li BIF    Low IF 0-rel to process
 * \li EIF    High IF 0-rel to process
 * \li BChan  Low channel 0-rel to process
 * \li EChan  High channel 0-rel to process
 * \li grid   Array of nSpec gridding arrays
 * \li uvw    Array of u,v,w values
 * \li ithread thread number, >0 -> no threading 
 * \li buffSize if >0 then the number of ofloats to copy from buffer to buffer on UVin
 * \li buffer   Data buffer to copy
 */
static gpointer ThreadUVGridMFBuffer (gpointer arg)
{
  /* Get arguments from structure */
  UVGridFuncArg *largs = (UVGridFuncArg*)arg;
  ObitUVGridMF *in = largs->in;
  ObitUV *UVin     = largs->UVin;
  olong BIF        = largs->BIF;
  olong EIF        = largs->EIF;
  olong BChan      = largs->BChan;
  olong EChan      = largs->EChan;
  ObitCArray *grid = largs->grid;
  ObitFArray  *uvw = largs->uvw;

 /* Make copy of u,v,w */
  copyUVW(largs);   
  
  /* prepare data */
  PrepBufferMF (in, UVin, BIF, EIF, BChan, EChan, uvw, grid);
  
  /* grid */
  GridBufferMF (in, UVin, BIF, EIF, BChan, EChan, uvw, grid);

  /* Indicate completion */
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadUVGridBuffer */

/** 
 * Reorders grid and does gridding correction.
 * NOTE: threading in FFTW apparently conflicts with Obit threads.
 * Arguments are given in the structure passed as arg
 * \param arg  Pointer to FFT2ImFuncArg argument with elements
 * \li thread    Thread with restart queue
 * \li in        Input ObitUVGrid object
 * \li array     Output ObitFArray Image array
 * \li ithread   thread number, >0 -> no threading 
 */
static gpointer ThreadFFT2ImMF (gpointer arg)
{
  /* Get arguments from structure */
  FFT2ImFuncArg *largs = (FFT2ImFuncArg*)arg;
  ObitUVGrid *in     = largs->in;
  ObitFArray *array  = largs->array;

  /* reorder to center at center */
  ObitFArray2DCenter (array);
  
  /* Do multiply to make griding correction */
  if (in->doBeam) 
    ObitFArrayMulColRow (array, in->xCorrBeam, in->yCorrBeam, array);
  else
    ObitFArrayMulColRow (array, in->xCorrImage, in->yCorrImage, array);
  goto finish;
  
  /* cleanup */
 finish:
  /* Indicate completion */
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadFFT2ImMF */

/** 
 * Copy u,v,w from uv buffer on UVin to uvw member
 */
static void copyUVW (UVGridFuncArg *arg)
{
  olong i, n, lrec;
  ObitUVDesc *desc = arg->UVin->myDesc;
  ofloat *u, *v, *w, *uvw;
  
  uvw  = arg->uvw->array;
  u    = arg->UVin->buffer + desc->ilocu;
  v    = arg->UVin->buffer + desc->ilocv;
  w    = arg->UVin->buffer + desc->ilocw;
  n    = MIN(desc->numVisBuff, arg->uvw->naxis[1]);
  lrec = desc->lrec;
  
  for (i=0; i<n; i++) {
    *(uvw)   = *u;
    *(uvw+1) = *v;
    *(uvw+2) = *w;
    u   += lrec;
    v   += lrec;
    w   += lrec;
    uvw += 3;
  } /* end copy loop */
} /* end copyUVW */
