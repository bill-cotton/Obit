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
#include <unistd.h>
#include "ObitFFT.h"
#include "ObitUVGridMF.h"
#include "ObitImageUtil.h"
#include "ObitThreadGrid.h"

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
/** gridding correction threaded function argument */
typedef struct {
  /* ObitThread with restart queue */
  ObitThread *thread;
  /* SkyModel with model components loaded (ObitSkyModelLoad) */
  ObitUVGrid *in;
  /* UV data set to model and subtract from current buffer */
  ObitFArray *array;
  /* plane number 0-rel  */
  olong       iplane;
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

/** Private: Threaded FFT/gridding correct */
static gpointer ThreadFFT2ImMF (gpointer arg);

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
 * Wideband, multi frequency version
 * Input data should be fully edited and calibrated, with any weighting applied 
 * and converted to the appropriate Stokes type.
 * The object UVin will be opened during this call if it is not already open.
 * image should describe the center, size and grid spacing of the desired
 * image.
 * The beams corresponding to each image should be made first using the
 * same ObitUVGridMF.
 * If imagee is of type ObitImageMF then the selection of IFs and channels per 
 * output channel are obtained from this object, otherwise, the imagee info member 
 * is expected to contain
 * \li "nSpec"     OBIT_long scalar = number of spectral channels, imagee should be 
 *                 adequately dimensioned.
 * \li "BIFSpec"   OBIT_long(nSpec)  = Beginning IF 0-rel in UVIn per spectral channel
 * \li "EIFSpec"   OBIT_long(nSpec)  = Endinging IF 0-rel in UVIn per spectral channel
 * \li "BChanSpec" OBIT_long(nSpec)  = Beginning channel 0-rel in UVIn per spectral channel
 * \li "EChanSpec" OBIT_long(nSpec)  = Endinging channel 0-rel in UVIn per spectral channel
 * \param in       Object to initialize
 * \param UVin     Uv data object to be gridded.
 *                 info entry "MFTaper" gives frequency dependent taper fudge factor
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
  ObitImage *image = (ObitImage*)imagee;
  ObitImage *myBeam;
  olong nx, ny, naxis[2], iSpec, size, nif, nfreq, nn;
  ofloat cellx, celly, dxyzc[3], xt, yt, zt, MFTape=0.8, BeamTaper=0.0;
  olong nSpec=1, maxOrder=-1, *BIFSpec=NULL, *EIFSpec=NULL, *BChanSpec=NULL, *EChanSpec=NULL;
  ofloat *ramp=NULL, *data=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean doCalSelect = FALSE;
  ObitIOAccess access;
  ofloat Beam[3] = {0.0,0.0,0.0}, tarBeam[3]={0.,0.,0.}, corrBeam[3]={0.,0.,0.};
  ofloat sigma2v, sigma2u, cpa, spa, taper;
  odouble tarFreq;
  gchar *routine="ObitUVGridMFSetup";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVGridMFIsA(in));
  g_assert (ObitUVIsA(UVin));
  g_assert (ObitImageIsA(image));
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

  /* Spectral information */
  if (ObitImageMFIsA(image)) {  /* MF Image */
      nSpec      = ((ObitImageMF*)image)->nSpec;     /* save number of coarse channels */;
      maxOrder   = ((ObitImageMF*)image)->maxOrder;  /* save max imaging order */
      BIFSpec    = ((ObitImageMF*)image)->BIFSpec;
      EIFSpec    = ((ObitImageMF*)image)->EIFSpec;
      BChanSpec  = ((ObitImageMF*)image)->BChanSpec;
      EChanSpec  = ((ObitImageMF*)image)->EChanSpec;
    } else {                   /* Spectral cube */
      ObitInfoListGetTest(image->info, "nSpec",     &type, dim, &nSpec);
      ObitInfoListGetP(image->info, "BIFSpec",   &type, dim, (gpointer*)&BIFSpec);
      ObitInfoListGetP(image->info, "EIFSpec",   &type, dim, (gpointer*)&EIFSpec);
      ObitInfoListGetP(image->info, "BChanSpec", &type, dim, (gpointer*)&BChanSpec);
      ObitInfoListGetP(image->info, "EChanSpec", &type, dim, (gpointer*)&EChanSpec);
      /* Check consistency with image size */
      Obit_return_if_fail((image->myDesc->inaxes[image->myDesc->jlocf]>=nSpec), err,
			  "%s: Image %s inadequately dimensioned for %d channels", 
			  routine, image->name, nSpec);
    }

    /* Spectral stuff defined? */
    Obit_return_if_fail(((BIFSpec!=NULL)   && (EIFSpec!=NULL) &&
			 (BChanSpec!=NULL) && (EChanSpec!=NULL)), err,
			"%s: Channel selection not given", 
			routine);
    in->nSpec    = nSpec;
    in->maxOrder = maxOrder;

  
  /* open uv data to fully instantiate if not already open */
  if (in->myStatus==OBIT_Inactive) {
    retCode = ObitUVOpen (UVin, access, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  uvDesc = UVin->myDesc;

  /* Get source position if it's not already in header */
  if ((uvDesc->crval[uvDesc->jlocr]==0.0) && 
      (uvDesc->crval[uvDesc->jlocd]==0.0)) {
    ObitUVGetRADec (UVin, &uvDesc->crval[uvDesc->jlocr], 
			&uvDesc->crval[uvDesc->jlocd], err);
    if (err->error) Obit_traceback_msg (err, routine, UVin->name);
  }

  /* Beam, image dependent stuff */
  in->nxBeam     = myBeam->myDesc->inaxes[0];
  in->nyBeam     = myBeam->myDesc->inaxes[1];
  in->icenxBeam  = in->nxBeam/2 + 1; 
  in->icenyBeam  = in->nyBeam/2 + 1;
  in->nxImage    = image->myDesc->inaxes[0];
  in->nyImage    = image->myDesc->inaxes[1];
  in->icenxImage = in->nxImage/2 + 1;
  in->icenyImage = in->nyImage/2 + 1;
   
  /* Frequency dependent fudge factor */
  ObitInfoListGetTest(UVin->info, "MFTaper", &type, dim, &MFTape);

  /* Any additional tapering (deg) */
  ObitInfoListGetTest(image->myDesc->info, "BeamTapr", &type, dim, &BeamTaper);
  if (BeamTaper>0.0) {
    /* DEBUG taper   = (1.0 / (BeamTaper*DG2RAD/2.35)/(G_PI));
       in->BeamTaperUV = log(0.3)/(taper*taper);*/
    taper = BeamTaper*DG2RAD*sqrt(2.0)*G_PI/2.35482;
    taper = BeamTaper*DG2RAD*sqrt(2.0)*G_PI / 1.17741022; /* DEBUG */
    in->BeamTaperUV = -taper*taper;
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
    in->BIFSpec[iSpec]   = BIFSpec[iSpec];
    in->EIFSpec[iSpec]   = EIFSpec[iSpec];
    in->BChanSpec[iSpec] = BChanSpec[iSpec];
    in->EChanSpec[iSpec] = EChanSpec[iSpec];
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

   /* Get target restoring beam */
  if (!ObitInfoListGetTest(UVin->info, "targBeam", &type, dim, Beam)) 
    ObitInfoListGetTest(UVin->info, "Beam", &type, dim, Beam);
  tarFreq = UVin->myDesc->freqIF[0];  /* Frequency of lowest IF */

  for (iSpec=0; iSpec<nn; iSpec++) {
    if (MFTape>0) {
      /* Target beam size scaled to this frequency */
      tarBeam[0] = Beam[0] * tarFreq/uvDesc->freqArr[iSpec];
      tarBeam[1] = Beam[1] * tarFreq/uvDesc->freqArr[iSpec];
      tarBeam[2] = Beam[2];
      /* Correction beam including any additional */
      corrBeam[0] = sqrt (MAX(1.0e-10,Beam[0]*Beam[0]-tarBeam[0]*tarBeam[0]) + BeamTaper);
      corrBeam[1] = sqrt (MAX(1.0e-10,Beam[1]*Beam[1]-tarBeam[1]*tarBeam[1]) + BeamTaper);
      corrBeam[2] = Beam[2];
      /* MFTape = beam fudge factor */
      taper   = (MFTape/ (((corrBeam[0]/2.35)/206265.))/(G_PI));
      sigma2u = log(0.3)/(taper*taper);
      taper   = (MFTape/ (((corrBeam[1]/2.35)/206265.))/(G_PI));
      sigma2v = log(0.3)/(taper*taper);
      cpa     = cos(corrBeam[2]*DG2RAD);
      spa     = sin(corrBeam[2]*DG2RAD);
      in->sigma1[iSpec]  = -fabs((cpa*cpa*sigma2v + spa*spa*sigma2u));
      in->sigma2[iSpec]  = -fabs((spa*spa*sigma2v + cpa*cpa*sigma2u));
      in->sigma3[iSpec]  = -fabs(2.0*cpa*spa*(sigma2v - sigma2u));
    } else {
      /* No frequency dependent taper */
      in->sigma1[iSpec]  = 0.0;
      in->sigma2[iSpec]  = 0.0;
      in->sigma3[iSpec]  = 0.0; 
    }
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
  olong   itemp;
  ObitThreadGrid *grids=NULL;
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

   /* Create threaded gridding object */
  grids = newObitThreadGrid("thread gridder");
  /* Init */
  ObitThreadGridSetup (grids, UVin, 1, &inn, in->nThreads, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* end initialize */

  /* loop gridding data */
  while (retCode == OBIT_IO_OK) {
    
    /* read buffer */
    if (doCalSelect) retCode = ObitUVReadSelect (UVin, NULL, err);
    else retCode = ObitUVRead (UVin, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    /* Do operation on buffer possibly with threads to grid all */
    ObitThreadGridGrid(grids);
  } /* end loop reading/gridding data */

  /* Shut down any threading */
  ObitThreadPoolFree (in->thread);

  /* fold negative u columns to conjugate cells */
  ObitThreadGridFlip(grids);
  
  /* Accumulate grids to "grid" members on in swapped for FFT */
  ObitThreadGridMerge(grids);
  
  /* release ThreadGrid object  */
  grids = ObitThreadGridUnref(grids);

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
 * \param UVin    UV data objects to be gridded.
 *                Should be the same as passed to previous call to 
 *                #ObitUVGridSetup for input in element.
 * \param err     ObitErr stack for reporting problems.
 */
void ObitUVGridMFReadUVPar (olong nPar, ObitUVGrid **inn, ObitUV *UVin, ObitErr *err)
{
  ObitUVGridMF **in = (ObitUVGridMF**)inn;
  ObitIOCode retCode = OBIT_IO_OK;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ofloat temp;
  olong ip, itemp;
  olong nTh,  nCopy, doCalib;
  ObitThreadGrid *grids=NULL;
  gboolean doCalSelect;
  gchar *routine="ObitUVGridMFReadUVPar";
  /* DEBUG 
  ObitFArray *dbgRArr=NULL, *dbgIArr=NULL;*/

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(UVin));
  g_assert (ObitUVDescIsA(UVin->myDesc));
  g_assert (UVin->myDesc->fscale!=NULL); /* frequency scaling table */
  if (nPar<=0) return;
  for (ip=0; ip<nPar; ip++) {
    g_assert (ObitUVGridIsA(in[ip]));
  }

  /*  ObitErrTimeLog(err, routine);  Add Timestamp */

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
  ObitInfoListGetTest(UVin->info, "doCalib", &type, dim, &doCalib);
  doCalSelect = FALSE;
  ObitInfoListGetTest(UVin->info, "doCalSelect", &type, dim, &doCalSelect);

  /* UVin[0] should have been opened in  ObitUVGridSetup */
  
  /* How many threads? */
  in[0]->nThreads = MAX (1, ObitThreadNumProc(in[0]->thread));

  /* How many threads? */
  nTh = in[0]->nThreads;

  /* Create threaded gridding object */
  grids = newObitThreadGrid("thread gridder");

  /* Init */
  ObitThreadGridSetup (grids, UVin, nPar, (ObitUVGrid**)inn, in[0]->nThreads, err);
  if (err->error) goto cleanup;

  ObitErrLog(err);

  /* loop gridding data */
  while (retCode == OBIT_IO_OK) {
    /* read buffer  */
    nCopy = 1;
    if (doCalSelect) retCode = ObitUVReadSelect (UVin, NULL, err);
    else             retCode = ObitUVRead (UVin, NULL, err);
    if (retCode==OBIT_IO_EOF) break;  /* Finished */
    if (err->error) goto cleanup;
    
    /* Do operation on buffer possibly with threads to grid all */
    ObitThreadGridGrid(grids);
  } /* end loop reading/gridding data */

  /* Shut down any threading */
  ObitThreadPoolFree (grids->GridInfo->thread);
  
  /* fold negative u columns to conjugate cells */
  ObitThreadGridFlip(grids);

  /* Accumulate grids to "grid" members on in swapped for FFT */
  ObitThreadGridMerge(grids);

  /* release ThreadGrid object  */
  grids = ObitThreadGridUnref(grids);

  /*ObitErrTimeLog(err, "Stop Grid Loop");  DEBUG */
  ObitErrLog(err);

  /* DEBUG - look at first complex grid
  dbgRArr = ObitCArrayakeF(in[0]->grids[0]);
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

  /* Close data */
  retCode = ObitUVClose (UVin, err);

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
  gboolean allBlank;
   /* DEBUG
  ObitFArray *dbgRArr=NULL, *dbgIArr=NULL; */
 
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
    in->BeamNorm = 0.0;
    allBlank = TRUE;  /* All planes blanked */
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
      in->BeamNorm += in->BeamNorms[j]; /* track sum of weights */

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
      pln = 2+in->maxOrder+j;                             /* MF channel/plane number */
      if ((in->nSpec>1)  && (in->maxOrder<0)) pln = j+1;  /* Spectral imaging */
      if (in->nSpec==1)                       pln = 1;    /* Normal imaging */
      plane[0] = pln;
      ObitImagePutPlane (out, array->array, plane, err);
      if (err->error) goto cleanup;
      
    } /* end loop over spectral planes */

    /* Determine combined beam unless a spectral cube */
    if ((in->nSpec>1) && (ObitImageMFIsA(out))) 
      ObitImageMFCombine ((ObitImageMF*)out, FALSE, err);
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

    /* DEBUG - look at first complex grid 
    dbgRArr = ObitCArrayMakeF(in->grids[0]);
    dbgIArr = ObitCArrayMakeF(in->grids[0]);
    ObitCArrayReal (in->grids[0], dbgRArr);
    ObitCArrayImag (in->grids[0], dbgIArr);
    ObitImageUtilArray2Image ("DbugGridReal1.fits", 0, dbgRArr, err);  
    ObitImageUtilArray2Image ("DbugGridImag1.fits", 0, dbgIArr, err);  
    dbgRArr = ObitFArrayUnref(dbgRArr);
    dbgIArr = ObitFArrayUnref(dbgIArr);  */
    /*   end DEBUG */

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
    if ((in->nSpec>1)  && (ObitImageMFIsA(out))) 
      ObitImageMFCombine ((ObitImageMF*)out, TRUE, err);
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
  olong narr, nThread, plane[5]={1,1,1,1,1};
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  FFT2ImFuncArg *args=NULL;
  ObitUVGridMF **in = (ObitUVGridMF**)inn;
  ObitThreadFunc func=(ObitThreadFunc)ThreadFFT2ImMF;
  ObitFArray **array=NULL;
  ObitImageClassInfo *imgClass;
  ObitImage *theBeam;
  gboolean OK, allBlank;
  ofloat *Corrp, fblank = ObitMagicF();
  gchar *routine = "ObitUVGridMFFFT2ImPar";
  /* DEBUG
  ObitFArray *dbgRArr=NULL, *dbgIArr=NULL; */

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

  if (err->prtLv>=5) {  /* Diagnostics */
    Obit_log_error(err, OBIT_InfoErr, "%s: start FFTs",routine);
    ObitErrLog(err); 
  }

  /* How many threads for FFT? 2 per 1K pixels in x */
  nThread = MAX (1, ObitThreadNumProc(in[0]->thread));
  nThread = MIN (MAX (1, 2*in[0]->nxImage/1024), nThread);
  ObitFFTNThreads (nThread);   /* Enable FFT with threading */

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
      /* DEBUG - look at first complex grid
      if (i==1) {
	dbgRArr = ObitCArrayMakeF(in[i]->grids[0]);
	dbgIArr = ObitCArrayMakeF(in[i]->grids[0]);
	ObitCArrayReal (in[i]->grids[0], dbgRArr);
	ObitCArrayImag (in[i]->grids[0], dbgIArr);
	ObitImageUtilArray2Image ("DbugGridRealPar.fits", 0, dbgRArr, err);  
	ObitImageUtilArray2Image ("DbugGridImagPar.fits", 0, dbgIArr, err);  
	dbgRArr = ObitFArrayUnref(dbgRArr);
	dbgIArr = ObitFArrayUnref(dbgIArr);
      }  *//* end DEBUG */
      /* Create array of pixel arrays */
      for (j=0; j<in[i]->nSpec; j++) {
	array[ip] = ObitFArrayCreate("Image", 2, xdim);
	/* do FFT */
	ObitFFTC2R (in[i]->FFTImage, in[i]->grids[j], array[ip++]);
      }
    }
    
  } /* end loop doing FFTs */

  /* DEBUG  
  ObitImageUtilArray2Image ("DbugRawBeam0.fits",  0, array[0], err); */
  /* END DEBUG */

  if (err->prtLv>=5) {  /* Diagnostics */
    Obit_log_error(err, OBIT_InfoErr, "%s: finished FFTs",routine);
    ObitErrLog(err); 
  }

  /* Get normalization loop looking for an in entry with doBeam member set an,
   the center peak is measured and used to normalize,  this peak is assumed
   to be the normalization for the subsequent image.
   if an image without corresponding beam is encountered, the BeamNorms
   member of it in[] is used to normalize */


  for (i=0; i<nPar; i++) {
    /* is this a beam? */
    if (in[i]->doBeam) {
      pos[0] = 0; pos[1] = 0; pos[2] = 0;
      /* Loop over planes */
      in[i]->BeamNorm = 0.0;
      allBlank = TRUE;  /* All planes blanked */
      for (j=0; j<in[i]->nSpec; j++) {
	Corrp = ObitFArrayIndex(array[i*in[i]->nSpec+j], pos);
	in[i]->BeamNorms[j] = *Corrp;
	in[i]->BeamNorm += in[i]->BeamNorms[j]; /* track sum of weights */
	/* Deal with blanked planes */
	if (in[i]->BeamNorms[j]==0.0) in[i]->BeamNorms[j] = fblank;
	allBlank = allBlank &&  (in[i]->BeamNorms[j]==fblank);
	/* Save normalization on in[i+1] */
	if (!in[i+1]->doBeam) in[i+1]->BeamNorms[j] = in[i]->BeamNorms[j];
	in[i+1]->BeamNorm = in[i]->BeamNorm;
      } /* end loop over planes */
      /* Save BeamNorms array on Beam */
      dim[0] = in[i]->nSpec;  dim[1] = dim[2] = dim[3] = dim[4] = 1;
      ObitInfoListAlwaysPut(out[i]->info, "BeamNorms", OBIT_float, dim, in[i]->BeamNorms);
    } /* end if beam */
    /* Fetch BeamNorms from beam if needed */
    if ((!allBlank) && ((in[i]->BeamNorms[0]==0.0) || (in[i]->BeamNorms[0]==fblank))) {
      imgClass  = (ObitImageClassInfo*)out[i]->ClassInfo;    /* Image class */
      if (!in[i]->doBeam) {   /* This is not a beam */
	theBeam   = imgClass->ObitImageGetBeam(out[i], 0, plane, err);
	ObitInfoListGetTest(theBeam->info, "BeamNorms", &type, dim, in[i]->BeamNorms);
      } else {  /* Is a beam */
 	ObitInfoListGetTest(out[i]->info, "BeamNorms", &type, dim, in[i]->BeamNorms);
     }
      if (err->error) goto cleanup;
    } /* end fetch Beam Norms */

  } /* end normalization factor loop */

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

  /* Gridding correction 
     Loop over images - do nSpec planes in parallel using threading */
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
	args->iplane = (ii+off) % in[0]->nSpec;
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
    } /* end loop over channels */
  } /* End loop gridding correcting images */

  /* DEBUG 
  ObitImageUtilArray2Image ("DbugRawBeam1.fits",  0, array[0], err);  */
  /* END DEBUG */

   if (err->prtLv>=5) {  /* Diagnostics */
    Obit_log_error(err, OBIT_InfoErr, "%s: finished Gridding corr",routine);
    ObitErrLog(err); 
  }

  /* Loop writing and combining */
  for (i=0; i<nPar; i++) {
    out[i]->image = ObitFArrayUnref(out[i]->image);  /* Free buffer */
    /* is this a beam? */
    if (in[i]->doBeam) {
      /* Loop over planes */
      for (j=0; j<in[i]->nSpec; j++) {
	/* Check */
	if (in[i]->BeamNorms[j]==0.0) {  /* No data? */
	  Obit_log_error(err, OBIT_InfoWarn, 
			 "%s Peak in beam  %d is zero for:%s, blank plane",
			 routine, j, in[i]->name);
	  in[i]->BeamNorms[j] = fblank;
	  /* DEBUG
	  ObitImageUtilArray2Image ("DbugRawBeamBad.fits",  0, array[i*in[i]->nSpec+j], err);   */
	  /* END DEBUG */
	} /* end no data */
	
	pln = 2+in[i]->maxOrder+j;      /* Get channel/plane number */
	if (in[i]->nSpec==1) pln = 1;   /* Normal imaging */
	plane[0] = pln;
	ObitImagePutPlane (out[i], array[i*in[i]->nSpec+j]->array, plane, err);
	if (err->error) goto cleanup;
	out[i]->image = ObitFArrayUnref(out[i]->image);  /* Free buffer */
      } /* end loop over planes */

      if (err->prtLv>=5) {  /* Diagnostics */
	Obit_log_error(err, OBIT_InfoErr, "Start Combine Beam field %d",i+1);
      }
      /* Determine combined beam */
      if ((in[i]->nSpec>1) && (ObitImageMFIsA(out[i]))) 
	ObitImageMFCombine ((ObitImageMF*)out[i], FALSE, err);
      if (err->error) goto cleanup;
      if (err->prtLv>=5) {  /* Diagnostics */
	Obit_log_error(err, OBIT_InfoErr, "End Combine Beam");
	ObitErrLog(err); 
      }

      i++;       /* Advance to image or next beam */
    } /* end if beam */
    /*  Now image */

    /* Loop over planes */
    for (j=0; j<in[i]->nSpec; j++) {
      /* Deal with blanked planes */
      if (in[i]->BeamNorms[j]==0.0) in[i]->BeamNorms[j] = fblank;
      /* Write output */
      pln = 2+j+in[i]->maxOrder;       /* Get channel/plane number */
      if (in[i]->nSpec==1) pln = 1;    /* Normal imaging */
      plane[0] = pln;
      ObitImagePutPlane (out[i], array[i*in[i]->nSpec+j]->array, plane, err);
      if (err->error) goto cleanup;
    } /* end loop over planes */

    if (err->prtLv>=6) {  /* Diagnostics */
      Obit_log_error(err, OBIT_InfoErr, "Start Combine Image field %d",i+1);
    }
     /* Determine combined image */
    if ((in[i]->nSpec>1) && (ObitImageMFIsA(out[i]))) 
      ObitImageMFCombine ((ObitImageMF*)out[i], TRUE, err);
    if (err->error) goto cleanup;
    if (err->prtLv>=6) {  /* Diagnostics */
      Obit_log_error(err, OBIT_InfoErr, "End Combine Image");
      ObitErrLog(err); 
    }
    
  } /* end normalization loop */
  
    /* DEBUG
    ObitImageUtilArray2Image ("DbugRawBeam1.fits",  0, array[0], err);   */
    /* END DEBUG */

   if (err->prtLv>=5) {  /* Diagnostics */
    Obit_log_error(err, OBIT_InfoErr, "%s: finished imaging",routine);
    ObitErrLog(err); 
  }
 /*  cleanup */
 cleanup:
  if (array) {
    for (i=0; i<narr; i++) {
      array[i] = ObitFArrayUnref( array[i]);
    }
    g_free(array);
  }

  /* Delete FFT objects */
  for (i=0; i<nPar; i++) {
   if (in[i]->doBeam) { /* Beam? */
     in[i]->FFTBeam = ObitFFTUnref(in[i]->FFTBeam);
   } else { /* Image */
     in[i]->FFTImage =  ObitFFTUnref(in[i]->FFTImage);
   }
  }

  /* Reset FFT threading */
  ObitFFTNThreads (1);

  /* Clear FFT Threads */
  ObitFFTClearThreads();

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
  ObitUVGridMF *in     = (ObitUVGridMF*)largs->in;
  ObitFArray *array    = largs->array;
  olong iplane         = largs->iplane;

  olong  i;
  ObitFArray *YCorr   = NULL;
  ofloat norm, fblank = ObitMagicF();

  if ((in->BeamNorms[iplane]<=0.0) || (in->BeamNorms[iplane]==fblank)){
    ObitFArrayFill (array, fblank); /* Blank fill */
    goto finish;
  } else norm = 1.0/in->BeamNorms[iplane]; /* Normalization factor */
 
  /* reorder to center at center */
  ObitFArray2DCenter (array);
  
  /* Do multiply to make griding correction */
  if (in->doBeam) {
    /* Temporary array with normalization */
    YCorr = ObitFArrayCreate ("tmp", in->yCorrBeam->ndim, in->yCorrBeam->naxis);
    for (i=0; i<in->yCorrBeam->arraySize; i++) 
      YCorr->array[i] = in->yCorrBeam->array[i] * norm;
    ObitFArrayMulColRow (array, in->xCorrBeam, YCorr, array);
  } else {
    /* Temporary array with normalization */
    YCorr = ObitFArrayCreate ("tmp", in->yCorrImage->ndim, in->yCorrImage->naxis);
    for (i=0; i<in->yCorrImage->arraySize; i++) 
      YCorr->array[i] = in->yCorrImage->array[i] * norm;
    ObitFArrayMulColRow (array, in->xCorrImage, YCorr, array);
  }
  YCorr = ObitFArrayUnref(YCorr);
  goto finish;
  
  /* cleanup */
 finish:
  /* Indicate completion */
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadFFT2ImMF */

