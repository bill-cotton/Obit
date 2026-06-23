/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2026                                               */
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

#include "ObitGPUBeamModel.h"
#include "ObitSkyModelVMBeamMF.h"
#include "ObitSpectrumMF.h"
#if HAVE_GPU==1  /* GPU? Real versions */
#include "ObitCUDAUtil.h"
#include "CUDAMath.h"
#include "ObitGPUBeamInterp.h"
#include "ObitCUDABeamModelInfoDef.h"
#include "CUDASkyGeom.h"
#include "ObitUtil.h"
#endif /* HAVE_GPU */

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUBeamModel.c
 * ObitGPUBeamModel class function definitions.
 * This class is derived from the Obit base class.
 * The ObitGPUBeamModel class uses GPUs to calculate sky models.
 * CUDA is used for the actual GPU routines but CUDA is REALLY primitive
 * so very basic inplementation is in ObitCUDABeamModel.cu and friends.
 * Portions of the class are in CUDA and are only implemented if the
 * compiler option -DHAVE_GPU=1 is used.  Some portions also need to 
 * variable IS_CUDA=1 to be set in the calling routines.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitGPUBeamModel";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitGPUBeamModelClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGPUBeamModelClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGPUBeamModelInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGPUBeamModelClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGPUBeamModelClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitGPUBeamModel* newObitGPUBeamModel (gchar* name)
{
  ObitGPUBeamModel* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGPUBeamModelClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitGPUBeamModel));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGPUBeamModelInit((gpointer)out);

 return out;
} /* end newObitGPUBeamModel */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitGPUBeamModelGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGPUBeamModelClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGPUBeamModelGetClass */

/**
 * Make a deep copy of an ObitGPUBeamModel.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitGPUBeamModel* ObitGPUBeamModelCopy  (ObitGPUBeamModel *in, ObitGPUBeamModel *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

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
    out = newObitGPUBeamModel(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */

  return out;
} /* end ObitGPUBeamModelCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an GPUBeamModel similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitGPUBeamModelClone  (ObitGPUBeamModel *in, ObitGPUBeamModel *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;

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

} /* end ObitGPUBeamModelClone */

/**
 * Creates an ObitGPUBeamModel 
 * Currently only DFT point supported
 * \param name  An optional name for the object.
 * \param type  mode comp. type 'DFT' or 'GRID'.  Only DFT supported.
 * \return the new object.
 */
ObitGPUBeamModel* ObitGPUBeamModelCreate (gchar* name, gchar *type)
{
  ObitGPUBeamModel* out;

  /* Create basic structure */
  out = newObitGPUBeamModel (name);

  return out;
} /* end ObitGPUBeamModelCreate */

/**
 * Initialize an ObitGPUBeamModel 
 * Currently only DFT  supported
 * \param in       Input object
 * \param skyModel Skymodel to implement on GPU
 *                 disguised as an Obit to keep from confusing the compiler
 * \param uvdata   UV data being operated on
 * \param delta_PA maximum range of parallactic angle (deg)
 * \param err      Obit error stack object.
 */
void ObitGPUBeamModelDFTInit (ObitGPUBeamModel *in,  Obit *SkyModel,
			      ObitUV *uvdata, ofloat delta_PA, ObitErr *err)
{
#if HAVE_GPU==1  /* GPU? Real versions */
  ObitUVDesc *inDesc = uvdata->myDesc;
  ObitSkyModelVMBeamMF *skyModel = (ObitSkyModelVMBeamMF*)SkyModel;
  ObitImage *img;
  int i, memsize, kincf, kincif, nchan, nbchan, nif, nant, nmodel;
  ObitIOAccess access;
  ObitInfoType type;
  ObitImageDesc *beamDesc = NULL;
  gint32 dim[MAXINFOELEMDIM];
  gboolean doCalSelect;
  olong ldevice[20], nVisPIO, nSpec;
  odouble *specFreq, f_0;
  ofloat delt_f, crp_f;
  gchar *routine = "ObitGPUBeamModelDFTInit";

  /* reset GPU if haven't started already */
  if (!in->gpuInfo) ObitCUDAResetGPU ();
  
  /* gpuInfo - host only */
  memsize = sizeof(GPUBeamInfo);  /* Basic structure in locked host memory */
  /* allocate if needed */
  if (!in->gpuInfo) in->gpuInfo = (GPUBeamInfo*)ObitCUDAUtilAllocHost(memsize);

  /* Initialize visibility range handling */
  in->gpuInfo->delta_ParAng = (float)delta_PA;
  in->gpuInfo->timeRange[0] = -1.0e5; in->gpuInfo->timeRange[1] = 1.0e5;
  in->gpuInfo->curParAng = -5000.;
  in->gpuInfo->visRange[0] =  in->gpuInfo->visRange[1] = 0;
  in->gpuInfo->doUnit     = 0;  /* Use unit rather than beam Jones matrices? */
  in->gpuInfo->haveUnit   = 0;  /* Have unit matrices been initialized? */
  in->antList = ObitAntennaListRef(skyModel->AntList[0]);
  in->curSource = ObitSourceRef(skyModel->curSource);

  /* Set GPU device, default 0 */
  in->gpuInfo->cuda_device = 0;
  ldevice[0] = in->gpuInfo->cuda_device;
  //DEBUGObitInfoListGetTest(uvdata->info, "GPU_no", &type, (gint32*)dim, ldevice);
  in->gpuInfo->cuda_device = (int)ldevice[0];
  ObitCUDASetGPU (in->gpuInfo->cuda_device);
  // Seems to deallocate locked host memory ObitCUDAResetGPU ();  /* Make sure reset */

  /* Print level */
  in->gpuInfo->prtLv = err->prtLv;

  /* reallocate data buffer on uvdata to locked memory
    first reset number of vis per IO */
  ObitUVClose (uvdata, err);
  ObitInfoListGetTest (uvdata->info, "nVisPIO", &type, dim,  
		       &in->gpuInfo->oldnVisPIO);
  nVisPIO = MIN(10000, inDesc->nvis);
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (uvdata->info, "nVisPIO", OBIT_long, dim,  &nVisPIO);
  doCalSelect = FALSE;
  ObitInfoListGetTest(uvdata->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;   
  ObitUVOpen (uvdata, access, err);
  if (err->error) Obit_traceback_msg (err, routine, uvdata->name);

  /* Fill gpuInfo */
  in->gpuInfo->bufferSize = uvdata->bufferSize;
  memsize = (uvdata->bufferSize+2*inDesc->lrec) * sizeof(ofloat);
  /* allocate locked host memory - some (2 vis ) padding */
  in->gpuInfo->h_data = (ofloat*)ObitCUDAUtilAllocHost(memsize);
  
  in->gpuInfo->nvis    = nVisPIO;
  in->gpuInfo->nchan   = inDesc->inaxes[inDesc->jlocf]*inDesc->inaxes[inDesc->jlocif];
  in->gpuInfo->nrparm  = inDesc->nrparm;
  in->gpuInfo->lenvis  = inDesc->lrec;
  in->gpuInfo->nstream = BEAM_STREAM_COUNT;

  /* Allocate resources */
  /* Device Frequency scaling array from uv descriptor */
  nchan   = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocif>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else                   nif = 1;
  memsize = nchan * nif * sizeof(float);
  in->gpuInfo->d_freq = ObitCUDAUtilAllocGPU (memsize);  /* host not locked - risky? */
  ObitCUDAUtilHost2GPU (in->gpuInfo->d_freq, inDesc->fscale, memsize, NULL);

  /* Increments in frequency tables */
  if (inDesc->jlocif>=0) {
    if (inDesc->jlocf<inDesc->jlocif) { /* freq before IF */
      kincf  = 1;
      kincif = inDesc->inaxes[inDesc->jlocf];
    } else { /* IF before freq  */
      kincif = 1;
      kincf  = inDesc->inaxes[inDesc->jlocif];
    } 
  } else {  /* NO IF axis */
      kincif = 1;
      kincf  = 1;
  }
  /****    Visibility info    ****/
  memsize = sizeof(GPUVisInfo);
  in->gpuInfo->d_visInfo = (GPUBeamVisInfo *)ObitCUDAUtilAllocGPU(memsize);
  in->gpuInfo->h_visInfo = (GPUBeamVisInfo *)ObitCUDAUtilAllocHost(memsize);
  in->gpuInfo->h_visInfo->d_freqScale = in->gpuInfo->d_freq;
  in->gpuInfo->h_visInfo->nvis   = in->gpuInfo->nvis;
  in->gpuInfo->h_visInfo->nprod  = inDesc->ncorr;
  in->gpuInfo->h_visInfo->nrparm = in->gpuInfo->nrparm;
  in->gpuInfo->h_visInfo->lenvis = in->gpuInfo->lenvis;
  in->gpuInfo->h_visInfo->nchan  = in->gpuInfo->nchan;
  in->gpuInfo->h_visInfo->ilocu  = inDesc->ilocu;
  in->gpuInfo->h_visInfo->ilocv  = inDesc->ilocv;
  in->gpuInfo->h_visInfo->ilocw  = inDesc->ilocw;
  in->gpuInfo->h_visInfo->iloct  = inDesc->iloct;
  in->gpuInfo->h_visInfo->ilocb  = inDesc->ilocb;
  in->gpuInfo->h_visInfo->ilocsu = inDesc->ilocsu;
  in->gpuInfo->h_visInfo->chanb  = uvdata->mySel->startChann-1;;
  in->gpuInfo->h_visInfo->chane  = uvdata->mySel->startChann+uvdata->mySel->numberChann-2;
  in->gpuInfo->h_visInfo->incf   = inDesc->incf; 
  if (inDesc->jlocif>0) in->gpuInfo->h_visInfo->nIF = inDesc->inaxes[inDesc->jlocif];
  else                  in->gpuInfo->h_visInfo->nIF = 1;
  in->gpuInfo->h_visInfo->nSpec = 1;  /* Default */
  ObitInfoListGetTest (skyModel->mosaic->images[0]->myDesc->info, "NSPEC",
		       &type, dim, &in->gpuInfo->h_visInfo->nSpec);
  in->gpuInfo->h_visInfo->IFb    = uvdata->mySel->startIF-1;
  in->gpuInfo->h_visInfo->IFe    = uvdata->mySel->startIF+uvdata->mySel->numberIF-2;
  in->gpuInfo->h_visInfo->incif  = inDesc->incif;
  in->gpuInfo->h_visInfo->nstok  = inDesc->inaxes[inDesc->jlocs];
  in->gpuInfo->h_visInfo->stokb  = 0;  /* ??? */
  in->gpuInfo->h_visInfo->stoke  = in->gpuInfo->h_visInfo->stokb + 
    MIN (2, in->gpuInfo->h_visInfo->nstok) - 1;   /* No more than 2 */
  in->gpuInfo->h_visInfo->incs   = inDesc->incs;
  in->gpuInfo->h_visInfo->kincif = kincif; /* IF increment in d_freqScale */
  in->gpuInfo->h_visInfo->kincf  = kincf;  /* Freq increment in d_freqScale */
  if (inDesc->crval[inDesc->jlocs]>0.0)       in->gpuInfo->h_visInfo->stype  = 0;
  else if (inDesc->crval[inDesc->jlocs]<-4.0) in->gpuInfo->h_visInfo->stype  = 2;
  else                                        in->gpuInfo->h_visInfo->stype  = 1;
  /* ratio of uv reference freq to image */
  in->gpuInfo->h_visInfo->h_freqRat = NULL;
  /* depends on type */
  if (ObitSkyModelVMBeamMFIsA(skyModel)) {  /* One per coarse channel */
    nSpec     = ((ObitSkyModelVMBeamMF*)skyModel)->nSpec;
    specFreq = ((ObitSkyModelVMBeamMF*)skyModel)->specFreq;
    memsize = nSpec * sizeof(float);
    in->gpuInfo->h_visInfo->h_freqRat = (float*)ObitCUDAUtilAllocHost(memsize);
    for (i=0; i<nSpec; i++) 
      in->gpuInfo->h_visInfo->h_freqRat[i] =  (float)(uvdata->myDesc->freq/specFreq[i]);
    /* to GPU */
    in->gpuInfo->h_visInfo->d_freqRat = (float*)ObitCUDAUtilAllocGPU(memsize);
    ObitCUDAUtilHost2GPU ((float*)in->gpuInfo->h_visInfo->d_freqRat, 
			  (float*)in->gpuInfo->h_visInfo->h_freqRat, 
			  memsize, NULL);
 } else {  /* Single spectrum */
    memsize = sizeof(float);
    in->gpuInfo->h_visInfo->h_freqRat = (float*)ObitCUDAUtilAllocHost(memsize);
    in->gpuInfo->h_visInfo->h_freqRat[0] = 1.0;   /* in case */
    if (((ObitSkyModel*)skyModel)->mosaic!=NULL) {
      img = ((ObitSkyModel*)skyModel)->mosaic->images[0];
      in->gpuInfo->h_visInfo->h_freqRat[0] =
	(float)(uvdata->myDesc->freq/img->myDesc->crval[img->myDesc->jlocf]);
    /* to GPU */
    in->gpuInfo->h_visInfo->d_freqRat = (float*)ObitCUDAUtilAllocGPU(memsize);
    ObitCUDAUtilHost2GPU ((float*)in->gpuInfo->h_visInfo->d_freqRat, 
			  (float*)in->gpuInfo->h_visInfo->h_freqRat, 
			  memsize, NULL);
    }
  }
  /* visInfo to GPU */
  memsize = sizeof(GPUVisInfo);
  ObitCUDAUtilHost2GPU ((float*)in->gpuInfo->d_visInfo, (float*)in->gpuInfo->h_visInfo, 
			memsize, NULL);

  /****    Model info    - much in DFTSetMod ****/
  nmodel = in->gpuInfo->nmodel;
  memsize = sizeof(GPUModelInfo);
  in->gpuInfo->d_modelInfo = (GPUBeamModelInfo*)ObitCUDAUtilAllocGPU(memsize);
  in->gpuInfo->h_modelInfo = (GPUBeamModelInfo*)ObitCUDAUtilAllocHost(memsize);
  in->gpuInfo->h_modelInfo->nmodel = nmodel;
  in->gpuInfo->h_modelInfo->nterm  = ((ObitSkyModel*)skyModel)->nSpecTerm;
  in->gpuInfo->h_modelInfo->nSpec  = ((ObitSkyModelVMBeamMF*)skyModel)->nSpec;
  in->gpuInfo->h_modelInfo->size   = in->gpuInfo->modelSize;
  /* Model Type (set for real in DFTSetMod) */
  in->gpuInfo->h_modelInfo->type    = GPUModTypePoint;
  in->gpuInfo->h_modelInfo->d_model = in->gpuInfo->d_model;
  in->gpuInfo->h_modelInfo->stokType= 1; /* Stokes I */

  /****     Antenna/Beam info    ****/
  /* gpuInfo->antInfo host and device */
  beamDesc = skyModel->RXBeam[0]->ImgDesc;  /* ImageDesc of one of the beam images */
  nant = skyModel->AntList[0]->number;      /* Doesn't handle subarrays */
  nbchan = beamDesc->inaxes[beamDesc->jlocf];  /* Number of beam channels */
  memsize = sizeof(GPUBeamAntInfo);
  in->gpuInfo->h_antInfo = (GPUBeamAntInfo*)ObitCUDAUtilAllocHost(memsize);
  in->gpuInfo->d_antInfo = (GPUBeamAntInfo*)ObitCUDAUtilAllocGPU(memsize);
  in->gpuInfo->h_antInfo->numAnt     = nant;
  in->gpuInfo->h_antInfo->numAntType = skyModel->numAntType;
  in->gpuInfo->h_antInfo->numJChan   = skyModel->numPlane[0]; /* Better all be the same */
  in->gpuInfo->h_antInfo->ilocb      = inDesc->ilocb;
  in->gpuInfo->h_antInfo->iloca1     = inDesc->iloca1;
  in->gpuInfo->h_antInfo->iloca2     = inDesc->iloca2;
  in->gpuInfo->h_antInfo->ilocsa     = inDesc->ilocsa;
  in->gpuInfo->h_antInfo->s2p        = 0.0;
  in->gpuInfo->h_antInfo->c2p        = 1.0;
  in->gpuInfo->h_antInfo->refantSize = skyModel->BeamShape->antSize;
  /* Factor for gain antSize*fudge*DG2RAD/VELIGHT, fudge=1.052 
     from ObitBeamShape.c: GetCos2Beam */
  in->gpuInfo->h_antInfo->refSizeFact =
    (float)(skyModel->BeamShape->antSize*1.052*0.017453292519943295/2.997924562e8);
  /* Arrays  Antenna type per antenna */
  memsize = nant*sizeof(int);
  in->gpuInfo->h_antInfo->h_AntType    = (int*)ObitCUDAUtilAllocHost(memsize);
  in->gpuInfo->h_antInfo->d_AntType    = (int*)ObitCUDAUtilAllocGPU(memsize);
  for (i=0;i<nant; i++) in->gpuInfo->h_antInfo->h_AntType[i] = (int)skyModel->AntType[i];
  ObitCUDAUtilHost2GPU ((float*)in->gpuInfo->h_antInfo->d_AntType, 
			(float*)in->gpuInfo->h_antInfo->h_AntType, 
			memsize, NULL);
  ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_antInfo->h_AntType);  /* Cleanup */
  in->gpuInfo->h_antInfo->h_AntType = NULL;

  /* Frequency per beam channel *h_BeamFreq, *d_BeamFreq; */
  memsize = nbchan*sizeof(float);
  in->gpuInfo->h_antInfo->h_BeamFreq  = (float*)ObitCUDAUtilAllocHost(memsize);
  in->gpuInfo->h_antInfo->d_BeamFreq  = (float*)ObitCUDAUtilAllocGPU(memsize);
  /* Number of vis channels per beam channel get from one of the beam images */
  beamDesc = skyModel->RXBeam[0]->ImgDesc;
  f_0 = beamDesc->crval[beamDesc->jlocf]; delt_f = beamDesc->cdelt[beamDesc->jlocf];
  crp_f = beamDesc->crpix[beamDesc->jlocf];
  for (i=0; i<nbchan; i++) 
    in->gpuInfo->h_antInfo->h_BeamFreq[i] = (float)(f_0+(i+crp_f-1)*delt_f);
  ObitCUDAUtilHost2GPU ((float*)in->gpuInfo->h_antInfo->d_BeamFreq,
			(float*)in->gpuInfo->h_antInfo->h_BeamFreq, 
			memsize, NULL);
  ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_antInfo->h_BeamFreq);  /* Cleanup */
  in->gpuInfo->h_antInfo->h_BeamFreq = NULL;

  /* beam channel per vis channel */
  memsize = skyModel->numUVChann*sizeof(int);
  in->gpuInfo->h_antInfo->h_visJChann  = (int*)ObitCUDAUtilAllocHost(memsize);
  in->gpuInfo->h_antInfo->d_visJChann  = (int*)ObitCUDAUtilAllocGPU(memsize);
  for (i=0; i<skyModel->numUVChann; i++)
    in->gpuInfo->h_antInfo->h_visJChann[i] = (int)skyModel->FreqPlane[i];
  ObitCUDAUtilHost2GPU ((float*)in->gpuInfo->h_antInfo->d_visJChann, 
			(float*)in->gpuInfo->h_antInfo->h_visJChann, 
			memsize, NULL);
  ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_antInfo->h_visJChann);  /* Cleanup */
  in->gpuInfo->h_antInfo->h_visJChann = NULL;

   /****    Allocate host/device vis data arrays    ****/
  in->gpuInfo->d_data  = g_malloc0(in->gpuInfo->nstream*sizeof(float*));
  memsize = ((in->gpuInfo->nvis+5) * in->gpuInfo->lenvis) * 
    sizeof(float)/in->gpuInfo->nstream;  // size of device vis buffer
  for (i=0; i<in->gpuInfo->nstream; i++) {
    in->gpuInfo->d_data[i]  = ObitCUDAUtilAllocGPU(memsize);
  }

  /****    Allocate stream stuff     ****/
  in->gpuInfo->stream    = g_malloc0(in->gpuInfo->nstream*sizeof(void*));
  in->gpuInfo->cycleDone = g_malloc0(in->gpuInfo->nstream*sizeof(void*));
  for (i=0; i<in->gpuInfo->nstream; ++i)
    {
      in->gpuInfo->stream[i]    = ObitCUDAStreamCreate();
      in->gpuInfo->cycleDone[i] = ObitCUDAEventCreate();
      ObitCUDAEventRecord (in->gpuInfo->cycleDone[i], in->gpuInfo->stream[i]);
    }

  /* CALL CUDA for setup in GPU */
  ObitCUDABeamModelDFTInit (in->gpuInfo);
#endif /* HAVE_GPU */
  } /* end  ObitGPUBeamModelDFTInit */


/**
 * Setup model for ObitGPUBeamModel 
 * Repackages model to include spectral index for each subband flux density.
 * \param in      Input object
 * \param skyModel Skymodel to implement on GPU, possibly a derived type
 *                 disguised as an Obit to keep from confusing the compiler
 * \param model   model components, w/ subband flux densities but w/o SI
 * \param err     Obit error stack object.
 */
void ObitGPUBeamModelDFTSetMod (ObitGPUBeamModel *in, Obit *skyModel,
			       ObitFArray *model, ObitErr *err)
{
#if HAVE_GPU==1  /* GPU? Real versions */
  int i, inc, high, memsize;
  ObitFArray *tmpmodel=NULL;
  ObitSpectrumMF *spEval=NULL;
  olong off, naxis[2], pos[2], imod, ispec, nSpec, j, ncopy, nterm, mterm, nplane, plnoff, *specIndex;
  gboolean isTSpec=FALSE, isGauss=FALSE;
  ofloat Threshold=0.0, *newmodel=NULL, *oldmodel=NULL, *SI=NULL, *PlnRMS=NULL;
  ofloat sumFlux=0.0;
  ObitSkyModel *skymod = (ObitSkyModel*)skyModel;
  ObitImageDesc *imDesc;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];

  /* Info from loaded skyModel */
  ObitInfoListGetTest (skymod->info, "Threshold",  &type, dim, &Threshold);
  nplane = skymod->mosaic->images[0]->myDesc->inaxes[2];  /* No. planes in image */
  in->gpuInfo->h_modelInfo->nterm  = skymod->nTerm;
  nterm = in->gpuInfo->h_modelInfo->nterm;  /* number of spectral terms in input model */
  in->gpuInfo->h_modelInfo->nSpec  = ((ObitSkyModelVMBeamMF*)skyModel)->nSpec;
  /* array of coarse frequency (subband) index per vis frequency channel */
  if (in->gpuInfo->h_modelInfo->h_specIndex) ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_modelInfo->h_specIndex);
  in->gpuInfo->h_modelInfo->h_specIndex = NULL;
   /* depends on type */
   if (ObitSkyModelVMBeamMFIsA(skyModel)) {  /* One per channel/IF  */
     nSpec    = in->gpuInfo->h_visInfo->nchan;
     specIndex = ((ObitSkyModelVMBeamMF*)skyModel)->specIndex;
     memsize = nSpec * sizeof(int);
     in->gpuInfo->h_modelInfo->h_specIndex = (int*)ObitCUDAUtilAllocHost(memsize);
     for (i=0; i<nSpec; i++) in->gpuInfo->h_modelInfo->h_specIndex[i] = (int)specIndex[i];
     /* to GPU */
     in->gpuInfo->h_modelInfo->d_specIndex = (int*)ObitCUDAUtilAllocGPU(memsize);
     ObitCUDAUtilHost2GPU ((float*)in->gpuInfo->h_modelInfo->d_specIndex, 
			   (float*)in->gpuInfo->h_modelInfo->h_specIndex, 
 			  memsize, NULL);
   } else in->gpuInfo->h_modelInfo->d_specIndex = NULL;
   /* Stokes type */
   imDesc = skymod->mosaic->images[0]->myDesc;
   in->gpuInfo->h_modelInfo->stokType=MIN (4, MAX(1,imDesc->crval[imDesc->jlocs]));

  
  /* Does this model require special handling?*/
  isTSpec = (((ObitSkyModel*)skyModel)->modType==OBIT_SkyModel_PointModTSpec) ||
    (((ObitSkyModel*)skyModel)->modType==OBIT_SkyModel_GaussModTSpec);
  isGauss = (((ObitSkyModel*)skyModel)->modType==OBIT_SkyModel_GaussModTSpec);

  /* Find last non zero entry */
  inc = model->naxis[0];
  if (isGauss) off = 10; /* Gaussian */
  else         off = 7;  /* Point */
  high = 1; sumFlux = 0.0;
  for (i=0; i<model->naxis[1]; i++)
    if (model->array[off+i*inc]!=0.0) {high = i+1; sumFlux += model->array[off+i*inc];}

  Obit_log_error(err, OBIT_InfoErr, "SkyModel by GPU with %d components, sum model %f", high,sumFlux);
  ObitErrLog(err);

  /* (re)build device model array  */
  if ((in->gpuInfo->nmodel!=high) || (in->gpuInfo->modelSize!=model->naxis[0])) {
    in->gpuInfo->nmodel = high;
    in->gpuInfo->h_modelInfo->nmodel = high;
    nSpec = ((ObitSkyModelVMBeamMF*)skyModel)->nSpec;
    if (isTSpec)
      in->gpuInfo->modelSize = model->naxis[0]-nterm+nSpec;  /* Add spectral indices */
    else
      in->gpuInfo->modelSize = model->naxis[0];  /* Shouldn't be a spectral fit */
    in->gpuInfo->h_modelInfo->size   = in->gpuInfo->modelSize;
    if (in->gpuInfo->d_model) ObitCUDAUtilFreeGPU(in->gpuInfo->d_model);
    /* new Device Model array */
    memsize = in->gpuInfo->nmodel * in->gpuInfo->modelSize * sizeof(float);
    in->gpuInfo->d_model  = ObitCUDAUtilAllocGPU(memsize);
    in->gpuInfo->h_modelInfo->nterm  =  ((ObitSkyModel*)skyModel)->nSpecTerm;
    in->gpuInfo->h_modelInfo->d_model  = in->gpuInfo->d_model;

    /* Temporary working host copy for TSpec models adding spectral indices */
    /* Spectrum from LoadComps has 7 terms giving geometry followed by a 
       nterm spectral fit. The spectrum is refitted  and the original fit discarded. */
    if (isTSpec) {
      /* Are the plane RMSes on in->skyModel? */
      ObitInfoListGetP(((ObitSkyModel*)skyModel)->info, "PlnRMS",
		       &type, dim, (gpointer*)&PlnRMS);
      /* Spectrum evaluator */
      mterm = MIN(4, nSpec-1);
      spEval = newObitSpectrumMFCreate ("Sp eval", nSpec,
					((ObitSkyModelVMBeamMF*)skyModel)->specFreq,
					((ObitSkyModelVMBeamMF*)skyModel)->refFreq, mterm);
      /* Set sigmas for weighting if available */
      plnoff = nplane-nSpec;  /* Offset in PlnRMS for subbands */
      if (PlnRMS)  ObitSpectrumMFSetSigma (spEval, &PlnRMS[plnoff]);
      /* Temporary array */
      naxis[0] = in->gpuInfo->modelSize;  naxis[1] = high;
      tmpmodel = ObitFArrayCreate ("Expand copy", 2, naxis);
      /* Copy adding spectral indices */
      SI = g_malloc0(nSpec*sizeof(ofloat));  /* Work array */
      pos[0] = 0; pos[1] = 0; sumFlux = 0.0;
      if (isGauss) ncopy = 10; /* Gaussian */
      else         ncopy = 7;  /* Point */
      for (imod=0; imod<high; imod++) {
  	pos[1] = imod;
	oldmodel = ObitFArrayIndex(model, pos);     /* Pointer in array */
	newmodel = ObitFArrayIndex(tmpmodel, pos);  /* Pointer in array */
	/* Missing fits will have flux=avg, spectral terms 0 */
	sumFlux += oldmodel[ncopy];  /* fooey, sum model flux density */
	for (j=0; j<ncopy; j++) newmodel[j] = oldmodel[j]; /* Copy initial values */
	/* Fit spectral indices wrong nterm, 2 in image */
	ObitSpectrumMFSI (spEval, &oldmodel[ncopy+nterm], SI, err);
	/* Copy tabulated spectra adding spectral indices */
	for (ispec=0; ispec<nSpec; ispec++){
	  newmodel[ncopy+ispec*2]   = oldmodel[ncopy+nterm+ispec]; /* Flux density */
	  if (fabs(SI[ispec])<5.0)
	    newmodel[ncopy+ispec*2+1] = SI[ispec];     /* Spectral index, 0 if out of range */
	  else newmodel[ncopy+ispec*2+1] = 0.0;
	  /* 0 SI if not IPol */
	  if (in->gpuInfo->h_modelInfo->stokType!=1) newmodel[ncopy+ispec*2+1] = 0.0;
	} /* end loop over tabulates spectrum */
      } /* end loop over models */
      if (SI) {g_free(SI); SI=NULL;} /* Cleanup */
      if (spEval) spEval =  ObitSpectrumMFUnref(spEval);
    } else { /* No repackaging needed */
      tmpmodel = model;
    }
    /* Model Type */
    if (((ObitSkyModel*)skyModel)->modType==OBIT_SkyModel_PointMod)
      in->gpuInfo->h_modelInfo->type = GPUModTypePoint;
    else if (((ObitSkyModel*)skyModel)->modType==OBIT_SkyModel_PointModSpec)  
      in->gpuInfo->h_modelInfo->type = GPUModTypePointSpec;
    else if (((ObitSkyModel*)skyModel)->modType==OBIT_SkyModel_PointModTSpec) 
      in->gpuInfo->h_modelInfo->type = GPUModTypePointTSpec;
    else if (((ObitSkyModel*)skyModel)->modType==OBIT_SkyModel_GaussMod)
      in->gpuInfo->h_modelInfo->type = GPUModTypeGauss;
    else if (((ObitSkyModel*)skyModel)->modType==OBIT_SkyModel_GaussModSpec)  
      in->gpuInfo->h_modelInfo->type = GPUModTypeGaussSpec;
    else if (((ObitSkyModel*)skyModel)->modType==OBIT_SkyModel_GaussModTSpec) 
      in->gpuInfo->h_modelInfo->type = GPUModTypeGaussTSpec;
    else g_error ("Model type not supported in GPU %d",((ObitSkyModel*)skyModel)->modType);  /* Bother */
    /* Operation type */
    in->gpuInfo->h_modelInfo->opType = GPUOpTypeSub;
    if (((ObitSkyModel*)skyModel)->doDivide)  in->gpuInfo->h_modelInfo->opType = GPUOpTypeDiv;
    if (((ObitSkyModel*)skyModel)->doReplace) in->gpuInfo->h_modelInfo->opType = GPUOpTypeRepl;
    /* Factor/doFlip */
    in->gpuInfo->h_modelInfo->stokFact[0] = ((ObitSkyModel*)skyModel)->factor;
    in->gpuInfo->h_modelInfo->stokFact[1] = ((ObitSkyModel*)skyModel)->factor;
    in->gpuInfo->h_modelInfo->stokFact[2] = ((ObitSkyModel*)skyModel)->factor;
    in->gpuInfo->h_modelInfo->stokFact[3] = ((ObitSkyModel*)skyModel)->factor;
    if (((ObitSkyModel*)skyModel)->doFlip) {
           in->gpuInfo->h_modelInfo->doSwap = 1;
    } else in->gpuInfo->h_modelInfo->doSwap = 0;
  } /* End rebuild model */
  /* Copy model to device NB: host memory NOT locked */
  memsize = in->gpuInfo->nmodel * in->gpuInfo->modelSize * sizeof(float);
  ObitCUDAUtilHost2GPU(in->gpuInfo->d_model, tmpmodel->array, memsize, NULL);
  if (isTSpec) tmpmodel = ObitFArrayUnref(tmpmodel); /* Free Temporary copy */

  /* Copy Model info to device */
  in->gpuInfo->h_modelInfo->nmodel = in->gpuInfo->nmodel;
  in->gpuInfo->h_modelInfo->size   = in->gpuInfo->modelSize;
  memsize = sizeof(GPUModelInfo);
  ObitCUDAUtilHost2GPU((float*)in->gpuInfo->d_modelInfo, 
		       (float*)in->gpuInfo->h_modelInfo, 
		       memsize, NULL);

   /* Allocate device Jones matrices */
  memsize = in->gpuInfo->h_antInfo->numJChan*in->gpuInfo->h_antInfo->numAntType*in->gpuInfo->h_modelInfo->nmodel*sizeof(cudaMatx);
  in->gpuInfo->h_antInfo->d_Jones      = (cudaMatx*)ObitCUDAUtilAllocGPU(memsize);
  /* d_BeamInterp managed in ObitGPUBeamModelDFTSetBeam */
  in->gpuInfo->h_antInfo->d_BeamInterp = NULL;
  /* Copy antInfo to device */
  memsize = sizeof(GPUBeamAntInfo);
  ObitCUDAUtilHost2GPU ((float*)in->gpuInfo->d_antInfo, 
			(float*)in->gpuInfo->h_antInfo, 
			memsize, NULL);

  /* Threshold reached? */
  if (fabs(sumFlux)<Threshold) {
    in->gpuInfo->doUnit = 1;
    Obit_log_error(err, OBIT_InfoErr, "Beam corrections turned off, model %g < Threshold %g",
		   sumFlux, Threshold);
  } else {
    Obit_log_error(err, OBIT_InfoErr, "Sum of model flux density %g with Beam Cor", sumFlux);
  }
    ObitErrLog(err);

  /* CALL CUDA setup*/
  ObitCUDABeamModelDFTSetMod (in->gpuInfo);
  return;
#endif /* HAVE_GPU */
} /* end ObitGPUBeamModelDFTSetMod */

/**
 * Setup./copy Beam interpolators to GPU for ObitGPUBeam
 * Ordering PP r, i, PQ r, i, QU r, i, PP r, i, antenna type (fastest to slowest).
 * P=X/R, Q=Y/L
 * \param in       Input object
 * \param skyModel SkyModel being evaluated as Obit to keep from confusing compiler
 * \param err     Obit error stack object.
 */
void ObitGPUBeamModelDFTSetBeam (ObitGPUBeamModel *in, Obit *SkyModel, ObitErr *err)
{
#if HAVE_GPU==1  /* GPU? Real versions */
  int hwidth = 2;
  olong itype, indx;
  olong numAntType=((ObitSkyModelVMBeam*)SkyModel)->numAntType;  /* Number of antenna types */
  int memsize;
  ObitSkyModelVMBeam *skyModel = (ObitSkyModelVMBeam*)SkyModel;
  CUDABeamInterp *h_fInt=NULL, **h_fIntArr = NULL;

  /* Array of FInterpolators, numAntType*8 */
  memsize = numAntType*8 * sizeof (CUDABeamInterp*);
  h_fIntArr = (CUDABeamInterp**)ObitCUDAUtilAllocHost(memsize);
  in->gpuInfo->h_antInfo->d_BeamInterp = (CUDABeamInterp**)ObitCUDAUtilAllocGPU(memsize);
  /* Copy antInfo to device to update d_BeamInterp */
  memsize = sizeof(GPUBeamAntInfo);
  ObitCUDAUtilHost2GPU ((float*)in->gpuInfo->d_antInfo, 
			(float*)in->gpuInfo->h_antInfo, 
			memsize, NULL);

  /* Copy array of beam interpolators to device */
  for (itype=0; itype<numAntType; itype++) {
    indx = itype*8;
    h_fInt = ObitGPUCUDABeamInterpCreate (skyModel->RXBeam[itype]->ImgPixels, skyModel->RXBeam[itype]->ImgDesc, hwidth);
    h_fIntArr[indx+0] = (CUDABeamInterp*)h_fInt->d_BeamInterp; ObitCUDAUtilFreeHost((float*)h_fInt->h_myArray); h_fInt->h_myArray=NULL;
    h_fInt = ObitGPUCUDABeamInterpCreate (skyModel->RXBeamIm[itype]->ImgPixels, skyModel->RXBeam[itype]->ImgDesc, hwidth);
    h_fIntArr[indx+1] =(CUDABeamInterp*) h_fInt->d_BeamInterp; ObitCUDAUtilFreeHost((float*)h_fInt->h_myArray);  h_fInt->h_myArray=NULL;
    h_fInt = ObitGPUCUDABeamInterpCreate (skyModel->LYBeam[itype]->ImgPixels, skyModel->LYBeam[itype]->ImgDesc, hwidth);
    h_fIntArr[indx+2] = (CUDABeamInterp*)h_fInt->d_BeamInterp; ObitCUDAUtilFreeHost((float*)h_fInt->h_myArray);  h_fInt->h_myArray=NULL;
    h_fInt = ObitGPUCUDABeamInterpCreate (skyModel->LYBeamIm[itype]->ImgPixels, skyModel->LYBeam[itype]->ImgDesc, hwidth);
    h_fIntArr[indx+3] = (CUDABeamInterp*)h_fInt->d_BeamInterp; ObitCUDAUtilFreeHost((float*)h_fInt->h_myArray); h_fInt->h_myArray=NULL;
    h_fInt = ObitGPUCUDABeamInterpCreate (skyModel->RLBeam[itype]->ImgPixels, skyModel->RLBeam[itype]->ImgDesc, hwidth);
    h_fIntArr[indx+4] = (CUDABeamInterp*)h_fInt->d_BeamInterp; ObitCUDAUtilFreeHost((float*)h_fInt->h_myArray); h_fInt->h_myArray=NULL;
    h_fInt = ObitGPUCUDABeamInterpCreate (skyModel->RLBeamIm[itype]->ImgPixels, skyModel->RLBeam[itype]->ImgDesc, hwidth);
    h_fIntArr[indx+5] =(CUDABeamInterp*) h_fInt->d_BeamInterp; ObitCUDAUtilFreeHost((float*)h_fInt->h_myArray); h_fInt->h_myArray=NULL;
    h_fInt = ObitGPUCUDABeamInterpCreate (skyModel->LRBeam[itype]->ImgPixels, skyModel->LRBeam[itype]->ImgDesc, hwidth);
    h_fIntArr[indx+6] = (CUDABeamInterp*)h_fInt->d_BeamInterp; ObitCUDAUtilFreeHost((float*)h_fInt->h_myArray); h_fInt->h_myArray=NULL;
    h_fInt = ObitGPUCUDABeamInterpCreate (skyModel->LRBeamIm[itype]->ImgPixels, skyModel->LRBeam[itype]->ImgDesc, hwidth);
    h_fIntArr[indx+7] = (CUDABeamInterp*)h_fInt->d_BeamInterp; ObitCUDAUtilFreeHost((float*)h_fInt->h_myArray); h_fInt->h_myArray=NULL;
    /* This CUDABeamInterpCreate is in ObitGPUBeamInterp.c
       pixel values are copied to the GPU returning the device address */
  } /* end loop over types */
  memsize = numAntType*8*sizeof(CUDABeamInterp*);
  ObitCUDAUtilHost2GPU ((float*)in->gpuInfo->h_antInfo->d_BeamInterp, (float*)h_fIntArr, memsize, NULL);
  /* Cleanup */
  //NO for (i=0; i<numAntType*8; i++) ObitCUDAUtilFreeHost((float*)h_fIntArr[i]);
  ObitCUDAUtilFreeHost((float*)h_fIntArr); h_fIntArr = NULL; 
#endif /* HAVE_GPU */
} /* end ObitGPUBeamModelDFTSetBeam */

/**
 * Calculate an ObitGPUBeamModel - process a buffer load of visibility data
 * Processing done in blocks of parallactic angle range with Jones matrices
 * computed at the beginning of each block.
 * Currently only DFT point supported
 * \param in      Iniput object
 * \param uvdata  UV data being operated on
 * \param err     Obit error stack object.
 */
void ObitGPUBeamModelDFTCalc (ObitGPUBeamModel *in, ObitUV *uvdata, ObitErr *err) {
/**************************************************************************************************
  This gets passed an arbitrary block of data, needs to be limited in Parallactic angle range
    could probably divide the data here, pass the data, remake the Jones matrices
    How does this affect using streams to overlap IO and computation?
    Calculating the Jones matrices at the beginning of a block doesn't have much I/O overhead as
    the beams will already be in the device How many streams???
    can calculate when the next parallactic angle update is needed without reading the data first.
**************************************************************************************************/
#if HAVE_GPU==1  /* GPU? Real versions */
  olong iloct, offset, nvis=uvdata->myDesc->numVisBuff, lenvis=uvdata->myDesc->lrec;
    size_t memsize;
    gboolean newScan=TRUE, isDone = FALSE;
    int doJones=1, doUnit=0;
    ObitUVDesc *desc = uvdata->myDesc;
    ofloat curPA,  bvTime, evTime, bTime, eTime, tTime;
    ofloat delta_ParAng=in->gpuInfo->delta_ParAng*DG2RAD;
    long ivis, visRange[2];
    ObitAntennaList *antList = in->antList;
    ObitSource *curSource    = in->curSource;

    // Progress report if in->prtLv>=5
    if (err->prtLv>=5){
      gchar timeStr[24];
      day2dhms(uvdata->buffer[uvdata->myDesc->iloct], timeStr); /* Human readable time */
      int numComp = in->gpuInfo->nmodel;
      Obit_log_error(err, OBIT_InfoErr,"Time %s nvis %6d, ncomp %6d", timeStr, nvis,numComp);
      ObitErrLog(err);
    }

    /* Copy data to locked buffer */
    memsize = nvis * lenvis * sizeof(ofloat);
    memcpy (in->gpuInfo->h_data, uvdata->buffer, memsize);
    
    /* In case all in the same Parallactic angle range */
    visRange[0] = 0; visRange[1] = nvis-1; doJones = 1; doUnit = 0;;
  
    /* Get start and end times in the buffer */
    iloct = desc->iloct;
    bTime = uvdata->buffer[iloct];
    offset = (nvis-1) * lenvis;
    eTime = uvdata->buffer[offset+iloct];
    newScan = bTime>in->gpuInfo->timeRange[1]; /* Is this a new scan */

    /* Need Parallactic angle and validity range*/
    curPA = ObitSkyModelVMBeamPAUpdate (1, 0, bTime, delta_ParAng,
					antList, curSource, &bvTime, &evTime);
    /* Does curPA need updating in GPU? */
    doJones = newScan || fabs(curPA-in->gpuInfo->curParAng)>delta_ParAng;  /* same parallactic range as end of last? */
    
    /* Loop dividing by parallactic angle range */
    while (!isDone) {
      /* May take several passes dividing the data by parallactic angle range. */
      /* If evTime>eTime then done */
      if (evTime>=eTime) {
	isDone=TRUE; visRange[1] = nvis-1;
      } else { /* have to divide */
	isDone = FALSE;
	/* Find first vis beyond evTime */
	for (ivis=visRange[0]; ivis<nvis; ivis++) {
	  if (uvdata->buffer[ivis*lenvis+iloct]>evTime) break;
	}
	visRange[1] = MIN(ivis, nvis-1);
      } /* end have to divide */
      
      /* CALL CUDA for this batch */
      /* Are unit matrices rather than Beam jones matrices been requested?
	 if so, have they been initialised?*/
      if (in->gpuInfo->doUnit) {
	doJones = 0; doUnit=1;
	if (in->gpuInfo->haveUnit) doUnit = 0;  /* Don't need again */
      }
      if (newScan) visRange[0]=0;
      ObitCUDABeamModelDFTCalc (in->gpuInfo, doJones, doUnit, curPA, visRange);
      newScan = FALSE;
      if (in->gpuInfo->doUnit) in->gpuInfo->haveUnit = 1;  /* now done */
      /* Save this parallactic angle, timerange */
      in->gpuInfo->curParAng = curPA;
      in->gpuInfo->timeRange[0] = bvTime;
      in->gpuInfo->timeRange[1] = evTime;
      /* Next batch, or Done */
      visRange[0] = visRange[1]+1;
      if (visRange[0]>=nvis) break;  /* Done? */
      tTime = uvdata->buffer[visRange[0]*lenvis+iloct];  /* may be gaps */
      curPA = ObitSkyModelVMBeamPAUpdate (1, 0, tTime, delta_ParAng,
					  antList, curSource, &bvTime, &evTime);
      doJones = 1;  /* New set of Jones matrices */
    } /* end not done */

   
    /* Copy data back to IO buffer */
    memsize = in->gpuInfo->nvis * in->gpuInfo->lenvis * sizeof(ofloat);
    memcpy (uvdata->buffer, in->gpuInfo->h_data, memsize);
    return;
#endif /* HAVE_GPU */
 } /* end ObitGPUBeamModelDFTCalc */

/**
 * Shutdown  ObitGPUBeamModel 
 * Currently only DFT point supported
 * \param in      Iniput object
 * \param uvdata  UV data being operated on
 * \param err     Obit error stack object.
 */
void ObitGPUBeamModelDFTShutdown (ObitGPUBeamModel *in, ObitUV *uvdata, ObitErr *err)
{
#if HAVE_GPU==1  /* GPU? Real versions */
  int i;
  gint32 dim[MAXINFOELEMDIM];
  /*fprintf(stderr,"DEBUG in ObitGPUBeamModelDFTShutdown\n");*/
  /* Something to do? */
  if (in==NULL) return;
  // Free resources
  if (in->gpuInfo) {
    if (in->gpuInfo->h_data) 
      {ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_data);    in->gpuInfo->h_data=NULL;}
    if (in->gpuInfo->h_visInfo) {
      if (in->gpuInfo->h_visInfo->d_freqRat) {
	ObitCUDAUtilFreeGPU((float*)in->gpuInfo->h_visInfo->d_freqRat);
	in->gpuInfo->h_visInfo->d_freqRat = NULL;
      }
      if (in->gpuInfo->h_visInfo->h_freqRat) {
	ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_visInfo->h_freqRat);
	in->gpuInfo->h_visInfo->h_freqRat = NULL;
      }
      ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_visInfo); in->gpuInfo->h_visInfo=NULL;
    }
    if (in->gpuInfo->d_visInfo)
      {ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_visInfo);  in->gpuInfo->d_visInfo=NULL;}
    if (in->gpuInfo->d_modelInfo)
      {ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_modelInfo);in->gpuInfo->d_modelInfo=NULL;}
    if (in->gpuInfo->d_model)
      {ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_model); in->gpuInfo->d_model=NULL;}
    for (i=0; i<in->gpuInfo->nstream; ++i)    {
      if (in->gpuInfo->stream)    ObitCUDAStreamDestroy((in->gpuInfo->stream)[i]);
      if (in->gpuInfo->cycleDone) ObitCUDAEventDestroy((in->gpuInfo->cycleDone)[i]);
      if (in->gpuInfo->d_data)    
	{ObitCUDAUtilFreeGPU(in->gpuInfo->d_data[i]); in->gpuInfo->d_data[i]=NULL;}
    }
    if (in->gpuInfo->stream)    {g_free(in->gpuInfo->stream);   in->gpuInfo->stream=NULL;}
    if (in->gpuInfo->cycleDone) {g_free(in->gpuInfo->cycleDone);in->gpuInfo->cycleDone=NULL;}
    if (in->gpuInfo->d_data)    {g_free(in->gpuInfo->d_data); in->gpuInfo->d_data=NULL;}
  
    /* Replace old nVisPIO */
    dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
    ObitInfoListAlwaysPut (uvdata->info, "nVisPIO", OBIT_long, dim,  
			   &in->gpuInfo->oldnVisPIO);

    /* CALL CUDA shutdown */
    ObitCUDABeamModelDFTShutdown (in->gpuInfo);
    in->gpuInfo = NULL;
  }
#endif /* HAVE_GPU */
} /*end ObitGPUBeamModelDFTShutdown */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitGPUBeamModelClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGPUBeamModelClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
} /* end ObitGPUBeamModelClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGPUBeamModelClassInfoDefFn (gpointer inClass)
{
  ObitGPUBeamModelClassInfo *theClass = (ObitGPUBeamModelClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGPUBeamModelClassInit;
  theClass->newObit       = (newObitFP)newObitGPUBeamModel;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGPUBeamModelClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGPUBeamModelGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitGPUBeamModelCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitGPUBeamModelClear;
  theClass->ObitInit      = (ObitInitFP)ObitGPUBeamModelInit;
#if HAVE_GPU==1  /* GPU? Real versions */
  theClass->ObitGPUBeamModelCreate = (ObitGPUBeamModelCreateFP)ObitGPUBeamModelCreate;
  theClass->ObitGPUBeamModelDFTInit = (ObitGPUBeamModelDFTInitFP)ObitGPUBeamModelDFTInit;
  theClass->ObitGPUBeamModelDFTSetMod = (ObitGPUBeamModelDFTSetModFP)ObitGPUBeamModelDFTSetMod;
  theClass->ObitGPUBeamModelDFTSetBeam= (ObitGPUBeamModelDFTSetBeamFP)ObitGPUBeamModelDFTSetBeam;
  theClass->ObitGPUBeamModelDFTCalc = (ObitGPUBeamModelDFTCalcFP)ObitGPUBeamModelDFTCalc;
  theClass->ObitGPUBeamModelDFTShutdown = (ObitGPUBeamModelDFTShutdownFP)ObitGPUBeamModelDFTShutdown;
#endif /* HAVE_GPU */
} /* end ObitGPUBeamModelClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitGPUBeamModelInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGPUBeamModel *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  /* gpuInfo - host only */
  in->gpuInfo   = NULL;
  in->antList   = NULL;
  in->curSource = NULL;
} /* end ObitGPUBeamModelInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitGPUBeamModel* cast to an Obit*.
 */
void ObitGPUBeamModelClear (gpointer inn)
{
#if HAVE_GPU==1  /* GPU? Real versions */
  ObitClassInfo *ParentClass;
  ObitGPUBeamModel *in = inn;
  olong i, j, indx, nAntType=0;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (!in) return;

  if (in->gpuInfo->h_antInfo) {
    nAntType = in->gpuInfo->h_antInfo->numAntType;
  }

  /* delete this class members */
  if (in->gpuInfo) {
    if (in->gpuInfo->d_visInfo) {
      ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_visInfo);
      in->gpuInfo->d_visInfo = NULL;
    }
    if (in->gpuInfo->h_visInfo) {
      if (in->gpuInfo->h_visInfo->d_freqRat) {
	ObitCUDAUtilFreeGPU((float*)in->gpuInfo->h_visInfo->d_freqRat);
	in->gpuInfo->h_visInfo->d_freqRat = NULL;
      }
      if (in->gpuInfo->h_visInfo->h_freqRat) {
	ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_visInfo->h_freqRat);
	in->gpuInfo->h_visInfo->h_freqRat = NULL;
      }
      ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_visInfo);
      in->gpuInfo->h_visInfo = NULL;
    }
    if (in->gpuInfo->d_freq) {
      ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_freq);
      in->gpuInfo->d_freq = NULL;
    }
    if (in->gpuInfo->d_model) {
      ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_model);
      in->gpuInfo->d_model = NULL;
    }
    if (in->gpuInfo->d_modelInfo) {
      ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_modelInfo);
      in->gpuInfo->d_modelInfo = NULL;
    }
    if (in->gpuInfo->h_modelInfo) {
      if (in->gpuInfo->h_modelInfo->d_specIndex) {
	ObitCUDAUtilFreeGPU((float*)in->gpuInfo->h_modelInfo->d_specIndex);
	in->gpuInfo->h_modelInfo->d_specIndex = NULL;
      }
      if (in->gpuInfo->h_modelInfo->h_specIndex) {
	ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_modelInfo->h_specIndex);
	in->gpuInfo->h_modelInfo->h_specIndex = NULL;
      }
      ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_modelInfo);
      in->gpuInfo->h_modelInfo = NULL;
    }

    /* Device antInfo */
    if (in->gpuInfo->d_antInfo) {
      if (in->gpuInfo->d_antInfo->d_BeamFreq) {
	ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_antInfo->d_BeamFreq);
	in->gpuInfo->d_antInfo->d_BeamFreq = NULL;
      }
      if (in->gpuInfo->d_antInfo->d_AntType) {
	ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_antInfo->d_AntType);
	in->gpuInfo->d_antInfo->d_AntType = NULL;
      }
      if (in->gpuInfo->d_antInfo->d_visJChann) {
	ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_antInfo->d_visJChann);
	in->gpuInfo->d_antInfo->d_visJChann = NULL;
      }
      if (in->gpuInfo->d_antInfo->d_Jones) {
	ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_antInfo->d_Jones);
	in->gpuInfo->d_antInfo->d_Jones = NULL;
      }
      if (in->gpuInfo->d_antInfo->d_BeamInterp) {
	/* Free interpolators */
	for (j=0; j<nAntType; j++) {
	  indx = j*nAntType;	  
	  for (i=0; i<8; i++) {
	    if (in->gpuInfo->d_antInfo->d_BeamInterp[indx+i]) {
	      CUDABeamInterpFree(in->gpuInfo->d_antInfo->d_BeamInterp[indx+i]);
	      in->gpuInfo->d_antInfo->d_BeamInterp[indx+i] =  NULL;
	    }
	  } /* end loop over Stokes product */
	} /* end loop over types */
	ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_antInfo->d_BeamInterp);
	in->gpuInfo->d_antInfo->d_BeamInterp = NULL;
      }
      ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_antInfo);
      in->gpuInfo->d_antInfo = NULL;
    } /* end if d_antInfo */
    
    /* Host antInfo */
    if (in->gpuInfo->h_antInfo) {
      if (in->gpuInfo->h_antInfo->h_BeamFreq) {
	ObitCUDAUtilFreeGPU((float*)in->gpuInfo->h_antInfo->h_BeamFreq);
	in->gpuInfo->h_antInfo->h_BeamFreq = NULL;
      }
      if (in->gpuInfo->h_antInfo->h_AntType) {
	ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_antInfo->h_AntType);
	in->gpuInfo->h_antInfo->h_AntType = NULL;
      }
      if (in->gpuInfo->h_antInfo->h_visJChann) {
	ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_antInfo->h_visJChann);
	in->gpuInfo->h_antInfo->h_visJChann = NULL;
      }
      ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_antInfo);
      in->gpuInfo->h_antInfo = NULL;
    } /* end if h_antInfo */

    /* Free antList */
    in->antList = ObitAntennaListUnref(in->antList);
    
    /* Free curSource */
    in->curSource = ObitAntennaListUnref(in->curSource);
    
    /* free top level */
    ObitCUDAUtilFreeHost((float*)in->gpuInfo); in->gpuInfo = NULL;
  } /* end if in->gpuinfo */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
#endif /* HAVE_GPU */
  
} /* end ObitGPUBeamModelClear */
