/* Still need 
1)  model to write to GPU
2)  set nvis per IO to appropriate for GPU
3)  reset nvis per IO
4)  Gaussian model 
*/
/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
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

#include "ObitGPUSkyModel.h"
#include "ObitCUDAUtil.h"
#include "ObitCUDASkyModelInfoDef.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUSkyModel.c
 * ObitGPUSkyModel class function definitions.
 * This class is derived from the Obit base class.
 * The ObitGPUSkyModel class uses GPUs to calculate sky models.
 * CUDA is used for the actual GPU routines but CUDA is REALLY primitive
 * so very basic inplementation is in ObitCUDASkyModel.cu and friends.
 * Portions of the class are in CUDA and are only implemented if the
 * compiler option -DHAVE_GPU=1 is used.  Some portions also need to 
 * variable IS_CUDA=1 to be set in the calling routines.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitGPUSkyModel";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitGPUSkyModelClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGPUSkyModelClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGPUSkyModelInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGPUSkyModelClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGPUSkyModelClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitGPUSkyModel* newObitGPUSkyModel (gchar* name)
{
  ObitGPUSkyModel* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGPUSkyModelClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitGPUSkyModel));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGPUSkyModelInit((gpointer)out);

 return out;
} /* end newObitGPUSkyModel */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitGPUSkyModelGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGPUSkyModelClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGPUSkyModelGetClass */

/**
 * Make a deep copy of an ObitGPUSkyModel.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitGPUSkyModel* ObitGPUSkyModelCopy  (ObitGPUSkyModel *in, ObitGPUSkyModel *out, ObitErr *err)
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
    out = newObitGPUSkyModel(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */

  return out;
} /* end ObitGPUSkyModelCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an GPUSkyModel similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitGPUSkyModelClone  (ObitGPUSkyModel *in, ObitGPUSkyModel *out, ObitErr *err)
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

} /* end ObitGPUSkyModelClone */

/**
 * Creates an ObitGPUSkyModel 
 * Currently only DFT point supported
 * \param name  An optional name for the object.
 * \param type  mode comp. type 'DFT' or 'GRID'.
 * \return the new object.
 */
ObitGPUSkyModel* ObitGPUSkyModelCreate (gchar* name, gchar *type)
{
  ObitGPUSkyModel* out;

  /* Create basic structure */
  out = newObitGPUSkyModel (name);

  return out;
} /* end ObitGPUSkyModelCreate */

/**
 * Initialize an ObitGPUSkyModel 
 * Currently only DFT point supported
 * \param in      Iniput object
 * \param uvdata  UV data being operated on
 * \param err     Obit error stack object.
 */
void ObitGPUSkyModelDFTInit (ObitGPUSkyModel *in, ObitUV *uvdata, ObitErr *err)
{
  ObitUVDesc *inDesc = uvdata->myDesc;
  int i, memsize, kincf, kincif, nchan, nif;
  ObitIOAccess access;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean doCalSelect;
  olong nVisPIO, oldnVisPIO;
  gchar *routine = "ObitGPUSkyModelDFTInit";

  /* Set GPU device, hard code 0 for now */
  in->gpuInfo->cuda_device = 0;
  ObitCUDASetGPU (in->gpuInfo->cuda_device);
  ObitCUDAResetGPU ();  /* Make sure reset */

  /* Print level */
  in->gpuInfo->prtLv = err->prtLv;

  /* reallocate data buffer on uvdata to locked memory
    first reset number of vis per IO */
  ObitUVClose (uvdata, err);
  oldnVisPIO = MIN(1000, inDesc->nvis);
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

  in->gpuInfo->bufferSize = uvdata->bufferSize;
  memsize = uvdata->bufferSize * sizeof(ofloat);
  in->gpuInfo->h_data = (ofloat*)ObitCUDAUtilAllocHost(memsize);

  in->gpuInfo->nvis    = nVisPIO;
  in->gpuInfo->nchan   = inDesc->inaxes[inDesc->jlocf];
  in->gpuInfo->nrparm  = inDesc->nrparm;
  in->gpuInfo->lenvis  = inDesc->lrec;
  in->gpuInfo->nstream = STREAM_COUNT;

  /* Allocate resources */
  /* Device Frequency scaling array from uv descriptor */
  nchan   = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocif>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else                   nif = 1;
  memsize = nchan * nif * sizeof(float);
  in->gpuInfo->d_freq = ObitCUDAUtilAllocGPU (memsize);
  ObitCUDAUtilHost2GPU (in->gpuInfo->d_freq, inDesc->fscale, memsize, NULL);
  /* Need locked version of fscale??? */

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
  /* Visibility info */
  memsize = sizeof(GPUVisInfo);
  in->gpuInfo->d_visInfo = (GPUVisInfo *)ObitCUDAUtilAllocGPU(memsize);
  in->gpuInfo->h_visInfo = (GPUVisInfo *)ObitCUDAUtilAllocHost(memsize);
  in->gpuInfo->h_visInfo->freqScale = in->gpuInfo->d_freq;
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
  in->gpuInfo->h_visInfo->incf   = inDesc->incf/3;  /* In units of correlations */
  in->gpuInfo->h_visInfo->nIF    = inDesc->inaxes[inDesc->jlocif];
  in->gpuInfo->h_visInfo->IFb    = uvdata->mySel->startIF-1;
  in->gpuInfo->h_visInfo->IFe    = uvdata->mySel->startIF+uvdata->mySel->numberIF-2;
  in->gpuInfo->h_visInfo->incif  = inDesc->incif/3;
  in->gpuInfo->h_visInfo->nstok  = inDesc->inaxes[inDesc->jlocs];
  in->gpuInfo->h_visInfo->stokb  = 0;  // ???
  in->gpuInfo->h_visInfo->stoke  = 1;
  in->gpuInfo->h_visInfo->incs   = inDesc->incs/3;
  in->gpuInfo->h_visInfo->kincif = kincif; // IF increment in freqScale
  in->gpuInfo->h_visInfo->kincf  = kincf;  // Freq increment in freqScale
  if (inDesc->crval[inDesc->jlocs]>0.0)       in->gpuInfo->h_visInfo->stype  = 0;
  else if (inDesc->crval[inDesc->jlocs]<-4.0) in->gpuInfo->h_visInfo->stype  = 2;
  else                                        in->gpuInfo->h_visInfo->stype  = 1;
  /* to GPU */
  ObitCUDAUtilHost2GPU ((float*)in->gpuInfo->d_visInfo, (float*)in->gpuInfo->h_visInfo, 
			memsize, NULL);

  /* Model info */
  memsize = sizeof(GPUModelInfo);
  in->gpuInfo->d_modelInfo = (GPUModelInfo*)ObitCUDAUtilAllocGPU(memsize);
  in->modelInfo = (GPUModelInfo*)ObitCUDAUtilAllocHost(memsize);
  in->modelInfo->nmodel = in->gpuInfo->nmodel;
  in->modelInfo->size   = in->gpuInfo->modelSize;
  in->modelInfo->type   = 0;  /* MORE HERE */
  in->modelInfo->model  = in->gpuInfo->d_model;

  /* Allocate host/device arrays */
  in->gpuInfo->d_data_in  = g_malloc0(in->gpuInfo->nstream*sizeof(float*));
  in->gpuInfo->d_data_out = g_malloc0(in->gpuInfo->nstream*sizeof(float*));
  memsize = in->gpuInfo->nvis * in->gpuInfo->lenvis * 
    sizeof(float)/in->gpuInfo->nstream;  // size of vis buffer
  for (i=0; i<in->gpuInfo->nstream; i++) {
    in->gpuInfo->d_data_in[i]  = ObitCUDAUtilAllocGPU(memsize);
    in->gpuInfo->d_data_out[i] = ObitCUDAUtilAllocGPU(memsize);
  }

  /* Allocate stream stuff  */
  in->gpuInfo->stream    = g_malloc0(in->gpuInfo->nstream*sizeof(void*));
  in->gpuInfo->cycleDone = g_malloc0(in->gpuInfo->nstream*sizeof(void*));
  for (i=0; i<in->gpuInfo->nstream; ++i)
    {
      in->gpuInfo->stream[i]    = ObitCUDAStreamCreate();
      in->gpuInfo->cycleDone[i] = ObitCUDAEventCreate();
      ObitCUDAEventRecord (in->gpuInfo->cycleDone[i], in->gpuInfo->stream[i]);
    }

  /* CALL CUDA; AND CURSE; */
  ObitCUDASkyModelDFTInit (in->gpuInfo, in->visInfo, in->modelInfo);
} /* end  ObitGPUSkyModelDFTInit */


/**
 * Setup model for ObitGPUSkyModel 
 * Currently only DFT point supported
 * \param in      Input object
 * \param model   model components
 * \param err     Obit error stack object.
 */
void ObitGPUSkyModelDFTSetMod (ObitGPUSkyModel *in, ObitFArray *model, ObitErr *err)
{
  int i, inc, high, memsize;

  /* Find last non zero entry */
  inc = model->naxis[0];
  high = 1;
  for (i=0; i<model->naxis[1]; i++) if (model->array[i+inc]!=0.0) high = i+1;

  /* (re)build model array */
  if ((in->gpuInfo->nmodel!=high) || (in->gpuInfo->modelSize!=model->naxis[0])) {
    in->gpuInfo->nmodel    = high;
    in->gpuInfo->modelSize = model->naxis[0];
    if (in->gpuInfo->d_model) ObitCUDAUtilFreeGPU(in->gpuInfo->d_model);
    /* new Device Model array */
    memsize = in->gpuInfo->nmodel * in->gpuInfo->modelSize * sizeof(float);
    in->gpuInfo->d_model = ObitCUDAUtilAllocGPU(memsize);
    in->modelInfo->model = in->gpuInfo->d_model;
  }
  /* Copy model to device */
  memsize = in->gpuInfo->nmodel * in->gpuInfo->modelSize * sizeof(float);
  ObitCUDAUtilHost2GPU(in->gpuInfo->d_model, model->array, memsize, NULL);

  /* Copy Model info to device */
  in->modelInfo->nmodel = in->gpuInfo->nmodel;
  in->modelInfo->size   = in->gpuInfo->modelSize;
  memsize = sizeof(GPUModelInfo);
  ObitCUDAUtilHost2GPU((float*)in->gpuInfo->d_modelInfo, (float*)in->modelInfo, 
		       memsize, NULL);

  /* CALL CUDA; AND CURSE; */
  ObitCUDASkyModelDFTSetMod (in->gpuInfo, in->visInfo, in->modelInfo);
  return;
} /* end ObitGPUSkyModelDFTSetMod */

/**
 * Calculate an ObitGPUSkyModel 
 * Currently only DFT point supported
 * \param in      Iniput object
 * \param uvdata  UV data being operated on
 * \param err     Obit error stack object.
 */
void ObitGPUSkyModelDFTCalc (ObitGPUSkyModel *in, ObitUV *uvdata, ObitErr *err)
{
  size_t memsize;
   /* Reset number of vis */
   in->gpuInfo->nvis = uvdata->myDesc->numVisBuff;
   in->visInfo->nvis = uvdata->myDesc->numVisBuff;
   in->visInfo->nprod = uvdata->myDesc->ncorr;

   /* Copy data to locked buffer */
   memsize = in->gpuInfo->nvis * in->gpuInfo->lenvis * sizeof(ofloat);
   memcpy (in->gpuInfo->h_data, uvdata->buffer, memsize);

  /* CALL CUDA; AND CURSE; */
  ObitCUDASkyModelDFTCalc (in->gpuInfo, in->visInfo, in->modelInfo);

  /* Copy data back */
  memcpy (uvdata->buffer, in->gpuInfo->h_data, memsize);
  return;
} /* end ObitGPUSkyModelDFTCalc */


/**
 * Shutdown  ObitGPUSkyModel 
 * Currently only DFT point supported
 * \param in      Iniput object
 * \param uvdata  UV data being operated on
 * \param err     Obit error stack object.
 */
void ObitGPUSkyModelDFTShutdown (ObitGPUSkyModel *in, ObitUV *uvdata, ObitErr *err)
{
  int i;
  gint32 dim[MAXINFOELEMDIM];
  // Free resources
  ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_data);
  ObitCUDAUtilFreeHost((float*)in->gpuInfo->h_visInfo);
  ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_visInfo);
  ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_modelInfo);
  ObitCUDAUtilFreeGPU((float*)in->gpuInfo->d_model);
  for (i =0; i<in->gpuInfo->nstream; ++i)    {
    ObitCUDAUtilFreeGPU(in->gpuInfo->d_data_in[i]);
    ObitCUDAUtilFreeGPU(in->gpuInfo->d_data_out[i]);
    ObitCUDAStreamDestroy((in->gpuInfo->stream)[i]);
    ObitCUDAEventDestroy((in->gpuInfo->cycleDone)[i]);
  }
  g_free(in->gpuInfo->stream);
  g_free(in->gpuInfo->cycleDone);
  
  ObitCUDAResetGPU ();

   /* Replace old nVisPIO */
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (uvdata->info, "nVisPIO", OBIT_long, dim,  
			 &in->gpuInfo->oldnVisPIO);

 /* CALL CUDA; AND CURSE; */
  ObitCUDASkyModelDFTShutdown (in->gpuInfo, in->visInfo, in->modelInfo);
} /*end ObitGPUSkyModelDFTShutdown */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitGPUSkyModelClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGPUSkyModelClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGPUSkyModelClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGPUSkyModelClassInfoDefFn (gpointer inClass)
{
  ObitGPUSkyModelClassInfo *theClass = (ObitGPUSkyModelClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGPUSkyModelClassInit;
  theClass->newObit       = (newObitFP)newObitGPUSkyModel;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGPUSkyModelClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGPUSkyModelGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitGPUSkyModelCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitGPUSkyModelClear;
  theClass->ObitInit      = (ObitInitFP)ObitGPUSkyModelInit;
  theClass->ObitGPUSkyModelCreate = (ObitGPUSkyModelCreateFP)ObitGPUSkyModelCreate;
  theClass->ObitGPUSkyModelDFTInit = (ObitGPUSkyModelDFTInitFP)ObitGPUSkyModelDFTInit;
  theClass->ObitGPUSkyModelDFTSetMod = (ObitGPUSkyModelDFTSetModFP)ObitGPUSkyModelDFTSetMod;
  theClass->ObitGPUSkyModelDFTCalc = (ObitGPUSkyModelDFTCalcFP)ObitGPUSkyModelDFTCalc;
  theClass->ObitGPUSkyModelDFTShutdown = (ObitGPUSkyModelDFTShutdownFP)ObitGPUSkyModelDFTShutdown;
} /* end ObitGPUSkyModelClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitGPUSkyModelInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGPUSkyModel *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->gpuInfo    = g_malloc0(sizeof(GPUInfo));
  in->visInfo    = g_malloc0(sizeof(GPUVisInfo));
  in->modelInfo  = g_malloc0(sizeof(GPUModelInfo));
} /* end ObitGPUSkyModelInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitGPUSkyModel* cast to an Obit*.
 */
void ObitGPUSkyModelClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGPUSkyModel *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->gpuInfo)   g_free(in->gpuInfo);
  if (in->visInfo)   g_free(in->visInfo);
  if (in->modelInfo) ObitCUDAUtilFreeHost((float*)in->modelInfo);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitGPUSkyModelClear */
