/* $Id: $    */
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

#include "ObitGPUFArray.h"
#if HAVE_GPU==1  /* GPU? Real versions */
#include "ObitCUDAUtil.h"
#endif /* HAVE_GPU */

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUFArray.c
 * ObitGPUFArray class function definitions.
 * This class is derived from the Obit base class.
 * ObitGPUFArrays consist of an array in locked host memory and corresponding
 * allocations in the GPU.  Data movement functions are provided.
 * Communications to the GPU is via CUDAFArray.
 * This class provide a limited subset of ObitFArray.
 * CUDA is used for the actual GPU routines but CUDA is REALLY primitive
 * so very basic inplementation is in CUDAFArray.cu and friends.
 * Portions of the class are in CUDA and are only implemented if the
 * compiler option -DHAVE_GPU=1 is used.  Some portions also need to 
 * variable IS_CUDA=1 to be set in the calling (CUDA) routines.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitGPUFArray";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitGPUFArrayClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGPUFArrayClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGPUFArrayInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGPUFArrayClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGPUFArrayClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitGPUFArray* newObitGPUFArray (gchar* name)
{
  ObitGPUFArray* out = NULL;
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  int memsize;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGPUFArrayClassInit();

  /* allocate/init structure in locked memory */
  memsize = sizeof(ObitGPUFArray);
  /*out = (ObitGPUFArray*)ObitCUDAUtilAllocHost(memsize);*/
  out = g_malloc0(memsize);

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGPUFArrayInit((gpointer)out);
#else /* Not GPU */
  g_warning ("Obit not built with GPU enabled");
#endif /* HAVE_GPU */
 return out;
} /* end newObitGPUFArray */

/**
 * Creates an GPUFArray of a specified geometry.
 * GPU resident structures created, initialized.
 * \param ndim  Number of dimensions desired, if <=0 data array not created.
 *              maximum value = MAXFARRAYDIM.
 * \param naxis Dimensionality along each axis. NULL => don't create array.
 * \return the new object.
 */
ObitGPUFArray* ObitGPUFArrayCreate (gchar *name, olong ndim, olong *naxis)
{
  ObitGPUFArray* out=NULL;
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  int i, size;

  /* Create basic structures - host */
  out =  newObitGPUFArray(name);
  /* CUDA version */
  out->FArray = CUDAFArrayCreate (ndim, naxis);

  /* copy geometry */
  out->ndim = ndim;
  size = 1; /* total size */
  for (i=0; i<ndim; i++) {
    out->naxis[i] = MAX (1, MIN(naxis[i],524288));  /* Not too big */
    size *= out->naxis[i]; /* total size */
  }
  out->size = size;
#else /* Not GPU */
  g_warning ("Obit not built with GPU enabled");
#endif /* HAVE_GPU */
  return out;
} /* end GPUFArrayCreate */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitGPUFArrayGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGPUFArrayClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGPUFArrayGetClass */

/**
 * Make a deep copy of an ObitGPUFArray.
 * NYI
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitGPUFArray* ObitGPUFArrayCopy  (ObitGPUFArray *in, ObitGPUFArray *out, ObitErr *err)
{
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
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
    out = newObitGPUFArray(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  g_error("ObitGPUFArrayCopy NOT implemented");

#endif /* HAVE_GPU */
  return out;
} /* end ObitGPUFArrayCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an GPUFArray similar to the input one.
 * NYI
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitGPUFArrayClone  (ObitGPUFArray *in, ObitGPUFArray *out, ObitErr *err)
{
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  const ObitClassInfo *ParentClass;

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  g_error("ObitGPUFArrayClone NOT implemented");
#endif /* HAVE_GPU */

} /* end ObitGPUFArrayClone */

/**
 * Copy host array to data
 * NB: Transfers to/from GPU need locked host memory
 * \param in   Object to copy
 * \param err  Obit error/message stack
 */
void ObitGPUFArrayToData (ObitGPUFArray *in, ofloat* data, 
			  ObitErr *err)
{
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  size_t memsize = in->size*sizeof(float);
  /* error checks */
  if (err->error) return;
  memmove ((float*)data, in->FArray->h_array,  memsize);
#endif /* HAVE_GPU */
} /* end ObitGPUFArrayToData */

/**
 * Copy from data to host array;
 * NB: Transfers to/from GPU need locked host memory
 * \param in   Object to copy
 * \param err  Obit error/message stack
 */
void ObitGPUFArrayFromData (ObitGPUFArray *in, ofloat* data, 
			    ObitErr *err)
{
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  size_t memsize = in->size*sizeof(float);
  /* error checks */
  if (err->error) return;
  memmove (in->FArray->h_array, (float*)data,  memsize);
#endif /* HAVE_GPU */
} /* end ObitGPUFArrayFromData */

/**
 * Copy data array from host to GPU
 * NB: Transfers to/from GPU need locked host memory
 * \param in   Object to copy
 * \param err  Obit error/message stack
 */
void ObitGPUFArrayToGPU (ObitGPUFArray *in, ObitErr *err)
{
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  int memsize = in->size*sizeof(float);
  /* error checks */
  if (err->error) return;
  ObitCUDAUtilHost2GPU (in->FArray->d_array, in->FArray->h_array, memsize, NULL);
#endif /* HAVE_GPU */
} /* end ObitGPUFArrayToGPU */

/**
 * Copy data array from GPU to Host
 * NB: Transfers to/from GPU need locked host memory
 * \param in   Object to copy
 * \param err  Obit error/message stack
 */
void ObitGPUFArrayToHost (ObitGPUFArray *in, ObitErr *err)
{
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  int memsize = in->size*sizeof(float);
  /* error checks */
  if (err->error) return;
  ObitCUDAUtilGPU2Host (in->FArray->h_array, in->FArray->d_array, memsize, NULL);
#endif /* HAVE_GPU */
} /* end ObitGPUFArrayToGPU */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitGPUFArrayClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGPUFArrayClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGPUFArrayClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGPUFArrayClassInfoDefFn (gpointer inClass)
{
  ObitGPUFArrayClassInfo *theClass = (ObitGPUFArrayClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGPUFArrayClassInit;
  theClass->newObit       = (newObitFP)newObitGPUFArray;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGPUFArrayClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGPUFArrayGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitGPUFArrayCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitGPUFArrayClear;
  theClass->ObitInit      = (ObitInitFP)ObitGPUFArrayInit;
  theClass->ObitGPUFArrayCreate   = (ObitGPUFArrayCreateFP)ObitGPUFArrayCreate;
  theClass->ObitGPUFArrayToData   = (ObitGPUFArrayToDataFP)ObitGPUFArrayToData;
  theClass->ObitGPUFArrayFromData = (ObitGPUFArrayFromDataFP)ObitGPUFArrayFromData;
  theClass->ObitGPUFArrayToGPU    = (ObitGPUFArrayToGPUFP)ObitGPUFArrayToGPU;
  theClass->ObitGPUFArrayToHost   = (ObitGPUFArrayToHostFP)ObitGPUFArrayToHost;
} /* end ObitGPUFArrayClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitGPUFArrayInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGPUFArray *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->ndim   = 0;
  in->naxis[0] = in->naxis[1] = in->naxis[2] = in->naxis[3] = in->naxis[4];
  in->size   = 0;
  in->FArray = NULL;
} /* end ObitGPUFArrayInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitGPUFArray* cast to an Obit*.
 */
void ObitGPUFArrayClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGPUFArray *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

#if HAVE_GPU==1  /* GPU? Real versions */
  /* delete this class members */
  if (in) {
    if (in->FArray) CUDAFArrayZap(in->FArray);
  }
#endif /* HAVE_GPU */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
     
  /* Delete basic object NO already done 
  ObitCUDAUtilFreeHost((float*)in);*/
 
} /* end ObitGPUFArrayClear */
