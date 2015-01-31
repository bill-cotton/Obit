/* $Id: $   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014-2015                                          */
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

#include "ObitGPUFInterpolate.h"
#if HAVE_GPU==1  /* GPU? Real versions */
#include "ObitCUDAUtil.h"
#include "CUDAFInterpolate.h"
#endif /* HAVE_GPU */

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUFInterpolate.c
 * ObitGPUFInterpolate class function definitions.
 * This class is derived from the Obit base class.
 * The ObitGPUFInterpolate class uses GPUs to calculate sky models.
 * CUDA is used for the actual GPU routines but CUDA is REALLY primitive
 * so very basic inplementation is in ObitCUDAFInterpolate.cu and friends.
 * Portions of the class are in CUDA and are only implemented if the
 * compiler option -DHAVE_GPU=1 is used.  Some portions also need to 
 * variable IS_CUDA=1 to be set in the calling routines.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitGPUFInterpolate";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitGPUFInterpolateClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGPUFInterpolateClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGPUFInterpolateInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGPUFInterpolateClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGPUFInterpolateClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitGPUFInterpolate* newObitGPUFInterpolate (gchar* name)
{
  ObitGPUFInterpolate* out=NULL;
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  int memsize;
/* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGPUFInterpolateClassInit();

  /* allocate/init structure */
  memsize = sizeof(ObitGPUFInterpolate);
  /*out = (ObitGPUFInterpolate*) ObitCUDAUtilAllocHost(memsize);*/
  out = g_malloc0(memsize);

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGPUFInterpolateInit((gpointer)out);

#else /* Not GPU */
  g_warning ("Obit not built with GPU enabled");
#endif /* HAVE_GPU */
 return out;
} /* end newObitGPUFInterpolate */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitGPUFInterpolateGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGPUFInterpolateClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGPUFInterpolateGetClass */

/**
 * Make a deep copy of an ObitGPUFInterpolate.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitGPUFInterpolate* ObitGPUFInterpolateCopy  (ObitGPUFInterpolate *in, 
					       ObitGPUFInterpolate *out, ObitErr *err)
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
    out = newObitGPUFInterpolate(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */

 #endif /* HAVE_GPU */
 return out;
} /* end ObitGPUFInterpolateCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an GPUFInterpolate similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitGPUFInterpolateClone  (ObitGPUFInterpolate *in, ObitGPUFInterpolate *out, ObitErr *err)
{
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
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
#endif /* HAVE_GPU */

} /* end ObitGPUFInterpolateClone */

/**
 * Creates an ObitGPUFInterpolate 
 * Currently only DFT point supported
 * \param name  An optional name for the object.
 * \param inArray  Array to be interpolated
 * \param xArray   Array of input x pixels for output
 *                 MUST be the same geometry as output array
 * \param yArray   Array of input y pixels for output
 * \param hwidth   Convolution kernal half width
 * \param err      Obit message/error stack
 * \return the new object.
 */
ObitGPUFInterpolate* ObitGPUFInterpolateCreate (gchar* name, ObitFArray *inArray, 
						ObitFArray *xArray, ObitFArray *yArray,
						int hwidth, ObitErr *err)
{
  ObitGPUFInterpolate* out=NULL;
  gchar *routine = "ObitGPUFInterpolateCreate";

  if (err->error) return out;  /* previous error */

  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  int memsize;
  out = newObitGPUFInterpolate(name);

  /* Make GPU versions of arrays */
  out->inArray = ObitGPUFArrayCreate("inArray", inArray->ndim, inArray->naxis);
  if (err->error) Obit_traceback_val (err, routine, inArray->name, out);

  out->xArray = ObitGPUFArrayCreate("xArray", xArray->ndim, xArray->naxis);
  ObitGPUFArrayFromData(out->xArray, xArray->array, err);
  ObitGPUFArrayToGPU (out->xArray, err);
  if (err->error) Obit_traceback_val (err, routine, inArray->name, out);

  out->yArray = ObitGPUFArrayCreate("yArray", yArray->ndim, yArray->naxis);
  ObitGPUFArrayFromData(out->yArray, yArray->array, err);
  ObitGPUFArrayToGPU (out->yArray, err);
  if (err->error) Obit_traceback_val (err, routine, inArray->name, out);

 /* GPU FInterpolate structure */
  memsize = sizeof(CUDAFInterpolate);
  out->FInterpolate = newCUDAFInterpolateCreate(out->inArray->FArray, hwidth, STREAM_COUNT);
  out->hwidth  = MAX (1, MIN (4, hwidth));
  out->nx      = inArray->naxis[0];
  out->ny      = inArray->naxis[1];

#else /* Not GPU */
  Obit_log_error(err, OBIT_Error, "%s: Obit not built with GPU enabled", 
		 routine);
#endif /* HAVE_GPU */
   return out;
} /* end ObitGPUFInterpolateCreate */

/**
 * Interpolate value at requested pixel in a plane of an n(>=2)-D array.
 * Interpolation between planes is not supported.
 * \param in       The object to interpolate
 * \param Array    Input pixel array
 * \param outArray Output pixel array
 * \param err      Obit message/error stack
 */
void ObitGPUFInterpolateImage (ObitGPUFInterpolate *in, ObitFArray *inArray,
			       ObitFArray *outArray, ObitErr *err)
{
  gchar *routine = "ObitGPUFInterpolateImage";
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  int i, memsize, nrowBuff, nx=outArray->naxis[0], ny=outArray->naxis[1]; 
  /* Compatibility checks */
  Obit_return_if_fail ((nx==in->xArray->naxis[0]) && (ny==in->xArray->naxis[1]), err, 
		       "outArray incompatable with xArray %d!=%d or %d!=%d ",
		       nx,in->xArray->naxis[0],ny,in->xArray->naxis[1]);
  Obit_return_if_fail ((nx==in->yArray->naxis[0]) && (ny==in->yArray->naxis[1]), err, 
		       "outArray incompatable with yArray %d!=%d or %d!=%d ",
		       nx,in->yArray->naxis[0],ny,in->yArray->naxis[1]);

  /* (re)create GPU buffers if needed */
  if ((in->FInterpolate->nxBuff!=nx) || (in->FInterpolate->nyBuff!=ny)) {
    /* Delete old */
    if (in->buffer) ObitCUDAUtilFreeHost((float*)in->buffer);
    for (i=0; i<in->FInterpolate->nstream; ++i)    {
      if ((in->FInterpolate->d_data)&&(in->FInterpolate->d_data[i])) 
	ObitCUDAUtilFreeGPU(in->FInterpolate->d_data[i]);
    }
    /* Create new */
    nrowBuff = 1 + (ny-1)/in->FInterpolate->nstream;
    memsize = nx*ny*sizeof(float);
    in->buffer = ObitCUDAUtilAllocHost(memsize);
    memsize = nrowBuff*nx*sizeof(float);
    for (i=0; i<in->FInterpolate->nstream; ++i)    {
      in->FInterpolate->d_data[i] = ObitCUDAUtilAllocGPU(memsize);
    }
  }
  
  /* Copy input data to GPU */
  ObitGPUFArrayFromData(in->inArray, inArray->array, err);
  ObitGPUFArrayToGPU (in->inArray, err);
  if (err->error) Obit_traceback_msg (err, routine, inArray->name);

  /* Do interpolation */
  CUDAFInterpolateImage (in->FInterpolate, in->xArray->FArray, 
			 in->yArray->FArray, in->buffer);

  /* Copy data to outArray */
  memsize = nx*ny*sizeof(float);
  memmove(outArray->array, in->buffer, memsize);
#else /* Not GPU */
  Obit_log_error(err, OBIT_Error, "%s: Obit not built with GPU enabled", 
		 routine);
#endif /* HAVE_GPU */
} /* end GPUFInterpolateImage */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitGPUFInterpolateClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGPUFInterpolateClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGPUFInterpolateClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGPUFInterpolateClassInfoDefFn (gpointer inClass)
{
  ObitGPUFInterpolateClassInfo *theClass = (ObitGPUFInterpolateClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGPUFInterpolateClassInit;
  theClass->newObit       = (newObitFP)newObitGPUFInterpolate;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGPUFInterpolateClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGPUFInterpolateGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitGPUFInterpolateCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitGPUFInterpolateClear;
  theClass->ObitInit      = (ObitInitFP)ObitGPUFInterpolateInit;
  theClass->ObitGPUFInterpolateCreate = (ObitGPUFInterpolateCreateFP)ObitGPUFInterpolateCreate;
  theClass->ObitGPUFInterpolateImage  = (ObitGPUFInterpolateImageFP)ObitGPUFInterpolateImage;
} /* end ObitGPUFInterpolateClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitGPUFInterpolateInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGPUFInterpolate *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->FInterpolate = NULL;
  in->inArray      = NULL;
  in->xArray       = NULL;
  in->yArray       = NULL;
  in->buffer       = NULL;
  in->nx = in->ny = 0;
  in->hwidth = 0;
} /* end ObitGPUFInterpolateInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitGPUFInterpolate* cast to an Obit*.
 */
void ObitGPUFInterpolateClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGPUFInterpolate *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

#if HAVE_GPU==1  /* GPU? Real versions */
  /* delete this class members */
  CUDAFInterpolateZap(in->FInterpolate);  /* This gets the input CUDAFArray */
  in->inArray->FArray = NULL;  /* Already Zapped */
  in->inArray = ObitGPUFArrayUnref(in->inArray);
  in->xArray  = ObitGPUFArrayUnref(in->xArray);
  in->yArray  = ObitGPUFArrayUnref(in->yArray);
  if (in->buffer) ObitCUDAUtilFreeHost((float*)in->buffer);
 #endif /* HAVE_GPU */
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

  /* Delete basic object NO already done
  ObitCUDAUtilFreeHost((float*)in); */
  
} /* end ObitGPUFInterpolateClear */
