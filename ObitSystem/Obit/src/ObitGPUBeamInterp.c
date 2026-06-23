/* $Id: $   */
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

#include "ObitGPUBeamInterp.h"
#if HAVE_GPU==1  /* GPU? Real versions */
#include "ObitCUDAUtil.h"
#include "CUDABeamInterp.h"
#endif /* HAVE_GPU */

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUBeamInterp.c
 * ObitGPUBeamInterp class function definitions.
 * This class is derived from the Obit base class.
 * The ObitGPUBeamInterp class uses GPUs to calculate sky models.
 * CUDA is used for the actual GPU routines but CUDA is REALLY primitive
 * so very basic inplementation is in ObitCUDABeamInterp.cu and friends.
 * Portions of the class are in CUDA and are only implemented if the
 * compiler option -DHAVE_GPU=1 is used.  Some portions also need to 
 * variable IS_CUDA=1 to be set in the calling routines.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitGPUBeamInterp";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitGPUBeamInterpClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGPUBeamInterpClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGPUBeamInterpInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGPUBeamInterpClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGPUBeamInterpClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitGPUBeamInterp* newObitGPUBeamInterp (gchar* name)
{
  ObitGPUBeamInterp* out=NULL;
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  int memsize;
/* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGPUBeamInterpClassInit();

  /* allocate/init structure */
  memsize = sizeof(ObitGPUBeamInterp);
  /*out = (ObitGPUBeamInterp*) ObitCUDAUtilAllocHost(memsize);*/
  out = g_malloc0(memsize);

  /* initialize other stuff */
  ObitGPUBeamInterpInit((gpointer)out);

#else /* Not GPU */
  g_warning ("Obit not built with GPU enabled");
#endif /* HAVE_GPU */
 return out;
} /* end newObitGPUBeamInterp */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitGPUBeamInterpGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGPUBeamInterpClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGPUBeamInterpGetClass */

/**
 * Make a deep copy of an ObitGPUBeamInterp.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitGPUBeamInterp* ObitGPUBeamInterpCopy  (ObitGPUBeamInterp *in, 
					   ObitGPUBeamInterp *out, ObitErr *err)
{
  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar* outName="BeamInterp";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    out = newObitGPUBeamInterp(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */

 #endif /* HAVE_GPU */
 return out;
} /* end ObitGPUBeamInterpCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an GPUBeamInterp similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitGPUBeamInterpClone  (ObitGPUBeamInterp *in, ObitGPUBeamInterp *out, ObitErr *err)
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

} /* end ObitGPUBeamInterpClone */

/**
 * Creates an ObitGPUBeamInterp and member CUDABeamInterp
 * \param name  An optional name for the object.
 * \param inArray  ObitFArray to be interpolated
 * \param Desc     Descriptor of image to be interpolated
 * \param hwidth   Convolution kernal half width
 * \param err      Obit message/error stack
 * \return the new object.
 */
ObitGPUBeamInterp* ObitGPUBeamInterpCreate (gchar* name, ObitFArray *inArray, 
					    ObitImageDesc *Desc, int hwidth,
					    ObitErr *err)
{
  ObitGPUBeamInterp* out=NULL;
  /*gchar *routine = "ObitGPUBeamInterpCreate";*/

  if (err->error) return out;  /* previous error */

  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* GPU? Real versions */
  out = newObitGPUBeamInterp(name);

  /* Create CUDABeamInterp */
  out->beamInterp =  ObitGPUCUDABeamInterpCreate (inArray, Desc, hwidth);

 return out;
#else /* Not GPU */
 gchar *routine = "ObitGPUBeamInterpCreate";
  Obit_log_error(err, OBIT_Error, "%s: Obit not built with GPU enabled", 
		 routine);
#endif /* HAVE_GPU */
   return out;
} /* end ObitGPUBeamInterpCreate */

  /* Only if building with GPU enabled */
#if HAVE_GPU==1  /* Have a GPU? */
/**
 * Create CUDABeamInterp and copy info from host to device
 * The usage for Beam correction differs from that in CUDABeamInterp
 * so different constructor and destructor is needed.
 * These are for interpolating the beam images.
 * This has to be in a pure c routine as the glib features used in Obit classes
 * don't play well in cuda.
 * \param  array    ObitFArray to be interpolated.
 * \param  Desc     Descriptor of image to be interpolated
 * \param  hwidth   half-width in pixels of interpolation kernal
 * \return host pointer for CUDABeamInterp, contents copied
 * should be Freeed using CUDABeamInterpFree
 */
CUDABeamInterp* ObitGPUCUDABeamInterpCreate (ObitFArray *array, ObitImageDesc *Desc, int hwidth)
{
  CUDABeamInterp* out=NULL;
  CUDAImageDesc *h_myDesc=NULL;
  int iwid, i, j, memsize;
  float prod;

  /* Create/init output structure */
  memsize = sizeof(CUDABeamInterp);  /* Basic structure in host */
  out = (CUDABeamInterp*) ObitCUDAUtilAllocHost(memsize);
  out->d_BeamInterp = NULL; /* device memory allocated */
  /* Copy FArray to device  */
  memsize = sizeof(CUDAFArray);  /* Basic structure in host */
  out->h_myArray = (CUDAFArray*) ObitCUDAUtilAllocHost(memsize); /* host */
  out->d_myArray = (CUDAFArray*) ObitCUDAUtilAllocGPU(memsize);  /* device */
  /* Fill in basic info */
  out->h_myArray->fblank = ObitMagicF();
  out->h_myArray->ndim   = array->ndim;
  for (i=0; i<array->ndim; i++)  out->h_myArray->naxis[i] = (int)array->naxis[i];
  /* pixel array in device */
  memsize = array->arraySize*sizeof(float);
  /* Need host copy in locked memory */
  out->h_myArray->h_array = ObitCUDAUtilAllocHost(memsize);
  for (i=0; i<array->arraySize; i++) out->h_myArray->h_array[i] = (float)array->array[i];
  out->h_myArray->d_array = (float*)ObitCUDAUtilAllocGPU(memsize);  /* allocate device */
  out->h_myArray->d_FArray = NULL;  /* Just in case */
  ObitCUDAUtilHost2GPU (out->h_myArray->d_array, (float*)out->h_myArray->h_array, memsize, NULL);
  memsize = sizeof(CUDAFArray);  /* Copy basic structure to device */
  ObitCUDAUtilHost2GPU ((float*)out->d_myArray, (float*)out->h_myArray, memsize, NULL);

  /* Image descriptor */
  memsize = sizeof(CUDAImageDesc);  /* Basic structure  */
  out->d_myDesc = (CUDAImageDesc*) ObitCUDAUtilAllocGPU(memsize);  /* device */
  h_myDesc = ObitGPUSkyGeomImageH2D (Desc);  /* Host copy in locked memory */
  ObitCUDAUtilHost2GPU ((float*)out->d_myDesc, (float*)h_myDesc, memsize, NULL);
  /* Free work memory */
  ObitGPUImageDescZap (h_myDesc);  h_myDesc = NULL;

  /* GPU structures for BeamInterp  */
  memsize = sizeof(CUDABeamInterp);  /* Device memory for interpolator */
  out->d_BeamInterp = ObitCUDAUtilAllocGPU(memsize);
  /* Set interpolator Kernal halfwidth */
  out->hwidth  = MAX (1, MIN (4, hwidth));

  /* Get array size info */
  out->nx     = array->naxis[0];
  out->ny     = array->naxis[1];

  /* Init Lagrangian denominators for hwidth */
  iwid = 1 + (2*out->hwidth);
  for (j= 1; j<=iwid; j++) {
    prod = 1.0;
    for (i= 1; i<=iwid; i++) {
      if (i != j) prod = prod * (j - i);
    } 
    out->denom[j-1] = 1.0 / prod;
  } 

  /* copy basic FInterpolator structure to device */
  memsize = sizeof(CUDABeamInterp);
  ObitCUDAUtilHost2GPU ((float*)out->d_BeamInterp, (float*)out, memsize, NULL);

  return out;
} /* end ObitGPUCUDABeamInterpCreate  */

/**
 * Free a device CUDABeamInterp
 * \param  in  device pointer to  CUDABeamInterp
 */
void CUDABeamInterpFree (CUDABeamInterp* in)
{
  if (!in) return;  /* Anybody home? */

  /* FArray */
  if (in->d_myArray) {
    if (in->d_myArray->d_array) {
      ObitCUDAUtilFreeGPU((float*)in->d_myArray->d_array); in->d_myArray->d_array = NULL;
    }
    ObitCUDAUtilFreeGPU((float*)in->d_myArray); in->d_myArray = NULL;
  }

  /* Image descriptor */
  if (in->d_myDesc) {
    ObitCUDAUtilFreeGPU((float*)in->d_myDesc); in->d_myDesc = NULL;
  }

  /* Basic structure */
  if (in->d_BeamInterp) {
    ObitCUDAUtilFreeGPU((float*)in->d_BeamInterp); in->d_BeamInterp = NULL;
  }
} /* end CUDABeamInterpFree  */
#endif /* HAVE_GPU */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitGPUBeamInterpClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGPUBeamInterpClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGPUBeamInterpClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGPUBeamInterpClassInfoDefFn (gpointer inClass)
{
  ObitGPUBeamInterpClassInfo *theClass = (ObitGPUBeamInterpClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGPUBeamInterpClassInit;
  theClass->newObit       = (newObitFP)newObitGPUBeamInterp;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGPUBeamInterpClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGPUBeamInterpGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitGPUBeamInterpCopy;
  theClass->ObitClone     = (ObitCloneFP)ObitGPUBeamInterpClone;
  theClass->ObitClear     = (ObitClearFP)ObitGPUBeamInterpClear;
  theClass->ObitInit      = (ObitInitFP)ObitGPUBeamInterpInit;
#if HAVE_GPU==1  /* GPU? Real versions */
  theClass->ObitGPUBeamInterpCreate = (ObitGPUBeamInterpCreateFP)ObitGPUBeamInterpCreate;
#endif /* HAVE_GPU */
} /* end ObitGPUBeamInterpClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitGPUBeamInterpInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGPUBeamInterp *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->beamInterp = NULL;
} /* end ObitGPUBeamInterpInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitGPUBeamInterp* cast to an Obit*.
 */
void ObitGPUBeamInterpClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGPUBeamInterp *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

#if HAVE_GPU==1  /* GPU? Real versions */
  /* delete this class members */
  /* CUDABeamInterp */
  if (in->beamInterp) {
    if (in->beamInterp->h_myArray) {
      CUDAFArrayZap(in->beamInterp->h_myArray);  /* This should get the device structures */
    } /* end CUDAFArray */
    if (in->beamInterp->d_myDesc) {
      ObitCUDAUtilFreeGPU((float*)in->beamInterp->d_myDesc);
    } /* end CUDAImageDesc */
    
    ObitCUDAUtilFreeHost((float*)in->beamInterp);  /* basic structure */
    } /* end zap CUDABeamInterp */
   #endif /* HAVE_GPU */
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

  /* Delete basic object NO already done
  ObitCUDAUtilFreeHost((float*)in); */
  
} /* end ObitGPUBeamInterpClear */
