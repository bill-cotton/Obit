/* $Id$         */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2025                                          */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#if HAVE_AVX512==1
#include <immintrin.h>
# define ALIGN32_BEG
# define ALIGN32_END __attribute__((packed,aligned(32)))
/* Union allowing c interface */
typedef __m512  V16SF; // vector of 16 float (avx512)
typedef __m512i V16SI; // vector of 16 int   (avx512)
typedef __mmask16 MASK16; // vector of 16 mask  (avx512)
typedef ALIGN32_BEG union {
  float f[16];
  int   i[16];
  V16SF   v;
} ALIGN32_END CV16SF;
typedef ALIGN32_BEG union {
  float f[16];
  int   i[16];
  V16SI   v;
} ALIGN32_END IV16SF;
#endif
#if HAVE_AVX==1
#include <immintrin.h>
# define ALIGN32_BEG
# define ALIGN32_END __attribute__((packed,aligned(32)))
/* Union allowing c interface */
typedef __m256  V8SF; // vector of 8 float (avx)
typedef __m256i V8SI; // vector of 8 int   (avx)
typedef ALIGN32_BEG union {
  float f[8];
  int   i[8];
  V8SF   v;
} ALIGN32_END CV8SF;
typedef ALIGN32_BEG union {
  float f[8];
  int   i[8];
  V8SI   v;
} ALIGN32_END IV8SF;
#endif
/** SSE implementation 4 floats in parallel */
#if HAVE_SSE==1
#include <xmmintrin.h>
# define ALIGN32_BEG
# define ALIGN32_END __attribute__((packed,aligned(32)))
typedef __m128 V4SF;  // vector of 4 float
typedef ALIGN32_BEG union {
  float f[4];
  int   i[4];
  V4SF  v;
} ALIGN32_END CV4SF;
#endif

#include "ObitCArray.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCArray.c
 * ObitCArray class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitCArray";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/*---------------Private structures----------------*/
/* Threaded function argument */
typedef struct {
  /* ObitThread to use */
  ObitThread *thread;
  /* ObitCArray to work on */
  ObitCArray *in;
  /* Second ObitCArray to work on */
  ObitCArray *in2;
  /* Output ObitCArray */
  ObitCArray *out;
  /* First ObitFArray */
  ObitFArray *FA_1;
  /* Second ObitFArray */
  ObitFArray *FA_2;
  /* Third ObitFArray */
  ObitFArray *FA_3;
  /* Fourth ObitFArray */
  ObitFArray *FA_4;
  /* First element (1-rel) number */
  olong        first;
  /* Highest element (1-rel) number */
  olong        last;
  /* Function dependent arguments */
  ofloat *arg1, *arg2,* arg3, *arg4, *arg5, *arg6, *arg7;
  /* Return value */
  gfloat value;
  /* Return position */
  olong        pos[MAXFARRAYDIM];
  /* thread number  */
  olong        ithread;
} CAFuncArg;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitCArrayClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitCArrayClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitCArrayInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitCArrayClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitCArrayClassInfoDefFn (gpointer inClass);

/** Private: Threaded Add */
static gpointer ThreadCAAdd (gpointer arg);

/** Private: Threaded Sub */
static gpointer ThreadCASub (gpointer arg);

/** Private: Threaded Multiply */
static gpointer ThreadCAMul (gpointer arg);

/** Private: Threaded Divide */
static gpointer ThreadCADiv (gpointer arg);

/** Private: Threaded Fill */
static gpointer ThreadCAFill (gpointer arg);

/** Private: Threaded Amplitude */
static gpointer ThreadCAAmp (gpointer arg);

/** Private: Threaded Phase */
static gpointer ThreadCAPhase (gpointer arg);

/** Private: Threaded MulAccum */
static gpointer ThreadCAMulAccum (gpointer arg);

/** Private: Threaded SMulAccum */
static gpointer ThreadCASMulAccum (gpointer arg);

/*-------------- Threading routines ----------------------------------*/
/** Private: Make Threaded args */
static olong MakeCAFuncArgs (ObitThread *thread, ObitCArray *in,
			     ObitCArray *in2, ObitCArray *out, 
			     ObitFArray *FA_1, ObitFArray *FA_2,
			     ObitFArray *FA_3, ObitFArray *FA_4,
			     olong larg1, olong larg2, olong larg3, 
			     olong larg4, olong larg5, 
			     olong larg6, olong larg7, 
			     CAFuncArg ***ThreadArgs);

/** Private: Delete Threaded args/free Thread pool */
static void KillCAFuncArgs (olong nargs, CAFuncArg **ThreadArgs);
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitCArray* newObitCArray (gchar* name)
{
  ObitCArray* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitCArrayClassInit();

  /* allocate/init structure */
  out = ObitMemAlloc0Name(sizeof(ObitCArray), "ObitCArray");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitCArrayInit((gpointer)out);

 return out;
} /* end newObitCArray */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitCArrayGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitCArrayClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitCArrayGetClass */

/**
 * Make a deep copy of an ObitCArray.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitCArray* ObitCArrayCopy  (ObitCArray *in, ObitCArray *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i;

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
    out = newObitCArray(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->info = ObitInfoListUnref(out->info);
  out->info = ObitInfoListCopy(in->info);

  /* arrays */
  out = ObitCArrayRealloc (out, in->ndim, in->naxis);
  /* copy data */
  for (i=0; i<2*in->arraySize; i++) out->array[i] = in->array[i];

  return out;
} /* end ObitCArrayCopy */

/**
 * Determine if the two input objects have a compatable geometry.
 * Must have same number of non degenerate dimensions and 
 * each dimension must be the same size.
 * \param in1  First object to test.
 * \param in2  Second object to test.
 * \return TRUE if compatable, else FALSE.
 */
gboolean ObitCArrayIsCompatable (ObitCArray *in1, ObitCArray *in2)
{
  olong i, ndim;

  /* test fails if either input is NULL */
  if (in1==NULL) return FALSE;
  if (in2==NULL) return FALSE;

   /* Same number of non degenerate dimensions? */
  ndim = in1->ndim;
  if (in1->ndim!=in2->ndim) {
    /* Don't bother if extra dimensions are length 1 */
    if (in1->ndim>in2->ndim) {  /* in1 has more */
      ndim = in2->ndim;
      for (i=in2->ndim; i<in1->ndim; i++)
	if (in1->naxis[i]>1) return FALSE;
    } else {  /* in2 has more */
      ndim = in1->ndim;
      for (i=in1->ndim; i<in2->ndim; i++)
	if (in2->naxis[i]>1) return FALSE;
    }
  }  /* end dimension check */

  /* check dimensions */
  for (i=0; i<ndim; i++) if (in1->naxis[i]!=in2->naxis[i]) return FALSE;

  return TRUE; /* must be OK */
} /* end ObitCArrayIsCompatable  */

/**
 * Creates an ObitCArray of a specified geometry.
 * \param name  An optional name for the object.
 * \param ndim  Number of dimensions desired, if <=0 data array not created.
 *              maximum value = 10.
 * \param naxis Dimensionality along each axis. NULL => don't create array.
 * \return the new object.
 */
ObitCArray* ObitCArrayCreate (gchar* name, olong ndim, olong *naxis)
{
  ObitCArray* out;
  olong i, size;

  /* Create basic structure */
  out = newObitCArray (name);

  /* create data array if wanted */
  if ((ndim<0) || (naxis==NULL)) return out;
  g_assert (ndim<11); /* sanity check */

  /* copy geometry */
  out->ndim = ndim;
  out->naxis = ObitMemAlloc0Name(ndim*sizeof(olong), "CArray naxis");
  size = 1;
  for (i=0; i<ndim; i++) {
    out->naxis[i] = MAX (1, naxis[i]);
    size *= out->naxis[i]; /* total size */
  }

  /* create array - add a bit extra, FFT seems to need it */
  out->array = ObitMemAlloc0Name(2*size*sizeof(ofloat) + 
				 2*out->naxis[0]*sizeof(ofloat),
				 "CArray array");
  out->arraySize = size;

  return out;
} /* end ObitCArrayCreate */

/**
 * Reallocate memory if needed, zero memory.
 * \param in    Object with structures to reallocate.
 * \param ndim  Number of dimensions desired, if <=0 data array not created.
 *              maximum value = 10.
 * \param naxis Dimensionality along each axis. NULL => don't create array.
 * \return the new object.
 */
ObitCArray* ObitCArrayRealloc (ObitCArray* in, olong ndim, olong *naxis)
{
  ObitCArray* out = in;
  olong i, size;

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (naxis != NULL);
  g_assert (ndim<11); /* sanity check */

  /* just deallocate old? */
  if ((ndim<0) || (naxis==NULL)) {
    /* just deallocate old */
    if (out->array) out->array = ObitMemFree(out->array); 
    out->arraySize = 0;
    if (out->naxis) out->naxis = ObitMemFree(out->naxis); 
    out->ndim = 0;
    return out;
  }
 
  /* different geometry? */
  if (out->ndim != ndim) {
    out->ndim = ndim;
    if (out->naxis!=NULL) out->naxis = ObitMemFree(out->naxis);
    out->naxis = ObitMemAlloc0Name(ndim*sizeof(olong), "CArray naxis");
  }

  /* set dimensions, find output size */
  size = 1;
  for (i=0; i<ndim; i++) {
    out->naxis[i] = MAX (1, naxis[i]);
    size *= out->naxis[i]; /* total size */
  }

  /* resize array if needed */
  if (size != out->arraySize) {
    out->array = 
      ObitMemRealloc(out->array, 2*size*sizeof(ofloat)+
		     2*out->naxis[0]*sizeof(ofloat));
    out->arraySize = size;
  }

  /* zero fill memory */
  memset (out->array, 0, 2*size*sizeof(ofloat));

  return out;
} /* end ObitCArrayRealloc */

/**
 * Calculate offset for a given pixel location and return pointer.
 * Subsequent data are stored in order of increasing dimension 
 * (rows, then columns...).
 * \param in      Object with data
 * \param pos     array of 0-rel pixel numbers on each axis
 * \return pointer to specified object; NULL if illegal pixel.
 */
ofloat*  ObitCArrayIndex (ObitCArray *in, olong *pos)
{
  ofloat *out = NULL;
  olong i, indx, previous;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (pos != NULL);
  g_assert (in->array != NULL);
  g_assert (in->naxis != NULL);

  /* Calculate offset */
  previous = 2;  /* size of previous dimension */
  indx = 0;
  for (i=0; i<in->ndim; i++) {
    /* check legality */
    if ((pos[i]<0) || (pos[i]>=in->naxis[i])) return NULL;
    indx += pos[i] * previous;
    previous *= in->naxis[i];
  }

  /* pointer to array */
  out = in->array + indx;

  return out;
} /* end ObitCArrayIndex  */

/**
 * Find maximum abs value pixel.
 * Return value and location in pos.
 * \param in      Object with data
 * \param pos     (out) array of 0-rel pixel numbers on each axis
 * \return maximum modulus value.
 */
ofloat ObitCArrayMaxAbs (ObitCArray *in, olong *pos)
{
  olong i, temp, maxCell;
  ofloat maxVal, *data, val, fblank = ObitMagicF();

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (pos != NULL);

  /* Loop over array */
  maxCell = -1;
  maxVal = -1.0E25;
  data = in->array;
  for (i=0; i<2*in->arraySize; i+=2) 
    {
      if ((data[i]!=fblank) && (data[i+1]!=fblank)) {
	val = data[i]*data[i] + data[i+1]*data[i+1];
	if ((val!=fblank) && (fabs(val)>maxVal)) {
	  maxCell = i;
	  maxVal  = val;
	}
      }
    }

  /* Convert cell to index */
  temp = maxCell;
  for (i=0; i<in->ndim; i++) {
    if (i==0) temp /= 2;  /*take account of it being complex */
    pos[i] = temp % in->naxis[i];
    temp = (temp - pos[i]) / in->naxis[i];
  }

  return sqrt(maxVal);
} /* end  ObitCArrayMaxAbs */

/**
 * Find minimum real or imaginary  value pixel.
 * Return value and location in pos.
 * \param in      Object with data
 * \param pos     (out) array of 0-rel pixel numbers on each axis
 * \return minimum real or imaginary  value.
 */
ofloat ObitCArrayMin (ObitCArray *in, olong *pos)
{
  olong i, temp, minCell;
  ofloat minVal, *data, val, fblank = ObitMagicF();

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (pos != NULL);

  /* Loop over array */
  minCell = -1;
  minVal = 1.0E25;
  data = in->array;
  for (i=0; i<2*in->arraySize; i+=2) 
    {
      if ((data[i]!=fblank) && (data[i+1]!=fblank)) {
	val = MIN (data[i], data[i+1]);
	if (val<minVal) {
	  minCell = i;
	  minVal  = val;
	} 
      }
    }

  /* Convert cell to index */
  temp = minCell;
  for (i=0; i<in->ndim; i++) {
    if (i==0) temp /= 2;  /*take account of it being complex */
    pos[i] = temp % in->naxis[i];
    temp = (temp - pos[i]) / in->naxis[i];
  }

  return minVal;
} /* end  ObitCArrayMin */

/**
 *  Negate each element of the array.
 * in = -in.
 * \param in Input object with data
 */
void ObitCArrayNeg (ObitCArray* in)
{
  olong i;

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<2*in->arraySize; i++) in->array[i] = -in->array[i];
} /* end  ObitCArrayNeg */

/**
 *  Conjugate each element of the array.
 * in = conjg(in).
 * \param in Input object with data
 */
void ObitCArrayConjg (ObitCArray* in)
{
  olong i;

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<2*in->arraySize; i += 2) in->array[i+1] = -in->array[i+1];
} /* end  ObitCArrayConjg */

/**
 * Fill the elements of an array with a complex scalar
 * in = cmpx
 * \param in      Input object with data
 * \param cmpx    Scalar value as  (real,imaginary)
 */
void ObitCArrayFillX (ObitCArray* in, ofloat cmpx[2])
{
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<2*in->arraySize; i += 2) {
    in->array[i]   = cmpx[0];
    in->array[i+1] = cmpx[1];
  }
} /* end ObitCArrayFillX */

/**
 * Fill the elements of an array with a complex scalar
 * in = cmpx
 * \param in      Input object with data
 * \param cmpx    Scalar value as  (real,imaginary)
 */
void ObitCArrayFill (ObitCArray* in, ofloat cmpx[2])
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  CAFuncArg **threadArgs;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  /* Initialize Threading */
  nThreads = MakeCAFuncArgs (in->thread, in, NULL, NULL, NULL, NULL, NULL, NULL, 2, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in->arraySize;
  /* At least 50,000 per thread */
  nTh = MAX (1, MIN((olong)(0.5+nElem/50000.),nThreads));
  nElemPerThread = nElem/nTh;
  if (nElem<50000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem;
    threadArgs[i]->last    = hiElem;
    threadArgs[i]->arg1[0] = cmpx[0];
    threadArgs[i]->arg1[1] = cmpx[1];
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  ObitThreadIterator (in->thread, nTh, 
		      (ObitThreadFunc)ThreadCAFill,
		      (gpointer**)threadArgs);

  /* Free local objects */
  KillCAFuncArgs(nThreads, threadArgs);
  
} /* end ObitCArrayFill */

/**
 *  Add a real scalar to each element of the array.
 * in = in + scalar
 * \param in      Input object with data
 * \param scalar  Scalar value
 */
void ObitCArraySAdd (ObitCArray* in, ofloat scalar)
{
  olong i;

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<2*in->arraySize; i += 2) {
    in->array[i] += scalar;
  }
} /* end ObitCArraySAdd */

/**
 *  Multiply each element of the array by a real scalar.
 * in = in * scalar
 * \param in      Input object with data
 * \param scalar  Scalar value
 */
void ObitCArraySMul (ObitCArray* in, ofloat scalar)
{
  olong i;

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<2*in->arraySize; i++) in->array[i]  *= scalar;
} /* end ObitCArraySMul */

/**
 *  Add a complex scalar to each element of the array.
 * in = in + scalar
 * \param in      Input object with data
 * \param scalar  Scalar value
 */
void ObitCArrayCSAdd (ObitCArray* in, ofloat scalar[2])
{
  olong i;
  ofloat sr, si, fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  sr = scalar[0];
  si = scalar[1];
  for (i=0; i<2*in->arraySize; i+=2) {
    if ((in->array[i]!=fblank) && (in->array[i+1]!=fblank)) {
      in->array[i]   += sr;
      in->array[i+1] += si;
    } else {
      in->array[i]   = fblank;
      in->array[i+1] = fblank;
    }
  }
} /* end ObitCArrayCSAdd */

/**
 *  Multiply each element of the array by a complex scalar.
 * in = in * scalar
 * \param in      Input object with data
 * \param scalar  Scalar value
 */
void ObitCArrayCSMul (ObitCArray* in, ofloat scalar[2])
{
  olong i;
  ofloat tr1, ti1, tr2, ti2, fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  tr2 = scalar[0];
  ti2 = scalar[1];
  for (i=0; i<2*in->arraySize; i+=2) {
    tr1 = in->array[i];
    ti1 = in->array[i+1];
    if ((tr1!=fblank) && (ti1!=fblank)) {
      in->array[i]   = tr1*tr2 - ti1*ti2;
      in->array[i+1] = ti1*tr2 + tr1*ti2;
    } else {
      in->array[i]   = fblank;
      in->array[i+1] = fblank;
    }
  }
} /* end ObitCArrayCSMul */

/**
 *  Add corresponding elements of two arrays.
 *  out = in1 + in2,  if either is blanked the result is blanked
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitCArrayAdd (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  CAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitCArrayIsCompatable(in1, in2));
  g_assert (ObitCArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeCAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work - pretend floats with 2 entries per cell*/
  nElem = 2*in1->arraySize;
  /* At least 100,000 per thread */
  nTh = MAX (1, MIN((olong)(0.5+nElem/100000.),nThreads));
  nElemPerThread = nElem/nTh;
  if (nElem<1000000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem;
    threadArgs[i]->last    = hiElem;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  ObitThreadIterator (in1->thread, nTh, 
		      (ObitThreadFunc)ThreadCAAdd,
		      (gpointer**)threadArgs);
  
  /* Free local objects */
  KillCAFuncArgs(nThreads, threadArgs);
  
} /* end ObitCArrayAdd */

/**
 *  Subtract corresponding elements of the arrays.
 *  out = in1 - in2, if either is blanked the result is blanked
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitCArraySub (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  CAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitCArrayIsCompatable(in1, in2));
  g_assert (ObitCArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeCAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work  - pretend floats with 2 entries per cell */
  nElem = 2*in1->arraySize;
  /* At least 100,000 per thread */
  nTh = MAX (1, MIN((olong)(0.5+nElem/100000.),nThreads));
  nElemPerThread = nElem/nTh;
  if (nElem<100000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem;
    threadArgs[i]->last    = hiElem;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  ObitThreadIterator (in1->thread, nTh, 
		      (ObitThreadFunc)ThreadCASub,
		      (gpointer**)threadArgs);

  /* Free local objects */
  KillCAFuncArgs(nThreads, threadArgs);
  
 } /* end ObitCArraySub */

/**
 *  Multiply corresponding elements of the arrays.
 *  out = in1 * in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitCArrayMul (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  CAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitCArrayIsCompatable(in1, in2));
  g_assert (ObitCArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeCAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in1->arraySize;
  /* At least 100,000 per thread */
  nTh = MAX (1, MIN((olong)(0.5+nElem/100000.),nThreads));
  nElemPerThread = nElem/nTh;
  if (nElem<100000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem;
    threadArgs[i]->last    = hiElem;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  ObitThreadIterator (in1->thread, nTh, 
		      (ObitThreadFunc)ThreadCAMul,
		      (gpointer**)threadArgs);
  
  /* Free local objects */
  KillCAFuncArgs(nThreads, threadArgs);
  
} /* end ObitCArrayMul */

/**
 *  Divide corresponding elements of the arrays.
 *  out = in1 / in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitCArrayDiv (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  CAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitCArrayIsCompatable(in1, in2));
  g_assert (ObitCArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeCAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in1->arraySize;
  /* At least 200,000 per thread */
  nTh = MAX (1, MIN((olong)(0.5+nElem/200000.),nThreads));
  nElemPerThread = nElem/nTh;
  if (nElem<100000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem;
    threadArgs[i]->last    = hiElem;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  ObitThreadIterator (in1->thread, nTh, 
		      (ObitThreadFunc)ThreadCADiv,
		      (gpointer**)threadArgs);

  /* Free local objects */
  KillCAFuncArgs(nThreads, threadArgs);
  
} /* end ObitCArrayDiv */

/**
 *  Add corresponding elements of the array.
 *  out = in1 + in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitCArrayAddX (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
{
  olong i;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitCArrayIsCompatable(in1, in2));
  g_assert (ObitCArrayIsCompatable(in1, out));

  for (i=0; i<2*in1->arraySize; i++)
    out->array[i] = in1->array[i] + in2->array[i];
} /* end ObitCArrayAddX */

/**
 *  Subtract corresponding elements of the arrays.
 *  out = in1 - in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitCArraySubX (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
{
  olong i;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitCArrayIsCompatable(in1, in2));
  g_assert (ObitCArrayIsCompatable(in1, out));

  for (i=0; i<2*in1->arraySize; i++) 
    out->array[i] = in1->array[i] - in2->array[i];
 } /* end ObitCArraySubX */

/**
 *  Multiply corresponding elements of the arrays.
 *  Output may be one of the inputs
 *  out = in1 * in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitCArrayMulX (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
{
  olong i;
  ofloat tr1, ti1, tr2, ti2;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitCArrayIsCompatable(in1, in2));
  g_assert (ObitCArrayIsCompatable(in1, out));

  for (i=0; i<2*in1->arraySize; i += 2) {
    tr1 = in1->array[i];
    ti1 = in1->array[i+1];
    tr2 = in2->array[i];
    ti2 = in2->array[i+1];
    out->array[i]   = tr1*tr2 - ti1*ti2;
    out->array[i+1] = ti1*tr2 + tr1*ti2;
  }
} /* end ObitCArrayMulX */

/**
 *  Divide corresponding elements of the arrays.
 *  out = in1 / in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitCArrayDivX (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
{
  olong i;
  ofloat denom;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitCArrayIsCompatable(in1, in2));
  g_assert (ObitCArrayIsCompatable(in1, out));

  for (i=0; i<2*in1->arraySize; i += 2) {
    denom = (in2->array[i]*in2->array[i]+ 
	     in2->array[i+1]*in2->array[i+1]);
    out->array[i]   = (in1->array[i] * in2->array[i] +
      in1->array[i+1] * in2->array[i+1]) / denom;
    out->array[i+1] = (in1->array[i+1] * in2->array[i] -
      in1->array[i] * in2->array[i+1]) / denom;
  }
} /* end ObitCArrayDivX */

/* ----------------  CArray - FArray functions  -------------------  */

/**
 *   Make an FArray with same geometry as a CArray
 *  out = copy((in))
 * \param in  Input object with data
 * \return the new object.
 */
ObitFArray* ObitCArrayMakeF  (ObitCArray *in)
{
  gchar *outName;
  ObitFArray* out = NULL;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Create  */
  /* derive object name */
  outName = g_strconcat ("Copy: ",in->name, NULL);
  out = ObitFArrayCreate(outName, in->ndim, in->naxis);
  g_free(outName);

  return out;
}  /* end ObitCArrayMakeF */

/**
 *  Make a CArray with same geometry as a FArray
 *  out = copy(in)
 * \param in  Input object with data
 * \return the new object.
 */
ObitCArray* ObitCArrayMakeC  (ObitFArray *in)
{
  gchar *outName;
  ObitCArray* out = NULL;

  /* error checks */
  g_assert (ObitFArrayIsA(in));

  /* Create  */
  /* derive object name */
  outName = g_strconcat ("Copy: ",in->name, NULL);
  out = ObitCArrayCreate(outName, in->ndim, in->naxis);
  g_free(outName);

  return out;
}  /* end ObitCArrayMakeC */

/**
 *  Test if the geometries of an FArray and a CArray are compatible
 *  true iff FArray compatable with CArray
 * \param Cin  Input CArray
 * \param Fin  Input FArray
 * \return TRUE if compatable, else FALSE.
 */
gboolean ObitCArrayIsFCompatable  (ObitCArray *Cin, ObitFArray *Fin)
{
  olong i, ndim;

  /* test fails if either input is NULL */
  if (Fin==NULL) return FALSE;
  if (Cin==NULL) return FALSE;

  /* Same number of non degenerate dimensions? */
  ndim = Cin->ndim;
  if (Cin->ndim!=Fin->ndim) {
    /* Don't bother if extra dimensions are length 1 */
    if (Cin->ndim>Fin->ndim) {  /* Cin has more */
      ndim = Fin->ndim;
      for (i=Fin->ndim; i<Cin->ndim; i++)
	if (Cin->naxis[i]>1) return FALSE;
    } else {  /* Fin has more */
      ndim = Cin->ndim;
      for (i=Cin->ndim; i<Fin->ndim; i++)
	if (Fin->naxis[i]>1) return FALSE;
    }
  }  /* end dimension check */

  /* Same number of dimensions? */
  if (Cin->ndim!=Fin->ndim) return FALSE;

  /* check dimensions */
  for (i=0; i<ndim; i++) if (Fin->naxis[i]!=Cin->naxis[i]) return FALSE;

  return TRUE; /* must be OK */
}  /* end ObitCArrayIsFCompatable */

/**
 *  Multiply the elements of a CArray by the elements of an FArray
 * Output may be one of the inputs.
 *  out = Cin * Fin
 * \param Cin  Input CArray
 * \param Fin  Input FArray
 * \param out  Output CArray
 */
void ObitCArrayFMul (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out)
{
  olong i, j;
  ofloat trc, tic, tr;

  /* error checks */
  g_assert (ObitCArrayIsA(Cin));
  g_assert (ObitFArrayIsA(Fin));
  g_assert (ObitCArrayIsCompatable(Cin, out));
  g_assert (ObitCArrayIsFCompatable(Cin, Fin));

  /* Multiply */
  j = 0;
  for (i=0; i<2*Cin->arraySize; i += 2) {
    trc = Cin->array[i];
    tic = Cin->array[i+1];
    tr  = Fin->array[j];
    out->array[i]   = trc * tr;
    out->array[i+1] = tic * tr;
    j++;
  }
}  /* end ObitCArrayFMul */

/**
 *  Divide the elements of a CArray by the elements of an FArray
 *  Output may be input
 *  out = Cin / Fin
 * \param Cin  Input CArray
 * \param Fin  Input FArray
 * \param out  Output CArray
 */
void ObitCArrayFDiv (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out)
{
  olong i, j;
  ofloat trc, tic, tr;

  /* error checks */
  g_assert (ObitCArrayIsA(Cin));
  g_assert (ObitFArrayIsA(Fin));
  g_assert (ObitCArrayIsCompatable(Cin, out));
  g_assert (ObitCArrayIsFCompatable(Cin, Fin));

  /* Divide */
  j = 0;
  for (i=0; i<2*Cin->arraySize; i += 2) {
    trc = Cin->array[i];
    tic = Cin->array[i+1];
    tr  = Fin->array[j];
    out->array[i]   = trc / tr;
    out->array[i+1] = tic / tr;
    j++;
  }
}  /* end ObitCArrayFDiv */

/**
 *  Add An FArray to the real elements of a CArray
 *  out = out + in
 * \param Cin  Input CArray
 * \param Fin  Input FArray
 * \param out  Output CArray
 */
void ObitCArrayFAdd (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out)
{
  olong i, j;

  /* error checks */
  g_assert (ObitCArrayIsA(Cin));
  g_assert (ObitFArrayIsA(Fin));
  g_assert (ObitCArrayIsCompatable(Cin, out));
  g_assert (ObitCArrayIsFCompatable(Cin, Fin));

  /* Add */
  j = 0;
  for (i=0; i<2*Cin->arraySize; i += 2) {
    out->array[i]   = Cin->array[i] + Fin->array[j];
    out->array[i+1] = Cin->array[i+1];
    j++;
  }
}  /* end ObitCArrayFAdd */

/**
 *  Rotate phase of the elements of a CArray by the elements of an FArray
 *  Output may be one of the inputs.
 *  out = Cin * exp(i*Fin)
 * \param Cin  Input CArray
 * \param Fin  Input FArray
 * \param out  Output CArray
 */
void ObitCArrayFRot (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out)
{
  olong i, j;
  ofloat trc, tic, tr, ti, fblank = ObitMagicF();

  /* error checks */
  g_assert (ObitCArrayIsA(Cin));
  g_assert (ObitFArrayIsA(Fin));
  g_assert (ObitCArrayIsCompatable(Cin, out));
  g_assert (ObitCArrayIsFCompatable(Cin, Fin));

  /* Rotate */
  j = 0;
  for (i=0; i<2*Cin->arraySize; i += 2) {
    trc = Cin->array[i];
    tic = Cin->array[i+1];
    if ((trc!=fblank) && (tic!=fblank) && (Fin->array[j]!=fblank)) {
      tr  = cos(Fin->array[j]);
      ti  = sin(Fin->array[j]);
      out->array[i]   = trc*tr - tic*ti;
      out->array[i+1] = tic*tr + trc*ti;
    } else {
      out->array[i]   = fblank;
      out->array[i+1] = fblank;
    }
    j++;
  }
}  /* end ObitCArrayFRot */

/**
 *  Combine two FArrays into a CArray
 *  out = complex(Fin1, Fin2)
 * \param Fin1 Input FArray for real part
 * \param Fin2 Input FArray for imaginary part
 * \param out  Output CArray
 */
void ObitCArrayComplex (ObitFArray* Fin1, ObitFArray* Fin2, ObitCArray* out)
{
  olong i, j;

  /* error checks */
  g_assert (ObitFArrayIsA(Fin1));
  g_assert (ObitFArrayIsA(Fin2));
  g_assert (ObitCArrayIsA(out));
  g_assert (ObitFArrayIsCompatable(Fin1, Fin2));
  g_assert (ObitCArrayIsFCompatable(out, Fin1));

  /* Combine */
  j = 0;
  for (i=0; i<2*Fin1->arraySize; i += 2) {
    out->array[i]   = Fin1->array[j];
    out->array[i+1] = Fin2->array[j];
    j++;
  }
}  /* end ObitCArrayComplex */

/**
 *  Return real part of CArray
 *  out = real(in)
 * \param in  Input CArray
 * \param out  Output CArray
 */
void ObitCArrayReal (ObitCArray* in, ObitFArray* out)
{
  olong i, j;

  /* error checks */
  g_assert (ObitCArrayIsA(in));
  g_assert (ObitFArrayIsA(out));
  g_assert (ObitCArrayIsFCompatable(in, out));

  /* Extract Real part */
  j = 0;
  for (i=0; i<2*in->arraySize; i += 2) {
    out->array[j]   = in->array[i];
    j++;
  }
}  /* end ObitCArrayReal */

/**
 *  Return imaginary part of CArray
 *  out = imag(in)
 * \param in  Input CArray
 * \param out Output FArray
 */
void ObitCArrayImag (ObitCArray* in, ObitFArray* out)
{
  olong i, j;

  /* error checks */
  g_assert (ObitCArrayIsA(in));
  g_assert (ObitFArrayIsA(out));
  g_assert (ObitCArrayIsFCompatable(in, out));

  /* Extract Imaginary part */
  j = 0;
  for (i=0; i<2*in->arraySize; i += 2) {
    out->array[j]   = in->array[i+1];
    j++;
  }
}  /* end ObitCArrayImag */

/**
 *  Return amplitude of elements of a CArray
 *  out = sqrt(real(in)^2 + imag(in)^2)
 * \param in  Input CArray
 * \param out Output FArray
 */
void ObitCArrayAmp (ObitCArray* in, ObitFArray* out)
{
  olong i, nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  CAFuncArg **threadArgs;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  /* Initialize Threading */
  nThreads = MakeCAFuncArgs (in->thread, in, NULL, NULL, out, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in->arraySize;
  /* At least 10,000 per thread */
  nTh = MAX (1, MIN((olong)(0.5+nElem/10000.),nThreads));
  nElemPerThread = nElem/nTh;
  if (nElem<10000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem;
    threadArgs[i]->last    = hiElem;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  ObitThreadIterator (in->thread, nTh, 
		      (ObitThreadFunc)ThreadCAAmp,
		      (gpointer**)threadArgs);

  /* Free local objects */
  KillCAFuncArgs(nThreads, threadArgs);
  
}  /* end ObitCArrayAmp */

/**
 *  Return amplitude^2 of elements of a CArray
 * Uses intrinsics but is not threaded on purpose.
 *  out = real(in)^2 + imag(in)^2
 * Supports blanking
 * \param in  Input CArray
 * \param out Output FArray
 */
void ObitCArrayAmp2 (ObitCArray* in, ObitFArray* out)
{
  olong i, ilast;
  ofloat fblank = ObitMagicF();
  ofloat *iArr = in->array;
#if HAVE_AVX512==1
  CV16SF vr, vi, vb, vamp;
  IV16SF vindxr, vindxi;
  MASK16 msk1, msk2;
#elif HAVE_AVX==1
  CV8SF  vb, vm1, vm2, tv1, tv2;
#endif

  /* error checks */
  g_assert (ObitCArrayIsA(in));
  g_assert (ObitFArrayIsA(out));
  g_assert (ObitCArrayIsFCompatable(in, out));

  /* Extract Amplitude^2 */
#if HAVE_AVX512==1  /* Vector length 16 complex - uses gather */
  vb.v     = _mm512_set1_ps(fblank);  /* vector of blanks */
  vindxr.v = _mm512_set_epi32 (30,28,26,24,22,20,18,16,14,12,10,8,6,4,2,0); /* Reals*/
  vindxi.v = _mm512_set_epi32 (31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1); /* Imags */
  /* Do blocks of 16 complex as vector */
  for (i=0; i<in->arraySize-16; i+=16) {
    vr.v  = _mm512_i32gather_ps (vindxr.v, &iArr[2*i], 4); /* Load Reals */
    msk1  = _mm512_cmp_ps_mask(vr.v, vb.v, _CMP_EQ_OQ);         /* find blanks */
    vi.v  = _mm512_i32gather_ps (vindxi.v, &iArr[2*i], 4); /* Load Imaginaries */
    msk2  = _mm512_cmp_ps_mask(vi.v, vb.v, _CMP_EQ_OQ);         /* find blanks */
    vr.v  = _mm512_mul_ps (vr.v, vr.v);                         /* Real squared */
    vi.v  = _mm512_mul_ps (vi.v, vi.v);                         /* Imag squared */
    vamp.v  = _mm512_add_ps(vr.v, vi.v);                        /* R^2 + i^2 */
    vamp.v  = _mm512_mask_blend_ps(msk1,vamp.v,vb.v);           /* replace blanks real */
    vamp.v  = _mm512_mask_blend_ps(msk2,vamp.v,vb.v);           /* replace blanks imag */
    _mm512_storeu_ps(&out->array[i], vamp.v);                   /* Save */
  } /* end vector loop */
  ilast = i; /* How far did I get? */
#elif HAVE_AVX==1  /* AVX Vector (8 complex) */ 
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  for (i=0; i<in->arraySize-8; i+=8) {
    /* Gather the hard way - reals then imag */
    tv1.v = _mm256_set_ps(iArr[14+i*2],iArr[12+i*2],iArr[10+i*2],iArr[8+i*2],
			  iArr[6+i*2],iArr[4+i*2],iArr[2+i*2],iArr[i*2]);
    vm1.v = _mm256_cmp_ps(tv1.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    tv1.v = _mm256_mul_ps (tv1.v, tv1.v);           /* real squared */
    tv2.v = _mm256_set_ps(iArr[15+i*2],iArr[13+i*2],iArr[11+i*2],iArr[9+i*2],
			  iArr[7+i*2],iArr[5+i*2],iArr[3+i*2],iArr[1+i*2]);
    vm2.v = _mm256_cmp_ps(tv2.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    tv2.v = _mm256_mul_ps (tv2.v, tv2.v);           /* imaginary squared */
    tv2.v = _mm256_add_ps (tv1.v, tv2.v);           /* now amplitude squared */
    tv2.v = _mm256_blendv_ps(tv2.v,vb.v, vm1.v);    /* replace blanks 1 */
    tv2.v = _mm256_blendv_ps(tv2.v,vb.v, vm2.v);    /* replace blanks 2 */
    _mm256_storeu_ps(&out->array[i], tv2.v);        /* Save */
  } /* end vector loop */
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  for (i=ilast; i<in->arraySize; i++) {
    if ((iArr[i*2]!=fblank) && (iArr[1+i*2]!=fblank)) {
      out->array[i] = iArr[i*2]*iArr[i*2]+iArr[1+i*2]*iArr[1+i*2];
    } else {
      out->array[i] = fblank;
    }
  } /* end scalar loop */
}  /* end ObitCArrayAmp2 */

/**
 *  Return phase (radians) of elements of a CArray
 *  out = atan2(imag(in), real(in))
 * \param in  Input CArray
 * \param out Output FArray
 */
void ObitCArrayPhase (ObitCArray* in, ObitFArray* out)
{
  olong i, nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  CAFuncArg **threadArgs;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  /* Initialize Threading */
  nThreads = MakeCAFuncArgs (in->thread, in, NULL, NULL, out, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in->arraySize;
  /* At least 10,000 per thread */
  nTh = MAX (1, MIN((olong)(0.5+nElem/10000.),nThreads));
  nElemPerThread = nElem/nTh;
  if (nElem<10000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem;
    threadArgs[i]->last    = hiElem;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  ObitThreadIterator (in->thread, nTh, 
		      (ObitThreadFunc)ThreadCAPhase,
		      (gpointer**)threadArgs);

  /* Free local objects */
  KillCAFuncArgs(nThreads, threadArgs);
  
}  /* end ObitCArrayPhase */

/**
 *  Return phase (radians) of elements of a CArray
 *  out = atan2(imag(in), real(in))
 * \param in  Input CArray
 * \param out Output FArray
 */
void ObitCArrayPhaseX (ObitCArray* in, ObitFArray* out)
{
  olong i, j;
  ofloat fblank = ObitMagicF();

  /* error checks */
  g_assert (ObitCArrayIsA(in));
  g_assert (ObitFArrayIsA(out));
  g_assert (ObitCArrayIsFCompatable(in, out));

  /* Extract Phase */
  j = 0;
  for (i=0; i<2*in->arraySize; i += 2) {
    if ((in->array[i]!=fblank) && (in->array[i+1]!=fblank)) {
      out->array[j] = atan2(in->array[i+1], in->array[i]+1.0e-20);
    } else {
      out->array[j] = fblank;
    }
    j++;
  }
}  /* end ObitCArrayPhaseX */

/**
 * In-place rearrangement of a half plane center-at-the edges array to 
 * center at the center, or the other way around.
 * The first and second halves of each column are swaped.
 * This is needed for the peculiar order of FFTs.
 * This uses the FFTW convention that half plane complex arrays
 * have the "short" (i.e. n/2+1) dimension first.  
 * \param in   Half plane complex 2D array to reorder
 * \return the new object.
 */
void ObitCArray2DCenter (ObitCArray* in)
{
  olong i, j, pos[2], nx, ny;
  ofloat temp, *inp, *outp;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* ... and God sayeth unto FFT: 
     "The middle shall be first and the last shall be the middle." */

  /* swap top and bottom halves */
  nx = in->naxis[0];
  ny = in->naxis[1];
  pos[0] = 0; pos[1] = 0;
  inp = ObitCArrayIndex (in, pos);
  pos[0] = 0; pos[1] = ny/2;
  outp = ObitCArrayIndex (in, pos);

  for (j=0; j<ny/2; j++) {
    for (i=0; i<2*nx; i += 2) {
      temp = outp[i];
      outp[i] = inp[i];
      inp[i] = temp;
      temp = outp[i+1];
      outp[i+1] = inp[i+1];
      inp[i+1] = temp;
    }
    /* next row */
    inp += nx*2;
    outp += nx*2;
  }
} /* end ObitCArray2DCenter */

/**
 * In-place rearrangement of a full plane center-at-the edges array to 
 * center at the center, or the other way around.
 * The first and second halves of each column are swaped.
 * This is needed for the peculiar order of FFTs.
 * This uses the FFTW convention that half plane complex arrays
 * have the "short" (i.e. n/2+1) dimension first.  
 * \param in   Full plane complex 2D array to reorder
 * \return the new object.
 */
void ObitCArray2DCenterFull (ObitCArray* in)
{
  olong i, j, pos[2], nx, ny;
  ofloat temp, *inp, *outp;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* ... and God sayeth unto FFT: 
     "The middle shall be first and the last shall be the middle." */

  /* swap first and fourth quadrents */
  nx = in->naxis[0];
  ny = in->naxis[1];
  pos[0] = 0; pos[1] = 0;
  inp = ObitCArrayIndex (in, pos);
  pos[0] = nx/2; pos[1] = ny/2;
  outp = ObitCArrayIndex (in, pos);

  for (j=0; j<ny/2; j++) {
    for (i=0; i<nx; i+=2) {
      temp = outp[i];
      outp[i] = inp[i];
      inp[i] = temp;
      temp = outp[i+1];
      outp[i+1] = inp[i+1];
      inp[i+1] = temp;
    }
    /* next rwo */
    inp  += nx*2;
    outp += nx*2;
  }

  /* swap second and third quadrents */
  nx = in->naxis[0];
  ny = in->naxis[1];
  pos[0] = nx/2; pos[1] = 0;
  inp = ObitCArrayIndex (in, pos);
  pos[0] = 0; pos[1] = ny/2;
  outp = ObitCArrayIndex (in, pos);

  for (j=0; j<ny/2; j++) {
    for (i=0; i<nx; i+=2) {
      temp = outp[i];
      outp[i] = inp[i];
      inp[i] = temp;
      temp = outp[i+1];
      outp[i+1] = inp[i+1];
      inp[i+1] = temp;
    }
    /* next row */
    inp  += nx*2;
    outp += nx*2;
  }
} /* end ObitCArray2DCenterFull */

/**
 * Create a new ObitCArray which is a copy of the input array
 * with numConjCol conjugate columns added.
 * Only does 2D half plane complex images.
 * This uses the FFTW convention that half plane complex arrays
 * have the "short" (i.e. n/2+1) dimension first.  
 * Thus rows must have added columns.
 * \param in          Input array
 * \param numConjCol  How many conjugate columns to add.
 * \return the new object.
 */
ObitCArray* ObitCArrayAddConjg (ObitCArray* in, olong numConjCol)
{
  olong ndim, i1, i2;
  olong naxis[MAXFARRAYDIM], ipos[2], opos[2];
  ObitCArray *out=NULL;
  gchar *outName;
  ofloat *inp, *outp;
  
  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->ndim == 2);  /* Only 2D */
  g_assert (numConjCol>0);
  g_assert (numConjCol<in->naxis[1]);
  
  /* Get size */
  ndim = in->ndim;
  naxis[0] = in->naxis[0] + numConjCol; naxis[1] = in->naxis[1];
  
  /* derive object name */
  outName = g_strconcat ("ConjugateAdded: ",in->name, NULL);
  
  /* Create output */
  out = ObitCArrayCreate (outName, ndim, naxis);
  g_free(outName); outName = NULL;
  
  /* Copy input array */
  opos[0] = numConjCol;
  ipos[0] = 0;
  opos[1] = 0;
  ipos[1] = 0;
   
  /* Array pointers */
  inp  = ObitCArrayIndex (in,  ipos);
  outp = ObitCArrayIndex (out, opos);
  
  /* Copy in to end of output rows */
  for (i2=0; i2<in->naxis[1]; i2++) { /* second dimension*/
    for (i1=0; i1<2*in->naxis[0]; i1++) outp[i1] = inp[i1];
    inp  += 2*in->naxis[0];
    outp += 2*out->naxis[0];
  } /* end loop over second dimension */

  /* Add conjugate columns,  Note: pixels 0, n/2+1 stay put */
   for (i2=1; i2<in->naxis[1]; i2++) { /*second dimension*/

    /* Get pointers */
    ipos[0] = 1; ipos[1] = i2; 
    inp = ObitCArrayIndex(in, ipos);
    opos[0] = numConjCol-1; opos[1] = in->naxis[1]-i2; 
    outp = ObitCArrayIndex(out, opos);
 
    for (i1=0; i1<2*numConjCol; i1 += 2) {
      /* Loop down row, flipping, conjugating */
      outp[-i1]   =  inp[i1];
      outp[-i1+1] = -inp[i1+1];
    }
  } /* end loop over second dimension */

  /* Now 0th row */
  ipos[0] = numConjCol+1; ipos[1] = 0; 
  inp = ObitCArrayIndex(in, ipos);
  opos[0] = numConjCol-1; opos[1] = 0; 
  outp = ObitCArrayIndex(out, opos);
 
  for (i1=0; i1<2*numConjCol; i1 += 2) {
    /* Loop down row, flipping, conjugating */
    outp[-i1]   =  inp[i1];
    outp[-i1+1] = -inp[i1+1];
  }
  return out;
} /* end ObitCArrayAddConjg */

/**
 *  Form complex from two FArrays, multiply by complex scalar and 
 *  complex accumulate, replace blanks w/ 0
 *  ObitCArrayMakeCAFuncArgs to initialize threading with (elem)
 * \param Fin_r   Input FArray for real part      (FA_1)
 * \param Fin_i   Input FArray for imaginary part (FA_2)
 * \param cscalar Complex Scalar value [r,i]      (arg1)
 * \param Accum   CArray Accumulator              (out)
 */
void ObitCArraySMulAccum (ObitFArray* Fin_r, ObitFArray* Fin_i, 
			  ofloat cscalar[2], ObitCArray* Accum) {
  olong nThreads;
  CAFuncArg **threadArgs;
 
  /* error checks */
  g_assert (ObitFArrayIsA(Fin_r));
  g_assert (ObitFArrayIsA(Fin_i));
  g_assert (ObitCArrayIsA(Accum));
  g_assert (ObitFArrayIsCompatable(Fin_r, Fin_i));
  g_assert (ObitCArrayIsFCompatable(Accum, Fin_i));
  
  /* Initialize Threading */
  nThreads = MakeCAFuncArgs (Fin_r->thread, NULL, NULL, Accum, Fin_r, Fin_i, NULL, NULL, 2, 0, 0, 0, 0, 0, 0,
			     &threadArgs);

  /* Call routine to do the work */
  ObitCArraySMulAccumTh (Fin_r, Fin_i, cscalar, Accum, nThreads, (gpointer)threadArgs);

  /* Free local objects */
  KillCAFuncArgs(nThreads, threadArgs);
} /* end  ObitCArraySMulAccum */

/**
 *  Form complex from two FArrays, multiply by complex scalar and 
 *  complex accumulate - threading initialized elsewhere.
 *  Uses thread on threadArgs[0]->thread which should be used by 
 * Blanking NOT supported, replace with 0
 *  MakeCAFuncArgs to initialize threading (elem)
 * \param Fin_r   Input FArray for real part,     (FA_1)
 * \param Fin_i   Input FArray for imaginary part (FA_2)
 * \param cscalar Complex Scalar value [r,i]      (arg1)
 * \param Accum   CArray Accumulator              (out)
 * \param nThreads   Number of threads allocated
 * \param threadArgs Thread arguments
 */
void ObitCArraySMulAccumTh (ObitFArray* Fin_r, ObitFArray* Fin_i, 
			    ofloat cscalar[2], ObitCArray* Accum,
			    olong nThreads, gpointer inArgs) {
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread;
  CAFuncArg **threadArgs = (CAFuncArg**)inArgs;
 
  /* Divide up work */
  nElem = Fin_r->arraySize;
  /* At least 50,000 per thread */
  nTh = MAX (1, MIN((olong)(0.5+nElem/50000.),nThreads));
  nElemPerThread = nElem/nTh;
  if (nElem<100000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem;
    threadArgs[i]->last    = hiElem;
    threadArgs[i]->arg1[0] = cscalar[0]; threadArgs[i]->arg1[1] = cscalar[1]; 
    /* Reset pointers */
    if (threadArgs[i]->out)  ObitCArrayUnref(threadArgs[i]->out);
    threadArgs[i]->out = ObitCArrayRef(Accum);
    if (threadArgs[i]->FA_1) ObitFArrayUnref(threadArgs[i]->FA_1);
    threadArgs[i]->FA_1 = ObitFArrayRef(Fin_r);
    if (threadArgs[i]->FA_2) ObitFArrayUnref(threadArgs[i]->FA_2);
    threadArgs[i]->FA_2 = ObitFArrayRef(Fin_i);
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  ObitThreadIterator (threadArgs[0]->thread, nTh, 
		      (ObitThreadFunc)ThreadCASMulAccum,
		      (gpointer**)threadArgs);

} /* end  ObitCArraySMulAccumTh */

/**
 *  Complex multiply two pairs of FArrays and accumulate
 * Accum += Fin1 * Fin2
 * Blanking NOT supported, replace with 0
 * \param Fin1_r   First FArray real part       (FA_1)
 * \param Fin1_i   First FArray imaginary part  (FA_2)
 * \param Fin2_r   Second FArray real part      (FA_3)
 * \param Fin2_i   Second FArray imaginary part (FA_4)
 * \param Accum    Complex Accumulator          (out)
 */
void ObitCArrayMulAccum (ObitFArray* Fin1_r, ObitFArray* Fin1_i, 
			 ObitFArray* Fin2_r, ObitFArray* Fin2_i, 
			 ObitCArray* Accum) {
  olong nThreads;
  CAFuncArg **threadArgs;

  /* error checks */
  g_assert (ObitFArrayIsCompatable(Fin1_r, Fin2_r));
  g_assert (ObitCArrayIsFCompatable(Accum, Fin1_r));
  
  /* Initialize Threading */
  nThreads = MakeCAFuncArgs (Fin1_r->thread, NULL, NULL, Accum, Fin1_r, Fin1_i, Fin2_r, Fin2_i,
			     2, 0, 0, 0, 0, 0, 0,
			     &threadArgs);

  /* Call routine to do the work */
  ObitCArrayMulAccumTh (Fin1_r, Fin1_i, Fin2_r, Fin2_i, Accum,  
			nThreads, (gpointer)threadArgs);

  /* Free local objects */
  KillCAFuncArgs(nThreads, threadArgs);
} /* end  ObitCArrayMulAccum */

/**
 * If elements in in are >= those in out1 then out1 is set to the value in in 
 * and the corresponding element in out2 is set to value, 
 * and the corresponding value in out3 is set to than in in2.
 * Blanking NOT supported
 * \param in     Input FArray 
 * \param out1   First output FArray
 * \param value  new value for out2
 * \param out2   Second FArray output whose elements are set to value when
 *               an element in out1 is updated.
 * \param in2    CArray whose elements are written to out3 on an update
 * \param out3   CArray output
 */
void ObitCArrayMaxSetValues (ObitFArray* in, ObitFArray* out1, ofloat value, 
			     ObitFArray* out2, ObitCArray* in2, ObitCArray* out3)
{
  olong i, ilast;
  ofloat *iArr   = in->array,  *oArr1 = out1->array, *oArr2 = out2->array;
  ofloat *iCArr  = in2->array, *oCArr = out3->array;
#if HAVE_AVX512==1
  CV16SF v1, v2, v3, vs;
  IV16SF vindxr, vindxi;
  MASK16 msk;
#elif HAVE_AVX==1
  olong j, k;
  CV8SF v1, v2, v3, vs, vm;
#endif
  
  /* error checks */
  g_assert (ObitFArrayIsCompatable(in, out1));
  g_assert (ObitFArrayIsCompatable(out1, out2));
  
#if HAVE_AVX512==1  /* AVX512 Vector */
  vs.v     = _mm512_set1_ps(value);  /* vector of value */
  vindxr.v = _mm512_set_epi32 (30,28,26,24,22,20,18,16,14,12,10,8,6,4,2,0); /* Reals*/
  vindxi.v = _mm512_set_epi32 (31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1); /* Imags */
  /* Do blocks of 16 as vector */
  for (i=0; i<in->arraySize-16; i+=16) {
    v1.v = _mm512_loadu_ps(&iArr[i]);            /* Load from in */
    v2.v = _mm512_loadu_ps(&oArr1[i]);           /* Load from out1 */
    msk  = _mm512_cmp_ps_mask(v1.v, v2.v, _CMP_GE_OQ); /* find new max */
    v2.v = _mm512_mask_blend_ps(msk,v2.v,v1.v);   /* Select out1 */
    _mm512_storeu_ps(&oArr1[i], v2.v);            /* Save out1 */
    v3.v = _mm512_loadu_ps(&oArr2[i]);            /* Load from out2 */
    v3.v = _mm512_mask_blend_ps(msk,v3.v,vs.v);   /* Select out2 */
    _mm512_storeu_ps(&oArr2[i], v3.v);            /* Save out2 */
     /* Select complex values using gather/scatter */
    v1.v = _mm512_i32gather_ps (vindxr.v, &iCArr[2*i], 4); /* Load Reals in */
    v2.v = _mm512_i32gather_ps (vindxr.v, &oCArr[2*i], 4); /* Load Reals out */
    v3.v = _mm512_mask_blend_ps(msk,v2.v,v1.v);            /* Select out3 real */
    _mm512_i32scatter_ps((void*)(&oCArr[i*2]), vindxr.v, v3.v, 4); /* Scatter reals */
    v1.v = _mm512_i32gather_ps (vindxi.v, &iCArr[2*i], 4); /* Load Imag in */
    v2.v = _mm512_i32gather_ps (vindxi.v, &oCArr[2*i], 4); /* Load Imag out */
    v3.v = _mm512_mask_blend_ps(msk,v2.v,v1.v);            /* Select out3 imag */
    _mm512_i32scatter_ps((void*)(&oCArr[i*2]), vindxi.v, v3.v, 4); /* Scatter imag */
  } /* end outer loop */
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* AVX Vector */
  vs.v   = _mm256_broadcast_ss(&value);  /* vector of value */
  /* Do blocks of 8 as vector */
  for (i=0; i<in->arraySize-8; i+=8) {
    v1.v = _mm256_loadu_ps(&iArr[i]);             /* Load from in */
    v2.v = _mm256_loadu_ps(&oArr1[i]);            /* Load from out1 */
    vm.v = _mm256_cmp_ps(v1.v, v2.v, _CMP_GE_OQ); /* find new max */
    v2.v = _mm256_blendv_ps(v2.v,v1.v, vm.v);     /* Select out1 */
    _mm256_storeu_ps(&oArr1[i], v2.v);            /* Save out1 */
    v3.v = _mm256_loadu_ps(&oArr2[i]);            /* Load from out2 */
    v3.v = _mm256_blendv_ps(v3.v,vs.v, vm.v);     /* Select out2*/
    _mm256_storeu_ps(&oArr2[i], v3.v);            /* Save out2 */
    /* Select complex values using mask vm */
    k = 2*i;
    v1.v = _mm256_set_ps(iCArr[14+k],iCArr[12+k],iCArr[10+k],iCArr[8+k],
                         iCArr[6+k],iCArr[4+k],iCArr[2+k],iCArr[k]);
    v2.v = _mm256_set_ps(oCArr[14+k],oCArr[12+k],oCArr[10+k],oCArr[8+k],
                         oCArr[6+k],oCArr[4+k],oCArr[2+k],oCArr[k]);
    v3.v = _mm256_blendv_ps(v2.v,v1.v, vm.v);     /* reals */
    /* Grumble */
    for (j=0; j<8; j++) oCArr[k+j*2] = v3.v[j];
    /* Now imaginary */
    v1.v = _mm256_set_ps(iCArr[15+k],iCArr[13+k],iCArr[11+k],iCArr[9+k],
                         iCArr[7+k],iCArr[5+k],iCArr[3+k],iCArr[1+k]);
    v2.v = _mm256_set_ps(oCArr[15+k],oCArr[13+k],oCArr[11+k],oCArr[9+k],
                         oCArr[7+k],oCArr[5+k],oCArr[3+k],oCArr[1+k]);
    v3.v = _mm256_blendv_ps(v2.v,v1.v, vm.v);     /* imaginaries */
    /* Grumble */
    for (j=0; j<8; j++) oCArr[1+k+j*2] = v3.v[j];
  } /* End outer loop */
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  for (i=ilast; i<in->arraySize; i++) {
    if (iArr[i]>oArr1[i]) {
      oArr1[i] = iArr[i];
      oArr2[i] = value;
      oCArr[2*i]   = iCArr[2*i];
      oCArr[1+2*i] = iCArr[1+2*i];
    }
  } /* end loop over array */
} /* end ObitCArrayMaxSetValues */

/**
 *  Complex multiply two pairs of FArrays and complex accumulate
 *  complex accumulate - threading managed elsewhere.
 * Accum += Fin1 * Fin2
 *  ObitCArrayMakeCAFuncArgs to initialize threading with (elem)
 * \param Fin1_r   First FArray real part         (FA_1)
 * \param Fin1_i   First FArray imaginary part    (FA_2)
 * \param Fin2_r   Second FArray real part        (FA_3)
 * \param Fin2_i   Second FArray imaginary part   (FA_4)
 * \param Accum    Complex Accumulator            (out)
 * \param nThreads   Number of thread args allocated
 * \param threadArgs Thread arguments
 */
void ObitCArrayMulAccumTh (ObitFArray* Fin1_r, ObitFArray* Fin1_i, 
			   ObitFArray* Fin2_r, ObitFArray* Fin2_i, 
			   ObitCArray* Accum,
			   olong nThreads, gpointer inArgs) {
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread;
  CAFuncArg **threadArgs = (CAFuncArg**)inArgs;

  /* Divide up work */
  nElem = Fin1_r->arraySize;
  /* At least 50,000 per thread */
  nTh = MAX (1, MIN((olong)(0.5+nElem/50000.),nThreads));
  nElemPerThread = nElem/nTh;
  if (nElem<100000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem;
    threadArgs[i]->last    = hiElem;
    /* Reset pointers */
    if (threadArgs[i]->FA_1) ObitFArrayUnref(threadArgs[i]->FA_1);
    threadArgs[i]->FA_1 = ObitFArrayRef(Fin1_r);
    if (threadArgs[i]->FA_2) ObitFArrayUnref(threadArgs[i]->FA_2);
    threadArgs[i]->FA_2 = ObitFArrayRef(Fin1_i);
    if (threadArgs[i]->FA_3) ObitFArrayUnref(threadArgs[i]->FA_3);
    threadArgs[i]->FA_3 = ObitFArrayRef(Fin2_r);
    if (threadArgs[i]->FA_4) ObitFArrayUnref(threadArgs[i]->FA_4);
    threadArgs[i]->FA_4 = ObitFArrayRef(Fin2_i);
    if (threadArgs[i]->out) ObitCArrayUnref(threadArgs[i]->out);
    threadArgs[i]->out = ObitCArrayRef(Accum);
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  ObitThreadIterator (threadArgs[0]->thread, nTh, 
		      (ObitThreadFunc)ThreadCAMulAccum,
		      (gpointer**)threadArgs);

} /* end  ObitCArrayMulAccumTh */

/**
 *  Matrix inner multiply
 *  out = in1 * in2
 * in1, in2 and out can be square matrices of the same dimension or 
 * in1 or in2 and be a row or column vector
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitCArrayMatrixMult (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
{
  olong ir, ic, ii, n, nrow1, nrow2, ncol1, ncol2, nrowo, ncolo;
  olong indx1, indx2, indxo;
  ofloat sumr, sumi;
  
  /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));
  ncol1 = in1->naxis[0];
  if (in1->ndim>1) nrow1 = in1->naxis[1];
  else nrow1 = 1;
  ncol2 = in2->naxis[0];
  if (in2->ndim>1) nrow2 = in2->naxis[1];
  else nrow2 = 1;
  ncolo = out->naxis[0];
  if (out->ndim>1) nrowo = out->naxis[1];
  else nrowo = 1;
  g_assert (nrow1==ncol2);
  g_assert (ncolo==MIN(nrow1, nrow2));
  g_assert (nrowo==MIN(ncol1, ncol2));

  n = MAX (ncol1, nrow2);
  n = MAX (n, ncol2);
  n = MAX (n, nrow1);
  /* Both matrices */
  if ((ncol1>1) && (nrow1>1) && (ncol2>1) && (nrow2>1)) {
    for (ir=0; ir<nrowo; ir++) {
      for (ic=0; ic<ncolo; ic++) {
	sumr = sumi = 0.0;
	indxo = 2*(ic + ir*ncolo);
	for (ii=0; ii<n; ii++) {
	  indx1 = 2*(ii*ncol1 + ic);
	  indx2 = 2*(ii + ir*ncol2);
	  sumr += in1->array[indx1]*in2->array[indx2] - 
	    in1->array[indx1+1]*in2->array[indx2+1];
	  sumi += in1->array[indx1]*in2->array[indx2+1] + 
	    in1->array[indx1+1]*in2->array[indx2];
	}
	out->array[indxo]   = sumr;
	out->array[indxo+1] = sumi;
      }
    } /* end row loop */
    /* First a column vector */
  } else if (nrow1==1) {
    for (ir=0; ir<nrowo; ir++) {
      for (ic=0; ic<ncolo; ic++) {
	sumr = sumi = 0.0;
	indxo = 2*(ic + ir*ncolo);
	for (ii=0; ii<n; ii++) {
	  indx1 = 2*ii;
	  indx2 = 2*(ii + ir*ncol2);
	  sumr += in1->array[indx1]*in2->array[indx2] - 
	    in1->array[indx1+1]*in2->array[indx2+1];
	  sumi += in1->array[indx1]*in2->array[indx2+1] + 
	    in1->array[indx1+1]*in2->array[indx2];
	}
	out->array[indxo]   = sumr;
	out->array[indxo+1] = sumi;
      }
    } /* end row loop */
  } else if (nrow2==1) {  /* Second a column vector */
    for (ir=0; ir<nrowo; ir++) {
      for (ic=0; ic<ncolo; ic++) {
	sumr = sumi = 0.0;
	indxo = 2*(ic + ir*ncolo);
	for (ii=0; ii<n; ii++) {
	  indx1 = 2*(ii*ncol1 + ir);
	  indx2 = 2*ii;
	  sumr += in1->array[indx1]*in2->array[indx2] - 
	    in1->array[indx1+1]*in2->array[indx2+1];
	  sumi += in1->array[indx1]*in2->array[indx2+1] + 
	    in1->array[indx1+1]*in2->array[indx2];
	}
	out->array[indxo]   = sumr;
	out->array[indxo+1] = sumi;
      }
    } /* end row loop */
     /* First a row vector */
  } else if (ncol1==1) {
    for (ir=0; ir<nrowo; ir++) {
      for (ic=0; ic<ncolo; ic++) {
	sumr = sumi = 0.0;
	indxo = 2*(ic + ir*ncolo);
	for (ii=0; ii<n; ii++) {
	  indx1 = 2*ii;
	  indx2 = 2*(ii + ic*ncol2);
	  sumr += in1->array[indx1]*in2->array[indx2] - 
	    in1->array[indx1+1]*in2->array[indx2+1];
	  sumi += in1->array[indx1]*in2->array[indx2+1] + 
	    in1->array[indx1+1]*in2->array[indx2];
	}
	out->array[indxo]   = sumr;
	out->array[indxo+1] = sumi;
      }
    } /* end row loop */
  } else if (ncol2==1) {  /* Second a row vector */
    for (ir=0; ir<nrowo; ir++) {
      for (ic=0; ic<ncolo; ic++) {
	sumr = sumi = 0.0;
	indxo = 2*(ic + ir*ncolo);
	for (ii=0; ii<n; ii++) {
	  indx1 = 2*(ii*ncol1 + ic);
	  indx2 = 2*ii;
	  sumr += in1->array[indx1]*in2->array[indx2] - 
	    in1->array[indx1+1]*in2->array[indx2+1];
	  sumi += in1->array[indx1]*in2->array[indx2+1] + 
	    in1->array[indx1+1]*in2->array[indx2];
	}
	out->array[indxo]   = sumr;
	out->array[indxo+1] = sumi;
      }
    } /* end row loop */
 } /* End by type */
  
} /* end ObitCArrayMatrixMult */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitCArrayClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitCArrayClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitCArrayClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitCArrayClassInfoDefFn (gpointer inClass)
{
  ObitCArrayClassInfo *theClass = (ObitCArrayClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitCArrayClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitCArrayClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitCArrayGetClass;
  theClass->newObit       = (newObitFP)newObitCArray;
  theClass->ObitCopy      = (ObitCopyFP)ObitCArrayCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitCArrayClear;
  theClass->ObitInit      = (ObitInitFP)ObitCArrayInit;
  theClass->ObitCArrayCreate = (ObitCArrayCreateFP)ObitCArrayCreate;
  theClass->ObitCArrayIsCompatable = 
    (ObitCArrayIsCompatableFP)ObitCArrayIsCompatable;
  theClass->ObitCArrayIndex  = (ObitCArrayIndexFP)ObitCArrayIndex;
  theClass->ObitCArrayMaxAbs = (ObitCArrayMaxAbsFP)ObitCArrayMaxAbs;
  theClass->ObitCArrayMin    = (ObitCArrayMinFP)ObitCArrayMin;
  theClass->ObitCArrayNeg    = (ObitCArrayNegFP)ObitCArrayNeg;
  theClass->ObitCArrayConjg  = (ObitCArrayConjgFP)ObitCArrayConjg;
  theClass->ObitCArrayFill   = (ObitCArrayFillFP)ObitCArrayFill;
  theClass->ObitCArraySAdd   = (ObitCArraySAddFP)ObitCArraySAdd;
  theClass->ObitCArraySMul   = (ObitCArraySMulFP)ObitCArraySMul;
  theClass->ObitCArrayCSAdd  = (ObitCArrayCSAddFP)ObitCArrayCSAdd;
  theClass->ObitCArrayCSMul  = (ObitCArrayCSMulFP)ObitCArrayCSMul;
  theClass->ObitCArrayAdd    = (ObitCArrayAddFP)ObitCArrayAdd;
  theClass->ObitCArraySub    = (ObitCArraySubFP)ObitCArraySub;
  theClass->ObitCArrayMul    = (ObitCArrayMulFP)ObitCArrayMul;
  theClass->ObitCArrayDiv    = (ObitCArrayDivFP)ObitCArrayDiv;
  theClass->ObitCArrayMakeF  = (ObitCArrayMakeFFP)ObitCArrayMakeF;
  theClass->ObitCArrayMakeC  = (ObitCArrayMakeCFP)ObitCArrayMakeC;
  theClass->ObitCArrayIsFCompatable = 
    (ObitCArrayIsFCompatableFP)ObitCArrayIsFCompatable;
  theClass->ObitCArrayComplex= (ObitCArrayComplexFP)ObitCArrayComplex;
  theClass->ObitCArrayFMul   = (ObitCArrayFMulFP)ObitCArrayFMul;
  theClass->ObitCArrayFDiv   = (ObitCArrayFDivFP)ObitCArrayFDiv;
  theClass->ObitCArrayFAdd   = (ObitCArrayFAddFP)ObitCArrayFAdd;
  theClass->ObitCArrayFRot   = (ObitCArrayFRotFP)ObitCArrayFRot;
  theClass->ObitCArrayReal   = (ObitCArrayRealFP)ObitCArrayReal;
  theClass->ObitCArrayImag   = (ObitCArrayImagFP)ObitCArrayImag;
  theClass->ObitCArrayAmp    = (ObitCArrayAmpFP)ObitCArrayAmp;
  theClass->ObitCArrayAmp2   = (ObitCArrayAmp2FP)ObitCArrayAmp2;
  theClass->ObitCArrayPhase  = (ObitCArrayPhaseFP)ObitCArrayPhase;
  theClass->ObitCArray2DCenter = (ObitCArray2DCenterFP)ObitCArray2DCenter;
  theClass->ObitCArrayAddConjg = (ObitCArrayAddConjgFP)ObitCArrayAddConjg;
  theClass->ObitCArrayMatrixMult  = (ObitCArrayMatrixMultFP)ObitCArrayMatrixMult ;
} /* end ObitCArrayClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitCArrayInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitCArray *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread       = newObitThread();
  in->info         = newObitInfoList(); 
  in->array        = NULL;
  in->arraySize    = 0;
  in->ndim         = 0;
  in->naxis        = NULL;

} /* end ObitCArrayInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitCArray* cast to an Obit*.
 */
void ObitCArrayClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitCArray *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  if (in->array)  in->array = ObitMemFree(in->array);
  if (in->naxis)  in->naxis = ObitMemFree(in->naxis);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitCArrayClear */

/**
 * Make arguments for a Threaded ThreadCAFunc?
 * For use outside the ObitCArray class
 * \param thread     ObitThread object to be used
 * \param in         CA to be operated on
 * \param in2        2nd CA to be operated on
 * \param out        output CA
 * \param FA_1       First FArray
 * \param FA_2       Second FArray
 * \param FA_3       Third FArray
 * \param FA_4       Fourth FArray
 * \param larg1      Length of function dependent arg1 in floats
 * \param larg2      Length of function dependent arg2 in floats
 * \param larg3      Length of function dependent arg3 in floats
 * \param larg4      Length of function dependent arg4 in floats
 * \param larg5      Length of function dependent arg5 in floats
 * \param larg6      Length of function dependent arg6 in floats
 * \param larg7      Length of function dependent arg7 in floats
 * \param ThreadArgs[out] Created array of CAFuncArg, 
 *                   delete with ObitCArrayKillCAFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
olong ObitCArrayMakeCAFuncArgs (ObitThread *thread, ObitCArray *in,
				ObitCArray *in2, ObitCArray *out, 
				ObitFArray *FA_1, ObitFArray *FA_2, 
				ObitFArray *FA_3, ObitFArray *FA_4, 
				olong larg1, olong larg2, olong larg3, 
				olong larg4, olong larg5,
				olong larg6, olong larg7, 
				gpointer *outArgs)

{
  olong nThreads;
  CAFuncArg **ThreadArgs=NULL;
  /* Use class internal version */
  nThreads = MakeCAFuncArgs(thread, in, in2, out, FA_1, FA_2, FA_3, FA_4, 
			    larg1,larg2,larg3,larg4,larg5,larg6,larg7,
			    &ThreadArgs);
  
  *outArgs = (gpointer)ThreadArgs;
  return nThreads;
} /*  end ObitCArrayMakeCAFuncArgs */

/**
 * Make arguments for a Threaded ThreadCAFunc?
 * For use outside the ObitCArray class
 * \param thread     ObitThread object to be used
 * \param in         CA to be operated on
 * \param in2        2nd CA to be operated on
 * \param out        output CA
 * \param FA_1       First FArray
 * \param FA_2       Second FArray
 * \param larg1      Length of function dependent arg1 in floats
 * \param larg2      Length of function dependent arg2 in floats
 * \param larg3      Length of function dependent arg3 in floats
 * \param larg4      Length of function dependent arg4 in floats
 * \param larg5      Length of function dependent arg5 in floats
 * \param larg6      Length of function dependent arg6 in floats
 * \param larg7      Length of function dependent arg7 in floats
 * \param ThreadArgs[out] Created array of CAFuncArg, 
 *                   delete with ObitCArrayKillCAFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
olong ObitCArrayMakeCAFuncArgsX (ObitThread *thread, ObitCArray *in,
				ObitCArray *in2, ObitCArray *out, 
				ObitFArray *FA_1, ObitFArray *FA_2, 
				olong larg1, olong larg2, olong larg3, 
				olong larg4, olong larg5,
				olong larg6, olong larg7, 
				gpointer *outArgs)

{
  olong i, j, nThreads;
  ollong size=1;
  CAFuncArg **ThreadArgs=NULL;
  
  /* Get size from one of the inputs */
  if (in)        size = in->arraySize;
  else if (FA_1) size = FA_1->arraySize;
  else if (out)  size = out->arraySize;
   
  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  ThreadArgs = g_malloc0(nThreads*sizeof(CAFuncArg*));
  for (i=0; i<nThreads; i++) 
    (ThreadArgs)[i] = g_malloc0(sizeof(CAFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (ThreadArgs)[i]->thread= ObitThreadRef(thread);
    (ThreadArgs)[i]->in    = ObitCArrayRef(in);
    if (in2) (ThreadArgs)[i]->in2   = ObitCArrayRef(in2);
    else (ThreadArgs)[i]->in2   = NULL;
    if (out) (ThreadArgs)[i]->out   = ObitCArrayRef(out);
    else (ThreadArgs)[i]->out = NULL;
    if (FA_1) (ThreadArgs)[i]->FA_1   = ObitFArrayRef(FA_1);
    else (ThreadArgs)[i]->FA_1 = NULL;
    if (FA_1) (ThreadArgs)[i]->FA_2   = ObitFArrayRef(FA_2);
    else (ThreadArgs)[i]->FA_2 = NULL;
    (ThreadArgs)[i]->first = 1;
    (ThreadArgs)[i]->last  = size;
    (ThreadArgs)[i]->value = 0.0;
    for (j=0; j<MAXCARRAYDIM; j++) (ThreadArgs)[i]->pos[j] = 0;
    if (larg1>0) (ThreadArgs)[i]->arg1 = g_malloc0(larg1*sizeof(float));
    else (ThreadArgs)[i]->arg1 = NULL;
    if (larg2>0) (ThreadArgs)[i]->arg2 = g_malloc0(larg2*sizeof(float));
    else (ThreadArgs)[i]->arg2 = NULL;
    if (larg3>0) (ThreadArgs)[i]->arg3 = g_malloc0(larg3*sizeof(float));
    else (ThreadArgs)[i]->arg3 = NULL;
    if (larg4>0) (ThreadArgs)[i]->arg4 = g_malloc0(larg4*sizeof(float));
    else (ThreadArgs)[i]->arg4 = NULL;
    if (larg5>0) (ThreadArgs)[i]->arg5 = g_malloc0(larg5*sizeof(float));
    else (ThreadArgs)[i]->arg5 = NULL;
    if (larg5>0) (ThreadArgs)[i]->arg6 = g_malloc0(larg6*sizeof(float));
    else (ThreadArgs)[i]->arg6 = NULL;
    if (larg5>0) (ThreadArgs)[i]->arg7 = g_malloc0(larg7*sizeof(float));
    else (ThreadArgs)[i]->arg7 = NULL;
    (ThreadArgs)[i]->ithread  = i;
  }
  *outArgs = (gpointer)ThreadArgs;
  return nThreads;
} /*  end ObitCArrayMakeCAFuncArgsX */

/**
 * Delete arguments for ThreadCAFunc
 * \param nargs      number of elements in ThreadArgs.
 * \param inArgs     pointer to array of CAFuncArg
 */
void ObitCArrayKillCAFuncArgs (olong nargs, gpointer inArgs)
{
  CAFuncArg **ThreadArgs = (CAFuncArg**)inArgs;

  if (ThreadArgs==NULL) return;

  /* Use class function */
  KillCAFuncArgs (nargs, ThreadArgs);
} /*  end ObitCArrayKillCAFuncArgs */

/**
 * Add portions of two CArrays, out = in + in2
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to CAFuncArg argument with elements:
 * \li in       ObitCArray to work on
 * \li in2      2nd ObitCArray to work on
 * \li out      Output ObitCArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadCAAdd (gpointer arg)
{
  /* Get arguments from structure */
  CAFuncArg *largs = (CAFuncArg*)arg;
  ObitCArray *in1       = largs->in;
  ObitCArray *in2       = largs->in2;
  ObitCArray *out       = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  /* local */
  olong   i, ilast;
  ofloat  fblank = ObitMagicF();
#if HAVE_AVX512==1
  CV16SF v1, v2, vb;
  MASK16 msk1, msk2;
#elif HAVE_AVX==1
  CV8SF v1, v2, vb, vm1, vm2;
#endif

  if (hiElem<loElem) goto finish;

  /* Loop over array */
#if HAVE_AVX512==1  /* Vector */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  /* Do blocks of 16 as vector */
  for (i=loElem; i<hiElem-16; i+=16) {
    v1.v  = _mm512_loadu_ps(&in1->array[i]);
    msk1  = _mm512_cmp_ps_mask(v1.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v2.v  = _mm512_loadu_ps(&in2->array[i]);
    msk2  = _mm512_cmp_ps_mask(v2.v, vb.v, _CMP_EQ_OQ);  /* find blanks */
    v2.v  = _mm512_add_ps(v1.v, v2.v);             /* do operation */
    v2.v  = _mm512_mask_blend_ps(msk1,v2.v,vb.v);    /* replace blanks 1 */
    v2.v  = _mm512_mask_blend_ps(msk2,v2.v,vb.v);    /* replace blanks 2 */
    _mm512_storeu_ps(&out->array[i], v2.v);        /* Save */
  }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* Vector */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  /* Do blocks of 8 as vector */
  for (i=loElem; i<hiElem-8; i+=8) {
    v1.v  = _mm256_loadu_ps(&in1->array[i]);
    vm1.v = _mm256_cmp_ps(v1.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v2.v  = _mm256_loadu_ps(&in2->array[i]);
    vm2.v = _mm256_cmp_ps(v2.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v2.v  = _mm256_add_ps(v1.v, v2.v);             /* do operation */
    v2.v  = _mm256_blendv_ps(v2.v,vb.v, vm1.v);    /* replace blanks 1 */
    v2.v  = _mm256_blendv_ps(v2.v,vb.v, vm2.v);    /* replace blanks 2 */
    _mm256_storeu_ps(&out->array[i], v2.v);        /* Save */
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = loElem;  /* Do all */
#endif
  for (i=ilast; i<hiElem; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = in1->array[i] + in2->array[i];
    else out->array[i] = fblank;
  }

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadCAAdd */

/**
 * Subtract portions of two CArrays, out = in1 - in2
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to CAFuncArg argument with elements:
 * \li in       ObitCArray to work on
 * \li in2      2nd ObitCArray to work on
 * \li out      Output ObitCArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadCASub (gpointer arg)
{
  /* Get arguments from structure */
  CAFuncArg *largs = (CAFuncArg*)arg;
  ObitCArray *in1       = largs->in;
  ObitCArray *in2       = largs->in2;
  ObitCArray *out       = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  /* local */
  olong   i, ilast;
  ofloat  fblank = ObitMagicF();
#if HAVE_AVX512==1
  CV16SF v1, v2, vb;
  MASK16 msk1, msk2;
#elif HAVE_AVX==1
  CV8SF v1, v2, vb, vm1, vm2;
#endif

  if (hiElem<loElem) goto finish;

  /* Loop over array */
#if HAVE_AVX512==1  /* Vector */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  /* Do blocks of 16 as vector */
  for (i=loElem; i<hiElem-16; i+=16) {
    v1.v  = _mm512_loadu_ps(&in1->array[i]);
    msk1  = _mm512_cmp_ps_mask(v1.v, vb.v, _CMP_EQ_OQ);   /* find blanks */
    v2.v  = _mm512_loadu_ps(&in2->array[i]);
    msk2  = _mm512_cmp_ps_mask(v2.v, vb.v, _CMP_EQ_OQ);   /* find blanks */
    v2.v  = _mm512_sub_ps(v1.v, v2.v);               /* do operation */
    v2.v  = _mm512_mask_blend_ps(msk1,v2.v,vb.v);    /* replace blanks 1 */
    v2.v  = _mm512_mask_blend_ps(msk2,v2.v,vb.v);    /* replace blanks 2 */
    _mm512_storeu_ps(&out->array[i], v2.v);
  }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* Vector */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  /* Do blocks of 8 as vector */
  for (i=loElem; i<hiElem-8; i+=8) {
    v1.v  = _mm256_loadu_ps(&in1->array[i]);
    vm1.v = _mm256_cmp_ps(v1.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v2.v  = _mm256_loadu_ps(&in2->array[i]);
    vm2.v = _mm256_cmp_ps(v2.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v2.v  = _mm256_sub_ps(v1.v, v2.v);             /* do operation */
    v2.v  = _mm256_blendv_ps(v2.v,vb.v, vm1.v);    /* replace blanks 1 */
    v2.v  = _mm256_blendv_ps(v2.v,vb.v, vm2.v);    /* replace blanks 2 */
    _mm256_storeu_ps(&out->array[i], v2.v);
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = loElem;  /* Do all */
#endif
  for (i=ilast; i<hiElem; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = in1->array[i] - in2->array[i];
    else out->array[i] = fblank;
  }

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadCASub */

/**
 * Multiply portions of two CArrays, out = in * in2
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to CAFuncArg argument with elements:
 * \li in       ObitCArray to work on
 * \li in2      2nd ObitCArray to work on
 * \li out      Output ObitCArray
 * \li first    First complex element (1-rel) number
 * \li last     Highest complex element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadCAMul (gpointer arg)
{
  /* Get arguments from structure */
  CAFuncArg *largs = (CAFuncArg*)arg;
  ObitCArray *in1       = largs->in;
  ObitCArray *in2       = largs->in2;
  ObitCArray *out       = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  olong i, iLast;
  ofloat tr1, ti1, tr2, ti2;
  ofloat  fblank = ObitMagicF();

#if HAVE_AVX512==1
  CV16SF v1r, v2r, v1i, v2i, vb, vt1, vt2, vr, vi, vm;
  MASK16 msk1, msk2, selMask;
#elif HAVE_AVX==1
  CV8SF vm1, vm2;
  CV8SF v1r, v1i, vb, v2r, v2i, vt1, vt2, vr, vi, vm;
#endif
  if (hiElem<loElem) goto finish;

  /* AVX?  some extra work here due to lack of good gather/scatter */
#if HAVE_AVX512==1  /* Vector AVX 512 */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  /* Do blocks of 16 as vector  = 8 complex **/
  for (i=loElem; i<hiElem-8; i+=8) {
    v1r.v = _mm512_loadu_ps(&in1->array[2*i]);           /* in1 complex */
    msk1  = _mm512_cmp_ps_mask(v1r.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v2r.v = _mm512_loadu_ps(&in2->array[2*i]);           /* in2 complex */
    msk2  = _mm512_cmp_ps_mask(v2r.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    /* convert to pairs of real/imag */
    v1i.v = _mm512_movehdup_ps (v1r.v);  /* Duplicate 1 Odds */
    v2i.v = _mm512_movehdup_ps (v2r.v);  /* Duplicate 2 Odds */
    v1r.v = _mm512_moveldup_ps (v1r.v);  /* Duplicate 1 Evens */
    v2r.v = _mm512_moveldup_ps (v2r.v);  /* Duplicate 2 Evens */
    /* Multiply */
    vt1.v = _mm512_mul_ps (v1r.v, v2r.v);  /* tr1*tr2 in pairs of entries */
    vt2.v = _mm512_mul_ps (v1i.v, v2i.v);  /* ti1*ti2 */
    vr.v  = _mm512_sub_ps (vt1.v, vt2.v);  /* Real part tr1*tr2 - ti1*ti2*/
    vt1.v = _mm512_mul_ps (v1i.v, v2r.v);  /* ti1*tr2 */
    vt2.v = _mm512_mul_ps (v1r.v, v2i.v);  /* tr1*ti2 */
    vi.v  = _mm512_add_ps (vt1.v, vt2.v);  /* Imaginary part ti1*tr2 +tr1*tri2 */
    selMask = _mm512_int2mask(0xaaaa);     /* alternating mask */
    vm.v  = _mm512_mask_blend_ps (selMask, vr.v, vi.v);  /* merge back to complex */
    vm.v  = _mm512_mask_blend_ps(msk1,vm.v,vb.v);    /* replace blanks 1 */
    vm.v  = _mm512_mask_blend_ps(msk2,vm.v,vb.v);    /* replace blanks 2 */
    _mm512_storeu_ps(&out->array[2*i], vm.v);        /* Save */
  }
  iLast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* Vector AVX */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  /* Do blocks of 8 as vector = 4 complex */
  for (i=loElem; i<hiElem-4; i+=4) {
    v1r.v = _mm256_loadu_ps(&in1->array[2*i]);       /* in1 complex */
    vm1.v = _mm256_cmp_ps(v1r.v, vb.v, _CMP_EQ_OQ);  /* find blanks */
    v2r.v = _mm256_loadu_ps(&in2->array[2*i]);       /* in2 complex */
    vm2.v = _mm256_cmp_ps(v2r.v, vb.v, _CMP_EQ_OQ);  /* find blanks */
    /* convert to pairs of real/imag */
    v1i.v = _mm256_movehdup_ps (v1r.v);  /* Duplicate 1 Odds */
    v2i.v = _mm256_movehdup_ps (v2r.v);  /* Duplicate 2 Odds */
    v1r.v = _mm256_moveldup_ps (v1r.v);  /* Duplicate 1 Evens */
    v2r.v = _mm256_moveldup_ps (v2r.v);  /* Duplicate 2 Evens */
    /* Multiply */
    vt1.v = _mm256_mul_ps (v1r.v, v2r.v);  /* tr1*tr2 in pairs of entries */
    vt2.v = _mm256_mul_ps (v1i.v, v2i.v);  /* ti1*ti2 */
    vr.v  = _mm256_sub_ps (vt1.v, vt2.v);  /* Real part tr1*tr2 - ti1*ti2 */
    vt1.v = _mm256_mul_ps (v1i.v, v2r.v);  /* ti1*tr2 */
    vt2.v = _mm256_mul_ps (v1r.v, v2i.v);  /* tr1*ti2 */
    vi.v  = _mm256_add_ps (vt1.v, vt2.v);  /* Imaginary part ti1*tr2 +tr1*tri2 */
    /* Merge (blend) AVX no Mask*/
    vm.v  = _mm256_blend_ps (vr.v, vi.v, 0xaa);  /* merge back to complex */
    vm.v  = _mm256_blendv_ps(vm.v,vb.v,vm1.v);   /* replace blanks 1 */
    vm.v  = _mm256_blendv_ps(vm.v,vb.v,vm2.v);   /* replace blanks 2 */
    _mm256_storeu_ps(&out->array[2*i], vm.v);    /* Save */
  }
  iLast = i;  /* How far did I get? */
#else /* Scalar */
  iLast = loElem;  /* Do all */
#endif
  /* Finish up or scalar */
  for (i=iLast; i<hiElem; i++) {
    tr1 = in1->array[2*i];
    ti1 = in1->array[2*i+1];
    tr2 = in2->array[2*i];
    ti2 = in2->array[2*i+1];
    if ((tr1!=fblank) && (tr2!=fblank) && (ti1!=fblank) && (ti2!=fblank)) {
      out->array[2*i]   = tr1*tr2 - ti1*ti2;
      out->array[2*i+1] = ti1*tr2 + tr1*ti2;
    } else out->array[2*i] = fblank;
  }

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadCAMul */

/**
 * Divide portions of two CArrays, out = in / in2
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to CAFuncArg argument with elements:
 * \li in       ObitCArray to work on
 * \li in2      2nd ObitCArray to work on
 * \li out      Output ObitCArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadCADiv (gpointer arg)
{
  /* Get arguments from structure */
  CAFuncArg *largs = (CAFuncArg*)arg;
  ObitCArray *in1       = largs->in;
  ObitCArray *in2       = largs->in2;
  ObitCArray *out       = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  olong i;
  ofloat tr1, ti1, tr2, ti2, denom;
  ofloat  fblank = ObitMagicF();
  if (hiElem<loElem) goto finish;


  for (i=loElem; i<hiElem; i++) {
    tr1 = in1->array[2*i];
    ti1 = in1->array[2*i+1];
    tr2 = in2->array[2*i];
    ti2 = in2->array[2*i+1];
    if ((tr1!=fblank) && (tr2!=fblank) && (ti1!=fblank) && (ti2!=fblank)) {
      denom = (tr2*tr2 + ti2*ti2);
      if (denom!=0.0) {
	out->array[2*i]   = (tr1*tr2 + ti1*ti2)/denom;
	out->array[2*i+1] = (ti1*tr2 - tr1*ti2)/denom;
      } else out->array[2*i] = fblank;
    } else out->array[2*i] = fblank;
  } /* end loop */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadCADiv */

/**
 * scalar fill portions a CArray
 * Callable as thread
 * \param arg Pointer to CAFuncArg argument with elements:
 * \li in       ObitCArray to work on
 * \li arg1     Complex scalar as real, imag
 * \li first    First complex element (1-rel) number
 * \li last     Highest complex element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadCAFill (gpointer arg)
{
  /* Get arguments from structure */
  CAFuncArg *largs = (CAFuncArg*)arg;
  ObitCArray *in        = largs->in;
  ofloat     *cmpx      = largs->arg1;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  olong i;
  if (hiElem<loElem) goto finish;

  for (i=loElem; i<hiElem; i++) {
    in->array[2*i]   = cmpx[0];
    in->array[2*i+1] = cmpx[1];

  }
  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadCAFill */

/**
 * Write amplitudes of a complex array as an FArray
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to CAFuncArg argument with elements:
 * \li in       ObitCArray to work on
 * \li FA_1     Output ObitFArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadCAAmp (gpointer arg)
{
  /* Get arguments from structure */
  CAFuncArg *largs = (CAFuncArg*)arg;
  ObitCArray *in1       = largs->in;
  ObitFArray *FA_1      = largs->FA_1;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  olong i, ilast;
  ofloat  tr1, ti1, fblank = ObitMagicF();
#if HAVE_AVX512==1
  CV16SF vr, vi, vb, vamp;
  IV16SF vindxr, vindxi;
  MASK16 msk1, msk2;
#endif

if (hiElem<loElem) goto finish;

#if HAVE_AVX512==1  /* Vector length 16 */
  vb.v     = _mm512_set1_ps(fblank);  /* vector of blanks */
  vindxr.v = _mm512_set_epi32 (30,28,26,24,22,20,18,16,14,12,10,8,6,4,2,0); /* Reals*/
  vindxi.v = _mm512_set_epi32 (31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1); /* Imags */
  /* Do blocks of 16 as vector */
  for (i=loElem; i<hiElem-16; i+=16) {
    vr.v  = _mm512_i32gather_ps (vindxr.v, &in1->array[2*i], 4);/* Load Reals */
    msk1  = _mm512_cmp_ps_mask(vr.v, vb.v, _CMP_EQ_OQ);         /* find blanks */
    vi.v  = _mm512_i32gather_ps (vindxi.v, &in1->array[2*i], 4);/* Load Imaginaries */
    msk2  = _mm512_cmp_ps_mask(vi.v, vb.v, _CMP_EQ_OQ);         /* find blanks */
    vr.v  = _mm512_mul_ps (vr.v, vr.v);                         /* Real squared */
    vi.v  = _mm512_mul_ps (vi.v, vi.v);                         /* Imag squared */
    vamp.v  = _mm512_add_ps(vr.v, vi.v);                        /* R^2 + i^2 */
    vamp.v  = _mm512_sqrt_ps(vamp.v);                            /* sqrt(R^2 + i^2) */
    vamp.v  = _mm512_mask_blend_ps(msk1,vamp.v,vb.v);           /* replace blanks real */
    vamp.v  = _mm512_mask_blend_ps(msk2,vamp.v,vb.v);           /* replace blanks imag */
    _mm512_storeu_ps(&FA_1->array[i], vamp.v);                   /* Save */
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = loElem;  /* Do all */
#endif
  for (i=ilast; i<hiElem; i++) {
    tr1 = in1->array[2*i];
    ti1 = in1->array[2*i+1];
    if ((tr1!=fblank) && (ti1!=fblank)) {
      FA_1->array[i] = sqrt (tr1*tr1 + ti1*ti1);
    } else FA_1->array[i] = fblank;
  } /* end loop */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadCAAmp */

/**
 * Write phases of a complex array as an FArray
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to CAFuncArg argument with elements:
 * \li in       ObitCArray to work on
 * \li FA_1     Output ObitFArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadCAPhase (gpointer arg)
{
  /* Get arguments from structure */
  CAFuncArg *largs = (CAFuncArg*)arg;
  ObitCArray *in1       = largs->in;
  ObitFArray *FA_1      = largs->FA_1;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  olong i;
  ofloat  tr1, ti1, fblank = ObitMagicF();
  if (hiElem<loElem) goto finish;

  for (i=loElem; i<hiElem; i++) {
    tr1 = in1->array[2*i];
    ti1 = in1->array[2*i+1];
    if ((tr1!=fblank) && (ti1!=fblank)) {
      FA_1->array[i] = atan2 (ti1, tr1+1.0e-20);
    } else FA_1->array[i] = fblank;
  } /* end loop */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadCAPhase */

/**
 * Form complex from two FArrays, multiply by complex scalar and 
 * complex accumulate. Magic value blanking not supported.
 * Callable as thread
 * \param arg Pointer to CAFuncArg argument with elements:
 * \li FA_1     Real part of input
 * \li FA_2     Imaginary part of input
 * \li arg1     Real/Imaginary parts of scalar
 * \li out      Output ObitCArray accumulator
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadCASMulAccum (gpointer arg)
{
  /* Get arguments from structure */
  CAFuncArg *largs = (CAFuncArg*)arg;
  ObitFArray *Fin_r     = largs->FA_1;
  ObitFArray *Fin_i     = largs->FA_2;
  ofloat     *cscalar   = largs->arg1;
  ObitCArray *Accum     = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  ofloat tr1, ti1, tr2, ti2;
  olong i, ilast;
  ofloat *rArr = Fin_r->array, *iArr = Fin_i->array, *oArr = Accum->array;
#if HAVE_AVX512==1
  CV16SF v1r, v1i, sr, si, tv1, tv2, tv3, tv4;
  V16SI vindex;
#elif HAVE_AVX==1
  olong j, k, l;
  CV8SF  v1r, v1i, sr, si, tv1, tv2, tv3, tv4;
  CV4SF  v4, vr4, vi4, vc4 ;
#endif
  
  if (hiElem<loElem) goto finish;
  tr2 = cscalar[0]; ti2 = cscalar[1];
  
#if HAVE_AVX512==1  /* AVX 512 Vector (16 float) */
  sr.v = _mm512_set1_ps(tr2); /* Scalar real */
  si.v = _mm512_set1_ps(ti2); /* Scalar imaginary */
  vindex = _mm512_set_epi32(30,28,26,24,22,20,18,16,14,12,10, 8, 6,4,2,0);
  for (i=loElem; i<hiElem-16; i+=16) {
    v1r.v = _mm512_loadu_ps(&rArr[i]); /* Input reals */
    v1i.v = _mm512_loadu_ps(&iArr[i]); /* Input imaginaries */
    tv1.v = _mm512_mul_ps (v1r.v, sr.v);  /* tr1*tr2 */
    tv2.v = _mm512_mul_ps (v1i.v, si.v);  /* ti1*ti2 */
    tv1.v = _mm512_sub_ps (tv1.v, tv2.v); /* now Real part */
    tv3.v = _mm512_i32gather_ps(vindex, (void const*)(&oArr[i*2]), 4); /* Gather reals */
    tv3.v = _mm512_add_ps (tv3.v, tv1.v);  /* Update reals */
    _mm512_i32scatter_ps((void*)(&oArr[i*2]), vindex, tv3.v, 4); /* Scatter reals */
    tv3.v = _mm512_mul_ps (v1i.v, sr.v);  /* ti1*tr2 */
    tv4.v = _mm512_mul_ps (v1r.v, si.v);  /* tr1*ti2 */
    tv2.v = _mm512_add_ps (tv3.v, tv4.v); /* now Imaginary part */
    tv3.v = _mm512_i32gather_ps(vindex, (void const*)(&oArr[1+(i*2)]), 4); /* Gather imag */
    tv3.v = _mm512_add_ps (tv3.v, tv2.v);  /* Update imag */
    _mm512_i32scatter_ps((void*)(&oArr[1+(i*2)]), vindex, tv3.v, 4); /* Scatter imag */
  } /* end outer loop */
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* AVX Vector (8 float) */ 
  sr.v = _mm256_broadcast_ss(&tr2); /* Scalar real */
  si.v = _mm256_broadcast_ss(&ti2); /* Scalar imaginary */
  for (i=loElem; i<hiElem-8; i+=8) {
    v1r.v = _mm256_loadu_ps(&rArr[i]); /* Input reals */
    v1i.v = _mm256_loadu_ps(&iArr[i]); /* Input imaginaries */
    tv1.v = _mm256_mul_ps (v1r.v, sr.v);  /* tr1*tr2 */
    tv2.v = _mm256_mul_ps (v1i.v, si.v);  /* ti1*ti2 */
    tv1.v = _mm256_sub_ps (tv1.v, tv2.v); /* now Real part */
    tv3.v = _mm256_mul_ps (v1i.v, sr.v);  /* ti1*tr2 */
    tv4.v = _mm256_mul_ps (v1r.v, si.v);  /* tr1*ti2 */
    tv2.v = _mm256_add_ps (tv3.v, tv4.v); /* now Imaginary part */
    /* use 4 float vectors to update output, they can be shuffled */
    k = 2*i;
    for (j=0; j<2; j++) {
      l = k + j*8;
      if (j==0) {
	vr4.v = _mm256_extractf128_ps(tv1.v, 0x0);  /* extract 4 floats from Real */
	vi4.v = _mm256_extractf128_ps(tv2.v, 0x0);  /* extract 4 floats from Imaginary */
      } else {	
   	vr4.v = _mm256_extractf128_ps(tv1.v, 0x1);  /* extract 4 floats from Real */
	vi4.v = _mm256_extractf128_ps(tv2.v, 0x1);  /* extract 4 floats from Imaginary */
      }
      v4.v  = _mm_loadu_ps(&oArr[l]);        /* Input next four output floats */
      vc4.v = _mm_unpacklo_ps(vr4.v, vi4.v); /* first 2 of 8 complex values */
      v4.v  = _mm_add_ps (v4.v, vc4.v);      /* update */
      _mm_storeu_ps (&oArr[l], v4.v);        /* rewrite */
      v4.v  = _mm_loadu_ps(&oArr[l+4]);      /* Next four output floats */
      vc4.v = _mm_unpackhi_ps(vr4.v, vi4.v); /* second 2 of 8 complex values */
      v4.v  = _mm_add_ps (v4.v, vc4.v);      /* update */
      _mm_storeu_ps (&oArr[l+4], v4.v);      /* rewrite */
    }
  } /* end outer loop */
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  /* Loop over whatever is left over */
  for (i=ilast; i<hiElem; i++) {
    tr1 = rArr[i];
    ti1 = iArr[i];
    oArr[i*2]   += tr1*tr2 - ti1*ti2;
    oArr[1+i*2] += ti1*tr2 + tr1*ti2;
  } /* end loop over array */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  } /* end ThreadCASMulAccum */

/**
 * Form complex from two FArrays, multiply by complex scalar and 
 * complex accumulate. Magic value blanking not supported.
 * Callable as thread
 * \param arg Pointer to CAFuncArg argument with elements:
 * \li FA_1     Real part of 1st input
 * \li FA_2     Imaginary part of 1st input
 * \li FA_3     Real part of 2nd input
 * \li FA_4     Imaginary part of 2nd input
 * \li out      Output ObitCArray accumulator
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadCAMulAccum (gpointer arg)
{
  /* Get arguments from structure */
  CAFuncArg *largs = (CAFuncArg*)arg;
  ObitFArray *Fin1_r    = largs->FA_1;
  ObitFArray *Fin1_i    = largs->FA_2;
  ObitFArray *Fin2_r    = largs->FA_3;
  ObitFArray *Fin2_i    = largs->FA_4;
  ObitCArray *Accum     = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  ofloat tr1, ti1, tr2, ti2;
  olong i, ilast;
  ofloat *rArr1 = Fin1_r->array, *iArr1 = Fin1_i->array, *oArr = Accum->array;
  ofloat *rArr2 = Fin2_r->array, *iArr2 = Fin2_i->array;
#if HAVE_AVX512==1
  CV16SF v1r, v1i, v2r, v2i, tv1, tv2, tv3, tv4;
  V16SI vindex;
#elif HAVE_AVX==1
  olong j, k, l;
  CV8SF  v1r, v1i, v2r, v2i, tv1, tv2, tv3, tv4;
  CV4SF  v4, vr4, vi4, vc4 ;
#endif
  
  if (hiElem<loElem) goto finish;
  
#if HAVE_AVX512==1  /* AVX 512 Vector (16 float) */
  vindex = _mm512_set_epi32(30,28,26,24,22,20,18,16,14,12,10, 8, 6,4,2,0);
  for (i=loElem; i<hiElem-16; i+=16) {
    v1r.v = _mm512_loadu_ps(&rArr1[i]);    /* Input 1 reals */
    v1i.v = _mm512_loadu_ps(&iArr1[i]);    /* Input 1 imaginaries */
    v2r.v = _mm512_loadu_ps(&rArr2[i]);    /* Input 2 reals */
    v2i.v = _mm512_loadu_ps(&iArr2[i]);    /* Input 2 imaginaries */
    tv1.v = _mm512_mul_ps (v1r.v, v2r.v);  /* tr1*tr2 */
    tv2.v = _mm512_mul_ps (v1i.v, v2i.v);  /* ti1*ti2 */
    tv1.v = _mm512_sub_ps (tv1.v, tv2.v);  /* now Real part */
    tv3.v = _mm512_i32gather_ps(vindex, (void const*)(&oArr[i*2]), 4); /* Gather reals */
    tv3.v = _mm512_add_ps (tv3.v, tv1.v);  /* Update reals */
    _mm512_i32scatter_ps((void*)(&oArr[i*2]), vindex, tv3.v, 4); /* Scatter reals */
    tv3.v = _mm512_mul_ps (v1i.v, v2r.v);  /* ti1*tr2 */
    tv4.v = _mm512_mul_ps (v1r.v, v2i.v);  /* tr1*ti2 */
    tv2.v = _mm512_add_ps (tv3.v, tv4.v);  /* now Imaginary part */
    tv3.v = _mm512_i32gather_ps(vindex, (void const*)(&oArr[1+(i*2)]), 4); /* Gather imag */
    tv3.v = _mm512_add_ps (tv3.v, tv2.v);  /* Update imag */
    _mm512_i32scatter_ps((void*)(&oArr[1+(i*2)]), vindex, tv3.v, 4); /* Scatter imag */
  } /* end outer loop */
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* AVX Vector (8 float) */ 
  for (i=loElem; i<hiElem-8; i+=8) {
    v1r.v = _mm256_loadu_ps(&rArr1[i]);    /* Input 1 reals */
    v1i.v = _mm256_loadu_ps(&iArr1[i]);    /* Input 1 imaginaries */
    v2r.v = _mm256_loadu_ps(&rArr2[i]);    /* Input 2 reals */
    v2i.v = _mm256_loadu_ps(&iArr2[i]);    /* Input 2 imaginaries */
    tv1.v = _mm256_mul_ps (v1r.v, v2r.v);  /* tr1*tr2 */
    tv2.v = _mm256_mul_ps (v1i.v, v2i.v);  /* ti1*ti2 */
    tv1.v = _mm256_sub_ps (tv1.v, tv2.v);  /* now Real part */
    tv3.v = _mm256_mul_ps (v1i.v, v2r.v);  /* ti1*tr2 */
    tv4.v = _mm256_mul_ps (v1r.v, v2i.v);  /* tr1*ti2 */
    tv2.v = _mm256_add_ps (tv3.v, tv4.v);  /* now Imaginary part */
    /* use 4 float vectors to update output, they can be shuffled */
    k = 2*i;
    for (j=0; j<2; j++) {
      l = k + j*8;
      if (j==0) {
	vr4.v = _mm256_extractf128_ps(tv1.v, 0x0);  /* extract 4 floats from Real */
	vi4.v = _mm256_extractf128_ps(tv2.v, 0x0);  /* extract 4 floats from Imaginary */
      } else {	
   	vr4.v = _mm256_extractf128_ps(tv1.v, 0x1);  /* extract 4 floats from Real */
	vi4.v = _mm256_extractf128_ps(tv2.v, 0x1);  /* extract 4 floats from Imaginary */
      }
      v4.v  = _mm_loadu_ps(&oArr[l]);        /* Input next four output floats */
      vc4.v = _mm_unpacklo_ps(vr4.v, vi4.v); /* first 2 of 8 complex values */
      v4.v  = _mm_add_ps (v4.v, vc4.v);      /* update */
      _mm_storeu_ps (&oArr[l], v4.v);        /* rewrite */
      v4.v  = _mm_loadu_ps(&oArr[l+4]);      /* Next four output floats */
      vc4.v = _mm_unpackhi_ps(vr4.v, vi4.v); /* first 2 of 8 complex values */
      v4.v  = _mm_add_ps (v4.v, vc4.v);      /* update */
      _mm_storeu_ps (&oArr[l+4], v4.v);      /* rewrite */
    }
  } /* end outer loop */
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  /* Loop over whatever is left over */
  for (i=ilast; i<hiElem; i++) {
    tr1 = rArr1[i];  ti1 = iArr1[i];
    tr2 = rArr2[i];  ti2 = iArr2[i];
    oArr[i*2]   += tr1*tr2 - ti1*ti2;
    oArr[1+i*2] += ti1*tr2 + tr1*ti2;
  } /* end loop over array */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  } /* end ThreadMulAccum */

/**
 * Make arguments for a Threaded ThreadCAFunc?
 * \param thread     ObitThread object to be used
 * \param in         CA to be operated on
 * \param in2        2nd CA to be operated on
 * \param out        output CA
 * \param FA_1       First FArray
 * \param FA_2       Second FArray
 * \param FA_3       Third FArray
 * \param FA_4       Fourth FArray
 * \param larg1      Length of function dependent arg1 in floats
 * \param larg2      Length of function dependent arg2 in floats
 * \param larg3      Length of function dependent arg3 in floats
 * \param larg4      Length of function dependent arg4 in floats
 * \param larg5      Length of function dependent arg5 in floats
 * \param larg6      Length of function dependent arg6 in floats
 * \param larg7      Length of function dependent arg7 in floats
 * \param ThreadArgs[out] Created array of CAFuncArg, 
 *                   delete with KillCAFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeCAFuncArgs (ObitThread *thread, ObitCArray *in,
			     ObitCArray *in2, ObitCArray *out, 
			     ObitFArray *FA_1, ObitFArray *FA_2, 
			     ObitFArray *FA_3, ObitFArray *FA_4, 
			     olong larg1, olong larg2, olong larg3, 
			     olong larg4, olong larg5,
			     olong larg6, olong larg7, 
			     CAFuncArg ***ThreadArgs)

{
  olong i, j, nThreads;
  ollong size=1;

  /* Get size from one of the inputs */
  if (in)        size = in->arraySize;
  else if (FA_1) size = FA_1->arraySize;
  else if (out)  size = out->arraySize;
   
  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(CAFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(CAFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->thread= ObitThreadRef(thread);
    (*ThreadArgs)[i]->in    = ObitCArrayRef(in);
    if (in2) (*ThreadArgs)[i]->in2   = ObitCArrayRef(in2);
    else (*ThreadArgs)[i]->in2   = NULL;
    if (out) (*ThreadArgs)[i]->out   = ObitCArrayRef(out);
    else (*ThreadArgs)[i]->out = NULL;
    if (FA_1) (*ThreadArgs)[i]->FA_1   = ObitFArrayRef(FA_1);
    else (*ThreadArgs)[i]->FA_1 = NULL;
    if (FA_2) (*ThreadArgs)[i]->FA_2   = ObitFArrayRef(FA_2);
    else (*ThreadArgs)[i]->FA_2 = NULL;
    if (FA_3) (*ThreadArgs)[i]->FA_3   = ObitFArrayRef(FA_3);
    else (*ThreadArgs)[i]->FA_3 = NULL;
    if (FA_4) (*ThreadArgs)[i]->FA_4   = ObitFArrayRef(FA_4);
    else (*ThreadArgs)[i]->FA_4 = NULL;
    (*ThreadArgs)[i]->first = 1;
    (*ThreadArgs)[i]->last  = size;
    (*ThreadArgs)[i]->value = 0.0;
    for (j=0; j<MAXCARRAYDIM; j++) (*ThreadArgs)[i]->pos[j] = 0;
    if (larg1>0) (*ThreadArgs)[i]->arg1 = g_malloc0(larg1*sizeof(float));
    else (*ThreadArgs)[i]->arg1 = NULL;
    if (larg2>0) (*ThreadArgs)[i]->arg2 = g_malloc0(larg2*sizeof(float));
    else (*ThreadArgs)[i]->arg2 = NULL;
    if (larg3>0) (*ThreadArgs)[i]->arg3 = g_malloc0(larg3*sizeof(float));
    else (*ThreadArgs)[i]->arg3 = NULL;
    if (larg4>0) (*ThreadArgs)[i]->arg4 = g_malloc0(larg4*sizeof(float));
    else (*ThreadArgs)[i]->arg4 = NULL;
    if (larg5>0) (*ThreadArgs)[i]->arg5 = g_malloc0(larg5*sizeof(float));
    else (*ThreadArgs)[i]->arg5 = NULL;
    if (larg5>0) (*ThreadArgs)[i]->arg6 = g_malloc0(larg6*sizeof(float));
    else (*ThreadArgs)[i]->arg6 = NULL;
    if (larg5>0) (*ThreadArgs)[i]->arg7 = g_malloc0(larg7*sizeof(float));
    else (*ThreadArgs)[i]->arg7 = NULL;
    (*ThreadArgs)[i]->ithread  = i;
  }

  return nThreads;
} /*  end MakeCAFuncArgs */

/**
 * Delete arguments for ThreadCAFunc
 * \param nargs      number of elements in ThreadArgs.
 * \param ThreadArgs Array of CAFuncArg
 */
static void KillCAFuncArgs (olong nargs, CAFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  ObitThreadPoolFree (ThreadArgs[0]->thread);  /* Free thread pool */
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      if (ThreadArgs[i]->thread) ObitThreadUnref(ThreadArgs[i]->thread);
      if (ThreadArgs[i]->in)   ObitCArrayUnref(ThreadArgs[i]->in);
      if (ThreadArgs[i]->in2)  ObitCArrayUnref(ThreadArgs[i]->in2);
      if (ThreadArgs[i]->out)  ObitCArrayUnref(ThreadArgs[i]->out);
      if (ThreadArgs[i]->FA_1) ObitFArrayUnref(ThreadArgs[i]->FA_1);
      if (ThreadArgs[i]->FA_2) ObitFArrayUnref(ThreadArgs[i]->FA_2);
      if (ThreadArgs[i]->FA_3) ObitFArrayUnref(ThreadArgs[i]->FA_3);
      if (ThreadArgs[i]->FA_4) ObitFArrayUnref(ThreadArgs[i]->FA_4);
      if (ThreadArgs[i]->arg1) g_free(ThreadArgs[i]->arg1);
      if (ThreadArgs[i]->arg2) g_free(ThreadArgs[i]->arg2);
      if (ThreadArgs[i]->arg3) g_free(ThreadArgs[i]->arg3);
      if (ThreadArgs[i]->arg4) g_free(ThreadArgs[i]->arg4);
      if (ThreadArgs[i]->arg5) g_free(ThreadArgs[i]->arg5);
      if (ThreadArgs[i]->arg6) g_free(ThreadArgs[i]->arg6);
      if (ThreadArgs[i]->arg7) g_free(ThreadArgs[i]->arg7);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillCAFuncArgs */
