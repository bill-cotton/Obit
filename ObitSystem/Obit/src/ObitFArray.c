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
/*;Correspondence about this software should be addressed as follows: */
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
typedef __m512  V16SF; // vector of 16 float (avx)
typedef __m512i V16SI; // vector of 16 int   (avx)
typedef __mmask16 MASK16; // vector of 16 mask  (avx)
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

#include "ObitThread.h"
#include "ObitFArray.h"
#include "ObitMem.h"
#include "ObitExp.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFArray.c
 * ObitFArray class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitFArray";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitFArrayClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitFArrayClassInfo myClassInfo = {FALSE};

/**
 * Init GSL random number generator
 */
#if HAVE_GSL==1  /* GSL stuff */
static gsl_rng *GSLran_gen=NULL;
#else
#include <stdlib.h>
#endif /* HAVE_GSL */

/*--------------- File Global Variables  ----------------*/


/*---------------Private structures----------------*/
/* Threaded function argument */
typedef struct {
  /* ObitThread to use */
  ObitThread *thread;
  /* ObitFArray to work on */
  ObitFArray *in;
  /* Second ObitFArray to work on */
  ObitFArray *in2;
  /* Output ObitFArray */
  ObitFArray *out;
  /* Additional ObitFArrays */
  ObitFArray *FA_3, *FA_4, *FA_5, *FA_6;
  /* First element (1-rel) number */
  olong        first;
  /* Highest element (1-rel) number */
  olong        last;
  /* Function dependent arguments */
  gpointer arg1, arg2, arg3, arg4, arg5, arg6, arg7;
  /* Return value */
  gfloat value;
  /* Return position */
  olong        pos[MAXFARRAYDIM];
  /* thread number  */
  olong        ithread;
} FAFuncArg;

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitFArrayInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitFArrayClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitFArrayClassInfoDefFn (gpointer inClass);

/** Private: Threaded Find max */
static gpointer ThreadFAMax (gpointer arg);

/** Private: Threaded Find min */
static gpointer ThreadFAMin (gpointer arg);

/** Private: Threaded Find abs max */
static gpointer ThreadFAAbsMax (gpointer arg);

/** Private: Threaded RMS  sums*/
static gpointer ThreadFARMSSum (gpointer arg);

/** Private: Threaded Accumulate histogram elements */
static gpointer ThreadFAHisto (gpointer arg);

/** Private: Threaded Convolve Gaussian */
static gpointer ThreadFAConvGaus (gpointer arg);

/** Private: Threaded Add */
static gpointer ThreadFAAdd (gpointer arg);

/** Private: Threaded Sub */
static gpointer ThreadFASub (gpointer arg);

/** Private: Threaded Multiply */
static gpointer ThreadFAMul (gpointer arg);

/** Private: Threaded Divide */
static gpointer ThreadFADiv (gpointer arg);

/** Private: Threaded SumArr */
static gpointer ThreadFASumArr (gpointer arg);

/** Private: Threaded AvgArr */
static gpointer ThreadFAAvgArr (gpointer arg);

/** Private: Threaded MaxArr */
static gpointer ThreadFAMaxArr (gpointer arg);

/** Private: Threaded MinArr */
static gpointer ThreadFAMinArr (gpointer arg);

/** Private: Threaded ExtArr */
static gpointer ThreadFAExtArr (gpointer arg);

/** Private: Threaded ShiftAdd */
static gpointer ThreadFAShAdd (gpointer arg);

/** Private: Threaded CplxSMulAccum */
static gpointer ThreadFACplxSMulAccum (gpointer arg);

/** Private: Threaded CplxMulAccum */
static gpointer ThreadFACplxMulAccum (gpointer arg);

/** Private: Make Threaded args */
static olong MakeFAFuncArgs (ObitThread *thread, ObitFArray *in,
			     ObitFArray *in2, ObitFArray *out,
			     ObitFArray *FA_3, ObitFArray *FA_4,
			     ObitFArray *FA_5, ObitFArray *FA_6,
			     olong larg1, olong larg2, olong larg3, 
			     olong larg4, olong larg5, 
			     olong larg6, olong larg7, 
			     FAFuncArg ***ThreadArgs);

/** Private: Delete Threaded args */
static void KillFAFuncArgs (olong nargs, FAFuncArg **ThreadArgs);


/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitFArray* newObitFArray (gchar* name)
{
  ObitFArray* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFArrayClassInit();

  /* allocate/init structure */
  out = ObitMemAlloc0Name (sizeof(ObitFArray), "ObitFArray");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitFArrayInit((gpointer)out);

 return out;
} /* end newObitFArray */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitFArrayGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFArrayClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitFArrayGetClass */

/**
 * Make a deep copy of an ObitFArray.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitFArray* ObitFArrayCopy  (ObitFArray *in, ObitFArray *out, ObitErr *err)
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
    out = newObitFArray(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */

  /* arrays */
  out = ObitFArrayRealloc (out, in->ndim, in->naxis);
  /* copy data */
  for (i=0; i<in->arraySize; i++) out->array[i] = in->array[i];

  return out;
} /* end ObitFArrayCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an FArray similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitFArrayClone  (ObitFArray *in, ObitFArray *out, ObitErr *err)
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

  /* arrays */
  out = ObitFArrayRealloc (out, in->ndim, in->naxis);

} /* end ObitFArrayClone */

/**
 * Determine if the two input objects have a compatable geometry.
 * Must have same number of non degenerate dimensions and 
 * each dimension must be the same size.
 * \param in1  First object to test.
 * \param in2  Second object to test.
 * \return TRUE if compatable, else FALSE.
 */
gboolean ObitFArrayIsCompatable (ObitFArray *in1, ObitFArray *in2)
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
} /* end ObitFArrayIsCompatable  */

/**
 * Creates an ObitFArray of a specified geometry.
 * \param name  An optional name for the object.
 * \param ndim  Number of dimensions desired, if <=0 data array not created.
 *              maximum value = MAXFARRAYDIM.
 * \param naxis Dimensionality along each axis. NULL => don't create array.
 * \return the new object.
 */
ObitFArray* ObitFArrayCreate (gchar* name, olong ndim, olong *naxis)
{
  ObitFArray* out;
  olong i, size;

  /* Create basic structure */
  out = newObitFArray (name);

  /* create data array if wanted */
  if ((ndim<0) || (naxis==NULL)) return out;
  g_assert (ndim<=MAXFARRAYDIM); /* sanity check */

  /* copy geometry */
  out->ndim = ndim;
  out->naxis = ObitMemAlloc0Name (ndim*sizeof(olong), "FArray naxis");
  if (ndim<=1) {  /* Single dimension */
    out->naxis[0] = MAX (1, naxis[0]);
    size = out->naxis[0]; /* total size */
  } else { /* Multi */
    size = 1; /* total size */
    for (i=0; i<ndim; i++) {
      /* WHAT??? out->naxis[i] = MAX (1, MIN(naxis[i],524288));  Not too big */
      out->naxis[i] = MAX (1, naxis[i]);
      size *= out->naxis[i]; /* total size */
    }
  }

  /* create array - add a bit extra, FFT seems to need it */
  out->array = ObitMemAlloc0Name (size*sizeof(ofloat) + 
				  out->naxis[0]*sizeof(ofloat),
				  "FArray array");
  out->arraySize = size;

  return out;
} /* end ObitFArrayCreate */

/**
 * Creates an ObitFArray of the specified subarray size of an extant
 * ObitFArray and copy the values.
 * \param in    Object with structures to subarray.
 * \param blc   (0-rel) lower index of first pixel to copy
 * \param trc   (0-rel) lower index of highest pixel to copy
 * \param err   Obit error stack object.
 * \return the new object.
 */
ObitFArray* ObitFArraySubArr (ObitFArray *in, olong *blc, olong *trc, 
			       ObitErr *err)
{
  olong ndim, idim, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10;
  olong naxis[MAXFARRAYDIM], ipos[MAXFARRAYDIM], opos[MAXFARRAYDIM];
  ObitFArray *out=NULL;
  gchar *outName;
  ofloat *inp, *outp;
  gchar *routine = "ObitFArraySubArr";

   /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (blc != NULL);
  g_assert (trc != NULL);

  /* Check window */
  Obit_retval_if_fail (((blc[0]>=0) && ((in->ndim==1) || (blc[1]>=0))), err, out,
		       "%s: Bad blc( %d,  %d)  in %s", 
		       routine, blc[0], blc[1], in->name);
  Obit_retval_if_fail (((trc[0]<=in->naxis[0]) && ((in->ndim==1) || (trc[1]<=in->naxis[1]))), 
		       err, out,
		       "%s: Bad trc( %d,  %d)  in %s", 
		       routine, trc[0], trc[1], in->name);

  /* Get size */
  ndim = in->ndim;
  for (idim=0; idim<MAXFARRAYDIM; idim++) naxis[idim] = 1;
  for (idim=0; idim<ndim; idim++) naxis[idim] = trc[idim]-blc[idim]+1;

  /* derive object name */
  outName = g_strconcat ("Subarray: ",in->name, NULL);
    
  /* Create output */
  out = ObitFArrayCreate (outName, ndim, naxis);
  g_free(outName); outName = NULL;

  /* Copy, loop over all possible dimensions */
  for (idim=0; idim<MAXFARRAYDIM; idim++) ipos[idim] = 0;
  for (i10=0; i10<naxis[9]; i10++) {
    opos[9] = i10;
    if(ndim>9) ipos[9] = blc[9] + i10;
    for (i9=0; i9<naxis[8]; i9++) {
      opos[8] = i9;
      if(ndim>8) ipos[8] = blc[8] + i9;
      for (i8=0; i8<naxis[7]; i8++) {
	opos[7] = i8;
	if(ndim>7) ipos[7] = blc[7] + i8;
	for (i7=0; i7<naxis[6]; i7++) {
	  opos[6] = i7;
	  if(ndim>6) ipos[6] = blc[6] + i7;
	  for (i6=0; i6<naxis[5]; i6++) {
	    opos[5] = i6;
	    if(ndim>5) ipos[5] = blc[5] + i6;
	    for (i5=0; i5<naxis[4]; i5++) {
	      opos[4] = i5;
	      if(ndim>4) ipos[4] = blc[4] + i5;
	      for (i4=0; i4<naxis[3]; i4++) {
		opos[3] = i4;
		if(ndim>3) ipos[3] = blc[3] + i4;
		for (i3=0; i3<naxis[2]; i3++) {
		  opos[2] = i3;
		  if(ndim>2) ipos[2] = blc[2] + i3;
		  for (i2=0; i2<naxis[1]; i2++) {
		    opos[1] = i2;
		    if(ndim>1) ipos[1] = blc[1] + i2;
		    opos[0] = 0;
		    ipos[0] = blc[0];
		    
		    /* Array pointers */
		    inp  = ObitFArrayIndex (in,  ipos);
		    outp = ObitFArrayIndex (out, opos);

		    /* Copy row */
		    for (i1=0; i1<naxis[0]; i1++) {
		      *outp++ = *inp++;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  return out;
} /* end  ObitFArraySubArr */

/**
 * Transpose an ObitFArray 
 * \param in    Object with structures to transpose
 * \param order output 1-rel order of the transposed axes, in storage order
 *              negative value = reverse order, 
 *              e,g, [2,1] = transpose 2D array
 * \param err   Obit error stack object.
 * \return the new object.
 */
ObitFArray* ObitFArrayTranspose (ObitFArray *in, olong *order,  ObitErr *err)
{
  olong ndim, idim, i2, i3, i4, i5, i6, i7, i8, i9, i10;
  olong incr, j1, j2, j3, j4, j5, j6, j7, j8, j9, j10;
  olong naxis[MAXFARRAYDIM], ipos[MAXFARRAYDIM], opos[MAXFARRAYDIM], 
    off[MAXFARRAYDIM], inAxis[MAXFARRAYDIM], i;
  gboolean flip[MAXFARRAYDIM];
  ObitFArray *out=NULL;
  gchar *outName;
  ofloat *inp, *outp;
  gchar *routine = "ObitFArrayTranspose";

   /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (order != NULL);


  /* Get size */
  ndim = in->ndim;
  for (idim=0; idim<MAXFARRAYDIM; idim++) {
    inAxis[idim] = 0;
    naxis[idim]  = -1;
    off[idim]    = 0;
  }
  for (idim=0; idim<ndim; idim++) {
    inAxis[idim] = abs(order[idim])-1;      /* 0-rel order */
    naxis[idim]  = in->naxis[inAxis[idim]];
  }
  for (idim=0; idim<ndim; idim++) {
    Obit_retval_if_fail ((naxis[idim]>0), err, out,
			 "%s: no order specified for dimension %d in %s", 
			 routine, idim+1, in->name);
  }
  for (idim=0; idim<MAXFARRAYDIM; idim++) naxis[idim] = MAX (naxis[idim], 1);

  /* calculate offsets */
  for (idim=0; idim<ndim; idim++) {
    flip[idim] = order[idim]<0;  /* Axis reversed in output? */
    off[idim] = 1;
    for (i=1; i<inAxis[idim]; i++) off[idim] *= in->naxis[i-1];
  }

  /* derive object name */
  outName = g_strconcat ("Transpose: ",in->name, NULL);
    
  /* Create output */
  out = ObitFArrayCreate (outName, ndim, naxis);
  g_free(outName); outName = NULL;

  /* Copy, loop over all possible dimensions - loop on output dimensions */
  for (idim=0; idim<MAXFARRAYDIM; idim++) ipos[idim] = 0;
  for (idim=0; idim<MAXFARRAYDIM; idim++) opos[idim] = -1;
  for (j10=0; j10<naxis[9]; j10++) {
    if (flip[9]) i10 = naxis[9]-j10-1;
    else         i10 = j10;
    opos[9] = j10;
    if(ndim>9) ipos[inAxis[9]] = i10;

    for (j9=0; j9<naxis[8]; j9++) {
      if (flip[8]) i9 = naxis[8]-j9-1;
      else         i9 = j9;
      opos[8] = j9;
      if(ndim>8) ipos[inAxis[8]] = i9;

      for (j8=0; j8<naxis[7]; j8++) {
	if (flip[7]) i8 = naxis[7]-j8-1;
	else         i8 = j8;
	opos[7] = j8;
	if(ndim>7) ipos[inAxis[7]] = i8;

	for (j7=0; j7<naxis[6]; j7++) {
	  if (flip[6]) i7 = naxis[6]-j7-1;
	  else         i7 = j7;
	  opos[6] = j7;
	  if(ndim>6) ipos[inAxis[6]] = i7;

	  for (j6=0; j6<naxis[5]; j6++) {
	    if (flip[5]) i6 = naxis[5]-j6-1;
	    else         i6 = j6;
	    opos[5] = j6;
	    if(ndim>5) ipos[inAxis[5]] = i6;

	    for (j5=0; j5<naxis[4]; j5++) {
	      if (flip[4]) i5 = naxis[4]-j5-1;
	      else         i5 = j5;
	      opos[4] = j5;
	      if(ndim>4) ipos[inAxis[4]] = i5;

	      for (j4=0; j4<naxis[3]; j4++) {
		if (flip[3]) i4 = naxis[3]-j4-1;
		else         i4 = j4;
		opos[3] = j4;
		if(ndim>3) ipos[inAxis[3]] = i4;

		for (j3=0; j3<naxis[2]; j3++) {
		  if (flip[2]) i3 = naxis[4]-j3-1;
		  else         i3 = j3;
		  opos[2] = j3;
		  if(ndim>2) ipos[inAxis[2]] = i3;

		  for (j2=0; j2<naxis[1]; j2++) {
		    if (flip[1]) i2 = naxis[1]-j2-1;
		    else         i2 = j2;
		    opos[1] = j2;
		    if(ndim>1) ipos[inAxis[1]] = i2;

		    /* where to start first axis */
		    opos[0] = 0;
		    if (flip[0]) {
		      incr = -off[0];
		      ipos[inAxis[0]] = naxis[0]-1;
		    } else {
		      incr = off[0];
		      ipos[inAxis[0]] = 0;
		    }

		    /* Output array pointera */
		    outp = ObitFArrayIndex (out, opos);
		    inp  = ObitFArrayIndex (in,  ipos);
		    /* DEBUG
		    fprintf (stderr, "ipos %d  %d %d opos %d  %d %d inAxis %d  %d  %d\n",
			     ipos[inAxis[0]], ipos[inAxis[1]], ipos[inAxis[2]], 
			     opos[0], opos[1],opos[2], inAxis[0], inAxis[1], inAxis[2]); */
		    
		    /* Copy row */
		    for (j1=0; j1<naxis[0]; j1++) {
		      *outp++ = *inp;
		      inp    += incr;  /* next element */
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  return out;
} /* end  ObitFArrayTranspose */

/**
 * Reallocate memory if needed, zero memory.
 * \param in    Object with structures to reallocate.
 * \param ndim  Number of dimensions desired, if <=0 data array not created.
 *              maximum value = #MAXFARRAYDIM.
 * \param naxis Dimensionality along each axis. NULL => don't create array.
 * \return the resized object.
 */
ObitFArray* ObitFArrayRealloc (ObitFArray* in, olong ndim, olong *naxis)
{
  ObitFArray* out = in;
  olong i, size;

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ndim<=MAXFARRAYDIM); /* sanity check */

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
    out->naxis = ObitMemAlloc0Name(ndim*sizeof(olong), "FArray naxis");
  }

  /* set dimensions, find output size */
  if (ndim<=1) {  /* Single dimension */
    out->naxis[0] = MAX (1, naxis[0]);
    size = out->naxis[0]; /* total size */
  } else { /* Multi */
    size = 1; /* total size */
    for (i=0; i<ndim; i++) {
      /* WHAT??? out->naxis[i] = MAX (1, MIN(naxis[i],524288));  Not too big */
      out->naxis[i] = MAX (1, naxis[i]);
      size *= out->naxis[i]; /* total size */
    }
  }
  
 /* resize array if needed */
  if (size != out->arraySize) {
    out->array = ObitMemRealloc(out->array, 
			   size*sizeof(ofloat)+out->naxis[0]*sizeof(ofloat));
    out->arraySize = size;
  }

  /* zero fill memory */
  memset (out->array, 0, size*sizeof(ofloat));

  return out;
} /* end ObitFArrayRealloc */

/**
 * Calculate offset for a given pixel location and return pointer.
 * Subsequent data are stored in order of increasing dimension 
 * (rows, then columns...).
 * \param in      Object with data
 * \param pos     array of 0-rel pixel numbers on each axis
 * \return pointer to specified cell; NULL if illegal pixel.
 */
ofloat*  ObitFArrayIndex (ObitFArray *in, olong *pos)
{
  ofloat *out = NULL;
  olong i, indx, previous;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (pos != NULL);
  g_assert (in->array != NULL);
  g_assert (in->naxis != NULL);

  /* Calculate offset */
  previous = 1;  /* size of previous dimension */
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
} /* end ObitFArrayIndex  */

/**
 * Find maximum pixel value.
 * Return value and location in pos.
 * \param in      Object with data
 * \param pos     (out) array of 0-rel pixel numbers on each axis
 * \return maximum value.
 */
ofloat ObitFArrayMax (ObitFArray *in, olong *pos)
{
  olong i, maxCell;
  ofloat maxVal, fblank = ObitMagicF();
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  gboolean OK;
  FAFuncArg **threadArgs;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (pos != NULL);

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in->thread, in, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
  OK = ObitThreadIterator (in->thread, nTh, 
			   (ObitThreadFunc)ThreadFAMax,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) return fblank;

  /* Find max */
  maxVal  = threadArgs[0]->value;
  maxCell = 0;
  for (i=1; i<nTh; i++) {
    if ((threadArgs[i]->value!=fblank) && (threadArgs[i]->value>maxVal)) {
	maxCell = i;
	maxVal  = threadArgs[i]->value;
      }
    }
  for (i=0; i<in->ndim; i++) pos[i] = threadArgs[maxCell]->pos[i];
    
  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
  
  return maxVal;
} /* end  ObitFArrayMax */

/**
 * Find maximum pixel value.
 * Return value and location in pos.
 * \param in      Object with data
 * \param pos     (out) array of 0-rel pixel numbers on each axis
 * \return maximum absolute (signed) value.
 */
ofloat ObitFArrayMaxAbs (ObitFArray *in, olong *pos)
{
  olong i, maxCell;
  ofloat maxAVal, fblank = ObitMagicF();
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  gboolean OK;
  FAFuncArg **threadArgs;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (pos != NULL);

   /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in->thread, in, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
  OK = ObitThreadIterator (in->thread, nTh, 
			   (ObitThreadFunc)ThreadFAAbsMax,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) return fblank;

  /* Find abs max */
  maxAVal = fabs(threadArgs[0]->value);
  maxCell = 0;
  for (i=1; i<nTh; i++) {
    if ((threadArgs[i]->value!=fblank) && (fabs(threadArgs[i]->value)>maxAVal)) {
	maxCell = i;
	maxAVal = fabs(threadArgs[i]->value);
      }
    }
  for (i=0; i<in->ndim; i++) pos[i] = threadArgs[maxCell]->pos[i];
    
  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);

  return MAX (0.0, maxAVal);
} /* end  ObitFArrayMaxAbs */

/**
 * Find minimum pixel value.
 * Return value and location in pos.
 * \param in      Object with data
 * \param pos     (out) array of 0-rel pixel numbers on each axis
 * \return minimum value.
 */
ofloat ObitFArrayMin (ObitFArray *in, olong *pos)
{
  olong i, minCell;
  ofloat minVal, fblank = ObitMagicF();
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  gboolean OK;
  FAFuncArg **threadArgs;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (pos != NULL);

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in->thread, in, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
  OK = ObitThreadIterator (in->thread, nTh, 
			   (ObitThreadFunc)ThreadFAMin,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) return fblank;

  /* Find min */
  minVal  = threadArgs[0]->value;
  minCell = 0;
  for (i=1; i<nTh; i++) {
    if ((threadArgs[i]->value!=fblank) && (threadArgs[i]->value<minVal)) {
      minCell = i;
      minVal  = threadArgs[i]->value;
    }
  }
  for (i=0; i<in->ndim; i++) pos[i] = threadArgs[minCell]->pos[i];
    

  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);

  return minVal;
} /* end  ObitFArrayMin */

/**
 * replace any magic value blanks with scalar
 * \param in      Object with data to deblank
 * \param scalar  Value to replace blanks.
 */
void ObitFArrayDeblank (ObitFArray *in, ofloat scalar)
{
  olong i, ilast;
  ofloat fblank = ObitMagicF();
#if HAVE_AVX512==1
  CV16SF v, vb, vs;
  MASK16 msk;
#elif HAVE_AVX==1
  CV8SF v, vb, vm, vs;
#endif

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Loop over array */
#if HAVE_AVX512==1  /* 16 float Vector */
  vb.v  = _mm512_set1_ps(fblank);   /* vector of blanks */
  vs.v  = _mm512_set1_ps(scalar);   /* vector of values to replace blank */
  /* Do blocks of 16 as vector */
  for (i=0; i<in->arraySize-16; i+=16) {
    v.v  = _mm512_loadu_ps(&in->array[i]);
    msk  = _mm512_cmp_ps_mask(v.v, vb.v, _CMP_EQ_OQ);  /* find blanks */
    v.v  = _mm512_mask_blend_ps(msk,v.v, vs.v);        /* replace blanks with scalar */
    _mm512_storeu_ps(&in->array[i], v.v);
  }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* 8 float Vector */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  vs.v  = _mm256_broadcast_ss(&scalar);   /* vector of values to replace blank */
  /* Do blocks of 8 as vector */
  for (i=0; i<in->arraySize-8; i+=8) {
    v.v  = _mm256_loadu_ps(&in->array[i]);
    vm.v = _mm256_cmp_ps(v.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v.v  = _mm256_blendv_ps(v.v,vs.v, vm.v);     /* replace blanks with scalar */
    _mm256_storeu_ps(&in->array[i], v.v);
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  for (i=ilast; i<in->arraySize; i++) 
      if (in->array[i]==fblank)in->array[i] = scalar;

} /* end  ObitFArrayDeblank */

/**
 *  Determine RMS noise in array.
 *  Value is based on a histogram analysis and is determined from 
 *  the width of the peak around the mode.
 *  out =  RMS (in.)
 * \param in Input object with data
 * \return rms of element distribution (-1 on error)
 */
ofloat ObitFArrayRMS (ObitFArray* in)
{
  olong i, j, modeCell=0, imHalf=0, ipHalf=0, numCell;
  olong i1, i2, ic, infcount;
  ofloat amax, amin, tmax, sum, sum2, x, count, mean, arg, cellFact=1.0;
  ofloat s, s2, half, *histo = NULL, *thist=NULL;
  ofloat rawRMS, rawMean, fiddle, out = -1.0, fblank = ObitMagicF();
  gboolean done = FALSE;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  ollong icount, c;
  gboolean OK;
  FAFuncArg **threadArgs;
 
  /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  /* Initialize Threading for initial values */
  nThreads = 
    MakeFAFuncArgs (in->thread, in, NULL, NULL, NULL, NULL, NULL, NULL, 
		    sizeof(ollong), sizeof(ofloat), sizeof(ofloat), 0, 0, 0, 0,
		    &threadArgs);
  
  /* Divide up work */
  nElem = in->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
  OK = ObitThreadIterator (in->thread, nTh, 
			   (ObitThreadFunc)ThreadFARMSSum,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) return fblank;

  /* sum parts */
  icount = 0; sum = sum2 = 0.0; 
  for (i=0; i<nTh; i++) {
    if (threadArgs[i]->value!=fblank) {
      c  = *(ollong*)(threadArgs[i]->arg1);
      s  = *(ofloat*)(threadArgs[i]->arg2);
      s2 = *(ofloat*)(threadArgs[i]->arg3);
      icount += c;
      sum    += s;
      sum2   += s2;
      /*fprintf (stderr,"%d %d %g %g %d %g %g \n",i,c,s,s2,count,sum,sum2); DEBUG*/
    }
  }

  /* cleanup */
  KillFAFuncArgs(nThreads, threadArgs);

  /* Initial values */
  count = (ofloat)icount;
  /* Better have something */
  if (count<5) return fblank;
  rawRMS = (sum2/icount) - ((sum / icount) * (sum / icount));
  if (rawRMS>0.0) rawRMS = sqrt(rawRMS);
  rawMean = sum / icount;
  /*fprintf(stderr,"Initial count %g rawRMS %g rawMean %g\n",count, rawRMS, rawMean);*/
  amin = rawMean - rawRMS;
  amax = rawMean + rawRMS;

  /* Make histogram size such that the average cell has at least 30 entries */
  numCell = icount / 30;
  numCell = MAX (100,  numCell);  /* but not too few */
  numCell = MIN (1000, numCell);  /* or too many */

  /* Initialize Threading for histogram */
  nThreads = 
    MakeFAFuncArgs (in->thread, in, NULL, NULL, NULL, NULL, NULL, NULL, 
		    sizeof(ollong), sizeof(ollong), numCell*sizeof(ofloat),sizeof(ofloat), sizeof(ofloat), 
		    sizeof(ollong), sizeof(ollong),
		    &threadArgs);
  
  /* Divide up work */
  nElem = in->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
  if (nElem<1000000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem;
    threadArgs[i]->last    = hiElem;
    memmove(threadArgs[i]->arg1, &numCell, sizeof(olong));
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Loop until a reasonable number of values in peak of histogram */
  infcount = 0;  /* Loop to check for endless loop */
  fiddle  = 1.0;  /* Factor to fiddle the half width */
  while (!done) {

    /* Don't do this forever */
    infcount++;
    fiddle *= 0.85;
    if (infcount>40) {
      KillFAFuncArgs(nThreads, threadArgs);
      return rawRMS;}  /* bag it */
    
    /* Set up thread arguments */
    for (i=0; i<nTh; i++) {
      memmove(threadArgs[i]->arg4, &amax, sizeof(ofloat));
      memmove(threadArgs[i]->arg5, &amin, sizeof(ofloat));
    }

    /* Do Form Histogram */
    OK = ObitThreadIterator (in->thread, nTh, 
			   (ObitThreadFunc)ThreadFAHisto,
			   (gpointer**)threadArgs);

    /* Check for problems */
    if (!OK) return fblank;

    /* Accumulate counts - histogram */
    icount = *(ollong*)(threadArgs[0])->arg2;
    histo = (ofloat*)(threadArgs[0]->arg3);
    for (i=1; i<nTh; i++) {\
      c     = *(ollong*)(threadArgs[i])->arg2;
     /*fprintf (stderr,"%d ncell %ld c %ld under %ld over %ld amax %g amin %g  \n",i, ncell, c, under, over, ttmax, ttmin);*/
      count += c;
      thist = (ofloat*)(threadArgs[i]->arg3);
      for (j=0; j<numCell; j++) histo[j] += thist[j];
    }
    /*fprintf (stderr,"DEBUG total count %ld  \n",icount);*/
    /* Find mode cell */
    cellFact =  numCell / (amax - amin + 1.0e-20);
    modeCell = -1;
    tmax = -1.0e20;
    for (i=1; i<numCell-1; i++) { /* ignore end cells */
      if (histo[i]>tmax) {
	tmax = histo[i];
	modeCell = i;
      }
    }
    
    /*fprintf (stderr,"DEBUG modeCell %d  cellFact %g\n",modeCell, cellFact);*/
    /*ttmax = histo[modeCell];
      fprintf (stderr,"%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f \n",
      histo[modeCell-5*2]/ttmax,histo[modeCell-4*2]/ttmax,histo[modeCell-3*2]/ttmax,histo[modeCell-2*2]/ttmax,histo[modeCell-1*2]/ttmax,
      histo[modeCell]/ttmax,
      histo[modeCell+1*2]/ttmax,histo[modeCell+2*2]/ttmax,histo[modeCell+3*2]/ttmax,histo[modeCell+4*2]/ttmax,histo[modeCell+5*2]/ttmax);
      
      fprintf (stderr,"%8.4g, %8.4g, %8.4g, %8.4g, %8.4g, %8.4g, %8.4g, %8.4g, %8.4g, %8.4g, %8.4g \n",
      histo[modeCell-5*2],histo[modeCell-4*2],histo[modeCell-3*2],histo[modeCell-2*2],histo[modeCell-1*2],
      histo[modeCell],
      histo[modeCell+1*2],histo[modeCell+2*2],histo[modeCell+3*2],histo[modeCell+4*2],histo[modeCell+5*2]);*/
    
    /* find half width of peak by finding number of cells positive 
       and negative from the mode the distribution stays above tmax/2.
       If the data are quantized then some (many) cells may have zero 
       occuputation (or nearly so) and should be ignored.
    */
    imHalf = modeCell;
    for (i=modeCell; i>=1; i--) {
      if (histo[i]>0.5*tmax)  imHalf = i;
      /*else if (histo[i]>0.2*tmax) break;*/
    }
    ipHalf = modeCell;
    for (i=modeCell; i<numCell-1; i++) {
      if (histo[i]>0.5*tmax) ipHalf = i;
      /*else if (histo[i]>0.2*tmax) break;*/
    }

    /*fprintf (stderr,"imHalf %d ipHalf %d diff %d\n",imHalf, ipHalf, ipHalf-imHalf); DEBUG*/
    /* acceptability tests */
    /* ipHalf - imHalf must be greater than 25 and less than 50 */
    if ((ipHalf-imHalf) < 25) {
      /* if peak completely unresolved */
      if ((ipHalf-imHalf)<=0) {  /* wild stab */
	half = fiddle*0.5 / cellFact; /* ~ halfwidth? 1/2 cell */
      } else { /* partly resolved */
	/*half = (fiddle*0.5 * (ipHalf-imHalf)) / cellFact;*/ /* ~ halfwidth */
 	half = (fiddle * (ipHalf-imHalf)) / cellFact; /* ~ halfwidth */
     }
      mean = amin + (modeCell-0.5) /  cellFact;
      /* don't spread over whole histogram */
      half *= numCell / 50.0; 
      /* Don't go below rawRMS */
      half = MAX(half,fiddle*rawRMS);
      /*fprintf (stderr,"DEBUG half %g fiddle %g rawRMS %g\n",half,fiddle,rawRMS);*/
      amax = mean + half;
      amin = mean - half;
      continue;  /* try it again */
    } else if ((ipHalf-imHalf) > 50) {  /* Spread too thinly? */
      mean = amin + (modeCell-0.5) /  cellFact;
      half = (0.5 * (ipHalf-imHalf)) / (fiddle*cellFact); /* ~ halfwidth */
      /* don't spread over whole histogram - try for 25 cells*/
      half *= numCell / 25.0;  
      /* Don't go below rawRMS */
      half = MAX(half,rawRMS);
      amax = mean + half;
      amin = mean - half;
      continue;  /* try it again */
    }
    done = TRUE;
    break;
  } /* end loop getting acceptable histogram */

  /* determine RMS */
  out = (ipHalf - imHalf) / cellFact; /* FWHM */
  out *= 2.35 / 2.0; /* to Sigma */

  /* get second moment around mode +/- 3 FWHM */
  i1 = modeCell - 3 * (modeCell-imHalf);
  i1 = MAX (0, i1);
  i2 = modeCell + 3 * (ipHalf - modeCell);
  i2 = MIN (numCell-1, i2);

  /* Get second moment */
  sum = sum2 = 0.0;
  count = 0.0; ic = 0;
  for (i=i1; i<=i2; i++) {
    if (histo[i]!=0) {
      x      = (ofloat)(i - modeCell);
      sum   += histo[i] * x;
      sum2  += histo[i] * x * x;
      count += histo[i];
      ic++;
    }
  }
    
  /* determine RMS */
  if (ic>2.0) {
    mean = sum / count;
    arg = (sum2/count) - mean*mean;
    if (arg>0.0) out = sqrt(arg);
    else out = 0.5 * fabs (mean);
    out /= cellFact;
  } else { /* trouble - histogram bad*/
    out = fblank;
  }
  /* debug fprintf (stderr,"infcount,%d\n",infcount);*/
  /* cleanup */
  KillFAFuncArgs(nThreads, threadArgs);

  return MIN (out, rawRMS);
} /* end  ObitFArrayRMS */

/**
 *  Determine RMS noise in array from average squared - average value
 *  out =  RMS (in.)
 * \param in Input object with data
 * \return rms of element distribution (-1 on error)
 */
ofloat ObitFArrayRawRMS (ObitFArray* in)
{
  ollong count;
  olong i, nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  ofloat sum, sum2, rawRMS, fblank = ObitMagicF();
  gboolean OK;
  FAFuncArg **threadArgs;
  ollong c;
  ofloat s,s2;

   /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  /* Initialize Threading */
  nThreads = 
    MakeFAFuncArgs (in->thread, in, NULL, NULL, NULL, NULL, NULL, NULL, 
		    sizeof(ollong), sizeof(ofloat), sizeof(ofloat), 0, 0, 0, 0,
		    &threadArgs);
  
  /* Divide up work */
  nElem = in->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
  OK = ObitThreadIterator (in->thread, nTh, 
			   (ObitThreadFunc)ThreadFARMSSum,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) return fblank;

  /* sum parts */
  count = 0; sum = sum2 = 0.0; 
  for (i=0; i<nTh; i++) {
    if (threadArgs[i]->value!=fblank) {
      c  = *(ollong*)(threadArgs[i]->arg1);
      s  = *(ofloat*)(threadArgs[i]->arg2);
      s2 = *(ofloat*)(threadArgs[i]->arg3);
      count += c;
      sum   += s;
      sum2  += s2;
      /*fprintf (stderr,"%d %d %g %g %d %g %g \n",i,c,s,s2,count,sum,sum2); DEBUG*/
    }
  }

  /* Get raw RMS */
  if (count>5)  {
    rawRMS = (sum2/count) - ((sum / count) * (sum / count));
    if (rawRMS>0.0) rawRMS = sqrt(rawRMS);
  } else rawRMS = fblank;

  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);

  return rawRMS;
} /* end  ObitFArrayRawRMS */

/**
 *  Determine RMS about zero in an array
 *  out =  RMS (in.)
 * \param in Input object with data
 * \return rms of element distribution (-1 on error)
 */
ofloat ObitFArrayRMS0 (ObitFArray* in)
{
  ollong count;
  olong i, nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  ofloat sum, sum2, rawRMS, fblank = ObitMagicF();
  gboolean OK;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  /* Initialize Threading */
  nThreads = 
    MakeFAFuncArgs (in->thread, in, NULL, NULL, NULL, NULL, NULL, NULL, 
		    sizeof(ollong), sizeof(ofloat), sizeof(ofloat), 0, 0, 0, 0,
		    &threadArgs);
  
  
  /* Divide up work */
  nElem = in->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
  OK = ObitThreadIterator (in->thread, nTh, 
			   (ObitThreadFunc)ThreadFARMSSum,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) return fblank;

  /* sum parts */
  count = 0; sum = sum2 = 0.0; 
  for (i=0; i<nTh; i++) {
    if (threadArgs[i]->value!=fblank) {
      count += *(ollong*)(threadArgs[i]->arg1);
      sum   += *(ofloat*)(threadArgs[i]->arg2);
      sum2  += *(ofloat*)(threadArgs[i]->arg3);
    }
  }

  /* Get RMS  about zero*/
  if (count>5)  {
    rawRMS = (sum2/count);
    if (rawRMS>0.0) rawRMS = sqrt(rawRMS);
  } else rawRMS = fblank;
  
  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);

  return rawRMS;
} /* end  ObitFArrayRMS0 */

/**
 *  Determine Histogram RMS in a possibly quantized image
 *  out =  RMS (in.)
 * If there are more than 100 quantization levels in the  RawRMS
 * then the distribution is assumed close enough to continuous and 
 * ObitFArrayRMS is called.
 * \param in Input object with data
 * \return rms of element distribution (-1 on error)
 */
ofloat ObitFArrayRMSQuant (ObitFArray* in)
{
  ofloat out = -1.0;
  olong i, ilast;
  ofloat quant, zero, center, rawRMS, mode, rangeFact, fblank = ObitMagicF();
  olong icell, modeCell=0, imHalf=0, ipHalf=0, numCell;
  olong i1, i2, ic, it;
  ofloat amin, tmax, sum, sum2, x, count, mean, arg, cellFact=1.0;
  ofloat *histo = NULL;
#if HAVE_AVX512==1
  CV16SF v, vb, vflag, vmin, vfac;
  IV16SF vcell;
  MASK16 msk;
  ofloat tmp;
  olong j;
#elif HAVE_AVX==1
  CV8SF v, vb, vm, vflag, vmin, vfac;
  IV8SF vcell;
  ofloat tmp;
  olong j;
#endif

   /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  /* Get quantization info */
  ObitFArrayQuant (in, &quant, &zero);

  /* Get Raw (unquantized ) RMS */
  rawRMS = ObitFArrayRMS(in);

  /* If more than 100 levels in Raw RMS assume close enough to continious
     distribution */
  if ((rawRMS/MAX(1.0e-20,quant))>100.0) return rawRMS;

  /* Find Mode */
  mode = ObitFArrayMode (in);

  /* Center the distribution of the quantization level closest to the mode */
  it = (mode - zero) / MAX(1.0e-20,quant);
  center = it * quant;

  /* Make histogram to cover +/- 2*rawRMS */
  rangeFact = 2.0;
  numCell = 2 * rangeFact * (rawRMS/MAX(1.0e-20,quant));
  if (numCell<10) return rawRMS;   /* Something to work with? */
  histo = ObitMemAllocName(numCell*sizeof(ofloat), "FArray histo");
  
  /* value range */
  /* amax = center + rangeFact*rawRMS; CHANGE amax not used */
  /* amin should be a quantization level */
  it = 0.5 + rangeFact*rawRMS/MAX(1.0e-20,quant);
  amin = center - it*quant;
  
  /* form histogram */
  cellFact =  1.0/quant;
  for (i=0; i<numCell; i++) histo[i] = 0.0;
#if HAVE_AVX512==1  /* 16 Vector */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  tmp = -cellFact*amin;
  vflag.v  = _mm512_set1_ps(tmp);  /* flagged values */
  tmp = amin;
  vmin.v  = _mm512_set1_ps(tmp);  /* amin */
  tmp = cellFact;
  vfac.v = _mm512_set1_ps(tmp);  /* fact */
  /* Do blocks of 16 as vector */
  for (i=0; i<in->arraySize-16; i+=16) {
    v.v    = _mm512_loadu_ps(&in->array[i]);         /* Data */
    msk    = _mm512_cmp_ps_mask(v.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v.v    = _mm512_mask_blend_ps(msk, v.v,vflag.v);    /* replace blanks */
    v.v    = _mm512_sub_ps(v.v,vmin.v);        /* -amin */
    v.v    = _mm512_mul_ps(v.v,vfac.v);        /* * fact */
    v.v    = _mm512_roundscale_ps(v.v,32);     /* round */
    vcell.v = _mm512_cvtps_epi32(v.v);
    for (j=0; j<16; j++) {
      if (vcell.i[j]<0.0) continue;
      icell = vcell.i[j];
      icell = MIN (numCell-1, MAX(0, icell));
      histo[icell]++;
    }
  }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* Vector */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  tmp = -cellFact*amin;
  vflag.v  = _mm256_broadcast_ss(&tmp);  /* flagged values */
  tmp = amin;
  vmin.v  = _mm256_broadcast_ss(&tmp);  /* amin */
  tmp = cellFact;
  vfac.v = _mm256_broadcast_ss(&tmp);  /* fact */
  /* Do blocks of 8 as vector */
  for (i=0; i<in->arraySize-8; i+=8) {
    v.v    = _mm256_loadu_ps(&in->array[i]);         /* Data */
    vm.v   = _mm256_cmp_ps(v.v, vb.v, _CMP_EQ_OQ);  /* find blanks */
    v.v    = _mm256_blendv_ps(v.v,vflag.v,vm.v);    /* replace blanks */
    v.v    = _mm256_sub_ps(v.v,vmin.v);        /* -amin */
    v.v    = _mm256_mul_ps(v.v,vfac.v);        /* * fact */
    v.v    = _mm256_round_ps(v.v,_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC);  /* round */
    vcell.v = _mm256_cvtps_epi32(v.v);
    for (j=0; j<8; j++) {
      if (vcell.i[j]<0.0) continue;
      icell = vcell.i[j];
      icell = MIN (numCell-1, MAX(0, icell));
      histo[icell]++;
    }
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar and rest */
  ilast = 0;
#endif
  for (i=ilast; i<in->arraySize; i++) {
    if (in->array[i]!=fblank){
      icell = 0.5 + cellFact * (in->array[i]-amin);
      icell = MIN (numCell-1, MAX(0, icell));
      histo[icell]++;
    }
  }
  
  /* Find mode cell */
  modeCell = -1;
  tmax = -1.0e20;
  for (i=1; i<numCell-1; i++) { /* ignore end cells */
    if (histo[i]>tmax) {
      tmax = histo[i];
      modeCell = i;
    }
  }
  
  
  /* find half width of peak by finding number of cells positive 
     and negative from the mode the distribution stays above tmax/2.
  */
  imHalf = modeCell;
  for (i=modeCell; i>=1; i--) {
    if (histo[i]>0.5*tmax)  imHalf = i;
    /*else if (histo[i]>0.2*tmax) break;*/
  }
  ipHalf = modeCell;
  for (i=modeCell; i<numCell-1; i++) {
    if (histo[i]>0.5*tmax) ipHalf = i;
    /*else if (histo[i]>0.2*tmax) break;*/
  }
  
  
  /* determine RMS */
  out = (ipHalf - imHalf) / cellFact; /* FWHM */
  out *= 2.35 / 2.0; /* Sigma */

  /* get second moment around mode +/- 3 sigma */
  i1 = modeCell - 3 * (modeCell-imHalf);
  i1 = MAX (0, i1);
  i2 = modeCell + 3 * (ipHalf - modeCell);
  i2 = MIN (numCell, i2);

  /* Get second moment */
  sum = sum2 = 0.0;
  count = 0.0; ic = 0;
  for (i=i1; i<=i2; i++) {
    if (histo[i]!=0) {
      x      = (ofloat)(i - modeCell);
      sum   += histo[i] * x;
      sum2  += histo[i] * x * x;
      count += histo[i];
      ic++;
    }
  }
    
  /* determine RMS */
  if (ic>2.0) {
    mean = sum / count;
    arg = (sum2/count) - mean*mean;
    if (arg>0.0) out = sqrt(arg);
    else out = 0.5 * fabs (mean);
    out /= cellFact;
  } else { /* trouble - histogram bad*/
    out = fblank;
  }
 
  /* cleanup */
  ObitMemFree(histo);

  return out;
} /* end  ObitFArrayRMSQuant */

/**
 *  Determine quantization from minimum difference between pixel values
 * and the value closest to zero;
 * \param in    Input object with data
 * \param quant [out] quantization level
 * \param zero  [out] closest level to zero
 */
void ObitFArrayQuant (ObitFArray* in, ofloat *quant, ofloat *zero)
{
  ofloat delta, small, fblank = ObitMagicF();
  olong i;
#if HAVE_AVX512==1
  CV16SF v[2], vbig[2], vzero, vt, vneg, vb, vacc, vsm;
  MASK16 msk;
  glong iv, io, is;
  ofloat big;
#elif HAVE_AVX==1
  CV8SF v[2], vbig[2], vzero, vt, vneg, vb, vm, vacc, vsm;
  glong iv, io, is;
  ofloat big;
#endif
   /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  delta = 1.0e10; small = 1.0e10;
#if HAVE_AVX512==1  /* Vector */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  vacc.v = _mm512_set1_ps(small);   /* min diff accumulator */
  vsm.v  = _mm512_set1_ps(small);   /* min accumulator */
  big = small;
  vbig[0].v = _mm512_set1_ps(big);   /* big values */
  big *= 1.05;
  vbig[1].v = _mm512_set1_ps(big);   /* bigger values */
  big = -small;
  vneg.v = _mm512_set1_ps(big);     /* big negative values */
  big = 0.0;
  vzero.v = _mm512_set1_ps(big);    /* zeroes */
  /* First vector load */
  iv = 0; is = 0;
  v[iv].v = _mm512_loadu_ps(&in->array[is]); is++;
  msk     = _mm512_cmp_ps_mask(v[iv].v, vb.v, _CMP_EQ_OQ);
  v[iv].v = _mm512_mask_blend_ps(msk, v[iv].v,vbig[iv].v); /* replace blanks */
  /* Loop over array differencing sliding blocks of 16 */
  while (is<in->arraySize-9) {
    iv = 1-iv; /* Swap */
    v[iv].v = _mm512_loadu_ps(&in->array[is]); is += 16;
    msk     = _mm512_cmp_ps_mask(v[iv].v, vb.v, _CMP_EQ_OQ);
    v[iv].v = _mm512_mask_blend_ps(msk,v[iv].v,vbig[iv].v); /* replace blanks */
    /* Max negative */
    msk     = _mm512_cmp_ps_mask(v[iv].v, vzero.v, _CMP_GT_OS);  /* mask of positives*/
    vt.v    = _mm512_mask_blend_ps(msk, v[iv].v, vneg.v);       /* set to large negative */
    vsm.v  = _mm512_max_ps(vsm.v, vt.v);            /* min accumulator */   
    /* Min positive */
    msk     = _mm512_cmp_ps_mask(v[iv].v, vzero.v, _CMP_LT_OS); /* mask of negatives */
    vt.v    = _mm512_mask_blend_ps(msk, v[iv].v, vbig[0].v);    /* set to large positives */
    vsm.v   = _mm512_min_ps(vsm.v, vt.v);           /* min accumulator */   
    /* Difference */
    io = 1-iv;  /* other */
    vt.v    = _mm512_sub_ps(v[iv].v,v[io].v);  /* Difference adjacent */
    vt.v    = _mm512_mul_ps(vt.v,vt.v);        /* No abs fn so square */
    vacc.v  = _mm512_min_ps(vacc.v, vt.v);     /* accumulate min */
  }
  /* parse results */
  delta = vacc.f[0]; 
  for (i=1; i<16; i++) if (vacc.f[i]<delta) delta = vacc.f[i];
  delta = sqrt(delta);  /* differences were squared */
  small = vsm.f[0];
  for (i=1; i<8; i++) if (fabs(vsm.f[i])<fabs(small)) small = vsm.f[i];
#elif HAVE_AVX==1  /* Vector */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  vacc.v = _mm256_broadcast_ss(&small);   /* min diff accumulator */
  vsm.v  = _mm256_broadcast_ss(&small);   /* min accumulator */
  big = small;
  vbig[0].v = _mm256_broadcast_ss(&big);   /* big values */
  big *= 1.05;
  vbig[1].v = _mm256_broadcast_ss(&big);   /* bigger values */
  big = -small;
  vneg.v = _mm256_broadcast_ss(&big);     /* big negative values */
  big = 0.0;
  vzero.v = _mm256_broadcast_ss(&big);    /* zeroes */
  /* First vector load */
  iv = 0; is = 0;
  v[iv].v = _mm256_loadu_ps(&in->array[is]); is++;
  vm .v   = _mm256_cmp_ps(v[iv].v, vb.v, _CMP_EQ_OQ);
  v[iv].v = _mm256_blendv_ps(v[iv].v,vbig[iv].v,vm.v); /* replace blanks */
  /* Loop over array differencing sliding blocks of 8 */
  while (is<in->arraySize-9) {
    iv = 1-iv; /* Swap */
    v[iv].v = _mm256_loadu_ps(&in->array[is]); is += 8;
    vm.v    = _mm256_cmp_ps(v[iv].v, vb.v, _CMP_EQ_OQ);
    v[iv].v = _mm256_blendv_ps(v[iv].v,vbig[iv].v,vm.v); /* replace blanks */
    /* Max negative */
    vm.v    = _mm256_cmp_ps(v[iv].v, vzero.v, _CMP_GT_OS);  /* mask of positives*/
    vt.v    = _mm256_blendv_ps(v[iv].v, vneg.v, vm.v); /* set to large negative */
    vsm.v  = _mm256_max_ps(vsm.v, vt.v);            /* min accumulator */   
    /* Min positive */
    vm.v    = _mm256_cmp_ps(v[iv].v, vzero.v, _CMP_LT_OS); /* mask of negatives */
    vt.v    = _mm256_blendv_ps(v[iv].v, vbig[0].v, vm.v);/* set to large positives */
    vsm.v   = _mm256_min_ps(vsm.v, vt.v);           /* min accumulator */   
    /* Difference */
    io = 1-iv;  /* other */
    vt.v    = _mm256_sub_ps(v[iv].v,v[io].v);  /* Difference adjacent */
    vt.v    = _mm256_mul_ps(vt.v,vt.v);        /* No abs fn so square */
    vacc.v  = _mm256_min_ps(vacc.v, vt.v);     /* accumulate min */
  }
  /* parse results */
  delta = vacc.f[0]; 
  for (i=1; i<8; i++) if (vacc.f[i]<delta) delta = vacc.f[i];
  delta = sqrt(delta);  /* differences were squared */
  small = vsm.f[0];
  for (i=1; i<8; i++) if (fabs(vsm.f[i])<fabs(small)) small = vsm.f[i];
# else /* scalar*/
  for (i=1; i<in->arraySize; i++) if (in->array[i]!=fblank) {
    if (in->array[i]!=in->array[i-1])
      delta = MIN (fabs(in->array[i]-in->array[i-1]), delta);
    if (fabs(in->array[i]) < fabs(small)) small = in->array[i];
  }
#endif

  /* Set output */
  *quant = delta;
  *zero  = small;
} /* end  ObitFArrayQuant */

/**
 *  Determine Mode of pixel value distribution in array.
 *  Value is based on a histogram analysis and is determined from 
 *  the peak in the distribution..
 *  out =  Mode (in.)
 * \param in Input object with data
 * \return mode of distribution
 */
ofloat ObitFArrayMode (ObitFArray* in)
{
  olong i, ilast, icell, modeCell, numCell, count, pos[MAXFARRAYDIM];
  ofloat amax, amin, tmax, cellFact, *histo;
  ofloat out = 0.0, fblank = ObitMagicF();
#if HAVE_AVX512==1
  CV16SF v, vb, vflag, vmin, vfac;
  IV16SF vcell;
  MASK16 msk;
  ofloat tmp;
  olong j;
#elif HAVE_AVX==1
  CV8SF v, vb, vm, vflag, vmin, vfac;
  IV8SF vcell;
  ofloat tmp;
  olong j;
#endif

   /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  /* get max/min */
  amax = ObitFArrayMax (in, pos);
  amin = ObitFArrayMin (in, pos);

   /* How many valid values? */
  count = 0;
  for (i=0; i<in->arraySize; i++) if (in->array[i]!=fblank) count++;
  if (count<5) return out; /* better have some */

 /* Make histogram size such that the average cell has 30 entries */
  numCell = count / 30;
  numCell = MAX (100, numCell);  /* but not too few */
  histo = ObitMemAllocName(numCell*sizeof(ofloat), "FArray histo");

  /* form histogram */
  cellFact =  numCell / (amax - amin + 1.0e-20);
  for (i=0; i<numCell; i++) histo[i] = 0.0;
#if HAVE_AVX512==1  /* 16 Vector */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  tmp = -cellFact*amin;
  vflag.v  = _mm512_set1_ps(tmp);  /* flagged values */
  tmp = amin;
  vmin.v  = _mm512_set1_ps(tmp);  /* amin */
  tmp = cellFact;
  vfac.v = _mm512_set1_ps(tmp);  /* fact */
  /* Do blocks of 16 as vector */
  for (i=0; i<in->arraySize-16; i+=16) {
    v.v    = _mm512_loadu_ps(&in->array[i]);        /* Data */
    msk    = _mm512_cmp_ps_mask(v.v, vb.v, _CMP_EQ_OQ);  /* find blanks */
    v.v    = _mm512_mask_blend_ps(msk, v.v,vflag.v);    /* replace blanks */
    v.v    = _mm512_sub_ps(v.v,vmin.v);        /* -amin */
    v.v    = _mm512_mul_ps(v.v,vfac.v);        /* * fact */
    v.v    = _mm512_roundscale_ps(v.v,32);     /* round */
    vcell.v = _mm512_cvtps_epi32(v.v);
    for (j=0; j<16; j++) {
      if (vcell.i[j]<0.0) continue;
      icell = vcell.i[j];
      icell = MIN (numCell-1, MAX(0, icell));
      histo[icell]++;
    }
  }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* Vector */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  tmp = -cellFact*amin;
  vflag.v  = _mm256_broadcast_ss(&tmp);  /* flagged values */
  tmp = amin;
  vmin.v  = _mm256_broadcast_ss(&tmp);  /* amin */
  tmp = cellFact;
  vfac.v = _mm256_broadcast_ss(&tmp);  /* fact */
  /* Do blocks of 8 as vector */
  for (i=0; i<in->arraySize-8; i+=8) {
    v.v    = _mm256_loadu_ps(&in->array[i]);         /* Data */
    vm.v   = _mm256_cmp_ps(v.v, vb.v, _CMP_EQ_OQ);  /* find blanks */
    v.v    = _mm256_blendv_ps(v.v,vflag.v,vm.v);    /* replace blanks */
    v.v    = _mm256_sub_ps(v.v,vmin.v);        /* -amin */
    v.v    = _mm256_mul_ps(v.v,vfac.v);        /* * fact */
    v.v    = _mm256_round_ps(v.v,_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC);  /* round */
    vcell.v = _mm256_cvtps_epi32(v.v);
    for (j=0; j<8; j++) {
      if (vcell.i[j]<0.0) continue;
      icell = vcell.i[j];
      icell = MIN (numCell-1, MAX(0, icell));
      histo[icell]++;
    }
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar and rest */
  ilast = 0;
#endif
  for (i=ilast; i<in->arraySize; i++) {
    if (in->array[i]!=fblank){
      icell = 0.5 + cellFact * (in->array[i]-amin);
      icell = MIN (numCell-1, MAX(0, icell));
      histo[icell]++;
    }
  }

  /* Find mode cell */
  modeCell = -1;
  tmax = -1.0e20;
  for (i=0; i<numCell; i++) {
    if (histo[i]>tmax) {
      tmax = histo[i];
      modeCell = i;
    }
  }

  /* convert from cell to value */
  out = amin + modeCell / cellFact;

  /* cleanup */
  ObitMemFree(histo);
  
  return out;
} /* end  ObitFArrayMode */

/**
 *  determine mean value in array
 *  out =  Mean (in.)
 * \param in Input object with data
 * \return mean of distribution
 */
ofloat ObitFArrayMean (ObitFArray* in)
{
  olong i, count, ilast;
  ofloat out = 0.0, fblank = ObitMagicF();
#if HAVE_AVX512==1
  CV16SF v, vb, vcnt,  vone, vzero, vt, vsum;
  MASK16 msk;
  ofloat tmp;
#elif HAVE_AVX==1
  CV8SF v, vb, vm, vcnt,  vone, vzero, vt, vsum;
  ofloat tmp;
#endif

   /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  out = 0.0;
  count = 0;
#if HAVE_AVX512==1  /* 16 Vector */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  tmp = 1.0;
  vone.v  = _mm512_set1_ps(tmp);   /* vector of ones */
  tmp = 0.0;
  vzero.v  = _mm512_set1_ps(tmp);  /* vector of zeroes */
  vcnt.v   = _mm512_set1_ps(tmp);  /* initialize counts */
  vsum.v   = _mm512_set1_ps(tmp);  /* initialize sum */
  /* Do blocks of 16 as vector */
  for (i=0; i<in->arraySize-16; i+=16) {
    v.v  = _mm512_loadu_ps(&in->array[i]);          /* load */
    msk  = _mm512_cmp_ps_mask(v.v, vb.v, _CMP_EQ_OQ);  /* find blanks */
    v.v  = _mm512_mask_blend_ps(msk,v.v,vzero.v);      /* replace blanks with zeroes */
    vt.v = _mm512_mask_blend_ps(msk,vone.v, vzero.v);  /* zero blank for counts */
    vcnt.v = _mm512_add_ps(vcnt.v,vt.v);  /* sum counts */
    vsum.v = _mm512_add_ps(vsum.v,v.v);   /* sum data */
  }
  ilast = i;  /* How far did I get? */
  vcnt.v = _mm512_roundscale_ps(vcnt.v,32);  /* round */
  for (i=0; i<16; i++) {
    count += (olong)vcnt.f[i];
    out   += vsum.f[i];
  }
#elif HAVE_AVX==1  /* Vector */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  tmp = 1.0;
  vone.v  = _mm256_broadcast_ss(&tmp);   /* vector of ones */
  tmp = 0.0;
  vzero.v  = _mm256_broadcast_ss(&tmp);  /* vector of zeroes */
  vcnt.v   = _mm256_broadcast_ss(&tmp);  /* initialize counts */
  vsum.v   = _mm256_broadcast_ss(&tmp);  /* initialize sum */
  /* Do blocks of 8 as vector */
  for (i=0; i<in->arraySize-8; i+=8) {
    v.v  = _mm256_loadu_ps(&in->array[i]);          /* load */
    vm.v = _mm256_cmp_ps(v.v, vb.v, _CMP_EQ_OQ);    /* find blanks */
    v.v  = _mm256_blendv_ps(v.v,vzero.v,vm.v);      /* replace blanks with zeroes */
    vt.v = _mm256_blendv_ps(vone.v, vzero.v, vm.v); /* zero blank for counts */
    vcnt.v = _mm256_add_ps(vcnt.v,vt.v);  /* sum counts */
    vsum.v = _mm256_add_ps(vsum.v,v.v);   /* sum data */
  }
  ilast = i;  /* How far did I get? */
  vcnt.v = _mm256_round_ps(vcnt.v,_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC);  /* round */
  for (i=0; i<8; i++) {
    count += (olong)vcnt.f[i];
    out   += vsum.f[i];
  }
 #else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  for (i=ilast; i<in->arraySize; i++) {
    if (in->array[i]!=fblank) {
      out += in->array[i];
      count++;
    }
  }
  
  if (count>0) out /= count;
  else out = fblank;

  return out;
} /* end  ObitFArrayMean */

/**
 *  Replace each element of the array with scalar.
 * in = scalar.
 * \param in Input object with data
 * \param scalar  Scalar value
 */
void ObitFArrayFill (ObitFArray* in, ofloat scalar)
{
  olong i, ilast;
#if HAVE_AVX512==1
  CV16SF vs;
#elif HAVE_AVX==1
  CV8SF vs;
#endif

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

#if HAVE_AVX512==1  /* Vector */
  vs.v   = _mm512_set1_ps(scalar);  /* vector of scalar */
  for (i=0; i<in->arraySize-16; i+=16) {
    _mm512_storeu_ps(&in->array[i], vs.v); /* Save */
 }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* Vector */
  vs.v   = _mm256_broadcast_ss(&scalar);  /* vector of scalar */
  for (i=0; i<in->arraySize-8; i+=8) {
    _mm256_storeu_ps(&in->array[i], vs.v); /* Save */
 }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  for (i=ilast; i<in->arraySize; i++) in->array[i] = scalar;
} /* end  ObitFArrayFill */

/**
 *  Negate each element of the array.
 * in = -in.
 * \param in Input object with data
 */
void ObitFArrayNeg (ObitFArray* in)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) 
    if (in->array[i]!=fblank) in->array[i] = -in->array[i];
} /* end  ObitFArrayNeg */

/**
 *  Takes abs value of each element of the array.
 * in = |in|.
 * \param in Input object with data
 */
void ObitFArrayAbs (ObitFArray* in)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) 
    if (in->array[i]!=fblank) in->array[i] = fabs(in->array[i]);
} /* end  ObitFArrayAbs */

/**
 *  Take sine of each element of the array.
 * in = sin(in).
 * \param in Input object with data
 */
void ObitFArraySin (ObitFArray* in)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) 
    if (in->array[i]!=fblank) in->array[i] = sin(in->array[i]);
} /* end  ObitFArraySin*/

/**
 *  Take cosine of each element of the array.
 * in = cos(in).
 * \param in Input object with data
 */
void ObitFArrayCos (ObitFArray* in)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) 
    if (in->array[i]!=fblank) in->array[i] = cos(in->array[i]);
} /* end  ObitFArrayCos */

/**
 *  Take square root of each element of the array.
 * in = sqrt(MAX(1.0e-20, in)).
 * \param in Input object with data
 */
void ObitFArraySqrt (ObitFArray* in)
{
  olong i, ilast;
  ofloat fblank = ObitMagicF();

#if HAVE_AVX512==1
  CV16SF vr, vb;
  MASK16 msk;
#endif

 /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

#if HAVE_AVX512==1  /* Vector length 16 */
  vb.v     = _mm512_set1_ps(fblank);  /* vector of blanks */
  for (i=0; i<in->arraySize-16; i+=16) {
    vr.v  = _mm512_loadu_ps (&in->array[i]);            /* Load Reals */
    msk   = _mm512_cmp_ps_mask(vr.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    vr.v  = _mm512_sqrt_ps(vr.v);                       /* sqrt */
    vr.v  = _mm512_mask_blend_ps(msk,vr.v,vb.v);        /* replace blanks */
    _mm512_storeu_ps(&in->array[i], vr.v);             /* Save */
  } /* end vector loop */
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  for (i=ilast; i<in->arraySize; i++) 
    if (in->array[i]!=fblank) 
      in->array[i] = sqrt(MAX(1.0e-20, in->array[i]));
} /* end  ObitFArraySqrt */

/**
 *  Sum each element of the array.
 * out =  Sum (in.)
 * \param in Input object with data
 * \return sum of elements 
 */
ofloat ObitFArraySum (ObitFArray* in)
{
  olong i, ilast=0;
  ofloat out = 0.0, fblank = ObitMagicF();

#if HAVE_AVX512==1
  CV16SF vr, vb, vz;
  MASK16 msk;
  ofloat tmp;
#endif
   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

#if HAVE_AVX512==1  /* Vector length 16 */
  vb.v     = _mm512_set1_ps(fblank);     /* vector of blanks */
  tmp = 0.0;
  vz.v  = _mm512_set1_ps(tmp);  /* vector of zeroes */
  for (ilast=0; ilast<in->arraySize-16; ilast+=16) {
    vr.v  = _mm512_loadu_ps (&in->array[ilast]);            /* Load Reals */
    msk   = _mm512_cmp_ps_mask(vr.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    vr.v  = _mm512_mask_blend_ps(msk,vr.v,vz.v);        /* replace blanks with zero */
    out += _mm512_reduce_add_ps(vr.v);                  /* Sum */
  } /* end vector loop */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  for (i=ilast; i<in->arraySize; i++) 
    if (in->array[i]!=fblank) out += in->array[i];

  return out;
} /* end  ObitFArraySum */

/**
 * How many valid elements are in in?
 * out =  Count of valid elements in (in.)
 * \param in Input object with data
 * \return count of valid elements 
 */
olong ObitFArrayCount (ObitFArray* in)
{
  olong i, ilast, out = 0;
  ofloat fblank = ObitMagicF();
  olong count=0;

#if HAVE_AVX512==1
  CV16SF v, vb, vzero, vone, vcnt, vt;
  MASK16 msk;
  ofloat tmp;
#elif HAVE_AVX==1
  CV8SF v, vb, vm, vzero, vone, vcnt, vt;
  ofloat tmp;
#endif
   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  count = 0;
#if HAVE_AVX512==1  /* Vector */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  tmp = 0.0;
  vzero.v  = _mm512_set1_ps(tmp);  /* vector of zeroes */
  vcnt.v   = _mm512_set1_ps(tmp);  /* initialize counts */
  tmp = 1.0;
  vone.v   = _mm512_set1_ps(tmp);  /* vector of ones */
  /* Do blocks of 16 as vector */
  for (i=0; i<in->arraySize-16; i+=16) {
    v.v  = _mm512_loadu_ps(&in->array[i]);
    msk  = _mm512_cmp_ps_mask(v.v, vb.v, _CMP_EQ_OQ);   /* find blanks */
    v.v  = _mm512_mask_blend_ps(msk,v.v,vzero.v);    /* replace blanks with zeroes */
    vt.v = _mm512_mask_blend_ps(msk,vone.v, vzero.v);
    vcnt.v = _mm512_add_ps(vcnt.v,vt.v);  /* sum counts */
  }
  ilast = i;  /* How far did I get? */
  vcnt.v = _mm512_roundscale_ps(vcnt.v,32);  /* round */
  for (i=0; i<16; i++) count += (olong)vcnt.f[i];
#elif HAVE_AVX==1  /* Vector */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  tmp = 0.0;
  vzero.v  = _mm256_broadcast_ss(&tmp);  /* vector of zeroes */
  vcnt.v   = _mm256_broadcast_ss(&tmp);  /* initialize counts */
  tmp = 1.0;
  vone.v   = _mm256_broadcast_ss(&tmp);  /* vector of ones */
  /* Do blocks of 8 as vector */
  for (i=0; i<in->arraySize-8; i+=8) {
    v.v  = _mm256_loadu_ps(&in->array[i]);
    vm.v = _mm256_cmp_ps(v.v, vb.v, _CMP_EQ_OQ);     /* find blanks */
    v.v  = _mm256_blendv_ps(v.v,vzero.v,vm.v); /* replace blanks with zeroes */
    vt.v = _mm256_blendv_ps(vone.v, vzero.v, vm.v);
    vcnt.v = _mm256_add_ps(vcnt.v,vt.v);  /* sum counts */
  }
  ilast = i;  /* How far did I get? */
  vcnt.v = _mm256_round_ps(vcnt.v,_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC);  /* round */
  for (i=0; i<8; i++) count += (olong)vcnt.f[i];
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  for (i=ilast; i<in->arraySize; i++) 
    if (in->array[i]!=fblank) count++;

  out = count;
  return out;
} /* end  ObitFArrayCount */

/**
 *  Add a scalar to each element of the array.
 * in = in + scalar
 * \param in      Input object with data
 * \param scalar  Scalar value
 */
void ObitFArraySAdd (ObitFArray* in, ofloat scalar)
{
  olong i, ilast;
  ofloat fblank = ObitMagicF();
#if HAVE_AVX512==1
  CV16SF v, vb, vs, vzero;
  MASK16 msk;
  ofloat tmp;
#elif HAVE_AVX==1
  CV8SF v, vb, vm, vs, vzero;
  ofloat tmp;
#endif

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

#if HAVE_AVX512==1  /* 16 Vector */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  vs.v   = _mm512_set1_ps(scalar);  /* vector of scalar */
  tmp = 0.0;
  vzero.v  = _mm512_set1_ps(tmp);  /* vector of zeroes */
  /* Do blocks of 16 as vector */
  for (i=0; i<in->arraySize-16; i+=16) {
    v.v  = _mm512_loadu_ps(&in->array[i]);       /* fetch */
    msk = _mm512_cmp_ps_mask(v.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v.v  = _mm512_mask_blend_ps(msk,v.v,vzero.v);   /* replace blanks with zeroes */
    v.v  = _mm512_add_ps(v.v,vs.v);                 /* do operation */
    v.v  = _mm512_mask_blend_ps(msk,v.v, vb.v);     /* reblank */
    _mm512_storeu_ps(&in->array[i], v.v);          /* Save */
 }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* Vector */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  vs.v   = _mm256_broadcast_ss(&scalar);  /* vector of scalar */
  tmp = 0.0;
  vzero.v  = _mm256_broadcast_ss(&tmp);  /* vector of zeroes */
  /* Do blocks of 8 as vector */
  for (i=0; i<in->arraySize-8; i+=8) {
    v.v  = _mm256_loadu_ps(&in->array[i]);       /* fetch */
    vm.v = _mm256_cmp_ps(v.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v.v  = _mm256_blendv_ps(v.v,vzero.v,vm.v);   /* replace blanks with zeroes */
    v.v  = _mm256_add_ps(v.v,vs.v);              /* do operation */
    v.v  = _mm256_blendv_ps(v.v, vb.v, vm.v);    /* reblank */
    _mm256_storeu_ps(&in->array[i], v.v);        /* Save */
 }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  for (i=ilast; i<in->arraySize; i++) 
    if (in->array[i]!=fblank) in->array[i] += scalar;
} /* end ObitFArraySAdd */

/**
 *  Multiply each element of the array by a scalar.
 * in = in * scalar
 * \param in      Input object with data
 * \param scalar  Scalar value
 */
void ObitFArraySMul (ObitFArray* in, ofloat scalar)
{
  olong i, ilast;
  ofloat fblank = ObitMagicF();
#if HAVE_AVX512==1
  CV16SF v, vb, vs, vone;
  MASK16 msk;
  ofloat tmp;
#elif HAVE_AVX==1
  CV8SF v, vb, vm, vs, vone;
  ofloat tmp;
#endif

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

#if HAVE_AVX512==1  /* Vector */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  vs.v   = _mm512_set1_ps(scalar);  /* vector of scalar */
  tmp = 1.0;
  vone.v  = _mm512_set1_ps(tmp);  /* vector of ones */
  /* Do blocks of 16 as vector */
  for (i=0; i<in->arraySize-16; i+=16) {
    v.v  = _mm512_loadu_ps(&in->array[i]);       /* fetch */
    msk  = _mm512_cmp_ps_mask(v.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v.v  = _mm512_mask_blend_ps(msk,v.v,vone.v);      /* replace blanks with ones */
    v.v  = _mm512_mul_ps(v.v,vs.v);                   /* do operation */
    v.v  = _mm512_mask_blend_ps(msk,v.v, vb.v);       /* reblank */
    _mm512_storeu_ps(&in->array[i], v.v);        /* Save */
 }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* Vector */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  vs.v   = _mm256_broadcast_ss(&scalar);  /* vector of scalar */
  tmp = 1.0;
  vone.v  = _mm256_broadcast_ss(&tmp);  /* vector of ones */
  /* Do blocks of 8 as vector */
  for (i=0; i<in->arraySize-8; i+=8) {
    v.v  = _mm256_loadu_ps(&in->array[i]);       /* fetch */
    vm.v = _mm256_cmp_ps(v.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v.v  = _mm256_blendv_ps(v.v,vone.v,vm.v);    /* replace blanks with ones */
    v.v  = _mm256_mul_ps(v.v,vs.v);              /* do operation */
    v.v  = _mm256_blendv_ps(v.v, vb.v, vm.v);    /* reblank */
    _mm256_storeu_ps(&in->array[i], v.v);        /* Save */
 }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  for (i=ilast; i<in->arraySize; i++) 
    if (in->array[i]!=fblank) in->array[i]  *= scalar;
} /* end ObitFArraySMul */

/**
 *  Divide each element of the array into a scalar.
 * No check for zeroes is made .
 * in = scalar / in
 * \param in      Input object with data
 * \param scalar  Scalar value
 */
void ObitFArraySDiv (ObitFArray* in, ofloat scalar)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) 
    if (in->array[i]!=fblank) in->array[i] = scalar / in->array[i];
} /* end ObitFArraySDiv */

/**
 *  Replace values outside of a given range with a new value
 * in = newVal where in < minVal or in > maxVal
 * \param in      Input object with data
 * \param minVal  Minimum allowed value
 * \param maxVal  Maximum allowed value
 * \param newVal  Value to use if out of range.
 */
void ObitFArrayClip (ObitFArray* in, ofloat minVal, ofloat maxVal, 
		     ofloat newVal)
{
  olong i;
  ofloat fblank = ObitMagicF();

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) 
    if ((in->array[i]!=fblank) && ((in->array[i]<minVal) || (in->array[i]>maxVal)))
      in->array[i] = newVal;
} /* end ObitFArrayClip */

/**
 *  Replace values inside of a given range with a new value
 * in = newVal where in >= minVal and in <= maxVal
 * \param in      Input object with data
 * \param minVal  Minimum allowed value
 * \param maxVal  Maximum allowed value
 * \param newVal  Value to use if out of range.
 */
void ObitFArrayInClip (ObitFArray* in, ofloat minVal, ofloat maxVal, 
		     ofloat newVal)
{
  olong i;
  ofloat fblank = ObitMagicF();

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) 
    if ((in->array[i]!=fblank) && ((in->array[i]>=minVal) && (in->array[i]<=maxVal)))
      in->array[i] = newVal;
} /* end ObitFArrayInClip */

/**
 *  Blank elements of array in1 where array in2 is blanked
 *  out = in1 or blank where in2 is blank
 * \param in1  Input object with data
 * \param in2  Input object with blanking
 * \param out  Output array (may be an input array).
 */
void ObitFArrayBlank (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i, ilast;
  ofloat fblank = ObitMagicF();
#if HAVE_AVX512==1
  CV16SF v, v2, vb;
  MASK16 msk;
#elif HAVE_AVX==1
  CV8SF v, v2, vb, vm;
#endif

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

#if HAVE_AVX512==1  /* 16 Vector */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  /* Do blocks of 16 as vector */
  for (i=0; i<in1->arraySize-16; i+=16) {
    v.v  = _mm512_loadu_ps(&in1->array[i]);       /* fetch 1 */
    v2.v = _mm512_loadu_ps(&in2->array[i]);       /* fetch 2 */
    msk  = _mm512_cmp_ps_mask(v2.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v.v  = _mm512_mask_blend_ps(msk,v.v, vb.v);     /* blank 1 */
    _mm512_storeu_ps(&out->array[i], v.v);        /* Save */
  }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* Vector */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  /* Do blocks of 8 as vector */
  for (i=0; i<in1->arraySize-8; i+=8) {
    v.v  = _mm256_loadu_ps(&in1->array[i]);       /* fetch 1 */
    v2.v = _mm256_loadu_ps(&in2->array[i]);       /* fetch 2 */
    vm.v = _mm256_cmp_ps(v2.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v.v  = _mm256_blendv_ps(v.v, vb.v, vm.v);     /* blank 1 */
    _mm256_storeu_ps(&out->array[i], v.v);        /* Save */
 }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  for (i=ilast; i<in1->arraySize; i++) 
    if (in2->array[i]!=fblank)
      out->array[i] = in1->array[i];
    else out->array[i] = fblank;
} /* end ObitFArrayBlank */

/**
 * Pick the larger nonblanked elements of two arrays.
 *  out = MAX (in1, in2) or whichever is not blanked
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArrayMaxArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in1->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
		       (ObitThreadFunc)ThreadFAMaxArr,
		       (gpointer**)threadArgs);
   
  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
  
} /* end ObitFArrayMaxArr */

/**
 * Pick the lesser nonblanked elements of two arrays.
 *  out = MIN (in1, in2) or whichever is not blanked
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArrayMinArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in1->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
		      (ObitThreadFunc)ThreadFAMinArr,
		      (gpointer**)threadArgs);

  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
} /* end ObitFArrayMinArr */

/**
 * Pick the more extreme (furthest from zero) nonblanked elements of two arrays.
 *  out =  Extreme (in1, in2) or whichever is not blanked
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArrayExtArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in1->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
		      (ObitThreadFunc)ThreadFAExtArr,
		      (gpointer**)threadArgs);

  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
} /* end ObitFArrayExtArr */

/**
 * Sum nonblanked elements of two arrays.
 *  out = (in1 + in2) or whichever is not blanked
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArraySumArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in1->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
		      (ObitThreadFunc)ThreadFASumArr,
		      (gpointer**)threadArgs);

  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
  
} /* end ObitFArraySumArr */

/**
 * Average nonblanked elements of two arrays.
 *  out = (in1 + in2)/2 or whichever is not blanked
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArrayAvgArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in1->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
		      (ObitThreadFunc)ThreadFAAvgArr,
		      (gpointer**)threadArgs);

  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
  
} /* end ObitFArrayAvgArr */

/**
 *  Add corresponding elements of two arrays.
 *  out = in1 + in2,  if either is blanked the result is blanked
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArrayAdd (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in1->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
		      (ObitThreadFunc)ThreadFAAdd,
		      (gpointer**)threadArgs);
  
  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
  
} /* end ObitFArrayAdd */

/**
 *  Abs add corresponding elements of two arrays.
 *  out = in1 + sign(1)*|in2|,  if either is blanked the result is blanked
 * \param in1  1st Input object with data
 * \param in2  2nd Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArrayAddAbs (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      if (out->array[i]>=0) 
	out->array[i] = in1->array[i] + fabs(in2->array[i]);
      else
	out->array[i] = in1->array[i] - fabs(in2->array[i]);
    else out->array[i] = fblank;
  }
} /* end ObitFArrayAddAbs */

/**
 *  Subtract corresponding elements of the arrays.
 *  out = in1 - in2, if either is blanked the result is blanked
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArraySub (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in1->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
		      (ObitThreadFunc)ThreadFASub,
		      (gpointer**)threadArgs);

  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
  
 } /* end ObitFArraySub */

/**
 *  Give the elements of one array the sign of another
 *  in2 = sign(in1)*|in2|, if either is blanked the result is blanked
 * \param in1  Input object with data
 * \param in2  Input object with data
 */
void ObitFArraySign (ObitFArray* in1, ObitFArray* in2)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) {
      if ((in1->array[i]>=0.0) && (in1->array[i]!=fblank))
	in2->array[i] =  fabs(in2->array[i]);
	else if (in2->array[i]!=fblank)
	  in2->array[i] = -fabs(in2->array[i]);
    } else in2->array[i] = fblank;
  }
 } /* end ObitFArraySign */

/**
 *  Multiply corresponding elements of the arrays.
 *  out = in1 * in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArrayMul (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in1->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
		      (ObitThreadFunc)ThreadFAMul,
		      (gpointer**)threadArgs);
  
  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
  
} /* end ObitFArrayMul */

/**
 *  Divide corresponding elements of the arrays.
 *  out = in1 / in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArrayDiv (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nElem = in1->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
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
		      (ObitThreadFunc)ThreadFADiv,
		      (gpointer**)threadArgs);

  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
  
} /* end ObitFArrayDiv */

/**
 *  Divide corresponding elements of the arrays with clipping.
 *  out = in1 / in2 where in2>minVal, else blanked
 * \param in1    Input object with data
 * \param in2    Input object with data
 * \param minVal minimum allowed value for in2
 * \param out    Output array (may be an input array).
 */
void ObitFArrayDivClip (ObitFArray* in1, ObitFArray* in2, ofloat minVal, ObitFArray* out)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank) && 
	(in2->array[i]>minVal)) {
      out->array[i] = in1->array[i] / in2->array[i];
      /* debug
      if (out->array[i] >1.0e5) fprintf (stderr,"bad div %d %g %g %g\n",
					  i,in1->array[i],in2->array[i],out->array[i]); */
    } else out->array[i] = fblank;
  }
} /* end ObitFArrayDiv */

/**
 *  Sum the products of the elements of two arrays
 *  out = Sum (in1 x in2)
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \return sum of product of elements 
 */
ofloat ObitFArrayDot (ObitFArray* in1, ObitFArray* in2)
{
  olong i;
  ofloat fblank = ObitMagicF();
  ofloat out = 0.0;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out += in1->array[i] * in2->array[i];
  }
  return out;
} /* end ObitFArrayDot */

/**
 *  Multiply the elements of a 2D array by the corresponding elements
 *  of a row and column vector.
 *  NOTE: this does not check for magic value blanking, this was causing 
 *  trouble in its major application - image formation - which should not 
 *  produce blanks.
 *  out[i,j] = in[i,j] * row[j] * col[i].
 * \param in   Input 2D array
 * \param row  Input row vector
 * \param col  Input column
 * \param out  Output array (may be an input array).
 */
void ObitFArrayMulColRow (ObitFArray* in, ObitFArray* row, ObitFArray* col,
			  ObitFArray* out)
{
  olong i, j, nx, ny;
  ofloat *inp, *outp, *rowp, *colp;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(row, &myClassInfo));
  g_assert (ObitIsA(col, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in, out));
  g_assert (in->ndim==2);
  g_assert (row->ndim==1);
  g_assert (col->ndim==1);
  g_assert (row->naxis[0]==in->naxis[0]);
  g_assert (col->naxis[0]==in->naxis[1]);

  /* Dimensions */
  nx = in->naxis[0];
  ny = in->naxis[1];

  /* init pointers */
  inp  = in->array;
  rowp = row->array;
  colp = col->array;
  outp = out->array;

  /* Double loop over array */
  for (j=0; j<ny; j++) {
    for (i=0; i<nx; i++) {
      outp[i] = inp[i] * rowp[j] * colp[i];
    } /* end row loop */
    inp  += nx; /* Pointers for next row */
    outp += nx;
  } /* end column loop */

} /* end ObitFArrayMulColRow  */

/**
 *  In-place rearrangement of a center-at-the edges array to 
 * center at the center, or the other way around.
 * This is needed for the peculiar order of FFTs.
 * FFTs don't like blanked values.
 * \param in   1D array to reorder
 */
void ObitFArray1DCenter (ObitFArray* in)
{
  olong i, pos[1], nx;
  ofloat temp, *inp, *outp;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* ... and God sayeth unto FFT: 
     "The middle shall be first and the last shall be the middle." */

  /* swap halves */
  nx = in->naxis[0];
  pos[0] = 0;
  inp = ObitFArrayIndex (in, pos);
  pos[0] = nx/2;
  outp = ObitFArrayIndex (in, pos);

  for (i=0; i<nx/2; i++) {
    temp    = outp[i];
    outp[i] = inp[i];
    inp[i]  = temp;
  }
} /* end ObitFArray1DCenter */

/**
 *  In-place rearrangement of a center-at-the edges array to 
 * center at the center, or the other way around.
 * This is needed for the peculiar order of FFTs.
 * FFTs don't like blanked values.
 * \param in   2D array to reorder
 */
void ObitFArray2DCenter (ObitFArray* in)
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
  inp = ObitFArrayIndex (in, pos);
  pos[0] = nx/2; pos[1] = ny/2;
  outp = ObitFArrayIndex (in, pos);

  for (j=0; j<ny/2; j++) {
    for (i=0; i<nx/2; i++) {
      temp = outp[i];
      outp[i] = inp[i];
      inp[i] = temp;
    }
    /* next rwo */
    inp += nx;
    outp += nx;
  }

  /* swap second and third quadrents */
  nx = in->naxis[0];
  ny = in->naxis[1];
  pos[0] = nx/2; pos[1] = 0;
  inp = ObitFArrayIndex (in, pos);
  pos[0] = 0; pos[1] = ny/2;
  outp = ObitFArrayIndex (in, pos);

  for (j=0; j<ny/2; j++) {
    for (i=0; i<nx/2; i++) {
      temp = outp[i];
      outp[i] = inp[i];
      inp[i] = temp;
    }
    /* next rwo */
    inp += nx;
    outp += nx;
  }
  
} /* end ObitFArray2DCenter */

/**
 *  In-place inversion of a symmetric matrix in an ObitFArray
 * Adopted from AIPS SYMINV.FOR
 * Magic blanking not supported
 * \param in   2D array to invers (max dim 50x50)
 * \param ierr return code, 0=>OK, else could not invert.
 */
void ObitFArray2DSymInv (ObitFArray* in, olong *ierr)
{
#define MAXSYMINVSIZE 50
  olong  j, k=0, l, m, n;
  ofloat  ab, big, p[MAXSYMINVSIZE], q[MAXSYMINVSIZE], r[MAXSYMINVSIZE];
  ofloat *a;
 

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  *ierr = 0;
  if (in->ndim!=2) *ierr = 1;
  if (in->naxis[0]>50) *ierr = 2;
  if (in->naxis[0]!=in->naxis[1]) *ierr = 3;
  if (*ierr != 0) return;

  /* array information */
  a = in->array;
  n = in->naxis[0];

  /* make sure symmetric */
  for (j= 1; j<=n; j++) { /* loop 10 */
    for (k= 1; k<=n; k++) { /* loop 5 */
      /* not symmetric */
      if (abs (a[(k-1)*n+(j-1)]-a[(j-1)*n+(k-1)]) > abs (a[(k-1)*n+(j-1)])*1.0e-4) {
	*ierr = 4;
	return;
      } 
    } /* end loop  L5:   */;
  } /* end loop  L10:  */;
  for (m= 1; m<=n; m++) { /* loop 15 */
    r[m-1] = 1.0;
  } /* end loop  L15:  */;
  for (m= 1; m<=n; m++) { /* loop 200 */
    big = 0.0;
    for (l= 1; l<=n; l++) { /* loop 20 */
      ab = abs (a[(l-1)*n+(l-1)]);
      if ((ab > big)  &&  (r[l-1] != 0)) {
	big = ab;
	k = l;
      } 
    } /* end loop  L20:  */;
    /* inverse indeterminant? */
    if (big == 0.0) {
      *ierr = 5;
      return;
    } 
    r[k-1] = 0.0;
    q[k-1] = 1.0 / a[(k-1)*n+(k-1)];
    p[k-1] = 1.0;
    a[(k-1)*n+(k-1)] = 0.0;
    for (l= 1; l<=k-1; l++) { /* loop 30 */
      p[l-1] = a[(k-1)*n+(l-1)];
      if (r[l-1] == 0.0) {
	q[l-1] = a[(k-1)*n+(l-1)] * q[k-1];
      } else {
	q[l-1] = -a[(k-1)*n+(l-1)] * q[k-1];
      } 
      a[(k-1)*n+(l-1)] = 0.0;
    } /* end loop  L30:  */;
    for (l= k+1; l<=n; l++) { /* loop 40 */
      if (r[l-1] != 0.0) {
	p[l-1] = a[(l-1)*n+(k-1)];
      } else {
	p[l-1] = -a[(l-1)*n+(k-1)];
      } 
      q[l-1] = -a[(l-1)*n+(k-1)] * q[k-1];
      a[(l-1)*n+(k-1)] = 0.0;
    } /* end loop  L40:  */;
    for (l= 1; l<=n; l++) { /* loop 60 */
      for (k= l; k<=n; k++) { /* loop 50 */
	a[(k-1)*n+(l-1)] = a[(k-1)*n+(l-1)] + p[l-1]*q[k-1];
      } /* end loop  L50:  */;
    } /* end loop  L60:  */;
  } /* end loop  L200: */;
  /* fill in symmetric half */
  m = n + 1;
  l = n;
  for (k= 2; k<=n; k++) { /* loop 80 */
    m = m - 1;
    l = l - 1;
    for (j= 1; j<=l; j++) { /* loop 70 */
      a[(j-1)*n+(m-1)] = a[(m-1)*n+(j-1)];
    } /* end loop  L70:  */;
  } /* end loop  L80:  */;
  
 } /* end ObitFArray2DSymInv */


/**
 * Make 2-D circular Gaussian in an FArray
 * \param array  Array to fill in
 * \param Cen    0-rel center pixel
 * \param Size   FWHM of Gaussian in pixels.
 */
void ObitFArray2DCGauss (ObitFArray *array, olong Cen[2], ofloat FWHM)
{
  olong ix, iy, indx, nx, ny, size[2];
  ofloat x, y, *data, sigma, factor, arg;

  /* error checks */
  g_assert(ObitFArrayIsA(array));

  nx = array->naxis[0];
  ny = array->naxis[1];
  size[0] = size[1] = 0;
  data = ObitFArrayIndex (array, size);

  /* Loop */
  sigma = FWHM / 2.3548;
  factor = 1.0 / (2.0 * sigma * sigma);
  for (iy=0; iy<ny; iy++) {
    for (ix=0; ix<nx; ix++) {
      indx = iy*nx + ix;
      x = (ofloat)(ix - Cen[0]);
      y = (ofloat)(iy - Cen[1]);
      arg = -(x*x + y*y) * factor;
      if (arg>-15.0) arg = ObitExpCalc (arg);
      else arg = 0.0;
      data[indx] = arg;
    }
  }

} /* end ObitFArray2DCGauss */

/**
 * Make 2-D elliptical Gaussian in an FArray
 * Model is added to previous contents of the FArray
 * \param array   Array to fill in
 * \param amp     Peak value of Gaussian
 * \param Cen     0-rel center pixel
 * \param GauMod  Gaussian parameters, Major axis, FWHM, minor axis 
 *        FWHM (both in pixels) and rotation angle wrt "X" axis (deg).
 */
void ObitFArray2DEGauss (ObitFArray *array, ofloat amp, ofloat Cen[2], ofloat GauMod[3])
{
  olong ix, iy, indx, nx, ny, size[2];
  ofloat *data, sigmaX, sigmaY, factorX, factorY, arg;
  ofloat CosPA, SinPA, xp, yp, x, y;

  /* error checks */
  g_assert(ObitFArrayIsA(array));

  nx = array->naxis[0];
  ny = array->naxis[1];
  size[0] = size[1] = 0;
  data = ObitFArrayIndex (array, size);

  /* Set values for calculating Gaussian */
  CosPA = cos ((GauMod[2]+90.0)*DG2RAD);
  SinPA = sin ((GauMod[2]+90.0)*DG2RAD);
  sigmaX = GauMod[1] / 2.3548;   /* X is minor */
  factorX = 1.0 / (2.0 * sigmaX * sigmaX);
  sigmaY = GauMod[0] / 2.3548;   /* Y is major */
  factorY = 1.0 / (2.0 * sigmaY * sigmaY);

  /* Loop */
  for (iy=0; iy<ny; iy++) {
    for (ix=0; ix<nx; ix++) {
      indx = iy*nx + ix;
      xp = (ofloat)(ix - Cen[0]);
      yp = (ofloat)(iy - Cen[1]);
      /* Rotate to elipse frame */
      x = CosPA*xp + yp*SinPA;
      y = CosPA*yp - xp*SinPA;
      arg = -(x*x*factorX + y*y*factorY);
      if (arg>-15.0) {
	arg = amp * ObitExpCalc (arg);
	data[indx] += arg;
      }
    }
  }

} /* end ObitFArray2DEGauss */

/**
 *  Shift and Add scaled arrays
 *  Two FArrays are aligned at specified pixels and the 
 *  corresponding pixels are added with a scalar multiplied 
 *  times the second.
 *  Only handles to 3 dimensions.
 *  If in1/out are 3D and in2 is 2D then the same plane in in2
 *  is used for all planes in in1/out.
 * NB: this works better if the alignment point is near the center of in2
 *  out = in1 + scalar x in2 in overlap, else in1
 * \param in1     First input object with data, may be blanked
 * \param pos1    Alignment pixel in in1 (0-rel)
 * \param in2     Second input object with data, blanked pixels ignored
 * \param pos2    Alignment pixel in in2  (0-rel)
 * \param scalar  factor to be multiplied times in2
 * \param out     Output array, may be an input array and MUST 
 *                have the same the same geometry.
 */
void ObitFArrayShiftAdd (ObitFArray* in1, olong *pos1, 
			 ObitFArray* in2, olong *pos2, 
			 ofloat scalar, ObitFArray* out)
{
  olong ip, np, hiy, loy, ny1, ny2, offy;
  olong i, nTh, nRow, loRow, hiRow, nRowPerThread, nThreads;
  gboolean OK;
  FAFuncArg **threadArgs;

  /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));
  g_assert(pos1!=NULL);
  g_assert(pos2!=NULL);
  /* Check dimensionality */
  g_assert (in1->ndim==in2->ndim);
  g_assert (in1->ndim==out->ndim);
  g_assert (in1->naxis[0]==out->naxis[0]);
  g_assert (in1->naxis[1]==out->naxis[1]);

  if (in1->ndim>2) np = in1->naxis[2];  /* NUmber of planes */
  else np = 1;

  /* determine regions of overlap in in1/out of in2 */
  offy = pos2[1] - pos1[1];
  ny1 = in1->naxis[1];
  ny2 = in2->naxis[1];
  loy = MAX (0, pos1[1] - in2->naxis[1]/2);
  hiy = MIN (ny1-1, pos1[1] + in2->naxis[1]/2);
  /* In case in2 not centered */
  if (loy+offy<0) loy -= loy+offy;
  if (hiy+offy>=ny2) hiy -= hiy+offy-ny2+1;
  /* Keep in range */
  loy = MAX (0, loy);
  hiy = MIN (hiy, ny1-1);

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in1->thread, in1, in2, out, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work */
  nRow = hiy - loy + 1;          /* How many rows? */
  nRowPerThread = nRow/nThreads;
  nTh = nThreads;
  if (nRow<5*nThreads) {nRowPerThread = nRow; nTh = 1;}
  loRow = loy;
  hiRow = loRow + nRowPerThread-1;
  hiRow = MIN (hiRow, hiy);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    threadArgs[i]->value = scalar; /* Scaling factor */
    threadArgs[i]->pos[0] = 0;     /* Plane */
    threadArgs[i]->pos[1] = pos1[0]; threadArgs[i]->pos[2] = pos1[1];   /* Alignment */
    threadArgs[i]->pos[3] = pos2[0]; threadArgs[i]->pos[4] = pos2[1];   /* Alignment */
    if (i==(nTh-1)) hiRow  = hiy;  /* Make sure do all */
    threadArgs[i]->first   = loRow;
    threadArgs[i]->last    = hiRow;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Row */
    loRow += nRowPerThread;
    hiRow += nRowPerThread;
    hiRow = MIN (hiRow, ny1);
  }

  /* Loop over planes */
  for (ip = 0; ip<np; ip++) {

    for (i=0; i<nTh; i++) threadArgs[i]->pos[0] = ip;    /* Plane */

    /* Do operation */
    OK = ObitThreadIterator (in1->thread, nTh, 
			     (ObitThreadFunc)ThreadFAShAdd,
			     (gpointer**)threadArgs);
    if (!OK) goto cleanup;
  } /* end loop over planes */

  /* Free local objects */
 cleanup:
  KillFAFuncArgs(nThreads, threadArgs);
  
} /* end ObitFArrayShiftAdd */

/**
 * Zero fills out and inserts in, centered and multiplied by factor.
 * Any blanks in in are replaced with zero.
 * This routine is intended for zero padding images before an FFT
 * May run threaded
 * to increase the resolution in the uv plane.
 * \param in      Object with structures to zero pad
 * \param out     Output object
 * \param factor  scaling factor for in
 */
void ObitFArrayPad (ObitFArray *in, ObitFArray* out, ofloat factor)
{
  olong ipos[2], opos[2];

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* Zero fill output */
  ObitFArrayFill (out, 0.0);

  /* Insert in into out - center as well as possible */
  ipos[0] = in->naxis[0]/2;  ipos[1] = in->naxis[1]/2;
  opos[0] = out->naxis[0]/2; opos[1] = out->naxis[1]/2;
  ObitFArrayShiftAdd (out, opos, in, ipos, factor, out);

  /* Deblank */
   ObitFArrayDeblank (out, 0.0);

} /* end  ObitFArrayPad */

/**
 * Convolves a list of points with a Gaussian and adds them to a 2-D array.
 * \param in      2-D array to add Gaussians to.
 * \param list    List of positions and fluxes of the Gaussians
 *                (x pixel, y pixel, flux)
 * \param ncomp   Number of components in list  
 *                (generally less than size of FArray).
 * \param gauss   Gaussian coefficients for (d_x*d_x, d_y*d_y, d_x*d_y)
 *                Gaussian maj = major axis FWHM, min=minor, pa = posn. angle
 *                cr=cos(pa+rotation), sr=sin(pa+rotation),
 *                cell_x, cell_y x, y cell spacing is same units as maj, min
 *                [0] = {(cr/min)^2 + ((sr/maj)^2)}*(cell_x^2)*4*log(2)
 *                [1] = {(sr/min)^2 + ((cr/maj)^2)}*(cell_y^2)*4*log(2)
 *                [2] = {(1/min)^2 - ((1/maj)^2)}*sr*cr*abs(cell_x*cell_y)*8*log(2)
 */
void  ObitFArrayConvGaus (ObitFArray* in, ObitFArray* list, olong ncomp, 
			  ofloat gauss[3])
{

  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  gboolean OK;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(list, &myClassInfo));

  /* anything to do? */
  if (ncomp<=0) return;
  ObitExpCalc(0.0);  /* Make sure initialized */


  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in->thread, in, NULL, NULL, NULL, NULL, NULL, NULL, 
			     0, sizeof(olong), 3*sizeof(ofloat), 0, 0, 0, 0,
			     &threadArgs);
  
  /* Divide up work by row - only thread if more than 100 rows */
  nElem = in->naxis[1];
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
  if (nElem<100) {nElemPerThread = nElem; nTh = 1;}
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
    threadArgs[i]->arg1 = (gpointer)list;
    memmove(threadArgs[i]->arg2, &ncomp, sizeof(olong));
    memmove(threadArgs[i]->arg3, gauss, 3*sizeof(ofloat));
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  OK = ObitThreadIterator (in->thread, nTh, 
			   (ObitThreadFunc)ThreadFAConvGaus,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) return;

  /* Reset arg1 - keep KillFAFuncArgs from zapping what it's pointing at */
  for (i=0; i<nTh; i++) {
    threadArgs[i]->arg1 = NULL;
  }

  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
} /* end ObitFArrayConvGaus */

/**
 * Select elements in an FArray by increment
 * \param in    Input Object
 * \param out   Output Object
 * \param blc   (0-rel) lower index of first pixel to copy
 * \param trc   (0-rel) lower index of highest pixel to copy
 * \param inc   increment on each axis
 * \param err   Obit error stack object.
 */
void ObitFArraySelInc (ObitFArray *in, ObitFArray *out, olong *blc, olong *trc, 
		       olong *inc, ObitErr *err)
{
  olong ndim, idim, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10;
  olong naxis[MAXFARRAYDIM], linc[MAXFARRAYDIM], ipos[MAXFARRAYDIM], opos[MAXFARRAYDIM];
  ofloat *inp, *outp, tmp;
  gchar *routine = "ObitFArraySelInc";

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));
  g_assert (blc != NULL);
  g_assert (trc != NULL);
  g_assert (inc != NULL);

  /* Check window */
  Obit_return_if_fail (((blc[0]>=0) && ((in->ndim==1) || (blc[1]>=0))), err,
		       "%s: Bad blc( %d,  %d)  in %s", 
		       routine, blc[0], blc[1], in->name);
  Obit_return_if_fail (((trc[0]<=in->naxis[0]) && ((in->ndim==1) || (trc[1]<=in->naxis[1]))), 
		       err,
		       "%s: Bad trc( %d,  %d)  in %s", 
		       routine, trc[0], trc[1], in->name);

  /* Get size */
  ndim = in->ndim;
  for (idim=0; idim<MAXFARRAYDIM; idim++) naxis[idim] = 1;
  for (idim=0; idim<MAXFARRAYDIM; idim++) linc[idim]  = 1;
  for (idim=0; idim<ndim; idim++) {
    naxis[idim] = in->naxis[idim];
    linc[idim]  = inc[idim];
    /* Check dimension */
    tmp = 0.99 + inc[idim] * ((ofloat)naxis[idim] / (ofloat)inc[idim]);
    Obit_return_if_fail (((olong)tmp ==naxis[idim]), err, 
		       "%s: incompatable dimension %d ( %d,  %d)  in %s", 
		       routine, idim, (olong)tmp, naxis[idim], in->name);
  }

  /* Copy, loop over all possible dimensions */
  for (idim=0; idim<MAXFARRAYDIM; idim++) ipos[idim] = 0;
  for (idim=0; idim<MAXFARRAYDIM; idim++) opos[idim] = -1;
  for (i10=0; i10<naxis[9]; i10+=linc[9]) {
    opos[9]++;
    if(ndim>9) ipos[9] = blc[9] + i10;
    if (ipos[9]>=trc[9]) continue;
    for (i9=0; i9<naxis[8]; i9+=linc[8]) {
      opos[8]++;
      if(ndim>8) ipos[8] = blc[8] + i9;
      if (ipos[8]>=trc[8]) continue;
      for (i8=0; i8<naxis[7]; i8+=linc[7]) {
	opos[7]++;
	if(ndim>7) ipos[7] = blc[7] + i8;
	if (ipos[7]>=trc[7]) continue;
	for (i7=0; i7<naxis[6]; i7+=linc[6]) {
	  opos[6]++;
	  if(ndim>6) ipos[6] = blc[6] + i7;
	  if (ipos[6]>=trc[6]) continue;
	  for (i6=0; i6<naxis[5]; i6+=linc[5]) {
	    opos[5]++;
	    if(ndim>5) ipos[5] = blc[5] + i6;
	    if (ipos[5]>=trc[5]) continue;
	    for (i5=0; i5<naxis[4]; i5+=linc[4]) {
	      opos[4]++;
	      if(ndim>4) ipos[4] = blc[4] + i5;
	      if (ipos[4]>=trc[4]) continue;
	      for (i4=0; i4<naxis[3]; i4+=linc[3]) {
		opos[3]++;
		if(ndim>3) ipos[3] = blc[3] + i4;
		if (ipos[3]>=trc[3]) continue;
		for (i3=0; i3<naxis[2]; i3+=linc[2]) {
		  opos[2]++;
		  if(ndim>2) ipos[2] = blc[2] + i3;
		  if (ipos[2]>=trc[2]) continue;
		  for (i2=0; i2<naxis[1]; i2+=linc[1]) {
		    opos[1]++;
		    if(ndim>1) ipos[1] = blc[1] + i2;
		    if (ipos[1]>=trc[1]) continue;
		    opos[0] = 0;
		    ipos[0] = blc[0];
		    
		    /* Array pointers */
		    inp  = ObitFArrayIndex (in,  ipos);
		    outp = ObitFArrayIndex (out, opos);

		    /* Copy row */
		    for (i1=0; i1<naxis[0]; i1+=linc[0]) {
		      if (blc[0]+i1>=trc[0]) continue;
		      *outp++ = inp[i1];
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  return;
} /* end  ObitFArraySelInc */

/**
 * Return histogram of elements in an FArray
 * \param in    Input Object
 * \param n     Number of elements in histogram
 * \param min   Min value in histogram
 * \param max   Max value in histogram
 * \return FArray with histogram, info has items:
 * \li nHisto OBIT_long scalar Number of elements in histogram
 * \li Min    OBIT_float scalar Minimum value in histogram
 * \li Max    OBIT_float scalar Maximum value in histogram
 * \li Total  OBIT_float scalar Total number of values in histogram
 * \li Under  OBIT_float scalar Number of underflows in histogram
 * \li Over   OBIT_float scalar Number of overflows in histogram
 */
ObitFArray*  ObitFArrayHisto (ObitFArray* in, olong n, ofloat min, ofloat max)
{
  ObitFArray *out=NULL;
  olong i, j, ndim, naxis[1];
  ollong cntUnder, cntOver, cntTotal;
  ofloat  *histo = NULL, *thist=NULL, ftemp;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  gboolean OK;
  gint32  dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  /* Create output */
  ndim = 1; naxis[0] = n;
  out = ObitFArrayCreate ("Histogram", ndim, naxis);

  /* Initialize Threading */
  nThreads = 
    MakeFAFuncArgs (in->thread, in, NULL, NULL, NULL, NULL, NULL, NULL, 
		    sizeof(ollong), sizeof(olong), n*sizeof(ofloat), sizeof(ofloat), sizeof(ofloat), 
		    sizeof(ollong), sizeof(ollong),
		    &threadArgs);
  
  /* Divide up work */
  nElem = in->arraySize;
  nElemPerThread = nElem/nThreads;
  nTh = nThreads;
  if (nElem<100000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);
  
  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem;
    threadArgs[i]->last    = hiElem;
    memmove(threadArgs[i]->arg1, &n, sizeof(olong));
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }
  
  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    memmove(threadArgs[i]->arg4, &max, sizeof(ofloat));
    memmove(threadArgs[i]->arg5, &min, sizeof(ofloat));
  }
  
  /* Do Form Histogram */
  OK = ObitThreadIterator (in->thread, nTh, 
			   (ObitThreadFunc)ThreadFAHisto,
			   (gpointer**)threadArgs);
  
  /* Check for problems */
  if (!OK) return out;
  
  /* Accumulate counts - histogram */
  cntTotal = *(olong*)(threadArgs[0])->arg2;
  cntUnder = *(olong*)(threadArgs[0])->arg6;
  cntOver  = *(olong*)(threadArgs[0])->arg7;
  histo = (ofloat*)out->array;
  memmove(histo, threadArgs[0]->arg3, n*sizeof(ofloat));
  for (i=1; i<nTh; i++) {
    cntTotal += *(olong*)(threadArgs[i])->arg2;
    cntUnder += *(olong*)(threadArgs[i])->arg6;
    cntOver  += *(olong*)(threadArgs[i])->arg7;
    thist = (ofloat*)(threadArgs[i]->arg3);
    for (j=0; j<n; j++) histo[j] += thist[j];
  }
  
  /* cleanup */
  KillFAFuncArgs(nThreads, threadArgs);
  
  /* Save info*/
  ObitInfoListAlwaysPut(out->info, "nHisto", OBIT_long,  dim, &n);
  ObitInfoListAlwaysPut(out->info, "Min",    OBIT_float, dim, &min);
  ObitInfoListAlwaysPut(out->info, "Max",    OBIT_float, dim, &max);
  ftemp = (ofloat)cntTotal;
  ObitInfoListAlwaysPut(out->info, "Total",  OBIT_float, dim, &ftemp);
  ftemp = (ofloat)cntUnder;
  ObitInfoListAlwaysPut(out->info, "Under",  OBIT_float, dim, &ftemp);
  ftemp = (ofloat)cntOver;
  ObitInfoListAlwaysPut(out->info, "Over",   OBIT_float, dim, &ftemp);
  
  return out;
} /* end ObitFArrayHisto  */

/**
 * Exponentiate each element of the array.
 * out = exp(in).
 * \param in  Input object with data
 * \param out Output object with data
 */
void ObitFArrayExp (ObitFArray* in, ObitFArray* out)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);
  g_assert (ObitFArrayIsCompatable(in, out));

  for (i=0; i<in->arraySize; i++) {
    if (in->array[i]!=fblank) out->array[i] = expf(in->array[i]);
    else                      out->array[i] = fblank;
  }
} /* end  ObitFArrayExp */

/**
 * Natural log of each element of the array.
 * out = ln(in). out blank where in==0
 * \param in  Input object with data
 * \param out Output object with data
 */
void ObitFArrayLog (ObitFArray* in, ObitFArray* out)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);
  g_assert (ObitFArrayIsCompatable(in, out));

  for (i=0; i<in->arraySize; i++) {
    if ((in->array[i]!=fblank) && (in->array[i]!=0.0)) 
      out->array[i] = logf(in->array[i]);
    else out->array[i] = fblank;
  }
} /* end  ObitFArrayLog */

/**
 * Raise in1 to the in2 power for each element of the array.
 * out = pow (in1, in2) = in1^in2. if in1<0 out is blanked
 * \param in1  1st input object with data, 
 * \param in2  2nd input object with data
 * \param out Output object with data
 */
void ObitFArrayPow (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (in1->array != NULL);
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank) && (in1->array[i]>0.))
      out->array[i] = powf(in1->array[i], in2->array[i]);
    else out->array[i] = fblank;
  }
} /* end  ObitFArrayPow */

/**
 * Generate Gaussian distributed random number
 * NB: this is stateful, DO NOT MULTITHREAD!
 * \param mean of distribution
 * \param sigma of distribution
 * \return random number
 */
ofloat ObitFArrayRandom (ofloat mean, ofloat sigma)
{
  ofloat val = mean;
#if HAVE_GSL==1  /* GSL stuff */
  odouble dsigma = sigma;
  if (GSLran_gen==NULL) GSLran_gen=gsl_rng_alloc(gsl_rng_default);
  val = mean + (ofloat)gsl_ran_gaussian (GSLran_gen, dsigma);
#else /* more primitive */
  int i;
  ofloat sum, norm;
  sum = 0.0;
  norm = 1.0 / RAND_MAX;
  for (i=0; i<12; i++) sum += (norm*rand());
  val = sigma * (sum - 6.0) + mean;
#endif /* HAVE_GSL */
  return val;
} /* end ObitFArrayRandom */

/**
 * Replace each element of the array with Gaussian random numbers.
 * \param in Input object 
 * \param mean of distribution
 * \param sigma of distribution
 */
void ObitFArrayRandomFill (ObitFArray* in, ofloat mean, ofloat sigma)
{
  olong i;

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) in->array[i] = ObitFArrayRandom(mean, sigma);
} /* end  ObitFArrayRandomFill */

/**
 * Replace each element of a rectangular subarray with value, 2D Only
 * \param in     Input 2D object 
 * \param win    blc,trc 0-rel
 * \param value  new value
 */
void ObitFArrayRectFill (ObitFArray* in, olong win[4], ofloat value)
{
  olong i, j, k, nx, ny, x_lo, x_hi, y_lo, y_hi, lwin[4];

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);
  nx = in->naxis[0]; ny = in->naxis[1];
  /* get in correct order */
  for (i=0; i<4; i++) lwin[i] = win[i];
  if ((win[0]>win[2]) || (win[1]>win[3])) {
    lwin[0]=win[2]; lwin[1]=win[3]; lwin[2]=win[0]; lwin[3]=win[1]; 
  }
  x_lo = MAX(0,lwin[0]); x_hi = MIN(nx,lwin[2]); 
  y_lo = MAX(0,lwin[1]); y_hi = MIN(ny,lwin[3]); 
  for (j=y_lo; j<y_hi; j++) {
    for (i=x_lo; i<x_hi; i++) {
      k = j*nx + i;
      in->array[k] = value;
    } /* end inner loop */
  } /* end outer loop */
} /* end  ObitFArrayRectFill */

/**
 * Replace each element of a round subarray with value, 2D Only
 * \param in     Input object 
 * \param win    radius, center 0-rel
 * \param value  new value
 */
void ObitFArrayRoundFill (ObitFArray* in, olong win[3], ofloat value)
{
  olong i, j, k, nx, ny, x_lo, x_hi, y_lo, y_hi;
  ofloat rad2, dist2, dx2, dy2;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);
  nx = in->naxis[0]; ny = in->naxis[1];
  x_lo = MAX(0,(win[1]-win[0])); x_hi = MIN(nx,(win[1]+win[0])); 
  y_lo = MAX(0,(win[2]-win[0])); y_hi = MIN(ny,(win[2]+win[0])); 
  /*fprintf(stderr,"RectFill x_lo %d x_hi %d y_lo %d y_hi %d \n", x_lo,x_hi,y_lo,y_hi);*/
  rad2 = (ofloat)(win[0])*(ofloat)(win[0]); /* max distance^2 */
  for (j=y_lo; j<y_hi; j++) {
    dy2 = (ofloat)(j-win[2])*(ofloat)(j-win[2]);
    for (i=x_lo; i<x_hi; i++) {
      dx2 = (ofloat)(i-win[1])*(ofloat)(i-win[1]);
      dist2 = dx2+dy2;  /* How far? */
      if (dist2<=rad2) {
	k = j*nx + i;
	in->array[k] = value;
      } /* end in circle */
    } /* end inner loop */
  } /* end outer loop */
} /* end  ObitFArrayRoundFill */

/**
 * If elements in in are >= those in out1 then out1 is set to the value in in 
 * and the corresponding element in out2 is set to value.
 * Blanking NOT supported
 * \param in     Input object 
 * \param out1   First output
 * \param value  new value for out2
 * \param out2   Second output whose elements are set to value when
 *               an element in out1 is updated.
 */
void ObitFArrayMaxSetValue (ObitFArray* in, ObitFArray* out1, ofloat value, 
			    ObitFArray* out2)
{
  olong i, ilast;
#if HAVE_AVX512==1
  CV16SF v1, v2, v3, vs;
  MASK16 msk;
#elif HAVE_AVX==1
  CV8SF v1, v2, v3, vs, vm;
#endif

  /* error checks */
  g_assert (ObitFArrayIsCompatable(in, out1));
  g_assert (ObitFArrayIsCompatable(out1, out2));

#if HAVE_AVX512==1  /* AVX512 Vector */
  vs.v   = _mm512_set1_ps(value);  /* vector of value */
  /* Do blocks of 16 as vector */
  for (i=0; i<in->arraySize-16; i+=16) {
    v1.v = _mm512_loadu_ps(&in->array[i]);       /* Load from in */
    v2.v = _mm512_loadu_ps(&out1->array[i]);     /* Load from out1 */
    msk  = _mm512_cmp_ps_mask(v1.v, v2.v, _CMP_GE_OQ); /* find new max */
    v2.v = _mm512_mask_blend_ps(msk,v2.v,v1.v);  /* Select out1 */
    _mm512_storeu_ps(&out1->array[i], v2.v);      /* Save out1 */
    v3.v = _mm512_loadu_ps(&out2->array[i]);     /* Load from out2 */
    v3.v = _mm512_mask_blend_ps(msk,v3.v,vs.v);  /* Select out2 */
    _mm512_storeu_ps(&out2->array[i], v3.v);      /* Save out2 */
  }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* AVX Vector */
  vs.v   = _mm256_broadcast_ss(&value);  /* vector of value */
  /* Do blocks of 8 as vector */
  for (i=0; i<in->arraySize-8; i+=8) {
    v1.v = _mm256_loadu_ps(&in->array[i]);
    v2.v = _mm256_loadu_ps(&out1->array[i]);
    vm.v = _mm256_cmp_ps(v1.v, v2.v, _CMP_GE_OQ); /* find new max */
    v2.v = _mm256_blendv_ps(v2.v,v1.v, vm.v);  /* Select out1 */
    _mm256_storeu_ps(&out1->array[i], v2.v);    /* Save out1 */
    v3.v = _mm256_loadu_ps(&out2->array[i]);
    v3.v = _mm256_blendv_ps(v3.v,vs.v, vm.v);  /* replace blanks 2 */
    _mm256_storeu_ps(&out2->array[i], v3.v);    /* Save out2 */
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  for (i=ilast; i<in->arraySize; i++) {
    if (in->array[i]>out1->array[i]) {
      out1->array[i] = in->array[i];
      out2->array[i] = value;
    }
  } /* end loop over array */
} /* end ObitFArrayMaxSetValue */

/*---------Complex functions using pairs of FArrays -----------------*/
/**
 * Amplitude squares from a pair of reals 
 * out = Fin_r**2 + Fin_i**2
 * Blanking NOT supported
 * \param Fin_r  Real part of input
 * \param Fin_i  Imaginary part of input
 * \param out    Output object
 */
void ObitFArrayCplxAmp2 (ObitFArray* Fin_r, ObitFArray* Fin_i, 
			 ObitFArray* out)
{
  olong i, ilast=0;
  ofloat *iArr_r = Fin_r->array, *iArr_i = Fin_i->array;
  ofloat *oArr = out->array;
#if HAVE_AVX512==1
  CV16SF v1r, v1i, tv1;
#elif HAVE_AVX==1
  CV8SF  v1r, v1i, tv1;
#endif

   /* error checks */
  g_assert (ObitFArrayIsCompatable(Fin_r, out));
  g_assert (ObitFArrayIsCompatable(Fin_i, out));

#if HAVE_AVX512==1  /* AVX512 Vector */
  /* Do blocks of 16 as vector */
  for (i=0; i<Fin_r->arraySize-16; i+=16) {
    v1r.v = _mm512_loadu_ps(&iArr_r[i]);   /* input real */
    v1r.v = _mm512_mul_ps(v1r.v, v1r.v);   /* square real*/
    v1i.v = _mm512_loadu_ps(&iArr_i[i]);   /* input imaginary */
    v1i.v = _mm512_mul_ps(v1i.v, v1i.v);   /* square imaginary*/
    tv1.v = _mm512_add_ps(v1r.v, v1i.v);   /* Sum */
    _mm512_storeu_ps(&oArr[i], tv1.v);     /* Save out */
  }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* AVX Vector */
  /* Do blocks of 8 as vector */
  for (i=0; i<Fin_r->arraySize-8; i+=8) {
    v1r.v = _mm256_loadu_ps(&iArr_r[i]);   /* input real */
    v1r.v = _mm256_mul_ps(v1r.v, v1r.v);   /* square real*/
    v1i.v = _mm256_loadu_ps(&iArr_i[i]);   /* input imaginary */
    v1i.v = _mm256_mul_ps(v1i.v, v1i.v);   /* square imaginary*/
    tv1.v = _mm256_add_ps(v1i.v, v1r.v);   /* Sum */
    _mm256_storeu_ps(&oArr[i], tv1.v);     /* Save out */
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  /* Whatever is left */
  for (i=ilast; i<Fin_r->arraySize; i++) 
      oArr[i] = iArr_r[i]*iArr_r[i]+iArr_i[i]*iArr_i[i];
} /* end  ObitFArrayCplxAmp2 */

/**
 * Amplitude from a pair of reals 
 * out = sqrt(Fin_r**2 + Fin_i**2)
 * Blanking NOT supported
 * \param Fin_r  Real part of input
 * \param Fin_i  Imaginary part of input
 * \param out    Output object
 */
void ObitFArrayCplxAmp (ObitFArray* Fin_r, ObitFArray* Fin_i, 
			 ObitFArray* out)
{
  olong i, ilast=0;
  ofloat *iArr_r = Fin_r->array, *iArr_i = Fin_i->array;
  ofloat *oArr = out->array;
#if HAVE_AVX512==1
  CV16SF v1r, v1i, tv1;
#elif HAVE_AVX==1
  CV8SF  v1r, v1i, tv1;
#endif

   /* error checks */
  g_assert (ObitFArrayIsCompatable(Fin_r, out));
  g_assert (ObitFArrayIsCompatable(Fin_i, out));

#if HAVE_AVX512==1  /* AVX512 Vector */
  /* Do blocks of 16 as vector */
  for (i=0; i<Fin_r->arraySize-16; i+=16) {
    v1r.v = _mm512_loadu_ps(&iArr_r[i]);   /* input real */
    v1r.v = _mm512_mul_ps(v1r.v, v1r.v);   /* square real*/
    v1i.v = _mm512_loadu_ps(&iArr_i[i]);   /* input imaginary */
    v1i.v = _mm512_mul_ps(v1i.v, v1i.v);   /* square imaginary*/
    tv1.v = _mm512_add_ps(v1r.v, v1i.v);   /* Sum */
    tv1.v = _mm512_sqrt_ps(tv1.v);         /* Sqrt */
    _mm512_storeu_ps(&oArr[i], tv1.v);   /* Save out2 */
  }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* AVX Vector */
  /* Do blocks of 8 as vector */
  for (i=0; i<Fin_r->arraySize-8; i+=8) {
    v1r.v = _mm256_loadu_ps(&iArr_r[i]);   /* input real */
    v1r.v = _mm256_mul_ps(v1r.v, v1r.v);   /* square real*/
    v1i.v = _mm256_loadu_ps(&iArr_i[i]);   /* input imaginary */
    v1i.v = _mm256_mul_ps(v1i.v, v1i.v);   /* square imaginary*/
    tv1.v = _mm256_add_ps(v1i.v, v1r.v);   /* Sum */
    tv1.v = _mm256_sqrt_ps(tv1.v);         /* Sqrt */
    _mm256_storeu_ps(&oArr[i], tv1.v);     /* Save out2 */
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  /* Whatever is left */
  for (i=ilast; i<Fin_r->arraySize; i++) 
    oArr[i] = sqrt(iArr_r[i]*iArr_r[i]+iArr_i[i]*iArr_i[i]);
} /* end  ObitFArrayCplxAmp */

/**
 * Phase from a pair of reals 
 * out = atan2(Fin_i, Fin_r)
 * Blanking NOT supported
 * \param Fin_r  Real part of input
 * \param Fin_i  Imaginary part of input
 * \param out    Output object
 */
void ObitFArrayCplxPhase (ObitFArray* Fin_r, ObitFArray* Fin_i, 
			  ObitFArray* out)
{
  olong i, ilast=0;
  ofloat *iArr_r = Fin_r->array, *iArr_i = Fin_i->array;
  ofloat *oArr = out->array;
#if HAVE_AVX512X==1
  CV16SF v1r, v1i, tv1;
#elif HAVE_AVXX==1
  CV8SF  v1r, v1i, tv1;
#endif

   /* error checks */
  g_assert (ObitFArrayIsCompatable(Fin_r, out));
  g_assert (ObitFArrayIsCompatable(Fin_i, out));

  /****************************************************************
   atan2 is implemented in the Short Vector Math Library (SVML)
   which doesn't seem to be supported by gcc
   AVX disabled
  ****************************************************************/

#if HAVE_AVX512X==1  /* AVX512 Vector */
  /* Do blocks of 16 as vector */
  for (i=0; i<Fin_r->arraySize-16; i+=16) {
    v1r.v = _mm512_loadu_ps(&iArr_r[i]);   /* input real */
    v1i.v = _mm512_loadu_ps(&iArr_i[i]);   /* input imaginary */
    tv1.v = _mm512_atan2_ps(v1i.v, v1r.v); /* atan2 */
   _mm512_storeu_ps(&oArr[i], tv1.v.v);    /* Save out */
  }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVXX==1  /* AVX Vector */
  /* Do blocks of 8 as vector */
  for (i=0; i<Fin_r->arraySize-8; i+=8) {
    v1r.v = _mm256_loadu_ps(&iArr_r[i]);    /* input real */
    v1i.v = _mm256_loadu_ps(&iArr_i[i]);    /* input imaginary */
    tv1.v = _mm256_atan2_ps(v1i.v, v1r.v);  /* atan2 */
    _mm256_storeu_ps(&oArr[i], tv1.v);      /* Save out */
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  /* Whatever is left */
  for (i=ilast; i<Fin_r->arraySize; i++) 
    oArr[i] = atan2(iArr_i[i], iArr_r[i]);
} /* end  ObitFArrayCplxPhase */

/**
 *  Form complex from two FArrays, multiply by complex scalar and 
 *  complex accumulate
 * MAY BE BUGGY
 * Accum += Fin1 * cscalar
 *  Blanking NOT supported, replace with 0
 * \param Fin_r   Input FArray for real part, replace blanks w/ 0
 * \param Fin_i   Input FArray for imaginary part
 * \param cscalar Complex Scalar value [r,i]
 * \param Accum_r Real Accumulator
 * \param Accum_i Imaginary Accumulator
 */
void ObitFArrayCplxSMulAccum (ObitFArray* Fin_r, ObitFArray* Fin_i, 
			      ofloat cscalar[2], 
			      ObitFArray* Accum_r, ObitFArray* Accum_i) {
  olong nThreads;
  FAFuncArg **threadArgs;
 
  /* error checks */
  g_assert (ObitFArrayIsCompatable(Fin_r, Fin_i));
  g_assert (ObitFArrayIsCompatable(Accum_r, Fin_i));
  
  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (Fin_r->thread, Fin_r, Fin_i, NULL, Accum_r, Accum_i, 
			     NULL, NULL, 2*sizeof(ofloat), 0, 0, 0, 0, 0, 0,
			     &threadArgs);

  /* Call routine to do the work */
  ObitFArrayCplxSMulAccumTh (Fin_r, Fin_i, cscalar, Accum_r, Accum_i, nThreads, (gpointer)threadArgs);

  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
} /* end  ObitFArrayCplxSMulAccum */

/**
 *  Form complex from two FArrays, multiply by complex scalar and 
 *  complex accumulate - threading initialized elsewhere.
 * MAY BE BUGGY
 * Accum += Fin1 * cscalar
 *  Blanking NOT supported, replace with 0
 *  ObitFArrayMakeFAFuncArgs to initialize threading with (elem)
 * \param Fin_r   Input FArray for real part,      (in)
 * \param Fin_i   Input FArray for imaginary part  (in2)
 * \param cscalar Complex Scalar value [r,i]
 * \param Accum_r Real Accumulator                 (FA_3)
 * \param Accum_i Imaginary Accumulator            (FA_4)
 * \param nThreads   Number of thread args allocated
 * \param threadArgs Thread arguments
 */
void ObitFArrayCplxSMulAccumTh (ObitFArray* Fin_r, ObitFArray* Fin_i, 
				ofloat cscalar[2], 
				ObitFArray* Accum_r, ObitFArray* Accum_i,
				olong nThreads, gpointer inArgs) {
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread;
  ofloat *farr;
  FAFuncArg **threadArgs = (FAFuncArg**)inArgs;

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
    farr = threadArgs[i]->arg1; farr[0] = cscalar[0]; farr[1] = cscalar[1];
    /* Reset pointers */
    if (threadArgs[i]->in)   ObitFArrayUnref(threadArgs[i]->in);
    threadArgs[i]->in  = ObitFArrayRef(Fin_r);
    if (threadArgs[i]->in2)  ObitFArrayUnref(threadArgs[i]->in2);
    threadArgs[i]->in2 = ObitFArrayRef(Fin_i);
    if (threadArgs[i]->FA_3) ObitFArrayUnref(threadArgs[i]->FA_3);
    threadArgs[i]->FA_3 = ObitFArrayRef(Accum_r);
    if (threadArgs[i]->FA_4) ObitFArrayUnref(threadArgs[i]->FA_4);
    threadArgs[i]->FA_4 = ObitFArrayRef(Accum_i);
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  ObitThreadIterator (threadArgs[0]->thread, nTh, 
		      (ObitThreadFunc)ThreadFACplxSMulAccum,
		      (gpointer**)threadArgs);

} /* end  ObitFArrayCplxSMulAccumTh */

/**
 *  Multiply two pairs of FArrays and complex accumulate
 * MAY BE BUGGY
 * Accum += Fin1 * Fin2
 * Blanking NOT supported, replace with 0
 * \param Fin1_r   First FArray real part
 * \param Fin1_i   First FArray imaginary part
 * \param Fin2_r   Second FArray real part
 * \param Fin2_i   Second FArray imaginary part
 * \param Accum_r  Real Accumulator
 * \param Accum_r  Imaginary Accumulator
 */
void ObitFArrayCplxMulAccum (ObitFArray* Fin1_r, ObitFArray* Fin1_i, 
			     ObitFArray* Fin2_r, ObitFArray* Fin2_i, 
			     ObitFArray* Accum_r, ObitFArray* Accum_i) {
  olong nThreads;
  FAFuncArg **threadArgs;

  /* error checks */
  g_assert (ObitFArrayIsCompatable(Fin1_r, Fin2_r));
  g_assert (ObitFArrayIsCompatable(Accum_r, Fin1_r));
  
  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (Fin1_r->thread, Fin1_r, Fin1_i, NULL, Fin2_r, Fin2_i, Accum_r, Accum_i, 
			     0, 0, 0, 0, 0, 0, 0,
			     &threadArgs);

  /* Call routine to do the work */
  ObitFArrayCplxMulAccumTh (Fin1_r, Fin1_i, Fin2_r, Fin2_i, Accum_r,  Accum_i, 
			    nThreads, (gpointer)threadArgs);

  /* Free local objects */
  KillFAFuncArgs(nThreads, threadArgs);
} /* end  ObitFArrayCplxMulAccum */

/**
 *  Form complex from two FArrays, multiply by complex scalar and 
 *  complex accumulate - threading managed elsewhere.
 * MAY BE BUGGY
 * Accum += Fin1 * Fin2
 *  ObitFArrayMakeFAFuncArgs to initialize threading with (elem)
 * \param Fin1_r   First FArray real part         (in)
 * \param Fin1_i   First FArray imaginary part    (in2)
 * \param Fin2_r   Second FArray real part        (FA_3)
 * \param Fin2_i   Second FArray imaginary part   (FA_4)
 * \param Accum_r  Real Accumulator               (FA_5)
 * \param Accum_r  Imaginary Accumulator          (FA_6)
 * \param nThreads   Number of thread args allocated
 * \param threadArgs Thread arguments
 */
void ObitFArrayCplxMulAccumTh (ObitFArray* Fin1_r, ObitFArray* Fin1_i, 
			       ObitFArray* Fin2_r, ObitFArray* Fin2_i, 
			       ObitFArray* Accum_r,ObitFArray* Accum_i,
			       olong nThreads, gpointer inArgs) {
  olong i;
  olong nTh, nElem, loElem, hiElem, nElemPerThread;
  FAFuncArg **threadArgs = (FAFuncArg**)inArgs;

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
    if (threadArgs[i]->in)   ObitFArrayUnref(threadArgs[i]->in);
    threadArgs[i]->in  = ObitFArrayRef(Fin1_r);
    if (threadArgs[i]->in2)  ObitFArrayUnref(threadArgs[i]->in2);
    threadArgs[i]->in2 = ObitFArrayRef(Fin1_i);
    if (threadArgs[i]->FA_3) ObitFArrayUnref(threadArgs[i]->FA_3);
    threadArgs[i]->FA_3 = ObitFArrayRef(Fin2_r);
    if (threadArgs[i]->FA_4) ObitFArrayUnref(threadArgs[i]->FA_4);
    threadArgs[i]->FA_4 = ObitFArrayRef(Fin2_i);
    if (threadArgs[i]->FA_5) ObitFArrayUnref(threadArgs[i]->FA_5);
    threadArgs[i]->FA_5 = ObitFArrayRef(Accum_r);
    if (threadArgs[i]->FA_6) ObitFArrayUnref(threadArgs[i]->FA_6);
    threadArgs[i]->FA_6 = ObitFArrayRef(Accum_i);
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elem */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  ObitThreadIterator (threadArgs[0]->thread, nTh, 
		      (ObitThreadFunc)ThreadFACplxMulAccum,
		      (gpointer**)threadArgs);

} /* end  ObitFArrayCplxMulAccumTh */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitFArrayClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitFArrayClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitFArrayClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitFArrayClassInfoDefFn (gpointer inClass)
{
  ObitFArrayClassInfo *theClass = (ObitFArrayClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitFArrayClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitFArrayClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitFArrayGetClass;
  theClass->newObit       = (newObitFP)newObitFArray;
  theClass->ObitCopy      = (ObitCopyFP)ObitFArrayCopy;
  theClass->ObitFArraySubArr = 
    (ObitFArraySubArrFP)ObitFArraySubArr;
  theClass->ObitFArrayTranspose = 
    (ObitFArrayTransposeFP)ObitFArrayTranspose;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitFArrayClear;
  theClass->ObitInit      = (ObitInitFP)ObitFArrayInit;
  theClass->ObitFArrayCreate = (ObitFArrayCreateFP)ObitFArrayCreate;
  theClass->ObitFArrayIsCompatable = 
    (ObitFArrayIsCompatableFP)ObitFArrayIsCompatable;
  theClass->ObitFArrayIndex  = (ObitFArrayIndexFP)ObitFArrayIndex;
  theClass->ObitFArrayFill   = (ObitFArrayFillFP)ObitFArrayFill;
  theClass->ObitFArrayNeg    = (ObitFArrayNegFP)ObitFArrayNeg;
  theClass->ObitFArrayAbs    = (ObitFArrayAbsFP)ObitFArrayAbs;
  theClass->ObitFArraySin    = (ObitFArraySinFP)ObitFArraySin;
  theClass->ObitFArrayCos    = (ObitFArrayCosFP)ObitFArrayCos;
  theClass->ObitFArraySqrt   = (ObitFArraySqrtFP)ObitFArraySqrt;
  theClass->ObitFArrayMax    = (ObitFArrayMaxFP)ObitFArrayMax;
  theClass->ObitFArrayMaxAbs = (ObitFArrayMaxAbsFP)ObitFArrayMaxAbs;
  theClass->ObitFArrayMin    = (ObitFArrayMinFP)ObitFArrayMin;
  theClass->ObitFArrayDeblank= (ObitFArrayDeblankFP)ObitFArrayDeblank;
  theClass->ObitFArraySum    = (ObitFArraySumFP)ObitFArraySum;
  theClass->ObitFArrayCount  = (ObitFArrayCountFP)ObitFArrayCount;
  theClass->ObitFArraySAdd   = (ObitFArraySAddFP)ObitFArraySAdd;
  theClass->ObitFArraySMul   = (ObitFArraySMulFP)ObitFArraySMul;
  theClass->ObitFArraySDiv   = (ObitFArraySDivFP)ObitFArraySDiv;
  theClass->ObitFArrayClip   = (ObitFArrayClipFP)ObitFArrayClip;
  theClass->ObitFArrayInClip = (ObitFArrayInClipFP)ObitFArrayInClip;
  theClass->ObitFArrayBlank  = (ObitFArrayBlankFP)ObitFArrayBlank;
  theClass->ObitFArrayMaxArr = (ObitFArrayMaxArrFP)ObitFArrayMaxArr;
  theClass->ObitFArrayMinArr = (ObitFArrayMinArrFP)ObitFArrayMinArr;
  theClass->ObitFArrayExtArr = (ObitFArrayExtArrFP)ObitFArrayExtArr;
  theClass->ObitFArraySumArr = (ObitFArraySumArrFP)ObitFArraySumArr;
  theClass->ObitFArrayAvgArr = (ObitFArrayAvgArrFP)ObitFArrayAvgArr;
  theClass->ObitFArrayAdd    = (ObitFArrayAddFP)ObitFArrayAdd;
  theClass->ObitFArrayAddAbs = (ObitFArrayAddAbsFP)ObitFArrayAddAbs;
  theClass->ObitFArraySub    = (ObitFArraySubFP)ObitFArraySub;
  theClass->ObitFArraySign   = (ObitFArraySignFP)ObitFArraySign;
  theClass->ObitFArrayMul    = (ObitFArrayMulFP)ObitFArrayMul;
  theClass->ObitFArrayDiv    = (ObitFArrayDivFP)ObitFArrayDiv;
  theClass->ObitFArrayDivClip= (ObitFArrayDivClipFP)ObitFArrayDivClip;
  theClass->ObitFArrayDot    = (ObitFArrayDotFP)ObitFArrayDot;
  theClass->ObitFArrayRMS    = (ObitFArrayRMSFP)ObitFArrayRMS;
  theClass->ObitFArrayRawRMS = (ObitFArrayRawRMSFP)ObitFArrayRawRMS;
  theClass->ObitFArrayRMS0   = (ObitFArrayRMS0FP)ObitFArrayRMS0;
  theClass->ObitFArrayRMSQuant = (ObitFArrayRMSQuantFP)ObitFArrayRMSQuant;
  theClass->ObitFArrayQuant  = (ObitFArrayQuantFP)ObitFArrayQuant;
  theClass->ObitFArrayMulColRow = 
    (ObitFArrayMulColRowFP)ObitFArrayMulColRow;
  theClass->ObitFArray1DCenter = 
    (ObitFArray1DCenterFP)ObitFArray1DCenter;
  theClass->ObitFArray2DCenter = 
    (ObitFArray2DCenterFP)ObitFArray2DCenter;
  theClass->ObitFArray2DSymInv = 
    (ObitFArray2DSymInvFP)ObitFArray2DSymInv;
  theClass->ObitFArray2DCGauss = 
    (ObitFArray2DCGaussFP)ObitFArray2DCGauss;
  theClass->ObitFArrayShiftAdd = 
    (ObitFArrayShiftAddFP)ObitFArrayShiftAdd;
  theClass->ObitFArrayPad = (ObitFArrayPadFP)ObitFArrayPad;
  theClass->ObitFArrayConvGaus = 
    (ObitFArrayConvGausFP)ObitFArrayConvGaus;
  theClass->ObitFArraySelInc = 
    (ObitFArraySelIncFP)ObitFArraySelInc;
  theClass->ObitFArrayHisto = 
    (ObitFArrayHistoFP)ObitFArrayHisto;
  theClass->ObitFArrayExp = 
    (ObitFArrayExpFP)ObitFArrayExp;
  theClass->ObitFArrayLog = 
    (ObitFArrayLogFP)ObitFArrayLog;
  theClass->ObitFArrayPow = 
    (ObitFArrayPowFP)ObitFArrayPow;
  theClass->ObitFArrayRandom = 
    (ObitFArrayRandomFP)ObitFArrayRandom;
  theClass->ObitFArrayRandomFill = 
    (ObitFArrayRandomFillFP)ObitFArrayRandomFill;
  theClass->ObitFArrayRectFill = 
    (ObitFArrayRectFillFP)ObitFArrayRectFill;
  theClass->ObitFArrayRoundFill = 
    (ObitFArrayRoundFillFP)ObitFArrayRoundFill;
  theClass->ObitFArrayMaxSetValue = 
    (ObitFArrayMaxSetValueFP)ObitFArrayMaxSetValue;

} /* end ObitFArrayClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitFArrayInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFArray *in = inn;

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

} /* end ObitFArrayInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitFArray* cast to an Obit*.
 */
void ObitFArrayClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFArray *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  if (in->array)  {ObitMemFree(in->array);  in->array = NULL;}
  if (in->naxis)  {ObitMemFree(in->naxis); in->naxis = NULL;}
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitFArrayClear */

/**
 * Make arguments for a Threaded ThreadFAFunc?
 * For use outside the ObitFArray class
 * Use NULL or 0 for values not used
 * \param thread     ObitThread object to be used
 * \param in         FA to be operated on
 * \param in2        2nd FA to be operated on
 * \param out        output FA
 * \param FA_3       Additional FArray
 * \param FA_4       Additional FArray
 * \param FA_5       Additional FArray
 * \param FA_6       Additional FArray
 * \param larg1      Length of function dependent arg1 in floats
 * \param larg2      Length of function dependent arg2 in floats
 * \param larg3      Length of function dependent arg3 in floats
 * \param larg4      Length of function dependent arg4 in floats
 * \param larg5      Length of function dependent arg5 in floats
 * \param larg6      Length of function dependent arg6 in floats
 * \param larg7      Length of function dependent arg7 in floats
 * \param ThreadArgs[out] Created array of CAFuncArg, 
 *                   delete with ObitFArrayKillCAFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
olong ObitFArrayMakeFAFuncArgs (ObitThread *thread, 
				ObitFArray *in,	ObitFArray *in2, 
				ObitFArray *out,
				ObitFArray *FA_3, ObitFArray *FA_4,
				ObitFArray *FA_5, ObitFArray *FA_6,
				olong larg1, olong larg2, olong larg3, 
				olong larg4, olong larg5,
				olong larg6, olong larg7, 
				gpointer *outArgs)
{
  olong nThreads;
  FAFuncArg **ThreadArgs=NULL;
  /* Use class internal version */
  nThreads = MakeFAFuncArgs(thread, in, in2, out, FA_3, FA_4, FA_5, FA_6, 
			    larg1,larg2,larg3,larg4,larg5,larg6,larg7,
			    &ThreadArgs);
  
  *outArgs = (gpointer)ThreadArgs;
  return nThreads;
} /*  end ObitFArrayMakeFAFuncArgs */

/**
 * Delete arguments for ThreadfAFunc, stop thread pool
 * \param nargs      number of elements in ThreadArgs.
 * \param inArgs     pointer to array of FAFuncArg
 */
void ObitFArrayKillFAFuncArgs (olong nargs, gpointer inArgs)
{
  FAFuncArg **ThreadArgs = (FAFuncArg**)inArgs;

  if (ThreadArgs==NULL) return;

  /* Use class function */
  KillFAFuncArgs (nargs, ThreadArgs);
} /*  end ObitFArrayKillFAFuncArgs */

/**
 * Find maximum value and position in a subset of an FArray, 
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li value    Return Value, blanked if no data
 * \li pos      Return position
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAMax (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in        = largs->in;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;
  olong      *pos       = largs->pos;


  /* local */
  olong      i, temp, maxCell;
  ofloat     maxVal, fblank = ObitMagicF();

  if (hiElem<loElem) goto finish;

  /* Loop over array */
  largs->value = fblank;  /* In case no valid data*/
  maxCell = -1;
  maxVal  = -1.0E25;
  for (i=loElem; i<hiElem; i++) 
    {
      if ((in->array[i]!=fblank) && (in->array[i]>maxVal)) {
	maxCell = i;
	maxVal  = in->array[i];
      }
    }
  /* Return Value */
  largs->value = maxVal;

  /* Convert cell to index */
  temp = maxCell;
  for (i=0; i<in->ndim; i++) {
    pos[i] = temp % in->naxis[i];
    temp = (temp - pos[i]) / in->naxis[i];
  }

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadFAMax */

/**
 * Find minimum value and position in a subset of an FArray, 
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li value    Return Value, blanked if no data
 * \li pos      Return position
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAMin (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in        = largs->in;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;
  olong      *pos       = largs->pos;

  /* local */
  olong      i, temp, maxCell;
  ofloat     minVal, fblank = ObitMagicF();

  largs->value = fblank;  /* In case no valid data*/
  if (hiElem<loElem) goto finish;

  /* Loop over array */
  maxCell = -1;
  minVal = 1.0E25;
  for (i=loElem; i<hiElem; i++) 
    {
      if ((in->array[i]!=fblank) && (in->array[i]<minVal)) {
	maxCell = i;
	minVal = in->array[i];
      }
    }
  /* Return Value */
  largs->value = minVal;

  /* Convert cell to index */
  temp = maxCell;
  for (i=0; i<in->ndim; i++) {
    pos[i] = temp % in->naxis[i];
    temp = (temp - pos[i]) / in->naxis[i];
  }

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /*  end ThreadFAMin */

/**
 * Add portions of two FArrays, out = in + in2
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li in2      2nd ObitFArray to work on
 * \li out      Output ObitFArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAAdd (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in1       = largs->in;
  ObitFArray *in2       = largs->in2;
  ObitFArray *out       = largs->out;
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
  
} /*  end ThreadFAAdd */

/**
 * Subtract portions of two FArrays, out = in1 - in2
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li in2      2nd ObitFArray to work on
 * \li out      Output ObitFArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFASub (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in1       = largs->in;
  ObitFArray *in2       = largs->in2;
  ObitFArray *out       = largs->out;
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
  
} /*  end ThreadFASub */

/**
 * Multiply portions of two FArrays, out = in * in2
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li in2      2nd ObitFArray to work on
 * \li out      Output ObitFArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAMul (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in1       = largs->in;
  ObitFArray *in2       = largs->in2;
  ObitFArray *out       = largs->out;
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
    v2.v  = _mm512_mul_ps(v1.v, v2.v);             /* do operation */
    v2.v  = _mm512_mask_blend_ps(msk1,v2.v,vb.v);  /* replace blanks 1 */
    v2.v  = _mm512_mask_blend_ps(msk2,v2.v,vb.v);  /* replace blanks 2 */
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
    v2.v  = _mm256_mul_ps(v1.v, v2.v);             /* do operation */
    v2.v  = _mm256_blendv_ps(v2.v,vb.v, vm1.v);    /* replace blanks 1 */
    v2.v  = _mm256_blendv_ps(v2.v,vb.v, vm2.v);    /* replace blanks 2 */
    _mm256_storeu_ps(&out->array[i], v2.v);        /* Save */
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = loElem;/* Do all */
#endif
  for (i=ilast; i<hiElem; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = in1->array[i] * in2->array[i];
    else out->array[i] = fblank;
  }

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadFAMul */

/**
 * Divide portions of two FArrays, out = in / in2
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li in2      2nd ObitFArray to work on
 * \li out      Output ObitFArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFADiv (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in1       = largs->in;
  ObitFArray *in2       = largs->in2;
  ObitFArray *out       = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  /* local */
  olong   i;
  ofloat  fblank = ObitMagicF();

  if (hiElem<loElem) goto finish;

  /* Loop over array */
  for (i=loElem; i<hiElem; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)
	&& (in2->array[i]!=0.0)) 
      out->array[i] = in1->array[i] / in2->array[i];
    else out->array[i] = fblank;
  }

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadFADiv */

/**
 * Sum nonblanked elements of two arrays.
 *  out = (in1 + in2) or whichever is not blanked
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li in2      2nd ObitFArray to work on
 * \li out      Output ObitFArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFASumArr (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in1       = largs->in;
  ObitFArray *in2       = largs->in2;
  ObitFArray *out       = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  /* local */
  olong   i, ilast;
  ofloat  fblank = ObitMagicF();
#if HAVE_AVX512==1
  CV16SF v1, v2, vb, vzero;
  MASK16 msk1, msk2;
  ofloat tmp;
#elif HAVE_AVX==1
  CV8SF v1, v2, vb, vm1, vm2, vzero;
  ofloat tmp;
#endif

  if (hiElem<loElem) goto finish;

  /* Loop over array */
 #if HAVE_AVX512==1  /* Vector */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  tmp = 0.0;
  vzero.v  = _mm512_set1_ps(tmp);  /* vector of zeroes */
  /* Do blocks of 16 as vector */
  for (i=loElem; i<hiElem-16; i+=16) {
    v1.v  = _mm512_loadu_ps(&in1->array[i]);
    msk1  = _mm512_cmp_ps_mask(v1.v, vb.v, _CMP_EQ_OQ);   /* find blanks */
    v1.v  = _mm512_mask_blend_ps(msk1,v1.v,vzero.v);  /* replace blanks with zeroes */
    v2.v  = _mm512_loadu_ps(&in2->array[i]);
    msk2  = _mm512_cmp_ps_mask(v2.v, vb.v, _CMP_EQ_OQ);  /* find blanks */
    v2.v  = _mm512_mask_blend_ps(msk2,v2.v,vzero.v);     /* replace blanks with zeroes */
    v2.v  = _mm512_add_ps(v1.v, v2.v);                   /* do operation */
    msk2  = _mm512_kand(msk1, msk2);                     /* both blanked? */
    v2.v  = _mm512_mask_blend_ps(msk2,v2.v,vb.v);        /* replace blanks 1 */
    _mm512_storeu_ps(&out->array[i], v2.v);        /* Save */
  }
  ilast = i;  /* How far did I get? */
 #elif HAVE_AVX==1  /* Vector */
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  tmp = 0.0;
  vzero.v  = _mm256_broadcast_ss(&tmp);  /* vector of zeroes */
  /* Do blocks of 8 as vector */
  for (i=loElem; i<hiElem-8; i+=8) {
    v1.v  = _mm256_loadu_ps(&in1->array[i]);
    vm1.v = _mm256_cmp_ps(v1.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v1.v  = _mm256_blendv_ps(v1.v,vzero.v,vm1.v);  /* replace blanks with zeroes */
    v2.v  = _mm256_loadu_ps(&in2->array[i]);
    vm2.v = _mm256_cmp_ps(v2.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v2.v  = _mm256_blendv_ps(v2.v,vzero.v,vm2.v);  /* replace blanks with zeroes */
    v2.v  = _mm256_add_ps(v1.v, v2.v);             /* do operation */
    vm2.v = _mm256_and_ps(vm1.v, vm2.v);           /* both blanked */
    v2.v  = _mm256_blendv_ps(v2.v,vb.v, vm2.v);    /* replace blanks 1 */
    _mm256_storeu_ps(&out->array[i], v2.v);        /* Save */
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = loElem;/* Do all */
#endif
  for (i=ilast; i<hiElem; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = (in1->array[i] + in2->array[i]);
    else if (in1->array[i]==fblank)  /* 1 blanked */
      out->array[i] = in2->array[i];
    else if (in2->array[i]==fblank)
      out->array[i] = in1->array[i];  /* 2 blanked */
    else out->array[i] = fblank;      /* both blanked */
  } /* End loop over selected elements */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadFASumArr */

/**
 * Average nonblanked elements of two arrays.
 *  out = (in1 + in2)/2 or whichever is not blanked
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li in2      2nd ObitFArray to work on
 * \li out      Output ObitFArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAAvgArr (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in1       = largs->in;
  ObitFArray *in2       = largs->in2;
  ObitFArray *out       = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;


  /* local */
  olong   i;
  ofloat  fblank = ObitMagicF();

  if (hiElem<loElem) goto finish;

  /* Loop over array */
  for (i=loElem; i<hiElem; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = 0.5 * (in1->array[i] + in2->array[i]);
    else if (in1->array[i]==fblank)  /* 1 blanked */
      out->array[i] = in2->array[i];
    else if (in2->array[i]==fblank)
      out->array[i] = in1->array[i];  /* 2 blanked */
    else out->array[i] = fblank;      /* both blanked */

  } /* End loop over selected elements */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadFAAvgArr */

/**
 * Pick the larger nonblanked elements of two arrays.
 *  out = MAX (in1, in2) or whichever is not blanked
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li in2      2nd ObitFArray to work on
 * \li out      Output ObitFArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAMaxArr (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in1       = largs->in;
  ObitFArray *in2       = largs->in2;
  ObitFArray *out       = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;


  /* local */
  olong   i;
  ofloat  fblank = ObitMagicF();

  if (hiElem<loElem) goto finish;

  /* Loop over array */
  for (i=loElem; i<hiElem; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = MAX (in1->array[i], in2->array[i]);
    else if (in1->array[i]==fblank)  /* 1 blanked */
      out->array[i] = in2->array[i];
    else if (in2->array[i]==fblank)
      out->array[i] = in1->array[i];  /* 2 blanked */
    else out->array[i] = fblank;      /* both blanked */
  } /* End loop over selected elements */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  return NULL;
  
} /*  end ThreadFAMaxArr */

/**
 * Pick the lesser nonblanked elements of two arrays.
 *  out = MIN (in1, in2) or whichever is not blanked
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li in2      2nd ObitFArray to work on
 * \li out      Output ObitFArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAMinArr (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in1       = largs->in;
  ObitFArray *in2       = largs->in2;
  ObitFArray *out       = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  /* local */
  olong   i;
  ofloat  fblank = ObitMagicF();

  if (hiElem<loElem) goto finish;

  /* Loop over array */
  for (i=loElem; i<hiElem; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = MIN (in1->array[i], in2->array[i]);
    else if (in1->array[i]==fblank)  /* 1 blanked */
      out->array[i] = in2->array[i];
    else if (in2->array[i]==fblank)
      out->array[i] = in1->array[i];  /* 2 blanked */
    else out->array[i] = fblank;      /* both blanked */
  } /* End loop over selected elements */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
     ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadFAMinArr */

/**
 * Pick the more extreme (furthest from zero) nonblanked elements of two arrays.
 *  out = Extreme (in1, in2) or whichever is not blanked
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li in2      2nd ObitFArray to work on
 * \li out      Output ObitFArray
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAExtArr (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in1       = largs->in;
  ObitFArray *in2       = largs->in2;
  ObitFArray *out       = largs->out;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  /* local */
  olong   i;
  ofloat  fblank = ObitMagicF();

  if (hiElem<loElem) goto finish;

  /* Loop over array */
  for (i=loElem; i<hiElem; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) {
      if (fabs(in1->array[i])>fabs(in2->array[i]))
	out->array[i] = in1->array[i];
      else
	out->array[i] = in2->array[i];
    } else if (in2->array[i]!=fblank)  /* 1 blanked */
      out->array[i] = in2->array[i];
    else if (in1->array[i]!=fblank)
      out->array[i] = in1->array[i];  /* 2 blanked */
    else out->array[i] = fblank;      /* both blanked */
  } /* End loop over selected elements */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  return NULL;
  
} /*  end ThreadFAExtArr */

/**
 * Inner shiftAdd
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li first    First row (0-rel) number
 * \li last     Highest row(0-rel) number
 * \li value    scaling value
 * \li pos      [0] = plane number
 *              [1,2] = pos1, [3,4] = pos2
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAShAdd (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in1       = largs->in;
  ObitFArray *in2       = largs->in2;
  ObitFArray *out       = largs->out;
  ofloat     scalar     = largs->value;
  olong      loRow      = largs->first;
  olong      hiRow      = largs->last;
  olong      *pos       = largs->pos;


  /* local */
  olong   ix, iy, lox, hix, ip, indx1, indx2, nx1, nx2, ny1, ny2, offx, offy, lenp1, lenp2;
  olong   pos1[2], pos2[2];
  gboolean areSame=FALSE;
  ofloat  fblank = ObitMagicF();

  if (hiRow<loRow) goto finish;

  /* IS the output one of the inputs? */
  areSame = (in1==out) || (in2==out);

  ip = pos[0];                          /* Plane number */
  pos1[0] = pos[1]; pos1[1] = pos[2];   /* Alignment pixels */
  pos2[0] = pos[3]; pos2[1] = pos[4];

  /* Size of in1/out */
  nx1 = in1->naxis[0];
  ny1 = in1->naxis[1];
  lenp1 = nx1*ny1;
 
  /* Size of in2 */
  nx2 = in2->naxis[0];
  ny2 = in2->naxis[1];
  if (in2->ndim>2) lenp2 = nx2*ny2;
  else lenp2 = 0;

  /* determine regions of overlap in in1/out of in2 */
  offx = pos2[0] - pos1[0];
  offy = pos2[1] - pos1[1];
  lox = MAX (0, pos1[0] - in2->naxis[0]/2);
  hix = MIN (nx1-1, pos1[0] + in2->naxis[0]/2);
  /* In case in2 not centered */
  if (lox+offx<0) lox -= lox+offx;
  if (hix+offx>=nx2) hix -= hix+offx-nx2+1;
  /* Keep in range */
  lox = MAX (0, lox);
  hix = MIN (hix, nx1-1);

  /* Loop over rows */
  for (iy=loRow; iy<=hiRow; iy++) {
    
    /* Loop over columns */
    for (ix=lox; ix<=hix; ix++) {
      
      /* indices in arrays */
      indx1 = ip*lenp1 + iy*nx1 + ix;
      indx2 = ip*lenp2 + (iy+offy) * nx2 + ix + offx;
      
      /* do operation */
      if (in1->array[indx1]!=fblank) {
	if (in2->array[indx2]!=fblank) {
	  out->array[indx1] = in1->array[indx1] + scalar * in2->array[indx2];
	}
      } else if (!areSame) {  /* Don't blank if the output is the accumulation of an input */
	out->array[indx1] = fblank;
      }
    } /* end x loop */
  } /* end y loop */ 

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadFADShAdd */

/**
 * Find maximum abs value and position in a subset of an FArray, 
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li value    Return Value, blanked if no data
 * \li pos      Return position
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAAbsMax (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in        = largs->in;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;
  olong      *pos       = largs->pos;

  /* local */
  olong  i, temp, maxCell;
  ofloat maxAVal=-1.0, maxVal, *data, val, fblank = ObitMagicF();

  largs->value = fblank;  /* In case no valid data*/
  if (hiElem<loElem) goto finish;

  /* Loop over array */
  maxCell = -1;
  maxVal  = 0.0;
  data    = in->array;
  for (i=loElem; i<hiElem; i++) 
    {
       val = data[i];
      if ((val!=fblank) && (fabs(val)>maxAVal)) {
	maxCell = i;
	maxAVal = fabs(val);
	maxVal  = val;
      }
    }
  /* Return Value */
  largs->value = maxVal;

  /* Convert cell to index */
  temp = maxCell;
  for (i=0; i<in->ndim; i++) {
    pos[i] = temp % in->naxis[i];
    temp = (temp - pos[i]) / in->naxis[i];
  }

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /*  end ThreadFAAbsMax */

/**
 * Determine RMS sums in a subset of an FArray, 
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li arg1     (olong) count of valid elements
 * \li arg2     (ofloat) sum of elements
 * \li arg3     (ofloat) sum of elements^2
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFARMSSum (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in        = largs->in;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  /* local */
  olong i, ilast;
  ofloat sum, sum2, fblank = ObitMagicF();
  ollong count;
#if HAVE_AVX512==1
  CV16SF v, vb, vzero, vone, vcnt, vsum, vsum2, vt;
  MASK16 msk;
  ofloat tmp;
#elif HAVE_AVX==1
  CV8SF v, vb, vm, vzero,  vsum, vsum2;
  /*CV8SF v, vb, vm, vzero, vone, vcnt, vsum, vsum2, vt;*/
  ofloat tmp;
  ollong cnt, j;
#endif

  if (hiElem<loElem) goto finish;

  /* Loop over array by vector type */
  count = 0; sum = sum2 = 0.0; 
#if HAVE_AVX512==1  /* Vector */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  tmp = 0.0;
  vzero.v  = _mm512_set1_ps(tmp);  /* vector of zeroes */
  vcnt.v   = _mm512_set1_ps(tmp);  /* initialize counts */
  vsum.v   = _mm512_set1_ps(tmp);  /* initialize sum */
  vsum2.v  = _mm512_set1_ps(tmp);  /* initialize sum**2 */
  tmp = 1.0;
  vone.v   = _mm512_set1_ps(tmp);  /* vector of ones */
  /* Do blocks of 16 as vector */
  for (i=loElem; i<hiElem-16; i+=16) {
    v.v  = _mm512_loadu_ps(&in->array[i]);
    msk  = _mm512_cmp_ps_mask(v.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v.v  = _mm512_mask_blend_ps(msk,v.v,vzero.v);     /* replace blanks with zeroes */
    vt.v = _mm512_mask_blend_ps(msk,vone.v, vzero.v); /* 1 for good, 0 blank */
    vcnt.v = _mm512_add_ps(vcnt.v,vt.v);  /* sum counts */
    vsum.v = _mm512_add_ps(vsum.v,v.v);   /* sum data */
    v.v    = _mm512_mul_ps(v.v, v.v);     /* square */
    vsum2.v = _mm512_add_ps(vsum2.v,v.v); /* sum data squared */
  }
  ilast = i;  /* How far did I get? */
  vcnt.v = _mm512_roundscale_ps(vcnt.v,32);  /* round */
  for (i=0; i<16; i++) {
    count += (olong)vcnt.f[i];
    sum   += vsum.f[i];
    sum2  += vsum2.f[i];
  }
#elif HAVE_AVX==1  /* Vector */
  cnt = 0;
  vb.v   = _mm256_broadcast_ss(&fblank);  /* vector of blanks */
  tmp = 0.0;
  vzero.v  = _mm256_broadcast_ss(&tmp);  /* vector of zeroes */
  /*vcnt.v   = _mm256_broadcast_ss(&tmp);   initialize counts */
  vsum.v   = _mm256_broadcast_ss(&tmp);  /* initialize sum */
  vsum2.v  = _mm256_broadcast_ss(&tmp);  /* initialize sum**2 */
  /*tmp = 1.0;
    vone.v   = _mm256_broadcast_ss(&tmp);   vector of ones */
  /* Do blocks of 8 as vector */
  for (i=loElem; i<hiElem-8; i+=8) {
    /* don't use AVX for counts */
    for (j=i; j<i+8; j++) if (in->array[j]!=fblank) cnt++;
    v.v  = _mm256_loadu_ps(&in->array[i]);
    vm.v = _mm256_cmp_ps(v.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v.v  = _mm256_blendv_ps(v.v,vzero.v,vm.v);   /* replace blanks with zeroes */
    /*vt.v = _mm256_blendv_ps(vone.v, vzero.v, vm.v);   1 for good, 0 blank */
    /* vcnt.v = _mm256_add_ps(vcnt.v,vt.v);   sum counts - doesn't work so well */
    vsum.v = _mm256_add_ps(vsum.v,v.v);   /* sum data */
    v.v    = _mm256_mul_ps(v.v, v.v);     /* square */
    vsum2.v = _mm256_add_ps(vsum2.v,v.v); /* sum data squared */
  }
  ilast = i;  /* How far did I get? */
  /*vcnt.v = _mm256_round_ps(vcnt.v,_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC);   round */
  count = cnt;
  for (i=0; i<8; i++) {
    /*count += (olong)vcnt.f[i];*/
    sum   += vsum.f[i];
    sum2  += vsum2.f[i];
  }
#else /* Scalar */
  ilast = loElem;  /* Do all */
#endif
  for (i=ilast; i<hiElem; i++) 
    {
      if (in->array[i]!=fblank) {
	count++;
	sum  += in->array[i];
	sum2 += in->array[i] * in->array[i];
      }
    }

  /* Return values */
  *(ollong*)largs->arg1 = count; 
  *(ofloat*)largs->arg2 = sum;
  *(ofloat*)largs->arg3 = sum2;

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /*  end ThreadFARMSSum */

/**
 * Accumulate histogram elements 
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li first    (olong) First element (1-rel) number
 * \li last     (olong) Highest element (1-rel) number
 * \li arg1     (olong) Number of cells in histogram
 * \li arg2     (ollong)  Count of pixels
 * \li arg3     (ofloat*) Histogram
 * \li arg4     (ofloat)  Max in histogram
 * \li arg5     (ofloat)  Min in histogram
 * \li arg6     (ollong)  Number of underflows
 * \li arg7     (ollong)  Number of overflows
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAHisto (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in        = largs->in;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;
  olong      numCell    = *(olong*)largs->arg1;
  ollong     count      = *(ollong*)largs->arg2;
  ofloat     *histo     = (ofloat*)largs->arg3;
  ofloat     amax       = *(ofloat*)largs->arg4;
  ofloat     amin       = *(ofloat*)largs->arg5;
  ollong     under      = *(ollong*)largs->arg6;
  ollong     over       = *(ollong*)largs->arg7;

  /* local */
  olong  i, icell, ilast;
  ofloat cellFact, fblank = ObitMagicF();
#if HAVE_AVX512XXX==1
  CV16SF v, vb, vflag, vmin, vfac;
  IV16SF vcell;
  MASK16 msk;
  ofloat tmp;
  olong j;
#elif HAVE_AVX==1
  CV8SF v, vb, vm, vflag, vmin, vfac;
  IV8SF vcell;
  ofloat tmp;
  olong j;
#endif

  if (hiElem<loElem) goto finish;

  /* Loop over array */
  cellFact =  numCell / (amax - amin + 1.0e-20);
  count = 0; under = 0; over = 0;
  for (i=0; i<numCell; i++) histo[i] = 0.0;
#if HAVE_AVX512XXX==1  /* 16 Vector - not sure it works right */
  vb.v   = _mm512_set1_ps(fblank);  /* vector of blanks */
  tmp = -cellFact*amin;
  vflag.v  = _mm512_set1_ps(tmp);  /* flagged values */
  tmp = amin;
  vmin.v  = _mm512_set1_ps(tmp);  /* amin */
  tmp = cellFact;
  vfac.v = _mm512_set1_ps(tmp);  /* fact */
  /* Do blocks of 16 as vector */
  for (i=loElem; i<hiElem-16; i+=16) {
    v.v    = _mm512_loadu_ps(&in->array[i]);            /* Data */
    msk    = _mm512_cmp_ps_mask(v.v, vb.v, _CMP_EQ_OQ); /* find blanks */
    v.v    = _mm512_mask_blend_ps(msk,v.v,vflag.v);     /* replace blanks */
    v.v    = _mm512_sub_ps(v.v,vmin.v);        /* -amin */
    v.v    = _mm512_mul_ps(v.v,vfac.v);        /* * fact */
    v.v    = _mm512_roundscale_ps(v.v,32);     /* round cell number */
    vcell.v = _mm512_cvtps_epi32(v.v);
    for (j=0; j<16; j++) {
      if (in->array[i+j] != fblank) {
	icell = vcell.i[j];
	if (icell<0) under++;
	else if (icell>=numCell) over++;
	else histo[icell]++;
	count++;
      }
    }
  }
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* Vector */
  vb.v   = _mm256_broadcast_ss(&fblank); /* vector of blanks */
  tmp = -cellFact*amin;
  vflag.v  = _mm256_broadcast_ss(&tmp);  /* flagged values */
  tmp = amin;
  vmin.v  = _mm256_broadcast_ss(&tmp);   /* amin */
  tmp = cellFact;
  vfac.v = _mm256_broadcast_ss(&tmp);    /* cell factor */
  /* Do blocks of 8 as vector */
  for (i=loElem; i<hiElem-8; i+=8) {
    v.v    = _mm256_loadu_ps(&in->array[i]);        /* Data */
    vm.v   = _mm256_cmp_ps(v.v, vb.v, _CMP_EQ_OQ);  /* find blanks */
    v.v    = _mm256_blendv_ps(v.v,vflag.v,vm.v);    /* replace blanks */
    v.v    = _mm256_sub_ps(v.v,vmin.v);        /* -amin */
    v.v    = _mm256_mul_ps(v.v,vfac.v);        /* * cell factor */
    v.v    = _mm256_round_ps(v.v,_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC);  /* round cell */
    vcell.v = _mm256_cvtps_epi32(v.v);
    for (j=0; j<8; j++) {
      if (in->array[i+j] != fblank) {
	icell = vcell.i[j];
	if (icell<0) under++;
	else if (icell>=numCell) over++;
	else histo[icell]++;
	count++;
      }
    }
  }
  ilast = i;  /* How far did I get? */
#else /* Scalar and rest */
  ilast = loElem;
#endif
  for (i=ilast; i<hiElem; i++) {
    if (in->array[i]!=fblank){
      icell = (olong)(0.49 + (cellFact * (in->array[i]-amin)));
      if (icell<0) under++;
      else if (icell>=numCell) over++;
      else histo[icell]++;
      count++;
    }
  }

  /* Save */
  *(ollong*)largs->arg2 = count;
  *(ollong*)largs->arg6 = under;
  *(ollong*)largs->arg7 = over;

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /*  end ThreadFAHisto */

/**
 * Thread convolve a list of points, 
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       ObitFArray to work on
 * \li first    First row (1-rel) row in in
 * \li last     Highest row (1-rel) in in
 * \li arg1     (ObitFArray*)list List of positions and fluxes of the Gaussians
 *              (x pixel, y pixel, flux)
 * \li arg2     (olong*)ncomp   Number of components in list  
 *              generally less than size of FArray).
 * \li arg3     (ofloat*) gauss[3]  Gaussian coefficients for (d_x*d_x, d_y*d_y, d_x*d_y)
 *                Gaussian maj = major axis FWHM, min=minor, pa = posn. angle
 *                cr=cos(pa+rotation), sr=sin(pa+rotation),
 *                cell_x, cell_y x, y cell spacing is same units as maj, min
 *                [0] = {(cr/min)^2 + ((sr/maj)^2)}*(cell_x^2)*4*log(2)
 *                [1] = {(sr/min)^2 + ((cr/maj)^2)}*(cell_y^2)*4*log(2)
 *                [2] = {(1/min)^2 - ((1/maj)^2)}*sr*cr*abs(cell_x*cell_y)*8*log(2)
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAConvGaus (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in        = largs->in;
  ObitFArray *list      = (ObitFArray*)largs->arg1;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last-1;
  olong      ncomp      = *(olong*)largs->arg2;
  ofloat     gauss[3];

  /* local */
  olong icomp, lrec, ix, iy, indx;
  ofloat *table, *image, dx, dy, farg, aa, bb, cc;
  ofloat minGaus=12.0;

  if (hiElem<loElem) goto finish;

  /* Get Gauss */
  memmove(gauss, largs->arg3, 3*sizeof(ofloat));

  /* Setup list to access */
  table = list->array;
  lrec  = list->naxis[0];

  /* access to in as image */
  image = in->array;

  /* Gaussian parameters */
  aa = gauss[0];
  bb = gauss[1];
  cc = gauss[2];

  /* Loop over elements in list */
  for (icomp=0; icomp<ncomp; icomp++) {
    if (table[2]==0.0) continue;  /* ignore zero flux */
    indx = loElem*in->naxis[0];   /* image array index */

    /* Loop over array convolving */
    for (iy = loElem; iy<=hiElem; iy++) {
      dy = iy - table[1];   /* y offset */
      if (bb*dy*dy>minGaus) {indx+=in->naxis[0]; continue;}
      for (ix = 0; ix<in->naxis[0]; ix++) {
	dx = ix - table[0];   /* x offset */
	if (aa*dx*dx>minGaus) {indx++; continue;}
	farg = -(aa*dx*dx + bb*dy*dy + cc*dx*dy);
	if (farg>-minGaus) {
	  image[indx] += table[2] * ObitExpCalc(farg);
	}
	indx++;
      }
    }
    table += lrec;
  } /* end loop over components list */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadFAConvGaus */

/**
 * Multiply complex from two FArrays by complex scalar and 
 * complex accumulate. Magic value blanking not supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       Real part of input
 * \li in2      Imaginary part of input
 * \li arg1     Real/Imaginary parts of scalar
 * \li FA_3     Output real accumulator
 * \li FA_4     Output imaginary accumulator
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFACplxSMulAccum (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *Fin_r     = largs->in;
  ObitFArray *Fin_i     = largs->in2;
  ofloat     *cscalar   = largs->arg1;
  ObitFArray *Accum_r   = largs->FA_3;
  ObitFArray *Accum_i   = largs->FA_4;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  ofloat tr1, ti1, tr2, ti2;
  olong i, ilast;
  ofloat *iArr_r = Fin_r->array,   *iArr_i = Fin_i->array;
  ofloat *oArr_r = Accum_r->array, *oArr_i = Accum_i->array;
#if HAVE_AVX512==1
  CV16SF v1r, v1i, sr, si, tv1, tv2, tv3;
#elif HAVE_AVXX==1
  CV8SF  v1r, v1i, sr, si, tv1, tv2, tv3;
#endif
  
  if (hiElem<loElem) goto finish;
  tr2 = cscalar[0]; ti2 = cscalar[1];
  
#if HAVE_AVX512==1  /* AVX 512 Vector (16 float) */
  sr.v = _mm512_set1_ps(tr2); /* Scalar real */
  si.v = _mm512_set1_ps(ti2); /* Scalar imaginary */
  for (i=loElem; i<hiElem-16; i+=16) { 
    v1r.v = _mm512_loadu_ps(&iArr_r[i]);  /* Input reals */
    tv1.v = _mm512_mul_ps (v1r.v, sr.v);  /* tr1*tr2 */
    v1i.v = _mm512_loadu_ps(&iArr_i[i]);  /* Input imaginaries */
    tv2.v = _mm512_mul_ps (v1i.v, si.v);  /* ti1*ti2 */
    tv1.v = _mm512_sub_ps (tv1.v, tv2.v); /* now Real part */
    tv3.v = _mm512_loadu_ps(&oArr_r[i]);  /* Output reals */
    tv3.v = _mm512_add_ps (tv3.v, tv1.v); /* Update reals */
    _mm512_storeu_ps(&oArr_r[i], tv3.v);  /* store reals */
    tv3.v = _mm512_mul_ps (v1i.v, sr.v);  /* ti1*tr2 */
    tv1.v = _mm512_mul_ps (v1r.v, si.v);  /* tr1*ti2 */
    tv2.v = _mm512_add_ps (tv3.v, tv1.v); /* now Imaginary part */
    tv3.v = _mm512_loadu_ps(&oArr_i[i]);  /* Output reals */
    tv3.v = _mm512_add_ps (tv3.v, tv2.v); /* Update imag */
    _mm512_storeu_ps(&oArr_i[i], tv3.v);  /* store reals */
  } /* end outer loop */
  ilast = i;  /* How far did I get? */
#elif HAVE_AVXXX==1  /* AVX Vector (8 float) */ 
  sr.v = _mm256_broadcast_ss(&tr2); /* Scalar real */
  si.v = _mm256_broadcast_ss(&ti2); /* Scalar imaginary */
  for (i=loElem; i<hiElem-8; i+=8) {
    v1r.v = _mm256_loadu_ps(&iArr_r[i]);  /* Input reals */
    tv1.v = _mm256_mul_ps (v1r.v, sr.v);  /* tr1*tr2 */
    v1i.v = _mm256_loadu_ps(&iArr_i[i]);  /* Input imaginaries */
    tv2.v = _mm256_mul_ps (v1i.v, si.v);  /* ti1*ti2 */
    tv1.v = _mm256_sub_ps (tv1.v, tv2.v); /* now Real part */
    tv3.v = _mm256_loadu_ps(&oArr_r[i]);  /* Output real */
    tv3.v = _mm256_add_ps (tv3.v, tv1.v); /* Update reals */
    _mm256_storeu_ps(&oArr_r[i], tv3.v);  /* Save real */
    tv3.v = _mm256_mul_ps (v1i.v, sr.v);  /* ti1*tr2 */
    tv1.v = _mm256_mul_ps (v1r.v, si.v);  /* tr1*ti2 */
    tv2.v = _mm256_add_ps (tv3.v, tv1.v); /* now Imaginary part */
    tv3.v = _mm256_loadu_ps(&oArr_i[i]);  /* Output imaginary */
    tv3.v = _mm256_add_ps (tv3.v, tv2.v); /* Update imag */
    _mm256_storeu_ps(&oArr_i[i], tv3.v);  /* Save imag */
  } /* end outer loop */
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  /* Loop over whatever is left over */
  for (i=ilast; i<hiElem; i++) {
    tr1 = iArr_r[i];
    ti1 = iArr_i[i];
    oArr_r[i] += tr1*tr2 - ti1*ti2;
    oArr_i[i] += ti1*tr2 + tr1*ti2;
  } /* end loop over array */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadFACplxSMulAccum */

/**
 * Complex multiply from pairs of floats, accumulate
 * complex accumulate. Magic value blanking not supported.
 * Callable as thread
 * \param arg Pointer to FAFuncArg argument with elements:
 * \li in       Real part of 1st input
 * \li in2      Imaginary part of 1st input
 * \li FA_3     Real part of 2nd input
 * \li FA_4     Imaginary part of 2nd input
 * \li FA_5     Output real accumulator
 * \li FA_6     Output imaginary accumulator
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFACplxMulAccum (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *Fin1_r    = largs->in;
  ObitFArray *Fin1_i    = largs->in2;
  ObitFArray *Fin2_r    = largs->FA_3;
  ObitFArray *Fin2_i    = largs->FA_4;
  ObitFArray *Accum_r   = largs->FA_5;
  ObitFArray *Accum_i   = largs->FA_6;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last;

  ofloat tr1, ti1, tr2, ti2;
  olong i, ilast;
  ofloat *iArr1_r = Fin1_r->array,  *iArr1_i = Fin1_i->array;
  ofloat *iArr2_r = Fin2_r->array,  *iArr2_i = Fin2_i->array;
  ofloat *oArr_r  = Accum_r->array, *oArr_i = Accum_i->array;
#if HAVE_AVX512==1
  CV16SF v1r, v1i, v2r, v2i, tv1, tv2, tv3;
#elif HAVE_AVX==1
  CV8SF  v1r, v1i, v2r, v2i, tv1, tv2, tv3;
#endif
  
  if (hiElem<loElem) goto finish;
  
#if HAVE_AVX512==1  /* AVX 512 Vector (16 float) */
  for (i=loElem; i<hiElem-16; i+=16) {
    v1r.v = _mm512_loadu_ps(&iArr1_r[i]);  /* Input 1st  reals */
    v2r.v = _mm512_loadu_ps(&iArr2_r[i]);  /* Input 2nd reals */
    tv1.v = _mm512_mul_ps (v1r.v, v2r.v);  /* tr1*tr2 */
    v1i.v = _mm512_loadu_ps(&iArr1_i[i]);  /* Input 1st imaginaries */
    v2i.v = _mm512_loadu_ps(&iArr2_i[i]);  /* Input 2nd imaginaries */
    tv2.v = _mm512_mul_ps (v1i.v, v2i.v);  /* ti1*ti2 */
    tv1.v = _mm512_sub_ps (tv1.v, tv2.v);  /* now Real part */
    tv3.v = _mm512_loadu_ps(&oArr_r[i]);   /* Output reals */
    tv3.v = _mm512_add_ps (tv3.v, tv1.v);  /* Update reals */
    _mm512_storeu_ps(&oArr_r[i], tv3.v);   /* store reals */
    tv3.v = _mm512_mul_ps (v1i.v, v2r.v);  /* ti1*tr2 */
    tv1.v = _mm512_mul_ps (v1r.v, v2i.v);  /* tr1*ti2 */
    tv2.v = _mm512_add_ps (tv3.v, tv1.v);  /* now Imaginary part */
    tv3.v = _mm512_loadu_ps(&oArr_i[i]);   /* Output reals */
    tv3.v = _mm512_add_ps (tv3.v, tv2.v);  /* Update imag */
    _mm512_storeu_ps(&oArr_i[i], tv3.v);   /* store imag */
  } /* end outer loop */
  ilast = i;  /* How far did I get? */
#elif HAVE_AVX==1  /* AVX Vector (8 float) */ 
  for (i=loElem; i<hiElem-8; i+=8) {
    v1r.v = _mm256_loadu_ps(&iArr1_r[i]);  /* Input 1st reals */
    v2r.v = _mm256_loadu_ps(&iArr2_r[i]);  /* Input 2nd reals */
    tv1.v = _mm256_mul_ps (v1r.v ,v2r.v);  /* tr1*tr2 */
    v1i.v = _mm256_loadu_ps(&iArr1_i[i]);  /* Input 1st imaginaries */
    v2i.v = _mm256_loadu_ps(&iArr2_i[i]);  /* Input 2nd imaginaries */
    tv2.v = _mm256_mul_ps (v1i.v, v2i.v);  /* ti1*ti2 */
    tv1.v = _mm256_sub_ps (tv1.v, tv2.v);  /* now Real part */
    tv3.v = _mm256_loadu_ps(&oArr_r[i]);   /* Output reals */
    tv3.v = _mm256_add_ps (tv3.v, tv1.v);  /* Update reals */
    _mm256_storeu_ps(&oArr_r[i], tv3.v);   /* Save real */
    tv3.v = _mm256_mul_ps (v1i.v, v2r.v);  /* ti1*tr2 */
    tv1.v = _mm256_mul_ps (v1r.v, v2i.v);  /* tr1*ti2 */
    tv2.v = _mm256_add_ps (tv3.v, tv1.v);  /* now Imaginary part */
    tv3.v = _mm256_loadu_ps(&oArr_i[i]);   /* Output imags */
    tv3.v = _mm256_add_ps (tv3.v, tv2.v);  /* Update imag */
    _mm256_storeu_ps(&oArr_i[i], tv3.v);   /* Save imag */
  } /* end outer loop */
  ilast = i;  /* How far did I get? */
#else /* Scalar */
  ilast = 0;  /* Do all */
#endif
  /* Loop over whatever is left over */
  for (i=ilast; i<hiElem; i++) {
    tr1 = iArr1_r[i]; ti1 = iArr1_i[i];
    tr2 = iArr2_r[i]; ti2 = iArr2_i[i];
    oArr_r[i] += tr1*tr2 - ti1*ti2;
    oArr_i[i] += ti1*tr2 + tr1*ti2;
  } /* end loop over array */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadCplxMulAccum */

/**
 * Make arguments for a Threaded ThreadFAFunc?
 * \param thread     ObitThread object to be used
 * \param in         FA to be operated on
 * \param in2        2nd FA to be operated on
 * \param out        output FA
 * \param FA_3       Additional FA
 * \param FA_4       Additional FA
 * \param FA_5       Additional FA
 * \param FA_6       Additional FA
 * \param larg1      Length of function dependent arg1 in bytes
 * \param larg2      Length of function dependent arg2 in bytes
 * \param larg3      Length of function dependent arg3 in bytes
 * \param larg4      Length of function dependent arg4 in bytes
 * \param larg5      Length of function dependent arg5 in bytes
 * \param larg6      Length of function dependent arg6 in bytes
 * \param larg7      Length of function dependent arg7 in bytes
 * \param ThreadArgs[out] Created array of FAFuncArg, 
 *                   delete with KillFAFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeFAFuncArgs (ObitThread *thread, ObitFArray *in,
			     ObitFArray *in2, ObitFArray *out,
			     ObitFArray *FA_3, ObitFArray *FA_4,
			     ObitFArray *FA_5, ObitFArray *FA_6,
			     olong larg1, olong larg2, olong larg3, 
			     olong larg4, olong larg5,
			     olong larg6, olong larg7, 
			     FAFuncArg ***ThreadArgs)

{
  olong i, j, nThreads;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(FAFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(FAFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->thread= ObitThreadRef(thread);
    (*ThreadArgs)[i]->in    = ObitFArrayRef(in);
    if (in2) (*ThreadArgs)[i]->in2   = ObitFArrayRef(in2);
    else (*ThreadArgs)[i]->in2   = NULL;
    if (out)  (*ThreadArgs)[i]->out   = ObitFArrayRef(out);
    else (*ThreadArgs)[i]->out = NULL;
    if (FA_3) (*ThreadArgs)[i]->FA_3  = ObitFArrayRef(FA_3);
    else (*ThreadArgs)[i]->FA_3 = NULL;
    if (FA_4) (*ThreadArgs)[i]->FA_4  = ObitFArrayRef(FA_4);
    else (*ThreadArgs)[i]->FA_4 = NULL;
    if (FA_5) (*ThreadArgs)[i]->FA_5  = ObitFArrayRef(FA_5);
    else (*ThreadArgs)[i]->FA_5 = NULL;
    if (FA_6) (*ThreadArgs)[i]->FA_6  = ObitFArrayRef(FA_6);
    else (*ThreadArgs)[i]->FA_6 = NULL;
    (*ThreadArgs)[i]->first = 1;
    (*ThreadArgs)[i]->last  = in->arraySize;
    (*ThreadArgs)[i]->value = 0.0;
    for (j=0; j<MAXFARRAYDIM; j++) (*ThreadArgs)[i]->pos[j] = 0;
    if (larg1>0) (*ThreadArgs)[i]->arg1 = g_malloc0(larg1);
    else (*ThreadArgs)[i]->arg1 = NULL;
    if (larg2>0) (*ThreadArgs)[i]->arg2 = g_malloc0(larg2);
    else (*ThreadArgs)[i]->arg2 = NULL;
    if (larg3>0) (*ThreadArgs)[i]->arg3 = g_malloc0(larg3);
    else (*ThreadArgs)[i]->arg3 = NULL;
    if (larg4>0) (*ThreadArgs)[i]->arg4 = g_malloc0(larg4);
    else (*ThreadArgs)[i]->arg4 = NULL;
    if (larg5>0) (*ThreadArgs)[i]->arg5 = g_malloc0(larg5);
    else (*ThreadArgs)[i]->arg5 = NULL;
    if (larg5>0) (*ThreadArgs)[i]->arg6 = g_malloc0(larg6);
    else (*ThreadArgs)[i]->arg6 = NULL;
    if (larg5>0) (*ThreadArgs)[i]->arg7 = g_malloc0(larg7);
    else (*ThreadArgs)[i]->arg7 = NULL;
    (*ThreadArgs)[i]->ithread  = i;
  }

  return nThreads;
} /*  end MakeInterpImageArgs */

/**
 * Delete arguments for ThreadFAFunc
 * \param nargs      number of elements in ThreadArgs.
 * \param ThreadArgs Array of FAFuncArg
 */
static void KillFAFuncArgs (olong nargs, FAFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  ObitThreadPoolFree (ThreadArgs[0]->thread);  /* Free thread pool */
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      if (ThreadArgs[i]->thread) ObitThreadUnref(ThreadArgs[i]->thread);
      if (ThreadArgs[i]->in)   ObitFArrayUnref(ThreadArgs[i]->in);
      if (ThreadArgs[i]->in2)  ObitFArrayUnref(ThreadArgs[i]->in2);
      if (ThreadArgs[i]->out)  ObitFArrayUnref(ThreadArgs[i]->out);
      if (ThreadArgs[i]->FA_3) ObitFArrayUnref(ThreadArgs[i]->FA_3);
      if (ThreadArgs[i]->FA_4) ObitFArrayUnref(ThreadArgs[i]->FA_4);
      if (ThreadArgs[i]->FA_5) ObitFArrayUnref(ThreadArgs[i]->FA_5);
      if (ThreadArgs[i]->FA_6) ObitFArrayUnref(ThreadArgs[i]->FA_6);
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
} /*  end KillFAFuncArgs */
