/* $Id$         */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2011                                          */
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

#include "ObitThread.h"
#include "ObitFArray.h"
#include "ObitMem.h"

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

/*--------------- File Global Variables  ----------------*/


/*---------------Private structures----------------*/
/* Threaded function argument */
typedef struct {
  /* ObitThread to use */
  ObitThread *thread;
  /* ObitFArray to work on */
  ObitFArray *in;
  /* First element (1-rel) number */
  olong        first;
  /* Highest element (1-rel) number */
  olong        last;
  /* Function dependent arguments */
  gpointer arg1, arg2, arg3, arg4, arg5;
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

/** Private: Make Threaded args */
static olong MakeFAFuncArgs (ObitThread *thread, ObitFArray *in,
			     olong larg1, olong larg2, olong larg3, 
			     olong larg4, olong larg5, 
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
      out->naxis[i] = MAX (1, MIN(naxis[i],524288));  /* Not too big */
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
      out->naxis[i] = MAX (1, MIN(naxis[i],524288));  /* Not too big */
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
  nThreads = MakeFAFuncArgs (in->thread, in, 0, 0, 0, 0, 0,
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
  nThreads = MakeFAFuncArgs (in->thread, in, 0, 0, 0, 0, 0,
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
  nThreads = MakeFAFuncArgs (in->thread, in, 0, 0, 0, 0, 0,
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
  olong i;
  ofloat fblank = ObitMagicF();

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Loop over array */
  for (i=0; i<in->arraySize; i++) 
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
  olong i, j, modeCell=0, imHalf=0, ipHalf=0, numCell, pos[MAXFARRAYDIM];
  olong i1, i2, ic, infcount;
  ofloat amax, amin, omax, omin, tmax, sum, sum2, x, count, mean, arg, cellFact=1.0;
  ofloat half, *histo = NULL, *thist=NULL;
  ofloat rawRMS, out = -1.0, fblank = ObitMagicF();
  gboolean done = FALSE;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  gboolean OK;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  /* get max/min */
  omax = ObitFArrayMax (in, pos);
  omin = ObitFArrayMin (in, pos);
  amax = omax;
  amin = omin;

  /* How many valid values? */
  count = 0; sum = sum2 = 0.0; 
  for (i=0; i<in->arraySize; i++) if (in->array[i]!=fblank) {
    count++;
    sum  += in->array[i];
    sum2 += in->array[i] * in->array[i];
  }
  if (count<5) return out; /* better have some */

  /* Get raw RMS */
  rawRMS = (sum2/count) - ((sum / count) * (sum / count));
  if (rawRMS>0.0) rawRMS = sqrt(rawRMS);
  if (rawRMS<(1.0e-5*MAX (fabs(omax), fabs(omin)))) return rawRMS;

  /* Make histogram size such that the average cell has 30 entries */
  numCell = count / 30;
  numCell = MAX (100, numCell);  /* but not too few */

  /* Initialize Threading */
  nThreads = 
    MakeFAFuncArgs (in->thread, in, 
		    sizeof(olong),  sizeof(olong), numCell*sizeof(ofloat),sizeof(ofloat), sizeof(ofloat), 
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
  while (!done) {

    /* Don't do this forever */
    infcount++;
    if (infcount>20) {
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
    count = *(olong*)(threadArgs[0])->arg2;
    histo = (ofloat*)(threadArgs[0]->arg3);
    for (i=1; i<nTh; i++) {
      count += *(olong*)(threadArgs[i])->arg2;
      thist = (ofloat*)(threadArgs[i]->arg3);
      for (j=0; j<numCell; j++) histo[j] += thist[j];
    }

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

    /* acceptability tests */
    /* ipHalf - imHalf must be greater than 10 and less than 50 */
    if ((ipHalf-imHalf) < 10) {
      /* if peak completely unresolved */
      if ((ipHalf-imHalf)<=0) {  /* wild stab */
	half = 0.5 / cellFact; /* ~ halfwidth? 1/2 cell */
      } else { /* partly resolved */
	half = (0.5 * (ipHalf-imHalf)) / cellFact; /* ~ halfwidth */
      }
      mean = amin + (modeCell-0.5) /  cellFact;
      /* don't spread over whole histogram - try for 25 cells*/
      half *= numCell / 25.0; 
      amax = mean + half;
      amin = mean - half;
      /* amin = MAX (amin, omin);
	 amax = MIN (amax, omax);*/
      continue;  /* try it again */
    } else if ((ipHalf-imHalf) > 50) {  /* Spread too thinly? */
      mean = amin + (modeCell-0.5) /  cellFact;
      half = (0.5 * (ipHalf-imHalf)) / cellFact; /* ~ halfwidth */
      /* don't spread over whole histogram - try for 25 cells*/
      half *= numCell / 25.0;  
      amax = mean + half;
      amin = mean - half;
      continue;  /* try it again */
    }
    done = TRUE;
    break;
  } /* end loop getting acceptable histogram */

  /* determine RMS */
  out = (ipHalf - imHalf) / cellFact; /* FWHM */
  out *= 2.35 / 2.0; /* Sigma */

  /* get second moment around mode +/- 3 sigma */
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
  olong i, count;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  ofloat sum, sum2, rawRMS, fblank = ObitMagicF();
  gboolean OK;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  /* Initialize Threading */
  nThreads = 
    MakeFAFuncArgs (in->thread, in, 
		    sizeof(olong), sizeof(ofloat), sizeof(ofloat), 0, 0,
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
      count += *(olong*)(threadArgs[i]->arg1);
      sum   += *(ofloat*)(threadArgs[i]->arg2);
      sum2  += *(ofloat*)(threadArgs[i]->arg3);
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
  olong i, count;
  olong nTh, nElem, loElem, hiElem, nElemPerThread, nThreads;
  ofloat sum, sum2, rawRMS, fblank = ObitMagicF();
  gboolean OK;
  FAFuncArg **threadArgs;

   /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  /* Initialize Threading */
  nThreads = 
    MakeFAFuncArgs (in->thread, in, 
		    sizeof(olong), sizeof(ofloat), sizeof(ofloat), 0, 0,
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
      count += *(olong*)(threadArgs[i]->arg1);
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
  olong i;
  ofloat quant, zero, center, rawRMS, mode, rangeFact, fblank = ObitMagicF();
  olong icell, modeCell=0, imHalf=0, ipHalf=0, numCell;
  olong i1, i2, ic, it;
  ofloat amax, amin, tmax, sum, sum2, x, count, mean, arg, cellFact=1.0;
  ofloat *histo = NULL;

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
  histo = ObitMemAllocName(numCell*sizeof(ofloat), "FArray histo");
  
  /* value range */
  amax = center + rangeFact*rawRMS;
  /* amin should be a quantization level */
  it = 0.5 + rangeFact*rawRMS/MAX(1.0e-20,quant);
  amin = center - it*quant;
  
  /* form histogram */
  cellFact =  1.0/quant;
  for (i=0; i<numCell; i++) histo[i] = 0.0;
  for (i=0; i<in->arraySize; i++) {
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
  olong i;
  ofloat delta, small, fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  delta = 1.0e20; small = 1.0e20;
  for (i=1; i<in->arraySize; i++) if (in->array[i]!=fblank) {
    if (in->array[i]!=in->array[i-1])
      delta = MIN (fabs(in->array[i]-in->array[i-1]), delta);
    if (fabs(in->array[i]) < fabs(small)) small = in->array[i];
  }

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
  olong i, icell, modeCell, numCell, count, pos[MAXFARRAYDIM];
  ofloat amax, amin, tmax, cellFact, *histo;
  ofloat out = 0.0, fblank = ObitMagicF();

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
  for (i=0; i<in->arraySize; i++) {
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
  olong i, count;
  ofloat out = 0.0, fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitFArrayIsA(in));
  g_assert (in->array != NULL);

  out = 0.0;
  count = 0;
  for (i=0; i<in->arraySize; i++) {
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
  olong i;

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) in->array[i] = scalar;
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
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) 
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
  olong i;
  ofloat out = 0.0, fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) 
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
  olong i, out = 0;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) 
    if (in->array[i]!=fblank) out++;

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
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) 
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
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<in->arraySize; i++) 
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
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if (in2->array[i]!=fblank)
      out->array[i] = in1->array[i];
    else out->array[i] = fblank;
  }
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
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = MAX (in1->array[i], in2->array[i]);
    else if (in1->array[i]==fblank)  /* 1 blanked */
      out->array[i] = in2->array[i];
    else if (in2->array[i]==fblank)
      out->array[i] = in1->array[i];  /* 2 blanked */
    else out->array[i] = fblank;      /* both blanked */
  }
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
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = MIN (in1->array[i], in2->array[i]);
    else if (in1->array[i]==fblank)  /* 1 blanked */
      out->array[i] = in2->array[i];
    else if (in2->array[i]==fblank)
      out->array[i] = in1->array[i];  /* 2 blanked */
    else out->array[i] = fblank;      /* both blanked */
  }
} /* end ObitFArrayMinArr */

/**
 * Pick the more extreme (furthest from zero) nonblanked elements of two arrays.
 *  out = MIN (in1, in2) or whichever is not blanked
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArrayExtArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) {
      if (fabs(in1->array[i])>fabs(in2->array[i]))
	out->array[i] = in1->array[i];
      else
	out->array[i] = in2->array[i];
    } else if (in1->array[i]==fblank)  /* 1 blanked */
      out->array[i] = in2->array[i];
    else if (in2->array[i]==fblank)
      out->array[i] = in1->array[i];  /* 2 blanked */
    else out->array[i] = fblank;      /* both blanked */
  }
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
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = (in1->array[i] + in2->array[i]);
    else if (in1->array[i]==fblank)  /* 1 blanked */
      out->array[i] = in2->array[i];
    else if (in2->array[i]==fblank)
      out->array[i] = in1->array[i];  /* 2 blanked */
    else out->array[i] = fblank;      /* both blanked */
  }
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
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = 0.5 * (in1->array[i] + in2->array[i]);
    else if (in1->array[i]==fblank)  /* 1 blanked */
      out->array[i] = in2->array[i];
    else if (in2->array[i]==fblank)
      out->array[i] = in1->array[i];  /* 2 blanked */
    else out->array[i] = fblank;      /* both blanked */
  }
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
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = in1->array[i] + in2->array[i];
    else out->array[i] = fblank;
  }
} /* end ObitFArrayAdd */

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
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = in1->array[i] - in2->array[i];
    else out->array[i] = fblank;
  }
 } /* end ObitFArraySub */

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
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)) 
      out->array[i] = in1->array[i] * in2->array[i];
    else out->array[i] = fblank;
  }
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
  ofloat fblank = ObitMagicF();

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitFArrayIsCompatable(in1, in2));
  g_assert (ObitFArrayIsCompatable(in1, out));

  for (i=0; i<in1->arraySize; i++) {
    if ((in1->array[i]!=fblank) && (in2->array[i]!=fblank)
	&& (in2->array[i]!=0.0)) 
      out->array[i] = in1->array[i] / in2->array[i];
    else out->array[i] = fblank;
  }
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
      arg = (x*x + y*y) * factor;
      if (arg<15.0) arg = exp (-arg);
      else arg = 0.0;
      data[indx] = arg;
    }
  }

} /* end ObitFArray2DCGauss */

/**
 * Make 2-D elliptical Gaussian in an FArray
 * Model is added to previoius contents of the FArray
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
      arg = x*x*factorX + y*y*factorY;
      if (arg<15.0) {
	arg = amp * exp (-arg);
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
  olong ix, iy, ip, np, hix, lox, hiy, loy, offx, offy;
  olong nx1, ny1, nx2, ny2, lenp1, lenp2, indx1, indx2;
  ofloat fblank = ObitMagicF();

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

  /* Size of in1/out */
  nx1 = in1->naxis[0];
  ny1 = in1->naxis[1];
  if (in1->ndim>2) np = in1->naxis[2];
  else np = 1;
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
  loy = MAX (0, pos1[1] - in2->naxis[1]/2);
  hiy = MIN (ny1-1, pos1[1] + in2->naxis[1]/2);
  /* In case in2 not centered */
  if (lox+offx<0) lox -= lox+offx;
  if (loy+offy<0) loy -= loy+offy;
  if (hix+offx>=nx2) hix -= hix+offx-nx2+1;
  if (hiy+offy>=ny2) hiy -= hiy+offy-ny2+1;
  /* Keep in range */
  lox = MAX (0, lox);
  hix = MIN (hix, nx1-1);
  loy = MAX (0, loy);
  hiy = MIN (hiy, ny1-1);

  /* Loop over planes */
  for (ip = 0; ip<np; ip++) {

    /* Loop over rows */
    for (iy=loy; iy<=hiy; iy++) {

      /* Loop over columns */
      for (ix=lox; ix<=hix; ix++) {

	/* indices in arrays */
	indx1 = ip*lenp1 + iy*nx1 + ix;
	indx2 = ip*lenp2 + (iy+offy) * nx2 + ix + offx;
	
	/* do operation */
	if (in1->array[indx1]!=fblank) {
	  if (in2->array[indx2]!=fblank) {
	    out->array[indx1] = in1->array[indx1] + scalar * in2->array[indx2];
	    /* debug
	    if (out->array[indx1] >1.0e5) fprintf (stderr,"bad accum %d  %d %g %g %g\n",
					       ix,iy,in1->array[indx1],in2->array[indx2],out->array[indx1]); */
	  }
	} else {
	  out->array[indx1] = fblank;
	}
      } 
    } 
  } /* end loop over planes */
 
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

  /* Initialize Threading */
  nThreads = MakeFAFuncArgs (in->thread, in, 
			     0, sizeof(olong), 3*sizeof(ofloat), 0, 0,
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
    for (i9=0; i9<naxis[8]; i9+=linc[8]) {
      opos[8]++;
      if(ndim>8) ipos[8] = blc[8] + i9;
      for (i8=0; i8<naxis[7]; i8+=linc[7]) {
	opos[7]++;
	if(ndim>7) ipos[7] = blc[7] + i8;
	for (i7=0; i7<naxis[6]; i7+=linc[6]) {
	  opos[6]++;
	  if(ndim>6) ipos[6] = blc[6] + i7;
	  for (i6=0; i6<naxis[5]; i6+=linc[5]) {
	    opos[5]++;
	    if(ndim>5) ipos[5] = blc[5] + i6;
	    for (i5=0; i5<naxis[4]; i5+=linc[4]) {
	      opos[4]++;
	      if(ndim>4) ipos[4] = blc[4] + i5;
	      for (i4=0; i4<naxis[3]; i4+=linc[3]) {
		opos[3]++;
		if(ndim>3) ipos[3] = blc[3] + i4;
		for (i3=0; i3<naxis[2]; i3+=linc[2]) {
		  opos[2]++;
		  if(ndim>2) ipos[2] = blc[2] + i3;
		  for (i2=0; i2<naxis[1]; i2+=linc[1]) {
		    opos[1]++;
		    if(ndim>1) ipos[1] = blc[1] + i2;
		    opos[0] = 0;
		    ipos[0] = blc[0];
		    
		    /* Array pointers */
		    inp  = ObitFArrayIndex (in,  ipos);
		    outp = ObitFArrayIndex (out, opos);

		    /* Copy row */
		    for (i1=0; i1<naxis[0]; i1+=linc[0]) {
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
  theClass->ObitFArraySub    = (ObitFArraySubFP)ObitFArraySub;
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
  if (in->array)  ObitMemFree(in->array);  in->array = NULL;
  if (in->naxis)  ObitMemFree(in->naxis); in->naxis = NULL;
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitFArrayClear */

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
  olong      hiElem     = largs->last-1;
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
  olong      hiElem     = largs->last-1;
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
  olong      hiElem     = largs->last-1;
  olong      *pos       = largs->pos;

  /* local */
  olong  i, temp, maxCell;
  ofloat maxAVal=0.0, maxVal, *data, val, fblank = ObitMagicF();

  largs->value = fblank;  /* In case no valid data*/
  if (hiElem<loElem) goto finish;

  /* Loop over array */
  maxCell = -1;
  maxVal = -1.0E25;
  data = in->array;
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
  olong      hiElem     = largs->last-1;

  /* local */
  olong i, count;
  ofloat sum, sum2, fblank = ObitMagicF();

  if (hiElem<loElem) goto finish;

  /* Loop over array */
  count = 0; sum = sum2 = 0.0; 
  for (i=loElem; i<hiElem; i++) 
    {
      if (in->array[i]!=fblank) {
	count++;
	sum  += in->array[i];
	sum2 += in->array[i] * in->array[i];
      }
    }

  /* Return values */
  *(olong*)largs->arg1  = count;
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
 * \li first    First element (1-rel) number
 * \li last     Highest element (1-rel) number
 * \li arg1     (olong) Number of cells in histogram
 * \li arg2     (olong) Count of pixels
 * \li arg3     (ofloat*) Histogram
 * \li arg4     (ofloat) Max in histogram
 * \li arg5     (ofloat) Min in histogram
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadFAHisto (gpointer arg)
{
  /* Get arguments from structure */
  FAFuncArg *largs = (FAFuncArg*)arg;
  ObitFArray *in        = largs->in;
  olong      loElem     = largs->first-1;
  olong      hiElem     = largs->last-1;
  olong      numCell    = *(olong*)largs->arg1;
  olong      count      = *(olong*)largs->arg2;
  ofloat     *histo     = (ofloat*)largs->arg3;
  ofloat     amax       = *(ofloat*)largs->arg4;
  ofloat     amin       = *(ofloat*)largs->arg5;

  /* local */
  olong  i, icell;
  ofloat cellFact, fblank = ObitMagicF();

  if (hiElem<loElem) goto finish;

  /* Loop over array */
  cellFact =  numCell / (amax - amin + 1.0e-20);
  count = 0;
  for (i=0; i<numCell; i++) histo[i] = 0.0;
  for (i=loElem; i<hiElem; i++) {
    if (in->array[i]!=fblank){
      icell = 0.5 + cellFact * (in->array[i]-amin);
      icell = MIN (numCell-1, MAX(0, icell));
      histo[icell]++;
      count ++;
    }
  }

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
    indx = loElem*in->naxis[0];  /* image array index */

    /* Loop over array convolving */
    for (iy = loElem; iy<=hiElem; iy++) {
      dy = iy - table[1];   /* y offset */
      for (ix = 0; ix<in->naxis[0]; ix++) {
	dx = ix - table[0];   /* x offset */
	farg = aa*dx*dx + bb*dy*dy + cc*dx*dy;
	if (farg<12.0) {
	  image[indx] += table[2] * exp(-farg);
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
 * Make arguments for a Threaded ThreadFAFunc?
 * \param thread     ObitThread object to be used
 * \param in         FA to be operated on
 * \param larg1      Length of function dependent arg1 in bytes
 * \param larg2      Length of function dependent arg2 in bytes
 * \param larg3      Length of function dependent arg3 in bytes
 * \param larg4      Length of function dependent arg4 in bytes
 * \param larg5      Length of function dependent arg5 in bytes
 * \param ThreadArgs[out] Created array of FAFuncArg, 
 *                   delete with KillFAFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeFAFuncArgs (ObitThread *thread, ObitFArray *in,
			     olong larg1, olong larg2, olong larg3, 
			     olong larg4, olong larg5, 
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
      if (ThreadArgs[i]->arg1) g_free(ThreadArgs[i]->arg1);
      if (ThreadArgs[i]->arg2) g_free(ThreadArgs[i]->arg2);
      if (ThreadArgs[i]->arg3) g_free(ThreadArgs[i]->arg3);
      if (ThreadArgs[i]->arg4) g_free(ThreadArgs[i]->arg4);
      if (ThreadArgs[i]->arg5) g_free(ThreadArgs[i]->arg5);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillFAFuncArgs */
