/* $Id$         */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
gfloat*  ObitCArrayIndex (ObitCArray *in, olong *pos)
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
  for (i=0; i<in->arraySize; i+=2) 
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
  for (i=0; i<in->arraySize; i+=2) 
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
void ObitCArrayFill (ObitCArray* in, ofloat cmpx[2])
{
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->array != NULL);

  for (i=0; i<2*in->arraySize; i += 2) {
    in->array[i]   = cmpx[0];
    in->array[i+1] = cmpx[1];
  }
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
 *  Add corresponding elements of the array.
 *  out = in1 + in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitCArrayAdd (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
{
  olong i;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitCArrayIsCompatable(in1, in2));
  g_assert (ObitCArrayIsCompatable(in1, out));

  for (i=0; i<2*in1->arraySize; i++)
    out->array[i] = in1->array[i] + in2->array[i];
} /* end ObitCArrayAdd */

/**
 *  Subtract corresponding elements of the arrays.
 *  out = in1 - in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitCArraySub (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
{
  olong i;

   /* error checks */
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitCArrayIsCompatable(in1, in2));
  g_assert (ObitCArrayIsCompatable(in1, out));

  for (i=0; i<2*in1->arraySize; i++) 
    out->array[i] = in1->array[i] - in2->array[i];
 } /* end ObitCArraySub */

/**
 *  Multiply corresponding elements of the arrays.
 *  Output may be one of the inputs
 *  out = in1 * in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitCArrayMul (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
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
} /* end ObitCArrayMul */

/**
 *  Divide corresponding elements of the arrays.
 *  out = in1 / in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitCArrayDiv (ObitCArray* in1, ObitCArray* in2, ObitCArray* out)
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
} /* end ObitCArrayDiv */

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
  if (Fin->ndim!=Fin->ndim) return FALSE;

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
  olong i, j;
  ofloat fblank = ObitMagicF();

  /* error checks */
  g_assert (ObitCArrayIsA(in));
  g_assert (ObitFArrayIsA(out));
  g_assert (ObitCArrayIsFCompatable(in, out));

  /* Extract Amplitude */
  j = 0;
  for (i=0; i<2*in->arraySize; i += 2) {
    if ((in->array[i]!=fblank) && (in->array[i+1]!=fblank)) {
      out->array[j] = sqrt (in->array[i]*in->array[i] + 
			    in->array[i+1]*in->array[i+1]);
    } else {
      out->array[j] = fblank;
    }
    j++;
  }
}  /* end ObitCArrayAmp */

/**
 *  Return phase (radians) of elements of a CArray
 *  out = atan2(imag(in), real(in))
 * \param in  Input CArray
 * \param out Output FArray
 */
void ObitCArrayPhase (ObitCArray* in, ObitFArray* out)
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
}  /* end ObitCArrayPhase */

/**
 * In-place rearrangement of a center-at-the edges array to 
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
  theClass->ObitCArrayMakeC  = (ObitCArrayMakeCFP)ObitCArrayMakeF;
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
  theClass->ObitCArrayPhase  = (ObitCArrayPhaseFP)ObitCArrayPhase;
  theClass->ObitCArray2DCenter = (ObitCArray2DCenterFP)ObitCArray2DCenter;
  theClass->ObitCArrayAddConjg = (ObitCArrayAddConjgFP)ObitCArrayAddConjg;
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


