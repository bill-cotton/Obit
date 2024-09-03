/* $Id$         */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2020-2024                                          */
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

#include "ObitMatx.h"
#include "ObitMem.h"
#include "ObitSinCos.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitMatx.c
 * ObitMatx class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitMatx";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitMatxClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitMatxClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitMatxInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitMatxClear (gpointer in);

/** Private: Inatilize Class  */
void ObitMatxClassInit (void);

/** Private: Set Class function pointers. */
static void ObitMatxClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitMatx* newObitMatx (gchar* name)
{
  ObitMatx* out;

  /* Class initialization if needed 
  if (!myClassInfo.initialized) ObitMatxClassInit();*/

  /* allocate/init structure */
  out = ObitMemAlloc0Name(sizeof(ObitMatx), "ObitMatx");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitMatxInit((gpointer)out);

 return out;
} /* end newObitMatx */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitMatxGetClass (void)
{
  /* Class initialization if needed 
  if (!myClassInfo.initialized) ObitMatxClassInit();*/

  return (gconstpointer)&myClassInfo;
} /* end ObitMatxGetClass */

/**
 * Make a deep copy of an ObitMatx.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitMatx* ObitMatxCopy  (ObitMatx *in, ObitMatx *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, j, n, m;
  ofloat val[3];

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  //g_assert (ObitIsA(in, &myClassInfo));
  //if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitMatx(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* arrays */
  n = in->naxis[0]; m = in->naxis[1];  /* Only 2D */
  for (j=0; j<n; j++) {
    for (i=0; j<m; j++) {
      ObitMatxGet(in, i, j, val);
      ObitMatxSet(out, val, i, j);
    }
  }
  return out;
} /* end ObitMatxCopy */

/**
 * Determine if the two input objects have a compatable geometry.
 * Must be same type and geometry
 * \param in1  First object to test.
 * \param in2  Second object to test.
 * \return TRUE if compatable, else FALSE.
 */
gboolean ObitMatxIsCompatible (ObitMatx *in1, ObitMatx *in2)
{
  olong i;

  /* test fails if either input is NULL */
  if (in1==NULL) return FALSE;
  if (in2==NULL) return FALSE;
  /* Same type? */
  if (in1->type!=in2->type) return FALSE;

   /* Same number of non degenerate dimensions? */
  if (in1->ndim!=in2->ndim) {
    /* Don't bother if extra dimensions are length 1 */
    if (in1->ndim>in2->ndim) {  /* in1 has more */
      for (i=in2->ndim; i<in1->ndim; i++)
	if (in1->naxis[i]>1) return FALSE;
    } else {  /* in2 has more */
      for (i=in1->ndim; i<in2->ndim; i++)
	if (in2->naxis[i]>1) return FALSE;
    }
  }  /* end dimension check */

  return TRUE; /* must be OK */
} /* end ObitMatxIsCompatible  */

/**
 * Creates an ObitMatx of a specified type and geometry.
 * \param type  type, OBIT_Real, OBIT_Complex, OBIT_Double
 * \param ndim  Number of dimensions desired, ONLY 2 currently
 * \param naxis Dimensionality along each axis. NULL => don't create array.
 * \return the new object.
 */
ObitMatx* ObitMatxCreate (ObitMatxType type, olong ndim, olong *naxis)
{
  ObitMatx* out;
  olong i, size;

  /* Initialize class if needed */
  if (!myClassInfo.initialized) ObitMatxClassInit();

  /* Create basic structure */
  out = newObitMatx ("Matx");

  /* create data array if wanted */
  if ((ndim<0) || (naxis==NULL)) return out;
  g_assert (ndim==2); /* sanity check */

  /* copy geometry */
  out->ndim = ndim;
  out->naxis = ObitMemAlloc0Name(ndim*sizeof(olong), "Matx naxis");
  size = 1;
  for (i=0; i<ndim; i++) {
    out->naxis[i] = MAX (1, naxis[i]);
    size *= out->naxis[i]; /* total size */
  }

  out->type = type;  /* data type */

  switch (type) {
  case OBIT_Real:
    out->array = ObitMemAlloc0Name(size*sizeof(ofloat) + 
					out->naxis[0]*sizeof(ofloat),
					"Matx array");
    out->flt = (ofloat*)out->array;
   break;
  case OBIT_Complex:
    size *= 2;
    out->array = ObitMemAlloc0Name(size*sizeof(ofloat) + 
					out->naxis[0]*sizeof(ofloat),
					"Matx array");
    out->cpx = (ocomplex*)out->array;
    break;
  case OBIT_Double:
    size *= 2;
    out->array = ObitMemAlloc0Name(size*sizeof(ofloat) + 
					out->naxis[0]*sizeof(ofloat),
					"Matx array");
    out->dbl = (odouble*)out->array;
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };
  out->arraySize = size;

  return out;
} /* end ObitMatxCreate */

/**
 *  Matrix inner multiply
 *  out = in1 * in2
 * in1, in2 and out are square matrices of the same dimension
 * \param in1  Input object
 * \param in2  Input object
 * \param out  Output object
 */
void ObitMatxMult (ObitMatx* in1, ObitMatx* in2, ObitMatx* out)
{
  olong ir, ic, ii, n;
  olong indx1, indx2, indxo;
  ofloat sumr;
  odouble sumd;
  ocomplex tcpx;
  
  /* error checks 
  g_assert (ObitMatxIsCompatible(in1, in2));
  g_assert (ObitMatxIsCompatible(in1, out));
  g_assert (in1->naxis[0]==in1->naxis[1]);  */ /* Square */

  n = in1->naxis[0];
  switch (in1->type) {
  case OBIT_Real:
    for (ir=0; ir<n; ir++) {
      for (ic=0; ic<n; ic++) {
	sumr = 0.0;
	indxo = ic + ir*n;
	for (ii=0; ii<n; ii++) {
	  indx1 = (ii + ir*n);
	  indx2 = (ii*n + ic);
	  sumr += in1->flt[indx1]*in2->flt[indx2];
	}
	out->flt[indxo]   = sumr;
      }
    }
    break;
  case OBIT_Complex:
    for (ir=0; ir<n; ir++) {
      for (ic=0; ic<n; ic++) {
	indxo = (ic + ir*n);
	COMPLEX_SET (out->cpx[indxo], 0.0, 0.0);
	for (ii=0; ii<n; ii++) {
	  indx1 = (ii + ir*n);
	  indx2 = (ii*n + ic);
	  COMPLEX_MUL2(tcpx, in1->cpx[indx1], in2->cpx[indx2]);
	  COMPLEX_ADD2 (out->cpx[indxo], out->cpx[indxo], tcpx);
	}
      }
    }
    break;
  case OBIT_Double:
    for (ir=0; ir<n; ir++) {
      for (ic=0; ic<n; ic++) {
	sumd = 0.0;
	indxo = ic + ir*n;
	for (ii=0; ii<n; ii++) {
	  indx1 = (ii + ir*n);
	  indx2 = (ii*n + ic);
	  sumd += in1->dbl[indx1]*in2->dbl[indx2];
	}
	out->dbl[indxo]   = sumd;
      }
    }
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };

 /* end row loop */
  
} /* end ObitMatxMatrixMult */

/**
 *  Matrix inner multiply of a matrix  by the conjugate transpose of another matrix
 *  out = in1 * in2^H
 * in1, in2 and out are square matrices of the same dimension
 * \param in1  Input object
 * \param in2  Input object, conjugate transposed used
 * \param out  Output object
 */
void ObitMatxMultCT (ObitMatx* in1, ObitMatx* in2, ObitMatx* out)
{
  olong ir, ic, ii, n;
  olong indx1, indx2, indxo;
  ofloat sumr;
  odouble sumd;
  ocomplex tcpx, tcpx2;
  
  /* error checks 
  g_assert (ObitMatxIsCompatible(in1, in2));
  g_assert (ObitMatxIsCompatible(in1, out));
  g_assert (in1->naxis[0]==in1->naxis[1]);*/  /* Square */

  n = in1->naxis[0];
  switch (in1->type) {
  case OBIT_Real:
    for (ir=0; ir<n; ir++) {
      for (ic=0; ic<n; ic++) {
	sumr = 0.0;
	indxo = ic + ir*n;
	for (ii=0; ii<n; ii++) {
	  indx1 = (ii + ir*n);
	  indx2 = (ii + ic*n);
	  sumr += in1->flt[indx1]*in2->flt[indx2];
	}
	out->flt[indxo]   = sumr;
      }
    }
    break;
  case OBIT_Complex:
    for (ir=0; ir<n; ir++) {
      for (ic=0; ic<n; ic++) {
	indxo = (ic + ir*n);
	COMPLEX_SET (out->cpx[indxo], 0.0, 0.0);
	for (ii=0; ii<n; ii++) {
	  indx1 = (ii + ir*n);
	  indx2 = (ii + ic*n);
	  COMPLEX_CONJUGATE(tcpx2, in2->cpx[indx2]);
	  COMPLEX_MUL2(tcpx, in1->cpx[indx1], tcpx2);
	  COMPLEX_ADD2 (out->cpx[indxo], out->cpx[indxo], tcpx);
	}
      }
    }
    break;
  case OBIT_Double:
    for (ir=0; ir<n; ir++) {
      for (ic=0; ic<n; ic++) {
	sumd = 0.0;
	indxo = ic + ir*n;
	for (ii=0; ii<n; ii++) {
	  indx1 = (ii + ir*n);
	  indx2 = (ii + ic*n);
	  sumd += in1->dbl[indx1]*in2->dbl[indx2];
	}
	out->dbl[indxo]   = sumd;
      }
    }
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };

 /* end row loop */
  
} /* end ObitMatxMultCT */

/**
 *  Add corresponding elements of the array.
 *  out = in1 + in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitMatxAdd (ObitMatx* in1, ObitMatx* in2, ObitMatx* out)
{
  olong i;

   /* error checks */
  g_assert (ObitMatxIsCompatible(in1, in2));
  g_assert (ObitMatxIsCompatible(in1, out));

  switch (in1->type) {
  case OBIT_Real:
   for (i=0; i<in1->arraySize; i++)
    out->flt[i] = in1->flt[i] + in2->flt[i];
   break;
  case OBIT_Complex:
   for (i=0; i<in1->arraySize/2; i++)
     COMPLEX_ADD2(out->cpx[i], in1->cpx[i], in2->cpx[i]);
    break;
  case OBIT_Double:
   for (i=0; i<in1->arraySize/2; i++)
     out->dbl[i] = in1->dbl[i] + in2->dbl[i];
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };

} /* end ObitMatxAdd */


/**
 *  Subtract corresponding elements of the arrays.
 *  out = in1 - in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitMatxSub (ObitMatx* in1, ObitMatx* in2, ObitMatx* out)
{
  olong i;

   /* error checks */
  g_assert (ObitMatxIsCompatible(in1, in2));
  g_assert (ObitMatxIsCompatible(in1, out));

  switch (in1->type) {
  case OBIT_Real:
   for (i=0; i<in1->arraySize; i++)
    out->flt[i] = in1->flt[i] - in2->flt[i];
   break;
  case OBIT_Complex:
   for (i=0; i<in1->arraySize/2; i++)
     COMPLEX_SUB(out->cpx[i], in1->cpx[i], in2->cpx[i]);
    break;
  case OBIT_Double:
   for (i=0; i<in1->arraySize/2; i++)
     out->dbl[i] = in1->dbl[i] - in2->dbl[i];
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };
 } /* end ObitMatxSub */

/**
 *   Multiply a scalar times an ObitMatx.
 *  out = in1 * scalar
 * \param in1  Input object with data
 * \param scalar to multiply each element with
 * \param out  Output object
 */
void ObitMatxSMult (ObitMatx* in1, void *scalar, ObitMatx* out)
{
  olong i;

   /* error checks */
  g_assert (ObitMatxIsCompatible(in1, out));

  switch (in1->type) {
  case OBIT_Real:
    for (i=0; i<in1->arraySize; i++)
      out->flt[i] = in1->flt[i] * ((ofloat*)scalar)[0];
    break;
  case OBIT_Complex:
    for (i=0; i<in1->arraySize/2; i++)
      COMPLEX_MUL2(out->cpx[i], in1->cpx[i], ((ocomplex*)scalar)[0]);
    break;
  case OBIT_Double:
    for (i=0; i<in1->arraySize/2; i++)
      out->dbl[i] = in1->dbl[i] * ((odouble*)scalar)[0];
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };
 } /* end ObitMatxSMult */


/**
 *  Transpose and conjugate a matrix.
 * \param in      Input object
 * \param out     Output object
 */
void ObitMatxCTrans (ObitMatx* in, ObitMatx* out)
{
  olong i, j, n, m;
  ofloat val[3];
  gboolean conjg = in->type==OBIT_Complex;
  /* error checks */
  g_assert (ObitMatxIsCompatible(in, out));
  n = in->naxis[0]; m = in->naxis[1];  /* Only 2D */
  for (j=0; j<n; j++) {
    for (i=0; i<m; i++) {
      ObitMatxGet(in, i, j, val);
      if (conjg) val[1] = -val[1];
      ObitMatxSet(out, val, j, i);
    }
  }
} /* end ObitMatxCTrans */


/**
 * Zero the elements of an array 
 * \param in      Input object to zero
 */
void ObitMatxZero (ObitMatx* in)
{
  olong i;
  ofloat val[3] = {0.0,0.0,0.0};

  /* error checks */
  //g_assert (ObitIsA(in, &myClassInfo));
  //g_assert (in->array != NULL);

  switch (in->type) {
  case OBIT_Real:
    for (i=0; i<in->arraySize; i++) in->flt[i] = 0.0;
   break;
  case OBIT_Complex:
    for (i=0; i<in->arraySize/2; i++) 
      COMPLEX_SET(in->cpx[i], val[0], val[1]); 
    break;
  case OBIT_Double:
    for (i=0; i<in->arraySize/2; i++) in->dbl[i] = 0.0;
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };
} /* end ObitMatxZero */

/**
 * Set the matrix to a unit matrix - must be square 2D matrix
 * \param in      Input object
 */
void ObitMatxUnit (ObitMatx* in)
{
  olong i, n;
  ofloat val[3] = {1.0,0.0,0.0};

  n = in->naxis[0]; /* dimension */
  switch (in->type) {
  case OBIT_Real:
    for (i=0; i<in->arraySize; i++) in->flt[i] = 0.0;
    /* Diagonal terms */
    for (i=0; i<n; i++) in->flt[i*n+i] = 1.0;
    break;
  case OBIT_Complex:
    /* Zero */
    for (i=0; i<in->arraySize/2; i++) COMPLEX_SET(in->cpx[i], val[1], val[2]); 
    /* Diagonal terms */
    for (i=0; i<n; i++) COMPLEX_SET(in->cpx[i*n+i], val[0], val[1]); 
    break;
  case OBIT_Double:
    for (i=0; i<in->arraySize; i++) in->dbl[i] = 0.0;
    /* Diagonal terms */
    for (i=0; i<n; i++) in->dbl[i*n+i] = 1.0;
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };
} /* end ObitMatxUnit */

/**
 * Fills values in a 2x2 complex ObitMatx as 1D
 * \param matx  ObitMatx to fill
 * \param R00   Real part of element (0,0)
 * \param I00   Imaginary part of element (0,0)
 * \param R10   Real part of element (1,0)
 * \param I10   Imaginary part of element (1,0)
 * \param R01   Real part of element (0,1)
 * \param I01   Imaginary part of element (0,1)
 * \param R11   Real part of element (1,1)
 * \param I11   Imaginary part of element (1,1)
 */
void ObitMatxSet2C(ObitMatx *matx, ofloat R00, ofloat I00, ofloat R10, ofloat I10, 
		   ofloat R01, ofloat I01, ofloat R11, ofloat I11)
{
  olong indx;

  indx = 0; COMPLEX_SET (matx->cpx[indx], R00, I00);
  indx = 1; COMPLEX_SET (matx->cpx[indx], R10, I10);
  indx = 2; COMPLEX_SET (matx->cpx[indx], R01, I01);
  indx = 3; COMPLEX_SET (matx->cpx[indx], R11, I11);
} /* end ObitMatxSet2C */

/**
 * Fetch values in a 2x2 complex ObitMatx as 1D
 * \param matx  ObitMatx to fetch values from
 * \param R00   Real part of element (0,0)
 * \param I00   Imaginary part of element (0,0)
 * \param R10   Real part of element (1,0)
 * \param I10   Imaginary part of element (1,0)
 * \param R01   Real part of element (0,1)
 * \param I01   Imaginary part of element (0,1)
 * \param R11   Real part of element (1,1)
 * \param I11   Imaginary part of element (1,1)
 */
void ObitMatxGet2C(ObitMatx *matx, ofloat* R00, ofloat *I00, ofloat *R10, ofloat *I10, 
		   ofloat *R01, ofloat *I01, ofloat *R11, ofloat *I11)
{
  olong indx;
  indx = 0; *R00 = matx->cpx[indx].real; *I00 = matx->cpx[indx].imag;
  indx = 1; *R10 = matx->cpx[indx].real; *I10 = matx->cpx[indx].imag;
  indx = 2; *R01 = matx->cpx[indx].real; *I01 = matx->cpx[indx].imag;
  indx = 3; *R11 = matx->cpx[indx].real; *I11 = matx->cpx[indx].imag;
} /* end ObitMatxGet2C */

/**
 * Invert 2x2 matrix 
 * \param in    ObitMatx to invert
 * \param out   Inverted ObitMatx
 */
/** Invert 2x2 matrix */
void ObitMatxInv2x2(ObitMatx *in, ObitMatx *out)
{
  ofloat fDet, Det[2], d, *matx;
  odouble dDet;
  /* error checks */
  g_assert (ObitMatxIsCompatible(in, out));
  g_assert (in->naxis[0]==2 && out->naxis[1]==2);  /* 2x2 */

  switch (in->type) {
  case OBIT_Real:
    fDet = in->flt[0]*in->flt[3] - in->flt[1]*in->flt[2];
    if (fDet>0.0) fDet = 1.0 / fDet;
    else          fDet = 1.0;
    out->flt[0] = +in->flt[3]*fDet;
    out->flt[1] = -in->flt[1]*fDet;
    out->flt[2] = -in->flt[2]*fDet;
    out->flt[3] = +in->flt[0]*fDet;
    break;
  case OBIT_Complex:
    /* inverse of determinant */
    matx = in->array;
    Det[0] = (matx[0]*matx[6] - matx[1]*matx[7]) - (matx[2]*matx[4] - matx[3]*matx[5]);
    Det[1] = (matx[0]*matx[7] + matx[1]*matx[6]) - (matx[2]*matx[5] + matx[3]*matx[4]);
    /* Inverse of determinant */
    d = Det[0]*Det[0] + Det[1]*Det[1];
    if (d!=0.0) d = 1.0 / d;
    else d = 1.0;
    Det[0] *=  d;
    Det[1] *= -d;
    /* Inverse matrix out */
    COMPLEX_SET (out->cpx[3],   matx[0]*Det[0]-matx[1]*Det[1],    matx[0]*Det[1]+matx[1]*Det[0]);
    COMPLEX_SET (out->cpx[2], -(matx[4]*Det[0]-matx[5]*Det[1]), -(matx[4]*Det[1]+matx[5]*Det[0]));
    COMPLEX_SET (out->cpx[1], -(matx[2]*Det[0]-matx[3]*Det[1]), -(matx[2]*Det[1]+matx[3]*Det[0]));
    COMPLEX_SET (out->cpx[0],   matx[6]*Det[0]-matx[7]*Det[1],    matx[6]*Det[1]+matx[7]*Det[0]);
    break;
  case OBIT_Double:
    dDet = in->dbl[0]*in->dbl[3] - in->dbl[1]*in->dbl[2];
    if (dDet>0.0) dDet = 1.0 / dDet;
    else          dDet = 1.0;
    out->dbl[0] = +in->dbl[3]*dDet;
    out->dbl[1] = -in->dbl[1]*dDet;
    out->dbl[2] = -in->dbl[2]*dDet;
    out->dbl[3] = +in->dbl[0]*dDet;
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };

}  /* end ObitMatxInv2x2 */

/** Inverse perfect linear feed Jones matrix complex 2x2 */
void ObitMatxIPerfLinJones(ObitMatx *out)
//void IJones(ofloat J[8])
{
  ofloat elp_x=0.0, elp_y=0.0, ori_x=0.0, ori_y=G_PI*0.5;
  ofloat angle[4], sina[4], cosa[4], Jones[8] ,Det[2], d;
  /*  Check that 2x2 complex? */

  angle[0] = G_PI*0.25+elp_x; angle[1] = G_PI*0.25-elp_y;
  angle[2] = ori_x;           angle[3] = ori_y;
  ObitSinCosVec(4, angle, sina, cosa);
  Jones[0] =  cosa[0]*cosa[2]; Jones[1] = -cosa[0]*sina[2];
  Jones[2] =  sina[0]*cosa[2]; Jones[3] =  sina[0]*sina[2];
  Jones[4] =  sina[1]*cosa[3]; Jones[5] = -sina[1]*sina[3];
  Jones[6] =  cosa[1]*cosa[3]; Jones[7] =  cosa[1]*sina[3];
  /* inverse of determinant */
  Det[0] = (Jones[0]*Jones[6] - Jones[1]*Jones[7]) - (Jones[2]*Jones[4] - Jones[3]*Jones[5]);
  Det[1] = (Jones[0]*Jones[7] + Jones[1]*Jones[6]) - (Jones[2]*Jones[5] + Jones[3]*Jones[4]);
  /* Inverse of determinant */
  d = Det[0]*Det[0] + Det[1]*Det[1];
  if (d!=0.0) d = 1.0 / d;
  else d = 1.0;
  Det[0] *=  d;
  Det[1] *= -d;
  /* Inverse matrix out */
  COMPLEX_SET (out->cpx[3],   Jones[0]*Det[0]-Jones[1]*Det[1],    Jones[0]*Det[1]+Jones[1]*Det[0]);
  COMPLEX_SET (out->cpx[2], -(Jones[4]*Det[0]-Jones[5]*Det[1]), -(Jones[4]*Det[1]+Jones[5]*Det[0]));
  COMPLEX_SET (out->cpx[1], -(Jones[2]*Det[0]-Jones[3]*Det[1]), -(Jones[2]*Det[1]+Jones[3]*Det[0]));
  COMPLEX_SET (out->cpx[0],   Jones[6]*Det[0]-Jones[7]*Det[1],    Jones[6]*Det[1]+Jones[7]*Det[0]);
} /* end ObitMatxIPerfLinJonesJones */

/** Inverse perfect circular feed Jones matrix  complex 2x2 
    Really just the perfect linear Jones matrix */
void ObitMatxIPerfCirJones(ObitMatx *out)
{
  ofloat elp_x=0.0, elp_y=0.0, ori_x=0.0, ori_y=G_PI*0.5;
  ofloat angle[4], sina[4], cosa[4], Jones[8];

  angle[0] = G_PI*0.25+elp_x; angle[1] = G_PI*0.25-elp_y;
  angle[2] = ori_x;           angle[3] = ori_y;
  ObitSinCosVec(4, angle, sina, cosa);

  Jones[0] =  cosa[0]*cosa[2]; Jones[1] = -cosa[0]*sina[2];
  Jones[2] =  sina[0]*cosa[2]; Jones[3] =  sina[0]*sina[2];
  Jones[4] =  sina[1]*cosa[3]; Jones[5] = -sina[1]*sina[3];
  Jones[6] =  cosa[1]*cosa[3]; Jones[7] =  cosa[1]*sina[3];
  
  /* Matrix out */
  COMPLEX_SET (out->cpx[0], Jones[0], Jones[1]);
  COMPLEX_SET (out->cpx[1], Jones[2], Jones[3]);
  COMPLEX_SET (out->cpx[2], Jones[4], Jones[5]);
  COMPLEX_SET (out->cpx[3], Jones[6], Jones[7]); 
} /* end ObitMatxIPerfCirJones */

/**
 *  Outer 2x2 complex multiply, original
 *  out = in1 (outer product) conjg(in2)
 * \param iin1  1st Input object
 * \param iin2  2nd Input object
 * \param oout  Output object
 */
void ObitMatxOuterMult2C(ObitMatx *iin1, ObitMatx *iin2, ObitMatx *oout)
//static void MatxOuter(ofloat in1[8], ofloat in2[8], ofloat out[32])
{
  ofloat *in1, *in2, *out;
  /* Local ofloat pointers */
  in1 = (ofloat*)iin1->array;
  in2 = (ofloat*)iin2->array;
  out = (ofloat*)oout->array;

  out[0]  =  in1[0] * in2[0] + in1[1] * in2[1];  out[1]  =  in1[1] * in2[0] - in1[0] * in2[1];
  out[2]  =  in1[0] * in2[2] + in1[1] * in2[3];  out[3]  =  in1[1] * in2[2] - in1[0] * in2[3];
  out[4]  =  in1[2] * in2[0] + in1[3] * in2[1];  out[5]  =  in1[3] * in2[0] - in1[2] * in2[1];
  out[6]  =  in1[2] * in2[2] + in1[3] * in2[3];  out[7]  =  in1[3] * in2[2] - in1[2] * in2[3];
  
  out[8]  =  in1[0] * in2[4] + in1[1] * in2[5];  out[9]  =  in1[1] * in2[4] - in1[0] * in2[5];
  out[10] =  in1[0] * in2[6] + in1[1] * in2[7];  out[11] =  in1[1] * in2[6] - in1[0] * in2[7];
  out[12] =  in1[2] * in2[4] + in1[3] * in2[5];  out[13] =  in1[3] * in2[4] - in1[2] * in2[5];
  out[14] =  in1[2] * in2[6] + in1[3] * in2[7];  out[15] =  in1[3] * in2[6] - in1[2] * in2[7];
  
  out[16] =  in1[4] * in2[0] + in1[5] * in2[1];  out[17] =  in1[5] * in2[0] - in1[4] * in2[1];
  out[18] =  in1[4] * in2[2] + in1[5] * in2[3];  out[19] =  in1[5] * in2[2] - in1[4] * in2[3];
  out[20] =  in1[6] * in2[0] + in1[7] * in2[1];  out[21] =  in1[7] * in2[0] - in1[6] * in2[1];
  out[22] =  in1[6] * in2[2] + in1[7] * in2[3];  out[23] =  in1[7] * in2[2] - in1[6] * in2[3];
  
  out[24] =  in1[4] * in2[4] + in1[5] * in2[5];  out[25] =  in1[5] * in2[4] - in1[4] * in2[5];
  out[26] =  in1[4] * in2[6] + in1[5] * in2[7];  out[27] =  in1[5] * in2[6] - in1[4] * in2[7];
  out[28] =  in1[6] * in2[4] + in1[7] * in2[5];  out[29] =  in1[7] * in2[4] - in1[6] * in2[5];
  out[30] =  in1[6] * in2[6] + in1[7] * in2[7];  out[31] =  in1[7] * in2[6] - in1[6] * in2[7];
  } /* end  MatxOuterMult2C*/

  /**
 *  4x4 complex matrix * 4x1 complex vector multiply 
 *  out = in1 x in2
 * \param iin1  1st Input object
 * \param iin2  2nd Input object
 * \param oout  Output object
 */
void ObitMatxVec4Mult(ObitMatx *iin1, ObitMatx *iin2, ObitMatx *oout) {
//  static void ObitMatxVec4Mult(ofloat in1[32], ofloat in2[8], ofloat out[8])
  olong ic, ii, n=4;
  ofloat sumr, sumi;
  ofloat *in1, *in2, *out;
  /* Local ofloat pointers */
  in1 = (ofloat*)iin1->array;
  in2 = (ofloat*)iin2->array;
  out = (ofloat*)oout->array;
  /* Check */
  g_assert ((iin1->naxis[0]==4) && (iin1->naxis[1]==4) && (iin1->type==OBIT_Complex));
  g_assert ((iin2->naxis[0]==4) && (iin2->naxis[1]==1) && (iin2->type==OBIT_Complex));
  g_assert ((oout->naxis[0]==4) && (oout->naxis[1]==1) && (oout->type==OBIT_Complex));

  for (ic=0; ic<4; ic++) {
    sumr = sumi = 0.0;
    for (ii=0; ii<4; ii++) {
      sumr += in1[(ic*n+ii)*2]*in2[ii*2]   - in1[(ic*n+ii)*2+1]*in2[ii*2+1];
      sumi += in1[(ic*n+ii)*2]*in2[ii*2+1] + in1[(ic*n+ii)*2+1]*in2[ii*2];
    }
    out[ic*2]   = sumr;
    out[ic*2+1] = sumi;
  }
} /* end ObitMatxVec4Mult */

/**
 * Print a 2x2 matrix to stderr
 * \param in      Input Matrix
 * \param label   Label
 */
void ObitMatxPrint (ObitMatx* in, gchar* label)
{
  ofloat R00, I00, R01, I01, R10, I10, R11, I11;
  odouble D00, D01, D10, D11;

  switch (in->type) {
  case OBIT_Real:
    ObitMatxGet(in, 0, 0, &R00); ObitMatxGet(in, 1, 0, &R10); 
    ObitMatxGet(in, 1, 0, &R01); ObitMatxGet(in, 1, 1, &R11); 
    fprintf (stderr, " %s \n   %10.6f %10.6f\n   %10.6f %10.6f\n", label, R00, R10, R01, R11);
   break;
  case OBIT_Complex:
    ObitMatxGet2C (in, &R00, &I00, &R10, &I10, &R01, &I01, &R11, &I11);
    fprintf (stderr, " %s \n   %10.6f %10.6f, %10.6f %10.6f\n   %10.6f %10.6f, %10.6f %10.6f\n", 
	     label, R00, I00, R10, I10, R01, I01, R11, I11);
    break;
  case OBIT_Double:
    ObitMatxGet(in, 0, 0, &D00); ObitMatxGet(in, 1, 0, &D10); 
    ObitMatxGet(in, 1, 0, &D01); ObitMatxGet(in, 1, 1, &D11); 
    fprintf (stderr, " %s \n   %12.8lf %12.8lf\n   %12.8lf %12.8lf\n", label, D00, D10, D01, D11);
    break;
  default:
    g_assert_not_reached(); /* unknown, barf */
  };
 
} /* end ObitMatxPrint */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitMatxClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitMatxClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitMatxClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitMatxClassInfoDefFn (gpointer inClass)
{
  ObitMatxClassInfo *theClass = (ObitMatxClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitMatxClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitMatxClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitMatxGetClass;
  theClass->newObit       = (newObitFP)newObitMatx;
  theClass->ObitCopy      = (ObitCopyFP)ObitMatxCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitMatxClear;
  theClass->ObitInit      = (ObitInitFP)ObitMatxInit;
  theClass->ObitMatxCreate = (ObitMatxCreateFP)ObitMatxCreate;
  theClass->ObitMatxIsCompatible = 
    (ObitMatxIsCompatibleFP)ObitMatxIsCompatible;
  theClass->ObitMatxMult   = (ObitMatxMultFP)ObitMatxMult;
  theClass->ObitMatxMultCT = (ObitMatxMultCTFP)ObitMatxMultCT;
  theClass->ObitMatxAdd    = (ObitMatxAddFP)ObitMatxAdd;
  theClass->ObitMatxSub    = (ObitMatxSubFP)ObitMatxSub;
  theClass->ObitMatxCTrans = (ObitMatxCTransFP)ObitMatxCTrans;
  theClass->ObitMatxZero   = (ObitMatxZeroFP)ObitMatxZero;
  theClass->ObitMatxUnit   = (ObitMatxUnitFP)ObitMatxUnit;
  theClass->ObitMatxSet2C  = (ObitMatxSet2CFP)ObitMatxSet2C;
  theClass->ObitMatxGet2C  = (ObitMatxGet2CFP)ObitMatxGet2C;
  theClass->ObitMatxIPerfLinJones = (ObitMatxIPerfLinJonesFP)ObitMatxIPerfLinJones;
  theClass->ObitMatxOuterMult2C   = (ObitMatxOuterMult2CFP)ObitMatxOuterMult2C;
} /* end ObitMatxClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitMatxInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitMatx *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->array        = NULL;
  in->arraySize    = 0;
  in->ndim         = 0;
  in->naxis        = NULL;

} /* end ObitMatxInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitMatx* cast to an Obit*.
 */
void ObitMatxClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitMatx *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->array)  in->array = ObitMemFree(in->array);
  if (in->naxis)  in->naxis = ObitMemFree(in->naxis);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitMatxClear */

