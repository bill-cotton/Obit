/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2020                                               */
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
/** 
 * Simple Matrix library
 */
#include "Obit.h"
#include "ObitErr.h"
#include "ObitComplex.h"

#ifndef OBITMATX_H 
#define OBITMATX_H 
/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitMatxType
 * enum for matrix data type.
 */
enum obitMatxType{
  /** Real*/
  OBIT_Real, 
  /** Complex */
  OBIT_Complex,  
  /** Double */
  OBIT_Double,
}; /* end enum obitMatxType */
typedef enum obitMatxType ObitMatxType;
/*--------------Class definitions-------------------------------------*/
/** ObitXML Class structure. */
typedef struct {
#include "ObitMatxDef.h"
} ObitMatx;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitMatx
 * returns a ObitMatx*.
 * in = object to unreference
 */
#define ObitMatxUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitMatx.
 * returns a ObitMatx*.
 * in = object to reference
 */
#define ObitMatxRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitMatxIsA(in) ObitIsA (in, ObitMatxGetClass())

/**
 * Routine for setting an element
 * \li matx = ObitMatx structure modify.
 * \li i1, i2 = 0 rel indices
 * \li val  = array of data (real, complex, double) for output
 */
static inline void ObitMatxGet(ObitMatx *matx, olong i1, olong i2, void *val) {
  if (matx->type==OBIT_Real) {
    *(ofloat*)val = matx->flt[i1+i2*matx->naxis[0]];
  } else if (matx->type==OBIT_Complex) {
    ((ofloat*)val)[0] = matx->cpx[i1+i2*matx->naxis[0]].real;
    ((ofloat*)val)[1] = matx->cpx[i1+i2*matx->naxis[0]].imag;
    COMPLEX_SET(matx->cpx[i1+i2*matx->naxis[0]], 
		((ofloat*)val)[0], ((ofloat*)val)[1]);
  } else { /* Double */
    *((odouble*)val) = matx->dbl[i1+i2*matx->naxis[0]];
  }
} /* end ObitMatxGet */

/**
 * Routine for setting an element
 * \li matx = ObitMatx structure modify.
 * \li val  = array of data (real, complex, double)
 * \li i1, i2 = 0 rel indices
 */
static inline void ObitMatxSet(ObitMatx *matx, void *val, olong i1, olong i2) {
  if (matx->type==OBIT_Real) {
    matx->flt[i1+i2*matx->naxis[0]] = ((float*)val)[0];
  } else if (matx->type==OBIT_Complex) {
    COMPLEX_SET(matx->cpx[i1+i2*matx->naxis[0]], 
		((ofloat*)val)[0], ((ofloat*)val)[1]); 
  } else { /* Double */
    matx->dbl[i1+i2*matx->naxis[0]] = *((odouble*)val);}
} /* end ObitMatxSet */

/** Create empty (0) matrix */
ObitMatx* ObitMatxCreate(ObitMatxType type, olong ndim, olong *naxis);
typedef ObitMatx* (*ObitMatxCreateFP) (ObitMatxType type, olong ndim, olong *naxis);
/** Copy matrix */
ObitMatx* ObitMatxCopy(ObitMatx *in, ObitMatx *out, ObitErr *err);
typedef ObitMatx* (*ObitMatxCopyFP) (ObitMatx *in, ObitMatx *out, ObitErr *err);
/** Are matrices compatible? */
gboolean ObitMatxIsCompatible(ObitMatx *in1, ObitMatx *in2);
typedef gboolean (*ObitMatxIsCompatibleFP) (ObitMatx *in1, ObitMatx *in2);
/** Multiply */
void ObitMatxMult(ObitMatx *in1, ObitMatx *in2, ObitMatx *out);
typedef void (*ObitMatxMultFP) (ObitMatx *in1, ObitMatx *in2, ObitMatx *out);
/** Add */
void ObitMatxAdd(ObitMatx *in1, ObitMatx *in2, ObitMatx *out);
typedef void (*ObitMatxAddFP) (ObitMatx *in1, ObitMatx *in2, ObitMatx *out);
/** Subtract in1-in2 */
void ObitMatxSub(ObitMatx *in1, ObitMatx *in2, ObitMatx *out);
typedef void (*ObitMatxSubFP) (ObitMatx *in1, ObitMatx *in2, ObitMatx *out);
/** Conjugate transpose */
void ObitMatxCTrans(ObitMatx *in, ObitMatx *out);
typedef void (*ObitMatxCTransFP) (ObitMatx *in, ObitMatx *out);
/** Zero values */
void ObitMatxZero(ObitMatx *out);
typedef void (*ObitMatxZeroFP) (ObitMatx *in);
/** Fill 2x2 */
void ObitMatxSet2C(ObitMatx *matx, ofloat R00, ofloat I00, ofloat R01, ofloat I01, 
		   ofloat R10, ofloat I10, ofloat R11, ofloat I11);
typedef void (*ObitMatxSet2CFP) (ObitMatx *matx, ofloat R00, ofloat I00, 
				 ofloat R01, ofloat I01,  ofloat R10, ofloat I10, 
				 ofloat R11, ofloat I11);
/** Inverse perfect linear feed Jones matrix */
void ObitMatxIPerfLinJones(ObitMatx *in);
typedef void (*ObitMatxIPerfLinJonesFP) (ObitMatx *out);
/** Outer 2x2 complex multiply*/
void ObitMatxOuterMult2C(ObitMatx *in1, ObitMatx *in2, ObitMatx *out);
typedef void (*ObitMatxOuterMult2CFP) (ObitMatx *in1, ObitMatx *in2, ObitMatx *out);
/** 4x4 complex matrix * 4x1 complex vector multiply */
void ObitMatxVec4Mult(ObitMatx *in1, ObitMatx *in2, ObitMatx *out);
typedef void (*ObitMatxVec4MultFP) (ObitMatx *in1, ObitMatx *in2, ObitMatx *out);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitMatxClassDef.h"
} ObitMatxClassInfo;
#endif /* OBITMATX_H */ 
