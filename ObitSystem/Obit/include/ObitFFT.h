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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITFFT_H 
#define OBITFFT_H 

/* FFTW3 ?  */
#if HAVE_FFTW3==1
#include <fftw3.h>
#elif HAVE_FFTW==1  /* else Try FFTW 2 */
#include <fftw.h>
#include <rfftw.h>
#elif HAVE_GSL==1  /* Else try GSL version */
#include <gsl/gsl_fft_complex_float.h>
#include <gsl/gsl_fft_halfcomplex_float.h>
#include <gsl/gsl_fft_real_float.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_real.h>
#endif /* HAVE_GSL */
#include "Obit.h"
#include "ObitThread.h"
#include "ObitFArray.h"
#include "ObitCArray.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFFT.h
 * ObitFFT Fast Fourier Transform class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This class is for performing FFT on memory resident data.
 * This implementation uses the FFTW package if available, else gsl
 *
 * \section FFTOrder  Data order
 * Data passed to and from ObitFFT routines are as an ObitFArray or ObitCArray
 * which are stored in column major (Fortran) order.
 * Data are passed and returned in "center-at-the-edge" (i.e. unnatural) order
 * and there is NO transpose of the array axes.
 * In the half complex form, the first axis is nx/2+1 where nx is the number of 
 * real elements.
 * Only even numbers of elements on each axis will work well.
 * 
 * \section ObitFFTaccess Creators and Destructors
 * An ObitFFT can be created using newObitFFT which allows specifying 
 * a name for the object, and the type, size and direction of the transform.
 *
 * A copy of a pointer to an ObitFFT should always be made using the
 * #ObitFFTRef function which updates the reference count in the object.
 * Then whenever freeing an ObitFFT or changing a pointer, the function
 * #ObitFFTUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum ObitFFTdir
 * enum for FFT direction
 * This specifies the status of the connection a disk resident data.
 */
enum obitFFTdir {
  /** Sign of exponent in transform = -1, real to complex */
  OBIT_FFT_Forward,
  /** Sign of exponent in transform = +1, complex to real  */
  OBIT_FFT_Reverse
}; /* end enum obitFFTdir */
/** typedef for enum for ObitIO object status. */
typedef enum obitFFTdir ObitFFTdir;

/**
 * \enum ObitFFTtype
 * enum for FFT type, complex to complex or real to complex/
 * This specifies the status of the connection a disk resident data.
 */
enum obitFFTtype {
  /** Full complex to complex transforms */
  OBIT_FFT_FullComplex,
  /** Real to half complex or reverse.  */
  OBIT_FFT_HalfComplex
}; /* end enum obitFFTtype */
/** typedef for enum for ObitIO object status. */
typedef enum obitFFTtype ObitFFTtype;

/*--------------Class definitions-------------------------------------*/
/** ObitFFT Class structure. */
typedef struct {
#include "ObitFFTDef.h"   /* this class definition */
} ObitFFT;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitFFT
 * returns a ObitFFT*.
 * in = object to unreference
 */
#define ObitFFTUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitFFT.
 * returns a ObitFFT*.
 * in = object to reference
 */
#define ObitFFTRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitFFTIsA(in) ObitIsA (in, ObitFFTGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitFFTClassInit (void);

/** Public: Constructor. */
ObitFFT* newObitFFT (gchar* name, ObitFFTdir dir, ObitFFTtype type, 
		     olong rank, olong *dim);
/** Typedef for definition of class pointer structure */
typedef ObitFFT* (*newObitFFTFP) (gchar *name, ObitFFTdir dir, 
				  ObitFFTtype type, 
				  olong rank, olong *dim);

/** Public: ClassInfo pointer */
gconstpointer ObitFFTGetClass (void);

/** Public: Suggest efficient size for a transform */
olong ObitFFTSuggestSize (gint length);

/** Public: Real to half Complex. */
void ObitFFTR2C (ObitFFT *in, ObitFArray *inArray, ObitCArray *outArray);
typedef void (*ObitFFTR2CFP) (ObitFFT *in, ObitFArray *inArray, 
			    ObitCArray *outArray);

/** Public: Half Complex to Real. */
void ObitFFTC2R (ObitFFT *in, ObitCArray *inArray, ObitFArray *outArray);
typedef void (*ObitFFTC2RFP) (ObitFFT *in, ObitCArray *inArray, 
			    ObitFArray *outArray);

/** Public: Full Complex to Complex. */
void ObitFFTC2C (ObitFFT *in, ObitCArray *inArray, ObitCArray *outArray);
typedef void (*ObitFFTC2CFP) (ObitFFT *in, ObitCArray *inArray, 
			    ObitCArray *outArray);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitFFTClassDef.h"
} ObitFFTClassInfo; 

#endif /* OBITFFT_H */ 
