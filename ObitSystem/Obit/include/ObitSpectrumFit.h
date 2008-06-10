/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008                                               */
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
#ifndef OBITSPECTRUMFIT_H 
#define OBITSPECTRUMFIT_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"
#include "ObitBeamShape.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#ifdef HAVE_GSL
#include <gsl/gsl_multifit.h>
#endif /* HAVE_GSL */ 

/*-------- Obit:  Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSpectrumFit.h
 *
 * ObitSpectrumFit Class for fitting spectra to image pixels
 *
 * This class does least squares fitting of log(s) as a polynomial in log($\nu$).
 * Either an image cube or a set of single plane images at arbitrary 
 * frequencies may be fitted.
 * The result is an image cube of Log(S) with multiples of powers of log($\nu$)
 * as the planes.
 * The function ObitSpectrumFitEval will evaluate this fit and return an image
 * with the flux densities at the desired frequencies.
 * 
 * \section ObitSpectrumFitaccess Creators and Destructors
 * An ObitSpectrumFit will usually be created using ObitSpectrumFitCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitSpectrumFit should always be made using the
 * #ObitSpectrumFitRef function which updates the reference count in the object.
 * Then whenever freeing an ObitSpectrumFit or changing a pointer, the function
 * #ObitSpectrumFitUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitSpectrumFit Class structure. */
typedef struct {
#include "ObitSpectrumFitDef.h"   /* this class definition */
} ObitSpectrumFit;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSpectrumFit
 * returns a ObitSpectrumFit*.
 * in = object to unreference
 */
#define ObitSpectrumFitUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSpectrumFit.
 * returns a ObitSpectrumFit*.
 * in = object to reference
 */
#define ObitSpectrumFitRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSpectrumFitIsA(in) ObitIsA (in, ObitSpectrumFitGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSpectrumFitClassInit (void);

/** Public: Default Constructor. */
ObitSpectrumFit* newObitSpectrumFit (gchar* name);

/** Public: Create/initialize ObitSpectrumFit structures */
ObitSpectrumFit* ObitSpectrumFitCreate (gchar* name, olong nterm);
/** Typedef for definition of class pointer structure */
typedef ObitSpectrumFit* (*ObitSpectrumFitCreateFP) (gchar* name, olong nterm);

/** Public: ClassInfo pointer */
gconstpointer ObitSpectrumFitGetClass (void);

/** Public: Copy (deep) constructor. */
ObitSpectrumFit* ObitSpectrumFitCopy  (ObitSpectrumFit *in, 
				       ObitSpectrumFit *out, ObitErr *err);

/** Public: Copy structure. */
void ObitSpectrumFitClone (ObitSpectrumFit *in, ObitSpectrumFit *out, 
			   ObitErr *err);

/** Public: Fit spectrum to an image cube */
void ObitSpectrumFitCube (ObitSpectrumFit* in, ObitImage *inImage, 
			  ObitImage *outImage, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void(*ObitSpectrumFitCubeFP) (ObitSpectrumFit* in, ObitImage *inImage, 
				      ObitImage *outImage, ObitErr *err);

/** Public: Fit spectrum to an array of images */
void ObitSpectrumFitImArr (ObitSpectrumFit* in, olong nimage, ObitImage **imArr, 
			   ObitImage *outImage, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void(*ObitSpectrumFitImArrFP) (ObitSpectrumFit* in, olong nimage, ObitImage **imArr, 
			   ObitImage *outImage, ObitErr *err);


/** Public: Evaluate spectrum */
void ObitSpectrumFitEval (ObitSpectrumFit* in, ObitImage *inImage, 
			  odouble outFreq, ObitImage *outImage, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void(*ObitSpectrumFitEvalFP) (ObitSpectrumFit* in, ObitImage *inImage, 
				      odouble outFreq, ObitImage *outImage, 
				      ObitErr *err);
/** Public: Fit single spectrum */
gfloat* ObitSpectrumFitSingle (gint nfreq, olong nterm, odouble *freq, 
			       ofloat *flux, ofloat *sigma, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef gfloat*(*ObitSpectrumFitSingleFP) (gint nfreq, olong nterm, odouble *freq, 
					   ofloat *flux, ofloat *sigma, 
					   ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitSpectrumFitClassDef.h"
} ObitSpectrumFitClassInfo; 

#endif /* OBITFSPECTRUMFIT_H */ 
