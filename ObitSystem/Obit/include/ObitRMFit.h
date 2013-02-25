/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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
#ifndef OBITRMFIT_H 
#define OBITRMFIT_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"
#include "ObitBeamShape.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#ifdef HAVE_GSL
#include <gsl/gsl_multifit_nlin.h>
#endif /* HAVE_GSL */ 

/*-------- Obit:  Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitRMFit.h
 *
 * ObitRMFit Class for fitting rotatioin measures (RM) to Stokes Q and U image pixels
 *
 * This class does least squares fitting of RM and EVPA at 0 lambda^2
 * Either an image cube or a set of single plane images at arbitrary 
 * lambda squares may be fitted.
 * The result is an image cube of RM and  EVPA at 0 lambda^2 as the planes.
 * The function ObitRMFitEval will evaluate this fit and return an image
 * with the flux densities at the desired frequencies.
 * 
 * \section ObitRMFitaccess Creators and Destructors
 * An ObitRMFit will usually be created using ObitRMFitCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitRMFit should always be made using the
 * #ObitRMFitRef function which updates the reference count in the object.
 * Then whenever freeing an ObitRMFit or changing a pointer, the function
 * #ObitRMFitUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitRMFit Class structure. */
typedef struct {
#include "ObitRMFitDef.h"   /* this class definition */
} ObitRMFit;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitRMFit
 * returns a ObitRMFit*.
 * in = object to unreference
 */
#define ObitRMFitUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitRMFit.
 * returns a ObitRMFit*.
 * in = object to reference
 */
#define ObitRMFitRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitRMFitIsA(in) ObitIsA (in, ObitRMFitGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitRMFitClassInit (void);

/** Public: Default Constructor. */
ObitRMFit* newObitRMFit (gchar* name);

/** Public: Create/initialize ObitRMFit structures */
ObitRMFit* ObitRMFitCreate (gchar* name, olong nterm);
/** Typedef for definition of class pointer structure */
typedef ObitRMFit* (*ObitRMFitCreateFP) (gchar* name, 
						     olong nterm);

/** Public: ClassInfo pointer */
gconstpointer ObitRMFitGetClass (void);

/** Public: Copy (deep) constructor. */
ObitRMFit* ObitRMFitCopy  (ObitRMFit *in, ObitRMFit *out, ObitErr *err);

/** Public: Copy structure. */
void ObitRMFitClone (ObitRMFit *in, ObitRMFit *out, ObitErr *err);

/** Public: Fit spectrum to an image cube */
void ObitRMFitCube (ObitRMFit* in, ObitImage *inQImage, ObitImage *inUImage, 
		    ObitImage *outImage, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void(*ObitRMFitCubeFP) (ObitRMFit* in, 
				ObitImage *inQImage, ObitImage *inUImage, 
				ObitImage *outImage, ObitErr *err);

/** Public: Fit spectrum to an array of images */
void ObitRMFitImArr (ObitRMFit* in, olong nimage, 
		     ObitImage **imQArr, ObitImage **imUArr, 
		     ObitImage *outImage, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void(*ObitRMFitImArrFP) (ObitRMFit* in, olong nimage, 
				 ObitImage **imQArr, ObitImage **imUArr, 
				 ObitImage *outImage, ObitErr *err);


/** Public: Fit single spectrum */
ofloat* ObitRMFitSingle (olong nlamb2, olong nterm, odouble refLamb2, odouble *lamb2, 
			 ofloat *qflux, ofloat *qsigma, ofloat *uflux, ofloat *usigma, 
			 ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ofloat*(*ObitRMFitSingleFP) (olong nlamb2, olong nterm, odouble refLamb2, odouble *lamb2,
				     ofloat *qflux, ofloat *qsigma, ofloat *uflux, ofloat *usigma,
				     ObitErr *err);

/** Public: Make fitting arg structure */
gpointer ObitRMFitMakeArg (olong nlamb2, olong nterm, 
			   odouble refLamb2, odouble *lamb2, 
			   ofloat **out, ObitErr *err);

/** Public: Fit single RM using arg */
void ObitRMFitSingleArg (gpointer arg,
			 ofloat *qflux, ofloat *qsigma, 
			 ofloat *uflux, ofloat *usigma,
			 ofloat *out);

/** Public: Kill fitting arg structure */
void ObitRMFitKillArg (gpointer arg);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitRMFitClassDef.h"
} ObitRMFitClassInfo; 

#endif /* OBITFRMFIT_H */ 
