/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2025                                               */
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
#ifndef OBITFARASYN_H 
#define OBITFARASYN_H 

#include "ObitCArray.h"
#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"
#include "ObitUtil.h"
#include "ObitFArray.h"
#include "ObitThread.h"
#include "ObitSinCos.h"
#include "ObitComplex.h"
#include "ObitSystem.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif /* VELIGHT */
#include <math.h>
void sincosf(float x, float *sin, float *cos); /* Shouldn't need */
#if HAVE_GSL==1  /* GSL stuff */
#include "gsl/gsl_blas.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_multifit_nlin.h"
#endif /* GSL stuff */

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFaraSyn.h
 *
 * ObitFaraSyn Class for Faraday Synthesis/Anaslysis functions
 *
 * Class documentation
 * This class is for Faraday Synthesis and analysis of Stokes Q and U
 * cubes and can be used to generate image cubes with Faraday depth as
 * the third axis or search/fit for peaks in the Faraday depth in each
 * pixel. 
 * 
 * A linearly polarized wave passing through a magnetized plasma will
 * experience a Faraday rotation of the angle of the polarized signal
 * \cite{Brentjens2005}:
 * \begin{equation}\label{rot}
 * \Delta \chi\ =\ \lambda^2\ 0.81 \int n_e B_\parallel dr,
 * \end{equation}
 * where $\lambda$ is the wavelength in m, $n_e$ is the electron density in
 * cm$^{-3}$, $B_\parallel$is the strength of the component of the
 * magnetic field along the line of sight in $\mu$Gauss and r is distance
 * in parsec.
 * from Obit Memo 76
 * 
 * \section ObitFaraSynaccess Creators and Destructors
 * An ObitFaraSyn will usually be created using ObitFaraSynCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitFaraSyn should always be made using the
 * #ObitFaraSynRef function which updates the reference count in the object.
 * Then whenever freeing an ObitFaraSyn or changing a pointer, the function
 * #ObitFaraSynUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitFaraSyn Class structure. */
typedef struct {
#include "ObitFaraSynDef.h"   /* this class definition */
} ObitFaraSyn;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitFaraSyn
 * returns a ObitFaraSyn*.
 * in = object to unreference
 */
#define ObitFaraSynUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitFaraSyn.
 * returns a ObitFaraSyn*.
 * in = object to reference
 */
#define ObitFaraSynRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitFaraSynIsA(in) ObitIsA (in, ObitFaraSynGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitFaraSynClassInit (void);

/** Public: Default Constructor. */
ObitFaraSyn* newObitFaraSyn (gchar* name);

/** Public: Create/initialize ObitFaraSyn structures */
ObitFaraSyn* ObitFaraSynCreate 
  (gchar* name,	
   olong nLamb2, odouble *lamb2, odouble refLamb2,
   ObitFArray **inQFArrays, ObitFArray **inUFArrays,
   ofloat *plnWt, ofloat *QRMS, ofloat *URMS,
   olong nOut, ofloat minRMSyn, ofloat maxRMSyn, ofloat delRMSyn,
   ObitCArray **RMPlanes, 
   gboolean doWrite, ObitImage* outAImage, ObitImage* outPImage,
   ObitImage *outFitRM, ofloat minQUSNR, ofloat minFrac, 
   gboolean doError, gboolean doRMSyn, ofloat maxChi2);
/** Typedef for definition of class pointer structure */
typedef ObitFaraSyn* (*ObitFaraSynCreateFP)
(gchar* name,	
 olong nLamb2, odouble *lamb2, odouble refLamb2,
 ObitFArray **inQFArrays, ObitFArray **inUFArrays,
 ofloat *plnWt, ofloat *QRMS, ofloat *URMS,
 olong nOut, ofloat minRMSyn, ofloat maxRMSyn, ofloat delRMSyn,
 ObitCArray **RMPlanes,
 gboolean doWrite, ObitImage* outAImage, ObitImage* outPImage,
 ObitImage *outFitRM, ofloat minQUSNR, ofloat minFrac, 
 gboolean doError, gboolean doRMSyn, ofloat maxChi2);

/** Public: Set Q & U data */
void ObitFaraSynSetQU (ObitFaraSyn* in, 
		       olong nlamb2, odouble *lamb2, odouble refLamb2,
		       ObitFArray **inQFArrays, ObitFArray **inUFArrays,
		       ofloat *plnWt, ofloat *QRMS, ofloat *URMS);

/** Public: Set RM Synthesis request/output */
void ObitFaraSynSetRMSyn 
  (ObitFaraSyn* in,
   olong nOut, ofloat minRMSyn, ofloat delRMSyn,
   ObitCArray **RMPlanes,
   gboolean doWrite, ObitImage* outAImage, ObitImage* outPImage);

/** Public: Set RM Analysis request/output */
void ObitFaraSynSetRMAna 
  (ObitFaraSyn* in,
   ofloat minRMSyn, float maxRMSyn, ofloat delRMSyn,
   ObitImage* outFitRM, ofloat minQUSNR, ofloat minFrac, 
   gboolean doError, gboolean doRMSyn, ofloat maxChi2);

/** Public: ClassInfo pointer */
gconstpointer ObitFaraSynGetClass (void);

/** Public: Copy (deep) constructor. */
ObitFaraSyn* ObitFaraSynCopy  (ObitFaraSyn *in, ObitFaraSyn *out, ObitErr *err);

/** Public: Copy structure. */
void ObitFaraSynClone (ObitFaraSyn *in, ObitFaraSyn *out, ObitErr *err);

/* Public: Faraday Synthesis. Faraday depth cube */
void ObitFaraSynRMSyn (ObitFaraSyn *in, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitFaraSynRMSynFP) (ObitFaraSyn *in, ObitErr *err);

/* Public: Faraday Analysis. Peak Faraday depth cube */
void ObitFaraSynRMAna (ObitFaraSyn *in, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitFaraSynRMAnaFP) (ObitFaraSyn *in, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitFaraSynClassDef.h"
} ObitFaraSynClassInfo; 

#endif /* OBITFFaraSyn_H */ 
