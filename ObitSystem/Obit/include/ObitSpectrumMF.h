/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2025,2026                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITSPECTRUMINTERP_H 
#define OBITSPECTRUMINTERP_H

#include "Obit.h"
#include "ObitInfoList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSpectrumMF.h
 * ObitSpectrumMF evaluates an MFImage-like tabulated spectrum
 *
 * This class is derived from the #Obit class.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitSpectrumMF Class structure. */
typedef struct {
#include "ObitSpectrumMFDef.h"   /* this class definition */
} ObitSpectrumMF;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSpectrumMF
 * returns a ObitSpectrumMF*.
 * in = object to unreference
 */
#define ObitSpectrumMFUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSpectrumMF.
 * returns a ObitSpectrumMF*.
 * in = object to reference
 */
#define ObitSpectrumMFRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSpectrumMFIsA(in) ObitIsA (in, ObitSpectrumMFGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSpectrumMFClassInit (void);

/** Public: Constructor. */
ObitSpectrumMF* newObitSpectrumMF (gchar* name);

/** Public: Constructor from value. */
ObitSpectrumMF* 
newObitSpectrumMFCreate (gchar* name, olong nFreq, odouble *Freqs,
			 odouble refFreq, olong nTerm);

/** Update info on an ObitSpectrumMF */
void ObitSpectrumMFUpdate (ObitSpectrumMF *in, olong nFreq, odouble *Freqs,
			   odouble refFreq, olong nTerm);

/** Public: ClassInfo pointer */
gconstpointer ObitSpectrumMFGetClass (void);

/** Public: Copy (deep) constructor. */
ObitSpectrumMF* ObitSpectrumMFCopy  (ObitSpectrumMF *in, ObitSpectrumMF *out, 
					     ObitErr *err);

/** Public: Copy (shallow) constructor. */
ObitSpectrumMF* ObitSpectrumMFClone (ObitSpectrumMF *in, ObitSpectrumMF *out);

/** Public: Get component flux density at a  Frequency */
ofloat ObitSpectrumMFEval (ObitSpectrumMF *in, odouble chFreq, 
			   olong sbno, ofloat *Spectrum, ofloat* spFit);

/** Public: Set channel sigmas for weighting */
void ObitSpectrumMFSetSigma (ObitSpectrumMF *in, ofloat* sigma);

/** Public: Get Spectral indices per spectral point */
void ObitSpectrumMFSI (ObitSpectrumMF *in,  ofloat *Spectrum,
		       ofloat *SI, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitSpectrumMFClassDef.h"
} ObitSpectrumMFClassInfo; 

#endif /* OBITSPECTRUMINTERP_H */ 
