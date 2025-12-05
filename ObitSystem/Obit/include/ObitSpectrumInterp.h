/* $Id$ */
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
 * \file ObitSpectrumInterp.h
 * ObitSpectrumInterp does Lagrangian interpolation of a spectrum 
 * at a given frequency
 *
 * This class is derived from the #Obit class.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitSpectrumInterp Class structure. */
typedef struct {
#include "ObitSpectrumInterpDef.h"   /* this class definition */
} ObitSpectrumInterp;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSpectrumInterp
 * returns a ObitSpectrumInterp*.
 * in = object to unreference
 */
#define ObitSpectrumInterpUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSpectrumInterp.
 * returns a ObitSpectrumInterp*.
 * in = object to reference
 */
#define ObitSpectrumInterpRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSpectrumInterpIsA(in) ObitIsA (in, ObitSpectrumInterpGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSpectrumInterpClassInit (void);

/** Public: Constructor. */
ObitSpectrumInterp* newObitSpectrumInterp (gchar* name);

/** Public: Constructor from value. */
ObitSpectrumInterp* 
newObitSpectrumInterpCreate (gchar* name, olong nFreq, odouble *Freqs,
			     ofloat *Spectrum, ofloat *chWt, olong hwidth);

/** Public: ClassInfo pointer */
gconstpointer ObitSpectrumInterpGetClass (void);

/** Public: Copy (deep) constructor. */
ObitSpectrumInterp* ObitSpectrumInterpCopy  (ObitSpectrumInterp *in, ObitSpectrumInterp *out, 
					     ObitErr *err);

/** Public: Copy (shallow) constructor. */
ObitSpectrumInterp* ObitSpectrumInterpClone (ObitSpectrumInterp *in, ObitSpectrumInterp *out);

/** Public: Replace member Spectrum */
void ObitSpectrumInterpReplace (ObitSpectrumInterp *in, ofloat *Spectrum);

/** Public: Set weights based on a  Spectrum */
void ObitSpectrumInterpSetWt (ObitSpectrumInterp *in, ofloat *Spectrum);

/** Public: Interpolate Frequency */
ofloat ObitSpectrumInterpFreq (ObitSpectrumInterp *in, odouble Freq);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitSpectrumInterpClassDef.h"
} ObitSpectrumInterpClassInfo; 

#endif /* OBITSPECTRUMINTERP_H */ 
