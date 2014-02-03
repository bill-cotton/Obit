/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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
 * Utility routine for fast sine/cosine calculation 
 * This routine uses a table lookup followed by a 1 term Taylor's 
 * series expansion.
 * This produces moderate accuracy but fast calculation of sine/cosine pairs.
 * Speed is ~3 times faster than standard sin/cos routines
 * Comparison accuracy in 10^8 trials over range of angles gives:
 * Avg difference 3.42813e-08, rms 1.51385e-06, max. difference 4.83e-6.
 * The rms difference corresponds to an angle error of 8.7e-5 deg.
 * This utility is useful for calculating instrumental responses to models 
 * or other applications in which errors do not seriously accumulate.
 */
#include "Obit.h"

#ifndef OBITSINCOS_H 
#define OBITSINCOS_H 
/** Init sine/cosine  */
void ObitSinCosInit();
/** Calculate sine/cosine of angle */
void ObitSinCosCalc(ofloat angle, ofloat *sin, ofloat *cos);
/** Calculate sine/cosine of vector of angles */
void ObitSinCosVec(olong n, ofloat *angle, ofloat *sin, ofloat *cos);
#endif /* OBITSINCOS_H */ 
