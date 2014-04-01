/* $Id$ 
hacked version for test*/
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2011-2014                                          */
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
 * Utility routine for fast exp(x) calculation 
 * This routine uses a MacLauren's series expansion.
 * This produces moderate accuracy but fast calculation of exp(x), 
 * This utility is useful for calculating instrumental responses to models 
 * or other applications in which errors do not seriously accumulate.
 */
#include "Obit.h"

#ifndef OBITEXP_H 
#define OBITEXP_H 
/** Calculate exponential of -arg */
ofloat ObitExpCalc(ofloat arg);
/** Calculate exp of negatives of a vector */
void ObitExpVec(olong n, ofloat *argarr, ofloat *exparr);
#endif /* OBITEXP_H */ 
