/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2012                                          */
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
#ifndef OBITFARRAYUTIL_H 
#define OBITFARRAYUTIL_H 

#include "ObitFArray.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFArrayUtil.h
 * ObitFArray numeric array utility module for ObitFArray class
 *
 * This class is for creating contains functions for manipulating a ObitFArrays.
 * Except as noted, magic value blanking is supported (OBIT_MAGIC).
 * There are no instances of this "class"
 * 
 */

/*---------------Public functions---------------------------*/
/** Public: Fit 2-D circular Gaussian. */
ofloat ObitFArrayUtilFitCGauss (ObitFArray *in, ofloat *FWHM, ofloat *center, 
				ofloat *peak, ObitErr *err);

/** Public: Fit 1-D Gaussian plus baseline. */
ofloat ObitFArrayUtilFit1DGauss (ObitFArray *in, ofloat *FWHM, ofloat *center, 
				 ofloat *peak, ofloat *a, ofloat *b, ObitErr *err);
/** Public: Fit multiple 1-D Gaussians plus baseline. */
ofloat ObitFArrayUtilFit1DGauss2 (ObitFArray *in, olong ngauss, ofloat *FWHM, 
				  ofloat *center,  ofloat *peak, ofloat *a, ofloat *b, 
				  ObitErr *err);
/** Public: Convolve two arrays */
ObitFArray* ObitFArrayUtilConvolve (ObitFArray *in1, ObitFArray *in2, 
				    ObitErr *err);
/** Public: Correlate two arrays */
ObitFArray* ObitFArrayUtilCorrel (ObitFArray *in1, ObitFArray *in2, 
				    ObitErr *err);
/** Public: Create Gaussian UV Taper*/
ObitFArray* ObitFArrayUtilUVGaus (olong *naxis, ofloat *cells, ofloat maprot,
				  ofloat Gaumaj, ofloat Gaumin, ofloat GauPA);
#endif /* OBITFARRAYUTIL_H */ 
