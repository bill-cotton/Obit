/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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
#ifndef OBITCONVUTIL_H 
#define OBITCONVUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"
#include "ObitFFT.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitConvUtil.h
 * ObitConvUtil utility module definition for the ObitImage class.
 *
 * Routines for image convolution
 */

/*---------------Public functions---------------------------*/
/*  Public: (de)Convolve all planes of an image with an FArray */
void ObitConvUtilConv (ObitImage *inImage, ObitFArray *convFn, 
		       gboolean doDivide, ofloat rescale,
		       ObitImage *outImage, ObitErr *err);

/*  Public: Create Gaussian array */
ObitFArray* ObitConvUtilGaus (ObitImage *inImage, ofloat Beam[3]);

/*  Public: Deconvolve two Gaussians */
void ObitConvUtilDeconv (ofloat fmaj, ofloat fmin, ofloat fpa, 
			 ofloat cmaj, ofloat cmin, ofloat cpa, 
			 ofloat *rmaj, ofloat *rmin, ofloat *rpa);

#endif /* OBITCONVUTIL_H */ 
