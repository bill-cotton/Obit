/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
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
#ifndef OBITFEATHERUTIL_H 
#define OBITFEATHERUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"
#include "ObitFFT.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFeatherUtil.h
 * Image feathering utility routine definition.
 *
 * Routines for feathering together images at different resolutions.
 */

/*---------------Public functions---------------------------*/
/**  Public: Create an FFT object suitable for FFTing an image */
ObitFFT* ObitFeatherUtilCreateFFT (ObitImage *in, ObitFFTdir dir);

/**  Public:  Return scratch Image the size needed for a padding */
ObitImage* ObitFeatherUtilCreatePadImage(ObitFFT *inFFT, 
					 ObitImage *inImage, ObitErr *err);
 
/**  Public:  Create a half plane CArray suitable for the output of 
     FFTing an image */
ObitCArray* ObitFeatherUtilCreateFFTArray(ObitFFT *inFFT);
 
/**  Public: Zero Pads an image as needed for an FFT */
void ObitFeatherUtilPad (ObitFFT *inFFT, ObitImage *inImage, 
			 ObitImage *outImage, ObitErr *err);

/**  Public: Increases the size of an image and zero pads */
void ObitFeatherUtilBigger (olong *naxis, ObitImage *inImage, 
			    ObitImage *outImage, ObitErr *err);

/**  Public:Zero Pads an array as needed for an FFT PadArray */
void ObitFeatherUtilPadArray (ObitFFT *inFFT, ObitFArray *inArray, 
			      ObitFArray *outArray);

/**  Public:  Extract a Real array from one padded for FFTs */
ObitFArray* ObitFeatherUtilExtract (ObitFFT *inFFT, ObitFArray *inArray, 
				    ObitErr *err);

/**  Public:Make uv plane weighting array  */
ObitFArray* ObitFeatherUtilMakeBeamMask (ObitImage *inImage, ObitFFT *inFFT, 
				  ObitErr *err);

/**  Public: Fill an FArray with a model the size and shape of the 
     resolution in an image */
void ObitFeatherUtilCreateModel (ObitImage *image, ObitFArray *outArray);

/**  Public: Accumulate the weighted FT of an FArray */
void ObitFeatherUtilAccumImage (ObitFFT *FFTfor, ObitImage *inImage, 
				ObitFArray *wtArray, ObitCArray *accArray, 
				ObitCArray *workArray, ObitErr *err);

/**  Public: HGEOM-like operation (Before EWG got to it) */
void ObitFeatherUtilInterpol (ObitImage *inImage, ObitImage *tmplImage, 
			      ObitImage *outImage, ObitErr *err);

/**  Public: Extract the subimage in inAray corresponding to outImage */
void ObitFeatherUtilSubImage (ObitImage *inImage, ObitFArray *inArray,
			      ObitImage *outImage, ObitErr *err);


#endif /* OBITFEATHERUTIL_H */ 
