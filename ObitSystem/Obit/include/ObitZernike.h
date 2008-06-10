/* $Id: ObitZernike.h,v 1.2 2007/08/31 17:24:49 bcotton Exp $  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
#ifndef OBITZERNIKE_H 
#define OBITZERNIKE_H 

#include "Obit.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitZernike.h
 * ObitZernike utility module definition.
 *
 * This utility package supplies Zernike terms and derivatives.
 * Zernike polynomials are used to describe phase errors across an aperature.
 */


/*---------------Public functions---------------------------*/
/** Return Zernike term N for X and Y */
ofloat ObitZernike (gint n, ofloat x, ofloat y);

/** Return Zernike term N gradient in X for X and Y */
ofloat ObitZernikeGradX (gint n, ofloat x, ofloat y);

/** Return Zernike term N gradient in Y for X and Y */
ofloat ObitZernikeGradY (gint n, ofloat x, ofloat y);

/** Return Zernike term N for Polar coordinates rho and phi */
ofloat ObitZernikePolar (gint n, ofloat rho, ofloat phi);

#endif /* ZERNIKE_H */ 
