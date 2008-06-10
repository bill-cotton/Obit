/* $Id: ObitUVWeight.h,v 1.8 2007/08/31 17:24:49 bcotton Exp $  */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITUVWEIGHT_H 
#define OBITUVWEIGHT_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitUV.h"
#include "ObitFArray.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVWeight.h
 * ObitUVWeight uv data class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This class is for weighting uv data in preparation for imaging.
 * Uniform, Natural and Briggs Robust weighting are supported.
 * No instances of this class should be needed outside of the 
 * ObitUVWeightData function.
 * 
 * ObitUVWeightData convolves UV weights onto a grid to determine weighting function.
 * Then the data weights are modified by the weighting function.
 * The control parameters are attached to the ObitInfoList member info
 * on uvdata.  The parameters are:
 * \li "nuGrid" OBIT_long scalar = Number of "U" pixels in weighting grid.
 *              [defaults to "nx"]
 * \li "nvGrid" OBIT_long scalar = Number of "V" pixels in weighting grid.
 * \li "WtBox"  OBIT_long scalar = Size of weighting box in cells [def 0]
 * \li "WtFunc" OBIT_long scalar = Weighting convolution function [def. 1]
 *              1=Pill box, 2=linear, 3=exponential, 4=Gaussian
 *              if positive, function is of radius, negative in u and v.
 * \li "xCells" OBIT_float scalar = Image cell spacing in X in asec.
 * \li "yCells" OBIT_float scalar = Image cell spacing in Y in asec.
 * \li "UVTaper"OBIT_float scalar = UV taper width in kilowavelengths. [def. no taper].
 *              NB: If the taper is applied her is should not also be applied
 *              in the imaging step as the taper will be applied to the
 *              output data.
 * \li "Robust" OBIT_float scalar = Briggs robust parameter. [def. 0.0]
 *              < -7 -> Pure Uniform weight, >7 -> Pure natural weight.
 *              Uses AIPS rather than Briggs definition of Robust.
 * \li "WtPower" OBIT_float scalar = Power to raise weights to.  [def = 1.0]
 *              Note: a power of 0.0 sets all the output weights to 1 as modified
 *              by uniform/Tapering weighting.  Applied in determing weights 
 *              as well as after.
 * \li "MaxBaseline" OBIT_float scalar = maximum baseline length in wavelengths.
 *              Default = 1.0e15.  Output data not flagged by this criteria.
 * \li "MinBaseline" OBIT_float scalar = minimum baseline length in wavelengths.
 *              Default = 1.0e15.Output data not flagged by this criteria.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitUVWeight Class structure. */
typedef struct {
#include "ObitUVWeightDef.h"   /* this class definition */
} ObitUVWeight;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVWeightClassInit (void);

/** Public: Constructor. */
ObitUVWeight* newObitUVWeight (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitUVWeightGetClass (void);

/** Public: Determine weighting and correct an ObitUV. */
void ObitUVWeightData (ObitUV *uvdata, ObitErr *err);
typedef void (*ObitUVWeightDataFP) (ObitUV *uvdata, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVWeightClassDef.h"
} ObitUVWeightClassInfo; 

#endif /* OBITUVWEIGHT_H */ 
