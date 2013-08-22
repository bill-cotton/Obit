/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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
#ifndef OBITVEGASUTIL_H 
#define OBITVEGASUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitIO.h"
#include "ObitImage.h"
#include "ObitOTF.h"
#include "ObitOTFArrayGeom.h"
#include "ObitOTFSkyModel.h"
#include "ObitFInterpolate.h"
#include "ObitTableCC.h"

/*-------- Obit:  Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitVEGASUtil.h
 * Utility routines for the Obit GBT/VEGAS utilities.
 *
 * \section ObitVEGASUtilparameters Control Parameters
 * The imaging control parameters are passed through the info object 
 * on the VEGAS data, these control both the output image files and the 
 * processing parameters.
 * VEGAS Data selection/calibration/editing control
 * \li  "doCalSelect" OBIT_bool (1,1,1) Select/calibrate/edit data?
 * \li  "doCalib" OBIT_int (1,1,1) >0 -> calibrate,
 * \li  "gainUse" OBIT_int (1,1,1) SN/CL table version number, 0-> use highest
 * \li  "flagVer" OBIT_int (1,1,1) Flag table version, 0-> use highest, <0-> none
 * \li  "Stokes" OBIT_string (4,1,1) Selected output Stokes parameters:
 *               "I", "V", " " -> "I"
 * \li  "BChan" OBIT_int (1,1,1) First spectral channel selected. [def all]
 * \li  "EChan" OBIT_int (1,1,1) Highest spectral channel selected. [def all]
 * \li  "Targets" OBIT_string (?,?,1) Target names selected. [def all]
 * \li  "timeRange" OBIT_float (2,1,1) Selected timerange in days. [def all]
 * \li  "Scans" OBIT_int (2,1,1) Lowest and highest selected scan numbers. [def all]
 * \li  "Feeds" OBIT_int (?,1,1) a list of selected feed numbers, [def all.]
 * 
 */

/*---------------Public functions---------------------------*/
/** Public: Average the frequencies in a GBT/VEGAS OTF. */
void ObitVEGASUtilAverage(ObitOTF *inOTF, ObitOTF *outOTF, olong chAvg, 
			  ObitErr *err);

#endif /* OBITVEGASUTIL_H */ 
