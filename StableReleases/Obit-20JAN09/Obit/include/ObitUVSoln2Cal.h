/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITUVSOLN2CAL_H 
#define OBITUVSOLN2CAL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitUV.h"
#include "ObitUVSoln.h"
#include "ObitTableSN.h"
#include "ObitTableCL.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVSoln2Cal.h
 * Routines to Apply a Soln (SN) table to a Cal (CL) table writing a 
 * new Cal table for UV class. 
 *
 * The following options can be entered onto the inUV uv data info list:
 * \li "interMode", OBIT_string (4,1,1)  Interpolation mode:
                    default or blank = "2PT "
 * \li "2PT " = linear vector interpolation with no SN smoothing.
 * \li "SELF" = Use only SN solution from same source which is closest in time.
 * \li "POLY" = Fit a polynomial to the SN rates and delays.
 *             Use the integral of the rate polynomial for the phases.
 * \li "SIMP" = Simple linear phase connection between SN phase
 *              entries, assumes phase difference less than 180 degrees.
 * \li "AMBG" = Linear phase connection using rates to resolve phase ambiguities.
 * \li "CUBE" = As AMBG but fit third order polynomial to phases and rates.
 * \li "MWF " = Median window filter of SN table before 2PT interpolation
 * \li "GAUS" = Gaussian smoothing of SN table before 2PT interpolation,
 * \li "BOX " = Boxcar smoothing of SN table before 2PT interpolation,
 *              boxcar width set by adverb INTPARM.
 * \li "interParm", OBIT_float (3,1,1) interpolation parameters
 *                 default = 0's
 * \li mode="BOX ", smoothing time in hours for amplitude, phase, delay/rate
 * \li mode="MWF ", window size in hours for amplitude, phase, delay/rate
 * 
 * \li "interNPoly", OBIT_int (1,1,1) number of terms in polynomial for
 *                   mode="POLY", default = 2
 *
 * \li "maxInter", OBIT_float (1,1,1) Max. time (min) over which to interpolate.
 *                 default = 1440.0;
 *
 * \li "allPass", OBIT_bool (1,1,1) If true copy unmodified entries as well.
 *                else only data calibrated.  Default = FALSE.
 *
 * \li "solnVer", OBIT_int (1,1,1) Solution (SN) table version to use.
 *                0=>highest, default 0;
 * \li "calIn",   OBIT_int (1,1,1) Input calibration (CL) table version to use.
 *                0=>highest, -1 => just convert solnVer to a CL table, default 0;
 * \li "calOut",  OBIT_int (1,1,1) Output calibration (CL) table version to use.
 *                0=>create new, default 0;
 * \li "subA",    OBIT_int (1,1,1) Subarray, default 1
 * \li "refAnt",  OBIT_int (1,1,1) Reference antenna, default 1
 *
 * Any data selection parameters on the input UV data info object will 
 * be applied.
 */

/*---------------Public functions---------------------------*/
/** Public: Apply a Soln table to a Cal table writing a new Cal table */
ObitTableCL* ObitUVSoln2Cal (ObitUV *in, ObitUV *out, ObitErr *err);

#endif /* OBITUVSOLN2CAL_H */ 
