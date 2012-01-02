/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2011                                               */
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
#ifndef OBITUVRLDELAY_H 
#define OBITUVRLDELAY_H 

#include "ObitUV.h"
#include "ObitTableSN.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVRLDelay.h
 * ObitUVRLDelay utility for cross polarized Delay/phase calibration
 *
 *
 * Routines determine calibration for an ObitUV and write an SN table.
 * Control parameters are on the info member.
 * \li "subA"    OBIT_int   (1,1,1) Selected subarray (default 1)
 * \li "minSNR"  OBIT_float (1,1,1) Minimum acceptable SNR (5)
 * \li "minNo"   OBIT_int   (1,1,1) Min. no. antennas. (default 4)
 * \li "antWt"   OBIT_float (*,1,1) Antenna weights. (default 1.0)
 * \li "UVR_Full"OBIT_float (2,1,1) Range of baseline lengths with full weight
 *                                  (kilolamda). If none is given then 
 *                                  derive one if possible.
 * \li "WtUV"    OBIT_float (1,1,1) Weight outside of UVRANG. (default 1.0)
 * \li "RLPhase" OBIT_float (1,1,1) Desired phase after correction (default 0.0)
 *                                  in deg at reference frequency
 * \li "RM"      OBIT_float (1,1,1) Rotation measure (rad/m^2) (default 0.0)
 * \li "refAnt"  OBIT_long  (1,1,1) reference antenna for SN table (default 1)
 * \li "prtLv"   OBIT_int   (1,1,1) Print level (default no print)
 * 
 */

/*---------------Private structures----------------*/
/*---------------Public functions---------------------------*/
/** Determine cross pol calibration from UV data divided by model. */
ObitTableSN* ObitUVRLDelayCal (ObitUV *inUV, ObitUV *outUV, ObitErr *err);

#endif /* OBITFUVRLDELAY_H */ 
