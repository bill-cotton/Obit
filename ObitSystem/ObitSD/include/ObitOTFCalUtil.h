/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2013                                          */
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
#ifndef OBITOTFCALUTIL_H 
#define OBITOTFCALUTIL_H 

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include <glibconfig.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitIO.h"
#include "ObitImage.h"
#include "ObitOTF.h"
#include "ObitOTFArrayGeom.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFCalUtil.h
 * Calibration utility routines for the OTF class.
 *
 */

/*---------------Public functions---------------------------*/
/** Public: Fit calibrator scans */
void ObitOTFCalUtilFitCal (ObitOTF *inOTF, olong detect, ObitErr *err);

/** Public: Fit calibrator on/off scans */
void ObitOTFCalUtilFitOnOff (ObitOTF *inOTF, olong detect, ObitErr *err);

/** Public: Fit bandpass cal from on/off scans */
void ObitOTFCalUtilFitBPOnOff (ObitOTF *inOTF, olong offScan, olong onScan, 
			       olong BPVer, ObitErr *err);

/** Public: Fit calibrator nod scan */
void ObitOTFCalUtilFitNod (ObitOTF *inOTF, olong detect, ObitErr *err);

/** Public: Fit tipping scan */
void ObitOTFCalUtilFitTip (ObitOTF *inOTF, ObitErr *err);

/** Public: Add flagging entry to flag table  */
ObitIOCode ObitOTFCalUtilFlag (ObitOTF *inOTF, ObitErr *err);

#endif /* OBITOTFCALUTIL_H */ 
