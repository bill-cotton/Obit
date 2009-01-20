/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
/*  Define the basic components of the ObitUVSelfCal structure        */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitUVSelfCalDef.h
 * ObitUVSelfCal structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** Scratch UV for divided data */
ObitUV *SCData;
/** Calibration solution object */
ObitUVGSolve *mySolver;
/** Display server */
ObitDisplay *display;;
/** Sky Model for subtracting Clean model from data */
ObitSkyModel *skyModel;
/** Model calculation mode for components */
ObitSkyModelMode modelMode;
/** Number of values in flux vs BL histogram */
olong numHist;
/** increment in baseline length (lambda) in hist */
ofloat HistInc;
/** Histogram of baseline length vs Average flux density */
ofloat *hist;
/** RMSes of values in hist */
ofloat *histRMS;
/** UV range of full weight */
ofloat UVFullRange[2];
/** Total Maximum image pixel value in clean window */
ofloat totalMax;
/** Sum of flux densitiy in SkyModel */
ofloat sumCC;
/** RMS of residuals in Field 1 */
ofloat RMSFld1;
/** Total quality measure (sum flux/ RMS field 1) */
ofloat totalQual;
/** How many values in "last?" arrays */
olong numLast;
/** Last up to 5 totalQual values */
ofloat lastQual[5];
/** SN table version of last up to 5 passes */
ofloat lastSNVer[5];
