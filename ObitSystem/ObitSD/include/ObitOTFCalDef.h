/* $Id: ObitOTFCalDef.h,v 1.5 2008/02/13 21:13:13 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
/* Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitOTFCal structure           */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitOTFCalDef.h
 * ObitOTFCal structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** I/O status */
ObitIOStatus myStatus;
/** ObitOTF data Descriptor */
ObitOTFDesc* myDesc;
/** ObitOTF data Selector */
ObitOTFSel* mySel;
/** data flagging structure */
ObitOTFCalFlagS *flag;
/* Flag table */
Obit *FlagTable;
/** Solution table (as Obit*) */
Obit *SolnTable;
/** Solution Table Row (as Obit*)*/
Obit *SolnTableRow;
/** Calibration  table (as Obit*) */
Obit *CalTable;
/** Calibration Table Row (as Obit*)*/
Obit *CalTableRow;
/** Last Row read */
olong LastRowRead;
/** Number of rows in table */
olong numRow;
/** Use Solution (else Calibration) table? */
gboolean doSolnTable;
/** Apply flagging? */
gboolean doFlag;
/** number of detectors (all feeds, poln, frequencies)*/
olong numDet;
/** number of polynomial terms */
olong numPoly;
/** Number of Feeds */
olong numFeed;
/** Number of Spectral channels */
olong numChan;
/** Start channel number (1-rel) */
olong bChan;
/** highest channel (1-rel) */
olong eChan;
/** Number of Stokes in data*/
olong numStok;
/** Stokes conversion type, 0=>none, 1=I, 2=V */
olong PolMode;
/** current Target ID in cal table */
olong CurTargID;
/** Prior Target ID in cal table */
olong PriorTargID;
/** Following Target ID in cal table */
olong FollowTargID;
/** Time if calibration */
ofloat CalTime;
/** Prior time in cal table */
ofloat PriorCalTime;
/** Following Time in cal table */
ofloat FollowCalTime;
/** RA, Dec offset to pointing */
ofloat CalApplyAzoff, CalApplyEloff;
/** Prior RA, Dec offset to pointing */
ofloat CalPriorAzoff, CalPriorEloff;
/** Following RA, Dec offset to pointing */
ofloat CalFollowAzoff, CalFollowEloff;
/** cal value to subtract if "cal on" (ndetect) */
ofloat *CalApplyCal;
/** additive term to apply (ndetect) */
ofloat *CalApplyAdd;
/** Multiplicative term to apply (ndetect)*/
ofloat *CalApplyMult;
/** Weight factor to apply (ndetect)*/
ofloat *CalApplyWt;
/** polynomial sky brightness to apply (npoly) */
ofloat *CalApplyPoly;
/** Prior cal value to subtract if "cal on" (ndetect) */
ofloat *CalPriorCal;
/** Prior additive term to apply (ndetect)*/
ofloat *CalPriorAdd;
/** Prior Multiplicative term to apply (ndetect)*/
ofloat *CalPriorMult;
/** Prior weight factor to apply (ndetect)*/
ofloat *CalPriorWt;
/** Prior polynomial sky brightness to apply (npoly)*/
ofloat *CalPriorPoly;
/** Following  cal value to subtract if "cal on" (ndetect) */
ofloat *CalFollowCal;
/** Following additive term to apply *(ndetect)*/
ofloat *CalFollowAdd;
/** Following Multiplicative term to apply (ndetect)*/
ofloat *CalFollowMult;
/** Following weight to apply (ndetect)*/
ofloat *CalFollowWt;
/** Following polynomial sky brightness to apply (npoly)*/
ofloat *CalFollowPoly;
/** Legendre polynomial terms, one for each coef and each detector
    index = (detector # -1)*ncoef + icoef
    az = "X", el = "Y" */
ofloat *poly;
/** Feed wanted detectors */
gboolean *WantDetect;
