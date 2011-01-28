/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitUVSVRFIXize structure      */
/*  This is intended to be included in a class structure definition   */
/* and to be used as the template for generating new classes derived  */
/* from Obit.                                                         */
/**
 * \file ObitUVRFIXizeDef.h
 * ObitUVSoln structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** UV data object to be corrected */
ObitUV *myUV;
/** Output, corrected UV data object */
ObitUV *outUV;
/** Residual UV data object from which correctsion are extracted */
ObitUV *residUV;
/** Counterrotated/averaged UV data object */
ObitUV *RFIUV;
/** Calibration object */
ObitUVSoln *SNSoln;
/** Zero fringe rate Calibration SN table */
ObitTableSN *SNTable;
/** SN Table Rows */
ObitTableSNRow *SNTableRow;
/** Number of antennas */
olong numAnt;
/** Number of baselines */
olong numBL;
/** Number of IFs */
olong numIF;
/** Number of vis in RFIUV */
olong numVis;
/** Last Vis processed in current buffer */
olong LastVisRead;
/** First vis to read in next call to NewTime */
olong NextVisRead;
/** length of visibility entry in VisApply/VisPrior/VisFollow */
olong lenVisArrayEntry;
/** current source ID in RFIUV */
olong CurSourID;
/** Prior source ID in RFIUV */
olong PriorSourID;
/** Following source ID RFIUV */
olong FollowSourID;
/** Data averaging time (days)in RFIUV */
ofloat AvgTime;
/** Prior time in RFIUV */
ofloat PriorTime;
/** Prior time per vis */
ofloat *PriorVisTime;
/** Prior vis number per vis */
olong *PriorVisNum;
/** Following Time in RFIUV */
ofloat FollowTime;
/** Following time per vis */
ofloat *FollowVisTime;
/** Following vis number per vis */
olong *FollowVisNum;
/** time of visibility in VisApply */
ofloat VisTime;
/** time of visibility in VisFollow/VisPrior */
ofloat ReadTime;
/** Apply time per vis */
ofloat *ApplyVisTime;
/** Apply vis number per vis */
olong *ApplyVisNum;
/** Visibility array for current time */
ofloat *VisApply;
/** Prior Visibility array */
ofloat *VisPrior;
/** Following Visibility array */
ofloat *VisFollow;
/** Antenna gain real per IF, per antenna */
ofloat *AntReal;
/** Antenna gain imaginary per IF, per antenna */
ofloat *AntImag;
/** Antenna group delay (sec) */
ofloat *AntDelay;
/** Antenna fringe rates (sec/sec) */
ofloat *AntRate;
/* Baseline index array */
olong *blLookup;

