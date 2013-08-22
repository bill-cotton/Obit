/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2013                                          */
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
/*  Define the basic components of the ObitOTFSel structure            */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitOTFSelDef.h
 * ObitOTFSel structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** File Type requested (FITS) */
ObitIOType FileType;
/** Max. number of records per read/write */
olong nRecPIO;
/** Number of records for next read */
olong numRecRead;
/** Which version of the calibrator table */
olong calVersion;
/** Calibrate/edit/select Data? */
gboolean doCalSelect;
/** Calibrate Data? */
gboolean doCal;
/** Bandpass calibrate Data? */
gboolean doBP;
/** Index Data? */
gboolean doIndex;
/** Translate Stokes? */
gboolean transPol;
/** Index table (as Obit*) */
Obit *IndexTable;
/** Index Table Row (as Obit*) */
Obit *IndexTableRow;
gboolean doFlag;
/** version number of FG table for flagging */
olong FGversion;
/** version number of OTFBP table for bandpass */
olong BPVer;
/** Number of rows in flag table */
olong numRow;
/** Last Row read */
olong LastRowRead;
/** First record in scan */
olong scanFirstRec;
/** Last record in scan */
olong scanLastRec;
/** Number of Stokes parameters in output */
olong numberPoln;
/** Increment in data array between Stokes */
olong jincs;
/** Start channel (1-rel) */
olong startChann;
/** Number of channels */
olong numberChann;
/** Increment in data array between Frequencies */
olong jincf;
/** Number of Feeds */
olong numberFeed;
/** Increment in data array between feeds */
olong jincfeed;
/** Number of entries in feeds 0=> all selected. */
olong numberFeedList;
/** List of selected feeds, NULL => all selected */
olong *feeds;
/** Number of entries in targets, 0=> all selected. */
olong numberTargetList;
/** List of selected target ids, NULL => all selected */
olong *targets;
/** Lowest and highest scan numbers */
olong scans[2];
/** Start and end times in days */
ofloat timeRange[2];
/** Selected Stokes parameter(s) */
gchar Stokes[5];
/** Keep cal-on data? */
gboolean keepCal;
/** Replace data with cal value */
gboolean replCal;
