/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2009                                          */
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
/*  Define the basic components of the ObitUV ClassInfo structure     */
/* This is intended to be included in a classInfo structure definition*/
#include "ObitDataClassDef.h"  /* Parent class ClassInfo definition file */
/** Function pointer to scratch copy constructor. */
newObitUVScratchFP newObitUVScratch;
/** Function pointer to are two UVs the same. */
ObitUVSameFP ObitUVSame;
/** Function pointer to Read. */
ObitUVReadFP ObitUVRead;
/** Function pointer to Read multi buffer. */
ObitUVReadMultiFP ObitUVReadMulti;
/** Function pointer to Reread multi buffer. */
ObitUVReReadMultiFP ObitUVReReadMulti;
/** Function pointer to Read with selection/calibration. */
ObitUVReadSelectFP ObitUVReadSelect;
/** Function pointer to Read with selection/calibration multi buffer. */
ObitUVReadMultiSelectFP ObitUVReadMultiSelect;
/** Function pointer to Reread with selection/calibration multi buffer. */
ObitUVReReadMultiSelectFP ObitUVReReadMultiSelect;
/** Function pointer to Write. */
ObitUVWriteFP ObitUVWrite;
/** Function pointer to Create an associated table. */
newObitUVTableFP newObitUVTable;
/** Destroy an associated Table(s) */
ObitUVZapTableFP ObitUVZapTable;
/** Fully instantiate. */
ObitUVFullInstantiateFP ObitUVFullInstantiate;
/** Copy associated Table(s) */
ObitUVCopyTablesFP ObitUVCopyTables;
/** Update disk resident tables information */
ObitUVUpdateTablesFP ObitUVUpdateTables;
/** Channel selection in FG table */
ObitUVChanSelFP  ObitUVChanSel;
