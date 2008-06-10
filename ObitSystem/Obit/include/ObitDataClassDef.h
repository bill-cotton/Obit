/* $Id$  */
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
/*  Define the basic components of the ObitData ClassInfo structure     */
/* This is intended to be included in a classInfo structure definition*/
#include "ObitClassDef.h"  /* Parent class ClassInfo definition file */
/** Function pointer to scratch copy constructor. */
newObitDataScratchFP newObitDataScratch;
/** Function pointer to are two ObitDatas the same. */
ObitDataSameFP ObitDataSame;
/** Function pointer to Copy. */
ObitDataCopyFP ObitDataCopy;
/** Function pointer to Clone. */
ObitDataCloneFP ObitDataClone;
/** Initialize IO member */
ObitDataSetupIOFP ObitDataSetupIO;
/** Function pointer to Rename. */
ObitDataRenameFP ObitDataRename;
/** Function pointer to Zap. */
ObitDataZapFP ObitDataZap;
/** Function pointer to Open. */
ObitDataOpenFP ObitDataOpen;
/** Function pointer to Close. */
ObitDataCloseFP ObitDataClose;
/** Function pointer to Reset IO. */
ObitDataIOSetFP ObitDataIOSet;
/** Function pointer to Create an associated table. */
newObitDataTableFP newObitDataTable;
/** Destroy an associated Table(s) */
ObitDataZapTableFP ObitDataZapTable;
/** Fully instantiate. */
ObitDataFullInstantiateFP ObitDataFullInstantiate;
/** Copy associated Table(s) */
ObitDataCopyTablesFP ObitDataCopyTables;
/** Update disk resident tables information */
ObitDataUpdateTablesFP ObitDataUpdateTables;
/** Copy a given table */
ObitDataCopyTableFP ObitDataCopyTable;
/** Function pointer to Create an associated History. */
newObitDataHistoryFP newObitDataHistory;
/** Function pointer to Write header keyword */
ObitDataWriteKeywordFP ObitDataWriteKeyword;
/** Function pointer to  Read header keyword*/
ObitDataReadKeywordFP ObitDataReadKeyword;
