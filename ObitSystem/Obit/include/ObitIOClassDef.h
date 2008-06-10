/* $Id$                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2008                                          */
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
/*  Define the basic components of the Obit InfoClass structure       */
/* This is intended to be included in a classInfo structure definition*/
#include "ObitClassDef.h"  /* Parent class ClassInfo definition file */
/** Function pointer to  Constructor */
newObitIOFP newObitIO;
/** Function pointer to Are underlying structures the same?. */
ObitIOSameFP ObitIOSame;
/** Function pointer to Rename. */
ObitIORenameFP ObitIORename;
/** Function pointer to Zap. */
ObitIOZapFP ObitIOZap;
/** Function pointer to  Open */
ObitIOOpenFP ObitIOOpen;
/** Function pointer to  Close */
ObitIOCloseFP ObitIOClose;
/** Function pointer to  Init I/O */
ObitIOSetFP ObitIOSet;
/** Function pointer to  Read */
ObitIOReadFP ObitIORead;
/** Function pointer to  Read specifying start row */
ObitIOReadRowFP ObitIOReadRow;
/** Function pointer to  Read with selection */
ObitIOReadSelectFP ObitIOReadSelect;
/** Function pointer to  Read with selection specifying start row */
ObitIOReadRowSelectFP ObitIOReadRowSelect;
/** Function pointer to  Write */
ObitIOWriteFP ObitIOWrite;
/** Function pointer to  Write specifying start row */
ObitIOWriteRowFP ObitIOWriteRow;
/** Function pointer to  Flush */
ObitIOFlushFP ObitIOFlush;
/** Public:  Read Descriptor */
ObitIOReadDescriptorFP ObitIOReadDescriptor;
/** Public:  Write Descriptor */
ObitIOWriteDescriptorFP ObitIOWriteDescriptor;
/** Function pointer to  Create buffer */
ObitIOCreateBufferFP ObitIOCreateBuffer;
/** Function pointer to Destroy buffer */
ObitIOFreeBufferFP ObitIOFreeBuffer;
/** Function pointer to create a table. */
newObitIOTableFP newObitIOTable;
/** Function pointer to update table info. */
ObitIOUpdateTablesFP ObitIOUpdateTables;
