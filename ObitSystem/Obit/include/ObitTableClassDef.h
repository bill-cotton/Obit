/* $Id$  */
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
/*  Define the basic components of the ObitTable  ClassInfo structure */
/* This is intended to be included in a classInfo structure definition*/
#include "ObitClassDef.h"  /* Parent class ClassInfo definition file */
/** Fully instantiate. */
ObitTableFullInstantiateFP ObitTableFullInstantiate;
/** Remove previous entries */
ObitTableClearRowsFP ObitTableClearRows;
/** Function pointer to Convert. */
ObitTableConvertFP ObitTableConvert;
/** Function pointer to Zap. */
ObitTableZapFP ObitTableZap;
/** Function pointer to Open. */
ObitTableOpenFP ObitTableOpen;
/** Function pointer to Close. */
ObitTableCloseFP ObitTableClose;
/** Function pointer to Read. */
ObitTableReadFP ObitTableRead;
/** Function pointer to Read with selection/calibration. */
ObitTableReadSelectFP ObitTableReadSelect;
/** Function pointer to Write. */
ObitTableWriteFP ObitTableWrite;
/** Function pointer to Read Row. */
ObitTableReadRowFP ObitTableReadRow;
/** Function pointer to Write Row. */
ObitTableWriteRowFP ObitTableWriteRow;
/** Function pointer to Attach Row. */
ObitTableSetRowFP ObitTableSetRow ;
/** Function pointer to Return table type. */
ObitTableGetTypeFP ObitTableGetType ;
/** Function pointer to Return table version. */
ObitTableGetVersionFP ObitTableGetVersion;
