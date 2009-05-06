/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2009                                          */
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
/*  Define the basic components of the ObitFile InfoClass structure.  */
/* This is intended to be included in a classInfo structure definition*/
#include "ObitClassDef.h"  /* Parent class ClassInfo definition file. */
/** Function pointer to Does file exist? */
ObitFileExistFP ObitFileExist;
/** Function pointer to file size */
ObitFileSizeFP ObitFileSize;
/** Function pointer to file name */
ObitFileNameFP ObitFileName;
/** Function pointer to Rename */
ObitFileRenameFP ObitFileRename;
/** Function pointer to Zap */
ObitFileZapFP ObitFileZap;
/** Function pointer to ZapFile */
ObitFileZapFileFP ObitFileZapFile;
/** Function pointer to  Open */
ObitFileOpenFP ObitFileOpen;
/** Function pointer to  Close */
ObitFileCloseFP ObitFileClose;
/** Function pointer to position */
ObitFileSetPosFP ObitFileSetPos;
/** Function pointer to position at end */
ObitFileEndFP ObitFileEnd;
/** Function pointer to  Read */
ObitFileReadFP ObitFileRead;
/** Function pointer to  Read text line*/
ObitFileReadLineFP ObitFileReadLine;
/** Function pointer to  Write */
ObitFileWriteFP ObitFileWrite;
/** Function pointer to  Write text line*/
ObitFileWriteLineFP ObitFileWriteLine;
/** Function pointer to Pad block. */
ObitFilePadFP ObitFilePad;
/** Function pointer to Pad end of file. */
ObitFilePadFileFP ObitFilePadFile;
/** Function pointer to  Flush */
ObitFileFlushFP ObitFileFlush;
/** Function pointer to Adds error message to err if errno non zero */
ObitFileErrMsgFP ObitFileErrMsg;
