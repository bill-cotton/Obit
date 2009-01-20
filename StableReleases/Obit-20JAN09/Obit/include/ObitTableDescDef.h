/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003                                               */
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
/*  Define the basic components of the ObitTableDesc structure        */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitTableDescDef.h
 * ObitImage structure members.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Type of access to file (read, write...) */
ObitIOAccess access;
/** Number of rows */
olong nrow;
/** Length of a row in bytes */
olong lrow;
/** Length of a row in bytes for I/O purposes */
olong lrowIO;
/** file position in bytes of the beginning of the table */
olong startData;
/** current beginning row (1-rel) read/write in buffer  */
olong firstRow;
/** number of rows in buffer */
olong numRowBuff;
/** Number of columns (fields) */
olong nfield;
/** Table version number */
olong version;
/** Number of keywords stored in info */
olong nkey;
/** Sort order, column by primary and secondary key, 
 * positive => ascending order, negative => descending 
 * An additional + or - 256 indicates sorting by absolute value */
olong sort[2];
/** Physical (1-rel) order number of each field - as stored in memory */
olong *order;
/** Offset from the first word in a row to the first element in the field
 * in units of the data type of the field */
olong *offset;
/** Offset from the first word in a row to the first element in the field
 * in units of bytes */
olong *byteOffset;
/** Total number of element (of column type) */
olong *repeat;
/** Field dimension [MAXINFOELEMDIM] */
gint32 **dim;
/** Field data types as #ObitInfoType */
ObitInfoType *type;
/** Table name */
gchar *TableName;
/** Field labels */
gchar **FieldName;
/** Field units */
gchar **FieldUnit;
/** InfoList for other keywords */
ObitInfoList *info;
/** Offset for "_status" column */
olong statusOff;
