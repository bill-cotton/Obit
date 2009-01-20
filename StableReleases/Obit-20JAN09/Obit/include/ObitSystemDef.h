/* $Id$ */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitSystem structure           */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitSystemDef.h
 * ObitSystem structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Threading info member object  */
ObitThread *thread;
/** Program name */
gchar *pgmName;
/** Program number */
olong pgmNumber;
/** Number of scratch files registered (not number current in list) */
olong numberScratch;
/** Number of AIPS disks */
olong numberAIPSdisk;
/** Number of FITS disks */
olong numberFITSdisk;
/** AIPS user number */
olong AIPSuser;
/** disk number of last scratch file assigned */
olong lastDisk;
 /** Number of entries in scratchList*/
olong number;
/** glib singly linked list for registered scratch files */
GSList* scratchList;
/** Error/message object */
ObitErr *err;
