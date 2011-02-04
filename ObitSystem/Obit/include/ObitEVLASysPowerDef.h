/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2011                                               */
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
/*  Define the basic components of the ObitEVLASysPower structure         */
/*  This is intended to be included in a class structure definition  */
/*  This class accesses data in the EVLA BDF format                  */
/**
 * \file ObitEVLASysPowerDef.h
 * ObitEVLASysPower structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** data file name */
gchar *DataFile;
/** Obit File for I/O */
ObitFile *file;
/** Size of file */
ObitFilePos fileSize;
/** how many rows in table */
olong nrow;
/** current row in table */
olong curRow;
/** Number of bytes in the buffer */
olong nBytesInBuffer;
/** I/O and parsing buffer */
gchar *buffer;
/* Current buffer pointer */
gchar *current;
/** Does data need byte flipping? */
gboolean byteFlip;
/** 0-rel Order number of antennaId */
olong ordantennaId;
/** 0-rel Order number of spectralWindowId */
olong ordspectralWindowId;
/** 0-rel Order number of feedId */
olong ordfeedId;
/** 0-rel Order number of timeInterval */
olong ordtimeInterval;
/** 0-rel Order number of numReceptor */
olong ordnumReceptor;
/** 0-rel Order number of switchedPowerDifference */
olong ordswitchedPowerDifference;
/** 0-rel Order number of  switchedPowerSum*/
olong ordswitchedPowerSum;
/** 0-rel Order number of requantizerGain */
olong ordrequantizerGain;
