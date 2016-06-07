/* $Id: ObitPointingDef.h 287 2011-02-04 08:35:31Z bill.cotton $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2016                                               */
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
/*  Define the basic components of the ObitPointing structure         */
/*  This is intended to be included in a class structure definition  */
/*  This class accesses data in the EVLA BDF format                  */
/**
 * \file ObitPointingDef.h
 * ObitPointing structure members for this and any derived classes.
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
/** how many elements in each row? */
olong nelem;
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
/** 0-rel Order number of pointingModelId */
olong ordpointingModelId;
/** 0-rel Order number of timeInterval */
olong ordtimeInterval;
/** 0-rel Order number of numSample */
olong ordnumSample;
/** 0-rel Order number of encoder */
olong ordencoder;
/** 0-rel Order number of pointingTracking */
olong ordpointingTracking;
/** 0-rel Order number of usePolynomials */
olong ordusePolynomials;
/** 0-rel Order number of timeOrigin */
olong ordtimeOrigin;
/** 0-rel Order number of numTerm */
olong ordnumTerm;
/** 0-rel Order number of pointingDirection */
olong ordpointingDirection;
/** 0-rel Order number of target */
olong ordtarget;
/** 0-rel Order number of offset */
olong ordoffset;
/** 0-rel Order number of overTheTop */
olong ordoverTheTop;
/** 0-rel Order number of sourceOffset */
olong ordsourceOffset;
/** 0-rel Order number of sourceOffsetReferenceCode */
olong ordsourceOffsetReferenceCode;
/** 0-rel Order number of sourceOffsetEquinox */
olong ordsourceOffsetEquinox;
/** 0-rel Order number of sampledTimeInterval */
olong ordsampledTimeInterval;
/** 0-rel Order number of atmosphericCorrection */
olong ordatmosphericCorrection;

