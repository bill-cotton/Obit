/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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
/*  Define the basic components of the ObitALMACalAtm structure       */
/*  This is intended to be included in a class structure definition  */
/*  This class accesses data in the EVLA BDF format                  */
/**
 * \file ObitALMACalAtmDef.h
 * ObitALMACalAtm structure members for this and any derived classes.
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
/** number of elements */
olong nelem;
/** 0-rel order of antennaName */
olong ordantennaName;
/** 0-rel order of receiverBand */
olong ordreceiverBand;
/** 0-rel order of basebandName */
olong ordbasebandName;
/** 0-rel order of calDataId */
olong ordcalDataId;
/** 0-rel order of calReductionId */
olong ordcalReductionId;
/** 0-rel order of startValidTime */
olong ordstartValidTime;
/** 0-rel order of endValidTime */
olong ordendValidTime;
/** 0-rel order of numFreq */
olong ordnumFreq;
/** 0-rel order of numLoad */
olong ordnumLoad;
/** 0-rel order of numReceptor */
olong ordnumReceptor;
/** 0-rel order of forwardEffSpectrum */
olong ordforwardEffSpectrum;
/** 0-rel order of frequencyRange */
olong ordfrequencyRange;
/** 0-rel order of groundPressure */
olong ordgroundPressure;
/** 0-rel order of groundRelHumidity */
olong ordgroundRelHumidity;
/** 0-rel order of frequencySpectrum */
olong ordfrequencySpectrum;
/** 0-rel order of groundTemperature */
olong ordgroundTemperature;
/** 0-rel order of polarizationTypes */
olong ordpolarizationTypes;
/** 0-rel order of powerSkySpectrum */
olong ordpowerSkySpectrum;
/** 0-rel order of powerLoadSpectrum */
olong ordpowerLoadSpectrum;
/** 0-rel order of syscalType */
olong ordsyscalType;
/** 0-rel order of tAtmSpectrum */
olong ordtAtmSpectrum;
/** 0-rel order of tRecSpectrum */
olong ordtRecSpectrum;
/** 0-rel order of tSysSpectrum */
olong ordtSysSpectrum;
/** 0-rel order of tauSpectrum */
olong ordtauSpectrum;
/** 0-rel order of tAtm */
olong ordtAtm;
/** 0-rel order of tRec */
olong ordtRec;
/** 0-rel order of tSys */
olong ordtSys;
/** 0-rel order of tau */
olong ordtau;
/** 0-rel order of water */
olong ordwater;
/** 0-rel order of waterError */
olong ordwaterError;
/** 0-rel order of alphaSpectrum */
olong ordalphaSpectrum;
/** 0-rel order of forwardEfficiency */
olong ordforwardEfficiency;
/** 0-rel order of forwardEfficiencyError */
olong ordforwardEfficiencyError;
/** 0-rel order of sbGain */
olong ordsbGain;
/** 0-rel order of sbGainError */
olong ordsbGainError;
/** 0-rel order of sbGainSpectrum */
olong ordsbGainSpectrum;
