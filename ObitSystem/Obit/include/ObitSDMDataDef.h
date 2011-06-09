/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010-2011                                          */
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
/*  Define the basic components of the ObitSDMData structure         */
/*  This is intended to be included in a class structure definition  */
/* This class accesses data in the EVLA SDM format                   */
/**
 * \file ObitSDMDataDef.h
 * ObitSDMData structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** Root of data directory */
gchar *DataRoot;
/** Schema version number */
olong schemaVersion;
/* Is this EVLA data? */
gboolean isEVLA;
/* Is this ALMA data? */
gboolean isALMA;
/** ASDM table */
ASDMTable* ASDMTab;
/** Main table */
ASDMMainTable* MainTab;
/** Antenna table */
ASDMAntennaTable* AntennaTab;
/** calData table */
ASDMcalDataTable* calDataTab;
/** calDevice table */
ASDMcalDeviceTable* calDeviceTab;
/** calPointing table */
ASDMcalPointingTable* calPointingTab;
/** CalReduction table */
ASDMCalReductionTable* CalReductionTab;
/** ConfigDescription table */
ASDMConfigDescriptionTable* ConfigDescriptionTab;
/** CorrelatorMode table */
ASDMCorrelatorModeTable* CorrelatorModeTab;
/** DataDescription table */
ASDMDataDescriptionTable* DataDescriptionTab;
/** Doppler table */
ASDMDopplerTable* DopplerTab;
/** ExecBlock table */
ASDMExecBlockTable* ExecBlockTab;
/** Feed table */
ASDMFeedTable* FeedTab;
/** Field table */
ASDMFieldTable* FieldTab;
/** Flag table */
ASDMFlagTable* FlagTab;
/** Pointing table */
ASDMPointingTable* PointingTab;
/** PointingModel table */
ASDMPointingModelTable* PointingModelTab;
/** Polarization table */
ASDMPolarizationTable* PolarizationTab;
/** Processor table */
ASDMProcessorTable* ProcessorTab;
/** Receiver table */
ASDMReceiverTable* ReceiverTab;
/** SBSummary table */
ASDMSBSummaryTable* SBSummaryTab;
/** Scan table */
ASDMScanTable* ScanTab;
/** Source table */
ASDMSourceTable* SourceTab;
/** SpectralWindow table */
ASDMSpectralWindowTable* SpectralWindowTab;
/** State table */
ASDMStateTable* StateTab;
/** Station table */
ASDMStationTable* StationTab;
/** Subscan table */
ASDMSubscanTable* SubscanTab;
/** SwitchCycle table */
ASDMSwitchCycleTable* SwitchCycleTab;
/** Weather table */
ASDMWeatherTable* WeatherTab;
/** SysCal table */
ASDMSysCalTable* SysCalTab;
/** SysPower table */
ASDMSysPowerTable* SysPowerTab;
/** XXXX table 
    ASDMXXXXTable* ASDMXXXXTab; */
/** Reference Frequency */
odouble refFreq;
/** Reference JD */
odouble refJD;
/** Integration time in days */
ofloat integTime;
