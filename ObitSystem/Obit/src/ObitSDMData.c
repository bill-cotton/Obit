/* $Id$        */
/*
Which ASDM tables have been coded (X)
X    ASDM.xml
X    Main.xml
X    Antenna.xml
X    CalAtmosphere.xml
     CalData.xml
X    CalDevice.xml
     CalPointing.xml
     CalReduction.xml
X    CalWVR.xml
X    ConfigDescription.xml
X    CorrelatorMode.xml
X    DataDescription.xml
     Doppler.xml
X    ExecBlock.xml
X    Feed.xml
X    Field.xml
X    Flag.xml
     PointingModel.xml
X    Pointing.xml
X    Polarization.xml
X    Processor.xml
     Receiver.xml
     SBSummary.xml
X    Scan.xml
X    Source.xml
X    SpectralWindow.xml
X    State.xml
X    Station.xml
X    Subscan.xml
X    SwitchCycle.xml
     SysCal.xml
X    SysPower.xml
X    Weather.xml
 */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010,2011                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
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

#include "ObitSDMData.h"
#include "ObitEVLASysPower.h"
#include "ObitFile.h"
#include "glib/gqsort.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSDMData.c
 * ObitSDMData class function definitions.
 * This class is derived from the Obit base class.
 * This class accesses data in the EVLA SDM format
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitSDMData";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitSDMDataClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSDMDataClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitSDMDataInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSDMDataClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSDMDataClassInfoDefFn (gpointer inClass);

/* ASDM Tables routines */
/** Private: Parse 32 bit integer from XML string  */
static olong ASDMparse_int(gchar *string, olong maxChar, 
			   gchar *prior, gchar **next);
/** Private: Parse 64 bit integer from XML string  */
static ollong ASDMparse_int64(gchar *string, olong maxChar, 
			      gchar *prior, gchar **next);
/** Private: Parse double from XML string  */
static odouble ASDMparse_dbl(gchar *string, olong maxChar, 
			     gchar *prior, gchar **next);
/** Private: Parse boolean from XML string  */
static gboolean ASDMparse_boo(gchar *string, olong maxChar, 
			      gchar *prior, gchar **next);
/** Private: Parse string from XML string  */
static gchar* ASDMparse_str(gchar *string, olong maxChar, 
			   gchar *prior, gchar **next);
/** Private: Parse quotedstring from XML string  */
static gchar* ASDMparse_quote_str(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next);
/** Private: Parse time from XML string  */
static odouble ASDMparse_time(gchar *string, olong maxChar, 
			     gchar *prior, gchar **next);

/** Private: Parse time range from XML string  */
static odouble* ASDMparse_timeRange(gchar *string, olong maxChar, 
				    gchar *prior, gchar **next);

/** Private: Parse time interval from XML string  */
static odouble ASDMparse_timeint(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next);

/** Private: Parse array of doubles from XML string  */
static odouble* ASDMparse_dblarray(gchar *string, olong maxChar, 
				   gchar *prior, gchar **next);
/** Private: Parse array of floats from XML string  */
static ofloat* ASDMparse_fltarray(gchar *string, olong maxChar, 
				   gchar *prior, gchar **next);
/** Private: Parse array of ints from XML string  */
static olong* ASDMparse_intarray(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next);
/** Private: Parse array of enums ending with the index from XML string  */
static olong* ASDMparse_enumarray(gchar *string, olong maxChar, 
				  gchar *prior, gchar **next);

/** Private: Parse string array from XML string  */
static gchar** ASDMparse_strarray(gchar *string, olong maxChar, 
				  gchar *prior, gchar **next);
/** Private: Parse boolean array from XML string  */
static gboolean* ASDMparse_booarray(gchar *string, olong maxChar, 
				    gchar *prior, gchar **next);

/** Private: Look up antennaMake enum */
static ObitASDMAntennaMake LookupAntennaMake(gchar *name);

/** Private: Look up baseband name enum */
static ObitASDMBasebandName LookupBasebandName(gchar *name);

/** Private: Look up axis name enum */
static ObitASDMAxisName LookupAxisName(gchar *name);

/** Private: Look up sideband rejection mode enum */
static ObitASDMSideBMode LookupSideBMode(gchar *name);

/** Private: Look up window function enum */
static ObitASDMWindowFn LookupWindowFn(gchar *name);

/** Private: Look up atmospheric correction enum */
static ObitASDMAtmPhCorr LookupAtmPhCorr(gchar *name);

/** Private: Parser for ASDM table from file */
static ASDMTable* ParseASDMTable(gchar *ASDMFile, 
				 ObitErr *err);
/** Private: Destructor for ASDMtable. */
static ASDMTable* KillASDMTable(ASDMTable* table);

/** Private: Destructor for Main table row. */
static ASDMMainRow* KillASDMMainRow(ASDMMainRow* row);
/** Private: Parser constructor for Main table from file */
static ASDMMainTable* ParseASDMMainTable(ObitSDMData *me,
					 gchar *MainFile, 
					 ObitErr *err);
/** Private: Destructor for Main table. */
static ASDMMainTable* KillASDMMainTable(ASDMMainTable* table);

/** Private: Destructor for Antenna table row. */
static ASDMAntennaRow* KillASDMAntennaRow(ASDMAntennaRow* row);
/** Private: Parser constructor for Antenna table from file */
static ASDMAntennaTable* ParseASDMAntennaTable(ObitSDMData *me,
					 gchar *AntennaFile, 
					 ObitErr *err);
/** Private: Destructor for Antenna table. */
static ASDMAntennaTable* KillASDMAntennaTable(ASDMAntennaTable* table);

/** Private: Destructor for calAtmosphere table row. */
static ASDMcalAtmosphereRow* KillASDMcalAtmosphereRow(ASDMcalAtmosphereRow* row);
/** Private: Parser constructor for calAtmosphere table from file */
static ASDMcalAtmosphereTable* ParseASDMcalAtmosphereTable(ObitSDMData *me,
					 gchar *calAtmosphereFile, 
					 ObitErr *err);
/** Private: Destructor for calAtmosphere table. */
static ASDMcalAtmosphereTable* KillASDMcalAtmosphereTable(ASDMcalAtmosphereTable* table);

/** Private: Destructor for calData table row. */
static ASDMcalDataRow* KillASDMcalDataRow(ASDMcalDataRow* row);
/** Private: Parser constructor for calData table from file */
static ASDMcalDataTable* ParseASDMcalDataTable(ObitSDMData *me,
					 gchar *calDataFile, 
					 ObitErr *err);
/** Private: Destructor for calData table. */
static ASDMcalDataTable* KillASDMcalDataTable(ASDMcalDataTable* table);

/** Private: Destructor for calDevice table row. */
static ASDMcalDeviceRow* KillASDMcalDeviceRow(ASDMcalDeviceRow* row);
/** Private: Parser constructor for calDevice table from file */
static ASDMcalDeviceTable* ParseASDMcalDeviceTable(ObitSDMData *me,
					 gchar *calDeviceFile, 
					 ObitErr *err);
/** Private: Destructor for calDevice table. */
static ASDMcalDeviceTable* KillASDMcalDeviceTable(ASDMcalDeviceTable* table);

/** Private: Destructor for calPointing table row. */
static ASDMcalPointingRow* KillASDMcalPointingRow(ASDMcalPointingRow* row);
/** Private: Parser constructor for calPointing table from file */
static ASDMcalPointingTable* ParseASDMcalPointingTable(ObitSDMData *me,
					 gchar *calPointingFile, 
					 ObitErr *err);
/** Private: Destructor for calPointing table. */
static ASDMcalPointingTable* KillASDMcalPointingTable(ASDMcalPointingTable* table);

/** Private: Destructor for CalReduction table row. */
static ASDMCalReductionRow* KillASDMCalReductionRow(ASDMCalReductionRow* row);
/** Private: Parser constructor for CalReduction table from file */
static ASDMCalReductionTable* ParseASDMCalReductionTable(ObitSDMData *me,
					 gchar *CalReductionFile, 
					 ObitErr *err);
/** Private: Destructor for CalReduction table. */
static ASDMCalReductionTable* KillASDMCalReductionTable(ASDMCalReductionTable* table);

/** Private: Destructor for CalWVR table row. */
static ASDMCalWVRRow* KillASDMCalWVRRow(ASDMCalWVRRow* row);
/** Private: Parser constructor for CalWVR table from file */
static ASDMCalWVRTable* ParseASDMCalWVRTable(ObitSDMData *me,
					 gchar *CalWVRFile, 
					 ObitErr *err);
/** Private: Destructor for CalWVR table. */
static ASDMCalWVRTable* KillASDMCalWVRTable(ASDMCalWVRTable* table);

/** Private: Destructor for ConfigDescription table row. */
static ASDMConfigDescriptionRow* KillASDMConfigDescriptionRow(ASDMConfigDescriptionRow* row);
/** Private: Parser constructor for ConfigDescription table from file */
static ASDMConfigDescriptionTable* ParseASDMConfigDescriptionTable(ObitSDMData *me,
					 gchar *ConfigDescriptionFile, 
					 ObitErr *err);
/** Private: Destructor for ConfigDescription table. */
static ASDMConfigDescriptionTable* 
KillASDMConfigDescriptionTable(ASDMConfigDescriptionTable* table);

/** Private: Destructor for CorrelatorMode table row. */
static ASDMCorrelatorModeRow* KillASDMCorrelatorModeRow(ASDMCorrelatorModeRow* row);
/** Private: Parser constructor for CorrelatorMode table from file */
static ASDMCorrelatorModeTable* ParseASDMCorrelatorModeTable(ObitSDMData *me,
					 gchar *CorrelatorModeFile, 
					 ObitErr *err);
/** Private: Destructor for CorrelatorMode table. */
static ASDMCorrelatorModeTable* KillASDMCorrelatorModeTable(ASDMCorrelatorModeTable* table);

/** Private: Destructor for DataDescription table row. */
static ASDMDataDescriptionRow* KillASDMDataDescriptionRow(ASDMDataDescriptionRow* row);
/** Private: Parser constructor for DataDescription table from file */
static ASDMDataDescriptionTable* ParseASDMDataDescriptionTable(ObitSDMData *me,
					 gchar *DataDescriptionFile, 
					 ObitErr *err);
/** Private: Destructor for DataDescription table. */
static ASDMDataDescriptionTable* KillASDMDataDescriptionTable(ASDMDataDescriptionTable* table);

/** Private: Destructor for Doppler table row. */
static ASDMDopplerRow* KillASDMDopplerRow(ASDMDopplerRow* row);
/** Private: Parser constructor for Doppler table from file */
static ASDMDopplerTable* ParseASDMDopplerTable(ObitSDMData *me,
					 gchar *DopplerFile, 
					 ObitErr *err);
/** Private: Destructor for Doppler table. */
static ASDMDopplerTable* KillASDMDopplerTable(ASDMDopplerTable* table);

/** Private: Destructor for ExecBlock table row. */
static ASDMExecBlockRow* KillASDMExecBlockRow(ASDMExecBlockRow* row);
/** Private: Parser constructor for ExecBlock table from file */
static ASDMExecBlockTable* ParseASDMExecBlockTable(ObitSDMData *me,
					 gchar *ExecBlockFile, 
					 ObitErr *err);
/** Private: Destructor for ExecBlock table. */
static ASDMExecBlockTable* KillASDMExecBlockTable(ASDMExecBlockTable* table);

/** Private: Destructor for Feed table row. */
static ASDMFeedRow* KillASDMFeedRow(ASDMFeedRow* row);
/** Private: Parser constructor for Feed table from file */
static ASDMFeedTable* ParseASDMFeedTable(ObitSDMData *me,
					 gchar *FeedFile, 
					 ObitErr *err);
/** Private: Destructor for Feed table. */
static ASDMFeedTable* KillASDMFeedTable(ASDMFeedTable* table);

/** Private: Destructor for Field table row. */
static ASDMFieldRow* KillASDMFieldRow(ASDMFieldRow* row);
/** Private: Parser constructor for Field table from file */
static ASDMFieldTable* ParseASDMFieldTable(ObitSDMData *me,
					 gchar *FieldFile, 
					 ObitErr *err);
/** Private: Destructor for Field table. */
static ASDMFieldTable* KillASDMFieldTable(ASDMFieldTable* table);

/** Private: Destructor for Flag table row. */
static ASDMFlagRow* KillASDMFlagRow(ASDMFlagRow* row);
/** Private: Parser constructor for Flag table from file */
static ASDMFlagTable* ParseASDMFlagTable(ObitSDMData *me,
					 gchar *FlagFile, 
					 ObitErr *err);
/** Private: Destructor for Flag table. */
static ASDMFlagTable* KillASDMFlagTable(ASDMFlagTable* table);

/** Private: Destructor for Pointing table row. */
static ASDMPointingRow* KillASDMPointingRow(ASDMPointingRow* row);
/** Private: Parser constructor for Pointing table from file */
static ASDMPointingTable* ParseASDMPointingTable(ObitSDMData *me,
					 gchar *PointingFile, 
					 ObitErr *err);
/** Private: Destructor for Pointing table. */
static ASDMPointingTable* KillASDMPointingTable(ASDMPointingTable* table);

/** Private: Destructor for PointingModel table row. */
static ASDMPointingModelRow* KillASDMPointingModelRow(ASDMPointingModelRow* row);
/** Private: Parser constructor for PointingModel table from file */
static ASDMPointingModelTable* ParseASDMPointingModelTable(ObitSDMData *me,
					 gchar *PointingModelFile, 
					 ObitErr *err);
/** Private: Destructor for PointingModel table. */
static ASDMPointingModelTable* KillASDMPointingModelTable(ASDMPointingModelTable* table);

/** Private: Destructor for Polarization table row. */
static ASDMPolarizationRow* KillASDMPolarizationRow(ASDMPolarizationRow* row);
/** Private: Parser constructor for Polarization table from file */
static ASDMPolarizationTable* ParseASDMPolarizationTable(ObitSDMData *me,
					 gchar *PolarizationFile, 
					 ObitErr *err);
/** Private: Destructor for Polarization table. */
static ASDMPolarizationTable* KillASDMPolarizationTable(ASDMPolarizationTable* table);

/** Private: Destructor for Processor table row. */
static ASDMProcessorRow* KillASDMProcessorRow(ASDMProcessorRow* row);
/** Private: Parser constructor for Processor table from file */
static ASDMProcessorTable* ParseASDMProcessorTable(ObitSDMData *me,
					 gchar *ProcessorFile, 
					 ObitErr *err);
/** Private: Destructor for Processor table. */
static ASDMProcessorTable* KillASDMProcessorTable(ASDMProcessorTable* table);

/** Private: Destructor for Receiver table row. */
static ASDMReceiverRow* KillASDMReceiverRow(ASDMReceiverRow* row);
/** Private: Parser constructor for Receiver table from file */
static ASDMReceiverTable* ParseASDMReceiverTable(ObitSDMData *me,
					 gchar *ReceiverFile, 
					 ObitErr *err);
/** Private: Destructor for Receiver table. */
static ASDMReceiverTable* KillASDMReceiverTable(ASDMReceiverTable* table);

/** Private: Destructor for SBSummary table row. */
static ASDMSBSummaryRow* KillASDMSBSummaryRow(ASDMSBSummaryRow* row);
/** Private: Parser constructor for SBSummary table from file */
static ASDMSBSummaryTable* ParseASDMSBSummaryTable(ObitSDMData *me,
					 gchar *SBSummaryFile, 
					 ObitErr *err);
/** Private: Destructor for SBSummary table. */
static ASDMSBSummaryTable* KillASDMSBSummaryTable(ASDMSBSummaryTable* table);

/** Private: Destructor for Scan table row. */
static ASDMScanRow* KillASDMScanRow(ASDMScanRow* row);
/** Private: Parser constructor for Scan table from file */
static ASDMScanTable* ParseASDMScanTable(ObitSDMData *me,
					 gchar *ScanFile, 
					 ObitErr *err);
/** Private: Destructor for Scan table. */
static ASDMScanTable* KillASDMScanTable(ASDMScanTable* table);

/** Private: Destructor for Source table row. */
static ASDMSourceRow* KillASDMSourceRow(ASDMSourceRow* row);
/** Private: Parser constructor for Source table from file */
static ASDMSourceTable* ParseASDMSourceTable(ObitSDMData *me,
					 gchar *SourceFile, 
					 ObitErr *err);
/** Private: Destructor for Source table. */
static ASDMSourceTable* KillASDMSourceTable(ASDMSourceTable* table);

/** Private: Destructor for SpectralWindow table row. */
static ASDMSpectralWindowRow* KillASDMSpectralWindowRow(ASDMSpectralWindowRow* row);
/** Private: Parser constructor for SpectralWindow table from file */
static ASDMSpectralWindowTable* ParseASDMSpectralWindowTable(ObitSDMData *me,
					 gchar *SpectralWindowFile, 
					 ObitErr *err);
/** Private: Destructor for SpectralWindow table. */
static ASDMSpectralWindowTable* KillASDMSpectralWindowTable(ASDMSpectralWindowTable* table);

/** Private: Destructor for State table row. */
static ASDMStateRow* KillASDMStateRow(ASDMStateRow* row);
/** Private: Parser constructor for State table from file */
static ASDMStateTable* ParseASDMStateTable(ObitSDMData *me,
					 gchar *StateFile, 
					 ObitErr *err);
/** Private: Destructor for State table. */
static ASDMStateTable* KillASDMStateTable(ASDMStateTable* table);

/** Private: Destructor for Station table row. */
static ASDMStationRow* KillASDMStationRow(ASDMStationRow* row);
/** Private: Parser constructor for Station table from file */
static ASDMStationTable* ParseASDMStationTable(ObitSDMData *me,
					 gchar *StationFile, 
					 ObitErr *err);
/** Private: Destructor for Station table. */
static ASDMStationTable* KillASDMStationTable(ASDMStationTable* table);

/** Private: Destructor for Subscan table row. */
static ASDMSubscanRow* KillASDMSubscanRow(ASDMSubscanRow* row);
/** Private: Parser constructor for Subscan table from file */
static ASDMSubscanTable* ParseASDMSubscanTable(ObitSDMData *me,
					 gchar *SubscanFile, 
					 ObitErr *err);
/** Private: Destructor for Subscan table. */
static ASDMSubscanTable* KillASDMSubscanTable(ASDMSubscanTable* table);

/** Private: Destructor for SwitchCycle table row. */
static ASDMSwitchCycleRow* KillASDMSwitchCycleRow(ASDMSwitchCycleRow* row);
/** Private: Parser constructor for SwitchCycle table from file */
static ASDMSwitchCycleTable* ParseASDMSwitchCycleTable(ObitSDMData *me,
					 gchar *SwitchCycleFile, 
					 ObitErr *err);
/** Private: Destructor for SwitchCycle table. */
static ASDMSwitchCycleTable* KillASDMSwitchCycleTable(ASDMSwitchCycleTable* table);

/** Private: Destructor for SysCal table row. */
static ASDMSysCalRow* KillASDMSysCalRow(ASDMSysCalRow* row);
/** Private: Parser constructor for SysCal table from file */
static ASDMSysCalTable* ParseASDMSysCalTable(ObitSDMData *me,
					 gchar *SysCalFile, 
					 ObitErr *err);
/** Private: Destructor for SysCal table. */
static ASDMSysCalTable* KillASDMSysCalTable(ASDMSysCalTable* table);

/** Private: Destructor for SysPower table row. */
static ASDMSysPowerRow* KillASDMSysPowerRow(ASDMSysPowerRow* row);
/** Private: Parser constructor for SysPower table from file */
static ASDMSysPowerTable* ParseASDMSysPowerTable(ObitSDMData *me,
					 gchar *SysPowerFile, 
					 ObitErr *err);
/** Private: Parser constructor for SysPowerl table from xml file */
static ASDMSysPowerTable* ParseASDMSysPowerTableXML(ObitSDMData *me,
						    gchar *SysPowerFile, 
						    ObitErr *err);
/** Private: Parser constructor for SysPower table from binary file */
static ASDMSysPowerTable* ParseASDMSysPowerTableBin(ObitSDMData *me,
						    gchar *SysPowerFile, 
						    ObitErr *err);
/** Private: Destructor for SysPower table. */
static ASDMSysPowerTable* KillASDMSysPowerTable(ASDMSysPowerTable* table);

/** Private: Destructor for Weather table row. */
static ASDMWeatherRow* KillASDMWeatherRow(ASDMWeatherRow* row);
/** Private: Parser constructor for Weather table from file */
static ASDMWeatherTable* ParseASDMWeatherTable(ObitSDMData *me,
					       gchar *WeatherFile, 
					       ObitErr *err);
/** Private: Destructor for Weather table. */
static ASDMWeatherTable* KillASDMWeatherTable(ASDMWeatherTable* table);
/** Private: Count rows in a table */
static olong CountTableRows(gchar *WeatherFile, ObitErr *err);

/** Private: Destructor for XXXX table row. */
static ASDMXXXXRow* KillASDMXXXXRow(ASDMXXXXRow* row);
/** Private: Parser constructor for XXXX table from file */
static ASDMXXXXTable* ParseASDMXXXXTable(ObitSDMData *me,
					 gchar *XXXXFile, 
					 ObitErr *err);
/** Private: Destructor for XXXX table. */
static ASDMXXXXTable* KillASDMXXXXTable(ASDMXXXXTable* table);
/** Private: Count rows in a table */
static olong CountTableRows(gchar *XXXXFile, ObitErr *err);

/** Private: Sort comparison function for frequencies */
static gint CompareFreq (gconstpointer in1, gconstpointer in2, 
			 gpointer ncomp);

/** Private: Squeeze all blanks out of a string */
static void Strip (gchar *s);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSDMData* newObitSDMData (gchar* name)
{
  ObitSDMData* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSDMDataClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSDMData));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSDMDataInit((gpointer)out);

 return out;
} /* end newObitSDMData */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSDMDataGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSDMDataClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitSDMDataGetClass */

/**
 * Make a deep copy of an ObitSDMData. NYI
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitSDMData* ObitSDMDataCopy  (ObitSDMData *in, ObitSDMData *out, ObitErr *err)
{
  /*const ObitClassInfo *ParentClass;*/
  /*gboolean oldExist;*/
  /*gchar *outName;*/

  /* error checks */
  if (err->error) return out;

  /* Stubbed */
  g_error("ObitSDMDataCopy: Stubbed");

  return out;
} /* end ObitSDMDataCopy */

/**
 * Make a copy of a object but do not copy the actual data NYI
 * This is useful to create an BDFData similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitSDMDataClone  (ObitSDMData *in, ObitSDMData *out, ObitErr *err)
{
  /*const ObitClassInfo *ParentClass;*/

  /* error checks */
  if (err->error) return;

  /* Stubbed */
  g_error("ObitSDMDataCopy: Stubbed");

} /* end ObitSDMDataClone */

/**
 * Creates an ObitSDMData 
 * Parses the ASMD XML tables and stores
 * \param name     An optional name for the object.
 * \param DataRoot Directory root of data
 * \return the new object.
 */
ObitSDMData* ObitSDMDataCreate (gchar* name, gchar *DataRoot, ObitErr *err)
{
  gchar *fullname;
  ObitSDMData* out;
  ASDMSpectralWindowArray *damn=NULL;
  gchar *routine="ObitSDMDataCreate";

  /* Create basic structure */
  out = newObitSDMData (name);

  /* set Values */
  out->DataRoot  = strdup(DataRoot);
  out->isEVLA    = FALSE;
  out->isALMA    = FALSE;
  out->selConfig = -1;
  out->selBand   = ASDMBand_Any;

  /* Parsing buffer */
  out->maxLine = 32768;
  out->line = g_malloc(out->maxLine+1);

  /* ASDM table - also get schemaVersion */
  fullname = g_strconcat (DataRoot,"/ASDM.xml", NULL);
  out->ASDMTab = ParseASDMTable(fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);
  out->schemaVersion = out->ASDMTab->schemaVersion;

  /* ExecBlock table - also determines which array */
  fullname = g_strconcat (DataRoot,"/ExecBlock.xml", NULL);
  out->ExecBlockTab = ParseASDMExecBlockTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Main table */
  fullname = g_strconcat (DataRoot,"/Main.xml", NULL);
  out->MainTab = ParseASDMMainTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Antenna table */
  fullname = g_strconcat (DataRoot,"/Antenna.xml", NULL);
  out->AntennaTab = ParseASDMAntennaTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* calAtmosphere table */
  fullname = g_strconcat (DataRoot,"/CalAtmosphere.xml", NULL);
  out->calAtmosphereTab = ParseASDMcalAtmosphereTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* calData table */
  fullname = g_strconcat (DataRoot,"/CalData.xml", NULL);
  out->calDataTab = ParseASDMcalDataTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* calDevice table */
  fullname = g_strconcat (DataRoot,"/CalDevice.xml", NULL);
  out->calDeviceTab = ParseASDMcalDeviceTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* calPointing table */
  fullname = g_strconcat (DataRoot,"/CalPointing.xml", NULL);
  out->calPointingTab = ParseASDMcalPointingTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* CalReduction table */
  fullname = g_strconcat (DataRoot,"/CalReduction.xml", NULL);
  out->CalReductionTab = ParseASDMCalReductionTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* CalWVR table */
  fullname = g_strconcat (DataRoot,"/CalWVR.xml", NULL);
  out->CalWVRTab = ParseASDMCalWVRTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* ConfigDescription table */
  fullname = g_strconcat (DataRoot,"/ConfigDescription.xml", NULL);
  out->ConfigDescriptionTab = ParseASDMConfigDescriptionTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* CorrelatorMode table */
  fullname = g_strconcat (DataRoot,"/CorrelatorMode.xml", NULL);
  out->CorrelatorModeTab = ParseASDMCorrelatorModeTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* DataDescription table */
  fullname = g_strconcat (DataRoot,"/DataDescription.xml", NULL);
  out->DataDescriptionTab = ParseASDMDataDescriptionTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Doppler table */
  fullname = g_strconcat (DataRoot,"/Doppler.xml", NULL);
  out->DopplerTab = ParseASDMDopplerTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Feed table */
  fullname = g_strconcat (DataRoot,"/Feed.xml", NULL);
  out->FeedTab = ParseASDMFeedTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Field table */
  fullname = g_strconcat (DataRoot,"/Field.xml", NULL);
  out->FieldTab = ParseASDMFieldTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Flag table */
  fullname = g_strconcat (DataRoot,"/Flag.xml", NULL);
  out->FlagTab = ParseASDMFlagTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Pointing table */
  fullname = g_strconcat (DataRoot,"/Pointing.xml", NULL);
  out->PointingTab = ParseASDMPointingTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* PointingModel table */
  fullname = g_strconcat (DataRoot,"/PointingModel.xml", NULL);
  out->PointingModelTab = ParseASDMPointingModelTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Polarization table */
  fullname = g_strconcat (DataRoot,"/Polarization.xml", NULL);
  out->PolarizationTab = ParseASDMPolarizationTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Processor table */
  fullname = g_strconcat (DataRoot,"/Processor.xml", NULL);
  out->ProcessorTab = ParseASDMProcessorTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Receiver table */
  fullname = g_strconcat (DataRoot,"/Receiver.xml", NULL);
  out->ReceiverTab = ParseASDMReceiverTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* SBSummary table */
  fullname = g_strconcat (DataRoot,"/SBSummary.xml", NULL);
  out->SBSummaryTab = ParseASDMSBSummaryTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Scan table */
  fullname = g_strconcat (DataRoot,"/Scan.xml", NULL);
  out->ScanTab = ParseASDMScanTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Source table */
  fullname = g_strconcat (DataRoot,"/Source.xml", NULL);
  out->SourceTab = ParseASDMSourceTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* SpectralWindow table */
  fullname = g_strconcat (DataRoot,"/SpectralWindow.xml", NULL);
  out->SpectralWindowTab = ParseASDMSpectralWindowTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* State table */
  fullname = g_strconcat (DataRoot,"/State.xml", NULL);
  out->StateTab = ParseASDMStateTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Station table */
  fullname = g_strconcat (DataRoot,"/Station.xml", NULL);
  out->StationTab = ParseASDMStationTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Subscan table */
  fullname = g_strconcat (DataRoot,"/Subscan.xml", NULL);
  out->SubscanTab = ParseASDMSubscanTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* SwitchCycle table */
  fullname = g_strconcat (DataRoot,"/SwitchCycle.xml", NULL);
  out->SwitchCycleTab = ParseASDMSwitchCycleTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* SysCal table */
  fullname = g_strconcat (DataRoot,"/SysCal.xml", NULL);
  out->SysCalTab = ParseASDMSysCalTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* SysPower table from either xml or binary */
  fullname = g_strconcat (DataRoot,"/SysPower", NULL);
  out->SysPowerTab = ParseASDMSysPowerTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Weather table */
  fullname = g_strconcat (DataRoot,"/Weather.xml", NULL);
  out->WeatherTab = ParseASDMWeatherTable(out, fullname, err);
  if (err->error) Obit_traceback_val (err, routine, fullname, out);
  g_free(fullname);

  /* Other info - what a piece of shit */
  damn = ObitSDMDataGetSWArray (out, 0, FALSE);
  Obit_retval_if_fail((damn!=NULL), err, out,
		      "%s: Failed to generate Spectral Window Array", 
		      routine);

  /* Reference JD from SW Array */
  out->refJD = damn->refJD;

  /* Reference Frequency */
  out->refFreq = damn->refFreq;

  damn = ObitSDMDataKillSWArray(damn);

  /* Get integration time from the first scan */
  out->integTime = out->MainTab->rows[0]->interval / 
    MAX(1, out->MainTab->rows[0]->numIntegration);

  /* release parsing buffer */
  out->maxLine = 0;
  g_free(out->line); out->line = NULL;

  return out;
} /* end ObitSDMDataCreate */

/**
 * Creates and fills n spectral window array
 * Parses the ASMD XML tables and stores
 * \param in       ASDM object to use
 * \param mainRow  Scan row number (in ASDMMain table)
 * \param SWOrder  If TRUE leave data is SW order
 * \return the new structure, NULL on error, delete using ObitSDMKillSWArray
 */
ASDMSpectralWindowArray* ObitSDMDataGetSWArray (ObitSDMData *in, olong mainRow, 
						gboolean SWOrder)
{ 
  ASDMSpectralWindowArray* out=NULL;
  olong configDescriptionId, dataDescriptionId, spectralWindowId, *dataDescriptions;
  olong jPoln, polarizationId;
  olong i, j, k, iMain, iConfig, jSW, iSW, jDD, numDD, iJD, ncomp;
  odouble JD;
  gboolean first = TRUE, isALMA;
  ASDMSpectralWindowArrayEntry **twinds=NULL;
  ofloat *sortStruct=NULL;
  gint number;
  size_t size;

  /* Find scan in Main table */
  iMain = mainRow;
  
  if (iMain>=in->MainTab->nrows) return out;

  /* Is this ALMA data? */
  isALMA = !strncmp(in->ExecBlockTab->rows[0]->telescopeName, "ALMA", 4);
  /* ALMA may be a bit confused */
  if (!isALMA) isALMA = !strncmp(in->ExecBlockTab->rows[0]->telescopeName, "OSF", 3);

  out = g_malloc0(sizeof(ASDMSpectralWindowArray));

  /* Assume no more spectral windows than total defined */
  out->winds = g_malloc0(in->SpectralWindowTab->nrows*sizeof(ASDMSpectralWindowArrayEntry));

  /* Find entry  in configDescription table */
  configDescriptionId = in->MainTab->rows[iMain]->configDescriptionId;
  for (iConfig=0; iConfig<in->ConfigDescriptionTab->nrows; iConfig++) {
    if (in->ConfigDescriptionTab->rows[iConfig]->configDescriptionId==configDescriptionId) break;
  }
  if (iConfig>=in->ConfigDescriptionTab->nrows) return NULL;

  /* Reference frequency - get from Spectral window on first dataDescription */
  /* Find Data description */
  dataDescriptionId = in->ConfigDescriptionTab->rows[iConfig]->dataDescriptionId[0];
  for (jDD=0; jDD<in->DataDescriptionTab->nrows; jDD++) {
    if (in->DataDescriptionTab->rows[jDD]->dataDescriptionId==dataDescriptionId) break;
  }
  if (jDD>=in->DataDescriptionTab->nrows) return NULL;

  /* Find spectralWindow */
  spectralWindowId = in->DataDescriptionTab->rows[jDD]->spectralWindowId;
  for (jSW=0; jSW<in->SpectralWindowTab->nrows; jSW++) {
    if (in->SpectralWindowTab->rows[jSW]->spectralWindowId==spectralWindowId) break;
  }
  if (jSW>=in->SpectralWindowTab->nrows) return NULL;
  /* Finally */
  out->refFreq = in->SpectralWindowTab->rows[jSW]->chanFreqStart;

  /* Reference JD (0 h ) */
  JD = in->MainTab->rows[iMain]->time;
  iJD = (olong)(JD-0.5);
  /* To 0 hours */
  out->refJD = iJD + 0.5;

  /* Data descriptions array */
  dataDescriptions = in->ConfigDescriptionTab->rows[iConfig]->dataDescriptionId;
  
  /* Count number of data descriptions */
  numDD = 0;
  i = 0;
  while(dataDescriptions[i++]>=0) {numDD++;}

  /* Create/copy structure/info */
  iSW = 0;
  for (i=0; i<numDD; i++) {
    /* Find Data descriptions */

    dataDescriptionId = dataDescriptions[iSW];
    for (jDD=0; jDD<in->DataDescriptionTab->nrows; jDD++) {
      if (in->DataDescriptionTab->rows[jDD]->dataDescriptionId==dataDescriptionId) break;
    }
    if (jDD>=in->DataDescriptionTab->nrows) return NULL;

    /* Find spectralWindow */
    spectralWindowId = in->DataDescriptionTab->rows[jDD]->spectralWindowId;
    for (jSW=0; jSW<in->SpectralWindowTab->nrows; jSW++) {
      if (in->SpectralWindowTab->rows[jSW]->spectralWindowId==spectralWindowId) break;
    }
    if (jSW>=in->SpectralWindowTab->nrows) return NULL;

    out->winds[iSW] = g_malloc(sizeof(ASDMSpectralWindowArrayEntry));
    out->winds[iSW]->spectralWindowId = in->SpectralWindowTab->rows[jSW]->spectralWindowId;
    out->winds[iSW]->selected         = TRUE;
    out->winds[iSW]->numChan          = in->SpectralWindowTab->rows[jSW]->numChan;
    out->winds[iSW]->netSideband      = g_strdup(in->SpectralWindowTab->rows[jSW]->netSideband);
    out->winds[iSW]->refFreq          = in->SpectralWindowTab->rows[jSW]->refFreq;
    out->winds[iSW]->totBandwidth     = in->SpectralWindowTab->rows[jSW]->totBandwidth;
    /* Which world view - single channel info or array of everything */
    if (in->SpectralWindowTab->rows[jSW]->chanFreqArray==NULL) { /* Single */
      out->winds[iSW]->chanFreqStart    = in->SpectralWindowTab->rows[jSW]->chanFreqStart;
      out->winds[iSW]->chanFreqStep     = in->SpectralWindowTab->rows[jSW]->chanFreqStep;
      out->winds[iSW]->chanWidth        = in->SpectralWindowTab->rows[jSW]->chanWidth;
      out->winds[iSW]->effectiveBw      = in->SpectralWindowTab->rows[jSW]->effectiveBw;
      out->winds[iSW]->resolution       = in->SpectralWindowTab->rows[jSW]->resolution;
    } else { /* Potential trouble - cannot deal with general case */
      out->winds[iSW]->chanFreqStart    = in->SpectralWindowTab->rows[jSW]->chanFreqArray[0];
      out->winds[iSW]->chanFreqStep     = in->SpectralWindowTab->rows[jSW]->chanFreqArray[1] - 
	in->SpectralWindowTab->rows[jSW]->chanFreqArray[0];
      out->winds[iSW]->chanWidth        = in->SpectralWindowTab->rows[jSW]->chanWidthArray[0];
      out->winds[iSW]->effectiveBw      = in->SpectralWindowTab->rows[jSW]->effectiveBwArray[0];
      out->winds[iSW]->resolution       = in->SpectralWindowTab->rows[jSW]->resolutionArray[0];
      /* Can I cope with/want this one??? */
      for (k=1; k<in->SpectralWindowTab->rows[jSW]->numChan; k++) {
	if (in->SpectralWindowTab->rows[jSW]->chanWidthArray[0] != 
	    in->SpectralWindowTab->rows[jSW]->chanWidthArray[k])  out->winds[iSW]->selected = FALSE;
	if (in->SpectralWindowTab->rows[jSW]->effectiveBwArray[0] != 
	    in->SpectralWindowTab->rows[jSW]->effectiveBwArray[k])  out->winds[iSW]->selected = FALSE;
	if (in->SpectralWindowTab->rows[jSW]->resolutionArray[0] != 
	    in->SpectralWindowTab->rows[jSW]->resolutionArray[k])  out->winds[iSW]->selected = FALSE;
      }
    }

    /* Fix up frequency for LSB - DEBUG STUB */
    if (out->winds[iSW]->netSideband[0]=='$') 
      out->winds[iSW]->chanFreqStart -= (out->winds[iSW]->numChan-1) * out->winds[iSW]->chanFreqStep;

    /* Set reference frequency */
    if (first) {
      out->refFreq = out->winds[iSW]->chanFreqStart;
      first = FALSE;
    }

    /* Find Polarization */
    polarizationId = in->DataDescriptionTab->rows[jDD]->polOrHoloId;
    for (jPoln=0; jPoln<in->PolarizationTab->nrows; jPoln++) {
      if (in->PolarizationTab->rows[jPoln]->polarizationId==polarizationId) break;
    }
    if (jPoln>=in->PolarizationTab->nrows) return NULL;

    out->winds[iSW]->nCPoln = in->PolarizationTab->rows[jPoln]->numCorr;
    out->winds[iSW]->nAPoln = MIN (3,in->PolarizationTab->rows[jPoln]->numCorr);

    /* Baseband - assume name starts with "BB_" and ends in number */
    out->winds[iSW]->basebandNum = 
      (olong)strtol(&in->SpectralWindowTab->rows[jSW]->basebandName[3], NULL, 10);

    /* Subband number - assume name starts with "Subband:" and ends in number */
    if (in->SpectralWindowTab->rows[jSW]->name) {
      out->winds[iSW]->subbandNum  = 
	(olong)strtol(&in->SpectralWindowTab->rows[jSW]->name[8], NULL, 10);
    } else out->winds[iSW]->subbandNum = out->winds[iSW]->basebandNum;

    /* Is this one selected? */
    if ((in->selChan>0) && (in->selChan!=out->winds[iSW]->numChan)) 
      out->winds[iSW]->selected = FALSE;
    if ((in->selBW>0.0) && (fabs(in->selBW-out->winds[iSW]->chanWidth)>100.0))
      out->winds[iSW]->selected = FALSE;

    iSW++;
  } /* end loop over DataDescriptions */

  out->nwinds = iSW;

  /* Sort into ascending reference frequencies - o
     float precision should be good enough */
  out->order = g_malloc0(out->nwinds*sizeof(olong));
  for (i=0; i<out->nwinds; i++) out->order[i] = i;

  /* Need to sort? */
  if (!SWOrder) {
    sortStruct = g_malloc0((2*out->nwinds+5)*sizeof(ofloat));
    for (i=0; i<out->nwinds; i++) {
      sortStruct[2*i]   = (ofloat)i;
      sortStruct[2*i+1] = (ofloat)out->winds[i]->chanFreqStart;
    }
    /* Sort */
    number = out->nwinds;
    size   = 2*sizeof(ofloat);
    ncomp  = 1;
    g_qsort_with_data (sortStruct, number, size, CompareFreq, &ncomp);
    
    /* Save sorted results */
    /* save initial windows to temporary array */
    twinds     = out->winds;
    out->winds = g_malloc0(in->SpectralWindowTab->nrows*sizeof(ASDMSpectralWindowArrayEntry));
    for (i=0; i<out->nwinds; i++) {
      j = (olong)(sortStruct[i*2]+0.5);
      out->winds[i] = twinds[j];
      out->order[i] = j;
    }
  } /* end sorting frequencies */

  /* Set reference frequency to first ordered Spectral window */
  out->refFreq = out->winds[0]->chanFreqStart;

  /* If there is only one ALMA window, it must be selected 
     this is needed for ALMA WVR data to find sources */
  if (isALMA && (out->nwinds==1)) out->winds[0]->selected = TRUE;

  /* Cleanup */
  if (sortStruct) g_free(sortStruct);
  if (twinds)     g_free(twinds);
  return out;
} /* end ObitSDMDataGetSWArray */

/**
 * Delete an ASDMSpectralWindowArray
 * \param in  Object to delete
 * \return NULL pointer
 */
ASDMSpectralWindowArray* ObitSDMDataKillSWArray (ASDMSpectralWindowArray *in)
{ 
  olong i;

  if (in==NULL) return NULL;
  /* delete entries */
  if (in->order) g_free(in->order);
  if (in->winds) {
    for (i=0; i<in->nwinds; i++) {
      if (in->winds[i]) {
	if (in->winds[i]->netSideband) g_free(in->winds[i]->netSideband);
	g_free(in->winds[i]);
      }
    }
    g_free(in->winds);
  }

  g_free(in);
  return NULL;
} /* end ObitSDMDataKillSWArray */

 /**
 * Select Spectral windows by number of channels/band
 * \param in       The structure to update
 * \param selChan selected number of channels
 * \param band    Selected band
 */
gboolean ObitSDMDataSelChan  (ASDMSpectralWindowArray *in, olong selChan,
			      olong selIF, ObitASDMBand band)
{
  gboolean out = FALSE, damn;
  olong iSW;
  /*gchar *routine = "ObitSDMDataSelChan";*/
  
  for (iSW=0; iSW<in->nwinds; iSW++) {
    damn = (in->winds[iSW]->numChan == selChan) &&
      ((ObitSDMDataFreq2Band (in->winds[iSW]->refFreq)==band) || 
       (band==ASDMBand_Any)) &&   (in->nwinds==selIF);
    in->winds[iSW]->selected = in->winds[iSW]->selected || damn;
    if (in->winds[iSW]->selected) out = TRUE;
  }
  return out;
} /* end ObitSDMDataSelChan */

/**
 * Creates and fills an antenna/station array
 * Parses the ASMD XML tables and stores
 * \param in       ASDM object to use
 * \param mainRow Main table number (in ASDMMain table)
 * \return the new structure, NULL on error, delete using ObitSDMKillAntArray
 */
ASDMAntennaArray* ObitSDMDataGetAntArray (ObitSDMData *in, olong mainRow)
{ 
  ASDMAntennaArray* out=NULL;
  olong  configDescriptionId, stationId, dataDescriptionId, spectralWindowId, execBlockId;
  olong i, iMain, iConfig, iAnt, jAnt, jDD, jSW, numAnt, iJD, iExec;
  odouble JD;

  /* Find scan in Main table */
  iMain = mainRow;
  if (iMain>=in->MainTab->nrows) return out;

  out = g_malloc0(sizeof(ASDMAntennaArray));

  /* Assume no more antennas than total defined */
  out->ants = g_malloc0(in->AntennaTab->nrows*sizeof(ASDMAntennaArrayEntry));

  /* How many? */
  numAnt = in->AntennaTab->nrows;
  out->nants  = numAnt;
  out->maxAnt = numAnt;

  /* Find entry  in configDescription table */
  configDescriptionId = in->MainTab->rows[iMain]->configDescriptionId;
  for (iConfig=0; iConfig<in->ConfigDescriptionTab->nrows; iConfig++) {
    if (in->ConfigDescriptionTab->rows[iConfig]->configDescriptionId==configDescriptionId) break;
  }
  if (iConfig>=in->ConfigDescriptionTab->nrows) return NULL;

  /* Find Exec block */
  execBlockId = in->MainTab->rows[iMain]->execBlockId;
  for (iExec=0; iExec<in->ExecBlockTab->nrows; iExec++) {
    if (in->ExecBlockTab->rows[iExec]->execBlockId==execBlockId) break;
  }
  if (iExec>=in->ExecBlockTab->nrows) return NULL;

  /* Array name */
  out->arrayName = g_strdup(in->ExecBlockTab->rows[iExec]->telescopeName); 

  /* Observer */
  if (strlen(in->ExecBlockTab->rows[iExec]->observerName)>1)
    out->obsName = g_strdup(in->ExecBlockTab->rows[iExec]->observerName); 
  else 
    out->obsName = g_strdup("Anon."); 

  /* Reference frequency - get from Spectral window on first dataDescription */
  /* Find Data description */
  dataDescriptionId = in->ConfigDescriptionTab->rows[iConfig]->dataDescriptionId[0];
  for (jDD=0; jDD<in->DataDescriptionTab->nrows; jDD++) {
    if (in->DataDescriptionTab->rows[jDD]->dataDescriptionId==dataDescriptionId) break;
  }
  if (jDD>=in->DataDescriptionTab->nrows) return NULL;

  /* Find spectralWindow */
  spectralWindowId = in->DataDescriptionTab->rows[jDD]->spectralWindowId;
  for (jSW=0; jSW<in->SpectralWindowTab->nrows; jSW++) {
    if (in->SpectralWindowTab->rows[jSW]->spectralWindowId==spectralWindowId) break;
  }
  if (jSW>=in->SpectralWindowTab->nrows) return NULL;
  /* Finally */
  out->refFreq = in->SpectralWindowTab->rows[jSW]->refFreq;

  /* Reference JD (0 h ) */
  JD = in->MainTab->rows[iMain]->time;
  iJD = (olong)(JD-0.5);
  /* To 0 hours */
  out->refJD = iJD + 0.5;

  /* Loop over antennas */
  for (iAnt=0; iAnt<numAnt; iAnt++) {
    out->ants[iAnt] = g_malloc0(sizeof(ASDMAntennaArrayEntry));
    
    /* Save antenna info */
    jAnt                          = iAnt;
    out->ants[iAnt]->antennaId    = in->AntennaTab->rows[jAnt]->antennaId;
    out->ants[iAnt]->antName      = g_strdup(in->AntennaTab->rows[jAnt]->name);
    out->ants[iAnt]->antennaMake  = in->AntennaTab->rows[jAnt]->antennaMake;
    out->ants[iAnt]->antennaType  = in->AntennaTab->rows[jAnt]->antennaType;
    out->ants[iAnt]->dishDiameter = in->AntennaTab->rows[jAnt]->dishDiameter;
    out->ants[iAnt]->antPosition  = g_malloc0(3*sizeof(odouble));
    for (i=0; i<3; i++) out->ants[iAnt]->antPosition[i] = in->AntennaTab->rows[jAnt]->position[i];
    out->ants[iAnt]->offset = g_malloc0(3*sizeof(odouble));
    for (i=0; i<3; i++) out->ants[iAnt]->offset[i] = in->AntennaTab->rows[jAnt]->offset[i];
    out->ants[iAnt]->time = in->AntennaTab->rows[jAnt]->time;
    out->ants[iAnt]->stationId    = in->AntennaTab->rows[jAnt]->stationId;
    /* Crack name to get number Assume EVLA starts with "ea" */
    if ((out->ants[iAnt]->antName[0]=='e') && (out->ants[iAnt]->antName[1]=='a'))
      out->ants[iAnt]->antennaNo = strtol(&out->ants[iAnt]->antName[2],NULL,10);
    else out->ants[iAnt]->antennaNo = out->ants[iAnt]->antennaId;
    out->maxAnt = MAX(out->maxAnt, out->ants[iAnt]->antennaNo);

    /* Find Station */
    stationId = out->ants[iAnt]->stationId;
    for (jAnt=0; jAnt<in->StationTab->nrows; jAnt++) {
      if (in->StationTab->rows[jAnt]->stationId==stationId) break;
    }
    if (jAnt>=in->StationTab->nrows) return NULL;

    /* Save station info */
    out->ants[iAnt]->staName     = g_strdup(in->StationTab->rows[jAnt]->name);
    out->ants[iAnt]->staPosition = g_malloc0(3*sizeof(odouble));
    for (i=0; i<3; i++) out->ants[iAnt]->staPosition[i] = in->StationTab->rows[jAnt]->position[i];
    out->ants[iAnt]->type        = in->StationTab->rows[jAnt]->type;

    /* Assume polarization for first PolarizationTab entry */
    out->ants[iAnt]->numPoln     = MIN(2, in->PolarizationTab->rows[0]->numCorr);
    out->ants[iAnt]->numPolnCorr = in->PolarizationTab->rows[0]->numCorr;
    out->ants[iAnt]->polnType[0] = in->PolarizationTab->rows[0]->corrType[0][0];
    if (out->ants[iAnt]->numPolnCorr>1)
      out->ants[iAnt]->polnType[1] = in->PolarizationTab->rows[0]->corrType[1][1];
  } /* end loop over antennas */
  
  return out;
} /* end ObitSDMDataGetAntArray */

/**
 * Delete an ASDMAntennaArray
 * \param in  Object to delete
 * \return NULL pointer
 */
ASDMAntennaArray* ObitSDMDataKillAntArray (ASDMAntennaArray *in)
{ 
  olong i;

  if (in==NULL) return NULL;

  /* delete entries */
  if (in->ants) {
    for (i=0; i<in->nants; i++) {
      if (in->ants[i]) {
	if (in->ants[i]->antName)     g_free(in->ants[i]->antName);
	if (in->ants[i]->antPosition) g_free(in->ants[i]->antPosition);
	if (in->ants[i]->offset)      g_free(in->ants[i]->offset);
	if (in->ants[i]->staName)     g_free(in->ants[i]->staName);
	if (in->ants[i]->staPosition) g_free(in->ants[i]->staPosition);
	g_free(in->ants[i]);
      }
    }
    g_free(in->ants);
  }
  if (in->arrayName) g_free(in->arrayName);
  if (in->obsName)   g_free(in->obsName);
  g_free(in);
  return NULL;
} /* end ObitSDMDataKillAntArray */

/**
 * Creates and fills a Source array
 * There is one entry per Spectral window
 * \param in   ASDM object to use
 * \return the new structure, NULL on error, 
 *         delete using ObitSDMDataKillSourceArray
 */
ASDMSourceArray* ObitSDMDataGetSourceArray (ObitSDMData *in)
{ 
  ASDMSourceArray* out=NULL;
  olong i, iSource, iField, sourceId, num;

  out = g_malloc0(sizeof(ASDMSourceArray));

  /* Assume no more sources than total Fields defined */
  out->sou  = g_malloc0(in->SourceTab->nrows*sizeof(ASDMSourceArrayEntry));
  out->nsou = 0;

  /* Loop over Sources */
  for (iSource=0; iSource<in->SourceTab->nrows; iSource++) {

    /* Create */
    out->sou[iSource] = g_malloc0(sizeof(ASDMSourceArrayEntry));
    out->nsou++;

    /* Copy source info */
    out->sou[iSource]->sourceId        = in->SourceTab->rows[iSource]->sourceId;
    out->sou[iSource]->sourceNo        = in->SourceTab->rows[iSource]->sourceNo;
    out->sou[iSource]->sourceName      = g_strdup(in->SourceTab->rows[iSource]->sourceName);
    out->sou[iSource]->timeInterval    = g_malloc0(2*sizeof(odouble));
    out->sou[iSource]->timeInterval[0] = in->SourceTab->rows[iSource]->timeInterval[0];
    out->sou[iSource]->timeInterval[1] = in->SourceTab->rows[iSource]->timeInterval[1];
    out->sou[iSource]->direction       = g_malloc0(2*sizeof(odouble));	       
    out->sou[iSource]->direction[0]    = in->SourceTab->rows[iSource]->direction[0];
    out->sou[iSource]->direction[1]    = in->SourceTab->rows[iSource]->direction[1];
    out->sou[iSource]->properMotion    = g_malloc0(2*sizeof(odouble));
    out->sou[iSource]->properMotion[0] = in->SourceTab->rows[iSource]->properMotion[0];
    out->sou[iSource]->properMotion[1] = in->SourceTab->rows[iSource]->properMotion[1];
    out->sou[iSource]->numLines        = in->SourceTab->rows[iSource]->numLines;
    num                                = out->sou[iSource]->numLines;
    out->sou[iSource]->restFrequency   = g_malloc0(num*sizeof(odouble));
    for (i=0; i<num; i++) out->sou[iSource]->restFrequency[i] = in->SourceTab->rows[iSource]->restFrequency[i];
    out->sou[iSource]->sysVel          = g_malloc0(num*sizeof(odouble));
    for (i=0; i<num; i++) out->sou[iSource]->sysVel[i] = in->SourceTab->rows[iSource]->sysVel[i];
    out->sou[iSource]->spectralWindowId = in->SourceTab->rows[iSource]->spectralWindowId;

    /* Find source in Fields */
    sourceId = out->sou[iSource]->sourceId;
    for (iField=0; iField<in->FieldTab->nrows; iField++) {
     if (in->FieldTab->rows[iField]->sourceId==sourceId) break;
       } /* end loop over fields */
    if (iField>=in->FieldTab->nrows) return NULL;
    
    /* Save code - ignore terminally stupid "NONE" */
    if (strcmp("NONE", in->FieldTab->rows[iField]->code) && 
	strcmp("none", in->FieldTab->rows[iField]->code))
	out->sou[iSource]->code = strdup(in->FieldTab->rows[iField]->code);
    else out->sou[iSource]->code = strdup(" ");
 } /* end loop over sources */

  return out;
} /* end ObitSDMDataGetSourceArray */

/**
 * Delete Source/field  info
 * \param in  Structure to delete
 * \return NULL pointer
 */
ASDMSourceArray* ObitSDMDataKillSourceArray (ASDMSourceArray *in)
{ 
  olong i;
  if (in==NULL) return NULL;
  
  /* delete entries */
  if (in->sou) {
    for (i=0; i<in->nsou; i++) {
      if (in->sou[i]) {
	if (in->sou[i]->timeInterval)  g_free(in->sou[i]->timeInterval);
	if (in->sou[i]->code)          g_free(in->sou[i]->code);
	if (in->sou[i]->direction)     g_free(in->sou[i]->direction);
	if (in->sou[i]->properMotion)  g_free(in->sou[i]->properMotion);
	if (in->sou[i]->sourceName)    g_free(in->sou[i]->sourceName);
	if (in->sou[i]->restFrequency) g_free(in->sou[i]->restFrequency);
	if (in->sou[i]->sysVel)        g_free(in->sou[i]->sysVel);
	g_free(in->sou[i]);
      }
    }
    g_free(in->sou);
  }
  g_free(in);
  return NULL;
} /* end ObitSDMDataKillSourceArray */

/**  Look up  Band code enum
 * \param  code  Band code to lookup
 * \return value
 */
ObitASDMBand ObitSDMDataBand2Band (gchar *code)
{
  ObitASDMBand out = ASDMBand_Any;
 
  if (!strncmp (code, "Any", 3)) return ASDMBand_Any;
  if (!strncmp (code, "4",   1)) return ASDMBand_4;
  if (!strncmp (code, "P",   1)) return ASDMBand_P;
  if (!strncmp (code, "L",   1)) return ASDMBand_L;
  if (!strncmp (code, "S",   1)) return ASDMBand_S;
  if (!strncmp (code, "C",   1)) return ASDMBand_C;
  if (!strncmp (code, "X",   1)) return ASDMBand_X;
  if (!strncmp (code, "Ku",  2)) return ASDMBand_Ku;
  if (!strncmp (code, "Ka",  2)) return ASDMBand_Ka;
  if (!strncmp (code, "K",   1)) return ASDMBand_K;
  if (!strncmp (code, "Q",   1)) return ASDMBand_Q;
  if (!strncmp (code, "A3",  1)) return ASDMBand_A3;
  if (!strncmp (code, "A4",  1)) return ASDMBand_A4;
  if (!strncmp (code, "A5",  1)) return ASDMBand_A5;
  if (!strncmp (code, "A6",  1)) return ASDMBand_A6;
  if (!strncmp (code, "A7",  1)) return ASDMBand_A7;
  if (!strncmp (code, "A8",  1)) return ASDMBand_A8;
  if (!strncmp (code, "A9",  1)) return ASDMBand_A9;
  if (!strncmp (code, "A10", 1)) return ASDMBand_A10;
  if (!strncmp (code, "A11", 1)) return ASDMBand_A11;
  return out;
} /* end ObitSDMDataBand2Band */

/**  Look up  Band code enum from frequency (Hz)
 * \param  freq  Frequency (Hz) to lookup
 * \return value, ASDMBand_Any if confused
 */
ObitASDMBand ObitSDMDataFreq2Band (odouble freq)
{
  ObitASDMBand out = ASDMBand_Any;
 
  if (freq<100.0e6) return ASDMBand_4;
  if (freq<900.0e6) return ASDMBand_P;
  if (freq<2.0e9)   return ASDMBand_L;
  if (freq<3.7e9)   return ASDMBand_S;
  if (freq<7.5e9)   return ASDMBand_C;
  if (freq<12.0e9)  return ASDMBand_X;
  if (freq<18.0e9)  return ASDMBand_Ku;
  if (freq<26.5e9)  return ASDMBand_K;
  if (freq<40.0e9)  return ASDMBand_Ka;
  if (freq<50.0e9)  return ASDMBand_Q;
  if (freq<117.0e9) return ASDMBand_A3;
  if (freq<163.0e9) return ASDMBand_A4;
  if (freq<211.0e9) return ASDMBand_A5;
  if (freq<275.0e9) return ASDMBand_A6;
  if (freq<375.0e9) return ASDMBand_A7;
  if (freq<510.0e9) return ASDMBand_A8;
  if (freq<730.0e9) return ASDMBand_A9;
  if (freq<960.0e9) return ASDMBand_A10;
  if (freq<2000.0e9) return ASDMBand_A11;
  return out;
} /* end ObitSDMDataFreq2Band */

/**  Find first selected Main table row, allows defaults 
 * \param  in      ASDM object
 * \param  selChan Number of selected channels [def 0]
 * \param  selIF   Number of selected IFs (spectral windows)[def 0]
 * \param  selBand   Selected band [def ASDMBand_Any]
 * \param  selConfig Selected configID, >=0 overrides selIF, selBand
 * \return 0-rel main row number, -1=> problem.
 */
olong ObitASDSelScan(ObitSDMData *in, olong selChan, olong selIF, 
		     ObitASDMBand selBand, olong selConfig)
{
  olong out = 1;
  olong configDescriptionId, dataDescriptionId, spectralWindowId;
  olong iMain, iConfig, iDD, jSW, jDD, numDD, i;
  gboolean doConfigId = selConfig>=0;
 
  /* Loop over scans in Main table */
  for (iMain=0; iMain<in->MainTab->nrows; iMain++) {
    out = iMain;  /* In case this is OK */
    
    /* Find entry  in configDescription table */
    configDescriptionId = in->MainTab->rows[iMain]->configDescriptionId;

    /* Specified config ID? */
    if (doConfigId && (selConfig==configDescriptionId)) return out;
    if (doConfigId) continue;  /* Others tests don't matter */

    for (iConfig=0; iConfig<in->ConfigDescriptionTab->nrows; iConfig++) {
      if (in->ConfigDescriptionTab->rows[iConfig]->configDescriptionId==configDescriptionId) break;
    }
    if (iConfig>=in->ConfigDescriptionTab->nrows) return -1;
    
      /* Count number of data descriptions = number of spectral windows */
      numDD = 0;
      i = 0;
      while(in->ConfigDescriptionTab->rows[iConfig]->dataDescriptionId[i++]>=0) {numDD++;}

    /* Reference frequency - get from Spectral window on any dataDescription */
    /* Find Data description */
    for (iDD=0; iDD<in->ConfigDescriptionTab->rows[iConfig]->numDataDescription; iDD++) {
      dataDescriptionId = in->ConfigDescriptionTab->rows[iConfig]->dataDescriptionId[iDD];
      for (jDD=0; jDD<in->DataDescriptionTab->nrows; jDD++) {
	if (in->DataDescriptionTab->rows[jDD]->dataDescriptionId==dataDescriptionId) break;
      }
      if (jDD>=in->DataDescriptionTab->nrows) return -1;
      
      /* Find spectralWindow */
      spectralWindowId = in->DataDescriptionTab->rows[jDD]->spectralWindowId;
      for (jSW=0; jSW<in->SpectralWindowTab->nrows; jSW++) {
	if (in->SpectralWindowTab->rows[jSW]->spectralWindowId==spectralWindowId) break;
      }
      if (jSW>=in->SpectralWindowTab->nrows) return -1;
      /* Finally, is this one a match? */
      if (((ObitSDMDataFreq2Band (in->SpectralWindowTab->rows[jSW]->refFreq)==selBand) || 
	   (selBand==ASDMBand_Any)) && 
	  ((in->SpectralWindowTab->rows[jSW]->numChan==selChan) || (selChan<=0)) && 
	  ((numDD==selIF) || (selIF<=0))) return out;
    } /* End loop over data descriptions */
  } /* End loop over main table */

  return -1;  /* Must not have found a match */
} /* end ObitASDSelScan */

/**
 * Give all sources with the same name the same source number
 * \param in  Structure with SourceTab to fix
 */
void ObitSDMSourceTabFix (ObitSDMData *in)
{ 
  olong i, j;
  ASDMSourceRow *irow, *jrow;

  if (in==NULL) return;
  /* fix entries */  
  for (i=0; i<in->SourceTab->nrows; i++) {
    if (in->SourceTab->rows[i]) {
      irow = in->SourceTab->rows[i];
	for (j=i+1; j<in->SourceTab->nrows; j++) {
	  if (in->SourceTab->rows[j]) {
	    jrow = in->SourceTab->rows[j];
	    if (!strcmp(irow->sourceName, jrow->sourceName))
	      jrow->sourceNo = irow->sourceNo;
	  }
	}
    }
  }
  return;
} /* end ObitSDMSourceTabFix */

/**
 * Give all sources with the same name and calcode the same source number
 * \param in  Structure with SourceTab to fix
 */
void ObitSDMSourceTabFixCode (ObitSDMData *in)
{ 
  olong i, j;
  ASDMSourceRow *irow, *jrow;

  if (in==NULL) return;
  /* fix entries */  
  for (i=0; i<in->SourceTab->nrows; i++) {
    if (in->SourceTab->rows[i]) {
      irow = in->SourceTab->rows[i];
	for (j=i+1; j<in->SourceTab->nrows; j++) {
	  if (in->SourceTab->rows[j]) {
	    jrow = in->SourceTab->rows[j];
	    if (!strcmp(irow->sourceName, jrow->sourceName) &&
		!strcmp(irow->code, jrow->code)) 
	      jrow->sourceNo = irow->sourceNo;
	  }
	}
    }
  }
  return;
} /* end ObitSDMSourceTabFixCode */

/**   Select Scan by source code
 * \param  in      ASDM object
 * \param  iMain   Main table number of scan in question
 * \param  selCode Desired Source code:
 *        "*" => Anything other than "NONE" or "none", "NONE"=>"NONE",
 *        No characters or starts with blank => any scan
 *        otherwise if any character matches the first character 
 *        in Field Table code, the scan is selected
 * \return True if scan desired
 */
gboolean ObitSDMDataSelCode  (ObitSDMData *in, olong iMain, gchar *selCode)
{
  gboolean want=TRUE;
  gboolean noCode;
  olong fieldId, iField, i, n;

  /* Checking anything? */
  if ((strlen(selCode)==0) || (selCode[0]==' ')) return want;

  /* Find Field table row */
  fieldId = in->MainTab->rows[iMain]->fieldId;
  for (iField=0; iField<in->FieldTab->nrows; iField++) {
    if (in->FieldTab->rows[iField]->fieldId==fieldId) break;
  }
  if (iField>=in->FieldTab->nrows) return want; /* Bother - not found */

  /* Selections */
  noCode = !strncmp(in->FieldTab->rows[iField]->code, "NONE", 4) ||
    !strncmp(in->FieldTab->rows[iField]->code, "none", 4);
  if ((selCode[0]=='*') && (!noCode)) return TRUE;
  if ((selCode[0]=='*') && (noCode))  return FALSE;
  if ((!strncmp(selCode, "NONE", 4)) && noCode) return TRUE;
  if ((!strncmp(selCode, "NONE", 4)) && 
      (in->FieldTab->rows[iField]->code[0]!=' '))  return FALSE;

  n    = (strlen(selCode));
  want = FALSE;
  for (i=0; i<n; i++) {
    if (selCode[i]== ' ') continue;
    if (selCode[i]==in->FieldTab->rows[iField]->code[0]) return TRUE;
  }

  return want;
} /* end ObitSDMDataSelCode */

/**   
 * Look up pointing for an antenna at a given time (JD) in PointingTab
 * \param  in      ASDM object
 * \param  JD      Time as JD in question
 * \param  ant     Antenna Id
 * \param  err     Obit error/message stack
 * \return pointer to selected table row entry or NULL if failed.
 */
ASDMPointingRow* ObitSDMDataPointingLookup  (ObitSDMData *in, odouble JD, olong ant, 
					     ObitErr *err)
{
  ASDMPointingRow* out=NULL;
  olong i, itab;
  odouble *timeI;
  gchar *routine = "PointingLookup";

  if (err->error) return out;

  /* Find matching entry */
  itab = -1;
  for (i=0; i<in->PointingTab->nrows; i++) {
    timeI = in->PointingTab->rows[i]->timeInterval;
    if ((in->PointingTab->rows[i]->antennaId==ant) &&
	(timeI[0]-0.5*timeI[1]<=JD) && (timeI[0]+0.5*timeI[1]>=JD))
      {itab=i; break;}
  } /* End loop looking for table */
  /* Find it? */
  Obit_retval_if_fail((itab>=0), err, out,
		      "%s: Failed to find pointing for ant %d JD %lf", 
		      routine, ant, JD);
  out = in->PointingTab->rows[itab];
  return out;
} /* end ObitSDMDataPointingLookup */

/**
 * Generate an SN table from data in an CalWVR SDM table
 * Delay calculated using:
 * T=270 K and eq. from Carilli, Carlstrom & Holdaway, Synthesis Imaging II
 * Output SN table will be associated with inUV.
 * \param inUV     Input uv data. 
 * \param SDM      ASDM srtucture
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created ObitTableSN
 */
ObitTableSN* ObitSDMDataWVR2SN (ObitUV *inUV, ObitSDMData *SDM,
				ObitErr *err)
{
  /* Speed of light */
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif
  /* CI = 1/speed of light */
#ifndef CI
#define CI 1.0 / VELIGHT
#endif
  ObitIOCode retCode;
  ObitTableSN    *outSN   = NULL;
  ObitTableSNRow *SNRow   = NULL;
  ASDMCalWVRRow  *WVRRow  = NULL;
  ASDMAntennaArray *antArray= NULL;
  olong iant, numPol, numIF, iif, mainRow=0;
  olong ver, iSNRow=0, iWVRRow;
  odouble delay, phase;
  gchar *tname;
  gchar *routine = "ObitSDMCalWVR2SN";
  
  /* error checks */
  if (err->error) return outSN;
  g_assert (ObitUVIsA(inUV));

  /* Must have a CalWVRTable */
  if ((SDM->CalWVRTab==NULL) || (SDM->CalWVRTab->nrows<1)) return outSN;
  
  /* create output table  */
  if (inUV->myDesc->jlocs>=0) numPol = MIN (2, inUV->myDesc->inaxes[inUV->myDesc->jlocs]);
  else numPol = 1;
  if (inUV->myDesc->jlocif>=0) numIF  = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
  else numIF = 1;
  if (outSN==NULL) {
    tname = g_strconcat ("Calibration for: ",inUV->name, NULL);
    ver = 0;
    outSN = newObitTableSNValue (tname, (ObitData*)inUV, &ver, OBIT_IO_WriteOnly,  
				 numPol, numIF, err);
    g_free (tname);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, outSN);
  }

  /* Tell about it */
  Obit_log_error(err, OBIT_InfoErr, "Converting CalWVR data to SN table %d", ver);
  
  /* Open SN table for write */
  retCode = ObitTableSNOpen (outSN, OBIT_IO_WriteOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, outSN->name, outSN);

  /* create row structure */
  SNRow = newObitTableSNRow(outSN);

   /* Attach row to output buffer */
  ObitTableSNSetRow (outSN, SNRow, err);
  if (err->error) Obit_traceback_val (err, routine, outSN->name, outSN);

  outSN->numAnt = 1;

 /* Initialize SN Row structure */
  SNRow->Time   = 0.0;
  SNRow->TimeI  = 0.0;
  SNRow->SourID = 0;
  SNRow->antNo  = 0;
  SNRow->SubA   = 0;
  SNRow->FreqID = 0;
  SNRow->IFR    = 0.0;
  SNRow->NodeNo = 1;
  SNRow->MBDelay1 = 0.0;
  SNRow->MBDelay2 = 0.0;
  SNRow->status = 1;
  for (iif = 0; iif < numIF; iif++) {
    SNRow->Real1[iif]   = 0.0;
    SNRow->Imag1[iif]   = 0.0;
    SNRow->Delay1[iif]  = 0.0;
    SNRow->Rate1[iif]   = 0.0;
    SNRow->Weight1[iif] = 0.0;
    SNRow->RefAnt1[iif] = 1;
    if (numPol>1) { /* Second polarization */
      SNRow->Real2[iif]   = 0.0;
      SNRow->Imag2[iif]   = 0.0;
      SNRow->Delay2[iif]  = 0.0;
      SNRow->Rate2[iif]   = 0.0;
      SNRow->Weight2[iif] = 0.0;
      SNRow->RefAnt2[iif] = 1;
    }
  }
  
  /* Get antenna information */
  antArray = ObitSDMDataGetAntArray (SDM, mainRow);

  /* Frequency info */
  if (inUV->myDesc->freqIF==NULL) {
    ObitUVGetFreq (inUV, err);
    if (err->error) goto cleanup;
  }

  /* Loop over WVR table converting */
  for (iWVRRow=0; iWVRRow<SDM->CalWVRTab->nrows; iWVRRow++) {
    WVRRow = SDM->CalWVRTab->rows[iWVRRow];
    /* Set antenna invariant values */
    SNRow->Time   = 0.5*(WVRRow->endValidTime + WVRRow->startValidTime) - antArray->refJD;
    SNRow->TimeI  = WVRRow->endValidTime - WVRRow->startValidTime;
    
    /* Lookup antenna number */
    for (iant=0; iant<antArray->nants; iant++) {
      if (!strcmp(WVRRow->antennaName, antArray->ants[iant]->antName)) {
	SNRow->antNo  = antArray->ants[iant]->antennaNo+1;
	break;
      }
    } /* end looking up antenna number */
    
    /* convert WVR model to SN */
    for (iif = 0; iif < numIF; iif++) {
      delay = -6.3 * WVRRow->water * CI;
      phase = 2*G_PI*inUV->myDesc->freqIF[iif]*delay;
      SNRow->Real1[iif]   = (ofloat)cos(phase);
      SNRow->Imag1[iif]   = (ofloat)sin(phase);
      SNRow->Delay1[iif]  = (ofloat)delay;
      SNRow->Weight1[iif] = 1.0;
      if (numPol>1) { /* Second polarization */
	SNRow->Real2[iif]   = (ofloat)cos(phase);
	SNRow->Imag2[iif]   = (ofloat)sin(phase);
	SNRow->Delay2[iif]  = (ofloat)delay;
	SNRow->Weight2[iif] = 1.0;
     }
    } /* end IF loop */

    /* write it */
    retCode = ObitTableSNWriteRow (outSN, iSNRow, SNRow, err);
    if (err->error) goto cleanup;
    
  } /* end loop over WVR Table */

  /* Close table */
  retCode = ObitTableSNClose (outSN, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* deallocate structures */
 cleanup:
  SNRow = ObitTableSNRowUnref(SNRow);
  antArray = ObitSDMDataKillAntArray (antArray);
  if (err->error) Obit_traceback_val (err, routine, outSN->name, outSN);

  return outSN;
} /* end ObitSDMDataWVR2SN  */

/**
 * Lookup further gain entries for a given antenna/time
 * \param calAtmosphereTab Atmospheric data table
 * \param UVDesc           UV descriptor
 * \param iAtmRow          Master table row
 * \param jAtmRow          [in/out] Next table row
 * \param BBIndex          Lookup table for IF corresponding to baseband number
 * \param gain1            [out] gains for first receptor (poln), 1/IF
 * \param gain2            [out] gains for second receptor (poln)
 * \return TRUE if there may be more entries
 */
static gboolean FindMoreAtmData(ASDMcalAtmosphereTable *calAtmosphereTab, ObitUVDesc *UVDesc, 
				olong iAtmRow, olong *jAtmRow, olong *BBIndex,
				ofloat *gain1, ofloat *gain2)
{
  gboolean more=FALSE;
  olong i, iif, numIF, numPol, iRow, oldRow=*jAtmRow;
  ofloat fblank = ObitMagicF();

  /* Finished */
  if (*jAtmRow>=calAtmosphereTab->nrows) return more;

  /* Numbers of things */
  if (UVDesc->jlocs>=0) numPol = MIN (2, UVDesc->inaxes[UVDesc->jlocs]);
  else numPol = 1;
  if (UVDesc->jlocif>=0) numIF  = UVDesc->inaxes[UVDesc->jlocif];
  else numIF = 1;

  /* Init on first */
  if (iAtmRow==oldRow) {
    for (i=0; i<numIF; i++) {
      gain1[i] = gain2[i] = fblank;
    }
    iRow = iAtmRow;
  } else {  /* Subsequent call - check if same antenna */
    if (!strcmp(calAtmosphereTab->rows[iAtmRow]->antennaName, 
		 calAtmosphereTab->rows[*jAtmRow]->antennaName)) {
      iRow = *jAtmRow;
    } else { /* Next antenna - stop looking */
      return FALSE;
    }
  }
  
  /* Find IF - set gains */
  iif = BBIndex[calAtmosphereTab->rows[iRow]->basebandId];
  if (iif>=0) {
    gain1[iif] = sqrt(calAtmosphereTab->rows[iRow]->tSys[0]);
    if (numPol>1) gain2[iif] = sqrt(calAtmosphereTab->rows[iRow]->tSys[1]);
  }
    
  /* More to find? */
  more = FALSE;
  for (i=0; i<numIF; i++) {
    if (gain1[i]==fblank) more = TRUE;
    if ((numPol>1) && (gain2[i]==fblank)) more = TRUE;
    if (more) break;
  }
  *jAtmRow = oldRow+1; /* next */

  return more;
} /* end FindMoreAtmData */

/**
 * Generate an SN table from data in an CalAtmosphere SDM table
 * Output SN table will be associated with inUV.
 * \param inUV       Input uv data. 
 * \param SDM        ASDM srtucture
 * \param SpWinArray Spectral window array used for UV data.
 * \param err        Error stack, returns if not empty.
 * \return Pointer to the newly created ObitTableSN
 */
ObitTableSN* ObitSDMDataAtm2SN (ObitUV *inUV, ObitSDMData *SDM,
				ASDMSpectralWindowArray* SpWinArray, 
				ObitErr *err)
{
  ObitIOCode retCode;
  ObitTableSN    *outSN   = NULL;
  ObitTableSNRow *SNRow   = NULL;
  ASDMcalAtmosphereRow  *AtmRow  = NULL;
  ASDMAntennaArray *antArray= NULL;
  olong i, iant, numPol, numIF, iif, mainRow=0;
  olong ver, iSNRow=0, iAtmRow, jAtmRow;
  olong BBIndex[33];
  ofloat *gain1=NULL, *gain2=NULL;
  ofloat fblank = ObitMagicF();
  gchar *tname;
  gchar *routine = "ObitSDMCalAtm2SN";
  
  /* error checks */
  if (err->error) return outSN;
  g_assert (ObitUVIsA(inUV));

  /* Must have a CalAtmTable */
  if ((SDM->calAtmosphereTab==NULL) || (SDM->calAtmosphereTab->nrows<1)) return outSN;

  /* Baseband lookup table 0-rel IF number */
  for (i=0; i<33; i++) BBIndex[i] = -1;
  for (i=0; i<SpWinArray->nwinds; i++) BBIndex[SpWinArray->winds[i]->basebandNum] = i;
  
  /* create output table  */
  if (inUV->myDesc->jlocs>=0) numPol = MIN (2, inUV->myDesc->inaxes[inUV->myDesc->jlocs]);
  else numPol = 1;
  if (inUV->myDesc->jlocif>=0) numIF  = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
  else numIF = 1;
  if (outSN==NULL) {
    tname = g_strconcat ("Calibration for: ",inUV->name, NULL);
    ver = 0;
    outSN = newObitTableSNValue (tname, (ObitData*)inUV, &ver, OBIT_IO_WriteOnly,  
				 numPol, numIF, err);
    g_free (tname);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, outSN);
  }

  /* Tell about it */
  Obit_log_error(err, OBIT_InfoErr, "Converting CalAtmosphere data to SN table %d", ver);
  
  /* Open SN table for write */
  retCode = ObitTableSNOpen (outSN, OBIT_IO_WriteOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, outSN->name, outSN);

  /* create row structure */
  SNRow = newObitTableSNRow(outSN);

   /* Attach row to output buffer */
  ObitTableSNSetRow (outSN, SNRow, err);
  if (err->error) Obit_traceback_val (err, routine, outSN->name, outSN);

  outSN->numAnt = 1;

 /* Initialize SN Row structure */
  SNRow->Time   = 0.0;
  SNRow->TimeI  = 0.0;
  SNRow->SourID = 0;
  SNRow->antNo  = 0;
  SNRow->SubA   = 0;
  SNRow->FreqID = 0;
  SNRow->IFR    = 0.0;
  SNRow->NodeNo = 1;
  SNRow->MBDelay1 = 0.0;
  SNRow->MBDelay2 = 0.0;
  SNRow->status = 1;
  for (iif = 0; iif < numIF; iif++) {
    SNRow->Real1[iif]   = 0.0;
    SNRow->Imag1[iif]   = 0.0;
    SNRow->Delay1[iif]  = 0.0;
    SNRow->Rate1[iif]   = 0.0;
    SNRow->Weight1[iif] = 0.0;
    SNRow->RefAnt1[iif] = 1;
    if (numPol>1) { /* Second polarization */
      SNRow->Real2[iif]   = 0.0;
      SNRow->Imag2[iif]   = 0.0;
      SNRow->Delay2[iif]  = 0.0;
      SNRow->Rate2[iif]   = 0.0;
      SNRow->Weight2[iif] = 0.0;
      SNRow->RefAnt2[iif] = 1;
    }
  }
  
  /* Get antenna information */
  antArray = ObitSDMDataGetAntArray (SDM, mainRow);

  /* Frequency info */
  if (inUV->myDesc->freqIF==NULL) {
    ObitUVGetFreq (inUV, err);
    if (err->error) goto cleanup;
  }

  /* Work arrays */
  gain1 = g_malloc0(numIF*sizeof(ofloat));
  gain2 = g_malloc0(numIF*sizeof(ofloat));

  /* Loop over Atmosphere table converting */
  for (iAtmRow=0; iAtmRow<SDM->calAtmosphereTab->nrows; ) {
    AtmRow = SDM->calAtmosphereTab->rows[iAtmRow];
    
    /* Routine to find further entries for this antenna and return array of gains */
    jAtmRow = iAtmRow;
    while (FindMoreAtmData(SDM->calAtmosphereTab, inUV->myDesc, iAtmRow, &jAtmRow, 
			   BBIndex, gain1, gain2));
    
    /* Set antenna invariant values */
    SNRow->Time   = AtmRow->startValidTime - antArray->refJD;
    SNRow->TimeI  = 0.5 * (AtmRow->endValidTime - AtmRow->startValidTime);
    
    /* Lookup antenna number */
    for (iant=0; iant<antArray->nants; iant++) {
      if (!strcmp(AtmRow->antennaName, antArray->ants[iant]->antName)) {
	SNRow->antNo  = antArray->ants[iant]->antennaNo+1;
	break;
      }
    } /* end looking up antenna number */
    
    /* convert Atm model to SN */
    for (iif = 0; iif < numIF; iif++) {
      if (gain1[iif]!=fblank) {
	SNRow->Real1[iif]   = gain1[iif];
	SNRow->Imag1[iif]   = 0.0;
	SNRow->Weight1[iif] = 1.0;
      } else {
	SNRow->Real1[iif]   = fblank;
	SNRow->Imag1[iif]   = fblank;
	SNRow->Weight1[iif] = 0.0;
      }
      if (numPol>1) { /* Second polarization */
	if (gain2[iif]!=fblank) {
	  SNRow->Real2[iif]   = gain2[iif];
	  SNRow->Imag2[iif]   = 0.0;
	  SNRow->Weight2[iif] = 1.0;
	} else {
	  SNRow->Real2[iif]   = fblank;
	  SNRow->Imag2[iif]   = fblank;
	  SNRow->Weight2[iif] = 0.0;
	}
      }
    } /* end IF loop */
    
    /* write it - first start time */
    retCode = ObitTableSNWriteRow (outSN, iSNRow, SNRow, err);
    /* then end time - no they're bogus
       SNRow->Time = AtmRow->endValidTime - antArray->refJD;
       retCode = ObitTableSNWriteRow (outSN, iSNRow, SNRow, err); */
    if (err->error) goto cleanup;
    iAtmRow = jAtmRow;  /* skip ones done */
  } /* end loop over Atm Table */

  /* Close table */
  retCode = ObitTableSNClose (outSN, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* deallocate structures */
 cleanup:
  if (gain1) g_free(gain1);
  if (gain2) g_free(gain2);
  SNRow = ObitTableSNRowUnref(SNRow);
  antArray = ObitSDMDataKillAntArray (antArray);
  if (err->error) Obit_traceback_val (err, routine, outSN->name, outSN);

  return outSN;
} /* end ObitSDMDataAtm2SN  */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSDMDataClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSDMDataClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSDMDataClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSDMDataClassInfoDefFn (gpointer inClass)
{
  ObitSDMDataClassInfo *theClass = (ObitSDMDataClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSDMDataClassInit;
  theClass->newObit       = (newObitFP)newObitSDMData;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSDMDataClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSDMDataGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitSDMDataCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitSDMDataClear;
  theClass->ObitInit      = (ObitInitFP)ObitSDMDataInit;
  theClass->ObitSDMDataCreate = (ObitSDMDataCreateFP)ObitSDMDataCreate;

} /* end ObitSDMDataClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitSDMDataInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSDMData *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->DataRoot             = NULL;
  in->ASDMTab              = NULL;
  in->MainTab              = NULL;
  in->AntennaTab           = NULL;
  in->calAtmosphereTab     = NULL;
  in->calDataTab           = NULL;
  in->calDeviceTab         = NULL;
  in->calPointingTab       = NULL;
  in->CalReductionTab      = NULL;
  in->CalWVRTab            = NULL;
  in->ConfigDescriptionTab = NULL;
  in->CorrelatorModeTab    = NULL;
  in->DataDescriptionTab   = NULL;
  in->DopplerTab           = NULL;
  in->ExecBlockTab         = NULL;
  in->FeedTab              = NULL;
  in->FieldTab             = NULL;
  in->FlagTab              = NULL;
  in->PointingTab          = NULL;
  in->PointingModelTab     = NULL;
  in->PolarizationTab      = NULL;
  in->ProcessorTab         = NULL;
  in->ReceiverTab          = NULL;
  in->SBSummaryTab         = NULL;
  in->ScanTab              = NULL;
  in->SourceTab            = NULL;
  in->SpectralWindowTab    = NULL;
  in->StateTab             = NULL;
  in->StationTab           = NULL;
  in->SubscanTab           = NULL;
  in->SwitchCycleTab       = NULL;
  in->SysCalTab            = NULL;
  in->SysPowerTab          = NULL;
  in->WeatherTab           = NULL;
  in->line                 = NULL;

} /* end ObitSDMDataInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitSDMData* cast to an Obit*.
 */
void ObitSDMDataClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSDMData *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->DataRoot) g_free(in->DataRoot); in->DataRoot = NULL;
  if (in->line) g_free(in->line);  in->line = NULL;
  in->ASDMTab              = KillASDMTable(in->ASDMTab);
  in->MainTab              = KillASDMMainTable(in->MainTab);
  in->AntennaTab           = KillASDMAntennaTable(in->AntennaTab);
  in->calAtmosphereTab     = KillASDMcalAtmosphereTable(in->calAtmosphereTab);
  in->calDataTab           = KillASDMcalDataTable(in->calDataTab);
  in->calDeviceTab         = KillASDMcalDeviceTable(in->calDeviceTab);
  in->calPointingTab       = KillASDMcalPointingTable(in->calPointingTab);
  in->CalReductionTab      = KillASDMCalReductionTable(in->CalReductionTab);
  in->CalWVRTab            = KillASDMCalWVRTable(in->CalWVRTab);
  in->ConfigDescriptionTab = KillASDMConfigDescriptionTable(in->ConfigDescriptionTab);
  in->CorrelatorModeTab    = KillASDMCorrelatorModeTable(in->CorrelatorModeTab);
  in->DataDescriptionTab   = KillASDMDataDescriptionTable(in->DataDescriptionTab);
  in->DopplerTab           = KillASDMDopplerTable(in->DopplerTab);
  in->ExecBlockTab         = KillASDMExecBlockTable(in->ExecBlockTab);
  in->FeedTab              = KillASDMFeedTable(in->FeedTab);
  in->FieldTab             = KillASDMFieldTable(in->FieldTab);
  in->FlagTab              = KillASDMFlagTable(in->FlagTab);
  in->PointingTab          = KillASDMPointingTable(in->PointingTab);
  in->PointingModelTab     = KillASDMPointingModelTable(in->PointingModelTab);
  in->PolarizationTab      = KillASDMPolarizationTable(in->PolarizationTab);
  in->ProcessorTab         = KillASDMProcessorTable(in->ProcessorTab);
  in->ReceiverTab          = KillASDMReceiverTable(in->ReceiverTab);
  in->SBSummaryTab         = KillASDMSBSummaryTable(in->SBSummaryTab);
  in->ScanTab              = KillASDMScanTable(in->ScanTab);
  in->SourceTab            = KillASDMSourceTable(in->SourceTab);
  in->SpectralWindowTab    = KillASDMSpectralWindowTable(in->SpectralWindowTab);
  in->StateTab             = KillASDMStateTable(in->StateTab);
  in->StationTab           = KillASDMStationTable(in->StationTab);
  in->SubscanTab           = KillASDMSubscanTable(in->SubscanTab);
  in->SwitchCycleTab       = KillASDMSwitchCycleTable(in->SwitchCycleTab);
  in->SysCalTab            = KillASDMSysCalTable(in->SysCalTab);
  in->SysPowerTab          = KillASDMSysPowerTable(in->SysPowerTab);
  in->WeatherTab           = KillASDMWeatherTable(in->WeatherTab);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSDMDataClear */

/* ASDM Tables routines */
/**  Parse 32 bit integer from XLM string 
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, 0 if problem
 */
static olong ASDMparse_int(gchar *string, olong maxChar, 
			   gchar *prior, gchar **next)
{
  olong out = 0;
  gchar *b;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  out = (olong)strtol(b, next, 10);
    
  return out;
} /* end ASDMparse_int */

/**  Parse 64 bit integer from XLM string 
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, 0 if problem
 */
static ollong ASDMparse_int64(gchar *string, olong maxChar, 
			     gchar *prior, gchar **next)
{
  ollong out = 0;
  gchar *b;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  out = (ollong)strtoll(b, next, 10);
    
  return out;
} /* end ASDMparse_int64 */

/**  Parse double from XLM string 
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, 0.0 if problem
 */
static odouble ASDMparse_dbl(gchar *string, olong maxChar, 
			     gchar *prior, gchar **next)
{
  odouble out = 0.0;
  gchar *b;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  out = (odouble)strtod(b, next);
    
  return out;
} /* end ASDMparse_dbl */

/**  Parse boolean from XLM string 
 * Either true or false
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, FALSE if problem
 */
static gboolean ASDMparse_boo(gchar *string, olong maxChar, 
			      gchar *prior, gchar **next)
{
  gboolean out = FALSE;
  gchar *b;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  if (!strncmp(b, "true",  4)) return TRUE;
  if (!strncmp(b, "false", 5)) return FALSE;
  *next = b;

  return out;
} /* end ASDMparse_boo */

/**  Parse unquoted string from XLM string 
 * All text from end of prior until next '<'
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 */
static gchar* ASDMparse_str(gchar *string, olong maxChar, 
			   gchar *prior, gchar **next)
{
  gchar *out = NULL;
  gchar *b;
  olong charLeft;
  olong i, n;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);

  /* count */
  charLeft = maxChar - (b-string);
  n = 0;
  for (i=0; i<charLeft; i++) {
    if (b[i]=='<') break;
    n++;
  }
  out = g_malloc(n+1);
  for (i=0; i<n; i++) out[i] = b[i]; out[i] = 0;
  *next = b + n;

  return out;
} /* end ASDMparse_str */

/**  Parse double quoted from XLM string 
 * All text from end of prior+1 until next '"'
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 */
static gchar* ASDMparse_quote_str(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next)
{
  gchar *out = NULL;
  gchar *b;
  olong charLeft;
  olong i, n;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  if (*b!='"') return out;  /* Make sure quote */
  b++;                      /* Skip quote */

  /* count */
  charLeft = maxChar - (b-string);
  n = 0;
  for (i=0; i<charLeft; i++) {
    if (b[i]=='"') break;
    n++;
  }
  out = g_malloc(n+1);
  for (i=0; i<n; i++) out[i] = b[i]; out[i] = 0;
  *next = b + n;

  return out;
} /* end ASDMparse_quote_str */

/**  Parse array of doubles from XML string 
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 */
static odouble* ASDMparse_dblarray(gchar *string, olong maxChar, 
				   gchar *prior, gchar **next)
{
  odouble *out = NULL;
  gchar *b;
  olong charLeft, ndim, naxis1=1, naxis2=1, naxis3=1, num;
  olong i;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);

  /* Get dimensionality - only can handle 2 */
  ndim = (olong)strtol(b, next, 10);
  b = *next;
  /* get number of values axis 1 */
  naxis1 = (olong)strtol(b, next, 10);
  b = *next;
  /* get number of values axis 2 */
  if (ndim==2) {
    naxis2 = (olong)strtol(b, next, 10);
    b = *next;
  }
  /* get number of values axis 3 */
  if (ndim==3) {
    naxis2 = (olong)strtol(b, next, 10);
    b = *next;
    naxis3 = (olong)strtol(b, next, 10);
    b = *next;
  }
  num = naxis1*naxis2*naxis3;
  out = g_malloc0(MAX(1,num)*sizeof(odouble));

  /* Loop parsing */
  for (i=0; i<num; i++) {
    charLeft = maxChar - (b-string);
    out[i] = (odouble)strtod(b, next);
    b = *next;
  } /* end loop parsing values */

  return out;
} /* end ASDMparse_dblarray */

/**  Parse array of floats from XML string 
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 */
static ofloat* ASDMparse_fltarray(gchar *string, olong maxChar, 
				   gchar *prior, gchar **next)
{
  ofloat *out = NULL;
  gchar *b;
  olong charLeft, ndim, naxis1=1, naxis2=1, num;
  olong i;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);

  /* Get dimensionality - only can handle 2 */
  ndim = (olong)strtol(b, next, 10);
  b = *next;
  /* get number of values axis 1 */
  naxis1 = (olong)strtol(b, next, 10);
  b = *next;
  /* get number of values axis 2 */
  if (ndim==2) {
    naxis2 = (olong)strtol(b, next, 10);
    b = *next;
  }
  num = naxis1*naxis2;
  out = g_malloc0(MAX(1,num)*sizeof(ofloat));

  /* Loop parsing */
  for (i=0; i<num; i++) {
    charLeft = maxChar - (b-string);
    out[i] = (ofloat)strtod(b, next);
    b = *next;
  } /* end loop parsing values */

  return out;
} /* end ASDMparse_fltarray */

/**  Parse array of doubles from XML string 
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 */
static olong* ASDMparse_intarray(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next)
{
  olong *out = NULL;
  gchar *b;
  olong charLeft, ndim, naxis1=1, naxis2=1, num;
  olong i;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);

  /* Get dimensionality - only can handle 2 */
  ndim = (olong)strtol(b, next, 10);
  b = *next;
  /* get number of values */
  naxis1 = (olong)strtol(b, next, 10);
  b = *next;
  /* get number of values axis 2 */
  if (ndim==2) {
    naxis2 = (olong)strtol(b, next, 10);
    b = *next;
  }
  num = naxis1*naxis2;
  out = g_malloc0(MAX(1,num)*sizeof(olong));

  /* Loop parsing */
  for (i=0; i<num; i++) {
    charLeft = maxChar - (b-string);
    out[i] = (olong)strtol(b, next, 10);
    b = *next;
  } /* end loop parsing values */

  return out;
} /* end ASDMparse_intarray */

/**  Parse array of enums ending with the index from XML string 
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 *         one more entry than size, last = -999;
 */
static olong* ASDMparse_enumarray(gchar *string, olong maxChar, 
				  gchar *prior, gchar **next)
{
  olong *out = NULL;
  gchar *b;
  olong charLeft, ndim, naxis1=1, naxis2=1, num;
  olong i, j;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);

  /* Get dimensionality - only can handle 2 */
  ndim = (olong)strtol(b, next, 10);
  b = *next;
  /* get number of values */
  naxis1 = (olong)strtol(b, next, 10);
  b = *next;
  /* get number of values axis 2 */
  if (ndim==2) {
    naxis2 = (olong)strtol(b, next, 10);
    b = *next;
  }
  num = naxis1*naxis2;
  out = g_malloc0(MAX(1,(num+1))*sizeof(olong));

  /* Loop parsing the integer after each '-' */
  for (i=0; i<num; i++) {
    charLeft = maxChar - (b-string);
    for (j=0; j<charLeft; j++) {
      if (*(b++)=='_') break;
    }
    out[i] = (olong)strtol(b, next, 10);
    b = *next;
  }  /* end loop parsing values */

  out[num] = -999;     /* Mark end */
  return out;
} /* end ASDMparse_enumarray */

/**  Parse array of unquoted strings from XLM string 
 * Each string text from end of prior until next ' ' or '<'
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 */
static gchar** ASDMparse_strarray(gchar *string, olong maxChar, 
				  gchar *prior, gchar **next)
{
  gchar **out = NULL;
  gchar *b;
  olong charLeft, ndim, naxis1=1, naxis2=1, num;
  olong i, j, n;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);

   /* Get dimensionality - only can handle 2 */
  ndim = (olong)strtol(b, next, 10);
  b = *next;
  /* get number of values */
  naxis1 = (olong)strtol(b, next, 10);
  b = *next+1;
  /* get number of values axis 2 */
  if (ndim==2) {
    naxis2 = (olong)strtol(b, next, 10);
    b = *next+1;
  }
  num = naxis1*naxis2;
  out = g_malloc0(MAX(1,(num+1))*sizeof(gchar*));

  /* Loop over strings */
  for (j=0; j<num; j++) {
    /* count */
    charLeft = maxChar - (b-string);
    n = 0;
    for (i=0; i<charLeft; i++) {
      if ((b[i]=='<') || (b[i]==' '))break;
      n++;
    }
    out[j] = g_malloc(n+1);
    for (i=0; i<n; i++) out[j][i] = b[i]; out[j][i] = 0;
    b += n+1;
    *next = b;
  } /* end loop over strings */

  return out;
} /* end ASDMparse_strarray */

/**  Parse array of booleans from XLM string 
 * Values of "true" or "false" terminated by ' ' or '<'
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 */
static gboolean* ASDMparse_booarray(gchar *string, olong maxChar, 
				  gchar *prior, gchar **next)
{
  gboolean *out = NULL;
  gchar *b;
  olong charLeft, ndim, naxis1=1, naxis2=1, num;
  olong i, j, n;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);

   /* Get dimensionality - only can handle 2 */
  ndim = (olong)strtol(b, next, 10);
  b = *next;
  /* get number of values */
  naxis1 = (olong)strtol(b, next, 10);
  b = *next+1;
  /* get number of values axis 2 */
  if (ndim==2) {
    naxis2 = (olong)strtol(b, next, 10);
    b = *next+1;
  }
  num = naxis1*naxis2;
  out = g_malloc0(MAX(1,num)*sizeof(gboolean));

  /* Loop over strings */
  for (j=0; j<num; j++) {
    /* count */
    charLeft = maxChar - (b-string);
    n = 0;
    for (i=0; i<charLeft; i++) {
      if ((b[i]=='<') || (b[i]==' '))break;
      n++;
    }
    /* Only check first character */
    out[j] = b[0]=='t';
    b += n+1;
    *next = b;
  } /* end loop over strings */

  return out;
} /* end ASDMparse_booarray */

/**  Parse time from XLM string  
 * Read time as a MJD nanoseconds, return as JD
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, 0.0 if problem
 */
static odouble ASDMparse_time(gchar *string, olong maxChar, 
			      gchar *prior, gchar **next)
{
  odouble out = 0.0;
  long long temp;
  gchar *b;
  odouble mjdJD0=2400000.5; /* JD of beginning of MJD time */

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  temp = strtoll(b, next, 10);
  out = (odouble)((temp*1.0e-9)/86400.0) + mjdJD0;
    
  return out;
} /* end ASDMparse_time */

/**  Parse time from XLM string  
 * Read time range as a mjd nanoseconds, return as pair of JD
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, [0.0,1.0e20] if problem, g_freed when done
 */
static odouble* ASDMparse_timeRange(gchar *string, olong maxChar, 
				   gchar *prior, gchar **next)
{
  odouble *out;
  long long temp;
  gchar *b;
  odouble mjdJD0=2400000.5; /* JD of beginning of MJD time */

  out = g_malloc(2*sizeof(odouble));
  out[0] = 0.0; out[1] = 1.0e20;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  temp = strtoll(b, next, 10);
  out[0] = (odouble)((temp*1.0e-9)/86400.0) + mjdJD0;
  b = *next;
  temp = strtoll(b, next, 10);
  out[1] = (odouble)((temp*1.0e-9)/86400.0) + mjdJD0;
   
  return out;
} /* end ASDMparse_timeRange */

/**  Parse time interval from XLM string  
 * Read time interval in nanoseconds, return as days
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, 0.0 if problem
 */
static odouble ASDMparse_timeint(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next)
{
  odouble out = 0.0;
  long long temp;
  gchar *b;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  temp = strtoll(b, next, 10);
  out = (odouble)((temp*1.0e-9)/86400.0);
    
  return out;
} /* end ASDMparse_timeint */

/**  Look up axis name enum 
 * \param  string  String to lookup
 * \return value
 */
static ObitASDMAxisName LookupAxisName(gchar *name)
{
  ObitASDMAxisName out = 0;
 
  if (!strncmp (name, "TIM", 3)) return ASDMAxisName_TIM;
  if (!strncmp (name, "BAL", 3)) return ASDMAxisName_BAL;
  if (!strncmp (name, "ANT", 3)) return ASDMAxisName_ANT;
  if (!strncmp (name, "BAB", 3)) return ASDMAxisName_BAB;
  if (!strncmp (name, "SPW", 3)) return ASDMAxisName_SPW;
  if (!strncmp (name, "SIB", 3)) return ASDMAxisName_SIB;
  if (!strncmp (name, "SUB", 3)) return ASDMAxisName_SUB;
  if (!strncmp (name, "BIN", 3)) return ASDMAxisName_BIN;
  if (!strncmp (name, "APC", 3)) return ASDMAxisName_APC;
  if (!strncmp (name, "SPP", 3)) return ASDMAxisName_SPP;
  if (!strncmp (name, "POL", 3)) return ASDMAxisName_POL;
  if (!strncmp (name, "STO", 3)) return ASDMAxisName_STO;
  if (!strncmp (name, "HOL", 3)) return ASDMAxisName_HOL;
  return out;
} /* end LookupAxisName */

/**  Look up baseband name enum 
 * \param  string  String to look up
 * \return value
 */
static ObitASDMBasebandName LookupBasebandName(gchar *name)
{
  ObitASDMBasebandName out = 0;
 
  if (!strncmp (name, "NOBB", 4))   return ASDMBasebandName_NOBB;
  if (!strncmp (name, "BB_1", 4))   return ASDMBasebandName_BB_1;
  if (!strncmp (name, "BB_2", 4))   return ASDMBasebandName_BB_2;
  if (!strncmp (name, "BB_3", 4))   return ASDMBasebandName_BB_3;
  if (!strncmp (name, "BB_4", 4))   return ASDMBasebandName_BB_4;
  if (!strncmp (name, "BB_5", 4))   return ASDMBasebandName_BB_5;
  if (!strncmp (name, "BB_6", 4))   return ASDMBasebandName_BB_6;
  if (!strncmp (name, "BB_7", 4))   return ASDMBasebandName_BB_7;
  if (!strncmp (name, "BB_8", 4))   return ASDMBasebandName_BB_8;
  if (!strncmp (name, "BB_ALL", 6)) return ASDMBasebandName_BB_ALL;
  return out;
} /* end LookupBasebandName */

/**  Look up antenna make enum 
 * \param  string  String to look up
 * \return value
 */
static ObitASDMAntennaMake LookupAntennaMake(gchar *name)
{
  ObitASDMAntennaMake out = 0;
 
  if (!strncmp (name, "UNKNOWN", 7))          return  ASDMAnt_UNKNOWN;
  if (!strncmp (name, "MITSUBISHI_12_A", 15)) return ASDMAnt_MITSUBISHI_12_A;
  if (!strncmp (name, "MITSUBISHI_12_B", 15)) return ASDMAnt_MITSUBISHI_12_B;
  if (!strncmp (name, "VERTEX_12_ATF", 13))   return ASDMAnt_VERTEX_12_ATF;
  if (!strncmp (name, "AEM_12_ATF", 10))      return ASDMAnt_AEM_12_ATF;
  if (!strncmp (name, "VERTEX_12", 9))        return ASDMAnt_VERTEX_12;
  if (!strncmp (name, "AEM_12", 6))           return ASDMAnt_AEM_12;
  if (!strncmp (name, "IRAM_15", 7))          return ASDMAnt_IRAM_15;
  return out;
} /* end LookupAntennaMake */

/**  Look up polarization type enum 
 * \param  string  String to look up
 * \return value
 */
static ObitASDMPolnType LookupPolnType(gchar *name)
{
  ObitASDMPolnType out = 0;
 
  if (name[0]=='R') return ASDMPolnType_R;
  if (name[0]=='L') return ASDMPolnType_L;
  if (name[0]=='X') return ASDMPolnType_X;
  if (name[0]=='Y') return ASDMPolnType_Y;
  return out;
} /* end LookupPolnType */

/**  Look up  sideband processing mode enum 
 * \param  string  String to look up
 * \return value
 */
static ObitASDMSideBMode LookupSideBMode(gchar *name)
{
  ObitASDMSideBMode out = 0;
 
  if (!strncmp (name, "NONE", 4))                         
    return ASDMSideBMode_NONE;
  if (!strncmp (name, "PHASE_SWITCH_SEPARATION", 23))     
    return ASDMSideBMode_PHASE_SWITCH_SEPARATION;
  if (!strncmp (name, "FREQUENCY_OFFSET_SEPARATION", 27)) 
    return ASDMSideBMode_FREQUENCY_OFFSET_SEPARATION;
  if (!strncmp (name, "PHASE_SWITCH_REJECTION", 22))      
    return ASDMSideBMode_PHASE_SWITCH_REJECTION;
  if (!strncmp (name, "FREQUENCY_OFFSET_REJECTION", 26))  
    return ASDMSideBMode_FREQUENCY_OFFSET_REJECTION;
  return out;
} /* end LookupSideBMode */

/**  Look up  sideband processing mode enum 
 * \param  string  String to look up
 * \return value
 */
static ObitASDMWindowFn LookupWindowFn(gchar *name)
{
  ObitASDMWindowFn out = 0;
  
  if (!strncmp (name, "UNIFORM", 7))           return ASDMWindowFn_UNIFORM;
  if (!strncmp (name, "HANNING", 7))           return ASDMWindowFn_HANNING;
  if (!strncmp (name, "HAMMING", 7))           return ASDMWindowFn_HAMMING;
  if (!strncmp (name, "BARTLETT", 8))          return ASDMWindowFn_BARTLETT;
  if (!strncmp (name, "BLACKMANN", 9))         return ASDMWindowFn_BLACKMANN;
  if (!strncmp (name, "BLACKMANN_HARRIS", 16)) return ASDMWindowFn_BLACKMANN_HARRIS;
  if (!strncmp (name, "WELCH", 5))             return ASDMWindowFn_WELCH;
  return out;
} /* end LookupWindowFn */

/**  Look up atmospheric phase correction enum 
 * \param  string  String to look up
 * \return value
 */
static ObitASDMAtmPhCorr LookupAtmPhCorr(gchar *name)
{
  ObitASDMAtmPhCorr out = 0;
  
  if (!strncmp (name, "AP_UNCORRECTED", 14))   return ASDMAtmPhCorr_AP_UNCORRECTED;
  if (!strncmp (name, "AP_CORRECTED", 12))     return ASDMAtmPhCorr_AP_CORRECTED;
  return out;
} /* end LookupAtmPhCorr */

/** Constructor for ASDM table parsing from file
 * Reads schema version number
 * Expect the number of rows in the line following the <Name> line.
 *       <Name>Doppler</Name>
 *       <NumberRows>0</NumberRows>
 * \param  ASDMFile filename containing table
 * \param  err  ObitErr for reporting errors.
 * \return number of rows in table
 */
static ASDMTable* ParseASDMTable(gchar *ASDMFile, 
				 ObitErr *err)
{
  ASDMTable* out = NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong  maxLine = 4098;
  gchar line[4099], *schemaVersion, *next, *prior, *tstr;
  gchar *routine = " ParseASDMTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMTable));

  /* Init rows count to -1 */
  out->MainRows             = -1;
  out->AntennaRows          = -1;
  out->calAtmosphereRows    = -1;
  out->calDataRows          = -1;
  out->calDeviceRows        = -1;
  out->calPointingRows      = -1;
  out->CalReductionRows     = -1;
  out->CalWVRRows           = -1;
  out->ConfigDescriptionRows= -1;
  out->CorrelatorModeRows   = -1;
  out->DataDescriptionRows  = -1;
  out->DopplerRows          = -1;
  out->ExecBlockRows        = -1;
  out->FeedRows             = -1;
  out->FieldRows            = -1;
  out->FlagRows             = -1;
  out->PointingRows         = -1;
  out->PointingModelRows    = -1;
  out->PolarizationRows     = -1;
  out->ProcessorRows        = -1;
  out->ReceiverRows         = -1;
  out->SBSummaryRows        = -1;
  out->ScanRows             = -1;
  out->SourceRows           = -1;
  out->SpectralWindowRows   = -1;
  out->StateRows            = -1;
  out->StationRows          = -1;
  out->SubscanRows          = -1;
  out->SwitchCycleRows      = -1;
  out->SysCalRows           = -1;
  out->SysPowerRows         = -1;
  out->WeatherRows          = -1;

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, ASDMFile,OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

   /* Loop over file */
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    /* schema version ALMA and EVLA differ here */
    if (g_strstr_len (line, maxLine, "<ASDM")!=NULL) {
      out->schemaVersion = -1;
      /* ALMA */
      schemaVersion = ASDMparse_quote_str(line, maxLine, "schemaVersion=", &next);
      if (schemaVersion) {
	out->schemaVersion = strtol(schemaVersion, &next, 10);
	g_free(schemaVersion);
      }
      if (out->schemaVersion<0) { /* EVLA? */
	retCode = ObitFileReadXML (file, line, maxLine, err);
	if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
	if (retCode==OBIT_IO_EOF) break;
	schemaVersion = ASDMparse_quote_str(line, maxLine, "schemaVersion=", &next);
	if (schemaVersion) {
	  out->schemaVersion = strtol(schemaVersion, &next, 10);
	  g_free(schemaVersion);
	}
      }
      if (out->schemaVersion<0) out->schemaVersion = 2;  /* Who knows??? */
      continue;

    }

    prior = "<Name>";
    tstr = ASDMparse_str(line, maxLine, prior, &next);
    /* If not found - skip */
    if (tstr==NULL) continue;
    Strip (tstr);  /* remove blanks */

    /* Number of Main rows */
    if (!strcmp(tstr,"Main")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->MainRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Antenna rows */
    if (!strcmp(tstr,"Antenna")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->AntennaRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of calAtmosphere rows */
    if (!strcmp(tstr,"CalAtmosphere")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->calAtmosphereRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of calData rows */
    if (!strcmp(tstr,"CalData")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->calDataRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of calDevice rows */
    if (!strcmp(tstr,"CalDevice")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->calDeviceRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of calPointing rows */
    if (!strcmp(tstr,"CalPointing")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->calPointingRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of CalReduction rows */
    if (!strcmp(tstr,"CalReduction")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->CalReductionRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of CalWVR rows */
    if (!strcmp(tstr,"CalWVR")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->CalWVRRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of ConfigDescription rows */
    if (!strcmp(tstr,"ConfigDescription")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->ConfigDescriptionRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of CorrelatorMode rows */
    if (!strcmp(tstr,"CorrelatorMode")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->CorrelatorModeRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of DataDescription rows */
    if (!strcmp(tstr,"DataDescription")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->DataDescriptionRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Doppler rows */
    if (!strcmp(tstr,"Doppler")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->DopplerRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of ExecBlock rows */
    if (!strcmp(tstr,"ExecBlock")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->ExecBlockRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Feed rows */
    if (!strcmp(tstr,"Feed")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->FeedRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Field rows */
    if (!strcmp(tstr,"Field")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->FieldRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Flag rows */
    if (!strcmp(tstr,"Flag")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->FlagRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Pointing rows */
    if (!strcmp(tstr,"Pointing")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->PointingRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of PointingModel rows */
    if (!strcmp(tstr,"PointingModel")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->PointingModelRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Polarization rows */
    if (!strcmp(tstr,"Polarization")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->PolarizationRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Processor rows */
    if (!strcmp(tstr,"Processor")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->ProcessorRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Receiver rows */
    if (!strcmp(tstr,"Receiver")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->ReceiverRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of SBSummary rows */
    if (!strcmp(tstr,"SBSummary")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->SBSummaryRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Scan rows */
    if (!strcmp(tstr,"Scan")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->ScanRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Source rows */
    if (!strcmp(tstr,"Source")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->SourceRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of SpectralWindow rows */
    if (!strcmp(tstr,"SpectralWindow")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->SpectralWindowRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of State rows */
    if (!strcmp(tstr,"State")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->StateRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Station rows */
    if (!strcmp(tstr,"Station")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->StationRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Subscan rows */
    if (!strcmp(tstr,"Subscan")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->SubscanRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of SwitchCycle rows */
    if (!strcmp(tstr,"SwitchCycle")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->SwitchCycleRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of SysCal rows */
    if (!strcmp(tstr,"SysCal")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->SysCalRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of SysPower rows */
    if (!strcmp(tstr,"SysPower")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->SysPowerRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    /* Number of Weather rows */
    if (!strcmp(tstr,"Weather")) {
      retCode = ObitFileReadXML (file, line, maxLine, err);
      if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
      if (retCode==OBIT_IO_EOF) break;
      out->WeatherRows = ASDMparse_int(line, maxLine, "<NumberRows>", &next);
      continue;
    }
    if (tstr); g_free(tstr); tstr = NULL;
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  if (out->schemaVersion<0) out->schemaVersion = 2;  /* Who knows??? */

  return out;
} /* end ParseASDMTable */

/** 
 * Destructor for ASDM table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMTable* KillASDMTable(ASDMTable* table)
{
  if (table==NULL) return NULL;
  g_free(table);
  return NULL;
} /* end KillASDMTable */

/* -----------------------------------------------------------------------------*/
/* ----------------------------  Main -----------------------------------------*/
/** 
 * Destructor for Main table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMMainRow* KillASDMMainRow(ASDMMainRow* row)
{
  if (row == NULL) return NULL;
  if (row->stateId)  g_free(row->stateId);
  if (row->entityId) g_free(row->entityId);
  g_free(row);
  return NULL;
} /* end   KillASDMMainRow */

/** 
 * Constructor for Main table parsing from file
 * \param  MainFile Name of file containing table
 * \param  err      ObitErr for reporting errors.
 * \return table structure,  use KillASDMMainTable to free
 */
static ASDMMainTable* ParseASDMMainTable(ObitSDMData *me, 
					 gchar *MainFile, 
					 ObitErr *err)
{
  ASDMMainTable* out=NULL;
  ObitIOCode retCode;
  olong i, irow, colon, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next, *tstr;
  ObitFile *file=NULL;
  gchar *routine = " ParseASDMMainTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMMainTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->MainRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMMainRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMMainRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, MainFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<time>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->time = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numAntenna>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numAntenna = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<timeSampling>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      tstr =  ASDMparse_str (line, maxLine, prior, &next);
      Strip(tstr);
      if (!strcmp(tstr, "INTEGRATION")) out->rows[irow]->timeSampling = ASDMINTEGRATION;
      else if (!strcmp(tstr, "SUBINTEGRATION")) out->rows[irow]->timeSampling = ASDMSUBINTEGRATION;
      g_free(tstr);
      continue;
    }
    prior = "<interval>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->interval = ASDMparse_timeint (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numIntegration>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numIntegration = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<scanNumber>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->scanNumber = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<subscanNumber>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->subscanNumber = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<dataSize>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      if (me->schemaVersion<=2) 
	out->rows[irow]->dataSize = ASDMparse_int (line, maxLine, prior, &next);
      else if (me->schemaVersion==3) 
	out->rows[irow]->dataSize = ASDMparse_int64 (line, maxLine, prior, &next);
      else g_error("Unsupported schemaVersion %d", me->schemaVersion);
      continue;
    }
    prior = "<dataOid>"; 
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      /* Parse off next line
	 retCode = ObitFileReadXML (file, line, maxLine, err);
	 if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
	 if (retCode==OBIT_IO_EOF) break;*/
      prior = "entityId=";
      tstr = ASDMparse_quote_str (line, maxLine, prior, &next);
      /* Replace slashes with underscore - drop up to colon */
      colon = -1;
      for (i=0; i<strlen(tstr); i++) {	
	if (tstr[i]=='/') tstr[i]='_';
	if (tstr[i]==':') colon = i;
      }
      out->rows[irow]->entityId = g_strdup(&tstr[colon+1]);
      g_free(tstr);
      continue;
    }
    /* ALMA version 3 */
    prior = "<dataUID>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      /* Parse off next line
	 retCode = ObitFileReadXML (file, line, maxLine, err);
	 if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
	 if (retCode==OBIT_IO_EOF) break; */
      prior = "entityId=";
      tstr = ASDMparse_quote_str (line, maxLine, prior, &next);
      /* Replace slashes with underscore - drop up to colon */
      colon = -1;
      for (i=0; i<strlen(tstr); i++) {	
	if (tstr[i]=='/') tstr[i]='_';
	if (tstr[i]==':') colon = i;
      }
      out->rows[irow]->entityId = g_strdup(&tstr[colon+1]);
      g_free(tstr);
      continue;
    }
    prior = "<flagRow>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->flagRow = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }
    prior = "<configDescriptionId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "ConfigDescription_";
      out->rows[irow]->configDescriptionId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<execBlockId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "ExecBlock_";
      out->rows[irow]->execBlockId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<fieldId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Field_";
      out->rows[irow]->fieldId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<stateId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->stateId = ASDMparse_enumarray (line, maxLine, prior, &next);
      continue;
    }

   /* Is this the end of a row? */
   if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
   Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMMainTable */

/** 
 * Destructor for Main table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMMainTable* KillASDMMainTable(ASDMMainTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMMainRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMMainTable */

/* ----------------------  Antenna ----------------------------------- */
/** 
 * Destructor for Antenna table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMAntennaRow* KillASDMAntennaRow(ASDMAntennaRow* row)
{
  if (row == NULL) return NULL;
  if (row->name)     g_free(row->name);
  if (row->position) g_free(row->position);
  if (row->offset)   g_free(row->offset);
  g_free(row);
  return NULL;
} /* end   KillASDMAntennaRow */

/** 
 * Constructor for Antenna table parsing from file
 * \param  AntennaFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMAntennaTable to free
 */
static ASDMAntennaTable* 
ParseASDMAntennaTable(ObitSDMData *me, 
		      gchar *AntennaFile, 
		      ObitErr *err)
{
  ASDMAntennaTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next, *tstr;
  gchar *routine = " ParseASDMAntennaTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMAntennaTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->AntennaRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMAntennaRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMAntennaRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, AntennaFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    /* Antenna ID */
    prior = "<antennaId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Antenna_";
      out->rows[irow]->antennaId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    /* name */
    prior = "<name>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->name = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->name);   /* Deblank */
      continue;
    }
    /* antenna Make enum */
    prior = "<antennaMake>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      tstr =  ASDMparse_str (line, maxLine, prior, &next);
      Strip(tstr);
      out->rows[irow]->antennaMake = LookupAntennaMake(tstr);
      g_free(tstr);
    }
    /* antenna Type enum */
    prior = "<antennaType>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      tstr =  ASDMparse_str (line, maxLine, prior, &next);
      Strip(tstr);
      if (!strcmp(tstr, "GROUND_BASED"))      out->rows[irow]->antennaType = ASDMAnt_GROUND_BASED;
      else if (!strcmp(tstr, "SPACE_BASED"))  out->rows[irow]->antennaType = ASDMAnt_SPACE_BASED;
      else if (!strcmp(tstr, "TRACKING_STN")) out->rows[irow]->antennaType = ASDMAnt_TRACKING_STN;
      g_free(tstr);
    }
    /* dish Diameter (m)*/
    prior = "<dishDiameter>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->dishDiameter = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    /* position array of doubles */
    prior = "<position>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->position = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    /* offset array of doubles */
    prior = "<offset>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->offset = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    /* time */
    prior = "<time>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->time = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    /* station Id */
    prior = "<stationId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Station_";
      out->rows[irow]->stationId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMAntennaTable */

/** 
 * Destructor for Antenna table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMAntennaTable* KillASDMAntennaTable(ASDMAntennaTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMAntennaRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMAntennaTable */

/* ----------------------  calAtmosphere ----------------------------------- */
/** 
 * Destructor for calAtmosphere table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMcalAtmosphereRow* KillASDMcalAtmosphereRow(ASDMcalAtmosphereRow* row)
{
  if (row == NULL) return NULL;
  if (row->antennaName)        g_free(row->antennaName);
  if (row->syscalType)         g_free(row->syscalType);
  if (row->polarizationTypes)  g_free(row->polarizationTypes);
  if (row->tAtm)               g_free(row->tAtm);
  if (row->tRec)               g_free(row->tRec);
  if (row->tSys)               g_free(row->tSys);
  if (row->tau)                g_free(row->tau);
  if (row->water)              g_free(row->water);
  if (row->waterError)         g_free(row->waterError);
  if (row->forwardEfficiency)  g_free(row->forwardEfficiency);
  if (row->sbGain)             g_free(row->sbGain);
  g_free(row);
  return NULL;
} /* end   KillASDMcalAtmosphereRow */

/** 
 * Constructor for calAtmosphere table parsing from file
 * \param  calAtmosphereFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMcalAtmosphereTable to free
 */
static ASDMcalAtmosphereTable* 
ParseASDMcalAtmosphereTable(ObitSDMData *me, 
		      gchar *calAtmosphereFile, 
		      ObitErr *err)
{
  ASDMcalAtmosphereTable* out=NULL;
  ObitIOCode retCode;
  olong i, j, charLeft, irow, maxLine = me->maxLine, ndim, naxis;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next, *b;
  ObitFile *file=NULL;
  gchar *routine = " ParseASDMcalAtmosphereTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMcalAtmosphereTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->calAtmosphereRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMcalAtmosphereRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMcalAtmosphereRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, calAtmosphereFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    /* receiverBand */
    prior = "<receiverBand>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "ALMA_RB_";
      out->rows[irow]->receiverBand = ASDMparse_int (line, maxLine, prior, &next);
    }
    /* antenna Name */
    prior = "<antennaName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->antennaName = ASDMparse_str (line, maxLine, prior, &next);
    }
    /* syscalType */
    prior = "<syscalType>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->syscalType = ASDMparse_str (line, maxLine, prior, &next);
    }
    /* baseband Id */
    prior = "<basebandName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "BB_";
      out->rows[irow]->basebandId = ASDMparse_int (line, maxLine, prior, &next);
    }
    /* number of frequencies */
    prior = "<numFreq>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numFreq = ASDMparse_int (line, maxLine, prior, &next);
    }
    /* number of Loads */
    prior = "<numLoad>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numLoad = ASDMparse_int (line, maxLine, prior, &next);
    }
    /* number of Receptors */
    prior = "<numReceptor>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numReceptor = ASDMparse_int (line, maxLine, prior, &next);
    }
    /* cal Id */
    prior = "<calDataId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "<calDataId>CalData_";
      out->rows[irow]->calDataId = ASDMparse_int (line, maxLine, prior, &next);
    }
    /*  cal reduction Id */
    prior = "<calReductionId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "<calReductionId>CalReduction_";
      out->rows[irow]->calReductionId = ASDMparse_int (line, maxLine, prior, &next);
    }
    /* polarization types, e.g. 'X', 'Y', convert to enums */
    prior = "<polarizationTypes>";
    b = g_strstr_len (line, maxLine, prior);
    if (b!=NULL) {
      b += strlen(prior);
      /* Parse array of strings */
      /* Get dimensionality - only can handle 1 */
      ndim = (olong)strtol(b, &next, 10);
      g_assert(ndim==1);  /* bother */
      b = next;
      /* get number of values */
      naxis = (olong)strtol(b, &next, 10);
      out->rows[irow]->polarizationTypes = g_malloc0(naxis*sizeof(ObitASDMPolnType));
      b = next;
      /* Find next non blank */
      while (*b==' ') b++;
      /* Cycle through values */
      for (i=0; i<naxis; i++) {
	/* Find next blank or '<' */
	charLeft =  maxLine - (b-line);
	for (j=0; j<charLeft; j++) {
	  if ((b[j]==' ') || (b[j]=='<')) {
	    b[j] = 0; break;
	  }
	}
	out->rows[irow]->polarizationTypes[i] =  LookupPolnType(b);
	b += j+1; /* Next value */
      }
      continue;
    }
    /* start time (days) */
    prior = "<startValidTime>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->startValidTime = ASDMparse_time (line, maxLine, prior, &next);
    }
    /* end time (days) */
    prior = "<endValidTime>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->endValidTime = ASDMparse_time (line, maxLine, prior, &next);
    }
    /* Ground pressure */
    prior = "<groundPressure>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->groundPressure = ASDMparse_dbl (line, maxLine, prior, &next);
    }
    /* Ground relative humidity */
    prior = "<groundRelHumidity>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->groundRelHumidity = ASDMparse_dbl (line, maxLine, prior, &next);
    }
    /* Ground temperature (K) */
    prior = "<groundTemperature>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->groundTemperature = ASDMparse_int (line, maxLine, prior, &next);
    }
    /* Frequency range (Hz) */
    prior = "<frequencyRange>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->frequencyRange = ASDMparse_dblarray (line, maxLine, prior, &next);
    }
    /* Atm. temp per receptor */
    prior = "<tAtm>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->tAtm = ASDMparse_dblarray (line, maxLine, prior, &next);
    }
    /* Receiver temp. per receptor */
    prior = "<tRec>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->tRec = ASDMparse_dblarray (line, maxLine, prior, &next);
    }
    /* System  temp. per receptor */
    prior = "<tSys>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->tSys = ASDMparse_dblarray (line, maxLine, prior, &next);
    }
    /* Opacity per receptor */
    prior = "<tau>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->tau = ASDMparse_dblarray (line, maxLine, prior, &next);
    }
    /* Water per receptor */
    prior = "<water>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->water = ASDMparse_dblarray (line, maxLine, prior, &next);
    }
    /* Error in water, per receptor */
    prior = "<waterError>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->waterError = ASDMparse_dblarray (line, maxLine, prior, &next);
    }
    /* forward efficiency, per receptor */
    prior = "<forwardEfficiency>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->forwardEfficiency = ASDMparse_dblarray (line, maxLine, prior, &next);
    }
    /* sb(?) gain, per receptor  */
    prior = "<sbGain>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->sbGain = ASDMparse_dblarray (line, maxLine, prior, &next);
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMcalAtmosphereTable */

/** 
 * Destructor for calAtmosphere table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMcalAtmosphereTable* KillASDMcalAtmosphereTable(ASDMcalAtmosphereTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMcalAtmosphereRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMcalAtmosphereTable */

/* ----------------------  calData ----------------------------------- */
/** 
 * Destructor for calData table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMcalDataRow* KillASDMcalDataRow(ASDMcalDataRow* row)
{
  if (row == NULL) return NULL;
  g_free(row);
  return NULL;
} /* end   KillASDMcalDataRow */

/** 
 * Constructor for calData table parsing from file
 * \param  calDataFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMcalDataTable to free
 */
static ASDMcalDataTable* 
ParseASDMcalDataTable(ObitSDMData *me, 
		      gchar *calDataFile, 
		      ObitErr *err)
{
  ASDMcalDataTable* out=NULL;
  ObitIOCode retCode;
  olong irow, maxLine =  me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  /*gchar *prior, *next;*/
  ObitFile *file=NULL;
  gchar *routine = " ParseASDMcalDataTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMcalDataTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->calDataRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMcalDataRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMcalDataRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, calDataFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMcalDataTable */

/** 
 * Destructor for calData table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMcalDataTable* KillASDMcalDataTable(ASDMcalDataTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMcalDataRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMcalDataTable */

/* ----------------------  calDevice ----------------------------------- */
/** 
 * Destructor for calDevice table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMcalDeviceRow* KillASDMcalDeviceRow(ASDMcalDeviceRow* row)
{
  olong i, n;
  if (row == NULL) return NULL;
  if (row->timeInterval)    g_free(row->timeInterval);
  if (row->calEff)          g_free(row->calEff);
  if (row->noiseCal)        g_free(row->noiseCal);
  if (row->coupledNoiseCal) g_free(row->coupledNoiseCal);
  if (row->temperatureLoad) g_free(row->temperatureLoad);
  if (row->calLoadNames) {
    n = row->numCalLoad;
    for (i=0; i<n; i++) if (row->calLoadNames[i]) g_free(row->calLoadNames[i]);
    g_free(row->calLoadNames);
  }
  g_free(row);
  return NULL;
} /* end   KillASDMcalDeviceRow */

/** 
 * Constructor for calDevice table parsing from file
 * \param  calDeviceFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMcalDeviceTable to free
 */
static ASDMcalDeviceTable* 
ParseASDMcalDeviceTable(ObitSDMData *me, 
		      gchar *calDeviceFile, 
		      ObitErr *err)
{
  ASDMcalDeviceTable* out=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  ObitFile *file=NULL;
  odouble mjdJD0=2400000.5; /* JD of beginning of MJD time */
  gchar *routine = " ParseASDMcalDeviceTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMcalDeviceTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->calDeviceRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMcalDeviceRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMcalDeviceRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, calDeviceFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<antennaId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Antenna_";
      out->rows[irow]->antennaId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<spectralWindowId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "SpectralWindow_";
      out->rows[irow]->spectralWindowId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<feedId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->feedId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numReceptor>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numReceptor = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<timeInterval>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->timeInterval = ASDMparse_timeRange(line, maxLine, prior, &next);
      /* Remove offset from second */
      if ((out->rows[irow]->timeInterval[1]<out->rows[irow]->timeInterval[0]) &&
	  (out->rows[irow]->timeInterval[1]>mjdJD0))
	out->rows[irow]->timeInterval[1] -= mjdJD0;
      continue;
    }
    prior = "<numCalload>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numCalLoad = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<calLoadNames>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->calLoadNames = ASDMparse_strarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<calEff>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->calEff = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<noiseCal>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->noiseCal = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<coupledNoiseCal>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->coupledNoiseCal = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<temperatureLoad>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->temperatureLoad = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMcalDeviceTable */

/** 
 * Destructor for calDevice table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMcalDeviceTable* KillASDMcalDeviceTable(ASDMcalDeviceTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMcalDeviceRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMcalDeviceTable */

/* ----------------------  calPointing ----------------------------------- */
/** 
 * Destructor for calPointing table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMcalPointingRow* KillASDMcalPointingRow(ASDMcalPointingRow* row)
{
  if (row == NULL) return NULL;
  g_free(row);
  return NULL;
} /* end   KillASDMcalPointingRow */

/** 
 * Constructor for calPointing table parsing from file
 * \param  calPointingFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMcalPointingTable to free
 */
static ASDMcalPointingTable* 
ParseASDMcalPointingTable(ObitSDMData *me, 
			  gchar *calPointingFile, 
			  ObitErr *err)
{
  ASDMcalPointingTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  /*gchar *prior, *next;*/
  gchar *routine = " ParseASDMcalPointingTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMcalPointingTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->calPointingRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMcalPointingRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMcalPointingRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, calPointingFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMcalPointingTable */

/** 
 * Destructor for calPointing table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMcalPointingTable* KillASDMcalPointingTable(ASDMcalPointingTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMcalPointingRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMcalPointingTable */

/* ----------------------  CalReduction ----------------------------------- */
/** 
 * Destructor for CalReduction table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMCalReductionRow* KillASDMCalReductionRow(ASDMCalReductionRow* row)
{
  if (row == NULL) return NULL;
  g_free(row);
  return NULL;
} /* end   KillASDMCalReductionRow */

/** 
 * Constructor for CalReduction table parsing from file
 * \param  CalReductionFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMCalReductionTable to free
 */
static ASDMCalReductionTable* 
ParseASDMCalReductionTable(ObitSDMData *me, 
			   gchar *CalReductionFile, 
			   ObitErr *err)
{
  ASDMCalReductionTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  /*gchar *prior, *next;*/
  gchar *routine = " ParseASDMCalReductionTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMCalReductionTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->CalReductionRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMCalReductionRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMCalReductionRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, CalReductionFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMCalReductionTable */

/** 
 * Destructor for CalReduction table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMCalReductionTable* KillASDMCalReductionTable(ASDMCalReductionTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMCalReductionRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMCalReductionTable */

/* ----------------------  CalWVR ----------------------------------- */
/** 
 * Destructor for CalWVR table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMCalWVRRow* KillASDMCalWVRRow(ASDMCalWVRRow* row)
{
  if (row == NULL) return NULL;
  if (row->wvrMethod)      g_free(row->wvrMethod);
  if (row->antennaName)    g_free(row->antennaName);
  if (row->chanFreq)       g_free(row->chanFreq);
  if (row->chanWidth)      g_free(row->chanWidth);
  if (row->refTemp)        g_free(row->refTemp);
  if (row->pathCoeff)      g_free(row->pathCoeff);
  if (row->polyFreqLimits) g_free(row->polyFreqLimits);
  g_free(row);
  return NULL;
} /* end   KillASDMCalWVRRow */

/** 
 * Constructor for CalWVR table parsing from file
 * \param  CalWVRFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMCalWVRTable to free
 */
static ASDMCalWVRTable* 
ParseASDMCalWVRTable(ObitSDMData *me, 
			   gchar *CalWVRFile, 
			   ObitErr *err)
{
  ASDMCalWVRTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  odouble *darry=NULL;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMCalWVRTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMCalWVRTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->CalWVRRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMCalWVRRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMCalWVRRow));

  file = newObitFile("CalWVR");
  retCode = ObitFileOpen(file, CalWVRFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<startValidTime>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->startValidTime = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    prior = "<endValidTime>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->endValidTime = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    /* method */
    prior = "<wvrMethod>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->wvrMethod = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->wvrMethod);   /* Deblank */
      continue;
    }
    /* antenna Name */
    prior = "<antennaName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->antennaName = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->antennaName);   /* Deblank */
      continue;
    }
     prior = "<numChan>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numChan = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
     prior = "<chanFreq>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->chanFreq = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<chanWidth>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->chanWidth = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<refTemp>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->refTemp = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
     prior = "<numPoly>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numPoly = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<pathCoeff>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->pathCoeff = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<polyFreqLimits>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->polyFreqLimits = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<wetPath>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      darry = ASDMparse_dblarray (line, maxLine, prior, &next);
      if (darry) {
	out->rows[irow]->wetPath = darry[0];
	g_free(darry);
      }
      continue;
    }
    prior = "<dryPath>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      darry = ASDMparse_dblarray (line, maxLine, prior, &next);
      if (darry) {
	out->rows[irow]->dryPath = darry[0];
	g_free(darry);
      }
      continue;
    }
    prior = "<water>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->water = ASDMparse_dbl (line, maxLine, prior, &next);
      /*fprintf (stdout, "wvr %s %lf %lf %lf \n",
	out->rows[irow]->antennaName, out->rows[irow]->wetPath, out->rows[irow]->dryPath, 
	out->rows[irow]->water);*/
      continue;
   }
    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMCalWVRTable */

/** 
 * Destructor for CalWVR table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMCalWVRTable* KillASDMCalWVRTable(ASDMCalWVRTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMCalWVRRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMCalWVRTable */

/* ---------------------- ConfigDescription ----------------------------------- */
/** 
 * Destructor for ConfigDescription table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMConfigDescriptionRow* KillASDMConfigDescriptionRow(ASDMConfigDescriptionRow* row)
{
  if (row == NULL) return NULL;
  if (row->atmPhaseCorrection) g_free(row->atmPhaseCorrection);
  if (row->antennaId)          g_free(row->antennaId);
  if (row->dataDescriptionId)  g_free(row->dataDescriptionId);
  if (row->feedId)             g_free(row->feedId);
  if (row->switchCycleId)      g_free(row->switchCycleId);
  g_free(row);
  return NULL;
} /* end   KillASDMConfigDescriptionRow */

/** 
 * Constructor for ConfigDescription table parsing from file
 * \param  ConfigDescriptionFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMConfigDescriptionTable to free
 */
static ASDMConfigDescriptionTable* 
ParseASDMConfigDescriptionTable(ObitSDMData *me, 
				gchar *ConfigDescriptionFile, 
				ObitErr *err)
{
  ASDMConfigDescriptionTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong i, j, irow, ndim, naxis, maxLine = me->maxLine, charLeft;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next, *tstr, *b;
  gchar *routine = " ParseASDMConfigDescriptionTable";
  
  /* error checks */
  if (err->error) return out;
  
  out = g_malloc0(sizeof(ASDMConfigDescriptionTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->ConfigDescriptionRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMConfigDescriptionRow*));
  for (irow=0; irow<out->nrows; irow++) 
    out->rows[irow] = g_malloc0(sizeof(ASDMConfigDescriptionRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, ConfigDescriptionFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  
  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<numAntenna>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numAntenna = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numDataDescription>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numDataDescription = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numFeed>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numFeed = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<correlationMode>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      tstr =  ASDMparse_str (line, maxLine, prior, &next);
      Strip(tstr);
      if (!strcmp(tstr, "CROSS_ONLY"))          
	out->rows[irow]->correlationMode = ASDMCorrMode_CROSS_ONLY;
      else if (!strcmp(tstr, "AUTO_ONLY"))      
	out->rows[irow]->correlationMode = ASDMCorrMode_AUTO_ONLY;
      else if (!strcmp(tstr, "CROSS_AND_AUTO")) 
	out->rows[irow]->correlationMode = ASDMCorrMode_CROSS_AND_AUTO;
      g_free(tstr);
    }
    prior = "<configDescriptionId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "ConfigDescription_";
      out->rows[irow]->configDescriptionId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numAtmPhaseCorrection>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numAtmPhaseCorrection = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<atmPhaseCorrection>";
    b = g_strstr_len (line, maxLine, prior);
    if (b!=NULL) {
      b += strlen(prior);
      /* Get dimensionality - only can handle 1 */
      ndim = (olong)strtol(b, &next, 10);
      g_assert(ndim==1);  /* bother */
      b = next;
      /* get number of values */
      naxis= (olong)strtol(b, &next, 10);
      b = next + 1;  /* Start of first value */
      /* Output */
      out->rows[irow]->atmPhaseCorrection = g_malloc0(naxis*sizeof(ObitASDMAtmPhCorr));
      /* Cycle through values */
      for (i=0; i<naxis; i++) {
	/* Find next blank or '<' */
	charLeft =  maxLine - (b-line);
	for (j=0; j<charLeft; j++) {
	  if ((b[j]==' ') || (b[j]=='<')) {
	    b[j] = 0; break;
	  }
	}
	if (!strcmp(b, "AP_UNCORRECTED"))        
	  out->rows[irow]->atmPhaseCorrection[i] = ASDMAtmPhCorr_AP_UNCORRECTED;
	else if (!strcmp(b, "AP_CORRECTED"))   
	  out->rows[irow]->atmPhaseCorrection[i] = ASDMAtmPhCorr_AP_CORRECTED;
	b += j+1; /* Next value */
      }
      continue;
    }
    prior = "<processorType>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      tstr =  ASDMparse_str (line, maxLine, prior, &next);
      Strip(tstr);
      if (!strcmp(tstr, "CORRELATOR"))        
	out->rows[irow]->processorType = ASDMProcrType_CORRELATOR;
      else if (!strcmp(tstr, "RADIOMETER"))   
	out->rows[irow]->processorType = ASDMProcrType_RADIOMETER;
      g_free(tstr);
    }
    prior = "<spectralType>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      tstr =  ASDMparse_str (line, maxLine, prior, &next);
      Strip(tstr);
      if (!strcmp(tstr, "CHANNEL_AVERAGE"))      
	out->rows[irow]->spectralType = ASDMSpecRes_CHANNEL_AVERAGE;
      else if (!strcmp(tstr, "BASEBAND_WIDE"))   
	out->rows[irow]->spectralType = ASDMSpecRes_BASEBAND_WIDE;
      else if (!strcmp(tstr, "FULL_RESOLUTION")) 
	out->rows[irow]->spectralType = ASDMSpecRes_FULL_RESOLUTION;
      g_free(tstr);
    }
    prior = "<antennaId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->antennaId = ASDMparse_enumarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<dataDescriptionId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->dataDescriptionId = ASDMparse_enumarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<feedId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->feedId = ASDMparse_intarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<processorId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->processorId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<switchCycleId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->switchCycleId = ASDMparse_enumarray (line, maxLine, prior, &next);
      continue;
    }
    
    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;
    
    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */
  
  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);
  
  return out;
} /* end ParseASDMConfigDescriptionTable */

/** 
 * Destructor for ConfigDescription table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMConfigDescriptionTable* 
KillASDMConfigDescriptionTable(ASDMConfigDescriptionTable* table)
{
  olong i;
  
  if (table==NULL) return NULL;  /* Anybody home? */
  
  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMConfigDescriptionRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMConfigDescriptionTable */

/* ----------------------  CorrelatorMode ----------------------------------- */
/** 
 * Destructor for CorrelatorMode table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMCorrelatorModeRow* KillASDMCorrelatorModeRow(ASDMCorrelatorModeRow* row)
{
  if (row == NULL) return NULL;
  if (row->basebandNames)  g_free(row->basebandNames);
  if (row->basebandConfig) g_free(row->basebandConfig);
  if (row->axesOrderArray) g_free(row->axesOrderArray);
  if (row->filterMode)     g_free(row->filterMode);
  if (row->correlatorName) g_free(row->correlatorName);
  g_free(row);
  return NULL;
} /* end   KillASDMCorrelatorModeRow */

/** 
 * Constructor for CorrelatorMode table parsing from file
 * \param  CorrelatorModeFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMCorrelatorModeTable to free
 */
static ASDMCorrelatorModeTable* 
ParseASDMCorrelatorModeTable(ObitSDMData *me, 
			     gchar *CorrelatorModeFile, 
			     ObitErr *err)
{
  ASDMCorrelatorModeTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong i, j, irow, charLeft, ndim, naxis, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next, *b, *tstr;
  gchar *routine = " ParseASDMCorrelatorModeTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMCorrelatorModeTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->CorrelatorModeRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMCorrelatorModeRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMCorrelatorModeRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, CorrelatorModeFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<correlatorModeId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "CorrelatorMode_";
      out->rows[irow]->correlatorModeId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numBaseband>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numBaseband = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

    prior = "<basebandNames>";
    b = g_strstr_len (line, maxLine, prior);
    if (b!=NULL) {
      b += strlen(prior);
      /* Parse array of strings */
      /* Get dimensionality - only can handle 1 */
      ndim = (olong)strtol(b, &next, 10);
      g_assert(ndim==1);  /* bother */
      b = next;
      /* get number of values */
      naxis = (olong)strtol(b, &next, 10);
      out->rows[irow]->basebandNames = g_malloc0(naxis*sizeof(ObitASDMBasebandName));
     /* Cycle through values */
      for (i=0; i<naxis; i++) {
	/* Find next blank or '<' */
	charLeft =  maxLine - (b-line);
	for (j=0; j<charLeft; j++) {
	  if ((b[j]==' ') || (b[j]=='<')) {
	    b[j] = 0; break;
	  }
	}
	out->rows[irow]->basebandNames[i] =  LookupBasebandName(b);
	b += j+1; /* Next value */
      }
      continue;
    }

    prior = "<basebandConfig>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->basebandConfig = ASDMparse_intarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<accumMode>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      tstr =  ASDMparse_str (line, maxLine, prior, &next);
      Strip(tstr);
      if (!strcmp(tstr, "FAST"))      
	 out->rows[irow]->accumMode = ASDMAccumMode_FAST;
      else if (!strcmp(tstr, "NORMAL"))   
	out->rows[irow]->accumMode = ASDMAccumMode_NORMAL;
      else if (!strcmp(tstr, "UNDEFINED")) 
	out->rows[irow]->accumMode = ASDMAccumMode_UNDEFINED;
      g_free(tstr);
      continue;
    }
    prior = "<binMode>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->binMode = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numAxes>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numAxes = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

    prior = "<axesOrderArray>";
    b = g_strstr_len (line, maxLine, prior);
    if (b!=NULL) {
      b += strlen(prior);
      /* Parse array of strings */
      /* Get dimensionality - only can handle 1 */
      ndim = (olong)strtol(b, &next, 10);
      g_assert(ndim==1);  /* bother */
      b = next;
      /* get number of values */
      naxis = (olong)strtol(b, &next, 10);
      out->rows[irow]->axesOrderArray = g_malloc0(naxis*sizeof(ObitASDMAxisName));
     /* Cycle through values */
      for (i=0; i<naxis; i++) {
	/* Find next blank or '<' */
	charLeft =  maxLine - (b-line);
	for (j=0; j<charLeft; j++) {
	  if ((b[j]==' ') || (b[j]=='<')) {
	    b[j] = 0; break;
	  }
	}
	out->rows[irow]->axesOrderArray[i] =  LookupAxisName(b);
	b += j+1; /* Next value */
      }
      continue;
    }

    prior = "<filterMode>";
    b = g_strstr_len (line, maxLine, prior);
    if (b!=NULL) {
      b += strlen(prior);
      /* Parse array of strings */
      /* Get dimensionality - only can handle 1 */
      ndim = (olong)strtol(b, &next, 10);
      g_assert(ndim==1);  /* bother */
      b = next;
      /* get number of values */
      naxis = (olong)strtol(b, &next, 10);
      out->rows[irow]->filterMode = g_malloc0(naxis*sizeof(ObitASDMFilterMode));
     /* Cycle through values */
      for (i=0; i<naxis; i++) {
	/* Find next blank or '<' */
	charLeft =  maxLine - (b-line);
	for (j=0; j<charLeft; j++) {
	  if ((b[j]==' ') || (b[j]=='<')) {
	    b[j] = 0; break;
	  }
	}
	if (!strcmp(b, "FILTER_NA"))      
	  out->rows[irow]->filterMode[i] = ASDMFilterMode_FILTER_NA;
	else if (!strcmp(b, "FILTER_TDM"))   
	  out->rows[irow]->filterMode[i] = ASDMFilterMode_FILTER_TDM;
	else if (!strcmp(b, "FILTER_TFB")) 
	  out->rows[irow]->filterMode[i] = ASDMFilterMode_FILTER_TFB;
	else if (!strcmp(b, "UNDEFINED")) 
	  out->rows[irow]->filterMode[i] = ASDMFilterMode_UNDEFINED;
	b += j+1; /* Next value */
      }
      continue;
    }

    prior = "<correlatorName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->correlatorName = ASDMparse_str (line, maxLine, prior, &next);
      Strip( out->rows[irow]->correlatorName);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;
    
    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMCorrelatorModeTable */

/** 
 * Destructor for CorrelatorMode table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMCorrelatorModeTable* 
KillASDMCorrelatorModeTable(ASDMCorrelatorModeTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMCorrelatorModeRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMCorrelatorModeTable */

/* ----------------------  DataDescription ----------------------------------- */
/** 
 * Destructor for DataDescription table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMDataDescriptionRow* 
KillASDMDataDescriptionRow(ASDMDataDescriptionRow* row)
{
  if (row == NULL) return NULL;
  g_free(row);
  return NULL;
} /* end   KillASDMDataDescriptionRow */

/** 
 * Constructor for DataDescription table parsing from file
 * \param  DataDescriptionFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMDataDescriptionTable to free
 */
static ASDMDataDescriptionTable* 
ParseASDMDataDescriptionTable(ObitSDMData *me, 
			      gchar *DataDescriptionFile, 
			      ObitErr *err)
{
  ASDMDataDescriptionTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMDataDescriptionTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMDataDescriptionTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->DataDescriptionRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMDataDescriptionRow*));
  for (irow=0; irow<out->nrows; irow++) 
    out->rows[irow] = g_malloc0(sizeof(ASDMDataDescriptionRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, DataDescriptionFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<dataDescriptionId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "DataDescription_";
      out->rows[irow]->dataDescriptionId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<polOrHoloId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Polarization_";
      out->rows[irow]->polOrHoloId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

    prior = "<spectralWindowId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "SpectralWindow_";
      out->rows[irow]->spectralWindowId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMDataDescriptionTable */

/** 
 * Destructor for DataDescription table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMDataDescriptionTable* 
KillASDMDataDescriptionTable(ASDMDataDescriptionTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMDataDescriptionRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMDataDescriptionTable */

/* ----------------------  Doppler ----------------------------------- */
/** 
 * Destructor for Doppler table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMDopplerRow* KillASDMDopplerRow(ASDMDopplerRow* row)
{
  if (row == NULL) return NULL;
  g_free(row);
  return NULL;
} /* end   KillASDMDopplerRow */

/** 
 * Constructor for Doppler table parsing from file
 * \param  DopplerFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMDopplerTable to free
 */
static ASDMDopplerTable* 
ParseASDMDopplerTable(ObitSDMData *me, 
		      gchar *DopplerFile, 
		      ObitErr *err)
{
  ASDMDopplerTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  /*gchar *prior, *next;*/
  gchar *routine = " ParseASDMDopplerTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMDopplerTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->DopplerRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMDopplerRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMDopplerRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, DopplerFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMDopplerTable */

/** 
 * Destructor for Doppler table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMDopplerTable* KillASDMDopplerTable(ASDMDopplerTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMDopplerRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMDopplerTable */

/* ----------------------  ExecBlock ----------------------------------- */
/** 
 * Destructor for ExecBlock table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMExecBlockRow* KillASDMExecBlockRow(ASDMExecBlockRow* row)
{
  olong i;
  if (row == NULL) return NULL;
  if (row->configName)       g_free(row->configName);
  if (row->telescopeName)    g_free(row->telescopeName);
  if (row->observerName)     g_free(row->observerName);
  if (row->observingLog)     {
    for (i=0; i<row->numObservingLog; i++ ) {
      if (row->observingLog[i]) g_free(row->observingLog[i]);
    }
    g_free(row->observingLog);
  }
  if (row->sessionReference) g_free(row->sessionReference);
  if (row->schedulerMode)    g_free(row->schedulerMode);
  if (row->antennaId)        g_free(row->antennaId);
  g_free(row);
  return NULL;
} /* end   KillASDMExecBlockRow */

/** 
 * Constructor for ExecBlock table parsing from file
 * Determines which array
 * \param  ExecBlockFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMExecBlockTable to free
 */
static ASDMExecBlockTable* 
ParseASDMExecBlockTable(ObitSDMData *me, 
			gchar *ExecBlockFile, 
			ObitErr *err)
{
  ASDMExecBlockTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMExecBlockTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMExecBlockTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->ExecBlockRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMExecBlockRow*));
  for (irow=0; irow<out->nrows; irow++) 
    out->rows[irow] = g_malloc0(sizeof(ASDMExecBlockRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, ExecBlockFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<execBlockId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "ExecBlock_";
      out->rows[irow]->execBlockId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<startTime>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->startTime = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    prior = "<endTime>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->endTime = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    prior = "<execBlockNum>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->execBlockNum = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<configName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->configName = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->configName);
      continue;
    }
    prior = "<telescopeName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->telescopeName = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->telescopeName);
      continue;
    }

    prior = "<observerName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->observerName = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->observerName);
      continue;
    }
    prior = "<observingLog>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      if (me->schemaVersion<=2) {
	out->rows[irow]->numObservingLog = 1;
	out->rows[irow]->observingLog    = g_malloc0(sizeof(gchar*));
	out->rows[irow]->observingLog[0] = ASDMparse_str (line, maxLine, prior, &next);
	Strip(out->rows[irow]->observingLog[0]);
      } else if (me->schemaVersion==3) {
	out->rows[irow]->observingLog = ASDMparse_strarray (line, maxLine, prior, &next);
	Strip(out->rows[irow]->observingLog[0]);
      } else g_error("Unsupported schemaVersion %d", me->schemaVersion);
      continue;
    }
    prior = "<numObservingLog>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numObservingLog = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<sessionReference>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->sessionReference = ASDMparse_str (line, maxLine, prior, &next);
      /*Strip(out->rows[irow]->sessionReference); Huh? */
      continue;
    }
    prior = "<schedulerMode>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->schedulerMode = ASDMparse_str (line, maxLine, prior, &next);
      /* Strip(out->rows[irow]->schedulerMode); curious */
      continue;
    }
    prior = "<baseRangeMin>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->baseRangeMin = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<baseRangeMax>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->baseRangeMax = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<baseRmsMinor>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->baseRmsMinor = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<baseRmsMajor>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->baseRmsMajor = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<basePa>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->basePa = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<siteAltitude>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->siteAltitude = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<siteLongitude>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->siteLongitude = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<siteLatitude>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->siteLatitude = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<aborted>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->aborted = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numAntenna>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numAntenna = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<flagRow>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->flagRow = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }
    prior = "<antennaId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->antennaId = ASDMparse_enumarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<sbSummaryId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->sbSummaryId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<ScaleId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->ScaleId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Get array */
  if (out->rows[0]->telescopeName) {
    /* Is this ALMA data? */
    me->isALMA = !strncmp(out->rows[0]->telescopeName, "ALMA", 4);
    /* ALMA may be a bit confused */
    if (!me->isALMA) me->isALMA = !strncmp(out->rows[0]->telescopeName, "OSF", 3);
    if (!me->isALMA) me->isALMA = !strncmp(out->rows[0]->telescopeName, "AOS", 3);
   
    /* Is this EVLA data? */
    me->isEVLA = !strncmp(out->rows[0]->telescopeName, "EVLA", 4);
  } else { /* If nothing probably ALMA */
    me->isALMA = TRUE;
    me->isEVLA = FALSE;
  }

 /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMExecBlockTable */

/** 
 * Destructor for ExecBloc table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMExecBlockTable* KillASDMExecBlockTable(ASDMExecBlockTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMExecBlockRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMExecBlockTable */

/* ----------------------    Feed ----------------------------------- */
/** 
 * Destructor for Feed table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMFeedRow* KillASDMFeedRow(ASDMFeedRow* row)
{
  if (row == NULL) return NULL;
  if (row->timeInterval)      g_free(row->timeInterval);
  if (row->beamOffset)        g_free(row->beamOffset);
  if (row->focusReference)    g_free(row->focusReference);
  if (row->polResponse)       g_free(row->polResponse);
  if (row->polarizationTypes) g_free(row->polarizationTypes);
  if (row->receptorAngle)     g_free(row->receptorAngle);
  if (row->receiverId)        g_free(row->receiverId);
  g_free(row);
  return NULL;
} /* end   KillASDMFeedRow */

/** 
 * Constructor for Feed table parsing from file
 * \param  FeedFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMFeedTable to free
 */
static ASDMFeedTable* ParseASDMFeedTable(ObitSDMData *me, 
					 gchar *FeedFile, 
					 ObitErr *err)
{
  ASDMFeedTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong i, j, charLeft, irow, maxLine = me->maxLine, ndim, naxis;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next, *b;
  gchar *routine = " ParseASDMFeedTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMFeedTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->FeedRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMFeedRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMFeedRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, FeedFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<feedId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->feedId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<timeInterval>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->timeInterval = ASDMparse_timeRange(line, maxLine, prior, &next);
      continue;
    }
    prior = "<numReceptor>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numReceptor = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<beamOffset>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->beamOffset = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<focusReference>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->focusReference = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<polarizationTypes>";
    b = g_strstr_len (line, maxLine, prior);
    if (b!=NULL) {
      b += strlen(prior);
      /* Parse array of strings */
      /* Get dimensionality - only can handle 1 */
      ndim = (olong)strtol(b, &next, 10);
      g_assert(ndim==1);  /* bother */
      b = next;
      /* get number of values */
      naxis = (olong)strtol(b, &next, 10);
      out->rows[irow]->polarizationTypes = g_malloc0(naxis*sizeof(ObitASDMPolnType));
      b = next;
      /* Find next non blank */
      while (*b==' ') b++;
     /* Cycle through values */
      for (i=0; i<naxis; i++) {
	/* Find next blank or '<' */
	charLeft =  maxLine - (b-line);
	for (j=0; j<charLeft; j++) {
	  if ((b[j]==' ') || (b[j]=='<')) {
	    b[j] = 0; break;
	  }
	}
	out->rows[irow]->polarizationTypes[i] =  LookupPolnType(b);
	b += j+1; /* Next value */
      }
      continue;
    }
    prior = "<polResponse>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->polResponse = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<receptorAngle>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->receptorAngle = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<antennaId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Antenna_";
      out->rows[irow]->antennaId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<receiverId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->receiverId = ASDMparse_intarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<spectralWindowId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "SpectralWindow_";
      out->rows[irow]->spectralWindowId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMFeedTable */

/** 
 * Destructor for Feed table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMFeedTable* KillASDMFeedTable(ASDMFeedTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMFeedRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMFeedTable */

/* ----------------------  Field ----------------------------------- */
/** 
 * Destructor for Field table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMFieldRow* KillASDMFieldRow(ASDMFieldRow* row)
{
  if (row == NULL) return NULL;
  if (row->fieldName)    g_free(row->fieldName);
  if (row->code)         g_free(row->code);
  if (row->delayDir)     g_free(row->delayDir);
  if (row->phaseDir)     g_free(row->phaseDir);
  if (row->referenceDir) g_free(row->referenceDir);
  g_free(row);
  return NULL;
} /* end   KillASDMFieldRow */

/** 
 * Constructor for Field table parsing from file
 * \param  FieldFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMFieldTable to free
 */
static ASDMFieldTable* ParseASDMFieldTable(ObitSDMData *me, 
					   gchar *FieldFile, 
					   ObitErr *err)
{
  ASDMFieldTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMFieldTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMFieldTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->FieldRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMFieldRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMFieldRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, FieldFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<fieldId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Field_";
      out->rows[irow]->fieldId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<fieldName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->fieldName = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->fieldName);  /* Deblank */
     continue;
    }
    prior = "<code>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->code = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->code);
      continue;
    }
    prior = "<numPoly>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numPoly = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<delayDir>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->delayDir = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<phaseDir>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->phaseDir = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<referenceDir>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->referenceDir = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<time>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->time = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    prior = "<sourceId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->sourceId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMFieldTable */

/** 
 * Destructor for Field table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMFieldTable* KillASDMFieldTable(ASDMFieldTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMFieldRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMFieldTable */

/* ----------------------  Flag ----------------------------------- */
/** 
 * Destructor for Flag table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMFlagRow* KillASDMFlagRow(ASDMFlagRow* row)
{
  if (row == NULL) return NULL;
  if (row->reason) g_free(row->reason);
  if (row->antennaId) {
    g_free(row->antennaId);
  }
  g_free(row);
  return NULL;
} /* end   KillASDMFlagRow */

/** 
 * Constructor for Flag table parsing from file
 * \param  FlagFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMFlagTable to free
 */
static ASDMFlagTable* ParseASDMFlagTable(ObitSDMData *me, 
					   gchar *FlagFile, 
					   ObitErr *err)
{
  ASDMFlagTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMFlagTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMFlagTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX (0, me->ASDMTab->FlagRows);
  /* Need to count? */
  if (me->ASDMTab->FlagRows<0)
    out->nrows = CountTableRows(FlagFile, err);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMFlagRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMFlagRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, FlagFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    /* old EVLA version */
    prior = "<antennaId>Antenna_";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Antenna_";
      out->rows[irow]->antennaId    = g_malloc0(sizeof(olong*));
      out->rows[irow]->antennaId[0] = ASDMparse_int (line, maxLine, prior, &next);
      out->rows[irow]->numAntenna   = 1;
      continue;
    }
      
    /* ALMA and new EVLA version */
    prior = "<antennaId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->antennaId = ASDMparse_enumarray (line, maxLine, prior, &next);
      continue;
    }

    prior = "<numAntenna>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numAntenna = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<reason>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->reason = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->reason);
      continue;
    }
    prior = "<startTime>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->startTime = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    prior = "<endTime>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->endTime = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }

    prior = "<numPolarizationType>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numPolarizationType = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

    prior = "<numSpectralWindow>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numSpectralWindow= ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

   /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMFlagTable */

/** 
 * Destructor for Flag table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMFlagTable* KillASDMFlagTable(ASDMFlagTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMFlagRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMFlagTable */

/* ----------------------  Pointing ----------------------------------- */
/** 
 * Destructor for Pointing table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMPointingRow* KillASDMPointingRow(ASDMPointingRow* row)
{
  if (row == NULL) return NULL;
  if (row->timeInterval)      g_free(row->timeInterval);
  if (row->encoder)           g_free(row->encoder);
  if (row->pointingDirection) g_free(row->pointingDirection);
  if (row->target)            g_free(row->target);
  if (row->offset)            g_free(row->offset);
  g_free(row);
  return NULL;
} /* end   KillASDMPointingRow */

/** 
 * Constructor for Pointing table parsing from file
 * \param  PointingFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMPointingTable to free
 */
static ASDMPointingTable* ParseASDMPointingTable(ObitSDMData *me, 
						 gchar *PointingFile, 
						 ObitErr *err)
{
  ASDMPointingTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  odouble mjdJD0=2400000.5; /* JD of beginning of MJD time */
  gchar *prior, *next;
  gchar *routine = " ParseASDMPointingTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMPointingTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->PointingRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMPointingRow*));
  for (irow=0; irow<out->nrows; irow++) 
    out->rows[irow] = g_malloc0(sizeof(ASDMPointingRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, PointingFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<timeInterval>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->timeInterval = ASDMparse_timeRange(line, maxLine, prior, &next);
      /* Remove offset from second */
      if ((out->rows[irow]->timeInterval[1]<out->rows[irow]->timeInterval[0]) &&
	  (out->rows[irow]->timeInterval[1]>mjdJD0))
	out->rows[irow]->timeInterval[1] -= mjdJD0;
      continue;
    }
    prior = "<numSample>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numSample = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<encoder>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->encoder = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<pointingTracking>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->pointingTracking = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }
    prior = "<usePolynomials>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->usePolynomials = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }
    prior = "<timeOrigin>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->timeOrigin = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numTerm>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numTerm = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<pointingDirection>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->pointingDirection = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<target>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->target = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<offset>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->offset = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<overTheTop>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->overTheTop = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }
    prior = "<antennaId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Antenna_";
      out->rows[irow]->antennaId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<pointingModelId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->pointingModelId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    if (irow>out->nrows) return out;  /* Needs more work */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMPointingTable */

/** 
 * Destructor for Pointing table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMPointingTable* KillASDMPointingTable(ASDMPointingTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMPointingRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMPointingTable */

/* ----------------------  PointingModel ----------------------------------- */
/** 
 * Destructor for PointingModel table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMPointingModelRow* 
KillASDMPointingModelRow(ASDMPointingModelRow* row)
{
  if (row == NULL) return NULL;
  g_free(row);
  return NULL;
} /* end   KillASDMPointingModelRow */

/** 
 * Constructor for PointingModel table parsing from file
 * \param  PointingModelFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMPointingModelTable to free
 */
static ASDMPointingModelTable* 
ParseASDMPointingModelTable(ObitSDMData *me, 
			    gchar *PointingModelFile, 
			    ObitErr *err)
{
  ASDMPointingModelTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  /*gchar *prior, *next;*/
  gchar *routine = " ParseASDMPointingModelTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMPointingModelTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->PointingModelRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMPointingModelRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMPointingModelRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, PointingModelFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMPointingModelTable */

/** 
 * Destructor for PointingModel table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMPointingModelTable* 
KillASDMPointingModelTable(ASDMPointingModelTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMPointingModelRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMPointingModelTable */

/* ----------------------  Polarization ----------------------------------- */
/** 
 * Destructor for Polarization table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMPolarizationRow* 
KillASDMPolarizationRow(ASDMPolarizationRow* row)
{
  olong i, n;
  if (row == NULL) return NULL;
  if (row->corrType) {
    n = row->numCorr;
    for (i=0; i<n; i++) if (row->corrType[i]) g_free(row->corrType[i]);
     g_free(row->corrType);
  }
  if (row->corrProduct) {
    n = row->numCorr*2;
    for (i=0; i<n; i++) {
      if (row->corrProduct[i]==NULL) break;
      if (row->corrProduct[i]) g_free(row->corrProduct[i]);
    }
     g_free(row->corrProduct);
  }
  g_free(row);
  return NULL;
} /* end   KillASDMPolarizationRow */

/** 
 * Constructor for Polarization table parsing from file
 * \param  PolarizationFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMPolarizationTable to free
 */
static ASDMPolarizationTable* 
ParseASDMPolarizationTable(ObitSDMData *me, 
			   gchar *PolarizationFile, 
			   ObitErr *err)
{
  ASDMPolarizationTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMPolarizationTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMPolarizationTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->PolarizationRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMPolarizationRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMPolarizationRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, PolarizationFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<polarizationId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Polarization_";
      out->rows[irow]->polarizationId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numCorr>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numCorr = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<corrType>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->corrType = ASDMparse_strarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<corrProduct>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->corrProduct = ASDMparse_strarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<flagRow>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->flagRow = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMPolarizationTable */

/** 
 * Destructor for Polarization table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMPolarizationTable* KillASDMPolarizationTable(ASDMPolarizationTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMPolarizationRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMPolarizationTable */

/* ----------------------  Processor ----------------------------------- */
/** 
 * Destructor for Processor table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMProcessorRow* KillASDMProcessorRow(ASDMProcessorRow* row)
{
  if (row == NULL) return NULL;
  if (row->processorType)    g_free(row->processorType);
  if (row->processorSubType) g_free(row->processorSubType);
  g_free(row);
  return NULL;
} /* end   KillASDMProcessorRow */

/** 
 * Constructor for Processor table parsing from file
 * \param  ProcessorFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMProcessorTable to free
 */
static ASDMProcessorTable* 
ParseASDMProcessorTable(ObitSDMData *me, 
			gchar *ProcessorFile, 
			ObitErr *err)
{
  ASDMProcessorTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMProcessorTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMProcessorTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->ProcessorRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMProcessorRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMProcessorRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, ProcessorFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<processorId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Processor_";
      out->rows[irow]->processorId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<modeId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "CorrelatorMode_";
      out->rows[irow]->modeId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<processorType>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->processorType = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->processorType);
      continue;
    }
    prior = "<processorSubType>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->processorSubType = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->processorSubType);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMProcessorTable */

/** 
 * Destructor for Processor table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMProcessorTable* KillASDMProcessorTable(ASDMProcessorTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMProcessorRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMProcessorTable */

/* ---------------------- Receiver ----------------------------------- */
/** 
 * Destructor for Receiver table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMReceiverRow* KillASDMReceiverRow(ASDMReceiverRow* row)
{
  if (row == NULL) return NULL;
  g_free(row);
  return NULL;
} /* end   KillASDMReceiverRow */

/** 
 * Constructor for Receiver table parsing from file
 * \param  ReceiverFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMReceiverTable to free
 */
static ASDMReceiverTable* 
ParseASDMReceiverTable(ObitSDMData *me, 
		       gchar *ReceiverFile, 
		       ObitErr *err)
{
  ASDMReceiverTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  /*gchar *prior, *next;*/
  gchar *routine = " ParseASDMReceiverTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMReceiverTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->ReceiverRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMReceiverRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMReceiverRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, ReceiverFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMReceiverTable */

/** 
 * Destructor for Receiver table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMReceiverTable* KillASDMReceiverTable(ASDMReceiverTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMReceiverRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMReceiverTable */

 /* ----------------------  SBSummary  ----------------------------------- */
/** 
 * Destructor for SBSummary table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMSBSummaryRow* KillASDMSBSummaryRow(ASDMSBSummaryRow* row)
{
  if (row == NULL) return NULL;
  g_free(row);
  return NULL;
} /* end   KillASDMSBSummaryRow */

/** 
 * Constructor for SBSummary table parsing from file
 * \param  SBSummaryFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMSBSummaryTable to free
 */
static ASDMSBSummaryTable* 
ParseASDMSBSummaryTable(ObitSDMData *me, 
			gchar *SBSummaryFile, 
			ObitErr *err)
{
  ASDMSBSummaryTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  /*gchar *prior, *next;*/
  gchar *routine = " ParseASDMSBSummaryTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMSBSummaryTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->SBSummaryRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMSBSummaryRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMSBSummaryRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, SBSummaryFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMSBSummaryTable */

/** 
 * Destructor for SBSummary table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMSBSummaryTable* KillASDMSBSummaryTable(ASDMSBSummaryTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMSBSummaryRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMSBSummaryTable */

/* ----------------------  Scan ----------------------------------- */
/** 
 * Destructor for Scan table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMScanRow* KillASDMScanRow(ASDMScanRow* row)
{
  olong i, n;

  if (row == NULL) return NULL;
  if (row->calibrationOnLine) g_free(row->calibrationOnLine);
  if (row->sourceName)        g_free(row->sourceName);
  if (row->scanIntent) {
    n = row->numIntent;
    for (i=0; i<n; i++) {
      if (row->scanIntent[i]==NULL) break;
      if (row->scanIntent[i]) g_free(row->scanIntent[i]);
    }
    g_free(row->scanIntent);
  }
  if (row->calDataType) {
    n = row->numIntent;
    for (i=0; i<n; i++) {
      if (row->calDataType[i]==NULL) break;
      if (row->calDataType[i]) g_free(row->calDataType[i]);
    }
    g_free(row->calDataType);
  }
  g_free(row);
  return NULL;
} /* end   KillASDMScanRow */

/** 
 * Constructor for Scan table parsing from file
 * \param  ScanFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMScanTable to free
 */
static ASDMScanTable* ParseASDMScanTable(ObitSDMData *me, 
					 gchar *ScanFile, 
					 ObitErr *err)
{
  ASDMScanTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMScanTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMScanTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->ScanRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMScanRow*));
  for (irow=0; irow<out->nrows; irow++) 
    out->rows[irow] = g_malloc0(sizeof(ASDMScanRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, ScanFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<scanNumber>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->scanNumber = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<startTime>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->startTime = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    prior = "<endTime>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->endTime = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numIntent>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numIntent = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
     prior = "<numSubScan>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numSubscan = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<scanIntent>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->scanIntent = ASDMparse_strarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<calDataType>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->calDataType = ASDMparse_strarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<calibrationOnLine>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]-> calibrationOnLine= ASDMparse_booarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<sourceName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->sourceName = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->sourceName);  /* Deblank */
     continue;
    }
    prior = "<flagRow>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->flagRow = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }
    prior = "<execBlockId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "ExecBlock_";
      out->rows[irow]->execBlockId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMScanTable */

/** 
 * Destructor for Scan table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMScanTable* KillASDMScanTable(ASDMScanTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMScanRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMScanTable */

/* ----------------------  Source ----------------------------------- */
/** 
 * Destructor for Source table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMSourceRow* KillASDMSourceRow(ASDMSourceRow* row)
{
  if (row == NULL) return NULL;
  if (row->timeInterval)  g_free(row->timeInterval);
  if (row->code)          g_free(row->code);
  if (row->direction)     g_free(row->direction);
  if (row->properMotion)  g_free(row->properMotion);
  if (row->sourceName)    g_free(row->sourceName);
  if (row->restFrequency) g_free(row->restFrequency);
  if (row->sysVel)        g_free(row->sysVel);
  g_free(row);
  return NULL;
} /* end   KillASDMSourceRow */

/** 
 * Constructor for Source table parsing from file
 * \param  SourceFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMSourceTable to free
 */
static ASDMSourceTable* ParseASDMSourceTable(ObitSDMData *me, 
					     gchar *SourceFile, 
					     ObitErr *err)
{
  ASDMSourceTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMSourceTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMSourceTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->SourceRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMSourceRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMSourceRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, SourceFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<sourceId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->sourceId = ASDMparse_int (line, maxLine, prior, &next);
      /* Initial number = ID + 1 */
      out->rows[irow]->sourceNo = out->rows[irow]->sourceId + 1;
      continue;
    }
    prior = "<timeInterval>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->timeInterval = ASDMparse_timeRange(line, maxLine, prior, &next);
      continue;
    }
    prior = "<code>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->code = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->code);
      continue;
    }
    prior = "<direction>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->direction = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<properMotion>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->properMotion = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<sourceName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->sourceName = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->sourceName);  /* Deblank */
      continue;
    }
    prior = "<numLines>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numLines = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<restFrequency>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->restFrequency = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<sysVel>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->sysVel = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<spectralWindowId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "SpectralWindow_";
      out->rows[irow]->spectralWindowId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

   /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMSourceTable */

/** 
 * Destructor for Source table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMSourceTable* KillASDMSourceTable(ASDMSourceTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMSourceRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMSourceTable */

/* ----------------------  SpectralWindow ----------------------------------- */
/** 
 * Destructor for SpectralWindow table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMSpectralWindowRow* 
KillASDMSpectralWindowRow(ASDMSpectralWindowRow* row)
{
  if (row == NULL) return NULL;
  if (row->basebandName)          g_free(row->basebandName);
  if (row->netSideband)           g_free(row->netSideband);
  if (row->name)                  g_free(row->name);
  if (row->correlationBit)        g_free(row->correlationBit);
  if (row->chanWidthArray)        g_free(row->chanWidthArray);
  if (row->effectiveBwArray)      g_free(row->effectiveBwArray);
  if (row->resolutionArray)       g_free(row->resolutionArray);
  if (row->SpecRes)               g_free(row->SpecRes);
  if (row->assocSpectralWindowId) g_free(row->assocSpectralWindowId);
  g_free(row);
  return NULL;
} /* end   KillASDMSpectralWindowRow */

/** 
 * Constructor for SpectralWindow table parsing from file
 * \param  SpectralWindowFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMSpectralWindowTable to free
 */
static ASDMSpectralWindowTable* 
ParseASDMSpectralWindowTable(ObitSDMData *me, 
			     gchar *SpectralWindowFile, 
			     ObitErr *err)
{
  ASDMSpectralWindowTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong i, irow, maxLine = me->maxLine;
  gchar *line=me->line, **assNat=NULL;
  gchar *endrow = "</row>";
  gchar *prior, *next, *tstr;
  gchar *routine = " ParseASDMSpectralWindowTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMSpectralWindowTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->SpectralWindowRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMSpectralWindowRow*));
  for (irow=0; irow<out->nrows; irow++) 
    out->rows[irow] = g_malloc0(sizeof(ASDMSpectralWindowRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, SpectralWindowFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<spectralWindowId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "SpectralWindow_";
      out->rows[irow]->spectralWindowId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<basebandName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->basebandName = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->basebandName);
      continue;
    }
    prior = "<netSideband>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->netSideband = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->netSideband);
      continue;
    }
    prior = "<numChan>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numChan = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<refFreq>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->refFreq = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<sidebandProcessingMode>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      tstr =  ASDMparse_str (line, maxLine, prior, &next);
      Strip(tstr);
      out->rows[irow]->sidebandProcessingMode = LookupSideBMode(tstr);
      g_free(tstr);
      continue;
    }
    prior = "<totBandwidth>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->totBandwidth = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<windowFunction>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      tstr =  ASDMparse_str (line, maxLine, prior, &next);
      Strip(tstr);
      out->rows[irow]->windowFunction = LookupWindowFn(tstr);
      g_free(tstr);
      continue;
    }
    prior = "<chanFreqStart>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->chanFreqStart = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<chanFreqStep>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->chanFreqStep = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<chanWidth>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->chanWidth = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<correlationBit>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->correlationBit = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->correlationBit);
      continue;
    }
    prior = "<effectiveBw>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->effectiveBw = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }
    prior = "<name>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->name = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->name);
      continue;
    }
    prior = "<oversampling>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->oversampling = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }
    prior = "<quantization>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->quantization = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }
    prior = "<resolution>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->resolution = ASDMparse_dbl (line, maxLine, prior, &next);
      continue;
    }

    prior = "<chanFreqArray>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->chanFreqArray = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<chanWidthArray>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->chanWidthArray = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<effectiveBwArray>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->effectiveBwArray = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<resolutionArray>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->resolutionArray = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numAssocValues>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numAssocValues = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<assocNature>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      assNat = ASDMparse_strarray (line, maxLine, prior, &next);
      out->rows[irow]->SpecRes = g_malloc0(MIN(1,out->rows[irow]->numAssocValues)*sizeof(ObitASDMSpecRes));
      i = 0;
      while(assNat[i]) {
	Strip(assNat[i]);
	if (i>=MIN(1,out->rows[irow]->numAssocValues)) break;  /* FULL? */
	if (!strcmp(assNat[i], "CHANNEL_AVERAGE"))      
	  out->rows[irow]->SpecRes[i] = (olong)ASDMSpecRes_CHANNEL_AVERAGE;
	else if (!strcmp(assNat[i], "BASEBAND_WIDE"))   
	  out->rows[irow]->SpecRes[i] = (olong)ASDMSpecRes_BASEBAND_WIDE;
	else if (!strcmp(assNat[i], "FULL_RESOLUTION")) 
	  out->rows[irow]->SpecRes[i] = (olong)ASDMSpecRes_FULL_RESOLUTION;
	g_free(assNat[i]);
	i++;
      }
      g_free(assNat);
      continue; 
    }
    
    prior = "<assocSpectralWindowId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->assocSpectralWindowId = ASDMparse_enumarray (line, maxLine, prior, &next);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMSpectralWindowTable */

/** 
 * Destructor for SpectralWindow table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMSpectralWindowTable* 
KillASDMSpectralWindowTable(ASDMSpectralWindowTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMSpectralWindowRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMSpectralWindowTable */

/* ----------------------  State ----------------------------------- */
/** 
 * Destructor for State table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMStateRow* KillASDMStateRow(ASDMStateRow* row)
{
  if (row == NULL) return NULL;
  if (row->calDeviceName) g_free(row->calDeviceName);
  g_free(row);
  return NULL;
} /* end   KillASDMStateRow */

/** 
 * Constructor for State table parsing from file
 * \param  StateFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMStateTable to free
 */
static ASDMStateTable* ParseASDMStateTable(ObitSDMData *me, 
					   gchar *StateFile, 
					   ObitErr *err)
{
  ASDMStateTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMStateTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMStateTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->StateRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMStateRow*));
  for (irow=0; irow<out->nrows; irow++) 
    out->rows[irow] = g_malloc0(sizeof(ASDMStateRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, StateFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
     prior = "<stateId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "State_";
      out->rows[irow]->stateId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<calDeviceName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->calDeviceName = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->calDeviceName);
      continue;
    }
    prior = "<sig>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->sig = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }
    prior = "<ref>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->ref = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }
    prior = "<onSky>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->onSky = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMStateTable */

/** 
 * Destructor for State table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMStateTable* KillASDMStateTable(ASDMStateTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMStateRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMStateTable */

/* ----------------------  Station ----------------------------------- */
/** 
 * Destructor for Station table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMStationRow* KillASDMStationRow(ASDMStationRow* row)
{
  if (row == NULL) return NULL;
  if(row->name)     g_free(row->name);
  if(row->position) g_free(row->position);
  g_free(row);
  return NULL;
} /* end   KillASDMStationRow */

/** 
 * Constructor for Station table parsing from file
 * \param  StationFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMStationTable to free
 */
static ASDMStationTable* ParseASDMStationTable(ObitSDMData *me, 
					       gchar *StationFile, 
					       ObitErr *err)
{
  ASDMStationTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next, *tstr;
  gchar *routine = " ParseASDMStationTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMStationTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->StationRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMStationRow*));
  for (irow=0; irow<out->nrows; irow++) 
    out->rows[irow] = g_malloc0(sizeof(ASDMStationRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, StationFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    /* Station ID */
    prior = "<stationId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Station_";
      out->rows[irow]->stationId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

    /* name */
    prior = "<name>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->name = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->name);
      continue;
    }

    /* position array of doubles */
    prior = "<position>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->position = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    
    /* station Type enum */
    prior = "<type>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      tstr =  ASDMparse_str (line, maxLine, prior, &next);
      Strip(tstr);
      if (!strcmp(tstr, "ANTENNA_PAD"))          
	out->rows[irow]->type = ASDMStn_ANTENNA_PAD;
      else if (!strcmp(tstr, "MAINTENANCE_PAD")) 
	out->rows[irow]->type = ASDMStn_MAINTENANCE_PAD;
      else if (!strcmp(tstr, "WEATHER_STATION")) 
	out->rows[irow]->type = ASDMStn_WEATHER_STATION;
      g_free(tstr);
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMStationTable */

/** 
 * Destructor for Station table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMStationTable* KillASDMStationTable(ASDMStationTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMStationRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMStationTable */

/* ---------------------- Subscan ----------------------------------- */
/** 
 * Destructor for Subscan table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMSubscanRow* KillASDMSubscanRow(ASDMSubscanRow* row)
{
  if (row == NULL) return NULL;
  if (row->fieldName)            g_free(row->fieldName);
  if (row->subscanIntent)        g_free(row->subscanIntent);
  if (row->numberSubintegration) g_free(row->numberSubintegration);
  g_free(row);
  return NULL;
} /* end   KillASDMSubscanRow */

/** 
 * Constructor for Subscan table parsing from file
 * \param  SubscanFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMSubscanTable to free
 */
static ASDMSubscanTable* 
ParseASDMSubscanTable(ObitSDMData *me, 
		      gchar *SubscanFile, 
		      ObitErr *err)
{
  ASDMSubscanTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMSubscanTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMSubscanTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->SubscanRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMSubscanRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMSubscanRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, SubscanFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<scanNumber>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->scanNumber = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<subscanNumber>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->subscanNumber = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<startTime>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->startTime = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    prior = "<endTime>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->endTime = ASDMparse_time (line, maxLine, prior, &next);
      continue;
    }
    prior = "<fieldName>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->fieldName = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->fieldName);  /* Deblank */
      continue;
    }
    prior = "<subscanIntent>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->subscanIntent = ASDMparse_str (line, maxLine, prior, &next);
      Strip(out->rows[irow]->subscanIntent);
      continue;
    }
    prior = "<numberIntegration>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numberIntegration = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numberSubintegration>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numberSubintegration = ASDMparse_intarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<flagRow>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->flagRow = ASDMparse_boo (line, maxLine, prior, &next);
      continue;
    }
    prior = "<execBlockId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "ExecBlock_";
      out->rows[irow]->execBlockId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMSubscanTable */

/** 
 * Destructor for Subscan table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMSubscanTable* KillASDMSubscanTable(ASDMSubscanTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMSubscanRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMSubscanTable */

/* ---------------------- SwitchCycle ----------------------------------- */
/** 
 * Destructor for SwitchCycle table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMSwitchCycleRow* KillASDMSwitchCycleRow(ASDMSwitchCycleRow* row)
{
  if (row == NULL) return NULL;
  if (row->weightArray)       g_free(row->weightArray);
  if (row->dirOffsetArray)    g_free(row->dirOffsetArray);
  if (row->freqOffsetArray)   g_free(row->freqOffsetArray);
  if (row->stepDurationArray) g_free(row->stepDurationArray);
  g_free(row);
  return NULL;
} /* end   KillASDMSwitchCycleRow */

/** 
 * Constructor for SwitchCycle table parsing from file
 * \param  SwitchCycleFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMSwitchCycleTable to free
 */
static ASDMSwitchCycleTable* 
ParseASDMSwitchCycleTable(ObitSDMData *me, 
			  gchar *SwitchCycleFile, 
			  ObitErr *err)
{
  ASDMSwitchCycleTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMSwitchCycleTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMSwitchCycleTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->SwitchCycleRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMSwitchCycleRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMSwitchCycleRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, SwitchCycleFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<switchCycleId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "SwitchCycle_";
      out->rows[irow]->switchCycleId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numStep>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numStep = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<weightArray>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->weightArray = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<dirOffsetArray>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->dirOffsetArray = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<freqOffsetArray>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->freqOffsetArray = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<stepDurationArray>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->stepDurationArray = ASDMparse_dblarray (line, maxLine, prior, &next);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMSwitchCycleTable */

/** 
 * Destructor for SwitchCycle table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMSwitchCycleTable* KillASDMSwitchCycleTable(ASDMSwitchCycleTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMSwitchCycleRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMSwitchCycleTable */

/* ---------------------- SysCal ----------------------------------- */
/** 
 * Destructor for SysCal table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMSysCalRow* KillASDMSysCalRow(ASDMSysCalRow* row)
{
  if (row == NULL) return NULL;
  g_free(row);
  return NULL;
} /* end   KillASDMSysCalRow */

/** 
 * Constructor for SysCal table parsing from file
 * \param  SysCalFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMSysCalTable to free
 */
static ASDMSysCalTable* 
ParseASDMSysCalTable(ObitSDMData *me, 
		     gchar *SysCalFile, 
		     ObitErr *err)
{
  ASDMSysCalTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  /*gchar *prior, *next;*/
  gchar *routine = " ParseASDMSysCalTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMSysCalTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->SysCalRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMSysCalRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMSysCalRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, SysCalFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMSysCalTable */

/** 
 * Destructor for SysCal table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMSysCalTable* KillASDMSysCalTable(ASDMSysCalTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMSysCalRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMSysCalTable */

/* ---------------------- SysPower ----------------------------------- */
/** 
 * Destructor for SysPower table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMSysPowerRow* KillASDMSysPowerRow(ASDMSysPowerRow* row)
{
  if (row == NULL) return NULL;
  if (row->timeInterval)            g_free(row->timeInterval);
  if (row->switchedPowerDifference) g_free(row->switchedPowerDifference);
  if (row->switchedPowerSum)        g_free(row->switchedPowerSum);
  if (row->requantizerGain)         g_free(row->requantizerGain);
  g_free(row);
  return NULL;
} /* end   KillASDMSysPowerRow */

/** 
 * Constructor for SysPower table parsing from file
 * \param  SysPowerFile Name root of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMSysPowerTable to free
 */
static ASDMSysPowerTable* 
ParseASDMSysPowerTable(ObitSDMData *me, 
		       gchar *SysPowerFile, 
		       ObitErr *err)
{
  ASDMSysPowerTable* out=NULL;
  gchar *fullname;
  gchar *routine = "ParseASDMSysPowerTable";

  if (err->error) return out;

  /* does either table exist? Try binary first */
  fullname = g_strconcat (SysPowerFile,".bin", NULL);
  if (ObitFileExist(fullname, err)) {
    out = ParseASDMSysPowerTableBin(me, fullname, err);
  }
  if (fullname) g_free(fullname);
  if (err->error) Obit_traceback_val (err, routine, me->name, out);
  if (out) return out;  /* Done? */

  /* does either table exist? Try binary first */
  fullname = g_strconcat (SysPowerFile,".xml", NULL);
  if (ObitFileExist(fullname, err)) {
    out = ParseASDMSysPowerTableXML(me, fullname, err);
  }
  if (fullname) g_free(fullname);
  if (err->error) Obit_traceback_val (err, routine, me->name, out);
  return out;
} /* end ParseASDMSysPowerTable */

/** 
 * Constructor for SysPower table parsing from xml file
 * \param  SysPowerFile Name of file containing table
 * \param  err          ObitErr for reporting errors.
 * \return table structure,  use KillASDMSysPowerTable to free
 */
static ASDMSysPowerTable* 
ParseASDMSysPowerTableXML(ObitSDMData *me, 
			  gchar *SysPowerFile, 
			  ObitErr *err)
{
  ASDMSysPowerTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  odouble mjdJD0=2400000.5; /* JD of beginning of MJD time */
  gchar *line=me->line;
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = "ParseASDMSysPowerTableXML";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMSysPowerTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->SysPowerRows);
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMSysPowerRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMSysPowerRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, SysPowerFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<antennaId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Antenna_";
      out->rows[irow]->antennaId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<spectralWindowId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "SpectralWindow_";
      out->rows[irow]->spectralWindowId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<feedId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->feedId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<numReceptor>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->numReceptor = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    prior = "<timeInterval>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->timeInterval = ASDMparse_timeRange(line, maxLine, prior, &next);
      /* Remove offset from second */
      if ((out->rows[irow]->timeInterval[1]<out->rows[irow]->timeInterval[0]) &&
	  (out->rows[irow]->timeInterval[1]>mjdJD0))
	out->rows[irow]->timeInterval[1] -= mjdJD0;
      continue;
    }
    prior = "<switchedPowerDifference>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->switchedPowerDifference = ASDMparse_fltarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<switchedPowerSum>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->switchedPowerSum = ASDMparse_fltarray (line, maxLine, prior, &next);
      continue;
    }
    prior = "<requantizerGain>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->requantizerGain = ASDMparse_fltarray (line, maxLine, prior, &next);
      continue;
    }

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMSysPowerTable */

/** 
 * Constructor for SysPower table parsing from binary file
 * \param  SysPowerFile Name of file containing table
 * \param  err          ObitErr for reporting errors.
 * \return table structure,  use KillASDMSysPowerTable to free
 */
static ASDMSysPowerTable* 
ParseASDMSysPowerTableBin(ObitSDMData *me, 
			  gchar *SysPowerFile, 
			  ObitErr *err)
{
  ASDMSysPowerTable* out=NULL;
  ObitEVLASysPower*  binTab=NULL;
  olong nrow, orow, irow;
  gchar *routine = "ParseASDMSysPowerTableBin";

  if (err->error) return out;

  /* Basic output */
  out = g_malloc0(sizeof(ASDMSysPowerTable));
  out->rows = NULL;

  /* Create binary table parsing structure */
  binTab = ObitEVLASysPowerCreate ("EVLASysPower", err);
  if (err->error) Obit_traceback_val (err, routine, me->name, out);

  /* Init file */
  ObitEVLASysPowerInitFile (binTab, SysPowerFile, err);
  if (err->error) goto cleanup;

  /* How big */
  nrow = ObitEVLASysPowerGetNrow (binTab);
  /* If not given use ASDM table */
  if (nrow<0) {
    nrow = me->ASDMTab->SysPowerRows;
    binTab->nrow = nrow;
  }
  out->nrows = nrow;
  if (nrow<1) return out;

  /* Finish building output */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMSysPowerRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMSysPowerRow));

  /* Loop reading rows */
  for (irow=0; irow<out->nrows; irow++) {
    orow = ObitEVLASysPowerGetRow (binTab, out->rows[irow], err);
    if (orow<=0) break;  /* Done? */
    if (err->error) goto cleanup;
    nrow = orow;
  } /* end loop reading */
  out->nrows = nrow;

  /* Cleanup */
  binTab = ObitEVLASysPowerUnref(binTab); 
  if (err->error) Obit_traceback_val (err, routine, me->name, out);
 cleanup:
  

  if (err->error) Obit_traceback_val (err, routine, me->name, out);
  return out;
} /* end ParseASDMSysPowerTableBin */
/** 
 * Destructor for SysPower table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMSysPowerTable* KillASDMSysPowerTable(ASDMSysPowerTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMSysPowerRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMSysPowerTable */

/* ---------------------- Weather ----------------------------------- */
/** 
 * Destructor for Weather table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMWeatherRow* KillASDMWeatherRow(ASDMWeatherRow* row)
{
  if (row == NULL) return NULL;
  if (row->timeInterval)  g_free(row->timeInterval);
  g_free(row);
  return NULL;
} /* end   KillASDMWeatherRow */

/** 
 * Constructor for Weather table parsing from file
 * Note: Fblank used to indicate invalid values.
 * \param  WeatherFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMWeatherTable to free
 */
static ASDMWeatherTable* ParseASDMWeatherTable(ObitSDMData *me, 
					       gchar *WeatherFile, 
					       ObitErr *err)
{
  ASDMWeatherTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = 1024;
  ofloat fblank = ObitMagicF();
  odouble dtemp;
  gboolean flagged;
  gchar line[1025];
  gchar *endrow = "</row>";
  gchar *prior, *next;
  gchar *routine = " ParseASDMWeatherTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMWeatherTable));
  out->rows = NULL;

  /* How many rows? */
  out->nrows = MAX(0, me->ASDMTab->WeatherRows); 
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMWeatherRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMWeatherRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, WeatherFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */
    prior = "<timeInterval>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      out->rows[irow]->timeInterval = ASDMparse_timeRange(line, maxLine, prior, &next);
      continue;
    }

    prior = "<pressure>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      dtemp = ASDMparse_dbl (line, maxLine, prior, &next);
      if (out->rows[irow]->pressure!=fblank)
	out->rows[irow]->pressure = (ofloat)dtemp;
      continue;
    }
    prior = "<pressureFlag>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      flagged = ASDMparse_boo (line, maxLine, prior, &next);
      if (flagged) out->rows[irow]->pressure = fblank;
      continue;
    }

    prior = "<relHumidity>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      dtemp = ASDMparse_dbl (line, maxLine, prior, &next);
      if (out->rows[irow]->relHumidity!=fblank)
	out->rows[irow]->relHumidity = (ofloat)dtemp;
      continue;
    }
    prior = "<relHumidityFlag>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      flagged = ASDMparse_boo (line, maxLine, prior, &next);
      if (flagged) out->rows[irow]->relHumidity = fblank;
      continue;
    }
 
    prior = "<temperature>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      dtemp = ASDMparse_dbl (line, maxLine, prior, &next);
      if (out->rows[irow]->temperature!=fblank)
	out->rows[irow]->temperature = (ofloat)dtemp;
      continue;
    }
    prior = "<temperatureFlag>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      flagged = ASDMparse_boo (line, maxLine, prior, &next);
      if (flagged) out->rows[irow]->temperature = fblank;
      continue;
    }

    prior = "<windDirection>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      dtemp = ASDMparse_dbl (line, maxLine, prior, &next);
      if (out->rows[irow]->windDirection!=fblank)
	out->rows[irow]->windDirection = (ofloat)dtemp;
      continue;
    }
    prior = "<windDirectionFlag>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      flagged = ASDMparse_boo (line, maxLine, prior, &next);
      if (flagged) out->rows[irow]->windDirection = fblank;
      continue;
    }

    prior = "<windSpeed>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      dtemp = ASDMparse_dbl (line, maxLine, prior, &next);
      if (out->rows[irow]->windSpeed!=fblank)
	out->rows[irow]->windSpeed = (ofloat)dtemp;
      continue;
    }
    prior = "<windSpeedFlag>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      flagged = ASDMparse_boo (line, maxLine, prior, &next);
      if (flagged) out->rows[irow]->windSpeed = fblank;
      continue;
    }

    prior = "<windMax>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      dtemp = ASDMparse_dbl (line, maxLine, prior, &next);
      if (out->rows[irow]->windMax!=fblank)
	out->rows[irow]->windMax = (ofloat)dtemp;
      continue;
    }
    prior = "<windMaxFlag>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      flagged = ASDMparse_boo (line, maxLine, prior, &next);
      if (flagged) out->rows[irow]->windMax = fblank;
      continue;
    }

    prior = "<dewPoint>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      dtemp = ASDMparse_dbl (line, maxLine, prior, &next);
      if (out->rows[irow]->dewPoint!=fblank)
	out->rows[irow]->dewPoint = (ofloat)dtemp;
      continue;
    }
    prior = "<dewPointFlag>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      flagged = ASDMparse_boo (line, maxLine, prior, &next);
      if (flagged) out->rows[irow]->dewPoint = fblank;
      continue;
    }

    prior = "<stationId>";
    if (g_strstr_len (line, maxLine, prior)!=NULL) {
      prior = "Station_";
      out->rows[irow]->stationId = ASDMparse_int (line, maxLine, prior, &next);
      continue;
    }
    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMWeatherTable */

/** 
 * Destructor for Weather table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMWeatherTable* KillASDMWeatherTable(ASDMWeatherTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMWeatherRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMWeatherTable */

/* ---------------------- XXXX ----------------------------------- */
/** 
 * Destructor for XXXX table row.
 * \param  structure to destroy
 * \return NULL row pointer
 */
static ASDMXXXXRow* KillASDMXXXXRow(ASDMXXXXRow* row)
{
  if (row == NULL) return NULL;
  g_free(row);
  return NULL;
} /* end   KillASDMXXXXRow */

/** 
 * Constructor for XXXX table parsing from file
 * \param  XXXXFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return table structure,  use KillASDMXXXXTable to free
 */
static ASDMXXXXTable* ParseASDMXXXXTable(ObitSDMData *me, 
					 gchar *XXXXFile, 
					 ObitErr *err)
{
  ASDMXXXXTable* out=NULL;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong irow, maxLine = me->maxLine;
  gchar *line=me->line;
  gchar *endrow = "</row>";
  /*gchar *prior, *next;*/
  gchar *routine = " ParseASDMXXXXTable";

  /* error checks */
  if (err->error) return out;

  out = g_malloc0(sizeof(ASDMXXXXTable));
  out->rows = NULL;

  /* How many rows? */
  /*out->nrows = MAX(0, me->ASDMTab->XXXXRows); */
  if (out->nrows<1) return out;

  /* Finish building it */
  out->rows = g_malloc0((out->nrows+1)*sizeof(ASDMXXXXRow*));
  for (irow=0; irow<out->nrows; irow++) out->rows[irow] = g_malloc0(sizeof(ASDMXXXXRow));

  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, XXXXFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  irow = 0;
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* Parse entries */

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) irow++;

    /* Check overflow */
    Obit_retval_if_fail((irow<=out->nrows), err, out,
			"%s: Found more rows than allocated (%d)", 
			routine, out->nrows);
  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end ParseASDMXXXXTable */

/** 
 * Destructor for XXXX table
 * \param  structure to destroy
 * \return NULL pointer
 */
static ASDMXXXXTable* KillASDMXXXXTable(ASDMXXXXTable* table)
{
  olong i;

  if (table==NULL) return NULL;  /* Anybody home? */

  /* Delete row structures */
  if (table->rows) {
    for (i=0; i<table->nrows; i++) 
      table->rows[i] = KillASDMXXXXRow(table->rows[i]);
    g_free(table->rows);
  }
  g_free(table);
  return NULL;
} /* end KillASDMXXXXTable */

/** 
 * Count the number of rows in the xml table in file CntFile 
 * \param  CntFile Name of file containing table
 * \param  err     ObitErr for reporting errors.
 * \return Count of the lines which contain "<\row>"
 */
static olong CountTableRows(gchar *CntFile, ObitErr *err)
{
  olong out = 0;
  ObitFile *file=NULL;
  ObitIOCode retCode;
  olong maxLine = 1028;
  gchar line[1032];
  gchar *endrow = "</row>";
  gchar *routine = "CountTableRows";

  /* error checks */
  if (err->error) return out;

  /* Check if it exists */
  if (!ObitFileExist(CntFile, err)) return out;

  /* Create file structure */
  file = newObitFile("ASDM");
  retCode = ObitFileOpen(file, CntFile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);

  /* Loop over file */
  while (retCode!=OBIT_IO_EOF) {

    retCode = ObitFileReadXML (file, line, maxLine, err);
    if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
    if (retCode==OBIT_IO_EOF) break;

    /* DEBUG 
    fprintf (stdout, "%s",line);*/

    /* Is this the end of a row? */
    if (g_strstr_len (line, maxLine, endrow)!=NULL) 
      out++;

  } /* end loop over table */

  /* Close up */
  retCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_val (err, routine, file->fileName, out);
  file = ObitFileUnref(file);

  return out;
} /* end CountTableRows */

/**
 * Compare frequencies, to give ascending order.
 * Conformant to function type GCompareDataFunc
 * \param in1   First list
 * \param in2   Second list
 * \param ncomp Number of values to compare (1)
 * \return <0 -> in1 > in2; =0 -> in1 == in2; >0 -> in1 < in2; 
 */
static gint CompareFreq (gconstpointer in1, gconstpointer in2, 
			 gpointer ncomp)
{
  gint out = 0;
  ofloat *float1, *float2;

  /* get correctly typed local values */
  float1 = (float*)(in1 + sizeof(ofloat));
  float2 = (float*)(in2 + sizeof(ofloat));
  if ((*float1)>(*float2))      out =  1;
  else if ((*float1)<(*float2)) out = -1;
  else                          out =  0;
  return out;
} /* end CompareFreq */

/**
 * Squeeze all blanks out of a string
 * \param String to squeeze
 */
static void Strip (gchar* s)
{
  olong n, i;

  if (s==NULL) return;  /* Uh oh */
  n = strlen(s);

  /* Leading blanks */
  while (s[0]==' ') {
    for (i=0; i<n-1; i++) s[i] = s[i+1];
    n--;
  }
  /* Trailing blanks */
   for (i=n-1; i>0; i--) {
    if (s[i]==' ') s[i] = 0;
    if (s[i]!=0) break;
  }
 
} /* end Strip */

