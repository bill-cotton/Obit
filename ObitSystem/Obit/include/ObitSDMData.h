/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010,2012                                          */
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
#ifndef OBITSDMDATA_H 
#define OBITSDMDATA_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitTableSN.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSDMData.h
 *
 * This class accesses data in the EVLA BDF format
 *
 * Class documentation should go here.
 * 
 * \section ObitSDMDataaccess Creators and Destructors
 * An ObitSDMData will usually be created using ObitSDMDataCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitSDMData should always be made using the
 * #ObitSDMDataRef function which updates the reference count in the object.
 * Then whenever freeing an ObitSDMData or changing a pointer, the function
 * #ObitSDMDataUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitBDFSampleType
 * enum for BDF Sampling type
 */
enum obitBDFSampleType {
  /** simple integration */
  ASDMINTEGRATION,
  /** sub integration */
  ASDMSUBINTEGRATION
}; /* end enum obitBDFSampleType */
/** typedef for enum for BDFSampleType. */
typedef enum obitBDFSampleType ObitSDMSampleType;

/**
 * \enum obitASDMAntennaMake
 * enum for Antenna Make
 */
enum obitASDMAntennaMake {
  ASDMAnt_UNKNOWN,
  ASDMAnt_MITSUBISHI_12_A,
  ASDMAnt_MITSUBISHI_12_B,
  ASDMAnt_VERTEX_12_ATF,
  ASDMAnt_AEM_12_ATF,
  ASDMAnt_VERTEX_12,
  ASDMAnt_AEM_12,
  ASDMAnt_IRAM_15
}; /* end enum obitASDMAntennaMake */
/** typedef for enum for ASDMAntennaMake. */
typedef enum obitASDMAntennaMake ObitASDMAntennaMake;

/**
 * \enum obitASDMAntennaType
 * enum for Antenna Type
 */
enum obitASDMAntennaType {
  ASDMAnt_GROUND_BASED,
  ASDMAnt_SPACE_BASED,
  ASDMAnt_TRACKING_STN
}; /* end enum obitASDMAntennaType */
/** typedef for enum for ASDMAntennaType. */
typedef enum obitASDMAntennaType ObitASDMAntennaType;

/**
 * \enum obitASDMStationType
 * enum for Station Type
 */
enum obitASDMStationType {
  ASDMStn_ANTENNA_PAD,
  ASDMStn_MAINTENANCE_PAD,
  ASDMStn_WEATHER_STATION
}; /* end enum obitASDMStationType */
/** typedef for enum for ASDMStationType. */
typedef enum obitASDMStationType ObitASDMStationType;

/**
 * \enum obitASDMCorrMode
 * Enum for Correlation Mode
 */
enum obitASDMCorrMode {
  ASDMCorrMode_CROSS_ONLY,
  ASDMCorrMode_AUTO_ONLY,
  ASDMCorrMode_CROSS_AND_AUTO
}; /* end enum obitASDMCorrMode */
/** typedef for enum for ASDMCorrMode. */
typedef enum obitASDMCorrMode ObitASDMCorrMode;

/**
 * \enum obitASDMAtmPhCorr
 * Enum for Atmospheric phase correction Mode
 */
enum obitASDMAtmPhCorr {
  ASDMAtmPhCorr_AP_UNCORRECTED,
  ASDMAtmPhCorr_AP_CORRECTED
}; /* end enum obitASDMAtmPhCorr */
/** typedef for enum for ASDMAtmPhCorr. */
typedef enum obitASDMAtmPhCorr ObitASDMAtmPhCorr;

/**
 * \enum obitASDMProcrType
 * Enum for Processor type
 */
enum obitASDMProcrType {
  ASDMProcrType_CORRELATOR,
  ASDMProcrType_RADIOMETER,
  ASDMProcrType_SPECTROMETER
}; /* end enum obitASDMProcrType */
/** typedef for enum for ASDMProcrType. */
typedef enum obitASDMProcrType ObitASDMProcrType;

/**
 * \enum obitASDMSpecRes
 * Enum for spectral resolution type
 */
enum obitASDMSpecRes {
  ASDMSpecRes_CHANNEL_AVERAGE,
  ASDMSpecRes_BASEBAND_WIDE,
  ASDMSpecRes_FULL_RESOLUTION
}; /* end enum obitASDMSpecRes */
/** typedef for enum for ASDMSpecRes. */
typedef enum obitASDMSpecRes ObitASDMSpecRes;

/**
 * \enum obitASDMAccumMode
 * Enum for Accumulation mode
 */
enum obitASDMAccumMode {
  ASDMAccumMode_FAST,
  ASDMAccumMode_NORMAL,
  ASDMAccumMode_UNDEFINED
}; /* end enum obitASDMAccumMode */
/** typedef for enum for ASDMAccumMode. */
typedef enum obitASDMAccumMode ObitASDMAccumMode;

/**
 * \enum obitASDMBasebandName
 * Enum for Accumulation mode
 */
enum obitASDMBasebandName {
  ASDMBasebandName_NOBB,
  ASDMBasebandName_BB_1,
  ASDMBasebandName_BB_2,
  ASDMBasebandName_BB_3,
  ASDMBasebandName_BB_4,
  ASDMBasebandName_BB_5,
  ASDMBasebandName_BB_6,
  ASDMBasebandName_BB_7,
  ASDMBasebandName_BB_8,
  ASDMBasebandName_BB_ALL
}; /* end enum obitASDMBasebandName */
/** typedef for enum for ASDMBasebandName. */
typedef enum obitASDMBasebandName ObitASDMBasebandName;

/**
 * \enum obitASDMAxisName
 * Enum for Axis name
 */
enum obitASDMAxisName {
  ASDMAxisName_TIM,
  ASDMAxisName_BAL,
  ASDMAxisName_ANT,
  ASDMAxisName_BAB,
  ASDMAxisName_SPW,
  ASDMAxisName_SIB,
  ASDMAxisName_SUB,
  ASDMAxisName_BIN,
  ASDMAxisName_APC,
  ASDMAxisName_SPP,
  ASDMAxisName_POL,
  ASDMAxisName_STO,
  ASDMAxisName_HOL
}; /* end enum obitASDMAxisName */
/** typedef for enum for ASDMAxisName. */
typedef enum obitASDMAxisName ObitASDMAxisName;

/**
 * \enum obitASDMFilterMode
 * Enum for Filterulation mode
 */
enum obitASDMFilterMode {
  ASDMFilterMode_FILTER_NA,
  ASDMFilterMode_FILTER_TDM,
  ASDMFilterMode_FILTER_TFB,
  ASDMFilterMode_UNDEFINED
}; /* end enum obitASDMFilterMode */
/** typedef for enum for ASDMFilterMode. */
typedef enum obitASDMFilterMode ObitASDMFilterMode;

/**
 * \enum obitASDMPolnType
 * Enum for Polarization type
 */
enum obitASDMPolnType {
  ASDMPolnType_R,
  ASDMPolnType_L,
  ASDMPolnType_X,
  ASDMPolnType_Y
}; /* end enum obitASDMPolnType */
/** typedef for enum for ASDMPolnType. */
typedef enum obitASDMPolnType ObitASDMPolnType;

/**
 * \enum obitASDMSideBMode
 * Enum for sideband processing mode
 */
enum obitASDMSideBMode {
  ASDMSideBMode_NONE,
  ASDMSideBMode_PHASE_SWITCH_SEPARATION,
  ASDMSideBMode_FREQUENCY_OFFSET_SEPARATION,
  ASDMSideBMode_PHASE_SWITCH_REJECTION,
  ASDMSideBMode_FREQUENCY_OFFSET_REJECTION
}; /* end enum obitASDMSideBMode */
/** typedef for enum for ASDMSideBMode. */
typedef enum obitASDMSideBMode ObitASDMSideBMode;

/**
 * \enum obitASDMWindowFn
 * Enum for windowing function
 */
enum obitASDMWindowFn {
  ASDMWindowFn_UNIFORM,
  ASDMWindowFn_HANNING,
  ASDMWindowFn_HAMMING,
  ASDMWindowFn_BARTLETT,
  ASDMWindowFn_BLACKMANN,
  ASDMWindowFn_BLACKMANN_HARRIS,
  ASDMWindowFn_WELCH
}; /* end enum obitASDMWindowFn */
/** typedef for enum for ASDMWindowFn. */
typedef enum obitASDMWindowFn ObitASDMWindowFn;

/**
 * \enum obitASDMBand
 * Enum for Band codes
 */
enum obitASDMBand {
  ASDMBand_Any,
  ASDMBand_4,
  ASDMBand_P,
  ASDMBand_L,
  ASDMBand_S,
  ASDMBand_C,
  ASDMBand_X,
  ASDMBand_Ku,
  ASDMBand_K,
  ASDMBand_Ka,
  ASDMBand_Q,
  ASDMBand_A3,
  ASDMBand_A4,
  ASDMBand_A5,
  ASDMBand_A6,
  ASDMBand_A7,
  ASDMBand_A8,
  ASDMBand_A9,
  ASDMBand_A10,
  ASDMBand_A11
}; /* end enum obitASDMBand */
/** typedef for enum for ASDMBand. */
typedef enum obitASDMBand ObitASDMBand;

/*--------------Class definitions-------------------------------------*/
/** ObitSDMData ASDM/BDF structures. */
#include "ObitSDMDataTypeDef.h" 

/** ObitSDMData Class structures. */
typedef struct {
#include "ObitSDMDataDef.h"   /* this class definition */
} ObitSDMData;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSDMData
 * returns a ObitSDMData*.
 * in = object to unreference
 */
#define ObitSDMDataUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSDMData.
 * returns a ObitSDMData*.
 * in = object to reference
 */
#define ObitSDMDataRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSDMDataIsA(in) ObitIsA (in, ObitSDMDataGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSDMDataClassInit (void);

/** Public: Default Constructor. */
ObitSDMData* newObitSDMData (gchar* name);

/** Public: Create/initialize ObitSDMData structures */
ObitSDMData* ObitSDMDataCreate (gchar* name, gchar *DataRoot, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitSDMData* (*ObitSDMDataCreateFP) (gchar* name, gchar *DataRoot, 
					     ObitErr *err);

/** Public: Create/initialize ObitSDMData with only Intent structures */
ObitSDMData* ObitSDMIntentCreate (gchar* name, gchar *DataRoot, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitSDMData* (*ObitSDMIntentCreateFP) (gchar* name, gchar *DataRoot, 
					       ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitSDMDataGetClass (void);

/** Public: Copy (deep) constructor. */
ObitSDMData* ObitSDMDataCopy  (ObitSDMData *in, ObitSDMData *out, ObitErr *err);

/** Public: Copy structure. */
void ObitSDMDataClone (ObitSDMData *in, ObitSDMData *out, ObitErr *err);

/** Public: Get spectral window array. */
ASDMSpectralWindowArray* ObitSDMDataGetSWArray (ObitSDMData *in, olong mainRow, 
						gboolean SWOrder);

/** Public: Delete spectral window array. */
ASDMSpectralWindowArray* ObitSDMDataKillSWArray (ASDMSpectralWindowArray *in);

/** Public: Select Spectral window by number of channels  */
gboolean ObitSDMDataSelChan  (ASDMSpectralWindowArray *in, olong selChan, 
			      olong selIF, ObitASDMBand band);

/** Public: Get antenna/station array. */
ASDMAntennaArray* ObitSDMDataGetAntArray (ObitSDMData *in, olong mainRow);

/** Public: Delete antenna/station array. */
ASDMAntennaArray* ObitSDMDataKillAntArray (ASDMAntennaArray *in);

/** Public: Get Source/field array */
ASDMSourceArray* ObitSDMDataGetSourceArray (ObitSDMData *in);

/** Public: Delete Source/field  array. */
ASDMSourceArray* ObitSDMDataKillSourceArray (ASDMSourceArray *in);

/** Public: Convert band code string to band code enum */
ObitASDMBand ObitSDMDataBand2Band (gchar *code);

/** Public: Convert frequency to band code enum */
ObitASDMBand ObitSDMDataFreq2Band (odouble freq);

/** Find first Main table row matching selection */
olong ObitASDSelScan(ObitSDMData *in, olong selChan, olong selIF, 
		     ObitASDMBand band, olong selConfig);

/** Fix source numbers in an ASDMSourceTable */
void ObitSDMSourceTabFix(ObitSDMData *in);

/** Fix source numbers in an ASDMSourceTable including the calcode */
void ObitSDMSourceTabFixCode(ObitSDMData *in);

/** Public: Select Scan by code  */
gboolean ObitSDMDataSelCode  (ObitSDMData *in, olong iMain, gchar *selCode);

/** Find a pointingTab row for an antenna/time */
ASDMPointingRow* ObitSDMDataPointingLookup(ObitSDMData *in, odouble JD, olong ant, 
					   ObitErr *err);

/** Public: Convert ALMA WVR data to SN table  */
ObitTableSN* ObitSDMDataWVR2SN (ObitUV *inUV, ObitSDMData *SDM,
				ObitErr *err);

/** Public: Convert ALMA Atmosphere data to SN table  */
ObitTableSN* ObitSDMDataAtm2SN (ObitUV *inUV, ObitSDMData *SDM,
				ASDMSpectralWindowArray* SpWinArray, 
				ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitSDMDataClassDef.h"
} ObitSDMDataClassInfo; 

#endif /* OBITFSDMDATA_H */ 
