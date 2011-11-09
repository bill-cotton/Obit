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
/* This file defines ASDM and BDF structures                          */
/**
 * \file ObitSDMDataTypeDef.h
 * ObitSDMData ASDM structures
 */

/* Structures interpreted from the ASDM */
/* SpectralWindow/Poln array */
typedef struct {
  /** Spectral window ID */
  olong spectralWindowId;
  /** Selected? */
  gboolean selected;
  /** Baseband number */
  olong basebandNum;
  /** subband number */
  olong subbandNum;
  /** Number of channels */
  olong numChan;
  /** Number of cross poln */
  olong nCPoln;
  /** Number of auto poln */
  olong nAPoln;
  /** net Sideband */
  gchar *netSideband;
  /** refFreq */
  odouble refFreq;
  /** total Bandwidth (Hz) */
  odouble totBandwidth;
  /** channel Frequency Start (Hz) */
  odouble chanFreqStart;
  /** channel Freq Step (Hz) */
  odouble chanFreqStep;
  /** channel Width (Hz) */
  odouble chanWidth;
  /** effectiveBw */
  odouble effectiveBw;
  /** spectral resolution in Hz */
  odouble resolution;
} ASDMSpectralWindowArrayEntry;
typedef struct {
  /* Reference Frequency */
  odouble refFreq;
  /* Reference JD */
  odouble refJD;
  /** Number of rows */
  olong nwinds;
  /** 0-rel order number of Spectral windows in data */
  olong *order;
  /** Array of ASDMSpectralWindowArray rows */
  ASDMSpectralWindowArrayEntry **winds;
} ASDMSpectralWindowArray;

/* Antenna/station array */
typedef struct {
  /** Antenna ID */
  olong antennaId;
  /** Antenna Number */
  olong antennaNo;
  /** name */
  gchar* antName;
  /** antenna Make enum */
  ObitASDMAntennaMake antennaMake;
  /** antenna Type enum */
  ObitASDMAntennaType antennaType;
  /** dish Diameter (m)*/
  odouble dishDiameter;
  /** position array of doubles */
  odouble *antPosition;
  /** offset array of doubles */
  odouble *offset;
  /** time */
  odouble time ;
  /** station Id */
  olong stationId;
  /** Station name */
  gchar *staName;
  /** position */
  odouble *staPosition;
  /** type */
  ObitASDMStationType type;
  /** Number of polarizations (1 or 2) */
  olong numPoln;
  /** Number of polarization orrelations (1-4) */
  olong numPolnCorr;
  /** Polarization types 'R', 'L', 'X', 'Y' */
  gchar polnType[2];
} ASDMAntennaArrayEntry;
typedef struct {
  /** Array Name */
  gchar *arrayName;
  /** Observer Name */
  gchar *obsName;
  /* Reference Frequency */
  odouble refFreq;
  /* Reference JD */
  odouble refJD;
  /** Number of rows */
  olong nants;
  /** Highest antenna number */
  olong maxAnt;
  /** Array of Antenna entries */
  ASDMAntennaArrayEntry **ants;
} ASDMAntennaArray;

/** Source array - one entry per Spectral Window */
typedef struct {
  /** source Id */
  olong sourceId;
  /** source Number (1-rel) */
  olong sourceNo;
  /** time Interval as JD */
  odouble *timeInterval;
  /** cal code */
  gchar *code;
  /** direction, array of RA, Dec in rad */
  odouble *direction;
  /** proper Motion, array of RA, Dec in rad */
  odouble *properMotion;
  /** sourceName */
  gchar *sourceName;
  /** numLines */
  olong numLines;
  /** rest Frequency, one per line */
  odouble *restFrequency;
  /** systemic Velocity, one per line */
  odouble *sysVel;
  /** spectral Window Id */
  olong spectralWindowId;
} ASDMSourceArrayEntry;
typedef struct {
  /** Number of rows */
  olong nsou;
  /** Array of Source entries */
  ASDMSourceArrayEntry **sou;
} ASDMSourceArray;

 /** EVLA ASDM Tables               */
 /* ASDM Table */
typedef struct {
  /** Schema version number */
  olong schemaVersion;
  /** Number of Main rows */
  olong MainRows;
  /** Number of Antenna rows */
  olong AntennaRows;
  /** Number of calAtmosphere rows */
  olong calAtmosphereRows;
  /** Number of calData rows */
  olong calDataRows;
  /** Number of calDevice rows */
  olong calDeviceRows;
  /** Number of calPointing rows */
  olong calPointingRows;
  /** Number of CalReduction rows */
  olong CalReductionRows;
  /** Number of CalWVR rows */
  olong CalWVRRows;
  /** Number of ConfigDescription rows */
  olong ConfigDescriptionRows;
  /** Number of CorrelatorMode rows */
  olong CorrelatorModeRows;
  /** Number of DataDescription rows */
  olong DataDescriptionRows;
  /** Number of Doppler rows */
  olong DopplerRows;
  /** Number of ExecBlock rows */
  olong ExecBlockRows;
  /** Number of Feed rows */
  olong FeedRows;
  /** Number of Field rows */
  olong FieldRows;
  /** Number of Flag rows */
  olong FlagRows;
   /** Number of Pointing rows */
  olong PointingRows;
  /** Number of PointingModel rows */
  olong PointingModelRows;
  /** Number of Polarization rows */
  olong PolarizationRows;
  /** Number of Processor rows */
  olong ProcessorRows;
  /** Number of Receiver rows */
  olong ReceiverRows;
  /** Number of SBSummary rows */
  olong SBSummaryRows;
  /** Number of Scan rows */
  olong ScanRows;
  /** Number of Source rows */
  olong SourceRows;
  /** Number of SpectralWindow rows */
  olong SpectralWindowRows;
  /** Number of State rows */
  olong StateRows;
  /** Number of Station rows */
  olong StationRows;
  /** Number of Subscan rows */
  olong SubscanRows;
  /** Number of SwitchCycle rows */
  olong SwitchCycleRows;
  /** Number of SysCal rows */
  olong SysCalRows;
  /** Number of SysPower rows */
  olong SysPowerRows;
  /** Number of Weather rows */
  olong WeatherRows;
} ASDMTable;


 /* Main Table */
typedef struct {
  odouble time;
  /** Number of antennas */
  olong numAntenna;
  /** Type of time sampling */
  ObitSDMSampleType timeSampling;
  /** Time interval in days covered */
  odouble interval;
  /** Number of integrations */
  olong numIntegration;
  /** Scan number */
  olong scanNumber;
  /** Subscan number */
  olong subscanNumber;
  /** Number of bytes in the BDF data file 64 bits */
  ollong dataSize;
  /** Name of associated BDF data file (replaced slash with underscore) */
  gchar *entityId;
  /** Is all data with this row flagged? */
  gboolean flagRow;
  /** Configuration description Identifier */
  olong configDescriptionId;
  /** execution block identifier */
  olong execBlockId;
  /** field identifier */
  olong fieldId;
  /** State identifier */
  olong *stateId;
} ASDMMainRow;

typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMMain rows */
  ASDMMainRow **rows;
} ASDMMainTable;


 /* Antenna Table */
typedef struct {
  /** Antenna ID */
  olong antennaId;
  /** name */
  gchar* name;
  /** antenna Make enum */
  ObitASDMAntennaMake antennaMake;
  /** antenna Type enum */
  ObitASDMAntennaType antennaType;
  /** dish Diameter (m)*/
  odouble dishDiameter;
  /** position array of doubles */
  odouble *position;
  /** offset array of doubles */
  odouble *offset;
  /** time */
  odouble time ;
  /** station Id */
  olong stationId;
} ASDMAntennaRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMAntenna rows */
  ASDMAntennaRow **rows;
} ASDMAntennaTable;

 /* ALMA CalAtmosphere Table, ignore spectra for now */
typedef struct {
  /** receiverBand */
  olong receiverBand;
  /** antenna name */
  gchar *antennaName;
  /** syscalType */
  gchar *syscalType;
  /** baseband Id */
  olong basebandId;
  /** number of frequencies */
  olong numFreq;
  /** number of Loads */
  olong numLoad;
  /** number of Receptors */
  olong numReceptor;
  /** calData Id */
  olong calDataId ;
  /**  cal reduction Id */
  olong calReductionId ;
  /** polarization types (1D array of poln enums) */
  olong *polarizationTypes;
  /** start time (days) */
  odouble startValidTime;
  /** end time (days) */
  odouble endValidTime;
  /** Ground pressure */
  odouble groundPressure;
  /** Ground relative humidity */
  odouble groundRelHumidity;
  /** Ground temperature (K) */
  odouble groundTemperature;
  /** Frequency range (Hz) */
  odouble *frequencyRange;
  /** Atm. temp per receptor */
  odouble *tAtm;
  /** Receiver temp. per receptor */
  odouble *tRec ;
  /** System  temp. at top of atmosphere per receptor */
  odouble *tSys;
  /** Opacity per receptor */
  odouble *tau;
  /** Water per receptor */
  odouble *water;
  /** Error in water, per receptor */
  odouble *waterError;
  /** forward efficiency, per receptor */
  odouble *forwardEfficiency;
  /** sb(?) gain, per receptor  */
  odouble *sbGain;
} ASDMcalAtmosphereRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMCalAtmosphere rows */
  ASDMcalAtmosphereRow **rows;
} ASDMcalAtmosphereTable;

 /* calData Table */
typedef struct {
  /** place holder */
  olong holder;
} ASDMcalDataRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMcalData rows */
  ASDMcalDataRow **rows;
} ASDMcalDataTable;

 /* calDevice Table */
typedef struct {
  /** antenna Id */
  olong antennaId;
  /** spectral Window Id */
  olong spectralWindowId;
  /** feedId */
  olong feedId;
  /** number of Cal loads (N_cal) */
  olong numCalLoad;
  /** number of Receptors (N_rec) */
  olong numReceptor;
  /** Names of calibration devices (1D array of strings per cal device) */
  gchar **calLoadNames;
  /** time Interval (days) */
  odouble *timeInterval;
  /**  Calibration efficiencies (2D array of double [N_rec][N_cal]) */
  odouble *calEff;
  /**  equivalent temp of noise cal (1D array of double, one per load) obsolete */
  odouble *noiseCal;
  /**  equivalent temp of noise cal (1D array of double, one per poln (R, L)) */
  odouble *coupledNoiseCal;
  /**  Physical temperature of loads (1D array of double, one per load) */
  odouble *temperatureLoad;
} ASDMcalDeviceRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMcalData rows */
  ASDMcalDeviceRow **rows;
} ASDMcalDeviceTable;

 /* calPointing Table */
typedef struct {
  /** place holder */
  olong holder;
} ASDMcalPointingRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMcalPointing rows */
  ASDMcalPointingRow **rows;
} ASDMcalPointingTable;

 /* CalReduction Table */
typedef struct {
  /** place holder */
  olong holder;
} ASDMCalReductionRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMCalReduction rows */
  ASDMCalReductionRow **rows;
} ASDMCalReductionTable;

 /* CalWVR Table */
typedef struct {
  /** start Time JD*/
  odouble startValidTime;
  /** end Time JD */
  odouble endValidTime;
  /**  wvr method */
  gchar *wvrMethod;
  /**  antenna Name */
  gchar *antennaName;
  /** numChan */
  olong numChan;
  /** Arrays of info - dim numChan
      channel frequency array */
  odouble *chanFreq;
  /* channel width array */
  odouble *chanWidth;
  /* channel temperatures */
  odouble *refTemp;
  /** number of Polynomials */
  olong numPoly;
  /* path coefficients */
  odouble *pathCoeff;
  /* polyFreqLimits */
  odouble *polyFreqLimits;
  /* wet Path */
  odouble wetPath;
  /* dry Path */
  odouble dryPath;
  /* water */
  odouble water;
} ASDMCalWVRRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMCalWVR rows */
  ASDMCalWVRRow **rows;
} ASDMCalWVRTable;

 /* ConfigDescription Table */
typedef struct {
  /** Number of antennas */
  olong numAntenna;
  /** Number of Data Descriptions */
  olong numDataDescription;
  /** Number of feeds */
  olong numFeed;
  /** correlation Mode */
  ObitASDMCorrMode correlationMode;
  /** config Description Id */
  olong configDescriptionId;
  /**  number of Atm Phase Correction */
  olong numAtmPhaseCorrection;
  /** Atm Phase Correction type enum array */
  ObitASDMAtmPhCorr *atmPhaseCorrection;
  /** processor Type enum */
  ObitASDMProcrType processorType;
  /** spectral Type enum */
  ObitASDMSpecRes spectralType;
  /** antenna Ids in configuration array */
  olong *antennaId;
  /** data Description Id array */
  olong *dataDescriptionId;
  /** feed Id array */
  olong *feedId;
  /** processor Id */
  olong processorId;
  /** switch Cycle Id */
  olong *switchCycleId;
} ASDMConfigDescriptionRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMConfigDescription rows */
  ASDMConfigDescriptionRow **rows;
} ASDMConfigDescriptionTable;

 /* CorrelatorMode Table */
typedef struct {
  /** correlator Mode Id */
  olong correlatorModeId;
  /** number of Basebands */
  olong numBaseband;
  /** baseband Name enums, per baseband */
  ObitASDMBasebandName *basebandNames;
  /** baseband Config array, per baseband*/
  olong *basebandConfig;
  /** accumMode enum*/
  ObitASDMAccumMode accumMode;
  /** bin Mode */
  olong binMode;
  /** number of Axes */
  olong numAxes;
  /** axes Order Array */
  ObitASDMAxisName *axesOrderArray;
  /** filterMode enum array per baseband */
  ObitASDMFilterMode *filterMode;
  /** correlator Name */
  gchar *correlatorName;
} ASDMCorrelatorModeRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMCorrelatorMode rows */
  ASDMCorrelatorModeRow **rows;
} ASDMCorrelatorModeTable;

 /* DataDescription Table */
typedef struct {
  /** data Description Id */
  olong dataDescriptionId;
  /** polarization Id */
  olong polOrHoloId;
  /** spectra lWindow Id */
  olong spectralWindowId;
} ASDMDataDescriptionRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMDataDescription rows */
  ASDMDataDescriptionRow **rows;
} ASDMDataDescriptionTable;

 /* Doppler Table */
typedef struct {
  /** place holder */
  olong holder;
} ASDMDopplerRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMDoppler rows */
  ASDMDopplerRow **rows;
} ASDMDopplerTable;

 /* ExecBlock Table */
typedef struct {
  /** exec Block Id */
  olong execBlockId;
  /** start Time (JD) */
  odouble startTime;
  /** end Time  (JD) */
  odouble endTime ;
  /** exec Block Number */
  olong execBlockNum;
  /** array config Name */
  gchar *configName;
  /** telescope Name */
  gchar *telescopeName;
  /** observer Name */
  gchar *observerName ;
  /** observing Log */
  gchar **observingLog ;
  /** session Reference */
  gchar *sessionReference;
  /** scheduler Mode */
  gchar *schedulerMode ;
  /** Number of observing Log */
  olong numObservingLog;
  /** baseline Range Min */
  odouble baseRangeMin;
  /** baseline Range Max */
  odouble baseRangeMax;
  /** baseline Rms Minor */
  odouble baseRmsMinor;
  /** basline Rms Major */
  odouble baseRmsMajor;
  /** baseline Pa */
  odouble basePa;
  /** site Altitude (m) */
  odouble siteAltitude;
  /** site Longitude (deg) */
  odouble siteLongitude;
  /** site Latitude (deg) */
  odouble siteLatitude;
  /** was scan aborted */
  gboolean aborted;
  /** number of Antenna */
  olong numAntenna;
  /** flag this Row */
  gboolean flagRow;
  /** antenna Id array */
  olong *antennaId;
  /** sbSummary Table ID */
  olong sbSummaryId;
  /** scale Id  */
  olong ScaleId;
} ASDMExecBlockRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMExecBlock rows */
  ASDMExecBlockRow **rows;
} ASDMExecBlockTable;

 /* Feed Table */
typedef struct {
  /** feedId */
  olong feedId;
  /** time Interval (days) */
  odouble *timeInterval;
  /** number of Receptors */
  olong numReceptor;
  /** beam Offset (2D array of double) */
  odouble *beamOffset;
  /** focus Reference  (2D array of double) */
  odouble *focusReference;
  /** polarizationTypes (1D array of poln enums) */
  olong *polarizationTypes;
  /** pol Response (2D array of double)*/
  odouble *polResponse;
  /** receptor Angle (1D array of double) */
  odouble *receptorAngle;
  /** antenna Id */
  olong antennaId;
  /** receiver Id (1D array of int) */
  olong *receiverId;
  /** spectral Window Id */
  olong spectralWindowId;
} ASDMFeedRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMFeed rows */
  ASDMFeedRow **rows;
} ASDMFeedTable;

 /* Field Table */
typedef struct {
  /** field Id */
  olong fieldId;
  /** field Name */
  gchar *fieldName;
  /** code */
  gchar *code;
  /** number of Polynomials */
  olong numPoly;
  /** delay  tracking center (rad) possibly also derivatives */
  odouble* delayDir;
  /**  phase tracking center (rad) possibly also derivatives */
  odouble *phaseDir;
  /** reference direction (rad) possibly also derivatives */
  odouble *referenceDir;
  /** time (JD)*/
  odouble time;
  /** source Id */
  olong sourceId;
} ASDMFieldRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMField rows */
  ASDMFieldRow **rows;
} ASDMFieldTable;

 /* Flag Table */
typedef struct {
  /** flag Id */
  olong flagId;
  /** antenna Id array */
  olong *antennaId;
  /** number of antennas */
  olong numAntenna;
  /** number of polarization types, 0=> all */
  olong numPolarizationType;
  /** number of spectral windows, 0=> all */
  olong numSpectralWindow;
  /** reason */
  gchar *reason;
  /** start Time (JD) */
  odouble startTime;
  /** end Time (JD) */
  odouble endTime;
} ASDMFlagRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMFlag rows */
  ASDMFlagRow **rows;
} ASDMFlagTable;

 /* Pointing Table */
typedef struct {
  /** time Interval (JD days) */
  odouble *timeInterval;
  /** number of samples */
  olong numSample;
  /** Encoder values */
  odouble *encoder;
  /** Was antenna in tracking mode */
  gboolean pointingTracking;
  /** using Polynomials */
  gboolean usePolynomials;
  /** Time origin JD of the polynomial expansion */
  odouble timeOrigin;
  /** number of terms in the polynomial */
  olong numTerm;
  /** Commanded pointing direction */
  odouble *pointingDirection;
  /** Direction of target */
  odouble *target;
  /** Horizion mapping offsets */
  odouble *offset;
  /** over the top observing? */
  gboolean overTheTop;
  /** antenna Id identifying row in Antenna table */
  olong antennaId;
  /** pointing model Id */
  olong pointingModelId;
} ASDMPointingRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMPointing rows */
  ASDMPointingRow **rows;
} ASDMPointingTable;

 /* PointingModel Table */
typedef struct {
  /** place holder */
  olong holder;
} ASDMPointingModelRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMPointingModel rows */
  ASDMPointingModelRow **rows;
} ASDMPointingModelTable;

 /* Polarization Table */
typedef struct {
  /** polarization Id */
  olong polarizationId;
  /** number of Correlations */
  olong numCorr;
  /** correlation polarization Type (1D array of numCorr poln codes) */
  gchar **corrType;
  /** corrProduct (2D array of numCorr x 2 poln codes, R, L, X, Y) */
  gchar **corrProduct;
  /** flagRow */
  gboolean flagRow;
} ASDMPolarizationRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMPolarization rows */
  ASDMPolarizationRow **rows;
} ASDMPolarizationTable;

 /* Processor Table */
typedef struct {
  /** processor Id */
  olong processorId;
  /** modeId */
  olong modeId;
  /** processorType */
  gchar *processorType;
  /** processorSubType */
  gchar *processorSubType;
} ASDMProcessorRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMProcessor rows */
  ASDMProcessorRow **rows;
} ASDMProcessorTable;

 /* Receiver Table */
typedef struct {
  /** place holder */
  olong holder;
} ASDMReceiverRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMReceiver rows */
  ASDMReceiverRow **rows;
} ASDMReceiverTable;

 /* SBSummary Table */
typedef struct {
  /** place holder */
  olong holder;
} ASDMSBSummaryRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMSBSummary rows */
  ASDMSBSummaryRow **rows;
} ASDMSBSummaryTable;

 /* Scan Table */
typedef struct {
  /** scan Number */
  olong scanNumber;
  /** start Time JD*/
  odouble startTime;
  /** end Time JD */
  odouble endTime;
  /** numIntent */
  olong numIntent;
  /** numSubscan */
  olong numSubscan;
  /** scan Intent array of strings */
  gchar **scanIntent;
  /** cal Data Type string array, one per subscan[?] */
  gchar **calDataType;
  /** calibration On Line? boolean array  */
  gboolean *calibrationOnLine;
  /** source Name */
  gchar *sourceName;
  /** Row flagged */
  gboolean flagRow;
  /** execBlockId */
  olong execBlockId;
} ASDMScanRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMScan rows */
  ASDMScanRow **rows;
} ASDMScanTable;

 /* Source Table */
typedef struct {
  /** source Id */
  olong sourceId;
  /** source Number (1-rel) */
  olong sourceNo;
  /** time Interval as JD */
  odouble *timeInterval;
  /** code */
  gchar *code;
  /** direction, array of RA, Dec in rad */
  odouble *direction;
  /** proper Motion, array of RA, Dec in rad */
  odouble *properMotion;
  /** sourceName */
  gchar *sourceName;
  /** numLines */
  olong numLines;
  /** rest Frequency, one per line */
  odouble *restFrequency;
  /** systemic Velocity, one per line */
  odouble *sysVel;
  /** spectral Window Id */
  olong spectralWindowId;
} ASDMSourceRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMSource rows */
  ASDMSourceRow **rows;
} ASDMSourceTable;

 /* SpectralWindow Table */
typedef struct {
  /** spectralWindowId */
  olong spectralWindowId;
  /**  basebandName */
  gchar *basebandName;
  /** net Sideband */
  gchar *netSideband;
  /** numChan */
  olong numChan;
  /** refFreq */
  odouble refFreq;
  /** sideband Processing Mode enum */
  ObitASDMSideBMode sidebandProcessingMode;
  /** total Bandwidth (Hz) */
  odouble totBandwidth;
  /** window Function enumeration */
  ObitASDMWindowFn windowFunction;
  /** channel Frequency Start (Hz) */
  odouble chanFreqStart;
  /** channel Freq Step (Hz) */
  odouble chanFreqStep;
  /** channel Width (Hz) */
  odouble chanWidth;
  /** correlationBit */
  gchar *correlationBit;
  /** effectiveBw */
  odouble effectiveBw;
  /** name */
  gchar *name;
  /** oversampling?  */
  gboolean oversampling ;
  /** quantization[?]  */
  gboolean quantization;
  /** spectral resolution in Hz */
  odouble resolution;
  /** Arrays of info - dim numChan
      channel frequency array */
  odouble *chanFreqArray;
  /* channel width array */
  odouble *chanWidthArray;
  /* effective bandwidth array */
  odouble *effectiveBwArray;
  /* resolution array */
  odouble *resolutionArray;
  /** Number of associated values */
  olong numAssocValues;
  /** Nature of associated values (ObitASDMSpecRes as olong) */
  olong *SpecRes;
  /** Associated Spectral window IDs */
  olong *assocSpectralWindowId;
} ASDMSpectralWindowRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMSpectralWindow rows */
  ASDMSpectralWindowRow **rows;
} ASDMSpectralWindowTable;

 /* State Table */
typedef struct {
  /** stateId */
  olong stateId;
  /** calibration Device Name */
  gchar *calDeviceName;
  /** Looking at sky? */
  gboolean sig;
  /** Looking at reference? */
  gboolean ref;
  /** On sky? */
  gboolean onSky;
} ASDMStateRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMState rows */
  ASDMStateRow **rows;
} ASDMStateTable;

 /* Station Table */
typedef struct {
  /** station Id */
  olong stationId;
  /** name */
  gchar *name;
  /** position */
  odouble *position;
  /** type */
  ObitASDMStationType type;
} ASDMStationRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMStation rows */
  ASDMStationRow **rows;
} ASDMStationTable;

 /* Subscan Table */
typedef struct {
  /** scanNumber */
  olong scanNumber;
  /** subscanNumber */
  olong subscanNumber;
  /** start Time JD */
  odouble startTime;
  /** end Time JD */
  odouble endTime;
  /** fieldName */
  gchar *fieldName;
  /** subscanIntent */
  gchar *subscanIntent;
  /** number of Integrations */
  olong numberIntegration;
  /** numberSubintegration 1 per integration */
  olong *numberSubintegration;
  /** flagRow */
  gboolean flagRow;
  /** execBlockId */
  olong execBlockId;
} ASDMSubscanRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMSubscan rows */
  ASDMSubscanRow **rows;
} ASDMSubscanTable;

 /* SwitchCycle Table */
typedef struct {
  /** switch Cycle Id */
  olong switchCycleId;
  /** number of steps Step */
  olong numStep;
  /** weight Array */
  odouble *weightArray;
  /** dir Offset Array */
  odouble *dirOffsetArray;
  /** freq Offset Array */
  odouble *freqOffsetArray;
  /** step Duration Array */
  odouble *stepDurationArray;
} ASDMSwitchCycleRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMSwitchCycle rows */
  ASDMSwitchCycleRow **rows;
} ASDMSwitchCycleTable;

 /* Weather Table */
typedef struct {
  /** time Interval (days) */
  odouble *timeInterval;
  /** atmospheric pressure (Pascal), fblank = invalid */
  ofloat  pressure;
  /** relative humidity, fblank = invalid */
  ofloat  relHumidity;
  /** temperature(K), fblank = invalid */
  ofloat temperature ;
  /** wind direction, Azimuth?, radians?, fblank = invalid */
  ofloat windDirection ;
  /** wind speed (m/s), fblank = invalid */
  ofloat  windSpeed;
  /** wind max - max gust(?) (m/s), fblank = invalid */
  ofloat  windMax;
  /** dew point (K), fblank = invalid */
  ofloat dewPoint;
  /** station Id */
  olong  stationId;
} ASDMWeatherRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMWeather rows */
  ASDMWeatherRow **rows;
} ASDMWeatherTable;

 /* SysCal Table */
typedef struct {
  /** place holder */
  olong holder;
} ASDMSysCalRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMSysCal rows */
  ASDMSysCalRow **rows;
} ASDMSysCalTable;

 /* SysPower Table */
typedef struct {
  /** antenna Id */
  olong antennaId;
  /** spectral Window Id */
  olong spectralWindowId;
  /** feedId */
  olong feedId;
  /** number of Receptors */
  olong numReceptor;
  /** time Interval (days) */
  odouble *timeInterval;
  /** switched Power Difference (1D array of double) */
  ofloat *switchedPowerDifference;
  /** switched Power Sum (1D array of double) */
  ofloat *switchedPowerSum;
  /** Requantizer gain Sum (1D array of double) */
  ofloat *requantizerGain;
} ASDMSysPowerRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMSysPower rows */
  ASDMSysPowerRow **rows;
} ASDMSysPowerTable;

 /* XXXX Table */
typedef struct {
  /** place holder */
  olong holder;
} ASDMXXXXRow;
typedef struct {
  /** Number of rows */
  olong nrows;
  /** Array of ASDMXXXX rows */
  ASDMXXXXRow **rows;
} ASDMXXXXTable;

