/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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
/*  Define the basic components of the ObitBDFData structure         */
/*  This is intended to be included in a class structure definition  */
/* This class accesses data in the EVLA BDF format                   */
/* This file defines BDF structures                          */

/** Spectral Window info from XML           */
typedef struct {
 /** spectral window number (NOT Id) */
  olong   spectralWindowNum;
 /** List of single dish (autocorrelation) products, e.g. RR, RL, LL */
  gchar   **sdPolProducts;
 /** Number of autocorrelation products */
  olong   numSdPolProducts;
 /**  List of crosscorrelation products, e.g. RR, RL, LR, LL */
  gchar   **crossPolProducts;
 /** Number of cross correlation products */
  olong   numCrossPolProducts;
 /** Number of spectral points */
  olong   numSpectralPoint;
 /** Number of data (e.g. pulsar) bins */
  olong   numBin;
  /* Scaling factor */
  odouble scaleFactor;
  /** Sideband */
  gchar *sideband;
} BDFSpecWindowInfo;

/** Baseband info from XML           */
typedef struct {
 /** Baseband name */
  gchar   *basebandName;
  /** Number of spectral windows */
  olong numSpectralWindow;
  /** Array of Spectral window information dimension(MAXBBSW) */
  BDFSpecWindowInfo** SWinds;
} BDFBasebandInfo;

/**
 * \file ObitBDFDataTypeDef.h
 * ObitBDFData BDF structures
 */
 /** Scan info from XML           */
typedef struct {
  /** Start time JD */
  odouble startTime;
  /** ASDM Main table row index */
  olong iMain;
  /** Number of antennas */
  olong numAntenna;
  /** Number of Times */
  olong numTime;
  /** correlation Mode */
  ObitBDFCorrMode correlationMode;
  /** spectral Resolution enum */
  ObitBDFSpecRes spectralResolution;
  /** Endian */
  ObitBDFEndian endian;
  /** Number of basebands */
  olong numBaseband;
  /** Array of baseband info structures */
  BDFBasebandInfo *BBinfo[MAXBBINFO];
  /** Flag array size */
  olong FlagSize;
  /* Order of flag array axes */
  ObitBDFAxisName *FlagAxes;
  /** actual times array size */
  olong actualTimesSize;
  /** Order of actual times axes */
  ObitBDFAxisName *actualTimesAxes;
  /**  actualDurations array size */
  olong actualDurationsSize;
  /** Order of actualDurations axes */
  ObitBDFAxisName *actualDurationsAxes;
  /**  crossData array size */
  olong crossDataSize;
  /** Order of crossData axes */
  ObitBDFAxisName *crossDataAxes;
  /**  autoData array size */
  olong autoDataSize;
  /** Order of autoData axes */
  ObitBDFAxisName *autoDataAxes;
  /**  weight array size */
  olong weightSize;
  /** Order of weight axes */
  ObitBDFAxisName *weightAxes;
  /**  zeroLag array size */
  olong zeroLagSize;
  /** Order of zeroLag axes */
  ObitBDFAxisName *zeroLagAxes;
} BDFScanInfo;

/** Integration info from XML           */
typedef struct {
  /** Time JD */
  odouble time;
  /** integation time days */
  odouble interval;
  /** Data type */
  ObitBDFDataType type;
} BDFIntegInfo;
