/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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
#ifndef OBITBDFDATA_H 
#define OBITBDFDATA_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFile.h"
#include "ObitUVDesc.h"
#include "ObitSDMData.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitBDFData.h
 *
 * This class accesses data in the EVLA BDF format
 *
 * Class documentation should go here.
 * 
 * \section ObitBDFDataaccess Creators and Destructors
 * An ObitBDFData will usually be created using ObitBDFDataCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitBDFData should always be made using the
 * #ObitBDFDataRef function which updates the reference count in the object.
 * Then whenever freeing an ObitBDFData or changing a pointer, the function
 * #ObitBDFDataUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitBDFMIMEType
 * Enum for xml mime type
 */
enum obitBDFMIMEType {
  BDFMIMEType_Unknown,
  BDFMIMEType_sdmDataHeader,       /* scan header */
  BDFMIMEType_desc,                /* integration header */
  BDFMIMEType_crossData,           /* cross correlation data */
  BDFMIMEType_autoData,            /* auto correlation data */
  BDFMIMEType_flags,               /* flag data */
  BDFMIMEType_actualTimes,         /* actual time data */
  BDFMIMEType_actualDurations,     /* actual duration data */
  BDFMIMEType_weights,             /* weight data */
  BDFMIMEType_zeroLags,            /* ??? */
  BDFMIMEType_EOF                  /* End of file */
}; /* end enum obitBDFMIMEType */
/** typedef for enum for BDFMIMEType. */
typedef enum obitBDFMIMEType ObitBDFMIMEType;

/**
 * \enum obitBDFBasebandName
 * Enum for Accumulation mode
 */
enum obitBDFBasebandName {
  BDFBasebandName_NOBB,
  BDFBasebandName_BB_1,
  BDFBasebandName_BB_2,
  BDFBasebandName_BB_3,
  BDFBasebandName_BB_4,
  BDFBasebandName_BB_5,
  BDFBasebandName_BB_6,
  BDFBasebandName_BB_7,
  BDFBasebandName_BB_8,
  BDFBasebandName_BB_ALL
}; /* end enum obitBDFBasebandName */
/** typedef for enum for BDFBasebandName. */
typedef enum obitBDFBasebandName ObitBDFBasebandName;

/**
 * \enum obitBDFCorrMode
 * Enum for Correlation Mode
 */
enum obitBDFCorrMode {
  BDFCorrMode_CROSS_ONLY,
  BDFCorrMode_AUTO_ONLY,
  BDFCorrMode_CROSS_AND_AUTO
}; /* end enum obitBDFCorrMode */
/** typedef for enum for BDFCorrMode. */
typedef enum obitBDFCorrMode ObitBDFCorrMode;

/**
 * \enum obitBDFAxisName
 * Enum for Axis name
 */
enum obitBDFAxisName {
  BDFAxisName_UNK,   /* Unknown */
  BDFAxisName_TIM,   /* Time axis */
  BDFAxisName_BAL,   /* Baseline axis */
  BDFAxisName_ANT,   /* Antenna axis */
  BDFAxisName_BAB,   /* Baseband axis */
  BDFAxisName_SPW,   /* Spectral window axis */
  BDFAxisName_SIB,   /* Sideband axis */
  BDFAxisName_SUB,   /* Subband axis */
  BDFAxisName_BIN,   /* binning axis */
  BDFAxisName_APC,   /* Atmospheric correction axis */
  BDFAxisName_SPP,   /* Spectral channel axis */
  BDFAxisName_POL,   /* Polarization axis */
  BDFAxisName_STO,   /* undocumented axis (possibly a synonym for POL) */
  BDFAxisName_HOL,   /* undocumented axis */
  BDFAxisName_END    /* last element in list */
}; /* end enum obitBDFAxisName */
/** typedef for enum for BDFAxisName. */
typedef enum obitBDFAxisName ObitBDFAxisName;

/**
 * \enum obitBDFSpecRes
 * Enum for spectral resolution type
 */
enum obitBDFSpecRes {
  BDFSpecRes_CHANNEL_AVERAGE,
  BDFSpecRes_BASEBAND_WIDE,
  BDFSpecRes_FULL_RESOLUTION
}; /* end enum obitBDFSpecRes */
/** typedef for enum for BDFSpecRes. */
typedef enum obitBDFSpecRes ObitBDFSpecRes;

/**
 * \enum obitBDFPolnType
 * Enum for Polarization type
 */
enum obitBDFPolnType {
  BDFPolnType_R,
  BDFPolnType_L,
  BDFPolnType_X,
  BDFPolnType_Y
}; /* end enum obitBDFPolnType */
/** typedef for enum for BDFPolnType. */
typedef enum obitBDFPolnType ObitBDFPolnType;

/**
 * \enum obitBDFDataType
 * Enum for data type
 */
enum obitBDFDataType {
  BDFDataType_Unknown,
  BDFDataType_INT16_TYPE,
  BDFDataType_INT32_TYPE,
  BDFDataType_INT64_TYPE,
  BDFDataType_FLOAT32_TYPE,
  BDFDataType_FLOAT64_TYPE
}; /* end enum obitBDFDataType */
/** typedef for enum for BDFDataType. */
typedef enum obitBDFDataType ObitBDFDataType;

/**
 * \enum obitBDFEndian
 * Enum for endianness of data
 */
enum obitBDFEndian {
  BDFEndian_Big,
  BDFEndian_Little
}; /* end enum obitBDFEndian */
/** typedef for enum for BDFEndian. */
typedef enum obitBDFEndian ObitBDFEndian;

/*----------------- Macroes ---------------------------*/
/** Granularity of buffer operations (frame size) */
#define BDFBUFFERSIZE 16384
/** Number of frames in buffer */
#define BDFBUFFERFRAMES 32
/** Maximum number of BasebandInfo structs */
#define MAXBBINFO 128
/** Maximum number of Spectral windows per BasebandInfo */
#define MAXBBSW 32
/** 
 * Macro to unreference (and possibly destroy) an ObitBDFData
 * returns a ObitBDFData*.
 * in = object to unreference
 */
#define ObitBDFDataUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitBDFData.
 * returns a ObitBDFData*.
 * in = object to reference
 */
#define ObitBDFDataRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitBDFDataIsA(in) ObitIsA (in, ObitBDFDataGetClass())

/*--------------Class definitions-------------------------------------*/
/** ObitBDFData BDF/BDF structures. */
#include "ObitBDFDataTypeDef.h" 

/** ObitBDFData Class structures. */
typedef struct {
#include "ObitBDFDataDef.h"   /* this class definition */
} ObitBDFData;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitBDFDataClassInit (void);

/** Public: Default Constructor. */
ObitBDFData* newObitBDFData (gchar* name);

/** Public: Create/initialize ObitBDFData structures */
ObitBDFData* ObitBDFDataCreate (gchar* name, ObitUVDesc *desc, 
				ObitSDMData *SDMData, 
				ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitBDFData* (*ObitBDFDataCreateFP) (gchar* name, ObitUVDesc *desc, 
					     ObitSDMData *SDMData, 
					     ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitBDFDataGetClass (void);

/** Public: Copy (deep) constructor. */
ObitBDFData* ObitBDFDataCopy  (ObitBDFData *in, ObitBDFData *out, ObitErr *err);

/** Public: Copy structure. */
void ObitBDFDataClone (ObitBDFData *in, ObitBDFData *out, ObitErr *err);

/** Public: Initialize file */
void ObitBDFDataInitFile (ObitBDFData *in, gchar *DataFile, ObitErr *err);

/** Public: Fill Buffer */
ObitIOCode ObitBDFDataFillBuffer (ObitBDFData *in, ObitErr *err);

/** Public: Initialize Scan */
void ObitBDFDataInitScan (ObitBDFData *in, olong iMain, gboolean SWOrder,
			  ObitErr *err);

/** Public: Select Spectral window by number of channels  */
gboolean ObitBDFDataSelChan  (ObitBDFData *in, olong selChan, 
			      olong selIF, ObitASDMBand band);

/** Public: Initialize Integration */
ObitIOCode ObitBDFDataInitInteg (ObitBDFData *in, ObitErr *err);

/** Public: Read Next Integration */
ObitIOCode ObitBDFDataReadInteg (ObitBDFData *in, ObitErr *err);

/** Public: Get visibility record */
ObitIOCode ObitBDFDataGetVis (ObitBDFData *in, ofloat *vis, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitBDFDataClassDef.h"
} ObitBDFDataClassInfo; 

#endif /* OBITFBDFDATA_H */ 
