/* $Id$     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITTIMEFILTER_H 
#define OBITTIMEFILTER_H 

#include "ObitThread.h"
#include "ObitFArray.h"
#include "ObitCArray.h"
#include "ObitFFT.h"
#include "ObitFInterpolate.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTimeFilter.h
 * ObitTimeFilter time filter class definition.
 *
 * This class is derived from the #Obit class.
 * This class is for performing TimeFilter on memory resident data.
 *
 * \section ObitTimeFilter Creators and Destructors
 * An ObitTimeFilter can be created using newObitTimeFilter which allows specifying 
 * a name for the object, and the type, size and direction of the transform.
 *
 * A copy of a pointer to an ObitTimeFilter should always be made using the
 * #ObitTimeFilterRef function which updates the reference count in the object.
 * Then whenever freeing an ObitTimeFilter or changing a pointer, the function
 * #ObitTimeFilterUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 * \section ObitTimeFilter usage
 * An ObitTimeFilter provides storage for one or more time series 
 * of floats in both the time and frequency (half complex) domains.
 * The class also contains tools for transforming from one domain to 
 * the other and for applying a filter in the frequency domain.  
 * Typcal usage sequence is:
 * \li Create using  #newObitTimeFilter
 * \li fill time series data using #timeData ofloat pointer member
 * \li transform to frequency using #ObitTimeFilter2Freq
 * \li apply filter using #ObitTimeFilterFilter
 * \li transform back to time domain using #ObitTimeFilter2Time
 * \li access filter time series using #timeData ofloat pointer member
 * \li unreference object using #ObitTimeFilterUnref
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitTimeFilterType
 * enum for type of ObitTimeFilter filter type.
 * This specifies the type of filtering to be performed.
 */
enum obitTimeFilterType {
  /** Low pass filter */
  OBIT_TimeFilter_LowPass, 
  /** High pass filter */
  OBIT_TimeFilter_HighPass, 
  /** Notch pass filter */
  OBIT_TimeFilter_NotchPass, 
  /** Notch block filter */
  OBIT_TimeFilter_NotchBlock, 
}; /* end enum obitTimeFilterType */

/** typedef for enum for obitTimeFilterType. */
typedef enum obitTimeFilterType ObitTimeFilterType;

/*--------------Class definitions-------------------------------------*/
/** ObitTimeFilter Class structure. */
typedef struct {
#include "ObitTimeFilterDef.h"   /* this class definition */
} ObitTimeFilter;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitTimeFilter
 * returns a ObitTimeFilter*.
 * in = object to unreference
 */
#define ObitTimeFilterUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitTimeFilter.
 * returns a ObitTimeFilter*.
 * in = object to reference
 */
#define ObitTimeFilterRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitTimeFilterIsA(in) ObitIsA (in, ObitTimeFilterGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitTimeFilterClassInit (void);

/** Public: Constructor. */
ObitTimeFilter* newObitTimeFilter (gchar* name, olong nTime, olong nSeries);

/** Typedef for definition of class pointer structure */
typedef ObitTimeFilter* (*newObitTimeFilterFP) (gchar *name, 
						gint nTime, olong nFilter);

/** Public: ClassInfo pointer */
gconstpointer ObitTimeFilterGetClass (void);

/** Public: Resize arrays. */
void ObitTimeFilterResize (ObitTimeFilter *in, olong nTime);

/** Public: Construct regular time series. */
void ObitTimeFilterGridTime (ObitTimeFilter *in, olong seriesNo,
			     ofloat dTime, olong nTime, ofloat *times, ofloat *data);
typedef void (*ObitTimeFilterGridTimeFP) (ObitTimeFilter *in, olong seriesNo,
			     ofloat dTime, olong nTime, ofloat *times, ofloat *data);

/** Public: Copy time series to external times. */
void ObitTimeFilterUngridTime (ObitTimeFilter *in, olong seriesNo,
			     olong nTime, ofloat *times, ofloat *data);
typedef void (*ObitTimeFilterUngridTimeFP) (ObitTimeFilter *in, olong seriesNo,
			     olong nTime, ofloat *times, ofloat *data);

/** Public: Compute frequency series. */
void ObitTimeFilter2Freq (ObitTimeFilter *in);
typedef void (*ObitTimeFilter2FreqFP) (ObitTimeFilter *in);

/** Public: Compute Time series. */
void ObitTimeFilter2Time (ObitTimeFilter *in);
typedef void (*ObitTimeFilter2TimeFP) (ObitTimeFilter *in);

/** Public: Apply Filter to Frequency series */
void ObitTimeFilterFilter (ObitTimeFilter *in, olong seriesNo,
			   ObitTimeFilterType type, ofloat *parms,
			   ObitErr *err);
typedef void (*ObitTimeFilterFilterFP) (ObitTimeFilter *in, olong seriesNo,
					ObitTimeFilterType type, ofloat *parms,
					ObitErr *err);

/** Public: Apply Filter to Frequency series with physical parameters */
void ObitTimeFilterDoFilter (ObitTimeFilter *in, olong seriesNo,
			     ObitTimeFilterType type, ofloat *freqs,
			     ObitErr *err);
typedef void (*ObitTimeFilterDoFilterFP) (ObitTimeFilter *in, olong seriesNo,
					  ObitTimeFilterType type, ofloat *freqs,
					  ObitErr *err);

/** Public: Plot power spectrum. */
void ObitTimeFilterPlotPower (ObitTimeFilter *in, olong seriesNo,
			      gchar *label, ObitErr *err);
typedef void (*ObitTimeFilterPlotPowerFP) (ObitTimeFilter *in, olong seriesNo,
					   gchar *label, ObitErr *err);

/** Public: Plot Time series. */
void ObitTimeFilterPlotTime (ObitTimeFilter *in, olong seriesNo,
			     gchar *label, ObitErr *err);
typedef void (*ObitTimeFilterPlotTimeFP) (ObitTimeFilter *in, olong seriesNo,
					  gchar *label, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitTimeFilterClassDef.h"
} ObitTimeFilterClassInfo; 

#endif /* OBITTIMEFILTER_H */ 
