/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
/*  Define the basic components of the ObitOTFDesc structure           */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitOTFDescDef.h
 * ObitOTFDesc structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Type of access to file */
ObitIOAccess access;
/** Source of data in file */
ObitGBTOTFType OTFType;
/** Number of records total */
olong nrecord;
/** Length (floats) of a record. */
olong lrec;
/** current beginning record (1-rel) read/write in buffer  */
olong firstRec;
/** number of records in buffer */
olong numRecBuff;
/** Number of columns in data */
olong ncol;
/** Number of descriptive parameters */
olong numDesc;
/** Column labels */
gchar   colType[OTF_MAX_COL][OTFLEN_KEYWORD];
/** Column units */
gchar   colUnit[OTF_MAX_COL][OTFLEN_VALUE];
/** Column repeat count */
olong    colRepeat[OTF_MAX_COL];
/** Name of object. */
gchar   object[OTFLEN_VALUE]; 
/** Name of telescope making observation. */
gchar   teles[OTFLEN_VALUE]; 
/** Observing date as yyyy-mm-dd */
gchar   obsdat[OTFLEN_VALUE];
/** File creation date as yyyy-mm-dd */
gchar   date[OTFLEN_VALUE];
/** Origin (software) of image. */
gchar   origin[OTFLEN_VALUE];
/** Units of data */
gchar   bunit[OTFLEN_VALUE];
/** Number of data axes. */
olong   naxis;
/** Dimensions of axes. */
olong   inaxes[OTF_MAXDIM];
/** WCS labels for each dimension of array. */
gchar   ctype[OTF_MAXDIM][OTFLEN_KEYWORD];
/** Axis coordinate increments. */
ofloat  cdelt[OTF_MAXDIM]; 
/** Axis reference pixels (1-rel) */
ofloat  crpix[OTF_MAXDIM];
/** Axis rotation angles (deg) */
ofloat  crota[OTF_MAXDIM];
/** Axis coordinate values at reference pixel. */
odouble crval[OTF_MAXDIM];
/** Alternate reference Pixel (frequency or velocity) */
ofloat altCrpix;
/** Epoch (years) of celestial coordinates,
 *  This is sometimes confused with equinox. */
ofloat  epoch;
/** Mean Equinox of celestial coordinates (e.g. 2000) */
ofloat  equinox;
/** Julian date of observations */
odouble JDObs;
/** Observed central RA (deg) */
odouble obsra; 
/** Observed central Dec (deg) */
odouble obsdec;
/** Beam size Gaussian FWHM in deg */
ofloat beamSize;
/** Telescope diameter in m */
ofloat diameter;
/** 0-rel index of time in record */
olong iloct;
/** 0-rel index of integration time in record */
olong ilocti;
/** 0-rel index of target index in record */
olong iloctar;
/** 0-rel index of scan index in record */
olong ilocscan;
/** 0-rel index of RA in record */
olong ilocra;
/** 0-rel index of Dec in record */
olong ilocdec;
/** 0-rel index of rotation in record */
olong ilocrot;
/** 0-rel index of cal in record */
olong iloccal;
/** 0-rel index of start of data array in record */
olong ilocdata;
/** 0-rel axis order: Data-Wt */
olong jlocdatawt;
/** 0-rel axis order: Feed */
olong jlocfeed;
/** 0-rel axis order: Stokes' parameters */
olong jlocs;
/** 0-rel axis order: Frequency */
olong jlocf;
/** 0-rel axis order: State */
olong jlocstate;
/** Increment in data: Data/Weight */
olong incdatawt;
/** Increment in data: Feed */
olong incfeed;
/** Increment in data: Stokes */
olong incs;
/** Increment in data: Frequency */
olong incf;
/** Increment in data: State */
olong incstate;
/** Sort order 1st 2 char meaningful. */
gchar isort[3];
/** InfoList for other keywords */
ObitInfoList *info;
