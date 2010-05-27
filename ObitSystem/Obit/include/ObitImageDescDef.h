/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2009                                          */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitImage structure            */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitImageDef.h
 * ObitImage structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
  /** Number of bits per pixel on disk. negative implies IEEE floating.*/
  olong    bitpix;
  /** Number of axes. */
  olong    naxis;
  /** Dimensions of axes. */
  olong   inaxes[IM_MAXDIM];
  /** WCS labels for each dimension of array. */
  gchar   ctype[IM_MAXDIM][IMLEN_KEYWORD];
  /** Axis coordinate increments. */
  ofloat  cdelt[IM_MAXDIM]; 
  /** Axis reference pixels (1-rel) */
  ofloat  crpix[IM_MAXDIM];
  /** Axis rotation angles (deg) */
  ofloat  crota[IM_MAXDIM];
  /** Alternate reference Pixel (frequency or velocity) */
  ofloat altCrpix;
  /** Offset in X (rotated RA) of phase center (deg)  */
  ofloat xshift;
  /** Offset in Y (rotated Dec) of phase center (deg)  */
  ofloat yshift;
  /** Convolving beam major axis in degrees  */
  ofloat beamMaj;
  /** Convolving beam minor axis in degrees   */
  ofloat beamMin;
  /** Convolving beam position angle in degrees */
  ofloat beamPA;
  /** Name of object. */
  gchar   object[IMLEN_VALUE]; 
  /** Name of telescope making observation. */
  gchar   teles[IMLEN_VALUE]; 
  /** Name of instrument making observation. */
  gchar   instrument[IMLEN_VALUE]; 
  /** Observer */
  gchar   observer[IMLEN_VALUE];
  /** Observing date as yyyy-mm-dd */
  gchar   obsdat[IMLEN_VALUE];
  /** Image creation date as yyyy-mm-dd */
  gchar   date[IMLEN_VALUE];
  /** Origin (software) of image. */
  gchar   origin[IMLEN_VALUE];
  /** Units of image */
  gchar   bunit[IMLEN_VALUE];
  /** Epoch (years) of celestial coordinates,
   *  This is sometimes confused with equinox.
   */
  ofloat  epoch;
  /** Mean Equinox of celestial coordinates (e.g. 2000) */
  ofloat  equinox;
  /** Observed RA (deg) */
  odouble obsra; 
  /** Observed Dec (deg) */
  odouble obsdec;
  /** Axis coordinate values at reference pixel. */
  odouble crval[IM_MAXDIM];
  /** Alternate reference value (frequency or velocity) */
  odouble altRef;
  /** Rest frequency of line (Hz)  */
  odouble restFreq;
  /** maximum value in image */
  ofloat  maxval;
  /** minimum value in image */
  ofloat  minval;
  /** Are there blanked pixels in the image? */
  gboolean  areBlanks;
  /** number of iterations of deconvolution */
  olong niter;
  /** Velocity reference frame: 1-3 (LSR, Helio, Observer) */
  olong VelReference;
  /** Velocity definition 0=>optical, 1=radio */
  olong VelDef;
  /** Coordinate type */
  ObitCoordType coordType;
  /** current row  (1-rel) read/write in buffer  */
  olong row;
  /** current plane (1-rel) read/write in buffer */
  olong plane;
  /** Access size code  */
  ObitIOSize IOsize;
  /** 0-rel axis order: RA (or longitude) , -1 means not present. */
  olong jlocr;
  /** 0-rel axis order: declination (or latitude), -1 means not present. */
  olong jlocd;
  /** 0-rel axis order: Stokes' parameters, -1 means not present. */
  olong jlocs;
  /** 0-rel axis order: Frequency. -1 means not present. */
  olong jlocf;
  /** 0-rel axis order: IF */
  olong jlocif;
  /** InfoList for other keywords */
  ObitInfoList *info;
  /** Is thi image using 3D imaging? */
  gboolean  do3D;
  /** x, y pixel offset if not do3D */
  ofloat  xPxOff, yPxOff;
