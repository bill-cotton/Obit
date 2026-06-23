/* $Id: $   */
/*HIDE c (esp.glib) structures from cuda*/
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2026                                               */
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
/*  Define the basic components of the CUDA ImageDesc and UVDesc structures */
/*  This is intended to be included in a class structure definition   */
/**
 * \file CUDAImageDescDef.h
 * Contents of CUDAImageDesc structure
 */

  /** Number of bits per pixel on disk. negative implies IEEE floating.*/
  long    bitpix;
  /** Number of axes. */
  long    naxis;
  /** Dimensions of axes. */
  long   inaxes[IM_MAXDIM];
  /** WCS labels for each dimension of array. */
  char   ctype[IM_MAXDIM][IMLEN_KEYWORD];
  /** Axis coordinate increments. */
  float  cdelt[IM_MAXDIM]; 
  /** Axis reference pixels (1-rel) */
  float  crpix[IM_MAXDIM];
  /** Axis rotation angles (deg) */
  float  crota[IM_MAXDIM];
  /** Alternate reference Pixel (frequency or velocity) */
  float altCrpix;
  /** Offset in X (rotated RA) of phase center (deg)  */
  float xshift;
  /** Offset in Y (rotated Dec) of phase center (deg)  */
  float yshift;
  /** Convolving beam major axis in degrees  */
  float beamMaj;
  /** Convolving beam minor axis in degrees   */
  float beamMin;
  /** Convolving beam position angle in degrees */
  float beamPA;
  /** Name of object. */
  char   object[IMLEN_VALUE]; 
  /** Name of telescope making observation. */
  char   teles[IMLEN_VALUE]; 
  /** Name of instrument making observation. */
  char   instrument[IMLEN_VALUE]; 
  /** Observer */
  char   observer[IMLEN_VALUE];
  /** Observing date as yyyy-mm-dd */
  char   obsdat[IMLEN_VALUE];
  /** Image creation date as yyyy-mm-dd */
  char   date[IMLEN_VALUE];
  /** Origin (software) of image. */
  char   origin[IMLEN_VALUE];
  /** Units of image */
  char   bunit[IMLEN_VALUE];
  /** Epoch (years) of celestial coordinates,
      This is sometimes confused with equinox. */
  float  epoch;
  /** Mean Equinox of celestial coordinates (e.g. 2000) */
  float  equinox;
  /** Observed RA (deg) */
  double obsra; 
  /** Observed Dec (deg) */
  double obsdec;
  /** Axis coordinate values at reference pixel. */
  double crval[IM_MAXDIM];
  /** Alternate reference value (frequency or velocity) */
  double altRef;
  /** Rest frequency of line (Hz)  */
  double restFreq;
  /** maximum value in image */
  float  maxval;
  /** minimum value in image */
  float  minval;
  /** Are there blanked pixels in the image? */
  int  areBlanks;
  /** number of iterations of deconvolution */
  long niter;
  /** Velocity reference frame: 1-3 (LSR, Helio, Observer) */
  long VelReference;
  /** Velocity definition 0=>optical, 1=radio */
  long VelDef;
  /** Coordinate type */
  //Fooey CUDACoordType coordType;
  /** current row  (1-rel) read/write in buffer  */
  long row;
  /** current plane (1-rel) plus higher order planes read/write in buffer */
  long plane, plane4, plane5, plane6, plane7;
  /** 0-rel axis order: RA (or longitude) , -1 means not present. */
  long jlocr;
  /** 0-rel axis order: declination (or latitude), -1 means not present. */
  long jlocd;
  /** 0-rel axis order: Stokes' parameters, -1 means not present. */
  long jlocs;
  /** 0-rel axis order: Frequency. -1 means not present. */
  long jlocf;
  /** 0-rel axis order: IF */
  long jlocif;
  /** Is this image using 3D imaging? */
  int  do3D;
  /** x, y pixel offset if not do3D */
  float  xPxOff, yPxOff;
