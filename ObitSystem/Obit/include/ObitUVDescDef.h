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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitUVDesc structure           */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitUVDescDef.h
 * ObitUVDesc structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Type of access to file */
  ObitIOAccess access;
  /** Number of visibilities total */
  olong nvis;
  /** Length (floats) of a vis record. */
  olong lrec;
  /** Number of random parameters.*/
  olong  nrparm;
  /** current beginning visibility (1-rel) read/write in buffer  */
  olong firstVis;
  /** number of visibilities in buffer */
  olong numVisBuff;
  /** Number of axes. */
  olong   naxis;
  /** Number of subarrays in the data */
  olong numSubA;
  /** Maximum antenna in the data */
  olong maxAnt;
  /** Dimensions of axes. */
  olong   inaxes[UV_MAXDIM];
  /** WCS labels for each dimension of array. */
  gchar   ctype[UV_MAXDIM][UVLEN_KEYWORD];
  /** WCS labels for random axes. */
  gchar   ptype[UV_MAX_RANP][UVLEN_KEYWORD];
  /** Axis coordinate increments. */
  ofloat  cdelt[UV_MAXDIM]; 
  /** Axis reference pixels (1-rel) */
  ofloat  crpix[UV_MAXDIM];
  /** Axis rotation angles (deg) */
  ofloat  crota[UV_MAXDIM];
  /** Alternate reference Pixel (frequency or velocity) */
  ofloat altCrpix;
  /** Offset in X (rotated RA) of phase center  (deg) */
  ofloat xshift;
  /** Offset in Y (rotated Dec) of phase center (deg)  */
  ofloat yshift;
  /** Data integration time in Days */
  ofloat DeltaTime;
  /** Name of object. */
  gchar   object[UVLEN_VALUE]; 
  /** Name of telescope making observation. */
  gchar   teles[UVLEN_VALUE]; 
  /** Name of instrument making observation. */
  gchar   instrument[UVLEN_VALUE]; 
  /** Observer */
  gchar   observer[UVLEN_VALUE];
  /** Observing date as yyyy-mm-dd */
  gchar   obsdat[UVLEN_VALUE];
  /** File creation date as yyyy-mm-dd */
  gchar   date[UVLEN_VALUE];
  /** Origin (software) of image. */
  gchar   origin[UVLEN_VALUE];
  /** Units of visibilities */
  gchar   bunit[UVLEN_VALUE];
  /** Epoch (years) of celestial coordinates,
   *  This is sometimes confused with equinox. */
  ofloat  epoch;
  /** Mean Equinox of celestial coordinates (e.g. 2000) */
  ofloat  equinox;
  /** Maximum baseline length if > 0 */
  ofloat maxBL;
  /** Maximum abs(W) if > 0 */
  ofloat maxW;
  /** Julian date of observations */
  odouble JDObs;
  /** Axis coordinate values at reference pixel. */
  odouble crval[UV_MAXDIM];
  /** Observed RA (deg) */
  odouble obsra; 
  /** Observed Dec (deg) */
  odouble obsdec;
  /** Alternate reference value (frequency or velocity) */
  odouble altRef;
  /** Rest frequency of line (Hz)  */
  odouble restFreq;
  /** Reference Frequency (Hz) for u,v,w */
  odouble freq; 
  /** Array of Frequencies (Hz) in the order of Freq, IF in data 
   *  dimension is nfreq*nif */
  odouble *freqArr;
  /** Array of differential frequency scaling factors from one frequency to the
   *  next, dimension is nfreq*nif */
  ofloat *fscale;
  /** Array of Frequencies per IF */
  odouble *freqIF;
  /** Array of channel frequency increments per IF */
  ofloat *chIncIF;
  /** Array of sidebands, one per IF */
  olong *sideband;
  /** Number of antennas (actually max ant number) per subarray */
  oint *numAnt;
  /** Velocity reference frame: 1-3 (LSR, Helio, Observer) */
  olong VelReference;
  /** Velocity definition 0=>optical, 1=radio */
  olong VelDef;
  /** Random parameter offset: U coordinate */
  olong ilocu;
  /** Random parameter offset: V coordinate */
  olong ilocv;
  /** Random parameter offset: W coordinate */
  olong ilocw;
  /** Random parameter offset:Time */
  olong iloct;
  /** Random parameter offset: Baseline */
  olong ilocb;
  /** Random parameter offset: Source id. */
  olong ilocsu;
  /** Random parameter offset:Freq id.  */
  olong ilocfq;
  /** Random parameter offset: Integration time */
  olong ilocit;
  /** Random parameter offset: Correlation id. */
  olong ilocid;
  /** 0-rel axis order: Weight-scale parameters for compressed data */
  olong ilocws;
  /** 0-rel axis order: Complex values */
  olong jlocc;
  /** 0-rel axis order: Stokes' parameters */
  olong jlocs;
  /** 0-rel axis order: Frequency */
  olong jlocf;
  /** 0-rel axis order: RA */
  olong jlocr;
  /**0-rel axis order: declination */
  olong jlocd;
  /** 0-rel axis order: IF */
  olong jlocif;
  /** Increment in data: Stokes */
  olong incs;
  /** Increment in data: Frequency */
  olong incf;
  /** Increment in data: IF */
  olong incif;
  /** Total number of complex visibilities */
  olong ncorr;
  /** Sort order 1st 2 char meaningful. */
  gchar isort[3];
  /** InfoList for other keywords */
  ObitInfoList *info;
