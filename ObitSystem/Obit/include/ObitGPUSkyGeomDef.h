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
 * \file ObitGPUSkyGeomDef.h
 * GPU (limited) version of ObitSkyGeom + ImageDesc and UVDesc structures
 */

#ifndef CUDASKYGEOMDEF_H 
#define CUDASKYGEOMDEF_H

#ifndef CUDAIMAGEDESCDEF_H // Prevent multiple definitions
#define CUDAIMAGEDESCDEF_H
/* Image Descriptor */
#define IM_MAXDIM 7       /* maximum array dimension */
#define IMLEN_VALUE 41    /* Maximum length of descriptor string value */
#define IMLEN_KEYWORD 21  /* Maximum length of descriptor keyword  */
typedef struct {
#include "CUDAImageDescDef.h"
} CUDAImageDesc;
#endif /* CUDAIMAGEDESCDEF_H */ 

/* Now UV descriptor */
#define UV_MAXDIM 7       /* maximum array dimension */
#define UV_MAX_RANP 14    /* maximum random parameter dimension */
#define UVLEN_VALUE 41    /* Maximum length of descriptor string value */
#define UVLEN_KEYWORD 21  /* Maximum length of descriptor keyword  */

typedef struct {
  /** Number of visibilities total */
  long nvis;
  /** Length (floats) of a vis record. */
  long lrec;
  /** Number of random parameters.*/
  long  nrparm;
  /** current beginning visibility (1-rel) read/write in buffer  */
  long firstVis;
  /** number of visibilities in buffer */
  long numVisBuff;
  /** Number of axes. */
  long   naxis;
  /** Number of subarrays in the data */
  long numSubA;
  /** Maximum antenna in the data */
  long maxAnt;
  /** Dimensions of axes. */
  long   inaxes[UV_MAXDIM];
  /** WCS labels for each dimension of array. */
  char   ctype[UV_MAXDIM][UVLEN_KEYWORD];
  /** WCS labels for random axes. */
  char   ptype[UV_MAX_RANP][UVLEN_KEYWORD];
  /** Axis coordinate increments. */
  float  cdelt[UV_MAXDIM]; 
  /** Axis reference pixels (1-rel) */
  float  crpix[UV_MAXDIM];
  /** Axis rotation angles (deg) */
  float  crota[UV_MAXDIM];
  /** Alternate reference Pixel (frequency or velocity) */
  float altCrpix;
  /** Offset in X (rotated RA) of phase center  (deg) */
  float xshift;
  /** Offset in Y (rotated Dec) of phase center (deg)  */
  float yshift;
  /** Convolving beam major axis in degrees  */
  float beamMaj;
  /** Convolving beam minor axis in degrees   */
  float beamMin;
  /** Convolving beam position angle in degrees */
  float beamPA;
  /** Data integration time in Days */
  float DeltaTime;
  /** Name of object. */
  char   object[UVLEN_VALUE]; 
  /** Name of telescope making observation. */
  char   teles[UVLEN_VALUE]; 
  /** Name of instrument making observation. */
  char   instrument[UVLEN_VALUE]; 
  /** Observer */
  char   observer[UVLEN_VALUE];
  /** Observing date as yyyy-mm-dd */
  char   obsdat[UVLEN_VALUE];
  /** File creation date as yyyy-mm-dd */
  char   date[UVLEN_VALUE];
  /** Origin (software) of image. */
  char   origin[UVLEN_VALUE];
  /** Units of visibilities */
  char   bunit[UVLEN_VALUE];
  /** Epoch (years) of celestial coordinates,
   *  This is sometimes confused with equinox. */
  float  epoch;
  /** Mean Equinox of celestial coordinates (e.g. 2000) */
  float  equinox;
  /** Maximum baseline length if > 0 */
  float maxBL;
  /** Maximum abs(W) if > 0 */
  float maxW;
  /** Julian date of observations */
  double JDObs;
  /** Axis coordinate values at reference pixel. */
  double crval[UV_MAXDIM];
  /** Observed RA (deg) */
  double obsra; 
  /** Observed Dec (deg) */
  double obsdec;
  /** Alternate reference value (frequency or velocity) */
  double altRef;
  /** Rest frequency of line (Hz)  */
  double restFreq;
  /** Reference Frequency (Hz) for u,v,w */
  double freq; 
  /** Velocity reference frame: 1-3 (LSR, Helio, Observer) */
  long VelReference;
  /** Velocity definition 0=>optical, 1=radio */
  long VelDef;
  /** Random parameter offset: U coordinate */
  long ilocu;
  /** Random parameter offset: V coordinate */
  long ilocv;
  /** Random parameter offset: W coordinate */
  long ilocw;
  /** Random parameter offset:Time */
  long iloct;
  /** Random parameter offset: Baseline */
  long ilocb;
  /** Random parameter offset: Antenna 1 */
  long iloca1;
  /** Random parameter offset: Antenna 2 + subarray*0.01 */
  long iloca2;
  /** Random parameter offset: subarray */
  long ilocsa;
  /** Random parameter offset: Source id. */
  long ilocsu;
  /** Random parameter offset:Freq id.  */
  long ilocfq;
  /** Random parameter offset: Integration time */
  long ilocit;
  /** Random parameter offset: Correlation id. */
  long ilocid;
  /** 0-rel axis order: Weight-scale parameters for compressed data */
  long ilocws;
  /** 0-rel axis order: Complex values */
  long jlocc;
  /** 0-rel axis order: Stokes' parameters */
  long jlocs;
  /** 0-rel axis order: Frequency */
  long jlocf;
  /** 0-rel axis order: RA */
  long jlocr;
  /**0-rel axis order: declination */
  long jlocd;
  /** 0-rel axis order: IF */
  long jlocif;
  /** Increment in data: Stokes */
  long incs;
  /** Increment in data: Frequency */
  long incf;
  /** Increment in data: IF */
  long incif;
  /** Total number of complex visibilities */
  long ncorr;
  /** Sort order 1st 2 char meaningful. */
  char isort[3];
} CUDAUVDesc;
#endif  /* CUDASKYGEOMDef_H */

