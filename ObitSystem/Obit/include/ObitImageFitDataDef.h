/* $Id: ObitImageFitDataDef.h,v 1.2 2007/09/20 03:08:13 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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
/*  Define the basic components of the ObitImageFitData structure     */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitImageFitDataDef.h
 * ObitImageFitData structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** Function pointer for evaluation */
ObitImageFitDataFuncFP fx;
/**  image pixel array */
ObitFArray *pixels;
/**  image pixel residual array */
ObitFArray *resids;
/** Number of components to fit */
olong ncomp; 
/**  Total number of parameters */
olong nparm;  
/**  Number of parameters to fit */
olong nvar;  
/** Number of model parameters per component  */
olong  *np;
/** Component number of fitted parameters  */
olong  *ivar;
/** Parameter number of fitted parameters */
olong  *jvar;
/** Model component types */
ObitFitModelCompType *type;
/** Parameters of model components [comp][param] */
odouble **p;
/** Small but significant differential in p [comp][param] */
odouble **dp;
/**  Fit this parameter? [comp][param] */
gboolean **pf;
/** Model component parameter upper bounds, fblank = undef */
odouble **pu; 
/** Model components parameter lower bounds, fblank = undef */
odouble **pl; 
/** model component parameter errors */
odouble **e;
/** Work array the length of the number of all parameters */
odouble *ptemp;
/** Amplitude scaling */
ofloat rscale;
/** Restoring beam area in pixels */
ofloat beamarea;
/** Restoring beam in pixels */
ofloat beam[3];
/** RMS in full image */
ofloat irms;
