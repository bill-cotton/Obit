/* $Id: ObitBeamShapeDef.h,v 1.1 2008/05/06 13:20:14 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008                                               */
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
/*  Define the basic components of the ObitBeamShape structure         */
/*  This is intended to be included in a class structure definition   */
/* and to be used as the template for generating new classes derived  */
/* from Obit.                                                         */
/**
 * \file ObitBeamShapeDef.h
 * ObitBeamShape structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** Image descriptor */
ObitImageDesc *myDesc;
/** Antenna pointing position (rad)*/
odouble raPnt, decPnt;
/** Reference Frequency */
odouble refFreq;
/** Minimum desired gain */
ofloat pbmin;
/** Antenna diameter (m) */
ofloat antSize;
/** Gain wanted? */
gboolean doGain;
/** Use Jinc or polynomial */
gboolean doJinc;
