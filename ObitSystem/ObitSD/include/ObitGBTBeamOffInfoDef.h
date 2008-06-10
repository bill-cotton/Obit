/* $Id: ObitGBTBeamOffInfoDef.h,v 1.1.1.1 2004/07/19 17:04:44 bcotton Exp $                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004                                               */
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
/*  Define the basic components of the ObitGBTBeamOffInfo structure.          */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitGBTBeamOffInfoListDef.h
 * ObitGBTBeamOffInfoDef structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
#define MAXNUMBEAMOFF 10
/** Number of BeamOffs in this structure */
olong nBeamOff;
/** Name */
gchar BeamName[MAXNUMBEAMOFF][8];
/** Azimuth*cos(el) offset in deg */
odouble xeloff[MAXNUMBEAMOFF];
/** elevation offset  in deg */
odouble eloff[MAXNUMBEAMOFF];
/** sig/ref feed 1 */
olong srfeed1[MAXNUMBEAMOFF];
/** sig/ref feed 2 */
olong srfeed2[MAXNUMBEAMOFF];
