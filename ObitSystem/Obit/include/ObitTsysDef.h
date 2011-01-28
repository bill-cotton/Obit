/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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
/*  Define the basic components of the ObitTsys structure             */
/*  This is intended to be included in a class structure definition   */
/* and to be used as the template for generating new classes derived  */
/* from Obit.                                                         */
/**
 * \file ObitTsysDef.h
 * ObitTsys structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** TY table */
ObitTableTY *TYTable;
/** TY table row */
ObitTableTYRow *TYRow;
/** Number of IFs */
olong numIF;
/** Number of Polns */
olong numPoln;
/** Number of subarrays */
olong numSubA;
/** Antenna list one per subarray */
ObitAntennaList **antList;
/** Current subarray */
olong SubA;
/** Number of antennas in current subarray */
olong numAnt;
/** Array of antenna prior row numbers */
olong *priorRowNo;
/** Array of antenna following row numbers */
olong *followRowNo;
/** Array of antenna prior TYRow */
ObitTableTYRow **priorRow;
/** Array of antenna follow TYRow */
ObitTableTYRow **followRow;
