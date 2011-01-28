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
/*  Define the basic components of the ObitSwPower structure             */
/**
 * \file ObitSwPowerDef.h
 * ObitSwPower structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** SY table */
ObitTableSY *SYTable;
/** SY table row */
ObitTableSYRow *SYRow;
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
/** Array of antenna prior SYRow */
ObitTableSYRow **priorRow;
/** Array of antenna follow SYRow */
ObitTableSYRow **followRow;
