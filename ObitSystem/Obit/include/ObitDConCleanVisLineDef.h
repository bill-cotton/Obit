/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013,2022                                          */
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
/*  Define the basic components of the ObitDConCleanVisLine structure */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitDConCleanVisLineDef.h
 * ObitDConCleanVisLine structure members for this and any derived classes.
 */
#include "ObitDConCleanVisDef.h"  /* Parent class definitions */
/** Number of parallel channels */
olong nPar;
/** Number of thread arguments */
olong nArgs;
/** Array of Channel thread arguments */
gpointer* chArgs;
/** Secret pointer to top level Clean object in channel specific versions 
    NEVER Ref or esp. Unref */
Obit *motherShip;
/** List of channel "done" flags, [nArgs] */
gboolean *chDone;
