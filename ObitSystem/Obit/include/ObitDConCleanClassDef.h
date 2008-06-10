/* $Id: ObitDConCleanClassDef.h,v 1.4 2005/04/19 11:46:30 bcotton Exp $  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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
/*  Define the basic components of the ObitDConClean ClassInfo structure */
/* This is intended to be included in a classInfo structure definition   */
#include "ObitDConClassDef.h"  /* Parent class ClassInfo definition file   */
/** Pointer to Set Default CLEAN windows.*/
ObitDConCleanDefWindowFP ObitDConCleanDefWindow;
/** Pointer to Prepare for minor cycle.*/
ObitDConCleanPixelStatsFP ObitDConCleanPixelStats;
/** Pointer to Determine image statistics.*/
ObitDConCleanImageStatsFP ObitDConCleanImageStats;
/** Pointer to Select components to be subtracted */
ObitDConCleanSelectFP ObitDConCleanSelect;
/** Pointer to Subtract components */
ObitDConCleanSubFP ObitDConCleanSub;
/** Pointer to Restore subtracted components.*/
ObitDConCleanRestoreFP ObitDConCleanRestore;
/** Pointer to Flatten multiple facets to one.*/
ObitDConCleanFlattenFP ObitDConCleanFlatten;
/** Pointer to Restore subtracted components from other fields.*/
ObitDConCleanXRestoreFP ObitDConCleanXRestore;
/** Pointer to .Automatically add window. */
ObitDConCleanAutoWindowFP ObitDConCleanAutoWindow;

