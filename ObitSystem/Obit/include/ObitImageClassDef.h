/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2015                                          */
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
/*  Define the basic components of the ObitImage ClassInfo structure  */
/* This is intended to be included in a classInfo structure definition*/
#include "ObitDataClassDef.h"  /* Parent class ClassInfo definition file */
/** Function pointer to scratch copy constructor. */
newObitImageScratchFP newObitImageScratch;
/** Function pointer to are two Images the same. */
ObitImageSameFP ObitImageSame;
/** Function pointer to read specified plane. */
ObitImageGetPlaneFP ObitImageGetPlane;
/** Function pointer to write specified plane. */
ObitImagePutPlaneFP ObitImagePutPlane;
/** Function pointer to Create an associated table. */
newObitImageTableFP newObitImageTable;
/** Destroy an associated Table(s) */
ObitImageZapTableFP ObitImageZapTable;
/** Fully instantiate. */
ObitImageFullInstantiateFP ObitImageFullInstantiate;
/** Copy associated Table(s) */
ObitImageCopyTablesFP ObitImageCopyTables;
/** Update disk resident tables information */
ObitImageUpdateTablesFP ObitImageUpdateTables;
/** Name beam */
ObitImageSetBeamNameFP ObitImageSetBeamName;
/** Get beam */
ObitImageGetBeamFP ObitImageGetBeam;
/** Get beam order */
ObitImageGetBeamOrderFP ObitImageGetBeamOrder;
/** Get plane Frequency */
ObitImageGetPlaneFreqFP ObitImageGetPlaneFreq;
