/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004,2010                                          */
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
/*  Define the basic components of the ObitDConCleanPxList ClassInfo structure */
/* This is intended to be included in a classInfo structure definition  */
#include "ObitClassDef.h"  /* Parent class ClassInfo definition file */
/** Function pointer to Create/initialize ObitDCon structures */
ObitDConCleanPxListCreateFP ObitDConCleanPxListCreate;
/** Function pointer to Get Parameters */
ObitDConCleanPxListGetParmsFP ObitDConCleanPxListGetParms;
/** Function pointer to  Reset Clean */
ObitDConCleanPxListResetFP ObitDConCleanPxListReset;
/** Function pointer to Resize Arrrays */
ObitDConCleanPxListResizeFP ObitDConCleanPxListResize;
/** Function pointer to Update with new image and window */
ObitDConCleanPxListUpdateFP ObitDConCleanPxListUpdate;
/** Function pointer to Do minor cycle BGC CLEANing. */
ObitDConCleanPxListCLEANFP ObitDConCleanPxListCLEAN;
/** Function pointer to Do SDI CLEANing */
ObitDConCleanPxListSDIFP ObitDConCleanPxListSDI;
/** Function pointer to Get results of CLEAN  */
ObitDConCleanPxListResultFP ObitDConCleanPxListResult;
