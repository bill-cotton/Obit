/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2010                                          */
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
/*  Define the basic components of the ObitDConCleanVis ClassInfo structure */
/* This is intended to be included in a classInfo structure definition   */
#include "ObitDConCleanClassDef.h"  /* Parent class ClassInfo definition file */
/** Pointer to Determine quality measure for field.*/
ObitDConCleanVisQualityFP ObitDConCleanVisQuality;
/** Pointer to Determine if reimaging needed.*/
ObitDConCleanVisReimageFP ObitDConCleanVisReimage;
/** Pointer to Determine if reimaging needed.*/
ObitDConCleanVisAddFieldFP ObitDConCleanVisAddField;
/** Pointer to Determine if recentering needed.*/
ObitDConCleanVisRecenterFP ObitDConCleanVisRecenter;
/** Pointer to Filter weak, isolated components.*/
ObitDConCleanVisFilterFP ObitDConCleanVisFilter;
/** Pointer to Set Default CLEAN windows .*/
ObitDConCleanVisDefWindowFP ObitDConCleanVisDefWindow;
/* Private functions for derived classes */
/** Pointer to (re)make residuals */
MakeResidualsFP MakeResiduals;
/** Pointer to  (re)make all residuals. */
MakeAllResidualsFP MakeAllResiduals;
/** Pointer to Low accuracy subtract CLEAN model.*/
SubNewCCsFP SubNewCCs;
/** Pointer to Create/init PxList.*/
NewPxListFP NewPxList;
/** Pointer to Create/init Pixarray.*/
NewPxArrayFP NewPxArray;
/** Pointer to Delete Pixarray.*/
KillPxArrayFP KillPxArray;
/** Pointer to Delete BeamPatches.*/
KillBeamPatchesFP KillBeamPatches;
