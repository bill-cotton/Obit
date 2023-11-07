/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2020,2023                                          */
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
/*  Define the basic components of the ObitMatx ClassInfo structure   */
/* This is intended to be included in a classInfo structure definition*/
#include "ObitClassDef.h"  /* Parent class ClassInfo definition file */
/** Function pointer to  Constructor. */
ObitMatxCreateFP ObitMatxCreate;
/** Function pointer to compatibility checker. */
ObitMatxIsCompatibleFP ObitMatxIsCompatible;
/** Function pointer to Multiply. */
ObitMatxMultFP ObitMatxMult;
/** Function pointer to Multiply by conjugate transpose. */
ObitMatxMultCTFP ObitMatxMultCT;
/** Function pointer to Add. */
ObitMatxAddFP ObitMatxAdd;
/** Function pointer to Subract. */
ObitMatxSubFP ObitMatxSub;
/** Function pointer to Conjugate transpose. */
ObitMatxCTransFP ObitMatxCTrans;
/** Function pointer to Zero values. */
ObitMatxZeroFP ObitMatxZero;
/** Function pointer to write unit matrix. */
ObitMatxUnitFP ObitMatxUnit;
/** Function pointer to set 2x2 complex values. */
ObitMatxSet2CFP ObitMatxSet2C;
/** Function pointer to get 2x2 complex values. */
ObitMatxGet2CFP ObitMatxGet2C;
/** Function pointer to Inverse perfect linear feed Jones matrix */
ObitMatxIPerfLinJonesFP ObitMatxIPerfLinJones;
/** Function pointer to Outer 2x2 complex multiply */
ObitMatxOuterMult2CFP ObitMatxOuterMult2C;
/** Function pointer to 4x4 complex matrix * 4x1 complex vector multiply */
ObitMatxVec4MultFP ObitMatxVec4Mult;
