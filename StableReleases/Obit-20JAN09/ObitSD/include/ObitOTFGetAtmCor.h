/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITOTFGETATMCOR_H 
#define OBITOTFGETATMCOR_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitOTF.h"
#include "ObitTableOTFSoln.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFGetAtmCor.h
 *
 * Routines to determine atmospheric opacity calibration for the OTF class.
 * Routines determine calibration for an ObitOTF and writes ObitTableOTFSoln
 * tables.
 */

/*---------------Public functions---------------------------*/
/** Public: Determine  atmospheric opacity calibration. */
ObitTableOTFSoln* ObitOTFGetAtmCor (ObitOTF *in, ObitOTF *out, ObitErr *err);

/** Public: Fill Soln table with atmospheric corrections. */
ObitTableOTFSoln* ObitOTFGetAtmEm (ObitOTF *inOTF, ObitOTF *outOTF, ObitErr *err);

#endif /* OBITOTFGETATMCOR_H */ 
