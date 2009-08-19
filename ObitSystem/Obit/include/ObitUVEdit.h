/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2009                                          */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITUVEDIT_H 
#define OBITUVEDIT_H 

#include "ObitErr.h"
#include "ObitUV.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVEdit.h
 * ObitUVEdit module definition for editing ObitUV data.
 *
 * This utility module contains utility functions for editing uv data.
 */

/*---------------Public functions---------------------------*/
/** Public:  Time domain editing, FG table out */
void ObitUVEditTD (ObitUV *inUV, ObitUV *outUV, ObitErr *err);

/** Public:  Time domain RMS/Avg editing, FG table out */
void ObitUVEditTDRMSAvg (ObitUV *inUV, ObitUV *outUV, ObitErr *err);

/** Public:  Time domain RMS/Avg vector editing, FG table out */
void ObitUVEditTDRMSAvgVec (ObitUV *inUV, ObitUV *outUV, ObitErr *err);

/** Public:  Frequency domain editing, FG table out */
void ObitUVEditFD (ObitUV *inUV, ObitUV *outUV, ObitErr *err);

/** Public: Edit by excessive amplitude in given Stokes, FG table out */
void ObitUVEditStokes (ObitUV *inUV, ObitUV *outUV, ObitErr *err);


/** Public: Clip raw visibilities, uv data out */
ObitUV* ObitUVEditClip (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			ObitErr *err);

/** Public: Flag visibilities by Stokes, uv data out */
ObitUV* ObitUVEditClipStokes (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			      ObitErr *err);

/** Public: Flag visibilities by running median */
void ObitUVEditMedian (ObitUV *inUV, ObitUV *outUV, ObitErr *err);

#endif /* OBITIUVEDIT_H */ 
