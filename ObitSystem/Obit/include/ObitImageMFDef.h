/* $Id: ObitImageMFDef.h 128 2009-09-23 14:48:29Z bill.cotton $  */
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
/*;  Correspondence concerning Obit should be addressed as follows:   */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitImageMF structure      */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitImageMFDef.h
 * ObitImageMF structure members for derived classes.
 */
#include "ObitImageDef.h"  /* Parent class definitions */
/** Maximum order of the imaging 
    Spectral index only = 1, plus curvature = 2 */
olong maxOrder;
/** Current order of the imaging 
    Spectral index only = 1, plus curvature = 2 */
olong curOrder;
/** Reference frequency */
odouble refFreq;
/** If TRUE, image has been formed but work 
    files/Decomposition not done */
gboolean fresh;
/** Number of coarse frequency planes */
olong nSpec;
/** Arrays of start and finish IFs (0-rel), per coarse channel */
olong *BIFSpec, *EIFSpec;
/** Arrays of start and finish Channels (0-rel), per coarse channel */
olong *BChanSpec, *EChanSpec;
/** Arrays of Center Frequency, per coarse channel */
odouble *specFreq;
/** Spectral index correction applied to data making image */
ofloat alpha;
