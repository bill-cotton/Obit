/* $Id: ObitSkyModelVMIonDef.h,v 1.1 2006/02/17 16:45:08 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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
/*  Define the basic components of the ObitSkyModelIon structure      */
/*  This class represents sky models and their Fourier transform      */
/*  This is intended to be included in a class structure definition.  */
/**
 * \file ObitSkyModelVMIonDef.h
 * ObitSkyModel structure members for this and any derived classes.
 */
#include "ObitSkyModelVMDef.h"  /* Parent class definitions */
/** Index table for fields with components */
olong *fieldIndex;
/** Table of rotation matrices for components 6 needed per row (field) */
ObitFArray *uRotTab;
/** Table of Zernike X gradient terms,  ncoef needed per row (field) */
ObitFArray *ZernY;
/** Table of Zernike Y gradient terms,  ncoef needed per row (field) */
ObitFArray *ZernX;
/** if true need 3D rotation multiply by matrix */
gboolean do3Dmul;
/** cosine, sine of rotation difference between uv, image */
ofloat ccrot, ssrot;
/** NI Table with ionospheric model */
ObitTableNI *NITable;
/** NI Table row */
ObitTableNIRow *NIRow;
/** Last row read */
olong lastNIRow;
/** Number of coefficients */
olong ncoef;
/** Time of preceeding coefficients */
ofloat priorTime;
/** Time of following coefficients */
ofloat followTime;
/** Weight of preceeding coefficients */
ofloat priorWeight;
/** Weight of following coefficients */
ofloat followWeight;
/** Array of coef for preceeding time */
ofloat *priorCoef;
/** Array of coef for following time */
ofloat *followCoef;
