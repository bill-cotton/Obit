/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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
/*  Define the basic components of the ObitRMFit structure      */
/**
 * \file ObitRMFitDef.h
 * ObitRMFit structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** Number of lambda squared data points to be fitted */
olong nlamb2;
/** Number of terms to be fitted */
olong nterm;
/** Do Error analysis: */
gboolean doError;
/** Size of planes in pixels */
olong nx, ny;
/** min Q/U SNR */
ofloat minQUSNR;
/** min fraction of planes included */
ofloat minFrac;
/** Output Image descriptor */
ObitImageDesc *outDesc;
/** Array of Q, U pixel arrays for input data (nlamb2) */
ObitFArray **inQFArrays, **inUFArrays;
/** Array of Q, U RMSes per inFArrays (nlamb2) */
ofloat *QRMS, *URMS;
/** reference lambda squared */
odouble refLamb2;
/** Array of lambda squares (nlamb2) */
odouble *lamb2;
/** Array of pixel arrays for output data (2*nterm) */
ObitFArray **outFArrays;
