/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITAIPSFORTRAN_H 
#define OBITAIPSFORTRAN_H 
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitAIPS.h"
#include "ObitAIPSObject.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitAIPSFortran.h
 * ObitAIPSFortran module definition.
 *
 * This is a Utility module with Fortran callable routines for AIPS.
 */

/*-------------- type definitions----------------------------------*/
/*-------------- enumerations -------------------------------------*/
/*---------------Public functions---------------------------*/
/** Public: Init Obit System */
void obintx_ (const gchar *pgmName, const oint *lenPgmName, const oint *pgmNumber, 
	      const oint *AIPSuser,
	      const oint *numberAIPSdisk, const gchar *AIPSdir, const oint *lenAIPSdir,
	      const oint *numberFITSdisk, const gchar *FITSdir, const oint *lenFITSdir,
	      const oint *F_TRUE, const oint *F_FALSE, oint *ierr);

/** Public: Shutdown Obit System */
void obshtx_ (void);

/** Public: Copy a uvdata set potentially applying calibration, editing and selection,*/
void obuvcp_ (AIPSObj uvin, AIPSObj uvout, oint *ierr);

/** Public: Determine and apply uniform weighting corrections to uv data */
void obufwt_ (AIPSObj uv, AIPSObj image, oint *ierr);

/** Public: Make beams and image from uv data */
void obimuv_ (AIPSObj uvdata, oint *ifield, oint *nfield, AIPSObj *image, AIPSObj *beam, 
	      oint *doBeam, oint *doCreate, oint *chan, oint* nchan,  oint* imchan, oint *ierr);

/** Public: Make beams and image with possible ionospheric corrections */
void obiuvi_ (AIPSObj uvdata, oint *ifield, oint *nfield, AIPSObj *image, AIPSObj *beam, 
	      oint *doBeam, oint *doCreate, oint *chan, oint* nchan,  oint* imchan, oint *ierr);


#endif /* OBITAIPSFORTRAN_H */ 

