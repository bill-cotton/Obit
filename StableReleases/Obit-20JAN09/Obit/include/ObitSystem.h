/* $Id$  */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITSYSTEM_H 
#define OBITSYSTEM_H 
#include "ObitErr.h"
#include "Obit.h"
#include "ObitIO.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSystem.h
 * ObitSystem manager class
 *
 * This class is derived from the #Obit class.
 *
 * This system wide information and controls creation and destruction
 * of Scratch files.
 * There should be one and only one instance of ObitSystem.
 *
 * \section ObitSystemScratch Scratch files.
 * The ObitSystem assists in creating and maintains a registry of scratch 
 * files. 
 * \li ObitSystemGetScratch assigns naming information to an info for a 
 * scratch object and assigns whatever resources are needed 
 * (e.g. assign AIPS catalog slot)
 * \li ObitSystemAddScratch registers a scratch file.  Any scratch files 
 * remaining in the list when ObitSystemShutdown is called will be Unreffed 
 * causing their destruction.
 * \li ObitSystemFreeScratch removes an entry from the scratch file list.
 *
 * \section ObitSystemUsage Usage
 * Instances can be obtained using the #ObitSystemStartup constructor
 * and deleted using the #ObitSystemShutdown function.
 * Obit function calls should only be made between these two ObitSystem calls.
 */

/*------------------- Macroes ----------------------------------------*/
/*--------------Class definitions-------------------------------------*/
/**
 * ObitSystem Class structure.
 *
 * This class contains Obit system-wide information
 */  
typedef struct {
#include "ObitSystemDef.h" /* actual definition */
} ObitSystem;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSystemClassInit (void);

/** Public: startup/Constructor. */
ObitSystem* 
ObitSystemStartup (gchar *pgmName, olong pgmNumber,
		   olong AIPSuser,
		   olong numberAIPSdisk, gchar* AIPSdir[], 
		   olong numberFITSdisk, gchar* FITSdir[], 
		   oint F_TRUE, oint F_FALSE, ObitErr *err);
/** define type for ClassInfo structure */
typedef ObitSystem* 
(*ObitSystemStartupFP) (gchar *pgmName, olong pgmNumber, 
			olong AIPSuser,
			olong numberAIPSdisk, gchar* AIPSdir[], 
			olong numberFITSdisk, gchar* FITSdir[], 
			ObitErr *err);

/** Public: Return class pointer. */
gconstpointer ObitSystemGetClass (void);

/** Public: Shutdown */
ObitSystem* ObitSystemShutdown (ObitSystem* in);
typedef ObitSystem* (*ObitSystemShutdownFP) (ObitSystem* in);

/** Public: Get Scratch file assignments */
void ObitSystemGetScratch (ObitIOType FileType, gchar *type,
			   ObitInfoList *info, ObitErr *err);
typedef void (*ObitSystemGetScratchFP) (ObitIOType FileType, gchar *type,
			   ObitInfoList *info, ObitErr *err);

/** Public: Add Scratch file assignments */
void ObitSystemAddScratch (Obit *in, ObitErr *err);
typedef void (*ObitSystemAddScratchFP) (Obit *in, ObitErr *err);

/** Public: Free Scratch file assignments */
void ObitSystemFreeScratch (Obit *in, ObitErr *err);
typedef void (*ObitSystemFreeScratchFP) (Obit *in, ObitErr *err);

/** Public: Tell If Obit is initialized */
gboolean ObitSystemIsInit (void);

/** Public: Tell Program Name */
gchar* ObitSystemGetPgmName (void);

/** Public: Reset Program Name */
void ObitSystemSetPgmName (gchar *pgmName);

/** Public: Tell Program Number */
olong ObitSystemGetPgmNumber (void);

/** Public: Reset Program Number */
void ObitSystemSetPgmNumber (olong pgmNumber);

/** Public: Tell AIPS user ID */
olong ObitSystemGetAIPSuser (void);

/** Public: Reset AIPS user ID */
void ObitSystemSetAIPSuser (olong AIPSuser);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitSystemClassDef.h" /* Actual definition */
} ObitSystemClassInfo; 

#endif /* OBITSYSTEM_H */ 

