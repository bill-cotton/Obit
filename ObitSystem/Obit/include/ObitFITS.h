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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITFITS_H 
#define OBITFITS_H 
#include <glib.h>
#include "ObitImageDesc.h"
#include "ObitFile.h"
#include "ObitIO.h"

/**
 * \file ObitFITS.h
 * ObitFITS  class definition.
 *
 * This is a Utility class to handle the interface with FITS.
 * This class is non-derivable and only one instance is allowed.
 * Information regarding FITS is stored in a file static structure
 * and is available from function calls.
 * The structure must be initialized by a call to #ObitFITSClassInit.
 */

/*--------------  FITS types  -------------------------------------*/
/** Typedef for INTEGER in FITS structures - should be size of float */
typedef gint32 FITSint;

/*-------------- enumerations -------------------------------------*/
/*------------------  Macros    -------------------------------------*/
/** The maximum number of FITS "disks" */ 
#define MAXFITSDISK 20

/*--------------Class definitions-------------------------------------*/
/** ObitFITS Class Structure. */  
typedef struct {
  /** class name for verification */
  gchar *className;
  /** Object reference count. */
  gboolean initialized; 
  /** Array of directory strings by "disk" */
  gchar *FITSdir[MAXFITSDISK];
  /** Number of actual disks */
  olong NumberDisks;
} ObitFITS;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitFITSClassInit (gint number, gchar* dir[]);

/** Public: Constructor. */
ObitFITS* newObitFITS (void);

/** Public: Shutdown. */
void ObitFITSShutdown (void);

/** Public: Add FITS data directory. */
olong ObitFITSAddDir (gchar* dir, ObitErr *err);

/** Public: Replace FITS data directory path. */
void ObitFITSSetDir (gchar* dir, gint disk, ObitErr *err);

/** Public: Generate name of FITS file */
gchar* 
ObitFITSFilename (gint disk, gchar* fileName, ObitErr *err);

/** Public: Get directory string */
gchar* ObitFITSDirname (gint disk, ObitErr *err);

/** Public: Assign a scratch file info */
void ObitFITSAssign(gchar *pgmName, olong pgmNumber, 
		    olong disk, olong scrNo, ObitInfoList *info, 
		    ObitErr *err);

/** Public: Rename a FITS file */
void ObitFITSRename (ObitIO *in, ObitInfoList *info, ObitErr *err);

#endif /* OBITFITS_H */ 

