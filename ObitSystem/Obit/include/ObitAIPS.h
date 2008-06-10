/* $Id: ObitAIPS.h,v 1.10 2007/09/11 12:38:38 bcotton Exp $  */
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
#ifndef OBITAIPS_H 
#define OBITAIPS_H 
#include <glib.h>
#include "Obit.h"
#include "ObitImageDesc.h"
#include "ObitFile.h"
#include "ObitTypes.h"
#include "ObitIO.h"

/**
 * \file ObitAIPS.h
 * ObitAIPS  class definition.
 *
 * This is a Utility class to handle the interface with AIPS.
 * This class is non-derivable and only one instance is allowed.
 * Information regarding AIPS is stored in a file static structure
 * and is available from function calls.
 * The structure must be initialized by a call to #ObitAIPSClassInit.
 */

/*--------------  AIPS types  -------------------------------------*/
/** Typedef for INTEGER in AIPS structures - should be size of float */
typedef oint AIPSint;

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitAIPSFileType.
 * enum for AIPS file types.
 */
enum obitAIPSFileType {
  /** Catalog file */
  OBIT_AIPS_Catalog = 0, 
  /** Header file */
  OBIT_AIPS_Header,
  /** Image file */
  OBIT_AIPS_Image,
  /** UV data */
  OBIT_AIPS_UVdata,
  /** Scratch file */
  OBIT_AIPS_Scratch,
  /** Table file */
  OBIT_AIPS_Table,
  /** History */
  OBIT_AIPS_History,
  /** Plot */
  OBIT_AIPS_Plot,
  /** Slice */
  OBIT_AIPS_Slice
}; 
/** typedef for enum for AIPS file types */
typedef enum obitAIPSFileType ObitAIPSFileType;

/*------------------  Macros    -------------------------------------*/
/** The maximum number of AIPS "disks" */ 
#define MAXAIPSDISK 35

/*--------------Class definitions-------------------------------------*/
/** ObitAIPS Class Structure. */  
typedef struct {
  /** class name for verification */
  gchar *className;
  /** Object reference count. */
  gboolean initialized; 
  /** Array of directory strings by "disk" */
  gchar *AIPSdir[MAXAIPSDISK];
  /** Array of flags telling if scratch files allowed */
  gboolean noScrat[MAXAIPSDISK];
  /** Number of actual disks */
  olong NumberDisks;
  /** FORTRAN values of True, False */
  oint F_TRUE, F_FALSE;
} ObitAIPS;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitAIPSClassInit (gint number, gchar* dir[], oint F_TRUE, oint F_FALSE);

/** Public: Shutdown. */
void ObitAIPSShutdown (void);

/** Public: Generate name of AIPS file */
gchar* 
ObitAIPSFilename (ObitAIPSFileType type, olong disk, olong cno, 
		  olong userid, gchar *tabType, olong tabVer, ObitErr *err);

/** Public: Set directory string */
olong ObitAIPSSetDirname (gint disk, gchar* dir, ObitErr *err);

/** Public: Get directory string */
gchar* ObitAIPSDirname (gint disk, ObitErr *err);

/** Public: Get number of AIPS disks */
olong ObitAIPSGetNumDisk (ObitErr *err);

/** Public: Determine file position offset in an image */
ObitFilePos ObitAIPSImageFileOffset (gint naxis, olong *naxes, olong *pos);

/** Public: Determine file position offset in a table */
ObitFilePos ObitAIPSTableFileOffset (ObitFilePos start, olong lrow, olong row);

/** Public: Determine file position of the end of an AIPS table */
ObitFilePos ObitAIPSTableEOF (ObitFilePos start, olong lrow, olong nrow);

/** Public: Wonky padding for end of UV data */
ObitFilePos ObitAIPSUVWonkyPad (ObitFilePos curPos);

/** Public: Convert a olong to EHex */
void ObitAIPSEHex (gint in, gchar *out);

/** Public: Assign a scratch file info */
void ObitAIPSAssign(gchar *pgmName, olong pgmNumber, gchar *type,
		    olong user, olong disk, olong scrNo, ObitInfoList *info, 
		    ObitErr *err);

/** Public: Rename an AIPS file */
void ObitAIPSRename (ObitIO *in, ObitInfoList *info, ObitErr *err);

/** Public: Convert Fortran LOGICAL to a gboolean */
gboolean ObitAIPSBooleanF2C (oint logical);

/** Public: Convert a gboolean to a Fortran LOGICAL  */
oint ObitAIPSBooleanC2F (gboolean bool);

/** Public: Mark/Unmark AIPS directory as noScrat */
void ObitAIPSnoScrat(gint disk, gboolean noScrat, ObitErr *err);

/** Public: Tell if AIPS directory is noScrat */
gboolean ObitAIPSisNoScrat(gint disk);

/** Public: Check for "noScrat" AIPS disks */
void ObitAIPSSetnoScrat(ObitInfoList *info, ObitErr *err);


#endif /* OBITAIPS_H */ 

