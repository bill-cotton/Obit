/* $Id$   */
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
#ifndef OBITAIPSDIR_H 
#define OBITAIPSDIR_H 
#include <glib.h>
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitAIPS.h"

/**
 * \file ObitAIPSDir.h
 * ObitAIPSDir  class definition.
 *
 * This is a Utility class to handle the interface with the AIPS catalog
 * directory system.
 * This class is non-derivable and has no public instances.
 */
/*------------------ Structures -----------------------------*/
/**
 * AIPS Catalog (directory) entry structure
 * A short description of the name and access info.
 */
typedef struct {
  /** User id or -1 if slot is unused */
  olong user;
  /** status of entry:
   * \li 0   => no programs are accessing,
   * \li n>0 => n programs have read access,
   * \li -1  =>one program with write access,
   * \li n<0 => one program writing, 1+n reading.
   */
  gint32 status;
  /** date and time of last access in packed format
   * (Packed years since 1900):
   * \li access[0] = 256 * (256 * YY + MM) + DD
   * \li access[1] = 256 * (256 * HH + MM) + SS
   */
  guint32 access[2];
  /** AIPS sequence number */
  gint32 seq;
  /** AIPS name (12 char, not NULL terminated) */
  gchar name[12];
  /** AIPS class (6 char, not NULL terminated) */
  gchar class[6];
  /** AIPS file type, "MA","UV", "SC" (2 char) */
  gchar type[2];
} ObitAIPSDirCatEntry;

/*----------------- enums ---------------------------*/
/** enum for Catalog status change codes. */
enum _ObitAIPSDirStatusCode {
  /** Add Write Status */
  OBIT_AIPS_Dir_AddWrite,
  /** Clear Write Status */
  OBIT_AIPS_Dir_ClearWrite,
  /** Increment Read Status */
  OBIT_AIPS_Dir_IncRead,
  /** Decrement Read Status */
  OBIT_AIPS_Dir_DecRead,
  /** Clear All Status */
  OBIT_AIPS_Dir_ClearAll
};/* end enum _ObitAIPSDirStatusCode */
/** typedef for enum for status codes */
typedef enum _ObitAIPSDirStatusCode ObitAIPSDirStatusCode;

/** enum for Catalog status change errors. */
enum _ObitAIPSDirStatusError {
  /** Operation OK */
  OBIT_AIPS_Dir_StatusOK,
  /** Input Specification Error */
  OBIT_AIPS_Dir_StatusSpecErr,
  /** I/O Error */
  OBIT_AIPS_Dir_StatusIOErr,
  /** Cannot add write due to read status */
  OBIT_AIPS_Dir_StatusRead,
  /** Cannot add read or write due to write status */
  OBIT_AIPS_Dir_StatusWrite
};/* end enum _ObitAIPSDirStatusError */
/** typedef for enum for error codes */
typedef enum _ObitAIPSDirStatusError ObitAIPSDirStatusError;

/*---------------Public functions---------------------------*/
/** Public: Find AIPS Catalog for a given AIPS name ... */
olong ObitAIPSDirFindCNO(olong disk, olong user, 
		     gchar Aname[13], gchar Aclass[7], gchar Atype[3], 
		     olong seq, ObitErr *err);

/** Public: Allocate AIPS directory slot an fill it in. */
olong ObitAIPSDirAlloc(olong disk, olong user, 
		     gchar Aname[13], gchar Aclass[7], gchar Atype[3], 
		     olong seq, gboolean *exist, ObitErr *err);

/** Public: Remove Catalog directory entry */
void ObitAIPSDirRemoveEntry(olong disk, olong user, olong cno, ObitErr *err);

/** Public: Determine maximum catalog slot number occupied */
olong ObitAIPSDirNumber(olong disk, olong user, ObitErr *err);

/** Public: Determine maximum sequence number used */
olong ObitAIPSDirHiSeq(olong disk, olong user, gchar Aname[13], 
		      gchar Aclass[7], gchar Atype[3], 
		      gboolean exist, ObitErr *err);

/** Public: Rename entry */
void ObitAIPSDirRename(olong disk, olong user,  olong cno, gchar *newName, 
		      gchar *newClass, olong newSeq, ObitErr *err);

/** Public: Get Catalog directory entry */
ObitAIPSDirCatEntry* 
ObitAIPSDirGetEntry(olong disk, olong user, olong cno, ObitErr *err);

/** Public: Get last access */
void ObitAIPSDirGetAccess(ObitAIPSDirCatEntry* entry, gchar *timeDate);

/** Public: Change status of Catalog directory entry */
ObitAIPSDirStatusError
ObitAIPSDirStatus(olong disk, olong user, olong cno, 
		  ObitAIPSDirStatusCode code, ObitErr *err);

#endif /* OBITAIPSDIR_H */ 

