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
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "ObitFile.h"
#include "ObitAIPSDir.h"
#include "ObitAIPSCat.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitAIPSDir.c
 * ObitAIPSDir class function definitions.
 */

/*-----------------File Globals ---------------------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitAIPSDir";

/** ObitThread to use to lock access to the directory */
static ObitThread *myLock = NULL;

/*------------------ Structures -----------------------------*/
/** ObitAIPSDir Class Structure. */  
typedef struct {
  /** class name for verification */
  gchar *className;
  /** Open I/O channel */
  ObitFile *myFile;
  /** File name */
  gchar      *CatFile;
  /** disk number */
  olong disk;
  /** User id */
  olong user;
  /** maximum slot */
  olong maxcno;
  /** does the buffer need to be flushed? */
  gboolean flush;
} ObitAIPSDir;

/**
 * AIPS Catalog directory structure
 * This is at the beginning of a catalog.
 */
typedef struct {
  /** disk number (not used?) */
  olong disk;
  olong dummy;
  /** number of catalog entries */
  olong ncat;
  /** date created (y, m, d) */
  olong date_created[3];
  /** time created (h, m, s) */
  olong time_created[3];
  /** date of last access (y, m, d) */
  olong date_access[3];
  /** time of last access (h, m, s) */
  olong time_access[3];
} ObitAIPSDirCatHead;

/*---------------Private function prototypes----------------*/
/** Private: Open catalog directory. */
static ObitAIPSDir* 
ObitAIPSDirOpen (gint disk, olong user, ObitErr *err);

/** Private: Close catalog directory. */
static void 
ObitAIPSDirClose (ObitAIPSDir* in, ObitErr *err);

/** Private: Lookup AIPS file in directory. */
static olong 
ObitAIPSDirFindEntry (ObitAIPSDir* in, gchar Aname[13], 
		      gchar Aclass[7], gchar Atype[3], 
		      olong seq, ObitErr *err);

/** Private: Find first free entry in directory. */
static olong 
ObitAIPSDirFindFree (ObitAIPSDir* in, gchar Aname[13], 
		      gchar Aclass[7], gchar Atype[3], 
		      olong seq, ObitErr *err);

/** Private: Read a Catalog directory entry. */
static void
ObitAIPSDirRead(ObitAIPSDir* in, olong cno,
		ObitAIPSDirCatEntry *entry, ObitErr *err);

/** Private: Write a Catalog directory entry. */
static void
ObitAIPSDirWrite(ObitAIPSDir* in, olong cno,
		 ObitAIPSDirCatEntry *entry, ObitErr *err);

/** Private: Copy a Catalog directory entry. */
static void
ObitAIPSDirCopy(ObitAIPSDirCatEntry *in, ObitAIPSDirCatEntry *out);

/** Private: Add a record block to the catalog directory. */
static void
ObitAIPSDirExtend(ObitAIPSDir* in, ObitErr *err);

/** Private: Initialize directory header */
static void
ObitAIPSDirInitHead(ObitAIPSDirCatHead *head, olong disk);

/** Private: Update access time on directory header */
static void
ObitAIPSDirUpdateHead(ObitAIPSDirCatHead *head);

/** Private: Update access time on directory entry */
static void
ObitAIPSDirUpdateEntry(ObitAIPSDirCatEntry *entry);

/** Private: Returns true if input is a  ObitAIPSDir* */
gboolean ObitAIPSDirIsA (ObitAIPSDir* in);

/** Private: Convert packed time word to its parts */
static void ObitAIPSDirUnpackTime (AIPSint pack, olong unpack[3]);

/** Private: Convert time triplet to packed. */
static void ObitAIPSDirPackTime (AIPSint *pack, olong unpack[3]);

/** Private: Convert packed date word to its parts */
static void ObitAIPSDirUnpackDate (AIPSint pack, olong unpack[3]);

/** Private: Convert date triplet to packed. */
static void ObitAIPSDirPackDate (AIPSint *pack, olong unpack[3]);

/** Private: Class initializer. */
void ObitAIPSDirClassInit (gint number, gchar* dir[]);

/*---------------Public functions---------------------------*/
/**
 * Look through AIPS catalog on a given disk to find a given
 * name, class, type and sequence
 * \param disk disk number to search.
 * \param user AIPS user number
 * \param Aname  AIPS name.
 * \param Aclass AIPS class.
 * \param Atype  AIPS file type ("MA", "UV", "SC", "  "=>any).
 * \param seq   AIPS sequence number.
 * \param err   Obit error stack
 * \return the catalog slot number or -1 if it was not found.
 */
olong ObitAIPSDirFindCNO(gint disk, olong user, 
		     gchar Aname[13], gchar Aclass[7], 
		     gchar Atype[3], olong seq, ObitErr *err)
{
  ObitAIPSDir         *myDir = NULL;
  olong cno=-1, ndisk, i;
  ObitAIPSDirCatEntry entry;
  gchar lAname[13], lAclass[7], lAtype[3];
  gchar *routine = "ObitAIPSDirFindCNO";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return -1;  /* previous error? */

  /* protect against bad strings */
  for (i=0; i<12; i++) lAname[i]  = Aname[i];  lAname[i] = 0;
  for (i=0; i<6; i++)  lAclass[i] = Aclass[i]; lAclass[i] = 0;
  for (i=0; i<2; i++)  lAtype[i]  = Atype[i];  lAtype[i] = 0;
  
  /* Check that disk legal */
  ndisk = ObitAIPSGetNumDisk(err);
  if ((disk <= 0) || (disk > ndisk)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: illegal AIPS disk number %d not in [1, %d]", 
      routine, disk, ndisk);
    return -1;
  }

 /* Open */
  myDir =  ObitAIPSDirOpen (disk, user, err);
  if (err->error) return -1;

  /* Find entry */
  cno = ObitAIPSDirFindEntry (myDir, lAname, lAclass, lAtype, 
			      seq, err);
  if (err->error) { /* attempt close on error */
    ObitAIPSDirClose (myDir, err); 
    return -1;
  }

  /* Log error and return */
  if (cno<1) {
    Obit_log_error(err, OBIT_InfoWarn, 
		   "Could not find AIPS %s file %s %s seq %d disk %d", 
		   lAtype, lAname, lAclass, seq, disk);
    ObitAIPSDirClose (myDir, err); 
    return -1;
  }

  /* Update last access */
  ObitAIPSDirRead(myDir, cno, &entry, err);
  if (err->error) { /* attempt close on error */
    ObitAIPSDirClose (myDir, err); 
    return -1;
  }

  /* change time */
  ObitAIPSDirUpdateEntry(&entry);

  /* write it back */
  ObitAIPSDirWrite(myDir, cno, &entry, err);
  if (err->error) { /* attempt close on error */
    ObitAIPSDirClose (myDir, err); 
    return -1;
  }
  
  /* close file */
  ObitAIPSDirClose (myDir, err);

  return cno;
} /* end ObitAIPSDirFindCNO */

/**
 * If the given entry already exist in the catalog is is
 * returned and the exist value is set TRUE.
 * If it doesn't exist find the first free slot and allocate
 * it to the new entry.
 * \param disk   disk number.
 * \param user   AIPS user number.
 * \param Aname  AIPS name.
 * \param Aclass AIPS class.
 * \param Atype  AIPS file type (MA, UV, SC).
 * \param seq    AIPS sequence number.  If <=0 then one is assigned
 *               if a matching file exists, the highest sequence number, 
 *               is used, else seq 1 is allocated.
 * \param exist  [out] TRUE iff the entry previously existed.
 * \param err    Obit error stack.
 * \return the catalog slot number or -1 if it was not found.
 */
olong ObitAIPSDirAlloc(gint disk, olong user, 
		     gchar Aname[13], gchar Aclass[7], gchar Atype[3], 
		     olong seq, gboolean *exist, ObitErr *err)
{
  ObitAIPSDir         *myDir = NULL;
  olong cno = -1, ndisk, i;
  ObitAIPSDirCatEntry entry;
  gchar lAname[13], lAclass[7], lAtype[3];
  gchar *routine = "ObitAIPSDirAlloc";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return -1;  /* previous error? */

  /* protect against bad strings */
  for (i=0; i<12; i++) lAname[i]  = Aname[i];  lAname[i] = 0;
  for (i=0; i<6; i++)  lAclass[i] = Aclass[i]; lAclass[i] = 0;
  for (i=0; i<2; i++)  lAtype[i]  = Atype[i];  lAtype[i] = 0;
  
  /* Check that disk legal */
  ndisk = ObitAIPSGetNumDisk(err);
  if ((disk <= 0) || (disk > ndisk)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: illegal AIPS disk number %d not in [1, %d]", 
      routine, disk, ndisk);
    return -1;
  }

  /* Check/fix sequence number */
  if (seq<1) {
    seq = ObitAIPSDirHiSeq (disk, user, lAname, lAclass, lAtype, TRUE, err);
    if (err->error) Obit_traceback_val (err, routine, "Catalog search", cno);
  }

  /* Open */
  myDir =  ObitAIPSDirOpen (disk, user, err);
  if (err->error) return -1;

  /* is it already there */
  cno = ObitAIPSDirFindEntry (myDir, lAname, lAclass, lAtype, 
			      seq, err);
  if (err->error) { /* attempt close on error */
    ObitAIPSDirClose (myDir, err); 
    return -1;
  }

  *exist = cno>0; /* find it? */

  /* If not - find a free one */
  if (!*exist) {
    cno = ObitAIPSDirFindFree (myDir, lAname, lAclass, lAtype, 
			       seq, err);
    if (err->error) { /* attempt close on error */
      ObitAIPSDirClose (myDir, err); 
      return -1;
    }
  }

  /* Enter values - read */
  ObitAIPSDirRead(myDir, cno, &entry, err);
  if (err->error) { /* attempt close on error */
    ObitAIPSDirClose (myDir, err); 
    return -1;
  }

  /* register information */
  if (!*exist) { /* only if new entry */
    entry.user  = user;
    entry.seq   = seq;
    g_memmove(entry.name,  lAname, 12);
    g_memmove(entry.class, lAclass, 6);
    g_memmove(entry.type,  lAtype, 2);
  }

  /* access time time */
  ObitAIPSDirUpdateEntry(&entry);

  /* write it back */
  ObitAIPSDirWrite(myDir, cno, &entry, err);
  if (err->error) { /* attempt close on error */
    ObitAIPSDirClose (myDir, err); 
    return -1;
  }
  
  /* close file */
  ObitAIPSDirClose (myDir, err);
  if (err->error) Obit_traceback_val (err, routine, "Catalog search", cno);

  /* Write Dummy AIPS header */
  if ((cno>0) && (!*exist)) {
    ObitAIPSCatDummy(disk, user, Aname, Aclass, Atype, seq, cno, err);
    if (err->error) Obit_traceback_val (err, routine, "Dummy header", cno);
  }

  return cno;
} /* end ObitAIPSDirAlloc */

/**
 * Mark specified catalog entry as unoccupied.
 * \param disk   disk number.
 * \param user   AIPS user number.
 * \param cno    Catalog slot number.
 * \param err    Obit error stack.
 */
void ObitAIPSDirRemoveEntry(gint disk, olong user, olong cno, ObitErr *err)
{
  ObitAIPSDir         *myDir = NULL;
  ObitAIPSDirCatEntry entry;
  olong ndisk;
  ObitAIPSDirStatusError retCode = OBIT_AIPS_Dir_StatusSpecErr;
  gchar *routine = "ObitAIPSDirRemoveEntry";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */

  /* Check that disk legal */
  ndisk = ObitAIPSGetNumDisk(err);
  if ((disk <= 0) || (disk > ndisk)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: illegal AIPS disk number %d not in [1, %d]", 
      routine, disk, ndisk);
    return;
  }

  /* Open */
  retCode = OBIT_AIPS_Dir_StatusIOErr;
  myDir =  ObitAIPSDirOpen (disk, user, err);
  if (err->error) return;

  /* Read entry */
  ObitAIPSDirRead(myDir, cno, &entry, err);
  if (err->error) { /* attempt close on error */
    ObitAIPSDirClose (myDir, err); 
    return;
  }

  /* Must have no read/write status */
  if (entry.status!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "Cannot remove Catalog entry - active status");
    return;
  }

  /* Mark as unoccupied */
  entry.user = -1;

  /* access time time */
  ObitAIPSDirUpdateEntry(&entry);

  /* write it back */
  ObitAIPSDirWrite(myDir, cno, &entry, err);
  if (err->error) { /* attempt close on error */
    retCode = OBIT_AIPS_Dir_StatusIOErr;
    ObitAIPSDirClose (myDir, err); 
    return;
  }
  
  /* close file */
  ObitAIPSDirClose (myDir, err);

} /* end ObitAIPSDirRemoveEntry */

/**
 * Return maximum allocated slot number
 * \param disk   disk number.
 * \param user   AIPS user number.
 * \param err    Obit error stack.
 * \return highest catalog slot number in directory
 */
olong ObitAIPSDirNumber(gint disk, olong user, ObitErr *err)
{
  olong                ndisk, out = 0;
  ObitAIPSDir         *myDir = NULL;
  ObitAIPSDirStatusError retCode = OBIT_AIPS_Dir_StatusSpecErr;
  gchar *routine = "ObitAIPSDirNumber";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;  /* previous error? */

  /* Check that disk legal */
  ndisk = ObitAIPSGetNumDisk(err);
  if ((disk <= 0) || (disk > ndisk)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: illegal AIPS disk number %d not in [1, %d]", 
      routine, disk, ndisk);
    return out;
  }

  /* Open */
  retCode = OBIT_AIPS_Dir_StatusIOErr;
  myDir =  ObitAIPSDirOpen (disk, user, err);
  if (err->error) return out;
  out = myDir->maxcno;  /* Get maximum allocate */

  /* close file */
  ObitAIPSDirClose (myDir, err);

  return out;
} /* end ObitAIPSDirNumber */

/**
 * Searches directory(ies) for matching file name to determine the highest
 * sequence number.  
 * If exist==TRUE and a match is found the highest value is returned.  
 * If If exist==FALSE and a match is found , the highest seq+1 is returned
 * If there are no matches, 1 is returned
 * \param disk   disk number., if 0, check all
 * \param user   AIPS user number.
 * \param Aname  AIPS name, NULLs assumed equivalent to blank
 * \param Aclass AIPS class. NULLs assumed equivalent to blank
 * \param Atype  AIPS file type (MA, UV, SC).
 * \param exist  TRUE iff the entry already exists (open for read)
 * \param err    Obit error stack.
 * \return the desired AIPS sequence number, -1 on error
 */
olong ObitAIPSDirHiSeq(gint disk, olong user, gchar Aname[13], gchar Aclass[7], 
		      gchar Atype[3], gboolean exist, ObitErr *err)
{
  olong i, outSeq = -1;
  olong loDisk, hiDisk, iDisk;
  ObitAIPSDir   *myDir = NULL;
  olong cno = -1, ndisk;
  ObitAIPSDirCatEntry entry;
  gboolean found;
  gchar *routine = "ObitAIPSDirHiSeq";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outSeq;  /* previous error? */

  /* Clean input strings - make sure no nulls in characters to be
     compared */
  for (i=0; i<12; i++) if (Aname[i]==0)  Aname[i] = ' ';
  for (i=0; i<6 ; i++) if (Aclass[i]==0) Aclass[i] = ' ';
  for (i=0; i<2 ; i++) if (Atype[i]==0)  Atype[i] = ' ';

  /* Loop over disks */
  ndisk = ObitAIPSGetNumDisk(err);
  if (disk<=0) {
    loDisk = 1;
    hiDisk = ndisk;
  } else {
    loDisk = disk;
    hiDisk = disk;
  }

  for (iDisk = loDisk; iDisk<= hiDisk; iDisk++) {
    /* Open */
    myDir =  ObitAIPSDirOpen (disk, user, err);
    if (err->error) Obit_traceback_val (err, routine, "Catalog search", -1);

    /* Loop over slots */
    for (cno=1; cno<=myDir->maxcno; cno++) {
      ObitAIPSDirRead(myDir, cno, &entry, err);
      if (err->error) { /* attempt close on error */
	ObitAIPSDirClose (myDir, err); 
	Obit_traceback_val (err, routine, "Catalog search", outSeq);
      }

      /* Does this match? */
      found  = (user == entry.user) &&
	(!strncmp (Aname,  entry.name, 12)) &&
	(!strncmp (Aclass, entry.class, 6)) &&
	(!strncmp (Atype,  entry.type,  2));
      if (found) outSeq = MAX (outSeq, entry.seq);
    }

    /* close directory */
    ObitAIPSDirClose (myDir, err);
    if (err->error) Obit_traceback_val (err, routine, "Catalog search", outSeq);

  } /* end loop over disk */

  /* One found? */
  if (outSeq<=0) {    /* Nope - use 1 */
    outSeq = 1;
  } else if (exist) { /* Want existing one, use highest found */
    outSeq = outSeq;
  } else {            /* Want new one, use highest found + 1 */
    outSeq++;
  }

  return outSeq;
} /* end ObitAIPSDirHiSeq */

/**
 * Rename directory entry
 * \param disk     Disk number
 * \param user     AIPS user number.
 * \param cno      Slot number to be renamed
 * \param newName  New AIPS name (12 characters)
 * \param newClass New AIPS Class (6 characters)
 * \param newSeq   New AIPS sequence
 * \param err      Obit error stack.
 */
void ObitAIPSDirRename(gint disk, olong user,  olong cno, gchar *newName, 
			gchar newClass[7], olong newSeq, ObitErr *err)
{
  ObitAIPSDir   *myDir = NULL;
  ObitAIPSDirCatEntry entry;
  olong i;
  gchar Aname[13], Aclass[7];
  gchar *routine = "ObitAIPSDirHiSeq";
  
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  
  /* Clean input strings - make sure no nulls in characters to be
     compared */
  strncpy (Aname, newName,  12);
  for (i=0; i<12; i++) if (Aname[i]==0)  Aname[i] = ' ';
  strncpy (Aclass, newClass, 6);
  for (i=0; i<6 ; i++) if (Aclass[i]==0) Aclass[i] = ' ';
  
  /* Open directory */
  myDir =  ObitAIPSDirOpen (disk, user, err);
  if (err->error) Obit_traceback_msg (err, routine, "Rename");

  ObitAIPSDirRead(myDir, cno, &entry, err);
  if (err->error) { /* attempt close on error */
    ObitAIPSDirClose (myDir, err); 
    Obit_traceback_msg (err, routine, "Rename");
  }

  /* Change info */
  entry.seq   = newSeq;
  g_memmove(entry.name,  Aname, 12);
  g_memmove(entry.class, Aclass, 6);
  
  /* access time time */
  ObitAIPSDirUpdateEntry(&entry);
  
  /* write it back */
  ObitAIPSDirWrite(myDir, cno, &entry, err);
  if (err->error) { /* attempt close on error */
    ObitAIPSDirClose (myDir, err); 
    return;
  }
  
  /* close file */
  ObitAIPSDirClose (myDir, err);
  if (err->error) Obit_traceback_msg (err, routine, "Rename");

} /* end ObitAIPSDirRename */

/** 
 * Returns the entry for a given catalog slot.
 * \param disk disk number.
 * \param user user id number.
 * \param cno  catalog slot number
 * \param  err ObitErr error stack.
 * \return pointer to a newly created ObitAIPSDirCatEntry, 
 *         NULL on failure, this must be freed (g_free).
 */
ObitAIPSDirCatEntry* 
ObitAIPSDirGetEntry(gint disk, olong user, olong cno, ObitErr *err)
{
  ObitAIPSDir         *myDir = NULL;
  olong ndisk;
  ObitAIPSDirCatEntry *entry;
  gchar *routine = "ObitAIPSDirGetEntry";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;  /* previous error? */

  /* Check that disk legal */
  ndisk = ObitAIPSGetNumDisk(err);
  if ((disk <= 0) || (disk > ndisk)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: illegal AIPS disk number %d not in [1, %d]", 
      routine, disk, ndisk);
    return NULL;
  }

  /* Open */
  myDir =  ObitAIPSDirOpen (disk, user, err);
  if (err->error) return NULL;

  /* create output */
  entry = g_malloc0(sizeof(ObitAIPSDirCatEntry));

  /* Read entry */
  ObitAIPSDirRead(myDir, cno, entry, err);
  if (err->error) { /* attempt close on error */
    ObitAIPSDirClose (myDir, err); 
    g_free(entry);
    return NULL;
  }

  /* access time time
     ObitAIPSDirUpdateEntry(entry); */

  /* write it back */
  ObitAIPSDirWrite(myDir, cno, entry, err);
  if (err->error) { /* attempt close on error */
    ObitAIPSDirClose (myDir, err); 
    g_free(entry);
    return NULL;
  }
  
  /* close file */
  ObitAIPSDirClose (myDir, err);

  return entry;
} /* end ObitAIPSDirGetEntry */

/** 
 * Gets time/date string from an AIPS directory entry
 * in the form "xx-Mon-yyyy hh:mm:ss"
 * \param entry    Structure whose access time is desired
 * \param timeDate String to accept date/time string
 *                 Must have at least 21 char allocated.
 */
void ObitAIPSDirGetAccess(ObitAIPSDirCatEntry* entry, gchar *timeDate)
{
  olong edate[3], etime[3];
  gchar *months[12] = {"Jan","Feb","Mar","Apr","May","Jun","Jul",
		       "Aug","Sep","Oct","Nov","Dec"};

  /* Unpack */
  ObitAIPSDirUnpackDate (entry->access[0], edate);
  ObitAIPSDirUnpackTime (entry->access[1], etime);

  /* encode */
  g_snprintf (timeDate, 21, "%2.2d-%s-%4d %2.2d:%2.2d:%2.2d",
	      edate[2], months[edate[1]-1], edate[0], 
	      etime[0], etime[1], etime[2]);

} /* end ObitAIPSDirGetAccess */

/** 
 * Change the status of an AIPS catalog directory entry.
 * Only one program/thread can mark an entry write and if the
 * entry is marked write no reads can be added. 
 * Errors are indicated in err and as the return code.
 * \param disk disk number.
 * \param user user id number.
 * \param cno  catalog slot number
 * \param code a status change code as an #ObitAIPSDirStatusCode
 *             defined in ObitAIPSDir.h.
 *             \li OBIT_AIPS_Dir_AddWrite   => Add Write Status
 *             \li OBIT_AIPS_Dir_ClearWrite => Clear Write Status
 *             \li OBIT_AIPS_Dir_IncRead    => Increment Read Status
 *             \li OBIT_AIPS_Dir_DecRead    => Decrement Read Status
 * \param  err ObitErr error stack.
 * \return a completion code as an #ObitAIPSDirStatusError.
 *         OBIT_AIPS_Dir_StatusOK on success.
 */
ObitAIPSDirStatusError
ObitAIPSDirStatus(gint disk, olong user, olong cno, 
		  ObitAIPSDirStatusCode code, ObitErr *err)
{
  ObitAIPSDir         *myDir = NULL;
  ObitAIPSDirCatEntry entry;
  ObitAIPSDirStatusError retCode = OBIT_AIPS_Dir_StatusSpecErr;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return retCode;  /* previous error? */

  /* Open */
  retCode = OBIT_AIPS_Dir_StatusIOErr;
  myDir =  ObitAIPSDirOpen (disk, user, err);
  if (err->error) return retCode;

  /* Read entry */
  ObitAIPSDirRead(myDir, cno, &entry, err);
  if (err->error) { /* attempt close on error */
    ObitAIPSDirClose (myDir, err); 
    return retCode;
  }

  /* Check if operation allowed */
  retCode = OBIT_AIPS_Dir_StatusOK;
  switch (code) { 
  case OBIT_AIPS_Dir_AddWrite:
    /* has a write status? */
    if (entry.status<0) {
      retCode = OBIT_AIPS_Dir_StatusWrite;
      Obit_log_error(err, OBIT_Error, 
		     "Cannot mark WRITE, Catalog Status already WRITE");
    }
    /* has a read status? */
    if (entry.status>0) {
      retCode = OBIT_AIPS_Dir_StatusRead;
      Obit_log_error(err, OBIT_Error, 
		     "Cannot mark WRITE, Catalog Status already Read");
    }
    break;
  case OBIT_AIPS_Dir_ClearWrite:
    /* Must be marked write */
    if (entry.status>=0) {
      retCode = OBIT_AIPS_Dir_StatusSpecErr;
      Obit_log_error(err, OBIT_Error, 
		     "Cannot clear WRITE, Catalog Status not WRITE");
    }
    break;
  case OBIT_AIPS_Dir_IncRead:
    /* has a write status? */
    if (entry.status<0) {
      retCode = OBIT_AIPS_Dir_StatusWrite;
      Obit_log_error(err, OBIT_Error, 
		     "Cannot mark WRITE, Catalog Status already WRITE");
    }
    break;
  case OBIT_AIPS_Dir_DecRead:
    /* Must be marked read */
    if ((entry.status==0) || (entry.status==-1)) {
      retCode = OBIT_AIPS_Dir_StatusSpecErr;
      Obit_log_error(err, OBIT_Error, 
		     "Cannot clear READ, Catalog Status not READ");
    }
    break;
  case OBIT_AIPS_Dir_ClearAll:
    /* Always succeed */
    break;
  default:
    retCode = OBIT_AIPS_Dir_StatusSpecErr;
    Obit_log_error(err, OBIT_Error, 
		   "Unknown Catalog Status change code %d", 
		   code);
    g_assert_not_reached(); /* unknown, barf */
  }; /* end switch to check if allowed */

  /* Do operation if allowed */
  if (retCode == OBIT_AIPS_Dir_StatusOK) {
    switch (code) { 
    case OBIT_AIPS_Dir_AddWrite:
      if (entry.status==0) /* no previous status */
	entry.status = -1;
      else if (entry.status > 0) /* previous write+read */
	entry.status = -(1+entry.status);
      break;
    case OBIT_AIPS_Dir_ClearWrite:
      if (entry.status==0) /* no previous status */
	/* Shouldn't happen (should have been caught earlier) */
	retCode = OBIT_AIPS_Dir_StatusSpecErr;
      else if (entry.status==-1) /* previous write */
	entry.status = 0;
      else if (entry.status < -1) /* previous write+read */
	entry.status = -(1+entry.status);
      break;
    case OBIT_AIPS_Dir_IncRead:
      if (entry.status>=0) /* no or only read previous status */
	entry.status++;
      else if (entry.status < 0) /* previous write+read? */
	entry.status--;
      break;
    case OBIT_AIPS_Dir_DecRead:
      if (entry.status==0) /* no previous status */
	/* Shouldn't happen */
	retCode = OBIT_AIPS_Dir_StatusSpecErr;
      else if (entry.status>0) /* previous read */
	entry.status--;
      else if (entry.status < -1) /* previous write+read */
	entry.status++;
      else if (entry.status == -1) /* previous write */
	/* should not occur - only marked as write */
	retCode = OBIT_AIPS_Dir_StatusWrite;
      break;
    case OBIT_AIPS_Dir_ClearAll:
      entry.status = 0;
      break;
    default:
      g_assert_not_reached(); /* unknown, barf */
    }; /* end switch */
  } /* end of do operation */

  /* access time time */
  ObitAIPSDirUpdateEntry(&entry);

  /* write it back */
  ObitAIPSDirWrite(myDir, cno, &entry, err);
  if (err->error) { /* attempt close on error */
    retCode = OBIT_AIPS_Dir_StatusIOErr;
    ObitAIPSDirClose (myDir, err); 
    return retCode;
  }
  
  /* close file */
  ObitAIPSDirClose (myDir, err);

  return retCode;
} /* end ObitAIPSDirStatus */

/*---------------Private functions---------------------------*/
/** 
 * Open catalog directory and lock against other thread access.
 * Updates last access time.
 * \param disk disk number.
 * \param user AIPS user number.
 * \param  err ObitErr error stack.
 * \return pointer to a newly created ObitAIPSDir, NULL on failure,
 *         this must be freed by a call to ObitAIPSDirClose.
 */
static ObitAIPSDir* 
ObitAIPSDirOpen (gint disk, olong user, ObitErr *err)
{
  ObitAIPSDir         *out = NULL;
  ObitAIPSDirCatHead  *head=NULL;
  AIPSint    buffer[256];
  ObitIOCode status;
  olong       ndisk;
  gboolean   exist;
  olong      size, wantPos;
  gchar *routine = "ObitAIPSDirOpen";

  /* create structure */
  out = g_malloc0(sizeof(ObitAIPSDir));
  out->className = myClassName;  /* class name pointer */
  out->disk      = disk;
  out->user      = user;
  out->maxcno    = 0;
  out->flush     = FALSE;
  out->CatFile   = NULL;
  out->myFile    = NULL;

  /* Check that disk legal */
  ndisk = ObitAIPSGetNumDisk(err);
  if ((disk <= 0) || (disk > ndisk)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: illegal AIPS disk number %d not in [1, %d]", 
      routine, disk, ndisk);
    return out;
  }

  /* Lock directory aginst other threads */
  if (myLock==NULL) myLock = newObitThread(); /* create lock first time */
  ObitThreadLock(myLock);

  /* Set file name */
  out->CatFile = 
    ObitAIPSFilename (OBIT_AIPS_Catalog, disk, 0, user, NULL, 0, err);
  if (err->error) {
    ObitThreadUnlock(myLock);
    Obit_traceback_val (err, routine, "Catalog search", NULL);
  }

  /* Does it currently exist? */
  if (out->myFile) ObitFileUnref(out->myFile);
  out->myFile = newObitFile("Catalog search");
  exist = ObitFileExist (out->CatFile, err);
  if (err->error) {
    ObitThreadUnlock(myLock);
    Obit_traceback_val (err, routine, "Catalog search", NULL);
  }

  /* open */
  size = 256 * sizeof(AIPSint);
  if (ObitFileOpen (out->myFile, out->CatFile, OBIT_IO_ReadWrite, 
		    OBIT_IO_Binary, size, err) || (err->error)) {
     Obit_log_error(err, OBIT_Error, 
		   "ERROR opening AIPS catalog file disk %d", out->disk);
    g_free(out->CatFile); /* going up in flames - clean up */
    g_free(out);  
    ObitThreadUnlock(myLock);
    Obit_traceback_val (err, routine, "Catalog search", NULL);
  }

  /* If it didn't previously exist - initialize buffer */
  if (!exist) {
    head = (ObitAIPSDirCatHead*)buffer;
    ObitAIPSDirInitHead (head, disk);  
    if (err->error) Obit_traceback_val (err, routine, "Catalog init", NULL);
  } else { /* exists - read header block */
    wantPos = 0;
    status = ObitFileRead (out->myFile, wantPos, size, (gchar*)buffer, err);
    if ((status!=OBIT_IO_OK) || (err->error)) {/* add traceback on error */
      Obit_log_error(err, OBIT_Error, 
		     "Status %d reading AIPS catalog file disk %d", status, out->disk);
      ObitThreadUnlock(myLock);
      Obit_traceback_val (err, routine, "Catalog search", NULL);
    }
  } /* end init/read header */

  /* update last access */
  /* Update header */
  head = (ObitAIPSDirCatHead*)buffer;
  ObitAIPSDirUpdateHead(head);

  /* How many entries */
  out->maxcno = head->ncat;

  /* rewrite */
  /* position file to beginning */
  wantPos = 0;

  /* write */
  status = ObitFileWrite (out->myFile, wantPos, size, (gchar*)buffer, err);
  if ((status!=OBIT_IO_OK) || (err->error)) {/* add traceback on error */
    ObitThreadUnlock(myLock);
    Obit_traceback_val (err, routine, "Catalog search", NULL);
  }

  /* If it didn't previously exist - add a block */
  if (!exist) ObitAIPSDirExtend (out, err);
  if (err->error)  {
    ObitThreadUnlock(myLock);
    Obit_traceback_val (err, routine, "Catalog extend", out);
  }

  out->flush = TRUE; /* something in buffer to write */
  return out;
} /* end ObitAIPSDirOpen  */

/**
 * Close catalog directory and release all resources.
 * \param in Pointer to catalog directory structure.
 *           This will be deallocated.
 * \param  err ObitErr error stack.
 */
static void 
ObitAIPSDirClose (ObitAIPSDir* in, ObitErr *err)
{
   ObitIOCode status;
   gchar *routine = "ObitAIPSDirClose";

  /* still unlock if in==NULL */
  if (in!=NULL) {
    /* close file */
    status = ObitFileClose (in->myFile, err);
    if ((status!=OBIT_IO_OK) || (err->error)) { /* add traceback on error */
      ObitThreadUnlock(myLock);
      Obit_traceback_msg (err, routine, "Catalog search");
    }

    /* delete */
    in->myFile = ObitFileUnref(in->myFile);

    if (in->CatFile) g_free (in->CatFile);

  } /* end of shutdown when in defined */

  /* Deallocate directory object */
  if (in) g_free (in);

  /* Unlock directory against other threads */
  ObitThreadUnlock(myLock);
} /* end ObitAIPSDirClose */

/**
 * Find a given entry in the open catalog.
 * \param  in Pointer to catalog directory structure info.
 * \param Aname  AIPS name.
 * \param Aclass AIPS class.
 * \param Atype  AIPS file type ("MA", "UV", "SC", "  "=>any).
 * \param seq    AIPS sequence number.
 * \param err    Obit error stack
 * \return the catalog slot number or -1 if it was not found.
 */
static olong 
ObitAIPSDirFindEntry (ObitAIPSDir* in, gchar Aname[13], 
		      gchar Aclass[7], gchar Atype[3], 
		      olong seq, ObitErr *err)
{
  olong       cno = 0;
  AIPSint    buffer[256];
  olong       i, nwpl, nlpr;
  ObitIOCode status;
  ObitFilePos size;
  gboolean   anyType, found = FALSE;
  olong      wantPos;
  ObitAIPSDirCatEntry *entry=NULL;
  gchar *routine = "ObitAIPSDirFindEntry";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return -1;  /* previous error? */
  g_assert (ObitAIPSDirIsA(in));

  /* Clean input strings - make sure no nulls in characters to be
     compared */
  for (i=0; i<12; i++) if (Aname[i]==0)  Aname[i] = ' ';
  for (i=0; i<6 ; i++) if (Aclass[i]==0) Aclass[i] = ' ';
  for (i=0; i<2 ; i++) if (Atype[i]==0)  Atype[i] = ' ';

  /* position file to second block - first is header */
  wantPos = 256 * sizeof(AIPSint);
  size    = 256 * sizeof(AIPSint);

  nwpl = 10;          /* number of AIPSint words per entry */
  nlpr = 256 / nwpl;  /* number of entries per "record" (256 words) */
  /* Blank type matches any */
  anyType = (Atype[0] == ' ') && (Atype[1] == ' ');
  /* Loop through file til found or EOF */
  while (!found) {

    /* read next block */
    status = ObitFileRead (in->myFile, wantPos, size, (gchar*)buffer, err);
    if (((status!=OBIT_IO_OK) && (status!=OBIT_IO_EOF)) || 
	(err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, "Catalog search", -1);
    wantPos = -1L; /* now sequential access */

    /* EOF? - quit loop  */
    if (status==OBIT_IO_EOF) break;

    entry = (ObitAIPSDirCatEntry*)buffer;  /* entry pointer into buffer */
    /* look through this block */
    for (i=0; i<nlpr; i++) {
      cno++;   /* count how many have gone by */

      /* check if this one matches */
      found  = (seq == entry->seq) &&
	(in->user == entry->user) &&
	(!strncmp (Aname,  entry->name, 12)) &&
	(!strncmp (Aclass, entry->class, 6)) &&
	((!strncmp (Atype,  entry->type,  2) || anyType));
      if (found) break; 

      entry++; /* advance pointer to next entry */
    }
  
  } /*  end loop through file */

  /* Was it found? */
  if (!found) cno = -1;

  return cno;
} /* end ObitAIPSDirFindEntry */

/**
 * Search through an open catalog and find first free slot.
 * Catalog directory will be extended if needed.
 * \param  in  Pointer to catalog directory structure info.
 * \param Aname  AIPS name.
 * \param Aclass AIPS class.
 * \param Atype  AIPS file type (MA, UV, SC).
 * \param seq    AIPS sequence number.
 * \param  err ObitErr error stack.
 * \return the catalog slot number.
 */
static olong ObitAIPSDirFindFree (ObitAIPSDir* in, gchar Aname[13], 
		      gchar Aclass[7], gchar Atype[3], 
		      olong seq, ObitErr *err)
{
  olong       cno = 0;
  AIPSint    buffer[256];
  olong       i, nwpl, nlpr;
  ObitIOCode status;
  ObitFilePos   wantPos, size;
  gboolean   found = FALSE;
  ObitAIPSDirCatEntry *entry=NULL;
  gchar *routine = "ObitAIPSDirFindFree";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return -1;  /* previous error? */
  g_assert (ObitAIPSDirIsA(in));

  /* position file to second block - first is header */
  wantPos = 256 * sizeof(AIPSint);
  size    = 256 * sizeof(AIPSint);

  nwpl = 10;          /* number of AIPSint words per entry */
  nlpr = 256 / nwpl;  /* number of entries per "record" (256 words) */
  /* A free entry has userid < 0 */
  /* Loop through file til found or EOF */
  while (!found) {

    /* read next block */
    status = ObitFileRead (in->myFile, wantPos, size, (gchar*)buffer, err);
    if (status==OBIT_IO_EOF) break;  /* EOF? - quit loop  */

    if ((status!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, "Catalog search", -1);
    wantPos = -1L; /* now sequential access */
    
    entry = (ObitAIPSDirCatEntry*)buffer;  /* entry pointer into buffer */
    /* look through this block */
    for (i=0; i<nlpr; i++) {
      cno++;   /* count how many have gone by */

      /* check if this one matches */
      found  =  (entry->user < 0);
      if (found) break; 
      entry++; /* advance pointer to next entry */
    }
  } /*  end loop through file */

  /* Was it found? */
  if (!found) {
    /* add another block and use first one */
    ObitAIPSDirExtend(in, err);
    if (err->error) return -1; /* it work? */
    cno++;
  }

  return cno;
} /* end ObitAIPSDirFindFree */

/**
 * Read the contents of the specified entry.
 * \param in    Pointer to catalog directory structure info.
 * \param cno   which slot number to read.
 * \param entry Catalog entry structure to accept data.
 * \param err   Obit error stack.
 */
static void
ObitAIPSDirRead(ObitAIPSDir* in, olong cno, 
		ObitAIPSDirCatEntry *entry, ObitErr *err)
{
  AIPSint    buffer[20];
  olong       nwpl, nlpr, ib, ir;
  ObitIOCode status;
  olong      size, wantPos;
  ObitAIPSDirCatEntry *fentry=NULL;
  gchar *routine = "ObitAIPSDirRead";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  g_assert (ObitAIPSDirIsA(in));
  g_assert(entry != NULL);

  /* Check CNO range */
  if ((cno<1) || (cno>in->maxcno)) {
    Obit_log_error(err, OBIT_Error, 
		   "Catalog slot %d out of range (1- %d) disk %d",
		   cno,in->maxcno,in->disk);
    return;
  }


  nwpl = 10;          /* number of AIPSint words per entry */
  nlpr = 256 / nwpl;  /* number of entries per "record" (256 words) */

  /* where is this entry */
  ib = 1 + ((cno - 1) / nlpr); /* block number */
  ir = 1 + ((cno-1) % nlpr);   /* record number */
  wantPos = sizeof(AIPSint) * ((ib * 256) + ((ir-1) * nwpl));
  size = nwpl * sizeof(AIPSint);

  /* read record */
  status = ObitFileRead (in->myFile, wantPos, size, (gchar*)buffer, err);
  if ((status!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, "Catalog search");

  /* Copy entry from buffer */
  fentry = (ObitAIPSDirCatEntry*)buffer;  /* entry pointer into buffer */
  ObitAIPSDirCopy (fentry, entry);

} /* end ObitAIPSDirRead */

/**
 * Write the contents of the specified entry.
 * \param in    Pointer to catalog directory structure info.
 * \param cno   which slot number to write.
 * \param entry Catalog entry structure to accept data.
 * \param err   Obit error stack.
 */
static void
ObitAIPSDirWrite(ObitAIPSDir* in, olong cno,
		 ObitAIPSDirCatEntry *entry, ObitErr *err)
{
  AIPSint    buffer[20];
  olong       nwpl, nlpr, ib, ir;
  ObitIOCode status;
  ObitFilePos size, wantPos;
  ObitAIPSDirCatEntry *fentry=NULL;
  gchar *routine = "ObitAIPSDirWrite";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  g_assert (ObitAIPSDirIsA(in));
  g_assert(entry != NULL);

  /* Check CNO range */
  if ((cno<1) || (cno>in->maxcno)) {
    Obit_log_error(err, OBIT_Error, 
		   "Catalog slot %d out of range (1- %d) disk %d", 
		   cno, in->maxcno, in->disk);
    return;
  }


  nwpl = 10;          /* number of AIPSint words per entry */
  nlpr = 256 / nwpl;  /* number of entries per "record" (256 words) */

  /* where is this entry */
  ib = 1 + ((cno - 1) / nlpr); /* block number */
  ir = 1 + ((cno-1) % nlpr);   /* record number */
  wantPos = sizeof(AIPSint) * ((ib * 256) + ((ir-1) * nwpl));
  size = nwpl * sizeof(AIPSint);

  /* Copy entry to buffer */
  fentry = (ObitAIPSDirCatEntry*)buffer;  /* entry pointer into buffer */
  ObitAIPSDirCopy (entry, fentry);

  /* write record */
  status = ObitFileWrite (in->myFile, wantPos, size, (gchar*)buffer, err);
  if ((status!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, "Catalog search");

  in->flush = TRUE; /* something in buffer to flush */
} /* end ObitAIPSDirWrite */

/**
 * Copy the contents from in to out.
 * \param in   Input Catalog entry structure.
 * \param out  Output Catalog entry structure to accept data.
 */
static void
ObitAIPSDirCopy(ObitAIPSDirCatEntry *in, ObitAIPSDirCatEntry *out)
{
  olong i;

  /* Error checks */
  g_assert(in  != NULL);
  g_assert(out != NULL);

  out->user      = in->user;
  out->status    = in->status;
  out->access[0] = in->access[0];
  out->access[1] = in->access[1];
  out->seq       = in->seq;
  g_memmove (out->name, in-> name, 20);

  /* replace any NULLs in the string with blanks */
  for (i=0; i<20; i++) 
    if (out->name[i]==0) out->name[i]=' ';

} /* end ObitAIPSDirCopy */

/**
 * Add a record block to the catalog directory and initialize.
 * NOTE: this should only be called with the file is positioned at the end
 * \param in  Pointer to catalog directory structure info.
 * \param err Obit error stack.
 */
static void
ObitAIPSDirExtend(ObitAIPSDir* in, ObitErr *err)
{
  AIPSint    buffer[256];
  olong       i, nwpl, nlpr;
  ObitIOCode status;
  olong      size, wantPos;
  gchar      blank[21] = "                    ";
  ObitAIPSDirCatEntry *entry=NULL;
  ObitAIPSDirCatHead  *head=NULL;
  gchar *routine = "ObitAIPSDirExtend";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  g_assert (ObitAIPSDirIsA(in));

  /* position file to end */
  /* This appears not to work and is extremely risky c IO really sucks 
     status = ObitFileEnd (in->myFile, err);*/

  /* initialize buffer */
  nwpl = 10;          /* number of AIPSint words per entry */
  nlpr = 256 / nwpl;  /* number of entries per "record" (256 words) */
  entry = (ObitAIPSDirCatEntry*)buffer;  /* entry pointer into buffer */
  for (i=0; i<nlpr; i++) {
    /* initialize record */
    entry->user      = -1;
    entry->status    = OBIT_IO_OK;
    entry->access[0] = 0;
    entry->access[1] = 0;
    entry->seq       = 0;
    g_memmove (entry->name, blank, 20);
    entry++; /* next */
  }

  size = 256 * sizeof(AIPSint);
  wantPos = -1L;
  status = ObitFileWrite (in->myFile, wantPos, size, (gchar*)buffer, err);
  if ((status!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, "Extend AIPS catalog");

  /* Flush it to disk */
  status = ObitFileFlush (in->myFile, err);
  if ((status!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, "Extend AIPS catalog");

  /* update header */
  /* position file to beginning */
  wantPos = 0;

  /* read */
  status = ObitFileRead (in->myFile, wantPos, size, (gchar*)buffer, err);
  if ((status!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, "Extend AIPS catalog");

  /* update */
  in->maxcno +=nlpr;
  head = (ObitAIPSDirCatHead*)buffer;
  head->ncat = in->maxcno;

  /* rewrite */
  status = ObitFileWrite (in->myFile, wantPos, size, (gchar*)buffer, err);
  if ((status!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, "Extend AIPS catalog");

  /* Flush it to disk */
  status = ObitFileFlush (in->myFile, err);
  if ((status!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, "Extend AIPS catalog");

  in->flush = TRUE; /* something in buffer to flush */
} /* end ObitAIPSDirExtend */

/**
 * Update access time on ObitAIPSDirCatHead structure.
 * \param head Header record to update.
 * \param disk AIPS "disk" number
 */
static void
ObitAIPSDirInitHead(ObitAIPSDirCatHead *head, olong disk)
{
  struct tm *lp;
  time_t clock;

  /* Header info */
  head->disk  = disk;  /* disk number */
  head->dummy = 0;
  head->ncat  = 0;     /* Number of entries */

  /* Init creation time */
  /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);

  /* Convert to  broken-down time. */
  lp = localtime (&clock);

  /* to output */
  head->date_access[0] = lp->tm_year;
  if (head->date_created[0]<1000)  head->date_access[0] += 1900; /* full year */
  head->date_created[1] = lp->tm_mon+1; /* For some bizzare reason, month is 0-rel */
  head->date_created[2] = lp->tm_mday;
  head->time_created[0] = lp->tm_hour;
  head->time_created[1] = lp->tm_min;
  head->time_created[2] = lp->tm_sec;

  /* Init access time */
  ObitAIPSDirUpdateHead (head);
} /* end ObitAIPSDirInitHead */

/**
 * Update access time on ObitAIPSDirCatHead structure.
 * \param head Header record to update.
 */
static void
ObitAIPSDirUpdateHead(ObitAIPSDirCatHead *head)
{
  struct tm *lp;
  time_t clock;

  /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);

  /* Convert to  broken-down time. */
  lp = localtime (&clock);

  /* to output */
  head->date_access[0] = lp->tm_year;
  if (head->date_access[0]<1000)  head->date_access[0] += 1900; /* full year */
  head->date_access[1] = lp->tm_mon+1; /* For some bizzare reason, month is 0-rel */
  head->date_access[2] = lp->tm_mday;
  head->time_access[0] = lp->tm_hour;
  head->time_access[1] = lp->tm_min;
  head->time_access[2] = lp->tm_sec;
} /* end ObitAIPSDirUpdateHead */

/**
 * Update access time on ObitAIPSDirCatEntry structure.
 * \param entry ObitAIPSDirCatEntry structure to update
 */
static void
ObitAIPSDirUpdateEntry(ObitAIPSDirCatEntry *entry)
{
  struct tm *lp;
  time_t clock;
  olong timea[3], datea[3];

  /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);

  /* Convert to  broken-down time. */
  lp = localtime (&clock);

  /* to local arrays */
  datea[0] = lp->tm_year;
  if (datea[0]<1000) datea[0] += 1900; /* full year */
  datea[1] = lp->tm_mon+1; /* For some bizzare reason, month is 0-rel */
  datea[2] = lp->tm_mday;
  timea[0] = lp->tm_hour;
  timea[1] = lp->tm_min;
  timea[2] = lp->tm_sec;

  /* update output structure */
  ObitAIPSDirPackDate (&entry->access[0], datea);
  ObitAIPSDirPackTime (&entry->access[1], timea);
} /* end ObitAIPSDirUpdateEntry */

/**
 * Check class  name string at beginning of structure.
 * \param  in Pointer to catalog directory structure info.
 * \return TRUE if the correct class else FALSE.
 */
gboolean ObitAIPSDirIsA (ObitAIPSDir* in)
{
  gboolean out;

  /* error checks */
  if (in == NULL) return FALSE;
  if (in->className == NULL) return FALSE;

  /* compare class name member */
  out = !strcmp(in->className, myClassName);

  return out;
} /* end ObitAIPSDirIsA */

/**
 * Unpack time triplet (h, m, s).
 * \param pack packed version of triplet
 * \param unpack unpacked version of triplet.
 */
void ObitAIPSDirUnpackTime (AIPSint pack, olong unpack[3])
{
  unpack[2] = pack % 256;
  pack /= 256;
  unpack[1] = pack % 256;
  unpack[0] = pack / 256;
} /* end ObitAIPSDirUnpackTime */

/**
 * Pack time triplet (h, m, s).
 * \param pack packed version of triplet
 * \param unpack unpacked version of triplet.
 */
void ObitAIPSDirPackTime (AIPSint *pack, olong unpack[3])
{
  olong i1, i2, i3;

  i1 = unpack[0]; i2 = unpack[1]; i3 = unpack[2];
  /* stuff 'em together */
  *pack = 256 * (256 * i1 + i2) + i3;
} /* end ObitAIPSDirPackTime */

/**
 * Unpack date triplet (y, m, d).
 * \param pack packed version of triplet
 * \param unpack unpacked version of triplet.
 */
void ObitAIPSDirUnpackDate (AIPSint pack, olong unpack[3])
{
  unpack[2] = pack % 256;
  pack /= 256;
  unpack[1] = pack % 256;
  unpack[0] = 1900 + pack / 256;
} /* end ObitAIPSDirUnpackDate */

/**
 * Pack date triplet (y, m, d).
 * \param pack packed version of triplet
 * \param unpack unpacked version of triplet.
 */
void ObitAIPSDirPackDate (AIPSint *pack, olong unpack[3])
{
  olong i1, i2, i3;

  i1 = unpack[0]; i2 = unpack[1]; i3 = unpack[2];
  /* packed version years since 1900 */
  if (i1 > 1000) i1 -= 1900;
  if (i2<1) i2 = 1;
  /* stuff 'em together */
  *pack = 256 * (256 * i1 + i2) + i3;
} /* end ObitAIPSDirPackDate */



