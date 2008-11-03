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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include "Obit.h"
#include "ObitFITS.h"
#include "ObitHistory.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFITS.c
 * ObitFITS class function definitions.
 */

/*-----------------File Globals ---------------------------*/
/** Class information structure */
static ObitFITS ObitFITSInfo = {"ObitFITS",FALSE};

/** Pointer to Class information structure */
static ObitFITS *myFITSInfo = &ObitFITSInfo;

/*------------------ Structures -----------------------------*/

/*---------------Private function prototypes----------------*/
/*---------------Public functions---------------------------*/
/**
 * Save names of the FITS directories.
 * \param number of disks defined [0,MAXFITSDISK], >=0 -> none
 * \param dir the names of the directories
 *    If NULL then look for environment variables FITS, FITS01, FITS02...
 */
void ObitFITSClassInit (gint number, gchar* dir[])
{
  olong i;
  gchar fitsxx[8], *ev;
  gchar *da[]={"01","02","03","04","05","06","07","08","09","10",
	       "11","12","13","14","15","16","17","18","19","20"};


  myFITSInfo->NumberDisks = 0;
  for (i=0; i<MAXFITSDISK; i++) {
    myFITSInfo->FITSdir[i] = NULL;
  }

  /* now initialized */
  myFITSInfo->initialized = TRUE;

  /* Are directories given or should I look for them? */
  if (dir==NULL) { /* Look in environment */
    myFITSInfo->NumberDisks = 0;

    /* First try $FITS */
    ev = getenv ("FITS");
    if (ev) {
      myFITSInfo->FITSdir[myFITSInfo->NumberDisks] =
	g_strconcat(ev, "/", NULL);
      myFITSInfo->NumberDisks++;
    }

    /* Look for FITS01, FITS02... */
    for (i=0; i<MAXFITSDISK; i++) {
      g_snprintf (fitsxx,7, "FITS%s",da[i]);
      ev = getenv (fitsxx);
      if (ev) {
	myFITSInfo->FITSdir[myFITSInfo->NumberDisks] =
	  g_strconcat(ev, "/", NULL);
	myFITSInfo->NumberDisks++;
      } else {
	break;
      }
    }
  } else { /* use what are given */
    
    /* error checks */
    g_assert (number<=MAXFITSDISK);
    
    if (number<=0) {
      myFITSInfo->NumberDisks = 0;
      myFITSInfo->FITSdir[0] = NULL;
      return;
    }
    
    /* save directory names */
    myFITSInfo->NumberDisks = number;
    for (i=0; i<number; i++) 
      myFITSInfo->FITSdir[i] = g_strconcat(dir[i], "/", NULL);
  } /* end of initialize data directories */
    
} /* end ObitFITSClassInit */

/**
 * Frees directory strings
 */
void ObitFITSShutdown (void)
{
  olong i;

  myFITSInfo->NumberDisks = 0;
  for (i=0; i<MAXFITSDISK; i++) {
    if (myFITSInfo->FITSdir[i]) g_free(myFITSInfo->FITSdir[i]);
    myFITSInfo->FITSdir[i] = NULL;
  }

  /* now uninitialized */
  myFITSInfo->initialized = FALSE;

} /*  end ObitFITSShutdown */

/**
 * Add a directory to the list of directories for FITS files
 * Limit of MAXFITSDISK (20) total disks 
 * #ObitFITSClassInit must have been used to initialize.
 * \param dir   names of the directories with terminal '/'
 * \param err   Error stack for any error messages.
 * \return new 1-rel disk number, -1 on failure
 */
olong ObitFITSAddDir (gchar* dir, ObitErr *err)
{
  olong out = -1;

  if (err->error) return out;
  if (myFITSInfo->NumberDisks>=MAXFITSDISK) {
    /* too many */
    Obit_log_error(err, OBIT_Error, "FITS directory list FULL");
    return out;
  }

  /* add to list */
  myFITSInfo->FITSdir[myFITSInfo->NumberDisks] =
    g_strconcat(dir, "/", NULL);
  out = ++myFITSInfo->NumberDisks;
  return out;
} /* end ObitFITSAddDir */

/**
 * Replace directory path
 * Limit of MAXFITSDISK (20) total disks 
 * #ObitFITSClassInit must have been used to initialize.
 * \param dir   name of the directory 
 * \param disk      FITS "disk" number. 1-rel, =0 => ignore directory
 * \param err   Error stack for any error messages.
 * \return new 1-rel disk number, -1 on failure
 */
void ObitFITSSetDir (gchar* dir, gint disk, ObitErr *err)
{

  if (err->error) return;
  if (myFITSInfo->NumberDisks<disk) {
    /* Must already be there */
    Obit_log_error(err, OBIT_Error, "FITS directory %d not yet defined", disk);
    return;
  }

  /* add to list */
  if (myFITSInfo->FITSdir[disk]) g_free(myFITSInfo->FITSdir[disk]);
  myFITSInfo->FITSdir[disk] = g_strconcat(dir, "/", NULL);
} /* end ObitFITSSetDir */

/**
 * Forms file name from the various parts.
 * #ObitFITSClassInit must have been used to initialize.
 * \param disk      FITS "disk" number. 1-rel, =0 => ignore directory
 * \param fileName  name
 * \param err       Error stack for any error messages.
 * \return full path name string, should be deallocated when done
 */
gchar* 
ObitFITSFilename (gint disk, gchar* fileName, ObitErr *err)
{
  gchar *out;
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  /* if disk <0 just return fileName */
  if (disk<=0) {
    out = g_strdup (fileName);
    ObitTrimTrail(out);  /* Trim any trailing blanks */
    return out;
  }
  if (!myFITSInfo->initialized) /* FITS directories uninitialized */
    Obit_log_error(err, OBIT_Error, 
		   "FITS directories uninitialized");
  if ((disk<0) || (disk>myFITSInfo->NumberDisks)) /* Disk number out of range */
    Obit_log_error(err, OBIT_Error, 
		   "FITS disk number %d out of range [%d,%d]", 
		   disk, 1, myFITSInfo->NumberDisks);
  if (err->error) return NULL;


  /* if fileName begins with '!', put it at the beginning */
  /* put it all together */
  if (fileName[0]=='!' )
    out = g_strconcat ("!", ObitFITSDirname(disk, err), &fileName[1], NULL);
  else
    out = g_strconcat (ObitFITSDirname(disk, err), fileName, NULL);
  ObitTrimTrail(out);  /* Trim any trailing blanks */
  
  return out;
} /* end ObitFITSFilename */

/**
 * Returns pointer to directory string by FITS disk.
 * \param disk FITS disk number.
 * \param err  Error stack for any error messages.
 * \return directory name string, this is a pointer into a global 
 *         class structure and should not be g_freeed.
 */
gchar* ObitFITSDirname (gint disk, ObitErr *err)
{

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  if (!myFITSInfo->initialized) /* FITS directories uninitialized */
    Obit_log_error(err, OBIT_Error, 
		   "FITS directories uninitialized");
  if ((disk<0) || (disk>myFITSInfo->NumberDisks)) /* Disk number out of range */
    Obit_log_error(err, OBIT_Error, 
		   "FITS disk number %d out of range [%d,%d]", 
		   disk, 1, myFITSInfo->NumberDisks);
  if (myFITSInfo->FITSdir[disk-1]==NULL) /* Directory not defined */
    Obit_log_error(err, OBIT_Error, 
		   "FITS directory %d not defined", disk);
  if (err->error) return NULL;

  return myFITSInfo->FITSdir[disk-1];
} /* ObitFITSDirname  */


/**
 * Assigns scratch file naming information and writes to the info.
 * Makes name "pgmName+pgmNumber+'Scr'+scrNo"
 * \param pgmName    Program name
 * \param pgmNumber  Program incarnation number
 * \param disk       FITS disk number.
 * \param scrNo      Which scratch file number
 * \param info       ObitInfoList to write to
 * \param err        Error stack for any error messages.
 */
/** Public: Assign a scratch file info */
void ObitFITSAssign(gchar *pgmName, olong pgmNumber, 
		    olong disk, olong scrNo, ObitInfoList *info, 
		    ObitErr *err)
{
  gchar name[121];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOType ft;

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(info));

  /* form name string */
  g_snprintf (name, 120, "%s%dScr%d", pgmName, pgmNumber, scrNo);

  /* write assignment to info */
  ft = OBIT_IO_FITS;
  ObitInfoListPut (info, "FileType", OBIT_long, dim, (gpointer)&ft,   err);
  ObitInfoListPut (info, "Disk",     OBIT_long, dim, (gpointer)&disk, err);

  dim[0] = strlen(name);
  ObitInfoListPut (info, "FileName", OBIT_string, dim, 
		   (gpointer)name, err);

} /* end ObitFITSAssign */

/**
 * Renames a FITS file
 * \param in   ObitIO on FITS file
 * \param info Associated ObitInfoList
 * \li "newFileName" OBIT_string (?,1,1) New Name of disk file.
 * \li "FileName" OBIT_string (?,1,1)    Old Name of disk file.
 * \li "Disk" OBIT_int (1,1,1)           Disk number
 * \param err  Error stack for any error messages.
 */
void ObitFITSRename(ObitIO *in, ObitInfoList *info, ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOType ft;
  ObitInfoType type;
  olong disk;
  gchar *newName, *oldName, *oldFull, *newFull;
  ObitHistory *inHist=NULL;
  gchar hiCard[73];
  gchar *routine = "ObitFITSRename";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error);
  g_assert (ObitIOIsA(in));

  /* check that FITS type */
  if (!ObitInfoListGet(info, "FileType", &type, dim, &ft, err)) {
    Obit_log_error(err, OBIT_Error, 
		"%s: entry FileType not in InfoList Object %s",	
		   routine, in->name);
    return;
  }
  Obit_return_if_fail ((ft == OBIT_IO_FITS), err,"%s: Object NOT FITS %s",
		       routine, in->name);
  /* Get file names */
  if (!ObitInfoListGet(info, "Disk", &type, dim, &disk, err))
    Obit_traceback_msg (err, routine, in->name);
    
  if(!ObitInfoListInfo(info, "newFileName", &type, dim, err)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry newFileName not in InfoList Object %s",
		   routine, in->name);
    return;
  }
  /* Allocate */
  newName = g_malloc0(dim[0]+1);
  if (!ObitInfoListGet(info, "newFileName", &type, dim, newName, err)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry newFileName not in InfoList Object %s",
		   routine, in->name);
    return;
  }
  /* Final NULL */
  newName[dim[0]] = 0;

  /* Full path */
  newFull = ObitFITSFilename (disk, newName, err);  
     
  /* Does new filename exist? */
  if (ObitFileExist (newFull, err) || err->error) {
    Obit_log_error(err, OBIT_Error, 
		   "%s:File %s already exists",
		   routine, newName);
    g_free(newName);
    return;
  }
  if(!ObitInfoListInfo(info, "FileName", &type, dim, err)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileName not in InfoList Object %s",
		   routine, in->name);
    return;
  }
  /* Allocate */
  oldName = g_malloc0(dim[0]+1);
  if(!ObitInfoListGet(info, "FileName", &type, dim, oldName, err)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileName not in InfoList Object %s",
		   routine, in->name);
    return;
  }
  /* Final NULL */
  oldName[dim[0]] = 0;
     
  /* Full path */
  oldFull = ObitFITSFilename (disk, oldName, err);

  /* Rename */
  ObitFileRename (oldFull, newFull, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Save new filename */
  dim[0] = strlen(newName); dim[1] = 1;
  ObitInfoListAlwaysPut(info, "FileName", OBIT_string, dim, newName);

  /* Add history entry */
  inHist = newObitHistoryValue ("History", info, err);
  ObitHistoryOpen (inHist, OBIT_IO_ReadWrite, err);
  ObitHistoryTimeStamp (inHist, "Obit Rename", err);
  g_snprintf ( hiCard, 72, "Obit / rename from %s",oldName);
  ObitHistoryWriteRec (inHist, -1, hiCard, err);
  g_snprintf ( hiCard, 72, "Obit /  to %s", newName);
  ObitHistoryWriteRec (inHist, -1, hiCard, err);
  ObitHistoryClose (inHist, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  inHist = ObitHistoryUnref(inHist);
  g_free(oldName);
  g_free(newName);
  g_free(newFull);
  g_free(oldFull);

} /* end ObitFITSRename */

/*---------------Private functions---------------------------*/
