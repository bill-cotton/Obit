/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2009                                          */
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
#include "ObitAIPS.h"
#include "ObitAIPSCat.h"
#include "ObitHistory.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitAIPS.c
 * ObitAIPS class function definitions.
 */

/*-----------------File Globals ---------------------------*/
/** number of bytes per "sector" in ancient aipsish */
static olong AIPS_NBPS = 256*sizeof(AIPSint); 

/** Current AIPS data revision code */
gchar AIPS_Version = 'D';

/** Class information structure */
static ObitAIPS ObitAIPSInfo = {"ObitAIPS", FALSE};

/** Pointer to Class information structure */
static ObitAIPS *myAIPSInfo = &ObitAIPSInfo;

/* AIPS magic blanking value */
static union FBLANKequiv {
  gchar string[4];
  ofloat fblank;
} FBLANK;

/*------------------ Structures -----------------------------*/

/*---------------Private function prototypes----------------*/
/*---------------Public functions---------------------------*/
/**
 * Save names of the AIPS directories.
 * \param number of disks defined [0,MAXAIPSDISK]], >=0 -> none
 * \param dir the names of the directories
 *        If NULL, check for environment variables $DA01...
 * \param F_TRUE   Value of Fortran TRUE (used in Fortran interface)
 * \param F_FALSE  Value of Fortran FALSE
 */
void ObitAIPSClassInit (olong number, gchar* dir[], oint F_TRUE, oint F_FALSE)
{
  olong i;
  gchar daxx[5], *ev;
  gchar *da[]={"01","02","03","04","05","06","07","08","09","0A",
	       "0B","0C","0D","0E","0F","0G","0H","0I","0J","0K",
	       "0L","0M","0N","0O","0P","0Q","0R","0S","0T","0U",
	       "0V","0W","0X","0Y","0Z"};

  myAIPSInfo->NumberDisks = 0;
  for (i=0; i<MAXAIPSDISK; i++) {
    myAIPSInfo->AIPSdir[i] = NULL;
    myAIPSInfo->noScrat[i] = FALSE;
  }

  /* now initialized */
  myAIPSInfo->initialized = TRUE;

  /* Are directories given or should I look for them? */
  if (dir==NULL) { /* Look in environment */
    for (i=0; i<MAXAIPSDISK; i++) {
      sprintf (daxx,"DA%s",da[i]);
      ev = getenv (daxx);
      if (ev) {
	myAIPSInfo->AIPSdir[myAIPSInfo->NumberDisks] =
	  /* strip DA?? = and add / */
	  g_strconcat(ev, "/", NULL);
	myAIPSInfo->NumberDisks++;
      } else {
	break;
      }
    }
  } else { /* use what are given */
    
    /* error checks */
    g_assert (number<=MAXAIPSDISK);
    
    if (number<=0) {
      myAIPSInfo->NumberDisks = 0;
      myAIPSInfo->AIPSdir[0] = NULL;
      return;
    }
    
    /* save directory names */
    myAIPSInfo->NumberDisks = number;
    for (i=0; i<number; i++) 
      myAIPSInfo->AIPSdir[i] =  g_strconcat(dir[i], "/", NULL);
    
  } /* end of initialize data directories */

  /* initialize catalog header structure info in ObitAIPSCat class */
  ObitAIPSCatInitDHDR();

  /* Save true and false */
  myAIPSInfo->F_TRUE = F_TRUE;
  myAIPSInfo->F_FALSE = F_FALSE;

  /* initialise AIPS magic blanking value float equiv of 'INDE' */
  FBLANK.fblank = ObitMagicF();
  
} /* end ObitAIPSClassInit */

/**
 * Frees directory strings
 */
void ObitAIPSShutdown (void)
{
  olong i;

  myAIPSInfo->NumberDisks = 0;
  for (i=0; i<MAXAIPSDISK; i++) {
    if (myAIPSInfo->AIPSdir[i]) g_free(myAIPSInfo->AIPSdir[i]);
    myAIPSInfo->AIPSdir[i] = NULL;
  }

  /* now uninitialized */
  myAIPSInfo->initialized = FALSE;

} /*  end ObitAIPSShutdown */

/**
 * Forms file name from the various parts.
 * #ObitAIPSClassInit must have been used to initialize.
 * \param type File type code
 * \param disk AIPS "disk" number. 1-rel
 * \param cno AIPS catalog slot number.
 * \param userid user number.
 * \param tabType two character code for table type, 
 *                NULL if not needed.
 * \param tabVer table version number.
 * \param err    Error stack for any error messages.
 * \return file name string, should be g_freeed when done.
 */
gchar* 
ObitAIPSFilename (ObitAIPSFileType type, olong disk, olong cno, 
		  olong userid, gchar *tabType, olong tabVer, ObitErr *err)
{
  gchar *out;
  gchar idEhex[8], cnoEhex[8], verEhex[8], file[20];
  gchar *types[] = {"CA","CB","MA","UV","SC","UK","HI","PL","SL"};
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  if (!myAIPSInfo->initialized) /* AIPS directories un initialized */
    Obit_log_error(err, OBIT_Error, 
		   "AIPS directories uninitialized");
  if ((disk<=0) || (disk>myAIPSInfo->NumberDisks)) /* Disk number out of range */
    Obit_log_error(err, OBIT_Error, 
		   "AIPS disk number %d out of range [%d,%d]", 
		   disk, 1, myAIPSInfo->NumberDisks);
  if (err->error) return NULL;

  /* put together basic file name */
  /* File type */
  file[0]=types[type][0]; file[1]=types[type][1];
  
  /* date revision code */
  file[2] = AIPS_Version;
  
  /* catalog slot */
  /* Convert to EHex */
  ObitAIPSEHex(cno, (gchar*)cnoEhex);
  file[3] = cnoEhex[0]; file[4] = cnoEhex[1]; file[5] = cnoEhex[2];
  
  /* version by type */
  if (type==OBIT_AIPS_Catalog) {
    verEhex[0] = '0'; verEhex[1] = '0';verEhex[2] = '0';
  } else if (type==OBIT_AIPS_Table) {
    ObitAIPSEHex(tabVer, (gchar*)verEhex);

    /* set table type */
    file[0] = tabType[0]; file[1] = tabType[1];
  } else { /* Everything else is one */
    verEhex[0] = '0'; verEhex[1] = '0';verEhex[2] = '1';
  }
  file[6] = verEhex[0]; file[7] = verEhex[1]; file[8] = verEhex[2];

  
  /* User id */
  /* Convert to EHex */
  ObitAIPSEHex(userid, (gchar*)idEhex);
  file[9] = '.';        file[10] = idEhex[0];
  file[11] = idEhex[1]; file[12] = idEhex[2];
  
  /* Top it all off with a VAXism */
  file[13] = ';';         file[14] = 0;
  
  /* put it all together */
  out = g_strconcat (ObitAIPSDirname(disk, err), file, NULL);
  
  return out;
} /* end ObitAIPSFilename */

/**
 * Change directory name
 * Note, the maximum number of AIPS disks is MAXAIPSDISK  new disks up to this 
 * max can be assigned.
 * \param disk  disk (1-rel) to be changed, must have been assigned at startup
 *              if <=0 then add new disk, up to MAXAIPSDISK
 * \param dir   new name of the directory
 * \return disk number assigned
 */
olong ObitAIPSSetDirname (olong disk, gchar* dir, ObitErr *err)
{
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return disk;
  if (!myAIPSInfo->initialized) /* AIPS directories un initialized */
    Obit_log_error(err, OBIT_Error, "AIPS directories uninitialized");
  if (disk<=0) disk = myAIPSInfo->NumberDisks+1;  /* Add directory? */
  if ((disk<1) || (disk>MAXAIPSDISK)) /* Disk number out of range */
    Obit_log_error(err, OBIT_Error, 
		   "AIPS disk number %d out of range [%d,%d]", 
		   disk, 1, myAIPSInfo->NumberDisks);
  if (err->error) return  disk;

  /* replace directory name */
  if (myAIPSInfo->AIPSdir[disk-1]) g_free(myAIPSInfo->AIPSdir[disk-1]);
  myAIPSInfo->AIPSdir[disk-1] =  g_strdup(dir);
  myAIPSInfo->noScrat[disk-1] = FALSE;
  myAIPSInfo->NumberDisks = MAX (myAIPSInfo->NumberDisks, disk);

  return disk;
} /* end ObitAIPSSetDirname */

/**
 * Lookup and possibly add a directory.
 * Searches the current list of directories for directory name dir,
 * if it exists, its disk number is returned.  Otherwise the directory
 * is added and its number is returned.
 * Note, the maximum number of AIPS disks is MAXAIPSDISK  new disks up to this 
 * max can be assigned.
 * \param disk  disk (1-rel) to be changed, must have been assigned at startup
 *              if <=0 then add new disk, up to MAXAIPSDISK
 * \param dir   name of the directory to be located or added
 * \return disk number assigned
 */
olong ObitAIPSFindDirname (gchar* dir, ObitErr *err)
{
  olong disk=0;

  /* error checks */
  if (err->error) return disk;

  if (!myAIPSInfo->initialized) { /* AIPS directories un initialized */
    Obit_log_error(err, OBIT_Error, "AIPS directories uninitialized");
    return disk;
  }

  /* See if it is currently defined */
  for (disk=1; disk<=myAIPSInfo->NumberDisks; disk++) {
    if (!strcmp(dir, myAIPSInfo->AIPSdir[disk-1])) return disk;
  }

    /* Nope - add */
    disk = myAIPSInfo->NumberDisks+1;  /* Add directory */
    if ((disk<1) || (disk>MAXAIPSDISK)) /* Disk number out of range? */
      Obit_log_error(err, OBIT_Error, 
		     "AIPS disk number %d out of range [%d,%d]", 
		     disk, 1, myAIPSInfo->NumberDisks);
    if (err->error) return  disk;

  /* Add directory name */
  if (myAIPSInfo->AIPSdir[disk-1]) g_free(myAIPSInfo->AIPSdir[disk-1]);
  myAIPSInfo->AIPSdir[disk-1] =  g_strdup(dir);
  myAIPSInfo->noScrat[disk-1] = FALSE;
  myAIPSInfo->NumberDisks = MAX (myAIPSInfo->NumberDisks, disk);

  return disk;
} /* end ObitAIPSFindDirname */

/**
 * Returns pointer to directory string by AIPS disk.
 * \param disk AIPS disk number.
 * \param err  Error stack for any error messages.
 * \return directory name string, this is a pointer into a global 
 *         class structure and should not be g_freeed.
 */
gchar* ObitAIPSDirname (olong disk, ObitErr *err)
{
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  if (!myAIPSInfo->initialized) /* AIPS directories uninitialized */
    Obit_log_error(err, OBIT_Error, 
		   "AIPS directories uninitialized");
  if ((disk<=0) || (disk>myAIPSInfo->NumberDisks)) /* Disk number out of range */
    Obit_log_error(err, OBIT_Error, 
		   "AIPS disk number %d out of range [%d,%d]", 
		   disk, 1, myAIPSInfo->NumberDisks);
  if (myAIPSInfo->AIPSdir[disk-1]==NULL) /* Directory not defined */
    Obit_log_error(err, OBIT_Error, 
		   "AIPS directory %d not defined", disk);
  if (err->error) return NULL;

  return myAIPSInfo->AIPSdir[disk-1];
} /* ObitAIPSDirname  */

/**
 * Returns number of defined AIPS disks
 * \param err  Error stack for any error messages.
 * \return number of disks
 */
olong ObitAIPSGetNumDisk (ObitErr *err)
{
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return 0;
  if (!myAIPSInfo->initialized) /* AIPS directories un initialized */
    Obit_log_error(err, OBIT_Error, 
		   "AIPS directories uninitialized");
  if (err->error) return 0;

  g_assert(myAIPSInfo->initialized);

  return myAIPSInfo->NumberDisks;
} /* ObitAIPSGetNumDisk  */

/**
 * Calculate first byte offset in an AIPS image file of a given pixel.
 * AIPS images have planes starting on an even #AIPS_NBPS byte boundry.
 * If the row size is less than AIPS_NBPS bytes, multiple rows are
 * stored in a sector but are not allowed to cross a sector boundry;
 * the end of a sector is padded and the next row starts at the beginning
 * of a new sector.
 * Rows longer than a sector have the last sector padded and the next 
 * row starts at the beginning of a new sector.
 * This was necessary back at the dawn of time when MODCOMP computers 
 * roamed the earth and could only start a transfer at the beginning 
 * of a disk sector.
 * This routine patterned after the aipsish $APLSUB/COMOFF.FOR.
 * (The original is incomprehensible).
 * \param naxis number of axes in image
 * \param naxes number of pixels on each axis.
 * \param pos   1-rel pixel position
 * \return byte offset from beginning of image file 
 */
ObitFilePos ObitAIPSImageFileOffset (olong naxis, olong *naxes, olong *pos)
{
  ObitFilePos nspp, nrps, nspr, np, ns, filePos=0;

  /* error checks */
  g_assert((naxis>0) && (naxis<=IM_MAXDIM));
  g_assert(naxes!=NULL);
  g_assert(pos!=NULL);

  /* how many sectors per plane? */
  nrps = AIPS_NBPS / (naxes[0] * sizeof(float)); /* # rows per sector */
  if (nrps > 0) { /* multiple rows per sector */
    nspp = (naxes[1] - 1) / nrps + 1;
  } else { /* multiple sectors per row */
    nspp = (naxes[0] * sizeof(float) - 1) / AIPS_NBPS + 1;
    nspp *= naxes[1];
  }

  /* how many planes down? */
  np = 0;
   /*if (naxis>2) np += (pos[2]-1);          third dimension */
  np += (pos[2]-1);  /* Always pad to end of plane */
  if (naxis>3) np += (pos[3]-1)*naxes[2]; /* 4th dim*/
  if (naxis>4) np += (pos[4]-1)*(naxes[3]*naxes[2]); /* 5th */
  if (naxis>5) np += (pos[5]-1)*(naxes[4]*naxes[3]*naxes[2]); /* 6th */
  if (naxis>6) np += (pos[6]-1)*(naxes[5]*naxes[4]*naxes[3]*naxes[2]); /* 7th */

  /* byte offset to the beginning of the plane */
  filePos = np * nspp * AIPS_NBPS;

  /* correction for row */
  if (nrps > 0) { /* multiple rows per sector */
    /* integral number of rows per sector */
    /* number of whole sectors since beginning of plane: */
    ns = (pos[1]-1)/nrps;  
    filePos += ns * AIPS_NBPS;
    /* fractional sector: */
    filePos += (pos[1]-ns*nrps-1) * (naxes[0] * sizeof(float));
  } else { /* multiple sectors per row */
    /* each row takes an integral number of sectors */
    /* # sectors per row */
    nspr = 1 + (((naxes[0] * sizeof(float))-1)/ AIPS_NBPS); 
    filePos += nspr * AIPS_NBPS * (pos[1]-1);
  }

  /* correction for column */
  filePos += (ObitFilePos)(pos[0]-1) * sizeof(float);

  return filePos;
}  /* end ObitAIPSImageFileOffset */


/**
 * Calculate first byte offset in an AIPS table file of the beginning of
 * a row.  Tables rows are written either an integral number of sectors
 * (#AIPS_NBPS bytes) per row, or an integral number of rows per sector.
 * The last "sector" is padded after the last datum.
 * \param start Byte offset of the beginning of the row data.
 * \param lrow  The length of a row in bytes.
 * \param row   1-rel row number desired.
 * \return byte offset from beginning of table file 
 */
ObitFilePos ObitAIPSTableFileOffset (ObitFilePos start, olong lrow, olong row)
{
  ObitFilePos nrps, nspr, ns, nr, filePos=0;
  olong rowoff;

  /* error checks */
  g_assert(start>0);
  g_assert(lrow>0);
  g_assert(row>0);

  /* how many sectors per row? */
  rowoff = row - 1;
  nrps = AIPS_NBPS / lrow; /* # rows per sector */
  if (nrps > 0) { /* multiple rows per sector */
    /* How many whole sectors */
    ns = rowoff / nrps;
    /* How many rows in this sector? */
    nr = rowoff % nrps;
  } else { /* multiple sectors per row */
    /* Number of sectors per row */
    nspr = 1 + (lrow-1) / AIPS_NBPS;
    /* How many whole sectors */
    ns = rowoff * nspr;
    /* How many rows in this sector? */
    nr = 0;
  }

  /* put the pieces together */
  filePos = start + ns * AIPS_NBPS + nr * lrow;

  return filePos;
}  /* end ObitAIPSTableFileOffset */

/**
 * Calculate where the end of the current "ModComp" AIPSish sector ends.
 * Tables rows are written either an integral number of sectors
 * (#AIPS_NBPS bytes) per row, or an integral number of rows per sector.
 * The last "sector" is padded after the last datum.
 * This routine calculates the location of the end of this sector.
 * \param start Byte offset of the beginning of the row data.
 * \param lrow  The length of a row in bytes.
 * \param nrow   Number of rows in the table.
 * \return byte offset from beginning of table file 
 */
ObitFilePos ObitAIPSTableEOF (ObitFilePos start, olong lrow, olong nrow)
{
  ObitFilePos nrps, nspr, ns, filePos=0;

  /* error checks */
  g_assert(start>0);
  g_assert(lrow>0);

  /* how many sectors per row? */
  nrps = AIPS_NBPS / lrow; /* # rows per sector */
  if (nrps > 0) { /* multiple rows per sector */
    /* How many whole sectors */
    ns = (olong)(((ofloat)nrow / (ofloat)nrps) + 0.9999);
  } else { /* multiple sectors per row */
    /* Number of sectors per row */
    nspr = 1 + (lrow-1) / AIPS_NBPS;
    /* How many whole sectors */
    ns = nrow * nspr;
  }

  /* put the pieces together */
  filePos = start + ns * AIPS_NBPS;

  return filePos;
}  /* end ObitAIPSTableEOF */

/**
 * Determine target padding size for wonky AIPS uvdata file
 * size rules.  Must be an integral number of AIPS blocks
 * \param curPos  The byte offset of the current end of the file.
 * \return target file position index in bytes to be padded to.
 */
ObitFilePos ObitAIPSUVWonkyPad (ObitFilePos curPos)
{
  ObitFilePos ns, out = curPos;

  /* If it's already an integral multiple of AIPS_NBPS, nothing more needed */
  if ((curPos %  AIPS_NBPS) == 0) return out;

  /* Need to pad the current AIPS block */
  /* number of whole blocks */
  ns = curPos / AIPS_NBPS;

  out = (ns + 1) * AIPS_NBPS;
  return out;
} /* end ObitAIPSUVWonkyPad */

/**
 * Converts an integer into an extended (base 36) Hex string
 * Only works up to 3 EHex digits.
 * \param in integer to be converted
 * \param out preexisting string into which to write value
 *            must be at least 4 characters.
 */
void ObitAIPSEHex (olong in, gchar *out)
{
  gchar *hexc = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  olong i, j, iv, jv, work = in;

  /* error tests */
  if (in<0) {out[0]='B'; out[1]='A';out[2]='D';out[3]=0;return;}
  g_assert(out!=NULL);

  /* init */
  out[0] = '0'; out[1] = '0'; out[2] = '0'; out[4] = 0; 

  /* loop doing converting */
  for (i=0; i<3; i++) {
    j = 2 - i;
    iv = work / 36;
    jv = work - iv*36;
    out[j] = hexc[jv];
    work = iv;
  }

} /* end ObitAIPSEHex */

/**
 * Assigns scratch file naming information and writes to the info.
 * Deliberately does not follow AIPS conventions.
 * name  = "OBIT SCRATCH"
 * class = "pgmName+pgmNumber"
 * seq   = scrNo
 * Makes name "pgmName+pgmNumber+'Scr'+scrNo"
 * \param pgmName    Program name
 * \param pgmNumber  Program incarnation number
 * \param type       AIPS file type ("MA", "UV").
 * \param user       AIPS user number
 * \param disk       AIPS disk number.
 * \param scrNo      Which scratch file number
 * \param info       ObitInfoList to write to
 * \param err        Error stack for any error messages.
 */
void ObitAIPSAssign(gchar *pgmName, olong pgmNumber, gchar *type,
		    olong user, olong disk, olong scrNo, ObitInfoList *info, 
		    ObitErr *err)
{
  gchar name[13]="OBIT SCRATCH", class[7];
  olong seq, cno;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean exist;
  ObitIOType ft;

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(info));

  /* form AIPSish class, seq */
  g_snprintf (class, 7, "%s%d", pgmName, pgmNumber);
  /* Make sure "POPS" number gets into class */
  if (strlen(pgmName)>=6) g_snprintf (&class[5], 2, "%d", pgmNumber);
  seq = scrNo;

  /* get assignment */
  cno = ObitAIPSDirAlloc(disk, user, name, class, type, seq, 
			  &exist, err);
  if (cno<0)
    Obit_log_error(err, OBIT_Error, 
		   "ERROR assigning AIPS scratch file CNO");

  /* write assignment to info */
  ft = OBIT_IO_AIPS;
  ObitInfoListPut (info, "FileType", OBIT_long, dim, (gpointer)&ft,   err);
  ObitInfoListPut (info, "Disk",     OBIT_long, dim, (gpointer)&disk, err);
  ObitInfoListPut (info, "CNO",      OBIT_long, dim, (gpointer)&cno,  err);
  ObitInfoListPut (info, "User",     OBIT_long, dim, (gpointer)&user, err);

} /* end ObitAIPSAssign */

/**
 * Renames a AIPS file
 * \param in   ObitIO on AIPS file
 * \param info Associated ObitInfoList
 * \li "Disk" OBIT_long (1,1,1)           Disk number
 * \li "CNO" OBIT_long (1,1,1)            Catalog slot number
 * \li "newName" OBIT_string (12,1,1) New AIPS Name 
 *      absent or Blank = don't change
 * \li "newClass" OBIT_string (6,1,1) New AIPS Class
 *      absent or Blank = don't changeO
 * \li "newSeq" OBIT_long (1,1,1) New AIPS Sequence
 *      0 => unique value
 * \param err  Error stack for any error messages.
 */
void ObitAIPSRename(ObitIO *in, ObitInfoList *info, ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOType ft;
  ObitInfoType type;
  olong disk, user, CNO, newSeq=-1, oldSeq, test;
  gchar newName[13], newClass[7], oldName[13], oldClass[7];
  ObitAIPSDirCatEntry* entry=NULL;
  ObitHistory *inHist=NULL;
  ObitErr *xerr=NULL;
  gchar hiCard[73];
  gchar *routine = "ObitAIPSRename";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error);
  g_assert (ObitIOIsA(in));

  /* check that AIPS type */
  if (!ObitInfoListGet(info, "FileType", &type, dim, &ft, err)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",	
		   routine, in->name);
    return;
  }
  Obit_return_if_fail ((ft == OBIT_IO_AIPS), err,"%s: Object NOT AIPS %s",
		       routine, in->name);
  /* Get naming info */
  if (!ObitInfoListGet(info, "Disk", &type, dim, &disk, err))
    Obit_traceback_msg (err, routine, in->name);
    
  if (!ObitInfoListGet(info, "User", &type, dim, &user, err))
    Obit_traceback_msg (err, routine, in->name);
    
  if (!ObitInfoListGet(info, "CNO", &type, dim, &CNO, err))
    Obit_traceback_msg (err, routine, in->name);
    
  if (!ObitInfoListGet(info, "newName", &type, dim, newName, err))
    Obit_traceback_msg (err, routine, in->name);
  newName[MIN(12,dim[0])] = 0;
    
  if (!ObitInfoListGet(info, "newClass", &type, dim, newClass, err))
    Obit_traceback_msg (err, routine, in->name);
  newClass[MIN(6,dim[0])] = 0;

  if (!ObitInfoListGet(info, "newSeq", &type, dim, &newSeq, err))
    Obit_traceback_msg (err, routine, in->name);
    
  /* Get current info */
  entry = ObitAIPSDirGetEntry (disk, user, CNO, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Get old values */
  oldSeq = entry->seq;
  g_memmove(oldName,  entry->name, 12); oldName[12] = 0;
  g_memmove(oldClass, entry->class, 6); oldClass[6] = 0;

  /* Default name to old */
  if (!strncmp (newName,  "    ", 4)) 
    {g_memmove(newName,  entry->name, 12); newName[12] = 0;}
  if (!strncmp (newClass, "    ", 4)) 
    {g_memmove(newClass, entry->class, 6); newClass[6] = 0;}

  /* Default sequence number */
  if (newSeq<=0) {
    newSeq = ObitAIPSDirHiSeq(disk, user, newName, newClass, 
			      entry->type, FALSE, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } else { /* make sure it doesn't already exist */
    /* See if already entry with name, and class,
       Note: unlike POPS, this will allow the same seq. no. 
       as exisits on another disk */
    /* Any failure assumed to be caused by nonexistance */
    xerr = newObitErr();
    test = ObitAIPSDirFindCNO(disk, user, newName, newClass, 
			      entry->type, newSeq, xerr);
    if (xerr->error) test=-1;
    xerr = ObitErrUnref(xerr);
    /* test */
    Obit_return_if_fail ((test==-1), err,
			 "%s: Entry %s.%s seq %d disk %d already exists",
			 routine, newName, newClass, newSeq, disk);
  }

  /* OK - change catalog entry */
  ObitAIPSDirRename(disk, user, CNO, newName, newClass, newSeq, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Update catalog entry */
  ObitAIPSCatRename(disk, user, CNO, newName, newClass, newSeq, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Add history entry */
  inHist = newObitHistoryValue ("History", info, err);
  ObitHistoryOpen (inHist, OBIT_IO_ReadWrite, err);
  ObitHistoryTimeStamp (inHist, "Obit Rename", err);
  g_snprintf ( hiCard, 72, "Obit / rename from %s.%s.%d",
	       oldName, oldClass, oldSeq);
  ObitHistoryWriteRec (inHist, -1, hiCard, err);
  g_snprintf ( hiCard, 72, "Obit /  to %s.%s.%d",
	       newName, newClass, newSeq);
  ObitHistoryWriteRec (inHist, -1, hiCard, err);
  ObitHistoryClose (inHist, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  inHist = ObitHistoryUnref(inHist);
    
} /* end ObitAIPSRename */

/**
 * Convert Fortran Logical to gboolean
 * \param logical Fortran LOGICAL value
 * \return gboolean equivalent (TRUE or FALSE)
 */
gboolean ObitAIPSBooleanF2C (oint logical)
{
  gboolean out;
  out = logical==myAIPSInfo->F_TRUE;
  return out;
} /* end ObitAIPSBooleanF2C */

/**
 * Convert gboolean to Fortran Logical
 * \param bool boolean (TRUE or FALSE)
 * \return Fortran LOGICAL value equivalent
 */
 oint ObitAIPSBooleanC2F (gboolean bool)
{
  oint out;

  if (bool) out = myAIPSInfo->F_TRUE;
  else      out = myAIPSInfo->F_FALSE;
  return out;
} /* end ObitAIPSBooleanC2F */

/**
 * Convert Mark AIPS data directory as "noScrat"
 * i.e. don't add scratch files on this data area
 * \param disk       AIPS disk number.
 * \param noScrat    If True don't put scratch files on disk
 *                   else allow.
 * \param err  Error stack for any error messages.
 */
void ObitAIPSnoScrat(olong disk, gboolean noScrat, ObitErr *err)
{
  gchar *routine = "ObitAIPSnoScrat";
  Obit_return_if_fail (((disk>0) && (disk<=myAIPSInfo->NumberDisks)), err,
			"%s: disk %d not in ranrfs [1,%d]",
		       routine, disk, myAIPSInfo->NumberDisks);
  myAIPSInfo->noScrat[disk-1] = noScrat;
} /* end ObitAIPSnoScrat */

/**
 * Tell is scratch files allowed on AIPS directory
 * \param disk       AIPS disk number.
 * \return True if scratch files disallowed or invalid disk
 */
gboolean ObitAIPSisNoScrat(olong disk)
{
  if ((disk<=0) || (disk>myAIPSInfo->NumberDisks)) return TRUE;
  return myAIPSInfo->noScrat[disk-1];
} /* end ObitAIPSisNoScrat */

/**
 * Check for "noScrat" entry in info and set all positive integer 
 * values (OBIT_int,OBIT_long or OBIT_oint) as noScratch disks.
 * \param info   List to check
 * \param err  Error stack for any error messages.
 */
void ObitAIPSSetnoScrat(ObitInfoList *info, ObitErr *err)
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, *iarr=NULL;
  gchar *routine = "ObitAIPSSetnoScrat";

  if (err->error) return;

  ObitInfoListGetP(info, "noScrat",  &type, dim, (gpointer)&iarr);
  if ((iarr==NULL) || (dim[0]<=0) || ((type!=OBIT_int) && (type!=OBIT_long) && 
				      (type!=OBIT_oint))) return;

  for (i=0; i<dim[0]; i++) {
    if ((iarr[i]>0) && (iarr[i]<=myAIPSInfo->NumberDisks))
      ObitAIPSnoScrat (iarr[i], TRUE, err);
    if (err->error) Obit_traceback_msg (err, routine, " ");
  }
  
} /* end ObitAIPSSetnoScrat */

/*---------------Private functions---------------------------*/

