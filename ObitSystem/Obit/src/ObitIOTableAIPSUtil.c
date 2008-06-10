/* $Id: ObitIOTableAIPSUtil.c,v 1.1 2006/06/26 16:46:43 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
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

#include "ObitIOTableAIPSUtil.h"
#include "ObitAIPS.h"
#include "ObitAIPSCat.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOTableAIPSUtil.c
 * ObitIOTableAIPSUtil -  AIPS table utilities
 */

/*----------------------Public functions---------------------------*/

/**
 * Read Table list from an AIPS catalog entry
 * \param in  Pointer to object to be read.
 * \param err ObitErr for reporting errors.
 */
void ObitIOTableAIPSUtilReadTableList(ObitData *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *HeaderFile=NULL;
  olong disk, CNO, UserId;
  AIPSint buffer[260];
  ObitFilePos wantPos;
  ObitFile *myFile=NULL;
  gsize size;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar *routine = "ObitIOTableAIPSUtilReadTableList";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(in));

  /* File info */
  if(!ObitInfoListGet(in->info, "Disk", &type, dim, &disk, err))
    Obit_traceback_msg (err, routine, in->name);
  if(!ObitInfoListGet(in->info, "User", &type, dim, &UserId, err))
    Obit_traceback_msg (err, routine, in->name);
  if(!ObitInfoListGet(in->info, "CNO", &type, dim, &CNO, err))
    Obit_traceback_msg (err, routine, in->name);
  
  /* Set file name */
  HeaderFile = ObitAIPSFilename (OBIT_AIPS_Header, disk, CNO, 
				 UserId, NULL, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* swallow file - open */
  myFile = newObitFile(in->name);
  size = 260 * sizeof(AIPSint);
  ObitFileOpen (myFile, HeaderFile, OBIT_IO_ReadOnly, 
		OBIT_IO_Binary, size, err);
  g_free(HeaderFile); HeaderFile = NULL;  /* cleanup */
  if  (err->error)  Obit_traceback_msg (err, routine, in->name);
  
  /* read */
  wantPos = 0;
  retCode = ObitFileRead (myFile, wantPos, size, (gchar*)buffer, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, in->name);

  /* Get table Info to Table List */
  ObitAIPSCatGetTable ((ObitTableList*)in->tableList, (gchar*)buffer,
		       UserId, disk, CNO, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

 /* close file */
  retCode = ObitFileClose (myFile, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, in->name);

  /* delete */
  myFile = ObitFileUnref(myFile);
}  /* End of ObitIOTableAIPSUtilReadTableList */

/**
 * Write Table list to an AIPS catalog entry
 * \param in  Pointer to object to be update
 * \param err ObitErr for reporting errors.
 */
void ObitIOTableAIPSUtilWriteTableList(ObitData *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *HeaderFile=NULL;
  olong disk, CNO, UserId;
  AIPSint buffer[260];
  ObitFilePos wantPos;
  ObitFile *myFile=NULL;
  gsize size;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar *routine = "ObitIOTableAIPSUtilReadTableList";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(in));

  /* File info */
  if(!ObitInfoListGet(in->info, "Disk", &type, dim, &disk, err))
    Obit_traceback_msg (err, routine, in->name);
  if(!ObitInfoListGet(in->info, "User", &type, dim, &UserId, err))
    Obit_traceback_msg (err, routine, in->name);
  if(!ObitInfoListGet(in->info, "CNO", &type, dim, &CNO, err))
    Obit_traceback_msg (err, routine, in->name);
  
  /* Set file name */
  HeaderFile = ObitAIPSFilename (OBIT_AIPS_Header, disk, CNO, 
				 UserId, NULL, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* swallow file - open */
  myFile = newObitFile(in->name);
  size = 260 * sizeof(AIPSint);
  ObitFileOpen (myFile, HeaderFile, OBIT_IO_ReadWrite, 
		OBIT_IO_Binary, size, err);
  g_free(HeaderFile); HeaderFile = NULL;  /* cleanup */
  if  (err->error)  Obit_traceback_msg (err, routine, in->name);
  
  /* read */
  wantPos = 0;
  retCode = ObitFileRead (myFile, wantPos, size, (gchar*)buffer, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, in->name);

  /* Update table Info from Table List */
  ObitAIPSCatSetTable (in->tableList, (gchar*)buffer, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

   /* Now (re)write it */
  wantPos = 0; /* File location */
  /* write it  */
  retCode = ObitFileWrite (myFile, wantPos, size, (gchar*)buffer, err);
  if ((retCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, in->name         );

  /* close file */
  retCode = ObitFileClose (myFile, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, in->name);

  /* delete */
  myFile = ObitFileUnref(myFile);
} /* End of ObitIOTableAIPSUtilWriteTableList */

/**
 * Get an associated table
 * \param in  Pointer to object with tables
 * \param err ObitErr for reporting errors.
 * \return Table or NULL on error
 */
ObitTable* ObitIOTableAIPSUtilGetTable(ObitData *in,ObitIOAccess access, 
				       gchar *tabType, olong *tabVer, 
				       ObitErr *err)
{
  ObitTable *out = NULL;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  olong version;
  olong disk, CNO, UserId;
  gboolean gotIt;
  gchar ttype[3], *outName, tabName[51];
  gchar *routine = "ObitIOTableFITSUtilGetTable";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitDataIsA(in));
  g_assert(tabType!=NULL);
  g_assert(tabVer!=NULL);

  /* the Tablelist object must be present */
  if (in->tableList==NULL) {
      Obit_log_error(err, OBIT_Error, 
		     "my tableList member is NULL, open %s first", 
		     in->name);
      return NULL;
  }

  /* Do we already have this one? */
  version = *tabVer;
  gotIt = ObitTableListGet ((ObitTableList*)in->tableList, tabType, &version, 
			    &out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);

  /* Check if we're forcing a new table */
  if ((access==OBIT_IO_WriteOnly) && (*tabVer <= 0)) {
    version++;
    out = ObitTableUnref(out);
  }

  /* Set output table version */
  *tabVer = version;
  
  if (gotIt && (out!=NULL)) return out; /* that was easy */
  
  /* If it doesn't exist and request is read only - return NULL */
  if ((!gotIt) && (access==OBIT_IO_ReadOnly)) return NULL;

   /* Create one - make descriptive name */
  g_snprintf (tabName, 50, "%s table %d for ",tabType, *tabVer);
  outName =  g_strconcat (tabName, in->name, NULL);
  out = newObitTable (outName);
  g_free(outName);

  /* File info */
  if(!ObitInfoListGet(in->info, "Disk", &type, dim, &disk, err))
    Obit_traceback_val (err, routine, in->name, NULL);
  if(!ObitInfoListGet(in->info, "User", &type, dim, &UserId, err))
    Obit_traceback_val (err, routine, in->name, NULL);
  if(!ObitInfoListGet(in->info, "CNO", &type, dim, &CNO, err))
    Obit_traceback_val (err, routine, in->name, NULL);

  /* Setup info needed for access */
  ttype[0] = tabType[5]; ttype[1] = tabType[6]; ttype[2] = 0;
  ObitTableSetAIPS(out, disk, CNO, ttype, version, UserId, 1, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
  
  /* register it in the TableList */
  ObitTableListPut (in->tableList, tabType, &version, out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
  
  /* Force Write to disk 
  ObitIOUVAIPSWriteDescriptor (in, err);*/
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
  
  return out;
} /* end ObitIOTableFITSUtilGetTable */

