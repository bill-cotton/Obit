/* $Id: ObitHistory.c,v 1.23 2007/08/23 14:50:48 bcotton Exp $     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include <sys/types.h>
#include <time.h>
#include "Obit.h"
#include "ObitFileFITS.h"
#include "ObitIOHistory.h"
#include "ObitIOHistoryFITS.h"
#include "ObitIOHistoryAIPS.h"
#include "ObitHistory.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitHistory.c
 * ObitHistory class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitHistory";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitHistoryClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitHistoryClassInfo myClassInfo = {FALSE};

/*----------------------Private functions---------------------------*/

/** Private: Initialize newly instantiated object. */
void  ObitHistoryInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitHistoryClear (gpointer in);

/** Private: Get selection */
static void 
ObitHistoryGetSelect (ObitHistory* in, ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitHistoryClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/

/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitHistory* newObitHistory (gchar* name)
{
  ObitHistory* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitHistoryClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitHistory));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitHistoryInit((gpointer)out);

 return out;
} /* end newObitHistory */

/**
 * Constructor from basic object ObitInfoList
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \param info Parent object list defining the underlying file
 *             e.g. FileType, disk, name for FITS, disk, user, cno for AIPS.
 * \param err  Error stack, returns if not empty.
 * \return the new object.
 */
ObitHistory* 
newObitHistoryValue (gchar* name, ObitInfoList *info, ObitErr *err)
{
  ObitHistory* out = NULL;
  gint32 dim[IM_MAXDIM];
  ObitInfoType type;
  gchar *routine = "newObitHistoryValue";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;

  /* Basic object */
  out = newObitHistory(name);

  /* get file type */
  if(!ObitInfoListGet(info, "FileType", &type, dim, &out->FileType, err))
    Obit_traceback_val (err, routine, name, out);

  /* Copy ObitInfoList to get relevant information */
  out->info = ObitInfoListUnref(out->info); /* out with the old - empty */
  out->info = ObitInfoListCopy(info);       /* in with the new - copy */
  
  return out;
} /* end newObitHistoryValue */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitHistoryGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitHistoryClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitHistoryGetClass */

/**
 * Delete underlying files and the basic object.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 * \return pointer for input object, NULL if deletion successful
 */
ObitHistory* ObitHistoryZap (ObitHistory *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitHistoryZap";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return in;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

 /* Open to fully instantiate */
  if ((in->myStatus!=OBIT_Active) & (in->myStatus!=OBIT_Modified)) {
    retCode = ObitHistoryOpen (in, OBIT_IO_ReadWrite, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, NULL);
  }
  
  /* Delete the file */
  ObitIOZap ((ObitIO*)in->myIO, err);
  if (err->error)
    Obit_traceback_val (err, routine, in->name, in);
  in->myIO = NULL;

  /* Get memory resident bits as well */
  in = ObitHistoryUnref(in);
  
  return in;
} /* end ObitHistoryZap */

/**
 * Make a deep copy of input object.
 * Both objects should be filly defined.
 * \param in  The object to copy, 
 * if underlying structures don't exist, it merely returns
 * without writing the out History.
 * \param out An existing object pointer for output
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitHistory* ObitHistoryCopy (ObitHistory *in, ObitHistory *out, 
			      ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOCode iretCode, oretCode;
  const ObitIOHistoryClassInfo *inClass, *outClass;
  gchar hiCard[73];
  gchar *routine = "ObitHistoryCopy";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* Copy data */
  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitHistoryOpen (in, OBIT_IO_ReadOnly, err);
  /* if it didn't work bail out */
  if ((iretCode!=OBIT_IO_OK) || (err->error)) {
    /* Remove any errors on err */
    ObitErrClearErr(err);
    return out;
  }

  /* test open output */
  oretCode = ObitHistoryOpen (out, OBIT_IO_ReadWrite, err);
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_val (err, routine,in->name, out);
  }

  /* we're in business, copy */
  inClass  = (ObitIOHistoryClassInfo*)in->myIO->ClassInfo;
  outClass = (ObitIOHistoryClassInfo*)out->myIO->ClassInfo;
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    iretCode = inClass->ObitIOHistoryReadRec (in->myIO, -1, hiCard, err);
    if (iretCode!=OBIT_IO_OK) break;
    oretCode = outClass->ObitIOHistoryWriteRec (out->myIO, -1, hiCard, err);
  }
  
  /* check for errors */
  if (((iretCode > OBIT_IO_EOF) && (iretCode!=OBIT_IO_NotFoundErr)) || 
      (oretCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine,in->name, out);
  
  /* close files */
  iretCode = ObitHistoryClose (in, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine,in->name, out);
  
  oretCode = ObitHistoryClose (out, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine,out->name, out);
  
  return out;
} /* end ObitHistoryCopy */

/**
 * Copy HISTORY cards from main file header of in to ObitHistory out,
 * Both objects should be fully defined.
 * \param in  The object to copy
 * \param out An existing object pointer for output 
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitIOCode ObitHistoryCopyHeader (ObitHistory *in, ObitHistory *out, 
				 ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  olong disk, i;
  ObitFileFITS *inFITS;
  gchar hiCardIn[81], hiCardOut[81], FileName[201];
  gint32 dim[IM_MAXDIM];
  ObitInfoType type;
  const ObitIOHistoryClassInfo *outClass;
  gchar *routine = "ObitHistoryCopyHeader";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return OBIT_IO_SpecErr;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* NOP if input isn't FITS */
  if (in->FileType!=OBIT_IO_FITS) return OBIT_IO_OK;

  /* get FITS file info */
  if(!ObitInfoListGet(in->info, "Disk", &type, dim, &disk, err))
    Obit_traceback_val (err, routine, in->name, OBIT_IO_ReadErr);
  
  if (!ObitInfoListGet(in->info, "FileName", &type, dim, FileName, err))
    Obit_traceback_val (err, routine, in->name, OBIT_IO_ReadErr);
  FileName[dim[0]] = 0;


  /* Open input FITS header */
  inFITS = newObitFileFITS(in->name);
  iretCode = ObitFileFITSOpen (inFITS, FileName, disk, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, iretCode);

  /* test open output */
  oretCode = ObitHistoryOpen (out, OBIT_IO_ReadWrite, err);
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_val (err, routine,in->name, oretCode);
  }

  /* we're in business, copy */
  outClass = (ObitIOHistoryClassInfo*)out->myIO->ClassInfo;
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    iretCode = ObitFileFITSReadHistory (inFITS, hiCardIn, err);
    if (iretCode!=OBIT_IO_OK) break;
    for (i=0; i<80; i++) hiCardOut[i] = ' '; hiCardOut[i] = 0;
    strncpy (hiCardOut, &hiCardIn[8], 70); hiCardOut[71] = 0;
    oretCode = outClass->ObitIOHistoryWriteRec (out->myIO, -1, hiCardOut, err);
  }
  
  /* check for errors */
  if (((iretCode > OBIT_IO_EOF) && (iretCode!=OBIT_IO_NotFoundErr)) || 
      (oretCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine,in->name, iretCode);
  
  /* close files */
  iretCode = ObitFileFITSClose (inFITS, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine,in->name, iretCode);
  
  oretCode = ObitHistoryClose (out, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine,out->name, oretCode);
  
  return OBIT_IO_OK;
} /* end ObitHistoryCopyHeader */

/**
 * Copy HISTORY cards from ObitHistory of in to main file header of out.
 * Both objects should be filly defined.
 * \param in  The object to copy
 * \param out Output object for HISTORY header entries.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitIOCode ObitHistoryCopy2Header (ObitHistory *in, ObitHistory *out, 
				   ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  olong disk, morekeys;
  ObitFileFITS *outFITS;
  gchar hiCardIn[81], hiCardOut[81], FileName[201];
  gint32 dim[IM_MAXDIM];
  ObitInfoType type;
  gchar *routine = "ObitHistoryCopy2Header";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return OBIT_IO_SpecErr;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* NOP if output isn't FITS */
  if (out->FileType!=OBIT_IO_FITS) return OBIT_IO_OK;

  /* get FITS file info */
  if(!ObitInfoListGet(out->info, "Disk", &type, dim, &disk, err))
    Obit_traceback_val (err, routine, out->name, OBIT_IO_ReadErr);
  
  if (!ObitInfoListGet(out->info, "FileName", &type, dim, FileName, err))
    Obit_traceback_val (err, routine, out->name, OBIT_IO_ReadErr);
  FileName[dim[0]] = 0;

  /* Open input FITS table */
  iretCode = ObitHistoryOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, iretCode);

  /* test open output */
  outFITS = newObitFileFITS(out->name);
  oretCode = ObitFileFITSOpen (outFITS, FileName, disk, OBIT_IO_ReadWrite, err);
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_val (err, routine,out->name, oretCode);
  }

  /* Expand header to accept new records */
  morekeys = ObitHistoryNumRec(in);
  oretCode = ObitFileFITSAddKeys (outFITS, morekeys, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_val (err, routine,out->name, oretCode);
  }
  
  /* we're in business, copy */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    iretCode = ObitHistoryReadRec (in, -1, hiCardIn, err);
    strncpy (hiCardOut, hiCardIn, 70);
    hiCardOut[70] = 0;  /* to be sure */
    if (iretCode!=OBIT_IO_OK) break;
    oretCode = ObitFileFITSWriteHistory (outFITS, hiCardOut, err);
   }
  
  /* check for errors */
  if (((iretCode > OBIT_IO_EOF) && (iretCode!=OBIT_IO_NotFoundErr)) || 
      (oretCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine,in->name, iretCode);
  
  /* close files */
  iretCode = ObitHistoryClose (in, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine,in->name, iretCode);
  
  oretCode = ObitFileFITSClose (outFITS, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, out->name, oretCode);
  
  return OBIT_IO_OK;
} /* end ObitHistoryCopy2Header */

/**
 * Copy HISTORY cards from main FITS header of in to main file header of out.
 * Both objects should be filly defined.
 * \param in  The object to copy
 * \param out Output object for HISTORY header entries.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitIOCode ObitHistoryHeader2Header (ObitHistory *in, ObitHistory *out, 
				     ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  olong disk;
  ObitFileFITS *outFITS, *inFITS;
  gchar hiCardIn[81], hiCardOut[81], FileName[201];
  gint32 dim[IM_MAXDIM];
  ObitInfoType type;
  gchar *routine = "ObitHistoryHeader2Header";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return OBIT_IO_SpecErr;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* NOP if both aren't FITS */
  if (in->FileType!=OBIT_IO_FITS)  return OBIT_IO_OK;
  if (out->FileType!=OBIT_IO_FITS) return OBIT_IO_OK;

  /* get input FITS file info */
  if(!ObitInfoListGet(in->info, "Disk", &type, dim, &disk, err))
    Obit_traceback_val (err, routine, out->name, OBIT_IO_ReadErr);
  
  if (!ObitInfoListGet(in->info, "FileName", &type, dim, FileName, err))
    Obit_traceback_val (err, routine, out->name, OBIT_IO_ReadErr);
  FileName[dim[0]] = 0;

  /* Open input FITS header */
  inFITS = newObitFileFITS(in->name);
  iretCode = ObitFileFITSOpen (inFITS, FileName, disk, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, iretCode);

  /* get output FITS file info */
  if(!ObitInfoListGet(out->info, "Disk", &type, dim, &disk, err))
    Obit_traceback_val (err, routine, out->name, OBIT_IO_ReadErr);
  
  if (!ObitInfoListGet(out->info, "FileName", &type, dim, FileName, err))
    Obit_traceback_val (err, routine, out->name, OBIT_IO_ReadErr);
  FileName[dim[0]] = 0;

  /* test open output */
  outFITS = newObitFileFITS(out->name);
  oretCode = ObitFileFITSOpen (outFITS, FileName, disk, OBIT_IO_ReadWrite, err);
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_val (err, routine,out->name, oretCode);
  }

  /* we're in business, copy */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    iretCode = ObitFileFITSReadHistory (inFITS, hiCardIn, err);
    strncpy (hiCardOut, hiCardIn, 70);
    if (iretCode!=OBIT_IO_OK) break;
    oretCode = ObitFileFITSWriteHistory (outFITS, hiCardOut, err);
   }
  
  /* check for errors */
  if (((iretCode > OBIT_IO_EOF) && (iretCode!=OBIT_IO_NotFoundErr)) || 
      (oretCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine,in->name, iretCode);
  
  /* close files */
  iretCode = ObitFileFITSClose (inFITS, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine,in->name, iretCode);
  
  oretCode = ObitFileFITSClose (outFITS, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, out->name, oretCode);
  
  return OBIT_IO_OK;
} /* end ObitHistoryHeader2Header */

/**
 * Initialize structures and open file.
 * The image descriptor is read if OBIT_IO_ReadOnly or 
 * OBIT_IO_ReadWrite and written to disk if opened OBIT_IO_WriteOnly.
 * After the file has been opened the member, buffer is initialized
 * for reading/storing the table unless member bufferSize is <0.
 * If the requested version ("Ver" in InfoList) is 0 then the highest
 * numbered table of the same type is opened on Read or Read/Write, 
 * or a new table is created on on Write.
 * The file etc. info should have been stored in the ObitInfoList:
 * \li "FileType" OBIT_long scalar = OBIT_IO_FITS or OBIT_IO_AIPS 
 *               for file type (see class documentation for details).
 * \param in     Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite,
 *               or OBIT_IO_WriteOnly).
 *               If OBIT_IO_WriteOnly any existing data in the output file
 *               will be lost.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitHistoryOpen (ObitHistory *in, ObitIOAccess access, 
			    ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitHistoryOpen";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitHistoryIsA(in));

  /* If the file is already open - close it  first */
  if ((in->myStatus==OBIT_Active) || (in->myStatus==OBIT_Modified)) {
    retCode = ObitIOHistoryClose (in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* get selection parameters */
  ObitHistoryGetSelect (in, err);
  if (err->error) Obit_traceback_val (err, routine,in->name, retCode);

  /* create appropriate ObitIOHistory */
  /* unlink any existing IO structure */
  in->myIO = ObitUnref (in->myIO);
  if (in->FileType==OBIT_IO_FITS) {
    in->myIO = (ObitIOHistory*)newObitIOHistoryFITS(in->name, in->info, err);

  } else if (in->FileType==OBIT_IO_AIPS) {
    in->myIO = (ObitIOHistory*)newObitIOHistoryAIPS(in->name, in->info, err);
  }

  in->myIO->access = access; /* save access type */

 /* most of the instructions for the I/O are in the ObitInfoList */
  retCode = ObitIOHistoryOpen (in->myIO, access, in->info, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, in->name, retCode);

  /* read or write Headers */
  if ((access == OBIT_IO_ReadOnly) || (access == OBIT_IO_ReadWrite)) {
    /* read header info */
    retCode = ObitIOHistoryReadDescriptor(in->myIO, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  } 
  if ((access == OBIT_IO_ReadWrite) || (access == OBIT_IO_WriteOnly)) {
    /* Write header info */
    retCode = ObitIOHistoryWriteDescriptor(in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* init I/O */
  retCode = ObitIOHistorySet (in->myIO, in->info, err);
  if ((retCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);

  /* set Status */
  in->myStatus = OBIT_Active;

  return retCode;
} /* end ObitHistoryOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitHistoryClose (ObitHistory *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitHistoryClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  /* Something going on? */
  if (in->myStatus == OBIT_Inactive) return OBIT_IO_OK;

  /* flush buffer if writing */
  if (((in->myIO->access==OBIT_IO_ReadWrite) || 
       (in->myIO->access==OBIT_IO_WriteOnly))) {

    retCode = ObitIOHistoryFlush (in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);

    /* Update header on disk if writing */
    retCode = ObitIOHistoryWriteDescriptor(in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);    
  }

  /* Close actual file */
  retCode = ObitIOHistoryClose (in->myIO, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* set Status */
  in->myStatus = OBIT_Inactive;

  return OBIT_IO_OK;;
} /* end ObitHistoryClose */

/**
 * Read one record of table data from disk.
 * \param in     Pointer to object to be read.
 * \param recno  Record number to read, -1 = next;
 * \param hiCard Char array to accept line
 * \param err    ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitHistoryReadRec (ObitHistory *in, olong recno, 
				    gchar hiCard[73], ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitHistoryReadRec";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* read record recno */
  retCode = ObitIOHistoryReadRec (in->myIO, recno, hiCard,  err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  return retCode;
} /* end ObitHistoryReadRec */

/**
 * Write a history record
 * \param in    Pointer to object to be written.
 * \param rowno Record number to write, -1 = next;
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitHistoryWriteRec (ObitHistory *in, olong recno, 
				     gchar hiCard[73], ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitHistoryWriteRec";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* write record recno */
  retCode = ObitIOHistoryWriteRec (in->myIO, recno, hiCard, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  return retCode;
} /* end ObitHistoryWriteRec */

/**
 * Add a time stamp and a label to the end of a History
 * \param in    Pointer to object to be written.
 * \param label Label string for record
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitHistoryTimeStamp (ObitHistory *in, 
				 gchar *label, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar line[80];
  struct tm *lp;
  time_t clock;
  olong timea[3], datea[3];
  gchar *routine = "ObitHistoryTimeStamp";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* date and time info */
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

  /* Compose line to write */
  g_snprintf (line,70, "        / %4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2d %s",
	   datea[0],datea[1],datea[2],timea[0], timea[1],timea[2],label);

  /* write row rowno */
  retCode = ObitIOHistoryWriteRec (in->myIO, in->myIO->CurrentNumber+1, line, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  return retCode;
} /* end ObitHistoryTimeStamp */

/**
 * How many history records?
 * \param in    Pointer to open object to be tested
 * \return number of records, <0 => problem
 */
olong ObitHistoryNumRec (ObitHistory *in)
{
  olong out = 0;

  /* error checks */
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  out = ObitIOHistoryNumRec (in->myIO);
  return out;
} /* end ObitHistoryNumRec */

/**
 * Copy values from a list of entries in an ObitInfoList to an open History
 * Only first 64 characters of string values copied
 * \param out  Output object for HISTORY header entries.
 * \param list NULL terminated list of entries in info
 * \param info ObitInfoList with values to copy
 * \param err  Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitIOCode 
ObitHistoryCopyInfoList (ObitHistory *out, gchar *pgmName, gchar *list[], 
			 ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gpointer     xdata;
  gboolean     found, *bdata;
  olong        i, j, more, indx, ltemp, lstr, *ldata, size;
  olong         *idata;
  oint         *odata;
  ofloat       *fdata;
  odouble      *ddata;
  gchar        hicard[81], bchar, *cdata, cstring[65];
  const ObitIOHistoryClassInfo *outClass;
  gchar *routine = "ObitHistoryyCopyInfoList";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(out, &myClassInfo));

  outClass = (ObitIOHistoryClassInfo*)out->myIO->ClassInfo;
  /* loop through list copying elements */
  i = 0;
  while (list[i]) {
    found = ObitInfoListGetP(info, list[i], &type, dim, &xdata);
    
    if (found) { /* copy by type to history */

      /* element count */
      size = dim[0];
      for (j=1; j<MAXINFOELEMDIM; j++) size *= MAX (1, dim[j]);

      switch (type) { 
      case OBIT_int:
	idata = (olong*)xdata;
	g_snprintf (hicard, 80, "%s %s = ", pgmName, list[i]);
	more = size;
	indx = strlen (hicard);
	while (more>0) {
	  for (j=0; j<8; j++) {
	    g_snprintf (&hicard[indx], 80-indx, "%d ", *idata++);
	    indx = strlen (hicard);
	    more--;                    /* finished? */
	    if (more<=0) break;
	  }
	  outClass->ObitIOHistoryWriteRec (out->myIO, -1, hicard, err);
	  if (err->error) Obit_traceback_val (err, routine, out->name, retCode);
	  g_snprintf (hicard, 80, "%s   ", pgmName);
	  indx = strlen (hicard);
	}

	break;
      case OBIT_oint:
	odata = (oint*)xdata;
	g_snprintf (hicard, 80, "%s %s = ",  pgmName, list[i]);
	more = size;
	indx = strlen (hicard);
	while (more>0) {
	  for (j=0; j<20; j++) {
	    ltemp = (olong)(*odata++);
	    g_snprintf (&hicard[indx], 80-indx, " %d ", ltemp);
	    indx = strlen (hicard);
	    more--;                    /* finished? */
	    if (more<=0) break;
	    if (indx>60) break;        /* Line full? */
	  }
	  outClass->ObitIOHistoryWriteRec (out->myIO, -1, hicard, err);
	  if (err->error) Obit_traceback_val (err, routine, out->name, retCode);
	  g_snprintf (hicard, 80, "%s   ", pgmName);
	  indx = strlen (hicard);
	}

	break;
      case OBIT_long:
	ldata = (olong*)xdata;
	g_snprintf (hicard, 80, "%s %s = ", pgmName, list[i]);
	indx = strlen (hicard);
	more = size;
	while (more>0) {
	  for (j=0; j<20; j++) {
	    ltemp = (olong)(*ldata++);
	    g_snprintf (&hicard[indx], 80-indx, " %d ", ltemp);
	    indx = strlen (hicard);
	    more--;                    /* finished? */
	    if (more<=0) break;
	    if (indx>60) break;        /* Line full? */
	  }
	  outClass->ObitIOHistoryWriteRec (out->myIO, -1, hicard, err);
	  if (err->error) Obit_traceback_val (err, routine, out->name, retCode);
	  g_snprintf (hicard, 80, "%s   ", pgmName);
	  indx = strlen (hicard);
	}

	break;
      case OBIT_float:
	fdata = (ofloat*)xdata;
	g_snprintf (hicard, 80, "%s %s = %15.5g ",
		    pgmName, list[i], *fdata++);
	more = size - 1;
	/* Can add second value in this line */
	if (more>=1) {
	  indx = strlen (hicard);
	  g_snprintf (&hicard[indx], 80-indx, "%15.5g ", *fdata++);
	  more--;                    /* finished? */
	}
	outClass->ObitIOHistoryWriteRec (out->myIO, -1, hicard, err);
	if (err->error) Obit_traceback_val (err, routine, out->name, retCode);
	g_snprintf (hicard, 80, "%s ", pgmName);
	indx = strlen (hicard);
	while (more>0) {
	  for (j=0; j<4; j++) {
	    g_snprintf (&hicard[indx], 80-indx, "%15.5g ", *fdata++);
	    indx = strlen (hicard);
	    more--;                    /* finished? */
	    if (more<=0) break;
	    if (indx>55) break;   /* Line full? */
	  }
	  outClass->ObitIOHistoryWriteRec (out->myIO, -1, hicard, err);
	  if (err->error) Obit_traceback_val (err, routine, out->name, retCode);
	  g_snprintf (hicard, 80, "%s   ", pgmName);
	  indx = strlen (hicard);
	}


	break;
      case OBIT_double:
	ddata = (odouble*)xdata;
	g_snprintf (hicard, 80, "%s %s = %25.12lg ",
		    pgmName, list[i], *ddata++);
	outClass->ObitIOHistoryWriteRec (out->myIO, -1, hicard, err);
	if (err->error) Obit_traceback_val (err, routine, out->name, retCode);
	more = size - 1;
	g_snprintf (hicard, 80, "%s ", pgmName);
	indx = strlen (hicard);
	while (more>0) {
	  for (j=0; j<2; j++) {
	    g_snprintf (&hicard[indx], 80-indx, "%25.12lg ", *ddata++);
	    indx = strlen (hicard);
	    more--;                    /* finished? */
	    if (more<=0) break;
	    if (indx>45) break;   /* Line full? */
	  }
	  outClass->ObitIOHistoryWriteRec (out->myIO, -1, hicard, err);
	  if (err->error) Obit_traceback_val (err, routine, out->name, retCode);
	  g_snprintf (hicard, 80, "%s   ", pgmName);
	  indx = strlen (hicard);
	}

	break;
      case OBIT_string:   /* only 64 char of string */
	cdata = (gchar*)xdata;
	lstr = dim[0];  /* length of string */
	strncpy (cstring, cdata, MIN (lstr, 64));
	cstring[MIN (lstr, 64)] = 0;  /* null terminate */
	cdata += lstr;         /* move down string array */
	g_snprintf (hicard, 80, "%s %s = '%s' ",
		    pgmName, list[i], cstring);
	outClass->ObitIOHistoryWriteRec (out->myIO, -1, hicard, err);
	if (err->error) Obit_traceback_val (err, routine, out->name, retCode);
	more = (size / lstr) - 1;
	g_snprintf (hicard, 80, "%s ", pgmName);
	indx = strlen (hicard);
	while (more>0) {
	  for (j=0; j<2; j++) {
	    strncpy (cstring, cdata, MIN (lstr, 64));
	    cstring[MIN (lstr, 64)] = 0;  /* null terminate */
	    cdata += lstr;         /* move down string array */
	    g_snprintf (&hicard[indx], 80-indx, "'%s' ", cstring);
	    indx = strlen (hicard);
	    more--;                    /* finished? */
	    if (more<=0) break;
	    if (indx>40) break;   /* Line full? */
	  }
	  outClass->ObitIOHistoryWriteRec (out->myIO, -1, hicard, err);
	  if (err->error) Obit_traceback_val (err, routine, out->name, retCode);
	  g_snprintf (hicard, 80, "%s   ", pgmName);
	  indx = strlen (hicard);
	}

	break;
      case OBIT_bool:
	bdata = (gboolean*)xdata;
	g_snprintf (hicard, 80, "%s %s = ", pgmName, list[i]);
	indx = strlen (hicard);
	more = size;
	while (more>0) {
	  for (j=0; j<30; j++) {
	    if (*bdata++) bchar = 'T';
	    else bchar = 'F';
	    g_snprintf (&hicard[indx], 80-indx, "%c ", bchar);
	    indx = strlen (hicard);
	    more--;                    /* finished? */
	    if (more<=0) break;
	    if (indx>60) break;   /* Line full? */
	  }
	  outClass->ObitIOHistoryWriteRec (out->myIO, -1, hicard, err);
	  if (err->error) Obit_traceback_val (err, routine, out->name, retCode);
	  g_snprintf (hicard, 80, "%s   ", pgmName);
	  indx = strlen (hicard);
	}

	break;
      default:
	g_assert_not_reached(); /* unknown, barf */
      }; /* end switch to copy by type */
    }
    /* loop */
    i++;  /* next item */
  } /* end loop over list */
  return retCode;
} /* end ObitHistoryCopyInfoList */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitHistoryClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();
  myClassInfo.hasScratch  = FALSE; /* Scratch files not allowed */

  /* Set function pointers */
  ObitHistoryClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitHistoryClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitHistoryClassInfoDefFn (gpointer inClass)
{
  ObitHistoryClassInfo *theClass = (ObitHistoryClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitHistoryClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitHistoryClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitHistoryGetClass;
  theClass->newObit       = (newObitFP)newObitHistory;
  theClass->ObitCopy      = (ObitCopyFP)ObitHistoryCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitHistoryClear;
  theClass->ObitInit      = (ObitInitFP)ObitHistoryInit;
  /* New this class */
  theClass->ObitHistoryOpen  = (ObitHistoryOpenFP)ObitHistoryOpen;
  theClass->ObitHistoryClose = (ObitHistoryCloseFP)ObitHistoryClose;
  theClass->ObitHistoryZap   = (ObitHistoryZapFP)ObitHistoryZap;

  /* *************** CHANGE HERE *********************************  */

} /* end ObitHistoryClassDefFn */

/*---------------Private functions--------------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitHistoryInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitHistory *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread    = newObitThread();
  in->info      = newObitInfoList();
  in->myIO      = NULL;

} /* end ObitHistoryInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitHistory* cast to an Obit*.
 */
void ObitHistoryClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitHistory *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  in->myIO      = ObitUnref(in->myIO);
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitHistoryClear */


/**
 * Get requested information from the ObitInfoList (FileType)
 * \param in   Pointer to Object
 * \param err  ObitErr for reporting errors.
 */
static void 
ObitHistoryGetSelect (ObitHistory* in, ObitErr *err)
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitHistoryIsA(in));

  /* what type of underlying file? */
  if (!ObitInfoListGet(in->info, "FileType", &type, dim, 
		       &in->FileType, err)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "ObitHistoryGetSelect: entry FileType not in InfoList Object %s",
		   in->name);
  }

} /* end ObitHistoryGetSelect */


