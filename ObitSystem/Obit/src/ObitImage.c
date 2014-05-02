/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2014                                          */
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

#include "ObitImageDesc.h"
#include "ObitImageSel.h"
#include "ObitImage.h"
#include "ObitIOImageFITS.h"
#include "ObitIOImageAIPS.h"
#include "ObitFITS.h"
#include "ObitAIPSDir.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitHistory.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImage.c
 * ObitImage class function definitions.
 * This class is derived from the ObitData base class.
 *
 * This class contains an astronomical image and allows access.
 * An ObitImage is the front end to a persistent disk resident structure.
 *
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitImage";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitDataGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitImageClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitImageClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitImageInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitImageClear (gpointer in);

/** Private: Read selection parameters from ObitInfoList. */
static void ObitImageGetSelect (ObitInfoList *info, ObitImageDesc* desc, 
				ObitImageSel *sel, ObitErr *err);

/** Private: Determine overall plane number. */
static olong PlaneNumber (olong plane[5], olong naxis, olong *inaxes);

/** Private: Assign myIO object */
static void ObitImageSetupIO (ObitImage *in, ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitImageClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitImage* newObitImage (gchar* name)
{
  ObitImage* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageClassInit();

  /* allocate/init structure */
  out = ObitMemAlloc0Name(sizeof(ObitImage), "ObitImage");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitImageInit((gpointer)out);

 return out;
} /* end newObitImage */

/**
 * Create an Image object with selection parameters set from an InfoList
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param inList InfoList to extract object information from
 * Following InfoList entries for AIPS files ("xxx" = prefix):
 * \li xxxName  OBIT_string  AIPS file name
 * \li xxxClass OBIT_string  AIPS file class
 * \li xxxDisk  OBIT_oint    AIPS file disk number
 * \li xxxSeq   OBIT_oint    AIPS file Sequence number
 * \li AIPSUser OBIT_oint    AIPS User number
 * \li xxxCNO   OBIT_oint    AIPS Catalog slot number
 * \li xxxDir   OBIT_string  Directory name for xxxDisk
 *
 * Following entries for FITS files ("xxx" = prefix):
 * \li xxxFile  OBIT_string  FITS file name
 * \li xxxDisk  OBIT_oint    FITS file disk number
 * \li xxxDir   OBIT_string  Directory name for xxxDisk
 *
 * For all File types:
 * \li xxxFileType OBIT_string "UV" = UV data, "MA"=>image, "Table"=Table, 
 *                "OTF"=OTF, etc
 * \li xxxDataType OBIT_string "AIPS", "FITS"
 *                 Defaults to value of "DataType"
 *    
 * For xxxDataType = "MA"
 * \li xxxBLC   OBIT_oint[7] (Images only) 1-rel bottom-left corner pixel
 * \li xxxTRC   OBIT_oint[7] (Images Only) 1-rel top-right corner pixel
 * \li xxxExist OBIT_bool (1,1,1) If True, file expected to exist (def TRUE)
 * \param err     ObitErr for reporting errors.
 * \return new data object with selection parameters set
 */
ObitImage* ObitImageFromFileInfo (gchar *prefix, ObitInfoList *inList, 
				  ObitErr *err)
{
  ObitImage    *out = NULL;
  ObitInfoType type;
  olong        Aseq, AIPSuser, disk, cno, i;
  gchar        *strTemp, inFile[129], stemp[256];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gboolean     exist;
  gchar        *keyword=NULL, *DataTypeKey = "DataType", *DataType=NULL;
  gchar        *routine = "ObitImageFromFileInfo";

  if (err->error) return out;  /* Previous error? */

  /* Create output */
  out = newObitImage (prefix);

  /* Is it expected to exist? */
  if (prefix) keyword = g_strconcat (prefix, "Exist", NULL);
  else        keyword = g_strdup("nVisPIO");
  exist = TRUE;
  ObitInfoListGetTest(inList, keyword, &type, dim, &exist);
  g_free(keyword);

  /* BLC */
  if (prefix) keyword = g_strconcat (prefix, "BLC", NULL);
  else        keyword = g_strdup("BLC");
  ObitInfoListGetTest(inList, keyword, &type, dim, blc); /* BLC */
  g_free(keyword);
  
  /* BLC */
  if (prefix) keyword = g_strconcat (prefix, "TRC", NULL);
  else        keyword = g_strdup("TRC");
  ObitInfoListGetTest(inList, keyword, &type, dim, trc); /* TRC */
  g_free(keyword);

  /* File type - could be either AIPS or FITS */
  if (prefix) keyword =  g_strconcat (prefix, DataTypeKey, NULL);
  else keyword =  g_strconcat (DataTypeKey, NULL);
  if (!ObitInfoListGetP (inList, keyword, &type, dim, (gpointer)&DataType)) {
    /* Try "DataType" */
    if (!ObitInfoListGetP(inList, "DataType", &type, dim, (gpointer)&DataType)) {
      /* couldn't find it - add message to err and return */
      Obit_log_error(err, OBIT_Error, 
		     "%s: entry %s not in InfoList", routine, keyword);
      g_free(keyword);
      return out;
    }
  }
  g_free(keyword);

  if (!strncmp (DataType, "AIPS", 4)) { /* AIPS */
    /* AIPS disk */
    if (prefix) keyword = g_strconcat (prefix, "Disk", NULL);
    else        keyword = g_strdup("Disk");
    ObitInfoListGet(inList, keyword, &type, dim, &disk, err);
    g_free(keyword);

    /* If prefixDir given, lookup disk number */
    if (prefix) keyword = g_strconcat (prefix, "Dir", NULL);
    else        keyword = g_strdup("Dir");
    if (ObitInfoListGetP(inList, keyword, &type, dim, (gpointer)&strTemp)) {
      /* make sure NULL terminated */
      strncpy (stemp, strTemp, MIN(255,dim[0])); stemp[MIN(255,dim[0])] = 0;
      disk = ObitAIPSFindDirname (stemp, err);
      if (err->error) Obit_traceback_val (err, routine, routine, out);
    }
    g_free(keyword);

    /* AIPS name */
    if (prefix) keyword = g_strconcat (prefix, "Name", NULL);
    else        keyword = g_strdup("Name");
    if (ObitInfoListGetP(inList, keyword, &type, dim, (gpointer)&strTemp)) {
      strncpy (Aname, strTemp, 13);
    } else { /* Didn't find */
      strncpy (Aname, "No Name ", 13);
    } 
    Aname[12] = 0;
    g_free(keyword);

    /* AIPS class */
    if (prefix) keyword = g_strconcat (prefix, "Class", NULL);
    else        keyword = g_strdup("Class");
    if  (ObitInfoListGetP(inList, keyword, &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    Aclass[6] = 0;
    g_free(keyword);

    /* input AIPS sequence */
    if (prefix) keyword = g_strconcat (prefix, "Seq", NULL);
    else        keyword = g_strdup("Seq");
    ObitInfoListGet(inList, keyword, &type, dim, &Aseq, err);
    g_free(keyword);

    /* if ASeq==0 want highest existing sequence */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
      if (err->error) Obit_traceback_val (err, routine, "inList", out);
      /* Save on inList*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(inList, "inSeq", OBIT_oint, dim, &Aseq);
    } 

    /* AIPS User no. */
    ObitInfoListGet(inList, "AIPSuser", &type, dim, &AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "inList", out);    

    /* Find/assign catalog number */
    if (exist) 
      cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    else { /* Create */
      cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
      Obit_log_error(err, OBIT_InfoErr, "Making AIPS image %s %s %d on disk %d cno %d",
		     Aname, Aclass, Aseq, disk, cno);
  }
    if (err->error) Obit_traceback_val (err, routine, "inList", out);
    
    /* define object */
    ObitImageSetAIPS (out, OBIT_IO_byPlane, disk, cno, AIPSuser, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "inList", out);
    
  } else if (!strncmp (DataType, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    if (prefix) keyword = g_strconcat (prefix, "File", NULL);
    else        keyword = g_strdup("File");
    if (ObitInfoListGetP(inList, keyword, &type, dim, (gpointer)&strTemp)) {
      strncpy (inFile, strTemp, 128);
    } else { 
      strncpy (inFile, "No_Filename_Given", 128);
    }
    g_free(keyword);
    
    /* input FITS disk */
    if (prefix) keyword = g_strconcat (prefix, "Disk", NULL);
    else        keyword = g_strdup("Disk");
    ObitInfoListGet(inList, keyword, &type, dim, &disk, err);
    g_free(keyword);

    /* If prefixDir given, lookup disk number */
    if (prefix) keyword = g_strconcat (prefix, "Dir", NULL);
    else        keyword = g_strdup("Dir");
    if (ObitInfoListGetP(inList, keyword, &type, dim, (gpointer)&strTemp)) {
      /* make sure NULL terminated */
      strncpy (stemp, strTemp, MIN(255,dim[0])); stemp[MIN(255,dim[0])] = 0;
      disk = ObitFITSFindDir (stemp, err);
      if (err->error) Obit_traceback_val (err, routine, routine, out);
    }
    g_free(keyword);

    /* define object */
    ObitImageSetFITS (out, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "inList", out);
    
    /* Tell about it if didn't exist */
    if (!exist)
      Obit_log_error(err, OBIT_InfoErr, "Making FITS image %s on disk %d",
		     inFile, disk);

  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   routine, DataType);
    return out;
  }

  /* Ensure out fully instantiated and OK if it exists */
  if (exist) ObitImageFullInstantiate (out, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "inList", out);

  /* Set defaults BLC, TRC - use size on myIO as blc, trc incorporated into myDesc */
  if (exist) {
    for (i=0; i<IM_MAXDIM; i++) {
      if (blc[i]<=0) blc[i] = 1;
      blc[i] = MAX (1,  blc[i]);
      if (trc[i]<=0) trc[i] = ((ObitImageDesc*)out->myIO->myDesc)->inaxes[i];
      trc[i] = MIN (trc[i], ((ObitImageDesc*)out->myIO->myDesc)->inaxes[i]);
    }
  }

  /* Save blc, trc */
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (out->info, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (out->info, "TRC", OBIT_long, dim, trc);

  /* Copy any InfoList Parameters */
  if (prefix) keyword = g_strconcat (prefix, "Info", NULL);
  else        keyword = g_strdup("Info");
  ObitInfoListCopyWithPrefix (inList, out->info, keyword, TRUE);
  
  return out;
} /* end ObitImageFromFileInfo */

/**
 * Create a scratch file suitable for accepting the data to be read from in.
 * A scratch Image is more or less the same as a normal Image except that it is
 * automatically deleted on the final unreference.
 * The output will have the underlying files of the same type as in already 
 * allocated.
 * \param in  The object to copy, info may have
 * \li ScrSize OBIT_int (?,1,1) Dimension of the desired scratch Image
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitImage* newObitImageScratch (ObitImage *in, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitImage *out=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, size[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *outName;
  gchar *routine = "newObitImageScratch";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Ensure in fully instantiated -assume OK if myIO exists */
  if (!in->myIO) ObitImageFullInstantiate (in, TRUE, err);
  if (err->error)Obit_traceback_val (err, routine, in->name, out);

  /* Create - derive object name */
  outName = g_strconcat ("Scratch Copy: ",in->name,NULL);
  out = newObitImage(outName);
  g_free(outName);

  /* Mark as scratch */
  out->isScratch = TRUE;

   /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* Copy descriptor */
  out->myDesc = (gpointer)ObitImageDescCopy(in->myDesc, out->myDesc, err);

  /* Check if different size needed */
  if (ObitInfoListGetTest(in->info, "ScrSize", &type, dim, size)) {
    for (i=0; i<MIN (dim[0], IM_MAXDIM); i++) 
      if (size[i]>0) out->myDesc->inaxes[i] = size[i];
  }
 
  /* Force to float pixels */
  out->myDesc->bitpix=-32;

  /* Allocate underlying file */
  ObitSystemGetScratch (in->mySel->FileType, "MA", out->info, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  
  /* Register in the scratch file list */
  ObitSystemAddScratch ((Obit*)out, err);
  
  /* same size IO as input */
  dim[0] = 1;
  ObitInfoListPut (out->info, "IOBy", OBIT_long, dim, &in->myDesc->IOsize, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /* Fully instantiate output */
  ObitImageFullInstantiate (out, FALSE, err);
  if (err->error)Obit_traceback_val (err, routine, out->name, out);
 
  return out;
} /* end newObitImageScratch */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitImageGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObiImagetGetClass */

/**
 * Test if two ObitImages have the same underlying structures.
 * This test is done using values entered into the #ObitInfoList
 * in case the object has not yet been opened.
 * \param in1 First object to compare
 * \param in2 Second object to compare
 * \param err ObitErr for reporting errors.
 * \return TRUE if to objects have the same underlying structures
 * else FALSE
 */
gboolean ObitImageSame (ObitImage *in1, ObitImage *in2, ObitErr *err )
{
  /* Call ObitData function */
  return ObitDataSame ((ObitData*)in1, (ObitData*)in2, err);
} /* end ObitImageSame */

/**
 * Delete underlying files and the basic object.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 * \return pointer for input object, NULL if deletion successful
 */
ObitImage* ObitImageZap (ObitImage *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitImageZap";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return in;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* Any actual I/O? If not just delete object */
  if (in->mySel->FileType==OBIT_IO_MEM) return ObitImageUnref(in);

  /* Close if still active */
  if ((in->myStatus == OBIT_Active) || (in->myStatus == OBIT_Modified)){
   retCode = ObitIOClose(in->myIO, err);
   if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
     Obit_traceback_val (err, routine, in->name, in);    
  }

  /* Ensure in fully instantiated */
  ObitErrLog(err); /* Show any pending messages as they may get lost */
  ObitImageFullInstantiate (in, TRUE, err);
  /* If this fails, clear errors and assume it doesn't exist */
  if (err->error) { 
    ObitErrClearErr(err); 
    return ObitImageUnref(in); 
  }

  /* Delete Image and all tables  */
  ObitIOZap (in->myIO, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, in);

  /* If it's scratch remove from list */
  if (in->isScratch) ObitSystemFreeScratch ((Obit*)in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, in);

  /* Delete object */
  in->isScratch = 0; /* Already deleted underlying structures */
  while (in) in = ObitImageUnref(in);
  
  return in;
} /* end ObitImageZap */

/**
 * Rename underlying files
 * New name information depends on the underlying file type and is
 * given on the info member.
 * For FITS files:
 * \li "newFileName" OBIT_string (?,1,1) New Name of disk file.
 *
 * For AIPS:
 * \li "newName" OBIT_string (12,1,1) New AIPS Name 
 *      Blank = don't change
 * \li "newClass" OBIT_string (6,1,1) New AIPS Class
 *     Blank = don't changeO
 * \li "newSeq" OBIT_int (1,1,1) New AIPS Sequence
 *      0 => unique value
 * \param in Pointer to object to be renamed.
 * \param err ObitErr for reporting errors.
 */
void ObitImageRename (ObitImage *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitImageRename";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* Any actual I/O? */
  if (in->mySel->FileType==OBIT_IO_MEM) return;

  /* Close if still active */
  if ((in->myStatus == OBIT_Active) || (in->myStatus == OBIT_Modified)){
   retCode = ObitIOClose(in->myIO, err);
   if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
     Obit_traceback_msg (err, routine, in->name);    
  }

  /* Ensure in fully instantiated */
  ObitErrLog(err); /* Show any pending messages as they may get lost */
  ObitImageFullInstantiate (in, TRUE, err);
  /* If this fails, clear errors and assume it doesn't exist */
  if (err->error) { 
    ObitErrClearErr(err); 
    return; 
  }

  /* Rename Image */
  ObitIORename (in->myIO, in->info, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  return;
} /* end ObitImageRename */

/**
 * Make a deep copy of input object.
 * Copies are made of complex members including disk files; these 
 * will be copied applying whatever selection is associated with the input.
 * Objects should be closed on input and will be closed on output.
 * In order for the disk file structures to be copied, the output file
 * must be sufficiently defined that it can be written.
 * The copy will be attempted but no errors will be logged until
 * both input and output have been successfully opened.
 * ObitInfoList and ObitThread members are only copied if the output object
 * didn't previously exist.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitImage* ObitImageCopy (ObitImage *in, ObitImage *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOCode iretCode, oretCode;
  gboolean oldExist;
  ObitHistory *inHist=NULL, *outHist=NULL;
  gchar *outName=NULL, *today=NULL;
  gchar *routine = "ObitImageCopy";

  /* error checks */
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Any actual I/O? Can't handle Memory only */
  if (in->mySel->FileType==OBIT_IO_MEM) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Image %s in memory only and I do not know how to copy it", 
		   routine, in->name);
    return NULL;
  }

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitImage(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes other additions only if out newly created */
  if (!oldExist) {
    /* copy */
    out->myDesc = ObitImageDescCopy(in->myDesc, out->myDesc, err);
    /* Don't copy selector */
    if (out->mySel) out->mySel = ObitUnref (out->mySel);
    out->mySel = newObitImageSel (out->name);
    /* Don't copy info */
    /*out->info = ObitInfoListUnref(out->info); */
    /*out->info = ObitInfoListRef(in->info); */
    /* Output will initially have no associated tables */
    out->tableList = ObitTableListUnref(out->tableList);
    out->tableList = newObitTableList(out->name);
    /* don't copy ObitThread  */
}

  /* If the output object was created this call it cannot be fully
     defined so we're done */
  if (!oldExist) return out;

  /* if input has file designated, copy data */
  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitImageOpen (in, OBIT_IO_ReadOnly, err);
  /* if it didn't work bail out */
  if ((iretCode!=OBIT_IO_OK) || (err->error)) {
    return out;
  }

  /* copy Descriptor - this time with full information */
  out->myDesc = ObitImageDescCopy(in->myDesc, out->myDesc, err);
  /* Float it */
  out->myDesc->bitpix = -32;

  /* Creation date today */
  today = ObitToday();
  strncpy (out->myDesc->date, today, IMLEN_VALUE);
  if (today) g_free(today);
 
  /* use same data buffer on input and output 
     so don't assign buffer for output */
  out->extBuffer = TRUE;

  /* test open output */
  ObitErrLog(err); /* Show any pending messages as they may get lost */
  oretCode = ObitImageOpen (out, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitImageOpen (out, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset external buffer */
    out->extBuffer = FALSE;
    return out;
  }

  /* Copy any history  unless Scratch */
  if (!in->isScratch && !out->isScratch) {
    inHist  = newObitDataHistory((ObitData*)in, OBIT_IO_ReadOnly, err);
    outHist = newObitDataHistory((ObitData*)out, OBIT_IO_WriteOnly, err);
    outHist = ObitHistoryCopy (inHist, outHist, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
    inHist  = ObitHistoryUnref(inHist);
    outHist = ObitHistoryUnref(outHist);
  }

  /* make sure the access sizes are the same */
  if (in->myDesc->IOsize != out->myDesc->IOsize) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Access sizes of two images differ", routine);
    Obit_log_error(err, OBIT_Error, 
		   "Objects %s %s", in->name, out->name);
    /* unset external buffer */
    out->extBuffer = FALSE;
    return out;
  }
  
  /* we're in business, copy */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    iretCode = ObitImageRead (in, in->image->array, err);
    if (iretCode!=OBIT_IO_OK) break;
    oretCode = ObitImageWrite (out, in->image->array, err);
  }
  
  /* unset external buffer */
  out->extBuffer = FALSE;
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, out);
  
  /* close files to be sure */
  iretCode = ObitImageClose (in, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, out);
  
  /* close files to be sure */
  oretCode = ObitImageClose (out, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, out->name, out);
  
  return out;
} /* end ObitImageCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an Image similar to the input one.
 * Output version set to floating pixels
 * \param in  The object to copy, info may have
 * \param out An existing object pointer for output, info may have
 * \li Size OBIT_int (?,1,1) Dimension of the desired  Image
 * \param err Error stack, returns if not empty.
 */
void ObitImageClone  (ObitImage *in, ObitImage *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOCode iretCode, oretCode;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, size[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ObitHistory *inHist=NULL, *outHist=NULL;
  gchar *today=NULL;
  gchar *exclude[]={"AIPS CC", "AIPS HI", "AIPS PL", "AIPS SL", NULL}; 
  gchar *routine = "ObitImageClone";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* Any actual I/O? Can't handle Memory only */
  if (in->mySel->FileType==OBIT_IO_MEM) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Image %s in memory only and I do not know how to clone it", 
		   routine, in->name);
    return;
  }

  /* Ensure in fully instantiated and OK  */
  ObitImageFullInstantiate (in, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes other additions */
  /* Don't copy selector */
  if (out->mySel) out->mySel = ObitUnref (out->mySel);
  out->mySel = newObitImageSel (out->name);
  /* Output will initially have no associated tables */
  out->tableList = ObitTableListUnref(out->tableList);
  out->tableList = newObitTableList(out->name);
  /* don't copy ObitThread  */

  /* Open to fully instantiate input and see if it's OK */
  iretCode = ObitImageOpen (in, OBIT_IO_ReadOnly, err);
  if ((iretCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, in->name);

  /* copy Descriptor */
  out->myDesc = ObitImageDescCopy(in->myDesc, out->myDesc, err);

  /* Check if different size needed */
  if (ObitInfoListGetTest(out->info, "Size", &type, dim, size)) {
    for (i=0; i<MIN (dim[0], IM_MAXDIM); i++) 
      if (size[i]>0) out->myDesc->inaxes[i] = size[i];
  }
 
  /* Creation date today */
  today = ObitToday();
  strncpy (out->myDesc->date, today, IMLEN_VALUE);
  if (today) g_free(today);
 
  /* Force to float pixels */
  out->myDesc->bitpix=-32;

  /* Open output */
  oretCode = ObitImageOpen (out, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitImageOpen (out, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset external buffer */
    out->extBuffer = FALSE;
    Obit_traceback_msg (err, routine, out->name);
  }

  /* Copy any history unless Scratch  */
  if (!in->isScratch && !out->isScratch) {
    inHist  = newObitDataHistory((ObitData*)in, OBIT_IO_ReadOnly, err);
    outHist = newObitDataHistory((ObitData*)out, OBIT_IO_WriteOnly, err);
    outHist = ObitHistoryCopy (inHist, outHist, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    inHist  = ObitHistoryUnref(inHist);
    outHist = ObitHistoryUnref(outHist);
  }

 /* Copy tables  */
  iretCode = ObitImageCopyTables (in, out, exclude, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Close files */
  ObitImageClose (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  ObitImageClose (out, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
} /* end ObitImageClone */

/**
 * Define image which has cells on the same grid as in2 but covering the 
 * region of in1.  
 * Only sets descriptor and data buffer members on out.
 * Suitable for memory only image, out->mySel->FileType set to OBIT_IO_MEM.
 * Returns without error if in1 cannot be projected onto in2.
 * Routine translated from the AIPSish 4MASS/SUB/FLATEN.FOR/FLTMSC
 * \param in1  Image whose region is to be covered, must have actual descriptor
 * \param in2  Image whose geometry is to be copied , must have actual descriptor
 * \param out  An existing object pointer for output.  
 *             If the image (FArray) member exists, it is resized if needed.
 *             Descriptor updated to reflect any resizing
 * \param err Error stack, returns if not empty.
 */
void ObitImageClone2  (ObitImage *in1, ObitImage *in2, ObitImage *out, 
		       ObitErr *err)
{
  olong   i, icx, icy;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong naxis[IM_MAXDIM], onaxis[IM_MAXDIM];
  ofloat      fact, pix[IM_MAXDIM], pix2[IM_MAXDIM], pix3[IM_MAXDIM];
  ofloat      crpix[IM_MAXDIM], crpix2[IM_MAXDIM];
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gchar *today=NULL;
  gchar *routine = "ObitImageClone2";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));
  /* Check Equinox */
  Obit_return_if_fail(((in1->myDesc->equinox==1950.0) || (in1->myDesc->equinox==2000.0)), err,
		      "%s: Bad equinox %f", routine, in1->myDesc->equinox);
  /* Check Equinox */
  Obit_return_if_fail(((in2->myDesc->equinox==1950.0) || (in2->myDesc->equinox==2000.0)), err,
		      "%s: Bad equinox %f", routine, in2->myDesc->equinox);

  /* Set output to full image, plane at a time */
  dim[0] = IM_MAXDIM;
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  ObitInfoListPut (out->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (out->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (out->info, "IOBy", OBIT_long, dim, &IOsize, err);

  /* Reference pixels */
  for (i=0; i<IM_MAXDIM; i++) crpix[i]  = 1.0;
  for (i=0; i<IM_MAXDIM; i++) crpix2[i] = 1.0;
  crpix[0]  = in1->myDesc->crpix[0]; crpix[1]  = in1->myDesc->crpix[1];
  crpix2[0] = in2->myDesc->crpix[0]; crpix2[1] = in2->myDesc->crpix[1];

  /* Creation date today */
  today = ObitToday();
  strncpy (out->myDesc->date, today, IMLEN_VALUE);
  if (today) g_free(today);
 
  /* Input image size */
  for (i=0; i<IM_MAXDIM; i++) naxis[i]  = in1->myDesc->inaxes[i];

  /* Center of subimage */
  /*pix[0] = crpix[0];
    pix[1] = crpix[1];*/
  pix[0] = in1->myDesc->inaxes[0]/2;
  pix[1] = in1->myDesc->inaxes[1]/2;
  pix[2] = 1;
  pix[3] = 1;
  pix[4] = 1;

  /* Find corresponding pixel in in2 */
  ObitImageDescCvtPixel (in1->myDesc, in2->myDesc, pix, pix2, err);
  if (err->error) {
    /* Ignore field if can't convert  coordinates  */
    ObitErrClearErr (err);
    return;
  } 

  /* Also 1 cell north for differential rotation */
  pix[1] = pix[1] + 1.0;
  ObitImageDescCvtPixel (in1->myDesc, in2->myDesc, pix, pix3, err);
  if (err->error) {
    /* Ignore field if can't convert  coordinates  */
    ObitErrClearErr (err);
    return;
  } 

  /* Need to enlarge image for relative rotation, the  enlargement factor is the RA 
     pixel shift between pix2 and  pix3. */
  fact = 2.0 * fabs (pix3[0] - pix2[0]);
  fact = 1.0 + MIN (1.0, fact);

  /* Take into account relative pixel sizes */
  fact *= (fabs(in1->myDesc->cdelt[1]) / fabs(in2->myDesc->cdelt[1]));

  /* Increase size for possible differential rotation. */
  naxis[0] = (olong) (0.5+naxis[0]*fact);
  naxis[1] = (olong) (0.5+naxis[1]*fact);

  /* Fudge size a bit */
  naxis[0] = naxis[0] * 1.1;
  naxis[1] = naxis[1] * 1.1;

  /* How big is any existant output array - don't make it smaller */
  for (i=0; i<IM_MAXDIM; i++) onaxis[i] = 0; 
  if (out->image) {
    for (i=0; i<out->image->ndim; i++) onaxis[i] = out->image->naxis[i];
    if (onaxis[0] > naxis[0]) naxis[0] = onaxis[0];
    if (onaxis[1] > naxis[1]) naxis[1] = onaxis[1];
  }

  /* Center of out image */
  icx = naxis[0] / 2;
  icy = naxis[1] / 2 + 1;
  
  /* Shift reference pixel. */
  crpix[0] = icx + (crpix2[0] - pix2[0]);
  crpix[1] = icy + (crpix2[1] - pix2[1]);

  /* Round reference pixel */
  crpix[0] = (ofloat)((olong)(crpix[0]+0.5));
  crpix[1] = (ofloat)((olong)(crpix[1]+0.5));

 /* If it already exists and is big  enough, use it, create if needed */
  if (out->image) out->image = ObitFArrayRealloc (out->image, 2, naxis);
  else            out->image = ObitFArrayCreate (in1->name, 2, naxis);

  /* Copy basic in2 descriptor to out */
  out->myDesc = ObitImageDescCopy(in2->myDesc, out->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, in1->name);

  /* Modify to reflect actual output */
  for (i=0; i<IM_MAXDIM; i++) out->myDesc->crpix[i] = crpix[i];  /* reference pixel. */
  for (i=0; i<IM_MAXDIM; i++) out->myDesc->inaxes[i] = naxis[i]; /* Size */

  /* set this up for use as a memory only image */
  out->mySel->FileType = OBIT_IO_MEM;
  dim[0] = 1;
  ObitInfoListPut (out->info, "FileType", OBIT_long, dim, &out->mySel->FileType, err);
  if (err->error) Obit_traceback_msg (err, routine, in1->name);

} /* end ObitImageClone2 */

/**
 * Copy the structure of an Image to a memory resident Image
 * Only sets descriptor and data buffer members on out.
 * out->mySel->FileType set to OBIT_IO_MEM.
 * \param in   Image to be duplicated
 * \param out  An existing object pointer for output.  
 *             If the image (FArray) member exists, it is resized if needed.
 *             Descriptor updated to reflect any resizing
 * \param err Error stack, returns if not empty.
 */
void ObitImageCloneMem  (ObitImage *in, ObitImage *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "ObitImageCloneMem";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

 /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* Copy basic in descriptor to out */
  out->myDesc = ObitImageDescCopy(in->myDesc, out->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  /* Force to float pixels */
  out->myDesc->bitpix=-32;

  /* Copy Image buffer */
  if (out->image==NULL) out->image = newObitFArray (NULL);
  ObitFArrayClone(in->image, out->image, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* set this up for use as a memory only image */
  out->mySel->FileType = OBIT_IO_MEM;
  dim[0] = 1;
  ObitInfoListPut (out->info, "FileType", OBIT_long, dim, &out->mySel->FileType, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitImageCloneMem */

/**
 * Initialize structures and open file.
 * The image descriptor is read if OBIT_IO_ReadOnly or
 * OBIT_IO_ReadWrite and written to disk if opened OBIT_IO_WriteOnly.
 * Update for selection (blc,trc) is only done on OBIT_IO_ReadOnly .
 * After the file has been opened, the member image is initialized
 * for storing the image unless member extBuffer is TRUE.
 * The file etc. info should have been stored in the ObitInfoList:
 * "FileType" OBIT_long scalar = OBIT_IO_FITS or OBIT_IO_AIPS 
 *    or OBIT_IO_MEM (no persistent form )    for file type.
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite or
 *               OBIT_IO_WriteOnly).
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitImageOpen (ObitImage *in, ObitIOAccess access, 
			  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitImageDesc *desc=NULL;
  gchar *routine = "ObitImageOpen";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* ReadCal is equivalent to ReadOnly */
  if (access == OBIT_IO_ReadCal) access = OBIT_IO_ReadOnly;

  /* If the file is already open - close it  first */
  if ((in->myStatus==OBIT_Active) || (in->myStatus==OBIT_Modified)) {
    if (in->myIO) retCode = ObitIOClose (in->myIO, err);
    else retCode = OBIT_IO_OK;
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* set Status */
  in->myStatus = OBIT_Active;

  /* get selection parameters */
  /* descriptor on IO, if it exists and not WriteOnly, is more accurate 
     reflection of reality */
  if (ObitIOIsA(in->myIO) && access!=OBIT_IO_WriteOnly) 
    desc = (ObitImageDesc*)in->myIO->myDesc;
  else desc = in->myDesc;
  if (desc==NULL) desc = in->myDesc;  /* Just in case */
  ObitImageGetSelect (in->info, desc, in->mySel, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
 
  /* create appropriate ObitIO */
  ObitImageSetupIO (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Save info is actually doing IO (not Mem_only) */
  if (ObitIOIsA(in->myIO)) {
    /* Add reference to tableList */
    in->myIO->tableList = (Obit*)ObitUnref(in->myIO->tableList);
    in->myIO->tableList = (Obit*)ObitRef(in->tableList);
    
    in->myIO->access = access; /* save access type */
    
    /* most of the instructions for the I/O are in the ObitInfoList */
    retCode = ObitIOOpen (in->myIO, access, in->info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
    
    /* read or write Headers */
    if ((access == OBIT_IO_ReadOnly)|| (access == OBIT_IO_ReadWrite)) {
      /* read header info */
      retCode = ObitIOReadDescriptor(in->myIO, err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
	Obit_traceback_val (err, routine, in->name, retCode);
      
    } else if (access == OBIT_IO_WriteOnly) {
      /* Write header info */
      retCode = ObitIOWriteDescriptor(in->myIO, err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
	Obit_traceback_val (err, routine, in->name, retCode);
    }
    
    /* if ReadOnly Set descriptors for the output on in to reflect the selection
       by in->mySel,  the descriptors on in->myIO will still describe
       the external representation */
    if (access == OBIT_IO_ReadOnly) {
      ObitImageSelSetDesc ((ObitImageDesc*)in->myIO->myDesc,
			  (ObitImageSel*)in->myIO->mySel, in->myDesc, err);
      /* Update IF selection */
      if (in->myDesc->jlocif>=0) 
	ObitImageSelSetIF (in->myDesc, (ObitImageSel*)in->myIO->mySel, (ObitData*)in, err);
    } else if (access == OBIT_IO_ReadWrite) 
    /* copy actual descriptor */
      in->myDesc = ObitImageDescCopy ((ObitImageDesc*)in->myIO->myDesc, in->myDesc, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    
    /* init I/O */
    retCode = ObitIOSet (in->myIO, in->info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);

    /* Allocate buffer - resize if necessary */
    if (!in->extBuffer) {
      in->image = ObitImageSelBuffer (in->image, (ObitImageDesc*)in->myIO->myDesc, 
				      in->mySel);
    } /* end buffer allocation */
    
    /* save current location */
    in->myDesc->plane   = ((ObitImageDesc*)in->myIO->myDesc)->plane;
    in->myDesc->row     = ((ObitImageDesc*)in->myIO->myDesc)->row;
  } else retCode = OBIT_IO_OK; /* end of if IO */
  
  return retCode;
} /* end ObitImageOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitImageClose (ObitImage *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitImageClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  /* Something going on? */
  if (in->myStatus == OBIT_Inactive) return OBIT_IO_OK;
  if (in->myIO == NULL) return OBIT_IO_OK;
  retCode = OBIT_IO_OK;

  /* Any actual I/O? */
  if (in->mySel->FileType!=OBIT_IO_MEM) {
    /* flush buffer if writing */
    if (((in->myIO->access==OBIT_IO_ReadWrite) || 
	 (in->myIO->access==OBIT_IO_WriteOnly)) &&
	(in->myStatus == OBIT_Modified)) {
      retCode = ObitIOFlush (in->myIO, err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
	Obit_traceback_val (err, routine, in->name, retCode);
      
      /* Save Max/min isBlank */
      in->myDesc->maxval = 
	MAX (in->myDesc->maxval, ((ObitImageDesc*)in->myIO->myDesc)->maxval);
      in->myDesc->minval = 
	MIN (in->myDesc->minval, ((ObitImageDesc*)in->myIO->myDesc)->minval);
      if (((ObitImageDesc*)in->myIO->myDesc)->areBlanks) 
	in->myDesc->areBlanks = TRUE;
      /* Update descriptor on myIO */
      ObitImageDescCopyDesc(in->myDesc, (ObitImageDesc*)in->myIO->myDesc, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    }
    
    if (((in->myIO->access==OBIT_IO_ReadWrite) || 
	 (in->myIO->access==OBIT_IO_WriteOnly))) {
      
      /* Update descriptor on IO (the one to be written) */
      in->myIO->myDesc = ObitImageDescCopy(in->myDesc, in->myIO->myDesc, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      
      if (in->myStatus == OBIT_Modified) {
	/* Update header on disk if writing */
	retCode = OBIT_IO_OK;
	if (in->myIO->myStatus != OBIT_Inactive)
	  retCode = ObitIOWriteDescriptor(in->myIO, err);
	if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
	  Obit_traceback_val (err, routine, in->name, retCode);
      }
    }
    
    retCode = ObitIOClose (in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  } /* end of if actual I/O */

  /* set Status */
  in->myStatus = OBIT_Inactive;

  return retCode;
} /* end ObitImageClose */

/**
 * Ensures full instantiation of object - basically open to read/write header
 * and verify or create file.
 * If object has previously been opened, as demonstrated by the existance
 * of its myIO member, this operation is a no-op.
 * Virtual - calls actual class member
 * \param in     Pointer to object
 * \param exist  TRUE if object should previously exist, else FALSE
 * \param err    ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
void ObitImageFullInstantiate (ObitImage *in, gboolean exist, ObitErr *err)
{
  ObitIOAccess access;
  gchar *routine = "ObitImageFullInstantiate";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  if (in->myIO) return;  /* is this needed */

  /* Open readonly if it should exist, else writeonly */
  if (exist) access = OBIT_IO_ReadOnly;
  else access = OBIT_IO_WriteOnly;
  in->extBuffer = TRUE;  /* Don't need to assign buffer here */

  /* Open and close */
  ObitImageOpen(in, access, err);
  ObitImageClose(in, err);
  if (err->error)Obit_traceback_msg (err, routine, in->name);
  in->extBuffer = FALSE;  /* May need buffer later */
} /* end ObitImageFullInstantiate */

/**
 * Read image data from disk.
 * Reads row in->myDesc->row + 1; plane in->myDesc->plane + 1
 * A series of calls will read sequential sections of the image,
 * either a row at a time or a plane at a time as specified to
 * ObitImageOpen.  
 * The ObitImageDesc maintains the current location in the image.
 * This is a NOP if in is a memory only image
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 *             if NULL, use the image member of in.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitImageRead (ObitImage *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  ofloat *buffer = data;
  gchar *routine = "ObitImageRead";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* This is a NOP if this is a memory only image */
  if (in->mySel->FileType==OBIT_IO_MEM) return OBIT_IO_OK;

  /* check and see if its open - if not attempt */
  if ((in->myStatus!=OBIT_Active) && (in->myStatus!=OBIT_Modified)) {
    access = OBIT_IO_ReadOnly;
    retCode = ObitIOOpen (in->myIO, access, in->info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback, return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* select internal or external buffer */
  if ((buffer==NULL) && (in->image!=NULL)) buffer = in->image->array;

  /* Check buffer */
  Obit_retval_if_fail((buffer != NULL), err, retCode,
 		      "%s: No buffer allocated for %s", routine, in->name);

  retCode = ObitIORead (in->myIO, buffer, err);
  if ((retCode > OBIT_IO_EOF) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* save current location */
  in->myDesc->plane   = ((ObitImageDesc*)in->myIO->myDesc)->plane;
  in->myDesc->row     = ((ObitImageDesc*)in->myIO->myDesc)->row;

  return retCode;
} /* end ObitImageRead */

/**
 * Write information to disk.
 * Writes row in->myDesc->row + 1; plane in->myDesc->plane + 1
 * A series of calls will write sequential sections of the image,
 * either a row at a time or a plane at a time as specified to
 * ObitImageOpen.
 * The ObitImageDesc maintains the current location in the image.
 * This is a NOP if in is a memory only image
 * \param in Pointer to object to be written.
 * \param data pointer to buffer containing input data.
 *             if NULL, use the image member of in.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitImageWrite (ObitImage *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  ofloat *buffer = data;
  gchar *routine = "ObitImageWrite";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* This is a NOP if this is a memory only image */
  if (in->mySel->FileType==OBIT_IO_MEM) return OBIT_IO_OK;
  
  /* check and see if its open - if not attempt */
  if ((in->myStatus!=OBIT_Modified) && (in->myStatus!=OBIT_Active)) {
    access = OBIT_IO_WriteOnly;
    retCode = ObitIOOpen (in->myIO, access, in->info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* select internal or external buffer */
  if ((buffer==NULL)&& (in->image!=NULL)) buffer = in->image->array;

  /* Check buffer */
  Obit_retval_if_fail((buffer != NULL), err, retCode,
 		      "%s: No buffer allocated for %s", routine, in->name);

  /* most of the instructions for the I/O are in the ObitInfoList */
  retCode = ObitIOWrite (in->myIO, buffer, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* set Status */
  in->myStatus = OBIT_Modified;

  /* save current location */
  in->myDesc->plane   = ((ObitImageDesc*)in->myIO->myDesc)->plane;
  in->myDesc->row     = ((ObitImageDesc*)in->myIO->myDesc)->row;
 
 /* save max/min/blanking */
  in->myDesc->maxval    = ((ObitImageDesc*)in->myIO->myDesc)->maxval;
  in->myDesc->minval    = ((ObitImageDesc*)in->myIO->myDesc)->minval;
  in->myDesc->areBlanks = ((ObitImageDesc*)in->myIO->myDesc)->areBlanks;

  return retCode;
} /* end ObitImageWrite */

/**
 * Read the specified plane of an image from disk.
 * If the object is open on call it is returned open, otherwise closed.
 * The ObitImageDesc maintains the current location in the image.
 * This is a NOP if in is a memory only image
 * Note: the underlying routines need more work for > 3 dimensions,
 *   descriptor plane needs to be turned into an array.
 * \param in    Pointer to object to be read.
 * \param data  Pointer to buffer to write results.
 *              if NULL, use the image member of in.
 * \param plane 5 element array giving pixel numbers (1-rel) on axes 3-7
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitImageGetPlane (ObitImage *in, ofloat *data, olong plane[5], ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  gboolean doOpen = FALSE;
  olong iPlane;
  ofloat *buffer = data;
  gchar *routine = "ObitImageGetPlane";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* This is a NOP if this is a memory only image */
  if (in->mySel->FileType==OBIT_IO_MEM) return OBIT_IO_OK;

  /* check and see if its open - if not attempt */
  if ((in->myStatus!=OBIT_Active) && (in->myStatus!=OBIT_Modified)) {
    access = OBIT_IO_ReadOnly;
    doOpen = TRUE;   /* will need to close */
    retCode = ObitImageOpen (in, access, err);
    if ((retCode!=OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* which plane? */
  iPlane = PlaneNumber(plane, in->myDesc->naxis, in->myDesc->inaxes);
  in->myDesc->plane = ((ObitImageDesc*)in->myIO->myDesc)->plane  = iPlane-1;
  in->mySel->blc[2] = ((ObitImageSel*)in->myIO->mySel)->blc[2]   = iPlane;
  in->mySel->trc[2] = ((ObitImageSel*)in->myIO->mySel)->trc[2]   = iPlane;
  in->myDesc->row   = ((ObitImageDesc*)in->myIO->myDesc)->row    = 0;
  in->myDesc->IOsize= ((ObitImageDesc*)in->myIO->myDesc)->IOsize = OBIT_IO_byPlane;

  /* select internal or external buffer */
  if ((buffer==NULL) && (in->image!=NULL)) buffer = in->image->array;

  /* Check buffer */
  Obit_retval_if_fail((buffer != NULL), err, retCode,
 		      "%s: No buffer allocated for %s", routine, in->name);

  retCode = ObitIORead (in->myIO, buffer, err);
  if ((retCode > OBIT_IO_EOF) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* save current location in case next read sequential */
  in->myDesc->plane   = ((ObitImageDesc*)in->myIO->myDesc)->plane;
  in->myDesc->row     = ((ObitImageDesc*)in->myIO->myDesc)->row;

  /* Close if needed */
  if (doOpen) {
    retCode = ObitImageClose (in, err);
    if ((retCode!=OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  return retCode;
} /* end ObitImageGetPlane */

/**
 * Write the specified plane of an image from disk.
 * If the object is open on call it is returned open, otherwise closed.
 * The ObitImageDesc maintains the current location in the image.
 * This is a NOP if in is a memory only image
 * Note: the underlying routines need more work for > 3 dimensions,
 *   descriptor plane needs to be turned into an array.
 * \param in    Pointer to object to be read.
 * \param data  Pointer to buffer with pixel data
 *              if NULL, use the image member of in (MUST exist).
 * \param plane 5 element array giving pixel numbers (1-rel) on axes 3-7
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitImagePutPlane (ObitImage *in, ofloat *data, olong plane[5], ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  gboolean saveExtBuffer, doOpen = FALSE;
  olong iPlane;
  ofloat *buffer = data;
  gchar *routine = "ObitImagePutPlane";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* This is a NOP if this is a memory only image */
  if (in->mySel->FileType==OBIT_IO_MEM) return OBIT_IO_OK;

  /* check and see if its open - if not attempt */
  if ((in->myStatus!=OBIT_Active) && (in->myStatus!=OBIT_Modified)) {
    access = OBIT_IO_ReadWrite;
    doOpen = TRUE;   /* will need to close */
    saveExtBuffer = in->extBuffer;
    in->extBuffer = TRUE;
    retCode = ObitImageOpen (in, access, err);
    if ((retCode!=OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);
    in->extBuffer = saveExtBuffer;
  }

  /* which plane? */
  iPlane = PlaneNumber(plane, in->myDesc->naxis, in->myDesc->inaxes);
  in->myDesc->plane = ((ObitImageDesc*)in->myIO->myDesc)->plane  = iPlane-1;
  in->mySel->blc[2] = ((ObitImageSel*)in->myIO->mySel)->blc[2]   = iPlane;
  in->mySel->trc[2] = ((ObitImageSel*)in->myIO->mySel)->trc[2]   = iPlane;
  in->myDesc->row   = ((ObitImageDesc*)in->myIO->myDesc)->row    = 0;
  in->myDesc->IOsize= ((ObitImageDesc*)in->myIO->myDesc)->IOsize = OBIT_IO_byPlane;

  /* select internal or external buffer */
  if ((buffer==NULL) && (in->image!=NULL)) buffer = in->image->array;

  /* Check buffer */
  Obit_retval_if_fail((buffer != NULL), err, retCode,
  		      "%s: No buffer allocated for %s", routine, in->name);

  retCode = ObitIOWrite (in->myIO, buffer, err);
  if ((retCode > OBIT_IO_EOF) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* save current location in case next read sequential */
  in->myDesc->plane   = ((ObitImageDesc*)in->myIO->myDesc)->plane;
  in->myDesc->row     = ((ObitImageDesc*)in->myIO->myDesc)->row;

  /* set Status */
  in->myStatus = OBIT_Modified;

  /* save max/min/blanking */
  in->myDesc->maxval    = ((ObitImageDesc*)in->myIO->myDesc)->maxval;
  in->myDesc->minval    = ((ObitImageDesc*)in->myIO->myDesc)->minval;
  in->myDesc->areBlanks = ((ObitImageDesc*)in->myIO->myDesc)->areBlanks;

  /* Close if needed */
  if (doOpen) {
    retCode = ObitImageClose (in, err);
    in->extBuffer = FALSE;
    if ((retCode!=OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);
    /* Free image buffer if allocated */
    if (data==NULL) in->image = ObitFArrayUnref(in->image);
  }

  return retCode;
} /* end ObitImagePutPlane */

/**
 * Return a ObitTable Object to a specified table associated with
 * the input ObitImage.  
 * If such an object exists, a reference to it is returned,
 * else a new object is created and entered in the ObitTableList.
 * \param in       Pointer to object with associated tables.
 *                 This MUST have been opened before this call.
 * \param access   access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite,
 *                 or OBIT_IO_WriteOnly).
 *                 This is used to determine defaulted version number
 *                 and a different value may be used for the actual 
 *                 Open.
 * \param tabType  The table type (e.g. "AIPS CC").
 * \param tabVer   Desired version number, may be zero in which case
 *                 the highest extant version is returned for read
 *                 and the highest+1 for write.
 * \param err      ObitErr for reporting errors.
 * \return pointer to created ObitTable, NULL on failure.
 */
ObitTable* 
newObitImageTable (ObitImage *in, ObitIOAccess access, 
		   gchar *tabType, olong *tabVer, ObitErr *err)
{
  /* Call ObitData function */
  return newObitDataTable ((ObitData*)in, access, tabType, tabVer, err);

} /* end newObitImageTable */

/**
 * Destroy a specified table(s) associated with the input ObitImage.  
 * The table is removed from the ObitTableList but the external form
 * may not be updated.
 * \param in       Pointer to object with associated tables.
 * \param tabType  The table type (e.g. "AIPS CC").
 * \param tabVer   Desired version number, may be zero in which case
 *                 the highest extant version is returned for read
 *                 and the highest+1 for write.
 *                 -1 => all versions of tabType
 * \param err      ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitImageZapTable (ObitImage *in, gchar *tabType, olong tabVer, 
			      ObitErr *err)
{
  /* Call ObitData function */
  return ObitDataZapTable ((ObitData*)in, tabType, tabVer, err);
} /* end ObitImageZapTable */

/**
 * Copies the associated tables from one ObitImage to another.
 * \param in      The ObitImage with tables to copy.
 * \param out     An ObitImage to copy the tables to, old ones replaced.
 * \param exclude a NULL termimated list of table types NOT to copy.
 *                If NULL, use include
 * \param include a NULL termimated list of table types to copy.
 *                ignored if exclude nonNULL.
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitImageCopyTables (ObitImage *in, ObitImage *out, gchar **exclude,
				gchar **include, ObitErr *err)
{
  /* Call ObitData function */
  return ObitDataCopyTables ((ObitData*)in, (ObitData*)out, 
			     exclude, include, err);
} /* end ObitImageCopyTables */

/**
 * Update any disk resident structures about the current tables.
 * \param in   Pointer to object to be updated.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitImageUpdateTables (ObitImage *in, ObitErr *err)
{
  /* Call ObitData function */
  return ObitDataUpdateTables ((ObitData*)in, err);
} /* end ObitImageUpdateTables */

/** Reposition IO to beginning of file
 * \param in   Pointer to object to be rewound.
 * \param err  ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitImageIOSet (ObitImage *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  return ObitIOSet (in->myIO, in->info, err);
} /* end ObitImageIOSet */

/**
 * Set external file representation of the myBeam member of an image
 * Make same type (FITS, AIPS) as image
 * \li FITS prepend "Beam"
 * \li AIPS class = "Beam", allocate catalog slot.
 * \param image  Image whose beam name is to be set 
 * \param err     ObitErr stack for reporting problems.
 */
void ObitImageSetBeamName (ObitImage *image, ObitErr *err) 
{
  ObitIOType FileType;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar tempStr[201], *filename;
  olong i, disk, user, cno, seq;
  ObitAIPSDirCatEntry *entry;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gboolean exist;
  ObitImage *theBeam;
  gchar AName[13], AClass[7] = "Beam  ", AType[3] = "MA";
  gchar *routine = "ObitImageSetBeamName";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageIsA(image));

  theBeam = (ObitImage*)image->myBeam;

  /* Get FileType */
  ObitInfoListGet(image->info, "FileType", &type, dim, &FileType, err);
  if (err->error)  {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",
		   routine, image->name);
    return;
  }

  /* Operate by type */
  if (FileType==OBIT_IO_FITS) {        /* FITS file */
    ObitInfoListGet(image->info, "Disk", &type, dim, &disk, err);
    if (err->error)  Obit_traceback_msg (err, routine, image->name);
    
    ObitInfoListGet(image->info, "FileName", &type, dim, tempStr, err);
    if (err->error)  Obit_traceback_msg (err, routine, image->name);
    tempStr[dim[0]] = 0;  /* NULL terminate */
    
    /* form file name for file */
    filename = g_strconcat ("Beam", tempStr, NULL);

    /* setup Obit object */
    ObitImageSetFITS (theBeam, OBIT_IO_byPlane, disk, filename, blc, trc, err);
    g_free(filename);
    if (err->error)  Obit_traceback_msg (err, routine, image->name);

  } else if (FileType==OBIT_IO_AIPS) { /* AIPS file */
    ObitInfoListGet(image->info, "Disk", &type, dim, &disk, err);
    if (err->error)  Obit_traceback_msg (err, routine, image->name);
    ObitInfoListGet(image->info, "User", &type, dim, &user, err);
    if (err->error)  Obit_traceback_msg (err, routine, image->name);
    ObitInfoListGet(image->info, "CNO", &type, dim, &cno, err);
    if (err->error)  Obit_traceback_msg (err, routine, image->name);

    /* Name etc from AIPS catalog entry */
    entry = ObitAIPSDirGetEntry (disk, user, cno, err);
    if (err->error)  Obit_traceback_msg (err, routine, image->name);

    /* Allocate new cno */
    seq = entry->seq;
    for (i=0; i<12;i++) AName[i] = entry->name[i]; AName[i] = 0;
    cno = ObitAIPSDirAlloc (disk, user, AName, AClass, AType, seq, &exist, err);
    g_free(entry);
    if (err->error)  Obit_traceback_msg (err, routine, image->myBeam->name);
 
    /* setup Obit object */
    ObitImageSetAIPS (theBeam, OBIT_IO_byPlane, disk, cno, user, blc, trc, err);
    if (err->error)  Obit_traceback_msg (err, routine, image->myBeam->name);

  } else if (FileType==OBIT_IO_MEM) {  /* Memory resident only */
    ObitInfoListPut(theBeam->info, "FileType", type, dim, &FileType, err);
    if (err->error)  Obit_traceback_msg (err, routine, image->name);
  }
} /* end ObitImageSetBeamName */

/**
 * Write header keyword/value
 * \param in   object to update, must be open during call with Write access
 * \param name The label (keyword) of the information. Max 8 char
 * \param type Data type of data element (enum defined in ObitInfoList class).
 * \param dim  Dimensionality of datum.  Only scalars and strings up to 8 char are allowed
 *             Note: for strings, the first element is the length in char.
 * \param data Pointer to the data. 
 * \param err  ObitErr for reporting errors.
 */
void ObitImageWriteKeyword (ObitImage *in, 
			   gchar* name, ObitInfoType type, gint32 *dim, 
			   gconstpointer data, ObitErr *err)
{
  ObitInfoListPut(in->myDesc->info, name, type, dim, data, err);
  ObitInfoListPut(((ObitImageDesc*)in->myIO->myDesc)->info, 
		  name, type, dim, data, err);
  in->myStatus = OBIT_Modified;
} /* end ObitImageWriteKeyword */

/**
 * Read header keyword/value
 * \param in   object to update, must be fully instantiated
 * \param name [out] The label (keyword) of the information. Max 8 char
 * \param type [out] Data type of data element (enum defined in ObitInfoList class).
 * \param dim  [out] Dimensionality of datum.  Only scalars and strings up to 8 char 
 *                   are supported
 *                   Note: for strings, the first element is the length in char.
 * \param data [out] Pointer to the data. 
 * \param err  ObitErr for reporting errors.
 */
void ObitImageReadKeyword (ObitImage *in, 
			  gchar* name, ObitInfoType *type, gint32 *dim, 
			  gpointer data, ObitErr *err)
{
  ObitInfoListGet(((ObitImageDesc*)in->myIO->myDesc)->info, 
		  name, type, dim, data, err);
} /* end ObitImageReadKeyword */


/**
 * Set size of data read/write and window in image
 * sets extBuffer to FALSE.
 * \param in   object to update, must be fully instantiated
 *             Will be opened and closed to get descriptor.
 * \param IOBy Size of I/O (OBIT_IO_byPlane or OBIT_IO_byRow).
 * \param blc  Array giving bottom left corner (1-rel)
 *             if zeroes, use descriptor to get actual value
 *             for whole image.  Actual values returned
 * \param trc  Array giving top right corner (1-rel)
 *             if zeroes, read descriptor to get actual value
 *             for whole image. Actual values returned
 * \param err  ObitErr for reporting errors.
 */
void ObitImageSetSelect (ObitImage *in, ObitIOSize IOBy, 
			 olong blc[IM_MAXDIM], olong trc[IM_MAXDIM],
			 ObitErr *err)
{
   gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
   olong i;
   olong tblc[IM_MAXDIM] = {1,1,1,1,1,1,1};
   olong ttrc[IM_MAXDIM] = {0,0,0,0,0,0,0};
   gchar *routine = "ObitImageSetSelect";
 
  /*  Open image ReadOnly to get proper descriptor */
  dim[0] = IM_MAXDIM;
  for (i=0; i<IM_MAXDIM; i++) {tblc[i] = 1; ttrc[i] = 0;}
  ObitInfoListPut (in->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (in->info, "TRC", OBIT_long, dim, trc, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* If it has an IO open and close */
  if (in->myIO) {
    if ((ObitImageOpen (in, OBIT_IO_ReadOnly, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, 
		     "ERROR opening image %s", in->name);
      return;
    }
    ObitImageClose (in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Set limits/real values */
  for (i=0; i<in->myDesc->naxis; i++) {
    if (blc[i]<=0) blc[i] = 1;
    blc[i] = MIN (blc[i], in->myDesc->inaxes[i]);
    if (trc[i]<=0) trc[i] = in->myDesc->inaxes[i];
    trc[i] = MAX (trc[i], blc[i]);
    trc[i] = MIN (trc[i], in->myDesc->inaxes[i]);
  }

  /* Set selection values */
  dim[0] = 1;
  ObitInfoListAlwaysPut (in->info, "IOBy", OBIT_long, dim,  &IOBy);
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (in->info, "BLC", OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (in->info, "TRC", OBIT_long, dim, trc);

  in->extBuffer = FALSE;  /* May need buffer */
  
} /* end ObitImageSetSelect */

/**
 * Get pointer to image beam
 * \param image  Image whose beam name is to be set 
 * \param beamNo Which order Beam, only 0 for base class [0-rel]
 * \param plane  [out] Plane number for beam [all 1s here]
 * \param err    Obit error structure
 * \return pointer to beam, NULL if not defined.
 */
ObitImage* ObitImageGetBeam (ObitImage *image, olong beamNo, 
			     olong plane[5], ObitErr *err) 
{
  olong i;
  ObitImage *theBeam;

  /* error checks */
  g_assert (ObitImageIsA(image));

  for (i=0; i<5; i++) plane[i] = 1;  /* Initialize plane */

  theBeam = (ObitImage*)image->myBeam;
  return theBeam;
} /* end ObitImageGetBeam */

/**
 * Get highest order of image beam
 * Used for Sault-Wieringa wideband imaging, 0 for base class
 * \param image  Image whose beam name is to be set 
 * \return order number
 */
olong ObitImageGetBeamOrder (ObitImage *image) 
{
  return 0;
} /* end ObitImageGetBeam */

/*-------Private functions called by ObitData class ------*/
/** Private:  Copy Constructor for scratch file*/
static ObitData* newObitDataImageScratch (ObitData *in, ObitErr *err)
{
  return (ObitData*) newObitImageScratch ((ObitImage*)in, err);
} /* end newObitDataImageScratch  */

/** Private: Copy (deep) constructor.  */
static ObitData* ObitDataImageCopy  (ObitData *in, ObitData *out, 
				  ObitErr *err)
{
  return (ObitData*) ObitImageCopy ((ObitImage*)in, (ObitImage*)out, err);
} /* end  ObitDataImageCopy*/

/** Private: Copy structure */
static void ObitDataImageClone (ObitData *in, ObitData *out, ObitErr *err)
{
  ObitImageClone ((ObitImage*)in, (ObitImage*)out, err);
} /* end ObitDataImageClone */

/** Private: Zap */
static ObitData* ObitDataImageZap (ObitData *in, ObitErr *err)
{
  return (ObitData*)ObitImageZap ((ObitImage*)in, err);
} /* end ObitDataImageZap */

/** Private: Rename */
static void ObitDataImageRename (ObitData *in, ObitErr *err)
{
  ObitImageRename ((ObitImage*)in, err);
} /* end ObitDataImageRename */

/** Private: Open */
static ObitIOCode ObitDataImageOpen (ObitData *in, ObitIOAccess access, 
				  ObitErr *err)
{
  return ObitImageOpen ((ObitImage*)in, access, err);
} /* end ObitUDataImageOpen */

/** Private: Close  */
static ObitIOCode ObitDataImageClose (ObitData *in, ObitErr *err)
{
  return ObitImageClose ((ObitImage*)in, err);
} /* end ObitDataImageClose */

/** Private:  Reset IO to start of file  */
static ObitIOCode ObitDataImageIOSet (ObitData *in, ObitErr *err)
{
  return ObitImageIOSet ((ObitImage*)in, err);
} /* end  ObitDataImageIOSet */

/** Private: Assign myIO object */
static void ObitDataImageSetupIO (ObitData *in, ObitErr *err)
{
  ObitImageSetupIO ((ObitImage*)in, err);
} /* end ObitDataImageSetupIO */

/** Private: full instantiation */
static void ObitDataImageFullInstantiate (ObitData *in, gboolean exist, 
					  ObitErr *err)
{
  ObitImageFullInstantiate ((ObitImage*)in, exist, err);
} /* end ObitDataImageFullInstantiate */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitImageClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();
  myClassInfo.hasScratch    = TRUE; /* Scratch files allowed */

  /* Set function pointers */
  ObitImageClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitImageClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitImageClassInfoDefFn ( gpointer inClass)
{
  ObitImageClassInfo *theClass = (ObitImageClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitImageClassInit;
  theClass->newObit       = (newObitFP)newObitImage;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitImageClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitImageGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitImageCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitImageClear;
  theClass->ObitInit      = (ObitInitFP)ObitImageInit;
  theClass->newObitImageScratch  = 
    (newObitImageScratchFP)newObitImageScratch;
  theClass->ObitImageSame = (ObitImageSameFP)ObitImageSame;
  theClass->ObitImageGetPlane = (ObitImageGetPlaneFP)ObitImageGetPlane;
  theClass->ObitImagePutPlane = (ObitImagePutPlaneFP)ObitImagePutPlane;
  theClass->newObitImageTable = (newObitImageTableFP)newObitImageTable;
  theClass->ObitImageZapTable= (ObitImageZapTableFP)ObitImageZapTable;
  theClass->ObitImageFullInstantiate= 
    (ObitImageFullInstantiateFP)ObitImageFullInstantiate;
  theClass->ObitImageCopyTables= 
    (ObitImageCopyTablesFP)ObitImageCopyTables;
  theClass->ObitImageUpdateTables= 
    (ObitImageUpdateTablesFP)ObitImageUpdateTables;
  theClass->ObitImageSetBeamName = 
    (ObitImageSetBeamNameFP)ObitImageSetBeamName;
  theClass->ObitImageGetBeam = 
    (ObitImageGetBeamFP)ObitImageGetBeam;
  theClass->ObitImageGetBeamOrder = 
    (ObitImageGetBeamOrderFP)ObitImageGetBeamOrder;

  /* Function pointers referenced from ObitData class */
  theClass->newObitDataScratch  = (newObitDataScratchFP)newObitDataImageScratch;
  theClass->ObitDataRename  = (ObitDataRenameFP)ObitDataImageRename;
  theClass->ObitDataZap     = (ObitDataZapFP)ObitDataImageZap;
  theClass->ObitDataClone   = (ObitDataCloneFP)ObitDataImageClone;
  theClass->ObitDataCopy    = (ObitDataCopyFP)ObitDataImageCopy;
  theClass->ObitDataOpen    = (ObitDataOpenFP)ObitDataImageOpen;
  theClass->ObitDataClose   = (ObitDataCloseFP)ObitDataImageClose;
  theClass->ObitDataIOSet   = (ObitDataIOSetFP)ObitDataImageIOSet;
  theClass->ObitDataSetupIO = (ObitDataSetupIOFP)ObitDataImageSetupIO;
  theClass->ObitDataFullInstantiate= 
    (ObitDataFullInstantiateFP)ObitDataImageFullInstantiate;
  theClass->ObitDataWriteKeyword= 
    (ObitDataWriteKeywordFP)ObitImageWriteKeyword;
  theClass->ObitDataReadKeyword= 
    (ObitDataReadKeywordFP)ObitImageReadKeyword;

} /* end ObitImageClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitImageInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImage *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->myIO      = NULL;
  in->myDesc    = newObitImageDesc(in->name);
  in->mySel     = newObitImageSel(in->name);
  in->myStatus  = OBIT_Inactive;
  in->image     = NULL;
  in->extBuffer = FALSE;
  in->myGrid    = NULL;
  in->myBeam    = NULL;
  in->isScratch = FALSE;

} /* end ObitImageInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitImage* cast to an Obit*.
 */
void ObitImageClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImage *in = inn;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Delete underlying files if isScratch */
  if (in->isScratch) {
    err = newObitErr(); /* for possible messages */
    /* Remove from ObitSystem list */
    ObitSystemFreeScratch ((Obit*)in, err);
    in->isScratch = FALSE;  /* avoid infinite recursion */
    ObitImageZap (in, err); /* delete files */
    ObitErrLog(err);
    err = ObitErrUnref(err);
  }

  /* delete this class members */
  in->tableList = ObitUnref(in->tableList);
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  in->myIO      = ObitUnref(in->myIO);
  in->myDesc    = ObitUnref(in->myDesc);
  in->mySel     = ObitUnref(in->mySel);
  in->image     = ObitUnref(in->image);
  in->myGrid    = ObitUnref(in->myGrid);
  in->myBeam    = ObitUnref(in->myBeam );
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitImageClear */


/**
 * get requested information from the ObitInfoList
 * \param info Pointer to InfoList
 * \param desc  pointer to image descriptor
 * \param sel   pointer to image selector to update.
 * \param err  ObitErr for reporting errors.
 */
static void ObitImageGetSelect (ObitInfoList *info, ObitImageDesc* desc, 
				ObitImageSel *sel, ObitErr *err)
{
  ObitInfoType type;
  gint32 i, dim[MAXINFOELEMDIM];
  gchar *routine = "ObitImageGetSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(info));
  g_assert (ObitIsA(sel, ObitImageSelGetClass()));

  /* what type of underlying file? */
  if (!ObitInfoListGet(info, "FileType", &type, dim, 
		       (gpointer)&sel->FileType, err)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		"%s: Image %s incompletely defined",	routine, sel->name);
  }

  /* set defaults */
  for (i=0; i<IM_MAXDIM; i++) sel->blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) sel->trc[i] = desc->inaxes[i];

  /* get BLC, TRC */
  ObitInfoListGetTest(info, "BLC", &type, dim, sel->blc);
  /* If defaulted use 1 */
  for (i=0; i<IM_MAXDIM; i++) if (sel->blc[i]<=0) sel->blc[i] = 1;

  ObitInfoListGetTest(info, "TRC", &type, dim, sel->trc);
  /* If defaulted use dim */
  for (i=0; i<IM_MAXDIM; i++) if (sel->trc[i]<=0) sel->trc[i] = desc->inaxes[i];
} /* end ObitImageGetSelect */

/**
 * Determine plane number from pixel indices
 * \param plane  Array of 1-rel pixel indices for planes 3-7
 * \param naxis  Number of axes to check
 * \param inaxes Array of axis dimentions
 * \return the 1-rel overall plane number
 */
static olong PlaneNumber (olong plane[5], olong naxis, olong *inaxes)
{
  int i, prev, plNumber = 1;

  prev = 1;
  for (i=3; i<=naxis; i++) {
    plNumber += (plane[i-3]-1) * prev;
    prev *= inaxes[i-1];  /* product of previous number of planes */
  }
  return plNumber;
} /* end PlaneNumber */


/**
 * Create myIO object depending on value of FileType in in->info.
 * This is the principle place where the underlying file type is known.
 * \param in   Image object to attach myIO
 * \param err  ObitErr for reporting errors.
 */
static void ObitImageSetupIO (ObitImage *in, ObitErr *err)
{
  ObitIOType FileType;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar *routine = "ObitImageSetupIO";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* Get FileType */
  if (!ObitInfoListGet(in->info, "FileType", &type, dim, 
		       (gpointer)&FileType, err)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",
		   routine, in->name);
    return;
  }

  /* unlink any existing IO structure */
  in->myIO = ObitUnref (in->myIO);
  if (FileType==OBIT_IO_FITS) {    /* FITS file */
    in->myIO = (ObitIO*)newObitIOImageFITS(in->name, in->info, err);
    /* copy selector */
    ((ObitIOImageFITS*)in->myIO)->mySel = 
      ObitImageSelCopy(in->mySel, 
		       ((ObitIOImageFITS*)in->myIO)->mySel, err);
    /* copy descriptor */
    ((ObitIOImageFITS*)in->myIO)->myDesc = 
      ObitImageDescCopy(in->myDesc, 
			((ObitIOImageFITS*)in->myIO)->myDesc, err);
    
  } else if (FileType==OBIT_IO_AIPS) { /* AIPS file */
    in->myIO = (ObitIO*)newObitIOImageAIPS(in->name, in->info, err);
    /* copy selector */
    ((ObitIOImageAIPS*)in->myIO)->mySel = 
      ObitImageSelCopy(in->mySel, 
		       ((ObitIOImageAIPS*)in->myIO)->mySel, err);
    /* copy descriptor */
    ((ObitIOImageAIPS*)in->myIO)->myDesc = 
      ObitImageDescCopy(in->myDesc, 
			((ObitIOImageAIPS*)in->myIO)->myDesc, err);
  } else if (in->mySel->FileType==OBIT_IO_MEM) {  /* Memory resident only */
  }
 
} /* end ObitImageSetupIO */

