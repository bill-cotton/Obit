/* $Id$         */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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

#include "ObitOTF.h"
#include "ObitOTFCal.h"
#include "ObitIOOTFFITS.h"
#include "ObitTableOTFTarget.h"
#include "ObitTableOTFTargetUtil.h"
#include "ObitTableOTFIndex.h"
#include "ObitTableOTFFlag.h"
#include "ObitSystem.h"
#include "ObitHistory.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTF.c
 * ObitOTF class function definitions.
 * 
 * This class is derived from the #ObitData base class.
 *
 * GBT On-the-Fly single dish imaging data type class
 */

/* Documentation for doxygen main page */
/**
 *  \mainpage Obit Classes
 * Obit uses a class derivation scheme that doxygen does not understand
 * so some care is needed in interpreting this documentation.
 * Class hierarchies are generally noted in the names of modules, i.e.
 * Obit is the base class from which (almost) all others are derived.
 * Obit class derivation is by means of nested include files; each class has an
 * include file for the data members and for the class function pointers.
 * These include files include the corresponding includes of their parent class.
 * 
 * The most effective use of this documentation page is to use the "File
 * List" function and select the header file (.h file) for the given
 * class, this gives the functional interface to the class.
 * Class data members can be viewed using the links to the class
 * structure on this page or use the "Class List" function on the main
 * page and select the desired class.
 * 
 * Each class has a "ClassInfo" structure with class specific
 * information, mostly, function pointers.
 * Each instance of a class has a number of data members, including a
 * pointer to the ClassInfo structure of its class.
 */
/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitOTF";

/** Function to obtain parent ClassInfo - ObitData */
static ObitGetClassFP ObitParentGetClass = ObitDataGetClass;

/**
 * ClassInfo structure ObitOTFClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitOTFClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitOTFInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitOTFClear (gpointer in);

/** Private: Read selection parameters from ObitInfoList. */
static void ObitOTFGetSelect (ObitOTF *in, ObitInfoList *info, ObitOTFSel *sel,
			     ObitErr *err);

/* Private: Setup for calibration */
static void ObitOTFSetupCal (ObitOTF *in, ObitErr *err);

/* Private: Average a buffer of data in frequency. */
void ObitOTFAverData(ObitOTFDesc *inDesc, ObitOTFDesc *outDesc, 
		  ofloat *inBuffer, ofloat *outBuffer);

/* Private: Determine target renumbering for merging data. */
void TarRenumber(ObitOTF *in, ObitOTF *out, olong **targetRenumber, ObitErr *err);

/** Private: Assign myIO object */
static void ObitOTFSetupIO (ObitOTF *in, ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitOTFClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitOTF* newObitOTF (gchar* name)
{
  ObitOTF* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOTFClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitOTF));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitOTFInit((gpointer)out);

 return out;
} /* end newObitOTF */

/**
 * Create an OTF object with selection parameters set from an InfoList
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param inList InfoList to extract object information from
 * Following entries for FITS files ("xxx" = prefix):
 * \li xxxFileName OBIT_string  FITS file name
 * \li xxxDisk     OBIT_oint    FITS file disk number
 * \li xxxDir      OBIT_string  Directory name for xxxDisk
 *
 * For all File types types:
 * \li xxxDataType OBIT_string "UV" = UV data, "MA"=>image, "Table"=Table, 
 *                "OTF"=OTF, etc
 * \li xxxFileType OBIT_oint File type as ObitIOType, OBIT_IO_FITS, OBIT_IO_AIPS
 *
 * For xxxDataType = "OTF"
 * \li xxxnRecPIO OBIT_int (1,1,1) Number of vis. records per IO call
 * \param err     ObitErr for reporting errors.
 * \return new data object with selection parameters set
 */
ObitOTF* ObitOTFFromFileInfo (gchar *prefix, ObitInfoList *inList, 
			      ObitErr *err)
{
  ObitOTF      *out = NULL;
  ObitInfoType type;
  olong        disk, i, nrec, nThreads;
  gchar        *strTemp, inFile[129];
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gpointer     listPnt;
  gchar        *keyword=NULL, *DataTypeKey = "DataType", *DataType=NULL;
  gchar        *parm[] = {"doCalSelect", "doCalib", "gainUse", "flagVer",
			  "BChan", "EChan", "Targets", "timeRange", "Scans",
			  "Feeds", "keepCal", "replCal",
			  NULL};
  gchar *routine = "ObiOTFFromFileInfo";


  if (err->error) return out;  /* Previous error? */

  /* Create output */
  out = newObitOTF (prefix);

  /* Number of Rec per IO  */
  nrec = 1000;
  nThreads = ObitThreadNumProc (out->thread);
  nrec *= MAX (1, nThreads);
  if (prefix) keyword = g_strconcat (prefix, "nRecPIO", NULL);
  else        keyword = g_strdup("nRecPIO");
  ObitInfoListGetTest(inList, keyword, &type, dim, &nrec);
  g_free(keyword);

  /* File type - could be FITS */
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

  if (!strncmp (DataType, "FITS", 4)) {  /* FITS input */
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

    /* define object */
    ObitOTFSetFITS (out, nrec, disk, inFile, err);
    if (err->error) Obit_traceback_val (err, routine, "inList", out);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   routine, DataType);
    return out;
  }

  /* Selection/calibration */
  i = 0;
  while (parm[i]) {
    if (prefix) keyword = g_strconcat (prefix, parm[i], NULL);
    else        keyword = g_strdup(parm[i]);
    if (ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&listPnt)) {
      ObitInfoListAlwaysPut(out->info, parm[i], type, dim, (gpointer*)&listPnt);
    }
    g_free(keyword);
  } /* end loop copying parameters */
  
  /* Ensure out fully instantiated and OK */
  ObitOTFFullInstantiate (out, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "inList", out);

  return out;
} /* end ObitOTFFromFileInfo */

/**
 * Create a scratch file suitable for accepting the data to be read from in.
 * A scratch OTF is more or less the same as a normal OTF except that it is
 * automatically deleted on the final unreference.
 * The output will have the underlying files of the same type as in already 
 * allocated.
 * \param in  The object to copy
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitOTF* newObitOTFScratch (ObitOTF *in, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitOTF *out=NULL;
  ObitIOCode iretCode;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  /* Don't copy Cal and Soln tables */
  gchar *exclude[]={"OTFScanData", "OTFSoln", 
		    "OTFCal", "OTFIndex", "OTFArrayGeom", NULL};
  gchar *outName;
  olong  i, numData, NPIO;
  gchar *routine = "newObitOTFScratch";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Ensure in fully instantiated -assume OK if myIO exists */
  if (!in->myIO) ObitOTFFullInstantiate (in, TRUE, err);
  if (err->error)Obit_traceback_val (err, routine, in->name, out);

  /* Create - derive object name */
  outName = g_strconcat ("Scratch Copy: ",in->name,NULL);
  out = newObitOTF(outName);
  g_free(outName);

  /* Mark as scratch */
  out->isScratch = TRUE;

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* Output will initially have no associated tables */
  out->tableList = ObitTableListUnref(out->tableList);
  out->tableList = newObitTableList(out->name);

  /* Copy descriptor */
  out->myDesc = (gpointer)ObitOTFDescCopy(in->myDesc, out->myDesc, err); 
  out->myDesc->nrecord = 0; 

  /* number of detectors + weights */
  numData = 1;
  for (i=0; i<out->myDesc->naxis; i++) numData *= MAX (1, out->myDesc->inaxes[i]);
  out->myDesc->colRepeat[out->myDesc->ncol-1] = numData;

  /* Copy number of records per IO */
  ObitInfoListGet (in->info, "nRecPIO", &type, dim,  (gpointer)&NPIO, err);
  ObitInfoListPut (out->info, "nRecPIO", type, dim,  (gpointer)&NPIO, err);

  /* Copy geometry */
  out->geom = ObitOTFArrayGeomAver(in->geom, in->myDesc, 
				   out->geom, out->myDesc, err);
 
  /* Allocate underlying file */
  ObitSystemGetScratch (in->mySel->FileType, "OTF", out->info, err);
  if (err->error)  Obit_traceback_val (err, routine, in->name, in);
  
  /* Register in the scratch file list */
  ObitSystemAddScratch ((Obit*)out, err);

  /* same file type */
  ObitInfoListPut(out->info, "FileType", OBIT_long, dim, 
		  (gpointer*)&in->mySel->FileType, err);
  out->mySel->FileType = in->mySel->FileType;
 
  /* Fully instantiate output */
  ObitOTFFullInstantiate (out, FALSE, err);
  if (err->error)Obit_traceback_val (err, routine, out->name, out);

 /* Copy tables   */
  ObitOTFOpen(out, OBIT_IO_ReadWrite, err);
  iretCode = ObitOTFCopyTables (in, out, exclude, NULL, err);
  ObitOTFClose(out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  return out;
} /* end newObitOTFScratch */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitOTFGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOTFClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitOTFGetClass */

/**
 * Test if two ObitUVs have the same underlying structures.
 * This test is done using values entered into the #ObitInfoList
 * in case the object has not yet been opened.
 * \param in1 First object to compare
 * \param in2 Second object to compare
 * \param err ObitErr for reporting errors.
 * \return TRUE if to objects have the same underlying structures
 * else FALSE
 */
gboolean ObitOTFSame (ObitOTF *in1, ObitOTF *in2, ObitErr *err )
{
  /* Call ObitData function */
  return ObitDataSame ((ObitData*)in1, (ObitData*)in2, err);
} /* end ObitOTFSame */

/**
 * Delete underlying files and the basic object.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 * \return pointer for input object, NULL if deletion successful
 */
ObitOTF* ObitOTFZap (ObitOTF *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitOTFZap";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return in;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* Close if still active */
  if ((in->myStatus == OBIT_Active) || (in->myStatus == OBIT_Modified)){
   retCode = ObitIOClose(in->myIO, err);
   if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
     Obit_traceback_val (err, routine, in->name, in);    
  }

  /* Ensure in fully instantiated */
  ObitErrLog(err); /* Show any pending messages as they may get lost */
  ObitOTFFullInstantiate (in, TRUE, err);
  /* If this fails, clear errors and assume it doesn't exist */
  if (err->error) { 
    ObitErrClearErr(err); 
    return ObitOTFUnref(in); 
  }

  /* Delete the file (if attached) */
  ObitIOZap (in->myIO, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, in);

  /* If it's scratch remove from list */
  if (in->isScratch) ObitSystemFreeScratch ((Obit*)in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, in);

  /* Delete object */
  in->isScratch = 0; /* Already deleted underlying structures */
  while (in) in = ObitOTFUnref(in);

  return in;
} /* end ObitOTFZap */

/**
 * Rename underlying files
 * New name information depends on the underlying file type and is
 * given on the info member.
 * For FITS files:
 * \li "newFileName" OBIT_string (?,1,1) New Name of disk file.
 * \param in Pointer to object to be renamed.
 * \param err ObitErr for reporting errors.
 */
void ObitOTFRename (ObitOTF *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitOTFRename";

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
  ObitOTFFullInstantiate (in, TRUE, err);
  /* If this fails, clear errors and assume it doesn't exist */
  if (err->error) { 
    ObitErrClearErr(err); 
    return; 
  }

  /* Rename OTF */
  ObitIORename (in->myIO, in->info, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  return;
} /* end ObitOTFRename */

/**
 * Make a deep copy of input object.
 * Copies are made of complex members including disk files; these 
 * will be copied applying whatever selection is associated with the input.
 * Objects should be closed on input and will be closed on output.
 * In order for the disk file structures to be copied, the output file
 * must be sufficiently defined that it can be written.
 * The copy will be attempted but no errors will be logged until
 * both input and output have been successfully opened.
 * If the contents of the data are copied, all associated tables are 
 * copied first.
 * ObitInfoList and ObitThread members are only copied if the output object
 * didn't previously exist.
 * Parent class members are included but any derived class info is ignored.
 * The file etc. info should have been stored in the ObitInfoList:
 * \li "doCalSelect" OBIT_boolean scalar if TRUE, calibrate/select/edit input data.
 * \li  "doCalib" OBIT_int (1,1,1) >0 -> calibrate,
 * \li  "gainUse" OBIT_int (1,1,1) SN/CL table version number, 0-> use highest
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitOTF* ObitOTFCopy (ObitOTF *in, ObitOTF *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOCode iretCode, oretCode;
  gboolean oldExist, doCalSelect;
  gchar *outName;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong  i, j, numData, NPIO;
  ofloat *inBuff, *outBuff;
  ObitHistory *inHist=NULL, *outHist=NULL;
  ObitIOAccess iaccess, oaccess;
  /* Don't copy Cal and Soln tables */
  gchar *exclude[]={"OTFScanData", "OTFSoln", "OTFArrayGeom",
		    "OTFCal", "OTFIndex", "History", NULL};
  gchar *routine = "ObitOTFCopy";
 
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitOTF(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes other additions only if out newly created */
  if (!oldExist) {
    out->mySel = newObitOTFSel (out->name);
    /* Don't copy info */
    /*out->info = ObitInfoListUnref(out->info); */
    /*out->info = ObitInfoListRef(in->info); */
    /* Output will initially have no associated tables */
    out->tableList = ObitTableListUnref(out->tableList);
    out->tableList = newObitTableList(out->name);
    /* don't copy ObitOTFSel, ObitThread */
  }

  /* If the output object was created this call it cannot be fully
     defined so we're done */
  if (!oldExist) return out;

  doCalSelect = FALSE;
  ObitInfoListGetTest(in->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) iaccess = OBIT_IO_ReadCal;
  else iaccess = OBIT_IO_ReadWrite;

  /* if input has file designated, copy data */
  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitOTFOpen (in, iaccess, err);
  if ((iretCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_val (err, routine, in->name, out);

  /* Copy number of records per IO */
  NPIO = 1000; dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListGetTest (in->info, "nRecPIO", &type, dim,  &NPIO);
  ObitInfoListAlwaysPut (out->info, "nRecPIO", OBIT_long, dim,  &NPIO);

  /* Copy descriptor */
  out->myDesc = ObitOTFDescCopy(in->myDesc, out->myDesc, err);
  out->myDesc->nrecord = 0;    /* may not copy all */

  /* number of detectors + weights */
  numData = 1;
  for (i=0; i<out->myDesc->naxis; i++) numData *= MAX (1, out->myDesc->inaxes[i]);
  out->myDesc->colRepeat[out->myDesc->ncol-1] = numData;

  /* copy Array Geometry - this time with full information */
  out->geom = ObitOTFArrayGeomAver(in->geom, in->myDesc, 
				   out->geom, out->myDesc, err);

  /* test open output */
  oaccess = OBIT_IO_WriteOnly;
  oretCode = ObitOTFOpen (out, oaccess, err);
  /* If this didn't work try OBIT_IO_ReadWritey */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
     ObitErrClear(err);
     oaccess = OBIT_IO_ReadWrite;
     oretCode = ObitOTFOpen (out, oaccess, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* ObitErrClear(err); probably not an error */
    Obit_traceback_val (err, routine, in->name, out);
  }

  /* Copy any history */
  inHist  = newObitHistoryValue("in history", in->info, err);
  outHist = newObitHistoryValue("out history", out->info, err);
  outHist = ObitHistoryCopy (inHist, outHist, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  inHist  = ObitHistoryUnref(inHist);
  outHist = ObitHistoryUnref(outHist);

  /* Copy tables before data */
  iretCode = ObitOTFCopyTables (in, out, exclude, NULL, err);
  if (err->error) Obit_traceback_val (err, routine,in->name, out);

  /* Close and reopen input to init calibration which will have been 
     disturbed by the table copy */
  ObitOTFClose (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  iretCode = ObitOTFOpen (in, iaccess, err);
  if ((iretCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_val (err, routine, in->name, out);

  /* we're in business, copy */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitOTFReadSelect (in, in->buffer, err);
    else iretCode = ObitOTFRead (in, in->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;

    /* How many */
    out->myDesc->numRecBuff = in->myDesc->numRecBuff;
    numData = in->myDesc->numRecBuff;

    /* Copy records from buffer */
    inBuff  = in->buffer;
    outBuff = out->buffer;
    for (i=0; i<numData; i++) {
      for (j=0; j<out->myDesc->lrec; j++) outBuff[j] = inBuff[j];
      inBuff  += in->myDesc->lrec;
      outBuff += out->myDesc->lrec;
    }
    
    oretCode = ObitOTFWrite (out, out->buffer, err);
  }  /* End loop copying file */
  
  /* close files */
  oretCode = ObitOTFClose (out, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine,out->name, out);
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) Obit_traceback_val (err, routine,in->name, out);
  
  iretCode = ObitOTFClose (in, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine,in->name, out);
  
  return out;
} /* end ObitOTFCopy */

/**
 * Make a deep copy of input object avering frequency channels
 * Copies are made of complex members including disk files; these 
 * will be copied applying whatever selection is associated with the input.
 * Objects should be closed on input and will be closed on output.
 * In order for the disk file structures to be copied, the output file
 * must be sufficiently defined that it can be written.
 * The copy will be attempted but no errors will be logged until
 * both input and output have been successfully opened.
 * If the contents of the data are copied, all associated tables are 
 * copied first.
 * ObitInfoList and ObitThread members are only copied if the output object
 * didn't previously exist.
 * Parent class members are included but any derived class info is ignored.
 * The file etc. info should have been stored in the ObitInfoList:
 * \li "doCalSelect" OBIT_boolean scalar if TRUE, calibrate/select/edit input data.
 * \li "doCalib" OBIT_int (1,1,1) >0 -> calibrate,
 * \li "gainUse" OBIT_int (1,1,1) SN/CL table version number, 0-> use highest
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitOTF* ObitOTFAver (ObitOTF *in, ObitOTF *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOCode iretCode, oretCode;
  ObitHistory *inHist=NULL, *outHist=NULL;
  gboolean oldExist, doCalSelect;
  gchar *outName;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  /* Don't copy Cal and Soln tables */
  gchar *exclude[]={"OTFScanData", "OTFSoln", "OTFCal", "OTFFlag", "OTFIndex", 
		    "OTFArrayGeom", NULL};
  gchar *routine = "ObitOTFAver";
 
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitOTFIsA(in));
  if (out) g_assert (ObitOTFIsA(out));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitOTF(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = ((ObitClassInfo*)in->ClassInfo)->ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes other additions only if out newly created */
  if (!oldExist) {
    out->mySel = newObitOTFSel (out->name);
    /* Don't copy info */
    /* out->info = ObitInfoListUnref(out->info); */
    /* out->info = ObitInfoListRef(in->info); */
    /* Output will initially have no associated tables */
    out->tableList = ObitTableListUnref(out->tableList);
    out->tableList = newObitTableList(out->name);
    /* don't copy ObitOTFSel, ObitThread */
  }

  /* If the output object was created this call it cannot be fully
     defined so we're done */
  if (!oldExist) return out;

  doCalSelect = FALSE;
  ObitInfoListGetTest(in->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* if input has file designated, copy data */
  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitOTFOpen (in, OBIT_IO_ReadWrite, err);
  if ((iretCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_val (err, routine, in->name, out);

  /* Copy descriptor */
  out->myDesc = ObitOTFDescCopy(in->myDesc, out->myDesc, err);
  out->myDesc->nrecord = 0;    /* may not copy all */

  /* Averaging in frequency */
  out->myDesc->colRepeat[out->myDesc->ncol-1] /= out->myDesc->inaxes[out->myDesc->jlocf];
  out->myDesc->cdelt[out->myDesc->jlocf]      *= out->myDesc->inaxes[out->myDesc->jlocf];
  out->myDesc->crpix[out->myDesc->jlocf]  = 1.0;
  out->myDesc->inaxes[out->myDesc->jlocf] = 1;
  ObitOTFDescIndex (out->myDesc);

  /* copy/average Array Geometry - this time with full information */
  out->geom = ObitOTFArrayGeomAver(in->geom, in->myDesc, out->geom, out->myDesc, err);

  /* test open output */
  oretCode = ObitOTFOpen (out, OBIT_IO_WriteOnly, err);
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* ObitErrClear(err); probably not an error */
    /* unset output buffer (may be multiply deallocated) */
    out->buffer = NULL;
    out->bufferSize = 0;
    return out;
  }

  /* Copy any history */
  inHist  = newObitHistoryValue("in history", in->info, err);
  outHist = newObitHistoryValue("out history", out->info, err);
  outHist = ObitHistoryCopy (inHist, outHist, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  inHist  = ObitHistoryUnref(inHist);
  outHist = ObitHistoryUnref(outHist);

  /* Copy tables before data */
  iretCode = ObitOTFCopyTables (in, out, exclude, NULL, err);
  if (err->error) {/* add traceback,return */
    out->buffer = NULL;
    out->bufferSize = 0;
    Obit_traceback_val (err, routine,in->name, out);
  }

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  ObitOTFClose (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  iretCode = ObitOTFOpen (in, access, err);
  if ((iretCode != OBIT_IO_OK) || (err->error>0)) Obit_traceback_val (err, routine, in->name, out);

  /* we're in business, copy */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitOTFReadSelect (in, in->buffer, err);
    else iretCode = ObitOTFRead (in, in->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;

    /* Average to output buffer */
    ObitOTFAverData (in->myDesc, out->myDesc, in->buffer, out->buffer);

    /* Write it back out */
    oretCode = ObitOTFWrite (out, out->buffer, err);
  } /* end loop copying/averaging */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,in->name, out);
  
  /* close files to be sure */
  iretCode = ObitOTFClose (in, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,in->name, out);
  
  /* close files to be sure */
  oretCode = ObitOTFClose (out, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,out->name, out);
  
  return out;
} /* end ObitOTFAver */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an OTF similar to the input one.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out The output object, must be previously created and associated files
 *            specified.
 * \param err Error stack, returns if not empty.
 */
void ObitOTFClone  (ObitOTF *in, ObitOTF *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOCode iretCode=OBIT_IO_OK, oretCode=OBIT_IO_OK;
  ObitHistory *inHist=NULL, *outHist=NULL;
  gboolean doOpenClose;
  /* Don't copy Cal and Soln tables */
  gchar *exclude[]={"OTFScanData", "OTFSoln", 
		    "OTFCal", "OTFIndex", NULL};
  gchar *routine = "ObitOTFClone";
 
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes other additions  */
  out->mySel = ObitOTFSelUnref(out->mySel);
  out->mySel = newObitOTFSel (out->name);
  /* Output will initially have no associated tables */
  out->tableList = ObitTableListUnref(out->tableList);
  out->tableList = newObitTableList(out->name);
  /* don't copy ObitOTFSel, ObitThread */

  /* Open to fully instantiate input and see if it's OK */
  doOpenClose = (in->myStatus!=OBIT_Active) && (in->myStatus!=OBIT_Modified);
  if (doOpenClose)
    iretCode = ObitOTFOpen (in, OBIT_IO_ReadWrite, err);
  if ((iretCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, in->name);

  /* Copy descriptor */
  out->myDesc = ObitOTFDescCopy(in->myDesc, out->myDesc, err);
  out->myDesc->nrecord = 0;    /* No data yet */

  /* copy Array Geometry  */
  out->geom = ObitOTFArrayGeomCopy(in->geom, out->geom, err);

 /* Open Output Data */
  oretCode = ObitOTFOpen (out, OBIT_IO_WriteOnly, err) ;
  if ((oretCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, out->name);

  /* Copy any history */
  inHist  = newObitHistoryValue("in history", in->info, err);
  outHist = newObitHistoryValue("out history", out->info, err);
  outHist = ObitHistoryCopy (inHist, outHist, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  inHist  = ObitHistoryUnref(inHist);
  outHist = ObitHistoryUnref(outHist);

  /* Copy tables  */
  iretCode = ObitOTFCopyTables (in, out, exclude, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Close files */
  if (doOpenClose) ObitOTFClose (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  ObitOTFClose (out, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);

  return;
} /* end ObitOTFClone */

/**
 * Append the data in in to the end of out
 * \param in  The object to copy
 * \param out An existing object pointer for output
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitIOCode ObitOTFConcat (ObitOTF *in, ObitOTF *out, ObitErr *err)
{
  ObitIOCode iretCode, oretCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean doCalSelect;
  ofloat dayOff, scanOffset;
  olong i, inTarget, outTarget, *targetRenumber=NULL;
  gchar *routine = "ObitOTFConcat";
 
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return oretCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* Check compatability */
  Obit_retval_if_fail(((in->myDesc->ncol==out->myDesc->ncol) &&
		       (in->myDesc->naxis==out->myDesc->naxis) &&
		       (in->myDesc->numDesc==out->myDesc->numDesc) &&
		       (in->myDesc->lrec==out->myDesc->lrec)), 
		      err, oretCode, "%s: InputOTFdatasets incompatable", routine);

  /* Open input */
  doCalSelect = FALSE;
  ObitInfoListGetTest(in->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;
  iretCode = ObitOTFOpen (in, OBIT_IO_ReadCal, err);
  if ((iretCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_val (err, routine, in->name, oretCode);

  /* Get translation information 
     - highest scan in output */
  scanOffset = (ofloat)ObitOTFHighScan(out, err); 
  if (err->error) Obit_traceback_val (err, routine, in->name, oretCode);

  /* target translation table */
  TarRenumber(in, out, &targetRenumber, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, oretCode);

  /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (out->buffer) ObitIOFreeBuffer(out->buffer); /* free existing */
  out->buffer = in->buffer;
  out->bufferSize = -1;

  /* Open output */
  oretCode = ObitOTFOpen (out, OBIT_IO_ReadWrite, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    out->buffer = NULL;
    out->bufferSize = 0;
    return oretCode;
  }

  /* append to end of the file */
  out->myDesc->firstRec = out->myDesc->nrecord+1;
 
  /* Day offset between the two datasets */
  dayOff = in->myDesc->JDObs - out->myDesc->JDObs;

 /* we're in business, copy */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitOTFReadSelect (in, in->buffer, err);
    else iretCode = ObitOTFRead (in, in->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
   /* How many */
    out->myDesc->numRecBuff = in->myDesc->numRecBuff;
    /* Update data */
    for (i=0; i<in->myDesc->numRecBuff; i++) {
      /* update time */
      in->buffer[in->myDesc->iloct+i*in->myDesc->lrec] += dayOff;
      /* update scan number */
      in->buffer[in->myDesc->ilocscan+i*in->myDesc->lrec] += scanOffset;
      /* Update target number */
      inTarget  = in->buffer[in->myDesc->iloctar+i*in->myDesc->lrec];
      outTarget = targetRenumber[inTarget];
      in->buffer[in->myDesc->iloctar+i*in->myDesc->lrec] = (ofloat)outTarget;
    }
    oretCode = ObitOTFWrite (out, in->buffer, err);
  }
  
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  out->buffer = NULL;
  out->bufferSize = 0;

  if (targetRenumber) g_free (targetRenumber); /* Cleanup */

  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,in->name, oretCode);
  
  /* close files to be sure */
  iretCode = ObitOTFClose (in, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, oretCode);
  
  /* close files to be sure */
  oretCode = ObitOTFClose (out, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, out->name, oretCode);
  
  return oretCode;
} /* end ObitOTFConcat */

/**
 * Initialize structures and open file.
 * The image descriptor is read if OBIT_IO_ReadOnly, OBIT_IO_ReadCal or 
 * OBIT_IO_ReadWrite and written to disk if opened OBIT_IO_WriteOnly.
 * If access is OBIT_IO_ReadCal then the calibration/selection/editing
 * needed is initialized.
 * See the #ObitOTFSel class for a description of the selection and 
 * calibration parameters.
 * After the file has been opened the member, buffer is initialized
 * for reading/storing the data unless member bufferSize is <0.
 * The file etc. info should have been stored in the ObitInfoList:
 * \li "FileType" OBIT_long scalar = OBIT_IO_FITS for file type.
 * \li "nRecPIO" OBIT_long scalar = Maximum number of visibilities
 *               per transfer, this is the target size for Reads (may be 
 *               fewer) and is used to create buffers.
 * \li "Compress" Obit_bool scalar = TRUE indicates output is to be 
 *               in compressed format. (access=OBIT_IO_WriteOnly only).
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite,
 *               OBIT_IO_ReadCal or OBIT_IO_WriteOnly).
 *               If OBIT_IO_WriteOnly any existing data in the output file
 *               will be lost.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitOTFOpen (ObitOTF *in, ObitIOAccess access, 
			ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong need, ver;
  ObitTableOTFArrayGeom* geomTable=NULL;
  ObitTableOTFTarget* targetTable=NULL;
  ObitTableOTFIndex*  indexTable=NULL;
  gchar *routine = "ObitOTFOpen";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  
  /* Same type of access on descriptor */
  in->myDesc->access = access;
  
  /* If the file is already open - close it  first */
  if ((in->myStatus==OBIT_Active) || (in->myStatus==OBIT_Modified)) {
    if(in->myIO) retCode = ObitOTFClose (in, err);
    else retCode = OBIT_IO_OK;
   if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }
  
  /* set Status */
  in->myStatus = OBIT_Active;
      
  
  /* create appropriate ObitIO */
  ObitOTFSetupIO (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    
  /* Add reference to tableList */
  in->myIO->tableList = (Obit*)ObitUnref(in->myIO->tableList);
  in->myIO->tableList = (Obit*)ObitRef(in->tableList);
  
  /* get selection parameters */
  ObitOTFGetSelect (in, in->info, in->mySel, err);
  if (err->error) /* add traceback,return on error */
    Obit_traceback_val (err, routine, in->name, retCode);
    
  in->myIO->access = access; /* save access type */
  /*+++++++++++++++++ Actual file open ++++++++++++++++++*/
  /* most of the instructions for the I/O are in the ObitInfoList */
  retCode = ObitIOOpen (in->myIO, access, in->info, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);
  
  /* read or write Headers */
  if ((access == OBIT_IO_ReadOnly) || (access == OBIT_IO_ReadCal) || (access == OBIT_IO_ReadWrite)) {
    /* read header info, array geometry */
    retCode = ObitIOReadDescriptor(in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) Obit_traceback_val (err, routine, in->name, retCode);
    
    /* Read Array Geometry table */
    ver = 1;
    geomTable = 
      newObitTableOTFArrayGeomValue ("GeomTable", (ObitData*)in, &ver, OBIT_IO_ReadWrite, err);
    retCode = ObitOTFArrayGeomRead (&in->geom, geomTable, err);
    geomTable = ObitTableOTFArrayGeomUnref(geomTable);
    if ((retCode != OBIT_IO_OK) || (err->error)) Obit_traceback_val (err, routine, in->name, retCode);
				       
    /* Set descriptors for the output on in to reflect the selection
       by in->mySel,  the descriptors on in->myIO will still describe
       the external representation */
    ObitOTFSelSetDesc ((ObitOTFDesc*)in->myIO->myDesc,(ObitOTFSel*)in->myIO->mySel, in->myDesc, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    
    /* Output */
  } else if (access == OBIT_IO_WriteOnly) {
    /* Set descriptors for the output on in to reflect the selection
       by in->mySel,  the descriptors on in->myIO will still describe
       the external representation */
    ObitOTFSelGetDesc (in->myDesc, (ObitOTFSel*)in->myIO->mySel, (ObitOTFDesc*)in->myIO->myDesc, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    
    /* Write geometry, header info */
    ver = 1;
    geomTable = 
      newObitTableOTFArrayGeomValue ("GeomTable", (ObitData*)in, &ver, OBIT_IO_ReadWrite, err);
    retCode = ObitOTFArrayGeomWrite (in->geom, geomTable, err);
    geomTable = ObitTableOTFArrayGeomUnref(geomTable);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);

    
    /* Init target table */
    ver = 1;
    targetTable = 
      newObitTableOTFTargetValue ("TargetTable", (ObitData*)in, &ver, OBIT_IO_ReadWrite, 
				  err);
    /* Open target table */
    if (targetTable) {
      retCode = ObitTableOTFTargetOpen (targetTable, OBIT_IO_ReadWrite, err);
      if ((retCode!= OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, in->name, retCode);
      /* Close  table */
      retCode = ObitTableOTFTargetClose (targetTable, err); 
      if ((retCode!= OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, in->name, retCode);
      targetTable = ObitTableOTFTargetUnref(targetTable);
    }
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);
    
     
    /* Init index table */
    ver = 1;
    indexTable = 
      newObitTableOTFIndexValue ("Index Table", (ObitData*)in, &ver, OBIT_IO_ReadWrite, err);
    if (indexTable) {
      /* Open index table */
      retCode = ObitTableOTFIndexOpen (indexTable, OBIT_IO_ReadWrite, err);
      if ((retCode!= OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, in->name, retCode);
      /* Close  table */
      retCode = ObitTableOTFIndexClose (indexTable, err); 
      if ((retCode!= OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, in->name, retCode);
      indexTable = ObitTableOTFIndexUnref(indexTable);
    }
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);
    
    retCode = ObitIOWriteDescriptor(in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }
  
  /* initialize any Calibration needed - complete output Descriptor, Selector*/
  if (access == OBIT_IO_ReadCal) {
    ObitOTFSetupCal (in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  }
  
  /* Allocate buffer - resize if necessary */
  /* buffer size < 0 => no buffer desired */
  if (in->bufferSize >= 0) {
    need = ObitOTFSelBufferSize(in->myDesc, in->mySel);
    /* is current one big enough? */
    if ((in->buffer!=NULL) && (need>in->bufferSize)) {
      /* no - deallocate */
      if (in->buffer) ObitIOFreeBuffer(in->buffer);
      in->buffer = NULL;
      in->bufferSize = 0;
    }
    /* new one if needed */
    if (in->buffer==NULL)  
      ObitIOCreateBuffer (&in->buffer, &in->bufferSize, in->myIO, 
			  in->info, err);
  } /* end buffer allocation */
  
  /* Set I/O to beginning of the file */
  ((ObitOTFDesc*)in->myIO->myDesc)->firstRec = 0;
  /* For WriteOnly the file is truncated.*/
  if (access == OBIT_IO_WriteOnly) {
    ((ObitOTFDesc*)in->myIO->myDesc)->firstRec = 1;
    ((ObitOTFDesc*)in->myIO->myDesc)->nrecord = 0;
    in->myDesc->nrecord = 0;
  }
  
  /* save current location */
  in->myDesc->firstRec   = ((ObitOTFDesc*)in->myIO->myDesc)->firstRec;
  in->myDesc->numRecBuff = ((ObitOTFDesc*)in->myIO->myDesc)->numRecBuff;
  
  return retCode;
} /* end ObitOTFOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitOTFClose (ObitOTF *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitOTFClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  /* Something going on? */
  if (in->myStatus == OBIT_Inactive) return OBIT_IO_OK;
  if (in->myIO == NULL) return OBIT_IO_OK;

  /* flush buffer if writing */
  if (((in->myIO->access==OBIT_IO_ReadWrite) || 
       (in->myIO->access==OBIT_IO_WriteOnly)) &&
      (in->myStatus == OBIT_Modified)) {
    retCode = ObitIOFlush (in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
 
   /* Update descriptor on myIO */
    ObitOTFDescCopyDesc(in->myDesc, (ObitOTFDesc*)in->myIO->myDesc, err);
    if (err->error)
      Obit_traceback_val (err, routine, in->name, retCode);

    /* Update header on disk if writing */
    retCode = OBIT_IO_OK;
    if (in->myIO->myStatus != OBIT_Inactive)
      retCode = ObitIOWriteDescriptor(in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);    
  }

  /* Close actual file */
  retCode = ObitIOClose (in->myIO, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* shutdown any calibration */
  if (in->myIO->access==OBIT_IO_ReadCal) {
    in->myIO->myCal = ObitOTFCalShutdown((ObitOTFCal*)in->myIO->myCal, err);
    if (err->error) /* add traceback,return on error */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* Shutdown any indexing */
  ObitOTFSelShutdown (in->mySel, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* set Status */
  in->myStatus = OBIT_Inactive;

  return retCode;
} /* end ObitOTFClose */

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
void ObitOTFFullInstantiate (ObitOTF *in, gboolean exist, ObitErr *err)
{
  ObitIOAccess access;
  gchar *routine = "ObitOTFFullInstantiate";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  if (in->myIO) return;  /* is this needed */

  /* Open readonly if it should exist, else writeonly */
  if (exist) {
    access = OBIT_IO_ReadOnly;
  } else access = OBIT_IO_WriteOnly;
  in->bufferSize = -1;  /* Don't need to assign buffer here */

  /* Open and close */
  ObitOTFOpen(in, access, err);
  ObitOTFClose(in, err);
  if (err->error)Obit_traceback_msg (err, routine, in->name);
  in->bufferSize = 0;  /* May need buffer later */
} /* end ObitOTFFullInstantiate */

/**
 * Read OTF data from disk.
 * The ObitOTFDesc maintains the current location in the file.
 * The number read will be mySel->nRecPIO (until the end of the selected
 * range of records in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstRec
 * and the number of visibilities is myDesc->numRecBuff.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 *             if NULL, use the buffer member of in.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitOTFRead (ObitOTF *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  ofloat *myBuf = data;
  olong need;
  gchar *routine = "ObitOTFRead";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* If calibration/selection is requested use ObitOTFReadSelect */
  if (in->mySel->doCalSelect) return ObitOTFReadSelect (in, data, err);

 /* check and see if its open - if not ,attempt */
  if ((in->myStatus!=OBIT_Active) && (in->myStatus!=OBIT_Modified)) {
    access = OBIT_IO_ReadWrite;
    retCode = ObitIOOpen (in->myIO, access, in->info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback, return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* select internal or external buffer */
  if (myBuf==NULL) {
    myBuf = in->buffer;
    /* Check that internal buffer ( defined in gfloats) large enough */
    need = in->mySel->nRecPIO*in->myDesc->lrec;
    if (need > in->bufferSize) {
      Obit_log_error(err, OBIT_Error, 
		     "IO buffer ( %d) too small, need %d for %s", 
		     in->bufferSize, need, in->name);
      return retCode;
    }
  } 
  g_assert (myBuf != NULL); /* check it */

  retCode = ObitIORead (in->myIO, myBuf, err);
  if ((retCode > OBIT_IO_EOF) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* save current location */
  in->myDesc->firstRec   = ((ObitOTFDesc*)in->myIO->myDesc)->firstRec;
  in->myDesc->numRecBuff = ((ObitOTFDesc*)in->myIO->myDesc)->numRecBuff;

  return retCode;
} /* end ObitOTFRead */

/**
 * Read data from disk applying selection.
 * The number read will be mySel->nRecPIO (until the end of the selected
 * range of record in which case it will be smaller).
 * The first record number after a read is myDesc->firstRec
 * and the number of records is myDesc->numRecBuff.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 *             if NULL, use the buffer member of in.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitOTFReadSelect (ObitOTF *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  ofloat *myBuf = data;
  olong need;
  gboolean done;
  gchar *routine = "ObitOTFReadSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

 /* check and see if its open - if not attempt */
  if ((in->myStatus!=OBIT_Active) && (in->myStatus!=OBIT_Modified)) {
    access = OBIT_IO_ReadCal;
    retCode = ObitIOOpen (in->myIO, access, in->info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback, return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* select internal or external buffer */
  if (myBuf==NULL) {
    myBuf = in->buffer;
    /* Check that internal buffer ( defined in gfloats) large enough */
    need = in->mySel->nRecPIO*in->myDesc->lrec;
    if (need > in->bufferSize) {
      Obit_log_error(err, OBIT_Error, 
		     "IO buffer ( %d) too small, need %d for %s", 
		     in->bufferSize, need, in->name);
      return retCode;
    }
  } 
  g_assert (myBuf != NULL); /* check it */

  /* Loop until something found */
  done = FALSE;
  while (!done) {
    retCode = ObitIOReadSelect (in->myIO, myBuf, err);
    if ((retCode > OBIT_IO_EOF) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
    done = ((ObitOTFDesc*)in->myIO->myDesc)->numRecBuff>0; /* Find something? */
    done = done || (retCode==OBIT_IO_EOF);
  }

  /* save current location */
  in->myDesc->firstRec   = ((ObitOTFDesc*)in->myIO->myDesc)->firstRec;
  in->myDesc->numRecBuff = ((ObitOTFDesc*)in->myIO->myDesc)->numRecBuff;

  return retCode;
} /* end ObitOTFReadSelect */

/**
 * Write information to disk.
 * The data in the buffer will be written starting at record
 * myDesc->firstRec and the number written will be myDesc->numRecBuff
 * which should not exceed mySel->nRecPIO if the internal buffer is used.
 * myDesc->firstRec will be maintained and need not be changed for
 * sequential writing.
 * \param in Pointer to object to be written.
 * \param data pointer to buffer containing input data.
 *             if NULL, use the buffer member of in.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitOTFWrite (ObitOTF *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  ofloat *myBuf = data;
  olong need;
  gchar *routine = "ObitOTFWrite";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* If nothing to write just return OK */
  if (in->myDesc->numRecBuff<=0) return OBIT_IO_OK;

  /* check and see if its open - if not attempt */
  if ((in->myStatus!=OBIT_Modified) && (in->myStatus!=OBIT_Active)) {
    access = OBIT_IO_WriteOnly;
    retCode = ObitIOOpen (in->myIO, access, in->info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* select internal or external buffer */
  if (myBuf==NULL) {
    myBuf = in->buffer;
    /* Check that internal buffer ( defined in gfloats) large enough */
    need = in->mySel->nRecPIO*in->myDesc->lrec;
    if (need > in->bufferSize) {
      Obit_log_error(err, OBIT_Error, 
		     "IO buffer ( %d) too small, need %d for %s", 
		     in->bufferSize, need, in->name);
      return retCode;
    }
  } 
  g_assert (myBuf != NULL); /* check it */

  /* set number and location to write on myIO descriptor */
  ((ObitOTFDesc*)in->myIO->myDesc)->firstRec   = in->myDesc->firstRec;
  ((ObitOTFDesc*)in->myIO->myDesc)->numRecBuff = in->myDesc->numRecBuff;

  /* most of the instructions for the I/O are in the ObitInfoList */
  retCode = ObitIOWrite (in->myIO, myBuf, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* set Status */
  in->myStatus = OBIT_Modified;

  /* save current location */
  in->myDesc->firstRec   = ((ObitOTFDesc*)in->myIO->myDesc)->firstRec;
  in->myDesc->numRecBuff = ((ObitOTFDesc*)in->myIO->myDesc)->numRecBuff;
  in->myDesc->nrecord    = ((ObitOTFDesc*)in->myIO->myDesc)->nrecord;

  return retCode;
} /* end ObitOTFWrite */

/**
 * Return a ObitTable Object to a specified table associated with
 * the input ObitOTF.  
 * If such an object exists, a reference to it is returned,
 * else a new object is created and entered in the ObitTableList.
 * \param in       Pointer to object with associated tables.
 *                 This MUST have been opened before this call.
 * \param access   access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite,
 *                 or OBIT_IO_WriteOnly).
 *                 This is used to determine defaulted version number
 *                 and a different value may be used for the actual 
 *                 Open.
 * \param tabType  The table type (e.g. "ArrayGeom").
 * \param tabVer   Desired version number, may be zero in which case
 *                 the highest extant version is returned for read
 *                 and the highest+1 for write.
 * \param err      ObitErr for reporting errors.
 * \return pointer to created ObitTable, NULL on failure.
 */
ObitTable* 
newObitOTFTable (ObitOTF *in, ObitIOAccess access, 
		gchar *tabType, olong *tabVer, ObitErr *err)
{
  /* Call ObitData function */
  return newObitDataTable ((ObitData*)in, access, tabType, tabVer, err);
} /* end newObitOTFTable */

/**
 * Destroy a specified table(a) associated with the input ObitOTF.  
 * The table is removed from the ObitTableList
 * \param in       Pointer to object with associated tables.
 * \param tabType  The table type (e.g. "AIPS CC").
 * \param tabVer   Desired version number, may be zero in which case
 *                 the highest extant version is returned for read
 *                 and the highest+1 for write.
 *                 -1 => all versions of tabType
 * \param err      ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitOTFZapTable (ObitOTF *in, gchar *tabType, olong tabVer, 
			   ObitErr *err)
{
  /* Call ObitData function */
  return ObitDataZapTable ((ObitData*)in, tabType, tabVer, err);
} /* end ObitOTFZapTable */

/**
 * Copies the associated tables from one ObitOTF to another.
 * \param in      The ObitOTF with tables to copy.
 * \param out     An ObitOTF to copy the tables to, old ones replaced.
 * \param exclude a NULL termimated list of table types NOT to copy.
 *                If NULL, use include
 * \param include a NULL termimated list of table types to copy.
 *                ignored if exclude nonNULL.
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitOTFCopyTables (ObitOTF *in, ObitOTF *out, gchar **exclude,
			     gchar **include, ObitErr *err)
{
  /* Call ObitData function */
  return ObitDataCopyTables ((ObitData*)in, (ObitData*)out, 
			     exclude, include, err);
} /* end ObitOTFCopyTables */

/**
 * Update any disk resident structures about the current tables.
 * \param in   Pointer to object to be updated.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitOTFUpdateTables (ObitOTF *in, ObitErr *err)
{
  /* Call ObitData function */
  return ObitDataUpdateTables ((ObitData*)in, err);
} /* end ObitOTFUpdateTables */

/**
 * Reposition IO to beginning of file
 * \param in   Pointer to object to be rewound.
 * \param err  ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitOTFIOSet (ObitOTF *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  in->myDesc->firstRec   = 1;
  in->myDesc->numRecBuff = 0;
  return ObitIOSet (in->myIO, in->info, err);
} /* end ObitOTFIOSet */

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
void ObitOTFWriteKeyword (ObitOTF *in, 
			   gchar* name, ObitInfoType type, gint32 *dim, 
			   gconstpointer data, ObitErr *err)
{
   ObitInfoListPut(in->myDesc->info, name, type, dim, data, err);
 ObitInfoListPut(((ObitOTFDesc*)in->myIO->myDesc)->info, 
		  name, type, dim, data, err);
  in->myStatus = OBIT_Modified;
} /* end ObitOTFWriteKeyword */

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
void ObitOTFReadKeyword (ObitOTF *in, 
			  gchar* name, ObitInfoType *type, gint32 *dim, 
			  gpointer data, ObitErr *err)
{
  ObitInfoListGet(((ObitOTFDesc*)in->myIO->myDesc)->info, 
		  name, type, dim, data, err);
} /* end ObitOTFReadKeyword */

/**
 * Determine the number of records in the current scan.
 * Only returns valid results after reading has begun.
 * If there is no selection in effect then the size of the file is returned.
 * \param in   object to query, should be open and reads started
 * \return number of records in the current scan
 */
olong ObitOTFNumRecScan (ObitOTF *inOTF)
{
  olong nrec;

  /* Use results from the selector if valid*/
  if (inOTF->mySel->scanLastRec>0)
    nrec = inOTF->mySel->scanLastRec - inOTF->mySel->scanFirstRec + 5;
  /* else tell about all */
  else  nrec = inOTF->myDesc->nrecord;

  return nrec;
} /* end ObitOTFNumRecScan */

/**
 * Determine the highest scan number in the data set.
 * \param in   object to query
 * \param err Error stack, returns if not empty.
 * \return highest scan number.
 */
olong ObitOTFHighScan (ObitOTF *in, ObitErr *err)
{
  olong hiScan=0;
  olong scanNo, i;
  gboolean done;
  ObitIOCode retCode;
  gchar *routine = "ObitOTFHighScan";

 /* Open Input Data */
  retCode = ObitOTFOpen (in, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_val (err, routine, in->name, hiScan);

  done = ((retCode != OBIT_IO_OK) || (in->myDesc->nrecord<=0));
  while (!done) {

    /* read buffer */
    retCode = ObitOTFRead (in, NULL, err);
    if ((retCode != OBIT_IO_OK) || (err->error>0)) goto cleanup;
    done = (retCode == OBIT_IO_EOF); /* done? */
    if (done) break;

     /* Loop over buffer */
    for (i=0; i<in->myDesc->numRecBuff; i++) {
      scanNo = (olong)in->buffer[in->myDesc->ilocscan+i*in->myDesc->lrec];
      hiScan = MAX (scanNo, hiScan);
    }
  } /* end loop over file */
  
  /* close file */
 cleanup:
  retCode = ObitOTFClose (in, err);
  if ((retCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, hiScan);
  
  return hiScan;
} /* end ObitOTFNumRecScan */

/*-------Private functions called by ObitData class ------*/
/** Private:  Copy Constructor for scratch file*/
static ObitData* newObitDataOTFScratch (ObitData *in, ObitErr *err)
{
  return (ObitData*) newObitOTFScratch ((ObitOTF*)in, err);
} /* end newObitDataOTFScratch  */

/** Private: Copy (deep) constructor.  */
static ObitData* ObitDataOTFCopy  (ObitData *in, ObitData *out, 
				  ObitErr *err)
{
  return (ObitData*) ObitOTFCopy ((ObitOTF*)in, (ObitOTF*)out, err);
} /* end  ObitDataOTFCopy*/

/** Private: Copy structure */
static void ObitDataOTFClone (ObitData *in, ObitData *out, ObitErr *err)
{
  ObitOTFClone ((ObitOTF*)in, (ObitOTF*)out, err);
} /* end ObitDataOTFClone */

/** Private: Zap */
static ObitData* ObitDataOTFZap (ObitData *in, ObitErr *err)
{
  return (ObitData*)ObitOTFZap ((ObitOTF*)in, err);
} /* end ObitDataOTFZap */

/** Private: Rename */
static void ObitDataOTFRename (ObitData *in, ObitErr *err)
{
  ObitOTFRename ((ObitOTF*)in, err);
} /* end ObitDataOTFRename */

/** Private: Open */
static ObitIOCode ObitDataOTFOpen (ObitData *in, ObitIOAccess access, 
				  ObitErr *err)
{
  return ObitOTFOpen ((ObitOTF*)in, access, err);
} /* end ObitUDataOTFOpen */

/** Private: Close  */
static ObitIOCode ObitDataOTFClose (ObitData *in, ObitErr *err)
{
  return ObitOTFClose ((ObitOTF*)in, err);
} /* end ObitDataOTFClose */

/** Private:  Reset IO to start of file  */
static ObitIOCode ObitDataOTFIOSet (ObitData *in, ObitErr *err)
{
  return ObitOTFIOSet ((ObitOTF*)in, err);
} /* end  ObitDataOTFIOSet */

/** Private: Assign myIO object */
static void ObitDataOTFSetupIO (ObitData *in, ObitErr *err)
{
  ObitOTFSetupIO ((ObitOTF*)in, err);
} /* end ObitOTFSetupIO */

/** Private: full instantiation */
static void ObitDataOTFFullInstantiate (ObitData *in, gboolean exist, 
					  ObitErr *err)
{
  ObitOTFFullInstantiate ((ObitOTF*)in, exist, err);
} /* end ObitDataOTFFullInstantiate */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitOTFClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();
  myClassInfo.hasScratch  = TRUE; /* Scratch files allowed */

  /* Set function pointers */
  ObitOTFClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitOTFClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitOTFClassInfoDefFn (gpointer inClass)
{
  ObitOTFClassInfo *theClass = (ObitOTFClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitOTFClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitOTFClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitOTFGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitOTFClear;
  theClass->ObitInit      = (ObitInitFP)ObitOTFInit;
  theClass->newObit       = (newObitFP)newObitOTF;
  theClass->newObitOTFScratch  = (newObitOTFScratchFP)newObitOTFScratch;
  theClass->ObitCopy      = (ObitCopyFP)ObitOTFCopy;
  theClass->ObitClone     = NULL; /* Different call */
  theClass->newObitOTFTable= (newObitOTFTableFP)newObitOTFTable;
  theClass->ObitOTFZapTable= (ObitOTFZapTableFP)ObitOTFZapTable;
  theClass->ObitOTFFullInstantiate= 
    (ObitOTFFullInstantiateFP)ObitOTFFullInstantiate;
  theClass->ObitOTFCopyTables= 
    (ObitOTFCopyTablesFP)ObitOTFCopyTables;
  theClass->ObitOTFUpdateTables= 
    (ObitOTFUpdateTablesFP)ObitOTFUpdateTables;
  /* Function pointers referenced from ObitData class */
  theClass->newObitDataScratch  = (newObitDataScratchFP)newObitDataOTFScratch;
  theClass->ObitDataZap     = (ObitDataZapFP)ObitDataOTFZap;
  theClass->ObitDataRename  = (ObitDataRenameFP)ObitDataOTFRename;
  theClass->ObitDataClone   = (ObitDataCloneFP)ObitDataOTFClone;
  theClass->ObitDataCopy    = (ObitDataCopyFP)ObitDataOTFCopy;
  theClass->ObitDataOpen    = (ObitDataOpenFP)ObitDataOTFOpen;
  theClass->ObitDataClose   = (ObitDataCloseFP)ObitDataOTFClose;
  theClass->ObitDataIOSet   = (ObitDataIOSetFP)ObitDataOTFIOSet;
  theClass->ObitDataSetupIO = (ObitDataSetupIOFP)ObitDataOTFSetupIO;
  theClass->ObitDataFullInstantiate= 
    (ObitDataFullInstantiateFP)ObitDataOTFFullInstantiate;
  theClass->ObitDataWriteKeyword= 
    (ObitDataWriteKeywordFP)ObitOTFWriteKeyword;
  theClass->ObitDataReadKeyword= 
    (ObitDataReadKeywordFP)ObitOTFReadKeyword;


} /* end ObitOTFClassDefFn */
/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitOTFInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTF *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->myIO      = NULL;
  in->myDesc    = newObitOTFDesc(in->name);
  in->mySel     = newObitOTFSel(in->name);
  in->geom      = NULL;
  in->myStatus  = OBIT_Inactive;
  in->buffer    = NULL;
  in->bufferSize= 0;
  in->isScratch = FALSE;

} /* end ObitOTFInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitOTF* cast to an Obit*.
 */
void ObitOTFClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTF *in = inn;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Delete underlying files if isScratch */
  if (in->isScratch) {
    err = newObitErr();     /* for possible messages */
    /* Remove from ObitSystem list */
    ObitSystemFreeScratch ((Obit*)in, err);
    in->isScratch = FALSE;  /* avoid infinite recursion */
    ObitOTFZap (in, err);    /* delete files */
    ObitErrLog(err);
    err = ObitErrUnref(err);
  }

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  in->myIO      = ObitUnref(in->myIO);
  in->myDesc    = ObitOTFDescUnref(in->myDesc);
  in->mySel     = ObitOTFSelUnref(in->mySel);
  in->geom      = ObitOTFArrayGeomUnref(in->geom);
  in->tableList = ObitUnref(in->tableList);
  if (in->buffer) ObitIOFreeBuffer(in->buffer); 
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitOTFClear */

/**
 * Get requested information from the ObitInfoList
 * \param in   OTFdata.
 * \param info Pointer to InfoList
 * \param sel  pointer to uvdata selector to update.
 * \param err  ObitErr for reporting errors.
 */
static void ObitOTFGetSelect (ObitOTF *in, ObitInfoList *info, ObitOTFSel *sel,
				ObitErr *err)
{
  ObitInfoType type;
  union ObitInfoListEquiv InfoReal; 
  oint itemp;
  gint32 i, *iptr, dim[MAXINFOELEMDIM];
  olong highVer, iver;
  ofloat ftempArr[10];
  odouble *dtempArr;
  olong itempArr[10];
  ObitTableOTFTarget *TarTable=NULL;
  gchar tempStr[5], *sptr;
  gchar *routine = "ObitOTFGetSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA(info));
  g_assert (ObitOTFSelIsA(sel));

  /* what type of underlying file? */
  if (!ObitInfoListGet(info, "FileType", &type, (gint32*)&dim, 
		       (gpointer)&sel->FileType, err)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		"ObitOTFGetSelect: entry FileType not in InfoList Object %s",
		sel->name);
  }

  /* Maximum number of records per read/write? [default 100] */
  sel->nRecPIO = 100;
  ObitInfoListGetTest(info, "nRecPIO", &type, (gint32*)dim, &sel->nRecPIO);
  sel->nRecPIO = MAX (1, sel->nRecPIO); /* no fewer than 1 */

  /* Following only needed for ReadCal */
  sel->doCalSelect = FALSE;
  sel->doCal = FALSE;
  if (in->myDesc->access != OBIT_IO_ReadCal) return; 

  /* Calibrate/select/edit output? */
  sel->doCalSelect = FALSE;
  ObitInfoListGetTest(info, "doCalSelect", &type, dim, &sel->doCalSelect);

  /* Want cal-on data? */
  sel->keepCal = TRUE;
  ObitInfoListGetTest(info, "keepCal", &type, dim, &sel->keepCal);

  /* Replace data with cal? */
  sel->replCal = FALSE;
  ObitInfoListGetTest(info, "replCal", &type, dim, &sel->replCal);

   /* Selection */
  InfoReal.itg = 0;type = OBIT_oint;
  ObitInfoListGetTest(info, "BChan", &type, dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  sel->startChann  = itemp;

  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "EChan", &type, dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  if (itemp>0) sel->numberChann = itemp - sel->startChann+1;
  else  sel->numberChann = 0;

  /* Time Range */
  ftempArr[0] = -1.0e20; ftempArr[1] = 1.0e20; type = OBIT_float;
  ObitInfoListGetTest(info, "timeRange", &type, dim, &ftempArr);
  if (type==OBIT_float)
    for (i=0; i<2; i++) sel->timeRange[i] = ftempArr[i];
  else if (type==OBIT_double) {
    dtempArr = (odouble*)ftempArr;
    for (i=0; i<2; i++) sel->timeRange[i] = dtempArr[i];
  }
  /* default */
  if ((sel->timeRange[0]==0.0) && (sel->timeRange[1]==0.0)) {
    sel->timeRange[0] = -1.0e20;
    sel->timeRange[1] =  1.0e20;
  }
    
  /* Scan Range */
  itempArr[0] = -1000000000; itempArr[1] = 1000000000; 
  ObitInfoListGetTest(info, "Scans", &type, (gint32*)dim, &itempArr);
  if ((itempArr[0]==0) && (itempArr[1]==0)) {itempArr[0]=1;itempArr[1]=10000000;}
  for (i=0; i<2; i++) sel->scans[i] = itempArr[i];
  
  for (i=0; i<4; i++) tempStr[i] = ' '; tempStr[4] = 0;
  ObitInfoListGetTest(info, "Stokes", &type, (gint32*)dim, &tempStr);
  for (i=0; i<4; i++) sel->Stokes[i] = tempStr[i]; sel->Stokes[4] = 0;
  
  /*  Calibration */
  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "doCalib", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  sel->doCal = itemp > 0;
  itemp = 0;
  ObitInfoListGetTest(info, "gainUse", &type, (gint32*)dim, &itemp);
  sel->calVersion = itemp;
  /* Make sure it actually exists */
  highVer = ObitTableListGetHigh (in->tableList, "OTFCal");
  if (highVer<=0) highVer = ObitTableListGetHigh (in->tableList, "OTFSoln");
  if (highVer<=0) {
    sel->calVersion = -1;
    /*???sel->doCalSelect = FALSE;*/
    sel->doCal = FALSE;
  }
 
  /* Flagging */
  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "flagVer", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  sel->doFlag = itemp >= 0;
  sel->FGversion = itemp;

   /* Selected feeds */
  if (ObitInfoListGetP(info, "Feeds", &type, (gint32*)dim, (gpointer)&iptr)) {
    sel->numberFeedList = dim[0];
    sel->feeds = g_realloc(sel->feeds, sel->numberFeedList*sizeof(olong));
    /* loop copying */
    for (i=0; i<sel->numberFeedList; i++) {
      sel->feeds[i] = abs (iptr[i]);
    }
  } else {
    sel->numberFeedList = 0;
  }

  /* Selected targets */
  if (ObitInfoListGetP(info, "Targets", &type, (gint32*)dim, (gpointer)&sptr)) {
    sel->numberTargetList = dim[1];
    sel->targets = g_realloc(sel->targets, sel->numberTargetList*sizeof(olong));
    /* have to lookup targets - need OTFTarget table for this. */
    iver = 1;
    TarTable = newObitTableOTFTargetValue (in->name, (ObitData*)in, &iver, 
					   OBIT_IO_ReadOnly, err);
    if (TarTable==NULL) {/* No target table - only one target and it is selected */
      sel->numberTargetList = 0;
    } else { /* Lookup targets to get numbers */
      /* Do lookup */
      ObitTableOTFTargetLookup (TarTable, dim, sptr, sel->targets, err);
      if(err->error)  Obit_traceback_msg (err, routine, in->name);
      TarTable = ObitTableOTFTargetUnref(TarTable); /* release table */
    }
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } else { /* no "Targets" specified */
    sel->numberTargetList = 0; /* everything selected */
  }

} /* end ObitOTFGetSelect */

/**
 * Do various operation that have to be done in the ObitOTF class
 * rather than the ObitOTFCal class.  This mostly involves setups
 * of tables since the ObitOTF cannot be visible from the ObitOTFCal class.
 * \param in   OTF object with OTFCal to prepare for calibration
 * \param err  ObitErr for reporting errors.
 */
static void ObitOTFSetupCal (ObitOTF *in, ObitErr *err)
{
  ObitOTFSel *sel = NULL;
  ObitOTFCal *cal = NULL;
  olong highVer, useVer;
  gchar *routine = "ObitOTFSetupCal";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(in));
  g_assert (ObitIOIsA(in->myIO));
  g_assert (ObitOTFSelIsA(in->mySel));

  /* Following only needed for ReadCal */
  if (in->myDesc->access != OBIT_IO_ReadCal) return; 

  /* Create/initialize Calibrator if needed */
  if (in->myIO->myCal==NULL) in->myIO->myCal = newObitOTFCal(in->name);
  cal = (ObitOTFCal*)in->myIO->myCal;
  
  /* Setup tables as needed for calibration */
  sel = in->mySel;
  
  /* Applyig calibration? */
  if (sel->doCal) {

    cal->doSolnTable = FALSE;
    /* Use OTFCal table if one exists */
    highVer = ObitTableListGetHigh (in->tableList, "OTFCal");
    if (highVer > 0) {
      /* if calVersion ==0 use highest */
      if (sel->calVersion==0) useVer = highVer;
      else useVer = sel->calVersion;
      /* Check that requested table exists */
      Obit_return_if_fail((highVer>=useVer), err,
			  "%s: OTFCal version %d does not exist in %s",
			  routine, useVer, in->name);
      cal->CalTable =
	(Obit*) newObitTableOTFCalValue (in->name, (ObitData*)in, &useVer, 
					 OBIT_IO_ReadWrite, 0, 0, err);
      if (cal->CalTable==NULL) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: NO calibration table %d for %s", routine, useVer, 
		       in->name);
	return;
      }
      /* Tell what you're about to do */
      Obit_log_error(err, OBIT_InfoErr,
		     "Applying OTFCal table version %d",  useVer);
   } else {

      /* No Cal table - Use Soln table */
      highVer = ObitTableListGetHigh (in->tableList, "OTFSoln");
      if (sel->calVersion==0) useVer = highVer;
      else useVer = sel->calVersion;
      /* Check that requested table exists */
      Obit_return_if_fail((highVer>=useVer), err,
			"%s: OTFSoln version %d does not exist in %s",
			routine, useVer, in->name);
      cal->doSolnTable = TRUE;
      cal->SolnTable =
	(Obit*) newObitTableOTFSolnValue (in->name, (ObitData*)in, &useVer, 
					  OBIT_IO_ReadWrite, 0, 0, err);
      if (cal->SolnTable==NULL) { /* Couldn't open table */
	Obit_log_error(err, OBIT_Error, 
		       "%s: NO calibration Soln table %d for %s", routine, 
		       useVer, in->name);
	return;
      }
      /* Tell what you're about to do */
      Obit_log_error(err, OBIT_InfoErr,
		     "Applying OTFSoln table version %d",  useVer);
    } /* end get cal table */
  }
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Flag table for Flagging */
  if (sel->doFlag) {
    /* if sel->FGversion ==0 use highest */
    highVer = ObitTableListGetHigh (in->tableList, "OTFFlag");
    if (sel->FGversion==0) useVer = highVer;
    else useVer = MIN (sel->FGversion, highVer);

    if (useVer>0) { /* have one - use */
      cal->FlagTable = 
	(Obit*)newObitTableOTFFlagValue (in->name, (ObitData*)in, &useVer, 
				    OBIT_IO_ReadWrite, err);
      if (cal->FlagTable==NULL) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: NO Flagging table %d for %s", routine, useVer, 
		       in->name);
	return;
      }
    } else {
      /* no flag table - ignore */
      sel->doFlag = FALSE;
      cal->FlagTable = NULL;
    }
  } /* end FG table setup */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

   /* Initialize indexing of the uv data if an index table exists */
  highVer = ObitTableListGetHigh (in->tableList, "OTFIndex");
  useVer = 1;
  sel->doIndex = FALSE;
  if (highVer>=1) {
    sel->IndexTable =
      (Obit*) newObitTableOTFIndexValue (in->name, (ObitData*)in, &useVer, 
					 OBIT_IO_ReadWrite, err);
    /* Check that it actually contains records */
    ObitTableOTFIndexOpen ((ObitTableOTFIndex*)sel->IndexTable, OBIT_IO_ReadWrite, err);
    ObitTableOTFIndexClose ((ObitTableOTFIndex*)sel->IndexTable, err); 
    if (((ObitTableOTFIndex*)(sel->IndexTable))->myDesc->nrow<=0) 
      sel->IndexTable = ObitTableOTFIndexUnref(sel->IndexTable);
 } else sel->IndexTable = NULL;
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  if ((highVer>=1) && (sel->IndexTable))
    ObitOTFSelNextInit (sel, (ObitOTFDesc*)in->myIO->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Start up calibration - finish output Descriptor, and Selector */
  ObitOTFCalStart ((ObitOTFCal*)in->myIO->myCal, in->mySel, 
		   (ObitOTFDesc*)in->myIO->myDesc, in->geom, in->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitOTFSetupCal */

/**
 * Average a buffer of data in frequency
 * \param inDesc     Input descriptor
 * \param outDesc    Output Descriptor
 * \param inBuffer   Input buffer
 * \param outBuffer  Output buffer
 */
void ObitOTFAverData(ObitOTFDesc *inDesc, ObitOTFDesc *outDesc, 
		  ofloat *inBuffer, ofloat *outBuffer)
{
  ofloat *inData, *outData, *inDes, *outDes, sum, sumWt, fblank=ObitMagicF();
  olong i, j, count, nstoke, nchan, nfeed;
  olong istoke, ifeed, ichan, isoff, ifoff, icoff, osoff, ofoff, ocoff;
  gboolean doDataWt;
  olong incdatawt;

  /* How many records to average? */
  outDesc->numRecBuff = inDesc->numRecBuff;
  nstoke = inDesc->inaxes[inDesc->jlocs];
  nchan  = inDesc->inaxes[inDesc->jlocf];
  nfeed  = inDesc->inaxes[inDesc->jlocfeed];
  incdatawt = inDesc->incdatawt; /* increment in data-wt axis */
  doDataWt = incdatawt>1;   /* Have Data-Wt axis? */

  inData  = inBuffer  + inDesc->ilocdata;
  outData = outBuffer + outDesc->ilocdata;
  inDes   = inBuffer;
  outDes  = outBuffer;

  /* Loop over buffer */
  for (i=0; i<inDesc->numRecBuff; i++) {

    /* Copy descriptive information */
    for (j=0; j<inDesc->ilocdata; j++) outDes[j] = inDes[j];

    /* loop over data - Stokes outer loop */
    for (istoke=0; istoke<nstoke; istoke++) {
      isoff = istoke * inDesc->incs; /* offset in data */
      osoff = istoke * outDesc->incs;
	
	/* Loop over feeds */
	for (ifeed=0; ifeed<nfeed; ifeed++) {
	  ifoff = isoff + ifeed * inDesc->incfeed;
	  ofoff = osoff + ifeed * outDesc->incfeed;
	    
	  /* Channel loop */
	  sum   = 0.0;
	  sumWt = 0.0;
	  count = 0;
	  for (ichan=0; ichan<nchan; ichan++) {
	    icoff = ifoff + ichan * inDesc->incf;
	    if (inData[icoff] != fblank) {
	      count++;
	      sum += inData[icoff];
	      /* Weights in data? */
	      if (doDataWt && (inData[icoff+1] > 0.0)) 
		sumWt += inData[icoff+1];
	    }
	  } /* end loop summing */
	  
	  /* Set average in output */
	  ocoff = ofoff;  /* Only one frequency */
	  if (count>0) {
	    outData[ofoff] = sum/count;
	    if (doDataWt) outData[ofoff+1] = sumWt/count;
	    } else {
	    outData[ocoff]   = fblank;
	    if (doDataWt) outData[ocoff+1] = 0.0;
	    }
	      
	} /* end Feed loop */
    } /* end Stokes loop */
    
    /* Update data pointers */
    inData  += inDesc->lrec;
    outData += outDesc->lrec;
    inDes   += inDesc->lrec;
    outDes  += outDesc->lrec;
  } /* end loop over buffer */
} /* end ObitOTFAverData */

/**
 * Return a table of target renumbering values
 * \param in             Input OTF
 * \param out            Output OTF
 * \param targetRenumber [out] table of renumbering values,
 *                       index of input target value gives output
 *                       Should be g_freeed when done
 * \param err            for reporting errors.
 */
void TarRenumber(ObitOTF *in, ObitOTF *out, olong **targetRenumber, ObitErr *err)
{
  ObitSourceList* SList=NULL;
  ObitTableOTFTarget* table;
  olong ver, i, maxid;
  gchar *routine = "ObitOTF:TarRenumber";

  /* Get source info from input data */
  /* create output Target table object */
  ver = 1;
  table = newObitTableOTFTargetValue ("Target table", (ObitData*)in, &ver, 
				      OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  SList = ObitTableOTFTargetGetList (table, err);
  if (err->error) goto cleanup;
  if (table) ObitTableOTFTargetUnref(table);

  if (SList->number<=0) goto cleanup; /* Anything to do */

  /* Find maximum ID */
  maxid = 0;
  for (i=0; i<SList->number; i++) {
    maxid = MAX (SList->SUlist[i]->SourID, maxid);
  }
  
  /* Create output array */
  *targetRenumber = g_malloc0((1+maxid)*sizeof(olong));

  /* create output Target table object */
  ver = 1;
  table = newObitTableOTFTargetValue ("Target table", (ObitData*)out, &ver, 
				      OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;

  /* Loop over input */
  for (i=0; i<SList->number; i++) {
    /* Update/get source renumber */
    (*targetRenumber)[SList->SUlist[i]->SourID] = 
      ObitTableOTFTargetGetAddSource(table, SList->SUlist[i]->SourceName,SList->SUlist[i]->Qual,
				     SList->SUlist[i]->RAMean, SList->SUlist[i]->DecMean,
				     SList->SUlist[i]->equinox, err);
   if (err->error) goto cleanup;
 }

  /* Cleanup */
 cleanup:
  if (SList) ObitSourceListUnref(SList);
  if (table) ObitTableOTFTargetUnref(table);
  if(err->error) Obit_traceback_msg (err, routine, in->name);

}  /* end TarRenumber */

/**
 * Create myIO object depending on value of FileType in in->info.
 * This is the principle place where the underlying file type is known.
 * \param in   OTF object to attach myIO
 * \param err  ObitErr for reporting errors.
 */
static void ObitOTFSetupIO (ObitOTF *in, ObitErr *err)
{
  ObitIOType FileType;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar *routine = "ObitOTFSetupIO";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* Get FileType */
  if (!ObitInfoListGet(in->info, "FileType", &type, (gint32*)&dim, 
		       (gpointer)&FileType, err)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",
		   routine, in->name);
    return;
  }

  /* unlink any existing IO structure */
  in->myIO = ObitUnref (in->myIO);
  if (FileType==OBIT_IO_FITS) {  /* FITS files */
    in->myIO = (ObitIO*)newObitIOOTFFITS(in->name, in->info, err);
    /* copy selector pointer - use same object for OTF and IO */
    ((ObitIOOTFFITS*)in->myIO)->mySel = ObitUnref(((ObitIOOTFFITS*)in->myIO)->mySel);
    ((ObitIOOTFFITS*)in->myIO)->mySel = ObitRef(in->mySel);
    /* copy descriptor */
    ((ObitIOOTFFITS*)in->myIO)->myDesc = ObitOTFDescCopy(in->myDesc, 
		     ((ObitIOOTFFITS*)in->myIO)->myDesc, err);

  } else if (FileType==OBIT_IO_AIPS) {  /* Oops - NO AIPS allowed */
    g_error("%s: NO AIPS support for OTF data",routine);
  }

} /* end ObitOTFSetupIO */


