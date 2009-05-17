/* $Id$          */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2009                                          */
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

#include "ObitUV.h"
#include "ObitIOUVFITS.h"
#include "ObitIOUVAIPS.h"
#include "ObitAIPSDir.h"
#include "ObitFITS.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"
#include "ObitTableFQ.h"
#include "ObitTableANUtil.h"
#include "ObitTableBP.h"
#include "ObitTableBL.h"
#include "ObitTableCL.h"
#include "ObitTableCQ.h"
#include "ObitTableSN.h"
#include "ObitTableFG.h"
#include "ObitTableNX.h"
#include "ObitTableSNUtil.h"
#include "ObitTableCLUtil.h"
#include "ObitTableFQUtil.h"
#include "ObitTableSUUtil.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitHistory.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUV.c
 * ObitUV class function definitions.
 * This class is derived from the ObitData base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUV";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitDataGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitUVClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVClear (gpointer in);

/** Private: Read selection parameters from ObitInfoList. */
static void ObitUVGetSelect (ObitUV *in, ObitInfoList *info, ObitUVSel *sel,
			     ObitErr *err);

/** Private: Setup for calibration */
static void ObitUVSetupCal (ObitUV *in, ObitErr *err);

/** Private: Assign myIO object */
static void ObitUVSetupIO (ObitUV *in, ObitErr *err);

/** Private: Copy tables with selection */
static ObitIOCode CopyTablesSelect (ObitUV *inUV, ObitUV *outUV, ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitUVClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUV* newObitUV (gchar* name)
{
  ObitUV* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUV));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVInit((gpointer)out);

 return out;
} /* end newObitUV */

/**
 * Create an UV object with selection parameters set from an InfoList
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param inList InfoList to extract object information from
 * Following InfoList entries for AIPS files ("xxx" = prefix):
 * \li xxxName  OBIT_string  AIPS file name
 * \li xxxClass OBIT_string  AIPS file class
 * \li xxxDisk  OBIT_oint    AIPS file disk number
 * \li xxxSeq   OBIT_oint    AIPS file Sequence number
 * \li xxxUser  OBIT_oint    AIPS User number
 * \li xxxCNO   OBIT_oint    AIPS Catalog slot number
 * \li xxxDir   OBIT_string  Directory name for xxxDisk
 *
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
 * For xxxDataType = "UV"
 * \li xxxnVisPIO  OBIT_int (1,1,1) Number of vis. records per IO call
 * \li xxxdoCalSelect OBIT_bool (1,1,1) Select/calibrate/edit data?
 * \li xxxStokes OBIT_string (4,1,1) Selected output Stokes parameters:
 *              "    "=> no translation,"I   ","V   ","Q   ", "U   ", 
 *              "IQU ", "IQUV",  "IV  ", "RR  ", "LL  ", "RL  ", "LR  ", 
 *              "HALF" = RR,LL, "FULL"=RR,LL,RL,LR. [default "    "]
 *              In the above 'F' can substitute for "formal" 'I' (both RR+LL).
 * \li xxxBChan OBIT_int (1,1,1) First spectral channel selected. [def all]
 * \li xxxEChan OBIT_int (1,1,1) Highest spectral channel selected. [def all]
 * \li xxxBIF   OBIT_int (1,1,1) First "IF" selected. [def all]
 * \li xxxEIF   OBIT_int (1,1,1) Highest "IF" selected. [def all]
 * \li xxxdoPol OBIT_int (1,1,1) >0 -> calibrate polarization.
 * \li xxxdoCalib OBIT_int (1,1,1) >0 -> calibrate, 2=> also calibrate Weights
 * \li xxxgainUse OBIT_int (1,1,1) SN/CL table version number, 0-> use highest
 * \li xxxflagVer OBIT_int (1,1,1) Flag table version, 0-> use highest, <0-> none
 * \li xxxBLVer   OBIT_int (1,1,1) BL table version, 0> use highest, <0-> none
 * \li xxxBPVer   OBIT_int (1,1,1) Band pass (BP) table version, 0-> use highest
 * \li xxxSubarray OBIT_int (1,1,1) Selected subarray, <=0->all [default all]
 * \li xxxdropSubA OBIT_bool (1,1,1) Drop subarray info?
 * \li xxxFreqID   OBIT_int (1,1,1) Selected Frequency ID, <=0->all [default all]
 * \li xxxtimeRange OBIT_float (2,1,1) Selected timerange in days.
 * \li xxxUVRange  OBIT_float (2,1,1) Selected UV range in kilowavelengths.
 * \li xxxInputAvgTime OBIT_float (1,1,1) Input data averaging time (sec).
 *               used for fringe rate decorrelation correction.
 * \li xxxSources OBIT_string (?,?,1) Source names selected unless any starts with
 *               a '-' in which case all are deselected (with '-' stripped).
 * \li xxxsouCode OBIT_string (4,1,1) Source Cal code desired, '    ' => any code selected
 *                                  '*   ' => any non blank code (calibrators only)
 *                                  '-CAL' => blank codes only (no calibrators)
 * \li xxxQual    Obit_int (1,1,1)  Source qualifier, -1 [default] = any
 * \li xxxAntennas OBIT_int (?,1,1) a list of selected antenna numbers, if any is negative
 *                 then the absolute values are used and the specified antennas are deselected.
 * \li xxxcorrType OBIT_int (1,1,1) Correlation type, 0=cross corr only, 1=both, 2=auto only.
 * \li xxxpassAl l OBIT_bool (1,1,1) If True, pass along all data when selecting/calibration
 *                                 even if it's all flagged, 
 *                                 data deselected by time, source, antenna etc. is not passed.
 * \li xxxdoBand  OBIT_int (1,1,1) Band pass application type <0-> none
 *     (1) if = 1 then all the bandpass data for each antenna
 *         will be averaged to form a composite bandpass
 *         spectrum, this will then be used to correct the data.
 *     (2) if = 2 the bandpass spectra nearest in time (in a weighted
 *         sense) to the uv data point will be used to correct the data.
 *     (3) if = 3 the bandpass data will be interpolated in time using
 *         the solution weights to form a composite bandpass spectrum,
 *         this interpolated spectrum will then be used to correct the
 *         data.
 *     (4) if = 4 the bandpass spectra nearest in time (neglecting
 *         weights) to the uv data point will be used to correct the
 *         data.
 *     (5) if = 5 the bandpass data will be interpolated in time ignoring
 *         weights to form a composite bandpass spectrum, this
 *         interpolated spectrum will then be used to correct the data.
 * \li xxxSmooth  OBIT_float (3,1,1) specifies the type of spectral smoothing
 *        Smooth(1) = type of smoothing to apply:
 *           0 => no smoothing
 *           1 => Hanning
 *           2 => Gaussian
 *           3 => Boxcar
 *           4 => Sinc (i.e. sin(x)/x)
 *         Smooth(2) = the "diameter" of the function, i.e.
 *           width between first nulls of Hanning triangle
 *           and sinc function, FWHM of Gaussian, width of
 *           Boxcar. Defaults (if < 0.1) are 4, 2, 2 and 3
 *           channels for Smooth(1) = 1 - 4.
 *         Smooth(3) = the diameter over which the convolving
 *           function has value - in channels.
 *           Defaults: 1, 3, 1, 4 times Smooth(2) used when
 * \li xxxSubScanTime Obit_float scalar [Optional] if given, this is the 
 *          desired time (days) of a sub scan.  This is used by the 
 *          selector to suggest a value close to this which will
 *          evenly divide the current scan.  See #ObitUVSelSubScan
 *          0 => Use scan average.
 *          This is only useful for ReadSelect operations on indexed ObitUVs.
 * \param err     ObitErr for reporting errors.
 * \return new data object with selection parameters set
 */
ObitUV* ObitUVFromFileInfo (gchar *prefix, ObitInfoList *inList, 
			    ObitErr *err)
{
  ObitUV       *out = NULL;
  ObitInfoType type;
  olong        Aseq, AIPSuser, disk, cno, i, nvis, nThreads;
  gchar        *strTemp, inFile[129], stemp[256];
  gchar        Aname[13], Aclass[7], *Atype = "UV";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gpointer     listPnt;
  gchar        *keyword=NULL, *DataTypeKey = "DataType", *DataType=NULL;
  gchar        *parm[] = {"DoCalSelect", "Stokes", "BChan", "EChan", "BIF", "EIF",
			  "doPol", "doCalib", "gainUse", "flagVer", "BLVer", "BPVer",
			  "Subarray", "dropSubA", "FreqID", "timeRange", "UVRange",
			  "InputAvgTime", "Sources", "souCode", "Qual", "Antennas",
			  "corrType", "passAll", "doBand", "Smooth", "SubScanTime",
			  NULL};
  gchar *routine = "ObiUVFromFileInfo";

  if (err->error) return out;  /* Previous error? */

  /* Create output */
  out = newObitUV (prefix);

  /* Number of Vis per IO  */
  nvis = 1000;
  nThreads = ObitThreadNumProc (out->thread);
  nvis *= MAX (1, nThreads);
  if (prefix) keyword = g_strconcat (prefix, "nVisPIO", NULL);
  else        keyword = g_strdup("nVisPIO");
  ObitInfoListGetTest(inList, keyword, &type, dim, &nvis);
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
  if (keyword) g_free(keyword); keyword = NULL;

  if (!strncmp (DataType, "AIPS", 4)) { /* AIPS */
    /* AIPS disk */
    if (prefix) keyword = g_strconcat (prefix, "Disk", NULL);
    else        keyword = g_strdup("Disk");
    ObitInfoListGet(inList, keyword, &type, dim, &disk, err);
    if (keyword) g_free(keyword); keyword = NULL;

    /* If prefixDir given, lookup disk number */
    if (prefix) keyword = g_strconcat (prefix, "Dir", NULL);
    else        keyword = g_strdup("Dir");
    if (ObitInfoListGetP(inList, keyword, &type, dim, (gpointer)&strTemp)) {
      /* make sure NULL terminated */
      strncpy (stemp, strTemp, MIN(255,dim[0])); stemp[MIN(255,dim[0])] = 0;
      disk = ObitAIPSFindDirname (stemp, err);
      if (err->error) Obit_traceback_val (err, routine, routine, out);
    }
    if (keyword) g_free(keyword); keyword = NULL;

    /* AIPS name */
    if (prefix) keyword = g_strconcat (prefix, "Name", NULL);
    else        keyword = g_strdup("Name");
    if (ObitInfoListGetP(inList, keyword, &type, dim, (gpointer)&strTemp)) {
      strncpy (Aname, strTemp, 13);
    } else { /* Didn't find */
      strncpy (Aname, "No Name ", 13);
    } 
    Aname[12] = 0;
    if (keyword) g_free(keyword); keyword = NULL;

    /* AIPS class */
    if (prefix) keyword = g_strconcat (prefix, "Class", NULL);
    else        keyword = g_strdup("Class");
    if  (ObitInfoListGetP(inList, keyword, &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    Aclass[6] = 0;
    if (keyword) g_free(keyword); keyword = NULL;

    /* input AIPS sequence */
    if (prefix) keyword = g_strconcat (prefix, "Seq", NULL);
    else        keyword = g_strdup("Seq");
    ObitInfoListGet(inList, keyword, &type, dim, &Aseq, err);
    if (keyword) g_free(keyword); keyword = NULL;

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

    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "inList", out);
    
    /* define object */
    ObitUVSetAIPS (out, nvis, disk, cno, AIPSuser, err);
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
    if (keyword) g_free(keyword); keyword = NULL;
    
    /* input FITS disk */
    if (prefix) keyword = g_strconcat (prefix, "Disk", NULL);
    else        keyword = g_strdup("Disk");
    ObitInfoListGet(inList, keyword, &type, dim, &disk, err);
    if (keyword) g_free(keyword); keyword = NULL;

    /* If prefixDir given, lookup disk number */
    if (prefix) keyword = g_strconcat (prefix, "Dir", NULL);
    else        keyword = g_strdup("Dir");
    if (ObitInfoListGetP(inList, keyword, &type, dim, (gpointer)&strTemp)) {
      /* make sure NULL terminated */
      strncpy (stemp, strTemp, MIN(255,dim[0])); stemp[MIN(255,dim[0])] = 0;
      disk = ObitFITSFindDir (stemp, err);
      if (err->error) Obit_traceback_val (err, routine, routine, out);
    }
    if (keyword) g_free(keyword); keyword = NULL;

    /* define object */
    ObitUVSetFITS (out, nvis, disk, inFile, err);
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
    i++;
  } /* end loop copying parameters */

  /* Copy any InfoList Parameters */
  if (prefix) keyword = g_strconcat (prefix, "Info", NULL);
  else        keyword = g_strdup("Info");
  ObitInfoListCopyWithPrefix (inList, out->info, keyword, TRUE);
  
  /* Ensure out fully instantiated and OK */
  ObitUVFullInstantiate (out, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "inList", out);

  return out;
} /* end ObitUVFromFileInfo */

/**
 * Create a scratch file suitable for accepting the data to be read from in.
 * A scratch UV is more or less the same as a normal UV except that it is
 * automatically deleted on the final unreference.
 * The output will have the underlying files of the same type as in already 
 * allocated.
 * The object is defined but the underlying structures are not created.
 * \param in  The object to copy
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitUV* newObitUVScratch (ObitUV *in, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitUV *out=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *outName;
  olong NPIO;
  /* Don't copy Cal and Soln tables */
  /*gchar *exclude[]={"AIPS UV", "AIPS CL", "AIPS NI",
    "AIPS SN", "AIPS NX", "AIPS HI", "AIPS PL", "AIPS SL", NULL};*/
  gchar *routine = "newObitUVScratch";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Ensure in fully instantiated -assume OK if myIO exists */
  if (!in->myIO) ObitUVFullInstantiate (in, TRUE, err);
  if (err->error)Obit_traceback_val (err, routine, in->name, out);

  /* Create - derive object name */
  outName = g_strconcat ("Scratch Copy: ",in->name,NULL);
  out = newObitUV(outName);
  g_free(outName);

  /* Mark as scratch */
  out->isScratch = TRUE;

   /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

 /* Copy descriptor */
  out->myDesc = (gpointer)ObitUVDescCopy(in->myDesc, out->myDesc, err);
  out->myDesc->nvis = 0;  /* no data yet */
 
   /* Copy number of records per IO */
  ObitInfoListGet (in->info, "nVisPIO", &type, dim,  (gpointer)&NPIO, err);
  ObitInfoListPut (out->info, "nVisPIO", type, dim,  (gpointer)&NPIO, err);

  /* Allocate underlying file */
  ObitSystemGetScratch (in->mySel->FileType, "UV", out->info, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, in);
  
  /* Register in the scratch file list */
  ObitSystemAddScratch ((Obit*)out, err);

  /* File may not be completely properly defined here, defer instantiation */
 return out;
} /* end newObitUVScratch */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVGetClass */

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
gboolean ObitUVSame (ObitUV *in1, ObitUV *in2, ObitErr *err )
{
  /* Call ObitData function */
  return ObitDataSame ((ObitData*)in1, (ObitData*)in2, err);
} /* end ObitUVSame */

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
void ObitUVRename (ObitUV *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitUVRename";

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
  ObitUVFullInstantiate (in, TRUE, err);
  /* If this fails, clear errors and assume it doesn't exist */
  if (err->error) { 
    ObitErrClearErr(err); 
    return; 
  }

  /* Rename UV */
  ObitIORename (in->myIO, in->info, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  return;
} /* end ObitUVRename */

/**
 * Delete underlying files and the basic object.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 * \return pointer for input object, NULL if deletion successful
 */
ObitUV* ObitUVZap (ObitUV *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitUVZap";

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
  ObitUVFullInstantiate (in, TRUE, err);
  /* If this fails, clear errors and assume it doesn't exist */
  if (err->error) { 
    ObitErrClearErr(err); 
    return ObitUVUnref(in); 
  }

  /* Delete the file and all tables */
  ObitIOZap (in->myIO, err);
  if (err->error)  Obit_traceback_val (err, routine, in->name, in);

  /* If it's scratch remove from list */
  if (in->isScratch) ObitSystemFreeScratch ((Obit*)in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, in);

  /* Delete object */
  in->isScratch = 0; /* Already deleted underlying structures */
  while (in) in = ObitUVUnref(in); 

  return in;
} /* end ObitUVZap */

/**
 * Make a deep copy of input object.
 * Copies are made of complex members including disk files; these 
 * will be copied applying whatever selection is associated with the input.
 * Objects should be closed on input and will be closed on output.
 * In order for the disk file structures to be copied, the output file
 * must be sufficiently defined that it can be written; the copy does
 * not apply any selection/calibration/translation.
 * The copy will be attempted but no errors will be logged until
 * both input and output have been successfully opened.
 * If the contents of the uv data are copied, all associated tables are 
 * copied first.
 * ObitInfoList and ObitThread members are only copied if the output object
 * didn't previously exist.
 * Parent class members are included but any derived class info is ignored.
 * The file etc. info should have been stored in the ObitInfoList:
 * \li "doCalSelect" OBIT_boolean scalar if TRUE, calibrate/select/edit input data.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitUV* ObitUVCopy (ObitUV *in, ObitUV *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOCode iretCode, oretCode;
  gboolean oldExist, doCalSelect;
  ObitInfoType type;
  olong NPIO;
  gint32 dim[MAXINFOELEMDIM];
  ObitHistory *inHist=NULL, *outHist=NULL;
  olong count;
  ObitIOAccess access;
  olong i;
  ofloat bl, maxbl, maxw, *u, *v, *w;
  gchar *outName=NULL, *today=NULL;
  gchar *routine = "ObitUVCopy";
 
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));
  Obit_retval_if_fail((in!=out), err, out,
 		      "%s: input and output are the same", routine);
  

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitUV(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes other additions only if out newly created */
  if (!oldExist) {
    /* copy */
    out->myDesc = (gpointer)ObitUVDescCopy(in->myDesc, out->myDesc, err);
    /* Don't copy selector */
    if (out->mySel) out->mySel = ObitUnref (out->mySel);
    out->mySel = newObitUVSel (out->name);
    /* Don't copy info */
    /*out->info = ObitInfoListUnref(out->info);*/
    /*out->info = ObitInfoListRef(in->info);*/
    /* Output will initially have no associated tables */
    out->tableList = ObitTableListUnref(out->tableList);
    out->tableList = newObitTableList(out->name);
    /* don't copy ObitUVSel, ObitThread */
  }

  /* If the output object was created this call it cannot be fully
     defined so we're done */
  if (!oldExist) return out;

  doCalSelect = FALSE;
  ObitInfoListGetTest(in->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* if output has file designated, copy data */
  /* test open to fully instantiate input and see if it's OK */
  in->bufferSize = MAX (0, in->bufferSize); /* need buffer */
  iretCode = ObitUVOpen (in, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,in->name, out);

  /* copy Descriptor - this time with full information */
  out->myDesc = ObitUVDescCopy(in->myDesc, out->myDesc, err);
  /* Creation date today */
  today = ObitToday();
  strncpy (out->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
 
  /* Copy number of records per IO to output */
  NPIO = 1000;
  type = OBIT_long;
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListGetTest   (in->info,  "nVisPIO", &type,      dim,  &NPIO);
  ObitInfoListAlwaysPut (out->info, "nVisPIO",  type, dim,  &NPIO);

  /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (out->buffer) ObitIOFreeBuffer(out->buffer); /* free existing */
  out->buffer = NULL;
  out->bufferSize = -1;

  /* test open output */
  oretCode = ObitUVOpen (out, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (out, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    out->buffer = NULL;
    out->bufferSize = 0;
    return out;
  }

  /* Copy tables before data */
  iretCode = CopyTablesSelect (in, out, err);
  if (err->error) {/* add traceback,return */
    out->buffer = NULL;
    out->bufferSize = 0;
    Obit_traceback_val (err, routine,in->name, out);
  }

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (in->myIO, in->info, err);
  oretCode = ObitIOSet (out->myIO, in->info, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /* No data in output yet */
  ((ObitUVDesc*)out->myIO->myDesc)->nvis = 0;
  out->myDesc->nvis = 0;

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (in, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,in->name, out);

  iretCode = ObitUVOpen (in, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,in->name, out);

  /* we're in business, copy */
  count = 0;
  maxbl = -1.0;
  maxw  = -1.0;
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (in, in->buffer, err);
    else iretCode = ObitUVRead (in, in->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;

    /* Baseline statistics */
    u   = in->buffer+in->myDesc->ilocu;
    v   = in->buffer+in->myDesc->ilocv;
    w   = in->buffer+in->myDesc->ilocw;
    for (i=0; i<in->myDesc->numVisBuff; i++) { /* loop over buffer */
      /* Get statistics */
      bl = ((*u)*(*u) + (*v)*(*v));
      maxbl = MAX (maxbl, bl);
      maxw = MAX (fabs(*w), maxw);
      
      /* update data pointers */
      u += in->myDesc->lrec;
      v += in->myDesc->lrec;
      w += in->myDesc->lrec;
    } /* end loop over buffer */

   /* How many */
    out->myDesc->numVisBuff = in->myDesc->numVisBuff;
    count += out->myDesc->numVisBuff;
    oretCode = ObitUVWrite (out, in->buffer, err);
  }
  
  /* Save baseline statistics in the descriptor */
  in->myDesc->maxBL  = sqrt(fabs(maxbl));
  in->myDesc->maxW   = maxw;
  out->myDesc->maxBL = sqrt(fabs(maxbl));
  out->myDesc->maxW  = maxw;

  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,in->name, out);
  
  /* Copy any history  unless Scratch */
  if (!in->isScratch && !out->isScratch) {
    inHist  = newObitDataHistory((ObitData*)in, OBIT_IO_ReadOnly, err);
    outHist = newObitDataHistory((ObitData*)out, OBIT_IO_WriteOnly, err);
    outHist = ObitHistoryCopy (inHist, outHist, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
    inHist  = ObitHistoryUnref(inHist);
    outHist = ObitHistoryUnref(outHist);
  }

  /* unset input buffer (may be multiply deallocated ;'{ ) */
  out->buffer = NULL;
  out->bufferSize = 0;

  /* close files */
  oretCode = ObitUVClose (out, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,out->name, out);

  /* close input */
  iretCode = ObitUVClose (in, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,in->name, out);
  
  /* Make sure something copied - if there is anything to copy */
  Obit_retval_if_fail(((count>0) || (in->myDesc->nvis<=0)), err, out,
 		      "%s: NO Data copied for %s", routine, in->name);
  
  return out;
} /* end ObitUVCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create a UV object similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 */
void ObitUVClone  (ObitUV *in, ObitUV *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOAccess access;
  ObitIOCode iretCode, oretCode;
  ObitHistory *inHist=NULL, *outHist=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong NPIO;
  gboolean doClose;
  gchar *today=NULL;
  gchar *routine = "ObitUVClone";
 
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
  /* Don't copy selector */
  if (out->mySel) out->mySel = ObitUnref (out->mySel);
  out->mySel = newObitUVSel (out->name);
  /* Output will initially have no associated tables */
  out->tableList = ObitTableListUnref(out->tableList);
  out->tableList = newObitTableList(out->name);
  /* don't copy ObitUVSel, ObitThread */

  /* Open if needed to to fully instantiate input and see if it's OK */
  if ((in->myStatus == OBIT_Defined) || (in->myStatus == OBIT_Inactive)) {
    doClose = TRUE;  /* Need to close later */
    iretCode = ObitUVOpen (in, OBIT_IO_ReadWrite, err);
    if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_msg (err, routine,in->name);
  } else doClose = FALSE;

  /* copy Descriptor - this time with full information */
  out->myDesc = ObitUVDescCopy(in->myDesc, out->myDesc, err);
  out->myDesc->nvis = 0;  /* no data yet */

  /* Creation date today */
  today = ObitToday();
  strncpy (out->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
 
  /* Copy number of records per IO to output */
  NPIO = 1000;
  type = OBIT_long;
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListGetTest   (in->info,  "nVisPIO", &type,      dim,  &NPIO);
  ObitInfoListAlwaysPut (out->info, "nVisPIO",  type, dim,  &NPIO);

  /* Open output */
  ObitErrLog(err); /* Show any pending messages as they may get lost */
  access = OBIT_IO_WriteOnly;
  oretCode = ObitUVOpen (out, access, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClearErr(err);
    access = OBIT_IO_ReadWrite;
    oretCode = ObitUVOpen (out, access, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_msg (err, routine, out->name);
  }
  /* Ignore any non error messages */
  if (!err->error) ObitErrClear(err); 

  /* Copy any history unless Scratch */
  if (!in->isScratch && !out->isScratch) {
    inHist  = newObitDataHistory((ObitData*)in, OBIT_IO_ReadOnly, err);
    outHist = newObitDataHistory((ObitData*)out, OBIT_IO_WriteOnly, err);
    outHist = ObitHistoryCopy (inHist, outHist, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    inHist  = ObitHistoryUnref(inHist);
    outHist = ObitHistoryUnref(outHist);
  }

  /* Copy tables */
  iretCode = CopyTablesSelect (in, out, err);
  if (err->error) Obit_traceback_msg (err, routine,in->name);

  /* Close files */
  if (doClose) ObitUVClose (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  ObitUVClose (out, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
} /* end ObitUVClone */

/**
 * Initialize structures and open file.
 * The image descriptor is read if OBIT_IO_ReadOnly, OBIT_IO_ReadCal or 
 * OBIT_IO_ReadWrite and written to disk if opened OBIT_IO_WriteOnly.
 * If access is OBIT_IO_ReadCal then the calibration/selection/editing
 * needed is initialized.
 * See the #ObitUVSel class for a description of the selection and 
 * calibration parameters.
 * After the file has been opened the member, buffer is initialized
 * for reading/storing the data unless member bufferSize is <0.
 * The file etc. info should have been stored in the ObitInfoList:
 * \li "FileType" OBIT_long scalar = OBIT_IO_FITS or OBIT_IO_AIPS 
 *               for file type.
 * \li "nVisPIO" OBIT_long scalar = Maximum number of visibilities
 *               per transfer, this is the target size for Reads (may be 
 *               fewer) and is used to create buffers.
 * \li "Compress" Obit_bool scalar = TRUE indicates output is to be 
 *               in compressed format. (access=OBIT_IO_WriteOnly only).
 * \li "SubScanTime" Obit_float scalar [Optional] if given, this is the 
 *               desired time (day) of a sub scan.  This is used by the 
 *               selector to suggest a value close to this which will
 *               evenly divided the current scan.  See #ObitUVSelSubScan
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite,
 *               OBIT_IO_ReadCal or OBIT_IO_WriteOnly).
 *               If OBIT_IO_WriteOnly any existing data in the output file
 *               will be lost.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitUVOpen (ObitUV *in, ObitIOAccess access, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong need;
  gchar *routine = "ObitUVOpen";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* Same type of access on descriptor */
  in->myDesc->access = access;

  /* If the file is already open - close it  first */
  if ((in->myStatus==OBIT_Active) || (in->myStatus==OBIT_Modified)) {
    if (in->myIO) retCode = ObitUVClose (in, err);
    else retCode = OBIT_IO_OK;
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* set Status */
  in->myStatus = OBIT_Active;

  ObitUVGetSelect (in, in->info, in->mySel, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* create appropriate ObitIO */
  ObitUVSetupIO (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Add reference to tableList */
  in->myIO->tableList = (Obit*)ObitUnref(in->myIO->tableList);
  in->myIO->tableList = (Obit*)ObitRef(in->tableList);

  in->myIO->access = access; /* save access type */

  /*+++++++++++++++++ Actual file open ++++++++++++++++++*/
  /* most of the instructions for the I/O are in the ObitInfoList */
  retCode = ObitIOOpen (in->myIO, access, in->info, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* read or write Headers */
  if ((access == OBIT_IO_ReadOnly) || (access == OBIT_IO_ReadCal) || 
      (access == OBIT_IO_ReadWrite)) {
    /* read header info */
    retCode = ObitIOReadDescriptor(in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) Obit_traceback_val (err, routine, in->name, retCode);

    /* Set descriptors for the output on in to reflect the selection
       by in->mySel,  the descriptors on in->myIO will still describe
       the external representation */
    ObitUVSelSetDesc ((ObitUVDesc*)in->myIO->myDesc,(ObitUVSel*)in->myIO->mySel, in->myDesc, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

    /* If single source selected, get source info in header */
    if (in->mySel->selectSources && (in->mySel->numberSourcesList==1) && 
	(access == OBIT_IO_ReadCal)) {
      ObitUVGetSouInfo (in, err);
    }

   /* Output */
  } else if (access == OBIT_IO_WriteOnly) {
    /* Set descriptors for the output on in to reflect the selection
       by in->mySel,  the descriptors on in->myIO will still describe
       the external representation */
    ObitUVSelGetDesc (in->myDesc, (ObitUVSel*)in->myIO->mySel, (ObitUVDesc*)in->myIO->myDesc, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    
    /* Write header info */
    retCode = ObitIOWriteDescriptor(in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* initialize any Calibration needed - complete output Descriptor, Selector*/
  /* get selection parameters for Read/cal/select */
  if ((access==OBIT_IO_ReadOnly) || (access==OBIT_IO_ReadWrite) || 
      (access==OBIT_IO_ReadCal)) {
    ObitUVSetupCal (in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* Initialize frequency arrays if not done */
  if (in->myDesc->freqArr==NULL) ObitUVGetFreq (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  
  /* Allocate buffer - resize if necessary */
  /* buffer size < 0 => no buffer desired */
  if (in->bufferSize >= 0) {
    need = ObitUVSelBufferSize(in->myDesc, in->mySel);
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

    /* Check if valid memory */
    if (!ObitMemValid (in->buffer)) {
       Obit_log_error(err, OBIT_Error, 
		     "%s: IO buffer not in Valid memory for %s", 
		     routine, in->name);
      return retCode;
    }
  } /* end buffer allocation */

  /* init I/O */
  retCode = ObitIOSet (in->myIO, in->info, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Set I/O to beginning of the file */
  ((ObitUVDesc*)in->myIO->myDesc)->firstVis = 0;
  /* For WriteOnly the file is truncated.*/
  if (access == OBIT_IO_WriteOnly) {
    ((ObitUVDesc*)in->myIO->myDesc)->firstVis = 1;
    ((ObitUVDesc*)in->myIO->myDesc)->nvis = 0;
    in->myDesc->nvis = 0;
  }

  /* save current location */
  in->myDesc->firstVis   = ((ObitUVDesc*)in->myIO->myDesc)->firstVis;
  in->myDesc->numVisBuff = ((ObitUVDesc*)in->myIO->myDesc)->numVisBuff;

  return retCode;
} /* end ObitUVOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitUVClose (ObitUV *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitUVClose";

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
    ObitUVDescCopyDesc(in->myDesc, (ObitUVDesc*)in->myIO->myDesc, err);
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

  /* shutdown any calibration, indexing */
  if (in->myIO->access==OBIT_IO_ReadCal) {
    ObitUVCalShutdown((ObitUVCal*)in->myIO->myCal, err);
    if (err->error)  Obit_traceback_val (err, routine, in->name, retCode);
    ObitUVSelShutdown(in->myIO->mySel, err);
    if (err->error)  Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* Free buffer */
  ObitIOFreeBuffer(in->buffer);
  in->buffer = NULL;

  /* set Status */
  in->myStatus = OBIT_Inactive;

  return retCode;
} /* end ObitUVClose */

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
void ObitUVFullInstantiate (ObitUV *in, gboolean exist, ObitErr *err)
{
  ObitIOAccess access;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean doCalSelect;
  gchar *routine = "ObitUVFullInstantiate";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  if (in->myIO) return;  /* is this needed? */

  /* Open readonly if it should exist, else writeonly */
  if (exist) {
    /* If doCalSelect set use ReadCal */
    doCalSelect = FALSE;
    ObitInfoListGetTest(in->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
    if (doCalSelect) access = OBIT_IO_ReadCal;
    else access = OBIT_IO_ReadOnly;
  }  else access = OBIT_IO_WriteOnly;

  in->bufferSize = -1;  /* Don't need to assign buffer here */
  
  /* Open and close */
  ObitUVOpen(in, access, err);
  ObitUVClose(in, err);
  if (err->error)Obit_traceback_msg (err, routine, in->name);
  in->bufferSize = 0;  /* May need buffer later */
} /* end ObitUVFullInstantiate */

/**
 * Read uv data data from disk.
 * The ObitUVDesc maintains the current location in the file.
 * The number read will be mySel->nVisPIO (until the end of the selected
 * range of visibilities in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities is myDesc->numVisBuff.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 *             if NULL, use the buffer member of in.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitUVRead (ObitUV *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  ofloat *myBuf = data;
  olong need;
  gchar *routine = "ObitUVRead";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

 /* check and see if its open - if not attempt */
  if ((in->myStatus!=OBIT_Active) && (in->myStatus!=OBIT_Modified)) {
    access = OBIT_IO_ReadOnly;
    retCode = ObitIOOpen (in->myIO, access, in->info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback, return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* select internal or external buffer */
  if (myBuf==NULL) {
    myBuf = in->buffer;
    /* Check that internal buffer ( defined in gfloats) large enough */
    need = in->mySel->nVisPIO*in->myDesc->lrec;
    if (need > in->bufferSize) {
      Obit_log_error(err, OBIT_Error, 
		     "IO buffer ( %d) too small, need %d for %s", 
		     in->bufferSize, need, in->name);
      return retCode;
    }
  } 

  /* Check buffer */
  Obit_retval_if_fail((myBuf != NULL), err, retCode,
 		      "%s: No buffer allocated for %s", routine, in->name);

  retCode = ObitIORead (in->myIO, myBuf, err);
  if ((retCode > OBIT_IO_EOF) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* save current location */
  in->myDesc->firstVis   = ((ObitUVDesc*)in->myIO->myDesc)->firstVis;
  in->myDesc->numVisBuff = ((ObitUVDesc*)in->myIO->myDesc)->numVisBuff;

  return retCode;
} /* end ObitUVRead */

/**
 * Read uv data data from disk to multiple buffers.
 * If amp/phase calibration being applied, it is done independently
 * for each buffer using the myCal in the associated in,
 * otherwise, the first buffer is processed and copied to the others.
 * All buffers must be the same size and the underlying dataset the same.
 * The ObitUVDesc maintains the current location in the file.
 * The number read will be mySel->nVisPIO (until the end of the selected
 * range of visibilities in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities is myDesc->numVisBuff.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data)
 * \param nBuff Number of buffers to be filled
 * \param in    Array of pointers to to object to be read; 
 *              must all be to same underlying data set.
 * \param data  Array of pointers to buffers to write results.
 *              If NULL, use buffers on in elements (must be open)
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitUVReadMulti (olong nBuff, ObitUV **in, ofloat **data, 
			    ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong need, ib;
  gchar *routine = "ObitUVReadMulti";

  /* error checks */
  if (err->error) return retCode;
  if (nBuff<=0)   return retCode;

  /* Setup for multiple buffers on first in */
  if (in[0]->nParallel!=nBuff) {
    /* Out with any old */
    if (in[0]->multiBuf) g_free(in[0]->multiBuf);
    if (in[0]->multiBufIO) {
      for (ib=0; ib<in[0]->nParallel; ib++) 
	in[0]->multiBufIO[ib] = ObitIOUnref(in[0]->multiBufIO[ib]);
    }
    /* In with the new */
    in[0]->nParallel  = nBuff;
    in[0]->multiBuf   = g_malloc0(nBuff*sizeof(ofloat*));
    in[0]->multiBufIO = g_malloc0(nBuff*sizeof(ObitIO*));
    /* Fill multi buff values */
    for (ib=0; ib<nBuff; ib++) {  
      in[0]->multiBuf[ib]   = in[ib]->buffer;
      in[0]->multiBufIO[ib] = ObitIORef(in[ib]->myIO);
    }
  } /* end allocate for multiple buffers */

  /* Loop over buffers */
  for (ib=0; ib<nBuff; ib++) {  
    /* select internal or external buffer */
    if ((data==NULL) || (data[ib]==NULL)) {
      in[0]->multiBuf[ib] = in[ib]->buffer;
      /* Check that internal buffer ( defined in gfloats) large enough */
      need = in[ib]->mySel->nVisPIO*in[0]->myDesc->lrec;
      if (need > in[ib]->bufferSize) {
	Obit_log_error(err, OBIT_Error, 
		       "IO buffer ( %d) too small, need %d for %s buff %d", 
		       in[ib]->bufferSize, need, in[ib]->name, ib);
	return retCode;
      }
    } else {  /* Use passed buffer */
      in[0]->multiBuf[ib] = data[ib];
    }
    
    /* Check buffer */
    Obit_retval_if_fail((in[0]->multiBuf[ib] != NULL), err, retCode,
			"%s: No buffer allocated for %s", routine, in[ib]->name);
  } /* End loop over buffers */

  /* Do IO */
  retCode = ObitIOReadMulti (nBuff, in[0]->multiBufIO, in[0]->multiBuf, err);
  if ((retCode > OBIT_IO_EOF) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in[0]->name, retCode);

  /* save current location */
  for (ib=0; ib<nBuff; ib++) {
    in[ib]->myDesc->firstVis   = ((ObitUVDesc*)in[0]->myIO->myDesc)->firstVis;
    in[ib]->myDesc->numVisBuff = ((ObitUVDesc*)in[0]->myIO->myDesc)->numVisBuff;
    ((ObitUVDesc*)in[ib]->myIO->myDesc)->firstVis = 
      ((ObitUVDesc*)in[0]->myIO->myDesc)->firstVis;
    ((ObitUVDesc*)in[ib]->myIO->myDesc)->numVisBuff = 
      ((ObitUVDesc*)in[0]->myIO->myDesc)->numVisBuff;
  }

  return retCode;
} /* end ObitUVReadMulti */

/**
 * Reread uv data data from disk to multiple buffers.
 * Retreives data read in a previous call to ObitUVReadMulti.
 * in which in[0] should be the same as in the call to ObitUVReadMulti
 * If amp/phase calibration being applied, it is done independently
 * for each buffer using the myCal in the associated in,
 * otherwise, data from the first buffer copied to the others.
 * NOTE: this depends on retreiving the data from the first element in 
 * All buffers must be the same size and the underlying dataset the same.
 * The ObitUVDesc maintains the current location in the file.
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities is myDesc->numVisBuff.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data)
 * \param nBuff Number of buffers to be filled
 * \param in    Array of pointers to to object to be read; 
 *              must all be to same underlying data set.
 * \param data  Array of pointers to buffers to write results.
 *              If NULL, use buffers on in elements (must be open)
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitUVReReadMulti (olong nBuff, ObitUV **in, ofloat **data, 
			      ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong need, ib;
  gchar *routine = "ObitUVReReadMulti";

  /* error checks */
  if (err->error) return retCode;
  if (nBuff<=0)   return retCode;

  /* Check buffer assignment */
  Obit_retval_if_fail((in[0]->nParallel>=nBuff), err, retCode,
 		      "%s: Buffer assignment error %d %d", 
		      routine, in[0]->nParallel, nBuff);

  /* Loop over buffers */
  for (ib=0; ib<nBuff; ib++) {  

    /* select internal or external buffer */
    if ((data==NULL) || (data[ib]==NULL)) {
      in[0]->multiBuf[ib] = in[ib]->buffer;
      /* Check that internal buffer ( defined in gfloats) large enough */
      need = in[ib]->mySel->nVisPIO*in[0]->myDesc->lrec;
      if (need > in[ib]->bufferSize) {
	Obit_log_error(err, OBIT_Error, 
		       "IO buffer ( %d) too small, need %d for %s buff %d", 
		       in[ib]->bufferSize, need, in[ib]->name, ib);
	return retCode;
      }
    } else {  /* Use passed buffer */
      in[0]->multiBuf[ib] = data[ib];
    }
    
    /* Check buffer */
    Obit_retval_if_fail((in[0]->multiBuf[ib] != NULL), err, retCode,
			"%s: No buffer allocated for %s", routine, in[ib]->name);
  } /* End loop over buffers */

  /* Retrieve old data */
  retCode = ObitIOReReadMulti (nBuff, in[0]->multiBufIO, in[0]->multiBuf, err);
  if ((retCode > OBIT_IO_EOF) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in[0]->name, retCode);

  /* save current location */
  for (ib=0; ib<nBuff; ib++) {
    in[ib]->myDesc->firstVis   = ((ObitUVDesc*)in[0]->myIO->myDesc)->firstVis;
    in[ib]->myDesc->numVisBuff = ((ObitUVDesc*)in[0]->myIO->myDesc)->numVisBuff;
    ((ObitUVDesc*)in[ib]->myIO->myDesc)->firstVis = 
      ((ObitUVDesc*)in[0]->myIO->myDesc)->firstVis;
    ((ObitUVDesc*)in[ib]->myIO->myDesc)->numVisBuff = 
      ((ObitUVDesc*)in[0]->myIO->myDesc)->numVisBuff;
  }

  return retCode;
} /* end ObitUVReReadMulti */

/**
 * Read data from disk applying selection.
 * The number read will be mySel->nVisPIO (until the end of the selected
 * range of visibilities in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities is myDesc->numVisBuff.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 *             if NULL, use the buffer member of in.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitUVReadSelect (ObitUV *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  ofloat *myBuf = data;
  olong need;
  gchar *routine = "ObitUVReadSelect";

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
    need = in->mySel->nVisPIO*in->myDesc->lrec;
    if (need > in->bufferSize) {
      Obit_log_error(err, OBIT_Error, 
		     "IO buffer ( %d) too small, need %d for %s", 
		     in->bufferSize, need, in->name);
      return retCode;
    }
  } 
  /* Check buffer */
  Obit_retval_if_fail((myBuf != NULL), err, retCode,
 		      "%s: No buffer allocated for %s", routine, in->name);


  retCode = ObitIOReadSelect (in->myIO, myBuf, err);
  if ((retCode > OBIT_IO_EOF) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* save current location */
  in->myDesc->firstVis   = ((ObitUVDesc*)in->myIO->myDesc)->firstVis;
  in->myDesc->numVisBuff = ((ObitUVDesc*)in->myIO->myDesc)->numVisBuff;

  return retCode;
} /* end ObitUVReadSelect */

/**
 * Read data from disk applying selection to multiple buffers.
 * The number read will be mySel->nVisPIO (until the end of the selected
 * range of visibilities in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities is myDesc->numVisBuff.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data)
 * \param nBuff Number of buffers to be filled
 * \param in    Array of pointers to to object to be read; 
 *              must all be to same underlying data set but with 
 *              possible independent calibration
 * \param data  array of pointers to buffers to write results.
 *              If NULL, use buffers on in elements (must be open)
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitUVReadMultiSelect (olong nBuff, ObitUV **in, ofloat **data, 
				  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong need, ib;
  gchar *routine = "ObitUVReadMultiSelect";

  /* error checks */
  if (err->error) return retCode;
  if (nBuff<=0)   return retCode;

  /* Setup for multiple buffers on first in */
  if (in[0]->nParallel!=nBuff) {
    /* Out with any old */
    if (in[0]->multiBuf) g_free(in[0]->multiBuf);
    if (in[0]->multiBufIO) {
      for (ib=0; ib<in[0]->nParallel; ib++) 
	in[0]->multiBufIO[ib] = ObitIOUnref(in[0]->multiBufIO[ib]);
    }
    /* In with the new */
    in[0]->nParallel  = nBuff;
    in[0]->multiBuf   = g_malloc0(nBuff*sizeof(ofloat*));
    in[0]->multiBufIO = g_malloc0(nBuff*sizeof(ObitIO*));
      
    /* Fill multi buff values */
    for (ib=0; ib<nBuff; ib++) {  
      in[0]->multiBuf[ib]   = in[ib]->buffer;
      in[0]->multiBufIO[ib] = ObitIORef(in[ib]->myIO);
    }
  } /* end allocate for multiple buffers */


  /* Loop over buffers */
  for (ib=0; ib<nBuff; ib++) {  

    /* select internal or external buffer */
    if ((data==NULL) || (data[ib]==NULL)) {
      in[0]->multiBuf[ib] = in[ib]->buffer;
      /* Check that internal buffer ( defined in gfloats) large enough */
      need = in[0]->mySel->nVisPIO*in[0]->myDesc->lrec;
      if (need > in[ib]->bufferSize) {
	Obit_log_error(err, OBIT_Error, 
		       "IO buffer ( %d) too small, need %d for %s", 
		       in[0]->bufferSize, need, in[ib]->name);
	return retCode;
      }
    } else {  /* Use passed buffer */
      in[0]->multiBuf[ib] = data[ib];
    }
    
    /* Check buffer */
    Obit_retval_if_fail((in[0]->multiBuf[ib] != NULL), err, retCode,
			"%s: No buffer allocated for %s %d", routine, in[ib]->name, ib);
  } /* End loop over buffers */
    
    retCode = ObitIOReadMultiSelect (nBuff, in[0]->multiBufIO, in[0]->multiBuf, err);
    if ((retCode > OBIT_IO_EOF) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in[0]->name, retCode);

  /* save current location */
  for (ib=0; ib<nBuff; ib++) {
    in[ib]->myDesc->firstVis   = ((ObitUVDesc*)in[0]->myIO->myDesc)->firstVis;
    in[ib]->myDesc->numVisBuff = ((ObitUVDesc*)in[0]->myIO->myDesc)->numVisBuff;
    ((ObitUVDesc*)in[ib]->myIO->myDesc)->firstVis = 
      ((ObitUVDesc*)in[0]->myIO->myDesc)->firstVis;
    ((ObitUVDesc*)in[ib]->myIO->myDesc)->numVisBuff = 
      ((ObitUVDesc*)in[0]->myIO->myDesc)->numVisBuff;
  }

  return retCode;
} /* end ObitUVReadMultiSelect */

/**
 * Reread data from disk applying selection to multiple buffers.
 * Retreives data read in a previous call to ObitUVReadMultiSelect
 * possibly applying new calibration.
 * NOTE: this depends on retreiving the data from the first element in 
 * in which should be the same as in the call to ObitUVReadMultiSelect
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities is myDesc->numVisBuff.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data)
 * \param nBuff Number of buffers to be filled
 * \param in    Array of pointers to to object to be read; 
 *              must all be to same underlying data set but with 
 *              possible independent calibration
 * \param data  array of pointers to buffers to write results.
 *              If NULL, use buffers on in elements (must be open)
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitUVReReadMultiSelect (olong nBuff, ObitUV **in, ofloat **data, 
				    ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong need, ib;
  gchar *routine = "ObitUVReReadMultiSelect";

  /* error checks */
  if (err->error) return retCode;
  if (nBuff<=0)   return retCode;

  /* Check buffer assignment */
  Obit_retval_if_fail((in[0]->nParallel>=nBuff), err, retCode,
 		      "%s: Buffer assignment error %d %d", 
		      routine, in[0]->nParallel, nBuff);

  /* Loop over buffers */
  for (ib=0; ib<nBuff; ib++) {  

    /* select internal or external buffer */
    if ((data==NULL) || (data[ib]==NULL)) {
      in[0]->multiBuf[ib] = in[ib]->buffer;
      /* Check that internal buffer ( defined in gfloats) large enough */
      need = in[0]->mySel->nVisPIO*in[0]->myDesc->lrec;
      if (need > in[ib]->bufferSize) {
	Obit_log_error(err, OBIT_Error, 
		       "IO buffer ( %d) too small, need %d for %s", 
		       in[0]->bufferSize, need, in[ib]->name);
	return retCode;
      }
    } else {  /* Use passed buffer */
      in[0]->multiBuf[ib] = data[ib];
    }
    
    /* Check buffer */
    Obit_retval_if_fail((in[0]->multiBuf[ib] != NULL), err, retCode,
			"%s: No buffer allocated for %s", routine, in[ib]->name);
  } /* End loop over buffers */
    
  /* Fetch data from buffer */
    retCode = ObitIOReReadMultiSelect (nBuff, in[0]->multiBufIO, in[0]->multiBuf, err);
    if ((retCode > OBIT_IO_EOF) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in[0]->name, retCode);

  /* save current location */
  for (ib=0; ib<nBuff; ib++) {
    in[ib]->myDesc->firstVis   = ((ObitUVDesc*)in[0]->myIO->myDesc)->firstVis;
    in[ib]->myDesc->numVisBuff = ((ObitUVDesc*)in[0]->myIO->myDesc)->numVisBuff;
    ((ObitUVDesc*)in[ib]->myIO->myDesc)->firstVis = 
      ((ObitUVDesc*)in[0]->myIO->myDesc)->firstVis;
    ((ObitUVDesc*)in[ib]->myIO->myDesc)->numVisBuff = 
      ((ObitUVDesc*)in[0]->myIO->myDesc)->numVisBuff;
  }

  return retCode;
} /* end ObitUVReReadMultiSelect */

/**
 * Write information to disk.
 * The data in the buffer will be written starting at visibility
 * myDesc->firstVis and the number written will be myDesc->numVisBuff
 * which should not exceed mySel->nVisPIO if the internal buffer is used.
 * myDesc->firstVis will be maintained and need not be changed for
 * sequential writing.
 * NB: If the same UV data is being both read and rewritten, 
 * use ObitUVRewrite instead of ObitUVWrite.
 * \param in Pointer to object to be written.
 * \param data pointer to buffer containing input data.
 *             if NULL, use the buffer member of in.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitUVWrite (ObitUV *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  ofloat *myBuf = data;
  olong need;
  gchar *routine = "ObitUVWrite";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

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
    need = in->mySel->nVisPIO*in->myDesc->lrec;
    if (need > in->bufferSize) {
      Obit_log_error(err, OBIT_Error, 
		     "IO buffer ( %d) too small, need %d for %s", 
		     in->bufferSize, need, in->name);
      return retCode;
    }
  } 
  /* Check buffer */
  Obit_retval_if_fail((myBuf != NULL), err, retCode,
 		      "%s: No buffer allocated for %s", routine, in->name);


  /* set number and location to write on myIO descriptor */
  in->myDesc->firstVis = MAX (1, in->myDesc->firstVis);
  ((ObitUVDesc*)in->myIO->myDesc)->firstVis   = in->myDesc->firstVis;
  ((ObitUVDesc*)in->myIO->myDesc)->numVisBuff = in->myDesc->numVisBuff;

  /* most of the instructions for the I/O are in the ObitInfoList */
  retCode = ObitIOWrite (in->myIO, myBuf, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* set Status */
  in->myStatus = OBIT_Modified;

  /* save current location */
  in->myDesc->firstVis   = ((ObitUVDesc*)in->myIO->myDesc)->firstVis;
  in->myDesc->numVisBuff = ((ObitUVDesc*)in->myIO->myDesc)->numVisBuff;
  in->myDesc->nvis       = ((ObitUVDesc*)in->myIO->myDesc)->nvis;

  return retCode;
} /* end ObitUVWrite */

/**
 * Rewrite information to disk.
 * This routines assumes that the input UV is also being read and so
 * the firstVis values in the UV descriptors on both in and its 
 * IO member are left unchanged.  Otherwise a call to this routine 
 * is equivalent to ObitUVWrite
 * The data in the buffer will be written starting at visibility
 * myDesc->firstVis and the number written will be myDesc->numVisBuff
 * which should not exceed mySel->nVisPIO if the internal buffer is used.
 * \param in Pointer to object to be written.
 * \param data pointer to buffer containing input data.
 *             if NULL, use the buffer member of in.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitUVRewrite (ObitUV *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong firstVis;
  gchar *routine = "ObitUVRewrite";

  /* error checks */
  if (err->error) return retCode;

  /* Save firstVis on input */
  firstVis = in->myDesc->firstVis;

  /* Actual write */
  retCode = ObitUVWrite(in, data, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Restore firstVis */ 
  in->myDesc->firstVis = firstVis;
  ((ObitUVDesc*)(in->myIO->myDesc))->firstVis = firstVis;

  return retCode;
} /* end ObitUVRewrite */

/**
 * Return a ObitTable Object to a specified table associated with
 * the input ObitUV.  
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
newObitUVTable (ObitUV *in, ObitIOAccess access, 
		gchar *tabType, olong *tabVer, ObitErr *err)
{
  /* Call ObitData function */
  return newObitDataTable ((ObitData*)in, access, tabType, tabVer, err);
} /* end newObitUVTable */

/**
 * Destroy a specified table(s) associated with the input ObitUV.  
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
ObitIOCode ObitUVZapTable (ObitUV *in, gchar *tabType, olong tabVer, 
			   ObitErr *err)
{
  /* Call ObitData function */
  return ObitDataZapTable ((ObitData*)in, tabType, tabVer, err);
} /* end ObitUVZapTable */

/**
 * Copies the associated tables from one ObitUV to another.
 * \param in      The ObitUV with tables to copy.
 * \param out     An ObitUV to copy the tables to, old ones replaced.
 * \param exclude a NULL termimated list of table types NOT to copy.
 *                If NULL, use include
 * \param include a NULL termimated list of table types to copy.
 *                ignored if exclude nonNULL.
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitUVCopyTables (ObitUV *in, ObitUV *out, gchar **exclude,
			     gchar **include, ObitErr *err)
{
  /* Call ObitData function */
  return ObitDataCopyTables ((ObitData*)in, (ObitData*)out, 
			     exclude, include, err);
} /* end ObitUVCopyTables */

/**
 * Update any disk resident structures about the current tables.
 * \param in   Pointer to object to be updated.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitUVUpdateTables (ObitUV *in, ObitErr *err)
{
  /* Call ObitData function */
  return ObitDataUpdateTables ((ObitData*)in, err);
} /* end ObitUVUpdateTables */

/**
 * Fills in frequency information in Descriptor from header and FQ table.
 * These are the myDesc->freqArr and myDesc->fscale array members.
 * Uses source dependent frequency info if available from in->info
 * \li "SouIFOff" OBIT_double (nif,1,1) Source frequency offset per IF
 * \li "SouBW"    OBIT_double (1,1,1)   Channel Bandwidth
 *
 * \param in      The ObitUV with descriptor to update.
 * \param err     ObitErr for reporting errors.
 */
void ObitUVGetFreq (ObitUV* in, ObitErr *err) 
{
  ObitUVDesc* desc = NULL;
  ObitTable *tab = NULL;
  ObitTableFQ *fqtab = NULL;
  odouble *SouIFOff=NULL, SouBW=0.0;
  olong ver;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "ObitUVGetFreq";
  
  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  desc = in->myDesc;
  /* Source info available? */
  ObitInfoListGetP    (in->info, "SouIFOff",  &type, dim, (gpointer)&SouIFOff);
  /* consistency check - expected number of values? */
  if ((desc->jlocif>=0) && desc->inaxes[desc->jlocif]>dim[0]) SouIFOff = NULL;
  ObitInfoListGetTest (in->info, "SouBW",     &type, dim, &SouBW);
  
  /* get FQ table if one exists */
  ver = 1;
  tab = newObitUVTable (in, OBIT_IO_ReadOnly, "AIPS FQ", &ver, err);
  if (tab!=NULL) {
    
    fqtab = ObitTableFQConvert(tab);
    tab = ObitTableUnref(tab);
  } else fqtab = NULL;
  
  /* update descriptors */
  ObitErrLog(err); /* Show any pending messages as they may get lost */
  ObitUVDescGetFreq (desc, (Obit*)fqtab, SouIFOff, err);
  if (err->error) goto cleanup;

  /* Also need it on the IO descriptor */
  if (in->myIO) {
    desc = in->myIO->myDesc;
    ObitUVDescGetFreq (desc, (Obit*)fqtab, SouIFOff, err);
    if (err->error) goto cleanup;
  }

  /* done with fq table */
 cleanup: fqtab = ObitTableUnref(fqtab);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end  ObitUVGetFreq */

/**
 * Determine subarray information for an ObitUV and add to the descriptors.
 * \param in    Object to obtain data from and with descriptors to update.
 *              Updated on in->myDesc and in->myIO->myDesc:
 * \li numSubA  Number of subarrays, always >0
 * \li numAnt   Array of maximum antenna numbers per subarray
 *              will be allocaded here and must be gfreeed when done
 *              NULL returned if no AN tables
 * \li maxAnt   Maximum antenna number in numAnt.  0 if no AN tables.
 * \param *err  ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitUVGetSubA (ObitUV *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong highANver, iANver;
  ObitTableAN *ANTable = NULL;
  ObitUVDesc  *desc = NULL, *IOdesc = NULL;
  gchar *routine = "ObitUVGetSubA";
 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitUVIsA(in));
  g_assert (ObitIOIsA(in->myIO));
  g_assert (ObitUVDescIsA(in->myIO->myDesc));

  /* pointer to descriptors */
  desc   = (ObitUVDesc*)in->myDesc;
  IOdesc = (ObitUVDesc*)in->myIO->myDesc;

  /* Default outputs */
  desc->numSubA = 1;
  desc->maxAnt  = 0;
  IOdesc->numSubA = 1;
  IOdesc->maxAnt  = 0;

  /* How many AN tables (subarrays) */
  highANver = ObitTableListGetHigh (in->tableList, "AIPS AN");

  /* are there any? */
  if (highANver <= 0) return OBIT_IO_OK;

  /* allocate arrays */
  desc->numSubA = highANver;
  if (desc->numAnt) g_free(desc->numAnt);
  desc->numAnt = g_malloc0(highANver*sizeof(oint*));
  IOdesc->numSubA = highANver;
  if (IOdesc->numAnt) g_free(IOdesc->numAnt);
  IOdesc->numAnt = g_malloc0(highANver*sizeof(oint*));
  ObitErrLog(err); /* Show any pending messages as they may get lost */

  /* Loop over AN tables */
  for (iANver=1; iANver<=highANver; iANver++) {

    /* Get table */
    ANTable = 
      newObitTableANValue (in->name, (ObitData*)in, &iANver, OBIT_IO_ReadOnly, 0, 0, err);
    /* Ignore missing table */
    if (ANTable==NULL) {
      ObitErrClear(err); 
      continue;
    }
    if (err->error) goto cleanup;

    /* Get info from table */
    retCode = ObitTableANGetInfo (ANTable, &(desc->numAnt[iANver-1]), NULL, NULL, err);
    /* Ignore missing table */
    if (retCode==OBIT_IO_OpenErr) {
      ObitErrClear(err); 
      retCode = OBIT_IO_OK;
      continue;
    }
    if ((err->error) || (retCode!=OBIT_IO_OK)) goto cleanup;
 
    /* save to I/O descriptor */
    IOdesc->numAnt[iANver-1] = desc->numAnt[iANver-1];

    /* max antenna number */
    desc->maxAnt = MAX (desc->maxAnt, desc->numAnt[iANver-1]);

    /* release table object */
    ANTable = ObitTableANUnref(ANTable);
  }

  /* save to I/O descriptor */
  IOdesc->maxAnt = desc->maxAnt;

 cleanup:ANTable = ObitTableANUnref(ANTable);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  return retCode;
} /* end ObitUVGetSubA  */

/**
 * Reposition IO to beginning of file
 * \param in   Pointer to object to be rewound.
 * \param err  ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitUVIOSet (ObitUV *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitUVIOSet";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  in->myDesc->firstVis   = 1;
  in->myDesc->numVisBuff = 0;

  /* Barf and die if inactive */
  Obit_retval_if_fail((in->myStatus!=OBIT_Inactive), err, retCode,
 		      "%s: IO inactive for %s", routine, in->name);
  
  return ObitIOSet (in->myIO, in->info, err);
} /* end ObitUVIOSet */

/**
 * Get source position.  
 * If single source file get from uvDesc, 
 * if multisource read from SU table
 * Checks that only one source selected.
 * Also fill in position like information in the descriptor 
 * for multi-source datasets
 * \param  uvdata  Data object from which position sought
 * \param  ra      [out] RA at mean epoch (deg)
 * \param  dec     [out] Dec at mean epoch (deg)
 * \param  err     Error stack, returns if not empty.
 */
void  ObitUVGetRADec (ObitUV *uvdata, odouble *ra, odouble *dec, 
		      ObitErr *err)
{
  ObitSourceList *sList=NULL;
  ObitSource     *source=NULL;
  ObitTable      *tmpTable;
  ObitTableSU    *SUTable=NULL;
  ObitUVSel      *sel;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong i, j, SourID, ver, count;
  olong  Qual;
  gboolean found;
  gchar *sptr, souCode[5], *tabType = "AIPS SU";
  gchar *routine = "ObitUVGetRADec";

  /* error checks */
  g_assert(ObitUVIsA(uvdata));
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* descriptor ilocsu > 0 indicates multisource */
  if (uvdata->myDesc->ilocsu >= 0) { /* multisource */
    sel = uvdata->mySel; /* Selector */

    /* Get SU table */
    ver = 1;
    tmpTable = newObitUVTable (uvdata, OBIT_IO_ReadWrite, tabType, &ver, err);
    SUTable = ObitTableSUConvert(tmpTable);
    tmpTable = ObitTableUnref(tmpTable);
    if (err->error) Obit_traceback_msg (err, routine, uvdata->name);

    /* Make sure selector initialized */
    if (!sel->sources) {

      /* Get selected sources */
      Qual = -1;   /* Qualifier - default = all */
      ObitInfoListGetTest(uvdata->info, "Qual", &type, dim, &Qual);
      /* Cal code */
      souCode[0] =  souCode[1] =  souCode[2] =  souCode[3] =  ' '; souCode[4] = 0; 
      ObitInfoListGetTest(uvdata->info, "souCode", &type, dim, souCode);
      if (ObitInfoListGetP(uvdata->info, "Sources", &type, dim, (gpointer)&sptr)) {
	/* Count actual entries in source list */
	count = 0;  j = 0;
	for (i=0; i<dim[1]; i++) {
	  if ((sptr[j]!=' ') || (sptr[j+1]!=' ')) count++;
	  j += dim[0];
	}
	/* In case all selected */
	ObitTableFullInstantiate ((ObitTable*)SUTable, TRUE, err);  /* Make sure fully defined */
	if(err->error)  Obit_traceback_msg (err, routine, SUTable->name);
	sel->numberSourcesList = SUTable->myDesc->nrow; 
	sel->sources = g_realloc(sel->sources, sel->numberSourcesList*sizeof(olong));
      } else { /* Trouble - no sources selected */
	Obit_log_error(err, OBIT_Error,"%s: Sources not specified on %s",
		       routine, uvdata->name);
	SUTable = ObitTableSUUnref(SUTable);   /* Done with table */
	return;
      }

      /* Do lookup */
      dim[1] = count;
      ObitTableSULookup (SUTable, dim, sptr, Qual, souCode,
			 sel->sources, &sel->selectSources,  
			 &sel->numberSourcesList, err);
      if(err->error)  Obit_traceback_msg (err, routine, SUTable->name);
    }

    /* There must be exactly one source selected */
    if (uvdata->mySel->numberSourcesList!=1) {
      Obit_log_error(err, OBIT_Error,"%s: Exactly one source must be selected in %s",
		   routine, uvdata->name);
      SUTable = ObitTableSUUnref(SUTable);   /* Done with table */
      return;
    }

    /* Which source is it? */
    SourID = uvdata->mySel->sources[0];

    /* Get source list */
    sList = ObitTableSUGetList (SUTable, err);
    SUTable = ObitTableSUUnref(SUTable);   /* Done with table */
    if (err->error) Obit_traceback_msg (err, routine, SUTable->name);

    /* Look up source in table */
    found = FALSE;
    for (i=0; i<sList->number; i++) {
      source = sList->SUlist[i];
      if (source->SourID == SourID) {
	*ra  = source->RAMean;
	*dec = source->DecMean;
	/* Set equinox on descriptor as well */
	uvdata->myDesc->equinox = source->equinox;
	/* Observed RA, dec if not given */
	if (uvdata->myDesc->obsra==0.0)  uvdata->myDesc->obsra  = source->RAMean;
	if (uvdata->myDesc->obsdec==0.0) uvdata->myDesc->obsdec = source->DecMean;
	/* Just in case */
	uvdata->myDesc->crval[uvdata->myDesc->jlocr] = source->RAMean;
	uvdata->myDesc->crval[uvdata->myDesc->jlocd] = source->DecMean;

	found = TRUE;
      }
    }
    sList = ObitSourceListUnref(sList);  /* Done with list */

    /* Check if found */
    if (!found) {
      Obit_log_error(err, OBIT_Error,"%s: Selected source %d not found in list %s",
		   routine, SourID, uvdata->name);
      return;
    }

  } else { /* single source */
    *ra  = uvdata->myDesc->crval[uvdata->myDesc->jlocr];
    *dec = uvdata->myDesc->crval[uvdata->myDesc->jlocd];
  }
} /* end ObitUVGetRADec */

/**
 * Get source information, position, velocity, update uv descriptor
 * If single source file get from uvDesc, 
 * if multisource read from SU table
 * Checks that only one source selected.
 * Also fill in position like information in the descriptor 
 * for multi-source datasets
 * If source table is read, the source dependent IF offsets and bandwidth
 * are written to the in->info object as 
 * \li "SouIFOff" OBIT_double (nif,1,1) Source frequency offset per IF
 * \li "SouBW"    OBIT_double (1,1,1)   Bandwidth
 *
 * \param  uvdata  Data object from which info sought
 * \param  err     Error stack, returns if not empty.
 */
void  ObitUVGetSouInfo (ObitUV *uvdata, ObitErr *err)
{
  ObitSourceList *sList=NULL;
  ObitSource     *source=NULL;
  ObitTable      *tmpTable;
  ObitTableSU    *SUTable=NULL;
  ObitUVSel      *sel;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1, 1, 1, 1, 1,};
  olong i, j, SourID, ver, count, iif;
  olong  Qual;
  gboolean found;
  odouble ra=0.0, dec=0.0, LSRVel=0.0, RestFreq=0.0;
  gchar *sptr, souCode[5], *tabType = "AIPS SU";
  gchar velRef[MAXKEYCHARTABLESU], velDef[MAXKEYCHARTABLESU];
  gchar *routine = "ObitUVGetRADec";

  /* error checks */
  g_assert(ObitUVIsA(uvdata));
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  
  /* Default velocity definitions */
  strncpy (velRef, "LSR",   3);
  strncpy (velDef, "RADIO", 5);

  /* descriptor ilocsu > 0 indicates multisource */
  if (uvdata->myDesc->ilocsu >= 0) { /* multisource */
    sel = uvdata->mySel; /* Selector */

    /* Get SU table */
    ver = 1;
    tmpTable = newObitUVTable (uvdata, OBIT_IO_ReadWrite, tabType, &ver, err);
    SUTable = ObitTableSUConvert(tmpTable);
    tmpTable = ObitTableUnref(tmpTable);
    if (err->error) Obit_traceback_msg (err, routine, uvdata->name);

    /* Make sure selector initialized */
    if (!sel->sources) {

      /* Get selected sources */
      Qual = -1;   /* Qualifier - default = all */
      ObitInfoListGetTest(uvdata->info, "Qual", &type, dim, &Qual);
      /* Cal code */
      souCode[0] =  souCode[1] =  souCode[2] =  souCode[3] =  ' '; souCode[4] = 0; 
      ObitInfoListGetTest(uvdata->info, "souCode", &type, dim, souCode);
      if (ObitInfoListGetP(uvdata->info, "Sources", &type, dim, (gpointer)&sptr)) {
	/* Count actual entries in source list */
	count = 0;  j = 0;
	for (i=0; i<dim[1]; i++) {
	  if ((sptr[j]!=' ') || (sptr[j+1]!=' ')) count++;
	  j += dim[0];
	}
	/* In case all selected */
	ObitTableFullInstantiate ((ObitTable*)SUTable, TRUE, err);  /* Make sure fully defined */
	if(err->error)  Obit_traceback_msg (err, routine, SUTable->name);
	sel->numberSourcesList = SUTable->myDesc->nrow; 
	sel->sources = g_realloc(sel->sources, sel->numberSourcesList*sizeof(olong));
      } else { /* Trouble - no sources selected */
	Obit_log_error(err, OBIT_Error,"%s: Sources not specified on %s",
		       routine, uvdata->name);
	SUTable = ObitTableSUUnref(SUTable);   /* Done with table */
	return;
      }

      /* Do lookup */
      dim[1] = count;
      ObitTableSULookup (SUTable, dim, sptr, Qual, souCode,
			 sel->sources, &sel->selectSources,  
			 &sel->numberSourcesList, err);
      if(err->error)  Obit_traceback_msg (err, routine, SUTable->name);
    }

    /* There must be exactly one source selected */
    if (uvdata->mySel->numberSourcesList!=1) {
      Obit_log_error(err, OBIT_Error,"%s: Exactly one source must be selected in %s",
		   routine, uvdata->name);
      SUTable = ObitTableSUUnref(SUTable);   /* Done with table */
      return;
    }

    /* Which source is it? */
    SourID = uvdata->mySel->sources[0];

    /* Get source list */
    sList = ObitTableSUGetList (SUTable, err);

    /* Actual velocity definitions */
    strncpy (velRef, SUTable->velType, MAXKEYCHARTABLESU-1);
    strncpy (velDef, SUTable->velDef,  MAXKEYCHARTABLESU-1);

    SUTable = ObitTableSUUnref(SUTable);   /* Done with table */
    if (err->error) Obit_traceback_msg (err, routine, SUTable->name);

    /* Look up source in table */
    found = FALSE;
    for (i=0; i<sList->number; i++) {
      source = sList->SUlist[i];
      if (source->SourID == SourID) {
	ra  = source->RAMean;
	dec = source->DecMean;
	/* Set equinox on descriptor as well */
	uvdata->myDesc->equinox = source->equinox;
	/* Observed RA, dec if not given */
	if (uvdata->myDesc->obsra==0.0)  uvdata->myDesc->obsra  = source->RAMean;
	if (uvdata->myDesc->obsdec==0.0) uvdata->myDesc->obsdec = source->DecMean;
	/* Just in case */
	uvdata->myDesc->crval[uvdata->myDesc->jlocr] = source->RAMean;
	uvdata->myDesc->crval[uvdata->myDesc->jlocd] = source->DecMean;

	/* Frequency/velocity info for first IF selected */
	iif = sel->startIF; /* Account for selection in IF */
	RestFreq  = source->RestFreq[iif-1];
	LSRVel    = source->LSRVel[iif-1];
	uvdata->myDesc->restFreq = RestFreq;
	uvdata->myDesc->altRef   = LSRVel;
	uvdata->myDesc->altCrpix = uvdata->myDesc->crpix[uvdata->myDesc->jlocf];
	if (!strncmp(velDef,"OPTICAL", 7))  uvdata->myDesc->VelDef = 0;
	if (!strncmp(velDef,"RADIO", 5))    uvdata->myDesc->VelDef = 1;	
	if (!strncmp(velRef,"LSR", 3))      uvdata->myDesc->VelReference = 1;
	if (!strncmp(velRef,"HELIO", 5))    uvdata->myDesc->VelReference = 2;
	if (!strncmp(velRef,"OBSERVER", 8)) uvdata->myDesc->VelReference = 3;

	iif = sel->startIF; /* Account for selection in IF */
	dim[0] = sel->numberIF;
	ObitInfoListAlwaysPut (uvdata->info, "SouIFOff", OBIT_double, dim, &source->FreqOff[iif-1]);
	dim[0] = 1;
	ObitInfoListAlwaysPut (uvdata->info, "SouBW", OBIT_double, dim, &source->Bandwidth);

	found = TRUE;
      }
    }
    sList = ObitSourceListUnref(sList);  /* Done with list */

    /* Check if found */
    if (!found) {
      Obit_log_error(err, OBIT_Error,"%s: Selected source %d not found in list %s",
		   routine, SourID, uvdata->name);
      return;
    }

  } else { /* single source */
    ra  = uvdata->myDesc->crval[uvdata->myDesc->jlocr];
    dec = uvdata->myDesc->crval[uvdata->myDesc->jlocd];
  }

  /* Update I/O headers */
  ((ObitUVDesc*)uvdata->myIO->myDesc)->equinox = uvdata->myDesc->equinox;
  ((ObitUVDesc*)uvdata->myIO->myDesc)->crval[((ObitUVDesc*)uvdata->myIO->myDesc)->jlocr] = ra;
  ((ObitUVDesc*)uvdata->myIO->myDesc)->crval[((ObitUVDesc*)uvdata->myIO->myDesc)->jlocd] = dec;
  ((ObitUVDesc*)uvdata->myIO->myDesc)->restFreq = RestFreq;
  ((ObitUVDesc*)uvdata->myIO->myDesc)->altRef   = LSRVel;
  ((ObitUVDesc*)uvdata->myIO->myDesc)->altCrpix = uvdata->myDesc->altCrpix;
  ((ObitUVDesc*)uvdata->myIO->myDesc)->VelReference = uvdata->myDesc->VelReference;
  ((ObitUVDesc*)uvdata->myIO->myDesc)->VelDef = uvdata->myDesc->VelDef;
} /* end ObitUVGetSouInfo */

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
void ObitUVWriteKeyword (ObitUV *in, 
			   gchar* name, ObitInfoType type, gint32 *dim, 
			   gconstpointer data, ObitErr *err)
{
  ObitInfoListPut(in->myDesc->info, name, type, dim, data, err);
  ObitInfoListPut(((ObitUVDesc*)in->myIO->myDesc)->info, 
		  name, type, dim, data, err);
  in->myStatus = OBIT_Modified;
} /* end ObitUVWriteKeyword */

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
void ObitUVReadKeyword (ObitUV *in, 
			  gchar* name, ObitInfoType *type, gint32 *dim, 
			  gpointer data, ObitErr *err)
{
  ObitInfoListGet(((ObitUVDesc*)in->myIO->myDesc)->info, 
		  name, type, dim, data, err);
} /* end ObitUVReadKeyword */


/**
 * Create channel selection flag table.
 * If a flag table is currently selected it is copied to a new AIPS FG
 * table on the uv data and channel selection  added.
 * If no flag table is selected then an new table is created.
 * The new flagging entries includ all channels and IFs in the range
 * specified in the UVSel values of BChan, EChan, BIF and EIF that are
 * NOT specified in the IChanSel array
 * \param in       UV with selection to which the new flag table is to be attached.
 * \param dim      Dimensionality of IChanSel
 * \param IChanSel Channel selection in groups of 4 values:
 *                 Sets of channels are selected by groups of 4 parameters:
 *                [0] 1-rel start channel number in block (def . 1)
 *                [1] 1-rel end channel number, (def. to end)
 *                [2] channel increment
 *                [3] 1-rel IF number (def all)
 *                [5,0,2,1] means every other channel from 5 to the highest
 *                in IF 1.  All zeroes means all channels and all IFs 
 *                selected by BChan, EChan, BIF, EIF.
 *                List terminated by group of 4 zeroes.
 * \param err  ObitErr for reporting errors.
 * \return the AIPS FG table version number, 
 *  -1 on no selection or failure.
 */
olong ObitUVChanSel (ObitUV *in, gint32 *dim, olong *IChanSel, 
		    ObitErr *err)
{
  olong BChan=0, EChan=0, BIF=0, EIF=0, nChan, nIF=1;
  olong i, j, k, is, ie, ii, iif;
  olong oldFGver, newFGver=-1, iFGRow;
  ObitTableFG *inFGTable=NULL, *outFGTable=NULL;
  ObitTableFGRow *row=NULL;
  gboolean **ChanIF=NULL;
  gchar reason[25];
  gchar *routine = "ObitUVChanSel";

  /* Anything requested? */
  if ((dim[1]<1) || ((IChanSel[0]==0) && (IChanSel[1]==0) &&
      (IChanSel[2]==0) && (IChanSel[3]==0))) return newFGver;

  /* error checks */
  if (err->error) return newFGver;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Fully instantiate UV data if needed */
  ObitUVFullInstantiate (in, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, newFGver);

  /* tell about it */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Continuum Channel selection via FG table");
  /* Get selection */
  if (in->mySel) {
    if (in->mySel->doFlag) oldFGver = in->mySel->FGversion;
    else oldFGver = -1;
    BChan    = in->mySel->startChann;
    EChan    = BChan + in->mySel->numberChann-1;
    BIF      = in->mySel->startIF;
    EIF      = BIF + in->mySel->numberIF-1;
  } else { /* Default selection */
    oldFGver = -1;
    BChan = 1;
    if (in->myDesc->jlocf>=0) EChan = in->myDesc->inaxes[in->myDesc->jlocf];
    else EChan = 1;
    BIF = 1;
    if (in->myDesc->jlocif>=0) EIF = in->myDesc->inaxes[in->myDesc->jlocif];
    else BIF = 1;
  }

  /* New FG table - make new one */
  newFGver = ObitTableListGetHigh (in->tableList, "AIPS FG") +1;
  /* Get table */
  outFGTable = 
    newObitTableFGValue (in->name, (ObitData*)in, &newFGver, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, newFGver);

  /* Copy old table or write new? */
  if (oldFGver>0) {  /* Copy */
    inFGTable = 
      newObitTableFGValue (in->name, (ObitData*)in, &oldFGver, OBIT_IO_ReadOnly, err);
    outFGTable =  ObitTableFGCopy (inFGTable, outFGTable, err);
    if (err->error)  goto cleanup;
    inFGTable = ObitTableFGUnref(inFGTable);
  }

  /* Add new selection - first make arrays of channels/IFs selected */
  ObitTableFGOpen (outFGTable, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;

  /* Create Row */
  row = newObitTableFGRow (outFGTable);
  
  /* Attach row to output buffer */
  ObitTableFGSetRow (outFGTable, row, err);
  if (err->error) goto cleanup;

  /* Initialize solution row */
  row->SourID  = 0; 
  row->SubA    = 0; 
  row->freqID  = 0; 
  row->ants[0] = 0; 
  row->ants[1] = 0; 
  row->TimeRange[0] = -1.0e20; 
  row->TimeRange[1] =  1.0e20; 
  row->ifs[0]    = BIF; 
  row->ifs[1]    = EIF; 
  row->chans[0]  = BChan; 
  row->chans[1]  = EChan; 
  row->pFlags[0] = 15; 
  row->pFlags[1] = 0; 
  row->pFlags[2] = 0; 
  row->pFlags[3] = 0; 
  g_snprintf (reason, 25, "Data selection");
  row->reason    = reason; 

  /* Create array giving selected channels/IFs */
  nChan = EChan - BChan + 1;
  nIF   = EIF - BIF + 1;
  ChanIF = g_malloc0(nIF*sizeof(gboolean*));
  for (i=0; i<nIF; i++) {
    ChanIF[i] = g_malloc0(nChan*sizeof(gboolean));
    /* Array of chann/IF to flag */
    for (j=0; j<nChan; j++)  ChanIF[i][j] = TRUE;
  }

  /* Loop over IChanSel */
  for (i=0; i<dim[1]; i++) {
    /* After first set of all zeroes - bag it */
    if ((i>0) && (IChanSel[i*4]==0) && (IChanSel[i*4+1]==0) &&
	(IChanSel[i*4+2]==0) && (IChanSel[i*4+3]==0)) break;
    is  = MAX (IChanSel[i*4]-1, 0);
    ie  = IChanSel[i*4+1];
    if (ie<=0) ie = nChan;
    ie = MIN (ie, nChan);
    ii  = IChanSel[i*4+2];
    if (ii<=0) ii = 1;
    iif = IChanSel[i*4+3];
    iif = MIN (iif, nIF);
    for (j=is; j<ie; j+=ii) {
      /* One or all IFs? */
      if (iif>0) {  /* One */
	ChanIF[iif-1][j] = FALSE;  /* Don't flag this one */
      } else {      /* All */
	for (k=0; k<nIF; k++) ChanIF[k][j] = FALSE;  /* Don't flag this one */
      }
    }
  } /* End loop over IChanSel */

  /* Loop through IF, channels writing flags */
  for (i=0; i<nIF; i++) {
    for (j=0; j<nChan; j++) {
      if (ChanIF[i][j]) {
	row->ifs[0]    = i+1;
	row->ifs[1]    = i+1;
	row->chans[0]  = j+1;
	row->chans[1]  = j+1;
	/* Write */
	iFGRow = -1;
	ObitTableFGWriteRow (outFGTable, iFGRow, row, err);
	if (err->error) goto cleanup;
      } /* End deselected IF/channel */
    } /* end loop over Channel */
  } /* end loop over IF */
  
  /* Close output table */
  ObitTableFGClose (outFGTable, err);
  
  /* Cleanup */
 cleanup:
  outFGTable = ObitTableFGUnref(outFGTable);
  row        = ObitTableFGRowUnref(row);
  if (ChanIF!=NULL) {
    for (i=0; i<nIF; i++) if (ChanIF[i]!=NULL) g_free(ChanIF[i]);
    g_free(ChanIF);
  }
  if (err->error) Obit_traceback_val (err, routine, in->name, newFGver);

  return newFGver;
} /* end ObitUVChanSel */


/*-------Private functions called by ObitData class ------*/
/** Private:  Copy Constructor for scratch file*/
static ObitData* newObitDataUVScratch (ObitData *in, ObitErr *err)
{
  return (ObitData*) newObitUVScratch ((ObitUV*)in, err);
} /* end newObitDataUVScratch  */

/** Private: Copy (deep) constructor.  */
static ObitData* ObitDataUVCopy  (ObitData *in, ObitData *out, 
				  ObitErr *err)
{
  return (ObitData*) ObitUVCopy ((ObitUV*)in, (ObitUV*)out, err);
} /* end  ObitDataUVCopy*/

/** Private: Copy structure */
static void ObitDataUVClone (ObitData *in, ObitData *out, ObitErr *err)
{
  ObitUVClone ((ObitUV*)in, (ObitUV*)out, err);
} /* end ObitDataUVClone */

/** Private: Zap */
static ObitData* ObitDataUVZap (ObitData *in, ObitErr *err)
{
  return (ObitData*)ObitUVZap ((ObitUV*)in, err);
} /* end ObitDataUVZap */

/** Private: Rename */
static void ObitDataUVRename (ObitData *in, ObitErr *err)
{
  ObitUVRename ((ObitUV*)in, err);
} /* end ObitDataUVRename */

/** Private: Open */
static ObitIOCode ObitDataUVOpen (ObitData *in, ObitIOAccess access, 
				  ObitErr *err)
{
  return ObitUVOpen ((ObitUV*)in, access, err);
} /* end ObitUDataUVOpen */

/** Private: Close  */
static ObitIOCode ObitDataUVClose (ObitData *in, ObitErr *err)
{
  return ObitUVClose ((ObitUV*)in, err);
} /* end ObitDataUVClose */

/** Private:  Reset IO to start of file  */
static ObitIOCode ObitDataUVIOSet (ObitData *in, ObitErr *err)
{
  return ObitUVIOSet ((ObitUV*)in, err);
} /* end  ObitDataUVIOSet */

/** Private: Assign myIO object */
static void ObitDataUVSetupIO (ObitData *in, ObitErr *err)
{
  ObitUVSetupIO ((ObitUV*)in, err);
} /* end ObitDataUVSetupIO */

/** Private: full instantiation */
static void ObitDataUVFullInstantiate (ObitData *in, gboolean exist, 
				       ObitErr *err)
{
  ObitUVFullInstantiate ((ObitUV*)in, exist, err);
} /* end ObitDataUCFullInstantiate */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();
  myClassInfo.hasScratch    = TRUE; /* Scratch files allowed */

  /* Set function pointers */
  ObitUVClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVClassInfoDefFn (gpointer inClass)
{
  ObitUVClassInfo *theClass = (ObitUVClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVClassInit;
  theClass->newObit       = (newObitFP)newObitUV;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVGetClass;
  theClass->newObitUVScratch  = (newObitUVScratchFP)newObitUVScratch;
  theClass->ObitUVSame    = (ObitUVSameFP)ObitUVSame;
  theClass->ObitCopy      = (ObitCopyFP)ObitUVCopy;
  theClass->ObitClone     = NULL;  /* Different call */
  theClass->ObitClear     = (ObitClearFP)ObitUVClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVInit;
  theClass->ObitUVRead    = (ObitUVReadFP)ObitUVRead;
  theClass->ObitUVReadMulti= (ObitUVReadMultiFP)ObitUVReadMulti;
  theClass->ObitUVReReadMulti= (ObitUVReReadMultiFP)ObitUVReReadMulti;
  theClass->ObitUVReadSelect = (ObitUVReadSelectFP)ObitUVReadSelect;
  theClass->ObitUVReadMultiSelect = 
    (ObitUVReadMultiSelectFP)ObitUVReadMultiSelect;
  theClass->ObitUVReReadMultiSelect = 
    (ObitUVReReadMultiSelectFP)ObitUVReReadMultiSelect;
  theClass->ObitUVWrite   = (ObitUVWriteFP)ObitUVWrite;
  theClass->newObitUVTable= (newObitUVTableFP)newObitUVTable;
  theClass->ObitUVZapTable= (ObitUVZapTableFP)ObitUVZapTable;
  theClass->ObitUVFullInstantiate= 
    (ObitUVFullInstantiateFP)ObitUVFullInstantiate;
  theClass->ObitUVCopyTables= 
    (ObitUVCopyTablesFP)ObitUVCopyTables;
  theClass->ObitUVUpdateTables= 
    (ObitUVUpdateTablesFP)ObitUVUpdateTables;
  /* Function pointers referenced from ObitData class */
  theClass->newObitDataScratch  = (newObitDataScratchFP)newObitDataUVScratch;
  theClass->ObitDataRename  = (ObitDataRenameFP)ObitDataUVRename;
  theClass->ObitDataZap     = (ObitDataZapFP)ObitDataUVZap;
  theClass->ObitDataClone   = (ObitDataCloneFP)ObitDataUVClone;
  theClass->ObitDataCopy    = (ObitDataCopyFP)ObitDataUVCopy;
  theClass->ObitDataOpen    = (ObitDataOpenFP)ObitDataUVOpen;
  theClass->ObitDataClose   = (ObitDataCloseFP)ObitDataUVClose;
  theClass->ObitDataIOSet   = (ObitDataIOSetFP)ObitDataUVIOSet;
  theClass->ObitDataSetupIO = (ObitDataSetupIOFP)ObitDataUVSetupIO;
  theClass->ObitDataFullInstantiate= 
    (ObitDataFullInstantiateFP)ObitDataUVFullInstantiate;
  theClass->ObitDataWriteKeyword= 
    (ObitDataWriteKeywordFP)ObitUVWriteKeyword;
  theClass->ObitDataReadKeyword= 
    (ObitDataReadKeywordFP)ObitUVReadKeyword;

} /* end ObitUVClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUV *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->myIO      = NULL;
  in->myDesc    = newObitUVDesc(in->name);
  in->mySel     = newObitUVSel(in->name);
  in->myStatus  = OBIT_Defined;
  in->buffer    = NULL;
  in->bufferSize= 0;
  in->nParallel = 0;
  in->multiBufIO= NULL;
  in->multiBuf  = NULL;
  in->isScratch = FALSE;

} /* end ObitUVInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUV* cast to an Obit*.
 */
void ObitUVClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUV *in = inn;
  olong ib;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Delete underlying files if isScratch */
  if (in->isScratch) {
    err = newObitErr();     /* for possible messages */
    /* Remove from ObitSystem list */
    ObitSystemFreeScratch ((Obit*)in, err);
    in->isScratch = FALSE;  /* avoid infinite recursion */
    ObitUVZap (in, err);    /* delete files */
    ObitErrLog(err);
    err = ObitErrUnref(err);
  }

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  in->myIO      = ObitUnref(in->myIO);
  in->myDesc    = ObitUVDescUnref(in->myDesc);
  in->mySel     = ObitUVSelUnref(in->mySel);
  in->tableList = ObitUnref(in->tableList);
  if (in->buffer) ObitIOFreeBuffer(in->buffer); 
  if (in->multiBuf) g_free(in->multiBuf);
  if ((in->nParallel>0) && (in->multiBufIO)) {
    for (ib=0; ib<in->nParallel; ib++) 
      in->multiBufIO[ib] = ObitIOUnref(in->multiBufIO[ib]);
  }
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVClear */

/**
 * Get requested information from the ObitInfoList
 * \param in   UVdata.
 * \param info Pointer to InfoList
 * \param sel  pointer to uvdata selector to update.
 * \param err  ObitErr for reporting errors.
 */
static void ObitUVGetSelect (ObitUV *in, ObitInfoList *info, ObitUVSel *sel,
			     ObitErr *err)
{
  ObitInfoType type;
  gint32 i, dim[MAXINFOELEMDIM];
  olong itemp, *iptr, Qual;
  olong iver, j, count=0;
  ofloat ftempArr[10];
  ObitTableSU *SUTable=NULL;
  union ObitInfoListEquiv InfoReal; 
  gchar tempStr[5], souCode[5], *sptr;
  gboolean badChar;
  ObitUVDesc *desc;
  gchar *routine = "ObitUVGetSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA(info));
  g_assert (ObitUVSelIsA(sel));

  /* what type of underlying file? */
  if (!ObitInfoListGet(info, "FileType", &type, dim, 
		       (gpointer)&sel->FileType, err)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		"ObitUVGetSelect: entry FileType not in InfoList Object %s",
		sel->name);
  }

  /* Use descriptor on IO */
  if (in->myIO && in->myIO->myDesc) desc = (ObitUVDesc*)in->myIO->myDesc;
  else  desc = in->myDesc;

  /* Maximum number of visibilities per read/write? [default 100] */
  sel->nVisPIO = 100;
  ObitInfoListGetTest(info, "nVisPIO", &type, dim, &sel->nVisPIO);
  sel->nVisPIO = MAX (1, sel->nVisPIO); /* no fewer than 1 */

  /* Compress output? */
  sel->Compress = FALSE;
  ObitInfoListGetTest(info, "Compress", &type, dim, &sel->Compress);

  /* Following only needed for ReadCal */
  if (in->myDesc->access != OBIT_IO_ReadCal) {
    /* Default selection */
    sel->SubA         = 0;
    sel->FreqID       = 0;
    sel->startChann   = 1;
    sel->numberChann  = desc->inaxes[desc->jlocf];
    sel->startIF      = 1;
    if (desc->jlocif>=0) sel->numberIF = desc->inaxes[desc->jlocif];
    else sel->numberIF = 1;
    sel->numberPoln   = MAX (1, desc->inaxes[desc->jlocs]);
    sel->doPolCal     = FALSE;
    sel->timeRange[0] = -1.0e20; sel->timeRange[1] = 1.0e20;
    return; 
  }

  /* Calibrate/select/edit output? */
  sel->doCalSelect = FALSE;
  ObitInfoListGetTest(info, "doCalSelect", &type, dim, &sel->doCalSelect);

  /* Selection */
  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "Subarray", &type, dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  sel->SubA  = itemp;

  InfoReal.itg = 0;type = OBIT_oint;
  ObitInfoListGetTest(info, "FreqID", &type, dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  sel->FreqID = itemp;

  InfoReal.itg = 0;type = OBIT_oint;
  ObitInfoListGetTest(info, "BChan", &type, dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  if (desc->jlocf>0) itemp = MIN (itemp, desc->inaxes[desc->jlocf]);
  sel->startChann  = MAX (1, itemp);

  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "EChan", &type, dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  if (desc->jlocf>0) itemp = MIN (itemp, desc->inaxes[desc->jlocf]);
  if (itemp>0) sel->numberChann = itemp - sel->startChann+1;
  else  sel->numberChann = desc->inaxes[desc->jlocf] - sel->startChann+1;
  sel->numberChann = MAX (1, MIN (sel->numberChann, desc->inaxes[desc->jlocf]));

  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "BIF", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  if (desc->jlocf>0) itemp = MIN (itemp, desc->inaxes[desc->jlocif]);
  sel->startIF  =  MAX (1, itemp);

  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "EIF", &type, dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  if (desc->jlocf>0) itemp = MIN (itemp, desc->inaxes[desc->jlocif]);
  if (itemp>0) sel->numberIF = itemp - sel->startIF+1;
  else sel->numberIF = desc->inaxes[desc->jlocif] - sel->startIF+1;
  sel->numberIF = MAX (1, MIN (sel->numberIF, desc->inaxes[desc->jlocif]));
  if (desc->jlocif<0) sel->numberIF = 1;

  for (i=0; i<4; i++) tempStr[i] = ' '; tempStr[4] = 0;
  ObitInfoListGetTest(info, "Stokes", &type, dim, &tempStr);
  /* Remove any trailing junk */
  badChar = FALSE;
  for (i=0; i<4; i++) {
    badChar = badChar || (tempStr[i]==0); 
    if (badChar) tempStr[i] = ' ';
  }
  for (i=0; i<4; i++) sel->Stokes[i] = tempStr[i]; sel->Stokes[4] = 0;

  /* Polarization calibration */
  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "doPol", &type, dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  sel->doPolCal = itemp > 0;

  /* amp/phase/delay/rate Calibration */
  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "doCalib", &type, dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  sel->doCal = itemp > 0;
  sel->doCalWt = itemp > 1;
  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "gainUse", &type, dim, &InfoReal);
  sel->calVersion = InfoReal.itg;

  /* Flagging */
  InfoReal.itg = -1; type = OBIT_oint;
  ObitInfoListGetTest(info, "flagVer", &type, dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  sel->doFlag = itemp >= 0;
  sel->FGversion = itemp;

  /* Correlation type */
  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "corrType", &type, dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  sel->corrType = itemp;

  /* Pass all data whether flagged etc */
  sel->passAll = FALSE;
  ObitInfoListGetTest(info, "passAll", &type, dim, &sel->passAll);
  if (sel->passAll) sel->corrType = 1; /* Both correlation types */

  /* Drop Subarray? */
  ObitInfoListGetTest(info, "dropSubA", &type, dim, &sel->dropSubA);

  /* Baseline dependent calibration */
  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "BLVer", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  sel->doBLCal = itemp > 0;
  sel->BPversion = itemp;

  /* Bandpass calibration */
  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(info, "doBand", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  sel->doBPCal = itemp > 0;
  sel->doBand = itemp;
  itemp = 0;
  ObitInfoListGetTest(info, "BPVer", &type, (gint32*)dim, &InfoReal);
  sel->BPversion = itemp;

  /* Spectral smoothing */
  for (i=0; i<3; i++) ftempArr[i] = 0.0; 
  ObitInfoListGetTest(info, "Smooth", &type, (gint32*)dim, &ftempArr);
  for (i=0; i<3; i++) sel->smooth[i] = ftempArr[i];

  /* Data averaging time */
  if (ObitInfoListGetTest(info, "InputAvgTime", &type, (gint32*)dim, 
			  &sel->InputAvgTime)) {
    /* Convert to days */
    sel->InputAvgTime /= 86400.0;
  }

  /* Time Range */
  ftempArr[0] = -1.0e20; ftempArr[1] = 1.0e20; 
  ObitInfoListGetTest(info, "timeRange", &type, (gint32*)dim, &ftempArr);
  for (i=0; i<2; i++) sel->timeRange[i] = ftempArr[i];
  /* If both=0.0 take all */
  if ((sel->timeRange[0]==0.0) && (sel->timeRange[1]==0.0)) {
    sel->timeRange[0] = -1.0e20; sel->timeRange[1] = 1.0e20;
  }

  /* UV Range */
  ftempArr[0] = 0.0; ftempArr[1] = 1.0e20; 
  ObitInfoListGetTest(info, "UVRange", &type, (gint32*)dim, &ftempArr);
  if ((ftempArr[0]<=0.0) && (ftempArr[1]<=0.0)) ftempArr[1] = 1.0e17;
  for (i=0; i<2; i++) sel->UVRange[i] = ftempArr[i]*1.0e3;

  /* Selected antennas */
  if (ObitInfoListGetP(info, "Antennas", &type, (gint32*)dim, (gpointer)&iptr)) {
    sel->numberAntList = dim[0];
    sel->ants = g_realloc(sel->ants, sel->numberAntList*sizeof(olong));
    /* loop copying, checking for deselection */
    sel->selectAnts = FALSE;
    for (i=0; i<sel->numberAntList; i++) {
      sel->ants[i] = abs (iptr[i]);
      sel->selectAnts = sel->selectAnts || (sel->ants[i]>0);
    }
    /* If not selecting by antenna free ants */
    if (!sel->selectAnts) {
      sel->numberAntList = 0;
      g_free(sel->ants); sel->ants = NULL;
    }
  } else {
    sel->selectAnts = FALSE;
    sel->numberAntList = 0;
  }

  /* Selected sources */
  Qual = -1;   /* Qualifier - default = all */
  ObitInfoListGetTest(info, "Qual", &type, dim, &Qual);
  /* Cal code */
  souCode[0] =  souCode[1] =  souCode[2] =  souCode[3] =  ' '; souCode[4] = 0; 
  ObitInfoListGetTest(info, "souCode", &type, dim, souCode);
  if (ObitInfoListGetP(info, "Sources", &type, dim, (gpointer)&sptr)) {
    sel->numberSourcesList = count;
    /* Count actual entries in source list */
    count = 0;  j = 0;
    for (i=0; i<dim[1]; i++) {
      if ((sptr[j]!=' ') || (sptr[j+1]!=' ')) count++;
      j += dim[0];
    }
    sel->numberSourcesList = count;
    if (count>0) {  /* Anything actually specified? */
      /* have to lookup sources - need SU table for this. */
      iver = 1;
      SUTable = newObitTableSUValue (in->name, (ObitData*)in, &iver, OBIT_IO_ReadOnly, 0, err);
      if (SUTable==NULL) {  /* No source table - only one source and it is selected */
	sel->numberSourcesList = 0;
      } else { /* Lookup sources to get numbers */
	/* In case all selected */
	/* Open table */
	ObitTableSUOpen (SUTable, OBIT_IO_ReadOnly, err);
	if (err->error)  Obit_traceback_msg (err, routine,in->name);
	sel->numberSourcesList = MAX (count, SUTable->myDesc->nrow); 
	ObitTableSUClose (SUTable, err);
	if (err->error)  Obit_traceback_msg (err, routine,in->name);
	if (sel->sources)
	  sel->sources = g_realloc(sel->sources, sel->numberSourcesList*sizeof(olong));
	else
	  sel->sources = g_malloc0(sel->numberSourcesList*sizeof(olong));
	
	/* Do lookup */
	dim[1] = count;
	ObitTableSULookup (SUTable, dim, sptr, Qual, souCode, 
			   sel->sources, &sel->selectSources, 
			   &sel->numberSourcesList, err); 
	if (err->error)  Obit_traceback_msg (err, routine, in->name);
	SUTable = ObitTableSUUnref(SUTable); /* release table */
      }
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    } else { /* end of sources specified */
      /* Deallocate sources if needed */
      if (sel->sources) {g_free(sel->sources); sel->sources=NULL;}
    }
  } else { /* no "Sources" specified */ 
   sel->numberSourcesList = 0; /* everything selected */
  }

  /* Subscan interval - default very large */
  sel->SubScanTime = 1.0e20;
  ObitInfoListGetTest(info, "SubScanTime", &type, dim, &sel->SubScanTime);
  sel->SubScanSuggest = sel->SubScanTime;  /* Default */

  /* Data selection parameters */
  sel->numberVis   = 1;
  sel->numberPoln  = MAX (1, MIN (2, desc->inaxes[desc->jlocs]));

} /* end ObitUVGetSelect */

/**
 * Do various operation that have to be done in the ObitUV class
 * rather than the ObitUVCal class.  This mostly involves setups
 * of tables since the ObitUV cannot be visible from the ObitUVCal class.
 * \param in   UV object with UVCal to prepare for calibration
 * \param err  ObitErr for reporting errors.
 */
static void ObitUVSetupCal (ObitUV *in, ObitErr *err)
{
  ObitUVSel *sel = NULL;
  ObitUVCal *cal = NULL;
  olong highVer, iVer, useVer;
  gchar *routine = "ObitUVSetupCal";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(in));
  g_assert (ObitIOIsA(in->myIO));
  g_assert (ObitUVSelIsA(in->mySel));

  if (err->error) Obit_traceback_msg (err, routine, in->name);
  /* Need array information */
  ObitUVGetSubA (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Following only needed for ReadCal */
  if (in->myDesc->access != OBIT_IO_ReadCal) return; 

  /* Create/initialize Calibrator if needed */
  if (in->myIO->myCal==NULL) in->myIO->myCal = (gpointer)newObitUVCal(in->name);
  cal = (ObitUVCal*)in->myIO->myCal;

  /* Setup tables as needed for calibration */
  sel = in->myIO->mySel;

  /* Need all AN tables for Polarization calibration and VLBI
  if (sel->doPolCal || sel->doCal) { ALWAYS NEED ? */
    /* How many AN tables (subarrays) */
    highVer = ObitTableListGetHigh (in->tableList, "AIPS AN");
    cal->numANTable = highVer;

    /* allocate */
    if (cal->ANTables) g_free(cal->ANTables);
    cal->ANTables = g_malloc0(highVer*sizeof(Obit*));

    /* Loop over AN tables */
    for (iVer=0; iVer<highVer; iVer++) {

      /* Get table */
      useVer = iVer+1;
      cal->ANTables[iVer] =
	(Obit*)newObitTableANValue (in->name, (ObitData*)in, &useVer, OBIT_IO_ReadOnly, 0, 0, err);
      if (err->error) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: NO Antenna AN table %d for %s", routine, useVer, in->name);
	return;
      }
    }
 /* }  end AN table setup */
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Get any SU table */
  highVer = ObitTableListGetHigh (in->tableList, "AIPS SU");
  useVer = 1;
  if (highVer>=1)
    cal->SUTable =
      (Obit*) newObitTableSUValue (in->name, (ObitData*)in, &useVer, OBIT_IO_ReadOnly, 0, err);
  else cal->SUTable = ObitTableSUUnref(cal->SUTable);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* BL table for Baseline dependent calibration */
  if (sel->doBLCal) {
    /* if sel->BLversion ==0 use highest */
    highVer = ObitTableListGetHigh (in->tableList, "AIPS BL");
    if (sel->BLversion==0) useVer = highVer;
    else useVer = sel->BLversion;
    cal->BLTable = 
      (Obit*)newObitTableBLValue (in->name, (ObitData*)in, &useVer, OBIT_IO_ReadOnly, 0, 0, err);
    if (cal->BLTable==NULL) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: NO Baseline BL table %d for %s", routine, useVer, in->name);
      return;
    }
  } /* end BL table setup */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* BP table for Bandpass calibration */
  if (sel->doBPCal) {
    /* if sel->BLversion ==0 use highest */
    highVer = ObitTableListGetHigh (in->tableList, "AIPS BP");
    if (sel->BPversion==0) useVer = highVer;
    else useVer = sel->BPversion;
    cal->BPTable = 
      (Obit*)newObitTableBPValue (in->name, (ObitData*)in, &useVer, OBIT_IO_ReadOnly, 0, 0, 0, err);
    if (cal->BPTable==NULL) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: NO Bandpass BP table %d for %s", routine, useVer, in->name);
	return;
    }
  } /* end BP table setup */
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* FG table for Flagging */
  if (sel->doFlag) {
    /* if sel->FGversion ==0 use highest */
    highVer = ObitTableListGetHigh (in->tableList, "AIPS FG");
    if (sel->FGversion==0) useVer = highVer;
    else useVer = MIN (sel->FGversion, highVer);

    if (useVer>0) { /* have one - use */
      cal->FGTable = 
	(Obit*)newObitTableFGValue (in->name, (ObitData*)in, &useVer, OBIT_IO_ReadOnly, err);
      if (cal->FGTable==NULL) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: NO Flagging FG table %d for %s", routine, useVer, in->name);
	return;
      }
    } else {
      /* no flag table - ignore */
      sel->doFlag = FALSE;
      cal->FGTable = ObitTableFGUnref(cal->FGTable);
    }
  } /* end FG table setup */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Applying amp/phase/delay/rate calibration */
  if (sel->doCal) {

    /* Use CL table if one exists */
    highVer = ObitTableListGetHigh (in->tableList, "AIPS CL");
    if (highVer > 0) {
      /* if calVersion ==0 use highest */
      if (sel->calVersion==0) useVer = highVer;
      else useVer = sel->calVersion;
      cal->CLTable =
	(Obit*) newObitTableCLValue (in->name, (ObitData*)in, &useVer, OBIT_IO_ReadOnly, 0, 0, 0, err);
      if (cal->CLTable==NULL) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: NO calibration CL table %d for %s", routine, useVer, in->name);
	return;
      }
    } else {

      /* No CL table - Use SN table */
      highVer = ObitTableListGetHigh (in->tableList, "AIPS SN");
      if (sel->calVersion==0) useVer = highVer;
      else useVer = sel->calVersion;
      cal->SNTable =
	(Obit*) newObitTableSNValue (in->name, (ObitData*)in, &useVer, OBIT_IO_ReadOnly, 0, 0, err);
      if (cal->SNTable==NULL) { /* Couldn't open table */
	    Obit_log_error(err, OBIT_Error, 
		"%s: NO calibration SN table %d for %s", routine, useVer, in->name);
	    return;
      }
    } /* end get cal table */
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Get any CQ table */
    highVer = ObitTableListGetHigh (in->tableList, "AIPS CQ");
    useVer = 1;
    if (highVer>=1)
      cal->CQTable =
	(Obit*) newObitTableCQValue (in->name, (ObitData*)in, &useVer, OBIT_IO_ReadOnly, 0, err);
    else cal->CQTable = ObitTableCQUnref(cal->CQTable);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
  } /* end of tables setup for amp/phase/delay/rate calibration */

  /* Initialize indexing of the uv data if an NX table exists */
  highVer = ObitTableListGetHigh (in->tableList, "AIPS NX");
  useVer = 1;
  if (highVer>=1)
    sel->NXTable =
      (Obit*) newObitTableNXValue (in->name, (ObitData*)in, &useVer, OBIT_IO_ReadOnly, err);
  else sel->NXTable = ObitTableNXUnref(sel->NXTable);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  if (highVer>=1)
    ObitUVSelNextInit (sel, (ObitUVDesc*)in->myIO->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Need frequency information */
  ObitUVGetFreq (in, err);

  /* Start up calibration - finish output Descriptor, and Selector */
  ObitUVCalStart ((ObitUVCal*)in->myIO->myCal, (ObitUVSel*)in->myIO->mySel, 
		  (ObitUVDesc*)in->myIO->myDesc, in->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitUVSetupCal */

/**
 * Create myIO object depending on value of FileType in in->info.
 * This is the principle place where the underlying file type is known.
 * \param in   UV object to attach myIO
 * \param err  ObitErr for reporting errors.
 */
static void ObitUVSetupIO (ObitUV *in, ObitErr *err)
{
  ObitIOType FileType;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar *routine = "ObitUVSetupIO";

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
  if (FileType==OBIT_IO_FITS) {
    in->myIO = (ObitIO*)newObitIOUVFITS(in->name, in->info, err);
    /* copy selector pointer - use same object for UV and IO */
    ((ObitIOUVFITS*)in->myIO)->mySel = ObitUnref(((ObitIOUVFITS*)in->myIO)->mySel);
    ((ObitIOUVFITS*)in->myIO)->mySel = ObitRef(in->mySel);
    /* copy descriptor */
    ((ObitIOUVFITS*)in->myIO)->myDesc = ObitUVDescCopy(in->myDesc, 
		     ((ObitIOUVFITS*)in->myIO)->myDesc, err);

  } else if (FileType==OBIT_IO_AIPS) {
    in->myIO = (ObitIO*)newObitIOUVAIPS(in->name, in->info, err);
    /* copy selector pointer - use same object for UV and IO */
    ((ObitIOUVAIPS*)in->myIO)->mySel = ObitUnref(((ObitIOUVAIPS*)in->myIO)->mySel);
    ((ObitIOUVAIPS*)in->myIO)->mySel = ObitRef(in->mySel);
    /* copy descriptor */
    ((ObitIOUVAIPS*)in->myIO)->myDesc = ObitUVDescCopy(in->myDesc, 
		     ((ObitIOUVAIPS*)in->myIO)->myDesc, err);
  }

} /* end ObitUVSetupIO */

/**
 * Copy standard UV data tables and then selected tables with selection 
 * specified on inUV.
 * Tables with selection:
 * \li FQ
 * The FQ table is always copied selecting by IF and FQID
 * Frequency offsets relative to first selected IF.
 * \li AN
 * All AN Tables are copied with selection by IF and antenna.
 * If polarization cal (doPol) is being applied then polarization
 * calibration information is reset
 * \li SN
 * If no calibration (doCalib) is applied then all SN tables are copied
 * selecting by time, source, antenna, freqID, IF, poln
 * \li CL
 * If no calibration (doCalib) is applied then all CL tables are copied
 * selecting by time, source, antenna, freqID, IF, poln
 *
 * \param inUV  The object to copy from, defines selection
 * Uses source dependent frequency info if available from inUV->info
 * \li "SouIFOff" OBIT_double (nif,1,1) Source frequency offset per IF
 * \li "SouBW"    OBIT_double (1,1,1)   Bandwidth
 *
 * \param outUV An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
static ObitIOCode CopyTablesSelect (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS FQ","AIPS SU",
		    "AIPS AN","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS CQ",
		    "AIPS WX","AIPS AT","AIPS NI","AIPS BP","AIPS OF",
		    "AIPS PS",
		    "AIPS HI","AIPS PL","AIPS SL", 
		    NULL};
  gboolean copySU;
  odouble *SouIFOff=NULL, SouBW=0.0;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  gchar *routine = "ObitUV:CopyTablesSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Source info available? */
  ObitInfoListGetP    (inUV->info, "SouIFOff",  &type, dim, (gpointer)&SouIFOff);
  ObitInfoListGetTest (inUV->info, "SouBW",     &type, dim, &SouBW);
  
  /* Copy standard tables */
  retCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inUV->name, retCode);

  /* If multisource out then copy SU table, multiple sources selected or
   sources deselected suggest MS out */
  copySU = (inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources);
  /* Or all sources in a MS file selected */
  copySU = copySU || 
    ((inUV->mySel->numberSourcesList==0) && (inUV->myDesc->ilocsu>=0));
  if (copySU) retCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inUV->name, retCode);

  /* Copy FQ Table */
  retCode =  ObitTableFQSelect (inUV, outUV, SouIFOff, SouBW, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inUV->name, retCode);

  /* Copy AN Tables */
  retCode =  ObitTableANSelect (inUV, outUV, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inUV->name, retCode);

  /* Copy SN Tables */
  retCode =  ObitTableSNSelect (inUV, outUV, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inUV->name, retCode);

  /* Copy CL Tables */
  retCode =  ObitTableCLSelect (inUV, outUV, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inUV->name, retCode);

  return retCode;
} /* end CopyTablesSelect */
