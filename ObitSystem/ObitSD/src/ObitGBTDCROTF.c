/* $Id: ObitGBTDCROTF.c,v 1.5 2006/10/26 13:49:24 bcotton Exp $ */
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

#include "ObitGBTDCROTF.h"
#include "ObitTableGBTDCRSTATE.h"
#include "ObitTableGBTSPDATA.h"
#include "ObitFileFITS.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitGBTDCROTF.c
 * ObitGBTDCROTF class function definitions.
 *
 * This class has the facilities to convert from a GBT DCR file and
 * associated tables to OTF format.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitGBTDCROTF";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitGBTDCROTFClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGBTDCROTFClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGBTDCROTFInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGBTDCROTFClear (gpointer in);

/** Private:  Get pointing times from Antenna file, Secondary focus */
static void GetAntennaGR (ObitGBTDCROTF *in, gchar *scanName, olong inDisk, ObitErr *err);

/** Private:  Get pointing times from Antenna file, Prime Focus */
static void GetAntennaPF (ObitGBTDCROTF *in, gchar *scanName, olong inDisk, ObitErr *err);

/** Private:  Get pointing for a given time */
static void GetPoint (ObitGBTDCROTF *in, odouble time, ofloat *ra, ofloat *dec);

/** Private:  Get file descriptor */
static gboolean GetHeader (ObitGBTDCROTF *in, gchar *scanName, olong inDisk, ObitErr *err);

/** Private:  Get data */
static void GetData (ObitGBTDCROTF *in, gchar *infile, olong inDisk, ObitErr *err);

/** Private:  Get target id */
static olong GetTarget (ObitGBTDCROTF *in, gboolean isNew, gchar *name, ObitErr *err);

/** Private:  Initialize Index table */
static void InitScan (ObitGBTDCROTF *in, gboolean isNew, ofloat scan, 
		      ofloat target, ObitErr *err);

/** Private:  Undate scan info */
static void SetScan (ObitGBTDCROTF *in, ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitGBTDCROTFClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \param name  A name for the object
 * \return the new object.
 */
ObitGBTDCROTF* newObitGBTDCROTF (gchar* name)
{
  ObitGBTDCROTF* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTDCROTFClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitGBTDCROTF));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGBTDCROTFInit((gpointer)out);

  return out;
} /* end newObitGBTDCROTF */

/**
 * Constructor from values.
 * \param name     A name for the object
 * \param outOTF   Extant OTF to append new data to
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitGBTDCROTF* 
newObitGBTDCROTFValue (gchar *name, ObitOTF *outOTF, ObitErr *err)
{
  ObitGBTDCROTF* out=NULL;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;

  /* Create basic object */
  out = newObitGBTDCROTF(name);

  /* attach data */
  out->outOTF = ObitOTFRef(outOTF);

  return out;
} /* end newObitGBTDCROTFValue */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitGBTDCROTFGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTDCROTFClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetIOClass */

/**
 * Convert a GBT DCR and associated tables to OTF
 * \param in       The conversion object
 * \param inDisk   input FITS disk number for base of GBT files
 * \param scanName Scan name part of FITS file names (e.g. "2003_06_27_01:32:01")
 * \param err      Obit error stack object.
 * \return pointer to the new object.
 */
ObitIOCode ObitGBTDCROTFConvert (ObitGBTDCROTF *in, olong inDisk, gchar *scanName, 
				 ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitGBTDCROTFConvert";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));


  olong  nrec;
  ObitIOAccess access;
  gboolean isNew;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  /* Get header info, array geometry, initialize output if necessary */
  isNew = GetHeader (in, scanName, inDisk, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  
  /* Get Telescope position array, Prime and secondary focus different */
  if (in->refFrequency>1.0e9)
    GetAntennaGR (in, scanName, inDisk, err);
  else  /* Low frequencies at prime focus */
    GetAntennaPF (in, scanName, inDisk, err);
  /* Say what went wrong if error */
  if (err->error) 
     Obit_log_error(err, OBIT_Error, "Error reading Antenna file for scan %s", scanName);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  
  Obit_log_error(err, OBIT_InfoErr, "Adding scan %s %6.0f to OTF file", scanName, in->scan);

  /* Define I/O size */
  dim[0] = 1;
  nrec = in->StateData->nDCRState;
  ObitInfoListPut (in->outOTF->info, "nRecPIO", OBIT_long, dim, &nrec, err);
  ObitInfoListPut (in->outOTF->info, "IOBy",    OBIT_long, dim, &nrec, err);

  /* Open output OTF */
  if (isNew) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  retCode = ObitOTFOpen (in->outOTF, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file");
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Get target id number */
  in->target = (ofloat)GetTarget (in, isNew, in->Name, err);

  /* Initialize scan in Index table */
  InitScan (in, isNew, in->scan, in->target, err);

  /* convert data  */
  GetData (in, scanName, inDisk, err);

  /* Update index table */
  SetScan (in, err);

  /* Close */
  retCode = ObitOTFClose (in->outOTF, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  
   /* cleanup */
   in->nAntTime = 0;
   if (in->AntDMJD) g_free(in->AntDMJD); in->AntDMJD = NULL;
   if (in->AntRA)   g_free(in->AntRA);   in->AntRA = NULL;
   if (in->AntDec)  g_free(in->AntDec);  in->AntDec = NULL;
   in->IFdata      = ObitGBTIFInfoUnref(in->IFdata); 
   in->BeamOffData = ObitGBTBeamOffInfoUnref(in->BeamOffData); 
   in->StateData   = ObitGBTDCRStateInfoUnref(in->StateData); 
   in->StateData   = ObitGBTDCRStateInfoUnref(in->StateData); 

   return OBIT_IO_OK;
} /* end ObitGBTDCROTFConvert */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitGBTDCROTFClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGBTDCROTFClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGBTDCROTFClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGBTDCROTFClassInfoDefFn (gpointer inClass)
{
  ObitGBTDCROTFClassInfo *theClass = (ObitGBTDCROTFClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGBTDCROTFClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGBTDCROTFClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGBTDCROTFGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitGBTDCROTFClear;
  theClass->ObitInit      = (ObitInitFP)ObitGBTDCROTFInit;
  theClass->newObit       = (newObitFP)newObitGBTDCROTF;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;

} /* end ObitGBTDCROTFClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitGBTDCROTFInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGBTDCROTF *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->outOTF = NULL;
  in->nAntTime = 0;
  in->AntDMJD = NULL;
  in->AntRA = NULL;
  in->AntDec = NULL;
  in->IFdata = NULL;
  in->StateData = NULL;
  in->BeamOffData = NULL;
  in->nfeed = 1;
  in->nchan = 1;
  in->nstok = 2;
  in->isBS = FALSE;

} /* end ObitGBTDCROTFInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitGBTDCROTFClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGBTDCROTF *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  if (in->AntDMJD) g_free(in->AntDMJD); in->AntDMJD = NULL;
  if (in->AntRA)   g_free(in->AntRA);   in->AntRA = NULL;
  if (in->AntDec)  g_free(in->AntDec);  in->AntDec = NULL;
  in->IFdata      = ObitGBTIFInfoUnref(in->IFdata); 
  in->BeamOffData = ObitGBTBeamOffInfoUnref(in->BeamOffData); 
  in->StateData   = ObitGBTDCRStateInfoUnref(in->StateData); 
  in->outOTF      = ObitOTFUnref(in->outOTF);
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitGBTDCROTFClear */

/**
 * Get antenna positions and leave in globals, Secondary focus (>1 GHz)
 * \param  in       The conversion object
 * \param  scanName   root of input file names 
 * \param  inDisk   FITS disk number for input
 * \param  err       Obit return error stack   
 */
static void GetAntennaGR (ObitGBTDCROTF *in, gchar *scanName, olong inDisk, ObitErr *err)
{
  gchar FullFile[128];
  olong irow;
  olong ver, nrow;
  gchar *tab;
  ObitTableGBTANTPOSGR* table;
  ObitTableGBTANTPOSGRRow* row;
  ObitIOCode retCode;
  gchar *routine = "GetAntennaGR";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(scanName!=NULL);

  /* get full file name */
  sprintf (FullFile,"Antenna/%s.fits", scanName);
  
  /* Create table structure */
  table = newObitTableGBTANTPOSGR("Antenna");

  /* Setup */
  tab = "ANTPOSGR";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(table,inDisk,FullFile,tab,ver,nrow,err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Open */
  retCode = ObitTableGBTANTPOSGROpen (table, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create Row structure */
  row = newObitTableGBTANTPOSGRRow (table);

  /* make sure there is data */
  if (table->myDesc->nrow<=0) {
     Obit_log_error(err, OBIT_Error, "No data in Antenna file for scan %s", scanName);
     return;
 }

  /* Create arrays */
  in->nAntTime = table->myDesc->nrow;
  in->AntDMJD = g_malloc0(in->nAntTime*sizeof(odouble));
  in->AntRA   = g_malloc0(in->nAntTime*sizeof(odouble));
  in->AntDec  = g_malloc0(in->nAntTime*sizeof(odouble));

  /* Loop over table */
  for (irow = 1; irow<=in->nAntTime; irow++) {
    retCode = ObitTableGBTANTPOSGRReadRow (table, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    in->AntDMJD[irow-1] = row->dmjd;
    in->AntRA[irow-1]   = row->raj2000;
    in->AntDec[irow-1]  = row->decj2000;
 } /* end loop over table */
  
  /* Close */
  retCode = ObitTableGBTANTPOSGRClose (table, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  table = ObitTableGBTANTPOSGRUnref(table);
  row = ObitTableGBTANTPOSGRRowUnref(row);
} /* end GetAntennaGR  */

/**
 * Get antenna positions and leave in globals, Prime focus (<1 GHz)
 * \param  in       The conversion object
 * \param  outData  Output OTF object
 * \param  scanName Name of scan
 * \param  inDisk   FITS disk number for input
 * \param  err      Obit return error stack
 */
static void GetAntennaPF (ObitGBTDCROTF *in, gchar *scanName, olong inDisk, ObitErr *err)
{
  gchar FullFile[128];
  olong irow;
  olong ver, nrow;
  gchar *tab;
  ObitTableGBTANTPOSPF* table;
  ObitTableGBTANTPOSPFRow* row;
  ObitIOCode retCode;
  gchar *routine = "GetAntennaPF";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(scanName!=NULL);

  /* get full file name */
  sprintf (FullFile,"Antenna/%s.fits", scanName);
  
  /* Create table structure */
  table = newObitTableGBTANTPOSPF("Antenna");

  /* Setup */
  tab = "ANTPOSPF";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(table,inDisk,FullFile,tab,ver,nrow,err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Open */
  retCode = ObitTableGBTANTPOSPFOpen (table, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create Row structure */
  row = newObitTableGBTANTPOSPFRow (table);

  /* make sure there is data */
  if (table->myDesc->nrow<=0) {
     Obit_log_error(err, OBIT_Error, "No data in Antenna file for scan %s", scanName);
     return;
 }

  /* Create arrays */
  in->nAntTime = table->myDesc->nrow;
  in->AntDMJD = g_malloc0(in->nAntTime*sizeof(odouble));
  in->AntRA   = g_malloc0(in->nAntTime*sizeof(odouble));
  in->AntDec  = g_malloc0(in->nAntTime*sizeof(odouble));

  /* Loop over table */
  for (irow = 1; irow<=in->nAntTime; irow++) {
    retCode = ObitTableGBTANTPOSPFReadRow (table, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    in->AntDMJD[irow-1] = row->dmjd;
    in->AntRA[irow-1]   = row->raj2000;
    in->AntDec[irow-1]  = row->decj2000;
 } /* end loop over table */
  
  /* Close */
  retCode = ObitTableGBTANTPOSPFClose (table, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  table = ObitTableGBTANTPOSPFUnref(table);
  row = ObitTableGBTANTPOSPFRowUnref(row);
} /* end GetAntennaPF  */

/**
 * Get header information from scan header files
 * \param  in       The conversion object
 * \param  outData  Output OTF object       
 * \param  scanName   root of input file names 
 * \param  inDisk   FITS disk number for input
 * \param  outDisk  FITS disk number for output
 * \param  err      Obit return error stack
 * \return  TRUE if the file is just created  
 */
static gboolean GetHeader (ObitGBTDCROTF *in, gchar *scanName, olong inDisk, ObitErr *err)
{
  gboolean out = FALSE;
  ObitOTFDesc *desc;
  ObitOTFArrayGeom *geom;
  ObitIOCode retCode;
  olong numberDetect, ndetec, ncol, ncopy, iscan;
  gchar Date[48], *fullname=NULL, *tempStr, FullFile[128], commnt[81];
  odouble SiteLat, SiteLong, SiteElev, JD, T, GMST0, GSTIAT0, nu, e2;
  odouble aEarth = 6378137.0; /* Semi major axis of Earth */
  odouble flat   = 1.0 / 298.257223563; /* inverse of flattening of Earth */ 
  ofloat UT1UTC, PolarX, PolarY;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  olong i, j, outDisk, ierr = 0;
  gchar backend[21];
  ObitFileFITS *FITS;
  gchar *routine = "GetHeader";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;
  g_assert(scanName!=NULL);

  desc = in->outOTF->myDesc;

  /* Get frequency setup */
  in->IFdata = newObitGBTIFInfoValue("IFdata", "DCR", inDisk, scanName, err);
  /* Save data */
  in->nchan        = in->IFdata->nchan;
  in->deltaFreq    = in->IFdata->delta[0];
  in->refPixel     = in->IFdata->refPixel[0];
  in->refFrequency = in->IFdata->refFrequency[0];
  in->poln[0]      = in->IFdata->poln[0];
  in->poln[1]      = in->IFdata->poln[1];

  /* Get Beam offset data */
  in->BeamOffData = newObitGBTBeamOffInfoValue("BeamOff", inDisk, scanName, err);

  /* Get State info data */
  in->StateData = newObitGBTDCRStateInfoValue("State", inDisk, scanName, err);

  /* Is this regular of beam switched? */
  in->isBS = (in->IFdata->srfeed1[0] > 0) || (in->IFdata->srfeed2[0] > 0);

  /* Set Feed Array Geometry */
  /* create Array geometry with 2 (R,L) elements */
  numberDetect = 2;
  if (in->isBS) numberDetect = 8; /* 2 real feeds x 2 poln x 2 states */
  geom = ObitOTFArrayGeomCreate (numberDetect);
  geom->azOffset = g_realloc(geom->azOffset, numberDetect*sizeof(ofloat));
  geom->elOffset = g_realloc(geom->elOffset, numberDetect*sizeof(ofloat));

  /* Init - for unswitched */
  in->azOff = 0.0; /* no position offset for single feed */
  in->elOff = 0.0; /* no position offset for single feed */
  for (i=0; i<numberDetect; i++) {
    geom->azOffset[i] = 0.0;
    geom->elOffset[i] = 0.0;
 }

  /* If beam switched get offsets from FITS file */
  if (in->isBS) {
    for (i=0; i<numberDetect; i+=2) {
      j = i/2;  /* feed number */
      if (i>=numberDetect/2) j = (i-numberDetect/2)/2;
      geom->azOffset[i] = -in->BeamOffData->xeloff[j];
      geom->elOffset[i] =  in->BeamOffData->eloff[j];
      /* Poln pair has the same offset */
      geom->azOffset[i+1] = -in->BeamOffData->xeloff[j];
      geom->elOffset[i+1] =  in->BeamOffData->eloff[j];
      /* Set polarizations */
      j *= 2;
      in->poln[i]   = in->IFdata->poln[j];
      in->poln[i+1] = in->IFdata->poln[j+1];
    }
    /* rerefer all the offsets to feed 1 */
    /* Correction to pointing */
    in->azOff = -geom->azOffset[0];
    in->elOff = -geom->elOffset[0];
    for (i=0; i<numberDetect; i++) {
      geom->azOffset[i] += in->azOff;
      geom->elOffset[i] += in->elOff;
    }
 }
  
  /* Other information - Get from Antenna file */
  /* get full file name */
  sprintf (FullFile,"Antenna/%s.fits", scanName);

  /* Get header keywords  */
  FITS = newObitFileFITS("FITS_file");
  retCode = ObitFileFITSOpen (FITS, FullFile, inDisk, OBIT_IO_ReadOnly,err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  
  /* Read main file keywords */
  UT1UTC = 0.0;
  retCode = ObitFileFITSReadKeyFlt (FITS, "DELTAUTC", &UT1UTC, commnt, err);
  if (retCode == OBIT_IO_NotFoundErr) retCode = OBIT_IO_OK;
  
  PolarX = 0.0;
  retCode = ObitFileFITSReadKeyFlt (FITS, "IERSPMX", &PolarX, commnt,  err);
  if (retCode == OBIT_IO_NotFoundErr) retCode = OBIT_IO_OK;
  
  PolarY = 0.0;
  retCode = ObitFileFITSReadKeyFlt (FITS,"IERSPMY", &PolarY, commnt,  err);
  if (retCode == OBIT_IO_NotFoundErr) retCode = OBIT_IO_OK;
  
  SiteLat = 3.8433119e+01;
  retCode = ObitFileFITSReadKeyDbl (FITS, "SITELAT", &SiteLat, commnt, err);
  if (retCode == OBIT_IO_NotFoundErr) retCode = OBIT_IO_OK;
  
  SiteLong = 7.9839833E+01 ;
  retCode = ObitFileFITSReadKeyDbl (FITS, "SITELONG", &SiteLong, commnt, err);
  if (retCode == OBIT_IO_NotFoundErr) retCode = OBIT_IO_OK;
  
  SiteElev = 8.24595E+02;
  retCode = ObitFileFITSReadKeyDbl (FITS, "SITEELEV", &SiteElev, commnt, err);
  if (retCode == OBIT_IO_NotFoundErr) retCode = OBIT_IO_OK;
  
  iscan = 1;
  retCode = ObitFileFITSReadKeyLng (FITS, "SCAN", &iscan, commnt, err);
  if (retCode == OBIT_IO_NotFoundErr) retCode = OBIT_IO_OK;
  
  retCode = ObitFileFITSReadKeyStr (FITS, "DATE-OBS", Date, commnt, err);
  if (retCode == OBIT_IO_NotFoundErr) retCode = OBIT_IO_OK;
  
  retCode = ObitFileFITSReadKeyStr (FITS, "OBJECT", in->Name, commnt, err);
  if (retCode == OBIT_IO_NotFoundErr) retCode = OBIT_IO_OK;
  
  /* Close FITS file */
  retCode = ObitFileFITSClose (FITS, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  
  /* Convert to JD */
  ObitOTFDescDate2JD (Date, &JD);
  
  /* GST at IAT=0 at 0h on reference date (deg)*/
  /* Tropical century from jan 0.5, 2000 */
  T = (JD - 2451545.0) / 36525.0;
  
  /* GMST at IAT=0 in radians */
  GMST0 = ((((((-6.2e-6 * T) + 0.093104) * T) + 8640184.812866) * T + 24110.54841) 
	   * 2.0 * G_PI / 86400.0);
  /* to degrees */
  GSTIAT0 = RAD2DG * fmod(GMST0, (2.0*G_PI));

   ncopy = strlen (Date);
   if (ncopy>10) ncopy = 10;
  for (i=0; i<ncopy; i++) geom->RefDate[i] = Date[i]; geom->RefDate[i]=0;
  geom->TimeSys[0] = 'U'; geom->TimeSys[1] = 'T';geom->TimeSys[2] = 'C';geom->TimeSys[3] = 0;

  /* Conversion of geocentric to geodetic */
  e2 = 2.0*flat - flat*flat;
  nu = aEarth / sqrt (1.0 - e2 * sin(DG2RAD*SiteLat) * sin(DG2RAD*SiteLat));
  geom->TeleX =  (nu + SiteElev) * cos(DG2RAD*SiteLat) * cos(DG2RAD*SiteLong);
  geom->TeleY = -(nu + SiteElev) * cos(DG2RAD*SiteLat) * sin(DG2RAD*SiteLong);
  geom->TeleZ =  ((1.0 - e2) * nu +  SiteElev) * sin(DG2RAD*SiteLat);


  geom->DegDay  = 3.6098564497330e+02;
  geom->GSTiat0 = GSTIAT0;
  geom->PolarX  = PolarX;
  geom->PolarY  = PolarY;
  geom->ut1Utc  = UT1UTC;
  geom->dataUtc = 0.0;
  geom->iatUtc  = 0.0;

  /* Compute some useful terms */
  /* telescope latitude in radians */
  geom->lat = SiteLat * RAD2DG;
  /* telescope longitude in radians */
  geom->lon = SiteLong * RAD2DG;
  /* LST at iat0 in radians */
  geom->LSTiat0 = geom->GSTiat0*1.74533e-2 + geom->lon;
  /* Earth rotation rate in rad/day */
  geom->RadDay = geom->DegDay*1.74533e-2;
  /* Data - IAT in days */
  geom->dataIat = (geom->dataUtc - geom->iatUtc) / 86400.0;

  /* Attach Array geometry to OTF */
  in->outOTF->geom = ObitOTFArrayGeomRef(geom);
  geom = ObitOTFArrayGeomUnref(geom);

  /* Other information from DCR file */
  /* get full file name */
  sprintf (FullFile,"DCR/%s.fits", scanName);

  /* Get header keywords */
  retCode = ObitFileFITSOpen (FITS, FullFile, inDisk, OBIT_IO_ReadOnly,err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Make sure data in "DCR" mode */
  retCode = ObitFileFITSReadKeyStr (FITS, "BACKEND", backend, commnt, err);
  if (retCode != OBIT_IO_OK) {
     Obit_log_error(err, OBIT_Error, "%s: No BACKEND keyword in DCR file", routine);
     return retCode;
  }
  if (strncmp (backend, "DCR",3)) {
    Obit_log_error(err, OBIT_Error, "ERROR backend %s NOT DCR", backend);
    return retCode;
  }
  
  /* Read main file keywords */
  ndetec = 1;
  retCode = ObitFileFITSReadKeyLng (FITS, "NRCVRS", &ndetec, commnt, err);
  if (retCode == OBIT_IO_NotFoundErr) retCode = OBIT_IO_OK;
  in->nfeed = ndetec/2;

  in->integTime = 0.0; /* Integration time in days */
  retCode = ObitFileFITSReadKeyDbl (FITS, "DURATION", &in->integTime, commnt, err);
  if (retCode == OBIT_IO_NotFoundErr) retCode = OBIT_IO_OK;
  if (ierr==KEY_NO_EXIST) ierr = 0;
  /* Convert from seconds to days */
  in->integTime /= 86400.0;
  
  /* Close FITS file */
  retCode = ObitFileFITSClose (FITS, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  
  /* initialize globals */
  in->refMJD = 0.0;
  in->target = 0.0;
  in->scan   = (ofloat)iscan;

   
  /* Does the output file exist?  BTW, I only do FITS */
  if(!ObitInfoListGet(in->outOTF->info, "Disk", &type, dim, &outDisk, err))
    Obit_traceback_val (err, routine, in->name, retCode);
  if(!ObitInfoListGet(in->outOTF->info, "FileName", &type, dim, FullFile, err))
    Obit_traceback_val (err, routine, in->name, out);
  FullFile[dim[0]] = 0; /* NULL terminate */
  /* Drop any initial "!" */
  tempStr = FullFile;
  if (FullFile[0] == '!') tempStr = FullFile+1;
  fullname = ObitFITSFilename (outDisk, tempStr, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /* check is output file exists, if not, initialize */
  if (!ObitFileExist (fullname, err)) {
    if (err->error>0) Obit_log_error(err, OBIT_Error, "ERROR testing file %s", fullname);

    /* Initialize  Descriptor */
    /* Approximate Beam size */
    desc->beamSize = 2.0e8 / in->refFrequency;

    /* diameter */
    desc->diameter = 100.0;  /* GBT diameter */

    /* new file */
    out = TRUE;
    
    strncpy (desc->object, "Sky", OTFLEN_VALUE);
    strncpy (desc->teles,  "GBT       ", OTFLEN_VALUE);
    strncpy (desc->origin, "Obit ", OTFLEN_VALUE);
    desc->isort[0] = 'T';  /* Time ordered */
    ncol = 0;
    
    desc->JDObs = 0.0;
    desc->epoch = 2000.0;
    desc->equinox = 2000.0;
    strncpy (desc->bunit,  "ADU     ", OTFLEN_VALUE);
    strncpy (desc->obsdat, Date, OTFLEN_VALUE);
    
    /* Time */
    strncpy (desc->colType[ncol], "TIME    ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DAYS    ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Integration time */
    strncpy (desc->colType[ncol], "TIME_INT", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DAYS    ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Target index */
    strncpy (desc->colType[ncol], "TARGET  ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "        ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Scan index */
    strncpy (desc->colType[ncol], "SCAN    ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "        ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Pointing RA */
    strncpy (desc->colType[ncol], "RA      ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DEGREE  ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Pointing Dec */
    strncpy (desc->colType[ncol], "DEC     ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DEGREE  ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Rotation of array on sky */
    strncpy (desc->colType[ncol], "ROTATE  ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DEGREE  ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Cal on? */
    strncpy (desc->colType[ncol], "CAL     ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "        ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Data - MUST be last column */
    desc->numDesc = ncol;
    strncpy (desc->colType[ncol], "DATA    ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "COUNTS  ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = numberDetect*2; /* With weights */
    ncol++;
    
    desc->ncol = ncol;
    
    /* Data array descriptors */
    desc->naxis = 0;

    /* Stokes axis */
    desc->inaxes[desc->naxis] = in->nstok;
    strncpy (desc->ctype[desc->naxis], "STOKES  ", OTFLEN_KEYWORD);
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    /* Set appropritate Stokes Type */
    if ((in->poln[0]=='R') || (in->poln[0]=='L')) { /* R, L */
      desc->cdelt[desc->naxis] = -1.0;
      desc->crval[desc->naxis] = -1.0;
    } else { /* Must be X, Y */
      desc->cdelt[desc->naxis] = -1.0;
      desc->crval[desc->naxis] = -5.0;
    }
    desc->naxis++;

    /* Data-Wt axis */
    desc->inaxes[desc->naxis] = 2;
    strncpy (desc->ctype[desc->naxis], "DATAWT", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = 1.0;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = 1.0;
    desc->naxis++;

    /* Feed axis */
    desc->inaxes[desc->naxis] = in->nfeed;
    if (in->isBS) desc->inaxes[desc->naxis] *= 2; /* Need extra "feeds" */
    strncpy (desc->ctype[desc->naxis], "FEED", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = 1.0;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = 1.0;
    desc->naxis++;

    /* Frequency axis */
    desc->inaxes[desc->naxis] = in->nchan;
    strncpy (desc->ctype[desc->naxis], "FREQ    ", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = 1.0;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = in->refFrequency;
    desc->naxis++;

    /* Reference Modified Julian date */
    in->refMJD = JD - 2400000.5;

  } /* end initialize descriptor */

  /* cleanup */
  if (fullname) g_free(fullname);
  FITS = ObitFileFITSUnref(FITS);

  /* Index the descriptor */
  ObitOTFDescIndex (desc);

  return out;
} /* end GetHeader */

/**
 * Read data from GB FITS files and write  
 * Assumptions:
 *    The following uses state information read from the FITS files to
 * decide what data is what. 
 * If the data is beam switched, a second set of "detectors is used 
 * for the "ref" state where the feeds are swapped and therefore have
 * different amplifier chains and therefore different gains
 * \param  in         The conversion object
 * \param  scanName   Scan part of input file name
 * \param  inDisk     FITS disk number for input
 * \param  err        Obit return error stack   
 */
static void GetData (ObitGBTDCROTF *in, gchar *scanName, olong inDisk, ObitErr *err)
{
  gchar FullFile[128];
  olong irow, numberDetect;
  olong i, ver, nrow, ip1off, ip2off, feed, stOff;
  olong nS, nIF, nPoln, iS, iIF, ii, jj, incdatawt;
  gchar *tab;
  ofloat *data, ra, dec, fblank = ObitMagicF();
  odouble timeCorrD, TimePhase[10], total, JD;
  ObitTableGBTDCRDATA* table;
  ObitTableGBTDCRDATARow* row;
  ObitOTFDesc *desc;
  ObitOTFArrayGeom *geom;
  ObitIOCode retCode;
  gchar *routine = "GetData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(scanName!=NULL);

  numberDetect = in->outOTF->geom->numberDetect; /* number of detectors */
  nS = in->StateData->nDCRState;   /* Number of phase states */
  nIF = in->IFdata->nIF;           /* Number of IFs */
  nPoln = 2;                       /* Number of polarizations */

  /* Stokes offsets in input */
  ip1off = 0;
  if (in->poln[0]=='R') ip1off = 0;
  else if (in->poln[0]=='L') ip1off = 1;
  else if (in->poln[0]=='X') ip1off = 0;
  else if (in->poln[0]=='Y') ip1off = 1;

  ip2off = 0;
  if (in->poln[1]=='R') ip2off = 0;
  else if (in->poln[1]=='L') ip2off = 1;
  else if (in->poln[1]=='X') ip2off = 0;
  else if (in->poln[1]=='Y') ip2off = 1;

  /* Set offsets (sec) to start time time in each phase */
  TimePhase[0] = 0.5 * in->StateData->phasetim[0];
  total = in->StateData->phasetim[0];
  for (iS=1; iS<nS; iS++) {
    TimePhase[iS] = TimePhase[iS-1] + 0.5 * in->StateData->phasetim[iS-1] +
      0.5 * in->StateData->phasetim[iS];
    total += in->StateData->phasetim[iS];  /* accumulate total time */
  }

  /* Convert time offsets to center time in days */
  for (iS=0; iS<nS; iS++) {
    TimePhase[iS] -=  0.5 * total;
    TimePhase[iS] /=  86400.0;
  }

  /* Get reference MJD, convert ref date  to JD */
  ObitOTFDescDate2JD (in->outOTF->geom->RefDate, &JD);
  /* Reference Modified Julian date */
  in->refMJD = JD - 2400000.5;

  /* get full file name */
  sprintf (FullFile,"DCR/%s.fits", scanName);
  
  /* Create table structure */
  table = newObitTableGBTDCRDATA("Data");
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Setup */
  tab = "DATA";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(table,inDisk,FullFile,tab,ver,nrow,err);
 if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Correction to time in days */
  timeCorrD = in->timeCorr / 86400.0;

  /* Open */
  retCode = ObitTableGBTDCRDATAOpen (table, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create Row structure */
  row = newObitTableGBTDCRDATARow (table);

  desc = in->outOTF->myDesc;
  geom = in->outOTF->geom;
  incdatawt = desc->incdatawt;  /* Increment on Data-Wt axis */

  /* Write at end */
  desc->firstRec = desc->nrecord+1;
  in->startRec = desc->firstRec;  /* First record in scan */
  in->startTime = -1.0e20;        /* dummy srart time */

  /* Loop over table */
  for (irow = 1; irow<=table->myDesc->nrow; irow++) {
    retCode = ObitTableGBTDCRDATAReadRow (table, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* first time in scan */
    if (in->startTime<-1000.0) in->startTime = row->timetag - in->refMJD;
    data = in->outOTF->buffer;     /* Output data array */

    /* Loop over state */
    stOff = 0;
    for (iS=0; iS<nS; iS++) {
      /* Fill record entry in data */
      data[desc->iloct]   = row->timetag - in->refMJD;    /* Time (days) */
      /* Correct for position in integration */
      data[desc->iloct] += TimePhase[iS];
      /* Correction for GBT timing error  */
      data[desc->iloct] += timeCorrD;
      /* Get position */
      GetPoint (in, row->timetag+TimePhase[iS]+timeCorrD, &ra, &dec);
      data[desc->ilocti]  = 0.0;          /* time interval (days) */
      data[desc->iloctar] = in->target;   /* target number */
      data[desc->ilocscan]= in->scan;     /* scan number */
      data[desc->iloccal] = 0.0;          /* Cal? */
      if (in->StateData->cal[iS]) data[desc->iloccal] = 1.0; 
      /* Parallactic angle (deg) */
      data[desc->ilocrot] = ObitOTFArrayGeomParAng(geom, data[desc->iloct], ra, dec);
      /* correction to position */
      ObitOTFArrayGeomCorrPoint(in->azOff, in->elOff, data[desc->ilocrot], &ra, &dec);
      data[desc->ilocra]  = ra;           /* RA  in deg. */
      data[desc->ilocdec] = dec;          /* Dec in deg.*/

      /* Init data to blank */	
      for (i=0; i<numberDetect; i++) {
	data[desc->ilocdata+i*incdatawt]   = fblank;
	data[desc->ilocdata+i*incdatawt+1] = 0.0;
      }

      /* sig or ref state? */      
      if (in->StateData->sigref[iS]) { 
	/* "ref" state - use second set of ficticious feeds */
	/* Loop over "IFs" */
	for (iIF=0; iIF<nIF; iIF+=2) {
	  feed = MAX (1, in->IFdata->feed[iIF]) - 1;
	  jj = nPoln * feed + numberDetect/2;
	  /* Swap 1 and 2 on input */
	  ii = ((nIF-iIF-2)+ip1off) * nS;
	  data[desc->ilocdata+(0+jj)*incdatawt]   =  row->data[iS+ii]; /* RCP or X*/
	  data[desc->ilocdata+(0+jj)*incdatawt+1] =  1.0; /* Initial weight */
	  ii = ((nIF-iIF-2)+ip2off) * nS;
	  data[desc->ilocdata+(1+jj)*incdatawt]   =  row->data[iS+ii]; /* LCP or Y */
	  data[desc->ilocdata+(1+jj)*incdatawt+1] =  1.0; /* Initial weight */
	}
      } else {
	/* "sig" state - everything is cool */
 	for (iIF=0; iIF<nIF; iIF+=2) { /* do poln pairs */
	  feed = MAX (1, in->IFdata->feed[iIF]) - 1;
	  ii = (iIF+ip1off) * nS;
	  jj = nPoln * feed;
	  data[desc->ilocdata+(0+jj*incdatawt)]   = row->data[iS+ii]; /* RCP or X*/
	  data[desc->ilocdata+(0+jj*incdatawt+1)] = 1.0; /* Initial weight */
	  ii = (iIF+ip2off) * nS;
	  data[desc->ilocdata+(1+jj)*incdatawt]   = row->data[iS+ii]; /* LCP or Y */
	  data[desc->ilocdata+(1+jj)*incdatawt+1] = 1.0; /* Initial weight */
	}
      }
      data += desc->lrec;  /* next record */
    } /* End loop over state */

    /* set number of records */
    desc->numRecBuff = in->StateData->nDCRState;
    
    /* Write output  */
    if ((ObitOTFWrite (in->outOTF, NULL, err) != OBIT_IO_OK) || (err->error>0))
      Obit_log_error(err, OBIT_Error, "ERROR writing output Table file");
    if (err->error) Obit_traceback_msg (err, routine, in->name);
 } /* end loop over table */
  
  /* Get end times and record numbers */
  desc->firstRec = desc->nrecord+1;
  in->endTime = row->timetag - in->refMJD;
  in->endRec  = desc->nrecord;  /* Last record in scan */

  /* Close */
  retCode = ObitTableGBTDCRDATAClose (table, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Cleanup */
  table = ObitTableGBTDCRDATAUnref(table);
  row = ObitTableGBTDCRDATARowUnref(row);

} /* end GetData  */

/**
 * Get antenna pointing at a given time, end points or linear interpolation   .
 * \param  in       The conversion object
 * \param  time   [in] The desired MJD 
 * \param  ra     [out] RA J2000 in degrees 
 * \param  dec    [out]Dec J2000 in degrees
 */
static void GetPoint (ObitGBTDCROTF *in, odouble time, ofloat *ra, ofloat *dec)
{
  olong i, best;
  odouble test, delta;
  ofloat w1, w2;
  
  /* Find closest */
  best = -1;
  delta = 1.0e20;
  for (i=0; i<in->nAntTime; i++) { /* loop over array */
    test = fabs (in->AntDMJD[i]-time);
      if (delta> test) {
	delta = test;
	best = i;
      } else { /* must be getting further, stop */
	break;
      }
  }

  /* end points */
  if ((best==0) || (best==(in->nAntTime-1))) {
    *ra = in->AntRA[best];
    *dec = in->AntDec[best];
  } else if (time<in->AntDMJD[best]){ /* interpolate with previous */
    w1 = (in->AntDMJD[best]-time) / (in->AntDMJD[best]-in->AntDMJD[best-1]);
    w2 = 1.0 - w1;
    *ra  = w1 * in->AntRA[best-1]  + w2 * in->AntRA[best];
    *dec = w1 * in->AntDec[best-1] + w2 * in->AntDec[best];
  } else { /* interpolate with following */
    w1 = (in->AntDMJD[best+1]-time) / (in->AntDMJD[best+1]-in->AntDMJD[best]);
    w2 = 1.0 - w1;
    *ra  = w1 * in->AntRA[best]  + w2 * in->AntRA[best+1];
    *dec = w1 * in->AntDec[best] + w2 * in->AntDec[best+1];
  }
} /* end GetPoint */

/**
 * Get target id, look through existing table create new entry if needed.
 * \param  in       The conversion object
 * \param  isNew    True if output file just created
 * \param  name     Name of target
 * \param  err       Obit return error stack
 * \return Target id 
 */
static olong GetTarget (ObitGBTDCROTF *in, gboolean isNew, gchar *name, ObitErr *err)
{
  olong targ = -1;
  ObitTableOTFTarget* table;
  ObitTableOTFTargetRow* row;
  olong iRow, ver;
  gboolean doWrite;
  ObitIOAccess access;
  gchar *routine = "GetTarget";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return targ;
  g_assert(name!=NULL);

  /* create Target table object */
  ver = 1;
  if (isNew) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  table = newObitTableOTFTargetValue ("Target table", (ObitData*)in->outOTF, &ver, access, err);
  if (err->error) Obit_traceback_val (err, routine, in->outOTF->name, targ);

  /* Open table */
  if ((ObitTableOTFTargetOpen (table, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFTarget table");
    return targ;
  }

  /* Create Row */
  row = newObitTableOTFTargetRow (table);

  /* attach to table buffer */
  ObitTableOTFTargetSetRow (table, row, err);
  if (err->error) Obit_traceback_val (err, routine, in->outOTF->name, targ);

  /* Newly created?  Just write new one */
  doWrite = FALSE;
  if (isNew) {
    targ = 1;
    row->TargID = targ;
    strncpy(row->Target, name, 16);
    doWrite = TRUE;
  } else { /* Existing, see if already exists? */

    /* loop through table */
    for (iRow = 1; iRow<=table->myDesc->nrow; iRow++) {
      if ((ObitTableOTFTargetReadRow (table, iRow, row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading OTFTarget Table file");
	return targ;
      }
      if (!strncmp (row->Target, name, 16)) {
	/* Found match */
	targ = row->TargID;
	break;
      }  
    } /* end loop over table */

    /* Add new entry? */
    if (targ<=0) {
      targ = table->myDesc->nrow + 1;
      row->TargID = targ;
      strncpy(row->Target, name, 16);
      doWrite = TRUE;
    }
  } /* end output table already exists */

  /* need to write new entry? */
  if (doWrite) {
    iRow = table->myDesc->nrow + 1;
    if ((ObitTableOTFTargetWriteRow (table, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR writing OTFTarget Table file");
      return targ;
    }
  }
  
 /* Close  table */
  if ((ObitTableOTFTargetClose (table, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFTarget Table file");
    return targ;
  }

  /* Cleanup */
  row = ObitTableOTFTargetRowUnref(row);
  table = ObitTableOTFTargetUnref(table);

  return targ;
} /* end  GetTarget */

/**
 * Initializes Index table for this scan creating an entry
 * \param  in       The conversion object
 * \param  isNew    True if output file just created 
 * \param  scan     Scan number
 * \param  target   Target ID number
 * \param  err      Obit return error stack    
 */
static void InitScan (ObitGBTDCROTF *in, gboolean isNew, ofloat scan, 
	       ofloat target, ObitErr *err)
{
  ObitTableOTFIndex* table;
  ObitTableOTFIndexRow* row;
  olong iRow, ver;
  olong scanID, targetID;
  ObitIOAccess access;
  gchar *routine = "InitScan";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* create Index table object */
  scanID = (olong)(scan+0.5);
  targetID = (olong)(target+0.5);
  ver = 1;
  if (isNew) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  table = newObitTableOTFIndexValue ("Index table", (ObitData*)in->outOTF, &ver, access, err);
  if (err->error) Obit_traceback_msg (err, routine, in->outOTF->name);

  /* Open table */
  if ((ObitTableOTFIndexOpen (table, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFIndex table");
    return;
  }

  /* Create Row */
  row = newObitTableOTFIndexRow (table);

  /* initialize row */
  row->ScanID = scanID;
  row->TargetID = targetID;
  row->Time = 0.0;
  row->TimeI = 0.0;
  row->StartRec = -1;
  row->EndRec = -1;

  /* attach to table buffer */
  ObitTableOTFIndexSetRow (table, row, err);
  if (err->error) Obit_traceback_msg (err, routine, in->outOTF->name);

  /* Write at end of table */
  iRow = table->myDesc->nrow + 1;
  if ((ObitTableOTFIndexWriteRow (table, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR writing OTFIndex Table file");
    return;
  }

 /* Close  table */
  if ((ObitTableOTFIndexClose (table, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFIndex Table file");
    return;
  }

  /* Cleanup */
  row = ObitTableOTFIndexRowUnref(row);
  table = ObitTableOTFIndexUnref(table);

} /* end  InitScan */

/**
 * Initializes Index table for this scan creating an entry 
 * \param  in       The conversion object
 * \param  outData   Output OTF object 
 * \param  isNew     True if output file just created   
 * \param  startTime Start time of scan in days
 * \param  endTime   End time of scan in days
 * \param  startRec  First record in scan
 * \param  endRec    Last record in scan
 * \param  err       Obit return error stac
 */
static void SetScan (ObitGBTDCROTF *in, ObitErr *err)
{
  ObitTableOTFIndex* table;
  ObitTableOTFIndexRow* row;
  olong iRow, ver;
  ObitIOAccess access;
  gchar *routine = "SetScan";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* create Index table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  table = newObitTableOTFIndexValue ("Index table", (ObitData*)in->outOTF, &ver, access, err);
  if (err->error) Obit_traceback_msg (err, routine, in->outOTF->name);

  /* Open table */
  if ((ObitTableOTFIndexOpen (table, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFIndex table");
    return;
  }

  /* Create Row */
  row = newObitTableOTFIndexRow (table);

  /* attach to table buffer */
  ObitTableOTFIndexSetRow (table, row, err);
  if (err->error) Obit_traceback_msg (err, routine, in->outOTF->name);

  /* Update last record */
  iRow = table->myDesc->nrow;
  if ((ObitTableOTFIndexReadRow (table, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR reading OTFIndex Table file");
    return;
  }

  /* upate row */
  row->Time = 0.5 * (in->startTime + in->endTime);
  row->TimeI = (in->endTime - in->startTime);
  row->StartRec = in->startRec;
  row->EndRec   = in->endRec;

  /* Rewrite at end of table */
  iRow = table->myDesc->nrow;
  if ((ObitTableOTFIndexWriteRow (table, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR writing OTFIndex Table file");
    return;
  }

 /* Close  table */
  if ((ObitTableOTFIndexClose (table, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFIndex Table file");
    return;
  }

  /* Cleanup */
  row = ObitTableOTFIndexRowUnref(row);
  table = ObitTableOTFIndexUnref(table);

} /* end  SetScan */

