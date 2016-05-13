/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010-2016                                          */
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

#include "ObitBDFData.h"
#include "ObitFile.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitBDFData.c
 * ObitBDFData class function definitions.
 * This class is derived from the Obit base class.
 * This class accesses data in the EVLA BDF format
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitBDFData";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitBDFDataClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitBDFDataClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitBDFDataInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitBDFDataClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitBDFDataClassInfoDefFn (gpointer inClass);

/* BDF XML routines */
/** Private: Parse integer from XML string  */
static olong BDFparse_int(gchar *string, olong maxChar, 
			   gchar *prior, gchar **next);
/** Private: Parse string from XML string  */
static gchar* BDFparse_str(gchar *string, olong maxChar, 
			   gchar *prior, gchar **next);
/** Private: Parse quotedstring from XML string  */
static gchar* BDFparse_quote_str(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next);
/** Private: Parse time from XML string  */
static odouble BDFparse_time(gchar *string, olong maxChar, 
			     gchar *prior, gchar **next);

/** Private: Parse time interval from XML string  */
static odouble BDFparse_timeint(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next);

/** Private: Parse string array from XML string  */
static gchar** BDFparse_strarray(gchar *string, olong maxChar, 
				  gchar *prior, gchar **next);
/** Private: Parse axis order array from XML string  */
static ObitBDFAxisName* BDFparse_axesarray(gchar *string, olong maxChar, 
				       gchar *prior, gchar **next);

/** Private: Look up axis name enum */
static ObitBDFAxisName LookupAxisName(gchar *name);

/** Private: Look up data type enum */
static ObitBDFDataType LookupDataType(gchar *name);

/** Private: Look up endian enum */
static ObitBDFEndian LookupEndian(gchar *name);

/** Private: delete BDFSpecWindowInfo */
static BDFSpecWindowInfo* KillBDFSpecWindowInfo(BDFSpecWindowInfo *info);

/** Private: delete BDFBasebandInfo */
static BDFBasebandInfo* KillBDFBasebandInfo(BDFBasebandInfo *info);

/** Private: delete BDFScanInfo */
static BDFScanInfo* KillBDFScanInfo(BDFScanInfo *info);

/** Private: delete BDFIntegInfo */
static BDFIntegInfo* KillBDFIntegInfo(BDFIntegInfo *info);

/* Get start of next MIME segment */
static ObitBDFMIMEType GetNextMIME(ObitBDFData *in, 
				   gchar *last, gchar **start, 
				   ObitErr *err);

/* Get start of next MIME segment */
static ObitBDFMIMEType GetNextMIMEType(ObitBDFData *in, 
				       gchar *last, gchar **start, 
				       ObitErr *err);
/* Copy binary float data */
static ObitIOCode CopyFloats (ObitBDFData *in, 
			      gchar *start, ofloat *target, olong n, 
			      gboolean byteFlip, ofloat scale, 
			      ObitBDFDataType Dtype, ObitErr *err);

/* Copy binary long integer data */
static ObitIOCode CopyLongs (ObitBDFData *in, 
			     gchar *start, olong *target, olong n, 
			     gboolean byteFlip, ObitBDFDataType Dtype, 
			     ObitErr *err);

/** Private: Squeeze all blanks out of a string */
static void Strip (gchar *s);

/** Private: Get axis order */
static olong GetAxisOrder (ObitBDFAxisName *axes, ObitBDFAxisName axis);

/** Private: find string */
static gchar* findString (gchar *start, olong maxStr, gchar *match);
/*----------------- Union definitions ----------------------*/
/** Used for byte swapping shorts */
 union sequiv { 
   gshort full; 
   gchar parts[2];
 }; 

/** Used for byte swapping floats */
 union fequiv { 
   ofloat full; 
   gchar parts[4]; 
 }; 

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitBDFData* newObitBDFData (gchar* name)
{
  ObitBDFData* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitBDFDataClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitBDFData));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitBDFDataInit((gpointer)out);

 return out;
} /* end newObitBDFData */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitBDFDataGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitBDFDataClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitBDFDataGetClass */

/**
 * Make a deep copy of an ObitBDFData. NYI
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitBDFData* ObitBDFDataCopy  (ObitBDFData *in, ObitBDFData *out, ObitErr *err)
{
  /*const ObitClassInfo *ParentClass;*/
  /*gboolean oldExist;*/
  /*gchar *outName;*/

  /* error checks */
  if (err->error) return out;

  /* Stubbed */
  g_error("ObitBDFDataCopy: Stubbed");

  return out;
} /* end ObitBDFDataCopy */

 /**
 * Make a copy of a object but do not copy the actual data NYI
 * This is useful to create an BDFData similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitBDFDataClone  (ObitBDFData *in, ObitBDFData *out, ObitErr *err)
{
  /*const ObitClassInfo *ParentClass;*/

  /* error checks */
  if (err->error) return;

  /* Stubbed */
  g_error("ObitBDFDataClone: Stubbed");

} /* end ObitBDFDataClone */

/**
 * Creates an ObitBDFData 
 * Parses the ASMD XML tables and stores
 * \param name      An optional name for the object.
 * \param desc      UV descriptor for data extracted
 * \param  SDMData  ASDM structure    
 * \param err       Obit error stack object.
 * \return the new object.
 */
ObitBDFData* ObitBDFDataCreate (gchar* name,  ObitUVDesc *desc, 
				ObitSDMData *SDMData,
				ObitErr *err)
{
  ObitBDFData* out;
  /*gchar *routine="ObitBDFDataCreate";*/

  /* Create basic structure */
  out = newObitBDFData (name);

  /* Save descriptor */
  out->desc = ObitUVDescRef(desc);

  /* Save the ASDM */
  out->SDMData = ObitSDMDataRef(SDMData);

  /* Init buffer */
  out->buffer = g_malloc0(BDFBUFFERSIZE*BDFBUFFERFRAMES+100);
  out->nBytesInBuffer = 0;
  out->current = out->buffer;

  return out;
} /* end ObitBDFDataCreate */

 /**
 * Initialize File
 * Initializes buffer, parses scan XML header
 * \param in       The object to fill
 * \param DataFile Name of Mime file with data
 * \param err      Obit error stack object.
 * \param err Obit error stack object.
 */
void ObitBDFDataInitFile  (ObitBDFData *in, gchar *DataFile, ObitErr *err)
{
  ObitIOCode retCode;
  gchar *routine = "ObitBDFDataInitFile";

  /* error checks */
  if (err->error) return;

   /* set Values - file name */
  in->DataFile = strdup(DataFile);

  /* Get size */
  in->fileSize = ObitFileSize(DataFile, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create file object */
  if (in->file==NULL) in->file = newObitFile ("BDF File");

  /* Open */
  ObitFileOpen (in->file, in->DataFile, OBIT_IO_ReadOnly, OBIT_IO_Binary, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Fill Buffer */
  retCode = ObitBDFDataFillBuffer (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitBDFDataInitFile */

 /**
 * Fill Buffer
 * If the buffer is filled to capacity, the bottom frame is copied to
 * to the top of the buffer and the remainder of the buffer filled.
 * Updates in->nBytesInBuffer, in->current.
 * \param in  The object to fill
 * \param err Obit error stack object.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_EOF = EOF.
 */
ObitIOCode ObitBDFDataFillBuffer (ObitBDFData *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  olong i, ncopy, size; 
  gchar *top, *copy;
  gchar *routine = "ObitBDFDataFillBuffer";

  /* error checks */
  if (err->error) return retCode;

  if ((in->fileSize-in->file->filePos) <= 0) return OBIT_IO_EOF;

   /* Is it already full? */
  if (in->nBytesInBuffer>=BDFBUFFERSIZE*BDFBUFFERFRAMES) { /* Yes - shift */
    top  = in->buffer;
    copy = in->buffer + (BDFBUFFERFRAMES-1) * BDFBUFFERSIZE;
    memmove (top, copy, (size_t)BDFBUFFERSIZE);
    top = &in->buffer[BDFBUFFERSIZE];
    ncopy = (BDFBUFFERFRAMES-1);
    in->nBytesInBuffer = BDFBUFFERSIZE;
    in->current -= (BDFBUFFERFRAMES-1) * BDFBUFFERSIZE; /* Current position */
  } else {  /* Nope - fill 'er up */
    top  = in->buffer;
    ncopy = BDFBUFFERFRAMES;
    in->nBytesInBuffer = 0;
    in->current = in->buffer; /* Current position */
  }

  /* Do reads */
  for (i=0; i<ncopy; i++) {
    /* No more than what's left */
    size = MIN ((olong)BDFBUFFERSIZE, (in->fileSize-in->file->filePos));
    if (size<=0) break;
    retCode = ObitFileRead (in->file, -1, size, top, err);
    /* DEBUG
    fprintf (stderr, "Read size %d pos %lld filesize %lld\n", size, in->file->filePos, in->fileSize); */
    if (err->error) {
      Obit_traceback_val (err, routine, in->name, retCode);
    }
    in->nBytesInBuffer += size;
    if (in->file->filePos>=in->fileSize) break;
    top += BDFBUFFERSIZE;
  }
  return retCode;
} /* end ObitBDFDataFillBuffer */

 /**
 * Parse scan info from buffer
 * Expects all of scan XML in buffer
 * \param in      The object to update
 * \param iMain   The ASDM Main table row
 * \param SWOrder If TRUE leave data is SW order
 * \param selChan Number of channels selected
 * \param selIF   Number of IFs selected
 * \param err     Obit error stack object.
 */
void ObitBDFDataInitScan  (ObitBDFData *in, olong iMain, gboolean SWOrder,
			   olong selChan, olong selIF, ObitErr *err)
{
  gchar *startInfo, *endInfo, *startBB, *endBB, *prior, *next, *xnext, *tstr;
  olong  configDescriptionId, fieldId, inext, ScanId=0, iScan, iIntent, iField;
  olong maxStr, maxStr2, i, j, count, *antIds, iConfig, iAnt, jAnt, iSW, jSW=0, jSource;
  olong blOrder, polnOrder, freqOrder, SPWOrder, BBOrder, APCOrder, binOrder;
  olong *SWoff=NULL, flagoff;
  gboolean done;
  gchar *aname;
  ObitIOCode retCode = OBIT_IO_OK;
  gchar *routine = "ObitBDFDataInitScan";

  /* error checks */
  if (err->error) return;

  /* Create info structure if not there */
  if (in->ScanInfo==NULL) in->ScanInfo = g_malloc0(sizeof(BDFScanInfo));

  /* Is this EVLA data? */
  in->isEVLA = !strncmp(in->SDMData->ExecBlockTab->rows[0]->telescopeName, "EVLA", 4);
  /* Is this ALMA data? */
  in->isALMA = !strncmp(in->SDMData->ExecBlockTab->rows[0]->telescopeName, "ALMA", 4);
  /* ALMA is a bit confused */
  if (!in->isALMA) in->isALMA = !strncmp(in->SDMData->ExecBlockTab->rows[0]->telescopeName, "OSF", 3);
  if (!in->isALMA) in->isALMA = !strncmp(in->SDMData->ExecBlockTab->rows[0]->telescopeName, "AOS", 3);

  /* Init */
  in->ScanInfo->iMain         = iMain;
  in->ScanInfo->scanNumber    = in->SDMData->MainTab->rows[iMain]->scanNumber;
  in->ScanInfo->subScanNumber = in->SDMData->MainTab->rows[iMain]->subscanNumber;
  in->ScanInfo->numAntenna  = -1;
  in->ScanInfo->numBaseband = 0;
  in->haveFlag = FALSE;
  if (in->ScanInfo->FlagAxes) 
    {g_free(in->ScanInfo->FlagAxes); in->ScanInfo->FlagAxes = NULL;}
  in->haveActualTimes = FALSE;
  if (in->ScanInfo->actualTimesAxes) 
    {g_free(in->ScanInfo->actualTimesAxes); in->ScanInfo->actualTimesAxes = NULL;}
  in->haveActualDurations = FALSE;
  if (in->ScanInfo->actualDurationsAxes) 
    {g_free(in->ScanInfo->actualDurationsAxes); in->ScanInfo->actualDurationsAxes = NULL;}
  in->haveCrossData = FALSE;
  if (in->ScanInfo->crossDataAxes) 
    {g_free(in->ScanInfo->crossDataAxes); in->ScanInfo->crossDataAxes = NULL;}
  in->haveAutoData = FALSE;
  if (in->ScanInfo->autoDataAxes) 
    {g_free(in->ScanInfo->autoDataAxes); in->ScanInfo->autoDataAxes = NULL;}
  in->haveWeight = FALSE;
  if (in->ScanInfo->weightAxes) 
    {g_free(in->ScanInfo->weightAxes); in->ScanInfo->weightAxes = NULL;}
  in->haveZeroLag = FALSE;
  if (in->ScanInfo->zeroLagAxes) 
    {g_free(in->ScanInfo->zeroLagAxes); in->ScanInfo->zeroLagAxes = NULL;}

  /* Lookup scan ID */
  for (iScan=0; iScan<in->SDMData->ScanTab->nrows; iScan++) {
    ScanId = iScan;
    if ( in->SDMData->MainTab->rows[iMain]->scanNumber == 
 	 in->SDMData->ScanTab->rows[ScanId]->scanNumber) break;
  }
  in->ScanInfo->scanId = ScanId;
  
  /* Is this a holography scan/subscan? */
  in->ScanInfo->isHolo = FALSE;
  for (iIntent=0; iIntent<in->SDMData->ScanTab->rows[ScanId]->numIntent; iIntent++) {
    if (in->SDMData->ScanTab->rows[ScanId]->scanIntent[iIntent]==NULL) break;
    in->ScanInfo->isHolo = in->ScanInfo->isHolo ||  
      (!strncmp(in->SDMData->ScanTab->rows[ScanId]->scanIntent[iIntent], 
 		"MAP_ANTENNA_SURFACE", 19));
    if (in->ScanInfo->isHolo) break;
  }
  
  /* If isHolo check subscan - may be phase reference = not holography */
  if (in->ScanInfo->isHolo) {
    for (iScan=0; iScan<in->SDMData->SubscanTab->nrows; iScan++) {
      ScanId = iScan;
      if ((in->ScanInfo->scanNumber==in->SDMData->SubscanTab->rows[ScanId]->scanNumber) &&
 	  (in->ScanInfo->subScanNumber==in->SDMData->SubscanTab->rows[ScanId]->subscanNumber))
 	break;
    }
    if ((in->SDMData->SubscanTab->rows[ScanId]->subscanIntent!=NULL) && 
	!strncmp(in->SDMData->SubscanTab->rows[ScanId]->subscanIntent,
 		 "REFERENCE", 9)) in->ScanInfo->isHolo = FALSE;
  }
  
  /* Parse scan header - first find limits */
  /* Find first non null */
  while (*in->current==0) (*in->current)++;  /* dangerous */
  maxStr    = in->nBytesInBuffer - (in->current-in->buffer);
  /*startInfo = g_strstr_len (in->current, maxStr, "<sdmDataHeader ");*/
  startInfo = findString (in->current, maxStr, "<sdmDataHeader ");
  /*endInfo   = g_strstr_len (in->current, maxStr, "</sdmDataHeader>");*/
  endInfo   = findString (in->current, maxStr, "</sdmDataHeader>");
  maxStr    = (olong)(endInfo-startInfo);
  
  /* Start Time */
  prior = "<startTime>";
  in->ScanInfo->startTime = BDFparse_time(startInfo, maxStr, prior, &next);
  
  /* Number of antennas */
  prior = "<numAntenna>";
  in->ScanInfo->numAntenna  = BDFparse_int(startInfo, maxStr, prior, &next);
  
  /* Number of times */
  prior = "<dimensionality axes=\"TIM\">";
  in->ScanInfo->numTime  = BDFparse_int(startInfo, maxStr, prior, &next);
  /* ALMA variant */
  prior = "<numTime>";
  in->ScanInfo->numTime  = BDFparse_int(startInfo, maxStr, prior, &next);
  in->numTime = in->ScanInfo->numTime;
  in->currTimeIndx = 0;
  
  /* correlation Mode */
  prior = "<correlationMode>";
  tstr = BDFparse_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    Strip(tstr);
    if (!strcmp(tstr, "CROSS_ONLY"))          
      in->ScanInfo->correlationMode = BDFCorrMode_CROSS_ONLY;
    else if (!strcmp(tstr, "AUTO_ONLY"))      
      in->ScanInfo->correlationMode = BDFCorrMode_AUTO_ONLY;
    else if (!strcmp(tstr, "CROSS_AND_AUTO")) 
      in->ScanInfo->correlationMode = BDFCorrMode_CROSS_AND_AUTO;
    g_free(tstr);
  }
  
  /* spectral Resolution */
  prior = "<spectralResolution>";
  tstr = BDFparse_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    Strip(tstr);
    if (!strcmp(tstr, "CHANNEL_AVERAGE"))      
      in->ScanInfo->spectralResolution = BDFSpecRes_CHANNEL_AVERAGE;
    else if (!strcmp(tstr, "BASEBAND_WIDE"))   
      in->ScanInfo->spectralResolution = BDFSpecRes_BASEBAND_WIDE;
    else if (!strcmp(tstr, "FULL_RESOLUTION")) 
      in->ScanInfo->spectralResolution = BDFSpecRes_FULL_RESOLUTION;
    g_free(tstr);
  }
  
  /* Endian */
  prior = "byteOrder=";
  tstr = BDFparse_quote_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    in->ScanInfo->endian = LookupEndian(tstr);
    g_free(tstr);
  }
  
  /* first find limits of dataStruct */
  maxStr    = in->nBytesInBuffer - (in->current-in->buffer);
  /*startInfo = g_strstr_len (in->current, maxStr, "<dataStruct ");*/
  startInfo = findString (in->current, maxStr, "<dataStruct ");
  /* May need next buffer */
  while (startInfo==NULL) {
    retCode = ObitBDFDataFillBuffer (in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    maxStr = in->nBytesInBuffer;
    /*startInfo = g_strstr_len (in->buffer, maxStr, "<sdmDataSubsetHeader ");*/
    startInfo = findString (in->buffer, maxStr, "<sdmDataSubsetHeader ");
   if (retCode==OBIT_IO_EOF) startInfo = in->buffer;
  }
  /* If the entire file was read and no start info - then something is wrong */
  if (retCode==OBIT_IO_EOF) return;
  maxStr    = in->nBytesInBuffer - (olong)(startInfo-in->buffer);
  /*endInfo   = g_strstr_len (startInfo, maxStr, "</dataStruct>");*/
  endInfo   = findString (startInfo, maxStr, "</dataStruct>");
  maxStr    = (olong)(endInfo-startInfo);

  /* flags */
   prior = "<flags size=";
  tstr = BDFparse_quote_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    in->ScanInfo->FlagSize = (olong)strtol(tstr, &xnext, 10);
    /* Axes types */
    prior = "axes=";
    next++;
    in->ScanInfo->FlagAxes = BDFparse_axesarray (next, maxStr, prior, &next);
    g_free(tstr);
  }
  
  /* actualTimes  */
  prior = "<actualTimes size=";
  tstr = BDFparse_quote_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    in->ScanInfo->actualTimesSize = (olong)strtol(tstr, &xnext, 10);
    /* Axes types */
    prior = "axes=";
    next++;
    in->ScanInfo->actualTimesAxes = BDFparse_axesarray (next, maxStr, prior, &next);
    g_free(tstr);
  }
  
  /* actualDurations  */
  prior = "<actualDurations size=";
  tstr = BDFparse_quote_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    in->ScanInfo->actualDurationsSize = (olong)strtol(tstr, &xnext, 10);
    /* Axes types */
    prior = "axes=";
    in->ScanInfo->actualDurationsAxes = BDFparse_axesarray (next, maxStr, prior, &next);
    g_free(tstr);
  }
  
  /*  crossData */
  prior = "<crossData size=";
  tstr = BDFparse_quote_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    in->ScanInfo->crossDataSize = (olong)strtol(tstr, &xnext, 10);
    /* Axes types */
    prior = "axes=";
    next++;
    in->ScanInfo->crossDataAxes = BDFparse_axesarray (next, maxStr, prior, &next);
    g_free(tstr);
  }
  
  /* autoData  */
  prior = "<autoData size=";
  tstr = BDFparse_quote_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    in->ScanInfo->autoDataSize = (olong)strtol(tstr, &xnext, 10);
    /* Axes types */
    prior = "axes=";
    next++;
    in->ScanInfo->autoDataAxes = BDFparse_axesarray (next, maxStr, prior, &next);
    g_free(tstr);
  }
  
  /* weight  */
  prior = "<weight size=";
  tstr = BDFparse_quote_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    in->ScanInfo->weightSize = (olong)strtol(tstr, &xnext, 10);
    /* Axes types */
    prior = "axes=";
    next++;
    in->ScanInfo->weightAxes = BDFparse_axesarray (next, maxStr, prior, &next);
    g_free(tstr);
  }
  
  /* zero lags  */
  prior = "<zeroLags size=";
  tstr = BDFparse_quote_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    in->ScanInfo->zeroLagSize = (olong)strtol(tstr, &xnext, 10);
    /* Axes types */
    prior = "axes=";
    next++;
    in->ScanInfo->zeroLagAxes = BDFparse_axesarray (next, maxStr, prior, &next);
    g_free(tstr);
  }
  
  /* Parse basebands */
  i = 0;
  in->ScanInfo->numBaseband = 0;
  in->ScanInfo->numBin = 0;
  next = startInfo;
  while (i<MAXBBINFO) {
    /* first find limits of next baseband info */
    maxStr = (olong)(endInfo-next);
    /*startBB = g_strstr_len (next, maxStr, "<baseband ");*/
    startBB = findString (next, maxStr, "<baseband ");
    /* More? */
    if (startBB==NULL) break;
    maxStr2 = (olong)(endInfo-next);
    endBB   = findString (startBB, maxStr, "</baseband>");
    maxStr2 = (olong)(endBB-startBB);
 
    in->ScanInfo->BBinfo[i] = KillBDFBasebandInfo(in->ScanInfo->BBinfo[i]);/* delete old */
    in->ScanInfo->BBinfo[i] = g_malloc0(sizeof(BDFBasebandInfo));          /* Allocate new */

    /* Allocate MAXBBSW Spectral window array */
    in->ScanInfo->BBinfo[i]->SWinds = g_malloc0(MAXBBSW*sizeof(BDFSpecWindowInfo*)); 
    in->ScanInfo->BBinfo[i]->numSpectralWindow = 0;

    /* Parse */
    /* basebandName */
    prior = "<baseband name=";
    in->ScanInfo->BBinfo[i]->basebandName = BDFparse_quote_str(startBB, maxStr2, prior, &next);
  
    /* loop over spectral window numbers */
    done = FALSE;
    j = 0;
    while (!done) {
      prior = "<spectralWindow sw=";
      tstr = BDFparse_quote_str (startBB, maxStr2, prior, &next);
      if (tstr==NULL) break;  /* All? */
      in->ScanInfo->BBinfo[i]->SWinds[j] = g_malloc0(sizeof(BDFSpecWindowInfo)); 
      if (tstr!=NULL) {
	in->ScanInfo->BBinfo[i]->SWinds[j]->spectralWindowNum = (olong)strtol(tstr, &next, 10);
	g_free (tstr);
      } else in->ScanInfo->BBinfo[i]->SWinds[j]->spectralWindowNum = 1;
      
      /* List of single dish (autocorrelation) products  */
      prior = "sdPolProducts=\"";
      in->ScanInfo->BBinfo[i]->SWinds[j]->sdPolProducts = BDFparse_strarray(startBB, maxStr2, prior, &next);
      count = 0;
      if (in->ScanInfo->BBinfo[i]->SWinds[j]->sdPolProducts) {
	while (in->ScanInfo->BBinfo[i]->SWinds[j]->sdPolProducts[count]) {count++;}
      }
      in->ScanInfo->BBinfo[i]->SWinds[j]->numSdPolProducts = count;
      
      /* List of crosscorrelation products */
      prior = "crossPolProducts=\"";
      in->ScanInfo->BBinfo[i]->SWinds[j]->crossPolProducts = BDFparse_strarray(startBB, maxStr2, prior, &next);
      /* Count 'em */
      count = 0;
      if (in->ScanInfo->BBinfo[i]->SWinds[j]->crossPolProducts) {
	while (in->ScanInfo->BBinfo[i]->SWinds[j]->crossPolProducts[count]) {count++;}
      }
      in->ScanInfo->BBinfo[i]->SWinds[j]->numCrossPolProducts = count;
      
      /* Number of spectral points  */
      prior = "numSpectralPoint=";
      tstr = BDFparse_quote_str (startBB, maxStr2, prior, &next);
      if (tstr!=NULL) {
	in->ScanInfo->BBinfo[i]->SWinds[j]->numSpectralPoint = (olong)strtol(tstr, &next, 10);
	g_free (tstr);
      } else in->ScanInfo->BBinfo[i]->SWinds[j]->numSpectralPoint = 1;
      
      /* Number of data (e.g. pulsar) bins  */
      prior = "numBin=";
      tstr = BDFparse_quote_str (startBB, maxStr2, prior, &next);
      if (tstr!=NULL) {
	in->ScanInfo->BBinfo[i]->SWinds[j]->numBin = (olong)strtol(tstr, &next, 10);
	g_free (tstr);      
      } else in->ScanInfo->BBinfo[i]->SWinds[j]->numBin = 1;
      in->ScanInfo->numBin = MAX (in->ScanInfo->numBin, in->ScanInfo->BBinfo[i]->SWinds[j]->numBin);
      
      /* Data scaling factor  */
      prior = "scaleFactor=";
      tstr = BDFparse_quote_str (startBB, maxStr2, prior, &next);
      if (tstr!=NULL) {
	in->ScanInfo->BBinfo[i]->SWinds[j]->scaleFactor = (odouble)strtod(tstr, &next);
	g_free (tstr);
      } else in->ScanInfo->BBinfo[i]->SWinds[j]->scaleFactor = 1.0;
     
      /* Sideband */
      prior = "sideband=";
      in->ScanInfo->BBinfo[i]->SWinds[j]->sideband = BDFparse_quote_str(startBB, maxStr2, prior, &next);

      /* Update text search range */
      maxStr2 -= next-startBB;
      startBB = next;
      
      in->ScanInfo->BBinfo[i]->numSpectralWindow++;
      j++;  /* Spectral window index */
    } /* End loop over baseband spectral windows */
    in->ScanInfo->numBaseband++;  /* Count basebands */
    i++;  /* baseband index */
 } /* end loop over baseband */
  
  in->current = next;  /* where in buffer */

  /* Only some orders of data axes are supported - check crossData */
  if (in->ScanInfo->crossDataAxes) {
    polnOrder = GetAxisOrder (in->ScanInfo->crossDataAxes, BDFAxisName_POL);
    if (polnOrder<0)
      polnOrder = GetAxisOrder (in->ScanInfo->crossDataAxes, BDFAxisName_STO);
    freqOrder =  GetAxisOrder (in->ScanInfo->crossDataAxes,  BDFAxisName_SPP);
    SPWOrder  =  GetAxisOrder (in->ScanInfo->crossDataAxes,  BDFAxisName_SPW);
    BBOrder   =  GetAxisOrder (in->ScanInfo->crossDataAxes,  BDFAxisName_BAB);
    APCOrder  =  GetAxisOrder (in->ScanInfo->crossDataAxes,  BDFAxisName_APC);
    blOrder   =  GetAxisOrder (in->ScanInfo->crossDataAxes,  BDFAxisName_BAL);
    binOrder  =  GetAxisOrder (in->ScanInfo->crossDataAxes,  BDFAxisName_BIN);

    /* Trap ALMA (WVR?) pathology */
    if ((polnOrder<0) && (freqOrder>=0) && (SPWOrder<0)) polnOrder = freqOrder+1;
   
    Obit_return_if_fail((blOrder==0), err,
			"%s: Only support baseline as first axis",  routine);
    Obit_return_if_fail(((binOrder<0) || (in->ScanInfo->numBin<=1)), err,
			"%s: Binned data not supported", routine);
    Obit_return_if_fail((freqOrder<polnOrder), err,
			"%s: Unsupported Freq/poln axis order", routine);
    Obit_return_if_fail(((APCOrder<0) || ((APCOrder<freqOrder) && (APCOrder<polnOrder))), err,
			"%s: Unsupported APC axis order", routine);
    Obit_return_if_fail(((BBOrder<0) || ((BBOrder<freqOrder) && (BBOrder<polnOrder))), err,
			"%s: Unsupported Baseband axis order", routine);
    Obit_return_if_fail(((SPWOrder<0) || ((SPWOrder<freqOrder) && (SPWOrder<polnOrder))), err,
			"%s: Unsupported Baseband axis order", routine);
  }


  /* Numbers of things */
  in->numBaseline       = (in->ScanInfo->numAntenna * (in->ScanInfo->numAntenna-1))/2;
  in->currBaseline      = 0;
  in->numBaseband       = in->ScanInfo->numBaseband;

  /* Spectral window array */
  in->SWArray = ObitSDMDataKillSWArray(in->SWArray);
  in->SWArray = ObitSDMDataGetSWArray (in->SDMData, iMain, SWOrder);

  /* Set selection */
  ObitBDFDataSelChan (in, selChan, selIF,  ASDMBand_Any);

  /* More reliable - find selected SW, count */
  in->numSpectralWindow = 0;
  for (iSW=0; iSW<in->SWArray->nwinds; iSW++) {
    if (in->SWArray->winds[iSW]->selected) {
      jSW = iSW;
      in->numSpectralWindow++;
    }
  }
  in->numSpectralChann  = in->SWArray->winds[jSW]->numChan;
  in->numCPoln          = in->SWArray->winds[jSW]->nCPoln;
  in->numAPoln          = in->SWArray->winds[jSW]->nAPoln;

  /* Init antennas */
  in->nant   = in->ScanInfo->numAntenna;
  if (in->antNo) g_free(in->antNo);
  in->antNo  = g_malloc0(in->nant*sizeof(olong));
  if (in->antId) g_free(in->antId);
  in->antId  = g_malloc0(200*sizeof(olong));

  /* Lookup table of antenna numbers corresponding to IDs */

  /* Find entry  in configDescription table */
  configDescriptionId = 
    in->SDMData->MainTab->rows[in->ScanInfo->iMain]->configDescriptionId;
  for (iConfig=0; iConfig<in->SDMData->ConfigDescriptionTab->nrows; iConfig++) {
    if (in->SDMData->ConfigDescriptionTab->rows[iConfig]->configDescriptionId==configDescriptionId) break;
  }
  Obit_return_if_fail((iConfig<in->SDMData->ConfigDescriptionTab->nrows), err,
		      "%s: Could not find configDescriptionId %d in ASDM", 
		      routine, configDescriptionId);

  /* Antenna id array */
  antIds = in->SDMData->ConfigDescriptionTab->rows[iConfig]->antennaId;
  /* Loop over antennas */
  for (iAnt=0; iAnt<in->nant; iAnt++) {
    /* See if past all on list */
    if (antIds[iAnt]<0) continue;
    /* Find antenna */
    for (jAnt=0; jAnt<in->SDMData->AntennaTab->nrows; jAnt++) {
      if (in->SDMData->AntennaTab->rows[jAnt]->antennaId==antIds[iAnt]) break;
    }
    Obit_return_if_fail((jAnt<in->SDMData->AntennaTab->nrows), err,
			"%s: Could not find antennaId %d in ASDM", 
			routine, antIds[iAnt]);

    /* Crack name to get number Assume EVLA starts with "ea" */
    aname = in->SDMData->AntennaTab->rows[jAnt]->name;
    if ((aname[0]=='e') && (aname[1]=='a'))
      in->antNo[iAnt]    = (olong)strtol(&aname[2],NULL,10);
    else in->antNo[iAnt] = in->SDMData->AntennaTab->rows[jAnt]->antennaNo;
    /* Inverse lookup */
    in->antId[in->antNo[iAnt]] = iAnt;
  } /* end loop over antennas */

  /* Atmospheric phase correction stuff */
  in->numAtmCorr = in->SDMData->ConfigDescriptionTab->rows[iConfig]->numAtmPhaseCorrection;
  in->numAtmCorr = MIN (2, MAX(1, in->numAtmCorr));
  if (in->numAtmCorr>1) {
    if ((in->selAtmCorr==FALSE) &&
	(in->SDMData->ConfigDescriptionTab->rows[iConfig]->atmPhaseCorrection[0]==ASDMAtmPhCorr_AP_UNCORRECTED))
      in->offAtmCorr = 0;
    else if ((in->selAtmCorr==TRUE) &&
	     (in->SDMData->ConfigDescriptionTab->rows[iConfig]->atmPhaseCorrection[0]==ASDMAtmPhCorr_AP_CORRECTED))
      in->offAtmCorr = 0;
    else if ((in->selAtmCorr==FALSE) &&
	     (in->SDMData->ConfigDescriptionTab->rows[iConfig]->atmPhaseCorrection[1]==ASDMAtmPhCorr_AP_UNCORRECTED))
      in->offAtmCorr = 1;
    else if ((in->selAtmCorr==TRUE) &&
	(in->SDMData->ConfigDescriptionTab->rows[iConfig]->atmPhaseCorrection[1]==ASDMAtmPhCorr_AP_CORRECTED))
      in->offAtmCorr = 1;
    else  in->offAtmCorr = 0; /* shouldn't happen */
  } else {
    in->offAtmCorr = 0;
  }

  /* Source No. - find in in->SDMData->SourceArray with key fieldId */
  fieldId = in->SDMData->MainTab->rows[in->ScanInfo->iMain]->fieldId;
  /* Depending of whether doQual (1 source per subscan) */
  if (in->SDMData->doQual) {
    for (jSource=0; jSource<in->SDMData->SourceArray->nsou; jSource++) {
      if (in->SDMData->SourceArray->sou[jSource]->fieldId==fieldId) break;
    }
    Obit_return_if_fail((jSource<in->SDMData->SourceArray->nsou), err,
			"%s: Could not find source Id %d in ASDM", 
			routine, fieldId);
  } else {  /* One source per scan - use Field Table */
    for (iField=0; iField<in->SDMData->FieldTab->nrows; iField++) {
      if (in->SDMData->FieldTab->rows[iField]->fieldId==fieldId) break;
    }
    Obit_return_if_fail((iField<in->SDMData->FieldTab->nrows), err,
			"%s: Could not find field Id %d in ASDM", 
			routine, fieldId);
    /* Now look in SourceArray for name */
    for (jSource=0; jSource<in->SDMData->SourceArray->nsou; jSource++) {
      if (!strcmp(in->SDMData->SourceArray->sou[jSource]->sourceName, 
		  in->SDMData->FieldTab->rows[iField]->fieldName)) break;
    }
    Obit_return_if_fail((jSource<in->SDMData->SourceArray->nsou), err,
			"%s: Could not find source %s in ASDM", 
			routine, in->SDMData->SourceTab->rows[iField]->sourceName);
     
  }
  in->sourceNo = in->SDMData->SourceArray->sou[jSource]->sourceNo;
  if (in->curSource) g_free(in->curSource);
  in->curSource = strdup(in->SDMData->SourceArray->sou[jSource]->sourceName);
  in->sourceQual = in->SDMData->SourceArray->sou[jSource]->sourceQual;

  /* Cross correlation frequency increment = no poln. x 2 */
  in->cincf  = in->numCPoln * 2;
  in->cfincf = in->numCPoln;   /* Flag frequency increment */

  /* Cross correlation IF increment = no poln.x no Chan x 2*/
  in->cincif  = in->numCPoln * in->numSpectralChann * 2;
  in->cfincif = in->numCPoln * in->numSpectralChann;   /* Flag if increment */

  /* Ordering of cross correlation polarizations - shuffle order */
  if (in->coffs) g_free(in->coffs);
  in->coffs = g_malloc0(in->numCPoln*sizeof(olong));
  if (in->cfoffs) g_free(in->cfoffs);
  in->cfoffs = g_malloc0(in->numCPoln*sizeof(olong));
  if (in->numCPoln==1) {
    in->coffs[0]  = 0;   /* Vis */
    in->cfoffs[0] = 0;   /* Flag */
  } else if (in->numCPoln==2) {
    in->coffs[0]  = 0;   /* Vis */
    in->coffs[1]  = 2;   /* Vis */
    in->cfoffs[0] = 0;   /* Flag */
    in->cfoffs[1] = 1;   /* Flag */
 } else if (in->numCPoln==4) {
    in->coffs[0]  = 0;   /* Vis */
    in->coffs[1]  = 6;   /* Vis */
    in->coffs[2]  = 2;   /* Vis */
    in->coffs[3]  = 4;   /* Vis */
    in->cfoffs[0] = 0;   /* Flag */
    in->cfoffs[1] = 3;   /* Flag */
    in->cfoffs[2] = 1;   /* Flag */
    in->cfoffs[3] = 2;   /* Flag */
 }

  /* Offsets of the crosscorrelation Spectral windows in input */
  SWoff = g_malloc0((in->SWArray->nwinds+2)*sizeof(olong));
  SWoff[0] = 0;
  for (iSW=0; iSW<in->SWArray->nwinds-1; iSW++) {
    SWoff[iSW+1] = SWoff[iSW] + 
      in->SWArray->winds[iSW]->numChan * in->SWArray->winds[iSW]->nCPoln * in->numAtmCorr * 2 ;
  }

  /* Notes on binary flags 
     no_flags = Napc * Nbl * Nbb * nCPoln + Nant * Nbb * nAPol
     where Napc = number of atmospheric corrected/uncorrected
     Nant = number antennas, Nbl = number of baselines, 
     Nbb = number of basebands, here assumed to be same as spectral windows
     nCPol = no. cross correlation polarizations (1,2,4)
     nAPol = number auto correlation poln (1,2,3)
     Axis order (sorta) defined in APCOrder (InitScan) but apparently in order 
     given in equation (first most slowly variable).
     NB: The bits in the flags have meaning and may need to be masked.
   */
  /* Which cross correlation IF/Spectral windows */
  if (in->coffif) g_free(in->coffif);                               /* Vis */
  in->coffif = g_malloc0((in->SWArray->nwinds+2)*sizeof(olong));    /* Vis */
  if (in->cfoffif) g_free(in->cfoffif);                             /* Flag */
  in->cfoffif = g_malloc0((in->SWArray->nwinds+2)*sizeof(olong));   /* Flag */
  /* Which ones selected? */
  inext          = 0;
  flagoff        = 0;
  in->coffif[0]  = 0;    /* Vis */
  in->cfoffif[0] = 0;   /* Flag */
  for (iSW=0; iSW<in->SWArray->nwinds; iSW++) {
    /* Use ordering */
    jSW = in->SWArray->order[iSW];
    if (in->SWArray->winds[jSW]->selected) {
      /* Correct for selection of Atm corr - 
	 this assumes Atm corr axis earlier than freq, poln*/
      in->coffif[inext] = (SWoff[jSW] +
			   in->offAtmCorr*in->SWArray->winds[jSW]->numChan*in->SWArray->winds[jSW]->nCPoln*2);
      in->cfoffif[inext] = jSW * in->numAtmCorr * in->SWArray->winds[jSW]->nCPoln; /* nAtm * nPoln */
      /*    SWoff[jSW]/(in->numSpectralChann * in->numAtmCorr * 2); */ 
      inext++;
    }  
    flagoff +=  in->numAtmCorr * in->SWArray->winds[jSW]->nCPoln; /* Offset in flagging array for next ACor */
  }
  if (SWoff) g_free(SWoff); SWoff = NULL; /* Cleanup*/
  
  /* Spectral window (IF) sidebands */
  if (in->isLSB) g_free(in->isLSB);
  in->isLSB = g_malloc0(in->SWArray->nwinds*sizeof(gboolean));
  inext = 0;
  for (iSW=0; iSW<in->SWArray->nwinds; iSW++) {
    if (in->SWArray->winds[iSW]->selected) {
      /*in->isLSB[iSW] = in->SWArray->winds[iSW]->netSideband[0]=='$';  DEBUG STUB */
      /*in->isLSB[iSW] = in->SWArray->winds[iSW]->chanFreqStep<0.0;  Do not know why */
      in->isLSB[iSW] = FALSE;  /* Do not flip spectra - thery will be labeled as lsb */
      inext++;
    }
  }

  /* Auto correlation frequency increment  RR, LL(XX,YY) are real RL (XY) complex */
  if (in->numAPoln<=2) in->aincf = in->numAPoln;    /* Vis */
  else in->aincf = 4;                               /* Vis */
  if (in->numAPoln<=2) in->afincf = in->numAPoln;   /* Flag */
  else in->afincf = 4;                              /* Flag */

  /* Auto correlation IF increment = no poln.values x no Chan*/
  in->aincif  = in->aincf * in->numSpectralChann;   /* Vis */
  in->afincif = in->aincf * in->numSpectralChann;   /* Flag */

  /* Ordering of auto correlation polarizations - shuffle order */
  if (in->aoffs) g_free(in->aoffs);                     /* Vis */
  in->aoffs = g_malloc0(in->numAPoln*sizeof(olong));    /* Vis */
  if (in->afoffs) g_free(in->afoffs);                   /* Flag */
  in->afoffs = g_malloc0(in->numAPoln*sizeof(olong));   /* Flag */
  if (in->numAPoln==1) {
    in->aoffs[0]  = 0;  /* Vis */
    in->afoffs[0] = 0;  /* Flag */
  } else if (in->numAPoln==2) {
    in->aoffs[0]  = 0;  /* Vis */
    in->aoffs[1]  = 1;  /* Vis */
    in->afoffs[0] = 0;  /* Flag */
    in->afoffs[1] = 1;  /* Flag */
  } else if ((in->numAPoln==3) || (in->numAPoln==4)) {
    in->aoffs[0] = 0;   /* Vis */
    in->aoffs[1] = 3;   /* Vis */
    in->aoffs[2] = 1;   /* Vis */
    in->afoffs[0] = 0;  /* Flag */
    in->afoffs[1] = 3;  /* Flag */
    in->afoffs[2] = 1;  /* Flag */
 }

  in->crossVisSize = 0;  /* Init vis size */
  in->autoVisSize  = 0;

  /* Offsets of the autocorrelation Spectral windows in input */
  SWoff    = g_malloc0((in->SWArray->nwinds+2)*sizeof(olong));
  SWoff[0] = 0;
  for (iSW=0; iSW<in->SWArray->nwinds-1; iSW++) {
    SWoff[iSW+1] = SWoff[iSW] + /* Note use of boolean */
      in->SWArray->winds[jSW]->selected*in->SWArray->winds[iSW]->numChan * in->aincf;
  }

  /* Which auto correlation IF/Spectral windows */
  if (in->aoffif) g_free(in->aoffif);                              /* Vis */
  in->aoffif = g_malloc0((in->SWArray->nwinds+2)*sizeof(olong));   /* Vis */
  if (in->afoffif) g_free(in->afoffif);                            /* Flag */
  in->afoffif = g_malloc0((in->SWArray->nwinds+2)*sizeof(olong));  /* Flag */
  /* Which ones selected? */
  inext         = 0;
  in->aoffif[0] = 0;   /* Flag */
  /* flagcnt = in->SWArray->winds[iSW]->nCPoln;  Offset in flagging array - initialized above */
  for (iSW=0; iSW<in->SWArray->nwinds; iSW++) {
    /* Use ordering */
    jSW = in->SWArray->order[iSW];
    if (in->SWArray->winds[jSW]->selected) {
      in->aoffif[inext]  = SWoff[jSW];  /* Vis */
      in->afoffif[inext] = flagoff + jSW * in->SWArray->winds[jSW]->nAPoln; /* nAtm * nPoln */
      inext++;
    }
   
    /* Count size of visibilities */
    in->crossVisSize += 2*in->SWArray->winds[iSW]->numChan * in->SWArray->winds[iSW]->nCPoln * in->numAtmCorr; 
    in->autoVisSize  +=   in->SWArray->winds[iSW]->numChan * in->SWArray->winds[iSW]->nAPoln;
  }
  if (SWoff) g_free(SWoff); SWoff = NULL;  /* Cleanup*/
} /* end ObitBDFDataInitScan */


 /**
 * Select Spectral windows by number of channels/band
 * Modifies in->SWArray
 * \param in       The structure to update
 * \param selChan selected number of channels
 * \param selIF   selected number of IFs (spectral windows)
 * \param band    Selected band
 * \return TRUE if some data selected
 */
gboolean ObitBDFDataSelChan  (ObitBDFData *in, olong selChan, 
			      olong selIF, ObitASDMBand band)
{
  /*gchar *routine = "ObitBDFDataSelChan";*/

  return ObitSDMDataSelChan (in->SWArray, selChan, selIF, band);
  
} /* end ObitBDFDataSelChan */

 /**
 * Parse integration info from buffer, ingest data for integration
 * May update buffer contents.
 * \param in  The object to update
 * \param err Obit error stack object.
 * \return I/O code, OBIT_IO_OK = OK, OBIT_IO_EOF=File read without finding data.
 */
ObitIOCode ObitBDFDataInitInteg  (ObitBDFData *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  gchar *startInfo, *endInfo, *prior, *next, *tstr;
  olong i, maxStr;
  gchar *routine = "ObitBDFDataInitInteg";

  /* error checks */
  if (err->error) return retCode;

  /* Create info structure if not there */
  if (in->IntegInfo==NULL) in->IntegInfo = g_malloc0(sizeof(BDFIntegInfo));

  /* Parse scan header - first find limits */
  maxStr    = in->nBytesInBuffer - (in->current-in->buffer);
  /*startInfo = g_strstr_len (in->current, maxStr, "<sdmDataSubsetHeader ");*/
  startInfo = findString (in->current, maxStr, "<sdmDataSubsetHeader ");
  /* May need next buffer */
  while (startInfo==NULL) {
    retCode = ObitBDFDataFillBuffer (in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    maxStr = in->nBytesInBuffer;
    /*startInfo = g_strstr_len (in->buffer, maxStr, "<sdmDataSubsetHeader ");*/
    startInfo = findString (in->buffer, maxStr, "<sdmDataSubsetHeader ");
   if (retCode==OBIT_IO_EOF) startInfo = in->buffer;
  }
  /* If the entire file was read and no start info - then something is wrong */
  if (retCode==OBIT_IO_EOF) return retCode;

  maxStr    = in->nBytesInBuffer - (olong)(startInfo-in->buffer);
  /*endInfo   = g_strstr_len (startInfo, maxStr, "</sdmDataSubsetHeader>");*/
  endInfo   = findString (startInfo, maxStr, "</sdmDataSubsetHeader>");
  maxStr    = (olong)(endInfo-startInfo);
  
  /* Start Time */
  prior = "<time>";
  in->IntegInfo->time = BDFparse_time(startInfo, maxStr, prior, &next);
  in->currTime = in->IntegInfo->time - in->SDMData->refJD;
  
  /* Integration Time (sec) */
  prior = "<interval>";
  in->IntegInfo->interval = BDFparse_timeint(startInfo, maxStr, prior, &next);
  in->currIntTime = in->IntegInfo->interval;
  
  /* Data type from dataStruct */
  /*startInfo = g_strstr_len (startInfo, maxStr, "<dataStruct ");*/
  startInfo = findString (startInfo, maxStr, "<dataStruct ");
  maxStr    = (olong)(endInfo-startInfo);
  prior = "type=";
  tstr = BDFparse_quote_str (startInfo, maxStr, prior, &next);
  in->IntegInfo->type = LookupDataType(tstr);
  g_free(tstr);
  in->current = next;  /* where in buffer */

  /* if multiple times, setup */
  if (in->numTime>1) {
      if (in->actualTimesData==NULL)
	in->actualTimesData = g_malloc0(sizeof(ofloat)*(in->numTime+10));
      else if (in->ScanInfo->actualTimesSize!=in->numTime)
	in->actualTimesData = g_realloc(in->actualTimesData, sizeof(ofloat)*(in->numTime+10));
      in->ScanInfo->actualTimesSize = in->numTime;
      for (i=0; i<in->numTime; i++) {
	in->actualTimesData[i] = (ofloat)(in->IntegInfo->time - in->SDMData->refJD);
	in->actualTimesData[i] +=
	  (ofloat)((i-in->numTime/2.)/(in->numTime))*((ofloat)(in->IntegInfo->interval/86400.0));
      }
      in->currTime     = in->actualTimesData[0];  /* First time */
      in->currIntTime /= in->numTime;             /* Actual integration time */
      in->currTimeIndx = 0;
  } /* End time array */

  /* Init antennas */
  in->ant1     = 0;
  in->ant2     = 1;
  in->topAnt   = 1;
  in->nextCVis = 0;   /* Cross vis index */
  in->nextAVis = 0;   /* Auto vis index */

  return retCode;
} /* end ObitBDFDataInitInteg */

 /**
 * Read integration info from buffer, ingest header, data for integration
 * May update buffer contents.
 * \param in  The object to update
 * \param err Obit error stack object.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_EOF = EOF.
 */
ObitIOCode ObitBDFDataReadInteg (ObitBDFData *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  ObitBDFMIMEType type;
  gboolean byteFlip;
  ofloat scale;
  ObitBDFDataType Dtype;
  gchar *last, *start;
  gchar *routine = "ObitBDFDataReadInteg";

  /* error checks */
  if (err->error) return retCode;

  /* Initialize */
  in->haveCrossData       = FALSE;
  in->haveAutoData        = FALSE;
  in->haveActualTimes     = FALSE;
  in->haveFlag            = FALSE;
  in->haveActualDurations = FALSE;
  in->haveWeight          = FALSE;
  in->haveZeroLag         = FALSE;
  
  /* Parse header */
  retCode = ObitBDFDataInitInteg (in, err);
  if (retCode==OBIT_IO_EOF) return retCode;
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Byte flip needed */
  byteFlip =  ((in->ScanInfo->endian==BDFEndian_Big) && (G_BYTE_ORDER==G_LITTLE_ENDIAN)) ||
    ((in->ScanInfo->endian==BDFEndian_Little) && (G_BYTE_ORDER==G_BIG_ENDIAN));

  /* scale - this may not be general enough */
  scale = 1.0;
  if (in->isALMA) scale = 1.0 / in->ScanInfo->BBinfo[0]->SWinds[0]->scaleFactor;

  /* Loop through data segments until next header */
  while (1) {
    /* Type and start of next segment */ 
    last = in->current;
    type = GetNextMIME (in, last, &start, err);
    in->current = start;
    if ((type==BDFMIMEType_sdmDataHeader) || (type==BDFMIMEType_desc)) break;
    /* Through? */
    if (type==BDFMIMEType_EOF) retCode = OBIT_IO_EOF;
    if ((type==BDFMIMEType_EOF) || (type==BDFMIMEType_Unknown)) break;
    if (retCode==OBIT_IO_EOF) break;

    /* Cross correlation data */
    if (type==BDFMIMEType_crossData) {
      /* Big enough? */
      if (in->crossCorr && (in->ScanInfo->crossDataSize>in->nCrossCorr)) {
	in->crossCorr = g_realloc(in->crossCorr, (in->ScanInfo->crossDataSize+10)*sizeof(ofloat));
      }
      in->nCrossCorr = in->ScanInfo->crossDataSize;
      in->haveCrossData = TRUE;
      /* Create if needed */
      if (in->crossCorr==NULL)
	in->crossCorr = g_malloc0((in->nCrossCorr+10)*sizeof(ofloat));
      Dtype = in->IntegInfo->type;
      scale = 1.0;
      if (in->isALMA) scale = 1.0 / in->ScanInfo->BBinfo[0]->SWinds[0]->scaleFactor;
      retCode = CopyFloats (in, start, in->crossCorr, in->nCrossCorr, byteFlip, scale, Dtype, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      if (retCode==OBIT_IO_EOF) return retCode;
      continue;
    } /* End cross correlation */
    
    /* Auto correlation data */
    if (type==BDFMIMEType_autoData) {
      /* Big enough? */
      if (in->autoCorr && (in->ScanInfo->autoDataSize>in->nAutoCorr)) {
	in->autoCorr = g_realloc(in->autoCorr, (in->ScanInfo->autoDataSize+10)*sizeof(ofloat));
      }
      in->nAutoCorr = in->ScanInfo->autoDataSize;
      in->haveAutoData = TRUE;
      /* Create if needed */
      if (in->autoCorr==NULL)
	in->autoCorr = g_malloc0((in->nAutoCorr+10)*sizeof(ofloat));
      if (in->isALMA) Dtype = BDFDataType_FLOAT32_TYPE;
      else Dtype = in->IntegInfo->type;
      scale = 1.0;
      if (in->isALMA) scale = 1.0 / in->ScanInfo->BBinfo[0]->SWinds[0]->scaleFactor;
      retCode = CopyFloats (in, start, in->autoCorr, in->nAutoCorr, byteFlip, scale, Dtype, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      if (retCode==OBIT_IO_EOF) return retCode;
      continue;
   } /* End auto correlation */
    
    /* flag data - these are really longs */
    if (type==BDFMIMEType_flags) {
      in->haveFlag = TRUE;
      /* Big enough? */
      if (in->flagData && (in->ScanInfo->FlagSize>in->nFlagData)) {
	in->flagData = g_realloc(in->flagData, (in->ScanInfo->FlagSize+10)*sizeof(olong));
      }
      in->nFlagData = in->ScanInfo->FlagSize;
      /* Size of XCor entry=  Napc * Nbb * nCPoln */
      in->FlagSize = in->numAtmCorr * in->numBaseband * in->numCPoln;
      /* Create if needed */
      if (in->flagData==NULL)
	/*in->flagData = g_malloc0(sizeof(olong)*(in->ScanInfo->FlagSize+10)/8);*/
	in->flagData = g_malloc0(sizeof(olong)*(in->ScanInfo->FlagSize+10));
      Dtype = BDFDataType_INT32_TYPE;
      retCode = CopyLongs (in, start, in->flagData, in->ScanInfo->FlagSize, 
			   byteFlip, Dtype, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      if (retCode==OBIT_IO_EOF) return retCode;
      continue;
  } /* End flag data */
    
    /* ActualTimes data */
    if (type==BDFMIMEType_actualTimes) {
      in->haveActualTimes = TRUE;
      /* Create if needed */
      if (in->actualTimesData==NULL)
	/*in->actualTimesData = g_malloc0(sizeof(ofloat)*(in->ScanInfo->actualTimesSize+10)/8);*/
	in->actualTimesData = g_malloc0(sizeof(ofloat)*(in->ScanInfo->actualTimesSize+10));
      Dtype = in->IntegInfo->type;
      scale = 1.0;
      if (in->isALMA) scale = 1.0 / in->ScanInfo->BBinfo[0]->SWinds[0]->scaleFactor;
      /*retCode = CopyFloats (in, start, (ofloat*)in->actualTimesData, in->ScanInfo->actualTimesSize/8, */
      retCode = CopyFloats (in, start, (ofloat*)in->actualTimesData, in->ScanInfo->actualTimesSize, 
			    byteFlip, scale, Dtype, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      if (retCode==OBIT_IO_EOF) return retCode;
      continue;
    } /* End actualTimes data */
     
    /* ActualDurations data */
    if (type==BDFMIMEType_actualDurations) {
      in->haveActualDurations = TRUE;
      /* Create if needed */
      if (in->actualDurationsData==NULL)
	in->actualDurationsData = g_malloc0(sizeof(olong)*(in->ScanInfo->actualDurationsSize+10));
      Dtype = in->IntegInfo->type;  scale = 1.0;
      scale = 1.0;
      if (in->isALMA) scale = 1.0 / in->ScanInfo->BBinfo[0]->SWinds[0]->scaleFactor;
      retCode = CopyFloats (in, start, (ofloat*)in->actualDurationsData, in->ScanInfo->actualDurationsSize, 
			    byteFlip, scale, Dtype, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      if (retCode==OBIT_IO_EOF) return retCode;
      continue;
    } /* End actualDurations data */
    
    /* weight data */
    if (type==BDFMIMEType_weights) {
      in->haveWeight = TRUE;
      /* Create if needed */
      if (in->weightData==NULL)
	in->weightData = g_malloc0(sizeof(olong)*(in->ScanInfo->weightSize+10));
      Dtype = in->IntegInfo->type; scale = 1.0;
      scale = 1.0;
      if (in->isALMA) scale = 1.0 / in->ScanInfo->BBinfo[0]->SWinds[0]->scaleFactor;
      retCode = CopyFloats (in, start, (ofloat*)in->weightData, in->ScanInfo->weightSize, byteFlip, 
			    scale, Dtype, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      if (retCode==OBIT_IO_EOF) return retCode;
      continue;
    } /* End weight data */

    /* zero lags data */
    if (type==BDFMIMEType_zeroLags) {
      in->haveZeroLag = TRUE;
      /* Create if needed */
      if (in->zeroLagData==NULL)
	in->zeroLagData = g_malloc0(sizeof(olong)*(in->ScanInfo->zeroLagSize+10));
      if (in->isALMA) Dtype = BDFDataType_FLOAT32_TYPE;
      else Dtype = in->IntegInfo->type;  scale = 1.0;
      scale = 1.0;
      if (in->isALMA) scale = 1.0 / in->ScanInfo->BBinfo[0]->SWinds[0]->scaleFactor;
      retCode = CopyFloats (in, start, (ofloat*)in->zeroLagData, in->ScanInfo->zeroLagSize, byteFlip, 
			    scale, Dtype, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      if (retCode==OBIT_IO_EOF) return retCode;
      continue;
    } /* End zero lags data */

 }  /* end loop over data segments */

  return retCode;
} /* end ObitBDFDataReadInteg */

/**
 * Return single visibility record from the BDF
 * \param in      Pointer to object to be read.
 * \param vis     [output] visibility record
 *                <0 => start at current position.
 * \param size    number of bytes to read.
 * \param buffer  pointer to buffer to accept results.
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_EOF = EOF.
 */
ObitIOCode ObitBDFDataGetVis (ObitBDFData *in, ofloat *vis, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_EOF;
  olong iStok, iChan, iIF, indx, ondx, kndx, voff, foff, jChan, ant1, ant2, suba=1;
  olong nChan, nIF, nStok;
  ofloat weight;
  /*gchar *routine = "ObitBDFDataGetVis";*/

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (vis != NULL);

  /* Sizes of things */
  nChan = in->desc->inaxes[in->desc->jlocf];
  nIF   = in->desc->inaxes[in->desc->jlocif];
  nStok = in->desc->inaxes[in->desc->jlocs];

  /* More cross data? */
  if (in->haveCrossData) {
    
    retCode = OBIT_IO_OK;  /* Have some */

    /* Random parameters */
    vis[in->desc->ilocu]  = 0.0;
    vis[in->desc->ilocv]  = 0.0;
    vis[in->desc->ilocw]  = 0.0;
    if (in->desc->ilocfq>=0) vis[in->desc->ilocfq] = 1.0;
    vis[in->desc->iloct]  = in->currTime;
    vis[in->desc->ilocit] = in->currIntTime;
    vis[in->desc->ilocsu] = (ofloat)in->sourceNo;
    ant1 = in->antNo[in->ant1];
    ant2 = in->antNo[in->ant2];
    weight = vis[in->desc->ilocit];  /* Use integration time as weight */

    /* Is the order of the baseline correct, ant1<ant2 ? */
    if (ant1<ant2) {
      ObitUVDescSetAnts (in->desc, vis, ant1, ant2, suba);
    
      /* Loop over visibilities */
      voff = in->nextCVis * in->crossVisSize;
      /* in->nextCVis x flagSize =  Napc * Nbb * nCPoln*/
      foff = in->nextCVis * in->FlagSize;
      for (iIF=0; iIF<nIF; iIF++) {
	for (iChan=0; iChan<nChan; iChan++) {
	  /* Need to reverse order for LSB? */
	  if (in->isLSB[iIF]) jChan = nChan-iChan-1;
	  else jChan = iChan;
	  for (iStok=0; iStok<nStok; iStok++) {
	    ondx = in->desc->nrparm +
	      iStok*in->desc->incs + jChan*in->desc->incf + iIF*in->desc->incif;
	    indx = voff + in->coffs[iStok] + iChan*in->cincf + in->coffif[iIF];
	    vis[ondx]   = in->crossCorr[indx];
	    vis[ondx+1] = in->crossCorr[indx+1];
	    vis[ondx+2] = weight;
	    /* Binary flagging per baseline/IF/Stokes (not chan)*/
	    if (in->haveFlag && in->binFlag) {
	      kndx = foff + in->cfoffs[iStok] + in->cfoffif[iIF];
	      /* May need to mask unwanted bits */
	      if (in->flagData[kndx]!=0) vis[ondx+2] = -fabs(vis[ondx+2]);
	    }
	  } /* end Stokes loop */
	} /* end Channel loop */
      } /* end IF loop */
    } else {   /* Flip baseline */
      ObitUVDescSetAnts (in->desc, vis, ant2, ant1, suba);
      
      /* Loop over visibilities */
      voff = in->nextCVis * in->crossVisSize;
      foff = in->nextCVis * (in->crossVisSize/(2*nChan));
      for (iIF=0; iIF<nIF; iIF++) {
	for (iChan=0; iChan<nChan; iChan++) {
	  /* Need to reverse order for LSB? */
	  if (in->isLSB[iIF]) jChan = nChan-iChan-1;
	  else jChan = iChan;
	  for (iStok=0; iStok<nStok; iStok++) {
	    ondx = in->desc->nrparm +
	      iStok*in->desc->incs + jChan*in->desc->incf + iIF*in->desc->incif;
	    indx = voff + iChan*in->cincf + in->coffif[iIF];
	    /* Use switch to deal with different polns */
	    switch (iStok) { 
	    case 0:     /* RR or XX */
	    case 1:     /* LL or YY */
	      vis[ondx]   =  in->crossCorr[in->coffs[iStok]+indx];
	      vis[ondx+1] = -in->crossCorr[in->coffs[iStok]+indx+1];   /* Conjugate */
	      vis[ondx+2] = weight;
	      break;
	    case 2:  /* Swap RL, LR */
	      vis[ondx]   =  in->crossCorr[in->coffs[iStok+1]+indx];
	      vis[ondx+1] = -in->crossCorr[in->coffs[iStok+1]+indx+1];   /* Conjugate */
	      vis[ondx+2] = weight;
	      break;
	    case 3:
	      vis[ondx]   =  in->crossCorr[in->coffs[iStok-1]+indx];
	      vis[ondx+1] = -in->crossCorr[in->coffs[iStok-1]+indx+1];   /* Conjugate */
	      vis[ondx+2] = weight;
	      break;
	    default:
	      g_assert_not_reached(); /* unknown, barf */
	    }; /* end switch on polarization */
	    /* Binary flagging */
	    if (in->haveFlag && in->binFlag) {
	      kndx = foff + in->cfoffs[iStok] + in->cfoffif[iIF];
	      if (in->flagData[kndx]!=0) vis[ondx+2] = -fabs(vis[ondx+2]);
	    }
	  } /* end Stokes loop */
	} /* end Channel loop */
      } /* end IF loop */
    } /* end reverse baseline */

    /* Update visibility and antenna numbers */
    in->nextCVis++;
    in->ant1++;
    if (in->ant1>=in->topAnt) {in->topAnt++; in->ant1=0; in->ant2 = in->topAnt;}
    /* Finished? */
    if (in->topAnt>=in->nant) in->haveCrossData = FALSE;
    if (!in->haveCrossData) {in->ant1 = 0; in->ant2 = 0; in->topAnt = 0;}
    return retCode;
  } /* end have cross data */

   
  /* More auto data? */
  if (in->haveAutoData) {
    
    retCode = OBIT_IO_OK;  /* Have some */

    /* Random parameters */
    vis[in->desc->ilocu]  = 0.0;
    vis[in->desc->ilocv]  = 0.0;
    vis[in->desc->ilocw]  = 0.0;
    vis[in->desc->iloct]  = in->currTime;
    vis[in->desc->ilocit] = in->currIntTime;
    ObitUVDescSetAnts (in->desc, vis,in->antNo[in->ant1], in->antNo[in->ant1], suba);
    vis[in->desc->ilocsu] = (ofloat)in->sourceNo;
    weight = vis[in->desc->ilocit];  /* Use integration time as weight */
    
    /* Loop over visibilities */
    voff = in->nextAVis * in->autoVisSize;
    foff = (in->nextCVis-1) * in->FlagSize + (voff/nChan);  /* Assumes cross flags before auto */
    for (iIF=0; iIF<nIF; iIF++) {
      for (iChan=0; iChan<nChan; iChan++) {
	/* Need to reverse order for LSB? */
	if (in->isLSB[iIF]) jChan = nChan-iChan-1;
	else jChan = iChan;
	for (iStok=0; iStok<nStok; iStok++) {
	  ondx = in->desc->nrparm +
	    iStok*in->desc->incs + jChan*in->desc->incf + iIF*in->desc->incif;
	  indx = voff + jChan*in->aincf + in->aoffif[iIF];
	  /* Use switch to deal with different polns */
	  switch (iStok) { 
	  case 0:
	    vis[ondx]   = in->autoCorr[indx+in->aoffs[0]];
	    vis[ondx+1] = 0.0;
	    vis[ondx+2] = weight;
	    break;
	  case 1:
	    vis[ondx] = in->autoCorr[indx+in->aoffs[1]];
	    vis[ondx+1] = 0.0;
	    vis[ondx+2] = weight;
	    break;
	  case 2:
	    vis[ondx]   = in->autoCorr[indx+in->aoffs[2]];
	    vis[ondx+1] = in->autoCorr[indx+in->aoffs[2]+1];
	    vis[ondx+2] = weight;
	    break;
	  case 3:
	    /* Use conjugate */
	    vis[ondx]   =  in->autoCorr[indx+in->aoffs[2]];
	    vis[ondx+1] = -in->autoCorr[indx+in->aoffs[2]+1];
	    vis[ondx+2] = weight;
	    break;
	  default:
	    g_assert_not_reached(); /* unknown, barf */
	  }; /* end switch on polarization */
	    /* Binary flagging */
	  if (in->haveFlag && in->binFlag) {
	      kndx = foff + in->afoffs[iStok] + in->afoffif[iIF];
	      if (in->flagData[kndx]!=0) vis[ondx+2] = -fabs(vis[ondx+2]);
	    }
	} /* end Stokes loop */
      } /* end Channel loop */
    } /* end IF loop */
    
    /* Update visibility and antenna numbers */
    in->nextAVis++;
    in->ant1++;
    /* Finished? */
    if (in->ant1>=in->nant) {
      in->currTimeIndx++;
      in->ant1 = 0;
      if (in->currTimeIndx>=in->numTime) in->haveAutoData = FALSE;
      else {
	/* Next time */
	in->currTime = in->actualTimesData[in->currTimeIndx];
      }
    }
    return retCode;
  } /* end have auto data */
    
  return retCode;
} /* end ObitBDFDataGetVis  */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitBDFDataClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitBDFDataClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitBDFDataClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitBDFDataClassInfoDefFn (gpointer inClass)
{
  ObitBDFDataClassInfo *theClass = (ObitBDFDataClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitBDFDataClassInit;
  theClass->newObit       = (newObitFP)newObitBDFData;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitBDFDataClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitBDFDataGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitBDFDataCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitBDFDataClear;
  theClass->ObitInit      = (ObitInitFP)ObitBDFDataInit;
  theClass->ObitBDFDataCreate = (ObitBDFDataCreateFP)ObitBDFDataCreate;

} /* end ObitBDFDataClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitBDFDataInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitBDFData *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->SDMData             = NULL;
  in->SWArray             = NULL;
  in->ScanInfo            = NULL;
  in->IntegInfo           = NULL;
  in->DataFile            = NULL;
  in->buffer              = NULL;
  in->crossCorr           = NULL;
  in->autoCorr            = NULL;
  in->flagData            = NULL;
  in->actualTimesData     = NULL;
  in->actualDurationsData = NULL;
  in->weightData          = NULL;
  in->zeroLagData         = NULL;
  in->antNo               = NULL;
  in->antId               = NULL;
  in->coffs               = NULL;
  in->coffif              = NULL;
  in->cfoffs              = NULL;
  in->cfoffif             = NULL;
  in->aoffs               = NULL;
  in->aoffif              = NULL;
  in->afoffs              = NULL;
  in->afoffif             = NULL;
  in->isLSB               = NULL;
  in->curSource           = NULL;
  in->isEVLA              = FALSE;
  in->isALMA              = FALSE;
  in->binFlag             = FALSE;
  in->selAtmCorr          = FALSE;
  in->numAtmCorr          = 1;
} /* end ObitBDFDataInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitBDFData* cast to an Obit*.
 */
void ObitBDFDataClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitBDFData *in = inn;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Close file */
  err = newObitErr();
  if (in->file) ObitFileClose (in->file, err);
  err = ObitErrUnref(err);

  /* delete this class members */
  in->file      = ObitFileUnref(in->file);
  in->ScanInfo  = KillBDFScanInfo(in->ScanInfo);
  in->IntegInfo = KillBDFIntegInfo(in->IntegInfo);
  in->desc      = ObitUVDescUnref(in->desc);
  in->SDMData   = ObitSDMDataUnref(in->SDMData);
  in->SWArray   = ObitSDMDataKillSWArray(in->SWArray);
  if (in->DataFile)            g_free(in->DataFile);
  if (in->buffer)              g_free(in->buffer);
  if (in->crossCorr)           g_free(in->crossCorr);
  if (in->autoCorr)            g_free(in->autoCorr);
  if (in->flagData)            g_free(in->flagData);
  if (in->actualTimesData)     g_free(in->actualTimesData);
  if (in->actualDurationsData) g_free(in->actualDurationsData);
  if (in->weightData)          g_free(in->weightData);
  if (in->zeroLagData)         g_free(in->zeroLagData);
  if (in->antNo)               g_free(in->antNo);
  if (in->antId)               g_free(in->antId);
  if (in->coffs)               g_free(in->coffs);
  if (in->coffif)              g_free(in->cfoffif);
  if (in->cfoffs)              g_free(in->cfoffs);
  if (in->cfoffif)             g_free(in->coffif);
  if (in->aoffs)               g_free(in->aoffs);
  if (in->aoffif)              g_free(in->aoffif);
  if (in->afoffs)              g_free(in->afoffs);
  if (in->afoffif)             g_free(in->afoffif);
  if (in->isLSB)               g_free(in->isLSB);
  if (in->curSource)           g_free(in->curSource);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitBDFDataClear */

/* BDF XML routines */
/**  Parse integer from XLM string 
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, 0 if problem
 */
static olong BDFparse_int(gchar *string, olong maxChar, 
			  gchar *prior, gchar **next)
{
  olong out = 0;
  gchar *b;

  *next = string;  /* if not found */
  /*b = g_strstr_len (string, maxChar, prior);*/
  b = findString (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  out = (olong)strtol(b, next, 10);
    
  return out;
} /* end BDFparse_int */

/**  Parse unquoted string from XLM string 
 * All text from end of prior until next '<'
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 */
static gchar* BDFparse_str(gchar *string, olong maxChar, 
			   gchar *prior, gchar **next)
{
  gchar *out = NULL;
  gchar *b;
  olong charLeft;
  olong i, n;

  *next = string;  /* if not found */
  /*b = g_strstr_len (string, maxChar, prior);*/
  b = findString (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);

  /* count */
  charLeft = maxChar - (b-string);
  n = 0;
  for (i=0; i<charLeft; i++) {
    if (b[i]=='<') break;
    n++;
  }
  out = g_malloc(n+1);
  for (i=0; i<n; i++) out[i] = b[i]; out[i] = 0;
  *next = b + n;

  return out;
} /* end BDFparse_str */

/**  Parse double quoted from XLM string 
 * All text from end of prior+1 until next '"'
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 */
static gchar* BDFparse_quote_str(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next)
{
  gchar *out = NULL;
  gchar *b;
  olong charLeft;
  olong i, n;

  *next = string;  /* if not found */
  /*b = g_strstr_len (string, maxChar, prior);*/
  b = findString (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  if (*b!='"') return out;  /* Make sure quote */
  b++;                      /* Skip quote */

  /* count */
  charLeft = maxChar - (b-string);
  n = 0;
  for (i=0; i<charLeft; i++) {
    if (b[i]=='"') break;
    n++;
  }
  out = g_malloc(n+1);
  for (i=0; i<n; i++) out[i] = b[i]; out[i] = 0;
  *next = b + n;

  return out;
} /* end BDFparse_quote_str */

/**  Parse array of unquoted strings from XLM string;
 * array of strings NULL terminated.
 * Each string text from end of prior until next ' ' or '<' or '"'
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done, last entry NULL
 */
static gchar** BDFparse_strarray(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next)
{
  gchar **out = NULL;
  gchar *b;
  olong charLeft, num;
  olong i, j, n;

  *next = string;  /* if not found */
  /*b = g_strstr_len (string, maxChar, prior);*/
  b = findString (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);

  /* Count blanks before next '"' (num val -1) */
  num = 1;  /* No blank after last */
  for (j=1; j<maxChar; j++) {
    if (b[j]==' ') num++;
    if (b[j]=='"') break;
  }
  out = g_malloc0(MAX(1,(num+1))*sizeof(gchar*));

  /* Loop over strings */
  for (j=0; j<num; j++) {
    /* count */
    charLeft = maxChar - (b-string);
    n = 0;
    for (i=0; i<charLeft; i++) {
      if ((b[i]=='<') || (b[i]==' ') || (b[i]=='"'))break;
      n++;
    }
    out[j] = g_malloc(n+1);
    for (i=0; i<n; i++) out[j][i] = b[i]; out[j][i] = 0;
    b += n+1;
    *next = b;
  } /* end loop over strings */

  /* NULL in highest */
  out[num] = NULL;
  return out;
} /* end BDFparse_strarray */

/**  Parse axis order array from XML string  
 * Each string text from end of prior until next ' ' or '"'
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 * one larger than actual, last = -999
 */
static ObitBDFAxisName* BDFparse_axesarray(gchar *string, olong maxChar, 
					   gchar *prior, gchar **next)
{
  ObitBDFAxisName *out = NULL;
  gchar *b;
  olong charLeft, num;
  olong i, j, n;
  
  *next = string;  /* if not found */
  /*b = g_strstr_len (string, maxChar, prior);*/
  b = findString (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  
  /* Count blanks before next '"' (num val -1) */
  num = 1;
  for (j=1; j<maxChar; j++) {
    if (b[j]==' ') num++;
    if (b[j]=='"') break;
  }
  out = g_malloc0(MAX(1,(num+1))*sizeof(ObitBDFAxisName));
  
  /* Loop over strings */
  b++;  /* Skip first quote */
  for (j=0; j<num; j++) {
    /* count */
    charLeft = maxChar - (b-string);
    n = 0;
    for (i=0; i<charLeft; i++) {
      if ((b[i]=='"') || (b[i]==' '))break;
      n++;
    }
    out[j] = LookupAxisName(b);
    b += n+1;
    *next = b;
  } /* end loop over strings */

  out[num] = BDFAxisName_END;  /* mark end */
  return out;
} /* end BDFparse_axesarray */

/**  Parse time from XLM string  
 * Read time as a mjd nanoseconds, return as JD
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, 0.0 if problem
 */
static odouble BDFparse_time(gchar *string, olong maxChar, 
			     gchar *prior, gchar **next)
{
  odouble out = 0.0;
  long long temp;
  gchar *b;
  odouble mjdJD0=2400000.5; /* JD of beginning of mjd time */

  *next = string;  /* if not found */
  /*b = g_strstr_len (string, maxChar, prior);*/
  b = findString (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  temp = strtoll(b, next, 10);
  out = (odouble)((temp*1.0e-9)/86400.0) + mjdJD0;
    
  return out;
} /* end BDFparse_time */

/**  Parse time interval from XLM string  
 * Read time interval in nanoseconds, return as seconds
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, 0.0 if problem
 */
static odouble BDFparse_timeint(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next)
{
  odouble out = 0.0;
  long long temp;
  gchar *b;

  *next = string;  /* if not found */
  /*b = g_strstr_len (string, maxChar, prior);*/
  b = findString (string, maxChar, prior);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prior);
  temp = strtoll(b, next, 10);
  out = (odouble)(temp*1.0e-9);
    
  return out;
} /* end BDFparse_timeint */

/**  Look up axis name enum 
 * \param  string  String to lookup
 * \return value
 */
static ObitBDFAxisName LookupAxisName(gchar *name)
{
  ObitBDFAxisName out = 0;
 
  if (!strncmp (name, "TIM", 3)) return BDFAxisName_TIM;
  if (!strncmp (name, "BAL", 3)) return BDFAxisName_BAL;
  if (!strncmp (name, "ANT", 3)) return BDFAxisName_ANT;
  if (!strncmp (name, "BAB", 3)) return BDFAxisName_BAB;
  if (!strncmp (name, "SPW", 3)) return BDFAxisName_SPW;
  if (!strncmp (name, "SIB", 3)) return BDFAxisName_SIB;
  if (!strncmp (name, "SUB", 3)) return BDFAxisName_SUB;
  if (!strncmp (name, "BIN", 3)) return BDFAxisName_BIN;
  if (!strncmp (name, "APC", 3)) return BDFAxisName_APC;
  if (!strncmp (name, "SPP", 3)) return BDFAxisName_SPP;
  if (!strncmp (name, "POL", 3)) return BDFAxisName_POL;
  if (!strncmp (name, "STO", 3)) return BDFAxisName_STO;
  if (!strncmp (name, "HOL", 3)) return BDFAxisName_HOL;
  return out;
} /* end LookupAxisName */

/**  Look up data type enum 
 * \param  name  String to lookup, NULL => BDFDataType_FLOAT32_TYPE
 * \return value
 */
static ObitBDFDataType LookupDataType(gchar *name)
{
  ObitBDFDataType out = BDFDataType_FLOAT32_TYPE;

  if (name==NULL) return out;
 
  if (!strncmp (name, "INT16_TYPE", 10)) return BDFDataType_INT16_TYPE;
  if (!strncmp (name, "INT32_TYPE", 10)) return BDFDataType_INT32_TYPE;
  if (!strncmp (name, "INT64_TYPE", 10)) return BDFDataType_INT64_TYPE;
  if (!strncmp (name, "FLOAT32_TYPE", 12)) return BDFDataType_FLOAT32_TYPE;
  if (!strncmp (name, "FLOAT64_TYPE", 12)) return BDFDataType_FLOAT64_TYPE;
  return out;
} /* end LookupDataType */

/**  Look up data type enum 
 * \param  name  String to lookup
 * \return value
 */
static ObitBDFEndian LookupEndian(gchar *name)
{
  ObitBDFEndian out = 0;
 
  if (!strncmp (name, "Little_Endian", 13)) return BDFEndian_Little;
  if (!strncmp (name, "Big_Endian",    10)) return BDFEndian_Big;
  return out;
} /* end LookupEndian */

/** 
 * Destructor for .BDFSpecWindowInfo
 * \param  structure to destroy
 * \return NULL pointer
 */
static BDFSpecWindowInfo* KillBDFSpecWindowInfo(BDFSpecWindowInfo *info)
{
  olong i;
  if (info == NULL) return NULL;
  if (info->sideband)     g_free(info->sideband);
  if (info->sdPolProducts) {
    i = 0;
    while (info->sdPolProducts[i]) {
      g_free(info->sdPolProducts[i++]);
    }
    g_free(info->sdPolProducts);
  }
  if (info->crossPolProducts) {
    i = 0;
    while (info->crossPolProducts[i]) {
      g_free(info->crossPolProducts[i++]);
    }
    g_free(info->crossPolProducts);
  }
  g_free(info);
  return NULL;
} /* end  KillBDFSpecWindowInfo */

/** 
 * Destructor for .BDFBasebandInfo
 * \param  structure to destroy
 * \return NULL pointer
 */
static BDFBasebandInfo* KillBDFBasebandInfo(BDFBasebandInfo *info)
{
  olong i;
  if (info == NULL) return NULL;
  if (info->basebandName) g_free(info->basebandName);
  for (i = 0; i<MAXBBSW; i++) {
    if (info->SWinds[i])  KillBDFSpecWindowInfo(info->SWinds[i]);
    else break;
  }
  g_free(info);
  return NULL;
} /* end  KillBDFBasebandInfo */

/** 
 * Destructor for BDFScanInfo
 * \param  structure to destroy
 * \return NULL pointer
 */
static BDFScanInfo* KillBDFScanInfo(BDFScanInfo *info)
{
  olong i;

  if (info == NULL) return NULL;
  if ((info->BBinfo)  && (info->numBaseband>0)) {
    /*for (i=0; i<info->numBaseband; i++) {*/
    for (i=0; i<MAXBBINFO; i++) {
      info->BBinfo[i] = KillBDFBasebandInfo(info->BBinfo[i]);
    }
  }
  if (info->FlagAxes)            g_free(info->FlagAxes);
  if (info->actualTimesAxes)     g_free(info->actualTimesAxes);
  if (info->actualDurationsAxes) g_free(info->actualDurationsAxes);
  if (info->crossDataAxes)       g_free(info->crossDataAxes);
  if (info->autoDataAxes)        g_free(info->autoDataAxes);
  if (info->weightAxes)          g_free(info->weightAxes);
  if (info->zeroLagAxes)         g_free(info->zeroLagAxes);
  g_free(info);
  return NULL;
} /* end KillBDFScanInfo */

/** 
 * Destructor for .BDFIntegInfo
 * \param  structure to destroy
 * \return NULL pointer
 */
static BDFIntegInfo* KillBDFIntegInfo(BDFIntegInfo *info)
{
  if (info == NULL) return NULL;
  g_free(info);
  return NULL;
} /* end KillBDFIntegInfo */


 /**
 * Copy binary and convert to float data 
 * May update buffer contents.
 * \param in       The object to update
 * \param start    input array
 * \param target   output array
 * \param n        Number of floats
 * \param byteFlip If TRUE flip bytes
 * \param scale    scaling factor
 * \param Dtype    Data type (short, float...)
 * \param err   Obit error stack object.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_EOF = EOF.
 */
static ObitIOCode CopyFloats (ObitBDFData *in, 
			      gchar *start, ofloat *target, olong n, 
			      gboolean byteFlip, ofloat scale, 
			      ObitBDFDataType Dtype, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  olong i, nleft, ncopy, ncopyb, nhere, shit;
  ofloat *out;
  short *sdata;
  gint32 *ldata;
  char  *bdata, btemp;
  gchar *lstart;
  union fequiv inu, outu;
  gchar *routine = "CopyFloats";

  /* error checks */
  if (err->error) return retCode;

  /* Copy data by type */
  if (Dtype==BDFDataType_INT16_TYPE) {  /* scaled shorts */
    /* How many in current buffer load */
    nhere = (in->nBytesInBuffer - (olong)(start-in->buffer))/2;
    nleft = n;
    out = target;
    lstart = start;
    /* All in buffer? */
    if (nhere>=n) {  /* All in current */
      sdata = (short*)lstart;
      ncopy = n*2;
      /* byte flip if necessary */
      if (byteFlip) {
	bdata = (char*)sdata;
	for (i=0; i<ncopy; i+=2) {
	  btemp      = bdata[i];
	  bdata[i]   = bdata[i+1];
	  bdata[i+1] = btemp;
	}
      } /* end byte flip */
      for (i=0; i<n; i++) out[i] = (ofloat)sdata[i];
      in->current = lstart + ncopy;
    } else {         /* Multiple buffers */
      while (nleft>0) {
	/* Copy what's here */
	ncopy = MIN (nleft, nhere);
	ncopyb = ncopy*2;  /* in bytes */
	sdata = (short*)lstart;
	/* byte flip if necessary */
	if (byteFlip) {
	  bdata = (char*)sdata;
	  for (i=0; i<ncopyb; i+=2) {
	    btemp      = bdata[i];
	    bdata[i]   = bdata[i+1];
	    bdata[i+1] = btemp;
	  }
	} /* end byte flip */
	for (i=0; i<ncopy; i++) out[i] = (ofloat)sdata[i];
	out += ncopy;
	shit = (olong)(out-target);
	in->current = lstart + ncopyb;
	nleft -= ncopy;
	if (nleft<=0) break;  /* Done? */
	/* Get more */
	retCode = ObitBDFDataFillBuffer (in, err);
	if (err->error) {
	  Obit_traceback_val (err, routine, in->name, retCode);
	}
	/* If EOF and not done - bail */
	if ((retCode==OBIT_IO_EOF) && (nleft>0)) return retCode;
	lstart = in->current;
	nhere = (in->nBytesInBuffer - (olong)(lstart-in->buffer))/2;
      }  /* end loop over buffers */
    } /* end multiple buffers */
    
    /* If in last segment of buffer, update */
    if ((olong)(in->current-in->buffer)>(BDFBUFFERSIZE*(BDFBUFFERFRAMES-1))) {
      retCode = ObitBDFDataFillBuffer (in, err);
      if (err->error) {
	Obit_traceback_val (err, routine, in->name, retCode);
      }
    }
    /* end 16 bit integer */
  /* Scaled 32 bit integers */
  } else if (Dtype==BDFDataType_INT32_TYPE) {  /* scaled 32 bit */
    /* How many in current buffer load */
    nhere   = (in->nBytesInBuffer - (olong)(start-in->buffer))/4;
    nleft  = n;
    out    = target;
    lstart = start;
    /* All in buffer? */
    if (nhere>=n) {  /* All in current */
      ldata = (gint32*)lstart;
      ncopy = n*4;
      /* byte flip if necessary */
      if (byteFlip) {
	bdata = (char*)ldata;
	for (i=0; i<ncopy; i+=4) {
	  btemp      = bdata[i];
	  bdata[i]   = bdata[i+1];
	  bdata[i+1] = btemp;
	}
      } /* end byte flip */
      for (i=0; i<n; i++) out[i] = (ofloat)ldata[i];
      in->current = lstart + ncopy;
    } else {         /* Multiple buffers */
      while (nleft>0) {
	/* Copy what's here */
	ncopy  = MIN (nleft, nhere);
	ncopyb = ncopy*4;  /* in bytes */
	ldata  = (gint32*)lstart;
	/* byte flip if necessary */
	if (byteFlip) {
	  bdata = (char*)ldata;
	  for (i=0; i<ncopyb; i+=4) {
	    btemp      = bdata[i];
	    bdata[i]   = bdata[i+1];
	    bdata[i+1] = btemp;
	  }
	} /* end byte flip */
	for (i=0; i<ncopy; i++) out[i] = (ofloat)ldata[i];
	out += ncopy;
	shit = (olong)(out-target);
	in->current = lstart + ncopyb;
	nleft -= ncopy;
	if (nleft<=0) break;  /* Done? */
	/* Get more */
	retCode = ObitBDFDataFillBuffer (in, err);
	if (err->error) {
	  Obit_traceback_val (err, routine, in->name, retCode);
	}
	/* If EOF and not done - bail */
	if ((retCode==OBIT_IO_EOF) && (nleft>0)) return retCode;
	lstart = in->current;
	nhere = (in->nBytesInBuffer - (olong)(lstart-in->buffer))/4;
      }  /* end loop over buffers */
    } /* end multiple buffers */
    
    /* If in last segment of buffer, update */
    if ((olong)(in->current-in->buffer)>(BDFBUFFERSIZE*(BDFBUFFERFRAMES-1))) {
      retCode = ObitBDFDataFillBuffer (in, err);
      if (err->error) {
	Obit_traceback_val (err, routine, in->name, retCode);
      }
    }
    /* end 32 bit integer */
    /* 32 bit floats */
  } else if (Dtype==BDFDataType_FLOAT32_TYPE) {
    /* How many in current buffer load */
    nhere = (in->nBytesInBuffer - (olong)(start-in->buffer))/sizeof(ofloat);
    nleft = n;
    out = target;
    lstart = start;
    /* All in buffer? */
    if (nhere>=n) {  /* All in current */
      ncopy = n*sizeof(ofloat);
      memmove (out, lstart, (size_t)ncopy); 
      in->current = lstart + ncopy;
    } else {         /* Multiple buffers */
      while (nleft>0) {
	/* Copy what's here */
	ncopy = MIN (nleft, nhere);
	ncopyb = ncopy*sizeof(ofloat);  /* in bytes */
	memmove (out, lstart, (size_t)ncopyb);
	out += ncopy;
	shit = (olong)(out-target);
	in->current = lstart + ncopyb;
	nleft -= ncopy;
	if (nleft<=0) break;  /* Done? */
	/* Get more */
	retCode = ObitBDFDataFillBuffer (in, err);
	if (err->error) {
	  Obit_traceback_val (err, routine, in->name, retCode);
	}
	/* If EOF and not done - bail */
	if ((retCode==OBIT_IO_EOF) && (nleft>0)) return retCode;
	lstart = in->current;
	nhere = (in->nBytesInBuffer - (olong)(lstart-in->buffer))/sizeof(ofloat);
      }  /* end loop over buffers */
    } /* end multiple buffers */
    
    /* If in last segment of buffer, update */
    if ((olong)(in->current-in->buffer)>(BDFBUFFERSIZE*(BDFBUFFERFRAMES-1))) {
      retCode = ObitBDFDataFillBuffer (in, err);
      if (err->error) {
	Obit_traceback_val (err, routine, in->name, retCode);
      }
    }
    /* byte flip if necessary */
    if (byteFlip) {
      out = target;
      for (i=0; i<n; i++) {
	inu.full = out[i];
	outu.parts[0] = inu.parts[3]; 
	outu.parts[1] = inu.parts[2]; 
	outu.parts[2] = inu.parts[1]; 
	outu.parts[3] = inu.parts[0]; 
	out[i] = outu.full;
      }
    } /* end byte flip */

    /* end 32 bit float */
  } else {  /* Unsupported type */
    Obit_log_error(err, OBIT_Error, 
		   "%s: Unsupported data type %d", routine, Dtype);
    return OBIT_IO_SpecErr;
  }

  /* scale if necessary */
  if (fabs(scale-1.000)>1.0e-5) {
    out = target;
    for (i=0; i<n; i++) out[i] *= scale;
  }  /* end scaling */
  return retCode;
} /* end CopyFloats */

 /**
 * Copy binary and convert to long data 
 * May update buffer contents.
 * \param in       The object to update
 * \param start    input array
 * \param target   output array
 * \param n        Number of floats
 * \param byteFlip If TRUE flip bytes
 * \param Dtype    Data type (short, float...)
 * \param err   Obit error stack object.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_EOF = EOF.
 */
static ObitIOCode CopyLongs (ObitBDFData *in, 
			     gchar *start, olong *target, olong n, 
			     gboolean byteFlip, ObitBDFDataType Dtype, 
			     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  olong i, nleft, ncopy, ncopyb, nhere, shit;
  olong  *out;
  short *sdata;
  gint32 *ldata;
  char  *bdata, btemp;
  gchar *lstart;
  union fequiv inu, outu;
  gchar *routine = "CopyLongs";

  /* error checks */
  if (err->error) return retCode;

  /* Copy data by type */
  if (Dtype==BDFDataType_INT16_TYPE) {  /* scaled shorts */
    /* How many in current buffer load */
    nhere = (in->nBytesInBuffer - (olong)(start-in->buffer))/2;
    nleft = n;
    out = target;
    lstart = start;
    /* All in buffer? */
    if (nhere>=n) {  /* All in current */
      sdata = (short*)lstart;
      ncopy = n*2;
      /* byte flip if necessary */
      if (byteFlip) {
	bdata = (char*)sdata;
	for (i=0; i<ncopy; i+=2) {
	  btemp      = bdata[i];
	  bdata[i]   = bdata[i+1];
	  bdata[i+1] = btemp;
	}
      } /* end byte flip */
      for (i=0; i<n; i++) out[i] = (olong)sdata[i];
      in->current = lstart + ncopy;
    } else {         /* Multiple buffers */
      while (nleft>0) {
	/* Copy what's here */
	ncopy = MIN (nleft, nhere);
	ncopyb = ncopy*2;  /* in bytes */
	sdata = (short*)lstart;
	/* byte flip if necessary */
	if (byteFlip) {
	  bdata = (char*)sdata;
	  for (i=0; i<ncopyb; i+=2) {
	    btemp      = bdata[i];
	    bdata[i]   = bdata[i+1];
	    bdata[i+1] = btemp;
	  }
	} /* end byte flip */
	for (i=0; i<ncopy; i++) out[i] = (olong)sdata[i];
	out += ncopy;
	shit = (olong)(out-target);
	in->current = lstart + ncopyb;
	nleft -= ncopy;
	if (nleft<=0) break;  /* Done? */
	/* Get more */
	retCode = ObitBDFDataFillBuffer (in, err);
	if (err->error) {
	  Obit_traceback_val (err, routine, in->name, retCode);
	}
	/* If EOF and not done - bail */
	if ((retCode==OBIT_IO_EOF) && (nleft>0)) return retCode;
	lstart = in->current;
	nhere = (in->nBytesInBuffer - (olong)(lstart-in->buffer))/2;
      }  /* end loop over buffers */
    } /* end multiple buffers */
    
    /* If in last segment of buffer, update */
    if ((olong)(in->current-in->buffer)>(BDFBUFFERSIZE*(BDFBUFFERFRAMES-1))) {
      retCode = ObitBDFDataFillBuffer (in, err);
      if (err->error) {
	Obit_traceback_val (err, routine, in->name, retCode);
      }
    }
    /* end 16 bit integer */
  /* Scaled 32 bit integers */
  } else if (Dtype==BDFDataType_INT32_TYPE) {  /* scaled 32 bit */
    /* How many in current buffer load */
    nhere   = (in->nBytesInBuffer - (olong)(start-in->buffer))/4;
    nleft  = n;
    out    = target;
    lstart = start;
    /* All in buffer? */
    if (nhere>=n) {  /* All in current */
      ldata = (gint32*)lstart;
      ncopy = n*4;
      /* byte flip if necessary */
      if (byteFlip) {
	bdata = (char*)ldata;
	for (i=0; i<ncopy; i+=4) {
	  btemp      = bdata[i];
	  bdata[i]   = bdata[i+1];
	  bdata[i+1] = btemp;
	}
      } /* end byte flip */
      for (i=0; i<n; i++) out[i] = (olong)ldata[i];
      in->current = lstart + ncopy;
    } else {         /* Multiple buffers */
      while (nleft>0) {
	/* Copy what's here */
	ncopy  = MIN (nleft, nhere);
	ncopyb = ncopy*4;  /* in bytes */
	ldata  = (gint32*)lstart;
	/* byte flip if necessary */
	if (byteFlip) {
	  bdata = (char*)ldata;
	  for (i=0; i<ncopyb; i+=4) {
	    btemp      = bdata[i];
	    bdata[i]   = bdata[i+1];
	    bdata[i+1] = btemp;
	  }
	} /* end byte flip */
	for (i=0; i<ncopy; i++) out[i] = (olong)ldata[i];
	out += ncopy;
	shit = (olong)(out-target);
	in->current = lstart + ncopyb;
	nleft -= ncopy;
	if (nleft<=0) break;  /* Done? */
	/* Get more */
	retCode = ObitBDFDataFillBuffer (in, err);
	if (err->error) {
	  Obit_traceback_val (err, routine, in->name, retCode);
	}
	/* If EOF and not done - bail */
	if ((retCode==OBIT_IO_EOF) && (nleft>0)) return retCode;
	lstart = in->current;
	nhere = (in->nBytesInBuffer - (olong)(lstart-in->buffer))/4;
      }  /* end loop over buffers */
    } /* end multiple buffers */
    
    /* If in last segment of buffer, update */
    if ((olong)(in->current-in->buffer)>(BDFBUFFERSIZE*(BDFBUFFERFRAMES-1))) {
      retCode = ObitBDFDataFillBuffer (in, err);
      if (err->error) {
	Obit_traceback_val (err, routine, in->name, retCode);
      }
    }
    /* end 32 bit integer */
    /* 32 bit floats */
  } else if (Dtype==BDFDataType_FLOAT32_TYPE) {
    /* How many in current buffer load */
    nhere = (in->nBytesInBuffer - (olong)(start-in->buffer))/sizeof(ofloat);
    nleft = n;
    out = target;
    lstart = start;
    /* All in buffer? */
    if (nhere>=n) {  /* All in current */
      for (i=0; i<n; i++) out[i] = (olong)lstart[i];
      in->current = lstart + n;
    } else {         /* Multiple buffers */
      while (nleft>0) {
	/* Copy what's here */
	ncopy = MIN (nleft, nhere);
	ncopyb = ncopy*4;  /* in bytes */
	for (i=0; i<ncopy; i++) out[i] = (olong)lstart[i];
	out += ncopy;
	shit = (olong)(out-target);
	in->current = lstart + ncopyb;
	nleft -= ncopy;
	if (nleft<=0) break;  /* Done? */
	/* Get more */
	retCode = ObitBDFDataFillBuffer (in, err);
	if (err->error) {
	  Obit_traceback_val (err, routine, in->name, retCode);
	}
	/* If EOF and not done - bail */
	if ((retCode==OBIT_IO_EOF) && (nleft>0)) return retCode;
	lstart = in->current;
	nhere = (in->nBytesInBuffer - (olong)(lstart-in->buffer))/sizeof(ofloat);
      }  /* end loop over buffers */
    } /* end multiple buffers */
    
    /* If in last segment of buffer, update */
    if ((olong)(in->current-in->buffer)>(BDFBUFFERSIZE*(BDFBUFFERFRAMES-1))) {
      retCode = ObitBDFDataFillBuffer (in, err);
      if (err->error) {
	Obit_traceback_val (err, routine, in->name, retCode);
      }
    }
    /* byte flip if necessary */
    if (byteFlip) {
      out = target;
      for (i=0; i<n; i++) {
	inu.full = out[i];
	outu.parts[0] = inu.parts[3]; 
	outu.parts[1] = inu.parts[2]; 
	outu.parts[2] = inu.parts[1]; 
	outu.parts[3] = inu.parts[0]; 
	out[i] = outu.full;
      }
    } /* end byte flip */

    /* end 32 bit float */
  } else {  /* Unsupported type */
    Obit_log_error(err, OBIT_Error, 
		   "%s: Unsupported data type %d", routine, Dtype);
    return OBIT_IO_SpecErr;
  }

  return retCode;
} /* end CopyLongs */

 /**
 * Find beginning of next Mime segment
 * May update buffer contents.
 * \param in    The object to update
 * \param last  Last byte in buffer of previous segment
 * \param start First byte in buffer of segment
 * \param err   Obit error stack object.
 */
static ObitBDFMIMEType GetNextMIME(ObitBDFData *in, 
				   gchar *last, gchar **start, 
				   ObitErr *err)
{
  ObitBDFMIMEType out = BDFMIMEType_Unknown;
  ObitIOCode retCode;
  gchar *prior, *tstr;
  olong maxStr;
  gboolean isData, isHeader, warn;
  gchar *dataType  = "\nContent-Type: application/octet-stream\n";
  gchar *headerType = "\nContent-Type: text/xml; charset=utf-8\n";
  gchar *routine = "GetNextMIME";

  /* error checks */
  if (err->error) return out;

  /* Look for MIME_boundary-2 */
  maxStr    = in->nBytesInBuffer - (olong)(last-in->buffer);
  prior = "--MIME_boundary-2";
  /*tstr = g_strstr_len (last, maxStr, prior);*/
  tstr = findString (last, maxStr, prior);
  /* Possibly junk in file */
  warn = TRUE;
  while (tstr==NULL) {
    /* Tell about it */
    if (warn) {
     Obit_log_error(err, OBIT_InfoWarn, 
		   "Confused or trashed file");
     warn = FALSE;
    }
    in->current = in->buffer+(BDFBUFFERSIZE-1)*BDFBUFFERFRAMES; /* to end */
    retCode = ObitBDFDataFillBuffer (in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
    /* Skip this buffer load */
    last = in->current;
    maxStr = in->nBytesInBuffer - (olong)(last-in->buffer);
    /* Has the entire file been read and no start info - then something is wrong */
    if (retCode==OBIT_IO_EOF) {
      Obit_log_error(err, OBIT_InfoWarn, 
		     "Hit EOF before next MIME block");
      return out;
    }
    /*tstr = g_strstr_len (last, maxStr, prior);*/
    tstr = findString (last, maxStr, prior);
  }
  *start = tstr + strlen(prior);

  /* if found - see what type */
  if (tstr!=NULL) {
    /* Data? */
    isData = !strncmp(*start, dataType, strlen(dataType));
    /* Header? */
    isHeader = !strncmp(*start, headerType, strlen(headerType));
    out = GetNextMIMEType (in, tstr, start, err);
    /* If it's incomplete, update buffer and try again */
    if (out==BDFMIMEType_Unknown) {
      retCode = ObitBDFDataFillBuffer (in, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, out);
      last = in->current;
      out = GetNextMIMEType (in, last, start, err);
    }
  } else { /* Not found, update buffer */
    retCode = ObitBDFDataFillBuffer (in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
    /* Try again */
    last = in->buffer;
    out = GetNextMIMEType (in, last, start, err);
 }
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  return out;
} /* end GetNextMIME */

 /**
 * Find type of next Mime segment
 * May update buffer contents.
 * \param in    The object to update
 * \param last  Last byte in buffer of previous segment
 * \param start First byte in buffer of segment
 * \param err   Obit error stack object.
 */
static ObitBDFMIMEType GetNextMIMEType(ObitBDFData *in, 
				     gchar *last, gchar **start, 
				     ObitErr *err)
{
  ObitBDFMIMEType out = BDFMIMEType_Unknown;
  gchar *s, *e, *slash, *dot, *Xpad, *prior, *tstr;
  olong i, maxStr;
  /*gchar *routine = "GetNextMIMEType";*/

  /* error checks */
  if (err->error) return out;

  /* Look for boundaries */
  maxStr    = in->nBytesInBuffer - (olong)(last-in->buffer);
  prior = "Content-Location:";
  /*tstr = g_strstr_len (last, maxStr, prior);*/
  tstr = findString (last, maxStr, prior);
  if (tstr==NULL) return out;  /* Cannot find? */
  s = tstr + strlen(prior);
  slash = s;  /* In case no slashing */
  prior = "X-pad: **";
  /*Xpad = g_strstr_len (last, maxStr, prior);*/
  Xpad = findString (last, maxStr, prior);
  if (Xpad==NULL) {;  /* Cannot find? ALMA doesn't have look for \n\n*/
    prior = "\n\n";
    /*Xpad = g_strstr_len (last, maxStr, prior);*/
    Xpad = findString (last, maxStr, prior);
  }
  if (Xpad==NULL) return out;  /* Still Cannot find? */
  e = Xpad;

  /* find last '/' and '.' */
  maxStr = (olong)(e-s);
  for (i=0; i<maxStr; i++) {
    if (s[i]=='/') slash = &s[i];
    if (s[i]=='.') dot   = &s[i];
  }

  /* What type? */
  slash++;  /* Skip slash or blank */
  if (!strncmp(slash, "sdmDataHeader.xml", 17))        out = BDFMIMEType_sdmDataHeader;
  else if (!strncmp(slash, "desc.xml", 8))             out = BDFMIMEType_desc;
  else if (!strncmp(slash, "crossData.bin", 13))       out = BDFMIMEType_crossData;
  else if (!strncmp(slash, "autoData.bin", 12))        out = BDFMIMEType_autoData;
  else if (!strncmp(slash, "flags.bin", 9))            out = BDFMIMEType_flags;
  else if (!strncmp(slash, "actualTimes.bin", 15))     out = BDFMIMEType_actualTimes;
  else if (!strncmp(slash, "actualDurations.bin", 19)) out = BDFMIMEType_actualDurations;
  else if (!strncmp(slash, "weights.bin", 11))         out = BDFMIMEType_weights;
  else if (!strncmp(slash, "zeroLags.bin", 12))        out = BDFMIMEType_zeroLags;
  else out =  BDFMIMEType_Unknown;

  /* Find start - after **\n\n */
  maxStr    = in->nBytesInBuffer - (olong)(Xpad-in->buffer);
  prior = "**\n\n";
  /*tstr = g_strstr_len (Xpad, maxStr, prior);*/
  tstr = findString (Xpad, maxStr, prior);
  /* If not found assume it's ALMA data and use the "\n\n" */
  if (tstr==NULL)  {prior = "\n\n"; tstr = Xpad;}
  *start = tstr + strlen(prior);  /* First byte of data */

  return out;
} /* end GetNextMIMEType */

/**
 * Squeeze all blanks out of a string
 * \param String to squeeze
 */
static void Strip (gchar* s)
{
  olong n, i;

  if (s==NULL) return;  /* Uh oh */
  n = strlen(s);

  /* Leading blanks */
  while (s[0]==' ') {
    for (i=0; i<n-1; i++) s[i] = s[i+1];
    n--;
  }
  /* Trailing blanks */
   for (i=n-1; i>0; i--) {
    if (s[i]==' ') s[i] = 0;
    if (s[i]!=0) break;
  }
 
} /* end Strip */

/**
 * Determine the order of an axis
 * \param axes  Array of axis type enums
 * \param axis  Axis type enum to search for
 * \return 0-rel order, -1 if not found
 */
static olong GetAxisOrder (ObitBDFAxisName *axes, ObitBDFAxisName axis)
{
  olong out = -1, i;

  /* undefined types */
  if ((axis==BDFAxisName_END) || (axis==BDFAxisName_UNK)) return out;
  for (i=0; i<15; i++) {
    if ((axes[i]==BDFAxisName_END) || (axes[i]==BDFAxisName_UNK)) break;
    if (axes[i]==axis) {out = i; break;}
  }
  return out;
} /* end GetAxisOrder */

/**
 * Find first occurance of a string in an array of gchar
 * \param start  beginning of array, may include NULLs
 * \param maxStr maximum length of search
 * \param match  Null terminated string to match, at least 2 characters
 * \return start of match or NULL if not found
 */
static gchar* findString (gchar *start, olong maxStr, gchar *match)
{
  gchar *out=NULL, *test, *end;
  gboolean OK;
  olong i, len;

  len  = strlen (match);
  end  = start + maxStr;
  test = start;
  while (test<end) {
    if ((test[0]==match[0]) && (test[1]==match[1])) {
      if ((test+len)>end) return out;   /* End? */
      /* Check if rest match */
      OK = TRUE;
      for (i=2; i<len; i++) {
	OK = OK && (test[i]==match[i]);
	if (!OK) break;
      }
      if (OK) return test;
    } /* end first two match */
    test++;
  }

  return out;
} /* end findString */
