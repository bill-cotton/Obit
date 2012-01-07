/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012                                               */
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

#include "ObitEVLASysPower.h"
#include "ObitFile.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitEVLASysPower.c
 * ObitEVLASysPower class function definitions.
 * This class is derived from the Obit base class.
 * This class accesses data in the EVLA BDF format
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitEVLASysPower";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitEVLASysPowerClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitEVLASysPowerClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitEVLASysPowerInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitEVLASysPowerClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitEVLASysPowerClassInfoDefFn (gpointer inClass);

/* EVLA SysPower (ESP) parsing routines */
/** Private: Parse quotedstring from XML string  */
static gchar* ESPparse_quote_str(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next);
/** Private: Parse 1 byte boolean  */
gboolean ESPparse_bool(gchar *buffer, 
		       gchar *prior, gchar **next);
/** Private: Parse integer  */
olong ESPparse_int(gchar *buffer, 
		   gchar *prior, gchar **next, 
		   gboolean byteFlip);
/** Private: Parse integer from string */
olong ESPparse_intstr(gchar *buffer, gchar *prefix,
		      gchar *prior, gchar **next, 
		      gboolean byteFlip);
/** Private: Parse long long integer  */
long long ESPparse_llong(gchar *buffer, 
			 gchar *prior, gchar **next, 
			 gboolean byteFlip);
/** Private: Parse float  */
ofloat ESPparse_flt(gchar *buffer, 
		    gchar *prior, gchar **next, 
		    gboolean byteFlip);
/** Private: Parse double  */
odouble ESPparse_dbl(gchar *buffer, 
		     gchar *prior, gchar **next, 
		     gboolean byteFlip);
/** Private: Parse string  */
gchar* ESPparse_str(gchar *buffer, 
		    gchar *prior, gchar **next, 
		    gboolean byteFlip);
/** Private: Parse time interval */
odouble* ESPparse_timeint(gchar *buffer, 
			  gchar *prior, gchar **next, 
			  gboolean byteFlip);

/** Parse 1D float  array */
ofloat* ESPparse_flt_arr(gchar *buffer, 
			 gchar *prior, gchar **next, 
			 gboolean byteFlip);

/** Parse 1D double  array */
odouble* ESPparse_dbl_arr(gchar *buffer, 
			  gchar *prior, gchar **next, 
			  gboolean byteFlip);

/** Private: Look up endian enum */
static ObitESPEndian LookupEndian(gchar *name);

/* Get start of next MIME segment */
static void GetNextMIME(ObitEVLASysPower *in, 
			gchar *last, gchar **start, 
			ObitErr *err);

/*----------------- Union definitions ----------------------*/
/** Used for byte swapping 32 bit ints */
 union lequiv { 
   olong full; 
   gchar parts[4];
 }; 

/** Used for byte swapping 64 bit ints */
 union llequiv { 
   long long full; 
   gchar parts[8];
 }; 

/** Used for byte swapping floats */
 union fequiv { 
   ofloat full; 
   gchar parts[4]; 
 }; 

/** Used for byte swapping doubles */
 union dequiv { 
   odouble full; 
   gchar parts[8]; 
 }; 

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitEVLASysPower* newObitEVLASysPower (gchar* name)
{
  ObitEVLASysPower* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitEVLASysPowerClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitEVLASysPower));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitEVLASysPowerInit((gpointer)out);

 return out;
} /* end newObitEVLASysPower */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitEVLASysPowerGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitEVLASysPowerClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitEVLASysPowerGetClass */

/**
 * Make a deep copy of an ObitEVLASysPower. NYI
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitEVLASysPower* ObitEVLASysPowerCopy  (ObitEVLASysPower *in, ObitEVLASysPower *out, ObitErr *err)
{
  /*const ObitClassInfo *ParentClass;*/
  /*gboolean oldExist;*/
  /*gchar *outName;*/

  /* error checks */
  if (err->error) return out;

  /* Stubbed */
  g_error("ObitEVLASysPowerCopy: Stubbed");

  return out;
} /* end ObitEVLASysPowerCopy */

 /**
 * Make a copy of a object but do not copy the actual data NYI
 * This is useful to create an EVLASysPower similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitEVLASysPowerClone  (ObitEVLASysPower *in, ObitEVLASysPower *out, ObitErr *err)
{
  /*const ObitClassInfo *ParentClass;*/

  /* error checks */
  if (err->error) return;

  /* Stubbed */
  g_error("ObitEVLASysPowerClone: Stubbed");

} /* end ObitEVLASysPowerClone */

/**
 * Creates an ObitEVLASysPower 
 * Parses the ASMD XML tables and stores
 * \param name         An optional name for the object.
 * \param err          Obit error stack object.
 * \return the new object.
 */
ObitEVLASysPower* ObitEVLASysPowerCreate (gchar* name, ObitErr *err)
{
  ObitEVLASysPower* out;
  /*gchar *routine="ObitEVLASysPowerCreate";*/

  /* Create basic structure */
  out = newObitEVLASysPower (name);

  /* Init buffer */
  out->buffer = g_malloc0(ESPBUFFERSIZE*ESPBUFFERFRAMES);
  out->nBytesInBuffer = 0;
  out->current = out->buffer;

  return out;
} /* end ObitEVLASysPowerCreate */

 /**
 * Initialize File
 * Initializes buffer, parses XML header
 * \param in       The object to fill
 * \param DataFile Name of Mime file with data
 * \param err      Obit error stack object.
 * \param err Obit error stack object.
 */
void ObitEVLASysPowerInitFile  (ObitEVLASysPower *in, gchar *DataFile, ObitErr *err)
{
  ObitIOCode retCode;
  olong i, count, maxStr, itemp;
  gchar *startInfo, *endInfo, *next, *start, *prior, *tstr;
  gchar *ord[10];
  /* gchar *tableUID=NULL, *containerUID=NULL; */
  ObitESPEndian endian=ESPEndian_Little;
  gchar *routine = "ObitEVLASysPowerInitFile";

  /* error checks */
  if (err->error) return;

   /* set Values - file name */
  in->DataFile = strdup(DataFile);

  /* Get size */
  in->fileSize = ObitFileSize(DataFile, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create file object */
  if (in->file==NULL) in->file = newObitFile ("ESP File");

  /* Open */
  ObitFileOpen (in->file, in->DataFile, OBIT_IO_ReadOnly, OBIT_IO_Binary, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Fill Buffer */
  retCode = ObitEVLASysPowerFillBuffer (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Parse xml header - first find limits */
  maxStr    = in->nBytesInBuffer - (in->current-in->buffer);
  startInfo = g_strstr_len (in->current, maxStr, "<SysPowerTable>");
  endInfo   = g_strstr_len (in->current, maxStr, "</SysPowerTable>");
  maxStr    = (olong)(endInfo-startInfo);
  
  /* Endian */
  prior = "byteOrder=";
  tstr = ESPparse_quote_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    endian = LookupEndian(tstr);
    g_free(tstr);
  }
  
  /* Is Byte flip needed? */
  in->byteFlip =  ((endian==ESPEndian_Big) && (G_BYTE_ORDER==G_LITTLE_ENDIAN)) ||
    ((endian==ESPEndian_Little) && (G_BYTE_ORDER==G_BIG_ENDIAN));

  /* Orders of entries - use buffer location to order */
  startInfo = g_strstr_len (in->current, maxStr, "<Attributes>");
  maxStr    = (olong)(endInfo-startInfo);
  ord[0]    = g_strstr_len (in->current, maxStr, "<antennaId/>");
  ord[1]    = g_strstr_len (in->current, maxStr, "<spectralWindowId/>");
  ord[2]    = g_strstr_len (in->current, maxStr, "<feedId/>");
  ord[3]    = g_strstr_len (in->current, maxStr, "<timeInterval/>");
  ord[4]    = g_strstr_len (in->current, maxStr, "<numReceptor/>");
  ord[5]    = g_strstr_len (in->current, maxStr, "<switchedPowerDifference/>");
  ord[6]    = g_strstr_len (in->current, maxStr, "<switchedPowerSum/>");
  ord[7]    = g_strstr_len (in->current, maxStr, "<requantizerGain/>");

  /* antennaId = 0 */
  count = 0;
  in->ordantennaId = 0;
  for (i=0; i<8; i++) {
    if ((i!=0) && (ord[i]<ord[0])) in->ordantennaId++;
  }

  /*spectralWindowId  = 1 */
  count = 0;
  in->ordspectralWindowId = 0;
  for (i=0; i<8; i++) {
    if ((i!=1) && (ord[i]<ord[1])) in->ordspectralWindowId++;
  }

  /* feedId = 2 */
  count = 0;
  in->ordfeedId = 0;
  for (i=0; i<8; i++) {
    if ((i!=2) && (ord[i]<ord[2])) in->ordfeedId++;
  }

  /* timeInterval = 3 */
  count = 0;
  in->ordtimeInterval = 0;
  for (i=0; i<8; i++) {
    if ((i!=3) && (ord[i]<ord[3])) in->ordtimeInterval++;
  }

  /* numReceptor = 4 */
  count = 0;
  in->ordnumReceptor = 0;
  for (i=0; i<8; i++) {
    if ((i!=4) && (ord[i]<ord[4])) in->ordnumReceptor++;
  }

  /* switchedPowerDifference = 5 */
  count =0;
  in->ordswitchedPowerDifference = 0;
  for (i=0; i<8; i++) {
    if ((i!=5) && (ord[i]<ord[5])) in->ordswitchedPowerDifference++;
  }

  /* switchedPowerSum = 6 */
  count = 0;
  in->ordswitchedPowerSum = 0;
  for (i=0; i<8; i++) {
    if ((i!=6) && (ord[i]<ord[6])) in->ordswitchedPowerSum++;
  }

  /* requantizerGain = 7 */
  count = 0;
  in->ordrequantizerGain = 0;
  for (i=0; i<8; i++) {
    if ((i!=7) && (ord[i]<ord[7])) in->ordrequantizerGain++;
  }
  in->current = endInfo;  /* where in buffer */

  /* Go to data table */
  GetNextMIME (in, endInfo, &start, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Parse stuff at beginning including number of rows 
  tableUID     = ESPparse_str(in->buffer, start, &next, in->byteFlip);
  start        = next;
  containerUID = ESPparse_str(in->buffer, start, &next, in->byteFlip);
  start        = next;
  in->nrow     = ESPparse_int(in->buffer, start, &next, in->byteFlip);
  in->curRow   = 0;
  if (tableUID) g_free(tableUID);
  if (containerUID) g_free(containerUID);*/

  /* What a BIZZARRE format - loop through strings until finding a prior
     "length" <=0 or > 1000 and presume that this is the number of rows */
  while (1) {
    itemp = ESPparse_int(in->buffer, start, &next, in->byteFlip);
    if ((itemp<=0) || (itemp>1000)) break;
    tstr    = ESPparse_str(in->buffer, start, &next, in->byteFlip);
    start   = next;
    if (tstr) g_free(tstr);
  }

  in->nrow     = ESPparse_int(in->buffer, start, &next, in->byteFlip);
  in->current = next;  /* where in buffer */


} /* end ObitEVLASysPowerInitFile */

 /**
 * Fill Buffer
 * If the buffer is filled to capacity, the bottom frame is copied to
 * to the top of the buffer and the remainder of the buffer filled.
 * Updates in->nBytesInBuffer, in->current.
 * \param in  The object to fill
 * \param err Obit error stack object.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_EOF = EOF.
 */
ObitIOCode ObitEVLASysPowerFillBuffer (ObitEVLASysPower *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  olong i, ncopy, size; 
  gchar *top, *copy;
  gchar *routine = "ObitEVLASysPowerFillBuffer";

  /* error checks */
  if (err->error) return retCode;

  if ((in->fileSize-in->file->filePos) <= 0) return OBIT_IO_EOF;

   /* Is it already full? */
  if (in->nBytesInBuffer>=ESPBUFFERSIZE*ESPBUFFERFRAMES) { /* Yes - shift */
    top  = in->buffer;
    copy = in->buffer + (ESPBUFFERFRAMES-1) * ESPBUFFERSIZE;
    memmove (top, copy, (size_t)ESPBUFFERSIZE);
    top = &in->buffer[ESPBUFFERSIZE];
    ncopy = (ESPBUFFERFRAMES-1);
    in->nBytesInBuffer = ESPBUFFERSIZE;
    in->current -= (ESPBUFFERFRAMES-1) * ESPBUFFERSIZE; /* Current position */
  } else {  /* Nope - fill 'er up */
    top  = in->buffer;
    ncopy = ESPBUFFERFRAMES;
    in->nBytesInBuffer = 0;
    in->current = in->buffer; /* Current position */
  }

  /* Do reads */
  for (i=0; i<ncopy; i++) {
    /* No more than what's left */
    size = MIN ((olong)ESPBUFFERSIZE, (in->fileSize-in->file->filePos));
    if (size<=0) break;
    retCode = ObitFileRead (in->file, -1, size, top, err);
    /* DEBUG
    fprintf (stderr, "Read size %d pos %lld filesize %lld\n", size, in->file->filePos, in->fileSize); */
    if (err->error) {
      Obit_traceback_val (err, routine, in->name, retCode);
    }
    in->nBytesInBuffer += size;
    if (in->file->filePos>=in->fileSize) break;
    top += ESPBUFFERSIZE;
  }
  return retCode;
} /* end ObitEVLASysPowerFillBuffer */

/**
 * Return number of rows in table
 * Should be called after ObitEVLASysPowerInitFile
 * \param in      Pointer to table
 * \return number of rows given by header
 */
olong ObitEVLASysPowerGetNrow (ObitEVLASysPower *in)
{
  return in->nrow;
} /* end ObitEVLASysPowerGetNrow*/

/**
 * Fill in values for next table row
 * \param in      Pointer to object to be read.
 * \param row     Row structure to fill in
 * \param err     ObitErr for reporting errors.
 * \return row number from table, -1 = done
 */
olong ObitEVLASysPowerGetRow (ObitEVLASysPower *in, 
			      ASDMSysPowerRow *row, ObitErr *err)
{
  olong out = in->curRow+1;
  olong maxStr, elem, nelem=8;
  gchar *start, *next, *done;
  gboolean gotIt;
  odouble mjdJD0=2400000.5; /* JD of beginning of MJD time */
  gchar *routine = "ObitEVLASysPowerGetRow";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (row != NULL);

  /* Done? */
  if (out>in->nrow) return -1;

  /* Look for --MIME_boundary-- marking end of data */
  maxStr = 30;
  done = g_strstr_len (in->current, maxStr, "--MIME_boundary--");
  if (done!=NULL) return -1;

  start = in->current;
  /* Loop over elements in row */
  for (elem=0; elem<nelem; elem++) {
    if (elem==in->ordantennaId) {
      /* antennaId */
      row->antennaId = ESPparse_intstr(in->buffer, "Antenna_", start, &next, in->byteFlip);
    } else if (elem==in->ordspectralWindowId) {
      /* spectralWindowId */  
      row->spectralWindowId = ESPparse_intstr(in->buffer, "SpectralWindow_", 
					      start, &next, in->byteFlip);
    } else if (elem==in->ordfeedId) {
      /* feedId */
      row->feedId = ESPparse_int(in->buffer, start, &next, in->byteFlip);
    } else if (elem==in->ordtimeInterval) {
      /* timeInterval */
      row->timeInterval = ESPparse_timeint(in->buffer, start, &next, in->byteFlip);
      /* Remove JD0 offset from second */
      if ((row->timeInterval[1]<row->timeInterval[0]) && (row->timeInterval[1]>mjdJD0))
	row->timeInterval[1] -= mjdJD0;
    } else if (elem==in->ordnumReceptor) {
      /* numReceptor */
      row->numReceptor = ESPparse_int(in->buffer, start, &next, in->byteFlip);
    } else if (elem==in->ordswitchedPowerDifference) {
      /* switchedPowerDifference */
      gotIt = ESPparse_bool(in->buffer, start, &next);
      if (gotIt) {
	start = next;
	row->switchedPowerDifference = ESPparse_flt_arr(in->buffer, start, &next, in->byteFlip);
      }
    } else if (elem==in->ordswitchedPowerSum) {
      /* switchedPowerSum */
      gotIt = ESPparse_bool(in->buffer, start, &next);
      if (gotIt) {
	start = next;
	row->switchedPowerSum = ESPparse_flt_arr(in->buffer, start, &next, in->byteFlip);
      }
    } else if (elem==in->ordrequantizerGain) {
      /* requantizerGain */
      gotIt = ESPparse_bool(in->buffer, start, &next);
      if (gotIt) {
	start = next;
	row->requantizerGain = ESPparse_flt_arr(in->buffer, start, &next, in->byteFlip);
      }
    } else {
      /* Bummer */
    }
    start = next;
  } /* end loop over elements */
  in->curRow++;          /* Keep track of row */
  in->current = start;   /* where in buffer */

  /* If out of first segment of buffer read more */
  if (((in->current-in->buffer)>ESPBUFFERSIZE) && (in->curRow<in->nrow)) {
    ObitEVLASysPowerFillBuffer (in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
  }
  
  return out;
} /* end ObitEVLASysPowerGetRow  */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitEVLASysPowerClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitEVLASysPowerClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitEVLASysPowerClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitEVLASysPowerClassInfoDefFn (gpointer inClass)
{
  ObitEVLASysPowerClassInfo *theClass = (ObitEVLASysPowerClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitEVLASysPowerClassInit;
  theClass->newObit       = (newObitFP)newObitEVLASysPower;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitEVLASysPowerClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitEVLASysPowerGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitEVLASysPowerCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitEVLASysPowerClear;
  theClass->ObitInit      = (ObitInitFP)ObitEVLASysPowerInit;
  theClass->ObitEVLASysPowerCreate = (ObitEVLASysPowerCreateFP)ObitEVLASysPowerCreate;

} /* end ObitEVLASysPowerClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitEVLASysPowerInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitEVLASysPower *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->DataFile                   = NULL;
  in->buffer                     = NULL;
  in->ordantennaId               = 0;
  in->ordspectralWindowId        = 1;
  in->ordfeedId                  = 2;
  in->ordtimeInterval            = 3;
  in->ordnumReceptor             = 4;
  in->ordswitchedPowerDifference = 5;
  in->ordswitchedPowerSum        = 6;
  in->ordrequantizerGain         = 7;
} /* end ObitEVLASysPowerInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitEVLASysPower* cast to an Obit*.
 */
void ObitEVLASysPowerClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitEVLASysPower *in = inn;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Close file */
  err = newObitErr();
  if (in->file) ObitFileClose (in->file, err);
  err = ObitErrUnref(err);

  /* delete this class members */
  in->file      = ObitFileUnref(in->file);
  if (in->DataFile)            g_free(in->DataFile);
  if (in->buffer)              g_free(in->buffer);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitEVLASysPowerClear */

/* EVLASysPower (ESP) parsing routines */
/**  Parse double quoted from XLM string 
 * All text from end of prior+1 until next '"'
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 */
static gchar* ESPparse_quote_str(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next)
{
  gchar *out = NULL;
  gchar *b;
  olong charLeft;
  olong i, n;

  *next = string;  /* if not found */
  b = g_strstr_len (string, maxChar, prior);
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
} /* end ESPparse_quote_str */

/**  Parse 8 bit boolean
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \return value
 */
gboolean ESPparse_bool(gchar *buffer, 
		       gchar *prior, gchar **next)
{
  gboolean out = FALSE;
  gchar *b = prior;

  *next = prior+1;  /* After byte */
  out = b[0] == 1;
    
  return out;
} /* end ESPparse_bool */

/**  Parse 32 bit integer byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
olong ESPparse_int(gchar *buffer, 
		   gchar *prior, gchar **next, 
		   gboolean byteFlip)
{
  olong out = 0;
  union lequiv equiv;
  gchar *b = prior;

  *next = prior+4;  /* After bytes */
  if (byteFlip) {
    equiv.parts[0] = b[3];
    equiv.parts[1] = b[2];
    equiv.parts[2] = b[1];
    equiv.parts[3] = b[0];
  } else {
    equiv.parts[0] = b[0];
    equiv.parts[1] = b[1];
    equiv.parts[2] = b[2];
    equiv.parts[3] = b[3];
  }
  out = equiv.full;
    
  return out;
} /* end ESPparse_int */

/**  Parse 32 bit integer from string
 * \param  buffer   buffer to parse
 * \param  prefix   Prefix for integer
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
olong ESPparse_intstr(gchar *buffer, gchar *prefix,
		      gchar *prior, gchar **next, 
		      gboolean byteFlip)
{
  olong out = 0;
  olong maxChar;
  gchar *temp, *b, *e;

  /* read string */
  temp = ESPparse_str(buffer, prior, next, byteFlip);

  maxChar = strlen(temp);
  b = g_strstr_len (temp, maxChar, prefix);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prefix);
  out = (olong)strtol(b, &e, 10);

  if (temp) g_free(temp);
  return out;
} /* end ESPparse_intstr */

/**  Parse 64 bit integer byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
long long ESPparse_llong(gchar *buffer, 
			 char *prior, gchar **next, 
			 gboolean byteFlip)
{
  long long out = 0;
  union llequiv equiv;
  gchar *b = prior;

  *next = prior+8;  /* After bytes */
  if (byteFlip) {
    equiv.parts[0] = b[7];
    equiv.parts[1] = b[6];
    equiv.parts[2] = b[5];
    equiv.parts[3] = b[4];
    equiv.parts[4] = b[3];
    equiv.parts[5] = b[2];
    equiv.parts[6] = b[1];
    equiv.parts[7] = b[0];
  } else {
    equiv.parts[0] = b[0];
    equiv.parts[1] = b[1];
    equiv.parts[2] = b[2];
    equiv.parts[3] = b[3];
    equiv.parts[4] = b[4];
    equiv.parts[5] = b[5];
    equiv.parts[6] = b[6];
    equiv.parts[7] = b[7];
  }
  out = equiv.full;
    
  return out;
} /* end ESPparse_llong */

/**  Parse 32 bit float byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
ofloat ESPparse_flt(gchar *buffer, 
		    gchar *prior, gchar **next, 
		    gboolean byteFlip)
{
  ofloat out = 0.0;
  union fequiv equiv;
  gchar *b = prior;

  *next = prior+4;  /* After bytes */
  if (byteFlip) {
    equiv.parts[0] = b[3];
    equiv.parts[1] = b[2];
    equiv.parts[2] = b[1];
    equiv.parts[3] = b[0];
  } else {
    equiv.parts[0] = b[0];
    equiv.parts[1] = b[1];
    equiv.parts[2] = b[2];
    equiv.parts[3] = b[3];
  }
  out = equiv.full;
    
  return out;
} /* end ESPparse_flt */

/**  Parse 64 bit float byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
odouble ESPparse_dbl(gchar *buffer, 
		     char *prior, gchar **next, 
		     gboolean byteFlip)
{
  odouble out = 0.0;
  union dequiv equiv;
  gchar *b = prior;

  *next = prior+8;  /* After bytes */
  if (byteFlip) {
    equiv.parts[0] = b[7];
    equiv.parts[1] = b[6];
    equiv.parts[2] = b[5];
    equiv.parts[3] = b[4];
    equiv.parts[4] = b[3];
    equiv.parts[5] = b[2];
    equiv.parts[6] = b[1];
    equiv.parts[7] = b[0];
  } else {
    equiv.parts[0] = b[0];
    equiv.parts[1] = b[1];
    equiv.parts[2] = b[2];
    equiv.parts[3] = b[3];
    equiv.parts[4] = b[4];
    equiv.parts[5] = b[5];
    equiv.parts[6] = b[6];
    equiv.parts[7] = b[7];
  }
  out = equiv.full;
    
  return out;
} /* end ESPparse_dbl */

/**  Parse string from buffer, length of string given by preceeding integer.
 * \param  buffer   buffer to parse
 * \param  maxChar  Maximum size of string
 * \param  prior    buffer prior to value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return string pointer, null terminated, should be g_freeed when done
 */
gchar* ESPparse_str(gchar *buffer, 
		    gchar *prior, gchar **next, 
		    gboolean byteFlip)
{
  gchar *out = NULL;
  gchar *b;
  olong i, len;

  /* How long */
  len = ESPparse_int(buffer, prior, next, byteFlip);

  /* output string */
  out = g_malloc(len+1);
  b = *next;
  for (i=0; i<len; i++) out[i] = b[i]; out[i] = 0;
  *next += len;

  return out;
} /* end ESPparse_str */

/**  Parse time interval
 * Read time interval (MJD) in nanoseconds, return as JD
 * \param  buffer  Bufferto parse
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value as pair of odouble
 */
odouble* ESPparse_timeint(gchar *buffer, 
			  gchar *prior, gchar **next,
			  gboolean byteFlip)
{
  odouble *out = NULL;
  long long temp;
  gchar *b;
  odouble mjdJD0=2400000.5; /* JD of beginning of MJD time */

  /* output */
  out = g_malloc0(2*sizeof(odouble));
  /* Time */
  temp = ESPparse_llong(buffer, prior, next, byteFlip);
  out[0] = (odouble)((temp*1.0e-9)/86400.0) + mjdJD0;
  /* Interval */
  b = *next;
  temp = ESPparse_llong(buffer, b, next, byteFlip);
  out[1] = (odouble)((temp*1.0e-9)/86400.0) + mjdJD0;

  return out;
} /* end ESPparse_timeint */

/**  Parse 1-D array of 32 bit float byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value array, should be g_freeed when done
 */
ofloat* ESPparse_flt_arr(gchar *buffer, 
			 gchar *prior, gchar **next, 
			 gboolean byteFlip)
{
  ofloat *out = NULL;
  olong ndim,  num, i;
  /* olong naxis[10]; */
  gchar *b;

  /* How many dimensions? */
  ndim = ESPparse_int(buffer, prior, next, byteFlip);
  b = *next;
  /* Size of each dimension
     NO - implementation NOT what documented 
  num = 1;
  for (i=0; i<ndim; i++) {
    naxis[i] = ESPparse_int(buffer, b, next, byteFlip);
    b = *next;
    num *= MAX(1, naxis[i]);
    }*/
  num = ndim;

  /* output */
  out = g_malloc0(num*sizeof(ofloat));

  /* Loop parsing */
  for (i=0; i<num; i++) {
    out[i] = ESPparse_flt(buffer, b, next, byteFlip);
    b = *next;
  }

  return out;
} /* end ESPparse_flt_arr */

/**  Parse 1-D array of 32 bit float byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value array, should be g_greeed when done
 */
odouble* ESPparse_dbl_arr(gchar *buffer, 
			  gchar *prior, gchar **next, 
			  gboolean byteFlip)
{
  odouble *out = NULL;
  olong ndim, num, i;
  /* olong naxis[10]; */
  gchar *b;

  /* How many dimensions? */
  ndim = ESPparse_int(buffer, prior, next, byteFlip);
  b = *next;
  /* Size of each dimension
     NO - implementation NOT what documented 
  num = 1;
  for (i=0; i<ndim; i++) {
    naxis[i] = ESPparse_int(buffer, b, next, byteFlip);
    b = *next;
    num *= MAX(1, naxis[i]);
    }*/
  num = ndim;

  /* output */
  out = g_malloc0(num*sizeof(odouble));

  /* Loop parsing */
  for (i=0; i<num; i++) {
    out[i] = ESPparse_dbl(buffer, b, next, byteFlip);
    b = *next;
  }

  return out;
} /* end ESPparse_dbl_arr */

 /**
 * Find beginning of next Mime segment
 * May update buffer contents.
 * \param in    The object to update
 * \param last  Last byte in buffer of previous segment
 * \param start First byte in buffer of data segment
 * \param err   Obit error stack object.
 */
static void GetNextMIME(ObitEVLASysPower *in, 
			gchar *last, gchar **start, 
			ObitErr *err)
{
  gchar *prior, *tstr;
  olong maxStr;
  gboolean isData;
  gchar *dataType  = "\nContent-Type: binary/octet-stream\n";
  gchar *routine = "GetNextMIME";

  /* error checks */
  if (err->error) return;

  /* Look for MIME_boundary */
  maxStr    = in->nBytesInBuffer - (olong)(last-in->buffer);
  prior = "--MIME_boundary";
  tstr = g_strstr_len (last, maxStr, prior);
  *start = tstr + strlen(prior);

  /* if found - see what type */
  if (tstr!=NULL) {
    /* Data? */
    isData = !strncmp(*start, dataType, strlen(dataType));
  } else { /* Not found, something is wrong */
    isData = FALSE;
 }
  /* Better be data */
  if (!isData) {
    Obit_log_error(err, OBIT_Error, "Could not find data table");
    return;
  }
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Find "\n\n" */
  prior = "\n\n";
  tstr = g_strstr_len (last, maxStr, prior);

  /* Space to first non blank */
  *start = tstr + strlen(prior);
  while ((*(*start))==' ') 
    {(*start)++;}

  return;
} /* end GetNextMIME */

/**  Look up data type enum 
 * \param  name  String to lookup
 * \return value
 */
static ObitESPEndian LookupEndian(gchar *name)
{
  ObitESPEndian out = 0;
 
  if (!strncmp (name, "Little_Endian", 13)) return ESPEndian_Little;
  if (!strncmp (name, "Big_Endian",    10)) return ESPEndian_Big;
  return out;
} /* end LookupEndian */

