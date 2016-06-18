/* $Id: ObitPointing.c 516 2015-08-09 00:21:40Z bill.cotton $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2016                                               */
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

#include "ObitPointing.h"
#include "ObitFile.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPointing.c
 * ObitPointing class function definitions.
 * This class is derived from the Obit base class.
 * This class accesses data in the ASDM binary Pointing table
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitPointing";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitPointingClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitPointingClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitPointingInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitPointingClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitPointingClassInfoDefFn (gpointer inClass);

/* Pointing binary table parsing routines */
/** Private: Parse quotedstring from XML string  */
static gchar* Pointingparse_quote_str(gchar *string, olong maxChar, 
				      gchar *prior, gchar **next);
/** Private: Parse 1 byte boolean  */
gboolean Pointingparse_bool(gchar *buffer, 
			    gchar *prior, gchar **next);
/** Private: Parse integer  */
olong Pointingparse_int(gchar *buffer, 
			gchar *prior, gchar **next, 
			gboolean byteFlip);
/** Private: Parse integer from string */
olong Pointingparse_intstr(gchar *buffer, gchar *prefix,
			   gchar *prior, gchar **next, 
			   gboolean byteFlip);
/** Private: Parse long long integer  */
long long Pointingparse_llong(gchar *buffer, 
			      gchar *prior, gchar **next, 
			      gboolean byteFlip);
/** Private: Parse float  */
ofloat Pointingparse_flt(gchar *buffer, 
			 gchar *prior, gchar **next, 
			 gboolean byteFlip);
/** Private: Parse double  */
odouble Pointingparse_dbl(gchar *buffer, 
			  gchar *prior, gchar **next, 
			  gboolean byteFlip);
/** Private: Parse string  */
gchar* Pointingparse_str(gchar *buffer, 
			 gchar *prior, gchar **next, 
			 gboolean byteFlip);
/** Private: Parse time */
odouble Pointingparse_time(gchar *buffer, 
			   gchar *prior, gchar **next, 
			   gboolean byteFlip);
/** Private: Parse time interval */
odouble* Pointingparse_timeint(gchar *buffer, 
			       gchar *prior, gchar **next, 
			       gboolean byteFlip);

/** Parse 1D float  array */
ofloat* Pointingparse_flt_arr(gchar *buffer, 
			      gchar *prior, gchar **next, 
			      gboolean byteFlip);

/** Parse 1D double  array */
odouble* Pointingparse_dbl_arr(gchar *buffer, 
			       gchar *prior, gchar **next, 
			       gboolean byteFlip);

/** Private: Find string in buffer */
static gchar* PointingFindString(gchar *buffer, olong maxChar, 
				 const gchar *string);

/** Private: Look up endian enum */
static ObitPointingEndian LookupEndian(gchar *name);

/* Get start of next MIME segment */
static void GetNextMIME(ObitPointing *in, 
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
ObitPointing* newObitPointing (gchar* name)
{
  ObitPointing* out;
  
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPointingClassInit();
  
  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitPointing));
  
  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");
  
  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;
  
  /* initialize other stuff */
  ObitPointingInit((gpointer)out);
  
 return out;
} /* end newObitPointing */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitPointingGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPointingClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitPointingGetClass */

/**
 * Make a deep copy of an ObitPointing. NYI
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitPointing* ObitPointingCopy  (ObitPointing *in, ObitPointing *out, ObitErr *err)
{
  /*const ObitClassInfo *ParentClass;*/
  /*gboolean oldExist;*/
  /*gchar *outName;*/

  /* error checks */
  if (err->error) return out;

  /* Stubbed */
  g_error("ObitPointingCopy: Stubbed");

  return out;
} /* end ObitPointingCopy */

 /**
 * Make a copy of a object but do not copy the actual data NYI
 * This is useful to create an Pointing similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitPointingClone  (ObitPointing *in, ObitPointing *out, ObitErr *err)
{
  /*const ObitClassInfo *ParentClass;*/

  /* error checks */
  if (err->error) return;

  /* Stubbed */
  g_error("ObitPointingClone: Stubbed");

} /* end ObitPointingClone */

/**
 * Creates an ObitPointing 
 * Parses the ASMD XML tables and stores
 * \param name         An optional name for the object.
 * \param err          Obit error stack object.
 * \return the new object.
 */
ObitPointing* ObitPointingCreate (gchar* name, ObitErr *err)
{
  ObitPointing* out;
  /*gchar *routine="ObitPointingCreate";*/

  /* Create basic structure */
  out = newObitPointing (name);

  /* Init buffer */
  out->buffer = g_malloc0(POINTINGBUFFERSIZE*POINTINGBUFFERFRAMES);
  out->nBytesInBuffer = 0;
  out->current = out->buffer;

  return out;
} /* end ObitPointingCreate */

 /**
 * Initialize File
 * Initializes buffer, parses XML header
 * \param in       The object to fill
 * \param DataFile Name of Mime file with data
 * \param err      Obit error stack object.
 * \param err Obit error stack object.
 */
void ObitPointingInitFile  (ObitPointing *in, gchar *DataFile, ObitErr *err)
{
  ObitIOCode retCode;
  olong i, io, maxStr, itemp;
  gchar *startInfo, *endInfo, *next, *start, *prior, *tstr;
  gchar *ord[20], nelem=18;
  /* gchar *tableUID=NULL, *containerUID=NULL; */
  ObitPointingEndian endian=PointingEndian_Little;
  gchar *routine = "ObitPointingInitFile";

  /* error checks */
  if (err->error) return;

   /* set Values - file name */
  in->DataFile = strdup(DataFile);

  /* Get size */
  in->fileSize = ObitFileSize(DataFile, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create file object */
  if (in->file==NULL) in->file = newObitFile ("Pointing File");

  /* Open */
  ObitFileOpen (in->file, in->DataFile, OBIT_IO_ReadOnly, OBIT_IO_Binary, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Fill Buffer */
  retCode = ObitPointingFillBuffer (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Parse xml header - first find limits */
  maxStr    = in->nBytesInBuffer - (in->current-in->buffer);
  startInfo = g_strstr_len (in->current, maxStr, "<PointingTable");
  endInfo   = g_strstr_len (in->current, maxStr, "</PointingTable>");
  maxStr    = (olong)(endInfo-startInfo);

  /* Find anything? */
  if (startInfo==NULL) return;
  
  /* Endian */
  prior = "byteOrder=";
  tstr = Pointingparse_quote_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    endian = LookupEndian(tstr);
    g_free(tstr);
  }
  
  /* Is Byte flip needed? */
  in->byteFlip =  ((endian==PointingEndian_Big) && (G_BYTE_ORDER==G_LITTLE_ENDIAN)) ||
    ((endian==PointingEndian_Little) && (G_BYTE_ORDER==G_BIG_ENDIAN));

  /* Orders of entries - use buffer location to order */
  startInfo = g_strstr_len (startInfo, maxStr, "<Attributes>");
  maxStr    = (olong)(endInfo-startInfo);
  io = 0;
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<antennaId/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<timeInterval/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<numSample/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<encoder/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<pointingTracking/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<usePolynomials/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<timeOrigin/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<numTerm/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<pointingDirection/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<target/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<offset/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<pointingModelId/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<overTheTop/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<sourceOffset/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<sourceOffsetReferenceCode/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<sourceOffsetEquinox/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<sampledTimeInterval/>");
  ord[io++]    = g_strstr_len (startInfo, maxStr, "<atmosphericCorrection/>");

  /* antennaId = 0 */
  io = 0; in->ordantennaId = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordantennaId++;
  }

  /* timeInterval = 1 */
  io = 1; in->ordtimeInterval = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordtimeInterval++;
  }

  /* numSample = 2 */
  io = 2; in->ordnumSample = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordnumSample++;
  }

  /* encoder = 3 */
  io = 3; in->ordencoder = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordencoder++;
  }

  /* pointingTracking = 4 */
  io = 4; in->ordpointingTracking = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordpointingTracking++;
  }

  /* usePolynomials = 5 */
  io = 5; in->ordusePolynomials = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordusePolynomials++;
  }

  /* timeOrigin = 6 */
  io = 6; in->ordtimeOrigin = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordtimeOrigin++;
  }

  /* numTerm = 7 */
  io = 7; in->ordnumTerm = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordnumTerm++;
  }

  /* pointingDirection = 8 */
  io = 8; in->ordpointingDirection = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordpointingDirection++;
  }

  /* target = 9 */
  io = 9; in->ordtarget = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordtarget++;
  }

  /* offset = 10 */
  io = 10; in->ordoffset = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordoffset++;
  }

  /* pointingModelId = 11 */
  io = 11; in->ordpointingModelId = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordpointingModelId++;
  }

  /* overTheTop = 12 */
  io = 12; in->ordoverTheTop = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordoverTheTop++;
  }

  /* sourceOffset = 13 */
  io = 13; in->ordsourceOffset = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordsourceOffset++;
  }

  /* sourceOffsetReferenceCode = 14 */
  io = 14; in->ordsourceOffsetReferenceCode = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordsourceOffsetReferenceCode++;
  }

  /* sourceOffsetEquinox = 15 */
  io = 15; in->ordsourceOffsetEquinox = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordsourceOffsetEquinox++;
  }

  /* sampledTimeInterval = 16 */
  io = 16; in->ordsampledTimeInterval = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordsampledTimeInterval++;
  }

  /* atmosphericCorrection = 17 */
  io = 17; in->ordatmosphericCorrection = 0;
  for (i=0; i<nelem; i++) {
    if ((i!=io) && (ord[i]<ord[io])) in->ordatmosphericCorrection++;
  }

  /* Need count of actual elements */
  nelem = 0;
  if (in->ordantennaId>=0) nelem++;
  if (in->ordpointingModelId>0) nelem++;
  if (in->ordtimeInterval>0) nelem++;
  if (in->ordnumSample>0) nelem++;
  if (in->ordencoder>0) nelem++;
  if (in->ordpointingTracking>0) nelem++;
  if (in->ordusePolynomials>0) nelem++;
  if (in->ordtimeOrigin>0) nelem++;
  if (in->ordnumTerm>0) nelem++;
  if (in->ordpointingDirection>0) nelem++;
  if (in->ordtarget>0) nelem++;
  if (in->ordoffset>0) nelem++;
  if (in->ordoverTheTop>0) nelem++;
  if (in->ordsourceOffset>0) nelem++;
  if (in->ordsourceOffsetReferenceCode>0) nelem++;
  if (in->ordsourceOffsetEquinox>0) nelem++;
  if (in->ordsampledTimeInterval>0) nelem++;
  if (in->ordatmosphericCorrection>0) nelem++;
  in->nelem = nelem;
  in->current = endInfo;  /* where in buffer */

  /* Go to data table */
  GetNextMIME (in, endInfo, &start, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Parse stuff at beginning including number of rows 
  tableUID     = Pointingparse_str(in->buffer, start, &next, in->byteFlip);
  start        = next;
  containerUID = Pointingparse_str(in->buffer, start, &next, in->byteFlip);
  start        = next;
  in->nrow     = Pointingparse_int(in->buffer, start, &next, in->byteFlip);
  in->curRow   = 0;
  if (tableUID) g_free(tableUID);
  if (containerUID) g_free(containerUID);*/

  /* What a BIZZARRE format - loop through strings until finding a prior
     "length" <=0 or > 1000 and presume that this is the number of rows */
  while (1) {
    itemp = Pointingparse_int(in->buffer, start, &next, in->byteFlip);
    if ((itemp<=0) || (itemp>1000)) break;
    tstr    = Pointingparse_str(in->buffer, start, &next, in->byteFlip);
    start   = next;
    if (tstr) g_free(tstr);
  }

  in->nrow     = Pointingparse_int(in->buffer, start, &next, in->byteFlip);
  in->current = next;  /* where in buffer */


} /* end ObitPointingInitFile */

 /**
 * Fill Buffer
 * If the buffer is filled to capacity, the bottom frame is copied to
 * to the top of the buffer and the remainder of the buffer filled.
 * Updates in->nBytesInBuffer, in->current.
 * \param in  The object to fill
 * \param err Obit error stack object.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_EOF = EOF.
 */
ObitIOCode ObitPointingFillBuffer (ObitPointing *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  olong i, ncopy, size; 
  gchar *top, *copy;
  gchar *routine = "ObitPointingFillBuffer";

  /* error checks */
  if (err->error) return retCode;

  if ((in->fileSize-in->file->filePos) <= 0) return OBIT_IO_EOF;

   /* Is it already full? */
  if (in->nBytesInBuffer>=POINTINGBUFFERSIZE*POINTINGBUFFERFRAMES) { /* Yes - shift */
    top  = in->buffer;
    copy = in->buffer + (POINTINGBUFFERFRAMES-1) * POINTINGBUFFERSIZE;
    memmove (top, copy, (size_t)POINTINGBUFFERSIZE);
    top = &in->buffer[POINTINGBUFFERSIZE];
    ncopy = (POINTINGBUFFERFRAMES-1);
    in->nBytesInBuffer = POINTINGBUFFERSIZE;
    in->current -= (POINTINGBUFFERFRAMES-1) * POINTINGBUFFERSIZE; /* Current position */
  } else {  /* Nope - fill 'er up */
    top  = in->buffer;
    ncopy = POINTINGBUFFERFRAMES;
    in->nBytesInBuffer = 0;
    in->current = in->buffer; /* Current position */
  }

  /* Do reads */
  for (i=0; i<ncopy; i++) {
    /* No more than what's left */
    size = MIN ((olong)POINTINGBUFFERSIZE, (in->fileSize-in->file->filePos));
    if (size<=0) break;
    retCode = ObitFileRead (in->file, -1, size, top, err);
    /* DEBUG
    fprintf (stderr, "Read size %d pos %lld filesize %lld\n", size, in->file->filePos, in->fileSize); */
    if (err->error) {
      Obit_traceback_val (err, routine, in->name, retCode);
    }
    in->nBytesInBuffer += size;
    if (in->file->filePos>=in->fileSize) break;
    top += POINTINGBUFFERSIZE;
  }
  return retCode;
} /* end ObitPointingFillBuffer */

/**
 * Return number of rows in table
 * Should be called after ObitPointingInitFile
 * \param in      Pointer to table
 * \return number of rows given by header
 */
olong ObitPointingGetNrow (ObitPointing *in)
{
  return in->nrow;
} /* end ObitPointingGetNrow*/

/**
 * Fill in values for next table row
 * \param in      Pointer to object to be read.
 * \param row     Row structure to fill in
 * \param err     ObitErr for reporting errors.
 * \return row number from table, -1 = done
 */
olong ObitPointingGetRow (ObitPointing *in, 
			      ASDMPointingRow *row, ObitErr *err)
{
  olong out = in->curRow+1;
  olong maxStr, elem, nterm = 1, nsamp = 1;
  gchar *start, *next, *done, *temp;
  gboolean anybodyHere;
  ObitIOCode retCode=OBIT_IO_OK;
  gchar *routine = "ObitPointingGetRow";
  
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
  for (elem=0; elem<in->nelem; elem++) {
    if (elem==in->ordantennaId) {
      row->antennaId = -1;
      /* antennaId 
	 Find next occurance of prefix */
      temp = g_strstr_len (start+4, in->nBytesInBuffer-(start-in->buffer), "Antenna_");
      if (temp) start = temp-4;
      else {
	temp = PointingFindString(start, in->nBytesInBuffer-(start-in->buffer), "Antenna_");
	if (temp) start = temp-4;
	else {
	  /* Something wrong force read new buffer */
	  start = in->buffer+in->nBytesInBuffer;
	  break;  
	}
      }
      row->antennaId = Pointingparse_intstr(in->buffer, "Antenna_", start, &next, in->byteFlip);
    } else if (elem==in->ordtimeInterval) {              /* timeInterval */
      row->timeInterval = Pointingparse_timeint(in->buffer, start, &next, in->byteFlip);
    } else if (elem==in->ordnumSample) {                 /* numSample */
      row->numSample = Pointingparse_int(in->buffer, start, &next, in->byteFlip);
      nsamp = row->numSample;
    } else if (elem==in->ordencoder) {                   /* encoder - NOT SURE OF TYPE */
      row->encoder = Pointingparse_dbl_arr(in->buffer, start, &next, in->byteFlip);
    } else if (elem==in->ordpointingTracking) {          /* pointingTracking */
      row->pointingTracking = Pointingparse_bool(in->buffer, start, &next);
    } else if (elem==in->ordusePolynomials) {            /* usePolynomials */
      row->usePolynomials = Pointingparse_bool(in->buffer, start, &next);
    } else if (elem==in->ordtimeOrigin) {                /* timeOrigin */
      row->timeOrigin = Pointingparse_time(in->buffer, start, &next, in->byteFlip);
    } else if (elem==in->ordnumTerm ) {                  /* numTerm */
      row->numTerm = Pointingparse_int(in->buffer, start, &next, in->byteFlip);
      nterm = row->numTerm;
    } else if (elem==in->ordpointingDirection) {         /* pointingDirection */
      row->pointingDirection = Pointingparse_dbl_arr(in->buffer, start, &next, in->byteFlip);
    } else if (elem==in->ordtarget) {                    /* target */
      row->target = Pointingparse_dbl_arr(in->buffer, start, &next, in->byteFlip);
    } else if (elem==in->ordoffset) {                    /* offset */
      row->offset = Pointingparse_dbl_arr(in->buffer, start, &next, in->byteFlip);
    } else if (elem==in->ordpointingModelId) {           /* pointingModelId */
      row->PointingModelId  = Pointingparse_int(in->buffer, start, &next, in->byteFlip);
      /* Optional - preceeded by a boolean saying if it is there */
    } else if (elem==in->ordoverTheTop) {                /* overTheTop */
      anybodyHere = Pointingparse_bool(in->buffer, start, &next);
      if (anybodyHere) row->OverTheTop = Pointingparse_bool(in->buffer, start+1, &next);
    } else if (elem==in->ordsourceOffset) {              /* sourceOffset */
      anybodyHere = Pointingparse_bool(in->buffer, start, &next);
      if (anybodyHere) row->SourceOffset = Pointingparse_dbl_arr(in->buffer, start+1, &next, in->byteFlip);
    } else if (elem==in->ordsourceOffsetReferenceCode) {  /* sourceOffsetReferenceCode */
       anybodyHere = Pointingparse_bool(in->buffer, start, &next);
      if (anybodyHere) row->SourceReferenceCodeOffset = Pointingparse_int(in->buffer, start+1, &next, in->byteFlip);
    } else if (elem==in->ordsourceOffsetEquinox) {      /* sourceOffsetEquinox */
      anybodyHere = Pointingparse_bool(in->buffer, start, &next);
      if (anybodyHere) row->SourceEquinoxOffset = Pointingparse_dbl(in->buffer, start+1, &next, in->byteFlip);
    } else if (elem==in->ordsampledTimeInterval) {       /* sampledTimeInterval */
      anybodyHere = Pointingparse_bool(in->buffer, start, &next);
      if (anybodyHere) row->SampledTimeInterval = Pointingparse_timeint(in->buffer, start+1, &next, in->byteFlip);
    } else if (elem==in->ordatmosphericCorrection) {     /*  atmosphericCorrection*/
      anybodyHere = Pointingparse_bool(in->buffer, start, &next);
      if (anybodyHere) row->AtmosphericCorrection = Pointingparse_dbl_arr(in->buffer, start+1, &next, in->byteFlip);
    } else {
                                                         /* Bummer */
    } /* end selection by order */

    start = next;
  } /* end loop over elements */
  in->curRow++;          /* Keep track of row */
  in->current = start;   /* where in buffer */

  /* If out of first segment of buffer read more */
  if (((in->current-in->buffer)>POINTINGBUFFERSIZE) && (in->curRow<in->nrow)) {
    retCode = ObitPointingFillBuffer (in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
  }

  /* Check if read all of file */
  if (retCode==OBIT_IO_EOF) out = -1;
  
  return out;
} /* end ObitPointingGetRow  */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitPointingClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitPointingClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitPointingClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitPointingClassInfoDefFn (gpointer inClass)
{
  ObitPointingClassInfo *theClass = (ObitPointingClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitPointingClassInit;
  theClass->newObit       = (newObitFP)newObitPointing;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitPointingClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitPointingGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitPointingCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitPointingClear;
  theClass->ObitInit      = (ObitInitFP)ObitPointingInit;
  theClass->ObitPointingCreate = (ObitPointingCreateFP)ObitPointingCreate;

} /* end ObitPointingClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitPointingInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPointing *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->DataFile                   = NULL;
  in->file                       = NULL;
  in->nelem                      = 0;
} /* end ObitPointingInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitPointing* cast to an Obit*.
 */
void ObitPointingClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPointing *in = inn;
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
  
} /* end ObitPointingClear */

/* Pointing (Pointing) parsing routines */
/**  Parse double quoted from XLM string 
 * All text from end of prior+1 until next '"'
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 */
static gchar* Pointingparse_quote_str(gchar *string, olong maxChar, 
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
  if ((*b!='"') && (*b!='\'')) return out;  /* Make sure quote */
  b++;                      /* Skip quote */

  /* count */
  charLeft = maxChar - (b-string);
  n = 0;
  for (i=0; i<charLeft; i++) {
    if ((b[i]=='"') || (b[i]=='\'')) break;
    n++;
  }
  out = g_malloc(n+1);
  for (i=0; i<n; i++) out[i] = b[i]; out[i] = 0;
  *next = b + n;

  return out;
} /* end Pointingparse_quote_str */

/**  Parse 8 bit boolean
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \return value
 */
gboolean Pointingparse_bool(gchar *buffer, 
		       gchar *prior, gchar **next)
{
  gboolean out = FALSE;
  gchar *b = prior;

  *next = prior+1;  /* After byte */
  out = b[0] == 1;
    
  return out;
} /* end Pointingparse_bool */

/**  Parse 32 bit integer byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
olong Pointingparse_int(gchar *buffer, 
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
} /* end Pointingparse_int */

/**  Parse 32 bit integer from string
 * \param  buffer   buffer to parse
 * \param  prefix   Prefix for integer
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
olong Pointingparse_intstr(gchar *buffer, gchar *prefix,
		      gchar *prior, gchar **next, 
		      gboolean byteFlip)
{
  olong out = 0;
  olong maxChar;
  gchar *temp, *b, *e;

   /* read string */
  temp = Pointingparse_str(buffer, prior, next, byteFlip);

  maxChar = strlen(temp);
  b = g_strstr_len (temp, maxChar, prefix);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prefix);
  out = (olong)strtol(b, &e, 10);

  if (temp) g_free(temp);
  return out;
} /* end Pointingparse_intstr */

/**  Parse 64 bit integer byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
long long Pointingparse_llong(gchar *buffer, 
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
} /* end Pointingparse_llong */

/**  Parse 32 bit float byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
ofloat Pointingparse_flt(gchar *buffer, 
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
} /* end Pointingparse_flt */

/**  Parse 64 bit float byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
odouble Pointingparse_dbl(gchar *buffer, 
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
} /* end Pointingparse_dbl */

/**  Parse string from buffer, length of string given by preceeding integer.
 * \param  buffer   buffer to parse
 * \param  maxChar  Maximum size of string
 * \param  prior    buffer prior to value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return string pointer, null terminated, should be g_freeed when done
 */
gchar* Pointingparse_str(gchar *buffer, 
		    gchar *prior, gchar **next, 
		    gboolean byteFlip)
{
  gchar *out = NULL;
  gchar *b;
  olong i, len;

  /* How long */
  len = Pointingparse_int(buffer, prior, next, byteFlip);

  /* output string */
  out = g_malloc(len+1);
  b = *next;
  for (i=0; i<len; i++) out[i] = b[i]; out[i] = 0;
  *next += len;

  return out;
} /* end Pointingparse_str */

/**  Parse time interval
 * Read time interval (MJD) in nanoseconds, return as JD
 * At least for ALMA data, the first values is a time and 
 * the second value is a time interval
 * \param  buffer  Bufferto parse
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value as pair of odouble
 */
odouble* Pointingparse_timeint(gchar *buffer, 
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
  temp = Pointingparse_llong(buffer, prior, next, byteFlip);
  out[0] = (odouble)((temp*1.0e-9)/86400.0) + mjdJD0;
  /* Interval */
  b = *next;
  temp = Pointingparse_llong(buffer, b, next, byteFlip);
  out[1] = (odouble)((temp*1.0e-9)/86400.0);

  return out;
} /* end Pointingparse_timeint */

/**  Parse time 
 * Read time  (MJD) in nanoseconds, return as JD
 * \param  buffer  Bufferto parse
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value as pair of odouble
 */
odouble Pointingparse_time(gchar *buffer, 
			   gchar *prior, gchar **next,
			   gboolean byteFlip)
{
  odouble out = 0.0;
  long long temp;
  odouble mjdJD0=2400000.5; /* JD of beginning of MJD time */

  /* Time */
  temp = Pointingparse_llong(buffer, prior, next, byteFlip);
  out = (odouble)((temp*1.0e-9)/86400.0) + mjdJD0;
  return out;
} /* end Pointingparse_time */

/**  Parse 2-D array of 32 bit float byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value array, should be g_freeed when done
 */
ofloat* Pointingparse_flt_arr(gchar *buffer, 
			      gchar *prior, gchar **next, 
			      gboolean byteFlip)
{
  ofloat *out = NULL;
  olong ndim,  num, i;
  olong naxis[10];
  gchar *b;

  /* How many dimensions? Always 2? without saying? 
  ndim = Pointingparse_int(buffer, prior, next, byteFlip);
  b = *next;*/
  num = 1;
  ndim = 2; b = prior;
  /* Size of each dimension */
  for (i=0; i<ndim; i++) {
    naxis[i] = Pointingparse_int(buffer, b, next, byteFlip);
    b = *next;
    num *= MAX(1, naxis[i]);
    }

  /* output */
  out = g_malloc0(num*sizeof(ofloat));

  /* Loop parsing */
  for (i=0; i<num; i++) {
    out[i] = Pointingparse_flt(buffer, b, next, byteFlip);
    b = *next;
  }

  return out;
} /* end Pointingparse_flt_arr */

/**  Parse 2-D array of 64 bit float byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value array, should be g_greeed when done
 */
odouble* Pointingparse_dbl_arr(gchar *buffer, 
			       gchar *prior, gchar **next, 
			       gboolean byteFlip)
{
  odouble *out = NULL;
  olong ndim, num, i;
  olong naxis[10];
  gchar *b;

  /* How many dimensions? Always 2" without saying? 
  ndim = Pointingparse_int(buffer, prior, next, byteFlip);
  b = *next;*/
  ndim = 2; b = prior;
  /* Size of each dimension */
  num = 1;
  for (i=0; i<ndim; i++) {
    naxis[i] = Pointingparse_int(buffer, b, next, byteFlip);
    b = *next;
    num *= MAX(1, naxis[i]);
    }

  /* output */
  out = g_malloc0(num*sizeof(odouble));

  /* Loop parsing */
  for (i=0; i<num; i++) {
    out[i] = Pointingparse_dbl(buffer, b, next, byteFlip);
    b = *next;
  }

  return out;
} /* end Pointingparse_dbl_arr */

/**  Find string in buffer, allows many nulls in buffer
 * \param  buffer   buffer start to search
 * \param  maxChar  Maximum number of characters to search
 * \param  string   String to locate
 * \return pointer to start of string if found, else NULL
 */
static gchar* PointingFindString (gchar *buffer, olong maxChar, const gchar *string)
{
  gchar *out = NULL, *tstr;
  olong nt, nchk, i;

  /* How far to check? */
  nt   = strlen(string);
  nchk = maxChar - nt;

  /* Loop */
  for (i=0; i<nchk; i++) {
    if (buffer[i]==string[0]) {  /* Maybe */
      tstr = g_strstr_len (&buffer[i], nt, string);
      if (tstr) return tstr;  /* Got it */
    }
  } /* end loop over buffer */

  return out;
} /* end PointingFindString */

 /**
 * Find beginning of next Mime segment
 * May update buffer contents.
 * \param in    The object to update
 * \param last  Last byte in buffer of previous segment
 * \param start First byte in buffer of data segment
 * \param err   Obit error stack object.
 */
static void GetNextMIME(ObitPointing *in, 
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
static ObitPointingEndian LookupEndian(gchar *name)
{
  ObitPointingEndian out = 0;
 
  if (!strncmp (name, "Little_Endian", 13)) return PointingEndian_Little;
  if (!strncmp (name, "Big_Endian",    10)) return PointingEndian_Big;
  return out;
} /* end LookupEndian */

