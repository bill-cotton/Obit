/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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

#include "ObitALMACalAtm.h"
#include "ObitFile.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitALMACalAtm.c
 * ObitALMACalAtm class function definitions.
 * This class is derived from the Obit base class.
 * This class accesses data in the ALMA binary calAtmosphere files
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitALMACalAtm";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitALMACalAtmClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitALMACalAtmClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitALMACalAtmInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitALMACalAtmClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitALMACalAtmClassInfoDefFn (gpointer inClass);

/* ALMA CalAtmosphere (ACA) parsing routines */
/** Private: Parse quoted string  */
static gchar* ACAparse_quote_str(gchar *string, olong maxChar, 
				 gchar *prior, gchar **next);
/** Private: Parse 1 byte boolean  */
gboolean ACAparse_bool(gchar *buffer, 
		       gchar *prior, gchar **next);
/** Private: Parse integer  */
olong ACAparse_int(gchar *buffer, 
		   gchar *prior, gchar **next, 
		   gboolean byteFlip);
/** Private: Parse integer from string */
olong ACAparse_intstr(gchar *buffer, gchar *prefix,
		      gchar *prior, gchar **next, 
		      gboolean byteFlip);
/** Private: Parse long long integer  */
long long ACAparse_llong(gchar *buffer, 
			 gchar *prior, gchar **next, 
			 gboolean byteFlip);
/** Private: Parse float  */
ofloat ACAparse_flt(gchar *buffer, 
		    gchar *prior, gchar **next, 
		    gboolean byteFlip);
/** Private: Parse double  */
odouble ACAparse_dbl(gchar *buffer, 
		     gchar *prior, gchar **next, 
		     gboolean byteFlip);
/** Private: Parse string  */
gchar* ACAparse_str(gchar *buffer, 
		    gchar *prior, gchar **next, 
		    gboolean byteFlip);
/** Private: Parse poln string  */
gchar* ACAparse_poln(gchar *buffer, 
		    gchar *prior, gchar **next, 
		    gboolean byteFlip);
/** Private: Parse time  */
odouble ACAparse_time(gchar *buffer, 
			  gchar *prior, gchar **next, 
			  gboolean byteFlip);

/** Private: Parse time interval */
odouble* ACAparse_timeint(gchar *buffer, 
			  gchar *prior, gchar **next, 
			  gboolean byteFlip);

/** Parse 1D float  array */
ofloat* ACAparse_flt_arr(gchar *buffer, 
			 gchar *prior, gchar **next, olong ncrap, olong nwords,
			 gboolean byteFlip);

/** Parse 1D double  array */
odouble* ACAparse_dbl_arr(gchar *buffer, 
			  gchar *prior, gchar **next, olong ncrap, olong nwords,
			  gboolean byteFlip);

/** Private: Look up endian enum */
static ObitACAEndian LookupEndian(gchar *name);

/* Get start of next MIME segment */
static void GetNextMIME(ObitALMACalAtm *in, 
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
ObitALMACalAtm* newObitALMACalAtm (gchar* name)
{
  ObitALMACalAtm* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitALMACalAtmClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitALMACalAtm));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitALMACalAtmInit((gpointer)out);

 return out;
} /* end newObitALMACalAtm */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitALMACalAtmGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitALMACalAtmClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitALMACalAtmGetClass */

/**
 * Make a deep copy of an ObitALMACalAtm. NYI
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitALMACalAtm* ObitALMACalAtmCopy  (ObitALMACalAtm *in, ObitALMACalAtm *out, ObitErr *err)
{
  /*const ObitClassInfo *ParentClass;*/
  /*gboolean oldExist;*/
  /*gchar *outName;*/

  /* error checks */
  if (err->error) return out;

  /* Stubbed */
  g_error("ObitALMACalAtmCopy: Stubbed");

  return out;
} /* end ObitALMACalAtmCopy */

 /**
 * Make a copy of a object but do not copy the actual data NYI
 * This is useful to create an ALMACalAtm similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitALMACalAtmClone  (ObitALMACalAtm *in, ObitALMACalAtm *out, ObitErr *err)
{
  /*const ObitClassInfo *ParentClass;*/

  /* error checks */
  if (err->error) return;

  /* Stubbed */
  g_error("ObitALMACalAtmClone: Stubbed");

} /* end ObitALMACalAtmClone */

/**
 * Creates an ObitALMACalAtm 
 * Parses the ASMD XML tables and stores
 * \param name         An optional name for the object.
 * \param err          Obit error stack object.
 * \return the new object.
 */
ObitALMACalAtm* ObitALMACalAtmCreate (gchar* name, ObitErr *err)
{
  ObitALMACalAtm* out;
  /*gchar *routine="ObitALMACalAtmCreate";*/

  /* Create basic structure */
  out = newObitALMACalAtm (name);

  /* Init buffer */
  out->buffer = g_malloc0(ACABUFFERSIZE*ACABUFFERFRAMES);
  out->nBytesInBuffer = 0;
  out->current = out->buffer;

  return out;
} /* end ObitALMACalAtmCreate */

 /**
 * Initialize File
 * Initializes buffer, parses XML header
 * \param in       The object to fill
 * \param DataFile Name of Mime file with data
 * \param err      Obit error stack object.
 * \param err Obit error stack object.
 */
void ObitALMACalAtmInitFile  (ObitALMACalAtm *in, gchar *DataFile, ObitErr *err)
{
  ObitIOCode retCode;
  olong i, count, maxStr, itemp, iord, nord;
  gchar *startInfo, *endInfo, *next, *start, *prior, *tstr;
  gchar *ord[40];
  gchar *tableUID=NULL, *containerUID=NULL, *junk=NULL;
  ObitACAEndian endian=ACAEndian_Little;
  gchar *routine = "ObitALMACalAtmInitFile";

  /* error checks */
  if (err->error) return;

   /* set Values - file name */
  in->DataFile = strdup(DataFile);

  /* Get size */
  in->fileSize = ObitFileSize(DataFile, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create file object */
  if (in->file==NULL) in->file = newObitFile ("ACA File");

  /* Open */
  ObitFileOpen (in->file, in->DataFile, OBIT_IO_ReadOnly, OBIT_IO_Binary, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Fill Buffer */
  retCode = ObitALMACalAtmFillBuffer (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Parse xml header - first find limits */
  maxStr    = in->nBytesInBuffer - (in->current-in->buffer);
  startInfo = g_strstr_len (in->current, maxStr, "<CalAtmosphereTable ");
  endInfo   = g_strstr_len (in->current, maxStr, "</CalAtmosphereTable>");
  maxStr    = (olong)(endInfo-startInfo);
  
  /* Endian */
  prior = "byteOrder=";
  tstr = ACAparse_quote_str (startInfo, maxStr, prior, &next);
  if (tstr) {
    endian = LookupEndian(tstr);
    g_free(tstr);
  }
  
  /* Is Byte flip needed? */
  in->byteFlip =  ((endian==ACAEndian_Big) && (G_BYTE_ORDER==G_LITTLE_ENDIAN)) ||
    ((endian==ACAEndian_Little) && (G_BYTE_ORDER==G_BIG_ENDIAN));

  /* Orders of entries - use buffer location to order */
  startInfo = g_strstr_len (in->current, maxStr, "<Attributes>");
  endInfo   = g_strstr_len (startInfo, maxStr, "</Attributes>");
  maxStr    = (olong)(endInfo-startInfo);
  iord = 0;
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<antennaName/>");             /* 0 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<receiverBand/>");            /* 1 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<basebandName/>");            /* 2 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<calDataId/>");               /* 3 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<calReductionId/>");          /* 4 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<startValidTime/>");          /* 5 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<endValidTime/>");            /* 6 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<numFreq/>");                 /* 7 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<numLoad/>");                 /* 8 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<numReceptor/>");             /* 9 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<forwardEffSpectrum/>");      /* 10 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<frequencyRange/>");          /* 11 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<groundPressure/>");          /* 12 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<groundRelHumidity/>");       /* 13 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<frequencySpectrum/>");       /* 14 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<groundTemperature/>");       /* 15 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<polarizationTypes/>");       /* 16 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<powerSkySpectrum/>");        /* 17 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<powerLoadSpectrum/>");       /* 18 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<syscalType/>");              /* 19 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<tAtmSpectrum/>");            /* 20 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<tRecSpectrum/>");            /* 21 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<tSysSpectrum/>");            /* 22 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<tauSpectrum/>");             /* 23 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<tAtm/>");                    /* 24 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<tRec/>");                    /* 25 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<tSys/>");                    /* 26 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<tau/>");                     /* 27 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<water/>");                   /* 28 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<waterError/>");              /* 29 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<alphaSpectrum/>");           /* 30 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<forwardEfficiency/>");       /* 31 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<forwardEfficiencyError/>");  /* 32 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<sbGain/>");                  /* 33 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<sbGainError/>");             /* 34 */
  ord[iord++] = g_strstr_len (startInfo, maxStr, "<sbGainSpectrum/>");          /* 35 */
  nord = iord;
  in->nelem = nord;   /* How many elements in table */

  /* Interprete orders */
  /* antennaName 0 */
  count = 0;
  in->ordantennaName = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[0])) in->ordantennaName++;
  }

  /* receiverBand 1 */
  count = 0;
  in->ordreceiverBand = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[1])) in->ordreceiverBand++;
  }

  /* basebandName 2 */
  count = 0;
  in->ordbasebandName = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[2])) in->ordbasebandName++;
  }

  /* calDataId 3 */
  count = 0;
  in->ordcalDataId = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[3])) in->ordcalDataId++;
  }

  /* calReductionId 4 */
  count = 0;
  in->ordcalReductionId = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[4])) in->ordcalReductionId++;
  }

  /* startValidTime 5 */
  count = 0;
  in->ordstartValidTime = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[5])) in->ordstartValidTime++;
  }

  /* endValidTime 6 */
  count = 0;
  in->ordendValidTime = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[6])) in->ordendValidTime++;
  }

  /* numFreq 7 */
  count = 0;
  in->ordnumFreq = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[7])) in->ordnumFreq++;
  }

  /* numLoad 8 */
  count = 0;
  in->ordnumLoad = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[8])) in->ordnumLoad++;
  }

  /* numReceptor 9 */
  count = 0;
  in->ordnumReceptor = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[9])) in->ordnumReceptor++;
  }

  /* forwardEffSpectrum 10 */
  count = 0;
  in->ordforwardEffSpectrum = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[10])) in->ordforwardEffSpectrum++;
  }

  /* frequencyRange 11 */
  count = 0;
  in->ordfrequencyRange = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[11])) in->ordfrequencyRange++;
  }

  /* groundPressure 12 */
  count = 0;
  in->ordgroundPressure = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[12])) in->ordgroundPressure++;
  }

  /* groundRelHumidity 13 */
  count = 0;
  in->ordgroundRelHumidity = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[13])) in->ordgroundRelHumidity++;
  }

  /* frequencySpectrum 14 */
  count = 0;
  in->ordfrequencySpectrum = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[14])) in->ordfrequencySpectrum++;
  }

  /* groundTemperature 15 */
  count = 0;
  in->ordgroundTemperature = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[15])) in->ordgroundTemperature++;
  }

  /* polarizationTypes 16 */
  count = 0;
  in->ordpolarizationTypes = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[16])) in->ordpolarizationTypes++;
  }

  /* powerSkySpectrum 17 */
  count = 0;
  in->ordpowerSkySpectrum = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[17])) in->ordpowerSkySpectrum++;
  }

  /* powerLoadSpectrum 18 */
  count = 0;
  in->ordpowerLoadSpectrum = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[18])) in->ordpowerLoadSpectrum++;
  }

  /* syscalType 19 */
  count = 0;
  in->ordsyscalType = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[19])) in->ordsyscalType++;
  }

  /* tAtmSpectrum 20 */
  count = 0;
  in->ordtAtmSpectrum = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[20])) in->ordtAtmSpectrum++;
  }

  /* tRecSpectrum 21 */
  count = 0;
  in->ordtRecSpectrum = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[21])) in->ordtRecSpectrum++;
  }

  /* tSysSpectrum 22 */
  count = 0;
  in->ordtSysSpectrum = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[22])) in->ordtSysSpectrum++;
  }

  /* tauSpectrum 23 */
  count = 0;
  in->ordtauSpectrum = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[23])) in->ordtauSpectrum++;
  }

  /* tAtm 24 */
  count = 0;
  in->ordtAtm = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[24])) in->ordtAtm++;
  }

  /* tRec 25 */
  count = 0;
  in->ordtRec = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[25])) in->ordtRec++;
  }

  /* tSys 26 */
  count = 0;
  in->ordtSys = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[26])) in->ordtSys++;
  }

  /* tau 27 */
  count = 0;
  in->ordtau = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[27])) in->ordtau++;
  }

  /* water 28 */
  count = 0;
  in->ordwater = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[28])) in->ordwater++;
  }

  /* waterError 29 */
  count = 0;
  in->ordwaterError = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[29])) in->ordwaterError++;
  }

  /* alphaSpectrum 30 */
  count = 0;
  in->ordalphaSpectrum = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[30])) in->ordalphaSpectrum++;
  }

  /* forwardEfficiency 31 */
  count = 0;
  in->ordforwardEfficiency = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[31])) in->ordforwardEfficiency++;
  }

  /* forwardEfficiencyError 32 */
  count = 0;
  in->ordforwardEfficiencyError = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[32])) in->ordforwardEfficiencyError++;
  }

  /* sbGain 33 */
  count = 0;
  in->ordsbGain = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[33])) in->ordsbGain++;
  }

  /* sbGainError 34 */
  count = 0;
  in->ordsbGainError = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[34])) in->ordsbGainError++;
  }

  /* sbGainSpectrum 35 */
  count = 0;
  in->ordsbGainSpectrum = 0;
  for (i=0; i<nord; i++) {
    if ((i!=0) && (ord[i]<=ord[35])) in->ordsbGainSpectrum++;
  }

  in->current = endInfo;  /* where in buffer */

  /* Go to data table */
  GetNextMIME (in, endInfo, &start, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Parse stuff at beginning including number of rows */
  /* Real amateur production here - this isn't documented but appears to work */
  tableUID     = ACAparse_str(in->buffer, start, &next, in->byteFlip);
  start        = next;
  junk         = ACAparse_str(in->buffer, start, &next, in->byteFlip);
  if (junk)    g_free(junk); junk = NULL;
  start        = next;
  junk         = ACAparse_str(in->buffer, start, &next, in->byteFlip);
  if (junk)    g_free(junk); junk = NULL;
  start        = next;
  junk         = ACAparse_str(in->buffer, start, &next, in->byteFlip);
  if (junk)    g_free(junk); junk = NULL;
  start        = next;
  junk         = ACAparse_str(in->buffer, start, &next, in->byteFlip);
  if (junk)    g_free(junk); junk = NULL;
  start        = next;
  containerUID = ACAparse_str(in->buffer, start, &next, in->byteFlip);
  start        = next;
  junk         = ACAparse_str(in->buffer, start, &next, in->byteFlip);
  if (junk)    g_free(junk); junk = NULL;
  start        = next;
  junk         = ACAparse_str(in->buffer, start, &next, in->byteFlip);
  if (junk)    g_free(junk); junk = NULL;
  start        = next;
  itemp        = ACAparse_int(in->buffer, start, &next, in->byteFlip);
  start        = next;
  itemp        = ACAparse_int(in->buffer, start, &next, in->byteFlip);
  start        = next+2;
  in->nrow     = ACAparse_int(in->buffer, start, &next, in->byteFlip);
  start        = next;
  in->curRow   = 0;
  if (tableUID) g_free(tableUID);
  if (containerUID) g_free(containerUID);

  /* What a BIZZARRE format - loop through strings until finding a prior
     "length" <=0 or > 1000 and presume that this is the number of rows */
  /*while (1) {
    itemp = ACAparse_int(in->buffer, start, &next, in->byteFlip);
    if ((itemp<=0) || (itemp>1000)) break;
    tstr    = ACAparse_str(in->buffer, start, &next, in->byteFlip);
    start   = next;
    if (tstr) g_free(tstr);
    }

  in->nrow     = ACAparse_int(in->buffer, start, &next, in->byteFlip);*/

  in->current = next;  /* where in buffer */

} /* end ObitALMACalAtmInitFile */

 /**
 * Fill Buffer
 * If the buffer is filled to capacity, the bottom frame is copied to
 * to the top of the buffer and the remainder of the buffer filled.
 * Updates in->nBytesInBuffer, in->current.
 * \param in  The object to fill
 * \param err Obit error stack object.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_EOF = EOF.
 */
ObitIOCode ObitALMACalAtmFillBuffer (ObitALMACalAtm *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  olong i, ncopy, size; 
  gchar *top, *copy;
  gchar *routine = "ObitALMACalAtmFillBuffer";

  /* error checks */
  if (err->error) return retCode;

  if ((in->fileSize-in->file->filePos) <= 0) return OBIT_IO_EOF;

   /* Is it already full? */
  if (in->nBytesInBuffer>=ACABUFFERSIZE*ACABUFFERFRAMES) { /* Yes - shift */
    top  = in->buffer;
    copy = in->buffer + (ACABUFFERFRAMES-1) * ACABUFFERSIZE;
    memmove (top, copy, (size_t)ACABUFFERSIZE);
    top = &in->buffer[ACABUFFERSIZE];
    ncopy = (ACABUFFERFRAMES-1);
    in->nBytesInBuffer = ACABUFFERSIZE;
    in->current -= (ACABUFFERFRAMES-1) * ACABUFFERSIZE; /* Current position */
  } else {  /* Nope - fill 'er up */
    top  = in->buffer;
    ncopy = ACABUFFERFRAMES;
    in->nBytesInBuffer = 0;
    in->current = in->buffer; /* Current position */
  }

  /* Do reads */
  for (i=0; i<ncopy; i++) {
    /* No more than what's left */
    size = MIN ((olong)ACABUFFERSIZE, (in->fileSize-in->file->filePos));
    if (size<=0) break;
    retCode = ObitFileRead (in->file, -1, size, top, err);
    /* DEBUG
    fprintf (stderr, "Read size %d pos %lld filesize %lld\n", size, in->file->filePos, in->fileSize); */
    if (err->error) {
      Obit_traceback_val (err, routine, in->name, retCode);
    }
    in->nBytesInBuffer += size;
    if (in->file->filePos>=in->fileSize) break;
    top += ACABUFFERSIZE;
  }
  return retCode;
} /* end ObitALMACalAtmFillBuffer */

/**
 * Return number of rows in table
 * Should be called after ObitALMACalAtmInitFile
 * \param in      Pointer to table
 * \return number of rows given by header
 */
olong ObitALMACalAtmGetNrow (ObitALMACalAtm *in)
{
  return in->nrow;
} /* end ObitALMACalAtmGetNrow*/

/**
 * Fill in values for next table row
 * \param in      Pointer to object to be read.
 * \param row     Row structure to fill in
 * \param err     ObitErr for reporting errors.
 * \return row number from table, -1 = done
 */
olong ObitALMACalAtmGetRow (ObitALMACalAtm *in, 
			    ASDMcalAtmosphereRow *row, ObitErr *err)
{
  olong out = in->curRow+1;
  olong maxStr, elem, nelem, i, n, ijunk;
  odouble *dbljunk=NULL;
  ofloat *fltjunk=NULL, fjunk;
  gchar *start, *next, *done, *strjunk=NULL;
  gboolean gotIt; 
  gchar *routine = "ObitALMACalAtmGetRow";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (row != NULL);

  /* Done? */
  if (out>in->nrow) return -1;
  nelem = in->nelem;

  /* Look for --MIME_boundary-- marking end of data */
  maxStr = 30;
  done = g_strstr_len (in->current, maxStr, "--MIME_boundary--");
  if (done!=NULL) return -1;

  start = in->current;
  /* Loop over elements in row */
  for (elem=0; elem<nelem; elem++) {

    /* antennaName */
    if (elem==in->ordantennaName) {
      row->antennaName =  ACAparse_str(in->buffer, start, &next, in->byteFlip);
      /* receiverBand */
    } else if (elem==in->ordreceiverBand) {
      row->receiverBand = ACAparse_intstr(in->buffer, "ALMA_RB_", start, &next, in->byteFlip);
      /* basebandName */
    } else if (elem==in->ordbasebandName) {
      row->basebandId = ACAparse_intstr(in->buffer, "BB_", start, &next, in->byteFlip);
      /* calDataId */
    } else if (elem==in->ordcalDataId) {
      row->calDataId = ACAparse_intstr(in->buffer, "CalData_", start, &next, in->byteFlip);
      /* calReductionId */
    } else if (elem==in->ordcalReductionId) {
      row->calReductionId = ACAparse_intstr(in->buffer, "CalReduction_", start, &next, in->byteFlip);
      /* startValidTime */
    } else if (elem==in->ordstartValidTime) {
      row->startValidTime = ACAparse_time(in->buffer, start, &next, in->byteFlip);
      /* endValidTime */
    } else if (elem==in->ordendValidTime) {
      row->endValidTime = ACAparse_time(in->buffer, start, &next, in->byteFlip);
      /* numFreq */
    } else if (elem==in->ordnumFreq) {
      row->numFreq = ACAparse_int(in->buffer, start, &next, in->byteFlip);
      /* numLoad */
    } else if (elem==in->ordnumLoad) {
      row->numLoad = ACAparse_int(in->buffer, start, &next, in->byteFlip);
      /* numReceptor */
    } else if (elem==in->ordnumReceptor) {
      row->numReceptor = ACAparse_int(in->buffer, start, &next, in->byteFlip);
      /* forwardEffSpectrum */
    } else if (elem==in->ordforwardEffSpectrum) {
      fltjunk = ACAparse_flt_arr(in->buffer, start, &next, 1, row->numReceptor*row->numFreq, in->byteFlip);
      if (fltjunk) g_free(fltjunk); fltjunk = NULL;
      /* frequencyRange */
    } else if (elem==in->ordfrequencyRange) {
      row->frequencyRange = g_malloc0(2*sizeof(odouble));
      ijunk = ACAparse_int(in->buffer, start, &next, in->byteFlip); 
      start = next; row->frequencyRange[0] = ACAparse_dbl(in->buffer, start, &next, in->byteFlip); 
      start = next; row->frequencyRange[1] = ACAparse_dbl(in->buffer, start, &next, in->byteFlip); 
      /* groundPressure */
    } else if (elem==in->ordgroundPressure) {
      row->groundPressure = ACAparse_dbl(in->buffer, start, &next, in->byteFlip);
      /* groundRelHumidity */
    } else if (elem==in->ordgroundRelHumidity) {
      row->groundRelHumidity = ACAparse_dbl(in->buffer, start, &next, in->byteFlip);
      /* frequencySpectrum */
    } else if (elem==in->ordfrequencySpectrum) {
      dbljunk = ACAparse_dbl_arr(in->buffer, start, &next, 0, row->numFreq, in->byteFlip);
      if (dbljunk) g_free(dbljunk); dbljunk = NULL;
      /* groundTemperature */
    } else if (elem==in->ordgroundTemperature) {
      row->groundTemperature = ACAparse_dbl(in->buffer, start, &next, in->byteFlip);
      /* polarizationTypes */
    } else if (elem==in->ordpolarizationTypes) {
      strjunk = ACAparse_poln(in->buffer, start, &next, in->byteFlip);
      if (strjunk) g_free(strjunk); strjunk = NULL;
      /* Not really what's needed here */
      /* powerSkySpectrum */
    } else if (elem==in->ordpowerSkySpectrum) {
      fltjunk = ACAparse_flt_arr(in->buffer, start, &next, 1, row->numReceptor*row->numFreq, in->byteFlip);
      if (fltjunk) g_free(fltjunk); fltjunk = NULL;
      /* powerLoadSpectrum */
    } else if (elem==in->ordpowerLoadSpectrum) {
      n = row->numReceptor*row->numFreq*MAX (1, row->numLoad);
      fltjunk = ACAparse_flt_arr(in->buffer, start, &next, 2, n, in->byteFlip);
      if (fltjunk) g_free(fltjunk); fltjunk = NULL;
      /* syscalType */
    } else if (elem==in->ordsyscalType) {
      row->syscalType = ACAparse_str(in->buffer, start, &next, in->byteFlip);
      /* tAtmSpectrum  */
    } else if (elem==in->ordtAtmSpectrum) {
      dbljunk = ACAparse_dbl_arr(in->buffer, start, &next, 1, row->numReceptor*row->numFreq, in->byteFlip);
      if (dbljunk) g_free(dbljunk); dbljunk = NULL;
      /* tRecSpectrum */
    } else if (elem==in->ordtRecSpectrum) {
      dbljunk = ACAparse_dbl_arr(in->buffer, start, &next, 1, row->numReceptor*row->numFreq, in->byteFlip);
      if (dbljunk) g_free(dbljunk); dbljunk = NULL;
      /* tSysSpectrum */
    } else if (elem==in->ordtSysSpectrum) {
      dbljunk = ACAparse_dbl_arr(in->buffer, start, &next, 1, row->numReceptor*row->numFreq, in->byteFlip);
      if (dbljunk) g_free(dbljunk); dbljunk = NULL;
      /* tauSpectrum */
    } else if (elem==in->ordtauSpectrum) {
      fltjunk = ACAparse_flt_arr(in->buffer, start, &next, 1, row->numReceptor*row->numFreq, in->byteFlip);
      if (fltjunk) g_free(fltjunk); fltjunk = NULL;
      /* tAtm */
    } else if (elem==in->ordtAtm) {
      row->tAtm = ACAparse_dbl_arr(in->buffer, start, &next, 0, row->numReceptor, in->byteFlip);
      /* tRec */
    } else if (elem==in->ordtRec) {
      row->tRec = ACAparse_dbl_arr(in->buffer, start, &next, 0, row->numReceptor, in->byteFlip);
      /* tSys */
    } else if (elem==in->ordtSys) {
      row->tSys = ACAparse_dbl_arr(in->buffer, start, &next, 0, row->numReceptor, in->byteFlip);
      /* tau */
    } else if (elem==in->ordtau) {
      fltjunk = ACAparse_flt_arr(in->buffer, start, &next, 0, row->numReceptor, in->byteFlip);
      row->tau = g_malloc0(row->numReceptor*sizeof(odouble));
      for (i=0; i<row->numReceptor; i++) row->tau[i] = fltjunk[i];
      if (fltjunk) g_free(fltjunk); fltjunk = NULL;
      /* water */
    } else if (elem==in->ordwater) {
      row->water = ACAparse_dbl_arr(in->buffer, start, &next, 0, row->numReceptor, in->byteFlip);
      /* waterError */
    } else if (elem==in->ordwaterError) {
      row->waterError = ACAparse_dbl_arr(in->buffer, start, &next, 0, row->numReceptor, in->byteFlip);
      /* alphaSpectrum optional */
    } else if (elem==in->ordalphaSpectrum) {
      gotIt = ACAparse_bool(in->buffer, start, &next);
      start = next;
      if (gotIt) {
	fltjunk = ACAparse_flt_arr(in->buffer, start, &next, 1, row->numReceptor*row->numFreq, in->byteFlip);
	if (fltjunk) g_free(fltjunk); fltjunk = NULL;
      }
      /* forwardEfficiency  optional */
    } else if (elem==in->ordforwardEfficiency) {
      gotIt = ACAparse_bool(in->buffer, start, &next);
      start = next;
      if (gotIt) {
	fltjunk = ACAparse_flt_arr(in->buffer, start, &next, 0, row->numReceptor, in->byteFlip);
	row->forwardEfficiency = g_malloc0(row->numReceptor*sizeof(odouble));
	for (i=0; i<row->numReceptor; i++) row->forwardEfficiency[i] = fltjunk[i];
	if (fltjunk) g_free(fltjunk); fltjunk = NULL;
      }
      /* forwardEfficiencyError optional */
    } else if (elem==in->ordforwardEfficiencyError) {
      gotIt = ACAparse_bool(in->buffer, start, &next);
      start = next;
      if (gotIt) {
	dbljunk = ACAparse_dbl_arr(in->buffer, start, &next, 0, row->numReceptor, in->byteFlip);
	if (dbljunk) g_free(dbljunk); dbljunk = NULL;
      }
      /* sbGain optional */
    } else if (elem==in->ordsbGain) {
      gotIt = ACAparse_bool(in->buffer, start, &next);
      start = next;
      if (gotIt) {
	fltjunk = ACAparse_flt_arr(in->buffer, start, &next, 0, row->numReceptor, in->byteFlip);
	row->sbGain = g_malloc0(row->numReceptor*sizeof(odouble));
	for (i=0; i<row->numReceptor; i++) row->sbGain[i] = fltjunk[i];
	if (fltjunk) g_free(fltjunk); fltjunk = NULL;
      }
      /* sbGainError optional */
    } else if (elem==in->ordsbGainError) {
      gotIt = ACAparse_bool(in->buffer, start, &next);
      start = next;
      if (gotIt) {
	fltjunk = ACAparse_flt_arr(in->buffer, start, &next, 0, row->numReceptor*row->numFreq, in->byteFlip);
	if (fltjunk) g_free(fltjunk); fltjunk = NULL;
      }
      /* sbGainSpectrum optional */
    } else if (elem==in->ordsbGainSpectrum) {
      gotIt = ACAparse_bool(in->buffer, start, &next);
      start = next;
      if (gotIt) {
	fltjunk = ACAparse_flt_arr(in->buffer, start, &next, 1, row->numReceptor*row->numFreq, in->byteFlip);
	if (fltjunk) g_free(fltjunk); fltjunk = NULL; 
      }
    } else {
      /* Bummer */
    }
    start = next;    /* for next */
  } /* end loop over elements */

  in->curRow++;          /* Keep track of row */
  in->current = start;   /* where in buffer */

  /* If out of first segment of buffer read more */
  if (((in->current-in->buffer)>ACABUFFERSIZE) && (in->curRow<in->nrow)) {
    ObitALMACalAtmFillBuffer (in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
  }
  
  return out;

} /* end ObitALMACalAtmGetRow  */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitALMACalAtmClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitALMACalAtmClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitALMACalAtmClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitALMACalAtmClassInfoDefFn (gpointer inClass)
{
  ObitALMACalAtmClassInfo *theClass = (ObitALMACalAtmClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitALMACalAtmClassInit;
  theClass->newObit       = (newObitFP)newObitALMACalAtm;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitALMACalAtmClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitALMACalAtmGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitALMACalAtmCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitALMACalAtmClear;
  theClass->ObitInit      = (ObitInitFP)ObitALMACalAtmInit;
  theClass->ObitALMACalAtmCreate = (ObitALMACalAtmCreateFP)ObitALMACalAtmCreate;

} /* end ObitALMACalAtmClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitALMACalAtmInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitALMACalAtm *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->DataFile                  = NULL;
  in->buffer                    = NULL;
  in->ordantennaName            =  0;
  in->ordreceiverBand           =  1;
  in->ordbasebandName           =  2;
  in->ordcalDataId              =  3;
  in->ordcalReductionId         =  4;
  in->ordstartValidTime         =  5;
  in->ordendValidTime           =  6;
  in->ordnumFreq                =  7;
  in->ordnumLoad                =  8;
  in->ordnumReceptor            =  9;
  in->ordforwardEffSpectrum     =  10;
  in->ordfrequencyRange         =  11;
  in->ordgroundPressure         =  12;
  in->ordgroundRelHumidity      =  13;
  in->ordfrequencySpectrum      =  14;
  in->ordgroundTemperature      =  15;
  in->ordpolarizationTypes      =  16;
  in->ordpowerSkySpectrum       =  17;
  in->ordpowerLoadSpectrum      =  18;
  in->ordsyscalType             =  19;
  in->ordtAtmSpectrum           =  20;
  in->ordtRecSpectrum           =  21;
  in->ordtSysSpectrum           =  22;
  in->ordtauSpectrum            =  23;
  in->ordtAtm                   =  24;
  in->ordtRec                   =  25;
  in->ordtSys                   =  26;
  in->ordtau                    =  27;
  in->ordwater                  =  28;
  in->ordwaterError             =  29;
  in->ordalphaSpectrum          =  30;
  in->ordforwardEfficiency      =  31;
  in->ordforwardEfficiencyError =  32;
  in->ordsbGain                 =  33;
  in->ordsbGainError            =  34;
  in->ordsbGainSpectrum         =  35;
  in->nelem                     =  36;
} /* end ObitALMACalAtmInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitALMACalAtm* cast to an Obit*.
 */
void ObitALMACalAtmClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitALMACalAtm *in = inn;
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
  
} /* end ObitALMACalAtmClear */

/* ALMACalAtm (ACA) parsing routines */


/**  Parse single quoted (ALMA quotes) from XLM string 
 * All text from end of prior+1 until next '''
 * \param  string  String to parse
 * \param  maxChar Maximum size of string
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \return value, NULL if problem, should be g_freeed when done
 */
static gchar* ACAparse_quote_str(gchar *string, olong maxChar, 
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
  if (*b!='\'') return out;  /* Make sure quote */
  b++;                      /* Skip quote */

  /* count */
  charLeft = maxChar - (b-string);
  n = 0;
  for (i=0; i<charLeft; i++) {
    if (b[i]=='\'') break;
    n++;
  }
  out = g_malloc(n+1);
  for (i=0; i<n; i++) out[i] = b[i]; out[i] = 0;
  *next = b + n;

  return out;
} /* end ACAparse_quote_str */

/**  Parse 8 bit boolean
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \return value
 */
gboolean ACAparse_bool(gchar *buffer, 
		       gchar *prior, gchar **next)
{
  gboolean out = FALSE;
  gchar *b = prior;

  *next = prior+1;  /* After byte */
  out = b[0] == 1;
    
  return out;
} /* end ACAparse_bool */

/**  Parse 32 bit integer byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
olong ACAparse_int(gchar *buffer, 
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
} /* end ACAparse_int */

/**  Parse 32 bit integer from string
 * \param  buffer   buffer to parse
 * \param  prefix   Prefix for integer
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
olong ACAparse_intstr(gchar *buffer, gchar *prefix,
		      gchar *prior, gchar **next, 
		      gboolean byteFlip)
{
  olong out = 0;
  olong maxChar;
  gchar *temp, *b, *e;

  /* read string */
  temp = ACAparse_str(buffer, prior, next, byteFlip);

  maxChar = strlen(temp);
  b = g_strstr_len (temp, maxChar, prefix);
  if (b==NULL) return out;  /* Found? */
  b += strlen(prefix);
  out = (olong)strtol(b, &e, 10);

  if (temp) g_free(temp);
  return out;
} /* end ACAparse_intstr */

/**  Parse 64 bit integer byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
long long ACAparse_llong(gchar *buffer, 
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
} /* end ACAparse_llong */

/**  Parse 32 bit float byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
ofloat ACAparse_flt(gchar *buffer, 
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
} /* end ACAparse_flt */

/**  Parse 64 bit float byte flipping as necessary
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value
 */
odouble ACAparse_dbl(gchar *buffer, 
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
} /* end ACAparse_dbl */

/**  Parse string from buffer, length of string given by preceeding integer.
 * \param  buffer   buffer to parse
 * \param  maxChar  Maximum size of string
 * \param  prior    buffer prior to value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return string pointer, null terminated, should be g_freeed when done
 */
gchar* ACAparse_str(gchar *buffer, 
		    gchar *prior, gchar **next, 
		    gboolean byteFlip)
{
  gchar *out = NULL;
  gchar *b;
  olong i, len;

  /* How long */
  len = ACAparse_int(buffer, prior, next, byteFlip);

  /* output string */
  out = g_malloc(len+1);
  b = *next;
  for (i=0; i<len; i++) out[i] = b[i]; out[i] = 0;
  *next += len;

  return out;
} /* end ACAparse_str */

/**  Parse poln strings from buffer, 
 * Assumes first integer number of strings
 * \param  buffer   buffer to parse
 * \param  maxChar  Maximum size of string
 * \param  prior    buffer prior to value
 * \param  next     pointer in buffer after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return string pointer, null terminated, should be g_freeed when done
 *  all poln in one string
 */
gchar* ACAparse_poln(gchar *buffer, 
		    gchar *prior, gchar **next, 
		    gboolean byteFlip)
{
  gchar *out = NULL, *t=NULL;
  gchar *b;
  olong i, len;

  /* How many */
  len = ACAparse_int(buffer, prior, next, byteFlip);

  /* output string */
  out = g_malloc(len+1);
  b = *next;
  for (i=0; i<len; i++) {
    t = ACAparse_str(buffer, b, next, byteFlip);
    out[i] = t[0]; 
    if (t) g_free(t); t = NULL;
    b = *next;
  }
  out[i] = 0;

  return out;
} /* end ACAparse_poln */

/**  Parse time interval
 * Read time interval (MJD) in nanoseconds, return as JD
 * \param  buffer  Bufferto parse
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value as pair of odouble
 */
odouble* ACAparse_timeint(gchar *buffer, 
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
  temp = ACAparse_llong(buffer, prior, next, byteFlip);
  out[0] = (odouble)((temp*1.0e-9)/86400.0) + mjdJD0;
  /* Interval */
  b = *next;
  temp = ACAparse_llong(buffer, b, next, byteFlip);
  out[1] = (odouble)((temp*1.0e-9)/86400.0) + mjdJD0;

  return out;
} /* end ACAparse_timeint */

/**  Parse time 
 * Read time (MJD) in nanoseconds, return as JD
 * \param  buffer  Bufferto parse
 * \param  prior string prior to value
 * \param  next  pointer in string after parsed value
 * \param  byteFlip True if bytes need flipping
 * \return value as pair of odouble
 */
odouble ACAparse_time(gchar *buffer, 
			  gchar *prior, gchar **next,
			  gboolean byteFlip)
{
  odouble out;
  long long temp;
  odouble mjdJD0=2400000.5; /* JD of beginning of MJD time */

  /* Time */
  temp = ACAparse_llong(buffer, prior, next, byteFlip);
  out = (odouble)((temp*1.0e-9)/86400.0) + mjdJD0;

  return out;
} /* end ACAparse_time */

/**  Parse 1-D array of 32 bit float byte flipping as necessary
 *  expects 4 byte integer of questionable value before data
 *  skip ncrap integers after questionable integer.
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  ncrap    how many extra useless integers before data?
 * \param  nwords   How many words of data actually there
 * \param  byteFlip True if bytes need flipping
 * \return value array, should be g_freeed when done
 */
ofloat* ACAparse_flt_arr(gchar *buffer, 
			 gchar *prior, gchar **next, olong ncrap, olong nwords,
			 gboolean byteFlip)
{
  ofloat *out = NULL;
  olong ndim,  i;
  /*olong num, naxis[10]; */
  gchar *b;
  
  /* How many dimensions? */
  ndim = ACAparse_int(buffer, prior, next, byteFlip);
  b = *next;
  /* Useless words?  */
  if (ncrap>0) {
    i = ACAparse_int(buffer, b, next, byteFlip);
    b = *next;
  }
  if (ncrap>1) {
    i = ACAparse_int(buffer, b, next, byteFlip);
    b = *next;
  }
  if (ncrap>2) {
    i = ACAparse_int(buffer, b, next, byteFlip);
    b = *next;
  }

  /* Size of each dimension - probably bad values - just skip
  documentation VERY WRONG
  num = 1;
  for (i=0; i<ndim; i++) {
    naxis[i] = ACAparse_int(buffer, b, next, byteFlip);
    b = *next;
    num *= MAX(1, naxis[i]);
    } */

  /* output */
  out = g_malloc0(nwords*sizeof(ofloat));

  /* Loop parsing */
  for (i=0; i<nwords; i++) {
    out[i] = ACAparse_flt(buffer, b, next, byteFlip);
    b = *next;
  }

  return out;
} /* end ACAparse_flt_arr */

/**  Parse 1-D array of 32 bit float byte flipping as necessary
 *  expects 4 byte integer of questionable value before data
 *  skip ncrap integers after questionable integer.
 * \param  buffer   buffer to parse
 * \param  prior    buffer start of value
 * \param  next     pointer in buffer after parsed value
 * \param  ncrap    how many extra useless integers before data?
 * \param  nwords   How many words of data actually there
 * \param  byteFlip True if bytes need flipping
 * \return value array, should be g_greeed when done
 */
odouble* ACAparse_dbl_arr(gchar *buffer, 
			  gchar *prior, gchar **next, olong ncrap, olong nwords,
			  gboolean byteFlip)
{
  odouble *out = NULL;
  olong ndim, num, i;
  /* olong naxis[10]; */
  gchar *b;

  /* How many dimensions? */
  ndim = ACAparse_int(buffer, prior, next, byteFlip);
  b = *next;
  /* Useless words?  */
  if (ncrap>0) {
    i = ACAparse_int(buffer, b, next, byteFlip);
    b = *next;
  }
  if (ncrap>1) {
    i = ACAparse_int(buffer, b, next, byteFlip);
    b = *next;
  }
  if (ncrap>2) {
    i = ACAparse_int(buffer, b, next, byteFlip);
    b = *next;
  }
  /* Size of each dimension
     NO - implementation NOT what documented 
  num = 1;
  for (i=0; i<ndim; i++) {
    naxis[i] = ACAparse_int(buffer, b, next, byteFlip);
    b = *next;
    num *= MAX(1, naxis[i]);
    }*/
  num = nwords;

  /* output */
  out = g_malloc0(num*sizeof(odouble));

  /* Loop parsing */
  for (i=0; i<num; i++) {
    out[i] = ACAparse_dbl(buffer, b, next, byteFlip);
    b = *next;
  }

  return out;
} /* end ACAparse_dbl_arr */

 /**
 * Find beginning of next Mime segment
 * May update buffer contents.
 * \param in    The object to update
 * \param last  Last byte in buffer of previous segment
 * \param start First byte in buffer of data segment
 * \param err   Obit error stack object.
 */
static void GetNextMIME(ObitALMACalAtm *in, 
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
static ObitACAEndian LookupEndian(gchar *name)
{
  ObitACAEndian out = 0;
 
  if (!strncmp (name, "Little_Endian", 13)) return ACAEndian_Little;
  if (!strncmp (name, "Big_Endian",    10)) return ACAEndian_Big;
  return out;
} /* end LookupEndian */

