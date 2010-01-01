/* $Id:  $  */
/* Convert Obit UV to FITS IDI format                                 */
/* to do
   1) Check Convertion of antenna positions to earth center in ARRAY_GEOMETRY 
      table if VLA
   3) FREQID not defined in context of AIPS AN table
   
 */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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
#include "ObitFITS.h"
#include "ObitSystem.h"
#include "ObitAIPSDir.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitUV.h"
#include "ObitFileFITS.h"
#include "ObitTableSU.h"
#include "ObitTableAN.h"
#include "ObitTableFQ.h"
#include "ObitTableFG.h"
#include "ObitTableCL.h"
#include "ObitTableBP.h"
#include "ObitTableTY.h"
#include "ObitTableWX.h"
#include "ObitTableIDI_ANTENNA.h"
#include "ObitTableIDI_ARRAY_GEOMETRY.h"
#include "ObitTableIDI_FREQUENCY.h"
#include "ObitTableIDI_SOURCE.h"
#include "ObitTableIDI_FLAG.h"
#include "ObitTableIDI_CALIBRATION.h"
#include "ObitTableIDI_BANDPASS.h"
#include "ObitTableIDI_SYSTEM_TEMPERATURE.h"
#include "ObitTableIDI_WEATHER.h"
#include "ObitTableIDI_UV_DATA.h"
#include "ObitHistory.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* IDIOutin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Create output uvdata */
ObitUV* setInputData (ObitInfoList *myInput, ObitErr *err);
/* Init output FITS IDI file */
void InitIDI (gchar *FITSfile, olong disk, ObitErr *err);
/* Put Frequency info */
void PutFrequencyInfo (ObitUV *inData, ObitData *outData, ObitErr *err);
/* Put Antenna info */
void PutAntennaInfo (ObitInfoList *myInput, ObitUV *inData, ObitData *outData, 
		     ObitErr *err);
/* Put Source info */
void PutSourceInfo (ObitUV *inData, ObitData *outData, ObitErr *err);
/* Write output visibility data */
void PutData (ObitUV *inData, ObitData *outData, ObitInfoList *myInput, 
	      ObitErr *err);
/* Copy any FG/FLAG tables */
void PutFlagInfo (ObitInfoList *myInput, ObitUV *inData, ObitData *outData, 
		  ObitErr *err);
/* Copy any CL/CALIBRATION tables */
void PutCalibrationInfo (ObitInfoList *myInput, ObitUV *inData, ObitData*outData, 
			 ObitErr *err);
/* Copy any BP/BANDPASS tables */
void PutBandpassInfo (ObitInfoList *myInput, ObitUV *inData, ObitData *outData, 
					     ObitErr *err);
/* Copy any TY/SYSTEM_TEMPERATURE tables */
void PutTSysInfo (ObitInfoList *myInput, ObitUV *inData, ObitData *outData, 
		  ObitErr *err);
/* Copy any WX/WEATHER tables */
void PutWeatherInfo (ObitInfoList *myInput, ObitUV *inData, ObitData *outData, 
		     ObitErr *err);
/* Write history */
void IDIOutHistory (ObitUV* inData, ObitInfoList* myInput, ObitData* outData, 
		    ObitErr* err);
/* Special version of routine to create IDI_UV_DATA table */
ObitTableIDI_UV_DATA* myObitTableIDI_UV_DATAValue (gchar* name, ObitData *file, olong *ver,
						   ObitUV *inData, ObitIOAccess access,
						   oint no_band, ObitErr *err);
/** copy table keywords to descriptor info list */
void myObitTableIDI_UV_DATADumpKey (ObitTableIDI_UV_DATA *in, ObitErr *err);
/** Write UV_DATA row */
ObitIOCode 
myObitTableIDI_UV_DATAWriteRow  (ObitTableIDI_UV_DATA *in, olong iIDI_UV_DATARow, 
				 ObitTableIDI_UV_DATARow *row,
				 ObitErr *err);
/** Close UV_DATA  */
ObitIOCode myObitTableIDI_UV_DATAClose (ObitTableIDI_UV_DATA *in, ObitErr *err);

/* Program globals */
gchar *pgmName = "IDIOut";     /* Program name */
gchar *infile  = "IDIOut.inp"; /* File with program inputs */
gchar *outfile = "IDIOut.out"; /* File to contain program outputs */
olong  pgmNumber;              /* Program number (like POPS no.) */
olong  AIPSuser;               /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;                /* Number of AIPS directories */
gchar **AIPSdirs=NULL;         /* List of AIPS data directories */
olong  nFITS=0;                /* Number of FITS directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
gchar **FITSdirs=NULL;         /* List of FITS data directories */
olong maxAnt;                  /* Maximum antenna number */


int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*    Convert Obit UV to FITS IDI format                                  */
/*----------------------------------------------------------------------- */
{
  olong  i, ierr=0;
  ObitSystem *mySystem= NULL;
  ObitData *outData= NULL;
  ObitUV   *inData=NULL;
  ObitErr *err= NULL;
  gchar outIDI[128];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong disk;

  err = newObitErr();  /* Obit error/message stack */

  /* Startup - parse command line */
  ierr = 0;
  myInput = IDIOutin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* output File name */
  for (i=0; i<128; i++) outIDI[i] = 0;
  ObitInfoListGet(myInput, "outFile", &type, dim, outIDI, err);
  outIDI[dim[0]] = 0;     /* null terminate */
  ObitTrimTrail(outIDI);  /* Trim trailing blanks */
  ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;  ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Output file  */
  outData = (ObitData*)newObitData("Output Data");
  ObitDataSetFITS(outData, disk, outIDI, err);

  /* Initialize Output - write file primary HDU */
  InitIDI (outIDI, disk, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Open output data */
  if ((ObitDataOpen (outData, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", outData->name);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Create ObitUV for input data */
  inData = setInputData (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Write Frequency info */
  PutFrequencyInfo (inData, outData, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
   
  /* Write Antenna info */
  PutAntennaInfo (myInput, inData, outData, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
   
  /* Write Source info */
  PutSourceInfo (inData, outData, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
   
  /* Copy other tables */
  PutFlagInfo (myInput, inData, outData, err);          /* FLAG tables */
  PutCalibrationInfo (myInput, inData, outData, err);   /* CALIBRATION tables */
  PutBandpassInfo (myInput, inData, outData, err);      /* BANDPASS tables */
  PutTSysInfo (myInput, inData, outData, err);          /* SYSTEM_TEMPERATURE tables */
  PutWeatherInfo (myInput, inData, outData, err);       /* WEATHER tables */
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* History as Obit History table */
  IDIOutHistory (inData, myInput, outData, err);
  
  /* convert data  */
  PutData (inData, outData, myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;
  
  /* Shutdown Obit */
 exit:
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  mySystem = ObitSystemShutdown (mySystem);
  
  /* cleanup */
  myInput  = ObitInfoListUnref(myInput);   /* delete input list */
  myOutput = ObitInfoListUnref(myOutput);  /* delete output list */
  inData   = ObitUnref(inData);
  outData  = ObitUnref(outData);
 
  return ierr;
} /* end of main */

ObitInfoList* IDIOutin (int argc, char **argv, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Parse control info from command line                                  */
/*   Input:                                                               */
/*      argc   Number of arguments from command line                      */
/*      argv   Array of strings from command line                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   return  parser list                                                  */
/*----------------------------------------------------------------------- */
{
  olong i, j, k, ax;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar *input_file="IDIOut.in", *arg;
  gboolean init=FALSE;
  oint itemp;
  gchar *strTemp;
  ObitInfoList* list;
  gchar *routine = "IDIOutin";

  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* command line arguments */
  if (argc<=1) Usage(); /* must have arguments */
  /* parse command line */
  for (ax=1; ax<argc; ax++) {
    arg = argv[ax];
    if (strcmp(arg, "-input") == 0){ /* input parameters */
      input_file = argv[++ax];
      /* parse input file */
      ObitParserParse (input_file, list, err);
      init = TRUE;

    } else if (strcmp(arg, "-output") == 0){ /* output results */
      outfile = argv[++ax];

    } else if (strcmp(arg, "-pgmNumber") == 0) { /*Program number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "pgmNumber", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-AIPSuser") == 0) { /* AIPS user number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "AIPSuser", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-inDisk") == 0) { /* Input disk */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "inDisk", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-inFile") == 0){ /* input FITS File name */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inFile", OBIT_string, dim, strTemp);

     } else if (strcmp(arg, "-inName") == 0) { /* AIPS UV inName */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inName", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-inClass") == 0) { /* AIPS UV inClass */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inClass", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-inSeq") == 0) { /* AIPS input UV sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "inSeq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-outFile") == 0){ /* Output FITS file */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outFile", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outDisk") == 0) { /* Output disk */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outDisk", OBIT_oint, dim, &itemp, err);
      
    } else { /* unknown argument */
      /* DEBUG fprintf (stderr,"DEBUG parameter %s \n",arg);*/
      Usage();
    }
  } /* end parsing input arguments */
  
  /* Read defaults if no file specified */
  if (!init) ObitParserParse (input_file, list, err);

  /* Extract basic information to program globals */
  ObitInfoListGet(list, "pgmNumber", &type, dim, &pgmNumber, err);
  ObitInfoListGet(list, "AIPSuser",  &type, dim, &AIPSuser,  err);
  ObitInfoListGet(list, "nAIPS",     &type, dim, &nAIPS,     err);
  ObitInfoListGet(list, "nFITS",     &type, dim, &nFITS,     err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  /* Directories more complicated */
  ObitInfoListGetP(list, "AIPSdirs",  &type, dim, (gpointer)&strTemp);
  if (strTemp) {  /* Found? */
    AIPSdirs = g_malloc0(dim[1]*sizeof(gchar*));
    for (i=0; i<dim[1]; i++) {
      AIPSdirs[i] =  g_malloc0(dim[0]*sizeof(gchar));
      k = 0;
      for (j=0; j<dim[0]; j++) { /* Don't copy blanks */
	if (strTemp[j]!=' ') {AIPSdirs[i][k] = strTemp[j]; k++;}
      }
      AIPSdirs[i][k] = 0;
      strTemp += dim[0];
    }
  }

  ObitInfoListGetP(list, "FITSdirs",  &type, dim, (gpointer)&strTemp);
  if (strTemp)   {  /* Found? */
    FITSdirs = g_malloc0(dim[1]*sizeof(gchar*));
    for (i=0; i<dim[1]; i++) {
      FITSdirs[i] =  g_malloc0(dim[0]*sizeof(gchar));
      k = 0;
      for (j=0; j<dim[0]; j++) { /* Don't copy blanks */
	if (strTemp[j]!=' ') {FITSdirs[i][k] = strTemp[j]; k++;}
      }
      FITSdirs[i][k] = 0;
      strTemp += dim[0];
    }
  }

  /* Initialize output */
  myOutput = defaultOutputs(err);
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

 return list;
} /* end IDIOutin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: IDIOut -input file -output ofile [args]\n");
    fprintf(stderr, "Convert an Obit/UV to FITS/IDI file format\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def IDIOut.in\n");
    fprintf(stderr, "  -output output result parameter file, def IDIOut.out\n");
    fprintf(stderr, "  -inFile output uv FITS  file\n");  
    fprintf(stderr, "  -inName output uv AIPS file name\n");
    fprintf(stderr, "  -inClass output uv AIPS file class\n");
    fprintf(stderr, "  -inSeq output uv AIPS file sequence\n");
    fprintf(stderr, "  -inDisk output uv (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -outFile output FITS IDI file names\n");    
    /*/exit(1);  bail out */
  }/* end Usage */

/*----------------------------------------------------------------------- */
/*  Create default input ObitInfoList                                     */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*  Values:                                                               */
/*     pgmNumber Int        Program number (like POPS number) def 1       */
/*     nFITS     Int        Number of FITS directories [def. 1]           */
/*     FITSdirs  Str [?,?]  FITS directories [def {"./"}]                 */
/*     AIPSuser  Int        AIPS user number [def 2}]                     */
/*     nAIPS     Int        Number of AIPS directories [def. 1]           */
/*     AIPSdirs  Str [?,?]  AIPS directories [def {"AIPSdata/"}]          */
/*     DataType  Str [4]    "AIPS" or "FITS" [def {"FITS"}]               */
/*     doCalSelect Boo      TRUE                                          */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  gchar *scan_name ="unspecified";
  ObitInfoList *out = newObitInfoList();
  gboolean doCalSelect;
  gchar *routine = "defaultInputs";

  /* add parser items */
  /* Program number */
  dim[0] = 1; dim[1] = 1;
  itemp = 1;
  ObitInfoListPut (out, "pgmNumber", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Default FITS directories - same directory */
  dim[0] = 1; dim[1] = 1;
  itemp = 2; /* number of FITS directories */
  ObitInfoListPut (out, "nFITS", OBIT_oint, dim, &itemp, err);

  /* AIPS user number */
  dim[0] = 1; dim[1] = 1;
  itemp = 2;
  ObitInfoListPut (out, "AIPSuser", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Default AIPS directories */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; /* number of AIPS directories */
  ObitInfoListPut (out, "nAIPS", OBIT_oint, dim, &itemp, err);

  /* Default type "FITS" */
  strTemp = "FITS";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "DataType", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* base of scan file names */
  dim[0] = strlen (scan_name);
  ObitInfoListPut (out, "inFile", OBIT_string, dim, scan_name, err);

  /* output FITS file name */
  strTemp = "uvData.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);

  /* doCalSelect */
  doCalSelect = TRUE;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (out, "doCalSelect", OBIT_bool, dim, &doCalSelect);

  return out;
} /* end defaultInputs */

/*----------------------------------------------------------------------- */
/*  Create default output ObitInfoList                                    */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultOutputs(ObitErr *err)
{
  /*gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};*/
  /*ofloat ftemp;*/
  ObitInfoList *out = newObitInfoList();
  /*gchar *routine = "defaultOutputs";*/

  /* No outputs */
  return out;
} /* end defaultOutputs */

/*----------------------------------------------------------------------- */
/*  Create input uv data                                                  */
/*   Input:                                                               */
/*      Source    Source name                                             */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/* Returns the output uv data                                             */
/*----------------------------------------------------------------------- */
ObitUV* setInputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitUV    *inUV = NULL;
  ObitInfoType type;
  olong      i, n, Aseq, disk, cno, lType;
  gchar     *Type, *strTemp, inFile[129];
  gchar     Aname[13], Aclass[7], *Atype = "UV";
  olong      nvis;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean  exist;
  gchar     tname[129], *fullname=NULL;
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "souCode", "Qual", "Stokes", "timeRange", 
    "BChan", "EChan", "BIF", "EIF", "FreqID",
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer", "doPol",
    "Smooth", "Antennas",  "subA", "Sources", "souCode", "Qual",
     NULL};
  gchar     *routine = "setInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inUV;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input UV Object */
  g_snprintf (tname, 100, "input UV data");
  inUV = newObitUV(tname);
    
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  lType = dim[0];
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */

    /* inName given? */
    ObitInfoListGetP (myInput, "inName", &type, dim, (gpointer)&strTemp);
    for (i=0; i<12; i++) Aname[i] = ' ';  Aname[i] = 0;
    for (i=0; i<MIN(12,dim[0]); i++) Aname[i] = strTemp[i];
    /* Save any defaulting on myInput */
    dim[0] = 12;
    ObitInfoListAlwaysPut (myInput, "inName", OBIT_string, dim, Aname);
      
    /* input AIPS class */
    if (ObitInfoListGetP(myInput, "inClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    /* Default in class is "IDIOut" */
    if (!strncmp(Aclass, "      ", 6)) strncpy (Aclass, "IDIOut", 7);

    /* input AIPS disk - default is inDisk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    if (disk<=0)
       ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "inSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", inUV);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "inSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Allocate catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inUV);

    /* define object */
    nvis = 1;
    ObitUVSetAIPS (inUV, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inUV);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */

    /* inFile given? */
    ObitInfoListGetP (myInput, "inFile", &type, dim, (gpointer)&strTemp);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) inFile[i] = strTemp[i]; inFile[i] = 0;
    ObitTrimTrail(inFile);  /* remove trailing blanks */

    /* Save any defaulting on myInput */
    dim[0] = strlen(inFile);
    ObitInfoListAlwaysPut (myInput, "inFile", OBIT_string, dim, inFile);

    /* input FITS disk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    if (disk<=0) /* defaults to inDisk */
      ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);

    /* Did it previously exist */
    fullname = ObitFITSFilename (disk, inFile, err);
    exist =  ObitFileExist (fullname, err);
    if (fullname) g_free(fullname);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inUV);
    
    /* define object */
    nvis = 1;
    ObitUVSetFITS (inUV, nvis, disk, inFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inUV);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return inUV;
  }
  
  /* Fully instantiate */
  ObitUVFullInstantiate (inUV, TRUE, err);

  /* ANtenna information */
  ObitUVGetSubA (inUV, err);
  maxAnt = inUV->myDesc->maxAnt;
  if (err->error) Obit_traceback_val (err, routine, "myInput", inUV);

  /* Get input parameters from myInput, copy to inUV */
  ObitInfoListCopyList (myInput, inUV->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inUV);

  return inUV;
} /* end setInputUV */

void InitIDI (gchar *FITSfile, olong disk, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Initialize output IDI header                                          */
/*   Input:                                                               */
/*      FITSfile  Name of output FITS file                                */
/*      disk      FITS disk number                                        */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitFileFITS* myFile=NULL;
  long naxes[2] = {0, 0};
  gchar *FullName=NULL;
  int status = 0;
  gchar *routine = "InitIDI";

  /* error checks */
  if (err->error) return;
  g_assert(FITSfile!=NULL);

  /* Get full file name - without '!' */
  if (FITSfile[0]=='!')
    FullName = ObitFITSFilename (disk, &FITSfile[1], err);
  else
    FullName = ObitFITSFilename (disk, FITSfile, err);

  /* If file exists and starts with '!' - Zap first */
  if (ObitFileExist(&FITSfile[1], err) && (FITSfile[0]=='!')) {
    ObitFileZapFile (FullName, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, "Output");

  myFile =  newObitFileFITS ("Output");
  ObitFileFITSOpen (myFile, FullName, 0, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, "Output");

  /* Dummy header */
  fits_write_imghdr (myFile->myFptr, 8, 0, naxes, &status);
  if (status) {             /* it went wrong */
    Obit_log_error(err, OBIT_Error, "%s: ERROR creating FITS file %s", 
		   routine, "Output");
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
  }

  /* Close */
  ObitFileFITSClose (myFile, err);
  if (err->error) Obit_traceback_msg (err, routine, "Output");		     
  
  myFile = ObitFileFITSUnref(myFile); /* cleanup */
  if (FullName) g_free(FullName);

  return;
} /* end InitIDI */

void PutFrequencyInfo (ObitUV *inData, ObitData *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Put Frequency info to IDI_FREQUENCY table                             */
/*  Write FREQUENCY table                                                 */
/*   Input:                                                               */
/*      inData   Input UV object                                          */
/*      outData  Output IDI FITS object                                   */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_FREQUENCY*    outTable=NULL;
  ObitTableIDI_FREQUENCYRow* outRow=NULL;
  ObitTableFQ*            inTable=NULL;
  ObitTableFQRow*         inRow=NULL;
  olong i, lim, iRow, oRow, ver;
  oint numIF;
  ObitIOAccess access;
  gchar *routine = "PutFrequencyInfo";
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(outData));
  g_assert (ObitUVIsA(inData));
  
  /* Create input FQ table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  numIF = 0;
  inTable = newObitTableFQValue ("Output table", (ObitData*)inData, 
				 &ver, access, numIF, err);
  if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with FQ table");
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Open table */
  if ((ObitTableFQOpen (inTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FQ table");
    return;
  }
  
  /* Create output FREQUENCY table object */
  ver = 1;
  access = OBIT_IO_WriteOnly;
  numIF = inTable->numIF;
  outTable = newObitTableIDI_FREQUENCYValue ("Input table", outData, 
					     &ver, access, numIF, err);
  if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_FREQUENCY table");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Open table */
  if ((ObitTableIDI_FREQUENCYOpen (outTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_FREQUENCY table");
    return;
  }
  
  /* Create output Row */
  outRow = newObitTableIDI_FREQUENCYRow (outTable);
  /* attach to table buffer */
  ObitTableIDI_FREQUENCYSetRow (outTable, outRow, err);
  
  /* Create input Row */
  inRow = newObitTableFQRow (inTable);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Set FREQUENCY table values */
  outTable->no_stkd  = inData->myDesc->inaxes[inData->myDesc->jlocs];
  outTable->stk_1    = (olong)inData->myDesc->crval[inData->myDesc->jlocs];
  outTable->no_band  = inData->myDesc->inaxes[inData->myDesc->jlocif];
  outTable->no_chan  = inData->myDesc->inaxes[inData->myDesc->jlocf];
  outTable->ref_freq = inData->myDesc->crval[inData->myDesc->jlocf];
  outTable->chan_bw  = inData->myDesc->cdelt[inData->myDesc->jlocf];
  outTable->ref_pixl = inData->myDesc->crpix[inData->myDesc->jlocf];
  lim = MIN (UVLEN_VALUE,MAXKEYCHARTABLEIDI_FREQUENCY);
  strncpy (outTable->obscode,   inData->myDesc->observer, 8);
  lim = MIN (MAXKEYCHARTABLEIDI_FREQUENCY, UVLEN_VALUE);
  strncpy (outTable->RefDate, inData->myDesc->obsdat, lim);
  outTable->myStatus = OBIT_Modified; /* Mark as modified */
  
  /* Initialize output row */
  outRow->fqid    = 0;
  for (i=0; i<numIF; i++) {
    outRow->bandfreq[i] = 0.0;
    outRow->chWidth[i]  = 0;
    outRow->totBW[i]    = 0;
    outRow->sideBand[i] = 0;
  }
  outRow->status    = 0;

  /* loop through input table */
  for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
    if ((ObitTableFQReadRow (inTable, iRow, inRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
     Obit_log_error(err, OBIT_Error, "ERROR reading FG Table");
      return;
    }

    /* Save to FREQUENCY table */
    outRow->fqid    = inRow->fqid;
    for (i=0; i<numIF; i++) {
      outRow->bandfreq[i] = inRow->freqOff[i];
      outRow->chWidth[i]  = inRow->chWidth[i];
      outRow->totBW[i]    = inRow->totBW[i];
      outRow->sideBand[i] = inRow->sideBand[i];
    }
    outRow->status    = 0;
    oRow = -1;
    if ((ObitTableIDI_FREQUENCYWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
       Obit_log_error(err, OBIT_Error, "ERROR writing IDI_FREQUENCY Table");
      return;
    }
  } /* end loop over input table */
  
  /* Close  tables */
  if ((ObitTableIDI_FREQUENCYClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_FREQUENCY Table file");
    return;
  }

  if ((ObitTableFQClose (inTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output FQ Table file");
    return;
  }

  /* Cleanup */
  outRow    = ObitTableIDI_FREQUENCYRowUnref(outRow);
  outTable  = ObitTableIDI_FREQUENCYUnref(inTable);
  inRow     = ObitTableFQRowUnref(outRow);
  inTable   = ObitTableFQUnref(outTable);

 } /* end  PutFrequencyInfo */

void PutAntennaInfo (ObitInfoList *myInput, ObitUV *inData, 
		     ObitData *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Writes IDI_ANTENNA and IDI_ARRAY_GEOMETRY tables                      */
/*  Does not copy poln cal. info if being applied                         */
/*   Input:                                                               */
/*      myInput  Input parser object                                      */
/*      inData   Input UV object                                          */
/*      outData  Output IDI FITS object                                   */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_ANTENNA         *outTable=NULL;
  ObitTableIDI_ANTENNARow      *outRow=NULL;
  ObitTableIDI_ARRAY_GEOMETRY  *out2Table=NULL;
  ObitTableIDI_ARRAY_GEOMETRYRow  *out2Row=NULL;
  ObitTableAN                  *inTable=NULL;
  ObitTableANRow               *inRow=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, j, lim, iRow, oRow, ver;
  olong nANver, iANver;
  oint numIF, numPCal, numOrb;
  ObitIOAccess access;
  gboolean isVLA, doPol;
  double VLALong, cosLong, sinLong;
  gchar *routine = "PutAntennaInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(outData));
  g_assert (ObitUVIsA(inData));

  /* How many AN tables are there? */
  nANver = ObitTableListGetHigh (inData->tableList, "AIPS AN");

  /* Is poln cal being applied? */
  doPol = -1;
  ObitInfoListGetTest(myInput, "doPol", &type, dim, &doPol);

  /* Array Geometry - loop over table */
  for (iANver=1; iANver<=nANver; iANver++) {
    /* Create input Antenna table object */
    ver = iANver;
    access  = OBIT_IO_ReadWrite;
    numOrb  = 0;
    numPCal = 0;
    inTable = newObitTableANValue ("Input table", (ObitData*)inData, 
				   &ver, access, numOrb, numPCal, err);
    if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableANOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output AN table");
      return;
    }

    /* Is this the VLA? */
    isVLA  = !strncmp(inTable->ArrName, "VLA     ", 8);
    VLALong = 1.878283678;
    cosLong = cos(VLALong);
    sinLong = sin(VLALong);

    /* Create Row */
    inRow = newObitTableANRow (inTable);
    
    /* Create output Array geometry table object */
    ver    = iANver;
    access = OBIT_IO_WriteOnly;
    numIF  = inData->myDesc->inaxes[inData->myDesc->jlocif];
    numOrb = inTable->numOrb;
    out2Table = newObitTableIDI_ARRAY_GEOMETRYValue ("Output table", outData, 
						     &ver, access, numIF, numOrb, err);
    if (out2Table==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_ARRAY_GEOMETRY table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableIDI_ARRAY_GEOMETRYOpen (out2Table, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output IDI_ARRAY_GEOMETRY table");
      return;
    }
    
    /* Create Row */
    out2Row = newObitTableIDI_ARRAY_GEOMETRYRow (out2Table);
    /* attach to table buffer */
    ObitTableIDI_ARRAY_GEOMETRYSetRow (out2Table, out2Row, err);
    
    /* Set ARRAY_GEOMETRY table values */
    out2Table->no_stkd  = inData->myDesc->inaxes[inData->myDesc->jlocs];
    out2Table->stk_1    = (olong)inData->myDesc->crval[inData->myDesc->jlocs];
    out2Table->no_band  = inData->myDesc->inaxes[inData->myDesc->jlocif];
    out2Table->no_chan  = inData->myDesc->inaxes[inData->myDesc->jlocf];
    out2Table->ref_freq = inData->myDesc->crval[inData->myDesc->jlocf];
    out2Table->chan_bw  = inData->myDesc->cdelt[inData->myDesc->jlocf];
    out2Table->ref_pixl = inData->myDesc->crpix[inData->myDesc->jlocf];
    out2Table->numOrb   = inTable->numOrb;
    out2Table->Freq     = inTable->Freq;
    out2Table->GSTiat0  = inTable->GSTiat0;
    out2Table->DegDay   = inTable->DegDay;
    out2Table->ut1Utc   = inTable->ut1Utc;
    out2Table->iatUtc   = inTable->iatUtc;
    out2Table->PolarX   = inTable->PolarX;
    out2Table->PolarY   = inTable->PolarY;
    out2Table->ArrayX   = inTable->ArrayX;
    out2Table->ArrayY   = inTable->ArrayY;
    out2Table->ArrayZ   = inTable->ArrayZ;
    lim = MIN (UVLEN_VALUE,MAXKEYCHARTABLEIDI_ARRAY_GEOMETRY);
    strncpy (out2Table->obscode, inData->myDesc->observer, 8);
    lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, MAXKEYCHARTABLEAN);
    strncpy (out2Table->RefDate, inTable->RefDate, lim);
    strncpy (out2Table->TimeSys, inTable->TimeSys, lim);
    strncpy (out2Table->ArrName, inTable->ArrName, lim);
    strncpy (out2Table->frame,   "GEOCENTRIC", lim);
    out2Table->myStatus = OBIT_Modified; /* Mark as modified */
    
    /* Initialize output row */
    out2Row->noSta     = 0;
    out2Row->mntSta    = 0;
    out2Row->diameter  = 0.0;
    out2Row->AntName[0]= 0;
    out2Row->StaXYZ[0] = 0.0;
    out2Row->StaXYZ[1] = 0.0;
    out2Row->StaXYZ[2] = 0.0;
    out2Row->derXYZ[0] = 0.0;
    out2Row->derXYZ[1] = 0.0;
    out2Row->derXYZ[2] = 0.0;
    out2Row->staXof[0] = 0.0;
    out2Row->staXof[1] = 0.0;
    out2Row->staXof[2] = 0.0;
    for (i=0; i<numOrb; i++) out2Row->OrbParm[0] = 0.0;
    out2Row->status      = 0;

    
    /* loop through input Antenna table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableANReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading AN Table");
	return;
      }

      /* If this is the VLA rotate to IERS system */
      if (isVLA) {
      /* ****** CHECK THIS  ***** */
	out2Row->StaXYZ[0] = cosLong*inRow->StaXYZ[0] - sinLong*inRow->StaXYZ[1];
	out2Row->StaXYZ[1] = cosLong*inRow->StaXYZ[1] + sinLong*inRow->StaXYZ[0];
	out2Row->StaXYZ[2] = inRow->StaXYZ[2];
      } else {
	/* Not VLA */
	out2Row->StaXYZ[0] = inRow->StaXYZ[0];
	out2Row->StaXYZ[1] = inRow->StaXYZ[1];
	out2Row->StaXYZ[2] = inRow->StaXYZ[2];
     }
      
      /* Set output Row */
      out2Row->noSta     = inRow->noSta;
      out2Row->mntSta    = inRow->mntSta;
      out2Row->staXof[0] = inRow->staXof;
      for (i=0; i<numOrb; i++) out2Row->OrbParm[0]= inRow->OrbParm[i];
      for (i=0; i<8;      i++) out2Row->AntName[i]= inRow->AntName[i];
      

      oRow = -1;
      if ((ObitTableIDI_ARRAY_GEOMETRYWriteRow (out2Table, oRow, out2Row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating ARRAY_GEOMETRY Table");
	return;
      }
    } /* end loop over input table */
    
    if ((ObitTableIDI_ARRAY_GEOMETRYClose (out2Table, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_ARRAY_GEOMETRY Table file");
      return;
    }
    
    if ((ObitTableANClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output Antenna Table file");
      return;
    }
    
    /* Cleanup */
    inRow     = ObitTableANRowUnref(inRow);
    inTable   = ObitTableANUnref(inTable);
    out2Row   = ObitTableIDI_ARRAY_GEOMETRYRowUnref(out2Row);
    out2Table = ObitTableIDI_ARRAY_GEOMETRYUnref(out2Table);
  } /* end loop converting to ARRAY_GEOMETRY tables */

  
  /* ANTENNA  - loop over table */
  for (iANver=1; iANver<=nANver; iANver++) {
    /* Create input Antenna table object */
    ver = iANver;
    access  = OBIT_IO_ReadWrite;
    numOrb  = 0;
    numPCal = 0;
    inTable = newObitTableANValue ("Input table", (ObitData*)inData, 
				   &ver, access, numOrb, numPCal, err);
    if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableANOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output AN table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableANRow (inTable);
    
    /* All output to same table - create/open,/init once */
    if (outTable==NULL) {
      /* Create output Array geometry table object */
      ver      = 1;
      access   = OBIT_IO_WriteOnly;
      numIF    = inData->myDesc->inaxes[inData->myDesc->jlocif];
      numPCal  = inTable->numPCal/numIF;  /*in AN table this is for all bands */
      if (doPol) numPCal = 0;
      outTable = newObitTableIDI_ANTENNAValue ("Output table", outData, 
					       &ver, access, numIF, numPCal, err);
      if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_ARRAY_GEOMETRY table");
      if (err->error) Obit_traceback_msg (err, routine, outData->name);
      
      /* Open table */
      if ((ObitTableIDI_ANTENNAOpen (outTable, access, err) 
	   != OBIT_IO_OK) || (err->error))  { /* error test */
	Obit_log_error(err, OBIT_Error, "ERROR opening output IDI_ARRAY_GEOMETRY table");
	return;
      }
      
      /* Create Row */
      outRow = newObitTableIDI_ANTENNARow (outTable);
      /* attach to table buffer */
      ObitTableIDI_ANTENNASetRow (outTable, outRow, err);
      
      /* Set ANTENNA table values */
      outTable->no_stkd  = inData->myDesc->inaxes[inData->myDesc->jlocs];
      outTable->stk_1    = (olong)inData->myDesc->crval[inData->myDesc->jlocs];
      outTable->no_band  = inData->myDesc->inaxes[inData->myDesc->jlocif];
      outTable->no_chan  = inData->myDesc->inaxes[inData->myDesc->jlocf];
      outTable->ref_freq = inData->myDesc->crval[inData->myDesc->jlocf];
      outTable->chan_bw  = inData->myDesc->cdelt[inData->myDesc->jlocf];
      outTable->ref_pixl = inData->myDesc->crpix[inData->myDesc->jlocf];
      lim = MIN (UVLEN_VALUE,MAXKEYCHARTABLEIDI_ARRAY_GEOMETRY);
      strncpy (outTable->obscode,   inData->myDesc->observer, 8);
      lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, MAXKEYCHARTABLEAN);
      strncpy (outTable->RefDate, inTable->RefDate, lim);
      strncpy (outTable->ArrName, inTable->ArrName, lim);
      if (doPol)  /* Applying cal? */
	strncpy (outTable->polType, "        ", lim);
      else
	strncpy (outTable->polType, inTable->polType, lim);
      outTable->myStatus = OBIT_Modified; /* Mark as modified */
      
      /* Initialize output row */
      outRow->time          = 0.0;
      outRow->time_interval = 1.0e20;
      outRow->antenna_no    = 0;
      outRow->array         = 0;
      outRow->freqid        = 0;
      outRow->no_levels     = 0;
      outRow->polTypeA      = 'R';
      outRow->polTypeB      = 'L';
      outRow->AntName[0]    = 0;
      for (i=0; i<numIF; i++) {
	outRow->PolAngA[i]  = 0.0;
	outRow->PolAngB[i]  = 0.0;
	outRow->PolCalA[i]  = 0.0;
	outRow->PolCalB[i]  = 0.0;
	outRow->BeamFWHM[i] = 0.0;
      }
      outRow->status      = 0;
    } /* end setup output first pass */
    
    
    /* loop through input Antenna table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableANReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading AN Table");
	return;
      }
      
      /* Set output Row */
      outRow->antenna_no    = inRow->noSta;
      outRow->array         = iANver;
      outRow->freqid        = 1; /* SHIT */
      outRow->polTypeA      = inRow->polTypeA[0];
      outRow->polTypeB      = inRow->polTypeB[0];
      for (i=0; i<8; i++) outRow->AntName[i]= inRow->AntName[i];
      for (i=0; i<numIF; i++) {
	outRow->PolAngA[i]  = inRow->PolAngA;
	outRow->PolAngB[i]  = inRow->PolAngB;
	/* numPCal is per band, in AIPS AN it is total */
	if (!doPol) {  /* Don't copy if applying */
	  for (j=0; j<numPCal; j++) {
	    outRow->PolCalA[j+i*numPCal]  = inRow->PolCalA[j+i*numPCal];
	    outRow->PolCalB[j+i*numPCal]  = inRow->PolCalB[j+i*numPCal];
	  }
	}
      }

      oRow = -1;
      if ((ObitTableIDI_ANTENNAWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating ANTENNA Table");
	return;
      }
    } /* end loop over input table */
    
    if ((ObitTableANClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output Antenna Table file");
      return;
    }
    
    /* Cleanup */
    inRow     = ObitTableANRowUnref(inRow);
    inTable   = ObitTableANUnref(inTable);
  } /* end loop converting to ANTENNA table */
  
  /* Close ANTENNA table */
  if ((ObitTableIDI_ANTENNAClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_ARRAY_GEOMETRY Table file");
    return;
  }
  
  /* Cleanup */
  outRow   = ObitTableIDI_ANTENNAUnref(outRow);
  outTable = ObitTableIDI_ANTENNAUnref(outTable);
} /* end  PutAntennaInfo */

void PutSourceInfo (ObitUV *inData, ObitData *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Put Source info to IDI_SOURCE table                                   */
/*   Input:                                                               */
/*      inData   Input UV object                                          */
/*      outData  Output IDI FITS object                                   */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/* ---------------------------------------------------------------------- */
{
  ObitTableSU*            inTable=NULL;
  ObitTableSURow*         inRow=NULL;
  ObitTableIDI_SOURCE*    outTable=NULL;
  ObitTableIDI_SOURCERow* outRow=NULL;
  olong i, lim, iRow, oRow, ver;
  oint numIF;
  ObitIOAccess access;
  gchar *routine = "PutSourceInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(outData));
  g_assert (ObitUVIsA(inData));

  /* Don't bother if no sourceID in data */
  if (inData->myDesc->ilocsu<0) return;

  /* Create input Source table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  numIF = 0;
  inTable = newObitTableSUValue ("Input table", (ObitData*)inData, 
				 &ver, access, numIF, err);
  if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with SU table");
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Open table */
  if ((ObitTableSUOpen (inTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input SU table");
    return;
  }

  /* Create output Source table object */
  ver = 1;
  access = OBIT_IO_WriteOnly;
  numIF = inTable->numIF;
  outTable = newObitTableIDI_SOURCEValue ("Output table", outData, 
					  &ver, access, numIF, err);
  if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_SOURCE table");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open table */
  if ((ObitTableIDI_SOURCEOpen (outTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
	Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_SOURCE table");
	return;
      }
  
  /* Create output Row */
  outRow = newObitTableIDI_SOURCERow (outTable);
  /* attach to table buffer */
  ObitTableIDI_SOURCESetRow (outTable, outRow, err);

  /* Create input Row */
  inRow = newObitTableSURow (inTable);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Set SOURCE table values */
  outTable->no_stkd  = inData->myDesc->inaxes[inData->myDesc->jlocs];
  outTable->stk_1    = (olong)inData->myDesc->crval[inData->myDesc->jlocs];
  outTable->no_band  = inData->myDesc->inaxes[inData->myDesc->jlocif];
  outTable->no_chan  = inData->myDesc->inaxes[inData->myDesc->jlocf];
  outTable->ref_freq = inData->myDesc->crval[inData->myDesc->jlocf];
  outTable->chan_bw  = inData->myDesc->cdelt[inData->myDesc->jlocf];
  outTable->ref_pixl = inData->myDesc->crpix[inData->myDesc->jlocf];
  lim = MIN (UVLEN_VALUE,MAXKEYCHARTABLEIDI_ARRAY_GEOMETRY);
  strncpy (outTable->obscode,   inData->myDesc->observer, 8);
  lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, UVLEN_VALUE);
  strncpy (outTable->RefDate, inData->myDesc->obsdat, lim);
  outTable->myStatus = OBIT_Modified; /* Mark as modified */

  /* Initialize output row */
  outRow->SourID    = 0;
  outRow->Qual      = 0;
  outRow->RAMean    = 0.0;
  outRow->DecMean   = 0.0;
  outRow->Epoch     = 0.0;
  outRow->RAApp     = 0.0;
  outRow->DecApp    = 0.0;
  outRow->PMRa      = 0.0;
  outRow->PMDec     = 0.0;
  outRow->Source[0] = 0;
  outRow->CalCode[0]= 0;
  for (i=0; i<inTable->numIF; i++) {
    outRow->IFlux[i]     = 0.0;
    outRow->QFlux[i]     = 0.0;
    outRow->UFlux[i]     = 0.0;
    outRow->VFlux[i]     = 0.0;
    outRow->FreqOff[i]   = 0.0;
    outRow->SysVel[i]    = 0.0;
    outRow->RestFreq[i]  = 0.0;
    outRow->alpha[i]     = 0.0;
  }
  outRow->status    = 0;

  /* loop through input table */
  for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
    if ((ObitTableSUReadRow (inTable, iRow, inRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR reading Source Table");
      return;
    }
    /* add to output */
    /* If first time update header */
    /* Set output row for end of table */
    outRow->SourID    = outTable->myDesc->nrow+1;
    outRow->Qual      = inRow->Qual;
    outRow->RAMean    = inRow->RAMean;
    outRow->DecMean   = inRow->DecMean;
    outRow->Epoch     = 0.0;
    outRow->RAApp     = inRow->RAApp;
    outRow->DecApp    = inRow->DecApp;
    outRow->PMRa      = inRow->PMRa;
    outRow->PMDec     = inRow->PMDec;
    strncpy(outRow->VelTyp, inTable->velType, 8);
    strncpy(outRow->VelDef, inTable->velDef, 8);
    if (inRow->Epoch>1975.) strncpy (outRow->Equinox, "J2000",8);
    else  strncpy (outRow->Equinox, "B1950",8);
    lim = MIN(inTable->myDesc->repeat[inTable->SourceCol], 
	      outTable->myDesc->repeat[outTable->SourceCol]);
    for (i=0; i<lim; i++) 
      outRow->Source[i] = inRow->Source[i];
    lim = MIN(inTable->myDesc->repeat[inTable->CalCodeCol], 
	      outTable->myDesc->repeat[outTable->CalCodeCol]);
    for (i=0; i<lim; i++) 
      outRow->CalCode[i]= inRow->CalCode[i];
    for (i=0; i<inTable->numIF; i++) {
      outRow->IFlux[i]    = inRow->IFlux[i];
      outRow->QFlux[i]    = inRow->QFlux[i];
      outRow->UFlux[i]    = inRow->UFlux[i];
      outRow->VFlux[i]    = inRow->VFlux[i];
      outRow->FreqOff[i]  = inRow->FreqOff[i];
      outRow->SysVel[i]   = inRow->LSRVel[i];
      outRow->RestFreq[i] = inRow->RestFreq[i];
    }

    oRow = -1;
    if ((ObitTableIDI_SOURCEWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating IDI_SOURCE Table");
      return;
    }
  } /* end loop over input table */
  
  /* Close  tables */
  if ((ObitTableIDI_SOURCEClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output IDI_SOURCE Table file");
    return;
  }

  if ((ObitTableSUClose (inTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input Source Table file");
    return;
  }

  /* Cleanup */
  outRow   = ObitTableIDI_SOURCERowUnref(outRow);
  outTable = ObitTableIDI_SOURCEUnref(outTable);
  inRow    = ObitTableSURowUnref(inRow);
  inTable  = ObitTableSUUnref(inTable);

} /* end  PutSourceInfo */

void PutFlagInfo (ObitInfoList *myInput, ObitUV *inData, ObitData *outData, 
		  ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any FLAG tables on inData to AIPS FG on outData               */
/*  Nothing is copied if flagging is being applied                        */
/*   Input:                                                               */
/*      myInput  Input parser object                                      */
/*      inData   Input UV object                                          */
/*      outData  Output IDI FITS object                                   */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableFG*          inTable=NULL;
  ObitTableFGRow*       inRow=NULL;
  ObitTableIDI_FLAG*    outTable=NULL;
  ObitTableIDI_FLAGRow* outRow=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, iver, hiVer, lim, iRow, oRow, ver;
  oint numIF, flagVer;
  ObitIOAccess access;
  gchar *routine = "PutFlagInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(outData));
  g_assert (ObitUVIsA(inData));

  /* Is flagging being applied? */
  flagVer = -1;
  ObitInfoListGetTest(myInput, "flagVer", &type, dim, &flagVer);
  if (flagVer>=0) return;

  /* Loop over plausible versions */
  hiVer = ObitTableListGetHigh (inData->tableList, "AIPS FG");
  for (iver=1; iver<=hiVer; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input Flag table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    inTable = newObitTableFGValue ("Input table", (ObitData*)inData, 
				   &ver, access, err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      continue;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    numIF = inData->myDesc->inaxes[inData->myDesc->jlocif];
    
    /* Open table */
    if ((ObitTableFGOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input FG table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableFGRow (inTable);
    
    /* Create output IDI_FLAG table object */
    ver = iver;
    access = OBIT_IO_WriteOnly;
    numIF = inData->myDesc->inaxes[inData->myDesc->jlocif];
    outTable = newObitTableIDI_FLAGValue ("Output table", outData, 
					  &ver, access, numIF, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with FG table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableIDI_FLAGOpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output FG table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTableIDI_FLAGRow (outTable);
    /* attach to table buffer */
    ObitTableIDI_FLAGSetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Set FLAG table values */
    outTable->no_stkd  = inData->myDesc->inaxes[inData->myDesc->jlocs];
    outTable->stk_1    = (olong)inData->myDesc->crval[inData->myDesc->jlocs];
    outTable->no_band  = inData->myDesc->inaxes[inData->myDesc->jlocif];
    outTable->no_chan  = inData->myDesc->inaxes[inData->myDesc->jlocf];
    outTable->ref_freq = inData->myDesc->crval[inData->myDesc->jlocf];
    outTable->chan_bw  = inData->myDesc->cdelt[inData->myDesc->jlocf];
    outTable->ref_pixl = inData->myDesc->crpix[inData->myDesc->jlocf];
    lim = MIN (UVLEN_VALUE,MAXKEYCHARTABLEIDI_ARRAY_GEOMETRY);
    strncpy (outTable->obscode,   inData->myDesc->observer, 8);
    lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, UVLEN_VALUE);
    strncpy (outTable->RefDate, inData->myDesc->obsdat, lim);
    outTable->myStatus = OBIT_Modified; /* Mark as modified */
    
    /* Initialize output row */
    outRow->status = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableFGReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading FG Table");
	return;
      }
      
      /* Loop over bands - write one at a time in FLAG table */
      outRow->SourID    = inRow->SourID;
      outRow->Array     = inRow->SubA;
      outRow->fqid      = inRow->freqID;
      outRow->ants[0]   = inRow->ants[0];
      outRow->ants[1]   = inRow->ants[1];
      outRow->timerange[0] = inRow->TimeRange[0];
      outRow->timerange[1] = inRow->TimeRange[1];
      outRow->chans[0]  = inRow->chans[0];
      outRow->chans[1]  = inRow->chans[1];
      outRow->pflags[0] = inRow->pFlags[0];
      outRow->pflags[1] = inRow->pFlags[1];
      outRow->pflags[2] = inRow->pFlags[2];
      outRow->pflags[3] = inRow->pFlags[3];
      strncpy (outRow->reason, inRow->reason, 24);
      for (i=0; i<numIF; i++) outRow->bands[i] = 0;
      for (i=inRow->ifs[0]; i<=inRow->ifs[i]; i++) 
	outRow->bands[i] = 1;
      
      /* Write */
      oRow = -1;
      if ((ObitTableIDI_FLAGWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating FG Table");
	return;
      }
    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_FLAGClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output IDI_FLAG Table file");
      return;
    }
    
    if ((ObitTableFGClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input FG Table file");
      return;
    }
    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied FLAG table %d", iver);

    /* Cleanup */
    inRow     = ObitTableFGRowUnref(inRow);
    inTable   = ObitTableFGUnref(inTable);
    outRow    = ObitTableIDI_FLAGRowUnref(outRow);
    outTable  = ObitTableIDI_FLAGUnref(outTable);
  } /* end loop over versions */


} /* end  PutFlagInfo */

void PutCalibrationInfo (ObitInfoList *myInput, ObitUV *inData, 
			 ObitData *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any CALIBRATION tables on inData to AIPS CL on outData        */
/*  Nothing is copied if calibration is being applied                     */
/*   Input:                                                               */
/*      myInput  parser object                                            */
/*      inData   Input UV object                                          */
/*      outData  Output IDI FITS object                                   */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableCL*          inTable=NULL;
  ObitTableCLRow*       inRow=NULL;
  ObitTableIDI_CALIBRATION*    outTable=NULL;
  ObitTableIDI_CALIBRATIONRow* outRow=NULL;
  olong i, iver, hiVer, lim, iRow, oRow, ver, doCalib;
  oint numIF, numAnt, numPol, numTerm;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOAccess access;
  gchar *routine = "PutCalibrationInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(outData));
  g_assert (ObitUVIsA(inData));

  /* Is calibration being applied? */
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doCalib", &type, dim, &doCalib);
  if (doCalib>0) return;

  /* Loop over plausible versions */
  hiVer = ObitTableListGetHigh (inData->tableList, "AIPS CL");
  for (iver=1; iver<=hiVer; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input CL table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    numIF = numAnt = numPol = numTerm, 0;
    inTable = newObitTableCLValue ("Input table", (ObitData*)inData, 
				   &ver, access, numPol, numIF, numTerm, err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      continue;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableCLOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input CL table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableCLRow (inTable);
    
    /* Create output CALIBRATION table object */
    ver = iver;
    access  = OBIT_IO_WriteOnly;
    numIF   = inTable->numIF;
    numAnt  = maxAnt;
    numPol  = inTable->numPol;
    numTerm = 0;
    outTable = newObitTableIDI_CALIBRATIONValue ("Output table", outData, 
						 &ver, access, numIF, numAnt, numPol, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with CALIBRATION table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableIDI_CALIBRATIONOpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output IDI_CALIBRATION table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTableIDI_CALIBRATIONRow (outTable);
    /* attach to table buffer */
    ObitTableIDI_CALIBRATIONSetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Set CALIBRATION table values */
    outTable->no_stkd  = inData->myDesc->inaxes[inData->myDesc->jlocs];
    outTable->stk_1    = (olong)inData->myDesc->crval[inData->myDesc->jlocs];
    outTable->no_band  = inData->myDesc->inaxes[inData->myDesc->jlocif];
    outTable->no_chan  = inData->myDesc->inaxes[inData->myDesc->jlocf];
    outTable->ref_freq = inData->myDesc->crval[inData->myDesc->jlocf];
    outTable->chan_bw  = inData->myDesc->cdelt[inData->myDesc->jlocf];
    outTable->ref_pixl = inData->myDesc->crpix[inData->myDesc->jlocf];
    lim = MIN (UVLEN_VALUE,MAXKEYCHARTABLEIDI_ARRAY_GEOMETRY);
    strncpy (outTable->obscode,   inData->myDesc->observer, 8);
    lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, UVLEN_VALUE);
    strncpy (outTable->RefDate, inData->myDesc->obsdat, lim);
    outTable->myStatus = OBIT_Modified; /* Mark as modified */
    
    /* Initialize output row */
    for (i=0; i<numIF; i++) {
      outRow->TSys1[i] = 0.0;
      outRow->TAnt1[i] = 0.0;
      outRow->sensitivity1[i] = 0.0;
    }
    if (numPol>1) {   /* 2 poln */
      outRow->TSys2[i] = 0.0;
      outRow->TAnt2[i] = 0.0;
      outRow->sensitivity2[i] = 0.0;
    }
    outRow->status      = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableCLReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading CL Table");
	return;
      }
      
      /* Save to CL table */
      outRow->Time        = inRow->Time;
      outRow->TimeI       = inRow->TimeI;
      outRow->SourID      = inRow->SourID;
      outRow->antNo       = inRow->antNo;
      outRow->Array       = inRow->SubA;
      outRow->fqid        = inRow->FreqID;
      for (i=0; i<numIF; i++) {
	outRow->phase1[i]      = atan2 (inRow->Imag1[i], inRow->Real1[i]);
	outRow->real1[i]      = inRow->Real1[i];
	outRow->imag1[i]      = inRow->Imag1[i];
	outRow->rate1[i]      = inRow->Rate1[i];
	outRow->delay1[i]     = inRow->Delay1[i];
	outRow->weight1[i]    = inRow->Weight1[i];
	outRow->refant1[i]    = inRow->RefAnt1[i];
      }
      if (numPol>1) {   /* 2 poln */
	for (i=0; i<numIF; i++) {
	  outRow->phase1[i]      = atan2 (inRow->Imag2[i], inRow->Real2[i]);
	  outRow->real2[i]   = inRow->Real2[i];
	  outRow->imag2[i]   = inRow->Imag2[i];
	  outRow->rate2[i]   = inRow->Rate2[i];
	  outRow->delay2[i]  = inRow->Delay2[i];
	  outRow->weight2[i] = inRow->Weight2[i];
	  outRow->refant2[i] = inRow->RefAnt2[i];
	}
      }
      /* Write */
      oRow = -1;
      if ((ObitTableIDI_CALIBRATIONWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating IDI_CALIBRATION Table");
	return;
      }

    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_CALIBRATIONClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_CALIBRATION Table file");
      return;
    }
    
    if ((ObitTableCLClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output CL Table file");
      return;
    }

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied CALIBRATION table %d", iver);

    /* Cleanup */
    inRow    = ObitTableIDI_CALIBRATIONRowUnref(inRow);
    inTable  = ObitTableIDI_CALIBRATIONUnref(inTable);
    outRow   = ObitTableCLRowUnref(outRow);
    outTable = ObitTableCLUnref(outTable);
  } /* end loop over versions */

} /* end  PutCalibrationInfo */

void PutBandpassInfo (ObitInfoList *myInput, ObitUV *inData, 
		      ObitData *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any BANDPASS tables on inData to AIPS BP on outData           */
/*  Nothing is copied if bandpass is being applied                        */
/*   Input:                                                               */
/*      myInput  parser object                                            */
/*      inData   Input UV object                                          */
/*      outData  Output  IDI FITS object                                  */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableBP*          inTable=NULL;
  ObitTableBPRow*       inRow=NULL;
  ObitTableIDI_BANDPASS*    outTable=NULL;
  ObitTableIDI_BANDPASSRow* outRow=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, iver, hiVer, lim, iRow, oRow, ver, doBand;
  oint numIF, numAnt, numPol, numBach, strtChn;
  ObitIOAccess access;
  gchar *routine = "PutBandpassInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(outData));
  g_assert (ObitUVIsA(inData));

  /* Is bandpass being applied? */
  doBand = -1;
  ObitInfoListGetTest(myInput, "doBand", &type, dim, &doBand);
  if (doBand>0) return;

  /* Loop over plausible versions */
  hiVer = ObitTableListGetHigh (inData->tableList, "AIPS BP");
  for (iver=1; iver<=hiVer; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input BP table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    numIF = numPol = numBach = 0;
    inTable = newObitTableBPValue ("Input table", (ObitData*)inData, 
				   &ver, access, 
				   numPol, numIF, numBach, 
				   err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      continue;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableBPOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input BP table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableBPRow (inTable);
    
    /* Create output BANDPASS table object */
    ver = iver;
    access  = OBIT_IO_WriteOnly;
    numIF   = inTable->numIF;
    numAnt  = maxAnt;
    numPol  = inTable->numPol;
    numBach = inTable->numChan;
    strtChn = 1;
    outTable = newObitTableIDI_BANDPASSValue ("Output table", (ObitData*)outData, 
					      &ver, access, numIF, numAnt, numPol, 
					      numBach, strtChn, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_BANDPASS table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableIDI_BANDPASSOpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output IDI_BANDPASS table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTableIDI_BANDPASSRow (outTable);
    /* attach to table buffer */
    ObitTableIDI_BANDPASSSetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Set BANDPASS table values */
    outTable->no_stkd  = inData->myDesc->inaxes[inData->myDesc->jlocs];
    outTable->stk_1    = (olong)inData->myDesc->crval[inData->myDesc->jlocs];
    outTable->no_band  = inData->myDesc->inaxes[inData->myDesc->jlocif];
    outTable->no_chan  = inData->myDesc->inaxes[inData->myDesc->jlocf];
    outTable->ref_freq = inData->myDesc->crval[inData->myDesc->jlocf];
    outTable->chan_bw  = inData->myDesc->cdelt[inData->myDesc->jlocf];
    outTable->ref_pixl = inData->myDesc->crpix[inData->myDesc->jlocf];
    lim = MIN (UVLEN_VALUE,MAXKEYCHARTABLEIDI_ARRAY_GEOMETRY);
    strncpy (outTable->obscode,   inData->myDesc->observer, 8);
    lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, UVLEN_VALUE);
    strncpy (outTable->RefDate, inData->myDesc->obsdat, lim);
    outTable->myStatus = OBIT_Modified; /* Mark as modified */
    
    /* Initialize output row */
    outRow->status      = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableBPReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading BP Table");
	return;
      }
      
      /* Save to BP table */
      outRow->Time        = inRow->Time;
      outRow->TimeI       = inRow->TimeI;
      outRow->SourID      = inRow->SourID;
      outRow->Array       = inRow->SubA;
      outRow->antNo       = inRow->antNo;
      outRow->fqid        = inRow->FreqID;
      for (i=0; i<numIF; i++) {
	outRow->refant1[i]   = inRow->RefAnt1;
	outRow->breal1[i]    = inRow->Real1[i];
	outRow->bimag1[i]    = inRow->Imag1[i];
      }
      if (numPol>1) {   /* 2 poln */
	for (i=0; i<numIF; i++) {
	  outRow->refant2[i]   = inRow->RefAnt2;
	  outRow->breal2[i]    = inRow->Real2[i];
	  outRow->bimag2[i]    = inRow->Imag2[i];
	}
      }
      /* Write */
      oRow = -1;
      if ((ObitTableIDI_BANDPASSWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating IDI_BANDPASS Table");
	return;
      }

    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_BANDPASSClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_BANDPASS Table file");
      return;
    }
    
    if ((ObitTableBPClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output BP Table file");
      return;
    }

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied BANDPASS table %d", iver);

    /* Cleanup */
    inRow   = ObitTableBPRowUnref(inRow);
    inTable = ObitTableBPUnref(inTable);
    outRow    = ObitTableIDI_BANDPASSRowUnref(outRow);
    outTable  = ObitTableIDI_BANDPASSUnref(outTable);
  } /* end loop over versions */

} /* end  PutBandpassInfo */

void PutTSysInfo (ObitInfoList *myInput, ObitUV *inData, ObitData *outData, 
		  ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any SYSTEM_TEMPERATURE tables on inData to AIPS TY on outData */
/*   Input:                                                               */
/*      myInput  parser object                                            */
/*      inData   Input UV object                                          */
/*      outData  Output IDI FITS object                                   */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableTY*          inTable=NULL;
  ObitTableTYRow*       inRow=NULL;
  ObitTableIDI_SYSTEM_TEMPERATURE*    outTable=NULL;
  ObitTableIDI_SYSTEM_TEMPERATURERow* outRow=NULL;
  olong i, iver, hiVer, lim, iRow, oRow, ver;
  oint numIF, numPol;
  ObitIOAccess access;
  gchar *routine = "PutTSysInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(outData));
  g_assert (ObitUVIsA(inData));

  /* Loop over plausible versions */
  hiVer = ObitTableListGetHigh (inData->tableList, "AIPS TY");
  for (iver=1; iver<=hiVer; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input TY table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    numIF = numPol = 0;
    inTable = newObitTableTYValue ("Input table", (ObitData*)inData, 
				   &ver, access, numPol, numIF, err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      continue;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableTYOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input TY table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableTYRow (inTable);
    
    /* Create output SYSTEM_TEMPERATURE table object */
    ver = iver;
    access  = OBIT_IO_WriteOnly;
    numIF   = inTable->numIF;
    numPol  = inTable->numPol;
    outTable = newObitTableIDI_SYSTEM_TEMPERATUREValue ("Output table", outData, 
				    &ver, access, numIF, numPol, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_SYSTEM_TEMPERATURE  table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableIDI_SYSTEM_TEMPERATUREOpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output IDI_SYSTEM_TEMPERATURE table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTableIDI_SYSTEM_TEMPERATURERow (outTable);
    /* attach to table buffer */
    ObitTableIDI_SYSTEM_TEMPERATURESetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Set SYSTEM_TEMPERATURE table values */
    outTable->no_stkd  = inData->myDesc->inaxes[inData->myDesc->jlocs];
    outTable->stk_1    = (olong)inData->myDesc->crval[inData->myDesc->jlocs];
    outTable->no_band  = inData->myDesc->inaxes[inData->myDesc->jlocif];
    outTable->no_chan  = inData->myDesc->inaxes[inData->myDesc->jlocf];
    outTable->ref_freq = inData->myDesc->crval[inData->myDesc->jlocf];
    outTable->chan_bw  = inData->myDesc->cdelt[inData->myDesc->jlocf];
    outTable->ref_pixl = inData->myDesc->crpix[inData->myDesc->jlocf];
    lim = MIN (UVLEN_VALUE,MAXKEYCHARTABLEIDI_ARRAY_GEOMETRY);
    strncpy (outTable->obscode,   inData->myDesc->observer, 8);
    lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, UVLEN_VALUE);
    strncpy (outTable->RefDate, inData->myDesc->obsdat, lim);
    outTable->myStatus = OBIT_Modified; /* Mark as modified */
    
    /* Initialize output row */
    outRow->status      = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableTYReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading TY Table");
	return;
      }
      
      /* Save to TY table */
      outRow->Time      = inRow->Time;
      outRow->TimeI     = inRow->TimeI;
      outRow->SourID    = inRow->SourID;
      outRow->antNo     = inRow->antennaNo;
      outRow->Array     = inRow->SubA ;
      outRow->fqid      = inRow->FreqID;
      for (i=0; i<numIF; i++) {
	outRow->TSys1[i] = inRow->Tsys1[i];
	outRow->TAnt1[i] = inRow->Tant1[i];
      }
      if (numPol>1) {   /* 2 poln */
	for (i=0; i<numIF; i++) {
	  outRow->TSys2[i] = inRow->Tsys2[i];
	  outRow->TAnt2[i] = inRow->Tant2[i];
	}
      }
      /* Write */
      oRow = -1;
      if ((ObitTableIDI_SYSTEM_TEMPERATUREWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating IDI_SYSTEM_TEMPERATURE Table");
	return;
      }

    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_SYSTEM_TEMPERATUREClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_SYSTEM_TEMPERATURE Table file");
      return;
    }
    
    if ((ObitTableTYClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output TY Table file");
      return;
    }

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied SYSTEM_TEMPERATURE table %d", iver);

    /* Cleanup */
    inRow   = ObitTableTYRowUnref(inRow);
    inTable = ObitTableTYUnref(inTable);
    outRow    = ObitTableIDI_SYSTEM_TEMPERATURERowUnref(outRow);
    outTable  = ObitTableIDI_SYSTEM_TEMPERATUREUnref(outTable);
  } /* end loop over versions */

} /* end  PutTSysInfo */

void PutWeatherInfo (ObitInfoList *myInput, ObitUV *inData, 
		     ObitData *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any WEATHER tables on inData to AIPS WX on outData            */
/*   Input:                                                               */
/*      myInput  Input parser object                                      */
/*      inData   Input UV object                                          */
/*      outData  Output IDI FITS object                                   */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableWX*          inTable=NULL;
  ObitTableWXRow*       inRow=NULL;
  ObitTableIDI_WEATHER*    outTable=NULL;
  ObitTableIDI_WEATHERRow* outRow=NULL;
  olong iver, hiVer, lim, iRow, oRow, ver;
  oint numIF;
  ObitIOAccess access;
  gchar *routine = "PutWeatherInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(outData));
  g_assert (ObitUVIsA(inData));

  /* Loop over plausible versions */
  hiVer = ObitTableListGetHigh (inData->tableList, "AIPS WX");
  for (iver=1; iver<=hiVer; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input WX table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    numIF =0;
    inTable = newObitTableWXValue ("Input table", (ObitData*)inData, 
				    &ver, access, err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      continue;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableWXOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input WX table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableWXRow (inTable);
    
    /* Create output IDI_WEATHER table object */
    ver = iver;
    access  = OBIT_IO_WriteOnly;
    outTable = newObitTableIDI_WEATHERValue ("Output table", (ObitData*)outData, 
				    &ver, access, numIF, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_WEATHER table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableIDI_WEATHEROpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output IDI_WEATHER table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTableIDI_WEATHERRow (outTable);
    /* attach to table buffer */
    ObitTableIDI_WEATHERSetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Set WEATHER table values */
    outTable->no_stkd  = inData->myDesc->inaxes[inData->myDesc->jlocs];
    outTable->stk_1    = (olong)inData->myDesc->crval[inData->myDesc->jlocs];
    outTable->no_band  = inData->myDesc->inaxes[inData->myDesc->jlocif];
    outTable->no_chan  = inData->myDesc->inaxes[inData->myDesc->jlocf];
    outTable->ref_freq = inData->myDesc->crval[inData->myDesc->jlocf];
    outTable->chan_bw  = inData->myDesc->cdelt[inData->myDesc->jlocf];
    outTable->ref_pixl = inData->myDesc->crpix[inData->myDesc->jlocf];
    lim = MIN (UVLEN_VALUE,MAXKEYCHARTABLEIDI_ARRAY_GEOMETRY);
    strncpy (outTable->obscode,   inData->myDesc->observer, 8);
    lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, UVLEN_VALUE);
    strncpy (outTable->RefDate, inData->myDesc->obsdat, lim);
    outTable->myStatus = OBIT_Modified; /* Mark as modified */
    
    /* Initialize output row */
    outRow->status      = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableWXReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading WX Table");
	return;
      }
      
      /* Save to WX table */
      outRow->Time           = inRow->Time;
      outRow->TimeI          = inRow->TimeI;
      outRow->antNo          = inRow->antNo;
      outRow->temperature    = inRow->temperature;
      outRow->pressure       = inRow->pressure;
      outRow->dewpoint       = inRow->dewpoint;
      outRow->wind_velocity  = inRow->windVelocity;
      outRow->wind_direction = inRow->windDirection;
      outRow->wvr_h2o        = inRow->wvrH2O;
      outRow->ionos_electron = inRow->onosElectron;

      /* Write */
      oRow = -1;
      if ((ObitTableIDI_WEATHERWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating IDI_WEATHER Table");
	return;
      }

    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_WEATHERClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_WEATHER Table file");
      return;
    }
    
    if ((ObitTableWXClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output WX Table file");
      return;
    }

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied WEATHER table %d", iver);

    /* Cleanup */
    inRow   = ObitTableWXRowUnref(inRow);
    inTable = ObitTableWXUnref(inTable);
    outRow    = ObitTableIDI_WEATHERRowUnref(outRow);
    outTable  = ObitTableIDI_WEATHERUnref(outTable);
  } /* end loop over versions */

} /* end  PutWeatherInfo */

void PutData (ObitUV *inData, ObitData *outData, 
	      ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Write visibility data to FITS IDI                                     */
/*   Input:                                                               */
/*      inData   Input UV object                                          */
/*      outData  Output FITS IDI file                                     */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitUVDesc    *desc=NULL;
  ObitTableIDI_UV_DATA    *outTable=NULL;
  ObitTableIDI_UV_DATARow *outRow=NULL;
  ObitIOCode retCode;
  olong lim, oRow, i, nwt;
  oint no_band=0;
  ObitIOAccess access;
  ofloat uvwFact;
  olong ver;
  gchar *routine = "PutData";
  int cnt=0;  /* DEBUG*/

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(outData));
  g_assert(myInput!=NULL);

  /* Define output IDI UV_DATA */
  desc = inData->myDesc;
  no_band = desc->inaxes[desc->jlocif];
  
  /* Create input UV_DATA table object */
  ver = 1;
  access = OBIT_IO_WriteOnly;
  outTable = myObitTableIDI_UV_DATAValue ("Input table", outData, 
					  &ver, inData, access, no_band, 
					  err);
  if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_UV_DATA table");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Open output table */
  if ((ObitTableIDI_UV_DATAOpen (outTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_UV_DATA table");
    return;
  }
  
  /* Set UV_DATA table values */
  outTable->no_stkd  = desc->inaxes[desc->jlocs];
  outTable->stk_1    = (olong)desc->crval[desc->jlocs];
  outTable->no_band  = desc->inaxes[desc->jlocif];
  outTable->no_chan  = desc->inaxes[desc->jlocf];
  outTable->ref_freq = desc->crval[desc->jlocf];
  outTable->chan_bw  = desc->cdelt[desc->jlocf];
  outTable->ref_pixl = desc->crpix[desc->jlocf];
  lim = MIN (UVLEN_VALUE,MAXKEYCHARTABLEIDI_ARRAY_GEOMETRY);
  strncpy (outTable->obscode,   desc->observer, 8);
  lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, UVLEN_VALUE);    
  strncpy (outTable->RefDate, desc->obsdat, lim);
  
  /* Create output row */
  outRow = newObitTableIDI_UV_DATARow (outTable);
  /* attach to table buffer */
  ObitTableIDI_UV_DATASetRow (outTable, outRow, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Open input data */
  retCode = ObitUVOpen (inData, OBIT_IO_ReadWrite, err);
  if (( retCode != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input file %s", inData->name);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Initialize output row */
  outRow->uu        = 0.0;
  outRow->vv        = 0.0;
  outRow->ww        = 0.0;
  outRow->date      = desc->JDObs;
  outRow->Time      = 0.0;
  outRow->Baseline  = 0;
  outRow->Array     = 0;
  outRow->Baseline  = 0;
  outRow->Source    = 0;
  outRow->FreqID    = 0;
  outRow->IntTim    = 0.0;
  outRow->GateID    = 0;
  nwt = outTable->no_stkd * outTable->no_band;
  for (i=0; i<nwt; i++) outRow->Weight[i] = 1.0;
  outRow->Flux[0]   = 0.0;
  outRow->status    = 0;

  /* Factor to convert from wavelengths at the reference frequency to sec */
  uvwFact = 1.0 / desc->crval[desc->jlocf];

  /* Loop over file - reading one visibility per call */\

  while (retCode==OBIT_IO_OK) {
    /* read 1 vis */
    retCode = ObitUVRead (inData, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);

    /* Convert */
    if (desc->ilocu>=0) outRow->uu = uvwFact*inData->buffer[desc->ilocu];
    if (desc->ilocv>=0) outRow->vv = uvwFact*inData->buffer[desc->ilocv];
    if (desc->ilocw>=0) outRow->ww = uvwFact*inData->buffer[desc->ilocw];
    if (desc->iloct>=0) outRow->Time = (odouble)inData->buffer[desc->iloct];
    if (desc->ilocb>=0) {
      outRow->Baseline = (olong)inData->buffer[desc->ilocb];
      outRow->Array = (olong)(1 + (inData->buffer[desc->ilocb]-outRow->Baseline)*100);
    }
    if (desc->ilocsu>=0) outRow->Source = (olong)inData->buffer[desc->ilocsu];
    if (desc->ilocfq>=0) outRow->FreqID = (olong)inData->buffer[desc->ilocfq];
    if (desc->ilocit>=0) outRow->IntTim = inData->buffer[desc->ilocit];
    /* Copy correlation data */
    for (i=0; i<desc->ncorr*3; i++) outRow->Flux[i] = inData->buffer[desc->nrparm+i];

    /* Write FITS IDI */
    oRow = -1;
    retCode = myObitTableIDI_UV_DATAWriteRow (outTable, oRow, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);    

    if (inData->myDesc->firstVis>=inData->myDesc->nvis) break; /* done? */
    cnt++; /* DEBUG */
     /*if (cnt>2000) break;   DEBUG */

  } /* end loop over input */

  /* Close input uv data */
  if ((ObitUVClose (inData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing input file");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Close output table */
  if ((myObitTableIDI_UV_DATAClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_UV_DATA Table file");
    return;
  }

  /* Cleanup */
  outRow = ObitTableIDI_UV_DATARowUnref(outRow);

  return;
} /* end PutData */

/*----------------------------------------------------------------------- */
/*  Write History for IDIOut                                              */
/*   Input:                                                               */
/*      inData   ObitUV to copy history from                              */
/*      myInput   Input parameters on InfoList                            */
/*      outData   FITS IDI to write history to                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void IDIOutHistory (ObitUV* inData, ObitInfoList* myInput, ObitData* outData, 
		    ObitErr* err)
{
  ObitHistory  *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "FreqID", "BChan", "EChan", "BIF", "EIF",  "Stokes", 
    "Sources",  "Qual", "souCode", "subA", "Antennas", 
    "doCalSelect", "doCalib", "gainUse", "doPol", "flagVer", 
    "doBand", "BPVer", "Smooth",  
    NULL};
  gchar *routine = "IDIOutHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));
  g_assert (ObitDataIsA(outData));

  /* Do history  */
  inHistory  = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory ((ObitData*)outData, OBIT_IO_WriteOnly, err);

  /* Copy to history table */
  ObitHistoryCopy (inHistory, outHistory, err);

  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end IDIOutHistory  */

/**
 * Constructor from values.
 * Creates a new table structure and attaches to the TableList of file.
 * If the specified table already exists then it is returned.
 * Initializes class if needed on first call.
 * Forces an update of any disk resident structures (e.g. AIPS header).
 * \param name   An optional name for the object.
 * \param file   ObitData which which the table is to be associated.
 * \param ver    Table version number. 0=> add higher, value used returned
 * \param inData ObitUV to be copied
 * \param access access (OBIT_IO_ReadOnly, means do not create if it doesn't exist.
 * \param no_band Number of frequency bands (IF)
 * \param err Error stack, returns if not empty.
 * \return the new object, NULL on failure.
 */
ObitTableIDI_UV_DATA* myObitTableIDI_UV_DATAValue (gchar* name, ObitData *file, olong *ver,
						   ObitUV *inData,
						   ObitIOAccess access,
						   oint no_band, ObitErr *err)
{
  ObitTableIDI_UV_DATA* out=NULL;
  ObitTable     *testTab=NULL;
  ObitTableDesc *desc=NULL;
  ObitTableList *list=NULL;
  ObitInfoList  *info=NULL;
  ObitUVDesc    *uvDesc=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean optional;
  olong colNo, i, ncol;
  ObitIOCode retCode;
  gchar *tabType = "UV_DATA";
  gchar *routine = "myObitTableIDI_UV_DATAValue";

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitDataIsA(file));

  /* Class initialization if needed */
  ObitTableIDI_UV_DATAClassInit();

  /* Get TableList */
  list = ((ObitData*)file)->tableList;
  info = ((ObitData*)file)->info;

  /* Create output table */
  /* create basal table */
  testTab = newObitDataTable ((ObitData*)file, access, tabType,
			       ver, err);
  if (err->error) Obit_traceback_val (err, routine,"", out);
  
  /* likely need to convert */
  if (ObitTableIDI_UV_DATAIsA(testTab)) { 
    out = ObitTableRef(testTab);
  } else { /* needs conversion */
    out = ObitTableIDI_UV_DATAConvert(testTab);
  }
  testTab = ObitTableUnref(testTab); /* remove reference */
  /*  out = newObitTableIDI_UV_DATA("UV_DATA");*/

  /* Update the TableList */
  *ver = 1;
  ObitTableListPut(list, tabType, ver, (ObitTable*)out, err);
  if (err->error) Obit_traceback_val (err, routine,"", out);

  /* Set values */
  uvDesc = inData->myDesc;
  out->no_band = MAX (0, no_band);
  out->tabrev = 1;
  out->no_stkd = 1;
  out->stk_1 = -1;
  out->no_chan = 1;
  out->ref_freq = 1.0;
  out->chan_bw = 1.0;
  out->ref_pixl = 1.0;
  strncpy (out->obscode, uvDesc->observer, MAXKEYCHARTABLEIDI_UV_DATA );
  strncpy (out->observer, uvDesc->observer, MAXKEYCHARTABLEIDI_UV_DATA );
  strncpy (out->teles, uvDesc->teles, MAXKEYCHARTABLEIDI_UV_DATA );
  strncpy (out->RefDate, uvDesc->obsdat, MAXKEYCHARTABLEIDI_UV_DATA );
  out->nmatrix = 1;
  out->maxis = uvDesc->naxis;
  strncpy (out->Equinox, "J2000   ", MAXKEYCHARTABLEIDI_UV_DATA );
  strncpy (out->WeighTyp, "NORMAL", MAXKEYCHARTABLEIDI_UV_DATA );
  strncpy (out->dateObs, uvDesc->obsdat, MAXKEYCHARTABLEIDI_UV_DATA );
  out->visScale = 1.0;

  /* initialize descriptor */
  desc   = out->myDesc;

  /* How many columns actually in table? */
  ncol = 14 + out->no_band*0;
  desc->FieldName = g_malloc0((ncol+1)*sizeof(gchar*));
  desc->FieldUnit = g_malloc0((ncol+1)*sizeof(gchar*));
  desc->type      = g_malloc0((ncol+1)*sizeof(ObitInfoType));
  desc->dim       = g_malloc0((ncol+1)*sizeof(gint32*));
  for (i=0; i<ncol+1; i++) 
    desc->dim[i] = g_malloc0(7*sizeof(gint32));

  desc->TableName = g_strdup(tabType);
  desc->sort[0] = uvDesc->isort[0];
  desc->sort[1] = uvDesc->isort[1];
  colNo = 0;

  /* Define Columns */
  desc->FieldName[colNo] = g_strdup("UU      ");
  desc->FieldUnit[colNo] = g_strdup("SECONDS");
  desc->type[colNo] = OBIT_float;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("VV      ");
  desc->FieldUnit[colNo] = g_strdup("SECONDS");
  desc->type[colNo] = OBIT_float;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("WW      ");
  desc->FieldUnit[colNo] = g_strdup("SECONDS");
  desc->type[colNo] = OBIT_float;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("DATE    ");
  desc->FieldUnit[colNo] = g_strdup("DAYS");
  desc->type[colNo] = OBIT_double;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("TIME    ");
  desc->FieldUnit[colNo] = g_strdup("DAYS");
  desc->type[colNo] = OBIT_double;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("BASELINE");
  desc->FieldUnit[colNo] = g_strdup("");
  desc->type[colNo] = OBIT_oint;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("ARRAY   ");
  desc->FieldUnit[colNo] = g_strdup("");
  desc->type[colNo] = OBIT_oint;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  if (uvDesc->ilocsu>=0) {
    desc->FieldName[colNo] = g_strdup("SOURCE  ");
    desc->FieldUnit[colNo] = g_strdup("");
    desc->type[colNo] = OBIT_oint;
    for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
    colNo++;
  }
  if (uvDesc->ilocfq>=0) {
    desc->FieldName[colNo] = g_strdup("FREQID  ");
    desc->FieldUnit[colNo] = g_strdup("");
    desc->type[colNo] = OBIT_oint;
    for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
    colNo++;
  }
  if (uvDesc->ilocit>=0) {
    desc->FieldName[colNo] = g_strdup("INTTIM");
    desc->FieldUnit[colNo] = g_strdup("SECONDS");
    desc->type[colNo] = OBIT_float;
    for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
    colNo++;
  }
  /* Weight */
  optional = FALSE;

  desc->FieldName[colNo] = g_strdup("WEIGHT");
  desc->FieldUnit[colNo] = g_strdup("");
  desc->type[colNo] = OBIT_float;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  /* Dimension noPoln x noBand */
  desc->dim[colNo][0] = uvDesc->inaxes[uvDesc->jlocs];
  desc->dim[colNo][1] = no_band;
  colNo++;

  optional = FALSE;
  desc->FieldName[colNo] = g_strdup("FLUX    ");
  desc->FieldUnit[colNo] = g_strdup("");
  desc->type[colNo] = OBIT_float;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;

  /* Regular axes */
  if (uvDesc->naxis>=1) {
    out->maxis1 = MAX (0, uvDesc->inaxes[0]);
    desc->dim[colNo][0] = out->maxis1;
    out->cdelt1 = uvDesc->cdelt[0];
    out->crpix1 = uvDesc->crpix[0];
    out->crval1 = uvDesc->crval[0];
    strncpy (out->ctype1, uvDesc->ctype[0], 8);
  }
  if (uvDesc->naxis>=2) {
    out->maxis2 = MAX (0, uvDesc->inaxes[1]);
    desc->dim[colNo][1] = out->maxis2;
    out->cdelt2 = uvDesc->cdelt[1];
    out->crpix2 = uvDesc->crpix[1];
    out->crval2 = uvDesc->crval[1];
    strncpy (out->ctype2, uvDesc->ctype[1], 8);
 }
  if (uvDesc->naxis>=3) {
    out->maxis3 = MAX (0, uvDesc->inaxes[2]);
    desc->dim[colNo][2] = out->maxis3;
    out->cdelt3 = uvDesc->cdelt[2];
    out->crpix3 = uvDesc->crpix[2];
    out->crval3 = uvDesc->crval[2];
    strncpy (out->ctype3, uvDesc->ctype[2], 8);
}
  if (uvDesc->naxis>=4) {
    out->maxis4 = MAX (0, uvDesc->inaxes[3]);
    desc->dim[colNo][3] = out->maxis4;
    out->cdelt4 = uvDesc->cdelt[3];
    out->crpix4 = uvDesc->crpix[3];
    out->crval4 = uvDesc->crval[3];
    strncpy (out->ctype4, uvDesc->ctype[3], 8);
  }
  if (uvDesc->naxis>=5) {
    out->maxis5 = MAX (0, uvDesc->inaxes[4]);
    desc->dim[colNo][4] = out->maxis5;
    out->cdelt5 = uvDesc->cdelt[4];
    out->crpix5 = uvDesc->crpix[4];
    out->crval5 = uvDesc->crval[4];
    strncpy (out->ctype5, uvDesc->ctype[4], 8);
 }
  if (uvDesc->naxis>=6) {
    out->maxis6 = MAX (0, uvDesc->inaxes[5]);
    desc->dim[colNo][5] = out->maxis6;
    out->cdelt6 = uvDesc->cdelt[5];
    out->crpix6 = uvDesc->crpix[5];
    out->crval6 = uvDesc->crval[5];
    strncpy (out->ctype6, uvDesc->ctype[5], 8);
 }
  if (uvDesc->naxis>=7) {
    out->maxis7 = MAX (0, uvDesc->inaxes[6]);
    desc->dim[colNo][6] = out->maxis7;
    out->cdelt7 = uvDesc->cdelt[6];
    out->crpix7 = uvDesc->crpix[6];
    out->crval7 = uvDesc->crval[6];
    strncpy (out->ctype7, uvDesc->ctype[6], 8);
 }
  colNo++;
  /* Mark data matrix */
  out->tmatx01 = colNo==1;
  out->tmatx02 = colNo==2;
  out->tmatx03 = colNo==3;
  out->tmatx04 = colNo==4;
  out->tmatx05 = colNo==5;
  out->tmatx06 = colNo==6;
  out->tmatx07 = colNo==7;
  out->tmatx08 = colNo==8;
  out->tmatx09 = colNo==9;
  out->tmatx10 = colNo==10;
  out->tmatx11 = colNo==11;
  out->tmatx12 = colNo==12;
  out->tmatx13 = colNo==13;
  out->tmatx14 = colNo==14;
  out->tmatx15 = colNo==15;

  /* Add _status column at end */
  desc->FieldName[colNo] = g_strdup("_status");
  desc->FieldUnit[colNo] = g_strdup("        ");
  desc->type[colNo] = OBIT_short;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  
  /* number of fields */
  desc->nfield = colNo + 1;

  /* Sort order */
  dim[0] = 2; dim[1] = dim[3] = 1;
  uvDesc->isort[2] = 0;
  ObitInfoListAlwaysPut (desc->info, "OBITSORT", OBIT_string, dim, uvDesc->isort);

  /* initialize descriptor keywords */
  myObitTableIDI_UV_DATADumpKey (out, err);
 
  /* index table descriptor */
  ObitTableDescIndex (desc);

  /* Open and Close to fully instantiate */
  retCode = ObitTableIDI_UV_DATAOpen(out, OBIT_IO_WriteOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, out->name, out);    
  
  retCode = ObitTableIDI_UV_DATAClose(out, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, out->name, out); 

  /* Force update of disk resident info */
  retCode = ObitIOUpdateTables (((ObitData*)file)->myIO, info, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, out->name, out); 
  
 return out;
} /* end myObitTableIDI_UV_DATAValue */

/**
 * Copy table specific (keyword) information  to infolist.
 * \param info Table to update
 * \param err  ObitErr for reporting errors.
 */
void myObitTableIDI_UV_DATADumpKey (ObitTableIDI_UV_DATA *in, ObitErr *err)
{
  ObitInfoList *info=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

 /* error checks */
   g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* Set Keywords */
  if (in->myIO!=NULL) info = ((ObitTableDesc*)(in->myIO->myDesc))->info;
  else info = in->myDesc->info;
  /* tabrev */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "TABREV", type, dim, 
		  (gpointer)&in->tabrev);
  /* no_stkd */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "NO_STKD", type, dim, 
		  (gpointer)&in->no_stkd);
  /* stk_1 */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "STK_1", type, dim, 
		  (gpointer)&in->stk_1);
  /* no_band */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "NO_BAND", type, dim, 
		  (gpointer)&in->no_band);
  /* no_chan */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "NO_CHAN", type, dim, 
		  (gpointer)&in->no_chan);
  /* ref_freq */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "REF_FREQ", type, dim, 
		  (gpointer)&in->ref_freq);
  /* chan_bw */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CHAN_BW", type, dim, 
		  (gpointer)&in->chan_bw);
  /* ref_pixl */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "REF_PIXL", type, dim, 
		  (gpointer)&in->ref_pixl);
  /* obscode */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "OBSCODE", type, dim, 
		  (gpointer)&in->obscode);
  /* ArrName */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "ARRNAM", type, dim, 
		  (gpointer)&in->ArrName);
  /* RefDate */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "RDATE", type, dim, 
		  (gpointer)&in->RefDate);
  /* Equinox */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "EQUINOX", type, dim, 
		  (gpointer)&in->Equinox);
  /* WeighTyp */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "WEIGHTYP", type, dim, 
		  (gpointer)&in->WeighTyp);
  /* nmatrix */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "NMATRIX", type, dim, 
		  (gpointer)&in->nmatrix);
  /* maxis */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "MAXIS", type, dim, 
		  (gpointer)&in->maxis);
  /* tmatx01 */
  if (in->tmatx01) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX01", type, dim, 
			  (gpointer)&in->tmatx01);
  }
  /* tmatx02 */
  if (in->tmatx02) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX02", type, dim, 
			  (gpointer)&in->tmatx02);
  }
  /* tmatx03 */
  if (in->tmatx03) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX03", type, dim, 
			  (gpointer)&in->tmatx03);
  }
  /* tmatx04 */
  if (in->tmatx04) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX04", type, dim, 
			  (gpointer)&in->tmatx04);
  }
  /* tmatx05 */
  if (in->tmatx05) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX05", type, dim, 
			  (gpointer)&in->tmatx05);
  }
  /* tmatx06 */
  if (in->tmatx06) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX06", type, dim, 
			  (gpointer)&in->tmatx06);
  }
  /* tmatx07 */
  if (in->tmatx07) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX07", type, dim, 
			  (gpointer)&in->tmatx07);
  }
  /* tmatx08 */
  if (in->tmatx08) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX08", type, dim, 
			  (gpointer)&in->tmatx08);
  }
  /* tmatx09 */
  if (in->tmatx09) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX09", type, dim, 
			  (gpointer)&in->tmatx09);
  }
  /* tmatx10 */
  if (in->tmatx10) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX10", type, dim, 
			  (gpointer)&in->tmatx10);
  }
  /* tmatx11 */
  if (in->tmatx11) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX11", type, dim, 
			  (gpointer)&in->tmatx11);
  }
  /* tmatx12 */
  if (in->tmatx12) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX12", type, dim, 
			  (gpointer)&in->tmatx12);
  }
  /* tmatx13 */
  if (in->tmatx13) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX13", type, dim, 
			  (gpointer)&in->tmatx13);
  }
  /* tmatx14 */
  if (in->tmatx14) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX14", type, dim, 
			  (gpointer)&in->tmatx14);
  }
  /* tmatx15 */
  if (in->tmatx15) {
    type  = OBIT_bool;
    dim[0] = 1;
    ObitInfoListAlwaysPut(info, "TMATX15", type, dim, 
			  (gpointer)&in->tmatx15);
  }
  /* maxis1 */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "MAXIS1", type, dim, 
			(gpointer)&in->maxis1);
  /* ctype1 */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "CTYPE1", type, dim, 
			(gpointer)&in->ctype1);
  /* cdelt1 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CDELT1", type, dim, 
			(gpointer)&in->cdelt1);
  /* crpix1 */
  type  = OBIT_float;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRPIX1", type, dim, 
			(gpointer)&in->crpix1);
  /* crval1 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRVAL1", type, dim, 
			(gpointer)&in->crval1);
  
  /* maxis2 */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "MAXIS2", type, dim, 
			(gpointer)&in->maxis2);
  /* ctype2 */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "CTYPE2", type, dim, 
			(gpointer)&in->ctype2);
  /* cdelt2 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CDELT2", type, dim, 
			(gpointer)&in->cdelt2);
  /* crpix2 */
  type  = OBIT_float;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRPIX2", type, dim, 
			(gpointer)&in->crpix2);
  /* crval2 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRVAL2", type, dim, 
			(gpointer)&in->crval2);
  
  /* maxis3 */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "MAXIS3", type, dim, 
			  (gpointer)&in->maxis3);
  /* ctype3 */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "CTYPE3", type, dim, 
			(gpointer)&in->ctype3);
  /* cdelt3 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CDELT3", type, dim, 
			(gpointer)&in->cdelt3);
  /* crpix3 */
  type  = OBIT_float;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRPIX3", type, dim, 
			(gpointer)&in->crpix3);
  /* crval3 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRVAL3", type, dim, 
			(gpointer)&in->crval3);
  
  /* maxis4 */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "MAXIS4", type, dim, 
			(gpointer)&in->maxis4);
  /* ctype4 */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "CTYPE4", type, dim, 
			(gpointer)&in->ctype4);
  /* cdelt4 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CDELT4", type, dim, 
			(gpointer)&in->cdelt4);
  /* crpix4 */
  type  = OBIT_float;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRPIX4", type, dim, 
			(gpointer)&in->crpix4);
  /* crval4 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRVAL4", type, dim, 
			(gpointer)&in->crval4);
  
  /* maxis5 */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "MAXIS5", type, dim, 
			(gpointer)&in->maxis5);
  /* ctype5 */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "CTYPE5", type, dim, 
			(gpointer)&in->ctype5);
  /* cdelt5 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CDELT5", type, dim, 
			(gpointer)&in->cdelt5);
  /* crpix5 */
  type  = OBIT_float;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRPIX5", type, dim, 
			(gpointer)&in->crpix5);
  /* crval5 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRVAL5", type, dim, 
			(gpointer)&in->crval5);
  /* maxis6 */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "MAXIS6", type, dim, 
			(gpointer)&in->maxis6);
  /* ctype6 */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "CTYPE6", type, dim, 
			(gpointer)&in->ctype6);
  /* cdelt6 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CDELT6", type, dim, 
			(gpointer)&in->cdelt6);
  /* crpix6 */
  type  = OBIT_float;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRPIX6", type, dim, 
			(gpointer)&in->crpix6);
  /* crval6 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRVAL6", type, dim, 
			(gpointer)&in->crval6);
  /* maxis7 */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "MAXIS7", type, dim, 
			(gpointer)&in->maxis7);
  /* ctype7 */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "CTYPE7", type, dim, 
			(gpointer)&in->ctype7);
  /* cdelt7 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CDELT7", type, dim, 
			(gpointer)&in->cdelt7);
  /* crpix7 */
  type  = OBIT_float;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRPIX7", type, dim, 
			(gpointer)&in->crpix7);
  /* crval7 */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "CRVAL7", type, dim, 
			(gpointer)&in->crval7);

  /* dateObs */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "DATE-OBS", type, dim, 
		  (gpointer)&in->dateObs);
  /* teles */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "TELESCOP", type, dim, 
		  (gpointer)&in->teles);
  /* observer */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "OBSERVER", type, dim, 
		  (gpointer)&in->observer);
  /* visScale */
  type  = OBIT_float;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "VIS_SCAL", type, dim, 
		  (gpointer)&in->visScale);
  /* sort */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEIDI_UV_DATA;
  ObitInfoListAlwaysPut(info, "SORT", type, dim, 
		  (gpointer)&in->sort);
   
} /* end myObitTableIDI_UV_DATADumpKey */

/**
 * Write a table row.
 * Before calling this routine, the row structure needs to be initialized
 * and filled with data. The array members of the row structure are  
 * pointers to independently allocated memory.  These pointers can be set to the 
 * correct table buffer locations using ObitTableIDI_UV_DATASetRow  
 * \param in       Table to read
 * \param iIDI_UV_DATARow   Row number, -1 -> next
 * \param row Table Row structure containing data
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode 
myObitTableIDI_UV_DATAWriteRow  (ObitTableIDI_UV_DATA *in, olong iIDI_UV_DATARow, 
				 ObitTableIDI_UV_DATARow *row,
				 ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  /*gshort    *siRow;*/
  odouble   *dRow;
  oint      *iRow, i;
  ofloat    *fRow;
  /*gchar     *cRow;
    gboolean  *lRow;
    guint8    *bRow;*/
  gchar *routine = "myObitTableIDI_UV_DATAWriteRow";
  

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;

  /* Typed pointers to row of data */  
  dRow  = (odouble*)in->buffer;
  /*siRow = (gshort*)in->buffer;*/
  iRow  = (oint*)in->buffer;
  fRow  = (ofloat*)in->buffer;
  /*cRow  = (gchar*)in->buffer;*/
  /*lRow  = (gboolean*)in->buffer;*/
  /*bRow  = (guint8*)in->buffer;*/
  
  /* Make full copy of all data */
  fRow[in->uuOff]   = row->uu;
  fRow[in->vvOff]   = row->vv;
  fRow[in->wwOff]   = row->ww;
  dRow[in->dateOff] = row->date;
  dRow[in->TimeOff] = row->Time;
  iRow[in->BaselineOff] = row->Baseline;
  iRow[in->ArrayOff]    = row->Array;
  if (in->SourceOff>=0) iRow[in->SourceOff] = row->Source;
  if (in->FreqIDOff>=0) iRow[in->FreqIDOff] = row->FreqID;
  if (in->IntTimOff>=0) fRow[in->IntTimOff] = row->IntTim;
  if (in->WeightCol >= 0) { 
    for (i=0; i<in->myDesc->repeat[in->WeightCol]; i++) {
      fRow[in->WeightOff+i] = row->Weight[i];
    }
  } 
  if (in->FluxCol >= 0) { 
    for (i=0; i<in->myDesc->repeat[in->FluxCol]; i++) 
      fRow[in->FluxOff+i] = row->Flux[i];
  } 

  /* copy status */
  iRow[in->myDesc->statusOff] = row->status;
   
  /* Write one row */
  in->myDesc->numRowBuff = 1;
 
  /* Write row iIDI_UV_DATARow */
  retCode = ObitTableWrite ((ObitTable*)in, iIDI_UV_DATARow, NULL,  err);
  if (err->error) 
    Obit_traceback_val (err, routine,in->name, retCode);

  return retCode;
} /*  end myObitTableIDI_UV_DATAWriteRow */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
ObitIOCode myObitTableIDI_UV_DATAClose (ObitTableIDI_UV_DATA *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitTableIDI_UV_DATAClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;

  /* Something going on? */
  if (in->myStatus == OBIT_Inactive) return OBIT_IO_OK;

  /* Update keywords on descriptor if not ReadOnly*/
  if (in->myDesc->access != OBIT_IO_ReadOnly) 
    myObitTableIDI_UV_DATADumpKey (in, err);
  if (err->error) 
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Close */
  retCode = ObitTableClose ((ObitTable*)in, err);
  if (err->error) 
    Obit_traceback_val (err, routine, in->name, retCode);

  return retCode;
} /* end myObitTableIDI_UV_DATAClose */
