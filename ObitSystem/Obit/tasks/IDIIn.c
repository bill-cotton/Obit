/* $Id$  */
/* Read IDI format data, convert to Obit UV                           */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2007-2009                                          */
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
#include "ObitTableSU.h"
#include "ObitTableAN.h"
#include "ObitTableFQ.h"
#include "ObitTableFG.h"
#include "ObitTableCL.h"
#include "ObitTableBP.h"
#include "ObitTableTY.h"
#include "ObitTableIM.h"
#include "ObitTableGC.h"
#include "ObitTablePC.h"
#include "ObitTableWX.h"
#include "ObitTableIDI_ANTENNA.h"
#include "ObitTableIDI_ARRAY_GEOMETRY.h"
#include "ObitTableIDI_FREQUENCY.h"
#include "ObitTableIDI_SOURCE.h"
#include "ObitTableIDI_FLAG.h"
#include "ObitTableIDI_CALIBRATION.h"
#include "ObitTableIDI_BANDPASS.h"
#include "ObitTableIDI_SYSTEM_TEMPERATURE.h"
#include "ObitTableIDI_GAIN_CURVE.h"
#include "ObitTableIDI_PHASE_CAL.h"
#include "ObitTableIDI_INTERFEROMETER_MODEL.h"
#include "ObitTableIDI_WEATHER.h"
#include "ObitTableIDI_UV_DATA.h"
#include "ObitHistory.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* IDIInin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Create output uvdata */
ObitUV* setOutputData (ObitInfoList *myInput, ObitErr *err);
/* Get file descriptor */
void GetHeader (ObitUV *outData, gchar *infile, ObitInfoList *myInput, 
		ObitErr *err);
/* Get data */
void GetData (ObitUV *outData, gchar *infile, ObitInfoList *myInput, 
	      ObitErr *err);
/* Get Antenna info */
void GetAntennaInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Get Frequency info */
void GetFrequencyInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Get Source info */
void GetSourceInfo (ObitData *inData, ObitUV *outData, gboolean isNew, 
		    ObitErr *err);
/* Copy any FLAG tables */
void GetFlagInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Copy any CALIBRATION tables */
void GetCalibrationInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Copy any BANDPASS tables */
void GetBandpassInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Copy any SYSTEM_TEMPERATURE tables */
void GetTSysInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Copy any GAIN_CURVEtables */
void GetGainCurveInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Copy any PHASE_CALtables */
void GetPhaseCalInfo (ObitData *inData, ObitUV *outData, 
				 ObitErr *err);
/* Copy any INTTERFEROMETER_MODELtables */
void GetInterferometerModelInfo (ObitData *inData, ObitUV *outData, 
				 ObitErr *err);
/* Copy any WEATHER tables */
void GetWeatherInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Read data */
void ProcessData (gchar *inscan, ofloat avgTime,
		  olong *ndetect, olong *ntime, ofloat *refDate,
		  ofloat** ATime, ofloat*** AData, ofloat** ACal, 
		  ObitErr *err);
/* Write history */
void IDIInHistory (ObitData* inData, ObitInfoList* myInput, ObitUV* outData, 
		   ObitErr* err);

/* Program globals */
gchar *pgmName = "IDIIn";       /* Program name */
gchar *infile  = "IDIIn.inp";   /* File with program inputs */
gchar *outfile = "IDIIn.out";   /* File to contain program outputs */
olong  pgmNumber;      /* Program number (like POPS no.) */
olong  AIPSuser;       /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;        /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;        /* Number of FITS directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
gchar **FITSdirs=NULL; /* List of FITS data directories */
gchar DataRoot[128];   /* Root directory of input data */
odouble refMJD;        /* reference Julian date */
odouble integTime;     /* Integration time in days */
ofloat *SourceID = NULL; /* Source number (1-rel) lookup table */
gboolean isNew=FALSE;  /* Is the output newly created */
olong  nchan=1;        /* Number of frequencies */
olong  nstok=1;        /* Number of Stokes */
olong  nIF=1;          /* Number of IFs */
ofloat deltaFreq;      /* Channel increment */
ofloat refPixel;       /* Frequency reference pixel */
odouble refFrequency;  /* reference frequency (Hz) */


int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*    Read IDI  data to an Obit UV dataset                                */
/*----------------------------------------------------------------------- */
{
  olong  i, ierr=0;
  ObitSystem *mySystem= NULL;
  ObitUV *outData= NULL;
  ObitData *inData=NULL;
  ObitErr *err= NULL;
  gchar inscan[128];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar FullFile[128];
  olong disk;

  err = newObitErr();  /* Obit error/message stack */

  /* Startup - parse command line */
  ierr = 0;
  myInput = IDIInin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input File name */
  for (i=0; i<128; i++) inscan[i] = 0;
  ObitInfoListGet(myInput, "Scan", &type, dim, inscan, err);
  inscan[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(inscan);  /* Trim trailing blanks */

  /* Get input data file root */
  for (i=0; i<128; i++) DataRoot[i] = 0;
  ObitInfoListGet(myInput, "DataRoot", &type, dim, DataRoot, err);
  DataRoot[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(DataRoot);  /* Trim trailing blanks */

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);

  if (err->error) ierr = 1;  ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Create ObitUV for data */
  outData = setOutputData (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
   
  /* Get header info, array geometry, initialize output if necessary */
  GetHeader (outData, inscan, myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Tell what's going on */   
  Obit_log_error(err, OBIT_InfoErr, "Adding file %s to output", 
		 inscan);
  ObitErrLog(err);

  /* show any errors */
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
   
  /* Open output data */
  if ((ObitUVOpen (outData, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output UV file %s", outData->name);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* convert data  */
  GetData (outData, inscan, myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Close output uv data */
  if ((ObitUVClose (outData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing output file");
  
  /* Input file for tables */
  inData = newObitData("Input Data");
  disk = 0;
  /* Full input file name */
  sprintf (FullFile,"%s%s.fits", DataRoot, inscan);
  ObitDataSetFITS(inData, disk, FullFile, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Copy tables */
  GetFlagInfo (inData, outData, err);          /* FLAG tables */
  GetCalibrationInfo (inData, outData, err);   /* CALIBRATION tables */
  GetBandpassInfo (inData, outData, err);      /* BANDPASS tables */
  GetTSysInfo (inData, outData, err);          /* SYSTEM_TEMPERATURE tables */
  GetGainCurveInfo (inData, outData, err);     /* GAIN_CURVE tables */
  GetInterferometerModelInfo (inData, outData, err); /* INTERFEROMETER_MODEL tables */
  GetPhaseCalInfo (inData, outData, err);      /* PHASE_CAL tables */
  GetWeatherInfo (inData, outData, err);       /* WEATHER tables */

  /* History */
  IDIInHistory (inData, myInput, outData, err);
  
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
  if (SourceID) g_free(SourceID);
 
  return ierr;
} /* end of main */

ObitInfoList* IDIInin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="IDIIn.in", *arg;
  gboolean init=FALSE;
  oint itemp;
  gchar *strTemp;
  ObitInfoList* list;
  gchar *routine = "IDIInin";

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
      
    } else if (strcmp(arg, "-Scan") == 0){ /* Scan name */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "Scan", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-DataRoot") == 0){ /* Data root directory */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "DataRoot", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outFile") == 0){ /* Output FITS file */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outFile", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outDisk") == 0) { /* Output disk */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outDisk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-outName") == 0) { /* AIPS UV outName */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outName", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-outClass") == 0) { /* AIPS UV outClass */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outClass", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outSeq") == 0) { /* AIPS output UV sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outSeq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-pgmNumber") == 0) { /*Program number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "pgmNumber", OBIT_oint, dim, &itemp, err);
      
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
} /* end IDIInin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: IDIIn -input file -output ofile [args]\n");
    fprintf(stderr, "Convert an IDI/IDI file format to Obit/UV\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def IDIIn.in\n");
    fprintf(stderr, "  -output output result file, def UVSub.out\n");
    fprintf(stderr, "  -scan of file root used to form scan FITS file names\n");
    fprintf(stderr, "  -DataRoot Directory name for input \n");
    fprintf(stderr, "  -outFile output uv FITS  file\n");  
    fprintf(stderr, "  -outName output uv AIPS file name\n");
    fprintf(stderr, "  -outClass output uv AIPS file class\n");
    fprintf(stderr, "  -outSeq output uv AIPS file sequence\n");
    fprintf(stderr, "  -outDisk output uv (AIPS or FITS) disk number (1-rel) \n");
    
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
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  gchar *scan_name ="unspecified";
  ObitInfoList *out = newObitInfoList();
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
  ObitInfoListPut (out, "Scan", OBIT_string, dim, scan_name, err);

  /* output FITS file name */
  strTemp = "uvData.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);

  /* root of data directory */
  strTemp = "dataRoot";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "DataRoot", OBIT_string, dim, strTemp, err);

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
/*  Create output uv data                                                 */
/*   Input:                                                               */
/*      Source    Source name                                             */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/* Returns the output uv data                                             */
/*----------------------------------------------------------------------- */
ObitUV* setOutputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitUV    *outUV = NULL;
  ObitInfoType type;
  olong      i, n, Aseq, disk, cno, lType;
  gchar     *Type, *strTemp, outFile[129];
  gchar     Aname[13], Aclass[7], *Atype = "UV";
  olong      nvis;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean  exist;
  gchar     tname[129], *fullname=NULL;
  gchar     *routine = "setOutputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic output UV Object */
  g_snprintf (tname, 100, "output UV data");
  outUV = newObitUV(tname);
    
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  lType = dim[0];
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */

    /* outName given? */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&strTemp);
    for (i=0; i<12; i++) Aname[i] = ' ';  Aname[i] = 0;
    for (i=0; i<MIN(12,dim[0]); i++) Aname[i] = strTemp[i];
    /* Save any defaulting on myInput */
    dim[0] = 12;
    ObitInfoListAlwaysPut (myInput, "outName", OBIT_string, dim, Aname);

      
    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    /* Default out class is "IDIIn" */
    if (!strncmp(Aclass, "      ", 6)) strncpy (Aclass, "IDIIn", 7);

    /* input AIPS disk - default is outDisk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (disk<=0)
       ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* output AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Allocate catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

    /* Did it previously exist? */
    isNew = !exist;
    
    /* define object */
    nvis = 1;
    ObitUVSetAIPS (outUV, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */

    /* outFile given? */
    ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&strTemp);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) outFile[i] = strTemp[i]; outFile[i] = 0;
    ObitTrimTrail(outFile);  /* remove trailing blanks */

    /* Save any defaulting on myInput */
    dim[0] = strlen(outFile);
    ObitInfoListAlwaysPut (myInput, "outFile", OBIT_string, dim, outFile);

    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (disk<=0) /* defaults to outDisk */
      ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);

    /* Did it previously exist */
    disk = 1;  /* input "disk" */
    fullname = ObitFITSFilename (disk, outFile, err);
    exist =  ObitFileExist (fullname, err);
    if (fullname) g_free(fullname);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    isNew = !exist;
    
    /* define object */
    nvis = 1;
    ObitUVSetFITS (outUV, nvis, disk, outFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outUV;
  }
  
  return outUV;
} /* end setOutputUV */

void GetHeader (ObitUV *outData, gchar *inscan, 
		ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get header information from scan header files                         */
/*  Returns TRUE if the file is just created                              */
/*   Input:                                                               */
/*      outData  Output UV object                                         */
/*      inscan   root of input file names                                 */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitUVDesc *desc;
  ObitTableIDI_UV_DATA *inTable=NULL;
  ObitData *inData=NULL;
  olong ncol;
  gchar FullFile[128], *today=NULL;
  olong disk, lim;
  oint no_band=0, maxis1=0, maxis2=0, maxis3=0, maxis4=0, maxis5=0;
  ObitIOAccess access;
  olong ver;
  gchar *routine = "GetHeader";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outData));
  g_assert(inscan!=NULL);
  g_assert(myInput!=NULL);

  /* Full input file name */
  sprintf (FullFile,"%s%s.fits", DataRoot, inscan);

  /* Create input Data from which to read tables */
  inData = newObitData("Input Data");
  disk = 0;
  ObitDataSetFITS(inData, disk, FullFile, err);
  /* Open and close to init TableList */
  ObitDataOpen (inData, OBIT_IO_ReadOnly, err);
  ObitDataClose (inData,  err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);


  /* Get source info, copy to output SU table, save lookup table in SourceID */
  ObitUVOpen (outData, OBIT_IO_ReadWrite, err);
  GetSourceInfo (inData, outData, isNew, err);
  ObitUVClose (outData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
 
  /* Define output descriptor if isNew */
  if (isNew) {
    desc = outData->myDesc;
    /* Create input Antenna table object */
    ver = 1;
    access = OBIT_IO_ReadOnly;
    inTable = newObitTableIDI_UV_DATAValue ("Input table", inData, 
					    &ver, access, no_band, 
					    maxis1, maxis2, maxis3, maxis4, maxis5,
					    err);
   if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_UV_DATA table");
   if (err->error) Obit_traceback_msg (err, routine, inData->name);

    /* Open table */
    if ((ObitTableIDI_UV_DATAOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_UV_DATA table");
      return;
    }

    /* Define header */
    desc->nvis = 0;
    strncpy (desc->origin, "Obit ", UVLEN_VALUE);
    desc->isort[0] = inTable->sort[0]; 
    desc->isort[1] = inTable->sort[1]; 
    desc->nrparm = 0;
    
    /* Creation date today */
    today = ObitToday();
    strncpy (desc->date, today, UVLEN_VALUE-1);
    if (today) g_free(today);
    desc->JDObs = 0.0;
    desc->epoch = 0.0;
    desc->equinox = 0.0;
    desc->obsra = 0.0;
    desc->obsdec = 0.0;
    desc->altCrpix = 0.0;
    desc->altRef = 0.0;
    desc->restFreq = 0.0;
    desc->xshift = 0.0;
    desc->yshift = 0.0;
    desc->VelDef = 0;
    desc->VelReference = 0;
    strncpy (desc->bunit, "        ", UVLEN_VALUE);
    lim = MIN (UVLEN_VALUE, MAXKEYCHARTABLEIDI_UV_DATA);
    strncpy (desc->obsdat,     inTable->RefDate, lim);
    strncpy (desc->teles,      inTable->teles, lim);
    strncpy (desc->observer,   inTable->observer, lim);
    strncpy (desc->instrument, inTable->ArrName, lim);
   
    /* Random parameters */
    ncol = 0;
    /* U */
    strncpy (desc->ptype[ncol], "UU-L-SIN", UVLEN_KEYWORD);
    ncol++;

    /* V */
    strncpy (desc->ptype[ncol], "VV-L-SIN", UVLEN_KEYWORD);
    ncol++;
    
    /* W */
    strncpy (desc->ptype[ncol], "WW-L-SIN", UVLEN_KEYWORD);
    ncol++;
    
    /* Baseline */
    strncpy (desc->ptype[ncol], "BASELINE", UVLEN_KEYWORD);
    ncol++;
    
    /* Time */
    strncpy (desc->ptype[ncol], "TIME1   ", UVLEN_KEYWORD);
    ncol++;
    
    /* Source - must have input SOURCE table */
    if ((inTable->myDesc->repeat[inTable->SourceCol]>0) && SourceID) {
      strncpy (desc->ptype[ncol], "SOURCE  ", UVLEN_KEYWORD);
      ncol++;
    }
    
    /* FreqID */
    if (inTable->myDesc->repeat[inTable->FreqIDCol]>0) {
      strncpy (desc->ptype[ncol], "FREQSEL ", UVLEN_KEYWORD);
      ncol++;
    }
    
    /* Integration time */
    if (inTable->myDesc->repeat[inTable->IntTimCol]>0) {
      strncpy (desc->ptype[ncol], "INTTIM  ", UVLEN_KEYWORD);
      ncol++;
    }
    desc->nrparm = ncol;  /* Number of random parameters */
    
    /* Data Matrix */
    ncol = 0;
    lim = MIN (UVLEN_KEYWORD, MAXKEYCHARTABLEIDI_UV_DATA);
    /* Dimension 1 */
    if (inTable->maxis>=1) {
      strncpy (desc->ctype[ncol], inTable->ctype1, lim);
      desc->inaxes[ncol] = MAX (inTable->maxis1, 3);
      desc->cdelt[ncol] = inTable->cdelt1;
      desc->crpix[ncol] = inTable->crpix1;
      desc->crval[ncol] = inTable->crval1;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    /* Dimension 2 */
    if (inTable->maxis>=2) {
      strncpy (desc->ctype[ncol], inTable->ctype2, lim);
      desc->inaxes[ncol] = inTable->maxis2;
      desc->cdelt[ncol] = inTable->cdelt2;
      desc->crpix[ncol] = inTable->crpix2;
      desc->crval[ncol] = inTable->crval2;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    /* Dimension 3 */
    if (inTable->maxis>=3) {
      strncpy (desc->ctype[ncol], inTable->ctype3, lim);
      desc->inaxes[ncol] = inTable->maxis3;
      desc->cdelt[ncol] = inTable->cdelt3;
      desc->crpix[ncol] = inTable->crpix3;
      desc->crval[ncol] = inTable->crval3;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    /* Dimension 4 */
    if (inTable->maxis>=4) {
      strncpy (desc->ctype[ncol], inTable->ctype4, lim);
      desc->inaxes[ncol] = inTable->maxis4;
      desc->cdelt[ncol] = inTable->cdelt4;
      desc->crpix[ncol] = inTable->crpix4;
      desc->crval[ncol] = inTable->crval4;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    /* Dimension 5 */
    if (inTable->maxis>=5) {
      strncpy (desc->ctype[ncol], inTable->ctype5, lim);
      desc->inaxes[ncol] = inTable->maxis5;
      desc->cdelt[ncol] = inTable->cdelt5;
      desc->crpix[ncol] = inTable->crpix5;
      desc->crval[ncol] = inTable->crval5;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    /* Dimension 6 */
    if (inTable->maxis>=6) {
      strncpy (desc->ctype[ncol], inTable->ctype6, lim);
      desc->inaxes[ncol] = inTable->maxis6;
      desc->cdelt[ncol] = inTable->cdelt6;
      desc->crpix[ncol] = inTable->crpix6;
      desc->crval[ncol] = inTable->crval6;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    /* Dimension 7 */
    if (inTable->maxis>=7) {
      strncpy (desc->ctype[ncol], inTable->ctype7, lim);
      desc->inaxes[ncol] = inTable->maxis7;
      desc->cdelt[ncol] = inTable->cdelt7;
      desc->crpix[ncol] = inTable->crpix7;
      desc->crval[ncol] = inTable->crval7;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    desc->naxis = ncol;  /* Number of dimensions */

    /* Close  table */
    if ((ObitTableIDI_UV_DATAClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_UV_DATA Table file");
      return;
    }

    /* index descriptor */
    ObitUVDescIndex (desc);

    /* Add Antenna and Frequency info */
    ObitUVOpen (outData, OBIT_IO_WriteOnly, err) ;
    GetFrequencyInfo (inData, outData, err);
    GetAntennaInfo (inData, outData, err);
    ObitUVClose (outData, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);

  } /* End define descriptor */

  /* Instantiate output Data */
  ObitUVFullInstantiate (outData, FALSE, err);
  if (err->error)Obit_traceback_msg (err, routine, outData->name);


  /* Cleanup */
  inData  = ObitDataUnref(inData);
  inTable = ObitTableIDI_UV_DATAUnref(inTable);

  return;
} /* end GetHeader */

void GetData (ObitUV *outData, gchar *inscan, ObitInfoList *myInput, 
	      ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Read data from IDI file, write outData                                */
/*      outData  Output OTF object, open on input                         */
/*      inscan   Scan part of input file name                             */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  olong iRow, kinc, lim, ver, ant1, ant2, tant1, tant2;
  olong disk, i;
  oint no_band=0, maxis1=0, maxis2=0, maxis3=0, maxis4=0, maxis5=0;
  ofloat lambdaPerSec, sid;
  ofloat *Buffer=NULL;
  odouble JD;
  ObitUVDesc *desc;
  ObitIOAccess access;
  ObitData *inData=NULL;
  ObitTableIDI_UV_DATA *inTable=NULL;
  ObitTableIDI_UV_DATARow *inRow=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar FullFile[128];
  gchar tstr1[32], tstr2[32];
  gchar *routine = "GetData";

  /* error checks */
  if (err->error) return;
  g_assert(inscan!=NULL);
  g_assert(myInput!=NULL);
  g_assert(ObitUVIsA(outData));

  /* Full input file name */
  sprintf (FullFile,"%s%s.fits", DataRoot, inscan);

  /* Create input Data from which to read tables */
  inData = newObitData("Input Data");
  disk = 0;
  ObitDataSetFITS(inData, disk, FullFile, err);
  /* Open and close to init TableList */
  ObitDataOpen (inData, OBIT_IO_ReadOnly, err);
  ObitDataClose (inData,  err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  ver = 1;
  access = OBIT_IO_ReadOnly;
  inTable = newObitTableIDI_UV_DATAValue ("Input table", inData, 
					  &ver, access, no_band, 
					  maxis1, maxis2, maxis3, maxis4, maxis5,
					  err);
  if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_UV_DATA table");
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Open table */
  if ((ObitTableIDI_UV_DATAOpen (inTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_UV_DATA table");
    return;
  }

  /* Create Row */
  inRow = newObitTableIDI_UV_DATARow (inTable);

  /* Consistency checks */
  desc = outData->myDesc;
  if (desc->jlocs>=0) nstok = desc->inaxes[desc->jlocs];
  if (desc->jlocf>=0) nchan = desc->inaxes[desc->jlocf];
  if (desc->jlocif>=0) nIF  = desc->inaxes[desc->jlocif];

  Obit_return_if_fail((nchan==inTable->no_chan), err,
		       "%s: Input number freq. incompatible %d != %d", 
		      routine, nchan, inTable->no_chan);
  Obit_return_if_fail((nstok==inTable->no_stkd), err,
		       "%s: Input number Poln incompatible %d != %d", 
		      routine, nstok, inTable->no_stkd);
  Obit_return_if_fail((nIF==inTable->no_band), err,
		       "%s: Input number Bands (IFs) incompatible %d != %d", 
		      routine, nIF, inTable->no_band);
  /* Truncate strings to last non blank */
  strncpy(tstr1, inTable->ctype1, 32);
  ObitTrimTrail(tstr1);
  strncpy(tstr2, desc->ctype[0], 32);
  ObitTrimTrail(tstr2);
  lim = MIN (UVLEN_KEYWORD, MAXKEYCHARTABLEIDI_UV_DATA);
  Obit_return_if_fail((!strncmp(tstr1, tstr2, lim)), err,
		       "%s: First axis different %s != %s", 
		      routine, tstr1, tstr2);  
  /* Truncate strings to last non blank */
  strncpy(tstr1, inTable->ctype2, 32);
  ObitTrimTrail(tstr1);
  strncpy(tstr2, desc->ctype[1], 32);
  ObitTrimTrail(tstr2);
  Obit_return_if_fail((!strncmp(tstr1, tstr2, lim)), err,
		      "%s: Second axis different %s != %s", 
		      routine, tstr1, tstr2);  
  if (inTable->maxis>=2) {
    /* Truncate strings to last non blank */
    strncpy(tstr1, inTable->ctype3, 32);
    ObitTrimTrail(tstr1);
    strncpy(tstr2, desc->ctype[2], 32);
    ObitTrimTrail(tstr2);
    Obit_return_if_fail((!strncmp(tstr1, tstr2, lim)), err,
			"%s: Third axis different %s != %s", 
			routine, tstr1, tstr2);  
  }
  if (inTable->maxis>=3) {
    /* Truncate strings to last non blank */
    strncpy(tstr1, inTable->ctype4, 32);
    ObitTrimTrail(tstr1);
    strncpy(tstr2, desc->ctype[3], 32);
    ObitTrimTrail(tstr2);
    Obit_return_if_fail((!strncmp(tstr1, tstr2, lim)), err,
			"%s: Fourth axis different %s != %s", 
			routine, tstr1, tstr2);  
  }
  if (inTable->maxis>=4) {
    /* Truncate strings to last non blank */
    strncpy(tstr1, inTable->ctype5, 32);
    ObitTrimTrail(tstr1);
    strncpy(tstr2, desc->ctype[4], 32);
    ObitTrimTrail(tstr2);
    Obit_return_if_fail((!strncmp(tstr1, tstr2, lim)), err,
			"%s: Fifth axis different %s != %s", 
			routine, tstr1, tstr2);  
  }
  if (inTable->maxis>=5) {
    /* Truncate strings to last non blank */
    strncpy(tstr1, inTable->ctype6, 32);
    ObitTrimTrail(tstr1);
    strncpy(tstr2, desc->ctype[5], 32);
    ObitTrimTrail(tstr2);
    Obit_return_if_fail((!strncmp(tstr1, tstr2, lim)), err,
			"%s: Sixth axis different %s != %s", 
			routine, tstr1, tstr2);  
  }
  if (inTable->maxis>=6) {
    /* Truncate strings to last non blank */
    strncpy(tstr1, inTable->ctype7, 32);
    ObitTrimTrail(tstr1);
    strncpy(tstr2, desc->ctype[6], 32);
    ObitTrimTrail(tstr2);
    Obit_return_if_fail((!strncmp(tstr1, tstr2, lim)), err,
			"%s: Seventh axis different %s != %s", 
			routine, tstr1, tstr2);  
  }
  Obit_return_if_fail((desc->ncorr*3==inTable->myDesc->repeat[inTable->FluxCol]), err,
		      "%s: Input and output data sizes different, %d !=  %d", 
		      routine, desc->ncorr*3, 
		      inTable->myDesc->repeat[inTable->FluxCol]);

  /* Prepare output */
  desc = outData->myDesc;
  Buffer = outData->buffer;
  desc->firstVis = desc->nvis+1;  /* Append to end of data */

  /* Get any prior Obit sort order */
  if (outData->myDesc->isort[0] == ' ') {
    ObitInfoListGetTest(myInput, "OBITSORT", &type, dim, tstr1);
    /* Only the first 2 characters */
    outData->myDesc->isort[0] = tstr1[0];
    outData->myDesc->isort[1] = tstr1[1];
    outData->myDesc->isort[2] = 0;
  }
  
  /* Get reference MJD, convert ref date  to JD */
  ObitUVDescDate2JD (inTable->RefDate, &JD);
  /* Reference Modified Julian date */
  refMJD = JD - 2400000.5;

  /* Number of values on axis 1 */
  kinc = inTable->maxis1;

  /* Number of wavelengths per second */
  lambdaPerSec = inTable->ref_freq;

  /* Loop over table */
  for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {

    if ((ObitTableIDI_UV_DATAReadRow (inTable, iRow, inRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR reading IDI_UV_DATA Table");
      return;
      }
    if (inRow->status==-1) continue;  /* Flagged entry */

    /* Antenna numbers in proper order */
    tant1 = (olong)(inRow->Baseline/256);
    tant2 = (olong)(inRow->Baseline - tant1*256);
    ant1 = MIN (tant1, tant2);
    ant2 = MAX (tant1, tant2);
 
    /* Convert table to UV data form */
    /* Convert U,V,W to lambda */
    Buffer[desc->ilocu] = (ofloat)inRow->uu*lambdaPerSec;
    Buffer[desc->ilocv] = (ofloat)inRow->vv*lambdaPerSec;
    Buffer[desc->ilocw] = (ofloat)inRow->ww*lambdaPerSec;
    Buffer[desc->ilocb] = (ofloat)(ant1*256+ant2);
    Buffer[desc->iloct] = (ofloat)(inRow->date+inRow->Time-JD);
    if (desc->ilocsu>=0) {
      if (SourceID) sid = SourceID[inRow->Source-1];
      else sid = inRow->Source;
      Buffer[desc->ilocsu] = sid;
    }
    if (desc->ilocfq>=0) Buffer[desc->ilocfq] = (ofloat)inRow->FreqID;
    if (desc->ilocit>=0) Buffer[desc->ilocit] = (ofloat)inRow->IntTim;
    for (i=0; i<desc->ncorr; i++) {
      Buffer[desc->nrparm+i*3]   = inRow->Flux[i*kinc];
      Buffer[desc->nrparm+i*3+1] = inRow->Flux[i*kinc+1];
      Buffer[desc->nrparm+i*3+2] = inRow->Weight[i];
    }

    /* set number of records */
    desc->numVisBuff = 1;
    
    /* Write output */
    if ((ObitUVWrite (outData, NULL, err) != OBIT_IO_OK) || (err->error))
      Obit_log_error(err, OBIT_Error, "ERROR writing output UV data"); 
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
  } /* end loop over table */  


  /* Close  tables */
  if ((ObitTableIDI_UV_DATAClose (inTable, err) 
       != OBIT_IO_OK) || (err->error)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_UV_DATA Table");
    return;
  }

  /* Cleanup */
  inData  = ObitDataUnref(inData);
  inTable = ObitTableIDI_UV_DATAUnref(inTable);
  inRow   = ObitTableIDI_UV_DATAUnref(inTable);

} /* end GetData  */

void GetAntennaInfo (ObitData *inData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get info from IDI_ANTENNA and IDI_ARRAY_GEOMETRY tables on indata     */
/*  Writes AN table output                                                */
/*  to output .                                                           */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_ANTENNA         *inTable=NULL;
  ObitTableIDI_ANTENNARow      *inRow=NULL;
  ObitTableIDI_ARRAY_GEOMETRY  *in2Table=NULL;
  ObitTableIDI_ARRAY_GEOMETRYRow  *in2Row=NULL;
  ObitTableAN                  *outTable=NULL;
  ObitTableANRow               *outRow=NULL;
  olong i, lim, iRow,i2Row, oRow, ver;
  oint numIF, numPCal, numOrb;
  gboolean found;
  ObitIOAccess access;
  gchar *routine = "GetAntennaInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Create input Antenna table object */
  ver = 1;
  access = OBIT_IO_ReadOnly;
  numIF =  numPCal = 0;
  inTable = newObitTableIDI_ANTENNAValue ("Input table", inData, 
					  &ver, access, numIF, numPCal, err);
  if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_ANTENNA table");
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Open table */
  if ((ObitTableIDI_ANTENNAOpen (inTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_ANTENNA table");
    return;
  }

  /* Create Row */
  inRow = newObitTableIDI_ANTENNARow (inTable);

  /* Create input Array geometry table object */
  ver = 1;
  access = OBIT_IO_ReadOnly;
  numIF = numOrb = 0;
  in2Table = newObitTableIDI_ARRAY_GEOMETRYValue ("Input table", inData, 
						  &ver, access, numIF, numOrb, err);
  if (in2Table==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_ARRAY_GEOMETRY table");
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Open table */
  if ((ObitTableIDI_ARRAY_GEOMETRYOpen (in2Table, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_ARRAY_GEOMETRY table");
    return;
  }

  /* Create Row */
  in2Row = newObitTableIDI_ARRAY_GEOMETRYRow (in2Table);

  /* Create output Antenna table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  numOrb  = in2Table->numOrb;
  numPCal = inTable->numPCal;
  outTable = newObitTableANValue ("Output table", (ObitData*)outData, 
				  &ver, access, numOrb, numPCal, err);
  if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open table */
  if ((ObitTableANOpen (outTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output AN table");
    return;
  }

  /* Create output Row */
  outRow = newObitTableANRow (outTable);
  /* attach to table buffer */
  ObitTableANSetRow (outTable, outRow, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Set AN table values */
  outTable->ArrayX = in2Table->ArrayX;
  outTable->ArrayY = in2Table->ArrayY;
  outTable->ArrayZ = in2Table->ArrayZ;
  outTable->GSTiat0 = in2Table->GSTiat0;
  outTable->DegDay  = in2Table->DegDay;
  outTable->Freq    = in2Table->ref_freq;
  lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, MAXKEYCHARTABLEAN);
  strncpy (outTable->RefDate, inTable->RefDate, lim);
  outTable->PolarX  = in2Table->PolarX;
  outTable->PolarY  = in2Table->PolarY;
  outTable->dataUtc = in2Table->ut1Utc;
  lim = MIN (MAXKEYCHARTABLEIDI_ARRAY_GEOMETRY, MAXKEYCHARTABLEAN);
  strncpy (outTable->TimeSys, in2Table->TimeSys, lim );
  outTable->FreqID = 0;
  outTable->iatUtc = in2Table->ut1Utc;
  lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, MAXKEYCHARTABLEAN);
  strncpy (outTable->polType, inTable->polType, lim );
  outTable->P_Refant = 0;
  outTable->P_Diff01 = 0.0;
  outTable->P_Diff02 = 0.0;
  outTable->P_Diff03 = 0.0;
  outTable->P_Diff04 = 0.0;
  outTable->P_Diff05 = 0.0;
  outTable->P_Diff06 = 0.0;
  outTable->P_Diff07 = 0.0;
  outTable->P_Diff08 = 0.0;

  /* Initialize output row */
  outRow->noSta     = 0;
  outRow->mntSta    = 0;
  outRow->staXof    = 0.0;
  outRow->PolAngA   = 0.0;
  outRow->PolAngB   = 0.0;
  outRow->AntName[0]= 0; 
  outRow->StaXYZ[0] = 0.0;
  outRow->StaXYZ[1 ]= 0.0;
  outRow->StaXYZ[2] = 0.0;
  outRow->OrbParm[0]= 0.0;
  outRow->polTypeA[0] = ' ';
  outRow->PolCalA[0]  = 0.0;
  outRow->polTypeB[0] =  ' ';
  outRow->PolCalB[0]  = 0.0;
  outRow->status      = 0;

  /* loop through input Antenna table */
  for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
    if ((ObitTableIDI_ANTENNAReadRow (inTable, iRow, inRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR reading IDI_ANTENNA Table");
      return;
    }

    /* Find in Array Geometry table */
    found = FALSE;
    for (i2Row = 1; i2Row<=in2Table->myDesc->nrow; i2Row++) {
      if ((ObitTableIDI_ARRAY_GEOMETRYReadRow (in2Table, i2Row, in2Row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading IDI_ARRAY_GEOMETRYR Table");
	return;
      }
    if (inRow->antenna_no==in2Row->noSta) {
      /* Found match */
      found = TRUE;
      break;
      }  
    } /* end loop over table */

    /* Set output Row */
    outRow->noSta     = inRow->antenna_no;
    outRow->PolAngA   = inRow->PolAngA[0];
    outRow->polTypeA[0]  = inRow->polTypeA;   
    for (i=0; i<numPCal; i++) 
	outRow->PolCalA[i] = inRow->PolCalA[i];
    if (inTable->no_stkd>1) {  /* Multiple polarizations */
      outRow->PolAngB   = inRow->PolAngB[0];   
      outRow->polTypeB[0]  = inRow->polTypeB;   
      for (i=0; i<numPCal; i++) 
	outRow->PolCalB[i] = inRow->PolCalB[i];
    }
    lim = MIN(inTable->myDesc->repeat[inTable->AntNameCol], 
	      outTable->myDesc->repeat[outTable->AntNameCol]);
    for (i=0; i<lim; i++) 
	outRow->AntName[i] = inRow->AntName[i];
    outRow->status      = 0;
    if (found) {
      outRow->mntSta    = in2Row->mntSta;
      outRow->staXof    = in2Row->staXof[0];
      outRow->StaXYZ[0] = in2Row->StaXYZ[0];
      outRow->StaXYZ[1] = in2Row->StaXYZ[1];
      outRow->StaXYZ[2] = in2Row->StaXYZ[2];
      for (i=0; i<numOrb; i++) outRow->OrbParm[i] = in2Row->OrbParm[i];
    }
    oRow = outRow->noSta;
    if ((ObitTableANWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating ANTENNA Table");
      return;
    }
  } /* end loop over input table */
  
  /* Close  tables */
  if ((ObitTableIDI_ANTENNAClose (inTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_ANTENNA Table file");
    return;
  }

  if ((ObitTableIDI_ARRAY_GEOMETRYClose (in2Table, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_ARRAY_GEOMETRY Table file");
    return;
  }

  if ((ObitTableANClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output Antenna Table file");
    return;
  }

  /* Cleanup */
  inRow    = ObitTableIDI_ANTENNARowUnref(inRow);
  inTable  = ObitTableIDI_ANTENNAUnref(inTable);
  in2Row   = ObitTableIDI_ARRAY_GEOMETRYRowUnref(in2Row);
  in2Table = ObitTableIDI_ARRAY_GEOMETRYUnref(in2Table);
  outRow   = ObitTableANRowUnref(outRow);
  outTable = ObitTableANUnref(outTable);

} /* end  GetAntennaInfo */

void GetFrequencyInfo (ObitData *inData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get Source info from IDI_FREQUENCY table on indata                    */
/*  Write FQ table                                                        */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_FREQUENCY*    inTable=NULL;
  ObitTableIDI_FREQUENCYRow* inRow=NULL;
  ObitTableFQ*            outTable=NULL;
  ObitTableFQRow*         outRow=NULL;
  olong i, iRow, oRow, ver;
  oint numIF;
  ObitIOAccess access;
  gchar *routine = "GetFrequencyInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Create input FREQUENCY table object */
  ver = 1;
  access = OBIT_IO_ReadOnly;
  numIF = 0;
  inTable = newObitTableIDI_FREQUENCYValue ("Input table", inData, 
					 &ver, access, numIF, err);
  if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_FREQUENCY table");
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Open table */
  if ((ObitTableIDI_FREQUENCYOpen (inTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_FREQUENCY table");
    return;
  }

  /* Create Row */
  inRow = newObitTableIDI_FREQUENCYRow (inTable);

  /* Create output FQ table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  numIF = inTable->no_band;
  outTable = newObitTableFQValue ("Output table", (ObitData*)outData, 
				  &ver, access, numIF, err);
  if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with FQ table");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open table */
  if ((ObitTableFQOpen (outTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FQ table");
    return;
  }

  /* Create output Row */
  outRow = newObitTableFQRow (outTable);
  /* attach to table buffer */
  ObitTableFQSetRow (outTable, outRow, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Initialize output row */
  outRow->fqid    = 0;
  for (i=0; i<numIF; i++) {
    outRow->freqOff[i]  = 0;
    outRow->chWidth[i]  = 0;
    outRow->totBW[i]    = 0;
    outRow->sideBand[i] = 0;
  }
  outRow->status    = 0;

  /* loop through input table */
  for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
    if ((ObitTableIDI_FREQUENCYReadRow (inTable, iRow, inRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR reading IDI_FREQUENCY Table");
      return;
    }

    /* Save to FQ table */
    outRow->fqid    = inRow->fqid;
    for (i=0; i<numIF; i++) {
      outRow->freqOff[i]  = inRow->bandfreq[i];
      outRow->chWidth[i]  = inRow->chWidth[i];
      outRow->totBW[i]    = inRow->totBW[i];
      outRow->sideBand[i] = inRow->sideBand[i];
    }
    outRow->status    = 0;
    oRow = outRow->fqid;
    if (oRow<1) oRow = -1;
    if ((ObitTableFQWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating FG Table");
      return;
    }
  } /* end loop over input table */
  
  /* Close  tables */
  if ((ObitTableIDI_FREQUENCYClose (inTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_FREQUENCY Table file");
    return;
  }

  if ((ObitTableFQClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output FQ Table file");
    return;
  }

  /* Cleanup */
  inRow    = ObitTableIDI_FREQUENCYRowUnref(inRow);
  inTable  = ObitTableIDI_FREQUENCYUnref(inTable);
  outRow   = ObitTableFQRowUnref(outRow);
  outTable = ObitTableFQUnref(outTable);

} /* end  GetFrequencyInfo */

void GetSourceInfo (ObitData *inData, ObitUV *outData, gboolean isNew, 
		    ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get Source info from IDI_SOURCE table on inData                       */
/*  Copies source info from intput to putput data and generated a lookup  */
/*  table, global, SourceID to give translation from input source ID      */
/*  to output.  If no SOURCE table leafe SourceID=NULL                    */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*      isNew    True if output file just created                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_SOURCE*    inTable=NULL;
  ObitTableIDI_SOURCERow* inRow=NULL;
  ObitTableSU*            outTable=NULL;
  ObitTableSURow*         outRow=NULL;
  olong i, lim, iRow, oRow, ver;
  oint numIF;
  gboolean found;
  ObitIOAccess access;
  gchar *routine = "GetSourceInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Print any messages */
  ObitErrLog(err);
  
  /* Create input Source table object */
  ver = 1;
  access = OBIT_IO_ReadOnly;
  numIF = 0;
  inTable = newObitTableIDI_SOURCEValue ("Input table", inData, 
					 &ver, access, numIF, err);
    /* Find it?  If not, just return */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      return;
    }
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Open table */
  if ((ObitTableIDI_SOURCEOpen (inTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_SOURCE table");
    return;
  }

  /* Create Row */
  inRow = newObitTableIDI_SOURCERow (inTable);

  /* Create global source ID lookup table */
  SourceID = g_malloc0((2+inTable->myDesc->nrow)*sizeof(ofloat));

  /* Create output Source table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  numIF = inTable->no_band;
  outTable = newObitTableSUValue ("Output table", (ObitData*)outData, 
				  &ver, access, numIF, err);
  if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with SU table");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open table */
  if ((ObitTableSUOpen (outTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output SU table");
    return;
  }

  /* Create output Row */
  outRow = newObitTableSURow (outTable);
  /* attach to table buffer */
  ObitTableSUSetRow (outTable, outRow, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Initialize output row */
  outRow->SourID    = 0;
  outRow->Qual      = 0;
  outRow->Bandwidth = 0.0;
  outRow->RAMean    = 0.0;
  outRow->DecMean   = 0.0;
  outRow->Epoch     = 0.0;
  outRow->RAApp     = 0.0;
  outRow->DecApp    = 0.0;
  outRow->PMRa      = 0.0;
  outRow->PMDec     = 0.0;
  outRow->Source[0] = 0;
  outRow->CalCode[0]= 0;
  for (i=0; i<inTable->no_band; i++) {
    outRow->IFlux[i]     = 0.0;
    outRow->QFlux[i]     = 0.0;
    outRow->UFlux[i]     = 0.0;
    outRow->VFlux[i]     = 0.0;
    outRow->FreqOff[i]   = 0.0;
    outRow->LSRVel[i]    = 0.0;
    outRow->RestFreq [i] = 0.0;
  }
  outRow->status    = 0;


  /* loop through input table */
  for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
    if ((ObitTableIDI_SOURCEReadRow (inTable, iRow, inRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR reading IDI_SOURCE Table");
      return;
    }

    /* See if source exists in output table */
    found = FALSE;
    for (oRow = 1; oRow<=outTable->myDesc->nrow; oRow++) {
      if ((ObitTableSUReadRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading Source Table");
	return;
      }
    if (!strncmp (inRow->Source, outRow->Source, 16)) {
      /* Found match */
      found = TRUE;
      break;
      }  
    } /* end loop over table */

    /* If found just remember output source ID */
    if (found) {
      /* Save source no. in globla lookup table */
      SourceID[inRow->SourID] = outRow->SourID;
    } else { /* Not found - add */
      /* If first time update header */
      if (isNew) {
	outTable->FreqID = inRow->FreqID;
	lim = MIN(inTable->myDesc->repeat[inTable->SysVelCol], MAXKEYCHARTABLESU);
	for (i=0; i<lim; i++) outTable->velDef[i] = inRow->SysVel[i];
	lim = MIN(inTable->myDesc->repeat[inTable->VelTypCol], MAXKEYCHARTABLESU);
	for (i=0; i<lim; i++) outTable->velType[i] = inRow->VelTyp[i];
      }
      /* Set output row for end of table */
      outRow->SourID    = outTable->myDesc->nrow+1;
      outRow->Qual      = inRow->Qual;
      outRow->RAMean    = inRow->RAMean;
      outRow->DecMean   = inRow->DecMean;
      outRow->Epoch     = inRow->Epoch;
      outRow->RAApp     = inRow->RAApp;
      outRow->DecApp    = inRow->DecApp;
      outRow->PMRa      = inRow->PMRa;
      outRow->PMDec     = inRow->PMDec;
      lim = MIN(inTable->myDesc->repeat[inTable->SourceCol], 
		outTable->myDesc->repeat[outTable->SourceCol]);
      for (i=0; i<lim; i++) 
	outRow->Source[i] = inRow->Source[i];
      lim = MIN(inTable->myDesc->repeat[inTable->CalCodeCol], 
		outTable->myDesc->repeat[outTable->CalCodeCol]);
      for (i=0; i<lim; i++) 
	outRow->CalCode[i]= inRow->CalCode[i];
      for (i=0; i<inTable->no_band; i++) {
	outRow->IFlux[i]    = inRow->IFlux[i];
	outRow->QFlux[i]    = inRow->QFlux[i];
	outRow->UFlux[i]    = inRow->UFlux[i];
	outRow->VFlux[i]    = inRow->VFlux[i];
	outRow->FreqOff[i]  = inRow->FreqOff[i];
	outRow->LSRVel[i]   = inRow->SysVel[i];
	outRow->RestFreq[i] = inRow->RestFreq[i];
      }
      outRow->status    = 0;
 
      /* Save source no. in global lookup table */
      SourceID[inRow->SourID] = outRow->SourID;

      oRow = outRow->SourID;
      if ((ObitTableSUWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating Source Table");
	return;
      }
   } /* End add new entry */
 } /* end loop over input table */
  
  /* Close  tables */
  if ((ObitTableIDI_SOURCEClose (inTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_SOURCE Table file");
    return;
  }

  if ((ObitTableSUClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output Source Table file");
    return;
  }

  /* Cleanup */
  inRow    = ObitTableIDI_SOURCERowUnref(inRow);
  inTable  = ObitTableIDI_SOURCEUnref(inTable);
  outRow   = ObitTableSURowUnref(outRow);
  outTable = ObitTableSUUnref(outTable);

} /* end  GetSourceInfo */

void GetFlagInfo (ObitData *inData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any FLAG tables on inData to AIPS FG on outData               */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_FLAG*    inTable=NULL;
  ObitTableIDI_FLAGRow* inRow=NULL;
  ObitTableFG*          outTable=NULL;
  ObitTableFGRow*       outRow=NULL;
  olong i, iver, iRow, oRow, ver;
  oint numIF;
  ObitIOAccess access;
  gchar *routine = "GetFlagInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Loop over plausible versions */
  for (iver=1; iver<1000; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input Flag table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    numIF = 0;
    inTable = newObitTableIDI_FLAGValue ("Input table", inData, 
					 &ver, access, numIF, err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      break;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    numIF = inTable->no_band;
    
    /* Open table */
    if ((ObitTableIDI_FLAGOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_FLAG table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableIDI_FLAGRow (inTable);
    
    /* Create output FG table object */
    ver = iver;
    access = OBIT_IO_ReadWrite;
    numIF = inTable->no_band;
    outTable = newObitTableFGValue ("Output table", (ObitData*)outData, 
				    &ver, access, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with FG table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableFGOpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output FG table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTableFGRow (outTable);
    /* attach to table buffer */
    ObitTableFGSetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Initialize output row */
    outRow->status       = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableIDI_FLAGReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading IDI_FLAG Table");
	return;
      }
      
      /* Loop over bands - write one at a time in FG table */
      for (i=0; i<numIF; i++) {
	/* Flagged? */
	if (!inRow->bands[i]) continue;
	/* Save to FG table */
	outRow->SourID       = inRow->SourID;
	outRow->SubA         = inRow->Array;
	outRow->freqID       = inRow->fqid;
	outRow->ants[0]      = inRow->ants[0];
	outRow->ants[1]      = inRow->ants[1];
	outRow->TimeRange[0] = inRow->timerange[0];
	outRow->TimeRange[1] = inRow->timerange[1];
	outRow->ifs[0]       = i+1;
	outRow->ifs[1]       = i+1;
	outRow->chans[0]     = inRow->chans[0];
	outRow->chans[1]     = inRow->chans[1];
	outRow->pFlags[0]    = inRow->pflags[0];
	outRow->pFlags[1]    = inRow->pflags[1];
	outRow->pFlags[2]    = inRow->pflags[2];
	outRow->pFlags[3]    = inRow->pflags[3];
	strncpy (outRow->reason, inRow->reason ,24);
	/* Write */
	oRow = -1;
	if ((ObitTableFGWriteRow (outTable, oRow, outRow, err)
	     != OBIT_IO_OK) || (err->error>0)) { 
	  Obit_log_error(err, OBIT_Error, "ERROR updating FG Table");
	  return;
	}
      }
    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_FLAGClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_FLAG Table file");
      return;
    }
    
    if ((ObitTableFGClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output FG Table file");
      return;
    }

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied FLAG table %d", iver);

    /* Cleanup */
    inRow    = ObitTableIDI_FLAGRowUnref(inRow);
    inTable  = ObitTableIDI_FLAGUnref(inTable);
    outRow   = ObitTableFGRowUnref(outRow);
    outTable = ObitTableFGUnref(outTable);
  } /* end loop over versions */

} /* end  GetFlagInfo */

void GetCalibrationInfo (ObitData *inData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any CALIBRATION tables on inData to AIPS CL on outData        */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_CALIBRATION*    inTable=NULL;
  ObitTableIDI_CALIBRATIONRow* inRow=NULL;
  ObitTableCL*          outTable=NULL;
  ObitTableCLRow*       outRow=NULL;
  olong i, iver, iRow, oRow, ver;
  oint numIF, numAnt, numPol, numTerm;
  ObitIOAccess access;
  gchar *routine = "GetCalibrationInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Loop over plausible versions */
  for (iver=1; iver<1000; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input Flag table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    numIF = numAnt = numPol = 0;
    inTable = newObitTableIDI_CALIBRATIONValue ("Input table", inData, 
						&ver, access, numIF, numAnt, numPol, err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      break;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableIDI_CALIBRATIONOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_CALIBRATION table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableIDI_CALIBRATIONRow (inTable);
    
    /* Create output CL table object */
    ver = iver;
    access  = OBIT_IO_ReadWrite;
    numIF   = inTable->no_band;
    numAnt  = inTable->numAnt;
    numPol  = inTable->numPol;
    numTerm = 0;
    outTable = newObitTableCLValue ("Output table", (ObitData*)outData, 
				    &ver, access, numPol, numIF, numTerm, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with CL table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableCLOpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output CL table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTableCLRow (outTable);
    /* attach to table buffer */
    ObitTableCLSetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Initialize output row */
    outRow->IFR         = 0.0;
    outRow->atmos       = 0.0;
    outRow->Datmos      = 0.0;
    outRow->GeoDelay[0] = 0.0;
    outRow->MBDelay1    = 0.0;
    outRow->clock1      = 0.0;
    outRow->Dclock1     = 0.0;
    outRow->dispers1    = 0.0;
    outRow->Ddispers1   = 0.0;
    for (i=0; i<numIF; i++) outRow->DopplerOff[i] = 0.0;
    if (numPol>1) {   /* 2 poln */
      outRow->MBDelay2  = 0.0;
      outRow->clock2    = 0.0;
      outRow->Dclock2   = 0.0;
      outRow->dispers2  = 0.0;
      outRow->Ddispers2 = 0.0;
    }
    outRow->status      = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableIDI_CALIBRATIONReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading IDI_CALIBRATION Table");
	return;
      }
      
      /* Save to CL table */
      outRow->Time        = inRow->Time;
      outRow->TimeI       = inRow->TimeI;
      outRow->SourID      = inRow->SourID;
      outRow->antNo       = inRow->antNo;
      outRow->SubA        = inRow->Array;
      outRow->FreqID      = inRow->fqid;
      for (i=0; i<numIF; i++) {
	outRow->Real1[i]      = inRow->real1[i];
	outRow->Imag1[i]      = inRow->imag1[i];
	outRow->Rate1[i]      = inRow->rate1[i];
	outRow->Delay1[i]     = inRow->delay1[i];
	outRow->Weight1[i]    = inRow->weight1[i];
	outRow->RefAnt1[i]    = inRow->refant1[i];
      }
      if (numPol>1) {   /* 2 poln */
	for (i=0; i<numIF; i++) {
	  outRow->Real2[i]   = inRow->real2[i];
	  outRow->Imag2[i]   = inRow->imag2[i];
	  outRow->Rate2[i]   = inRow->rate2[i];
	  outRow->Delay2[i]  = inRow->delay2[i];
	  outRow->Weight2[i] = inRow->weight2[i];
	  outRow->RefAnt2[i] = inRow->refant2[i];
	}
      }
      /* Write */
      oRow = -1;
      if ((ObitTableCLWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating CL Table");
	return;
      }

    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_CALIBRATIONClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_CALIBRATION Table file");
      return;
    }
    
    if ((ObitTableCLClose (outTable, err) 
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

} /* end  GetCalibrationInfo */

void GetBandpassInfo (ObitData *inData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any BANDPASS tables on inData to AIPS BP on outData           */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_BANDPASS*    inTable=NULL;
  ObitTableIDI_BANDPASSRow* inRow=NULL;
  ObitTableBP*          outTable=NULL;
  ObitTableBPRow*       outRow=NULL;
  olong i, iver, iRow, oRow, ver;
  oint numIF, numAnt, numPol, numBach, strtChn;
  ObitIOAccess access;
  gchar *routine = "GetBandpassInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Loop over plausible versions */
  for (iver=1; iver<1000; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input Flag table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    numIF = numAnt = numPol = numBach =  strtChn = 0;
    inTable = newObitTableIDI_BANDPASSValue ("Input table", inData, 
					     &ver, access, 
					     numIF, numAnt, numPol, numBach, strtChn,
					     err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      break;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableIDI_BANDPASSOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_BANDPASS table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableIDI_BANDPASSRow (inTable);
    
    /* Create output BP table object */
    ver = iver;
    access  = OBIT_IO_ReadWrite;
    numIF   = inTable->no_band;
    numAnt  = inTable->numAnt;
    numPol  = inTable->numPol;
    numBach = inTable->numBach;
    outTable = newObitTableBPValue ("Output table", (ObitData*)outData, 
				    &ver, access, numPol, numIF, numBach, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with BP table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableBPOpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output BP table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTableBPRow (outTable);
    /* attach to table buffer */
    ObitTableBPSetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Initialize output row */
    for (i=0; i<numIF; i++) outRow->ChanShift[i] = 0.0;
    outRow->status      = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableIDI_BANDPASSReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading IDI_BANDPASS Table");
	return;
      }
      
      /* Save to BP table */
      outRow->Time        = inRow->Time;
      outRow->TimeI       = inRow->TimeI;
      outRow->SourID      = inRow->SourID;
      outRow->SubA        = inRow->Array;
      outRow->antNo       = inRow->antNo;
      outRow->FreqID      = inRow->fqid;
      outRow->RefAnt1    = inRow->refant1[0];
      for (i=0; i<numIF; i++) {
	outRow->Real1[i]      = inRow->breal1[i];
	outRow->Imag1[i]      = inRow->bimag1[i];
	outRow->Weight1[i]    = 1.0;
      }
      if (numPol>1) {   /* 2 poln */
	outRow->RefAnt2 = inRow->refant2[0];
	for (i=0; i<numIF; i++) {
	  outRow->Real2[i]   = inRow->breal2[i];
	  outRow->Imag2[i]   = inRow->bimag2[i];
	  outRow->Weight2[i] = 1.0;
	}
      }
      /* Write */
      oRow = -1;
      if ((ObitTableBPWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating BP Table");
	return;
      }

    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_BANDPASSClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_BANDPASS Table file");
      return;
    }
    
    if ((ObitTableBPClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output BP Table file");
      return;
    }

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied BANDPASS table %d", iver);

    /* Cleanup */
    inRow    = ObitTableIDI_BANDPASSRowUnref(inRow);
    inTable  = ObitTableIDI_BANDPASSUnref(inTable);
    outRow   = ObitTableBPRowUnref(outRow);
    outTable = ObitTableBPUnref(outTable);
  } /* end loop over versions */

} /* end  GetBandpassInfo */

void GetTSysInfo (ObitData *inData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any SYSTEM_TEMPERATURE tables on inData to AIPS TY on outData */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_SYSTEM_TEMPERATURE*    inTable=NULL;
  ObitTableIDI_SYSTEM_TEMPERATURERow* inRow=NULL;
  ObitTableTY*          outTable=NULL;
  ObitTableTYRow*       outRow=NULL;
  olong i, iver, iRow, oRow, ver;
  oint numIF, numPol;
  ObitIOAccess access;
  gchar *routine = "GetTSysInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Loop over plausible versions */
  for (iver=1; iver<1000; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input Flag table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    numIF = numPol = 0;
    inTable = newObitTableIDI_SYSTEM_TEMPERATUREValue ("Input table", inData, 
						       &ver, access, numIF, numPol, err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      break;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableIDI_SYSTEM_TEMPERATUREOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_SYSTEM_TEMPERATURE table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableIDI_SYSTEM_TEMPERATURERow (inTable);
    
    /* Create output TY table object */
    ver = iver;
    access  = OBIT_IO_ReadWrite;
    numIF   = inTable->no_band;
    numPol  = inTable->numPol;
    outTable = newObitTableTYValue ("Output table", (ObitData*)outData, 
				    &ver, access, numPol, numIF, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with TY table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableTYOpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output TY table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTableTYRow (outTable);
    /* attach to table buffer */
    ObitTableTYSetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Initialize output row */
    outRow->status      = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableIDI_SYSTEM_TEMPERATUREReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading IDI_SYSTEM_TEMPERATURE Table");
	return;
      }
      
      /* Save to TY table */
      outRow->Time      = inRow->Time;
      outRow->TimeI     = inRow->TimeI;
      outRow->SourID    = inRow->SourID;
      outRow->antennaNo = inRow->antNo;
      outRow->SubA      = inRow->Array;
      outRow->FreqID    = inRow->fqid;
      for (i=0; i<numIF; i++) {
	outRow->Tsys1[i] = inRow->TSys1[i];
	outRow->Tant1[i] = inRow->TAnt1[i];
      }
      if (numPol>1) {   /* 2 poln */
	for (i=0; i<numIF; i++) {
	  outRow->Tsys2[i] = inRow->TSys2[i];
	  outRow->Tant2[i] = inRow->TAnt2[i];
	}
      }
      /* Write */
      oRow = -1;
      if ((ObitTableTYWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating TY Table");
	return;
      }

    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_SYSTEM_TEMPERATUREClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_SYSTEM_TEMPERATURE Table file");
      return;
    }
    
    if ((ObitTableTYClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output TY Table file");
      return;
    }

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied SYSTEM_TEMPERATURE table %d", iver);

    /* Cleanup */
    inRow    = ObitTableIDI_SYSTEM_TEMPERATURERowUnref(inRow);
    inTable  = ObitTableIDI_SYSTEM_TEMPERATUREUnref(inTable);
    outRow   = ObitTableTYRowUnref(outRow);
    outTable = ObitTableTYUnref(outTable);
  } /* end loop over versions */

} /* end  GetTSysInfo */

void GetGainCurveInfo (ObitData *inData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any GAIN_CURVE tables on inData to AIPS GC on outData         */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_GAIN_CURVE*    inTable=NULL;
  ObitTableIDI_GAIN_CURVERow* inRow=NULL;
  ObitTableGC*          outTable=NULL;
  ObitTableGCRow*       outRow=NULL;
  olong i, j, iver, iRow, oRow, ver;
  oint numIF, numPol, numTabs;
  ObitIOAccess access;
  gchar *routine = "GetGainCurveInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Loop over plausible versions */
  for (iver=1; iver<1000; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input Flag table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    numIF = numPol = numTabs = 0;
    inTable = newObitTableIDI_GAIN_CURVEValue ("Input table", inData, &ver, access, 
					       numPol, numIF, numTabs, 
					       err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      break;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableIDI_GAIN_CURVEOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_GAIN_CURVE table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableIDI_GAIN_CURVERow (inTable);
    
    /* Create output GC table object */
    ver = iver;
    access  = OBIT_IO_ReadWrite;
    numIF   = inTable->no_band;
    numPol  = inTable->numPol;
    numTabs =  inTable->numTabs;
    outTable = newObitTableGCValue ("Output table", (ObitData*)outData, 
				    &ver, access, numIF, numPol, numTabs, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with GC table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableGCOpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output GC table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTableGCRow (outTable);
    /* attach to table buffer */
    ObitTableGCSetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Initialize output row */
    outRow->status = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableIDI_GAIN_CURVEReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading IDI_GAIN_CURVE Table");
	return;
      }
      
      /* Save to GC table */
      outRow->antennaNo  = inRow->antNo;
      outRow->SubArray   = inRow->Array;
      outRow->FreqID     = inRow->fqid;
      for (i=0; i<numIF; i++) {
	outRow->Type1[i] = inRow->type1[i];
	outRow->NTerm1[i]= inRow->nterm1[i];
	outRow->XTyp1[i] = inRow->x_typ1[i];
	outRow->YTyp1[i] = inRow->y_typ1[i];
	outRow->XVal1[i] = inRow->x_val1[i];
	outRow->sens1[i] = inRow->sens1[i];
	for (j=0; j<numTabs; j++) {
	  outRow->YVal1[i*numTabs+j] = inRow->y_val1[i*numTabs+j];
	  outRow->gain1[i*numTabs+j] = inRow->gain1[i*numTabs+j];
	}
      }
      if (numPol>1) {   /* 2 poln */
	for (i=0; i<numIF; i++) {
	  outRow->Type2[i] = inRow->type2[i];
	  outRow->NTerm2[i]= inRow->nterm2[i];
	  outRow->XTyp2[i] = inRow->x_typ2[i];
	  outRow->YTyp2[i] = inRow->y_typ2[i];
	  outRow->XVal2[i] = inRow->x_val2[i];
	  outRow->sens2[i] = inRow->sens2[i];
	  for (j=0; j<numTabs; j++) {
	    outRow->YVal2[i*numTabs+j] = inRow->y_val2[i*numTabs+j];
	    outRow->gain2[i*numTabs+j]  = inRow->gain2[i*numTabs+j];
	  }
	}
      }

      /* Write */
      oRow = -1;
      if ((ObitTableGCWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating GC Table");
	return;
      }

    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_GAIN_CURVEClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_GAIN_CURVE Table file");
      return;
    }
    
    if ((ObitTableGCClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output GC Table file");
      return;
    }

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied GAIN_CURVE table %d", iver);

    /* Cleanup */
    inRow    = ObitTableIDI_GAIN_CURVERowUnref(inRow);
    inTable  = ObitTableIDI_GAIN_CURVEUnref(inTable);
    outRow   = ObitTableGCRowUnref(outRow);
    outTable = ObitTableGCUnref(outTable);

  } /* end loop over versions */

} /* end  GetGainCurveInfo */

void GetPhaseCalInfo (ObitData *inData, ObitUV *outData, 
				 ObitErr *err)
/*----------------------------------------------------------------------- */
/* Convert any PHASE_CAL tables on inData to AIPS PC on outData           */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_PHASE_CAL*    inTable=NULL;
  ObitTableIDI_PHASE_CALRow* inRow=NULL;
  ObitTablePC*          outTable=NULL;
  ObitTablePCRow*       outRow=NULL;
  olong i, iver, iRow, oRow, ver;
  oint numIF, numPol, numTones;
  ObitIOAccess access;
  gchar *routine = "GetPhaseCalInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Loop over plausible versions */
  for (iver=1; iver<1000; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input Flag table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    numIF = numPol = numTones = 0;
    inTable = newObitTableIDI_PHASE_CALValue ("Input table", inData, &ver, access, 
					      numPol, numIF, numTones, 
					      err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      break;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableIDI_PHASE_CALOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_PHASE_CAL table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableIDI_PHASE_CALRow (inTable);
    
    /* Create output PC table object */
    ver = iver;
    access   = OBIT_IO_ReadWrite;
    numIF    = inTable->no_band;
    numPol   = inTable->numPol;
    numTones =  inTable->numTones;
    outTable = newObitTablePCValue ("Output table", (ObitData*)outData, 
				    &ver, access, numPol, numIF, numTones, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with PC table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTablePCOpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output PC table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTablePCRow (outTable);
    /* attach to table buffer */
    ObitTablePCSetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Initialize output row */
    outRow->status      = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableIDI_PHASE_CALReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading IDI_PHASE_CAL Table");
	return;
      }
      
      /* Save to PC table */
      outRow->Time     = inRow->Time;
      outRow->SourID   = inRow->SourID;
      outRow->antennaNo= inRow->antennaNo;
      outRow->Array    = inRow->Array;
      outRow->FreqID   = inRow->FreqID;
      outRow->CableCal = inRow->CableCal;
      for (i=0; i<numTones; i++) {
	outRow->State1[i]   = inRow->State1[i];
	outRow->PCFreq1[i]  = inRow->PCFreq1[i];
	outRow->PCReal1[i]  = inRow->PCReal1[i];
	outRow->PCImag1[i]  = inRow->PCImag1[i];
	outRow->PCRate1[i]  = inRow->PCRate1[i];
	
	if (numPol>1) {   /* 2 poln */
	  outRow->PCFreq2[i]  = inRow->PCFreq2[i];
	  outRow->PCReal2[i]  = inRow->PCReal2[i];
	  outRow->PCImag2[i]  = inRow->PCImag2[i];
	  outRow->PCRate2[i]  = inRow->PCRate2[i];
	  outRow->State2[i]   = inRow->State2[i];
	}
      }

      /* Write */
      oRow = -1;
      if ((ObitTablePCWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating PC Table");
	return;
      }

    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_PHASE_CALClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_PHASE_CAL Table file");
      return;
    }
    
    if ((ObitTablePCClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output PC Table file");
      return;
    }

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied PHASE_CAL table %d", iver);

    /* Cleanup */
    inRow    = ObitTableIDI_PHASE_CALRowUnref(inRow);
    inTable  = ObitTableIDI_PHASE_CALUnref(inTable);
    outRow   = ObitTablePCRowUnref(outRow);
    outTable = ObitTablePCUnref(outTable);

  } /* end loop over versions */

} /* end  GetPhaseCalInfo */

void GetInterferometerModelInfo (ObitData *inData, ObitUV *outData, 
				 ObitErr *err)
/*----------------------------------------------------------------------- */
/* Convert any INTERFEROMETER_MODEL tables on inData to AIPS IM on outData*/
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_INTERFEROMETER_MODEL*    inTable=NULL;
  ObitTableIDI_INTERFEROMETER_MODELRow* inRow=NULL;
  ObitTableIM*          outTable=NULL;
  ObitTableIMRow*       outRow=NULL;
  olong i, j, iver, iRow, oRow, ver;
  oint numIF, numPol, npoly;
  ObitIOAccess access;
  gchar *routine = "GetInterferometerModelInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Loop over plausible versions */
  for (iver=1; iver<1000; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input Flag table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    numIF = numPol = npoly = 0;
    inTable = newObitTableIDI_INTERFEROMETER_MODELValue ("Input table", inData, &ver, access, 
							 numPol, npoly, numIF, 
							 err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      break;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableIDI_INTERFEROMETER_MODELOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_INTERFEROMETER_MODEL table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableIDI_INTERFEROMETER_MODELRow (inTable);
    
    /* Create output IM table object */
    ver = iver;
    access  = OBIT_IO_ReadWrite;
    numIF   = inTable->no_band;
    numPol  = inTable->numPol;
    npoly   =  inTable->npoly;
    outTable= newObitTableIMValue ("Output table", (ObitData*)outData, 
				   &ver, access, numPol, numIF, npoly, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IM table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableIMOpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output IM table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTableIMRow (outTable);
    /* attach to table buffer */
    ObitTableIMSetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Initialize output row */
    outRow->status      = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableIDI_INTERFEROMETER_MODELReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading IDI_INTERFEROMETER_MODEL Table");
	return;
      }
      
      /* Save to IM table */
      outRow->Time      = inRow->Time;
      outRow->TimeI     = inRow->TimeI;
      outRow->SourID    = inRow->SourID;
      outRow->antennaNo = inRow->antennaNo;
      outRow->Array     = inRow->Array;
      outRow->FreqID    = inRow->FreqID;
      outRow->IFR       = inRow->IFR;
      outRow->Disp1     = inRow->Disp1;
      outRow->DRate1    = inRow->DRate1;
    for (i=0; i<numIF; i++) {
	outRow->FreqVar[i]   = inRow->FreqVar[i];
	for (j=0; j<npoly; j++) {
	  outRow->PDelay1[i*npoly+j]   = inRow->PDelay1[i*npoly+j];
	  outRow->PRate1[i*npoly+j]    = inRow->PRate1[i*npoly+j];
	  outRow->Disp1     = inRow->Disp1;
	}
      }
      /* There are more in the IDI table */
      for (j=0; j<npoly; j++) {
	outRow->GDelay1[j]   = inRow->GDelay1[j];
	outRow->GRate1[j]    = inRow->GRate1[j];
      }
      if (numPol>1) {   /* 2 poln */
	outRow->Disp2  = inRow->Disp2;
	outRow->DRate2 = inRow->DRate2;
	for (i=0; i<numIF; i++) {
	  for (j=0; j<npoly; j++) {
	    outRow->PDelay2[i*npoly+j]   = inRow->PDelay2[i*npoly+j];
	    outRow->PRate2[i*npoly+j]    = inRow->PRate2[i*npoly+j];
	  }
	}
	for (j=0; j<npoly; j++) {
	  /* There are more in the IDI table */
	  outRow->GDelay2[j]   = inRow->GDelay2[j];
	  outRow->GRate2[j]    = inRow->GRate2[j];
	}
      } /* end poln 2 */

      /* Write */
      oRow = -1;
      if ((ObitTableIMWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating IM Table");
	return;
      }

    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_INTERFEROMETER_MODELClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_INTERFEROMETER_MODEL Table file");
      return;
    }
    
    if ((ObitTableIMClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output IM Table file");
      return;
    }

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied INTERFEROMETER_MODEL table %d", iver);

    /* Cleanup */
    inRow    = ObitTableIDI_INTERFEROMETER_MODELRowUnref(inRow);
    inTable  = ObitTableIDI_INTERFEROMETER_MODELUnref(inTable);
    outRow   = ObitTableIMRowUnref(outRow);
    outTable = ObitTableIMUnref(outTable);

  } /* end loop over versions */

} /* end  GetInterferometerModelInfo */

void GetWeatherInfo (ObitData *inData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any WEATHER tables on inData to AIPS WX on outData            */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_WEATHER*    inTable=NULL;
  ObitTableIDI_WEATHERRow* inRow=NULL;
  ObitTableWX*          outTable=NULL;
  ObitTableWXRow*       outRow=NULL;
  olong iver, iRow, oRow, ver;
  oint numIF;
  ObitIOAccess access;
  gchar *routine = "GetWeatherInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Loop over plausible versions */
  for (iver=1; iver<1000; iver++) {
    
    /* Print any messages */
    ObitErrLog(err);

    /* Create input Flag table object */
    ver = iver;
    access = OBIT_IO_ReadOnly;
    numIF =0;
    inTable = newObitTableIDI_WEATHERValue ("Input table", inData, 
					    &ver, access, numIF,err);
    /* Find it? */
    if (inTable==NULL) {
      ObitErrClearErr (err);
      break;
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Open table */
    if ((ObitTableIDI_WEATHEROpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_WEATHER table");
      return;
    }
    
    /* Create Row */
    inRow = newObitTableIDI_WEATHERRow (inTable);
    
    /* Create output WX table object */
    ver = iver;
    access  = OBIT_IO_ReadWrite;
    outTable = newObitTableWXValue ("Output table", (ObitData*)outData, 
				    &ver, access, err);
    if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with WX table");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Open table */
    if ((ObitTableWXOpen (outTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output WX table");
      return;
    }
    
    /* Create output Row */
    outRow = newObitTableWXRow (outTable);
    /* attach to table buffer */
    ObitTableWXSetRow (outTable, outRow, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Initialize output row */
    outRow->status      = 0;
    
    /* loop through input table */
    for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
      if ((ObitTableIDI_WEATHERReadRow (inTable, iRow, inRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading IDI_WEATHER Table");
	return;
      }
      
      /* Save to WX table */
      outRow->Time          = inRow->Time;
      outRow->TimeI         = inRow->TimeI;
      outRow->antNo         = inRow->antNo;
      outRow->temperature   = inRow->temperature;
      outRow->pressure      = inRow->pressure;
      outRow->dewpoint      = inRow->dewpoint;
      outRow->windVelocity  = inRow->wind_velocity;
      outRow->windDirection = inRow->wind_direction;
      outRow->wvrH2O        = inRow->wvr_h2o;
      outRow->onosElectron  = inRow->ionos_electron;
      /* Write */
      oRow = -1;
      if ((ObitTableWXWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating WX Table");
	return;
      }

    } /* end loop over input table */
    
    /* Close  tables */
    if ((ObitTableIDI_WEATHERClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_WEATHER Table file");
      return;
    }
    
    if ((ObitTableWXClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output WX Table file");
      return;
    }

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied WEATHER table %d", iver);

    /* Cleanup */
    inRow    = ObitTableIDI_WEATHERRowUnref(inRow);
    inTable  = ObitTableIDI_WEATHERUnref(inTable);
    outRow   = ObitTableWXRowUnref(outRow);
    outTable = ObitTableWXUnref(outTable);
  } /* end loop over versions */

} /* end  GetWeatherInfo */

/*----------------------------------------------------------------------- */
/*  Write History for IDIIn                                               */
/*   Input:                                                               */
/*      inData    FITS IDI to copy history from                           */
/*      myInput   Input parameters on InfoList                            */
/*      outData   ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void IDIInHistory (ObitData* inData, ObitInfoList* myInput, ObitUV* outData, 
		   ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "Scan", "DataRoot",
    NULL};
  gchar *routine = "IDIInHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(outData));

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

  outHistory = ObitHistoryUnref(outHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end IDIInHistory  */
