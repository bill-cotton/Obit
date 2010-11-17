/* $Id:  $  */
/* Read BDF format data, convert to Obit UV                           */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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
#include "ObitUVUtil.h"
#include "ObitTableSU.h"
#include "ObitTableSUUtil.h"
#include "ObitTableAN.h"
#include "ObitTableANUtil.h"
#include "ObitTableFQ.h"
#include "ObitTableFG.h"
#include "ObitTableCD.h"
#include "ObitTableCL.h"
#include "ObitTableCLUtil.h"
#include "ObitTableBP.h"
#include "ObitTableTY.h"
#include "ObitTableIM.h"
#include "ObitTableGC.h"
#include "ObitTablePC.h"
#include "ObitTableSY.h"
#include "ObitTableWX.h"
#include "ObitTableNX.h"
#include "ObitTableOT.h"
#include "ObitSDMData.h"
#include "ObitBDFData.h"
#include "ObitHistory.h"
#include "ObitPrecess.h"
#include "ObitVLAGain.h"
#include "ObitGainCal.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif

/* internal prototypes */
/* Get inputs */
ObitInfoList* BDFInin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Create output uvdata */
ObitUV* setOutputData (ObitInfoList *myInput, ObitErr *err);
/* Get file descriptor */
void GetHeader (ObitUV *outData, ObitSDMData *SDMData, ObitInfoList *myInput, 
		ObitErr *err);
/* Get data */
void GetData (ObitSDMData *SDMData, ObitInfoList *myInput, ObitUV *outData, 
	      ObitErr *err);

/* Get Antenna info */
void GetAntennaInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err);
/* Get Frequency info */
void GetFrequencyInfo (ObitSDMData *SDMData, ObitUV *outData, 
		       ASDMSpectralWindowArray* SpWinArray, ObitErr *err);
/* Get Source info */
void GetSourceInfo (ObitSDMData *SDMData, ObitUV *outData, olong iScan, 
		    ObitErr *err);
/* Copy any FLAG tables */
void GetFlagInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err);
/* Copy  CalDevice to AIPS CD table */
void GetCalDeviceInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err);
/* Copy  SysPower to AIPS SY table */
void GetSysPowerInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err);
/* Copy  Over the top from Pointing table to AIPS OT table */
void GetOTTInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err);
/* Copy any CALIBRATION tables */
void GetCalibrationInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Copy any BANDPASS tables */
void GetBandpassInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Copy any SYSTEM_TEMPERATURE tables */
void GetTSysInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Get Gain curve  (GC) table */
void GetGainCurveInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err);
/* Copy any PHASE_CALtables */
void GetPhaseCalInfo (ObitData *inData, ObitUV *outData, 
				 ObitErr *err);
/* Copy any INTERFEROMETER_MODELtables */
void GetInterferometerModelInfo (ObitData *inData, ObitUV *outData, 
				 ObitErr *err);
/* Copy any WEATHER tables */
void GetWeatherInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err);
/* Write history */
void BDFInHistory (ObitInfoList* myInput, ObitSDMData *SDMData, ObitUV* outData, 
		   ObitErr* err);
/* Update AN tables */
void UpdateAntennaInfo (ObitUV *outData, olong arrno, ObitErr *err);
/* Calculate UVW */
void CalcUVW (ObitUV *outData, ObitBDFData *BDFData, ofloat *Buffer, ObitErr *err);
/* Days to human string */
void day2dhms(ofloat time, gchar *timeString);
/* Check for zero visibilities */
gboolean CheckAllZero(ObitUVDesc *desc, ofloat *Buffer);
/*  Write FG table entries for scans intended for online calibration */
void FlagIntent (ObitSDMData *SDMData, ObitUV* outData, ObitErr* err);
/* Is a given time in next scan? */
gboolean nextScan(ObitSDMData *SDMData, olong curScan, odouble time,  
		  olong *curScanI, olong *nextScanNo, olong *SourID);
/* Nominal VLA sensitivity */
ofloat nomSen(ASDMAntennaArray*  AntArray);

/* Program globals */
gchar *pgmName = "BDFIn";       /* Program name */
gchar *infile  = "BDFIn.inp";   /* File with program inputs */
gchar *outfile = "BDFIn.out";   /* File to contain program outputs */
olong  pgmNumber;      /* Program number (like POPS no.) */
olong  AIPSuser;       /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;        /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;        /* Number of FITS directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
gchar **FITSdirs=NULL; /* List of FITS data directories */
odouble refJD = 0.0;   /* reference Julian date */
odouble refMJD = 0.0;  /* reference Modified Julian date */
odouble integTime;     /* Integration time in days */
olong  nchan=1;        /* Number of frequencies */
olong  nstok=1;        /* Number of Stokes */
olong  nIF=1;          /* Number of IFs */
ofloat deltaFreq;      /* Channel increment */
ofloat refPixel;       /* Frequency reference pixel */
odouble refFrequency;  /* reference frequency (Hz) */
olong numArray;        /* Number of input arrays */
odouble *dataRefJDs =NULL; /* Array of reference JD from data per input array */
odouble *arrayRefJDs=NULL; /* Array of reference JD from array geometry per input array */
odouble *arrRefJDCor=NULL; /* Array of input array ref day corrections */
ObitAntennaList **antennaLists=NULL;  /* Array of antenna lists for uvw calc */
ObitSourceList *uvwSourceList=NULL;   /* Source List for uvw calc */
olong uvwSourID=-1;                   /* Source ID for uvw calc */
olong uvwcurSourID=-1;                /* Current source ID for uvw calc */
ofloat uvrot=-0.0;                    /* Current source rotation of u-v */
ofloat **antLongs=NULL;               /* Array of Antenna longitudes per subarray */
ofloat **antLats=NULL;                /* Array of Antenna latitudes per subarray */
olong selScan=-1;                     /* First selected scan number */
olong selChan=-1;                     /* Selected number of channels */
olong selIF=-1;                       /* Selected number of IFs (SpWin) */
gboolean isEVLA;                      /* Is this EVLA data? */
gboolean SWOrder=FALSE;               /* Leave in SW Order? */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*    Read BDF  data to an Obit UV dataset                                */
/*----------------------------------------------------------------------- */
{
  olong  i, ierr=0;
  ObitSystem *mySystem= NULL;
  ObitUV *outData= NULL;
  ObitErr *err= NULL;
  gboolean doOnline=FALSE;
  gchar dataroot[132];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitSDMData *SDMData=NULL;
  ObitBDFData *BDFData=NULL;
  ObitTable   *CalTab=NULL;

  err = newObitErr();  /* Obit error/message stack */

  /* Startup - parse command line */
  ierr = 0;
  myInput = BDFInin (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input DataRoot name */
  for (i=0; i<132; i++) dataroot[i] = 0;
  ObitInfoListGet(myInput, "DataRoot", &type, dim, dataroot, err);
  dataroot[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(dataroot);  /* Trim trailing blanks */

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);

  if (err->error) ierr = 1;  ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Swallow SDM */
  SDMData = ObitSDMDataCreate ("SDM", dataroot, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Create ObitUV for data */
  outData = setOutputData (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Get header info, array geometry, initialize output */
  GetHeader (outData, SDMData, myInput, err); 
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit; 

  /* Open output data */
  if ((ObitUVOpen (outData, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output UV file %s", outData->name);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
  
  /* convert data  */
  GetData (SDMData, myInput, outData, err); 
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Close output uv data */
  if ((ObitUVClose (outData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing output file");
  
  /* Copy EVLA tables */
  if (isEVLA) {
    GetFlagInfo (SDMData, outData, err);       /*   FLAG tables */
    /* GetCalibrationInfo (inData, outData, err);   CALIBRATION tables */
    /* GetBandpassInfo (inData, outData, err);      BANDPASS tables */
    /* GetTSysInfo (inData, outData, err);          SYSTEM_TEMPERATURE tables */
    /* GetInterferometerModelInfo (inData, outData, err); INTERFEROMETER_MODEL tables */
    /* GetPhaseCalInfo (inData, outData, err);      PHASE_CAL tables */
    GetCalDeviceInfo (SDMData, outData, err);  /*   calDevice table */
    GetSysPowerInfo  (SDMData, outData, err);  /*   SysPower table */
    GetOTTInfo  (SDMData, outData, err);       /*   Over the top table */
    GetGainCurveInfo (SDMData, outData, err);  /*   gain curve (GC) table */
  }
  GetWeatherInfo   (SDMData, outData, err);  /*   Weather table */
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Check scan intents for online only calibrations */
  ObitInfoListGetTest(myInput, "doOnline", &type, dim, &doOnline);
  if (isEVLA && !doOnline) FlagIntent (SDMData, outData, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Update An tables with correct ref. date */
  /*for (i=1; i<=numArray; i++)  UpdateAntennaInfo (outData, i, err);*/
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Create CL table - interval etc,  set in setOutputData */
  if (isEVLA) {
    Obit_log_error(err, OBIT_InfoErr, "Creating CL Table");
    ObitErrLog(err);
    /*ObitTableCLGetDummy (outData, outData, 1, err); */
    CalTab = ObitGainCalCalc (outData, FALSE, err);
    if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
  }

  /* History */
  BDFInHistory (myInput, SDMData, outData, err);
  
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;
  
  /* Shutdown Obit */
 exit:
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  mySystem = ObitSystemShutdown (mySystem);
  
  /* cleanup */
  myInput  = ObitInfoListUnref(myInput);   /* delete input list */
  myOutput = ObitInfoListUnref(myOutput);  /* delete output list */
  SDMData  = ObitSDMDataUnref (SDMData);
  BDFData  = ObitBDFDataUnref (BDFData);
  outData  = ObitUnref(outData);
  CalTab   = ObitTableUnref(CalTab);
  uvwSourceList = ObitUnref(uvwSourceList);
  if (dataRefJDs)  g_free(dataRefJDs);
  if (arrayRefJDs) g_free(arrayRefJDs);
  if (arrRefJDCor) g_free(arrRefJDCor);
  if (antennaLists) {
    for (i=0; i<numArray; i++) antennaLists[i] = ObitUnref(antennaLists[i]);
    g_free(antennaLists);
  }
 
  return ierr;
} /* end of main */

ObitInfoList* BDFInin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="BDFIn.in", *arg;
  gboolean init=FALSE;
  oint itemp;
  gchar *strTemp;
  ObitInfoList* list;
  gchar *routine = "BDFInin";

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
      
    } else if (strcmp(arg, "-File") == 0){ /* Scan name */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "File", OBIT_string, dim, strTemp);

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
} /* end BDFInin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: BDFIn -input file -output ofile [args]\n");
    fprintf(stderr, "Convert an BDF/BDF file format to Obit/UV\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def BDFIn.in\n");
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
  ObitInfoListPut (out, "File", OBIT_string, dim, scan_name, err);

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
  olong      nvis, itemp;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean  exist, btemp;
  gchar     tname[129], *fullname=NULL;
  ofloat    calInt = 0.5, ftemp;
  gchar     *outParms[] = {  /* Parameters for output data */
    "Compress", "maxGap", "maxScan", "doSwPwr", NULL};
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
    /* Default out class is "BDFIn" */
    if (!strncmp(Aclass, "      ", 6)) strncpy (Aclass, "BDFIn", 7);

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

   
    /* define object */
    nvis = 1;
    ObitUVSetAIPS (outUV, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output AIPS UV data %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);
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

    /* Did it previously exist */
    fullname = ObitFITSFilename (disk, outFile, err);
    exist =  ObitFileExist (fullname, err);
    if (fullname) g_free(fullname);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
    /* define object */
    nvis = 1;
    ObitUVSetFITS (outUV, nvis, disk, outFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output FITS UV data %s on disk %d", outFile, disk);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outUV;
  }

  /* Set CL table interval */
  ObitInfoListGetTest(myInput, "calInt", &type, dim, &calInt);
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  /* to seconds */
  calInt *= 60.0;
  ObitInfoListAlwaysPut(outUV->info, "calInt", OBIT_float, dim, &calInt);
  
  /* Other CL Table controls */
  itemp = 1;   /* Output CL table version */
  ObitInfoListAlwaysPut(outUV->info, "calVer", OBIT_float, dim, &itemp);
  btemp = TRUE;
  ObitInfoListAlwaysPut(outUV->info, "doGain",  OBIT_bool, dim, &btemp);
  ObitInfoListAlwaysPut(outUV->info, "doOpac",  OBIT_bool, dim, &btemp);
  ftemp = 0.5;   /* Weight of weather in opacity */
  ObitInfoListAlwaysPut(outUV->info, "WXWeight", OBIT_float, dim, &ftemp);
  ObitInfoListAlwaysPut(outUV->info, "solInt", OBIT_float, dim, &calInt);
  
  /* Get input parameters from myInput, copy to outUV */
  ObitInfoListCopyList (myInput, outUV->info, outParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

  return outUV;
} /* end setOutputUV */

void GetHeader (ObitUV *outData, ObitSDMData *SDMData, ObitInfoList *myInput, 
		ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get header information from scan header files                         */
/*  Returns TRUE if the file is just created                              */
/*   Input:                                                               */
/*      outData  Output UV object                                         */
/*      SDMData  ASDM structure                                           */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitUVDesc *desc;
  olong ncol;
  gchar *today=NULL;
  olong i, lim, selConfig, iScan;
  ofloat epoch=2000.0, equinox=2000.0;
  olong nchan=1, npoln=1, nIF=1;
  odouble refFreq, startFreq=1.0;
  ofloat freqStep=1.0;
  ASDMSpectralWindowArray* SpWinArray;
  ASDMAntennaArray*  AntArray;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar selBand[12];
  ObitASDMBand band;
  gchar *bandCodes[] = {"Any", "4","P","L","S","C","X","Ku","K","Ka","Q","W"};
  gchar *routine = "GetHeader";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outData));
  g_assert(SDMData!=NULL);
  g_assert(myInput!=NULL);

  /* Get total number of arrays - only one */
  numArray = 1;
  if (numArray>=1) {
    /* Create array of reference JDs */
    dataRefJDs  = g_malloc(numArray*sizeof(odouble));
    arrayRefJDs = g_malloc(numArray*sizeof(odouble));
    arrRefJDCor = g_malloc(numArray*sizeof(odouble));
    for (i=0; i<numArray; i++) {
      dataRefJDs[i]  = -1.0;
      arrayRefJDs[i] = -1.0;
      arrRefJDCor[i] =  0.0;
    }
  }

  /* ConfigID selection */
  selConfig = -1;
  ObitInfoListGetTest(myInput, "selConfig", &type, dim, &selConfig);
  
  /* Channel selection */
  selChan = 0;
  ObitInfoListGetTest(myInput, "selChan", &type, dim, &selChan);

  /* Band selection */
  for (i=0; i<12; i++) selBand[i] = 0;
  ObitInfoListGetTest(myInput, "selBand", &type, dim, selBand);
  band = ObitSDMDataBand2Band (selBand);

  /* IF selection */
  selIF = 0;
  ObitInfoListGetTest(myInput, "selIF", &type, dim, &selIF);

  /* Check if Spectral window order desired */
  ObitInfoListGetTest(myInput, "SWOrder", &type, dim, &SWOrder);
  if (SWOrder)
    Obit_log_error(err, OBIT_InfoWarn, 
		   "Frequencies in Spectral Window Possibly NOT Frequency order");

  /* Find first selected scan */
  iScan = ObitASDSelScan (SDMData, selChan, selIF, band, selConfig);
  Obit_return_if_fail((iScan>=0), err,
		      "%s: No scans found matching selection criteria", 
		      routine);
  /* Extract ASDM Spectral windows data  */
  selScan = iScan;    /* Save to global */
  SpWinArray  = ObitSDMDataGetSWArray (SDMData, iScan, SWOrder);
  Obit_return_if_fail((SpWinArray), err,
		      "%s: Could not extract Spectral Windows from ASDM", 
		      routine);

  /* selConfig overrides selBand, selIF */
  if (selConfig>=0) {
    /* Set selIF, selBand, default selChan */
    selIF = SpWinArray->nwinds;
    band  = ObitSDMDataFreq2Band(SpWinArray->refFreq);
    /* Default selChans? use first Spectral window */
    if (selChan<=0) selChan = SpWinArray->winds[0]->numChan;

  } else { /* Selection including selBand, selIF */
    /* Default selChans? use first Spectral window */
    if (selChan<=0) selChan = SpWinArray->winds[0]->numChan;
    if (selIF<=0)   selIF   = SpWinArray->nwinds;
    if (band==ASDMBand_Any)  band  = ObitSDMDataFreq2Band(SpWinArray->refFreq);
    /* Set default configID */
    if (selConfig<0) selConfig = SDMData->MainTab->rows[iScan]->configDescriptionId;
  }

  /* Save selection for history */
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(myInput, "selConfig", OBIT_long, dim, &selConfig);
  Obit_log_error(err, OBIT_InfoErr, "Selecting scans with ConfigID %d", selConfig);
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  strcpy(selBand,  bandCodes[band]);
  dim[0] = strlen(selBand);
  ObitInfoListAlwaysPut(myInput, "selBand", OBIT_string, dim, selBand);
  Obit_log_error(err, OBIT_InfoErr, "Selecting %s band", bandCodes[band]);
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(myInput, "selIF", OBIT_long, dim, &selIF);
  Obit_log_error(err, OBIT_InfoErr, "Selecting scans with %d Spectral Windows", selIF);
  ObitInfoListAlwaysPut(myInput, "selChan", OBIT_long, dim, &selChan);
  Obit_log_error(err, OBIT_InfoErr, "Selecting spectral windows with %d channels", selChan);
  ObitErrLog(err);
 
 /* Define output descriptor */
  desc = outData->myDesc;
  /* Extract ASDM data  */
  SpWinArray  = ObitSDMDataGetSWArray (SDMData, iScan, SWOrder);
  Obit_return_if_fail((SpWinArray), err,
		      "%s: Could not extract Spectral Windows from ASDM", 
		      routine);
  AntArray    = ObitSDMDataGetAntArray(SDMData, iScan);
  Obit_return_if_fail((AntArray), err,
		      "%s: Could not extract Antenna info from ASDM", 
		      routine);
  /* Select Spectral windows */
  ObitSDMDataSelChan (SpWinArray, selChan, selIF, ASDMBand_Any);

  /* Frequency info from first selected Spectral window */
  for (i=0; i<SpWinArray->nwinds; i++) {
    if (SpWinArray->winds[i]->selected) {
      nchan     = SpWinArray->winds[i]->numChan;
      npoln     = SpWinArray->winds[i]->nCPoln;
      refFreq   = SpWinArray->winds[i]->refFreq;
      startFreq = SpWinArray->winds[i]->chanFreqStart;
      freqStep  = (ofloat)SpWinArray->winds[i]->chanFreqStep;
      break;
    }
  }

  /* EVLA always effectively USB */
  if (!strncmp(AntArray->arrayName, "EVLA", 4)) freqStep = fabs(freqStep );
  
  /* Count selected Spectral windows */
  nIF = 0;
  for (i=0; i<SpWinArray->nwinds; i++) {
    if (SpWinArray->winds[i]->selected) nIF++;
  }
  Obit_log_error(err, OBIT_InfoErr, "Selecting %d spectral windows (IFs)", nIF);
  ObitErrLog(err);
  
  /* Define header */
  desc->nvis = 0;
  strncpy (desc->origin, "Obit ", UVLEN_VALUE);
  desc->isort[0] = 'T'; 
  desc->isort[1] = 'B'; 
  desc->nrparm = 0;
  
  /* Creation date today */
  today = ObitToday();
  strncpy (desc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  desc->freq    = startFreq;
  desc->JDObs   = SDMData->refJD;
  refJD         = SDMData->refJD;
  desc->epoch   = epoch;
  desc->equinox = equinox;
  desc->obsra   = 0.0;
  desc->obsdec  = 0.0;
  desc->altCrpix = 0.0;
  desc->altRef = 0.0;
  desc->restFreq = 0.0;
  desc->xshift = 0.0;
  desc->yshift = 0.0;
  desc->VelDef = 0;
  desc->VelReference = 0;
  strncpy (desc->bunit, "        ", UVLEN_VALUE);
  lim = UVLEN_VALUE;
  ObitUVDescJD2Date (SDMData->refJD, desc->obsdat);
  strncpy (desc->teles,      AntArray->arrayName, lim);
  strncpy (desc->observer,   AntArray->obsName, lim);
  strncpy (desc->instrument, AntArray->arrayName, lim);
  strncpy (desc->object,   "MULTI   ", lim);
  
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
  
  /* Source  */
  strncpy (desc->ptype[ncol], "SOURCE  ", UVLEN_KEYWORD);
  ncol++;
  
  /* Integration time */
  strncpy (desc->ptype[ncol], "INTTIM  ", UVLEN_KEYWORD);
  ncol++;
  
  /* FreqID - NYI */
  /* strncpy (desc->ptype[ncol], "FREQSEL ", UVLEN_KEYWORD);
     ncol++; */
  
  desc->nrparm = ncol;  /* Number of random parameters */
  
  /* Data Matrix */
  ncol = 0;
  lim = UVLEN_KEYWORD;
  /* Dimension 1 = COMPLEX */
  strncpy (desc->ctype[ncol], "COMPLEX ", lim);
  desc->inaxes[ncol] = 3;
  desc->cdelt[ncol]  = 1.0;
  desc->crpix[ncol]  = 1.0;
  desc->crval[ncol]  = 1.0;
  desc->crota[ncol]  = 0.0;
  ncol++;
  
  /* Dimension 2 = STOKES */
  strncpy (desc->ctype[ncol], "STOKES  ", lim);
  desc->inaxes[ncol] = npoln;
  desc->crota[ncol]  =  0.0;
  desc->cdelt[ncol]  = -1.0;
  desc->crpix[ncol]  =  1.0;
  /* Circular or linear */
  if (AntArray->ants[0]->polnType[0]=='R') { /* Circular */
    desc->crval[ncol]  = -1.0;
  } else { /* Assume linear */
    desc->crval[ncol]  = -5.0;
  }
  ncol++;
 
  /* Dimension 3 = FREQ */
  strncpy (desc->ctype[ncol], "FREQ    ", lim);
  desc->inaxes[ncol] = nchan;
  desc->cdelt[ncol]  = freqStep;
  desc->crpix[ncol]  = 1.0;
  desc->crval[ncol]  = startFreq;
  desc->crota[ncol]  = 0.0;
  ncol++;
  
  /* Dimension 4 = IF (SpectralWindow in BDF) */
  strncpy (desc->ctype[ncol], "IF      ", lim);
  desc->inaxes[ncol] = nIF;
  desc->cdelt[ncol]  = 1.0;
  desc->crpix[ncol]  = 1.0;
  desc->crval[ncol]  = 1.0;
  desc->crota[ncol]  = 0.0;
  ncol++;  
  
  /* Dimension 5 Dummy RA */
  strncpy (desc->ctype[ncol], "RA      ", lim);
  desc->inaxes[ncol] = 1;
  desc->cdelt[ncol]  = 1.0;
  desc->crpix[ncol]  = 1.0;
  desc->crval[ncol]  = 0.0;
  desc->crota[ncol]  = 0.0;
  ncol++;  
  
  /* Dimension 6 Dummy Dec */
  strncpy (desc->ctype[ncol], "DEC     ", lim);
  desc->inaxes[ncol] = 1;
  desc->cdelt[ncol]  = 1.0;
  desc->crpix[ncol]  = 1.0;
  desc->crval[ncol]  = 0.0;
  desc->crota[ncol]  = 0.0;
  ncol++;  
  
  desc->naxis = ncol;  /* Number of dimensions */
  
  /* index descriptor */
  ObitUVDescIndex (desc);
  
  /* Get source info, copy to output SU table, save lookup table in SourceID */
  ObitUVOpen (outData, OBIT_IO_WriteOnly, err);
  GetSourceInfo (SDMData, outData, iScan, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Save any equinox information */
  outData->myDesc->epoch   = 2000;   /* ASDM Lacking */
  outData->myDesc->equinox = 2000;
  epoch   = outData->myDesc->epoch;
  equinox = outData->myDesc->equinox;
 
  /* Add Antenna and Frequency info */
  GetAntennaInfo (SDMData, outData,err);
  GetFrequencyInfo (SDMData, outData,  SpWinArray, err);
  ObitUVClose (outData, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Instantiate output Data */
  ObitUVFullInstantiate (outData, FALSE, err);
  if (err->error)Obit_traceback_msg (err, routine, outData->name);

  /* Cleanup */
  SpWinArray = ObitSDMDataKillSWArray (SpWinArray);
  AntArray   = ObitSDMDataKillAntArray(AntArray);

  return;
} /* end GetHeader */

/*----------------------------------------------------------------------- */
/*  Write History for BDFIn                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      SDMData   ASDM structure                                          */
/*      outData   ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void BDFInHistory (ObitInfoList* myInput, ObitSDMData *SDMData, 
		   ObitUV* outData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  ASDMScanTable* ScanTab;
  olong          iScan, iIntent;
  gchar          hicard[81], begString[17], endString[17];
  gchar         *hiEntries[] = {
    "DataRoot", "selChan", "selIF", "selBand", "selConfig", "dropZero", 
    "doCode", "calInt", "doSwPwr", "doOnline", "SWOrder", 
    NULL};
  gchar *routine = "BDFInHistory";
  
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(outData));

  /* Do history  */
  outHistory = newObitDataHistory ((ObitData*)outData, OBIT_IO_WriteOnly, err);

  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Add intent information from ASDM */
  ScanTab = SDMData->ScanTab;
  for (iScan=0; iScan<ScanTab->nrows; iScan++) {

    /* Timerange in human form */
    day2dhms(ScanTab->rows[iScan]->startTime-refJD, begString);
    day2dhms(ScanTab->rows[iScan]->endTime-refJD,   endString);
 
    g_snprintf (hicard, 80, "%s Scan=%d Source='%s' time= %s - %s", 
		pgmName, ScanTab->rows[iScan]->scanNumber,
		ScanTab->rows[iScan]->sourceName, begString, endString);
    ObitHistoryWriteRec (outHistory, -1, hicard, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);

    for (iIntent=0; iIntent<ScanTab->rows[iScan]->numIntent; iIntent++) {
      g_snprintf (hicard, 80, "%s   Intent[%d]='%s'", 
		  pgmName, iIntent+1,ScanTab->rows[iScan]->scanIntent[iIntent]);
      ObitHistoryWriteRec (outHistory, -1, hicard, err);
      if (err->error) Obit_traceback_msg (err, routine, outData->name);
    } /* end intent loop */
  } /* End scan look */
  
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  outHistory = ObitHistoryUnref(outHistory);  /* cleanup */
} /* End BDFInHistory */ 

/*----------------------------------------------------------------------- */
/* Write FG table entries for scans intended for online calibration       */
/*  e.g. pointing                                                         */
/*   Input:                                                               */
/*      SDMData   ASDM structure                                          */
/*      outData   ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void FlagIntent (ObitSDMData *SDMData, ObitUV* outData, ObitErr* err)
{
  ObitTableFG*     outTable=NULL;
  ObitTableFGRow*  outRow=NULL;
  olong            i, oRow, ver;
  ObitIOAccess     access;
  gboolean         flagIt;
  ASDMScanTable*   ScanTab;
  olong            iScan, iIntent, count=0;
  gchar            *flagEntries[] = {
    "CALIBRATE_POINTING",
    NULL};
  gchar *routine = "FlagIntent";
  
  /* error checks */
  if (err->error) return;
  g_assert (SDMData);
  g_assert (ObitUVIsA(outData));

  /* Create output FG table object */
  ver      = 1;
  access   = OBIT_IO_ReadWrite;
  outTable = newObitTableFGValue ("Output FG table", (ObitData*)outData, 
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
  outRow->SourID    = 0;
  outRow->SubA      = 0;
  outRow->freqID    = 0;
  outRow->ifs[0]    = 1;
  outRow->ifs[1]    = 0;
  outRow->chans[0]  = 1;
  outRow->chans[1]  = 0;
  outRow->pFlags[0] = -1;
  outRow->pFlags[1] = -1;
  outRow->pFlags[2] = -1;
  outRow->pFlags[3] = -1;
  outRow->status    = 0;
  
  /* Check intent information from ASDM */
  ScanTab = SDMData->ScanTab;
  for (iScan=0; iScan<ScanTab->nrows; iScan++) {

    flagIt = FALSE;
    i = 0;
    /* Is an intent for this scan on the flagEntries listed in the list? */
     for (iIntent=0; iIntent<ScanTab->rows[iScan]->numIntent; iIntent++) {
       while (flagEntries[i]) {
	 flagIt = flagIt || !strcmp(ScanTab->rows[iScan]->scanIntent[iIntent], flagEntries[i]);
	 strncpy (outRow->reason, ScanTab->rows[iScan]->scanIntent[iIntent], 24);
	 i++;
       }
     } /* end intent loop */
     
     if (flagIt) {
       /* Time range */
       outRow->TimeRange[0] =  ScanTab->rows[iScan]->startTime-refJD;
       outRow->TimeRange[1] =  ScanTab->rows[iScan]->endTime-refJD;
       
       
       /* Write */
       oRow = -1;
       count++;
       if ((ObitTableFGWriteRow (outTable, oRow, outRow, err)
	    != OBIT_IO_OK) || (err->error>0)) { 
	 Obit_log_error(err, OBIT_Error, "ERROR updating FG Table");
	 return;
       }
     } /* end write flag entry */
  } /* End scan loop */
  
  /* Close  table */
  if ((ObitTableFGClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output FG Table file");
    return;
  }
  
  /* Tell about it */
  Obit_log_error(err, OBIT_InfoErr, "Flagged %d online cal only scans", count);

  /* Cleanup */
  outRow   = ObitTableFGRowUnref(outRow);
  outTable = ObitTableFGUnref(outTable);

} /* end FlagIntent  */

void GetAntennaInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get info from ASDM structure                                          */
/*  Writes AN table output                                                */
/*  Missing antennas dummied                                              */
/*   Input:                                                               */
/*      SDMData  ASDM structure                                           */
/*      outData  Output UV object                                         */
/*      arrno     Array number                                            */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ASDMAntennaArray*  AntArray;
  ObitTableAN        *outTable=NULL;
  ObitTableANRow     *outRow=NULL;
  olong i, lim, iRow, oRow, ver, iarr;
  oint numIF, numPCal, numOrb;
  odouble JD, GASTM, Rate;
  ObitIOAccess access;
  gchar *out = "OUT";
  gchar *routine = "GetAntennaInfo";
  
  /* error checks */
  if (err->error) return;
  g_assert (SDMData!=NULL);
  g_assert (ObitUVIsA(outData));
  
  /* Extract antenna info */
  AntArray    = ObitSDMDataGetAntArray(SDMData, selScan);

  /* Create output Antenna table object */
  ver      = 1;
  access   = OBIT_IO_ReadWrite;
  numOrb   = 0;
  numPCal  = 0;
  numIF    = outData->myDesc->inaxes[outData->myDesc->jlocif];
  outTable = newObitTableANValue ("Output table", (ObitData*)outData, 
				  &ver, access, numIF, numOrb, numPCal, err);
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

  /* Is this the EVLA? */
  isEVLA = !strncmp(AntArray->arrayName, "EVLA", 4);

  /* Set AN table values */
  outTable->ArrayX  = 0.0;   /* Earth centered */
  outTable->ArrayY  = 0.0;
  outTable->ArrayZ  = 0.0;
  
  outTable->Freq    = outData->myDesc->crval[outData->myDesc->jlocf];
  /*MORE;*/
  outTable->PolarX  = 0.0;    /* ASDM Lacks */
  outTable->PolarY  = 0.0;
  outTable->dataUtc = 0.0;
  lim = MAXKEYCHARTABLEAN;
  strncpy (outTable->FRAME,   "ITRF    ", lim);        /* ASDM doesn't say */
  for (i=0; i<lim; i++) if (outTable->FRAME[i]==0) outTable->FRAME[i]=' ';
  if (isEVLA)  
    strncpy (outTable->TimeSys, "UTC     ", lim);        /* ASDM doesn't say */
  else
    strncpy (outTable->TimeSys, "IAT     ", lim);        /* ASDM doesn't say */
  for (i=0; i<lim; i++) if (outTable->TimeSys[i]==0) outTable->TimeSys[i]=' ';
  strncpy (outTable->ArrName, AntArray->arrayName, lim);
  for (i=0;i<lim;i++) if (outTable->ArrName[i]==0) outTable->ArrName[i] = ' ';
  outTable->FreqID = 0;
  outTable->iatUtc = 0.0;                        /* ASDM Lacks */
  strncpy (outTable->polType, "        ", lim);  /* ASDM Lacks */
  for (i=0;i<lim;i++) if (outTable->polType[i]==0) outTable->polType[i] = ' ';
  outTable->P_Refant = 0;
  outTable->P_Diff01 = 0.0;
  outTable->P_Diff02 = 0.0;
  outTable->P_Diff03 = 0.0;
  outTable->P_Diff04 = 0.0;
  outTable->P_Diff05 = 0.0;
  outTable->P_Diff06 = 0.0;
  outTable->P_Diff07 = 0.0;
  outTable->P_Diff08 = 0.0;
  /* Get reference date for array  */
  ObitUVDescJD2Date (SDMData->refJD, outTable->RefDate);

  /* Calculate earth rotation rate, GMST at UT midnight if not given */
  JD = SDMData->refJD;
  iarr = 0;
  arrayRefJDs[iarr] = JD;
  ObitUVDescJD2Date (JD, outTable->RefDate);
  ObitPrecessGST0 (JD, &GASTM, &Rate);
  outTable->DegDay  = Rate * 360.0;
  outTable->GSTiat0 = GASTM * 15.0;
  if (outTable->GSTiat0 < 0.0)  outTable->GSTiat0 += 360.0;;

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
  oRow = 1;
  for (iRow = 1; iRow<=AntArray->nants; iRow++) {

    /* Dummy missing antennas  */
    while (oRow<AntArray->ants[iRow-1]->antennaNo) {
      for (i=0; i<3; i++) outRow->AntName[i] = out[i];
      outRow->noSta     = oRow;
      outRow->status    = 0;
      outRow->mntSta    = 0;   /* ASDM lacks information */
      outRow->staXof    = 0.0;
      outRow->StaXYZ[0] = outTable->ArrayX;
      outRow->StaXYZ[1] = outTable->ArrayY;
      outRow->StaXYZ[2] = outTable->ArrayZ;
      for (i=0; i<numOrb; i++) outRow->OrbParm[i] = 0.0;
      
      if ((ObitTableANWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating ANTENNA Table");
	return;
      }
      oRow++;
    } /* end dummy loop */
    /* Set output Row */
    outRow->noSta          = AntArray->ants[iRow-1]->antennaNo;
    outRow->PolAngA        = 0.0;
    outRow->polTypeA       = g_strdup("    "); 
    outRow->polTypeA[0]    = AntArray->ants[iRow-1]->polnType[0];
    for (i=0; i<numPCal; i++) 
	outRow->PolCalA[i] = 0.0;
    if (AntArray->ants[iRow-1]->numPoln>1) {  /* Multiple polarizations */
      outRow->PolAngB      = 0.0;   
      outRow->polTypeB     = g_strdup("    "); 
      outRow->polTypeB[0]  = AntArray->ants[iRow-1]->polnType[1];
      for (i=0; i<numPCal; i++) 
	outRow->PolCalB[i] = 0.0;
    }
    for (i=0; i<outTable->myDesc->repeat[outTable->AntNameCol]; i++) 
      outRow->AntName[i] = ' ';
    lim = MIN(strlen(AntArray->ants[iRow-1]->staName), 
	      outTable->myDesc->repeat[outTable->AntNameCol]);
    for (i=0; i<lim; i++) 
      outRow->AntName[i] = AntArray->ants[iRow-1]->staName[i];
    outRow->status    = 0;
    outRow->mntSta    = 0;   /* ASDM lacks information */
    outRow->staXof    = AntArray->ants[iRow-1]->offset[0];  /* Not sure */
    outRow->StaXYZ[0] = AntArray->ants[iRow-1]->staPosition[0];
    outRow->StaXYZ[1] = AntArray->ants[iRow-1]->staPosition[1];
    outRow->StaXYZ[2] = AntArray->ants[iRow-1]->staPosition[2];
    for (i=0; i<numOrb; i++) outRow->OrbParm[i] = 0.0;

    if ((ObitTableANWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating ANTENNA Table");
      return;
    }
    oRow++;
  } /* end loop over input table */
  
  /* Close  table */
  if ((ObitTableANClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output Antenna Table file");
    return;
  }

  /* Cleanup */
  AntArray  = ObitSDMDataKillAntArray(AntArray);
  outRow    = ObitTableANRowUnref(outRow);
  outTable  = ObitTableANUnref(outTable);

} /* end GetAntennaInfo */

void GetFrequencyInfo (ObitSDMData *SDMData, ObitUV *outData, 
		       ASDMSpectralWindowArray* SpWinArray, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get Source info from ASDM                                             */
/*  Write FQ table                                                        */
/*   Input:                                                               */
/*      SDMData    ASDM structure                                         */
/*      outData    Output UV object                                       */
/*      SpWinArray Selected Spectral window array                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableFQ*            outTable=NULL;
  ObitTableFQRow*         outRow=NULL;
  olong i, oRow, ver, numIF, iIF;
  odouble refFreq=0.0;
  ObitIOAccess access;
  gchar *routine = "GetFrequencyInfo";

  /* error checks */
  if (err->error) return;
  g_assert (SDMData!=NULL);
  g_assert (ObitUVIsA(outData));

  /* Frequency info from first selected Spectral window */
  for (i=0; i<SpWinArray->nwinds; i++) {
    if (SpWinArray->winds[i]->selected) {
      refFreq   = SpWinArray->winds[i]->chanFreqStart;
      break;
    }
  }

  /* Count elected Spectral window */
  numIF = 0;
  for (i=0; i<SpWinArray->nwinds; i++) {
    if (SpWinArray->winds[i]->selected) numIF++;
  }
  /* Better be some */
  Obit_return_if_fail((numIF>=1), err,
		      "%s: NO SpectralWindows selected", routine);

  /* Create output FQ table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
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

  /* Save to FQ table */
  outRow->fqid    = 1;
  iIF  = 0;
  for (i=0; i<SpWinArray->nwinds; i++) {
    if (!SpWinArray->winds[i]->selected) continue; /* Want this one? */
    outRow->freqOff[iIF]  = SpWinArray->winds[i]->chanFreqStart - refFreq;
    outRow->chWidth[iIF]  = SpWinArray->winds[i]->chanWidth;
    outRow->totBW[iIF]    = SpWinArray->winds[i]->totBandwidth;
    if (SpWinArray->winds[i]->netSideband[0]=='U') outRow->sideBand[iIF] = 1;
    else outRow->sideBand[iIF] = -1;
    /* FOR EVLA everything is upper even if it's not */
    if (isEVLA) outRow->sideBand[iIF] = 1;
    iIF++;
  }
  outRow->status    = 0;
  oRow = outRow->fqid;
  if (oRow<1) oRow = -1;
  if ((ObitTableFQWriteRow (outTable, oRow, outRow, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR updating FG Table");
    return;
  }
  
  /* Close  table */
  if ((ObitTableFQClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output FQ Table file");
    return;
  }

  /* Cleanup */
  outRow     = ObitTableFQRowUnref(outRow);
  outTable   = ObitTableFQUnref(outTable);

} /* end  GetFrequencyInfo */

void GetSourceInfo (ObitSDMData *SDMData, ObitUV *outData, olong iScan,
		    ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get Source info from ASDM                                             */
/*  ASDM structure has one source entry per spectral window(=IF)          */
/*   Input:                                                               */
/*      SDMData  ASDM structure                                           */
/*      outData  Output UV object                                         */
/*      iScan    Scan number to use for Spectral Window array             */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ASDMSourceArray*   SourceArray;
  ASDMSpectralWindowArray* SpWinArray;
  ObitTableSU*       outTable=NULL;
  ObitTableSURow*    outRow=NULL;
  ObitSource *source=NULL;
  olong i, jSU, lim, iRow, oRow, ver, SourceID, lastSID=-1, iSW, SWId;
  oint numIF;
  ObitIOAccess access;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean *isDone=NULL, doCode=TRUE;
  gchar *blank = "        ";
  gchar *routine = "GetSourceInfo";

  /* error checks */
  if (err->error) return;
  g_assert (SDMData);
  g_assert (ObitUVIsA(outData));

  /* Print any messages */
  ObitErrLog(err);
  
  /* Give each source in source table by name a constant source number 
     include calcode or not? */
  ObitInfoListGetTest(myInput, "doCode", &type, dim, &doCode);
  if (doCode) 
    ObitSDMSourceTabFixCode(SDMData);
  else
    ObitSDMSourceTabFix(SDMData);

  /* Extract info */
  SourceArray = ObitSDMDataGetSourceArray(SDMData);
  Obit_return_if_fail((SourceArray), err,
		      "%s: Could not extract Source info from ASDM", 
		      routine);
  SpWinArray  = ObitSDMDataGetSWArray (SDMData, iScan, SWOrder);
  Obit_return_if_fail((SpWinArray), err,
		      "%s: Could not extract Spectral Windows from ASDM", 
		      routine);

  /* Create output Source table object */
  ver      = 1;
  access   = OBIT_IO_ReadWrite;
  numIF    = outData->myDesc->inaxes[outData->myDesc->jlocif];
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

  /* Init table header */
  outTable->FreqID = 1;
  lim = MIN(8,MAXKEYCHARTABLESU);
  for (i=0; i<lim; i++) outTable->velDef[i]  = blank[i];
  for (i=0; i<lim; i++) outTable->velType[i] = blank[i];

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
  for (i=0; i<numIF; i++) {
    outRow->IFlux[i]     = 0.0;
    outRow->QFlux[i]     = 0.0;
    outRow->UFlux[i]     = 0.0;
    outRow->VFlux[i]     = 0.0;
    outRow->FreqOff[i]   = 0.0;
    outRow->LSRVel[i]    = 0.0;
    outRow->RestFreq [i] = 0.0;
  }
  outRow->status    = 0;

  isDone = g_malloc0(SourceArray->nsou*sizeof(gboolean));
  for (iRow=0; iRow<SourceArray->nsou; iRow++) isDone[iRow] = FALSE;

  /* loop through input table */
  for (iRow=0; iRow<SourceArray->nsou; iRow++) {

    /* Done this one? */
    /*if ((isDone[iRow]) || (SourceArray->sou[iRow]->sourceNo<=lastSID)) continue;*/
    if (isDone[iRow]) continue;

    /* Is this one selected? - check in Spectral Window array */
    SWId = SourceArray->sou[iRow]->spectralWindowId;
    for (iSW=0; iSW<SpWinArray->nwinds; iSW++) {
      if (SpWinArray->winds[iSW]->spectralWindowId==SWId) break;
    }
    /* Not found or not selected? */
    if ((iSW>=SpWinArray->nwinds) || (!SpWinArray->winds[iSW]->selected)) 
      continue; 

    /* Set output row  */
    outRow->SourID    = SourceArray->sou[iRow]->sourceNo;
    isDone[iRow]      = TRUE;   /* Mark as done */
    outRow->RAMean    = SourceArray->sou[iRow]->direction[0]*RAD2DG;
    outRow->DecMean   = SourceArray->sou[iRow]->direction[1]*RAD2DG;
    outRow->PMRa      = SourceArray->sou[iRow]->properMotion[0]*RAD2DG*365.25;
    outRow->PMDec     = SourceArray->sou[iRow]->properMotion[1]*RAD2DG*365.25;
    outRow->Epoch     = 2000.0;   /* ASDM Lacking - AIPS naming is wrong */
    lastSID           = SourceArray->sou[iRow]->sourceNo;
    /* Precess */
    source = newObitSource("Temp");
    source->equinox = outRow->Epoch;
    source->RAMean  = outRow->RAMean;
    source->DecMean = outRow->DecMean;
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (outData->myDesc, source);
    outRow->RAApp  = source->RAApp;
    outRow->DecApp = source->DecApp;
    source = ObitSourceUnref(source);

    /* blank fill source name */
    lim = outTable->myDesc->repeat[outTable->SourceCol];
    for (i=0; i<lim; i++) outRow->Source[i] = ' ';
    lim = MIN(strlen(SourceArray->sou[iRow]->sourceName), 
	      outTable->myDesc->repeat[outTable->SourceCol]);
    for (i=0; i<lim; i++) outRow->Source[i] = SourceArray->sou[iRow]->sourceName[i];
    if (SourceArray->sou[iRow]->code) {
      lim = outTable->myDesc->repeat[outTable->CalCodeCol];
      for (i=0; i<lim; i++) outRow->CalCode[i] = ' ';
      lim = MIN(strlen(SourceArray->sou[iRow]->code), 
		outTable->myDesc->repeat[outTable->CalCodeCol]);
      for (i=0; i<lim; i++) outRow->CalCode[i]= SourceArray->sou[iRow]->code[i];
    } else {outRow->CalCode[0]=outRow->CalCode[1]=outRow->CalCode[2]=outRow->CalCode[3] = ' ';}

    /* Lines in first IF */
    if (SourceArray->sou[iRow]->sysVel)
      outRow->LSRVel[0]   = SourceArray->sou[iRow]->sysVel[0];
    else
      outRow->LSRVel[0]   = 0.0;
    if (SourceArray->sou[iRow]->restFrequency)
      outRow->RestFreq[0] = SourceArray->sou[iRow]->restFrequency[0];
    else
      outRow->RestFreq[0] = 0.0;
    outRow->FreqOff[0]  = 0.0;   /* ASDM lacking */

    /* Grumble, grumble, lookup velocity info for first line in other SWs */
    for (i=1; i<numIF; i++) {
      SourceID = SourceArray->sou[iRow]->sourceNo;
      for (jSU=iRow; jSU<SourceArray->nsou; jSU++) {
	if (SourceArray->sou[jSU]->sourceId==SourceID) break;
      }
      if (jSU<SourceArray->nsou) {
	isDone[iRow] = TRUE;    /* Mark as done */
	if (SourceArray->sou[jSU]->sysVel)
	  outRow->LSRVel[i]   = SourceArray->sou[jSU]->sysVel[0];
	else
	  outRow->LSRVel[i]   = 0.0;
	if (SourceArray->sou[jSU]->restFrequency)
	  outRow->RestFreq[i] = SourceArray->sou[jSU]->restFrequency[0];
	else
	  outRow->RestFreq[i] = 0.0;
	outRow->FreqOff[i]  = 0.0;   /* ASDM lacking */
      } else {  /* No more found */
	  outRow->LSRVel[i]   = 0.0;
  	  outRow->RestFreq[i] = 0.0;
	  outRow->FreqOff[i]  = 0.0;   /* ASDM lacking */
      }
    }
    outRow->status    = 0;
    
    oRow = -1; /*outRow->SourID;*/
    if ((ObitTableSUWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating Source Table");
      return;
    }
  } /* end loop over input table */
  
  if ((ObitTableSUClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output Source Table file");
    return;
  }

  /* Cleanup */
  SourceArray = ObitSDMDataKillSourceArray(SourceArray);
  SpWinArray  = ObitSDMDataKillSWArray (SpWinArray);
  outRow      = ObitTableSURowUnref(outRow);
  outTable    = ObitTableSUUnref(outTable);
  if (isDone) g_free(isDone);

} /* end  GetSourceInfo */

void GetData (ObitSDMData *SDMData, ObitInfoList *myInput, ObitUV *outData, 
	      ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Read data from BDF file, write outData                                */
/*      SDMData  ASDM structure                                           */
/*      myInput  parser object                                            */
/*      outData  Output UV object, open on input                          */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitBDFData *BDFData=NULL;
  ObitIOCode retCode;
  olong iMain, iInteg, ScanId=0, i, j, iBB, selChan, selIF, selConfig, 
    iSW, jSW, kBB, nIFsel, cntDrop=0, ver, iRow, sourId=0, iIntent, iScan;
  ofloat *Buffer=NULL, tlast=-1.0e20, startTime, endTime;
  ObitUVDesc *desc;
  ObitTableNX* NXtable;
  ObitTableNXRow* NXrow=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar selBand[12], begString[17], endString[17];
  ObitASDMBand band;
  ASDMSpectralWindowArray* SpWinArray=NULL;
  gboolean dropZero=TRUE, found=FALSE, doOnline=FALSE, drop;
  gchar dataroot[132];
  gchar *filename;
  gchar *routine = "GetData";

  /* error checks */
  if (err->error) return;
  g_assert(myInput!=NULL);
  g_assert(ObitUVIsA(outData));

  /* info */
  desc = outData->myDesc;
  if (desc->jlocs>=0) nstok = desc->inaxes[desc->jlocs];
  if (desc->jlocf>=0) nchan = desc->inaxes[desc->jlocf];
  if (desc->jlocif>=0) nIF  = desc->inaxes[desc->jlocif];
 
  /* Prepare output */
  desc   = outData->myDesc;
  Buffer = outData->buffer;
  desc->firstVis = desc->nvis+1;  /* Append to end of data */

  /* Set Obit sort order */
  if (outData->myDesc->isort[0]==' ') outData->myDesc->isort[0] = 'T';
  if (outData->myDesc->isort[1]==' ') outData->myDesc->isort[1] = 'B';
  
  /* Create BDF Structure  */
  BDFData = ObitBDFDataCreate ("BDF", outData->myDesc, SDMData, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Channel/IF/config selection - should have been completely specified in GetHeader */
  selChan = 0;
  ObitInfoListGetTest(myInput, "selChan", &type, dim, &selChan);
  selIF   = 0;
  ObitInfoListGetTest(myInput, "selIF", &type, dim, &selIF);
  selConfig = -1;
  ObitInfoListGetTest(myInput, "selConfig", &type, dim, &selConfig);

  /* Drop Zero vis? */
  ObitInfoListGetTest(myInput, "dropZero", &type, dim, &dropZero);

  /* Want Online scans? */
  ObitInfoListGetTest(myInput, "doOnline", &type, dim, &doOnline);
  
  /* Band selection */
  for (i=0; i<12; i++) selBand[i] = 0;
  ObitInfoListGetTest(myInput, "selBand", &type, dim, selBand);
  band = ObitSDMDataBand2Band (selBand);

  /* input DataRoot name */
  for (i=0; i<132; i++) dataroot[i] = 0;
  ObitInfoListGet(myInput, "DataRoot", &type, dim, dataroot, err);
  dataroot[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(dataroot);  /* Trim trailing blanks */

  /* Create Index table object */
  ver = 1;
  NXtable = newObitTableNXValue ("Index table", (ObitData*)outData, &ver, 
				 OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Clear existing rows */
  ObitTableClearRows ((ObitTable*)NXtable, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open Index table */
  if ((ObitTableNXOpen (NXtable, OBIT_IO_WriteOnly, err)
       != OBIT_IO_OK) || (err->error)) goto done;

  /* Create Index Row */
  NXrow = newObitTableNXRow (NXtable);

  /* initialize NX row */
  NXrow->SourID   = 0;
  NXrow->SubA     = 0;
  NXrow->Time     = 0.0;
  NXrow->TimeI    = 0.0;
  NXrow->StartVis = -1;
  NXrow->EndVis   = -1;
  NXrow->FreqID   = 1;

  /* attach to table buffer */
  ObitTableNXSetRow (NXtable, NXrow, err);
  if (err->error)  goto done;

  /* Write at beginning of iNdeX Table */
  iRow = 0;
  NXtable->myDesc->nrow = 0; /* ignore any previous entries */

  /* Loop over Main Table in ASDM */
  for (iMain=0; iMain<SDMData->MainTab->nrows; iMain++) {

    /* Selected? */
    if (selConfig != SDMData->MainTab->rows[iMain]->configDescriptionId) continue;
    if (SpWinArray) SpWinArray  = ObitSDMDataKillSWArray (SpWinArray);
    SpWinArray  = 
      ObitSDMDataGetSWArray (SDMData, SDMData->MainTab->rows[iMain]->scanNumber, 
			     SWOrder);
    Obit_return_if_fail((SpWinArray), err,
			"%s: Could not extract Spectral Windows from ASDM", 
			routine);
    /* Selection here mostly by ConfigID */
    if (!ObitSDMDataSelChan (SpWinArray, selChan, selIF, ASDMBand_Any)) continue;

    /* Ignore online cal scans unless doOnline */
    if (!doOnline) {
      drop = FALSE;
      /* DAMN ASDM - lookup ScanID */
      for (iScan=0; iScan<SDMData->ScanTab->nrows; iScan++) {
	ScanId = iScan;
	if ( SDMData->MainTab->rows[iMain]->scanNumber == 
	     SDMData->ScanTab->rows[ScanId]->scanNumber) break;
      }
      for (iIntent=0; iIntent<SDMData->ScanTab->rows[ScanId]->numIntent; iIntent++) {
	drop = drop ||  (!strncmp(SDMData->ScanTab->rows[ScanId]->scanIntent[iIntent], 
				  "CALIBRATE_POINTING", 18));
      }
      if (drop) {
	Obit_log_error(err, OBIT_InfoErr, "Drop online cal scan %3.3d subscan %3.3d", 
		       SDMData->MainTab->rows[iMain]->scanNumber, SDMData->MainTab->rows[iMain]->subscanNumber);
	ObitErrLog(err);
	continue;
      }
    }
    
    /* Get filename */
    filename = g_strconcat (dataroot, "/ASDMBinary/uid_", SDMData->MainTab->rows[iMain]->entityId, NULL);
    /* Complain and bail if file doesn't exist */
    if (!ObitFileExist(filename, err)) {
      /* Try adding .txt to end */
      if(filename) g_free(filename); 
      filename = g_strconcat (dataroot, "/ASDMBinary/uid_", SDMData->MainTab->rows[iMain]->entityId, ".txt",NULL);
      if (!ObitFileExist(filename, err)) {
	Obit_log_error(err, OBIT_InfoWarn, "Skipping, BDF file %s not found", filename);
	continue;
      }
    }
 
   /* File initialization */
    ObitBDFDataInitFile (BDFData, filename, err);
    g_free(filename);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* Init Scan */
    ObitBDFDataInitScan (BDFData, iMain, SWOrder, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
 
    /* Set selection - here mostly by ConfigID */
    ObitBDFDataSelChan (BDFData, selChan, selIF, ASDMBand_Any);

    /* Consistency check - loop over selected Spectral windows */
    nIFsel = 0;   /* Number of selected IFs */
    iSW    = 0;   /* input Spectral window index */
    kBB    = 0;   /* Baseband index */
    iBB    = 0;   /* SW Index in baseband */
    while (iSW<BDFData->SWArray->nwinds) {
      /* Get from ordered list */
      jSW = BDFData->SWArray->order[iSW];
      if (BDFData->SWArray->winds[jSW]->selected) {
	Obit_return_if_fail((nchan==BDFData->SWArray->winds[jSW]->numChan), err,
			    "%s: Input number freq. incompatible %d != %d", 
			    routine, nchan, BDFData->SWArray->winds[jSW]->numChan);
	Obit_return_if_fail((nstok==BDFData->SWArray->winds[jSW]->nCPoln), err,
			    "%s: Input number Poln incompatible %d != %d", 
			    routine, nstok, BDFData->SWArray->winds[jSW]->nCPoln);
	/* Baseband - assume name starts with "BB_" and ends in number 
	   BBNum = (olong)strtol(&BDFData->ScanInfo->BBinfo[kBB]->basebandName[3], NULL, 10);
	   Obit_return_if_fail((BBNum==BDFData->SWArray->winds[jSW]->basebandNum), err,
	   "%s: Input basebands inconsistent %d != %d, IF %d", 
	   routine, BBNum, BDFData->SWArray->winds[jSW]->basebandNum, nIFsel);*/
	/* Test frequency */
	Obit_return_if_fail((fabs(BDFData->SWArray->winds[iSW]->chanFreqStart-
				  outData->myDesc->freqIF[iSW]) < 1.0e3), err,
			    "%s: Frequencies inconsistent %lf != %lf, IF %d", 
			    routine, BDFData->SWArray->winds[iSW]->chanFreqStart, 
			    outData->myDesc->freqIF[iSW], nIFsel);
  	nIFsel++;
      } 
      iSW++;
      /* All of this baseband? */
      iBB++;
      if (iBB>=BDFData->ScanInfo->BBinfo[kBB]->numSpectralWindow) {kBB++; iBB=0;} /* Next baseband */
    } /* End loop over basebands/spectral windows */
    /* Same number of IFs? */
    Obit_return_if_fail((nIF==nIFsel), err,
			"%s: Input number Bands (IFs) incompatible %d != %d", 
			routine, nIF, nIFsel);

    /* Tell about scan */
    ScanId = SDMData->MainTab->rows[iMain]->scanNumber;
    for (j=0; j<SDMData->ScanTab->nrows; j++) {
      if (SDMData->ScanTab->rows[j]->scanNumber==ScanId) break;
    }
    if (j<SDMData->ScanTab->nrows) {
      /* Timerange in human form */
      day2dhms(SDMData->ScanTab->rows[j]->startTime-refJD, begString);
      day2dhms(SDMData->ScanTab->rows[j]->endTime-refJD,   endString);
      Obit_log_error(err, OBIT_InfoErr, "Scan %3.3d %s time %s  - %s", 
		     ScanId, SDMData->ScanTab->rows[j]->sourceName, begString, endString);
      ObitErrLog(err);
  }

    /* Initialize index */
    startTime = -1.0e20;
    endTime   = startTime;
    NXrow->StartVis = outData->myDesc->nvis+1;
    NXrow->EndVis   = NXrow->StartVis;

    /* Loop over integrations */
    for (iInteg=0; iInteg<SDMData->MainTab->rows[iMain]->numIntegration; iInteg++) {

      /* Read integration */
      retCode =  ObitBDFDataReadInteg(BDFData, err);
      if (retCode == OBIT_IO_EOF) break;
      if (err->error) Obit_traceback_msg (err, routine, outData->name);
      
      /* Loop over data */
      while (1) {
	/* fetch visibility */
	retCode =  ObitBDFDataGetVis (BDFData, Buffer, err);
	/* Done? */
	if (retCode == OBIT_IO_EOF) break;
	if (err->error) Obit_traceback_msg (err, routine, outData->name);
	
	/* Calculate uvw */
	CalcUVW (outData, BDFData, Buffer, err);
	if (err->error) Obit_traceback_msg (err, routine, outData->name);

	/* Check sort order */
	if (Buffer[desc->iloct]<tlast) 
	  if (Buffer[desc->iloct]<tlast) {
	    if (desc->isort[0]=='T') {
	      day2dhms(tlast, begString);
	      day2dhms(Buffer[desc->iloct], endString);
	      Obit_log_error(err, OBIT_InfoWarn, "Lost Time ordering at scan %d %s>%s", 
			     ScanId, begString, endString);
	    }
	    {desc->isort[0]=' '; desc->isort[1]=' ';}
	  }
	tlast = Buffer[desc->iloct];
	
	/* set number of records */
	desc->numVisBuff = 1;

	if (CheckAllZero(desc, Buffer) && dropZero) {
	  cntDrop++;
	  continue;
	}

	/* Get indexing information on first vis written */
	if (startTime<-1.0e10) {
	  startTime = Buffer[desc->iloct];
	  sourId    = (olong)(Buffer[desc->ilocsu]+0.5);
	}
	
	/* Write output */
	found   = TRUE;
	endTime = Buffer[desc->iloct];
	if ((ObitUVWrite (outData, NULL, err) != OBIT_IO_OK) || (err->error))
	  Obit_log_error(err, OBIT_Error, "ERROR writing output UV data"); 
	if (err->error) Obit_traceback_msg (err, routine, outData->name);
      } /* End loop over data in integration */
    } /* end loop over integrations */
    
    /* Write Index table */
    if (startTime>-1.0e10) {
      NXrow->Time     = 0.5 * (startTime + endTime);
      NXrow->TimeI    = (endTime - startTime);
      NXrow->EndVis   = desc->nvis;
      NXrow->SourID   = sourId;
      iRow++;
      if ((ObitTableNXWriteRow (NXtable, iRow, NXrow, err)
	   != OBIT_IO_OK) || (err->error>0)) goto done; 
    }
  } /* End loop over scans */

  /* Tell results */
 done:
  Obit_log_error(err, OBIT_InfoErr, 
		 "Wrote %d visibilities, drop %d all zero", 
		 outData->myDesc->nvis, cntDrop);
  ObitErrLog(err);

  /* Cleanup */
  if ((ObitTableNXClose (NXtable, err) != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, NXtable->name);

  BDFData  = ObitBDFDataUnref (BDFData);
  if (SpWinArray) SpWinArray  = ObitSDMDataKillSWArray (SpWinArray);
  NXrow    = ObitTableNXRowUnref(NXrow);
  NXtable   = ObitTableNXUnref(NXtable);

} /* end GetData  */

void CalcUVW (ObitUV *outData, ObitBDFData *BDFData, ofloat *Buffer, 
	      ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Calculates u,v,w (lambda) for a uv data row                           */
/*  Does nothing for autocorrelations                                     */
/*  Short baseline approximation                                          */
/*   Input:                                                               */
/*      outData  Output UV object                                         */
/*      BDFData  BDF data structure                                       */
/*      Buffer   Input uv data row, u,v,w will be modified                */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  olong souId, ant1, ant2, i, iANver, iarr, cnt;
  ObitTableAN *ANTable=NULL;
  ObitTableSU *SUTable=NULL;
  ObitSource *source=NULL;
  ObitAntennaList *AntList;
  ofloat time, uvw[3], bl[3], tmp, u, v;
  odouble arrJD, DecR, RAR, AntLst, HrAng=0.0, ArrLong, ArrLat;
  odouble sum, xx, yy, zz, lambda, DecOff, RAOff;
  gchar *routine = "CalcUVW";

  /* error checks */
  if (err->error) return;

  /* Which antennas? */
  tmp = outData->buffer[outData->myDesc->ilocb];
  ant1 = (olong)(tmp/256.0);
  ant2 = (olong)(tmp - ant1*256 +0.00005);
  
  /* NOP for autocorrelations  */
  if (ant1==ant2) return;

  /* Get antenna lists if they don't already exist */
  if (antennaLists==NULL) {

    /* Array of subarrays for antenna longitudes/latitudes */
    antLongs = g_malloc0(numArray*sizeof(ofloat*));
    antLats  = g_malloc0(numArray*sizeof(ofloat*));

    /* Create AntennaLists */
    antennaLists = g_malloc0(numArray*sizeof(ObitAntennaList*));
    for (i=0; i<numArray; i++) {
      iANver = i+1;
      ANTable = newObitTableANValue ("AN table", (ObitData*)outData, &iANver, 
				     OBIT_IO_ReadOnly, 0, 0, 0, err);
      antennaLists[i] = ObitTableANGetList (ANTable, err);
      if (err->error) Obit_traceback_msg (err, routine, outData->name);
      /* release table object */
      ANTable = ObitTableANUnref(ANTable);

      /* Subarray list of antenna longitudes/latitudes */
      antLongs[i] = g_malloc0(antennaLists[i]->number*sizeof(ofloat));
      antLats[i]  = g_malloc0(antennaLists[i]->number*sizeof(ofloat));

      /* Get average longitude */
      AntList = antennaLists[i];
      sum = 0.0; cnt = 0;
      for (iarr=0; iarr<AntList->number; iarr++) {
	antLongs[i][iarr] = atan2 (AntList->ANlist[iarr]->AntXYZ[1], 
				  AntList->ANlist[iarr]->AntXYZ[0]);
	antLats[i][iarr]  = AntList->ANlist[iarr]->AntLat;
	if (fabs(AntList->ANlist[iarr]->AntXYZ[1])>1.0) {
	  sum += antLongs[i][iarr];
	  cnt++;
	}
      }
      ArrLong = sum / cnt;
      AntList->ANlist[0]->AntLong = ArrLong;

      lambda = VELIGHT/BDFData->desc->freq;  /* Reference wavelength */

      ArrLong = -ArrLong;   /* Other way for rotation */
      
      /* Convert positions to lambda - rotate to frame of the array */
      for (iarr=0; iarr<AntList->number; iarr++) {
	xx = AntList->ANlist[iarr]->AntXYZ[0] / lambda;
	yy = AntList->ANlist[iarr]->AntXYZ[1] / lambda;
	zz = AntList->ANlist[iarr]->AntXYZ[2] / lambda;
	AntList->ANlist[iarr]->AntXYZ[0] = xx*cos(ArrLong) - yy*sin(ArrLong);
	AntList->ANlist[iarr]->AntXYZ[1] = xx*sin(ArrLong) + yy*cos(ArrLong);
	AntList->ANlist[iarr]->AntXYZ[2] = zz;
      }
    }
  } /* end create antenna lists */

  /* Get source information if necessary */
  if (uvwSourceList==NULL) {
    SUTable = newObitTableSUValue ("SU table", (ObitData*)outData, &iANver, 
				   OBIT_IO_ReadOnly, 0, err);
    uvwSourceList = ObitTableSUGetList (SUTable, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
   /* release table object */
    SUTable = ObitTableSUUnref(SUTable);
    /* Make sure all have precessed positions */
    for (i=0; i<uvwSourceList->number; i++) {
      ObitPrecessUVJPrecessApp (outData->myDesc, uvwSourceList->SUlist[i]);
    }
  } /* end create source list */

  /* New source? */
  souId = (olong)Buffer[BDFData->desc->ilocsu];
  if (uvwSourID!=souId ) {
    uvwSourID = souId;
    /* Find in Source List */
    for (i=0; i<uvwSourceList->number; i++) {
      uvwcurSourID = i;
      if (uvwSourceList->SUlist[i]->SourID==uvwSourID) break;
    }
    
    /* Get rotation to get v north at standard epoch - 
       precess posn. 10" north */
    source = newObitSource("Temp");
    source->RAMean  = uvwSourceList->SUlist[uvwcurSourID]->RAMean;
    source->DecMean = uvwSourceList->SUlist[uvwcurSourID]->DecMean + 10.0/3600.0;
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (outData->myDesc, source);
    RAOff  = source->RAApp;
    DecOff = source->DecApp;
    source = ObitSourceUnref(source);
    /* uvrot global = rotation to north */
    uvrot = -(ofloat)atan2(RAOff-uvwSourceList->SUlist[uvwcurSourID]->RAApp,
			   DecOff-uvwSourceList->SUlist[uvwcurSourID]->DecApp);
  } /* end new source */
  
    /* Array number (0-rel)*/
  iarr = 0;
  /* Array reference JD */
  arrJD = arrayRefJDs[iarr];

  /* Array geometry  */
  AntList = antennaLists[iarr];
  ArrLong = 0.5*(antLongs[iarr][ant1-1] + antLongs[iarr][ant2-1]);
  ArrLat  = 0.5*(antLats[iarr][ant1-1]  + antLats[iarr][ant2-1]);
  
  bl[0] =  AntList->ANlist[ant1-1]->AntXYZ[0] - AntList->ANlist[ant2-1]->AntXYZ[0];
  bl[1] =  AntList->ANlist[ant1-1]->AntXYZ[1] - AntList->ANlist[ant2-1]->AntXYZ[1];
  bl[2] =  AntList->ANlist[ant1-1]->AntXYZ[2] - AntList->ANlist[ant2-1]->AntXYZ[2];

  /* Position in radians */
  RAR  = uvwSourceList->SUlist[uvwcurSourID]->RAApp*DG2RAD;
  DecR = uvwSourceList->SUlist[uvwcurSourID]->DecApp*DG2RAD;

  /* LST and hour angle (radians) */
  time   = Buffer[BDFData->desc->iloct];
  AntLst = AntList->GSTIAT0 + ArrLong + time*AntList->RotRate;
  HrAng  = AntLst - RAR;

  /* Compute uvw - short baseline approximation */
  ObitUVUtilUVW (bl, DecR, (ofloat)HrAng, uvw);

  /* Rotate in u-v plane to north of standard epoch */
  u = uvw[0];
  v = uvw[1];
  uvw[0] = u*cos(uvrot) - v*sin(uvrot);
  uvw[1] = v*cos(uvrot) + u*sin(uvrot);

  Buffer[BDFData->desc->ilocu] = uvw[0];
  Buffer[BDFData->desc->ilocv] = uvw[1];
  Buffer[BDFData->desc->ilocw] = uvw[2];

} /* end CalcUVW */

void GetFlagInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any flag info from ASDM to AIPS FG on outData                 */
/*  ASDM only has antenna/time flagging                                   */
/*   Input:                                                               */
/*      SDMData  ASDM structure                                           */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableFG*       outTable=NULL;
  ObitTableFGRow*    outRow=NULL;
  olong              iRow, oRow, ver, iarr, antId, antNo, iAnt, i;
  ObitIOAccess       access;
  ASDMAntennaArray*  AntArray;
  gchar              *routine = "GetFlagInfo";

  /* error checks */
  if (err->error) return;
  g_assert (SDMData);
  g_assert (ObitUVIsA(outData));

  /* Any entries? */
  if (SDMData->FlagTab->nrows<=0) return;

  /* Print any messages */
  ObitErrLog(err);
  
  /* Create output FG table object */
  ver      = 1;
  access   = OBIT_IO_ReadWrite;
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
  
  /* Antenna array to lookup antenna numbers */
  AntArray    = ObitSDMDataGetAntArray(SDMData, selScan);
  Obit_return_if_fail((AntArray), err,
		      "%s: Could not extract Antenna info from ASDM", 
		      routine);

  /* Create output Row */
  outRow = newObitTableFGRow (outTable);
  /* attach to table buffer */
  ObitTableFGSetRow (outTable, outRow, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Initialize output row */
  outRow->SourID    = 0;
  outRow->SubA      = 0;
  outRow->freqID    = 0;
  outRow->ifs[0]    = 1;
  outRow->ifs[1]    = 0;
  outRow->chans[0]  = 1;
  outRow->chans[1]  = 0;
  outRow->pFlags[0] = -1;
  outRow->pFlags[1] = -1;
  outRow->pFlags[2] = -1;
  outRow->pFlags[3] = -1;
  outRow->status    = 0;
  
  /* Array number 0-rel */
  iarr = 0;

  /* loop through input table */
  for (iRow =0; iRow<SDMData->FlagTab->nrows; iRow++) {

     /* Make sure valid */
    if (SDMData->FlagTab->rows[iRow]->reason==NULL) continue;

   /* Look up antenna number from Id */
    antId = SDMData->FlagTab->rows[iRow]->antennaId;
    antNo = antId;
    for (iAnt=0; iAnt<AntArray->nants; iAnt++) {
      if (AntArray->ants[iAnt]->antennaId==antId) 
	{antNo = AntArray->ants[iAnt]->antennaNo; break;}
    }
    
    outRow->ants[0]      = antNo;
    outRow->ants[1]      = 0;
	 outRow->TimeRange[0] = SDMData->FlagTab->rows[iRow]->startTime - SDMData->refJD;
    outRow->TimeRange[1] = SDMData->FlagTab->rows[iRow]->endTime   - SDMData->refJD;
    strncpy (outRow->reason, SDMData->FlagTab->rows[iRow]->reason, 24);
    /* Patch for archaic software */
    for (i=0; i<24; i++) if (outRow->reason[i]==0) outRow->reason[i] = ' ';

    /* Write */
    oRow = -1;
    if ((ObitTableFGWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating FG Table");
      return;
    }
  } /* end loop over input table */
  
  /* Close  table */
  if ((ObitTableFGClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output FG Table file");
    return;
  }
  
  /* Tell about it */
  Obit_log_error(err, OBIT_InfoErr, "Copied %d flag records", outTable->myDesc->nrow);

  /* Cleanup */
  outRow   = ObitTableFGRowUnref(outRow);
  outTable = ObitTableFGUnref(outTable);
  AntArray = ObitSDMDataKillAntArray(AntArray);

} /* end  GetFlagInfo */

void GetWeatherInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert any WEATHER table from ASDM  to AIPS WX on outData            */
/*   Input:                                                               */
/*      SDMData  ASDM structure                                           */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableWX*          outTable=NULL;
  ObitTableWXRow*       outRow=NULL;
  ASDMWeatherTable*     inTab=SDMData->WeatherTab;
  olong iRow, oRow, ver;
  odouble mjdJD0=2400000.5; /* JD of beginning of MJD time */
  odouble K0 = 273.15;      /* Zero point of Centigrade in Kelvin */
  ofloat fblank = ObitMagicF();
  ObitIOAccess access;
  gchar *routine = "GetWeatherInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outData));

  /* Any entries? */
  if (inTab->nrows<=0) return;

  /* Print any messages */
  ObitErrLog(err);
  
  /* Create output WX table object */
  ver      = 1;
  access   = OBIT_IO_ReadWrite;
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
  for (iRow=0; iRow<inTab->nrows; iRow++) {
    
    /* Make sure valid */
    if (inTab->rows[iRow]->timeInterval==NULL) continue;

    /* Save to WX table */
    outRow->Time          = inTab->rows[iRow]->timeInterval[0]-refJD;
    outRow->TimeI         = inTab->rows[iRow]->timeInterval[1] -
                            inTab->rows[iRow]->timeInterval[0];
    /* Correct interval if necessary */
    if (inTab->rows[iRow]->timeInterval[1]<inTab->rows[iRow]->timeInterval[0])
      outRow->TimeI = inTab->rows[iRow]->timeInterval[1] - mjdJD0;
    outRow->antNo         = inTab->rows[iRow]->stationId;
    outRow->temperature   = inTab->rows[iRow]->temperature - K0;  /* K -> C */
    outRow->pressure      = inTab->rows[iRow]->pressure*0.01;     /* Pascal to millibar */
    outRow->dewpoint      = inTab->rows[iRow]->dewPoint - K0;     /* K -> C */
    outRow->windVelocity  = inTab->rows[iRow]->windSpeed;
    outRow->windDirection = inTab->rows[iRow]->windDirection*RAD2DG; /* rad to deg */
    outRow->wvrH2O        = fblank;
    outRow->onosElectron  = fblank;
    /* Write */
    oRow = -1;
    if ((ObitTableWXWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating WX Table");
      return;
    }
    
  } /* end loop over input table */
  
    /* Close  table */
  if ((ObitTableWXClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output WX Table file");
    return;
  }
  
  /* Tell about it */
  Obit_log_error(err, OBIT_InfoErr, "Copied WEATHER table");
  
  /* Cleanup */
  outRow   = ObitTableWXRowUnref(outRow);
  outTable = ObitTableWXUnref(outTable);

} /* end  GetWeatherInfo */

void GetCalDeviceInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert the CalDevice table from ASDM  to AIPS CD on outData          */
/*  Entries are accumulated into AIPS CD table                            */
/*   Input:                                                               */
/*      SDMData  ASDM structure                                           */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableCD*          outTable=NULL;
  ObitTableCDRow*       outRow=NULL;
  ASDMcalDeviceTable*   inTab=SDMData->calDeviceTab;
  ASDMAntennaArray*    AntArray;
  ASDMSpectralWindowArray* SpWinArray;
  olong i, j, iRow, oRow, ver, maxAnt, iAnt, jAnt, IFno;
  olong *antLookup, *SpWinLookup;
  oint numIF, numPol;
  ofloat fblank = ObitMagicF();
  gboolean want;
  ObitIOAccess access;
  gchar *routine = "GetCalDeviceInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outData));

  /* Any entries? */
  if (inTab->nrows<=0) return;

  /* Print any messages */
  ObitErrLog(err);
  
  /* Extract ASDM SpWin data  - selScan global */
  SpWinArray  = ObitSDMDataGetSWArray (SDMData, selScan, SWOrder);
  Obit_return_if_fail((SpWinArray), err,
		      "%s: Could not extract Spectral Windows from ASDM", 
		      routine);
  /* Extract antenna info */
  AntArray    = ObitSDMDataGetAntArray(SDMData, selScan);
  Obit_return_if_fail((AntArray), err,
		      "%s: Could not extract Antenna info from ASDM", 
		      routine);

  /* Select Spectral windows - use global selChan, selIF */
  ObitSDMDataSelChan (SpWinArray, selChan, selIF, ASDMBand_Any);

  /* Highest antenna number? */
  maxAnt = AntArray->maxAnt;
  numIF  = outData->myDesc->inaxes[outData->myDesc->jlocif];
  numPol = MIN (2, outData->myDesc->inaxes[outData->myDesc->jlocs]);
  
  /* Antenna number lookup table */
  antLookup = g_malloc(maxAnt*sizeof(olong));
  for (i=0; i<maxAnt; i++) antLookup[i] = -1;  /* For missing ants */
  for (i=0; i<AntArray->nants; i++) {
    if ((AntArray->ants[i]->antennaId>=0) && (AntArray->ants[i]->antennaId<maxAnt))
      antLookup[AntArray->ants[i]->antennaId] = AntArray->ants[i]->antennaNo;
  }

  /* Spectral window/IF lookup table */ 
  SpWinLookup  = g_malloc0(SDMData->SpectralWindowTab->nrows*sizeof(olong));
  for (i=0; i<SDMData->SpectralWindowTab->nrows; i++) SpWinLookup[i] = -1;  /* For deselected */
  for (i=0; i<SpWinArray->nwinds; i++) {
    if ((SpWinArray->winds[i]->spectralWindowId>=0) && 
	(SpWinArray->winds[i]->spectralWindowId<SDMData->SpectralWindowTab->nrows) &&
	SpWinArray->winds[i]->selected)
      SpWinLookup[SpWinArray->winds[i]->spectralWindowId] = SpWinArray->order[i];
  }

  /* Create output CD table object */
  ver      = 1;
  access   = OBIT_IO_ReadWrite;
  outTable = newObitTableCDValue ("Output table", (ObitData*)outData, 
				  &ver, access, numIF, numPol, err);
  if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with CD table");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Open table */
  if ((ObitTableCDOpen (outTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output CD table");
    return;
  }
  
  /* Table Header */
  ObitUVDescJD2Date (SDMData->refJD, outTable->RefDate);

  /* Create output Row */
  outRow = newObitTableCDRow (outTable);
  /* attach to table buffer */
  ObitTableCDSetRow (outTable, outRow, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Initialize output row */
  outRow->status      = 0;
  outRow->SubA        = 1;
  outRow->FreqID      = 1;

  /* Loop over antennas collecting SpWin/IFs for each antenna */
  for (iAnt=1; iAnt<=maxAnt; iAnt++) {
    outRow->antennaNo  = iAnt;

    /* Blank data in case missing */
    for (i=0; i<numIF; i++) {
      outRow->TCal1[i] = fblank;
      if (numPol>1) 
	outRow->TCal2[i] = fblank;
  }
    
    /* loop through input table */
    for (iRow=0; iRow<inTab->nrows; iRow++) {

      /* Make sure valid */
      if (inTab->rows[iRow]->timeInterval==NULL) continue;

      /* Is this the desired antenna? */
      if ((inTab->rows[iRow]->antennaId>=0) && 
	  (inTab->rows[iRow]->antennaId<AntArray->nants)) 
	jAnt = antLookup[inTab->rows[iRow]->antennaId];
      else jAnt = -1;
      if (jAnt!=iAnt) continue;

      /* Get 0-rel IF number - Sp Win must be in current config and selected */
      want = FALSE; 
      for (j=0; j<SpWinArray->nwinds; j++) {
	if ((inTab->rows[iRow]->spectralWindowId==SpWinArray->winds[j]->spectralWindowId) &&
	  SpWinArray->winds[j]->selected) {want=TRUE; break;}
      }

      if (want) {
	IFno = SpWinLookup[inTab->rows[iRow]->spectralWindowId];
	IFno = MAX (0, MIN(IFno, (numIF-1)));
      } else continue;
      
      /* Save to CD table row */
      /* Apply efficiencies to cal temp */
      if (inTab->rows[iRow]->calEff) {
	outRow->TCal1[IFno] = (ofloat)(inTab->rows[iRow]->noiseCal[0]*inTab->rows[iRow]->calEff[0]);
	if (numPol>1) 
	  outRow->TCal2[IFno] = (ofloat)(inTab->rows[iRow]->noiseCal[0]*inTab->rows[iRow]->calEff[1]);
      } else { /* - No efficiency - just cal */
 	outRow->TCal1[IFno] = (ofloat)(inTab->rows[iRow]->noiseCal[0]);
	if (numPol>1) 
	  outRow->TCal2[IFno] = (ofloat)(inTab->rows[iRow]->noiseCal[0]);
      }
    } /* end loop over input table */

    /* Write */
    oRow = -1;
    if ((ObitTableCDWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating CD Table");
      return;
    }
    
  } /* end Antenna loop */

  /* Close  table */
  if ((ObitTableCDClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output CD Table file");
    return;
  }
  
  /* Tell about it */
  Obit_log_error(err, OBIT_InfoErr, "Copied calDevice table %d rows",
		 outTable->myDesc->nrow);
  
  /* Cleanup */
  SpWinArray = ObitSDMDataKillSWArray (SpWinArray);
  AntArray   = ObitSDMDataKillAntArray(AntArray);
  outRow     = ObitTableCDRowUnref(outRow);
  outTable   = ObitTableCDUnref(outTable);
  if (antLookup)   g_free(antLookup);
  if (SpWinLookup) g_free(SpWinLookup);

} /* end  GetCalDeviceInfo */

void GetSysPowerInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert the SysPower table from ASDM  to AIPS SY on outData           */
/*  Entries are accumulated into AIPS SY table                            */
/*  Input SysPower table is modified (antennaId)                          */
/*   Input:                                                               */
/*      SDMData  ASDM structure                                           */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableSY*          outTable=NULL;
  ObitTableSYRow*       outRow=NULL;
  ASDMSysPowerTable*     inTab=SDMData->SysPowerTab;
  ASDMAntennaArray*    AntArray;
  ASDMSpectralWindowArray* SpWinArray;
  olong i, j, iRow, jRow, oRow, ver, maxAnt, IFno, SourNo, SWId;
  olong *antLookup=NULL, *SpWinLookup=NULL, *SpWinLookup2=NULL;
  olong curScan, curScanI, nextScanNo;
  oint numIF, numPol;
  ofloat fblank = ObitMagicF();
  gboolean want;
  ObitIOAccess access;
  gchar *routine = "GetSysPowerInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outData));

  /* Any entries? */
  if (inTab->nrows<=0) return;

  /* Print any prior messages */
  ObitErrLog(err);
  
  /* Extract ASDM SpWin data  - selScan global */
  curScan    = selScan;
  curScanI   = -1;  /* Force init */
  SourNo     = 0;
  nextScanNo = 0;
  SpWinArray  = ObitSDMDataGetSWArray (SDMData, selScan, SWOrder);
  Obit_return_if_fail((SpWinArray), err,
		      "%s: Could not extract Spectral Windows from ASDM", 
		      routine);
  /* Extract antenna info */
  AntArray    = ObitSDMDataGetAntArray(SDMData, selScan);
  Obit_return_if_fail((AntArray), err,
		      "%s: Could not extract Antenna info from ASDM", 
		      routine);

  /* Select Spectral windows - use global selChan, selIF */
  ObitSDMDataSelChan (SpWinArray, selChan, selIF, ASDMBand_Any);

  /* Highest antenna number? */
  maxAnt = AntArray->maxAnt;
  numIF  = outData->myDesc->inaxes[outData->myDesc->jlocif];
  numPol = MIN (2, outData->myDesc->inaxes[outData->myDesc->jlocs]);

  /* Antenna number lookup table */
  antLookup = g_malloc(maxAnt*sizeof(olong));
  for (i=0; i<maxAnt; i++) antLookup[i] = -1;  /* For missing ants */
  for (i=0; i<AntArray->nants; i++) {
    if ((AntArray->ants[i]->antennaId>=0) && (AntArray->ants[i]->antennaId<maxAnt))
      antLookup[AntArray->ants[i]->antennaId] = AntArray->ants[i]->antennaNo;
  }

  /* Spectral window/IF lookup tables */ 
  /* Lookup[SWId] = if order number, -1=>unused  */
  SpWinLookup  = g_malloc0(SDMData->SpectralWindowTab->nrows*sizeof(olong));
  /* Lookup2[SWId] = SpWinArray element for that SW, -1=not */
  SpWinLookup2 = g_malloc0(SDMData->SpectralWindowTab->nrows*sizeof(olong));
  for (i=0; i<SDMData->SpectralWindowTab->nrows; i++) SpWinLookup[i]  = -1;  /* For deselected */
  for (i=0; i<SDMData->SpectralWindowTab->nrows; i++) SpWinLookup2[i] = -1;
  for (i=0; i<SpWinArray->nwinds; i++) {
    SpWinLookup2[SpWinArray->winds[i]->spectralWindowId] = i;
    if ((SpWinArray->winds[i]->spectralWindowId>=0) && 
	(SpWinArray->winds[i]->spectralWindowId<SDMData->SpectralWindowTab->nrows) &&
	SpWinArray->winds[i]->selected)
      SpWinLookup[SpWinArray->winds[i]->spectralWindowId] = SpWinArray->order[i];
  }

  /* Create output SY table object */
  ver      = 1;
  access   = OBIT_IO_ReadWrite;
  outTable = newObitTableSYValue ("Output table", (ObitData*)outData, 
				  &ver, access, numIF, numPol, err);
  if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with SY table");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Open table */
  if ((ObitTableSYOpen (outTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output SY table");
    return;
  }

  /* Header values */
  outTable->nAnt = maxAnt;
  
  /* Create output Row */
  outRow = newObitTableSYRow (outTable);
  /* attach to table buffer */
  ObitTableSYSetRow (outTable, outRow, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Have to collect multiple SpWin into one record */
  
  /* Initialize output row */
  outRow->status      = 0;
  outRow->SubA        = 1;
  outRow->FreqID      = 1;
  /* Blank data in case missing */
  for (i=0; i<numIF; i++) {
    outRow->PwrDif1[i] = fblank;
    outRow->PwrSum1[i] = fblank;
    outRow->Gain1[i]   = fblank;
    if (numPol>1) {
      outRow->PwrDif2[i] = fblank;
      outRow->PwrSum2[i] = fblank;
      outRow->Gain2[i]   = fblank;
    }
  }
    
  /* loop through input table */
  for (iRow=0; iRow<inTab->nrows; iRow++) {

    /* Make sure valid */
    if (inTab->rows[iRow]->timeInterval==NULL) continue;

    /* Look for previously handled data (antennaId=-10) */
    if (inTab->rows[iRow]->antennaId<=-10) continue;

    /* Is this in the same scan?  Antennas may change but not SpWin */
    if (nextScan(SDMData, curScan, inTab->rows[iRow]->timeInterval[0], 
		&curScanI, &nextScanNo, &SourNo)) {
      curScan = nextScanNo;

      /* Extract antenna info */
      AntArray = ObitSDMDataKillAntArray (AntArray);  /* Delete old */
      AntArray = ObitSDMDataGetAntArray(SDMData, curScan);
      Obit_return_if_fail((AntArray), err,
			  "%s: Could not extract Antenna info from ASDM", 
			  routine);

      /* Antenna number lookup table */
      antLookup = g_malloc(maxAnt*sizeof(olong));
      for (i=0; i<maxAnt; i++) antLookup[i] = -1;  /* For missing ants */
      for (i=0; i<AntArray->nants; i++) {
	if ((AntArray->ants[i]->antennaId>=0) && (AntArray->ants[i]->antennaId<maxAnt))
	  antLookup[AntArray->ants[i]->antennaId] = AntArray->ants[i]->antennaNo;
      }
    } /* End new scan */
      
    /* Save info to SY table row */
    outRow->Time          = inTab->rows[iRow]->timeInterval[0]-refJD;
    outRow->TimeI         = inTab->rows[iRow]->timeInterval[1];
    outRow->SourID        = SourNo;
    /* Convert antennaID to antenna number */
    if ((inTab->rows[iRow]->antennaId>=0) && 
	(inTab->rows[iRow]->antennaId<AntArray->nants)) 
      outRow->antennaNo = antLookup[inTab->rows[iRow]->antennaId];
    else continue;  /* ignore if antennaId bad */

    /* Get 0-rel IF number - Sp Win must be in current config and selected */
    want = FALSE; 
    for (j=0; j<SpWinArray->nwinds; j++) {
      if ((inTab->rows[iRow]->spectralWindowId==SpWinArray->winds[j]->spectralWindowId) &&
	  SpWinArray->winds[j]->selected) {want=TRUE; break;}
    }
    
    if (want) {
      IFno = SpWinLookup[inTab->rows[iRow]->spectralWindowId];
      IFno = MAX (0, MIN(IFno, (numIF-1)));
    } else continue;
    
    /* snatch data */
    outRow->PwrDif1[IFno] = inTab->rows[iRow]->switchedPowerDifference[0];
    outRow->PwrSum1[IFno] = inTab->rows[iRow]->switchedPowerSum[0];
    outRow->Gain1[IFno]   = inTab->rows[iRow]->requantizerGain[0];
    if (numPol>1) {
      outRow->PwrDif2[IFno] = inTab->rows[iRow]->switchedPowerDifference[1];
      outRow->PwrSum2[IFno] = inTab->rows[iRow]->switchedPowerSum[1];
      outRow->Gain2[IFno]   = inTab->rows[iRow]->requantizerGain[1];
    }

    /* Search table until find entry for same antennaId with a later time */
    for (jRow=iRow+1; jRow<inTab->nrows; jRow++) {
      /* Gone far enough? */
      if ((inTab->rows[jRow]->timeInterval[0]>inTab->rows[iRow]->timeInterval[0])) break;
      /* DEBUG 
	 if ((inTab->rows[jRow]->timeInterval[0]>inTab->rows[iRow]->timeInterval[0]) &&
	 (inTab->rows[jRow]->antennaId==inTab->rows[iRow]->antennaId)) break; */

      /* Look for previously handled data (antennaId=-10) */
      if (inTab->rows[jRow]->antennaId<=-10) continue;

      /* Desired antenna? */
      if (inTab->rows[jRow]->antennaId!=inTab->rows[iRow]->antennaId) continue;

      /* Must want this one - work out IF number - must be valid and selected */
      SWId = SpWinLookup2[inTab->rows[jRow]->spectralWindowId];
      if ((SWId>=0) && (SWId<SpWinArray->nwinds) &&  SpWinArray->winds[SWId]->selected) {
	IFno = SpWinLookup[inTab->rows[jRow]->spectralWindowId];
	IFno = MAX (0, MIN(IFno, (numIF-1)));
      } else continue;
      
      /* snatch data */
      outRow->PwrDif1[IFno] = inTab->rows[jRow]->switchedPowerDifference[0];
      outRow->PwrSum1[IFno] = inTab->rows[jRow]->switchedPowerSum[0];
      outRow->Gain1[IFno]   = inTab->rows[jRow]->requantizerGain[0];
      if (numPol>1) {
	outRow->PwrDif2[IFno] = inTab->rows[jRow]->switchedPowerDifference[1];
	outRow->PwrSum2[IFno] = inTab->rows[jRow]->switchedPowerSum[1];
	outRow->Gain2[IFno]   = inTab->rows[jRow]->requantizerGain[1];
      }
      /* Mark table row as done */
      inTab->rows[jRow]->antennaId = -10;
    } /* end loop looking for rest of the IFs */
    /* Write */
    oRow = -1;
    if ((ObitTableSYWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating SY Table");
      return;
    }
    /* Blank data in case missing */
    for (i=0; i<numIF; i++) {
      outRow->PwrDif1[i] = fblank;
      outRow->PwrSum1[i] = fblank;
      outRow->Gain1[i]   = fblank;
      if (numPol>1) {
	outRow->PwrDif2[i] = fblank;
	outRow->PwrSum2[i] = fblank;
	outRow->Gain2[i]   = fblank;
      }
    }
    
  } /* end loop over input table */
  
    /* Close  table */
  if ((ObitTableSYClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output SY Table file");
    return;
  }
  
  /* Tell about it */
  Obit_log_error(err, OBIT_InfoErr, "Copied SysPower table %d rows",
		 outTable->myDesc->nrow);
  
  /* Cleanup */
  outRow   = ObitTableSYRowUnref(outRow);
  outTable = ObitTableSYUnref(outTable);
  ObitSDMDataKillSWArray (SpWinArray);
  ObitSDMDataKillAntArray (AntArray);
  if (antLookup)    g_free(antLookup);
  if (SpWinLookup)  g_free(SpWinLookup);
  if (SpWinLookup2) g_free(SpWinLookup2);

} /* end  GetSysPowerInfo */

void GetOTTInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert the Pointing table from ASDM  to AIPS OT on outData           */
/*   Input:                                                               */
/*      SDMData  ASDM structure                                           */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableOT*          outTable=NULL;
  ObitTableOTRow*       outRow=NULL;
  ASDMPointingTable*    inTab=SDMData->PointingTab;
  ASDMAntennaArray*     AntArray;
  olong i, iRow, oRow, ver, maxAnt, SourNo;
  olong *antLookup, curScan, curScanI, nextScanNo;
  ObitIOAccess access;
  gchar *routine = "GetSysPowerInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outData));

  /* Any entries? */
  if (inTab->nrows<=0) return;

  /* Print any prior messages */
  ObitErrLog(err);
  
  /* Extract ASDM Antenna data  - selScan global */
  curScan    = selScan;
  curScanI   = -1;  /* Force init */
  SourNo     = 0;
  nextScanNo = 0;
  AntArray    = ObitSDMDataGetAntArray(SDMData, selScan);
  Obit_return_if_fail((AntArray), err,
		      "%s: Could not extract Antenna info from ASDM", 
		      routine);

  /* Highest antenna number? */
  maxAnt = AntArray->maxAnt;

  /* Antenna number lookup table */
  antLookup = g_malloc(maxAnt*sizeof(olong));
  for (i=0; i<maxAnt; i++) antLookup[i] = -1;  /* For missing ants */
  for (i=0; i<AntArray->nants; i++) {
    if ((AntArray->ants[i]->antennaId>=0) && (AntArray->ants[i]->antennaId<maxAnt))
      antLookup[AntArray->ants[i]->antennaId] = AntArray->ants[i]->antennaNo;
  }

  /* Create output OT table object */
  ver      = 1;
  access   = OBIT_IO_ReadWrite;
  outTable = newObitTableOTValue ("Output table", (ObitData*)outData, 
				  &ver, access, err);
  if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with OT table");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Open table */
  if ((ObitTableOTOpen (outTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OT table");
    return;
  }
  
  /* Create output Row */
  outRow = newObitTableOTRow (outTable);
  /* attach to table buffer */
  ObitTableOTSetRow (outTable, outRow, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Have to collect multiple SpWin into one record */
  
  /* Initialize output row */
  outRow->status      = 0;
    
  /* loop through input table */
  for (iRow=0; iRow<inTab->nrows; iRow++) {

    /* Make sure valid */
    if (inTab->rows[iRow]->timeInterval==NULL) continue;

    /* Look for previously handled data (antennaId=-10) */
    if (inTab->rows[iRow]->antennaId<=-10) continue;

    /* Is this in the same scan?  Antennas may change but not SpWin */
    if (nextScan(SDMData, curScan, inTab->rows[iRow]->timeInterval[0], 
		&curScanI, &nextScanNo, &SourNo)) {
      curScan = nextScanNo;

      /* Extract antenna info */
      AntArray = ObitSDMDataKillAntArray (AntArray);  /* Delete old */
      AntArray = ObitSDMDataGetAntArray(SDMData, curScan);
      Obit_return_if_fail((AntArray), err,
			  "%s: Could not extract Antenna info from ASDM", 
			  routine);

      /* Antenna number lookup table */
      antLookup = g_malloc(maxAnt*sizeof(olong));
      for (i=0; i<maxAnt; i++) antLookup[i] = -1;  /* For missing ants */
      for (i=0; i<AntArray->nants; i++) {
	if ((AntArray->ants[i]->antennaId>=0) && (AntArray->ants[i]->antennaId<maxAnt))
	  antLookup[AntArray->ants[i]->antennaId] = AntArray->ants[i]->antennaNo;
      }
    } /* End new scan */
      
    /* Save info to OT table row */
    outRow->Time          = inTab->rows[iRow]->timeInterval[0]-refJD;
    outRow->TimeI         = inTab->rows[iRow]->timeInterval[1];
    outRow->SourID        = SourNo;
    /* Convert antennaID to antenna number */
    if ((inTab->rows[iRow]->antennaId>=0) && 
	(inTab->rows[iRow]->antennaId<AntArray->nants)) 
      outRow->Antenna = antLookup[inTab->rows[iRow]->antennaId];
    else continue;  /* ignore if antennaId bad */
    /* Over the top? */
    outRow->OverTop = inTab->rows[iRow]->overTheTop;

    /* Write */
    oRow = -1;
    if ((ObitTableOTWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating OT Table");
      return;
    }
  } /* end loop over input table */
  
    /* Close  table */
  if ((ObitTableOTClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OT Table file");
    return;
  }
  
  /* Tell about it */
  Obit_log_error(err, OBIT_InfoErr, "Copied OT table");
  
  /* Cleanup */
  outRow   = ObitTableOTRowUnref(outRow);
  outTable = ObitTableOTUnref(outTable);
  ObitSDMDataKillAntArray (AntArray);
  if (antLookup)   g_free(antLookup);

} /* end  GetSysPowerInfo */

void GetGainCurveInfo (ObitSDMData *SDMData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Generate AIPS GC on outData                                           */
/*  Gives band average Stokes I gain curves                               */
/*   Input:                                                               */
/*      SDMData  ASDM structure                                           */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableGC*       outTable=NULL;
  ObitTableGCRow*    outRow=NULL;
  ASDMAntennaArray*  AntArray;
  olong i, j, iAnt, Ant, oRow, ver;
  ofloat gains[10], sens, fblank = ObitMagicF();
;
  odouble refJD, Freq;
  oint numIF, numPol, numTabs;
  ObitIOAccess access;
  gchar *routine = "GetGainCurveInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outData));

   /* Extract antenna info */
  AntArray    = ObitSDMDataGetAntArray(SDMData, selScan);
  refJD       = SDMData->refJD;
  Freq        = outData->myDesc->crval[outData->myDesc->jlocf];
  sens        = nomSen(AntArray); /* get nominal sensitivity */

 /* Create output GC table object */
  ver      = 1;
  access   = OBIT_IO_ReadWrite;
  numIF    = outData->myDesc->inaxes[outData->myDesc->jlocif];
  numPol   = MIN (2, outData->myDesc->inaxes[outData->myDesc->jlocs]);
  numTabs  = 4;   /* 4 values in gain curves */
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
  
  /* over antennas */
  for (iAnt=1; iAnt<=AntArray->nants; iAnt++) {

    Ant = AntArray->ants[iAnt-1]->antennaNo;  /* Actual antenna number */

    /* get curve */
    ObitVLAGainParseGain (Ant, refJD, Freq, gains);
    
    /* Save to GC table */
    outRow->antennaNo  = Ant;
    outRow->SubArray   = 1;
    outRow->FreqID     = 1;
    for (i=0; i<numIF; i++) {
      outRow->Type1[i] = 2;
      outRow->NTerm1[i]= numTabs;
      outRow->XTyp1[i] = 2;
      outRow->YTyp1[i] = 2;
      outRow->XVal1[i] = fblank;
      outRow->sens1[i] = sens;
      for (j=0; j<numTabs; j++) {
	outRow->YVal1[i*numTabs+j] = fblank;
	outRow->gain1[i*numTabs+j] = gains[j];
      }
    }
    if (numPol>1) {   /* 2 poln */
      for (i=0; i<numIF; i++) {
	outRow->Type2[i] = 2;
	outRow->NTerm2[i]= numTabs;
	outRow->XTyp2[i] = 2;
	outRow->YTyp2[i] = 2;
	outRow->XVal2[i] = fblank;
	outRow->sens2[i] = sens;
	for (j=0; j<numTabs; j++) {
	  outRow->YVal2[i*numTabs+j] = fblank;
	  outRow->gain2[i*numTabs+j] = gains[j];
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
    
    /* Close  table */
    if ((ObitTableGCClose (outTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing output GC Table file");
      return;
    }

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Generated GC table %d", ver);

    /* Cleanup */
    outRow   = ObitTableGCRowUnref(outRow);
    outTable = ObitTableGCUnref(outTable);
    ObitSDMDataKillAntArray (AntArray);
} /* end  GetGainCurveInfo */

/** 
 * Convert Time in days to a human readable form "dd/hh:mm:ss.s"
 * \param time  Time in days
 * \param timeString [out] time as string, should be >16 char
 */
void day2dhms(ofloat time, gchar *timeString)
{
  olong day, thour, tmin;
  ofloat ttim, ssec;

  /* Trap bad times */
  if ((time<-100.0) || (time>1000.0)) {
    sprintf (timeString, "Bad time");
    return;
  }

  day   = (olong)(time);
  ttim  = 24.0*(time - day);
  thour = MIN ((olong)(ttim), 23);
  ttim  = 60.0*(ttim - thour);
  tmin  = MIN ((olong)(ttim), 59);
  ssec  = 60.0*(ttim - tmin);
  /* avoid silliness */
  if (ssec>59.951) {
    tmin++;
    ssec = 0.0;
  }
  if (tmin>=60) {
    thour++;
    tmin -= 60;
  }
  if (thour>23) {
    day++;
    thour -= 24;
  }
  sprintf (timeString, "%2.2d/%2.2d:%2.2d:%4.1f", day, thour, tmin, ssec);
  /* Zero fill seconds */
  if (timeString[9]==' ') timeString[9] = '0';
} /* end day2dhms */

/** 
 * Check uv buffer for zero visibilities
 * Any visibilities which are exactly zero will have weights set to zero
 * Returns TRUE if all visibilities are exactly zero
 * \param desc
 * \param Buffer
 * \return TRUE if all visibilities are exactly zero
 */
gboolean CheckAllZero (ObitUVDesc *desc, ofloat *Buffer)
{
  gboolean out=TRUE;
  olong i, off;

  for (i=0; i<desc->ncorr; i++) {
    off = desc->nrparm + i*3;
    if ((Buffer[off]==0.0) && (Buffer[off+1]==0.0)) Buffer[off+2] = 0.0;
    else out = FALSE;
  }

  return out;
} /* end CheckAllZero */

/*----------------------------------------------------------------------- */
/*  Is a given time in the current scan?                                  */
/*   Input:                                                               */
/*      SDMData   ASDM structure                                          */
/*      curScan   Current scan number                                     */
/*       time     Time to test in JD                                      */
/*       curScanI Current scan index in ScanTab, updated on new scan      */
/*       nextScan Scan number of next scan, updated on new scan           */
/*       SourNo   Source number for next Scan, updated on new scan        */
/*   Return:                                                              */
/*       TRUE if time in another scan                                     */
/*----------------------------------------------------------------------- */
gboolean nextScan(ObitSDMData *SDMData, olong curScan, odouble time,  
		  olong *curScanI, olong *nextScanNo, olong *SourNo)
{
  gboolean out = FALSE;
  olong iScan = *curScanI;
  olong i, iMain, fieldId, jField, sourceId, jSource;

  /* New scan? */
  if ((iScan>=0) && (time<=SDMData->ScanTab->rows[iScan]->endTime)) return out;
  out = TRUE;

  /* Find scan whose end time greater than time - use last if not found */
  iScan = MAX(0, iScan);
  for (i=iScan; i<SDMData->ScanTab->nrows; i++) {
    if (time<=SDMData->ScanTab->rows[i]->endTime) break;
  }
  iScan = MIN (i, SDMData->ScanTab->nrows-1);
  *curScanI   = iScan;
  *nextScanNo = SDMData->ScanTab->rows[iScan]->scanNumber;

  /* Now lookup source */
  /* Find scan number in MainTab */
  for (iMain=0; iMain<SDMData->MainTab->nrows; iMain++) {
    if (SDMData->MainTab->rows[iMain]->scanNumber==(*nextScanNo)) break;
  }

  /* Source Id - have to look down goddamn tree */
  fieldId = SDMData->MainTab->rows[iMain]->fieldId;
  for (jField=0; jField<SDMData->FieldTab->nrows; jField++) {
    if (SDMData->FieldTab->rows[jField]->fieldId==fieldId) break;
  }
  sourceId = SDMData->FieldTab->rows[jField]->sourceId;
  for (jSource=0; jSource<SDMData->SourceTab->nrows; jSource++) {
    if (SDMData->SourceTab->rows[jSource]->sourceId==sourceId) break;
  }
  *SourNo = SDMData->SourceTab->rows[jSource]->sourceNo;

  return out;
} /* end nextScan */

/*----------------------------------------------------------------------- */
/*  Get nominal sensitivity if EVLA                                       */
/*  Values as per EVLA web site Sep. 2010 (some interpolation)            */
/*   Input:                                                               */
/*     AntArray  Antenna information array                                */
/*   Return:                                                              */
/*       Nominal sensitivity in K/Jy if EVLA, else fblank                 */
/*----------------------------------------------------------------------- */
ofloat nomSen(ASDMAntennaArray*  AntArray)
{
  ofloat out = ObitMagicF();
  olong i, band;
   /* Number of bands tabulated */
  olong nband=10;
 /* Upper frequency of bands */
  odouble upFreq[] = 
    /*  4        P      L       S      C       X      U        K       A      Q */
    {100.0e6, 900.0e6, 2.0e9, 3.7e9, 7.5e9, 12.0e9, 18.0e9, 26.5e9, 40.0e9, 60.0e9};
  /* Band sensitivity EVLA web site Sep 2010 */
  ofloat bandSen[] = 
    /*  4(Ha)     P      L     S(?)    C       X      U        K       A(?)    Q */
    {0.071,   0.071,  0.098,  0.12,  0.123, 0.112,  0.103,   0.071,  0.066, 0.062};

  /* Only for EVLA */
  if (strncmp(AntArray->arrayName, "EVLA", 4)) return out;
 
  /* Lookup band */
  band = 0;
  for (i=0; i<nband; i++) {
    band = i;
    if (AntArray->refFreq<upFreq[i]) break;
  }
  
  return bandSen[band];
} /* end nomSen  */
