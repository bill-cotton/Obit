/* $Id$  */
/* Read GBT Quadrant detector data                                */
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

#include "ObitUtil.h"
#include "ObitOTF.h"
#include "ObitTableOTFSoln.h"
#include "ObitIOOTFFITS.h"
#include "ObitFITS.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitTimeFilter.h"
#include "ObitTableOTFTarget.h"
#include "ObitTableOTFIndex.h"
#include "ObitTableGBTANTPOSGR.h"
#include "ObitGBTIFInfo.h"
#include "ObitTableGBTQUADDETECTOR.h"
#include "ObitTableGBTPARSENSOR.h"
#include "ObitTableGBTPARSENSORUtil.h"
#define NUMDETECT 72   /* Maximum number of detectors */
/* internal prototypes */
/* Get inputs */
ObitInfoList* QuadDetin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Get pointing times from Antenna file, Secondary focus */
void GetAntennaGR (gchar *infile, ObitInfoList *myInput, ObitErr *err);
/* Get pointing for a given time */
void GetElev (odouble time, ofloat *el);
/* Get file descriptor */
gboolean GetHeader (ObitOTF *outData, char *outfile, gchar *infile, 
		    ObitInfoList *myInput, ObitErr *err);
/* Get Quadrant detector data */
void GetQuadData(gchar *inscan, ObitInfoList *myInput, ObitErr *err);
/* Convert data and write */
void CnvrtData (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, 
		ObitErr *err);
/* Get target id */
olong GetTarget (ObitOTF *outData, gboolean isNew, gchar *name, 
		odouble ra, odouble dec, odouble equinox, ObitErr *err);
/* Initialize Index table */
void InitScan (ObitOTF *outData, gboolean isNew, ofloat scan, ofloat target, 
	       ObitErr *err);
/* Undate scan info */
void SetScan (ObitOTF *outData, odouble startTime, odouble endTime, 
	      olong startRec, olong endRec, ObitErr *err);
/* Convert measurements to offsets */
void ProcessQDData (ObitErr *err);
/* Get pointing offsetfor a given time */
void GetPointOffset (odouble time, ofloat *dAz, ofloat *dEl);

/* Program globals */
gchar *pgmName = "QuadDet2";      /* Program name */
gchar *infile  = "QuadDet2.inp";  /* File with program inputs */
gchar *outfile = "QuadDet2.out";  /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
/* gchar **FITSdirs=NULL; List of FITS data directories */
gchar DataRoot[128]; /* Root directory of input data */
gchar Name[48];   /* Target name */
odouble targRA;   /* Target RA in deg */
odouble targDec;  /* Target Dec in deg */
odouble Equinox;  /* Equinox of targRA, targDec */
olong nAntTime;        /* number of antenna time samples */
odouble *AntDMJD=NULL; /* Array of antenna times */
odouble *AntAz=NULL;   /* Array of Antenna Azimuth  values */
odouble *AntEl=NULL;   /* Array of Antenna Elevation values */
olong nQDetTime;       /* number of Quadrant detector time samples */
odouble *QDetDMJD=NULL;/* Array of  Quadrant detector times */
ofloat *QDetV1=NULL;   /* Array of  Quadrant detector voltage 1 */
ofloat *QDetV3=NULL;   /* Array of  Quadrant detector voltage 3 */
ofloat *QDetV4=NULL;   /* Array of  Quadrant detector voltage 4 */
ofloat *QDetV5=NULL;   /* Array of  Quadrant detector voltage 5 */
ofloat *QDetdAz=NULL;  /* Array of  Quadrant detector Az offset, in arcseconds */
ofloat *QDetdEl=NULL;  /* Array of  Quadrant detector El offset, in arcseconds */
odouble refMJD;   /* reference Julian date */
odouble integTime;/* Integration time in days */
odouble startTime;/* Start time of Scan in days */
odouble endTime;  /* End time of Scan in days */
ofloat target;    /* target number */
ofloat scan;      /* scan number */
olong  startRec;  /* First record number (1-rel) in scan */
olong  endRec;    /* End record number (1-rel) in scan */
olong  nfeed=NUMDETECT;  /* Number of feeds */
olong  nchan=1;   /* Number of frequencies */
olong  nstok=1;   /* Number of Stokes */
olong   ndetect=NUMDETECT;/* Number of detectors */
ofloat deltaFreq; /* Channel increment */
ofloat refPixel;  /* Frequency reference pixel */
odouble refFrequency; /* reference frequency (Hz) */
ofloat azOff, elOff; /* az and el offsets to apply to pointing positions */

olong  *xlate=NULL;     /* Translation from pixel offset array to data array */
/* Note: must flip sign of el offset */
/* Beam offset, (cross el,elevation, col, row),  asec -9999.0 => no data*/
  /*               col     row  comment */
ofloat tbeamOffset[] = {  
  -9999.,  -9999.,   0.00,    8.00, /*  00  Dark pixel */
  -9999.,  -9999.,   1.00,    8.00, /*  01  Dark pixel */
  -9999.,  -9999.,   2.00,    8.00, /*  02  Dark pixel */
  -9999.,  -9999.,   3.00,    8.00, /*  03  Dark pixel */
  -9999.,  -9999.,   4.00,    8.00, /*  04  Dark pixel */
  -9999.,  -9999.,   5.00,    8.00, /*  05  Dark pixel */
  -9999.,  -9999.,   6.00,    8.00, /*  06  Dark pixel */
  -9999.,  -9999.,   7.00,    8.00, /*  07  Dark pixel */
  14.30,  -3.69,   0.00,    0.00, /*  08    */
  -0.74,  -2.86,   1.00,    0.00, /*  09    */
  -3.69, -14.30,   2.00,    0.00, /*  10    */
 -10.05,  -3.50,   3.00,    0.00, /*  11    */
 -20.01,   5.17,   4.00,    0.00, /*  12    */
   0.74,   2.86,   5.00,    0.00, /*  13    */
   3.69,  14.30,   6.00,    0.00, /*  14    */
   2.86,  -0.74,   7.00,    0.00, /*  15    */
  16.42,  -7.29,   0.00,    1.00, /*  16    */
   1.38,  -6.46,   1.00,    1.00, /*  17    */
  -7.29, -16.42,   2.00,    1.00, /*  18    */
 -13.65,  -5.63,   3.00,    1.00, /*  19    */
 -16.42,   7.29,   4.00,    1.00, /*  20    */
  -1.38,   6.46,   5.00,    1.00, /*  21    */
   7.29,  16.42,   6.00,    1.00, /*  22    */
   6.46,   1.38,   7.00,    1.00, /*  23    */
   8.58,  -2.21,   0.00,    2.00, /*  24    */
   5.63, -13.65,   1.00,    2.00, /*  25    */
  -2.21,  -8.58,   2.00,    2.00, /*  26    */
  -6.46,  -1.38,   3.00,    2.00, /*  27    */
 -10.70,   5.81,   4.00,    2.00, /*  28    */
  -5.63,  13.65,   5.00,    2.00, /*  29    */
   2.21,   8.58,   6.00,    2.00, /*  30    */
  13.65,   5.63,   7.00,    2.00, /*  31    */
  10.70,  -5.81,   0.00,    3.00, /*  32    */
   3.50, -10.05,   1.00,    3.00, /*  33    */
  -5.81, -10.70,   2.00,    3.00, /*  34    */
 -15.77,  -2.03,   3.00,    3.00, /*  35    */
 -12.82,   9.41,   4.00,    3.00, /*  36    */
  -3.50,  10.05,   5.00,    3.00, /*  37    */
   5.81,  10.70,   6.00,    3.00, /*  38    */
  10.05,   3.50,   7.00,    3.00, /*  39    */
  12.82,  -9.41,   0.00,    4.00, /*  40    */
   2.03, -15.77,   1.00,    4.00, /*  41    */
  -9.41, -12.82,   2.00,    4.00, /*  42    */
 -12.17,   0.09,   3.00,    4.00, /*  43    */
  -7.10,   7.93,   4.00,    4.00, /*  44    */
   0.09,  12.17,   5.00,    4.00, /*  45    */
   9.41,  12.82,   6.00,    4.00, /*  46    */
  15.77,   2.03,   7.00,    4.00, /*  47    */
   4.98,  -4.34,   0.00,    5.00, /*  48    */
  -0.09, -12.17,   1.00,    5.00, /*  49    */
  -4.34,  -4.98,   2.00,    5.00, /*  50    */
  -8.58,   2.21,   3.00,    5.00, /*  51    */
  -9.22,  11.53,   4.00,    5.00, /*  52    */
  -2.03,  15.77,   5.00,    5.00, /*  53    */
   4.34,   4.98,   6.00,    5.00, /*  54    */
  12.17,  -0.09,   7.00,    5.00, /*  55    */
   9.22, -11.53,   0.00,    6.00, /*  56    */
  -1.57, -17.89,   1.00,    6.00, /*  57    */
 -11.53,  -9.22,   2.00,    6.00, /*  58    */
 -17.89,   1.57,   3.00,    6.00, /*  59    */
  -4.98,   4.34,   4.00,    6.00, /*  60    */
   1.57,  17.89,   5.00,    6.00, /*  61    */
  11.53,   9.22,   6.00,    6.00, /*  62    */
  17.89,  -1.57,   7.00,    6.00, /*  63    */
   7.10,  -7.93,   0.00,    7.00, /*  64    */
  -5.17, -20.01,   1.00,    7.00, /*  65    */
  -7.93,  -7.10,   2.00,    7.00, /*  66    */
 -14.30,   3.69,   3.00,    7.00, /*  67    */
  -2.86,   0.74,   4.00,    7.00, /*  68    */
   5.17,  20.01,   5.00,    7.00, /*  69    */
   7.93,   7.10,   6.00,    7.00, /*  70    */
  20.01,  -5.17,   7.00,    7.00}; /*  71    */

ofloat *beamOffset = tbeamOffset;
odouble BTime = 2.50e-6;  /* blank 220 msec after cal switch */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*    Read GBT Quadrant detector files and convert to OTFSoln file        */
/*----------------------------------------------------------------------- */
{
  olong  i, disk, nrec, ierr=0;
  ObitSystem *mySystem= NULL;
  ObitOTF *outData= NULL;
  ObitIOAccess access;
  gchar *FITSdirs[] = {"./", NULL};
  ObitErr *err= NULL;
  gchar inscan[128], outOTF[128];
  ObitInfoType type;
  gboolean isNew;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  err = newObitErr();

  /* Startup - parse command line */
  ierr = 0;
  myInput = QuadDetin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input Scan name */
  for (i=0; i<128; i++) inscan[i] = 0;
  ObitInfoListGet(myInput, "Scan", &type, dim, inscan, err);
  inscan[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(inscan);  /* Trim trailing blanks */

  /* output FITS file name */
  for (i=0; i<128; i++) outOTF[i] = 0;
  ObitInfoListGet(myInput, "outOTF", &type, dim, outOTF, err);
  outOTF[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(outOTF);  /* Trim trailing blanks */

  /* Get input data file root */
  for (i=0; i<128; i++) DataRoot[i] = 0;
  ObitInfoListGet(myInput, "DataRoot", &type, dim, DataRoot, err);
  DataRoot[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(DataRoot);  /* Trim trailing blanks */

  /* Initialize Obit */
  FITSdirs[1] =  DataRoot;  /* Input FITS data directory */
  nFITS = 2;
  mySystem = ObitSystemStartup ("QuadDet", 1, 0, 0, NULL, nFITS, FITSdirs, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;  ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Create ObitOTF for data */
  outData = newObitOTF("Output data");

  /* show any errors */
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
   
  /* Get header info, array geometry, initialize output if necessary */
  isNew = GetHeader (outData, outOTF, inscan, myInput, err);
   /* Say what went wrong if error */
  if (err->error) 
     Obit_log_error(err, OBIT_Error, "Error reading Header file for scan %s", inscan);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
 
  /* Get Telescope position array */
  GetAntennaGR (inscan, myInput, err);
  /* Say what went wrong if error */
  if (err->error) 
     Obit_log_error(err, OBIT_Error, "Error reading Antenna file for scan %s", inscan);
  
  /* Get Quadrant detector data arrays */
  GetQuadData (inscan, myInput, err);
  /* Say what went wrong if error */
  if (err->error) 
     Obit_log_error(err, OBIT_Error, "Error reading Quadrant detector file for scan %s", inscan);
  
  /* Informative message 
     show any errors */
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
   
  /* Define output, I/O size */
  disk = 1;
  nrec = 1;
  ObitOTFSetFITS(outData,nrec,disk,outOTF,err);
  
  /* show any errors */
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
   
  /* Open output OTF */
  if (isNew) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  if ((ObitOTFOpen (outData, access, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", outOTF);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Get target id number */
  target = (ofloat)GetTarget (outData, isNew, Name, targRA, targDec, Equinox, err);

  /* convert data  */
  CnvrtData (outData, inscan, myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Initialize scan in Index table - 
     do not make entry before data successfully converted */
  InitScan (outData, isNew, scan, target, err);

  /* Update index table */
  SetScan (outData, startTime, endTime, startRec, endRec, err);

  /* show any errors */
   if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;
   
  /* Close */
  if ((ObitOTFClose (outData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing output file");
  
  /* show any errors */
   if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;
   
   /* Shutdown Obit */
 exit:
   ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
   mySystem = ObitSystemShutdown (mySystem);
   
   /* cleanup */
   myInput = ObitInfoListUnref(myInput);  /* delete input list */
   outData = ObitUnref(outData);
   if (AntDMJD)  g_free(AntDMJD);
   if (AntAz)    g_free(AntAz);
   if (AntEl)    g_free(AntEl);
   if (QDetDMJD) g_free(QDetDMJD);
   if (QDetV1)   g_free(QDetV1);
   if (QDetV3)   g_free(QDetV3);
   if (QDetV4)   g_free(QDetV4);
   if (QDetV5)   g_free(QDetV5);
   if (QDetdAz)  g_free(QDetdAz);
   if (QDetdEl)  g_free(QDetdEl);
 
   return ierr;
} /* end of main */

ObitInfoList* QuadDetin (int argc, char **argv, ObitErr *err)
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
  olong ax;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar *input_file="QuadDet2.in", *arg;
  gboolean init=FALSE;
  oint itemp;
  gchar *strTemp;
  ObitInfoList* list;
  gchar *routine = "QuadDetIn";

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

    } else if (strcmp(arg, "-Scan") == 0){ /* Scan name */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "Scan", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-DataRoot") == 0){ /* Data root directory */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "DataRoot", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outOTF") == 0){ /* Output OTF file */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outOTF", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outDisk") == 0) { /* Output FITS disk */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outDisk", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-pgmNumber") == 0) { /*Program number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "pgmNumber", OBIT_oint, dim, &itemp, err);

    } else if (strcmp(arg, "-AIPSuser") == 0) { /* AIPS User */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outDisk", OBIT_oint, dim, &itemp, err);
      
      
    } else { /* unknown argument */
      fprintf(stderr,"Unknown parameter %s\n",arg);
      Usage();
    }
  } /* end parsing input arguments */
  
  /* Read defaults if no file specified */
  if (!init) ObitParserParse (input_file, list, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  /* Extract basic information to program globals */
  ObitInfoListGet(list, "pgmNumber", &type, dim, &pgmNumber, err);
  /*ObitInfoListGet(list, "nFITS",     &type, dim, &nFITS,     err);*/
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  /* Initialize output */
  myOutput = defaultOutputs(err);
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

 return list;
} /* end QuadDetin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: QuadDet2 -input file [-Scan date/time]\n");
    fprintf(stderr, "Convert a GBT Quadrant format to Obit/OTF OTFSoln table\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def QuadDet.in\n");
    fprintf(stderr, "  -scan date/time used for form scan FITS file names\n");
    fprintf(stderr, "  -config configuration file \n");
    
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
  strTemp = "DataOTF.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "outOTF", OBIT_string, dim, strTemp, err);

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
/*  Values:                                                               */
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

void GetAntennaGR (gchar *inscan, ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get antenna positions and leave in globals, Secondary focus (>1 GHz)  */
/*      outData  Output OTF object                                        */
/*      inscan   root of input file names                                 */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  gchar FullFile[128];
  olong irow;
  olong disk, ver, nrow;
  gchar *tab;
  ObitTableGBTANTPOSGR* table;
  ObitTableGBTANTPOSGRRow* row;
  ObitIOCode retCode;

  /* error checks */
  if (err->error) return;
  g_assert(inscan!=NULL);
  g_assert(myInput!=NULL);

  /* get full file name */
  sprintf (FullFile,"Antenna/%s.fits", inscan);
  
  /* Create table structure */
  table = newObitTableGBTANTPOSGR("Antenna");
  if (err->error) return;

  /* Setup */
  disk = 2;  /* Input data directory */
  tab = "ANTPOSGR";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(table,disk,FullFile,tab,ver,nrow,err);

  /* Open */
  retCode = ObitTableGBTANTPOSGROpen (table, OBIT_IO_ReadOnly, err);
  if (err->error) return;

  /* Create Row structure */
  row = newObitTableGBTANTPOSGRRow (table);

  /* make sure there is data */
  if (table->myDesc->nrow<=0) {
     Obit_log_error(err, OBIT_Error, "No data in Antenna file for scan %s", inscan);
     return;
 }

  /* Create arrays */
  nAntTime = table->myDesc->nrow;
  AntDMJD  = g_malloc0(nAntTime*sizeof(odouble));
  AntAz    = g_malloc0(nAntTime*sizeof(odouble));
  AntEl    = g_malloc0(nAntTime*sizeof(odouble));

  /* Loop over table */
  for (irow = 1; irow<=nAntTime; irow++) {
    retCode = ObitTableGBTANTPOSGRReadRow (table, irow, row, err);
    if (err->error) return;
    AntDMJD[irow-1] = row->dmjd;
    AntAz[irow-1]   = row->mntAaz;
    AntEl[irow-1]   = row->mntEl;
 } /* end loop over table */
  
  /* Close */
  retCode = ObitTableGBTANTPOSGRClose (table, err);
  if (err->error) return;

  /* Cleanup */
  table = ObitTableGBTANTPOSGRUnref(table);
  row = ObitTableGBTANTPOSGRRowUnref(row);
} /* end GetAntennaGR  */

void GetQuadData(gchar *inscan, ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get Quadrant detector data, leave in globals                          */
/*      outData  Output OTF object                                        */
/*      inscan   root of input file names                                 */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  gchar FullFile[128];
  olong irow;
  olong disk, ver, nrow;
  gchar *tab;
  ObitTableGBTQUADDETECTOR* table;
  ObitTableGBTQUADDETECTORRow* row;
  ObitIOCode retCode;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(inscan!=NULL);
  g_assert(myInput!=NULL);

  /* get full file name */
  sprintf (FullFile,"QuadrantDetector-QuadrantDetector-QuadrantDetectorData/%s.fits", inscan);
  
  /* Create table structure */
  table = newObitTableGBTQUADDETECTOR("Antenna");
  if (err->error) return;

  /* Setup */
  disk = 2;  /* Input data directory */
  tab  = "QuadrantDetectorData"; 
  ver  = 1; 
  nrow = 1;
  ObitTableSetFITS(table,disk,FullFile,tab,ver,nrow,err);

  /* Open */
  retCode = ObitTableGBTQUADDETECTOROpen (table, OBIT_IO_ReadOnly, err);
  if (err->error) return;

  /* Create Row structure */
  row = newObitTableGBTQUADDETECTORRow (table);

  /* make sure there is data */
  if (table->myDesc->nrow<=0) {
     Obit_log_error(err, OBIT_Error, "No data in Antenna file for scan %s", inscan);
     return;
 }

  /* Create arrays */
  nQDetTime = table->myDesc->nrow;
  QDetDMJD = g_malloc0(nQDetTime*sizeof(odouble));
  QDetV1   = g_malloc0(nQDetTime*sizeof(ofloat));
  QDetV3   = g_malloc0(nQDetTime*sizeof(ofloat));
  QDetV4   = g_malloc0(nQDetTime*sizeof(ofloat));
  QDetV5   = g_malloc0(nQDetTime*sizeof(ofloat));
  QDetdAz  = g_malloc0(nQDetTime*sizeof(ofloat));
  QDetdEl  = g_malloc0(nQDetTime*sizeof(ofloat));

  /* Loop over table */
  for (irow = 1; irow<=nQDetTime; irow++) {
    retCode = ObitTableGBTQUADDETECTORReadRow (table, irow, row, err);
    if (err->error) return;
    QDetDMJD[irow-1]= row->dmjd;
    QDetV1[irow-1]  = row->ch1Voltage;
    QDetV3[irow-1]  = row->ch3Voltage;
    QDetV4[irow-1]  = row->ch4Voltage;
    QDetV5[irow-1]  = row->ch5Voltage;
    QDetdAz[irow-1] = row->X_Axis;
    QDetdEl[irow-1] = row->Z_Axis;
 } /* end loop over table */
  
  /* Close */
  retCode = ObitTableGBTQUADDETECTORClose (table, err);
  if (err->error) return;

  /* Cleanup */
  table = ObitTableGBTQUADDETECTORUnref(table);
  row = ObitTableGBTQUADDETECTORRowUnref(row);
} /* end GetQuadData  */

gboolean GetHeader (ObitOTF *outData, gchar *outOTF, gchar *inscan, 
		    ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get header information from scan header files                         */
/*  Returns TRUE if the file is just created                              */
/*   Input:                                                               */
/*      outData  Output OTF object                                        */
/*      outOTF   Output file name                                         */
/*      inscan   root of input file names                                 */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  gboolean out = FALSE;
  ObitOTFDesc *desc;
  ObitOTFArrayGeom *geom;
  ObitTableGBTPARSENSOR *sensorTable=NULL;
  olong numberDetect, ncol, ncopy, ver, indx;
  gchar *tab;
  long iscan;
  gchar Date[48], *fullname=NULL, FullFile[128], commnt[81];
  odouble SiteLat, SiteLong, SiteElev, JD, T, GMST0, GSTIAT0, nu, e2;
  odouble aEarth = 6378137.0; /* Semi major axis of Earth */
  odouble flat   = 1.0 / 298.257223563; /* inverse of flattening of Earth */ 
  ofloat UT1UTC, PolarX, PolarY, adcsampt;
  long smpperst, intgdp, aisxl;
  olong i, disk, nrow, ierr = 0;
  fitsfile *fptr;
  double dtemp;
  float ftemp;
  olong tcol[NUMDETECT], trow[NUMDETECT];
  gchar *routine = "GetHeader";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitOTFIsA(outData));
  g_assert(outOTF!=NULL);
  g_assert(inscan!=NULL);
  g_assert(myInput!=NULL);

  desc = outData->myDesc;

  /* Fixed frequency setup */
  nchan        = 1;
  deltaFreq    = 8.0e9;
  refPixel     = 1.0;
  refFrequency = 90.0e9;

  /* Sensor order table */
  /* get full file name */
  sprintf (FullFile,"Rcvr_PAR/%s.fits", inscan);
  
  /* Create table structure */
  sensorTable = newObitTableGBTPARSENSOR("Data");
  if (err->error) return out;

  /* Setup */
  disk = 2;
  tab  = "Sensor";
  ver  = 1; 
  nrow = 1;
  ObitTableSetFITS(sensorTable,disk,FullFile,tab,ver,nrow,err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, out);

  numberDetect =  NUMDETECT;
  for (i=0; i<numberDetect; i++) {
    tcol[i] = (olong)beamOffset[i*4+2];
    trow[i] = (olong)beamOffset[i*4+3];
  }
  xlate = ObitTableGBTPARSENSORXlate (sensorTable, numberDetect, trow, tcol, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, out);
  sensorTable = ObitTableGBTPARSENSORUnref(sensorTable);

  /* Set Feed Array Geometry 
     create Array geometry with enough elements */
  numberDetect =  NUMDETECT;
  geom = ObitOTFArrayGeomCreate (numberDetect);
  geom->azOffset = g_realloc(geom->azOffset, numberDetect*sizeof(ofloat));
  geom->elOffset = g_realloc(geom->elOffset, numberDetect*sizeof(ofloat));

  /* Create detector table - each output data stream has one  */
  for (i=0; i<numberDetect; i++) {
    /* Index of data into beamOffset */
    indx = xlate[i];
    if ((indx>=0) && (beamOffset[indx*4]>-9000.0)) {
      geom->azOffset[i] =  beamOffset[indx*4]   / 3600.0;
      /* NOTE: sign flip for Brian's measurements*/
      geom->elOffset[i] = -beamOffset[indx*4+1] / 3600.0;
    } else {  /* No valid data */
      geom->azOffset[i] =  0.0;
      geom->elOffset[i] =  0.0;
    }
  }

  /* rerefer all the offsets to feed 1 */
  /* Don't do for PAR
  azOff = -geom->azOffset[0];
  elOff = -geom->elOffset[0];
  for (i=0; i<numberDetect; i++) {
    geom->azOffset[i] += azOff;
    geom->elOffset[i] += elOff;
  }
  */
  
  /* Source info from GO file 
     get full file name */
  sprintf (FullFile,"%s/GO/%s.fits", DataRoot, inscan);

  /* Get header keywords from naked cfitsio */
   if ( fits_open_file(&fptr, FullFile, READONLY, &ierr) ) {
      Obit_log_error(err, OBIT_Error, "ERROR %d opening input FITS file %s", 
		     ierr, FullFile);
      return out;
   }

   /* Read main file keywords */
   iscan = 1;
   fits_read_key_lng (fptr, "SCAN", &iscan, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   /* Target name */
   fits_read_key_str (fptr, "OBJECT", Name, commnt, &ierr);
   /* Make sure at least 16 characters */
   for (i=0; i<16; i++) if (Name[i]==0) Name[i]=' ';

   /* Target RA */
   targRA = 0.0;
   fits_read_key_dbl (fptr, "RA      ", &targRA, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   /* Now in degrees targRA *= 15.0;   to degrees */

   /* Target Dec */
   targDec = 0.0;
   fits_read_key_dbl (fptr, "DEC     ", &targDec, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   /* Target Equinox */
   Equinox = 2000.0;
   fits_read_key_dbl (fptr, "EQUINOX ", &Equinox, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   /* Close FITS file */
   fits_close_file (fptr, &ierr);
   if (ierr!=0) {
     Obit_log_error(err, OBIT_Error, "ERROR reading input FITS file %s", FullFile);
     return out;
   }

   /* Other information - Get from Antenna file */
   /* get full file name */
   sprintf (FullFile,"%s/Antenna/%s.fits", DataRoot, inscan);
   
   /* Get header keywords from naked cfitsio */
   if ( fits_open_file(&fptr, FullFile, READONLY, &ierr) ) {
     Obit_log_error(err, OBIT_Error, "ERROR %d opening input FITS file %s", 
		    ierr, FullFile);
     return out;
   }

   /* Read main file keywords */
   ftemp = 0.0;
   fits_read_key_flt (fptr, "DELTAUTC", &ftemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   UT1UTC = (ofloat)ftemp;

   ftemp = 0.0;
   fits_read_key_flt (fptr, "IERSPMX", &ftemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   PolarX = (ofloat)ftemp;

   ftemp = 0.0;
   fits_read_key_flt (fptr, "IERSPMY", &ftemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   PolarY = (ofloat)ftemp;

   dtemp = 0.0;
   fits_read_key_dbl (fptr, "SITELAT", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   SiteLat = (odouble)dtemp;

   dtemp = 0.0;
   fits_read_key_dbl (fptr, "SITELONG", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   SiteLong = (odouble)dtemp;
   
   dtemp = 0.0;
   fits_read_key_dbl (fptr, "SITEELEV", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   SiteElev = (odouble)dtemp;

   fits_read_key_str (fptr, "DATE-OBS", Date, commnt, &ierr);
 
   /* Close FITS file */
   fits_close_file (fptr, &ierr);
   if (ierr!=0) {
     Obit_log_error(err, OBIT_Error, "ERROR reading input FITS file %s", FullFile);
     return out;
   }

   /* Convert to JD */
  ObitOTFDescDate2JD (Date, &JD);

  /* GST at IAT=0 at 0h on reference date (deg)
     Tropical century from jan 0.5, 2000 */
  T = (JD - 2451545.0) / 36525.0;

  /* GMST at IAT=0 in radians */
  GMST0 = ((((((-6.2e-6 * T) + 0.093104) * T) + 8640184.812866) * T + 24110.54841) 
	   * 2.0 * G_PI / 86400.0);
  /* to degrees */
  GSTIAT0 = RAD2DG * fmod(GMST0, (2.0*G_PI));
  
  ncopy = strlen (Date);
  if (ncopy>10) ncopy = 10;
  for (i=0; i<ncopy; i++) geom->RefDate[i] = Date[i]; geom->RefDate[i]=0;
  geom->TimeSys[0] = 'U'; geom->TimeSys[1] = 'T';
  geom->TimeSys[2] = 'C'; geom->TimeSys[3] = 0;

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

  /* Compute some useful terms 
     telescope latitude in radians */
  geom->lat = SiteLat * DG2RAD;
  /* telescope East longitude in radians */
  geom->lon = -SiteLong * DG2RAD;
  /* LST at iat0 in radians */
  geom->LSTiat0 = geom->GSTiat0*1.74533e-2 + geom->lon;
  /* Earth rotation rate in rad/day */
  geom->RadDay = geom->DegDay*1.74533e-2;
  /* Data - IAT in days */
  geom->dataIat = (geom->dataUtc - geom->iatUtc) / 86400.0;

  /* Attach Array geometry to OTF */
  outData->geom = ObitOTFArrayGeomRef(geom);
  geom = ObitOTFArrayGeomUnref(geom);
  
  /* Other information from PAR file 
     get full file name */
  sprintf (FullFile,"%s/Rcvr_PAR/%s.fits", DataRoot, inscan);
  
  /* Get header keywords from naked cfitsio */
   if ( fits_open_file(&fptr, FullFile, READONLY, &ierr) ) {
      Obit_log_error(err, OBIT_Error, "ERROR %d opening input FITS file %s", 
		     ierr, FullFile);
      return out;
   }

   /* Read main file keywords */
   /* PAR defined to have 64 feeds x 1 poln + 8 "dark" detectors */

   /* Have to work out integration time */
   adcsampt = 1.0e-07;
   fits_read_key_flt (fptr, "ADCSAMPT", &adcsampt, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   smpperst = 1250;
   fits_read_key_lng (fptr, "SMPPERST", &smpperst, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   intgdp   = 10;
   fits_read_key_lng (fptr, "INTGPD", &intgdp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   /* Cal flag */
   aisxl    = 0;
   fits_read_key_lng (fptr, "AISXL", &aisxl, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   /* Integration time in seconds */
   integTime = 0.0;
   /* Convert from seconds to days */
   integTime /= 86400.0;
   
   /* Close FITS file */
   fits_close_file (fptr, &ierr);
   if (ierr!=0) {
     Obit_log_error(err, OBIT_Error, "ERROR reading input FITS file %s", FullFile);
     return out;
   }

  /* initialize globals */
  refMJD = 0.0;
  target = 0.0;
  scan   = (ofloat)iscan;

   
  /* Does the input file exist?  BTW, I only do FITS */
  disk = 1;  /* input "disk" */
  fullname = ObitFITSFilename (disk, outOTF, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, out);


  /* check if output file exists, if not, initialize */
  if (!ObitFileExist (fullname, err)) {
    if (err->error>0) Obit_log_error(err, OBIT_Error, "ERROR testing file %s", fullname);

    /* Initialize  Descriptor */
    /* Approximate Beam size */
    desc->beamSize = 2.0e8 / refFrequency;

    /* diameter */
    desc->diameter = 90.0;  /* GBT diameter illuminated */

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
    desc->colRepeat[ncol] = numberDetect*2;  /* Data + weight */
    ncol++;
    
    desc->ncol = ncol;
    
    /* Data array descriptors */
    desc->naxis = 0;

    /* Data-Wt axis */
    desc->inaxes[desc->naxis] = 2;
    strncpy (desc->ctype[desc->naxis], "DATAWT", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = 1.0;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = 1.0;
    desc->naxis++;

    /* Feed axis */
    desc->inaxes[desc->naxis] = nfeed;
    strncpy (desc->ctype[desc->naxis], "FEED", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = 1.0;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = 1.0;
    desc->naxis++;

    /* Stokes axis */
    desc->inaxes[desc->naxis] = 1;
    strncpy (desc->ctype[desc->naxis], "STOKES  ", OTFLEN_KEYWORD);
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    /* Set appropritate Stokes Type (I) */
    desc->cdelt[desc->naxis] = 1.0;
    desc->crval[desc->naxis] = 1.0;
    desc->naxis++;

    /* Frequency axis */
    desc->inaxes[desc->naxis] = 1;
    strncpy (desc->ctype[desc->naxis], "FREQ    ", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = deltaFreq;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = refFrequency;
    desc->naxis++;

    /* Reference Modified Julian date */
    refMJD = JD - 2400000.5;

    /* Mark as PAR data */
    desc->OTFType = OBIT_GBTOTF_PAR;
    /* end initialize descriptor */
  }

  /* cleanup */
  if (fullname) g_free(fullname);

  /* Index the descriptor */
  ObitOTFDescIndex (desc);

  return out;
} /* end GetHeader */

void CnvrtData (ObitOTF *outData, gchar *inscan, ObitInfoList *myInput, 
		ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert data and write as OTFSoln table                               */
/*      outData  Output OTF object                                        */
/*      inscan   Scan part of input file name                             */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  olong irow, iRow;
  olong i, j, ndetect, npoly, ntime, ver, solnVer;
  ofloat refDay, corDay;
  odouble JD;
  ObitOTFDesc *desc;
  ObitTableOTFSoln *outSoln=NULL;
  ObitTableOTFSolnRow *row=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar *tname;
  ofloat fblank = ObitMagicF();  
  gchar *routine = "CnvrtData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(inscan!=NULL);
  g_assert(myInput!=NULL);

  /* Consistency checks */
  desc = outData->myDesc;
  Obit_return_if_fail((nfeed==desc->inaxes[desc->jlocfeed]), err,
		       "%s: Request number feeds incompatible %d !=  %d", 
		      routine, nfeed, desc->inaxes[desc->jlocfeed]);
  Obit_return_if_fail((nstok==desc->inaxes[desc->jlocs]), err,
		       "%s: Request number Poln incompatible %d !=  %d", 
		      routine, nstok, desc->inaxes[desc->jlocs]);
  Obit_return_if_fail((nchan==desc->inaxes[desc->jlocf]), err,
		       "%s: Request number freq. incompatible %d !=  %d", 
		      routine, nchan, desc->inaxes[desc->jlocf]);
  
  /* Get reference MJD, convert ref date  to JD */
  ObitOTFDescDate2JD (outData->geom->RefDate, &JD);
  /* Reference Modified Julian date */
  refMJD = JD - 2400000.5;

  /* Convert measurements to offsets */
  ProcessQDData (err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Output OTFSoln table version */
  solnVer = 1;
  ObitInfoListGetTest(myInput, "solnVer",  &type, dim, &solnVer);

  /* Create output OTFSoln table */
  tname = g_strconcat ("Calibration for: ",outData->name, NULL);
  ver     = solnVer;
  ndetect = 1;  /* only dummy */
  npoly = 1;
  outSoln = newObitTableOTFSolnValue(tname, (ObitData*)outData, &ver, OBIT_IO_ReadWrite,  
				     ndetect, npoly, err);
  g_free (tname);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open output table */
  if ((ObitTableOTFSolnOpen (outSoln, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR opening input OTFSoln table", routine);
    return;
  }

  /* Create Row */
  row = newObitTableOTFSolnRow (outSoln);

  /* Attach row to output buffer */
  ObitTableOTFSolnSetRow (outSoln, row, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Initialize solution row */
  row->dAz = 0.0;
  row->dEl = 0.0;
  row->Target = 0;
  row->TimeI  = (QDetDMJD[9]-QDetDMJD[0])*0.1;
  for (j=0; j<ndetect; j++) row->mult[j] = 1.0;
  for (j=0; j<ndetect; j++) row->wt[j]   = 1.0;
  for (j=0; j<ndetect; j++) row->cal[j]  = 0.0;
  for (j=0; j<ndetect; j++) row->add[j]  = 0.0;
  for (i=0; i<npoly; i++)   row->poly[i] = 0.0;

  /* Write at end */
  startRec  = MAX(1,outSoln->myDesc->nrow);  /* First record in scan */
  startTime = -1.0e20;                       /* dummy start time */
  /* Initialize end times and record numbers in case of no data */
  endTime = startTime;
  endRec  = startRec; 

  /* Get difference in output reference day and first day in input */
  corDay = refMJD - refDay;
  
  /* Anything to do? */
  ntime = nQDetTime;
  if (ntime<=0) return;

  /* Loop over data */
  for (irow=0; irow<ntime; irow++) {
    if ((QDetdAz[irow]==fblank) || (QDetdEl[irow]==fblank)) continue;

    /* first time in scan */
    if (startTime<-1000.0) startTime = QDetDMJD[irow-1] - corDay;

    /* Get values from array */
    row->Time = QDetDMJD[irow] - corDay;
    row->dAz  = QDetdAz[irow] / 3600.0;  /* Deg */
    row->dEl  = QDetdEl[irow] / 3600.0;

    /* Write Soln table */
    iRow = -1;
    if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table", routine);
      return;
    }
	
  } /* end loop over table */  

  /* Get end times and record numbers */
  endTime = row->Time;
  endRec  = outSoln->myDesc->nrow;  /* Last record in scan */

  /* Close output cal table */
  if ((ObitTableOTFSolnClose (outSoln, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR closing output OTFSoln Table", routine);
    return;
  }
  
  /* Cleanup */
  row = ObitTableOTFSolnUnref(row);
  outSoln = ObitTableOTFSolnUnref(outSoln);

} /* end CnvrtData  */

void GetElev (odouble time, ofloat *el)
/*----------------------------------------------------------------------- */
/*  Get antenna elevation at a given time,                                */
/*  End points or linear interpolation                                    */
/*   Input:                                                               */
/*      time   The desired MJD                                            */
/*   Output:                                                              */
/*      el     Elevation in degrees                                       */
/*----------------------------------------------------------------------- */
{
  olong i, best;
  odouble test, delta;
  ofloat w1, w2;

  /* Initial value */
  *el  = 0.0;
  
  /* Find closest */
  best = -1;
  delta = 1.0e20;
  for (i=0; i<nAntTime; i++) { /* loop over array */
    test = fabs (AntDMJD[i]-time);
      if (delta> test) {
	delta = test;
	best = i;
      } else { /* must be getting further, stop */
	break;
      }
  }

  /* end points */
  if ((best==0) || (best==(nAntTime-1))) {
    /* Drop data outside of range - put at SCP */
    if ((time<AntDMJD[0]) || (time>AntDMJD[nAntTime-1])) {
      *el  = -90.0;
    } else {
      *el = AntEl[best];
    }
  } else if (time<AntDMJD[best]){ /* interpolate with previous */
    w1 = (AntDMJD[best]-time) / (AntDMJD[best]-AntDMJD[best-1]);
    w2 = 1.0 - w1;
    *el  = w1 * AntEl[best-1] + w2 * AntEl[best];
  } else { /* interpolate with following */
    w1 = (AntDMJD[best+1]-time) / (AntDMJD[best+1]-AntDMJD[best]);
    w2 = 1.0 - w1;
    *el  = w1 * AntEl[best] + w2 * AntEl[best+1];
  }
} /* end GetElev */

olong GetTarget (ObitOTF *outData, gboolean isNew, gchar *name, 
		odouble ra, odouble dec, odouble equinox, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get target id, look through existing table create new entry if needed.*/
/*  Returns Target id                                                     */
/*   Input:                                                               */
/*      outData  Output OTF object                                        */
/*      isNew    True if output file just created                         */
/*      name     Name of target                                           */
/*      ra       RA of target                                             */
/*      dec      Dec of target                                            */
/*      equinox  equinox of ra, dec                                       */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*   Return:                                                              */
/*      Target id                                                         */
/*----------------------------------------------------------------------- */
{
  olong targ = -1;
  ObitTableOTFTarget* table;
  ObitTableOTFTargetRow* row;
  olong iRow, ver;
  gboolean doWrite;
  ObitIOAccess access;
  gchar tName[20];
  gchar *routine = "GetTarget";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return targ;
  g_assert (ObitOTFIsA(outData));
  g_assert(name!=NULL);

  ObitTrimTrail(name);  /* Trim trailing blanks */

  /* create Target table object */
  ver = 1;
  if (isNew) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  table = newObitTableOTFTargetValue ("Target table", (ObitData*)outData, &ver, access, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, targ);

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
  if (err->error) Obit_traceback_val (err, routine, outData->name, targ);

  /* Newly created?  Just write new one */
  doWrite = FALSE;
  if (isNew) {
    targ = 1;
    row->TargID = targ;
    strncpy(row->Target, name, 16);
    row->RAMean  = ra;
    row->DecMean = dec;
    row->Epoch   = equinox;
    doWrite = TRUE;
  } else { /* Existing, see if already exists? */

    /* loop through table */
    for (iRow = 1; iRow<=table->myDesc->nrow; iRow++) {
      if ((ObitTableOTFTargetReadRow (table, iRow, row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading OTFTarget Table file");
	return targ;
      }
      strncpy (tName, row->Target, 16); tName[16]=0;
      ObitTrimTrail(tName);
      if (!strncmp (tName, name, 16)) {
	/* Found match */
	targ = row->TargID;
	break;
      }  
    } /* end loop over table */

    /* Add new entry? */
    if (targ<=0) {
      targ = table->myDesc->nrow + 1;
      row->TargID = targ;
      row->RAMean  = ra;
      row->DecMean = dec;
      row->Epoch   = equinox;
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

void InitScan (ObitOTF *outData, gboolean isNew, ofloat scan, ofloat target,
	       ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Initializes Index table for this scan creating an entry               */
/*   Input:                                                               */
/*      outData  Output OTF object                                        */
/*      isNew    True if output file just created                         */
/*      scan     Scan number                                              */
/*      target   Target ID number                                         */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*   Output global:                                                       */
/*       lastTime  End time or prior scan                                 */
/*----------------------------------------------------------------------- */
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
  g_assert (ObitOTFIsA(outData));

  /* create Index table object */
  scanID = (olong)(scan+0.5);
  targetID = (olong)(target+0.5);
  ver = 1;
  if (isNew) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  table = newObitTableOTFIndexValue ("Index table", (ObitData*)outData, &ver, access, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

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
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

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

void SetScan (ObitOTF *outData, odouble startTime, odouble endTime, 
	      olong startRec, olong endRec, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Updates Index table for this scan creating an entry                   */
/*   Input:                                                               */
/*      outData   Output OTF object                                       */
/*      isNew     True if output file just created                        */
/*      startTime Start time of scan in days                              */
/*      endTime   End time of scan in days                                */
/*      startRec  First record in scan                                    */
/*      endRec    Last record in scan                                     */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitTableOTFIndex* table;
  ObitTableOTFIndexRow* row;
  olong iRow, ver;
  ObitIOAccess access;
  gchar *routine = "SetScan";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(outData));

  /* create Index table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  table = newObitTableOTFIndexValue ("Index table", (ObitData*)outData, &ver, access, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

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
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Read next to the last record if one and times not specified */
  if ((table->myDesc->nrow>0) && (endTime<=0.0)) {
    iRow = table->myDesc->nrow-1;
    if ((ObitTableOTFIndexReadRow (table, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR reading OTFIndex Table file");
      return;
    }
    startTime = row->Time + 0.51*row->TimeI;
    endTime   = row->Time + 0.52*row->TimeI;
    startRec  = row->EndRec;
    endRec    = row->EndRec;
  }

  /* Update last record */
  iRow = table->myDesc->nrow;
  if ((ObitTableOTFIndexReadRow (table, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR reading OTFIndex Table file");
    return;
  }

  /* upate row */
  row->Time = 0.5 * (startTime + endTime);
  row->TimeI = (endTime - startTime);
  row->StartRec = startRec;
  row->EndRec   = endRec;

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

void ProcessQDData (ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert measurements to offsets                                       */
/*  Based on GBT Memo PTCS/PN/64.3 (Ries, Hunter and Ghigo) and IDL       */
/*  routines by Paul Ries                                                 */
/*   Input in global                                                      */
/*     nQDetTime  Number of Quadrant detector time samples                */
/*     QDetDMJD   Array of  Quadrant detector times                       */
/*     QDetV1     Array of  Quadrant detector voltage 1                   */
/*     QDetV3     Array of  Quadrant detector voltage 3                   */
/*     QDetV4     Array of  Quadrant detector voltage 4                   */
/*     QDetV5     Array of  Quadrant detector voltage 5                   */
/*   Output in global:                                                    */
/*     QDetdAz    Array of  Quadrant detector Az offset, in arcseconds    */
/*     QDetdEl    Array of  Quadrant detector El offset, in arcseconds    */
/*     err        Obit return error stack                                 */
/*----------------------------------------------------------------------- */
{
  ofloat *el=NULL, median, med1, med3, med4, med5;
  odouble sum, *C1Pos=NULL, *C2Pos=NULL, **g=NULL; 
  odouble **C_comps=NULL, **C_err=NULL; 
  /* C matrix comes from August 8, 2008 */
  odouble ESM_C[2][3]={{-1.01742164436251, 0.341134251429814, 0.390820531900982},
		       {-0.240294958082037, 1.09469317346167, -0.822147144704265}};
  odouble ESM_D[2][2] =  {{-549.9392, -156.4483},{-570.0, 468.0}};
  olong i, j, k, n;
  ofloat fblank = ObitMagicF();  
  gchar *routine = "ProcessQDData";

  /* error checks */
  if (err->error) return;

  /* Allocate arrays */
  el    = g_malloc0(nQDetTime*sizeof(ofloat));
  C1Pos = g_malloc0(nQDetTime*sizeof(odouble));
  C2Pos = g_malloc0(nQDetTime*sizeof(odouble));
  g     = g_malloc0(3*sizeof(odouble*));
  g[0]  = g_malloc0(nQDetTime*sizeof(odouble));
  g[1]  = g_malloc0(nQDetTime*sizeof(odouble));
  g[2]  = g_malloc0(nQDetTime*sizeof(odouble));
  C_comps    = g_malloc0(2*sizeof(odouble*));
  C_comps[0] = g_malloc0(nQDetTime*sizeof(odouble));
  C_comps[1] = g_malloc0(nQDetTime*sizeof(odouble));
  C_err      = g_malloc0(2*sizeof(odouble*));
  C_err[0]   = g_malloc0(nQDetTime*sizeof(odouble));
  C_err[1]   = g_malloc0(nQDetTime*sizeof(odouble));

  /* Get elevations, -90 = no value */
  for (i=0; i<nQDetTime; i++) GetElev (QDetDMJD[i], &el[i]);

  /* QD position */
  for (i=0; i<nQDetTime; i++) {
    C1Pos[i] = (QDetV1[i] - QDetV3[i]) / (QDetV1[i] + QDetV3[i]);
    C2Pos[i] = (QDetV4[i] - QDetV5[i]) / (QDetV4[i] + QDetV5[i]);
  }

  /* Fill g matrix */
  for (i=0; i<nQDetTime; i++) {
    if (el[i]>0.0) {
      g[0][i] = cos(el[i]*DG2RAD);
      g[1][i] = sin(el[i]*DG2RAD);
      g[2][i] = 1.0;
    } else {  /* No valid el */
      g[0][i] = sqrt(2.0);
      g[1][i] = sqrt(2.0);
      g[2][i] = 1.0;
    }
  }

  /* determine structure model
     C_comps = ESM_C*g */
  n = nQDetTime;
  for (j=0; j<2; j++) {
    for (i=0; i<n; i++) {
      sum = 0.0;
      for (k=0; k<3; k++) {
	sum += ESM_C[j][k] * g[k][i];
      }
      C_comps[j][i] = sum;
    }
  }
  
  /* subtract engineering model  
     C_err=[[C1Posi],[C2Posi]]-C_comps */
  for (i=0; i<n; i++) {
    C_err[0][i] = C1Pos[i] -  C_comps[0][i];
    C_err[1][i] = C2Pos[i] -  C_comps[1][i];
  }

  /*  convert from QD units to pointing error
      qd_xz=ESM_D##C_err (use QDetdAz, QDetdEl )*/
  for (i=0; i<n; i++) {
    QDetdEl[i] = +(ofloat)(ESM_D[0][1] * C_err[0][i] + ESM_D[1][1] * C_err[1][i]);
    QDetdAz[i] = -(ofloat)(ESM_D[0][0] * C_err[0][i] + ESM_D[1][0] * C_err[1][i]); 
  }

  /* Blank values with no el */
  for (i=0; i<n; i++) {
    if (el[i]<0.0) {
      QDetdAz[i] = fblank;
      QDetdEl[i] = fblank;
    }
  }


  /* Filter - remove median */
  median = medianValue(QDetdAz,1,n);
  for (i=0; i<n; i++) 
    if (QDetdAz[i]!=fblank) QDetdAz[i] -= median;
  median = medianValue(QDetdEl,1,n);
  for (i=0; i<n; i++) 
    if (QDetdEl[i]!=fblank) QDetdEl[i] -= median;

  /* Blank values with no el - again, medianValue unblanks them */
  for (i=0; i<n; i++) {
    if (el[i]<0.0) {
      QDetdAz[i] = fblank;
      QDetdEl[i] = fblank;
    }
  }
  /* Data quality checks
     fog warning */
  med1 = medianValue(QDetV1,1,n);
  med3 = medianValue(QDetV3,1,n);
  med4 = medianValue(QDetV4,1,n);
  med5 = medianValue(QDetV5,1,n);
  if (((med1+med3+med4+med5)/4.0) < 0.08) 
    Obit_log_error(err, OBIT_InfoWarn, 
		   "%s: Low QD voltages, possible fog. Dubious QD correction", 
		   routine);

  /* possible hysteresis warning */
  med1 = medianValue(QDetdAz,1,n/4);
  med3 = medianValue(&QDetdAz[n-n/4],1,n/4);
  if (fabs(med1-med3) > (10.0/3600.0)) 
    Obit_log_error(err, OBIT_InfoWarn, 
		   "%s: Possible hysteresis event in quadrant detector", 
		   routine);
  
  /* Cleanup */
  if (el)    g_free(el);
  if (C1Pos) g_free(C1Pos);
  if (C2Pos) g_free(C2Pos);
  if (g) {
    if(g[0]) g_free(g[0]);
    if(g[1]) g_free(g[1]);
    if(g[2]) g_free(g[2]);
    g_free(g);
  }
  if (C_comps) {
    if(C_comps[0]) g_free(C_comps[0]);
    if(C_comps[1]) g_free(C_comps[1]);
    g_free(C_comps);
  }
  if (C_err) {
    if(C_err[0]) g_free(C_err[0]);
    if(C_err[1]) g_free(C_err[1]);
    g_free(C_err);
  }
} /* end ProcessQDData  */

void GetPointOffset (odouble time, ofloat *dAz, ofloat *dEl)
/*----------------------------------------------------------------------- */
/*  Get antenna pointing error at a given time form quadrant detector     */
/*  End points or linear interpolation                                    */
/*   Input:                                                               */
/*      time   The desired MJD                                            */
/*   Output:                                                              */
/*      dAz    Azimuth in degrees                                         */
/*      dEl    Elevation in degrees                                       */
/*----------------------------------------------------------------------- */
{
  olong i, best;
  odouble test, delta;
  ofloat w1, w2;

  /* Initial values */
  *dAz  = 0.0;
  *dEl  = 0.0;
  
  /* Find closest */
  best = -1;
  delta = 1.0e20;
  for (i=0; i<nQDetTime; i++) { /* loop over array */
    test = fabs (QDetDMJD[i]-time);
      if (delta> test) {
	delta = test;
	best = i;
      } else { /* must be getting further, stop */
	break;
      }
  }

  /* end points */
  if ((best==0) || (best==(nQDetTime-1))) {
    /* Drop data outside of range - put at SCP */
    if ((time<QDetDMJD[0]) || (time>QDetDMJD[nQDetTime-1])) {
      *dAz  = 0.0;
      *dEl  = -90.0;
    } else {
      *dAz = QDetdAz[best];
      *dEl = QDetdEl[best];
    }
  } else if (time<QDetDMJD[best]){ /* interpolate with previous */
    w1 = (QDetDMJD[best]-time) / (QDetDMJD[best]-QDetDMJD[best-1]);
    w2 = 1.0 - w1;
    *dAz  = w1 * QDetdAz[best-1] + w2 * QDetdAz[best];
    *dEl  = w1 * QDetdEl[best-1] + w2 * QDetdEl[best];
  } else { /* interpolate with following */
    w1 = (QDetDMJD[best+1]-time) / (QDetDMJD[best+1]-QDetDMJD[best]);
    w2 = 1.0 - w1;
    *dAz  = w1 * QDetdAz[best] + w2 * QDetdAz[best+1];
    *dEl  = w1 * QDetdEl[best] + w2 * QDetdEl[best+1];
  }
} /* end GetPointOffset */

