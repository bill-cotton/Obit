/* $Id$  */
/* Read Penn Array (Mustang) data                                */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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
#include "ObitGBTBeamOffInfo.h"
#include "ObitTableGBTPARDATA.h"
#include "ObitTableGBTPARSENSOR.h"
#define NUMDETECT 72   /* Maximum number of detectors */
/* internal prototypes */
/* Get inputs */
ObitInfoList* PAROTFin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Get pointing times from Antenna file, Secondary focus */
void GetAntennaGR (gchar *infile, ObitInfoList *myInput, ObitErr *err);
/* Get pointing for a given time */
void GetPoint (odouble time, ofloat *ra, ofloat *dec);
/* Get file descriptor */
gboolean GetHeader (ObitOTF *outData, char *outfile, gchar *infile, 
		    ObitInfoList *myInput, ObitErr *err);
/* Get data */
void GetData (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, 
	      ofloat avgTime, ofloat offTime, olong scanNo, ObitErr *err);
/* Get target id */
olong GetTarget (ObitOTF *outData, gboolean isNew, gchar *name, 
		odouble ra, odouble dec, odouble equinox, ObitErr *err);
/* Initialize Index table */
void InitScan (ObitOTF *outData, gboolean isNew, ofloat scan, ofloat target, 
	       ObitErr *err);
/* Undate scan info */
void SetScan (ObitOTF *outData, odouble startTime, odouble endTime, 
	      olong startRec, olong endRec, ObitErr *err);
/* Read, filter, average data */
void ProcessData (gchar *inscan, ofloat avgTime,
		  olong *ndetect, olong *ntime, ofloat *refDate,
		  gfloat** ATime, gfloat*** AData, gfloat** ACal, 
		  ObitErr *err);
/* Read data for a single detector */
void ReadDatum (ObitTableGBTPARDATA* PARtable, olong ndetect, gboolean *bad, 
		gint *ntime, gdouble** time, gfloat** data, gfloat** cal, 
		ObitErr *err);
/* Average data a single detector */
void AverageDatum (gboolean bad, olong detector, olong nraw, odouble *time, 
		   ofloat *data, ofloat *cal, ofloat avgTime, olong *ntime, 
		   ofloat *refDate, ofloat **ATime, ofloat **AData, ofloat **ACal);
/* Remove baseline jumps for a single detector */
void deJumpDatum (gint ndata, ofloat *data, ofloat *cal);
/* IIR Filter for 1.4 Hz refrigerator hum, 10 Hz sampling */
void IIRFilter10 (gint ndata, ofloat *data);
/* IIR Filter for 1.4 Hz refrigerator hum, 20 Hz sampling */
void IIRFilter20 (gint ndata, ofloat *data);
/* Fit a given frequency to a set of time series */
void FitFreq (ofloat freq, olong ntime, ofloat *Atime, 
	      olong ndetect, ofloat **AData, gboolean *bad,  
	      ObitErr *err);
/* Average phase over detector */
ofloat AveragePhase (olong ndetect, ofloat *amp, ofloat *phase, gboolean *bad);

/* Program globals */
gchar *pgmName = "PAROTF";       /* Program name */
gchar *infile  = "PAROTF.inp";   /* File with program inputs */
gchar *outfile = "PAROTF.out";   /* File to contain program outputs */
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
olong nAntTime;   /* number of antenna time samples */
odouble *AntDMJD; /* Array of antenna times */
odouble *AntRA;   /* Array of Antenna RA J2000 values */
odouble *AntDec;  /* Array of Antenna RA J2000 values */
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
ObitGBTIFInfo*       IFdata     = NULL; /* IF information structure */
ObitGBTBeamOffInfo* BeamOffData = NULL; /* Beam offset information structure */
ofloat azOff, elOff; /* az and el offsets to apply to pointing positions */

/* Note: must flip sign of el offset */
/* Beam offset, (cross el,elevation),  asec -9999.0 => no data*/
ofloat tbeamOffset[] = {  
  /*               index row     col  comment */
  -9999.,  -9999., /*  00  0.00    8.00 Dark pixel */
  -9999.,  -9999., /*  01  1.00    8.00 Dark pixel */
  -9999.,  -9999., /*  02  2.00    8.00 Dark pixel */
  -9999.,  -9999., /*  03  3.00    8.00 Dark pixel */
  -9999.,  -9999., /*  04  4.00    8.00 Dark pixel */
  -9999.,  -9999., /*  05  5.00    8.00 Dark pixel */
  -9999.,  -9999., /*  06  6.00    8.00 Dark pixel */
  -9999.,  -9999., /*  07  7.00    8.00 Dark pixel */
  14.30,  -3.69, /*  08  0.00    0.00   */
  -0.74,  -2.86, /*  09  1.00    0.00   */
  -3.69, -14.30, /*  10  2.00    0.00   */
 -10.05,  -3.50, /*  11  3.00    0.00   */
 -20.01,   5.17, /*  12  4.00    0.00   */
   0.74,   2.86, /*  13  5.00    0.00   */
   3.69,  14.30, /*  14  6.00    0.00   */
   2.86,  -0.74, /*  15  7.00    0.00   */
  16.42,  -7.29, /*  16  0.00    1.00   */
   1.38,  -6.46, /*  17  1.00    1.00   */
  -7.29, -16.42, /*  18  2.00    1.00   */
 -13.65,  -5.63, /*  19  3.00    1.00   */
 -16.42,   7.29, /*  20  4.00    1.00   */
  -1.38,   6.46, /*  21  5.00    1.00   */
   7.29,  16.42, /*  22  6.00    1.00   */
   6.46,   1.38, /*  23  7.00    1.00   */
   8.58,  -2.21, /*  24  0.00    2.00   */
   5.63, -13.65, /*  25  1.00    2.00   */
  -2.21,  -8.58, /*  26  2.00    2.00   */
  -6.46,  -1.38, /*  27  3.00    2.00   */
 -10.70,   5.81, /*  28  4.00    2.00   */
  -5.63,  13.65, /*  29  5.00    2.00   */
   2.21,   8.58, /*  30  6.00    2.00   */
  13.65,   5.63, /*  31  7.00    2.00   */
  10.70,  -5.81, /*  32  0.00    3.00   */
   3.50, -10.05, /*  33  1.00    3.00   */
  -5.81, -10.70, /*  34  2.00    3.00   */
 -15.77,  -2.03, /*  35  3.00    3.00   */
 -12.82,   9.41, /*  36  4.00    3.00   */
  -3.50,  10.05, /*  37  5.00    3.00   */
   5.81,  10.70, /*  38  6.00    3.00   */
  10.05,   3.50, /*  39  7.00    3.00   */
  12.82,  -9.41, /*  40  0.00    4.00   */
   2.03, -15.77, /*  41  1.00    4.00   */
  -9.41, -12.82, /*  42  2.00    4.00   */
 -12.17,   0.09, /*  43  3.00    4.00   */
  -7.10,   7.93, /*  44  4.00    4.00   */
   0.09,  12.17, /*  45  5.00    4.00   */
   9.41,  12.82, /*  46  6.00    4.00   */
  15.77,   2.03, /*  47  7.00    4.00   */
   4.98,  -4.34, /*  48  0.00    5.00   */
  -0.09, -12.17, /*  49  1.00    5.00   */
  -4.34,  -4.98, /*  50  2.00    5.00   */
  -8.58,   2.21, /*  51  3.00    5.00   */
  -9.22,  11.53, /*  52  4.00    5.00   */
  -2.03,  15.77, /*  53  5.00    5.00   */
   4.34,   4.98, /*  54  6.00    5.00   */
  12.17,  -0.09, /*  55  7.00    5.00   */
   9.22, -11.53, /*  56  0.00    6.00   */
  -1.57, -17.89, /*  57  1.00    6.00   */
 -11.53,  -9.22, /*  58  2.00    6.00   */
 -17.89,   1.57, /*  59  3.00    6.00   */
  -4.98,   4.34, /*  60  4.00    6.00   */
   1.57,  17.89, /*  61  5.00    6.00   */
  11.53,   9.22, /*  62  6.00    6.00   */
  17.89,  -1.57, /*  63  7.00    6.00   */
   7.10,  -7.93, /*  64  0.00    7.00   */
  -5.17, -20.01, /*  65  1.00    7.00   */
  -7.93,  -7.10, /*  66  2.00    7.00   */
 -14.30,   3.69, /*  67  3.00    7.00   */
  -2.86,   0.74, /*  68  4.00    7.00   */
   5.17,  20.01, /*  69  5.00    7.00   */
   7.93,   7.10, /*  70  6.00    7.00   */
  20.01,  -5.17};/*  71  7.00    7.00   */

ofloat *beamOffset = tbeamOffset;
odouble BTime = 2.50e-6;  /* blank 220 msec after cal switch */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*    Read GBT Penn Array  data to an OTF dataset                         */
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
  ofloat avgTime, offTime;
  olong scanNo;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  err = newObitErr();

  /* Startup - parse command line */
  ierr = 0;
  myInput = PAROTFin (argc, argv, err);
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
  mySystem = ObitSystemStartup ("PAROTF", 1, 0, 0, NULL, nFITS, FITSdirs, 
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

  /* Averaging time */
  avgTime = 0.05;
  ObitInfoListGetTest(myInput, "avgTime",  &type, dim, &avgTime);
  if (avgTime<=0.0) avgTime = 0.05;
  avgTime /= 86400.0;  /* To days */

  /* Time offset */
  offTime = 0.0;
  ObitInfoListGetTest(myInput, "offTime",  &type, dim, &offTime);
  offTime /= 86400.0;  /* To days */

  /* Scan no. */
  scanNo = 0;
  ObitInfoListGetTest(myInput, "scanNo",  &type, dim, &scanNo);
  if (scanNo!=0) scan = (ofloat)scanNo;
  Obit_log_error(err, OBIT_InfoErr, "Adding scan %s %6.0f to OTF file %s", 
		 inscan, scan, outOTF);
  ObitErrLog(err);

  /* convert data  */
  GetData (outData, inscan, myInput, avgTime, offTime, scanNo, err);
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
   if (AntDMJD) g_free(AntDMJD);
   if (AntRA)   g_free(AntRA);
   if (AntDec)  g_free(AntDec);
   IFdata = ObitGBTIFInfoUnref(IFdata); 
   BeamOffData = ObitGBTBeamOffInfoUnref(BeamOffData); 
 
   return ierr;
} /* end of main */

ObitInfoList* PAROTFin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="PAROTF.in", *arg;
  gboolean init=FALSE;
  oint itemp;
  gchar *strTemp, config[128];
  ObitInfoList* list;
  ofloat *p;
  gchar *routine = "PAROTFIn";

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

    } else if (strcmp(arg, "-config") == 0){ /* configuration file */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "config", OBIT_string, dim, strTemp);

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

  /* Parse any configuration file */
  if (ObitInfoListGetTest(list, "config", &type, dim, config)) {
    if ((strncmp(config,"None",4)) && (strncmp(config,"    ",4)) && 
	(strlen(config)>=4)) {
      ObitTrimTrail(config);  /* Trim trailing blanks */
      ObitParserParse (config, list, err);
      if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

      /* Get values */
      if (ObitInfoListGetP (list, "beamOffset", &type, dim, (gpointer)&p)) 
	beamOffset = p;
      ObitInfoListGetTest (list, "BTime", &type, dim, &BTime);
    }
  } /* End parse configuration file */

  /* Extract basic information to program globals */
  ObitInfoListGet(list, "pgmNumber", &type, dim, &pgmNumber, err);
  /*ObitInfoListGet(list, "nFITS",     &type, dim, &nFITS,     err);*/
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  /* Initialize output */
  myOutput = defaultOutputs(err);
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

 return list;
} /* end PAROTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: PAROTF -input file [-Scan date/time]\n");
    fprintf(stderr, "Convert an GBT PAR file format to Obit/OTF\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def PAROTF.in\n");
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
  /*gfloat ftemp;*/
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
  g_assert(ObitErrIsA(err));
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
  AntDMJD = g_malloc0(nAntTime*sizeof(odouble));
  AntRA   = g_malloc0(nAntTime*sizeof(odouble));
  AntDec  = g_malloc0(nAntTime*sizeof(odouble));

  /* Loop over table */
  for (irow = 1; irow<=nAntTime; irow++) {
    retCode = ObitTableGBTANTPOSGRReadRow (table, irow, row, err);
    if (err->error) return;
    AntDMJD[irow-1] = row->dmjd;
    AntRA[irow-1]   = row->raj2000;
    AntDec[irow-1]  = row->decj2000;
 } /* end loop over table */
  
  /* Close */
  retCode = ObitTableGBTANTPOSGRClose (table, err);
  if (err->error) return;

  /* Cleanup */
  table = ObitTableGBTANTPOSGRUnref(table);
  row = ObitTableGBTANTPOSGRRowUnref(row);
} /* end GetAntennaGR  */

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
  olong numberDetect, ncol, ncopy;
  long iscan;
  gchar Date[48], *fullname=NULL, FullFile[128], commnt[81];
  odouble SiteLat, SiteLong, SiteElev, JD, T, GMST0, GSTIAT0, nu, e2;
  odouble aEarth = 6378137.0; /* Semi major axis of Earth */
  odouble flat   = 1.0 / 298.257223563; /* inverse of flattening of Earth */ 
  ofloat UT1UTC, PolarX, PolarY, adcsampt;
  long smpperst, intgdp, aisxl;
  olong i, disk, ierr = 0;
  fitsfile *fptr;
  double dtemp;
  float ftemp;
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

  /* Get Beam offset data
  BeamOffData = newObitGBTBeamOffInfoValue("BeamOff", disk, inscan, err); */

  /* Set Feed Array Geometry 
     create Array geometry with enough elements */
  numberDetect =  NUMDETECT;
  geom = ObitOTFArrayGeomCreate (numberDetect);
  geom->azOffset = g_realloc(geom->azOffset, numberDetect*sizeof(ofloat));
  geom->elOffset = g_realloc(geom->elOffset, numberDetect*sizeof(ofloat));

  /* Create detector table - each output data stream has one  */
  for (i=0; i<numberDetect; i++) {
    if (beamOffset[i*2]>-9000.0) {
      geom->azOffset[i] =  beamOffset[i*2]   / 3600.0;
      /* NOTE: sign flip for Brian's measurements*/
      geom->elOffset[i] = -beamOffset[i*2+1] / 3600.0;
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
  /* MORE WORK HERE */
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
   /* MORE WORK HERE */
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

void GetData (ObitOTF *outData, gchar *inscan, ObitInfoList *myInput, 
	      ofloat avgTime, ofloat offTime, olong scanNo, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Read data from GB FITS files and write                                */
/*      outData  Output OTF object                                        */
/*      inscan   Scan part of input file name                             */
/*      myInput  parser object                                            */
/*      avgTime  Data averaging time in days                              */
/*      offTime  Time offset for data in days                             */
/*      scanNo   Scan number, 0=> use GBT scan                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  olong irow;
  olong i, ntime, incdatawt;
  ofloat *ATime=NULL, **AData=NULL, *ACal=NULL;
  ofloat *data, ra, dec, refDay, corDay;
  ofloat fblank = ObitMagicF();
  odouble JD, dmjd=0.0;
  ObitOTFDesc *desc;
  ObitOTFArrayGeom *geom;
  gchar *routine = "GetData";

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

  /* Read, filter, average data from FITS file */
  ProcessData (inscan, avgTime, &ndetect, &ntime, &refDay, 
	       &ATime, &AData, &ACal, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  desc = outData->myDesc;
  geom = outData->geom;
  incdatawt = desc->incdatawt; /* increment in data-wt axis */
  
  /* Write at end */
  desc->firstRec = desc->nrecord+1;
  startRec = desc->firstRec;  /* First record in scan */
  startTime = -1.0e20;        /* dummy start time */
  /* Initialize end times and record numbers in case of no data */
  endTime = startTime;
  endRec  = startRec; 

  /* Get difference in output reference day and first day in input */
  corDay = refMJD - refDay;
  
  /* Anything to do? */
  if (ntime<=0) return;

  /* Loop over data */
  for (irow=0; irow<ntime; irow++) {

    /* Time in mjd relative to reference day */
    dmjd = ATime[irow];

    /* first time in scan */
    if (startTime<-1000.0) startTime = dmjd - corDay;
    data = outData->buffer;     /* Output data array */
    
    /* Fill record entry in data */
    data[desc->iloct] = dmjd - corDay + offTime;    /* Time (days) */

    /* Get position */
    GetPoint (dmjd+offTime+refDay, &ra, &dec);
    /* Drop data without good pointint value */
    if (dec<-85.0) continue;
    data[desc->ilocti]  = avgTime;      /* time interval (days) */
    data[desc->iloctar] = target;       /* target number */
    data[desc->ilocscan]= scan;         /* scan number */
    if (scanNo!=0) data[desc->ilocscan] = (ofloat)scanNo; 
    data[desc->iloccal] = ACal[irow];   /* Cal? */
    
    /* Parallactic angle (deg) */
    data[desc->ilocrot] = ObitOTFArrayGeomParAng(geom, data[desc->iloct], ra, dec);
    
    /* correction to position */
    ObitOTFArrayGeomCorrPoint(azOff, elOff, data[desc->ilocrot], &ra, &dec);
    data[desc->ilocra]  = ra;           /* RA  in deg. */
    data[desc->ilocdec] = dec;          /* Dec in deg.*/
    
    /* Copy data */	
    for (i=0; i<ndetect; i++) {
      if (beamOffset[i*2]>-9000.0) {/* Any valid data? */
	data[desc->ilocdata+i*incdatawt]   = AData[i][irow];
	data[desc->ilocdata+i*incdatawt+1] = 1.0; /* Initial Weight */
	if (AData[i][irow]==fblank) data[desc->ilocdata+i*incdatawt+1] = 0.0;
      } else {
	data[desc->ilocdata+i*incdatawt]   = fblank;
	data[desc->ilocdata+i*incdatawt+1] = 0.0; /* Initial Weight */
      }
    }
    
    /* set number of records */
    desc->numRecBuff = 1;
    
    /* Write output  */
    if ((ObitOTFWrite (outData, NULL, err) != OBIT_IO_OK) || (err->error>0))
      Obit_log_error(err, OBIT_Error, "ERROR writing output Table file");
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
  } /* end loop over table */  

  /* Get end times and record numbers */
  desc->firstRec = desc->nrecord+1;
  endTime = dmjd - corDay;
  endRec  = desc->nrecord;  /* Last record in scan */

  /* Cleanup */
  if (ATime) g_free(ATime);
  if (ACal)  g_free(ACal);
  if (AData) {
    for (i=0; i<ndetect; i++) g_free(AData[i]);
    g_free(AData);
  }


} /* end GetData  */

void GetPoint (odouble time, ofloat *ra, ofloat *dec)
/*----------------------------------------------------------------------- */
/*  Get antenna pointing at a given time,                                 */
/*  End points or linear interpolation                                    */
/*   Input:                                                               */
/*      time   The desired MJD                                            */
/*   Output:                                                              */
/*      ra     RA J2000 in degrees                                        */
/*      dec    Dec J2000 in degrees                                       */
/*      myInput  parser object                                            */
/*----------------------------------------------------------------------- */
{
  olong i, best;
  odouble test, delta;
  ofloat w1, w2;

  /* Initial values */
  *ra  = 0.0;
  *dec = 0.0;
  
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
      *ra  = 0.0;
      *dec = -90.0;
    } else {
      *ra = AntRA[best];
      *dec = AntDec[best];
    }
  } else if (time<AntDMJD[best]){ /* interpolate with previous */
    w1 = (AntDMJD[best]-time) / (AntDMJD[best]-AntDMJD[best-1]);
    w2 = 1.0 - w1;
    *ra  = w1 * AntRA[best-1]  + w2 * AntRA[best];
    *dec = w1 * AntDec[best-1] + w2 * AntDec[best];
  } else { /* interpolate with following */
    w1 = (AntDMJD[best+1]-time) / (AntDMJD[best+1]-AntDMJD[best]);
    w2 = 1.0 - w1;
    *ra  = w1 * AntRA[best]  + w2 * AntRA[best+1];
    *dec = w1 * AntDec[best] + w2 * AntDec[best+1];
  }
} /* end GetPoint */

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

void ProcessData (gchar *inscan, ofloat avgTime,
		  olong *ndetect, olong *ntime, ofloat *refDate,
		  gfloat** ATime, gfloat*** AData, gfloat** ACal, 
		  ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Read data for detectors from GB FITS and return arrays of:            */
/*    time = timetag for each datum                                       */
/*    data = datum as float for selected entry in data                    */
/*    cal  = 1.0 if cal on, else 0.0                                      */
/*  All rows in table read                                                */
/*  Inputs:                                                               */
/*     inscan   Scan part of input file name                              */
/*     avgTime  Averaging time in days                                    */
/*  Inputs/output:                                                        */
/*     ndetect    Number of detectors in AData                            */
/*     ntime      Length of arrays ATime, AData, ACal                     */
/*     refDate    Day number of reference day                             */
/*     ATime      Pointer to be set to array filled with time values      */
/*                Time is relative to 0 UT on refDay                      */
/*                Should be g_freeed when done                            */
/*     AData      Pointer to be set to array of selected data entries     */
/*                [detector][time], detector in order in FITS table       */
/*                Should be g_freeed when done                            */
/*     ACal       Pointer to be set to array of cal on flags              */
/*                Should be g_freeed when done                            */
/*   Output:                                                              */
/*      err       Obit return error stack                                 */
/*----------------------------------------------------------------------- */
{
  gchar FullFile[128];
  ObitTableGBTPARDATA* PARtable;
  ObitTimeFilter *tFilt=NULL;
  ofloat freq, *data[NUMDETECT], *cal, fblank=ObitMagicF();
  odouble *time;
  gchar *tab;
  olong i, id, disk, ver, nrow, nraw;
  gboolean calOn=FALSE, calOff=FALSE, bad[NUMDETECT];
  gchar *routine = "ProcessData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(data!=NULL);
  /*g_assert(time!=NULL);*/
  /*g_assert(cal!=NULL);*/

  /* get full file name */
  sprintf (FullFile,"Rcvr_PAR/%s.fits", inscan);
  
  /* Create table structure */
  PARtable = newObitTableGBTPARDATA("Data");
  if (err->error) return;

  /* Setup */
  disk = 2;
  tab  = "DATA";
  ver  = 1; 
  nrow = 1;
  ObitTableSetFITS(PARtable,disk,FullFile,tab,ver,nrow,err);

  /* Open */
  ObitTableGBTPARDATAOpen (PARtable, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, PARtable->name);

  /* How many detectors? */
  *ndetect = PARtable->myDesc->repeat[PARtable->daccountsCol];

  /* Allocate output detector array */
  *AData = g_malloc0((*ndetect)*sizeof(gfloat**));

  /* Init detectors arrays */
  for (id=0; id<(*ndetect); id++) {
    bad[id] = (beamOffset[id*2]<-9000.0);
    data[id] = NULL;
  }

  /* Get all detector data from table */
  ReadDatum (PARtable, *ndetect, bad, &nraw, &time, data, &cal, err);
  if (err->error) Obit_traceback_msg (err, routine, PARtable->name);

  /* Anything to do? */
  if (nraw<=0) {
    *ntime=0; 
    Obit_log_error(err, OBIT_InfoErr, "No data found");
    ObitErrLog(err);
    return;
  }

  /* Loop over detectors */
  for (id=0; id<(*ndetect); id++) {

    /* Time averaging  */
    AverageDatum (bad[id], id, nraw, time, data[id], cal, avgTime, 
		  ntime, refDate, ATime, &(*AData)[id], ACal);

    /* Free raw arrays */
    g_free(data[id]);

    /* remove baseline jumps */
    if (!bad[id]) deJumpDatum (*ntime, (*AData)[id], *ACal);

    /* Time domain filtering 
     First grid data onto a regular grid */

    /* Is cal switching? If so use no filtering */
    if ((!bad[id]) && (!calOn && !calOff)) {
      for (i=0; i<*ntime; i++) {
	if ((*AData)[id][i]!=fblank) {
	  if ((*ACal)[i]) calOn = TRUE;
	  else calOff = TRUE;
	}
      }
    }

    /* Only save time and cal flag on the last detector */
    if (id<((*ndetect)-1)) {
      g_free(*ATime); *ATime = NULL;
      g_free(*ACal);  *ACal = NULL;
    }	
  } /* end loop over table */  
  
  /* Close */
  ObitTableGBTPARDATAClose (PARtable, err);
  if (err->error) Obit_traceback_msg (err, routine, PARtable->name);
  
  /* Cleanup */
  PARtable = ObitTableGBTPARDATAUnref(PARtable);
  tFilt    = ObitTimeFilterUnref(tFilt);

  /* Fit/remove power at freq if not using cal */
  if (!(calOn&&calOff)) {
    freq = 1.41170; /* Hz   */
    FitFreq (freq, *ntime, *ATime, *ndetect, *AData, bad, err);
  }

  /* Free raw arrays */
  if (time)  g_free(time);
  if (cal)   g_free(cal);


} /* end ProcessData  */

  void ReadDatum (ObitTableGBTPARDATA* PARtable, olong ndetect, gboolean *bad, 
		  olong *ntime, gdouble** time, gfloat** data, gfloat** cal, 
		  ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Read data for a single detector from GB FITS and return arrays of:    */
/*    time = timetag for each datum                                       */
/*    data = datum as float for selected entry in data                    */
/*    cal  = 1.0 if cal on, else 0.0                                      */
/*    Blanks BTime data after a cal state switch                          */
/*  If the mean sae value exceeds 5 sigma then all the data is flagged    */
/*  Also blanks BTime after a maxSigma (20) sigma point                   */
/*  All rows in table read                                                */
/*  Inputs:                                                               */
/*     PARtable   Open table with PAR data                                */
/*     ndetect    Number of detectors                                     */
/*  Inputs/output:                                                        */
/*     bad        If true there is no valid data per detector             */
/*     ntime      Length of arrays time, data, cal                        */
/*     time       Pointer to be set to array filled with time values      */
/*                Time is relative to 0 UT on first data in data          */
/*                Should be g_freeed when done                            */
/*     data       Array of pointers to raw data arrays                    */
/*                Each should be g_freeed when done                       */
/*     cal        Pointer to be set to array of cal on flags              */
/*                Should be g_freeed when done                            */
/*   Output:                                                              */
/*      err       Obit return error stack                                 */
/*----------------------------------------------------------------------- */
{
  ObitTableGBTPARDATARow* row;
  ObitIOCode retCode;
  olong irow, i, j;
  ofloat val, mean, sigma, timeErr=0.0, *timeErrSum=NULL;
  odouble saeSum[NUMDETECT], saeSum2[NUMDETECT], saeCnt[NUMDETECT];
  odouble dacSum[NUMDETECT], dacSum2[NUMDETECT], dacCnt[NUMDETECT];
  ofloat maxSigma = 20.0;  /* How many sigma to flag */
  ofloat fblank = ObitMagicF();
  gboolean doBlank=FALSE, calOn=FALSE, lastCalOn=FALSE;
  odouble endBlank=-1000.0, tErr;
  olong calOnCnt=0, calOffCnt=0, timeErrCnt=0;
  gchar *routine = "ReadDatum";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(time!=NULL);
  g_assert(data!=NULL);
  g_assert(cal!=NULL);

  /* Allocate output arrays */
  *ntime = PARtable->myDesc->nrow;
  *time = g_malloc0((*ntime)*sizeof(odouble));
  *cal  = g_malloc0((*ntime)*sizeof(ofloat));
  timeErrSum = g_malloc0((*ntime)*sizeof(ofloat));

  /* Per detector setup */
  for (i=0; i<ndetect; i++) {
    data[i] = g_malloc0((*ntime)*sizeof(ofloat));
    saeSum[i] = saeSum2[i] = saeCnt[i] = 0.0;
    dacSum[i] = dacSum2[i] = dacCnt[i] = 0.0;
  }

  /* Create Row structure */
  row = newObitTableGBTPARDATARow (PARtable);

  /* Loop over table */
  for (irow=1; irow<=PARtable->myDesc->nrow; irow++) {

    /* Read row */
    retCode = ObitTableGBTPARDATAReadRow (PARtable, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, PARtable->name);

    /* Save selected data */	
    (*time)[irow-1] = row->TimeStamp;
    /*(*cal)[irow-1]  = (ofloat)(1 - row->DigInput[0]); ??? logic reversed? */
    (*cal)[irow-1]  = (ofloat)row->DigInput[0];/*??? logic reversed? */
    /* Count cal states */
    if (fabs ((*cal)[irow-1]) > 0.1) calOnCnt++;
    else calOffCnt++;

    /* Accumulate time offset using the 1 pps */
    if (row->DigInput[1]) {
      tErr = row->TimeStamp - (olong)row->TimeStamp;
      tErr *= 86400.0;
      tErr -= (olong)tErr;
      /* Assume within 0.5 sec */
      if (tErr>0.5) tErr = -(1.0 - tErr);
      timeErrSum[timeErrCnt++] = tErr;
    }

    /* Per detector section */
    for (i=0; i<ndetect; i++) {
       data[i][irow-1] = fblank;  /* Initial value */
      if (bad[i]) continue;

      /* Sae statistics */
      val = (ofloat)row->saecounts[i];
      saeSum[i]  += val;
      saeSum2[i] += val*val;
      saeCnt[i]++;
      
      /* Dac statistics */
      val = (ofloat)row->daccounts[i];
      dacSum[i]  += val;
      dacSum2[i] += val*val;
      dacCnt[i]++;
      
      data[i][irow-1] = val;  /* data value */
    }

  } /* end loop over table */  

  /* Cleanup */
  row   = ObitTableGBTPARDATARowUnref(row);

  /* Blank data after a cal change - if this really has switching */
  if ((calOnCnt>((*ntime)/10)) && (calOffCnt>((*ntime)/10))) {
    endBlank=-1000.0;
    for (j=0; j<*ntime; j++) {
      /* Cal on? */
      calOn = fabs ((*cal)[j]) > 0.1;
      if (irow==1) lastCalOn = calOn;
      
      /* Cal state switched? */
      if (lastCalOn != calOn) { /* Cal state switched? */
	endBlank = (*time)[j] + BTime;
      }
      /* In blanking period? */
      doBlank = ((*time)[j]<endBlank);
      if (doBlank) for (i=0; i<ndetect; i++) data[i][j] = fblank;

      lastCalOn = calOn;  /* last cal state */
    } /* end loop over times */
    /* end cal modulated */
  } else { /* No cal modulation - set cal off for everything */
    for (j=0; j<*ntime; j++) (*cal)[j] = 0.0;
  }

  /* Fix Time errors > 1 msec */
  if (timeErrCnt>0) {
    timeErr = medianAvg(timeErrSum, 1, timeErrCnt/4, FALSE, timeErrCnt);
    if (timeErrSum) g_free(timeErrSum);
  }
  if (fabs(timeErr)>(0.001)) {
    Obit_log_error(err, OBIT_InfoErr, "Correcting timestamps by %5.1lf msec", 
		   timeErr*1000.0);
    ObitErrLog(err);
    timeErr /= 86400.0;   /* to days */
    for (i=0; i<(*ntime); i++) {
      (*time)[i] -= timeErr;
    }
  }

  /* Per detector checks  */
  for (i=0; i<ndetect; i++) {
    if (bad[i]) continue;

    /* Check for unlocked detector (sae) */
    if (saeCnt[i]>2.0) {
      mean = fabs (saeSum[i] / saeCnt[i]);
      sigma = sqrt ((saeSum2[i]/saeCnt[i]) - mean*mean);
      if (mean>750.0)     bad[i] = TRUE;
      if (mean>5.0*sigma) bad[i] = TRUE;
      if (mean==0.0)      bad[i] = TRUE;
      if (sigma==0.0)     bad[i] = TRUE;
    }
    if (bad[i]) continue;  /* If bad, check no further */

    /* Check for DAC in range, big outliers */
    if (dacCnt[i]>2.0) {
      mean = dacSum[i] / dacCnt[i];
      sigma = MAX (10.0, sqrt ((dacSum2[i]/dacCnt[i]) - mean*mean));
      if (mean==0.0) bad[i] = TRUE;
      if (sigma==0.0) bad[i] = TRUE;
      endBlank=-1000.0;
      for (j=0; j<*ntime; j++) {
	/* DAC in range [1,16383] */
	if ((data[i][j]<1.0) || (data[i][j]>16383.0)) data[i][j] = fblank;
	/* Big outlier */
	if ((data[i][j]!=fblank) && (fabs(data[i][j]-mean) > maxSigma*sigma)) {
	  /*  Kill next BTime */
	  endBlank = (*time)[i] + BTime;
	}
	doBlank = (*time)[i] < endBlank;
	if (doBlank) data[i][j] = fblank;
      }
    }
  } /* End per detector checks */ 
} /* end ReadDatum  */

void AverageDatum (gboolean bad,  olong detector, olong nraw, odouble *time, 
		   ofloat *data, ofloat *cal, ofloat avgTime, olong *ntime, 
		   ofloat *refDate, ofloat **ATime, ofloat **AData, ofloat **ACal)
/*----------------------------------------------------------------------- */
/*  Average time filtered PAR data for one detector                       */
/*  Inputs:                                                               */
/*     bad        If true there is no valid data                          */
/*     detector   0-rel entry number of data array                        */
/*     nraw       Length of arrays time, data, cal                        */
/*     time       Pointer to be set to array filled with time values      */
/*                Time is relative to 0 UT on first data in data          */
/*     data       Pointer to array of selected data entries               */
/*     cal        Pointer to array of cal on flags                        */
/*     avgTime    Averaging time in days                                  */
/*  Output:                                                               */
/*     ntime      Length of arrays ATime, AData, ACal                     */
/*     ATime      Pointer to be set to array filled with time values      */
/*                Time is relative to 0 UT on first data in data          */
/*                Should be g_freeed when done                            */
/*     AData      Pointer to be set to array of selected data entries     */
/*                Should be g_freeed when done                            */
/*     ACal       Pointer to be set to array of cal on flags              */
/*                Should be g_freeed when done                            */
/*----------------------------------------------------------------------- */
{
  olong i, count, tcount, next;
  olong itemp;
  ofloat sumData, sumCal, defaultCal=1.0, fblank = ObitMagicF();
  odouble sumTime, tnext, tlast;

  /* error checks */
  g_assert(time!=NULL);
  g_assert(data!=NULL);
  g_assert(cal!=NULL);
  g_assert(ATime!=NULL);
  g_assert(AData!=NULL);
  g_assert(ACal!=NULL);

  /* How many times? Fudge a bit */
  *ntime = 10 + (time[nraw-1] - time[0]) / avgTime;

  /* Get reference day number */
  itemp = time[0];
  *refDate = (ofloat)(itemp) ;

  /* Allocate output arrays */
  *ATime = g_malloc0((*ntime)*sizeof(ofloat));
  *AData = g_malloc0((*ntime)*sizeof(ofloat));
  *ACal  = g_malloc0((*ntime)*sizeof(ofloat));

  /* This one all bad? */
  if (bad && (detector<(ndetect-1))) {
    for (i=0; i<*ntime; i++) {
      (*ATime)[i]  = 0.0;
      (*AData)[i]  = fblank;
      (*ACal)[i]   = 0.0;
    }
    return;
  }

  /* Loop over input data */
  tlast = time[0];
  tnext = time[0] + avgTime;
  next = 0;
  sumTime = 0.0; sumData = sumCal = 0.0; count = 0; tcount = 0;
  for (i=0; i<nraw; i++) {
    /* Filled arrays? */
    if (next>=(*ntime)) break;
    /* Ignore data before current integration - not a totally satifactory fix */
    if (time[i]<(tnext-1.01*avgTime)) continue;
    /* Average finished? */
    if (time[i]>tnext) {
      if (tcount>0) (*ATime)[next]  = (ofloat)(sumTime / tcount);
      else (*ATime)[next]  = tnext - 0.5*avgTime; /* Middle of accumulation */
      if (sumCal>0.0) (*ACal)[next] = 1.0;
      else if (count>0) (*ACal)[next] = 0.0; /* data but no cal */
      else (*ACal)[next] = defaultCal;
      /* Flag any data which is not all cal on or cal off */
      if ((count>0) && ((fabs(sumCal-count)<0.1) || (sumCal==0.0))) {
	(*AData)[next]  = sumData / count;
      } else { /* no data this interval  */
	(*AData)[next] = fblank;
      }
      tlast = time[i-1];
      tnext += avgTime;
      sumTime = 0.0; sumData = sumCal = 0.0; count = 0; tcount = 0;
      if ((*AData)[next]!=fblank) defaultCal = (*ACal)[next]; /* Default cal */
      next++;
      /* Is time still in current accumulation? */
      while (time[i]>tnext) {
	(*ATime)[next]   = tlast - (*refDate);
	(*AData)[next++] = fblank;
	tnext += avgTime;
	tlast += avgTime;
      }
    } /* End if average finished */

    /* Accumulate */
    tcount++;
    sumTime += (time[i] - (*refDate));
    sumCal  += cal[i];
    if (data[i]!=fblank) {
      count++;
      sumData += data[i];
    }
  } /* end loop over data */

  /* Actual number of times output */
  *ntime = MIN ((*ntime), next);

} /* end AverageDatum  */

void deJumpDatum (gint ndata, ofloat *data, ofloat *cal)
/*----------------------------------------------------------------------- */
/*  Remove jumps in data baseline                                         */
/*  Determines a 9 point running median and when there is a persistent    */
/*  jump not associated with a cal  state change, then the following      */
/*  data is adjusted by the difference.                                   */
/*  Multiple jumps may be detected                                        */
/*  Inputs:                                                               */
/*     ndata      Length of arrays data, cal                              */
/*     cal        Array of cal on flags                                   */
/*  Inputs/output:                                                        */
/*     data       Array of selected data entries                          */
/*----------------------------------------------------------------------- */
{
  olong i, j, count, p, f, ncopy, nrun=9, half=nrun/2;
  ofloat tmed, *tdata=NULL, *twork=NULL, maxDelt, delta;
  ofloat val, sum, sum2, fblank = ObitMagicF();

  if (ndata<=nrun) return;  /* Enough data to bother with? */

  /* error checks */
  g_assert(data!=NULL);
  g_assert(cal!=NULL);

  /* Allocate work arrays */
  tdata = g_malloc0((ndata)*sizeof(ofloat));
  ncopy = nrun*sizeof(ofloat);
  twork = g_malloc0(ncopy);

  /* Fill with running median */
  memcpy (twork, data, ncopy);  /* medianValue reorders the data */
  tmed = medianValue(twork, 1, nrun);
  for (i=0; i<half; i++) tdata[i] = tmed;  /* start */

  for (i=half; i<ndata-half; i++) { /* middle */
    memcpy (twork, &data[i-half], ncopy);
    tmed = medianValue(twork, 1, nrun);
    tdata[i] = tmed;
  }
  
  for (i=ndata-half; i<ndata; i++) tdata[i] = tmed;  /* end */

  /* RMS of difference values from median */
  sum = sum2 = 0.0; count = 0;
  for (i=0; i<ndata; i++) {
    if ((data[i]!=fblank) && (tdata[i]!=fblank)) {
      val = data[i]-tdata[i];  /* Difference from running median */
      sum  += val;
      sum2 += val*val;
      count++;
    }
  }
  /* maximum jump allowed 100 sigma */
  sum  /= MAX (1, count);
  sum2 /= MAX (1, count);
  maxDelt = 100.0 * sqrt(sum2 - sum*sum);

  /* loop looking for jump */
  for (i=0; i<ndata-1; i++) {
    if ((tdata[i]!=fblank) && (tdata[i+1]!=fblank) && fabs(tdata[i+1]-tdata[i])>maxDelt) {
      /* require also a jump in data but no change in cal */
      if ((data[i]!=fblank) && (data[i+1]!=fblank) && 
	   (fabs(data[i+1]-data[i])>(0.75*maxDelt)) && (cal[i+1]==cal[i])) {
	p = MAX (0, i-half);
	f = MIN (ndata-1, i+half);
	delta = tdata[p]-tdata[f];
	/* offset all following values by delta */
	for (j=i+1; j<ndata; j++) {
	  data[j]  += delta;
	  tdata[j] += delta;
	}
	/* Blank this pair for good measure */
	data[i] = fblank;
	data[i+1] = fblank;
      } /* end jump in data */
    } /* end jump in median */
  }  /* end loop looking for jumps */

  /* cleanup */
  if (tdata) g_free(tdata);
  if (twork) g_free(twork);
} /* end deJumpDatum  */

/**
 * IIF filter for 1.4 Hz oscillations in 10 Hz sampled data
 * See 
 * http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
 * Order 3 bessel bandstop filter 1.38-1.45 Hz
 * The first valid data point will be subtracted from the data
 * \param ndata    Length of data
 * \param data     data to be filtered, may be magic value blanked
 */
void IIRFilter10 (gint ndata, ofloat *data)
{
#define NZEROES 6
#define NPOLES  6
  ofloat GAIN = 1.039215778;
  ofloat xv[NZEROES+1], yv[NPOLES+1], last, subx, dd;
  ofloat fblank = ObitMagicF();
  olong i;
  
  /* Subtract first valid point from data */
  subx = data[0];
  last = fblank;
  for (i=0; i<ndata; i++) {
    if (data[i] != fblank) {
      if (last == fblank) { /* Save first good value */
	last = data[i];
	subx = last;
      }
      data[i] -= subx;
    }
  }

  /* Any good data? */
  if (last==fblank) return;

  /* Init internal registers */
  for (i=0; i<=NZEROES; i++) xv[i] = 0.0;
  for (i=0; i<=NPOLES;  i++) yv[i] = 0.0;

  /* Loop over data */
  last = 0.0;
  for (i=0; i<ndata; i++) {
    if (data[i]!=fblank) {  /* Good datum? */
      dd = data[i];
      last = dd;
    } else {  /* No use last valid value */
      dd = last;
    }
    xv[0] = xv[1]; xv[1] = xv[2]; xv[2] = xv[3]; xv[3] = xv[4]; xv[4] = xv[5]; xv[5] = xv[6];
    xv[6] = dd/GAIN;
    yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3]; yv[3] = yv[4]; yv[4] = yv[5]; yv[5] = yv[6];
    yv[6] =   (xv[0] + xv[6]) -   3.7817176297 * (xv[1] + xv[5]) +   7.7671294101 * (xv[2] + xv[4]) 
      -   9.5665394074 * xv[3] 
      + ( -0.9256760402 * yv[0]) + (  3.5460530638 * yv[1]) 
      + ( -7.3772622509 * yv[2]) + (  9.2040900878 * yv[3]) 
      + ( -7.5696488979 * yv[4]) + (  3.7334159240 * yv[5]);
    if  (data[i]!=fblank) data[i] = yv[6];
  } /* end loop over data */

} /* end IIRFilter10 */

/**
 * IIF filter for 1.4 Hz oscillations in 20 Hz sampled data
 * See 
 * http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
 * Order 3 bessel bandstop filter 1.38-1.45 Hz
 * The first valid data point will be subtracted from the data
 * \param ndata    Length of data
 * \param data     data to be filtered, may be magic value blanked
 */
void IIRFilter20 (gint ndata, ofloat *data)
{
#define NZEROES2 6
#define NPOLES2  6
  ofloat GAIN = 1.019454963;
  ofloat xv[NZEROES2+1], yv[NPOLES2+1], last, subx, dd;
  ofloat fblank = ObitMagicF();
  olong i;
  
  /* Subtract first valid point from data */
  subx = data[0];
  last = fblank;
  for (i=0; i<ndata; i++) {
    if (data[i] != fblank) {
      if (last == fblank) { /* Save first good value */
	last = data[i];
	subx = last;
      }
      data[i] -= subx;
    }
  }

  /* Any good data? */
  if (last==fblank) return;

  /* Init internal registers */
  for (i=0; i<=NZEROES2; i++) xv[i] = 0.0;
  for (i=0; i<=NPOLES2;  i++) yv[i] = 0.0;

  /* Loop over data */
  last = 0.0;
  for (i=0; i<ndata; i++) {
    if (data[i]!=fblank) {  /* Good datum? */
      dd = data[i];
      last = dd;
    } else {  /* No use last valid value */
      dd = last;
    }
    xv[0] = xv[1]; xv[1] = xv[2]; xv[2] = xv[3]; xv[3] = xv[4]; xv[4] = xv[5]; xv[5] = xv[6]; 
    xv[6] = dd/GAIN;
    yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3]; yv[3] = yv[4]; yv[4] = yv[5]; yv[5] = yv[6]; 
    yv[6] =   (xv[0] + xv[6]) -   5.4171909273 * (xv[1] + xv[5]) +  12.7819858470 * (xv[2] + xv[4])
      -  16.7222579640 * xv[3]
      + ( -0.9621250877 * yv[0]) + (  5.2456848138 * yv[1])
      + (-12.4572294670 * yv[2]) + ( 16.4026074470 * yv[3])
      + (-12.6185948320 * yv[4]) + (  5.3824651695 * yv[5]);
    
    if  (data[i]!=fblank) data[i] = yv[6];
  } /* end loop over data */

} /* end IIRFilter20 */

/**
 * Estimate amplitude and phase of a given frequency for each detector 
 * and remove from data stream.
 * Filter 10 sigma deviations from an alpha=0.5 50 cell running median average
 * \param  freq       Frequency to remove (Hz)
 * \param  ntime      Length of arrays ATime, AData
 * \param  ATime      Pointer to  array filled with time values (day)
 * \param  ndetect    Number of detectors in AData
 * \param  AData      Pointer to  array of data entries
 *                    [detector][time], detector in order in FITS table.
 *                    Modified on output
 * \param  bad        Array of detector valid flags 
 * \param  err        Obit Error/message stack  
 */
void FitFreq (ofloat freq, olong ntime, ofloat *Atime, 
	      olong ndetect, ofloat **AData, gboolean *bad, ObitErr *err)
{
  ofloat fblank = ObitMagicF();
  olong i, j, k, wind;
  ofloat time, phase, dphase, avgPhs, *fdata;
  ofloat ss, cs, cr, cnt, alpha, RMS;
  ofloat *sine=NULL, *cosine=NULL, *amp=NULL, *phs=NULL, *work1=NULL, *work2=NULL;

  work1 = g_malloc0(ntime*sizeof(ofloat));
  work2 = g_malloc0(ntime*sizeof(ofloat));

  /* Generate sine, cosine patterns */
  sine   = g_malloc0(ntime*sizeof(ofloat));
  cosine = g_malloc0(ntime*sizeof(ofloat));
  dphase = 2.0 * G_PI * freq * 86400.0;
  for (i=0; i<ntime; i++) {
    time = Atime[i] - Atime[0];
    phase = time*dphase;
    sine[i]   = sin(phase);
    cosine[i] = cos(phase);
  } /* end loop filling arrays */

  /* Subtract first point from each series */
  for (j=0; j<ndetect; j++) {
    fdata = AData[j];
    for (i=1; i<ntime; i++) {
      if (fdata[i]!=fblank) {
	fdata[i] -= fdata[0];
      }
    } /* end loop over time */
    fdata[0] = 0.0;
  } /* end loop over detector */

  wind  = 50;
  alpha = 0.5;

  /* Loop over detectors getting dot product */
  amp = g_malloc0(ndetect*sizeof(ofloat));
  phs = g_malloc0(ndetect*sizeof(ofloat));
  for (j=0; j<ndetect; j++) {
    if (bad[j]) continue;  /* Known to be bad? */
    fdata = AData[j];

    /* Filter data stream - Running median in work1 */
    RunningMedian (ntime, wind, fdata, alpha, &RMS, work1, work2);
    /* Filter data to work2 */
    for (k=0; k<ntime; k++) {
      /* Flagged or deviant? */
      if ((fdata[k]==fblank) || (work1[k]==fblank) || 
	  (fabs(fdata[k]-work1[k])>10.0*RMS)) {
	work2[k] = fblank;
	/* DEBUG
	if ((fdata[k]!=fblank) && (work1[k]!=fblank) &&
	    (fabs(fdata[k]-work1[k])>10.0*RMS))
	  fprintf (stdout,"det %d pt %d %f %f %f\n",j,k,fdata[k],work1[k],10.0*RMS); */
	/* DEBUG */
      } else 
	work2[k] = fdata[k];  /* OK */
    }

    ss = cs = cnt = 0.0;
    for (i=0; i<ntime; i++) {
      if (work2[i]!=fblank) {
	cnt += 0.5;
	ss  += work2[i] * sine[i];
	cs  += work2[i] * cosine[i];
      }
    } /* end loop over time */
    
    /* Normalize */
    if (cnt!=0.0) {
      ss /= cnt;
      cs /= cnt;
      amp[j] = sqrt(ss*ss + cs*cs);
      phs[j] = atan2(ss, cs+1.0e-20);
      /* Tell
      fprintf (stdout, "detector %d amp %f phase %f ss %f cs %f\n",
	       j, amp[j], 57.296*phs[j], ss, cs); */
    }
  } /* end loop over detector */

  avgPhs = AveragePhase (ndetect, amp, phs, bad); 
  Obit_log_error(err, OBIT_InfoErr, "Avgerage PT phase %f", avgPhs*57.296);
  ObitErrLog(err);
 
  /* Correct using average phase */
  for (j=0; j<ndetect; j++) {
    /* DEBUG
    amp[j] = fabs(amp[j]); avgPhs = phs[j];
    fprintf (stdout, "det %d amp %f avgPhas %f\n", j+1, amp[j], avgPhs*57.296); */
    /* DEBUG */
    cs = amp[j] * cos(avgPhs);
    ss = amp[j] * sin(avgPhs);
    fdata = AData[j];
    for (i=0; i<ntime; i++) {
      if (fdata[i]!=fblank) {
	cr = cs * cosine[i] + ss * sine[i];
	fdata[i] -= cr;
      }
    } /* end time loop */
  } /* end loop over detector */
 
  /* Cleanup*/
  if (sine)   g_free(sine);
  if (cosine) g_free(cosine);
  if (amp)    g_free(amp);
  if (phs)    g_free(phs);
  if (work1)  g_free(work1);
  if (work2)  g_free(work2);
  
} /* end FitFreq */

/**
 * Average phases and adjust sign of the amplitude
 * \param ndetect  Number of detectors in amp, phase
 * \param amp      Amplitudes, on output may be negative
 * \param phase    Phases (radians)
 * \param bad      Array of detector valid flags 
 * \return average phase (radians)
 */
ofloat AveragePhase (olong ndetect, ofloat *amp, ofloat *phase, gboolean *bad)
{
  olong i, id, best, dist[30], ndist=30;
  ofloat tphs, roughPhase, sum, sumwt, avgPhase;

  /* Get distribution in range 0-pi */
  for (i=0; i<ndist; i++) dist[i] = 0;
  for (i=8; i<ndetect; i++) {
    if ((amp[i]>1.0) && !bad[i]) {  /* ignore poor detectors */
      tphs = phase[i];
      if (tphs<0.0) tphs += G_PI; /* Fold to positive */
      id = (olong)(0.5 + ndist * tphs/G_PI);
      id = MAX (0, MIN (ndist-1, id));
      dist[id]++;
    } /* end if good detector */
  } /* end loop over detectors */

  /* Find mode of distribution */
  best = dist[0];
  id = 0;
  for (i=1; i<ndist; i++) {
    if (dist[i]>best) {
      best = dist[i];
      id = i;
    }
  }
  roughPhase = id * G_PI/ndist;

  /* Average phases flipping amplitudes and phases depending 
     on whether the phase is within pi of roughPhase */
  sum = sumwt = 0.0;
  for (i=0; i<ndetect; i++) {
    if ((amp[i]>1.0) && !bad[i]) {  /* ignore poor detectors */
      tphs = phase[i];
      /* Work out relative sign of gain - more than half turn? */
      if (fabs(tphs-roughPhase)>0.5*G_PI) { 
	/* Is the difference closer to a turn than half? */
	if ((tphs-roughPhase)>1.5*G_PI) {
	  tphs -= 2.0*G_PI;
	} else if ((tphs-roughPhase)<(-1.5*G_PI)) {
	  tphs += 2.0*G_PI;
	} else {  /* other sign */
	  amp[i] = -fabs(amp[i]);
	  if (tphs>roughPhase) tphs -= G_PI;
	  else tphs += G_PI;
	}
      }
      sum   += fabs(amp[i])*tphs;
      sumwt += fabs(amp[i]);
    } /* end if good detector */
  } /* end loop over detectors */
  
  /* Get average */
  if (sumwt>0.0) avgPhase = sum / sumwt;
  else avgPhase = 0.0;
  
  return avgPhase;
} /* end AveragePhase */
