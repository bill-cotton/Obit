/* $Id:  $  */
/* List contents of ASDM SysPower Table                              */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
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

#include <unistd.h>
#include "ObitSystem.h"
#include "ObitAIPSDir.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitSDMData.h"
#include "ObitPrinter.h"


/* internal prototypes */
/* Get inputs */
ObitInfoList* SysPowerViewin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Days to human string */
void day2dhms(ofloat time, gchar *timeString);
/*  Dump SysPower Table */
void Dump (ObitInfoList *myInput, ObitSDMData *SDMData, ObitErr* err);

/* Program globals */
gchar *pgmName = "SysPowerView";       /* Program name */
gchar *infile  = NULL; /* File with program inputs */
gchar *outfile = NULL; /* File to contain program outputs */
olong  pgmNumber=1;    /* Program number (like POPS no.) */
olong  AIPSuser=2;     /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;        /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;        /* Number of FITS directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
gchar **FITSdirs=NULL; /* List of FITS data directories */
odouble refJD = 0.0;   /* reference Julian date */
odouble refMJD = 0.0;  /* reference Modified Julian date */
gboolean isEVLA;                      /* Is this EVLA data? */
gboolean isALMA;                      /* Is this ALMA data? */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*    Summarize contents of an ASDM                                       */
/*----------------------------------------------------------------------- */
{
  olong  i, ierr=0;
  ObitSystem *mySystem= NULL;
  ObitErr *err= NULL;
  gchar dataroot[132];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitSDMData *SDMData=NULL;

  err = newObitErr();  /* Obit error/message stack */

  /* Startup - parse command line */
  ierr = 0;
  myInput = SysPowerViewin (argc, argv, err);
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

  /* List it */
  Dump (myInput, SDMData, err);
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;
  
  /* Shutdown Obit */
 exit:
  if (outfile)
    ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  mySystem = ObitSystemShutdown (mySystem);
  
  /* cleanup */
  myInput  = ObitInfoListUnref(myInput);   /* delete input list */
  myOutput = ObitInfoListUnref(myOutput);  /* delete output list */
  SDMData  = ObitSDMDataUnref (SDMData);
 
  return ierr;
} /* end of main */

ObitInfoList* SysPowerViewin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file=NULL, *arg;
  gboolean init=FALSE;
  oint itemp, iarray[2];
  ofloat farray[2];;
  gchar *strTemp;
  ObitInfoList* list;
  gchar *routine = "SysPowerViewin";

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

    } else if (strcmp(arg, "-prtFile") == 0){ /* Output text file */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "prtFile", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-doCrt") == 0) { /* terminal/file output  */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "doCrt", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-pgmNumber") == 0) { /* Program number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "pgmNumber", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-timeRange") == 0) { /* Time range */
      dim[0] = 2;
      farray[0] = (ofloat)strtof(argv[++ax], NULL);
      farray[1] = (ofloat)strtof(argv[++ax], NULL);
      ObitInfoListAlwaysPut (list, "timeRange", OBIT_float, dim, farray);
    } else if (strcmp(arg, "-SWId") == 0) { /* SW Id range */
      iarray[0] = (oint)strtol(argv[++ax], NULL, 0);
      iarray[1] = (oint)strtol(argv[++ax], NULL, 0);
      ObitInfoListAlwaysPut (list, "SWId", OBIT_float, dim, iarray);
    } else if (strcmp(arg, "-AntId") == 0) { /* Antenna Id range */
      iarray[0] = (oint)strtol(argv[++ax], NULL, 0);
      iarray[1] = (oint)strtol(argv[++ax], NULL, 0);
      ObitInfoListAlwaysPut (list, "AntId", OBIT_float, dim, iarray);
    } else { /* unknown argument */
      /* DEBUG fprintf (stderr,"DEBUG parameter %s \n",arg);*/
      Usage();
    }
  } /* end parsing input arguments */

  /* Initialize output */
  myOutput = defaultOutputs(err);
  if (outfile)
    ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

 return list;
} /* end SysPowerViewin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: SysPowerView -input file -output ofile [args]\n");
    fprintf(stderr, "List an ASDM SysPower Table \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input optional input parameter file\n");
    fprintf(stderr, "  -output output result file, def SysPowerView.out\n");
    fprintf(stderr, "  -DataRoot Directory name for root of ASDM \n");
    fprintf(stderr, "  -doCrt terminal/file output [def 132]\n");
    fprintf(stderr, "  -prtFile output text file [def. print.list]\n");  
    fprintf(stderr, "  -timeRange start end [time range in days]\n");
    fprintf(stderr, "  -AntId start end [ant. Id range]\n");
    fprintf(stderr, "  -SWId start end [SW Id range]\n");
    
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
  oint   itemp, iarray[2];
  ofloat farray[2];
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

  /* root of data directory */
  strTemp = "dataRoot";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "DataRoot", OBIT_string, dim, strTemp, err);

  /* doCrt */
  dim[0] = 1; dim[1] = 1;
  itemp = 132;
  ObitInfoListPut (out, "doCrt", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* print file */
  strTemp = "print.list";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "prtFile", OBIT_string, dim, strTemp, err);

  /* Timerange */
  farray[0] = farray[1] = 0.0;
  dim[0] = 2; dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut (out, "timeRange", OBIT_float, dim, farray);

  /* SW Id range */
  iarray[0] = iarray[1] = 0;
  dim[0] = 2; dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut (out, "SWId", OBIT_oint, dim, iarray);

  /* Ant Id range */
  iarray[0] = iarray[1] = 0;
  dim[0] = 2; dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut (out, "AntId", OBIT_oint, dim, iarray);

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

/**
 *  Dump contents of ASDM SysPower Table
 * \param myInput Input parameters on InfoList 
 * \param SDMData ASDM data structure
 * \param err     Obit Error stack 
*/
void Dump (ObitInfoList *myInput, ObitSDMData *SDMData, ObitErr *err)
{
  ObitPrinter  *myPrint = NULL;
  FILE         *outStream = NULL;
  gboolean     isInteractive = FALSE, quit = FALSE;
  gchar        line[1024], Title1[1024], Title2[1024];
  gchar        DataRoot[256], obsdat[25];
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *prtFile=NULL, begString[25], endString[25];
  ASDMSysPowerTable *SyPwrTab;
  ASDMAntennaArray*  AntArray;
  ASDMSourceArray*   SourceArray;
  olong        pLimit=1000000;  /* Page limit */
  olong        i, doCrt=0, LinesPerPage=0,  nants;
  olong        iTab, iScan, SWId[2], AntId[2], swid, antid, ia;
  gchar        *ants[100];
  gboolean     want;
  odouble      refJD = SDMData->refJD;
  ofloat       time, TR[2];
  gchar        *routine = "Dump";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitSDMDataIsA(SDMData));
  
  /* Get parameters */
  /* input DataRoot name */
  for (i=0; i<132; i++) DataRoot[i] = 0;
  ObitInfoListGet(myInput, "DataRoot", &type, dim, DataRoot, err);
  DataRoot[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(DataRoot);  /* Trim trailing blanks */
  
  /* timeRange */
  ObitInfoListGet(myInput, "timeRange",&type, dim, TR,  err);
  if (TR[1]<=TR[0]) TR[1] = 10.0;

  /* SWId range */
  ObitInfoListGet(myInput, "SWId",&type, dim, SWId,  err);
  if (SWId[1]==0) SWId[1] = 10000;

  /* AntId range */
  ObitInfoListGet(myInput, "AntId",&type, dim, AntId,  err);
  if (AntId[1]==0) AntId[1] = 100;

  ObitInfoListGet(myInput, "doCrt", &type, dim, &doCrt, err);
  isInteractive = doCrt>0;
  if (isInteractive) { /* interactive, write to stdout */
    outStream = stdout;
  } else { /* write to file */
    ObitInfoListGetP(myInput, "prtFile",  &type, dim, (gpointer)&prtFile);
    if (prtFile) prtFile[MIN (127,dim[0])] = 0;
    ObitTrimTrail(prtFile);
    /* Make sure file named */
    Obit_return_if_fail(((strlen(prtFile)>2) && 
			 ((prtFile[0]!=' ') && (prtFile[1]!=' '))), err,
			"%s: Printer file not specified", routine);
    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Writing output to %s",prtFile);
    ObitErrLog(err);
  }
  if (err->error) Obit_traceback_msg (err, routine, "myInput");

  /* Reference date */
  iScan       = SDMData->MainTab->rows[0]->scanNumber;
  AntArray    = ObitSDMDataGetAntArray(SDMData, iScan);
  refJD       = AntArray->refJD;
  ObitUVDescJD2Date (AntArray->refJD, obsdat);

  /* List of antenna names */
  for (i=0; i<100; i++) ants[i] = "unkn";
  nants = AntArray->nants;
  for (i=0; i<AntArray->nants; i++) {
    ia = AntArray->ants[i]->antennaId;
    ants[ia] = AntArray->ants[i]->antName;
  }

  /* Source Info */
  SourceArray = ObitSDMDataGetSourceArray(SDMData);
 
  /* Is this the EVLA? */
  isEVLA = !strncmp(AntArray->arrayName, "EVLA", 4);
  /* Is this the ALMA? */
  isALMA = !strncmp(AntArray->arrayName, "ALMA", 4);
  /* ALMA is a bit confused */
  if (!isALMA) isALMA = !strncmp(AntArray->arrayName, "OSF", 3);

  /* Titles */
  sprintf (Title1, "Dump of ASDM SysPower table in %s",DataRoot);
  sprintf (Title2, "  ");
  
  /* Create/Open Printer */
  myPrint = ObitPrinterCreate (pgmName, isInteractive, outStream, prtFile);
  /* page limit */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myPrint->info, "pageLimit", OBIT_long, dim, &pLimit);
  ObitPrinterOpen (myPrint, LinesPerPage, Title1, Title2, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
  
  /* Titles */
  ObitPrinterWrite (myPrint, Title1, &quit, err);
  ObitPrinterWrite (myPrint, Title2, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
  
  sprintf(line,"Observed %s, Array = %s, Observer = %s", 
	  obsdat, AntArray->arrayName, AntArray->obsName);
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Selection */
  day2dhms(TR[0], begString);
  day2dhms(TR[1], endString);
  sprintf(line,"Selected timerange: %s - %s", begString, endString);
  ObitPrinterWrite (myPrint, line, &quit, err);
  sprintf(line,"Selected ant Ids:   %3d - %3d", AntId[0], AntId[1]);
  ObitPrinterWrite (myPrint, line, &quit, err);
  sprintf(line,"Selected SW Ids:    %3d - %3d", SWId[0],  SWId[1]);
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Loop over SysPower Table */
  SyPwrTab = SDMData->SysPowerTab;
  for (iTab=0; iTab<SyPwrTab->nrows; iTab++) {
    time  = SyPwrTab->rows[iTab]->timeInterval[0]-refJD;
    swid  = SyPwrTab->rows[iTab]->spectralWindowId;
    antid = SyPwrTab->rows[iTab]->antennaId;
    want = ((time>=TR[0])     && (time<=TR[1])) &&
           ((swid>=SWId[0])   && (swid<=SWId[1])) &&
           ((antid>=AntId[0]) && (antid<=AntId[1]));
    if (want) {
      /* Timerange in human form */
      day2dhms(time, begString);
      g_snprintf (line, 130, "time %s AntID %2d (%s)  SWId %3d SPDif %6.3f SPSum %6.1f RG %6.3f             ", 
		  begString, antid, ants[antid], swid,
		  SyPwrTab->rows[iTab]->switchedPowerDifference[0],
		  SyPwrTab->rows[iTab]->switchedPowerSum[0],
		  SyPwrTab->rows[iTab]->requantizerGain[0]);
      /* Second Poln? */
      if (SyPwrTab->rows[iTab]->numReceptor>1) {
	g_snprintf (&line[79], 130, "  pol 2: SPDif %6.3f SPSum %6.1f RG %6.3f ", 
		    SyPwrTab->rows[iTab]->switchedPowerDifference[1],
		    SyPwrTab->rows[iTab]->switchedPowerSum[1],
		    SyPwrTab->rows[iTab]->requantizerGain[1]);
      }
      ObitPrinterWrite (myPrint, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    }  /* end wanted */
  } /* end SysPower Loop */
  /* Done - Close Printer */
 Quit:
  ObitPrinterClose (myPrint, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Cleanup */
  AntArray    = ObitSDMDataKillAntArray(AntArray);
  SourceArray = ObitSDMDataKillSourceArray(SourceArray);

} /* end Dump */

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

