/* $Id$  */
/* Summarize contents of ASDM                                        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010-2014                                          */
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
ObitInfoList* ASDMListin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Days to human string */
void day2dhms(ofloat time, gchar *timeString);
/*  Suymmarize ASDM */
void Summary (ObitInfoList *myInput, ObitSDMData *SDMData, ObitErr* err);

/* Program globals */
gchar *pgmName = "ASDMList";       /* Program name */
gchar *infile  = "ASDMList.inp";   /* File with program inputs */
gchar *outfile = "ASDMList.out";   /* File to contain program outputs */
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
  myInput = ASDMListin (argc, argv, err);
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
  Summary (myInput, SDMData, err);
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
 
  return ierr;
} /* end of main */

ObitInfoList* ASDMListin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="ASDMList.in", *arg;
  gboolean init=FALSE;
  oint itemp;
  gchar *strTemp;
  ObitInfoList* list;
  gchar *routine = "ASDMListin";

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
} /* end ASDMListin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: ASDMList -input file -output ofile [args]\n");
    fprintf(stderr, "Summarize an ASDM data archive \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def ASDMList.in\n");
    fprintf(stderr, "  -output output result file, def ASDMList.out\n");
    fprintf(stderr, "  -DataRoot Directory name for input \n");
    fprintf(stderr, "  -doCrt terminal/file output\n");
    fprintf(stderr, "  -prtFile output text file\n");  
    fprintf(stderr, "  -doCrt terminal/file output\n");
    
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
 *  Summarize contents of ASDM
 * \param myInput Input parameters on InfoList 
 * \param SDMData ASDM data structure
 * \param err     Obit Error stack 
*/
void Summary (ObitInfoList *myInput, ObitSDMData *SDMData, ObitErr *err)
{
  ObitPrinter  *myPrint = NULL;
  FILE         *outStream = NULL;
  gboolean     isInteractive = FALSE, quit = FALSE;
  gchar        line[1024], Title1[1024], Title2[1024];
  gchar        DataRoot[256], obsdat[25];
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *prtFile=NULL, begString[25], endString[25];
  gchar        calCode[24];
  ASDMMainTable *MainTab;
  ASDMScanTable *ScanTab;
  ASDMDataDescriptionTable *DataDescTab;
  ASDMConfigDescriptionTable *ConfigTab;
  ASDMSpectralWindowTable *SpectralWindowTab;
  ASDMAntennaArray*  AntArray;
  ASDMSourceArray*   SourceArray;
  olong        pLimit=1000000;  /* Page limit */
  olong        iScan, iIntent, iConfig, configID, dataDescriptionId;
  olong        i, ii, doCrt, LinesPerPage=0, iDD, jDD, jSW;
  olong        spectralWindowId, ScanID, iMain, iSource;
  gchar        bcode[10];
  gchar        *routine = "Summary";

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
  
  ObitInfoListGet(myInput, "doCrt",     &type, dim, &doCrt,  err);
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

  /* Source Info */
  SourceArray = ObitSDMDataGetSourceArray(SDMData);
 
  /* Is this the EVLA? */
  isEVLA = !strncmp(AntArray->arrayName, "EVLA", 4);
  /* Is this the ALMA? */
  isALMA = !strncmp(AntArray->arrayName, "ALMA", 4);
  /* ALMA is a bit confused */
  if (!isALMA) isALMA = !strncmp(AntArray->arrayName, "OSF", 3);

  /* Titles */
  sprintf (Title1, "Summary of ASDM in %s",DataRoot);
  sprintf (Title2, "  Configuration Summary");
  
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
  
  sprintf(line,"Configurations");
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  sprintf(line,"      ");
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Loop over configurations */
  ConfigTab         = SDMData->ConfigDescriptionTab;
  DataDescTab       = SDMData->DataDescriptionTab;
  SpectralWindowTab = SDMData->SpectralWindowTab;
  
  for (iConfig=0; iConfig<ConfigTab->nrows; iConfig++) {
    sprintf(line,"  Configuration  %d, no. SpWin = %d  ", 
	    iConfig, ConfigTab->rows[iConfig]->numDataDescription);
    ObitPrinterWrite (myPrint, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    /* Loop over data descriptions to get spectral windows */
    for (iDD=0; iDD<ConfigTab->rows[iConfig]->numDataDescription;  iDD++) {
      
      /* Find it in table */
      dataDescriptionId = ConfigTab->rows[iConfig]->dataDescriptionId[iDD];
      for (jDD=0; jDD<DataDescTab->nrows; jDD++) {
	if (DataDescTab->rows[jDD]->dataDescriptionId==dataDescriptionId) break;
      }
      if (jDD>=DataDescTab->nrows) return;  /* Shouldn't need this */
      
      /* Now find spectralWindow */
      spectralWindowId = DataDescTab->rows[jDD]->spectralWindowId;
      for (jSW=0; jSW<SpectralWindowTab->nrows; jSW++) {
	if (SpectralWindowTab->rows[jSW]->spectralWindowId==spectralWindowId) break;
      }
      if (jSW>=SpectralWindowTab->nrows) return;  /* Shouldn't need this */
      
      /* Finally give spectral window info */
      if (SpectralWindowTab->rows[jSW]->bandcode) 
	strncpy (bcode, SpectralWindowTab->rows[jSW]->bandcode,9);
      else
	strcpy (bcode, "Unknown");
      sprintf(line,"     SpWin= %2d Freq= %7.3lf GHz, %d chan of BW=%9.3lf kHz, tot BW=%8.3lf MHz, %s Band=%s", 
	      SpectralWindowTab->rows[jSW]->spectralWindowId, 
	      SpectralWindowTab->rows[jSW]->refFreq*1.0e-9,
	      SpectralWindowTab->rows[jSW]->numChan,
	      SpectralWindowTab->rows[jSW]->chanWidth*1.0e-3,
	      SpectralWindowTab->rows[jSW]->totBandwidth*1.0e-6,
	      SpectralWindowTab->rows[jSW]->netSideband,
	      bcode);
      ObitPrinterWrite (myPrint, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    } /* end loop over DataDescriptions */
    /* ALMA specific stuff */
    if (isALMA) {
      if (ConfigTab->rows[iConfig]->numAtmPhaseCorrection==1) {
	if (ConfigTab->rows[iConfig]->atmPhaseCorrection[0]==ASDMAtmPhCorr_AP_CORRECTED)
	  sprintf(line,"       ALMA Data with atmospheric correction");
	if (ConfigTab->rows[iConfig]->atmPhaseCorrection[0]==ASDMAtmPhCorr_AP_UNCORRECTED)
	  sprintf(line,"       ALMA Data without atmospheric correction");
	ObitPrinterWrite (myPrint, line, &quit, err);
	if (quit) goto Quit;
      } else if (ConfigTab->rows[iConfig]->numAtmPhaseCorrection==2) {
	sprintf(line,"       ALMA Data with and without atmospheric correction");
	ObitPrinterWrite (myPrint, line, &quit, err);
	if (quit) goto Quit;
      }
      if (ConfigTab->rows[iConfig]->spectralType==ASDMSpecRes_CHANNEL_AVERAGE) {
	sprintf(line,"       Spec Res = CHANNEL_AVERAGE");
      } else if (ConfigTab->rows[iConfig]->spectralType==ASDMSpecRes_BASEBAND_WIDE) {
	sprintf(line,"       Spec Res = BASEBAND_WIDE");
      } else if (ConfigTab->rows[iConfig]->spectralType==ASDMSpecRes_FULL_RESOLUTION) {
	sprintf(line,"       Spec Res = FULL_RESOLUTION");
      }
    ObitPrinterWrite (myPrint, line, &quit, err);
    if (quit) goto Quit;
    } /* End ALMA */
     
  } /* end loop over configs */
  
  sprintf(line,"    ");
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (quit) goto Quit;
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Source info - new page */
  ObitPrinterSetTitle (myPrint, NULL, NULL, err);  /* No page titles at top */
  ObitPrinterNewPage (myPrint, &quit, err);
  if (quit) goto Quit;
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    
  /* Titles */
  sprintf (Title1, "Summary of ASDM in %s",DataRoot);
  sprintf (Title2, "  Scan Summary");
  ObitPrinterWrite (myPrint, Title1, &quit, err);
  if (quit) goto Quit;
  ObitPrinterWrite (myPrint, Title2, &quit, err);
  if (quit) goto Quit;
  ObitPrinterSetTitle (myPrint, Title1, Title2, err);  /* No page titles at top */
  if (quit) goto Quit;
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Loop over scans */

  /* Add intent information from ASDM */
  ScanTab = SDMData->ScanTab;
  MainTab = SDMData->MainTab;
  for (iScan=0; iScan<ScanTab->nrows; iScan++) {

    /* Get configID from MainTab */
    ScanID = ScanTab->rows[iScan]->scanNumber;
    for (iMain=0; iMain<MainTab->nrows; iMain++) {
      if (MainTab->rows[iMain]->scanNumber==ScanID) break;
    }
    if (iMain>=MainTab->nrows) return;  /* Shouldn't need this */
    configID = MainTab->rows[iMain]->configDescriptionId;
    
    /* Timerange in human form */
    day2dhms(ScanTab->rows[iScan]->startTime-refJD, begString);
    day2dhms(ScanTab->rows[iScan]->endTime-refJD,   endString);

    /* Get Calcode */
    calCode[0] = ' '; calCode[1] = 0;
    for (iSource=0; iSource<SourceArray->nsou; iSource++) {
      if (!strcmp(SourceArray->sou[iSource]->fieldName, 
		  ScanTab->rows[iScan]->sourceName)) {
	g_snprintf (calCode, 20, "%s", SourceArray->sou[iSource]->code);
	break;
      }
    }
    
    g_snprintf (line, 90, "Scan=%d config=%d Source=%s Code='%s' time= %s-%s", 
		ScanTab->rows[iScan]->scanNumber,configID,
		ScanTab->rows[iScan]->sourceName, calCode, begString, endString);
    ObitPrinterWrite (myPrint, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    
    for (iIntent=0; iIntent<ScanTab->rows[iScan]->numIntent; iIntent++) {
      g_snprintf (line, 80, "   Intent[%d]='%s'", 
		  iIntent+1,ScanTab->rows[iScan]->scanIntent[iIntent]);
      ObitPrinterWrite (myPrint, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    } /* end intent loop */
    /* ALMA specific */
    if (isALMA) {
      ii = 0;
      sprintf(line,"     ALMA Data types: ");
      while (ScanTab->rows[iScan]->calDataType[ii]) {
	strncat(line, ScanTab->rows[iScan]->calDataType[ii++], 1000);
	strncat(line, " ", 1000);
      }
      ObitPrinterWrite (myPrint, line, &quit, err);
      if (quit) goto Quit;
    }
  } /* End scan look */
  
  /* Done - Close Printer */
 Quit:
  ObitPrinterClose (myPrint, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Cleanup */
  AntArray    = ObitSDMDataKillAntArray(AntArray);
  SourceArray = ObitSDMDataKillSourceArray(SourceArray);

} /* end Summary */

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

