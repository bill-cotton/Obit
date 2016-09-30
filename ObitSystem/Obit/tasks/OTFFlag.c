/* $Id:  $  */
/* Obit task to Assure constant editing for EVLA OTF data       */
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

#include "ObitUVDesc.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitUVUtil.h"
#include "ObitUVEdit.h"
#include "ObitTableFG.h"
#include "ObitTableAN.h"
#include "ObitTableANUtil.h"
#include "ObitHistory.h"
#include "ObitData.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* OTFFlagIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void OTFFlagOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Digest inputs */
void digestInputs(ObitInfoList *myInput, ObitErr *err);
/* Get input data */
ObitUV* getInputData (ObitInfoList *myInput, ObitErr *err);
/* Flag all selected data */
void  FlagData(ObitUV* inData, ObitErr* err);
/* Write history */
void OTFFlagHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);

/* Program globals */
gchar *pgmName = "OTFFlag";       /* Program name */
gchar *infile  = "OTFFlag.in" ;   /* File with program inputs */
gchar *outfile = "OTFFlag.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL;  /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL;  /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   "Machine gun the lifeboats!"                                         */
/*   Obit task to flag survivors if too few left                          */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitTableFG  *FlagTab=NULL;
  olong        hiTab, flagTab=0;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = OTFFlagIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get input uvdata */
  inData = getInputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Make sure flagTab exists, assign if passed 0 */
  ObitInfoListGetTest(myInput, "flagTab",  &type, dim, &flagTab);
  FlagTab = newObitTableFGValue("tmpFG", (ObitData*)inData, &flagTab, OBIT_IO_ReadWrite, 
				err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  /* Open/close table */
  ObitTableFGOpen (FlagTab, OBIT_IO_ReadWrite, err);
  ObitTableFGClose (FlagTab, err);
  FlagTab = ObitTableFGUnref(FlagTab);
  /* Save for history */
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(myInput, "flagTab",  OBIT_long, dim, &flagTab);

  /* Get (new) highest table number */
  hiTab = ObitTableListGetHigh (inData->tableList, "AIPS FG")+1;
  /* Write flags to temporary, highest numbered table */
  ObitInfoListAlwaysPut (inData->info , "flagTmp", OBIT_long,  dim, &hiTab);

   /* Do editing */
  FlagData(inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Copy flags to flagTab, delete temp highest */
  ObitUVEditAppendFG (inData, flagTab, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Write history */
  OTFFlagHistory (myInput, inData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  FlagTab   = ObitUnref(FlagTab);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput  = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem  = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* OTFFlagIn (int argc, char **argv, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Parse control info from command line                                  */
/*   Input:                                                               */
/*      argc   Number of arguments from command line                      */
/*      argv   Array of strings from command line                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   return  ObitInfoList with defaults/parsed values                     */
/*----------------------------------------------------------------------- */
{
  olong ax;
  gchar *arg;
  gboolean init=FALSE;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint    itemp, i, j, k;
  ObitInfoList* list=NULL;
  gchar *routine = "OTFFlagIn";

  /* error checks */
  if (err->error) return list;

  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* Initialize output */
  myOutput = defaultOutputs(err);
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  /* command line arguments */
  if (argc<=1) Usage(); /* must have arguments */
  /* parse command line */
  for (ax=1; ax<argc; ax++) {

    arg = argv[ax];
    if (strcmp(arg, "-input") == 0){ /* input parameters */
      infile = argv[++ax];
      /* parse input file */
      ObitParserParse (infile, list, err);
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
      
    } else if (strcmp(arg, "-inSeq") == 0) { /* AIPS sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "inSeq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-inDisk") == 0) { /* input disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "inDisk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-DataType") == 0) { /* Image type AIPS or FITS */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "DataType", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-inName") == 0) { /* AIPS inName*/
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inName", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-inClass") == 0) { /* AIPS inClass*/
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inClass", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-inFile") == 0) { /*inFile */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inFile", OBIT_string, dim, strTemp);
      
   } else { /* unknown argument */
      Usage();
    }
    if (err->error) Obit_traceback_val (err, routine, "GetInput", list);
  } /* end parsing input arguments */
  
  /* Read defaults if no file specified */
  if (!init) ObitParserParse (infile, list, err);

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

  return list;
} /* end OTFFlagIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: OTFFlag -input file -output ofile [args]\n");
    fprintf(stderr, "OTFFlag: Assure constant flagging in scan\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def OTFFlag.in\n");
    fprintf(stderr, "  -output parameter file, def OTFFlag.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input image\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2\n");
    fprintf(stderr, "  -inFile input FITS uvdata file\n");
    fprintf(stderr, "  -inName input AIPS uvdata file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input (AIPS or FITS) disk number (1-rel) \n");
    /*/exit(1);  bail out */
  }/* end Usage */

/*----------------------------------------------------------------------- */
/*  Create default input ObitInfoList                                     */
/*  Note: Other parameters may be passed through the input text file      */
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
/*     inFile    Str [?]    input FITS image file name [def "Image.fits"] */
/*     inName    Str [12]   input AIPS image name  [no def]               */
/*     inClass   Str [6]    input AIPS image class  [no def]              */
/*     inSeq     Int        input AIPS image sequence no  [no def]        */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat farray[5];
  gboolean btemp;
  ObitInfoList *out = newObitInfoList();
  /*gchar *routine = "defaultInputs";*/

  /* error checks */
  if (err->error) return out;

  /* add parser items */
  /* Program number */
  dim[0] = 1; dim[1] = 1;
  itemp = 1;
  ObitInfoListAlwaysPut (out, "pgmNumber", OBIT_oint, dim, &itemp);

  /* Default FITS directories - same directory */
  dim[0] = 1; dim[1] = 1;
  itemp = 0; /* number of FITS directories */
  ObitInfoListAlwaysPut (out, "nFITS", OBIT_oint, dim, &itemp);

  /* AIPS user number */
  dim[0] = 1; dim[1] = 1;
  itemp = 2;
  ObitInfoListAlwaysPut (out, "AIPSuser", OBIT_oint, dim, &itemp);

  /* Default AIPS directories */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; /* number of AIPS directories */
  ObitInfoListAlwaysPut (out, "nAIPS", OBIT_oint, dim, &itemp);

  /* Default type "FITS" */
  strTemp = "FITS";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "DataType", OBIT_string, dim, strTemp);

  /* input FITS file name */
  strTemp = "OTFFlag.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "inFile", OBIT_string, dim, strTemp);

  /* input AIPS file name */
  strTemp = "OTFFlagName";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "inName", OBIT_string, dim, strTemp);

  /* input AIPS file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "inClass", OBIT_string, dim, strTemp);

  /* AIPS sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListAlwaysPut (out, "inSeq", OBIT_oint, dim, &itemp);

  /* AIPS or FITS disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListAlwaysPut (out, "inDisk", OBIT_oint, dim, &itemp);

  /* Timerange in days */
  dim[0] = 2;dim[1] = 1;
  farray[0] = -1.0e20; farray[1] = 1.0e20;
  ObitInfoListAlwaysPut (out, "timeRange", OBIT_float, dim, farray);

  /*  Apply calibration/selection?, def=True */
  dim[0] = 1; dim[1] = 1;
  btemp = TRUE;
  ObitInfoListAlwaysPut (out, "doCalSelect", OBIT_bool, dim, &btemp);

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
  ObitInfoList *out = newObitInfoList();
  /*  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
      ofloat ftemp;
      gchar *routine = "defaultOutputs";"*/

  /* error checks */
  if (err->error) return out;

  /* add parser items - nothing */
  return out;
} /* end defaultOutputs */

/*----------------------------------------------------------------------- */
/*  Digest inputs                                                         */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void digestInputs(ObitInfoList *myInput, ObitErr *err)
{
  /*ObitInfoType type;*/
  /*gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};*/
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

  /* Initialize Threading */
  ObitThreadInit (myInput);

} /* end digestInputs */

/*----------------------------------------------------------------------- */
/*  Get input data                                                        */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitUV with input data                                           */
/*----------------------------------------------------------------------- */
ObitUV* getInputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitUV       *inData = NULL;
  olong        nvis=1000;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "souCode", "Qual", "timeRange",  "FreqID", "subA", "Antennas", 
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer", 
    "doPol", "PDVer", "flagTab", 
   NULL};
  gchar *routine = "getInputData";

  /* error checks */
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Build basic input UV data Object */
  inData = ObitUVFromFileInfo ("in", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Set buffer size */
  nvis = 1000;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO",  OBIT_long, dim,  &nvis);

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Ensure inData fully instantiated and OK and selector set */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Flag data product for scan if previously there is flagged and         */
/*  unflagged data                                                        */
/*   Input:                                                               */
/*      inData    ObitUV to flag                                          */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void  FlagData(ObitUV* inData, ObitErr* err)
{
  ObitIOCode iretCode;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        flagTab=0;
  ObitTableFG  *FlagTab=NULL;
  ObitTableFGRow *row=NULL;
  ObitUVDesc *inDesc;
  ObitAntennaList *AntList  = NULL;
  ObitTableAN     *ANTable  = NULL;
  olong i, j, ant1, ant2, lastSubA, iBL, jBL, iVis, iIF, iStok, iChan, iFGRow=-1;
  olong nchan, nif, nstok, ncorr, incs, incf, incif, soff;
  olong maxAnt, numBL, numIF;
  olong subA, numPCal, numOrb, ver;
  ollong count=0, total = 0;
  ofloat *Buffer, sec, lastSou, timeBeg, timeEnd;
  gboolean *OKFlags=NULL, *badFlags=NULL;
  gchar reason[25];
  struct tm *lp;
  time_t clock;
  gchar *routine = "FlagData";
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inData));

  /* Get temporary output flag table */
  ObitInfoListGetTest(inData->info, "flagTmp",  &type, dim, &flagTab);
  FlagTab = newObitTableFGValue("outFG", (ObitData*)inData, &flagTab, OBIT_IO_ReadWrite, 
				err);
  ObitTableFGOpen (FlagTab, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Create Row */
  row = newObitTableFGRow (FlagTab);
  
  /* Attach row to output buffer */
  ObitTableFGSetRow (FlagTab, row, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Initialize flag row */
  row->SourID  = 0; 
  row->SubA    = 0; 
  row->freqID  = 0; 
  row->ants[0] =  row->ants[1] = 0; 
  row->TimeRange[0] = -1.0e20; 
  row->TimeRange[1] =  1.0e20; 
  row->ifs[0]    = row->ifs[1]    = 0; 
  row->chans[0]  = row->chans[1]  = 0; 
  row->pFlags[0] = 1<<0 | 1<<1 | 1<<2 | 1<<3; 
  row->pFlags[0] = row->pFlags[1] = row->pFlags[2] = row->pFlags[3] = 0;
  /* Reason includes time/date */
  /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);
  /* Convert to  broken-down time. */
  lp = localtime (&clock);
  lp->tm_mon++; /* For some bizzare reason, month is 0-rel */
  if (lp->tm_year<1000)  lp->tm_year += 1900; /* full year */
  sec = (ofloat)lp->tm_sec;
  g_snprintf (reason, 25, "OTFFlag %d/%d/%d %d:%d:%f", 
	      lp->tm_year, lp->tm_mon, lp->tm_mday, 
	      lp->tm_hour, lp->tm_min, sec);
  row->reason    = reason; /* Unique string */

  /* Get maximum antenna number */
  /* Convert AN table into AntennaList */
  ver = 1;  numPCal =  numOrb  = numIF = 0;
  ANTable = newObitTableANValue (inData->name, (ObitData*)inData, &ver, 
				 numIF, numOrb, numPCal, OBIT_IO_ReadOnly, err);
  if (ANTable) AntList = ObitTableANGetList (ANTable, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  maxAnt = 0;
  for (i=0; i<AntList->number; i++) maxAnt = MAX (maxAnt, AntList->ANlist[i]->AntID);
  ANTable = ObitTableANUnref(ANTable);
  AntList = ObitAntennaListUnref(AntList);
  numBL = maxAnt*maxAnt; /* inefficient but simple */
 
  /* Open visibility */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  inDesc  = inData->myDesc;  /* Get descriptor */

  /* Set up for parsing data */
  lastSou = -1;
  timeBeg = timeEnd = -1.0e20;
  lastSubA = 0;
  /* get data product increments */
  ncorr = inDesc->ncorr;
  nchan = inDesc->inaxes[inDesc->jlocf];
  incf  = inDesc->incf/3;
  if (inDesc->jlocf>=0) {
    nif = inDesc->inaxes[inDesc->jlocif];
    incif = inDesc->incif/3;
  } else {nif = 1; incif = 1;}
  if (inDesc->jlocs>=0) {
    nstok = inDesc->inaxes[inDesc->jlocs];
    incs  = inDesc->incs/3;
  } else {nstok = 1; incs = 1;}

  /* Work arrays */
  OKFlags  = g_malloc0(ncorr*numBL*sizeof(gboolean));
  badFlags = g_malloc0(ncorr*numBL*sizeof(gboolean));
  for (i=0; i<ncorr*numBL; i++) OKFlags[i]  = FALSE;
  for (i=0; i<ncorr*numBL; i++) badFlags[i] = FALSE;

  /* Loop over Data */
  while (1) {
    /* Read next buffer */
    iretCode = ObitUVReadSelect (inData, inData->buffer, err);
    if (err->error) goto cleanup;
    if (iretCode==OBIT_IO_EOF) break;  /* Done? */

    /* Loop over data in buffer */
    Buffer = inData->buffer;
    if (timeBeg<-1.0e5) {
      timeBeg = Buffer[inDesc->iloct];
      timeEnd = Buffer[inDesc->iloct];
    }
    for (iVis=0; iVis<inDesc->numVisBuff; iVis++) {
      ObitUVDescGetAnts(inData->myDesc, Buffer, &ant1, &ant2, &subA);
      iBL = (ant1-1)*maxAnt + ant2-1;
      /* New scan? i.e. new Source */
      if ((inDesc->ilocsu>=0) && (Buffer[inDesc->ilocsu] != lastSou)) {
	/* See if flagging needed */
	/* In case something bad */
	if (inDesc->ilocsu>=0) row->SourID = (olong)(lastSou+0.5); 
	/*if (inDesc->ilocfq>=0) row->freqID = (olong)(Buffer[inDesc->ilocfq]+0.5);  */
	row->TimeRange[0] = timeBeg-1.0e-6; /* +/-100 msec */
	row->TimeRange[1] = timeEnd+1.0e-6; 
	row->SubA         = lastSubA; 

	/* Loop over baselines */
	for (jBL=0; jBL<numBL; jBL++) {
	  /* Loop over data products */
	  for (i=0; i<ncorr; i++) {
	    total++;   /* Total possibilities */
	    if (OKFlags[i+jBL*ncorr] && badFlags[i+jBL*ncorr]) {
	      /* Data to flag */
	      ant1 = jBL/maxAnt; ant2 = jBL-ant1*maxAnt;
	      row->ants[0] = ant1+1; row->ants[1] = ant2+1;
	      /* Which data product is this? */
	      j = i;
	      iIF   = j/incif;   j -= iIF*incif;
	      iChan = j/incf;    j -= iChan*incf;
	      iStok = j/incs;
	      row->ifs[0]    = row->ifs[1]    = iIF+1; 
	      row->chans[0]  = row->chans[1] = iChan+1; 
	      row->pFlags[0] = 1<<iStok;
	      /* Write flag */
	      iFGRow = -1;
	      ObitTableFGWriteRow (FlagTab, iFGRow, row, err);
	      if (err->error) goto cleanup;
	      count++;  /* count flags  */
	    } /* end flag this one */
	  } /* end data products loop */
	} /* end baseline loop */
	/* Set up for next scan */
	timeBeg = MAX (timeBeg, Buffer[inDesc->iloct]);
	timeEnd = MAX (timeEnd, Buffer[inDesc->iloct]);
	lastSubA = subA;
	lastSou = Buffer[inDesc->ilocsu]; /* Save source */
	for (i=0; i<ncorr*numBL; i++) OKFlags[i]  = FALSE;
	for (i=0; i<ncorr*numBL; i++) badFlags[i] = FALSE;
      } /* end if new scan */

      /* Check each data product */
      timeEnd = MAX (timeEnd, Buffer[inDesc->iloct]);
      soff = inDesc->nrparm+2;
      for (i=0; i<ncorr; i++) {
	if (Buffer[soff]>0.0) OKFlags[i+iBL*ncorr]  = TRUE;
	else                  badFlags[i+iBL*ncorr] = TRUE;
	soff += 3;
      }
      Buffer += inDesc->lrec;  /* Update buffer pointer */
    } /* end loop over record */
  } /* end loop over data file */

  /* See if flagging needed in last scan */
  /* In case something bad */
  if (inDesc->ilocsu>=0) row->SourID = (olong)(lastSou+0.5); 
  /*if (inDesc->ilocfq>=0) row->freqID = (olong)(Buffer[inDesc->ilocfq]+0.5); */
  row->TimeRange[0] = timeBeg-1.0e-6; /* +/-100 msec */
  row->TimeRange[1] = timeEnd+1.0e-6; 
  row->SubA         = lastSubA; 
  
  /* Loop over baselines */
  for (jBL=0; jBL<numBL; jBL++) {
    /* Loop over data products */
    for (i=0; i<ncorr; i++) {
      total++;   /* Total possibilities */
      if (OKFlags[i+jBL*ncorr] && badFlags[i+jBL*ncorr]) {
	/* Data to flag */
	ant1 = jBL/numBL; ant2 = jBL-ant1*numBL;
	row->ants[0] = ant1+1; row->ants[1] = ant2+1;
	/* Which data product is this? */
	iIF   = i/incif;
	iChan = i/incf;
	iStok = i/incs;
	row->ifs[0]    = row->ifs[1]    = iIF+1; 
	row->chans[0]  = row->chans[1] = iChan+1; 
	row->pFlags[0] = 1<<iStok;
	/* Write flag */
	iFGRow = -1;
	ObitTableFGWriteRow (FlagTab, iFGRow, row, err);
	if (err->error) goto cleanup;
	count++;  /* count flags  */
      } /* end flag this one */
    } /* end data products loop */
  } /* end baseline loop */
  
   /* Give report */
  if (err->prtLv>=1) 
    Obit_log_error(err, OBIT_InfoErr, "%s: flag %8.1lf of %8.1lf vis/IFs= %6.3lf percent",
		   routine, (odouble)count, (odouble)total, 
		   100.0*(odouble)count/((odouble)total));
 /* close uv file */
 cleanup:
  ObitUVClose (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Close output table */
  ObitTableFGClose (FlagTab, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Cleanup */
  FlagTab = ObitTableFGUnref(FlagTab);
  row     = ObitTableFGRowUnref(row);
  if (OKFlags)  g_free(OKFlags);
  if (badFlags) g_free(badFlags);
} /* end FlagData */

/*----------------------------------------------------------------------- */
/*  Write History for OTFFlag                                             */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to update history                                */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void OTFFlagHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "FreqID", "Qual", "Sources", "timeRange", "subA", "Antennas", 
    "doCalib", "gainUse", "doPol", "PDVer", "flagVer", "doBand", "BPVer", 
    "flagTab", 
    NULL};
  gchar *routine = "OTFFlagHistory";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

   /* Do history */
  outHistory = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadWrite, err);

  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  outHistory = ObitHistoryUnref(outHistory); /* cleanup */
 
} /* end OTFFlagHistory  */

