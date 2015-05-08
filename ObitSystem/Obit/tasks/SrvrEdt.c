/* $Id$  */
/* Obit task to Flag selected UV-data       */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014-2015                                          */
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
#include "ObitHistory.h"
#include "ObitData.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* SrvrEdtIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void SrvrEdtOut (ObitInfoList* outList, ObitErr *err);
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
void SrvrEdtHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);

/* Program globals */
gchar *pgmName = "SrvrEdt";       /* Program name */
gchar *infile  = "SrvrEdt.in" ;   /* File with program inputs */
gchar *outfile = "SrvrEdt.out";   /* File to contain program outputs */
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
  myInput = SrvrEdtIn (argc, argv, err);
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
  SrvrEdtHistory (myInput, inData, err); 
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

ObitInfoList* SrvrEdtIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "SrvrEdtIn";

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
} /* end SrvrEdtIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: SrvrEdt -input file -output ofile [args]\n");
    fprintf(stderr, "SrvrEdt: Survivor flagging\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def SrvrEdt.in\n");
    fprintf(stderr, "  -output parameter file, def SrvrEdt.out\n");
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
  strTemp = "SrvrEdt.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "inFile", OBIT_string, dim, strTemp);

  /* input AIPS file name */
  strTemp = "SrvrEdtName";
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

  /* min OK fractions (IF, record) */
  dim[0] = 5;dim[1] = 1;
  farray[0] = 0.10; farray[1] = 0.10; farray[2]=0; farray[3]=0; farray[4]=0; 
  ObitInfoListAlwaysPut (out, "minOK", OBIT_float, dim, farray);

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
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

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
    "doPol", "PDVer", "minOK", "flagTab", 
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
/*  Flag all selected data                                                */
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
  olong iVis, iif, istok, ichan, iFGRow=-1;
  olong nchan, nif, nstok, incs, incf, incif, soff, ifoff, choff;
  olong recGood, recBad, *IFGood=NULL, *IFBad=NULL;
  olong ant1, ant2, subA;
  ollong count=0, total = 0;
  ofloat *Buffer, fracOK,  minOK[5]={0.1,0.1,0.,0.,0.};
  ofloat sec;
  gchar reason[25];
  struct tm *lp;
  time_t clock;
  gchar *routine = "FlagData";
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inData));

  /* Control parameters */
  ObitInfoListGetTest(inData->info, "minOK",  &type, dim, minOK);
 
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
  row->ants[0] = 0; 
  row->ants[1] = 0; 
  row->TimeRange[0] = -1.0e20; 
  row->TimeRange[1] =  1.0e20; 
  row->ifs[0]    = 1; 
  row->ifs[1]    = 0; 
  row->chans[0]  = 1; 
  row->chans[1]  = 0; 
  row->pFlags[0] = 1<<0 | 1<<1 | 1<<2 | 1<<3; 
  row->pFlags[1] = row->pFlags[2] = row->pFlags[3] = 0;
  /* Reason includes time/date */
  /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);
  /* Convert to  broken-down time. */
  lp = localtime (&clock);
  lp->tm_mon++; /* For some bizzare reason, month is 0-rel */
  if (lp->tm_year<1000)  lp->tm_year += 1900; /* full year */
  sec = (ofloat)lp->tm_sec;
  g_snprintf (reason, 25, "SrvrEdt %d/%d/%d %d:%d:%f", 
	      lp->tm_year, lp->tm_mon, lp->tm_mday, 
	      lp->tm_hour, lp->tm_min, sec);
  row->reason    = reason; /* Unique string */
  
  /* Open visibility */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  inDesc  = inData->myDesc;  /* Get descriptor */

  /* get visibility increments */
  incs  = inDesc->incs;
  incf  = inDesc->incf;
  incif = inDesc->incif;

  /* Set up for parsing data */
  nchan = inDesc->inaxes[inDesc->jlocf];
  incf  = inDesc->incf;
  if (inDesc->jlocf>=0) {
    nif = inDesc->inaxes[inDesc->jlocif];
    incif = inDesc->incif;
  } else {nif = 1; incif = 1;}
  if (inDesc->jlocs>=0) {
    nstok = inDesc->inaxes[inDesc->jlocs];
    incs  = inDesc->incs;
  } else {nstok = 1; incs = 1;}

  /* Work arrays */
  IFGood = g_malloc0(nif*sizeof(olong));
  IFBad  = g_malloc0(nif*sizeof(olong));

  /* Loop over Data */
  while (1) {
    /* Read next buffer */
    iretCode = ObitUVReadSelect (inData, inData->buffer, err);
    if (err->error) goto cleanup;
    if (iretCode==OBIT_IO_EOF) break;  /* Done? */

    /* Loop over data in buffer */
    Buffer = inData->buffer;
    for (iVis=0; iVis<inDesc->numVisBuff; iVis++) {
      /* Count good, bad and ugly */
      recGood = recBad = 0;
      for (iif=0; iif<nif; iif++) IFGood[iif] = IFBad[iif] = 0;
      /* loop over IF */
      for (iif=0; iif<nif; iif++) {
	ifoff = inDesc->nrparm + iif * incif;
	/* Loop over frequency channel */
	for (ichan=0; ichan<nchan; ichan++) {
	  choff = ifoff + ichan * incf;
	  soff  = choff;
	  /* Loop over polarization */
	  for (istok=0; istok<nstok; istok++) {
	    /* If this OK? */
	    if (Buffer[soff+2]>0.0) {
	      recGood++;
	      IFGood[iif]++;
	    } else {  /* flagged */
	      recBad++;
	      IFBad[iif]++;
	    }
	    soff += incs;
	  } /* end loop over stokes */
	} /* end loop over channel */
      } /* end loop over IF */
 
      total += nif;  /* count of total record/IFs */
      /* In case something bad */
      if (inDesc->ilocsu>=0) row->SourID = (olong)(Buffer[inDesc->ilocsu]+0.5); 
      if (inDesc->ilocfq>=0) row->freqID = (olong)(Buffer[inDesc->ilocfq]+0.5); 
      row->TimeRange[0] = Buffer[inDesc->iloct]-1.0e-6; /* +/-100 msec */
      row->TimeRange[1] = Buffer[inDesc->iloct]+1.0e-6; 
      ObitUVDescGetAnts(inDesc, Buffer, &ant1, &ant2, &subA);
      row->ants[0]      = ant1; 
      row->ants[1]      = ant2; 
      row->SubA         = subA; 
      row->ifs[0]       = 1; 
      row->ifs[1]       = 0; 

      /* Record flagging trumps IF */
      if (recGood+recBad>0) fracOK = recGood / ((ofloat)(recGood+recBad));
      else fracOK = 1.0;  /* Don't bother */
      if ((fracOK<minOK[1]) && (recGood>0)) {  /* axe it? */
	/* Write flags */
	iFGRow = -1;
	ObitTableFGWriteRow (FlagTab, iFGRow, row, err);
	if (err->error) goto cleanup;
	count += nif;  /* count flags */
      } else { /* Check for bad IFs */
	for (iif=0; iif<nif; iif++) {
	  if (IFGood[iif]+IFBad[iif]>0) 
	    fracOK = IFGood[iif] / ((ofloat)(IFGood[iif]+IFBad[iif]));
	  else fracOK = 1.0;  /* Don't bother */
	  if ((fracOK<minOK[0]) && (IFGood[iif]>0)) {  /* axe it? */
	    /* Write flags */
	    row->ifs[0] = row->ifs[1] = iif+1; 
	    iFGRow = -1;
	    ObitTableFGWriteRow (FlagTab, iFGRow, row, err);
	    if (err->error) goto cleanup;
	    count++;  /* count flags - need total as well */
	  }
	}
      } /* end IF check */
      Buffer += inDesc->lrec;  /* Update buffer pointer */
    } /* end loop over record */
  } /* end loop over data file */

   /* Give report */
  if (err->prtLv>=1) 
    Obit_log_error(err, OBIT_InfoErr, "%s: flag %8.1lf of %8.1lf vis/IFs= %5.1lf percent",
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
  if (IFBad)  g_free(IFBad);
  if (IFGood) g_free(IFGood);
} /* end FlagData */

/*----------------------------------------------------------------------- */
/*  Write History for SrvrEdt                                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to update history                                */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void SrvrEdtHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "FreqID", "Qual", "Sources", "timeRange", "subA", "Antennas", 
    "doCalib", "gainUse", "doPol", "PDVer", "flagVer", "doBand", "BPVer", 
    "flagTab", "minOK", 
    NULL};
  gchar *routine = "SrvrEdtHistory";

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
 
} /* end SrvrEdtHistory  */

