/* $Id$  */
/* Obit Radio interferometry calibration software                     */
/* applies user-selected corrections to the calibration SN table      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2014                                          */
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

#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitUV.h"
#include "ObitTableSN.h"

/* Control information structure */
typedef struct { 
  /**  UV descriptor */
  ObitUVDesc *desc;;
  /** Selected Stokes */
  gchar Stokes[5];
  /** First Stokes */
  olong bStoke;
  /** Highest Stokes */
  olong eStoke;
  /** First IF */
  olong BIF;
  /** Highest IF */
  olong EIF;
   /**  blanking value */
  ofloat fblank;
 /**  SNCorParm */
  ofloat SNCorParm[30];
  /**  PhasParm */
  ofloat PhasParm[30];
  /**  XFER last phase Pol1 */
  ofloat XFERPhas1[50];
  /**  XFER last phase Pol2 */
  ofloat XFERPhas2[50];
}  ControlInfo; 
  
/* internal prototypes */
/* Get inputs */
ObitInfoList* SNCorIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void SNCorOut (ObitInfoList* outList, ObitErr *err);
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
/* Write history */
void SNCorHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Modify SN table */
void SNCorDoCor (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Modification functions */
void doAVRT (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doCLPA (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doCLPP (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doCLPD (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doCLPR (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doCLPW (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doXFER (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doZPHS (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doZRAT (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doZDEL (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doMULA (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doREFP (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doCPRT (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doCPSN (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doPCOP (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doPNEG (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doNORM (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);
void doRSET (ControlInfo *control, ObitTableSNRow *row, ObitErr *err);


/* Program globals */
gchar *pgmName = "SNCor";       /* Program name */
gchar *infile  = "SNCor.in" ;   /* File with program inputs */
gchar *outfile = "SNCor.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */

int main ( int argc, char **argv )
/*------------------------------------------------------------------------- */
/*   Obit task applies user-selected corrections to the calibration SN table */
/*------------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitUV       *inData=NULL;
  ObitSystem   *mySystem= NULL;
  ObitErr      *err= NULL;

  /* Startup - parse command line */
  err = newObitErr();
  myInput = SNCorIn (argc, argv, err);
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

  /* Modify */
   SNCorDoCor (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Write history */
  SNCorHistory (myInput, inData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUVUnref(inData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput  = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* SNCorIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "SNCorIn";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return list;

  /* Make default inputs InfoList */
  list = defaultInputs(err);
  myOutput = defaultOutputs(err);

  /* command line arguments */
  /* fprintf (stderr,"DEBUG arg %d %s\n",argc,argv[0]); DEBUG */
  if (argc<=1) Usage(); /* must have arguments */
  /* parse command line */
  for (ax=1; ax<argc; ax++) {

     /*fprintf (stderr,"DEBUG next arg %s %s\n",argv[ax],argv[ax+1]); DEBUG */
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

  /* Initialize output */
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  return list;
} /* end SNCorIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: SNCor -input file -output ofile [args]\n");
    fprintf(stderr, "Obit task to modify an AIPS SN table \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def SNCor.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def SNCor.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input uv data\n");
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
/*     inFile    Str [?]    input FITS image file name                    */
/*     inName    Str [12]   input AIPS image name  [no def]               */
/*     inClass   Str [6]    input AIPS image class  [no def]              */
/*     inSeq     Int        input AIPS image sequence no  [no def]        */
/*     Sources   Str (16,1) Sources selected, blank = all                 */
/*     Stokes    Str (4)    Stokes parameter to image, def='    '         */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat farray[3];
  ObitInfoList *out = newObitInfoList();
  gchar *routine = "defaultInputs";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;

  /* add parser items */
  /* Program number */
  dim[0] = 1; dim[1] = 1;
  itemp = 1;
  ObitInfoListPut (out, "pgmNumber", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Default FITS directories - same directory */
  dim[0] = 1; dim[1] = 1;
  itemp = 0; /* number of FITS directories */
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

  /* input FITS file name */
  strTemp = "SNCor.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "SNCorName";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* AIPS sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "inSeq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* AIPS or FITS disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "inDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Sources selected, blank = all */
  strTemp = "                ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Sources", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* Timerange in days */
  dim[0] = 2;dim[1] = 1;
  farray[0] = -1.0e20; farray[1] = 1.0e20;
  ObitInfoListPut (out, "timeRange", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Subarray */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "subA", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

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
  g_assert(ObitErrIsA(err));
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
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doCalSelect;
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* Make sure doCalSelect set properly */
  doCalSelect = TRUE;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);
 
} /* end digestInputs */

/*----------------------------------------------------------------------- */
/*  Get input data                                                        */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitUV with SN table to modify                                   */
/*----------------------------------------------------------------------- */
ObitUV* getInputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitUV       *inData = NULL;
  ObitInfoType type;
  olong         Aseq, disk, cno, nvis=1000;
  gchar        *Type, *strTemp, inFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "UV";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "Stokes", "timeRange", "BIF", "EIF", "subA",
    "doCalSelect", "Antennas", "FreqID", "souCode", "Qual", 
     NULL};
  gchar *routine = "getInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input UV data Object */
  inData = newObitUV("input UV data");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
    /* input AIPS disk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    /* input AIPS name */
    if (ObitInfoListGetP(myInput, "inName", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aname, strTemp, 13);
    } else { /* Didn't find */
      strncpy (Aname, "No Name ", 13);
    } 
    Aname[12] = 0;
    /* input AIPS class */
    if  (ObitInfoListGetP(myInput, "inClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    Aclass[6] = 0;
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "inSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 want highest existing sequence */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "inSeq", OBIT_oint, dim, &Aseq);
    }

    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
    /* define object */
    ObitUVSetAIPS (inData, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    if (ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp)) {
      strncpy (inFile, strTemp, 128);
    } else { 
      strncpy (inFile, "No_Filename_Given", 128);
    }
    
    /* input FITS disk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);

    /* define object */
    ObitUVSetFITS (inData, nvis, disk, inFile,  err); 
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return inData;
  }


  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

 /* Ensure inData fully instantiated and OK */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Write History for SNCor                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void SNCorHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "souCode", "Qual", "EditStokes", "BIF", "EIF", 
    "FreqID", "timeRange",  "subA", "Antennas", 
    "solnVer", "corMode", "SNCParm", "PhasParm",
    NULL};
  gchar *routine = "SNCorHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
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
 
} /* end SNCorHistory  */

/*----------------------------------------------------------------------- */
/*  Apply corrections to SN table                                         */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV with SN table to modify                          */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void SNCorDoCor (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  gchar corMode[5];
  oint numpol, numif;
  ObitInfoType type;
  olong SNver, inSNRow;
  ofloat t1, t2, fblank = ObitMagicF();
  odouble freq;
  olong i, itemp, nif, nstok;
  gint32 dim[MAXINFOELEMDIM];
  ObitTableSN *SNTable=NULL;
  ObitTableSNRow *inRow=NULL;
  ObitIOCode retCode;
  ControlInfo control;
  ObitUVSel *sel;
  gboolean wanted;
  gchar *routine = "SNCorDoCor";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Get Operation type */
  ObitInfoListGet(myInput, "corMode", &type, dim, corMode, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Get SN table */
  itemp = 0;
  ObitInfoListGetTest(myInput, "solnVer", &type, dim, &itemp);
  SNver = itemp;
  numpol = 0;
  numif  = 0;
  SNTable = newObitTableSNValue("SN table", (ObitData*)inData, &SNver, 
				OBIT_IO_ReadWrite, numpol, numif, err);
  /* Open */
  ObitTableSNOpen (SNTable, OBIT_IO_ReadWrite,  err);
  if (err->error) goto cleanup;

  /* Set row */
  inRow  = newObitTableSNRow (SNTable);

  /* Set control information */
  control.fblank = fblank;
  if (inData->myDesc->jlocif>=0) nif = inData->myDesc->inaxes[inData->myDesc->jlocif];
  else nif = 1;
  if (inData->myDesc->jlocs>=0) nstok = inData->myDesc->inaxes[inData->myDesc->jlocs];
  else nstok = 1;
  nstok = MIN (nstok, 2);
  if (inData->myDesc->jlocf>=0) freq = inData->myDesc->crval[inData->myDesc->jlocf];
  else freq = 1.0;
  control.desc = inData->myDesc;
  ObitInfoListGet(myInput,  "EditStokes", &type, dim, control.Stokes, err);
  if (err->error) goto cleanup;
  if (control.Stokes[0]=='R') {  /* R = first poln */
    control.bStoke = 1;
    control.eStoke = 1;
  } else if (control.Stokes[0]=='L') {  /* L  */
    if ((inData->myDesc->jlocs>=0) && 
	(inData->myDesc->crval[inData->myDesc->jlocs]<-1.1)) { /* Only L in data */
      control.bStoke = 1;
      control.eStoke = 1;
    } else { /* L second poln */
      control.bStoke = 2;
      control.eStoke = 2;
    }
  } else {   /* All Stokes */
    control.bStoke = 1;
    control.eStoke = nstok;
  }
  control.BIF = 1;
  ObitInfoListGetTest(myInput,  "BIF", &type, dim, &control.BIF);
  control.BIF = MAX (1, control.BIF);
  control.EIF = 0;
  ObitInfoListGetTest(myInput,  "EIF", &type, dim, &control.EIF);
  if (control.EIF<=0) control.EIF = nif;
  control.EIF = MAX (control.BIF, MIN (nif, control.EIF));
  for (i=0; i<30; i++) control.SNCorParm[i] = 0.0;
  ObitInfoListGetTest(myInput,  "SNCParm", &type, dim, control.SNCorParm);
  for (i=0; i<30; i++) control.PhasParm[i] = 0.0;
  ObitInfoListGetTest(myInput,  "PhasParm",  &type, dim, control.PhasParm);

  /* Open and close UV data to update selector */
  ObitUVOpen(inData, OBIT_IO_ReadCal, err);
  ObitUVClose(inData, err);
  if (err->error) goto cleanup;
  sel  = inData->mySel;

  /* Any corMode specific checks or processing */
  if (!strncmp (corMode, "PCOP", 4)) {
    /* Check reference IF */
    if ((control.SNCorParm[0]!=1.0) && (control.SNCorParm[0]!=2.0)) {
      Obit_log_error(err, OBIT_Error,
		     "%s: Unknown control parameter %f ", 
		     routine, control.SNCorParm[0]);
      return;
    }
  } else  if (!strncmp (corMode, "CPSN", 4)) {
    /* Check reference IF */
    if ((control.SNCorParm[0]<1.0) || (control.SNCorParm[0]>nif)) {
      Obit_log_error(err, OBIT_Error,
		     "%s: Invalid IF number %f ", 
		     routine, control.SNCorParm[0]);
      return;
    }
  } else  if (!strncmp (corMode, "CPRT", 4)) {
    /* Check reference IF */
    if ((control.SNCorParm[0]<1.0) || (control.SNCorParm[0]>nif)) {
      Obit_log_error(err, OBIT_Error,
		     "%s: Invalid IF number %f ", 
		     routine, control.SNCorParm[0]);
      return;
    }
  } else  if (!strncmp (corMode, "REFP", 4)) {
    /* Check reference IF */
    if ((control.SNCorParm[0]<1.0) || (control.SNCorParm[0]>nif)) {
      Obit_log_error(err, OBIT_Error,
		     "%s: Invalid IF number %f ", 
		     routine, control.SNCorParm[0]);
      return;
    }
  } else  if (!strncmp (corMode, "XFER", 4)) {
    /* Initialize last phases to 0 */
    for (i=0; i<50; i++) control.XFERPhas1[i] = control.XFERPhas2[i] = 0.0;

    /* Check upper and lower frequencies */
    if ((control.PhasParm[28]<100.0) || (control.PhasParm[29]==100.0)) {
      Obit_log_error(err, OBIT_Error,
		     "%s: Invalid upper/lower frequencies %f %f", 
		     routine, control.SNCorParm[29], control.SNCorParm[29]);
      return;
    }
  } else  if (!strncmp (corMode, "CLPR", 4)) {
    /* Convert milliHz to sec/sec */
    t1 = control.SNCorParm[0];
    t2 = control.SNCorParm[1];
    control.SNCorParm[0] = MIN (t1, t2) * 1.0e3 / freq;
    control.SNCorParm[1] = MAX (t1, t2) * 1.0e3 / freq;
  } else  if (!strncmp (corMode, "CLPD", 4)) {
    /* Convert nsec to sec */
    t1 = control.SNCorParm[0];
    t2 = control.SNCorParm[1];
    control.SNCorParm[0] = MIN (t1, t2) * 1.0e-9;
    control.SNCorParm[1] = MAX (t1, t2) * 1.0e-9;
  } else  if (!strncmp (corMode, "CLPP", 4)) {
    /* Convert phases to radians */
    t1 = control.SNCorParm[0];
    t2 = control.SNCorParm[1];
    control.SNCorParm[0] = MIN (t1, t2) * DG2RAD;
    control.SNCorParm[1] = MAX (t1, t2) * DG2RAD;
  }

  /* Loop over table */
  for (inSNRow=1; inSNRow<=SNTable->myDesc->nrow; inSNRow++) {
    retCode = ObitTableSNReadRow (SNTable, inSNRow, inRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    if (inRow->status==-1) continue;

    /* This one selected? */
    wanted = (inRow->Time>=sel->timeRange[0]) && (inRow->Time<=sel->timeRange[1]);
    wanted = wanted && ObitUVSelWantSour (sel, inRow->SourID);
    wanted = wanted && ObitUVSelWantAnt (sel, inRow->antNo);
    wanted = wanted && ((inRow->FreqID==sel->FreqID) || (sel->FreqID<=0));
    wanted = wanted && ((inRow->SubA==sel->SubA) || (sel->SubA<=0));
    if (!wanted) continue;

    /* Do operation */
    if (!strncmp (corMode, "AVRT", 4)) {
      doAVRT (&control, inRow, err);
    } else if (!strncmp (corMode, "CLPA", 4)) {
      doCLPA (&control, inRow, err);
    } else if (!strncmp (corMode, "CLPP", 4)) {
      doCLPP (&control, inRow, err);
    } else if (!strncmp (corMode, "CLPD", 4)) {
      doCLPD (&control, inRow, err);
    } else if (!strncmp (corMode, "CLPR", 4)) {
      doCLPR (&control, inRow, err);
    } else if (!strncmp (corMode, "CLPW", 4)) {
      doCLPW (&control, inRow, err);
    } else if (!strncmp (corMode, "XFER", 4)) {
      doXFER (&control, inRow, err);
    } else if (!strncmp (corMode, "ZPHS", 4)) {
      doZPHS (&control, inRow, err);
    } else if (!strncmp (corMode, "ZRAT", 4)) {
      doZRAT (&control, inRow, err);
    } else if (!strncmp (corMode, "ZDEL", 4)) {
      doZDEL (&control, inRow, err);
    } else if (!strncmp (corMode, "MULA", 4)) {
      doMULA (&control, inRow, err);
    } else if (!strncmp (corMode, "REFP", 4)) {
      doREFP (&control, inRow, err);
    } else if (!strncmp (corMode, "CPRT", 4)) {
      doCPRT (&control, inRow, err);
    } else if (!strncmp (corMode, "CPSN", 4)) {
      doCPSN (&control, inRow, err);
    } else if (!strncmp (corMode, "PCOP", 4)) {
      doPCOP (&control, inRow, err);
    } else if (!strncmp (corMode, "PNEG", 4)) {
      doPNEG (&control, inRow, err);
    } else if (!strncmp (corMode, "NORM", 4)) {
      doNORM (&control, inRow, err);
    } else if (!strncmp (corMode, "RSET", 4)) {
      doRSET (&control, inRow, err);
    } else {  /* Unknown */
      Obit_log_error(err, OBIT_Error,
		     "%s: Unknown corMode: %s", routine, corMode);
      return;
    } /* End of cal by mode */
    if (err->error) goto cleanup;

    /* Rewrite row */
    retCode = ObitTableSNWriteRow (SNTable, inSNRow, inRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  } /* end loop over table */

  /* Close table */
  ObitTableSNClose (SNTable, err);

  /* Cleanup */
 cleanup:
  SNTable = ObitTableSNUnref(SNTable);
  inRow = ObitTableSNRowUnref(inRow);
  if (err->error)	Obit_traceback_msg (err, routine, SNTable->name);

} /* end SNCorDoCor  */

/* Modification functions */
/*----------------------------------------------------------------------- */
/*  Apply AVRT corrections to SN table                                    */
/*  Average selected fringe rates in a record                             */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doAVRT (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat sum1, sum2, rate1, rate2, fblank;
  olong count1, count2;
  if (err->error) return;

  fblank = control->fblank;

  /* Loop over IF summing rates */
  sum1 = sum2 = 0.0;
  count1 = count2 = 0;
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  if ((row->Real1[jif]!=fblank) && (row->Imag1[jif]!=fblank) &&
	      (row->Rate1[jif]!=fblank) && (row->Weight1[jif]>0.0)) {
	    sum1 += row->Rate1[jif];
	    count1++;
	  }
	} else {       /* second pol */
	  if ((row->Real2[jif]!=fblank) && (row->Imag2[jif]!=fblank) &&
	      (row->Rate2[jif]!=fblank) && (row->Weight2[jif]>0.0)) {
	    sum2 += row->Rate2[jif];
	    count2++;
	  }
	}
    } /* end loop over Stokes */
  } /* end loop over IF */

  /* Average rates */
  if (count1>0) rate1 = sum1/count1;
  else rate1 = fblank;
  if (count2>0) rate2 = sum2/count2;
  else rate2 = fblank;

  /* Loop over IF replacing rate with average */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	   if ((row->Real1[jif]!=fblank) && (row->Imag1[jif]!=fblank) &&
	       (row->Rate1[jif]!=fblank) && (row->Weight1[jif]>0.0)) {
	    row->Rate1[jif]   = rate1;
	   } else {
	    row->Rate1[jif]   = fblank;
	   }
	} else {       /* second pol */
	   if ((row->Real2[jif]!=fblank) && (row->Imag2[jif]!=fblank) &&
	       (row->Rate2[jif]!=fblank) && (row->Weight2[jif]>0.0)) {
	    row->Rate2[jif]   = rate2;
	   } else {
	    row->Rate2[jif]   = fblank;
	   }
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doAVRT */

/*----------------------------------------------------------------------- */
/*  Apply CLPA corrections to SN table                                    */
/*  Flag amplitudes outside of the range SNCorParm[0] (min) to            */
/*  SNCorParm[1] (max).                                                   */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doCLPA (ControlInfo *control, ObitTableSNRow *row, 
	     ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat amax, amin, ampl, fblank;
  gboolean bad;
  if (err->error) return;

  fblank = control->fblank;
  amin   = control->SNCorParm[0];
  amax   = control->SNCorParm[1];

  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  bad = (row->Real1[jif]==fblank) || (row->Imag1[jif]==fblank);
	  if (bad) ampl = 1.0;
	  else ampl = sqrt (row->Real1[jif]*row->Real1[jif] +
			    row->Imag1[jif]*row->Imag1[jif]);
	  if (bad || (ampl<amin) || (ampl>amax)) {
	    row->Real1[jif]   = fblank;
	    row->Imag1[jif]   = fblank;
	    row->Delay1[jif]  = fblank;
	    row->Rate1[jif]   = fblank;
	    row->Weight1[jif] = 0.0;
	  } /* end blanked */
	} else {       /* second pol */
	  bad = (row->Real2[jif]==fblank) || (row->Imag2[jif]==fblank);
	  if (bad) ampl = 1.0;
	  else ampl = sqrt (row->Real2[jif]*row->Real2[jif] +
			    row->Imag2[jif]*row->Imag2[jif]);
	  if (bad || (ampl<amin) || (ampl>amax)) {
	    row->Real2[jif]   = fblank;
	    row->Imag2[jif]   = fblank;
	    row->Delay2[jif]  = fblank;
	    row->Rate2[jif]   = fblank;
	    row->Weight2[jif] = 0.0;
	  } /* end blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doCLPA */

/*----------------------------------------------------------------------- */
/*  Apply CLPP corrections to SN table                                    */
/*  Flags any gain entries with a phase outside of the range given by     */
/*  SNCorParm[0] (min) to SNCorParm[1] (max).                             */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doCLPP (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat pmax, pmin, phas, fblank;
  gboolean bad;
  if (err->error) return;

  fblank = control->fblank;
  pmin   = control->SNCorParm[0];
  pmax   = control->SNCorParm[1];

  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  bad = (row->Real1[jif]==fblank) || (row->Imag1[jif]==fblank);
	  if (bad) phas = 0.0;
	  else phas = atan2 (row->Real1[jif], row->Imag1[jif]+1.0e-20);
	  if (bad || (phas<pmin) || (phas>pmax)) {
	    row->Real1[jif]   = fblank;
	    row->Imag1[jif]   = fblank;
	    row->Delay1[jif]  = fblank;
	    row->Rate1[jif]   = fblank;
	    row->Weight1[jif] = 0.0;
	  } /* end blanked */
	} else {       /* second pol */
	  bad = (row->Real2[jif]==fblank) || (row->Imag2[jif]==fblank);
	  if (bad) phas = 0.0;
	  else phas = atan2 (row->Real2[jif], row->Imag2[jif]+1.0e-20);
	  if (bad || (phas<pmin) || (phas>pmax)) {
	    row->Real2[jif]   = fblank;
	    row->Imag2[jif]   = fblank;
	    row->Delay2[jif]  = fblank;
	    row->Rate2[jif]   = fblank;
	    row->Weight2[jif] = 0.0;
	  } /* end blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doCLPP */

/*----------------------------------------------------------------------- */
/*  Apply CLPD corrections to SN table                                    */
/*  Flags any  entries with a delay whose value in nanoseconds is outside */
/*  of the range given by SNCorParm[0] (min) to SNCorParm[1] (max)        */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doCLPD (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat dmax, dmin, fblank;
  if (err->error) return;

  fblank = control->fblank;
  dmin   = control->SNCorParm[0];
  dmax   = control->SNCorParm[1];

  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  if ((row->Real1[jif]==fblank) || (row->Weight1[jif]<=0.0) ||
	      (row->Delay1[jif]<dmin) || (row->Delay1[jif]>dmax)) {
	    row->Real1[jif]   = fblank;
	    row->Imag1[jif]   = fblank;
	    row->Delay1[jif]  = fblank;
	    row->Rate1[jif]   = fblank;
	    row->Weight1[jif] = 0.0;
	  } /* end blanked */
	} else {       /* second pol */
	  if ((row->Real2[jif]==fblank) || (row->Weight2[jif]<=0.0) ||
	      (row->Delay2[jif]<dmin) || (row->Delay2[jif]>dmax)) {
	    row->Real2[jif]   = fblank;
	    row->Imag2[jif]   = fblank;
	    row->Delay2[jif]  = fblank;
	    row->Rate2[jif]   = fblank;
	    row->Weight2[jif] = 0.0;
	  } /* end blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doCLPD */

/*----------------------------------------------------------------------- */
/*  Apply CLPR corrections to SN table                                    */
/*  Flags any  entries with a fringe rate whose value in milliHz is       */
/*  outside of the range given by SNCorParm[0] (min) to SNCorParm[1] (max)*/
/*  The given fringe rate in mHz is at the reference frequency.           */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doCLPR (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat rmax, rmin, fblank;
  if (err->error) return;

  fblank = control->fblank;
  rmin   = control->SNCorParm[0];
  rmax   = control->SNCorParm[1];

  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  if ((row->Real1[jif]==fblank) || (row->Weight1[jif]<=0.0) ||
	      (row->Rate1[jif]<rmin) || (row->Rate1[jif]>rmax)) {
	    row->Real1[jif]   = fblank;
	    row->Imag1[jif]   = fblank;
	    row->Delay1[jif]  = fblank;
	    row->Rate1[jif]   = fblank;
	    row->Weight1[jif] = 0.0;
	  } /* end blanked */
	} else {       /* second pol */
	  if ((row->Real2[jif]==fblank) || (row->Weight2[jif]<=0.0) ||
	      (row->Rate2[jif]<rmin) || (row->Rate2[jif]>rmax)) {
	    row->Real2[jif]   = fblank;
	    row->Imag2[jif]   = fblank;
	    row->Delay2[jif]  = fblank;
	    row->Rate2[jif]   = fblank;
	    row->Weight2[jif] = 0.0;
	  } /* end blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doCLPR */

/*----------------------------------------------------------------------- */
/*  Apply CLPW corrections to SN table                                    */
/*  Flags any gain entries whose weight lies outside the range            */
/*  SNCorParm[0] (min) to SNCorParm[1] (max).                             */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doCLPW (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat wmax, wmin, fblank;
  if (err->error) return;

  fblank = control->fblank;
  wmin   = control->SNCorParm[0];
  wmax   = control->SNCorParm[1];

  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  if ((row->Weight1[jif]<wmin) || (row->Weight1[jif]>wmax) ) {
	    row->Real1[jif]   = fblank;
	    row->Imag1[jif]   = fblank;
	    row->Delay1[jif]  = fblank;
	    row->Rate1[jif]   = fblank;
	    row->Weight1[jif] = 0.0;
	  } /* end blanked */
	} else {       /* second pol */
	  if ((row->Weight2[jif]<wmin) || (row->Weight2[jif]>wmax) ) {
	    row->Real2[jif]   = fblank;
	    row->Imag2[jif]   = fblank;
	    row->Delay2[jif]  = fblank;
	    row->Rate2[jif]   = fblank;
	    row->Weight2[jif] = 0.0;
	  } /* end blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doCLPW */

/*----------------------------------------------------------------------- */
/*  Apply XFER corrections to SN table                                    */
/*      This option transfers phases between frequencies and also sets    */
/*      antenna amplitudes.                                               */
/*  Antenna amplitudes: SNCorParm[30]:                                    */
/*      If left blank, the antenna amplitude gains will be unchanged.     */
/*      To manually set the antenna amplitude gains, enter                */
/*      the desired amplitudes strictly in order 1 to 28 into             */
/*      SNCorParm[0] to SNCorParm[27] and set SNCorParm[29] to -1.        */
/*      To set all antenna amplitude gains to 1.0 set SNCorParm[29]       */
/*      to +1.                                                            */
/*  Phase transfer parameters: PHASParm[30]                               */
/*      PhasParm[0] - PhasParm[27] are the instrumental phase offsets     */
/*      at the two frequencies for Antennas 1 to 28 (in radians).         */
/*      [They must be strictly in order 1 to 28.]                         */
/*      PhasParm[28] = The higher frequency in MHz.                       */
/*      PhasParm[29] = The lower frequency in MHz.                        */
/*  Adopted from the AIPSish SNCOR/XFER                                   */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doXFER (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, jif, kant;
  ofloat ampl, amp, phas, ratio, fblank;
  gboolean bad;
  if (err->error) return;

  fblank = control->fblank;
  kant = MIN (0, row->antNo-1);
  ratio = control->PhasParm[28] /  control->PhasParm[29];
  if (control->SNCorParm[29]<0.0) ampl = control->SNCorParm[kant];
  else if  (control->SNCorParm[29]>0.0) ampl = 1.0;
  else ampl = -1.0;

  /* First poln */
  if (control->bStoke==1) {
    /* Loop over IF */
    for (iif=control->BIF; iif<=control->EIF; iif++) {
      jif = iif - 1;
      bad = (row->Real1[jif]==fblank) || (row->Imag1[jif]==fblank);
      if (!bad) {
	if (ampl<0.0) amp = sqrt (row->Real1[jif]*row->Real1[jif] +
				  row->Imag1[jif]*row->Imag1[jif]);
	else amp = ampl;
	phas = atan2 (row->Real1[jif], row->Imag1[jif]+1.0e-20);
	/* Closest turn to last value */
	while (fabs(phas-control->XFERPhas1[kant])>=G_PI) {
	  if ((phas-control->XFERPhas1[kant])>0.0) phas -= 2.0*G_PI;
	  else phas += 2.0*G_PI;
	}
	/* Save last  for first IF */
	if (iif==control->BIF) control->XFERPhas1[kant] = phas;  
	/* Scale and offset phase */
	phas = ratio * phas + control->PhasParm[kant];
	row->Real1[jif]  = amp * cos(phas);
	row->Imag1[jif]  = amp * sin(phas);
      } /* end not blanked */
    } /* end IF loop */
  } /* end first poln */

  /* Second poln */
  if ((control->bStoke==2) || (control->eStoke==2))  {
    /* Loop over IF */
    for (iif=control->BIF; iif<=control->EIF; iif++) {
      jif = iif - 1;
      bad = (row->Real2[jif]==fblank) || (row->Imag2[jif]==fblank);
      if (!bad) {
	if (ampl<0.0) amp = sqrt (row->Real2[jif]*row->Real2[jif] +
				  row->Imag2[jif]*row->Imag2[jif]);
	else amp = ampl;
	phas = atan2 (row->Real2[jif], row->Imag2[jif]+1.0e-20);
	/* Closest turn to last value */
	while (fabs(phas-control->XFERPhas2[kant])>=G_PI) {
	  if ((phas-control->XFERPhas2[kant])>0.0) phas -= 2.0*G_PI;
	  else phas += 2.0*G_PI;
	}
	/* Save last  for first IF */
	if (iif==control->BIF) control->XFERPhas2[kant] = phas;  
	/* Scale and offset phase */
	phas = ratio * phas + control->PhasParm[kant];
	row->Real2[jif]  = amp * cos(phas);
	row->Imag2[jif]  = amp * sin(phas);
      } /* end not blanked */
    } /* end IF loop */
  } /* end second poln */

} /* end doXFER */

/*----------------------------------------------------------------------- */
/*  Apply ZPHS corrections to SN table                                    */
/*   Cause the phases of selected gaines to be set to zero.               */
/*   This is useful for the calibration of amplitudes only.               */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doZPHS (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat ampl, fblank;
  if (err->error) return;

  fblank = control->fblank;
  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  if ((row->Real1[jif]!=fblank) && (row->Imag1[jif]!=fblank) &&
	      (row->Weight1[jif]>0.0)) {
	    ampl = sqrt (row->Real1[jif]*row->Real1[jif] + 
			 row->Imag1[jif]*row->Imag1[jif]);
	    row->Real1[jif] = ampl;
	    row->Imag1[jif] = 0.0;
	  } /* end not blanked */
	} else {       /* second pol */
	  if ((row->Real2[jif]!=fblank) && (row->Imag2[jif]!=fblank) &&
	      (row->Weight2[jif]>0.0)) {
	    ampl = sqrt (row->Real2[jif]*row->Real2[jif] + 
			 row->Imag2[jif]*row->Imag2[jif]);
	    row->Real2[jif] = ampl;
	    row->Imag2[jif] = 0.0;
	  } /* end not blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doZPHS */

/*----------------------------------------------------------------------- */
/*  Apply ZRAT corrections to SN table                                    */
/*  Cause the residual fringe-rates of the selected solutions to be       */
/*  set to zero.                                                          */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doZRAT (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat fblank;
  if (err->error) return;

  fblank = control->fblank;
  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  if ((row->Real1[jif]!=fblank) && (row->Weight1[jif]>0.0)) {
	    row->Rate1[jif] = 0.0;
	  } /* end not blanked */
	} else {       /* second pol */
	  if ((row->Real2[jif]!=fblank) && (row->Weight2[jif]>0.0)) {
	    row->Rate2[jif] = 0.0;
	  } /* end not blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doZRAT */

/*----------------------------------------------------------------------- */
/*  Apply ZDEL corrections to SN table                                    */
/*  Cause the residual delays of the selected solutions to be set to zero.*/
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doZDEL (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat fblank;
  if (err->error) return;

  fblank = control->fblank;
  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  if ((row->Real1[jif]!=fblank) && (row->Weight1[jif]>0.0)) {
	    row->Delay1[jif] = 0.0;
	  } /* end not blanked */
	} else {       /* second pol */
	  if ((row->Real2[jif]!=fblank) && (row->Weight2[jif]>0.0)) {
	    row->Delay2[jif] = 0.0;
	  } /* end not blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doZDEL */

/*----------------------------------------------------------------------- */
/*  Apply MULA corrections to SN table                                    */
/*  Cause the amplitudes of the selected complex gains to be multiplied   */
/*  by SNCorParm[0].                                                      */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doMULA (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat fblank;
  if (err->error) return;

  fblank = control->fblank;
  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  if ((row->Real1[jif]!=fblank) && (row->Imag1[jif]!=fblank)) {
	    row->Real1[jif] *= control->SNCorParm[0];
	    row->Imag1[jif] *= control->SNCorParm[0];
	  } /* end not blanked */
	} else {       /* second pol */
	  if ((row->Real2[jif]!=fblank) && (row->Imag2[jif]!=fblank)) {
	    row->Real2[jif] *= control->SNCorParm[0];
	    row->Imag2[jif] *= control->SNCorParm[0];
	  } /* end not blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doMULA */

/*----------------------------------------------------------------------- */
/*  Apply REFP corrections to SN table                                    */
/*  Reference all the phases to the IF specified by SNCorParm[0].         */
/*  If the phase for this IF is flagged all phases will be flagged.       */
/*  Default = IF 1.                                                       */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doREFP (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif, kif;
  ofloat norm, rotR1=0.0, rotI1=0.0, rotR2=0.0, rotI2=0.0, tre, tim, fblank;
  gboolean bad1=FALSE, bad2=FALSE;
  if (err->error) return;

  fblank = control->fblank;
  kif = control->SNCorParm[0]-0.5;

  /* Phase rotations */
  if (control->bStoke==1) {
    bad1 = (row->Real1[kif]!=fblank) || (row->Imag1[kif]!=fblank) || 
      (row->Weight1[kif]<=0.0);
    if (!bad1) {
      rotR1 =  row->Real1[kif];
      rotI1 = -row->Imag1[kif];
      norm = rotR1*rotR1 + rotI1*rotI1;
      if (norm>1.0e-25) norm = 1.0 / sqrt(norm);
      else norm = 1.0;
      rotR1 *= norm;
      rotI1 *= norm;
    }
  }
  if ((control->eStoke==2) || (control->bStoke==2)) {
    bad2 = (row->Real2[kif]!=fblank) || (row->Imag2[kif]!=fblank) || 
      (row->Weight2[kif]<=0.0);
    if (!bad2) {
      rotR2 =  row->Real2[kif];
      rotI2 = -row->Imag2[kif];
      norm = rotR2*rotR2 + rotI2*rotI2;
      if (norm>1.0e-25) norm = 1.0 / sqrt(norm);
      else norm = 1.0;
      rotR2 *= norm;
      rotI2 *= norm;
    }
  }
  
  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
      jif = iif - 1;
      if (istoke==1) {  /* First pol */
	if ((row->Rate1[jif]!=fblank) && (row->Real1[jif]!=fblank) &&  (!bad1)) {
	  tre = row->Real1[jif];
	  tim = row->Imag1[jif];
	  row->Real1[jif]   = tre*rotR1 - tim*rotI1;
	  row->Imag1[jif]   = tre*rotI1 + tim*rotR1;
	  /* end not blanked */
	} else {
	  row->Real1[jif]   = fblank;
	  row->Imag1[jif]   = fblank;
	  row->Weight1[jif] = 0.0;
	}
      } else {       /* second pol */
	if ((row->Rate2[jif]!=fblank) && (row->Real2[jif]!=fblank) && (!bad2)) {
	  tre = row->Real2[jif];
	  tim = row->Imag2[jif];
	  row->Real2[jif]   = tre*rotR2 - tim*rotI2;
	  row->Imag2[jif]   = tre*rotI2 + tim*rotR2;
	  /* end not blanked */
	} else {
	  row->Real2[jif]   = fblank;
	  row->Imag2[jif]   = fblank;
	  row->Weight2[jif] = 0.0;
	}
      }
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doREFP */

/*----------------------------------------------------------------------- */
/*  Apply CPRT corrections to SN table                                    */
/*  Copy the residual fringe-rate from the IF denoted by SNCorParm[0] to  */
/*  all specified IFs.  Obviously If the rate for the reference IF is     */
/*  flagged all IFS will be flagged.  Default = IF 1.                     */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doCPRT (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif, kif;
  ofloat fblank;
  gboolean bad1=FALSE, bad2=FALSE;
  if (err->error) return;

  fblank = control->fblank;
  kif = control->SNCorParm[0]-0.5;
  if (control->bStoke==1) {
    bad1 = (row->Real1[kif]!=fblank) || (row->Imag1[kif]!=fblank) || 
      (row->Weight1[kif]<=0.0);
  }
  
  if ((control->eStoke==2) || (control->bStoke==2)) {
    bad2 = (row->Real2[kif]!=fblank) || (row->Imag2[kif]!=fblank) || 
      (row->Weight2[kif]<=0.0);
  }

  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
      jif = iif - 1;
      if (istoke==1) {  /* First pol */
	if ((row->Rate1[jif]!=fblank) && (row->Real1[jif]!=fblank) && (!bad1)) {
	  row->Rate1[jif]   = row->Rate1[kif];
	  /* end not blanked */
	} else {
	  row->Rate1[jif]   = fblank;
	  row->Weight1[jif] = 0.0;
	}
      } else {       /* second pol */
	if ((row->Rate2[jif]!=fblank) && (row->Real2[jif]!=fblank) && (!bad2)) {
	  row->Rate2[jif]   = row->Rate2[kif];
	  /* end not blanked */
	} else {
	  row->Rate2[jif]   = fblank;
	  row->Weight2[jif] = 0.0;
	}
      }
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doCPRT */

/*----------------------------------------------------------------------- */
/*  Apply CPSN corrections to SN table                                    */
/*  Copy the whole solution from the IF denoted by SNCorParm[0] to all    */
/*  specified IFs.  Obviously If the solution for the reference IF is     */
/*  flagged all IFS will be flagged. Default = IF 1.                      */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doCPSN (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif, kif;
  ofloat fblank;
  if (err->error) return;

  fblank = control->fblank;
  kif = control->SNCorParm[0]-0.5;
  
  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
      jif = iif - 1;
      if (istoke==1) {  /* First pol */
	if ((row->Real1[jif]!=fblank) && (row->Imag1[jif]!=fblank) &&
	    (row->Real1[kif]!=fblank) && (row->Imag1[kif]!=fblank)) {
	  row->Real1[jif]   = row->Real1[kif];
	  row->Imag1[jif]   = row->Imag1[kif];
	  row->Delay1[jif]  = row->Delay1[kif];
	  row->Rate1[jif]   = row->Rate1[kif];
	  row->Weight1[jif] = row->Weight1[kif];
	  row->RefAnt1[jif] = row->RefAnt1[kif];
	  /* end not blanked */
	} else {
	  row->Real1[jif]   = fblank;
	  row->Imag1[jif]   = fblank;
	  row->Weight1[jif] = 0.0;
	}
      } else {       /* second pol */
	if ((row->Real2[jif]!=fblank) && (row->Imag2[jif]!=fblank) &&
	   (row->Real2[kif]!=fblank) && (row->Imag2[kif]!=fblank)) {
	  row->Real2[jif]   = row->Real2[kif];
	  row->Imag2[jif]   = row->Imag2[kif];
	  row->Delay2[jif]  = row->Delay2[kif];
	  row->Rate2[jif]   = row->Rate2[kif];
	  row->Weight2[jif] = row->Weight2[kif];
	  row->RefAnt2[jif] = row->RefAnt2[kif];
	  /* end not blanked */
	} else {
	  row->Real2[jif]   = fblank;
	  row->Imag2[jif]   = fblank;
	  row->Weight2[jif] = 0.0;
	}
      }
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doCPSN */

/*----------------------------------------------------------------------- */
/*  Apply PCOP corrections to SN table                                    */
/*  Copy the whole solution from one polarization to the other, for all   */
/*  specified IF's. The copy direction is specified by SNCorParm[0]:      */
/*      SNCorParm[0] = 1   => Copy polzn. 1 (R) to 2 (L)                  */
/*                   = 2   => Copy polzn. 2 (L) to 1 (R)                  */
/*  No action is taken if there is only one polarization in the SN table. */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doPCOP (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, jif;
  gboolean R2L=FALSE;

  if (err->error) return;

  /* What wanted */
  if (control->SNCorParm[0]==1.0) R2L = TRUE;
  else if (control->SNCorParm[0]==2.0) R2L = FALSE;

  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    if (R2L) { /* Copy Poln 1 to Poln 2 */
      row->MBDelay2     = row->MBDelay1;
      row->Real2[jif]   = row->Real1[jif];
      row->Imag2[jif]   = row->Imag1[jif];
      row->Delay2[jif]  = row->Delay1[jif];
      row->Rate2[jif]   = row->Rate1[jif];
      row->Weight2[jif] = row->Weight1[jif];
      row->RefAnt2[jif] = row->RefAnt1[jif];
    } else {       /* Copy Poln 2 to Poln 1 */
      row->MBDelay1     = row->MBDelay2;
      row->Real1[jif]   = row->Real2[jif];
      row->Imag1[jif]   = row->Imag2[jif];
      row->Delay1[jif]  = row->Delay2[jif];
      row->Rate1[jif]   = row->Rate2[jif];
      row->Weight1[jif] = row->Weight2[jif];
      row->RefAnt1[jif] = row->RefAnt2[jif];
    }
  } /* end loop over IF */
 
} /* end doPCOP */

/*----------------------------------------------------------------------- */
/*  Apply PNEG corrections to SN table                                    */
/*  Flip the sign of the gain phase for all selected SN solutions.        */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doPNEG (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat fblank;
  if (err->error) return;

  fblank = control->fblank;
  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  if (row->Imag1[jif]!=fblank) {
	    row->Imag1[jif] = -row->Imag1[jif];
	  } /* end not blanked */
	} else {       /* second pol */
	  if (row->Imag2[jif]!=fblank) {
	    row->Imag2[jif] = -row->Imag2[jif];
	  } /* end not blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doPNEG */

/*----------------------------------------------------------------------- */
/*  Apply NORM corrections to SN table                                    */
/*  This option normalizes amplitudes to 1.                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV with SN table to modify                          */
/*      SNTable   SN table to modify                                      */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doNORM (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat norm, fblank;
  if (err->error) return;

  fblank = control->fblank;
  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  if ((row->Real1[jif]!=fblank) && (row->Imag1[jif]!=fblank)) {
	    norm = row->Real1[jif]*row->Real1[jif] + 
	           row->Imag1[jif]*row->Imag1[jif];
	    if (norm>1.0e-25) norm = 1.0 / sqrt(norm);
	    else norm = 1.0;
	    row->Real1[jif] *= norm;
	    row->Imag1[jif] *= norm;
	  } /* end not blanked */
	} else {       /* second pol */
	  if ((row->Real2[jif]!=fblank) && (row->Imag2[jif]!=fblank)) {
	    norm = sqrt (row->Real2[jif]*row->Real2[jif] + 
			 row->Imag2[jif]*row->Imag2[jif]);
	    if (norm>1.0e-25) norm = 1.0 / norm;
	    else norm = 1.0;
	    row->Real2[jif] *= norm;
	    row->Imag2[jif] *= norm;
	  } /* end not blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doNORM */

/*----------------------------------------------------------------------- */
/*  Apply RSET corrections to SN table                                    */
/*  Resets gains to (1,0) zero delay, rate and weight = 1.0; even if      */
/*  previously blanked                                                    */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      row       SN table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doRSET (ControlInfo *control, ObitTableSNRow *row, ObitErr *err)
{
  olong iif, istoke, jif;
  if (err->error) return;

  /* Loop over IF */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  row->IFR = 0.0;
	  row->MBDelay1     = 0.0;
	  row->Real1[jif]   = 1.0;
	  row->Imag1[jif]   = 0.0;
	  row->Delay1[jif]  = 0.0;
	  row->Rate1[jif]   = 0.0;
	  row->Weight1[jif] = 1.0;
	  row->RefAnt1[jif] = 1;
	} else {       /* second pol */
	  row->MBDelay2     = 0.0;
	  row->Real2[jif]   = 1.0;
	  row->Imag2[jif]   = 0.0;
	  row->Delay2[jif]  = 0.0;
	  row->Rate2[jif]   = 0.0;
	  row->Weight2[jif] = 1.0;
	  row->RefAnt2[jif] = 1;
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
 
} /* end doRSET */

