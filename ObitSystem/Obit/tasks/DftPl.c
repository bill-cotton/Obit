/* $Id$  */
/* Obit Task to Plot average uv data v time          .                */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2018                                               */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitUV.h"
#include "ObitPlot.h"
#include "ObitAIPSDir.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* DftPlIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void DftPlOut (ObitInfoList* outList, ObitErr *err);
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
/* Average data */
void AvgData (ObitUV* inData, olong *nplot, ofloat* plotS, ofloat* plotSerr, 
	      ofloat* plotTime, ObitErr *err);
/* Plot data */
void PlotData (ObitUV* inData, olong nplot, ofloat* plotS, ofloat* plotSerr, 
	       ofloat* plotTime, ObitErr *err);

/* Program globals */
gchar *pgmName = "DftPl";       /* Program name */
gchar *infile  = "DftPl.in" ;   /* File with program inputs */
gchar *outfile = "DftPl.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL;  /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Obit test uv data rutines                                            */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem=NULL;
  ObitUV       *inData=NULL;
  ObitErr      *err= NULL;
#define MAXPLOT  10000 /* Maximum number of samples to plot */
  ofloat plotS[MAXPLOT], plotSerr[MAXPLOT], plotTime[MAXPLOT]; 
  olong        nplot=MAXPLOT;

  /* Startup - parse command line, read inputs */
  err = newObitErr();

  myInput = DftPlIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Digest inputs */
  digestInputs(myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get input uvdata */
  inData = getInputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Average data */
  AvgData(inData, &nplot, plotS, plotSerr, plotTime, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Plot data */
  PlotData(inData, nplot, plotS, plotSerr, plotTime, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
 
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* DftPlIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "DftPlIn";

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
      
     } else if (strcmp(arg, "-DataType") == 0) { /* Data type AIPS or FITS */
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

    } else if (strcmp(arg, "-BChan") == 0) { /* BChan */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "BChan", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-EChan") == 0) { /* EChan */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "EChan", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-BIF") == 0) { /* BIF */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "BIF", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-EIF") == 0) { /* EIF */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "EIF", OBIT_oint, dim, &itemp, err);
      
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
} /* end DftPlIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: DftPl -input file -output ofile [args]\n");
    fprintf(stderr, "DftPl - plots average data v time\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def DftPl.in\n");
    fprintf(stderr, "  -output output result file, def DftPl.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input image\n");
    fprintf(stderr, "  -inFile input FITS UV file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input image (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -BChan first channel to copy\n");
    fprintf(stderr, "  -EChan highest channel to copy\n");
    fprintf(stderr, "  -BIF first IF to copy\n");
    fprintf(stderr, "  -EIF highest IF to copy\n");
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
/*     inFile    Str [?]    input FITS uv file name [no def]              */
/*     inDisk    Int        input AIPS or FITS uv disk no  [def 1]        */
/*     inName    Str [12]   input AIPS uv name  [no def]                  */
/*     inClass   Str [6]    input AIPS uv class  [no def]                 */
/*     inSeq     Int        input AIPS uv sequence no  [no def]           */
/*     BChan     Int [1]    channel number, 0=>all                        */
/*     EChan     Int [1]    channel number, 0=>all                        */
/*     BIF       Int [1]    first IF to process, 0=>1                     */
/*     EIF       Int [1]    highest IF to process, 0=>all                 */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
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
  ObitInfoListPut (out, "outDType", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input FITS file name */
  strTemp = "DftPl.intab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS input uv file name */
  strTemp = "DftPlName";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS input uv file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* AIPS input uv sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "inSeq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* AIPS or FITS input uv disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "inDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* channels  */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "BChan", OBIT_oint, dim, &itemp, err);
  itemp = 0; 
  ObitInfoListPut (out, "EChan", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* IFs  */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "BIF", OBIT_oint, dim, &itemp, err);
  itemp = 0; 
  ObitInfoListPut (out, "EIF", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Only want cross correlations */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "corrType", OBIT_string, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

   /* Stokes always "   " */
  strTemp = "    ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Stokes", OBIT_string, dim, strTemp, err);
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
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar Stokes[4]; 
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  ObitInfoListGetTest(myInput, "Stokes",  &type, dim, Stokes);
  if (Stokes[0]==' ') Stokes[0] = 'I';
  ObitInfoListAlwaysPut (myInput, "Stokes", type, dim, Stokes);

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

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
  ObitInfoType type;
  olong         Aseq, disk, cno, nvis=1000;
  gchar        *Type, *strTemp, inFile[129];
  oint         doCalib;
  gchar        Aname[13], Aclass[7], *Atype = "UV";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doCalSelect;
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "souCode", "Qual", "Stokes", "timeRange", "UVRange",
    "BChan", "EChan", "chanInc", "BIF", "EIF", "IFInc", "FreqID", "corrType", 
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer", 
    "doPol", "PDVer", "Smooth", "Antennas",  "subA", "Sources", "souCode", "Qual",
    "timeAvg", "Shift", "title", "Range", "nplot",
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
    inFile[128] = 0;
    ObitTrimTrail(inFile);  /* remove trailing blanks */
    
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

  /* Make sure doCalSelect set properly */
  doCalSelect = TRUE;
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doCalib",  &type, dim, &doCalib);
  doCalSelect = doCalSelect || (doCalib>0);
  /* Always */
  doCalSelect = TRUE;
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);
 

  /* Ensure inData fully instantiated and OK */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/**
 * Time Average data for plotting
 * \param inData   Input uv data to average, 
 *                 Any request for calibration, editing and selection honored
 * \param nplot    [in] maximum number of points to plot
 *                 [out] Number of points actually plotted.
 * \param plotS    Average flux density per time bin
 * \param plotSerr RMS error of points in plotS
 * \param plotTime Time(hours) of points in plotS
 * \param err      Error stack, returns if not empty.
 */
void AvgData (ObitUV* inData, olong *nplot, ofloat* plotS, ofloat* plotSerr, 
	      ofloat* plotTime, ObitErr *err)
{
  ObitIOCode iretCode;
  gboolean doCalSelect;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong ncorr, nrparm;
  ollong i, indx, iindx=0;
  ObitIOAccess access;
  ObitUVDesc *inDesc;
  olong lastSourceID, curSourceID;
  olong ilocu, ilocv;
  ofloat u, v, curTime, startTime, endTime;
  ofloat *inBuffer;
  odouble irefFreq;
  olong ivis=0, cnt, iplot,maxPlot=*nplot;
  gboolean done, gotOne, doShift;
  ofloat timeAvg, shift[2], sumWt, sumTime, sumRe, sumRe2, phase;
  gchar *routine = "AvgData";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inData));

  /* Selection/calibration/editing of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inData->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* Other control */
  timeAvg = 0.25;
  ObitInfoListGetTest(inData->info, "timeAvg", &type, dim, &timeAvg);
  if (timeAvg<=0.0) timeAvg = 0.25;
  timeAvg /= 60.*24.;  /* to days */
  shift[0] = shift[1] = 0.0;
  ObitInfoListGetTest(inData->info, "Shift", &type, dim, shift);
  /* to radians * 2 pi */
  shift[0] *= AS2RAD * 2.0 * G_PI;
  shift[1] *= AS2RAD * 2.0 * G_PI;
  doShift = ((shift[0]!=0.0) || (shift[1]!=0.0));
  
  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inData, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inData->name);
  /* Get descriptors */
  inDesc  = inData->myDesc;
  ilocu = inDesc->ilocu; ilocv = inDesc->ilocv;
  irefFreq = 1.0 / inDesc->freq;
  /* data info */
  ncorr   = inData->myDesc->ncorr;
  nrparm  = inData->myDesc->nrparm;

  /* Initialize things */
  startTime = -1.0e20;
  endTime   =  1.0e20;
  lastSourceID = -1;
  curSourceID  = 0;
  inBuffer  = inData->buffer;
  done   = FALSE;
  gotOne = FALSE;
  cnt = 0; sumWt = sumTime = sumRe = sumRe2 = 0.0;

  /* we're in business, average data */
  iplot = 0;
  while (iretCode==OBIT_IO_OK) {
    if ((!gotOne) || (inData->myDesc->numVisBuff<=0)) { /* need to read new record? */
      if (doCalSelect) iretCode = ObitUVReadSelect (inData, inData->buffer, err);
      else iretCode = ObitUVRead (inData, inData->buffer, err);
    }

    /* Are we there yet??? */
    done = (inDesc->firstVis >= inDesc->nvis) || (iretCode==OBIT_IO_EOF);
    if (done && (startTime>0.0)) goto process; /* Final? */

    /* Make sure valid data found */
    if (inData->myDesc->numVisBuff<=0) continue;

    /* loop over visibilities */
    for (ivis=0; ivis<inDesc->numVisBuff; ivis++) { 
      iindx = ivis*inDesc->lrec;
      gotOne = FALSE;
      
      curTime = inBuffer[iindx+inDesc->iloct]; /* Time */
      if (inDesc->ilocsu>=0) curSourceID = inBuffer[iindx+inDesc->ilocsu];
      if (startTime < -1000.0) {  /* Set time window etc. if needed */
	startTime = curTime;
	endTime   = startTime + timeAvg;
	lastSourceID = curSourceID;
      }

      /* Still in current interval/source? */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (inDesc->firstVis<=inDesc->nvis) && (iretCode==OBIT_IO_OK)) {
	
	/* Accumulate*/
	indx = iindx+nrparm; /* offset of start of vis data */
	for (i=0; i<ncorr; i++) {
	  if (inBuffer[indx+2] > 0.0) {
	    cnt ++;
	    sumWt   += inBuffer[indx+2];
	    sumTime += curTime;
	    /* Shift? */
	    if (doShift) {
	      u = (ofloat)inBuffer[ilocu]*inDesc->freqArr[i]*irefFreq;  /* Scale to freq */
	      v = (ofloat)inBuffer[ilocv]*inDesc->freqArr[i]*irefFreq;
	      phase = shift[0]*u + shift[1]*v;
	      inBuffer[indx] = inBuffer[indx]*cos(phase) - inBuffer[indx+1]*sin(phase);
	    }
	    sumRe   += inBuffer[indx]*inBuffer[indx+2];
	    sumRe2  += inBuffer[indx]*inBuffer[indx]*inBuffer[indx+2];
	  } 
	  indx += 3;
	} /* end loop over correlations */;
      } else {  /* process interval */
	
      process:
	/* Time bin statistics */
	/* Check that there is room */
	if (iplot>=maxPlot) {
	  Obit_log_error(err, OBIT_InfoWarn,"%s plot length truncated at %d",
			 routine, maxPlot);
	  goto done;
	}
	if ((cnt>1) && (sumWt>0.0)) {
	  plotTime[iplot] = 24.0 * sumTime / cnt;  /* hours */
	  plotS[iplot]    = sumRe/sumWt;
	  plotSerr[iplot] = sqrtf(MAX(0.0, sumRe2/sumWt -plotS[iplot]*plotS[iplot])) ;
	  plotSerr[iplot] /= sqrtf((ofloat)(cnt-1));  /* of mean */
	  iplot++;
	} /* End any data this baseline */
	/* Reinitialize things */
	cnt = 0; sumWt = sumTime = sumRe = sumRe2 = 0.0;
	startTime = -1.0e20;
	endTime   =  1.0e20;
      } /* end process interval */
    } /* end loop over buffer */
    if (done) goto done;
  } /* End loop over input file */
  
  /* End of processing */
 done:
  /* In case anything left
     Check that there is room */
  if ((cnt>0) && (iplot>=maxPlot)) {
    Obit_log_error(err, OBIT_InfoWarn,"%s plot length truncated at %d",
		   routine, maxPlot);
    cnt = 0;
  }
  if ((cnt>1) && (sumWt>0.0)) {
    plotTime[iplot] = 24.0 * sumTime / cnt;  /* hours */
    plotS[iplot]    = sumRe/sumWt;
    plotSerr[iplot] = sqrtf(MAX(0.0, sumRe2/sumWt -plotS[iplot]*plotS[iplot]));
    plotSerr[iplot] /= sqrtf((ofloat)(cnt-1));  /* of mean */
    iplot++;
  } /* End any data left */
  *nplot = iplot;   /* size of plot */
  iretCode = ObitUVClose (inData, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, inData->name);

} /* end AvgData */

/**
 * Plot time series of measurements
 * \param inData   Input uv data to plot - used for labeling 
 * \param nplot    Number of points to plot
 * \param plotS    Average flux density per time bin
 * \param plotSerr RMS error of points in plotS
 * \param plotTime Time(hours) of points in plotS
 * \param err      Error stack, returns if not empty.
 */
void PlotData (ObitUV* inData, olong nplot, ofloat* plotS, ofloat* plotSerr, 
	       ofloat* plotTime, ObitErr *err)
{
  ObitPlot *plot=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *defPlotFile="Plot.png", *defFormat="png", *defTitle="   ";
  gchar *outputfile=NULL, *plotFile=NULL, *format=NULL, label[500], *title=NULL;
  ofloat dmin, dmax, sum1, sum2, sig, ymin, ymax, Range[2];
  olong i, nplotpp=1;
  gchar *routine = "PlotData";

  if (err->error) return; /* previous error? */
  
  /* Better be some */
  Obit_return_if_fail ((nplot>=3), err, "%s: Only %d point(s) in plot", routine, nplot); 

  if (!ObitInfoListGetP(myInput, "format", &type, dim, (gpointer)&format)) {
    format = defFormat;
  }
  ObitTrimTrail (format);
  if (!ObitInfoListGetP(myInput, "plotFile", &type, dim, (gpointer)&plotFile)) {
    plotFile = defPlotFile;
  }
  ObitTrimTrail (plotFile);
  if (!ObitInfoListGetP(myInput, "title", &type, dim, (gpointer)&title)) {
    title = defTitle;
  } 
  if (title[0]!=' ') ObitTrimTrail (title);

  Range[0] = Range[1] = 0.0;
  ObitInfoListGetTest(myInput, "Range", &type, dim, Range);
  if ((Range[0]!=0.0) || (Range[1]!=0.0)) {
    ymin = Range[0]; ymax = Range[1];
  } else {
    /* PLot statistics for range */
    sum1 = sum2 = 0.0; dmin=1.0e10; dmax = -1.0e10;
    for (i=0; i<nplot; i++) {
      sum1 += plotS[i];    sum2 += plotSerr[i];
      sig = sum2 / nplot;
      dmin = MIN(plotS[i], dmin);
      dmax = MAX(plotS[i], dmax);
    }
    ymin = dmin-3*sig;
    ymax = dmax+3*sig;
  }

  ObitInfoListGetTest(myInput, "nplot", &type, dim, &nplotpp);

  /* Output file name */
  outputfile = g_strconcat (plotFile, "/", format, NULL);

  /* Make plot White background */
  plot = newObitPlot ("Plot");
  ObitPlotInitPlot (plot, outputfile, 15, 1, nplotpp, err);
  /* Labeling */
  dim[0] = strlen(title); dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(plot->info,"TITLE", OBIT_string, dim, title);
  strcpy (label, "Time (hours)");
  dim[0] = strlen(label); dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(plot->info,"XLABEL", OBIT_string, dim, label);
  strcpy (label, "Flux density (Jy)");
  dim[0] = strlen(label); dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(plot->info,"YLABEL", OBIT_string, dim, label);
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(plot->info,"YMIN", OBIT_float, dim, &ymin);
  ObitInfoListAlwaysPut(plot->info,"YMAX", OBIT_float, dim, &ymax);
  /* Plot */
  ObitPlotXYErr (plot, 2, nplot, plotTime, plotS, plotSerr, err);
  /* finalize plot */
  ObitPlotFinishPlot (plot, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  /* Cleanup */
  if (outputfile) g_free(outputfile);
  ObitPlotUnref(plot);
} /* end PlotData */
