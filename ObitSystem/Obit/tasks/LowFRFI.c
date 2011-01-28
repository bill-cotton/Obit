/* $Id$  */
/* Low Frequency Radio Interferometry RFI removal                     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2010                                          */
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
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitUVRFIXize.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* LowFRFIIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void LowFRFIOut (ObitInfoList* outList, ObitErr *err);
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
/* Get residual data */
ObitUV* getResidualData (ObitInfoList *myInput, ObitErr *err);
/* Get output data */
ObitUV* getOutputData (ObitInfoList *myInput, ObitUV* inData, ObitErr *err);
/* Write history */
void LowFRFIHistory (ObitInfoList* myInput, ObitUV* inData, ObitUVRFIXize *myRFI,
		     ObitErr* err);

/* Program globals */
gchar *pgmName = "LowFRFI";       /* Program name */
gchar *infile  = "LowFRFI.in" ;   /* File with program inputs */
gchar *outfile = "LowFRFI.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Obit task to determine gains for a uv data set                       */
/*----------------------------------------------------------------------- */
{
  oint          ierr = 0;
  ObitSystem    *mySystem= NULL;
  ObitUV        *inData=NULL, *resData=NULL, *outData=NULL;
  ObitUVRFIXize *myRFI=NULL;
  ObitErr       *err= NULL;
  gchar        *RFIParms[] = {  /* Control parameters */
    "timeAvg", "minRFI", "timeInt", "solInt", "doInvert", "minRot", "maxRot",
    NULL};

  /* Startup - parse command line */
  err = newObitErr();
  myInput = LowFRFIIn (argc, argv, err);
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

  /* Get input uvdata */
  inData = getInputData (myInput, err);

  /* Get Residual uvdata */
  resData = getResidualData (myInput, err);

  /* Get output uvdata */
  outData = getOutputData (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Create object */
  myRFI = ObitUVRFIXizeCreate ("RFI", inData, resData, outData);
  
  /* Get control parameters from myInput, copy to myRFI */
  ObitInfoListCopyList (myInput, myRFI->info, RFIParms);
  
  /* Do it */
  ObitUVRFIXizeCounterRot (myRFI, err);
  ObitUVRFIXizeFilter (myRFI, err);
  ObitUVRFIXizeCorrect (myRFI, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Write history */
  LowFRFIHistory (myInput, outData, myRFI, err); 

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  resData   = ObitUnref(resData);
  outData   = ObitUnref(outData);
  myRFI = ObitUVRFIXizeUnref(myRFI);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput  = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  err   = ObitErrUnref(err);
 
  return ierr;
} /* end of main */

ObitInfoList* LowFRFIIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "LowFRFIIn";

  /* error checks */
  if (err->error) return list;

  /* Make default inputs InfoList */
  list = defaultInputs(err);

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
  myOutput = defaultOutputs(err);
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  return list;
} /* end LowFRFIIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: LowFRFI -input file -output ofile [args]\n");
    fprintf(stderr, "LowFRFI estimate RFI and remove from UV data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def LowFRFI.in\n");
    fprintf(stderr, "  -output parameter file, def LowFRFI.out\n");
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
/*     Sources   Str (16,1) Sources selected, blank = all                 */
/*     Stokes    Str (4)    Stokes parameter to image, def=I              */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*     doCalSelect Boo (1)  Apply calibration/selection?  def=False       */
/*     doCalib   Int (1)    >0 => apply calibration, 2=> cal. wt, def=-1  */
/*     gainUse   Int (1)    Gain table (CL/SN) table to apply, 0=> highest*/
/*     doBand    Int (1)    If >0.5 apply bandpass cal.                   */
/*     flagVer   Int (1)    Flagging table version, def=0                 */
/*     BPVer     Int (1)    Bandpass table version, 0=highest, def=0      */
/*     doPol     Boo (1)    Apply polarization calibration?, def=False    */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat farray[3];
  gboolean btemp;
  ObitInfoList *out = newObitInfoList();
  gchar *routine = "defaultInputs";

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
  strTemp = "LowFRFI.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "inFile", OBIT_string, dim, strTemp);

  /* input AIPS file name */
  strTemp = "            ";
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

  /* Sources selected, blank = all */
  strTemp = "                ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "Sources", OBIT_string, dim, strTemp);
    
  /* Timerange in days */
  dim[0] = 2;dim[1] = 1;
  farray[0] = -1.0e20; farray[1] = 1.0e20;
  ObitInfoListAlwaysPut (out, "timeRange", OBIT_float, dim, farray);

  /*  Apply calibration/selection?, def=False */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListAlwaysPut (out, "doCalSelect", OBIT_bool, dim, &btemp);

  /*  >0 => apply gain calibration, 2=> cal. wt, def=no cal. */
  dim[0] = 1;dim[1] = 1;
  itemp = -1; 
  ObitInfoListAlwaysPut (out, "doCalib", OBIT_oint, dim, &itemp);

  /*  Gain table (CL/SN) table to apply, 0=> highest, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListAlwaysPut (out, "gainUse", OBIT_oint, dim, &itemp);

  /*  If >0.5 apply bandpass cal, def = no BP cal. */
  dim[0] = 1;dim[1] = 1;
  itemp = -1; 
  ObitInfoListAlwaysPut (out, "doBand", OBIT_oint, dim, &itemp);

  /*  Bandpass table version, 0=highest, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListAlwaysPut (out, "BPVer", OBIT_oint, dim, &itemp);

  /* Flagging table version, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListAlwaysPut (out, "flagVer", OBIT_oint, dim, &itemp);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* Apply polarization calibration?, def= False */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListAlwaysPut (out, "doPol", OBIT_bool, dim, &btemp);
  
  /* Subarray */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListAlwaysPut (out, "subA", OBIT_oint, dim, &itemp);

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
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *strTemp;
  gchar        outFile[129];
  gchar        Aname[13], Aclass[7];
  gboolean     doCalSelect, btemp;
  oint         doCalib;
  olong        i, n, disk, seq;
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

  /* Make sure doCalSelect set properly */
  doCalSelect = TRUE;
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doCalib",  &type, dim, &doCalib);
  doCalSelect = doCalSelect || (doCalib>0);
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);

  /* Output name -
   FITS */
  ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&strTemp);
  /* if not use inName */
  if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12))) {
    ObitInfoListGetP (myInput, "inFile", &type, dim, (gpointer)&strTemp);
    if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12))) {
      n = MIN (128, dim[0]);
      for (i=0; i<n; i++) outFile[i] = strTemp[i]; outFile[i] = 0;
      ObitTrimTrail(outFile);  /* remove trailing blanks */
      dim[0] = strlen(outFile);
      ObitInfoListAlwaysPut (myInput, "outFile", OBIT_bool, dim, outFile);
    }
  }

  /* AIPS */
  /* outName given? */
  ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&strTemp);
  /* if not, use inName */
  if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12))) {
    ObitInfoListGetP (myInput, "inName", &type, dim, (gpointer)&strTemp);
    for (i=0; i<12; i++) Aname[i] = ' ';  Aname[i] = 0;
    for (i=0; i<MIN(12,dim[0]); i++) Aname[i] = strTemp[i];
    /* Save any defaulting on myInput */
    dim[0] = 12;
    ObitInfoListAlwaysPut (myInput, "outName", OBIT_string, dim, Aname);
  }

  /* output AIPS class,  Default is "LowFRF"*/
  ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp);
  if ((strTemp==NULL) || (!strncmp(strTemp, "      ", 6))) {
    strncpy (Aclass, "LowFRF", 7);
    dim[0] = 6;
    ObitInfoListAlwaysPut (myInput, "outClass", OBIT_string, dim, Aclass);
  }

  /* output AIPS disk - default is inDisk */
  disk = 0;
  ObitInfoListGetTest(myInput, "outDisk", &type, dim, &disk);
  if (disk<=0) {
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    dim[0] = 1;
    ObitInfoListAlwaysPut (myInput, "outDisk", OBIT_oint, dim, &disk);
  }

  /* output AIPS sequence - default is new (0) */
  seq = 0;
  if (!ObitInfoListGetTest(myInput, "outSeq", &type, dim, &seq)) {
    dim[0] = 1;
    ObitInfoListAlwaysPut (myInput, "outSeq", OBIT_oint, dim, &seq);
  }

  /*  Output file not expected to exist */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListAlwaysPut (myInput, "outExist", OBIT_bool, dim, &btemp);

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
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "Stokes", "timeRange", "BChan", "EChan",   "BIF", "EIF", "subA",
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer", "doPol",
    "Antennas",
     NULL};
  gchar *routine = "getInputData";

  /* error checks */
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Build basic input UV data Object */
  inData = ObitUVFromFileInfo ("in", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Ensure inData fully instantiated and OK */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Get residual data                                                     */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitUV with residual data                                        */
/*----------------------------------------------------------------------- */
ObitUV* getResidualData (ObitInfoList *myInput, ObitErr *err)
{
  ObitUV       *resData = NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        flagVer;
  gchar *routine = "getResidualData";

  /* error checks */
  if (err->error) return resData;
  g_assert (ObitInfoListIsA(myInput));

  /* Build basic input UV data Object */
  resData = ObitUVFromFileInfo ("in2", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", resData);
  
  /* copy flag version */
  flagVer = -1;
  ObitInfoListGetTest (myInput, "flagRes", &type, dim, &flagVer); 
  ObitInfoListAlwaysPut (resData->info, "flagVer", OBIT_oint, dim, &flagVer);

  return resData;
} /* end getResidualData */

/*----------------------------------------------------------------------- */
/*  Get output data                                                       */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV from which to clone output                 */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitUV with residual data                                        */
/*----------------------------------------------------------------------- */
ObitUV* getOutputData (ObitInfoList *myInput, ObitUV* inData, ObitErr *err)
{
  ObitUV  *outData = NULL;
  gchar *routine = "getOutputData";

  /* error checks */
  if (err->error) return outData;
  g_assert (ObitInfoListIsA(myInput));

  /* Build basic input UV data Object */
  outData = ObitUVFromFileInfo ("out", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);

  /* Clone from input */
  ObitUVClone (inData, outData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);
  
  /* Ensure outData fully instantiated and OK */
  ObitUVFullInstantiate (outData, FALSE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);

  return outData;
} /* end getOutputData */

/*----------------------------------------------------------------------- */
/*  Write History for LowFRFI                                             */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*      myRFI     RFI Object                                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void LowFRFIHistory (ObitInfoList* myInput, ObitUV* inData, 
		     ObitUVRFIXize *myRFI, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "in2File",  "in2Disk", "in2Name", "in2Class", "in2Seq", 
    "Sources", "Qual", "souCode", "timeRange",  "subA",
    "selBand", "selFreq", "FreqID", "BChan", "EChan", 
    "doCalSelect",  "doCalib",  "gainUse",  "doBand ",  "BPVer",  "flagVer", 
    "doPol", "Antennas", 
    "timeAvg", "minRFI", "timeInt", "solInt", "doInvert", "minRot", "maxRot",
    NULL};
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat fractFlag=0.0, fractMod=0.0;
  gchar *routine = "LowFRFIHistory";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

   /* Do history */
  outHistory = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadWrite, err);

  /* Add this programs history 
     Copy selected values from myInput */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Results */
  ObitInfoListGetTest(myRFI->info, "fractFlag",  &type, dim, &fractFlag);
  ObitInfoListGetTest(myRFI->info, "fractMod",   &type, dim, &fractMod);
  g_snprintf (hicard, 80, "%s / fraction flagged %8.5f",pgmName, fractFlag);
  ObitHistoryWriteRec (outHistory, -1, hicard, err);
  g_snprintf (hicard, 80, "%s / fraction modified %8.5f",pgmName, fractMod);
  ObitHistoryWriteRec (outHistory, -1, hicard, err);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  outHistory = ObitHistoryUnref(outHistory); /* cleanup */
 
} /* end LowFRFIHistory  */

