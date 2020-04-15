/* $Id$  */
/* Obit Radio interferometry calibration software                     */
/* Copy calibration solution data                                     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2020                                               */
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
  /* Input and output number of IFs */
  olong inNumIF, outNumIF;
  /* Input and output number of Stokes */
  olong inNumStok, outNumStok;
  /* Array of input IF's (0-rel) closest to output IFs */
  olong *outCloseIn;
  /* Array of frequency differences between input and output IFs */
  ofloat *FreqDif;
  /**  blanking value */
  ofloat fblank;
 /**  SNCpyParm */
  ofloat SNCpyParm[64];
}  ControlInfo; 
  
/* internal prototypes */
/* Get inputs */
ObitInfoList* SNCpyIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void SNCpyOut (ObitInfoList* outList, ObitErr *err);
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
/* Create output uvdata */
ObitUV* setOutputUV (ObitInfoList *myInput, ObitUV* inData, ObitErr *err);
/* Write history */
void SNCpyHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Modify SN table */
void SNCpyDoCpy (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		 ObitErr* err);
/* Modification functions */
void doClosControl (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		    ControlInfo *control, ObitErr* err);
void doClos (ControlInfo *control, ObitTableSNRow *inrow, 
	     ObitTableSNRow *outrow, ObitErr *err);


/* Program globals */
gchar *pgmName = "SNCpy";       /* Program name */
gchar *infile  = "SNCpy.in" ;   /* File with program inputs */
gchar *outfile = "SNCpy.out";   /* File to contain program outputs */
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
/*  Copy selected data from one SN table to another                         */
/*------------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitUV       *inData=NULL, *outData=NULL;
  ObitSystem   *mySystem= NULL;
  ObitErr      *err= NULL;

  /* Startup - parse command line */
  err = newObitErr();
  myInput = SNCpyIn (argc, argv, err);
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

  /* Get output uvdata  */
  outData = setOutputUV (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Modify */
  SNCpyDoCpy (myInput, inData, outData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Write history */
  SNCpyHistory (myInput, outData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUVUnref(inData);
  outData   = ObitUVUnref(outData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput  = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* SNCpyIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "SNCpyIn";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return list;

  /* Make default inputs InfoList */
  list = defaultInputs(err);
  myOutput = defaultOutputs(err);

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

  /* Initialize output */
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  return list;
} /* end SNCpyIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: SNCpy -input file -output ofile [args]\n");
    fprintf(stderr, "Obit task to SN row data \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def SNCpy.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def SNCpy.out\n");
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
  strTemp = "SNCpy.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "SNCpyName";
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
  ObitInfoType type;
  gchar        str[20];
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

  /* cpyMode default = 'Clos' */
  strcpy(str, "    "); dim[0] = 4; type = OBIT_string;
  ObitInfoListGetTest(myInput, "cpyMode", &type, dim, str);
  if (!strncmp(str, "    ",4)) {
    strcpy(str, "Clos"); dim[0] = 4; type = OBIT_string;
    ObitInfoListAlwaysPut (myInput, "cpyMode", type, dim, str);
  }

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
/*  Create output uv data, defaults to input                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV from which to clone output                 */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/* Returns the output uv data                                             */
/*----------------------------------------------------------------------- */
ObitUV* setOutputUV (ObitInfoList *myInput, ObitUV* inData, ObitErr *err)
{
  ObitUV    *outUV = NULL;
  ObitInfoType type;
  olong      i, n, cno;
  oint      disk, Aseq;
  gchar     *Type, *strTemp=NULL;
  gchar     Aname[13], Aclass[7], *Atype = "UV";
  olong      nvis;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar     tname[129], outFile[129];
  gchar     *routine = "setOutputUV";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Create basic output UV Object */
  g_snprintf (tname, 100, "output UV data");
  outUV = newObitUV(tname);
    
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS output */

    /* outName given? */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&strTemp);
    /* if not use inName */
    if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12)))
      ObitInfoListGetP (myInput, "inName", &type, dim, (gpointer)&strTemp);
    for (i=0; i<12; i++) Aname[i] = ' ';  Aname[i] = 0;
    for (i=0; i<MIN(12,dim[0]); i++) Aname[i] = strTemp[i];
    /* Save any defaulting on myInput */
    dim[0] = 12;
    ObitInfoListAlwaysPut (myInput, "outName", OBIT_string, dim, Aname);

    /* outClass given? */
    ObitInfoListGetP (myInput, "outClass", &type, dim, (gpointer)&strTemp);
    /* if not use inName */
    if ((strTemp==NULL) || (!strncmp(strTemp, "      ", 6)))
      ObitInfoListGetP (myInput, "inClass", &type, dim, (gpointer)&strTemp);
    for (i=0; i<6; i++) Aclass[i] = ' ';  Aclass[i] = 0;
    for (i=0; i<MIN(6,dim[0]); i++) Aclass[i] = strTemp[i];
    /* Save any defaulting on myInput */
    dim[0] = 6;
    ObitInfoListAlwaysPut (myInput, "outClass", OBIT_string, dim, Aclass);

    /* outSeq given? */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);
    if (Aseq<=0) 
      ObitInfoListGet(myInput, "inSeq", &type, dim, &Aseq, err);
    /* Save any defaulting on myInput */
    dim[0] = 1;
    ObitInfoListAlwaysPut (myInput, "outSeq", OBIT_oint, dim, &Aseq);

    /* outDisk given? */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (disk<=0) 
      ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    /* Save any defaulting on myInput */
    dim[0] = 1;
    ObitInfoListAlwaysPut (myInput, "outDisk", OBIT_oint, dim, &disk);

    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
    /* define object */
    nvis = 1;
    ObitUVSetAIPS (outUV, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */

    /* outFile given? */
    ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&strTemp);
    /* if not use inName */
    if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12)))
      ObitInfoListGetP (myInput, "inFile", &type, dim, (gpointer)&strTemp);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) outFile[i] = strTemp[i]; outFile[i] = 0;
    ObitTrimTrail(outFile);  /* remove trailing blanks */

    /* Save any defaulting on myInput */
    dim[0] = strlen(outFile);
    ObitInfoListAlwaysPut (myInput, "outFile", OBIT_string, dim, outFile);

    /* outDisk given? */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (disk<=0) 
      ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    /* Save any defaulting on myInput */
    dim[0] = 1;
    ObitInfoListAlwaysPut (myInput, "outDisk", OBIT_oint, dim, &disk);

    /* define object */
    nvis = 1;
    ObitUVSetFITS (outUV, nvis, disk, outFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outUV;
  }
  
  /* Ensure outData fully instantiated and OK */
  ObitUVFullInstantiate (outUV, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

 return outUV;
} /* end setOutputUV */

/*----------------------------------------------------------------------- */
/*  Write History for SNCpy                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void SNCpyHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "souCode", "Qual", "EditStokes", "BIF", "EIF", 
    "FreqID", "timeRange",  "subA", "Antennas", 
    "inVer", "outVer", "cpyMode", "SNCParm", 
    NULL};
  gchar *routine = "SNCpyHistory";

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
 
} /* end SNCpyHistory  */

/*----------------------------------------------------------------------- */
/*  Copy SN data                                                          */
/*   Input:                                                               */
/*      myInput    Input parameters on InfoList                           */
/*      inData     ObitUV with input SN table                             */
/*      outData    ObitUV with output SN table                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void SNCpyDoCpy (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		 ObitErr* err)
{
  gchar cpyMode[5];
  oint numpol, numif;
  ObitInfoType type;
  olong inSNver, inSNRow, outSNver, outSNRow;
  ofloat fblank = ObitMagicF();
  olong i, itemp, nif, nstok;
  gint32 dim[MAXINFOELEMDIM];
  ObitTableSN *inSNTable=NULL, *outSNTable=NULL;
  ObitTableSNRow *inRow=NULL, *outRow=NULL;
  ObitIOCode retCode;
  ControlInfo control;
  ObitUVSel *sel;
  gboolean wanted;
  gchar *routine = "SNCpyDoCpy";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Get Operation type */
  ObitInfoListGet(myInput, "cpyMode", &type, dim, cpyMode, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Get SN tables */
  itemp = 0;
  ObitInfoListGetTest(myInput, "inVer", &type, dim, &itemp);
  inSNver = itemp;
  numpol = inData->myDesc->inaxes[inData->myDesc->jlocs];
  numif  = inData->myDesc->inaxes[inData->myDesc->jlocif];
  inSNTable = newObitTableSNValue("in SN table", (ObitData*)inData, &inSNver, 
				OBIT_IO_ReadOnly, numpol, numif, err);
  itemp = 0;
  ObitInfoListGetTest(myInput, "outVer", &type, dim, &itemp);
  outSNver = itemp;
  numpol = outData->myDesc->inaxes[outData->myDesc->jlocs];
  numif  = outData->myDesc->inaxes[outData->myDesc->jlocif];
  outSNTable = newObitTableSNValue("out SN table", (ObitData*)outData, &outSNver, 
				   OBIT_IO_ReadWrite, numpol, numif, err);
  /* Clear any extant rows */
  ObitTableClearRows ((ObitTable*)outSNTable, err);
  /* Open */
  ObitTableSNOpen (inSNTable,  OBIT_IO_ReadOnly,  err);
  ObitTableSNOpen (outSNTable, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;

  /* Copy other info */
  outSNTable->numAnt    = inSNTable->numAnt;
  outSNTable->numNodes  = inSNTable->numNodes;
  outSNTable->mGMod     = inSNTable->mGMod;
  outSNTable->isApplied = inSNTable->isApplied;

  /* Set rows */
  inRow   = newObitTableSNRow (inSNTable);
  outRow  = newObitTableSNRow (outSNTable);

  /* Set control information */
  control.fblank = fblank;
  control.outCloseIn = NULL;
  control.FreqDif    = NULL;
  control.inNumIF    = inData->myDesc->inaxes[inData->myDesc->jlocif];
  control.inNumStok  = inData->myDesc->inaxes[inData->myDesc->jlocs];
  control.outNumIF   = outData->myDesc->inaxes[inData->myDesc->jlocif];
  control.outNumStok = outData->myDesc->inaxes[inData->myDesc->jlocs];
  if (inData->myDesc->jlocif>=0) nif = inData->myDesc->inaxes[inData->myDesc->jlocif];
  else nif = 1;
  if (inData->myDesc->jlocs>=0) nstok = inData->myDesc->inaxes[inData->myDesc->jlocs];
  else nstok = 1;
  nstok = MIN (nstok, 2);
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
  for (i=0; i<30; i++) control.SNCpyParm[i] = 0.0;
  ObitInfoListGetTest(myInput,  "SNCParm", &type, dim, control.SNCpyParm);

  /* Open and close input UV data to update selector */
  ObitUVOpen(inData, OBIT_IO_ReadCal, err);
  ObitUVClose(inData, err);
  if (err->error) goto cleanup;
  sel  = inData->mySel;

  /* Any cpyMode specific checks or processing */
  if (!strncmp (cpyMode, "Clos", 4)) {
    /* Find closest IF, freq offsets */
    doClosControl (myInput, inData, outData, &control, err);
  } /* End initial setup */
  if (err->error) goto cleanup;

  /* Loop over table */
  for (inSNRow=1; inSNRow<=inSNTable->myDesc->nrow; inSNRow++) {
    retCode = ObitTableSNReadRow (inSNTable, inSNRow, inRow, err);
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
    if (!strncmp (cpyMode, "Clos", 4)) {
      doClos (&control, inRow, outRow, err);
    } else {  /* Unknown */
      Obit_log_error(err, OBIT_Error,
		     "%s: Unknown cpyMode: %s", routine, cpyMode);
      return;
    } /* End of cal by mode */
    if (err->error) goto cleanup;

    /* Write row */
    outSNRow = -1;
    retCode = ObitTableSNWriteRow (outSNTable, outSNRow, outRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  } /* end loop over table */

  /* Close tables */
  ObitTableSNClose (inSNTable, err);
  ObitTableSNClose (outSNTable, err);

  /* Cleanup */
 cleanup:
  inSNTable = ObitTableSNUnref(inSNTable);
  inRow = ObitTableSNRowUnref(inRow);
  outSNTable = ObitTableSNUnref(outSNTable);
  outRow = ObitTableSNRowUnref(outRow);
  if (err->error) Obit_traceback_msg (err, routine, outSNTable->name);

} /* end SNCpyDoCpy  */

/* Modification functions */
/*----------------------------------------------------------------------- */
/*  Setup for doClos                                                      */
/*   Input:                                                               */
/*      myInput    Input parameters on InfoList                           */
/*      inData     ObitUV with input SN table                             */
/*      outData    ObitUV with output SN table                            */
/*   In/Output:                                                           */
/*      control    control structure                                      */
/*       outCloseIn input index (0-rel) for each output IF                */
/*       FreqDif    Frequency difference (Hz) between in and out (out-in) */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doClosControl (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		    ControlInfo *control, ObitErr* err)
{
  olong iif, oif, bestIF;
  odouble *inFreqIF  = inData->myDesc->freqIF;
  odouble *outFreqIF = outData->myDesc->freqIF;
  odouble delt, best;

  /* Create work arrays */
  control->outCloseIn = g_malloc0(control->outNumIF*sizeof(olong));
  control->FreqDif    = g_malloc0(control->outNumIF*sizeof(ofloat));

  /* Loop over output array finding closest input freq */
  for (oif=0; oif<control->outNumIF; oif++) {
    bestIF = -1; best = 1.0e20;
    for (iif=0; iif<control->inNumIF; iif++) {
      delt = fabs(outFreqIF[oif] - inFreqIF[iif]);
      if (delt<best) {
	best   = delt;
	bestIF = iif;
      }
    } /* end input loop */
    control->outCloseIn[oif] = bestIF;
    control->FreqDif[oif]    = best;
  } /* endi output loop */
} /* end doClosControL */

/*----------------------------------------------------------------------- */
/*  Copy closest input to output with phase update                        */
/*  Average selected fringe rates in a record                             */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      inRow     input SN table row                                      */
/*      outRow    output SN table row                                     */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doClos (ControlInfo *control, ObitTableSNRow *inRow, 
	     ObitTableSNRow *outRow, ObitErr *err) {
  olong iif, oif;
  ofloat dfq, dphs, xr, xi, fblank;
  if (err->error) return;

  fblank = control->fblank;

  /* copy scalars */
  outRow->Time     = inRow->Time;
  outRow->TimeI    = inRow->TimeI;
  outRow->SourID   = inRow->SourID;
  outRow->antNo    = inRow->antNo;
  outRow->SubA     = inRow->SubA;
  outRow->FreqID   = inRow->FreqID;
  outRow->IFR      = inRow->IFR;
  outRow->NodeNo   = inRow->NodeNo;
  outRow->MBDelay1 = inRow->MBDelay1;
  if ((control->inNumStok>1) && (control->outNumStok>1)) {
    outRow->MBDelay2 = inRow->MBDelay2;
  }
  /* Ignore dispersive delay for now */

  /* Loop over output IFs */
  for (oif=0; oif<control->outNumIF; oif++) {
    iif = control->outCloseIn[oif];
    dfq = control->FreqDif[oif];
    outRow->Delay1[oif]  = inRow->Delay1[iif];
    outRow->Rate1[oif]   = inRow->Rate1[iif];
    outRow->Weight1[oif] = inRow->Weight1[iif];
    outRow->RefAnt1[oif] = inRow->RefAnt1[iif];
    /* correct for delay if given*/
    if ((inRow->Real1[iif]!=fblank) && (inRow->Weight1[iif]>0.0) &&
	(fabs(inRow->Delay1[iif])>0.0))  {
          dphs = 2*G_PI* dfq * inRow->Delay1[iif];
	  xr = inRow->Real1[iif]*cos(dphs) - inRow->Imag1[iif]*sin(dphs);
	  xi = inRow->Real1[iif]*sin(dphs) + inRow->Imag1[iif]*cos(dphs);
	  outRow->Real1[oif] = xr; outRow->Imag1[oif] = xi; 
    } else {outRow->Real1[oif] = inRow->Real1[iif]; 
            outRow->Imag1[oif] = inRow->Imag1[iif];}
    /* Two Stokes? */
    if ((control->inNumStok>1) && (control->outNumStok>1)) {
      outRow->Delay2[oif]  = inRow->Delay2[iif];
      outRow->Rate2[oif]   = inRow->Rate2[iif];
      outRow->Weight2[oif] = inRow->Weight2[iif];
      if ((inRow->Real2[iif]!=fblank) && (inRow->Weight2[iif]>0.0) &&
	  (fabs(inRow->Delay2[iif])>0.0)) {
	dphs = 2*G_PI* dfq * inRow->Delay2[iif];
	xr = inRow->Real2[iif]*cos(dphs) - inRow->Imag2[iif]*sin(dphs);
	xi = inRow->Real2[iif]*sin(dphs) + inRow->Imag2[iif]*cos(dphs);
	outRow->Real2[oif] = xr; outRow->Imag2[oif] = xi; 
      } else {outRow->Real2[oif] = inRow->Real2[iif]; 
              outRow->Imag2[oif] = inRow->Imag2[iif];} 
    }
  } /* end IF loop */

} /* end doClos */
