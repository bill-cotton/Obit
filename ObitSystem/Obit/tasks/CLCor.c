/* $Id$  */
/* Obit Radio interferometry calibration software                     */
/* applies user-selected corrections to the calibration CL table      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012                                               */
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
#include "ObitTableCL.h"
#include "ObitTableSUUtil.h"
#include "ObitTableANUtil.h"

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
  /** Subarray */
  olong subA;
   /**  blanking value */
  ofloat fblank;
 /**  CLCorParm */
  ofloat CLCorParm[30];
  /* Antenna list */
  ObitAntennaList *AntList;
  /* Source list */
  ObitSourceList *SrcList;
  /* index in Source list */
  olong SrcIndex;
}  ControlInfo; 
  
/* internal prototypes */
/* Get inputs */
ObitInfoList* CLCorIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void CLCorOut (ObitInfoList* outList, ObitErr *err);
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
void CLCorHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Modify CL table */
void CLCorDoCor (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Modification functions */
void doPANG (ControlInfo *control, ObitTableCLRow *row, ObitErr *err);
void doMULA (ControlInfo *control, ObitTableCLRow *row, ObitErr *err);
void doPNEG (ControlInfo *control, ObitTableCLRow *row, ObitErr *err);


/* Program globals */
gchar *pgmName = "CLCor";       /* Program name */
gchar *infile  = "CLCor.in" ;   /* File with program inputs */
gchar *outfile = "CLCor.out";   /* File to contain program outputs */
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
/*   Obit task applies user-selected corrections to the calibration CL table */
/*------------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitUV       *inData=NULL;
  ObitSystem   *mySystem= NULL;
  ObitErr      *err= NULL;

  /* Startup - parse command line */
  err = newObitErr();
  myInput = CLCorIn (argc, argv, err);
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
   CLCorDoCor (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Write history */
  CLCorHistory (myInput, inData, err); 
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

ObitInfoList* CLCorIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "CLCorIn";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return list;

  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* Initialize output */
  myOutput = defaultOutputs(err);
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

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

  return list;
} /* end CLCorIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: CLCor -input file -output ofile [args]\n");
    fprintf(stderr, "Obit task to modify an AIPS CL table \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def CLCor.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def CLCor.out\n");
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
  strTemp = "CLCor.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "CLCorName";
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
/*       ObitUV with CL table to modify                                   */
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
/*  Write History for CLCor                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void CLCorHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "souCode", "Qual", "EditStokes", "BIF", "EIF", 
    "FreqID", "timeRange",  "subA", "Antennas", 
    "calIn", "calOut", "corMode", "CLCParm", 
    NULL};
  gchar *routine = "CLCorHistory";

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
 
} /* end CLCorHistory  */

/*----------------------------------------------------------------------- */
/*  Apply corrections to CL table                                         */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV with CL table to modify                          */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void CLCorDoCor (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitTableAN     *ANTable  = NULL;
  ObitTableSU     *SUTable  = NULL;
  gchar corMode[5];
  oint numpol, numif, numterm;
  ObitInfoType type;
  olong iCLver, oCLver, inCLRow, outCLRow, ver, numPCal, numOrb;
  ofloat fblank = ObitMagicF();
  odouble freq;
  olong i, j, itemp, nif, nstok, count;
  gint32 dim[MAXINFOELEMDIM];
  ObitTableCL *iCLTable=NULL, *oCLTable=NULL;
  ObitTableCLRow *inRow=NULL, *outRow=NULL;
  ObitIOCode retCode;
  ControlInfo control;
  ObitUVSel *sel;
  gboolean wanted;
  gchar *routine = "CLCorDoCor";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Get Operation type */
  ObitInfoListGet(myInput, "corMode", &type, dim, corMode, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Get CL tables */
  itemp = 0;
  ObitInfoListGetTest(myInput, "calIn", &type, dim, &itemp);
  iCLver = itemp;
  itemp = 0;
  ObitInfoListGetTest(myInput, "calOut", &type, dim, &itemp);
  oCLver = itemp;
  if (oCLver<=0) oCLver = MAX(2, iCLver);
  numpol  = 0;
  numif   = 0;
  numterm = 0;
  iCLTable = newObitTableCLValue("in CL table", (ObitData*)inData, &iCLver, 
				OBIT_IO_ReadWrite, numpol, numif, numterm, err);
  /* Tell user */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Updating CL %d corMode %s to CL %d",iCLver, corMode, oCLver);

  /* Open */
  ObitTableCLOpen (iCLTable, OBIT_IO_ReadWrite,  err);
  if (err->error) goto cleanup;

  /* Set rows */
  inRow   = newObitTableCLRow (iCLTable);

  /* Input different from output? */
  if (iCLver!=oCLver) {
    numpol   = iCLTable->numPol;
    numif    = iCLTable->numIF;
    numterm  = iCLTable->numTerm;
    oCLTable = newObitTableCLValue("out CL table", (ObitData*)inData, &oCLver, 
				   OBIT_IO_ReadWrite, numpol, numif,numterm,  err);
    ObitTableCLOpen (oCLTable, OBIT_IO_ReadWrite,  err);
    if (err->error) goto cleanup;
    outRow  = newObitTableCLRow (oCLTable);
    ObitTableCLSetRow (oCLTable, outRow, err);
    if (err->error) goto cleanup;
  } else {  /* Same */
    oCLTable = ObitTableCLRef(iCLTable);
    outRow   = ObitTableCLRowRef(inRow);
  }

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
  for (i=0; i<30; i++) control.CLCorParm[i] = 0.0;
  ObitInfoListGetTest(myInput,  "CLCParm", &type, dim, control.CLCorParm);
  control.subA = 1;
  ObitInfoListGetTest(myInput,  "subA", &type, dim, &control.subA);
  control.SrcList  = NULL;
  control.SrcIndex = -999;

  /* Open and close UV data to update selector */
  ObitUVOpen(inData, OBIT_IO_ReadCal, err);
  ObitUVClose(inData, err);
  if (err->error) goto cleanup;
  sel  = inData->mySel;

  /* Any corMode specific checks or processing */
  if (!strncmp (corMode, "PANG", 4)) {
    /* Get antenna list */
    ver = control.subA;
    if (ver<=0) ver = 1;
    numPCal  = 0;
    numOrb   = 0;
    ANTable = newObitTableANValue (inData->name, (ObitData*)inData, &ver, 
				   numif, numOrb, numPCal, OBIT_IO_ReadOnly, err);
    if (ANTable) control.AntList = ObitTableANGetList (ANTable, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    ANTable = ObitTableANUnref(ANTable);
    if ((control.AntList==NULL)) {
      Obit_log_error(err, OBIT_Error,
		     "%s: Problem with antenna table %d ", 
		     routine, control.subA);;
      return;
    }
    /* Get antenna list */
    ver = 1;
    numPCal  = 0;
    numOrb   = 0;
    SUTable = newObitTableSUValue (inData->name, (ObitData*)inData, &ver, 
				   numif, OBIT_IO_ReadOnly, err);
    if (SUTable) control.SrcList = ObitTableSUGetList (SUTable, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    SUTable = ObitTableSUUnref(SUTable);
    if ((control.SrcList==NULL)) {
      Obit_log_error(err, OBIT_Error,
		     "%s: Problem with source list", routine);
      return;
    }
  } /* End setup for PANG */

  /* Loop over table */
  count = 0;
  for (inCLRow=1; inCLRow<=iCLTable->myDesc->nrow; inCLRow++) {
    retCode = ObitTableCLReadRow (iCLTable, inCLRow, inRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    if (inRow->status==-1) continue;

    /* Need new source info? */
    if ((control.SrcList!=NULL) && 
	((control.SrcIndex<0) || 
	 (control.SrcList->SUlist[control.SrcIndex]->SourID!=inRow->SourID))) {
      control.SrcIndex = -999;
      for (j=0; j<control.SrcList->number; j++) {
	if (control.SrcList->SUlist[j]->SourID==inRow->SourID) {
	  control.SrcIndex = j;
	  break;
	}
      }
      /* Better have it */
      if (control.SrcIndex<0) {
	Obit_log_error(err, OBIT_Error,
		       "%s: Could not locate source ID %d", routine, inRow->SourID);
	return;
      }
    } /* end get new source */

    /* This one selected? */
    wanted = (inRow->Time>=sel->timeRange[0]) && (inRow->Time<=sel->timeRange[1]);
    wanted = wanted && ObitUVSelWantSour (sel, inRow->SourID);
    wanted = wanted && ObitUVSelWantAnt (sel, inRow->antNo);
    wanted = wanted && ((inRow->FreqID==sel->FreqID) || (sel->FreqID<=0));
    wanted = wanted && ((inRow->SubA==sel->SubA) || (sel->SubA<=0));
    if (!wanted) continue;

    /* Copy row data */
    if (iCLver!=oCLver) {
      outRow->Time   = inRow->Time;
      outRow->TimeI  = inRow->TimeI;
      outRow->SourID = inRow->SourID;
      outRow->antNo  = inRow->antNo;
      outRow->SubA   = inRow->SubA;
      outRow->FreqID = inRow->FreqID;
      outRow->IFR    = inRow->IFR;
      outRow->atmos  = inRow->atmos;
      outRow->Datmos = inRow->Datmos;
      outRow->status = inRow->status;
      for (j=0; j<iCLTable->numTerm; j++) 
	outRow->GeoDelay[j] = inRow->GeoDelay[j];
      outRow->MBDelay1  = inRow->MBDelay1;
      outRow->clock1    = inRow->clock1;
      outRow->Dclock1   = inRow->Dclock1;
      outRow->dispers1  = inRow->dispers1;
      outRow->Ddispers1 = inRow->Ddispers1;
      for (j=0; j<iCLTable->numIF; j++) {
	outRow->Real1[j]   = inRow->Real1[j];
	outRow->Imag1[j]   = inRow->Imag1[j];
	outRow->Rate1[j]   = inRow->Rate1[j];
	outRow->Delay1[j]  = inRow->Delay1[j];
	outRow->Weight1[j] = inRow->Weight1[j];
	outRow->RefAnt1[j] = inRow->RefAnt1[j];
      }
      if (iCLTable->numPol>1) {
	outRow->MBDelay2  = inRow->MBDelay2;
	outRow->clock2    = inRow->clock2;
	outRow->Dclock2   = inRow->Dclock2;
	outRow->dispers2  = inRow->dispers2;
	outRow->Ddispers2 = inRow->Ddispers2;
	for (j=0; j<iCLTable->numIF; j++) {
	  outRow->Real2[j]   = inRow->Real2[j];
	  outRow->Imag2[j]   = inRow->Imag2[j];
	  outRow->Rate2[j]   = inRow->Rate2[j];
	  outRow->Delay2[j]  = inRow->Delay2[j];
	  outRow->Weight2[j] = inRow->Weight2[j];
	  outRow->RefAnt2[j] = inRow->RefAnt2[j];
	}
      } /* end 2 poln copy */
    } /* End copy if needed */

    /* Do operation */
    if (!strncmp (corMode, "PANG", 4)) {
      doPANG (&control, outRow, err);
    } else if (!strncmp (corMode, "MULA", 4)) {
      doMULA (&control, outRow, err);
    } else if (!strncmp (corMode, "PNEG", 4)) {
      doPNEG (&control, outRow, err);
    } else {  /* Unknown */
      Obit_log_error(err, OBIT_Error,
		     "%s: Unknown corMode: %s", routine, corMode);
      return;
    } /* End of cal by mode */
    if (err->error) goto cleanup;

    /* Rewrite row */
    if (iCLver!=oCLver) outCLRow = -1;
    else                outCLRow = inCLRow;
    count++;
    retCode = ObitTableCLWriteRow (oCLTable, outCLRow, outRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  } /* end loop over table */
  
  /* Tell user */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Modified %d of %d records", count, iCLTable->myDesc->nrow);
  /* Close table(s) */
  ObitTableCLClose (iCLTable, err);
  if (iCLver!=oCLver) ObitTableCLClose (oCLTable, err);
  
    /* Cleanup */
 cleanup:
  iCLTable = ObitTableCLUnref(iCLTable);
  oCLTable = ObitTableCLUnref(oCLTable);
  inRow = ObitTableCLRowUnref(inRow);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
} /* end CLCorDoCor  */

/* Modification functions */
/*----------------------------------------------------------------------- */
/*  Apply PANG corrections to CL table                                    */
/*  Average selected fringe rates in a record                             */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      CLRow       CL table row to modify                                */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doPANG (ControlInfo *control, ObitTableCLRow *CLRow, ObitErr *err)
{
  olong iif, istoke, jif;
  ofloat ParAng, rotR, rotI, tr, ti, fblank;
  if (err->error) return;

  fblank = control->fblank;

  /* Get parallactic angle */
  ParAng = ObitAntennaListParAng(control->AntList, CLRow->antNo, CLRow->Time, 
				 control->SrcList->SUlist[control->SrcIndex]);

  if (control->CLCorParm[0]>0.0) {
    /* Add corrections */
    rotR = cos(ParAng);
    rotI = sin(-ParAng);
  } else {
    /* remove corrections */
    rotR = cos(ParAng);
    rotI = sin(ParAng);
  }

  /* Loop over IFs modifying */
  for (iif=control->BIF; iif<=control->EIF; iif++) {
    jif = iif - 1;
    for (istoke=control->bStoke; istoke<=control->eStoke; istoke++) {
	if (istoke==1) {  /* First pol */
	  if ((CLRow->Real1[jif]!=fblank) && (CLRow->Imag1[jif]!=fblank) &&
	      CLRow->Weight1[jif]>0.0) {
	    tr = CLRow->Real1[jif];
	    ti = CLRow->Imag1[jif];
	    CLRow->Real1[jif] = tr*rotR - ti*rotI;
	    CLRow->Imag1[jif] = tr*rotI + ti*rotR;
	  }
	} else {       /* second pol - opposite phase */
	  if ((CLRow->Real2[jif]!=fblank) && (CLRow->Imag2[jif]!=fblank) &&
	      (CLRow->Weight2[jif]>0.0)) {
	    tr = CLRow->Real2[jif];
	    ti = CLRow->Imag2[jif];
	    CLRow->Real2[jif] = tr*rotR + ti*rotI;
	    CLRow->Imag2[jif] = -tr*rotI + ti*rotR;
	  }
	}
    } /* end loop over Stokes */
  } /* end loop over IF */

} /* end doPANG */

/*----------------------------------------------------------------------- */
/*  Apply MULA corrections to CL table                                    */
/*  Cause the amplitudes of the selected complex gains to be multiplied   */
/*  by CLCorParm[0].                                                      */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      CLRow     CL table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doMULA (ControlInfo *control, ObitTableCLRow *CLRow, ObitErr *err)
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
	  if ((CLRow->Real1[jif]!=fblank) && (CLRow->Imag1[jif]!=fblank)) {
	    CLRow->Real1[jif] *= control->CLCorParm[0];
	    CLRow->Imag1[jif] *= control->CLCorParm[0];
	  } /* end not blanked */
	} else {       /* second pol */
	  if ((CLRow->Real2[jif]!=fblank) && (CLRow->Imag2[jif]!=fblank)) {
	    CLRow->Real2[jif] *= control->CLCorParm[0];
	    CLRow->Imag2[jif] *= control->CLCorParm[0];
	  } /* end not blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doMULA */

/*----------------------------------------------------------------------- */
/*  Apply PNEG corrections to CL table                                    */
/*  Flip the sign of the gain phase for all selected CL solutions.        */
/*   Input:                                                               */
/*      control   Contol information                                      */
/*      CLRow     CL table row to modify                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doPNEG (ControlInfo *control, ObitTableCLRow *CLRow, ObitErr *err)
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
	  if (CLRow->Imag1[jif]!=fblank) {
	    CLRow->Imag1[jif] = -CLRow->Imag1[jif];
	  } /* end not blanked */
	} else {       /* second pol */
	  if (CLRow->Imag2[jif]!=fblank) {
	    CLRow->Imag2[jif] = -CLRow->Imag2[jif];
	  } /* end not blanked */
	}
    } /* end loop over Stokes */
  } /* end loop over IF */
} /* end doPNEG */

