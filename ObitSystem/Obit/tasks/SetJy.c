/* $Id$  */
/* Obit Radio interferometry calibration software                     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2007-2025                                          */
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
#include "ObitUV.h"
#include "ObitData.h"
#include "ObitTableSU.h"
#include "ObitTableSUUtil.h"
#include "ObitTableFQ.h"
#include "ObitTableFQUtil.h"
#define VELIGHT 2.997924562e8

/* internal prototypes */
/* Get inputs */
ObitInfoList* SetJyIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void SetJyOut (ObitInfoList* outList, ObitErr *err);
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
/* Update AIPS SU table */
ObitTableSU* SetJyUpdate (ObitUV* inData, ObitErr* err);
/* Write history */
void SetJyHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Calculate source flux density */
void CalcFlux (ObitTableSURow* row, ofloat *Parms, odouble Freq,
	       ofloat Flux[4], ObitErr* err);
/* Calculate source flux density for Perley & Butler 2017 */
void CalcFluxPB17 (ObitTableSURow* row, ofloat *Parms, odouble Freq,
		   ofloat Flux[4], ObitErr* err);


/* Program globals */
gchar *pgmName = "SetJy";       /* Program name */
gchar *infile  = "SetJy.in" ;   /* File with program inputs */
gchar *outfile = "SetJy.out";   /* File to contain program outputs */
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
/*   Obit task to apply set values in an AIPS SU table                    */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitErr      *err= NULL;
  ObitTableSU  *SUTable = NULL;
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "Qual", "calCode", "OPType", "BIF", "EIF", "ZeroFlux", 
    "Alpha", "FreqID", "RestFreq", "SysVel", "VelType", "VelDef", "Parms",
    NULL};

  /* Startup - parse command line */
  err = newObitErr();
  myInput = SetJyIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Get input uvdata */
  inData = getInputData (myInput, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Copy selection/calibration info to data */
  ObitInfoListCopyList (myInput, inData->info, dataParms);

  /* Update SU table */
  SUTable =  SetJyUpdate (inData, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Write history */
  SetJyHistory (myInput, inData, err); 
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* show any messages and errors */
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* cleanup */
  SUTable   = ObitTableSUUnref(SUTable);
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput  = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* SetJyIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "SetJyIn";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) {return list;}

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
} /* end SetJyIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: SetJy -input file -output ofile [args]\n");
    fprintf(stderr, "Obit task to modify AIPS SU table \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def SetJy.in\n");
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
/*     Alpha     float      Spectral index=0                              */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat Alpha = 0.0;
  /*gfloat farray[3];*/
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
  strTemp = "SetJy.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "SetJyName";
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
    
   /* Spectral index */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "Alpha", OBIT_float, dim, &Alpha, err);
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
  /*ObitInfoType type;
    gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
    gboolean     doCalSelect;*/
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
 
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
  gchar        Aname[13], Aclass[7], *Atype = "UV";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "FreqID", "Qual", "BIF", "EIF",
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
/*  Write History for SetJy                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void SetJyHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "calCode", "Qual", "OPType", "BIF", "EIF",
    "FreqID", "ZeroFlux",  "Alpha", "SysVel", "RestFreq", 
    "VelType", "VelDef",  "Parms", 
    NULL};
  gchar *routine = "SetJyHistory";

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
 
} /* end SetJyHistory  */

/*----------------------------------------------------------------------- */
/*  Update AIPS SU table                                                  */
/*   Input:                                                               */
/*      inData    ObitUV with SU table to update                          */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Returns AIPS SU table                                                */
/*----------------------------------------------------------------------- */
ObitTableSU* SetJyUpdate (ObitUV* inData, ObitErr* err)
{
  ObitTableSU *outSU=NULL;
  ObitTableFQ *FQTab=NULL;
  ObitTableSURow *row=NULL;
  ObitUVDesc *desc;
  olong ver;
  oint fqid, numIF, *sideBand=NULL;
  ofloat *chBandw=NULL, Alpha=0.0, SIFact;
  odouble Freq, nux, *freqOff=NULL;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong i, j, irow, lsou, nsou, iif, bif, eif;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *Sources=NULL, *calCode=NULL, *OPType=NULL, *VelType=NULL, *VelDef=NULL;
  ofloat SysVel=0.0, ZeroFlux[4]={0.0,0.0,0.0,0.0}, *Parms=NULL;
  ofloat altrfp, refpix, centpix, velinc, IFlux=0.0;
  olong Qual, lsign, ncheck, FreqID=0, BIF=1, EIF=0;
  odouble RestFreq=0.0, refFreq=0.0;
  gboolean wanted, match=FALSE;
  gchar tempName[101]; /* should always be big enough */
  gchar *blank = "        ";
  gchar *routine = "SetJyUpdate";

  /* error checks */
  if (err->error) return outSU;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Get control parameters */
  desc    = inData->myDesc;
  refpix  = desc->crpix[desc->jlocf];
  refFreq = desc->crval[desc->jlocf];
  centpix = 1.0 + desc->inaxes[desc->jlocf]/2.0;
  ObitInfoListGetP(inData->info, "Sources",  &type, dim, (gpointer)&Sources);
  /* Count number of actual sources */
  lsou = dim[0];
  nsou = 0;
  for (i=0; i<dim[1]; i++) {
    if (Sources[i*lsou]!=' ') nsou = i+1;
    else break;
  }
  ObitInfoListGetP(inData->info, "calCode",  &type, dim, (gpointer)&calCode);
  ObitInfoListGetP(inData->info, "OPType",   &type, dim, (gpointer)&OPType);
  ObitInfoListGetP(inData->info, "VelType",  &type, dim, (gpointer)&VelType);
  ObitInfoListGetP(inData->info, "VelDef",   &type, dim, (gpointer)&VelDef);
  ObitInfoListGetP(inData->info, "Parms",    &type, dim, (gpointer)&Parms);
  ObitInfoListGetTest(inData->info, "Qual",      &type, dim, &Qual);
  ObitInfoListGetTest(inData->info, "BIF",       &type, dim, &BIF);
  ObitInfoListGetTest(inData->info, "EIF",       &type, dim, &EIF);
  ObitInfoListGetTest(inData->info, "SysVel",    &type, dim, &SysVel);
  ObitInfoListGetTest(inData->info, "ZeroFlux",  &type, dim, ZeroFlux);
  ObitInfoListGetTest(inData->info, "Alpha",     &type, dim, &Alpha);
  ObitInfoListGetTest(inData->info, "FreqID",    &type, dim, &FreqID);
  ObitInfoListGetTest(inData->info, "RestFreq",  &type, dim, &RestFreq);
  if (calCode==NULL) calCode = blank;
  if (OPType==NULL)  OPType  = blank;
  if (VelType==NULL) VelType = blank;
  if (VelDef==NULL)  VelDef  = blank;
  IFlux = ZeroFlux[0];   /* Save */

  /* Digest velocity definition */
  lsign = 1;
  if (!strncmp (VelDef, "RADI", 4)) lsign = -1;
  if (!strncmp (VelDef, "OPTI", 4)) lsign =  1;
  if (Parms[0]==0.0) Parms[0] = 1.0;  /* Reference pixel */
  if (Parms[2]==0.0) Parms[2] = 1.0;  /* Flux correction factor */

  /* Get frequency Info */
  ver = 1;
  numIF = 0;
  FQTab = newObitTableFQValue (inData->name, (ObitData*)inData, &ver, 
			       OBIT_IO_ReadOnly, numIF, err);
  fqid = MAX (1, FreqID);
  ObitTableFQGetInfo (FQTab, fqid, &numIF, &freqOff, &sideBand, &chBandw, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outSU);

  /* Get AIPS SU table */
  outSU = newObitTableSUValue (inData->name, (ObitData*)inData, &ver, 
			       OBIT_IO_ReadWrite, numIF, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outSU);
 
  /* Open table */
  retCode = ObitTableSUOpen (outSU, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* IF range */
  numIF = outSU->numIF;
  bif = MAX (1, BIF);
  eif = MIN (EIF, numIF);
  if (eif<=0) eif = numIF;

  /* Create table row */
  row = newObitTableSURow (outSU);
  ObitTableSUSetRow (outSU, row, err);
  if (err->error) goto cleanup;

  /* Velocity type */
  if (strncmp (VelType, "        ", 8)) strncpy (outSU->velType, VelType, 8);
  
  /* Velocity definition */
  if (strncmp (VelDef, "        ", 8)) strncpy (outSU->velDef, VelDef, 8);

  /* FreqID */
  outSU->FreqID = fqid;
  
  /* Loop over table */
  ncheck = MIN (lsou, outSU->myDesc->repeat[outSU->SourceCol]);
  irow = 0;
  retCode = OBIT_IO_OK;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableSUReadRow (outSU, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

    /* Is this one wanted */
    wanted = FALSE;
    if (nsou<=0) wanted = (Qual<0) || (Qual==row->Qual);  /* All selected */
    else { /* Check list */

      for (i=0; i<nsou; i++) {
	/* Is this a match? */
	match = ObitStrCmp(&Sources[i*lsou], row->Source, ncheck);
	if (match) break;
      } /* end check is source desired */
      wanted = match && ((Qual<0) || (Qual==row->Qual));
    } /* end Source names given */

    /* Update this one? */
    if (wanted) {

      for (j=0; j<lsou; j++) {tempName[j] = ' ';} tempName[j] = 0;
	/* get blank padded name */
	for (j=0; j<lsou; j++) {
	  if (row->Source[j]==0) break;  /* only values in string */
	  tempName[j] = row->Source[j]; 
	}

      /* Update flux */
      for (iif=bif; iif<=eif; iif++) {
	Freq = desc->freq + freqOff[iif-1] + row->FreqOff[iif-1];
	/* Spectral index factor at IF ref Freq*/
	if (Alpha==0.0) SIFact = 1.0;
	else            SIFact = pow(Freq/refFreq, Alpha);
	/* Offset to IF center */
	Freq += (centpix-refpix) * chBandw[iif] * sideBand[iif];
	/* Calculate flux density? */
	if (!strncmp (OPType, "CALC", 4)) {
	  CalcFlux (row, Parms, Freq, ZeroFlux, err);
	  if (err->error) goto cleanup;
	} else if (IFlux>0.0) {  /* Calculate from IFlux */
	  ZeroFlux[0] = IFlux * SIFact;
	}
	if (ZeroFlux[0]>=0.0) row->IFlux[iif-1] = ZeroFlux[0];
	row->QFlux[iif-1] = ZeroFlux[1]*SIFact;
	row->UFlux[iif-1] = ZeroFlux[2]*SIFact;
	row->VFlux[iif-1] = ZeroFlux[3]*SIFact;
	/* Message */
	Obit_log_error(err, OBIT_InfoErr,
		       "%s %c IF %d Flux %8.3f %8.3f %8.3f %8.3f",
		       tempName, calCode[0], iif,  ZeroFlux[0], 
		       ZeroFlux[1]*SIFact, ZeroFlux[2]*SIFact, 
		       ZeroFlux[3]*SIFact);
      }
      /* Or perhaps just reset everything */
      if (!strncmp (OPType, "REJY", 4)) {
	for (iif=bif; iif<=eif; iif++) {
	  if (ZeroFlux[0]>=0.0) row->IFlux[iif-1] = ZeroFlux[0];
	  row->QFlux[iif-1] = ZeroFlux[1];
	  row->UFlux[iif-1] = ZeroFlux[2];
	  row->VFlux[iif-1] = ZeroFlux[3];
	}
      }

      /* CalCode */
      if (strncmp (calCode, "    ", 4))  strncpy (row->CalCode, calCode, 4);
      if (!strncmp (calCode, "----", 4)) strncpy (row->CalCode, blank, 4);

      /* RestFreq */
      if (RestFreq>0.0) {
	for (iif=bif; iif<=eif; iif++) 
	  row->RestFreq[iif-1] = RestFreq;
      }

      /* Systemic velocity */
      if (SysVel!=0.0) {
	for (iif=bif; iif<=eif; iif++) {
	  altrfp = Parms[0];
	  /* Calculate signed freq. inc */
	  if (sideBand[iif]==0) sideBand[iif] = 1;
	  nux = desc->freq + freqOff[iif-1] + (altrfp-refpix) * chBandw[iif-1]*sideBand[iif-1] + 
	    row->FreqOff[iif-1];
	  velinc = -((VELIGHT + lsign*SysVel*1.0e3) * chBandw[iif-1]*sideBand[iif-1]) / nux;
	  row->LSRVel[iif-1] = SysVel*1.0e3 + velinc * (refpix-altrfp);
	}
      }
      if ((!strncmp (OPType, "REVL", 4)) || (!strncmp (OPType, "RESE", 4)))
	for (iif=bif; iif<=eif; iif++) row->LSRVel[iif-1] = 0.0;

      /* Update in table */
      retCode = ObitTableSUWriteRow (outSU, irow, row, err);
      if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    } /* end update source entry */

  } /* end loop over table */


  /* Cleanup */
 cleanup:
  /* Close table */
  retCode = ObitTableSUClose (outSU, err);
  row = ObitTableSURowUnref (row);
  FQTab = ObitTableFQUnref (FQTab);
  if (sideBand) g_free(sideBand);
  if (chBandw)  g_free(chBandw);
  if (freqOff)  g_free(freqOff);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outSU);
  return outSU;
} /* end SetJyUpdate  */

/*----------------------------------------------------------------------- */
/*  Calculate calibrator flux density                                     */
/*  Recognizes 3C48,3C138,3C147,3C295,3C286 and 1934-638 (from R. Duncan) */
/*  0408-65 (From T.Mauch/Ben Hugo                                        */
/*  Flux density, S, at frequency, nu, is calculated assuming             */
/*                                                                        */
/*   Log  (S) = F  and Log  (nu) = x, and                                 */
/*      10                10                                              */
/*                                                                        */
/*   F  = [ a + b * x + c * x**2 + d * x**3]                              */
/*                                                                        */
/*  where  a, b, c and d are observed parameters.                         */
/* The error in F as a function of errors in a,b,c and d are              */
/*                                                                        */
/*    2      dF 2  2     dF 2  2     dF 2  2     dF 2  2                  */
/*   E   =  (--)  E   + (--)  E  +  (--)  E  +  (--)  E                   */
/*    F      da    a     db    b     dc    c     dd    d                  */
/*                                                                        */
/* where the dF/da are partial derivatives.  So                           */
/*                                                                        */
/*    2       2     2  2     4  2     6  2                                */
/*   E   =   E   + x  E  +  x  E  +  x  E                                 */
/*    F       a        b        c        d                                */
/*                                                                        */
/* The error in S is                                                      */
/*                                    F                                   */
/*           dS             dS      10                F dF                */
/*   E   =  (--)  E   and  (--) = d --  =  Log (10) 10  -- = S Log (10)   */
/*    S      dF    F        dF      dF        e         dF        e       */
/* so,                                                                    */
/*   E   = S Log (10) * E                                                 */
/*    S         e        F                                                */
/* Last Documented by Glen Langston on 93 March 2003                      */
/*  Adopted from the AIPSish SETJY/GETFLX                                 */
/*   Input:                                                               */
/*      row     SU table row describing source                            */
/*      Parms   User supplied parameters                                  */
/*       [1]: Only for 'CALC' option:                                     */
/*            <= 0 => use latest VLA values (1999.2) or, for 1934-638, the*/
/*               ATCA value of 30Jul94.                                   */
/*            1 => use Baars values or old ATCA/PKS values for 1934-638   */
/*            2 => use VLA 1995.2 values or for 1934-638 the ATCA value of*/
/*                 30Jul94.                                               */
/*            3 => use oldest VLA values (1990) or,for 1934-638, the      */
/*                 ATCA value of 30Jul94.                                 */
/*            4 => Perley & Butler 2012 3C123, 3C196, 3C286, 3C295        */
/*          >=5 => Perley & Butler 2017                                   */
/*        [2]: Only for 'CALC' option:                                    */
/*             multiply the calculated fluxes by Parms[2]                 */
/*      Freq    Frequency                                                 */
/*   Output:                                                              */
/*      Flux[0] Flux density for source at Freq                           */
/*      err     Obit Error/message stack                                  */
/*----------------------------------------------------------------------- */
void CalcFlux (ObitTableSURow* row, ofloat *Parms, odouble Freq,
	       ofloat Flux[4], ObitErr* err)
{
  olong ictype, isrc, i;
  odouble temp2=0.0, dt, ferror, freqm;
  /* number of recognized source names */
  olong xnsou = 45;

  /* Source lists, Baars et al. */
  ofloat coeff[9][4] = { 
    {1.480, 0.292, -0.124, 0.000},       /* 3C286 */
    {2.345, 0.071, -0.138, 0.000},       /* 3C48 */
    {1.766, 0.447, -0.184, 0.000},       /* 3C147 */
    {2.009, -0.07176, -0.0862, 0.000},   /* 3C138 */
    {-23.839, 19.569, -4.8168, 0.35836}, /* 1934-638 */
    {1.485, 0.759, -0.255, 0.000},        /* 3C295 */
    {1.8077,   -0.8018, -0.1157,   0.000},     /* 3C123 P&B 2012*/
    {1.2969,   -0.8690, -0.1788,   0.0305},     /* 3C196 P&B 2012*/
    {-113.51820049,  110.7808751,   -35.26819126,    3.70123262} /*0408-65 Mauch/Hugo */
  };

  /* Baars et al errors */
  ofloat coerr[9][4] = {
    {0.018, 0.006, 0.001, 0.000},  /* 3C286 */
    {0.030, 0.001, 0.001, 0.000},  /* 3C48 */
    {0.017, 0.006, 0.001, 0.000},  /* 3C147 */
    {0.000, 0.000, 0.000, 0.000},  /* 3C138 */
    {0.000, 0.000, 0.000, 0.000},  /* 1934-638 */
    {0.013, 0.009, 0.001, 0.000},  /* 3C295 */
    {0.000, 0.000, 0.000, 0.000},  /* 3C123 P&B 2012*/
    {0.000, 0.000, 0.000, 0.000},  /* 3C196 P&B 2012*/
    {0.000, 0.000, 0.000, 0.000}   /* 0408-65 Mauch/Hugo  */
  };

  /* Source lists, Taylor 1999.2
     Note that the VLA coefficients are for freq. in GHz, not MHz,
     not true of 1934-638.  This requires some modifications in
     the calculation loop for rcoeff */
  ofloat rcoeff[9][4] = {
    {1.23734, -0.43276, -0.14223,  0.00345},   /* 3C286 */
    {1.31752, -0.74090, -0.16708,  0.01525},   /* 3C48 */
    {1.44856, -0.67252, -0.21124,  0.04077},   /* 3C147 */
    {1.00761, -0.55629, -0.11134, -0.01460},   /* 3C138 */
    {-30.7667,  26.4908,  -7.0977, 0.605334},  /* 1934-638 (Reynolds, 02/Jul/94) */
    {1.46744, -0.77350, -0.25912,  0.00752},   /* 3C295 */
    {1.8077,   -0.8018, -0.1157,   0.000},     /* 3C123 P&B 2012*/
    {1.2969,   -0.8690, -0.1788,   0.0305},    /* 3C196 P&B 2012*/
    /*old {-113.51820049,  110.7808751,   -35.26819126,    3.70123262} 0408-65 Mauch/Hugo */
    {-41.59292326, 44.81497595, -15.18277173, 1.67166122} /*0408-65 Mauch/Hugo Mar21*/
  };
  /* Source lists, Perley 1990 */
  ofloat ncoeff[9][4] = {
    {1.35899, 0.35990, -0.13338, 0.000},     /* 3C286, */
    {2.0868, 0.20889, -0.15498, 0.000},      /* 3C48 */
    {1.92641, 0.36072, -0.17389, 0.000},     /* 3C147 */
    {2.009, -0.07176, -0.0862, 0.000},       /* 3C138 (Baars again) */
    {-30.7667, 26.4908, -7.0977, 0.605334},  /* 1934-638 (Reynolds, 02/Jul/94) */
    {1.485, 0.759, -0.255, 0.000},           /* 3C295 (Baars again) */
    {1.8077,   -0.8018, -0.1157,   0.000},   /* 3C123 P&B 2012*/
    {1.2969,   -0.8690, -0.1788,   0.0305},   /* 3C196 P&B 2012*/
    /*old {-113.51820049,  110.7808751,   -35.26819126,    3.70123262} 0408-65 Mauch/Hugo */
    {-41.59292326, 44.81497595, -15.18277173, 1.67166122} /*0408-65 Mauch/Hugo Mar21*/
 };
  /* Source lists, Perley 1995.2 */
  ofloat pcoeff[9][4] = {  /*  */
    {0.50344,  1.05026, -0.31666,  0.01602},   /* 3C286 */
    {1.16801,  1.07526, -0.42254,  0.02699},   /* 3C48 */
    {0.05702,  2.09340, -0.70760,  0.05477},   /* 3C147 */
    {1.97498, -0.23918,  0.01333, -0.01389},   /* 3C138 */
    {-30.7667, 26.4908, -7.0977,   0.605334},  /* 1934-638 (Reynolds, 02/Jul/94) */
    {1.28872,  0.94172, -0.31113,  0.00569},   /* 3C295  */
    {1.8077,   -0.8018, -0.1157,   0.000},     /* 3C123 P&B 2012*/
    {1.2969,   -0.8690, -0.1788,   0.0305},     /* 3C196 P&B 2012*/
    /*old??? {0.42884, 1.9395659, -0.66243187, 0.03926736} 0408-65 Mauch/Hugo */
    {-41.59292326, 44.81497595, -15.18277173, 1.67166122} /*0408-65 Mauch/Hugo Mar21*/
 };
  /* Source lists, Perley & Butler 2012 */
  ofloat pbcoeff[9][4] = {  /*  */
    {1.2515,  -0.4605,  -0.1715,   0.0336},    /* 3C286 P&B 2012 */
    {1.16801,  1.07526, -0.42254,  0.02699},   /* 3C48 Perley 1995.2 */
    {0.05702,  2.09340, -0.70760,  0.05477},   /* 3C147  Perley 1995.2 */
    {1.97498, -0.23918,  0.01333, -0.01389},   /* 3C138 Perley 1995.2 */
    {-30.7667,  26.4908,  -7.0977, 0.605334},  /* 1934-638 (Reynolds, 02/Jul/94) */
    {1.4866,   -0.7871, -0.3440,   0.0749},    /* 3C295 P&B 2012 */
    {1.8077,   -0.8018, -0.1157,   0.000},     /* 3C123 P&B 2012 */
    {1.2969,   -0.8690, -0.1788,   0.0305},    /* 3C196 P&B 2012 */
    /* old {-113.51820049,  110.7808751,   -35.26819126,    3.70123262} 0408-65 Mauch/Hugo */
    {-41.59292326, 44.81497595, -15.18277173, 1.67166122} /*0408-65 Mauch/Hugo Mar21*/
  };
  /* Source list - aliases */
  gchar *knosou[] = {
    "3C286",   "1328+307", "1331+305", "J1331+3030", "3c286", "3C 286",
    "3C48",    "0134+329", "0137+331", "J0137+3309", "3c48", "3C 48",
    "3C147",   "0538+498", "0542+498", "J0542+4951", "3c147", "3C 147",
    "3C138",   "0518+165", "0521+166", "J0521+1638", "3c138", "3C 138", "3C_138",
    "1934-638","1934-638", "1934-638", "J1939-6342", "PKS1934-638",
    "3C295",   "1409+524", "1411+522", "J1411+5212", "3c295", "3C 295", 
    "3C123",   "0433+295", "0437+296", "J0437+2940", "3c123", "3C 123", 
    "3C196",   "0809+483", "0813+482", "J0813+4822", "3c196", "3C 196",
    "0408-65" , "0408-65", "0408-65",  "0407-658" ,  "J0407-658",
  };
  /* Number of characters to check */
  olong lenchk[] = {
    5,8,8,10,5,6,   4,8,8,10,4,5, 5,8,8,10,5,6,
    5,8,8,10,5,6,6, 8,8,8,10,11,  5,8,8,10,5,6,
    5,8,8,10,5,6,   5,8,8,10,5,6, 7,7,7,8,8,
  };
  /* Source number in coef table */
  olong ksouno[] = {
    1,1,1,1,1,1,   2,2,2,2,2,2, 3,3,3,3,3,3,
    4,4,4,4,4,4,4, 5,5,5,5,5,   6,6,6,6,6,6,
    7,7,7,7,7,7,   8,8,8,8,8,8, 9,9,9,9,9,
  };
  /*  Frequency break points for bands 
  ofloat fband[] = {
    0.15e3, 0.7e3, 2.0e3, 6.0e3, 11.5e3, 18.e3, 28.e3
  };*/
  olong j;
  gchar tempName[101];

  /* Perley & Butler 2017 has more coefficients */
  if (Parms[1]>=4.9) {
    CalcFluxPB17 (row, Parms, Freq, Flux, err);
    return;
  }

  /* Name */
  for (j=0; j<16; j++) {tempName[j] = row->Source[j];} tempName[j] = 0;
  if (strlen(tempName)<16) {
    for (j=strlen(tempName); j<16; j++) tempName[j] = ' '; 
    tempName[j] = 0;
  }

  /* Type of calculation */
  ictype = (olong)(Parms[1]+0.5);

  /* Work out band 
  iband = 8;*/
  freqm = Freq*1.0e-6;  /* Frequency in MHz */
  /*for (i=0; i<7; i++) {
    if (freqm < fband[7-i]) iband = 7-i;
    }*/

  /* Lookup source */
  isrc = -1;
  for (i=0; i<xnsou; i++) {
    if (!strncmp(tempName, knosou[i], lenchk[i])) isrc = ksouno[i]-1;
  }

  /* Complain if not found */
  if (isrc<0) {
    Obit_log_error(err, OBIT_Error, "Unknown Flux calibrator %s", tempName);
    return;
  }

  /*  Compute flux */
  dt = log10 (freqm);
  /* All entries for 3C123 and 3C196 are from P&B 2012 with coefficients in GHz */
  if ((isrc==6) || (isrc==7)) dt -= 3.0e0;


  /* Assume no error */
  ferror = 0.0;
  /* Taylor 1999.2 &Reynolds (1934-638) */
  if ((ictype<=0) && (isrc!=4) && (isrc!=8)) {
    dt -= 3.0e0;
    temp2 = rcoeff[isrc][0] + dt * (rcoeff[isrc][1] + dt * (rcoeff[isrc][2] + dt * rcoeff[isrc][3]));
  }  
  if ((ictype<=0) && ((isrc==4)||(isrc==8))) {
    temp2 = rcoeff[isrc][0] + dt * (rcoeff[isrc][1] + dt * (rcoeff[isrc][2] + dt * rcoeff[isrc][3]));
  }

  /* Baars scale */
  if (ictype==1) {
    temp2 = coeff[isrc][0] + dt * (coeff[isrc][1] + dt *(coeff[isrc][2] + dt * coeff[isrc][3]));
    
    /* Baars error, sum of squares */
    ferror =  (coerr[isrc][0]*coerr[isrc][0]) +
      ((coerr[isrc][1]*dt)*(coerr[isrc][1]*dt)) +
      ((coerr[isrc][2]*dt*dt)*(coerr[isrc][2]*dt*dt)) +
      ((coerr[isrc][3]*dt*dt*dt)*(coerr[isrc][3]*dt*dt*dt));
  }

  /* Perley 1995.2 & Reynolds (1934-638) */
  if (ictype==2) {
    temp2 = pcoeff[isrc][0] + dt * (pcoeff[isrc][1] + dt * (pcoeff[isrc][2] + dt * pcoeff[isrc][3]));
  }

  /* Perley 1990 & Reynolds (1934-638) */
  if (ictype==3) {
    temp2 = ncoeff[isrc][0] + dt * (ncoeff[isrc][1] + dt * (ncoeff[isrc][2] + dt * ncoeff[isrc][3]));
  }
  /* Perley & Butler 2012 with others for completeness */
  if (ictype>=4) {
    /* coefficients for P&B 2012 for GHz, rest MHz, 3C123, 3C196 already changed. */
    if ((isrc==0) || (isrc==5)) dt -= 3.0e0;
    temp2 = pbcoeff[isrc][0] + dt * (pbcoeff[isrc][1] + dt * (pbcoeff[isrc][2] + dt * pbcoeff[isrc][3]));
  }
  
  Flux[0] = pow (10.0, temp2);
  
  /*  If non-zero error exponent */
  if (ferror>0.0) {
    /* sqrt sum of squares * factor */
    ferror = log(10.) * Flux[0] * sqrt(ferror);
  }
  
Flux[0] *= Parms[2];  /* Apply fudge factor */
} /* end CalcFlux */

/*----------------------------------------------------------------------- */
/*  Calculate calibrator flux density from Perley & Butler 2017           */
/*  Flux density, S, at frequency, nu, is calculated assuming             */
/*                                                                        */
/*   Log  (S) = F  and Log  (nu) = x, and                                 */
/*      10                10                                              */
/*                                                                        */
/*   F  = [ a + b * x + c * x**2 + d * x**3 + e * x**4 + f * x**5         */
/*                                                                        */
/*  where  a, b, c, d, e, f are observed parameters.                      */
/* The error in F as a function of errors in a,b,c and d are              */
/*                                                                        */
/*    2      dF 2  2     dF 2  2     dF 2  2     dF 2  2                  */
/*   E   =  (--)  E   + (--)  E  +  (--)  E  +  (--)  E                   */
/*    F      da    a     db    b     dc    c     dd    d                  */
/*                                                                        */
/* where the dF/da are partial derivatives.  So                           */
/*                                                                        */
/*    2       2     2  2     4  2     6  2                                */
/*   E   =   E   + x  E  +  x  E  +  x  E                                 */
/*    F       a        b        c        d                                */
/*                                                                        */
/* The error in S is                                                      */
/*                                    F                                   */
/*           dS             dS      10                F dF                */
/*   E   =  (--)  E   and  (--) = d --  =  Log (10) 10  -- = S Log (10)   */
/*    S      dF    F        dF      dF        e         dF        e       */
/* so,                                                                    */
/*   E   = S Log (10) * E                                                 */
/*    S         e        F                                                */
/* Last Documented by Glen Langston on 93 March 2003                      */
/*  Adopted from the AIPSish SETJY/GETFLX                                 */
/*   Input:                                                               */
/*      row     SU table row describing source                            */
/*      Parms   User supplied parameters                                  */
/*       [1]: Only for 'CALC' option:                                     */
/*            <= 0 => use latest VLA values (1999.2) or, for 1934-638, the*/
/*               ATCA value of 30Jul94.                                   */
/*            1 => use Baars values or old ATCA/PKS values for 1934-638   */
/*            2 => use VLA 1995.2 values or for 1934-638 the ATCA value of*/
/*                 30Jul94.                                               */
/*            >= 3 => use oldest VLA values (1990) or,for 1934-638, the   */
/*                 ATCA value of 30Jul94.                                 */
/*        [2]: Only for 'CALC' option:                                    */
/*             multiply the calculated fluxes by Parms[2]                 */
/*      Freq    Frequency                                                 */
/*   Output:                                                              */
/*      Flux[0] Flux density for source at Freq                           */
/*      err     Obit Error/message stack                                  */
/*
  source       alp_0     alp_1    alp_2     alp_3   alp_4   alp_5 chi^2   Freq range (GHz)
  J0133-3629  1.0440   -0.6619   -0.2252                            267     0.2-4
  3C48        1.3253   -0.7553   -0.1914   0.0498                   3.1    0.05 - 50
  Fornax A    2.2175   -0.6606                                       17    0.2 - 0.5
  3C123       1.8017   -0.7884   -0.1035  -0.0248   0.0090          1.9    0.05 - 50
  J0444+2809  0.9710   -0.8938   -0.1176                            3.3    0.2 - 2.0
  3C138       1.0088   -0.4981   -0.1552  -0.0102   0.0223          1.5    0.2 - 50
  Pictor A    1.9380   -0.7470   -0.0739  -0.0739                   8.1    0.2 - 4.0
  Taurus A    2.9516   -0.2173   -0.0473  -0.0074                   1.9    0.05 - 4.0
  3C147       1.4516   -0.6961   -0.2007   0.0640  -0.0464  0.0289  2.2    0.05 - 50
  3C196       1.2872   -0.8530   -0.1534  -0.0200   0.0201          1.6    0.05 - 50
  Hydra A     1.7795   -0.9176   -0.0843  -0.0139   0.0295          3.6    0.05 - 12
  Virgo A     2.4466   -0.8116   -0.0483                            2.0    0.05 - 50
  3C286       1.2482   -0.4507   -0.1798   0.0357                   1.9    0.05 - 50
  3C295       1.4701   -0.7658   -0.2780  -0.0347   0.0399          1.6    0.05 - 50
  Hercules A  1.8298   -1.0247   -0.0951                            2.3    0.2 - 12
  3C353       1.8627   -0.6938   -0.0998  -0.0732                   2.2    0.2 - 4
  3C380       1.2320   -0.7909    0.0947   0.0976  -0.1794  -0.1566 2.9    0.05 - 50
  Cygnus A    3.3498   -1.0022   -0.2246   0.0227   0.0425          1.9    0.05 - 12
  3C444       1.1064   -1.0052   -0.0750  -0.0767                   5.7    0.2 - 12
  Cassiopeia A 3.3584  -0.7518   -0.0347  -0.0705                   2.1    0.2 -4.
Errors
source         0       1      2       3      4       5
J0133-3629  0.0010  0.0018  0.0063
3C48        0.0005  0.0009  0.0011  0.0009
Fornax A    0.0030  0.0063
3C123       0.0007  0.0012  0.0023  0.0013  0.0013
J0444+2809  0.0011  0.0040  0.0096
3C138       0.0009  0.0022  0.0030  0.0065  0.0031
Pictor A    0.0010  0.0013  0.0055
Taurus A    0.0010  0.0032  0.0052  0.0132
3C147       0.0010  0.0017  0.0050  0.0044  0.0035  0.0025
3C196       0.0007  0.0012  0.0023  0.0013  0.0013
Hydra A     0.0009  0.0012  0.0041  0.0014  0.0028
Virgo A     0.0007  0.0020  0.0027
3C286       0.0005  0.0009  0.0011  0.0009
3C295       0.0007  0.0012  0.0023  0.0013  0.0013
Hercules A  0.0007  0.0012  0.0020
3C353       0.0010  0.0014  0.0052  0.0047
3C380       0.0016  0.0037  0.0218  0.0218  0.0597  0.0499
Cygnus A    0.0010  0.0014  0.0055  0.0021  0.0045
3C444       0.0009  0.0020  0.0039  0.0051
CassiopeiaA 0.0010  0.0014  0.0052  0.0047
*/
/*----------------------------------------------------------------------- */
void CalcFluxPB17 (ObitTableSURow* row, ofloat *Parms, odouble Freq,
		   ofloat Flux[4], ObitErr* err)
{
  olong isrc, i;
  odouble temp2=0.0, dt, ferror, freqg;

  /* Source lists */
 ofloat coeff[21][6] = { 
   {1.0440,  -0.6619,  -0.2252,  0.0   ,  0.0   , 0.0}, /*   J0133-3629  */
   {1.3253,  -0.7553,  -0.1914,  0.0498,  0.0   , 0.0}, /*   3C48        */
   {2.2175,  -0.6606,    0.0  ,  0.0   ,  0.0   , 0.0}, /*   Fornax A    */
   {1.8017,  -0.7884,  -0.1035, -0.0248,  0.0090, 0.0}, /*   3C123       */
   {0.9710,  -0.8938,  -0.1176,  0.0   ,  0.0   , 0.0}, /*   J0444+2809  */
   {1.0088,  -0.4981,  -0.1552, -0.0102,  0.0223, 0.0}, /*   3C138       */
   {1.9380,  -0.7470,  -0.0739, -0.0739,  0.0   , 0.0}, /*   Pictor A    */
   {2.9516,  -0.2173,  -0.0473, -0.0074,  0.0   , 0.0}, /*   Taurus A    */
   {1.4516,  -0.6961,  -0.2007,  0.0640, -0.0464, 0.0}, /*   3C147       */
   {1.2872,  -0.8530,  -0.1534, -0.0200,  0.0201, 0.0}, /*   3C196       */
   {1.7795,  -0.9176,  -0.0843, -0.0139,  0.0295, 0.0}, /*   Hydra A     */
   {2.4466,  -0.8116,  -0.0483,  0.0   ,  0.0   , 0.0}, /*   Virgo A     */
   {1.2482,  -0.4507,  -0.1798,  0.0357,  0.0   , 0.0}, /*   3C286       */
   {1.4701,  -0.7658,  -0.2780, -0.0347,  0.0399, 0.0}, /*   3C295       */
   {1.8298,  -1.0247,  -0.0951,  0.0   ,  0.0   , 0.0}, /*   Hercules A  */
   {1.8627,  -0.6938,  -0.0998, -0.0732,  0.0   , 0.0}, /*   3C353       */
   {1.2320,  -0.7909,   0.0947,  0.0976, -0.1794,-0.1}, /*   3C380       */
   {3.3498,  -1.0022,  -0.2246,  0.0227,  0.0425, 0.0}, /*   Cygnus A    */
   {1.1064,  -1.0052,  -0.0750, -0.0767,  0.0   , 0.0}, /*   3C444       */
   {3.3584,  -0.7518,  -0.0347  -0.0705,  0.0   , 0.0}, /*   Cassiopeia A*/
   /*old {0.42884, 1.9395659, -0.66243187, 0.03926736}  0408-65 Mauch/Hugo */
   {-41.59292326, 44.81497595, -15.18277173, 1.67166122} /*0408-65 Mauch/Hugo Mar21*/
};

  /* errors */
  ofloat coerr[21][6] = {
    {0.0010, 0.0018, 0.0063, 0.0   , 0.0   , 0.0},    /*  J0133-3629 */
    {0.0005, 0.0009, 0.0011, 0.0009, 0.0   , 0.0},    /*  3C48       */
    {0.0030, 0.0063, 0.0   , 0.0   , 0.0   , 0.0},    /*  Fornax A   */
    {0.0007, 0.0012, 0.0023, 0.0013, 0.0013, 0.0},    /*  3C123      */
    {0.0011, 0.0040, 0.0096, 0.0   , 0.0   , 0.0},    /*  J0444+2809 */
    {0.0009, 0.0022, 0.0030, 0.0065, 0.0031, 0.0},    /*  3C138      */
    {0.0010, 0.0013, 0.0055, 0.0   , 0.0   , 0.0},    /*  Pictor A   */
    {0.0010, 0.0032, 0.0052, 0.0132, 0.0   , 0.0},    /*  Taurus A   */
    {0.0010, 0.0017, 0.0050, 0.0044, 0.0035, 0.0025}, /*  3C147      */
    {0.0007, 0.0012, 0.0023, 0.0013, 0.0013, 0.0},    /*  3C196      */
    {0.0009, 0.0012, 0.0041, 0.0014, 0.0028, 0.0},    /*  Hydra A    */
    {0.0007, 0.0020, 0.0027, 0.0   , 0.0   , 0.0},    /*  Virgo A    */
    {0.0005, 0.0009, 0.0011, 0.0009, 0.0   , 0.0},    /*  3C286      */
    {0.0007, 0.0012, 0.0023, 0.0013, 0.0013},         /*  3C295      */
    {0.0007, 0.0012, 0.0020, 0.0   , 0.0},            /*  Hercules A */
    {0.0010, 0.0014, 0.0052, 0.0047, 0.0   , 0.0},    /*  3C353      */
    {0.0016, 0.0037, 0.0218, 0.0218, 0.0597, 0.0499}, /*  3C380      */
    {0.0010, 0.0014, 0.0055, 0.0021, 0.0045, 0.0},    /*  Cygnus A   */
    {0.0009, 0.0020, 0.0039, 0.0051, 0.0   , 0.0},    /*  3C444      */
    {0.0010, 0.0014, 0.0052, 0.0047, 0.0   , 0.0},    /*  Cassiopeia A */
    {0.000,  0.000,  0.000,  0.000,  0.0,    0.0}     /* 0408-65 Mauch/Hugo  */
 };

  /* Validity Frequency range */
  ofloat range[21][2] = { 
    {0.2,   4.0}, /*   J0133-3629    */
    {0.05, 50.0}, /*   3C48          */
    {0.2,   0.5}, /*   Fornax A      */
    {0.05, 50.0}, /*   3C123         */
    {0.2,   2.0}, /*   J0444+2809    */
    {0.2,  50.0}, /*   3C138         */
    {0.2,   4.0}, /*   Pictor A      */
    {0.05,  4.0}, /*   Taurus A      */
    {0.05, 50.0}, /*   3C147         */
    {0.05, 50.0}, /*   3C196         */
    {0.05, 12.0}, /*   Hydra A       */
    {0.05, 50.0}, /*   Virgo A       */
    {0.05, 50.0}, /*   3C286         */
    {0.05, 50.0}, /*   3C295         */
    {0.2,  12.0}, /*   Hercules A    */
    {0.2,   4.0}, /*   3C353         */
    {0.05, 50.0}, /*   3C380         */
    {0.05, 12.0}, /*   Cygnus A      */
    {0.2,  12.0}, /*   3C444         */
    {0.2,  4.0},  /*   Cassiopeia A  */
    {0.5,  2.0}   /*   0408-65 Mauch/Hugo  */
  };

  /* Source list - aliases */
  gchar *knosou[] = {
    "J0133-3629",                                               /* 0 */
    "3C48",       "0134+329", "0137+331", "J0137+3309", "3c48", /* 1 */
    "Fornax A",                                                 /* 2 */
    "3C123",      "0433+295", "0437+296", "J0437+2940", "3c123",/* 3 */
    "J0444+2809",                                               /* 4 */
    "3C138",      "0518+165", "0521+166", "J0521+1638", "3c138",/* 5 */
    "Pictor A",                                                 /* 6 */
    "Taurus A",                                                 /* 7 */
    "3C147",      "0538+498", "0542+498", "J0542+4951", "3c147",/* 8 */
    "3C196",      "0809+483", "0813+482", "J0813+4822", "3c196",/* 9 */
    "Hydra A",                                                  /* 10 */
    "Virgo A",    "M87", "3C274", "3c274",                      /* 11 */
    "3C286",      "1328+307", "1331+305", "J1331+3030", "3c286", "3C 286",/* 12 */
    "3C295",      "1409+524", "1411+522", "J1411+5212", "3c295",/* 13 */
    "Hercules A",                                               /* 14 */
    "3C353  ",    "3c353",                                      /* 15 */
    "3C380",      "3c380",                                      /* 16 */
    "Cygnus A",                                                 /* 17 */
    "3C444",      "3c444",                                      /* 18 */
    "Cassiopeia A",                                             /* 19 */ 
    "0408-65",   "0408-658", "J0407-658",                      /* 20 */ 
    NULL
  };
  
  /* Number of characters to check */
  olong lenchk[] = {
  /* 0      1       2        3       4      5      */
    10,  4,8,8,10,4, 8,  5,8,8,10,5, 10,  5,8,8,10,5,
  /* 6  7     8           9       10    11           12  */
    8,  8, 5,8,8,10,5, 5,8,8,10,5, 8,  8,3,5,5,  5,8,8,10,5,5, 
  /*   13       14  15    16   17 18   19     20*/
    5,8,8,10,5, 10, 5,5, 5,5,  8, 5,5,  12,  7,8,9
  };
  
  /* Source number 0-rel in coef table */
  olong ksouno[] = {
    0,  1,1,1,1,1,  2,  3,3,3,3,3, 4,  5,5,5,5,5,
    6,  7,  8,8,8,8,8,  9,9,9,9,9, 10, 11,11,11,11, 12,12,12,12,12,12,
    13,13,13,13,13,  14, 15,15,  16,16,  17,  18,18, 19, 20,20,20
  };
  
  olong j;
  gchar tempName[101];
  
  /* Name */
  for (j=0; j<16; j++) {tempName[j] = row->Source[j];} tempName[j] = 0;
  if (strlen(tempName)<16) {
    for (j=strlen(tempName); j<16; j++) tempName[j] = ' '; 
    tempName[j] = 0;
  }
  
  freqg = Freq*1.0e-9;  /* Frequency in GHz */
  
  /* Lookup source NULL terminated list */
  isrc = -1; i = 0;
  while (knosou[i]) {
    if (!strncmp(tempName, knosou[i], lenchk[i])) isrc = ksouno[i];
    i++;
  }
  
  /* Complain if not found */
  if (isrc<0) {
    Obit_log_error(err, OBIT_Error, "Unknown Flux calibrator %s", tempName);
    return;
  }
  /* Warn if out of range */
  if ((freqg<range[isrc][0]) || (freqg>range[isrc][1]))
     Obit_log_error(err, OBIT_InfoWarn, "Frequency %lf out of range [%f,%f] GHz for %s", 
		    freqg, range[isrc][0], range[isrc][1], tempName);

  /*  Compute flux */
  dt = log10 (freqg);
  temp2 = coeff[isrc][0] + dt * (coeff[isrc][1] + dt * (coeff[isrc][2] + 
	    dt * (coeff[isrc][3] + dt * (coeff[isrc][4] + dt * (coeff[isrc][5])))));
  Flux[0] = pow (10.0, temp2);

  /* Error - sum of squares */
  ferror =  (coerr[isrc][0]*coerr[isrc][0]) +
    ((coerr[isrc][1]*dt)*(coerr[isrc][1]*dt)) +
    ((coerr[isrc][2]*dt*dt)*(coerr[isrc][2]*dt*dt)) +
    ((coerr[isrc][3]*dt*dt*dt)*(coerr[isrc][3]*dt*dt*dt)) +
    ((coerr[isrc][4]*dt*dt*dt*dt)*(coerr[isrc][4]*dt*dt*dt*dt)) +
    ((coerr[isrc][5]*dt*dt*dt*dt*dt)*(coerr[isrc][5]*dt*dt*dt*dt*dt));
  /*  If non-zero error exponent */
  if (ferror>0.0) {
    /* sqrt sum of squares * factor */
    ferror = log(10.) * Flux[0] * sqrt(ferror);
  }
  
  Flux[0] *= Parms[2];  /* Apply fudge factor */
} /* end CalcFluxPB17 */
