/* $Id$  */
/* Obit Radio interferometry calibration software                     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2025                                          */
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
#include "ObitUV.h"
#include "ObitTableSN.h"
#include "ObitTableSU.h"
#include "ObitTableSUUtil.h"
#include "ObitSpectrumFit.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* GetJyIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void GetJyOut (ObitInfoList* outList, ObitErr *err);
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
void GetJyHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Update AIPS SN/SU tables */
ObitTableSN* GetJyUpdate (ObitUV* inData, ObitErr* err);
/* Swallow selected SN table entries */
void  ReadSN (ObitTableSN* SNTab, ObitUV *inData, olong CalSrc,
	      olong maxIF, olong maxAnt, olong maxSou, ofloat *oldFlux,
	      gchar *Sources, olong lsou, olong nsou, gchar *souCode, 
	      olong subA, olong Qual, olong FreqID, olong BIF, olong EIF, 
	      ofloat timeRange[2], olong *Antennas, olong nant, 
	      ofloat *offCnt, ofloat *offCnt2, ofloat *offSum, ofloat *offSum2, 
	      olong *offCalSrc, ObitErr* err);
/* Determine flux densities */
void DetFlux (olong maxIF, olong maxAnt, olong maxSou, ofloat *offCnt, ofloat *offCnt2, 
	      ofloat *offSum, ofloat *offSum2, olong *offCalSou, ofloat *oldFlux, 
	      ofloat *souFlux, ofloat *souErr, ObitErr* err);
/* Fit calibrator spectra */
void FitSpec (ObitUV* inData, gchar **source, olong maxIF, olong maxSou, olong nterm, 
	      ofloat *souFlux, ofloat *souErr, 
	      ofloat *newFlux, ObitErr* err);
/* Update SN table */
void UpdateSN (ObitTableSN* SNTab, olong maxIF, olong maxAnt, olong maxSou, 
	       olong BIF, olong EIF,
	       olong *offCalSou, ofloat *oldFlux, ofloat *souFlux,ObitErr* err);
/* Determine alpha median average  */
static ofloat MedianAvg (ollong n, ofloat *value, ofloat alpha);
/* Determine straight average  */
static ofloat StraightAvg (ollong n, ofloat *value);

/* Determine sigma for Median */
static ofloat MedianSigma (ollong n, ofloat *value, ofloat mean, ofloat alpha);
/* Weighted linear regression */
static void LinRegWt (olong n, odouble *nu, ofloat *S, ofloat *sigma,
		      ofloat *S_0, ofloat* alpha);
 
/* Program globals */
gchar *pgmName = "GetJy";       /* Program name */
gchar *infile  = "GetJy.in" ;   /* File with program inputs */
gchar *outfile = "GetJy.out";   /* File to contain program outputs */
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
/*   Obit task to determine calibrator flux densities and update cal.     */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitErr      *err= NULL;
  ObitTableSN  *SNTable = NULL;
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "souCode", "Qual", "calSour", "calCode",
    "FreqID", "timeRange",  "subA", "Antennas", "BIF", "EIF",
    "solnVer", 
    NULL};

  /* Startup - parse command line */
  err = newObitErr();
  myInput = GetJyIn (argc, argv, err);
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

  /* Solve for flux densities and update SN/SU table */
  SNTable =  GetJyUpdate (inData, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Write history */
  GetJyHistory (myInput, inData, err); 
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* show any messages and errors */
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* cleanup */
  SNTable   = ObitTableSNUnref(SNTable);
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput  = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* GetJyIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "GetJyIn";

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
} /* end GetJyIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: GetJy -input file -output ofile [args]\n");
    fprintf(stderr, "Obit task to setewrmine flux densities and apply \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def GetJy.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def GetJy.out\n");
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
  strTemp = "GetJy.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "GetJyName";
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
    "Sources", "timeRange",  "BIF", "EIF", "subA", "Antennas",
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
/*  Write History for GetJy                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void GetJyHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "souCode", "Qual", "calSour", "calCode", "BIF", "EIF",
    "FreqID", "timeRange",  "subA", "Antennas",  "solnVer", "nTerm",
    NULL};
  gchar *routine = "GetJyHistory";

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
 
} /* end GetJyHistory  */

/*----------------------------------------------------------------------- */
/*  Determine Source flux densities, update SU table, apply to SN         */
/*   Input:                                                               */
/*      inData    ObitUV with SU/SN tables to update                      */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Returns AIPS SU table                                                */
/*----------------------------------------------------------------------- */
ObitTableSN* GetJyUpdate (ObitUV* inData, ObitErr* err)
{
  ObitTableSU *SUTab=NULL;
  ObitTableSN *SNTab=NULL;
  ObitTableSURow *row=NULL;
  ofloat *offCnt=NULL, *offCnt2=NULL, *offSum=NULL, *offSum2=NULL;
  ofloat *oldFlux=NULL, *souFlux=NULL, *newFlux=NULL, *souErr=NULL;
  olong *offCalSou=NULL;
  olong ver, offset, nterm=2;
  oint numIF=1, numPol=1;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong i, j, isou, irow, lsou, nsou, lcal, ncal,iif,  bif, eif;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *Sources=NULL, *souCode=NULL, *calSour=NULL, *calCode=NULL;
  olong *Antennas=NULL;
  ofloat timeRange[2];
  olong Qual, maxAnt, maxSou=0, maxIF, nant, FreqID=0, BIF=1, EIF=0, subA=0, solnVer=0;
  gchar tempName[101]; /* should always be big enough */
  gchar *blank = "        ", **source=NULL;
  gchar *routine = "GetJyUpdate";

  /* error checks */
  if (err->error) return SNTab;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Get control parameters */
  ObitInfoListGetP(inData->info, "Sources",  &type, dim, (gpointer)&Sources);
  /* Count number of actual sources */
  lsou = dim[0];
  nsou = 0;
  for (i=0; i<dim[1]; i++) {
    if (Sources[i*lsou]!=' ') nsou = i+1;
    else break;
  }
  ObitInfoListGetP(inData->info, "calSour",  &type, dim, (gpointer)&calSour);
  /* Count number of actual calibrators */
  lcal = dim[0];
  ncal = 0;
  for (i=0; i<dim[1]; i++) {
    if (calSour[i*lcal]!=' ') ncal = i+1;
    else break;
  }
  ObitInfoListGetP(inData->info, "souCode",  &type, dim, (gpointer)&souCode);
  ObitInfoListGetP(inData->info, "calCode",  &type, dim, (gpointer)&calCode);
  ObitInfoListGetP(inData->info, "Antennas",  &type, dim,(gpointer) &Antennas);
  nant = dim[1];
  ObitInfoListGetTest(inData->info, "timeRange",  &type, dim, timeRange);
  if (timeRange[1]<=0.0) timeRange[1] = 1.0e20;
  ObitInfoListGetTest(inData->info, "Qual",       &type, dim, &Qual);
  ObitInfoListGetTest(inData->info, "BIF",        &type, dim, &BIF);
  ObitInfoListGetTest(inData->info, "EIF",        &type, dim, &EIF);
  ObitInfoListGetTest(inData->info, "subA",       &type, dim, &subA);
  ObitInfoListGetTest(inData->info, "FreqID",     &type, dim, &FreqID);
  ObitInfoListGetTest(inData->info, "solnVer",    &type, dim, &solnVer);
  ObitInfoListGetTest(inData->info, "nTerm",      &type, dim, &nterm);
  if (calCode==NULL) calCode = blank;
  if (souCode==NULL) souCode = blank;

  /* Get AIPS SN table */
  ver = solnVer;
  numPol = 0;
  SNTab = newObitTableSNValue (inData->name, (ObitData*)inData, &ver, 
			       OBIT_IO_ReadWrite, numPol, numIF, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, SNTab);
 
  /* Get AIPS SU table */
  ver = 1;
  SUTab = newObitTableSUValue (inData->name, (ObitData*)inData, &ver, 
			       OBIT_IO_ReadWrite, numIF, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, SNTab);
 
  /* Open SU table */
  retCode = ObitTableSUOpen (SUTab, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Create table row */
  row = newObitTableSURow (SUTab);
  ObitTableSUSetRow (SUTab, row, err);
  if (err->error) goto cleanup;

  /* Lookup maximum source ID, get previous fluxes */
  maxSou = 0;
  irow = 0;
  retCode = OBIT_IO_OK;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableSUReadRow (SUTab, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    maxSou = MAX (maxSou, row->SourID);
  } /* end loop looking up max source ID */

  /* IF range */
  numIF = SUTab->numIF;
  bif = MAX (1, BIF);
  eif = MIN (EIF, numIF);
  if (eif<=0) eif = numIF;

  /* Create accumulation arrays big enough to swallow everything 
   data in order IF, ant, source */

  maxAnt = inData->myDesc->maxAnt;
  maxIF  = numIF;
  offCnt  = g_malloc0(maxAnt*maxIF*maxSou*sizeof(ofloat));
  offCnt2 = g_malloc0(maxAnt*maxIF*maxSou*sizeof(ofloat));
  offSum  = g_malloc0(maxAnt*maxIF*maxSou*sizeof(ofloat));
  offSum2 = g_malloc0(maxAnt*maxIF*maxSou*sizeof(ofloat));
  oldFlux = g_malloc0(maxIF*maxSou*sizeof(ofloat));
  souFlux = g_malloc0(maxIF*maxSou*sizeof(ofloat));
  newFlux = g_malloc0(maxIF*maxSou*sizeof(ofloat));
  souErr  = g_malloc0(maxIF*maxSou*sizeof(ofloat));
  offCalSou= g_malloc0(maxSou*sizeof(olong));
  source   = g_malloc0(maxSou*sizeof(gchar*));

  /* Get previous fluxes, source names */
  irow = 0;
  retCode = OBIT_IO_OK;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableSUReadRow (SUTab, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    /* Save old flux */
    for (iif=bif-1; iif<eif; iif++) {
      offset = MAX (0, row->SourID-1) * maxIF;
      /* Default flux density = 1.0 (as used in Calib) */
      if (row->IFlux[iif]==0.0)  oldFlux[offset+iif] = 1.0;
      else oldFlux[offset+iif] = row->IFlux[iif];
    }
    offset = MAX (0, row->SourID-1);
    source[offset] = g_strdup(row->Source);
    source[offset][15] = 0;  /* Terminate */
  } /* end loop looking up max source ID */

  /* Close SU table - avoid conflict in ReadSN */
  retCode = ObitTableSUClose (SUTab, err);

  /* Swallow selected SN table entries */
  /* Source list */
  ReadSN (SNTab, inData, 1,
	  maxIF, maxAnt, maxSou, oldFlux, Sources, lsou, nsou, 
	  souCode, subA, Qual, FreqID, bif, eif, timeRange, 
	  Antennas, nant, 
	  offCnt, offCnt2, offSum, offSum2, offCalSou, err);
  if (err->error) goto cleanup;

  /* Calibrator list */
  ReadSN (SNTab, inData, 2,
	  maxIF, maxAnt, maxSou, oldFlux, calSour, lcal, ncal, 
	  calCode, subA, Qual, FreqID, bif, eif, timeRange, 
	  Antennas, nant, 
	  offCnt, offCnt2, offSum, offSum2, offCalSou, err);
  if (err->error) goto cleanup;

  /* Determine flux densities */
  DetFlux (maxIF, maxAnt, maxSou, offCnt, offCnt2, offSum, offSum2, 
	   offCalSou, oldFlux, souFlux, souErr, err);
  if (err->error) goto cleanup;

  /* Fit calibrator spectra */
  FitSpec (inData, source, maxIF, maxSou, nterm, souFlux, souErr, newFlux, err);
  if (err->error) goto cleanup;

  /* Open SU table */
  retCode = ObitTableSUOpen (SUTab, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Update SU table */
  irow = 0;
  retCode = OBIT_IO_OK;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableSUReadRow (SUTab, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    isou =  MAX (0, row->SourID-1);
    if (offCalSou[isou]!=1) continue;  /* Only fitted sources */

    /* Update flux density  */
    for (iif=bif-1; iif<eif; iif++) {
      offset = MAX (0, row->SourID-1) * maxIF;
      row->IFlux[iif] = newFlux[offset+iif];

      /* Message */
      /* Get source name */
      for (j=0; j<16; j++) {tempName[j] = ' ';} tempName[j] = 0;
      /* get blank padded name */
      for (j=0; j<lsou; j++) {
	if (row->Source[j]==0) break;  /* only values in string */
	tempName[j] = row->Source[j]; 
      }
      if (souFlux[offset+iif]>0.0) {
	Obit_log_error(err, OBIT_InfoErr,
		       "%s IF %d Flux %8.3f +/- %8.3f Jy, Smoothed %8.3f Jy",
		       tempName, iif+1, souFlux[offset+iif], souErr[offset+iif], 
		       newFlux[offset+iif]);
      } else { /* No solution */
	Obit_log_error(err, OBIT_InfoWarn,
		       "%s IF %d No Flux density measured,      Smoothed %8.3f Jy", 
		       tempName, iif+1, newFlux[offset+iif]);
      }
    } /* end IF loop */

    /* Update SU row */
    retCode = ObitTableSUWriteRow (SUTab, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  } /* end loop updating SU table */

  /* Update SN table */
  UpdateSN (SNTab, maxIF, maxAnt, maxSou, bif, eif, offCalSou, 
	    oldFlux, newFlux, err);
  if (err->error) goto cleanup;

  /* Cleanup */
 cleanup:
  /* Close SU table */
  retCode = ObitTableSUClose (SUTab, err);
  row = ObitTableSURowUnref (row);
  SUTab = ObitTableSUUnref (SUTab);
  if (offCnt)    g_free(offCnt);
  if (offCnt2)   g_free(offCnt2);
  if (offSum)    g_free(offSum);
  if (offSum2)   g_free(offSum2);
  if (oldFlux)   g_free(oldFlux);
  if (souFlux)   g_free(souFlux);
  if (newFlux)   g_free(newFlux);
  if (souErr)    g_free(souErr);
  if (offCalSou) g_free(offCalSou);
  if (source) {
    for (isou=0; isou<maxSou; isou++) {
      if (source[isou]) g_free (source[isou]);
    }
    g_free (source);
  }
  if (err->error) Obit_traceback_val (err, routine, inData->name, SNTab);
  return SNTab;
} /* end GetJyUpdate  */

/*----------------------------------------------------------------------- */
/*  Accumulate selected amplitude entries in SN table                     */
/*  Uses MedianAvg w/ alpha=0.5 for robust averaging                      */
/*   Input:                                                               */
/*     SNTab     SNTable to process                                       */
/*     inData    UV data with which SNTab associated, will modify selector*/
/*     CalSou    1 or 2 indicating source or calibrator                   */
/*     maxIF     Maximum IF number in accumulators                        */
/*     maxAnt    Maximum antenna number in accumulators                   */
/*     maxSou    Maximum source ID in accumulators                        */
/*     Sources   Selected source names, [0] blank=> any                   */
/*     lsou      length of source name in Sources                         */
/*     nsou      maximum number of entries in Sources                     */
/*     souCode   selection of Source by Calcode,if not specified in Source*/
/*                  '    ' => any calibrator code selected                */
/*                  '*   ' => any non blank code (cal. only)              */
/*                  '-CAL' => blank codes only (no calibrators)           */
/*                  anything else = calibrator code to select.            */
/*               NB: The souCode test is applied in addition to the       */
/*               other tests, i.e. Sources and Qual, in the               */
/*               selection of sources to process.                         */
/*     subA      Subarray 0-> any                                         */
/*     Qual      Source qualifier, -1 => any                              */
/*     FreqID    Selected Frequency ID                                    */
/*     BIF       First IF to process                                      */
/*     EIF       Highest IF to process                                    */
/*     timeRange time range (days) to consider, 0's -> all                */
/*     Antennas  List of antennas to include. If any number is negative   */
/*               then all antennas listed  are NOT desired and all others */
/*               are. All 0 => use all.                                   */
/*     nant      maximum number of entries in Antennas                    */
/*   Output:                                                              */
/*     offCnt   Counts of samples (source, ant, IF) Pol 1.                */
/*     offCnt2  Counts of samples (source, ant, IF) Pol 2                 */
/*     offSum   Sums of gains (source,ant,IF) Pol 1.                      */
/*     offSum2  Sums of gains (source,ant,IF) Pol 2                       */
/*     offCalSou array of 1 (source) or 2 (Calibrator) per SourceID       */
/*      err     Obit Error stack                                          */
/*----------------------------------------------------------------------- */
void  ReadSN (ObitTableSN* SNTab, ObitUV *inData, olong CalSou,
	      olong maxIF, olong maxAnt, olong maxSou, ofloat *oldFlux,
	      gchar *Sources, olong lsou, olong nsou, gchar *souCode, 
	      olong subA, olong Qual, olong FreqID, olong BIF, olong EIF, 
	      ofloat timeRange[2], olong *Antennas, olong nant, 
	      ofloat *offCnt, ofloat *offCnt2, ofloat *offSum, ofloat *offSum2, 
	      olong *offCalSou, ObitErr* err)
{
  ObitIOCode retCode;
  ObitTableSNRow *SNrow=NULL;
  ObitUVSel *sel;
  olong iif, irow, offset, offset2, numPol, ii, jj, isou, iant;
  ollong nmedn, i, size;
  ofloat **store1=NULL, **store2=NULL;
  ofloat alpha=0.5; /* median alpha */
  olong maxCells = 1000;  /* Maximum number of time cells */
  gboolean want, exceeded=FALSE;
  ofloat amp, fblank = ObitMagicF();
  gchar *routine = "ReadSN";
 
  /* error checks */
  if (err->error) return;

  /* Setup selector */
  sel  = inData->mySel; 
  ObitUVSelSetAnt (sel, Antennas, nant);
  ObitUVSelSetSour (sel, (ObitData*)inData, Qual, souCode, Sources, lsou, nsou, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);

  /* Open table */
  retCode = ObitTableSNOpen (SNTab, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Create table row */
  SNrow = newObitTableSNRow (SNTab);
  ObitTableSNSetRow (SNTab, SNrow, err);
  numPol  = SNTab->numPol;  /* Number of polarizations */

  /* Create working storage arrays, maxCells cells of maxSou*maxIF*maxAnt */
  size = maxSou*maxIF*maxAnt;
  store1 = g_malloc0(size*sizeof(ofloat*)); 
  if (numPol>1) store2 = g_malloc0(size*sizeof(ofloat*));
  for (i=0; i<size; i++) {
    store1[i] = g_malloc0(maxCells*sizeof(ofloat));
    if (numPol>1) store2[i] = g_malloc0(maxCells*sizeof(ofloat));
  }

  /* Loop over table */
  for (irow=1; irow<=SNTab->myDesc->nrow; irow++) {
    retCode = ObitTableSNReadRow (SNTab, irow, SNrow, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    if (SNrow->status < 0) continue; /* entry flagged? */

    /* Want this one?  check subarray */
    want = ((SNrow->SubA == subA)  ||  (SNrow->SubA <= 0) || (subA <= 0));
    
    /* check frqsel */
    want = want &&
      ((SNrow->FreqID == FreqID) || (SNrow->FreqID <= 0) || (FreqID <= 0));

    /* Want source? */
    want = want &&  ObitUVSelWantSour (sel, SNrow->SourID);

    /* Want Antenna? */
    want = want &&  ObitUVSelWantAnt (sel, SNrow->antNo);

    /* Want time */
    want = want &&  (SNrow->Time>=timeRange[0]) && (SNrow->Time<=timeRange[1]);

    /* skip if not wanted */
    if (!want) continue;

    /* Save src/cal flag */
    offCalSou[SNrow->SourID-1] = CalSou;

    /* Offset in accumulators */
    offset = MAX (0, SNrow->SourID-1) * maxAnt * maxIF +
      MAX (0, SNrow->antNo-1) * maxIF;
    offset2 = MAX (0, SNrow->SourID-1) * maxIF;
    for (iif=BIF-1; iif<EIF; iif++) {
      if ((SNrow->Weight1[iif]>0.0) && (SNrow->Real1[iif]!=fblank)) {
	if (offCnt[offset+iif]<maxCells) { /* Is there room? */
	  amp = 1.0 / (SNrow->Real1[iif]*SNrow->Real1[iif] + SNrow->Imag1[iif]*SNrow->Imag1[iif]);
	  if (oldFlux[offset2+iif]>0.0) amp *= oldFlux[offset2+iif];
	  ii = offset+iif;         /* source, IF, ant index */
	  jj = offCnt[offset+iif]; /* time cell index */
	  offCnt[offset+iif]++;
	  store1[ii][jj] = amp;
	} else exceeded = TRUE;  /* exceeded arrays */
      }
      /* Second poln? */
      if ((numPol>1) && (SNrow->Weight2[iif]>0.0) && (SNrow->Real2[iif]!=fblank)) {
	if (offCnt2[offset+iif]<maxCells) { /* Is there room? */
	  amp = 1.0 / (SNrow->Real2[iif]*SNrow->Real2[iif] + SNrow->Imag2[iif]*SNrow->Imag2[iif]);
	  if (oldFlux[offset2+iif]>0.0) amp *= oldFlux[offset2+iif];
	  ii = offset+iif;          /* source, IF, ant index */
	  jj = offCnt2[offset+iif]; /* time cell index */
	  offCnt2[offset+iif]++;
	  store2[ii][jj] = amp;
	} else exceeded = TRUE;  /* exceeded arrays */
      }
    } /* end IF loop */
   
  } /* end loop over table */

  /* Median alpha averaging */
  for (isou=0; isou<maxSou; isou++) {
    if (!ObitUVSelWantSour (sel, isou+1)) continue;    
    for (iant=0; iant<maxAnt; iant++) {
      for (iif=BIF-1; iif<EIF; iif++) {
	offset = isou*maxAnt*maxIF + iant*maxIF;
	if (offCnt[offset+iif]>3) {
	  nmedn = offCnt[offset+iif];
	  offSum[offset+iif] = MedianAvg(nmedn, store1[offset+iif], alpha);
	  offCnt[offset+iif] = 1;
	} else {
	  nmedn = offCnt[offset+iif];
	  offSum[offset+iif] = StraightAvg(nmedn, store1[offset+iif]);
	  offCnt[offset+iif] = MIN(1,offCnt[offset+iif]);
	}
	if (numPol>1) {
	  if (offCnt2[offset+iif]>3) {
	    nmedn = offCnt2[offset+iif];
	    offSum2[offset+iif] = MedianAvg(nmedn, store2[offset+iif], alpha);
	    offCnt2[offset+iif] = 1;
	  } else  {
	    nmedn = offCnt2[offset+iif];
	    offSum2[offset+iif] = StraightAvg(nmedn, store2[offset+iif]);
	    offCnt2[offset+iif] =  MIN(1,offCnt2[offset+iif]);
	  }
	}
      } /* end IF loop */
    } /* end antenna loop */
  } /* end source loop */

  /* Did I blow core? */
  if (exceeded) {
	Obit_log_error(err, OBIT_InfoErr,
		       "%s Exceeded maximum number of time samples", routine);
  }

 cleanup:
  SNrow = ObitTableSNRowUnref(SNrow);
  retCode = ObitTableSNClose (SNTab, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
  /* Free working storage arrays*/
  size = maxSou*maxIF*maxAnt;
  for (i=0; i<size; i++) {
    if (store1[i]) g_free(store1[i]);
    if (store2 && store2[i]) g_free(store2[i]);
  }
  if (store1) g_free(store1);
  if (store2) g_free(store2);

} /* end ReadSN */

/*----------------------------------------------------------------------- */
/*  Determine flux densities and errors                                   */
/* Adapted from AIPSish GetJY.FOR/GJYFLX                                  */
/*  Uses MedianAvg w/ alpha=0.5 for robust averaging                      */
/*   Input:                                                               */
/*     maxIF     Maximum IF number in accumulators                        */
/*     maxAnt    Maximum antenna number in accumulators                   */
/*     maxSou    Maximum source ID in accumulators                        */
/*     offCnt    Counts of samples (source, ant, IF) Pol 1.               */
/*     offCnt2   Counts of samples (source, ant, IF) Pol 2                */
/*     offSum    Sums of gains (source,ant,IF) Pol 1.                     */
/*     offSum2   Sums of gains (source,ant,IF) Pol 2                      */
/*     offCalSou array of 1 (source) or 2 (Calibrator) per SourceID       */
/*     oldFlux   Flux density in SU table (source, IF)                    */
/*   Output:                                                              */
/*     souFlux   Derived flux density  (source, IF)                       */
/*     souErr    Error in derived flux density  (source, IF)              */
/*      err      Obit Error stack                                         */
/*----------------------------------------------------------------------- */
void DetFlux (olong maxIF, olong maxAnt, olong maxSou, ofloat *offCnt, ofloat *offCnt2, 
	      ofloat *offSum, ofloat *offSum2, olong *offCalSou, ofloat *oldFlux, 
	      ofloat *souFlux, ofloat *souErr, ObitErr* err)
{
  olong isou, iant, iif;
  olong offset, offset2, offset3;
  ofloat *rcal1=NULL, *rcal2=NULL, *ical1=NULL, *ical2=NULL;
  ofloat fluxd, rms;
  ollong count, size;
  ofloat *store=NULL;
  ofloat alpha=0.5; /* median alpha */
  gchar *routine = "DetFlux";

  /* error checks */
  if (err->error) return;

  /* Local accumulators */
  rcal1 = g_malloc0(maxAnt*maxIF*maxSou*sizeof(ofloat));
  ical1 = g_malloc0(maxAnt*maxIF*maxSou*sizeof(ofloat));
  rcal2 = g_malloc0(maxAnt*maxIF*maxSou*sizeof(ofloat));
  ical2 = g_malloc0(maxAnt*maxIF*maxSou*sizeof(ofloat));

  /* Loop over calibrators, compute receiver cal factors in Jy/whatever */
  count = 0;
  for (isou=0; isou<maxSou; isou++) {
    if (offCalSou[isou]!=2) continue;
    for (iif=0; iif<maxIF; iif++) {
      for (iant=0; iant<maxAnt; iant++) {
	offset  = isou * maxAnt * maxIF + iant*maxIF;
	offset2 = isou * maxIF;
	offset3 = iant * maxIF;
	/* First polarization */
	if ((offCnt[offset+iif]>0) && (oldFlux[offset2+iif]>1.0e-20)) {
	  rcal1[offset3+iif] += oldFlux[offset2+iif] * offCnt[offset+iif] / offSum[offset+iif];
	  ical1[offset3+iif]++;
	  count++;
	}
	/* Second polarization */
	if ((offCnt2[offset+iif]>0) && (oldFlux[offset2+iif]>1.0e-20)) {
	  rcal2[offset3+iif] += oldFlux[offset2+iif] * offCnt2[offset+iif] / offSum2[offset+iif];
	  ical2[offset3+iif]++;
	  count++;
	}
      } /* end loop over antenna */
    } /* end loop over IF */
  } /* end loop over source */

  /* Check that got some */
  if (count<=0) {
    Obit_log_error(err, OBIT_Error, "NO Valid Calibration data found");
    goto cleanup;
  }

  /* Average receiver calibration factors - 
     ical should be 0 or 1 here so not really necessary any more */
  for (iif=0; iif<maxIF; iif++) {
    for (iant=0; iant<maxAnt; iant++) {
      offset2 = iant*maxIF;
      if (ical1[offset2+iif]>0.0) rcal1[offset2+iif] /= ical1[offset2+iif];
      if (ical2[offset2+iif]>0.0) rcal2[offset2+iif] /= ical2[offset2+iif];
    } /* end loop over antenna */
  } /* end loop over IF */
 
  /* Storage arrays for robust average/RMS */
  size = maxAnt*2;
  store = g_malloc0(size*sizeof(ofloat)); 

  /* Loop over sources - compute source fluxes and RMSes */
  for (isou=0; isou<maxSou; isou++) {
    if (offCalSou[isou]!=1) continue;
    for (iif=0; iif<maxIF; iif++) {
      count = 0;
      offset2 = isou*maxIF;
      for (iant=0; iant<maxAnt; iant++) {
	offset  = isou*maxAnt*maxIF + iant*maxIF;
	offset3 = iant*maxIF;
	/* First poln */
	if ((offCnt[offset+iif]>0) && (ical1[offset3+iif]>0.0)) {
	  fluxd = rcal1[offset3+iif] * (offSum[offset+iif] / offCnt[offset+iif]);
	  store[count] = fluxd;
	  count++;
	}
	/* Second poln */
	if ((offCnt2[offset+iif]>0) && (ical2[offset3+iif]>0.0)) {
	  fluxd = rcal2[offset3+iif] * (offSum2[offset+iif] / offCnt2[offset+iif]);
	  store[count] = fluxd;
	  count++;
	}
      } /* end loop over antenna */

      /* Compute  source average, RMS */
      fluxd = 0.0;
      rms   = -1.0;
      if (count>3) fluxd = MedianAvg(count, store, alpha);
      else fluxd = StraightAvg(count, store);
      if (count>5) rms   = MedianSigma(count, store, fluxd, alpha);
      /* Save values */
      souFlux[offset2+iif] = fluxd;
      souErr[offset2+iif]  = rms;
    } /* end loop over IF */
  } /* end loop over source */

  /* Cleanup */
 cleanup:
  if (store) g_free(store);
  if (rcal1) g_free(rcal1);
  if (ical1) g_free(ical1);
  if (rcal2) g_free(rcal2);
  if (ical2) g_free(ical2);
  if (err->error) Obit_traceback_msg (err, routine, routine);
} /* end DetFlux */

/*----------------------------------------------------------------------- */
/*  Least Squares fit to calibrator spectra                               */
/*  Spectra are fitted to each calibrator and the IF values in newFlux    */
/*  are those determined from the fit.                                    */
/*  If nterm<1 then no fitting is done and souFlux values are copied      */
/*  to newFlux                                                            */
/*   Input:                                                               */
/*     inData    ObitUV with SU/SN tables to update                       */
/*     source    Array of source names                                    */
/*     maxIF     Maximum IF number in accumulators                        */
/*     maxSou    Maximum source ID in accumulators                        */
/*     ntern     Number of terms in the spectral fit                      */
/*     souFlux   Derived flux density  (source, IF)                       */
/*     souErr    Error in derived flux density  (source, IF)              */
/*   Output:                                                              */
/*     newFlux   Derived flux density  (source, IF)                       */
/*      err      Obit Error stack                                         */
/*----------------------------------------------------------------------- */
void FitSpec (ObitUV* inData, gchar **source, 
	      olong maxIF, olong maxSou, olong nterm, 
	      ofloat *souFlux, ofloat *souErr, float *newFlux, 
	      ObitErr* err)
{
  olong isou, iif, iterm, offset2, count;
  ofloat *flux=NULL, *sigma=NULL, fitVal[5]={0.,0.,0.,0.,0.};
  odouble *freq=NULL, lnu, lll, arg;

  /* error checks */
  if (err->error) return;
  if (maxIF<=0) return;
  if (maxSou<=0) return;
  /* Copy souFlux to newFlux in case */
  for (isou=0; isou<maxSou; isou++) {
    for (iif=0; iif<maxIF; iif++) {
      offset2 = isou*maxIF;
      newFlux[offset2+iif] = souFlux[offset2+iif];
    }
  }
  /* Anything to do? */
  if (nterm<=0) return;

  /* Work arrays */
  flux  = g_malloc0(maxIF*sizeof(ofloat));
  sigma = g_malloc0(maxIF*sizeof(ofloat));
  freq  = g_malloc0(maxIF*sizeof(odouble));

  /* Loop over sources */
  for (isou=0; isou<maxSou; isou++) {
    /* Collect data to fit */
    count = 0;
    for (iif=0; iif<maxIF; iif++) {
      offset2 = isou*maxIF;
      if (souFlux[offset2+iif]>0.0) {
	freq[count]  = inData->myDesc->freqIF[iif];
	flux[count]  = souFlux[offset2+iif];
	sigma[count] = souErr[offset2+iif];
	/* Weight by 1/std dev. rather than 1/variance 
	sigma[count] = 0.1*sqrtf(souErr[offset2+iif]);*/
	//DEBUGif (count==0) sigma[count] *=0.01;
	count++;
      }
    }

    /* Anything? */
    if (count<1) continue;

    nterm = MIN (nterm, count);
    Obit_log_error(err, OBIT_InfoErr, "Smoothing by fitting %d term spectrum", nterm);

    /* Fit 
    fitVal = ObitSpectrumFitSingle(count, nterm, freq[0], freq, flux, sigma, FALSE, FALSE, NULL, err);*/

    /* Tell results 
    if (nterm==1) {
      Obit_log_error(err, OBIT_InfoErr,
		     "Source %s Flux density=%8.3f(%8.3f) Jy Chi^2= %8.3f", 
		     source[isou], fitVal[0], fitVal[nterm], fitVal[nterm*2]);
    }

    if (nterm>=2) {
      Obit_log_error(err, OBIT_InfoErr,
		     "Source %s Flux density=%8.3f(%5.3f) Jy, SI=%5.3f(%5.3f) Jy Chi^2=%5.1f", 
		     source[isou], fitVal[0], fitVal[nterm], fitVal[1], 
		     fitVal[nterm+1], fitVal[nterm*2]);
    }*/

    /* Try linear regression */
    LinRegWt (count, freq, flux, sigma, &fitVal[0], &fitVal[1]);
    Obit_log_error(err, OBIT_InfoErr,
		   "Source %s Fitted Flux density=%8.3f Jy, SI=%5.3f", 
		   source[isou], fitVal[0], fitVal[1]);
    
    /* Calculate new */
    for (iif=0; iif<maxIF; iif++) {
      offset2 = isou*maxIF;
      lnu = log (inData->myDesc->freqIF[iif]/inData->myDesc->freqIF[0]);
      lll = lnu; arg = 0.0;
      for (iterm=1; iterm<nterm; iterm++) {
	arg += fitVal[iterm] * lll;
	lll *= lnu;
      }
      newFlux[offset2+iif] =  exp(arg) * fitVal[0];
    } /* end IF loop */
    /*if (fitVal) g_free(fitVal);*/
  } /* end of source loop */
  /* Cleanup */
  if (flux)  g_free(flux);
  if (sigma) g_free(sigma);
  if (freq)  g_free(freq);
} /* end FitSpec  */
/*----------------------------------------------------------------------- */
/*  Update SN table for new fluxes                                        */
/*   Input:                                                               */
/*     SNTab     SNTable to process                                       */
/*     maxIF     Maximum IF number in accumulators                        */
/*     maxAnt    Maximum antenna number in accumulators                   */
/*     maxSou    Maximum source ID in accumulators                        */
/*     offCalSou array of 1 (source) or 2 (Calibrator) per SourceID       */
/*     oldFlux   Flux density in SU table (source, IF)                    */
/*     souFlux   Derived flux density  (source, IF)                       */
/*   Output:                                                              */
/*      err     Obit Error stack                                          */
/*----------------------------------------------------------------------- */
void UpdateSN (ObitTableSN* SNTab, olong maxIF, olong maxAnt, olong maxSou, 
	       olong BIF, olong EIF, 
	       olong *offCalSou, ofloat *oldFlux, ofloat *souFlux, ObitErr* err)
{
  ObitIOCode retCode;
  ObitTableSNRow *SNrow=NULL;
  olong iif, isou, irow, offset, numPol;
  ofloat amp, fblank = ObitMagicF();
  gchar *routine = "UpdateSN";
 
  /* error checks */
  if (err->error) return;

  /* Open table */
  retCode = ObitTableSNOpen (SNTab, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Create table row */
  SNrow = newObitTableSNRow (SNTab);
  ObitTableSNSetRow (SNTab, SNrow, err);
  numPol  = SNTab->numPol;  /* Number of polarizations */

  /* Loop over table */
  for (irow=1; irow<=SNTab->myDesc->nrow; irow++) {
    retCode = ObitTableSNReadRow (SNTab, irow, SNrow, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    if (SNrow->status < 0) continue; /* entry flagged? */

    /* Is this a source to update? */
    isou = MAX (0, SNrow->SourID-1);
    if (offCalSou[isou]!=1) continue; /* Source to be updated? */

    /* Offset in accumulators */
    offset = isou * maxIF;
    for (iif=BIF-1; iif<EIF; iif++) {
      if ((oldFlux[offset+iif]<=0.0) || (souFlux[offset+iif]<=0.0)) continue;
      amp = sqrt(souFlux[offset+iif]/oldFlux[offset+iif]);
      if ((SNrow->Weight1[iif]>0.0) && (SNrow->Real1[iif]!=fblank)) {
	SNrow->Real1[iif] *= amp;
	SNrow->Imag1[iif] *= amp;
      }
      /* Second poln? */
      if ((numPol>1) && (SNrow->Weight2[iif]>0.0) && (SNrow->Real2[iif]!=fblank)) {
 	SNrow->Real2[iif] *= amp;
	SNrow->Imag2[iif] *= amp;
     }
    } /* end IF loop */

    /* Rewrite */   
    retCode = ObitTableSNWriteRow (SNTab, irow, SNrow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  } /* end loop over table */

 cleanup:
  SNrow = ObitTableSNRowUnref(SNrow);
  retCode = ObitTableSNClose (SNTab, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
} /* end UpdateSN */
/**
 * ofloat comparison of two arguments
 * \param arg1 first value to compare
 * \param arg2 second value to compare
 * \return negative if arg1 is less than arg2, zero if equal
 *  and positive if arg1 is greater than arg2.
 */
static int compare_ofloat  (const void* arg1,  const void* arg2)
{
  int out = 0;
  ofloat larg1, larg2;

  larg1 = *(ofloat*)arg1;
  larg2 = *(ofloat*)arg2;
  if (larg1<larg2)      out = -1;
  else if (larg1>larg2) out = 1;
  return out;
} /* end compare_ofloat */

/**
 * Determine alpha median/average of an ofloat array
 * Use center 1-alpha of points, excluding at least one point from each end
 * \param n       Number of points
 * \param value   Array of values, sorted on return
 * \param alpha   0 -> 1 = central fraction of ordered list to use
 *                0.5 is generally a good compromise 
 * \return alpha median average, fblank if cannot determine
 */
static ofloat MedianAvg (ollong n, ofloat *value, ofloat alpha)
{
  ofloat out=0.0;
  ofloat fblank = ObitMagicF();
  ofloat beta, sum;
  olong i, i1, i2, count;

  if (n<=0) return out;

  /* Sort to ascending order */
  qsort ((void*)value, n, sizeof(ofloat), compare_ofloat);

  out = value[n/2];

  beta = MAX (0.05, MIN (0.95, 1.0-alpha)) / 2.0; /*  Average around median factor */

  /* Average around the center */
  i1 = MAX (0, (n/2)-(olong)(beta*n+0.5));
  i2 = MIN (n, (n/2)+(olong)(beta*n+0.5));

  if (i2>i1) {
    sum = 0.0;
    count = 0;
    for (i=i1; i<i2; i++) {
      if (value[i]!=fblank) {
	sum += value[i];
	count++;
      }
    }
    if (count>0) out = sum / count;
  }
   
  return out;
} /* end MedianAvg */

/**
 * Determine arithmetic average of an ofloat array
 * \param n       Number of points
 * \param value   Array of values, sorted on return
 * \return alpha median average
 */
static ofloat StraightAvg (ollong n, ofloat *value)
{
  ofloat out=0.0;
  ofloat sum;
  olong i, count;

  if (n<=0) return out;

  sum = 0.0; count = 0;
  for (i=0; i<n; i++) {
    sum += value[i]; count++;
  }

  if (count>0) out = sum / count;
   
  return out;
} /* end StraightAvg */

/**
 * Determine robust RMS value of a ofloat array about mean
 * Use center 1-alpha of points, excluding at least one point from each end
 * \param n       Number of points, needs at least 4
 * \param value   Array of values assumed sorted
 * \param mean    Mean value of value
 * \param alpha   0 -> 1 = central fraction of ordered list to use
 *                0.5 is generally a good compromise 
 * \return RMS value, fblank if cannot determine
 */
static ofloat MedianSigma (ollong n, ofloat *value, ofloat mean, ofloat alpha)
{
  ofloat fblank = ObitMagicF();
  ofloat out;
  ofloat sum, beta;
  olong i, i1, i2, count;

  out = fblank;
  if (n<=4) return out;
  if (mean==fblank) return out;

  beta = MAX (0.05, MIN (0.95, 1.0-alpha)) / 2.0; /*  Average around median factor */

  /* Get RMS around the center 1-alpha */
  i1 = MAX (1,   (n/2)-(olong)(beta*n+0.5));
  i2 = MIN (n-1, (n/2)+(olong)(beta*n+0.5));

  if (i2>i1) {
    sum = 0.0;
    count = 0;
    for (i=i1; i<i2; i++) {
      if (value[i]!=fblank) {
	sum += (value[i]-mean)*(value[i]-mean);
	count++;
      }
    }
    if (count>1) out = sqrt(sum / (count-1));
  }
   
  return out;
} /* end MedianSigma */

/**
 * Determine Weighted linear regression of a spectrum
 * Fitting in log(s)/log(nu) space
 * Negative and/or blanked values ignored
 * \param n       Number of spectral points
 * \param nu      Array of frequencies (Hz)
 * \param S       Array of flux densities, possibly blanked
 * \param sigma   Array of uncertainties for S
 * \param S_0     [out] flux density @ nu[0]
 * \param alpha   [out] spectral index
 */
static void LinRegWt (olong n, odouble *nu, ofloat *S, ofloat *sigma,
		      ofloat *S_0, ofloat* alpha) {
  ofloat fblank = ObitMagicF();
  odouble a, b, nu_0, s_w, s_x, s_y, s_xx, s_yy, s_xy, lnu, ls, w;
  olong i;
  s_w=0.; s_x=0.; s_y=0.; s_xx=0.; s_yy=0.; s_xy=0.; 
  nu_0 = nu[0]; 
  for (i=0; i<n; i++) {
    if ((S[i]>0.0) && (S[i]!=fblank)) {
      lnu = log(nu[i]/nu_0);
      ls  = log(S[i]);
      w   = 1/sigma[i];
      s_x += w*lnu; s_y += w*ls; s_xy += w*lnu*ls;
      s_xx += w*lnu*lnu; s_yy += w*ls*ls; s_w+=w;
    }
  }
  a   = (s_w*s_xy-s_x*s_y)/(s_w*s_xx-s_x*s_x);
  b = (s_y-a*s_x)/s_w;
  *S_0   = exp(b);
  *alpha = a;
} /* end linRegWt */
