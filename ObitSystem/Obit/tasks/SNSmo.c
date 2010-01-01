/* $Id$  */
/* Obit Radio interferometry calibration software                     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2009                                          */
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
#include "ObitUV.h"
#include "ObitData.h"
#include "ObitTableUtil.h"
#include "ObitTableSN.h"
#include "ObitUVSoln.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* SNSmoIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void SNSmoOut (ObitInfoList* outList, ObitErr *err);
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
void SNSmoHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Copy selected data */
void SNSmoCopy (ObitInfoList* myInput, ObitUV* inData,  ObitErr* err);
/* Do Any Clipping */
void SNSmoClip (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* reference phases */
void SNSmoRef (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Smooth */
void SNSmoSmooth (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Clipping routines */
void clipMBDly (ObitTableSN *SNTable, ObitUVSel *sel, 
		ofloat stdela, ofloat mxdela, olong sub,	
		olong maxtim, olong *wrkrec, ofloat *wrktim,
		ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
		ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
		ObitErr* err);
void clipDly (ObitTableSN *SNTable, ObitUVSel *sel, olong iif, 
	      ofloat stdela, ofloat mxdela, olong sub, 
	      olong maxtim, olong *wrkrec, ofloat *wrktim,
	      ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
	      ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
	      ObitErr*  err);
void clipRat (ObitTableSN *SNTable, ObitUVSel *sel, ofloat strate, ofloat mxrate, 
	      olong sub, 
	      olong maxtim, olong *wrkrec, ofloat *wrktim,
	      ofloat *work1, ofloat *work2, ofloat *work4, 
	      ofloat *work5, ofloat *work7, ofloat *work8, 
	      ObitErr* err);
void clipAmp (ObitTableSN *SNTable, ObitUVSel *sel, olong iif, 
	      ofloat stamp, ofloat mxamp, olong sub, 
	      olong maxtim, olong *wrkrec, ofloat *wrktim,
	      ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
	      ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
	      ObitErr* err);

void clipPh (ObitTableSN *SNTable,  ObitUVSel *sel, olong iif, 
	     ofloat stph, ofloat mxph, olong sub, 
	     olong maxtim, olong *wrkrec, ofloat *wrktim,
	     ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
	     ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
	     ObitErr* err);

/* Program globals */
gchar *pgmName = "SNSmo";       /* Program name */
gchar *infile  = "SNSmo.in" ;   /* File with program inputs */
gchar *outfile = "SNSmo.out";   /* File to contain program outputs */
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
/*   Obit task which smooths and filters Solution(SN) tables.             */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitErr      *err= NULL;

  /* Startup - parse command line */
  err = newObitErr();
  myInput = SNSmoIn (argc, argv, err);
  if (err->error) ierr = 1;
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

  /* Copy selected table */

  SNSmoCopy (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Clipping */
  SNSmoClip (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* re reference phases */
  SNSmoRef (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* smooth solutions */
  SNSmoSmooth (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Write history */
  SNSmoHistory (myInput, inData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* SNSmoIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "SNSmoIn";

  /* error checks */
  g_assert(ObitErrIsA(err));
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
} /* end SNSmoIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: SNSmo -input file -output ofile [args]\n");
    fprintf(stderr, "Clip and smooth an SN table \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def SNSmo.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def SNSmo.out\n");
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
  strTemp = "SNSmo.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "SNSmoName";
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
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doCalSelect;
  oint         doSNSmo;
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  /* Make sure doCalSelect set properly */
  doCalSelect = FALSE;
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doSNSmo = -1;
  ObitInfoListGetTest(myInput, "doSNSmo",  &type, dim, &doSNSmo);
  doCalSelect = doCalSelect || (doSNSmo>0);
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);
 


} /* end digestInputs */

/*----------------------------------------------------------------------- */
/*  Get input data, builds selector                                       */
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
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         Aseq, disk, cno, nvis=1000;
  gchar        *Type, *strTemp, inFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "UV";
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "timeRange", "BChan", "EChan",   "BIF", "EIF", "subA",
    "Antennas", "FredID", 
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

  /* Full instantiation - needed to initialize tableList */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Open and close UV data to set selection */
  ObitUVOpen(inData, OBIT_IO_ReadCal, err);
  ObitUVClose(inData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Write History for SNSmo                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void SNSmoHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "FreqID", "timeRange",  "subA", "Antennas", 
    "solnIn", "solnOut", "smoFunc", "smoParm", "doBlank", "smoType", 
    "refAnt",
    NULL};
  gchar *routine = "SNSmoHistory";

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
 
} /* end SNSmoHistory  */


/**
 * Routine to deselect records in an SN table if they match a given
 * subarray, have a selected FQ id, appear on a list of antennas and
 * are in a given timerange.
 * \param SNTab      SN table object 
 * \param isub       Subarray number, <=0 -> any
 * \param fqid       Selected FQ id, <=0 -> any
 * \param nantf      Number of antennas in ants
 * \param ants       List of antennas, NULL or 0 in first -> flag all
 * \param nsou       Number of source ids in sources
 * \param sources    List of sources, NULL or 0 in first -> flag all
 * \param timerange  Timerange to flag, 0s -> all
 * \param err        Error/message stack, returns if error.
 */
void SNSmoCopy (ObitInfoList* myInput, ObitUV* inData,  ObitErr* err)
{
  ObitIOCode retCode;
  ObitTableSN *inTab=NULL, *outTab=NULL;
  ObitTableSNRow *row=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong   itemp, isuba,  fqid;
  olong solnIn, solnOut, highVer;
  gboolean dropIt;
  ofloat tr[2];
  olong  loop, outrow;
  gchar  *routine = "SNSmoCopy";
  
  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  g_assert(ObitUVIsA(inData));
 
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS SN");
  itemp = 0;
  ObitInfoListGetTest(myInput, "solnIn", &type, dim, &itemp);
  solnIn = itemp;
  if (solnIn==0) solnIn = highVer;
  /* Save to inputs */
  dim[0] = dim[1] = dim[2] = 1; itemp = solnIn;
  ObitInfoListAlwaysPut(myInput, "solnIn", OBIT_long, dim, &itemp);

  ObitInfoListGetTest(myInput, "solnOut", &type, dim, &itemp);
  solnOut = itemp;
  if (solnOut==0) solnOut = highVer+1;
  /* Save to inputs */
  dim[0] = dim[1] = dim[2] = 1; itemp = solnOut;
  ObitInfoListAlwaysPut(myInput, "solnOut", OBIT_long, dim, &itemp);


  inTab = newObitTableSNValue (inData->name, (ObitData*)inData, &solnIn, 
			       OBIT_IO_ReadWrite, 0, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  /* Open input table */
  retCode = ObitTableSNOpen (inTab, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, inTab->name);

  outTab = newObitTableSNValue (inData->name, (ObitData*)inData, &solnOut, 
				OBIT_IO_ReadWrite, inTab->numPol, inTab->numIF, 
				err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

   /* If SN table previously existed, deselect values about to be replaced. 
     get information from selector on inData */
  ObitUVSolnDeselSN (outTab, inData->mySel->SubA, inData->mySel->FreqID, 
				inData->mySel->numberAntList, inData->mySel->ants, 
				inData->mySel->numberSourcesList, inData->mySel->sources, 
				inData->mySel->timeRange, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
 
  /* Open output table */
  retCode = ObitTableSNOpen (outTab, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, outTab->name);
  outTab->numAnt = inTab->numAnt;  /* Save number of antennas */
  /* If there are already entries, mark as unsorted */
  if (outTab->myDesc->nrow>0) {
    outTab->myDesc->sort[0] = outTab->myDesc->sort[1] = 0;
    ((ObitTableDesc*)outTab->myIO->myDesc)->sort[0] = 0;
    ((ObitTableDesc*)outTab->myIO->myDesc)->sort[1] = 0;
  }
 
  /* Get selection criteria */
  isuba        = inData->mySel->SubA;
  fqid         = inData->mySel->FreqID;

  /* Timerange */
  tr[0] = inData->mySel->timeRange[0];
  tr[1] = inData->mySel->timeRange[1];
  if ((tr[0]==0.0) && (tr[1]==0.0)) {
    tr[0] = -1.0e20;
    tr[1] =  1.0e20;
  }

  /* Create Row */
  row = newObitTableSNRow (outTab);
  /* Attach row to output buffer */
  ObitTableSNSetRow (outTab, row, err);
  if (err->error) Obit_traceback_msg (err, routine, outTab->name);

  /* Loop through table */
  for (loop=1; loop<=inTab->myDesc->nrow; loop++) { /* loop 20 */

    retCode = ObitTableSNReadRow (inTab, loop, row, err);
    if (err->error) Obit_traceback_msg (err, routine, outTab->name);
    if (row->status<0) continue;  /* Skip deselected records */

    /*  Drop this one? */
    dropIt = (row->SubA!=isuba) && (isuba>0);                    /* by subarray */
    dropIt = dropIt || ((row->FreqID!=fqid) && (fqid>0));        /* by FQ id */
    dropIt = dropIt || ((row->Time<tr[0]) || (row->Time>tr[1])); /* by time */
    if (dropIt) continue; /* Drop */

    if (!ObitUVSelWantAnt (inData->mySel,  row->antNo))  continue;  /* Check Antenna */
    if (!ObitUVSelWantSour (inData->mySel, row->SourID)) continue;  /* Check Source */

    /* Write record to output */
    outrow = -1;
    retCode = ObitTableSNWriteRow (outTab, outrow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, outTab->name);
  } /* End loop over table */

  /* Close tables */
  retCode = ObitTableSNClose (inTab, err);
  if (err->error) Obit_traceback_msg (err, routine, inTab->name);
  retCode = ObitTableSNClose (outTab, err);
  if (err->error) Obit_traceback_msg (err, routine, outTab->name);
  row = ObitTableSNRowUnref (row);  /* Cleanup */
} /* end SNSmoCopy */

void SNSmoClip (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
/*----------------------------------------------------------------------- */
/*  Do Any Clipping                                                       */
/*  Data selection based on selector on inData                            */
/*  Adopted from the AIPSish:  SNSMO.FOR/SNCLIP                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
{
  ObitTableSN *SNTable =  NULL;
  olong SNver, highVer;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat clipSmo[5]  = {0.0, 0.0, 0.0, 0.0, 0.0};
  ofloat stdela, stmbde, strate, stamp, stph;
  ofloat clipParm[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  ofloat mxdela, mxmbde, mxrate, mxamp, mxph;
  gboolean dodela, dombde, dorate, doamp, doph;
  olong itemp, numSub, iSub, subA, bif, eif, i;
  olong mxtim, *wrkrec=NULL;
  ofloat *wrktim=NULL, *work1=NULL, *work2=NULL, *work3=NULL, *work4=NULL,
    *work5=NULL, *work6=NULL, *work7=NULL, *work8=NULL;
  gchar *routine = "SNSmoClip";

  /* Get SN table to clip */
  itemp = 0;
  ObitInfoListGetTest(myInput, "solnOut", &type, dim, &itemp);
  SNver = itemp;
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS SN");
  if (SNver==0) SNver = highVer;
  SNTable = newObitTableSNValue (inData->name, (ObitData*)inData, &SNver, 
				 OBIT_IO_ReadWrite, 0, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Get clipping parameters */
  ObitInfoListGetTest(myInput, "clipSmo",  &type, dim, clipSmo);
  stamp  = clipSmo[0] / 24.0;
  stph   = clipSmo[1] / 24.0;
  strate = clipSmo[2] / 24.0;
  stdela = clipSmo[3] / 24.0;
  stmbde = clipSmo[4] / 24.0;
  ObitInfoListGetTest(myInput, "clipParm", &type, dim, clipParm);
  mxamp  = fabs (clipParm[0]);
  mxrate = fabs (clipParm[2]) / (inData->myDesc->freqIF[0] * 1.0e3);
  mxph   = fabs (clipParm[1]) / 57.296;
  mxdela = fabs (clipParm[3]) * 1.0e-9;
  mxmbde = fabs (clipParm[4]) * 1.0e-9;
  ObitInfoListGetTest(myInput, "BIF",      &type, dim, &bif);
  if (bif==0) bif = 1;
  if (inData->myDesc->jlocif>=0) 
    bif = MAX (1, MIN (bif, inData->myDesc->inaxes[inData->myDesc->jlocif]));
  else bif = 1;
  ObitInfoListGetTest(myInput, "EIF",      &type, dim, &eif);
  if (inData->myDesc->jlocif>=0) {
    if (eif<=0) eif = inData->myDesc->inaxes[inData->myDesc->jlocif];
    eif = MAX (1, MIN (eif, inData->myDesc->inaxes[inData->myDesc->jlocif]));
  } else eif = 1;
  ObitInfoListGetTest(myInput, "subA",     &type, dim, &subA);

  /* What's to be clipped? */
  dodela = mxdela >=  1.0e-10;
  dombde = mxmbde >=  1.0e-10;
  dorate = mxrate >=  1.0e-18;
  doamp  = mxamp  >=  1.0e-10;
  doph   = mxph   >=  1.0e-10;

  /* Anything to do? */
  if (!doamp && !doph && !dombde && !dodela && !dorate)
    {ObitTableSNUnref(SNTable); return;}

  /* Zero smoothing time actually very large */
  if (mxdela < 1.0e-10) mxdela = 1.0e20;
  if (mxmbde < 1.0e-10) mxmbde = 1.0e20;
  if (mxrate < 1.0e-18) mxrate = 1.0e20;
  if (mxamp  < 1.0e-10) mxamp  = 1.0e20;
  if (mxph   < 1.0e-10) mxph   = 1.0e20;

  /* Number of subarrays = number AN tables */
  numSub = ObitTableListGetHigh (inData->tableList, "AIPS AN");

  /* Sort to antenna time order */
  ObitTableUtilSort2f ((ObitTable*)SNTable, "ANTENNA", 1, FALSE, "TIME  ", 
		       1, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Open table */
  ObitTableSNOpen (SNTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Create work arrays - make big enough to swallow whole table */
  mxtim = SNTable->myDesc->nrow;
  wrkrec = g_malloc(mxtim*sizeof(olong));
  wrktim = g_malloc(mxtim*sizeof(ofloat));
  work1  = g_malloc(mxtim*sizeof(ofloat));
  work2  = g_malloc(mxtim*sizeof(ofloat));
  work3  = g_malloc(mxtim*sizeof(ofloat));
  work4  = g_malloc(mxtim*sizeof(ofloat));
  work5  = g_malloc(mxtim*sizeof(ofloat));
  work6  = g_malloc(mxtim*sizeof(ofloat));
  work7  = g_malloc(mxtim*sizeof(ofloat));
  work8  = g_malloc(mxtim*sizeof(ofloat));

  /* Loop over subarrays */
  for (iSub=1; iSub<=numSub; iSub++) {
    /* This one wanted? */
    if ((subA > 0)  &&  (iSub != subA)) continue;
    
    /* Multiband delay */
    if (dombde) clipMBDly (SNTable, inData->mySel, stmbde, mxmbde, iSub, 
			   mxtim, wrkrec, wrktim, work1, work2, work3,
			   work4, work5, work6, work7, work8,
			   err);
    
    /* Single band delay */
    if (dodela) {
      for (i= bif; i<=eif; i++) { /* loop 60 */
	clipDly (SNTable, inData->mySel, i, stdela, mxdela, iSub, 
		 mxtim, wrkrec, wrktim, work1, work2, work3,
		 work4, work5, work6, work7, work8,
		 err);
	if (err->error) goto cleanup;
     } /* end loop  L60:  */;
    } 
    
    /* Clip rates */
    if (dorate) {
      clipRat (SNTable, inData->mySel, strate, mxrate, iSub, 
	       mxtim, wrkrec, wrktim, work1, work2, 
	       work4, work6, work7, work8,
	       err);
      if (err->error) goto cleanup;
    } 
    
    /* Clip amp */
    if (doamp) {
      /* Loop over IF */
      for (i= bif; i<=eif; i++) { /* loop 100 */
	clipAmp (SNTable, inData->mySel, i, stamp, mxamp, iSub, 
		 mxtim, wrkrec, wrktim, work1, work2, work3,
		 work4, work5, work6, work7, work8,
		 err);
	if (err->error) goto cleanup;
      } /* end loop  L100: */;
    }

    /* Clip phase */
    if (doph) {
      /* Loop over IF */
      for (i= bif; i<=eif; i++) { /* loop 100 */
	clipPh (SNTable, inData->mySel, i, stph, mxph, iSub, 
		mxtim, wrkrec, wrktim, work1, work2, work3,
		work4, work5, work6, work7, work8,
		err);
	if (err->error) goto cleanup;
      } /* end loop  L100: */;
    }
  } /* end loop over subarrays */
    
 cleanup:
  /* Close SN Table */
  ObitTableSNClose (SNTable, err);
  SNTable = ObitTableSNUnref(SNTable);
  if (wrkrec) g_free (wrkrec);
  if (wrktim) g_free (wrktim);
  if (work1)  g_free (work1);
  if (work2)  g_free (work2);
  if (work3)  g_free (work3);
  if (work4)  g_free (work4);
  if (work5)  g_free (work5);
  if (work6)  g_free (work6);
  if (work7)  g_free (work7);
  if (work8)  g_free (work8);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
 
} /* end SNSmoClip */

void SNSmoRef (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
/*----------------------------------------------------------------------- */
/*  reference phases                                                      */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
{
  ObitTableSN *SNTable =  NULL;
  olong SNver, highVer;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong itemp, numSub, iSub, subA, refAnt;
  gchar *routine = "SNSmoRef";

  /* Get SN table to rereference */
  itemp = 0;
  ObitInfoListGetTest(myInput, "solnOut", &type, dim, &itemp);
  SNver = itemp;
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS SN");
  if (SNver==0) SNver = highVer;
  SNTable = newObitTableSNValue (inData->name, (ObitData*)inData, &SNver, 
				 OBIT_IO_ReadWrite, 0, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* reference antenna */
  refAnt = 0;
  ObitInfoListGetTest(myInput, "refAnt", &type, dim, &refAnt);

  /* desired subarray */
  numSub = (olong)ObitTableListGetHigh (inData->tableList, "AIPS AN");
  subA = 0;
  ObitInfoListGetTest(myInput, "subA",     &type, dim, &subA);

  /* Sort to time-antenna  order */
  ObitTableUtilSort2f ((ObitTable*)SNTable, "ANTENNA", 1, FALSE, "TIME  ", 
		       1, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Loop over subarrays */
  for (iSub=1; iSub<=numSub; iSub++) {
    /* This one wanted? */
    if ((subA > 0)  &&  (iSub != subA)) continue;
    
    ObitUVSolnRefAnt (SNTable, subA, &refAnt, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
  } /* end loop over subarrays */

  SNTable = ObitTableSNUnref(SNTable);
} /*end SNSmoRef */

/* Smooth */
void SNSmoSmooth (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
/*----------------------------------------------------------------------- */
/*   Smooth                                                               */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
{
  ObitTableSN *SNTable =  NULL;
  olong SNver, highVer;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong itemp, numSub, iSub, subA;
  gchar        *smoParms[] = {
    "smoFunc", "smoParm", "doBlank", "smoType", 
    NULL};
  gchar *routine = "SNSmoSmooth";

  /* Get SN table to rereference */
  itemp = 0;
  ObitInfoListGetTest(myInput, "solnOut", &type, dim, &itemp);
  SNver = itemp;
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS SN");
  if (SNver==0) SNver = highVer;
  SNTable = newObitTableSNValue (inData->name, (ObitData*)inData, &SNver, 
				 OBIT_IO_ReadWrite, 0, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* desired subarray */
  subA = 0;
  numSub = (olong)ObitTableListGetHigh (inData->tableList, "AIPS AN");
  ObitInfoListGetTest(myInput, "subA",  &type, dim, &subA);

  /* copy smoothing parameters to SNTable */
  ObitInfoListCopyList (myInput, SNTable->info, smoParms);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Loop over subarrays */
  for (iSub=1; iSub<=numSub; iSub++) {
    /* This one wanted? */
    if ((subA > 0)  &&  (iSub != subA)) continue;
    
    ObitUVSolnSNSmo (SNTable, subA, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
  } /* end loop over subarrays */

  /* Sort to time-antenna  order */
  ObitTableUtilSort2f ((ObitTable*)SNTable,"TIME  " , 1, FALSE, "ANTENNA", 
		       1, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  SNTable = ObitTableSNUnref(SNTable);
} /* end SNSmoSmooth */

/*----------------------------------------------------------------------- */
/*  Clip multiband delay                                                  */
/* Routine adopted from the AIPSish SNSMO.FOR/CLPDLY                      */
/*   Input:                                                               */
/*      SNTable Table to clip                                             */
/*      sel      Selector telling which data wanted                       */
/*      stdela   Smoothing time (day)                                     */
/*      mxdela   Max. excursion from smoothed value                       */
/*      sub      Subarray (1-rel)                                         */
/*      mxtim    Number of elements in wrkrec, wrktim, work1-8            */
/*      wrktim, wrkrec, work* work arrays                                 */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void clipMBDly (ObitTableSN *SNTable, ObitUVSel *sel, 
		ofloat stdela, ofloat mxdela, olong sub,	
		olong maxtim, olong *wrkrec, ofloat *wrktim,
		ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
		ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
		ObitErr* err)
{
  ObitTableSNRow *row=NULL;
  ofloat fblank =  ObitMagicF();
  olong   loopr, loopa, numtim, ant, fstrec, nleft, save, itime;
  olong  isnrno=0.0;
  gboolean   bad, bad2, want;
  ofloat     diff;
  odouble timoff=0.0;
  gchar *routine = "clipMBDly";

  /* Create Row */
  row = newObitTableSNRow (SNTable);
  /* Attach row to output buffer */
  ObitTableSNSetRow (SNTable, row, err);
  fstrec = 0;  /* Record number read in table */

  /* Loop over antenna */
  for (loopa= 1; loopa<=SNTable->numAnt; loopa++) { /* loop 600 */
    ant = loopa;
    /* Want this antenna? */
    if (!ObitUVSelWantAnt (sel, ant)) continue;

    /* Set pointers, counters */
    numtim = 0;
    nleft = SNTable->myDesc->nrow - fstrec;  /* How many rows? */

    /* Loop in time, reading */
    for (loopr= 1; loopr<=nleft; loopr++) { /* loop 100 */
      isnrno = fstrec + loopr;
      ObitTableSNReadRow (SNTable, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTable->name);
      if (row->status<0) continue;  /* Skip deselected records */
    
      /* Finished antenna? */
      if (row->antNo > ant) break;

      /* See if wanted. */
      want = ObitUVSelWantSour (sel, row->SourID);

      /* Check subarray */
      want = want  &&  (row->SubA == sub);

      /* Not all antennas wanted */
      want = want  &&  (row->antNo == ant);
      if (want) {
	/* See if flagged value */
	bad = (row->MBDelay1 == fblank) || (row->Weight1[0]<=0.0);
	bad2 = (SNTable->numPol <= 1)  ||  (row->MBDelay2 == fblank) ||
	  (row->Weight2[0]<=0.0);
	if (numtim < maxtim) {
	  if (numtim == 0) timoff = row->Time;
	  wrktim[numtim] = row->Time - timoff;
	  wrkrec[numtim] = isnrno;
	  if (bad) {
	    work2[numtim] = fblank;
	    work4[numtim] = 0.0;
	  } else {
	    work2[numtim] = row->MBDelay1;
	    work4[numtim] = row->Weight1[0];
	  } 
	  if (bad2) {
	    work3[numtim] = fblank;
	    work5[numtim] = 0.0;
	  } else {
	    work3[numtim] = row->MBDelay2;
	    work5[numtim] = row->Weight2[0];
	  } 
	  numtim++;
	} 
      } 
    } /* end loop  L100: */;
    save = isnrno - 1; /* How far did we get? */

    if (numtim <= 0) goto endAnt;  /* Catch anything? */

    /* Smooth as requested */
    ObitUVSolnSmooMWF (stdela, 1.0, wrktim, work2, work4, numtim, work1, 
		       work6, work7, work8, FALSE);
    /* Second Poln? */
    if (SNTable->numPol > 1) {
      ObitUVSolnSmooMWF (stdela, 1.0, wrktim, work3, work5, numtim, work2, 
			 work4, work7, work8, FALSE);
    } 

    /* Clip */
    for (itime=0; itime<numtim; itime++) { /* loop 200 */
      isnrno = wrkrec[itime];
      ObitTableSNReadRow (SNTable, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTable->name);

      bad = FALSE;
      diff = fabs (row->MBDelay1 - work1[itime]);
      if (diff > mxdela) {
	row->MBDelay1 = fblank;
	bad = TRUE;
      } 
      /* Second polarization? */
      if (SNTable->numPol >= 2) {
	diff = fabs (row->MBDelay2 - work2[itime]);
	if (diff > mxdela) {
	  row->MBDelay2 = fblank;
	  bad = TRUE;
	} 
      } 
      /* Rewrite record */
      if (bad) ObitTableSNWriteRow (SNTable, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTable->name);
    } /* end loop clipping L200: */;
  endAnt: fstrec = save;
  } /* end loop  of antenna loop L600: */;

  /* Cleanup */
  row = ObitTableSNRowUnref (row);  

} /* end clipMBDly */

/*----------------------------------------------------------------------- */
/*  Clip single band delay                                                */
/* Routine adopted from the AIPSish SNSMO.FOR/CLPDLY                      */
/*   Input:                                                               */
/*      SNTable Table to clip                                             */
/*      sel      Selector telling which data wanted                       */
/*      iif      IF to operate on, 1-rel                                  */
/*      stdela   Smoothing time (day)                                     */
/*      mxdela   Max. excursion from smoothed value  (sec)                */
/*      sub      Subarray (1-rel)                                         */
/*      mxtim    Number of elements in wrkrec, wrktim, work1-8            */
/*      wrktim, wrkrec, work* work arrays                                 */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void clipDly (ObitTableSN *SNTable, ObitUVSel *sel, olong iif, 
	      ofloat stdela, ofloat mxdela, olong sub, 
	      olong maxtim, olong *wrkrec, ofloat *wrktim,
	      ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
	      ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
	      ObitErr*  err)
{
  ObitTableSNRow *row=NULL;
  ofloat fblank =  ObitMagicF();
  olong   loopr, loopa, numtim, ant, fstrec, nleft, save, itime;
  olong  isnrno=0;
  gboolean   bad, bad2, want;
  ofloat     diff;
  odouble timoff=0.0;
  gchar *routine = "clipDly";

  /* Create Row */
  row = newObitTableSNRow (SNTable);
  /* Attach row to output buffer */
  ObitTableSNSetRow (SNTable, row, err);
  fstrec = 0;  /* Record number read in table */

  /* Loop over antenna */
  for (loopa= 1; loopa<=SNTable->numAnt; loopa++) { /* loop 600 */
    ant = loopa;
    /* Want this antenna? */
    if (!ObitUVSelWantAnt (sel, ant)) continue;

    /* Set pointers, counters */
    numtim = 0;
    nleft = SNTable->myDesc->nrow - fstrec;  /* How many rows? */

    /* Loop in time, reading */
    for (loopr= 1; loopr<=nleft; loopr++) { /* loop 100 */
      isnrno = fstrec + loopr;
      ObitTableSNReadRow (SNTable, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTable->name);
      if (row->status<0) continue;  /* Skip deselected records */
    
      /* Finished antenna? */
      if (row->antNo > ant) break;

      /* See if wanted. */
      want = ObitUVSelWantSour (sel, row->SourID);

      /* Check subarray */
      want = want  &&  (row->SubA == sub);

      /* Not all antennas wanted */
      want = want  &&  (row->antNo == ant);
      if (want) {
	/* See if flagged value */
	bad = (row->Delay1[iif-1] == fblank) || (row->Weight1[iif-1]<=0.0);
	bad2 = (SNTable->numPol <= 1)  ||  (row->Delay2[iif-1] == fblank)
	  || (row->Weight2[iif-1]<=0.0);
	if (numtim < maxtim) {
	  if (numtim == 0) timoff = row->Time;
	  wrktim[numtim] = row->Time - timoff;
	  wrkrec[numtim] = isnrno;
	  if (bad) {
	    work2[numtim] = fblank;
	    work4[numtim] = 0.0;
	  } else {
	    work2[numtim] = row->Delay1[iif-1];
	    work4[numtim] = row->Weight1[iif-1];
	  } 
	  if (bad2) {
	    work3[numtim] = fblank;
	    work5[numtim] = 0.0;
	  } else {
	    work3[numtim] = row->Delay2[iif-1];
	    work5[numtim] = row->Weight2[iif-1];
	  } 
	  numtim++;
	} 
      } 
    } /* end loop  L100: */;
    save = isnrno - 1; /* How far did we get? */

    if (numtim <= 0) goto endAnt;  /* Catch anything? */

    /* Smooth as requested */
    ObitUVSolnSmooMWF (stdela, 1.0, wrktim, work2, work4, numtim, work1, 
		       work6, work7, work8, FALSE);
    /* Second Poln? */
    if (SNTable->numPol > 1) {
      ObitUVSolnSmooMWF (stdela, 1.0, wrktim, work3, work5, numtim, work2, 
			 work4, work7, work8, FALSE);
    } 

    /* Clip */
    for (itime=0; itime<numtim; itime++) { /* loop 200 */
      isnrno = wrkrec[itime];
      ObitTableSNReadRow (SNTable, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTable->name);

      bad = FALSE;
      diff = fabs (row->Delay1[iif-1] - work1[itime]);
      if (diff > mxdela) {
	row->Delay1[iif-1] = fblank;
	bad = TRUE;
      } 
      /* Second polarization? */
      if (SNTable->numPol >= 2) {
	diff = fabs (row->Delay2[iif-1] - work2[itime]);
	if (diff > mxdela) {
	  row->Delay2[iif-1] = fblank;
	  bad = TRUE;
	} 
      } 
      /* Rewrite record */
      if (bad) ObitTableSNWriteRow (SNTable, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTable->name);
    } /* end loop clipping L200: */;
  endAnt: fstrec = save;
  } /* end loop  of antenna loop L600: */;

  /* Cleanup */
  row = ObitTableSNRowUnref (row);  

} /* end clipDly */

/*----------------------------------------------------------------------- */
/*  Clip rates, rates averaged before smoothing                           */
/* Routine adopted from the AIPSish SNSMO.FOR/CLPRAT                      */
/*   Input:                                                               */
/*      SNTable Table to clip                                             */
/*      sel      Selector telling which data wanted                       */
/*      strate   Smoothing time (day)                                     */
/*      mxrate   Max. excursion from smoothed value                       */
/*      sub      Subarray (1-rel)                                         */
/*      mxtim    Number of elements in wrkrec, wrktim, work1-8            */
/*      wrktim, wrkrec, work* work arrays                                 */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void clipRat (ObitTableSN *SNTable, ObitUVSel *sel, 
	      ofloat strate, ofloat mxrate, olong sub, 
	      olong maxtim, olong *wrkrec, ofloat *wrktim,
	      ofloat *work1, ofloat *work2, ofloat *work4, 
	      ofloat *work6, ofloat *work7, ofloat *work8, 
	      ObitErr* err)
{
  ObitTableSNRow *row=NULL;
  ofloat fblank =  ObitMagicF();
  olong   loopr, loopa, numtim, ant, fstrec, nleft, save, itime;
  olong  isnrno=0, i;
  gboolean   bad, want;
  ofloat     sum, count, diff;
  odouble timoff=0.0;
  gchar *routine = "clipRat";

  /* Create Row */
  row = newObitTableSNRow (SNTable);
  /* Attach row to output buffer */
  ObitTableSNSetRow (SNTable, row, err);
  fstrec = 0;  /* Record number read in table */

  /* Loop over antenna */
  for (loopa= 1; loopa<=SNTable->numAnt; loopa++) { /* loop 600 */
    ant = loopa;
    /* Want this antenna? */
    if (!ObitUVSelWantAnt (sel, ant)) continue;

    /* Set pointers, counters */
    numtim = 0;
    nleft = SNTable->myDesc->nrow - fstrec;  /* How many rows? */

    /* Loop in time, reading */
    for (loopr= 1; loopr<=nleft; loopr++) { /* loop 100 */
      isnrno = fstrec + loopr;
      ObitTableSNReadRow (SNTable, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTable->name);
      if (row->status<0) continue;  /* Skip deselected records */
    
      /* Finished antenna? */
      if (row->antNo > ant) break;

      /* See if wanted. */
      want = ObitUVSelWantSour (sel, row->SourID);

      /* Check subarray */
      want = want  &&  (row->SubA == sub);

      /* Not all antennas wanted */
      want = want  &&  (row->antNo == ant);
      if (want) {
	/* Average all rates */
	sum   = 0.0;
	count = 0.0;
	for (i=0; i<SNTable->numIF; i++) { /* loop 20 */
	  if ((row->Weight1[i] > 0.0) && (row->Rate1[i] != fblank)) {
	    sum   += row->Rate1[i];
	    count += 1.0;
	  }
	  if ((SNTable->numPol > 1)  &&  (row->Weight2[i] > 0.0) && (row->Rate2[i] != fblank)) {
	    sum   += row->Rate2[i];
	    count += 1.0;
	  } 
	} /* end loop  L20:  */;

	/* See if flagged value */
	bad = count  <=  0.1;

	/* add to work arrays */	
	if (numtim < maxtim) {
	  if (numtim == 0) timoff = row->Time;
	  wrktim[numtim] = row->Time - timoff;
	  wrkrec[numtim] = isnrno;
	  if (bad) {
	    work2[numtim] = fblank;
	    work4[numtim] = 0.0;
	  } else {
	    work2[numtim] = sum/count;
	    work4[numtim] = row->Weight1[0];
	  } 
	  numtim++;
	} 
      } 
    } /* end loop  L100: */;
    save = isnrno - 1; /* How far did we get? */

    if (numtim <= 0) goto endAnt;  /* Catch anything? */

    /* Smooth as requested */
    ObitUVSolnSmooMWF (strate, 1.0, wrktim, work2, work4, numtim, work1, 
		       work6, work7, work8, FALSE);

    /* Clip */
    for (itime=0; itime<numtim; itime++) { /* loop 200 */
      isnrno = wrkrec[itime];
      ObitTableSNReadRow (SNTable, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTable->name);

      bad = FALSE;
      for (i=0; i<SNTable->numIF; i++) { /* loop 20 */

	diff = fabs (row->Rate1[i] - work1[itime]);
	if (diff > mxrate) {
	  row->Rate1[i] = fblank;
	  bad = TRUE;
	} 
	/* Second polarization? */
	if (SNTable->numPol >= 2) {
	  diff = fabs (row->Rate2[i] - work2[itime]);
	  if (diff > mxrate) {
	    row->Rate2[i] = fblank;
	    bad = TRUE;
	  } 
	} 
      } /* end loop over IF */

      /* Rewrite record */
      if (bad) ObitTableSNWriteRow (SNTable, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTable->name);
    } /* end loop clipping L200: */;
  endAnt: fstrec = save;
  } /* end loop  of antenna loop L600: */;

  /* Cleanup */
  row = ObitTableSNRowUnref (row);  

} /* end clipRat */

/*----------------------------------------------------------------------- */
/*  Clip amp                                                              */
/* Routine adopted from the AIPSish SNSMO.FOR/ CLPAPH                     */
/*   Input:                                                               */
/*      SNTable  Table to clip                                            */
/*      iif      IF to operate on, 1-rel                                  */
/*      sel      Selector telling which data wanted                       */
/*      stamp    Smoothing time (day) for amplitude                       */
/*      mxamp    Max. amplitude excursion from smoothed value             */
/*      sub      Subarray (1-rel)                                         */
/*      mxtim    Number of elements in wrkrec, wrktim, work1-8            */
/*      wrktim, wrkrec, work* work arrays                                 */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
 void clipAmp (ObitTableSN *SNTable, ObitUVSel *sel, olong iif, 
	       ofloat stamp, ofloat mxamp, olong sub, 
	       olong maxtim, olong *wrkrec, ofloat *wrktim,
	       ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
	       ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
	       ObitErr* err)
   {
  ObitTableSNRow *row=NULL;
  ofloat fblank =  ObitMagicF();
  olong   loopr, loopa, numtim, ant, fstrec, nleft, save, itime;
  olong  isnrno=0;
  gboolean   bad, bad2, want;
  ofloat     diff;
  odouble timoff=0.0;
  gchar *routine = "clipAmp";

  /* Create Row */
  row = newObitTableSNRow (SNTable);
  /* Attach row to output buffer */
  ObitTableSNSetRow (SNTable, row, err);
  fstrec = 0;  /* Record number read in table */

  /* Loop over antenna */
  for (loopa= 1; loopa<=SNTable->numAnt; loopa++) { /* loop 600 */
    ant = loopa;
    /* Want this antenna? */
    if (!ObitUVSelWantAnt (sel, ant)) continue;

    /* Set pointers, counters */
    numtim = 0;
    nleft = SNTable->myDesc->nrow - fstrec;  /* How many rows? */

    /* Loop in time, reading */
    for (loopr= 1; loopr<=nleft; loopr++) { /* loop 100 */
      isnrno = fstrec + loopr;
      ObitTableSNReadRow (SNTable, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTable->name);
      if (row->status<0) continue;  /* Skip deselected records */
    
      /* Finished antenna? */
      if (row->antNo < ant) continue;
      if (row->antNo > ant) break;

      /* See if wanted. */
      want = ObitUVSelWantSour (sel, row->SourID);

      /* Check subarray */
      want = want  &&  (row->SubA == sub);

      /* Not all antennas wanted */
      want = want  &&  (row->antNo == ant);
      if (want) {
	/* See if flagged value */
	bad = (row->Real1[iif-1] == fblank) || (row->Weight1[iif-1]<=0.0);
	bad2 = (SNTable->numPol <= 1)  ||  (row->Real2[iif-1] == fblank)
	  || (row->Weight2[iif-1]<=0.0);
	if (numtim < maxtim) {
	  if (numtim == 0) timoff = row->Time;
	  wrktim[numtim] = row->Time - timoff;
	  wrkrec[numtim] = isnrno;
	  if (bad) {
	    work2[numtim] = fblank;
	    work4[numtim] = 0.0;
	  } else {
	    work2[numtim] = sqrt (row->Real1[iif-1]*row->Real1[iif-1] + 
				  row->Imag1[iif-1]*row->Imag1[iif-1]);
	    work4[numtim] = row->Weight1[iif-1];
	  } 
	  if (bad2) {
	    work3[numtim] = fblank;
	    work5[numtim] = 0.0;
	  } else {
	    work3[numtim] = sqrt (row->Real2[iif-1]*row->Real2[iif-1] + 
				  row->Imag2[iif-1]*row->Imag2[iif-1]);
	    work5[numtim] = row->Weight2[iif-1];
	  } 
	  numtim++;
	} 
      } 
    } /* end loop  L100: */;
    save = isnrno - 1; /* How far did we get? */
    
    if (numtim <= 0) goto endAnt;  /* Catch anything? */

    /* Smooth as requested */
    ObitUVSolnSmooMWF (stamp, 1.0, wrktim, work2, work4, numtim, work1, 
		       work6, work7, work8, FALSE);
    /* Second Poln? */
    if (SNTable->numPol > 1) {
      ObitUVSolnSmooMWF (stamp, 1.0, wrktim, work3, work5, numtim, work2, 
			 work4, work7, work8, FALSE);
    } 

    /* Clip */
    for (itime=0; itime<numtim; itime++) { /* loop 200 */
      isnrno = wrkrec[itime];
      ObitTableSNReadRow (SNTable, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTable->name);

      bad = FALSE;
      diff = fabs (sqrt (row->Real1[iif-1]*row->Real1[iif-1] + 
			 row->Imag1[iif-1]*row->Imag1[iif-1]) 
		   - work1[itime]);
      if (diff > mxamp) {
	row->Real1[iif-1] = fblank;
	row->Imag1[iif-1] = fblank;
	bad = TRUE;
      } 
      /* Second polarization? */
      if (SNTable->numPol >= 2) {
	diff = fabs (sqrt (row->Real2[iif-1]*row->Real2[iif-1] + 
			   row->Imag2[iif-1]*row->Imag2[iif-1]) 
		     - work2[itime]);
	if (diff > mxamp) {
	  row->Real2[iif-1] = fblank;
	  row->Imag2[iif-1] = fblank;
	  bad = TRUE;
	} 
      } 
      /* Rewrite record */
      if (bad) ObitTableSNWriteRow (SNTable, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTable->name);
    } /* end loop clipping L200: */;
  endAnt: fstrec = save;
  } /* end loop  of antenna loop L600: */;
  
  /* Cleanup */
  row = ObitTableSNRowUnref (row);  
  
   } /* end clipAmp */
 
/*----------------------------------------------------------------------- */
/*  Clip phase                                                            */
/* Not yet implemented - the AIPSish version doesn't look like it works   */
/* since DOPH is never tested                                             */
/* Routine adopted from the AIPSish SNSMO.FOR/CLPAPH                      */
/*   Input:                                                               */
/*      SNTable  Table to clip                                            */
/*      iif      IF to operate on, 1-rel                                  */
/*      sel      Selector telling which data wanted                       */
/*      stamp    Smoothing time (day) for amplitude                       */
/*      mxamp    Max. amplitude excursion from smoothed value             */
/*      doamp    Clip by Amplitude?                                       */
/*      stph     Smoothing time (day) for phases                          */
/*      mxph     Max. excursion from smoothed value                       */
/*      doph     Clip by Phase?                                           */
/*      sub      Subarray (1-rel)                                         */
/*      mxtim    Number of elements in wrkrec, wrktim, work1-8            */
/*      wrktim, wrkrec, work* work arrays                                 */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
  void clipPh (ObitTableSN *SNTable, ObitUVSel *sel, olong iif, 
	       ofloat stph, ofloat mxph, olong sub, 
	       olong maxtim, olong *wrkrec, ofloat *wrktim,
	       ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
	       ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
	       ObitErr* err)
{
  gchar *routine = "clipPh";
  
  /* Give excuses and leave */
  Obit_log_error(err, OBIT_InfoWarn, "%s: Clipping Phases not yet implemented", 
		 routine);
} /* end clipPh */


