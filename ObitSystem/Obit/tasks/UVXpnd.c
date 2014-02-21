/* $Id$  */
/* Obit Task to Task to expand uv data in frequency                   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
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

#include "ObitUV.h"
#include "ObitUVUtil.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitTableANUtil.h"
#include "ObitTableCLUtil.h"
#include "ObitTableSNUtil.h"
#include "ObitTableSUUtil.h"
#include "ObitTableFQUtil.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* UVXpndIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void UVXpndOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Digest inputs */
void digestInputs(ObitInfoList *myInput, ObitErr *err);
/* Get input data */
ObitUV* getTargetData (ObitInfoList *myInput, ObitErr *err);
/* Get input data */
ObitUV* getInputData (ObitInfoList *myInput, ObitErr *err);
/* Create output uvdata */
ObitUV* setOutputUV (ObitInfoList *myInput, ObitUV* inData, 
		     ObitUV* targData, ObitErr *err);
/* process */
void doUVXpnd (ObitInfoList* myInput, ObitUV* inData, 
	       ObitUV* targData, ObitErr* err);
void dataUVXpnd (ObitUV* inData, ObitUV* targData, olong nChXpnd, 
		ObitErr* err);
/* Write history */
void UVXpndHistory (ObitInfoList* myInput, ObitUV* inData, 
		    ObitUV* outData, ObitErr* err);
/* Fix FQ table */
static void FQSel (ObitUV *inUV, ObitUV *targUV, olong chXpn, 
		   olong fqid, ObitErr *err);

/* Program globals */
gchar *pgmName = "UVXpnd";       /* Program name */
gchar *infile  = "UVXpnd.in" ;   /* File with program inputs */
gchar *outfile = "UVXpnd.out";   /* File to contain program outputs */
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
/*   Obit Task to expand uv data in frequency                             */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL, *targData = NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = UVXpndIn (argc, argv, err);
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

  /* Get target uvdata */
  targData = getTargetData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Process */
  doUVXpnd (myInput, inData, targData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  targData  = ObitUnref(targData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* UVXpndIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "UVXpndIn";

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
      
    } else if (strcmp(arg, "-outSeq") == 0) { /* AIPS image sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outSeq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-outDisk") == 0) { /* output image disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outDisk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-outName") == 0) { /* AIPS image outName */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outName", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-outClass") == 0) { /* AIPS image outClass */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outClass", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-outFile") == 0) { /*outFile */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outFile", OBIT_string, dim, strTemp);

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
} /* end UVXpndIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: UVXpnd -input file -output ofile [args]\n");
    fprintf(stderr, "UVXpnd expand uv data in frequency\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def UVXpnd.in\n");
    fprintf(stderr, "  -output output result file, def UVXpnd.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input image\n");
    fprintf(stderr, "  -inFile input FITS UV file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input image (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -outFile output image (FITS Image file\n");  
    fprintf(stderr, "  -outName output image (AIPS file name\n");
    fprintf(stderr, "  -outClass output image (AIPS file class\n");
    fprintf(stderr, "  -outSeq output image (AIPS file sequence\n");
    fprintf(stderr, "  -outDisk output image ((AIPS or FITS) disk number (1-rel) \n");
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
/*     inName    Str [12]   input AIPS uv name  [no def]                  */
/*     inClass   Str [6]    input AIPS uv class  [no def]                 */
/*     inSeq     Int        input AIPS uv sequence no  [no def]           */
/*     outDisk   Int        output AIPS or FITS image disk no  [def 1]    */
/*     outFile   Str [?]    output FITS image file name [def "Image.fits" */
/*     outName   Str [12]   output AIPS image name  [no def]              */
/*     outClass  Str [6]    output AIPS image class  [sClean]             */
/*     outSeq    Int        output AIPS image sequence no  [new]          */
/*     Sources   Str (16,1) Sources selected, blank = all                 */
/*     Stokes    Str (4)    Stokes parameter to image, def=I              */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*     doCalSelect Boo (1)  Apply calibration/selection?  def=True        */
/*     doCalib   Int (1)    >0 => apply calibration, 2=> cal. wt, def=-1  */
/*     gainUse   Int (1)    Gain table (CL/SN) table to apply, 0=> highest*/
/*     doBand    Int (1)    If >0.5 apply bandpass cal.                   */
/*     flagVer   Int (1)    Flagging table version, def=0                 */
/*     BPVer     Int (1)    Bandpass table version, 0=highest, def=0      */
/*     doPol     Boo (1)    Apply polarization calibration?, def=False    */
/*     doFull    Boo (1)    Make full field (flattened) image? def=True   */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat ftemp, farray[3];
  gboolean btemp;
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
  strTemp = "UVXpnd.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "UVXpndName";
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

  /* output FITS file name */
  strTemp = "UVXpndOut.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS  file name */
  strTemp = "            ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS  file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "outSeq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output AIPS or FITS disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "outDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Sources selected, blank = all */
  strTemp = "                ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Sources", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* Stokes parameter to image */
  strTemp = "I   ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Stokes", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Field of view in deg def = 0.0 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "FOV", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Range in uv plane in klambda, 0=>all */
  dim[0] = 2;dim[1] = 1;
  farray[0] = 0.0; farray[1] = 0.0;
  ObitInfoListPut (out, "UVRange", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Timerange in days */
  dim[0] = 2;dim[1] = 1;
  farray[0] = -1.0e20; farray[1] = 1.0e20;
  ObitInfoListPut (out, "timeRange", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Apply calibration/selection?, def=False */
  dim[0] = 1; dim[1] = 1;
  btemp = TRUE;
  ObitInfoListPut (out, "doCalSelect", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  >0 => apply gain calibration, 2=> cal. wt, def=no cal. */
  dim[0] = 1;dim[1] = 1;
  itemp = -1; 
  ObitInfoListPut (out, "doCalib", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Gain table (CL/SN) table to apply, 0=> highest, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "gainUse", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  If >0.5 apply bandpass cal, def = no BP cal. */
  dim[0] = 1;dim[1] = 1;
  itemp = -1; 
  ObitInfoListPut (out, "doBand", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Bandpass table version, 0=highest, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "BPVer", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Flagging table version, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "flagVer", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* Apply polarization calibration?, def= False */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListPut (out, "doPol", OBIT_bool, dim, &btemp, err);
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
  gboolean doCalWt;
  olong doCalib;
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* If doCalWt set, make sure doCalib OK */
  doCalWt = FALSE;
  ObitInfoListGetTest(myInput, "doCalWt",  &type, dim, &doCalWt);
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doCalib",  &type, dim, &doCalib);
  if (doCalWt && (doCalib>0)) doCalib = 2;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "doCalib", OBIT_long, dim, &doCalib);

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
  olong        nvis=1000;
  oint         doCalib;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doCalSelect;
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "timeRange", "BIF", "EIF", "BChan", "EChan", "subA", "Sources", "Qual",
    "souCode", "Stokes", "FredID", "doCalSelect", 
    "doCalib", "gainUse", "doBand", "BPVer", "flagVer", "doPol", "PDVer", 
    NULL
  };
  gchar *routine = "getInputData";
  
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Build basic input UV data Object */
  inData = ObitUVFromFileInfo ("in", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Set buffer size */
  nvis = 1;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO",  OBIT_long, dim,  &nvis);
    
  /* Make sure doCalSelect set properly */
  doCalSelect = FALSE;
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doCalib",  &type, dim, &doCalib);
  doCalSelect = doCalSelect || (doCalib>0);
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);
 
  /* Set number of vis per IO */
  nvis = 1000;  /* How many vis per I/O? */
  nvis =  ObitUVDescSetNVis (inData->myDesc, myInput, nvis);
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO", OBIT_long, dim,  &nvis);

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Ensure inData fully instantiated, with selection and OK */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Get target data                                                       */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitUV with input data                                           */
/*----------------------------------------------------------------------- */
ObitUV* getTargetData (ObitInfoList *myInput, ObitErr *err)
{
  ObitUV       *inData = NULL;
  ObitInfoType type;
  olong        nvis=1000;
  oint         doCalib;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doCalSelect;
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    NULL
  };
  gchar *routine = "getTargetData";
  
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Build basic input UV data Object */
  inData = ObitUVFromFileInfo ("in2", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Set buffer size */
  nvis = 1;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO",  OBIT_long, dim,  &nvis);
    
  /* Make sure doCalSelect set properly */
  doCalSelect = FALSE;
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doCalib",  &type, dim, &doCalib);
  doCalSelect = doCalSelect || (doCalib>0);
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);
 
  /* Set number of vis per IO */
  nvis = 1000;  /* How many vis per I/O? */
  nvis =  ObitUVDescSetNVis (inData->myDesc, myInput, nvis);
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO", OBIT_long, dim,  &nvis);

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Ensure inData fully instantiated, with selection and OK */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getTargetData */

/*----------------------------------------------------------------------- */
/*  Create output uv data                                                 */
/*   Input:                                                               */
/*      Source    Source name                                             */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV from which to clone output                 */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/* Returns the output uv data                                             */
/*----------------------------------------------------------------------- */
ObitUV* setOutputUV (ObitInfoList *myInput, ObitUV* inData, ObitUV* targData, 
		     ObitErr *err)
{
  ObitUV    *outUV = NULL;
  ObitInfoType type;
  olong      i, n, Aseq, disk, cno;
  gchar     *Type, *strTemp, outFile[129];
  gchar     Aname[13], Aclass[7], *Atype = "UV";
  olong     nvis;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean  exist;
  gchar     tname[129];
  gchar     *today=NULL;
  gchar     *FITS = "FITS";
  gchar     *outParms[] = {  /* Parameters for output data */
    "Compress", NULL};
  gchar     *routine = "setOutputUV";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Create basic output UV Object */
  g_snprintf (tname, 100, "output UV data");
  outUV = newObitUV(tname);
    
  /* Output File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4)))
    ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */

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

    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    /* Default out class is "UVXpnd" */
    if (!strncmp(Aclass, "      ", 6)) strncpy (Aclass, "UVXpnd", 7);

    /* input AIPS disk - default is outDisk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (disk<=0)
       ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* output AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Allocate catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
    /* define object */
    nvis = 1000;
    ObitUVSetAIPS (outUV, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output AIPS UV data %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);
    
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

    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    
    /* define object */
    nvis = 1000;
    ObitUVSetFITS (outUV, nvis, disk, outFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output FITS UV data %s on disk %d", outFile, disk);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outUV;
  }

  /* Modify descriptor for number of output channels */
  /* copy input Descriptor */
  outUV->myDesc = ObitUVDescCopy(inData->myDesc, outUV->myDesc, err);
  /* Same number of channels, ref, channel width as target */
  outUV->myDesc->inaxes[outUV->myDesc->jlocf] =
    targData->myDesc->inaxes[targData->myDesc->jlocf];
  outUV->myDesc->crval[outUV->myDesc->jlocf] =
    targData->myDesc->crval[targData->myDesc->jlocf];
  outUV->myDesc->cdelt[outUV->myDesc->jlocf] =
    targData->myDesc->cdelt[targData->myDesc->jlocf];
  
  /* Creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);

  /* Get input parameters from myInput, copy to outUV */
  ObitInfoListCopyList (myInput, outUV->info, outParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

   /* Ensure inData fully instantiated, with selection and OK */
  ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  ObitUVClose (outUV, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

  ObitErrLog(err); /* Show messages */
  return outUV;
} /* end setOutputUV */

/*----------------------------------------------------------------------- */
/*  Upper level expand data                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV                                            */
/*      targData  Target ObitUV                                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doUVXpnd (ObitInfoList* myInput, ObitUV* inData, ObitUV* targData, 
	       ObitErr* err)
{
  ObitUV       *outData=NULL;
  ObitUVDesc *inDesc, *outDesc, *targDesc;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        doCalib=-1, doBand=-1, flagVer=-1, nChXpnd;
  ofloat       ChXpnd;
  gchar *FGInclude[] = {"AIPS FG", NULL};
  gchar *BPInclude[] = {"AIPS BP", NULL};
  gchar        *outParms[] = {  /* Parameters for output data */
    "Compress",
    NULL };
  gchar *routine = "doUVXpnd";

   /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Check compatability */
  inDesc   = inData->myDesc;
  targDesc = targData->myDesc;
  Obit_return_if_fail ((inDesc->inaxes[inDesc->jlocif]==targDesc->inaxes[targDesc->jlocif]), err, 
			"%s Input and target different no. IFs %d %d", 
			routine, inDesc->inaxes[inDesc->jlocif], targDesc->inaxes[targDesc->jlocif]);
  /* Check that integral number of channels in output */
  ChXpnd  = inDesc->cdelt[inDesc->jlocf] / targDesc->cdelt[targDesc->jlocf];
  nChXpnd = (olong)(ChXpnd+0.01);
  Obit_return_if_fail ((fabs(ChXpnd-nChXpnd)<0.01), err, 
			"%s Input and target channel widths incompatable %f %f", 
			routine, inDesc->cdelt[inDesc->jlocf], targDesc->cdelt[targDesc->jlocf]);

  /* Create output file for data */
  outData = setOutputUV (myInput, inData, targData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Check compatability */
  outDesc   = outData->myDesc;
  targDesc = targData->myDesc;
    Obit_return_if_fail (((outDesc->lrec==targDesc->lrec) || (outData==NULL)), err, 
		       "%s Input and output data records incompatable", routine);  

  /* Get output parameters from myInput, copy to outData */
  ObitInfoListCopyList (myInput, outData->info, outParms);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Expand to output */
  dataUVXpnd (inData, outData, nChXpnd, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Copy calibration tables not applied */
  ObitInfoListGetTest (myInput, "doCalib", &type, dim, &doCalib); 
  if (doCalib<=0) {
    ObitTableCLSelect (inData, outData, err);
    ObitTableSNSelect (inData, outData, err);
  }
  ObitInfoListGetTest (myInput, "doBand", &type, dim, &doBand); 
  if (doBand<=0) 
    ObitUVCopyTables (inData, outData, NULL, BPInclude, err);
  ObitInfoListGetTest (myInput, "flagVer", &type, dim, &flagVer); 
  if (flagVer<0) 
    ObitUVCopyTables (inData, outData, NULL, FGInclude, err);

  /* Antenna table(s) */
  ObitTableANSelect (inData, outData, err);

  /* FQ tables */
  ObitTableFQSelect (inData, outData, NULL, 0.0, err);
  FQSel (outData, targData, nChXpnd, 1, err);  /* Fix up chan widths */

  /* Multisource out? Copy Source table */
  if (outData->myDesc->ilocsu>=0) {
     ObitTableSUSelect (inData, outData, err);

     /* Index output  */
     ObitUVUtilIndex (outData, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, inData->name);  

  /* Do history */
  UVXpndHistory (myInput, inData, outData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);  
  
  /* Cleanup */
  outData  = ObitUVUnref(outData);

}  /* end doUVXpnd */

/*----------------------------------------------------------------------- */
/*  Write History for UVXpnd                                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outData   ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void UVXpndHistory (ObitInfoList* myInput, ObitUV* inData, 
		   ObitUV* outData, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "in2File",  "in2Disk", "in2Name", "in2Class", "in2Seq",
    "FreqID", "Sources", "souCode", "Qual", 
    "outDType", "outFile",  "outDisk", "outName", "outClass", "outSeq",
    "BIF", "EIF", "BChan", "EChan", "UVRange",  "timeRange",  "Compress", 
    "doCalSelect",  "doCalib",  "gainUse",  "doBand ", "BPVer", "flagVer", 
    "doPol", "PDVer", "doWtCal",  "corrType", 
    NULL};
  gchar *routine = "UVXpndHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Do history  */
  inHistory  = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory ((ObitData*)outData, OBIT_IO_WriteOnly, err);

  /* If FITS copy header */
  if (inHistory->FileType==OBIT_IO_FITS) {
    ObitHistoryCopyHeader (inHistory, outHistory, err);
  } else { /* simply copy history */
    ObitHistoryCopy (inHistory, outHistory, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end UVXpndHistory  */

/*----------------------------------------------------------------------- */
/*  Lower level expand data                                               */
/*   Input:                                                               */
/*      inData    Input ObitUV                                            */
/*      outData   Output ObitUV                                           */
/*      nChXpnd   Number of channels to expand from input to output       */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void dataUVXpnd (ObitUV* inData, ObitUV* outData, olong nChXpnd,
	        ObitErr* err)
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean doCalSelect;
  ObitIOCode iretCode, oretCode;
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  ofloat *inBuffer, *outBuffer, wtFact;
  olong i, j, k, indx, jndx, jf, jif, js, ojf;
  olong nstok, nchan, nif, iincs, iincf, iincif, oincs, oincf, oincif;
  gchar *routine = "dataUVXpnd";

   /* error checks */
  if (err->error) return;

  doCalSelect = FALSE;
  ObitInfoListGetTest(inData->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;
  nChXpnd = MAX (1, nChXpnd);
  wtFact = 1.0 / nChXpnd;
  
  /* Open files */
  iretCode = ObitUVOpen (inData, access, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* test open output */
  oretCode = ObitUVOpen (outData, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outData, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    outData->buffer = NULL;
    outData->bufferSize = 0;
    Obit_traceback_msg (err, routine, outData->name);
  }

  /* reset to beginning of uv data */
  ObitIOSet (inData->myIO,  inData->info, err);
  ObitIOSet (outData->myIO, outData->info, err);
  if (err->error) Obit_traceback_msg (err, routine,inData->name);

  /* Get descriptors */
  inDesc  = inData->myDesc;
  outDesc = outData->myDesc;
  /* Data increments */
  nstok  = outDesc->inaxes[outDesc->jlocs];
  nchan  = inDesc->inaxes[outDesc->jlocf];
  nif    = inDesc->inaxes[outDesc->jlocif];
  iincs  = inDesc->incs;
  iincf  = inDesc->incf;
  iincif = inDesc->incif;
  oincs  = outDesc->incs;
  oincf  = outDesc->incf;
  oincif = outDesc->incif;

  /* we're in business, expand data */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inData, inData->buffer, err);
    else             iretCode = ObitUVRead (inData, inData->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
   /* How many */
    outDesc->numVisBuff = inDesc->numVisBuff;

    /* Modify data */
    inBuffer  = inData->buffer;
    outBuffer = outData->buffer;
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      /* Zero output */
      for (j=inDesc->nrparm; j<outDesc->lrec; j++) outBuffer[j] = 0.0;
      /* Copy random parameters */
      for (j=0; j<inDesc->nrparm; j++) outBuffer[j] = inBuffer[j];
      /* Expand data */
      for (js=0; js<nstok; js++) {  /* Stokes loop */
	for (jif=0; jif<nif; jif++) {  /* IF loop */
	  for (jf=0; jf<nchan; jf++) {  /* Frequency loop */
	    jndx = inDesc->nrparm + js*iincs + jif*iincif + jf*iincf;
	    /* Expansion loop */
	    for (k=0; k<nChXpnd; k++) {
	      ojf  = jf*nChXpnd + k;
	      indx = outDesc->nrparm + js*oincs + jif*oincif + ojf*oincf;
	      outBuffer[indx]   = inBuffer[jndx];
	      outBuffer[indx+1] = inBuffer[jndx+1];
	      outBuffer[indx+2] = inBuffer[jndx+2] * wtFact;
	    } /* end expansion loop */
	  } /* end Frequency loop */
	} /* end IF loop */
      } /* end Stokes loop */
      inBuffer  += inDesc->lrec;  /* Update buffer pointers */
      outBuffer += outDesc->lrec;
    } /* end loop over visibilities */

    /* Write */
    oretCode = ObitUVWrite (outData, outData->buffer, err);
  } /* end loop processing data */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine,inData->name);
  
  /* close files */
  iretCode = ObitUVClose (inData, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inData->name);
  
  oretCode = ObitUVClose (outData, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, outData->name);
  
}  /* end dataUVXpnd */

/** 
 * Update FQ table for expanding frequency 
 * \param inUV   Input UV data, to be modified
 * \param targUV Target UV data
 * \param chXpn  Number of channels averaged
 * \param fqid   Desired FQ ID to update
 * \param err    Error stack, returns if not empty.
 */
static void FQSel (ObitUV *inUV, ObitUV *targUV, olong chXpn, 
		   olong fqid, ObitErr *err)
{
  ObitTableFQ    *inTab=NULL, *targTab=NULL;
  olong iFQver, highFQver;
  oint numIF;
  olong nif, nchan;
  odouble *freqOff=NULL;
  ofloat *chBandw=NULL;
  oint *sideBand=NULL;
  gchar *FQType = "AIPS FQ";
  gchar *routine = "ObitUVUtil:FQSel";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));

  /* How many target FQ tables  */
  highFQver = ObitTableListGetHigh (targUV->tableList, FQType);

  /* Are there any? */
  if (highFQver <= 0) return;

  /* Should only be one FQ table */
  iFQver = 1;
  nchan = inUV->myDesc->inaxes[inUV->myDesc->jlocf];
  if (inUV->myDesc->jlocif>=0) 
    nif = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
  else
    nif = 1;

  /* Get target table */
  numIF = 0;
  targTab = 
    newObitTableFQValue (targUV->name, (ObitData*)targUV, &iFQver, OBIT_IO_ReadOnly, 
			 numIF, err);
  if (err->error) Obit_traceback_msg (err, routine, targUV->name);
  /* Find it? */
  Obit_return_if_fail(((targTab!=NULL) || (nif<=1)), err,
		      "%s: Could not find FQ table for %s %d IFs", 
		      routine, targUV->name, nif);
  
  /* Fetch values */
  ObitTableFQGetInfo (targTab, fqid, &nif, &freqOff, &sideBand, &chBandw, err);
  if (err->error) Obit_traceback_msg (err, routine, inTab->name);
  
  /* Update channel widths
  for (i=0; i<nif; i++) chBandw[i] /= (ofloat)chXpn; */
  
  /* Get input file table */
  numIF = 0;
  inTab = 
    newObitTableFQValue (inUV->name, (ObitData*)inUV, &iFQver, OBIT_IO_ReadOnly, 
			 numIF, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  /* Find it? */
  Obit_return_if_fail(((inTab!=NULL) || (nif<=1)), err,
		      "%s: Could not find FQ table for %s %d IFs", 
		      routine, inUV->name, nif);
  
  /* Save values */
  ObitTableFQPutInfo (inTab, fqid, nif, freqOff, sideBand, chBandw, err);
  if (err->error) Obit_traceback_msg (err, routine, inTab->name);
  
   /* Cleanup */
  if (freqOff)  g_free(freqOff);
  if (sideBand) g_free(sideBand);
  if (chBandw)  g_free(chBandw);
  inTab   = ObitTableFQUnref(inTab);
  targTab = ObitTableFQUnref(targTab);
  
} /* end FQSel */
