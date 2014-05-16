/* $Id$  */
/* Obit task to Map beam polarization                                 */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009-2013                                          */
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

#include "ObitImageUtil.h"
#include "ObitUVUtil.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitTableSUUtil.h"
#include "ObitTableANUtil.h"
#include "ObitPrecess.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* MapBeamIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void MapBeamOut (ObitInfoList* outList, ObitErr *err);
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

/* Create output image */
ObitImage* setOutput (gchar *Source, olong iStoke, olong ant, 
		      gboolean doRMS, gboolean doPhase,
		      ObitInfoList *myInput, ObitUV* inData, 
		      ObitErr *err);

/* Loop over sources */
void doSources (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);

/* Loop over Channels/Poln */
void doChanPoln (gchar *Source, olong ant, ObitInfoList* myInput, 
		 ObitUV* inData, ObitErr* err);

/* Image */
ObitFArray** doImage (gboolean doRMS, gboolean doPhase, olong ant, 
		      ObitInfoList* myInput, 
		      ObitUV* inData, olong *nchan, olong *nIF, olong *npoln, 
		      ObitErr* err);
/* Write history */
void MapBeamHistory (gchar *Source, gchar Stok, ObitInfoList* myInput, 
		    ObitUV* inData, ObitImage* outImage, ObitErr* err);

/* Average  data */
ObitUV* doAvgData (ObitInfoList *myInput, ObitUV* inData, ObitErr *err);

/* Get list of antennas */
olong* getAntList (ObitInfoList* myInput, ObitErr* err);

/* Accumulate data into lists */
void  accumData (ObitUV* inData, ObitInfoList* myInput, olong ant,
		 olong nchan, olong nIF, olong selem, olong *nelem,
		 ofloat *SumIr, ofloat *SumIi, ofloat *SumII, ofloat *SumIWt,
		 ofloat *SumQr, ofloat *SumQi, ofloat *SumQQ, ofloat *SumQWt,
		 ofloat *SumUr, ofloat *SumUi, ofloat *SumUU, ofloat *SumUWt,
		 ofloat *SumVr, ofloat *SumVi, ofloat *SumVV, ofloat *SumVWt,
		 ofloat *SumAzCell, ofloat *SumElCell, ofloat *SumPACell, 
		 olong *CntCell, 
		 ofloat *avgAz, ofloat *avgEl, ofloat *avgPA, 
		 ObitErr* err);

/* Grid data into cells */
void  gridData (ObitInfoList* myInput, olong nchan, olong nIF, olong npoln,
		olong selem, olong nelem, gboolean doRMS, gboolean doPhase,
		ofloat *SumIr, ofloat *SumIi, ofloat *SumII, ofloat *SumIWt,
		ofloat *SumQr, ofloat *SumQi, ofloat *SumQQ, ofloat *SumQWt,
		ofloat *SumUr, ofloat *SumUi, ofloat *SumUU, ofloat *SumUWt,
		ofloat *SumVr, ofloat *SumVi, ofloat *SumVV, ofloat *SumVWt,
		ofloat *SumAzCell, ofloat *SumElCell, ofloat *SumPACell, 
		olong *CntCell, ObitFArray **grids);

/* Lagrangian interpolation coefficients */
void lagrange(ofloat x, ofloat y, olong n, olong hwid, 
	      ofloat *xlist, ofloat *ylist, ofloat *coef);

/* Program globals */
gchar *pgmName = "MapBeam";      /* Program name */
gchar *infile  = "MapBeam.in" ;  /* File with program inputs */
gchar *outfile = "MapBeam.out";  /* File to contain program outputs */
olong  pgmNumber;                /* Program number (like POPS no.) */
olong  AIPSuser;                 /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;                  /* Number of AIPS directories */
gchar **AIPSdirs=NULL;           /* List of AIPS data directories */
olong  nFITS=0;                  /* Number of FITS directories */
gchar **FITSdirs=NULL;           /* List of FITS data directories */
ObitInfoList *myInput  = NULL;   /* Input parameter list */
ObitInfoList *myOutput = NULL;   /* Output parameter list */
ofloat avgAz=0.0, avgEl=0.0, avgPA=0.0; /* Average observing Az, El, par Ang (deg) */
odouble RAMean=0.0, DecMean=0.0; /* Mean position of current source */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Obit task to make beam poln images from quasi holography data        */
/*----------------------------------------------------------------------- */
{
  oint         ierr      = 0;
  ObitSystem   *mySystem = NULL;
  ObitUV       *inData   = NULL;
  ObitErr      *err      = NULL;
 
   /* Startup - parse command line */
  err = newObitErr();
  myInput = MapBeamIn (argc, argv, err);
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

  /* Process */
  doSources (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);                /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* MapBeamIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "MapBeamIn";

  /* error checks */
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
      
     } else if (strcmp(arg, "-outDType") == 0) { /* Image type AIPS or FITS */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "DataType", OBIT_string, dim, strTemp);
      
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
      
     } else if (strcmp(arg, "-outFile") == 0) { /*outFile */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outFile", OBIT_string, dim, strTemp);

     } else if (strcmp(arg, "-AIPSdir") == 0) { /* Single AIPS Directory */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "AIPSdirs", OBIT_string, dim, strTemp);
      /* Only one AIPS Directory */
      dim[0] = 1;dim[1] = 1;
      itemp = 1; /* number of AIPS directories (1) */
      ObitInfoListPut (list, "nAIPS", OBIT_oint, dim, &itemp, err);

     } else if (strcmp(arg, "-FITSdir") == 0) { /* Single FITS Directory */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "FITSdirs", OBIT_string, dim, strTemp);
      /* Only one FITS Directory */
      dim[0] = 1;dim[1] = 1;
      itemp = 1; /* number of FITS directories (1) */
      ObitInfoListPut (list, "nFITS", OBIT_oint, dim, &itemp, err);

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
} /* end MapBeamIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program                                         */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: MapBeam -input file -output ofile [args]\n");
    fprintf(stderr, "MapBeam Obit task to make beam poln images\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def MapBeam.in\n");
    fprintf(stderr, "  -output output result file, def MapBeam.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -AIPSuser AIPS user number, def 2 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input \n");
    fprintf(stderr, "  -outDType AIPS or FITS type for output\n");
    fprintf(stderr, "  -inFile input FITS UV file\n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input image (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -BChan first channel to image\n");
    fprintf(stderr, "  -EChan highest channel to image\n");
    fprintf(stderr, "  -BIF first IF to image\n");
    fprintf(stderr, "  -EIF highest IF to image\n");
    fprintf(stderr, "  -outFile output image (FITS Image file\n");  
    fprintf(stderr, "  -outName output image (AIPS file name\n");
    fprintf(stderr, "  -outSeq output image (AIPS file sequence\n");
    fprintf(stderr, "  -outDisk output image ((AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -AIPSdir single AIPS data directory\n");
    fprintf(stderr, "  -FITSdir single FITS data directory\n");
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
/*     outSeq    Int        output AIPS image sequence no  [new]          */
/*     Sources   Str (16,1) Sources selected, blank = all                 */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*     doCalSelect Boo (1)  Apply calibration/selection?  def=True        */
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
  ofloat farray[3];
  gboolean btemp;
  olong  itemp, iarr[3];
  ObitInfoList *out = newObitInfoList();
  gchar *routine = "defaultInputs";

  /* error checks */
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
  strTemp = "MapBeam.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "MapBeamName";
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

  /* output FITS Image file name */
  strTemp = "MapBeamOut.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS Image file name */
  strTemp = "MapBeamOut";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS Image sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "outSeq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output AIPS or FITS Image disk number */
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

  /* Timerange in days */
  dim[0] = 2;dim[1] = 1;
  farray[0] = -1.0e20; farray[1] = 1.0e20;
  ObitInfoListPut (out, "timeRange", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Apply calibration/selection?, Always True */
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
  
  /* Reference antennas */
  dim[0] = 3; dim[1] = 1;
  iarr[0] = iarr[1] = iarr[2] = 0;
  ObitInfoListAlwaysPut (out, "RefAnts", OBIT_oint, dim, iarr);

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
  /* ObitInfoType type;
     gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1}; */
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

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
  ObitInfoType type;
  olong         Aseq, disk, cno, nvis;
  gchar        *Type, *strTemp, inFile[129];
  oint         doCalib;
  gchar        Aname[13], Aclass[7], *Atype = "UV";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doCalSelect;
  gchar        Stokes[5];
  gchar *routine = "getInputData";

  /* error checks */
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
    
    /* define object  */
    nvis = 1000;
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
    nvis = 1;
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
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);
 

  /* Convert to IQUV */
  strcpy (Stokes, "IQUV");
  dim[0] = strlen(Stokes);
  ObitInfoListAlwaysPut (inData->info, "Stokes", OBIT_string, dim, Stokes);

  /* Ensure inData fully instantiated and OK */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Set number of vis per IO */
  nvis = 1000;  /* How many vis per I/O? */
  nvis =  ObitUVDescSetNVis (inData->myDesc, myInput, nvis);
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO", OBIT_long, dim,  &nvis);

  return inData;
} /* end getInputData */


/*----------------------------------------------------------------------- */
/*  Create output image                                                   */
/*  One output image cube for requested Antenna/Stokes                    */
/*  output axes, Azimuth, Elevation, Channel, IF                          */
/*   Input:                                                               */
/*      Source    Source name                                             */
/*      iStoke    Stokes number (0-rel), I, Q, U, V                       */
/*      ant       Antenna number of image, 0=>all averaged                */
/*      doRMS     If TRUE, image is RMS                                   */
/*      doPhase   If TRUE, image is phase                                 */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV defining data                                    */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/*   Return  Output image depending on Stokes request                     */
/*----------------------------------------------------------------------- */
ObitImage* setOutput (gchar *Source, olong iStoke, olong ant, 
		      gboolean doRMS, gboolean doPhase, 
		      ObitInfoList *myInput, ObitUV* inData, 
		      ObitErr *err)
{
  ObitImage *outImage=NULL;
  ObitInfoType type;
  ObitIOType IOType;
  olong      i, n, Aseq, disk, cno, axNum, axFNum, axIFNum, nx=11, ny=11;
  gchar     *Type, *strTemp, outFile[129], *outName, *outF, stemp[32];;
  gchar     Aname[13], Aclass[7], *Atype = "MA";
  gchar     *today=NULL;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong     blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong     trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gboolean  exist;
  gchar     tname[129], *chStokes="IQUV";
  odouble   StokCrval[4] = {1.0,2.0,3.0,4.0};
  gfloat    xCells = 1.0, yCells=1.0;
  gchar     *FITS = "FITS";
  ObitImageDesc *outDesc;
  gchar     *routine = "setOutput";

  /* error checks */
  if (err->error) return outImage;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Notes:
     1) AIPS Name = Source +outName
            Class = "Santno"/"Sall" S=Stokes, specific or average ant.
        FITS Name = Source+"Santno"/"Sall"+outFile S=Stokes, specific or average an
        for RMS add "RM" to end of class or equivalent
        for Phase add "Ph" to end of class or equivalent
  */

  /* Create basic output Image Object */
  g_snprintf (tname, 100, "output Image %cPol",chStokes[iStoke-1]);
  outImage = newObitImage(tname);
    
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4)))
    ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS output */
    /* Generate output name from Source, outName */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&outName);
    /* Something in source name? */
    if ((Source[0]==' ') || (Source[0]==0)) g_snprintf (tname, 100, "%s", outName);
    else g_snprintf (tname, 100, "%s%s", Source, outName);
    /* If no name use input name */
    if ((tname[0]==' ') || (tname[0]==0)) {
      ObitInfoListGetP (myInput, "inName", &type, dim, (gpointer)&strTemp);
      g_snprintf (tname, 100, "%s", strTemp);
    }
      
    IOType = OBIT_IO_AIPS;  /* Save file type */
    /* input AIPS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);
    for (i=0; i<12; i++) Aname[i] = ' '; Aname[i] = 0;
    strncpy (Aname, tname, 13); 
    Aname[12] = 0;

    /* output AIPS class */
    /* Stokes or RMS? */
    if (doRMS) {  /* RMS */
      /* One or all antennas? */
      if (ant>0) {   /* One */
	g_snprintf (tname, 100, "%c%3.3dRM", chStokes[iStoke],ant);
      } else {       /* all */
	g_snprintf (tname, 100, "%cAllRM", chStokes[iStoke]);
      }
    }  else if (doPhase) {  /* Phase */
      /* One or all antennas? */
      if (ant>0) {   /* One */
	g_snprintf (tname, 100, "%c%3.3dPH", chStokes[iStoke],ant);
      } else {       /* all */
	g_snprintf (tname, 100, "%cAllPH", chStokes[iStoke]);
      }
    } else {      /* Stokes */
      if (ant>0) {   /* One */
	g_snprintf (tname, 100, "%c%3.3d  ", chStokes[iStoke],ant);
      } else {       /* all */
	g_snprintf (tname, 100, "%cAll  ", chStokes[iStoke]);
      }
    } /* end of branch by type */
    strncpy (Aclass, tname, 7); 
    Aclass[6] = 0;

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Find catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Output AIPS image %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);

    /* define object */
    ObitImageSetAIPS (outImage, OBIT_IO_byPlane, disk, cno, AIPSuser, 
		      blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);

  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */
    /* Generate output name from Source, outName */
    ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&outF);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) tname[i] = outF[i]; tname[i] = 0;
    /* If blank use ".fits" */
    if ((tname[0]==' ') || (tname[0]==0)) g_snprintf (tname, 128, ".fits");
    /* Something in source name? */
    if ((Source[0]==' ') || (Source[0]==0)) 
      g_snprintf (stemp, 30, "Beam");
    else g_snprintf (stemp, 30, Source);
    ObitTrimTrail(stemp);  /* remove trailing blanks */
	   
    IOType = OBIT_IO_FITS;  /* Save file type */

    /* Set output file name */
    /* Stokes or RMS? */
    if (doRMS) {  /* RMS */
      /* One or all antennas? */
      if (ant>0) {   /* One */
	g_snprintf (outFile, 128, "%s%c%3.3dRMS%s", stemp,chStokes[iStoke],ant,tname);
      } else {       /* all */
	g_snprintf (outFile, 128, "%s%cAllRMS%s", stemp,chStokes[iStoke],tname);
      }
     } else if (doPhase) {  /* Phase */
      /* One or all antennas? */
      if (ant>0) {   /* One */
	g_snprintf (outFile, 128, "%s%c%3.3dPh%s", stemp,chStokes[iStoke],ant,tname);
      } else {       /* all */
	g_snprintf (outFile, 128, "%s%cAllPh%s", stemp,chStokes[iStoke],tname);
      }
    } else {      /* Stokes */
      if (ant>0) {   /* One */
	g_snprintf (outFile, 128, "%s%c%3.3d%s", stemp,chStokes[iStoke],ant,tname);
      } else {       /* all */
	g_snprintf (outFile, 128, "%s%cAll%s", stemp,chStokes[iStoke],tname);
      }
    }  /* end branch by type */

    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    
    /* Give output Image name */
    Obit_log_error(err, OBIT_InfoErr, "Output FITS image %s on disk %d ",
		   outFile, disk);

    /* define object */
    ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, outFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outImage;
  } /* end of branch by output type */
  if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
  
  /* Size of image */
  ObitInfoListGetTest(myInput, "nx", &type, dim, &nx);
  ObitInfoListGetTest(myInput, "ny", &type, dim, &ny);
  ObitInfoListGetTest(myInput, "xCells", &type, dim, &xCells);
  ObitInfoListGetTest(myInput, "yCells", &type, dim, &yCells);

  /* Define header */
  /* Most info from UV Data descriptor */
  outDesc = outImage->myDesc;
  ObitImageUtilUV2ImageDesc (inData->myDesc, outDesc, TRUE, 1);

  /* Creation date today */
  today = ObitToday();
  strncpy (outDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);

  /* units */
  if (doPhase) {
    strcpy (outImage->myDesc->bunit, "DEGREE  ");
  }else {
    strcpy (outImage->myDesc->bunit, "JY/BEAM ");
  }

  /* Define axes */
  axNum = 0;

  /* Azimuth */
  strcpy (outImage->myDesc->ctype[axNum], "AZIMUTH");
  outImage->myDesc->inaxes[axNum] = nx;
  outImage->myDesc->crpix[axNum]  = nx/2 + 1.0;
  outImage->myDesc->cdelt[axNum]  = xCells/3600.;
  outImage->myDesc->crval[axNum]  = 0.0;
  axNum++;

  /* Elevation */
  strcpy (outImage->myDesc->ctype[axNum], "ELEVATIO");
  outImage->myDesc->inaxes[axNum] = ny;
  outImage->myDesc->crpix[axNum]  = ny/2 + 1.0;
  outImage->myDesc->cdelt[axNum]  = yCells/3600.;
  outImage->myDesc->crval[axNum]  = 0.0;
  axNum++;

  /* If there is only one channel and multiple IFs, IF axis is first - else Freq */
  if ((inData->myDesc->inaxes[inData->myDesc->jlocf]==1) && 
      ((inData->myDesc->jlocif>=0) && (inData->myDesc->inaxes[inData->myDesc->jlocif]>1))) {
    axFNum  = axNum+1;
    axIFNum = axNum;
  } else {
    axFNum  = axNum;
    axIFNum = axNum+1;
  }

  /* Channel */
  for (i=0; i<8; i++)
    outImage->myDesc->ctype[axFNum][i] = inData->myDesc->ctype[inData->myDesc->jlocf][i];
  outImage->myDesc->inaxes[axFNum] = inData->myDesc->inaxes[inData->myDesc->jlocf];
  outImage->myDesc->crpix[axFNum]  = inData->myDesc->crpix[inData->myDesc->jlocf];
  outImage->myDesc->cdelt[axFNum]  = inData->myDesc->cdelt[inData->myDesc->jlocf];
  outImage->myDesc->crval[axFNum]  = inData->myDesc->crval[inData->myDesc->jlocf];
  axNum++;


  /* IF */
  for (i=0; i<8; i++)
    outImage->myDesc->ctype[axIFNum][i] = inData->myDesc->ctype[inData->myDesc->jlocif][i];
  outImage->myDesc->inaxes[axIFNum] = inData->myDesc->inaxes[inData->myDesc->jlocif];
  outImage->myDesc->crpix[axIFNum]  = inData->myDesc->crpix[inData->myDesc->jlocif];
  outImage->myDesc->cdelt[axIFNum]  = inData->myDesc->cdelt[inData->myDesc->jlocif];
  outImage->myDesc->crval[axIFNum]  = inData->myDesc->crval[inData->myDesc->jlocif];
  axNum++;
  
   /* Stokes */
  for (i=0; i<8; i++)
    outImage->myDesc->ctype[axNum][i] = inData->myDesc->ctype[inData->myDesc->jlocs][i];
  outImage->myDesc->inaxes[axNum] = 1;
  outImage->myDesc->crpix[axNum]  = 1.0;
  outImage->myDesc->cdelt[axNum]  = 1.0;
  outImage->myDesc->crval[axNum]  = StokCrval[iStoke];
  axNum++;

  outImage->myDesc->naxis = axNum;  /* Total number of axes */

  /* reset image max/min */
  outDesc->maxval    = -1.0e20;
  outDesc->minval    =  1.0e20;

  /* Instantiate */
  ObitImageFullInstantiate (outImage, FALSE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);


  return outImage;
} /* end setOutput */

/*----------------------------------------------------------------------- */
/*  Loop over selected sources, these are all sources in the source table */
/*  with selection by souCode, Sources, timeRange etc.                    */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to image                                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doSources  (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  gchar        Source[17], lastSource[17];
  ObitSourceList* doList;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         maxlen, isource, failed=0, good=0;
  olong         iant, *antList=NULL, ant;
  gboolean     isBad = FALSE;
  gchar        *Fail="Failed  ", *Done="Done    ";
  gchar        *dataParms[] = {  /* Source selection*/
    "Sources", "souCode", "Qual", "timeRange", "UVRange", "FreqID",
    NULL
  };
  gchar *routine = "doSources";

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Make sure selector set on inData */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  
  /* Get source list to do */
  doList = ObitUVUtilWhichSources (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Get list of antennas */
  antList = getAntList(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Loop over list of sources */
  strncpy (lastSource, "None such       ", 16);
  for (isource = 0; isource<doList->number; isource++) {
    if (!doList->SUlist[isource]) continue; /* removed? */
    maxlen = MIN (16, strlen(doList->SUlist[isource]->SourceName));
    strncpy (Source, doList->SUlist[isource]->SourceName, maxlen);
    Source[maxlen] = 0;
    /* Ignore if same source name as last - just different qualifier */
    if (!strncmp(lastSource, Source, 16)) continue;
    strncpy (lastSource, Source, 16);
    /* Save position in global */
    RAMean  = doList->SUlist[isource]->RAMean;
    DecMean = doList->SUlist[isource]->DecMean;

    Obit_log_error(err, OBIT_InfoErr, " ******  Source %s ******", Source);
    ObitTrimTrail(Source);  /* remove trailing blanks */

    /* Save field name */
    dim[0] = 16; dim[1] = 1;
    ObitInfoListAlwaysPut (myInput, "FieldName", OBIT_string, dim, Source);

    /* Loop over antennas */
    iant = 0;
    while (antList[iant]>=0) {
      ant = antList[iant];
      iant++;

      /* Process source */
      doChanPoln (Source, ant, myInput, inData, err);
      /* Allow up to 10 failures before first success or up to 10% of large run */
      if (err->error) {
	ObitErrLog(err); /* Show failure messages */
	failed++;
	isBad = TRUE;
	if (((failed>=10) && (good<=0)) || 
	    (((failed>=10)&&(failed>0.1*doList->number)))) {
	  /* This isn't working - Give up */
	  Obit_log_error(err, OBIT_Error, "%s: Too many failures, giving up", 
			 routine);
	  return;
	}
      } else {
	isBad = FALSE;
	good++;
      } /* OK */
    } /* end antenna loop */


    /* Save processing summary - success or failure? */
    dim[0] = 8; dim[1] = 1;
    if (isBad) 
      ObitInfoListAlwaysPut (myInput, "Status", OBIT_string, dim, Fail);
    else
      ObitInfoListAlwaysPut (myInput, "Status", OBIT_string, dim, Done);
  } /* end source loop */

  doList = ObitSourceListUnref(doList);
  if (antList) g_free(antList);

}  /* end doSources */

/*----------------------------------------------------------------------- */
/*  Loop over frequencies and polarizations for a single source/antenna   */
/*   Input:                                                               */
/*      Source    Name of source being imaged                             */
/*      ant       Antenna number of image, 0=>all averaged                */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to image                                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doChanPoln (gchar *Source, olong ant, ObitInfoList* myInput, 
		 ObitUV* inData, ObitErr* err)
{
  ObitUV       *avgData = NULL;
  ObitImage    *outImage[4]={NULL,NULL,NULL,NULL};
  ObitFArray   **imgList = NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        i, loop, nloop, nchan, nIF, npoln, ichan, iIF, ipoln, plane[5];
  gboolean     doRMS=FALSE, doPhase=FALSE, isRMS[3] = {FALSE, TRUE, TRUE};
  gboolean     isPhase[3] = {FALSE, TRUE, TRUE};
  gchar        *chStokes=" IQUV";
  gchar        *tableList[]={"AIPS FQ","AIPS AN", NULL};
  gchar        *routine = "doChanPoln";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Get averaged data */
  avgData = doAvgData (myInput, inData, err);
  if (err->error) goto cleanup;

  /* Want RMS or phase images as well? */
  ObitInfoListGetTest(myInput, "doRMS",   &type, dim, &doRMS);
  ObitInfoListGetTest(myInput, "doPhase", &type, dim, &doPhase);
  nloop = 1;
  /* Set up loops for, ampl, ampl+phase, amp+rms and amp+phase+rms */
  if (doRMS)   nloop++;
  if (doPhase) nloop++;
  if (doRMS && doPhase) {isPhase[2] = isRMS[1] = FALSE; }

  /* Loop over poln, phase, RMS */
  for (loop=0; loop<nloop; loop++) {

    /* Image data in order (fastest to slowest, channel, IF, Stokes ) */
    imgList = doImage (isRMS[loop], isPhase[loop], ant, myInput, avgData, &nchan, &nIF, &npoln,  
		       err);
    if (err->error) goto cleanup;
    
    /* Write images */
    for (i=0; i<5; i++) plane[i] = 1;
    i = 0;
    for (ipoln=0; ipoln<npoln; ipoln++) {
      /* Create output for each poln */
      outImage[ipoln] = setOutput (Source, ipoln, ant, isRMS[loop], isPhase[loop], 
				   myInput, avgData, err);
      if (err->error) goto cleanup;
      ObitImageOpen (outImage[ipoln], OBIT_IO_WriteOnly, err);
      /* Save average geometry */
      dim[0] = dim[1] = dim[2] = 1;
      ObitInfoListAlwaysPut (outImage[ipoln]->myDesc->info, "avgAz", OBIT_float, dim, &avgAz);
      ObitInfoListAlwaysPut (outImage[ipoln]->myDesc->info, "avgEl", OBIT_float, dim, &avgEl);
      ObitInfoListAlwaysPut (outImage[ipoln]->myDesc->info, "avgPA", OBIT_float, dim, &avgPA);

      /* If there is only one channel and multiple IFs, IF axis is first - else Freq */
      if ((nchan==1) && (nIF>1)) { /* IF first */
	for (ichan=0; ichan<nchan; ichan++) {
	  plane[1] = ichan+1;
	  for (iIF=0; iIF<nIF; iIF++) {
	    plane[0] = iIF+1;
	    ObitImagePutPlane (outImage[ipoln], imgList[i++]->array, plane,  err);
	    if (err->error) goto cleanup;
	  }
	}
      } else { /* Freq first */
	for (iIF=0; iIF<nIF; iIF++) {
	  plane[1] = iIF+1;
	  for (ichan=0; ichan<nchan; ichan++) {
	    plane[0] = ichan+1;
	    ObitImagePutPlane (outImage[ipoln], imgList[i++]->array, plane,  err);
	    if (err->error) goto cleanup;
	  } /* end Channel loop */
	} /* end IF loop */
      }
      ObitImageClose (outImage[ipoln], err);
      if (err->error) goto cleanup;
      /* Copy FQ & AN tables */
      ObitDataCopyTables ((ObitData*)inData, (ObitData*)outImage[ipoln], NULL,
			  tableList, err);
    } /* end stokes loop */
    
    /* History */
    for (ipoln=0; ipoln<npoln; ipoln++) {
      MapBeamHistory (Source, chStokes[ipoln], myInput, inData, 
		      outImage[ipoln], err);
      if (err->error) goto cleanup;
    }  
    
    /* Cleanup */
    if (imgList) {
      i = 0;
      for (ipoln=0; ipoln<npoln; ipoln++) {
	outImage[ipoln] = ObitImageUnref(outImage[ipoln]);
	for (iIF=0; iIF<nIF; iIF++) {
	  for (ichan=0; ichan<nchan; ichan++) {
	    ObitFArrayUnref(imgList[i++]);
	  }
	}
      }
      g_free(imgList); imgList = NULL;
    } /* end cleanup image list */
  } /* end Stokes/RMS loop */
  
    /* Cleanup */
 cleanup:
  if (imgList) {
    i = 0;
    for (ipoln=0; ipoln<npoln; ipoln++) {
      outImage[ipoln] = ObitImageUnref(outImage[ipoln]);
      for (iIF=0; iIF<nIF; iIF++) {
	for (ichan=0; ichan<nchan; ichan++) {
	  ObitFArrayUnref(imgList[i++]);
	}
      }
    }
    g_free(imgList); imgList = NULL;
  } /* end cleanup image list */

  avgData  = ObitUVUnref(avgData); /* Scratch averaged data */
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
}  /* end doChanPoln */

/*----------------------------------------------------------------------- */
/*  Image all selected data                                               */
/*  Returns a list of ObitFArrays containing pixel data                   */
/*  fastest to slowest, channel, IF, Stokes                               */
/*   Input:                                                               */
/*      doRMS     If true make images of RMS, else Poln                   */
/*      doPhase   If TRUE, image is Phase, else Amplitude                 */
/*      ant       Antenna number of image, 0=>all                         */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to image                                         */
/*   Output:                                                              */
/*      nchan     Number of channels in output                            */
/*      nIF       Number of IFs in output                                 */
/*      npoln     Number of Stokes in output                              */
/*      err       Obit Error stack                                        */
/*   Output in globals                                                    */
/*      avgAz     Average Azimuth (deg)                                   */
/*      avgEl     Average Elevation (deg)                                 */
/*      avgPA     Average Parallactic angle (deg)                         */
/*   Return: List of FArrays, Unref/g_free when done                      */
/*----------------------------------------------------------------------- */
ObitFArray** doImage (gboolean doRMS, gboolean doPhase, olong ant, 
		      ObitInfoList* myInput, 
		      ObitUV* inData, olong *nchan, olong *nIF, olong *npoln, 
		      ObitErr* err)
{
  ObitFArray   **out = NULL;
  ObitFArray   *fatemp = NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        nx=11, ny=11, naxis[2], ipos[2], i, ipoln, ichan, iIF;
  olong        selem, nelem, size, indx, qndx, undx, vndx;
  ofloat       *SumIr=NULL, *SumIi=NULL, *SumII=NULL, *SumIWt=NULL;
  ofloat       *SumQr=NULL, *SumQi=NULL, *SumQQ=NULL, *SumQWt=NULL;
  ofloat       *SumUr=NULL, *SumUi=NULL, *SumUU=NULL, *SumUWt=NULL;
  ofloat       *SumVr=NULL, *SumVi=NULL, *SumVV=NULL, *SumVWt=NULL;
  ofloat       *SumAzCell=NULL, *SumElCell=NULL, *SumPACell=NULL;
  ofloat       *Center, *ICenter, value;
  olong        *CntCell=NULL;
  gchar *routine = "doImage";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /*  Number of channels */
  *nchan = inData->myDesc->inaxes[inData->myDesc->jlocf];
  /* Number of IFs */
  if (inData->myDesc->jlocif>=0)
    *nIF  = inData->myDesc->inaxes[inData->myDesc->jlocif];
  else *nIF = 1;
  /* Number of polarizations */
  *npoln = inData->myDesc->inaxes[inData->myDesc->jlocs];

  /* How big an image? */
  ObitInfoListGetTest(myInput, "nx", &type, dim, &nx);
  ObitInfoListGetTest(myInput, "ny", &type, dim, &ny);
  naxis[0] = nx; naxis[1] = ny;

  /* Make accumulation lists */
  selem = (*nIF) * (*nchan); /* size of list element */
  nelem = nx * ny;           /* Number of list elements */
  nelem *= 2;                /* To be sure big enough */
  size = selem * nelem;      /* size of list */
  SumIr     = g_malloc0(size*sizeof(ofloat));
  SumIi     = g_malloc0(size*sizeof(ofloat));
  SumII     = g_malloc0(size*sizeof(ofloat));
  SumIWt    = g_malloc0(size*sizeof(ofloat));
  SumQr     = g_malloc0(size*sizeof(ofloat));
  SumQi     = g_malloc0(size*sizeof(ofloat));
  SumQQ     = g_malloc0(size*sizeof(ofloat));
  SumQWt    = g_malloc0(size*sizeof(ofloat));
  SumUr     = g_malloc0(size*sizeof(ofloat));
  SumUi     = g_malloc0(size*sizeof(ofloat));
  SumUU     = g_malloc0(size*sizeof(ofloat));
  SumUWt    = g_malloc0(size*sizeof(ofloat));
  SumVr     = g_malloc0(size*sizeof(ofloat));
  SumVi     = g_malloc0(size*sizeof(ofloat));
  SumVWt    = g_malloc0(size*sizeof(ofloat));
  SumVV     = g_malloc0(size*sizeof(ofloat));
  SumAzCell = g_malloc0(nelem*sizeof(ofloat));  
  SumElCell = g_malloc0(nelem*sizeof(ofloat));
  SumPACell = g_malloc0(nelem*sizeof(ofloat));
  CntCell   = g_malloc0(nelem*sizeof(olong));

  /* Accumulate data onto lists */
  accumData (inData, myInput, ant, *nchan, *nIF, selem, &nelem,
	     SumIr, SumIi, SumII, SumIWt, SumQr, SumQi, SumQQ, SumQWt,
	     SumUr, SumUi, SumUU, SumUWt, SumVr, SumVi, SumVV, SumVWt,
	     SumAzCell, SumElCell, SumPACell, CntCell, 
	     &avgAz, &avgEl, &avgPA, err);
  if (err->error) goto cleanup;

  /* Create output */
  out = g_malloc0((*nchan)*(*nIF)*(*npoln)*sizeof(ObitFArray*));
  i = 0;
  for (ipoln=0; ipoln<(*npoln); ipoln++) {
    for (iIF=0; iIF<(*nIF); iIF++) {
      for (ichan=0; ichan<(*nchan); ichan++) {
	out[i++] = ObitFArrayCreate (NULL, 2, naxis);
      }
    }
  }

  /* Grid data */
  gridData (myInput, *nchan, *nIF, *npoln, selem, nelem, doRMS, doPhase,
	    SumIr, SumIi, SumII, SumIWt, SumQr, SumQi, SumQQ, SumQWt,
	    SumUr, SumUi, SumUU, SumUWt, SumVr, SumVi, SumVV, SumVWt,
	    SumAzCell, SumElCell, SumPACell, CntCell, 
	    out);

  /* No normalization for RMS, phase */
  if (doRMS|| doPhase) goto cleanup;

  /* Subtract center Q, U (source poln) * Ipol beam from rest */
  fatemp = ObitFArrayCreate (NULL, 2, naxis);  /* Work array */
  ipos[0] = nx/2; ipos[1] = ny/2;
  for (iIF=0; iIF<(*nIF); iIF++) {
    for (ichan=0; ichan<(*nchan); ichan++) {
      /* Q */ 
      indx = 0*selem + iIF*(*nchan) + ichan;
      qndx = 1*selem + iIF*(*nchan) + ichan;
      ICenter = ObitFArrayIndex (out[indx],  ipos);
      Center  = ObitFArrayIndex (out[qndx],  ipos);
      value = -(*Center) / (*ICenter);  /* Normalize IPol beam */
      fatemp = ObitFArrayCopy (out[indx], fatemp, err);
      if (err->error) goto cleanup;
      ObitFArraySMul (fatemp, value);  /* Q * IPol beam */
      ObitFArrayAdd (out[qndx], fatemp, out[qndx]);
      /* U */ 
      undx = 2*selem + iIF*(*nchan) + ichan;
      Center = ObitFArrayIndex (out[undx],  ipos);
      value = -(*Center) / (*ICenter);
      fatemp = ObitFArrayCopy (out[indx], fatemp, err);
      if (err->error) goto cleanup;
      ObitFArraySMul (fatemp, value);  /* U * IPol beam */
      ObitFArrayAdd (out[undx], fatemp, out[undx]);
    }
  }

  /* Divide Q,U,V by I  */
  for (iIF=0; iIF<(*nIF); iIF++) {
    for (ichan=0; ichan<(*nchan); ichan++) {
      indx = 0*selem + iIF*(*nchan) + ichan;
      qndx = 1*selem + iIF*(*nchan) + ichan;
      ObitFArrayDiv (out[qndx], out[indx], out[qndx]);
      undx = 2*selem + iIF*(*nchan) + ichan;
      ObitFArrayDiv (out[undx], out[indx], out[undx]);
      vndx = 3*selem + iIF*(*nchan) + ichan;
      ObitFArrayDiv (out[vndx], out[indx], out[vndx]);
    }
  }

  /* Normalize I */
  for (iIF=0; iIF<(*nIF); iIF++) {
    for (ichan=0; ichan<(*nchan); ichan++) {
      indx = 0*selem + iIF*(*nchan) + ichan;
      Center = ObitFArrayIndex (out[indx],  ipos);
      value = (*Center);
      if (value!=0.0) value = 1.0 / value;
      ObitFArraySMul (out[indx], value);
    }
  }

  /* Cleanup */
 cleanup:
  if (SumIr)     g_free(SumIr);
  if (SumIi)     g_free(SumIi);
  if (SumII)     g_free(SumII);
  if (SumIWt)    g_free(SumIWt);
  if (SumQr)     g_free(SumQr);
  if (SumQi)     g_free(SumQi);
  if (SumQQ)     g_free(SumQQ);
  if (SumQWt)    g_free(SumQWt);
  if (SumUr)     g_free(SumUr);
  if (SumUi)     g_free(SumUi);
  if (SumUU)     g_free(SumUU);
  if (SumUWt)    g_free(SumUWt);
  if (SumVr)     g_free(SumVr);
  if (SumVi)     g_free(SumVi);
  if (SumVV)     g_free(SumVV);
  if (SumVWt)    g_free(SumVWt);
  if (SumAzCell) g_free(SumAzCell);
  if (SumElCell) g_free(SumElCell);
  if (CntCell)   g_free(CntCell);
  fatemp = ObitFArrayUnref(fatemp);
  if (err->error) Obit_traceback_val (err, routine, inData->name, out);

  return out;
} /* end doImage */

/*----------------------------------------------------------------------- */
/*  Write History for MapBeam                                             */
/*   Input:                                                               */
/*      Source    Name of source being imaged                             */
/*      Stoke     Stokes's parameter imaged I, Q, U, V                    */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outImage  ObitImage to write history to                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void MapBeamHistory (gchar *Source, gchar Stoke, ObitInfoList* myInput, 
		    ObitUV* inData, ObitImage* outImage, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "Sources", "Qual", "souCode", "timeRange",  
    "FreqID", "BIF", "EIF", "BChan", "EChan",  
    "doCalib", "gainUse", "doPol",  "PDVer", "flagVer", "doBand ",  "BPVer", "Smooth",
    "outDType", "outFile",  "outDisk", "outName", "outSeq",
    "nx", "ny", "xCells", "yCells", "hwid", "avgTime", "avgFreq", "chAvg", "ChanSel",
    "blnkTime", "avgAnt", "doRMS", "doPhase", "doPolNorm", "Antennas", "RefAnts",
    NULL};
  gchar *routine = "MapBeamHistory";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));
  g_assert (ObitImageIsA(outImage));

  /* Do Image history  */
  inHistory  = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory ((ObitData*)outImage, OBIT_IO_WriteOnly, err);

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
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Write source and poln */
  if (Stoke==' ') Stoke = 'I';
  g_snprintf (hicard, 80, "%s Source = '%s', Stokes= '%c'", pgmName, Source, Stoke);
  ObitHistoryWriteRec (outHistory, -1, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Average geometry */
  g_snprintf (hicard, 80, "%s / Average observing Azimuth %f deg", pgmName, avgAz);
  ObitHistoryWriteRec (outHistory, -1, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  g_snprintf (hicard, 80, "%s / Average observing Elevation %f deg", pgmName, avgEl);
  ObitHistoryWriteRec (outHistory, -1, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  g_snprintf (hicard, 80, "%s / Average observing Para. Angle %f deg", pgmName, avgPA);
  ObitHistoryWriteRec (outHistory, -1, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end MapBeamHistory  */

/*----------------------------------------------------------------------- */
/*   Average data as requested                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList, uses:                     */
/*      inData    ObitUV to average from                                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return:                                                              */
/*     Averaged scratch data file, Unref when done                        */
/*----------------------------------------------------------------------- */
ObitUV* doAvgData (ObitInfoList *myInput, ObitUV* inData, ObitErr *err)
{
  ObitUV       *avgData=NULL;
  ObitUV       *scrData=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM]= {1,1,1,1,1};
  olong        avgFreq, nchAvg;
  ofloat       timeAvg;
  gboolean     isScratch, doAvgAll;
  gchar        Stokes[5];
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "UVRange", "timeRange", "doCalSelect", 
    "BIF", "EIF", "BChan", "EChan","subA", "Antennas",
    "doCalib", "gainUse", "doBand", "BPVer", "Smooth", "flagVer", 
    "doPol", "PDVer",  "avgTime", "avgFreq", "ChanSel", 
    NULL
  };
  gchar        *routine= "doAvgData";

  if (err->error) return avgData;  /* Prior error? */

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, inData->name, avgData);
  
  /* Convert to IQUV */
  strcpy (Stokes, "IQUV");
  dim[0] = strlen(Stokes);
  ObitInfoListAlwaysPut (inData->info, "Stokes", OBIT_string, dim, Stokes);

  /* Make scratch file, possibly with time and freq averaging */
  avgFreq = 0;
  ObitInfoListGetTest(myInput, "avgFreq",  &type, dim, &avgFreq);
  nchAvg = 1;
  ObitInfoListGetTest(myInput, "chAvg",  &type, dim, &nchAvg);
  timeAvg = 0.0;
  ObitInfoListGetTest(myInput, "avgTime",  &type, dim, &timeAvg);
  timeAvg /= 60.0;  /* Convert to min. */

  /* Average all channels/IFs? */
  doAvgAll = (avgFreq==3);

  /* If both temporal and frequency averaging, frequency average to scratch */
  if ((avgFreq>0) && (timeAvg>0.0)) {
    /* First frequency */
    dim[0] = dim[1] = 1;
    ObitInfoListAlwaysPut (inData->info, "NumChAvg", OBIT_long, dim, &nchAvg);
    dim[0] = dim[1] = 1;
    ObitInfoListAlwaysPut (inData->info, "doAvgAll", OBIT_bool, dim, &doAvgAll);
    isScratch = TRUE;
    scrData = ObitUVUtilAvgF (inData, isScratch, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, avgData);
    /* Then time */
    dim[0] = 1;
    ObitInfoListAlwaysPut (scrData->info, "timeAvg", OBIT_float, dim, &timeAvg);
    isScratch = TRUE;
    avgData = ObitUVUtilAvgT (scrData, isScratch, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, avgData);
    scrData = ObitUVUnref(scrData);

  } else if (avgFreq>0) {    /* Freq averaging only */
    dim[0] = dim[1] = 1;
    ObitInfoListAlwaysPut (inData->info, "NumChAvg", OBIT_long, dim, &nchAvg);
    dim[0] = dim[1] = 1;
    ObitInfoListAlwaysPut (inData->info, "doAvgAll", OBIT_bool, dim, &doAvgAll);
    isScratch = TRUE;
    avgData = ObitUVUtilAvgF (inData, isScratch, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, avgData);

  } else if (timeAvg>0.0) {  /* Time averaging only */
    dim[0] = 1;
    ObitInfoListAlwaysPut (inData->info, "timeAvg", OBIT_float, dim, &timeAvg);
    isScratch = TRUE;
    avgData = ObitUVUtilAvgT (inData, isScratch, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, avgData);

  } else { /* No averaging - straight copy */
    /* Scratch file for copy */
    avgData = newObitUVScratch (inData, err);
    /* Calibrate/edit/copy data to scratch file */
    avgData = ObitUVCopy (inData, avgData, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, avgData);
  }

  return avgData;
} /* end avgData */

/*----------------------------------------------------------------------- */
/*  Get list of antennas                                                  */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList, uses:                     */
/*        avgAnt    If TRUE average all antennas                          */
/*        Antennas  Desired antennas                                      */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return:                                                              */
/*     -1 terminated list, 0=> average antennas, g_free when done         */
/*----------------------------------------------------------------------- */
olong* getAntList (ObitInfoList* myInput, ObitErr* err)
{
  olong        *out=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM];
  gboolean     avgAnt=TRUE;
  olong        i, count, *Antennas;
  gchar        *routine= "getAntList";

  if (err->error) return out;  /* Prior error? */

  ObitInfoListGetTest(myInput, "avgAnt", &type, dim, &avgAnt);
  if (avgAnt) { /* Average all selected */
    out = g_malloc0(2*sizeof(olong));
    out[0] = 0;
    out[1] = -1;
    return out;
  } else { /* individual antenna output */
    /* Antenna list */
    ObitInfoListGetP(myInput, "Antennas",  &type, dim, (gpointer)&Antennas);
    if (Antennas==NULL) {  /* Not given */
      Obit_log_error(err, OBIT_Error, "%s: Antenna list not given", routine);
      return out;
    }
    /* Count */
    count = 0;
    for (i=0; i<dim[0]; i++) {
      if (Antennas[i]==0) break;
      count++;
    }
    /* If none selected all wanted */
    if (count<=0) {
      out = g_malloc0((dim[0]+1)*sizeof(olong));
      for (i=0; i<dim[0]; i++) out[i] = i+1;
      out[i] = -1;
    } else { /* explicit list */
      out = g_malloc0((count+1)*sizeof(olong));
      for (i=0; i<count; i++) out[i] = Antennas[i];
      out[i] = -1;
    }
  } /* end individual antenna output */

  return out;
} /* end getAntList */

/*----------------------------------------------------------------------- */
/*  Accumulate data into lists                                            */
/*  Accumulate real part of correlations for all Stokes/freq/IF           */
/*  Sums are in a list with each set of entries corresponding to a given  */
/*  pointing.                                                             */
/*  Q and U corrected for parallactic angle                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to image                                         */
/*      ant       Antenna number of image, 0=>all                         */
/*      nchan     Number of channels in output                            */
/*      nIF       Number of IFs in output                                 */
/*      selem     Size (floats) of list element                           */
/*   In/Output:                                                           */
/*      nelem     (in) Max number of list elements                        */
/*                (out) Actual number of list elements                    */
/*   Output:                                                              */
/*      SumIr     Real Stokes I accumulation list                         */
/*      SumIi     Imag Stokes I accumulation list                         */
/*      SumII     Stokes I*I accumulation list                            */
/*      SumIWt    I Weight accumulation list                              */
/*      SumQr     Real Stokes Q accumulation list                         */
/*      SumQi     Imag Stokes Q accumulation list                         */
/*      SumQQ     Stokes Q*Q accumulation list                            */
/*      SumQWt    Q Weight accumulation list                              */
/*      SumUr     Real Stokes U accumulation list                         */
/*      SumUi     Imag Stokes U accumulation list                         */
/*      SumUU     Stokes U*U accumulation list                            */
/*      SumUWt    U Weight accumulation list                              */
/*      SumVr     Real Stokes V accumulation list                         */
/*      SumVi     Imag Stokes V accumulation list                         */
/*      SumVV     Stokes V*V accumulation list                            */
/*      SumVWt    V Weight accumulation list                              */
/*      SumAzCell Azimuth offset accumulation list                        */
/*      SumElCell Elevation offset accumulation list                      */
/*      SumPACell Parallactic angle accumulation list                     */
/*      CntCell   Counts in geometry accumulation                         */
/*      avgAz     Average Azimuth (deg)                                   */
/*      avgEl     Average Elevation (deg)                                 */
/*      avgPA     Average Parallactic angle (deg)                         */
/*      err       Obit Error stack                                        */
/*----------------------------------------------------------------------- */
void  accumData (ObitUV* inData, ObitInfoList* myInput, olong ant,
		 olong nchan, olong nIF, olong selem, olong *nelem,
		 ofloat *SumIr, ofloat *SumIi, ofloat *SumII, ofloat *SumIWt,
		 ofloat *SumQr, ofloat *SumQi, ofloat *SumQQ, ofloat *SumQWt,
		 ofloat *SumUr, ofloat *SumUi, ofloat *SumUU, ofloat *SumUWt,
		 ofloat *SumVr, ofloat *SumVi, ofloat *SumVV, ofloat *SumVWt,
		 ofloat *SumAzCell, ofloat *SumElCell, ofloat *SumPACell, 
		 olong *CntCell, 
		 ofloat *avgAz, ofloat *avgEl, ofloat *avgPA, 
		 ObitErr* err)
{
  ObitAntennaList *AList=NULL;
  ObitSource *Source=NULL;
  ObitTableAN *ANTable=NULL;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitInfoType type;
  gint32   dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat   xCells=1.0, yCells=1.0, blnkTime=0.0, tblank=-1.0e20;
  ofloat   u, v, time, base, ulast, vlast, tol, Az, El, PA;
  ofloat   ss, cc, xr, xi, fblank =  ObitMagicF();
  odouble  sumAz, sumEl, sumPA;
  olong    count, maxElem=*nelem, iElem, indx, iant, ant1, ant2, off=0, iver;
  olong    i, j, jlocs, jlocf, jlocif, incs, incf, incif, doff, ddoff;
  olong    nx, ny, iIF, ichan, *refAnts, nRefAnt, ix, iy, prtLv=0;
  gboolean OK1, OK2;
  gchar    *routine = "accumData";

  /* error checks */
  if (err->error) return;

  /* Get control parameters */
  ObitInfoListGetTest(myInput, "xCells", &type,   dim, &xCells);
  ObitInfoListGetTest(myInput, "yCells", &type,   dim, &yCells);
  ObitInfoListGetTest(myInput, "blnkTime", &type, dim, &blnkTime);
  ObitInfoListGetTest(myInput, "prtLv",  &type, dim, &prtLv);
  /* How big an image? */
  ObitInfoListGetTest(myInput, "nx", &type, dim, &nx);
  ObitInfoListGetTest(myInput, "ny", &type, dim, &ny);
  /* Cell spacing to radians */
  xCells = (xCells / 3600.0) * DG2RAD;
  yCells = (yCells / 3600.0) * DG2RAD;
  /* Blanking time to days */
  blnkTime /= 86400.0;

  /* Ref antennas */
  nRefAnt = 0;
  ObitInfoListGetP(myInput, "RefAnts",  &type, dim, (gpointer)&refAnts);
  for (j=0; j<dim[0]; j++) {
    if (refAnts[j]<=0) break;
    nRefAnt++;
  }
  /* Check that a ref antenna given */
  Obit_return_if_fail ((nRefAnt>0), err, 
		       "%s NO reference antennas given",  routine);  

  /* Initialize */
  *avgAz = 0.0;
  *avgEl = 0.0;
  *avgPA = 0.0;
  sumAz  = sumEl = sumPA = 0.0;
  count  = 0;
  iElem  = -1;
  ulast  = vlast = 1.0e20;
  tol    = 0.3 * xCells;  /* Tolerance 0.3 cells  */
  iant = MAX (1, ant);

  /* Get Antenna List */
  iver = 1;  /* Or Subarray no. */
  ANTable = newObitTableANValue (inData->name, (ObitData*)inData, &iver, 
				 OBIT_IO_ReadOnly, 0, 0, 0, err);
  AList   = ObitTableANGetList (ANTable, err);
  ANTable = ObitTableANUnref(ANTable);   /* Done with table */
  if (err->error) Obit_traceback_msg (err, routine, ANTable->name);

  /* Source - get position from global */
  Source = newObitSource("Temp Source");
  Source->equinox = inData->myDesc->equinox;
  Source->RAMean  = RAMean;
  Source->DecMean = DecMean;
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  /*  Source->RAMean  = inData->myDesc->crval[inData->myDesc->jlocr];
      Source->DecMean = inData->myDesc->crval[inData->myDesc->jlocd];*/
  /* Compute apparent position */
  ObitPrecessUVJPrecessApp (inData->myDesc, Source);
  
  jlocs = inData->myDesc->jlocs;
  incs  = inData->myDesc->incs;
  jlocf = inData->myDesc->jlocf;
  incf  = inData->myDesc->incf;
  if (inData->myDesc->jlocif>=0) {
    jlocif = inData->myDesc->jlocif;
    incif  = inData->myDesc->incif;
  } else {
    jlocif = 0;
    incif  = 0;
  }

  /* Open uv data if not already open */
  if (inData->myStatus==OBIT_Inactive) {
    retCode = ObitUVOpen (inData, OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
  }

  /* Loop through data */
  while (retCode==OBIT_IO_OK) {
    /* read buffer full */
    retCode = ObitUVRead (inData, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    indx = 0;

    /* First time? */
    if (iElem<0) {
      iElem  = 0;
      tblank = inData->buffer[indx+inData->myDesc->iloct] + blnkTime;
      ulast  = inData->buffer[indx+inData->myDesc->ilocu];
      vlast  = inData->buffer[indx+inData->myDesc->ilocv];
    }

    /* loop over buffer */
    for (i=0; i<inData->myDesc->numVisBuff; i++) { 
      /* where are we? */
      time =  inData->buffer[indx+inData->myDesc->iloct];
      /* fix online system bug - flip sign of u */
      u    = -inData->buffer[indx+inData->myDesc->ilocu];
      v    =  inData->buffer[indx+inData->myDesc->ilocv];

      /* In blanked time? */
      if (time<tblank) goto next;

      /* Want antennas? */
      base = inData->buffer[indx+inData->myDesc->ilocb];
      /* crack Baseline */
      ant1 = (base / 256.0) + 0.001;
      ant2 = (base - ant1 * 256) + 0.001;
      if (ant>0) {
	if ((ant1!=ant) && (ant2!=ant)) goto next;
      }
      /* One and only one must be a reference antenna */
      OK1 = FALSE;
      OK2 = FALSE;
      for (j=0; j<nRefAnt; j++) {
	OK1 = OK1 || (ant1==refAnts[j]);
	OK2 = OK2 || (ant2==refAnts[j]);
      }
      if (!(OK1 || OK2)) goto next;
      if (OK1 && OK2)    goto next;

      /* First time in pointing? */
      if (CntCell[iElem]<=0) {
 	ulast  = u;
	vlast  = v;
     }

      /* New pointing? */
      if ((fabs(u-ulast)>tol) || (fabs(v-vlast)>tol)) {
	ulast  = u;
	vlast  = v;
	tblank = time + blnkTime;
	iElem++;
	/* Check if in bounds */
	Obit_return_if_fail ((iElem<maxElem), err, 
			     "%s Too many pointings %d >= %d",  
			     routine, iElem, maxElem);  
      }

      /* Observing geometry (radians) */
      Az = ObitAntennaListAz (AList, ant1, time, Source);
      El = ObitAntennaListElev (AList, ant1, time, Source);
      PA = ObitAntennaListParAng (AList, ant1, time, Source);

      /* Accumulate */
      sumAz += Az;
      sumEl += El;
      sumPA += PA;
      count++;

      SumAzCell[iElem] += u/xCells;
      SumElCell[iElem] += v/yCells;
      SumPACell[iElem] += PA;
      CntCell[iElem]++;
      /* Loop over freq, IF */
      for (iIF=0; iIF<nIF; iIF++) {
	for (ichan=0; ichan<nchan; ichan++) {
	  off  = iElem*selem + iIF*nchan + ichan;
	  doff = indx + inData->myDesc->nrparm + iIF*incif + ichan*incf;
	  ddoff = doff;
	  /* DEBUG 
	  if ((iIF==1) && (fabs((u/xCells)+5.0)<2.1) && (fabs((v/yCells)-2.0)<2.1)) {
	    fprintf (stderr,"cx %f cy %f wt %f\n",u/xCells-10.0, v/yCells-10.0, inData->buffer[ddoff+2]);
	  }*/
	  /* I */
	  if (inData->buffer[ddoff+2]>0.0) {
	    SumIr[off]  += inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumIi[off]  += inData->buffer[ddoff+1]*inData->buffer[ddoff+2];
	    SumII[off]  += inData->buffer[ddoff]*inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumIWt[off] += inData->buffer[ddoff+2]; 
	  }
	  /* Q */
	  ddoff = doff + incs;
	  if (inData->buffer[ddoff+2]>0.0) {
	    SumQr[off]  += inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumQi[off]  += inData->buffer[ddoff+1]*inData->buffer[ddoff+2];
	    SumQQ[off]  += inData->buffer[ddoff]*inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumQWt[off] += inData->buffer[ddoff+2]; 
	  }
	  /* U */
	  ddoff = doff + 2*incs;
	  if (inData->buffer[ddoff+2]>0.0) {
	    SumUr[off]  += inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumUi[off]  += inData->buffer[ddoff+1]*inData->buffer[ddoff+2];
	    SumUU[off]  += inData->buffer[ddoff]*inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumUWt[off] += inData->buffer[ddoff+2]; 
	  }
	  /* V */
	  ddoff = doff + 3*incs;
	  if (inData->buffer[ddoff+2]>0.0) {
	    SumVr[off]   += inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumVi[off]   += inData->buffer[ddoff+1]*inData->buffer[ddoff+2];
	    SumVV[off]   += inData->buffer[ddoff]*inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumVWt[off]  += inData->buffer[ddoff+2]; 
	  }
	}
      }

      /* update data pointers */
    next:
      indx  += inData->myDesc->lrec;
    } /* end loop over buffer */
  } /* end loop over file */
  
    /* Close */
  retCode = ObitUVClose (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Better have something */
  Obit_return_if_fail ((iElem>=2), err, 
		       "%s: Insufficient data included or bad cell spacing",  
		       routine);  

  /* Number of actual elements */
  *nelem = iElem+1;

  /* Normalize things */
  *avgAz = (sumAz/count)*RAD2DG;
  *avgEl = (sumEl/count)*RAD2DG;
  *avgPA = (sumPA/count)*RAD2DG;

  for (i=0; i<(*nelem); i++) {  /* Loop over elements */
    if (CntCell[i]>0) {
      SumAzCell[i] /= CntCell[i];
      SumElCell[i] /= CntCell[i];
      SumPACell[i] /= CntCell[i];
    } else {  /* Big value so they will be ignored */
      SumAzCell[i] = 1.0e10;
      SumElCell[i] = 1.0e10;
      SumPACell[i] = 1.0e10;
    }

    /* Loop over freq, IF */
    for (iIF=0; iIF<nIF; iIF++) {
      for (ichan=0; ichan<nchan; ichan++) {
	off  = i*selem + iIF*nchan + ichan;
	if (SumIWt[off]>0) {
	  SumIr[off] /= SumIWt[off];
	  SumIi[off] /= SumIWt[off];
	  SumII[off] = (((SumII[off]/SumIWt[off]) - SumIr[off]*SumIr[off]))/SumIWt[off];
	  if (SumII[off]>0.0) SumII[off] = sqrt(SumII[off]);
	  SumII[off] /= SumIr[off];  /* Normalize variance by I */
	} else {
	  SumIr[off] = fblank;
	  SumIi[off] = fblank;
	  SumII[off] = fblank;
	}
	if (SumQWt[off]>0) {
	  SumQr[off] /= SumQWt[off];
	  SumQi[off] /= SumQWt[off];
	  SumQQ[off] = (((SumQQ[off]/SumQWt[off]) - SumQr[off]*SumQr[off]))/SumQWt[off];
	  if (SumQQ[off]>0.0) SumQQ[off] = sqrt(SumQQ[off]);
	  SumQQ[off] /= SumIr[off];  /* Normalize variance by I */
	} else {
	  SumQr[off]  = fblank;
	  SumQi[off]  = fblank;
	  SumQQ[off] = fblank;
	}
	if (SumUWt[off]>0) {
	  SumUr[off] /= SumUWt[off];
	  SumUi[off] /= SumUWt[off];
	  SumUU[off] = (((SumUU[off]/SumUWt[off]) - SumUr[off]*SumUr[off]))/SumUWt[off];
	  if (SumUU[off]>0.0) SumUU[off] = sqrt(SumUU[off]);
	  SumUU[off] /= SumIr[off];  /* Normalize variance by I */
	} else {
	  SumUr[off] = fblank;
	  SumUi[off] = fblank;
	  SumUU[off] = fblank;
	}
	if (SumVWt[off]>0) {
	  SumVr[off] /= SumVWt[off];
	  SumVi[off] /= SumVWt[off];
	  SumVV[off] = (((SumVV[off]/SumVWt[off]) - SumVr[off]*SumVr[off]))/SumVWt[off];
	  if (SumVV[off]>0.0) SumVV[off] = sqrt(SumVV[off]);
	  SumVV[off] /= SumIr[off];  /* Normalize variance by I */
	} else {
	  SumVr[off] = fblank;
	  SumVi[off] = fblank;
	  SumVV[off] = fblank;
	}
	/* Counter rotate (Q+iU) for parallactic angle */
	if ((SumQr[off]!=fblank) && (SumUr[off]!=fblank)) {
	  cc =  cos(2.0*SumPACell[i]);
	  ss = -sin(2.0*SumPACell[i]);
	  /* Not sure this works */
	  xr = SumQr[off];
	  xi = SumUr[off];
	  SumQr[off] = cc*xr - ss*xi;
	  SumUr[off] = cc*xi + ss*xr;
	  xr = SumQi[off];
	  xi = SumUi[off];
	  SumQi[off] = cc*xr - ss*xi;
	  SumUi[off] = cc*xi + ss*xr;
	} /* end counter rotate */
      } /* end channel loop */
    } /* end IF loop */


    /* Add diagnostics */
    if (prtLv>=2) {
      ix = (olong) (SumAzCell[i] + nx/2 + 1.5);
      iy = (olong) (SumElCell[i] + ny/2 + 1.5);
      if ((SumAzCell[i]>1000.) || (SumElCell[i]>1000.)) continue;
      Obit_log_error(err, OBIT_InfoErr, 
		     "%3.3d Cell %3d %3d Az %8.1f cell, El %8.1f cell, I %6.3f %6.3f Q %6.3f %6.3f U %6.3f %6.3f V %6.3f %6.3f Jy",
		     i, ix,iy, 
		     /*SumAzCell[i]*xCells*206265., SumElCell[i]*yCells*206265., offset in asec */
		     SumAzCell[i], SumElCell[i],   /* offset in cells */
		     SumIr[i*selem],SumIi[i*selem], SumQr[i*selem],SumQi[i*selem],
		     SumUr[i*selem],SumUi[i*selem], SumVr[i*selem],SumVi[i*selem]);
    }
  } /* End loop normalizing list */

} /* end accumData  */

/*----------------------------------------------------------------------- */
/*  Interpolate quasi regular lists of points onto a regular grid         */
/*  Accumulate real part of correlations for all Stokes/freq/IF           */
/*  Sums are in a list with each set of entries corresponding to a given  */
/*  pointing.                                                             */
/*  Adapted from AIPS/MAPBM.FOR                                           */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      nchan     Number of channels in output                            */
/*      nIF       Number of IFs in output                                 */
/*      selem     Size (floats) of list element                           */
/*      nelem     Number of list elements                                 */
/*      doRMS     If TRUE, image is RMS                                   */
/*      doPhase   If TRUE, image is Phase, else Amplitude                 */
/*      SumIr     Real Stokes I accumulation list                         */
/*      SumIi     Imag Stokes I accumulation list                         */
/*      SumII     Stokes I*I accumulation list                            */
/*      SumIWt    I Weight accumulation list                              */
/*      SumQr     Real Stokes Q accumulation list                         */
/*      SumQi     Imag Stokes Q accumulation list                         */
/*      SumQQ     Stokes Q*Q accumulation list                            */
/*      SumQWt    Q Weight accumulation list                              */
/*      SumUr     Real Stokes U accumulation list                         */
/*      SumUi     Imag Stokes U accumulation list                         */
/*      SumUU     Stokes U*U accumulation list                            */
/*      SumUWt    U Weight accumulation list                              */
/*      SumVr     Real Stokes V accumulation list                         */
/*      SumVi     Imag Stokes V accumulation list                         */
/*      SumVV     Stokes V*V accumulation list                            */
/*      SumVWt    V Weight accumulation list                              */
/*      SumAzCell Azimuth offset accumulation list                        */
/*      SumElCell Elevation offset  accumulation list                     */
/*      SumPACell Parallactic angle accumulation list                     */
/*                Contents destroyed                                      */
/*      CntCell   Counts in geometry accumulation                         */
/*   Output:                                                              */
/*      grids     Array of ObitFArrays                                    */
/*                fastest to slowest, channel, IF, Stokes                 */
/*      err       Obit Error stack                                        */
/*----------------------------------------------------------------------- */
void  gridData (ObitInfoList* myInput, olong nchan, olong nIF, olong npoln,
		olong selem, olong nelem, gboolean doRMS, gboolean doPhase,
		ofloat *SumIr, ofloat *SumIi, ofloat *SumII, ofloat *SumIWt,
		ofloat *SumQr, ofloat *SumQi, ofloat *SumQQ, ofloat *SumQWt,
		ofloat *SumUr, ofloat *SumUi, ofloat *SumUU, ofloat *SumUWt,
		ofloat *SumVr, ofloat *SumVi, ofloat *SumVV, ofloat *SumVWt,
		ofloat *SumAzCell, ofloat *SumElCell, ofloat *SumPACell, 
		olong *CntCell, ObitFArray **grids)
{
  ofloat x, y, xcen, ycen, closest, xCells=1.0, yCells=1.0;
  ofloat *coef, amp, ph, fblank = ObitMagicF();
  odouble sumIWt , sumQWt, sumUWt, sumVWt;
  odouble valIr,  valIi, valII, valQr, valQi, valQQ;
  odouble valUr,  valUi, valUU, valVr, valVi, valVV;
  olong   i, iIF, ichan, nx=1, ny=1, hwid=1, ix, iy, indx, jndx, off;
  ObitFArray *array;
  ObitInfoType type;
  gint32   dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  /* Use SumPACell as work space */
  coef = SumPACell;

  /* How big an image? */
  ObitInfoListGetTest(myInput, "nx", &type, dim, &nx);
  ObitInfoListGetTest(myInput, "ny", &type, dim, &ny);
  ObitInfoListGetTest(myInput, "xCells", &type, dim, &xCells);
  ObitInfoListGetTest(myInput, "yCells", &type, dim, &yCells);
  ObitInfoListGetTest(myInput, "hwid",   &type, dim, &hwid);
  xcen = (ofloat)(nx/2);
  ycen = (ofloat)(ny/2);

  /* Loop over y (El) */
  for (iy=0; iy<ny; iy++) {
      y = iy - ycen;

    /* Loop over x (Az) */
    for (ix=0; ix<nx; ix++) {
      x = ix - xcen;
      
      /* Get interpolation coefficients */
      lagrange (x, y, nelem, hwid, SumAzCell, SumElCell, coef);

      /* Loop over IFs */
      for (iIF=0; iIF<nIF; iIF++) {
	/* Loop over chans */
	for (ichan=0; ichan<nchan; ichan++) {
	  valIr = valII = valQr = valQQ = valUr = valUU = valVr = valVV = 0.0;
	  valIi = valQi = valUi = valVi = 0.0;
	  sumIWt = sumQWt = sumUWt = sumVWt = 0.0;
	  off = ichan + iIF*nchan;
	  closest = 1000.0; /* Closest element */
	  /* Loop over lists summing */
	  for (i=0; i<nelem; i++) {
	    closest = MIN (closest, MAX (fabs(SumAzCell[i]-x), fabs(SumElCell[i]-y)));
	    if (coef[i]!=0.0) {
	      /* DEBUG */
	      if ((iIF==1) && (abs(ix-21)<=2) && (abs(iy-28)<=2)) {
		fprintf (stderr,"i %d x %d y %d flx %f coef %f sum %f\n",i, ix, iy, SumIr[i*selem+off], coef[i],sumIWt);
	      }
	      if (SumIr[i*selem+off]!=fblank) {
		valIr  += coef[i]*SumIr[i*selem+off];
		valIi  += coef[i]*SumIi[i*selem+off];
		valII  += coef[i]*SumII[i*selem+off];
		sumIWt += coef[i];
	      }
	      if (SumQr[i*selem+off]!=fblank) {
		valQr  += coef[i]*SumQr[i*selem+off];
		valQi  += coef[i]*SumQi[i*selem+off];
		valQQ  += coef[i]*SumQQ[i*selem+off];
		sumQWt += coef[i];
	      }
	      if (SumUr[i*selem+off]!=fblank) {
		valUr += coef[i]*SumUr[i*selem+off];
		valUi += coef[i]*SumUi[i*selem+off];
		valUU += coef[i]*SumUU[i*selem+off];
		sumUWt += coef[i];
	      }
	      if (SumVr[i*selem+off]!=fblank) {
		valVr += coef[i]*SumVr[i*selem+off];
		valVi += coef[i]*SumVi[i*selem+off];
		valVV += coef[i]*SumVV[i*selem+off];
		sumVWt += coef[i];
	      }
	    }
	  } /* end loop over lists */
	  /* DEBUG 
	  if (iIF==1) {
	    fprintf (stderr, "x %f y %f c %f s %f v %f\n", x, y, closest, sumIWt, valIr/sumIWt );
	  }*/
	  /* Better be something within 0.5 cells */
	  if (closest>0.5) {
	    sumIWt = sumQWt = sumUWt = sumVWt = 0.0;
	  }
	  if (fabs(sumIWt)>0.1) {
	    valIr /= sumIWt;
	    valIi /= sumIWt;
	    valII /= sumIWt;
	  } else {
	    valIr = fblank;
	    valIi = fblank;
	    valII = fblank;
	  }
	  if (fabs(sumQWt)>0.1) {
	    valQr /= sumQWt;
	    valQi /= sumQWt;
	    valQQ /= sumQWt;
	  } else {
	    valQr  = fblank;
	    valQi  = fblank;
	    valQQ = fblank;
	  }
	      if (fabs(sumUWt)>0.1) {
	    valUr /= sumUWt;
	    valUi /= sumUWt;
	    valUU /= sumUWt;
	  } else {
	    valUr = fblank;
	    valUi = fblank;
	    valUU = fblank;
	  }
	      if (fabs(sumVWt)>0.1) {
	    valVr /= sumVWt;
	    valVi /= sumVWt;
	    valVV /= sumVWt;
	  } else {
	    valVr = fblank;
	    valVi = fblank;
	    valVV = fblank;
	  }
	  
	  /* Insert into FArrays */
	  jndx  = iy*nx + ix;
	  /* I */
	  indx  = ichan + iIF*nchan;
	  /* DEBUG 
	  if ((iIF==1) && (ix==21) && (iy==28)) {
	    fprintf (stderr,"x %d y %d r %f i %f jndx %d indx %d\n",ix,iy, valIr,valIi, jndx,indx );
	  }*/
	  array = grids[indx];
	  if ((valIr==fblank) || (valIi==fblank)) {
	    amp = ph = fblank;
	  } else {
	    /* Amplitude and phase within +/- 90 deg */
	    amp = sqrt (valIr*valIr + valIi*valIi);
	    ph  = RAD2DG*atan2(valIi, valIr);
	    if (ph>90.0) {
	      amp = -amp;
	      ph -= 180.0;
	    } else if (ph<-90.0) {
	      amp = -amp;
	      ph += 180.0;
	    }
	  }
	  if (doRMS) {  /* RMS */
	    array->array[jndx] = valII;
	  } else if (doPhase) {  /* phase */
	    array->array[jndx] = ph;
	  } else {      /* Stokes ampl */
	    array->array[jndx] = amp;
	  }
	  /* Q */
	  indx  = ichan + iIF*nchan + nchan*nIF;
	  array = grids[indx];
	  if ((valQr==fblank) || (valQi==fblank)) {
	    amp = ph = fblank;
	  } else {
	    /* Amplitude and phase within +/- 90 deg */
	    amp = sqrt (valQr*valQr + valQi*valQi);
	    ph  = RAD2DG*atan2(valQi, valQr);
	    if (ph>90.0) {
	      amp = -amp;
	      ph -= 180.0;
	    } else if (ph<-90.0) {
	      amp = -amp;
	      ph += 180.0;
	    }
	  }
	  if (doRMS) {  /* RMS */
	    array->array[jndx] = valQQ;
	  } else if (doPhase) {  /* phase */
	    array->array[jndx] = ph;
	  } else {      /* Stokes ampl */
	    array->array[jndx] = amp;
	  }
	  /* U */
	  indx  = ichan + iIF*nchan + 2*nchan*nIF;
	  array = grids[indx];
	  if ((valUr==fblank) || (valUi==fblank)) {
	    amp = ph = fblank;
	  } else {
	    /* Amplitude and phase within +/- 90 deg */
	    amp = sqrt (valUr*valUr + valUi*valUi);
	    ph  = RAD2DG*atan2(valUi, valUr);
	    if (ph>90.0) {
	      amp = -amp;
	      ph -= 180.0;
	    } else if (ph<-90.0) {
	      amp = -amp;
	      ph += 180.0;
	    }
	  }
	  if (doRMS) {  /* RMS */
	    array->array[jndx] = valUU;
	  } else if (doPhase) {  /* phase */
	    array->array[jndx] = ph;
	  } else {      /* Stokes ampl */
	    array->array[jndx] = amp;
	  }
	  /* V */
	  indx  = ichan + iIF*nchan + 3*nchan*nIF;
	  array = grids[indx];
	  if ((valVr==fblank) || (valVi==fblank)) {
	    amp = ph = fblank;
	  } else {
	    /* Amplitude and phase within +/- 90 deg */
	    amp = sqrt (valVr*valVr + valVi*valVi);
	    ph  = RAD2DG*atan2(valVi, valVr);
	    if (ph>90.0) {
	      amp = -amp;
	      ph -= 180.0;
	    } else if (ph<-90.0) {
	      amp = -amp;
	      ph += 180.0;
	    }
	  }
	  if (doRMS) {  /* RMS */
	    array->array[jndx] = valVV;
	  } else if (doPhase) {  /* phase */
	    array->array[jndx] = ph;
	  } else {      /* Stokes ampl */
	    array->array[jndx] = amp;
	  }
	} /* end channel Loop */
      } /* end IF Loop */
    } /* end x loop */
  } /* end y loop */

} /* end gridData */

/*----------------------------------------------------------------------- */
/*  Lagrangian interpolation coefficients                                 */
/*  For interpolating in a quasi regular grid represented by lists        */
/*  Determine coefficients for elements in lists to interpolate to (x,y)  */
/*   Input:                                                               */
/*      x         Coordinate on first axis                                */
/*      y         Coordinate on second axis                               */
/*      n         Length of lists                                         */
/*      hwid      Halfwidth of interpolation kernal                       */
/*      xlist     List of coordinates on  first axis                      */
/*      ylist     List coordinate on second axis                          */
/*   Output:                                                              */
/*      coef      Array of interpolation coefficients for xlist,ylist     */
/*----------------------------------------------------------------------- */
void lagrange(ofloat x, ofloat y, olong n, olong hwid, 
	      ofloat *xlist, ofloat *ylist, ofloat *coef)
{
  ofloat xhwid = (ofloat)hwid, sum;
  odouble prodx, prodxd, prody, prodyd;
  olong  i, j, countx, county;

  /* DEBUG - closest cell 
  for (j=0; j<n; j++) {
    coef[j] = 0.0;
    if ((fabs(x-xlist[j])<=0.5) && (fabs(y-ylist[j])<=0.5)) coef[j] = 1.0;
  }
  return; */
  /* end DEBUG */

  /* Loop over list */
  sum = 0.0;
  for (j=0; j<n; j++) {
    prodx = prodxd = prody = prodyd = 1.0;
    countx = county = 0;
    coef[j] = 0.0;

    /* Within hwid? and i!=j */
    if ((fabs(x-xlist[j])<=xhwid) && (fabs(y-ylist[j])<=xhwid)) {
      coef[j] = 1.0;  /* In case nothing else within hwid */
      countx++;
      prodx  *= (odouble)(x - xlist[j]);
      prodxd *= (odouble)(xlist[j] - xlist[j]);
      county++;
      prody  *= (odouble)(y - ylist[j]);
      prodyd *= (odouble)(ylist[j] - ylist[j]);
      
      /* Inner loop over list */
      for (i=0; i<n; i++) {
	/* X Within hwid? and i!=j */
	/*if (fabs(x-xlist[i])<=xhwid) {*/
	if ((fabs(x-xlist[i])<=xhwid) && 
	    (i!=j) && (fabs(xlist[j]-xlist[i])>0.3)) {
	  countx++;
	  prodx  *= (odouble)(x - xlist[i]);
	  prodxd *= (odouble)(xlist[j] - xlist[i]);
	}
	/* Y Within hwid? and i!=j */
	/*if (fabs(y-ylist[i])<=xhwid) {*/
	if ((fabs(y-ylist[i])<=xhwid) &&
	    (i!=j) &&  (fabs(ylist[j]-ylist[i])>0.3)) {
	  county++;
	  prody  *= (odouble)(y - ylist[i]);
	  prodyd *= (odouble)(ylist[j] - ylist[i]);
	}
      }
    }
    /* put it together */
    if ((countx>=1)  || (county>=1)) {
      if ((prodxd!=0.0) && (prodyd!=0.0))
	coef[j] = (ofloat)(prodx*prody / (prodxd*prodyd));
      else
	coef[j] = 1.0;
      sum += coef[j];
    }
  } /* end loop over list */

  /* Normalize if anything found */
  if (fabs(sum)<0.01) return;
  prodx = 1.0 / sum;
  for (j=0; j<n; j++) coef[j] *= prodx;

  /* DEBUG */
  fprintf(stderr,"lagrange x %f, y %f, sum %f\n", x, y, sum);
  
} /* end lagrange */

