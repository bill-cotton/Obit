/* $Id$  */
/* Obit Filter SN table phases to remove peculiar phases  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2007-2013                                          */
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
#include "ObitTableANUtil.h"
#include "ObitAntennaList.h"
#include "ObitUVSoln.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* SNFiltIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void SNFiltOut (ObitInfoList* outList, ObitErr *err);
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
void SNFiltHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Filter SN table */
void SNFilter (ObitInfoList* myInput, ObitUV* inData,  ObitErr* err);
/* Setup for fitting */
void SNFiltSet (ObitInfoList* myInput, ObitUV* inData, ObitTableSN *SNTab, 
		ObitErr* err);
/* Count times in SN table */
olong timeCount (ObitTableSN *SNTab, ObitErr* err);
/* Ingest SN table */
void ReadSNTab (ObitTableSN *SNTab, ObitErr* err);
/* Fit model */
void FitSNData (gboolean init, ofloat maxRMS, ObitErr* err);
/* Initial Gradient from FT */
void GradFT (olong itime, ofloat *phase, ofloat *gradX, ofloat *gradY, ofloat *resid);
/* Gradient from Least Squares fit */
void GradFit (olong itime, ofloat *phase, ofloat *gradX, ofloat *gradY, ofloat *resid, ofloat *rms);
/* Edit data */
gboolean EditSNData (ofloat maxRMS, gboolean *allBad, ObitErr *err);
/* Write output */
void  WriteSNTab (ObitTableSN *inTab, ObitTableSN *outTab, 
		  gboolean doGrad, gboolean doRes, gboolean doBlank, ObitErr *err);
/* Show model, residuals */
void  ShowModel (ObitErr* err);
/* Smooth gradients */
void  SmooGrad (ObitInfoList* myInput);
/* Unwrap VLA phases down arms */
void VLAUnwrap (ObitErr* err);


/* Program globals */
gchar *pgmName = "SNFilt";       /* Program name */
gchar *infile  = "SNFilt.in" ;   /* File with program inputs */
gchar *outfile = "SNFilt.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
olong prtLv=0;          /* Diagnostic print level */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */

/* Fitting data */
olong numtim;           /* Number of solution times */
olong numant;           /* Number of antennas */
olong numIF;            /* Number of IFs */
olong numpol;           /* Number of polarizations */
olong refant;           /* reference antenna */
olong timinc;           /* Time increment in phase */
olong antinc;           /* Antenna increment in phase, insphs */
olong IFinc;            /* IF increment in phase, insphs */
olong polinc;           /* Poln increment in phase, insphs */
olong count;            /* Total number of valid phases measured */
gboolean *badtim=NULL; /* Valid flags for each time */
odouble uvfreq;        /* Reference frequency in Hz */
odouble *freqs=NULL;   /* Frequencies (Hz) per IF */
odouble *sntime=NULL;  /* Times of SN data (day) */
ofloat  *phase=NULL;   /* SN phases as (time, ant, IF, poln) */
ofloat  *insphs=NULL;  /* Instrumental phase (IF,ant,poln?)*/
ofloat  *gradX=NULL;   /* Gradient in X per time */
ofloat  *gradY=NULL;   /* Gradient in Y per time */
ofloat  timerange[2];  /* Reference time range */
ofloat  maxdis;        /* Maximum antenna difference from reference antenna */
ObitAntennaList *Ant=NULL;     /* Antenna info reprojected to local coordinates */
gboolean doUnwrap=FALSE;      /* Unwrap VLA phases? */

/* Work arrays */
ofloat *vis=NULL;      /* Visibilities for FT gradient fit */
ofloat *uv=NULL;       /* UVs for FT gradient fit */
olong grdsiz = 101;     /* Size of FT grid */
ofloat *grid=NULL;     /* FT grid for FT gradient fit */
olong *xres=NULL;       /* LSQ work array [maxant] */
olong *yres=NULL;       /* LSQ work array [maxant] */
ofloat *resid=NULL;    /* LSQ work array */
ofloat *phz=NULL;      /* LSQ work array */
ofloat *sum1p1=NULL;   /* LSQ work array */
ofloat *sum1p2=NULL;   /* LSQ work array */
ofloat *sum2p1=NULL;   /* LSQ work array */
ofloat *sum2p2=NULL;   /* LSQ work array */
ofloat *sum3p1=NULL;   /* LSQ work array */
ofloat *sum3p2=NULL;   /* LSQ work array */
/*----------------- Macroes ---------------------------*/
/** 
 * Macro to index phase, insphs
 * adds errMsg/traceback to err for errors
 * \li time   0 rel time index, use 0 for insphs
 * \li ant    0 rel antenna index
 * \li ifno   0 rel IF index
 * \li polno  0 rel Poln index
 */
#define ARR_INDEX(time, ant, ifno, polno) time*timinc + ant*antinc + ifno*IFinc + polno*polinc


int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/* Obit task which  filters Solution(SN) tables to remove peculiar phases */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitErr      *err= NULL;

  /* Startup - parse command line */
  err = newObitErr();
  myInput = SNFiltIn (argc, argv, err);
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

  /* Filter Table */
  SNFilter (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Write history */
  SNFiltHistory (myInput, inData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  
  /* Shutdown Obit */
 exit:
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput  = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* SNFiltIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "SNFiltIn";

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
} /* end SNFiltIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: SNFilt -input file -output ofile [args]\n");
    fprintf(stderr, "Filter SN table to remove peculiar phases \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def SNFilt.in\n");
    fprintf(stderr, "  -output output parameter file, def SNFilt.out\n");
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
  strTemp = "SNFilt.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "SNFiltName";
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
/*  Set size of FT search window                                          */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void digestInputs(ObitInfoList *myInput, ObitErr *err)
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat       ftemp;
  gboolean     doCalSelect;
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* Make sure doCalSelect set properly */
  doCalSelect = FALSE;
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);

  /* Set size of FT gradient search */
  ftemp = 10.0;
  ObitInfoListGetTest(myInput, "search", &type, dim, &ftemp);
  ftemp = MIN (5.0, ftemp);
  grdsiz = (olong)(2.0*ftemp + 0.5);
  /* Make Odd */
  grdsiz = 2*grdsiz + 1;

  /* Unwrap wanted? */
  doUnwrap = FALSE;
  ObitInfoListGetTest(myInput, "doUnwrap", &type, dim, &doUnwrap);

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

  /* Full instantiation - needed to initialize tableList */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Write History for SNFilt                                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void SNFiltHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "solnIn", "solnOut", "maxRMS", "doGrad", "doRes", "doBlank", "doUnwrap",
    "timeRange", "refAnt", "width", "alpha", "search", 
    NULL};
  gchar *routine = "SNFiltHistory";

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
 
} /* end SNFiltHistory  */

/*----------------------------------------------------------------------- */
/* Routine to fit gradient + instrumental phases to an SN table and       */
/* replace phases with this model                                         */
/*   Really only works for the VLA or arrays with the same geometry.      */
/* Routine translated from the AIPSish SNFLT.FOR/SNFILT                   */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV with input and output SN tables                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void SNFilter (ObitInfoList* myInput, ObitUV* inData,  ObitErr* err)
{
  ObitTableSN *inTab=NULL, *outTab=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong solnIn, solnOut, highVer;
  olong itemp;
  gboolean drop,allBad,  doGrad, doRes=FALSE, doBlank=FALSE;
  ofloat maxRMS;
  gchar  *routine = "SNFiltCopy";
  
  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  g_assert(ObitUVIsA(inData));

  /* Get input parameters */
  ObitInfoListGetTest(myInput, "prtLv", &type, dim, &prtLv);
  doGrad = FALSE;
  ObitInfoListGetTest(myInput, "doGrad", &type, dim, &doGrad);
  doRes = FALSE;
  ObitInfoListGetTest(myInput, "doRes",  &type, dim, &doRes);
  doBlank = FALSE;
  ObitInfoListGetTest(myInput, "doblank", &type, dim, &doBlank);
  maxRMS = 1.0e20;
  ObitInfoListGetTest(myInput, "maxRMS", &type, dim, &maxRMS);
  if (maxRMS<=1.0e-10) maxRMS = 1.0e20;
 
  /* Table version numbers */
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

  /* Sort to antenna time order */
  ObitTableUtilSort2f ((ObitTable*)inTab, 
		       "TIME  ", 1,  FALSE, 
		       "ANTENNA" , 1, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, inTab->name);

  /* Setup for fitting */
  SNFiltSet(myInput, inData, inTab, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Read Input SN table */
  ReadSNTab (inTab, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Output table - Zap  to be sure */
  if (solnOut>0) ObitDataZapTable ((ObitData*)inData, "AIPS SN", solnOut, err);
  outTab = newObitTableSNValue (inData->name, (ObitData*)inData, &solnOut, 
				OBIT_IO_WriteOnly, inTab->numPol, inTab->numIF, 
				err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Fit SN data */
  FitSNData (TRUE, maxRMS, err);
  if (err->error) Obit_traceback_msg (err, routine, inTab->name);
  ObitErrLog(err); /* show any error messages on err */

  /* Edit Data - any data dropped? */
  drop = EditSNData (maxRMS, &allBad, err);
  if (err->error) Obit_traceback_msg (err, routine, inTab->name);
  ObitErrLog(err); /* show any error messages on err */
  /* Check if all failed */
  Obit_return_if_fail((!allBad), err,
		      "%s: All solutions rejected", routine);

  /* Fit again if data dropped */
  if (drop) FitSNData (FALSE, maxRMS, err);
  if (err->error) Obit_traceback_msg (err, routine, inTab->name);
  ObitErrLog(err); /* show any error messages on err */

  /* Smooth gradients */
  SmooGrad (myInput);

  /* Show fitted model */
  if (prtLv>=2) ShowModel(err);
  ObitErrLog(err); /* show any error messages on err */

  /* Write output SN table */
  WriteSNTab (inTab, outTab, doGrad, doRes, doBlank, err);
  if (err->error) Obit_traceback_msg (err, routine, inTab->name);
} /* end SNFilter */

/**
 * Create global arrays needed for fitting,
 * Reproject array onto plane tangent at reference antenna
 * Routine translated from the AIPSish SNFLT.FOR/SNFSET 
 * \param myInput  Input parameters on InfoList 
 * \param inData   ObitUV with Tables
 * \param SNTab    SN table object 
 * \param Ant      Antenna list
 * \param err      Error/message stack, returns if error.
 */
void SNFiltSet (ObitInfoList* myInput, ObitUV* inData, ObitTableSN *SNTab, 
		ObitErr* err) 
{
  olong i;
  ObitTableAN  *ANTab=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong iANver;
  ofloat dis;
  odouble ArrLong, ArrLat, clat, slat, tx, ty, tz,  stnxr, stnyr, stnzr;
  gchar *routine = "SNFiltSet";

  /* Error checks */
  if (err->error) return;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));
  g_assert (ObitUVIsA(inData));

  /* Parameters */
  timerange[0] = 0.0; timerange[1] = 1.0e20;
  ObitInfoListGetTest(myInput, "timeRange", &type, dim, timerange);
  refant = 1;
  ObitInfoListGetTest(myInput, "refAnt", &type, dim, &refant);

  /* Get number of times */
  numtim = timeCount (SNTab, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);

  /* Numbers of other things */
  numIF  = SNTab->numIF;
  numpol = SNTab->numPol;
  numant = SNTab->numAnt;

  /* Values for indexing */
  polinc = 1;
  IFinc  = numpol*polinc;
  antinc = numIF * IFinc;
  timinc = numant * antinc;

  /* Create data arrays */
  badtim = g_malloc0(numtim*sizeof(gboolean));
  freqs  = g_malloc0(numIF*sizeof(odouble));
  sntime = g_malloc0(numtim*sizeof(odouble));
  phase  = g_malloc0(numtim*numant*numIF*numpol*sizeof(ofloat));
  insphs = g_malloc0(numant*numIF*numpol*sizeof(ofloat));
  gradX  = g_malloc0(numtim*sizeof(ofloat));
  gradY  = g_malloc0(numtim*sizeof(ofloat));

  /* FT gradient fit arrays */
  vis  = g_malloc0(numant*numant*sizeof(ofloat));
  uv   = g_malloc0(numant*numant*sizeof(ofloat));
  grid = g_malloc0(grdsiz*grdsiz*sizeof(ofloat));

  /* Fitting arrays */
  xres   = g_malloc0(numant*sizeof(olong));
  yres   = g_malloc0(numant*sizeof(olong));
  sum1p1 = g_malloc0(numIF*numant*numpol*sizeof(ofloat));
  sum1p2 = g_malloc0(numIF*numant*numpol*sizeof(ofloat));
  sum2p1 = g_malloc0(numtim*sizeof(ofloat));
  sum2p2 = g_malloc0(numtim*sizeof(ofloat));
  sum3p1 = g_malloc0(numtim*sizeof(ofloat));
  sum3p2 = g_malloc0(numtim*sizeof(ofloat));
  phz    = g_malloc0(numtim*numant*numIF*numpol*sizeof(ofloat));
  resid  = g_malloc0(numant*numIF*numpol*sizeof(ofloat));

  /* No bad times yet */
  for (i=0; i<numtim; i++) badtim[i] = FALSE;

  /* Save frequency info */
  ObitUVGetFreq (inData, err);
  if (err->error)  goto cleanup;
  uvfreq = inData->myDesc->freq;
  for (i=0; i<numIF; i++) freqs[i] = inData->myDesc->freqArr[i];

  iANver = 1;
  /* Get table */
  ANTab = newObitTableANValue (inData->name, (ObitData*)inData, &iANver, 
			       OBIT_IO_ReadOnly, 0, 0, 0, err);
  /* Get Antenna list */
  Ant = ObitTableANGetList (ANTab, err);
  if ((err->error) || (ANTab==NULL) || (Ant==NULL)) goto cleanup;

  /* Reprojected location of the  reference antenna */
  ArrLong = Ant->ANlist[refant-1]->AntLong;
  ArrLat  = Ant->ANlist[refant-1]->AntLat;
  clat = cos (ArrLat*DG2RAD);
  slat = sin (ArrLat*DG2RAD);
  tx = Ant->ANlist[refant-1]->AntXYZ[0];
  ty = Ant->ANlist[refant-1]->AntXYZ[1];
  tz = Ant->ANlist[refant-1]->AntXYZ[2];
  stnxr = ty;
  stnyr = tz * clat - tx * slat;
  stnzr = tx * clat + tz * slat;
  /* Reproject onto Plains of San  Augustin  new Y=N, X=E, Z=up */
  for (i=0; i<numant; i++) { /* loop 20 */
    tx = Ant->ANlist[i]->AntXYZ[0];
    ty = Ant->ANlist[i]->AntXYZ[1];
    tz = Ant->ANlist[i]->AntXYZ[2];
    Ant->ANlist[i]->AntXYZ[0] = ty - stnxr;
    Ant->ANlist[i]->AntXYZ[1] = tz * clat - tx * slat - stnyr;
    Ant->ANlist[i]->AntXYZ[2] = tx * clat + tz * slat - stnzr;
  } /* end loop  L20:  */;

  /* Find maximum antenna distance from center */
  maxdis = 0.0;
  for (i=0; i<numant; i++) { /* loop 210 */
    dis = sqrt (Ant->ANlist[i]->AntXYZ[0]*Ant->ANlist[i]->AntXYZ[0] +
		Ant->ANlist[i]->AntXYZ[1]*Ant->ANlist[i]->AntXYZ[1]);
    maxdis = MAX (maxdis, dis);
  } /* end loop  L210: */;

  /* deallocate structures */
 cleanup:
  ANTab = ObitTableANUnref(ANTab); 
} /* end of routine SNFiltSet */ 

/**
 * Routine to count number of distinct times in SN table
 * All valid entries are included.  
 * Routine translated from the AIPSish REFCNT.FOR/REFCNT  
 * \param SNTab    SN table object 
 * \param err      Error/message stack, returns if error.
 * \return  number of distinct times in table. 
 */
olong timeCount (ObitTableSN *SNTab, ObitErr* err) 
{
  olong   outCount = 0;
  ObitIOCode retCode;
  olong  loop;
  ofloat lastTime;
  ObitTableSNRow *row=NULL;
  gchar *routine = "timeCount";

  /* Error checks */
  if (err->error) return outCount;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));

  /* Open table */
  retCode = ObitTableSNOpen (SNTab, OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, SNTab->name, outCount);

  lastTime = -1.0e20;
  count = 0;
  
  /* Create Row */
  row = newObitTableSNRow (SNTab);

  /* Loop through table */
  for (loop=1; loop<=SNTab->myDesc->nrow; loop++) { /* loop 20 */

    retCode = ObitTableSNReadRow (SNTab, loop, row, err);
    if (err->error) break;
    if (row->status<0) continue;  /* Skip deselected record */

    /* Count times */
    if (row->Time>(lastTime+0.5*row->TimeI)) {
      lastTime = row->Time;
      outCount ++;
    }
  } /* end loop over table */
  /* Close SN table */
  ObitTableSNClose (SNTab, err);
  if (err->error) Obit_traceback_val (err, routine, SNTab->name, outCount);

  row = ObitTableSNRowUnref(row); /* delete row object */
  if (err->error) Obit_traceback_val (err, routine, SNTab->name, outCount);

  return outCount;
} /* end of routine timeCount */ 

void ReadSNTab (ObitTableSN *SNTab, ObitErr* err)
/*----------------------------------------------------------------------- */
/*  Read SN table saving values                                           */
/*  Routine translated from the AIPSish SNFLT.FOR/SNFRED                  */
/*   Input:                                                               */
/*      SNTab  SN table to swallow                                        */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
{
  ObitIOCode retCode;
  ObitTableSNRow *row=NULL;
  olong i, numb, loop, itime, iif, iant, indx;
  ofloat fblank = ObitMagicF();
  odouble time;
  gchar *routine = "ReadSNTab";

  /* Error checks */
  if (err->error) return ;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));

  /* Blank fill phase */
  numb = numant*numtim*numIF*numpol;
  for (i=0; i<numb; i++) phase[i] = fblank;

  /* Open table */
  retCode = ObitTableSNOpen (SNTab, OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, SNTab->name);

  /* Create Row */
  row = newObitTableSNRow (SNTab);

  time = -1.0e20;
  count = 0;
  itime = -1;

  /* Loop through table */
  for (loop=1; loop<=SNTab->myDesc->nrow; loop++) { /* loop L600: */

    retCode = ObitTableSNReadRow (SNTab, loop, row, err);
    if (err->error) break;
    if (row->status<0) continue;  /* Skip deselected record */

    /* Is this a new time? */
    if (row->Time > time) {
      /* If hit limit, stop reading here, complain  and use what you've got */
      if (itime >= numtim) {
	Obit_log_error(err, OBIT_InfoWarn, "%s: exceed max times %d", 
		       routine, numtim);
	break;
      }
      itime++;
      sntime[itime] = row->Time;
      time          = row->Time;
    } /* end new time */ 

    /* Save phases */
    iant = row->antNo-1;
    for (iif=0; iif<numIF; iif++) { /* loop 420 */
      indx = ARR_INDEX(itime,iant,iif,0);
      if ((row->Weight1[iif] > 0.0)  &&  (row->Imag1[iif] != fblank)  &&  
	  (row->Real1[iif] != fblank)) {
	/* Phase OK */
	phase[indx] = atan2 (row->Imag1[iif], row->Real1[iif]+1.0e-20);
	count++;
      } else {
	/* Phase bad */
	phase[indx] = fblank;
      } 
      if (numpol>1) { /* Second polarization? */
      indx = ARR_INDEX(itime,iant,iif,1);
      if ((row->Weight2[iif] > 0.0)  &&  (row->Imag2[iif] != fblank)  &&  
	    (row->Real2[iif] != fblank)) {
	  /* Phase OK */
	  phase[indx] = atan2 (row->Imag2[iif], row->Real2[iif]+1.0e-20);
	  count++;
	} else {
	  /* Phase bad */
	  phase[indx] = fblank;
	} 
      } /* end second poln */
    } /* end loop  L420: */

 } /* end loop  L600: over table */
  
  /* Close SN table */
  ObitTableSNClose (SNTab, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);

  row = ObitTableSNRowUnref(row); /* delete row object */
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
} /* end ReadSNTab */

/**
 * Fit model of instrumental phase plus gradients.  
 * Solution uses a relaxation method from Fred Schwab:  
 * Pn+1 = Pn + atan {(dChi2/dP) / (d2Chi2/dP2)}  
 * for each parameter P where n or n+1 indicates a given iteration,   
 * dChi2/dP is the first partial derivative of Chi squared wrt P,  
 * d2Chi2/d2P is the second partial derivative of Chi squared wrt P,  
 * Chi2 = Sum (Abs(I + gxX + gyY - sn_phase)**2),  
 * I is the instrumental phase, gx the x gradient, X the X array  
 * coordinate, gy the y gradient, Y the Y array coordinate and  
 * sn_phase is the measured phase for this time, antenna, IF and  
 * poln.   
 * Routine translated from the AIPSish SNFLT.FOR/SNFITT  
 * \param init    If true initialize the solutions 
 * \param maxRMS  Maximum allowable RMS in deg. 
 * \param err     Obit Error stack     
 */
void FitSNData (gboolean init, ofloat maxRMS, ObitErr* err)   {
  olong i, gcnt, indx, jndx;
  olong   ipol, itime, iif, iant, count, iter, rmscnt, itmp;
  gboolean   convgd, allBad=FALSE;
  ofloat      sumr, sumi, w, arg, delta, tol, test, rms, normg, 
    normi, delx, dely;
  ofloat  stnx, stny;
  ofloat fblank = ObitMagicF();
  gchar *routine = "FitSNData";

  /* Need to Initialize? */
  if (init) {
    /* Initial guess,  Zero gradients */
    for (i=0; i<numtim; i++) {gradX[i] = 0.0; gradY[i] = 0.0;}
    
    /* Average phases in specified  timerange  */
    gcnt = 0;
    for (ipol=0; ipol<numpol; ipol++) { /* loop 40 */
      for (iif=0; iif<numIF; iif++) { /* loop 30 */
	for (iant=0; iant<numant; iant++) { /* loop 20 */
	  jndx = ARR_INDEX(0,iant,iif,ipol);
	  insphs[jndx] = 0.0;
	  count = 0;
	  sumr = 0.0;
	  sumi = 0.0;
	  for (itime=0; itime<numtim; itime++) { /* loop 10 */
	    indx = ARR_INDEX(itime,iant,iif,ipol);
	    if ((sntime[itime] >= timerange[0])  &&  
		(sntime[itime] <= timerange[1])) {
	      if (phase[indx] != fblank) {
		count++;
		sumr += cos (phase[indx]);
		sumi += sin (phase[indx]);
	      } 
	    } 
	  } /* end loop  L10:  */;
	  if (count > 0) {
	    gcnt++;
	    sumr = sumr / count;
	    sumi = sumi / count;
	    insphs[jndx] = atan2 (sumi, sumr+1.0e-20);
	  } else {
	    insphs[jndx] = fblank;
	  } 
	} /* end loop  L20:  */;
      } /* end loop  L30:  */;
    } /* end loop  L40:  */
    
    /* Check if no data */
    Obit_return_if_fail((gcnt > 0), err,
			"%s: NO Valid data in timerange", routine);
    
    /* If VLA Unwrap phases */
    if (doUnwrap) VLAUnwrap (err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
 
    /* Make copy of data with average  phase removed. */
    for (iant=0; iant<numant; iant++) { /* loop 60 */
      for (itime=0; itime<numtim; itime++) { /* loop 58 */
	for (iif=0; iif<numIF; iif++) { /* loop 56 */
	  for (ipol=0; ipol<numpol; ipol++) { /* loop 54 */
	    indx = ARR_INDEX(itime,iant,iif,ipol);
	    jndx = ARR_INDEX(0,iant,iif,ipol);
	    if ((phase[indx] != fblank)  
		&& (insphs[jndx] != fblank))  {
	      phz[indx] = phase[indx] - insphs[jndx];
	    } else {
	      phz[indx] = fblank;
	    } 
	  } /* end loop  L54:  */
	} /* end loop  L56:  */
      } /* end loop  L58:  */
    } /* end loop  L60:  */
    
    /* Fit Gradients */
    allBad = TRUE;
    for (itime=0; itime<numtim; itime++) { /* loop 95 */
      /* Initial value with FT */
      GradFT (itime, phz, &gradX[itime], &gradY[itime], resid) ;
      
      /* Refine fit with least squares */
      GradFit (itime, phz, &gradX[itime], &gradY[itime], resid, &rms) ;
      
      /* Weed out bad times */
      if (rms > maxRMS) badtim[itime] = TRUE;
      else allBad = FALSE;

      /*   Diagnostics */
      if (prtLv>=2)
	Obit_log_error(err, OBIT_InfoErr, 
		       "time %5d refraction %10.1f %10.1f rms %10.1f",
		       itime+1, 1.330e5*gradX[itime], 1.330e5*gradY[itime], rms);
      ObitErrLog(err); /* show any error messages on err */
    } /* end loop  L95:  */;
    
  } /* End of initialization */
  ObitErrLog(err); /* show any error messages on err */

  /* Check if all bad */
  if (allBad) {
    Obit_log_error(err, OBIT_Error, "ALL solutions rejected");
    return;
  }

  normg = 0.0;
  normi = 0.0;

  /* Loop over iterations */
  for (iter= 1; iter<=300; iter++) { /* loop 600 */
    
    convgd = TRUE;  /* Can decide otherwise: */

    /* Instrumental phase sums */
    for (ipol=0; ipol<numpol; ipol++) { /* loop 140 */
      for (iif=0; iif<numIF; iif++) { /* loop 130 */
	for (iant=0; iant<numant; iant++) { /* loop 120 */
	  stnx = Ant->ANlist[iant]->AntXYZ[0];
	  stny = Ant->ANlist[iant]->AntXYZ[1] ;
	  jndx = ARR_INDEX(0,iant,iif,ipol);
	  sum1p1[jndx] = 0.0;
	  sum1p2[jndx] = 1.0e-20;
	  for (itime=0; itime<numtim; itime++) { /* loop 110 */
	    indx = ARR_INDEX(itime,iant,iif,ipol);
	    if (badtim[itime]) continue;
	    if ((phase[indx] != fblank)  && (insphs[jndx] != fblank)) {
	      arg = insphs[jndx] +  
		gradX[itime] * stnx + gradY[itime] * stny - 
		phase[indx];
	      
	      /* Keep in range */
	      if (!doUnwrap) {
		itmp = arg / (2.0*G_PI);
		arg  -= itmp *  2.0*G_PI;
		if (arg > G_PI)  arg -= 2.0*G_PI;
		if (arg < -G_PI) arg += 2.0*G_PI;
	      }
	      sum1p1[jndx] += arg;
	      sum1p2[jndx] += 1.0;
	    } 
	  } /* end loop  L110: */
	} /* end loop  L120: */
      } /* end loop  L130: */
    } /* end loop  L140: */
  
    /* Gradient sums */
    rmscnt = 0;
    rms = 0.0;
    for (itime=0; itime<numtim; itime++) { /* loop 190 */
      if (badtim[itime]) continue;
      sum2p1[itime] = 0.0;
      sum2p2[itime] = 1.0e-20;
      sum3p1[itime] = 0.0;
      sum3p2[itime] = 1.0e-20;
      for (ipol=0; ipol<numpol; ipol++) { /* loop 180 */
	for (iif=0; iif<numIF; iif++) { /* loop 170 */
	  for (iant=0; iant<numant; iant++) { /* loop 160 */
	    jndx = ARR_INDEX(0,iant,iif,ipol);
	    indx = ARR_INDEX(itime,iant,iif,ipol);
	    stnx = Ant->ANlist[iant]->AntXYZ[0];
	    stny = Ant->ANlist[iant]->AntXYZ[1];
	    if ((phase[indx] != fblank)  &&  
		(insphs[jndx] != fblank)) {
	      arg = insphs[jndx] +  
		gradX[itime] * stnx + gradY[itime] * stny - 
		phase[indx];
	      /* Keep in range */
	      if (!doUnwrap) {
		itmp = arg / (2.0*G_PI);
		arg -= itmp * 2.0*G_PI;
		if (arg > G_PI)  arg -= 2.0*G_PI;
		if (arg < -G_PI) arg += 2.0*G_PI;
	      }
	      sum2p1[itime] += stnx * arg ;
	      sum2p2[itime] += stnx * stnx;
	      sum3p1[itime] += stny * arg;
	      sum3p2[itime] += stny * stny;
	      rmscnt++;
	      rms += arg*arg;
	    } /* end if valid data */ 
	  } /* end loop  L160: */
	} /* end loop  L170: */
       } /* end loop  L180: */
    } /* end loop  L190: */

    /* Diagnostics */
    if (prtLv>=3) {
      Obit_log_error(err, OBIT_InfoErr, 
		     "Iteration %d grad. norm, %g %g RMS %f",
		     iter, normg, normi,57.296*sqrt (rms/rmscnt));
      ObitErrLog(err); /* show any error messages on err */
    }

    /* Update solutions  Damping factor for solutions */
    w = MIN (0.5, 0.01 * iter);

    /* Convergence criterion - lower the bar  */
    tol = 5.0e-6 + iter * 1.0e-5;
    normg = 0.0;
    normi = 0.0;

    /* Gradients */
    for (itime=0; itime<numtim; itime++) { /* loop 210 */
      if (badtim[itime]) continue;
      delx = atan (sum2p1[itime] / sum2p2[itime]);
      test = tol * fabs (gradX[itime]);
      gradX[itime] -= w * delx;

      /* Convergence test */
      convgd = convgd  &&  (fabs(delx) <= test);
      normg +=  delx*delx;
      dely = atan (sum3p1[itime] / sum3p2[itime]);
      test = tol * fabs (gradY[itime]);
      gradY[itime] -= w * dely;
      /* Convergence test */
      convgd = convgd  &&  (fabs(dely) <= test);
      normg += dely*dely;

    } /* end loop  L210: updating gradients */

    if (iter < 20) continue; /* No instrumental terms before 20 iterations */

    /* Instrumental phases */
    for (ipol=0; ipol<numpol; ipol++) { /* loop 240 */
      for (iif=0; iif<numIF; iif++) { /* loop 230 */
	for (iant=0; iant<numant; iant++) { /* loop 220 */
	  jndx = ARR_INDEX(0,iant,iif,ipol);
	  if (insphs[jndx] != fblank) {
	    delta = atan (sum1p1[jndx] / sum1p2[jndx]) ;
	    insphs[jndx] -= w * delta;
	    /* Convergence test */
	    convgd = convgd  &&  (fabs(delta) <= tol);
	    normi += delta*delta;
	    /* Keep in range */
	    itmp = insphs[jndx] / (2.0*G_PI);
	    insphs[jndx] -= itmp * 2.0*G_PI;
	  } 
	} /* end loop  L220: */;
      } /* end loop  L230: */;
    } /* end loop  L240: */;

    /* Converged? Need at least 22  iterations */
    if (convgd  &&  (iter >= 22)) break;
    } /* end iteration loop  L600: */

    /* Set RMS */
    if (rmscnt > 2) {
      rms = 57.296*sqrt (rms/rmscnt);
    } else {
      rms = -1.0;
    } 

    /* Diagnostics */
    if (prtLv>=1)
      Obit_log_error(err, OBIT_InfoErr, 
		     "after %d iterations the rms residual = %f",
		     iter, rms);
    ObitErrLog(err); /* show any error messages on err */
} /* end of FitSNData */ 

/**
 * Fit model of gradients Using Fourier transform , Average in IF, Stokes  
 * Routine translated from the AIPSish SNFLT.FOR/GRDFT  
 * \param itime     Which time in phase 
 * \param phase     Phases in radians, fblank if absent 
 * \param gradX     X Phase gradient rad/m, FBLANK if none 
 * \param gradY     Y Phase gradient rad/m, FBLANK if none 
 * \param resid     Residuals (ant,IF,poln) from model (not used)
 */
void GradFT (olong itime, ofloat *phase, ofloat *gradX, ofloat *gradY, ofloat *resid)
{
  olong   i, j, iant, jant, ipol, iif, ibase, nbase, xmax, ymax, vcount, indx, jndx;
  gboolean   good;
  ofloat bl, mxbl, cells, x, y, xcen, ycen, phas, yphas, peak;
  ofloat  tvis[2]; 
  ofloat fblank = ObitMagicF();

  /* Initial values */
  *gradX = 0.0;
  *gradY = 0.0;
  
  /* Make pseudo visibilities by  differencing antenna phases  */
  nbase = -1;
  mxbl  = 0.0;
  for (jant=0; jant<numant-1; jant++) { /* loop 40 */
    for (iant=jant+1; iant<numant; iant++) { /* loop 30 */
      /* Any data for this pair? */
      good    = FALSE;
      vcount  = 0;
      tvis[0] = 0.0;
      tvis[1] = 0.0;

      for (ipol=0; ipol<numpol; ipol++) { /* loop 20 */
	for (iif=0; iif<numIF; iif++) { /* loop 10 */
	  indx = ARR_INDEX(itime,iant,iif,ipol);
	  jndx = ARR_INDEX(itime,jant,iif,ipol);
	  if ((phase[indx] != fblank)  &&   
	      (phase[jndx] != fblank)) {
	    good = TRUE;
	    vcount++;
	    tvis[0] += cos (phase[indx] - phase[jndx]);
	    tvis[1] += sin (phase[indx] - phase[jndx]);
	  } 
	} /* end loop  L10:  */;
      } /* end loop  L20:  */;
      
      if (good) {
	nbase++;
	vis[nbase*2]   = tvis[0] / vcount;
	vis[nbase*2+1] = tvis[1] / vcount;
	uv[nbase*2]    = Ant->ANlist[iant]->AntXYZ[0] - Ant->ANlist[jant]->AntXYZ[0];
	uv[nbase*2+1]  = Ant->ANlist[iant]->AntXYZ[1] - Ant->ANlist[jant]->AntXYZ[1];
	bl   = uv[nbase*2]*uv[nbase*2] + uv[nbase*2+1]*uv[nbase*2+1];
	mxbl = MAX (mxbl, bl);
	uv[nbase*2]   *= 2.0*G_PI;
	uv[nbase*2+1] *= 2.0*G_PI;
      } 
    } /* end loop  L30:  */
  } /* end loop  L40:  */

  /* Better be some */
  if (nbase < 10) return;

  /* Zero Grid */
  for (j=0; j<grdsiz; j++) { /* loop 60 */
    for (i=0; i<grdsiz; i++) { /* loop 50 */
      grid[j*grdsiz+i] = 0.0;
    } /* end loop  L50:  */
  } /* end loop  L60:  */

  /* Cellspacing 1/2 turn on longest  baseline */
  cells = 1.0 / (2.0 * sqrt (mxbl));
  xcen = grdsiz / 2;
  ycen = grdsiz / 2;

  /* Loop over pseudo visibilities - making cosine transform */
  for (ibase=0; ibase<nbase; ibase++) { /* loop 200 */
    for (j=0; j<grdsiz; j++) { /* loop 120 */
      y = (ycen - j) * cells;
      yphas =  uv[ibase*2+1]*y;
      for (i=0; i<grdsiz; i++) { /* loop 110 */
	x = (xcen - i) * cells;
	phas = uv[ibase*2]*x + yphas;
	grid[j*grdsiz+i] +=  cos (phas)*vis[ibase*2] - sin (phas)*vis[ibase*2+1];
      } /* end loop  L110: */
    } /* end loop  L120: */
  } /* end loop  L200: */

  /* Find Peak */
  peak = -1.0e20;
  xmax = xcen + 0.5;
  ymax = ycen + 0.5;
  for (j=0; j<grdsiz; j++) { /* loop 220 */
    for (i=0; i<grdsiz; i++) { /* loop 210 */
      if (grid[j*grdsiz+i] > peak) {
	peak = grid[j*grdsiz+i];
	xmax = i;
	ymax = j;
      } 
    } /* end loop  L210: */
  } /* end loop  L220: */

  /* Set peak as gradient */
  *gradX = (xmax - xcen) * 2.0*G_PI * cells;
  *gradY = (ymax - ycen) * 2.0*G_PI * cells;
} /* end of GradFT */ 

/**
 * Fit model of gradients.  
 * Solution uses a relaxation method from Fred Schwab:  
 * Pn+1 = Pn + atan {(dChi2/dP) / (d2Chi2/dP2)}  
 * for each parameter P where n or n+1 indicates a given iteration,   
 * dChi2/dP is the first partial derivative of Chi squared wrt P,  
 * d2Chi2/d2P is the second partial derivative of Chi squared wrt P,  
 * Chi2 = Sum (Abs(gxX + gyY - sn_phase)**2),  
 * gx the x gradient, X the X array coordinate,   
 * gy the y gradient, Y the Y array coordinate and  
 * sn_phase is the measured phase.  
 * Routine translated from the AIPSish SNFLT.FOR/GRDFIT  
 * \param itime     Which time in phase 
 * \param phase     Phases in radians, fblank if absent 
 * \param gradX     X Phase gradient rad/m, FBLANK if none 
 * \param gradY     Y Phase gradient rad/m, FBLANK if none 
 * \param resid     [out] Residuals (ant,IF,poln) from model
 */
void GradFit (olong itime, ofloat *phase, ofloat *gradX, ofloat *gradY, 
	      ofloat *resid, ofloat *rms)
{
  gboolean   convgd;
  olong   iter, iant, ipol, iif, rmscnt, itmp, indx;
  ofloat w, arg, tol, test, normg, normi, delx, dely;
  ofloat stnx, stny;
  ofloat sumg2p1, sumg2p2,  sumg3p1, sumg3p2;
  ofloat fblank = ObitMagicF();

  /* Initialize residuals to data */
  for (iant=0; iant<numant; iant++) { /* loop 30 */
    for (ipol=0; ipol<numpol; ipol++) { /* loop 20 */
      for (iif=1; iif<numIF; iif++) { /* loop 10 */
	indx = ARR_INDEX(itime,iant,iif,ipol);
	resid[ARR_INDEX(0,iant,iif,ipol)] =  phase[indx];
      } /* end loop  L10:  */
    } /* end loop  L20:  */
  } /* end loop  L30:  */
  
  
  
  /* Loop over iterations */
  for (iter= 1; iter<=250; iter++) { /* loop 600 */
	
    /* Can decide otherwise: */
    convgd = TRUE;
    /* Gradient sums */
    rmscnt = 0;
    *rms = 0.0;
    sumg2p1 = 0.0;
    sumg2p2 = 1.0e-20;
    sumg3p1 = 0.0;
    sumg3p2 = 1.0e-20;
    for (ipol=0; ipol<numpol; ipol++) { /* loop 180 */
      for (iif=0; iif<numIF; iif++) { /* loop 170 */
	for (iant=0; iant<numant; iant++) { /* loop 160 */
	  stnx = Ant->ANlist[iant]->AntXYZ[0];
	  stny = Ant->ANlist[iant]->AntXYZ[1] ;
	  indx = ARR_INDEX(itime,iant,iif,ipol);
	  if (phase[indx] != fblank) {
	    arg = (*gradX) * stnx + (*gradY) * stny - phase[indx];
	    /* Keep in range */
	    if (!doUnwrap) {
	      itmp = arg / (2.0*G_PI);
	      arg  -= itmp *  2.0*G_PI;
	      if (arg > G_PI)  arg -= 2.0*G_PI;
	      if (arg < -G_PI) arg += 2.0*G_PI;
	    }
	    sumg2p1 += stnx * arg;
	    sumg2p2 += stnx * stnx;
	    sumg3p1 += stny * arg;
	    sumg3p2 += stny * stny;
	    rmscnt++;
	    (*rms) += arg*arg;
	  } 
	} /* end loop  L160: */
      } /* end loop  L170: */
    } /* end loop  L180: */

    /* Check that there is some valid  data */
    if (rmscnt < 5) return;

    /* Update solutions */
    w = MIN (0.5, 0.01 * iter);

    /* Convergence criterion - lower  the bar  */
    tol = 5.0e-6 + iter * 1.0e-5;
    normg = 0.0;
    normi = 0.0;

    /* Gradients */
    delx = atan (sumg2p1 / sumg2p2);
    test = tol * fabs (*gradX);
    (*gradX) -=  w * delx;

    /* Convergence test */
    convgd = convgd  &&  (fabs(delx) <= test);
    normg += delx*delx;
    dely = atan (sumg3p1 / sumg3p2);
    test = tol * fabs (*gradY);
    (*gradY) -= w * dely;
    
    /* Convergence test */
    convgd = convgd  &&  (fabs(dely) <= test);
    normg += dely*dely;

    if (iter < 20) continue;
    /* Converged? Need at least 22 iterations */
    if (convgd  &&  (iter >= 22)) break;
  } /* end loop iteration loop L600: */

  /*  Set RMS residual */
  if (rmscnt > 2) {
    (*rms) = RAD2DG*sqrt ((*rms)/rmscnt);
  } else {
    (*rms) = -1.0;
  } 
} /* end of GradFit */ 

/**
 * Edit Phase data based on time RMS residual from model  
 * Routine translated from the AIPSish SNFLT.FOR/SNFEDT  
 * \param maxRMS  Maximum allowable RMS in deg. 
 * \param allBad  [out] if true all data rejected
 * \param err     Obit Error stack     
 * \return drop   Some times dropped 
 */
gboolean EditSNData (ofloat maxRMS, gboolean *allBad, ObitErr *err) 
{
  gboolean drop = FALSE;
  olong   ipol, itime, iif, iant,itmp, count, drpcnt, it1, it2;
  ofloat res, sum2, rms, rtemp, rt1, stnx, stny;
  ofloat fblank = ObitMagicF();

  if (prtLv>=1)
    Obit_log_error(err, OBIT_InfoErr, 
		   "removing times with residual rms > %f deg", maxRMS);
  *allBad = FALSE;

  drpcnt = 0;
  /* Loop over times */
  for (itime=0; itime<numtim; itime++) { /* loop 600 */
    if ((gradX[itime] == fblank)  ||  (gradY[itime] == fblank)  
	|| badtim[itime]) continue;

    /* Clear accumulators */
    count = 0;
    sum2 = 0.0;

    /* Accumulate time residuals */
    for (ipol=0; ipol<numpol; ipol++) { /* loop 180 */
      for (iif=0; iif<numIF; iif++) { /* loop 170 */
	for (iant=0; iant<numant; iant++) { /* loop 160 */
	  stnx = Ant->ANlist[iant]->AntXYZ[0];
	  stny = Ant->ANlist[iant]->AntXYZ[1] ;
	  if ((phase[ARR_INDEX(itime,iant,iif,ipol)] != fblank)  &&   
	      (insphs[ARR_INDEX(0,iant,iif,ipol)] != fblank)) {
	    res = insphs[ARR_INDEX(0,iant,iif,ipol)] +  
	      gradX[itime] * stnx + gradY[itime] * stny - 
	      phase[ARR_INDEX(itime,iant,iif,ipol)];
	    /* Keep in range */
	    if (!doUnwrap) {
	      itmp = res / (2.0*G_PI);
	      res -= itmp * 2.0*G_PI;
	      if (res > G_PI)  res -= 2.0*G_PI;
	      if (res < -G_PI) res += 2.0*G_PI;
	    }
	    count++;
	    sum2 += res*res;
	  } 
	} /* end loop  L160: */
      } /* end loop  L170: */
    } /* end loop  L180: */

    /* Time RMS */
    if (count > 6) {
      rms = sqrt (sum2 / count);
      /* To degrees */
      rms *= RAD2DG;
    } else {
      rms = 0.0;
    } 
    /* Tell  */
    rtemp = 24.0*sntime[itime];
    it1 = rtemp;
    rtemp = (rtemp - it1)*60.0;
    it2 = rtemp;
    rt1 = (rtemp - it2)*60.0;
    if (prtLv>=1)
      Obit_log_error(err, OBIT_InfoErr, 
		     "%5d time: %3d %3d %6.2f rms=%7.1f no obs=%4d", 
		   itime+1, it1, it2, rt1, rms, count);
    /* Gradient in asec at zenith. */
    if (prtLv>=1)
      Obit_log_error(err, OBIT_InfoErr, 
		     " refraction in x =%10.2f in y =%10.2f",
		     1.330e5*gradX[itime], 1.330e5*gradY[itime]);
  
   /* Do we chuck it? */
    if (rms > maxRMS) {
      drpcnt = drpcnt + 1;
      if (prtLv>=1)
	Obit_log_error(err, OBIT_InfoErr, 
		       "remove preceeding time");

      badtim[itime] = TRUE;
      for (ipol=0; ipol<numpol; ipol++) { /* loop 280 */
	for (iif=0; iif<numIF; iif++) { /* loop 270 */
	  for (iant=0; iant<numant; iant++) { /* loop 260 */
	    phase[ARR_INDEX(itime,iant,iif,ipol)] = fblank;
	  } /* end loop  L260: */
	} /* end loop  L270: */
      } /* end loop  L280: */
    } /* toss it out */ 
    ObitErrLog(err); /* show any error messages on err */
  } /* end time loop  L600: */

  /* Final result */
  drop = drpcnt>0;
  *allBad = (drpcnt >= numtim);

  if (prtLv>=1)
    Obit_log_error(err, OBIT_InfoErr, 
		   "drop  %5d of %5d times", drpcnt, numtim);
  return drop;
} /* end  EditSNData  */ 

/**
 * Read input SN table and replace phases writing new output 
 * Routine translated from the AIPSish SNFLT.FOR/SNFWRI  
 * \param inTab   Input SN table to swallow 
 * \param outTab  Output SN table to swallow 
 * \param doGrad  True if ionospheric gradients are to be added
 * \param doRes   True if residuals to model requested.
 * \param doBlank True if to ignore prior blanking.
 * \param err     Obit Error stack     
 */
void WriteSNTab (ObitTableSN *inTab, ObitTableSN *outTab, 
		 gboolean doGrad, gboolean doRes, gboolean doBlank, 
		 ObitErr *err)
{
  olong   loop, i, ipol, iif, itime, osnrow, antno;
  ofloat amp, faz, refiph;
  gboolean  bad, OK;
  olong inSNRow, outSNRow, indx, jndx;
  ObitIOCode retCode;
  ObitTableSNRow *inRow=NULL, *outRow=NULL;
  odouble ttol;
  ofloat stnx, stny, fblank = ObitMagicF();
  gchar *routine = "SWriteSNTab";

  /* Error checks */
  if (err->error) return ;  /* previous error? */
  g_assert(ObitTableSNIsA(inTab));
  g_assert(ObitTableSNIsA(outTab));

  /* Subtract ref Ant. phase if instrumental is what is going into output table */
  if (!doRes) {
    for (ipol=0; ipol<numpol; ipol++) { /* loop 30 */
      for (iif=0; iif<numIF; iif++) { /* loop 25 */
	antno = refant-1;
	jndx = ARR_INDEX(0,antno,iif,ipol);
	refiph = insphs[jndx];
	for (i=0; i<numant; i++) { /* loop 20 */
	  jndx = ARR_INDEX(0,i,iif,ipol);
	  if (insphs[jndx] != fblank) {
	    insphs[jndx] -= refiph;
	  } 
	} /* end loop  L20:  */;
      } /* end loop  L25:  */;
    } /* end loop  L30:  */;
  } /* end doRes */

  /* Two second tolerance. for matchex s in tables  */
  ttol = 2*1.157407e-5;

  /* Open tables */
  retCode = ObitTableSNOpen (inTab, OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, inTab->name);
  inRow = newObitTableSNRow (inTab);  /* Create Row */

  retCode = ObitTableSNOpen (outTab, OBIT_IO_WriteOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, outTab->name);
  outTab->numAnt = inTab->numAnt;  /* Number of antennas */
  outRow = newObitTableSNRow (outTab);     /* Create Row */
  ObitTableSNSetRow (outTab, outRow, err);  /* Attach */
  if (err->error) goto cleanup;

  /* Loop over table */
  itime = 1;
  osnrow = 0;
  for (loop= 1; loop<=inTab->myDesc->nrow; loop++) { /* loop 600 */

    /* Read */
    inSNRow = loop;
    retCode = ObitTableSNReadRow (inTab, inSNRow, inRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    if (inRow->status==-1) continue;

    /* Update output row */
    outRow->Time     = inRow->Time;
    outRow->TimeI    = inRow->TimeI;
    outRow->SourID   = inRow->SourID;
    outRow->antNo    = inRow->antNo;
    outRow->SubA     = inRow->SubA;
    outRow->FreqID   = inRow->FreqID;
    outRow->IFR      = inRow->IFR;
    outRow->NodeNo   = inRow->NodeNo;
    outRow->MBDelay1 = inRow->MBDelay1;
    if (numpol>1) {
      outRow->MBDelay2 = inRow->MBDelay2;
    }

    /* Get time number - two second tolerance.  */
    if (fabs (inRow->Time - sntime[itime])  >  ttol) {
      OK = FALSE;
      for (i=0; i<numtim; i++) { /* loop 40 */
	itime = i;
	if (fabs (inRow->Time-sntime[itime]) <=  ttol) {OK = TRUE; break;}
      } /* end loop  L40:  */
      /* Trouble if not OK */
      if (!OK) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: cannot find time %10.6f in internal table",
		       routine, inRow->Time);
	bad = TRUE;
	goto cleanup;
      }
    }
    /* Is this time OK? */
    if (badtim[itime]) continue;
 
    antno = inRow->antNo-1;
    stnx  = Ant->ANlist[antno]->AntXYZ[0];
    stny  = Ant->ANlist[antno]->AntXYZ[1] ;

    /* Update phases */
    for (iif=0; iif<numIF; iif++) { /* loop 60 */
      /* Copy stuff for poln 1 */
      outRow->Delay1[iif]  = inRow->Delay1[iif];
      outRow->Rate1[iif]   = inRow->Rate1[iif];
      outRow->Weight1[iif] = inRow->Weight1[iif];
      outRow->RefAnt1[iif] = inRow->RefAnt1[iif];
      /* in case no solution - blank */
      jndx = ARR_INDEX(0,antno,iif,0);
      if (!doBlank && !doRes) {  /*replace blanking unless giving residuals */
	outRow->Real1[iif]   = cos(insphs[jndx]);
	outRow->Imag1[iif]   = sin(insphs[jndx]);
	outRow->Weight1[iif] = 1.0;
      } else { /* blanking */
	outRow->Real1[iif] = fblank;
	outRow->Imag1[iif] = fblank;
      }
      
      /* Is previous solution for pol 1 good? */
      if ((inRow->Weight1[iif] > 0.0)  &&  
	  (inRow->Imag1[iif] != fblank)  &&  (inRow->Real1[iif] != fblank)) {
	/* Have solution? */
	if ((insphs[jndx] != fblank)  &&  (gradX[itime] != fblank)) {
	  amp = sqrt (inRow->Real1[iif]*inRow->Real1[iif] + 
		      inRow->Imag1[iif]*inRow->Imag1[iif]);
	  
	  /* Replace phase with model? */
	  faz = insphs[jndx];
	  if (doGrad) faz +=  stnx * gradX[itime] +  stny * gradY[itime];
	  
	  /* Replace phase with residuals? */
	  if (doRes) {
	    indx = ARR_INDEX(itime,antno,iif,0);
	    if (phase[indx] != fblank) {
	      /* Add gradient? */
	      if (doGrad) { /* yes */
		faz = phase[indx] - insphs[jndx];
	      } else {     /* No */
	      faz = phase[indx] - 
		(insphs[jndx] +  gradX[itime] * stnx + gradY[itime] * stny);
	      }
	    } else {
	      faz = fblank;
	    } 
	  } 
	  if (faz != fblank) {
	    outRow->Real1[iif] = amp * cos (faz);
	    outRow->Imag1[iif] = amp * sin (faz);
	  } 
	}
      } /* end if pol1 valid */ 
      
      /* Is previous solution for pol 2 good? */
      if (numpol>1) {
	/* Copy stuff for poln 2 */
	outRow->Delay2[iif]  = inRow->Delay2[iif];
	outRow->Rate2[iif]   = inRow->Rate2[iif];
	outRow->Weight2[iif] = inRow->Weight2[iif];
	outRow->RefAnt2[iif] = inRow->RefAnt2[iif];
	
	/* in case no solution - blank */
	jndx = ARR_INDEX(0,antno,iif,1);
	if (!doBlank && !doRes) {  /*replace blanking unless giving residuals */
	  outRow->Real2[iif]   = cos(insphs[jndx]);
	  outRow->Imag2[iif]   = sin(insphs[jndx]);
	  outRow->Weight2[iif] = 1.0;
	} else { /* blanking */
	  outRow->Real2[iif] = fblank;
	  outRow->Imag2[iif] = fblank;
	}
	
	if ((inRow->Weight2[iif] > 0.0)  &&  
	    (inRow->Imag2[iif] != fblank)  &&  (inRow->Real2[iif] != fblank)) {
	  /* Have solution? */
	  if ((insphs[jndx] != fblank)  &&  (gradX[itime] != fblank)) {
	    amp = sqrt (inRow->Real2[iif]*inRow->Real2[iif] + inRow->Imag2[iif]*inRow->Imag2[iif]);
	    
	    /* Replace phase with model? */
	    faz = insphs[jndx];
	    if (doGrad) faz +=  stnx * gradX[itime] +  stny * gradY[itime];
	    
	    /* Replace phase with residuals? */
	    if (doRes) {
	      indx = ARR_INDEX(itime,antno,iif,1);
	      if (phase[indx] != fblank) {
		/* Add gradient? */
		if (doGrad) { /* yes */
		  faz = phase[indx] - insphs[jndx];
		} else {     /* No */
		  faz = phase[indx] - 
		    (insphs[jndx] +  gradX[itime] * stnx + gradY[itime] * stny);
		}
	      } else {
		faz = fblank;
	      } 
	    } /* end do res */ 
	    if (faz != fblank) {
	      outRow->Real2[iif] = amp * cos (faz);
	      outRow->Imag2[iif] = amp * sin (faz);
	    } 
	  } 
	} /* end if pol2 valid */ 
      } /* end if 2 poln */
    } /* end loop  L60:  */
    
    /* Write */
    outSNRow = -1;
    retCode = ObitTableSNWriteRow (outTab, outSNRow, outRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  } /* end loop  L600: */
  
 cleanup:  /* shutdown */
 /* Close SN tables */
  ObitTableSNClose (inTab, err);
  if (err->error) Obit_traceback_msg (err, routine, inTab->name);
  ObitTableSNClose (outTab, err);
  if (err->error) Obit_traceback_msg (err, routine, outTab->name);

  inRow  = ObitTableSNRowUnref(inRow);  /* delete row object */
  outRow = ObitTableSNRowUnref(outRow); /* delete row object */
} /* end WriteSNTab */

/**
 * Show model and final RMSes
 * \param err     Obit Error stack     
 */
void  ShowModel (ObitErr* err)
{
  olong   ipol, itime, iif, iant, itmp, count, it1, it2;
  ofloat insp, res, sum2, rms, rtemp, rt1, stnx, stny;
  ofloat fblank = ObitMagicF();

  /* Loop over antennas giving instrumental phase */
  Obit_log_error(err, OBIT_InfoErr,"Final Instrumental phases:"); 
  for (iant=0; iant<numant; iant++) { 
    for (iif=0; iif<numIF; iif++) { 
      for (ipol=0; ipol<numpol; ipol++) { 
	insp = insphs[ARR_INDEX(0,iant,iif,ipol)];
	if (insp==fblank) {  /* invalid */
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Ant: %3d IF %3d Poln %2d ins. phase = invalid", 
			 iant+1, iif+1, ipol+1);
	} else {             /* valid */
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Ant: %3d IF %3d Poln %2d ins. phase = %6.2f", 
			 iant+1, iif+1, ipol+1, insp*RAD2DG);
	}
      }
    }
  } /* end antenna loop */
  ObitErrLog(err); /* show any error messages on err */
 
  /* Loop over times determining RMS and showing gradients */
  for (itime=0; itime<numtim; itime++) { /* loop 600 */
    if ((gradX[itime] == fblank)  ||  (gradY[itime] == fblank)  
	|| badtim[itime]) continue;

    /* Clear accumulators */
    count = 0;
    sum2 = 0.0;

    /* Accumulate time residuals */
    for (ipol=0; ipol<numpol; ipol++) { 
      for (iif=0; iif<numIF; iif++) { 
	for (iant=0; iant<numant; iant++) { 
	  stnx = Ant->ANlist[iant]->AntXYZ[0];
	  stny = Ant->ANlist[iant]->AntXYZ[1] ;
	  if ((phase[ARR_INDEX(itime,iant,iif,ipol)] != fblank)  &&   
	      (insphs[ARR_INDEX(0,iant,iif,ipol)] != fblank)) {
	    res = insphs[ARR_INDEX(0,iant,iif,ipol)] +  
	      gradX[itime] * stnx + gradY[itime] * stny - 
	      phase[ARR_INDEX(itime,iant,iif,ipol)];
	    /* Keep in range */
	    if (!doUnwrap) {
	      itmp = res / (2.0*G_PI);
	      res -= itmp * 2.0*G_PI;
	      if (res > G_PI)  res -= 2.0*G_PI;
	      if (res < -G_PI) res += 2.0*G_PI;
	    }
	    /* To degrees */
	    res *= RAD2DG;
	    count++;
	    sum2 += res*res;
	  } 
	} 
      } 
    } 

    /* Time RMS */
    if (count > 6) {
      rms = sqrt (sum2 / count);
    } else {
      rms = 0.0;
    } 
    /* Tell  */
    rtemp = 24.0*sntime[itime];
    it1 = rtemp;
    rtemp = (rtemp - it1)*60.0;
    it2 = rtemp;
    rt1 = (rtemp - it2)*60.0;
    Obit_log_error(err, OBIT_InfoErr, 
		   "%5d time: %3d %3d %6.2f rms=%7.1f grad= %10.2f, %10.2f", 
		   itime+1, it1, it2, rt1, rms, 
		   1.330e5*gradX[itime], 1.330e5*gradY[itime]);
    ObitErrLog(err); /* show any error messages on err */
  } /* end time loop  L600: */
  
} /* end ShowModel */

/**
 * Does a median window smoothing of Gradients
 * \param     myInput   Input parameters on InfoList :
 * \li width Width in min of smoothing of fitted gradients.
 *            This is only useful if doRes=True and doGrad=False.
 *            If width>0.0 the the fitted gradients are smoothed.
 *            The parameter controls the type of smoothing.
 * \li alpha. 0 -> 1 = pure boxcar -> pure MWF (ALPHA of the 
 *            data samples are discarded and the rest averaged). 
 *            This is only useful if doRes=True and doGrad=False.
 */
void  SmooGrad (ObitInfoList* myInput)
{
  ofloat *t=NULL, *wor=NULL, *yor=NULL, *w=NULL, *ys=NULL, *ws=NULL;
  ofloat width=2.0, alpha=0.8, fblank = ObitMagicF();
  olong i, n;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  /* Control parameters */
  ObitInfoListGetTest(myInput, "width", &type, dim, &width);
  width /= 1440.0; /* to days */
  ObitInfoListGetTest(myInput, "alpha", &type, dim, &alpha);

  n = numtim;

  /* Allocate work arrays */
  t   = g_malloc0(n*sizeof(ofloat));
  wor = g_malloc0(n*sizeof(ofloat));
  yor = g_malloc0(n*sizeof(ofloat));
  w   = g_malloc0(n*sizeof(ofloat));
  ys  = g_malloc0(n*sizeof(ofloat));
  ws  = g_malloc0(n*sizeof(ofloat));

  /* Setup for gradX */
  for (i=0; i<n; i++) {
    t[i] = (ofloat)sntime[i];
    if (badtim[i]) w[i] = 0.0;
    else w[i] = 1.0;
  }

  /* Smooth gradX */
  ObitUVSolnSmooMWF (width, alpha, t, gradX, w, n, ys, ws, yor, wor, TRUE);

  /* Save smoothed values for gradX */
  for (i=0; i<n; i++) {
    gradX[i] = ys[i];
    if ((w[i]>0) || (ys[i]==fblank)) badtim[i] = FALSE;
    else badtim[i] = TRUE;
  }

  /* Setup for gradY */
  for (i=0; i<n; i++) {
    if (badtim[i]) w[i] = 0.0;
    else w[i] = 1.0;
  }

  /* Smooth gradY */
  ObitUVSolnSmooMWF (width, alpha, t, gradY, w, n, ys, ws, yor, wor, TRUE);

  /* Save smoothed values for gradY */
  for (i=0; i<n; i++) {
    gradY[i] = ys[i];
    if ((w[i]>0) || (ys[i]==fblank)) badtim[i] = FALSE;
    else badtim[i] = TRUE;
  }

  /* Cleanup */
  if (t)   g_free(t);
  if (wor) g_free(wor);
  if (yor) g_free(yor);
  if (w)   g_free(w);
  if (ys)  g_free(ys);
  if (ws)  g_free(ws);
} /* end SmooGrad */

/**
 * If the VLA, unwrap phases down the VLA arms
 * \param err     Obit Error stack     
 */
void VLAUnwrap (ObitErr* err)
{
  olong i, j, last, prev, ipad=0, fpad=0, nread=0, itime, iif, ipol;
  olong nNorth, nEast, nWest, VLAN[10], VLAE[10], VLAW[10];
  olong iph1=0, iph2=0, ipht=0, jndx1=0, jndx2=0;
  ofloat sum, delta, deltalast, fblank = ObitMagicF();
  gchar pad[9];
  gchar *routine = "VLAUnwrap";

  if (err->error) return;

  /* Ignore if not VLA/EVLA */
  if (!Ant->isVLA) return;

  if (prtLv>=2) Obit_log_error(err, OBIT_InfoErr, "Unwrapping VLA phases");
  ObitErrLog(err); /* show any error messages on err */

  /* Find order on arms */
  for (i=0; i<10; i++) VLAN[i] = VLAE[i] = VLAW[i] = 0;

  /* North Arm */
  last = 1000; prev = 0; nNorth = 0; 
  for (j=0; j<10; j++) { /* outer loop over north arm */
    for (i=0; i<numant; i++) {
      if ((!strncmp(Ant->ANlist[i]->AntName,"VLA:N",5)) ||
	  (!strncmp(Ant->ANlist[i]->AntName,"EVLA:N",6))) { 
	/* On North arm - get number */
	if (!strncmp(Ant->ANlist[i]->AntName,"VLA:N",5)) {
	  strncpy (pad, &Ant->ANlist[i]->AntName[5], 3); pad[3] = 0;
	  nread = sscanf (pad, "%d", &ipad);
	} else if (!strncmp(Ant->ANlist[i]->AntName,"EVLA:N",6)) {
	  strncpy (pad, &Ant->ANlist[i]->AntName[6], 2); pad[2] = 0;
	  nread = sscanf (pad, "%d", &ipad);
	}
	if (nread!=1) {
	  Obit_log_error(err, OBIT_Error, "%s: Error antenna pad number from %s", 
			 routine, Ant->ANlist[i]->AntName);
	  return;
	}
	/* Want smallest pad number greater than prev */
	if ((ipad<last) && (ipad>prev)) {
	  last = ipad;
	  VLAN[nNorth] = i;
	}
      }
    } /* end loop over North arm */;
    prev = last;
    last = 1000; 
    if (prev<100) nNorth++;
  } /* end outer loop over north arm */

  /* East Arm */
  last = 1000; prev = 0; nEast = 0;
  for (j=0; j<10; j++) { /* outer loop over north arm */
    for (i=0; i<numant; i++) {
      if ((!strncmp(Ant->ANlist[i]->AntName,"VLA:E",5)) ||
	  (!strncmp(Ant->ANlist[i]->AntName,"EVLA:E",6))) { 
	/* On North arm - get number */
	if (!strncmp(Ant->ANlist[i]->AntName,"VLA:E",5)) {
	  strncpy (pad, &Ant->ANlist[i]->AntName[5], 3); pad[3] = 0;
	  nread = sscanf (pad, "%d", &ipad);
	} else if (!strncmp(Ant->ANlist[i]->AntName,"EVLA:E",6)) {
	  strncpy (pad, &Ant->ANlist[i]->AntName[6], 2); pad[2] = 0;
	  nread = sscanf (pad, "%d", &ipad);
	}
	if (nread!=1) {
	  Obit_log_error(err, OBIT_Error, "%s: Error antenna pad number from %s", 
			 routine, Ant->ANlist[i]->AntName);
	  return;
	}
	/* Want smallest pad number greater than prev */
	if ((ipad<last) && (ipad>prev)) {
	  last = ipad;
	  VLAE[nEast] = i;
	}
      }
    } /* end loop over east arm */;
    prev = last;
    last = 1000; 
    if (prev<100) nEast++;
  } /* end outer loop over east arm */

  /* West Arm */
  last = 1000; prev = 0; nWest = 0;
  for (j=0; j<10; j++) { /* outer loop over wes6 arm */
    for (i=0; i<numant; i++) {
      if ((!strncmp(Ant->ANlist[i]->AntName,"VLA:W",5)) ||
	  (!strncmp(Ant->ANlist[i]->AntName,"EVLA:W",6))) { 
	/* On North arm - get number */
	if (!strncmp(Ant->ANlist[i]->AntName,"VLA:W",5)) {
	  strncpy (pad, &Ant->ANlist[i]->AntName[5], 3); pad[3] = 0;
	  nread = sscanf (pad, "%d", &ipad);
	} else if (!strncmp(Ant->ANlist[i]->AntName,"EVLA:W",6)) {
	  strncpy (pad, &Ant->ANlist[i]->AntName[6], 2); pad[2] = 0;
	  nread = sscanf (pad, "%d", &ipad);
	}
	if (nread!=1) {
	  Obit_log_error(err, OBIT_Error, "%s: Error antenna pad number from %s", 
			 routine, Ant->ANlist[i]->AntName);
	  return;
	}
	/* Want smallest pad number greater than prev */
	if ((ipad<last) && (ipad>prev)) {
	  last = ipad;
	  VLAW[nWest] = i;
	}
      }
    } /* end loop over west arm */;
    prev = last;
    last = 1000; 
    if (prev<100) nWest++;
  } /* end outer loop over west arm */

  /* Loop over times determining RMS and showing gradients */
  for (ipol=0; ipol<numpol; ipol++) {        /* loop over polarization */
    for (iif=0; iif<numIF; iif++) {          /* loop over IF */
      for (itime=0; itime<numtim; itime++) { /* loop over time */

	/* Loop down north arm unwrapping */
	deltalast = 0.0;
	/* Find first pad with valid data */
	for (ipad=0; ipad<nNorth; ipad++) {
	  ipht = ARR_INDEX(itime,VLAN[ipad],iif,ipol);
	  if (phase[ipht]!=fblank) {
	    fpad = ipad;
	    iph1 = ipht;
	    jndx1 = ARR_INDEX(0,VLAN[ipad],iif,ipol);
	    break;
	  }
	} /* end finding first good pad */
	
	/* Unwrap big turns in the wrong direction */
	sum = 0.0;
	for (ipad=fpad+1; ipad<nNorth; ipad++) {
	  iph2 = ARR_INDEX(itime,VLAN[ipad],iif,ipol);
	  if (phase[iph2]==fblank) continue;  /* This one OK? */
	  if ((phase[iph1]!=fblank) && (phase[iph2]!=fblank)) {
	    jndx2 = ARR_INDEX(0,VLAN[ipad],iif,ipol);
	    delta = phase[iph2] - phase[iph1] - insphs[jndx2] + insphs[jndx1] + sum;
	    /* Assume that the phase gradient has the same sign as the 
	       last one if the differences are > 0.25 rad */
	    if ((deltalast>0.25) && (delta<-0.25)) {
	      sum += 2.0*G_PI;
	    } else if ((deltalast<-0.25) && (delta>0.25)) {
	      sum -= 2.0*G_PI;
	    } else if (delta>G_PI) sum -= 2.0*G_PI;
	    else if (delta<-G_PI  )sum += 2.0*G_PI;
	    phase[iph2] += sum;
	    delta = phase[iph2] - phase[iph1] - insphs[jndx2] + insphs[jndx1];
	    deltalast = delta;
	    iph1  = iph2;
	    jndx1 = jndx2;
	  }
	} /* end loop over north pads */


	/* Loop down East arm unwrapping */
	deltalast = 0.0;
	/* Find first pad with valid data */
	for (ipad=0; ipad<nNorth; ipad++) {
	  ipht = ARR_INDEX(itime,VLAE[ipad],iif,ipol);
	  if (phase[ipht]!=fblank) {
	    fpad = ipad;
	    iph1 = ipht;
	    jndx1 = ARR_INDEX(0,VLAE[ipad],iif,ipol);
	    break;
	  }
	} /* end finding first good pad */
	
	/* Unwrap big turns in the wrong direction */
	sum = 0.0;
	for (ipad=fpad+1; ipad<nNorth; ipad++) {
	  iph2 = ARR_INDEX(itime,VLAE[ipad],iif,ipol);
	  if (phase[iph2]==fblank) continue;  /* This one OK? */
	  if ((phase[iph1]!=fblank) && (phase[iph2]!=fblank)) {
	    jndx2 = ARR_INDEX(0,VLAE[ipad],iif,ipol);
	    delta = phase[iph2] - phase[iph1] - insphs[jndx2] + insphs[jndx1] + sum;
	    /* Assume that the phase gradient has the same sign as the 
	       last one if the differences are > 0.25 rad */
	    if ((deltalast>0.25) && (delta<-0.25)) {
	      sum += 2.0*G_PI;
	    } else if ((deltalast<-0.25) && (delta>0.25)) {
	      sum -= 2.0*G_PI;
	    } else if (delta>G_PI) sum -= 2.0*G_PI;
	    else if (delta<-G_PI  )sum += 2.0*G_PI;
	    phase[iph2] += sum;
	    delta = phase[iph2] - phase[iph1] - insphs[jndx2] + insphs[jndx1];
	    deltalast = delta;
	    iph1  = iph2;
	    jndx1 = jndx2;
	  }
	} /* end loop over east pads */

	/* Loop down West arm unwrapping */
	deltalast = 0.0;
	/* Find first pad with valid data */
	for (ipad=0; ipad<nNorth; ipad++) {
	  ipht = ARR_INDEX(itime,VLAW[ipad],iif,ipol);
	  if (phase[ipht]!=fblank) {
	    fpad = ipad;
	    iph1 = ipht;
	    jndx1 = ARR_INDEX(0,VLAW[ipad],iif,ipol);
	    break;
	  }
	} /* end finding first good pad */
	
	/* Unwrap big turns in the wrong direction */
	sum = 0.0;
	for (ipad=fpad+1; ipad<nNorth; ipad++) {
	  iph2 = ARR_INDEX(itime,VLAW[ipad],iif,ipol);
	  if (phase[iph2]==fblank) continue;  /* This one OK? */
	  if ((phase[iph1]!=fblank) && (phase[iph2]!=fblank)) {
	    jndx2 = ARR_INDEX(0,VLAW[ipad],iif,ipol);
	    delta = phase[iph2] - phase[iph1] - insphs[jndx2] + insphs[jndx1] + sum;
	    /* Assume that the phase gradient has the same sign as the 
	       last one if the differences are > 0.25 rad */
	    if ((deltalast>0.25) && (delta<-0.25)) {
	      sum += 2.0*G_PI;
	    } else if ((deltalast<-0.25) && (delta>0.25)) {
	      sum -= 2.0*G_PI;
	    } else if (delta>G_PI) sum -= 2.0*G_PI;
	    else if (delta<-G_PI  )sum += 2.0*G_PI;
	    phase[iph2] += sum;
	    delta = phase[iph2] - phase[iph1] - insphs[jndx2] + insphs[jndx1];
	    deltalast = delta;
	    iph1  = iph2;
	    jndx1 = jndx2;
	  }
	} /* end loop over west pads */

      } /* End loop over time */
    } /* End loop over IF */
  } /* End loop over poln */
} /* end VLAUnwrap */
