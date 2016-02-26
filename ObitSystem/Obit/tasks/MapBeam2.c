/* $Id: MapBeam2.c  $  */
/* Obit task to Map beam polarization in correlations (RR,LL,LR,RL)   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2015-2016                                          */
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

#include "ObitUVDesc.h"
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
void doChanPoln (gchar *Source, olong *antList, ObitInfoList* myInput, 
		 ObitUV* inData, ObitUV* avgData, ObitErr* err);

/* Image */
void doImage (gchar *Source, gboolean doRMS, gboolean doPhase, olong *antList, 
	      ObitInfoList* myInput, ObitUV* inData, 
	      olong *nchan, olong *nIF, olong *npoln, ObitErr* err);

/* Write history */
void MapBeamHistory (gchar *Source, gchar* Stok, ObitInfoList* myInput, 
		    ObitUV* inData, ObitImage* outImage, ObitErr* err);

/* Average  data */
ObitUV* doAvgData (ObitInfoList *myInput, ObitUV* inData, ObitErr *err);

/* Get list of antennas */
olong* getAntList (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);

/* Accumulate data into lists */
void  accumData (ObitUV* inData, ObitInfoList* myInput, olong ant,
		 olong nchan, olong nIF, olong selem, olong *nelem,
		 ofloat *SumRRr, ofloat *SumRRi, ofloat *SumRRI, ofloat *SumRRWt,
		 ofloat *SumRLr, ofloat *SumRLi, ofloat *SumRLQ, ofloat *SumRLWt,
		 ofloat *SumLRr, ofloat *SumLRi, ofloat *SumLRU, ofloat *SumLRWt,
		 ofloat *SumLLr, ofloat *SumLLi, ofloat *SumLLV, ofloat *SumLLWt,
		 ofloat *SumAzCell, ofloat *SumElCell, ofloat *SumPACell, 
		 olong *CntCell, 
		 ofloat *avgAz, ofloat *avgEl, ofloat *avgPA, 
		 ObitErr* err);

/* Grid data into cells */
void  gridData (ObitInfoList* myInput, olong nchan, olong nIF, olong npoln,
		olong selem, olong nelem,
		ofloat *SumRRr, ofloat *SumRRi, ofloat *SumRRI, ofloat *SumRRWt,
		ofloat *SumRLr, ofloat *SumRLi, ofloat *SumRLQ, ofloat *SumRLWt,
		ofloat *SumLRr, ofloat *SumLRi, ofloat *SumLRU, ofloat *SumLRWt,
		ofloat *SumLLr, ofloat *SumLLi, ofloat *SumLLV, ofloat *SumLLWt,
		ofloat *SumAzCell, ofloat *SumElCell, ObitFArray **grids);

/* Lagrangian interpolation coefficients */
void lagrange(ofloat x, ofloat y, olong n, olong hwid, olong inc,
	      ofloat *xlist, ofloat *ylist, ofloat *wlist, ofloat *coef);

/* Other interpolation coefficients */
void interpolate(ofloat x, ofloat y, olong n, olong hwid, olong inc,
		 ofloat *xlist, ofloat *ylist, ofloat *wlist, ofloat *coef);

/* Normalize and convert to type */
static void normArray(gboolean doRMS, gboolean doPhase, 
		      ObitFArray *arrayr, ObitFArray *arrayi, ObitFArray *arrayw, 
		      ObitFArray *array2, ObitFArray *out);
/* Program globals */
gchar *pgmName = "MapBeam2";     /* Program name */
gchar *infile  = "MapBeam2.in" ; /* File with program inputs */
gchar *outfile = "MapBeam2.out"; /* File to contain program outputs */
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

/*---------------Private structures----------------*/
/* Threaded function argument */
typedef struct {
  ObitThread *thread;  /* ObitThread to use */
  olong loy;           /* first (0-rel) y*/
  olong hiy;           /* laast (0-rel) y*/
  olong nx;            /* Number of columns */
  olong ny;            /* Number of rows */
  olong hwid;          /* Half width of interpolation */
  olong nchan;         /* Number of channels in output */
  olong nIF;	       /* Number of IFs in output  */
  olong npoln;	       /* Number of polarizations in output  */
  olong selem;         /* Size (floats) of list element */
  olong nelem; 	       /* Number of list elements  */
  ofloat xcen;         /* X center cell */
  ofloat ycen;         /* Y center cell */
  ofloat *SumRRr;      /* Real Stokes RR (or XX) accumulation list */
  ofloat *SumRRi;      /* Imag Stokes RR (or XX) accumulation list  */
  ofloat *SumRR2;      /* Stokes RR**2 accumulation list  */
  ofloat *SumRRw;      /* Stokes RR weight accumulation list  */
  ofloat *SumRLr;      /* Real Stokes RL (or XY) accumulation list */
  ofloat *SumRLi;      /* Imag Stokes RL (or XY) accumulation list */
  ofloat *SumRL2;      /* Stokes RL**2 (or XY**2) accumulation list */
  ofloat *SumRLw;      /* Stokes RL weight accumulation list  */
  ofloat *SumLRr;      /* Real Stokes LR (or YX) accumulation list */
  ofloat *SumLRi;      /* Imag Stokes LR (or YX) accumulation list */
  ofloat *SumLR2;      /* Stokes LR**2 (or YX**2) accumulation list  */
  ofloat *SumLRw;      /* Stokes LR weight accumulation list  */
  ofloat *SumLLr;      /* Real Stokes LL (or YY) accumulation list */
  ofloat *SumLLi;      /* Imag Stokes LL (or YY) accumulation list */
  ofloat *SumLL2;      /* Stokes  LL**2 (or YY*2) accumulation list */
  ofloat *SumLLw;      /* Stokes LL weight accumulation list  */
  ofloat *SumAzCell;   /* Azimuth offset accumulation list */
  ofloat *SumElCell;   /* Elevation offset  accumulation list */
  ofloat *coef[4];     /* Work space array */
  ObitFArray **grids;  /* Array of output ObitFArrays
			  fastest to slowest, channel, IF, Stokes 
                          Stokes in order pp,qq,pq,qp */
   olong      ithread;  /* Thread number  */
} MBFuncArg;
		     
/** Private: Make Threaded args */
static olong MakeMBFuncArgs (ObitThread *thread, 
			     olong nchan, olong nIF, olong npoln, olong selem, olong nelem, 
			     olong nx, olong ny, olong hwid, ofloat xcen, ofloat ycen,
			     ofloat *SumRRr, ofloat *SumRRi, ofloat *SumRR2, ofloat *SumRRw,
			     ofloat *SumRLr, ofloat *SumRLi, ofloat *SumRL2, ofloat *SumRLw,
			     ofloat *SumLRr, ofloat *SumLRi, ofloat *SumLR2, ofloat *SumLRw,
			     ofloat *SumLLr, ofloat *SumLLi, ofloat *SumLL2, ofloat *SumLLw,
			     ofloat *SumAzCell, ofloat *SumElCell, ObitFArray **grids,
			     MBFuncArg ***ThreadArgs);

/** Private: Delete Threaded args */
static void KillMBFuncArgs (olong nargs, MBFuncArg **ThreadArgs);

/** Private: Threaded Gridding */
static gpointer ThreadMB2Grid (gpointer arg);

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
/*      iStoke    Stokes number (0-rel), pp, qq, pq, qp (p&q ortho. feeds)*/
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
  olong     i, n, Aseq, disk, cno, axNum, axFNum, axIFNum, nx=11, ny=11;
  olong     jStoke;
  gchar     *Type, *strTemp, outFile[129], *outName, *outF, stemp[32];
  gchar     Aname[13], Aclass[7], *Atype = "MA";
  gchar     *today=NULL;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong     blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong     trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gboolean  exist;
  gchar     tname[129], **chStokes;
  gchar     *chCcor[]={"RR","LL","RL","LR"}, *chLcor[]={"XX","YY","XY","YX"};
  odouble   *StokCrval;
  odouble   CFCrval[4] = {-1.0, -2.0,-3.0,-4.0}, LFCrval[4] = {-5.0, -6.0,-7.0,-8.0};
  gfloat    xCells, yCells;
  gchar     *FITS = "FITS";
  ObitImageDesc *outDesc;
  gchar     *routine = "setOutput";

  /* error checks */
  if (err->error) return outImage;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* correlator types */
  if (inData->myDesc->crval[inData->myDesc->jlocs]<=-3.0) chStokes = chLcor;
  else                                                    chStokes = chCcor;
  if (inData->myDesc->crval[inData->myDesc->jlocs]<=-7.0) StokCrval = LFCrval;
  else                                                    StokCrval = CFCrval;
  /* Which poln? */
  jStoke = iStoke;

  /* Notes:
     1) AIPS Name = Source +outName
            Class = "Santno"/"Sall" S=Stokes, specific or average ant.
        FITS Name = Source+"Santno"/"Sall"+outFile S=Stokes, specific or average an
        for RMS add "R" to end of class or equivalent
        for Phase add "P" to end of class or equivalent
  */

  /* Create basic output Image Object */
  g_snprintf (tname, 100, "output Image %sPol",chStokes[jStoke]);
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
	g_snprintf (tname, 100, "%s%3.3dR", chStokes[jStoke],ant);
      } else {       /* all */
	g_snprintf (tname, 100, "%sAllR", chStokes[jStoke]);
      }
    }  else if (doPhase) {  /* Phase */
      /* One or all antennas? */
      if (ant>0) {   /* One */
	g_snprintf (tname, 100, "%s%3.3dP", chStokes[jStoke],ant);
      } else {       /* all */
	g_snprintf (tname, 100, "%sAllP", chStokes[jStoke]);
      }
    } else {      /* Stokes */
      if (ant>0) {   /* One */
	g_snprintf (tname, 100, "%s%3.3d  ", chStokes[jStoke],ant);
      } else {       /* all */
	g_snprintf (tname, 100, "%sAll  ", chStokes[jStoke]);
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
	g_snprintf (outFile, 128, "%s%s%3.3dRMS%s", stemp,chStokes[jStoke],ant,tname);
      } else {       /* all */
	g_snprintf (outFile, 128, "%s%sAllRMS%s", stemp,chStokes[jStoke],tname);
      }
     } else if (doPhase) {  /* Phase */
      /* One or all antennas? */
      if (ant>0) {   /* One */
	g_snprintf (outFile, 128, "%s%s%3.3dPh%s", stemp,chStokes[jStoke],ant,tname);
      } else {       /* all */
	g_snprintf (outFile, 128, "%s%sAllPh%s", stemp,chStokes[jStoke],tname);
      }
    } else {      /* Stokes */
      if (ant>0) {   /* One */
	g_snprintf (outFile, 128, "%s%s%3.3d%s", stemp,chStokes[jStoke],ant,tname);
      } else {       /* all */
	g_snprintf (outFile, 128, "%s%sAll%s", stemp,chStokes[jStoke],tname);
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
  outImage->myDesc->crval[axNum]  = StokCrval[jStoke];
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
  ObitUV*      avgData = NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        maxlen, isource, failed=0, good=0;
  olong        *antList=NULL;
  gboolean     isBad = FALSE;
  gchar        *Fail="Failed  ", *Done="Done    ";
  gchar        *dataParms[] = {  /* Source selection*/
    "Sources", "souCode", "Qual", "Stokes", "timeRange", "UVRange", "FreqID",
    NULL
  };
  gchar *routine = "doSources";

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Make sure selector set on inData */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Get source list to do */
  doList = ObitUVUtilWhichSources (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Get list of antennas */
  antList = getAntList(myInput, inData, err);
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

    /* Get averaged data */
    avgData = doAvgData (myInput, inData, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    if (err->prtLv>=3)  /* Timing info */
      Obit_log_error(err, OBIT_InfoErr, "Averaged Data");

    /* Process source - loop overr antennas */
    doChanPoln (Source, antList, myInput, inData, avgData,  err);
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


    /* Save processing summary - success or failure? */
    dim[0] = 8; dim[1] = 1;
    if (isBad) 
      ObitInfoListAlwaysPut (myInput, "Status", OBIT_string, dim, Fail);
    else
      ObitInfoListAlwaysPut (myInput, "Status", OBIT_string, dim, Done);
    avgData  = ObitUVUnref(avgData); /* Scratch averaged data */
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
  } /* end source loop */

  doList = ObitSourceListUnref(doList);
  if (antList) g_free(antList);

}  /* end doSources */

/*----------------------------------------------------------------------- */
/*  Loop over frequencies and polarizations for a single source/antenna   */
/*  Accumulates individual antennas if averaging                          */
/*   Input:                                                               */
/*      Source    Name of source being imaged                             */
/*      antList   -1 terminated list of antennas                          */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to image                                         */
/*      avgData   Averaged ObitUV to image                                */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doChanPoln (gchar *Source, olong *antList, ObitInfoList* myInput, 
		 ObitUV* inData, ObitUV* avgData, ObitErr* err)
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        nchan, nIF, npoln;
  gboolean     doRMS=FALSE, doPhase=FALSE;
  gchar        *routine = "doChanPoln";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  if (err->prtLv>=3)  /* Timing info */
    Obit_log_error(err, OBIT_InfoErr, "In %s",routine);

  /* Want RMS or phase images as well? */
  ObitInfoListGetTest(myInput, "doRMS",   &type, dim, &doRMS);
  ObitInfoListGetTest(myInput, "doPhase", &type, dim, &doPhase);

  /* Image data in order (fastest to slowest, channel, IF, Stokes ) */
  doImage (Source, doRMS, doPhase, antList, myInput, avgData, &nchan, &nIF, &npoln,  
	   err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
}  /* end doChanPoln */

/*----------------------------------------------------------------------- */
/*  Image all selected data                                               */
/*  Returns a list of ObitFArrays containing pixel data                   */
/*  fastest to slowest, channel, IF, Stokes                               */
/*   Input:                                                               */
/*      doRMS     If true make images of RMS, else Poln                   */
/*      doPhase   If TRUE, image is Phase, else Amplitude                 */
/*      antList   -1 terminated list of antennas                          */
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
void doImage (gchar *Source, gboolean ldoRMS, gboolean ldoPhase, 
	      olong *antList, ObitInfoList* myInput, ObitUV* inData, 
	      olong *nchan, olong *nIF, olong *npoln, ObitErr* err)
{
  ObitImage    *outImage[4]={NULL,NULL,NULL,NULL};
  ObitFArray   **imgList = NULL;
  ObitFArray   *fatemp = NULL, *work[4]={NULL,NULL,NULL,NULL};
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        nx=11, ny=11, naxis[2], ipos[2], i, ipoln,jpoln,  ichan, iIF;
  olong        selem, nelem, maxelem, size;
  ofloat       *SumRRr=NULL, *SumRRi=NULL, *SumRR2=NULL, *SumRRWt=NULL;
  ofloat       *SumRLr=NULL, *SumRLi=NULL, *SumRL2=NULL, *SumRLWt=NULL;
  ofloat       *SumLRr=NULL, *SumLRi=NULL, *SumLR2=NULL, *SumLRWt=NULL;
  ofloat       *SumLLr=NULL, *SumLLi=NULL, *SumLL2=NULL, *SumLLWt=NULL;
  ofloat       *SumAzCell=NULL, *SumElCell=NULL, *SumPACell=NULL;
  ofloat       *pqCenter, *qpCenter, *ppCenter, *qqCenter, value, val1, val2;
  ofloat       fblank =  ObitMagicF();
  olong        nant, iant, jant, ant, tloop, ntloop, indx, plane[5], *CntCell=NULL;
  gboolean     doRMS, doPhase, avgAnt=TRUE, doPolNorm=FALSE;
  gboolean     isRMS[3] = {FALSE, TRUE, TRUE},isPhase[3] = {FALSE, TRUE, TRUE};
  gchar        *tableList[]={"AIPS FQ","AIPS AN", NULL};
  gchar        *chStokes[]={"  ","RR","LL","RL","LR"};
  gchar *routine = "doImage";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  if (err->prtLv>=3)  /* Timing info */
    Obit_log_error(err, OBIT_InfoErr, "In %s",routine);

  /* Averaging antennas? */
  ObitInfoListGetTest(myInput, "avgAnt", &type, dim, &avgAnt);
  /* How many antennas in antList */
  nant = 0;
  for (iant=0; iant<256; iant++) {
    if (antList[iant]<0) break;
    nant++;
  }

  /* Set up loops for, ampl, ampl+phase, amp+rms and amp+phase+rms */
  ntloop = 1;
  if (ldoRMS)   ntloop++;
  if (ldoPhase) ntloop++;
  if (ldoRMS && ldoPhase)  {isPhase[2] = isRMS[1] = FALSE; }
  if (!ldoRMS)   {isRMS[1]   = isRMS[2]   = FALSE;}
  if (!ldoPhase) {isPhase[1] = isPhase[2] = FALSE;}

  /*  Number of channels */
  *nchan = inData->myDesc->inaxes[inData->myDesc->jlocf];
  /* Number of IFs */
  if (inData->myDesc->jlocif>=0)
    *nIF  = inData->myDesc->inaxes[inData->myDesc->jlocif];
  else *nIF = 1;
  /* Number of polarizations */
  *npoln = inData->myDesc->inaxes[inData->myDesc->jlocs];
  /* Normalize cross pos by I Pol? */
  ObitInfoListGetTest(myInput, "doPolNorm", &type, dim, &doPolNorm);

  /* How big an image? */
  ObitInfoListGetTest(myInput, "nx", &type, dim, &nx);
  ObitInfoListGetTest(myInput, "ny", &type, dim, &ny);
  naxis[0] = nx; naxis[1] = ny;

  /* Make accumulation lists */
  selem = (*nIF) * (*nchan); /* size of list element */
  nelem = nx * ny;           /* Number of list elements */
  nelem *= 2;                /* To be sure big enough */
  maxelem = nelem;
  size = selem * nelem;      /* size of list */
  SumRRr     = g_malloc0(size*sizeof(ofloat));
  SumRRi     = g_malloc0(size*sizeof(ofloat));
  SumRR2     = g_malloc0(size*sizeof(ofloat));
  SumRRWt    = g_malloc0(size*sizeof(ofloat));
  if (*npoln>=2) {
    SumLLr     = g_malloc0(size*sizeof(ofloat));
    SumLLi     = g_malloc0(size*sizeof(ofloat));
    SumLLWt    = g_malloc0(size*sizeof(ofloat));
    SumLL2     = g_malloc0(size*sizeof(ofloat));
  }
  if (*npoln>=3) {
    SumRLr     = g_malloc0(size*sizeof(ofloat));
    SumRLi     = g_malloc0(size*sizeof(ofloat));
    SumRL2     = g_malloc0(size*sizeof(ofloat));
    SumRLWt    = g_malloc0(size*sizeof(ofloat));
  }
  if (*npoln>=4) {
    SumLRr     = g_malloc0(size*sizeof(ofloat));
    SumLRi     = g_malloc0(size*sizeof(ofloat));
    SumLR2     = g_malloc0(size*sizeof(ofloat));
    SumLRWt    = g_malloc0(size*sizeof(ofloat));
  }
  SumAzCell = g_malloc0(nelem*sizeof(ofloat));  
  SumElCell = g_malloc0(nelem*sizeof(ofloat));
  SumPACell = g_malloc0(nelem*sizeof(ofloat));
  CntCell   = g_malloc0(nelem*sizeof(olong));

  /* Loop over antennas */
  for (iant=0; iant<nant; iant++) {
    ant = antList[iant];  /* antenna to accumulate/grid */
    nelem = maxelem; 

    /* Create beam image real, imaginary, wt, ^2 sets */
    if (imgList==NULL) {
      imgList = g_malloc0(4*(*nchan)*(*nIF)*(*npoln)*sizeof(ObitFArray*));
      i = 0;
      for (ipoln=0; ipoln<(*npoln); ipoln++) {
	work[ipoln] = ObitFArrayCreate (NULL, 2, naxis);  /* poln work */
	for (iIF=0; iIF<(*nIF); iIF++) {
	  for (ichan=0; ichan<(*nchan); ichan++) {
	    imgList[i++] = ObitFArrayCreate (NULL, 2, naxis);  /* real */
	    imgList[i++] = ObitFArrayCreate (NULL, 2, naxis);  /* Imag */
	    imgList[i++] = ObitFArrayCreate (NULL, 2, naxis);  /* Wt   */
	    imgList[i++] = ObitFArrayCreate (NULL, 2, naxis);  /* Amp^2  */
	  }
	}
      }
    } /* end create accumulation images */

    /* Zero antenna accumulators */
    size = selem * nelem;      /* size of list */
    memset (SumRRr,    0, size*sizeof(ofloat));
    memset (SumRRi,    0, size*sizeof(ofloat));
    memset (SumRR2,    0, size*sizeof(ofloat));
    memset (SumRRWt,   0, size*sizeof(ofloat));
    memset (SumAzCell, 0, nelem*sizeof(ofloat));
    memset (SumElCell, 0, nelem*sizeof(ofloat));
    memset (SumPACell, 0, nelem*sizeof(ofloat));
    memset (CntCell,   0, nelem*sizeof(olong));
    if (*npoln>=2) {
      memset (SumLLr , 0, size*sizeof(ofloat));
      memset (SumLLi,  0, size*sizeof(ofloat));
      memset (SumLL2,  0, size*sizeof(ofloat));
      memset (SumLLWt, 0, size*sizeof(ofloat));
    }
    if (*npoln>=3) {
      memset (SumRLr , 0, size*sizeof(ofloat));
      memset (SumRLi,  0, size*sizeof(ofloat));
      memset (SumRL2,  0, size*sizeof(ofloat));
      memset (SumRLWt, 0, size*sizeof(ofloat));
      memset (SumLRr , 0, size*sizeof(ofloat));
      memset (SumLRi,  0, size*sizeof(ofloat));
      memset (SumLR2,  0, size*sizeof(ofloat));
      memset (SumLRWt, 0, size*sizeof(ofloat));
    }
    /* Accumulate data onto lists */
    accumData (inData, myInput, ant, *nchan, *nIF, selem, &nelem,
	       SumRRr, SumRRi, SumRR2, SumRRWt, SumRLr, SumRLi, SumRL2, SumRLWt,
	       SumLRr, SumLRi, SumLR2, SumLRWt, SumLLr, SumLLi, SumLL2, SumLLWt,
	       SumAzCell, SumElCell, SumPACell, CntCell, 
	       &avgAz, &avgEl, &avgPA, err);
    if (err->error) goto cleanup;
    
    /* Grid data */
    if (err->prtLv>=3)  /* Timing info */
      Obit_log_error(err, OBIT_InfoErr, "Grid data antenna %d",ant);
    gridData (myInput, *nchan, *nIF, *npoln, selem, nelem, 
	      SumRRr, SumRRi, SumRR2, SumRRWt, SumRLr, SumRLi, SumRL2, SumRLWt,
	      SumLRr, SumLRi, SumLR2, SumLRWt, SumLLr, SumLLi, SumLL2, SumLLWt,
	      SumAzCell, SumElCell, imgList);

    /* If individual antennas or average on last */
    if (!avgAnt || (iant==(nant-1))) {
 
      /* Loop over amp, phase, RMS */
      for (tloop=0; tloop<ntloop; tloop++) {
	doRMS   = isRMS[tloop];
	doPhase = isPhase[tloop];
	/* Create output for each poln */
	for (ipoln=0; ipoln<(*npoln); ipoln++) {
	  if (avgAnt) jant = 0;
	  else        jant = ant;
	  /* Trap LL or YY only */
	  if ((ipoln==0) && 
	      ((fabs(inData->myDesc->crval[inData->myDesc->jlocs]+2.0)<0.1) ||
	       ((fabs(inData->myDesc->crval[inData->myDesc->jlocs]+6.0)<0.1)))) 
	        jpoln = 1;
	  else  jpoln = ipoln;
	  outImage[ipoln] = setOutput (Source, jpoln, jant, doRMS, doPhase, 
				       myInput, inData, err);
	  if (err->error) goto cleanup;
	  ObitImageOpen (outImage[ipoln], OBIT_IO_WriteOnly, err);
	  /* Save average geometry */
	  dim[0] = dim[1] = dim[2] = 1;
	  ObitInfoListAlwaysPut (outImage[ipoln]->myDesc->info, "avgAz", OBIT_float, dim, &avgAz);
	  ObitInfoListAlwaysPut (outImage[ipoln]->myDesc->info, "avgEl", OBIT_float, dim, &avgEl);
	  ObitInfoListAlwaysPut (outImage[ipoln]->myDesc->info, "avgPA", OBIT_float, dim, &avgPA);
	  ObitImageClose (outImage[ipoln], err);
	  if (err->error) goto cleanup;
	} /* end loop creating output images */
      
	/* Looping over frequency */
	for (iIF=0; iIF<(*nIF); iIF++) {
	  for (ichan=0; ichan<(*nchan); ichan++) {
	    /* Normalize, convert */
	    indx = 4*(0*selem + iIF*(*nchan) + ichan);  /* pp */
	    normArray (doRMS, doPhase, imgList[indx], imgList[indx+1], 
		       imgList[indx+2], imgList[indx+3], work[0]);
	    if (*npoln>=2) {
	      indx = 4*(1*selem + iIF*(*nchan) + ichan);  /* qq */
	      normArray (doRMS, doPhase, imgList[indx], imgList[indx+1], 
			 imgList[indx+2], imgList[indx+3], work[1]);
	    }
	    if (*npoln>=3) {
	      indx = 4*(2*selem + iIF*(*nchan) + ichan);  /* pq */
	      normArray (doRMS, doPhase, imgList[indx], imgList[indx+1], 
			 imgList[indx+2], imgList[indx+3], work[2]);
	    }
	    if (*npoln>=4) {
	      indx = 4*(3*selem + iIF*(*nchan) + ichan);  /* qp */
	      normArray (doRMS, doPhase, imgList[indx], imgList[indx+1], 
			 imgList[indx+2], imgList[indx+3], work[3]);
	    }
	    
	    /* amplitude cross hands: Subtract center RL/XY, LR/YX (source poln) * 0.5(RR/XX+LL/YY) beam from rest */
	    if ((*npoln>=3) && (tloop==0)) {
	      fatemp = ObitFArrayCreate (NULL, 2, naxis);  /* Work array */
	      ipos[0] = nx/2; ipos[1] = ny/2;
	      /* pq (RL or XY) */ 
	      ppCenter = ObitFArrayIndex (work[0], ipos);
	      qqCenter = ObitFArrayIndex (work[1], ipos);
	      pqCenter = ObitFArrayIndex (work[2], ipos);
	      qpCenter = ObitFArrayIndex (work[3], ipos);
	      ObitFArrayAdd (work[1], work[2], fatemp);            /* Ipol = 0.5*(qq+pp) */
	      value = -(*pqCenter) / (0.5*(*ppCenter + *qqCenter));  /* Normalize by IPol beam */
	      ObitFArraySMul (fatemp, value);  /* pq * IPol beam */
	      ObitFArrayAdd (work[2], fatemp, work[2]);
	      /* qp  (LR or YX) */ 
	      ObitFArrayAdd (work[1], work[2], fatemp);            /* Ipol = 0.5*(qq+pp) */
	      value = -(*qpCenter) / (0.5*(*ppCenter + *qqCenter));  /* Normalize by IPol beam */
	      ObitFArraySMul (fatemp, value);  /* qp * IPol beam */
	      ObitFArrayAdd (work[3], fatemp, work[3]);
	      
	      /* Divide pq, qp by I?  */
	      if (doPolNorm) {
		/* Make Stokes I from RR/XX + LL/YY */
		ObitFArrayAdd (work[1], work[2], fatemp);            /* Ipol = 0.5*(qq+pp) */
		value = 0.5;
		ObitFArraySMul (fatemp, value);
		ObitFArrayDiv (work[2], fatemp, work[2]);
		ObitFArrayDiv (work[3], fatemp, work[3]);
	      } else {  /* Otherwise normalize by center I */
		ipos[0] = nx/2; ipos[1] = ny/2;
		ppCenter = ObitFArrayIndex (work[0],  ipos);
		qqCenter = ObitFArrayIndex (work[1],  ipos);
		value = 0.5*(*ppCenter + *qqCenter);
		if ((*ppCenter!=fblank) && (*qqCenter!=fblank)) {
		  if (value!=0.0) value = 1.0 / value;
		  ObitFArraySMul (work[2], value);
		  ObitFArraySMul (work[3], value);
		}
	      } /* End normalize x pol by center */
	      fatemp = ObitFArrayUnref (fatemp);  /* Work array */
	    } /* end if xpol */
	    
	    /* Normalize parallel hand correlations by center */
	    if (tloop==0) {  /* Only amplitude */
	      ipos[0] = nx/2; ipos[1] = ny/2;
	      ppCenter = ObitFArrayIndex (work[0],  ipos);
	      val1 = *ppCenter;
	      if (val1!=fblank) {
		if (val1!=0.0) val1 = 1.0 / val1;
		ObitFArraySMul (work[0], val1);
	      }
	      if (*npoln>=2) {
		qqCenter = ObitFArrayIndex (work[1],  ipos);
		val2 = *qqCenter;
		if (val2!=fblank) {
		  if (val2!=0.0) val2 = 1.0 / val2;
		  ObitFArraySMul (work[1], val2);
		}
	      }
	      /* DEBUG IF 12 chan 2  pixel 31,33
	      if ((iIF==11) && (ichan==1)) {
		fprintf (stderr,"write pp norm %9.6f 31,33 value %9.5f\n",
			 val1, work[0]->array[32*nx+30] );
		fprintf (stderr,"write qq norm %9.5f 31,33 value %9.5f\n",
			 val2, work[1]->array[32*nx+30] );
	      } */
	    } /* end if amp */
	    
	    /* Write images */
	    for (i=0; i<5; i++) plane[i] = 1;
	    i = 0;
	    for (ipoln=0; ipoln<(*npoln); ipoln++) {
	      /* If there is only one channel and multiple IFs, IF axis is first - else Freq */
	      if ((*nchan==1) && (*nIF>1)) { /* IF first */
		plane[1] = ichan+1;
		plane[0] = iIF+1;
		ObitImagePutPlane (outImage[ipoln], work[ipoln]->array, plane,  err);
		if (err->error) goto cleanup;
	      } else { /* Freq first */
		plane[1] = iIF+1;
		plane[0] = ichan+1;
		ObitImagePutPlane (outImage[ipoln], work[ipoln]->array, plane,  err);
		if (err->error) goto cleanup;
	      }
	    } /* end stokes loop */
	  } /* end Channel loop */
	} /* end IF loop */
	for (ipoln=0; ipoln<(*npoln); ipoln++) {
	  /* Copy FQ & AN tables */
	  ObitDataCopyTables ((ObitData*)inData, (ObitData*)outImage[ipoln], NULL,
			      tableList, err);
	  /* History */
	  MapBeamHistory (Source, chStokes[ipoln+1], myInput, inData, 
			  outImage[ipoln], err);
	  if (err->error) goto cleanup;
	} /* end poln loop */
      } /* end data type loop */
      
      /* Cleanup */
      if (imgList) {
	i = 0;
	for (ipoln=0; ipoln<(*npoln); ipoln++) {
	  work[ipoln] = ObitFArrayUnref(work[ipoln]);  /* poln work */
	  outImage[ipoln] = ObitImageUnref(outImage[ipoln]);
	  for (iIF=0; iIF<*nIF; iIF++) {
	    for (ichan=0; ichan<*nchan; ichan++) {
	      ObitFArrayUnref(imgList[i++]);  /* Real */
	      ObitFArrayUnref(imgList[i++]);  /* Imag */
	      ObitFArrayUnref(imgList[i++]);  /* Wt  */
	      ObitFArrayUnref(imgList[i++]);  /* Amp^2  */
	    }
	  }
	}
	g_free(imgList); imgList = NULL;
      } /* end cleanup image list */
    
    } /* end if last */
  } /* end antenna loop */
  /* Cleanup */
 cleanup:
  if (SumRRr)     g_free(SumRRr);
  if (SumRRi)     g_free(SumRRi);
  if (SumRR2)     g_free(SumRR2);
  if (SumRRWt)    g_free(SumRRWt);
  if (SumRLr)     g_free(SumRLr);
  if (SumRLi)     g_free(SumRLi);
  if (SumRL2)     g_free(SumRL2);
  if (SumRLWt)    g_free(SumRLWt);
  if (SumLRr)     g_free(SumLRr);
  if (SumLRi)     g_free(SumLRi);
  if (SumLR2)     g_free(SumLR2);
  if (SumLRWt)    g_free(SumLRWt);
  if (SumLLr)     g_free(SumLLr);
  if (SumLLi)     g_free(SumLLi);
  if (SumLL2)     g_free(SumLL2);
  if (SumLLWt)    g_free(SumLLWt);
  if (SumAzCell)  g_free(SumAzCell);
  if (SumElCell)  g_free(SumElCell);
  if (CntCell)    g_free(CntCell);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  return;
} /* end doImage */

/*----------------------------------------------------------------------- */
/*  Write History for MapBeam                                             */
/*   Input:                                                               */
/*      Source    Name of source being imaged                             */
/*      Stoke     Stokes's parameter imaged RR,LL,RL,LR                   */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outImage  ObitImage to write history to                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void MapBeamHistory (gchar *Source, gchar* Stoke, ObitInfoList* myInput, 
		    ObitUV* inData, ObitImage* outImage, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "Sources", "Qual", "souCode", "timeRange",  "Stokes",
    "FreqID", "BIF", "EIF", "BChan", "EChan",  
    "doCalib", "gainUse", "doPol",  "PDVer", "flagVer", "doBand ",  "BPVer", "Smooth",
    "outDType", "outFile",  "outDisk", "outName", "outSeq", "sclAzEl",
    "nx", "ny", "xCells", "yCells", "hwid", "avgTime", "avgFreq", "chAvg", "ChanSel",
    "blnkTime", "avgAnt", "doRMS", "doPhase", "doPolNorm", "Antennas", "RefAnts",
    "OffAz", "OffEl",
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
  g_snprintf (hicard, 80, "%s Source = '%s', Stokes= '%s'", pgmName, Source, Stoke);
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
  gboolean     isScratch, doAvgAll, noScale=TRUE;
  /*gchar        Stokes[5];*/
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "UVRange", "timeRange", "doCalSelect", 
    "BIF", "EIF", "BChan", "EChan","subA", "Antennas", "Stokes",
    "doCalib", "gainUse", "doBand", "BPVer", "Smooth", "flagVer", 
    "doPol", "PDVer",  "avgTime", "avgFreq", "ChanSel", 
    NULL
  };
  gchar        *routine= "doAvgData";

  if (err->error) return avgData;  /* Prior error? */

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, inData->name, avgData);
  
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
    dim[0] = dim[1] = 1;
    ObitInfoListAlwaysPut (inData->info, "noScale", OBIT_bool, dim, &noScale);
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
    dim[0] = dim[1] = 1;
    ObitInfoListAlwaysPut (inData->info, "noScale", OBIT_bool, dim, &noScale);
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
/*        Antennas  Desired antennas                                      */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return:                                                              */
/*     -1 terminated list, 0=> average antennas, g_free when done         */
/*----------------------------------------------------------------------- */
olong* getAntList (ObitInfoList* myInput, ObitUV *inData, ObitErr* err)
{
  olong        *out=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM];
  gboolean     want;
  olong        i, j, iout, count, *Antennas, *RefAnts, nant, nref;
  gchar        *routine= "getAntList";

  if (err->error) return out;  /* Prior error? */
  /* reference antennas */
  ObitInfoListGetP(myInput, "RefAnts",  &type, dim, (gpointer)&RefAnts);
  nref = dim[0];

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
    out  = g_malloc0((dim[0]+1)*sizeof(olong));
    iout = 0;
    nant = MIN (dim[0], inData->myDesc->numAnt[0]);
    for (i=0; i<nant; i++) {
      /* Want? */
      want = TRUE;
      for (j=0; j<nref; j++) if ((i+1)==RefAnts[j]) want = FALSE;
      if (want) out[iout++] = i+1;
    }
    out[iout] = -1;
  } else { /* explicit list */
    out  = g_malloc0((count+1)*sizeof(olong));
    iout = 0;
    for (i=0; i<count; i++) {
      /* Want? */
      want = TRUE;
      for (j=0; j<nref; j++) if (Antennas[i]==RefAnts[j]) want = FALSE;
      if (want) out[iout++] = Antennas[i];
    }
    out[iout] = -1;
  }

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
/*      SumRRr     Real Stokes RR (or XX) accumulation list               */
/*      SumRRi     Imag Stokes  RR (or XX) accumulation list              */
/*      SumRR2     Stokes  RR (or XX) accumulation list                   */
/*      SumRRWt    RR (or XX) Weight accumulation list                   */
/*      SumRLr     Real Stokes RL (or XY) accumulation list               */
/*      SumRLi     Imag Stokes RL (or XY) accumulation list               */
/*      SumRL2     Stokes RL**2 (or XY**2) accumulation list              */
/*      SumRLWt    RL (or XY) Weight accumulation list                    */
/*      SumLRr     Real Stokes LR (or YX) accumulation list               */
/*      SumLRi     Imag Stokes LR (or YX) accumulation list               */
/*      SumLR2     Stokes LR**2 (or YX**2) accumulation list              */
/*      SumLRWt    LR (or YX) Weight accumulation list                    */
/*      SumLLr     Real Stokes LL (or YY) accumulation list               */
/*      SumLLi     Imag Stokes LL (or YY) accumulation list               */
/*      SumLL2     Stokes LL**2 (or YY**2) accumulation list              */
/*      SumLLWt    LL (or YY) Weight accumulation list                    */
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
		 ofloat *SumRRr, ofloat *SumRRi, ofloat *SumRR2, ofloat *SumRRWt,
		 ofloat *SumRLr, ofloat *SumRLi, ofloat *SumRL2, ofloat *SumRLWt,
		 ofloat *SumLRr, ofloat *SumLRi, ofloat *SumLR2, ofloat *SumLRWt,
		 ofloat *SumLLr, ofloat *SumLLi, ofloat *SumLL2, ofloat *SumLLWt,
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
  ofloat   xCells=1.0, yCells=1.0, sclAzEl=1.0, blnkTime=0.0, tblank=-1.0e20;
  ofloat   u, v, time, ulast, vlast, tol, Az, El, PA, *farr;
  ofloat   ss, cc, xr, xi, *OffAz=NULL, *OffEl=NULL, fblank =  ObitMagicF();
  odouble  sumAz, sumEl, sumPA;
  olong    count, maxElem=*nelem, iElem, indx, iant, ant1, ant2, suba, off=0, iver;
  olong    i, j, jlocs, jlocf, jlocif, incs, incf, incif, doff, ddoff, hwid=0;
  olong    nx, ny, iIF, ichan, *refAnts, nRefAnt, ix, iy, prtLv=0;
  gboolean OK1, OK2, isCirc, doPolNorm=FALSE;
  gchar    *routine = "accumData";
  /* olong soff, uc, vc;  DEBUG */

  /* error checks */
  if (err->error) return;

  if (err->prtLv>=3)  /* Timing info */
    Obit_log_error(err, OBIT_InfoErr, "In %s",routine);

  /* Get control parameters */
  ObitInfoListGetTest(myInput, "xCells", &type,   dim, &xCells);
  ObitInfoListGetTest(myInput, "yCells", &type,   dim, &yCells);
  if (yCells==0.0) yCells = xCells;
  /* Cell spacing to radians */
  xCells = (xCells / 3600.0) * DG2RAD;
  yCells = (yCells / 3600.0) * DG2RAD;
  ObitInfoListGetTest(myInput, "sclAzEl", &type,   dim, &sclAzEl);
  if (fabs(sclAzEl)<0.001) sclAzEl = 1.0;
  ObitInfoListGetTest(myInput, "blnkTime", &type, dim, &blnkTime);
  ObitInfoListGetTest(myInput, "prtLv",  &type, dim, &prtLv);
  /* How big an image? */
  ObitInfoListGetTest(myInput, "nx", &type, dim, &nx);
  ObitInfoListGetTest(myInput, "ny", &type, dim, &ny);
  ObitInfoListGetTest(myInput, "hwid", &type, dim, &hwid);
  /* Normalize cross pos by I Pol? */
  ObitInfoListGetTest(myInput, "doPolNorm", &type, dim, &doPolNorm);
  /* Blanking time to days */
  blnkTime /= 86400.0;
  isCirc = (inData->myDesc->crval[inData->myDesc->jlocs]>=-4.0); /* Circular feeds? */
  /* Pointing offsets */
  if (ObitInfoListGetP(myInput, "OffAz",  &type, dim, (gpointer)&farr)) {
    OffAz = g_malloc0(dim[0]*sizeof(ofloat));
    for (i=0; i<dim[0]; i++) OffAz[i] = DG2RAD * farr[i]/60.0;  /* to radians */
  }
  if (ObitInfoListGetP(myInput, "OffEl",  &type, dim, (gpointer)&farr)) {
    OffEl = g_malloc0(dim[0]*sizeof(ofloat));
    for (i=0; i<dim[0]; i++) OffEl[i] = DG2RAD * farr[i]/60.0;  /* to radians */
  }
  
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
  /* Set tolerance if using closest cell or interpolation */
  tol = 0.7 * xCells;  /* Tolerance 0.7 cells  */

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
      /* Sign of u seems OK */
      u    = sclAzEl*inData->buffer[indx+inData->myDesc->ilocu];
      v    = sclAzEl*inData->buffer[indx+inData->myDesc->ilocv];

      /* In blanked time? */
      if (time<tblank) goto next;

      /* Want antennas? */
      ObitUVDescGetAnts(inData->myDesc, &inData->buffer[indx], &ant1, &ant2, &suba);
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

      /* Add pointing offsets */
      if (OffAz) {
	u += OffAz[ant1-1];
	u += OffAz[ant2-1];
      }
      if (OffEl) {
	v += OffEl[ant1-1];
	v += OffEl[ant2-1];
      }

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
	  /* RR/XX */
	  if (inData->buffer[ddoff+2]>0.0) {
	    SumRRr[off]  += inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumRRi[off]  += inData->buffer[ddoff+1]*inData->buffer[ddoff+2];
	    SumRR2[off]  += inData->buffer[ddoff]*inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumRRWt[off] += inData->buffer[ddoff+2]; 
	    /* DEBUG show offset x=-2, y=0
	    if (err->prtLv>4) {
	      if (u>0.0) uc = (olong)(u/xCells+0.5); 
	      else       uc = (olong)(u/xCells-0.5);
	      if (v>0.0) vc = (olong)(v/yCells+0.5);
	      else       vc = (olong)(v/yCells-0.5);
	      if ((uc>=-3) && (uc<=-1) && (vc>=-1) && (vc<=1) && (iIF==11) && (ichan==1)) {
		soff = off;
		fprintf (stderr,"bl %3d-%3d elem %4d  cx %7.2f cy %7.2f vis %9.6f %9.6f wt %5.2f accum %9.6f\n",
			 ant1, ant2, iElem, u/xCells, v/yCells, 
			 inData->buffer[ddoff], inData->buffer[ddoff+1], inData->buffer[ddoff+2], 
			 SumRRr[soff]/SumRRWt[soff]);
	      } 
	    } *//* end debug */
	  }
	  /* RL/XY */
	  ddoff = doff + 2*incs;
	  if (SumRLr && inData->buffer[ddoff+2]>0.0) {
	    SumRLr[off]  += inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumRLi[off]  += inData->buffer[ddoff+1]*inData->buffer[ddoff+2];
	    SumRL2[off]  += inData->buffer[ddoff]*inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumRLWt[off] += inData->buffer[ddoff+2]; 
	  }
	  /* LR/YX */
	  ddoff = doff + 3*incs;
	  if (SumLRr && inData->buffer[ddoff+2]>0.0) {
	    SumLRr[off]  += inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumLRi[off]  += inData->buffer[ddoff+1]*inData->buffer[ddoff+2];
	    SumLR2[off]  += inData->buffer[ddoff]*inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumLRWt[off] += inData->buffer[ddoff+2]; 
	  }
	  /* LL/YY */
	  ddoff = doff + incs;
	  if (SumLLr && inData->buffer[ddoff+2]>0.0) {
	    SumLLr[off]   += inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumLLi[off]   += inData->buffer[ddoff+1]*inData->buffer[ddoff+2];
	    SumLL2[off]   += inData->buffer[ddoff]*inData->buffer[ddoff]*inData->buffer[ddoff+2];
	    SumLLWt[off]  += inData->buffer[ddoff+2]; 
	  }
	  /* DEBUG show offset x=0, y=-3
	  if ((iIF==11) && (ichan==0) && (fabs((u/xCells))<0.8) && (fabs((v/yCells)-3.0)<0.8)) {
	    fprintf (stderr,"bl %2d-%2d elem %4d cx %7.2f cy %7.2f RR %8.5f %8.5f %7.5f LL %8.5f %8.5f %7.5f V %8.5f %8.5f sum R %8.5f sum L %8.5f \n",
		     ant1, ant2,iElem, u/xCells, v/yCells, 
		     inData->buffer[doff], inData->buffer[doff+1], inData->buffer[doff+2],
		     inData->buffer[doff+incs], inData->buffer[doff+incs+1],inData->buffer[doff+incs+2],
		     inData->buffer[doff]-inData->buffer[doff+incs],
		     inData->buffer[doff+1]-inData->buffer[doff+incs+1],
		     SumRRr[off]/SumRRWt[off], SumLLr[off]/SumLLWt[off]);
	  }  *//* end DEBUG */
	} /* end channel loop */
      } /* end IF loop */

      /* update data pointers */
    next:
      indx  += inData->myDesc->lrec;
    } /* end loop over buffer */
  } /* end loop over file */
  
    /* Close */
  retCode = ObitUVClose (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* DEBUG 
  if (err->prtLv>4) {
    fprintf (stderr," RR cell % d accum %9.6f %9.6f wt %7.3f\n",
	     soff, SumRRr[soff]/ SumRRWt[soff],  SumRRi[soff]/ SumRRWt[soff],  SumRRWt[soff]);
  }*/
    
  /* Number of actual elements */
  *nelem = iElem+1;
  Obit_log_error(err, OBIT_InfoErr, "Found antenna %d data in %d cells",ant,*nelem);

  /* Better have something */
  if (iElem<2) return;
 
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
	if (SumRRWt[off]>0) {
	  xr = SumRRr[off] / SumRRWt[off];
	  SumRRr[off] /= SumRRWt[off];
	  SumRRi[off] /= SumRRWt[off];
	  SumRR2[off] = (((SumRR2[off]/SumRRWt[off]) - xr*xr))/SumRRWt[off];
	  if (SumRR2[off]>0.0) SumRR2[off] = sqrt(SumRR2[off]);
	  else                 SumRR2[off] = 0.0;
	  SumRR2[off] /= fabs(xr);  /* Normalize variance by RR/XX */
	} else {
	  SumRRr[off] = fblank;
	  SumRRi[off] = fblank;
	  SumRR2[off] = fblank;
	}
	if (SumRLWt) {
	  if (SumRLWt[off]>0) {
	    xi = SumRLr[off] / SumRLWt[off];
	    SumRLr[off] /= SumRLWt[off];
	    SumRLi[off] /= SumRLWt[off];
	    SumRL2[off] = (((SumRL2[off]/SumRLWt[off]) - xi*xi))/SumRLWt[off];
	    if (SumRL2[off]>0.0) SumRL2[off] = sqrt(SumRL2[off]);
	    else                 SumRL2[off] = 0.0;
	    if (doPolNorm) SumRL2[off] /= fabs(xr);  /* Normalize variance by RR/XX */
	  } else {
	    SumRLr[off] = fblank;
	    SumRLi[off] = fblank;
	    SumRL2[off] = fblank;
	  }
	}
	if (SumLRWt) {
	  if (SumLRWt[off]>0) {
	    xi = SumLRr[off] / SumLRWt[off];
	    SumLRr[off] /= SumLRWt[off];
	    SumLRi[off] /= SumLRWt[off];
	    SumLR2[off] = (((SumLR2[off]/SumLRWt[off]) - xi*xi))/SumLRWt[off];
	    if (SumLR2[off]>0.0) SumLR2[off] = sqrt(SumLR2[off]);
	    else                 SumLR2[off] = 0.0;
	    if (doPolNorm) SumLR2[off] /= fabs(xr);  /* Normalize variance by RR/XX */
	  } else {
	    SumLRr[off] = fblank;
	    SumLRi[off] = fblank;
	    SumLR2[off] = fblank;
	  }
	}
	if (SumLLWt) {
	  if (SumLLWt[off]>0) {
	    xi = SumLLr[off] / SumLLWt[off];
	    SumLLr[off] /= SumLLWt[off];
	    SumLLi[off] /= SumLLWt[off];
	    SumLL2[off] = (((SumLL2[off]/SumLLWt[off]) - xi*xi))/SumLLWt[off];
	    if (SumLL2[off]>0.0) SumLL2[off] = sqrt(SumLL2[off]);
	    else                 SumLL2[off] = 0.0;
	    SumLL2[off] /= fabs(xr);  /* Normalize variance by RR/XX */
	  } else {
	    SumLLr[off] = fblank;
	    SumLLi[off] = fblank;
	    SumLL2[off] = fblank;
	  }
	}
	/* Counter rotate pq, qp if circular feeds for parallactic angle */
	if (isCirc && SumRLr && (SumRLr[off]!=fblank) && (SumLRr[off]!=fblank)) {
	  cc = cos(SumPACell[i]);
	  ss = sin(SumPACell[i]);
	  /* Not sure this right */
	  xr = SumRLr[off];
	  xi = SumRLi[off];
	  SumRLr[off] = cc*xr - ss*xi;
	  SumLRr[off] = cc*xi + ss*xr;
	  xr = SumLRr[off];
	  xi = SumLRi[off];
	  SumRLi[off] = cc*xr + ss*xi;
	  SumLRi[off] = cc*xi - ss*xr;
	} /* end counter rotate */
      } /* end channel loop */
    } /* end IF loop */


    /* Add diagnostics IF 1, ch 1*/
    if (prtLv>=2) {
      ix = (olong) (SumAzCell[i] + nx/2 + 1.5);
      iy = (olong) (SumElCell[i] + ny/2 + 1.5);
      if ((SumAzCell[i]>1000.) || (SumElCell[i]>1000.)) continue;
      ichan = 0; iIF = 0; 
      ichan = 1; iIF =11; /* DEBUG IF 12, ch 2 */
      off = i*selem + iIF*nchan + ichan;
      if (i==0)
	Obit_log_error(err, OBIT_InfoErr, "Data for IF %d, channel %d labeled as circular",
		       iIF+1, ichan+1);
      if (SumRLr) { /* Poln data? */
	Obit_log_error(err, OBIT_InfoErr, 
		       "%3.3d Cell %3d %3d Az%8.1f cell, El %6.1f cell, RR %6.3f %6.3f RL %6.3f %6.3f LR %6.3f %6.3f LL %6.3f %6.3f Jy",
		       i, ix,iy, 
		       /*SumAzCell[i]*xCells*206265., SumElCell[i]*yCells*206265., offset in asec */
		       SumAzCell[i], SumElCell[i],   /* offset in cells */
		       SumRRr[off],SumRRi[off], SumRLr[off],SumRLi[off], SumLRr[off],SumLRi[off], SumLLr[off],SumLLi[off]);
		       /* SumRRr[off],SumRRWt[off], SumRLr[off],SumRLWt[off], SumLRr[off],SumLRWt[off], SumLLr[off],SumLLWt[off]); DEBUG */
      } else {
	Obit_log_error(err, OBIT_InfoErr, 
		       "%3.3d Cell %3d %3d Az %8.1f cell, El %6.1f cell, RR %6.3f %6.3f Jy, LL %6.3f %6.3f Jy",
		       i, ix,iy, 
		       /*SumAzCell[i]*xCells*206265., SumElCell[i]*yCells*206265., offset in asec */
		       SumAzCell[i], SumElCell[i],   /* offset in cells */
		       SumRRr[off],SumRRi[off], SumLLr[off],SumLLi[off]);
      }
    }
    ObitErrLog(err);
  } /* End loop normalizing list */

  /* Cleanup */
  if (OffAz) g_free(OffAz);
  if (OffEl) g_free(OffEl);

} /* end accumData  */

/*----------------------------------------------------------------------- */
/*  Interpolate quasi regular lists of points onto a regular grid         */
/*  Accumulate real part of correlations for all Stokes/freq/IF           */
/*  Sums are in a list with each set of entries corresponding to a given  */
/*  pointing.                                                             */
/*  May use threading.                                                    */
/*  Adapted from AIPS/MAPBM.FOR                                           */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      nchan     Number of channels in output                            */
/*      nIF       Number of IFs in output                                 */
/*      npoln     Number of polarizations in output                       */
/*      selem     Size (floats) of list element                           */
/*      nelem     Number of list elements                                 */
/*      SumRRr     Real Stokes RR (or XX) accumulation list               */
/*      SumRRi     Imag Stokes RR (or XX) accumulation list               */
/*      SumRR2     Stokes RR**2 (or XX**2) accumulation list              */
/*      SumRRWt     RR (or XX) Weight accumulation list                   */
/*      SumRLr     Real Stokes RL (or XY) accumulation list               */
/*      SumRLi     Imag Stokes RL (or XY) accumulation list               */
/*      SumRL2     Stokes RL**2 (or X**2Y) accumulation list              */
/*      SumRLWt    RL (or XY) Weight accumulation list                    */
/*      SumLRr     Real Stokes LR (or YX) accumulation list               */
/*      SumLRi     Imag Stokes LR (or YX) accumulation list               */
/*      SumLR2     Stokes LR**2 (or YX**2) accumulation list              */
/*      SumLRWt    LR (or YX) Weight accumulation list                    */
/*      SumLLr     Real Stokes LL (or YY) accumulation list               */
/*      SumLLi     Imag Stokes  LL (or YY) accumulation list              */
/*      SumLL2     Stokes  LL**2 (or YY**2) accumulation list             */
/*      SumLLWt    LL (or YY) Weight accumulation list                    */
/*      SumAzCell Azimuth offset accumulation list                        */
/*      SumElCell Elevation offset  accumulation list                     */
/*   Output:                                                              */
/*      grids     Array of ObitFArrays                                    */
/*                fastest to slowest, channel, IF, Stokes                 */
/*                Stokes in order pp,qq,pq,qp                             */
/*      err       Obit Error stack                                        */
/*----------------------------------------------------------------------- */
void  gridData (ObitInfoList* myInput, olong nchan, olong nIF, olong npoln,
		olong selem, olong nelem, 
		ofloat *SumRRr, ofloat *SumRRi, ofloat *SumRR2, ofloat *SumRRWt,
		ofloat *SumRLr, ofloat *SumRLi, ofloat *SumRL2, ofloat *SumRLWt,
		ofloat *SumLRr, ofloat *SumLRi, ofloat *SumLR2, ofloat *SumLRWt,
		ofloat *SumLLr, ofloat *SumLLi, ofloat *SumLL2, ofloat *SumLLWt,
		ofloat *SumAzCell, ofloat *SumElCell, ObitFArray **grids)
{
  ObitInfoType type;
  gint32   dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong nThreads, nx, ny, hwid;
  ofloat xcen, ycen;
  gboolean OK;
  ObitThread *thread;
  MBFuncArg **threadArgs;

  /* Any thing to do? */
  if (nelem<2) return;

  /* How big an image? */
  ObitInfoListGetTest(myInput, "nx", &type, dim, &nx);
  ObitInfoListGetTest(myInput, "ny", &type, dim, &ny);
  ObitInfoListGetTest(myInput, "hwid",   &type, dim, &hwid);
  xcen = (ofloat)(nx/2);
  ycen = (ofloat)(ny/2);

  /* Initialize Threading */
  thread = newObitThread ();
  nThreads = 
    MakeMBFuncArgs (thread, 
		    nchan, nIF, npoln,selem, nelem, 
		    nx, ny, hwid, xcen, ycen,
		    SumRRr, SumRRi, SumRR2, SumRRWt, SumRLr, SumRLi, SumRL2, SumRLWt,
		    SumLRr, SumLRi, SumLR2, SumLRWt, SumLLr, SumLLi, SumLL2, SumLLWt,
		    SumAzCell, SumElCell, grids,
		    &threadArgs);
  /* Do operation */
  OK = ObitThreadIterator (thread, nThreads, 
			   (ObitThreadFunc)ThreadMB2Grid,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) return;

  /* Free local objects */
  KillMBFuncArgs(nThreads, threadArgs);
  ObitThreadUnref(thread);
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
/*      inc       stride in wlist                                         */
/*      xlist     List of coordinates on  first axis                      */
/*      ylist     List coordinate on second axis                          */
/*      wlist     List element weight                                     */
/*   Output:                                                              */
/*      coef      Array of interpolation coefficients for xlist,ylist     */
/*----------------------------------------------------------------------- */
void lagrange(ofloat x, ofloat y, olong n, olong hwid, olong inc,
	      ofloat *xlist, ofloat *ylist, ofloat *wlist, ofloat *coef)
{
  ofloat xhwid = (ofloat)hwid, sum;
  odouble prodx, prodxd, prody, prodyd;
  olong  i, j, countx, county;

  /* Closest cell?  */
  if (hwid<1) {
    for (j=0; j<n; j++) {
      coef[j] = 0.0;
      /*if ((fabs(x-xlist[j])<=0.5) && (fabs(y-ylist[j])<=0.5)) coef[j] = wlist[j*inc];*/
      if ((fabs(x-xlist[j])<=0.5) && (fabs(y-ylist[j])<=0.5) && 
	  (wlist[j*inc]>0.0)) coef[j] = 1.0;
    }
    return;
  }

  /* Lagrangian doesn't seem to work as intended */
  /* Loop over list */
  sum = 0.0;
  for (j=0; j<n; j++) {
    if (wlist[j*inc]<=0.0) continue;  /* Valid data? */
    prodx = prodxd = prody = prodyd = 1.0;
    countx = county = 0;
    coef[j] = 0.0;

    /* Within hwid? and i!=j */
    if ((fabs(x-xlist[j])<=xhwid) && (fabs(y-ylist[j])<=xhwid)) {
      coef[j] = 1.0;  /* In case nothing else within hwid */
     
      /* Inner loop over list */
      for (i=0; i<n; i++) {
	if (i==j) continue;                    /* i!=j */
	if (fabs(x-xlist[i])>xhwid) continue;  /* X within halfwidth */
	if (fabs(y-ylist[i])>xhwid) continue;  /* Y within halfwidth */
	if (fabs(xlist[j]-xlist[i])>0.3) {
	  countx++;
	  prodx  *= (odouble)(x - xlist[i]);
	  prodxd *= (odouble)(xlist[j] - xlist[i]);
	}
	if (fabs(ylist[j]-ylist[i])>0.3) {
	  county++;
	  prody  *= (odouble)(y - ylist[i]);
	  prodyd *= (odouble)(ylist[j] - ylist[i]);
	}
      } /* end inner loop */
    } /* end j within half width */
    /* put it together */
    if ((countx>=1)  || (county>=1)) {
      if ((prodxd!=0.0) && (prodyd!=0.0))
	coef[j] = (ofloat)(prodx*prody / (prodxd*prodyd));
      else
	coef[j] = 1.0;
      sum += coef[j];
    }
  } /* end loop over list */

  /* DEBUG 
  if ((x=1.0) && (y==1.0) && (sum!=0.0))
    fprintf(stderr,"lagrange x %f, y %f, sum %f\n", x, y, sum);
  if ((x=1.0) && (y==0.0) && (sum!=0.0))
    fprintf(stderr,"lagrange x %f, y %f, sum %f\n", x, y, sum);*/
  
  /* Normalize if anything found */
  if (fabs(sum)<0.01) return;
  prodx = 1.0 / sum;
  for (j=0; j<n; j++) coef[j] *= prodx;

  } /* end lagrange */

/*----------------------------------------------------------------------- */
/*  Interpolation coefficients                                            */
/*  For interpolating in a quasi regular grid represented by lists        */
/*  Determine coefficients for elements in lists to interpolate to (x,y)  */
/*  weight ~ exp(-distance**2)                                            */
/*   Input:                                                               */
/*      x         Coordinate on first axis                                */
/*      y         Coordinate on second axis                               */
/*      n         Length of lists                                         */
/*      hwid      Halfwidth of interpolation kernal                       */
/*      inc       stride in wlist                                         */
/*      xlist     List of coordinates on  first axis                      */
/*      ylist     List coordinate on second axis                          */
/*      wlist     List element weight                                     */
/*   Output:                                                              */
/*      coef      Array of interpolation coefficients for xlist,ylist     */
/*----------------------------------------------------------------------- */
void interpolate(ofloat x, ofloat y, olong n, olong hwid, olong inc,
		 ofloat *xlist, ofloat *ylist, ofloat *wlist, ofloat *coef)
{
  ofloat xhwid = (ofloat)hwid, d2, sum, prod;
  olong  j;

  /* Closest cell?  */
  if (hwid<1) {
    sum = 0;
    for (j=0; j<n; j++) {
      coef[j] = 0.0;
      if ((fabs(x-xlist[j])<=0.5) && (fabs(y-ylist[j])<=0.5)) 
	{coef[j] = wlist[j*inc]; sum += coef[j];}
     }
  } else {

    /* weight by exp(-dist**2) - Loop over list */
    sum = 0.0;
    for (j=0; j<n; j++) {
      coef[j] = 0.0;
      /* Within hwid?  */
      if ((fabs(x-xlist[j])<=xhwid) && (fabs(y-ylist[j])<=xhwid)) {
	/* exp(-(d2)) outside 0.1 cells */
	d2 = (x-xlist[j])*(x-xlist[j]) + (y-ylist[j])*(y-ylist[j]);
	if (d2<=0.01) coef[j] = wlist[j*inc];
	else          coef[j] = wlist[j*inc]*exp(-(d2-0.01));
	sum += coef[j];
      } /* end if within hwid of (x,y) */
    } /* end loop over list */
  } /* end weighted interpolation */

  /* Normalize if anything found */
  if (fabs(sum)<0.01) return;
  prod = 1.0 / sum;
  for (j=0; j<n; j++) coef[j] *= prod;
  } /* end interpolate */

/**
 * Make arguments for a Threaded ThreadMBFunc?
 * \param thread     ObitThread object to be used
 * \param  nchan     Number of channels in output
 * \param  npoln     Number of polarizations in output
 * \param  nIF       Number of IFs in output
 * \param  npoln     Number of polarizations in output
 * \param  selem     Size (floats) of list element
 * \param  nelem     Number of list elements
 * \param  nx;       Number of columns
 * \param  ny;       Number of rows
 * \param  hwid;     Half width of interpolation
 * \param  xcen;     X center cell
 * \param  ycen;     Y center cell
 * \param  SumRRr     Real Stokes RR (or XX) accumulation list
 * \param  SumRRi     Imag Stokes RR (or XX) accumulation list
 * \param  SumRR2     Stokes RR**2 (or XX**2) accumulation list
 * \param  SumRRw     Stokes RR weight list
 * \param  SumRLr     Real Stokes RL (or XY) accumulation list
 * \param  SumRLi     Imag Stokes RL (or XY) accumulation list
 * \param  SumRL2     Stokes RL**2 (or XY**2) accumulation list
 * \param  SumRlw     Stokes RL weight list
 * \param  SumLRr     Real Stokes LR (or YX) accumulation list
 * \param  SumLRi     Imag Stokes LR (or YX) accumulation list
 * \param  SumLR2     Stokes LR**2 (or YX**2) accumulation list
 * \param  SumLRw     Stokes LR weight list
 * \param  SumLLr     Real Stokes LL (or YY) accumulation list
 * \param  SumLLi     Imag Stokes LL (or YY) accumulation list
 * \param  SumLL2     Stokes LL**2 (or YY**2) accumulation list
 * \param  SumLLw     Stokes LL weight list
 * \param  SumAzCell Azimuth offset accumulation list
 * \param  SumElCell Elevation offset  accumulation list
 * \param  grids     [out] Array of ObitFArrays
 *                   fastest to slowest, channel, IF, Stokes
 * \param ThreadArgs[out] Created array of MBFuncArg, 
 *                   delete with KillMBFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeMBFuncArgs (ObitThread *thread, 
			     olong nchan, olong nIF, olong npoln, olong selem, olong nelem, 
			     olong nx, olong ny, olong hwid, ofloat xcen, ofloat ycen,
			     ofloat *SumRRr, ofloat *SumRRi, ofloat *SumRR2, ofloat *SumRRw,
			     ofloat *SumRLr, ofloat *SumRLi, ofloat *SumRL2, ofloat *SumRLw,
			     ofloat *SumLRr, ofloat *SumLRi, ofloat *SumLR2, ofloat *SumLRw,
			     ofloat *SumLLr, ofloat *SumLLi, ofloat *SumLL2, ofloat *SumLLw,
			     ofloat *SumAzCell, ofloat *SumElCell, ObitFArray **grids,
			     MBFuncArg ***ThreadArgs)

{
  olong i, iy, nThreads, nrow;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* How many rows per thread? */
  nrow = ny/nThreads;
  iy   = 0;

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(MBFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(MBFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->thread    = ObitThreadRef(thread);
    (*ThreadArgs)[i]->loy       = iy;     /* Divvy up processing */
    (*ThreadArgs)[i]->hiy       = MIN (iy+nrow, ny);
    if (i==(nThreads-1)) (*ThreadArgs)[i]->hiy = ny;  /* be sure to do all*/
    iy += nrow;
    (*ThreadArgs)[i]->nchan     = nchan;
    (*ThreadArgs)[i]->npoln     = npoln;
    (*ThreadArgs)[i]->nIF       = nIF;
    (*ThreadArgs)[i]->npoln     = npoln;
    (*ThreadArgs)[i]->selem     = selem;
    (*ThreadArgs)[i]->nelem     = nelem;
    (*ThreadArgs)[i]->nx        = nx;
    (*ThreadArgs)[i]->ny        = ny;
    (*ThreadArgs)[i]->hwid      = hwid;
    (*ThreadArgs)[i]->xcen      = xcen;
    (*ThreadArgs)[i]->ycen      = ycen;
    (*ThreadArgs)[i]->SumRRr    = SumRRr;
    (*ThreadArgs)[i]->SumRRi    = SumRRi;
    (*ThreadArgs)[i]->SumRR2    = SumRR2;
    (*ThreadArgs)[i]->SumRRw    = SumRRw;
    (*ThreadArgs)[i]->SumRLr    = SumRLr;
    (*ThreadArgs)[i]->SumRLi    = SumRLi;
    (*ThreadArgs)[i]->SumRL2    = SumRL2;
    (*ThreadArgs)[i]->SumRLw    = SumRLw;
    (*ThreadArgs)[i]->SumLRr    = SumLRr;
    (*ThreadArgs)[i]->SumLRi    = SumLRi;
    (*ThreadArgs)[i]->SumLR2    = SumLR2;
    (*ThreadArgs)[i]->SumLRw    = SumLRw;
    (*ThreadArgs)[i]->SumLLr    = SumLLr;
    (*ThreadArgs)[i]->SumLLi    = SumLLi;
    (*ThreadArgs)[i]->SumLL2    = SumLL2;
    (*ThreadArgs)[i]->SumLLw    = SumLLw;
    (*ThreadArgs)[i]->SumAzCell = SumAzCell;
    (*ThreadArgs)[i]->SumElCell = SumElCell;
    (*ThreadArgs)[i]->grids     = grids;
    (*ThreadArgs)[i]->coef[0]   = g_malloc0(nelem*sizeof(ofloat));
    (*ThreadArgs)[i]->coef[1]   = g_malloc0(nelem*sizeof(ofloat));
    (*ThreadArgs)[i]->coef[2]   = g_malloc0(nelem*sizeof(ofloat));
    (*ThreadArgs)[i]->coef[3]   = g_malloc0(nelem*sizeof(ofloat));
    (*ThreadArgs)[i]->ithread   = i;
  }

  return nThreads;
} /*  end MakeMBFuncArgs */

/**
 * Delete arguments for ThreadMBFunc
 * \param nargs      number of elements in ThreadArgs.
 * \param ThreadArgs Array of MBFuncArg
 */
static void KillMBFuncArgs (olong nargs, MBFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  ObitThreadPoolFree (ThreadArgs[0]->thread);  /* Free thread pool */
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      if (ThreadArgs[i]->thread) ObitThreadUnref(ThreadArgs[i]->thread);
      if (ThreadArgs[i]->coef[0])  g_free(ThreadArgs[i]->coef[0]);
      if (ThreadArgs[i]->coef[1])  g_free(ThreadArgs[i]->coef[1]);
      if (ThreadArgs[i]->coef[2])  g_free(ThreadArgs[i]->coef[2]);
      if (ThreadArgs[i]->coef[3])  g_free(ThreadArgs[i]->coef[3]);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillMBFuncArgs */

/**
 * Thread accumulate a set of rows in the image
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to MBFuncArg argument with elements:
 * \li  loy       First row (0-rel) 
 * \li  hiy       Highest row (0-rel) 
 * \li  nchan     Number of channels in output
 * \li  npoln     Number of polarizations in output
 * \li  nIF       Number of IFs in output
 * \li  npoln     Number of polarizations in output
 * \li  selem     Size (floats) of list element
 * \li  nelem     Number of list elements
 * \li  nx;       Number of columns
 * \li  ny;       Number of rows
 * \li  hwid;     Half width of interpolation
 * \li  xcen;     X center cell
 * \li  ycen;     Y center cell
 * \li  SumRRr     Real Stokes RR (or XX) accumulation list
 * \li  SumRRi     Imag Stokes RR (or XX) accumulation list
 * \li  SumRR2     Stokes RR**2 (or XX**2) accumulation list
 * \li  SumRRw     Stokes RR weight list
 * \li  SumRLr     Real Stokes RL (or XY) accumulation list
 * \li  SumRLi     Imag Stokes RL (or XY) accumulation list
 * \li  SumRL2     Stokes RL**2 (or XY**2) accumulation list
 * \li  SumRLw     Stokes RL weight list
 * \li  SumLRr     Real Stokes LR (or YX) accumulation list
 * \li  SumLRi     Imag Stokes LR (or YX) accumulation list
 * \li  SumLR2     Stokes LR**2 (or YX**2) accumulation list
 * \li  SumLRw     Stokes LR weight list
 * \li  SumLLr     Real Stokes LL (or YY) accumulation list
 * \li  SumLLi     Imag Stokes  LL (or YY) accumulation list
 * \li  SumLL2     Stokes  LL**2 (or YY**2) accumulation list
 * \li  SumLLw     Stokes LL weight list
 * \li  SumAzCell Azimuth offset accumulation list
 * \li  SumElCell Elevation offset  accumulation list
 * \li  grids     [out] Array of ObitFArrays
 *                   fastest to slowest, (real,imag, wt, ^2), channel, 
 *                   IF, Stokes in order pp,qq,pq,qp
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadMB2Grid (gpointer arg)
{
  /* Get arguments from structure */
  MBFuncArg *largs = (MBFuncArg*)arg;
  olong   loy        = largs->loy;
  olong   hiy        = largs->hiy;
  olong   nIF        = largs->nIF;
  olong   nchan      = largs->nchan;
  olong   npoln      = largs->npoln;
  olong   selem      = largs->selem;
  olong   nelem      = largs->nelem;
  olong  nx          = largs->nx; 
  /*olong  ny          = largs->ny;*/
  olong  hwid        = largs->hwid;
  ofloat xcen        = largs->xcen;  
  ofloat ycen        = largs->ycen;  
  ofloat* SumRRr     = largs->SumRRr;
  ofloat* SumRRi     = largs->SumRRi;
  ofloat* SumRR2     = largs->SumRR2;
  ofloat* SumRRw     = largs->SumRRw;
  ofloat* SumRLr     = largs->SumRLr;
  ofloat* SumRLi     = largs->SumRLi;
  ofloat* SumRL2     = largs->SumRL2;
  ofloat* SumRLw     = largs->SumRLw;
  ofloat* SumLRr     = largs->SumLRr;
  ofloat* SumLRi     = largs->SumLRi;
  ofloat* SumLR2     = largs->SumLR2;
  ofloat* SumLRw     = largs->SumLRw;
  ofloat* SumLLr     = largs->SumLLr;
  ofloat* SumLLi     = largs->SumLLi;
  ofloat* SumLL2     = largs->SumLL2;
  ofloat* SumLLw     = largs->SumLLw;
  ofloat* SumAzCell  = largs->SumAzCell;
  ofloat* SumElCell  = largs->SumElCell;
  ofloat* coef1      = largs->coef[0];
  ofloat* coef2      = largs->coef[1];
  ofloat* coef3      = largs->coef[2];
  ofloat* coef4      = largs->coef[3];
  ObitFArray** grids  = largs->grids;
  
  /* local */
  olong i, ix, iy, iIF, ichan, indx, jndx, off;
  ofloat x, y, closest;
  ofloat fblank = ObitMagicF();
  odouble sumRRWt , sumRLWt, sumLRWt, sumLLWt;
  odouble valRRr,  valRRi, valRR2, valRLr, valRLi, valRL2;
  odouble valLRr,  valLRi, valLR2, valLLr, valLLi, valLL2;
  ObitFArray *arrayr, *arrayi, *arrayw, *array2;

  if (hiy<loy) goto finish;
  
  /* Loop over y (El) */
  for (iy=loy; iy<hiy; iy++) {
    y = iy - ycen;
    
    /* Loop over x (Az) */
    for (ix=0; ix<nx; ix++) {
      x = ix - xcen;
      
      /* Loop over IFs */
      for (iIF=0; iIF<nIF; iIF++) {
	/* Loop over chans */
	for (ichan=0; ichan<nchan; ichan++) {
	  valRRr = valRR2 = valRLr = valRL2 = valLRr = valLR2 = valLLr = valLL2 = 0.0;
	  valRRi = valRLi = valLRi = valLLi = 0.0;
	  sumRRWt = sumRLWt = sumLRWt = sumLLWt = 0.0;
	  off = ichan + iIF*nchan;
	  /* Get interpolation coefficients */
	  if (SumRRw) lagrange (x, y, nelem, hwid, selem, SumAzCell, SumElCell, &SumRRw[off], coef1);
	  if (SumRLw) lagrange (x, y, nelem, hwid, selem, SumAzCell, SumElCell, &SumRLw[off], coef2);
	  if (SumLRw) lagrange (x, y, nelem, hwid, selem, SumAzCell, SumElCell, &SumLRw[off], coef3);
	  if (SumLLw) lagrange (x, y, nelem, hwid, selem, SumAzCell, SumElCell, &SumLLw[off], coef4);
	  /*if (SumRRw) interpolate (x, y, nelem, hwid, selem, SumAzCell, SumElCell, &SumRRw[off], coef1);*/
	  /*if (SumRLw) interpolate (x, y, nelem, hwid, selem, SumAzCell, SumElCell, &SumRLw[off], coef2);*/
	  /*if (SumLRw) interpolate (x, y, nelem, hwid, selem, SumAzCell, SumElCell, &SumLRw[off], coef3);*/
	  /*if (SumLLw) interpolate (x, y, nelem, hwid, selem, SumAzCell, SumElCell, &SumLLw[off], coef4);*/
	  closest = 1000.0; /* Closest element */
	  /* Loop over lists summing */
	  for (i=0; i<nelem; i++) {
	    closest = MIN (closest, MAX (fabs(SumAzCell[i]-x), fabs(SumElCell[i]-y)));
	    if ((SumRRr[i*selem+off]!=fblank) && (coef1[i]!=0.0)) {
	      valRRr  += coef1[i]*SumRRr[i*selem+off];
	      valRRi  += coef1[i]*SumRRi[i*selem+off];
	      valRR2  += coef1[i]*SumRR2[i*selem+off];
	      sumRRWt += coef1[i];
	    }
	    if (SumRLr && (SumRLr[i*selem+off]!=fblank) && (coef2[i]!=0.0)) {
	      valRLr  += coef2[i]*SumRLr[i*selem+off];
	      valRLi  += coef2[i]*SumRLi[i*selem+off];
	      valRL2  += coef2[i]*SumRL2[i*selem+off];
	      sumRLWt += coef2[i];
	    }
	    if (SumLRr && (SumLRr[i*selem+off]!=fblank) && (coef3[i]!=0.0)) {
	      valLRr  += coef3[i]*SumLRr[i*selem+off];
	      valLRi  += coef3[i]*SumLRi[i*selem+off];
	      valLR2  += coef3[i]*SumLR2[i*selem+off];
	      sumLRWt += coef3[i];
	    }
	    if ((SumLLr && SumLLr[i*selem+off]!=fblank) && (coef4[i]!=0.0)) {
	      valLLr  += coef4[i]*SumLLr[i*selem+off];
	      valLLi  += coef4[i]*SumLLi[i*selem+off];
	      valLL2  += coef4[i]*SumLL2[i*selem+off];
	      sumLLWt += coef4[i];
	    }
	  } /* end loop over lists */
	  /* Better be something within 0.5 cells */
	  if (closest>0.5) {
	    sumRRWt = sumRLWt = sumLRWt = sumLLWt = 0.0;
	  }
	  
	  /* Accumulate into FArrays */
	  jndx  = iy*nx + ix;  /* Cell in grid */
	  /* qq */
	  indx  = 4*(ichan + iIF*nchan);
	  arrayr = grids[indx];
	  arrayi = grids[indx+1];
	  arrayw = grids[indx+2];
	  array2 = grids[indx+3];
	  if (fabs(sumRRWt)>0.0) {
	    arrayr->array[jndx] += valRRr;
	    arrayi->array[jndx] += valRRi;
	    arrayw->array[jndx] += sumRRWt;
	    array2->array[jndx] += valRR2;
	  }

	  /* pq */
	  if (npoln>=3) {
	    indx  = 4*(ichan + iIF*nchan + 2*nchan*nIF);
	    arrayr = grids[indx];
	    arrayi = grids[indx+1];
	    arrayw = grids[indx+2];
	    array2 = grids[indx+3];
	    if (fabs(sumRRWt)>0.0) {
	      arrayr->array[jndx] += valRLr;
	      arrayi->array[jndx] += valRLi;
	      arrayw->array[jndx] += sumRLWt;
	      array2->array[jndx] += valRL2;
	    }
	  }

	  /* qp */
	  if (npoln>=4) {
	    indx  = 4*(ichan + iIF*nchan + 3*nchan*nIF);
	    arrayr = grids[indx];
	    arrayi = grids[indx+1];
	    arrayw = grids[indx+2];
	    array2 = grids[indx+3];
	    if (fabs(sumRRWt)>0.0) {
	      arrayr->array[jndx] += valLRr;
	      arrayi->array[jndx] += valLRi;
	      arrayw->array[jndx] += sumLRWt;
	      array2->array[jndx] += valLR2;
	    }
	  }

	  /* qq */
	  if (npoln>=2) {
	    indx  = 4*(ichan + iIF*nchan + nchan*nIF);
	    arrayr = grids[indx];
	    arrayi = grids[indx+1];
	    arrayw = grids[indx+2];
	    array2 = grids[indx+3];
	    if (fabs(sumRRWt)>0.0) {
	      arrayr->array[jndx] += valLLr;
	      arrayi->array[jndx] += valLLi;
	      arrayw->array[jndx] += sumLLWt;
	      array2->array[jndx] += valLL2;
	    }
	  }
	} /* end channel Loop */
      } /* end IF Loop */
    } /* end x loop */
  } /* end y loop */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadMB2Grid */


/**
 * Normalize image and convert to type
 * \param doRMS     If TRUE, image is RMS
 * \param doPhase   If TRUE, image is Phase, else Amplitude
 * \param arrayr    Weighted Real part, returns if NULL
 * \param arrayi    Weighted Imaginary part
 * \param arrayw    Sum of weights Weight
 * \param array2    Weighted Amp^2
 * \param out       output array
 */
static void normArray(gboolean doRMS, gboolean doPhase, 
		      ObitFArray *arrayr, ObitFArray *arrayi, ObitFArray *arrayw, 
		      ObitFArray *array2, ObitFArray *out) 
{
  ofloat re, im, amp2, amp, ph, fblank = ObitMagicF();
  olong i;

  /* Anything to do? */
  if (arrayr==NULL) return;

  /* Loop over elements */
  for (i=0; i<arrayr->arraySize; i++) {
    if (fabs(arrayw->array[i])==0.0) { /* no data, blank */
      out->array[i] = fblank;
    } else { /* normalize */
      re   = arrayr->array[i]/arrayw->array[i];
      im   = arrayi->array[i]/arrayw->array[i];
      amp2 = array2->array[i]/arrayw->array[i];
      /* Amplitude and phase within +/- 90 deg */
      amp = sqrt (re*re + im*im);
      ph  = RAD2DG*atan2(im, re);
      if (ph>90.0) {
	amp = -amp;
	ph -= 180.0;
      } else if (ph<-90.0) {
	amp = -amp;
	ph += 180.0;
      }
      /* Save type requested */
      if (doRMS) {  /* RMS */
	out->array[i] = amp2;
      } else if (doPhase) {  /* phase */
	out->array[i] = ph;
      } else {      /* Stokes ampl */
	out->array[i] = amp;
      }
      /* DEBUG/
	 out->array[i] = arrayw->array[i];  */
    } /* end OK */
  } /* end loop over array */
} /* end normArray */
