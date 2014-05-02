/* $Id$  */
/* Makes images of catalog sources for initial calibration            */
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

#include "ObitImageUtil.h"
#include "ObitUVUtil.h"
#include "ObitSystem.h"
#include "ObitInfoList.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitFITS.h"
#include "ObitTableVL.h"
#include "ObitTableCC.h"
#include "ObitPBUtil.h"

/* Source list (GSList) stuff */
/* Catalog list management */
/** Source list element structure */
typedef struct { 
  /**  catalog celestial position */
  odouble ra, dec;
  /** shift from reference position */
  ofloat shift[3];
  /** deconvolved gaussian parameters (maj, min, PA) (FWHM, deg) */
  ofloat gparm[3];
  /** Estimated catalog flux density */
  ofloat flux;
}  SouListElem; 

/** SouList structure. */
typedef struct {
  /** Number of entries */
  olong number;
  /** glib singly linked list */
  GSList* list;
} SouList;

  
/* internal prototypes */
/* Get inputs */
ObitInfoList* FacesIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void FacesOut (ObitInfoList* outList, ObitErr *err);
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
/* Loop over sources */
void doSources (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Get output image */
ObitImage* getOutputImage (ObitInfoList *myInput, ObitUV* inData,
			   gchar *Source, ObitErr *err);
/* Write history */
void FacesHistory (gchar *Source, ObitInfoList* myInput, 
		    ObitUV* inData, ObitImage* outImage, ObitErr* err);
/* Get source list from catalog */
SouList* ParseCat (ObitInfoList* myInput, ObitImage* outImage,  
		   ObitErr *err);
/* Add source list to CC Table */
void SouList2TableCC (SouList* slist, ObitImage* outImage, olong CCVer, 
		      ObitErr *err);
/* Add source list to Image */
void SouList2Image (SouList* slist, ObitImage* outImage, ObitErr *err);


/* Program globals */
gchar *pgmName = "Faces";       /* Program name */
gchar *infile  = "Faces.in" ;   /* File with program inputs */
gchar *outfile = "Faces.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */

/** Private: SouListElem Constructor  */ 
 SouListElem* 
newSouListElem (odouble ra, odouble dec, ofloat shift[3],
		ofloat gparm[3], ofloat flux); 

/**  Private: Update contents of an SouListElem */
void
SouListElemUpdate (SouListElem *elem, odouble ra, odouble dec, 
		   ofloat shift[3], ofloat gparm[3], ofloat flux);

/**  Private: Print contents of an SouListElem */
void
SouListElemPrint (SouListElem *elem, FILE *file);

/** Private: SouListElem Destructor  */ 
void freeSouListElem(SouListElem *me); 
  
/** Private: SouList structure.constructor. */
SouList* newSouList (void);

/** Private: Add an SouListElem to the SouList  */
void SouListAdd (SouList *in, SouListElem *elem);

/** Private: Remove a SouListElem from the list. */
void SouListRemove (SouList *in, SouListElem *elem);

/** Private: Remove all items from the list. */
void SouListClear (SouList *in);

/** Private: Remove all items from the list. */
void SouListPrint (SouList *in, FILE *file);

/** Private: destructor. */
void freeSouList (SouList *in);

/* Private: Deconvolve Gaussian  */
olong deconv (ofloat fmaj, ofloat fmin, ofloat fpa, 
	      ofloat cmaj, ofloat cmin, ofloat cpa, 
	      ofloat peak, ofloat irms,
	      ofloat *rmaj, ofloat *rmin, ofloat *rpa);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*  Makes images of catalog sources for initial calibration               */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitErr      *err= NULL;
 
   /* Startup - parse command line */
  err = newObitErr();
  myInput = FacesIn (argc, argv, err);
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

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* FacesIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "FacesIn";

  /* error checks */
  if (err->error) return list;

  /* Make default inputs/outputs InfoList */
  list = defaultInputs(err);
  myOutput = defaultOutputs(err);

   /* Initialize output */
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
      
     } else if (strcmp(arg, "-outDType") == 0) { /* Image type AIPS or FITS */
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

  return list;
} /* end FacesIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: Faces -input file -output ofile [args]\n");
    fprintf(stderr, "Create model of field from catalog\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def Faces.in\n");
    fprintf(stderr, "  -output output result file, def Faces.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -AIPSuser AIPS user number, def 2 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input \n");
    fprintf(stderr, "  -outDType AIPS or FITS type for output\n");
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
/*     outClass  Str [6]    output AIPS image class  [sClean]             */
/*     outSeq    Int        output AIPS image sequence no  [new]          */
/*     Sources   Str (16,1) Sources selected, blank = all                 */
/*     FOV       Flt (1)    Field of view in deg , NO DEFAULT (0.0)       */
/*     Catalog   Str (48)   Outlier catalog name, def 'NVSSVZ.FIT'        */
/*     maxDist   Flt (1)    Maximum distance to add sources (deg), def=0  */
/*     minFlux   Flt (1)    Min. estimated source flux den. (Jy), def=0   */
/*     SI        Flt (1)    Spectral index to est. flux den., def=-0.7    */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat ftemp;
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
  strTemp = "Faces.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "FacesName";
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
  strTemp = "FacesOut.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS Image file name */
  strTemp = "FacesOut";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS Image file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outClass", OBIT_string, dim, strTemp, err);
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

  /* Field of view in deg def = 0.0 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "FOV", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Outlier catalog name, def = NVSSVL.FIT" */
  strTemp = "NVSSVL.FIT";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Catalog", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Maximum distance to add sources (deg),, def=0.0 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "maxDist", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Minimum estimated source flux density (Jy), def = 0.0 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "minFlux", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Spectral index to estimate flux density, def = -0.7 */
  dim[0] = 1;dim[1] = 1;
  ftemp = -0.7; 
  ObitInfoListPut (out, "OutlierSI", OBIT_float, dim, &ftemp, err);
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
  gchar        outClass[12], *strTemp;
  /*gchar *routine = "digestInputs";*/
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* Default outclass = "Faces" */
  if  (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
    strncpy (outClass, strTemp, 7);
  } else { /* Didn't find */
    strncpy (outClass, "Faces", 7);
  }
  /* Replace blank */
  if (!strncmp("      ", outClass, 6))  strncpy (outClass, "Faces", 7);
  outClass[6] = 0;
  /* input AIPS sequence */
  dim[0] = 6;
  ObitInfoListAlwaysPut(myInput, "outClass", OBIT_string, dim, outClass);
  

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
  olong        nvis=1000;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", 
   NULL};
  gchar *routine = "getInputData";

  /* error checks */
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Build basic input UV data Object */
  inData = ObitUVFromFileInfo ("in", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Set buffer size */
  nvis = 1000;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO",  OBIT_long, dim,  &nvis);

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Ensure inData fully instantiated and OK and selector set */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Loop over selected sources, these are all sources in the source table */
/*  with selection by Sources.                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to image                                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doSources  (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  gchar        Source[17];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitImage *outImage=NULL;
  ObitSourceList* doList;
  olong maxlen, isource, CCVer=1;
  SouList *slist=NULL;
  gchar *dataParms[] = {  /* Source selection*/
    "Sources", 
    NULL
  };
  gchar *routine = "doSources";


  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* CC Table version */
  ObitInfoListGetTest(myInput, "CCVer", &type, dim, &CCVer);

  /* Make sure selector set on inData */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  
  /* Get source list to do */
  doList = ObitUVUtilWhichSources (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Loop over list of sources */
  for (isource = 0; isource<doList->number; isource++) {
    if (!doList->SUlist[isource]) continue; /* removed? */
    maxlen = MIN (16, strlen(doList->SUlist[isource]->SourceName));
    strncpy (Source, doList->SUlist[isource]->SourceName, maxlen);
    Source[maxlen] = 0;
    ObitTrimTrail (Source);   /* Get rid of any trailing blanks */

    Obit_log_error(err, OBIT_InfoErr, " ******  Source %s ******", Source);

    /* Generate image and CC Table */
    outImage = getOutputImage (myInput, inData, Source, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);

    /* Get source List for outImage */
    slist =  ParseCat (myInput, outImage, err);

    /* Attach to output */
    SouList2TableCC (slist, outImage, CCVer, err);
    SouList2Image (slist, outImage, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);

    /* History */
    FacesHistory (Source, myInput, inData, outImage, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);

    /* cleanup */
    outImage = ObitImageUnref(outImage);
    freeSouList(slist);
  } /* end source loop */

  doList = ObitSourceListUnref(doList);

}  /* end doSources */


/*----------------------------------------------------------------------- */
/*  Write History for Faces                                               */
/*   Input:                                                               */
/*      Source    Name of source being imaged                             */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outImage  ObitImage to write history to                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void FacesHistory (gchar *Source, ObitInfoList* myInput, 
		   ObitUV* inData, ObitImage* outImage, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "outDType", "outFile",  "outDisk", "outName", "outClass", "outSeq",
    "Catalog", "CatDisk", "maxDist",  "minFlux", "SI", "CCVer", 
    "FOV", "xCells", "yCells", "antSize", "usegauss",
    NULL};
  gchar *routine = "FacesHistory";

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

  /* Write source  */
  g_snprintf (hicard, 80, "%s Source = '%s', ", pgmName, Source);
  ObitHistoryWriteRec (outHistory, -1, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end FacesHistory  */

/*----------------------------------------------------------------------- */
/*  Create output image                                                   */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy from                                     */
/*      Source    Name of source being imaged                             */
/*      outImage  ObitImage to write history to                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return:                                                              */
/*     output image                                                       */
/*----------------------------------------------------------------------- */
ObitImage* getOutputImage (ObitInfoList *myInput, ObitUV* inData,
			   gchar *Source, ObitErr *err)
{
  ObitImage *outImage=NULL;
  ObitTableCC *outCC=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat FOV, xCells, yCells;
  olong CCVer, highVer;
  gboolean exist=FALSE;
  gchar oldName[32], tname[129];
  gchar *routine = "getOutputImage";

  /* error checks */
  if (err->error) return outImage;

  /* Get old output name */
  ObitInfoListGet(myInput, "outName", &type, dim, oldName, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outImage);

  /* Add source name */
  g_snprintf (tname, 100, "%s%s", Source, oldName);
  dim[0] = 12;
  ObitInfoListAlwaysPut(myInput, "outName",  OBIT_string, dim, tname);

  /* Shouldn't exist */
  dim[0] = 1;
  ObitInfoListAlwaysPut(myInput, "outExist",  OBIT_bool, dim, &exist);
 
  /* Create image */
  outImage = ObitImageFromFileInfo("out", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outImage);

  /* Reset old output name */
  dim[0] = 12;
  ObitInfoListAlwaysPut(myInput, "outName",  OBIT_string, dim, oldName);

  /* Set selected source */
  dim[0] = 16; dim[1] = 1;
  ObitInfoListAlwaysPut(inData->info, "Sources",  OBIT_string, dim, Source);
  /* Open UV for selected source */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outImage);

  /* Get info from UV Data */
  ObitImageUtilUV2ImageDesc(inData->myDesc, outImage->myDesc, TRUE, 100000);

  /* Set size, cell spacing */
  ObitInfoListGet(myInput, "FOV",     &type, dim, &FOV,     err);
  ObitInfoListGet(myInput, "xCells",  &type, dim, &xCells,  err);
  ObitInfoListGet(myInput, "yCells",  &type, dim, &yCells,  err);
  ObitInfoListGet(myInput, "CCVer",   &type, dim, &CCVer,   err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outImage);
  outImage->myDesc->inaxes[outImage->myDesc->jlocr] = 
    (olong)(0.5 + 2. * FOV / (fabs(xCells)/3600.));
  outImage->myDesc->inaxes[outImage->myDesc->jlocd] = 
    (olong)(0.5 + 2. * FOV / (fabs(yCells)/3600.));
  outImage->myDesc->crpix[outImage->myDesc->jlocr] = 
    outImage->myDesc->inaxes[outImage->myDesc->jlocr]/2;
  outImage->myDesc->crpix[outImage->myDesc->jlocd] = 
    outImage->myDesc->inaxes[outImage->myDesc->jlocd]/2;
  outImage->myDesc->cdelt[outImage->myDesc->jlocr]  = -fabs(xCells)/3600.;
  outImage->myDesc->cdelt[outImage->myDesc->jlocd]  =  fabs(yCells)/3600.;
  ObitImageFullInstantiate (outImage, FALSE, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outImage);

  /* Zero fill */
  ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outImage);
  ObitFArrayFill(outImage->image, 0.0);
  ObitImageWrite (outImage, NULL, err);
  ObitImageClose (outImage, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outImage);

   /* Delete old CC Table */
  highVer = ObitTableListGetHigh (outImage->tableList, "AIPS CC");
  if ((CCVer>0) && (CCVer<=highVer))
    ObitImageZapTable (outImage, "AIPS CC", CCVer, err);
  if (err->error) Obit_traceback_val (err, routine, outImage->name, outImage);

 /* Make CC Table */
  outCC = newObitTableCCValue ("CC", (ObitData*)outImage,
			       &CCVer, OBIT_IO_WriteOnly, 4, err);
  /* Open/close to generate */
  ObitTableCCOpen (outCC, OBIT_IO_WriteOnly, err);
  ObitTableCCClose (outCC, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outImage);
  outCC = ObitTableCCUnref(outCC);

  return outImage;
} /* end getOutputImage */

/**
 * Searches Catalog for sources in the desired ra, dec, and flux 
 * range taking into account the estimated single-dish beam 
 * \param myInput      Task parameter list
 * \param outImage     Image centered on position
 * \param err          Error stack, returns if not empty.
 * \return SouList with entries
 */
SouList* ParseCat (ObitInfoList* myInput, ObitImage* outImage,  
		   ObitErr *err) 
{
  SouList *outList = NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong count;
  odouble ra0, dec0, Freq, ra, dec, ra2000, dc2000, dra;
  odouble xx, yy, zz, dist, refreq;
  ofloat radius, minflux, diam, alpha, pbf, sumFlux=0.0;
  ofloat flux, scale, xsh[3], gparm[3], fparm[3], beam[3];
  gfloat rotate=0.0, dxyzc[3];
  gboolean wanted, doJinc, useGauss=FALSE;
  olong blc[IM_MAXDIM] = {1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0};
  olong ver, nrows, irow, CatDisk;
  gboolean doJ2B=FALSE;
  ObitIOCode retCode;
  ObitImage *VLImage=NULL;
  ObitTableVL *VLTable=NULL;
  ObitTableVLRow *VLRow=NULL;
  gchar Catalog[128];
  gchar *routine = "ParseCat";
  
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outList;

  /* get control parameters */
  ObitInfoListGet(myInput, "Catalog",  &type, dim, Catalog,  err);
  ObitInfoListGet(myInput, "CatDisk",  &type, dim, &CatDisk, err);
  ObitInfoListGet(myInput, "maxDist",  &type, dim, &radius,  err);
  ObitInfoListGet(myInput, "minFlux",  &type, dim, &minflux, err);
  ObitInfoListGet(myInput, "SI",       &type, dim, &alpha,   err);
  ObitInfoListGet(myInput, "antSize",  &type, dim, &diam,    err);
  ObitInfoListGet(myInput, "useGauss", &type, dim, &useGauss, err);
  if (err->error) Obit_traceback_val (err, routine, outImage->name, outList);

  /* Really needed? */
  if (radius<=1.0e-10) return outList;

  /* Info from image */
  ra0   = outImage->myDesc->crval[outImage->myDesc->jlocr];
  dec0  = outImage->myDesc->crval[outImage->myDesc->jlocd];
  Freq  = outImage->myDesc->crval[outImage->myDesc->jlocf];
  doJ2B = outImage->myDesc->equinox==1950.0;
  /* get j2000 position to lookup in Catalog in radians */
  ra2000 = ra0;
  dc2000 = dec0;
  if (doJ2B) ObitSkyGeomBtoJ (&ra2000, &dc2000);
  ra2000 *= DG2RAD;
  dc2000 *= DG2RAD;

  /* set defaults. */
  if (radius <= 0.0) radius = 15.0;
  if (diam   <= 0.0) diam   = 25.0;
  if (alpha  == 0.0) alpha  = -0.75;

  /* which beam model to use */
  doJinc = (Freq >= 1.0e9);

  /* Open Catalog (VL table on an image) */
  VLImage = newObitImage("Catalog image");
  ObitImageSetFITS(VLImage, OBIT_IO_byPlane, CatDisk, Catalog, blc, trc, err);

  /* Open to fully instantiate */
  ObitImageOpen(VLImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, VLImage->name, outList);

  /* Now get VL table */
  ver = 1;
  VLTable =  newObitTableVLValue("Catalog table", (ObitData*)VLImage, &ver, 
				 OBIT_IO_ReadOnly, err);
  ObitTableVLOpen(VLTable, OBIT_IO_ReadOnly, err);
  VLRow =  newObitTableVLRow (VLTable);  /* Table row */
  if (err->error) Obit_traceback_val (err, routine, VLTable->name, outList);

  /* Get table info */
  /*refreq  = VLTable->refFreq; shit */
  refreq  = 1.4e9;  
  nrows   = VLTable->myDesc->nrow;
  beam[0] = VLTable->BeamMajor;
  beam[1] = VLTable->BeamMinor;
  beam[2] = VLTable->BeamPA;

  /* frequency scaling */
  scale = pow ((Freq / refreq), alpha);

  /* Create output list */
  outList = newSouList ();

  /* loop through table */
  count = 0;
  for (irow= 1; irow<=nrows; irow++) { /* loop 500 */
    /* read */
    retCode = ObitTableVLReadRow (VLTable, irow, VLRow, err);
    if (err->error) Obit_traceback_val (err, routine, VLTable->name, outList);
   
    /* spectral scaling of flux density */
    flux = VLRow->PeakInt * scale;

    /* position, etc */
    ra   = VLRow->Ra2000;
    dec  = VLRow->Dec2000;

    /* select (crude) */
    xx = DG2RAD * ra;
    yy = DG2RAD * dec;
    dra = fabs (ra0-ra);
    if (dra>180.0) dra -= 360.0;
    if ((fabs(dec0-dec) <= 1.2*radius) && 
	(fabs(dra)*cos(yy) <= 1.2*radius) &&  
	(flux >= minflux)) {
      /* separation from pointing center */
      zz = sin (yy) * sin (dc2000) + cos (yy) * cos (dc2000) * cos (xx-ra2000);
      zz = MIN (zz, 1.000);
      dist = acos (zz) * RAD2DG;
      wanted =  (dist <= radius);

      /* primary beam correction to flux density */
      if (wanted) {
	if (doJinc) {
	  pbf = ObitPBUtilJinc (dist, Freq, diam, 0.05);
	} else {
	  pbf = ObitPBUtilPoly (dist, Freq, 0.05);
	} 
	flux *= MAX (0.05, pbf); /* Don't trust below 5% */
      }
      
      /* select (fine) */
      wanted = ((flux >= minflux)  && wanted);
      if (wanted) {
	if (doJ2B) {  /* precess to 1950 if necessary */
	  ObitSkyGeomJtoB (&ra, &dec);
	} 

	/* get shift needed */
	if (!strncmp(&outImage->myDesc->ctype[0][4], "-SIN", 4)) 
	  ObitSkyGeomShiftSIN (ra0, dec0, rotate, ra, dec, dxyzc);
	else if (!strncmp(&outImage->myDesc->ctype[0][4], "-NCP", 4)) 
	  ObitSkyGeomShiftNCP (ra0, dec0, rotate, ra, dec, dxyzc);
	else
	  ObitSkyGeomShiftSIN (ra0, dec0, rotate, ra, dec, dxyzc);

	xsh[0] =  dxyzc[0] / (DG2RAD * 2.0 * G_PI);
	xsh[1] =  dxyzc[1] / (DG2RAD * 2.0 * G_PI);
	xsh[2] =  dxyzc[2] / (DG2RAD * 2.0 * G_PI);

	/* add it to linked list */
	fparm[0] = VLRow->MajorAxis;
	fparm[1] = VLRow->MinorAxis;
	fparm[2] = VLRow->PosAngle;
	if (useGauss) {
	  /* Deconvolve */
	  deconv (fparm[0], fparm[1], fparm[2], beam[0], beam[1], beam[2], 
		  VLRow->PeakInt, VLRow->IRMS, 
		  &gparm[0], &gparm[1], &gparm[2]);
	} else {  /* Use as point */
	  gparm[0] = gparm[1] = gparm[2] = 0.0;
	}
	SouListAdd(outList, newSouListElem(ra, dec, xsh, gparm, flux));
	sumFlux += flux;   /* How much? */
	count ++;
      } /* end if wanted */ 
    } /* end crude selection */
  } /* end loop over table */

  /* tell how many */
  if (err->prtLv>=1) Obit_log_error(err, OBIT_InfoErr, 
				    "%s: Included %d sources with %f Jy", 
				    routine, count, sumFlux);
  /* Dump list */
  if (err->prtLv>2) SouListPrint (outList, stderr);

  /* Close up */
  ObitImageClose(VLImage, err);
  retCode = ObitTableVLClose(VLTable, err);
  if (err->error) Obit_traceback_val (err, routine, VLTable->name, outList);
  
  /* clean up */
  VLImage = ObitImageUnref(VLImage);
  VLTable = ObitTableUnref(VLTable);
  VLRow   = ObitTableRowUnref(VLRow);
  
  return outList;
} /* end of routine ParseCat */ 

/**
 * Add source list to AIPS CC table 
 * \param slist        List of sources to add
 * \param outImage     Image centered on position
 * \param err          Error stack, returns if not empty.
 */
void SouList2TableCC (SouList* slist, ObitImage* outImage, olong CCVer, 
		      ObitErr *err)
{
  GSList *tmp;
  SouListElem *elem;
  ObitTableCC *CCTable=NULL;
  ObitTableCCRow *CCRow=NULL;
  olong orow, noParms;
  gchar *routine = "SouList2TableCC";
  
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* Now get CC table */
  noParms = 4;
  CCTable =  newObitTableCCValue("Catalog table", (ObitData*)outImage, &CCVer, 
				 OBIT_IO_ReadWrite, noParms, err);
  ObitTableCCOpen(CCTable, OBIT_IO_ReadWrite, err);
  CCRow =  newObitTableCCRow (CCTable);  /* Table row */
  if (err->error) Obit_traceback_msg (err, routine, CCTable->name);

  /* loop through source list processing elements */
  tmp = slist->list;
  while (tmp!=NULL) {
    if (tmp->data) {
      /*  Append to table */
      elem            =  (SouListElem*)tmp->data;
      CCRow->Flux     = elem->flux;
      CCRow->DeltaX   = elem->shift[0];
      CCRow->DeltaY   = elem->shift[1];
      CCRow->DeltaZ   = elem->shift[2];
      CCRow->parms[0] = elem->gparm[0];
      CCRow->parms[1] = elem->gparm[1];
      CCRow->parms[2] = elem->gparm[2];
      CCRow->parms[3] = 1.0;
      /* Write row */
      orow = -1;
      ObitTableCCWriteRow (CCTable, orow, CCRow, err);
      if (err->error) Obit_traceback_msg (err, routine, CCTable->name);
    }
    tmp = g_slist_next(tmp);
  }
  /* Close/cleanup */
  ObitTableCCClose (CCTable, err);
  CCRow = ObitTableCCRowUnref(CCRow);
  CCTable = ObitTableCCUnref(CCTable);
  if (err->error) Obit_traceback_msg (err, routine, CCTable->name);

 } /* end SouList2TableCC */

/**
 * Add source list to Image 
 * \param slist        List of sources to add
 * \param outImage     Image centered on position
 * \param err          Error stack, returns if not empty.
 */
void SouList2Image (SouList* slist, ObitImage* outImage, ObitErr *err)
{
  GSList *tmp;
  SouListElem *elem;
  ObitFArray *pixels;
  ObitImageDesc *outDesc, *tmpDesc=NULL;
  ObitFArray *plist=NULL;
  olong pln[5], naxis[3], pos[2];
  ofloat *ppix, gauss[3], cellx, celly, bmaj, bmin, bpa, sr, cr, inPixl[2], outPixl[2];
  gchar *routine = "SouList2Image";
  
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* Now get image pixel array */
  ObitImageOpen(outImage, OBIT_IO_ReadWrite, err);
  pln[0] =  pln[1] = pln[2] = pln[3] = pln[4] = 1;
  ObitImageGetPlane (outImage, NULL, pln, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  pixels = outImage->image;

  /* Need fake ImageDescriptor for sources to add */
  outDesc = outImage->myDesc;
  tmpDesc = ObitImageDescCopy (outDesc, tmpDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  tmpDesc->inaxes[tmpDesc->jlocr] = 15;  /* Pretend it's small */
  tmpDesc->inaxes[tmpDesc->jlocd] = 15;
  inPixl[0] = tmpDesc->inaxes[tmpDesc->jlocr]/2;
  inPixl[1] = tmpDesc->inaxes[tmpDesc->jlocd]/2;
  tmpDesc->crpix[tmpDesc->jlocr] = inPixl[0];
  tmpDesc->crpix[tmpDesc->jlocd] = inPixl[1];

  /* Image grid info */
  cellx = outDesc->cdelt[outDesc->jlocr];
  celly = outDesc->cdelt[outDesc->jlocd];

  /* Create plist for positions, fluxes */
  naxis[0] = 3;
  plist = ObitFArrayCreate ("plist", 1, naxis);

  /* loop through source list processing elements */
  tmp = slist->list;
  while (tmp!=NULL) {
    if (tmp->data) {
      elem            =  (SouListElem*)tmp->data;
      /*  Add to image if needed */
      tmpDesc->crval[tmpDesc->jlocr] = elem->ra;
      tmpDesc->crval[tmpDesc->jlocd] = elem->dec;
      if (ObitImageDescCvtPixel(tmpDesc, outDesc, inPixl, outPixl, err)) {
	/* Cell on exact pixel, flux to plist */
	plist->array[0] = outPixl[0];
	plist->array[1] = outPixl[1];
	plist->array[2] = elem->flux;
	/* Gaussian in terms of sigmas and cells */
	bmaj = elem->gparm[0]/fabs(cellx);
	bmin = elem->gparm[1]/fabs(cellx);
	bpa  = elem->gparm[2];
	/* Hack for "semi" resolved */
	if ((bmaj>0.0) && (bmin<=0.0)) bmin = cellx;
	/* Point or Gaussian? Normalize Gaussian to unit area */
	if (bmaj>0.0) {
	    cr = cos ((bpa + outDesc->crota[outDesc->jlocd])*DG2RAD);
	    sr = sin ((bpa + outDesc->crota[outDesc->jlocd])*DG2RAD);
	    gauss[0] = ((cr*cr)/(bmin*bmin) + (sr*sr)/(bmaj*bmaj)) *
	      4.0*log(2.0);
	    gauss[1] =  ((sr*sr)/(bmin*bmin) + (cr*cr)/(bmaj*bmaj)) *
	      4.0*log(2.0);
	    gauss[2] = (1.0/(bmin*bmin) - 1.0/(bmaj*bmaj)) *
	      sr*cr*8.0*log(2.0);
	    plist->array[2] /= (2*G_PI*bmaj*bmin);
	}
	/* Add */
	if (bmaj>0.0) {
	  ObitFArrayConvGaus (pixels, plist, 1, gauss);
	} else { /* Point, just set pixel */
	  pos[0] = (olong)(outPixl[0]+0.5);
	  pos[1] = (olong)(outPixl[1]+0.5);
	  ppix   = ObitFArrayIndex (pixels, pos);
	  if (ppix!=NULL) *ppix  = elem->flux;
	}
	if (err->error) Obit_traceback_msg (err, routine, outImage->name);
      } /* end if overlap */
    }
    tmp = g_slist_next(tmp);
  }
  /* Write image */
  ObitImagePutPlane (outImage, NULL, pln, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Close/cleanup */
  ObitImageClose (outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
} /* end SouList2Image */


/*                 SouListElem functions             */
/**
 * SouListElem Constructor
 * \param ra      Catalog celestial position RA (deg)
 * \param dec     Catalog celestial position RA (deg)
 * \param shift   Shift from reference position
 * \param gparm   Gaussian parameters
 * \param flux    Estimated catalog flux density
 * \param epoch   Epoch number
 * \return the new  object.
 */
 SouListElem* 
newSouListElem (odouble ra, odouble dec, ofloat shift[3],
		ofloat gparm[3], ofloat flux)
{
  SouListElem *out=NULL;

  out = g_malloc0(sizeof(SouListElem));
  out->ra        = ra;
  out->dec       = dec;
  out->shift[0]  = shift[0];
  out->shift[1]  = shift[1];
  out->shift[2]  = shift[2];
  out->gparm[0]  = gparm[0];
  out->gparm[1]  = gparm[1];
  out->gparm[2]  = gparm[2];
  out->flux      = flux;
  return out;
} /* end newSouListElem */

/**
 * Update SouListElem
 * \param elem    SouListElem to update
 * \param ra      Catalog celestial position RA (deg)
 * \param dec     Catalog celestial position Dec (deg)
 * \param shift   Shift from reference position
 * \param gparm   Expected gparm in reference image
 * \param flux    Estimated catalog flux density
 * \return the new  object.
 */
 void
SouListElemUpdate (SouListElem *elem, odouble ra, odouble dec, 
		   ofloat shift[3], ofloat gparm[3], ofloat flux)
{
  elem->ra        = ra;
  elem->dec       = dec;
  elem->shift[0]  = shift[0];
  elem->shift[1]  = shift[1];
  elem->shift[2]  = shift[2];
  elem->gparm[0]  = gparm[0];
  elem->gparm[1]  = gparm[1];
  elem->gparm[2]  = gparm[2];
  elem->flux      = flux;
} /* end SouListElemUpdate */

/**
 * Print contents 
 * \param in Object to print
 * \param file  FILE* to write to
 */
 void SouListElemPrint (SouListElem *in, FILE *file)
{
  if (!in) return;

  fprintf (file, "RA=%lf, Dec=%lf\n",in->ra, in->dec);
  fprintf (file, "Position shift (deg) [%f,%f,%f]\n",
	   in->shift[0],in->shift[1],in->shift[2]);
  fprintf (file, "Gaussian parms [%f,%f,%f]\n",
	   in->gparm[0],in->gparm[1],in->gparm[2]);
  fprintf (file, "Flux density %f\n",  in->flux);
} /* end SouListElemPrint */

/**
 * Destructor 
 * \param in Object to delete
 */
 void freeSouListElem (SouListElem *in)
{
  if (in) g_free(in);
} /* end freeSouListElem */



/*  SouList functions */
/**
 * SouList Constructor
 * \return the new  SouList structure.
 */
 SouList* newSouList (void)
{
  SouList *out=NULL;

  out = g_malloc0(sizeof(SouList));
  out->number = 0;
  out->list   = NULL;
  return out;
} /* end newSouListElem */

/**
 * Attach elem to list in
 * \param in   list to add elem to
 * \param elem the element to add.
 */
 void SouListAdd (SouList *in, SouListElem *elem)
{
  /* link to list */
  in->list = g_slist_append (in->list, elem);
  in->number++;
} /* end SouListAdd */

/**
 * Remove elem from list in
 * \param in   list to remove elem from
 * \param elem the element to remove.
 */
 void SouListRemove (SouList *in, SouListElem *elem)
{
  /* remove from list */
  in->list = g_slist_remove(in->list, elem);
  in->number--; /* keep count */  
} /* end SouListRemove  */

/**
 * Remove all elements from list in
 * \param in   list to remove elem from
 * \param elem the element to remove.
 */
 void SouListClear (SouList *in)
{
  GSList *tmp;

  if (in==NULL) return;  /* Does it exist? */
  if (in->list==NULL) return;  /* Anything in it? */

  /* loop through list deleting elements */
  tmp = in->list;
  while (tmp!=NULL) {
    if (tmp->data) freeSouListElem(tmp->data);
    tmp = g_slist_next(tmp);
  }

  /* delete members  */
  g_slist_free(in->list);
  in->list = NULL;
  in->number = 0;

} /* end SouListClear  */

/**
 * Print all elements in list in to file
 * \param in   list to remove elem from
 * \param elem the element to remove.
 */
 void SouListPrint (SouList *in, FILE *file)
{
  GSList *tmp;

  if (in==NULL) return;  /* Does it exist? */
  if (in->list==NULL) return;  /* Anything in it? */

  fprintf (file, "Listing Of SouList\n");

  /* loop through list printing elements */
  tmp = in->list;
  while (tmp!=NULL) {
    if (tmp->data) SouListElemPrint(tmp->data, file);
    tmp = g_slist_next(tmp);
  }

} /* end SouListPrint  */

/**
 * Destructor 
 * \param in Object to delete
 */
 void freeSouList (SouList *in)
{
  /* Clear List */
  SouListClear(in);

  /* delete object */
  g_free(in);
} /* end freeSouListElem */

/**
 * Deconvolves a Gaussian "beam" from a gaussian component.  
 * If the signifigance of the axis fits is <98% or the size is
 * less than a third of the psf, the component is zeroed.
 * \param fmaj   Fitted major axis 
 * \param fmin   Fitted minor axis 
 * \param fpa    Fitted position angle of major axis 
 * \param cmaj   Point source major axis 
 * \param cmin   Point source minor axis 
 * \param cpa    Point source position angle of major axis 
 * \param rmaj   [out] Real major axis; = 0 => unable to fit 
 * \param rmin   [out] Real minor axis; = 0 => unable to fit 
 * \param rpa    [out] Real position angle of major axis 
 * \return Error return: 0 => ok 
 *         1,2-> # components unable to deconvolve 
 */
olong deconv (ofloat fmaj, ofloat fmin, ofloat fpa, 
	      ofloat cmaj, ofloat cmin, ofloat cpa, 
	      ofloat peak, ofloat irms,
	      ofloat *rmaj, ofloat *rmin, ofloat *rpa)
{
  olong   ierr=0;
  ofloat cmj2, cmn2, fmj2, fmn2, sinc, cosc, rhoc, sigic2, 
    det, rhoa, lfpa, lcpa, konst = 28.647888;
  ofloat snr, snrmaj, snrmin, errmaj, errmin;
  /* Confidence-level parameters for "significant" resolution; 
     2.33 => 98% confidence */
  ofloat sigmax = 2.33;

  /* fitted sizes must be at least psf */
  fmaj = MAX (fmaj, cmaj);
  fmin = MAX (fmin, cmin);

  /* Effective SNRs^2 to account for correlated noise. */
  snr = peak / irms;

  /* SNR**2 for major axis error */
  snrmaj = ((fmaj)*(fmin)/(4.0*cmaj*cmaj)) * 
    pow((1.0+(cmaj/(fmaj))*(cmaj/(fmaj))), 2.5) *
    pow((1.0+(cmin/(fmin))*(cmin/(fmin))), 0.5) *
    snr * snr;

  /* SNR**2 for minor axis/PA errors */
  snrmin = ((fmaj)*(fmin)/(4.0*cmaj*cmaj)) * 
    pow((1.0+(cmaj/(fmaj))*(cmaj/(fmaj))), 0.5) *
    pow((1.0+(cmin/(fmin))*(cmin/(fmin))), 2.5) *
    snr * snr;
  
  /* Axis sizes errors include 2% calibration error. */
  errmaj = sqrt (((2.0 * fmaj*fmaj) / snrmaj) +  (0.02*cmaj)*(0.02*cmaj));
  errmin = sqrt (((2.0 * fmin*fmin) / snrmin) +  (0.02*cmin)*(0.02*cmin));

  /* Get useful constants */
  lfpa = fmodf (fpa+900.0, 180.0);
  lcpa = fmodf (cpa+900.0, 180.0);
  cmj2 = cmaj * cmaj;
  cmn2 = cmin * cmin;
  fmj2 = fmaj * fmaj;
  fmn2 = fmin * fmin;
  sinc = (lfpa - lcpa) / konst;
  cosc = cos (sinc);
  sinc = sin (sinc);

  /* Trigonometry now */
  rhoc = (fmj2 - fmn2) * cosc - (cmj2 - cmn2);
  if (rhoc == 0.0) {
    sigic2 = 0.0;
    rhoa = 0.0;
  } else {
    sigic2 = atan((fmj2 - fmn2) * sinc / rhoc);
    rhoa = ((cmj2 - cmn2) - (fmj2 - fmn2) * cosc) / (2.0 * cos (sigic2));
  } 

  (*rpa) = sigic2 * konst + lcpa;
  det = ((fmj2 + fmn2) -(cmj2 + cmn2)) / 2.0;
  (*rmaj) = det - rhoa;
  (*rmin) = det + rhoa;
  ierr = 0;
  if (*rmaj < 0.0) ierr++;
  if (*rmin < 0.0) ierr++;

  /* Swap to get major > minor */
  (*rmaj) = MAX (0.0, *rmaj);
  (*rmin) = MAX (0.0, *rmin);
  (*rmaj) = sqrt (fabs (*rmaj));
  (*rmin) = sqrt (fabs (*rmin));
  if (*rmaj < *rmin) {
    sinc = (*rmaj);
    (*rmaj) = (*rmin);
    (*rmin) = sinc;
    (*rpa) = (*rpa)+90.0;
  } 

  /* Fix up PA */
  (*rpa) = fmodf (*rpa+900.0, 180.0);
  if (*rmaj == 0.0) {
    (*rpa) = 0.0;
  } else if (*rmin == 0.0) {
    if ((fabs(*rpa-lfpa) > 45.0)  &&  (fabs(*rpa-lfpa) < 135.0)) 
      (*rpa) = fmodf (*rpa+450.0, 180.0);
  } 

  /* Is resolution significant? */
  if ((*rmin<sigmax*errmin) || (*rmin<0.33*cmin)) *rmin = 0.0;
  if ((*rmaj<sigmax*errmaj) || (*rmaj<0.33*cmaj)) {*rmaj = 0.0; *rpa = 0.0;}
  return ierr;
} /* end of routine deconv */ 
