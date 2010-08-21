/* $Id$  */
/* FndSou Obit task - generate source list from image                 */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2010                                          */
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

#include <math.h>
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitImage.h"
#include "ObitFArray.h"
#include "ObitImageUtil.h"
#include "ObitHistory.h"
#include "ObitFitRegionList.h"
#include "ObitImageFit.h"
#include "ObitTableUtil.h"
#include "ObitTableMFUtil.h"
#include "ObitTableVLUtil.h"

/* Private structures */
/** IslandElem */
typedef struct { 
  /**  key/index */
  olong index;
  /** Max value */
  ofloat peak;
  /** bottom left corner (0-rel) */
  olong blc[2];
  /** top right corner (0-rel) */
  olong trc[2];
}  IslandElem; 

/** IslandList structure. */  
typedef struct  {
  /** glib singly linked list */
  GSList* list;
  /** How many entries */
  gulong number;
  /** Maximum index */
  gulong maxIndex;
} IslandList; 

/* internal prototypes */
/* Get inputs */
ObitInfoList* FndSouIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void FndSouOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Digest inputs */
void digestInputs(ObitInfoList *myInput, ObitErr *err);
/* Get image from inputs */
void FndSouGetImage(ObitInfoList *myInput, ObitErr *err);
/* Find and fit sources */
void doFndSou (ObitInfoList *myInput, ObitImage *inImage, 
	       ObitImage *outImage, ObitErr *err);
/* Write History */
void doHistory (ObitInfoList *myInput, ObitImage *inImage, 
		ObitImage *outImage, ObitErr *err);
/* Find islands */
IslandList* Islands(ObitFArray *data, ofloat cutt);
/* Convert Islands to Fit regions and fit */
ObitFitRegionList* Island2Region 
(ObitInfoList *myInput, IslandList* island, ObitImage *image, ObitErr *err);
/* Fit a region */
void FitRegion (ObitInfoList *myInput, ObitFitRegion *reg,
		ObitImage *image, olong indx, ObitErr *err);
/* Get initial model for an Island/Region */
 void  GetInitialModel (IslandElem* isElem,  ObitImage *image, gboolean doMulti,
			ofloat cutt, gboolean doPoint, gboolean doPA, 
			olong *nmodel, ObitFitModel *models[], 
			ObitErr *err);
/* Initial single model using moments */
void snglDef (ObitFArray *pixels, gboolean doPoint, gboolean doPA,
	      ofloat dmax, olong ixmax, olong iymax, 
	      ofloat cb[3], ObitFitModel **model) ;
/* Convert moments to Gaussian */
void scdmom (ofloat dmax, ofloat* sum2, ofloat* sumd2, ofloat* sum4, 
	     ofloat* a, ofloat* b, ofloat* theta, gboolean* singul) ;
/* Initial multiple component model */
void multDef (olong peakNo, gboolean doPA, olong* xpk, olong* ypk, 
	      ofloat* spk, ofloat cb[3], ObitFitModel *models[]);
/* Split island in two */
void fitTwo (ObitFArray *pixels, ObitFitRegion *reg, ObitFitRegion *oldreg, 
	     ofloat cb[3], gboolean doPoint, ObitErr* err);

/*------------  IslandElem Function protptypes  ----------*/
/** Private: Create a IslandElem. */
static IslandElem* newIslandElem (olong index, ofloat peak, 
				  olong blc[2], olong trc[2]);

/** Private: Delete a IslandElem. */
static void freeIslandElem (IslandElem *in);

/** Private: Delete a IslandElem. */
static void freeIslandElem (IslandElem *in);

/** Private: Print a IslandElem. */
static void IslandElemPrint (IslandElem *elem, FILE *file);

/*------------  IslandList Function protptypes  ----------*/
/** Private: Create a IslandList. */
static IslandList* newIslandList (void);

/** Private: Delete a IslandList. */
static void freeIslandList (IslandList *in);

/** Private: Print a IslandList. */
static void IslandListPrint (IslandList *in, FILE *file);

/** Private: Append an IslandElem to the list. */
static void IslandListAppend(IslandList *in, IslandElem *elem);

/** Private: Remove an IslandElem from the list. */
static void IslandListRemove (IslandList *in, IslandElem *elem);

/** Private: Find item in a list */
static IslandElem*  IslandListFind(IslandList *in, olong index);

/** Private: Find item in a list with highest peak */
static IslandElem*  IslandListFindHi(IslandList *in);

/** Private: Merge Islands */
static void IslandListMerge(IslandList *in, olong in1, olong in2, 
			    olong nx, olong *prev, olong *curr);
/** Private: Expand Island */
static void IslandListExpand(IslandList *in, olong in1, 
			     ofloat val, olong i, olong j);

/* Program globals */
gchar *pgmName = "FndSou";       /* Program name */
gchar *infile  = "FndSou.in" ;   /* File with program inputs */
gchar *outfile = "FndSou.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitImage *inImage;    /* Input image */
ObitImage *outImage=NULL; /* output image */
olong  prtLv=0;         /* Message level desired */
FILE  *prtFile;        /* Output file */
gboolean doResid=FALSE;/* save Residuals? */
olong nGood=0;         /* Number of fitted components */
olong breakIsland=0;   /* Number of islands broken into multiple */
olong failBreak=0;     /* number of islands failing to break into multiple */
olong rejectLoFlux=0;  /* number of componnts rejected due to low flux */
olong iterLimit=0;     /* number of fits hitting iteration limit */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */


int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/* FndSou Obit program to generate source list from image                 */
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = FndSouIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}
 
  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get input image and output image */
  FndSouGetImage(myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Munge image */
  doFndSou (myInput, inImage, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* History */
  doHistory (myInput, inImage, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* cleanup */
  inImage     = ObitImageUnref(inImage);
  outImage    = ObitImageUnref(outImage);
  myInput     = ObitInfoListUnref(myInput); 
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  if (prtFile && (prtFile!=stdout)) fclose (prtFile);
  
  return ierr;
} /* end of main */

ObitInfoList* FndSouIn (int argc, char **argv, ObitErr *err)
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
  gchar *strTemp=NULL;
  oint    itemp, i, j, k;
  ObitInfoList* list;
  gchar *routine = "FndSouIn";

  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* command line arguments */
  /* fprintf (stderr,"DEBUG arg %d %s\n",argc,argv[0]); DEBUG */
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
      
     } else if (strcmp(arg, "-DataType") == 0) { /* Image type AIPS or FITS */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "DataType", OBIT_string, dim, strTemp);
      
    } else { /* unknown argument */
      Usage();
    }
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
} /* end FndSouIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: FndSou -input file -output ofile [args]\n");
    fprintf(stderr, "FndSou generates a source list from an image\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def FndSou.in\n");
    fprintf(stderr, "  -output output result file, def FndSou.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input image\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    
    /*/exit(1);  bail out */
  }/* end Usage */

/*----------------------------------------------------------------------- */
/*  Create default input ObitInfoList                                     */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*  Values:                                                               */
/*     pgmNumber Int        Program number (like POPS number) def 1       */
/*     nFITS     Int        Number of FITS directories [def. 1]           */
/*     FITSdirs  Str [?,?]  FITS directories [def {"./"}]                 */
/*     AIPSuser  Int        AIPS user number [def 2}]                     */
/*     nAIPS     Int        Number of AIPS directories [def. 1]           */
/*     AIPSdirs  Str [?,?]  AIPS directories [def std. AIPS]              */
/*     DataType  Str [4]    "AIPS" or "FITS" [def {"FITS"}]               */
/*     inFile    Str [?]    input FITS image file name [no def]           */
/*     inName    Str [12]   input AIPS image name  [no def]               */
/*     inClass   Str [6]    input AIPS image class  [no def]              */
/*     inSeq     Int        input AIPS image sequence no  [no def]        */
/*     inDisk    Int        input AIPS or FITS image disk no  [def 1]     */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp=NULL;
  ofloat ftemp;
  oint   itemp;
  ObitInfoList *out = newObitInfoList();
  gchar tname[51];
  gchar *routine = "defaultInputs";

  /* add parser items */
  /* Program number */
  dim[0] = 1; dim[1] = 1;
  itemp = 1;
  ObitInfoListPut (out, "pgmNumber", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Default number of FITS directories  */
  dim[0] = 1; dim[1] = 1;
  itemp = 0; /* number of FITS directories */
  ObitInfoListPut (out, "nFITS", OBIT_oint, dim, &itemp, err);
  /* If fitsdirs is not defined then $FITS, $FITS01... will be used */

  /* AIPS user number */
  dim[0] = 1; dim[1] = 1;
  itemp = 2;
  ObitInfoListPut (out, "AIPSuser", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Default number of AIPS directories */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; /* number of AIPS directories */
  ObitInfoListPut (out, "nAIPS", OBIT_oint, dim, &itemp, err);
  /* If aipsdirs is not defined then $DA01, $DA02... will be used */

  /* Default type "FITS" */
  strTemp = "FITS";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "DataType", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  strTemp = "undefn";

  /* input FITS file name */
  g_snprintf (tname, 50, "inFile");
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, tname, OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* input AIPS file name */
  g_snprintf (tname, 50, "inName");
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, tname, OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* input AIPS file class */
  g_snprintf (tname, 50, "inClass");
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* AIPS sequence */
  g_snprintf (tname, 50, "inSeq");
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, tname, OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* AIPS or FITS disk number */
  g_snprintf (tname, 50, "inDisk");
  dim[0] = 1;dim[1] = 1;
  ftemp = 1; 
  ObitInfoListPut (out, tname, OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Default output */
  /* FITS file name */
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* AIPS file name */
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* AIPS file class */
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* AIPS sequence */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "outSeq", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* AIPS or FITS disk number */
  dim[0] = 1;dim[1] = 1;
  ftemp = 1; 
  ObitInfoListPut (out, "outDisk", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  return out;
} /* end defaultInputs */

/*----------------------------------------------------------------------- */
/*  Create default output ObitInfoList                                    */
/*  Nothing for this program                                              */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*  Values:                                                               */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultOutputs(ObitErr *err)
{
  /*gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};*/
  /* ofloat ftemp;*/
  ObitInfoList *out = newObitInfoList();
  gchar *routine = "defaultOutputs";

  if (err->error) return out;  /* existing error */

  /* add parser items */
  /* Image mean */
  /* dim[0] = 1; dim[1] = 1; */
  /* ftemp = 0.0; */
  /* ObitInfoListPut (out, "mean", OBIT_float, dim, &ftemp, err); */
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  return out;
} /* end defaultOutputs */

/*----------------------------------------------------------------------- */
/*  Digest inputs                                                         */
/*   Sets globals doResid, prtLv and prtFile                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void digestInputs(ObitInfoList *myInput, ObitErr *err)
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar pfile[256];
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* prtLv */
  ObitInfoListGetTest(myInput, "prtLv", &type, dim, &prtLv);

  /* OutPrint/prtFile */
  sprintf (pfile, "          ");
  ObitInfoListGetTest(myInput, "OutPrint", &type, dim, pfile);
  if (strncmp (pfile, "        ", 8)) {
    ObitTrimTrail(pfile);  /* Trim any trailing blanks */
    prtFile = fopen (pfile, "a");
  } else prtFile=stdout;

  /* doResid */
  ObitInfoListGetTest(myInput, "doResid", &type, dim, &doResid);

} /* end digestInputs */

/*----------------------------------------------------------------------- */
/*  Get images from myInput                                               */
/*  Values:                                                               */
/*      myImage   InfoList with inputs                                    */
/*      inImage   [out] ObitImage pointer                                 */
/*                 Passed as global                                       */
/*      outImage  [out] Output ObitImage pointer                          */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void FndSouGetImage(ObitInfoList *myInput, ObitErr *err)
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong         j, k, Aseq, disk, cno;
  gboolean     exist;
  gchar        *strTemp=NULL, inFile[128];
  gchar        iAname[13], iAclass[7];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gchar        tname[101];
  gchar *routine = "FndSouGetImage";

  if (err->error) return;  /* existing error? */

  /* Get region from myInput */
  ObitInfoListGetTest(myInput, "BLC", &type, dim, blc); /* BLC */
  ObitInfoListGetTest(myInput, "TRC", &type, dim, trc); /* TRC */

  /* File type - could be either AIPS or FITS */
  ObitInfoListGet (myInput, "DataType", &type, dim, tname, err);
  if (err->error) Obit_traceback_msg (err, routine, routine);
  if (!strncmp (tname, "AIPS", 4)) { /* AIPS input */

    /* input AIPS disk */
    g_snprintf (tname, 100, "inDisk");
    ObitInfoListGet(myInput, tname, &type, dim, &disk, err);

    /* input AIPS name */
    g_snprintf (tname, 50, "inName");
    for (k=0; k<12; k++) iAname[k] = ' '; iAname[k] = 0;
    ObitInfoListGet(myInput, tname, &type, dim, iAname, err);

    /* input AIPS class */
    g_snprintf (tname, 50, "inClass");
    for (k=0; k<6; k++) iAclass[k] = ' '; iAclass[k] = 0;
    ObitInfoListGet(myInput, tname, &type, dim, iAclass, err);

    /* input AIPS sequence */
    g_snprintf (tname, 50, "inSeq");
    ObitInfoListGet(myInput, tname, &type, dim, &Aseq, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    /* if ASeq==0 want highest existing sequence */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, iAname, iAclass, Atype, TRUE, err);
      if (err->error) Obit_traceback_msg (err, routine, "myInput");
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, tname, OBIT_oint, dim, &Aseq);
    } 
    
    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, iAname, iAclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /* Generate Object name from AIPS name */
    g_snprintf (tname, 100, "%s.%s:%d.%d", iAname, iAclass, Aseq, disk);
    inImage = newObitImage(tname);
    
    /* define image */
    ObitImageSetAIPS (inImage, OBIT_IO_byPlane, disk, cno, AIPSuser, 
		      blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /* Make sure it's OK */
    ObitImageFullInstantiate (inImage, TRUE, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);

  } else if (!strncmp (tname, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    for (j=0; j<128; j++) inFile[j] = 0;
        g_snprintf (tname, 100, "inFile");
    ObitInfoListGet(myInput, tname, &type, dim, inFile, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /* input FITS disk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /*  Object name from FITS name */
    inImage = newObitImage(inFile);
    
    /* define image */
    ObitImageSetFITS (inImage, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /* Make sure it's OK */
    ObitImageFullInstantiate (inImage, TRUE, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
  
  } /* end FITS input */
    
    /* Output image only for doResid */
  if (doResid) {
    /* File type - could be either AIPS or FITS */
    ObitInfoListGet (myInput, "DataType", &type, dim, tname, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    if (!strncmp (tname, "AIPS", 4)) { /* AIPS input */
      /* AIPS disk */
      ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
      /* AIPS name */
      for (k=0; k<12; k++) Aname[k] = ' '; Aname[k] = 0;
      ObitInfoListGet(myInput, "outName", &type, dim, Aname, err);
      Aname[dim[0]] = 0;
      /* Default */
      if (!strncmp(Aname,"     ",5)) strcpy (Aname, iAname);
      
      /* AIPS class */
      for (k=0; k<6; k++) Aclass[k] = ' '; Aclass[k] = 0;
      ObitInfoListGet(myInput, "outClass", &type, dim, Aclass, err);
      Aclass[dim[0]] = 0;
      /* Default */
      if (!strncmp(Aclass,"      ",6)) strcpy (Aclass, "Resid ");
      
      /* AIPS sequence */
      ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
      
      /* if ASeq==0 create new, high+1 */
      if (Aseq<=0) {
	Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
	if (err->error) Obit_traceback_msg (err, routine, "myInput");
	/* Save on myInput*/
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
      } 
      
      /* Find catalog number */
      cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, 
			     &exist, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
      
      /* Generate Object name from AIPS name */
      g_snprintf (tname, 100, "%s.%s:%d.%d", Aname, Aclass, Aseq, disk);
      outImage = newObitImage(tname);
      
      /* reset BLC, TRC */
      for (j=0; j<IM_MAXDIM; j++)  {blc[j] = 1; trc[j] = 0;}
      
      /* Tell about it */
      Obit_log_error(err, OBIT_InfoErr, "Output AIPS image %s %s %d on disk %d cno %d",
		     Aname, Aclass, Aseq, disk, cno);

      /* define image */
      ObitImageSetAIPS (outImage, OBIT_IO_byPlane, disk, cno, AIPSuser, 
			blc, trc, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
      
    } else if (!strncmp (tname, "FITS", 4)) {  /* FITS input */
      
      /* Output image */ 
      /* FITS file name */
      for (j=0; j<128; j++) inFile[j] = 0;
      ObitInfoListGet(myInput, "outFile", &type, dim, inFile, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
      ObitTrimTrail(inFile);  /* remove trailing blanks */
     
      /*  FITS disk */
      ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
      
      /*  Object name from FITS name */
      outImage = newObitImage(inFile);
      
      /* reset BLC, TRC */
      for (j=0; j<IM_MAXDIM; j++)  {blc[j] = 1; trc[j] = 0;}
      
      /* Give output Image name */
      Obit_log_error(err, OBIT_InfoErr, "Output FITS image %s on disk %d ",
		     inFile, disk);

     /* define image */
      ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
      
    } else { /* Unknown type - barf and bail */
      Obit_log_error(err, OBIT_Error, "%s: Unknown Image type %s", 
		     pgmName, strTemp);
    } /* end of output image */
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /* Clone output from input if needed */
    ObitImageClone ( inImage, outImage, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
    /* Ensure outImage fully instantiated and OK */
    ObitImageFullInstantiate (outImage, FALSE, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  } /* end if output needed */
} /* end FndSouGetImage */

/*----------------------------------------------------------------------- */
/*  Write history                                                         */
/*  Values:                                                               */
/*      myInputs  input parameter array                                   */
/*      inImage   Array of ObitImage pointers                             */
/*      outImage  Output ObitImage pointer                                */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void doHistory (ObitInfoList *myInput, ObitImage *inImage, 
		ObitImage *outImage, ObitErr *err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", "inFile",
    "BLC", "TRC", "doVL", "doResid", "NGauss", "CutOff",
    "Retry", "Sort", "doMult", "doWidth", "Gain", "Parms", "RMSsize",
    NULL};
  gchar *routine = "doHistory";

  if (err->error) return;  /* existing error? */

  /* Add history to input */
  inHistory  = newObitDataHistory ((ObitData*)inImage, OBIT_IO_ReadWrite, err);
  /* Add this programs history */
  ObitHistoryOpen (inHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (inHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (inHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  ObitHistoryClose (inHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  
  /* add history to output if any */
  if (!doResid || !outImage) return;

  /* Open outImage to make sure History gets registered in header if created */
  ObitImageOpen (outImage, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Do history  */
  inHistory  = newObitDataHistory ((ObitData*)inImage, OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory ((ObitData*)outImage, OBIT_IO_WriteOnly, err);
  ObitHistoryCopyHeader (inHistory, outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  
  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  
  /* Close outImage to update header */
  outImage->myStatus = OBIT_Modified;
  ObitImageClose (outImage, err);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end doHistory */
/*----------------------------------------------------------------------- */
/*  Generate source list, make residual image                             */
/*  Values:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inImage   input ObitImage                                         */
/*      convFn    FndSouving Function                                     */
/*      outImage  Output ObitImage pointer                                */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void doFndSou (ObitInfoList *myInput, ObitImage *inImage, 
	       ObitImage *outImage, ObitErr *err)
{
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitFArray *data=NULL;
  olong  i, plane[5] = {1,1,1,1,1};
  olong ipass, npass;
  ofloat cutt, rms;
  gboolean doVL;
  IslandList* islands=NULL;
  ObitFitRegionList *regList=NULL;
  ObitTableMF *TableMF=NULL;
  ObitTableVL *TableVL=NULL;
  olong iver, oldnGood=0, oldbreakIsland=0, oldfailBreak=0, oldrejectLoFlux=0, 
    olditerLimit=0, nIslands=0;
  gchar Sort[3];
  gchar *colName[3]  = {"DELTAX", "DELTAY", "FLUX"};
  gchar *ImgParms[] = {  
    "doPBCorr", "asize",  /* Primary beam correction */
    "RMSsize",            /* RMS box half width */
    NULL};
  gchar      *routine = "doFndSou";
      
  if (err->error) return;  /* existing error? */
  
  /* Get/Create output MF table */
  iver = 1;
  TableMF = newObitTableMFValue (inImage->name, (ObitData*)inImage, &iver, 
				 OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;

  /* Get/Create output VL table if requested  */
  doVL = FALSE;
  ObitInfoListGetTest(myInput, "doVL", &type, dim, &doVL);
  if (doVL) {
    iver = 1;
    TableVL = newObitTableVLValue (inImage->name, (ObitData*)inImage, &iver, 
				   OBIT_IO_ReadWrite, err);
    if (err->error) goto cleanup;
  }

  /* Read and extract image data*/
  ObitImageGetPlane (inImage, NULL, plane, err);
  data = ObitFArrayRef(inImage->image);
  if (err->error) goto cleanup;

  /* Get RMS from histogram/ save */
  rms = ObitFArrayRMS(data);
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inImage->info, "RMS", OBIT_float, dim, &rms);

  /* Get minimum for defining islands */
  cutt = 0.0;
  ObitInfoListGetTest(myInput, "CutOff", &type, dim, &cutt);
  if (cutt<=0.0) {
    /* Default to 5 sigma */
    cutt = 5.0 * rms;
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut(myInput, "CutOff", OBIT_float, dim, &cutt);
  }

  /* Loop npass times */
  npass = 1;
  ObitInfoListGetTest(myInput, "NPass", &type, dim, &npass);
  npass = MAX (1, npass);
  for (ipass = 1; ipass<=npass; ipass++) {
    /* Locate islands */
    islands = Islands(data, cutt);
    Obit_log_error(err, OBIT_InfoErr, "Found %d islands pass %d", 
		   (olong)islands->number, ipass);
    ObitErrLog(err); 
    nIslands += islands->number;  /* total number of islands */
    
    /* Listing of islands requested? */
    if (prtLv>=2)  IslandListPrint (islands, prtFile);
    
    /* Loop converting  Islands to regions with initial model the fit
       After fitting region, if the residuals are too high (>RESMAX) 
       then attempt to break in two and refit.
       Take the better of the original or second fit.
       Then finally subtract region fit from the image attached to inImage */
    /* Convert Islands to regions with initial model */
    regList = Island2Region (myInput, islands, inImage, err);
    if (err->error) goto cleanup;
    
    /* Tell about things */
    if (npass==1) {
      Obit_log_error(err, OBIT_InfoErr, "Successfully fitted %d components",  nGood);
      Obit_log_error(err, OBIT_InfoErr, "Attempt to break %d islands into multiple", 
		     breakIsland);
      Obit_log_error(err, OBIT_InfoErr, " %d Attempts to break islands failed", 
		     failBreak);
      Obit_log_error(err, OBIT_InfoErr, " %d components rejected for low peak", 
		     rejectLoFlux);
      Obit_log_error(err, OBIT_InfoErr, " %d fits hit iteration limit", iterLimit);
    } else {  /* Multiple passes */
      Obit_log_error(err, OBIT_InfoErr, "Pass %d: Successfully fitted %d components",  
		     ipass, nGood-oldnGood);
      Obit_log_error(err, OBIT_InfoErr, "   Attempt to break %d islands into multiple", 
		     breakIsland-oldbreakIsland);
      Obit_log_error(err, OBIT_InfoErr, "   %d Attempts to break islands failed", 
		     failBreak-oldfailBreak);
      Obit_log_error(err, OBIT_InfoErr, "   %d components rejected for low peak", 
		     rejectLoFlux-oldrejectLoFlux);
      Obit_log_error(err, OBIT_InfoErr, "   %d fits hit iteration limit", 
		     iterLimit-olditerLimit);
      oldnGood        = nGood;
      oldbreakIsland  = breakIsland;
      oldfailBreak    = failBreak;
      oldrejectLoFlux = rejectLoFlux;
      olditerLimit    = iterLimit;
    }
    ObitErrLog(err); 

    /* Convert RegionFitList to MF table entries */
    ObitTableMFRegions2MF (TableMF, regList, inImage, err);
    if (err->error) goto cleanup;

    regList = ObitFitRegionListUnref(regList);
    freeIslandList(islands);
  } /* End pass loop */

  /* Write residuals if requested */
  if (doResid && outImage) {
    ObitImagePutPlane (outImage, data->array, plane, err);
    if (err->error) goto cleanup;
  }

  /* Multipass summary */
  if (npass>1) {
    Obit_log_error(err, OBIT_InfoErr, "Total: Successfully fitted %d components",  nGood);
    Obit_log_error(err, OBIT_InfoErr, "Total: Attempt to break %d islands into multiple", 
		   breakIsland);
    Obit_log_error(err, OBIT_InfoErr, "Total: %d Attempts to break islands failed", 
		   failBreak);
    Obit_log_error(err, OBIT_InfoErr, "Total: %d components rejected for low peak", 
		   rejectLoFlux);
    Obit_log_error(err, OBIT_InfoErr, "Total: %d fits hit iteration limit", iterLimit);
    ObitErrLog(err); 
  }

 /* Sort MF table*/
  Sort[0] = ' ';
  ObitInfoListGetTest(myInput, "Sort", &type, dim, Sort);
  if (Sort[0]=='Y') i = 1;
  else if (Sort[0]=='S')  i = 2;
  else i=0;
  ObitTableUtilSort ((ObitTable*)TableMF, colName[i], FALSE, err);
  if (err->error) goto cleanup;

  if (doVL) {
    /* Primary beam correction? conversion parameters from myInput */
    ObitInfoListCopyList (myInput, inImage->info, ImgParms);
    if (err->error) goto cleanup;

    /* Convert MF to VL table */
    ObitTableMF2VL (TableMF, TableVL, inImage, err);

    /* Sort/index VL table */
    ObitTableVLIndex(TableVL, err);
    if (err->error) goto cleanup;
    
    /* List VL table if requested */
    if (prtLv>=1) {
      fprintf (prtFile, "Found %d islands in %d passes\n", nIslands, npass);
      fprintf (prtFile, "Successfully fitted %d components\n", nGood);
      fprintf (prtFile, "Attempt to break %d islands into multiple\n",  breakIsland);
      fprintf (prtFile, " %d Attempts to break islands failed\n", failBreak);
      fprintf (prtFile, " %d components rejected for low peak\n", rejectLoFlux);
      fprintf (prtFile, " %d fits hit iteration limit\n\n", iterLimit);
      ObitTableVLPrint (TableVL, inImage, prtFile, err);
      if (err->error) goto cleanup;
    }
    
  } else {
    /* List MF table if requested */
    if (prtLv>=1) {
      fprintf (prtFile, "Found %d islands\n", nIslands);
      fprintf (prtFile, "Successfully fitted %d components\n", nGood);
      fprintf (prtFile, "Attempt to break %d islands into multiple\n",  breakIsland);
      fprintf (prtFile, " %d Attempts to break islands failed\n", failBreak);
      fprintf (prtFile, " %d components rejected for low peak\n", rejectLoFlux);
      fprintf (prtFile, " %d fits hit iteration limit\n\n", iterLimit);
      ObitTableMFPrint (TableMF, inImage, prtFile, err);
      if (err->error) goto cleanup;
    }
  }

 cleanup:
  data = ObitFArrayUnref(data);
  regList = ObitFitRegionListUnref(regList);
  TableMF = ObitTableMFUnref(TableMF);
  TableVL = ObitTableVLUnref(TableVL);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  
} /* end doFndSou */

/**
 * Locate contigious areas in data above cutt and write into island list
 * Adopted from Walter Jaffes's version in AIPS task SAD
 * \param data  Image pixel array
 * \param cutt  Minimum value in data
 * \return the IslandList
 */
/* Find islands */
IslandList* Islands(ObitFArray *data, ofloat cutt)
{
  IslandList* islands = NULL;
  IslandElem* elem=NULL;
  olong i, j, k, nx, ny, *prev=NULL, *curr=NULL;
  olong blc[2], trc[2];
  ofloat *row, fblank = ObitMagicF();
  olong pos[2];

  /* Create output */
  islands = newIslandList();

  /* Arrays of previous and current island numbers */
  nx = data->naxis[0];
  ny = data->naxis[1];
  prev = g_malloc0 (nx*sizeof(olong));
  curr = g_malloc0 (nx*sizeof(olong));

  /* Loop over array */
  for (j=0; j<ny; j++) {

    /* This row's data */
    pos[0] = 0; pos[1] = j;
    row = ObitFArrayIndex (data, pos);
  
    /* Loop over row */
    for (i=0; i<nx; i++) {
      /* check if even close to fblank - rounding can screw up 32 bit integers */
      if (fabs(row[i]-fblank)<0.00001*fblank) continue;

      /* Is point above cutoff? */
      if ((row[i] < cutt)  ||  (row[i] == fblank)) {
	curr[i] = 0;  /* not interesting */
	/* Are any adjacent points, on  currnt Line or previous line 
	   already marked? */
      } else if ((i > 0)  &&  (curr[i-1] > 0)) {
	curr[i] = curr[i-1];
	
	/* Bridge to previous entry? */
	if ((i < (nx-1))  &&  (prev[i+1] > 0)) {
	  IslandListMerge (islands, curr[i], prev[i+1], nx, prev, curr);
	}
      } else if ((i > 0)  &&  (prev[i-1] != 0)) {
	curr[i] = prev[i-1];
	/* Is this a link between two  previous distinct islands? */
	if ((i < (nx-1))  &&  (prev[i+1] > 0)) {
	  IslandListMerge (islands, curr[i], prev[i+1], nx, prev, curr);
	}
      } else if (prev[i] != 0) {
	curr[i] = prev[i];
      } else if ((i < (nx-1))  &&  (prev[i+1] != 0)) {
	curr[i] = prev[i+1];
	
	/* Totally new island */
      } else {
	if (row[i]!=fblank) {
	  curr[i] = islands->maxIndex + 1;
	  blc[0] = i; blc[1] = j;
	  trc[0] = i; trc[1] = j;
	  elem = newIslandElem (curr[i], row[i], blc, trc);
	  IslandListAppend (islands, elem); /* add to list*/
	}
      }

      /* New addition to old island */
      if ((curr[i] != 0) && (row[i]!=fblank)) {
	IslandListExpand (islands, curr[i], row[i], i, j);
      }
    } /* end loop over row */

    /* Get ready for next line. */
    for (k=0; k<nx; k++) prev[k] = curr[k];
  
  } /* end loop over array */

  /* cleanup */
  /* cleanup:*/
  if (prev) g_free(prev);
  if (curr) g_free(curr);

  return islands;
} /* end Islands */

/**
 * Convert an IslandList to a FitRegionList
 * Each region has the name "regnnnnnn" where nnnnnn is the 1-rel
 * number of the entry with leading zeroes.
 * Estimates initial model from a moment analysis of the image.
 * Output list in descending order of peak flux density.and all 
 * entries will be removed from the input list.
 * \param myInput Inputs object
 * \param island  List to convert 
 * \param image   ObitImage being described
 * \param err     Obit error/message stack object.
 * \return the new  FitRegionList object.
 */
ObitFitRegionList* Island2Region (ObitInfoList *myInput, IslandList* island, 
				  ObitImage *image, ObitErr *err)
{
  ObitFitRegionList* out=NULL;
  IslandElem*  isElem=NULL;
  ObitFitRegion* reg=NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat cutt, peakResid=0.0, RMSResid=0.0, fluxResid=0.0;
  gboolean doPoint, doPA, doMult;
  olong nisland, corner[2], adim[2], k, nx, ny, i, bsize, d[2], nmodel=1;
  gchar *regname=NULL;
  ObitFitModel *models[10];
  gchar *routine = "Island2Region";

  if (err->error) return out;

  /* Create output list */
  out = ObitFitRegionListCreate ("RegionList", image);
  nx = image->image->naxis[0];
  ny = image->image->naxis[1];
  d[0] = nx-1; d[1] = ny-1;

  /* Allow multiple components? */
  doMult = FALSE;
  ObitInfoListGetTest(myInput, "doMult", &type, dim, &doMult);

  /* Maximum numper of islands? */
  nisland = 100;
  ObitInfoListGetTest(myInput, "NGauss", &type, dim, &nisland);

  /* Get minimum for defining islands */
  cutt = 0.0;
  ObitInfoListGetTest(myInput, "CutOff", &type, dim, &cutt);

   /* Not fitting size? */
  doPoint = TRUE;
  ObitInfoListGetTest(myInput, "doWidth", &type, dim, &doPoint);
  doPoint = !doPoint;

   /* Fitting PA? */
  doPA = doPoint;

  /* Beam size in pixels */
  bsize = (olong)((image->myDesc->beamMaj / 
		  fabs (image->myDesc->cdelt[0]))+0.5);
  bsize = MAX (bsize, 3);

  /* Loop over Island list converting */
  k = 0;
  isElem = IslandListFindHi (island);
  while (isElem) {

    k++;  /* region number */

    /* Done enough? */
    if (k>nisland) break;

    /* Expand window */
    for (i=0; i<2; i++) {
      isElem->blc[i] = MAX (0, isElem->blc[i]-bsize);
      isElem->trc[i] = MIN (d[i],  isElem->trc[i]+bsize);
   }
    
    corner[0] = isElem->blc[0]; corner[1] = isElem->blc[1];
    adim[0] = isElem->trc[0]-isElem->blc[0]+1; 
    adim[1] = isElem->trc[1]-isElem->blc[1]+1; 

    /* Get initial model */
    nmodel = 10;
    GetInitialModel (isElem, image, doMult, cutt, doPoint, doPA, 
		     &nmodel, models, err);
    if (err->error) Obit_traceback_val (err, routine, image->name, out);

    /* Find anything? */
    if (models[0]->Peak!=0.0) {
    
      /* Create region */
      regname   = ObitFitRegionName(k);
      reg = ObitFitRegionCreate (regname, corner, adim, isElem->peak,
				 peakResid, RMSResid, fluxResid, nmodel, models);
      if (regname) g_free(regname);
      
      /* Fit */
      FitRegion (myInput, reg, image, isElem->index, err);
      if (err->error) Obit_traceback_val (err, routine, image->name, out);
      
      /* Add to output list */
      ObitFitRegionListAppend (out, reg);
    } /* End found something to fit */
 
   /* Remove from IslandList */
    IslandListRemove (island, isElem);
    /* Get next */
    isElem = IslandListFindHi (island);
  } /* end loop over island list */

  return out;
} /*  end  Island2Region */

/**
 * Determine initial model for a region
 * Each region has the name "regnnnnnn" where nnnnnn is the 1-rel
 * number of the entry with leading zeroes.
 * Estimates initial model from a moment analysis of the image.
 * Output list in descending order of peak flux density.and all 
 * entries will be removed from the input list.
 * Adapted from the AIPSish VSAD.FOR/SADDAT
 * \param isElem  Island element to be converted
 * \param image   ObitImage with image data array attached
 * \param doMult  If true allow multiple components
 * \param cutt    Lower bound of flux density
 * \param doPoint True if not fitting size
 * \param doPA    True if fitting position angle
 * \param nmodel  Number of models. on input the max allowed
 *                on output, the actual,
 * \param models  Array of ObitFitModels
 *                If peak flux density==0 then nothing found.
 * \param err    Obit error/message stack object.
 */
 void  GetInitialModel (IslandElem* isElem,  ObitImage *image, gboolean doMult,
			ofloat cutt, gboolean doPoint, gboolean doPA, 
			olong *nmodel, ObitFitModel *models[], 
			ObitErr *err)
{
  olong blc[2], trc[2], pos[2];
  ObitFArray *pixels=NULL;
  ofloat *pixData;
  ofloat fblank = ObitMagicF();
  olong  maxPk=*nmodel, mmaxPk, nmpk;
  olong i, j, ind, ix, iy, idx, idy, nx, ny, peakNo, ipts, ipt, kmpk;
  gboolean noGo;
  olong *xxPk=NULL, *yyPk=NULL, *xPk=NULL, *yPk=NULL;
  ofloat *ssPk=NULL,*sPk=NULL, smax, cbeam[3];
  gchar *routine = "GetInitialModel";

  /* Get pixels and get initial model */
  blc[0] = isElem->blc[0]; blc[1] = isElem->blc[1]; 
  trc[0] = isElem->trc[0]; trc[1] = isElem->trc[1]; 
  pixels = ObitFArraySubArr (image->image, blc, trc, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);

 /* Pointer to data */
  pos[0] = pos[1] = 0;
  pixData = ObitFArrayIndex (pixels, pos);
  nx = pixels->naxis[0]; 
  ny = pixels->naxis[1]; 

  /* Allocate work arrays */
  mmaxPk = 10 * maxPk;  /* Should always be big enough */
  xPk  = g_malloc0(mmaxPk*sizeof(olong));
  yPk  = g_malloc0(mmaxPk*sizeof(olong));
  sPk  = g_malloc0(mmaxPk*sizeof(ofloat));
  xxPk = g_malloc0(mmaxPk*sizeof(olong));
  yyPk = g_malloc0(mmaxPk*sizeof(olong));
  ssPk = g_malloc0(mmaxPk*sizeof(ofloat));
  peakNo = 0;

  /* Find Peak if single */
  if (!doMult) {
    ssPk[peakNo] = ObitFArrayMax (pixels, pos);
    xxPk[peakNo] = 1.0 + pos[0];
    yyPk[peakNo] = 1.0 + pos[1];
    peakNo++;
  }

  /* Multiple peaks allowed */
  if (doMult) {
    
    /*  Loop over pixels excluding edges */
    for (iy= 2; iy<= ny-1; iy++) { /* loop 80 */
      for (ix= 2; ix<= nx - 1; ix++) { /* loop 70 */
	/* Quit if there are already maxPk */
	if (peakNo < mmaxPk) {
	  /*  Position in pixData */
	  ipts = (iy-1)*nx + ix -1;
	  /* Only count points above cutt */
	  if (pixData[ipts] <  cutt)   continue;
	  if (pixData[ipts] == fblank) continue;
		  
	  /*  Bigger than surrounding points? */
	  noGo = FALSE;
	  for (idy= -nx; idy<= nx; idy+=nx) { /* loop 60 */
	    for (idx= -1; idx<= 1; idx++) { /* loop 50 */
	      ipt = ipts + idy + idx;
	      /* Jump out if not a maximum */
	      if (pixData[ipt] != fblank) 
		if (pixData[ipts] < pixData[ipt]) noGo = TRUE; /*goto L68;*/
	      if (noGo) break;  /* This one OK? */
	    } /* end loop  L50:  */
	      if (noGo) break;  /* This one OK? */
	  } /* end loop  L60:  */
	} else {
	  /*  Quit if there are already mmaxPk */
	  goto foundEnough; /* L90 */
	} 
	if (noGo) continue;  /* This one OK? */

	/* This is an acceptable maximum but check that there's no adjacent 
	   point (can happen with 2 equal points) */
	for (kmpk=0; kmpk<peakNo; kmpk++) { /* loop 65 */
	  if ((abs(xxPk[kmpk]-ix) <= 1)  ||  (abs(yyPk[kmpk]-iy) <= 1)) 
	    noGo = TRUE; /*goto L68;*/
	} /* end loop  L65:  */
	if (noGo) continue;  /* This one OK? */

	/* Record position and flux */
	xxPk[peakNo] = ix;
	yyPk[peakNo] = iy;
	ssPk[peakNo] = pixData[ipts];
	peakNo++;
	/* L68: Jump here if you have good reason to believe this isn't a maximum */
      } /* end loop  L70:  */
    } /* end loop  L80:  */
  } /* end if multiple peaks */
			 
  /* Jump here if we're done */
foundEnough:  /* L90 */

  /* Take up to maxPk of the brightest  */
  nmpk = MIN (peakNo, maxPk);
  peakNo = 0;
  for (i=0; i<nmpk; i++) { /* loop 150 */
    /* Find brightest remaining */
    smax = fabs (ssPk[0]);
    ind = 0;
    for (j= 1; j<nmpk; j++) { /* loop 120 */
      if (fabs(ssPk[j])  >  smax) {
	ind = j;
	smax = fabs (ssPk[j]);
      } 
    } /* end loop  L120: */
    xPk[peakNo] = xxPk[ind];
    yPk[peakNo] = yyPk[ind];
    sPk[peakNo] = ssPk[ind];
    peakNo++;
    /* Drop this one */
    ssPk[ind] = 0.0;
  } /* end loop  L150: */

  /* Clean beam in pixels */
  cbeam[0] = image->myDesc->beamMaj / fabs (image->myDesc->cdelt[0]);
  cbeam[1] = image->myDesc->beamMin / fabs (image->myDesc->cdelt[0]);
  cbeam[2] = image->myDesc->beamPA * DG2RAD;

  /* For single sources use 2nd moments for starting estimates */
  if ((!doMult)  ||  (peakNo <= 1)) {
    *nmodel = 1;
    snglDef (pixels, doPoint, doPA, sPk[0], xPk[0], yPk[0], cbeam, &models[0]);

    /* For mutiple peaks use several points */
  } else {
    *nmodel = peakNo;
    multDef (peakNo, doPA, xPk, yPk, sPk, cbeam, models);
  } 

  pixels = ObitFArrayUnref(pixels); /* Don't need further */
  

  /* cleanup: */
  if (xPk) g_free(xPk);
  if (yPk) g_free(yPk);
  if (sPk) g_free(sPk);
  if (xxPk) g_free(xxPk);
  if (yyPk) g_free(yyPk);
  if (ssPk) g_free(ssPk);

} /* end GetInitialModel */

/**
 *  Loop through region list fitting.  After fitting region, if the residuals
 *  are too high (>RESMAX) then attempt to break in two and refit.
 *  Take the better of the original or second fit.
 *  Then finally subtract region fit from the image attached to inImage
 * \param myInput Inputs object
 * \param reg     Region to fit
 * \param image   ObitImage being described
 * \param indx    Island number
 * \param err     Obit error/message stack object.
 */
void FitRegion (ObitInfoList *myInput, ObitFitRegion *reg,
		 ObitImage *image, olong indx, ObitErr *err)
{
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, j, ier;
  ObitInfoType type;
  ObitFitRegion *oldReg=NULL;
  ObitImageFit* fitter = NULL;
  gboolean doMult, doPoint, revert;
  odouble dtemp;
  ofloat gain, parms[20], icut, tcut, rcut, xcut, cbeam[3], oldRMS, oldPeak;
  olong blc[2], trc[2];
  ObitFArray *pixels=NULL;
  gchar *routine = "FitRegions";

  if (err->error) return;

  /* Allow multiple components? */
  doMult = FALSE;
  ObitInfoListGetTest(myInput, "doMult", &type, dim, &doMult);

   /* Not fitting size? */
  doPoint = TRUE;
  ObitInfoListGetTest(myInput, "doWidth", &type, dim, &doPoint);
  doPoint = !doPoint;

  /* Control */
  xcut = 1.0;
  ObitInfoListGetTest(myInput, "CutOff", &type, dim, &xcut);
  icut = xcut;
  ObitInfoListGetTest(myInput, "Retry", &type, dim, &icut);
  if (icut==0) icut = xcut;
  for (i=0; i<20; i++) parms[i] = 0;
  ObitInfoListGetTest(myInput, "Parms", &type, dim, parms);
  rcut = parms[0];
  gain = 0.05;
  ObitInfoListGetTest(myInput, "Gain", &type, dim, &gain);

  /* Create fitter */
  fitter = ObitImageFitCreate ("my Fitter");

  /* Print level */
  i = prtLv-2;
  ObitInfoListAlwaysPut (fitter->info, "prtLv",  OBIT_long, dim, &i);

  /* Add bounds */
  if (parms[0]>0.0) {
    dtemp = 0.0;
    ObitInfoListAlwaysPut (fitter->info, "FluxLow",  OBIT_double, dim, &dtemp);
  }
  if (parms[1]>0.0) {
    dtemp = parms[1];
    ObitInfoListAlwaysPut (fitter->info, "GMajUp",   OBIT_double, dim, &dtemp);
    ObitInfoListAlwaysPut (fitter->info, "GMinUp",   OBIT_double, dim, &dtemp);
  }
  if (parms[2]>0.0) {
    dtemp = -parms[2];
    ObitInfoListAlwaysPut (fitter->info, "PosGuard", OBIT_double, dim, &dtemp);
  }
  if (parms[3]>0.0) {
    dtemp = image->myDesc->beamMaj / fabs (image->myDesc->cdelt[0]);
    ObitInfoListAlwaysPut (fitter->info, "GMajLow",  OBIT_double, dim, &dtemp);
    dtemp = image->myDesc->beamMin / fabs (image->myDesc->cdelt[0]);
    ObitInfoListAlwaysPut (fitter->info, "GMinLow",  OBIT_double, dim, &dtemp);
  }

  /* Clean beam in pixels */
  cbeam[0] = image->myDesc->beamMaj / fabs (image->myDesc->cdelt[0]);
  cbeam[1] = image->myDesc->beamMin / fabs (image->myDesc->cdelt[0]);
  cbeam[2] = image->myDesc->beamPA * DG2RAD;

  /* Fit */
  ier = ObitImageFitFit (fitter, image, reg, err);
  if (ier==1) iterLimit++;  /* Hit iteration limit? */
  
  /* Are residuals acceptable? If only 1 and not, then try 2 */
  tcut = gain*reg->models[0]->Peak;
  tcut = sqrt (icut*icut + tcut*tcut);
  if ((reg->nmodel==1) && (fabs(reg->RMSResid)>tcut)  && 
      (fabs(reg->peakResid)>tcut) && doMult) {
    breakIsland++;  /* count */
    /* Save old solution */
    oldReg =  ObitFitRegionCopy (reg, oldReg, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
    oldRMS  = reg->RMSResid;
    oldPeak = fabs(reg->peakResid);
    
    /* Get pixels data */
    blc[0] = reg->corner[0]; blc[1] = reg->corner[1];
    trc[0] = reg->corner[0] + reg->dim[0] - 1; 
    trc[1] = reg->corner[1] + reg->dim[1] - 1; 
    pixels = ObitFArraySubArr (image->image, blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
    
    /* reset region */
    fitTwo (pixels, reg, oldReg, cbeam, doPoint, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
    
    pixels = ObitFArrayUnref(pixels); /* Don't need further */
    
    /* Refit */
    ObitImageFitFit (fitter, image, reg, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
    
    /* Is this one any better? If not restore old 
       Use new one unless both peak and RMS worse */
    revert = (oldRMS < reg->RMSResid) && (oldPeak <  fabs(reg->peakResid));

    /* Must have significant difference in position or size */
    revert = revert || 
      ((fabs(reg->models[0]->DeltaX - oldReg->models[0]->DeltaX)<1.0) &&
       (fabs(reg->models[0]->DeltaY - oldReg->models[0]->DeltaY)<1.0) &&
       (fabs((reg->models[0]->parms[0]/oldReg->models[0]->parms[0])-1.0)<0.1) &&
       (fabs((reg->models[0]->parms[1]/oldReg->models[0]->parms[1])-1.0)<0.1));

    /* Both must be above threshold */
    revert = revert || 
      ((reg->models[0]->Peak<rcut) || (oldReg->models[0]->Peak<rcut));

    /* revert to previous model? */
    if (revert) {
      failBreak++; /* count */
      ObitFitRegionResize(reg, 1);
      reg->models[0] = ObitFitModelCopy (oldReg->models[0], reg->models[0], err);
      if (err->error) Obit_traceback_msg (err, routine, image->name);
    }
  } /* end retry */

  /* Reject (zero) components with peaks below threshold, count good */
  for (j=0; j<reg->nmodel; j++) {
    if (fabs(reg->models[j]->Peak)<rcut) {
      reg->models[j]->Peak = 0.0;
      rejectLoFlux++; /* count */
    } else nGood++;
  }
  
  /* Subtract fitted model from residual */
  ObitFitRegionSubtract (reg, image, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
  oldReg = ObitFitRegionUnref(oldReg);
  fitter = ObitImageFitUnref(fitter);
} /* end FitRegion */

/**
 * Convert FitRegionList to entries in an MF table on image
 * Sorts table
 * \param myInput Inputs object
 * \param regList List to convert 
 * \param image   ObitImage being described
 * \param err     Obit error/message stack object.
 * \return TableMF created on image
 */
/**
 * Set the beginning parameters for the fitting routine, using the  
 * estimates based on moment fits about the peak of the region.  
 * Routine translated from the AIPSish VSAD.FOR/SADDEF  
 * \param pixels  Contains data we're trying to fit 
 * \param doPoint True if not fitting size
 * \param doPA    True if fitting position angle
 * \param dmax    Flux density at peak
 * \param ixmax   X-coord of max in Array 
 * \param iymax   Y-coord of max in Array 
 * \param cb      Clean beam parameters (pix, pix, rad)
 * \param model   Model parameters, created
 */
void snglDef (ObitFArray *pixels, gboolean doPoint, gboolean doPA,
	      ofloat dmax, olong ixmax, olong iymax, 
	      ofloat cb[3], ObitFitModel **model) 
{
  olong   i, j, k, pts, sumpts, nx, ny;
  olong pos[2];
  ofloat sum2[3], sumd2[3], sum4[5], a, b, theta, slit, x, y, temp;
  ofloat g[6], *pixData;
  gboolean sing;
  ofloat cratio = 0.6, fblank = ObitMagicF();

  /* Pointer to data */
  pos[0] = pos[1] = 0;
  pixData = ObitFArrayIndex (pixels, pos);
  nx = pixels->naxis[0]; 
  ny = pixels->naxis[1]; 

  /* For point source use clean beam */
  if (doPoint) {
    g[3] = cb[0];
    g[4] = cb[1];
    g[5] = cb[2];

  } else {
    /* Loop around, find a lot of moments */
    for (i=0; i<3; i++) sum2[i]  = 0.0;
    for (i=0; i<3; i++) sumd2[i] = 0.0;
    for (i=0; i<5; i++) sum4[i]  = 0.0;
    slit = cratio * dmax;

    sumpts = 0;
    for (j=0; j<ny; j++) { /* loop 60 */
      for (i=0; i<nx; i++) { /* loop 50 */
	pts = j*ny + i;
	if ((pixData[pts] != fblank)  &&  (pixData[pts] >= slit)) {
	  x = i - ixmax;
	  y = j - iymax;
	  sumpts++;
	  for (k= 0; k<=2; k++) { /* loop 20 */
	    temp = 1.;
	    if (k > 0) temp *= powf (x, (ofloat)k);
	    if (k < 2) temp *= powf (y, (ofloat)(2-k));
	    sum2[k]  += temp;
	    sumd2[k] += pixData[pts]*temp;
	  } /* end loop  L20:  */

	  for (k= 0; k<=4; k++) { /* loop 30 */
	    temp = 1.0;
	    if (k > 0) temp *= powf (x, (ofloat)k);
	    if (k < 4) temp *= powf (y, (ofloat)(4-k));
	    sum4[k] += temp;
	  } /* end loop  L30:  */
	} 
      } /* end loop  L50:  */
    } /* end loop  L60:  */

    /* Convert moments into bmaj, bmin,  bpa */
    sing = sumpts <= 3;
    if (!sing) scdmom (dmax, sum2, sumd2, sum4, &a, &b, &theta, &sing);
    if (!sing) {
      g[3] = MAX (cb[0], MIN (4.0*cb[0], a));
      g[4] = MAX (cb[1], MIN (4.0*cb[1], b));
      g[5] = theta;
      /* Didn't work, use point spread  fn. */
    } else {
      g[3] = cb[0];
      g[4] = cb[1];
      g[5] = cb[2];
    } 
  } 
  /* Fill in Peak flux and position */
  g[0] = dmax;
  g[1] = ixmax;
  g[2] = iymax;

  if (!doPA) g[5] = 0.0;

  *model = ObitFitModelCreate("single", OBIT_FitModel_GaussMod, 
			      g[0], g[1], g[2], 3, &g[3]);
} /* end of routine snglDef */ 

/**
 * From the various 2nd moments find the best least-squares quadratic  
 * fit to the values near the peak.  The assumed form of the fit is  
 * I = dmax - a*x*x - b*x*y - c*y*y  
 * If there are too few points to make a fit, or if there are other  
 * problems, SINGUL will be set to be .TRUE.  
 * Routine translated from the AIPSish VSAD.FOR/SCDMOM  
 * \param dmax     Maximum in array 
 * \param sum2     Vector of moments, SUM(I) is the sum over 
 *                 valid pixels of X**I * Y **(2-I) 
 * \param sumd2    Sum of Flux * X**I * Y**(2-I) 
 * \param sum4     Sum of X**I * Y**(4-I) 
 * \param a        [out] Estimate of major axis 
 * \param b        [out] Estimate of minor axis 
 * \param theta    [out] Estimate of position angle (radians) 
 * \param singul   [out] If .FALSE. couldn't find a decent solution 
 */
void scdmom (ofloat dmax, ofloat* sum2, ofloat* sumd2, ofloat* sum4, 
	     ofloat* a, ofloat* b, ofloat* theta, gboolean* singul) 
{
  ofloat mat[4][3], sol[3];
  ofloat bmin, bplus, x,temp, w, denom;
  olong   i, j, k, pivot[3];
  
  /* Set up matrix for least square  solution */
  for (i= 0; i<=2; i++) { /* loop 10 */
    mat[2][i] = sumd2[i] - dmax * sum2[i];
    for (j= 0; j<=2; j++) { /* loop 5 */
      mat[j][i] = sum4[i-j+1];
    } /* end loop   L5:  */
  } /* end loop  L10:  */

  /* Pivoted Gaussian elimination */
  *singul = FALSE;

  /* Reduce to Right triangular */
  for (i= 0; i<=2; i++) { /* loop 40 */
    x = 0.0;
    /* Find pivot */
    for (j= 0; j<=2; j++) { /* loop 20 */
      if (fabs(mat[j][i]) > x) {
	x = fabs (mat[j][i]);
	pivot[i] = j;
      } 
    } /* end loop  L20:  */

    *singul = (x == 0.0);
    if (!singul) {
      for (j= i+1; j<=2; j++) { /* loop 30 */
	x = mat[pivot[i]][j] /  mat[pivot[i]][i];
	for (k= 0; k<=3; k++) { /* loop 25 */
	  temp = mat[k][j];
	  mat[k][j] = temp - mat[k][i] * x;
	  
	  /* Anything that eliminates too  well is zero */
	  if ((k<3)  &&  (fabs(mat[k][j]) < 1.0e-4*fabs (temp))) mat[k][j] = 0.;
	} /* end loop  L25:  */
      } /* end loop  L30:  */

    } else {  /* No can do */
      return;
    } 
  } /* end loop  L40:  */

  /* From Right triangle find  solutions */
  for (i= 2; i<=0; i++) { /* loop 60 */
    for (j= 2; j<=i+1; j++) { /* loop 50 */
      mat[2][i] -= mat[2][j] * mat[pivot[j]][i];
    } /* end loop  L50:  */
    mat[2][i] /= mat[pivot[i]][i];
  } /* end loop  L60:  */
  
  /* Unpivot */
  for (i= 0; i<=2; i++) { /* loop 70 */
    sol[pivot[i]] = -mat[2][i] / dmax;
  } /* end loop  L70:  */

    /*  Does solution for a,b,c make sense? */
    if ((sol[0] < 0.)  ||  (sol[2] < 0.)  ||  (sol[1]*sol[1] >= 4.*sol[0]*sol[2])) {
      *singul = TRUE;
      return;
    } 

    /*  Convert to bmaj, bmin, bpa */
      bmin  = sol[2] - sol[0];
      bplus = sol[2] + sol[0];
      w = sqrt (bmin*bmin + sol[1]*sol[1]);

      /* empirical fudge factor */
      x = 1.6 * logf(16.0);
      denom = (bplus + w);
      if (fabs (denom) <= 1.0e-20) denom = 1.0;
      *b = sqrt (x/denom);
      denom = (bplus - w);
      if (fabs (denom) <= 1.0e-20) denom = 1.0;
      *a = sqrt (x/denom);
      *theta = atan2 (sol[1], -bmin);
  } /* end of routine scdmom */ 

/**
 * There are peakNo > 1 peaks above the cutoff.  Set the initial  guesses to be 
 * points at each peak.  
 * Routine translated from the AIPSish VSAD.FOR/MULDEF  
 * \param peakNo  Number of peaks 
 * \param doPA    True if fitting position angle
 * \param xpk     X-coord of peaks relative to window 
 * \param ypk     Y-coord of peaks relative to window 
 * \param spk     Flux at each peak 
 * \param cb      Clean beam parameters (pix, pix, rad)
 * \param model   Array of model parameters, entries created
 */
void multDef (olong peakNo, gboolean doPA, olong* xpk, olong* ypk, ofloat* spk, ofloat cb[3], 
	      ObitFitModel *models[]) 
{
  olong i;
  ofloat g[6];
  gchar label[20];

  for (i=0; i<peakNo; i++) { /* loop 100 */
    g[0] = spk[i];
    g[1] = xpk[i];
    g[2] = ypk[i];
    g[3] = cb[0];
    g[4] = cb[1];
    g[5] = cb[2];
    if (!doPA) g[5] = 0.0;
    sprintf (label, "comp%d", i);
    models[i] = ObitFitModelCreate(label, OBIT_FitModel_GaussMod, 
				  g[0], g[1], g[2], 3, &g[3]);
  } /* end loop  L100: */
} /* end of routine multDef */ 

/**
 * fitTwo prepares for redoing a fit trying 2 overlapping Gaussians  
 * rather than a single one.  
 * Routine translated from the AIPSish VSAD.FOR/SADRED  
 * \param pixels  Region pixel subarray
 * \param reg     Region to update
 * \param oldReg  Previous region fit
 * \param cb      Clean beam (pix, pix, rad)
 * \param doPoint True if not fitting size
 * \param err     Obit error/message stack object.
 */
void fitTwo (ObitFArray *pixels, ObitFitRegion *reg, ObitFitRegion *oldReg, 
	     ofloat cb[3], gboolean doPoint, ObitErr* err) 
{
  olong   i, nx, ny, ngauss;
  olong pos[2];
  ofloat datmax, *pixData, g[3][6];
  /*gchar *routine = "fitTwo";*/

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */

  /* Pointer to data */
  pos[0] = pos[1] = 0;
  pixData = ObitFArrayIndex (pixels, pos);
  nx = pixels->naxis[0]; 
  ny = pixels->naxis[1]; 

  /* Resize output */
  ngauss = 2;
  ObitFitRegionResize(reg, ngauss);

  /* Old solution */
  g[2][0] = oldReg->models[0]->Peak;
  g[2][1] = oldReg->models[0]->DeltaX;
  g[2][2] = oldReg->models[0]->DeltaY;
  g[2][3] = oldReg->models[0]->parms[0];
  g[2][4] = oldReg->models[0]->parms[1];
  g[2][5] = oldReg->models[0]->parms[2];
 
  /* Find maximum */
  datmax = ObitFArrayMax (pixels, pos);

  /* 1st = CB at max */
  g[0][0] = 0.8 * datmax;
  g[0][1] = pos[0];
  g[0][2] = pos[1];
  g[0][3] = cb[0];
  g[0][4] = cb[1];
  g[0][5] = cb[2];
  g[1][3] = cb[0];
  g[1][4] = cb[1];
  g[1][5] = cb[2];

  /* separated */
  if (fabs(g[0][1]- g[2][1]) + fabs (g[0][2]-g[2][2]) > 1.5) {
    g[1][0] = 0.8 * g[2][0];
    g[1][1] = 2.  * g[2][1] - g[0][1];
    g[1][2] = 2.  * g[2][2] - g[0][2];
    if (!doPoint) {
      g[1][3] = (cb[0] + g[2][3]) / 2.0;
      g[1][4] = (cb[1] + g[2][4]) / 2.0;
      g[1][5] = (cb[2] + g[2][5]) / 2.0;
    } 

    /* core-halo ? */
  } else {
    g[0][0] = 0.8 * g[0][0];
    g[1][0] = datmax - g[0][0];
    g[1][1] = 2. * g[2][1] - g[0][1];
    g[1][2] = 2. * g[2][2] - g[0][2];
    if (!doPoint) {
      g[1][3] = 2 * cb[0];
      g[1][4] = 2 * cb[1];
      g[1][5] = cb[2];
    } 
  }
  /* Save to res */
  for (i=0; i<2; i++) {
    reg->models[i] = 
      ObitFitModelCreate ("model", OBIT_FitModel_GaussMod,
			  g[i][0], g[i][1], g[i][2], 3, &g[i][3]);
			  
  }

} /* end of routine fitTwo */ 

/*------------  IslandElem Functions         ----------*/
/**
 * IslandElem Constructor
 * \param index Island number
 * \param peak  Peak value in Island
 * \param blc   BLC corner of island (0-rel)
 * \param trc   TRC corner of island (0-rel)
 * \return the new  object.
 */
static IslandElem* newIslandElem (olong index, ofloat peak, 
				  olong blc[2], olong trc[2])
{
  IslandElem *out=NULL;

  out = g_malloc0(sizeof(IslandElem));
  out->index   = index;  
  out->peak    = peak;  
  out->blc[0]  = blc[0];  
  out->blc[1]  = blc[1];  
  out->trc[0]  = trc[0];  
  out->trc[1]  = trc[1];  

  return out;
} /* end newIslandElem */

/**
 * IslandElem Destructor 
 * \param in Object to delete
 */
static void freeIslandElem (IslandElem *in)
{
  if (in) g_free(in);
} /* end freeIslandElem */

/**
 * Print IslandElem 
 * \param in    Object to print
 * \param file  FILE* to write to
 */
static void IslandElemPrint (IslandElem *elem, FILE *file)
{
  fprintf (file, "%d peak=%f blc= %d %d trc=%d% d\n", 
	   elem->index, elem->peak, elem->blc[0], elem->blc[1], 
	   elem->trc[0], elem->trc[1]);
} /* end IslandElemPrint */

/*------------  IslandList Functions         ----------*/
/**
 * IslandList Constructor
 * \return the new  object.
 */
static IslandList* newIslandList (void)
{
  IslandList *out=NULL;

  out = g_malloc0(sizeof(IslandList));
  out->list    = NULL;
  out->number  = 0;
  out->maxIndex = 0;

  return out;
} /* end newIslandList */

/**
 * IslandList Destructor 
 * \param in Object to delete
 */
static void freeIslandList (IslandList *in)
{
  if (!in) return;
  GSList *tmp;

  /* loop through list deleting elements */
  tmp = in->list;
  while (tmp!=NULL) {
    if (tmp->data) freeIslandElem(tmp->data);
    tmp = g_slist_next(tmp);
  }

  /* delete members  */
  g_slist_free(in->list);

  /* delete object */
  g_free (in);
    
} /* end freeIslandElem */

/**
 * Print IslandList 
 * \param in    Object to print
 * \param file  FILE* to write to
 */
static void IslandListPrint (IslandList *in, FILE *file)
{
  GSList *tmp;
  IslandElem *elem;

  fprintf (file, "Listing of IslandList\n");
  /* loop through list printing elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (IslandElem*)tmp->data;
    IslandElemPrint(elem, file);
    tmp = g_slist_next(tmp);
  }
} /* end IslandElemPrint */

/**
 * Append elem to list in
 * \param in   Object with table to add elem to
 * \param elem the element to add. MUST NOT have an index already in list
 */
static void IslandListAppend (IslandList *in, IslandElem *elem)
{
  IslandElem *tmp = NULL;

  /* Make sure it's not already in list */
  tmp =  IslandListFind (in, elem->index);
  if (tmp!=NULL) { /* trouble - die in a noisy fashion */
    g_error ("IslandAppend: trying to add redundant island: %d new %d",
	     tmp->index, elem->index);
  }

  /* append to list */
  in->list = g_slist_append  (in->list, elem);
  in->number++;
  in->maxIndex = MAX (in->maxIndex, elem->index);

} /* end IslandListAppend */

/**
 * Remove elem from list in
 * \param in   Object with table to remove elem from
 * \param elem the element to remove. MUST be in list, freeed
 */
static void IslandListRemove (IslandList *in, IslandElem *elem)
{
  /* remove from table */
  in->list = g_slist_remove (in->list, elem);

  in->number--; /* keep count */
  if (elem->index==in->maxIndex) in->maxIndex--;
  freeIslandElem (elem);

} /* end IslandListRemove  */

/**
 * Find a pointer is in list in
 * \param in    Object with table to search
 * \param index island number of item to search for
 * \return pointer to element containing item, NULL if not found.
 */
static IslandElem* IslandListFind (IslandList *in, olong index)
{
  GSList *tmp;
  IslandElem *elem;

  /* loop through list testing elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (IslandElem*)tmp->data;
    /* check if this is a match */
    if (elem->index==index) return elem;
    tmp = g_slist_next(tmp);
  }
  return NULL; /* didn't find */
} /* end IslandListFind */

/**
 * Find element with highest peak
 * \param in    Object with table to search
 * \return pointer to element containing item, NULL if not found.
 */
static IslandElem* IslandListFindHi (IslandList *in)
{
  GSList *tmp;
  ofloat max = -1.0e20;
  IslandElem *elem, *theOne=NULL;

  if (!in->list) return NULL;  /* empty? */

  /* loop through list testing elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (IslandElem*)tmp->data;
    /* check if this peak is higher */
    if (elem->peak>max) {
      max = elem->peak;
      theOne = elem;
    }
    tmp = g_slist_next(tmp);
  }
  return theOne;
} /* end IslandListFindHi */

/**
 * Merge two IslandElements and update island lists
 * \param in    Island list with elements, in1, in2
 * \param in1   index of first element
 * \param in2   index of second element
 * \param nx    length of prev, cur
 * \param prev  previous island numbers on image row
 * \param curr  current island numbers on image row
 */
static void IslandListMerge (IslandList *in, olong in1, olong in2, 
			     olong nx, olong *prev, olong *curr)
{
  IslandElem *elem1, *elem2, *elemLo, *elemHi;
  olong i, iLo, iHi;

  /* Dont bother if they're the same */
  if (in1==in2) return;

  /* Look them up */
  elem1 = IslandListFind(in, in1);
  elem2 = IslandListFind(in, in2);
  iLo = MIN (elem1->index, elem2->index);
  iHi = MAX (elem1->index, elem2->index);

  /* Which one is which? */
  if (elem1->index == iLo) {elemLo = elem1; elemHi = elem2;}
  else {elemLo = elem2; elemHi = elem1;}
  /* Merge blc, trc to elemLo */
  elemLo->blc[0] = MIN (elemLo->blc[0], elemHi->blc[0]);
  elemLo->blc[1] = MIN (elemLo->blc[1], elemHi->blc[1]);
  elemLo->trc[0] = MAX (elemLo->trc[0], elemHi->trc[0]);
  elemLo->trc[1] = MAX (elemLo->trc[1], elemHi->trc[1]);
  /* total peak */
  elemLo->peak = MAX (elemLo->peak, elemHi->peak);

  /* delete elemHi from list */
  IslandListRemove (in, elemHi);

  /* Update island numbers */
  for (i=0; i<nx; i++) {
    if (prev[i]==iHi) prev[i] = iLo;
    if (curr[i]==iHi) curr[i] = iLo;
  }
} /* end IslandListMerge */

/**
 * expand island in1 to include pixwl (i,j)
 * \param in    Island list with element in1
 * \param in1   index of element
 * \param val   data value of new element
 * \param i     x pixel
 * \param j     y pixel
 */
static void IslandListExpand (IslandList *in, olong in1, 
			      ofloat val, olong i,  olong j)
{
  IslandElem *elem1;

  /* Look it up */
  elem1 = IslandListFind(in, in1);

  /* expand */
  elem1->blc[0] = MIN (elem1->blc[0], i);
  elem1->blc[1] = MIN (elem1->blc[1], j);
  elem1->trc[0] = MAX (elem1->trc[0], i);
  elem1->trc[1] = MAX (elem1->trc[1], j);

  elem1->peak = MAX (elem1->peak, val); /* new max? */ 
} /* end IslandListExpand */

