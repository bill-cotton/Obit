/* $Id$  */
/* Obit task - interpolates image to different geometry               */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
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
#include "ObitImageUtil.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* HGeomIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void HGeomOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Get input image */
ObitImage* getInputImage (ObitInfoList *myInput, ObitErr *err);
/* Get template image */
ObitImage* getTemplImage (ObitInfoList *myInput, ObitErr *err);
/* Define output image */
ObitImage* getOutputImage (ObitInfoList *myInput, ObitImage *tmplImage, 
			   ObitImage *inImage, ObitErr *err);
/* Interpolate */
void HGeomInterp (ObitInfoList* myInput, ObitImage* inImage, 
		  ObitImage* outImage, ObitErr* err);
/* Write history */
void HGeomHistory (ObitInfoList* myInput, ObitImage* inImage, 
		   ObitImage* outImage, ObitErr* err);

/* Program globals */
gchar *pgmName = "HGeom";       /* Program name */
gchar *infile  = "HGeom.inp";   /* File with program inputs */
gchar *outfile = "HGeom.out";   /* File to contain program outputs */
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
/*   Obit program - extract a subregion of an image                       */
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitImage    *inImage = NULL, *tmplImage = NULL, *outImage = NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = HGeomIn (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) goto exit;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Get input Image Object */
  inImage = getInputImage (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Get template image object */
  tmplImage = getTemplImage (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Define output Image Object */
  outImage = getOutputImage (myInput, tmplImage, inImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Copy */
  HGeomInterp (myInput, inImage, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Do history */
  HGeomHistory ( myInput, inImage, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inImage   = ObitUnref(inImage);
  tmplImage = ObitUnref(tmplImage);
  outImage  = ObitUnref(outImage);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* HGeomIn (int argc, char **argv, ObitErr *err)
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
  ObitInfoList* list;
  gchar *routine = "HGeomIn";

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
      
     } else if (strcmp(arg, "-DataType") == 0) { /* Image type AIPS or FITS */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "DataType", OBIT_string, dim, strTemp);
      
    } else if (strcmp(arg, "-inSeq") == 0) { /* input AIPS sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "inSeq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-inDisk") == 0) { /* input disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "inDisk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-inName") == 0) { /* input AIPS Name*/
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inName", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-inClass") == 0) { /* input AIPS Class*/
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inClass", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-inFile") == 0) { /*input  FITSFile */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inFile", OBIT_string, dim, strTemp);
      
    } else if (strcmp(arg, "-in2Seq") == 0) { /* Template  AIPS sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "in2Seq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-in2Disk") == 0) { /* Template disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "in2Disk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-in2Name") == 0) { /* Template  AIPS Name*/
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "in2Name", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-in2Class") == 0) { /* Template  AIPS Class*/
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "in2Class", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-in2File") == 0) { /* Template  FITSFile */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "in2File", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-outFile") == 0) { /* output  FITSFile */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outFile", OBIT_string, dim, strTemp);
      
      } else if (strcmp(arg, "-outName") == 0) { /* output AIPS Name*/
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outName", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-outClass") == 0) { /* output AIPS Class*/
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outClass", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outSeq") == 0) { /* output AIPS sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outSeq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-outDisk") == 0) { /* output disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outDisk", OBIT_oint, dim, &itemp, err);
      
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
} /* end HGeomIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: HGeom -input file -output ofile [args]\n");
    fprintf(stderr, "Interpolate image to grid defined by template image\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def HGeom.inp\n");
    fprintf(stderr, "  -output output result file, def HGeom.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType 'AIPS' or 'FITS' type for input image\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inFile input FITS Image file\n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -in2File template FITS Image file\n");
    fprintf(stderr, "  -in2Name template  AIPS file name\n");
    fprintf(stderr, "  -in2Class template AIPS file class\n");
    fprintf(stderr, "  -in2Seq itemplate AIPS file sequence\n");
    fprintf(stderr, "  -in2Disk template (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -outFile output FITS Image file\n");
    fprintf(stderr, "  -outName output AIPS file name [def inName]\n");
    fprintf(stderr, "  -outClass output AIPS file class [def 'SUBIM']\n");
    fprintf(stderr, "  -outSeq output AIPS file sequence\n");
    fprintf(stderr, "  -outDisk output (AIPS or FITS) disk number (1-rel) \n");
    
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
/*     AIPSdirs  Str [?,?]  AIPS directories [def {"AIPSdata/"}]          */
/*     DataType  Str [4]    "AIPS" or "FITS" [def {"FITS"}]               */
/*     inFile    Str [?]    input FITS image file name [def "Image.fits"] */
/*     inName    Str [12]   input AIPS image name  [no def]               */
/*     inClass   Str [6]    input AIPS image class  [no def]              */
/*     inSeq     Int        input AIPS image sequence no  [no def]        */
/*     inDisk    Int        input AIPS or FITS image disk no  [def 1]     */
/*     in2File   Str [?]    tmpl FITS image file name [def "Image.fits"]  */
/*     in2Name   Str [12]   tmpl AIPS image name  [no def]                */
/*     in2Class  Str [6]    tmpl AIPS image class  [no def]               */
/*     in2Seq    Int        tmpl AIPS image sequence no  [no def]         */
/*     in2Disk   Int        tmpl AIPS or FITS image disk no  [def 1]      */
/*     BLC       Int  [7]   bottom-left (1-rel) corner[def {1,1,...)]     */
/*     TRC       Int  [7]   top-right (1-rel) corner[def {0,0,...)]       */
/*     size      Int  [2]   optional output size                          */
/*     crpix     Flt  [2]   optional output reference pixel (2 axes)      */
/*     hwidth    Int        interpolation kernal halfwidth                */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  ofloat farray[2];
  oint   itemp, iarray[2];
  olong   blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong   trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ObitInfoList *out = newObitInfoList();
  gchar *routine = "defaultInputs";

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
  strTemp = "Image.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "HGeomName";
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

  /* BLC, TRC */
  dim[0] = IM_MAXDIM;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "BLC", OBIT_oint, dim, blc, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  ObitInfoListPut (out, "TRC", OBIT_oint, dim, trc, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Template FITS file name */
  strTemp = "Template.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "in2File", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Template AIPS file name */
  strTemp = "Template";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "in2Name", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Template AIPS file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "in2Class", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Template AIPS sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "in2Seq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Template AIPS or FITS disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "in2Disk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output FITS file name */
  strTemp = "HGeomOut.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS file name */
  strTemp = "HGeomName";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "outSeq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* OutputAIPS or FITS disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "outDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* size */
  dim[0] = 2;
  iarray[0] = iarray[1] = 0;
  ObitInfoListPut (out, "size", OBIT_oint, dim, iarray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* crpix */
  dim[0] = 2;
  farray[0] = farray[1] = 0.0;
  ObitInfoListPut (out, "crpix", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* hwidth */
  dim[0] = 1;
  itemp = 2;
  ObitInfoListPut (out, "hwidth", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  return out;
} /* end defaultInputs */

/*----------------------------------------------------------------------- */
/*  Create default output ObitInfoList                                    */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*  Values:                                                               */
/*     Mean     Int        Image pixel mean  [0.0]                        */
/*     RMS      Int        Image pixel rms   [0.0]                        */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultOutputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat ftemp;
  ObitInfoList *out = newObitInfoList();
  gchar *routine = "defaultOutputs";

  /* add parser items */
  /* Image mean */
  dim[0] = 1; dim[1] = 1;
  ftemp = 0.0;
  ObitInfoListPut (out, "Mean", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Image mean */
  dim[0] = 1; dim[1] = 1;
  ftemp = 0.0;
  ObitInfoListPut (out, "RMS", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  return out;
} /* end defaultOutputs */

/*----------------------------------------------------------------------- */
/*  Get input image                                                       */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitImage with input image                                       */
/*----------------------------------------------------------------------- */
ObitImage* getInputImage (ObitInfoList *myInput, ObitErr *err)
{
  ObitImage    *inImage = NULL;
  ObitInfoType type;
  olong         i, Aseq, disk, cno;
  gchar        *Type, *strTemp, inFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *routine = "getInputImage";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inImage;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input UV data Object */
  inImage = newObitImage("input Image");
  
  /* Get region from myInput */
  ObitInfoListGetTest(myInput, "BLC", &type, dim, blc); /* BLC */
  ObitInfoListGetTest(myInput, "TRC", &type, dim, trc); /* TRC */

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
      if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "inSeq", OBIT_oint, dim, &Aseq);
    } 


    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);
    
    /* define object */
    ObitImageSetAIPS (inImage, OBIT_IO_byPlane, disk, cno, AIPSuser, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);
    
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
    ObitImageSetFITS (inImage, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return inImage;
  }

  /* Ensure inImage fully instantiated and OK */
  ObitImageFullInstantiate (inImage, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);

  /* Set defaults BLC, TRC */
  for (i=0; i<IM_MAXDIM; i++) {
    if (blc[i]<=0) blc[i] = 1;
    blc[i] = MAX (1,  blc[i]);
    if (trc[i]<=0) trc[i] = inImage->myDesc->inaxes[i];
    trc[i] = MIN (trc[i], inImage->myDesc->inaxes[i]);
  }

  /* Save blc, trc */
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (myInput, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (inImage->info, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (myInput, "TRC", OBIT_long, dim, trc);
  ObitInfoListAlwaysPut (inImage->info, "TRC", OBIT_long, dim, trc);

  return inImage;
} /* end getInputImage */

/*----------------------------------------------------------------------- */
/*  Get template image                                                    */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitImage with template image                                    */
/*----------------------------------------------------------------------- */
ObitImage* getTemplImage (ObitInfoList *myInput, ObitErr *err)
{
  ObitImage    *tmplImage = NULL;
  ObitInfoType type;
  olong         Aseq, disk, cno;
  gchar        *Type, *strTemp, inFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *routine = "getTemplImage";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return tmplImage;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input UV data Object */
  tmplImage = newObitImage("input Image");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
    /* input AIPS disk */
    ObitInfoListGet(myInput, "in2Disk", &type, dim, &disk, err);
    /* input AIPS name */
    if (ObitInfoListGetP(myInput, "in2Name", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aname, strTemp, 13);
    } else { /* Didn't find */
      strncpy (Aname, "No Name ", 13);
    } 
    Aname[12] = 0;
    /* input AIPS class */
    if  (ObitInfoListGetP(myInput, "in2Class", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    Aclass[6] = 0;
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "in2Seq", &type, dim, &Aseq, err);

    /* if ASeq==0 want highest existing sequence */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", tmplImage);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "in2Seq", OBIT_oint, dim, &Aseq);
    } 


    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", tmplImage);
    
    /* define object */
    ObitImageSetAIPS (tmplImage, OBIT_IO_byPlane, disk, cno, AIPSuser, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", tmplImage);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    if (ObitInfoListGetP(myInput, "in2File", &type, dim, (gpointer)&strTemp)) {
      strncpy (inFile, strTemp, 128);
    } else { 
      strncpy (inFile, "No_Filename_Given", 128);
    }
    
    /* input FITS disk */
    ObitInfoListGet(myInput, "in2Disk", &type, dim, &disk, err);

    /* define object */
    ObitImageSetFITS (tmplImage, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", tmplImage);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return tmplImage;
  }

  /* Ensure tmplImage fully instantiated and OK */
  ObitImageFullInstantiate (tmplImage, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", tmplImage);

  return tmplImage;
} /* end getTemplImage */

/*----------------------------------------------------------------------- */
/*  Define output image, cloned from template with size override          */
/*  Image is created with program parameters applied                      */
/*  Output image will have geometry defined by the first two axes of the  */
/*  Template image, possibly modified by parameters size and crpix but    */
/*  the geometry on the other axes of the input image.                    */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      tmplImage Template image to copy                                  */
/*      inImage   Image to interpolate                                    */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitImage for output image                                       */
/*----------------------------------------------------------------------- */
ObitImage* getOutputImage (ObitInfoList *myInput, ObitImage *tmplImage, 
			   ObitImage *inImage, ObitErr *err)
{
  ObitImage    *outImage = NULL;
  ObitInfoType type;
  gboolean     exist;
  olong         Aseq, disk, idisk, cno;
  gchar        *Type, *strTemp, *strTemp2, outFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong         i, j, size[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ofloat       crpix[2];
  gchar *routine = "getOutputImage";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outImage;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input UV data Object */
  outImage = newObitImage("output Image");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
    /* output AIPS disk default = inDisk*/
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    ObitInfoListGetTest(myInput, "outDisk", &type, dim, &disk);
    /* output AIPS name - default = input */
    ObitInfoListGetP(myInput, "inName", &type, dim, (gpointer)&strTemp2);
    if (ObitInfoListGetP(myInput, "outName", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aname, strTemp, 13);
    } else { /* Didn't find - use inName */
      strncpy (Aname, strTemp, 13);
    } 
    if (!strncmp (Aname, "    ", 4)) strncpy (Aname, strTemp2, 13);
    Aname[12] = 0;
    /* output AIPS class */
    if  (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find - use "SUBIM"*/
      strncpy (Aclass, "SUBIM ", 7);
    }
    /* Default if output class blank */
    if ((Aclass[0]==' ') && (Aclass[1]==' ') && (Aclass[2]==' '))
      strncpy (Aclass, "SUBIM ", 7);
    Aclass[6] = 0;
    /* output AIPS sequence */
    Aseq = 1;
    ObitInfoListGetTest(myInput, "outSeq", &type, dim, &Aseq);

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);

      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Allocate catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
     /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Making AIPS image %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);
   
   /* define object */
    ObitImageSetAIPS (outImage, OBIT_IO_byPlane, disk, cno, AIPSuser, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */
    /* output FITS file name */
    ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp2);
    if (ObitInfoListGetP(myInput, "outFile", &type, dim, (gpointer)&strTemp)) {
      strncpy (outFile, strTemp, 128);
    } else { 
      g_snprintf (outFile, 129, "HGeom%s", strTemp2);
    }
    /* If blank use HGeom+inFile */
    if (!strncmp (outFile, "    ", 4)) {
      g_snprintf (outFile, 129, "HGeom%s", strTemp2);
    }
    
    /* output FITS disk default = inDisk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    idisk = disk;
    ObitInfoListGetTest(myInput, "outDisk", &type, dim, &disk);
    if (disk<=0) disk = idisk;

    /* Tell about it */
    ObitTrimTrail(outFile);
    Obit_log_error(err, OBIT_InfoErr, "Making FITS image %s on disk %d",
		   outFile, disk);

    /* define object */
    ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, outFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return outImage;
  }

  /* Clone from template image */
  ObitInfoListGetTest(myInput, "size", &type, dim, size);
  /* get higher planes from input Image */
  for (i=2; i<IM_MAXDIM; i++) size[i] =  MAX (1, inImage->myDesc->inaxes[i]);
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut(outImage->info, "Size", type, dim, size);  
  ObitImageClone (tmplImage, outImage, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);

  /* open image to update header */
  ObitImageOpen (outImage, OBIT_IO_ReadWrite, err);

  /* Copy descriptive information from the input image for dimensions 
     higher than 2 */
  outImage->myDesc->naxis = inImage->myDesc->naxis;
  for (i=2; i<inImage->myDesc->naxis; i++) {
    for (j=0; j<IMLEN_KEYWORD; j++) 
      outImage->myDesc->ctype[i][j] = inImage->myDesc->ctype[i][j];
    outImage->myDesc->cdelt[i] = inImage->myDesc->cdelt[i];
    outImage->myDesc->crpix[i] = inImage->myDesc->crpix[i];
    outImage->myDesc->crota[i] = inImage->myDesc->crota[i];
  }

  /* Adjust reference pixel if needed */
  crpix[0] = crpix[1] = 0.0;
  ObitInfoListGetTest(myInput, "crpix", &type, dim, crpix);
  if ((crpix[0]!=0.0) || (crpix[1]!=0.0) || 
      (size[0]!=0.0) || (size[1]!=0.0)) {
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    /* if crpix specified use it */
    if ((crpix[0]!=0.0) || (crpix[1]!=0.0)) {
      /* just use crpix as is */
    } else { /* default for resize - use center of output */
      crpix[0] = 1.0 + outImage->myDesc->inaxes[0] * 0.5;
      crpix[1] = 1.0 + outImage->myDesc->inaxes[1] * 0.5;
    }
      outImage->myDesc->crpix[0] = crpix[0];
      outImage->myDesc->crpix[1] = crpix[1];
  } /* end ajust crpix */

  outImage->myStatus = OBIT_Modified;  /* Force header update */

  ObitImageClose (outImage, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
  
  return outImage;
} /* end getOutputImage */

/*----------------------------------------------------------------------- */
/*  Write History for HGeom                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inImage   Image to copy history from                              */
/*      outImage  Image to write history to                               */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void HGeomHistory (ObitInfoList* myInput, ObitImage* inImage, 
		    ObitImage* outImage, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "in2File",  "in2Disk", "in2Name", "in2Class", "in2Seq",
    "BLC",  "TRC",  "size", "crpix", "hwidth", 
    NULL};
  gchar *routine = "HGeomHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  /* Do history  */
  inHistory  = newObitDataHistory ((ObitData*)inImage,  OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory ((ObitData*)outImage, OBIT_IO_WriteOnly, err);

  /* If FITS copy header */
  if (inHistory->FileType==OBIT_IO_FITS) {
    ObitHistoryCopyHeader (inHistory, outHistory, err);
  } else { /* simply copy history */
     ObitHistoryCopy (inHistory, outHistory, err);
  }
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

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end HGeomHistory  */

/*----------------------------------------------------------------------- */
/*  Interpolate image and copy tables                                     */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inImage   Image to copy history from, shold be fully instantiated */
/*      outImage  Image to write history to                               */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void HGeomInterp (ObitInfoList* myInput, ObitImage* inImage, 
		  ObitImage* outImage, ObitErr* err)
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         inPlane[5], outPlane[5], ip1, ip2, ip3, ip4, ip5;
  olong         i, itemp, naxis[IM_MAXDIM];
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar        *exCClist[] = {"AIPS HI", "AIPS PL", "AIPS SL",
			      NULL};
  olong        hwidth;
  gchar        *routine = "HGeomInterp";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  /* Control parameters */
  itemp = 2;
  ObitInfoListGetTest(myInput, "hwidth", &type, dim, &itemp);
  hwidth = itemp; 

  /* Get region from myInput */
  ObitInfoListGetTest(myInput, "BLC", &type, dim, blc); /* BLC */
  ObitInfoListGetTest(myInput, "TRC", &type, dim, trc); /* TRC */

  /* Array to loop over */
  for (i=2; i<IM_MAXDIM; i++) naxis[i-2] = MAX (1, inImage->myDesc->inaxes[i]);

  /* Loop oer all possible planed */
  for (ip5=1; ip5<=naxis[4]; ip5++) {
    inPlane[4] = outPlane[4] = ip5+ MAX (blc[6],1) - 1;
    for (ip4=1; ip4<=naxis[3]; ip4++) {
      inPlane[3] = outPlane[3] = ip4+ MAX (blc[5],1) - 1;
      for (ip3=1; ip3<=naxis[2]; ip3++) {
	inPlane[2] = outPlane[2] = ip3+ MAX (blc[4],1) - 1;
	for (ip2=1; ip2<=naxis[1]; ip2++) {
	  inPlane[1] = outPlane[1] = ip2+ MAX (blc[3],1) - 1;
	  for (ip1=1; ip1<=naxis[0]; ip1++) {
	    inPlane[0] = outPlane[0] = ip1 + MAX (blc[2],1) - 1;
	    /* Interpolate plane */
	    ObitImageUtilInterpolateImage (inImage, outImage, inPlane, outPlane,
					   hwidth, err);
	    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
	    
	  } /* end loop over axis 1 */
	} /* end loop over axis 2 */
      } /* end loop over axis 3 */
    } /* end loop over axis 4 */
  } /* end loop over axis 5 */

  /* Copy  tables */
  ObitDataCopyTables ((ObitData*)inImage, (ObitData*)outImage, exCClist, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);


} /* end HGeomInterp  */


