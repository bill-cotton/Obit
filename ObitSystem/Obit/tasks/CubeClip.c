/* $Id$  */
/* Obit task - Clip insignificant pixels                              */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitImage.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitTableCCUtil.h"
#include "ObitFArrayUtil.h"
#include "ObitFeatherUtil.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* CubeClipIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void CubeClipOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Get input image */
ObitImage* getInputImage (ObitInfoList *myInput, ObitErr *err);
/* Define output image */
ObitImage* getOutputImage (ObitInfoList *myInput, ObitErr *err);
/* Clip */
void CubeClipClip (ObitInfoList* myInput, ObitImage* inImage, 
		   ObitImage* outImage, ObitErr* err);
/* Write history */
void CubeClipHistory (ObitInfoList* myInput, ObitImage* inImage, 
		      ObitImage* outImage, ObitErr* err);

/* Program globals */
gchar *pgmName = "CubeClip";       /* Program name */
gchar *infile  = "CubeClip.inp";   /* File with program inputs */
gchar *outfile = "CubeClip.out";   /* File to contain program outputs */
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
/*   Obit program - Clip a subregion of an image                          */
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitImage    *inImage = NULL, *outImage = NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = CubeClipIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) goto done;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Get input Image Object */
  inImage = getInputImage (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto done;

  /* Define output Image Object */
  outImage = getOutputImage (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto done;

  /* Clip */
  CubeClipClip (myInput, inImage, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto done;

  /* Do history */
  CubeClipHistory ( myInput, inImage, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto done;

 done:
  /* show any messages and errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inImage   = ObitUnref(inImage);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* CubeClipIn (int argc, char **argv, ObitErr *err)
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
  oint    itemp, i, j, k, iarr[IM_MAXDIM];
  ObitInfoList* list;
  gchar *routine = "CubeClipIn";

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
      
    } else if (strcmp(arg, "-inSeq") == 0) { /* input AIPS sequence number */
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
      
     } else if (strcmp(arg, "-BLC") == 0) { /* BLC */
      dim[0] = 1;
      /* read until something starts with "-" of hit end */
      i = 0;
      while ( ((ax+1)<argc) && (argv[ax+1][0]!='-')) {
	iarr[i++] = strtol(argv[++ax], NULL, 0);
      }
      dim[0] = i;
      ObitInfoListAlwaysPut (list, "BLC", OBIT_oint, dim, &iarr);
      
     } else if (strcmp(arg, "-TRC") == 0) { /* TRC */
      dim[0] = 1;
      /* read until something starts with "-" of hit end */
      i = 0;
      while ( ((ax+1)<argc) && (argv[ax+1][0]!='-')) {
	iarr[i++] = strtol(argv[++ax], NULL, 0);
      }
      dim[0] = i;
      ObitInfoListAlwaysPut (list, "TRC", OBIT_oint, dim, &iarr);
      
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
} /* end CubeClipIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: CubeClip -input file -output ofile [args]\n");
    fprintf(stderr, "Clip insignificant pixels in 3D image\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def CubeClip.inp\n");
    fprintf(stderr, "  -output output result file, def CubeClip.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType 'AIPS' or 'FITS' type for input image\n");
    fprintf(stderr, "  -inFile input FITS Image file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -BLC bottom-left (1-rel) pixel def. [1,1,..] \n");
    fprintf(stderr, "  -TRC top-right (1-rel) pixel def. [nx,ny,..] \n");
    fprintf(stderr, "  -outFile output FITS Image file\n");
    fprintf(stderr, "  -outName output AIPS file name [def inName]\n");
    fprintf(stderr, "  -outClass output AIPS file class [def '?CubCl']\n");
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
/*     blc       Int  [7]   bottom-left (1-rel) corner[def {1,1,...)]     */
/*     trc       Int  [7]   top-right (1-rel) corner[def {0,0,...)]       */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   i, itemp;
  olong   blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong   trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ObitInfoList *out = newObitInfoList();
  ofloat ftemp[5] = {0.0,0.0,0.0,0.0,0.0};
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
  strTemp = "CubeClipName";
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

  /* Parms */
  dim[0] = 4;
  for (i=0; i<4; i++) ftemp[i] = 0.0; ftemp[3] = 5.0;
  ObitInfoListPut (out, "Parms", OBIT_float, dim, ftemp, err);
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
  /*gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
    ofloat ftemp;*/
  ObitInfoList *out = newObitInfoList();
  /*gchar *routine = "defaultOutputs";*/

  /* add parser items */
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
/*  Define output image                                                   */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitImage for output image                                       */
/*----------------------------------------------------------------------- */
ObitImage* getOutputImage (ObitInfoList *myInput, ObitErr *err)
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
    ObitInfoListGetP(myInput, "inClass", &type, dim, (gpointer)&strTemp2);
    if  (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find - use "xCubCl" x=first char of inClass */
      strncpy (Aclass, "XSqish", 7);
      Aclass[0] = strTemp2[0];
    }
    /* If blank use inClass char 1 + 'Sqish' */
    if (!strncmp (Aclass, "    ", 4)) {
      strncpy (Aclass, "XCubCl", 7);
      Aclass[0] = strTemp2[0];
    }
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
    ObitImageSetAIPS (outImage, OBIT_IO_byPlane, disk, cno, AIPSuser, 
		      blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */
    /* output FITS file name */
    ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp2);
    if (ObitInfoListGetP(myInput, "outFile", &type, dim, (gpointer)&strTemp)) {
      strncpy (outFile, strTemp, 128);
    } else { 
      g_snprintf (outFile, 129, "CubeClip%s", strTemp2);
    }
    /* If blank use CubClp+inFile */
    if (!strncmp (outFile, "    ", 4)) {
      g_snprintf (outFile, 129, "CubClp%s", strTemp2);
    }
    
    /* output FITS disk default = inDisk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    idisk = disk;
    ObitInfoListGetTest(myInput, "outDisk", &type, dim, &disk);
    if (disk<=0) disk = idisk;

    /* Tell about it */
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

  return outImage;
} /* end getOutputImage */

/*----------------------------------------------------------------------- */
/*  Write History for CubeClip                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inImage   Image to copy history from                              */
/*      outImage  Image to write history to                               */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void CubeClipHistory (ObitInfoList* myInput, ObitImage* inImage, 
		    ObitImage* outImage, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "BLC",  "TRC",  "Parms", 
    NULL};
  gchar *routine = "CubeClipHistory";

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
 
} /* end CubeClipHistory  */

/*----------------------------------------------------------------------- */
/*  Clip image and copy tables                                            */
/*  Parms:                                                                */
/*      [0] min. RMS (convolved image)                                    */
/*      [1] min. fraction of  peak                                        */
/*      [2] >0.5 => blank fill, else zero fill                            */
/*      [3] Convolution size (0->5)                                       */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inImage   Image to copy history from                              */
/*      outImage  Image to write history to                               */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void CubeClipClip (ObitInfoList* myInput, ObitImage* inImage, 
		    ObitImage* outImage, ObitErr* err)
{
  oint         noParms;
  ObitIOCode   iretCode, oretCode;
  ofloat       RMS, maxF, newVal, *Parms=NULL, fblank =  ObitMagicF();
  ofloat       minAllow, FWHM, rescale;
  ObitFFT      *FFTfor=NULL, *FFTrev=NULL;
  ObitFArray   *ConvFn=NULL, *convlArray=NULL, *padConvFn=NULL, *padImage=NULL;
  ObitCArray   *wtArray=NULL, *FTArray=NULL;
  ObitInfoType type;
  ObitTableCC *inCC=NULL, *outCC=NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM], blc0[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM], trc0[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar        *tabType = "AIPS CC", *today=NULL;
  olong        i, inVer, outVer, highVer, lo, hi, pos[IM_MAXDIM], Cen[2], temp;
  olong        ndim=2, naxis[2], tblc[2], ttrc[2], cen[2];
  gchar        *exclude[]={"AIPS CC", "AIPS HI", "AIPS PL", "AIPS SL", 
			   NULL}; 
  gboolean     odd;
  gchar        *routine = "CubeClipClip";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  /* Control parameters */
  ObitInfoListGetP(myInput, "Parms", &type, dim, (gpointer)&Parms); 
  /* Get region from myInput */
  ObitInfoListGetTest(myInput, "BLC", &type, dim, blc);
  /* Defaults */
  for (i=0; i<dim[0]; i++) blc[i] = MAX (1, blc[i]);
  ObitInfoListGetTest(myInput, "TRC", &type, dim, trc);

  /* Open and close to get full size */
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (inImage->info, "BLC", OBIT_long, dim, blc0);
  ObitInfoListAlwaysPut (inImage->info, "TRC", OBIT_long, dim, trc0);
  iretCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  iretCode = ObitImageClose (inImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Set actual requested BLC, TRC */
  ObitInfoListAlwaysPut (inImage->info, "BLC", OBIT_long, dim, blc);
  for (i=0; i<3; i++) if (trc[i]<=0.0) trc[i] = inImage->myDesc->inaxes[i];
  /* First two must be even for convolution to work */
  for (i=0; i<2; i++) {
    temp = trc[i] - blc[i] + 1;
    odd = (2*(temp/2)) != temp;
    if (odd) { /* make size even - fiddle TRC is needed */
      if (trc[i]<inImage->myDesc->inaxes[i]) trc[i] += 1.0;
      else trc[i] -= 1.0;
    }
  }
  ObitInfoListAlwaysPut (inImage->info, "TRC", OBIT_long, dim, trc);

  /* Blank or zero for out of range points? */
  if (Parms[2]>=0.5) newVal = fblank;
  else newVal = 0.0;

  /* Clone output - will copy non CC tables 
     ObitImageClone (inImage, outImage, err);
     if (err->error) Obit_traceback_msg (err, routine, inImage->name);*/

  /* Open input image */
  iretCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Copy descriptor */
  outImage->myDesc = ObitImageDescCopy(inImage->myDesc, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Create FFTs */
  FFTfor = ObitFeatherUtilCreateFFT(inImage, OBIT_FFT_Forward);
  FFTrev = ObitFeatherUtilCreateFFT(inImage, OBIT_FFT_Reverse);

  /* Normalize rescale for FFT */
  rescale = 1.0 / (ofloat)(FFTfor->dim[0] * FFTfor->dim[1]);

  /* Setup for convolutions - make convolving function */
  ConvFn = ObitFArrayCreate("ConvFn", 2, inImage->myDesc->inaxes);
  if (Parms[3]<=0.0) FWHM = 5.0;
  else FWHM   = Parms[3];
  Cen[0] = inImage->myDesc->inaxes[0]/2; 
  Cen[1] = inImage->myDesc->inaxes[1]/2; 
  ObitFArray2DCGauss (ConvFn, Cen, FWHM);

   /* Padded convolving function */
  naxis[0] = FFTfor->dim[0];  naxis[1] = FFTfor->dim[1]; 
  padConvFn = ObitFArrayCreate("Pad Conv Fn", ndim, naxis);
  ObitFeatherUtilPadArray (FFTfor, ConvFn, padConvFn);

  /* Arrays for padded image/ FFT */
  naxis[0] = FFTfor->dim[0];  naxis[1] = FFTfor->dim[1]; 
  padImage = ObitFArrayCreate("Pad Image", ndim, naxis);
  FTArray  = ObitFeatherUtilCreateFFTArray (FFTfor);

  /* FFT Convolving function to wtArray */
  wtArray = ObitFeatherUtilCreateFFTArray (FFTfor);
  ObitFArray2DCenter (padConvFn); /* Swaparoonie to FFT order */
  ObitFFTR2C (FFTfor, padConvFn, wtArray);
  if (ConvFn)    ConvFn    = ObitFArrayUnref(ConvFn);
  if (padConvFn) padConvFn = ObitFArrayUnref(padConvFn);

  /* Creation date today */
  today = ObitToday();
  strncpy (outImage->myDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
 
  /* Force to float pixels */
  outImage->myDesc->bitpix=-32;

  /* Open output image */
  /* Use external buffer for writing output */
  outImage->extBuffer = TRUE;
  oretCode = ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Loop over planes until hitting EOF */
  while (iretCode!= OBIT_IO_EOF) {
    iretCode = ObitImageRead (inImage, NULL, err);
    if (iretCode == OBIT_IO_EOF) break;
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);

    /* Get plane statistics */
    maxF = fabs (ObitFArrayMaxAbs(inImage->image, pos));

    /* Convolve plane to get mask 
       First, pad image */
    ObitFeatherUtilPadArray (FFTfor, inImage->image, padImage);
    
    /* FFT padded image to FTArray */
    ObitFArray2DCenter (padImage); /* Swaparoonie to FFT order */
    ObitFFTR2C (FFTfor, padImage, FTArray);

    /* Multiply by transfer function */
    ObitCArrayMul (FTArray, wtArray, FTArray);

    /* Back FFT */
    ObitFFTC2R(FFTrev, FTArray, padImage);
    ObitFArray2DCenter (padImage);/* Swaparoonie */

    /* Get window to extract */
    cen[0] = FFTfor->dim[0]/2;  cen[1] = FFTfor->dim[1]/2; 
    tblc[0] = cen[0] - inImage->image->naxis[0] / 2; 
    ttrc[0] = cen[0] + inImage->image->naxis[0] / 2;
    ttrc[0] -= (ttrc[0]-tblc[0]+1) - inImage->image->naxis[0];
    tblc[1] = cen[1] - inImage->image->naxis[1] / 2; 
    ttrc[1] = cen[1] + inImage->image->naxis[1] / 2;
    ttrc[1] -= (ttrc[1]-tblc[1]+1) - inImage->image->naxis[1];
    
    /* Extract convolved image */
    convlArray = ObitFArraySubArr(padImage, tblc, ttrc, err);
    if (err->error) goto cleanup;

    /* rescale units */
    ObitFArraySMul (convlArray, rescale);

    /* Blank convolved image the same as Image */
    ObitFArrayBlank (convlArray, inImage->image, convlArray);
    
    /* Blank out of range pixels in convolved image */
    RMS  = ObitFArrayRMS(convlArray);
    minAllow = MAX (Parms[0]*RMS, Parms[1]*maxF);
    ObitFArrayInClip (convlArray, -minAllow, minAllow, fblank);

    /* Blank Image plane */
    ObitFArrayBlank (inImage->image, convlArray, inImage->image);

    /* newVal Replace blanks if needed */
    if (newVal != fblank) ObitFArrayDeblank (inImage->image, newVal);

    /* Write plane */
    oretCode = ObitImageWrite(outImage, inImage->image->array, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);

    /* Loop cleanup */
    if (convlArray) convlArray = ObitFArrayUnref(convlArray);
  } /* End loop over input image */

  /* Close input */
  iretCode = ObitImageClose (inImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Close output */
  oretCode = ObitImageClose (outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Free image buffer */
  inImage->image = ObitFArrayUnref(inImage->image);

  /* Unset external buffer for writing */
  outImage->extBuffer = FALSE;

  /* Cleanup */
 cleanup:
  if (padImage) padImage = ObitFArrayUnref(padImage);
  if (convlArray) convlArray = ObitFArrayUnref(convlArray);
  if (wtArray) wtArray = ObitCArrayUnref(wtArray);
  if (FTArray) FTArray = ObitCArrayUnref(FTArray);
  if (FFTfor) FFTfor = ObitFFTUnref(FFTfor);
  if (FFTrev) FFTrev = ObitFFTUnref(FFTrev);
 
  /* Copy nonCC tables */
  iretCode = ObitImageCopyTables (inImage, outImage, exclude, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Copy any relevant CCTables  */
  /* Number of Tables */
  highVer = ObitTableListGetHigh (inImage->tableList, tabType);
  lo = MAX (blc[2], 1);
  hi = MIN (trc[2], highVer);
  if (hi <= 0) hi = highVer;
  for (inVer=lo; inVer<=hi; inVer++) {
    noParms = 0;
    inCC = newObitTableCCValue ("inCC", (ObitData*)inImage, &inVer, OBIT_IO_ReadOnly, 
				noParms, err);
    noParms = inCC->noParms;
    outVer = inVer-blc[2]+1;
    outCC = newObitTableCCValue ("outCC", (ObitData*)outImage, &outVer, OBIT_IO_WriteOnly, 
				 noParms, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    outCC = ObitTableCCCopy (inCC, outCC, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    inCC = ObitTableCCUnref(inCC);
    outCC = ObitTableCCUnref(outCC);
  } /* end loop over tables */
  
} /* end CubeClipClip  */


