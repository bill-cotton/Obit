/* $Id$  */
/*  MCube: put together images into a cube                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2007-2011                                          */
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
#include "ObitImage.h"
#include "ObitFArray.h"
#include "ObitImageUtil.h"
#include "ObitHistory.h"
#include "ObitTableCC.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* MCubeIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void MCubeOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Digest inputs */
void digestInputs(ObitInfoList *myInput, ObitErr *err);
/* Get image from inputs */
void GetImages(ObitInfoList *myInput, gboolean *isNew, ObitErr *err);
/* MCube images together to outImage */
void doMCube (ObitInfoList *myInput, ObitImage *inImage, 
	       ObitImage *outImage, ObitErr *err);
/* Write History */
void doHistory (ObitInfoList *myInput, ObitImage *inImage, 
		ObitImage *outImage, gboolean isNew, ObitErr *err);
/* Generate output Image descriptor */
void MakeCubeDesc (ObitImageDesc *inDesc, ObitImageDesc *outDesc, 
		   olong axNum, olong axDim, 
		   ofloat axCRPix, ofloat axCDelt,
		   odouble axCRVal, ObitErr *err);
/* Insert input plane into output image */
void InsertPlane (ObitImage *in, ObitImage *out, olong *plane, 
		  ObitErr *err);
/* Combine header keyword frequency information */
void UpdateFreq(ObitInfoList *myInput, 
		ObitInfoList *oldInfo, ObitInfoList *newInfo, 
		ObitImage *in, ObitImage *out, 
		ObitErr *err);

/* Program globals */
gchar *pgmName = "MCube";       /* Program name */
gchar *infile  = "MCube.in" ;   /* File with program inputs */
gchar *outfile = "MCube.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitImage *inImage=NULL;           /* Input image */
ObitImage *outImage=NULL;          /* output image */
ObitInfoList *newInfo=NULL;        /* Input Desc.List */
ObitInfoList *oldInfo=NULL;        /* Old output Desc.List */
ObitInfoList *myInput  = NULL;     /* Input parameter list */
ObitInfoList *myOutput = NULL;     /* Output parameter list */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/* MCube:  put together images into a cube                                */
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitErr      *err= NULL;
  gboolean     isNew=FALSE;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = MCubeIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Get input image and output image */
  GetImages(myInput, &isNew, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* MCube them together */
  doMCube (myInput, inImage, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Update any header keyword frequency info */
  UpdateFreq (myInput, newInfo, oldInfo, inImage, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* History */
  doHistory (myInput, inImage, outImage, isNew, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* cleanup */
  inImage = ObitImageUnref(inImage);
  outImage= ObitImageUnref(outImage);
  myInput = ObitInfoListUnref(myInput); 
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput = ObitInfoListUnref(myOutput);
  newInfo  = ObitInfoListUnref(newInfo);
  oldInfo  = ObitInfoListUnref(oldInfo);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* MCubeIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "MCubeIn";

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
      
    } else if (strcmp(arg, "-axPix") == 0) { /* output pixel number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "axPix", OBIT_oint, dim, &itemp, err);
      
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
} /* end MCubeIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: MCube -input file -output ofile [args]\n");
    fprintf(stderr, "MCube Obit task = combine images into a cube\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def MCube.in\n");
    fprintf(stderr, "  -output output result file, def MCube.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input image\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inFile input FITS UV file\n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input image (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -axPix 1-rel pixel number on output image for input\n");
    fprintf(stderr, "  -outFile output image (FITS Image file\n");  
    fprintf(stderr, "  -outName output image (AIPS file name\n");
    fprintf(stderr, "  -outClass output image (AIPS file class\n");
    fprintf(stderr, "  -outSeq output image (AIPS file sequence\n");
    fprintf(stderr, "  -outDisk output image ((AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -out2File output uv FITS Image file\n");  
    fprintf(stderr, "  -AIPSdir single AIPS data directory\n");
    fprintf(stderr, "  -FITSdir single FITS data directory\n");
    
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
/*     inFile?   Str [?]    input FITS image file name [no def]           */
/*     inName?   Str [12]   input AIPS image name  [no def]               */
/*     inClass?  Str [6]    input AIPS image class  [no def]              */
/*     inSeq?    Int        input AIPS image sequence no  [no def]        */
/*     inDisk?   Int        input AIPS or FITS image disk no  [def 1]     */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  ofloat ftemp;
  oint   itemp;
  olong   jtemp;
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
  dim[0] = 12; dim[1] = 1;
  ObitInfoListPut (out, "inName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* input AIPS file class */
  dim[0] = 6; dim[1] = 1;
  ObitInfoListPut (out, "inClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* AIPS sequence */
  dim[0] = 1;dim[1] = 1;
  jtemp = 0; 
  ObitInfoListPut (out, "inSeq", OBIT_long, dim, &jtemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* AIPS or FITS disk number */
  dim[0] = 1;dim[1] = 1;
  jtemp = 1; 
  ObitInfoListPut (out, "inDisk", OBIT_long, dim, &jtemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Default output */
  /* FITS file name */
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* AIPS file name */
  dim[0] = 12; dim[1] = 1;
  ObitInfoListPut (out, "outName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* AIPS file class */
  dim[0] = 6; dim[1] = 1;
  ObitInfoListPut (out, "outClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* AIPS sequence */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "outSeq", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* AIPS or FITS disk number */
  dim[0] = 1;dim[1] = 1;
  jtemp = 1; 
  ObitInfoListPut (out, "outDisk", OBIT_long, dim, &jtemp, err);
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
  /*gchar *routine = "defaultOutputs";*/

  if (err->error) return out;  /* existing error */

  return out;
} /* end defaultOutputs */

/*----------------------------------------------------------------------- */
/*  Digest inputs                                                         */
/*   Set default input BLC, TRC = all                                     */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void digestInputs(ObitInfoList *myInput, ObitErr *err)
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ObitImageDesc *desc ;
  olong         i;
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* Get region from myInput */
  ObitInfoListGetTest(myInput, "BLC", &type, dim, blc); /* BLC */
  ObitInfoListGetTest(myInput, "TRC", &type, dim, trc); /* TRC */
  
  /* Set default input BLC, TRC */
  desc = (ObitImageDesc*)inImage->myIO->myDesc;
  for (i=0; i<IM_MAXDIM; i++) {
    blc[i] = MAX (1,  blc[i]);
    if (trc[i]<=0) trc[i] = desc->inaxes[i];
    trc[i] = MIN (trc[i], desc->inaxes[i]);
  }

  /* Save blc, trc */
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (myInput, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (inImage->info, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (myInput, "TRC", OBIT_long, dim, trc);
  ObitInfoListAlwaysPut (inImage->info, "TRC", OBIT_long, dim, trc);
} /* end digestInputs */

/*----------------------------------------------------------------------- */
/*  Get images                                                            */
/*  Values:                                                               */
/*      myInput   InfoList with inputs                                    */
/*      inImage   [out] ObitImage pointers                                */
/*                 Passed as global                                       */
/*      outImage  [out] Output ObitImage pointer                          */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void GetImages(ObitInfoList *myInput, gboolean *isNew, ObitErr *err)
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong        j, k, Aseq, disk, cno;
  gboolean     exist;
  gchar        *strTemp=NULL, inFile[128], *FITS="FITS";
  gchar        iAname[13], iAclass[7];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gchar        tname[101], *Type = NULL;
  olong        axNum=3, axDim=1, axPix=1;
  ofloat       axCRPix=1.0, axCDelt=-1000.0;
  odouble      axCRVal=-1000.0;
  gchar *routine = "GetImage";

  if (err->error) return;  /* existing error? */

  /* Get input region from myInput */
  ObitInfoListGetTest(myInput, "BLC", &type, dim, blc); /* BLC */
  ObitInfoListGetTest(myInput, "TRC", &type, dim, trc); /* TRC */

  /* Output file info */
  ObitInfoListGetTest(myInput, "axNum", &type, dim, &axNum); 
  ObitInfoListGetTest(myInput, "axDim", &type, dim, &axDim); 
  ObitInfoListGetTest(myInput, "axPix", &type, dim, &axPix);
  ObitInfoListGetTest(myInput, "axCRPix", &type, dim, &axCRPix);
  ObitInfoListGetTest(myInput, "axCDelt", &type, dim, &axCDelt);
  ObitInfoListGetTest(myInput, "axCRVal", &type, dim, &axCRVal);
  /* Default axis to expand = 3 */
  if (axNum<1) {
    axDim = 1;
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "axNum", OBIT_float, dim, &axNum);
  }

  /* Default dimension = 1 */
  if (axDim<1) {
    axDim = 1;
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "axDim", OBIT_float, dim, &axDim);
  }

  /* Default pixel = 1 */
  if (axPix<1) {
    axPix = 1;
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "axPix", OBIT_float, dim, &axPix);
  }

  /* Is this the first pixel? */
  *isNew = *isNew || (axPix==1);

  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */

    /* input AIPS disk */
    ObitInfoListGet(myInput,"inDisk", &type, dim, &disk, err);

    /* input AIPS name */
    for (k=0; k<12; k++) iAname[k] = ' '; iAname[k] = 0;
    ObitInfoListGet(myInput, "inName", &type, dim, iAname, err);

    /* input AIPS class */
    for (k=0; k<6; k++) iAclass[k] = ' '; iAclass[k] = 0;
    ObitInfoListGet(myInput, "inClass", &type, dim, iAclass, err);

    /* input AIPS sequence */
    ObitInfoListGet(myInput, "inSeq", &type, dim, &Aseq, err);
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
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);
    
    /* Make sure it's OK */
    ObitImageFullInstantiate (inImage, TRUE, err);
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */

    /* input FITS file name */
    for (j=0; j<128; j++) inFile[j] = 0;
    ObitInfoListGet(myInput, "inFile", &type, dim, inFile, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /* input FITS disk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /*  Object name from FITS name */
    inImage = newObitImage(inFile);
    
    /* define image */
    ObitImageSetFITS (inImage, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);
    
    /* Make sure it's OK */
    ObitImageFullInstantiate (inImage, TRUE, err);
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Image type %s", 
                   pgmName, strTemp);
  }
  if (err->error) Obit_traceback_msg (err, routine, routine);

  /* Save Descriptor list */
  newInfo = ObitInfoListCopy(inImage->myDesc->info);
   
  /* Output image */

  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4)))
    ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
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
    if (!strncmp(Aclass,"      ",6)) strcpy (Aclass, iAclass);

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
    *isNew = *isNew || (!exist);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /* Generate Object name from AIPS name */
    g_snprintf (tname, 100, "%s.%s:%d.%d", Aname, Aclass, Aseq, disk);
    outImage = newObitImage(tname);

    /* reset BLC, TRC */
    for (j=0; j<IM_MAXDIM; j++)  {blc[j] = 1; trc[j] = 0;}
    
    /* define image */
    ObitImageSetAIPS (outImage, OBIT_IO_byPlane, disk, cno, AIPSuser, 
		      blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output AIPS image %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */

    /* FITS file name */
    for (j=0; j<128; j++) inFile[j] = 0;
    ObitInfoListGet(myInput, "outFile", &type, dim, inFile, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /*  FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /*  Object name from FITS name */
    outImage = newObitImage(inFile);
    
    /* reset BLC, TRC */
    for (j=0; j<IM_MAXDIM; j++)  {blc[j] = 1; trc[j] = 0;}

    /* define image */
    ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output FITS Image %s on disk %d", inFile, disk);
    /* Check if file exists */
    if (!ObitFileExist(inFile, err)) *isNew = TRUE;

  } /* end FITS output */

  /* Default reference pixel the same as input+BLC-1 */
  if ((axCRPix<-999.9) && (axCRPix>-1001.0)) {
    axCRPix = inImage->myDesc->crpix[axNum-1]+blc[axNum-1]-1.0;
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "axCRPix", OBIT_float, dim, &axCRPix);
  }

  /* Default increment the same as input */
  if ((axCDelt<-999.9) && (axCDelt>-1001.0)) {
    axCDelt = inImage->myDesc->cdelt[axNum-1];
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "axCDelt", OBIT_float, dim, &axCDelt);
  }

  /* Default reference value the same as input */
  if ((axCRVal<-999.9) && (axCRVal>-1001.0)) {
    axCRVal = inImage->myDesc->crval[axNum-1];
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "axCRVal", OBIT_double, dim, &axCRVal);
  }

  /* Initialize output image */
  if (*isNew) MakeCubeDesc (inImage->myDesc, outImage->myDesc, 
			    axNum, axDim, axCRPix, axCDelt,
			    axCRVal, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
   
  /* Make sure it's OK */
  if (*isNew) { /* Open and close */
    ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
    ObitImageClose (outImage, err);
  } else { /* test */
    ObitImageOpen (outImage, OBIT_IO_ReadWrite, err);
    ObitImageClose (outImage, err);
    oldInfo = ObitInfoListCopy(outImage->myDesc->info);
  }
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
} /* end GetImages */

/*----------------------------------------------------------------------- */
/*  Insert inImage into outImage                                          */
/*  Values:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inImage   input ObitImage                                         */
/*      outImage  Output ObitImage pointer                                */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void doMCube (ObitInfoList *myInput, ObitImage *inImage, ObitImage *outImage, 
	      ObitErr *err)
{
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  oint         noParms=0;
  gboolean     copyCC=FALSE;
  ObitTableCC  *inCC=NULL, *outCC=NULL;
  ObitInfoType type;
  olong        plane[5] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        inVer, outVer, highVer;
  gchar        *tabType = "AIPS CC";
  olong axNum=3, axPix=1;
  gchar *routine = "doMCube";

  if (err->error) return;  /* existing error? */

  /* Which plane? */
  ObitInfoListGetTest(myInput, "axNum", &type, dim, &axNum); 
  ObitInfoListGetTest(myInput, "axPix", &type, dim, &axPix);
  plane[axNum-3] = axPix;

  /* Insert */
  ObitImageUtilInsertCube (inImage, outImage, plane, axNum, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Copy CC Tables? */
  ObitInfoListGetTest(myInput, "copyCC", &type, dim, &copyCC); 
  if (copyCC) {
    ObitInfoListGetTest(myInput, "BLC", &type, dim, blc);
    
    inVer   = MAX (1, blc[axNum-1]);
    outVer  = MAX (1, axPix);
    highVer = ObitTableListGetHigh (inImage->tableList, tabType);
    while (inVer<=highVer) {
      inCC = newObitTableCCValue ("inCC", (ObitData*)inImage, &inVer, 
				  OBIT_IO_ReadOnly, noParms, err);
      if (inCC==NULL) return;
      noParms = inCC->noParms;
      outCC = newObitTableCCValue ("outCC", (ObitData*)outImage, &outVer, 
				   OBIT_IO_WriteOnly, noParms, err);
      if (err->error) Obit_traceback_msg (err, routine, outImage->name);
      if (outCC==NULL) return;

      outCC = ObitTableCCCopy (inCC, outCC, err);
      if (err->error) Obit_traceback_msg (err, routine, outImage->name);
      inCC = ObitTableCCUnref(inCC);
      outCC = ObitTableCCUnref(outCC);

      /* Tell about it */
      Obit_log_error(err, OBIT_InfoErr, "Copied CC table %d to  %d",
		     inVer, outVer);
      inVer++;
      outVer++;
    } /* end copy loop */
  } /* End of copy CC table */
} /* end doMCube */

/*----------------------------------------------------------------------- */
/*  Write history                                                         */
/*  Values:                                                               */
/*      myInputs  input parameter array                                   */
/*      inImage   Array of ObitImage pointers                             */
/*      outImage  Output ObitImage pointer                                */
/*      isNew     True if file newly created and input history is copied. */
/*                Otherwise only new history written                      */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void doHistory (ObitInfoList *myInput, ObitImage *inImage, 
		ObitImage *outImage, gboolean isNew, ObitErr *err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "BLC", "TRC", "copyCC",
    "axNum", "axDim", "axPix", "axCRPix", "axCDelt", "axCRVal", 
    NULL};
  gchar *routine = "doHistory";

  if (err->error) return;  /* existing error? */

  /* Open outImage to make sure History gets registered in header if created */
  ObitImageOpen (outImage, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Do history  */
  if (isNew) {  /* Copy old history */
    inHistory  = newObitDataHistory ((ObitData*)inImage, OBIT_IO_ReadOnly, err);
    outHistory = newObitDataHistory ((ObitData*)outImage, OBIT_IO_WriteOnly, err);
    if (inHistory->FileType==OBIT_IO_FITS) {
      ObitHistoryCopyHeader (inHistory, outHistory, err);
    } else { /* simply copy history */
      ObitHistoryCopy (inHistory, outHistory, err);
    }
    inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  } else { /* Add to old */
    outHistory = newObitDataHistory ((ObitData*)outImage, OBIT_IO_ReadWrite, err);
  }
  
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

  outHistory = ObitHistoryUnref(outHistory); /* cleanup */
 
} /* end doHistory */

/**
 * Set ImageDesc values for a cube, assumes that one axis of 
 * input is to be expanded.
 * This should be called before the image is Opened or instantiated.
 * \param inDesc     Input Image Descriptor.
 * \param outDesc    Output Image Descriptor 
 * \param axNum      1-rel axis number of input image which is to 
 *                   be expanded to form the cube,
 * \param axDim      The output dimensionality of axis axNum
 * \param axCRPix    Reference pixel on axis axNum.
 * \param axCDelt    Coordinate increment for axis axNum.
 * \param axCRVal    Coordinate reference value on axis axNum
 * \param err        Error stack, returns if error.
 */
void MakeCubeDesc (ObitImageDesc *inDesc, ObitImageDesc *outDesc, 
		   olong axNum, olong axDim, 
		   ofloat axCRPix, ofloat axCDelt,
		   odouble axCRVal,
		   ObitErr *err)
{
  gchar *name, *today=NULL;
  gchar *routine = "ObitImageUtilMakeCube";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageDescIsA(inDesc));

  /* Save output name */
  if (outDesc->name) name = g_strdup (outDesc->name);
  else  name = g_strdup ("Descriptor");

  /* Most info from inDesc */
  outDesc = ObitImageDescCopy (inDesc, outDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inDesc->name);

  /* Creation date today */
  today = ObitToday();
  strncpy (outDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
 
  /* restore name */
  if (outDesc->name) g_free(outDesc->name);
  outDesc->name = name;

  /* Set output axis */
  outImage->myDesc->inaxes[axNum-1] = axDim;
  outImage->myDesc->crpix[axNum-1]  = axCRPix;
  outImage->myDesc->cdelt[axNum-1]  = axCDelt;
  outImage->myDesc->crval[axNum-1]  = axCRVal;


  /* reset image max/min */
  outDesc->maxval    = -1.0e20;
  outDesc->minval    =  1.0e20;

  /* 32 bit float */
  outDesc->bitpix = -32;

  return;
} /* end MakeCubeDesc */

/*----------------------------------------------------------------------- */
/*  Combine header keyword frequency information                          */
/*  If "NSPEC" not found on input images descriptor list, return          */
/*  Values:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      newInfo   Input Desc List                                         */
/*      oldInfo   Existing output Desc List                               */
/*      in        Input ObitImage                                         */
/*      out       Output ObitImage                                        */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void UpdateFreq(ObitInfoList *myInput, 
		ObitInfoList *newInfo, ObitInfoList *oldInfo, 
		ObitImage *in, ObitImage *out, 
		ObitErr *err)
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat farr[10], inAlpha, outAlpha;
  odouble *inSpecFreq=NULL, *outSpecFreq=NULL, *combSpecFreq=NULL;
  odouble darr[10], inAlphaRefF, outAlphaRefF;
  gchar keyword[12];
  olong i, j, nTerm, axPix, inSpec, inTerm, onSpec, onTerm, cnSpec, cnTerm;
  gchar *routine = "UpdateFreq";

  if (err->error) return;  /* existing error? */

  /* Number of spectral channels */
  inSpec = 0;
  ObitInfoListGetTest (newInfo, "NSPEC", &type, dim, &inSpec);
  if (inSpec<=0) return;

  /* Open image */
  ObitImageOpen (out, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
  out->myStatus = OBIT_Modified;  /* Force update on disk */

  /* Number of input spectral terms */
  inTerm = 0;
  ObitInfoListGetTest (newInfo, "NTERM", &type, dim, &inTerm);

  /* Number of spectral terms  */
  nTerm = inTerm;
  ObitInfoListGetTest (myInput, "nTerm", &type, dim, &nTerm);
  if (nTerm<=0) return;  /* Want some? */
  if (nTerm>0) cnTerm = nTerm;

  /* Where inserted  */
  axPix = nTerm+1;
  ObitInfoListGetTest (myInput, "axPix", &type, dim, &axPix);

  /* Alpha */
  farr[0] = 0.0;
  ObitInfoListGetTest (newInfo, "ALPHA", &type, dim, farr);
  inAlpha = farr[0];

  /* Alpha reference frequency - default to reference frequency */
  darr[0] = in->myDesc->crval[in->myDesc->jlocf];
  ObitInfoListGetTest (newInfo, "ALPHARF", &type, dim, darr);
  inAlphaRefF = darr[0];

  /* Create frequency array */
  inSpecFreq  = g_malloc0(inSpec*sizeof(odouble));

  /* Fetch frequencies */
  for (i=0; i<inSpec; i++) {
    inSpecFreq[i] = 1.0;
    sprintf (keyword, "FREQ%4.4d",i+1);
    ObitInfoListGetTest (newInfo, keyword, &type, 
			 dim, &inSpecFreq[i]);
  }
  
  /* See if they are already on output */
  onSpec = 0;
  if (oldInfo) {
    /* Number of spectral channels */
    ObitInfoListGetTest (oldInfo, "NSPEC", &type, dim, &onSpec);
    if (onSpec>0) {
      
      /* Number of input spectral terms */
      onTerm = 0;
      ObitInfoListGetTest (oldInfo, "NTERM", &type, dim, &onTerm);
      cnTerm = onTerm;
      
      /* Alpha */
      farr[0] = 0.0;
      ObitInfoListGetTest (oldInfo, "ALPHA", &type, dim, farr);
      outAlpha = farr[0];
      
      /* Alpha reference frequency - default to reference frequency */
      darr[0] = out->myDesc->crval[out->myDesc->jlocf];
      ObitInfoListGetTest (oldInfo, "ALPHARF", &type, dim, darr);
      outAlphaRefF = darr[0];
      
      /* Create frequency array */
      outSpecFreq  = g_malloc0(onSpec*sizeof(odouble));
      
      /* Fetch frequencies */
      for (i=0; i<onSpec; i++) {
	outSpecFreq[i] = 1.0;
	sprintf (keyword, "FREQ%4.4d",i+1);
	ObitInfoListGetTest (oldInfo, keyword, &type, 
			     dim, &outSpecFreq[i]);
      }
      
      /* IF both alphas not 0.0 then AlphaRefFs must be the same */
      if (inAlpha!=outAlpha){
	Obit_log_error(err, OBIT_InfoWarn, 
		       "%s: Prior ALPHAs unequal, %f != %f", 
		       routine, inAlpha, outAlpha);
      }
      if ((inAlpha==outAlpha) && (inAlpha!=0.0) && (inAlphaRefF!=outAlphaRefF)) {
	Obit_log_error(err, OBIT_InfoWarn, 
		       "%s: Prior ALPHAref freq  unequal, %lf != %lf", 
		       routine, inAlphaRefF, outAlphaRefF);
      }
    }
  } /* end if exist already on output */

  /* Make combined output */
  cnSpec = out->myDesc->inaxes[out->myDesc->jlocf] - nTerm;
  combSpecFreq  = g_malloc0(cnSpec*sizeof(odouble));
  if ((onSpec>0) && outSpecFreq) {
    j = 0;
    for (i=0; i<onSpec; i++) combSpecFreq[j++] = outSpecFreq[i];
  }
  j = axPix-cnTerm-1;
  for (i=0; i<inSpec; i++) combSpecFreq[j++] = inSpecFreq[i];
  
  /* Update header list */
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut (out->myDesc->info, "NSPEC", OBIT_long, dim, &cnSpec);
  ObitInfoListAlwaysPut (out->myDesc->info, "NTERM", OBIT_long, dim, &cnTerm);
  
  /* Add frequency info to descriptor */
  for (i=0; i<cnSpec; i++) {
    sprintf (keyword, "FREQ%4.4d",i+1);
    ObitInfoListAlwaysPut (out->myDesc->info, keyword, OBIT_double, 
			   dim, &combSpecFreq[i]);
  }

  /* New info output? */
  if (onSpec<=0) {
    /* Save Alpha */
    ObitInfoListAlwaysPut (out->myDesc->info, "ALPHA", OBIT_float, dim, &inAlpha);
    
    /* Save Alpha reference frequency */
    ObitInfoListAlwaysPut (out->myDesc->info, "ALPHARF", OBIT_double, dim, &inAlphaRefF);
  }

  /* Open image */
  ObitImageClose (out, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);

  /* Cleanup */
  if (inSpecFreq)   g_free(inSpecFreq);
  if (outSpecFreq)  g_free(outSpecFreq);
  if (combSpecFreq) g_free(combSpecFreq);
 
} /* end  UpdateFreq */
