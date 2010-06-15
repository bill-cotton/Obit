/* $Id$  */
/* Convert Ionospheric movie to 2D structure func (distance,time)     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2007,2010                                          */
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

#include "ObitImage.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitFArray.h"


/* internal prototypes */
/* Get inputs */
ObitInfoList* IonSFIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void IonSFOut (ObitInfoList* outList, ObitErr *err);
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
/* Collapse */
void IonSFCollapse (ObitInfoList* myInput, ObitImage* inImage, 
		    ObitImage* outImage, ObitErr* err);
/* Write history */
void IonSFHistory (ObitInfoList* myInput, ObitImage* inImage, 
		   ObitImage* outImage, ObitErr* err);
/* determine structure function */
void ObitFArrayPStrFunc (ObitFArray* in, olong ncell, ofloat *sf, 
			 ofloat *work);

/* Program globals */
gchar *pgmName = "IonSF";       /* Program name */
gchar *infile  = "IonSF.inp";   /* File with program inputs */
gchar *outfile = "IonSF.out";   /* File to contain program outputs */
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
/*  onvert Ionospheric movie to 2D structure func (distance,time)         */
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitImage    *inImage = NULL, *outImage = NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = IonSFIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) goto exit;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Get input Image Object */
  inImage = getInputImage (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Define output Image Object */
  outImage = getOutputImage (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Collapse into 2D */
  IonSFCollapse (myInput, inImage, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Do history */
  IonSFHistory ( myInput, inImage, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inImage   = ObitUnref(inImage);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* IonSFIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "IonSFIn";

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
} /* end IonSFIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: IonSF -input file -output ofile [args]\n");
    fprintf(stderr, "Ionospheric structure function image from movie \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def IonSF.inp\n");
    fprintf(stderr, "  -output output result file, def IonSF.out\n");
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
    fprintf(stderr, "  -outClass output AIPS file class [def 'IonSF']\n");
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
/*     subCC     Boo        Select CC tables                              */
/*     mrgCC     Boo        Merge CC tables                               */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean btemp;
  gchar *strTemp;
  oint   itemp;
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
  strTemp = "IonSFName";
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
  dim[0] = 5;
  btemp = FALSE;
  ObitInfoListPut (out, "Parms", OBIT_float, dim, ftemp, err);
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
  olong         Aseq, disk, cno;
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
    ObitTrimTrail(inFile);  /* remove trailing blanks */
   
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
  gchar        *FITS = "FITS";
  gchar        *routine = "getOutputImage";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outImage;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input UV data Object */
  outImage = newObitImage("output Image");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4)))
    ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;
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
    if (!strncmp (Aname, "    ", 4)) strncpy (Aname, "No Name", 13);
    Aname[12] = 0;

    /* output AIPS class */
    ObitInfoListGetP(myInput, "inClass", &type, dim, (gpointer)&strTemp2);
    if  (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find - use  IonSF */
      strncpy (Aclass, "IonSF", 7);
    }
    if (!strncmp (Aclass, "    ", 4)) {
      strncpy (Aclass, "IonSF", 7);
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
    ObitImageSetAIPS (outImage, OBIT_IO_byPlane, disk, cno, AIPSuser, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */
    /* output FITS file name */
    ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp2);
    if (ObitInfoListGetP(myInput, "outFile", &type, dim, (gpointer)&strTemp)) {
      strncpy (outFile, strTemp, 128);
    } else { 
      g_snprintf (outFile, 129, "IonSF%s", strTemp2);
    }
    /* If blank use IonSF+inFile */
    if (!strncmp (outFile, "    ", 4)) {
      g_snprintf (outFile, 129, "IonSF%s", strTemp2);
    }
    ObitTrimTrail(outFile);  /* remove trailing blanks */
  
    /* output FITS disk default = inDisk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    idisk = disk;
    ObitInfoListGetTest(myInput, "outDisk", &type, dim, &disk);
    if (disk<0) disk = idisk;

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
/*  Write History for IonSF                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inImage   Image to copy history from                              */
/*      outImage  Image to write history to                               */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void IonSFHistory (ObitInfoList* myInput, ObitImage* inImage, 
		    ObitImage* outImage, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "BLC",  "TRC",
    NULL};
  gchar *routine = "IonSFHistory";

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
  } else if (outHistory->FileType==OBIT_IO_FITS) {
    ObitHistoryCopy2Header (inHistory, outHistory, err);
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

  /* If output FITS - copy to header */
  if (outHistory->FileType==OBIT_IO_FITS) {
    ObitHistoryCopy2Header (outHistory, outHistory, err);
  }

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end IonSFHistory  */

/*----------------------------------------------------------------------- */
/*  Determine structure function and collapse to 2D                       */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inImage   Input cube                                              */
/*      outImage  output image                                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void IonSFCollapse (ObitInfoList* myInput, ObitImage* inImage, 
		    ObitImage* outImage, ObitErr* err)
{
  ObitIOCode   iretCode, oretCode;
  ofloat       time, *work=NULL, *sf=NULL, fblank = ObitMagicF();
  olong         nx, ny;
  olong        i, naxis[2], pos[2], iplane;
  ObitFArray   *accum=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar        *today=NULL;
  gchar        *routine = "IonSFCollapse";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  /* Control parameters */
  /* Get region from myInput */
  ObitInfoListGetTest(myInput, "BLC", &type, dim, blc);
  /* Defaults */
  for (i=0; i<dim[0]; i++) blc[i] = MAX (1, blc[i]);
  ObitInfoListGetTest(myInput, "TRC", &type, dim, trc);

  /* Open input image */
  iretCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* how big can the structure function be? */
  nx = 1.414 * MAX (inImage->myDesc->inaxes[0], inImage->myDesc->inaxes[1]);

  /* How many rows? */
  ny = inImage->myDesc->inaxes[2];

  /* Copy descriptor */
  outImage->myDesc = ObitImageDescCopy(inImage->myDesc, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Creation date today */
  today = ObitToday();
  strncpy (outImage->myDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
 
  /* Output image size */
  outImage->myDesc->inaxes[0] = nx;
  outImage->myDesc->inaxes[1] = ny;

  outImage->myDesc->naxis = 2;
  /* Structure function axis */
  g_snprintf (outImage->myDesc->ctype[0], 9, "Offset  ");
  outImage->myDesc->crval[0]  = 0.0;
  outImage->myDesc->crpix[0]  = 1.0;
  outImage->myDesc->crota[0]  = 0.0;
  /* Time axis */
  g_snprintf (outImage->myDesc->ctype[1], 9, "Time    ");
  outImage->myDesc->crval[1]  = inImage->myDesc->crval[2];
  outImage->myDesc->cdelt[1]  = inImage->myDesc->cdelt[2];
  outImage->myDesc->crpix[1]  = 1.0;
  outImage->myDesc->crota[1]  = 0.0;
  /* Clear third axis */
  outImage->myDesc->inaxes[2] = 0;
  g_snprintf (outImage->myDesc->ctype[2], 9, "        ");
  outImage->myDesc->crval[2]  = 0.0;
  outImage->myDesc->cdelt[2]  = 0.0;
  outImage->myDesc->crpix[2]  = 0.0;
  outImage->myDesc->crota[2]  = 0.0;

 /* Force to float pixels */
  outImage->myDesc->bitpix=-32;
  g_snprintf (outImage->myDesc->origin, 9, "Obit    ");

  /* Make accumulation array */
  naxis[0] = nx; naxis[1] = ny;
  accum = ObitFArrayCreate("Accumulator", 2, naxis);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  ObitFArrayFill (accum, 0.0);

  /* Work array */
  work = g_malloc0(nx*sizeof(ofloat));

  /* Loop over planes */
  iplane = -1;
  pos[0] = 0;
  for (iplane=0; iplane<ny; iplane++) {
    iretCode = ObitImageRead (inImage, NULL, err);
    if (iretCode == OBIT_IO_EOF) break;
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);

    pos[1] = iplane;
    sf = ObitFArrayIndex (accum, pos); /* pointer in output buffer */

    /* save time (first pixel) */
    time =  inImage->image->array[0];

    /* Ignore first pixel (time) */
    inImage->image->array[0] = fblank;

    /* Get plane structure function */
    ObitFArrayPStrFunc (inImage->image, nx, sf, work);

    /* Save time in cell nx */
    sf[nx-1] = time;

  } /* End loop over input image */

   /* Close input */
  iretCode = ObitImageClose (inImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  /* Free image buffer */
  inImage->image = ObitFArrayUnref(inImage->image);

  /* Ensure outImage fully instantiated and OK */
  ObitImageFullInstantiate (outImage, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
   
  /* Open output image */
  /* Use external buffer for writing output */
  outImage->extBuffer = TRUE;
  oretCode = ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  /* Write plane */
  oretCode = ObitImageWrite(outImage, accum->array, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  oretCode = ObitImageClose (outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  /* Unset external buffer for writing */
  outImage->extBuffer = FALSE;

  if (work) g_free(work);  /* Cleanup */
  accum = ObitFArrayUnref(accum);
} /* end IonSFCollapse  */

/**
 * Structure function (actually square root of ) of a 2D array of phases
 * \param in    Input object with phases in deg
 * \param ncell Number of cells in sf, work
 * \param sf    [out] Structure function, unpopulated cells are blanked.
 * \param work  [out] Work array the size of sf
 */
void ObitFArrayPStrFunc (ObitFArray* in, olong ncell, ofloat *sf, 
			 ofloat *work)
{
  olong i, j, ix, iy, jx, jy, icell;
  ofloat dx, dy, cell, dph, fblank = ObitMagicF();
  
  /* error checks */
  g_assert (ObitFArrayIsA(in));
  
  /* Initialize */
  for (i=0; i<ncell; i++) {
    sf[i] = work[i] = 0.0;
  }
  
  /* Loop over array looking for non blank value */
  for (i=0; i<in->arraySize; i++) {
    if (in->array[i]!=fblank) {
      /* where am I? */
      iy = i / in->naxis[0];
      ix = i - iy*in->naxis[0];
      /* Search remainder */
      for (j=i+1; j<in->arraySize; j++) {
	if (in->array[j]!=fblank) {
	  /* difference phases */
	  dph = in->array[i]-in->array[j];
	  /* No more than half a turn */
	  if (dph> 180.0) dph -= 360.0;
	  if (dph<-180.0) dph += 360.0;
	  /* where am I? */
	  jy = j / in->naxis[0];
	  jx = j - jy*in->naxis[0];
	  dx = ix - jx;
	  dy = iy - jy;
	  cell = sqrt (dx*dx+dy*dy);
	  icell = (olong)(cell+0.5);
	  icell = MIN (icell, ncell-1);
	  sf[icell] += dph*dph;
	  work[icell]++;
	}
      } /* end inner loop over array */
    }
  } /* end outer loop over array */

  /* Normalize */
  for (i=0; i<ncell; i++) {
    if (work[i]>0.0) {
      sf[i] /= work[i];
      sf[i] = sqrt (sf[i]);
    } else {
      sf[i] = fblank;
    }
  }
} /* end  ObitFArrayPStrFunc */

