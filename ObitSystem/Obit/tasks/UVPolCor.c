/* $Id$  */
/* Task to correct off-axis instrumental polarization in UV data      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009-2025                                          */
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

#include "ObitImageMosaic.h"
#include "ObitSkyModel.h"
#include "ObitSkyModelVMBeamMF.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitTableCCUtil.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitFITS.h"
#include "ObitImageInterp.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* UVPoCoIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void UVPoCoOut (ObitInfoList* outList, ObitErr *err);
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
/* Get input sky model */
ObitSkyModelVMBeamMF* getInputSkyModel (ObitInfoList *myInput,  ObitUV* inData, 
					ObitErr *err);
/* Create output uvdata */
ObitUV* setOutputData (ObitInfoList *myInput, ObitUV* inData, ObitErr *err);
/* Write history */
void UVPoCoHistory (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		   ObitErr* err);
/* Get beam images*/
void getBeam (ObitInfoList *myInput, gboolean doCmplx, ofloat Stokes0, olong *numAntType, 
	      ObitImage ***RXpol, ObitImage ***LYpol, ObitImage ***RLpol, ObitImage ***LRpol, 
	      ObitImage ***RXpolIm, ObitImage ***LYpolIm, ObitImage ***RLpolIm, ObitImage ***LRpolIm, 
	      ofloat **Diams, ObitErr *err);

/* Program globals */
gchar *pgmName = "UVPolCor";       /* Program name */
gchar *infile  = "UVPolCor.in" ;   /* File with program inputs */
gchar *outfile = "UVPolCor.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------------*/
/*   Obit Task to correct for instrumental polarization and remove Ipol model */
/*----------------------------------------------------------------------------*/
{
  oint         ierr = 0;
  ObitSystem   *mySystem=NULL;
  ObitUV       *inData=NULL, *outData=NULL;
  ObitSkyModelVMBeamMF *skyModel=NULL;
  gchar        *opcode=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitErr      *err= NULL;

   /* Startup - parse command line, read inputs */
  err = newObitErr();
  myInput = UVPoCoIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Digest inputs */
  digestInputs(myInput, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Get input uvdata */
  inData = getInputData (myInput, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Get input sky model */
  skyModel = getInputSkyModel (myInput, inData, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Get output uvdata */
  outData = setOutputData (myInput, inData, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Process */
  ObitInfoListGetP(myInput, "Opcode", &type, dim, (gpointer)&opcode);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  if (!strncmp (opcode, "DIV", 3)) {
    /* Divide */
       ObitSkyModelDivUV ((ObitSkyModel*)skyModel, inData, outData, err);
       if (err->error) Obit_log_error(err, OBIT_Error, "ERROR dividing");

    /* Anything else is subtract ("MODL" trapped elsewhere */
  } else {
    /* Subtract */
    ObitSkyModelSubUV ((ObitSkyModel*)skyModel, inData, outData, err);
    if (err->error) Obit_log_error(err, OBIT_Error, "ERROR subtracting");
  }
  /* anything in output? */
  if (outData->myDesc->nvis<10) 
    Obit_log_error(err, OBIT_Error, "NO data written to output");
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* History */
  UVPoCoHistory (myInput, inData, outData, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* show any messages and errors */
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  outData   = ObitUnref(outData);
  skyModel  = ObitUnref(skyModel);
 
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* UVPoCoIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "UVPoCoIn";

  /* error checks */
  g_assert(ObitErrIsA(err));
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
      
     } else if (strcmp(arg, "-DataType") == 0) { /* Data type AIPS or FITS */
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
      
    } else if (strcmp(arg, "-in2Seq") == 0) { /* AIPS sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "in2Seq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-in2Disk") == 0) { /* output uv disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "in2Disk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-in2Name") == 0) { /* AIPS  in2Name */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "in2Name", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-in2Class") == 0) { /* AIPS uin2Class */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "in2Class", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-in2File") == 0) { /* in2File (FITS) */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "in2File", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-in3Seq") == 0) { /* AIPS beam sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "in3Seq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-in3Disk") == 0) { /* beam disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "in3Disk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-in3Name") == 0) { /* AIPS beam in3Name */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "in3Name", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-in3Class") == 0) { /* AIPS beam in3Class */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "in3Class", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-in3File") == 0) { /* in3File (FITS beam) */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "in3File", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outSeq") == 0) { /* AIPS output UV sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outSeq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-outDisk") == 0) { /* output UV disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outDisk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-outName") == 0) { /* AIPS UV outName */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outName", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-outClass") == 0) { /* AIPS UV outClass */
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
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  return list;
} /* end UVPoCoIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: UVPolCor -input file -output ofile [args]\n");
    fprintf(stderr, "UVPolCor: Correct off-axis instrumental polarization\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def UVPolCor.in\n");
    fprintf(stderr, "  -output output result file, def UVPolCor.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input \n");
    fprintf(stderr, "  -inFile input FITS UV file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input uvdata (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -in2File input FITS Image file\n");  
    fprintf(stderr, "  -in2Name input AIPS Image file name\n");
    fprintf(stderr, "  -in2Class input AIPS Image file class\n");
    fprintf(stderr, "  -in2Seq input AIPS Image sequence\n");
    fprintf(stderr, "  -in2Disk input Image disk number (1-rel) \n");
    fprintf(stderr, "  -in3File input FITS Beam file\n");  
    fprintf(stderr, "  -in3Name input AIPS Beam file name\n");
    fprintf(stderr, "  -in3Class input AIPS Beam file class\n");
    fprintf(stderr, "  -in3Seq input AIPS Beam sequence\n");
    fprintf(stderr, "  -in3Disk input Beam disk number (1-rel) \n");
    fprintf(stderr, "  -outFile output uv FITS  file\n");  
    fprintf(stderr, "  -outName output uv AIPS file name\n");
    fprintf(stderr, "  -outClass output uv AIPS file class\n");
    fprintf(stderr, "  -outSeq output uv AIPS file sequence\n");
    fprintf(stderr, "  -outDisk output uv (AIPS or FITS) disk number (1-rel) \n");
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
/*     inDisk    Int        input AIPS or FITS uv disk no  [def 1]        */
/*     inName    Str [12]   input AIPS uv name  [no def]                  */
/*     inClass   Str [6]    input AIPS uv class  [no def]                 */
/*     inSeq     Int        input AIPS uv sequence no  [no def]           */
/*     channel   Int [1]    channel number, 0=>all                        */
/*     BIF       Int [1]    first IF to process, 0=>1                     */
/*     EIF       Int [1]    highest IF to process, 0=>all                 */
/*     in2Disk   Int        input AIPS or FITS image disk no  [def 1]     */
/*     in2File   Str [?]    input FITS image file name [def "Image.fits"  */
/*     in2Name   Str [12]   input AIPS image name  [no def]               */
/*     in2Class  Str [6]    input AIPS image class  [no def]              */
/*     in2Seq    Int        input AIPS image  sequence no  [no def]       */
/*     in2Disk   Int        input AIPS or FITS image disk no  [def 1]     */
/*     nmaps     Int [1]    number of fields in sky model                 */
/*     CCVer     Int [1]    CC file ver. number.  0 => highest.           */
/*     BComp     Int [64]   First clean component to process per field    */
/*     EComp     Int [64]   Highest clean component to process per field  */
/*     Flux      float [1]  Only components > Flux are used in the model. */
/*     outDisk   Int        output AIPS or FITS image disk no  [def 1]    */
/*     outFile   Str [?]    output FITS image file name [def "Image.fits" */
/*     outName   Str [12]   output AIPS image name  [no def]              */
/*     outClass  Str [6]    output AIPS image class  [no def]             */
/*     outSeq    Int        output AIPS image sequence no  [no def]       */
/*     Cmethod   Str [4]    "DFT ", "GRID", "    "=> fastest              */
/*     Cmodel    Str [4]    "COMP"=Clean comps, "IMAG"=use image          */
/*     Factor    Flt [1]    Model factor, -1 => add                     . */
/*     Opcode    Str [4]    "DIV " => divide "MODL" =>replace data        */
/*     mrgCC     Boo [1]    Merge CC tables                               */
/*     PBCor     Boo [1]    pri. beam corr?, def=True                     */
/*     antSize   Flt [1]    effective diameter (m) of primary, def=25     */
/*     Alpha     Flt (1)    default spectral index (0)                    */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   i, itemp, iarray[64];
  ofloat ftemp;
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
  strTemp = "UVPolCor.intab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS input uv file name */
  strTemp = "UVPolCorName";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS input uv file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* AIPS input uv sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "inSeq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* AIPS or FITS input uv disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "inDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* channel  */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "channel", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* BIF  */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "BIF", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* EIF  */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "EIF", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Input FITS image file name root*/
  strTemp = "UVPolCorModel.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "in2File", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Input AIPS image file name */
  strTemp = "UV";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "in2Name", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Input AIPS image file class root */
  strTemp = "IMAGER";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "in2Class", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Input AIPS image sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "in2Seq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Input AIPS or FITS image disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "in2Disk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* nmaps */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "nmaps", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* CCVer */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "CCVer", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* BComp */
  dim[0] = 64;dim[1] = 1;
  for (i=0; i<64; i++) iarray[i] = 1;
  ObitInfoListPut (out, "BComp", OBIT_oint, dim, &iarray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* EComp */
  dim[0] = 64;dim[1] = 1;
  for (i=0; i<64; i++) iarray[i] = 0;
  ObitInfoListPut (out, "EComp", OBIT_oint, dim, &iarray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* Flux */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "Flux", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output FITS UV file name */
  strTemp = "UVPolCorOut.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS UV file name */
  strTemp = "UVPolCorOut";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS UV file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS UV sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "outSeq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output AIPS or FITS UV disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "outDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Cmethod */
  strTemp = "    ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Cmethod", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Cmodel */
  strTemp = "    ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Cmodel", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Factor  */
  dim[0] = 1;dim[1] = 1;
  ftemp = 1.0; 
  ObitInfoListPut (out, "Factor", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Opcode */
  strTemp = "    ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Opcode", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* mrgCC */
  dim[0] = 1;
  btemp = FALSE;
  ObitInfoListPut (out, "mrgCC", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* PBCor, def=True */
  dim[0] = 1; dim[1] = 1;
  btemp = TRUE;
  ObitInfoListPut (out, "PBCor", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* antSize  */
  dim[0] = 1;dim[1] = 1;
  ftemp = 25.0; 
  ObitInfoListPut (out, "antSize", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* do3D, Always FALSE */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListPut (out, "do3D", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

   /* Stokes always "   " */
  strTemp = "    ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Stokes", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* default Spectral index */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "Alpha", OBIT_float, dim, &ftemp, err);
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
  gchar *strTemp, *opcode=NULL, outFile[129], Aname[13], Aclass[7];
  gboolean replace, doCalSelect, btemp;
  olong doCalib, i, n, disk, seq;
  ObitSkyModelMode modelMode;
  ObitSkyModelType modelType;
  ofloat modelFlux, Factor;
  gchar *routine = "digestInputs";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

  /* Convert test Cmethod to enum  Mode */
  ObitInfoListGetP (myInput, "Cmethod", &type, dim, (gpointer)&strTemp);
  if (!strncmp (strTemp, "GRID", 4)) modelMode = OBIT_SkyModel_Grid;
  else if (!strncmp (strTemp, "DFT", 3)) modelMode = OBIT_SkyModel_DFT;
  else modelMode = OBIT_SkyModel_Fastest;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "Mode", OBIT_long, dim, &modelMode);

  /* Convert test Cmodel to enum  ModelType */
  ObitInfoListGetP (myInput, "Cmodel", &type, dim, (gpointer)&strTemp);
  modelFlux = 0.0;
  ObitInfoListGetTest (myInput, "modelFlux", &type, dim, &modelFlux); 
  if (!strncmp (strTemp, "COMP", 4)) modelType = OBIT_SkyModel_Comps;
  else if (!strncmp (strTemp, "IMAG", 3)) modelType = OBIT_SkyModel_Image;
  else modelType = OBIT_SkyModel_Comps;
  /* Is a model given in the parameters? */
  if (modelFlux!=0.0)  modelType = OBIT_SkyModel_Point;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "ModelType", OBIT_long, dim, &modelType);

  /* replace data with model? */
  ObitInfoListGetP(myInput, "Opcode", &type, dim, (gpointer)&opcode);
  if ((opcode!=NULL) && !strncmp (opcode, "MODL", 4)) replace = TRUE;
  else replace = FALSE;
  dim[0] = dim[1] = 1;
  ObitInfoListAlwaysPut (myInput, "REPLACE", OBIT_bool, dim, &replace);

  /* if Factor==0.0 replace with 1.0 */
  Factor = 1.0;
  ObitInfoListGetTest(myInput, "Factor",  &type, dim, &Factor);
  if (Factor==0.0) Factor = 1.0;
  dim[0] = 1;dim[1] = 1;
  ObitInfoListAlwaysPut (myInput, "Factor", OBIT_float, dim, &Factor);

  /* Initialize Threading */
  ObitThreadInit (myInput);
 
  /* Make sure doCalSelect set properly */
  doCalSelect = TRUE;
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doCalib",  &type, dim, &doCalib);
  doCalSelect = doCalSelect || (doCalib>0);
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);

  /* Default data types */
  ObitInfoListGetP (myInput, "in3DType", &type, dim, (gpointer)&strTemp);
  /* if not, use DataType */
  if ((strTemp==NULL) || (!strncmp(strTemp, "    ", 4))) {
    ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&strTemp);
  }
  ObitInfoListAlwaysPut (myInput, "in3DataType", type, dim, strTemp);

  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&strTemp);
  /* if not, use DataType */
  if ((strTemp==NULL) || (!strncmp(strTemp, "    ", 4))) {
    ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&strTemp);
  }
  ObitInfoListAlwaysPut (myInput, "outDataType", type, dim, strTemp);

  /* Output name -
   FITS - defaults to "UVPolCor"+inFile */
  ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&strTemp);
  /* if not use "UVPolCor"+inFile */
  if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12))) {
    ObitInfoListGetP (myInput, "inFile", &type, dim, (gpointer)&strTemp);
    strncpy (outFile, "UVPolCor", 9);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) {outFile[i+9] = strTemp[i];} outFile[i+9] = 0;
    ObitTrimTrail(outFile);  /* remove trailing blanks */
    dim[0] = strlen(outFile);
    ObitInfoListAlwaysPut (myInput, "outFile", OBIT_string, dim, outFile);
  }

  /* AIPS */
  /* outName given? */
  ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&strTemp);
  /* if not, use inName */
  if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12))) {
    ObitInfoListGetP (myInput, "inName", &type, dim, (gpointer)&strTemp);
    for (i=0; i<12; i++) {Aname[i] = ' ';}  Aname[i] = 0;
    for (i=0; i<MIN(12,dim[0]); i++) Aname[i] = strTemp[i];
    /* Save any defaulting on myInput */
    dim[0] = 12;
    ObitInfoListAlwaysPut (myInput, "outName", OBIT_string, dim, Aname);
  }

  /* output AIPS class,  Default is "UVPoCo"*/
  ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp);
  if ((strTemp==NULL) || (!strncmp(strTemp, "      ", 6))) {
    strncpy (Aclass, "UVPoCo", 7);
    dim[0] = 6;
    ObitInfoListAlwaysPut (myInput, "outClass", OBIT_string, dim, Aclass);
  }

  /* output AIPS disk - default is inDisk */
  disk = 0;
  ObitInfoListGetTest(myInput, "outDisk", &type, dim, &disk);
  if (disk<=0) {
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    dim[0] = 1;
    ObitInfoListAlwaysPut (myInput, "outDisk", OBIT_oint, dim, &disk);
  }

  /* output AIPS sequence - default is new (0) */
  seq = 0;
  if (!ObitInfoListGetTest(myInput, "outSeq", &type, dim, &seq)) {
    dim[0] = 1;
    ObitInfoListAlwaysPut (myInput, "outSeq", OBIT_oint, dim, &seq);
  }

  /*  Output file not expected to exist */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListAlwaysPut (myInput, "outExist", OBIT_bool, dim, &btemp);

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
  olong         nvis;
  oint          doCalib;
  gint32        dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean      doCalSelect;
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "Stokes", "timeRange", "BChan", "EChan",  "BIF", "EIF", "subA",
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer",
    "doPol", "PDVer", "keepLin",
    "Smooth", "Antennas",  "Sources",  "souCode", "Qual", "Alpha", 
     NULL};
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
 
  /* Ensure inData fully instantiated and OK */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Set number of vis per IO */
  nvis = 1000;  /* How many vis per I/O? */
  nvis =  ObitUVDescSetNVis (inData->myDesc, myInput, nvis);
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO", OBIT_long, dim,  &nvis);

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Get input sky model                                                   */
/*  Does CC table merge if requested.                                     */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      uvdata    Input uv data                                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*      Sky Model to be used                                              */
/*----------------------------------------------------------------------- */
ObitSkyModelVMBeamMF* getInputSkyModel (ObitInfoList *myInput, ObitUV *inData, 
				        ObitErr *err)
{
  ObitSkyModelVMBeamMF *skyModel=NULL;
  ObitImageMosaic *mosaic=NULL;
  ObitImage    **RXBeam=NULL,   **LYBeam=NULL,   **RLBeam=NULL,   **LRBeam=NULL;
  ObitImage    **RXBeamIm=NULL, **LYBeamIm=NULL, **RLBeamIm=NULL, **LRBeamIm=NULL;
  ObitImage    **image=NULL;
  ObitInfoType type;
  ObitTableCC *inCC=NULL;
  gboolean     mrgCC=FALSE, doCmplx=FALSE, do3D=FALSE;
  oint         noParms, CCVer;
  olong        Aseq, disk, cno, i, nparm, nmaps, channel, numAntType;
  gchar        *Type, *strTemp, inFile[129], inRoot[129];
  gchar        Aname[13], Aclass[8], Aroot[8], *Atype = "MA";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ofloat       modelFlux, modelPos[2], *modelParm=NULL, *Diams=NULL;
  ofloat       modptflx,  modptxof, modptyof, modptypm[8];
  olong        inVer;
  gchar        name[101];
  gchar        *skyModelParms[] = {  /* skyModel parameters to save*/
    "prtLv", "noNeg",
    NULL};
  gchar        *dataParms[] = {  /* Control parameters */
    "CCVer",  "BComp",  "EComp",  "Flux", "PBCor", "antSize", "Factor", 
    "minFlux", "Mode", "ModelType", "REPLACE", "Stokes", 
    "BIF", "EIF", "BCHAN", "ECHAN",
    "MODPTFLX", "MODPTXOF", "MODPTYOF", "MODPTYPM", 
    NULL};
  gchar *routine = "getInputSkyModel";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return skyModel;
  g_assert (ObitInfoListIsA(myInput));

  /* Want Imag. images? */
  ObitInfoListGetTest(myInput, "doCmplx", &type, dim, &doCmplx); 

  /* Image model or model in parameters? */
  modelFlux = 0.0;
  ObitInfoListGetTest (myInput, "modelFlux", &type, dim, &modelFlux);
  modelPos[0] = modelPos[1] = 0.0;
  ObitInfoListGetTest (myInput, "modelPos", &type, dim, &modelPos);
  ObitInfoListGetP (myInput, "modelParm", &type, dim, (gpointer)&modelParm);
  nparm = dim[0];
  if (modelFlux!=0.0) {
    /* Model passed - switch sign of shift */
    modptflx = modelFlux;
    modptxof = -modelPos[0] / 3600.0;
    modptyof = -modelPos[1] / 3600.0;
    if (modelParm!=NULL) {
      for (i=0; i<nparm; i++) modptypm[i] = modelParm[i];
    } else {
      for (i=0; i<8; i++) modptypm[i] = 0.0;
    }
    dim[0] = dim[1] = 1;
    ObitInfoListAlwaysPut (myInput, "MODPTFLX", OBIT_float, dim, &modptflx);
    ObitInfoListAlwaysPut (myInput, "MODPTXOF", OBIT_float, dim, &modptxof);
    ObitInfoListAlwaysPut (myInput, "MODPTYOF", OBIT_float, dim, &modptyof);
    dim[0] = nparm;
    ObitInfoListAlwaysPut (myInput, "MODPTYPM", OBIT_float, dim, modptypm);

    /* Create Sky Model */
    skyModel = newObitSkyModelVMBeamMF ("Sky Model");

  } else {
    /* image or components model */
    
    /* How many fields? */
    nmaps = 1;
    ObitInfoListGetTest(myInput, "nmaps", &type, dim, &nmaps);

    /* Allocate Image array */
    image = g_malloc0(nmaps*sizeof(ObitImage));
    
    /* Create image mosaic */
    mosaic = newObitImageMosaic ("Mosaic", nmaps);
    
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
	strncpy (Aroot, strTemp, 7);
      } else { /* Didn't find */
	strncpy (Aroot, "NoClas", 7);
      }
      if (nmaps>1) Aroot[2] = 0;  /* Multiple facets */

      /* input AIPS sequence */
      ObitInfoListGet(myInput, "in2Seq", &type, dim, &Aseq, err);
      
      /* if ASeq==0 want highest existing sequence */
      if (Aseq<=0) {
	g_snprintf (Aclass, 7, "%s%4.4d",Aroot,1);
	Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	/* Save on myInput*/
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(myInput, "inSeq", OBIT_oint, dim, &Aseq);
      }
      
      /* Loop over fields */
      for (i=0; i<nmaps; i++) {
	g_snprintf (name, 100, "Input image %d",i+1);
	if (nmaps>1)  g_snprintf (Aclass, 7, "%s%4.4d",Aroot,i+1); /* Multiple facets */
	else strncpy (Aclass, Aroot, 7);

	/* Find catalog number */
	cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
	if (cno<0) Obit_log_error(err, OBIT_Error, "Failure looking up %s", name);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	
 	/* define object */
	image[i] = newObitImage(name);
	ObitImageSetAIPS(image[i], OBIT_IO_byPlane, disk, cno, AIPSuser,  blc, trc, err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);

	/* Attach Image */
	ObitImageMosaicSetImage (mosaic, i, image[i], err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
      } /* end loop over fields */
      
      /* get do3D from first image */
      do3D = image[0]->myDesc->do3D;
      
    } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
      /* input FITS file name */
      if (ObitInfoListGetP(myInput, "in2File", &type, dim, (gpointer)&strTemp)) {
	strncpy (inRoot, strTemp, 128);
      } else { 
	strncpy (inRoot, "No_Filename_Given", 128);
      }
      ObitTrimTrail(inRoot);  /* remove trailing blanks */
      
      /* input FITS disk */
      ObitInfoListGet(myInput, "in2Disk", &type, dim, &disk, err);
      
      /* Loop over fields */
      for (i=0; i<nmaps; i++) {
	/* Set file name */
	if (nmaps>1) g_snprintf (inFile, 128, "%s%4.4d",inRoot,i);
	else g_snprintf (inFile, 128, "%s",inRoot);

 	/* define object */
	g_snprintf (name, 100, "Input image %d",i+1);
	image[i] = newObitImage(name);
	ObitImageSetFITS(image[i], OBIT_IO_byPlane, disk, inFile, blc, trc, err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);

	/* Attach Image */
	ObitImageMosaicSetImage (mosaic, i, image[i], err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
      } /* end loop over fields */
      
      /* Save parameters */
      ObitInfoListCopyList (myInput, skyModel->info, skyModelParms);
      
    } else { /* Unknown type - barf and bail */
      Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		     pgmName, Type);
      return skyModel;
    }

    /* Merge CC tables? */
    ObitInfoListGetTest(myInput, "mrgCC", &type, dim, &mrgCC); 
    CCVer = 1;
    ObitInfoListGetTest(myInput, "CCVer", &type, dim, &CCVer); 
    if (mrgCC) {
      noParms = 0;
      for (i=0; i<nmaps; i++) {

	/* Open Image to get table List */
	ObitImageOpen (image[i], OBIT_IO_ReadWrite, err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	inVer = CCVer;
	inCC = newObitTableCCValue ("inCC", (ObitData*)image[i], &inVer, OBIT_IO_ReadOnly, 
				    noParms, err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	if (inCC) {  /* May not exist */
	  ObitTableCCUtilMerge (inCC, inCC, err);
	  if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	  inCC = ObitTableCCUnref(inCC);
	  ObitImageClose (image[i],  err);
	  if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	}
	image[i] = ObitImageUnref(image[i]);  /* deallocate image */
      }
      /* give messages */
      ObitErrLog(err); 
    } /* end merge CC tables */
    else {
      /* deallocate images */
      for (i=0; i<nmaps; i++) image[i] = ObitImageUnref(image[i]);
    }
    
    g_free(image);  /* Deallocate array */

  } /* End image or components model */
  
  /* Get Beam */
  getBeam (myInput, doCmplx, inData->myDesc->crval[inData->myDesc->jlocs], &numAntType,
	   &RXBeam,   &LYBeam,   &RLBeam,   &LRBeam, 
	   &RXBeamIm, &LYBeamIm, &RLBeamIm, &LRBeamIm, &Diams, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);

  /* Create Sky model */
  skyModel = ObitSkyModelVMBeamMFCreate("Sky Model", mosaic, inData, numAntType,
					RXBeam, LYBeam, RLBeam, LRBeam,
					RXBeamIm, LYBeamIm, RLBeamIm, LRBeamIm,
					Diams, err);

 /* Create sky model */
  if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);

  /* Get input parameters from myInput, copy to skyModel */
  ObitInfoListCopyList (myInput, skyModel->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, skyModel->name, skyModel);

  /* Save do3D */
  dim[0] = 1; dim[1] = 1;
  ObitInfoListAlwaysPut (skyModel->info, "do3D", OBIT_bool, dim, &do3D);
  
  /* If channel given, select by channel */
  channel = 0;
  ObitInfoListGetTest (myInput, "channel", &type, dim, &channel);
  if (channel>0) {
    dim[0] = dim[1] = dim[2] = dim[3] = 1;
    ObitInfoListAlwaysPut (skyModel->info, "BChan", OBIT_oint, dim, &channel);
    ObitInfoListAlwaysPut (skyModel->info, "EChan", OBIT_oint, dim, &channel);
  }
  
  return skyModel;
} /* end getInputSkyModel */
  
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
ObitUV* setOutputData (ObitInfoList *myInput, ObitUV* inData, ObitErr *err)
{
  ObitUV    *outUV = NULL;
  ObitInfoType type;
  olong     nvis;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar     *routine = "setOutputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Build basic input UV data Object */
  outUV = ObitUVFromFileInfo ("out", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

  /* Set buffer size */
  nvis = 1000; type = OBIT_long;
  ObitInfoListGetTest(inData->info, "nVisPIO", &type, dim, &nvis);
  ObitInfoListAlwaysPut (outUV->info, "nVisPIO",  type, dim,  &nvis);
    
  /* Clone from input */
  ObitUVClone (inData, outUV, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
  
  /* Ensure outUV fully instantiated and OK */
  ObitUVFullInstantiate (outUV, FALSE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

  ObitErrLog(err); /* Show messages */
  return outUV;
} /* end setOutputUV */

/*----------------------------------------------------------------------- */
/*  Write History for UVPolCor                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outData   ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void UVPoCoHistory (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		   ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", 
    "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "channel", "BIF", "EIF",   "Sources",  "Qual", 
    "doCalSelect", "doCalib", "gainUse", "doPol", "PDVer", "keepLin", 
    "flagVer", "doBand", "BPVer", "Smooth", 
    "in2File",  "in2Disk", "in2Name", "in2Class", "in2Seq", "doCmplx", 
    "in3Diam", "in3DType", "in3File",  "in3Disk", "in3Name", "in3Class", "in3Seq",
    "in4File",  "in4Disk", "in4Name", "in4Class", "in4Seq",
    "in5Diam", "in5DType", "in5File",  "in5Disk", "in5Name", "in5Class", "in5Seq",
    "in6File",  "in6Disk", "in6Name", "in6Class", "in6Seq",
    "nmaps", "CCVer", "BComp",  "EComp", "Flux",
    "outDType", "outFile",  "outDisk",  "outName", "outClass", "outSeq",
    "Cmethod", "Cmodel", "Factor",  "Opcode", 
    "modelFlux", "modelPos", "modelParm", "noNeg",
    "mrgCC", "PBCor", "antSize", "Alpha",
    "nThreads",
    NULL};
  gchar *routine = "UVPoCoHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));
  g_assert (ObitUVIsA(outData));

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
 
} /* end UVPoCoHistory  */

/*----------------------------------------------------------------------- */
/*  Get Beam images                                                       */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      doCmplx   If TRUE, also get imag. images                          */
/*   Output:                                                              */
/*      numAntType number of antenna types in RXpol... arrays             */
/*      RXpol     R/X pol image array                                     */
/*      LYpol     L/Y pol image array                                     */
/*      RLpol     RL/XY pol image array  or NULL if none                  */
/*      LRpol     LR/YX pol image array  or NULL if none                  */
/*      RXpolIm   R/X pol imag. image array or NULL if none               */
/*      LYpolIm   L/Y  pol imag. image array or NULL if none              */
/*      RLpolIm   RL/XY pol imag. image array or NULL if none             */
/*      LRpolIm   LR/YX pol imag. image array or NULL if none             */
/*      Diams     Diameters of antenna types (m)                          */
/*----------------------------------------------------------------------- */
void getBeam (ObitInfoList *myInput, gboolean doCmplx, ofloat Stokes0,
	      olong *numAntType, 
	      ObitImage ***RXpol, ObitImage ***LYpol, 
	      ObitImage ***RLpol, ObitImage ***LRpol, 
	      ObitImage ***RXpolIm, ObitImage ***LYpolIm, 
	      ObitImage ***RLpolIm, ObitImage ***LRpolIm, 
	      ofloat **Diams, ObitErr *err)
{
  ObitInfoType type;
  ObitImage *RXpol1=NULL, *LYpol1=NULL, *RLpol1=NULL, *LRpol1=NULL;
  ObitImage *RXpol2=NULL, *LYpol2=NULL, *RLpol2=NULL, *LRpol2=NULL;
  ObitImage *RXpolIm1=NULL, *LYpolIm1=NULL, *RLpolIm1=NULL, *LRpolIm1=NULL;
  ObitImage *RXpolIm2=NULL, *LYpolIm2=NULL, *RLpolIm2=NULL, *LRpolIm2=NULL;
  gchar     Aclass[20], *strTemp, *Type, inFile[129];
  ofloat    Diam1=0.0, Diam2=0.0;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean  doRRLL;
  olong     i;
  gchar *routine = "getBeam";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  doRRLL = (Stokes0>-4.5);

  /* Antenna diameters */
  Diam1 = 0.0; Diam2 = 0.0; 
  ObitInfoListGetTest(myInput, "in3Diam", &type, dim, &Diam1);
  ObitInfoListGetTest(myInput, "in5Diam", &type, dim, &Diam2);

  /* First type Beam */
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "in3DType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || ((Type[0]==' ')&&(Type[1]==' ')&&(Type[2]==' ')))
    ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  dim[0] = 4; dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "in3DataType", OBIT_string, dim, Type);
  ObitInfoListAlwaysPut (myInput, "in4DataType", OBIT_string, dim, Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
    /* AIPS Class */
    if (doRRLL) strncpy (Aclass, "RR    ", 7);
    else        strncpy (Aclass, "XX    ", 7);
    ObitInfoListGetTest(myInput, "in3Class", &type, dim, Aclass);
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    if (ObitInfoListGetP(myInput, "in3File", &type, dim, (gpointer)&strTemp)) {
      strncpy (inFile, strTemp, 128);
    } else { 
      strncpy (inFile, "No_Filename_Given", 128);
    }
    ObitTrimTrail(inFile);  /* remove trailing blanks */
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return;
  }
  
  /* Stokes R/X */
  if (doRRLL) {inFile[0]='R';inFile[1]='R';}
  else        {inFile[0]='X';inFile[1]='X';}
  dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "in3File", OBIT_string, dim, inFile);
  if (doRRLL) strncpy (Aclass, "RR    ", 7);
  else        strncpy (Aclass, "XX    ", 7);
  Aclass[6] = 0;
  dim[0] = 6; dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "in3Class", OBIT_string, dim, Aclass);
  RXpol1 = ObitImageFromFileInfo ("in3", myInput, err);
   /* Set name */
  if (RXpol1) {
    if (RXpol1->name) g_free(RXpol1->name);
    RXpol1->name = g_strdup("RXBeam1");
  }
  if (err->error) Obit_traceback_msg (err, routine, "myInput");
  ObitErrLog(err); /* Show messages */

  /* Stokes L/Y */
  if (doRRLL) {inFile[0]='L';inFile[1]='L';}
  else        {inFile[0]='Y';inFile[1]='Y';}
  dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "in3File", OBIT_string, dim, inFile);
  if (doRRLL) strncpy (Aclass, "LL    ", 7);
  else        strncpy (Aclass, "YY    ", 7);
  dim[0] = 6; dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "in3Class", OBIT_string, dim, Aclass);
  LYpol1 = ObitImageFromFileInfo ("in3", myInput, err);
  /* Set name */
  if (LYpol1) {
    if ((LYpol1)->name) g_free((LYpol1)->name);
    (LYpol1)->name = g_strdup("LYBeam1");
  }
  if (err->error) Obit_traceback_msg (err, routine, "myInput");
  ObitErrLog(err); /* Show messages */

  /* Stokes  RL/XY */
  if (doRRLL) {inFile[0]='R';inFile[1]='L';}
  else        {inFile[0]='X';inFile[1]='Y';}
  dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "in3File", OBIT_string, dim, inFile);
  Aclass[0]='Q';
  if (doRRLL) strncpy (Aclass, "RL    ", 7);
  else        strncpy (Aclass, "XY    ", 7);
  dim[0] = 6; dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "in3Class", OBIT_string, dim, Aclass);
  RLpol1 = ObitImageFromFileInfo ("in3", myInput, err);
  /* Set name */
  if (RLpol1) {
    if ((RLpol1)->name) g_free((RLpol1)->name);
    (RLpol1)->name = g_strdup("RLBeam");
  }
  ObitErrClear(err);  /* Suppress failure messages */

  /* Stokes U if present  - debug LR/YX */
  inFile[0]='U';
  if (doRRLL) {inFile[0]='L';inFile[1]='R';}
  else        {inFile[0]='Y';inFile[1]='X';}
  dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "in3File", OBIT_string, dim, inFile);
  Aclass[0]='U';
  if (doRRLL) strncpy (Aclass, "LR    ", 7);
  else        strncpy (Aclass, "YX    ", 7);
  dim[0] = 6; dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "in3Class", OBIT_string, dim, Aclass);
  LRpol1 = ObitImageFromFileInfo ("in3", myInput, err);
  /* Set name */
  if (LRpol1) {
    if ((LRpol1)->name) g_free((LRpol1)->name);
    (LRpol1)->name = g_strdup("LRBeam");
  }
  ObitErrClear(err);  /* Suppress failure messages */

  /* Also imag.? */
  if (doCmplx) {
    ObitInfoListGetP (myInput, "in4Type", &type, dim, (gpointer)&Type);
    if ((Type==NULL) || ((Type[0]==' ')&&(Type[1]==' ')&&(Type[2]==' ')))
      ObitInfoListGetP (myInput, "in3DType", &type, dim, (gpointer)&Type);
    /* Get base parts of name */
    if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
      /* AIPS Class */
      strncpy (Aclass, "I     ", 7);
      ObitInfoListGetTest(myInput, "in4Class", &type, dim, Aclass);
    } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
      /* input FITS file name */
      if (ObitInfoListGetP(myInput, "in4File", &type, dim, (gpointer)&strTemp)) {
	strncpy (inFile, strTemp, 128);
      } else { 
	strncpy (inFile, "No_Filename_Given", 128);
      }
      ObitTrimTrail(inFile);  /* remove trailing blanks */
    } 

    /* Stokes R/X */
    if (doRRLL) {inFile[0]='R';inFile[1]='R';}
    else        {inFile[0]='X';inFile[1]='X';}
    dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in4File", OBIT_string, dim, inFile);
    if (doRRLL) strncpy (Aclass, "RR    ", 7);
    else        strncpy (Aclass, "XX    ", 7);
    dim[0] = 6; dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in4Class", OBIT_string, dim, Aclass);
    RXpolIm1 = ObitImageFromFileInfo ("in4", myInput, err);
    /* Set name */
    if (RXpolIm1) {
      if ((RXpolIm1)->name) g_free((RXpolIm1)->name);
      (RXpolIm1)->name = g_strdup("RXBeam1 imag.");
    }
    if (err->error) Obit_traceback_msg (err, routine, "myInput");
    ObitErrLog(err); /* Show messages */

    /* Stokes L/Y */
    if (doRRLL) {inFile[0]='L';inFile[1]='L';}
    else        {inFile[0]='Y';inFile[1]='Y';}
    dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in4File", OBIT_string, dim, inFile);
    if (doRRLL) strncpy (Aclass, "LL    ", 7);
    else        strncpy (Aclass, "YY    ", 7);
    dim[0] = 6; dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in4Class", OBIT_string, dim, Aclass);
    LYpolIm1 = ObitImageFromFileInfo ("in4", myInput, err);
    /* Set name */
    if (LYpolIm1) {
      if ((LYpolIm1)->name) g_free((LYpolIm1)->name);
      (LYpolIm1)->name = g_strdup("LYBeam1 imag.");
    }
    if (err->error) Obit_traceback_msg (err, routine, "myInput");
    ObitErrLog(err); /* Show messages */

    /* Stokes RL/XY if present */
    if (RLpol1) {
      if (doRRLL) {inFile[0]='R';inFile[1]='L';}
      else        {inFile[0]='X';inFile[1]='Y';}
      dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
      ObitInfoListAlwaysPut (myInput, "in4File", OBIT_string, dim, inFile);
      if (doRRLL) strncpy (Aclass, "RL    ", 7);
      else        strncpy (Aclass, "XY    ", 7);
      dim[0] = 6; dim[1] = dim[2] = 1;
      ObitInfoListAlwaysPut (myInput, "in4Class", OBIT_string, dim, Aclass);
      RLpolIm1 = ObitImageFromFileInfo ("in4", myInput, err);
      /* Set name */
      if (RLpolIm1) {
	if (RLpolIm1->name) g_free(RLpolIm1->name);
	RLpolIm1->name = g_strdup("RLBeam1 imag.");
      }
      if (err->error) Obit_traceback_msg (err, routine, "myInput");
      ObitErrClear(err);  /* Suppress failure messages */
    } /* end if RLPol */

    /* Stokes LR/YX if present */
    if (LRpol1) {
      if (doRRLL) {inFile[0]='L';inFile[1]='R';}
      else        {inFile[0]='Y';inFile[1]='X';}
      dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
      ObitInfoListAlwaysPut (myInput, "in4File", OBIT_string, dim, inFile);
      if (doRRLL) strncpy (Aclass, "LR    ", 7);
      else        strncpy (Aclass, "YX    ", 7);
      dim[0] = 6; dim[1] = dim[2] = 1;
      ObitInfoListAlwaysPut (myInput, "in4Class", OBIT_string, dim, Aclass);
      LRpolIm1 = ObitImageFromFileInfo ("in4", myInput, err);
      /* Set name */
      if (LRpolIm1) {
	if ((LRpolIm1)->name) g_free((LRpolIm1)->name);
	(LRpolIm1)->name = g_strdup("LRBeam1 imag.");
      }
      if (err->error) Obit_traceback_msg (err, routine, "myInput");
      ObitErrClear(err);  /* Suppress failure messages */
    } /* end if LRpol1 */

  } /* End also imag. */
  /****************** Second type Beam ************************/
  if (Diam2>0.0) {
    /* File type - could be either AIPS or FITS */
    ObitInfoListGetP (myInput, "in5DType", &type, dim, (gpointer)&Type);
    if ((Type==NULL) || ((Type[0]==' ')&&(Type[1]==' ')&&(Type[2]==' ')))
      ObitInfoListGetP (myInput, "in3DType", &type, dim, (gpointer)&Type);
    dim[0] = 4; dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in5DataType", OBIT_string, dim, Type);
    ObitInfoListAlwaysPut (myInput, "in6DataType", OBIT_string, dim, Type);
    if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
      /* AIPS Class */
      if (doRRLL) strncpy (Aclass, "RR    ", 7);
      else        strncpy (Aclass, "XX    ", 7);
      ObitInfoListGetTest(myInput, "in3Class", &type, dim, Aclass);
    } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
      /* input FITS file name */
      if (ObitInfoListGetP(myInput, "in5File", &type, dim, (gpointer)&strTemp)) {
	strncpy (inFile, strTemp, 128);
      } else { 
	strncpy (inFile, "No_Filename_Given", 128);
      }
      ObitTrimTrail(inFile);  /* remove trailing blanks */
    } else { /* Unknown type - barf and bail */
      Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		     pgmName, Type);
      return;
    }
    
    /* Stokes R/X */
    if (doRRLL) {inFile[0]='R';inFile[1]='R';}
    else        {inFile[0]='X';inFile[1]='X';}
    dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in5File", OBIT_string, dim, inFile);
    if (doRRLL) strncpy (Aclass, "RR    ", 7);
    else        strncpy (Aclass, "XX    ", 7);
    Aclass[6] = 0;
    dim[0] = 6; dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in5Class", OBIT_string, dim, Aclass);
    RXpol2 = ObitImageFromFileInfo ("in5", myInput, err);
    /* Set name */
    if (RXpol2) {
      if (RXpol2->name) g_free(RXpol2->name);
      RXpol2->name = g_strdup("RXBeam2");
    }
    if (err->error) Obit_traceback_msg (err, routine, "myInput");
    ObitErrLog(err); /* Show messages */
    
    /* Stokes L/Y */
    if (doRRLL) {inFile[0]='L';inFile[1]='L';}
    else        {inFile[0]='Y';inFile[1]='Y';}
    dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in5File", OBIT_string, dim, inFile);
    if (doRRLL) strncpy (Aclass, "LL    ", 7);
    else        strncpy (Aclass, "YY    ", 7);
    dim[0] = 6; dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in5Class", OBIT_string, dim, Aclass);
    LYpol2 = ObitImageFromFileInfo ("in5", myInput, err);
    /* Set name */
    if (LYpol2) {
      if ((LYpol2)->name) g_free((LYpol2)->name);
      (LYpol2)->name = g_strdup("LYBeam2");
    }
    if (err->error) Obit_traceback_msg (err, routine, "myInput");
    ObitErrLog(err); /* Show messages */
    
    /* Stokes  RL/XY */
    if (doRRLL) strncpy (Aclass, "RL    ", 7);
    else        strncpy (Aclass, "XY    ", 7);
    dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in5File", OBIT_string, dim, inFile);
    Aclass[0]='Q';
    if (doRRLL) {Aclass[0]='R';Aclass[1]='L';}
    else        {Aclass[0]='X';Aclass[1]='Y';}
    dim[0] = 6; dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in5Class", OBIT_string, dim, Aclass);
    RLpol2 = ObitImageFromFileInfo ("in5", myInput, err);
    /* Set name */
    if (RLpol2) {
      if ((RLpol2)->name) g_free((RLpol2)->name);
      (RLpol2)->name = g_strdup("RLBeam2");
    }
    ObitErrClear(err);  /* Suppress failure messages */
    
    /* Stokes U if present  - debug LR/YX */
    inFile[0]='U';
    if (doRRLL) {inFile[0]='L';inFile[1]='R';}
    else        {inFile[0]='Y';inFile[1]='X';}
    dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in5File", OBIT_string, dim, inFile);
    Aclass[0]='U';
    if (doRRLL) strncpy (Aclass, "LR    ", 7);
    else        strncpy (Aclass, "YX    ", 7);
    dim[0] = 6; dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "in5Class", OBIT_string, dim, Aclass);
    LRpol2 = ObitImageFromFileInfo ("in5", myInput, err);
    /* Set name */
    if (LRpol2) {
      if ((LRpol2)->name) g_free((LRpol2)->name);
      (LRpol2)->name = g_strdup("LRBeam2");
    }
    ObitErrClear(err);  /* Suppress failure messages */
    
    /* Also imag.? */
    if (doCmplx) {
      ObitInfoListGetP (myInput, "in6DType", &type, dim, (gpointer)&Type);
      if ((Type==NULL) || ((Type[0]==' ')&&(Type[1]==' ')&&(Type[2]==' ')))
	ObitInfoListGetP (myInput, "in5DType", &type, dim, (gpointer)&Type);
      /* Get base parts of name */
      if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
	/* AIPS Class */
	strncpy (Aclass, "I     ", 7);
	ObitInfoListGetTest(myInput, "in6Class", &type, dim, Aclass);
      } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
	/* input FITS file name */
	if (ObitInfoListGetP(myInput, "in6File", &type, dim, (gpointer)&strTemp)) {
	  strncpy (inFile, strTemp, 128);
	} else { 
	  strncpy (inFile, "No_Filename_Given", 128);
	}
	ObitTrimTrail(inFile);  /* remove trailing blanks */
      } 
      
      /* Stokes R/X */
      if (doRRLL) {inFile[0]='R';inFile[1]='R';}
      else        {inFile[0]='X';inFile[1]='X';}
      dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
      ObitInfoListAlwaysPut (myInput, "in6File", OBIT_string, dim, inFile);
      if (doRRLL) strncpy (Aclass, "RR    ", 7);
      else        strncpy (Aclass, "XX    ", 7);
      dim[0] = 6; dim[1] = dim[2] = 1;
      ObitInfoListAlwaysPut (myInput, "in6Class", OBIT_string, dim, Aclass);
      RXpolIm2 = ObitImageFromFileInfo ("in6", myInput, err);
      /* Set name */
      if (RXpolIm2) {
	if ((RXpolIm2)->name) g_free((RXpolIm2)->name);
	(RXpolIm2)->name = g_strdup("RXBeam2 imag.");
      }
      if (err->error) Obit_traceback_msg (err, routine, "myInput");
      ObitErrLog(err); /* Show messages */
      
      /* Stokes L/Y */
      if (doRRLL) {inFile[0]='L';inFile[1]='L';}
      else        {inFile[0]='Y';inFile[1]='Y';}
      dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
      ObitInfoListAlwaysPut (myInput, "in6File", OBIT_string, dim, inFile);
      if (doRRLL) strncpy (Aclass, "LL    ", 7);
      else        strncpy (Aclass, "YY    ", 7);
      dim[0] = 6; dim[1] = dim[2] = 1;
      ObitInfoListAlwaysPut (myInput, "in6Class", OBIT_string, dim, Aclass);
      LYpolIm2 = ObitImageFromFileInfo ("in6", myInput, err);
      /* Set name */
      if (LYpolIm2) {
	if ((LYpolIm2)->name) g_free((LYpolIm2)->name);
	(LYpolIm2)->name = g_strdup("LYBeam2 imag.");
      }
      if (err->error) Obit_traceback_msg (err, routine, "myInput");
      ObitErrLog(err); /* Show messages */
      
      /* Stokes RL/XY if present */
      if (RLpol2) {
	if (doRRLL) {inFile[0]='R';inFile[1]='L';}
	else        {inFile[0]='X';inFile[1]='Y';}
	dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
	ObitInfoListAlwaysPut (myInput, "in6File", OBIT_string, dim, inFile);
	if (doRRLL) strncpy (Aclass, "RL    ", 7);
	else        strncpy (Aclass, "XY    ", 7);
	dim[0] = 6; dim[1] = dim[2] = 1;
	ObitInfoListAlwaysPut (myInput, "in6Class", OBIT_string, dim, Aclass);
	RLpolIm2 = ObitImageFromFileInfo ("in6", myInput, err);
	/* Set name */
	if (RLpolIm2) {
	  if (RLpolIm2->name) g_free(RLpolIm2->name);
	  RLpolIm2->name = g_strdup("RLBeam2 imag.");
	}
	if (err->error) Obit_traceback_msg (err, routine, "myInput");
	ObitErrClear(err);  /* Suppress failure messages */
      } /* end if RLPol */
      
      /* Stokes LR/YX if present */
      if (LRpol2) {
	if (doRRLL) {inFile[0]='L';inFile[1]='R';}
	else        {inFile[0]='Y';inFile[1]='X';}
	dim[0] = strlen(inFile); dim[1] = dim[2] = 1;
	ObitInfoListAlwaysPut (myInput, "in6File", OBIT_string, dim, inFile);
	if (doRRLL) strncpy (Aclass, "LR    ", 7);
	else        strncpy (Aclass, "YX    ", 7);
	dim[0] = 6; dim[1] = dim[2] = 1;
	ObitInfoListAlwaysPut (myInput, "in6Class", OBIT_string, dim, Aclass);
	LRpolIm2 = ObitImageFromFileInfo ("in6", myInput, err);
	/* Set name */
	if (LRpolIm2) {
	  if ((LRpolIm2)->name) g_free((LRpolIm2)->name);
	  (LRpolIm2)->name = g_strdup("LRBeam2 imag.");
	}
	if (err->error) Obit_traceback_msg (err, routine, "myInput");
	ObitErrClear(err);  /* Suppress failure messages */
      } /* end if LRpol2 */
    } /* End also imag. */
  }  /* End of second antenna type */

  /* How many antenna types */
  if (Diam2>0.0) *numAntType = 2;
  else *numAntType = 1;
  /* Create output */
  *Diams = g_malloc(*numAntType * sizeof(ofloat));
  *RXpol = g_malloc(*numAntType * sizeof(ObitImage*));
  *LYpol = g_malloc(*numAntType * sizeof(ObitImage*));
  *RLpol = g_malloc(*numAntType * sizeof(ObitImage*));
  *LRpol = g_malloc(*numAntType * sizeof(ObitImage*));
  *RXpolIm = g_malloc(*numAntType * sizeof(ObitImage*));
  *LYpolIm = g_malloc(*numAntType * sizeof(ObitImage*));
  *RLpolIm = g_malloc(*numAntType * sizeof(ObitImage*));
  *LRpolIm = g_malloc(*numAntType * sizeof(ObitImage*));
  for (i=0; i<*numAntType; i++) {
    (*Diams)[i]   = 0.0;
    (*RXpol)[i]   = NULL; (*LYpol)[i]   = NULL; (*RLpol)[i]   = NULL; (*LRpol)[i]   = NULL;
    (*RXpolIm)[i] = NULL; (*LYpolIm)[i] = NULL; (*RLpolIm)[i] = NULL; (*LRpolIm)[i] = NULL;
  }
  /* Fill in arrays*/
  (*Diams)[0]   = Diam1;
  (*RXpol)[0]   = RXpol1;   (*LYpol)[0]   = LYpol1;   (*RLpol)[0]   = RLpol1;   (*LRpol)[0]   = LRpol1;
  (*RXpolIm)[0] = RXpolIm1; (*LYpolIm)[0] = LYpolIm1; (*RLpolIm)[0] = RLpolIm1; (*LRpolIm)[0] = LRpolIm1;
  if (*numAntType==2) {
    (*Diams)[1]   = Diam2;
    (*RXpol)[1]   = RXpol2;   (*LYpol)[1]   = LYpol2;   (*RLpol)[1]   = RLpol2;   (*LRpol)[1]   = LRpol2;
    (*RXpolIm)[1] = RXpolIm2; (*LYpolIm)[1] = LYpolIm2; (*RLpolIm)[1] = RLpolIm2; (*LRpolIm)[1] = LRpolIm2;
 }
  
} /* end getBeam */

