/* $Id: Convol.c,v 1.9 2007/07/13 20:32:34 bcotton Exp $  */
/* Convol Obit task convolve an image with another image or a model   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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
#include "ObitCArray.h"
#include "ObitFFT.h"
#include "ObitFeatherUtil.h"
#include "ObitConvUtil.h"
#include "ObitImageUtil.h"
#include "ObitHistory.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* ConvolIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void ConvolOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Digest inputs */
void digestInputs(ObitInfoList *myInput, ObitErr *err);
/* Get image from inputs */
void convolGetImage(ObitInfoList *myInput, ObitErr *err);
/* Get convolving Fn from inputs */
ObitFArray* convolGetConvFn(ObitInfoList *myInput, ObitErr *err);
/* Convol images together to outImage */
void doConvol (ObitInfoList *myInput, ObitImage *inImage, ObitFArray *convFn, 
	       ObitImage *outImage, ObitErr *err);
/* Write History */
void doHistory (ObitInfoList *myInput, ObitImage *inImage, 
		ObitImage *outImage, ObitErr *err);

/* Program globals */
gchar *pgmName = "Convol";       /* Program name */
gchar *infile  = "Convol.in" ;   /* File with program inputs */
gchar *outfile = "Convol.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitImage *inImage;    /* Input image */
ObitImage *outImage;          /* output image */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/* Convol Obit program to convol together images of different resolution*/
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitErr      *err= NULL;
  ObitFArray   *convFn=NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = ConvolIn (argc, argv, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Get input image and output image */
  convolGetImage(myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get convolving function as ObitFArray */
  convFn = convolGetConvFn(myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Convol them together */
  doConvol (myInput, inImage, convFn, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* History */
  doHistory (myInput, inImage, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* cleanup */
  inImage = ObitImageUnref(inImage);
  outImage    = ObitImageUnref(outImage);
  myInput     = ObitInfoListUnref(myInput); 
  convFn = ObitFArrayUnref(convFn);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* ConvolIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "ConvolIn";

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
} /* end ConvolIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: Convol -input file -output ofile [args]\n");
    fprintf(stderr, "Convol Obit task = convol together up to 10 images\n");
    fprintf(stderr, "Images must be given in order of decreasing resolution\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def Convol.in\n");
    fprintf(stderr, "  -output output result file, def Convol.out\n");
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
/*   - Determines convolving beam and writes in useBeam                   */
/*   - Determines scaling and write in rescale                            */
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
  olong         i;
  ofloat  ftemp, Beam[3]={0.0,0.0,0.0}, useBeam[3]={0.0,0.0,0.0}, 
    oldBeam[3]={0.0,0.0,0.0};
  gboolean found;
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

  /* Deconvolve beam */
  ObitInfoListGet (myInput, "Beam", &type, dim, Beam, err);
  if (err->error) Obit_traceback_msg (err, routine, routine);

  /* Beam in angular units or pixels? */
  if (fabs(inImage->myDesc->cdelt[0])>0.0) {
    /* To degrees */
    Beam[0] /= 3600.0;
    Beam[1] /= 3600.0;
    /* Beam from Header */
    oldBeam[0] = inImage->myDesc->beamMaj;
    oldBeam[1] = inImage->myDesc->beamMin;
    oldBeam[2] = inImage->myDesc->beamPA;
    
    /* Deconvolve to get beam to convolve with */
    if ((oldBeam[0]>0.0) && (oldBeam[1]>0.0))
      ObitConvUtilDeconv (Beam[0], Beam[1], Beam[2], 
			  oldBeam[0], oldBeam[1], oldBeam[2],
			  &useBeam[0], &useBeam[1], &useBeam[2]);
    else {useBeam[0]=Beam[0]; useBeam[1]=Beam[1]; useBeam[2]=Beam[2];}
    /* No pixel spacing - use pixels */
    Obit_log_error(err, OBIT_InfoErr,"Convolving with Gaussian %f x %f @ %f",
		   useBeam[0]*3600.0,useBeam[1]*3600.0,useBeam[2] );
  } else {
    useBeam[0]=Beam[0]; useBeam[1]=Beam[1]; useBeam[2]=Beam[2];
    Obit_log_error(err, OBIT_InfoErr,"Convolving with Gaussian %f x %f pixels @ %f",
		   useBeam[0],useBeam[1],useBeam[2] );
  }
  ObitErrLog(err);
  
  /* Save */
  dim[0] = 3; dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "useBeam", OBIT_float, dim, useBeam);
 
  /* Unit rescaling factor */
  ftemp = 0.0;
  found = ObitInfoListGetTest(myInput, "Factor", &type, dim, &ftemp);
  if ((!found) || (ftemp==0.0)) {
    /* Work it out from change in beam size */
    if ((oldBeam[0]>0.0) && (oldBeam[1]>0.0) && (Beam[0]>0.0) && (Beam[1]>0.0))
      ftemp = (Beam[0]/oldBeam[0]) * (Beam[1]/oldBeam[1]);
      else  /* Oh bother - assume units per pixel */
	if (fabs(inImage->myDesc->cdelt[0])>0.0)
	  ftemp = (1.1331 * Beam[0]*Beam[1]) / 
	    (fabs(inImage->myDesc->cdelt[0])*fabs(inImage->myDesc->cdelt[1]));
	else ftemp = 1.0; /* No pixel spacing either */
  }

  /* Save */
  dim[0] = 1; dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "rescale", OBIT_float, dim, &ftemp);
  Obit_log_error(err, OBIT_InfoErr,"Scaling image by %f to preserve units", 
		 ftemp);
  ObitErrLog(err);

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
} /* end digestInputs */

/*----------------------------------------------------------------------- */
/*  Get images from myInput                                               */
/*  Values:                                                               */
/*      myImage   InfoList with inputs                                    */
/*      inImage   [out] ObitImage pointers                                */
/*                 Passed as global                                       */
/*      outImage  [out] Output ObitImage pointer                          */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void convolGetImage(ObitInfoList *myInput, ObitErr *err)
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
  gchar *routine = "convolGetImage";

  if (err->error) return;  /* existing error? */

  /* Get region from myInput */
  ObitInfoListGetTest(myInput, "BLC", &type, dim, blc); /* BLC */
  ObitInfoListGetTest(myInput, "TRC", &type, dim, trc); /* TRC */

  /* File type - could be either AIPS or FITS */
  ObitInfoListGet (myInput, "DataType", &type, dim, tname, err);
  if (err->error) Obit_traceback_msg (err, routine, routine);
  if (!strncmp (tname, "AIPS", 4)) { /* AIPS input */

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
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /* Make sure it's OK */
    ObitImageFullInstantiate (inImage, TRUE, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /* Output image */
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
    
  } else if (!strncmp (tname, "FITS", 4)) {  /* FITS input */

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
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /* Make sure it's OK */
    ObitImageFullInstantiate (inImage, TRUE, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
  
    /* Output image */ 
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
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Image type %s", 
                   pgmName, strTemp);
  }
  if (err->error) Obit_traceback_msg (err, routine, routine);

} /* end convolGetImage */

/*----------------------------------------------------------------------- */
/*  Get convolving function, either from the specified image of the model */
/*  Values:                                                               */
/*      myImage   InfoList with inputs                                    */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
ObitFArray* convolGetConvFn(ObitInfoList *myInput, ObitErr *err)
{
  ObitInfoType type;
  ObitFArray   *outArray=NULL;
  ObitImage    *outImage=NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong         plane[5] = {1,1,1,1,1};
  olong         j, k, Aseq, disk, cno;
  olong        ndim=2, naxis[2];
  ofloat       useBeam[3];
  gboolean     found;
  gchar        *strTemp=NULL, Opcode[5], inFile[128];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gchar        tname[101];
  gchar *routine = "convolGetConvFn";

  if (err->error) return outArray;  /* existing error? */

  /* Which operation wanted? GAUS, IMAG, DCON, DGAU */
  strcpy (Opcode, "    ");
  found = ObitInfoListGetTest(myInput, "Opcode", &type, dim, Opcode);
  if (!strcmp (Opcode, "    ")) strcpy (Opcode, "GAUS");
  dim[0] = 4; dim[1] = 1;
  ObitInfoListAlwaysPut(myInput, "Opcode", OBIT_string, dim, Opcode);

  if (!strcmp (Opcode, "GAUS") || !strcmp (Opcode, "DCON") ) {
    /* Use Gaussian */
    useBeam[0] = useBeam[1] = fabs(inImage->myDesc->cdelt[0]); 
    useBeam[2] = 0.0;
    ObitInfoListGetTest (myInput, "useBeam", &type, dim, useBeam);

    /* Create convolving Gaussian */
    outArray = ObitConvUtilGaus (inImage, useBeam);
    /* DEBUG 
    ObitImageUtilArray2Image ("ConvolDebug1.fits",1,outArray, err);*/
  } else if (!strcmp (Opcode, "IMAG")) {
    /* Convolve with an image */
    /* File type - could be either AIPS or FITS */
    ObitInfoListGet (myInput, "DataType", &type, dim, tname, err);
    if (err->error) Obit_traceback_val (err, routine, routine, outArray);
    if (!strncmp (tname, "AIPS", 4)) { /* AIPS input */
      
      /* input AIPS disk */
      ObitInfoListGet(myInput, "in2Disk", &type, dim, &disk, err);

      /* input AIPS name */
      for (k=0; k<12; k++) Aname[k] = ' '; Aname[k] = 0;
      ObitInfoListGet(myInput, "in2Name", &type, dim, Aname, err);

      /* input AIPS class */
      for (k=0; k<6; k++) Aclass[k] = ' '; Aclass[k] = 0;
      ObitInfoListGet(myInput, "in2Class", &type, dim, Aclass, err);

      /* input AIPS sequence */
      ObitInfoListGet(myInput, "in2Seq", &type, dim, &Aseq, err);
      if (err->error) Obit_traceback_val (err, routine, routine, outArray);
      /* if ASeq==0 want highest existing sequence */
      if (Aseq<=0) {
	Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
	if (err->error) Obit_traceback_val (err, routine, routine, outArray);
	/* Save on myInput*/
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(myInput, "in2Seq", OBIT_oint, dim, &Aseq);
      } 
      
      /* Find catalog number */
      cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
      if (err->error) Obit_traceback_val (err, routine, routine, outArray);
      
      /* Generate Object name from AIPS name */
      g_snprintf (tname, 100, "%s.%s:%d.%d", Aname, Aclass, Aseq, disk);
      outImage = newObitImage(tname);
      
      /* define image */
      ObitImageSetAIPS (outImage, OBIT_IO_byPlane, disk, cno, AIPSuser, 
			blc, trc, err);
      if (err->error) Obit_traceback_val (err, routine, routine, outArray);
       
    } else if (!strncmp (tname, "FITS", 4)) {  /* FITS input */
      
      /* input FITS file name */
      for (j=0; j<128; j++) inFile[j] = 0;
      g_snprintf (tname, 100, "inFile");
      ObitInfoListGet(myInput, tname, &type, dim, inFile, err);
      if (err->error) Obit_traceback_val (err, routine, routine, outArray);
      
      /* input FITS disk */
      ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
     if (err->error) Obit_traceback_val (err, routine, routine, outArray);
      
      /*  Object name from FITS name */
      outImage = newObitImage(inFile);
      
      /* define image */
      ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
      if (err->error) Obit_traceback_val (err, routine, routine, outArray);
       
    } else { /* Unknown type - barf and bail */
      Obit_log_error(err, OBIT_Error, "%s: Unknown Image type %s", 
		     pgmName, strTemp);
    }
    if (err->error) Obit_traceback_val (err, routine, routine, outArray);
 
    /* Make sure image OK */
    ObitImageFullInstantiate (outImage, TRUE, err);
    if (err->error) Obit_traceback_val (err, routine, routine, outArray);
    
    /* Read and extract image data*/
    naxis[0] = outImage->myDesc->inaxes[0];  
    naxis[1] = outImage->myDesc->inaxes[1]; 
    outArray = ObitFArrayCreate("Convolving function", ndim, naxis);
    ObitImageGetPlane (outImage, outArray->array, plane, err);
    outImage = ObitImageUnref(outImage);
    if (err->error) Obit_traceback_val (err, routine, routine, outArray);
   } else if (!strcmp (Opcode, "DGAU")) {
    Obit_log_error(err, OBIT_Error, "%s: Opcode DGAU not yet implemented", 
                   routine);
    return outArray;
  } /* End of Opcode type */

  return outArray;
} /* end convolGetConvFn */

/*----------------------------------------------------------------------- */
/*  Convolve images together to outImage                                  */
/*  Values:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inImage   input ObitImage                                         */
/*      convFn    Convolving Function                                     */
/*      outImage  Output ObitImage pointer                                */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void doConvol (ObitInfoList *myInput, ObitImage *inImage, ObitFArray *convFn, 
	       ObitImage *outImage, ObitErr *err)
{
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat     rescale, Beam[3];
  gboolean doDivide;
  gchar    Opcode[5];
  gchar *routine = "doConvol";

  if (err->error) return;  /* existing error? */

  /* Multiply or divide by FT of convolving function */
  ObitInfoListGetTest(myInput, "Opcode", &type, dim, Opcode);
  doDivide = (!strcmp (Opcode, "DCON")) || (!strcmp (Opcode, "DGAU"));

  /* Unit scaling */
  rescale = 1.0;
  ObitInfoListGetTest(myInput, "rescale", &type, dim, &rescale);

  /* Set output image resolution */
  if (ObitInfoListGetTest (myInput, "Beam", &type, dim, Beam)) {
    Beam[0] /= 3600.0;
    Beam[1] /= 3600.0;
  } else {  /* Use what's in input */
    Beam[0] = inImage->myDesc->beamMaj;
    Beam[1] = inImage->myDesc->beamMin;
    Beam[2] = inImage->myDesc->beamPA;
  }
  outImage->myDesc->beamMaj = Beam[0];
  outImage->myDesc->beamMin = Beam[1];
  outImage->myDesc->beamPA  = Beam[2];


  ObitConvUtilConv (inImage, convFn, doDivide, rescale, outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

} /* end doConvol */

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
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "in2File",  "in2Disk", "in2Name", "in2Class", "inS2eq",
    "outFile",  "outDisk", "outName", "outClass", "outSeq",
    "Beam", "BLC", "TRC", "Factor", "Opcode",
    NULL};
  gchar *routine = "doHistory";

  if (err->error) return;  /* existing error? */

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

