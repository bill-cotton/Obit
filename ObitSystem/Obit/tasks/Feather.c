/* $Id$  */
/* Feather Obit task - Feathers together images            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2015                                          */
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
#include "ObitFeatherUtil.h"
#include "ObitImageUtil.h"
#include "ObitHistory.h"

#define MAXINPUT 10  /* Maximum number of input images */
/* internal prototypes */
/* Get inputs */
ObitInfoList* FeatherIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void FeatherOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Get images from inputs */
void featherGetImage(ObitInfoList *myInput, olong *numImage, 
		     ObitErr *err);
/* Feather images together to outImage */
void doFeather (ObitInfoList *myInput, olong numImage, ObitImage *inImage[], 
		ObitImage *outImage, ObitErr *err);
/* Write History */
void doHistory (olong numImage, ObitImage *inImage[], ObitImage *outImage, 
		ObitErr *err);

/* Program globals */
gchar *pgmName = "Feather";       /* Program name */
gchar *infile  = "Feather.in" ;   /* File with program inputs */
gchar *outfile = "Feather.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitImage *inImage[MAXINPUT]; /* Input images */
ObitImage *outImage;          /* output image */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/* Feather Obit program to feather together images of different resolution*/
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitErr      *err= NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {IM_MAXDIM,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong        i, j, numImage, nplane, iplane, BChan=1, EChan=0;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = FeatherIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Get list of input images and output image */
  featherGetImage(myInput, &numImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Feather them together - loop over planes */
  nplane = inImage[0]->myDesc->inaxes[2];
  ObitInfoListGetTest(myInput, "BChan", &type, dim, &BChan);
  ObitInfoListGetTest(myInput, "EChan", &type, dim, &EChan);
  dim[0] = IM_MAXDIM;
  if (EChan<BChan) EChan = nplane;
  for (iplane=BChan; iplane<=EChan; iplane++) {
    /* Tell if more than one */
    if (nplane>1) {
      Obit_log_error(err, OBIT_InfoErr, "Plane %d", iplane);
      ObitErrLog(err);
      blc[2] = trc[2] = iplane;
      for (j=0; j<numImage; j++) {
	ObitInfoListAlwaysPut (inImage[j]->info, "BLC", OBIT_long, dim, blc); 
	ObitInfoListAlwaysPut (inImage[j]->info, "TRC", OBIT_long, dim, trc);
      }
      ObitInfoListAlwaysPut (outImage->info, "BLC", OBIT_long, dim, blc); 
      ObitInfoListAlwaysPut (outImage->info, "TRC", OBIT_long, dim, trc);
    }
    /* doit */
    doFeather (myInput, numImage, inImage, outImage, err);
    if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
  } /* end loop over planes */

  /* History */
  doHistory (numImage, inImage, outImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* cleanup */
  for (i=0; i<numImage; i++) {
    inImage[i] = ObitImageUnref(inImage[i]);
  }
  outImage    = ObitImageUnref(outImage);
  myInput     = ObitInfoListUnref(myInput); 
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* FeatherIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "FeatherIn";

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
      
     } else if (strcmp(arg, "-outDType") == 0) { /* Image type AIPS or FITS */
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
} /* end FeatherIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: Feather -input file -output ofile [args]\n");
    fprintf(stderr, "Feather Obit task = feather together up to 10 images\n");
    fprintf(stderr, "Images must be given in order of decreasing resolution\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def Feather.in\n");
    fprintf(stderr, "  -output output result file, def Feather.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input image\n");
    fprintf(stderr, "  -outDType AIPS or FITS type for output\n");
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
/*     DataType  Str [4]    input "AIPS" or "FITS" [def {"FITS"}]         */
/*     outDType  Str [4]    output "AIPS" or "FITS" [def {"FITS"}]        */
/*     in?File   Str [?]    input FITS image file name [no def]           */
/*     in?Name   Str [12]   input AIPS image name  [no def]               */
/*     in?Class  Str [6]    input AIPS image class  [no def]              */
/*     in?Seq    Int        input AIPS image sequence no  [no def]        */
/*     in?Disk   Int        input AIPS or FITS image disk no  [def 1]     */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  ofloat ftemp;
  oint   itemp, i;
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
  ObitInfoListPut (out, "outDType", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  strTemp = "undefn";
  /* Loop over potential inputs */
  for (i=0; i<MAXINPUT; i++) {

    /* input FITS file name */
    if (i==0) g_snprintf (tname, 50, "inFITS");
    else g_snprintf (tname, 50, "in%dFITS", i+1);
    dim[0] = strlen (strTemp); dim[1] = 1;
    ObitInfoListPut (out, tname, OBIT_string, dim, strTemp, err);
    if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
    /* input AIPS file name */
    if (i==0) g_snprintf (tname, 50, "inName");
    else g_snprintf (tname, 50, "in%dName", i+1);
    dim[0] = strlen (strTemp); dim[1] = 1;
    ObitInfoListPut (out, tname, OBIT_string, dim, strTemp, err);
    if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
    /* input AIPS file class */
    if (i==0) g_snprintf (tname, 50, "inClass");
    else g_snprintf (tname, 50, "in%dClass", i+1);
    dim[0] = strlen (strTemp); dim[1] = 1;
    ObitInfoListPut (out, "temp", OBIT_string, dim, strTemp, err);
    if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
    /* AIPS sequence */
    if (i==0) g_snprintf (tname, 50, "inSeq");
    else g_snprintf (tname, 50, "in%dSeq", i+1);
    dim[0] = 1;dim[1] = 1;
    ftemp = 0.0; 
    ObitInfoListPut (out, tname, OBIT_float, dim, &ftemp, err);
    if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
    /* AIPS or FITS disk number */
    if (i==0) g_snprintf (tname, 50, "inDisk");
    else g_snprintf (tname, 50, "in%dDisk", i+1);
    dim[0] = 1;dim[1] = 1;
    jtemp = 1; 
    ObitInfoListPut (out, tname, OBIT_long, dim, &jtemp, err);
    if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  }
  /* Default output */
  /*FITS file name */
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFITS", OBIT_string, dim, strTemp, err);
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
  jtemp = 0.0; 
  ObitInfoListPut (out, "outSeq", OBIT_long, dim, &jtemp, err);
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
/*  Get images from myInput                                               */
/*  Values:                                                               */
/*      myImage   InfoList with inputs                                    */
/*      numImage  [out] Number of images in inImage                       */
/*      inImage   [out] Array of ObitImage pointers                       */
/*                 Passed as global as there doesn't seem to be a way to  */
/*                 get it through a call sequence                         */
/*      outImage  [out] Output ObitImage pointer                          */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void featherGetImage(ObitInfoList *myInput, olong *numImage, 
		     ObitErr *err)
{
  ObitInfoType type;
  ObitImage    *tImage=NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong         i, j, k, Aseq, disk, cno;
  ofloat       ftemp;
  olong        number, nplane;
  gboolean     exist;
  gchar        *strTemp=NULL, inFile[128];
  gchar        Aname[13], Aclass[9], *Atype = "MA";
  gchar        tname[101], *Type, *FITS="FITS";
  gchar        *feathr = "Feathr";
  gchar *routine = "featherGetImage";

  if (err->error) return;  /* existing error? */

  /* Get number of input images */
  ObitInfoListGet(myInput, "numImage", &type, dim, &number, err);
  if (err->error) Obit_traceback_msg (err, routine, routine);
  *numImage = number;

  /* In range? */
  Obit_return_if_fail(((number>1)&&(number<=MAXINPUT)), err,
		      "%s: Number of input images,  %d, not in range [2,%d]", 
		      routine, number, MAXINPUT);
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGet (myInput, "DataType", &type, dim, tname, err);
  if (err->error) Obit_traceback_msg (err, routine, routine);
  if (!strncmp (tname, "AIPS", 4)) { /* AIPS input */

    /* Loop over input images */
    for (i=0; i<number; i++) {
      /* input AIPS disk */
      if (i==0) g_snprintf (tname, 50, "inDisk");
      else g_snprintf (tname, 100, "in%dDisk", i+1);
      ObitInfoListGet(myInput, tname, &type, dim, &disk, err);
      /* input AIPS name */
      if (i==0) g_snprintf (tname, 50, "inName");
      else g_snprintf (tname, 100, "in%dName", i+1);
      for (k=0; k<12; k++) Aname[k] = ' '; Aname[k] = 0;
      ObitInfoListGet(myInput, tname, &type, dim, Aname, err);
      /* input AIPS class */
      if (i==0) g_snprintf (tname, 50, "inClass");
      else g_snprintf (tname, 100, "in%dClass", i+1);
      for (k=0; k<6; k++) Aclass[k] = ' '; Aclass[k] = 0;
      ObitInfoListGet(myInput, tname, &type, dim, Aclass, err);
      /* input AIPS sequence */
      if (i==0) g_snprintf (tname, 50, "inSeq");
      else g_snprintf (tname, 100, "in%dSeq", i+1);
      ObitInfoListGet(myInput, tname, &type, dim, &Aseq, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
      /* if ASeq==0 want highest existing sequence */
      if (Aseq<=0) {
	Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
	if (err->error) Obit_traceback_msg (err, routine, "myInput");
	/* Save on myInput*/
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(myInput, tname, OBIT_oint, dim, &Aseq);
      } 
      
      /* Find catalog number */
      cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);

      /* Generate Object name from AIPS name */
      g_snprintf (tname, 100, "%s.%s:%d.%d", Aname, Aclass, Aseq, disk);
      inImage[i] = newObitImage(tname);
    
      /* define image */
      tImage = inImage[i];
      ObitImageSetAIPS (tImage, OBIT_IO_byPlane, disk, cno, AIPSuser, 
			blc, trc, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);

      /* Make sure it's OK */
      ObitImageFullInstantiate (inImage[i], TRUE, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
    } /* End loop over inputs */

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
    /* Loop over input images */
    for (i=0; i<number; i++) {
      /* input FITS file name */
      for (j=0; j<128; j++) inFile[j] = 0;
      if (i==0) g_snprintf (tname, 50, "inFile");
      else g_snprintf (tname, 100, "in%dFile", i+1);
      ObitInfoListGet(myInput, tname, &type, dim, inFile, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
      
      /* input FITS disk */
      if (i==0) g_snprintf (tname, 50, "inDisk");
      else g_snprintf (tname, 100, "in%dDisk", i+1);
      ObitInfoListGet(myInput, tname, &type, dim, &ftemp, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
      disk = ftemp + 0.5;
      
      /*  Object name from FITS name */
      inImage[i] = newObitImage(inFile);
      
      /* define image */
      tImage = inImage[i];
      ObitImageSetFITS (tImage, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);

      /* Make sure it's OK */
      ObitImageFullInstantiate (inImage[i], TRUE, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
    } /* end loop over inputs */
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Image type %s", 
                   pgmName, strTemp);
  }
  if (err->error) Obit_traceback_msg (err, routine, routine);

  /* Output image */
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4)))
    ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;

  if (!strncmp (Type, "AIPS", 4)) { 
    /* AIPS file */
    /* AIPS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* AIPS name */
    for (k=0; k<12; k++) Aname[k] = ' '; Aname[k] = 0;
    ObitInfoListGet(myInput, "outName", &type, dim, Aname, err);
    Aname[dim[0]] = 0;
    /* AIPS class */
    for (k=0; k<6; k++) Aclass[k] = feathr[k]; Aclass[k] = 0;
    ObitInfoListGet(myInput, "outClass", &type, dim, Aclass, err);
    Aclass[dim[0]] = 0;
    if (!strncmp (Aclass, "      ", 6)) 
      for (k=0; k<6; k++) Aclass[k] = feathr[k]; Aclass[k] = 0;
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
    
  } else if (!strncmp (Type, "FITS", 4)){ 

    /* FITS file */  
    /* Output image */ 
    /* FITS file name */
    for (i=0; i<128; i++) inFile[i] = 0;
    ObitInfoListGet(myInput, "outFile", &type, dim, inFile, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
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
  }
  if (err->error) Obit_traceback_msg (err, routine, routine);

  /* Check that images all have the same number of planes */
  nplane = inImage[0]->myDesc->inaxes[2];
  for (i=1; i<*numImage; i++) {
    Obit_return_if_fail((nplane==inImage[i]->myDesc->inaxes[2]), err,
			 "%s: Number of planes image %d = %d not %d", 
			 routine, i+1, inImage[i]->myDesc->inaxes[2], nplane);
   }

  /* Clone outImage from first input image */
  ObitImageClone (inImage[0], outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage[0]->name);
  
} /* end featherGetImage */

/*----------------------------------------------------------------------- */
/*   Feather images together to outImage                                  */
/*  Values:                                                               */
/*      numImage  Number of images in inImage                             */
/*      inImage   Array of ObitImage pointers                             */
/*      outImage  Output ObitImage pointer                                */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void doFeather (ObitInfoList *myInput, olong numImage, ObitImage *inImage[], 
		ObitImage *outImage, ObitErr *err)
{
  ObitInfoType type;
  gint32     dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong      i;
  olong      ndim=2, naxis[2], pos[2];
  ofloat     cmplx[2], peak, norm, *weights, tweight[100];
  ObitFFT    *FFTfor=NULL, *FFTrev=NULL;
  ObitImage  *padImage[MAXINPUT];
  ObitImage  *tmplImage=NULL, *tmpImage=NULL;
  ObitFArray *wtArray[MAXINPUT], *resultArray=NULL, *workArray2=NULL;
  ObitCArray *accArray=NULL, *workArray=NULL;
  /* ObitFArray *tArr;       DEBUG */
  /* gchar debugFile[200];   DEBUG */
  gchar *routine = "doFeather";

  if (err->error) return;  /* existing error? */

  /* Get weight array - default = all 1.0 */
  if (!ObitInfoListGetP (myInput, "weight", &type, dim, (gpointer)&weights)) {
    for (i=0; i<100; i++) tweight[i] = 1.0;
    weights = tweight;
  }

  /* Create FFTs */
  FFTfor = ObitFeatherUtilCreateFFT(inImage[0], OBIT_FFT_Forward);
  FFTrev = ObitFeatherUtilCreateFFT(inImage[0], OBIT_FFT_Reverse);

  /* Create padded images for FFT size */
  Obit_log_error(err, OBIT_InfoErr,"Pad/interpolate Images to same grid");
  ObitErrLog(err);

  for (i=0; i<numImage; i++) {
    Obit_log_error(err, OBIT_InfoErr,
		   "Pad Loop %d %s", i+1, inImage[i]->name);
    ObitErrLog(err);
    padImage[i] = ObitFeatherUtilCreatePadImage (FFTfor, inImage[0], err);
    if (err->error) Obit_traceback_msg (err, routine, inImage[i]->name);
    /* Merely pad first image */
    if (i==0) {
      ObitFeatherUtilPad (FFTfor, inImage[i], padImage[i], err);
    } else { /* interpolate and pad rest to same grid as first */
      ObitFeatherUtilInterpol (inImage[i], padImage[0], padImage[i], err);
    }
    if (err->error) Obit_traceback_msg (err, routine, inImage[i]->name);
  } /* end loop creating padded images */
  
  /* Create masks in FArrays, first get weights from restoring beams/resolution */
  Obit_log_error(err, OBIT_InfoErr,"Create weighting masks");
  ObitErrLog(err);

  for (i=0; i<numImage; i++) {
    wtArray[i] =  ObitFeatherUtilMakeBeamMask (padImage[i],  FFTfor, err);
    if (err->error) Obit_traceback_msg (err, routine, inImage[0]->name);
    /* DEBUG
       ObitImageOpen (padImage[i], OBIT_IO_ReadWrite, err);
       ObitImageRead (padImage[i], NULL, err);
       ObitImageClose (padImage[i],err);
       if (err->error) Obit_traceback_msg (err, routine, padImage[i]->name);
       g_snprintf (debugFile,200,"debugPad%d.fits",i+1);
       ObitImageUtilArray2Image (debugFile,0,padImage[i]->image, err); */
  } /* end loop creating weighting masks */

  /* derive weight masks from FT of beams, Weights are 1 with a Gaussian 
     hole in the middle representing the uv coverage of the next smallest 
     array/telescope */
  for (i=0; i<numImage; i++) {
    ObitFArrayFill (wtArray[i], 1.0);  /* 1 fill */
    /* If this is not the lowest resolution image, subtract next lowest */
    if (i<numImage-1) 
      ObitFArraySub (wtArray[i], wtArray[i+1], wtArray[i]);
    /* Scale by weights */
    ObitFArraySMul (wtArray[i], weights[i]);
    /* DEBUG 
       g_snprintf (debugFile,200,"debugWeight%d.fits",i+1);
       ObitImageUtilArray2Image (debugFile,0,wtArray[i], err);*/
  }  

  /* Make accumulation array and work array */
  Obit_log_error(err, OBIT_InfoErr,"Accumulate Weighted FFTs of images");
  ObitErrLog(err);
  accArray  = ObitFeatherUtilCreateFFTArray(FFTfor);
  cmplx[0] = 0.0; cmplx[1] = 0.0;
  ObitCArrayFill(accArray, cmplx);  /* Zero fill accumulation */
  workArray = ObitFeatherUtilCreateFFTArray(FFTfor);

  /* Loop accumulating images */
  for (i=0; i<numImage; i++) {
    ObitFeatherUtilAccumImage(FFTfor, padImage[i], wtArray[i], 
			      accArray, workArray, err);
    if (err->error) Obit_traceback_msg (err, routine, inImage[0]->name);
    /* DEBUG
       tArr = ObitCArrayMakeF(accArray);
       ObitCArrayAmp (workArray, tArr);
       g_snprintf (debugFile,200,"debugFT%d.fits",i+1);
       ObitImageUtilArray2Image (debugFile,0,tArr, err);
       tArr = ObitFArrayUnref(tArr); */
  } /* end loop accumulating images */

  /* DEBUG
     tArr = ObitCArrayMakeF(accArray);
     ObitCArrayReal (accArray, tArr);
     ObitImageUtilArray2Image ("FeatherDebug1.fits",0,tArr, err);
     tArr = ObitFArrayUnref(tArr); */

  /* FFT back to image domain */
  Obit_log_error(err, OBIT_InfoErr,"FFT back to image domain");
  ObitErrLog(err);
  
  naxis[0] = padImage[0]->myDesc->inaxes[0]; 
  naxis[1] = padImage[0]->myDesc->inaxes[1]; 
  resultArray = ObitFArrayCreate("Result Array", ndim, naxis);
  ObitFFTC2R(FFTrev, accArray, resultArray);
  ObitFArray2DCenter (resultArray);/* Swaparoonie */

  /* Get normalization by repeating but using the padded images 
     replaced by the beam */
  Obit_log_error(err, OBIT_InfoErr,"Get normalization using point models");
  ObitErrLog(err);
  
  cmplx[0] = 0.0; cmplx[1] = 0.0;
  ObitCArrayFill(accArray, cmplx);  /* Zero fill accumulation */
  
  /* Loop accumulating normalization images */
  for (i=0; i<numImage; i++) {
    /* replace array on padImage with model */
    ObitImageOpen (padImage[i], OBIT_IO_ReadWrite, err);
    if (err->error) Obit_traceback_msg (err, routine, padImage[i]->name);
    ObitFeatherUtilCreateModel(padImage[i], padImage[i]->image);
    ObitImageWrite (padImage[i], NULL, err);
    ObitImageClose (padImage[i],err);
    if (err->error) Obit_traceback_msg (err, routine, padImage[i]->name);
    ObitFeatherUtilAccumImage(FFTfor, padImage[i], wtArray[i], 
			      accArray, workArray, err);
    if (err->error) Obit_traceback_msg (err, routine, inImage[0]->name);
  } /* end loop accumulating normalization images */

  /* FFT normalization image back to image domain */
  workArray2 = newObitFArray("Scratch Array");
  ObitFArrayClone(resultArray, workArray2, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage[0]->name);
  ObitFFTC2R(FFTrev, accArray, workArray2);
  ObitFArray2DCenter (workArray2);/* Swaparoonie */
   
  /* Do normalization from peak in workArray2 */
  pos[0] = 0; pos[1] = 0;
  peak = ObitFArrayMax(workArray2, pos);
  Obit_log_error(err, OBIT_InfoErr,"peak in normalization image %f",peak);
  ObitErrLog(err);
  if (peak!=0.0)  norm = 1.0 / peak;
  else norm = 1.0;
  ObitFArraySMul(resultArray, norm);
  
  /* Extract to outImage from resultArray */
  ObitFeatherUtilSubImage (padImage[0], resultArray, outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage[0]->name);
  
  /* Delete scratch files */
  for (i=0; i<numImage; i++) {
    ObitImageZap(padImage[i], err);
    if (err->error) Obit_traceback_msg (err, routine, inImage[0]->name);
  } 
    
  /* cleanup */
  for (i=0; i<<numImage; i++) {
    wtArray[i] = ObitFArrayUnref(wtArray[i]);
  }
  tmplImage   = ObitImageUnref(tmplImage);
  tmpImage    = ObitImageUnref(tmpImage);
  accArray    = ObitCArrayUnref(accArray);
  workArray   = ObitCArrayUnref(workArray);
  resultArray = ObitFArrayUnref(resultArray);
  workArray2  = ObitFArrayUnref(workArray2);
  FFTfor      = ObitFFTUnref(FFTfor);
  FFTrev      = ObitFFTUnref(FFTrev);
} /* end doFeather */

/*----------------------------------------------------------------------- */
/*  Write history                                                         */
/*  Values:                                                               */
/*      numImage  Number of images in inImage                             */
/*      inImage   Array of ObitImage pointers                             */
/*      outImage  Output ObitImage pointer                                */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void doHistory (olong numImage, ObitImage *inImage[], ObitImage *outImage, 
		ObitErr *err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  olong         i;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "weight", "BChan", "EChan",
    NULL};
  gchar *routine = "doHistory";

  if (err->error) return;  /* existing error? */

  /* Open outImage to make sure History gets registered in header if created */
  ObitImageOpen (outImage, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage[0]->name);

  /* Do history  */
  inHistory  = newObitDataHistory ((ObitData*)inImage[0], OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory ((ObitData*)outImage, OBIT_IO_WriteOnly, err);
  ObitHistoryCopyHeader (inHistory, outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage[0]->name);
  
  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  for (i=0; i<numImage; i++) {
    g_snprintf (hicard, 80, "%s / input %d = %s",pgmName, i+1, inImage[i]->name);
    ObitHistoryWriteRec (outHistory, -1, hicard, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  } /* end loop adding input file names */

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Close outImage to update header */
  outImage->myStatus = OBIT_Modified;
  ObitImageClose (outImage, err);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end doHistory */
