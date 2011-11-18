/* $Id$  */
/* Combine beam images           */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010,2011                                          */
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
#include "ObitImageInterp.h"
#include "ObitImageUtil.h"
#include "ObitHistory.h"
#include "ObitTableFQUtil.h"

#define MAXINPUT 10  /* Maximum number of input images */
/* internal prototypes */
/* Get inputs */
ObitInfoList* ComBeamIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void ComBeamOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Get images from inputs */
void ComBeamGetImage(ObitInfoList *myInput, ObitErr *err);
/* ComBeam images together to outImage */
void doComBeam (ObitInfoList *myInput, ObitErr *err);
/* Write History */
void doHistory (ObitErr *err);

/* Program globals */
gchar *pgmName = "ComBeam";         /* Program name */
gchar *infile  = "ComBeam.in" ;     /* File with program inputs */
gchar *outfile = "ComBeam.out";     /* File to contain program outputs */
olong  pgmNumber;                   /* Program number (like POPS no.) */
olong  AIPSuser;                    /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;                     /* Number of AIPS directories */
gchar **AIPSdirs=NULL;              /* List of AIPS data directories */
olong  nFITS=0;                     /* Number of FITS directories */
gchar **FITSdirs=NULL;              /* List of FITS data directories */
ObitImage *outImage;                /* output image */
ObitImage *inImage[MAXINPUT];       /* Input images */
ObitImageInterp *inInterp[MAXINPUT]; /* Input interpolators */
olong nImage;                       /* Number of input images */
olong nPlane[MAXINPUT];             /* Number of planes per image */
odouble *freq[MAXINPUT];            /* Freq. (Hz) per plane per image */
ObitInfoList *myInput  = NULL;      /* Input parameter list */
ObitInfoList *myOutput = NULL;      /* Output parameter list */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/* Combine beam images                                                    */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitErr      *err= NULL;
  olong         i;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = ComBeamIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Get list of input images and output image */
  ComBeamGetImage(myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Combine them together */
  doComBeam (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* History */
  doHistory (err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* cleanup */
  for (i=0; i<nImage; i++) {
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

ObitInfoList* ComBeamIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "ComBeamIn";

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

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(list, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  return list;
} /* end ComBeamIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: ComBeam -input file -output ofile [args]\n");
    fprintf(stderr, "Combines beam images\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def ComBeam.in\n");
    fprintf(stderr, "  -output output result file, def ComBeam.out\n");
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
/*  Get images from myInput, convert to ImageInterp                       */
/*  Values:                                                               */
/*      myImage   InfoList with inputs                                    */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void ComBeamGetImage(ObitInfoList *myInput, ObitErr *err)
{
  ObitInfoType type;
  ObitImage    *tImage=NULL;
  ObitTableFQ  *FQTab=NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong        i, j, k, Aseq, disk, cno;
  ofloat       ftemp, deltFreq, xCells, yCells, *chBandw=NULL;
  odouble      refFreq, *freqOff=NULL;
  olong        number, nfreq, nx, ny, FQver;
  oint         fqid, nif, *sideBand=NULL;
  gboolean     exist;
  gchar        *strTemp=NULL, inFile[128];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gchar        tname[101], *Type, *FITS="FITS";
  gchar *routine = "ComBeamGetImage";

  if (err->error) return;  /* existing error? */

  /* Get number of input images */
  ObitInfoListGet(myInput, "numImage", &type, dim, &number, err);
  if (err->error) Obit_traceback_msg (err, routine, routine);
  nImage = number;

   /* In range? */
  Obit_return_if_fail(((number>=1)&&(number<=MAXINPUT)), err,
		      "%s: Number of input images,  %d, not in range [1,%d]", 
		      routine, number, MAXINPUT);
  
  /* Output info */
  ObitInfoListGet(myInput, "nfreq",    &type, dim, &nfreq, err);
  ObitInfoListGet(myInput, "deltFreq", &type, dim, &deltFreq, err);
  ObitInfoListGet(myInput, "refFreq",  &type, dim, &refFreq, err);
  ObitInfoListGet(myInput, "nx",       &type, dim, &nx, err);
  ObitInfoListGet(myInput, "ny",       &type, dim, &ny, err);
  ObitInfoListGet(myInput, "xCells",   &type, dim, &xCells, err);
  ObitInfoListGet(myInput, "yCells",   &type, dim, &yCells, err);
  if (err->error) Obit_traceback_msg (err, routine, routine);
  /* to Deg */
  xCells /= 3600.0;
  yCells /= 3600.0;

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

      /* Make interpolator */
      inInterp[i] = ObitImageInterpCreate(inImage[i]->name, inImage[i], 2, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
      /* Interpolator info */
      nPlane[i] = inInterp[i]->nplanes;
      freq[i]   = inInterp[i]->freqs;  /* just pointer*/
      
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

      /* Make interpolator */
      inInterp[i] = ObitImageInterpCreate(inImage[i]->name, inImage[i], 2, err);
      if (err->error) Obit_traceback_msg (err, routine, routine);
      /* Interpolator info */
      nPlane[i] = inInterp[i]->nplanes;
      freq[i]   = inInterp[i]->freqs;  /* just pointer*/

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
    for (k=0; k<6; k++) Aclass[k] = ' '; Aclass[k] = 0;
    ObitInfoListGet(myInput, "outClass", &type, dim, Aclass, err);
    Aclass[dim[0]] = 0;
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

  /* Axes for output image */
  outImage->myDesc = ObitImageDescCopy(inInterp[0]->ImgDesc, outImage->myDesc, 
				       err);
  if (err->error) Obit_traceback_msg (err, routine, inInterp[0]->name);
  outImage->myDesc->inaxes[0] = nx;
  outImage->myDesc->crpix[0]  = nx/2 + 1;
  outImage->myDesc->cdelt[0]  = xCells;
  outImage->myDesc->inaxes[1] = ny;
  outImage->myDesc->crpix[1]  = ny/2 + 1;
  outImage->myDesc->cdelt[1]  = yCells;
  outImage->myDesc->crval[outImage->myDesc->jlocf]   = refFreq;
  outImage->myDesc->cdelt[outImage->myDesc->jlocf]   = deltFreq;
  outImage->myDesc->inaxes[outImage->myDesc->jlocif] = nfreq;

  /* Make sure it's OK */
  ObitImageFullInstantiate (outImage, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Frequency to FQ table */
  nif   = nfreq;
  FQver = 1;
  fqid  = 1;
  chBandw  = g_malloc0(nfreq*sizeof(ofloat));
  sideBand = g_malloc0(nfreq*sizeof(oint));
  freqOff  = g_malloc0(nfreq*sizeof(odouble));
  for (i=0; i<nif; i++) {
    sideBand[i] = 1;
    chBandw[i]  = deltFreq;
    freqOff[i]  = i * deltFreq;
  }

  /* FQ Table */
  FQTab = 
    newObitTableFQValue (outImage->name, (ObitData*)outImage, &FQver, 
			 OBIT_IO_WriteOnly, nif, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);


  ObitTableFQPutInfo (FQTab, fqid, nif, freqOff, sideBand, chBandw, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Cleanup */
  if (freqOff)  g_free(freqOff);
  if (sideBand) g_free(sideBand);
  if (chBandw)  g_free(chBandw);
  FQTab = ObitTableFQUnref(FQTab);
} /* end ComBeamGetImage */

/*----------------------------------------------------------------------- */
/*   ComBeam merges images together to outImage                           */
/*  Values:                                                               */
/*      myImage   InfoList with inputs                                    */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void doComBeam (ObitInfoList *myInput, ObitErr *err)
{
  ObitInfoType type;
  gint32     dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitFArray *outData=NULL;
  olong      i, ix, iy, iz, ip,  totPlanes, adim[2], indx;
  olong      nfreq, nx, ny, iplane, jplane, plane[5] = {1,1,1,1,1};;
  ofloat     deltFreq, xCells, yCells, value, *wt=NULL;
  ofloat     fblank = ObitMagicF();
  odouble    refFreq, freq, x, y, xs, ys, dif, sum, sumwt;
  gchar *routine = "doComBeam";

  if (err->error) return;  /* existing error? */

   /* Output info */
  ObitInfoListGet(myInput, "nfreq",    &type, dim, &nfreq, err);
  ObitInfoListGet(myInput, "deltFreq", &type, dim, &deltFreq, err);
  ObitInfoListGet(myInput, "refFreq",  &type, dim, &refFreq, err);
  ObitInfoListGet(myInput, "nx",       &type, dim, &nx, err);
  ObitInfoListGet(myInput, "ny",       &type, dim, &ny, err);
  ObitInfoListGet(myInput, "xCells",   &type, dim, &xCells, err);
  ObitInfoListGet(myInput, "yCells",   &type, dim, &yCells, err);
  if (err->error) Obit_traceback_msg (err, routine, routine);
  /* to Deg */
  xCells /= 3600.0;
  yCells /= 3600.0;

  /* Total number of input planes */
  totPlanes = 0;
  for (i=0; i<nImage; i++) totPlanes += nPlane[i];

  /* Create work FArrays */
  adim[0] = nx; adim[1] = ny;
  outData = ObitFArrayCreate ("out", 2, adim);
  wt =  g_malloc0(totPlanes*sizeof(ofloat));

  /* Loop over output planes */
  for (iplane = 0; iplane<nfreq; iplane++) {
    freq = refFreq + iplane*deltFreq;  /* Frequency for plane */

    /* Determine weights */
    jplane = 0;
    for (iz=0; iz<nImage; iz++) {
      /* Loop over input planes */
      for (ip=0; ip<nPlane[iz]; ip++) {
	dif = 1.0e-9 * fabs (freq - inInterp[iz]->freqs[ip]);
	dif = MAX (dif, 0.001);
	wt[jplane++] = 1.0 / dif;
      }
    }

    /* Loop over image */
    for (iy=0; iy<ny; iy++) {
      y = (iy-ny/2) * yCells;
      for (ix=0; ix<nx; ix++) {
	x = (ix-nx/2) * xCells;

	/* Loop over input images */
	sum = sumwt = 0.0;
	jplane = 0;
	for (iz=0; iz<nImage; iz++) {
	  /* Loop over input planes */
	  for (ip=0; ip<nPlane[iz]; ip++) {
	    /* Scale location for difference in frequency */
	    xs = x * freq / inInterp[iz]->freqs[ip];
	    ys = y * freq / inInterp[iz]->freqs[ip];
	    /* Interpolate */
	    value = ObitImageInterpValue (inInterp[iz], xs, ys, 0.0, ip, err);
	    if (value!=fblank) {
	      /*DEBUG if ((ix==26)&&(iy==26))
		fprintf (stderr, "%d %d %f %f %f %f\n", iz, ip, xs, ys, value, wt[jplane]);*/
	      sum   += value * wt[jplane];
	      sumwt += wt[jplane++];
	    }
	  } /* end plane loop */
	} /* end image loop */

	/* save value */
	indx = ix + iy*nx;
	if (sumwt>0.0) {
	  outData->array[indx] = sum / sumwt;
	  /*DEBUG if ((ix==26)&&(iy==26))
	    fprintf (stderr, "%d %d %d %f %f %f %lf %lf\n", 
	    iplane, ix, iy, x, y, outData->array[indx], sum, sumwt);*/
	} else {  /* No value */
	  outData->array[indx] = fblank;
	}
	
      } /* end loop over y */
    } /* end loop over y */

    /* Write plane */
    plane[0] = iplane + 1;
    ObitImagePutPlane (outImage, outData->array, plane, err);

  } /* end loop over planes */

  /* Cleanup */
  if (wt) g_free(wt);
  if (outData) ObitFArrayUnref(outData);

 } /* end doComBeam */

/*----------------------------------------------------------------------- */
/*  Write history                                                         */
/*  Values:                                                               */
/*      nImage    Number of images in inImage                             */
/*      inImage   Array of ObitImage pointers                             */
/*      outImage  Output ObitImage pointer                                */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void doHistory (ObitErr *err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  olong         i;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "nx", "ny", "nfreq", "xCells", "yCells", "deltFreq", "refFreq", 
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
  for (i=0; i<nImage; i++) {
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
