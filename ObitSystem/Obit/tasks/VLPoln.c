/* $Id:  */
/* Extract polarization info for VL table                             */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2017,2018                                          */
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
#include "ObitHistory.h"
#include "ObitImage.h"
#include "ObitData.h"
#include "ObitTableVL.h"
#include "ObitFInterpolate.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* VLPolnIn (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Get input data */
ObitImage* getInputData (ObitInfoList *myInput, ObitErr *err);
/* Get input Q data */
ObitImage* getInputQData (ObitInfoList *myInput, ObitErr *err);
/* Get input U data */
ObitImage* getInputUData (ObitInfoList *myInput, ObitErr *err);
/* Get input RM data */
ObitImage* getInputRMData (ObitInfoList *myInput, ObitErr *err);
/* Extract Poln data */
void ExtractPoln (ObitImage* inData, ObitImage* QData, ObitImage* UData, ObitErr *err);
/* Extract RM data */
void ExtractRM (ObitImage* inData, ObitImage* RMData, ObitErr *err);
/* Write history */
void VLPolnHistory (ObitInfoList* myInput, ObitImage* inData, ObitErr* err);
/* Debias poln amp */
void PolnDeBias (ofloat *p, ofloat rms);

/* Program globals */
gchar *pgmName = "VLPoln";       /* Program name */
gchar *infile  = "VLPoln.in" ;   /* File with program inputs */
gchar *outfile = "VLPoln.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
olong prtLv=0;          /* Diagnostic print level */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
gboolean haveRM      = FALSE;  /* Have RM cube? */

/*----------------- Macroes ---------------------------*/

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/* Obit task which Convert VL table to FS table, filters                  */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitImage    *inData = NULL, *QData = NULL, *UData = NULL, *RMData = NULL;
  ObitErr      *err= NULL;
  gchar        *dataParms[] = {  /* Control parameter */
    "RMSsize",  "inVL", NULL};

  /* Startup - parse command line */
  err = newObitErr();
  myInput = VLPolnIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get input images */
  inData = getInputData (myInput, err);
  QData  = getInputQData (myInput, err);
  UData  = getInputUData (myInput, err);
  RMData = getInputRMData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Copy control parameters to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);

  /* Extract Poln info to VL Table */
  ExtractPoln (inData, QData, UData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Extract RM if Given */
  if (haveRM) {
    ExtractRM (inData, RMData, err);
    if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  }

  /* Write history */
  VLPolnHistory (myInput, inData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  
  /* Shutdown Obit */
 exit:
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput  = ObitInfoListUnref(myOutput);
  inData    = ObitImageUnref(inData);
  QData     = ObitImageUnref(QData);
  UData     = ObitImageUnref(UData);
  RMData    = ObitImageUnref(RMData);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* VLPolnIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "VLPolnIn";

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
} /* end VLPolnIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: VLPoln -input file -output ofile [args]\n");
    fprintf(stderr, "Convert VL table to FS \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def VLPoln.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def VLPoln.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input image\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2\n");
    fprintf(stderr, "  -inFile input FITS uvdata file\n");
    fprintf(stderr, "  -inName input AIPS uvdata file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input (AIPS or FITS) disk number (1-rel) \n");
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
/*     inFile    Str [?]    input FITS image file name [def "Image.fits"] */
/*     inName    Str [12]   input AIPS image name  [no def]               */
/*     inClass   Str [6]    input AIPS image class  [no def]              */
/*     inSeq     Int        input AIPS image sequence no  [no def]        */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
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
  strTemp = "VLPoln.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "VLPolnName";
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
/*  Get input data,                                                       */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitImage with input data                                        */
/*----------------------------------------------------------------------- */
ObitImage* getInputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitImage    *inData = NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong         Aseq, disk, cno;
  gchar        *Type, *strTemp, inFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gchar *routine = "getInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input Image data Object */
  inData = newObitImage("input Image data");
  
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
      if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "inSeq", OBIT_oint, dim, &Aseq);
    }

    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
    /* define object */
    ObitImageSetAIPS (inData, OBIT_IO_byPlane, disk, cno, AIPSuser, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    if (ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp)) {
      ObitTrimTrail(strTemp);
      strncpy (inFile, strTemp, 127);
    } else { 
      strncpy (inFile, "No_Filename_Given", 127);
    }
    
    /* input FITS disk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);

    /* define object */
    ObitImageSetFITS (inData, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return inData;
  }

  /* Full instantiation - needed to initialize tableList */
  ObitImageFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Get input Q Pol data,                                                 */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitImage with input data                                        */
/*----------------------------------------------------------------------- */
ObitImage* getInputQData (ObitInfoList *myInput, ObitErr *err)
{
  ObitImage    *QData = NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong        Aseq, disk, cno;
  gchar        *Type, *strTemp, QFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gchar *routine = "getInputQData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return QData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input Image data Object */
  QData = newObitImage("input Q Image");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
    /* input AIPS disk */
    ObitInfoListGet(myInput, "QDisk", &type, dim, &disk, err);
    /* input AIPS name */
    if (ObitInfoListGetP(myInput, "QName", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aname, strTemp, 13);
    } else { /* Didn't find */
      strncpy (Aname, "No Name ", 13);
    } 
    Aname[12] = 0;
    /* input AIPS class */
    if  (ObitInfoListGetP(myInput, "QClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    Aclass[6] = 0;
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "QSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 want highest existing sequence */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", QData);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "QSeq", OBIT_oint, dim, &Aseq);
    }

    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", QData);
    
    /* define object */
    ObitImageSetAIPS (QData, OBIT_IO_byPlane, disk, cno, AIPSuser, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", QData);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    if (ObitInfoListGetP(myInput, "QFile", &type, dim, (gpointer)&strTemp)) {
      ObitTrimTrail(strTemp);
      strncpy (QFile, strTemp, 127);
    } else { 
      strncpy (QFile, "No_Filename_Given", 127);
    }
    
    /* input FITS disk */
    ObitInfoListGet(myInput, "QDisk", &type, dim, &disk, err);

    /* define object */
    ObitImageSetFITS (QData, OBIT_IO_byPlane, disk, QFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", QData);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return QData;
  }

  /* Full instantiation */
  ObitImageFullInstantiate (QData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", QData);

  return QData;
} /* end getInputQData */

/*----------------------------------------------------------------------- */
/*  Get input U Pol data,                                                 */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitImage with input data                                        */
/*----------------------------------------------------------------------- */
ObitImage* getInputUData (ObitInfoList *myInput, ObitErr *err)
{
  ObitImage    *UData = NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong        Aseq, disk, cno;
  gchar        *Type, *strTemp, UFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gchar *routine = "getInputUData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return UData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input Image data Object */
  UData = newObitImage("input U Image");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
    /* input AIPS disk */
    ObitInfoListGet(myInput, "UDisk", &type, dim, &disk, err);
    /* input AIPS name */
    if (ObitInfoListGetP(myInput, "UName", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aname, strTemp, 13);
    } else { /* Didn't find */
      strncpy (Aname, "No Name ", 13);
    } 
    Aname[12] = 0;
    /* input AIPS class */
    if  (ObitInfoListGetP(myInput, "UClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    Aclass[6] = 0;
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "USeq", &type, dim, &Aseq, err);

    /* if ASeq==0 want highest existing sequence */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", UData);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "USeq", OBIT_oint, dim, &Aseq);
    }

    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", UData);
    
    /* define object */
    ObitImageSetAIPS (UData, OBIT_IO_byPlane, disk, cno, AIPSuser, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", UData);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    if (ObitInfoListGetP(myInput, "UFile", &type, dim, (gpointer)&strTemp)) {
      ObitTrimTrail(strTemp);
      strncpy (UFile, strTemp, 127);
    } else { 
      strncpy (UFile, "No_Filename_Given", 127);
    }
    
    /* input FITS disk */
    ObitInfoListGet(myInput, "UDisk", &type, dim, &disk, err);

    /* define object */
    ObitImageSetFITS (UData, OBIT_IO_byPlane, disk, UFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", UData);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return UData;
  }

  /* Full instantiation */
  ObitImageFullInstantiate (UData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", UData);

  return UData;
} /* end getInputUData */

/*----------------------------------------------------------------------- */
/*  Get input RM Pol data, returns NULL of not given and sets haveRM=F    */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitImage with input data                                        */
/*----------------------------------------------------------------------- */
ObitImage* getInputRMData (ObitInfoList *myInput, ObitErr *err)
{
  ObitImage    *RMData = NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong        Aseq, disk, cno;
  gchar        *Type, *strTemp, RMFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gchar *routine = "getInputRMData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return RMData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input Image data Object */
  RMData = newObitImage("input RM Image");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
   /* input AIPS disk */
    ObitInfoListGet(myInput, "RMDisk", &type, dim, &disk, err);
    /* input AIPS name */
    if (ObitInfoListGetP(myInput, "RMName", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aname, strTemp, 13);
      /* Check if given */
      if (!strncmp (Aname, "            ", 12)) {
	haveRM = FALSE;
	return NULL;
      }
    } else { /* Didn't find */
      haveRM = FALSE;
      return NULL;
    } 
    Aname[12] = 0;
    /* input AIPS class */
    if  (ObitInfoListGetP(myInput, "RMClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    Aclass[6] = 0;
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "RMSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 want highest existing sequence */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", RMData);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "RMSeq", OBIT_oint, dim, &Aseq);
    }

    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", RMData);
    
    /* define object */
    ObitImageSetAIPS (RMData, OBIT_IO_byPlane, disk, cno, AIPSuser, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", RMData);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    if (ObitInfoListGetP(myInput, "RMFile", &type, dim, (gpointer)&strTemp)) {
      /* Check if given */
      if (!strncmp (strTemp, "            ", 12)) {
	haveRM = FALSE;
	return NULL;
      }
      ObitTrimTrail(strTemp);
      strncpy (RMFile, strTemp, 127);
    } else { 
      haveRM = FALSE;
      return NULL;
    }
    
    /* input FITS disk */
    ObitInfoListGet(myInput, "RMDisk", &type, dim, &disk, err);

    /* define object */
    ObitImageSetFITS (RMData, OBIT_IO_byPlane, disk, RMFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", RMData);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return RMData;
  }

  /* Full instantiation  */
  ObitImageFullInstantiate (RMData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", RMData);

  haveRM = TRUE;   /* Must be OK if it gets here */
  return RMData;
} /* end getInputRMData */

/*----------------------------------------------------------------------- */
/*  Write History for VLPoln                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitImage to write history to                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void VLPolnHistory (ObitInfoList* myInput, ObitImage* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", "inVL", 
    "QFile",  "QDisk", "QName", "QClass", "QSeq", 
    "UFile",  "UDisk", "UName", "UClass", "USeq", 
    "RMFile",  "RMDisk", "RMName", "RMClass", "RMSeq",
    "RMSsize",
    NULL};
  gchar *routine = "VLPolnHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitImageIsA(inData));

   /* Do history */
  outHistory = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadWrite, err);

  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  outHistory = ObitHistoryUnref(outHistory); /* cleanup */
 
} /* end VLPolnHistory  */

void ExtractPoln (ObitImage* inData, ObitImage* QData, ObitImage* UData, ObitErr *err);
/**
 * Extract polarization  info from images, update VL table
 * \param inData    Input file with VL Table FS table
 *                  Control parameter on info
 * \li "RMSsize"    OBIT_float (1,1,1) halfwidth of region to determine RMS [def 50]
 *                  <=0 ->use full image
 * \param QData     Q Poln image
 * \param UData     U Poln image
 * \param err       ObitErr error stack.
 */
void ExtractPoln (ObitImage* inData, ObitImage* QData, ObitImage* UData, ObitErr *err)
{
  ObitTableVL *inVL=NULL;
  ObitTableVLRow *VLrow=NULL;
  ObitFArray *local=NULL;
  ObitIOSize IOBy;
  ObitFInterpolate *Qinterp=NULL, *Uinterp=NULL;
  ObitInfoType type;
  ofloat pflux, RMS, QRMS, URMS, lQRMS, lURMS, Qval, Uval, Qpixel[2], Upixel[2];
  ofloat fblank = ObitMagicF();
  odouble coord[2];
  olong irow, RMSsize, iPlane, blc[2], trc[2], Plane[5]={1,1,1,1,1};
  olong  nx, ny, hwidth, iVer, count = 0;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean toGal;
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitTableFSGetSpectrum";

  /* error checks */
  if (err->error) return;

/* Get parameters */
  RMSsize = 50;
  InfoReal.otg = RMSsize; type = OBIT_long;
  ObitInfoListGetTest(inData->info, "RMSsize",  &type, dim,  &InfoReal);
  if (type==OBIT_float) RMSsize = InfoReal.flt + 0.5;
  else if (type==OBIT_long)  RMSsize = InfoReal.otg;
  else if (type==OBIT_oint)  RMSsize = InfoReal.itg;
  /* VL table */
  iVer = 0;
  ObitInfoListGetTest(inData->info, "inVL",  &type, dim,  &iVer);
		    
   /* Get VL Table */
  inVL = newObitTableVLValue (inData->name, (ObitData*)inData, &iVer, 
				 OBIT_IO_ReadWrite, err);
  /* Open table */
  ObitTableVLOpen (inVL, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  VLrow  = newObitTableVLRow (inVL);
  ObitTableSetRow ((ObitTable*)inVL, (ObitTableRow*)VLrow, err);
  if (err->error) goto cleanup;

  /* Anything to work on? */
  if (inVL->myDesc->nrow<=0) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: NO entries in VL table", 
		   routine);
    goto cleanup;
  }

  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (QData->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  ObitInfoListPut (UData->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);

  /* Open images */
  if ((ObitImageOpen (QData, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, QData->name);
    goto cleanup;
  }
  if ((ObitImageOpen (UData, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, UData->name);
    goto cleanup;
  }

  /* Numbers of things */
  nx  = QData->myDesc->inaxes[QData->myDesc->jlocr];
  ny  = QData->myDesc->inaxes[QData->myDesc->jlocd];

  /* Need to convert position to Galactic? */
  toGal = (QData->myDesc->ctype[0][0]=='G' && QData->myDesc->ctype[0][1]=='L') &&
    !(inData->myDesc->ctype[0][0]=='G' && inData->myDesc->ctype[0][1]=='L');

  /* Make interpolators */
  hwidth = 2;
  Qinterp = newObitFInterpolateCreate ("Interpolator", QData->image, 
				      QData->myDesc, hwidth);
  Uinterp = newObitFInterpolateCreate ("Interpolator", UData->image, 
				      UData->myDesc, hwidth);
  /* Only interested in one plane */
  iPlane = 1;
  /* Read input plane */
  Plane[0] = iPlane;
  if ((ObitImageGetPlane (QData, NULL, Plane, err)
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		   routine, QData->name);
    goto cleanup;
  }
  if ((ObitImageGetPlane (UData, NULL, Plane, err)
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		   routine, UData->name);
    goto cleanup;
  }

  /* Plane statistics */
  QRMS    = ObitFArrayRMS (QData->image);
  URMS    = ObitFArrayRMS (UData->image);
  
  /* Set interpolators */
  ObitFInterpolateReplace (Qinterp, QData->image);
  ObitFInterpolateReplace (Uinterp, UData->image);
  
  /* Loop over VL Table */
  for (irow=1; irow<=inVL->myDesc->nrow; irow++) {
    ObitTableVLReadRow (inVL, irow, VLrow, err);
    if (err->error) goto cleanup;
      
    /* Want this one? */
    if (VLrow->status<0) continue;  /* Skip deselected record */
      
    /* Is this in QData? */
    coord[0] = VLrow->Ra2000;
    coord[1] = VLrow->Dec2000;

    /* Convert to Galactic? */
    if (toGal) {
      ObitSkyGeomJtoB  (&coord[0], &coord[1]);
      ObitSkyGeomEq2Gal(&coord[0], &coord[1]);
    }

    ObitImageDescGetPixel(QData->myDesc, coord, Qpixel, err);
    ObitImageDescGetPixel(UData->myDesc, coord, Upixel, err);
    if (err->error) {
      /* Can't deal with this one */
      ObitErrClear(err);
      continue;
    }
    if ((Qpixel[0]<0) || (Qpixel[1]<0) || (Qpixel[0]>=nx) || (Qpixel[1]>=ny)) continue;
    if ((Upixel[0]<0) || (Upixel[1]<0) || (Upixel[0]>=nx) || (Upixel[1]>=ny)) continue;
      
    Qval = ObitFInterpolatePosition(Qinterp, coord, err);
    Uval = ObitFInterpolatePosition(Uinterp, coord, err);
    if (err->error) goto cleanup;

    /* Local RMS? */
    if (RMSsize>0) {
      /* Q */
      blc[0] = MAX (1, (olong)(Qpixel[0] - RMSsize + 0.5));
      blc[1] = MAX (1, (olong)(Qpixel[1] - RMSsize + 0.5));
      trc[0] = MIN (nx-2, (olong)(Qpixel[0] + RMSsize + 0.5));
      trc[1] = MIN (ny-2, (olong)(Qpixel[1] + RMSsize + 0.5));
      local = ObitFArraySubArr (QData->image, blc, trc, err);
      if (err->error) goto cleanup;
      lQRMS   = ObitFArrayRMS (local);
      local   = ObitFArrayUnref(local);
      /* u */
      blc[0] = MAX (1, (olong)(Upixel[0] - RMSsize + 0.5));
      blc[1] = MAX (1, (olong)(Upixel[1] - RMSsize + 0.5));
      trc[0] = MIN (nx-2, (olong)(Upixel[0] + RMSsize + 0.5));
      trc[1] = MIN (ny-2, (olong)(Upixel[1] + RMSsize + 0.5));
      local = ObitFArraySubArr (UData->image, blc, trc, err);
      if (err->error) goto cleanup;
      lURMS   = ObitFArrayRMS (local);
      local   = ObitFArrayUnref(local);
    } else {  /* Use plane */
      lQRMS    = QRMS;
      lURMS    = URMS;
    }

    /* Update */
    RMS = (lQRMS + lURMS)*0.5;  /* Average RMS */
    VLrow->PolRMS  = RMS;
    VLrow->QCenter = Qval;
    VLrow->UCenter = Uval;
    if ((Qval!=fblank) && (Qval!=fblank)) {
      pflux = sqrt(Qval*Qval+Uval*Uval);
      /* Debias */
      PolnDeBias(&pflux, RMS);
    } else pflux = fblank;
    VLrow->PFlux   = pflux;
    
    /* reWrite row */
    count++;
    ObitTableVLWriteRow (inVL, irow, VLrow, err);
    if (err->error) goto cleanup;
  } /* end loop over  VL table */

cleanup:
  Obit_log_error(err, OBIT_InfoErr, "Updated %d polarization entries", count);
  /* Close */
  ObitImageClose (QData, err);
  ObitImageClose (UData, err);
  ObitTableVLClose (inVL,  err);
  /* Free image buffers if not memory resident */
  if (QData->mySel->FileType!=OBIT_IO_MEM) 
    QData->image = ObitFArrayUnref(QData->image);
  if (UData->mySel->FileType!=OBIT_IO_MEM) 
    UData->image = ObitFArrayUnref(UData->image);
  /* release objects */
  Qinterp = ObitFInterpolateUnref(Qinterp);
  Uinterp = ObitFInterpolateUnref(Uinterp);
  VLrow = ObitTableVLRowUnref(VLrow);
  inVL = ObitTableVLUnref(inVL);
  if (err->error) Obit_traceback_msg (err, routine, QData->name);
} /* end ExtractPoln */

/**
 * Extract Rotation measure  info from cube, update VL table
 * Also copy error if in cube
 * \param inData    Input file with VL Table FS table
 * \param QData     Q Poln image
 * \param UData     U Poln image
 * \param err       ObitErr error stack.
 */
void ExtractRM (ObitImage* inData, ObitImage* RMData, ObitErr *err)
{
  ObitTableVL *inVL=NULL;
  ObitTableVLRow *VLrow=NULL;
  ObitIOSize IOBy;
  ObitFInterpolate *Interp=NULL;
  ObitInfoType type;
  ofloat val, pixel[2];
  ofloat fblank = ObitMagicF();
  odouble coord[2];
  olong irow, iPlane, Plane[5]={1,1,1,1,1};
  olong nx, ny, hwidth, iVer, nPlane, count=0;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean toGal, doRMSyn;
  gchar *routine = "ObitTableFSGetSpectrum";

  /* error checks */
  if (err->error) return;

  /* Use RM synthesis or least squares cubes? */
  if (!strncmp(RMData->myDesc->ctype[2],"MaxRMSyn",8))
    doRMSyn = TRUE;
  else doRMSyn = FALSE;
  if (doRMSyn) Obit_log_error(err, OBIT_InfoErr, "Using RM Synthesis cube.");
  ObitErrLog(err); /* show any error messages on err */

  /* Get parameters */
  /* VL table */
  iVer = 0;
  ObitInfoListGetTest(inData->info, "inVL",  &type, dim,  &iVer);
  
  /* Get VL Table */
  inVL = newObitTableVLValue (inData->name, (ObitData*)inData, &iVer, 
			      OBIT_IO_ReadWrite, err);
  /* Open table */
  ObitTableVLOpen (inVL, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  VLrow  = newObitTableVLRow (inVL);
  ObitTableSetRow ((ObitTable*)inVL, (ObitTableRow*)VLrow, err);
  if (err->error) goto cleanup;
  
  /* Anything to work on? */
  if (inVL->myDesc->nrow<=0) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: NO entries in VL table", 
		   routine);
    goto cleanup;
  }
  
  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (RMData->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);

  /* Open image */
  if ((ObitImageOpen (RMData, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, RMData->name);
    goto cleanup;
  }

  /* Numbers of things */
  nx  = RMData->myDesc->inaxes[RMData->myDesc->jlocr];
  ny  = RMData->myDesc->inaxes[RMData->myDesc->jlocd];
  nPlane = MIN(2,RMData->myDesc->inaxes[2]) ; 
  if (RMData->myDesc->inaxes[2]>=4) nPlane = 4;

  /* Need to convert position to Galactic? */
  toGal = (RMData->myDesc->ctype[0][0]=='G' && RMData->myDesc->ctype[0][1]=='L') &&
    !(inData->myDesc->ctype[0][0]=='G' && inData->myDesc->ctype[0][1]=='L');

  /* Make interpolator */
  hwidth = 3;
  Interp = newObitFInterpolateCreate ("Interpolator", RMData->image, 
				      RMData->myDesc, hwidth);
  /* Loop over planes */
  for (iPlane=1; iPlane<=nPlane; iPlane++) {
    /* Read input plane */
    Plane[0] = iPlane;
    if ((ObitImageGetPlane (RMData, NULL, Plane, err)
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		     routine, RMData->name);
      goto cleanup;
    }

  /* Set interpolator */
  ObitFInterpolateReplace (Interp, RMData->image);
  
  /* Loop over VL Table */
  for (irow=1; irow<=inVL->myDesc->nrow; irow++) {
    ObitTableVLReadRow (inVL, irow, VLrow, err);
    if (err->error) goto cleanup;
    
    /* Want this one? */
    if (VLrow->status<0) continue;  /* Skip deselected record */
    
    /* Is this in RMData? */
    coord[0] = VLrow->Ra2000;
    coord[1] = VLrow->Dec2000;
    
    /* Convert to Galactic? */
    if (toGal) {
      ObitSkyGeomJtoB  (&coord[0], &coord[1]);
      ObitSkyGeomEq2Gal(&coord[0], &coord[1]);
    }
    
    ObitImageDescGetPixel(RMData->myDesc, coord, pixel, err);
    if (err->error) {
      /* Can't deal with this one */
      ObitErrClear(err);
      continue;
    }
    if ((pixel[0]<0) || (pixel[1]<0) || (pixel[0]>=nx) || (pixel[1]>=ny)) continue;
    
    val = ObitFInterpolatePosition(Interp, coord, err);
    if (err->error) goto cleanup;
    
    /* Update */
    if (doRMSyn) {  /* RM Synthesis -> get pamp from plane 3, no errors */
      if (iPlane==1) VLrow->RM      = val; VLrow->RMerr   = -1.0;
      if (iPlane==3) {
	if (VLrow->RM==fblank) val = fblank;
	if (val!=fblank) PolnDeBias(&val, VLrow->PolRMS); 	
	VLrow->PFlux = val;}
      if (iPlane==2) VLrow->EVPA    = val*RAD2DG;   /* In deg */
      if (iPlane==4) VLrow->EVPAerr = -1.0;
    } else {       /* Least squares RM with errors */
      if (iPlane==1) VLrow->RM      = val;
      if (iPlane==3) VLrow->RMerr   = val;
      if (iPlane==2) VLrow->EVPA    = val*RAD2DG;   /* In deg */
      if (iPlane==4) VLrow->EVPAerr = val*RAD2DG;
    }
    
    /* reWrite row */
    count++;
    ObitTableVLWriteRow (inVL, irow, VLrow, err);
    if (err->error) goto cleanup;
    } /* end loop over  VL table */
  } /* end loop over planes */

cleanup:
  Obit_log_error(err, OBIT_InfoErr, "Updated %d RM entries", count/nPlane);
  /* Close */
  ObitImageClose (RMData, err);
  ObitTableVLClose (inVL,  err);
  /* Free image buffers if not memory resident */
  if (RMData->mySel->FileType!=OBIT_IO_MEM) 
    RMData->image = ObitFArrayUnref(RMData->image);
  /* release objects */
  Interp = ObitFInterpolateUnref(Interp);
  VLrow = ObitTableVLRowUnref(VLrow);
  inVL = ObitTableVLUnref(inVL);
  if (err->error) Obit_traceback_msg (err, routine, RMData->name);
} /* end ExtractRM */

/**
 * Estimates the polarization bias in a polarization amplitude, p,  
 * measured in the presence of Q and U RMS noise, rms.  Returns the  
 * corrected value.  
 * The bias correction is such that the average bias is removed;  
 * thus the average in the absence of a signal is zero.  Does table  
 * lookup of values calculated by J. Condon in the range of p/rms of  
 * 1.253 (the no signal limit) and 4 (the high SNR regime).  Does  
 * second order Lagrange interpolation.  At lower values of P/RMS the  
 * bias is a constant 1.253*RMS. Above a signal-to-noise ratio of 4,  
 * use the formula:  
 * normalized bias = 1 / (2 * s) + 1 / (8 * s**3),  
 * where s is the true normalized flux density, iterating once to  
 * estimate s from the normalized map flux density.  "Normalized" means  
 * divided by the rms noise in the q and u maps.  
 * Routine translated from the AIPSish corerr.FOR/PDBIAS  
 *
 * \param p   On output, p is the estimated intrinsic total 
 *             polarized intensity. 
 * \param rms  The standard deviation of the (assumed equal) 
 *             Gaussian distributions of the Stokes Q or U maps. 
 */
void  PolnDeBias (ofloat *p, ofloat rms) 
{
  olong   i, index, i1, i2, i3;
  ofloat  pnorm, bias, d1, d2, d3, wt1, wt2, wt3, sum, sumwt;
  /* (map_flux,map_bias) pairs split into table1 and table2 */
  ofloat table1[] = {
      1.253, 1.256, 1.266, 1.281, 1.303, 1.330, 1.364, 1.402, 1.446, 
      1.495, 1.549, 1.606, 1.668, 1.734, 1.803, 1.875, 1.950, 2.027, 
      2.107, 2.189, 2.272, 2.358, 2.444, 2.532, 2.621, 2.711, 2.802, 
      2.894, 2.986, 3.079, 3.173, 3.267, 3.361, 3.456, 3.551, 3.646, 
      3.742, 3.838, 3.934, 4.031};
    ofloat table2[] = {
      1.253,  1.156,  1.066,  0.9814, 0.9030, 0.8304, 0.7636, 0.7023, 
      0.6462, 0.5951, 0.5486, 0.5064, 0.4683, 0.4339, 0.4028, 0.3749, 
      0.3498, 0.3273, 0.3070, 0.2888, 0.2724, 0.2576, 0.2442, 0.2321, 
      0.2212, 0.2112, 0.2021, 0.1938, 0.1861, 0.1791, 0.1726, 0.1666, 
      0.1610, 0.1557, 0.1509, 0.1463, 0.1420, 0.1380, 0.1342, 0.1306};

    /* Check RMS */
    if (rms <= 0.0) return;
    pnorm = (*p) / rms;
    
    /* Which regime? */
    if (pnorm <= table1[0]) {
      
      /* Low (no) SNR case */
      (*p) -= table2[0] * rms;
    } else if (pnorm >= table1[39]) {
      /* High SNR */
      bias = 1.0 / (2.0 * pnorm) + 1.0 / (8.0 * pnorm* pnorm* pnorm);
      pnorm = pnorm - bias;
      bias = 1.0 / (2.0 * pnorm) + 1.0 / (8.0 * pnorm*pnorm*pnorm);
      
      /* Correct for bias */
      *p -= bias * rms;
    } else {
      /* Middle, interpolate in table */
      index = 2;
      for (i= 3; i<=39; i++) {
	if (pnorm < table1[i-1]) break;
	index = i;
      } 
      /* Lagrange interpolation */
      i1 = index - 1;
      i2 = index;
      i3 = index + 1;
      d1 = (table1[i1-1] - table1[i2-1]) * (table1[i1-1] - table1[i3-1]);
      d2 = (table1[i2-1] - table1[i1-1]) * (table1[i2-1] - table1[i3-1]);
      d3 = (table1[i3-1] - table1[i1-1]) * (table1[i3-1] - table1[i2-1]);
      wt1 = (pnorm - table1[i2-1]) * (pnorm - table1[i3-1]) / d1;
      wt2 = (pnorm - table1[i1-1]) * (pnorm - table1[i3-1]) / d2;
      wt3 = (pnorm - table1[i1-1]) * (pnorm - table1[i2-1]) / d3;
      sum = table2[i1-1] * wt1 + table2[i2-1] * wt2 + table2[i3-1] * wt3;
      sumwt = wt1 + wt2 + wt3;
      if (sumwt > 0.0) {
	bias = sum / sumwt;
      } else {
	/* Shouldn't ever get here but do something reasonable. */
	bias = table2[i2-1];
      } 
      /* Correct for bias */
      *p -= bias * rms;
    } 
} /* end of routine PolnDeBias */ 

