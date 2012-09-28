/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2012                                          */
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
#include "ObitData.h"
#include "ObitTable.h"
#include "ObitTableList.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* TabCopyIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void TabCopyOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Digest inputs */
void digestInputs(ObitInfoList *myInput, ObitErr *err);
/* Get input data */
ObitData* getInputData (ObitInfoList *myInput, ObitErr *err);
/* Get output data */
ObitData* getOutputData (ObitInfoList *myInput, ObitErr *err);
/* Write history */
void TabCopyHistory (ObitInfoList* myInput, ObitData* inData, 
		     ObitData* outData, ObitErr* err);
/* Copy selected tables */
void TabCopyCopy (ObitInfoList* myInput, ObitData* inData,  
		  ObitData* outData, ObitErr* err);
/* Is table wanted? */
gboolean WantTable (ObitTable *inTab, gchar *KeyWord, ofloat KeyRange[2], 
		    gchar *KeyString, ObitErr* err);

/* Program globals */
gchar *pgmName = "TabCopy";       /* Program name */
gchar *infile  = "TabCopy.in" ;   /* File with program inputs */
gchar *outfile = "TabCopy.out";   /* File to contain program outputs */
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
/*   Obit task to copy tables arrached to ObitData objects                */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitData     *inData = NULL, *outData = NULL;
  ObitErr      *err= NULL;

  /* Startup - parse command line */
  err = newObitErr();
  myInput = TabCopyIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get input data */
  inData = getInputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* Get output data */
  outData = getOutputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Copy selected tables */
  TabCopyCopy (myInput, inData, outData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Write history */
  TabCopyHistory (myInput, inData, outData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  outData   = ObitUnref(outData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* TabCopyIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "TabCopyIn";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return list;

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
      
    } else if (strcmp(arg, "-inSeq") == 0) { /* AIPS sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "inSeq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-inDisk") == 0) { /* input disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "inDisk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-DataType") == 0) { /* File type AIPS or FITS */
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
  myOutput = defaultOutputs(err);
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  return list;
} /* end TabCopyIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: TabCopy -input file -output ofile [args]\n");
    fprintf(stderr, "Clip and smooth an SN table \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def TabCopy.in\n");
    fprintf(stderr, "  -output parameters (none) def TabCopy.out\n");
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
/*     inFile    Str [?]    input FITS image file name [def "TabCopy.fits"] */
/*     inName    Str [12]   input AIPS image name  [no def]               */
/*     inClass   Str [6]    input AIPS image class  [no def]              */
/*     inSeq     Int        input AIPS image sequence no  [no def]        */
/*     KeyWord   Str [?]    Keyword to test [blank = none'                */
/*     KeyRange  Flt (2)    Keyword value range, def=all                  */
/*     KeyString Str [?]    Keyword string to test, def=blank             */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat farray[3];
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
  strTemp = "TabCopy.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "TabCopyName";
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

  /* Keyword to test */
  strTemp = "        ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "KeyWord", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Key range  */
  dim[0] = 2;dim[1] = 1;
  farray[0] = -1.0e20; farray[1] = 1.0e20;
  ObitInfoListPut (out, "keyRange", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

   /* Key string to test */
  strTemp = "        ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "KeyString", OBIT_string, dim, strTemp, err);
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
  gboolean     doCalSelect;
  oint         doTabCopy;
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  /* Make sure doCalSelect set properly */
  doCalSelect = FALSE;
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doTabCopy = -1;
  ObitInfoListGetTest(myInput, "doTabCopy",  &type, dim, &doTabCopy);
  doCalSelect = doCalSelect || (doTabCopy>0);
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);

} /* end digestInputs */

/*----------------------------------------------------------------------- */
/*  Get input data, builds selector                                       */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitData with input data                                         */
/*----------------------------------------------------------------------- */
ObitData* getInputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitData       *inData = NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         Aseq, disk, cno;
  gchar        *Type, *strTemp, inFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "  ";
  gchar *routine = "getInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input data Object */
  inData = newObitData("input data");
  
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
    ObitDataSetAIPS (inData, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
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
    ObitDataSetFITS (inData, disk, inFile,  err); 
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return inData;
  }

  /* Full instantiation - needed to initialize tableList */
  ObitDataOpen (inData, OBIT_IO_ReadOnly, err);
  ObitDataClose (inData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Define output data object, output must exist                         */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitData for output image                                        */
/*----------------------------------------------------------------------- */
ObitData* getOutputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitData    *outData = NULL;
  ObitInfoType type;
  olong         Aseq, disk, idisk, cno;
  gchar        *Type, *strTemp, *strTemp2, outFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "  ";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "getOutputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input data Object */
  outData = newObitData("output Data");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4)))
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
    } else { /* Didn't find - use default */
      strncpy (Aclass, "noClas", 7);
      Aclass[0] = strTemp2[0];
    }
    /* If blank use noClas */
    if (!strncmp (Aclass, "    ", 4)) {
      strncpy (Aclass, "noClas", 7);
      Aclass[0] = strTemp2[0];
    }
    Aclass[6] = 0;

    /* output AIPS sequence */
    Aseq = 1;
    ObitInfoListGetTest(myInput, "outSeq", &type, dim, &Aseq);

    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outData);
    
    /* define object */
    ObitDataSetAIPS (outData, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outData);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */
    /* output FITS file name */
    ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp2);
    if (ObitInfoListGetP(myInput, "outFile", &type, dim, (gpointer)&strTemp)) {
      strncpy (outFile, strTemp, 128);
    } else { 
      g_snprintf (outFile, 129, "Squish%s", strTemp2);
    }
    /* If blank use Squish+inFile */
    if (!strncmp (outFile, "    ", 4)) {
      g_snprintf (outFile, 129, "Squish%s", strTemp2);
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
    ObitDataSetFITS (outData, disk, outFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outData);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return outData;
  }

  /* Full instantiation - needed to initialize tableList */
  ObitDataOpen (outData, OBIT_IO_ReadWrite, err);
  ObitDataClose (outData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);

  return outData;
} /* end getOutputData */

/*----------------------------------------------------------------------- */
/*  Write History for TabCopy                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitData to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void TabCopyHistory (ObitInfoList* myInput, ObitData* inData, ObitData* outData,
		     ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "inTab", "inVer",  "nCopy", 
    "outFile",  "outDisk", "outName", "outClass", "outSeq", "outVer",  
    "KeyWord",  "KeyRange", "KeyString", 
    NULL};
  gchar *routine = "TabCopyHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitDataIsA(inData));
  g_assert (ObitDataIsA(outData));

   /* Do history */
  outHistory = newObitDataHistory ((ObitData*)outData, OBIT_IO_ReadWrite, err);

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

  outHistory = ObitHistoryUnref(outHistory); /* cleanup */
 
} /* end TabCopyHistory  */


/**
 * Routine Copy selected Tables from inData to outData
 * \param inData     Input Data object
 * \param outData    Output Data object
 * \param err        Error/message stack, returns if error.
 */
void TabCopyCopy (ObitInfoList* myInput, ObitData* inData,  ObitData* outData,  
		  ObitErr* err)
{
  ObitTable *inTab=NULL, *outTab=NULL;
  olong i, inVer, outVer, nCopy;
  gchar *TabType=NULL, *KeyWord=NULL, *KeyString=NULL;
  ofloat KeyRange[3];
  olong inHigh, outHigh, iV, oV;
  ObitInfoType type;
  gint32   dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitTableList *inTL, *outTL;
  gboolean want=TRUE;
  gchar  *routine = "TabCopyCopy";
  
  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  g_assert(ObitDataIsA(inData));
  g_assert(ObitDataIsA(outData));

  /* Parameters used locally */
  inVer = 0;
  ObitInfoListGetTest(myInput, "inVer",      &type, dim, &inVer);
  nCopy = 1;
  ObitInfoListGetTest(myInput, "nCopy",      &type, dim, &nCopy);
  outVer = 0;
  ObitInfoListGetTest(myInput, "outVer",     &type, dim, &outVer);
  KeyRange[0]=-1.0e20; KeyRange[1]=1.0e20; 
  ObitInfoListGetTest(myInput, "KeyRange",   &type, dim, KeyRange);
  ObitInfoListGetP(myInput, "inTab",      &type, dim, (gpointer)&TabType);
  ObitInfoListGetP(myInput, "KeyWord",    &type, dim, (gpointer)&KeyWord);
  ObitInfoListGetP(myInput, "KeyString",  &type, dim, (gpointer)&KeyString);
  ObitTrimTrail (TabType);  /* No trailing blanks */

  /* Check */
  Obit_return_if_fail((TabType!=NULL), err, "%s: no table type specified ", routine);

  /* Table lists for source and destination */
  inTL    = inData->tableList;
  outTL   = outData->tableList;
  inHigh  = ObitTableListGetHigh (inTL, TabType);
  outHigh = ObitTableListGetHigh (outTL, TabType);
  /* Check */
  Obit_return_if_fail((inHigh>0), err, "%s: no tables of  type %s found in %s ", 
		      routine, TabType, inData->name);

  if (inVer>0) iV = inVer;
  else iV = inHigh;
  if (outVer>0) oV = outVer;
  else oV = outHigh+1;

  /* Loop over tables */
  for (i=0; i<nCopy; i++) {

    /* Want this input table? */
    if (!ObitTableListGet(inTL, TabType, &iV, &inTab, err)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: Table %s %d not found on on %s", 
		     routine, TabType, iV, inData->name);
    } else {

      /* May have to instantiate table */
      if (inTab==NULL) {
	inTab = newObitDataTable (inData, OBIT_IO_ReadOnly,  TabType, &iV, err);
	if (err->error) Obit_traceback_msg (err, routine, inData->name);
      }
      want  = WantTable(inTab, KeyWord, KeyRange, KeyString, err);
      inTab = ObitTableUnref(inTab);  /* Unreference table */
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);

    if (!want) continue;

    /* Output Table must not exist */
    if (ObitTableListGet(outTL, TabType, &oV, &outTab, err)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: Table %s %d already exists on %s", 
		     routine, TabType, oV, outData->name);
    } else {
      outTab = ObitTableUnref(outTab);  /* Unreference table */
    }
    if (err->error) Obit_traceback_msg (err, routine, inData->name);

    /* Copy table */
    ObitDataCopyTable (inData, outData, TabType, &iV, &oV, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);

    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Copied table %s version %d to output  %d",
		   TabType, iV, oV);
    ObitErrLog(err);

    /* Update next tables */
    iV++;
    oV++;
  } /* End loop over tables */

} /* end TabCopyCopy */

/**
 * Determine if table wanted
 * \param inTab      Input Table
 * \param KeyWord    Header keyword to check, blank or NULL -> OK
 * \param KeyRange   Range of allowed values if keyword is numeric
 *                   if KeyRange[1]<KeyRange[0] then only values OUTSIDE
 *                   [KeyRange[0],KeyRange[1] are allowed.
 * \param KeyString  String to match if KeyWord String
 * \param err        Error/message stack, returns if error.
 * \return TRUE is table wanted, else FALSE
 */
/* Is table wanted? */
gboolean WantTable (ObitTable *inTab, gchar *KeyWord, ofloat KeyRange[2], 
		    gchar *KeyString, ObitErr* err)
{
  gboolean want = TRUE;
  gpointer P=NULL;
  gshort *shortP;
  olong   *intP;
  olong  *longP;
  oint   *ointP;
  ofloat *floatP;
  odouble *doubleP, test;
  gchar *tstring;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar *routine = "WantTable";

  if (err->error) return want;  /* previous error? */

  /* Anything to test? */
  if ((KeyWord==NULL) || !strncmp(KeyWord, "        ", 8)) return want;

  /* Instantiate */
  ObitTableFullInstantiate (inTab, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, inTab->name, FALSE);

  /* Look up KeyWord */
  if (ObitInfoListGetP(inTab->myDesc->info, KeyWord, &type, dim, &P)) {
    switch (type) { 
    case OBIT_short:
      shortP = (gshort*)P;
      test = *shortP;
      break;
    case OBIT_int:
      intP = (olong*)P;
      test = *intP;
      break;
    case OBIT_oint:
      ointP = (oint*)P;
      test = *ointP;
      break;
    case OBIT_long:
      longP = (olong*)P;
      test = *longP;
      break;
    case OBIT_float:
      floatP = (ofloat*)P;
      test = *floatP;
      break;
    case OBIT_double:
      doubleP = (odouble*)P;
      test = *doubleP;
      break;
    case OBIT_string:
      /* Do string comparison here */
      if (KeyString==NULL) return FALSE;
      ObitTrimTrail (KeyString);
      tstring = g_strndup (P, dim[0]);
      ObitTrimTrail (tstring);
      want = !strcmp(tstring, KeyString);
      g_free(tstring);
      return want;
    default:
      /* No test */
      return TRUE;
    }; /* end switch  */
  } else { /* not found */
    Obit_log_error(err, OBIT_Error, "%s: Keyword %s not found in %s", 
                   routine, KeyWord, inTab->name);
    return FALSE;
  }
  
  /* Numeric test */
  if (KeyRange[0]<KeyRange[1]) {  /* in range wanted **/
    want = (test>=KeyRange[0]) && (test<=KeyRange[1]);
  } else { /* out of range wanted */
    want = (test<=KeyRange[0]) || (test>=KeyRange[1]);
  }

  return want;
} /* end WantTable */
