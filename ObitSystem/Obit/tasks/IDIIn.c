/* $Id$  */
/* Read IDI format data                               */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2007                                               */
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

#include "ObitUV.h"
#include "ObitFITS.h"
#include "ObitSystem.h"
#include "ObitAIPSDir.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitUV.h"
#include "ObitTableSU.h"
#include "ObitTableAN.h"
#include "ObitTableFQ.h"
#include "ObitTableIDI_ANTENNA.h"
#include "ObitTableIDI_ARRAY_GEOMETRY.h"
#include "ObitTableIDI_FREQUENCY.h"
#include "ObitTableIDI_SOURCE.h"
#include "ObitTableIDI_UV_DATA.h"
#include "ObitHistory.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif

/* internal prototypes */
/* Get inputs */
ObitInfoList* IDIInin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Create output uvdata */
ObitUV* setOutputData (ObitInfoList *myInput, ObitErr *err);
/* Get file descriptor */
void GetHeader (ObitUV *outData, gchar *infile, ObitInfoList *myInput, 
		ObitErr *err);
/* Get data */
void GetData (ObitUV *outData, gchar *infile, ObitInfoList *myInput, 
	      ObitErr *err);
/* Get Antenna info */
void GetAntennaInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Get Frequency info */
void GetFrequencyInfo (ObitData *inData, ObitUV *outData, ObitErr *err);
/* Get Source info */
void GetSourceInfo (ObitData *inData, ObitUV *outData, gboolean isNew, 
		    ObitErr *err);
/* Read data */
void ProcessData (gchar *inscan, ofloat avgTime,
		  olong *ndetect, olong *ntime, ofloat *refDate,
		  gfloat** ATime, gfloat*** AData, gfloat** ACal, 
		  ObitErr *err);
/* Write history */
void IDIInHistory (ObitInfoList* myInput, ObitUV* outData, ObitErr* err);

/* Program globals */
gchar *pgmName = "IDIIn";       /* Program name */
gchar *infile  = "IDIIn.inp";   /* File with program inputs */
gchar *outfile = "IDIIn.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
gchar **FITSdirs=NULL; /* List of FITS data directories */
gchar DataRoot[128]; /* Root directory of input data */
odouble refMJD;   /* reference Julian date */
odouble integTime;/* Integration time in days */
ofloat *SourceID = NULL; /* Source number (1-rel) lookup table */
gboolean isNew=FALSE;  /* Is the output newly created */
olong  nchan=1;   /* Number of frequencies */
olong  nstok=1;   /* Number of Stokes */
olong  nIF=1;     /* Number of IFs */
ofloat deltaFreq; /* Channel increment */
ofloat refPixel;  /* Frequency reference pixel */
odouble refFrequency; /* reference frequency (Hz) */


int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*    Read IDI IDI  data to a UV dataset                                  */
/*----------------------------------------------------------------------- */
{
  olong  i, ierr=0;
  ObitSystem *mySystem= NULL;
  ObitUV *outData= NULL;
  ObitErr *err= NULL;
  gchar inscan[128];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  err = newObitErr();

  /* Startup - parse command line */
  ierr = 0;
  myInput = IDIInin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input File name */
  for (i=0; i<128; i++) inscan[i] = 0;
  ObitInfoListGet(myInput, "Scan", &type, dim, inscan, err);
  inscan[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(inscan);  /* Trim trailing blanks */

  /* Get input data file root */
  for (i=0; i<128; i++) DataRoot[i] = 0;
  ObitInfoListGet(myInput, "DataRoot", &type, dim, DataRoot, err);
  DataRoot[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(DataRoot);  /* Trim trailing blanks */

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);

  if (err->error) ierr = 1;  ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Create ObitUV for data */
  outData = setOutputData (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
   
  /* Get header info, array geometry, initialize output if necessary */
  GetHeader (outData, inscan, myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Tell what's going on */   
  Obit_log_error(err, OBIT_InfoErr, "Adding file %s to output", 
		 inscan);
  ObitErrLog(err);

  /* show any errors */
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
   
  /* Open output data */
  if ((ObitUVOpen (outData, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", outData->name);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* convert data  */
  GetData (outData, inscan, myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Close */
  if ((ObitUVClose (outData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing output file");
  
  /* History */
  IDIInHistory (myInput, outData, err);
  
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;
  
  /* Shutdown Obit */
 exit:
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  mySystem = ObitSystemShutdown (mySystem);
  
  /* cleanup */
  myInput = ObitInfoListUnref(myInput);   /* delete input list */
  myInput = ObitInfoListUnref(myOutput);  /* delete output list */
  outData = ObitUnref(outData);
  if (SourceID) g_free(SourceID);
 
  return ierr;
} /* end of main */

ObitInfoList* IDIInin (int argc, char **argv, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Parse control info from command line                                  */
/*   Input:                                                               */
/*      argc   Number of arguments from command line                      */
/*      argv   Array of strings from command line                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   return  parser list                                                  */
/*----------------------------------------------------------------------- */
{
  olong i, j, k, ax;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar *input_file="IDIIn.in", *arg;
  gboolean init=FALSE;
  oint itemp;
  gchar *strTemp;
  ObitInfoList* list;
  gchar *routine = "IDIInin";

  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* command line arguments */
  if (argc<=1) Usage(); /* must have arguments */
  /* parse command line */
  for (ax=1; ax<argc; ax++) {
    arg = argv[ax];
    if (strcmp(arg, "-input") == 0){ /* input parameters */
      input_file = argv[++ax];
      /* parse input file */
      ObitParserParse (input_file, list, err);
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
      
    } else if (strcmp(arg, "-Scan") == 0){ /* Scan name */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "Scan", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-DataRoot") == 0){ /* Data root directory */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "DataRoot", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outFile") == 0){ /* Output FITS file */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outFile", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outDisk") == 0) { /* Output disk */
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

    } else if (strcmp(arg, "-outSeq") == 0) { /* AIPS output UV sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outSeq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-pgmNumber") == 0) { /*Program number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "pgmNumber", OBIT_oint, dim, &itemp, err);
      
    } else { /* unknown argument */
      /* DEBUG fprintf (stderr,"DEBUG parameter %s \n",arg);*/
      Usage();
    }
  } /* end parsing input arguments */
  
  /* Read defaults if no file specified */
  if (!init) ObitParserParse (input_file, list, err);

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
} /* end IDIInin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: IDIIn -input file -output ofile [args]\n");
    fprintf(stderr, "Convert an IDI/IDI file format to Obit/UV\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def IDIIn.in\n");
    fprintf(stderr, "  -output output result file, def UVSub.out\n");
    fprintf(stderr, "  -scan date/time used for form scan FITS file names\n");
    fprintf(stderr, "  -DataRoot Directory name for input \n");
    fprintf(stderr, "  -outFile output uv FITS  file\n");  
    fprintf(stderr, "  -outName output uv AIPS file name\n");
    fprintf(stderr, "  -outClass output uv AIPS file class\n");
    fprintf(stderr, "  -outSeq output uv AIPS file sequence\n");
    fprintf(stderr, "  -outDisk output uv (AIPS or FITS) disk number (1-rel) \n");
    
    /*/exit(1);  bail out */
  }/* end Usage */

/*----------------------------------------------------------------------- */
/*  Create default input ObitInfoList                                     */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
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
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  gchar *scan_name ="unspecified";
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
  itemp = 2; /* number of FITS directories */
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

  /* base of scan file names */
  dim[0] = strlen (scan_name);
  ObitInfoListPut (out, "Scan", OBIT_string, dim, scan_name, err);

  /* output FITS file name */
  strTemp = "uvData.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);

  /* root of data directory */
  strTemp = "dataRoot";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "DataRoot", OBIT_string, dim, strTemp, err);

  return out;
} /* end defaultInputs */

/*----------------------------------------------------------------------- */
/*  Create default output ObitInfoList                                    */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultOutputs(ObitErr *err)
{
  /*gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};*/
  /*gfloat ftemp;*/
  ObitInfoList *out = newObitInfoList();
  /*gchar *routine = "defaultOutputs";*/

  /* No outputs */
  return out;
} /* end defaultOutputs */

/*----------------------------------------------------------------------- */
/*  Create output uv data                                                 */
/*   Input:                                                               */
/*      Source    Source name                                             */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/* Returns the output uv data                                             */
/*----------------------------------------------------------------------- */
ObitUV* setOutputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitUV    *outUV = NULL;
  ObitInfoType type;
  olong      i, n, Aseq, disk, cno, lType;
  gchar     *Type, *strTemp, outFile[129];
  gchar     Aname[13], Aclass[7], *Atype = "UV";
  olong      nvis;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean  exist;
  gchar     tname[129], *fullname=NULL;
  gchar     *routine = "setOutputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic output UV Object */
  g_snprintf (tname, 100, "output UV data");
  outUV = newObitUV(tname);
    
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  lType = dim[0];
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */

    /* outName given? */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&strTemp);
    /* if not use inName */
    if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12)))
      ObitInfoListGetP (myInput, "inName", &type, dim, (gpointer)&strTemp);
    for (i=0; i<12; i++) Aname[i] = ' ';  Aname[i] = 0;
    for (i=0; i<MIN(12,dim[0]); i++) Aname[i] = strTemp[i];
    /* Save any defaulting on myInput */
    dim[0] = 12;
    ObitInfoListAlwaysPut (myInput, "outName", OBIT_string, dim, Aname);

      
    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    /* Default out class is "UVSub" */
    if (!strncmp(Aclass, "      ", 6)) strncpy (Aclass, "UVSub", 7);

    /* input AIPS disk - default is outDisk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (disk<=0)
       ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* output AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Allocate catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

    /* Did it previously exist? */
    isNew = !exist;
    
    /* define object */
    nvis = 1;
    ObitUVSetAIPS (outUV, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */

    /* outFile given? */
    ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&strTemp);
    /* if not use inName */
    if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12)))
      ObitInfoListGetP (myInput, "inFile", &type, dim, (gpointer)&strTemp);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) outFile[i] = strTemp[i]; outFile[i] = 0;
    ObitTrimTrail(outFile);  /* remove trailing blanks */

    /* Save any defaulting on myInput */
    dim[0] = strlen(outFile);
    ObitInfoListAlwaysPut (myInput, "outFile", OBIT_string, dim, outFile);

    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (disk<=0) /* defaults to outDisk */
      ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);

    /* Did it previously exist */
    disk = 1;  /* input "disk" */
    fullname = ObitFITSFilename (disk, outFile, err);
    exist =  ObitFileExist (fullname, err);
    if (fullname) g_free(fullname);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    isNew = !exist;
    
    /* define object */
    nvis = 1;
    ObitUVSetFITS (outUV, nvis, disk, outFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outUV;
  }
  
  return outUV;
} /* end setOutputUV */

void GetHeader (ObitUV *outData, gchar *inscan, 
		ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get header information from scan header files                         */
/*  Returns TRUE if the file is just created                              */
/*   Input:                                                               */
/*      outData  Output UV object                                         */
/*      inscan   root of input file names                                 */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitUVDesc *desc;
  ObitTableIDI_UV_DATA *inTable=NULL;
  ObitData *inData=NULL;
  olong ncol;
  gchar FullFile[128], *today=NULL;
  olong disk, lim;
  oint no_band=0, maxis1=0, maxis2=0, maxis3=0, maxis4=0, maxis5=0;
  ObitIOAccess access;
  olong ver;
  gchar *routine = "GetHeader";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outData));
  g_assert(inscan!=NULL);
  g_assert(myInput!=NULL);

  /* Full input file name */
  sprintf (FullFile,"%s/%s.fits", DataRoot, inscan);

  /* Create input Data from which to read tables */
  inData = newObitData("Input Data");
  disk = 0;
  ObitDataSetFITS(inData, disk, FullFile, err);
  /* Open and close to init TableList */
  ObitDataOpen (inData, OBIT_IO_ReadOnly, err);
  ObitDataClose (inData,  err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);


  /* Define output descriptor if isNew */
  if (isNew) {
    desc = outData->myDesc;
    /* Create input Antenna table object */
    ver = 1;
    access = OBIT_IO_ReadOnly;
    inTable = newObitTableIDI_UV_DATAValue ("Input table", inData, 
					    &ver, access, no_band, 
					    maxis1, maxis2, maxis3, maxis4, maxis5,
					    err);
   if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_UV_DATA table");
   if (err->error) Obit_traceback_msg (err, routine, inData->name);

    /* Open table */
    if ((ObitTableIDI_UV_DATAOpen (inTable, access, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_UV_DATA table");
      return;
    }

    /* Define header */
    desc->nvis = 0;
    strncpy (desc->origin, "Obit ", UVLEN_VALUE);
    desc->isort[0] = inTable->sort[0]; 
    desc->isort[1] = inTable->sort[1]; 
    desc->nrparm = 0;
    
    /* Creation date today */
    today = ObitToday();
    strncpy (desc->date, today, UVLEN_VALUE-1);
    if (today) g_free(today);
    desc->JDObs = 0.0;
    desc->epoch = 0.0;
    desc->equinox = 0.0;
    desc->obsra = 0.0;
    desc->obsdec = 0.0;
    desc->altCrpix = 0.0;
    desc->altRef = 0.0;
    desc->restFreq = 0.0;
    desc->xshift = 0.0;
    desc->yshift = 0.0;
    desc->VelDef = 0;
    desc->VelReference = 0;
    strncpy (desc->bunit, "        ", UVLEN_VALUE);
    lim = MIN (UVLEN_VALUE, MAXKEYCHARTABLEIDI_UV_DATA);
    strncpy (desc->obsdat,     inTable->RefDate, lim);
    strncpy (desc->teles,      inTable->teles, lim);
    strncpy (desc->observer,   inTable->observer, lim);
    strncpy (desc->instrument, inTable->ArrName, lim);
   
    /* Random parameters */
    ncol = 0;
    /* U */
    strncpy (desc->ptype[ncol], "UU-L-SIN", UVLEN_KEYWORD);
    ncol++;

    /* V */
    strncpy (desc->ptype[ncol], "VV-L-SIN", UVLEN_KEYWORD);
    ncol++;
    
    /* W */
    strncpy (desc->ptype[ncol], "WW-L-SIN", UVLEN_KEYWORD);
    ncol++;
    
    /* Baseline */
    strncpy (desc->ptype[ncol], "BASELINE", UVLEN_KEYWORD);
    ncol++;
    
    /* Time */
    strncpy (desc->ptype[ncol], "TIME1   ", UVLEN_KEYWORD);
    ncol++;
    
    /* Source */
    if (inTable->myDesc->repeat[inTable->SourceCol]>0) {
      strncpy (desc->ptype[ncol], "SOURCE  ", UVLEN_KEYWORD);
      ncol++;
    }
    
    /* FreqID */
    if (inTable->myDesc->repeat[inTable->FreqIDCol]>0) {
      strncpy (desc->ptype[ncol], "FREQSEL ", UVLEN_KEYWORD);
      ncol++;
    }
    
    /* Integration time */
    if (inTable->myDesc->repeat[inTable->IntTimCol]>0) {
      strncpy (desc->ptype[ncol], "INTTIM  ", UVLEN_KEYWORD);
      ncol++;
    }
    desc->nrparm = ncol;  /* Number of random parameters */
    
    /* Data Matrix */
    ncol = 0;
    lim = MIN (UVLEN_KEYWORD, MAXKEYCHARTABLEIDI_UV_DATA);
    /* Dimension 1 */
    if (inTable->maxis>=1) {
      strncpy (desc->ctype[ncol], inTable->ctype1, lim);
      desc->inaxes[ncol] = MAX (inTable->maxis1, 3);
      desc->cdelt[ncol] = inTable->cdelt1;
      desc->crpix[ncol] = inTable->crpix1;
      desc->crval[ncol] = inTable->crval1;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    /* Dimension 2 */
    if (inTable->maxis>=2) {
      strncpy (desc->ctype[ncol], inTable->ctype2, lim);
      desc->inaxes[ncol] = inTable->maxis2;
      desc->cdelt[ncol] = inTable->cdelt2;
      desc->crpix[ncol] = inTable->crpix2;
      desc->crval[ncol] = inTable->crval2;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    /* Dimension 3 */
    if (inTable->maxis>=3) {
      strncpy (desc->ctype[ncol], inTable->ctype3, lim);
      desc->inaxes[ncol] = inTable->maxis3;
      desc->cdelt[ncol] = inTable->cdelt3;
      desc->crpix[ncol] = inTable->crpix3;
      desc->crval[ncol] = inTable->crval3;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    /* Dimension 4 */
    if (inTable->maxis>=4) {
      strncpy (desc->ctype[ncol], inTable->ctype4, lim);
      desc->inaxes[ncol] = inTable->maxis4;
      desc->cdelt[ncol] = inTable->cdelt4;
      desc->crpix[ncol] = inTable->crpix4;
      desc->crval[ncol] = inTable->crval4;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    /* Dimension 5 */
    if (inTable->maxis>=5) {
      strncpy (desc->ctype[ncol], inTable->ctype5, lim);
      desc->inaxes[ncol] = inTable->maxis5;
      desc->cdelt[ncol] = inTable->cdelt5;
      desc->crpix[ncol] = inTable->crpix5;
      desc->crval[ncol] = inTable->crval5;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    /* Dimension 6 */
    if (inTable->maxis>=6) {
      strncpy (desc->ctype[ncol], inTable->ctype6, lim);
      desc->inaxes[ncol] = inTable->maxis6;
      desc->cdelt[ncol] = inTable->cdelt6;
      desc->crpix[ncol] = inTable->crpix6;
      desc->crval[ncol] = inTable->crval6;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    /* Dimension 7 */
    if (inTable->maxis>=7) {
      strncpy (desc->ctype[ncol], inTable->ctype7, lim);
      desc->inaxes[ncol] = inTable->maxis7;
      desc->cdelt[ncol] = inTable->cdelt7;
      desc->crpix[ncol] = inTable->crpix7;
      desc->crval[ncol] = inTable->crval7;
      desc->crota[ncol] = 0.0;
      ncol++;
    }
    desc->naxis = ncol;  /* Number of dimensions */

    /* Close  table */
    if ((ObitTableIDI_UV_DATAClose (inTable, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_UV_DATA Table file");
      return;
    }

    /* index descriptor */
    ObitUVDescIndex (desc);

    /* Add Antenna and Frequency info */
    ObitUVOpen (outData, OBIT_IO_WriteOnly, err) ;
    GetFrequencyInfo (inData, outData, err);
    GetAntennaInfo (inData, outData, err);
    ObitUVClose (outData, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);

  } /* End define descriptor */

  /* Get source info, copy to output SU table, save lookup table in SourceID */
  ObitUVOpen (outData, OBIT_IO_ReadWrite, err);
  GetSourceInfo (inData, outData, isNew, err);
  ObitUVClose (outData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
 
  /* Instantiate output Data */
  ObitUVFullInstantiate (outData, FALSE, err);
  if (err->error)Obit_traceback_msg (err, routine, outData->name);


  /* Cleanup */
  inData  = ObitDataUnref(inData);
  inTable = ObitTableIDI_UV_DATAUnref(inTable);

  return;
} /* end GetHeader */

void GetData (ObitUV *outData, gchar *inscan, ObitInfoList *myInput, 
	      ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Read data from IDI file, write outData                                */
/*      outData  Output OTF object, open on input                         */
/*      inscan   Scan part of input file name                             */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  olong iRow, kinc, lim, ver, ant1, ant2, tant1, tant2;
  olong disk, i;
  oint no_band=0, maxis1=0, maxis2=0, maxis3=0, maxis4=0, maxis5=0;
  ofloat lambdaPerSec, sid;
  ofloat *Buffer=NULL;
  odouble JD;
  ObitUVDesc *desc;
  ObitIOAccess access;
  ObitData *inData=NULL;
  ObitTableIDI_UV_DATA *inTable=NULL;
  ObitTableIDI_UV_DATARow *inRow=NULL;
  gchar FullFile[128];
  gchar *routine = "GetData";

  /* error checks */
  if (err->error) return;
  g_assert(inscan!=NULL);
  g_assert(myInput!=NULL);
  g_assert(ObitUVIsA(outData));

  /* Full input file name */
  sprintf (FullFile,"%s/%s.fits", DataRoot, inscan);

  /* Create input Data from which to read tables */
  inData = newObitData("Input Data");
  disk = 0;
  ObitDataSetFITS(inData, disk, FullFile, err);
  /* Open and close to init TableList */
  ObitDataOpen (inData, OBIT_IO_ReadOnly, err);
  ObitDataClose (inData,  err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  ver = 1;
  access = OBIT_IO_ReadOnly;
  inTable = newObitTableIDI_UV_DATAValue ("Input table", inData, 
					  &ver, access, no_band, 
					  maxis1, maxis2, maxis3, maxis4, maxis5,
					  err);
  if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_UV_DATA table");
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Open table */
  if ((ObitTableIDI_UV_DATAOpen (inTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_UV_DATA table");
    return;
  }

  /* Create Row */
  inRow = newObitTableIDI_UV_DATARow (inTable);

  /* Consistency checks */
  desc = outData->myDesc;
  if (desc->jlocs>=0) nstok = desc->inaxes[desc->jlocs];
  if (desc->jlocf>=0) nchan = desc->inaxes[desc->jlocf];
  if (desc->jlocif>=0) nIF  = desc->inaxes[desc->jlocif];

  Obit_return_if_fail((nchan==inTable->no_chan), err,
		       "%s: Input number freq. incompatible %d != %d", 
		      routine, nchan, inTable->no_chan);
  Obit_return_if_fail((nstok==inTable->no_stkd), err,
		       "%s: Input number Poln incompatible %d != %d", 
		      routine, nstok, inTable->no_stkd);
  Obit_return_if_fail((nIF==inTable->no_band), err,
		       "%s: Input number Bands (IFs) incompatible %d != %d", 
		      routine, nIF, inTable->no_band);
  lim = MIN (UVLEN_KEYWORD, MAXKEYCHARTABLEIDI_UV_DATA);
  Obit_return_if_fail((!strncmp(inTable->ctype1,desc->ctype[0],lim)), err,
		       "%s: First axis different %s != %s", 
		      routine, inTable->ctype1,desc->ctype[0]);  
  Obit_return_if_fail((!strncmp(inTable->ctype2,desc->ctype[1],lim)), err,
		      "%s: Second axis different %s != %s", 
		      routine, inTable->ctype2,desc->ctype[1]);
  if (inTable->maxis>=2) {
    Obit_return_if_fail((!strncmp(inTable->ctype3,desc->ctype[2],lim)), err,
			"%s: Second axis different %s != %s", 
			routine, inTable->ctype3,desc->ctype[2]);
  }
  if (inTable->maxis>=3) {
    Obit_return_if_fail((!strncmp(inTable->ctype3,desc->ctype[2],lim)), err,
			"%s: Third axis different %s != %s", 
			routine, inTable->ctype3,desc->ctype[2]);
  }
  if (inTable->maxis>=4) {
    Obit_return_if_fail((!strncmp(inTable->ctype4,desc->ctype[3],lim)), err,
			"%s: Fourth axis different %s != %s", 
			routine, inTable->ctype4,desc->ctype[3]);
  }
  if (inTable->maxis>=5) {
    Obit_return_if_fail((!strncmp(inTable->ctype5,desc->ctype[4],lim)), err,
			"%s: Fifth axis different %s != %s", 
			routine, inTable->ctype5,desc->ctype[4]);
  }
  if (inTable->maxis>=6) {
    Obit_return_if_fail((!strncmp(inTable->ctype6,desc->ctype[5],lim)), err,
			"%s: Sixth axis different %s != %s", 
			routine, inTable->ctype6,desc->ctype[5]);
  }
  Obit_return_if_fail((desc->ncorr==inTable->myDesc->repeat[inTable->WeightCol]), err,
		      "%s: Input and output data sizes different, %d !=  %d", 
		      routine, desc->ncorr, 
		      inTable->myDesc->repeat[inTable->WeightCol]);

  /* Prepare output */
  desc = outData->myDesc;
  Buffer = outData->buffer;
  desc->firstVis = desc->nvis+1;  /* Append to end of data */
  
  /* Get reference MJD, convert ref date  to JD */
  ObitUVDescDate2JD (inTable->RefDate, &JD);
  /* Reference Modified Julian date */
  refMJD = JD - 2400000.5;

  /* Number of values on axis 1 */
  kinc = inTable->maxis1;

  /* Number of wavelengths per second */
  lambdaPerSec = inTable->ref_freq;

  /* Loop over table */
  for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {

    if ((ObitTableIDI_UV_DATAReadRow (inTable, iRow, inRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR reading IDI_UV_DATA Table");
      return;
      }
    if (inRow->status==-1) continue;  /* Flagged entry */

    /* Antenna numbers in proper order */
    tant1 = (olong)(inRow->Baseline/256);
    tant2 = (olong)(inRow->Baseline - tant1*256);
    ant1 = MIN (tant1, tant2);
    ant2 = MAX (tant1, tant2);
 
    /* Convert table to UV data form */
    /* Convert U,V,W to lambda */
    Buffer[desc->ilocu] = (ofloat)inRow->uu*lambdaPerSec;
    Buffer[desc->ilocv] = (ofloat)inRow->vv*lambdaPerSec;
    Buffer[desc->ilocw] = (ofloat)inRow->ww*lambdaPerSec;
    Buffer[desc->ilocb] = (ofloat)(ant1*256+ant2);
    Buffer[desc->iloct] = (ofloat)(inRow->date+inRow->Time-JD);
    if (desc->ilocsu>=0) {
      sid = SourceID[inRow->Source-1];
      Buffer[desc->ilocsu] = sid;
    }
    if (desc->ilocfq>=0) Buffer[desc->ilocfq] = (ofloat)inRow->FreqID;
    if (desc->ilocit>=0) Buffer[desc->ilocit] = (ofloat)inRow->IntTim;
    for (i=0; i<desc->ncorr; i++) {
      Buffer[desc->nrparm+i*3]   = inRow->Flux[i*kinc];
      Buffer[desc->nrparm+i*3+1] = inRow->Flux[i*kinc+1];
      Buffer[desc->nrparm+i*3+2] = inRow->Weight[i];
    }

    /* set number of records */
    desc->numVisBuff = 1;
    
    /* Write output */
    if ((ObitUVWrite (outData, NULL, err) != OBIT_IO_OK) || (err->error))
      Obit_log_error(err, OBIT_Error, "ERROR writing output UV data"); 
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
  } /* end loop over table */  


  /* Close  tables */
  if ((ObitTableIDI_UV_DATAClose (inTable, err) 
       != OBIT_IO_OK) || (err->error)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_UV_DATA Table");
    return;
  }

  /* Cleanup */
  inData  = ObitDataUnref(inData);
  inTable = ObitTableIDI_UV_DATAUnref(inTable);
  inRow   = ObitTableIDI_UV_DATAUnref(inTable);

} /* end GetData  */

void GetAntennaInfo (ObitData *inData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get info from IDI_ANTENNA and IDI_ARRAY_GEOMETRY tables on indata     */
/*  Writes AN table output                                                */
/*  to output .                                                           */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_ANTENNA         *inTable=NULL;
  ObitTableIDI_ANTENNARow      *inRow=NULL;
  ObitTableIDI_ARRAY_GEOMETRY  *in2Table=NULL;
  ObitTableIDI_ARRAY_GEOMETRYRow  *in2Row=NULL;
  ObitTableAN                  *outTable=NULL;
  ObitTableANRow               *outRow=NULL;
  olong i, lim, iRow,i2Row, oRow, ver;
  oint numIF, numPCal, numOrb;
  gboolean found;
  ObitIOAccess access;
  gchar *routine = "GetAntennaInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Create input Antenna table object */
  ver = 1;
  access = OBIT_IO_ReadOnly;
  numIF =  numPCal = 0;
  inTable = newObitTableIDI_ANTENNAValue ("Input table", inData, 
					  &ver, access, numIF, numPCal, err);
  if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_ANTENNA table");
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Open table */
  if ((ObitTableIDI_ANTENNAOpen (inTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_ANTENNA table");
    return;
  }

  /* Create Row */
  inRow = newObitTableIDI_ANTENNARow (inTable);

  /* Create input Array geometry table object */
  ver = 1;
  access = OBIT_IO_ReadOnly;
  numIF = numOrb = 0;
  in2Table = newObitTableIDI_ARRAY_GEOMETRYValue ("Input table", inData, 
						  &ver, access, numIF, numOrb, err);
  if (in2Table==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_ARRAY_GEOMETRY table");
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Open table */
  if ((ObitTableIDI_ARRAY_GEOMETRYOpen (in2Table, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_ARRAY_GEOMETRY table");
    return;
  }

  /* Create Row */
  in2Row = newObitTableIDI_ARRAY_GEOMETRYRow (in2Table);

  /* Create output Antenna table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  numOrb  = in2Table->numOrb;
  numPCal = inTable->numPCal;
  outTable = newObitTableANValue ("Output table", (ObitData*)outData, 
				  &ver, access, numOrb, numPCal, err);
  if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open table */
  if ((ObitTableANOpen (outTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output AN table");
    return;
  }

  /* Create output Row */
  outRow = newObitTableANRow (outTable);
  /* attach to table buffer */
  ObitTableANSetRow (outTable, outRow, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Set AN table values */
  outTable->ArrayX = in2Table->ArrayX;
  outTable->ArrayY = in2Table->ArrayY;
  outTable->ArrayZ = in2Table->ArrayZ;
  outTable->GSTiat0 = in2Table->GSTiat0;
  outTable->DegDay  = in2Table->DegDay;
  outTable->Freq    = in2Table->ref_freq;
  lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, MAXKEYCHARTABLEAN);
  strncpy (outTable->RefDate, inTable->RefDate, lim);
  outTable->PolarX  = in2Table->PolarX;
  outTable->PolarY  = in2Table->PolarY;
  outTable->dataUtc = in2Table->ut1Utc;
  lim = MIN (MAXKEYCHARTABLEIDI_ARRAY_GEOMETRY, MAXKEYCHARTABLEAN);
  strncpy (outTable->TimeSys, in2Table->TimeSys, lim );
  outTable->FreqID = 0;
  outTable->iatUtc = in2Table->ut1Utc;
  lim = MIN (MAXKEYCHARTABLEIDI_ANTENNA, MAXKEYCHARTABLEAN);
  strncpy (outTable->polType, inTable->polType, lim );
  outTable->P_Refant = 0;
  outTable->P_Diff01 = 0.0;
  outTable->P_Diff02 = 0.0;
  outTable->P_Diff03 = 0.0;
  outTable->P_Diff04 = 0.0;
  outTable->P_Diff05 = 0.0;
  outTable->P_Diff06 = 0.0;
  outTable->P_Diff07 = 0.0;
  outTable->P_Diff08 = 0.0;

  /* Initialize output row */
  outRow->noSta     = 0;
  outRow->mntSta    = 0;
  outRow->staXof    = 0.0;
  outRow->PolAngA   = 0.0;
  outRow->PolAngB   = 0.0;
  outRow->AntName[0]= 0; 
  outRow->StaXYZ[0] = 0.0;
  outRow->StaXYZ[1 ]= 0.0;
  outRow->StaXYZ[2] = 0.0;
  outRow->OrbParm[0]= 0.0;
  outRow->polTypeA[0] = ' ';
  outRow->PolCalA[0]  = 0.0;
  outRow->polTypeB[0] =  ' ';
  outRow->PolCalB[0]  = 0.0;
  outRow->status      = 0;

  /* loop through input Antenna table */
  for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
    if ((ObitTableIDI_ANTENNAReadRow (inTable, iRow, inRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR reading IDI_ANTENNA Table");
      return;
    }

    /* Find in Array Geometry table */
    found = FALSE;
    for (i2Row = 1; i2Row<=in2Table->myDesc->nrow; i2Row++) {
      if ((ObitTableIDI_ARRAY_GEOMETRYReadRow (in2Table, i2Row, in2Row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading IDI_ARRAY_GEOMETRYR Table");
	return;
      }
    if (inRow->antenna_no==in2Row->noSta) {
      /* Found match */
      found = TRUE;
      break;
      }  
    } /* end loop over table */

    /* Set output Row */
    outRow->noSta     = inRow->antenna_no;
    outRow->PolAngA   = inRow->PolAngA;
    outRow->polTypeA  = inRow->polTypeA;   
    for (i=0; i<numPCal; i++) 
	outRow->PolCalA[i] = inRow->PolCalA[i];
    if (inTable->no_stkd>1) {  /* Multiple polarizations */
      outRow->PolAngB   = inRow->PolAngB;   
      outRow->polTypeB  = inRow->polTypeB;   
      for (i=0; i<numPCal; i++) 
	outRow->PolCalB[i] = inRow->PolCalB[i];
    }
    lim = MIN(inTable->myDesc->repeat[inTable->AntNameCol], 
	      outTable->myDesc->repeat[outTable->AntNameCol]);
    for (i=0; i<lim; i++) 
	outRow->AntName[i] = inRow->AntName[i];
    outRow->status      = 0;
    if (found) {
      outRow->mntSta    = in2Row->mntSta;
      outRow->staXof    = in2Row->staXof[0];
      outRow->StaXYZ[0] = in2Row->StaXYZ[0];
      outRow->StaXYZ[1] = in2Row->StaXYZ[1];
      outRow->StaXYZ[2] = in2Row->StaXYZ[2];
      for (i=0; i<numOrb; i++) outRow->OrbParm[i] = in2Row->OrbParm[i];
    }
    oRow = outRow->noSta;
    if ((ObitTableANWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating ANTENNA Table");
      return;
    }
  } /* end loop over input table */
  
  /* Close  tables */
  if ((ObitTableIDI_ANTENNAClose (inTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_ANTENNA Table file");
    return;
  }

  if ((ObitTableIDI_ARRAY_GEOMETRYClose (in2Table, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_ARRAY_GEOMETRY Table file");
    return;
  }

  if ((ObitTableANClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output Antenna Table file");
    return;
  }

  /* Cleanup */
  inRow    = ObitTableIDI_ANTENNARowUnref(inRow);
  inTable  = ObitTableIDI_ANTENNAUnref(inTable);
  in2Row   = ObitTableIDI_ARRAY_GEOMETRYRowUnref(in2Row);
  in2Table = ObitTableIDI_ARRAY_GEOMETRYUnref(in2Table);
  outRow   = ObitTableANRowUnref(outRow);
  outTable = ObitTableANUnref(outTable);

} /* end  GetAntennaInfo */

void GetFrequencyInfo (ObitData *inData, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get Source info from IDI_FREQUENCY table on indata                    */
/*  Write FQ table                                                        */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_FREQUENCY*    inTable=NULL;
  ObitTableIDI_FREQUENCYRow* inRow=NULL;
  ObitTableFQ*            outTable=NULL;
  ObitTableFQRow*         outRow=NULL;
  olong i, iRow, oRow, ver;
  oint numIF;
  ObitIOAccess access;
  gchar *routine = "GetFrequencyInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Create input Source table object */
  ver = 1;
  access = OBIT_IO_ReadOnly;
  numIF = 0;
  inTable = newObitTableIDI_FREQUENCYValue ("Input table", inData, 
					 &ver, access, numIF, err);
  if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_FREQUENCY table");
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Open table */
  if ((ObitTableIDI_FREQUENCYOpen (inTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_FREQUENCY table");
    return;
  }

  /* Create Row */
  inRow = newObitTableIDI_FREQUENCYRow (inTable);

  /* Create output FQ table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  numIF = inTable->no_band;
  outTable = newObitTableFQValue ("Output table", (ObitData*)outData, 
				  &ver, access, numIF, err);
  if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with FQ table");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open table */
  if ((ObitTableFQOpen (outTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FQ table");
    return;
  }

  /* Create output Row */
  outRow = newObitTableFQRow (outTable);
  /* attach to table buffer */
  ObitTableFQSetRow (outTable, outRow, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Initialize output row */
  outRow->fqid    = 0;
  for (i=0; i<numIF; i++) {
    outRow->freqOff[i]  = 0;
    outRow->chWidth[i]  = 0;
    outRow->totBW[i]    = 0;
    outRow->sideBand[i] = 0;
  }
  outRow->status    = 0;

  /* loop through input table */
  for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
    if ((ObitTableIDI_FREQUENCYReadRow (inTable, iRow, inRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR reading IDI_FREQUENCY Table");
      return;
    }

    /* Save to FQ table */
    outRow->fqid    = inRow->fqid;
    for (i=0; i<numIF; i++) {
      outRow->freqOff[i]  = inRow->bandfreq[i];
      outRow->chWidth[i]  = inRow->chWidth[i];
      outRow->totBW[i]    = inRow->totBW[i];
      outRow->sideBand[i] = inRow->sideBand[i];
    }
    outRow->status    = 0;
    oRow = outRow->fqid;
    if (oRow<1) oRow = -1;
    if ((ObitTableFQWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating Source Table");
      return;
    }
  } /* end loop over input table */
  
  /* Close  tables */
  if ((ObitTableIDI_FREQUENCYClose (inTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_FREQUENCY Table file");
    return;
  }

  if ((ObitTableFQClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output Source Table file");
    return;
  }

  /* Cleanup */
  inRow    = ObitTableIDI_FREQUENCYRowUnref(inRow);
  inTable  = ObitTableIDI_FREQUENCYUnref(inTable);
  outRow   = ObitTableFQRowUnref(outRow);
  outTable = ObitTableFQUnref(outTable);

} /* end  GetFrequencyInfo */

void GetSourceInfo (ObitData *inData, ObitUV *outData, gboolean isNew, 
		    ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get Source info from IDI_SOURCE table on indata                       */
/*  Copies source info from intput to putput data and generated a lookup  */
/*  table, global, SourceID to give translation from input source ID      */
/*  to output .                                                           */
/*   Input:                                                               */
/*      inData   Input IDI FITS object                                    */
/*      outData  Output UV object                                         */
/*      isNew    True if output file just created                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableIDI_SOURCE*    inTable=NULL;
  ObitTableIDI_SOURCERow* inRow=NULL;
  ObitTableSU*            outTable=NULL;
  ObitTableSURow*         outRow=NULL;
  olong i, lim, iRow, oRow, ver;
  oint numIF;
  gboolean found;
  ObitIOAccess access;
  gchar *routine = "GetSourceInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Create input Source table object */
  ver = 1;
  access = OBIT_IO_ReadOnly;
  numIF = 0;
  inTable = newObitTableIDI_SOURCEValue ("Input table", inData, 
					 &ver, access, numIF, err);
  if (inTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with IDI_SOURCE table");
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Open table */
  if ((ObitTableIDI_SOURCEOpen (inTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input IDI_SOURCE table");
    return;
  }

  /* Create Row */
  inRow = newObitTableIDI_SOURCERow (inTable);

  /* Create global source ID lookup table */
  SourceID = g_malloc0((2+inTable->myDesc->nrow)*sizeof(ofloat));

  /* Create output Source table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  numIF = inTable->no_band;
  outTable = newObitTableSUValue ("Output table", (ObitData*)outData, 
				  &ver, access, numIF, err);
  if (outTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with SU table");
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open table */
  if ((ObitTableSUOpen (outTable, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output SU table");
    return;
  }

  /* Create output Row */
  outRow = newObitTableSURow (outTable);
  /* attach to table buffer */
  ObitTableSUSetRow (outTable, outRow, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Initialize output row */
  outRow->SourID    = 0;
  outRow->Qual      = 0;
  outRow->Bandwidth = 0.0;
  outRow->RAMean    = 0.0;
  outRow->DecMean   = 0.0;
  outRow->Epoch     = 0.0;
  outRow->RAApp     = 0.0;
  outRow->DecApp    = 0.0;
  outRow->PMRa      = 0.0;
  outRow->PMDec     = 0.0;
  outRow->Source[0] = 0;
  outRow->CalCode[0]= 0;
  for (i=0; i<inTable->no_band; i++) {
    outRow->IFlux[i]     = 0.0;
    outRow->QFlux[i]     = 0.0;
    outRow->UFlux[i]     = 0.0;
    outRow->VFlux[i]     = 0.0;
    outRow->FreqOff[i]   = 0.0;
    outRow->LSRVel[i]    = 0.0;
    outRow->RestFreq [i] = 0.0;
  }
  outRow->status    = 0;


  /* loop through input table */
  for (iRow = 1; iRow<=inTable->myDesc->nrow; iRow++) {
    if ((ObitTableIDI_SOURCEReadRow (inTable, iRow, inRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR reading IDI_SOURCE Table");
      return;
    }

    /* See if source exists in output table */
    found = FALSE;
    for (oRow = 1; oRow<=outTable->myDesc->nrow; oRow++) {
      if ((ObitTableSUReadRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading Source Table");
	return;
      }
    if (!strncmp (inRow->Source, outRow->Source, 16)) {
      /* Found match */
      found = TRUE;
      break;
      }  
    } /* end loop over table */

    /* If found just remember output source ID */
    if (found) {
      /* Save source no. in globla lookup table */
      SourceID[inRow->SourID] = outRow->SourID;
    } else { /* Not found - add */
      /* If first time update header */
      if (isNew) {
	outTable->FreqID = inRow->FreqID;
	lim = MIN(inTable->myDesc->repeat[inTable->SysVelCol], MAXKEYCHARTABLESU);
	for (i=0; i<lim; i++) outTable->velDef[i] = inRow->SysVel[i];
	lim = MIN(inTable->myDesc->repeat[inTable->VelTypCol], MAXKEYCHARTABLESU);
	for (i=0; i<lim; i++) outTable->velType[i] = inRow->VelTyp[i];
      }
      /* Set output row for end of table */
      outRow->SourID    = outTable->myDesc->nrow+1;
      outRow->Qual      = inRow->Qual;
      outRow->RAMean    = inRow->RAMean;
      outRow->DecMean   = inRow->DecMean;
      outRow->Epoch     = inRow->Epoch;
      outRow->RAApp     = inRow->RAApp;
      outRow->DecApp    = inRow->DecApp;
      outRow->PMRa      = inRow->PMRa;
      outRow->PMDec     = inRow->PMDec;
      lim = MIN(inTable->myDesc->repeat[inTable->SourceCol], 
		outTable->myDesc->repeat[outTable->SourceCol]);
      for (i=0; i<lim; i++) 
	outRow->Source[i] = inRow->Source[i];
      lim = MIN(inTable->myDesc->repeat[inTable->CalCodeCol], 
		outTable->myDesc->repeat[outTable->CalCodeCol]);
      for (i=0; i<lim; i++) 
	outRow->CalCode[i]= inRow->CalCode[i];
      for (i=0; i<inTable->no_band; i++) {
	outRow->IFlux[i]    = inRow->IFlux[i];
	outRow->QFlux[i]    = inRow->QFlux[i];
	outRow->UFlux[i]    = inRow->UFlux[i];
	outRow->VFlux[i]    = inRow->VFlux[i];
	outRow->FreqOff[i]  = inRow->FreqOff[i];
	outRow->LSRVel[i]   = inRow->SysVel[i];
	outRow->RestFreq[i] = inRow->RestFreq[i];
      }
      outRow->status    = 0;
 
      /* Save source no. in global lookup table */
      SourceID[inRow->SourID] = outRow->SourID;

      oRow = outRow->SourID;
      if ((ObitTableSUWriteRow (outTable, oRow, outRow, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR updating Source Table");
	return;
      }
   } /* End add new entry */
 } /* end loop over input table */
  
  /* Close  tables */
  if ((ObitTableIDI_SOURCEClose (inTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input IDI_SOURCE Table file");
    return;
  }

  if ((ObitTableSUClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output Source Table file");
    return;
  }

  /* Cleanup */
  inRow    = ObitTableIDI_SOURCERowUnref(inRow);
  inTable  = ObitTableIDI_SOURCEUnref(inTable);
  outRow   = ObitTableSURowUnref(outRow);
  outTable = ObitTableSUUnref(outTable);

} /* end  GetSourceInfo */

/*----------------------------------------------------------------------- */
/*  Write History for IDIIn                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      outData   ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void IDIInHistory (ObitInfoList* myInput, ObitUV* outData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "Scan", 
    NULL};
  gchar *routine = "IDIInHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(outData));

  /* Do history  */
  outHistory = newObitDataHistory ((ObitData*)outData, OBIT_IO_WriteOnly, err);

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

  outHistory = ObitHistoryUnref(outHistory);
 
} /* end IDIInHistory  */
