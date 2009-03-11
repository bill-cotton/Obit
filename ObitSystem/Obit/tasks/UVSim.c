/* $Id:  $  */
/* Simulate UV data                                                   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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
#include "ObitUVUtil.h"
#include "ObitTableSU.h"
#include "ObitTableAN.h"
#include "ObitTableANUtil.h"
#include "ObitTableFQ.h"
#include "ObitPrecess.h"
#include "ObitAntennaList.h"
#include "ObitHistory.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif

/* internal prototypes */
/* Get inputs */
ObitInfoList* UVSimin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Create output uvdata */
ObitUV* setOutputData (ObitInfoList *myInput, ObitErr *err);
/* Get file descriptor */
void GetHeader (ObitUV *outData, ObitInfoList *myInput, ObitErr *err);
/* Get data */
void GetData (ObitUV *outData, ObitInfoList *myInput, ObitErr *err);
/* Get Antenna info */
void GetAntennaInfo (ObitInfoList *myInput, ObitUV *outData, ObitErr *err);
/* Get Frequency info */
void GetFrequencyInfo (ObitInfoList *myInput, ObitUV *outData, ObitErr *err);
/* Get Source info */
void GetSourceInfo (ObitInfoList *myInput, ObitUV *outData, gboolean isNew, 
		    ObitErr *err);
/* Write history */
void UVSimHistory (ObitInfoList* myInput, ObitUV* outData, ObitErr* err);

/* Program globals */
gchar *pgmName = "UVSim";       /* Program name */
gchar *infile  = "UVSim.inp";   /* File with program inputs */
gchar *outfile = "UVSim.out";   /* File to contain program outputs */
olong  pgmNumber;           /* Program number (like POPS no.) */
olong  AIPSuser;            /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;             /* Number of AIPS directories */
gchar **AIPSdirs=NULL;      /* List of AIPS data directories */
olong  nFITS=0;             /* Number of FITS directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
gchar **FITSdirs=NULL;         /* List of FITS data directories */
odouble refJD;         /* reference Julian date */
odouble integTime;     /* Integration time in days */
ofloat SourceID = 1.0; /* Source number (1-rel) lookup table */
gboolean isNew=FALSE;  /* Is the output newly created */
olong  nchan=1;        /* Number of frequencies */
olong  nstok=4;        /* Number of Stokes */
olong  nIF  =1;        /* Number of IFs */
olong  nAnt =1;        /* Number of IFs */
ofloat deltaFreq;      /* Channel increment */
ofloat deltaIF;        /* IF Freq increment */
ofloat refPixel=1.0;   /* Frequency reference pixel */
odouble refFrequency;  /* reference frequency (Hz) */
odouble RA, Dec;       /* Mean position (deg) */
odouble RAApp, DecApp; /* Apparent position (deg) */
odouble arrayXYZ[3];   /* Array center, earth centered coordinates */
odouble *antXYZ=NULL;  /* Array coordinates (x,y,z in m) wrt arrayXYZ, 
			  VLA convention */
odouble GSTiat0;       /* GST at iat 0 */
odouble DegDay;        /* Earth rotation rate in deg per day */
ObitAntennaList *AntList=NULL; /* Antenna List */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*    Simulate a UV dataset                                               */
/*----------------------------------------------------------------------- */
{
  olong        ierr   = 0;
  ObitSystem   *mySystem = NULL;
  ObitUV       *outData  = NULL;
  ObitErr      *err      = NULL;

  err = newObitErr();

  /* Startup - parse command line */
  ierr = 0;
  myInput = UVSimin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;  ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Create ObitUV for data */
  outData = setOutputData (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
   
  /* Get header info, array geometry, initialize output if necessary */
  GetHeader (outData, myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Open output data */
  if ((ObitUVOpen (outData, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", 
		   outData->name);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Simulate data  */
  GetData (outData, myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Close */
  if ((ObitUVClose (outData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing output file");
  
  /* History */
  UVSimHistory (myInput, outData, err);
  
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
 
  return ierr;
} /* end of main */

ObitInfoList* UVSimin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="UVSim.in", *arg;
  gboolean init=FALSE;
  oint itemp;
  gchar *strTemp;
  ObitInfoList* list;
  gchar *routine = "UVSimin";

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
} /* end UVSimin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: UVSim -input file -output ofile [args]\n");
    fprintf(stderr, "Simulate UV data and write Obit/UV\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def UVSim.in\n");
    fprintf(stderr, "  -output output result file, def UVSim.out\n");
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
/*     outDType  Str [4]    "AIPS" or "FITS" [def {"FITS"}]               */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
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
  ObitInfoListPut (out, "ODType", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output FITS file name */
  strTemp = "uvData.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);

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
  /*ofloat ftemp;*/
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
    /* Default out class is "UVSim" */
    if (!strncmp(Aclass, "      ", 6)) strncpy (Aclass, "UVSim", 7);

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

void GetHeader (ObitUV *outData, ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get header information from myInput                                   */
/*  Returns TRUE if the file is just created                              */
/*   Input:                                                               */
/*      outData  Output UV object                                         */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitUVDesc *desc;
  ObitTableAN *ANTable=NULL;
  olong ncol;
  gchar *today=NULL;
  olong numOrb,  numPCal;
  ObitIOAccess access;
  olong ver;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar refDate[12];
  gchar *routine = "GetHeader";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outData));
  g_assert(myInput!=NULL);

  /* Define output descriptor if isNew */
  if (isNew) {
    desc = outData->myDesc;

    /* Define header */
    desc->nvis = 0;
    strncpy (desc->origin, "Obit ", UVLEN_VALUE);
    desc->isort[0] = 'T'; 
    desc->isort[1] = 'B';
    desc->nrparm = 0;
    
    /* Creation date today */
    today = ObitToday();
    strncpy (desc->date, today, UVLEN_VALUE-1);
    if (today) g_free(today);
    desc->JDObs    = 0.0;
    desc->epoch    = 2000.0;
    desc->equinox  = 2000.0;
    desc->obsra    = 0.0;
    desc->obsdec   = 0.0;
    desc->altCrpix = 0.0;
    desc->altRef   = 0.0;
    desc->restFreq = 0.0;
    desc->xshift   = 0.0;
    desc->yshift   = 0.0;
    desc->VelDef   = 0;
    desc->VelReference = 0;
    strncpy (desc->bunit, "        ", UVLEN_VALUE);
    ObitInfoListGet(myInput, "refDate", &type, dim, refDate, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    strncpy (desc->obsdat,     refDate, UVLEN_VALUE);
    strncpy (desc->teles,      "Simulate", UVLEN_VALUE);
    strncpy (desc->observer,   "Simulate", UVLEN_VALUE);
    strncpy (desc->instrument, "Simulate", UVLEN_VALUE);
   
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
    strncpy (desc->ptype[ncol], "SOURCE  ", UVLEN_KEYWORD);
    ncol++;
    
    /* FreqID */
    strncpy (desc->ptype[ncol], "FREQSEL ", UVLEN_KEYWORD);
    ncol++;
    
    /* Integration time */
    strncpy (desc->ptype[ncol], "INTTIM  ", UVLEN_KEYWORD);
    ncol++;
    desc->nrparm = ncol;  /* Number of random parameters */
    
    /* Data Matrix */
    ncol = 0;
    /* COMPLEX */
    strncpy (desc->ctype[ncol], "COMPLEX ", UVLEN_KEYWORD);
    desc->inaxes[ncol] = 3;
    desc->cdelt[ncol]  = 1.0;
    desc->crpix[ncol]  = 1.0;
    desc->crval[ncol]  = 1.0;
    desc->crota[ncol]  = 0.0;
    ncol++;

    /* STOKES */
    strncpy (desc->ctype[ncol], "STOKES  ", UVLEN_KEYWORD);
    desc->inaxes[ncol] = nstok;
    desc->cdelt[ncol]  = -1.0;
    desc->crpix[ncol]  =  1.0;
    desc->crval[ncol]  = -1.0;
    desc->crota[ncol]  =  0.0;
    ncol++;

    /* FREQ */
    ObitInfoListGet(myInput, "nFreq",   &type, dim, &nchan,   err);
    ObitInfoListGet(myInput, "refFreq", &type, dim, &refFrequency, err);
    ObitInfoListGet(myInput, "delFreq", &type, dim, &deltaFreq, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    strncpy (desc->ctype[ncol], "FREQ    ", UVLEN_KEYWORD);
    desc->inaxes[ncol] = nchan;
    desc->cdelt[ncol]  = deltaFreq;
    desc->crpix[ncol]  = 1.0;
    desc->crval[ncol]  = refFrequency;
    desc->crota[ncol]  = 0.0;
    ncol++;

    /* IF */
    ObitInfoListGet(myInput, "nIF",   &type, dim, &nIF,   err);
    ObitInfoListGet(myInput, "delIF", &type, dim, &deltaIF, err);
    strncpy (desc->ctype[ncol], "IF      ", UVLEN_KEYWORD);
    desc->inaxes[ncol] = nIF;
    desc->cdelt[ncol]  = 1.0;
    desc->crpix[ncol]  = 1.0;
    desc->crval[ncol]  = 1.0;
    desc->crota[ncol]  = 0.0;
    ncol++;

    /* RA */
    ObitInfoListGet(myInput, "RA", &type, dim, &RA, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    strncpy (desc->ctype[ncol], "RA      ", UVLEN_KEYWORD);
    desc->inaxes[ncol] = 1;
    desc->cdelt[ncol]  = 1.0;
    desc->crpix[ncol]  = 1.0;
    desc->crval[ncol]  = RA;
    desc->crota[ncol]  = 0.0;
    ncol++;

    /* Dec */
    ObitInfoListGet(myInput, "Dec", &type, dim, &Dec, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    strncpy (desc->ctype[ncol], "DEC", UVLEN_KEYWORD);
    desc->inaxes[ncol] = 1;
    desc->cdelt[ncol]  = 1.0;
    desc->crpix[ncol]  = 1.0;
    desc->crval[ncol]  = Dec;
    ncol++;

    desc->naxis = ncol;  /* Number of dimensions */

    /* index descriptor */
    ObitUVDescIndex (desc);

    /* Add Antenna and Frequency info */
    ObitUVOpen (outData, OBIT_IO_WriteOnly, err) ;
    GetFrequencyInfo (myInput, outData, err);
    GetAntennaInfo (myInput, outData, err);
    ObitUVClose (outData, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);

  } /* End define descriptor */

  /* Get source info, copy to output SU table, save lookup table in SourceID */
  ObitUVOpen (outData, OBIT_IO_ReadWrite, err);
  GetSourceInfo (myInput, outData, isNew, err);
  ObitUVClose (outData, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
 
  /* Instantiate output Data */
  ObitUVFullInstantiate (outData, FALSE, err);
  if (err->error)Obit_traceback_msg (err, routine, outData->name);

  /* Get AntennaList in Global */
  /* Create output Antenna table object */
  ver      = 1;
  access   = OBIT_IO_ReadOnly;
  numOrb   = 0;
  numPCal  = 0;
  ANTable = newObitTableANValue ("AN table", (ObitData*)outData, 
				 &ver, access, numOrb, numPCal, err);
  if (ANTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
  AntList = ObitTableANGetList (ANTable, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Cleanup */
  ANTable = ObitTableANUnref(ANTable);

  return;
} /* end GetHeader */

void GetData (ObitUV *outData,  ObitInfoList *myInput,   ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Read data from myInput, write outData                                 */
/*      outData  Output UV Data object, open on input                     */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  olong ant1, ant2,  i, count=0;
  ofloat *Buffer=NULL;
  odouble DecR, RAR, ArrLong, ArrLat, AntLst, HrAng, cosdec, sindec, darg;
  ObitUVDesc *desc;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat time, timeRange[2], delTime, uvw[3], bl[3], minEl, el, chad;
  gchar *routine = "GetData";

  /* error checks */
  if (err->error) return;
  g_assert(myInput!=NULL);
  g_assert(ObitUVIsA(outData));

  /* Get Scan information */
  ObitInfoListGet(myInput, "timeRange", &type, dim, timeRange, err);
  ObitInfoListGet(myInput, "delTime",   &type, dim, &delTime,  err);
  ObitInfoListGet(myInput, "minEl",     &type, dim, &minEl,  err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  delTime /= 86400.0; /* Time increment to days */
  minEl   *= DG2RAD;  /* Min. el to radians */

  /* How many of what? */
  desc = outData->myDesc;
  nchan = desc->inaxes[desc->jlocf];
  nstok = desc->inaxes[desc->jlocs];
  nIF   = desc->inaxes[desc->jlocif];
  nAnt  = AntList->number;
  refFrequency = desc->crval[desc->jlocf];

  /* Prepare output */
  Buffer = outData->buffer;
  desc->firstVis = desc->nvis+1;  /* Append to end of data */
  
  /* Antenna coordinates to wavelengths at reference frequency */
  for (i=0; i<nAnt; i++) {
    AntList->ANlist[i]->AntXYZ[0] *= refFrequency/VELIGHT;
    AntList->ANlist[i]->AntXYZ[1] *= refFrequency/VELIGHT;
    AntList->ANlist[i]->AntXYZ[2] *= refFrequency/VELIGHT;
  }

  /* Position in radians */
  RAR      = RA*DG2RAD;
  DecR     = Dec*DG2RAD;

  /* Array geometry - assume first for all */
  ArrLong = AntList->ANlist[0]->AntLong;
  ArrLat  = AntList->ANlist[0]->AntLat;

  /* Loop over scan */
  time = timeRange[0];
  while (time<=timeRange[1]) {

    /* Min el (same for array) */
    cosdec = cos (DecR);
    sindec = sin (DecR);
    chad = cos (HrAng);
    darg = sin (ArrLat) * sindec + cos (ArrLat) * cosdec * chad;
    el = (1.570796327 - acos (MIN (darg, 1.000)));
    if (el<minEl) continue;

    /* LST and hour angle (radians) */
    AntLst = AntList->GSTIAT0 + ArrLong + time*AntList->RotRate;
    HrAng  = AntLst - RAR;

    /* Loop over antennas  */
    for (ant1=1; ant1<nAnt; ant1++) {
      for (ant2=ant1+1; ant2<=nAnt; ant2++) {
	bl[0] = AntList->ANlist[ant1-1]->AntXYZ[0] - AntList->ANlist[ant2-1]->AntXYZ[0];
	bl[1] = AntList->ANlist[ant1-1]->AntXYZ[1] - AntList->ANlist[ant2-1]->AntXYZ[1];
	bl[2] = AntList->ANlist[ant1-1]->AntXYZ[2] - AntList->ANlist[ant2-1]->AntXYZ[2];
	/* Compute uvw - short baseline approximation */
	ObitUVUtilUVW (bl, DecR, (ofloat)HrAng, uvw);
	/* Set vis in buffer */
	Buffer[desc->ilocu] = uvw[0];
	Buffer[desc->ilocv] = uvw[1];
	Buffer[desc->ilocw] = uvw[2];
	Buffer[desc->ilocb] = (ofloat)(ant1*256+ant2);
	Buffer[desc->iloct] = time;
	Buffer[desc->ilocsu] = SourceID;
	Buffer[desc->ilocfq] = 1;
	Buffer[desc->ilocit] = delTime;
	for (i=0; i<desc->ncorr; i++) { /* Visibilities */
	  Buffer[desc->nrparm+i*3]   = 0.0;
	  Buffer[desc->nrparm+i*3+1] = 0.0;
	  Buffer[desc->nrparm+i*3+2] = 1.0;
	}

	/* set number of records */
	desc->numVisBuff = 1;	
	count += desc->numVisBuff;
	/* Write output one vis at a time */
	if ((ObitUVWrite (outData, NULL, err) != OBIT_IO_OK) || (err->error))
	  Obit_log_error(err, OBIT_Error, "ERROR writing output UV data"); 
	if (err->error) Obit_traceback_msg (err, routine, outData->name);
	
      } /* end loop over second antenna */
    } /* end loop over first antenna */

    /* Next time */
    time += delTime;
  } /* end loop over time */

  /* Tell how many added */
  Obit_log_error(err, OBIT_InfoErr, "Added %d visibilities", count);
  ObitErrLog(err);

} /* end GetData  */

void GetAntennaInfo (ObitInfoList *myInput, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get info from user set parameters                                     */
/*  Writes AN table output                                                */
/*  to output .                                                           */
/*   Input:                                                               */
/*      myInput  Inputs object with user parameters                       */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableAN    *outTable=NULL;
  ObitTableANRow *outRow  =NULL;
  olong          iRow, oRow, ver;
  oint           numPCal, numOrb;
  odouble        JD, GASTM, Rate;
  ObitIOAccess   access;
  ObitInfoType   type;
  gint32         dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar          charTemp[32];
  gchar *routine = "GetAntennaInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outData));

  /* Get antenna info from myInput */
  ObitInfoListGet(myInput, "nAnts",    &type, dim, &nAnt,    err);
  ObitInfoListGet(myInput, "arrayXYZ", &type, dim, arrayXYZ, err);
  ObitInfoListGetP(myInput, "antXYZ",  &type, dim, (gpointer)&antXYZ);

  /* Default array location is VLA */
  if ((fabs(arrayXYZ[0])<0.1) && (fabs(arrayXYZ[1])<0.1) && (fabs(arrayXYZ[2])<0.1)) {
    arrayXYZ[0] = -1.601185365e+06;
    arrayXYZ[1] = -5.041977547e+06;
    arrayXYZ[2] =  3.554875870e+06;
    dim[0] = 3;
    ObitInfoListAlwaysPut(myInput, "arrayXYZ", OBIT_double, dim, arrayXYZ);
  }

  /* Create output Antenna table object */
  ver      = 1;
  access   = OBIT_IO_ReadWrite;
  numOrb   = 0;
  numPCal  = 0;
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
  outTable->ArrayX  = arrayXYZ[0];
  outTable->ArrayY  = arrayXYZ[1];
  outTable->ArrayZ  = arrayXYZ[2];

  /* Compute earth rotation rate and GMST at UT midnight  */
  JD = outData->myDesc->JDObs;
  ObitPrecessGST0 (JD, &GASTM, &Rate);
  outTable->GSTiat0 = GASTM * 15.0;
  if (outTable->GSTiat0<0.0) outTable->GSTiat0 += 360.0;
  GSTiat0 = outTable->GSTiat0;  /* save as global */
  outTable->DegDay  = Rate * 360.0;
  DegDay = outTable->DegDay;   /* save as global */

  outTable->Freq    = refFrequency;
  strncpy (outTable->ArrName, "Simulate", MAXKEYCHARTABLEAN);
  strncpy (outTable->RefDate, outData->myDesc->obsdat, MAXKEYCHARTABLEAN);
  outTable->PolarX  = 0.0;  /* Nothing wobbles on our planet */
  outTable->PolarY  = 0.0;
  outTable->dataUtc = 0.0;  /* In a perfect world */
  strncpy (outTable->TimeSys, "IAT     ", MAXKEYCHARTABLEAN);
  outTable->FreqID = -1;
  outTable->iatUtc = 0.0;   /* In a perfect world */
  strncpy (outTable->polType, "    ", MAXKEYCHARTABLEAN);
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
  outRow->polTypeA[0] = 'R';
  outRow->polTypeB[0] = 'L';
  outRow->status      = 0;
  
  /* loop through input Antenna table */
  for (iRow = 1; iRow<=nAnt; iRow++) {
    /* Set output Row */
    outRow->noSta     = iRow;
    outRow->StaXYZ[0] = antXYZ[(iRow-1)*3];
    outRow->StaXYZ[1] = antXYZ[(iRow-1)*3+1];
    outRow->StaXYZ[2] = antXYZ[(iRow-1)*3+2];
    sprintf (charTemp, "Ant%2.2d          ", iRow);
    strncpy (outRow->AntName, charTemp, outTable->myDesc->repeat[outTable->AntNameCol]);
    outRow->status      = 0;
    
    oRow = outRow->noSta;
    if ((ObitTableANWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating ANTENNA Table");
      return;
    }
  } /* end loop over antennas */
  
  /* Close  table */
  if ((ObitTableANClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output Antenna Table file");
    return;
  }

  /* Cleanup */
  outRow   = ObitTableANRowUnref(outRow);
  outTable = ObitTableANUnref(outTable);

} /* end  GetAntennaInfo */

void GetFrequencyInfo (ObitInfoList *myInput, ObitUV *outData, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get Source info from myInput                                          */
/*  Write FQ table                                                        */
/*   Input:                                                               */
/*      myInput  Inputs object with user parameters                       */
/*      outData  Output UV object                                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableFQ*            outTable=NULL;
  ObitTableFQRow*         outRow=NULL;
  olong i, oRow, ver;
  oint numIF;
  ObitIOAccess access;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "GetFrequencyInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outData));

  /* Get frequency information */
  ObitInfoListGet(myInput, "nIF",     &type, dim, &nIF,   err);
  ObitInfoListGet(myInput, "delIF",   &type, dim, &deltaIF, err);
  ObitInfoListGet(myInput, "nFreq",   &type, dim, &nchan,   err);
  ObitInfoListGet(myInput, "refFreq", &type, dim, &refFrequency, err);
  ObitInfoListGet(myInput, "delFreq", &type, dim, &deltaFreq, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Create output FQ table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  numIF = nIF;
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
  outRow->fqid    = 1;
  for (i=0; i<numIF; i++) {
    outRow->freqOff[i]  = i*deltaIF;
    outRow->chWidth[i]  = deltaFreq;
    outRow->totBW[i]    = deltaFreq*nchan;
    outRow->sideBand[i] = 0;
  }
  outRow->status    = 0;

  /* Write */
  oRow = outRow->fqid;
  if (oRow<1) oRow = -1;
  if ((ObitTableFQWriteRow (outTable, oRow, outRow, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR updating Source Table");
    return;
  }

  /* Close  table */
  if ((ObitTableFQClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output Source Table file");
    return;
  }

  /* Cleanup */
  outRow   = ObitTableFQRowUnref(outRow);
  outTable = ObitTableFQUnref(outTable);

} /* end  GetFrequencyInfo */

void GetSourceInfo (ObitInfoList *myInput, ObitUV *outData, gboolean isNew, 
		    ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get info from user set parameters                                     */
/*  Copies source info from myInput to output data and generated a lookup */
/*  table, global, SourceID to give translation from input source ID      */
/*  to output .                                                           */
/*   Input:                                                               */
/*      myInput  Inputs object with user parameters                       */
/*      outData  Output UV object                                         */
/*      isNew    True if output file just created                         */
/*   Output:                                                              */
/*       err     Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  ObitTableSU    *outTable=NULL;
  ObitTableSURow *outRow=NULL;
  ObitSource *source=NULL;
  olong i, oRow, ver;
  oint numIF;
  gboolean found;
  ObitIOAccess access;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar Source[20];
  gchar *routine = "GetSourceInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outData));

  /* Get Source information */
  ObitInfoListGet(myInput, "Source",  &type, dim, Source,   err);
  Source[dim[0]] = 0;
  ObitInfoListGet(myInput, "RA", &type, dim, &RA, err);
  ObitInfoListGet(myInput, "Dec", &type, dim, &Dec, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Create output Source table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  numIF  = nIF;
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
  for (i=0; i<nIF; i++) {
    outRow->IFlux[i]     = 0.0;
    outRow->QFlux[i]     = 0.0;
    outRow->UFlux[i]     = 0.0;
    outRow->VFlux[i]     = 0.0;
    outRow->FreqOff[i]   = 0.0;
    outRow->LSRVel[i]    = 0.0;
    outRow->RestFreq [i] = 0.0;
  }
  outRow->status    = 0;


  /* See if source exists in output table */
  found = FALSE;
  for (oRow = 1; oRow<=outTable->myDesc->nrow; oRow++) {
    if ((ObitTableSUReadRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR reading Source Table");
      return;
    }
    if (!strncmp (Source, outRow->Source, 16)) {
      /* Found match */
      found = TRUE;
      break;
    }  
  } /* end loop over table */
  
  /* If found just remember output source ID */
  if (found) {
    /* Save source no. in global lookup table */
    SourceID = outRow->SourID;
  } else { /* Not found - add */
    /* If first time update header */
    if (isNew) {
      outTable->FreqID = 1;
      strncpy (outTable->velDef,  "TOPOCENT", MAXKEYCHARTABLESU);
      strncpy (outTable->velType, "        ", MAXKEYCHARTABLESU);
    }
    /* Set output row for end of table */
    outRow->SourID    = outTable->myDesc->nrow+1;
    outRow->Qual      = 0;
    outRow->RAMean    = RA;
    outRow->DecMean   = Dec;
    outRow->Epoch     = 2000.0;
    outRow->PMRa      = 0.0;
    outRow->PMDec     = 0.0;
    strncpy (outRow->Source, Source, 16);
    strncpy (outRow->CalCode, "    ", 4);
    /* Precess */
    source = newObitSource(NULL);
    source->equinox= outRow->Epoch;
    source->RAMean  = outRow->RAMean;
    source->DecMean = outRow->DecMean;
    ObitPrecessUVJPrecessApp (outData->myDesc, source);
    outRow->RAApp  = source->RAApp;
    outRow->DecApp = source->DecApp;
    RAApp  = source->RAApp;          /* Save to global */
    DecApp = source->DecApp;         /* Save to global */
    source = ObitSourceUnref(source);
    for (i=0; i<nIF; i++) {
      outRow->IFlux[i]    = 0.0;
      outRow->QFlux[i]    = 0.0;
      outRow->UFlux[i]    = 0.0;
      outRow->VFlux[i]    = 0.0;
      outRow->FreqOff[i]  = 0.0;
      outRow->LSRVel[i]   = 0.0;
      outRow->RestFreq[i] = 0.0;
    }
    outRow->status    = 0;
    
    /* Save source no. in global */
    SourceID = outRow->SourID;
    
    oRow = outRow->SourID;
    if ((ObitTableSUWriteRow (outTable, oRow, outRow, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR updating Source Table");
      return;
    }
  } /* End add new entry */  
  /* Close  table */
  if ((ObitTableSUClose (outTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output Source Table file");
    return;
  }

  /* Cleanup */
  outRow   = ObitTableSURowUnref(outRow);
  outTable = ObitTableSUUnref(outTable);
} /* end  GetSourceInfo */

/*----------------------------------------------------------------------- */
/*  Write History for UVSim                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      outData   ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void UVSimHistory (ObitInfoList* myInput, ObitUV* outData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntriesNew[] = {
    "refDate", "timeRange", "delTime", 
    "Source", "RA", "Dec", "minEl",
    "refFreq", "nFreq", "delFreq", "nIF", "delIF",
    "arrayXYZ", "nAnts", "antXYZ",
    NULL};
  gchar        *hiEntriesOld[] = {
    "timeRange", "delTime",  "Source", "RA", "Dec", "minEl",
    NULL};
  gchar *routine = "UVSimHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(outData));

  /* Do history  */
  if (isNew)
    outHistory = newObitDataHistory ((ObitData*)outData, OBIT_IO_WriteOnly, err);
  else
    outHistory = newObitDataHistory ((ObitData*)outData, OBIT_IO_ReadWrite, err);

  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Copy selected values from myInput - only Scan info for existing file */
  if (isNew)
    ObitHistoryCopyInfoList (outHistory, pgmName, hiEntriesNew, myInput, err);
  else
    ObitHistoryCopyInfoList (outHistory, pgmName, hiEntriesOld, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  outHistory = ObitHistoryUnref(outHistory);
 
} /* end UVSimHistory  */
