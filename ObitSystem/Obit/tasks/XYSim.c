/* $Id$  */
/* Simulate UV data for antennas with linear (XY) feeds               */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2025                                               */
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

#include "ObitUVDesc.h"
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
#include "ObitTableCCUtil.h"
#include "ObitPrecess.h"
#include "ObitAntennaList.h"
#include "ObitSourceList.h"
#include "ObitTableSU.h"
#include "ObitTableSUUtil.h"
#include "ObitHistory.h"
#include "ObitSkyModel.h"
#include "ObitSkyModelMF.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif

/* internal prototypes */
/* Get inputs */
ObitInfoList* XYSimin (int argc, char **argv, ObitErr *err);
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
/* Digest inputs */
void digestInputs(ObitInfoList *myInput, ObitErr *err);
/* Get input sky model */
ObitSkyModel* getInputSkyModel (ObitInfoList *myInput, ObitErr *err);
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
void XYSimHistory (ObitInfoList* myInput, ObitUV* outData, ObitErr* err);

/* Program globals */
gchar *pgmName = "XYSim";       /* Program name */
gchar *infile  = "XYSim.inp";   /* File with program inputs */
gchar *outfile = "XYSim.out";   /* File to contain program outputs */
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
gchar   *antName=NULL; /* Antenna names (8 char) */
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
  ObitSkyModel *skyModel = NULL;
  ObitUV       *outData  = NULL, *scrData=NULL;
  ObitErr      *err      = NULL;
  ofloat       noise, scale=1.0, IQUV[4]={0.,0.,0.,0.};
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM];

  err = newObitErr();

  /* Startup - parse command line */
  ierr = 0;
  myInput = XYSimin (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);   if (ierr!=0) goto exit;}

  /* Digest inputs */
  digestInputs(myInput, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Create ObitUV for data */
  outData = setOutputData (myInput, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;}
   
  /* Get source, header info, array geometry, initialize output if necessary */
  GetHeader (outData, myInput, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;}

  /* Add IQUV? */
  ObitInfoListGetTest(myInput, "IQUV", &type, dim, IQUV);
  if (IQUV[0]>0.0) {
    Obit_log_error(err, OBIT_InfoErr, "Adding pt. src at origin IQUV=%f %f %f %f", 
		   IQUV[0], IQUV[1], IQUV[2], IQUV[3]);
    ObitErrLog(err);
    
  }
  if (err->error) {ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;}

  /* Open output data - use scrData */
  scrData = newObitUVScratch (outData, err);

  /* Get frequency information */
  ObitUVOpen (scrData, OBIT_IO_WriteOnly, err);
  GetFrequencyInfo (myInput, scrData, err);
  ObitUVClose (scrData, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;}

  if ((ObitUVOpen (scrData, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", 
		   outData->name);
  if (err->error) {ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;}

  /* Simulate data to scratch file */
  GetData (scrData, myInput, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;}

  /* Close */
  if ((ObitUVClose (scrData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing output file");

  /* Add model? */
  skyModel = getInputSkyModel(myInput, err);
  if (skyModel) {
    Obit_log_error(err, OBIT_InfoErr, "Adding sky model");
    ObitErrLog(err);
    ObitSkyModelSubUV (skyModel, scrData, scrData, err);
  }
  if (err->error) {ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;}
    
  /* Add noise? */
  noise = 0.0;
  ObitInfoListGetTest(myInput, "Noise", &type, dim, &noise);
  if (noise>0.0) {
    Obit_log_error(err, OBIT_InfoErr, "Adding %f noise", noise);
    ObitErrLog(err);
    
    noise /= sqrt(2.0);  /* value per real/imag */
    ObitUVUtilNoise(scrData, scrData, scale, noise, err);
  }
  if (err->error) {ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;}

  /* Append scratch file to output */
  Obit_log_error(err, OBIT_InfoErr, "Appending to output");
  ObitErrLog(err);
  ObitUVUtilAppend (scrData, outData, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;}

  /* History */
  XYSimHistory (myInput, outData, err);
  
  /* show any errors */
  if (err->error) {ierr = 1;   ObitErrLog(err);  if (ierr!=0) goto exit;}
  
  /* Shutdown Obit */
 exit:
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  mySystem = ObitSystemShutdown (mySystem);
  
  /* cleanup */
  myInput  = ObitInfoListUnref(myInput);   /* delete input list */
  myInput  = ObitInfoListUnref(myOutput);  /* delete output list */
  outData  = ObitUnref(outData);
  skyModel = ObitUnref(skyModel);
 
  return ierr;
} /* end of main */

ObitInfoList* XYSimin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="XYSim.in", *arg;
  gboolean init=FALSE;
  oint itemp;
  gchar *strTemp;
  ObitInfoList* list;
  gchar *routine = "XYSimin";

  /* Make default inputs InfoList */
  list = defaultInputs(err);
  myOutput = defaultOutputs(err);

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
  if (err->error) {Obit_traceback_val (err, routine, "GetInput", list);}

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
} /* end XYSimin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: XYSim -input file -output ofile [args]\n");
    fprintf(stderr, "Simulate UV data and write Obit/UV\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def XYSim.in\n");
    fprintf(stderr, "  -output output result file, def XYSim.out\n");
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
  olong      i, n, Aseq, disk, cno;
  gchar     *Type, *strTemp, outFile[129];
  gchar     Aname[13], Aclass[7], *Atype = "UV";
  olong      nvis;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean  exist;
  gchar     tname[129], *fullname=NULL;
  gchar     *outParms[] = {  /* Parameters for output data */
    "Compress", NULL};
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
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */

    /* outName given? */
    for (i=0; i<12; i++) {Aname[i] = ' ';}  Aname[i] = 0;
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&strTemp);
    /* if not use "No Name" */
    if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12))) {
      strncpy (Aname, "No Name     ", 13);
    } else {for (i=0; i<MIN(12,dim[0]); i++) Aname[i] = strTemp[i];}
    /* Save any defaulting on myInput */
    dim[0] = 12;
    ObitInfoListAlwaysPut (myInput, "outName", OBIT_string, dim, Aname);

      
    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    /* Default out class is "XYSim" */
    if (!strncmp(Aclass, "      ", 6)) strncpy (Aclass, "XYSim", 7);

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
    nvis = 1000;
    ObitUVSetAIPS (outUV, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */

    /* outFile given? */
    ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&strTemp);
    /* if not use inName */
    if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12)))
      ObitInfoListGetP (myInput, "in2File", &type, dim, (gpointer)&strTemp);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) {outFile[i] = strTemp[i];} outFile[i] = 0;
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
    nvis = 1000;
    ObitUVSetFITS (outUV, nvis, disk, outFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outUV;
  }
  
  /* Copy control parameters if this is a new file */
  if (isNew) ObitInfoListCopyList (myInput, outUV->info, outParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

  return outUV;
} /* end setOutputUV */

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
  gchar *strTemp;
  gboolean replace;
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
  replace = TRUE;
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
 
} /* end digestInputs */

/*----------------------------------------------------------------------- */
/*  Get input sky model                                                   */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*      Sky Model to be used, NULL if none given                          */
/*----------------------------------------------------------------------- */
ObitSkyModel* getInputSkyModel (ObitInfoList *myInput, ObitErr *err)
{
  ObitSkyModel *skyModel=NULL;
  ObitImageMosaic *mosaic=NULL;
  ObitImage    **image=NULL;
  ObitCCCompType CCType;
  ObitInfoType type;
  gboolean     do3D=TRUE;
  olong        Aseq, disk, cno,i=0, ver, nmaps;
  gchar        *Type, *Type2, *strTemp, inFile[129], inRoot[129];
  gchar        Aname[13], Aclass[7], Aroot[7], *Atype = "MA";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ofloat       smodel[20], modptflx,  modptxof, modptyof, modptypm[8];
  gchar        name[101];
  gchar        *dataParms[] = {  /* Control parameters */
    "CCVer",  "BComp",  "EComp",  "Flux", "PBCor", "antSize", "Factor", 
    "minFlux", "Mode", "ModelType", "REPLACE", "Stokes", 
    "MODPTFLX", "MODPTXOF", "MODPTYOF", "MODPTYPM", "doGPU",
    NULL};
  gchar *routine = "getInputSkyModel";

  /* error checks */
  if (err->error) return skyModel;
  g_assert (ObitInfoListIsA(myInput));

  /* How many fields? */
  nmaps = 1;
  ObitInfoListGetTest(myInput, "nmaps", &type, dim, &nmaps);

  /* Image model of model in parameters? */
  smodel[0] = 0.0;
  ObitInfoListGetTest (myInput, "modelFlux", &type, dim, &smodel[0]);
  if ((smodel!=NULL) && (smodel[0]!=0.0)) {
    /* Model passed - get rest */
    for (i=1; i<11; i++) smodel[i] = 0.0;
    ObitInfoListGetTest (myInput, "modelPos",  &type, dim, &smodel[1]);
    ObitInfoListGetTest (myInput, "modelParm", &type, dim, &smodel[3]);
    modptflx = smodel[0];
    modptxof = smodel[1] / 3600.0;
    modptyof = smodel[2] / 3600.0;
    for (i=0; i<8; i++) modptypm[i] = smodel[3+i];
    dim[0] = dim[1] = 1;
    ObitInfoListAlwaysPut (myInput, "MODPTFLX", OBIT_float, dim, &modptflx);
    ObitInfoListAlwaysPut (myInput, "MODPTXOF", OBIT_float, dim, &modptxof);
    ObitInfoListAlwaysPut (myInput, "MODPTYOF", OBIT_float, dim, &modptyof);
    dim[0] = 6;
    ObitInfoListAlwaysPut (myInput, "MODPTYPM", OBIT_float, dim, modptypm);

    /* Create Sky Model */
    skyModel = newObitSkyModel ("Sky Model");

    Obit_log_error(err, OBIT_InfoErr, "Using input model parameters");
  } else if (nmaps>0) {  /* Image given */
    
    Obit_log_error(err, OBIT_InfoErr, "Using image model");

    /* Allocate Image array */
    image = g_malloc0(nmaps*sizeof(ObitImage));
    
    /* Create image mosaic */
    mosaic = newObitImageMosaic ("Mosaic", nmaps);
    
    /* File type - could be either AIPS or FITS use DataType2 (default DataType) */
    ObitInfoListGetP (myInput, "DataType",  &type, dim, (gpointer)&Type);
    ObitInfoListGetP (myInput, "DataType2", &type, dim, (gpointer)&Type2);
    if (!strncmp (Type2, "    ", 4)) Type2 = Type;
    if (!strncmp (Type2, "AIPS", 4)) { /* AIPS input */
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

      /* input AIPS sequence */
      ObitInfoListGet(myInput, "in2Seq", &type, dim, &Aseq, err);
      
      /* if ASeq==0 want highest existing sequence */
      if (Aseq<=0) {
	/* If only one field use class given */
	if ((nmaps==1) && (i==0)) {
	} else { /* derive class from field number */
	  Aroot[2] = 0;
	  g_snprintf (Aclass, 7, "%s%4.4d",Aroot,i+1);
	} /* end one or many fields */
	Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	/* Save on myInput*/
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(myInput, "inSeq", OBIT_oint, dim, &Aseq);
      }
      
      /* Loop over fields */
      for (i=0; i<nmaps; i++) {
	g_snprintf (name, 100, "Input image %d",i+1);
	/* If only one field use class given */
	if ((nmaps==1) && (i==0)) {
	  g_snprintf (Aclass, 7, "%s",Aroot);
	} else { /* derive class from field number */
	  Aroot[2] = 0;
	  g_snprintf (Aclass, 7, "%s%4.4d",Aroot,i+1);
	} /* end one or many fields */
	
	  /* Find catalog number */
	cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
	if (cno<0) Obit_log_error(err, OBIT_Error, "Failure looking up %s", name);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	
	/* define object */
	image[i] = newObitImage(name);
	ObitImageSetAIPS(image[i], OBIT_IO_byPlane, disk, cno, AIPSuser,  blc, trc, err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	
	/* Ensure image fully instantiated and OK */
	ObitImageFullInstantiate (image[i], TRUE, err);

	/* Attach Image */
	ObitImageMosaicSetImage (mosaic, i, image[i], err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
      } /* end loop over fields */
      
    } else if (!strncmp (Type2, "FITS", 4)) {  /* FITS input */
      /* input FITS file name */
      if (ObitInfoListGetP(myInput, "in2File", &type, dim, (gpointer)&strTemp)) {
	strncpy (inRoot, strTemp, 128);
      } else { 
	strncpy (inRoot, "No_Filename_Given", 128);
      }
      ObitTrimTrail(inRoot);  /* remove trailing blanks */
   
      /* input FITS disk */
      ObitInfoListGet(myInput, "in2Disk", &type, dim, &disk, err);

      if (nmaps==1) {
 	
	/* Set file name */
	g_snprintf (inFile, 128, "%s",inRoot);
 
	/* define object */
	g_snprintf (name, 100, "Input image");
	image[0] = newObitImage(name);
	ObitImageSetFITS(image[0], OBIT_IO_byPlane, disk, inFile, blc, trc, err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	
	/* Ensure image fully instantiated and OK */
	ObitImageFullInstantiate (image[0], TRUE, err);

	/* Attach Image */
	ObitImageMosaicSetImage (mosaic, 0, image[0], err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
     } else { /* Multiple fields */
	
	/* Loop over fields */
	for (i=0; i<nmaps; i++) {
	  /* Set file name */
	  g_snprintf (inFile, 128, "%s%d",inRoot,i);
	  
	  /* define object */
	  g_snprintf (name, 100, "Input image %d",i+1);
	  image[i] = newObitImage(name);
	  ObitImageSetFITS(image[i], OBIT_IO_byPlane, disk, inFile, blc, trc, err);
	  if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	  
	  /* Ensure image fully instantiated and OK */
	  ObitImageFullInstantiate (image[i], TRUE, err);

	  /* Attach Image */
	  ObitImageMosaicSetImage (mosaic, i, image[i], err);
	  if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	} /* end loop over fields */
      }

   } else { /* Unknown type - barf and bail */
      Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		     pgmName, Type2);
      return skyModel;
    }

    /* Create Sky Model for appropriate type */
    ver = 0;
    ObitInfoListGetTest(myInput, "CCVer", &type, dim, &ver);
    if (ver>0) CCType = ObitTableCCUtilGetType ((ObitData*)mosaic->images[0], ver, err);
    else       CCType = OBIT_CC_Unknown;
    if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
    if ((CCType==OBIT_CC_PointModTSpec)|| (CCType==OBIT_CC_GaussModTSpec) ||
	(CCType==OBIT_CC_CGaussModTSpec) || (CCType==OBIT_CC_USphereModTSpec)) {
      skyModel = (ObitSkyModel*)ObitSkyModelMFCreate ("Sky Model", mosaic);
      Obit_log_error(err, OBIT_InfoErr, "Using tabulated spectrum sky model");
    } else
      skyModel = ObitSkyModelCreate ("Sky Model", mosaic);

    /* deallocate images */
    for (i=0; i<nmaps; i++) image[i] = ObitImageUnref(image[i]);
    g_free(image);  /* Deallocate array */

   /* End image or components model */
  } else {  /* Neither model given just return */
    Obit_log_error(err, OBIT_InfoErr, 
		   "Using source fluxes from SU table for point model");
    return skyModel;
  }
 
  /* Get input parameters from myInput, copy to skyModel */
  ObitInfoListCopyList (myInput, skyModel->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, skyModel->name, skyModel);
  
  /* get do3D from first image */
  if ((skyModel->mosaic) && (skyModel->mosaic->images[0]))
    do3D = skyModel->mosaic->images[0]->myDesc->do3D;
  else do3D = FALSE;
      
  /* Save do3D */
  dim[0] = 1; dim[1] = 1;
  ObitInfoListAlwaysPut (skyModel->info, "do3D", OBIT_bool, dim, &do3D);
  
  return skyModel;
} /* end getInputSkyModel */
  
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
  olong numIF=1, numOrb,  numPCal;
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

  /* Get position of this target */
  ObitInfoListGet(myInput, "RA", &type, dim, &RA, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  ObitInfoListGet(myInput, "Dec", &type, dim, &Dec, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
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
    desc->crval[ncol]  = -5.0;  /* XX,YY,XY,YX */
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
    numIF = nIF;
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
    strncpy (desc->ctype[ncol], "RA      ", UVLEN_KEYWORD);
    desc->inaxes[ncol] = 1;
    desc->cdelt[ncol]  = 1.0;
    desc->crpix[ncol]  = 1.0;
    desc->crval[ncol]  = RA;
    desc->crota[ncol]  = 0.0;
    ncol++;

    /* Dec */
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
  access   = OBIT_IO_ReadWrite;
  numOrb   = 0;
  numPCal  = 0;
  ANTable = newObitTableANValue ("AN table", (ObitData*)outData, 
				 &ver, access, numIF, numOrb, numPCal, err);
  if (ANTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
  AntList = ObitTableANGetList (ANTable, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Cleanup */
  ANTable = ObitTableANUnref(ANTable);

  return;
} /* end GetHeader */

void GetData (ObitUV *outData, ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Read description from myInput, write scratch uv data                  */
/*      outData  Output UV Data object, open on input                     */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  olong ant1, ant2,  suba=1, i, count=0, suId=1;
  ofloat *Buffer=NULL, IQUV[4], ipol=0.0, qpol=0.0, upol=0.0, vpol=0.0;
  ofloat PA, XYVis[12];  /* Polarized model visibility */
  gboolean doPol=FALSE;
  odouble DecR, RAR, ArrLong, ArrLat, AntLst, HrAng=0.0, cosdec, sindec, darg;
  ObitUVDesc *desc;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat time, timeRange[2], UVRange[2]={0.0,1.0e10}, delTime;
  ofloat uvw[3], bl[3], blen2, minEl, el, chad;
  ObitSourceList *SList=NULL;
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
  /* Square of UV Range */
  ObitInfoListGetTest(myInput, "UVRange", &type, dim, UVRange);
  if (UVRange[1]<=0.0) UVRange[1] = 1.0e10;  /* Sky's the limit */
  UVRange[0]  *= 1000.0;     UVRange[1]  *= 1000.0;  /* to wavelengths */
  UVRange[0]  *= UVRange[0]; UVRange[1]  *= UVRange[1];
  /* Polarized model */
  ObitInfoListGetTest(myInput, "IQUV", &type, dim, IQUV);
  if (IQUV[0]>0.0) {
    doPol = TRUE;
    ipol = IQUV[0]; qpol = IQUV[1]; upol = IQUV[2]; vpol = IQUV[3]; 
  }
 
 
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

  /* Set position in header */
  outData->myDesc->crval[outData->myDesc->jlocr] = RA;
  outData->myDesc->crval[outData->myDesc->jlocd] = Dec;

  /* Position in radians */
  RAR      = RA*DG2RAD;
  DecR     = Dec*DG2RAD;

  /* Array geometry - assume first for all */
  ArrLong = AntList->ANlist[0]->AntLong;
  ArrLat  = AntList->ANlist[0]->AntLat;

  /* Get Source List */
  SList = ObitUVGetSourceList (outData, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  if (outData->myDesc->ilocsu<0) {
    /* May have selected one source */
    if (outData->mySel->selectSources && (outData->mySel->numberSourcesList==1)) 
         suId = outData->mySel->sources[0];
    else suId = 1;
  }
 
    /* Loop over scan */
  time = timeRange[0];
  while (time<=timeRange[1]) {

    /* LST and hour angle (radians) */
    AntLst = AntList->GSTIAT0 + ArrLong + time*AntList->RotRate;
    HrAng  = AntLst - RAR;

    /* Min el (same for array) */
    cosdec = cos (DecR);
    sindec = sin (DecR);
    chad = cos (HrAng);
    darg = sin (ArrLat) * sindec + cos (ArrLat) * cosdec * chad;
    el = (1.570796327 - acos (MIN (darg, 1.000)));
    if (el<minEl) goto skip;

    /* Polarized visibility */
    if (doPol) {
      /* Need parallactic angle */
      PA = ObitAntennaListParAng (AntList, 1, time,  SList->SUlist[suId-1]);
      /* XX */
      XYVis[0] = ipol + qpol*cos(2*PA) + upol*sin(2*PA); XYVis[1] = 0.0; XYVis[2] = 1.0;
      /* YY */
      XYVis[3] = ipol - qpol*cos(2*PA) - upol*sin(2*PA); XYVis[4] = 0.0; XYVis[5] = 1.0;
      /* XY */
      XYVis[6] = - qpol*sin(2*PA) + upol*cos(2*PA); XYVis[7]  = +vpol; XYVis[8] = 1.0;
      /* YX */
      XYVis[9] = - qpol*sin(2*PA) + upol*cos(2*PA); XYVis[10] = -vpol; XYVis[11] = 1.0;
    }

    /* Loop over antennas  */
    for (ant1=1; ant1<nAnt; ant1++) {
      for (ant2=ant1+1; ant2<=nAnt; ant2++) {
	bl[0] = AntList->ANlist[ant1-1]->AntXYZ[0] - AntList->ANlist[ant2-1]->AntXYZ[0];
	bl[1] = AntList->ANlist[ant1-1]->AntXYZ[1] - AntList->ANlist[ant2-1]->AntXYZ[1];
	bl[2] = AntList->ANlist[ant1-1]->AntXYZ[2] - AntList->ANlist[ant2-1]->AntXYZ[2];
	/* In UV Range? */
	blen2 = bl[0]*bl[0] + bl[1]*bl[1];
	if ((blen2<UVRange[0]) || (blen2>UVRange[1])) continue;
	/* Compute uvw - short baseline approximation */
	ObitUVUtilUVW (bl, DecR, (ofloat)HrAng, uvw);
	/* Set vis in buffer */
	Buffer[desc->ilocu] = uvw[0];
	Buffer[desc->ilocv] = uvw[1];
	Buffer[desc->ilocw] = uvw[2];
	ObitUVDescSetAnts(desc, Buffer, ant1, ant2, suba);
	Buffer[desc->iloct] = time;
	Buffer[desc->ilocsu] = SourceID;
	if (doPol) {
	  for (i=0; i<desc->ncorr; i+=4) { /* model visibilities */
	    Buffer[desc->nrparm+i*3+0] = XYVis[0]; /* XX */
	    Buffer[desc->nrparm+i*3+1] = XYVis[1];
	    Buffer[desc->nrparm+i*3+2] = XYVis[2];
	    Buffer[desc->nrparm+i*3+3] = XYVis[3]; /* YY */
	    Buffer[desc->nrparm+i*3+4] = XYVis[4];
	    Buffer[desc->nrparm+i*3+5] = XYVis[5];
	    Buffer[desc->nrparm+i*3+6] = XYVis[6]; /* XY */
	    Buffer[desc->nrparm+i*3+7] = XYVis[7];
	    Buffer[desc->nrparm+i*3+8] = XYVis[8];
	    Buffer[desc->nrparm+i*3+9] = XYVis[9]; /* YX */
	    Buffer[desc->nrparm+i*3+10] = XYVis[10];
	    Buffer[desc->nrparm+i*3+11] = XYVis[11];
	  }
	} else { /* zeroes */
	  for (i=0; i<desc->ncorr; i++) { /* Visibilities */
	    Buffer[desc->nrparm+i*3]   = 0.0;
	    Buffer[desc->nrparm+i*3+1] = 0.0;
	    Buffer[desc->nrparm+i*3+2] = 1.0;
	  }
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
  skip:
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
  oint           numIF, numPCal, numOrb;
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
  ObitInfoListGetP(myInput, "antName", &type, dim, (gpointer)&antName);

  /* Default array location is Earth Center -  0's OK*/

  /* Create output Antenna table object */
  ver      = 1;
  access   = OBIT_IO_ReadWrite;
  numOrb   = 0;
  numPCal  = 0;
  if (outData->myDesc->jlocif>=0)
    numIF    = outData->myDesc->inaxes[outData->myDesc->jlocif];
  else
    numIF = 1;
  outTable = newObitTableANValue ("Output table", (ObitData*)outData, 
				  &ver, access, numIF, numOrb, numPCal, err);
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
  outTable->PolarX  = 90.0; 
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
  outRow->polTypeA[0] = 'X';
  outRow->polTypeB[0] = 'Y';
  outRow->status      = 0;
  
  /* loop through input Antenna table */
  for (iRow = 1; iRow<=nAnt; iRow++) {
    /* Set output Row */
    outRow->noSta     = iRow;
    outRow->StaXYZ[0] = antXYZ[(iRow-1)*3];
    outRow->StaXYZ[1] = antXYZ[(iRow-1)*3+1];
    outRow->StaXYZ[2] = antXYZ[(iRow-1)*3+2];
    if (!strncmp("        ",&antName[(iRow-1)*8],8)) {
	sprintf (charTemp, "Ant%2.2d          ", iRow);
	strncpy (outRow->AntName, charTemp, outTable->myDesc->repeat[outTable->AntNameCol]);
      } else {
	strncpy (outRow->AntName, &antName[(iRow-1)*8], 8);
      }
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
    strncpy (outRow->CalCode, "   ", 4);
    /* Precess */
    source = newObitSource(NULL);
    source->equinox = outRow->Epoch;
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
/*  Write History for XYSim                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      outData   ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void XYSimHistory (ObitInfoList* myInput, ObitUV* outData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntriesNew[] = {
    "refDate", "timeRange", "delTime", "UVRange", 
    "Source", "RA", "Dec", "minEl",
    "refFreq", "nFreq", "delFreq", "nIF", "delIF",
    "arrayXYZ", "nAnts", "antXYZ", "antName", "Noise", "Compress", "IQUV",
    "DataType", "in2File", "in2Name", "in2Class", "in2Disk", "in2Seq",
    "nmaps", "CCVer", "BComp", "Ecomp", "Flux", "Factor", "Cmethod", "Cmodel",
    "modelFlux", "modelPos", "modelParm", "mrgCC', PBCor", "antSize",
    "outDType", "outFile", "outName", "outClass", "outSeq", "outDisk", 
    NULL};
  gchar        *hiEntriesOld[] = {
    "timeRange", "delTime",  "Source", "RA", "Dec", "minEl", "Noise",
    NULL};
  gchar *routine = "XYSimHistory";

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
 
} /* end XYSimHistory  */
