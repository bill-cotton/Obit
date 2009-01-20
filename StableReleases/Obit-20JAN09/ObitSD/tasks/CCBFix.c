/* $Id$  */
/* Fix screwup in  CalTech Continuum Backend data                     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008                                               */
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

#include "ObitOTF.h"
#include "ObitIOOTFFITS.h"
#include "ObitFITS.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitTableOTFArrayGeom.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* CCBFixin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(gchar *scan_name, ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Get input data */
ObitOTF* getInputData (ObitInfoList *myInput, ObitErr *err);
/* Create output data */
ObitOTF* setOutputData (ObitInfoList *myInput, ObitOTF* inData, gboolean *exist,
			ObitErr *err);
/* Copy Array Geometry table */
void CopyData (ObitInfoList *myInput, ObitOTF* inData, ObitOTF* outData, ObitErr *err);


/* Program globals */
gchar *pgmName = "CCBFix";       /* Program name */
gchar *input_file  = "CCBFix.inp";   /* File with program inputs */
gchar *output_file = "CCBFix.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
/* gchar **FITSdirs=NULL; List of FITS data directories */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*    Fix screwup in CalTech Continuum Backend data                       */
/*----------------------------------------------------------------------- */
{
  olong  ierr=0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *inData=NULL, *outData= NULL;
  gboolean exist;
  ObitErr      *err= NULL;

  err = newObitErr();

  /* Startup - parse command line */
  ierr = 0;
  myInput = CCBFixin (argc, argv, err);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;


  /* Get input data */
  inData =  getInputData (myInput, err);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Create output file - copy some tables */
  outData = setOutputData(myInput, inData, &exist, err);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Copy data */
  CopyData(myInput, inData, outData, err);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;

  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;
   
   /* Shutdown Obit */
 exit:
   ObitReturnDumpRetCode (ierr, output_file, myOutput, err);  /* Final output */
   mySystem = ObitSystemShutdown (mySystem);
   
   /* cleanup */
   myInput = ObitInfoListUnref(myInput);  /* delete input list */
   outData = ObitUnref(outData);
   inData = ObitUnref(inData);
 
   return ierr;
} /* end of main */

ObitInfoList* CCBFixin (int argc, char **argv, ObitErr *err)
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
  olong i, j, k, ax, itemp;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *arg;
  gchar *scan_name ="unspecified";
  gboolean init=FALSE;
  gchar *strTemp;
  ObitInfoList* list;
  gchar *routine = "CCBFixin";

  /* Make default inputs & output InfoList */
  list = defaultInputs(scan_name, err);
  myOutput = defaultOutputs(err);
  if (err->error) Obit_traceback_val (err, routine, routine, list);

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
      output_file = argv[++ax];

    } else if (strcmp(arg, "-Scan") == 0){ /* scan name */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "Scan", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outOTF") == 0){ /* Output OTF file */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outOTF", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outDisk") == 0) { /* Output FITS disk */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outDisk", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-pgmNumber") == 0) { /*Program number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "pgmNumber", OBIT_oint, dim, &itemp, err);

    } else if (strcmp(arg, "-AIPSuser") == 0) { /* AIPS User */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "AIPSuser", OBIT_oint, dim, &itemp, err);
      
    } else { /* unknown argument */
      fprintf(stderr,"Unknown parameter %s\n",arg);
      Usage();
    }
  }
  
   /* Read defaults if no file specified */
  if (!init) ObitParserParse (input_file, list, err);
  if (err->error) Obit_traceback_val (err, routine, routine, list);

  /* Extract basic information to program globals */
  pgmNumber = 1;
  ObitInfoListGetTest(list, "pgmNumber", &type, dim, &pgmNumber);
  AIPSuser=1;
  ObitInfoListGetTest(list, "AIPSuser",  &type, dim, &AIPSuser);
  nAIPS = 0;
  ObitInfoListGetTest(list, "nAIPS",     &type, dim, &nAIPS);
  nFITS = 0;
  ObitInfoListGetTest(list, "nFITS",     &type, dim, &nFITS);
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
  ObitReturnDumpRetCode (-999, output_file, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, routine, list);

 return list;
} /* end CCBFixin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: CCBFix -input file [-scan date/time]\n");
    fprintf(stderr, "Convert an GBT CCB file format to Obit/OTF\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def CCBFix.in\n");
    fprintf(stderr, "  -output output parameter file, def CCBFix.out\n");
    fprintf(stderr, "  -outOTF Output OTF file\n");
    fprintf(stderr, "  -outDisk Output FITS disk\n");
    
    /*/exit(1);  bail out */
  }/* end Usage */

/*----------------------------------------------------------------------- */
/*  Create default input ObitInfoList                                     */
/*   Input:                                                               */
/*       scan_name Date/time base name for GBT scan FITS files            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(gchar *scan_name, ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  ofloat ftemp;
  gboolean doBS;
  ObitInfoList *out = newObitInfoList();

  /* add parser items */
  /* base of scan file names */
  dim[0] = strlen (scan_name);
  ObitInfoListPut (out, "Scan", OBIT_string, dim, scan_name, err);

  /* output FITS file name */
  strTemp = "DataOTF.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "outOTF", OBIT_string, dim, strTemp, err);

  /* root of data directory */
  strTemp = "dataRoot";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "DataRoot", OBIT_string, dim, strTemp, err);

  /* Correction in sec to GBT time labels */
  dim[0] = 1;
  ftemp = 0.0;
  ObitInfoListPut (out, "offTime", OBIT_float, dim, &ftemp, err);

  /* Beamswitch? */
  dim[0] = 1;
  doBS = TRUE;
  ObitInfoListPut (out, "doBS", OBIT_bool, dim, &doBS, err);

  /* Data normalization 1 for all */
  dim[0] = 1;
  ftemp = 1.0;
  ObitInfoListPut (out, "dataNorm", OBIT_bool, dim, &ftemp, err);

  return out;
} /* end defaultInputs */

/*----------------------------------------------------------------------- */
/*  Create default output ObitInfoList                                    */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*  Values:  Nothing returned                                             */
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
/*  Get input data                                                        */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitOTF with input data                                          */
/*----------------------------------------------------------------------- */
ObitOTF* getInputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitOTF       *inData = NULL;
  ObitInfoType type;
  olong         disk, nrec=1000;
  gchar        *strTemp, inFile[129];
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "getInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input OTF data Object */
  inData = newObitOTF("input OTF data");
  
  /* input FITS file name */
  if (ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp)) {
    strncpy (inFile, strTemp, 128);
  } else { 
    strncpy (inFile, "No_Filename_Given", 128);
  }
  ObitTrimTrail(inFile);  /* Trim trailing blanks */
 
  /* input FITS disk */
  ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);

  /* define object */
  ObitOTFSetFITS (inData, nrec, disk, inFile,  err); 
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Ensure inData fully instantiated and OK */
  ObitOTFFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Set output data                                                       */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input data                                              */
/*   Output:                                                              */
/*      exist     If true output file already exists                      */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitOTF with output data                                         */
/*----------------------------------------------------------------------- */
ObitOTF* setOutputData (ObitInfoList *myInput, ObitOTF* inData, 
			gboolean *exist, ObitErr *err)
{
  ObitOTF       *outData = NULL;
  ObitInfoType type;
  olong         disk, nrec=1000;
  gchar        *strTemp, outFile[129];
  ofloat       *azOffset=NULL, *elOffset=NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *newFull;
  gchar *exclude[]={"OTFScanData", "OTFSoln", "OTFCal", "OTFArrayGeom", 
		    NULL};
  gchar *routine = "setOutputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input OTF data Object */
  outData = newObitOTF("input OTF data");
  
  /* input FITS file name */
  if (ObitInfoListGetP(myInput, "outFile", &type, dim, (gpointer)&strTemp)) {
    strncpy (outFile, strTemp, 128);
  } else { 
    strncpy (outFile, "No_Filename_Given", 128);
  }
  ObitTrimTrail(outFile);  /* Trim trailing blanks */
 
  /* output FITS disk */
  ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);

  /* Full path */
  newFull = ObitFITSFilename (disk, outFile, err);  

  /* define object */
  ObitOTFSetFITS (outData, nrec, disk, outFile,  err); 
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);

  /* Does it exist */
  *exist = ObitFileExist (newFull, err);
  g_free(newFull);

  if (*exist) {  /* If old exists just verify */
    ObitOTFFullInstantiate (outData, TRUE, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outData);
    return outData;
  }
  
  /* Copy descriptor */
  outData->myDesc = ObitOTFDescCopy (inData->myDesc, outData->myDesc, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);
  outData->myDesc->inaxes[1] = 1;   /* Drop one poln */
  outData->myDesc->colRepeat[outData->myDesc->ilocdata] /= 2;   /* Drop one poln */
  /* Index the descriptor */
  ObitOTFDescIndex (outData->myDesc);

  /* Copy Array Geometry */
  outData->geom = ObitOTFArrayGeomCopy (inData->geom, outData->geom, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);
  
  /* Resize Array Geometry */
  outData->geom->numberDetect = 8;
  azOffset = outData->geom->azOffset;
  elOffset = outData->geom->elOffset;
  outData->geom->azOffset = g_malloc0(outData->geom->numberDetect*sizeof(ofloat));
  outData->geom->elOffset = g_malloc0(outData->geom->numberDetect*sizeof(ofloat));
  outData->geom->azOffset[0] = azOffset[0];
  outData->geom->azOffset[1] = azOffset[3];
  outData->geom->azOffset[2] = azOffset[4];
  outData->geom->azOffset[3] = azOffset[7];
  outData->geom->azOffset[4] = azOffset[8];
  outData->geom->azOffset[5] = azOffset[11];
  outData->geom->azOffset[6] = azOffset[12];
  outData->geom->azOffset[7] = azOffset[15];
  outData->geom->elOffset[0] = elOffset[0];
  outData->geom->elOffset[1] = elOffset[3];
  outData->geom->elOffset[2] = elOffset[4];
  outData->geom->elOffset[3] = elOffset[7];
  outData->geom->elOffset[4] = elOffset[8];
  outData->geom->elOffset[5] = elOffset[11];
  outData->geom->elOffset[6] = elOffset[12];
  outData->geom->elOffset[7] = elOffset[15];
  if (azOffset) g_free(azOffset);
  if (elOffset) g_free(elOffset);

  /* Ensure outData fully instantiated and OK */
  ObitOTFOpen (outData, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);

 /* Copy tables   */
  ObitOTFOpen(outData, OBIT_IO_ReadWrite, err);
  ObitOTFCopyTables (inData, outData, exclude, NULL, err);
  ObitOTFClose(outData, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outData);

  return outData;
} /* end setOutputData */

/*----------------------------------------------------------------------- */
/* Copy Array Geometry table                                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input data                                              */
/*      outData   Output data                                             */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void CopyData (ObitInfoList *myInput, ObitOTF* inData, ObitOTF* outData, 
		ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  olong indx, r, i, j, ncopy=8, numData;
  ofloat *inBuff, *outBuff;
  olong copy[8] = {0,3,4,7,8,11,12,15};
  gchar *routine = "CopyData";
  
  iretCode = ObitOTFOpen (inData, OBIT_IO_ReadOnly, err);
  if ((iretCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, inData->name);
  oretCode = ObitOTFOpen (outData, OBIT_IO_ReadWrite, err);
  if ((oretCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, outData->name);
  outData->myDesc->firstRec = outData->myDesc->nrecord+1;  /* At end */

  /* we're in business, copy */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    iretCode = ObitOTFRead (inData, inData->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;

    /* How many */
    outData->myDesc->numRecBuff = inData->myDesc->numRecBuff;
    numData = inData->myDesc->numRecBuff;

    /* Copy records from buffer */
    inBuff  = inData->buffer;
    outBuff = outData->buffer;
    for (r=0; r<numData; r++) {
      for (j=0; j<inData->myDesc->numDesc; j++) outBuff[j] = inBuff[j];
      j = 0;
      for (i=0; i<ncopy; i++) {
        indx = copy[i]*inData->myDesc->incdatawt;
	outBuff[outData->myDesc->numDesc+j*outData->myDesc->incdatawt]  = 
	  inBuff[inData->myDesc->numDesc+indx];
	if (outData->myDesc->incdatawt>1) {
	  outBuff[outData->myDesc->numDesc+j*outData->myDesc->incdatawt+1]  = 
	    inBuff[inData->myDesc->numDesc+indx+1];
	}
	j++;
	indx++;
      }
      inBuff  += inData->myDesc->lrec;
      outBuff += outData->myDesc->lrec;
    }
    
    oretCode = ObitOTFWrite (outData, outData->buffer, err);
  }  /* End loop copying file */
  
  /* close files */
  oretCode = ObitOTFClose (outData, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, outData->name);
  
  iretCode = ObitOTFClose (inData, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inData->name);
} /* end CopyArray */
