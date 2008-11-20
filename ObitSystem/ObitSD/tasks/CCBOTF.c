/* $Id$  */
/* Read CalTech Continuum Backend data                                */
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

#include "ObitOTF.h"
#include "ObitIOOTFFITS.h"
#include "ObitFITS.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitTableOTFTarget.h"
#include "ObitTableOTFIndex.h"
#include "ObitTableGBTANTPOSGR.h"
#include "ObitTableGBTANTPOSPF.h"
#include "ObitTableGBTCCBDATA.h"
#include "ObitGBTIFInfo.h"
#include "ObitGBTCCBStateInfo.h"
#include "ObitGBTCCBPortInfo.h"
#include "ObitGBTBeamOffInfo.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* CCBOTFin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(gchar *scan_name, ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Get pointing times from Antenna file, Secondary focus */
void GetAntennaGR (gchar *infile, ObitInfoList *myInput, ObitErr *err);
/* Get pointing times from Antenna file, Prime Focus */
void GetAntennaPF (gchar *infile, ObitInfoList *myInput, ObitErr *err);
/* Get pointing for a given time */
void GetPoint (odouble time, ofloat *ra, ofloat *dec);
/* New or old version of the CCB? */
gboolean isNewOld (gchar *infile, ObitErr *err);
/* Get new) file descriptor */
gboolean GetHeaderNew (ObitOTF *outData, char *outfile, gchar *infile, 
		       ObitInfoList *myInput, ObitErr *err);
/* Get (new) data */
void GetDataNew (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, ObitErr *err);
/* Get (old) file descriptor */
gboolean GetHeaderOld (ObitOTF *outData, char *outfile, gchar *infile, 
		       ObitInfoList *myInput, ObitErr *err);
/* Get (old) data */
void GetDataOld (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, ObitErr *err);
/* Get target id */
olong GetTarget (ObitOTF *outData, gboolean isNewFile, gchar *name, 
		odouble ra, odouble dec, odouble equinox, ObitErr *err);
/* Initialize Index table */
void InitScan (ObitOTF *outData, gboolean isNewFile, ofloat scan, ofloat target, 
	       ObitErr *err);
/* Undate scan info */
void SetScan (ObitOTF *outData, odouble startTime, odouble endTime, 
	      olong startRec, olong endRec, ObitErr *err);
/* Reordering indices */
void SetReorderNew (olong *reorder);
void SetReorderOld (olong *reorder);


/* Program globals */
gchar *pgmName = "CCBOTF";       /* Program name */
gchar *input_file  = "CCBOTF.inp";   /* File with program inputs */
gchar *output_file = "CCBOTF.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
/* gchar **FITSdirs=NULL; List of FITS data directories */
gchar DataRoot[128]; /* Root directory of input data */
gchar Name[48];   /* Target name */
odouble targRA;   /* Target RA in deg */
odouble targDec;  /* Target Dec in deg */
odouble Equinox;  /* Equinox of targRA, targDec */
olong nAntTime;   /* number of antenna time samples */
odouble *AntDMJD; /* Array of antenna times */
odouble *AntRA;   /* Array of Antenna RA J2000 values */
odouble *AntDec;  /* Array of Antenna RA J2000 values */
odouble refMJD;   /* reference Julian date */
odouble integTime;/* Integration time in days */
odouble startTime;/* Start time of Scan in days */
odouble endTime;  /* End time of Scan in days */
ofloat target;    /* target number */
ofloat scan;      /* scan number */
olong  startRec;  /* First record number (1-rel) in scan */
olong  endRec;    /* End record number (1-rel) in scan */
olong  nfeed=2;   /* Number of feeds */
olong  nchan=4;   /* Number of frequencies */
olong  nstok=2;   /* Number of Stokes */
olong  nstate=4;  /* Number of Switching states */
olong  ndetect=64;/* Number of detectors */
olong  nport=16;  /* Number of ports in use */
gboolean isBS=TRUE; /* True if data beam switched */
ofloat dataNorm[64];/* Normalization per detector/state */
gchar  poln[64];  /* Polarization states per detector, 'R', 'L', 'X', 'Y' */
ofloat deltaFreq; /* Channel increment */
ofloat refPixel;  /* Frequency reference pixel */
odouble refFrequency; /* reference frequency (Hz) */
ObitGBTIFInfo*       IFdata     = NULL; /* IF information structure */
ObitGBTCCBStateInfo* StateData  = NULL; /* State information structure */
ObitGBTCCBPortInfo* PortData    = NULL; /* Port information structure */
ObitGBTBeamOffInfo* BeamOffData = NULL; /* Beam offset information structure */
ofloat azOff, elOff; /* az and el offsets to apply to pointing positions */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*    Read CalTech Continuum Backend data to an OTF dataset               */
/*----------------------------------------------------------------------- */
{
  olong  i, disk, nrec, ierr=0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *outData= NULL;
  ObitIOAccess access;
  ObitErr *err= NULL;
  gchar *FITSdir[] = {"./", NULL};
  gchar infile[128], outOTF[128];
  ObitInfoType type;
  gboolean isNewFile, isNew;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  err = newObitErr();

  /* Startup - parse command line */
  ierr = 0;
  myInput = CCBOTFin (argc, argv, err);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Get inputs */
  /* input FITS file name */
  for (i=0; i<128; i++) infile[i] = 0;
  ObitInfoListGet(myInput, "Scan", &type, dim, infile, err);
  infile[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(infile);  /* Trim trailing blanks */

  /* output FITS file name */
  for (i=0; i<128; i++) outOTF[i] = 0;
  ObitInfoListGet(myInput, "outOTF", &type, dim, outOTF, err);
  outOTF[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(outOTF);  /* Trim trailing blanks */

  /* Get input data file root */
  for (i=0; i<128; i++) DataRoot[i] = 0;
  ObitInfoListGet(myInput, "DataRoot", &type, dim, DataRoot, err);
  DataRoot[dim[0]] = 0;  /* null terminate */
  ObitTrimTrail(DataRoot);  /* Trim trailing blanks */

  /* Initialize Obit */
  FITSdir[1] =  DataRoot;  /* Input FITS data directory */
  nFITS = 2;
  mySystem = ObitSystemStartup (pgmName, 1, 0, 0, NULL, nFITS, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Create ObitOTF for data */
  outData = newObitOTF("Output data");
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Which version of CCB? */
  isNew = isNewOld (infile, err);

  /* show any errors */
  /* Get header info, array geometry, initialize output if necessary */
  if (isNew) 
    isNewFile = GetHeaderNew (outData, outOTF, infile, myInput, err);
  else
    isNewFile = GetHeaderOld (outData, outOTF, infile, myInput, err);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;
  
  /* Get Telescope position array, Prime and secondary focus different */
  if (refFrequency>1.0e9)
    GetAntennaGR (infile, myInput, err);
  else  /* Low frequencies at prime focus */
    GetAntennaPF (infile, myInput, err);
  /* Say what went wrong if error */
  if (err->error) 
     Obit_log_error(err, OBIT_Error, "Error reading Antenna file for scan %s", infile);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Define output, I/O size */
  disk = 1;
  nrec = StateData->nCCBState;
  ObitOTFSetFITS(outData, nrec, disk, outOTF,err);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;
   
  /* Open output OTF */
  if (isNewFile) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  if ((ObitOTFOpen (outData, access, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", outOTF);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Get target id number */
  target = (ofloat)GetTarget (outData, isNewFile, Name, targRA, targDec, Equinox, err);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Initialize scan in Index table */
  InitScan (outData, isNewFile, scan, target, err);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;

  /* convert data  */
  if (isNew) 
    GetDataNew (outData, infile, myInput, err);
  else
    GetDataOld (outData, infile, myInput, err);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;

  /* Update index table */
  SetScan (outData, startTime, endTime, startRec, endRec, err);
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;
   
  /* Close */
  if ((ObitOTFClose (outData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing output file");
  
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);   if (ierr!=0) goto exit;
   
   /* Shutdown Obit */
 exit:
   ObitReturnDumpRetCode (ierr, output_file, myOutput, err);  /* Final output */
   mySystem = ObitSystemShutdown (mySystem);
   
   /* cleanup */
   myInput = ObitInfoListUnref(myInput);  /* delete input list */
   outData = ObitUnref(outData);
   if (AntDMJD) g_free(AntDMJD);
   if (AntRA)   g_free(AntRA);
   if (AntDec)  g_free(AntDec);
   IFdata      = ObitGBTIFInfoUnref(IFdata); 
   BeamOffData = ObitGBTBeamOffInfoUnref(BeamOffData); 
   StateData   = ObitGBTCCBStateInfoUnref(StateData); 
   PortData    = ObitGBTCCBPortInfoUnref(PortData); 
 
   return ierr;
} /* end of main */

ObitInfoList* CCBOTFin (int argc, char **argv, ObitErr *err)
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
  olong ax, itemp;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *arg;
  gchar *scan_name ="unspecified";
  gboolean init=FALSE;
  gchar *strTemp;
  ObitInfoList* list;
  gchar *routine = "CCBOTFin";

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

    } else if (strcmp(arg, "-DataRoot") == 0){ /* Data root directory */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "DataRoot", OBIT_string, dim, strTemp);

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
  /*ObitInfoListGet(list, "pgmNumber", &type, dim, &pgmNumber, err);*/
  /*ObitInfoListGet(list, "nFITS",     &type, dim, &nFITS,     err);*/
  if (err->error) Obit_traceback_val (err, routine, routine, list);

  /* Initialize output */
  ObitReturnDumpRetCode (-999, output_file, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, routine, list);

 return list;
} /* end CCBOTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: CCBOTF -input file [-scan date/time]\n");
    fprintf(stderr, "Convert an GBT CCB file format to Obit/OTF\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def CCBOTF.in\n");
    fprintf(stderr, "  -output output parameter file, def CCBOTF.out\n");
    fprintf(stderr, "  -Scan date/time used for form scan FITS file names\n");
    fprintf(stderr, "  -DataRoot GBT Archive Data root directory\n");
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

void GetAntennaGR (gchar *infile, ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get antenna positions and leave in globals, Secondary focus (>1 GHz)  */
/*      outData  Output OTF object                                        */
/*      infile   root of input file names                                 */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  gchar FullFile[128];
  olong irow;
  olong disk, ver, nrow;
  gchar *tab;
  ObitTableGBTANTPOSGR* table;
  ObitTableGBTANTPOSGRRow* row;
  ObitIOCode retCode;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(infile!=NULL);
  g_assert(myInput!=NULL);

  /* get full file name */
  sprintf (FullFile,"Antenna/%s.fits", infile);
  
  /* Create table structure */
  table = newObitTableGBTANTPOSGR("Antenna");
  if (err->error) return;

  /* Setup */
  disk = 2;  /* Input data directory */
  tab = "ANTPOSGR";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(table,disk,FullFile,tab,ver,nrow,err);

  /* Open */
  retCode = ObitTableGBTANTPOSGROpen (table, OBIT_IO_ReadOnly, err);
  if (err->error) return;

  /* Create Row structure */
  row = newObitTableGBTANTPOSGRRow (table);

  /* make sure there is data */
  if (table->myDesc->nrow<=0) {
     Obit_log_error(err, OBIT_Error, "No data in Antenna file for scan %s", infile);
     return;
 }

  /* Create arrays */
  nAntTime = table->myDesc->nrow;
  AntDMJD = g_malloc0(nAntTime*sizeof(odouble));
  AntRA   = g_malloc0(nAntTime*sizeof(odouble));
  AntDec  = g_malloc0(nAntTime*sizeof(odouble));

  /* Loop over table */
  for (irow = 1; irow<=nAntTime; irow++) {
    retCode = ObitTableGBTANTPOSGRReadRow (table, irow, row, err);
    if (err->error) return;
    AntDMJD[irow-1] = row->dmjd;
    AntRA[irow-1]   = row->raj2000;
    AntDec[irow-1]  = row->decj2000;
 } /* end loop over table */
  
  /* Close */
  retCode = ObitTableGBTANTPOSGRClose (table, err);
  if (err->error) return;

  /* Cleanup */
  table = ObitTableGBTANTPOSGRUnref(table);
  row = ObitTableGBTANTPOSGRRowUnref(row);
} /* end GetAntennaGR  */

void GetAntennaPF (gchar *infile, ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get antenna positions and leave in globals, Prime focus (<1 GHz)      */
/*      outData  Output OTF object                                        */
/*      infile   root of input file names                                 */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  gchar FullFile[128];
  olong irow;
  olong disk, ver, nrow;
  gchar *tab;
  ObitTableGBTANTPOSPF* table;
  ObitTableGBTANTPOSPFRow* row;
  ObitIOCode retCode;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(infile!=NULL);
  g_assert(myInput!=NULL);

  /* get full file name */
  sprintf (FullFile,"Antenna/%s.fits", infile);
  
  /* Create table structure */
  table = newObitTableGBTANTPOSPF("Antenna");
  if (err->error) return;

  /* Setup */
  disk = 2;  /* Input data directory */
  tab = "ANTPOSPF";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(table,disk,FullFile,tab,ver,nrow,err);

  /* Open */
  retCode = ObitTableGBTANTPOSPFOpen (table, OBIT_IO_ReadOnly, err);
  if (err->error) return;

  /* Create Row structure */
  row = newObitTableGBTANTPOSPFRow (table);

  /* make sure there is data */
  if (table->myDesc->nrow<=0) {
     Obit_log_error(err, OBIT_Error, "No data in Antenna file for scan %s", infile);
     return;
 }

  /* Create arrays */
  nAntTime = table->myDesc->nrow;
  AntDMJD = g_malloc0(nAntTime*sizeof(odouble));
  AntRA   = g_malloc0(nAntTime*sizeof(odouble));
  AntDec  = g_malloc0(nAntTime*sizeof(odouble));

  /* Loop over table */
  for (irow = 1; irow<=nAntTime; irow++) {
    retCode = ObitTableGBTANTPOSPFReadRow (table, irow, row, err);
    if (err->error) return;
    AntDMJD[irow-1] = row->dmjd;
    AntRA[irow-1]   = row->raj2000;
    AntDec[irow-1]  = row->decj2000;
 } /* end loop over table */
  
  /* Close */
  retCode = ObitTableGBTANTPOSPFClose (table, err);
  if (err->error) return;

  /* Cleanup */
  table = ObitTableGBTANTPOSPFUnref(table);
  row = ObitTableGBTANTPOSPFRowUnref(row);
} /* end GetAntennaPF  */

gboolean isNewOld (gchar *infile, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Determine if file is new or old CCB                                   */
/*   Reads GBTCCBPortInfo and if there are fewer than 16 ports            */
/*   then this is new CCB data                                            */
/*   Input:                                                               */
/*      infile   root of input file names                                 */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  gboolean out = FALSE;
  olong disk;
  gchar *routine = "isNewOld";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;
  g_assert(infile!=NULL);

  /* Get Port info data */
  disk = 2;  /* input "disk" */
  PortData = newObitGBTCCBPortInfoValue("Port", disk, infile, err);
  if (err->error) Obit_traceback_val (err, routine, infile, out);
  out =  PortData->nCCBPort<12;  /* New data only has 8 ports, old 16 */

  /* Which kind of data? */
  if (out)
    Obit_log_error(err, OBIT_InfoErr, "Data from New CCB");
  else
    Obit_log_error(err, OBIT_InfoErr, "Data from Old CCB");
  ObitErrLog(err); /* Tell it */

  return out;
} /* end isNewOld */

gboolean GetHeaderNew (ObitOTF *outData, gchar *outfile, gchar *infile, 
		       ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get header information from scan header files                         */
/*  New (Post 2006?) CCB                                                  */
/*  Returns TRUE if the file is just created                              */
/*   Input:                                                               */
/*      outData  Output OTF object                                        */
/*      outfile  Output file name                                         */
/*      infile   root of input file names                                 */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  gboolean out = FALSE;
  ObitOTFDesc *desc;
  ObitOTFArrayGeom *geom;
  olong numberDetect, ncol, ncopy;
  long intgdp, iscan, smpperst, aisxl;
  double dtemp;
  float ftemp;
  gchar Date[48], *fullname=NULL, FullFile[128], commnt[81];
  odouble SiteLat, SiteLong, SiteElev, JD, T, GMST0, GSTIAT0, nu, e2;
  odouble aEarth = 6378137.0; /* Semi major axis of Earth */
  odouble flat   = 1.0 / 298.257223563; /* inverse of flattening of Earth */ 
  ofloat UT1UTC, PolarX, PolarY, adcsampt;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, j, disk, ierr = 0;
  fitsfile *fptr;
  gchar *routine = "GetHeaderNew";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitOTFIsA(outData));
  g_assert(outfile!=NULL);
  g_assert(infile!=NULL);
  g_assert(myInput!=NULL);

  desc = outData->myDesc;

  /* Data normalization factors */
  for (i=1; i<64; i++) dataNorm[i] = 0.0;
  ObitInfoListGetTest(myInput, "dataNorm", &type, dim, dataNorm);
  if (dataNorm[0]==0.0) dataNorm[0] = 1.0;
  for (i=1; i<64; i++) {
    if (dataNorm[i]==0.0) dataNorm[i] = dataNorm[0];
  }

  /* Get frequency setup */
  disk = 2;  /* input "disk" */
  IFdata = newObitGBTIFInfoValue("IFdata", "CCB26_40", disk, infile, err);
  /* Save data */
  nchan        = IFdata->nchan;
  deltaFreq    = IFdata->delta[0];
  refPixel     = IFdata->refPixel[0];
  refFrequency = 1.0e20; /* Order seems reversed - find min */
  for (i=0; i<IFdata->nIF; i++) {
    refFrequency = MIN (IFdata->refFrequency[i], refFrequency);
  }
  poln[0]      = 'X';  /* IF table incorrect? */
  poln[1]      = 'Y';

  /* Get Beam offset data */
  BeamOffData = newObitGBTBeamOffInfoValue("BeamOff", disk, infile, err);

  /* Get State info data */
  StateData = newObitGBTCCBStateInfoValue("State", disk, infile, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, out);
  nstate = StateData->nCCBState;

  /* Get Port info data */
  nport =  PortData->nCCBPort;

  /* Is this regular or beam switched? */
  ObitInfoListGetTest(myInput, "doBS",   &type, dim, &isBS);
  if (isBS) nstate = 1;  /* Will combine states */

  /* Set Feed Array Geometry */
  /* create Array geometry with enough (R,L) elements */
  numberDetect =  nport * nstate;
  geom = ObitOTFArrayGeomCreate (numberDetect);
  geom->azOffset = g_realloc(geom->azOffset, numberDetect*sizeof(ofloat));
  geom->elOffset = g_realloc(geom->elOffset, numberDetect*sizeof(ofloat));

  /* Create detector table - each output data stream has one  */
  j = 1;
  nstok = 1;
  for (i=0; i<numberDetect; i+=2) {
    /* Only one poln per feed */
    j = 1-j;  /* real feed number */
    geom->azOffset[i] =  BeamOffData->xeloff[j];
    geom->elOffset[i] =  BeamOffData->eloff[j];
    j = 1-j;  /* real feed number */
    geom->azOffset[i+1] =  BeamOffData->xeloff[j];
    geom->elOffset[i+1] =  BeamOffData->eloff[j];
    /* Set polarizations */
    poln[i]   = 'X';  /* IF table incorrect? */
    poln[i+1] = 'Y';
  }

  /* rerefer all the offsets to feed 1 */
  /* Correction to pointing */
  azOff = -geom->azOffset[0];
  elOff = -geom->elOffset[0];
  for (i=0; i<numberDetect; i++) {
    geom->azOffset[i] += azOff;
    geom->elOffset[i] += elOff;
  }

  
  /* Source info from GO file */
  /* get full file name */
  sprintf (FullFile,"%sGO/%s.fits", DataRoot, infile);

  /* Get header keywords from naked cfitsio */
   if ( fits_open_file(&fptr, FullFile, READONLY, &ierr) ) {
      Obit_log_error(err, OBIT_Error, "ERROR %d opening input FITS file %s", 
		     ierr, FullFile);
      return out;
   }

   /* Read main file keywords */
   iscan = 1;
   fits_read_key_lng (fptr, "SCAN", &iscan, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   /* Target name */
   fits_read_key_str (fptr, "OBJECT", Name, commnt, &ierr);

   /* Target RA */
   dtemp = 0.0;
   fits_read_key_dbl (fptr, "RA      ", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   /* Now in degrees targRA *= 15.0;   to degrees */
   targRA = (odouble)dtemp;

   dtemp = 0.0;
   fits_read_key_dbl (fptr, "DEC     ", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   targDec = (odouble)dtemp;

   /* Target Equinox */
   dtemp = 2000.0;
   fits_read_key_dbl (fptr, "EQUINOX ", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   Equinox = (odouble)dtemp;

   /* Close FITS file */
   fits_close_file (fptr, &ierr);
   if (ierr!=0) {
     Obit_log_error(err, OBIT_Error, "ERROR reading input FITS file %s", FullFile);
     return out;
   }

   /* Other information - Get from Antenna file */
  /* get full file name */
  sprintf (FullFile,"%sAntenna/%s.fits", DataRoot, infile);

  /* Get header keywords from naked cfitsio */
   if ( fits_open_file(&fptr, FullFile, READONLY, &ierr) ) {
      Obit_log_error(err, OBIT_Error, "ERROR %d opening input FITS file %s", 
		     ierr, FullFile);
      return out;
   }

   /* Read main file keywords */
   ftemp = 0.0;
   fits_read_key_flt (fptr, "DELTAUTC", &ftemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   UT1UTC = (ofloat)ftemp;

   ftemp = 0.0;
   fits_read_key_flt (fptr, "IERSPMX", &ftemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   PolarX = (ofloat)ftemp;

   ftemp = 0.0;
   fits_read_key_flt (fptr, "IERSPMY", &ftemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   PolarY = (ofloat)ftemp;

   dtemp = 0.0;
   fits_read_key_dbl (fptr, "SITELAT", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   SiteLat = (odouble)dtemp;

   dtemp = 0.0;
   fits_read_key_dbl (fptr, "SITELONG", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   SiteLong = (odouble)dtemp;
   
   dtemp = 0.0;
   fits_read_key_dbl (fptr, "SITEELEV", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   SiteElev = (odouble)dtemp;

   fits_read_key_str (fptr, "DATE-OBS", Date, commnt, &ierr);
 
   /* Close FITS file */
   fits_close_file (fptr, &ierr);
   if (ierr!=0) {
     Obit_log_error(err, OBIT_Error, "ERROR reading input FITS file %s", FullFile);
     return out;
   }

   /* Convert to JD */
  ObitOTFDescDate2JD (Date, &JD);

  /* GST at IAT=0 at 0h on reference date (deg)*/
  /* Tropical century from jan 0.5, 2000 */
  T = (JD - 2451545.0) / 36525.0;

  /* GMST at IAT=0 in radians */
  GMST0 = ((((((-6.2e-6 * T) + 0.093104) * T) + 8640184.812866) * T + 24110.54841) 
	   * 2.0 * G_PI / 86400.0);
  /* to degrees */
  GSTIAT0 = RAD2DG * fmod(GMST0, (2.0*G_PI));

   ncopy = strlen (Date);
   if (ncopy>10) ncopy = 10;
  for (i=0; i<ncopy; i++) geom->RefDate[i] = Date[i]; geom->RefDate[i]=0;
  geom->TimeSys[0] = 'U'; geom->TimeSys[1] = 'T';geom->TimeSys[2] = 'C';geom->TimeSys[3] = 0;

  /* Conversion of geocentric to geodetic */
  e2 = 2.0*flat - flat*flat;
  nu = aEarth / sqrt (1.0 - e2 * sin(DG2RAD*SiteLat) * sin(DG2RAD*SiteLat));
  geom->TeleX =  (nu + SiteElev) * cos(DG2RAD*SiteLat) * cos(DG2RAD*SiteLong);
  geom->TeleY = -(nu + SiteElev) * cos(DG2RAD*SiteLat) * sin(DG2RAD*SiteLong);
  geom->TeleZ =  ((1.0 - e2) * nu +  SiteElev) * sin(DG2RAD*SiteLat);


  geom->DegDay  = 3.6098564497330e+02;
  geom->GSTiat0 = GSTIAT0;
  geom->PolarX  = PolarX;
  geom->PolarY  = PolarY;
  geom->ut1Utc  = UT1UTC;
  geom->dataUtc = 0.0;
  geom->iatUtc  = 0.0;

  /* Compute some useful terms */
  /* telescope latitude in radians */
  geom->lat = SiteLat * DG2RAD;
  /* telescope East longitude in radians */
  geom->lon = -SiteLong * DG2RAD;
  /* LST at iat0 in radians */
  geom->LSTiat0 = geom->GSTiat0*1.74533e-2 + geom->lon;
  /* Earth rotation rate in rad/day */
  geom->RadDay = geom->DegDay*1.74533e-2;
  /* Data - IAT in days */
  geom->dataIat = (geom->dataUtc - geom->iatUtc) / 86400.0;

  /* Attach Array geometry to OTF */
  outData->geom = ObitOTFArrayGeomRef(geom);
  geom = ObitOTFArrayGeomUnref(geom);

  /* Other information from CCB file */
  /* get full file name */
  sprintf (FullFile,"%sCCB26_40/%s.fits", DataRoot, infile);

  /* Get header keywords from naked cfitsio */
   if ( fits_open_file(&fptr, FullFile, READONLY, &ierr) ) {
      Obit_log_error(err, OBIT_Error, "ERROR %d opening input FITS file %s", 
		     ierr, FullFile);
      return out;
   }

   /* Make sure data in "CCB" mode */

   /* The header is current defective and doesn't have this 
      fits_read_key_str (fptr, "BACKEND", backend, commnt, &ierr);
      if (strncmp (backend, "CCB26_40",3)) {
      Obit_log_error(err, OBIT_Error, "ERROR backend %s NOT CCB", backend);
      return out;
      } */

   /* Read main file keywords */
   /* New CCB defined to have 2 feeds x 1 poln, with four sets of detectors, 
      beamswitching will combine the four states */
   ndetect = nport * nstate;

   /* Have to work out integration time */
   ftemp = 1.0e-07;
   fits_read_key_flt (fptr, "ADCSAMPT", &ftemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   adcsampt = (ofloat)ftemp;

   smpperst = 1250;
   fits_read_key_lng (fptr, "SMPPERST", &smpperst, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   intgdp   = 10;
   fits_read_key_lng (fptr, "INTGPD", &intgdp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   /* Cal flag */
   aisxl    = 0;
   fits_read_key_lng (fptr, "AISXL", &aisxl, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   /* Integration time in seconds */
   integTime = adcsampt * smpperst * intgdp * nstate;
   /* Convert from seconds to days */
   integTime /= 86400.0;

   /* Close FITS file */
   fits_close_file (fptr, &ierr);
   if (ierr!=0) {
     Obit_log_error(err, OBIT_Error, "ERROR reading input FITS file %s", FullFile);
     return out;
   }

  /* initialize globals */
  refMJD = 0.0;
  target = 0.0;
  scan   = (ofloat)iscan;

   
  /* Does the input file exist?  BTW, I only do FITS */
  disk = 1;  /* input "disk" */
  fullname = ObitFITSFilename (disk, outfile, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, out);


  /* check is output file exists, if not, initialize */
  if (!ObitFileExist (fullname, err)) {
    if (err->error>0) Obit_log_error(err, OBIT_Error, "ERROR testing file %s", fullname);

    /* Initialize  Descriptor */
    /* Approximate Beam size */
    desc->beamSize = 2.0e8 / refFrequency;

    /* diameter */
    desc->diameter = 100.0;  /* GBT diameter */

    /* new file */
    out = TRUE;
    
    strncpy (desc->object, "Sky", OTFLEN_VALUE);
    strncpy (desc->teles,  "GBT       ", OTFLEN_VALUE);
    strncpy (desc->origin, "Obit ", OTFLEN_VALUE);
    desc->isort[0] = 'T';  /* Time ordered */
    ncol = 0;
    
    desc->JDObs = 0.0;
    desc->epoch = 2000.0;
    desc->equinox = 2000.0;
    strncpy (desc->bunit,  "ADU     ", OTFLEN_VALUE);
    strncpy (desc->obsdat, Date, OTFLEN_VALUE);
    
    /* Time */
    strncpy (desc->colType[ncol], "TIME    ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DAYS    ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Integration time */
    strncpy (desc->colType[ncol], "TIME_INT", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DAYS    ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Target index */
    strncpy (desc->colType[ncol], "TARGET  ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "        ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Scan index */
    strncpy (desc->colType[ncol], "SCAN    ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "        ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Pointing RA */
    strncpy (desc->colType[ncol], "RA      ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DEGREE  ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Pointing Dec */
    strncpy (desc->colType[ncol], "DEC     ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DEGREE  ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Rotation of array on sky */
    strncpy (desc->colType[ncol], "ROTATE  ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DEGREE  ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Cal on? */
    strncpy (desc->colType[ncol], "CAL     ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "        ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Data - MUST be last column */
    desc->numDesc = ncol;
    strncpy (desc->colType[ncol], "DATA    ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "COUNTS  ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = numberDetect*2;  /* With weights */
    ncol++;
    
    desc->ncol = ncol;
    
    /* Data array descriptors */
    desc->naxis = 0;

     /* Data-Wt axis */
    desc->inaxes[desc->naxis] = 2;
    strncpy (desc->ctype[desc->naxis], "DATAWT", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = 1.0;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = 1.0;
    desc->naxis++;

   /* Stokes axis */
    desc->inaxes[desc->naxis] = nstok;
    strncpy (desc->ctype[desc->naxis], "STOKES  ", OTFLEN_KEYWORD);
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    /* Set appropritate Stokes Type */
    if ((poln[0]=='R') || (poln[0]=='L')) { /* R, L */
      desc->cdelt[desc->naxis] = -1.0;
      desc->crval[desc->naxis] = -1.0;
    } else { /* Must be X, Y */
      desc->cdelt[desc->naxis] = -1.0;
      desc->crval[desc->naxis] = -5.0;
    }
    desc->naxis++;

    /* Feed axis */
    desc->inaxes[desc->naxis] = nfeed;
    strncpy (desc->ctype[desc->naxis], "FEED", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = 1.0;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = 1.0;
    desc->naxis++;

    /* Frequency axis */
    desc->inaxes[desc->naxis] = nchan;
    strncpy (desc->ctype[desc->naxis], "FREQ    ", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = deltaFreq;
    desc->crpix[desc->naxis] = refPixel;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = refFrequency;
    desc->naxis++;

    /* State axis */
    desc->inaxes[desc->naxis] = nstate;
    strncpy (desc->ctype[desc->naxis], "STATE   ", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = 1.0;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = 1.0;
    desc->naxis++;

    /* Reference Modified Julian date */
    refMJD = JD - 2400000.5;

    /* Mark as CCB data */
    desc->OTFType = OBIT_GBTOTF_CCB;
    /* end initialize descriptor */
  }

  /* cleanup */
  if (fullname) g_free(fullname);

  /* Index the descriptor */
  ObitOTFDescIndex (desc);

  return out;
} /* end GetHeaderNew */

void GetDataNew (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, 
		 ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Read data from GB FITS files and write                                */
/*  New (Post 2006?) CCB                                                  */
/*  Assumptions:                                                          */
/*    The following uses state information read from the FITS files to    */
/*  decide what data is what.                                             */
/*  If the data is beam switched, a second set of "detectors is used      */
/*  for the "ref" state where the feeds are swapped and therefore have    */
/*  different amplifier chains and therefore different gains              */
/*  Cal column = positive if Cal A on, negative if cal B on               */
/*      outData  Output OTF object                                        */
/*      infile   Scan part of input file name                             */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  gchar FullFile[128];
  olong irow;
  olong i, j, k, state, port, disk, ver, nrow;
  olong iD, reorder[64], incdatawt, scanNo;
  gchar *tab;
  ofloat *data, ra, dec, cTime, sum, fblank = ObitMagicF();
  ofloat offTime, ldata[64];
  odouble offTimeD, JD;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gboolean valid;
  ObitTableGBTCCBDATA* table;
  ObitTableGBTCCBDATARow* row;
  ObitOTFDesc *desc;
  ObitOTFArrayGeom *geom;
  ObitIOCode retCode;
  gchar *routine = "GetDataNew";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(infile!=NULL);
  g_assert(myInput!=NULL);

  /* Consistency checks */
  desc = outData->myDesc;
  Obit_return_if_fail((nstate==desc->inaxes[desc->jlocstate]), err,
		       "%s: Request beamswitch states incompatible %d !=  %d", 
		      routine, nstate, desc->inaxes[desc->jlocstate]);
  Obit_return_if_fail((nfeed==desc->inaxes[desc->jlocfeed]), err,
		       "%s: Request number feeds incompatible %d !=  %d", 
		      routine, nfeed, desc->inaxes[desc->jlocfeed]);
  Obit_return_if_fail((nstok==desc->inaxes[desc->jlocs]), err,
		       "%s: Request number Poln incompatible %d !=  %d", 
		      routine, nstok, desc->inaxes[desc->jlocs]);
  Obit_return_if_fail((nchan==desc->inaxes[desc->jlocf]), err,
		       "%s: Request number freq. incompatible %d !=  %d", 
		      routine, nchan, desc->inaxes[desc->jlocf]);
  
  /* Get reference MJD, convert ref date  to JD */
  ObitOTFDescDate2JD (outData->geom->RefDate, &JD);
  /* Reference Modified Julian date */
  refMJD = JD - 2400000.5;

  /* get full file name */
  sprintf (FullFile,"CCB26_40/%s.fits", infile);
  
  /* Create table structure */
  table = newObitTableGBTCCBDATA("Data");
  if (err->error) return;

  /* Setup */
  disk = 2;
  tab = "DATA";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(table,disk,FullFile,tab,ver,nrow,err);

  /* Set reordering coefficients to change order in GBT table to OTF */
  SetReorderNew(reorder);

  /* Time offset */
  offTime = 0.0;
  ObitInfoListGetTest(myInput, "offTime",  &type, dim, &offTime);

  /* Correction to time in days */
  offTimeD = offTime / 86400.0;

  /* Get correction from start time to center */
  cTime = 0.5 * integTime;

  /* Scan no. */
  scanNo = 0;
  ObitInfoListGetTest(myInput, "scanNo",  &type, dim, &scanNo);
  if (scanNo!=0) scan = (ofloat)scanNo;
  Obit_log_error(err, OBIT_InfoErr, "Adding scan %s %6.0f to OTF file", 
		 infile, scan);
  ObitErrLog(err);

  /* Open */
  retCode = ObitTableGBTCCBDATAOpen (table, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Create Row structure */
  row = newObitTableGBTCCBDATARow (table);

  desc = outData->myDesc;
  geom = outData->geom;
  incdatawt = desc->incdatawt; /* increment in data-wt axis */
  
  /* Write at end */
  desc->firstRec = desc->nrecord+1;
  startRec = desc->firstRec;  /* First record in scan */
  startTime = -1.0e20;        /* dummy srart time */
  
  /* Loop over table */
  for (irow = 1; irow<=table->myDesc->nrow; irow++) {
    retCode = ObitTableGBTCCBDATAReadRow (table, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* first time in scan */
    if (startTime<-1000.0) startTime = row->dmjd + cTime - refMJD;
    data = outData->buffer;     /* Output data array */
    
    /* Check validity */
    if (row->usable) {
    
      /* Fill record entry in data */
      data[desc->iloct]   = row->dmjd + cTime  - refMJD;    /* Time (days) */
      /* Correction for GBT timing error  */
      data[desc->iloct] += offTimeD;
      /* Get position */
      GetPoint (row->dmjd+offTimeD+cTime, &ra, &dec);
      data[desc->ilocti]  = cTime*2.0;    /* time interval (days) */
      data[desc->iloctar] = target;       /* target number */
      data[desc->ilocscan]= scan;         /* scan number */
      data[desc->iloccal] = 0.0;          /* Cal? */
      /* Two possible cal states, sig and ref */
      if (row->cala) data[desc->iloccal] =  1.0;
      if (row->calb) data[desc->iloccal] = -1.0;
      /* If both on bag it */
      if ((row->cala) && (row->calb)) continue;

      /* Parallactic angle (deg) */
      data[desc->ilocrot] = ObitOTFArrayGeomParAng(geom, data[desc->iloct], ra, dec);
      /* correction to position */
      ObitOTFArrayGeomCorrPoint(azOff, elOff, data[desc->ilocrot], &ra, &dec);
      data[desc->ilocra]  = ra;           /* RA  in deg. */
      data[desc->ilocdec] = dec;          /* Dec in deg.*/
      
      /* Init data to blank */	
      for (i=0; i<ndetect; i++) {
	data[desc->ilocdata+i*incdatawt]   = fblank;
	data[desc->ilocdata+i*incdatawt+1] = 0.0;
      }

      /* Data normalization */
      for (iD=0; iD<table->myDesc->repeat[table->dataCol]; iD++) {
	/* Check validity */
	state = iD / nport;
	port  = iD - state*nport;
	j = PortData->slave[port];
	valid = (!row->ovrflow[iD]) && (row->slaveok[j]);
	if (valid) ldata[iD] = row->data[iD] * dataNorm[reorder[iD]];
	else  ldata[iD] = 0.0;
      }
      
      /* Beamswitching? */
      if (isBS) {
	/* Loop over data */
	for (iD=0; iD<ndetect; iD++) {
	  j = PortData->slave[iD];
	  sum = 0.0;
	  valid =  row->slaveok[j];
	  /* Loop over states */
	  for (k=0; k<StateData->nCCBState; k++) {
	    /* Check validity */
	    valid = valid && (!row->ovrflow[iD+k*nport]);
	    /* Accumulate */
	    if (StateData->phia[k]==StateData->phib[k]) 
	      sum += ldata[iD+k*nport];
	    else sum -= ldata[iD+k*nport];
	  }
	  if (valid) {
	    data[desc->ilocdata+reorder[iD]*incdatawt]   = 0.5 * sum;
	    data[desc->ilocdata+reorder[iD]*incdatawt+1] = 1.0;
	  }
	} /* end detector loop */
      } else { /* Pass all states */
	/* Loop over data */
	for (iD=0; iD<ndetect; iD++) {
	  /* Check validity */
	  state = iD / nport;
	  port  = iD - state*nport;
	  j = PortData->slave[port];
	  valid = (!row->ovrflow[iD]) && (row->slaveok[j]);
	  if (valid) {
	    data[desc->ilocdata+reorder[iD]*incdatawt]   = ldata[iD];
	    data[desc->ilocdata+reorder[iD]*incdatawt+1] = 1.0;
	  }
	}
      } /* end if beamswitch */
      
      /* set number of records */
      desc->numRecBuff = 1;
      
      /* Write output  */
      if ((ObitOTFWrite (outData, NULL, err) != OBIT_IO_OK) || (err->error>0))
	Obit_log_error(err, OBIT_Error, "ERROR writing output Table file");
      if (err->error) Obit_traceback_msg (err, routine, outData->name);
    } /* End valid row */
  } /* end loop over table */
  
  /* Get end times and record numbers */
  desc->firstRec = desc->nrecord+1;
  endTime = row->dmjd + cTime - refMJD;
  endRec  = desc->nrecord;  /* Last record in scan */

  /* Close */
  retCode = ObitTableGBTCCBDATAClose (table, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Cleanup */
  table = ObitTableGBTCCBDATAUnref(table);
  row = ObitTableGBTCCBDATARowUnref(row);


} /* end GetDataNew  */

gboolean GetHeaderOld (ObitOTF *outData, gchar *outfile, gchar *infile, 
		       ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get header information from scan header files                         */
/*  Old (2006?) CCB                                                       */
/*  Returns TRUE if the file is just created                              */
/*   Input:                                                               */
/*      outData  Output OTF object                                        */
/*      outfile  Output file name                                         */
/*      infile   root of input file names                                 */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  gboolean out = FALSE;
  ObitOTFDesc *desc;
  ObitOTFArrayGeom *geom;
  olong numberDetect, ncol, ncopy;
  gchar Date[48], *fullname=NULL, FullFile[128], commnt[81];
  odouble SiteLat, SiteLong, SiteElev, JD, T, GMST0, GSTIAT0, nu, e2;
  odouble aEarth = 6378137.0; /* Semi major axis of Earth */
  odouble flat   = 1.0 / 298.257223563; /* inverse of flattening of Earth */ 
  ofloat UT1UTC, PolarX, PolarY, adcsampt;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  long intgdp, iscan, smpperst, aisxl;
  double dtemp;
  float ftemp;
  olong i, j, disk, ierr = 0;
  fitsfile *fptr;
  gchar *routine = "GetHeaderOld";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitOTFIsA(outData));
  g_assert(outfile!=NULL);
  g_assert(infile!=NULL);
  g_assert(myInput!=NULL);

  desc = outData->myDesc;

  /* Data normalization factors */
  for (i=1; i<64; i++) dataNorm[i] = 0.0;
  ObitInfoListGetTest(myInput, "dataNorm", &type, dim, dataNorm);
  if (dataNorm[0]==0.0) dataNorm[0] = 1.0;
  for (i=1; i<64; i++) {
    if (dataNorm[i]==0.0) dataNorm[i] = dataNorm[0];
  }

  /* Get frequency setup */
  disk = 2;  /* input "disk" */
  IFdata = newObitGBTIFInfoValue("IFdata", "CCB26_40", disk, infile, err);
  /* Save data */
  nchan        = IFdata->nchan;
  deltaFreq    = IFdata->delta[0];
  refPixel     = IFdata->refPixel[0];
  refFrequency = 1.0e20; /* Order seems reversed - find min */
  for (i=0; i<IFdata->nIF; i++) {
    refFrequency = MIN (IFdata->refFrequency[i], refFrequency);
  }
  poln[0]      = IFdata->poln[0];
  poln[1]      = IFdata->poln[1]; 

  /* Get Beam offset data */
  BeamOffData = newObitGBTBeamOffInfoValue("BeamOff", disk, infile, err);

  /* Get State info data */
  StateData = newObitGBTCCBStateInfoValue("State", disk, infile, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, out);
  nstate = StateData->nCCBState;

  /* Get Port info data */
  nport =  PortData->nCCBPort;

  /* Is this regular or beam switched? */
  ObitInfoListGetTest(myInput, "doBS",   &type, dim, &isBS);
  if (isBS) nstate = 1;  /* Will combine states */

  /* Set Feed Array Geometry */
  /* create Array geometry with enough (R,L) elements */
  numberDetect =  nport * nstate;
  geom = ObitOTFArrayGeomCreate (numberDetect);
  geom->azOffset = g_realloc(geom->azOffset, numberDetect*sizeof(ofloat));
  geom->elOffset = g_realloc(geom->elOffset, numberDetect*sizeof(ofloat));

  /* Create detector table - each output data stream has one  */
  j = 1;
  for (i=0; i<numberDetect; i+=2) {
    j = 1-j;  /* real feed number */
    geom->azOffset[i] =  BeamOffData->xeloff[j];
    geom->elOffset[i] =  BeamOffData->eloff[j];
    /* Poln pair has the same offset */
    geom->azOffset[i+1] =  BeamOffData->xeloff[j];
    geom->elOffset[i+1] =  BeamOffData->eloff[j];
    /* Set polarizations */
    poln[i]   = IFdata->poln[j*2];
    poln[i+1] = IFdata->poln[j*2+1];
  }

  /* rerefer all the offsets to feed 1 */
  /* Correction to pointing */
  azOff = -geom->azOffset[0];
  elOff = -geom->elOffset[0];
  for (i=0; i<numberDetect; i++) {
    geom->azOffset[i] += azOff;
    geom->elOffset[i] += elOff;
  }

  
  /* Source info from GO file */
  /* get full file name */
  sprintf (FullFile,"%sGO/%s.fits", DataRoot, infile);

  /* Get header keywords from naked cfitsio */
   if ( fits_open_file(&fptr, FullFile, READONLY, &ierr) ) {
      Obit_log_error(err, OBIT_Error, "ERROR %d opening input FITS file %s", 
		     ierr, FullFile);
      return out;
   }

   /* Read main file keywords */
   iscan = 1;
   fits_read_key_lng (fptr, "SCAN", &iscan, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   /* Target name */
   fits_read_key_str (fptr, "OBJECT", Name, commnt, &ierr);

  /* Target RA */
   dtemp = 0.0;
   fits_read_key_dbl (fptr, "RA      ", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   /* Now in degrees targRA *= 15.0;   to degrees */
   targRA = (odouble)dtemp;

   dtemp = 0.0;
   fits_read_key_dbl (fptr, "DEC     ", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   targDec = (odouble)dtemp;

   /* Target Equinox */
   dtemp = 2000.0;
   fits_read_key_dbl (fptr, "EQUINOX ", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   Equinox = (odouble)dtemp;

   /* Close FITS file */
   fits_close_file (fptr, &ierr);
   if (ierr!=0) {
     Obit_log_error(err, OBIT_Error, "ERROR reading input FITS file %s", FullFile);
     return out;
   }

   /* Other information - Get from Antenna file */
  /* get full file name */
  sprintf (FullFile,"%sAntenna/%s.fits", DataRoot, infile);

  /* Get header keywords from naked cfitsio */
   if ( fits_open_file(&fptr, FullFile, READONLY, &ierr) ) {
      Obit_log_error(err, OBIT_Error, "ERROR %d opening input FITS file %s", 
		     ierr, FullFile);
      return out;
   }

   /* Read main file keywords */
   ftemp = 0.0;
   fits_read_key_flt (fptr, "DELTAUTC", &ftemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   UT1UTC = (ofloat)ftemp;

   ftemp = 0.0;
   fits_read_key_flt (fptr, "IERSPMX", &ftemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   PolarX = (ofloat)ftemp;

   ftemp = 0.0;
   fits_read_key_flt (fptr, "IERSPMY", &ftemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   PolarY = (ofloat)ftemp;

   dtemp = 0.0;
   fits_read_key_dbl (fptr, "SITELAT", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   SiteLat = (odouble)dtemp;

   dtemp = 0.0;
   fits_read_key_dbl (fptr, "SITELONG", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   SiteLong = (odouble)dtemp;
   
   dtemp = 0.0;
   fits_read_key_dbl (fptr, "SITEELEV", &dtemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   SiteElev = (odouble)dtemp;

   fits_read_key_str (fptr, "DATE-OBS", Date, commnt, &ierr);
 
   /* Close FITS file */
   fits_close_file (fptr, &ierr);
   if (ierr!=0) {
     Obit_log_error(err, OBIT_Error, "ERROR reading input FITS file %s", FullFile);
     return out;
   }

   /* Convert to JD */
  ObitOTFDescDate2JD (Date, &JD);

  /* GST at IAT=0 at 0h on reference date (deg)*/
  /* Tropical century from jan 0.5, 2000 */
  T = (JD - 2451545.0) / 36525.0;

  /* GMST at IAT=0 in radians */
  GMST0 = ((((((-6.2e-6 * T) + 0.093104) * T) + 8640184.812866) * T + 24110.54841) 
	   * 2.0 * G_PI / 86400.0);
  /* to degrees */
  GSTIAT0 = RAD2DG * fmod(GMST0, (2.0*G_PI));

  ncopy = strlen (Date);
  if (ncopy>10) ncopy = 10;
  for (i=0; i<ncopy; i++) geom->RefDate[i] = Date[i]; geom->RefDate[i]=0;
  geom->TimeSys[0] = 'U'; geom->TimeSys[1] = 'T';geom->TimeSys[2] = 'C';geom->TimeSys[3] = 0;

  /* Conversion of geocentric to geodetic */
  e2 = 2.0*flat - flat*flat;
  nu = aEarth / sqrt (1.0 - e2 * sin(DG2RAD*SiteLat) * sin(DG2RAD*SiteLat));
  geom->TeleX =  (nu + SiteElev) * cos(DG2RAD*SiteLat) * cos(DG2RAD*SiteLong);
  geom->TeleY = -(nu + SiteElev) * cos(DG2RAD*SiteLat) * sin(DG2RAD*SiteLong);
  geom->TeleZ =  ((1.0 - e2) * nu +  SiteElev) * sin(DG2RAD*SiteLat);


  geom->DegDay  = 3.6098564497330e+02;
  geom->GSTiat0 = GSTIAT0;
  geom->PolarX  = PolarX;
  geom->PolarY  = PolarY;
  geom->ut1Utc  = UT1UTC;
  geom->dataUtc = 0.0;
  geom->iatUtc  = 0.0;

  /* Compute some useful terms */
  /* telescope latitude in radians */
  geom->lat = SiteLat * DG2RAD;
  /* telescope East longitude in radians */
  geom->lon = -SiteLong * DG2RAD;
  /* LST at iat0 in radians */
  geom->LSTiat0 = geom->GSTiat0*1.74533e-2 + geom->lon;
  /* Earth rotation rate in rad/day */
  geom->RadDay = geom->DegDay*1.74533e-2;
  /* Data - IAT in days */
  geom->dataIat = (geom->dataUtc - geom->iatUtc) / 86400.0;

  /* Attach Array geometry to OTF */
  outData->geom = ObitOTFArrayGeomRef(geom);
  geom = ObitOTFArrayGeomUnref(geom);

  /* Other information from CCB file */
  /* get full file name */
  sprintf (FullFile,"%sCCB26_40/%s.fits", DataRoot, infile);

  /* Get header keywords from naked cfitsio */
   if ( fits_open_file(&fptr, FullFile, READONLY, &ierr) ) {
      Obit_log_error(err, OBIT_Error, "ERROR %d opening input FITS file %s", 
		     ierr, FullFile);
      return out;
   }

   /* Make sure data in "CCB" mode */

   /* The header is current defective and doesn't have this 
      fits_read_key_str (fptr, "BACKEND", backend, commnt, &ierr);
      if (strncmp (backend, "CCB26_40",3)) {
      Obit_log_error(err, OBIT_Error, "ERROR backend %s NOT CCB", backend);
      return out;
      } */

   /* Read main file keywords */
   /* CCB defined to have 2 feeds x 2 poln, with four sets of detectors, 
      beamswitching will combine the four states */
   ndetect = nport * nstate;

   /* Have to work out integration time */
   ftemp = 1.0e-07;
   fits_read_key_flt (fptr, "ADCSAMPT", &ftemp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;
   adcsampt = (ofloat)ftemp;

   smpperst = 1250;
   fits_read_key_lng (fptr, "SMPPERST", &smpperst, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   intgdp   = 10;
   fits_read_key_lng (fptr, "INTGPD", &intgdp, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   /* Cal flag */
   aisxl    = 0;
   fits_read_key_lng (fptr, "AISXL", &aisxl, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   /* Integration time in seconds */
   integTime = adcsampt * smpperst * intgdp * nstate;
   /* Convert from seconds to days */
   integTime /= 86400.0;

   /* Close FITS file */
   fits_close_file (fptr, &ierr);
   if (ierr!=0) {
     Obit_log_error(err, OBIT_Error, "ERROR reading input FITS file %s", FullFile);
     return out;
   }

  /* initialize globals */
  refMJD = 0.0;
  target = 0.0;
  scan   = (ofloat)iscan;

    /* Does the input file exist?  BTW, I only do FITS */
  disk = 1;  /* input "disk" */
  fullname = ObitFITSFilename (disk, outfile, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, out);


  /* check is output file exists, if not, initialize */
  if (!ObitFileExist (fullname, err)) {
    if (err->error>0) Obit_log_error(err, OBIT_Error, "ERROR testing file %s", fullname);

    /* Initialize  Descriptor */
    /* Approximate Beam size */
    desc->beamSize = 2.0e8 / refFrequency;

    /* diameter */
    desc->diameter = 100.0;  /* GBT diameter */

    /* new file */
    out = TRUE;
    
    strncpy (desc->object, "Sky", OTFLEN_VALUE);
    strncpy (desc->teles,  "GBT       ", OTFLEN_VALUE);
    strncpy (desc->origin, "Obit ", OTFLEN_VALUE);
    desc->isort[0] = 'T';  /* Time ordered */
    ncol = 0;
    
    desc->JDObs = 0.0;
    desc->epoch = 2000.0;
    desc->equinox = 2000.0;
    strncpy (desc->bunit,  "ADU     ", OTFLEN_VALUE);
    strncpy (desc->obsdat, Date, OTFLEN_VALUE);
    
    /* Time */
    strncpy (desc->colType[ncol], "TIME    ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DAYS    ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Integration time */
    strncpy (desc->colType[ncol], "TIME_INT", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DAYS    ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Target index */
    strncpy (desc->colType[ncol], "TARGET  ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "        ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Scan index */
    strncpy (desc->colType[ncol], "SCAN    ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "        ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Pointing RA */
    strncpy (desc->colType[ncol], "RA      ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DEGREE  ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Pointing Dec */
    strncpy (desc->colType[ncol], "DEC     ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DEGREE  ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Rotation of array on sky */
    strncpy (desc->colType[ncol], "ROTATE  ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "DEGREE  ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Cal on? */
    strncpy (desc->colType[ncol], "CAL     ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "        ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = 1;
    ncol++;
    
    /* Data - MUST be last column */
    desc->numDesc = ncol;
    strncpy (desc->colType[ncol], "DATA    ", OTFLEN_KEYWORD);
    strncpy (desc->colUnit[ncol], "COUNTS  ", OTFLEN_VALUE);
    desc->colRepeat[ncol] = numberDetect*2;  /* With weights */
    ncol++;
    
    desc->ncol = ncol;
    
    /* Data array descriptors */
    desc->naxis = 0;

     /* Data-Wt axis */
    desc->inaxes[desc->naxis] = 2;
    strncpy (desc->ctype[desc->naxis], "DATAWT", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = 1.0;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = 1.0;
    desc->naxis++;

   /* Stokes axis */
    desc->inaxes[desc->naxis] = nstok;
    strncpy (desc->ctype[desc->naxis], "STOKES  ", OTFLEN_KEYWORD);
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    /* Set appropritate Stokes Type */
    if ((poln[0]=='R') || (poln[0]=='L')) { /* R, L */
      desc->cdelt[desc->naxis] = -1.0;
      desc->crval[desc->naxis] = -1.0;
    } else { /* Must be X, Y */
      desc->cdelt[desc->naxis] = -1.0;
      desc->crval[desc->naxis] = -5.0;
    }
    desc->naxis++;

    /* Feed axis */
    desc->inaxes[desc->naxis] = nfeed;
    strncpy (desc->ctype[desc->naxis], "FEED", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = 1.0;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = 1.0;
    desc->naxis++;

    /* Frequency axis */
    desc->inaxes[desc->naxis] = nchan;
    strncpy (desc->ctype[desc->naxis], "FREQ    ", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = deltaFreq;
    desc->crpix[desc->naxis] = refPixel;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = refFrequency;
    desc->naxis++;

    /* State axis */
    desc->inaxes[desc->naxis] = nstate;
    strncpy (desc->ctype[desc->naxis], "STATE   ", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = 1.0;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = 1.0;
    desc->naxis++;

    /* Reference Modified Julian date */
    refMJD = JD - 2400000.5;

    /* Mark as CCB data */
    desc->OTFType = OBIT_GBTOTF_CCB;
    /* end initialize descriptor */
  }

  /* cleanup */
  if (fullname) g_free(fullname);

  /* Index the descriptor */
  ObitOTFDescIndex (desc);

  return out;
} /* end GetHeaderOld */

void GetDataOld (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, 
		 ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Read data from GB FITS files and write                                */
/*  Old (2006?) CCB                                                       */
/*  Assumptions:                                                          */
/*    The following uses state information read from the FITS files to    */
/*  decide what data is what.                                             */
/*  If the data is beam switched, a second set of "detectors is used      */
/*  for the "ref" state where the feeds are swapped and therefore have    */
/*  different amplifier chains and therefore different gains              */
/*  Cla column = positive if Cal A on, negative if cal B on               */
/*      outData  Output OTF object                                        */
/*      infile   Scan part of input file name                             */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  gchar FullFile[128];
  olong irow;
  olong i, j, k, state, port, disk, ver, nrow;
  olong iD, reorder[64], incdatawt, scanNo;
  gchar *tab;
  ofloat *data, ra, dec, cTime, sum, fblank = ObitMagicF();
  ofloat offTime, ldata[64];
  odouble offTimeD, JD;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gboolean valid;
  ObitTableGBTCCBDATA* table;
  ObitTableGBTCCBDATARow* row;
  ObitOTFDesc *desc;
  ObitOTFArrayGeom *geom;
  ObitIOCode retCode;
  gchar *routine = "GetDataOld";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(infile!=NULL);
  g_assert(myInput!=NULL);

  /* Consistency checks */
  desc = outData->myDesc;
  Obit_return_if_fail((nstate==desc->inaxes[desc->jlocstate]), err,
		       "%s: Request beamswitch states incompatible %d !=  %d", 
		      routine, nstate, desc->inaxes[desc->jlocstate]);
  Obit_return_if_fail((nfeed==desc->inaxes[desc->jlocfeed]), err,
		       "%s: Request number feeds incompatible %d !=  %d", 
		      routine, nfeed, desc->inaxes[desc->jlocfeed]);
  Obit_return_if_fail((nstok==desc->inaxes[desc->jlocs]), err,
		       "%s: Request number Poln incompatible %d !=  %d", 
		      routine, nstok, desc->inaxes[desc->jlocs]);
  Obit_return_if_fail((nchan==desc->inaxes[desc->jlocf]), err,
		       "%s: Request number freq. incompatible %d !=  %d", 
		      routine, nchan, desc->inaxes[desc->jlocf]);
  
  /* Get reference MJD, convert ref date  to JD */
  ObitOTFDescDate2JD (outData->geom->RefDate, &JD);
  /* Reference Modified Julian date */
  refMJD = JD - 2400000.5;

  /* get full file name */
  sprintf (FullFile,"CCB26_40/%s.fits", infile);
  
  /* Create table structure */
  table = newObitTableGBTCCBDATA("Data");
  if (err->error) return;

  /* Setup */
  disk = 2;
  tab = "DATA";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(table,disk,FullFile,tab,ver,nrow,err);

  /* Set reordering coefficients to change order in GBT table to OTF */
  SetReorderOld(reorder);

  /* Time offset */
  offTime = 0.0;
  ObitInfoListGetTest(myInput, "offTime",  &type, dim, &offTime);

  /* Correction to time in days */
  offTimeD = offTime / 86400.0;

  /* Get correction from start time to center */
  cTime = 0.5 * integTime;

  /* Scan no. */
  scanNo = 0;
  ObitInfoListGetTest(myInput, "scanNo",  &type, dim, &scanNo);
  if (scanNo!=0) scan = (ofloat)scanNo;
  Obit_log_error(err, OBIT_InfoErr, "Adding scan %s %6.0f to OTF file", 
		 infile, scan);
  ObitErrLog(err);

  /* Open */
  retCode = ObitTableGBTCCBDATAOpen (table, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Create Row structure */
  row = newObitTableGBTCCBDATARow (table);

  desc = outData->myDesc;
  geom = outData->geom;
  incdatawt = desc->incdatawt; /* increment in data-wt axis */
  
  /* Write at end */
  desc->firstRec = desc->nrecord+1;
  startRec = desc->firstRec;  /* First record in scan */
  startTime = -1.0e20;        /* dummy srart time */
  
  /* Loop over table */
  for (irow = 1; irow<=table->myDesc->nrow; irow++) {
    retCode = ObitTableGBTCCBDATAReadRow (table, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, outData->name);
    
    /* first time in scan */
    if (startTime<-1000.0) startTime = row->dmjd + cTime - refMJD;
    data = outData->buffer;     /* Output data array */
    
    /* Check validity */
    if (row->usable) {
    
      /* Fill record entry in data */
      data[desc->iloct]   = row->dmjd + cTime  - refMJD;    /* Time (days) */
      /* Correction for GBT timing error  */
      data[desc->iloct] += offTimeD;
      /* Get position */
      GetPoint (row->dmjd+offTimeD+cTime, &ra, &dec);
      data[desc->ilocti]  = cTime*2.0;    /* time interval (days) */
      data[desc->iloctar] = target;       /* target number */
      data[desc->ilocscan]= scan;         /* scan number */
      data[desc->iloccal] = 0.0;          /* Cal? */
      /* Two possible cal states, sig and ref */
      if (row->cala) data[desc->iloccal] =  1.0;
      if (row->calb) data[desc->iloccal] = -1.0;

      /* Parallactic angle (deg) */
      data[desc->ilocrot] = ObitOTFArrayGeomParAng(geom, data[desc->iloct], ra, dec);
      /* correction to position */
      ObitOTFArrayGeomCorrPoint(azOff, elOff, data[desc->ilocrot], &ra, &dec);
      data[desc->ilocra]  = ra;           /* RA  in deg. */
      data[desc->ilocdec] = dec;          /* Dec in deg.*/
      
      /* Init data to blank */	
      for (i=0; i<ndetect; i++) {
	data[desc->ilocdata+i*incdatawt]   = fblank;
	data[desc->ilocdata+i*incdatawt+1] = 0.0;
      }

      /* Data normalization */
      for (iD=0; iD<table->myDesc->repeat[table->dataCol]; iD++) {
	/* Check validity */
	state = iD / nport;
	port  = iD - state*nport;
	j = PortData->slave[port];
	valid = (!row->ovrflow[iD]) && (row->slaveok[j]);
	if (valid) ldata[iD] = row->data[iD] * dataNorm[reorder[iD]];
	else  ldata[iD] = 0.0;
      }
      
      /* Beamswitching? */
      if (isBS) {
	/* Loop over data */
	for (iD=0; iD<ndetect; iD++) {
	  j = PortData->slave[iD];
	  sum = 0.0;
	  valid =  row->slaveok[j];
	  /* Loop over states */
	  for (k=0; k<StateData->nCCBState; k++) {
	    /* Check validity */
	    valid = valid && (!row->ovrflow[iD+k*nport]);
	    /* Accumulate */
	    if (StateData->phia[k]==StateData->phib[k]) 
	      sum += ldata[iD+k*nport];
	    else sum -= ldata[iD+k*nport];
	  }
	  if (valid) {
	    data[desc->ilocdata+reorder[iD]*incdatawt]   = 0.5 * sum;
	    data[desc->ilocdata+reorder[iD]*incdatawt+1] = 1.0;
	  }
	} /* end detector loop */
      } else { /* Pass all states */
	/* Loop over data */
	for (iD=0; iD<ndetect; iD++) {
	  /* Check validity */
	  state = iD / nport;
	  port  = iD - state*nport;
	  j = PortData->slave[port];
	  valid = (!row->ovrflow[iD]) && (row->slaveok[j]);
	  if (valid) {
	    data[desc->ilocdata+reorder[iD]*incdatawt]   = ldata[iD];
	    data[desc->ilocdata+reorder[iD]*incdatawt+1] = 1.0;
	  }
	}
      } /* end if beamswitch */
      
      /* set number of records */
      desc->numRecBuff = 1;
      
      /* Write output  */
      if ((ObitOTFWrite (outData, NULL, err) != OBIT_IO_OK) || (err->error>0))
	Obit_log_error(err, OBIT_Error, "ERROR writing output Table file");
      if (err->error) Obit_traceback_msg (err, routine, outData->name);
    } /* End valid row */
  } /* end loop over table */
  
  /* Get end times and record numbers */
  desc->firstRec = desc->nrecord+1;
  endTime = row->dmjd + cTime - refMJD;
  endRec  = desc->nrecord;  /* Last record in scan */

  /* Close */
  retCode = ObitTableGBTCCBDATAClose (table, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
  /* Cleanup */
  table = ObitTableGBTCCBDATAUnref(table);
  row = ObitTableGBTCCBDATARowUnref(row);


} /* end GetDataOld  */

void GetPoint (odouble time, ofloat *ra, ofloat *dec)
/*----------------------------------------------------------------------- */
/*  Get antenna pointing at a given time,                                 */
/*  End points or linear interpolation                                    */
/*   Input:                                                               */
/*      time   The desired MJD                                            */
/*   Output:                                                              */
/*      ra     RA J2000 in degrees                                        */
/*      dec    Dec J2000 in degrees                                       */
/*      myInput  parser object                                            */
/*----------------------------------------------------------------------- */
{
  olong i, best;
  odouble test, delta;
  ofloat w1, w2;
  
  /* Find closest */
  best = -1;
  delta = 1.0e20;
  for (i=0; i<nAntTime; i++) { /* loop over array */
    test = fabs (AntDMJD[i]-time);
      if (delta> test) {
	delta = test;
	best = i;
      } else { /* must be getting further, stop */
	break;
      }
  }

  /* end points */
  if ((best==0) || (best==(nAntTime-1))) {
    *ra = AntRA[best];
    *dec = AntDec[best];
  } else if (time<AntDMJD[best]){ /* interpolate with previous */
    w1 = (AntDMJD[best]-time) / (AntDMJD[best]-AntDMJD[best-1]);
    w2 = 1.0 - w1;
    *ra  = w1 * AntRA[best-1]  + w2 * AntRA[best];
    *dec = w1 * AntDec[best-1] + w2 * AntDec[best];
  } else { /* interpolate with following */
    w1 = (AntDMJD[best+1]-time) / (AntDMJD[best+1]-AntDMJD[best]);
    w2 = 1.0 - w1;
    *ra  = w1 * AntRA[best]  + w2 * AntRA[best+1];
    *dec = w1 * AntDec[best] + w2 * AntDec[best+1];
  }
} /* end GetPoint */

olong GetTarget (ObitOTF *outData, gboolean isNewFile, gchar *name, 
		odouble ra, odouble dec, odouble equinox, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get target id, look through existing table create new entry if needed.*/
/*  Returns Target id                                                     */
/*   Input:                                                               */
/*      outData   Output OTF object                                       */
/*      isNewFile True if output file just created                        */
/*      name      Name of target                                          */
/*      ra        RA of target                                            */
/*      dec       Dec of target                                           */
/*      equinox   equinox of ra, dec                                      */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*   Return:                                                              */
/*      Target id                                                         */
/*----------------------------------------------------------------------- */
{
  olong targ = -1;
  ObitTableOTFTarget* table;
  ObitTableOTFTargetRow* row;
  olong iRow, ver;
  gboolean doWrite;
  ObitIOAccess access;
  gchar tName[20];
  gchar *routine = "GetTarget";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return targ;
  g_assert (ObitOTFIsA(outData));
  g_assert(name!=NULL);

  /* create Target table object */
  ver = 1;
  if (isNewFile) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  table = newObitTableOTFTargetValue ("Target table", (ObitData*)outData, &ver, access, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, targ);

  /* Open table */
  if ((ObitTableOTFTargetOpen (table, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFTarget table");
    return targ;
  }

  /* Create Row */
  row = newObitTableOTFTargetRow (table);

  /* attach to table buffer */
  ObitTableOTFTargetSetRow (table, row, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, targ);

  /* Newly created?  Just write new one */
  doWrite = FALSE;
  if (isNewFile) {
    targ = 1;
    row->TargID = targ;
    strncpy(row->Target, name, 16);
    row->RAMean  = ra;
    row->DecMean = dec;
    row->Epoch   = equinox;
    doWrite = TRUE;
  } else { /* Existing, see if already exists? */

    /* loop through table */
    for (iRow = 1; iRow<=table->myDesc->nrow; iRow++) {
      if ((ObitTableOTFTargetReadRow (table, iRow, row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading OTFTarget Table file");
	return targ;
      }
      strncpy (tName, row->Target, 16); tName[16]=0;
      ObitTrimTrail(tName);
      if (!strncmp (tName, name, 16)) {
	/* Found match */
	targ = row->TargID;
	break;
      }  
    } /* end loop over table */

    /* Add new entry? */
    if (targ<=0) {
      targ = table->myDesc->nrow + 1;
      row->TargID = targ;
      row->RAMean  = ra;
      row->DecMean = dec;
      row->Epoch   = equinox;
      strncpy(row->Target, name, 16);
      doWrite = TRUE;
    }
  } /* end output table already exists */

  /* need to write new entry? */
  if (doWrite) {
    iRow = table->myDesc->nrow + 1;
    if ((ObitTableOTFTargetWriteRow (table, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR writing OTFTarget Table file");
      return targ;
    }
  }
  
 /* Close  table */
  if ((ObitTableOTFTargetClose (table, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFTarget Table file");
    return targ;
  }

  /* Cleanup */
  row = ObitTableOTFTargetRowUnref(row);
  table = ObitTableOTFTargetUnref(table);

  return targ;
} /* end  GetTarget */

void InitScan (ObitOTF *outData, gboolean isNewFile, ofloat scan, 
	       ofloat target, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Initializes Index table for this scan creating an entry               */
/*   Input:                                                               */
/*      outData   Output OTF object                                       */
/*      isNewFile True if output file just created                        */
/*      scan      Scan number                                             */
/*      target    Target ID number                                        */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitTableOTFIndex* table;
  ObitTableOTFIndexRow* row;
  olong iRow, ver;
  olong scanID, targetID;
  ObitIOAccess access;
  gchar *routine = "InitScan";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(outData));

  /* create Index table object */
  scanID = (olong)(scan+0.5);
  targetID = (olong)(target+0.5);
  ver = 1;
  if (isNewFile) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  table = newObitTableOTFIndexValue ("Index table", (ObitData*)outData, &ver, access, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open table */
  if ((ObitTableOTFIndexOpen (table, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFIndex table");
    return;
  }

  /* Create Row */
  row = newObitTableOTFIndexRow (table);

  /* initialize row */
  row->ScanID = scanID;
  row->TargetID = targetID;
  row->Time = 0.0;
  row->TimeI = 0.0;
  row->StartRec = -1;
  row->EndRec = -1;

  /* attach to table buffer */
  ObitTableOTFIndexSetRow (table, row, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Write at end of table */
  iRow = table->myDesc->nrow + 1;
  if ((ObitTableOTFIndexWriteRow (table, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR writing OTFIndex Table file");
    return;
  }

 /* Close  table */
  if ((ObitTableOTFIndexClose (table, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFIndex Table file");
    return;
  }

  /* Cleanup */
  row = ObitTableOTFIndexRowUnref(row);
  table = ObitTableOTFIndexUnref(table);

} /* end  InitScan */

void SetScan (ObitOTF *outData, odouble startTime, odouble endTime, 
	      olong startRec, olong endRec, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Initializes Index table for this scan creating an entry               */
/*   Input:                                                               */
/*      outData   Output OTF object                                       */
/*      startTime Start time of scan in days                              */
/*      endTime   End time of scan in days                                */
/*      startRec  First record in scan                                    */
/*      endRec    Last record in scan                                     */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitTableOTFIndex* table;
  ObitTableOTFIndexRow* row;
  olong iRow, ver;
  ObitIOAccess access;
  gchar *routine = "SetScan";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(outData));

  /* create Index table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  table = newObitTableOTFIndexValue ("Index table", (ObitData*)outData, &ver, access, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open table */
  if ((ObitTableOTFIndexOpen (table, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFIndex table");
    return;
  }

  /* Create Row */
  row = newObitTableOTFIndexRow (table);

  /* attach to table buffer */
  ObitTableOTFIndexSetRow (table, row, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Update last record */
  iRow = table->myDesc->nrow;
  if ((ObitTableOTFIndexReadRow (table, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR reading OTFIndex Table file");
    return;
  }

  /* upate row */
  row->Time = 0.5 * (startTime + endTime);
  row->TimeI = (endTime - startTime);
  row->StartRec = startRec;
  row->EndRec   = endRec;

  /* Rewrite at end of table */
  iRow = table->myDesc->nrow;
  if ((ObitTableOTFIndexWriteRow (table, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR writing OTFIndex Table file");
    return;
  }

 /* Close  table */
  if ((ObitTableOTFIndexClose (table, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFIndex Table file");
    return;
  }

  /* Cleanup */
  row = ObitTableOTFIndexRowUnref(row);
  table = ObitTableOTFIndexUnref(table);

} /* end  SetScan */

void SetReorderNew (olong *reorder)
/*----------------------------------------------------------------------- */
/*  Fills index array for reordering of data from GBT FITS to OTF         */
/*   Output order:                                                        */
/*     Stokes axis R, L                                                   */
/*     Feed axis 1, 2                                                     */
/*     Freq axis  2.775E+10, 3.125E+10, 3.475E+10,3.825E+10               */ 
/*     State axis (0,0), (0,1), (1,0), (1,1)                              */
/*   Output:                                                              */
/*       reorder  Reordering array,value is output 0-rel index            */
/*----------------------------------------------------------------------- */
{
  olong i, j, istate, ipol, ifeed, ifreq;
  /* Set reordering coefficients to change order in GBT table to OTF */

  istate = 0;
  /* Loop over state */
  for (j=0; j<StateData->nCCBState; j++) {
    if ((StateData->phia[j]==0) && (StateData->phib[j]==0)) istate = 0;
    else if ((StateData->phia[j]==0) && (StateData->phib[j]==1)) istate = 1;
    else if ((StateData->phia[j]==1) && (StateData->phib[j]==0)) istate = 2;
    else istate = 3;
    /* Loop over port */
    for (i=0; i<nport; i++) {
      if (IFdata->poln[i]=='R') ipol = 0;
      else ipol = 1;
      ifeed = IFdata->feed[i]-1;
      if ((IFdata->refFrequency[i]>=2.5e10) && (IFdata->refFrequency[i]<=3.0e10)) ifreq = 0;
      else if ((IFdata->refFrequency[i]>=3.0e10) && (IFdata->refFrequency[i]<=3.3e10)) ifreq = 1;
      else if ((IFdata->refFrequency[i]>=3.3e10) && (IFdata->refFrequency[i]<=3.6e10)) ifreq = 2;
      else ifreq = 3;
      /* State dependence 0 & 4 are sig, 1 & 2 are ref */
      if ((istate==1) || (istate==2)) { /* switch feed  and poln */
	ifeed = 1-ifeed;
	ipol  = 1-ipol;
      }
      reorder[i+j*nport] = ifeed*nstok + ifreq*nstok*nfeed + istate*nport;
    }  /* End loop over port */
  }  /* End loop over state */

} /* end  SetReorderNew */

void SetReorderOld (olong *reorder)
/*----------------------------------------------------------------------- */
/*  Fills index array for reordering of data from GBT FITS to OTF         */
/*   Output order:                                                        */
/*     Stokes axis R, L                                                   */
/*     Feed axis 1, 2                                                     */
/*     Freq axis  2.775E+10, 3.125E+10, 3.475E+10,3.825E+10               */ 
/*     State axis (0,0), (0,1), (1,0), (1,1)                              */
/*   Output:                                                              */
/*       reorder  Reordering array,value is output 0-rel index            */
/*----------------------------------------------------------------------- */
{
  olong i, j, istate, ipol, ifeed, ifreq;
  /* Set reordering coefficients to change order in GBT table to OTF */

  istate = 0;
  /* Loop over state */
  for (j=0; j<StateData->nCCBState; j++) {
    if ((StateData->phia[j]==0) && (StateData->phib[j]==0)) istate = 0;
    else if ((StateData->phia[j]==0) && (StateData->phib[j]==1)) istate = 1;
    else if ((StateData->phia[j]==1) && (StateData->phib[j]==0)) istate = 2;
    else istate = 3;
    /* Loop over port */
    for (i=0; i<nport; i++) {
      if (IFdata->poln[i]=='R') ipol = 0;
      else ipol = 1;
      ifeed = IFdata->feed[i]-1;
      if ((IFdata->refFrequency[i]>=2.5e10) && (IFdata->refFrequency[i]<=3.0e10)) ifreq = 0;
      else if ((IFdata->refFrequency[i]>=3.0e10) && (IFdata->refFrequency[i]<=3.3e10)) ifreq = 1;
      else if ((IFdata->refFrequency[i]>=3.3e10) && (IFdata->refFrequency[i]<=3.6e10)) ifreq = 2;
      else ifreq = 3;
      /* State dependence 0 & 4 are sig, 1 & 2 are ref */
      if ((istate==1) || (istate==2)) { /* switch feed  and poln */
	ifeed = 1-ifeed;
	ipol  = 1-ipol;
      }
      reorder[i+j*nport] = ipol + ifeed*nstok + ifreq*nstok*nfeed + istate*nport;
    }  /* End loop over port */
  }  /* End loop over state */

} /* end  SetReorderOld */

