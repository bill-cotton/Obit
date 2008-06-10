/* $Id: PASim2OTF.c,v 1.6 2006/09/29 16:31:04 bcotton Exp $                            */
/*--------------------------------------------------------------------*/
/* Convert Penn Array IDL simulator data to OTF format                */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008-2008                                     */
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
#include "ObitOTFCal.h"
#include "ObitOTFUtil.h"
#include "ObitIOOTFFITS.h"
#include "ObitFITS.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitTableOTFTarget.h"
#include "ObitTableOTFIndex.h"
#include <errno.h>

/* internal prototypes */
/* Get inputs */
ObitInfoList* PASim2OTFin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Get file descriptor */
void GetHeader (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, 
		ObitErr *err);
/* Get next record */
gboolean GetRecord (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, 
		    ObitErr *err);
/* Read pointing file */
void GetAntenna (gchar *infile, ObitInfoList *myInput, ObitErr *err);
/* Get antenna ponting for a time */
void GetPoint (gint integ, odouble *ra, odouble *dec, odouble *mjd);
/* Get target id */
olong GetTarget (ObitOTF *outData, gboolean isNew, gchar *name, ObitErr *err);
/* Initializes Index table for this scan creating an entry   */
void InitScan (ObitOTF *outData, gboolean isNew, ofloat scan, ofloat target,
	       ObitErr *err);
/* Undate scan info */
void SetScan (ObitOTF *outData, odouble startTime, odouble endTime, 
	      olong startRec, olong endRec, ObitErr *err);

/* Program globals */
/* input geometry file */
gchar *geom_file = NULL;
/* input data file */
gchar *data_file = NULL;
/* input file pointer */
FILE *myFile = NULL;
/* Input file line buffer */
gchar myLine[4196];
/* Cal equivalent in Jy */
odouble calValue;
/* Reference date as yyyy-mm-dd */
gchar refDate[20];
/* Reference mjd */
odouble refMJD;
/* Time offset (days) to add to label time */
odouble TimeOff;
/* Last time (days) */
odouble lastTime=100000.0;
/* Time offset to add */
ofloat scan;      /* scan number */
olong  startRec;  /* First record number (1-rel) in scan */
olong  endRec;    /* End record number (1-rel) in scan */
odouble timeOffset = 0.0;
odouble startTime;/* Start time of Scan in days */
odouble endTime;  /* End time of Scan in days */
olong ncount = 0;  /* count of how many left to blank after a cal. */
olong numInteg=0;  /* number of integrations read */
olong nAntTime;   /* number of antenna time samples */
odouble *AntDMJD; /* Array of antenna times */
odouble *AntRA;   /* Array of Antenna RA J2000 values */
odouble *AntDec;  /* Array of Antenna RA J2000 values */
gchar Name[48];   /* Target/scan type name */
ofloat target;    /* target number */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Convert to Obit/OTF FITS tables format                               */
/*----------------------------------------------------------------------- */
{
  olong  disk, nrec, iscan, i, ierr=0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *outData= NULL;
  ObitErr *err= NULL;
  ObitIOAccess access;
  gchar *FITSdir[] = {"FITSdata/"}, *fullname=NULL;
  gchar geomfile[128], datafile[128], outfile[128];
  ObitInfoType type;
  gboolean isNew;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *errMsg, FullFile[128];
  gboolean done;

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("PASim2OTF", 1, 0, 0, NULL, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  ierr = 0;
  myInput = PASim2OTFin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input file names */
  for (i=0; i<128; i++) geomfile[i] = 0;
  ObitInfoListGet(myInput, "geom", &type, dim, geomfile, err);
  geomfile[dim[0]] = 0;  /* null terminate */

  for (i=0; i<128; i++) datafile[i] = 0;
  ObitInfoListGet(myInput, "data", &type, dim, datafile, err);
  datafile[dim[0]] = 0;  /* null terminate */

  /* output FITS file name */
  for (i=0; i<128; i++) outfile[i] = 0;
  ObitInfoListGet(myInput, "outfile", &type, dim, outfile, err);
  outfile[dim[0]] = 0;  /* null terminate */

  /* Target name/scantype */
  ObitInfoListGet(myInput, "Target", &type, dim, &Name, err);
  Name[dim[0]] = 0;  /* null terminate */

  /* Scan time offset in days */
  ObitInfoListGet(myInput, "TimeOff", &type, dim, &TimeOff, err);

  /* Scan number */
  ObitInfoListGet(myInput, "ScanNo", &type, dim, &iscan, err);
  scan = (ofloat)iscan;

  /* Create ObitOTF for data */
  outData = newObitOTF("Output data");

  /* Get header info, array geometry */
  GetHeader (outData, geomfile, myInput, err);
  
  /* Get Telescope position array */
  GetAntenna (datafile, myInput, err);

  /* Define output, I/O size */
  disk = 1;
  nrec = 1;
  ObitOTFSetFITS(outData,nrec,disk,outfile,err);
  
  /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;

  /* See if out exists, if so append */
  fullname = ObitFITSFilename (disk, outfile, err);
  isNew = !ObitFileExist (fullname, err);
  if (err->error>0) Obit_log_error(err, OBIT_Error, "ERROR testing file %s", fullname);
 
  /* Open output OTF */
  if (isNew) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  if ((ObitOTFOpen (outData, access, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", outfile);

  /* Write at end */
  outData->myDesc->firstRec = outData->myDesc->nrecord+1;
  startRec = outData->myDesc->firstRec;  /* First record in scan */
  startTime = -1.0e20;                   /* dummy srart time */
   
  /* Get target id number */
  target = (ofloat)GetTarget (outData, isNew, Name, err);

  /* Initialize scan in Index table */
  InitScan (outData, isNew, scan, target, err);

  /* Open input data text file */
  /* get full file name */
  sprintf (FullFile,"%s.raw", datafile);
  myFile = fopen (FullFile, "rt");
  if (myFile==NULL) {
    Obit_log_error(err, OBIT_Error, "ERROR opening file %s", FullFile);
    errMsg = strerror(errno);
    Obit_log_error(err, OBIT_Error, "%s", errMsg);
  }
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
  
  /* Loop translating file */
  done = FALSE;
  while (!done) {
    done = GetRecord (outData, datafile, myInput, err);
    if (err->error) done = TRUE;
    /* Write record */
    if (!done) {
     if ((ObitOTFWrite (outData, NULL, err) != OBIT_IO_OK) || (err->error>0))
      Obit_log_error(err, OBIT_Error, "ERROR writing output Table file");
    }
    if (err->error) done = TRUE;
  } /* end loop translating */
  
   /* show any errors */
   if (err->error) ierr = 1;
   ObitErrLog(err);
   if (ierr!=0) return ierr;

   /* record number */
   endRec  = outData->myDesc->nrecord;
  
  /* Update index table */
  SetScan (outData, startTime, endTime, startRec, endRec, err);

  /* Close */
  if ((ObitOTFClose (outData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing output file");
  
  /* show any errors */
   if (err->error) ierr = 1;
   ObitErrLog(err);
   if (ierr!=0) return ierr;
   
   /* Close input file */
   fclose(myFile);
  
   /* Shutdown Obit */
   mySystem = ObitSystemShutdown (mySystem);
   
   /* cleanup */
   myInput = ObitInfoListUnref(myInput);  /* delete input list */
   outData = ObitUnref(outData);
   if (fullname) g_free(fullname);
 
   return ierr;
} /* end of main */

ObitInfoList* PASim2OTFin (int argc, char **argv, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Parse control info from command line                                  */
/*   Input:                                                               */
/*      argc   Number of arguments from command line                      */
/*      argv   Array of strings from command line                         */
/*   Output:                                                              */
/*      geom_file   input ascii array geometry file                       */
/*      data_file   input ascii data file                                 */
/*      calValue    Cal equivalent in Jy                                  */
/*      refDate     Reference date as yyyy-mm-dd                          */
/*      TimeOff     Time offset (days) to add to label time               */
/*      scan        Current scan number                                   */
/*      err    Obit Error stack                                           */
/*   return  parser list                                                  */
/*----------------------------------------------------------------------- */
{
  olong ax;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *input_file="PASim2OTF.in", *arg;
  gchar *strTemp;
  oint    itemp;
  odouble dtemp;
  ObitInfoList* list;
  gboolean init=FALSE;

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

    } else if (strcmp(arg, "-geom") == 0){ /* array geometry */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "geom", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-data") == 0){ /* data */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "data", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-outfile") == 0){ /* output FITS file name */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outfile", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-refDate") == 0){ /* reference date */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "data", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-Target") == 0){ /* array geometry */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "Target", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-ScanNo") == 0) { /* Scan number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "ScanNo", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-TimeOff") == 0) { /* Scan time offset in days */
      dim[0] = 1;
      dtemp = strtod(argv[++ax], NULL);
      ObitInfoListPut (list, "TimeOff", OBIT_double, dim, &dtemp, err);

    } else if (strcmp(arg, "-calValue") == 0) { /* Jy equivalent of cal */
      dim[0] = 1;
      dtemp = strtod(argv[++ax], NULL);
      ObitInfoListPut (list, "calValue", OBIT_double, dim, &dtemp, err);

    } else { /* unknown argument */
      Usage();
    }
  }
  
  /* parse input file if specified */
  if (!init) ObitParserParse (input_file, list, err);

  return list;
} /* end PASim2OTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: PASim2OTF -input file [-geom file -data file...\n");
    fprintf(stderr, "Convert an external file format to Obit/OTF\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def PASim2OTF.in\n");
    fprintf(stderr, "  -geom array geometry file, def PennArray.geom\n");
    fprintf(stderr, "  -data input data file, def PennArray.dat\n");
    fprintf(stderr, "  -outfile output FITS file name, def PASimOTF.fits\n");
    fprintf(stderr, "  -refDate reference date,  def 20030500\n");
    fprintf(stderr, "  -ScanNo Scan number (better set) \n");
    fprintf(stderr, "  -TimeOff Scan time offset in days (better set) \n");
    fprintf(stderr, "  -calValue Jy equivalent of cal, def 1.0\n");
    fprintf(stderr, "  -Target Target/scan type name def Unknown\n");
    
    /*/exit(1);  bail out */
  }/* end Usage */

/*----------------------------------------------------------------------- */
/*  Create default input ObitInfoList                                     */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  odouble dtemp;
  oint itemp;
  ObitInfoList *out = newObitInfoList();

  /* add parser items */
  /* output FITS file name */
  strTemp = "PASimOTF.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "outfile", OBIT_string, dim, strTemp, err);

  /* input geometry file name */
  strTemp = "PennArray.geom";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "geom", OBIT_string, dim, strTemp, err);

  /* input data file name */
  strTemp = "PennArray.dat";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "data", OBIT_string, dim, strTemp, err);

  /* input Target name */
  strTemp = "Unknown";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "Target", OBIT_string, dim, strTemp, err);

  /* reference date */
  strTemp = "20030500";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "refDate", OBIT_string, dim, strTemp, err);

  /* Jy equivalent of cal */
  dim[0] = 1;
  dtemp = 1.0;
  ObitInfoListPut (out, "calValue", OBIT_double, dim, &dtemp, err);

  /* Scan time offset in days */
  dim[0] = 1;
  dtemp = 0.0;
  ObitInfoListPut (out, "TimeOff", OBIT_double, dim, &dtemp, err);

  /* Scan number */
  dim[0] = 1;
  itemp = 1;
   ObitInfoListPut (out, "ScanNo", OBIT_oint, dim, &itemp, err);
 

  return out;
} /* end defaultInputs */

void GetHeader (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, 
		ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get header information from text file/parse list                      */
/*   Input:                                                               */
/*      outData  Output OTF object                                        */
/*      infile   input geometry file name                                 */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitOTFDesc *desc;
  ObitOTFArrayGeom *geom;
  olong numberDetect, ncol, i, ncopy;
  gchar *errMsg;
  odouble JD, T, GMST0, GSTIAT0;
  ofloat azOff, elOff;
  olong nscan;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(outData));
  g_assert(infile!=NULL);
  g_assert(myInput!=NULL);

  /* get values from myInput */
  /* reference date */
  for (i=0; i<20; i++) refDate[i] = 0;
  ObitInfoListGet(myInput, "refDate", &type, dim, refDate, err);
  refDate[dim[0]] = 0;  /* null terminate */

  /* Convert to JD */
  ObitOTFDescDate2JD (refDate, &JD);
  /* Modified Julian date */
  refMJD = JD - 2400000.5;

  /* GST at IAT=0 at 0h on reference date (deg)*/
  /* Tropical century from jan 0.5, 1900 */
  T = (JD - 2451545.0) / 36525.0;

  /* GMST at IAT=0 in radians */
  GMST0 = ((((((-6.2e-6 * T) + 0.093104) * T) + 8640184.812866) * T + 24110.54841) 
	   * 2.0 * G_PI / 86400.0);
  /* to degrees */
  GSTIAT0 = RAD2DG * fmod(GMST0, (2.0*G_PI));

  /* Jy equivalent of cal */
  ObitInfoListGet(myInput, "calValue", &type, dim, &calValue, err);

  desc = outData->myDesc;

  /* Set Feed Array Geometry */
  /* create Array geometry with 64 elements */
  numberDetect = 64;
  geom = ObitOTFArrayGeomCreate (numberDetect);
  geom->azOffset = g_realloc(geom->azOffset, numberDetect*sizeof(ofloat));
  geom->elOffset = g_realloc(geom->elOffset, numberDetect*sizeof(ofloat));
 
  /* Other information - time info mostly for 6 May 2003 */
  ncopy = strlen (refDate);
  for (i=0; i<ncopy; i++) geom->RefDate[i] = refDate[i]; geom->RefDate[i]=0;
  geom->TimeSys[0] = 'U'; geom->TimeSys[1] = 'T';geom->TimeSys[2] = 'C';geom->TimeSys[3] = 0;
  geom->TeleX   =  882879.8949; /* 140 ft */
  geom->TeleY   = -4924482.3088;
  geom->TeleZ   =  3944130.6875;
  geom->DegDay  = 3.6098564497330e+02;
  geom->GSTiat0 = GSTIAT0;
  geom->PolarX  = -1.4314E-05;
  geom->PolarY  = 1.4333E-04;
  geom->ut1Utc  = -3.6465E-01;
  geom->dataUtc = 0.0;
  geom->iatUtc  = 0.0;

  /* Compute some useful terms */
  /* telescope latitude in radians */
  if (fabs(geom->TeleX)<1.0) geom->TeleX = 1.0;
  geom->lat = asin (geom->TeleZ / 
		    sqrt(geom->TeleX*geom->TeleX + 
			 geom->TeleY*geom->TeleY + 
			 geom->TeleZ*geom->TeleZ));
  /* telescope longitude in radians */
  geom->lon = atan2 (geom->TeleY,  geom->TeleX);
  /* LST at iat0 in radians */
  geom->LSTiat0 = geom->GSTiat0*1.74533e-2 + geom->lon;
  /* Earth rotation rate in rad/day */
  geom->RadDay = geom->DegDay*1.74533e-2;
  /* Data - IAT in days */
  geom->dataIat = (geom->dataUtc - geom->iatUtc) / 86400.0;

  /* Attach Array geometry to OTF */
  outData->geom = ObitOTFArrayGeomRef(geom);
  
  /* read geometry from table */ 
  /* Open input text file */
  myFile = fopen (infile, "r");
  if (myFile==NULL) {
    Obit_log_error(err, OBIT_Error, "ERROR opening file %s", infile);
    errMsg = strerror(errno);
    Obit_log_error(err, OBIT_Error, "%s", errMsg);
    return;
  }
  /* Loop over file */
  for (i=0; i<numberDetect; i++) {
    /* read data */
    fgets(myLine, 4195, myFile);
    if (ferror(myFile)) {
      Obit_log_error(err, OBIT_Error, "ERROR opening file %s", infile);
      errMsg = strerror(errno);
      Obit_log_error(err, OBIT_Error, "%s", errMsg);
      return;
    }
    
    /* Parse text line */
    nscan = sscanf (myLine,"%f %f", &azOff, &elOff);
    if (nscan<2) {
      Obit_log_error(err, OBIT_Error, "ERROR parsing line %s", myLine);
      errMsg = strerror(errno);
      Obit_log_error(err, OBIT_Error, "%s", errMsg);
      return;
    }
    geom->azOffset[i] = azOff;
    geom->elOffset[i] = elOff;
    /* DEBUG */
    /*geom->azOffset[i] = -azOff;*/
    /*geom->elOffset[i] = -elOff;*/
  } /* end loop over file */

  /* Close input file */
  fclose(myFile);
  
  /* delete object */
  geom = ObitOTFArrayGeomUnref(geom);

  /* Initialize  Descriptor */
  /* Beam size (GBT at 90 GHz) */
  desc->beamSize = 2.0e8 / 90.0e9;

  strncpy (desc->object, "Sky", OTFLEN_VALUE);
  strncpy (desc->teles,  "GBT       ", OTFLEN_VALUE);
  strncpy (desc->origin, "Obit ", OTFLEN_VALUE);
  desc->isort[0] = 'T';  /* Time ordered */
  ncol = 0;

  desc->JDObs = 0.0;
  desc->epoch = 2000.0;
  desc->equinox = 2000.0;
  strncpy (desc->bunit,  "ADU     ", OTFLEN_VALUE);
  strncpy (desc->obsdat, refDate, OTFLEN_VALUE);

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
  desc->colRepeat[ncol] = numberDetect;
  ncol++;

  desc->ncol = ncol;

  /* Data array descriptors */
  desc->naxis = 0;
  
  /* Stokes axis (Stokes I) */
  desc->inaxes[desc->naxis] = 1;
  strncpy (desc->ctype[desc->naxis], "STOKES  ", OTFLEN_KEYWORD);
  desc->crpix[desc->naxis] = 1.0;
  desc->crota[desc->naxis] = 0.0;
  desc->cdelt[desc->naxis] = 1.0;
  desc->crval[desc->naxis] = 1.0; /* Only "I" */
  desc->naxis++;
  
  /* Feed axis */
  desc->inaxes[desc->naxis] = numberDetect;
  strncpy (desc->ctype[desc->naxis], "FEED", OTFLEN_KEYWORD);
  desc->cdelt[desc->naxis] = 1.0;
  desc->crpix[desc->naxis] = 1.0;
  desc->crota[desc->naxis] = 0.0;
  desc->crval[desc->naxis] = 1.0;
  desc->naxis++;
  
  /* Frequency axis (90 GHz) */
  desc->inaxes[desc->naxis] = 1;
  strncpy (desc->ctype[desc->naxis], "FREQ    ", OTFLEN_KEYWORD);
  desc->cdelt[desc->naxis] = 1.0;
  desc->crpix[desc->naxis] = 1.0;
  desc->crota[desc->naxis] = 0.0;
  desc->crval[desc->naxis] = 90.0e9;
  desc->naxis++;

  /* Index the descriptor */
  ObitOTFDescIndex (desc);

} /* end GetHeader */

void GetAntenna (gchar *infile, ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get antenna positions and leave in globals                            */
/*      infile   root of input file names                                 */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  gchar FullFile[128], *errMsg;
  FILE *thisFile = NULL;
  olong nscan, ierr = 0;
  gboolean ImDone = FALSE;
  olong count = 0, irow;
  odouble dmjd, az, el, raj2000, decj2000;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(infile!=NULL);
  g_assert(myInput!=NULL);

  /* get full file name */
  sprintf (FullFile,"%s.att", infile);
  
  /* Pass throuth the file to count entries */
  /* Open input data text file */
  thisFile = fopen (FullFile, "rt");
  if (thisFile==NULL) {
    Obit_log_error(err, OBIT_Error, "ERROR opening file %s", FullFile);
    errMsg = strerror(errno);
    Obit_log_error(err, OBIT_Error, "%s", errMsg);
  }
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return;

  /* Loop over file */
  while (!ImDone) {
  
    fgets(myLine, 4195, thisFile);
    if (ferror(thisFile)) {
      Obit_log_error(err, OBIT_Error, "ERROR reading file %s", infile);
      errMsg = strerror(errno);
      Obit_log_error(err, OBIT_Error, "%s", errMsg);
      return;
    }
    ImDone = feof(thisFile); /* hit EOF? */
    if (ImDone) break;

    count++;   /* count 'em */
  } /* end loop over file */
  /* Close input file */
  fclose(thisFile);

  /* make sure there is data */
  if (count<=0) {
     Obit_log_error(err, OBIT_Error, "No data in Antenna file for scan %s", infile);
     return;
 }

  /* Allocate arrays */
  nAntTime = count;
  AntDMJD = g_malloc0(count*sizeof(odouble));
  AntRA   = g_malloc0(count*sizeof(odouble));
  AntDec  = g_malloc0(count*sizeof(odouble));

  /* Another pass to read the data */
  /* Open input data text file */
  thisFile = fopen (FullFile, "rt");
  if (thisFile==NULL) {
    Obit_log_error(err, OBIT_Error, "ERROR opening file %s", FullFile);
    errMsg = strerror(errno);
    Obit_log_error(err, OBIT_Error, "%s", errMsg);
  }
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return;

  /* Loop over file */
  irow = 0; ImDone = FALSE;
  while (!ImDone) {
  
    fgets(myLine, 4195, thisFile);
    if (ferror(thisFile)) {
      Obit_log_error(err, OBIT_Error, "ERROR reading file %s", infile);
      errMsg = strerror(errno);
      Obit_log_error(err, OBIT_Error, "%s", errMsg);
      return;
    }
    ImDone = feof(thisFile); /* hit EOF? */
    if (ImDone) break;

    /* parse data */
    nscan = sscanf (myLine,"%lf %lf %lf %lf %lf", &dmjd, &az, &el, &raj2000, &decj2000);
    if (nscan<5) {
      Obit_log_error(err, OBIT_Error, "ERROR parsing line %s", myLine);
      errMsg = strerror(errno);
      Obit_log_error(err, OBIT_Error, "%s", errMsg);
      return;
    }
    irow++;
    AntDMJD[irow-1] = dmjd;
    AntRA[irow-1]   = raj2000;
    AntDec[irow-1]  = decj2000;
  } /* end loop over file */
  /* Close input file */
  fclose(thisFile);

} /* end GetAntenna  */

gboolean GetRecord (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, 
		    ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get next record from input and copy to outData buffer                 */
/*  Uses global file IO pointer myFile                                    */
/*  Fills in two records, one call off, one cal on                        */
/*   Input:                                                               */
/*      outData  Output OTF object                                        */
/*      infile   input text file                                          */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*   return  parser list                                                  */
/*----------------------------------------------------------------------- */
{
  gboolean done=TRUE;
  ObitOTFDesc *desc;
  ObitOTFArrayGeom *geom;
  gchar *cdata, *errMsg;
  olong i, nskip=3;
  ofloat *data, ra, dec;
  odouble mjd, t[64], ira, idec, ical, itime, irot;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitOTFIsA(outData));

  desc = outData->myDesc;
  geom = outData->geom;

  /* read data */
  fgets(myLine, 4195, myFile);
  if (ferror(myFile)) {
    Obit_log_error(err, OBIT_Error, "ERROR reading file %s", infile);
    errMsg = strerror(errno);
    Obit_log_error(err, OBIT_Error, "%s", errMsg);
    return TRUE;
  }
  done = feof(myFile); /* hit EOF? */
  if (done) return done;

  /* Parse text line */
  /* Data */
  cdata = myLine;
  for (i=0; i<geom->numberDetect; i++) t[i] = strtod (cdata, &cdata);
  ical = strtod (cdata, &cdata);
  itime = strtod (cdata, &cdata);
  itime /= 86400.0;  /* time to days */
  irot = 0.0;   /* No info */

  /* get pointing position */
  GetPoint (numInteg++, &ira, &idec, &mjd);
 
  lastTime = mjd;
  
  /* Fill record entry in data */
  /* have to fiddle times to get them to increase */
  ra  = (ofloat)ira;
  if (ra<0.0)   ra += 360.0;
  if (ra>360.0) ra -= 360.0;
  dec = (ofloat)idec;
  if (ical>0.0)  {
    ncount = nskip;
    ical = calValue;
  }

  data = outData->buffer;
  data[desc->iloct]   = mjd-refMJD+timeOffset;/* Time (days) */
  data[desc->ilocti]  = 0.0;           /* time interval (days) */
  data[desc->iloctar] = target;        /* target number */
  data[desc->ilocra]  = ra;            /* RA  in deg. */
  data[desc->ilocdec] = dec;           /* Dec in deg.*/
  data[desc->iloccal] = (ofloat)ical;  /* cal? */
  data[desc->ilocrot] = (ofloat)irot;  /* Parallactic angle (deg) */
  /*data[desc->ilocrot] = (ofloat)(irot+90.0); debug */
  data[desc->ilocscan]= scan;  /* scan number */
  
  /* set data values */
  for (i=0; i<geom->numberDetect; i++) data[desc->ilocdata+i] = (ofloat)t[i];

  /* set number of records */
  desc->numRecBuff = 1;

  /* ignore some after cal */
     if ((ncount>0) & (ncount<nskip)) desc->numRecBuff = 0;
  ncount--;

  /* Append to end of file */
  desc->firstRec = desc->nrecord + 1;

  /* first time in scan */
  if (startTime<-1000.0) startTime = mjd - refMJD;

  /* Get end time */
  endTime = mjd - refMJD;
 
  return done;
} /* end GetRecord */

void GetPoint (gint integ, odouble *ra, odouble *dec, odouble *mjd)
/*----------------------------------------------------------------------- */
/*  Get antenna pointing for a given integration                          */
/*   Input:                                                               */
/*      integ  Integration number 0-rel                                   */
/*   Output:                                                              */
/*      ra     RA J2000 in degrees                                        */
/*      dec    Dec J2000 in degrees                                       */
/*      mjd    MJD                                                        */
/*----------------------------------------------------------------------- */
{

  *ra  = AntRA[integ];
  *dec = AntDec[integ];
  *mjd = AntDMJD[integ];
} /* end GetPoint */

olong GetTarget (ObitOTF *outData, gboolean isNew, gchar *name, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get target id, look through existing table create new entry if needed.*/
/*  Returns Target id                                                     */
/*   Input:                                                               */
/*      outData  Output OTF object                                        */
/*      isNew    True if output file just created                         */
/*      name     Name of target                                           */
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
  gchar *routine = "GetTarget";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return targ;
  g_assert (ObitOTFIsA(outData));
  g_assert(name!=NULL);

  /* create Target table object */
  ver = 1;
  if (isNew) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  table = newObitTableOTFTargetValue ("Target table", (Obit*)outData, &ver, access, err);
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
  if (isNew) {
    targ = 1;
    row->TargID = targ;
    strncpy(row->Target, name, 16);
    doWrite = TRUE;
  } else { /* Existing, see if already exists? */

    /* loop through table */
    for (iRow = 1; iRow<=table->myDesc->nrow; iRow++) {
      if ((ObitTableOTFTargetReadRow (table, iRow, row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading OTFTarget Table file");
	return targ;
      }
      if (!strncmp (row->Target, name, 16)) {
	/* Found match */
	targ = row->TargID;
	break;
      }  
    } /* end loop over table */

    /* Add new entry? */
    if (targ<=0) {
      targ = table->myDesc->nrow + 1;
      row->TargID = targ;
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

void InitScan (ObitOTF *outData, gboolean isNew, ofloat scan, ofloat target,
	       ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Initializes Index table for this scan creating an entry               */
/*   Input:                                                               */
/*      outData  Output OTF object                                        */
/*      isNew    True if output file just created                         */
/*      scan     Scan number                                              */
/*      target   Target ID number                                         */
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
  if (isNew) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  table = newObitTableOTFIndexValue ("Index table", (Obit*)outData, &ver, access, err);
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
/*      isNew     True if output file just created                        */
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
  table = newObitTableOTFIndexValue ("Index table", (Obit*)outData, &ver, access, err);
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

