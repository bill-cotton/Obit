/* To do:
  0) 
  1) ??? 1 detector and 2 poln???, get poln from IF/IF table,
     can also get Center Sky freq.
  6) Generic Obit read/Write FITS keywords - replace explicit cfitsio
  7) History
*/
/* $Id$                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003                                               */
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
/*;  Correspondence this software should be addressed as follows:     */
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
#include "ObitTableOTFTarget.h"
#include "ObitTableOTFIndex.h"
#include "ObitTableGBTANTPOSGR.h"
#include "ObitTableGBTANTPOSPF.h"
#include "ObitTableGBTSPDATA.h"
#include "ObitGBTIFInfo.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* SP2OTFin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(gchar *scan_name, ObitErr *err);
/* Get pointing times from Antenna file, Secondary focus */
void GetAntennaGR (gchar *infile, ObitInfoList *myInput, ObitErr *err);
/* Get pointing times from Antenna file, Prime Focus */
void GetAntennaPF (gchar *infile, ObitInfoList *myInput, ObitErr *err);
/* Get pointing for a given time */
void GetPoint (odouble time, ofloat *ra, ofloat *dec);
/* Get file descriptor */
gboolean GetHeader (ObitOTF *outData, char *outfile, gchar *infile, 
		    ObitInfoList *myInput, ObitErr *err);
/* Get data */
void GetData (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, ObitErr *err);
/* Get target id */
olong GetTarget (ObitOTF *outData, gboolean isNew, gchar *name, ObitErr *err);
/* Initialize Index table */
void InitScan (ObitOTF *outData, gboolean isNew, ofloat scan, ofloat target, 
	       ObitErr *err);
/* Undate scan info */
void SetScan (ObitOTF *outData, odouble startTime, odouble endTime, 
	      olong startRec, olong endRec, ObitErr *err);


/* Program globals */
gchar DataRoot[128]; /* Root directory of input data */
gchar Name[48];   /* Target name */
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
olong  nfeed=1;   /* Number of feeds */
olong  nchan=1024;/* Number of frequencies */
olong  nstok=2;   /* Number of Stokes */
gchar  poln[2];   /* Polarization states, 'R', 'L', 'X', 'Y' */
ofloat deltaFreq; /* Channel increment */
ofloat refPixel;  /* Frequency reference pixel */
odouble refFrequency; /* reference frequency (Hz) */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Convert to GBT SpectralProcessor data to Obit/OTF FITS tables format */
/*----------------------------------------------------------------------- */
{
  olong  i, disk, nrec, ierr=0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *outData= NULL;
  ObitIOAccess access;
  ObitErr *err= NULL;
  gchar *FITSdir[] = {"FITSdata/", NULL};
  gchar infile[128], outfile[128];
  ObitInfoType type;
  gboolean isNew;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  err = newObitErr();

  /* Startup - parse command line */
  ierr = 0;
  myInput = SP2OTFin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input FITS file name */
  for (i=0; i<128; i++) infile[i] = 0;
  ObitInfoListGet(myInput, "scan", &type, dim, infile, err);
  infile[dim[0]] = 0;  /* null terminate */

  /* output FITS file name */
  for (i=0; i<128; i++) outfile[i] = 0;
  ObitInfoListGet(myInput, "outfile", &type, dim, outfile, err);
  outfile[dim[0]] = 0;  /* null terminate */

  /* Get input data file root */
  for (i=0; i<128; i++) DataRoot[i] = 0;
  ObitInfoListGet(myInput, "DataRoot", &type, dim, DataRoot, err);
  DataRoot[dim[0]] = 0;  /* null terminate */

  /* Initialize Obit */
  FITSdir[1] =  DataRoot;  /* Input FITS data directory */
  mySystem = ObitSystemStartup ("SP2OTF", 1, 0, 0, NULL, 2, FITSdir, (oint)TRUE, 
				(oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Create ObitOTF for data */
  outData = newObitOTF("Output data");

  /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
   
  /* Get header info, array geometry, initialize output if necessary */
  isNew = GetHeader (outData, outfile, infile, myInput, err);
  
  /* Get Telescope position array, Prime and secondary focus different */
  if (refFrequency>1.0e9)
    GetAntennaGR (infile, myInput, err);
  else  /* Low frequencies at prime focus */
    GetAntennaPF (infile, myInput, err);
  /* Say what went wrong if error */
  if (err->error) 
     Obit_log_error(err, OBIT_Error, "Error reading Antenna file for scan %s", infile);
  
  /* Informative message */
  /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
   
  Obit_log_error(err, OBIT_InfoErr, "Adding scan %s %6.0f to OTF file %s", infile, scan, 
		 outfile);
  ObitErrLog(err);

  /* Define output, I/O size */
  disk = 1;
  nrec = 2;
  ObitOTFSetFITS(outData,nrec,disk,outfile,err);
  
  /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
   
  /* Open output OTF */
  if (isNew) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  if ((ObitOTFOpen (outData, access, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", outfile);
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;

  /* Get target id number */
  target = (ofloat)GetTarget (outData, isNew, Name, err);

  /* Initialize scan in Index table */
  InitScan (outData, isNew, scan, target, err);

  /* convert data  */
  GetData (outData, infile, myInput, err);
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;

  /* Update index table */
  SetScan (outData, startTime, endTime, startRec, endRec, err);

  /* show any errors */
   if (err->error) ierr = 1;
   ObitErrLog(err);
   if (ierr!=0) return ierr;
   
  /* Close */
  if ((ObitOTFClose (outData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing output file");
  
  /* show any errors */
   if (err->error) ierr = 1;
   ObitErrLog(err);
   if (ierr!=0) return ierr;
   
  
   /* Shutdown Obit */
   mySystem = ObitSystemShutdown (mySystem);
   
   /* cleanup */
   myInput = ObitInfoListUnref(myInput);  /* delete input list */
   outData = ObitUnref(outData);
   if (AntDMJD) g_free(AntDMJD);
   if (AntRA)   g_free(AntRA);
   if (AntDec)  g_free(AntDec);
  
   return ierr;
} /* end of main */

ObitInfoList* SP2OTFin (int argc, char **argv, ObitErr *err)
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
  olong ax;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *input_file="SP2OTF.in", *arg;
  gchar *scan_name ="unspecified";
  gboolean init=FALSE;
  gchar *strTemp;
  ObitInfoList* list;

  /* Make default inputs InfoList */
  list = defaultInputs(scan_name, err);

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

    } else if (strcmp(arg, "-scan") == 0){ /* scan name */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "scan", OBIT_string, dim, strTemp);

    } else { /* unknown argument */
      Usage();
    }
  }
  
  /* Read defaults if no file specified */
  if (!init) ObitParserParse (input_file, list, err);

  return list;
} /* end SP2OTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: SP2OTF -input file [-scan date/time]\n");
    fprintf(stderr, "Convert an GBT SpectralProcessor file format to Obit/OTF\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def SP2OTF.in\n");
    fprintf(stderr, "  -scan date/time used for form scan FITS file names\n");
    
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
  ObitInfoList *out = newObitInfoList();

  /* add parser items */
  /* base of scan file names */
  dim[0] = strlen (scan_name);
  ObitInfoListPut (out, "scan", OBIT_string, dim, scan_name, err);

  /* output FITS file name */
  strTemp = "DataOTF.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "outfile", OBIT_string, dim, strTemp, err);

  /* root of data directory */
  strTemp = "dataRoot";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "DataRoot", OBIT_string, dim, strTemp, err);

  return out;
} /* end defaultInputs */

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

gboolean GetHeader (ObitOTF *outData, gchar *outfile, gchar *infile, 
		    ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get header information from scan header files                         */
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
  ObitGBTIFInfo* IFdata=NULL;
  olong numberDetect, ndetec, ncol, i, ncopy, iscan;
  gchar Date[48], *fullname=NULL, FullFile[128], commnt[81];
  odouble SiteLat, SiteLong, SiteElev, JD, T, GMST0, GSTIAT0, nu, e2;
  odouble aEarth = 6378137.0; /* Semi major axis of Earth */
  odouble flat   = 1.0 / 298.257223563; /* inverse of flattening of Earth */ 
  ofloat UT1UTC, PolarX, PolarY;
  olong disk, ierr = 0;
  gchar backend[21];
  fitsfile *fptr;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitOTFIsA(outData));
  g_assert(outfile!=NULL);
  g_assert(infile!=NULL);
  g_assert(myInput!=NULL);

  desc = outData->myDesc;

  /* Get frequency setup */
  disk = 2;  /* input "disk" */
  IFdata = newObitGBTIFInfoValue("IFdata", "SpectralProcessor", disk,
				 DataRoot, infile, err);
  /* Save data */
  nchan        = IFdata->nchan;
  poln[0]      = IFdata->poln[0];
  poln[1]      = IFdata->poln[1];
  deltaFreq    = IFdata->delta[0];
  refPixel     = IFdata->refPixel[0];
  refFrequency = IFdata->refFrequency[0];

  IFdata = ObitGBTIFInfoUnref(IFdata); /* Cleanup */

  /* Set Feed Array Geometry */
  /* create Array geometry with 2*nchan (R,L) elements */
  numberDetect = nstok * nfeed * nchan;
  geom = ObitOTFArrayGeomCreate (numberDetect);
  geom->azOffset = g_realloc(geom->azOffset, numberDetect*sizeof(ofloat));
  geom->elOffset = g_realloc(geom->elOffset, numberDetect*sizeof(ofloat));
 
  /* Only one feed but each poln/freq has an offset in the table */
  for (i=0; i<numberDetect; i++) {
    geom->azOffset[i] = 0.0;
    geom->elOffset[i] = 0.0;
  }

  /* Other information - Get from Antenna file */
  /* get full file name */
  sprintf (FullFile,"%sAntenna/%s.fits", DataRoot, infile);

  /* Get header keywords from naked cfitsio */
   if ( fits_open_file(&fptr, FullFile, READONLY, &ierr) ) {
      Obit_log_error(err, OBIT_Error, "ERROR %d opening input FITS file %s", ierr, FullFile);
      return out;
   }

   /* Read main file keywords */
   UT1UTC = 0.0;
   fits_read_key_flt (fptr, "DELTAUTC", &UT1UTC, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   PolarX = 0.0;
   fits_read_key_flt (fptr, "IERSPMX", &PolarX, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   PolarY = 0.0;
   fits_read_key_flt (fptr, "IERSPMY", &PolarY, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   SiteLat = 3.8433119e+01;
   fits_read_key_dbl (fptr, "SITELAT", &SiteLat, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   SiteLong = 7.9839833E+01 ;
   fits_read_key_dbl (fptr, "SITELONG", &SiteLong, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   SiteElev = 8.24595E+02;
   fits_read_key_dbl (fptr, "SITEELEV", &SiteElev, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   iscan = 1;
   fits_read_key_lng (fptr, "SCAN", &iscan, commnt, &ierr);
   if (ierr==KEY_NO_EXIST) ierr = 0;

   fits_read_key_str (fptr, "DATE-OBS", Date, commnt, &ierr);
   fits_read_key_str (fptr, "OBJECT", Name, commnt, &ierr);

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
  geom->lat = SiteLat * RAD2DG;
  /* telescope longitude in radians */
  geom->lon = SiteLong * RAD2DG;
  /* LST at iat0 in radians */
  geom->LSTiat0 = geom->GSTiat0*1.74533e-2 + geom->lon;
  /* Earth rotation rate in rad/day */
  geom->RadDay = geom->DegDay*1.74533e-2;
  /* Data - IAT in days */
  geom->dataIat = (geom->dataUtc - geom->iatUtc) / 86400.0;

  /* Attach Array geometry to OTF */
  outData->geom = ObitOTFArrayGeomRef(geom);
  geom = ObitOTFArrayGeomUnref(geom);

  /* Other information from SP file */
  /* get full file name */
  sprintf (FullFile,"%sSpectralProcessor/%s.fits", DataRoot, infile);

  /* Get header keywords from naked cfitsio */
   if ( fits_open_file(&fptr, FullFile, READONLY, &ierr) ) {
      Obit_log_error(err, OBIT_Error, "ERROR %d opening input FITS file %s", ierr, FullFile);
      return out;
   }

   /* Make sure data in "SpectralProcessor" mode */
   fits_read_key_str (fptr, "BACKEND", backend, commnt, &ierr);
   if (strncmp (backend, "SPAB",17)) {
     Obit_log_error(err, OBIT_Error, "ERROR backend %s NOT SpectralProcessor", backend);
     return out;
   }

   /* Read main file keywords */
   ndetec = 2;
   /* Spectral processor can only handle 2*/
   nfeed = ndetec/2;

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

   
  /* Does the input file exist?  */
  disk = 1;  /* input "disk" */
  fullname = ObitFITSFilename (disk, outfile, err);

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

    /* Frequency axis */
    desc->inaxes[desc->naxis] = nchan;
    strncpy (desc->ctype[desc->naxis], "FREQ    ", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = deltaFreq;
    desc->crpix[desc->naxis] = refPixel;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = refFrequency;
    desc->naxis++;

    /* Feed axis */
    desc->inaxes[desc->naxis] = nfeed;
    strncpy (desc->ctype[desc->naxis], "FEED", OTFLEN_KEYWORD);
    desc->cdelt[desc->naxis] = 1.0;
    desc->crpix[desc->naxis] = 1.0;
    desc->crota[desc->naxis] = 0.0;
    desc->crval[desc->naxis] = 1.0;
    desc->naxis++;

    /* Reference Modified Julian date */
    refMJD = JD - 2400000.5;

  } /* end initialize descriptor */

  /* cleanup */
  if (fullname) g_free(fullname);

  /* Index the descriptor */
  ObitOTFDescIndex (desc);

  return out;
} /* end GetHeader */

void GetData (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Read data from GB FITS files and write                                */
/*  Assumptions:                                                          */
/*    The following assumes that the data were recorded with two states   */
/*    per integration, the firs being call off and the second cal on.     */
/*    The details are contained in the STATE table in the data FILE.      */
/*    There is one row per state and the CAL row is 1 if the cal is on    */
/*    else 0.                                                             */
/*      outData  Output OTF object                                        */
/*      myInput  parser object                                            */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  gchar FullFile[128];
  olong irow, ichan;
  olong disk, ver, nrow, ip1off, ip2off, incdatawt;
  gchar *tab;
  ofloat *data, ra, dec;
  odouble JD=0.0, timetag=0.0;
  ObitTableGBTSPDATA* table;
  ObitTableGBTSPDATARow* row;
  ObitOTFDesc *desc;
  ObitOTFArrayGeom *geom;
  ObitIOCode retCode;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(infile!=NULL);
  g_assert(myInput!=NULL);

  /* Stokes offsets in input */
  ip1off = 0;
  if (poln[0]=='R') ip1off = 0;
  else if (poln[0]=='L') ip1off = 2;
  else if (poln[0]=='X') ip1off = 0;
  else if (poln[0]=='Y') ip1off = 2;

  ip2off = 0;
  if (poln[1]=='R') ip2off = 0;
  else if (poln[1]=='L') ip2off = 2;
  else if (poln[1]=='X') ip2off = 0;
  else if (poln[1]=='Y') ip2off = 2;

  /* Get reference MJD, convert ref date  to JD */
  ObitOTFDescDate2JD (outData->geom->RefDate, &JD);
  /* Reference Modified Julian date */
  refMJD = JD - 2400000.5;


  /* get full file name */
  sprintf (FullFile,"SpectralProcessor/%s.fits", infile);
  
  /* Create table structure */
  table = newObitTableGBTSPDATA("Data");
  if (err->error) return;

  /* Setup */
  disk = 2;
  tab = "DATA";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(table,disk,FullFile,tab,ver,nrow,err);

  /* Open */
  retCode = ObitTableGBTSPDATAOpen (table, OBIT_IO_ReadOnly, err);
  if (err->error) return;

  integTime = table->inttime ; /* Integration time in sec */
  /* Convert from seconds to days */
  integTime /= 86400.0;

  /* Create Row structure */
  row = newObitTableGBTSPDATARow (table);

  desc = outData->myDesc;
  geom = outData->geom;
  incdatawt = desc->incdatawt; /* increment in data-wt axis */

  /* Write at end */
  desc->firstRec = desc->nrecord+1;
  startRec = desc->firstRec;  /* First record in scan */
  startTime = -1.0e20;        /* dummy srart time */

  /* Loop over table */
  for (irow = 1; irow<=table->myDesc->nrow; irow++) {
    retCode = ObitTableGBTSPDATAReadRow (table, irow, row, err);
    if (err->error) return;

    /* first time in scan */
    timetag = row->utdate + row->utcstart/86400.0; /* Time in days */
    if (startTime<-1000.0) startTime = timetag - refMJD;

    /* Fill record entry in data */
    data = outData->buffer;
    data[desc->iloct]   = timetag - refMJD;    /* Time (days) */
    /* Correct for position in integration */
    data[desc->iloct] -= 0.25 * integTime;
    /* Get position */
    GetPoint (timetag-0.25*integTime, &ra, &dec);
    data[desc->ilocti]  = 0.0;          /* time interval (days) */
    data[desc->iloctar] = target;       /* target number */
    data[desc->ilocscan]= scan;         /* scan number */
    data[desc->ilocra]  = ra;           /* RA  in deg. */
    data[desc->ilocdec] = dec;          /* Dec in deg.*/
    data[desc->iloccal] = 0.0;          /* No cal yet */
    data[desc->ilocrot] = 0.0;          /* Parallactic angle (deg) */
    
    /* The following makes assumptions about the order of the data */
    /* set data values, cal off */
    /* Loop over channels */
    for (ichan=0; ichan<nchan; ichan++) {
      data[desc->ilocdata+ichan*2]   = row->data[ip1off*nchan+ichan]; /* RCP or X */
      data[desc->ilocdata+1+ichan*2] = row->data[ip2off*nchan+ichan]; /* LCP or Y */
    }
    
    /* Cal on */
    data += desc->lrec;  /* next record */
    data[desc->iloct]   = timetag - refMJD; /* Time (days) */
    /* Correct for position in integration */
    data[desc->iloct] += 0.25 * integTime;
    /* Get position */
    GetPoint (timetag+0.25*integTime, &ra, &dec);
    data[desc->ilocti]  = 0.0;          /* time interval (days) */
    data[desc->iloctar] = target;       /* target number */
    data[desc->ilocra]  = ra;           /* RA  in deg. */
    data[desc->ilocdec] = dec;          /* Dec in deg.*/
    data[desc->iloccal] = 1.0;          /* cal */
    data[desc->ilocrot] = 0.0;          /* Parallactic angle (deg) */
    data[desc->ilocscan]= scan;         /* scan number */
    
    /* set data values */
    /* Loop over channels */
    for (ichan=0; ichan<nchan; ichan++) {
      data[desc->ilocdata+(ichan*2)*incdatawt]     = row->data[(1+ip1off)*nchan+ichan]; /* RCP or X */
      data[desc->ilocdata+(ichan*2)*incdatawt+1]   = 1.0; /* Initial weight */
      data[desc->ilocdata+(1+ichan*2)*incdatawt]   = row->data[(1+ip2off)*nchan+ichan]; /* LCP or Y */
      data[desc->ilocdata+(1+ichan*2)*incdatawt+1] = 1.0; /* Initial weight */
    }
  
    /* set number of records */
    desc->numRecBuff = 2;
    
    /* Write output  */
   if ((ObitOTFWrite (outData, NULL, err) != OBIT_IO_OK) || (err->error>0))
      Obit_log_error(err, OBIT_Error, "ERROR writing output Table file");
   } /* end loop over table */
  
  /* Get end times and record numbers */
  desc->firstRec = desc->nrecord+1;
  endTime = timetag - refMJD;
  endRec  = desc->nrecord;  /* Last record in scan */

  /* Close */
  retCode = ObitTableGBTSPDATAClose (table, err);
  if (err->error) return;
  
  /* Cleanup */
  table = ObitTableGBTSPDATAUnref(table);
  row = ObitTableGBTSPDATARowUnref(row);


} /* end GetData  */

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
    w1 = (time-AntDMJD[best+1]) / (AntDMJD[best+1]-AntDMJD[best]);
    w2 = 1.0 - w1;
    *ra  = w1 * AntRA[best]  + w2 * AntRA[best+1];
    *dec = w1 * AntDec[best] + w2 * AntDec[best+1];
  }
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

