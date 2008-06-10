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
#include "ObitOTFCal.h"
#include "ObitOTFUtil.h"
#include "ObitIOOTFFITS.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include <errno.h>

/* internal prototypes */
/* Get inputs */
ObitInfoList* X2OTFin (int argc, char **argv, ObitErr *err);
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

/* Program globals */
/* input file pointer */
FILE *myFile = NULL;
/* Input file line buffer */
gchar myLine[120];
/* Cal equivalent in Jy */
odouble calValue;
/* Reference date as yyyy-mm-dd */
gchar refDate[20];
/* Reference mjd */
odouble refMJD;

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Convert to Obit/OTF FITS tables format                               */
/*----------------------------------------------------------------------- */
{
  olong  disk, nrec, i, ierr=0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *outData= NULL;
  ObitErr *err= NULL;
  ObitIOAccess access;
 gchar *FITSdir[] = {"FITSdata/"};
  gchar infile[128], outfile[128];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean done;

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("X2OTF", 1, 0, 0, NULL, 1, FITSdir, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  ierr = 0;
  myInput = X2OTFin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input file name */
  for (i=0; i<128; i++) infile[i] = 0;
  ObitInfoListGet(myInput, "infile", &type, dim, infile, err);
  infile[dim[0]] = 0;  /* null terminate */

  /* output FITS file name */
  for (i=0; i<128; i++) outfile[i] = 0;
  ObitInfoListGet(myInput, "outfile", &type, dim, outfile, err);
  outfile[dim[0]] = 0;  /* null terminate */

  /* Cal value in Jy */
  ObitInfoListGet(myInput, "calValue", &type, dim, &calValue, err);

  /* Create ObitOTF for data */
  outData = newObitOTF("Output data");

  /* Get header info, array geometry */
  GetHeader (outData, infile, myInput, err);
  
  /* Define output, I/O size */
  disk = 1;
  nrec = 2; /* cal on and cal off separate */
  ObitOTFSetFITS(outData,nrec,disk,outfile,err);
  
   /* show any errors */
   if (err->error) ierr = 1;
   ObitErrLog(err);
   if (ierr!=0) return ierr;
   
   /* Open output OTF - try Read/Write and if that fails Write only */
   access = OBIT_IO_ReadWrite;
   if ((ObitOTFOpen (outData, access, err) 
       != OBIT_IO_OK) || (err->error>0)) {
    ObitErrClear(err);  /* don't report these */
    access = OBIT_IO_WriteOnly;  /* try again */
    if ((ObitOTFOpen (outData, access, err) 
	 != OBIT_IO_OK) || (err->error>0))  /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", outfile);
  }
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;

  /* Loop translating file */
  done = FALSE;
  while (!done) {
    done = GetRecord (outData, infile, myInput, err);
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
   
   return ierr;
} /* end of main */

ObitInfoList* X2OTFin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="X2OTF.in", *arg;
  ObitInfoList* list;

  /* command line arguments */
  if (argc<=1) Usage(); /* must have arguments */
  /* parse command line */
  for (ax=1; ax<argc; ax++) {
    arg = argv[ax];
    if (strcmp(arg, "-input") == 0){ /* input parameters */
      input_file = argv[++ax];
    } else { /* unknown argument */
      Usage();
    }
  }
  
  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* parse input file */
  ObitParserParse (input_file, list, err);

  return list;
} /* end X2OTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: X2OTF -input file\n");
    fprintf(stderr, "Convert an external file format to Obit/OTF\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def X2OTF.in\n");
    
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
  ObitInfoList *out = newObitInfoList();

  /* add parser items */
  /* output FITS file name */
  strTemp = "OTFdata.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "infile", OBIT_string, dim, strTemp, err);

  /* output FITS file name */
  strTemp = "DataOTF.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "outfile", OBIT_string, dim, strTemp, err);

  /* reference date */
  strTemp = "20030500";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "refDate", OBIT_string, dim, strTemp, err);

  /* Jy equivalent of cal */
  dim[0] = 1;
  dtemp = 1.0;
  ObitInfoListPut (out, "calValue", OBIT_double, dim, &dtemp, err);

  return out;
} /* end defaultInputs */

void GetHeader (ObitOTF *outData, gchar *infile, ObitInfoList *myInput, 
		ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get header information from text file/parse list                      */
/*   Input:                                                               */
/*      outData  Output OTF object                                        */
/*      infile   root of input file names                                 */
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

 /* Open input text file */
  myFile = fopen (infile, "rt");
  if (myFile==NULL) {
    Obit_log_error(err, OBIT_Error, "ERROR opening file %s", infile);
    errMsg = strerror(errno);
    Obit_log_error(err, OBIT_Error, "%s", errMsg);
    return;
  }

  desc = outData->myDesc;

  /* Set Feed Array Geometry */
  /* create Array geometry with 2 (X, Y pol) elements */
  numberDetect = 2;
  geom = ObitOTFArrayGeomCreate (numberDetect);
  geom->azOffset = g_realloc(geom->azOffset, numberDetect*sizeof(ofloat));
  geom->elOffset = g_realloc(geom->elOffset, numberDetect*sizeof(ofloat));
 
  /* Set feed offsets */
  for (i=0; i<numberDetect; i++) geom->azOffset[i] = 0.0;
  for (i=0; i<numberDetect; i++) geom->elOffset[i] = 0.0;

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
  geom = ObitOTFArrayGeomUnref(geom);

   
  /* Initialize  Descriptor */
  /* Beam size (GBT at 20 cm) */
  desc->beamSize = 9.0/60.0;

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

  /* Index the descriptor */
  ObitOTFDescIndex (desc);

} /* end GetHeader */

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
  gchar *errMsg;
  olong scan, sample, Xcaloff, Xcalon, Ycaloff, Ycalon;
  ofloat *data, ra, dec, az, el, nscan;
  odouble mjd;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitOTFIsA(outData));

  desc = outData->myDesc;
  geom = outData->geom;

  /* read data */
  fgets(myLine, 120, myFile);
  if (ferror(myFile)) {
    Obit_log_error(err, OBIT_Error, "ERROR opening file %s", infile);
    errMsg = strerror(errno);
    Obit_log_error(err, OBIT_Error, "%s", errMsg);
    return TRUE;
  }
  done = feof(myFile); /* hit EOF? */

  /* Parse text line */
  nscan = sscanf (myLine,"%d %d %lf %f %f %f %f %d %d %d %d",
		  &scan, &sample, &mjd, &ra, &dec, &az, &el, 
		  &Xcaloff, &Xcalon, &Ycaloff, &Ycalon);
  if (nscan<11) {
    Obit_log_error(err, OBIT_Error, "ERROR parsing line %s", myLine);
    errMsg = strerror(errno);
    Obit_log_error(err, OBIT_Error, "%s", errMsg);
    return TRUE;
  }

  /* convert units */
  
  /* Fill record entry in data */
  /* Cal off */
  mjd /= 86400.0;  /* MJD to days */
  ra  *= RAD2DG;   /* RA to deg */
  if (ra<0.0)   ra += 360.0;
  if (ra>360.0) ra -= 360.0;
  dec *= RAD2DG;   /* dec to deg */

  data = outData->buffer;
  data[desc->iloct]   = mjd - refMJD; /* Time (days) */
  data[desc->ilocti]  = 0.0;          /* time interval (days) */
  data[desc->iloctar] = 1.0;          /* target number */
  /* debug - write el as target */
  data[desc->iloctar] = el*RAD2DG;          /* target number */
  data[desc->ilocra]  = ra;           /* RA  in deg. */
  data[desc->ilocdec] = dec;          /* Dec in deg.*/
  data[desc->iloccal] = 0.0;          /* No cal */
  data[desc->ilocrot] = 0.0;          /* Parallactic angle (deg) */
  data[desc->ilocscan]= (ofloat)scan; /* scan number */
  
  /* set data values */
  data[desc->ilocdata+0] = Xcaloff;
  data[desc->ilocdata+1] = Ycaloff;
  
  /* Cal on */
  data += desc->lrec;  /* next record */
  data[desc->iloct]   = mjd - refMJD; /* Time (days) */
  data[desc->ilocti]  = 0.0;          /* time interval (days) */
  data[desc->iloctar] = 1.0;          /* target number */
  /* debug - write el as target */
  data[desc->iloctar] = el*RAD2DG;          /* target number */
  data[desc->ilocra]  = ra;           /* RA  in deg. */
  data[desc->ilocdec] = dec;          /* Dec in deg.*/
  data[desc->iloccal] = calValue;     /* cal */
  data[desc->ilocrot] = 0.0;          /* Parallactic angle (deg) */
  data[desc->ilocscan]= (ofloat)scan; /* scan number */
  
  /* set data values */
  data[desc->ilocdata+0] = Xcalon;
  data[desc->ilocdata+1] = Ycalon;
  
  /* set number of records */
  desc->numRecBuff = 2;

  /* Append to end of file */
  desc->firstRec = desc->nrecord + 1;

  return done;
} /* end GetRecord */

