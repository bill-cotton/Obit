#include <stdio.h>
#include <stdlib.h>
#include "ObitAll.h"
#include "ObitOTF.h"
#include "ObitOTFUtil.h"
#include "ObitOTFGetSoln.h"
#include "ObitOTFGetAtmCor.h"
#include "ObitOTFSoln2Cal.h"

/* Testbed for OTF self calibration*/

int main ( int argc, char **argv )
{
  ObitSystem *mySystem;
  ObitOTF *inData=NULL, *tmpData=NULL;
  ObitImage *Dirty=NULL, *outImage=NULL, *masterImage=NULL;
  ObitTableOTFSoln *soln=NULL;
  ObitErr *err;
  olong ierr=0;
  ObitTableOTFSoln *solnTable;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong itemp, order, ver, disk = 1 ;
  gchar *inFile     = "testQBandOTF.fits";      /* input OTF data */
  gchar *outFile    = "!testQBandCalOTF.fits";  /* output OTF data */
  gchar *tmpFile    = "tmpQBandOut.fits";       /* scratch output image */
  gchar *dirtFile   = "testQBandOut.fits";      /* output dirty image file */
  gchar *masterFile = "testQBandMaster.fits";   /* Standard lie */
  gchar *target[] = {"1256-057"};
  olong scans[] = {108,133};
  gchar *calType="MultiBeam";
  ofloat minFlux, clip, minEl, tau0, xtemp;
  ofloat aTemp[] = {1.643, 1.526, 2.999, 2.779, 1.686, 1.514, 2.964, 2.872};
  ofloat tRx[]   = {11.51, 8.601, 18.21, 14.31, 9.461, 8.278, 21.68, 15.38};
  ofloat calJy[] = {9.9, 10.1, 6.17, 5.97, 10.0, 11.18, 5.75, 5.55};
  gchar *FITSdir[] = {"../testIt"};
  oint ncoef;
  olong  lver;
  olong nx, ny, iter, niter, solnuse, calin, calout;
  gboolean doCalSelect, bad, bail;
  ofloat xCells, yCells, SI, minWt, rTemp;

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("testSelfCal", 1, 0, 0, NULL, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR initializing Obit");
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Create ObitOTF for data */
  inData = newObitOTF("input OTF");

  /* Define output, I/O size */
  ObitOTFSetFITS(inData,1000,disk,inFile,err);
  ObitOTFFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR creating input");
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* delete any Soln/Cal Tables */

  ObitOTFZapTable (inData, "OTFSoln", -1, err);
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;
  ObitOTFZapTable (inData, "OTFCal",  -1, err);
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Initial (dummy) calibration table */
  dim[0] = dim[1] = dim[2] = 1;
  SI = 1.25  / 86400.0; /* Solution interval  in sec. */
  ObitInfoListAlwaysPut(inData->info, "SOLINT", OBIT_float, dim, &SI);
  lver = 1;
  ncoef=1;
  ObitOTFGetDummyCal (inData, inData, lver, ncoef, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR Dummy Cal table");
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Atmospheric solution parameters */
  SI = 120.0 / 86400.0;/* Solution interval 120.0  sec */
  ObitInfoListAlwaysPut(inData->info, "SOLINT", OBIT_float, dim, &SI);
  tau0 = 0.053;
  ObitInfoListAlwaysPut(inData->info, "TAU0",  OBIT_float, dim, &tau0);
  ObitInfoListAlwaysPut(inData->info, "MINEL", OBIT_float, dim, &minEl);
  dim[0] = 8;
  ObitInfoListAlwaysPut(inData->info, "ATEMP", OBIT_float, dim, aTemp);
  ObitInfoListAlwaysPut(inData->info, "TRX",   OBIT_float, dim, tRx);
  ObitInfoListAlwaysPut(inData->info, "CALJY", OBIT_float, dim, calJy);

  /* Atmospheric calibration */
  fprintf (stderr, "Atmospheric calibration\n");
  soln = ObitOTFGetAtmCor(inData, inData, err);
  soln = ObitTableOTFSolnUnref(soln); 
  dim[0] = 1;
  ver = 1;
  ObitInfoListAlwaysPut(inData->info, "SOLNUSE", OBIT_long, dim, &ver);
  ver = 1;
  ObitInfoListAlwaysPut(inData->info, "CALIN",   OBIT_long, dim, &ver);
  ver = 0;
  ObitInfoListAlwaysPut(inData->info, "CALOUT",  OBIT_long, dim, &ver);
  ObitOTFSoln2Cal(inData, inData, err);  
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR Atmospheric calibration");
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;
  
  /* Linear baseline */
  fprintf (stderr, "Baseline removal\n");
  dim[0] = dim[1] = dim[2] = 1;
  SI = 10.0 / 86400.0;/* Solution interval 10.0  sec */
  ObitInfoListAlwaysPut(inData->info, "SOLINT", OBIT_float, dim, &SI);
  order = 1;
  ObitInfoListAlwaysPut(inData->info, "ORDER", OBIT_long, dim, &order);
  ver = 1;
  ObitInfoListAlwaysPut(inData->info, "GAINUSE", OBIT_long, dim, &ver);
  ver = 1;
  ObitInfoListAlwaysPut(inData->info, "DOCAL",  OBIT_long, dim, &ver);
  doCalSelect = TRUE;
  ObitInfoListAlwaysPut(inData->info, "doCalSelect",  OBIT_bool, dim, &doCalSelect);
  ver = 1;
  ObitInfoListAlwaysPut(inData->info, "FLAGVER",  OBIT_long, dim, &ver);
  dim[0] = 2;
  ObitInfoListAlwaysPut(inData->info, "SCANS",  OBIT_long, dim, scans);
  soln = ObitOTFGetSolnPolyBL (inData, inData, err);
  soln = ObitTableOTFSolnUnref(soln); 
  /* Soln2Cal parameters  */
  dim[0] = 1;
  ver = 2;
  ObitInfoListAlwaysPut(inData->info, "SOLNUSE", OBIT_long, dim, &ver);
  ver = 2;
  ObitInfoListAlwaysPut(inData->info, "CALIN",   OBIT_long, dim, &ver);
  ver = 3;
  ObitInfoListAlwaysPut(inData->info, "CALOUT",  OBIT_long, dim, &ver);
  ObitOTFSoln2Cal(inData, inData, err);  
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR in Baseline calibration");
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Imaging parameters */
  dim[0] = dim[1] = dim[2] = 1;
  ver = 0;
  ObitInfoListAlwaysPut(inData->info, "GAINUSE", OBIT_long, dim, &ver);
  xtemp =  194.04646;
  ObitInfoListAlwaysPut(inData->info, "RA",  OBIT_float, dim, &xtemp);
  xtemp = -5.78892;
  ObitInfoListAlwaysPut(inData->info, "Dec",  OBIT_float, dim, &xtemp);
  nx = 100;
  ObitInfoListAlwaysPut(inData->info, "nx",  OBIT_long, dim, &nx);
  ny = 100;
  ObitInfoListAlwaysPut(inData->info, "nx",  OBIT_long, dim, &ny);
  xCells = -4.0 / 3600.0;
  ObitInfoListAlwaysPut(inData->info, "xCells",  OBIT_float, dim, &xCells);
  yCells = 4.0 / 3600.0;
  ObitInfoListAlwaysPut(inData->info, "yCells",  OBIT_float, dim, &yCells);
  xtemp = 0.9;
  ObitInfoListAlwaysPut(inData->info, "ConvFactor",  OBIT_float, dim, &xtemp);
  xtemp = 1.0;
  ObitInfoListAlwaysPut(inData->info, "minWt",  OBIT_float, dim, &xtemp);

  /* Create initial image */
  Dirty = newObitImage("Dirty Image");
  ObitImageSetFITS(Dirty, OBIT_IO_byPlane, disk, dirtFile, blc, trc, err);
  /*ObitImageFullInstantiate (Dirty, FALSE, err);*/
  ObitOTFUtilMakeImage (inData, Dirty, TRUE, NULL, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR forming Dirty Image");
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Clean up */
  inData = ObitUnref(inData);
  Dirty  = ObitUnref(Dirty);

  /* show any errors */
  ObitErrLog(err);

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return 0;
} /* end of main */

