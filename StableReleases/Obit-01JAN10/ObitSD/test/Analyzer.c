#include <stdio.h>
#include <stdlib.h>
#include "ObitAll.h"
#include "ObitOTF.h"
#include "ObitOTFUtil.h"
#include "ObitOTFGetSoln.h"
#include "ObitPennArrayUtil.h"
#include "ObitOTFSoln2Cal.h"

/* Testbed data analyser for GBT Penn Array */

int main ( int argc, char **argv )
{
  ObitSystem *mySystem;
  ObitOTF *simData, *scrData;
  ObitImage  *outImage=NULL;
  ObitErr *err;
  ObitTableOTFSoln *solnTable;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong convType;
  /*  gchar *Filename="GBTSimulate.fits";*/
  gchar *Filename="PASimOTF.fits";
  gchar *FITSdir[] = {"FITSdata/"};
  oint itemp, ncoef;
  olong  ver;
  olong nx, ny, iter, niter, solnuse, calin, calout;
  gboolean doCalSelect, bad, bail;
  ofloat pMode, pRMS, xCells, yCells, SI, minWt, rTemp, clip;
  gchar *calType="Offset", *outFilename="!Analyzer.fits";

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("Analyzer", 1, 0, 0, NULL, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  ObitErrLog(err); /* show any error messages on err */

  /* Create ObitOTF for data */
  simData = newObitOTF("Penn Array Simulation");

  /* Define output, I/O size */
  ObitOTFSetFITS(simData,1000,1,Filename,err);

  /* Open/close input OTF to fully instantiate */
  if ((ObitOTFOpen (simData, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input FITS file %s", Filename);
  
  /* delete any Soln/Cal Tables */
  ObitOTFZapTable (simData, "OTFSoln", -1, err);
  ObitOTFZapTable (simData, "OTFCal", -1, err);
  ObitErrLog(err);

  /* Close */
  if ((ObitOTFClose (simData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing input file");
  
  /* Initial (dummy) calibration table 
     Solution interval  */
  SI = 2.0 / 86400.0;
  ObitInfoListPut(simData->info, "SOLINT", OBIT_float, dim, (gpointer*)&SI, err);
  calin = 1;
  ver = calin;
  ncoef=1;
  ObitOTFGetDummyCal (simData, simData, ver, ncoef, err);
  if (err->error) 
    Obit_log_error(err, OBIT_Error, "ERROR calibrating FITS file %s", Filename);

  /* Solution interval 30.0  sec */
  SI = 30.0 / 86400.0;
  ObitInfoListPut(simData->info, "SOLINT", OBIT_float, dim, (gpointer*)&SI, err);

  /* Determine instrumental calibration */
  solnTable = ObitOTFGetInstCal (simData, simData, err);
  if (err->error) 
    Obit_log_error(err, OBIT_Error, "ERROR calibrating FITS file %s", Filename);
  
  /* apply calibration to data */
  doCalSelect = TRUE;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (simData->info, "doCalSelect", OBIT_bool, dim, (gpointer)&doCalSelect, err);

  /* Apply solution to calibration
     Specify desired tables */
  dim[0] = 1;
  solnuse = 0;
  ObitInfoListPut(simData->info, "SOLNUSE",  OBIT_long, dim, (gpointer*)&solnuse, err);
  ObitInfoListPut(simData->info, "CALIN",    OBIT_long, dim, (gpointer*)&calin,   err);
  calout = calin+1;
  ObitInfoListPut(simData->info, "CALOUT",   OBIT_long, dim, (gpointer*)&calout,  err);
  calin = calout;
 
   /* Update calibration */
  ObitOTFSoln2Cal (simData, simData, err);
  if (err->error) 
    Obit_log_error(err, OBIT_Error, "ERROR calibrating FITS file %s", Filename);

  /* show any errors */
  bail = err->error;
  ObitErrLog(err);
  if (bail) return 1;

  itemp = 1;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (simData->info, "DOCALIB", OBIT_long, dim, 
		   (gpointer)&itemp, err);

  /* Most recent table */
  itemp = 0;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (simData->info, "GAINUSE", OBIT_long, dim, 
		   (gpointer)&itemp, err);

  /* show any errors */
  ObitErrLog(err);

  /* Remove Median over array */
  solnTable = ObitOTFGetSolnMBBase (simData, simData, err);
  if (err->error) 
    Obit_log_error(err, OBIT_Error, "ERROR calibrating FITS file %s", Filename);

  /* Apply solution to calibration
     Specify desired tables */
  dim[0] = 1;
  solnuse = 0;
  ObitInfoListPut(simData->info, "SOLNUSE",  OBIT_long, dim, (gpointer*)&solnuse, err);
  ObitInfoListPut(simData->info, "CALIN",    OBIT_long, dim, (gpointer*)&calin,   err);
  calout = calin+1;
  ObitInfoListPut(simData->info, "CALOUT",   OBIT_long, dim, (gpointer*)&calout,  err);
  calin = calout;
 
   /* Update calibration */
  ObitOTFSoln2Cal (simData, simData, err);
  if (err->error) 
    Obit_log_error(err, OBIT_Error, "ERROR calibrating FITS file %s", Filename);

  /* show any errors */
  bail = err->error;
  ObitErrLog(err);
  if (bail) return 1;

  /* Set imaging parameters */
  nx = 1000;
  ObitInfoListPut(simData->info, "nx", OBIT_long, dim, (gpointer*)&nx, err);
  
  ny = 1000;
  ObitInfoListPut(simData->info, "ny", OBIT_long, dim, (gpointer*)&ny, err);
  
  xCells = -(8.0/5.0) / 3600.0; /* ~5 samples per beam */
  xCells = -8.0 / 5.0;
  ObitInfoListPut(simData->info, "xCells", OBIT_float, dim, (gpointer*)&xCells, err);
  
  yCells = (8.0/5.0) / 3600.0; /* ~5 samples per beam */
  yCells = 8.0 / 5.0;
  ObitInfoListPut(simData->info, "yCells", OBIT_float, dim, (gpointer*)&yCells, err);
  
  /* Center */
  rTemp = 22.0; /* RA in deg */
  ObitInfoListPut(simData->info, "RA", OBIT_float, dim, (gpointer*)&rTemp, err);
  rTemp = 0.0; /* dec in deg. */
  ObitInfoListPut(simData->info, "Dec", OBIT_float, dim, (gpointer*)&rTemp, err);

  minWt = 0.3;
  ObitInfoListPut(simData->info, "minWt", OBIT_float, dim, (gpointer*)&minWt, err);
  
  convType = 5; /* Sph wave fn convolving fn. */
  ObitInfoListPut(simData->info, "ConvType", OBIT_long, dim, (gpointer*)&convType, err);
 
  /* Create output image object from OTF */
  outImage = ObitOTFUtilCreateImage (simData, err);

  /* define output */
  ObitImageSetFITS(outImage,OBIT_IO_byPlane,1,outFilename,blc,trc,err);
    
  /* Open and close to fully instantiate */
  ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  ObitImageClose (outImage, err);

   /* show any errors */
  ObitErrLog(err);

  /* Form image */
  ObitOTFUtilMakeImage (simData, outImage, TRUE, NULL, err);

  /* show any errors */
  ObitErrLog(err);

  /* Scratch file for calibration */
  scrData = newObitOTFScratch (simData, err);
  
  /* Calibration type */
  dim[0] = strlen(calType);
  ObitInfoListPut(scrData->info, "calType",   OBIT_string, dim, calType, err);

  /* show any errors */
  bail = err->error;
  ObitErrLog(err);
  if (bail) return 1;

  /* Iterate a few times */
  niter = 5;
  bad = FALSE;
  /* debug
  return 0 ;*/
  for (iter = 0; iter<niter; iter++) {
    Obit_log_error(err, OBIT_InfoErr, "Begin iteration %d ",iter+1);

    /* Get statistics on image */
    pMode = ObitFArrayMode(outImage->image);
    pRMS  = ObitFArrayRMS(outImage->image);
    Obit_log_error(err, OBIT_InfoErr, "Image mode %g  RMS %g ",pMode, pRMS);
    /* show any errors/messages */
    ObitErrLog(err);

    /* clip image below 3 sigma after first round  */
    /*clip = (5.0 + 2*(10-iter))*pRMS;*/
    clip = (3.0)*pRMS;
    clip = MAX (clip, 1.0e-5);
    if (iter<2) clip = 0.5e20;
    Obit_log_error(err, OBIT_InfoErr, "Clip image below %g ",clip);
    ObitFArrayClip (outImage->image, clip, 1.0e20, 0.0);
    /* show any errors */
    bad = bad || err->error;
    ObitErrLog(err);

    /* Calibrate before subtraction */
    itemp = 1;
    dim[0] = 1; dim[1] = 1;
    ObitInfoListPut (simData->info, "DOCALIB", OBIT_long, dim, (gpointer)&itemp, err);
    doCalSelect = TRUE;
    ObitInfoListPut (simData->info, "doCalSelect", OBIT_bool, dim, (gpointer)&doCalSelect, err);
    itemp = calout;
    /* Don't use new calibration before decent model
    if (iter<=1) itemp = 2; */
    ObitInfoListPut (simData->info, "GAINUSE", OBIT_long, dim, (gpointer)&itemp, err);


    /* Subtract image from simData to scrData */
    ObitOTFUtilSubImage(simData, scrData, outImage->image, outImage->myDesc, err);
    if (err->error && (!bad))
      Obit_log_error(err, OBIT_Error, "ERROR subtracting model");
    bad = bad || err->error;
 
    /* Determine calibration */
    if ((iter==1) || (iter==3)) {
      SI = 5.0 / 86400.0;
      if (iter>=3) SI =  5.0*0.25 / 86400.0;
      if (iter>=4) SI = 5.0*0.5 / 86400.0;
      if (iter>=5) SI = 5.0*0.25 / 86400.0;
      Obit_log_error(err, OBIT_InfoErr, "Do additive calibration %f sec",SI*86400.);
      ObitInfoListPut(scrData->info, "SOLINT", OBIT_float, dim, (gpointer*)&SI, err);
      solnTable = ObitOTFGetSolnGain (scrData, simData, err);
    } else {
      SI = MAX (0.0, 0.6*(0+iter)*0.01 / 86400.0);
      SI = 10.0 /  86400.0;
      ObitInfoListPut(scrData->info, "SOLINT", OBIT_float, dim, (gpointer*)&SI, err);
      Obit_log_error(err, OBIT_InfoErr, "Do polynomial calibration %f sec",SI*86400.);
      solnTable = ObitOTFGetSolnCal (scrData, simData, err);
    }
    if (err->error && (!bad)) 
      Obit_log_error(err, OBIT_Error, "ERROR calibrating FITS file %s", Filename);
    bad = bad || err->error;
    
    /* Apply solution to calibration */
    /* Specify desired tables */
    dim[0] = 1;
    solnuse = 0;
    ObitInfoListPut(simData->info, "SOLNUSE",  OBIT_long, dim, (gpointer*)&solnuse, err);
    ObitInfoListPut(simData->info, "CALIN",    OBIT_long, dim, (gpointer*)&calin,   err);
    calout = calout+1;
    ObitInfoListPut(simData->info, "CALOUT",   OBIT_long, dim, (gpointer*)&calout,  err);
    calin = calout;
    /* if (iter>1) calin = calout; Don't use new calibration before decent model */
 
    /* Update calibration */
    ObitOTFSoln2Cal (simData, simData, err);
    if (err->error) 
      Obit_log_error(err, OBIT_Error, "ERROR calibrating FITS file %s", Filename);

    /* apply calibration to data */
    doCalSelect = TRUE;
    dim[0] = 1; dim[1] = 1;
    ObitInfoListPut (simData->info, "doCalSelect", OBIT_bool, dim, (gpointer)&doCalSelect, err);
    
    itemp = 1;
    ObitInfoListPut (simData->info, "DOCALIB", OBIT_long, dim, (gpointer)&itemp, err);

    /* Most recent table */
    itemp = calout;
    dim[0] = 1; dim[1] = 1;
    ObitInfoListPut (simData->info, "GAINUSE", OBIT_long, dim, (gpointer)&itemp, err);

    /* Form image */
    ObitOTFUtilMakeImage (simData, outImage, TRUE, NULL, err);
    if (err->error && (!bad)) 
      Obit_log_error(err, OBIT_Error, "ERROR forming image");
    bad = bad || err->error;
    
    /* show any errors */
    bad = bad || err->error;
    ObitErrLog(err);

    if (bad) break; /* stop if in error condition */
    
 } /* end iteration loop */

  /* Clean up*/
  simData = ObitUnref(simData);

  /* show any errors */
  ObitErrLog(err);

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return 0;
} /* end of main */

