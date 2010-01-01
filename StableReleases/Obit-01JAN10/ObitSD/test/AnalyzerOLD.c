#include <stdio.h>
#include <stdlib.h>
#include "ObitAll.h"
#include "ObitOTF.h"
#include "ObitOTFUtil.h"
#include "ObitOTFGetSoln.h"
#include "ObitPennArrayUtil.h"

/* Subtract image */
void ObitOTFUtilSubImage(ObitOTF *inOTF, ObitOTF *outOTF, ObitFArray *image, 
		     ObitImageDesc *desc, ObitErr *err);

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
  gchar *Filename="GBTSimulate.fits";
  gchar *FITSdir[] = {"FITSdata/"};
  oint itemp;
  olong nx, ny, iter, niter;
  gboolean doCalSelect, bad;
  ofloat pMode, pRMS, xCells, yCells, SI, minWt;
  gchar *outFilename="!Analyzer.fits";

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("Simulator", 1, 0, 0, NULL, 1, FITSdir, (oint)TRUE, (oint)FALSE, err);
  ObitErrLog(err); /* show any error messages on err */

  /* Create ObitOTF for data */
  simData = newObitOTF("GBT Simulation");

  /* Define output, I/O size */
  ObitOTFSetFITS(simData,1000,1,Filename,err);

  /* Open/close input OTF to fully instantiate */
  if ((ObitOTFOpen (simData, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input FITS file %s", Filename);
  
  /* Close */
  if ((ObitOTFClose (simData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing input file");
  
  /* Solution interval 1 sec */
  SI = 1.0 / 86400.0;
  SI = 0.4 / 86400.0;
  ObitInfoListPut(simData->info, "SOLINT", OBIT_float, dim, (gpointer*)&SI, err);

  /* Determine calibration */
  solnTable = ObitOTFGetSolnCal (simData, simData, err);
  if (err->error) 
   Obit_log_error(err, OBIT_Error, "ERROR calibrating FITS file %s", Filename);

   /* apply calibration to data */
  doCalSelect = TRUE;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (simData->info, "doCalSelect", OBIT_bool, dim, (gpointer)&doCalSelect, err);

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

  /* Set imaging parameters */
  nx = 500;
  ObitInfoListPut(simData->info, "nx", OBIT_long, dim, (gpointer*)&nx, err);
  
  ny = 500;
  ObitInfoListPut(simData->info, "ny", OBIT_long, dim, (gpointer*)&ny, err);
  
  xCells = -(4.0/5.0) / 3600.0; /* ~5 samples per beam */
  ObitInfoListPut(simData->info, "xCells", OBIT_float, dim, (gpointer*)&xCells, err);
  
  yCells = (4.0/5.0) / 3600.0; /* ~5 samples per beam */
  ObitInfoListPut(simData->info, "yCells", OBIT_float, dim, (gpointer*)&yCells, err);
  
  minWt = 0.5;
  ObitInfoListPut(simData->info, "minWt", OBIT_float, dim, (gpointer*)&minWt, err);
  
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
  ObitOTFUtilMakeImage (simData, outImage, TRUE, err);

  /* show any errors */
  ObitErrLog(err);

  /* Scratch file for calibration */
  scrData = newObitOTFScratch (simData, err);
  
  /* Solution interval  */
  SI = 1.0 / 86400.0;
  SI = 0.4 / 86400.0;
  ObitInfoListPut(scrData->info, "SOLINT", OBIT_float, dim, (gpointer*)&SI, err);

  /* show any errors */
  ObitErrLog(err);

  /* Iterate a few times */
  niter = 4;
  bad = FALSE;
  for (iter = 0; iter<niter; iter++) {
    Obit_log_error(err, OBIT_InfoErr, "Begin iteration %d ",iter+1);

    /* Get statistics on image */
    pMode = ObitFArrayMode(outImage->image);
    pRMS  = ObitFArrayRMS(outImage->image);
    Obit_log_error(err, OBIT_InfoErr, "Image mode %f  RMS %f ",pMode, pRMS);
    /* show any errors/messages */
    ObitErrLog(err);

    /* clip image below 5 sigma */
    ObitFArrayClip (outImage->image, 5.0*pRMS, 1.0e20, 0.0);
    /* show any errors */
    bad = bad || err->error;
    ObitErrLog(err);

    /* Don't calibrate before subtraction */
    itemp = -1;
    dim[0] = 1; dim[1] = 1;
    ObitInfoListPut (simData->info, "DOCALIB", OBIT_long, dim, (gpointer)&itemp, err);
    doCalSelect = FALSE;
    ObitInfoListPut (simData->info, "doCalSelect", OBIT_bool, dim, (gpointer)&doCalSelect, err);

    /* Subtract image from simData to scrData */
    ObitOTFUtilSubImage(simData, scrData, outImage->image, outImage->myDesc, err);
    if (err->error && (!bad))
      Obit_log_error(err, OBIT_Error, "ERROR subtracting model");
    bad = bad || err->error;
 
    /* Determine calibration */
    solnTable = ObitOTFGetSolnCal (scrData, simData, err);
    if (err->error && (!bad)) 
      Obit_log_error(err, OBIT_Error, "ERROR calibrating FITS file %s", Filename);
    bad = bad || err->error;
    
    /* apply calibration to data */
    doCalSelect = TRUE;
    dim[0] = 1; dim[1] = 1;
    ObitInfoListPut (simData->info, "doCalSelect", OBIT_bool, dim, (gpointer)&doCalSelect, err);
    
    itemp = 1;
    ObitInfoListPut (simData->info, "DOCALIB", OBIT_long, dim, (gpointer)&itemp, err);

    /* Form image */
    ObitOTFUtilMakeImage (simData, outImage, TRUE, err);
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

