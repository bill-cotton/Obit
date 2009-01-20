#include <stdio.h>
#include <stdlib.h>
#include "ObitSystem.h"
#include "ObitUV.h"
#include "ObitImage.h"
#include "ObitImageUtil.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanVis.h"

/* program to test CleanVis functionality */
int main ( int argc, char **argv )
{
  ObitSystem *mySystem;
  ObitErr *err;
  olong ierr, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *AIPSdir[] = {"../AIPSdata/"};
  gchar *FITSdir[] = {"../testIt/"};
  ObitUV *uvdata=NULL, *outdata=NULL;
  ObitImage *model=NULL;
  ObitImageMosaic *mosaic=NULL;
  ObitSkyModel *SkyModel=NULL; 
  /* data */
  olong inDisk = 1;
  gchar *inFile  = "PModel.uvtab";
  gchar *inModel = "PModelI.fits";
  olong outDisk  = 1;
  gchar *outFile = "!UVSubTestOutC.uvtab";
  olong user=103;
  /* ObitIOType type = OBIT_IO_FITS;  FOR FITS */
  /* ObitIOType type = OBIT_IO_AIPS;  FOR AIPS */
  /* olong Adisk=1, Acno;
     gchar Aname[13] = {"WideField20 "};
     gchar Aclass[7] = {"uvdata"};
     gchar Atype[3] = {"UV"};
     olong  Aseq = 2; 
     olong  masterDisk = 1;
     gchar *masterFile  = "CleanVisMaster.fits";
     olong i;
     ofloat bmaj, bmin, bpa; */
  /* Control */
  gboolean Tr=TRUE, Fl=FALSE;
  olong  ccVer[1];
  ofloat antSize;

  /* Initialize Obit */
  err = newObitErr();
  ierr = 0;
  mySystem = ObitSystemStartup ("testUVSub", 1, user, 1, AIPSdir, 1, FITSdir, 
				(oint)Tr, (oint)Fl, err);
  ObitErrLog(err); /* show any error messages on err */

  /* setup input */
  uvdata = newObitUV("input FITS UV data");
  /* FOR FITS  */
  ObitUVSetFITS(uvdata, 1000, inDisk,inFile, err);

  /* setup output */
  outdata = newObitUV("output FITS UV data");
  ObitUVSetFITS(outdata, 1000, outDisk, outFile, err);

  /* setup model */
  model = newObitImage("Model image");
  ObitImageSetFITS(model, OBIT_IO_byPlane, inDisk, inModel, blc, trc, err);

  /* Make mosaic object */
  mosaic = newObitImageMosaic ("Mosaic", 1);

  /* Attach Image */
  ObitImageMosaicSetImage (mosaic, 0, model, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR creating mosaic object");
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Create Sky Model */
  SkyModel = ObitSkyModelCreate ("Sky Model", mosaic);

 /* Set inputs on sky model */
  dim[0] = dim[1] = 1;
  ObitInfoListAlwaysPut(SkyModel->info, "PBCor",   OBIT_bool,   dim, &Tr);
  /*ObitInfoListAlwaysPut(SkyModel->info, "PBCor",   OBIT_bool,   dim, &Fl);*/
  antSize = 25.0;
  ObitInfoListAlwaysPut(SkyModel->info, "antSize", OBIT_float,  dim, &antSize);
  ccVer[0] = 2;
  ObitInfoListAlwaysPut(SkyModel->info, "CCVer",   OBIT_long,    dim, ccVer);

  /* Subtract
     ObitSkyModelSubUV (SkyModel, uvdata, outdata, err); */
  /* Divide */
  ObitSkyModelDivUV (SkyModel, uvdata, outdata, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR subtracting");
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);

  fprintf (stdout,"Ends successfully\n");
  
  return 0;
} /* end of main */

