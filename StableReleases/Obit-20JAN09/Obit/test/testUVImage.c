#include <stdio.h>
#include <stdlib.h>
#include "ObitSystem.h"
#include "ObitUV.h"
#include "ObitImage.h"
#include "ObitImageUtil.h"
#include "ObitImageMosaic.h"
#include "ObitUVImager.h"

/* program to test UVImage functionality */
int main ( int argc, char **argv )
{
  ObitSystem *mySystem;
  ObitErr *err;
  olong ierr, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *AIPSdir[] = {"../AIPSdata/"};
  gchar *FITSdir[] = {"../testIt/"};
  ObitUV *uvdata=NULL;
  ObitImage *fullField=NULL, *outImage=NULL;
  /* data */
  olong inDisk = 1;
  /*gchar *inFile = "UVImageTestIn.uvtab";*/
  gchar *inFile = "WideField20UC.uvtab";
  gchar *sources[] = {"C346R424"};
  olong user=103;
  ObitIOType type = OBIT_IO_FITS;  /* FOR FITS */
  /* ObitIOType type = OBIT_IO_AIPS;  FOR AIPS */
  /* olong Adisk=1, Acno;
     gchar Aname[13] = {"WideField20 "};
     gchar Aclass[7] = {"uvdata"};
     gchar Atype[3] = {"UV"};
     olong  Aseq = 2; */
  olong outDisk  = 1;
  gchar *outFile = "!testUVImageOut.fits";
 /*  olong  masterDisk = 1;
     gchar *masterFile  = "CleanVisMaster.fits";
     olong i;
     ofloat bmaj, bmin, bpa; */
  /* Control */
  gchar *Stokes="I   ";
  ofloat FOV = 0.4167; /* 25/60 = 25' */
  /*gfloat FOV = 0.1; Small for testing */
  ofloat UVRange[] = {0.0,0.0};
  ofloat TimeRange[] = {0.0,10.0};
  ofloat Robust = 0.0;
  ofloat UVTaper[] = {0.0,0.0};
  gboolean doCalSelect = TRUE;
  olong docalib = 2;
  olong gainuse = 2;
  gboolean doFull = TRUE;
  gboolean doBeam = TRUE;
  /* Outlier */
  gchar *Catalog = "NVSSVZ.FIT";
  ofloat OutlierDist = 1.0;  /* Maximum distance to add outlyers (deg) */
  ofloat OutlierFlux = 0.01;/* Minimum estimated outlier flux density (Jy) */
  ofloat OutlierSI   = -1.0; /* Spectral index to estimate flux density */
  olong   OutlierSize = 50;   /* Size of outlyer field */
  gboolean Tr=TRUE, Fl=FALSE;
  gchar *name = "TstName", *class = "Class";
  olong seq=1;
  ObitUVImager *myImager=NULL;

  /* Initialize Obit */
  err = newObitErr();
  ierr = 0;
  mySystem = ObitSystemStartup ("testUVImage", 1, user, 1, AIPSdir, 1, FITSdir, 
				(oint)Tr, (oint)Fl, err);
  ObitErrLog(err); /* show any error messages on err */

  uvdata = newObitUV("input FITS UV data");

  /* setup input */
  /* FOR FITS  */
  ObitUVSetFITS(uvdata, 1000, inDisk,inFile, err);

   /* FOR AIPS
      Adisk = 1;
      Acno = ObitAIPSDirFindCNO(Adisk, user, Aname, Aclass, Atype, Aseq, err);
      if (Acno<0) Obit_log_error(err, OBIT_Error, 
      "Failure looking up input uvdata file");
      ObitUVSetAIPS(uvdata,1000,Adisk,Acno,user,err); */
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;
  

 /* Set inputs */
  dim[0] = dim[1] = 1;
  ObitInfoListAlwaysPut(uvdata->info, "doCalSelect", OBIT_bool,   dim, &doCalSelect);
  ObitInfoListAlwaysPut(uvdata->info, "doCalib",     OBIT_long,    dim, &docalib);
  ObitInfoListAlwaysPut(uvdata->info, "gainuse",     OBIT_long,    dim, &gainuse);
  ObitInfoListAlwaysPut(uvdata->info, "doBeam",      OBIT_bool,   dim, &doBeam);
  ObitInfoListAlwaysPut(uvdata->info, "imFileType",  OBIT_long,    dim, &type);
  ObitInfoListAlwaysPut(uvdata->info, "imSeq",       OBIT_long,    dim, &seq);
  ObitInfoListAlwaysPut(uvdata->info, "imDisk",      OBIT_long,    dim, &inDisk);
  ObitInfoListAlwaysPut(uvdata->info, "FOV",         OBIT_float,  dim, &FOV);
  ObitInfoListAlwaysPut(uvdata->info, "doFull",      OBIT_bool,   dim, &doFull);
  ObitInfoListAlwaysPut(uvdata->info, "OutlierDist", OBIT_float,  dim, &OutlierDist);
  ObitInfoListAlwaysPut(uvdata->info, "OutlierFlux", OBIT_float,  dim, &OutlierFlux);
  ObitInfoListAlwaysPut(uvdata->info, "OutlierSI",   OBIT_float,  dim, &OutlierSI);
  ObitInfoListAlwaysPut(uvdata->info, "OutlierSize", OBIT_long,    dim, &OutlierSize);
  ObitInfoListAlwaysPut(uvdata->info, "Robust",      OBIT_float,  dim, &Robust);
  dim[0] = 4;
  ObitInfoListAlwaysPut(uvdata->info, "Stokes",      OBIT_string, dim, Stokes);
  dim[0] = 12;
  ObitInfoListAlwaysPut(uvdata->info, "imName",      OBIT_string, dim, name);
  dim[0] = 6;
  ObitInfoListAlwaysPut(uvdata->info, "imClass",     OBIT_string, dim, class);
  dim[0] = strlen(Catalog);
  ObitInfoListAlwaysPut(uvdata->info, "Catalog",     OBIT_string, dim, Catalog);
  dim[0] = 2;
  ObitInfoListAlwaysPut(uvdata->info, "UVRange",     OBIT_float,  dim, UVRange);
  dim[0] = 2;
  ObitInfoListAlwaysPut(uvdata->info, "UVTaper",     OBIT_float,  dim, UVTaper);
  dim[0] = 2;
  ObitInfoListAlwaysPut(uvdata->info, "timeRange",   OBIT_float,  dim, TimeRange);
  dim[0] = 16;
  ObitInfoListAlwaysPut(uvdata->info, "Sources",      OBIT_string, dim, sources[0]);

  /* Make uvImager object */
  myImager = ObitUVImagerCreate ("Imager Object", uvdata, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR creating Imager object");
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Weight */
  ObitUVImagerWeight (myImager, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR Weighting");
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Image */
  ObitUVImagerImage (myImager, 0, Fl, Tr, Fl, err);
  /*ObitUVImagerImage (myImager, 0, Tr, Tr, Fl, err);*/
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR Imaging");
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Flatten */
  ObitUVImagerFlatten (myImager,  err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR Flattening");
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Quantize output */
  fullField = ObitImageMosaicGetFullImage (myImager->mosaic, err);
  outImage  = ObitImageUtilQuanFITS(fullField, outFile, outDisk, 0.25, err);
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Cleanup */
  myImager = ObitUVImagerUnref(myImager);

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return 0;
} /* end of main */

