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
  gchar *AIPSdir[] = {"../AIPSdata/"};
  gchar *FITSdir[] = {"../testIt/"};
  ObitUV *uvdata=NULL;
  ObitImage *fullField=NULL, *outImage=NULL;
  /* data */
  olong inDisk = 1;
  /*gchar *inFile = "UVImageTestIn.uvtab";*/
  /*gchar *inFile = "FLS3minBADASS.uvtab";*/
  gchar *inFile = "WideField20.uvtab";
  /*gchar *sources[] = {"C346R422"};*/
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
  gchar *outFile = "!testCleanVisOut.fits";
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
  ofloat OutlierFlux = 0.001;/* Minimum estimated outlier flux density (Jy) */
  ofloat OutlierSI   = -1.0; /* Spectral index to estimate flux density */
  olong   OutlierSize = 50;   /* Size of outlyer field */
  /* CLEAN parameters */
  olong niter = 1000;
  olong minpatch = 200;
  ofloat minflux = 0.0001;
  ofloat gain = 0.1;
  olong ccver = 1;
  gboolean Tr=TRUE, Fl=FALSE;
  gchar *name = "TstName", *class = "Class";
  olong seq=1;
  ObitDConCleanVis *Clean=NULL;

  /* Initialize Obit */
  err = newObitErr();
  ierr = 0;
  mySystem = ObitSystemStartup ("CleanVis", 1, user, 1, AIPSdir, 1, FITSdir, 
				(oint)Tr, (oint)Fl, err);
  ObitErrLog(err); /* show any error messages on err */

  uvdata = newObitUV("input FITS UV data");

  /* setup input */
  /* FOR FITS  */
  ObitUVSetFITS(uvdata ,1000, inDisk,inFile, err);

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
  ObitInfoListAlwaysPut(uvdata->info, "gainUse",     OBIT_long,    dim, &gainuse);
  ObitInfoListAlwaysPut(uvdata->info, "doBeam",      OBIT_bool,   dim, &doBeam);
  ObitInfoListAlwaysPut(uvdata->info, "imFileType",    OBIT_long,    dim, &type);
  ObitInfoListAlwaysPut(uvdata->info, "imSeq",         OBIT_long,    dim, &seq);
  ObitInfoListAlwaysPut(uvdata->info, "imDisk",        OBIT_long,    dim, &inDisk);
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
  ObitInfoListAlwaysPut(uvdata->info, "imName",        OBIT_string, dim, name);
  dim[0] = 6;
  ObitInfoListAlwaysPut(uvdata->info, "imClass",       OBIT_string, dim, class);
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

  /* Make clean object */
  Clean = ObitDConCleanVisCreate ("Clean Object", uvdata, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR creating CLEAN object");
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  dim[0] = dim[1] = 1;
  ObitInfoListAlwaysPut(Clean->info, "Niter",       OBIT_long,    dim, &niter);
  ObitInfoListAlwaysPut(Clean->info, "minFlux",     OBIT_float,  dim, &minflux);
  ObitInfoListAlwaysPut(Clean->info, "gain",        OBIT_float,  dim, &gain);
  ObitInfoListAlwaysPut(Clean->info, "minpatch",    OBIT_long,    dim, &minpatch);
  ObitInfoListAlwaysPut(Clean->info, "CCVer",       OBIT_long,    dim, &ccver);
  ObitInfoListAlwaysPut(Clean->info, "doRestore",   OBIT_bool,   dim, &Tr);
  ObitInfoListAlwaysPut(Clean->info, "doFlatten",   OBIT_bool,   dim, &Tr);
  /* bmaj = 10.0/3600.0; bmin = 3.0/3600.0; bpa = 30.0; */
  /* ObitInfoListAlwaysPut(Clean->info, "BMAJ",      OBIT_float,  dim, &bmaj); */
  /* ObitInfoListAlwaysPut(Clean->info, "BMIN",      OBIT_float,  dim, &bmin); */
  /* ObitInfoListAlwaysPut(Clean->info, "BPA",       OBIT_float,  dim, &bpa); */

  /* Set default window */
  ObitDConCleanVisDefWindow ((ObitDConClean*)Clean, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR adding windows");
 
  /* CLEAN */
  ObitDConCleanVisDeconvolve ((ObitDCon*)Clean, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR CLEANing");
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Quantize output */
  fullField = ObitImageMosaicGetFullImage (Clean->mosaic, err);
  outImage  = ObitImageUtilQuanFITS(fullField, outFile, outDisk, err);
  /* show any errors */
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return 0;
} /* end of main */

