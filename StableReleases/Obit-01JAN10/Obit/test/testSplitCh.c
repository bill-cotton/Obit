#include <stdio.h>
#include <stdlib.h>
#include "ObitSystem.h"
#include "ObitUV.h"
#include "ObitImage.h"
#include "ObitImageUtil.h"
#include "ObitImageMosaic.h"
#include "ObitUVImager.h"

/* program to test Splitting data by channel into multiple files */
int main ( int argc, char **argv )
{
  ObitSystem *mySystem;
  ObitErr *err;
  olong ierr, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *AIPSdir[] = {"../AIPSdata/"};
  gchar *FITSdir[] = {"../testIt/"};
  ObitUV *uvdata=NULL, *outdata[3];
  /* data */
  olong inDisk = 1;
  /*gchar *inFile = "UVImageTestIn.uvtab";*/
  gchar *inFile = "WideField20.uvtab";
  gchar *sources[] = {"C346R424"};
  olong user=103, itemp;
  /* ObitIOType type = OBIT_IO_FITS;  FOR FITS */
  /* ObitIOType type = OBIT_IO_AIPS;  FOR AIPS */
  /* olong Adisk=1, Acno;
     gchar Aname[13] = {"WideField20 "};
     gchar Aclass[7] = {"uvdata"};
     gchar Atype[3] = {"UV"};
     olong  Aseq = 2; */
  olong outDisk  = 1;
  gchar *outFile[3] = {"!testSplit1.uvtab","!testSplit2.uvtab","!testSplit3.uvtab"};
  /* Control */
  gchar *Stokes="I   ";
  gboolean doCalSelect = TRUE;
  olong docalib = 2;
  olong gainuse = 2;
  gboolean Tr=TRUE, Fl=FALSE;

  /* Initialize Obit */
  err = newObitErr();
  ierr = 0;
  mySystem = ObitSystemStartup ("testSplitCh", 1, user, 1, AIPSdir, 1, FITSdir, 
				(oint)Tr, (oint)Fl, err);
  ObitErrLog(err); /* show any error messages on err */

  uvdata = newObitUV("input FITS UV data");

  /* setup input */
  /* FOR FITS  */
  ObitUVSetFITS(uvdata, 1000, inDisk, inFile, err);

  /* Output files */
  outdata[0] =  newObitUV("output FITS UV data 1");
  ObitUVSetFITS(outdata[0], 1000, outDisk, outFile[0], err);
  outdata[1] =  newObitUV("output FITS UV data 1");
  ObitUVSetFITS(outdata[1], 1000, outDisk, outFile[1], err);
  outdata[2] =  newObitUV("output FITS UV data 1");
  ObitUVSetFITS(outdata[2], 1000, outDisk, outFile[2], err);

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
  itemp = 1;
  ObitInfoListAlwaysPut(uvdata->info, "BChan",     OBIT_long,    dim, &itemp);
  ObitInfoListAlwaysPut(uvdata->info, "BIF",     OBIT_long,    dim, &itemp);
  ObitInfoListAlwaysPut(uvdata->info, "EIF",     OBIT_long,    dim, &itemp);
  itemp = 7;
  ObitInfoListAlwaysPut(uvdata->info, "EChan",     OBIT_long,    dim, &itemp);
  dim[0] = 4;
  ObitInfoListAlwaysPut(uvdata->info, "Stokes",      OBIT_string, dim, Stokes);
  dim[0] = 16;
  ObitInfoListAlwaysPut(uvdata->info, "Sources",      OBIT_string, dim, sources[0]);

  /* split channels */
  ObitUVUtilSplitCh (uvdata, 3, outdata, err);
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return 0;
} /* end of main */

