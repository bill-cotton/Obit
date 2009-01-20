#include <stdio.h>
#include <stdlib.h>
#include "ObitAll.h"

/* Second program to test functionality */
/* Test uvdata and imaging  */
void TestTrace (ObitErr *err);
int main ( int argc, char **argv )
{
  ObitSystem *mySystem;
  ObitImage *outImage, *beamImage;
  ObitUV *uvdata, *outuv;
  ObitTable *tabTmp;
  ObitTableSN *SNin, *SNout;
  ObitTableSNRow *SNrow;
  ObitErr *err;
  olong blc1[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc1[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *AIPSdir[] = {"AIPSdata/"};
  gchar *FITSdir[] = {"FITSdata/"};
  olong disk, cno, user;
  gchar Cname[13] = {"1331+305"};
  gchar Cclass[7] = {"AIPS"};
  gchar Ctype[3] = {"UV"};
  olong  Cseq = 1;
  /* gchar Ename[13] = {"0319+415"}; UV imaging test
     gchar Eclass[7] = {"UVIPol"};
     gchar Etype[3] = {"UV"};
     olong  Eseq = 1; */
  gchar Fname[13] = {"1331+305"};/* UV imaging test beam */
  gchar Fclass[7] = {"IBeam "};
  gchar Ftype[3] = {"MA"};
  olong  Fseq = 2;
  gchar Gname[13] = {"1331+305"};/* UV imaging test image */
  gchar Gclass[7] = {"IMap  "};
  gchar Gtype[3] = {"MA"};
  olong  Gseq = 2;

  gchar *Stokes = "I   ";
  olong ver;
  gboolean exist, doCalSelect;
  olong nVisPIO, nfield, nChAvg, nx[1], ny[1], nxBeam[1], nyBeam[1];
  ofloat rotate, xCells[1], yCells[1], xShift[1], yShift[1];
  olong iSNRow;

  /* Initialize Obit */
  err = newObitErr();
  user = 100;
  mySystem = ObitSystemStartup ("test2", 1, user, 1, AIPSdir, 1, FITSdir,  
				(oint)TRUE, (oint)FALSE, err);
  ObitErrLog(err); /* show any error messages on err */

  /*++++++++++++++ Test imaging uvdata ++++++++++++++++++++++++++++*/

  /* uv data */
  uvdata = newObitUV("Test Imaging UV data");

 
  /* setup input */
  disk = 1;
  cno = ObitAIPSDirFindCNO(disk, user, Cname, Cclass, Ctype, Cseq, err);
  if (cno<0) Obit_log_error(err, OBIT_Error, 
			    "Failure looking up input uvdata file");

  /* debug - use cno 14
  cno = 14; */

  ObitUVSetAIPS(uvdata,3000,disk,cno,user,err);

  /* Select Stokes */
  dim[0] = 4; dim[1] = 1;
  ObitInfoListPut (uvdata->info, "Stokes", OBIT_string, dim, 
		 (gpointer)Stokes, err);

  /* Open, close uvdata to fully instantiate */
  /* Open */
  if ((ObitUVOpen (uvdata, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) {
    Obit_log_error(err, OBIT_Error, "ERROR opening input UVdata file");
  } 
  /* Close */
  if ((ObitUVClose (uvdata, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input UVdata file");
  }
  ObitErrLog(err);

  /* create scratch file */
  outuv = newObitUVScratch(uvdata, err);
  ObitErrLog(err);

  /* apply calibration/selection to uvdata*/
  doCalSelect = TRUE;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (uvdata->info, "doCalSelect", OBIT_bool, dim, 
		   (gpointer)&doCalSelect, err);

  /* copy data */
  outuv =  ObitUVCopy(uvdata, outuv, err);
  ObitErrLog(err);

  /* Set imaging parameters */
  nVisPIO = 3000;
  nfield    = 1;
  nChAvg    = 0;
  nx[0]     = 128;
  nxBeam[0] = 128;
  ny[0]     = 128;
  nyBeam[0] = 128;
  xCells[0] = -0.010 / 3600.0;
  yCells[0] = +0.010 / 3600.0;
  xShift[0] = 0.0;
  yShift[0] = 0.0;
  ObitImageUtilSet(outuv,nVisPIO,nChAvg,rotate,nfield,nx,nxBeam,ny,nyBeam,
		   xCells,yCells,xShift,yShift,err);

  /* Output image - create from outuv */
  outImage =  ObitImageUtilCreateImage (outuv, 1, TRUE, err);
  ObitErrLog(err);

  disk = 1;
  user = 100;
  cno = ObitAIPSDirAlloc(disk, user, Gname, Gclass, Gtype, Gseq, 
			 &exist, err);
  if (cno<0) g_error ("Failure setting up output file");
  ObitImageSetAIPS(outImage,OBIT_IO_byPlane,disk,cno,user,blc1,trc1,err);

  /* Get Beam */
  beamImage = ObitImageRef((ObitImage*)outImage->myBeam);
  disk = 1;
  user = 100;
  cno = ObitAIPSDirAlloc(disk, user, Fname, Fclass, Ftype, Fseq, 
			 &exist, err);
  if (cno<0) g_error ("Failure setting up output file");
  ObitImageSetAIPS(beamImage,OBIT_IO_byPlane,disk,cno,user,blc1,trc1,err);

  /* make beam and image, uniform weight */ 
  ObitImageUtilMakeImage (outuv, outImage, 0, TRUE, TRUE, err);

  /* test image destruction with beam */
  outImage->myBeam = ObitImageUnref(outImage->myBeam); /*no dangling references*/
  /* beamImage = ObitImageZap (beamImage, err);*/
  /*beamImage->isScratch = TRUE;*/
 
  /* show any errors */
  ObitErrLog(err);

  /*++++++++++++++ Test SN Table I/O ++++++++++++++++++++++++++++*/
  ver = 1;
  tabTmp = newObitUVTable(uvdata, OBIT_IO_ReadOnly, "AIPS SN", &ver, err);
  if ((tabTmp==NULL) || (err->error)) {
    Obit_log_error(err, OBIT_Error, "ERROR creating SN table");
  }
  SNin = ObitTableSNConvert (tabTmp);
  tabTmp = ObitTableUnref(tabTmp);

  /* Create Table row structure */
  SNrow = newObitTableSNRow (SNin);
 
  /* output table - generate from scratch */
  ver = 2;
  SNin->numIF = 2; /* Input not opened */
  SNin->numPol = 2;
  SNout = newObitTableSNValue ("Test Output SN table", (Obit*)uvdata, &ver,
			       OBIT_IO_ReadWrite, SNin->numIF, SNin->numPol, err);
  if ((SNout==NULL) || (err->error)) {
    Obit_log_error(err, OBIT_Error, "ERROR creating SN table");
  };

  if ((ObitTableSNOpen (SNin, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input SN table table");
  
  if ((ObitTableSNOpen (SNout, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output SN table table");
  
  /* show any errors */
  ObitErrLog(err);

  /* Read a row - use internal buffer  */
  iSNRow = 9;
  if ((ObitTableSNReadRow (SNin, iSNRow, SNrow, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR reading SN Table file");
  }
  /* show any errors */
  ObitErrLog(err);
  
  /* write the chosen row */
  if (!err->error) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "SN: time %lf ant %d source %d R12 %f I22 %f ",
		   SNrow->Time,SNrow->antNo,SNrow->SourID,SNrow->Real1[1],SNrow->Imag2[1]);
  }

  /* show any errors */
  ObitErrLog(err);
  
  /* Write a row   */
  iSNRow = -1;
  if ((ObitTableSNWriteRow (SNout, iSNRow, SNrow, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR writing SN Table file");
  }
  /* show any errors */
  ObitErrLog(err);

/* Close */
  if ((ObitTableSNClose (SNin, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input SN Table file");
  }
  if ((ObitTableSNClose (SNout, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output SN Table file");
  }
  /* show any errors */
  ObitErrLog(err);

 /* cleanup */
  outImage = ObitImageUnref(outImage);
  beamImage = ObitImageUnref(beamImage);
  uvdata = ObitUVUnref(uvdata);
  SNin = ObitTableSNUnref(SNin);
  SNout = ObitTableSNUnref(SNout);

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return 0;
} /* end of main */

/* routine to test ObitErr traceback and return */
void TestTrace (ObitErr *err)
{
  /* test logging */
  Obit_log_error(err, OBIT_MildError,"Simulated error condition");

  /* test Obit_return_if_fail */
  Obit_return_if_fail ((1==2), err, "testing");

  /* test traceback */
  Obit_traceback_msg (err, "xxx", "Issac"); 
}
