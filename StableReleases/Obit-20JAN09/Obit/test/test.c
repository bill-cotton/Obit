#include <stdio.h>
#include <stdlib.h>
#include "ObitAll.h"
#include "ObitTableCC.h"

/* program to test functionality */
void TestTrace (ObitErr *err);
int main ( int argc, char **argv )
{
  ObitSystem *mySystem;
  ObitImage *obitImage, *outImage;
  ObitUV *uvdata, *outuv;
  ObitTable *CCOut, *CCTab2;
  ObitTableCC *CCTab;
  ObitTableCCRow *row=NULL;
  Obit *obit, *obit2;
  ObitErr *err;
  oint nparm;
  olong blc1[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc1[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong blc2[IM_MAXDIM] = {27,24,1,1,1,1,1};
  olong trc2[IM_MAXDIM] = {59,50,2,0,0,0,0};
  olong dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  odouble data, odata;
  gconstpointer  ObitClass;
  ObitInfoType type;
  gchar *Filename="3C138VLA.PCube.fits.gz";
  /*gchar *Filename="testout.fits";*/
  gchar *Outname="!testout3.fits";
  gchar *Inname2 = "3C48.uvtab.gz";
  gchar *Outname2="!testUVout.uvfits";
  gchar *AIPSdir[] = {"AIPSdata/"};
  gchar *FITSdir[] = {"FITSdata/"};
  olong disk, cno, cnoo, user;
  gchar Aname[13] = {"3C138       "};
  gchar Aclass[7] = {"PCube "};
  gchar Atype[3] = {"MA"};
  olong  Aseq = 3;
  gchar Bname[13] = {"3C138  copy "};
  gchar Bclass[7] = {"Obit "};
  gchar Btype[3] = {"MA"};
  olong  Bseq = 2;
  gchar Cname[13] = {"1331+305"};
  gchar Cclass[7] = {"SPLIT"};
  gchar Ctype[3] = {"UV"};
  olong  Cseq = 1;
  gchar Dname[13] = {"0137+331 tst"};
  gchar Dclass[7] = {"Obit "};
  gchar Dtype[3] = {"UV"};
  olong  Dseq = 1;
  gboolean exist, compress;
  olong ver, rowno;

  /* Initialize Obit */
  err = newObitErr();
  user = 100;
  mySystem = ObitSystemStartup ("test", 1, user, 1, AIPSdir, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  ObitErrLog(err); /* show any error messages on err */

  /* test create base class */
  obit = newObit("Nobody");

  /* test Referencing */
  obit2 = ObitRef(obit);

  /* unreference the original */
  obit = ObitUnref(obit);

  /* is the copy still OK? */
  ObitClass =  ObitGetClass();
  g_assert (ObitIsA(obit2, ObitClass));

  /* Now the unreference the copy */
  obit2 = ObitUnref(obit2);

  /* test ObitErr */
  TestTrace(err);
  ObitErrLog(err);

  /*++++++++++++++++++ Test AIPS Images +++++++++++++++++++++++++++*/
  /* Copy an AIPS image to a FITS file */
  obitImage = newObitImage("AIPS");
  outImage =  newObitImage("FITS");

  /* test catalog entry lookup */
  disk = 1;
  cno = ObitAIPSDirFindCNO(disk, user, Aname, Aclass, Atype, Aseq, err);
  if (cnoo<0) Obit_log_error(err, OBIT_Error, 
			     "Failure looking up input file");
  ObitErrLog(err);
  /* setup input */
  ObitImageSetAIPS(obitImage,OBIT_IO_byPlane,disk,cno,user,blc2,trc2,err);
 
  /* Setup output */
  cnoo = ObitAIPSDirAlloc(disk, user, Bname, Bclass, Btype, Bseq, 
			  &exist, err);
  if (cnoo<0) g_error ("Failure setting up output file");
  ObitImageSetAIPS(outImage,OBIT_IO_byPlane,disk,cnoo,user,blc1,trc1,err);

  /* copy */
  outImage =  ObitImageCopy(obitImage, outImage, err);

  /* show any errors */
  ObitErrLog(err);
  ObitErrClear(err);

  /*++++++++++++++++++ Test AIPS tables +++++++++++++++++++++++++++*/
  ver = 1;
  nparm = 0;
  CCTab = newObitTableCCValue ("AIPS CC", (Obit*)obitImage, &ver, 
			       OBIT_IO_ReadWrite, nparm, err);
  if ((CCTab==NULL) || (err->error)) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR finding input AIPS table");
  }

  /* Open */
  if ((ObitTableCCOpen (CCTab, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, 
		   "ERROR opening input AIPS table");
		   }

  row =  newObitTableCCRow (CCTab);
  ObitTableCCSetRow (CCTab, row, err);

 /* Read a row - use internal buffer */
  rowno = 1;
  if ((ObitTableCCReadRow (CCTab, rowno, row, err)   
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR reading input Table file");
  }

  /* write new row */
  rowno = 10;
  ObitTableCCWriteRow (CCTab, rowno, row, err);
  if (!err->error) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "CC: %f %f %f ",
		   CCTab->buffer[CCTab->myDesc->offset[0]], 
		   CCTab->buffer[CCTab->myDesc->offset[1]], 
		   CCTab->buffer[CCTab->myDesc->offset[2]]);
  }

  /* Close */
  if ((ObitTableCCClose (CCTab, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input Table file");
  }
  ver = 1;
  CCOut = newObitImageTable(outImage, OBIT_IO_WriteOnly, "AIPS CC", &ver, err);

  /* copy
  CCOut =  ObitTableCopy(CCTab, CCOut, err); */

  /* release objects */
  CCTab     = ObitTableUnref(CCTab);
  obitImage = ObitImageUnref(obitImage);
  CCOut     = ObitTableUnref(CCOut);
  outImage  = ObitImageUnref(outImage);

  /* show any errors */
  ObitErrLog(err);
  ObitErrClear(err);

  /*++++++++++++++++++ Test FITS Images +++++++++++++++++++++++++++*/
  obitImage = newObitImage("Alphonse");
  g_message ("Hello World, I am %s",obitImage->name);

  /* test storing a value and retrieving it */
  data = 5.7890;
  ObitInfoListPut (obitImage->info, "Value", OBIT_double, dim, 
		 (gpointer)&data, err);
  if (ObitInfoListGet (obitImage->info, "Value", &type, 
		     (gint32*)&dim, (gpointer)&odata, err)) {
  Obit_log_error(err, OBIT_InfoErr, 
		 "InfoList test, wrote %f, retrieved %f",data, odata);
   }

  /* test setup macro */
  ObitImageSetFITS(obitImage,OBIT_IO_byPlane,1,Filename,blc2,trc2,err);
    
  /* test image opening */
  if ((ObitImageOpen (obitImage, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR opening input FITS file %s", Filename);
  }

  /* show any errors */
  ObitErrLog(err);
  ObitErrClear(err);

  /* copy image to output */
  /* create output */
  outImage =  newObitImage("Output");
  ObitImageSetFITS(outImage,OBIT_IO_byPlane,1,Outname,blc1,trc1,err);

  /* copy */
  outImage =  ObitImageCopy(obitImage, outImage, err);

  /* show any errors */
  ObitErrLog(err);
  ObitErrClear(err);

  /*++++++++++++++++++ Test FITS tables +++++++++++++++++++++++++++*/
  ver = 1;
  CCTab2 = newObitImageTable(obitImage, OBIT_IO_ReadOnly, "AIPS CC", 
			    &ver, err);
  if ((CCTab2==NULL) || (err->error)) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR finding input FITS table");
  }

  if ((ObitTableOpen (CCTab2, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR opening input FITS table");
  }

  /* Read a row - use internal buffer  */
  if ((ObitTableRead (CCTab2, 3, NULL, err) 
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR reading FITS Table file");
  }

  /* write the first row */
  if (!err->error) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "CC: %f %f %f ",
		   CCTab2->buffer[CCTab2->myDesc->offset[0]], 
		   CCTab2->buffer[CCTab2->myDesc->offset[1]], 
		   CCTab2->buffer[CCTab2->myDesc->offset[2]]);
  }

  /* Close */
  if ((ObitTableClose (CCTab2, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input Table file");
  }
  ver = 1;
  CCOut = newObitImageTable(outImage, OBIT_IO_WriteOnly, "AIPS CC", &ver, err);

  /* copy */
  CCOut =  ObitTableCopy(CCTab2, CCOut, err);

  /* release objects */
  obitImage = ObitImageUnref(obitImage);
  outImage  = ObitImageUnref(outImage);
  CCTab2     = ObitTableUnref(CCTab2);
  CCOut     = ObitTableUnref(CCOut);

  /* show any errors */
  ObitErrLog(err);
  ObitErrClear(err);

  /*++++++++++++++ Test AIPS uvdata ++++++++++++++++++++++++++++*/
  uvdata = newObitUV("Test AIPS UV data");

  /* setup input */
  disk = 1;
  cno = ObitAIPSDirFindCNO(disk, user, Cname, Cclass, Ctype, Cseq, err);
  if (cno<0) Obit_log_error(err, OBIT_Error, 
			    "Failure looking up input uvdata file");

  ObitUVSetAIPS(uvdata,3000,disk,cno,user,err);

  /* Open */
  if ((ObitUVOpen (uvdata, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) {
    Obit_log_error(err, OBIT_Error, "ERROR opening input UVdata file");
  } 

  /* Read a vis - use internal buffer  */
  if ((ObitUVRead (uvdata, NULL, err) 
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR reading input UVdata file");
  }

  /* Close */
  if ((ObitUVClose (uvdata, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input UVdata file");
  }

  /* write the first 5 random parameters (u,v,w,b,t) */
  if (!err->error) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "Vis: %f %f %f %f %f",
		   uvdata->buffer[0], uvdata->buffer[1], uvdata->buffer[2], 
		   uvdata->buffer[3], uvdata->buffer[4]*24.0);
  }
  
  /* Setup output */
  cnoo = ObitAIPSDirAlloc(disk, user, Dname, Dclass, Dtype, Dseq, 
			  &exist, err);
  if (cnoo<0) Obit_log_error(err, OBIT_Error, 
			     "Failure creating output AIPS uv file");

  outuv = newObitUV("Test output AIPS UV data");
  ObitUVSetAIPS(outuv,3000,disk,cnoo,user,err);

  /* Compress output */
  compress = TRUE;
  ObitInfoListPut (outuv->info, "Compress", OBIT_bool, dim, 
		 (gpointer)&compress, err);

  /* copy */
  outuv =  ObitUVCopy(uvdata, outuv, err);

  /* show any errors */
  ObitErrLog(err);
  ObitErrClear(err);

  /* delete objects */
  uvdata = ObitUVUnref(uvdata);
  outuv = ObitUVUnref(outuv);

  /*++++++++++++++ Test FITS uvdata ++++++++++++++++++++++++++++*/
  uvdata = newObitUV("Test FITS UV data");

  /* setup input */
  ObitUVSetFITS(uvdata,3000,1,Inname2,err);

  /* Open */
  if ((ObitUVOpen (uvdata, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) {
    Obit_log_error(err, OBIT_Error, "ERROR opening input UVdata file");
  } 

  /* Read a vis - use internal buffer  */
  if ((ObitUVRead (uvdata, NULL, err) 
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR reading input UVdata file");
  }

  /* Close */
  if ((ObitUVClose (uvdata, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input UVdata file");
  }

  /* write the first 5 random parameters (u,v,w,b,t) */
  if (!err->error) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "Vis: %f %f %f %f %f",
		   uvdata->buffer[0], uvdata->buffer[1], uvdata->buffer[2], 
		   uvdata->buffer[3], uvdata->buffer[4]*24.0);
  }
  
  outuv = newObitUV("Test output FITS UV data");
  ObitUVSetFITS(outuv,3000,1,Outname2,err);

  /* Compress output */
  compress = TRUE;
  ObitInfoListPut (outuv->info, "Compress", OBIT_bool, dim, 
		 (gpointer)&compress, err);

  /* copy */
  outuv =  ObitUVCopy(uvdata, outuv, err);

  /* show any errors */
  ObitErrLog(err);

  /* zap 'em */
  uvdata = ObitUnref(uvdata);
  outuv  = ObitUnref(outuv);

  /* show any errors */
  ObitErrLog(err);


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
