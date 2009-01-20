#include <stdio.h>
#include <stdlib.h>
#include "ObitAll.h"

/* Third program to test functionality */
/* Test multisource uv data selection and calibration */
void TestTrace (ObitErr *err);
int main ( int argc, char **argv )
{
  ObitSystem *mySystem;
  ObitUV *uvdata, *outuv;
  ObitErr *err;
  olong dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *AIPSdir[] = {"AIPSdata/"};
  gchar *FITSdir[] = {"FITSdata/"};
  olong disk, cno, user, itemp;
  /* Multisource uv data */
  gchar Aname[13] = {"BC105short "};
  gchar Aclass[7] = {"AIPS"};
  gchar Atype[3] = {"UV"};
  olong  Aseq = 1;
  gchar Bname[13] = {"BC105short "};
  gchar Bclass[7] = {"Obit"};
  gchar Btype[3] = {"UV"};
  olong  Bseq = 1;
  gchar *sources[] = {"BL_LAC   ","xxx"};
  gchar Fname[13] = {"NI Table"};
  gchar Fclass[7] = {"AIPS "};
  gchar Ftype[3] = {"UV"};
  olong  Fseq = 1;
  /*gchar Gname[13] = {"NI Table"};
    gchar Gclass[7] = {"Obit  "};
    gchar Gtype[3] = {"UV"};
    olong  Gseq = 1;*/

  ObitTableSN *TableSN = NULL;
  ObitTableNI *TableNI = NULL;
  odouble off[2] = {0.0, 0.0};
  olong ver;
  oint numCoef;
  gchar *Stokes = "    ";
  gboolean doCalSelect, exist;

  /* Initialize Obit */
  err = newObitErr();
  user = 100;
  mySystem = ObitSystemStartup ("test3", 1, user, 1, AIPSdir, 1, FITSdir,  
				(oint)TRUE, (oint)FALSE, err);
  ObitErrLog(err); /* show any error messages on err */

  /*++++++++++++++ Test convert NI to SN table ++++++++++++++++++++++++++++*/

  /* uv data */
  uvdata = newObitUV("Input UV for Test NI to SN");

  /* setup input */
  disk = 1;
  cno = ObitAIPSDirFindCNO(disk, user, Fname, Fclass, Ftype, Fseq, err);
  if (cno<0) Obit_log_error(err, OBIT_Error, 
			    "Failure looking up input multisource uvdata file");

  ObitUVSetAIPS(uvdata,300,disk,cno,user,err);

  /* Open, close uvdata to fully instantiate */
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

  /* get input table */
  ver = 1;
  numCoef = 0;
  TableNI = newObitTableNIValue ("Test Input NI table", (Obit*)uvdata, &ver,
			       OBIT_IO_ReadOnly, numCoef, err);
  if ((TableNI==NULL) || (err->error)) {
    Obit_log_error(err, OBIT_Error, "ERROR getting NI table");
  };

  /* Convert NI #1 to new SN table */
  TableSN = ObitIoN2SolNTableConvert (uvdata, TableNI, TableSN, off, err);
  ObitErrLog(err);

 /* cleanup */
  uvdata  = ObitUVUnref(uvdata);
  TableNI = ObitTableNIUnref(TableNI);
  TableSN = ObitTableSNUnref(TableSN);

  /*++++++++++++++ Test calibrating uvdata ++++++++++++++++++++++++++++*/

  /* uv data */
  uvdata = newObitUV("Test calibrating UV data");

  /* setup input */
  disk = 1;
  cno = ObitAIPSDirFindCNO(disk, user, Aname, Aclass, Atype, Aseq, err);
  if (cno<0) Obit_log_error(err, OBIT_Error, 
			    "Failure looking up input multisource uvdata file");

  ObitUVSetAIPS(uvdata,300,disk,cno,user,err);

  /* Select Stokes */
  dim[0] = 4; dim[1] = 1;
  ObitInfoListPut (uvdata->info, "Stokes", OBIT_string, dim, 
		 (gpointer)Stokes, err);

  /* Open, close uvdata to fully instantiate */
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

  disk = 1;
  user = 100;
  cno = ObitAIPSDirAlloc(disk, user, Bname, Bclass, Btype, Bseq, 
			 &exist, err);
  ObitErrLog(err);
  if (cno<0) g_error ("Failure setting up output file");

 /* create output file */
  outuv = newObitUV("Calibrated test data");
  ObitUVSetAIPS(outuv,300,disk,cno,user,err);
  
  ObitErrLog(err);

  /* apply calibration/selection to uvdata */
  doCalSelect = TRUE;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (uvdata->info, "doCalSelect", OBIT_bool, dim, 
		   (gpointer)&doCalSelect, err);

  itemp = 2;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (uvdata->info, "doCalib", OBIT_long, dim, 
		   (gpointer)&itemp, err);

    itemp = 0;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (uvdata->info, "gainUse", OBIT_long, dim, 
		   (gpointer)&itemp, err);

  itemp = 1;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (uvdata->info, "doBand", OBIT_long, dim, 
		   (gpointer)&itemp, err);

  itemp = 1;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (uvdata->info, "BPVer", OBIT_long, dim, 
		   (gpointer)&itemp, err);

  itemp = -1;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (uvdata->info, "doPol", OBIT_long, dim, 
		   (gpointer)&itemp, err);

  itemp = 1;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (uvdata->info, "flagVer", OBIT_long, dim, 
		   (gpointer)&itemp, err);

  itemp = 0;
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (uvdata->info, "corrType", OBIT_long, dim, 
		   (gpointer)&itemp, err);

  dim[0] = strlen(sources[0]); dim[1] = 1;
  ObitInfoListPut (uvdata->info, "Sources", OBIT_string, dim, 
    (gpointer)sources[0], err);


  /* copy data */
  outuv =  ObitUVCopy(uvdata, outuv, err);
  ObitErrLog(err);

  /* show any errors */
  ObitErrLog(err);

 /* cleanup */
  uvdata = ObitUVUnref(uvdata);

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return 0;
} /* end of main */

