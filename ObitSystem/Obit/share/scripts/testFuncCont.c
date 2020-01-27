#include "Obit.h"
#include "ObitData.h"
#include "ObitImage.h"
#include "ObitXML.h"
#include "ObitRPC.h"
#include "ObitMultiProc.h"
#include "ObitThread.h"
#include "ObitSystem.h"

/* Local version of function */
ObitInfoList* localFunc (ObitInfoList *myInputs, ObitErr *err);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   test Obit Function container - send RPC requests and execute         */
/*----------------------------------------------------------------------- */
{
  ObitErr *err = newObitErr();
  ObitSystem   *mySystem= NULL;
  ObitImage    *inData = NULL;
  ObitInfoList *taskList=NULL, *reply=NULL, *argList=NULL;
  gint32 dim[MAXINFOELEMDIM]={1,1,1,1,1};
  olong i, iarr[5] = {1,2,3,4,5};
  olong njob=8, nstream=2, ports[2] = {8770, 8771};
  gchar *Aname="DemoTest", *Aclass="I", *DataType="FITS";
  gchar *Adir[]={"/export/ssd/bcotton/SSD"};
  gchar *Fdir[]={"./"};
  gchar *pgmName="testFuncCont";
  olong pgmNumber = 1;
  olong Fdisk=0, Adisk=1, Auser=100, Aseq=1, nAIPS=1, nFITS=1;
  gchar *Fname="aaaSomeFile.fits";
  gchar *argument="Call arguments";
  ObitMultiProc *MP=NULL;
  ports[0] = 8770; ports[1] = 8771;
  gchar *URLs[] = {"http://localhost:8770/RPC2",
		   "http://localhost:8771/RPC2",};

  /* Number of parallel jobs */
  njob  = 2;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, Auser, nAIPS, Adir, 
				nFITS, Fdir, (oint)TRUE, (oint)FALSE, err);
  if (err->error) goto done;

  /* Create "task" list */
  taskList = newObitInfoList();
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  dim[0] = strlen(URLs[0]); dim[1] = njob;
  ObitInfoListAlwaysPut (taskList, "RPCURL", OBIT_string, dim, URLs);
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  dim[0] = njob;
  ObitInfoListAlwaysPut (taskList, "RPCport", OBIT_long, dim, ports);
  dim[0] = 1;
  ObitInfoListAlwaysPut (taskList, "nThreads", OBIT_long, dim, &nstream);
  /* Test file info */
  dim[0] = 1;
  if (!strcmp(DataType,"AIPS"))
    ObitInfoListAlwaysPut (taskList, "inDisk",   OBIT_long, dim, &Adisk);
  else
    ObitInfoListAlwaysPut (taskList, "inDisk",   OBIT_long, dim, &Fdisk);
  ObitInfoListAlwaysPut (taskList, "inSeq",    OBIT_long, dim, &Aseq);
  ObitInfoListAlwaysPut (taskList, "AIPSuser", OBIT_long, dim, &Auser);
  ObitInfoListAlwaysPut (taskList, "nAIPS", OBIT_long, dim, &nAIPS);
  dim[0] = strlen(Fname);
  ObitInfoListAlwaysPut (taskList, "inFile", OBIT_string, dim, Fname);
  dim[0] = strlen(Fdir[Fdisk]);
  ObitInfoListAlwaysPut (taskList, "inDir", OBIT_string, dim, Fdir[Fdisk]);
  dim[0] = strlen(Aname);
  ObitInfoListAlwaysPut (taskList, "inName", OBIT_string, dim, Aname);
  dim[0] = strlen(Aclass);
  ObitInfoListAlwaysPut (taskList, "inClass", OBIT_string, dim, Aclass);
  dim[0] = strlen(Adir[0]);
  ObitInfoListAlwaysPut (taskList, "AIPSdirs", OBIT_string, dim, Adir[0]);
  dim[0] = strlen(DataType);
  ObitInfoListAlwaysPut (taskList, "inDataType", OBIT_string, dim, DataType);
  ObitInfoListAlwaysPut (taskList, "outDType", OBIT_string, dim, DataType);

  /* Initialize threading/multiprocessing */
  ObitThreadInit(taskList);
  ObitMultiProcStart(taskList, err);
  if (err->error) goto done;

  /* Make test Image object */
  inData = ObitImageFromFileInfo("in", taskList, err);
  if (err->error) goto done;

  /* Create MultiProc */
  MP = ObitMultiProcCreate("Multi", njob, "test", (ObitMultiProcFunc)localFunc, err);
   if (err->error) goto done;

  /* Create argument */
  argList = newObitInfoList();
  dim[0] = strlen(argument);
  ObitInfoListAlwaysPut (argList, "Status", OBIT_string, dim, argument);
  dim[0] = 5;
  ObitInfoListAlwaysPut (argList, "iarray", OBIT_int, dim, iarr);
  if (err->error) goto done;

  /* attach test image as "out" */
  ObitDataGetFileInfo ((ObitData*)inData, "out", argList, err);
  if (err->error) goto done;

  /* Set argument */
  for (i=0; i<njob; i++) {
    dim[0] = 1;
    ObitInfoListAlwaysPut (argList, "job", OBIT_long, dim, &i);
    ObitMultiProcSetFuncArg (MP, i, argList);
 }

  /* Do call */
  ObitMultiProcExecute (MP, 0.0, err);
  if (err->error) goto done;

  /* Show results */
   for (i=0; i<njob; i++) {
     reply = ObitMultiProcGetFuncRet(MP, i);
     if (err->error) goto done;
     ObitInfoListPrint (reply, stdout);
     ObitInfoListUnref(reply);
   }

done:
  ObitErrLog(err); /* show any error messages on err */

  /* Shutdown multiprocessing remote procedures
  ObitMultiProcShutdown(err); */
  ObitErrLog(err); /* show any error messages on err */
  return 0;
} /* end main */

/**
 * Local test function
 * \param myInput input list
 * \param err     error structure
 * \return output list
 */
ObitInfoList* localFunc (ObitInfoList *myInput, ObitErr *err)
{  
  ObitInfoList *outList=NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *OK = "OK", *failed="Failed";
  ObitImage *testImg=NULL;
  olong plane[] = {1,1,1,1,1};
  ofloat rms;
  gchar *routine = "test";

   /* Initialize Return results */
  outList = newObitInfoList();
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  dim[0] = strlen(failed);
  ObitInfoListAlwaysPut (outList, "Status", OBIT_string, dim, failed);

  /* error? */
  if (err->error) goto finish;

  /* Do operation here */
  /* DEBUG */
  ObitInfoListPrint (myInput, stdout);

  /* Test Logging */
  Obit_log_error(err, OBIT_InfoErr, "%s: Started in localhost",routine);
  Obit_log_error(err, OBIT_InfoErr, "%s: some logging message",routine);
  Obit_log_error(err, OBIT_InfoErr, "%s: another logging message",routine);
  Obit_log_error(err, OBIT_InfoErr, "%s: yet another logging message",routine);
  /*Obit_log_error(err, OBIT_Error,   "%s: Give a test error message",routine);*/
  /* Read Image, give RMS */
  testImg = ObitImageFromFileInfo("out", myInput, err);
  ObitImageGetPlane(testImg, NULL, plane, err);
  if (err->error) Obit_traceback_val (err, routine, routine, outList);
  rms = ObitFArrayRMS(testImg->image);
  testImg = ObitImageUnref(testImg);
  Obit_log_error(err, OBIT_InfoErr, "%s: Image RMS %f",routine,rms);
  Obit_log_error(err, OBIT_InfoErr, "%s: Finished",routine);

  /* Check if OK */
  if (!err->error) {
    dim[0] = strlen(OK);
    ObitInfoListAlwaysPut (outList, "Status", OBIT_string, dim, OK);
  }

  /* Done */
 finish:
  ObitErrLog(err); /* show any error messages on err */
  
  return outList;
} /* end test */
