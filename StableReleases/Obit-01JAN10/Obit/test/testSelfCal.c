#include <stdio.h>
#include <stdlib.h>
#include "ObitSystem.h"
#include "ObitUV.h"
#include "ObitImage.h"
#include "ObitImageUtil.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanVis.h"
#include "ObitUVSelfCal.h"
#include "ObitAIPSDir.h"
#include "ObitMem.h"

/* program to test Self Cal functionality */
int main ( int argc, char **argv )
{
  ObitSystem *mySystem;
  ObitErr *err;
  olong ierr, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  /* gchar *AIPSdir[] = {"../AIPSdata/"};*/
  gchar *AIPSdir[] = {
    "/usr/AIPS/DATA/GOLLUM_1",
    "/usr/AIPS/DATA/GOLLUM_2",
    "/usr/AIPS/DATA/GOLLUM_3",
    "/usr/AIPS/DATA/GOLLUM_4",
    "/usr/AIPS/DATA/GOLLUM_5",
    "/usr/AIPS/DATA/GOLLUM_6",
    "/usr/AIPS/DATA/GOLLUM_7"
  };
  gchar *FITSdir[] = {"../testIt/"};
  ObitUV *uvdata=NULL, *outdata=NULL;
  ObitImage *model=NULL;
  ObitImageMosaic *mosaic=NULL;
  ObitSkyModel *SkyModel=NULL; 
  ObitUVSelfCal *selfCal=NULL; 
  ObitTableSN *SNTab=NULL;
  /* data */
  olong inDisk = 1;
  gchar *inFile  = "SelfCalTest.uvtab";
  /*gchar *inFile  = "PModelZeroUC.uvtab";*/
  gchar *inModel = "SelfCalTest.fits";
  olong outDisk  = 1;
  gchar *outFile = "!SelfCalOutC.uvtab";
  olong user=105;
  /* ObitIOType type = OBIT_IO_FITS;  FOR FITS */
  /* ObitIOType type = OBIT_IO_AIPS;  FOR AIPS */
  /* olong  Adisk=1, Acno; */
  /* gchar Aname[13] = {"G29CENTR "}; */
  /* gchar Aclass[7] = {"SPLIT"}; */
  /* gchar Atype[3] = {"UV"}; */
  /* olong  Aseq = 1;  */
  /* olong  Bdisk=1, Bcno; */
  /* gchar Bname[13] = {"G29CENTR "}; */
  /* gchar Bclass[7] = {"ICL001"}; */
  /* gchar Btype[3] = {"MA"}; */
  /* olong  Bseq = 600;  */
  olong  Adisk=5, Acno;
  gchar Aname[13] = {"J1613+3412"};
  gchar Aclass[7] = {"SPLIT"};
  gchar Atype[3] = {"UV"};
  olong  Aseq = 667; 
  olong  Bdisk=5, Bcno;
  /* AIPS Model */
  /*gchar Bname[13] = {"AIPSSCMAP"};
    gchar Bclass[7] = {"ICL001"};
    gchar Btype[3] = {"MA"};
    olong  Bseq = 1; 
  */
  /* Obit clean model */
  gchar Bname[13] = {"J1613+3412"};
  gchar Bclass[7] = {"IM0001"};
  gchar Btype[3] = {"MA"};
  olong  Bseq = 1; 
     /* olong  masterDisk = 1;
	gchar *masterFile  = "CleanVisMaster.fits";
	gint i;
	ofloat bmaj, bmin, bpa; */
  /* Control */
  gboolean Tr=TRUE, Fl=FALSE;
  olong  refant, isuba, prtlv, ccVer[1], itemp;
  ofloat antSize, solint, snrmin, uvfull[2], wtuv;

  /* Initialize Obit */
  err = newObitErr();
  ierr = 0;
  mySystem = ObitSystemStartup ("testSelfCal", 1, user, 7, AIPSdir, 1, FITSdir, 
				(oint)Tr, (oint)Fl, err);
  ObitErrLog(err); /* show any error messages on err */

  /* setup input */
  uvdata = newObitUV("input UV data");
  /* FOR FITS 
  ObitUVSetFITS(uvdata, 1000, inDisk, inFile, err); */

  /* FOR AIPS  */
  Acno = ObitAIPSDirFindCNO(Adisk, user, Aname, Aclass, Atype, Aseq, err);
  if (Acno<0) Obit_log_error(err, OBIT_Error, 
			     "Failure looking up input uvdata file");
  ObitUVSetAIPS(uvdata, 1000, Adisk, Acno, user, err); 
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* setup output scratch file */
  /* Save
  outdata = newObitUV("output UV data");
  ObitUVSetFITS(outdata, 1000, outDisk, outFile, err); */
  /* Scratch */
  outdata = newObitUVScratch(uvdata, err);

  /* setup model */
  model = newObitImage("Model image");
  /* FOR FITS 
     ObitImageSetFITS(model, OBIT_IO_byPlane, inDisk, inModel, blc, trc, err); */
  /* FOR AIPS  */
  Bcno = ObitAIPSDirFindCNO(Bdisk, user, Bname, Bclass, Btype, Bseq, err);
  if (Bcno<0) Obit_log_error(err, OBIT_Error, 
			     "Failure looking up input model file");
  ObitImageSetAIPS(model, OBIT_IO_byPlane, Bdisk, Bcno, user, blc, trc, err); 
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

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
  /*ObitInfoListAlwaysPut(SkyModel->info, "PBCor",   OBIT_bool,   dim, &Tr);*/
  ObitInfoListAlwaysPut(SkyModel->info, "PBCor",   OBIT_bool,   dim, &Fl);
  antSize = 25.0;
  ObitInfoListAlwaysPut(SkyModel->info, "antSize", OBIT_float,  dim, &antSize);
  ccVer[0] = 2;
  ccVer[0] = 3;
  ccVer[0] = 1;
  ObitInfoListAlwaysPut(SkyModel->info, "CCVer",   OBIT_long,    dim, ccVer);
  itemp = 20;  /* Number of components */
  itemp = 19;  /* Number of components */
  itemp = 22;  /* Number of components */
  ObitInfoListAlwaysPut(SkyModel->info, "EComp", OBIT_long,  dim, &itemp);

  /* Create Self Cal */
  selfCal = ObitUVSelfCalCreate ("SelfCal", SkyModel);
 
  /* Self cal parms */
  dim[0] = dim[1] = 1;
  refant = 1;
  ObitInfoListAlwaysPut(selfCal->info, "refAnt", OBIT_long,  dim, &refant);
  ObitInfoListAlwaysPut(selfCal->info, "avgPol",   OBIT_bool,   dim, &Tr);
  ObitInfoListAlwaysPut(selfCal->info, "avgIF",   OBIT_bool,   dim, &Tr);
  snrmin = 3.0;
  ObitInfoListAlwaysPut(selfCal->info, "minSNR", OBIT_float, dim, &snrmin);
  prtlv = 0;
  ObitInfoListAlwaysPut(selfCal->info, "prtLv",  OBIT_long,  dim, &prtlv);
  solint = 0.3333;
  ObitInfoListAlwaysPut(selfCal->info, "solInt", OBIT_float,  dim, &solint);
  wtuv = 0.01;
  ObitInfoListAlwaysPut(selfCal->info, "WtUV", OBIT_float, dim, &wtuv);
  dim[0] = 2;
  uvfull[0] = 1.032e+07; uvfull[1] = 1.224e+09;
  ObitInfoListAlwaysPut(selfCal->info, "UVR_Full", OBIT_float, dim, &uvfull);
  dim[0] = 4;
  ObitInfoListAlwaysPut(selfCal->info, "solType", OBIT_string,  dim, "L1  ");
  ObitInfoListAlwaysPut(selfCal->info, "solMode", OBIT_string,  dim, "P   ");
 
  /* Divide */
  ObitSkyModelDivUV (SkyModel, uvdata, outdata, err);
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;
  fprintf (stdout,"Divided model\n");

   /* Self Cal */
  SNTab = ObitUVSelfCalCal (selfCal, outdata, uvdata, err);
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* rereference phases */
  refant = 1;
  isuba = 1;
  ObitUVSelfCalRefAnt (selfCal, SNTab, isuba, &refant, err);
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) return ierr;

  /* DEBUG
  ObitMemPrint (stdout); */

  /* Cleanup */
  uvdata  = ObitUVUnref(uvdata);
  outdata = ObitUVUnref(outdata);
  model   = ObitImageUnref(model);
  mosaic  = ObitImageMosaicUnref(mosaic);
  selfCal = ObitUVSelfCalUnref(selfCal);
  SNTab   = ObitTableSNUnref(SNTab);

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);

  fprintf (stdout,"Ends successfully\n");
  
  return 0;
} /* end of main */

