/* Include defining Obit structures and prototypes */
#include "ObitAll.h"

/* Example program to copy an AIPS image to a FITS file  */
int main ( int argc, char **argv )
{
  ObitSystem *mySystem;
  ObitErr *err;
  ObitImage *obitImage, *outImage;
  olong disk, cno, user;

  /* relative paths to AIPS data directories */
  gchar *AIPSdir[] = {"AIPSdata/"}; 

  /* relative path to FITS data directories */
  gchar *FITSdir[] = {"FITSdata/"};

  /* bottom left corner of image specified */ 
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};

  /* top right corner of image specified */ 
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};

  /*  Define input AIPS file */
  gchar Aname[13] = {"3C138       "};
  gchar Aclass[7] = {"PCube "};
  gchar Atype[3] = {"MA"};
  olong  Aseq = 2;

  /* Define output FITS file name - allow overwrite */
  gchar *Filename="!3C138.PCube.fits";

  /* Initialize Obit, define where disk directories are */
  err = newObitErr();  /* Create error/message object */
  user = 100;     /* User ID for AIPS */
  mySystem = ObitSystemStartup ("AIPS2FITS", 1, user, 1, AIPSdir, 1, FITSdir,  
                                (oint)TRUE, (oint)FALSE, err);
  ObitErrLog(err); /* show any error messages on err */

  /* Define image objects */
  obitImage = newObitImage("AIPS");
  outImage =  newObitImage("FITS");

  /* Associate input object with AIPS file */
  disk = 1;      /* AIPS and or FITS disk number */
  cno = ObitAIPSDirFindCNO(disk, user, Aname, Aclass, Atype, Aseq, err);
  if (cno<0) Obit_log_error(err, OBIT_Error, "Failure looking up input file");
  ObitErrLog(err);  /* display any errors */

  /* setup input */
  ObitImageSetAIPS(obitImage,OBIT_IO_byPlane,disk,cno,user,blc,trc,err);
 
   /* Associate output object with FITS file */
  ObitImageSetFITS(outImage,OBIT_IO_byPlane,disk,Filename,blc,trc,err);
    
  /* copy */
  outImage =  ObitImageCopy(obitImage, outImage, err);

  /* show any errors */
  ObitErrLog(err);

  /* Unreference image objects  - this will destroy the object but not
     the disk files */
  obitImage = ObitImageUnref(obitImage);
  outImage =  ObitImageUnref(outImage);

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  err = ObitErrUnref(err);  /* delete error/message object */
  
  return 0;
} /* end of main */
