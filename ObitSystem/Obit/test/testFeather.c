/* $Id: testFeather.c,v 1.2 2006/12/28 16:10:04 bcotton Exp $  */
/*Debugging tool            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005                                               */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
/*;  This program is free software; you can redistribute it and/or    */
/*;  modify it under the terms of the GNU General Public License as   */
/*;  published by the Free Software Foundation; either version 2 of   */
/*;  the License, or (at your option) any later version.              */
/*;                                                                   */
/*;  This program is distributed in the hope that it will be useful,  */
/*;  but WITHOUT ANY WARRANTY; without even the implied warranty of   */
/*;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    */
/*;  GNU General Public License for more details.                     */
/*;                                                                   */
/*;  You should have received a copy of the GNU General Public        */
/*;  License along with this program; if not, write to the Free       */
/*;  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,     */
/*;  MA 02139, USA.                                                   */
/*;                                                                   */
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitImage.h"
#include "ObitSystem.h"
#include "ObitAIPSDir.h"
#include "ObitImageUtil.h"
#include "ObitHistory.h"
#include "ObitFeatherUtil.h"

/* Program globals */
gchar *pgmName = "testFeather";       /* Program name */
gchar *infile  = "testFeather.inp";   /* File with program inputs */
gchar *outfile = "testFeather.out";   /* File to contain program outputs */
olong  pgmNumber=1;     /* Program number (like POPS no.) */
olong  AIPSuser=100;    /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=1;         /* Number of FITS directories */
gchar *FITSdirs[]={"../testIt"}; /* List of FITS data directories */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Template Obit program - compute mean and RMS of an image             */
/*----------------------------------------------------------------------- */
{
#define MAXINPUT 10  /* Maximum number of input images */
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  /*ObitInfoType type;*/
  ObitErr      *err= NULL;
  /*gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};*/
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ofloat       mean, rms, cmplx[2], peak, norm;
  olong  inDisk = 1;
  gchar *inFile[] = {"FeatherTestIn1.fits", "FeatherTestIn2.fits"};
  /*gchar *inFile[] = {"FeatherTestIn1.fits", "ComaAPGBT1950.fits"};*/
  olong  tmplDisk = 1;
  gchar *tmplFile = "FeatherTestTmpl.fits";
  /*gint  masterDisk = 1;*/
  /*gchar *masterFile  = "FeatherTestMaster.fits";*/
  olong outDisk = 1;
  gchar *outFile  = "!FeatherTestOut.fits";
  ObitImage *inImage[MAXINPUT], *padImage[MAXINPUT];
  ObitImage *tmplImage=NULL, *outImage=NULL, *tmpImage=NULL;
  ObitFArray *wtArray[MAXINPUT], *resultArray=NULL, *workArray2=NULL;
  ObitCArray *accArray=NULL, *workArray=NULL;
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  olong inPlane[] = {1,1,1,1,1};
  olong outPlane[] = {1,1,1,1,1};
  olong i, numImage, hwidth=3, pos[2], naxis[2], ndim=2;
  ObitFFT *FFTfor=NULL, *FFTrev=NULL;
  gchar *name=NULL, tname[50], hicard[81];
  
  /* Initialize Obit */
  err = newObitErr();
  ierr = 0;
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); return ierr;}
  
  /* Set input data */
  numImage = 2;  /* Number of images */
  for (i=0; i<numImage; i++) {
    g_snprintf (tname, 49, "Input Image  %d", i);
    inImage[i] = newObitImage(tname);
    ObitImageSetFITS (inImage[i], OBIT_IO_byPlane, inDisk, inFile[i], blc, trc, err);
    ObitImageFullInstantiate (inImage[i], TRUE, err);
    if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  }
  
  /* Create FFTs */
  FFTfor = ObitFeatherUtilCreateFFT(inImage[0], OBIT_FFT_Forward);
  FFTrev = ObitFeatherUtilCreateFFT(inImage[0], OBIT_FFT_Reverse);

  /* Create padded images for FFT size */
  fprintf (stderr, "Pad/interpolate Images to same grid\n");

  for (i=0; i<numImage; i++) {
    fprintf (stderr, "Pad Loop %d %s\n", i+1, inImage[i]->name);
    name = g_strconcat ("Pad", inFile[i], NULL);
    padImage[i] = newObitImage(name);
    ObitImageSetFITS(padImage[i],OBIT_IO_byPlane,inDisk,name,blc,trc,err);
    /*ObitImageFullInstantiate (padImage[i], FALSE, err);*/
    g_free(name);
    if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
    /* Merely pad first image */
    if (i==0) {
      ObitFeatherUtilPad (FFTfor, inImage[i], padImage[i], err);
    } else { /* interpolate and pad rest to same grid as first */
      ObitFeatherUtilInterpol (inImage[i], padImage[0], padImage[i], err);
    }
    if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  } /* end loop creating padded images */
  
  /* SHIT DEBUG */
  /* tArr = ObitCArrayMakeF(accArray); */
  /* ObitCArrayReal (accArray, tArr); */
  /* ObitImageUtilArray2Image ("FeatherDebug.fits",1,padImage[0]->image, err); */
  /* tArr = ObitFArrayUnref(tArr); */

  /* Create masks in FArrays, first get weights from restoring beams/resolution */
  fprintf (stderr, "Create weighting masks\n");

  for (i=0; i<numImage; i++) {
    wtArray[i] =  ObitFeatherUtilMakeBeamMask (padImage[i],  FFTfor, err);
    if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  } /* end loop creating weighting masks */

  /* derive weight masks from FT of beams, Weights are 1 with a Gaussian 
     hole in the middle representing the uv coverage of the next smallest 
     array/telescope */
  for (i=0; i<numImage; i++) {
    ObitFArrayFill (wtArray[i], 1.0);  /* 1 fill */
    /* If this is not the lowest resolution image, subtract next lowest */
    if (i<numImage-1) 
      ObitFArraySub (wtArray[i], wtArray[i+1], wtArray[i]);
  }  

  /* Make accumulation array and work array */
  fprintf (stderr, "Accumulate Weighted FFTs of images\n");
  accArray  = ObitFeatherUtilCreateFFTArray(FFTfor);
  cmplx[0] = 0.0; cmplx[1] = 0.0;
  ObitCArrayFill(accArray, cmplx);  /* Zero fill accumulation */
  workArray = ObitFeatherUtilCreateFFTArray(FFTfor);

  /* Loop accumulating images */
  for (i=0; i<numImage; i++) {
    /* Read input */
    ObitImageOpen  (padImage[i], OBIT_IO_ReadOnly, err);
    ObitImageRead (padImage[i], NULL, err);
    ObitImageClose (padImage[i],err);
    if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
    ObitFeatherUtilAccumImage(FFTfor, padImage[i], wtArray[i], 
			      accArray, workArray, err);
    if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  } /* end loop accumulating images */

  /* FFT back to image domain */
  fprintf (stderr, "FFT back to image domain\n");
  
  naxis[0] = padImage[0]->myDesc->inaxes[0]; 
  naxis[1] = padImage[0]->myDesc->inaxes[1]; 
  resultArray = ObitFArrayCreate("Result Array", ndim, naxis);
  ObitFFTC2R(FFTrev, accArray, resultArray);
  ObitFArray2DCenter (resultArray);/* Swaparoonie */

  /* Get normalization by repeating but using the padded images 
     replaced by the beam */
  fprintf (stderr, "Get normalization using point models\n");
  
  cmplx[0] = 0.0; cmplx[1] = 0.0;
  ObitCArrayFill(accArray, cmplx);  /* Zero fill accumulation */
  
  /* Loop accumulating normalization images */
  for (i=0; i<numImage; i++) {
    /* replace array on padImage with model */
    /* Read input */
    ObitImageOpen (padImage[i], OBIT_IO_ReadOnly, err);
    ObitImageRead (padImage[i], NULL, err);
    ObitImageClose (padImage[i],err);
    if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
    /* Replace data with model */
    ObitFeatherUtilCreateModel(padImage[i], padImage[i]->image);
    ObitFeatherUtilAccumImage(FFTfor, padImage[i], wtArray[i], 
			      accArray, workArray, err);
    if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  } /* end loop accumulating normalization images */

  /* FFT normalization image back to image domain */
  workArray2 = newObitFArray("Scratch Array");
  ObitFArrayClone(resultArray, workArray2, err);
  if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  ObitFFTC2R(FFTrev, accArray, workArray2);
  ObitFArray2DCenter (workArray2);/* Swaparoonie */
   
  /* Do normalization from peak in workArray2 */
  pos[0] = 0; pos[1] = 0;
  peak = ObitFArrayMax(workArray2, pos);
  fprintf(stderr, "peak in normalization image %f\n",peak);
  if (peak!=0.0)  norm = 1.0 / peak;
  else norm = 1.0;
  ObitFArraySMul(resultArray, norm);
  
  /* Generate scratch file from inImage[0] */
  tmpImage  = newObitImageScratch (inImage[0], err);
  ObitImageOpen (tmpImage, OBIT_IO_WriteOnly, err);  /* Open */
  if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  
  /* Extract to tmpImage from resultArray */
  ObitFeatherUtilSubImage (padImage[0], resultArray, tmpImage, err);
  if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  
  /* Do history to scratch image as table */
  inHistory  = newObitHistoryValue ("Input History", inImage[0]->info, err);
  outHistory = newObitHistoryValue ("Output History", tmpImage->info, err);
  ObitHistoryCopyHeader (inHistory, outHistory, err);
  if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  
  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  for (i=0; i<numImage; i++) {
    g_snprintf (hicard, 80, "%s input %d = %s",pgmName, i+1, inFile[i]);
    ObitHistoryWriteRec (outHistory, -1, hicard, err);
    if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  } /* end loop adding input file names */
  ObitHistoryClose (outHistory, err);
  if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
  /* Copy to quantized integer image with same geometry as inImage[0] */
  fprintf (stderr, "Write output image\n");
  outImage = ObitImageUtilQuanFITS (tmpImage, outFile, outDisk, err);
  if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;

  /* Copy history */
  inHistory  = newObitHistoryValue ("Input History", tmpImage->info, err);
  outHistory = newObitHistoryValue ("Output History", outImage->info, err);
  ObitHistoryCopy2Header (inHistory, outHistory, err);
  if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);

  /* Compare with master lie [rms diff, max abs diff, max. master]
     masterImage  = Image.newPFImage("Master image",   masterFile,  masterDisk,  True, err)
     diff = Image.PCompare(masterImage, outImage, err);
     print "Comparison, rel. max. residual",diff[1]/diff[0], 
     " rel RMS residual",diff[2]/diff[0] */

  /* Say something */
  fprintf (stderr, "Wrote feathered output to %s\n",outFile);
  
  /* Delete scratch files */
  for (i=0; i<numImage; i++) {
    ObitImageZap(padImage[i], err);
    if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) return ierr;
  } 
    
  /* cleanup */
  for (i=0; i<<numImage; i++) {
    inImage[i] = ObitImageUnref(inImage[i]);
    wtArray[i] = ObitFArrayUnref(wtArray[i]);
  }
  outImage    = ObitImageUnref(outImage);
  tmplImage   = ObitImageUnref(tmplImage);
  tmpImage    = ObitImageUnref(tmpImage);
  accArray    = ObitCArrayUnref(accArray);
  workArray   = ObitCArrayUnref(workArray);
  resultArray = ObitFArrayUnref(resultArray);
  workArray2  = ObitFArrayUnref(workArray2);

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;

} /* end main */
