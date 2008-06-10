/* $Id$  */
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

/* Program globals */
gchar *pgmName = "testHGeom";       /* Program name */
gchar *infile  = "testHGeom.inp";   /* File with program inputs */
gchar *outfile = "testHGeom.out";   /* File to contain program outputs */
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
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  /*ObitInfoType type;*/
  ObitErr      *err= NULL;
  /*gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};;*/
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ofloat       mean, rms;
  olong  inDisk = 1;
  gchar *inFile   = "HGeomTestIn.fits";
  olong  tmplDisk = 1;
  gchar *tmplFile = "HGeomTestTmpl.fits";
  /*gint  masterDisk = 1;;*/
  /*gchar *masterFile  = "HGeomTestMaster.fits";;*/
  olong outDisk = 1;
  gchar *outFile  = "!HGeomTestOut.fits";
  ObitImage *inImage=NULL, *tmplImage=NULL, *outImage=NULL, *tmpImage=NULL;
  olong inPlane[] = {1,1,1,1,1};
  olong outPlane[] = {1,1,1,1,1};
  olong hwidth=3;

  /* Initialize Obit */
  err = newObitErr();
  ierr = 0;
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); return ierr;}

  /* Set data */
  inImage = newObitImage("input Image");
  ObitImageSetFITS (inImage, OBIT_IO_byPlane, inDisk, inFile, blc, trc, err);
  tmplImage = newObitImage("Template image");
  ObitImageSetFITS (tmplImage, OBIT_IO_byPlane, tmplDisk, tmplFile, blc, trc, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); return ierr;}

  /* Generate scratch file from tmplFile*/
  tmpImage  = newObitImageScratch(tmplImage, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); return ierr;}

  /*  ObitImageOpen(tmpImage, OBIT_IO_WriteOnly, err);
      if (err->error) {ierr = 1;  ObitErrLog(err); return ierr;}*/

  /* Interpolate */
  ObitImageUtilInterpolateImage(inImage, tmpImage, inPlane, outPlane, hwidth, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); return ierr;}

  /* Quantize to output */
  outImage = ObitImageUtilQuanFITS (tmpImage, outFile, outDisk, 0.25, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); return ierr;}

  /* Open and read Image, image on member image, an ObitFArray */
  ObitImageOpen (outImage, OBIT_IO_ReadOnly, err);
  ObitImageRead (outImage, NULL, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); return ierr;}

  /* Get statistics from inImage FArray */
  mean = ObitFArrayMean(outImage->image);
  rms = ObitFArrayRMS(outImage->image);

  ObitImageClose (outImage, err); /* Close */
  if (err->error) {ierr = 1;  ObitErrLog(err); return ierr;}

  /* Tell results */
   Obit_log_error(err, OBIT_InfoErr, 
		  "%s: mean %f RMS %f", pgmName, mean, rms);

   /* show any messages and errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
  
/* cleanup */
  inImage     = ObitUnref(inImage);
  outImage    = ObitUnref(outImage);
  tmplImage   = ObitUnref(tmplImage);
  
  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

