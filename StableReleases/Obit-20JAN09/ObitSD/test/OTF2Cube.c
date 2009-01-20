/* Convert a Penn Array OTF data set into an image cube */
/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004                                               */
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

#include "ObitOTF.h"
#include "ObitImage.h"
#include "ObitSystem.h"
#include "ObitParser.h"

/* internal prototypes */
/* Get inputs */
void OTF2Cubein (int argc, char **argv, 
		 gchar **input_file, gchar **output_file, ObitErr *err);
/* Give basic usage on error */
void Usage(void);

/* do conversion */
void OTF2Cube (ObitOTF *inOTF, ObitImage *outImage,  ObitErr *err);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Concert a Penn Array OTF into a series of Frames in qan image        */
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *inData= NULL;
  ObitImage *outImage=NULL;
  ObitErr *err= NULL;
  gchar *FITSdir[] = {"FITSdata/"};
  olong nrec, disk;
  gchar *infile, *outfile;
  olong blc[IM_MAXDIM] = {1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0};

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("OTF2Cube", 1, 0, 0, NULL, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  OTF2Cubein (argc, argv, &infile, &outfile, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Create input ObitOTF for data */
  inData = newObitOTF("Input data");
  
  /* Define input, I/O size */
  disk = 1;
  nrec = 100;
  ObitOTFSetFITS(inData,nrec,disk,infile,err);

  /* Create basic output Image */
  outImage = newObitImage("Output Image");
 
 /* Define output, I/O size */
  disk = 1;
  ObitImageSetFITS(outImage,OBIT_IO_byPlane,disk,outfile,blc,trc,err);

  /* Convert */
  OTF2Cube (inData, outImage, err);
  
   /* show any errors */
  ObitErrLog(err);

  /* cleanup */
  myInput = ObitInfoListUnref(myInput);  /* delete input list */
  inData  = ObitUnref(inData);
  outImage = ObitUnref(outImage);

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
 
  return ierr;
} /* end of main */

void OTF2Cubein (int argc, char **argv, 
		 gchar **input_file, gchar **output_file, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Parse control info from command line                                  */
/*   Input:                                                               */
/*      argc   Number of arguments from command line                      */
/*      argv   Array of strings from command line                         */
/*   Output:                                                              */
/*      input_file  Input OTF file name                                   */
/*      output_file Output OTF file name                                  */
/*      err         Obit Error stack                                      */
/*----------------------------------------------------------------------- */
{
  olong ax;
  gchar *arg;

  /* defaults */
  *input_file ="OTF2Cube.input";
  *output_file="OTF2Cube.out";

  /* command line arguments */
  if (argc<=1) Usage(); /* must have arguments */
  /* parse command line */
  for (ax=1; ax<argc; ax++) {
    arg = argv[ax];
    if (strcmp(arg, "-input") == 0){ /* input OTF file */
      *input_file = argv[++ax];
    } else if (strcmp(arg, "-output") == 0){ /* output Image file */
      *output_file = argv[++ax];
    } else { /* unknown argument */
      Usage();
    }
  }
  
} /* end OTF2Cubein */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: OTF2Cube -input file\n");
    fprintf(stderr, "Convert an Obit Penn Array OTF data file to image cube\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input OTF file, def OTF2Cube.input\n");
    fprintf(stderr, "  -output input Image file, def OTF2Cube.out\n");
    
    /*/exit(1);  bail out */
  }/* end Usage */


 void OTF2Cube (ObitOTF *inData, ObitImage *outImage, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Convert a Penn Array OTF to an image cube                             */
/*  Each row of the OTF is converted into an image plane                  */
/*   Input:                                                               */
/*      inOTF    OTF to convert                                           */
/*      outImage Output Image, created but not fullt instantiated         */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/*----------------------------------------------------------------------- */
{
  ObitIOCode retCode, oretCode;
  ofloat *array;
  ofloat *rec;
  olong ibuf, i;
  gchar *routine = "OTF2Cube";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inData));
  g_assert (ObitImageIsA(outImage));

  /* Open input OTF to fully instantiate */
  if ((ObitOTFOpen (inData, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_traceback_msg (err, routine, inData->name);
  
  /* Fill in descriptor - know this is Penn Array data = 8x8 */
  outImage->myDesc->naxis  = 3;
  outImage->myDesc->bitpix = -32;   /* Floating */

  strncpy (outImage->myDesc->ctype[0], "AZIMUTH ", IMLEN_KEYWORD-1);
  strncpy (outImage->myDesc->ctype[1], "ELEVATIO", IMLEN_KEYWORD-1);
  strncpy (outImage->myDesc->ctype[2], "FRAME   ", IMLEN_KEYWORD-1);
  outImage->myDesc->crval[0]  = 0.0;
  outImage->myDesc->crval[1]  = 0.0;
  outImage->myDesc->crval[2]  = 1.0;
  outImage->myDesc->crota[0]  = 0.0;
  outImage->myDesc->crota[1]  = 0.0;
  outImage->myDesc->crota[2]  = 0.0;
  outImage->myDesc->crpix[0]  = 5.0;
  outImage->myDesc->crpix[1]  = 5.0;
  outImage->myDesc->crpix[2]  = 1.0;
  outImage->myDesc->cdelt[0]  = 4.0 / 3600.0;
  outImage->myDesc->cdelt[1]  = 4.0 / 3600.0;
  outImage->myDesc->cdelt[2]  = 1.0;
  outImage->myDesc->inaxes[0] = 8;
  outImage->myDesc->inaxes[1] = 8;
  outImage->myDesc->inaxes[2] = inData->myDesc->nrecord;
  outImage->myDesc->xshift    = 0.0;
  outImage->myDesc->yshift    = 0.0;

  /* Open image */
  if ((ObitImageOpen (outImage, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_traceback_msg (err, routine, outImage->name);

  /* get data pointer */
  array = outImage->image->array;
 
  /* Loop over file converting rows to image frames */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {
    /* read buffer */
    retCode = ObitOTFRead (inData, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    if (retCode==OBIT_IO_EOF) break; /* done? */

    /* Record pointer */
    rec = inData->buffer;

    /* Loop over buffer */
    for (ibuf=0; ibuf<inData->myDesc->numRecBuff; ibuf++) {

      /* Convert to image plane */
      for (i=0; i<64; i++) array[i] = rec[inData->myDesc->ilocdata+i];

      /* write Image Plane */
      oretCode = ObitImageWrite (outImage, outImage->image->array, err);
      if (err->error) Obit_traceback_msg (err, routine, outImage->name);

      rec += inData->myDesc->lrec; /* Update data record pointer */
    } /* end loop over buffer */ 
  } /* end loop over OTF */
  
  /* Close */
  if ((ObitOTFClose (inData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_traceback_msg (err, routine, inData->name);
  if ((ObitImageClose (outImage, err) != OBIT_IO_OK) || (err->error>0))
    Obit_traceback_msg (err, routine, outImage->name);

  
} /* end OTF2Cube */
