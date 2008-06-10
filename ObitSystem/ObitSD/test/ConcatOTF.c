/* $Id: ConcatOTF.c,v 1.2 2004/08/23 15:26:14 bcotton Exp $                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003                                               */
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
#include "ObitOTFCal.h"
#include "ObitOTFUtil.h"
#include "ObitOTFGetSoln.h"
#include "ObitIOOTFFITS.h"
#include "ObitFITS.h"
#include "ObitSystem.h"
#include "ObitParser.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* ConcatOTFin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Concatenate OTF datasets                                             */
/* At present this does not concatenate any calibration tables            */
/*----------------------------------------------------------------------- */
{
  oint i, ierr = 0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *inOTF= NULL, *outOTF=NULL;
  ObitInfoType type;
  ObitErr *err= NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *FITSdir[] = {"FITSdata/"};
  olong nrec, disk;
  ofloat dayOff;
  gchar infile[128], outfile[128], *fullname;
  /* Don't copy Cal and Soln tables */
  ObitIOCode iretCode, oretCode;
 
  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("ConcatOTF", 1, 0, 0, NULL, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  myInput = ConcatOTFin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input OTF FITS file name */
  for (i=0; i<128; i++) infile[i] = 0;
  ObitInfoListGet(myInput, "infile", &type, dim, infile, err);

  /* output FITS file name */
  for (i=0; i<128; i++) outfile[i] = 0;
  ObitInfoListGet(myInput, "outfile", &type, dim, outfile, err);

  /* error check */
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Create input ObitOTF for data */
  inOTF = newObitOTF("Input data");
  
  /* Define input, I/O size */
  disk = 1;
  nrec = 1000;
  ObitOTFSetFITS(inOTF,nrec,disk,infile,err);
  
  /* Create output ObitOTF for data */
  outOTF = newObitOTF("Output data");
  
  /* Define output, I/O size */
  disk = 1;
  nrec = 1000;
  ObitOTFSetFITS(outOTF,nrec,disk,outfile,err);
  
  /* Open input OTF */
  if ((ObitOTFOpen (inOTF, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input FITS file %s", infile);

  /* I only do FITS */
  fullname = ObitFITSFilename (disk, outfile, err);

  /* check is output file exists, if not, copy input */
  if (!ObitFileExist (fullname, err)) {
    if (err->error>0) Obit_log_error(err, OBIT_Error, "ERROR testing file %s", fullname);

    /* Copy from inOTF to outOTF */
    ObitOTFCopy(inOTF, outOTF, err);

  } else { /* end of initialize new file */
  
    /* use same data buffer on input and output so don't assign buffer for output */
    if (outOTF->buffer) ObitIOFreeBuffer(outOTF->buffer); /* free existing */
    outOTF->buffer = inOTF->buffer;
    outOTF->bufferSize = -1;
    
    /* Open output OTF */
    if ((ObitOTFOpen (outOTF, OBIT_IO_ReadWrite, err) 
	 != OBIT_IO_OK) || (err->error>0))  /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", outfile);
    
    /* append to end of the file */
    outOTF->myDesc->firstRec = outOTF->myDesc->nrecord+1;
    
    /* show any errors */
    ObitErrLog(err);
    
    /* Day offset between the two datasets */
    dayOff = inOTF->myDesc->JDObs - outOTF->myDesc->JDObs;

    /* Copy data from inOTF to outOTF */
    iretCode = OBIT_IO_OK;
    oretCode = OBIT_IO_OK;
    while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
      iretCode = ObitOTFRead (inOTF, inOTF->buffer, err);
      if (iretCode!=OBIT_IO_OK) break;
      /* How many */
      outOTF->myDesc->numRecBuff = inOTF->myDesc->numRecBuff;
      /* Update time */
      for (i=0; i<inOTF->myDesc->numRecBuff; i++) {
	inOTF->buffer[inOTF->myDesc->iloct+i*inOTF->myDesc->lrec] += dayOff;
      }
      oretCode = ObitOTFWrite (outOTF, inOTF->buffer, err);
    }
    
    /* error? */
    if (err->error) Obit_log_error(err, OBIT_Error, "ERROR copying data");
    
    /* unset output buffer (may be multiply deallocated ;'{ ) */
    outOTF->buffer = NULL;
    outOTF->bufferSize = 0;
    
    /* Close */
    if ((ObitOTFClose (inOTF, err) != OBIT_IO_OK) || (err->error>0))
      Obit_log_error(err, OBIT_Error, "ERROR closing input file");
    if (err->error) ierr = 1;
    
    if ((ObitOTFClose (outOTF, err) != OBIT_IO_OK) || (err->error>0))
      Obit_log_error(err, OBIT_Error, "ERROR closing output file");
    if (err->error) ierr = 1;
    
  } /* End copy to existing file */
  
  /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
  
  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  /* cleanup */
  myInput = ObitInfoListUnref(myInput);  /* delete input list */
  inOTF   = ObitUnref(inOTF);
  outOTF  = ObitUnref(outOTF);
  if (fullname) g_free(fullname);
 
  return ierr;
} /* end of main */

ObitInfoList* ConcatOTFin (int argc, char **argv, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Parse control info from command line                                  */
/*   Input:                                                               */
/*      argc   Number of arguments from command line                      */
/*      argv   Array of strings from command line                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   return  ObitInfoList with defaults/parsed values                     */
/*----------------------------------------------------------------------- */
{
  olong ax;
  gchar *input_file="ConcatOTF.in", *arg;
  ObitInfoList* list;

  /* command line arguments */
  if (argc<=1) Usage(); /* must have arguments */
  /* parse command line */
  for (ax=1; ax<argc; ax++) {
    arg = argv[ax];
    if (strcmp(arg, "-input") == 0){ /* input parameters */
      input_file = argv[++ax];
    } else { /* unknown argument */
      Usage();
    }
  }
  
  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* parse input file */
  ObitParserParse (input_file, list, err);

  return list;
} /* end ConcatOTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: ConcatOTF -input file\n");
    fprintf(stderr, "Concat Obit/OTF data files \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def ConcatOTF.in\n");
    
    /*/exit(1);  bail out */
  }/* end Usage */

/*----------------------------------------------------------------------- */
/*  Create default input ObitInfoList                                     */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  ObitInfoList *out = newObitInfoList();

  /* add parser items */
  /* input OTF FITS file name */
  strTemp = "OTFdata.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "infile", OBIT_string, dim, strTemp, err);

  /* output OTF FITS file name */
  strTemp = "outputOTF.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "outfile", OBIT_string, dim, strTemp, err);

  return out;
} /* end defaultInputs */

