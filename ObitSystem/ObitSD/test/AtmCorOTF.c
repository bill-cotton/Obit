/* $Id: AtmCorOTF.c,v 1.1.1.1 2004/07/19 17:04:43 bcotton Exp $                            */ 
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
#include "ObitOTFGetAtmCor.h"
#include "ObitIOOTFFITS.h"
#include "ObitSystem.h"
#include "ObitParser.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* AtmCorOTFin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Corrections for atmospheric model                                    */
/*----------------------------------------------------------------------- */
{
  oint i, ierr = 0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *inData= NULL;
  ObitInfoType type;
  ObitErr *err= NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *FITSdir[] = {"FITSdata/"};
  ObitTableOTFSoln *solnTable;
  olong nrec, disk, ndetect, j;
  odouble dtemp, tarr[1000];
  ofloat solint, tau0, minEl, RAoff, Decoff, *atemp=NULL, *trx=NULL, *calJy=NULL;
  gchar infile[128];

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("AtmCorOTF", 1, 0, 0, NULL, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error message`<s on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  myInput = AtmCorOTFin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input FITS file name */
  for (i=0; i<128; i++) infile[i] = 0;
  ObitInfoListGet(myInput, "infile", &type, dim, infile, err);

  /* Solution interval */
  ObitInfoListGet(myInput, "SOLINT", &type, dim, &dtemp, err);
  solint = dtemp / 86400.0; /* convert to days */

  /* Zenith opacity */
  ObitInfoListGet(myInput, "TAU0", &type, dim, &dtemp, err);
  tau0 = dtemp;

  /* Effective atmospheric temperature */
  ObitInfoListGet(myInput, "ATEMP", &type, dim, &tarr, err);
  ndetect = dim[0];  /* one per detector */
  atemp = g_malloc0(ndetect*sizeof(ofloat));
  for (j=0; j<ndetect; j++) atemp[j] = tarr[j];

  /* Minimum elevation */
  ObitInfoListGet(myInput, "MINEL", &type, dim, &dtemp, err);
  minEl = dtemp;

   /* RA offset */
  ObitInfoListGet(myInput, "RAoff", &type, dim, &dtemp, err);
  RAoff = dtemp;

   /* Dec offset */
  ObitInfoListGet(myInput, "Decoff", &type, dim, &dtemp, err);
  Decoff = dtemp;

 /* Receiver temperature per detector */
  ObitInfoListGet(myInput, "TRX", &type, dim, tarr, err);
  ndetect = dim[0];  /* one per detector */
  trx = g_malloc0(ndetect*sizeof(ofloat));
  for (j=0; j<ndetect; j++) trx[j] = tarr[j];

  /* Cal value in Jy per detector */
  ObitInfoListGet(myInput, "CALJY", &type, dim, tarr, err);
  ndetect = dim[0];  /* one per detector */
  calJy = g_malloc0(ndetect*sizeof(ofloat));
  for (j=0; j<ndetect; j++) calJy[j] = tarr[j];

  /* error check */
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Create ObitOTF for data */
  inData = newObitOTF("Input data");
  
  /* Define output, I/O size */
  disk = 1;
  nrec = 1000;
  ObitOTFSetFITS(inData,nrec,disk,infile,err);
  
  /* Open/close input OTF to fully instantiate */
  if ((ObitOTFOpen (inData, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input FITS file %s", infile);
  
  /* Close */
  if ((ObitOTFClose (inData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing input file");
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
  
  /* Set calibration parameters */
  dim[0] = 1;
  ObitInfoListPut(inData->info, "SOLINT",  OBIT_float, dim, (gpointer*)&solint, err);
  ObitInfoListPut(inData->info, "TAU0",    OBIT_float, dim, (gpointer*)&tau0,   err);
  ObitInfoListPut(inData->info, "MINEL",   OBIT_float, dim, (gpointer*)&minEl,  err);
  ObitInfoListPut(inData->info, "RAoff",   OBIT_float, dim, (gpointer*)&RAoff,  err);
  ObitInfoListPut(inData->info, "Decoff",  OBIT_float, dim, (gpointer*)&Decoff,  err);
  dim[0] = ndetect;
  ObitInfoListPut(inData->info, "ATEMP",   OBIT_float, dim, (gpointer*)atemp,  err);
  ObitInfoListPut(inData->info, "TRX",     OBIT_float, dim, (gpointer*)trx, err);
  ObitInfoListPut(inData->info, "CALJY",   OBIT_float, dim, (gpointer*)calJy, err);
 
   /* Determine calibration */
  solnTable = ObitOTFGetAtmCor (inData, inData, err);
  if (err->error) 
    Obit_log_error(err, OBIT_Error, "ERROR gain calibrating FITS file %s", infile);

  /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
  
  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  /* cleanup */
  myInput = ObitInfoListUnref(myInput);  /* delete input list */
  inData = ObitUnref(inData);
  if (trx) g_free(trx);
  if (atemp) g_free(atemp);
  if (calJy) g_free(calJy);
  
  return ierr;
} /* end of main */

ObitInfoList* AtmCorOTFin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="AtmCorOTF.in", *arg;
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
} /* end AtmCorOTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: AtmCorOTF -input file\n");
    fprintf(stderr, "Calculate atmospheric corrections for an Obit/OTF data file\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def AtmCorOTF.in\n");
    
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
  odouble dtemp;
  ObitInfoList *out = newObitInfoList();

  /* add parser items */
  /* input FITS file name */
  strTemp = "OTFdata.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "infile", OBIT_string, dim, strTemp, err);

  /* Solution interval 30 sec */
  dim[0] = 1;
  dtemp = 30.0;
  ObitInfoListPut (out, "SOLINT", OBIT_float, dim, &dtemp, err);

  /* Zenith Opacity */
  dim[0] = 1;
  dtemp = 0.1;
  ObitInfoListPut (out, "TAU0", OBIT_float, dim, &dtemp, err);

  /* Minimum elevation */
  dim[0] = 1;
  dtemp = 1.0;
  ObitInfoListPut (out, "MINEL", OBIT_float, dim, &dtemp, err);

  /* RA Offset in degrees, to be added */
  dim[0] = 1;
  dtemp = 0.0;
  ObitInfoListPut (out, "RAoff", OBIT_float, dim, &dtemp, err);

  /* Declination Offset in degrees, to be added */
  dim[0] = 1;
  dtemp = 0.0;
  ObitInfoListPut (out, "Decoff", OBIT_float, dim, &dtemp, err);

  return out;
} /* end defaultInputs */

