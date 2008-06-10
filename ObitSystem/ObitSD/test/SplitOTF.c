/* To do:
1) Either calculate a proper index table or delete the zero length one,
   this is really fouling up the indexing.  This is no longer screwing up
   but one of these fixes should be implemented.
*/
/* $Id: SplitOTF.c,v 1.1.1.1 2004/07/19 17:04:43 bcotton Exp $                            */
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
#include "ObitSystem.h"
#include "ObitParser.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* SplitOTFin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/*  */
void SplitCopyArrayGeom (ObitOTF *in, ObitOTF *out, ObitErr *err);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Apply calibration and write a new OTF                                */
/*----------------------------------------------------------------------- */
{
  oint i, ierr = 0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *inData= NULL, *outData=NULL;
  ObitInfoType type;
  ObitErr *err= NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *FITSdir[] = {"FITSdata/"};
  olong nrec, disk, itemp, gainuse, flagver;
  gboolean doCalSelect, average;
  gchar infile[128], outfile[128];

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("SplitOTF", 1, 0, 0, NULL, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  myInput = SplitOTFin (argc, argv, err);
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

  /* Calibration */
  ObitInfoListGet(myInput, "GAINUSE", &type, dim, &itemp, err);
  gainuse = itemp;

  /* Editing */
  ObitInfoListGet(myInput, "FLAGVER", &type, dim, &itemp, err);
  flagver = itemp;

  /* Averaging */
  ObitInfoListGet(myInput, "Average", &type, dim, &itemp, err);
  average = itemp != 0;

  /* error check */
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Create input ObitOTF for data */
  inData = newObitOTF("Input data");
  
  /* Define input, I/O size */
  disk = 1;
  nrec = 100;
  ObitOTFSetFITS(inData,nrec,disk,infile,err);
  
  /* Create output ObitOTF for data */
  outData = newObitOTF("Output data");
  
  /* Define output, I/O size */
  disk = 1;
  nrec = 100;
  ObitOTFSetFITS(outData,nrec,disk,outfile,err);
  
  /* Open/close input OTF to fully instantiate */
  if ((ObitOTFOpen (inData, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input FITS file %s", infile);
  
  /* Close */
  if ((ObitOTFClose (inData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing input file");
  if (err->error) ierr = 1;
  
   /* show any errors */
  ObitErrLog(err);

  /* Apply prior calibration as requested */
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (inData->info, "GAINUSE", OBIT_long, dim, (gpointer)&gainuse, err);
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (inData->info, "FLAGVER", OBIT_long, dim, (gpointer)&flagver, err);
  dim[0] = 1; dim[1] = 1;
  if (gainuse>=0) itemp = 1;
  else itemp = 0;
  ObitInfoListPut (inData->info, "DOCALIB", OBIT_long, dim, (gpointer)&itemp, err);
  doCalSelect = (gainuse >= 0) || (flagver >= 0);
  ObitInfoListPut (inData->info, "doCalSelect", OBIT_bool, dim, (gpointer)&doCalSelect, err);
  
  /* Copy image from inData to outData */
  if (average) {
    ObitOTFAver  (inData, outData, err);          /* Copy/average data */
    /*    SplitCopyArrayGeom(inData, outData, err); *//* copy/select array geometry */
  } else ObitOTFCopy(inData, outData, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR copy/calibration data");
  
  /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
  
  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  /* cleanup */
  myInput = ObitInfoListUnref(myInput);  /* delete input list */
  inData  = ObitUnref(inData);
  outData = ObitUnref(outData);
 
  return ierr;
} /* end of main */

ObitInfoList* SplitOTFin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="SplitOTF.in", *arg;
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
} /* end SplitOTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: SplitOTF -input file\n");
    fprintf(stderr, "Split an Obit/OTF data file applying calibration\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def SplitOTF.in\n");
    
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
  olong    itemp;
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

  /* Cal/Soln table to apply, -1 => none, 0 => highest numbered */
  dim[0] = 1;
  itemp = -1;
  ObitInfoListPut (out, "GAINUSE", OBIT_long, dim, &itemp, err);

  /* Flag table to apply, <=0 => none */
  dim[0] = 1;
  itemp = -1;
  ObitInfoListPut (out, "FLAGVER", OBIT_long, dim, &itemp, err);

  /* Minimum elevation */
  dim[0] = 1;
  dtemp = 1.0;
  ObitInfoListPut (out, "MINEL", OBIT_float, dim, &dtemp, err);

  /* Average Spectrum> != 0 -> average spectrum */
  dim[0] = 1;
  itemp = 0;
  ObitInfoListPut (out, "Average", OBIT_long, dim, &itemp, err);

  return out;
} /* end defaultInputs */

