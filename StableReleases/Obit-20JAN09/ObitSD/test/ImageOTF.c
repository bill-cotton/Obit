/* $Id$     */
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
#include "ObitIOOTFFITS.h"
#include "ObitSystem.h"
#include "ObitParser.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* ImageOTFin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Image OTF data                                                       */
/*----------------------------------------------------------------------- */
{
  oint i, itemp, ierr = 0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *inData= NULL;
  ObitImage  *outImage=NULL;
  ObitInfoType type;
  ObitErr *err= NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *FITSdir[] = {"FITSdata/"};
  olong nx, ny, nrec, disk, gainuse, flagver;
  gboolean doCalSelect;
  ofloat RA, Dec, xCells, yCells, minWt;
  odouble dtemp;
  gchar infile[128], outfile[128];

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("ImageOTF", 1, 0, 0, NULL, 1, FITSdir, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  myInput = ImageOTFin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input FITS file name */
  for (i=0; i<128; i++) infile[i] = 0;
  ObitInfoListGet(myInput, "infile", &type, dim, infile, err);

  /* output FITS file name */
  for (i=0; i<128; i++) outfile[i] = 0;
  ObitInfoListGet(myInput, "outfile", &type, dim, outfile, err);

  /* Center position */
  ObitInfoListGet(myInput, "RA", &type, dim, &dtemp, err);
  RA = dtemp;
  ObitInfoListGet(myInput, "Dec", &type, dim, &dtemp, err);
  Dec = dtemp;

  /* cell spacing */
  ObitInfoListGet(myInput, "xSpace", &type, dim, &dtemp, err);
  xCells = dtemp;
  ObitInfoListGet(myInput, "ySpace", &type, dim, &dtemp, err);
  yCells = dtemp;

  /* Image size */
  ObitInfoListGet(myInput, "nx", &type, dim, &itemp, err);
  nx = itemp;
  ObitInfoListGet(myInput, "ny", &type, dim, &itemp, err);
  ny = itemp;

  /* Calbration table if any, -1 => no cal */
  ObitInfoListGet(myInput, "GAINUSE", &type, dim, &itemp, err);
  gainuse = itemp;

   /* Flagging */
  ObitInfoListGet(myInput, "FLAGVER", &type, dim, &itemp, err);
  flagver = itemp;

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
  
  /* Set imaging parameters */
   ObitInfoListPut(inData->info, "RA",  OBIT_float, dim, &RA, err);
   ObitInfoListPut(inData->info, "Dec", OBIT_float, dim, &Dec, err);
   ObitInfoListPut(inData->info, "nx", OBIT_long, dim, &nx, err);
   ObitInfoListPut(inData->info, "ny", OBIT_long, dim, &ny, err);
   xCells = -fabs(xCells); /* RA goes backwards */
   ObitInfoListPut(inData->info, "xCells", OBIT_float, dim, &xCells, err);
   ObitInfoListPut(inData->info, "yCells", OBIT_float, dim, &yCells, err);
   minWt = 0.1; /* minimum sum of weights */
   ObitInfoListPut(inData->info, "minWt", OBIT_float, dim, &minWt, err);

   /* Calibration if requested */
   if (gainuse>=0) {
     /* apply calibration to data */
     itemp = 1;
     dim[0] = 1; dim[1] = 1;
     ObitInfoListPut (inData->info, "DOCALIB", OBIT_long, dim, &itemp, err);
     
     /* Most recent table */
     itemp = gainuse;
     dim[0] = 1; dim[1] = 1;
     ObitInfoListPut (inData->info, "GAINUSE", OBIT_long, dim, &itemp, err);
     
   } /* end set calibration */
  
   /* selection? */
   doCalSelect = (gainuse >= 0) || (flagver>0);
   dim[0] = 1; dim[1] = 1;
   ObitInfoListPut (inData->info, "doCalSelect", OBIT_bool, dim, &doCalSelect, err);
   
   /* flagging */
   dim[0] = 1; dim[1] = 1;
   ObitInfoListPut (inData->info, "FLAGVER", OBIT_long, dim, (gpointer)&flagver, err);

  /* Create output image object from OTF */
   outImage = ObitOTFUtilCreateImage (inData, err);
   
   /* define output */
   disk = 1;
   ObitImageSetFITS(outImage,OBIT_IO_byPlane,disk,outfile,blc,trc,err);
   
   /* Open and close to fully instantiate */
   ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
   ObitImageClose (outImage, err);
   
   /* show any errors */
   if (err->error) ierr = 1;
   ObitErrLog(err);
   if (ierr!=0) return ierr;
   
   /* Form image */
   ObitOTFUtilMakeImage (inData, outImage, TRUE, NULL, err);
   
   /* show any errors */
   if (err->error) ierr = 1;
   ObitErrLog(err);
   if (ierr!=0) return ierr;
   
   /* Shutdown Obit */
   mySystem = ObitSystemShutdown (mySystem);
   
   /* cleanup */
   myInput = ObitInfoListUnref(myInput);  /* delete input list */
   inData = ObitUnref(inData);
   
   return ierr;
} /* end of main */

ObitInfoList* ImageOTFin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="ImageOTF.in", *arg;
  gboolean init=FALSE;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint    itemp;
  ObitInfoList* list;

  /* Make default inputs InfoList */
  list = defaultInputs(err);

 /* command line arguments */
  if (argc<=1) Usage(); /* must have arguments */
  /* parse command line */
  for (ax=1; ax<argc; ax++) {
    arg = argv[ax];
    if (strcmp(arg, "-input") == 0){ /* input parameters */
      input_file = argv[++ax];
      /* parse input file */
      ObitParserParse (input_file, list, err);
      init = TRUE;

    } else if (strcmp(arg, "-outfile") == 0) { /* output image */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outfile", OBIT_string, dim, strTemp);
      
    } else if (strcmp(arg, "-GAINUSE") == 0) { /* which cal table? */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "GAINUSE", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-FLAGVER") == 0) { /* which flag table? */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "FLAGVER", OBIT_oint, dim, &itemp, err);
      
    } else { /* unknown argument */
      Usage();
    }
  } /* end parsing input arguments */
  
  /* Read defaults if no file specified */
  if (!init) ObitParserParse (input_file, list, err);

  return list;
} /* end ImageOTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: ImageOTF -input file [args]\n");
    fprintf(stderr, "Image an Obit/OTF data file\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def ImageOTF.in\n");
    fprintf(stderr, "  -outfile output OTF file, def OTFdata.fits\n");
    fprintf(stderr, "  -GAINUSE which cal table, 0=>high, -1=> none\n");
    fprintf(stderr, "  -FLAGVER which flag table?, -1=> none\n");
    
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
  oint   itemp;
  ObitInfoList *out = newObitInfoList();

  /* add parser items */
  /* output FITS file name */
  strTemp = "OTFdata.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "infile", OBIT_string, dim, strTemp, err);

  /* output FITS file name */
  strTemp = "Image.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "outfile", OBIT_string, dim, strTemp, err);

  /* RA of center in degrees */
  dim[0] = 1;
  dtemp = 0.0;
  ObitInfoListPut (out, "RA", OBIT_double, dim, &dtemp, err);

  /* Dec of center in degrees */
  dim[0] = 1;
  dtemp = 0.0;
  ObitInfoListPut (out, "Dec", OBIT_double, dim, &dtemp, err);

  /* RA cell spacing in degrees (def 1 arcmin)*/
  dim[0] = 1;
  dtemp = 0.016666667;
  ObitInfoListPut (out, "xSpace", OBIT_double, dim, &dtemp, err);

  /* Dec of center in degrees */
  dim[0] = 1;
  dtemp = 0.016666667;
  ObitInfoListPut (out, "ySpace", OBIT_double, dim, &dtemp, err);

  /* number of samples in RA */
  dim[0] = 1;
  itemp = 500;
  ObitInfoListPut (out, "nx", OBIT_oint, dim, &itemp, err);

  /* number of samples in Dec */
  dim[0] = 1;
  itemp = 500;
  ObitInfoListPut (out, "ny", OBIT_oint, dim, &itemp, err);

  /* Default no calibration */
  dim[0] = 1;
  itemp = -1;
  ObitInfoListPut (out, "GAINUSE", OBIT_oint, dim, &itemp, err);

  /* FLAG table to apply, -1 => none */
  dim[0] = 1;
  itemp = -1;
  ObitInfoListPut (out, "FLAGVER", OBIT_oint, dim, &itemp, err);

  return out;
} /* end defaultInputs */

