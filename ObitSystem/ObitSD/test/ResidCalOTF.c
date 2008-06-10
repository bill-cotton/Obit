/* $Id$                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
ObitInfoList* ResidCalOTFin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Residual calibration from a sky model                                */
/*   Solution type depends on value of parameter calType:                 */
/*     "Gain" solve for multiplicative term from "cals" in data           */
/*        SOLINT  Solution interval in sec.                               */
/*        MINRMS  Minimum RMS residual to solution                        */
/*        MINEL   Minimum elevation deg                                   */
/*        CALJY   Noise cal value in Jy per detector                      */
/*     "Offset" Solve for additive terms from residuals to the model.     */
/*        SOLINT  Solution interval in sec.                               */
/*        MINEL   Minimum elevation deg                                   */
/*     "GainOffset" Solve both gain and offset                            */
/*        SOLINT  Solution interval in sec.                               */
/*        MINRMS  Minimum RMS residual to solution                        */
/*        MINEL   Minimum elevation deg                                   */
/*        CALJY   Noise cal value in Jy per detector                      */
/*     "Filter"  Additive terms from filteres residuals to the model.     */
/*        SOLINT  Time constant in sec., shorter timescalles filtered out */
/*        MINEL   Minimum elevation deg                                   */
/*        CLIP    Clipping level for residuals                            */
/*----------------------------------------------------------------------- */
{
  oint i, ierr = 0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *inData= NULL, *scrData=NULL;
  ObitImage  *Image=NULL;
  ObitInfoType type;
  ObitErr *err= NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *FITSdir[] = {"FITSdata/"};
  ObitTableOTFSoln *solnTable;
  olong nrec, disk, ndetect, j, gainuse, flagver;
  oint itemp;
  odouble dtemp, tarr[1000];
  gboolean doCalSelect;
  ofloat clip, solint, minrms, minEl, minFlux, *calJy=NULL;
  gchar infile[128], imagefile[128], calType[50];

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("ResidCalOTF", 1, 0, 0, NULL, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  myInput = ResidCalOTFin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input OTF FITS file name */
  for (i=0; i<128; i++) infile[i] = 0;
  ObitInfoListGet(myInput, "infile", &type, dim, infile, err);

  /* image FITS file name */
  for (i=0; i<128; i++) imagefile[i] = 0;
  ObitInfoListGet(myInput, "image", &type, dim, imagefile, err);

  /* calibration type */
  for (i=0; i<50; i++) calType[i] = 0;
  ObitInfoListGet(myInput, "calType", &type, dim, calType, err);

  /* Solution interval  */
  ObitInfoListGet(myInput, "SOLINT", &type, dim, &dtemp, err);
  solint = dtemp / 86400.0; /* convert to days */

  /* Minimum elevation (deg) */
  ObitInfoListGet(myInput, "MINEL", &type, dim, &dtemp, err);
  minEl = dtemp;

  /* Clipping level */
  ObitInfoListGet(myInput, "CLIP", &type, dim, &dtemp, err);
  clip = dtemp;

  /* Minimum brightness in model */
  ObitInfoListGet(myInput, "MINFLUX", &type, dim, &dtemp, err);
  minFlux = dtemp;

  /* Prior Calibration */
  ObitInfoListGet(myInput, "GAINUSE", &type, dim, &itemp, err);
  gainuse = itemp;

  /* Flagging */
  ObitInfoListGet(myInput, "FLAGVER", &type, dim, &itemp, err);
  flagver = itemp;

  /* minimum RMS */
  ObitInfoListGet(myInput, "MINRMS", &type, dim, &dtemp, err);
  minrms = dtemp;

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
  if ((ObitOTFOpen (inData, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input FITS file %s", infile);
  
  /* Close */
  if ((ObitOTFClose (inData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing input file");
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
  
  /* Create model image  */
  Image = newObitImage("Model image");

  /* Specify image and I/O */
  ObitImageSetFITS(Image,OBIT_IO_byPlane,1,imagefile,blc,trc,err);
    
  /* Open and read image */
  ObitImageOpen (Image, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR opening image");

  if (!err->error) ObitImageRead (Image, Image->image->array, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR reading image");
  
  ObitImageClose (Image, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR closing image");

   /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
  
  /* Scratch file for residual data */
  scrData = newObitOTFScratch (inData, err);
  
  /* Apply prior calibration as requested */
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (inData->info, "FLAGVER", OBIT_long, dim, (gpointer)&flagver, err);
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (inData->info, "GAINUSE", OBIT_long, dim, (gpointer)&gainuse, err);
  dim[0] = 1; dim[1] = 1;
  if (gainuse>=0) itemp = 1;
  else itemp = 0;
  ObitInfoListPut (inData->info, "DOCALIB", OBIT_long, dim, (gpointer)&itemp, err);
  doCalSelect = (gainuse >= 0) || (flagver>0);
  ObitInfoListPut (inData->info, "doCalSelect", OBIT_bool, dim, (gpointer)&doCalSelect, err);

  /* clip image below minFlux  */
  Obit_log_error(err, OBIT_InfoErr, "Clip image below %g ", minFlux);
  ObitFArrayClip (Image->image, minFlux, 1.0e20, 0.0);

  
  /* Subtract image from inData to scrData */
  ObitOTFUtilSubImage(inData, scrData, Image->image, Image->myDesc, err);
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR subtracting image");
  
  /* Set calibration parameters */
  dim[0] = 1;
  ObitInfoListPut(scrData->info, "SOLINT",  OBIT_float, dim, (gpointer*)&solint, err);
  ObitInfoListPut(scrData->info, "MINRMS",  OBIT_float, dim, (gpointer*)&minrms, err);
  ObitInfoListPut(scrData->info, "MINEL",   OBIT_float, dim, (gpointer*)&minEl,  err);
  ObitInfoListPut(scrData->info, "CLIP",    OBIT_float, dim, (gpointer*)&clip,  err);
  dim[0] = ndetect;
  ObitInfoListPut(scrData->info, "CALJY",   OBIT_float, dim, (gpointer*)calJy, err);
  dim[0] = strlen(calType);
  ObitInfoListPut(scrData->info, "calType",   OBIT_string, dim, calType, err);
 
  /* Determine calibration by type */
  if (!strncmp(calType, "Filter",6)) { /* Time filtering */
    solnTable = ObitOTFGetSolnFilter (scrData, inData, err);
  } else { /* offset of gain */
    solnTable = ObitOTFGetSolnGain (scrData, inData, err);
  }
  if (err->error) Obit_log_error(err, OBIT_Error, "ERROR calibrating FITS file %s", infile);

  /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
  
  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  /* cleanup */
  myInput = ObitInfoListUnref(myInput);  /* delete input list */
  inData = ObitUnref(inData);
  if (calJy) g_free(calJy);
 
  return ierr;
} /* end of main */

ObitInfoList* ResidCalOTFin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="ResidCalOTF.in", *arg;
  gboolean init=FALSE;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  odouble dtemp;
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

    } else if (strcmp(arg, "-image") == 0) { /* reference image */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "image", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-calType") == 0) { /* Calibration type */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "calType", OBIT_string, dim, strTemp);
      
    } else if (strcmp(arg, "-SOLINT") == 0) { /* solution interval */
      dim[0] = 1;
      dtemp = strtod(argv[++ax], NULL);
      ObitInfoListPut (list, "SOLINT", OBIT_double, dim, &dtemp, err);
      
    } else if (strcmp(arg, "-MINEL") == 0) { /* min elevation */
      dim[0] = 1;
      dtemp = strtod(argv[++ax], NULL);
      ObitInfoListPut (list, "MINEL", OBIT_double, dim, &dtemp, err);
      
    } else if (strcmp(arg, "-CLIP") == 0) { /* residual clipping level */
      dim[0] = 1;
      dtemp = strtod(argv[++ax], NULL);
      ObitInfoListPut (list, "CLIP", OBIT_double, dim, &dtemp, err);
      
    } else if (strcmp(arg, "-MINFLUX") == 0) { /* min brightness to use in model */
      dim[0] = 1;
      dtemp = strtod(argv[++ax], NULL);
      ObitInfoListPut (list, "MINFLUX", OBIT_double, dim, &dtemp, err);
      
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
} /* end ResidCalOTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: ResidCalOTF -input file [-args]\n");
    fprintf(stderr, "ResidCal an Obit/OTF data file\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def ResidCalOTF.in\n");
    fprintf(stderr, "  -image    Name of reference image \n");
    fprintf(stderr, "  -SOLINT   Solution interval (sec) \n");
    fprintf(stderr, "  -MINFLUX  Min. brightness in model to use \n");
    fprintf(stderr, "  -MINEL    Minimum elevation deg \n");
    fprintf(stderr, "  -CLIP     Residual clipping level \n");
    fprintf(stderr, "  -GAINUSE  Which gain table to apply (-1=none)\n");
    fprintf(stderr, "  -FLAGVER  Flag table version (-1=none)\n");
   
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
  oint    itemp;
  ObitInfoList *out = newObitInfoList();

  /* add parser items */
  /* input OTF FITS file name */
  strTemp = "OTFdata.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "infile", OBIT_string, dim, strTemp, err);

  /* input image FITS file name */
  strTemp = "image.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "image", OBIT_string, dim, strTemp, err);

  /* calibration type */
  strTemp = "Offset";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "calType", OBIT_string, dim, strTemp, err);

  /* Solution interval 30 sec */
  dim[0] = 1;
  dtemp = 30.0;
  ObitInfoListPut (out, "SOLINT", OBIT_double, dim, &dtemp, err);

  /* Cal/Soln table to apply, -1 => none, 0 => highest numbered */
  dim[0] = 1;
  itemp = -1;
  ObitInfoListPut (out, "GAINUSE", OBIT_oint, dim, &itemp, err);

  /* FLAG table to apply, -1 => none */
  dim[0] = 1;
  itemp = -1;
  ObitInfoListPut (out, "FLAGVER", OBIT_oint, dim, &itemp, err);

  /* minimum fractional solution RMS */
  dim[0] = 1;
  dtemp = 0.1;
  ObitInfoListPut (out, "MINRMS", OBIT_double, dim, &dtemp, err);

  /* Minimum elevation */
  dim[0] = 1;
  dtemp = 1.0;
  ObitInfoListPut (out, "MINEL", OBIT_double, dim, &dtemp, err);

  /* Clipping level */
  dim[0] = 1;
  dtemp = 1.0e20;
  ObitInfoListPut (out, "CLIP", OBIT_double, dim, &dtemp, err);

  /* Minimum brightness in model */
  dim[0] = 1;
  dtemp = 0.0;
  ObitInfoListPut (out, "MINFLUX", OBIT_double, dim, &dtemp, err);

  return out;
} /* end defaultInputs */

