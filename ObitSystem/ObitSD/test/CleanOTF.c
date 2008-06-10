/* $Id$                            */
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
#include "ObitOTFCLEAN.h"
#include "ObitSystem.h"
#include "ObitParser.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* CleanOTFin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Clean Image created from OTF data                                    */
/*----------------------------------------------------------------------- */
{
  oint i, itemp, ierr = 0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitImage *Dirty=NULL, *Beam=NULL, *Clean=NULL;
  ObitOTFCLEAN  *myClean=NULL;
  ObitInfoType type;
  ObitErr *err= NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *FITSdir[] = {"FITSdata/"};
  olong disk, niter;
  ofloat gain, flux, beamsize;
  odouble dtemp;
  gchar dirtyfile[128], beamfile[128], cleanfile[128];

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("CleanOTF", 1, 0, 0, NULL, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  myInput = CleanOTFin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input dirty file name */
  for (i=0; i<128; i++) dirtyfile[i] = 0;
  ObitInfoListGet(myInput, "dirtyfile", &type, dim, dirtyfile, err);

  /* input beam file name */
  for (i=0; i<128; i++) beamfile[i] = 0;
  ObitInfoListGet(myInput, "beamfile", &type, dim, beamfile, err);

  /* output CLEAN file name */
  for (i=0; i<128; i++) cleanfile[i] = 0;
  ObitInfoListGet(myInput, "cleanfile", &type, dim, cleanfile, err);

  /* Number of iterations */
  ObitInfoListGet(myInput, "NITER", &type, dim, &itemp, err);
  niter = (olong)itemp;

  /* CLEAN loop gain */
  ObitInfoListGet(myInput, "GAIN", &type, dim, &dtemp, err);
  gain = (ofloat)dtemp;

  /* Minimum flux density for CLEAN */
  ObitInfoListGet(myInput, "FLUX", &type, dim, &dtemp, err);
  flux = (ofloat)dtemp;

  /* Restoring beam size in asec. */
  ObitInfoListGet(myInput, "BEAMSIZE", &type, dim, &dtemp, err);
  beamsize = (ofloat)dtemp;

  /* error check */
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Create/setup images */
  disk = 1;
  Dirty = newObitImage("Dirty Map");
  ObitImageSetFITS(Dirty,OBIT_IO_byPlane,disk,dirtyfile,blc,trc,err);
  Beam = newObitImage("Dirty Beam");
  ObitImageSetFITS(Beam,OBIT_IO_byPlane,disk,beamfile,blc,trc,err);
  Clean = newObitImage("Clean Map");
  ObitImageSetFITS(Clean,OBIT_IO_byPlane,disk,cleanfile,blc,trc,err);

  /* error check */
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;
  
  /* Create Clean object */
  myClean = newObitOTFCLEAN ("CLEAN object", OBIT_OTFCLEAN_Hogbom,
			     Dirty, Beam, Clean, err);
  /* Set CLEANing parameters */
  ObitInfoListPut(myClean->info, "NITER",    OBIT_long,   dim, &niter,    err);
  ObitInfoListPut(myClean->info, "GAIN",     OBIT_float, dim, &gain,     err);
  ObitInfoListPut(myClean->info, "FLUX",     OBIT_float, dim, &flux,     err);
  ObitInfoListPut(myClean->info, "BEAMSIZE", OBIT_float, dim, &beamsize, err);

  /* error check */
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;
  
  /* Do CLEAN */
  ObitOTFCLEANClean (myClean, err);

  /* error check */
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;
  
  /* Restore if requested */
  ObitOTFCLEANRestore (myClean, err);

  /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
   
   /* cleanup */
   myInput = ObitInfoListUnref(myInput);  /* delete input list */
   myClean = ObitUnref(myClean);
   Dirty   = ObitUnref(Dirty);
   Beam    = ObitUnref(Beam);
   Clean   = ObitUnref(Clean);
   
   /* Shutdown Obit */
   mySystem = ObitSystemShutdown (mySystem);
   
   return ierr;
} /* end of main */

ObitInfoList* CleanOTFin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="CleanOTF.in", *arg;
  gboolean init=FALSE;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint    itemp;
  odouble dtemp;
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

    } else if (strcmp(arg, "-dirtyfile") == 0) { /* dirty image */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "dirtyfile", OBIT_string, dim, strTemp);
      
    } else if (strcmp(arg, "-beamfile") == 0) { /* beam image */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "beamfile", OBIT_string, dim, strTemp);
      
    } else if (strcmp(arg, "-cleanfile") == 0) { /* clean image */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "cleanfile", OBIT_string, dim, strTemp);
      
    } else if (strcmp(arg, "-NITER") == 0) { /* How man CLEAN cycles? */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "NITER", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-GAIN") == 0) { /* CLEAN loop gain? */
      dim[0] = 1;
      dtemp = strtod(argv[++ax], NULL);
      ObitInfoListPut (list, "GAIN", OBIT_double, dim, &dtemp, err);
      
    } else if (strcmp(arg, "-FLUX") == 0) { /* Min. CLEAN flux? */
      dim[0] = 1;
      dtemp = strtod(argv[++ax], NULL);
      ObitInfoListPut (list, "FLUX", OBIT_double, dim, &dtemp, err);
      
    } else if (strcmp(arg, "-BEAMSIZE") == 0) { /* Restoring beam size? */
      dim[0] = 1;
      dtemp = strtod(argv[++ax], NULL);
      ObitInfoListPut (list, "BEAMSIZE", OBIT_double, dim, &dtemp, err);
      
    } else { /* unknown argument */
      Usage();
    }
  } /* end parsing input arguments */
  
  /* Read defaults if no file specified */
  if (!init) ObitParserParse (input_file, list, err);

  return list;
} /* end CleanOTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: CleanOTF -input file [args]\n");
    fprintf(stderr, "Do CLEAN deconvolution to an image \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file,  def CleanOTF.in\n");
    fprintf(stderr, "  -dirtyfile input dirty image, def DirtyMap.fits\n");
    fprintf(stderr, "  -beamfile  input dirty beam,  def DirtyBeam.fits\n");
    fprintf(stderr, "  -cleanfile output clean image def CleanMap.fits\n");
    fprintf(stderr, "  -NITER number of iterations, def. 200\n");
    fprintf(stderr, "  -GAIN  loop gain, def. 0.1\n");
    fprintf(stderr, "  -FLUX  min. flux density, def 0.0\n");
    fprintf(stderr, "  -BEAMSIZE restoring beam in asec, <=0 => no restore\n");
    
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
  /* Dirty image FITS file name */
  strTemp = "DirtyMap.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "dirtyfile", OBIT_string, dim, strTemp, err);

  /* Dirty Beam FITS file name */
  strTemp = "DirtyBeam.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "beamfile", OBIT_string, dim, strTemp, err);

  /* output CLEAN FITS file name */
  strTemp = "!CleanMap.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "cleanfile", OBIT_string, dim, strTemp, err);

  /* Maximum number of iterations */
  dim[0] = 1;
  itemp = 200;
  ObitInfoListPut (out, "NITER", OBIT_oint, dim, &dtemp, err);

  /* CLEAN loop gain */
  dim[0] = 1;
  dtemp = 0.1;
  ObitInfoListPut (out, "GAIN", OBIT_double, dim, &dtemp, err);

  /* Minimum flux density for CLEAN (map units) */
  dim[0] = 1;
  dtemp = 0.0;
  ObitInfoListPut (out, "FLUX", OBIT_double, dim, &dtemp, err);

  /* Restoring beam size in asec */
  dim[0] = 1;
  dtemp = 10.0;
  ObitInfoListPut (out, "BEAMSIZE", OBIT_double, dim, &dtemp, err);

  return out;
} /* end defaultInputs */

