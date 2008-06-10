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
#include "ObitOTFCal.h"
#include "ObitOTFUtil.h"
#include "ObitOTFSoln2Cal.h"
#include "ObitIOOTFFITS.h"
#include "ObitSystem.h"
#include "ObitParser.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* Soln2CalOTFin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Apply a Soln table to a Cal table and write a new Cal table          */
/*----------------------------------------------------------------------- */
{
  oint i, itemp, ierr = 0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *inData= NULL;
  ObitInfoType type;
  ObitErr *err= NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *FITSdir[] = {"FITSdata/"};
  ObitTableOTFCal *calTable;
  olong nrec, disk, solnuse, calin, calout;
  gchar infile[128];

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("Soln2CalOTF", 1, 0, 0, NULL, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  myInput = Soln2CalOTFin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input FITS file name */
  for (i=0; i<128; i++) infile[i] = 0;
  ObitInfoListGet(myInput, "infile", &type, dim, infile, err);

  /* Input Solution table version */
  ObitInfoListGet(myInput, "SOLNUSE", &type, dim, &itemp, err);
  solnuse = itemp;

  /* Input Cal table version */
  ObitInfoListGet(myInput, "CALIN", &type, dim, &itemp, err);
  calin = itemp;

  /* Output Calibration table version */
  ObitInfoListGet(myInput, "CALOUT", &type, dim, &itemp, err);
  calout = itemp;

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
  
  /* Specify desired tables */
  dim[0] = 1;
  ObitInfoListPut(inData->info, "SOLNUSE",  OBIT_long, dim, (gpointer*)&solnuse, err);
  ObitInfoListPut(inData->info, "CALIN",    OBIT_long, dim, (gpointer*)&calin,   err);
  ObitInfoListPut(inData->info, "CALOUT",   OBIT_long, dim, (gpointer*)&calout,  err);
 
   /* Update calibration */
  calTable = ObitOTFSoln2Cal (inData, inData, err);
  if (err->error) 
    Obit_log_error(err, OBIT_Error, "ERROR calibrating FITS file %s", infile);

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

ObitInfoList* Soln2CalOTFin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="Soln2Cal.in", *arg;
  gboolean init=FALSE;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
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

    } else if (strcmp(arg, "-SOLNUSE") == 0) { /* which soln table? */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "SOLNUSE", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-CALIN") == 0) { /* which input cal table? */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "CALIN", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-CALOUT") == 0) { /* which output cal table? */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "CALOUT", OBIT_oint, dim, &itemp, err);
      
    } else { /* unknown argument */
      Usage();
    }
  } /* end parsing input arguments */
  
  /* Read defaults if no file specified */
  if (!init) ObitParserParse (input_file, list, err);

  return list;
} /* end Soln2CalOTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: Soln2Cal -input file [-args]\n");
    fprintf(stderr, "Apply a Soln table to a Cal Table for an Obit/OTF data file\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def Soln2Cal.in\n");
    fprintf(stderr, "  -SOLNUSE  Which Soln table to apply (0=>high)\n");
    fprintf(stderr, "  -CALIN    Which input Cal table (-1=none)\n");
    fprintf(stderr, "  -CALOUT   Which output Cal table (0=>high)\n");
    
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
  oint    itemp;
  ObitInfoList *out = newObitInfoList();

  /* add parser items */
  /* input FITS file name */
  strTemp = "OTFdata.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "infile", OBIT_string, dim, strTemp, err);

  /* Input Solution table version  */
  dim[0] = 1;
  itemp = 0;
  ObitInfoListPut (out, "SOLNUSE", OBIT_oint, dim, &itemp, err);

  /* Input Cal table version */
  dim[0] = 1;
  itemp = 0;
  ObitInfoListPut (out, "CALIN", OBIT_oint, dim, &itemp, err);

  /* Output Calibration table version  */
  dim[0] = 1;
  itemp = 0;
  ObitInfoListPut (out, "CALOUT", OBIT_oint, dim, &itemp, err);

  return out;
} /* end defaultInputs */

