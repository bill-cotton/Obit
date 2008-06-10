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
#include "ObitOTFGetSoln.h"
#include "ObitIOOTFFITS.h"
#include "ObitTableOTFFlag.h"
#include "ObitTableOTFTarget.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitTableOTFTargetUtil.h"
/* internal prototypes */
/* Get inputs */
ObitInfoList* EditOTFin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Make entries in OTF Flag table                                       */
/*----------------------------------------------------------------------- */
{
  oint i, itemp, iarr[2], ierr = 0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitOTF *inData= NULL;
  ObitInfoType type;
  ObitTableOTFFlag    *FlagTable=NULL;
  ObitTableOTFFlagRow *FlagRow=NULL;
  ObitTableOTFTarget  *TargetTable=NULL;
  ObitErr *err= NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *FITSdir[] = {"FITSdata/"};
  olong nrec, disk, bchan, echan, FeedFlag, TargID;
  odouble darr[2], TimeRange[2];
  gchar infile[128], stokes[5], Target[20], reason[33], *tname;
  olong ver, iRow, flagVer;

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("EditOTF", 1, 0, 0, NULL, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  myInput = EditOTFin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input FITS file name */
  for (i=0; i<128; i++) infile[i] = 0;
  ObitInfoListGet(myInput, "infile", &type, dim, infile, err);

  /* Flag table version */
  ObitInfoListGet(myInput, "FlagVer", &type, dim, &itemp, err);
  flagVer = itemp;

  /* Target name to flag "    " => all */
  for (i=0; i<20; i++) Target[i] = 0;
  ObitInfoListGet(myInput, "Target", &type, dim, Target, err);

  /* Time range to flag (days) */
  ObitInfoListGet(myInput, "TimeRange", &type, dim, darr, err);
  TimeRange[0] = darr[0];
  TimeRange[1] = darr[1];

  /* Range of Channels to flag */
  ObitInfoListGet(myInput, "ChanRange", &type, dim, iarr, err);
  bchan = iarr[0];
  echan = iarr[1];

  /* Polarization to flag "    "=> all */
  for (i=0; i<5; i++) stokes[i] = 0;
  ObitInfoListGet(myInput, "StokesFlag", &type, dim, stokes, err);

  /* Feed to flag, 0=> all */
  ObitInfoListGet(myInput, "FeedFlag", &type, dim, &itemp, err);
  FeedFlag = itemp;

   /* Reason */
  for (i=0; i<33; i++) reason[i] = 0;
  ObitInfoListGet(myInput, "Reason", &type, dim, reason, err);

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
  
  /* Look up Target number */
  /* Instantiate/Create output Flag Table */
  tname = g_strconcat ("OTFFlag table for: ",inData->name, NULL);
  ver = 1;
  TargetTable = newObitTableOTFTargetValue(tname, (Obit*)inData, &ver, 
				       OBIT_IO_ReadOnly, err);
  dim[0] = strlen(Target);
  dim[1] = 1;
  ObitTableOTFTargetLookup (TargetTable, dim, Target, iarr, err);
  TargID = iarr[0];  /* Target ID */
  g_free (tname);

  /* Instantiate/Create output Flag Table */
  tname = g_strconcat ("OTFFlag table for: ",inData->name, NULL);
  FlagTable = newObitTableOTFFlagValue(tname, (Obit*)inData, &flagVer, 
				       OBIT_IO_ReadWrite, err);
  g_free (tname);

  /* Open table */
  if ((ObitTableOTFFlagOpen (FlagTable, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output Flag table");
  }

  /* Create Table Row */
  FlagRow = newObitTableOTFFlagRow (FlagTable);
  
  /* Attach  row to output buffer */
  ObitTableOTFFlagSetRow (FlagTable, FlagRow, err);

  /* Fill in Flag row */
  FlagRow->TargetID = TargID;
  FlagRow->Feed     = FeedFlag;
  FlagRow->TimeRange[0] = TimeRange[0];
  FlagRow->TimeRange[1] = TimeRange[1];
  FlagRow->chans[0] = bchan;
  FlagRow->chans[1] = echan;
  FlagRow->pFlags[0] = 0;
  if (stokes[0]!=' ') FlagRow->pFlags[0] += 1;
  if (stokes[1]!=' ') FlagRow->pFlags[0] += 2;
  if (stokes[2]!=' ') FlagRow->pFlags[0] += 4;
  if (stokes[3]!=' ') FlagRow->pFlags[0] += 8;
  strncpy (FlagRow->reason, reason, 24);
  
  /* write row */
  iRow = FlagTable->myDesc->nrow+1;
  ObitTableOTFFlagWriteRow (FlagTable, iRow, FlagRow, err);
   
  /* Close Flag table */
  if ((ObitTableOTFFlagClose (FlagTable, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing Flag Table file");
  }
  
 /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
  
  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  /* cleanup */
  myInput = ObitInfoListUnref(myInput);  /* delete input list */
  inData = ObitUnref(inData);
  TargetTable = ObitUnref(TargetTable);
  FlagTable = ObitUnref(FlagTable);
  FlagRow   = ObitUnref(FlagRow);

 
  return ierr;
} /* end of main */

ObitInfoList* EditOTFin (int argc, char **argv, ObitErr *err)
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
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *input_file="EditOTF.in", *arg;
  gchar *strTemp;
  gboolean init=FALSE;
  odouble darray[2];
  olong    itemp, iarray[2];
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

    } else if (strcmp(arg, "-Target") == 0) { /* Target */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "Target", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-TimeRange") == 0) { /* Time Range */
      dim[0] = 2;
      darray[0] = strtod(argv[++ax], NULL);
      darray[1] = strtod(argv[++ax], NULL);
      ObitInfoListPut (list, "TimeRange", OBIT_double, dim, darray, err);


    } else if (strcmp(arg, "-ChanRange") == 0) { /* Channel Range */
      dim[0] = 2;
      iarray[0] = strtol(argv[++ax], NULL, 0);
      iarray[1] = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "ChanRange", OBIT_long, dim, iarray, err);

    } else if (strcmp(arg, "-StokesFlag") == 0) { /* Stokes Flag */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "StokesFlag", OBIT_string, dim, strTemp);

    } else if (strcmp(arg, "-FeedFlag") == 0) { /* Feed flag */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "FeedFlag", OBIT_long, dim, &itemp, err);

    } else if (strcmp(arg, "-Reason") == 0) {  /* Stokes Flag */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "Reason", OBIT_string, dim, strTemp);

    } else { /* unknown argument */
      Usage();
    }
  }
  
  /* Read defaults if no file specified */
  if (!init) ObitParserParse (input_file, list, err);

  return list;
} /* end EditOTFin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: EditOTF -input file [-args]\n");
    fprintf(stderr, "Edit an Obit/OTF data file\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def EditOTF.in\n");
    fprintf(stderr, "  -arg: any data selection item in input file \n");
    
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
  odouble darray[2];
  olong    itemp, iarray[2];
  ObitInfoList *out = newObitInfoList();

  /* add parser items */
  /* input FITS file name */
  strTemp = "OTFdata.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "infile", OBIT_string, dim, strTemp, err);

  /* Flag table version */
  dim[0] = 1;
  itemp = 1;
  ObitInfoListPut (out, "FlagVer", OBIT_long, dim, &itemp, err);

  /* Target to Flag */
  strTemp = "                ";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "Target", OBIT_string, dim, strTemp, err);

  /* Time Range */
  dim[0] = 2;
  darray[0] = 0.0;
  darray[1] = 1.0e20;
  ObitInfoListPut (out, "TimeRange", OBIT_double, dim, darray, err);

  /* Channel Range */
  dim[0] = 2;
  iarray[0] = 1;
  iarray[1] = 0;
  ObitInfoListPut (out, "ChanRange", OBIT_long, dim, iarray, err);

  /* Stokes Flag */
  strTemp = "    ";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "StokesFlag", OBIT_string, dim, strTemp, err);

  /* Feed flag */
  dim[0] = 1;
  itemp = 0;
  ObitInfoListPut (out, "FeedFlag", OBIT_long, dim, &itemp, err);

  /* Stokes Flag */
  strTemp = "No reason given     ";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "Reason", OBIT_string, dim, strTemp, err);

  return out;
} /* end defaultInputs */

