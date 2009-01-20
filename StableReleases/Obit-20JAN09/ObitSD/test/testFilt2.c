/* $Id$                            */
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

#include "ObitTimeFilter.h"
#include "ObitPlot.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitOTF.h"

#define MAXSAMPLE 5000   /* Maximum number of samples in scan */

/* internal prototypes */
/* Get inputs */
ObitInfoList* testFiltin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);

static void ReadScan (ObitOTF *inOTF, olong detect, olong scan, olong *Count,
		       ofloat Time[MAXSAMPLE], ofloat RA[MAXSAMPLE], 
		       ofloat Dec[MAXSAMPLE], ofloat *data, 
		      ObitErr *err);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Testbed program                                         */
/*----------------------------------------------------------------------- */
{
  oint i, nTime, ndetect, ierr = 0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitErr *err= NULL;
  ObitOTF *inOTF=NULL;
  ofloat Time[MAXSAMPLE], RA[MAXSAMPLE], Dec[MAXSAMPLE], data[2*MAXSAMPLE];
  olong Count;
  ofloat parms[1], fblank = ObitMagicF();
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar strTemp[121], *FITSdir[] = {"PythonData/"};
  gchar *inFile = "GCXbandCalSubOTF.fits";
  olong scan, nrec, detect, disk;
  ObitTimeFilter *filter = NULL;
  ObitPlot *plot = NULL;

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("testFilt", 1, 0, 0, NULL, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  myInput = testFiltin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

 /* error check */
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

   /* Set data */
  disk = 1;
  nrec = 1000;
  inOTF  = newObitOTF("Input data");
  ObitOTFSetFITS(inOTF, nrec, disk, inFile,err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Open/close input OTF to fully instantiate */
  if ((ObitOTFOpen (inOTF, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input FITS file %s", inFile);
  
  /* Close */
  if ((ObitOTFClose (inOTF, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing input file");
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
  
 /* read data */
  scan = 63;
  detect = -1;
  ReadScan (inOTF, detect, scan, &Count, Time, RA, Dec, data, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Clip data */
  for (i=0; i<Count; i++)  if (data[i]>0.1)  data[i]=0.1;
  for (i=0; i<Count; i++)  if (data[i]<-0.1) data[i]=-0.1;

 /* Test filter */
  nTime = Count;
  ndetect = 1;
  filter = newObitTimeFilter("testFilter", nTime, ndetect);

  /* copy to filter */
  for (i=0; i<nTime; i++)  filter->timeData[0][i] = data[i];

 /* test plot */
  plot = newObitPlot ("myPlot");

  /* Transform to frequency */
  ObitTimeFilter2Freq (filter);

  /* Plot labels */
  strcpy (strTemp, "Time Series Before");
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "TITLE", OBIT_string, dim, strTemp);
  /* plot time sequence */
  ObitPlotInitPlot (plot, NULL, err);
  ObitPlotXYPlot (plot, 2, nTime, NULL, filter->timeData[0], err);


  for (i=0; i<10; i++) filter->freqData[0][i] = 0.0;
  ObitPlotInitPlot (plot, NULL, err);
  ObitPlotXYPlot (plot, 2, nTime, NULL, filter->freqData[0], err);
  /* Apply filter */
  parms[0] = 0.25;
  ObitTimeFilterFilter (filter, -1, OBIT_TimeFilter_LowPass, parms, err);
  
  ObitPlotInitPlot (plot, NULL, err);
  ObitPlotXYPlot (plot, 2, nTime, NULL, filter->freqData[0], err);

  /* FT back to time */
  ObitTimeFilter2Time (filter);

  /* Plot labels */
  strcpy (strTemp, "Time Series After");
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "TITLE", OBIT_string, dim, strTemp);

  ObitPlotInitPlot (plot, NULL, err);
  /*  ObitPlotXYPlot (plot, 2, nTime, NULL, filter->freqData[0], err);*/
  ObitPlotXYPlot (plot, 2, nTime, NULL, filter->timeData[0], err);

 /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
  
  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  /* cleanup */
  myInput = ObitInfoListUnref(myInput);  /* delete input list */
  filter = ObitTimeFilterUnref(filter);
  plot = ObitPlotUnref(plot);
 
  return ierr;
} /* end of main */

ObitInfoList* testFiltin (int argc, char **argv, ObitErr *err)
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
  gchar *input_file="testFilt.in", *arg;
  gchar *strTemp;
  gboolean init=FALSE;
  odouble darray[2];
  olong    itemp, iarray[2];
  ObitInfoList* list;

  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* debug - no command line input yet */
  if (list) return list;

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
} /* end testFiltin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: testFilt -input file [-args]\n");
    fprintf(stderr, "Edit an Obit/OTF data file\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def testFilt.in\n");
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

/**
 * Reads single scan.
 * Gets calibrator information from the Target table (position, flux density)
 * \param inOTF    Input OTF data. 
 * \param detect   detector number (0-rel), -1 => all
 * \param scan     which scan to process
 * \param Count    Number of samples read
 * \param Time     Time tags of data
 * \param RA       RAs of data
 * \param Dec      Decs of data
 * \param data     Data array MAXSAMPLE*ndetect in size
 * \param err      Error stack, returns if not empty.
 */
static void ReadScan (ObitOTF *inOTF, olong detect, olong scan, olong *Count,
		       ofloat Time[MAXSAMPLE], ofloat RA[MAXSAMPLE], 
		       ofloat Dec[MAXSAMPLE], ofloat *data, 
		       ObitErr *err)
  {
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOCode retCode;
  gboolean doCalSelect;
  olong i, doCal, scans[2], ndetect,  mdetect, loDet, hiDet, iDet;
  ofloat *rec;
  gchar *routine = "ReadScan";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));

  /* How many detectors? */
  ndetect = inOTF->geom->numberDetect;
  /* How many to use */
  if (detect<0) {
    mdetect = ndetect;
    loDet = 0;
    hiDet = ndetect-1;
  } else {
    mdetect = 1;
    loDet = detect;
    hiDet = detect;
  }

  /* Select scan on input */
  doCalSelect = TRUE;
  dim[0] = 1;
  ObitInfoListPut(inOTF->info, "doCalSelect", OBIT_bool, dim, &doCalSelect, err);
  doCal = 1;
  ObitInfoListPut(inOTF->info, "DOCAL", OBIT_bool, dim, &doCal, err);
  scans[0] = scan; scans[1] = scan;
  dim[0] = 2;
  ObitInfoListPut(inOTF->info, "SCANS", OBIT_long, dim, scans, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

   /* open OTF data to fully instantiate  */
  retCode = ObitOTFOpen (inOTF, OBIT_IO_ReadCal, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Read data into storage */
  (*Count) = 0;      /* no. values in arrays */

  /* loop reading data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
    if (retCode==OBIT_IO_EOF) break; /* done? */

    /* Record pointer */
    rec = inOTF->buffer;
  
    /* Loop over buffer */
       for (i=0; i<inOTF->myDesc->numRecBuff; i++) {

	 /* Check -debug
	 if (scan!=((olong)(rec[inOTF->myDesc->ilocscan]+0.5))) {
	   fprintf (stderr,"Holy shit Batman this %d is not scan %d\n",
		     ((olong)(rec[inOTF->myDesc->ilocscan]+0.5)), scan);
	   exit(9);
	 } */
	   
	 /* Accumulate values  */
	 if (*Count<MAXSAMPLE) {	   
	   Time[*Count]  = rec[inOTF->myDesc->iloct];
	   RA[*Count]    = rec[inOTF->myDesc->ilocra];
	   Dec[*Count]   = rec[inOTF->myDesc->ilocdec];
	   for (iDet=loDet; iDet<=hiDet; iDet++) 
	     data[(iDet-loDet)*MAXSAMPLE+(*Count)] = rec[inOTF->myDesc->ilocdata+iDet];
	   (*Count)++;
	 } /* end  accumulate */
	 rec += inOTF->myDesc->lrec; /* Data record pointer */	 
       } /* end loop over buffer load */
  } /* end loop reading data */ 
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

}  /* end ReadScan */

