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

#include "ObitTimeFilter.h"
#include "ObitPlot.h"
#include "ObitSystem.h"
#include "ObitParser.h"
/* internal prototypes */
/* Get inputs */
ObitInfoList* testFiltin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Make entries in OTF Flag table                                       */
/*----------------------------------------------------------------------- */
{
#define NDATA 64
  oint i, nTime, ndetect, ierr = 0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitErr *err= NULL;
  ofloat x[NDATA], y[NDATA], parms[1];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar strTemp[121], *FITSdir[] = {"FITSdata/"};
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

  /* Test filter */
  nTime = NDATA;
  ndetect = 1;
  filter = newObitTimeFilter("testFilter", nTime, ndetect);

  /* Set data */
  for (i=0; i<nTime; i++) {
    x[i] = (ofloat)i;
    y[i] = sin (2*G_PI*i/nTime) + cos (20*G_PI*i/nTime);
  }

  /* copy to filter */
  for (i=0; i<nTime; i++)  filter->timeData[0][i] = y[i];

 /* test plot */
  plot = newObitPlot ("myPlot");

  /* Transform to frequency */
  ObitTimeFilter2Freq (filter);

  /* Plot labels */
  strcpy (strTemp, "Time Series Before");
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "TITLE", OBIT_string, dim, strTemp);

  ObitPlotInitPlot (plot, NULL, err);
  /*  ObitPlotXYPlot (plot, 2, nTime, NULL, filter->freqData[0], err);*/
  ObitPlotXYPlot (plot, 2, nTime, NULL, filter->timeData[0], err);
  /* Apply filter */
  parms[0] = 0.25;
  ObitTimeFilterFilter (filter, -1, OBIT_TimeFilter_LowPass, parms);
  
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

