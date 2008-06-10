/* $Id: FitTip.c,v 1.1.1.1 2004/07/19 17:04:43 bcotton Exp $                            */
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
#include "ObitIOOTFFITS.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include <errno.h>

/* internal prototypes */
/* Get inputs */
ObitInfoList* FitTipin (int argc, char **argv, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Get Data */
void GetData (gchar *infile, olong maxdata, ofloat *el, 
	      ofloat *T1Cal, ofloat *T2Cal, ofloat *T1Off, ofloat *T2Off, 
	      olong *ndata, ObitErr *err);
/* Average cal signals */
void AvgCal (gint ndata, ofloat *el, 
	      ofloat *T1Cal, ofloat *T2Cal, ofloat *T1Off, ofloat *T2Off, 
	     ofloat *cal);
/* Fit tipping curves */
void FitData (gint ndata, ofloat *el, 
	      ofloat *T1Cal, ofloat *T2Cal, ofloat *T1Off, ofloat *T2Off, 
	      ofloat *Tatm, ofloat *Trx);

/* Program globals */
/* Cal equivalent in Jy */
odouble calValue;

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Fit Tipping scan(s)                                                  */
/*----------------------------------------------------------------------- */
{
  olong  i, ierr=0;
  ObitInfoList *myInput = NULL;
  ObitSystem *mySystem= NULL;
  ObitErr *err= NULL;
  gchar *FITSdir[] = {"FITSdata/"};
  gchar infile[128];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong ndata;
#define MAXDATA 20000  /* Maximum number of data points */
  ofloat el[MAXDATA], T1Cal[MAXDATA], T2Cal[MAXDATA], T1Off[MAXDATA], T2Off[MAXDATA];
  ofloat cal[2], Tatm[2], Trx[2];

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("FitTip", 1, 0, 0, NULL, 1, FITSdir, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Startup - parse command line */
  ierr = 0;
  myInput = FitTipin (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return ierr;

  /* Get inputs */
  /* input file name */
  for (i=0; i<128; i++) infile[i] = 0;
  ObitInfoListGet(myInput, "infile", &type, dim, infile, err);
  infile[dim[0]] = 0;  /* null terminate */

  /* Cal value in Jy */
  ObitInfoListGet(myInput, "calValue", &type, dim, &calValue, err);

  /* Get data from input file */
  GetData (infile, MAXDATA, el, T1Cal, T2Cal, T1Off, T2Off, &ndata, err);
  
  /* Average cals */
  AvgCal (ndata, el, T1Cal, T2Cal, T1Off, T2Off, cal);

  /* Fit tipping curves */
  FitData (ndata, el, T1Cal, T2Cal, T1Off, T2Off, Tatm, Trx);

  /* show any errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) return ierr;
   
   /* Shutdown Obit */
   mySystem = ObitSystemShutdown (mySystem);
   
   /* cleanup */
   myInput = ObitInfoListUnref(myInput);  /* delete input list */
   
   return ierr;
} /* end of main */

ObitInfoList* FitTipin (int argc, char **argv, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Parse control info from command line                                  */
/*   Input:                                                               */
/*      argc   Number of arguments from command line                      */
/*      argv   Array of strings from command line                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   return  parser list                                                  */
/*----------------------------------------------------------------------- */
{
  olong ax;
  gchar *input_file="FitTip.in", *arg;
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
} /* end FitTipin */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: FitTip -input file\n");
    fprintf(stderr, "Convert an external file format to Obit/OTF\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def FitTip.in\n");
    
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
  /* input file name */
  strTemp = "OTFdata.fits";
  dim[0] = strlen (strTemp);
  ObitInfoListPut (out, "infile", OBIT_string, dim, strTemp, err);

  /* Jy equivalent of cal */
  dim[0] = 1;
  dtemp = 1.0;
  ObitInfoListPut (out, "calValue", OBIT_double, dim, &dtemp, err);

  return out;
} /* end defaultInputs */

void GetData (gchar *infile, olong maxdata, ofloat *el, 
	      ofloat *T1Cal, ofloat *T2Cal, ofloat *T1Off, ofloat *T2Off, 
	      olong *ndata, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get Data  from text file                                              */
/*   Input:                                                               */
/*      infile   input file name                                          */
/*      maxdata  maximum number of data points                            */
/*      el       Array of elevations in rad                               */
/*      T1Cal    Array of T1 Cal on                                       */
/*      T2Cal    Array of T2 Cal on                                       */
/*      T1Off    Array of T1 Cal off                                      */
/*      T1Off    Array of T2 Cal off                                      */
/*   Output:                                                              */
/*      ndata    number of data points                                    */
/*      err      Obit return error stack                                  */
/*----------------------------------------------------------------------- */
{
  olong scan, sample, Xcaloff, Xcalon, Ycaloff, Ycalon;
  ofloat ra, dec, azx, elx, nscan;
  odouble mjd;
  gboolean done = FALSE;
  gchar *errMsg;
  FILE *myFile = NULL;
  gchar myLine[120];

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(infile!=NULL);
  g_assert(maxdata>10);
  g_assert(el!=NULL);
  g_assert(T1Cal!=NULL);
  g_assert(T2Cal!=NULL);
  g_assert(T1Off!=NULL);
  g_assert(T1Off!=NULL);

 /* Open input text file */
  myFile = fopen (infile, "rt");
  if (myFile==NULL) {
    Obit_log_error(err, OBIT_Error, "ERROR opening file %s", infile);
    errMsg = strerror(errno);
    Obit_log_error(err, OBIT_Error, "%s", errMsg);
    return;
  }

  /* Loop over file */
  *ndata = 0;
  while (!done) {
    /* read data */
    fgets(myLine, 120, myFile);
    if (ferror(myFile)) {
      Obit_log_error(err, OBIT_Error, "ERROR opening file %s", infile);
      errMsg = strerror(errno);
      Obit_log_error(err, OBIT_Error, "%s", errMsg);
      return;
    }
    done = feof(myFile); /* hit EOF? */
    if (done) break;
    
    /* Parse text line */
    nscan = sscanf (myLine,"%d %d %lf %f %f %f %f %d %d %d %d",
		    &scan, &sample, &mjd, &ra, &dec, &azx, &elx, 
		    &Xcaloff, &Xcalon, &Ycaloff, &Ycalon);
    if (nscan<11) {
      Obit_log_error(err, OBIT_Error, "ERROR parsing line %s", myLine);
      errMsg = strerror(errno);
      Obit_log_error(err, OBIT_Error, "%s", errMsg);
      return;
    }

    /* Save data in arrays */
    if (*ndata<maxdata) {
      el[*ndata] = elx;
      T1Off[*ndata] = Xcaloff;
      T1Cal[*ndata] = Xcalon;
      T2Off[*ndata] = Ycaloff;
      T2Cal[*ndata] = Ycalon;
      (*ndata)++;
    } else { /* blew core */
      Obit_log_error(err, OBIT_Error, "Too many data points");
      return;
    }


  } /* end loop over file */

   /* Close input file */
   fclose(myFile);
  
} /* end GetData */

void AvgCal (gint ndata, ofloat *el, 
	      ofloat *T1Cal, ofloat *T2Cal, ofloat *T1Off, ofloat *T2Off, 
	      ofloat *cal)
/*----------------------------------------------------------------------- */
/*  Get Data  from text file                                              */
/*   Input:                                                               */
/*      ndata    number of data points                                    */
/*      el       Array of elevations in rad                               */
/*      T1Cal    Array of T1 Cal on                                       */
/*      T2Cal    Array of T2 Cal on                                       */
/*      T1Off    Array of T1 Cal off                                      */
/*      T1Off    Array of T2 Cal off                                      */
/*   Output:                                                              */
/*      cal      Average cal value for T1 and T2                          */
/*----------------------------------------------------------------------- */
{
  ofloat sum[2]={0.0,0.0}, fblank=ObitMagicF();
  olong i, count[2]={0,0};

  for (i=0; i<ndata; i++) {
    if ((T1Off[i]!=fblank) && (T1Cal[i]!=fblank)) {
      sum[0] += T1Cal[i] - T1Off[i];
      count[0]++;
    }
    if ((T2Off[i]!=fblank) && (T2Cal[i]!=fblank)) {
      sum[1] += T2Cal[i] - T2Off[i];
      count[1]++;
    }
  } /* end loop over data */
  
  /* average */
  if (count[0]>0) cal[0] = sum[0]/count[0];
  else cal[0] = fblank;
  if (count[1]>0) cal[1] = sum[1]/count[1];
  else cal[1] = fblank;

  /* debug */
  fprintf (stdout, "Cal averages = %f %f\n", cal[0], cal[1]);
} /* end AvgCal */

void FitData (gint ndata, ofloat *el, 
	      ofloat *T1Cal, ofloat *T2Cal, ofloat *T1Off, ofloat *T2Off, 
	      ofloat *Tatm, ofloat *Trx)
/*----------------------------------------------------------------------- */
/*  Fit each T1 and T2 data in terms of rx and atm components             */
/*  T = Trx + Tatm * airmass, where airmass = 1.0 / cos (za)              */
/*   Input:                                                               */
/*      ndata    number of data points                                    */
/*      el       Array of elevations in rad                               */
/*      T1Cal    Array of T1 Cal on                                       */
/*      T2Cal    Array of T2 Cal on                                       */
/*      T1Off    Array of T1 Cal off                                      */
/*      T1Off    Array of T2 Cal off                                      */
/*   Output:                                                              */
/*     Tatm      The Atmospheric contribution per airmass, T1, T2         */
/*     Trx       The receiver contribution per airmass, T1, T2            */
/*------------------------------------------------------------------------ */
{
  odouble x, y, Sxx, Sxy;
  odouble sumxT1=0.0, sumyT1=0.0, sumxxT1=0.0, sumyyT1=0.0, sumxyT1=0.0;
  odouble sumxT2=0.0, sumyT2=0.0, sumxxT2=0.0, sumyyT2=0.0, sumxyT2=0.0;
  ofloat fblank=ObitMagicF();
  olong i, countT1=0, countT2=0;

  /* Do sums */
  for (i=0; i<ndata; i++) {
    x = 1.0 / cos (1.5708-el[i]);
    if (T1Off[i]!=fblank) {
      y = T1Off[i];
      countT1++;
      sumxT1 += x;
      sumyT1 += y;
      sumxxT1 += x*x;
      sumyyT1 += y*y;
      sumxyT1 += x*y;
    }

    if (T2Off[i]!=fblank) {
      y = T2Off[i];
      countT2++;
      sumxT2 += x;
      sumyT2 += y;
      sumxxT2 += x*x;
      sumyyT2 += y*y;
      sumxyT2 += x*y;
    }
  } /* end loop over data */
  
  /* Determine coefficients */
  if (countT1>0) {
    Sxx = sumxxT1 - sumxT1*sumxT1/countT1;
    Sxy = sumxyT1 - sumxT1*sumyT1/countT1;
    Tatm[0] = Sxy/Sxx;
    Trx[0]  = (sumyT1-Tatm[0]*sumxT1)/countT1;
  } else {
    Tatm[0] = fblank;
    Trx[0]  = fblank;
  }

  if (countT2>0) {
    Sxx = sumxxT2 - sumxT2*sumxT2/countT2;
    Sxy = sumxyT2 - sumxT2*sumyT2/countT2;
    Tatm[1] = Sxy/Sxx;
    Trx[1]  = (sumyT2-Tatm[0]*sumxT2)/countT2;
  } else {
    Tatm[1] = fblank;
    Trx[1]  = fblank;
  }


  /* debug */
  fprintf (stdout, "TRx  = %f %f\n", Trx[0],  Trx[1]);
  fprintf (stdout, "Tatm = %f %f\n", Tatm[0], Tatm[1]);
  
} /* end FitData */
