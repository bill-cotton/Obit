/* $Id: ObitTimeFilter.c,v 1.9 2008/02/20 15:12:03 bcotton Exp $  */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include <math.h>
#include "ObitTimeFilter.h"
#include "ObitPlot.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTimeFilter.c
 * ObitTimeFilter class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitTimeFilter";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitTimeFilterClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitTimeFilterClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitTimeFilterInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitTimeFilterClear (gpointer in);

/** Private: interpolate to replace blanks in a 1-D ObitFarray. */
static void  ObitTimeFilterDeblank (ObitTimeFilter *in, ObitFArray *array);

/** Private: Set Class function pointers. */
static void ObitTimeFilterClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name    An optional name for the object.
 * \param nTime   Number of times in arrays to be filtered
 *                It is best to add some extra padding (10%) to allow a smooth
 *                transition from the end of the sequence back to the beginning.
 *                Remember the FFT algorithm assumes the function is periodic.
 * \param nSeries Number of time sequences to be filtered
 * \return the new object.
 */
ObitTimeFilter* newObitTimeFilter (gchar* name, olong nTime, olong nSeries)
{
  ObitTimeFilter* out;
  olong i, rank, dim[1];
  olong ndim, naxisr[1], naxisc[1], pos[1] = {0};

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTimeFilterClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitTimeFilter));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitTimeFilterInit((gpointer)out);

  /* save inputs */
  out->nTime   = nTime;
  out->nFreq   = 1 + nTime / 2;
  out->nSeries = nSeries;

  /* Create arrays */
  out->timeData   = g_malloc0 (nSeries*sizeof(ofloat*));
  out->freqData   = g_malloc0 (nSeries*sizeof(ofloat*));
  out->times      = g_malloc0 (out->nTime*sizeof(ofloat));
  out->freqs      = g_malloc0 (out->nFreq*sizeof(ofloat));
  out->timeSeries = g_malloc0 (nSeries*sizeof(ObitFArray*));
  out->freqSeries = g_malloc0 (nSeries*sizeof(ObitCArray*));
  ndim = 1;
  naxisr[0] = out->nTime;
  naxisc[0] = out->nFreq;
  for (i=0; i<nSeries; i++) {
    out->timeSeries[i] = ObitFArrayCreate("Time", ndim, naxisr);
    out->timeData[i]   = ObitFArrayIndex (out->timeSeries[i], pos); /* save pointer to data */
    out->freqSeries[i] = ObitCArrayCreate("Freq", ndim, naxisc);
    out->freqData[i]   = ObitCArrayIndex (out->freqSeries[i], pos); /* save pointer to data */
  }

  /* Create FFT objects */
  rank = 1;
  dim[0] = out->nTime;
  out->FFTFor = newObitFFT ("Time2FreqFFT", OBIT_FFT_Forward, 
			   OBIT_FFT_HalfComplex, rank, dim);
  out->FFTRev = newObitFFT ("Freq2TimeFFT", OBIT_FFT_Reverse, 
			   OBIT_FFT_HalfComplex, rank, dim);

  /* Interpolator object */
  out->interp = newObitFInterpolateCreate ("TimeSeriesInterpolator", 
					  out->timeSeries[0], NULL, 3);

  return out;
} /* end newObitTimeFilter */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitTimeFilterGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTimeFilterClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitTimeFilterGetClass */

/**
 * Change the sizes of the arrays on the input object to nTime.
 * \param in    ObitTimeFilter to resize.
 * \param nTime Number of times in arrays to be filtered
 *              It is best to add some extra padding (10%) to allow a smooth
 *              transition from the end of the sequence back to the beginning.
 *              Remember the FFT algorithm assumes the function is periodic.
 */
void ObitTimeFilterResize (ObitTimeFilter *in, olong nTime)
{
  olong i, rank, dim[1];
  olong ndim, naxisr[1], naxisc[1], pos[1] = {0};

  /* error checks */
  g_assert (ObitTimeFilterIsA(in));

  /* Don't bother if it's the same size */
  if (nTime==in->nTime) return;

  /* Lock ObitObjects against other threads */
  ObitThreadLock(in->thread);

  /* save inputs */
  in->nTime   = nTime;
  in->nFreq   = 1 + nTime / 2;

  /* Realloc float arrays */
  in->times  = g_realloc (in->times, in->nTime*sizeof(ofloat));
  in->freqs  = g_realloc (in->freqs, in->nFreq*sizeof(ofloat));
  
  /* Loop over series reallocating */
  ndim = 1;
  naxisr[0] = in->nTime;
  naxisc[0] = in->nFreq;
  for (i=0; i<in->nSeries; i++) {
    in->timeSeries[i] = ObitFArrayRealloc(in->timeSeries[i], ndim, naxisr);
    in->timeData[i]   = ObitFArrayIndex (in->timeSeries[i], pos); /* save pointer to data */
    in->freqSeries[i] = ObitCArrayRealloc(in->freqSeries[i], ndim, naxisc);
    in->freqData[i]   = ObitCArrayIndex (in->freqSeries[i], pos); /* save pointer to data */
  } /* end loop over series */

  /* Recreate FFT objects */
  in->FFTFor = ObitFFTUnref(in->FFTFor);  /* out with the old */
  rank = 1;
  dim[0] = in->nTime;
  in->FFTFor = newObitFFT ("Time2FreqFFT", OBIT_FFT_Forward, 
			   OBIT_FFT_HalfComplex, rank, dim);
  in->FFTRev = ObitFFTUnref(in->FFTRev);  /* out with the old */
  in->FFTRev = newObitFFT ("Freq2TimeFFT", OBIT_FFT_Reverse, 
			   OBIT_FFT_HalfComplex, rank, dim);

  /* Replace array on interpolator member */
  ObitFInterpolateReplace (in->interp, in->timeSeries[0]);

  /* Unlock ObitObjects */
  ObitThreadUnlock(in->thread);
 
} /* end ObitTimeFilterResize */

/**
 * Convert time-ordered data stream into regular time grid.
 * Will resize in if needed.
 * Data will be averaging into time bins.
 * \param in       Object with TimeFilter structures.
 * \param seriesNo Which time/frequency series to apply to (0-rel)
 * \param dTime    Increment of desired time grid (days)
 * \param nTime    Number of times in times, data
 * \param times    Array of times (days)
 * \param data     Array of data elements corresponding to times.
 */
void ObitTimeFilterGridTime (ObitTimeFilter *in, olong seriesNo,
			     ofloat dTime, olong nTime, ofloat *times, ofloat *data)
{
  olong i, j, nt, ntt, icell, delta, half, *count=NULL;
  ofloat idt, fblank = ObitMagicF();

  /* How many times? */
  nt = (times[nTime-1] - times[0]) / MAX( 1.0e-20, dTime);
  ntt = ObitFFTSuggestSize (nt);

  /* Resize if needed */
  ObitTimeFilterResize (in, ntt);

  /* Allocate work array */
  count = g_malloc0(10+ntt*sizeof(olong));

  /* Lock ObitObjects against other threads */
  ObitThreadLock(in->thread);

  /* Time and frequency increments */
  in->dTime = dTime;
  idt = 1.0 / dTime;
  in->dFreq = 1.0 / (in->nTime*dTime*86400.0); 

  /* Loop initializing arrays */
  half = in->nTime/2;  /* in->nTime should always be even */
  for (i=0; i<in->nFreq; i++) in->freqs[i] = i*in->dFreq;
  for (i=0; i<in->nTime; i++) {
    in->timeData[seriesNo][i] = 0.0;
    count[i] = 0;
    /* Time - center at the edges */
    if (i<half) j = half+i;
    else j = i-half;
    in->times[j] = i*dTime; 
  }

  /* Loop gridding - center at edges */
  delta = in->nTime/2 + (in->nTime-nTime)/2;
  for (i=0; i<nTime/2; i++) {
    icell = (olong)((times[i]-times[0])*idt + 0.5);
    icell = MAX (0, MIN (half-1, icell));
    if (data[i]!=fblank) {
      in->timeData[seriesNo][icell+delta] += data[i];
      count[icell+delta]++;
    }
  }
  half = nTime/2;
  for (i=nTime/2; i<nTime; i++) {
    icell = (olong)((times[i]-times[0])*idt + 0.5);
    icell = MAX (half, MIN (in->nTime-1, icell));
    if (data[i]!=fblank) {
      in->timeData[seriesNo][icell-half] += data[i];
	count[icell-half]++;
    }
  }
  /* Loop normalizing */
  for (i=0; i<in->nTime; i++) {
    if (count[i]>0) {
      in->timeData[seriesNo][i] /= count[i];
    } else {  /* No data this time */
      in->timeData[seriesNo][i] = fblank;
    }
  }

  if (count) g_free(count); /* Cleanup */

  /* Unlock ObitObjects */
  ObitThreadUnlock(in->thread);
} /* end  ObitTimeFilterGridTime*/


/**
 * Copy time series to external form.
 * Time series data are copied by nearest time stamp
 * \param in       Object with TimeFilter structures.
 * \param seriesNo Which time/frequency series to apply to (0-rel)
 * \param nTime    Number of times in times, data
 * \param times    [in] Array of times (days)
 * \param data     [out] Array of date elements corresponding to times.
 */
void ObitTimeFilterUngridTime (ObitTimeFilter *in, olong seriesNo,
			       olong nTime, ofloat *times, ofloat *data)
{
  olong i, half, delta, icell;
  ofloat idt, fblank = ObitMagicF();

  /* Lock ObitObjects against other threads */
  ObitThreadLock(in->thread);

  /* loop over output data - center at the edges */
  idt = 1.0 / in->dTime;
  half = nTime/2;
  delta = in->nTime/2 + (in->nTime-nTime)/2;
  for (i=0; i<nTime/2; i++) {
    icell = (olong)((times[i]-times[0])*idt + 0.5);
    icell = MAX (0, MIN (half-1, icell));
    if (data[i]!=fblank) data[i] = in->timeData[seriesNo][icell+delta];
  }

  for (i=nTime/2; i<nTime; i++) {
    icell = (olong)((times[i]-times[0])*idt + 0.5);
    icell = MAX (half, MIN (in->nTime-1, icell));
    if (data[i]!=fblank) data[i] = in->timeData[seriesNo][icell-half];
  }

  /* Unlock ObitObjects */
  ObitThreadUnlock(in->thread);
} /* end  ObitTimeFilterUngridTime*/


/**
 * Fourier transform all time series to frequency series
 * Any blanked values in the time series are interpolated and then the time series 
 * is FFTed to the frequency domain.  A linear interpolation between the 
 * last valid point and the first valid point is made to reduce the wraparound edge
 * effects.
 * \param in  Object with TimeFilter structures.
 */
void ObitTimeFilter2Freq (ObitTimeFilter *in)
{
  olong i;

  /* error checks */
  g_assert (ObitTimeFilterIsA(in));

  /* Lock ObitObjects against other threads */
  ObitThreadLock(in->thread);

  /* Loop over series */
  for (i=0; i<in->nSeries; i++) {

    /* Interpolated any blanked values, ensure smooth wrap around */
    ObitTimeFilterDeblank(in, in->timeSeries[i]);
    
    /* Transform */
    ObitFFTR2C (in->FFTFor, in->timeSeries[i], in->freqSeries[i]);

  } /* end loop over series */

  /* Unlock ObitObjects */
  ObitThreadUnlock(in->thread);
 
} /* end ObitTimeFilter2Freq */

/**
 * Fourier transform filtered frequency series to time series
 * \param in  Object with TimeFilter structures.
 */
void ObitTimeFilter2Time (ObitTimeFilter *in)
{
  olong i, j;
  float norm;

  /* error checks */
  g_assert (ObitTimeFilterIsA(in));

  /* Lock ObitObjects against other threads */
  ObitThreadLock(in->thread);

  /* Loop over series */
  for (i=0; i<in->nSeries; i++) {

    /* Transform */
    ObitFFTC2R (in->FFTRev, in->freqSeries[i], in->timeSeries[i]);

    /* Normalize by n */
    norm = 1.0 / in->nTime;
    if (in->nTime>1) for (j=0; j<in->nTime; j++) in->timeData[i][j] *= norm;

  } /* end loop over series */

  /* Unlock ObitObjects */
  ObitThreadUnlock(in->thread);
 
} /* end ObitTimeFilter2Time */

/**
 * Apply specified filter to specified time series.
 *
 * Following Filters are supported:
 * \li OBIT_TimeFilter_LowPass - 
 *     Zeroes frequencies above a fraction, parm[0], of the highest.
 * \li OBIT_TimeFilter_HighPass - 
 *     Zeroes frequencies below a fraction, parm[0], of the highest.
 * \li OBIT_TimeFilter_NotchPass - 
 *     Zeroes frequencies not in frequency range parm[0]->parm[1]
 * \li OBIT_TimeFilter_NotchBlock - 
 *     Zeroes frequencies in frequency range parm[0]->parm[1]
 *
 * \param in       Object with TimeFilter structures.
 * \param seriesNo Which time/frequency series to apply to (0-rel), <0 => all
 * \param type     Filter type to apply
 * \param *parm    Parameters for filter, meaning depends on type.
 * \param err      Error stack
 */
void ObitTimeFilterFilter (ObitTimeFilter *in, olong seriesNo,
			   ObitTimeFilterType type, ofloat *parms, ObitErr *err)
{
  olong hi, lo, fLo, fHi, fCen, iSeries, iFreq;
  gchar *routine = "ObitTimeFilterFilter";

  /* error checks */
  g_assert (ObitTimeFilterIsA(in));
 
  /* Lock ObitObjects against other threads */
  ObitThreadLock(in->thread);

  /* Range of series to use */
  if (seriesNo>0) {
    lo = seriesNo;
    hi = seriesNo;
  } else {
    lo = 0;
    hi = in->nSeries-1;
  }
   
  /* Apply filter by type */
  switch (type) {
  case OBIT_TimeFilter_LowPass:   /* Low pass filter */
    /* Check parameters */
    if (!((parms[0]>=0.0) && (parms[0]<=1.0))) { /* Check range */
      Obit_log_error(err, OBIT_Error, "%s: parm[0] out of range [0,1] %g", 
		     routine, parms[0]);
      return;
    }

    /* Set range of frequency values to zero */
    fLo = (olong)((parms[0] * in->nFreq) + 0.99999);
    fLo *= 2; /* As complex */
    fHi = 2 * in->nFreq;

    /* Loop over series */
    for (iSeries=lo; iSeries<=hi; iSeries++) {
      for (iFreq = fLo; iFreq<fHi; iFreq++) in->freqData[iSeries][iFreq] = 0.0;
    } /* end loop over series */
    break;

  case OBIT_TimeFilter_HighPass:  /* High pass filter */
    /* Check parameters */
    if (!((parms[0]>=0.0) && (parms[0]<=1.0))) { /* Check range */
      Obit_log_error(err, OBIT_Error, "%s: parm[0] out of range [0,1] %g", 
		     routine, parms[0]);
      return;
    }

    /* Set range of frequency values to zero */
    fLo = 0;
    fHi = (olong)((parms[0] * in->nFreq) + 0.5) - 1;
    fHi *= 2; /* As complex */

    /* Loop over series */
    for (iSeries=lo; iSeries<=hi; iSeries++) {
      for (iFreq = fLo; iFreq<fHi; iFreq++) in->freqData[iSeries][iFreq] = 0.0;
    } /* end loop over series */
    break;

  case OBIT_TimeFilter_NotchPass:    /* Notch pass filter */
    /* Check parameters */
    if (!((parms[0]>=0.0) && (parms[0]<=1.0))) { /* Check range */
      Obit_log_error(err, OBIT_Error, "%s: parm[0] out of range [0,1] %g", 
		     routine, parms[0]);
      return;
    }
    if (!((parms[1]>=0.0) && (parms[1]<=1.0))) { /* Check range */
      Obit_log_error(err, OBIT_Error, "%s: parm[1] out of range [0,1] %g", 
		     routine, parms[1]);
      return;
    }

    /* Set range of frequency values to zero */
    fLo = (olong)((parms[0] * in->nFreq) + 0.5) - 1;
    fHi = (olong)((parms[1] * in->nFreq) + 0.99999);
    /* Make sure somthing in the middle */
    fCen = (olong)(0.5*(fLo+fHi+1) + 0.5) - 1;
    fLo = MAX(0, MIN (fLo, fCen));
    fHi = MIN(in->nFreq-1, MAX (fHi, fCen+1));
    fLo *= 2; /* As complex */
    fHi *= 2;

    /* Loop over low frequency part of series */
    for (iSeries=lo; iSeries<=hi; iSeries++) {
      for (iFreq = 0; iFreq<fLo; iFreq++) in->freqData[iSeries][iFreq] = 0.0;
    } /* end loop over lo series */

    /* Loop over high frequency part of series */
    for (iSeries=lo; iSeries<=hi; iSeries++) {
      for (iFreq = fHi; iFreq<2*in->nFreq; iFreq++) in->freqData[iSeries][iFreq] = 0.0;
    } /* end loop over hi series */
    break;

  case OBIT_TimeFilter_NotchBlock:   /* Notch block filter */
    /* Check parameters */
    if (!((parms[0]>=0.0) && (parms[0]<=1.0))) { /* Check range */
      Obit_log_error(err, OBIT_Error, "%s: parm[0] out of range [0,1] %g", 
		     routine, parms[0]);
      return;
    }
    if (!((parms[1]>=0.0) && (parms[1]<=1.0))) { /* Check range */
      Obit_log_error(err, OBIT_Error, "%s: parm[0] out of range [0,1] %g", 
		     routine, parms[1]);
      return;
    }

    /* Set range of frequency values to zero */
    fLo = (olong)((parms[0] * in->nFreq) + 0.99999);
    fHi = (olong)((parms[1] * in->nFreq) + 0.5);
    /* Make sure somthing in the middle */
    fCen = (olong)(0.5*(fLo+fHi+1) + 0.5) - 1;
    fLo = MAX(0, MIN (fLo, fCen));
    fHi = MIN(in->nFreq-1, MAX (fHi, fCen+1));
    fLo *= 2;  /* As complex */
    fHi *= 2;

    /* DEBUG print power spectrum
    iSeries = 0;
    for (iFreq=0; iFreq<in->nFreq/2; iFreq++) {
      fprintf (stderr,"%d %f %f\n",iFreq, iFreq*in->dFreq,
	       sqrt(in->freqData[iSeries][iFreq*2]*in->freqData[iSeries][iFreq*2] +
		    in->freqData[iSeries][iFreq*2+1]*in->freqData[iSeries][iFreq*2+1]));
    } */

    /* Loop over series */
    for (iSeries=lo; iSeries<=hi; iSeries++) {
      for (iFreq=fLo; iFreq<fHi; iFreq++) in->freqData[iSeries][iFreq] = 0.0;
    } /* end loop over series */

    /* DEBUG print power spectrum after
    iSeries = 0;
    fprintf (stderr,"zeroes %d - %d \n",fLo/2, fHi/2);
    for (iFreq=0; iFreq<in->nFreq/2; iFreq++) {
      fprintf (stderr,"%d %f %f\n",iFreq, iFreq*in->dFreq,
	       sqrt(in->freqData[iSeries][iFreq*2]*in->freqData[iSeries][iFreq*2] +
		    in->freqData[iSeries][iFreq*2+1]*in->freqData[iSeries][iFreq*2+1]));
    }  */
    break;

  default:                        /* Unknown - barf and die */
    g_assert_not_reached(); 
  }; /* end switch on filter type */
  
  
  /* Unlock ObitObjects */
  ObitThreadUnlock(in->thread);
} /* end ObitTimeFilterFilter */

/**
 * Apply specified filter to specified time series.
 *
 * Following Filters are supported:
 * \li OBIT_TimeFilter_LowPass - 
 *     Zeroes frequencies above freq[0] (Hz)
 * \li OBIT_TimeFilter_HighPass - 
 *     Zeroes frequencies below freq[0] (Hz)
 * \li OBIT_TimeFilter_NotchPass - 
 *     Zeroes frequencies not in frequency range freq[0]->freq[1]
 * \li OBIT_TimeFilter_NotchBlock - 
 *     Zeroes frequencies in frequency range freq[0]->freq[1]
 *
 * \param in       Object with TimeFilter structures.
 * \param seriesNo Which time/frequency series to apply to (0-rel), <0 => all
 * \param type     Filter type to apply
 * \param freq     Frequencies (Hz) for filter, meaning depends on type.
 * \param err      Error stack
 */
void ObitTimeFilterDoFilter (ObitTimeFilter *in, olong seriesNo,
			     ObitTimeFilterType type, ofloat *freq, ObitErr *err)
{
  ofloat parm[5];
  gchar *routine = "ObitTimeFilterDoFilter";

  /* error checks */
  g_assert (ObitTimeFilterIsA(in));

  /* Convert physical parameters to parm */
  parm[0] = freq[0] / (in->nFreq*in->dFreq);
  parm[1] = freq[1] / (in->nFreq*in->dFreq);

  /* Filter */
  ObitTimeFilterFilter (in, seriesNo, type, parm, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitTimeFilterDoFilter */

/**
 * Plot the power spectrum of a selected frequency series
 * 
 * \param in       ObitTimeFilter to plot
 * \param seriesNo Which frequency series to plot (0-rel)
 * \param label    If nonNULL, a label for the plot
 * \param err      Error stack, returns if not empty.
 */
void ObitTimeFilterPlotPower (ObitTimeFilter *in, olong seriesNo,
			      gchar *label, ObitErr *err)
{
  olong i, nplot;
  ofloat *pfreq, *pdata, ymax, ymin;
  ObitPlot *plot = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar strTemp[121];
  gchar *routine = "ObitTimeFilterPlotPower";

  /* error checks */
  g_assert (ObitTimeFilterIsA(in));
  if (err->error) return;

  /* Asked for something? */
  if (in->nFreq<=0) return;

  /* initialize plotting */
  plot = newObitPlot (routine);

  /* Plot labels */
  if (label) {
    dim[0] = strlen (label);
    ObitInfoListAlwaysPut (plot->info, "TITLE", OBIT_string, dim, label);
  } else {
    g_snprintf (strTemp, 120, "Power Spectrum, series %d", seriesNo);
    dim[0] = strlen (strTemp);
    ObitInfoListAlwaysPut (plot->info, "TITLE", OBIT_string, dim, strTemp);
  }
  strncpy (strTemp, "Freq (Hz)", 120);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "XLABEL", OBIT_string, dim, strTemp);
  strncpy (strTemp, "Power", 120);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "YLABEL", OBIT_string, dim, strTemp);

  /* initialize plotting */
  ObitPlotInitPlot (plot, NULL, 15,1,1, err);
  if (err->error) {
    plot = ObitPlotUnref(plot);
    Obit_traceback_msg (err, routine, plot->name);
  }

  /* Allocate arrays */
  pfreq   = g_malloc0(in->nFreq*sizeof(ofloat));
  pdata   = g_malloc0(in->nFreq*sizeof(ofloat));

  /* Power spectrum to plot */
  nplot = 0;
  ymax = -1.0e20;
  ymin =  1.0e20;
  for (i=0; i<in->nFreq; i++) {
    pfreq[nplot] = i*in->dFreq;
    pdata[nplot] = sqrt(in->freqData[seriesNo][i*2]*in->freqData[seriesNo][i*2] +
			in->freqData[seriesNo][i*2+1]*in->freqData[seriesNo][i*2+1]);
    /* Find extrema */
    ymax = MAX (ymax, pdata[nplot]);
    ymin = MIN (ymin, pdata[nplot]);
   nplot++;
  }

  /* Set extrema */
  dim[0] = 1;
  ObitInfoListAlwaysPut(plot->info, "YMAX", OBIT_float, dim, (gpointer*)&ymax);
  ObitInfoListAlwaysPut(plot->info, "YMIN", OBIT_float, dim, (gpointer*)&ymin);
  ymin = 0.0; ymax = in->nFreq*in->dFreq;
  ObitInfoListAlwaysPut(plot->info, "XMAX", OBIT_float, dim, (gpointer*)&ymax);
  ObitInfoListAlwaysPut(plot->info, "XMIN", OBIT_float, dim, (gpointer*)&ymin);

  /* plot it */
  ObitPlotXYPlot (plot, -1, nplot, pfreq, pdata, err);
  ObitPlotFinishPlot (plot, err);
  if (err->error) {
    plot = ObitPlotUnref(plot);
    if (pfreq)  g_free(pfreq);
    if (pdata)  g_free(pdata);
    Obit_traceback_msg (err, routine, routine);
  }

  /* Deallocate arrays */
  plot = ObitPlotUnref(plot);
  if (pfreq)  g_free(pfreq);
  if (pdata)  g_free(pdata);
} /* end ObitTimeFilterPlotPower */

/**
 * Plot a selected time series
 * 
 * \param in       ObitTimeFilter to plot
 * \param seriesNo Which time series to plot (0-rel)
 * \param label    If nonNULL, a label for the plot
 * \param err      Error stack, returns if not empty.
 */
void ObitTimeFilterPlotTime (ObitTimeFilter *in, olong seriesNo,
			     gchar *label, ObitErr *err)
{
  olong i, half, nplot;
  ofloat *ptime, *pdata, ymax, ymin, t0, fblank = ObitMagicF();
  ObitPlot *plot = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar strTemp[121];
  gchar *routine = "ObitTimeFilterPlotTime";

  /* error checks */
  g_assert (ObitTimeFilterIsA(in));
  if (err->error) return;

  /* Asked for something? */
  if (in->nTime<=0) return;

  /* initialize plotting */
  plot = newObitPlot (routine);

  /* Plot labels */
  if (label) {
    dim[0] = strlen (label);
    ObitInfoListAlwaysPut (plot->info, "TITLE", OBIT_string, dim, label);
  } else {
    g_snprintf (strTemp, 120, "Filter time series %d", seriesNo);
    dim[0] = strlen (strTemp);
    ObitInfoListAlwaysPut (plot->info, "TITLE", OBIT_string, dim, strTemp);
  }
  strncpy (strTemp, "Time (seconds)", 120);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "XLABEL", OBIT_string, dim, strTemp);
  strncpy (strTemp, "Value", 120);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "YLABEL", OBIT_string, dim, strTemp);

  /* initialize plotting */
  ObitPlotInitPlot (plot, NULL, 15,1,1, err);
  if (err->error) {
    plot = ObitPlotUnref(plot);
    Obit_traceback_msg (err, routine, plot->name);
  }

  /* Allocate arrays */
  ptime   = g_malloc0(in->nTime*sizeof(ofloat));
  pdata   = g_malloc0(in->nTime*sizeof(ofloat));

  /* Time series to plot */
  nplot = 0;
  ymax = -1.0e20;
  ymin =  1.0e20;
  half = in-> nTime/2;
  t0   = in->times[half];
  /* Data stored as center at the edges */
  for (i=0; i<half; i++) {
    if (in->timeData[seriesNo][i+half]!=fblank) {
      ptime[nplot] = (in->times[i+half]-t0)*86400.0;
      pdata[nplot] = in->timeData[seriesNo][i+half];
      /* Find extrema */
      ymax = MAX (ymax, pdata[nplot]);
      ymin = MIN (ymin, pdata[nplot]);
      nplot++;
    }
  }
  for (i=half; i<in->nTime; i++) {
    if (in->timeData[seriesNo][i-half]!=fblank) {
      ptime[nplot] = (in->times[i-half]-t0)*86400.0;
      pdata[nplot] = in->timeData[seriesNo][i-half];
      /* Find extrema */
      ymax = MAX (ymax, pdata[nplot]);
      ymin = MIN (ymin, pdata[nplot]);
      nplot++;
    }
  }

  /* Set extrema */
  dim[0] = 1;
  ObitInfoListAlwaysPut(plot->info, "YMAX", OBIT_float, dim, (gpointer*)&ymax);
  ObitInfoListAlwaysPut(plot->info, "YMIN", OBIT_float, dim, (gpointer*)&ymin);
  ymin = 0.0; ymax = in->nTime*in->dTime*86400.0;
  ObitInfoListAlwaysPut(plot->info, "XMAX", OBIT_float, dim, (gpointer*)&ymax);
  ObitInfoListAlwaysPut(plot->info, "XMIN", OBIT_float, dim, (gpointer*)&ymin);

  /* plot it */
  ObitPlotXYPlot (plot, -1, nplot, ptime, pdata, err);
  ObitPlotFinishPlot (plot, err);
  if (err->error) {
    plot = ObitPlotUnref(plot);
    if (ptime)  g_free(ptime);
    if (pdata)  g_free(pdata);
    Obit_traceback_msg (err, routine, routine);
  }

  /* Deallocate arrays */
  plot = ObitPlotUnref(plot);
  if (ptime)  g_free(ptime);
  if (pdata)  g_free(pdata);
} /* end ObitTimeFilterPlotTime */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitTimeFilterClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitTimeFilterClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitTimeFilterClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitTimeFilterClassInfoDefFn (gpointer inClass)
{
  ObitTimeFilterClassInfo *theClass = (ObitTimeFilterClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitTimeFilterClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTimeFilterClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitTimeFilterGetClass;
  theClass->newObit       = (newObitFP)newObitTimeFilter;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitTimeFilterClear;
  theClass->ObitInit      = (ObitInitFP)ObitTimeFilterInit;
  /* New to this class */
  theClass->newObitTimeFilter  = 
    (newObitTimeFilterFP)newObitTimeFilter;
  theClass->ObitTimeFilter2Freq  = 
    (ObitTimeFilter2FreqFP)ObitTimeFilter2Freq;
  theClass->ObitTimeFilter2Time  = 
    (ObitTimeFilter2TimeFP)ObitTimeFilter2Time;
  theClass->ObitTimeFilterFilter  = 
    (ObitTimeFilterFilterFP)ObitTimeFilterFilter;

  theClass->ObitTimeFilterPlotPower  = 
    (ObitTimeFilterPlotPowerFP)ObitTimeFilterPlotPower;
  theClass->ObitTimeFilterPlotTime  = 
    (ObitTimeFilterPlotTimeFP)ObitTimeFilterPlotTime;
} /* end ObitTimeFilterClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitTimeFilterInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTimeFilter *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread     = newObitThread();
  in->nSeries    = 0;
  in->nTime      = 0;
  in->nFreq      = 0;
  in->times      = NULL;
  in->freqs      = NULL;
  in->timeData   = NULL;
  in->freqData   = NULL;
  in->timeSeries = NULL;
  in->freqSeries = NULL;
  in->FFTFor     = NULL;
  in->FFTRev     = NULL;
  in->interp     = NULL;

} /* end ObitTimeFilterInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *         Actually it should be an ObitTimeFilter* cast to an Obit*.
 */
void ObitTimeFilterClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTimeFilter *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->FFTFor = ObitFFTUnref(in->FFTFor);
  in->FFTRev = ObitFFTUnref(in->FFTRev);
  in->interp = ObitFInterpolateUnref(in->interp);
  if (in->times) g_free(in->times);
  if (in->freqs) g_free(in->freqs);
  for (i=0; i < in->nSeries; i++) {
    if (in->timeSeries) in->timeSeries[i] = ObitFArrayUnref(in->timeSeries[i]);
    if (in->freqSeries) in->freqSeries[i] = ObitCArrayUnref(in->freqSeries[i]);
  }
  if (in->timeData) g_free(in->timeData);
  if (in->freqData) g_free(in->freqData);
  if (in->timeSeries) g_free(in->timeSeries);
  if (in->freqSeries) g_free(in->freqSeries);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitTimeFilterClear */


/**
 * Replace blanks in an array by interpolation.
 * Valid values are extended into blanked regions that cannot be interpolated.
 * \param in    Object with Interpolation structures.
 * \param array FArray to deblank.
 */
static void  ObitTimeFilterDeblank (ObitTimeFilter *in, ObitFArray *array)
{
  ofloat *data, good, value, first, last, w1, w2, iDelta;
  ofloat fblank = ObitMagicF();
  olong i, n, iFirst, iLast, pos[1] = {0};
  gboolean blanked = FALSE;

  /* error checks */
  g_assert (ObitTimeFilterIsA(in));
  g_assert (ObitFArrayIsA(array));

  /* Get array pointer */
  data = ObitFArrayIndex (array, pos);
  n = array->naxis[0];

  /* Is there anything to do? */
  good = fblank;
  iFirst = -1;
  first = data[0];
  iLast = -1;
  last = data[n-1];
  for (i=0; i<n; i++) {
    blanked = blanked || (data[i] == fblank);
    if (data[i] != fblank) good = data[i];
    /* Find first and last good points */
    if ((iFirst<0) && (data[i]!=fblank)) {
      iFirst = i;
      first = data[i];
    }
    if (data[i]!=fblank) {
      iLast = i;
      last = data[i];
    }
  }
  if (!blanked) return;

  /* If it's all blanked, zero fill and return */
  if (good==fblank) {
    for (i=0; i<n; i++) data[i] = 0.0;
    return;
  }

  /* replace array on interpolator */
  ObitFInterpolateReplace (in->interp, array);

  /* Working value needed for smooth wrap around */
  iDelta = 1.0 / ((ofloat)MAX(1, (n-iLast+iFirst)));

  /* interpolate if possible */
  for (i=0; i<n; i++) {
    if (data[i]==fblank) { /* need to interpolate or end wrap? */
      if (i<iFirst) { /* beginning - enforce smooth transition */
	w1 = (iFirst - i) * iDelta;
	w2 = 1.0 - w1;
	value = w1 * last + w2 * first;
      } else if (i>iLast) { /* end - enforce smooth transition */
	w2 = (i-iLast) * iDelta;
	w1 = 1.0 - w2;
	value = w1 * last + w2 * first;
      } else { /* in middle - interpolate */
	value = ObitFInterpolate1D (in->interp, (ofloat)(i+1.0));
      }

      /* if it's good - use it */
      if (value!=fblank) { /* good */
	data[i] = value;
	good = value;
      } else { /* bad - use last good value */
	data[i] = good;
      }
    } else { /* this one OK on input, save as last good */
      good = data[i];
    }

  } /* end loop interpolating */

} /* end ObitTimeFilterDeblank */

