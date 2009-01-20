/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitTimeFilter type                    */
/*                                                                    */
/*;  Copyright (C) 2008                                               */
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

%{
#include "ObitTimeFilter.h"
%}


%inline %{

/** Public: Create time filter. */
extern ObitTimeFilter* TimeFilterCreate(char* name, int nTime, int nSeries) {
   return  newObitTimeFilter ((gchar*)name, (olong)nTime, (olong)nSeries);
}

/** Public: Resize filter. */
extern void TimeFilterResize (ObitTimeFilter* in, int nTime) {
  ObitTimeFilterResize (in, (olong)nTime);
} // end TimeFilterResize

/** Public: Construct regular time series. */
extern void TimeFilterGridTime (ObitTimeFilter* in, int seriesNo, float dTime, int nTime, 
                        float *times, float *data) {
  ObitTimeFilterGridTime (in, (olong)seriesNo, (ofloat)dTime, (olong)nTime, 
                        (ofloat*)times, (ofloat*)data);
} // end TimeFilterGridTime

/** Public: Return time series. */
extern PyObject* TimeFilterUngridTime (ObitTimeFilter* in, int seriesNo,
			         int nTime, float *times) {
  olong i;
  PyObject *out;
  ofloat *tdata=NULL;

  tdata = g_malloc0(nTime*sizeof(ofloat));
  ObitTimeFilterUngridTime (in, (olong)seriesNo, (olong)nTime, 
                            (ofloat*)times, tdata);
  out = PyList_New(nTime);
  for (i=0; i<nTime; i++) 
    PyList_SetItem(out, i, PyFloat_FromDouble((double)tdata[i]));
  if (tdata) g_free(tdata);
  return out;
} // end TimeFilterUngridTime

/** Public: Compute frequency series. */
extern void TimeFilter2Freq (ObitTimeFilter* in) {
  ObitTimeFilter2Freq (in);
} // end TimeFilter2Freq

/** Public: Compute Time series. */
extern void TimeFilter2Time (ObitTimeFilter* in) {
  ObitTimeFilter2Time (in);
} // end TimeFilter2Time

/** Public: Apply Filter to Frequency series. */
extern void TimeFilterFilter (ObitTimeFilter* in, int seriesNo, int type, 
                              float *parms, ObitErr *err) {
  ObitTimeFilterType ltype=OBIT_TimeFilter_LowPass;
  if (type==1) ltype = OBIT_TimeFilter_HighPass;
  if (type==2) ltype = OBIT_TimeFilter_NotchPass;
  if (type==3) ltype = OBIT_TimeFilter_NotchBlock;

  ObitTimeFilterFilter (in, (olong)seriesNo, ltype, (ofloat*)parms, err);
} // end TimeFilterFilter

/** Public: Apply Filter to Frequency series with physical parameters. */
extern void TimeFilterDoFilter (ObitTimeFilter* in, int seriesNo, int type, 
                                float *parms, ObitErr *err) {
  ObitTimeFilterType ltype=OBIT_TimeFilter_LowPass;
  if (type==1) ltype = OBIT_TimeFilter_HighPass;
  if (type==2) ltype = OBIT_TimeFilter_NotchPass;
  if (type==3) ltype = OBIT_TimeFilter_NotchBlock;

  ObitTimeFilterDoFilter (in, (olong)seriesNo, ltype, (ofloat*)parms, err);
} // end TimeFilterDoFilter

/** Public:  Plot power spectrum. */
extern void TimeFilterPlotPower (ObitTimeFilter* in, int seriesNo, char *label, ObitErr *err) {
  ObitTimeFilterPlotPower (in, (olong)seriesNo, (gchar*)label, err);
} // end TimeFilter

/** Public: Plot Time series.. */
extern void TimeFilterPlotTime (ObitTimeFilter* in, int seriesNo, char *label, ObitErr *err) {
  ObitTimeFilterPlotTime (in, (olong)seriesNo, (gchar*)label, err);
} // end TimeFilter

extern int TimeFilterIsA (ObitTimeFilter* in) {
  return ObitTimeFilterIsA(in);
} // end TimeFilterIsA 

extern ObitTimeFilter* TimeFilterRef (ObitTimeFilter* in) {
  return ObitTimeFilterRef (in);
} // end TimeFilterRef

extern ObitTimeFilter* TimeFilterUnref (ObitTimeFilter* in) {
  if (!ObitTimeFilterIsA(in)) return NULL;
  return ObitTimeFilterUnref (in);
} // end TimeFilterUnref

extern char* TimeFilterGetName (ObitTimeFilter* in) {
  return in->name;
} // end  TimeFilterGetName

// return dict with {"dTime", "time", "data"}
extern PyObject* TimeFilterGetTime (ObitTimeFilter* in, int seriesNo) {
  PyObject *outDict = PyDict_New();
  PyObject *time, *data;
  olong i;

  PyDict_SetItemString(outDict, "dTime", PyFloat_FromDouble((double)in->dTime));
  time = PyList_New(in->nTime);
  data = PyList_New(in->nTime);
  for (i=0; i<in->nTime; i++) {
    PyList_SetItem(time, i, PyFloat_FromDouble((double)in->times[i]));
    PyList_SetItem(data, i, PyFloat_FromDouble((double)in->timeData[seriesNo][i]));
  }
  PyDict_SetItemString(outDict, "time", time);
  PyDict_SetItemString(outDict, "data", data);
   
  return outDict;
} // end  TimeFilterGetTime

// return dict with {"dFreq", "freq", "data"}
extern PyObject* TimeFilterGetFreq (ObitTimeFilter* in, int seriesNo) {
  PyObject *outDict = PyDict_New();
  PyObject *freq, *data;
  olong i;

  PyDict_SetItemString(outDict, "dFreq", PyFloat_FromDouble((double)in->dFreq));
  freq = PyList_New(in->nFreq);
  data = PyList_New(in->nFreq);

  for (i=0; i<in->nFreq; i++) {
    PyList_SetItem(freq, i, PyFloat_FromDouble((double)(24.0*in->freqs[i])));
    PyList_SetItem(data, i, PyComplex_FromDoubles((double)in->freqData[seriesNo][i*2], 
                                                   (double)in->freqData[seriesNo][i*2+1]));
  }
  PyDict_SetItemString(outDict, "freq", freq);
  PyDict_SetItemString(outDict, "data", data);
   
  return outDict;
} // end  TimeFilterGetFreq

// return dict with {"dFreq", "freq", "data"}
// Power spectra
extern PyObject* TimeFilterGetPower (ObitTimeFilter* in, int seriesNo) {
  PyObject *outDict = PyDict_New();
  PyObject *freq, *data;
  double val, vr, vi;
  olong i;

  PyDict_SetItemString(outDict, "dFreq", PyFloat_FromDouble((double)in->dFreq));
  freq = PyList_New(in->nFreq);
  data = PyList_New(in->nFreq);

  for (i=0; i<in->nFreq; i++) {
    PyList_SetItem(freq, i, PyFloat_FromDouble((double)(24.0*in->freqs[i])));
    vr = (double)in->freqData[seriesNo][i*2];
    vi = (double)in->freqData[seriesNo][i*2+1];
    val = vr*vr + vi*vi;
    PyList_SetItem(data, i, PyFloat_FromDouble(val));
  }
  PyDict_SetItemString(outDict, "freq", freq);
  PyDict_SetItemString(outDict, "data", data);
   
  return outDict;
} // end  TimeFilterGetPower

// expect dict with {"dTime", "time", "data"}
extern void TimeFilterSetTime (ObitTimeFilter* in, int seriesNo,
                               PyObject *inDict) {
  PyObject *time, *data;
  olong len, i;

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return;
  }

  in->dTime = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "dTime"));
  time = PyDict_GetItemString(inDict, "time");
  data = PyDict_GetItemString(inDict, "data");
  len = PyList_Size(time);
  if (len>in->nTime) {
    PyErr_SetString(PyExc_TypeError,"Data length incompatible with filter");
    return;
  }

  for (i=0; i<in->nTime; i++) {
    in->times[i]              = (ofloat)PyFloat_AsDouble(PyList_GetItem(time, i));
    in->timeData[seriesNo][i] = (ofloat)PyFloat_AsDouble(PyList_GetItem(data, i));
  }
} // end  TimeFilterSetTime

// expect dict with {"dFreq", "freq", "data"}
extern void TimeFilterSetFreq (ObitTimeFilter* in, int seriesNo,
                               PyObject *inDict) {
  PyObject *freq, *data, *cx;
  olong len, i;

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return;
  }

  in->dFreq = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "dFreq"));
  freq = PyDict_GetItemString(inDict, "freq");
  data = PyDict_GetItemString(inDict, "data");
  len = PyList_Size(freq);
  if (len>in->nFreq) {
    PyErr_SetString(PyExc_TypeError,"Data length incompatible with filter");
    return;
  }

  for (i=0; i<in->nTime; i++) {
    in->freqs[i] = (ofloat)PyFloat_AsDouble(PyList_GetItem(freq, i));
    cx = PyList_GetItem(data, i);
    in->freqData[seriesNo][i*2]   = (ofloat)PyComplex_RealAsDouble(cx);
    in->freqData[seriesNo][i*2+1] = (ofloat)PyComplex_ImagAsDouble(cx);
  }
   
} // end  TimeFilterSetFreq

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitTimeFilter *me;
} TimeFilter;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitTimeFilter *me;
} TimeFilter;

%addmethods TimeFilter { 
  TimeFilter(char* name, int nTime, int nSeries) {
     TimeFilter *out;
     out = (TimeFilter*) malloc(sizeof(TimeFilter));
     if (strcmp(name, "None")) out->me = TimeFilterCreate(name, nTime, nSeries);
     else out->me = NULL;
     return out;
   }
  ~TimeFilter() {
    self->me = TimeFilterUnref(self->me);
    free(self);
  }
};

