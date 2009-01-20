/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitInfoList type                      */
/*                                                                    */
/*   Copyright (C) 2004-2008                                          */
/*   Associated Universities, Inc. Washington DC, USA.                */
/*                                                                    */
/*   This program is free software; you can redistribute it and/or    */
/*   modify it under the terms of the GNU General Public License as   */
/*   published by the Free Software Foundation; either version 2 of   */
/*   the License, or (at your option) any later version.              */
/*                                                                    */
/*   This program is distributed in the hope that it will be useful,  */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of   */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    */
/*   GNU General Public License for more details.                     */
/*                                                                    */
/*   You should have received a copy of the GNU General Public        */
/*   License along with this program; if not, write to the Free       */
/*   Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,     */
/*   MA 02139, USA.                                                   */
/*                                                                    */
/* Correspondence about this software should be addressed as follows: */
/*          Internet email: bcotton@nrao.edu.                         */
/*          Postal address: William Cotton                            */
/*                          National Radio Astronomy Observatory      */
/*                          520 Edgemont Road                         */
/*                          Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{
#include "ObitInfoList.h"
%}

// This cleans up the InfoListBlob we malloc'd before the function call
%typemap(python,freearg) InfoListBlob * {
   InfoListBlob *me = (InfoListBlob *)$source;
   if (me!=NULL) {
      if (me->name) free(me->name); me->name = NULL;
      if (me->data) free(me->data); me->data = NULL;
      free(me);
   }
}

// This tells SWIG to treat a InfoListBlob * argument with name 'outBlob' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.
// A type of -999 indicates failure

%typemap(python,argout) InfoListBlob *outBlob {
  PyObject *o, *t;
  int i, j, count;
  char tstring[1000];

 // Check for failure
 if ($source->type==-999) {
     o = PyString_FromString("Item not Found");
     PyList_Append($target,o);
     Py_XDECREF(o);
     return NULL;
 }
 o = PyString_FromString($source->name);
  if ((!$target) || ($target == Py_None)) {
      $target = PyList_New(0);
      PyList_Append($target,o);
      Py_XDECREF(o);
  } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    PyList_Append($target,o);
    Py_XDECREF(o);
  }

  // Type
    o = PyInt_FromLong((long)$source->type);
    PyList_Append($target,o);
    Py_XDECREF(o);

  // Dim (5)
    t = PyList_New(0);
    count = 1;
    for (i=0; i<5; i++) {
      if ($source->dim[i]>0) count *= $source->dim[i];
      o = PyInt_FromLong((long)$source->dim[i]);
      PyList_Append(t, o);
      Py_XDECREF(o);
    }
    PyList_Append($target,t);
    Py_XDECREF(t);

    if ($source->type == OBIT_string)
      count /= $source->dim[0];  /* now number of strings rather than total */

   // Data
    t = PyList_New(0);
    for (i=0; i<count; i++) {
      // a few for now 
      switch ($source->type) {
      case   OBIT_int: 
	o = PyInt_FromLong((long)((int*)$source->data)[i]);
	break;
      case   OBIT_oint:
	o = PyInt_FromLong((long)((oint*)$source->data)[i]);
        break;
      case   OBIT_long:  
	o = PyInt_FromLong(((long*)$source->data)[i]);
        break;
      case   OBIT_float: 
	o = PyFloat_FromDouble((double)((float*)$source->data)[i]);
        break;
      case   OBIT_string:
        for (j=0; j<$source->dim[0]; j++) 
	tstring[j] = ((char*)$source->data)[i*$source->dim[0]+j];
	tstring[j] = 0;
	o = PyString_FromString(tstring);
        break;
      case   OBIT_bool:
	o = PyInt_FromLong((long)((int*)$source->data)[i]);
        break;
      case   OBIT_double:
	o = PyFloat_FromDouble(((double*)$source->data)[i]);
        break;
      default:
        PyErr_SetString(PyExc_TypeError,"Unknown InfoList data type");
        return NULL;
      }; // end switch
      PyList_Append(t, o);
      Py_XDECREF(o);
    }
    PyList_Append($target,t);
    Py_XDECREF(t);
}



// Functions not implemented
//void ObitInfoListPut(ObitInfoList *in, 
//		      char* name, ObitInfoType type, gint32 *dim, 
//		      gconstpointer data, ObitErr *err);
//void ObitInfoListAlwaysPut(ObitInfoList *in, 
//		      char* name, ObitInfoType type, gint32 *dim, 
//		      gconstpointer data);
//int ObitInfoListInfo(ObitInfoList *in, 
//		     char *name, ObitInfoType *type, gint32 *dim, 
//		      ObitErr *err);
//int ObitInfoListGet(ObitInfoList *in, 
//		      char *name, ObitInfoType *type, gint32 *dim, 
//		      gpointer data, ObitErr *err);
//int ObitInfoListGetP(ObitInfoList *in, 
//		     char *name, ObitInfoType *type, gint32 *dim, 
//		      pointer *data);
//int ObitInfoListGetTest(ObitInfoList *in, 
//		      char *name, ObitInfoType *type, gint32 *dim, 
//		      pointer data);
//int ObitInfoListGetNumber(ObitInfoList *in, olong number,
//		      char *name, ObitInfoType *type, gint32 *dim, 
//		      pointer data, ObitErr *err);
//int ObitInfoListGetNumberP(ObitInfoList *in, olong number,
//		       char *name, ObitInfoType *type, gint32 *dim, 
//		       pointer *data);


%inline %{
// Structure for InfoList entry
typedef struct {
  char *name;
  int  type;
  long dim[5];
  void *data;
} InfoListBlob;

extern InfoListBlob* makeInfoListBlob() {
   InfoListBlob *out;
   out = malloc(sizeof(InfoListBlob));
   out->name = NULL;
   out->data = NULL;
   return out;
}

extern void freeInfoListBlob(InfoListBlob* blob) {
   if (!blob) return;
   if (blob->name) free(blob->name);
   if (blob->data) free(blob->data);
   free (blob);
}

ObitInfoList* InfoListCreate (void) {
  return newObitInfoList();
} // end InfoListCreate

ObitInfoList* freeInfoList (ObitInfoList *in) {
   return freeObitInfoList(in);
} // end freeList

ObitInfoList* InfoListCopy (ObitInfoList* in) {
  return ObitInfoListCopy (in);
} // end InfoListCopy

ObitInfoList* InfoListRef (ObitInfoList* in) {
  return ObitInfoListRef (in);
} // end InfoListRef

ObitInfoList* InfoListUnref (ObitInfoList* in) {
  if (!ObitInfoListIsA(in)) return NULL;
  return ObitInfoListUnref (in);
} // end InfoListUnref

ObitInfoList* InfoListCopyData (ObitInfoList* in, ObitInfoList* out) {
  return ObitInfoListCopyData (in, out);
} // end InfoListCopyData

void InfoListRemove (ObitInfoList *in, char *name) {
  return ObitInfoListRemove (in, name);
} // end InfoListRemove

void InfoListItemResize(ObitInfoList *in, 
			    char *name, int type, long *dim) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
  ObitInfoListResize (in, name, type, ddim);
} // end InfoListItemResize

int InfoListIsA (ObitInfoList* in) {
  return  ObitInfoListIsA(in);
} // end InfoListListIsA


// Put by type
extern void InfoListPutInt(ObitInfoList *in, char *name, long *dim, 
	             int* data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_long, ddim, (gpointer)data, err);
} // end InfoListPutInt

extern void InfoListAlwaysPutInt(ObitInfoList *in, char *name, 
	             long *dim, int* data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_long, ddim, (gpointer)data);
} // end InfoListAlwaysPutInp

extern void InfoListPutLong(ObitInfoList *in, char *name, long *dim, 
	             long* data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_long, ddim, (gpointer)data, err);
} // end InfoListPutLong

extern void InfoListAlwaysPutLong(ObitInfoList *in, char *name, 
	             long *dim, long* data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_long, ddim, (gpointer)data);
} // end InfoListAlwaysPutLong

extern void InfoListPutFloat(ObitInfoList *in, char *name, long *dim, 
	             float* data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_float, ddim, (gpointer)data, err);
} // end InfoListPutFloat

extern void InfoListAlwaysPutFloat(ObitInfoList *in, char *name, 
	             long *dim, float* data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_float, ddim, (gpointer)data);
} // end InfoListAlwaysPutFloat

extern void InfoListPutDouble(ObitInfoList *in, char *name, long *dim, 
	             double* data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_double, ddim, (gpointer)data, err);
} // end InfoListPutDouble

extern void InfoListAlwaysPutDouble(ObitInfoList *in, char *name, 
	             long *dim, double* data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_double, ddim, (gpointer)data);
} // end InfoListAlwaysPutDouble

extern void InfoListPutBoolean(ObitInfoList *in, char *name, long *dim, 
	             int* data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_bool, ddim, (gpointer)data, err);
} // end InfoListPutBoolean

extern void InfoListAlwaysPutBoolean(ObitInfoList *in, char *name, 
	             long *dim, int* data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_bool, ddim, (gpointer)data);
} // end InfoListAlwaysPutBoolean

extern void InfoListPutString(ObitInfoList *in, char *name, long *dim, 
	             char** data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_string, ddim, (gpointer)data, err);
} // end InfoListPutString

extern void InfoListAlwaysPutString(ObitInfoList *in, char *name, 
	             long *dim, char** data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_string, ddim, (gpointer)data);
} // end InfoListAlwayPutString

// single string version 
extern void InfoListPutSString(ObitInfoList *in, char *name, long *dim, 
	             char** data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_string, ddim, data[0], err);
} // end InfoListPutSString

// single string version 
extern void InfoListAlwaysPutSString(ObitInfoList *in, char *name, 
	             long *dim, char** data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_string, ddim, data[0]);
} // end InfoListAlwayPutsString

static int InfoListGetHelper(ObitInfoList *in, char *name, ObitInfoType type,
	gint32 dim[5], void *data, InfoListBlob *outBlob)
{
  gint32 i, j, count;
  float  *fdata;
  long   *ldata;
  int    *idata;
  oint   *odata;
  double *ddata;
  gboolean *bdata;
  char   *cdata, *codata;

  outBlob->name = strdup(name);
  outBlob->type = type;

  count = 1;
  for (i=0; i<5; i++) {
    outBlob->dim[i]  = dim[i];
    if (dim[i]>0) count *= dim[i];
  } 
  switch (type) {
    case   OBIT_int: 
      idata = (int*)data;
      outBlob->data = (void*)malloc(count*sizeof(int));
      for (i=0; i<count; i++) ((int*)outBlob->data)[i] = idata[i];
      break;
    case   OBIT_oint:
      odata = (oint*)data;
      outBlob->data = (void*)malloc(count*sizeof(oint));
      for (i=0; i<count; i++) ((oint*)outBlob->data)[i] = odata[i];
      break;
    case   OBIT_long:  
      ldata = (long*)data;
      outBlob->data = (void*)malloc(count*sizeof(long));
      for (i=0; i<count; i++) ((long*)outBlob->data)[i] = ldata[i];
      break;
    case   OBIT_float: 
      fdata = (float*)data;
      outBlob->data = (void*)malloc(count*sizeof(float));
      for (i=0; i<count; i++) ((float*)outBlob->data)[i] = fdata[i];
      break;
    case   OBIT_double:
      ddata = (double*)data;
      outBlob->data = (void*)malloc(count*sizeof(double));
      for (i=0; i<count; i++) ((double*)outBlob->data)[i] = ddata[i];
      break;
    case   OBIT_bool:
      bdata = (gboolean*)data;
      outBlob->data = (void*)malloc(count*sizeof(gboolean));
      for (i=0; i<count; i++) ((gboolean*)outBlob->data)[i] = bdata[i];
      break;
    case   OBIT_string:
      cdata = (char*)data;
      outBlob->data = (void*)malloc(count+1);
      codata = (char*)outBlob->data;
      count /= dim[0];  /* now number of strings rather than total */
      for (i=0; i<count; i++) {  // add nulls at end of strings
	for (j=0; j<dim[0]; j++)  *(codata++) = *(cdata++);
	*(codata++) = 0;
      }
      break;
    default:
      PyErr_SetString(PyExc_TypeError,"Unknown InfoList data type");
      return 0;
  }; // end switch

  return 0;
} // end InfoListGetHelper

// Convert InfoList data to a python list
static PyObject* Info2List(ObitInfoType type, gint32 dim[5], void *data)
{
  PyObject *outList=NULL, *outList2=NULL, *outList3=NULL, *o=NULL;
  gint32 i, j, k, ii, count[2], ndim;
  float  *fdata;
  long   *ldata;
  int    *idata;
  oint   *odata;
  double *ddata;
  gboolean *bdata;
  char   *tstring=NULL, *cdata;

  // how much data? - do three dimensions for now
  count[0] = dim[0]; ndim = 1;
  count[1] = dim[1]; if (dim[1]>1) ndim=2;
  count[2] = dim[2]; if (dim[2]>1) ndim=3;

  outList = PyList_New(count[0]);  // Output list
  ii = 0;

  switch (type) {
    case   OBIT_int: 
      outList = PyList_New(count[0]);  // Output list
      idata = (int*)data;
      for (i=0; i<count[0]; i++) {
  	if (ndim==1) {  // single dimension
          o = PyInt_FromLong((long)idata[ii++]);
          PyList_SetItem(outList, i, o);
        } else { // 2 or more dimensions
          outList2 = PyList_New(count[1]);
          for (j=0; j<count[1]; j++) {  // loop over second dimension
            if (ndim==2) {  // two dimensions
              o = PyInt_FromLong((long)idata[ii++]);
              PyList_SetItem(outList2, j, o);
            } else { // 3 (or more) dimensions
              outList3 = PyList_New(count[2]);
              for (k=0; k<count[2]; k++) {  // loop over third dimension
                o = PyInt_FromLong((long)idata[ii++]);
                PyList_SetItem(outList3, k, o);
              } // end loop over third dimension
              PyList_SetItem(outList2, j, outList3);
            }  
          }  // end loop over second dimension
          PyList_SetItem(outList, i, outList2);
        }
      } // end loop over first dimension
      break;
    case   OBIT_oint:
      outList = PyList_New(count[0]);  // Output list
      odata = (oint*)data;
      for (i=0; i<count[0]; i++) {
  	if (ndim==1) {  // single dimension
          o = PyInt_FromLong((long)odata[ii++]);
          PyList_SetItem(outList, i, o);
        } else { // 2 or more dimensions
          outList2 = PyList_New(count[1]);
          for (j=0; j<count[1]; j++) {  // loop over second dimension
            if (ndim==2) {  // two dimensions
              o = PyInt_FromLong((long)odata[ii++]);
              PyList_SetItem(outList2, j, o);
            } else { // 3 (or more) dimensions
              outList3 = PyList_New(count[2]);
              for (k=0; k<count[2]; k++) {  // loop over third dimension
                o = PyInt_FromLong((long)odata[ii++]);
                PyList_SetItem(outList3, k, o);
              } // end loop over third dimension
              PyList_SetItem(outList2, j, outList3);
            }  
          }  // end loop over second dimension
          PyList_SetItem(outList, i, outList2);
        }
      } // end loop over first dimension
      break;
    case   OBIT_long:  
      outList = PyList_New(count[0]);  // Output list
      ldata = (long*)data;
      for (i=0; i<count[0]; i++) {
  	if (ndim==1) {  // single dimension
          o = PyInt_FromLong((long)ldata[ii++]);
          PyList_SetItem(outList, i, o);
        } else { // 2 or more dimensions
          outList2 = PyList_New(count[1]);
          for (j=0; j<count[1]; j++) {  // loop over second dimension
            if (ndim==2) {  // two dimensions
              o = PyInt_FromLong((long)ldata[ii++]);
              PyList_SetItem(outList2, j, o);
            } else { // 3 (or more) dimensions
              outList3 = PyList_New(count[2]);
              for (k=0; k<count[2]; k++) {  // loop over third dimension
                o = PyInt_FromLong((long)ldata[ii++]);
                PyList_SetItem(outList3, k, o);
              } // end loop over third dimension
              PyList_SetItem(outList2, j, outList3);
            }  
          }  // end loop over second dimension
          PyList_SetItem(outList, i, outList2);
        }
      } // end loop over first dimension
      break;
    case   OBIT_float: 
      outList = PyList_New(count[0]);  // Output list
      fdata = (float*)data;
      for (i=0; i<count[0]; i++) {
  	if (ndim==1) {  // single dimension
          o = PyFloat_FromDouble((double)fdata[ii++]);
          PyList_SetItem(outList, i, o);
        } else { // 2 or more dimensions
          outList2 = PyList_New(count[1]);
          for (j=0; j<count[1]; j++) {  // loop over second dimension
            if (ndim==2) {  // two dimensions
              o = PyFloat_FromDouble((double)fdata[ii++]);
              PyList_SetItem(outList2, j, o);
            } else { // 3 (or more) dimensions
              outList3 = PyList_New(count[2]);
              for (k=0; k<count[2]; k++) {  // loop over third dimension
                o = PyFloat_FromDouble((double)fdata[ii++]);
                PyList_SetItem(outList3, k, o);
              } // end loop over third dimension
              PyList_SetItem(outList2, j, outList3);
            }  
          }  // end loop over second dimension
          PyList_SetItem(outList, i, outList2);
        }
      } // end loop over first dimension
      break;
    case   OBIT_double:
      outList = PyList_New(count[0]);  // Output list
      ddata = (double*)data;
      for (i=0; i<count[0]; i++) {
  	if (ndim==1) {  // single dimension
          o = PyFloat_FromDouble((double)ddata[ii++]);
          PyList_SetItem(outList, i, o);
        } else { // 2 or more dimensions
          outList2 = PyList_New(count[1]);
          for (j=0; j<count[1]; j++) {  // loop over second dimension
            if (ndim==2) {  // two dimensions
              o = PyFloat_FromDouble((double)ddata[ii++]);
              PyList_SetItem(outList2, j, o);
            } else { // 3 (or more) dimensions
              outList3 = PyList_New(count[2]);
              for (k=0; k<count[2]; k++) {  // loop over third dimension
                o = PyFloat_FromDouble((double)ddata[ii++]);
                PyList_SetItem(outList3, k, o);
              } // end loop over third dimension
              PyList_SetItem(outList2, j, outList3);
            }  
          }  // end loop over second dimension
          PyList_SetItem(outList, i, outList2);
        }
      } // end loop over first dimension
      break;
    case   OBIT_bool:
      outList = PyList_New(count[0]);  // Output list
      bdata = (gboolean*)data;
      for (i=0; i<count[0]; i++) {
  	if (ndim==1) {  // single dimension
          o = PyInt_FromLong((long)bdata[ii++]);
          PyList_SetItem(outList, i, o);
        } else { // 2 or more dimensions
          outList2 = PyList_New(count[1]);
          for (j=0; j<count[1]; j++) {  // loop over second dimension
            if (ndim==2) {  // two dimensions
              o = PyInt_FromLong((long)bdata[ii++]);
              PyList_SetItem(outList2, j, o);
            } else { // 3 (or more) dimensions
              outList3 = PyList_New(count[2]);
              for (k=0; k<count[2]; k++) {  // loop over third dimension
                o = PyInt_FromLong((long)bdata[ii++]);
                PyList_SetItem(outList3, k, o);
              } // end loop over third dimension
              PyList_SetItem(outList2, j, outList3);
            }  
          }  // end loop over second dimension
          PyList_SetItem(outList, i, outList2);
        }
      } // end loop over first dimension
      break;
    case   OBIT_string:
      cdata = (char*)data;
      outList = PyList_New(MAX (1, count[1]));  // Output list
      // First dimension is the length of strings
      tstring	= g_malloc0(count[0]+1);
      if (ndim==1) {  // single dimension - one string
        for (i=0; i<count[0]; i++) tstring[i]=cdata[ii++]; tstring[i]=0;
        o = PyString_InternFromString(tstring);
        PyList_SetItem(outList, 0, o);
      } else { // 2 or more dimensions
        for (j=0; j<count[1]; j++) {  // loop over second dimension
          if (ndim==2) {  // two dimensions
            for (i=0; i<count[0]; i++) tstring[i]=cdata[ii++]; tstring[i]=0;
            o = PyString_InternFromString(tstring);
            PyList_SetItem(outList, j, o);
          } else { // 3 (or more) dimensions
            outList2 = PyList_New(count[2]);
            for (k=0; k<count[2]; k++) {  // loop over third dimension
              for (i=0; i<count[0]; i++) tstring[i]=cdata[ii++]; tstring[i]=0;
              o = PyString_InternFromString(tstring);
              PyList_SetItem(outList2, k, o);
            } // end loop over third dimension
            PyList_SetItem(outList, j, outList2);
          }  
        }  // end loop over second dimension
      } // end multiple dimensions
      if (tstring) g_free(tstring);
      break;
    default:
      PyErr_SetString(PyExc_TypeError,"Unknown InfoList data type");
      return outList;
  }; // end switch

  return outList;
} // end Info2list

// Retrieve value as python list
//       0 - return code, 0=OK else failed
//       1 - name
//       2 - type
//       3 - dimension array
//       4 - data array
extern PyObject* InfoListGet(ObitInfoList *in, char *name) {
  gint32 i, dim[5];
  ObitInfoType type;
  void  *data;
  PyObject *dlist=NULL;
  PyObject *outList=NULL, *tl, *o;

  // output list
  outList = PyList_New(5);

  if (ObitInfoListGetP(in, name, &type, dim, &data)) {
    dlist = Info2List(type, dim, data);
    o = PyInt_FromLong((long)0);
    PyList_SetItem(outList, 0, o);
    o = PyString_InternFromString(name);
    PyList_SetItem(outList, 1, o);	
    o = PyInt_FromLong((long)type);
    PyList_SetItem(outList, 2, o);
    tl = PyList_New(5);
    for (i=0; i<5; i++) {
      o = PyInt_FromLong((long)dim[i]);
      PyList_SetItem(tl, i, o);
    }	
    PyList_SetItem(outList, 3, tl);	
    PyList_SetItem(outList, 4, dlist);	
    return outList;
   } /* end worked OK */

  // failed 
  //PyErr_SetString(PyExc_TypeError,"Item not found");
  o = PyInt_FromLong((long)999);
  PyList_SetItem(outList, 0, o);
  o = PyString_InternFromString(name);
  PyList_SetItem(outList, 1, o);	
  PyList_SetItem(outList, 2, NULL);	
  PyList_SetItem(outList, 3, NULL);	
  PyList_SetItem(outList, 4, NULL);	
  return outList;
} // end InfoListGet

// Retrieve numbered value as python list
//       0 - return code, 0=OK else failed
//       1 - name
//       2 - type
//       3 - dimension array
//       4 - data array
extern PyObject* InfoListGetNumber(ObitInfoList *in, int number)
{
  gint32 i, dim[5];
  char *name;
  ObitInfoType type;
  void  *data;
  PyObject *dlist=NULL;
  PyObject *outList=NULL, *tl, *o;

  // output list
  outList = PyList_New(5);

  if (ObitInfoListGetNumberP(in, number, &name, &type, dim, &data)) {
    dlist = Info2List(type, dim, data);
    o = PyInt_FromLong((long)0);
    PyList_SetItem(outList, 0, o);
    o = PyString_InternFromString(name);
    PyList_SetItem(outList, 1, o);	
    o = PyInt_FromLong((long)type);
    PyList_SetItem(outList, 2, o);
    tl = PyList_New(5);
    for (i=0; i<5; i++) {
      o = PyInt_FromLong((long)dim[i]);
      PyList_SetItem(tl, i, o);
    }	
    PyList_SetItem(outList, 3, tl);	
    PyList_SetItem(outList, 4, dlist);	
    return outList;
   } /* end worked OK */

  // failed 
  PyErr_SetString(PyExc_TypeError,"Item not found");
  o = PyInt_FromLong((long)999);
  PyList_SetItem(outList, 0, o);
  o = PyString_InternFromString(name);
  PyList_SetItem(outList, 1, o);	
  PyList_SetItem(outList, 2, NULL);	
  PyList_SetItem(outList, 3, NULL);	
  PyList_SetItem(outList, 4, NULL);	
  return outList;
} // end InfoListGetNumber

// Retrieve whole infolist as a python dict
// each entry as python list
//       0 - type
//       1 - dimension array
//       2 - data array
extern PyObject* InfoListGetDict(ObitInfoList *in)
{
  PyObject *outDict=PyDict_New();
  gint32 i, iel, dim[5];
  char *name;
  ObitInfoType type;
  void  *data;
  PyObject *dlist=NULL;
  PyObject *tList=NULL, *tl, *o;

  // loop over elements
  for (iel = 1; iel<=in->number; iel++) {
    // output list
    tList = PyList_New(3);
  
    if (ObitInfoListGetNumberP(in, iel, &name, &type, dim, &data)) {
      dlist = Info2List(type, dim, data);
      o = PyInt_FromLong((long)type);
      PyList_SetItem(tList, 0, o);
      tl = PyList_New(5);
      for (i=0; i<5; i++) {
        o = PyInt_FromLong((long)dim[i]);
        PyList_SetItem(tl, i, o);
      }	
      PyList_SetItem(tList, 1, tl);	
      PyList_SetItem(tList, 2, dlist);	
      // add to output
      PyDict_SetItemString(outDict, name, tList);
    } /* end worked OK */
  }	  // end element loop
  return outDict;
} // end InfoListGetDict

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitInfoList *me;
} InfoList;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitInfoList *me;
} InfoList;

%addmethods InfoList { 
  InfoList(void) {
     InfoList *out;
     out = (InfoList *) malloc(sizeof(InfoList));
     out->me = InfoListCreate();
    return out;
   }
  ~InfoList() {
    self->me = InfoListUnref(self->me);
    free(self);
  }
};

