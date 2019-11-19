/* $Id$   */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageDesc type                         */
/*                                                                    */
/*;  Copyright (C)2005-2019                                           */
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

%{
#include "ObitUVDesc.h"
#include "ObitPrecess.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUVDesc *me;
} UVDesc;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUVDesc *me;
} UVDesc;

%inline %{
extern ObitUVDesc* UVDescCreate (char *name) {
  return newObitUVDesc (name);
} // end UVDescCreate

ObitUVDesc* UVDescRef (ObitUVDesc* in) {
  return ObitUVDescRef (in);
} // end UVDescRef

ObitUVDesc* UVDescUnref (ObitUVDesc* in) {
  if (!ObitUVDescIsA(in)) return NULL;
  return ObitUVDescUnref (in);
} // end UVDescUnref

void CreateUVDesc(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVDesc, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateUVDesc: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((UVDesc*)ptr)->me = UVDescCreate(name);
   else  ((UVDesc*)ptr)->me = NULL;
   return;
 } // end CreateUVDesc

void DeleteUVDesc(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVDesc, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteUVDesc: could not recover c struct");
	return;
   }
   ((UVDesc*)ptr)->me = ObitUVDescUnref(((UVDesc*)ptr)->me);
}// end DeleteUVDesc

ObitUVDesc* UVDesc_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVDesc, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"UVDesc_Get_me: could not recover c struct");
	return NULL;
   }
   return ((UVDesc*)ptr)->me;
}  // end UVDesc_Get_me

void UVDesc_Set_me(PyObject *self, ObitUVDesc *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVDesc, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"UVDesc_Set_me: could not recover c struct");
	return;
   }
   ((UVDesc*)ptr)->me =  UVDescRef(me);
}  // end UVDesc_Set_me

extern ObitUVDesc* UVDescCopy (ObitUVDesc* in, 
		              ObitUVDesc* out, ObitErr *err) {
  return ObitUVDescCopy (in, out, err);
} // end UVDescCopy

extern void UVDescCopyDesc (ObitUVDesc* in, ObitUVDesc* out,
			ObitErr *err) {
  ObitUVDescCopyDesc  (in, out, err);
} // end UVDescCopyDesc

extern void UVDescIndex (ObitUVDesc* in) {
  ObitUVDescIndex (in);
} // end UVDescIndex

extern ObitInfoList* UVDescGetList (ObitUVDesc* in) {
  return ObitInfoListRef(in->info);
}
 
extern double UVDescDate2JD (char *date) {
  double JD;
  ObitUVDescDate2JD ((const gchar*)date, &JD);
  return JD;
}
 
extern double UVDescGST0 (double JD) {
  double GSTUTC0, Rate;
  ObitPrecessGST0 (JD, &GSTUTC0, &Rate);
  return GSTUTC0;
}
 
extern double UVDescERate (double JD) {
  double GSTUTC0, Rate;
  ObitPrecessGST0 (JD, &GSTUTC0, &Rate);
  return Rate;
}
 
extern PyObject *UVDescGetDict(ObitUVDesc* in) {
  PyObject *outDict = PyDict_New();
  PyObject *list1, *list2, *list3, *list4, *list5, *list6, *list7;
  PyObject *value;
  Py_ssize_t doh;
  int i, pos = 0;

  PyDict_SetItemString(outDict, "name",    PyString_InternFromString(in->name));
  PyDict_SetItemString(outDict, "object",  PyString_InternFromString(in->object));
  PyDict_SetItemString(outDict, "teles",   PyString_InternFromString(in->teles));
  PyDict_SetItemString(outDict, "instrume",PyString_InternFromString(in->instrument));
  PyDict_SetItemString(outDict, "observer",PyString_InternFromString(in->observer));
  PyDict_SetItemString(outDict, "origin",  PyString_InternFromString(in->origin));
  PyDict_SetItemString(outDict, "date",    PyString_InternFromString(in->date));
  PyDict_SetItemString(outDict, "obsdat",  PyString_InternFromString(in->obsdat));
  PyDict_SetItemString(outDict, "isort",   PyString_InternFromString(in->isort));
  PyDict_SetItemString(outDict, "bunit",   PyString_InternFromString(in->bunit));
  PyDict_SetItemString(outDict, "obsra",   PyFloat_FromDouble((double)in->obsra));
  PyDict_SetItemString(outDict, "obsdec",  PyFloat_FromDouble((double)in->obsdec));
  PyDict_SetItemString(outDict, "epoch",   PyFloat_FromDouble((double)in->epoch));
  PyDict_SetItemString(outDict, "equinox", PyFloat_FromDouble((double)in->equinox));
  PyDict_SetItemString(outDict, "beamMaj", PyFloat_FromDouble((double)in->beamMaj));
  PyDict_SetItemString(outDict, "beamMin", PyFloat_FromDouble((double)in->beamMin));
  PyDict_SetItemString(outDict, "beamPA",  PyFloat_FromDouble((double)in->beamPA));
  PyDict_SetItemString(outDict, "xshift",  PyFloat_FromDouble((double)in->xshift));
  PyDict_SetItemString(outDict, "yshift",  PyFloat_FromDouble((double)in->yshift));
  PyDict_SetItemString(outDict, "altCrpix",PyFloat_FromDouble((double)in->altCrpix));
  PyDict_SetItemString(outDict, "altRef",  PyFloat_FromDouble((double)in->altRef));
  PyDict_SetItemString(outDict, "restFreq",PyFloat_FromDouble((double)in->restFreq));
  PyDict_SetItemString(outDict, "JDObs",   PyFloat_FromDouble((double)in->JDObs));
  PyDict_SetItemString(outDict, "naxis",   PyLong_FromLong((long)in->naxis));
  PyDict_SetItemString(outDict, "nvis",    PyLong_FromLong((long)in->nvis));
  PyDict_SetItemString(outDict, "nrparm",  PyLong_FromLong((long)in->nrparm));
  PyDict_SetItemString(outDict, "lrec",    PyLong_FromLong((long)in->lrec));
  PyDict_SetItemString(outDict, "ncorr",   PyLong_FromLong((long)in->ncorr));
  PyDict_SetItemString(outDict, "VelReference", PyLong_FromLong((long)in->VelReference));
  PyDict_SetItemString(outDict, "VelDef", PyLong_FromLong((long)in->VelDef));
  list1 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list1, i, PyString_InternFromString(in->ctype[i]));
  PyDict_SetItemString(outDict, "ctype", list1);
  list2 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list2, i, PyFloat_FromDouble((double)in->crval[i]));
  PyDict_SetItemString(outDict, "crval", list2);
  list3 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list3, i, PyLong_FromLong((long)in->inaxes[i]));
  PyDict_SetItemString(outDict, "inaxes", list3);
  list4 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list4, i, PyFloat_FromDouble((double)in->cdelt[i]));
  PyDict_SetItemString(outDict, "cdelt", list4);
  list5 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list5, i, PyFloat_FromDouble((double)in->crpix[i]));
  PyDict_SetItemString(outDict, "crpix", list5);
  list6 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list6, i, PyFloat_FromDouble((double)in->crota[i]));
  PyDict_SetItemString(outDict, "crota", list6);
  list7 = PyList_New(UV_MAX_RANP);
  for (i=0; i<UV_MAX_RANP; i++) PyList_SetItem(list7, i, PyString_InternFromString(in->ptype[i]));
  PyDict_SetItemString(outDict, "ptype", list7);
  PyDict_SetItemString(outDict, "firstVis", PyLong_FromLong((long)in->firstVis));
  PyDict_SetItemString(outDict, "numVisBuff", PyLong_FromLong((long)in->numVisBuff));
/* Structural members - mostly read only */
  PyDict_SetItemString(outDict, "ilocu",    PyLong_FromLong((long)in->ilocu));
  PyDict_SetItemString(outDict, "ilocv",    PyLong_FromLong((long)in->ilocv));
  PyDict_SetItemString(outDict, "ilocw",    PyLong_FromLong((long)in->ilocw));
  PyDict_SetItemString(outDict, "iloct",    PyLong_FromLong((long)in->iloct));
  PyDict_SetItemString(outDict, "ilocb",    PyLong_FromLong((long)in->ilocb));
  PyDict_SetItemString(outDict, "ilocsu",   PyLong_FromLong((long)in->ilocsu));
  PyDict_SetItemString(outDict, "ilocfq",   PyLong_FromLong((long)in->ilocfq));
  PyDict_SetItemString(outDict, "ilocit",   PyLong_FromLong((long)in->ilocit));
  PyDict_SetItemString(outDict, "ilocid",   PyLong_FromLong((long)in->ilocid));
  PyDict_SetItemString(outDict, "ilocws",   PyLong_FromLong((long)in->ilocws));
  PyDict_SetItemString(outDict, "jlocc",    PyLong_FromLong((long)in->jlocc));
  PyDict_SetItemString(outDict, "jlocs",    PyLong_FromLong((long)in->jlocs));
  PyDict_SetItemString(outDict, "jlocf",    PyLong_FromLong((long)in->jlocf));
  PyDict_SetItemString(outDict, "jlocr",    PyLong_FromLong((long)in->jlocr));
  PyDict_SetItemString(outDict, "jlocd",    PyLong_FromLong((long)in->jlocd));
  PyDict_SetItemString(outDict, "jlocif",   PyLong_FromLong((long)in->jlocif));
  PyDict_SetItemString(outDict, "incs",     PyLong_FromLong((long)in->incs));
  PyDict_SetItemString(outDict, "incf",     PyLong_FromLong((long)in->incf));
  PyDict_SetItemString(outDict, "incif",    PyLong_FromLong((long)in->incif));
 
  /* Discard references to newly created objects. */
  doh = pos;
  while (PyDict_Next(outDict, &doh, NULL, &value))
    Py_DECREF(value);

  return outDict;
} // end UVDescGetDict

extern void UVDescSetDict(ObitUVDesc* in, PyObject *inDict) {
  PyObject *list1, *list2, *list3, *list4, *list5, *list6, *list7;
  char *tstr;
  int i, number;

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return;
  }

  tstr = PyString_AsString(PyDict_GetItemString(inDict, "object"));
  strncpy (in->object, tstr, UVLEN_VALUE); in->object[UVLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obsdat"));
  strncpy (in->obsdat, tstr, UVLEN_VALUE); in->obsdat[UVLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "bunit"));
  strncpy (in->bunit, tstr, UVLEN_VALUE); in->bunit[UVLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "teles"));
  strncpy (in->teles, tstr, UVLEN_VALUE); in->teles[UVLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "instrume"));
  strncpy (in->instrument, tstr, IMLEN_VALUE); in->instrument[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "observer"));
  strncpy (in->observer, tstr, IMLEN_VALUE); in->observer[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "origin"));
  strncpy (in->origin, tstr, UVLEN_VALUE); in->origin[UVLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "date"));
  strncpy (in->date, tstr, UVLEN_VALUE); in->date[UVLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "isort"));
  strncpy (in->isort, tstr, 2); in->isort[2]=0;
  in->epoch   = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "epoch"));
  in->equinox = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "equinox"));
  in->obsra   = PyFloat_AsDouble(PyDict_GetItemString(inDict, "obsra"));
  in->obsdec  = PyFloat_AsDouble(PyDict_GetItemString(inDict, "obsdec"));
  in->beamMaj = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "beamMaj"));
  in->beamMin = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "beamMin"));
  in->beamPA  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "beamPA"));
  in->restFreq= (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "restFreq"));
  in->JDObs   = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "JDObs"));
  in->xshift  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "xshift"));
  in->yshift  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "yshift"));
  in->altCrpix= (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "altCrpix"));
  in->altRef  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "altRef"));
  in->naxis   = (int)PyLong_AsLong(PyDict_GetItemString(inDict, "naxis"));
  in->nrparm  = PyLong_AsLong(PyDict_GetItemString(inDict, "nrparm"));
  in->nvis    = PyLong_AsLong(PyDict_GetItemString(inDict, "nvis"));
  in->numVisBuff = PyLong_AsLong(PyDict_GetItemString(inDict, "numVisBuff"));
  in->VelReference  = (int)PyLong_AsLong(PyDict_GetItemString(inDict, "VelReference"));
  in->VelDef  = (int)PyLong_AsLong(PyDict_GetItemString(inDict, "VelDef"));
  list1 = PyDict_GetItemString(inDict, "ctype");
  number = MIN (IM_MAXDIM, PyList_Size(list1));
  for (i=0; i<number; i++) {
    tstr = PyString_AsString(PyList_GetItem(list1, i));
    strncpy (in->ctype[i], tstr, UVLEN_KEYWORD);
  }
  list2 = PyDict_GetItemString(inDict, "crval");
  number = MIN (IM_MAXDIM, PyList_Size(list2));
  for (i=0; i<number; i++) in->crval[i] = PyFloat_AsDouble(PyList_GetItem(list2, i));
  list3 = PyDict_GetItemString(inDict, "inaxes");
  number = MIN (IM_MAXDIM, PyList_Size(list3));
  for (i=0; i<number; i++) in->inaxes[i] = (long)PyLong_AsLong(PyList_GetItem(list3, i));
  list4 = PyDict_GetItemString(inDict, "cdelt");
  number = MIN (IM_MAXDIM, PyList_Size(list4));
  for (i=0; i<number; i++) in->cdelt[i] = (float)PyFloat_AsDouble(PyList_GetItem(list4, i));
  list5 = PyDict_GetItemString(inDict, "crpix");
  number = MIN (IM_MAXDIM, PyList_Size(list5));
  for (i=0; i<number; i++) in->crpix[i] = (float)PyFloat_AsDouble(PyList_GetItem(list5, i));
  list6 = PyDict_GetItemString(inDict, "crota");
  number = MIN (IM_MAXDIM, PyList_Size(list6));
  for (i=0; i<number; i++) in->crota[i] = (float)PyFloat_AsDouble(PyList_GetItem(list6, i));
  list7 = PyDict_GetItemString(inDict, "ptype");
  number = MIN (UV_MAX_RANP, PyList_Size(list7));
  for (i=0; i<number; i++) {
    tstr = PyString_AsString(PyList_GetItem(list7, i));
    strncpy (in->ptype[i], tstr, UVLEN_KEYWORD);
  }
  in->firstVis = (int)PyLong_AsLong(PyDict_GetItemString(inDict, "firstVis"));
  /* Reindex just to be sure */
  ObitUVDescIndex (in);
} // end UVDescSetDict

extern int UVDescIsA (ObitUVDesc* in) {
  return ObitUVDescIsA(in);
}
%}

