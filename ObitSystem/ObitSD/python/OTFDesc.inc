/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig module description for OTFDesc type                           */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitOTFDesc.h"
%}

%inline %{
ObitOTFDesc* OTFDescCreate (gchar *name) {
  return newObitOTFDesc (name);
} // end OTFDescCreate

ObitOTFDesc* OTFDescCopy (ObitOTFDesc* in, 
		              ObitOTFDesc* out, ObitErr *err) {
  return ObitOTFDescCopy (in, out, err);
} // end OTFDescCopy

void OTFDescCopyDesc (ObitOTFDesc* in, ObitOTFDesc* out,
			ObitErr *err) {
  ObitOTFDescCopyDesc  (in, out, err);
} // end OTFDescCopyDesc

void OTFDescIndex (ObitOTFDesc* in) {
  ObitOTFDescIndex (in);
} // end OTFDescIndex

extern ObitInfoList* OTFDescGetList (ObitOTFDesc* in) {
  return ObitInfoListRef(in->info);
}

extern PyObject *OTFDescGetDict(ObitOTFDesc* in) {
  PyObject *outDict = PyDict_New();
  PyObject *list1, *list2, *list3, *list4, *list5, *list6;
  int i;

  PyDict_SetItemString(outDict, "name",    PyString_InternFromString(in->name));
  PyDict_SetItemString(outDict, "object",  PyString_InternFromString(in->object));
  PyDict_SetItemString(outDict, "teles",   PyString_InternFromString(in->teles));
  PyDict_SetItemString(outDict, "origin",  PyString_InternFromString(in->origin));
  PyDict_SetItemString(outDict, "date",    PyString_InternFromString(in->date));
  PyDict_SetItemString(outDict, "obsdat",  PyString_InternFromString(in->obsdat));
  PyDict_SetItemString(outDict, "isort",   PyString_InternFromString(in->isort));
  PyDict_SetItemString(outDict, "bunit",   PyString_InternFromString(in->bunit));
  PyDict_SetItemString(outDict, "JDObs",   PyFloat_FromDouble((double)in->JDObs));
  PyDict_SetItemString(outDict, "obsra",   PyFloat_FromDouble((double)in->obsra));
  PyDict_SetItemString(outDict, "obsdec",  PyFloat_FromDouble((double)in->obsdec));
  PyDict_SetItemString(outDict, "epoch",   PyFloat_FromDouble((double)in->epoch));
  PyDict_SetItemString(outDict, "equinox", PyFloat_FromDouble((double)in->equinox));
  PyDict_SetItemString(outDict, "beamSize",PyFloat_FromDouble((double)in->beamSize));
  PyDict_SetItemString(outDict, "diameter",PyFloat_FromDouble((double)in->diameter));
  PyDict_SetItemString(outDict, "naxis",   PyInt_FromLong((long)in->naxis));
  list1 = PyList_New(OTF_MAXDIM);
  for (i=0; i<OTF_MAXDIM; i++) PyList_SetItem(list1, i, PyString_InternFromString(in->ctype[i]));
  PyDict_SetItemString(outDict, "ctype", list1);
  list2 = PyList_New(OTF_MAXDIM);
  for (i=0; i<OTF_MAXDIM; i++) PyList_SetItem(list2, i, PyFloat_FromDouble((double)in->crval[i]));
  PyDict_SetItemString(outDict, "crval", list2);
  list3 = PyList_New(OTF_MAXDIM);
  for (i=0; i<OTF_MAXDIM; i++) PyList_SetItem(list3, i, PyInt_FromLong((long)in->inaxes[i]));
  PyDict_SetItemString(outDict, "inaxes", list3);
  list4 = PyList_New(OTF_MAXDIM);
  for (i=0; i<OTF_MAXDIM; i++) PyList_SetItem(list4, i, PyFloat_FromDouble((double)in->cdelt[i]));
  PyDict_SetItemString(outDict, "cdelt", list4);
  list5 = PyList_New(OTF_MAXDIM);
  for (i=0; i<OTF_MAXDIM; i++) PyList_SetItem(list5, i, PyFloat_FromDouble((double)in->crpix[i]));
  PyDict_SetItemString(outDict, "crpix", list5);
  list6 = PyList_New(OTF_MAXDIM);
  for (i=0; i<OTF_MAXDIM; i++) PyList_SetItem(list6, i, PyFloat_FromDouble((double)in->crota[i]));
  PyDict_SetItemString(outDict, "crota", list6);
/* Structural members - mostly read only */
  PyDict_SetItemString(outDict, "nrecord",  PyInt_FromLong((long)in->nrecord));
  PyDict_SetItemString(outDict, "firstRec", PyInt_FromLong((long)in->firstRec));
  PyDict_SetItemString(outDict, "numRecBuff", PyInt_FromLong((long)in->numRecBuff));
  PyDict_SetItemString(outDict, "ncol",     PyInt_FromLong((long)in->ncol));
  PyDict_SetItemString(outDict, "numDesc",  PyInt_FromLong((long)in->numDesc));
  PyDict_SetItemString(outDict, "lrec",     PyInt_FromLong((long)in->lrec));
  PyDict_SetItemString(outDict, "iloct",    PyInt_FromLong((long)in->iloct));
  PyDict_SetItemString(outDict, "ilocti",   PyInt_FromLong((long)in->ilocti));
  PyDict_SetItemString(outDict, "iloctar",  PyInt_FromLong((long)in->iloctar));
  PyDict_SetItemString(outDict, "ilocscan", PyInt_FromLong((long)in->ilocscan));
  PyDict_SetItemString(outDict, "ilocra",   PyInt_FromLong((long)in->ilocra));
  PyDict_SetItemString(outDict, "ilocdec",  PyInt_FromLong((long)in->ilocdec));
  PyDict_SetItemString(outDict, "ilocrot",  PyInt_FromLong((long)in->ilocrot));
  PyDict_SetItemString(outDict, "iloccal",  PyInt_FromLong((long)in->iloccal));
  PyDict_SetItemString(outDict, "ilocdata", PyInt_FromLong((long)in->ilocdata));
  PyDict_SetItemString(outDict, "jlocfeed", PyInt_FromLong((long)in->jlocfeed));
  PyDict_SetItemString(outDict, "jlocs",    PyInt_FromLong((long)in->jlocs));
  PyDict_SetItemString(outDict, "jlocf",    PyInt_FromLong((long)in->jlocf));
  PyDict_SetItemString(outDict, "incfeed" , PyInt_FromLong((long)in->incfeed));
  PyDict_SetItemString(outDict, "incs",     PyInt_FromLong((long)in->incs));
  PyDict_SetItemString(outDict, "incf",     PyInt_FromLong((long)in->incf));
  return outDict;
}

extern void OTFDescSetDict(ObitOTFDesc* in, PyObject *inDict) {
  PyObject *list1, *list2, *list4, *list5, *list6;
  char *tstr;
  int i;

  if (!PyDict_Check(inDict)) PyErr_SetString(PyExc_TypeError,"Input not a Dict");

  tstr = PyString_AsString(PyDict_GetItemString(inDict, "object"));
  strncpy (in->object, tstr, OTFLEN_VALUE); in->object[OTFLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obsdat"));
  strncpy (in->obsdat, tstr, OTFLEN_VALUE); in->obsdat[OTFLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "bunit"));
  strncpy (in->bunit, tstr, OTFLEN_VALUE); in->bunit[OTFLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "teles"));
  strncpy (in->teles, tstr, OTFLEN_VALUE); in->teles[OTFLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "origin"));
  strncpy (in->origin, tstr, OTFLEN_VALUE); in->origin[OTFLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "date"));
  strncpy (in->date, tstr, OTFLEN_VALUE); in->date[OTFLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "isort"));
  strncpy (in->isort, tstr, 2); in->isort[2]=0;
  in->JDObs   = PyFloat_AsDouble(PyDict_GetItemString(inDict, "JDObs"));
  in->obsra   = PyFloat_AsDouble(PyDict_GetItemString(inDict, "obsra"));
  in->obsdec  = PyFloat_AsDouble(PyDict_GetItemString(inDict, "obsdec"));
  in->epoch   = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "epoch"));
  in->equinox = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "equinox"));
  in->beamSize= (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "beamSize"));
  in->diameter= (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "diameter"));
  list1 = PyDict_GetItemString(inDict, "ctype");
  for (i=0; i<OTF_MAXDIM; i++) {
    tstr = PyString_AsString(PyList_GetItem(list1, i)); 
    strncpy (in->ctype[i], tstr, OTFLEN_KEYWORD); in->ctype[i][OTFLEN_KEYWORD]=0;
  }
  list2 = PyDict_GetItemString(inDict, "crval");
  for (i=0; i<OTF_MAXDIM; i++) in->crval[i] = PyFloat_AsDouble(PyList_GetItem(list2, i));
  list4 = PyDict_GetItemString(inDict, "cdelt");
  for (i=0; i<OTF_MAXDIM; i++) in->cdelt[i] = (float)PyFloat_AsDouble(PyList_GetItem(list4, i));
  list5 = PyDict_GetItemString(inDict, "crpix");
  for (i=0; i<OTF_MAXDIM; i++) in->crpix[i] = (float)PyFloat_AsDouble(PyList_GetItem(list5, i));
  list6 = PyDict_GetItemString(inDict, "crota");
  for (i=0; i<OTF_MAXDIM; i++) in->crota[i] = (float)PyFloat_AsDouble(PyList_GetItem(list6, i));
  /* Structural members - mostly read only */
  in->nrecord  = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "nrecord"));
  in->firstRec = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "firstRec"));
  in->numRecBuff = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "numRecBuff"));
}

ObitOTFDesc* OTFDescRef (ObitOTFDesc* in) {
  return ObitOTFDescRef (in);
} // end OTFDescRef

ObitOTFDesc* OTFDescUnref (ObitOTFDesc* in) {
  if (!ObitOTFDescIsA(in)) return NULL;
  return ObitOTFDescUnref (in);
} // end OTFDescUnref

extern int OTFDescIsA (ObitOTFDesc* in) {
  return ObitOTFDescIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitOTFDesc *me;
} OTFDesc;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitOTFDesc *me;
} OTFDesc;

%addmethods OTFDesc { 
  OTFDesc(char *name) {
     OTFDesc *out;
     out = (OTFDesc *) malloc(sizeof(OTFDesc));
     if (strcmp(name, "None")) out->me = OTFDescCreate(name);
     else  out->me = NULL;
     return out;
   }
  ~OTFDesc() {
    self->me = OTFDescUnref(self->me);
    free(self);
  }
};

