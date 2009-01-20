/* $Id:    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageDesc type                         */
/*                                                                    */
/*;  Copyright (C) 2005,2008                                          */
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
#include "ObitTableDesc.h"
%}

%inline %{
// Routine to remove trailing blanks from a string 
static void TableDescDeeBlank (gchar *in) 
{
  olong i;

  for (i=strlen(in)-1; i>=0; i--) {
     if (in[i]==' ') in[i] = 0;
     else if (in[i]!=' ') break;
  }

} // end TableDescDeeBlank

extern ObitTableDesc* TableDescCreate (char *name) {
  return newObitTableDesc (name);
} // end TableDescCreate

extern ObitTableDesc* TableDescCopy (ObitTableDesc* in, 
		              ObitTableDesc* out, ObitErr *err) {
  return ObitTableDescCopy (in, out, err);
} // end TableDescCopy

extern void TableDescCopyDesc (ObitTableDesc* in, ObitTableDesc* out,
			ObitErr *err) {
  ObitTableDescCopyDesc  (in, out, err);
} // end TableDescCopyDesc

extern void TableDescIndex (ObitTableDesc* in) {
  ObitTableDescIndex (in);
} // end TableDescIndex

extern ObitInfoList* TableDescGetList (ObitTableDesc* in) {
  return ObitInfoListRef(in->info);
}
 
extern PyObject *TableDescGetDict(ObitTableDesc* in) {
  PyObject *outDict = PyDict_New();
  PyObject *list, *value;
  gchar *ctemp;
  int i, pos = 0;

  /* test validity */
  if (in->nfield<=0) {
    PyErr_SetString(PyExc_TypeError,"Input not fully defined");
    return outDict;
  }

  PyDict_SetItemString(outDict, "Table name", PyString_InternFromString(in->TableName));
  PyDict_SetItemString(outDict, "version",PyInt_FromLong((long)in->version));
  PyDict_SetItemString(outDict, "nrow",   PyInt_FromLong((long)in->nrow));
  PyDict_SetItemString(outDict, "lrow",   PyInt_FromLong((long)in->lrow));
  PyDict_SetItemString(outDict, "sortOrder1", PyInt_FromLong((long)in->sort[0]));
  PyDict_SetItemString(outDict, "sortOrder2", PyInt_FromLong((long)in->sort[1]));

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) PyList_SetItem(list, i, PyInt_FromLong((long)in->repeat[i]));
  PyDict_SetItemString(outDict, "repeat", list);

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) PyList_SetItem(list, i, PyInt_FromLong((long)in->dim[i][0]));
  PyDict_SetItemString(outDict, "dim0", list);

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) PyList_SetItem(list, i, PyInt_FromLong((long)in->dim[i][1]));
  PyDict_SetItemString(outDict, "dim1", list);

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) PyList_SetItem(list, i, PyInt_FromLong((long)in->dim[i][2]));
  PyDict_SetItemString(outDict, "dim2", list);

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) PyList_SetItem(list, i, PyInt_FromLong((long)in->type[i]));
  PyDict_SetItemString(outDict, "type", list);

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) {
    ctemp =  g_strdup(in->FieldName[i]);
    TableDescDeeBlank(ctemp);
    PyList_SetItem(list, i, PyString_InternFromString(ctemp));
    g_free (ctemp);
  }
  PyDict_SetItemString(outDict, "FieldName", list);

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) PyList_SetItem(list, i, PyString_InternFromString(in->FieldUnit[i]));
  for (i=0; i<in->nfield; i++) {
    ctemp =  g_strdup(in->FieldUnit[i]);
    TableDescDeeBlank(ctemp);
    PyList_SetItem(list, i, PyString_InternFromString(ctemp));
    g_free (ctemp);
  }
  PyDict_SetItemString(outDict, "FieldUnit", list);

  /* Discard references to newly created objects. */
  while (PyDict_Next(outDict, &pos, NULL, &value))
    Py_DECREF(value);

  return outDict;
} // end TableDescGetDict

extern void TableDescSetDict(ObitTableDesc* in, PyObject *inDict) {
  PyObject *list;
  char *tstr;
  int i, number;

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return;
  }

  in->version  = PyInt_AsLong(PyDict_GetItemString(inDict, "version"));
  in->sort[0]  = PyInt_AsLong(PyDict_GetItemString(inDict, "sortOrder1"));
  in->sort[1]  = PyInt_AsLong(PyDict_GetItemString(inDict, "sortOrder2"));

  list = PyDict_GetItemString(inDict, "FieldName");
  number = MIN (in->nfield, PyList_Size(list));
  for (i=0; i<number; i++) {
    tstr = PyString_AsString(PyList_GetItem(list, i));
    if (in->FieldName[i]) g_free(in->FieldName[i]);
    in->FieldName[i] = g_strdup(tstr);
  }

  list = PyDict_GetItemString(inDict, "FieldUnit");
  number = MIN (in->nfield, PyList_Size(list));
  for (i=0; i<number; i++) {
    tstr = PyString_AsString(PyList_GetItem(list, i));
    if (in->FieldUnit[i]) g_free(in->FieldUnit[i]);
    in->FieldUnit[i] = g_strdup(tstr);
  }

} // end TableDescSetDict

//  Define a descriptor from a python dict
extern ObitTableDesc* TableDescDef(PyObject *inDict) {
  ObitTableDesc *out=NULL;
  PyObject *fieldname, *fieldunit, *dim0, *dim1, *dim2, *repeat, *type;
  PyObject *thing;
  gchar *tstr;
  olong i, nfield;

 
  // Define output
  out = newObitTableDesc(NULL);

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return out;
  }

  // Check input
  fieldname = PyDict_GetItemString(inDict, "FieldName");
  if (!fieldname) {
    PyErr_SetString(PyExc_TypeError,"FieldName Array not found");
    return out;
  }

  nfield = PyList_Size(fieldname);

  fieldunit = PyDict_GetItemString(inDict, "FieldUnit");
  if (!fieldunit) {
    PyErr_SetString(PyExc_TypeError,"FieldUnit Array not found");
    return out;
  }
  if (PyList_Size(fieldunit)!=nfield) {
    PyErr_SetString(PyExc_TypeError,"FieldUnit Array wrong dimension");
    return out;
  }

  repeat = PyDict_GetItemString(inDict, "repeat");
  if (!repeat) {
    PyErr_SetString(PyExc_TypeError,"repeat Array not found");
    return;
  }
  if (PyList_Size(repeat)!=nfield) {
    PyErr_SetString(PyExc_TypeError,"repeat Array wrong dimension");
    return out;
  }

  dim0   = PyDict_GetItemString(inDict, "dim0");
  if (!dim0) {
    PyErr_SetString(PyExc_TypeError,"dim0 Array not found");
    return out;
  }
  if (PyList_Size(dim0)!=nfield) {
    PyErr_SetString(PyExc_TypeError,"dim0 Array wrong dimension");
    return out;
  }

  dim1   = PyDict_GetItemString(inDict, "dim1");
  if (!dim1) {
    PyErr_SetString(PyExc_TypeError,"dim1 Array not found");
    return out;
  }
  if (PyList_Size(dim1)!=nfield) {
    PyErr_SetString(PyExc_TypeError,"dim1 Array wrong dimension");
    return out;
  }

  dim2   = PyDict_GetItemString(inDict, "dim2");
  if (!dim2) {
    PyErr_SetString(PyExc_TypeError,"dim2 Array not found");
    return out;
  }
  if (PyList_Size(dim2)!=nfield) {
    PyErr_SetString(PyExc_TypeError,"dim2 Array wrong dimension");
    return out;
  }

  type   = PyDict_GetItemString(inDict, "type");
  if (!type) {
    PyErr_SetString(PyExc_TypeError,"type Array not found");
    return out;
  }
  if (PyList_Size(type)!=nfield) {
    PyErr_SetString(PyExc_TypeError,"type Array wrong dimension");
    return out;
  }

  // Resize output
  ObitTableDescRealloc (out, nfield);

  thing = PyDict_GetItemString(inDict, "version");
  if (!thing) {
    PyErr_SetString(PyExc_TypeError,"version not found");
    return out;
  }
  out->version  = PyInt_AsLong(thing);
  thing = PyDict_GetItemString(inDict, "sortOrder1");
  if (!thing) {
    PyErr_SetString(PyExc_TypeError,"sortOrder1 not found");
    return out;
  }
  out->sort[0]  = PyInt_AsLong(thing);
  thing = PyDict_GetItemString(inDict, "sortOrder2");
  if (!thing) {
    PyErr_SetString(PyExc_TypeError,"sortOrder2 not found");
    return out;
  }
  out->sort[1]  = PyInt_AsLong(thing);
  thing = PyDict_GetItemString(inDict, "Table name");
  if (!thing) {
    PyErr_SetString(PyExc_TypeError,"Table name not found");
    return out;
  }
  tstr = PyString_AsString(thing);
  if (out->TableName) g_free(out->TableName);
  out->TableName = g_strdup(tstr);

  // field names
  for (i=0; i<nfield; i++) {
    tstr = PyString_AsString(PyList_GetItem(fieldname, i));
    if (out->FieldName[i]) g_free(out->FieldName[i]);
    out->FieldName[i] = g_strdup(tstr);
  }

  // field units
  for (i=0; i<nfield; i++) {
    tstr = PyString_AsString(PyList_GetItem(fieldunit, i));
    if (out->FieldUnit[i]) g_free(out->FieldUnit[i]);
    out->FieldUnit[i] = g_strdup(tstr);
  }

  // other field stuff
  for (i=0; i<nfield; i++) {
    out->repeat[i] = (olong)PyInt_AsLong(PyList_GetItem(repeat, i));
    out->type[i]   = (olong)PyInt_AsLong(PyList_GetItem(type, i));
    out->dim[i][0] = (olong)PyInt_AsLong(PyList_GetItem(dim0, i));
    out->dim[i][1] = (olong)PyInt_AsLong(PyList_GetItem(dim1, i));
    out->dim[i][2] = (olong)PyInt_AsLong(PyList_GetItem(dim2, i));
  }

  ObitTableDescIndex (out);

  return out;

} // end TableDescDef

ObitTableDesc* TableDescRef (ObitTableDesc* in) {
  return ObitTableDescRef (in);
} // end TableDescRef

ObitTableDesc* TableDescUnref (ObitTableDesc* in) {
  if (!ObitTableDescIsA(in)) return NULL;
  return ObitTableDescUnref (in);
} // end TableDescUnref

extern int TableDescIsA (ObitTableDesc* in) {
  return ObitTableDescIsA(in);
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitTableDesc *me;
} TableDesc;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitTableDesc *me;
} TableDesc;

%addmethods TableDesc { 
  TableDesc(char *name) {
     TableDesc *out;
     out = (TableDesc *) malloc(sizeof(TableDesc));
     if (strcmp(name, "None")) out->me = TableDescCreate (name);
     else out->me = NULL;
     return out;
   }
  ~TableDesc() {
    self->me = TableDescUnref(self->me);
    free(self);
  }
};

