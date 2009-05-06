/* $Id:  $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for OTFArrayGeom type                      */
/*                                                                    */
/*;  Copyright (C) 2009                                               */
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
#include "ObitOTFArrayGeom.h"
#include "ObitOTFSkyModel.h"
%}

%inline %{
ObitOTFArrayGeom* OTFArrayGeomCreate (gchar *name) {
  return newObitOTFArrayGeom (name);
} // end OTFArrayGeomCreate

ObitOTFArrayGeom* OTFArrayGeomCopy (ObitOTFArrayGeom* in, 
		              ObitOTFArrayGeom* out, ObitErr *err) {
  return ObitOTFArrayGeomCopy (in, out, err);
} // end OTFArrayGeomCopy

extern ObitInfoList* OTFArrayGeomGetList (ObitOTFArrayGeom* in) {
  return ObitInfoListRef(in->info);
}

extern PyObject *OTFArrayGeomGetDict (ObitOTFArrayGeom* in) {
  PyObject *outDict = PyDict_New();
  PyObject *list1, *list2;
  int i;

  PyDict_SetItemString(outDict, "numberDetect", PyInt_FromLong((long)in->numberDetect));
  PyDict_SetItemString(outDict, "RefDate",      PyString_InternFromString(in->RefDate));
  PyDict_SetItemString(outDict, "TimeSys",      PyString_InternFromString(in->TimeSys));
  PyDict_SetItemString(outDict, "TeleX",        PyFloat_FromDouble((double)in->TeleX));
  PyDict_SetItemString(outDict, "TeleY",        PyFloat_FromDouble((double)in->TeleY));
  PyDict_SetItemString(outDict, "TeleZ",        PyFloat_FromDouble((double)in->TeleZ));
  PyDict_SetItemString(outDict, "DegDay",       PyFloat_FromDouble((double)in->DegDay));
  PyDict_SetItemString(outDict, "GSTiat0",      PyFloat_FromDouble((double)in->GSTiat0));
  PyDict_SetItemString(outDict, "PolarX",       PyFloat_FromDouble((double)in->PolarX));
  PyDict_SetItemString(outDict, "PolarY",       PyFloat_FromDouble((double)in->PolarY));
  PyDict_SetItemString(outDict, "ut1Utc",       PyFloat_FromDouble((double)in->ut1Utc));
  PyDict_SetItemString(outDict, "dataUtc",      PyFloat_FromDouble((double)in->dataUtc));
  PyDict_SetItemString(outDict, "iatUtc",       PyFloat_FromDouble((double)in->iatUtc));
  PyDict_SetItemString(outDict, "lat",          PyFloat_FromDouble((double)in->lat));
  PyDict_SetItemString(outDict, "lon",          PyFloat_FromDouble((double)in->lon));
  PyDict_SetItemString(outDict, "LSTiat0",      PyFloat_FromDouble((double)in->LSTiat0));
  PyDict_SetItemString(outDict, "RadDay",       PyFloat_FromDouble((double)in->RadDay));
  PyDict_SetItemString(outDict, "dataIat",      PyFloat_FromDouble((double)in->dataIat));

  // Offsets
  list1 = PyList_New((long)in->numberDetect);
  for (i=0; i<in->numberDetect; i++) PyList_SetItem(list1, i, PyFloat_FromDouble((double)in->azOffset[i]));
  PyDict_SetItemString(outDict, "azOffset", list1);
  list2 = PyList_New((long)in->numberDetect);
  for (i=0; i<in->numberDetect; i++) PyList_SetItem(list2, i, PyFloat_FromDouble((double)in->elOffset[i]));
  PyDict_SetItemString(outDict, "elOffset", list2);
  return outDict;
} // end OTFArrayGeomGetDict

extern void OTFArrayGeomSetDict (ObitOTFArrayGeom* in, PyObject *inDict) {
  PyObject *list1, *list2;
  char *tstr;
  int i;

  if (!PyDict_Check(inDict)) PyErr_SetString(PyExc_TypeError,"Input not a Dict");

  // Strings
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (in->RefDate, tstr, 12); in->RefDate[11]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "TimeSys"));
  strncpy (in->TimeSys, tstr, 4); in->TimeSys[3]=0;

   // offsets
  list1 = PyDict_GetItemString(inDict, "azOffset");
  for (i=0; i<in->numberDetect; i++) {
    in->azOffset[i] = (ofloat)PyFloat_AsDouble(PyList_GetItem(list1, i));
  }
  Py_DECREF(list1);
  list2 = PyDict_GetItemString(inDict, "elOffset");
  for (i=0; i<in->numberDetect; i++) {
    in->elOffset[i] = (ofloat)PyFloat_AsDouble(PyList_GetItem(list2, i));
  }
  Py_DECREF(list2);

  // other stuff
  in->TeleX =   (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "TeleX"));
  in->TeleY =   (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "TeleY"));
  in->TeleZ =   (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "TeleZ"));
  in->DegDay =  (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "DegDay"));
  in->GSTiat0 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "GSTiat0"));
  in->PolarX =  (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "PolarX"));
  in->PolarY =  (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "PolarY"));
  in->ut1Utc =  (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ut1Utc"));
  in->dataUtc = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "dataUtc"));
  in->iatUtc =  (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "iatUtc"));
  in->lat =     (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "lat"));
  in->lon =     (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "lon"));
  in->LSTiat0 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "LSTiat0"));
  in->RadDay =  (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "RadDay"));
  in->dataIat = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "dataIat"));
} // end OTFArrayGeomSetDict

/** Public: Get Parallactic angle for a given time and direction */
extern float OTFArrayGeomParAng (ObitOTFArrayGeom *in, float time, float ra, float dec) {
  return (float)ObitOTFArrayGeomParAng (in, (ofloat)time, (ofloat)ra, (ofloat)dec);
} // end OTFArrayGeomParAng 

/** Public: Get Elevation for a given time and direction */
extern float OTFArrayGeomElev (ObitOTFArrayGeom *in, float time, float ra, float dec) {
  return (float)ObitOTFArrayGeomElev (in, (ofloat)time, (ofloat)ra, (ofloat)dec);
} // end OTFArrayGeomElev 

/** Public: Get detector coordinates on sky, return (x,y) as list of lists */
extern PyObject *OTFArrayGeomCoord (ObitOTFArrayGeom *in, float raPoint, float decPoint, float rot) {
  ofloat *x=NULL, *y=NULL;
  PyObject *list1, *list2, *outList;
  int i;

  x = g_malloc0(in->numberDetect*sizeof(ofloat));
  y = g_malloc0(in->numberDetect*sizeof(ofloat));
  ObitOTFArrayGeomCoord (in, (ofloat)raPoint, (ofloat)decPoint, (ofloat)rot, x, y);
  list1   = PyList_New(in->numberDetect);
  list2   = PyList_New(in->numberDetect);
  for (i=0; i<in->numberDetect; i++) {
    PyList_SetItem(list1, i, PyFloat_FromDouble((double)x[i]));
    PyList_SetItem(list2, i, PyFloat_FromDouble((double)y[i]));
  }
  outList = PyList_New(2);
  PyList_SetItem(outList, 0, list1);
  PyList_SetItem(outList, 1, list2);
  if (x)  g_free(x);
  if (y)  g_free(y);
  return outList;
} // end  OTFArrayGeomCoord

/** Public: Get detector locations projected onto a plane, return (x,y) as list of lists */
// ProjCode 0=-SIN , 1=-ARC, 2=-TAN
extern PyObject *OTFArrayGeomProj (ObitOTFArrayGeom *in, float raPoint, float decPoint, float rot,
		      float raProj, float decProj, int ProjCode) {
  ofloat *x=NULL, *y=NULL;
  PyObject *list1, *list2, *outList;
  ObitOTFProj Proj = OBIT_OTF_SIN;
  int i;
  
  if (ProjCode==0)      Proj = OBIT_OTF_SIN;
  else if (ProjCode==1) Proj = OBIT_OTF_ARC;
  else if (ProjCode==2) Proj = OBIT_OTF_TAN;
  else                  Proj = OBIT_OTF_SIN;

  x = g_malloc0(in->numberDetect*sizeof(ofloat));
  y = g_malloc0(in->numberDetect*sizeof(ofloat));
  ObitOTFArrayGeomProj (in, (ofloat)raPoint, (ofloat)decPoint, (ofloat)rot, 
      (ofloat)raProj, (ofloat)decProj,  Proj, x, y);
  list1   = PyList_New(in->numberDetect);
  list2   = PyList_New(in->numberDetect);
  for (i=0; i<in->numberDetect; i++) {
    PyList_SetItem(list1, i, PyFloat_FromDouble((double)x[i]));
    PyList_SetItem(list2, i, PyFloat_FromDouble((double)y[i]));
  }
  outList = PyList_New(2);
  PyList_SetItem(outList, 0, list1);
  PyList_SetItem(outList, 1, list2);
  if (x)  g_free(x);
  if (y)  g_free(y);
  return outList;
} // end  OTFArrayGeomProj

/** Public: Offset a celestial position in az, el, return (raPoint, decPoint) as list */
extern PyObject *OTFArrayGeomCorrPoint (float azOff, float elOff, float pa) {
  ofloat x, y;
  PyObject *list1;

  ObitOTFArrayGeomCorrPoint ((ofloat)azOff, (ofloat)elOff, (ofloat)pa, &x, &y);
  list1 = PyList_New(2);
  PyList_SetItem(list1, 0, PyFloat_FromDouble((double)x));
  PyList_SetItem(list1, 1, PyFloat_FromDouble((double)y));
  return list1;
} // end  OTFArrayGeomCorrPoint

ObitOTFArrayGeom* OTFArrayGeomRef (ObitOTFArrayGeom* in) {
  return ObitOTFArrayGeomRef (in);
} // end OTFArrayGeomRef

ObitOTFArrayGeom* OTFArrayGeomUnref (ObitOTFArrayGeom* in) {
  if (!ObitOTFArrayGeomIsA(in)) return NULL;
  return ObitOTFArrayGeomUnref (in);
} // end OTFArrayGeomUnref

extern int OTFArrayGeomIsA (ObitOTFArrayGeom* in) {
  return ObitOTFArrayGeomIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitOTFArrayGeom *me;
} OTFArrayGeom;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitOTFArrayGeom *me;
} OTFArrayGeom;

%addmethods OTFArrayGeom { 
  OTFArrayGeom(char *name) {
     OTFArrayGeom *out;
     out = (OTFArrayGeom *) malloc(sizeof(OTFArrayGeom));
     if (strcmp(name, "None")) out->me = OTFArrayGeomCreate(name);
     else  out->me = NULL;
     return out;
   }
  ~OTFArrayGeom() {
    self->me = OTFArrayGeomUnref(self->me);
    free(self);
  }
};

