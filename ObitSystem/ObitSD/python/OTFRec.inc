/* $Id: OTFRec.inc,v 1.1 2008/01/14 12:15:40 bcotton Exp $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for OTFRec type                            */
/*                                                                    */
/*;  Copyright (C)2008                                                */
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
#include "ObitOTF.h"
#include "ObitOTFDesc.h"
#include "ObitOTFArrayGeom.h"
%}

%inline %{
extern PyObject* OTFRecGet (ObitOTF* inOTF, ObitErr *err) {
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean doCalSelect;
  ObitIOCode iretCode;
  ofloat *recp=NULL;
  olong ndata, i;
  PyObject *rec;
  PyObject *vList, *tup;

   rec = PyDict_New();

   doCalSelect = FALSE;
   ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, dim, &doCalSelect);

    if (doCalSelect) iretCode = ObitOTFReadSelect (inOTF, inOTF->buffer, err);
    else iretCode = ObitOTFRead (inOTF, inOTF->buffer, err);

    if (iretCode==OBIT_IO_EOF)  {
      PyDict_SetItemString(rec, "EOF",PyInt_FromLong((long)1));
      return rec;
    }
    if (iretCode!=OBIT_IO_OK) return rec;

    PyDict_SetItemString(rec, "recNo", PyInt_FromLong((long)inOTF->myDesc->firstRec));
    PyDict_SetItemString(rec, "time",  PyFloat_FromDouble((double)inOTF->buffer[inOTF->myDesc->iloct]));
    PyDict_SetItemString(rec, "timei", PyFloat_FromDouble((double)inOTF->buffer[inOTF->myDesc->ilocti]));
    PyDict_SetItemString(rec, "scan",  PyInt_FromLong((long)(inOTF->buffer[inOTF->myDesc->ilocscan]+0.5)));
    PyDict_SetItemString(rec, "target",PyInt_FromLong((long)(inOTF->buffer[inOTF->myDesc->iloctar]+0.5)));
    PyDict_SetItemString(rec, "ra",    PyFloat_FromDouble((double)inOTF->buffer[inOTF->myDesc->ilocra]));
    PyDict_SetItemString(rec, "dec",   PyFloat_FromDouble((double)inOTF->buffer[inOTF->myDesc->ilocdec]));
    PyDict_SetItemString(rec, "rot",   PyFloat_FromDouble((double)inOTF->buffer[inOTF->myDesc->ilocrot]));
    PyDict_SetItemString(rec, "cal",   PyFloat_FromDouble((double)inOTF->buffer[inOTF->myDesc->iloccal]));
    //  How much data?
    ndata = 1;
    for (i=0; i<inOTF->myDesc->naxis; i++) ndata *= MAX (1, inOTF->myDesc->inaxes[i]);
    ndata /= inOTF->myDesc->incdatawt;
    vList = PyList_New(ndata);
    recp = &inOTF->buffer[inOTF->myDesc->numDesc];
    for (i=0; i<ndata; i++) {
      tup = PyTuple_New(2);
      // Does data have weight?
      if (inOTF->myDesc->incdatawt>1) { // yes
        PyTuple_SetItem (tup, 0, PyFloat_FromDouble((double)*(recp+0)));
        PyTuple_SetItem (tup, 1, PyFloat_FromDouble((double)*(recp+1)));
        recp += 2;
      } else { // no
        PyTuple_SetItem (tup, 0, PyFloat_FromDouble((double)*(recp)));
        PyTuple_SetItem (tup, 1, PyFloat_FromDouble((double)(1.0)));
        recp++;
      }
      PyList_SetItem(vList, i, tup);
    }

    PyDict_SetItemString(rec, "data", vList);
    Py_DECREF(vList);

   return rec;
} // end OTFRecGet

extern void OTFRecSet (PyObject* rec, ObitOTF* outOTF, ObitErr *err) {
  ObitIOCode oretCode;
  ofloat *recp=NULL;
  olong ndata, i, len;
  PyObject *vList, *tup;

  outOTF->buffer[outOTF->myDesc->iloct]    = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(rec, "time"));
  outOTF->buffer[outOTF->myDesc->ilocti]   = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(rec, "timei"));
  outOTF->buffer[outOTF->myDesc->ilocscan] = (ofloat)PyInt_AsLong(PyDict_GetItemString(rec, "scan"));
  outOTF->buffer[outOTF->myDesc->iloctar]  = (ofloat)PyInt_AsLong(PyDict_GetItemString(rec, "target"));
  outOTF->buffer[outOTF->myDesc->ilocra]   = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(rec, "ra"));
  outOTF->buffer[outOTF->myDesc->ilocdec]  = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(rec, "dec"));
  outOTF->buffer[outOTF->myDesc->ilocrot]  = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(rec, "rot"));
  outOTF->buffer[outOTF->myDesc->iloccal]  = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(rec, "cal"));

   vList = PyDict_GetItemString(rec, "data");
   len = PyList_Size(vList);
   //  How much data?
   ndata = 1;
   for (i=0; i<outOTF->myDesc->naxis; i++) ndata *= MAX (1, outOTF->myDesc->inaxes[i]);
   ndata /= outOTF->myDesc->incdatawt;
   if (len!=ndata) {
     PyErr_SetString(PyExc_TypeError,"OTFRec incompatible with ObitOTF");
     return;
   }
    recp = &outOTF->buffer[outOTF->myDesc->numDesc];
    for (i=0; i<len; i++) {
      tup = PyList_GetItem(vList, i);
      *(recp++) = (ofloat)PyFloat_AsDouble(PyTuple_GetItem (tup, 0));
      // Does data have weight?
      if (outOTF->myDesc->incdatawt>1) { // yes
          *(recp++) = (ofloat)PyFloat_AsDouble(PyTuple_GetItem (tup, 1));
      }
    }

    outOTF->myDesc->numRecBuff = 1;
    oretCode = ObitOTFWrite (outOTF, outOTF->buffer, err);

} // end OTFRecSet

// return dict with "xpos" = array of RAs (deg) of detectors
//  "ypos" = array of decs (deg) of detectors
extern PyObject* OTFRecGeomProj (PyObject* rec, ObitOTF* inOTF) {
  ofloat ra, dec, rot, *xpos=NULL, *ypos=NULL;
  olong ndetect, i;
  PyObject *out;
  PyObject *pList, *d;

  out = PyDict_New();

  ndetect = inOTF->geom->numberDetect;
  xpos = g_malloc0(ndetect*sizeof(ofloat));
  ypos = g_malloc0(ndetect*sizeof(ofloat));
  ra  = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(rec, "ra"));
  dec = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(rec, "dec"));
  rot = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(rec, "rot"));

  ObitOTFArrayGeomCoord (inOTF->geom, ra, dec, rot, xpos, ypos);

  pList = PyList_New(ndetect);
  for (i=0; i<ndetect; i++) {
    d = PyFloat_FromDouble((double)xpos[i]);
    PyList_SetItem(pList, i, d);
  }
  PyDict_SetItemString(out, "xpos", pList);
  Py_DECREF(pList);

  pList = PyList_New(ndetect);
  for (i=0; i<ndetect; i++) {
    d = PyFloat_FromDouble((double)ypos[i]);
    PyList_SetItem(pList, i, d);
  }
  PyDict_SetItemString(out, "ypos", pList);
  Py_DECREF(pList);

  if (xpos) g_free(xpos);
  if (ypos) g_free(ypos);
  return out;
} // end  OTFRecGeomProj
%}

