/* $Id$   */  
/*--------------------------------------------------------------------*/
/* Swig module description for UVVis type                             */
/*                                                                    */
/*;  Copyright (C)2007                                                */
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
#include "ObitUV.h"
#include "ObitUVDesc.h"
%}

%inline %{
extern PyObject* UVVisGet (ObitUV* inUV, ObitErr *err) {
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean doCalSelect;
  ObitIOCode iretCode;
  ofloat cbase, *visp=NULL;
  olong ant1, ant2, suid, fqid, icorr;
  PyObject *vis;
  PyObject *vList, *cvis, *cx, *tup;

   vis = PyDict_New();

   doCalSelect = FALSE;
   ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);

    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);

    if (iretCode==OBIT_IO_EOF)  {
      PyDict_SetItemString(vis, "EOF",PyInt_FromLong((long)1));
      return vis;
    }
    if (iretCode!=OBIT_IO_OK) return vis;

    PyDict_SetItemString(vis, "visNo",  PyInt_FromLong((long)inUV->myDesc->firstVis));
    cbase = inUV->buffer[inUV->myDesc->ilocb]; /* Baseline */
    ant1 = (cbase / 256.0) + 0.001;
    ant2 = (cbase - ant1 * 256) + 0.001;
    PyDict_SetItemString(vis, "ant1",    PyInt_FromLong((long)ant1));
    PyDict_SetItemString(vis, "ant2",    PyInt_FromLong((long)ant2));
    PyDict_SetItemString(vis, "u", PyFloat_FromDouble((double)inUV->buffer[inUV->myDesc->ilocu]));
    PyDict_SetItemString(vis, "v", PyFloat_FromDouble((double)inUV->buffer[inUV->myDesc->ilocv]));
    PyDict_SetItemString(vis, "w", PyFloat_FromDouble((double)inUV->buffer[inUV->myDesc->ilocw]));
    PyDict_SetItemString(vis, "time", PyFloat_FromDouble((double)inUV->buffer[inUV->myDesc->iloct]));
    if (inUV->myDesc->ilocsu>=0) {
        suid = (olong)(inUV->buffer[inUV->myDesc->ilocsu]);
        PyDict_SetItemString(vis, "suid", PyInt_FromLong((long)suid));
    }
    if (inUV->myDesc->ilocfq>=0) {
        fqid = (olong)(inUV->buffer[inUV->myDesc->ilocfq]);
        PyDict_SetItemString(vis, "fqid", PyInt_FromLong((long)fqid));
    }
    vList = PyList_New(inUV->myDesc->ncorr);
    visp = &inUV->buffer[inUV->myDesc->nrparm];
    for (icorr=0; icorr<inUV->myDesc->ncorr; icorr++) {
      tup = PyTuple_New(2);
      cx = PyComplex_FromDoubles((double)*visp, (double)*(visp+1));
      PyTuple_SetItem (tup, 0, cx);
      PyTuple_SetItem (tup, 1, PyFloat_FromDouble((double)*(visp+2)));
      PyList_SetItem(vList, icorr, tup);
      /*Py_DECREF(tup);*/
      visp += 3;
    }

    PyDict_SetItemString(vis, "vis", vList);
    Py_DECREF(vList);

   return vis;
} // end UVVisGet

extern void UVVisSet (PyObject* vis, ObitUV* outUV, ObitErr *err) {
  ObitIOCode oretCode;
  ofloat *visp=NULL;
  olong ant1, ant2, icorr, len;
  PyObject *vList, *cvis, *cx, *tup;


    ant1 = (olong)PyInt_AsLong(PyDict_GetItemString(vis, "ant1"));
    ant2 = (olong)PyInt_AsLong(PyDict_GetItemString(vis, "ant2"));
    outUV->buffer[outUV->myDesc->ilocb] = (ofloat)(ant1*256 + ant2);
    outUV->buffer[outUV->myDesc->ilocu] = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(vis, "u"));
    outUV->buffer[outUV->myDesc->ilocv] = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(vis, "v"));
    outUV->buffer[outUV->myDesc->ilocw] = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(vis, "w"));
    outUV->buffer[outUV->myDesc->iloct] = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(vis, "time"));
    if (outUV->myDesc->ilocsu>=0) {
        outUV->buffer[outUV->myDesc->ilocsu] = (ofloat)PyInt_AsLong(PyDict_GetItemString(vis, "suid"));
    }
    if (outUV->myDesc->ilocfq>=0) {
        outUV->buffer[outUV->myDesc->ilocfq] = (ofloat)PyInt_AsLong(PyDict_GetItemString(vis, "fqid"));
    }

    vList = PyDict_GetItemString(vis, "vis");
    len = PyList_Size(vList);
    if (len!=outUV->myDesc->ncorr) {
      PyErr_SetString(PyExc_TypeError,"UVVis incompatible with ObitUV");
      return;
    }
    visp = &outUV->buffer[outUV->myDesc->nrparm];
    for (icorr=0; icorr<len; icorr++) {
      tup = PyList_GetItem(vList, icorr);
      cx  = PyTuple_GetItem (tup, 0);
      *(visp++) = (ofloat)PyComplex_RealAsDouble(cx);
      *(visp++) = (ofloat)PyComplex_ImagAsDouble(cx);
      *(visp++) = (ofloat)PyFloat_AsDouble(PyTuple_GetItem (tup, 1));
    }

    outUV->myDesc->numVisBuff = 1;
    oretCode = ObitUVWrite (outUV, outUV->buffer, err);

} // end UVVisSet

%}

