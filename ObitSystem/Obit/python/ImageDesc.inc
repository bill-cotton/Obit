/* $Id$   */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageDesc type                         */
/*                                                                    */
/*;  Copyright (C) 2004-2013                                          */
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
#include "ObitImageDesc.h"
%}

%inline %{
extern ObitImageDesc* ImageDescCreate (char *name) {
  return newObitImageDesc (name);
} // end ImageDescCreate

extern ObitImageDesc* ImageDescCopy (ObitImageDesc* in, 
		              ObitImageDesc* out, ObitErr *err) {
  return ObitImageDescCopy (in, out, err);
} // end ImageDescCopy

extern void ImageDescCopyDesc (ObitImageDesc* in, ObitImageDesc* out,
			ObitErr *err) {
  ObitImageDescCopyDesc  (in, out, err);
} // end ImageDescCopyDesc

extern ObitImageDesc* ImageDescDefault (char *name) {
  return ObitImageDescDefault(name);
} // end ImageDescDefault

extern void ImageDescIndex (ObitImageDesc* in) {
  ObitImageDescIndex (in);
} // end ImageDescIndex

extern int ImageDescOverlap (ObitImageDesc* in1, 
		             ObitImageDesc* in2, ObitErr *err) {
  int tout;
  tout = ObitImageDescOverlap (in1, in2, err);
  return tout;
} // end ImageDescOverlap

// returns an array of 2 floats
extern PyObject* ImageDescCvtPixel (ObitImageDesc* in, ObitImageDesc* out, 
		      float *inPixel, ObitErr *err) {
  ofloat outPixel[7];
  PyObject *outList, *o;
  ObitImageDescCvtPixel (in, out, inPixel, outPixel, err);
  outList = PyList_New(2);
  o = PyFloat_FromDouble((double)outPixel[0]);
  PyList_SetItem(outList, 0, o);
  o = PyFloat_FromDouble((double)outPixel[1]);
  PyList_SetItem(outList, 1, o);	
  return outList;
} // end ImageDescCvtPixel

// returns an array of 2 floats
extern PyObject* ImageDescGetPos (ObitImageDesc* in, 
		      float *inPixel, ObitErr *err) {
  odouble outPos[2];
  PyObject *outList, *o;
  ObitImageDescGetPos (in, (ofloat*)inPixel, outPos, err);
  outList = PyList_New(2);
  o = PyFloat_FromDouble((double)outPos[0]);
  PyList_SetItem(outList, 0, o);
  o = PyFloat_FromDouble((double)outPos[1]);
  PyList_SetItem(outList, 1, o);	
  return outList;
} // end ImageDescGetPos

// returns an array of 2 floats
extern PyObject* ImageDescGetPixel (ObitImageDesc* in, 
		      double *inPos, ObitErr *err) {
  ofloat outPixel[7];
  PyObject *outList, *o;
  ObitImageDescGetPixel (in, (odouble*)inPos, outPixel, err);
  outList = PyList_New(2);
  o = PyFloat_FromDouble((double)outPixel[0]);
  PyList_SetItem(outList, 0, o);
  o = PyFloat_FromDouble((double)outPixel[1]);
  PyList_SetItem(outList, 1, o);	
  return outList;
} // end ImageDescGetPixel

extern ObitInfoList* ImageDescGetList (ObitImageDesc* in) {
  return ObitInfoListRef(in->info);
}
 
extern PyObject *ImageDescGetDict(ObitImageDesc* in) {
  PyObject *outDict = PyDict_New();
  PyObject *list1, *list2, *list3, *list4, *list5, *list6;
  PyObject *value;
  int doh;
  int i, pos = 0;

  PyDict_SetItemString(outDict, "name",    PyString_InternFromString(in->name));
  PyDict_SetItemString(outDict, "object",  PyString_InternFromString(in->object));
  PyDict_SetItemString(outDict, "teles",   PyString_InternFromString(in->teles));
  PyDict_SetItemString(outDict, "instrume",PyString_InternFromString(in->instrument));
  PyDict_SetItemString(outDict, "observer",PyString_InternFromString(in->observer));
  PyDict_SetItemString(outDict, "origin",  PyString_InternFromString(in->origin));
  PyDict_SetItemString(outDict, "date",    PyString_InternFromString(in->date));
  PyDict_SetItemString(outDict, "obsdat",  PyString_InternFromString(in->obsdat));
  PyDict_SetItemString(outDict, "bunit",   PyString_InternFromString(in->bunit));
  PyDict_SetItemString(outDict, "obsra",   PyFloat_FromDouble((double)in->obsra));
  PyDict_SetItemString(outDict, "obsdec",  PyFloat_FromDouble((double)in->obsdec));
  PyDict_SetItemString(outDict, "epoch",   PyFloat_FromDouble((double)in->epoch));
  PyDict_SetItemString(outDict, "equinox", PyFloat_FromDouble((double)in->equinox));
  PyDict_SetItemString(outDict, "beamMaj", PyFloat_FromDouble((double)in->beamMaj));
  PyDict_SetItemString(outDict, "beamMin", PyFloat_FromDouble((double)in->beamMin));
  PyDict_SetItemString(outDict, "beamPA",  PyFloat_FromDouble((double)in->beamPA));
  PyDict_SetItemString(outDict, "minval",  PyFloat_FromDouble((double)in->minval));
  PyDict_SetItemString(outDict, "maxval",  PyFloat_FromDouble((double)in->maxval));
  PyDict_SetItemString(outDict, "xshift",  PyFloat_FromDouble((double)in->xshift));
  PyDict_SetItemString(outDict, "yshift",  PyFloat_FromDouble((double)in->yshift));
  PyDict_SetItemString(outDict, "altCrpix",PyFloat_FromDouble((double)in->altCrpix));
  PyDict_SetItemString(outDict, "altRef",  PyFloat_FromDouble((double)in->altRef));
  PyDict_SetItemString(outDict, "restFreq",PyFloat_FromDouble((double)in->restFreq));
  PyDict_SetItemString(outDict, "xpixoff", PyFloat_FromDouble((double)in->xPxOff));
  PyDict_SetItemString(outDict, "ypixoff", PyFloat_FromDouble((double)in->yPxOff));
  PyDict_SetItemString(outDict, "do3D",    PyInt_FromLong((long)in->do3D));
  PyDict_SetItemString(outDict, "areBlanks",PyInt_FromLong((long)in->areBlanks));
  PyDict_SetItemString(outDict, "niter",   PyInt_FromLong((long)in->niter));
  PyDict_SetItemString(outDict, "naxis",   PyInt_FromLong((long)in->naxis));
  PyDict_SetItemString(outDict, "bitpix",  PyInt_FromLong((long)in->bitpix));
  PyDict_SetItemString(outDict, "IOsize",  PyInt_FromLong((long)in->IOsize));
  PyDict_SetItemString(outDict, "VelReference", PyInt_FromLong((long)in->VelReference));
  PyDict_SetItemString(outDict, "VelDef", PyInt_FromLong((long)in->VelDef));
  list1 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list1, i, PyString_InternFromString(in->ctype[i]));
  PyDict_SetItemString(outDict, "ctype", list1);
  list2 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list2, i, PyFloat_FromDouble((double)in->crval[i]));
  PyDict_SetItemString(outDict, "crval", list2);
  list3 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list3, i, PyInt_FromLong((long)in->inaxes[i]));
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
  /* Structural members - mostly read only */
  PyDict_SetItemString(outDict, "jlocr",    PyInt_FromLong((long)in->jlocr));
  PyDict_SetItemString(outDict, "jlocd",    PyInt_FromLong((long)in->jlocd));
  PyDict_SetItemString(outDict, "jlocs",    PyInt_FromLong((long)in->jlocs));
  PyDict_SetItemString(outDict, "jlocf",    PyInt_FromLong((long)in->jlocf));

  /* Discard references to newly created objects. */
  doh = pos;
  while (PyDict_Next(outDict, &doh, NULL, &value))
    Py_DECREF(value);

  return outDict;
} // end ImageDescGetDict

extern void ImageDescSetDict(ObitImageDesc* in, PyObject *inDict) {
  PyObject *list1, *list2, *list3, *list4, *list5, *list6;
  char *tstr;
  int i, number;

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return;
  }

  tstr = PyString_AsString(PyDict_GetItemString(inDict, "object"));
  strncpy (in->object, tstr, IMLEN_VALUE); in->object[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obsdat"));
  strncpy (in->obsdat, tstr, IMLEN_VALUE); in->obsdat[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "bunit"));
  strncpy (in->bunit, tstr, IMLEN_VALUE); in->bunit[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "teles"));
  strncpy (in->teles, tstr, IMLEN_VALUE); in->teles[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "instrume"));
  strncpy (in->instrument, tstr, IMLEN_VALUE); in->instrument[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "observer"));
  strncpy (in->observer, tstr, IMLEN_VALUE); in->observer[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "origin"));
  strncpy (in->origin, tstr, IMLEN_VALUE); in->origin[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "date"));
  strncpy (in->date, tstr, IMLEN_VALUE); in->date[IMLEN_VALUE-1]=0;
  in->epoch   = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "epoch"));
  in->equinox = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "equinox"));
  in->obsra   = PyFloat_AsDouble(PyDict_GetItemString(inDict, "obsra"));
  in->obsdec  = PyFloat_AsDouble(PyDict_GetItemString(inDict, "obsdec"));
  in->beamMaj = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "beamMaj"));
  in->beamMin = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "beamMin"));
  in->beamPA  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "beamPA"));
  in->restFreq= (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "restFreq"));
  in->minval  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "minval"));
  in->maxval  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "maxval"));
  in->xshift  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "xshift"));
  in->yshift  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "yshift"));
  in->altCrpix= (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "altCrpix"));
  in->altRef  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "altRef"));
  in->xPxOff  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "xpixoff"));
  in->yPxOff  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ypixoff"));
  in->niter   = PyInt_AsLong(PyDict_GetItemString(inDict, "niter"));
  in->naxis   = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "naxis"));
  in->bitpix  = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "bitpix"));
  in->IOsize  = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "IOsize"));
  in->VelReference  = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "VelReference"));
  in->VelDef  = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "VelDef"));
  in->areBlanks = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "areBlanks"));
  in->do3D      = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "do3D"));
  list1 = PyDict_GetItemString(inDict, "ctype");
  number = MIN (IM_MAXDIM, PyList_Size(list1));
  for (i=0; i<number; i++) {
    tstr = PyString_AsString(PyList_GetItem(list1, i));
    strncpy (in->ctype[i], tstr, IMLEN_KEYWORD);
  }
  list2 = PyDict_GetItemString(inDict, "crval");
  number = MIN (IM_MAXDIM, PyList_Size(list2));
  for (i=0; i<number; i++) in->crval[i] = PyFloat_AsDouble(PyList_GetItem(list2, i));
  list3 = PyDict_GetItemString(inDict, "inaxes");
  number = MIN (IM_MAXDIM, PyList_Size(list3));
  for (i=0; i<number; i++) in->inaxes[i] = (long)PyInt_AsLong(PyList_GetItem(list3, i));
  list4 = PyDict_GetItemString(inDict, "cdelt");
  number = MIN (IM_MAXDIM, PyList_Size(list4));
  for (i=0; i<number; i++) in->cdelt[i] = (float)PyFloat_AsDouble(PyList_GetItem(list4, i));
  list5 = PyDict_GetItemString(inDict, "crpix");
  number = MIN (IM_MAXDIM, PyList_Size(list5));
  for (i=0; i<number; i++) in->crpix[i] = (float)PyFloat_AsDouble(PyList_GetItem(list5, i));
  list6 = PyDict_GetItemString(inDict, "crota");
  number = MIN (IM_MAXDIM, PyList_Size(list6));
  for (i=0; i<number; i++) in->crota[i] = (float)PyFloat_AsDouble(PyList_GetItem(list6, i));
  /* Reindex just to be sure */
  ObitImageDescIndex (in);
} // end ImageDescSetDict

ObitImageDesc* ImageDescRef (ObitImageDesc* in) {
  return ObitImageDescRef (in);
} // end ImageDescRef

ObitImageDesc* ImageDescUnref (ObitImageDesc* in) {
  if (!ObitImageDescIsA(in)) return NULL;
  return ObitImageDescUnref (in);
} // end ImageDescUnref

extern int ImageDescIsA (ObitImageDesc* in) {
  return ObitImageDescIsA(in);
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImageDesc *me;
} ImageDesc;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImageDesc *me;
} ImageDesc;

%addmethods ImageDesc { 
  ImageDesc(char *name) {
     ImageDesc *out;
     out = (ImageDesc *) malloc(sizeof(ImageDesc));
     if (strcmp(name, "None")) out->me = ImageDescCreate (name);
     else out->me = NULL;
     return out;
   }
  ~ImageDesc() {
    self->me = ImageDescUnref(self->me);
    free(self);
  }
};

