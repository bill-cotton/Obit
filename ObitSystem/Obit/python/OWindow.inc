/* $Id: OWindow.inc,v 1.2 2005/12/19 00:12:10 bcotton Exp $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitDConCleanWindow type               */
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
#include "ObitErr.h"
#include "ObitDConClean.h"
#include "ObitDConCleanWindow.h"
%}

%inline %{
extern ObitDConCleanWindow* OWindowCreate (char *name,
                                           ObitImageMosaic *mosaic, 
                                           ObitErr *err) {
  return ObitDConCleanWindowCreate (name, mosaic, err);
} // end OWindowCreate

extern ObitDConCleanWindow* OWindowCreate1 (char *name,
                                            long naxis[2], 
                                            ObitErr *err) {
  olong lnaxis[2];
  lnaxis[0] = (olong)naxis[0];
  lnaxis[1] = (olong)naxis[1];
  return ObitDConCleanWindowCreate1 (name, lnaxis, err);
} // end OWindowCreate1

extern ObitDConCleanWindow* OWindowCopy (ObitDConCleanWindow* in, 
		              ObitDConCleanWindow* out, ObitErr *err) {
  return ObitDConCleanWindowCopy (in, out, err);
} // end OWindowCopy

// Convert list of lists into windows for field 
// list for each window is ID, type, followed by parameters
// previous entries removed
extern void OWindowSetList(ObitDConCleanWindow* in, PyObject *inList,
                      long field, ObitErr *err) {
  PyObject *list=NULL;
  ObitDConCleanWindowType type;
  olong lfield=field, i, n, iD, id, window[10];
  gchar *routine = "OWindowSetList";

  if  (err->error) return;

  /* Clear previous entries */
  ObitDConCleanWindowDel (in, lfield, -1, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  n = PyList_Size(inList);
  for (iD=0; iD<n; iD++) {
    list = PyList_GetItem(inList, iD);
    type = (ObitDConCleanWindowType)PyInt_AsLong(PyList_GetItem(list, 1));
    /* Ignore iD <0 */
    id = (ObitDConCleanWindowType)PyInt_AsLong(PyList_GetItem(list, 0));
    if (id<0) continue;
    switch (type) {
    case OBIT_DConCleanWindow_rectangle:
      for (i=0; i<4; i++) 
         window[i] = (olong)PyInt_AsLong(PyList_GetItem(list, i+2));
      break;
    case OBIT_DConCleanWindow_round:
      for (i=0; i<3; i++) 
         window[i] = (olong)PyInt_AsLong(PyList_GetItem(list, i+2));
      break;
    default:
      g_error ("Undefined Clean window type");
      return;
    }; /* end switch by window type */

    ObitDConCleanWindowAdd (in, lfield, type, window, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }  /* end loop */
} // end OWindowSetList

// Convert windows for field into a list of lists
// list for each window is ID, type, followed by parameters
extern PyObject* OWindowGetList(ObitDConCleanWindow* in, long field, 
                                ObitErr *err) {
  PyObject *outList=NULL;
  PyObject *list=NULL;
  ObitDConCleanWindowType type;
  olong i, n, iD, *window;

  if  (err->error) return outList;

  outList = PyList_New(in->maxId[field-1]);
  n = 0;
  for (iD=1; iD<=in->maxId[field-1]; iD++) {
    if (ObitDConCleanWindowInfo(in, field, iD, &type, &window, err)) {
      switch (type) {
      case OBIT_DConCleanWindow_rectangle:
        list = PyList_New(6);
        PyList_SetItem(list, 0, PyInt_FromLong((long)iD));
        PyList_SetItem(list, 1, PyInt_FromLong((long)type));
        for (i=0; i<4; i++) 
          PyList_SetItem(list, i+2, PyInt_FromLong((long)window[i]));
        break;
      case OBIT_DConCleanWindow_round:
        list = PyList_New(5);
        PyList_SetItem(list, 0, PyInt_FromLong((long)iD));
        PyList_SetItem(list, 1, PyInt_FromLong((long)type));
        for (i=0; i<3; i++) 
          PyList_SetItem(list, i+2, PyInt_FromLong((long)window[i]));
        break;
      default:
        g_error ("Undefined Clean window type");
        return outList;
      }; /* end switch by window type */
      PyList_SetItem(outList, n++, list);
    }  /* end if iD exists */
  } /* end loop over iD */

  /* if any more list entries fill with -1 */
  for (iD=n; iD<in->maxId[field-1]; iD++) {
    list = PyList_New(2);
    PyList_SetItem(list, 0, PyInt_FromLong((long)-1));
    PyList_SetItem(list, 1, PyInt_FromLong((long)-1));
    PyList_SetItem(outList, iD, list);
  }

  return outList;
} // end OWindowGetList

// return iD of new
extern long OWindowAdd (ObitDConCleanWindow* in, long field, long type, long *window, 
                        ObitErr *err) {
  olong lwindow[4] = {(olong)window[0],(olong)window[1],(olong)window[2],(olong)window[3]}; 
  return (long)ObitDConCleanWindowAdd (in, (olong)field, 
                                       (ObitDConCleanWindowType)type,
                                       lwindow, err);
} // end OWindowAdd

extern void OWindowUpdate (ObitDConCleanWindow* in, long field, long iD, 
                           long type, long *window, ObitErr *err) {
  olong lwindow[4] = {(olong)window[0],(olong)window[1],(olong)window[2],(olong)window[3]}; 
  ObitDConCleanWindowUpdate (in, (olong)field, (olong)iD,
                             (ObitDConCleanWindowType)type,
                             lwindow, err);
} // end OWindowUpdate

extern void OWindowDel (ObitDConCleanWindow* in, long field, long iD, ObitErr *err) {
  ObitDConCleanWindowDel (in, (olong)field, (olong)iD, err);
} // end OWindowDel

extern long OWindowGetMaxID (ObitDConCleanWindow* in, long field) {
  return (long)in->maxId[field-1];
} // end OWindowMaxID

extern ObitDConCleanWindow* OWindowRef (ObitDConCleanWindow* in) {
  return ObitDConCleanWindowRef (in);
} // end OWindowRef

extern ObitDConCleanWindow* OWindowUnref (ObitDConCleanWindow* in) {
  if (!ObitDConCleanWindowIsA(in)) return NULL;
  return ObitDConCleanWindowUnref (in);
} // end OWindowUnref

extern int OWindowIsA (ObitDConCleanWindow* in) {
  return ObitDConCleanWindowIsA(in);
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitDConCleanWindow *me;
} OWindow;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDConCleanWindow *me;
} OWindow;

%addmethods OWindow { 
  OWindow(void) {
     OWindow *out;
     out = (OWindow *) malloc(sizeof(OWindow));
     /* Defer creating Obit Object */
     out->me = NULL;
     return out;
   }
  ~OWindow() {
    self->me = OWindowUnref(self->me);
    free(self);
  }
};

