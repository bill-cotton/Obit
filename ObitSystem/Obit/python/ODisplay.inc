/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitDisplay type                       */
/*                                                                    */
/*;  Copyright (C) 2005-2019                                          */
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
#include "ObitDisplay.h"
#include "ObitImage.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanWindow.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitDisplay *me;
} ODisplay;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDisplay *me;
} ODisplay;


%inline %{
extern ObitDisplay* ODisplayCreate(char* name, char* serverURL, ObitErr *err) {
   return  ObitDisplayCreate ((gchar*)name, (gchar*)serverURL, err );
}

ObitDisplay* ODisplayRef (ObitDisplay* in) {
  return ObitDisplayRef (in);
} // end ODisplayRef

ObitDisplay* ODisplayUnref (ObitDisplay* in) {
  if (!ObitDisplayIsA(in)) return NULL;
  return ObitDisplayUnref (in);
} // end ODisplayUnref

void CreateODisplay(PyObject *self, char* name, char* serverURL, ObitErr *err) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ODisplay, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Create: could not recover c struct");
	return;
   }
   char *server;
   if ((serverURL==NULL) || (!strncmp(serverURL, "ObitView", 8))) 
      server = "http://localhost:8765/RPC2";
   else server = serverURL;
   if (strcmp(name, "None")) ((ODisplay*)ptr)->me = ODisplayCreate(name, server, err);
   else ((ODisplay*)ptr)->me = NULL;
   return;
 } // end CreateODisplay

void DeleteODisplay(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ODisplay, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Delete: could not recover c struct");
	return;
   }
   ((ODisplay*)ptr)->me = ODisplayUnref(((ODisplay*)ptr)->me);
} // end DeleteODisplay

ObitDisplay* ODisplay_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ODisplay, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Get_me: could not recover c struct");
	return NULL;
   }
   return ((ODisplay*)ptr)->me;
} // end ODisplay_Get_me

void ODisplay_Set_me(PyObject *self, ObitDisplay *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ODisplay, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"ODisplay_Set_me: could not recover c struct");
	return;
   }
   ((ODisplay*)ptr)->me =  ODisplayRef(me);
}  // end ODisplay_Set_me

extern long ODisplayImage  (ObitDisplay *display, ObitImage *image, ObitErr *err) {
  return (long) ObitDisplayShow (display, (Obit*)image, NULL, 1, err);
} // end ODisplayImage

extern long ODisplayMosaic  (ObitDisplay *display, ObitImageMosaic *mosaic, long field, 
                           ObitErr *err) {
  return (long) ObitDisplayShow (display, (Obit*)mosaic, NULL, (olong)field, err);
} // end ODisplayMosaic

extern long ODisplayImageEdit  (ObitDisplay *display, ObitImage *image, 
                              ObitDConCleanWindow *window, ObitErr *err) {
  return (long) ObitDisplayShow (display, (Obit*)image, window, 1, err);
} // end ODisplayImageEdit

extern long ODisplayMosaicEdit  (ObitDisplay *display, ObitImageMosaic *mosaic, long field, 
                               ObitDConCleanWindow *window, ObitErr *err) {
  return (long) ObitDisplayShow (display, (Obit*)mosaic, window, (olong)field, err);
} // end ODisplayMosaicEdit

extern long ODisplayMarkPos  (ObitDisplay *display, char *pos, ObitErr *err) {
  return (long) ObitDisplayMarkPos (display, (gchar*)pos, err);
} // end ODisplayMarkPos

extern long ODisplayIsA (ObitDisplay* in) {
  return (long)ObitDisplayIsA(in);
} // end  ODisplayIsA 

extern char* ODisplayGetName (ObitDisplay* in) {
  return in->name;
} // end  ODisplayGetName

%}


