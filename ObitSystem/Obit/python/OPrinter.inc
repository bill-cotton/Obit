/* $Id: OPrinter.inc  $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitPrinter type                       */
/*                                                                    */
/*;  Copyright (C) 2012,2019                                          */
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
#include "ObitPrinter.h"
#include "ObitImage.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanWindow.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitPrinter *me;
} OPrinter;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitPrinter *me;
} OPrinter;


%inline %{
extern ObitPrinter* OPrinterCreate(char* name, long isInteractive, 
                                   char *streamname, char* fileName) {
  gboolean lisInteractive = isInteractive!=0;
gchar *fn=NULL;
  FILE *outStream=NULL;
  if (!strncmp(streamname, "stdout", 6)) outStream = stdout;
  if (!strncmp(streamname, "stderr", 6)) outStream = stderr;
  if (strncmp(fileName, "None", 4))     fn = (gchar*)fileName;
  return  ObitPrinterCreate ((gchar*)name, lisInteractive, outStream, fn);
} // end OPrinterCreate

ObitPrinter* OPrinterRef (ObitPrinter* in) {
  return ObitPrinterRef (in);
} // end OPrinterRef

ObitPrinter* OPrinterUnref (ObitPrinter* in) {
  if (!ObitPrinterIsA(in)) return NULL;
  return ObitPrinterUnref (in);
} // end OPrinterUnref

void CreateOPrinter(PyObject *self, char* name, long isInteractive, 
           char *streamname, char* fileName) {
   void *ptr;
   int ret, flags=0;
   gboolean lisInteractive = isInteractive!=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_OPrinter, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Create: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((OPrinter*)ptr)->me = OPrinterCreate(name,
          lisInteractive, streamname, fileName);
   else ((OPrinter*)ptr)->me = NULL;
   return;
 } // end CreateOPrinter

void DeleteOPrinter(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_OPrinter, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Delete: could not recover c struct");
	return;
   }
   ((OPrinter*)ptr)->me = OPrinterUnref(((OPrinter*)ptr)->me);
} // end DeleteOPrinter

ObitPrinter* OPrinter_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_OPrinter, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Get_me: could not recover c struct");
	return NULL;
   }
   return ((OPrinter*)ptr)->me;
} // end OPrinter_Get_me

void OPrinter_Set_me(PyObject *self, ObitPrinter *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_OPrinter, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"OPrinter_Set_me: could not recover c struct");
	return;
   }
   ((OPrinter*)ptr)->me =  OPrinterRef(me);
}  // end OPrinter_Set_me

extern void OPrinterOpen  (ObitPrinter *printer, long LinesPerPage, 
                          char *Title1, char *Title2, ObitErr *err) {
  ObitPrinterOpen (printer, (olong)LinesPerPage,
                   (gchar*)Title1, (gchar*)Title2, err);
} // end OPrinterOpen

// Returns quit as long
extern long OPrinterWrite  (ObitPrinter *printer, char *line, ObitErr *err) {
   gboolean quit;
   ObitPrinterWrite (printer, (gchar*)line, &quit, err);
  return (long) quit;
} // end OPrinterWrite

extern void OPrinterClose  (ObitPrinter *printer, ObitErr *err) {
  ObitPrinterClose(printer, err);
} // end OPrinterClose

extern long OPrinterNewPage  (ObitPrinter *printer, ObitErr *err) {
  gboolean quit;
  ObitPrinterNewPage (printer, &quit, err);
  return (long) quit;
} // end OPrinterSetNewPage

extern ObitInfoList* OPrinterGetList (ObitPrinter* in) {
  return ObitInfoListRef(in->info);
}
 
extern long OPrinterIsA (ObitPrinter* in) {
  return (long)ObitPrinterIsA(in);
} // end  OPrinterIsA 

extern char* OPrinterGetName (ObitPrinter* in) {
  return in->name;
} // end  OPrinterGetName



%}

