/* $Id$  */  
/*--------------------------------------------------------------------*/
/* Swig module description for History utilities                      */
/*                                                                    */
/*;  Copyright (C) 2004-2019                                          */
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
#include "ObitHistory.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitHistory *me;
} History;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitHistory *me;
} History;

%inline %{

extern ObitHistory* HistoryCreate (gchar* name) {
  return newObitHistory (name);
} // end  HistoryCreate

extern ObitHistory* HistoryUnref (ObitHistory* in) {
  if (!ObitHistoryIsA(in)) return NULL;
  return ObitHistoryUnref(in);
}

extern ObitHistory*  HistoryRef (ObitHistory* in) {
  return ObitHistoryRef(in);
}

void CreateHistory(PyObject *self, char* name, ObitInfoList *info, ObitErr *err) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_History, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Create: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((History*)ptr)->me = newObitHistoryValue(name, info, err);
   else ((History*)ptr)->me = NULL;
   return;
 } // end CreateHistory

void DeleteHistory(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_History, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Delete: could not recover c struct");
	return;
   }
   ((History*)ptr)->me = HistoryUnref(((History*)ptr)->me);
} // end DeleteHistory

ObitHistory* History_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_History, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Get_me: could not recover c struct");
	return NULL;
   }
   return ((History*)ptr)->me;
} // end History_Get_me

void History_Set_me(PyObject *self, ObitHistory *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_History, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"History_Set_me: could not recover c struct");
	return;
   }
   ((History*)ptr)->me =  HistoryRef(me);
}  // end History_Set_me


extern ObitHistory* HistoryZap  (ObitHistory *in, ObitErr *err) {
  return ObitHistoryZap (in, err);
} // end HistoryZap

extern ObitHistory* HistoryCopy  (ObitHistory *in, ObitHistory *out, 
			      ObitErr *err) {

  return ObitHistoryCopy (in, out, err);
} // end  HistoryCopy

extern long HistoryCopyHeader  (ObitHistory *in, ObitHistory *out, 
		 	       ObitErr *err) {

  ObitIOCode ret;
  ret = ObitHistoryCopyHeader (in, out, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  HistoryCopy

extern long HistoryCopy2Header  (ObitHistory *in, ObitHistory *out, 
		 	       ObitErr *err) {

  ObitIOCode ret;
  ret = ObitHistoryCopy2Header (in, out, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  HistoryCopy2Header

extern long HistoryHeader2Header  (ObitHistory *in, ObitHistory *out, 
		 	          ObitErr *err) {

  ObitIOCode ret;
  ret = ObitHistoryHeader2Header (in, out, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  HistoryHeader2Header

extern long HistoryOpen (ObitHistory* in, long access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitHistoryOpen (in, laccess, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
}

extern long HistoryClose (ObitHistory* in, ObitErr *err) {
  ObitIOCode ret;
  ret =  ObitHistoryClose(in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
}

extern PyObject* HistoryReadRec (ObitHistory* in, long recno,  ObitErr *err) {
  ObitIOCode ret;
  gchar hiCard[73];
  ret = ObitHistoryReadRec(in, recno, hiCard, err);
  if (ret==OBIT_IO_OK) return PyString_FromString(hiCard);
  return PyString_FromString("");
}

extern long HistoryWriteRec (ObitHistory* in, long recno, char hiCard[73],
		           ObitErr *err) {
  ObitIOCode ret;
  ret = ObitHistoryWriteRec(in, recno, hiCard, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
}

extern long HistoryEdit (ObitHistory* in, long startr, long endr,
		         ObitErr *err) {
  ObitIOCode ret;
  ret = ObitHistoryEdit(in, startr, endr, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
}

extern long HistoryTimeStamp (ObitHistory* in, char *label, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitHistoryTimeStamp(in, label, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
}

extern ObitInfoList* HistoryGetList (ObitHistory* in) {
  return ObitInfoListRef(in->info);
}

extern long HistoryIsA (ObitHistory* in) {
  return (long)ObitHistoryIsA(in);
}

extern char* HistoryGetName (ObitHistory* in) {
  if (ObitHistoryIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}
%}

