/* $Id$   */  
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
#include "ObitTableList.h"
%}

%inline %{
extern ObitTableList* TableListCreate (char *name) {
  return newObitTableList (name);
} // end TableListCreate

extern ObitTableList* TableListCopy (ObitTableList* in, 
		              ObitTableList* out, ObitErr *err) {
  return ObitTableListCopy (in, out, err);
} // end TableListCopy

extern PyObject *TableListGetList(ObitTableList* in, ObitErr *err) {
  PyObject *outList=NULL;
  PyObject *list=NULL;
  ObitTable *table=NULL;
  gchar *tabName=NULL;
  olong tabVer, i;

  if  (err->error) return outList;

  //ObitTableListPrint(in, err);
  outList = PyList_New(in->number);
  for (i=0; i<in->number; i++) {
    ObitTableListGetNumber(in, i+1, &tabName, &tabVer, &table, err);
    list = PyList_New(2);
    PyList_SetItem(list, 0, PyInt_FromLong((long)tabVer));
    PyList_SetItem(list, 1, PyString_InternFromString(tabName));
    PyList_SetItem(outList, i, list);
    if (tabName) g_free(tabName);
    ObitTableUnref(table);
  }
  return outList;
} // end TableListGetList

long TableListGetHigh (ObitTableList* in, char *tabType) {
  return ObitTableListGetHigh (in, tabType);
} // end TableListGetHigh

void TableListPutHi (ObitTableList* in, ObitErr *err) {
  olong lversion = 1;
  gchar *tabType = "AIPS HI";
  ObitTableListPut (in, tabType, &lversion, NULL, err);
} // end TableListPutHi

ObitTableList* TableListRef (ObitTableList* in) {
  return ObitTableListRef (in);
} // end TableListRef

ObitTableList* TableListUnref (ObitTableList* in) {
  if (!ObitTableListIsA(in)) return NULL;
  return ObitTableListUnref (in);
} // end TableListUnref

extern int TableListIsA (ObitTableList* in) {
  return ObitTableListIsA(in);
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitTableList *me;
} TableList;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitTableList *me;
} TableList;

%addmethods TableList { 
  TableList(char *name) {
     TableList *out;
     out = (TableList *) malloc(sizeof(TableList));
     if (strcmp(name, "None")) out->me = TableListCreate (name);
     else out->me = NULL;
     return out;
   }
  ~TableList() {
    self->me = TableListUnref(self->me);
    free(self);
  }
};

