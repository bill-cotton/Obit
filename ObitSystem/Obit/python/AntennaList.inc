/* $Id: $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for AntennaList type                       */
/*                                                                    */
/*;  Copyright (C)2013,2019                                           */
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
#include "ObitAntennaList.h"
#include "ObitTableANUtil.h"
#include "ObitUV.h"
%}

/* A copy of the struct for c */
%{
typedef struct {
  ObitAntennaList *me;
} AntennaList;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitAntennaList *me;
} AntennaList;

%inline %{
extern ObitAntennaList* AntennaListCreate (char *name, ObitUV *inUV, long subA, 
    ObitErr *err) {
  ObitAntennaList *out=NULL;
  ObitTableAN *ANTable=NULL;
  olong iver;
  gchar *routine = "AntennaListCreate";

  // get Antenna List
  iver = subA;
  ANTable = newObitTableANValue (inUV->name, (ObitData*)inUV, &iver, 
				 OBIT_IO_ReadOnly, 0, 0, 0, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, out);
  if (ANTable) out = ObitTableANGetList (ANTable, err);
  ANTable = ObitTableANUnref(ANTable);   // Done with table 
  if (out==NULL) {
     Obit_log_error(err, OBIT_Error,"%s Problem generating AntennaList",
        routine);
     Obit_traceback_val (err, routine, inUV->name, out);
    }
  return out;
} // end AntennaListCreate

ObitAntennaList* AntennaListRef (ObitAntennaList* in) {
  return ObitAntennaListRef (in);
} // end AntennaListRef

ObitAntennaList* AntennaListUnref (ObitAntennaList* in) {
  return ObitAntennaListUnref (in);
} // end AntennaListUnref

void CreateAntennaList(PyObject *self, char* name, ObitUV *inUV, long subA, ObitErr *err) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_AntennaList, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateAntennaList: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((AntennaList*)ptr)->me = AntennaListCreate(name, inUV, subA, err);
   else  ((AntennaList*)ptr)->me = NULL;
   return;
 } // end CreateAntennaList

void DeleteAntennaList(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_AntennaList, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteAntennaList: could not recover c struct");
	return;
   }
   ((AntennaList*)ptr)->me = ObitAntennaListUnref(((AntennaList*)ptr)->me);
}// end DeleteAntennaList

ObitAntennaList* AntennaList_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_AntennaList, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"AntennaList_Get_me: could not recover c struct");
	return NULL;
   }
   return ((AntennaList*)ptr)->me;
}  // end AntennaList_Get_me

void AntennaList_Set_me(PyObject *self, ObitAntennaList *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_AntennaList, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"AntennaList_Set_me: could not recover c struct");
	return;
   }
   ((AntennaList*)ptr)->me =  ObitAntennaListRef(me);
}  // end AntennaList_Set_me


extern ObitAntennaList* AntennaListCopy (ObitAntennaList* in, 
		              ObitAntennaList* out, ObitErr *err) {
  return ObitAntennaListCopy (in, out, err);
} // end AntennaListCopy

extern float AntennaListGetElev (ObitAntennaList* in, long ant, float time, 
    ObitSource *Source) {
  return (float)ObitAntennaListElev (in, (olong)ant, (ofloat)time, Source);
}

extern float AntennaListGetAz (ObitAntennaList* in, long ant, float time, 
    ObitSource *Source) {
  return (float)ObitAntennaListAz (in, (olong)ant, (ofloat)time, Source);
}

extern float AntennaListGetParAng (ObitAntennaList* in, long ant, float time, 
    ObitSource *Source) {
  return (float)ObitAntennaListParAng (in, (olong)ant, (ofloat)time, Source);
}

extern double AntennaListGetRefJD (ObitAntennaList* in) {
  return (double)in->JD;
}

extern gchar* AntennaListGetArrName (ObitAntennaList* in) {
  return (gchar*)in->ArrName;
}
%}

