/* $Id:  $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV Doppler correction utilities        */
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
#include "ObitDoppler.h"
#include "ObitUVDesc.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the st<ruct for c */
%{
typedef struct {
  ObitDoppler *me;
} Doppler;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDoppler *me;
} Doppler;

%inline %{
extern ObitDoppler* newDoppler (char* name) {
  return newObitDoppler (name);
} // end  newDoppler

extern ObitDoppler* DopplerUnref (ObitDoppler* in) {
  if (!ObitDopplerIsA(in)) return NULL;
  return ObitDopplerUnref(in);
}

extern ObitDoppler*  DopplerRef (ObitDoppler* in) {
  return ObitDopplerRef(in);
}

void CreateDoppler(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_Doppler, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateDoppler: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((Doppler*)ptr)->me = newDoppler(name);
   else  ((Doppler*)ptr)->me = NULL;
   return;
 } // end CreateDoppler

void DeleteDoppler(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_Doppler, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteDoppler: could not recover c struct");
	return;
   }
   ((Doppler*)ptr)->me = ObitDopplerUnref(((Doppler*)ptr)->me);
}// end DeleteDoppler

ObitDoppler* Doppler_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_Doppler, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Doppler_Get_me: could not recover c struct");
	return NULL;
   }
   return ((Doppler*)ptr)->me;
}  // end Doppler_Get_me

void Doppler_Set_me(PyObject *self, ObitDoppler *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_Doppler, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Doppler_Set_me: could not recover c struct");
	return;
   }
   ((Doppler*)ptr)->me =  ObitDopplerRef(me);
}  // end Doppler_Set_me

extern ObitDoppler* DopplerCopy  (ObitDoppler *in, ObitDoppler *out, 
				    ObitErr *err) {
  return ObitDopplerCopy (in, out, err);
} // end  DopplerCopy

extern ObitInfoList* DopplerGetList (ObitDoppler* in) {
  return ObitInfoListRef(in->info);
}

extern long DopplerIsA (ObitDoppler* in) {
  return (long)ObitDopplerIsA(in);
}

extern float DopplerFreqLSR (double rest, float vlsr, 
			    double ra, double dec, 
			    char *date, double ut, 
			    double x, double y, double z) {
  odouble out, JD;
  olong year, doy;
  ofloat temp;
  // munge date
  ObitUVDescDate2JD ((gchar*)date, &JD);
  ObitDopplerJD2Date (JD, &year, &doy, &temp);
  out = ObitDopplerFreqLSR (rest, vlsr, ra, dec, year, doy, ut,
                            x, y, z);
  return (float)out;
}

extern ObitUV* DopplerCVel (ObitUV *inUV, long scratch, ObitUV *outUV, 
		   	    ObitErr *err) {
  gboolean lscratch = scratch!=0;
  return ObitDopplerCVel (inUV, lscratch, outUV, err);
}

extern char* DopplerGetName (ObitDoppler* in) {
  if (ObitDopplerIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

%}

