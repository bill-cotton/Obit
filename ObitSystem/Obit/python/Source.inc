/* $Id: $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for Source type                            */
/*                                                                    */
/*;  Copyright (C )2013,2019                                          */
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
#include "ObitSource.h"
#include "ObitTableSUUtil.h"
#include "ObitSourceList.h"
#include "ObitUV.h"
#include "ObitPrecess.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitSource *me;
} Source;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitSource *me;
} Source;

%inline %{
extern ObitSource* SourceCreate (char *name) {
  return newObitSource (name);
} // end SourceCreate

ObitSource* SourceRef (ObitSource* in) {
  return ObitSourceRef (in);
} // end SourceRef

ObitSource* SourceUnref (ObitSource* in) {
  if (!ObitSourceIsA(in)) return NULL;
  return ObitSourceUnref (in);
} // end SourceUnref

void CreateSource(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_Source, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateSource: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((Source*)ptr)->me = SourceCreate(name);
   else  ((Source*)ptr)->me = NULL;
   return;
 } // end CreateSource

void DeleteSource(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_Source, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteSource: could not recover c struct");
	return;
   }
   ((Source*)ptr)->me = ObitSourceUnref(((Source*)ptr)->me);
}// end DeleteSource

ObitSource* Source_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_Source, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Source_Get_me: could not recover c struct");
	return NULL;
   }
   return ((Source*)ptr)->me;
}  // end Source_Get_me

void Source_Set_me(PyObject *self, ObitSource *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_Source, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Source_Set_me: could not recover c struct");
	return;
   }
   ((Source*)ptr)->me =  ObitSourceRef(me);
}  // end Source_Set_me

extern ObitSource* SourceCreateByNumber (char *name, ObitUV *inUV, long Number, ObitErr *err) {
  ObitSource *out=NULL;
  ObitSourceList *tsList=NULL;
  olong iver, i;
  ObitTableSU *SUTable=NULL;
  gchar *routine = "SourceCreateByNumber";

  // Full Source list if available
  iver = 1;
  SUTable = newObitTableSUValue (inUV->name, (ObitData*)inUV, &iver, 
				 OBIT_IO_ReadOnly, 0, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, out);
  if (SUTable) {
    tsList = ObitTableSUGetList (SUTable, err);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, out);
    // Find source and copy reference to output
    for (i=0; i<tsList->number; i++) {
      if (tsList->SUlist[i]->SourID==Number) {
        out = ObitSourceRef(tsList->SUlist[i]);
        break;
      }
    }
    tsList = ObitSourceListUnref(tsList);
    if (out==NULL) {
       Obit_log_error(err, OBIT_Error,"%s Source number %ld not found",
		   routine, Number);
       Obit_traceback_val (err, routine, inUV->name, out);
    }
  } else {  /* Use position /name from header */
    out = newObitSource (name);
    strncpy (out->SourceName, inUV->myDesc->object, MIN(20,UVLEN_VALUE));
    out->equinox = inUV->myDesc->equinox;
    out->RAMean  = inUV->myDesc->crval[inUV->myDesc->jlocr];
    out->DecMean = inUV->myDesc->crval[inUV->myDesc->jlocd];
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (inUV->myDesc, out);
    return out;
  }

  SUTable = ObitTableSUUnref(SUTable);   /* Done with table */

  return out;
} // end SourceCreateByNumber

extern ObitSource* SourceCreateByName (char *name, ObitUV *inUV, char *Name, long Qual, 
                                       ObitErr *err) {
  ObitSource *out=NULL;
  ObitSourceList *tsList=NULL;
  olong iver, i;
  ObitTableSU *SUTable=NULL;
  gchar *routine = "SourceCreateByName";

  // Full Source list if available
  iver = 1;
  SUTable = newObitTableSUValue (inUV->name, (ObitData*)inUV, &iver, 
				 OBIT_IO_ReadOnly, 0, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, out);
  if (SUTable) {
    tsList = ObitTableSUGetList (SUTable, err);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, out);
    // Find source and copy reference to output
    for (i=0; i<tsList->number; i++) {
      if ((!strncmp(Name, tsList->SUlist[i]->SourceName, MIN(20,strlen(Name))))
        && (tsList->SUlist[i]->Qual==Qual)) {
        out = ObitSourceRef(tsList->SUlist[i]);
        break;
      }
    }
    tsList = ObitSourceListUnref(tsList);
    if (out==NULL) {
       Obit_log_error(err, OBIT_Error,"%s Source %s Qual %ld not found",
		   routine, Name, Qual);
       Obit_traceback_val (err, routine, inUV->name, out);
    }
  } else {  /* Use position /name from header */
    out = newObitSource (name);
    strncpy (out->SourceName, inUV->myDesc->object, MIN(20,UVLEN_VALUE));
    out->equinox = inUV->myDesc->equinox;
    out->RAMean  = inUV->myDesc->crval[inUV->myDesc->jlocr];
    out->DecMean = inUV->myDesc->crval[inUV->myDesc->jlocd];
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (inUV->myDesc, out);
    return out;
  }

  SUTable = ObitTableSUUnref(SUTable);   /* Done with table */

  return out;
} // end SourceCreateByName

extern ObitSource* SourceCopy (ObitSource* in, 
		              ObitSource* out, ObitErr *err) {
  return ObitSourceCopy (in, out, err);
} // end SourceCopy


extern char* SourceGetName (ObitSource* in) {
  return (char*)in->SourceName;
}

extern long SourceGetSID (ObitSource* in) {
  return (long)in->SourID;
}

extern long SourceGetQual (ObitSource* in) {
  return (long)in->Qual;
}

extern float SourceGetEquinox (ObitSource* in) {
  return (float)in->equinox;
}

extern double SourceGetRAMean (ObitSource* in) {
  return (double)in->RAMean;
}

extern double SourceGetDecMean (ObitSource* in) {
  return (double)in->DecMean;
}

extern double SourceGetRAApp (ObitSource* in) {
  return (double)in->RAApp;
}

extern double SourceGetDecApp (ObitSource* in) {
  return (double)in->DecApp;
}

extern void SourceSetRAMean (ObitSource* in, double value) {
  in->RAMean = (odouble)value;
}

extern void SourceSetDecMean (ObitSource* in, double value) {
  in->DecMean = (odouble)value;
}

extern void SourceSetRAApp (ObitSource* in, double value) {
  in->RAApp = (odouble)value;
}

extern void SourceSetDecApp (ObitSource* in, double value) {
  in->DecApp = (odouble)value;
}

extern long SourceIsA (ObitSource* in) {
  return (long)ObitSourceIsA(in);
}
%}

