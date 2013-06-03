/* $Id: $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for Source type                            */
/*                                                                    */
/*;  Copyright (C )2013                                               */
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

%inline %{
extern ObitSource* SourceCreate (char *name) {
  return newObitSource (name);
} // end SourceCreate

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


ObitSource* SourceRef (ObitSource* in) {
  return ObitSourceRef (in);
} // end SourceRef

extern char* SourceGetName (ObitSource* in) {
  return (char*)in->SourceName;
}

extern int SourceGetSID (ObitSource* in) {
  return (int)in->SourID;
}

extern int SourceGetQual (ObitSource* in) {
  return (int)in->Qual;
}

extern int SourceGetEquinox (ObitSource* in) {
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

ObitSource* SourceUnref (ObitSource* in) {
  if (!ObitSourceIsA(in)) return NULL;
  return ObitSourceUnref (in);
} // end SourceUnref

extern int SourceIsA (ObitSource* in) {
  return ObitSourceIsA(in);
}
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

%addmethods Source { 
  Source(char *name) {
     Source *out;
     out = (Source *) malloc(sizeof(Source));
     if (strcmp(name, "None")) out->me = SourceCreate (name);
     else out->me = NULL;
     return out;
   }
  ~Source() {
    self->me = SourceUnref(self->me);
    free(self);
  }
};

