/* $Id$  */  
/*--------------------------------------------------------------------*/
/* Swig module description for OTF type                               */
/*                                                                    */
/*;  Copyright (C) 2004-2009                                          */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{
#include "ObitOTF.h"
#include "ObitTableOTFTargetUtil.h"
#include "ObitIOOTFFITS.h"
%}


%inline %{
extern void OTFSetFITS(ObitOTF *in, int nsamp, int disk, char *file, 
                         ObitErr *err) {

  ObitOTFSetFITS(in, nsamp, disk, file, err);
 }

extern ObitOTF* OTFCreate (char* name) {
  return newObitOTF (name);
} // end  OTFCreate

extern ObitOTF* OTFScratch (ObitOTF *in, ObitErr *err) {
  return newObitOTFScratch (in, err);
} // end  OTFScratch

extern PyObject* OTFInfo (ObitOTF *in, ObitErr *err) {
  ObitIOOTFFITS *FITSIO=NULL;
  PyObject *outDict=PyDict_New();
  PyObject *o=NULL;

  if (err->error) return outDict;

  // Ensure in fully instantiated -assume OK if myIO exists 
  if (!in->myIO) ObitOTFFullInstantiate (in, TRUE, err);
  if (err->error) return outDict;

  // Get details and save in dict
  if (ObitIOOTFFITSIsA(in->myIO)) {  // FITS
    o = PyString_InternFromString("FITS");
    PyDict_SetItemString(outDict, "type", o);
    FITSIO = (ObitIOOTFFITS*)in->myIO;
    o = PyInt_FromLong((long)FITSIO->disk);
    PyDict_SetItemString(outDict, "disk", o);
    o = PyString_InternFromString((char*)FITSIO->FileName);
    PyDict_SetItemString(outDict, "filename", o);
  } else {  // Don't know this one
    o = PyString_InternFromString("UNKNOWN");
    PyDict_SetItemString(outDict, "type", o);
  }
  return outDict;
} // end  OTFInfo

extern ObitOTF* OTFZap  (ObitOTF *in, ObitErr *err) {
  return ObitOTFZap (in, err);
} // end OTFZap

extern void OTFRename  (ObitOTF *in, ObitErr *err) {
  ObitOTFRename (in, err);
} // end OTFRename

extern ObitOTF* OTFCopy  (ObitOTF *in, ObitOTF *out, 
			  ObitErr *err) {
  return ObitOTFCopy (in, out, err);
} // end OTF Copy

extern void OTFClone (ObitOTF *in, ObitOTF *out, ObitErr *err) {
   ObitOTFClone (in, out, err);
} // end  OTFClone

extern int OTFConcat (ObitOTF *in, ObitOTF *out, ObitErr *err) {
  ObitIOCode ret;
   ret = ObitOTFConcat (in, out, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  OTFClone

// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
extern int OTFOpen (ObitOTF *in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean doCalSelect;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;

  doCalSelect = FALSE;
  ObitInfoListGetTest(in->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect &&((access==1))) laccess = OBIT_IO_ReadCal;

  ret = ObitOTFOpen (in, laccess, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end OTFOpen

extern PyObject *UVGetRecBuf (ObitUV *in) {
  return PyBuffer_FromReadWriteMemory(in->buffer,
				      in->bufferSize * sizeof(ofloat));
}

// force header update 
extern void OTFDirty (ObitOTF *in) {
  in->myStatus = OBIT_Modified;
} // end Dirty

extern int OTFClose (ObitOTF *in, ObitErr *err) {
  ObitIOCode ret;
  ret =  ObitOTFClose (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end OTFClose

// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
// Table verion returned as outValue1
extern ObitTable* newOTFTable (ObitOTF *in, int access, 
			       char *tabType, long *outValue1, 
                               ObitErr *err) {
  ObitIOAccess laccess;
  olong ltabVer = (olong)*outValue1;
  ObitTable *outTable=NULL;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  outTable = newObitOTFTable (in, laccess, tabType, &ltabVer, err);
  *outValue1 = (long)ltabVer;
  return outTable;
} // end  newTable

extern int OTFZapTable (ObitOTF *in, char *tabType, long tabVer, 
			    ObitErr *err) {
  ObitIOCode ret;
  ret = ObitOTFZapTable (in, tabType, tabVer, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  ZapTable

extern int OTFCopyTables (ObitOTF *in, ObitOTF *out, char **exclude,
		  	        char **include, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitOTFCopyTables  (in, out, exclude, include, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  CopyTables

extern int OTFUpdateTables (ObitOTF *in, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitOTFUpdateTables (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  UpdateTables

// Open and close to fully instantiate
// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
extern int OTFfullInstantiate (ObitOTF* in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitOTFOpen (in, laccess, err);
  ret = ObitOTFClose (in, err);
  if ((err->error) || (ret!=OBIT_IO_OK)) return 1;
  else return 0;
} // end OTFfullInstantiate


extern void OTFSetTarget (ObitOTF *in, char *Target, 
                         float Flux, double RA, double Dec, 
                         ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess=OBIT_IO_ReadWrite;
  ObitTableOTFTarget *Table=NULL;
  olong targID[2] = {-1,-1};
  gint32 dim[MAXINFOELEMDIM];
  olong ver = 1;

  // Get Target Table
  Table = newObitTableOTFTargetValue ("TargetTable", (ObitData*)in, &ver, laccess, 
				err);
  if (err->error) return;

  // Get target ID from Target
  dim[0] = strlen(Target);
  dim[1] = dim[2] = dim[3] = 1; dim[4]=0;
  ret = ObitTableOTFTargetLookup (Table, dim, Target, targID, err);
  if ((ret!=OBIT_IO_OK) || (err->error)) return;

  // Set info
  ret = ObitTableOTFTargetSetSource (Table, targID[0], RA, Dec, Flux, err);
  if ((ret!=OBIT_IO_OK) || (err->error)) return;

  // Delete table object
  Table = ObitTableOTFTargetUnref(Table);
} // end  OTFSetTarget

extern ObitOTF* OTFUnref (ObitOTF* in) {
  if (!ObitOTFIsA(in)) return NULL;
  return ObitOTFUnref(in);
}

extern ObitOTF*  OTFRef (ObitOTF* in) {
  return ObitOTFRef(in);
}

extern ObitOTFDesc* OTFGetDesc (ObitOTF* in) {
  return ObitOTFDescRef(in->myDesc);
}

extern ObitOTFArrayGeom* OTFGetArrayGeom (ObitOTF* in) {
  return ObitOTFArrayGeomRef(in->geom);
}

extern ObitInfoList* OTFGetList (ObitOTF* in) {
  return ObitInfoListRef(in->info);
}

extern ObitTableList* OTFGetTableList (ObitOTF* in) {
  return ObitTableListRef(in->tableList);
}

extern long OTFGetHighVer (ObitOTF* in, char *tabType) {
  return ObitTableListGetHigh(in->tableList, tabType);
}

extern int OTFisScratch (ObitOTF* in) {
  return (int)in->isScratch;
}

extern int OTFIsA (ObitOTF* in) {
  return ObitOTFIsA(in);
}

extern char* OTFGetName (ObitOTF* in) {
  if (ObitOTFIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitOTF *me;
} OTF;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitOTF *me;
} OTF;

%addmethods OTF { 
  OTF(char *name) {
     OTF *out;
     out = (OTF *) malloc(sizeof(OTF));
     if (strcmp(name, "None")) out->me = OTFCreate(name);
     else  out->me = NULL;
     return out;
   }
  ~OTF() { /* Scratch files may be deleted separately*/
    if (self->me-> ReferenceCount>0) 
       self->me = OTFUnref(self->me);	
    free(self);
  }
};

