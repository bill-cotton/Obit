/* $Id$  */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV  type                               */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{
#include "ObitUV.h"
#include "ObitIOUVAIPS.h"
#include "ObitIOUVFITS.h"
#include "ObitData.h"
#include "ObitUVUtil.h"
#include "ObitUVEdit.h"
#include "ObitUVWeight.h"
#include "ObitTableSNUtil.h"
#include "ObitTableNXUtil.h"
#include "ObitTableFG.h"
#include "ObitTableCLUtil.h"
%}


/* A copy of the struct for c */
%{
typedef struct {
  ObitUV *me;
} UV;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUV *me;
} UV;

%inline %{

extern void UVSetFITS(ObitUV *in, long nvis, long disk, char *file, 
                      ObitErr *err) {
  ObitUVSetFITS(in, nvis, disk, file, err);
 } // end UVSetFITS

extern void UVSetAIPS(ObitUV *in, long nvis, long disk,long cno, long user, 
                      ObitErr *err) {
  ObitUVSetAIPS(in, nvis, disk, cno, user, err);
 }

// Cast a UV data to an OData
extern ObitData* UVCastData (ObitUV* inUV) {
  return (ObitData*)inUV;
} // end  UVCastData

extern ObitUV* UVCreate (char* name) {
  return newObitUV (name);
} // end  UVCreate

void CreateUV(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UV, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateUV: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((UV*)ptr)->me = UVCreate(name);
   else  ((UV*)ptr)->me = NULL;
   return;
 } // end CreateUV

void DeleteUV(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UV, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteUV: could not recover c struct");
	return;
   }
   ((UV*)ptr)->me = ObitUVUnref(((UV*)ptr)->me);
}// end DeleteUV

ObitUV* UV_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UV, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"UV_Get_me: could not recover c struct");
	return NULL;
   }
   return ((UV*)ptr)->me;
}  // end UV_Get_me

void UV_Set_me(PyObject *self, ObitUV *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UV, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"UV_Set_me: could not recover c struct");
	return;
   }
   ((UV*)ptr)->me =  ObitUVRef(me);
}  // end UV_Set_me

extern ObitUV* UVScratch (ObitUV *in, ObitErr *err) {
  return newObitUVScratch (in, err);
} // end  UVScratch

extern PyObject* UVInfo (ObitUV *in, ObitErr *err) {
  ObitIOUVAIPS *AIPSIO=NULL;
  ObitIOUVFITS *FITSIO=NULL;
  PyObject *outDict=PyDict_New();
  PyObject *o=NULL;

  if (err->error) return outDict;

  // Ensure in fully instantiated -assume OK if myIO exists 
  if (!in->myIO) ObitUVFullInstantiate (in, TRUE, err);
  if (err->error) return outDict;

  // Get details and save in dict
  if (ObitIOUVAIPSIsA(in->myIO)) {  // AIPS
    o = PyString_InternFromString("AIPS");
    PyDict_SetItemString(outDict, "type", o);
    AIPSIO = (ObitIOUVAIPS*)in->myIO;
    o = PyLong_FromLong((long)AIPSIO->disk);
    PyDict_SetItemString(outDict, "disk", o);
    o = PyLong_FromLong((long)AIPSIO->CNO);
    PyDict_SetItemString(outDict, "CNO", o);
    o = PyLong_FromLong((long)AIPSIO->UserId);
    PyDict_SetItemString(outDict, "user", o);
  } else if (ObitIOUVFITSIsA(in->myIO)) {  // FITS
    o = PyString_InternFromString("FITS");
    PyDict_SetItemString(outDict, "type", o);
    FITSIO = (ObitIOUVFITS*)in->myIO;
    o = PyLong_FromLong((long)FITSIO->disk);
    PyDict_SetItemString(outDict, "disk", o);
    o = PyString_InternFromString((char*)FITSIO->FileName);
    PyDict_SetItemString(outDict, "filename", o);
  } else {  // Don't know this one
    o = PyString_InternFromString("UNKNOWN");
    PyDict_SetItemString(outDict, "type", o);
  }
  return outDict;
} // end  UVInfo

extern ObitUV* UVZap  (ObitUV *in, ObitErr *err) {
  return ObitUVZap (in, err);
} // end UVZap

extern void UVRename  (ObitUV *in, ObitErr *err) {
  ObitUVRename (in, err);
} // end UVRename

extern ObitUV* UVCopy  (ObitUV *in, ObitUV *out, 
			         ObitErr *err) {
  return ObitUVCopy (in, out, err);
} // end  UVCopy

extern void UVClone (ObitUV *in, ObitUV *out, ObitErr *err) {
   return  ObitUVClone (in, out, err);
} // end  UVClone

// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
// Table verion returned as outValue1
extern ObitTable* newUVTable (ObitUV *in, long access, 
			      char *tabType, long *outValue1, ObitErr *err) {
  ObitIOAccess laccess;
  olong ltabVer = (olong)*outValue1;
  ObitTable *outTable=NULL;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  outTable = newObitUVTable (in, laccess, tabType, &ltabVer, err);
  *outValue1 = (long)ltabVer;
  return outTable;
} // end  newUVTable

extern long UVZapTable (ObitUV *in, char *tabType, long tabVer, 
			ObitErr *err) {
  ObitIOCode ret;
  ret = ObitUVZapTable (in, tabType, tabVer, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  UVZapTable

extern long UVCopyTables (ObitUV *in, ObitUV *out, char **exclude,
		  	 char **include, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitUVCopyTables  (in, out, exclude, include, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  UVCopyTables

extern long UVUpdateTables (ObitUV *in, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitUVUpdateTables (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  UVUpdateTables

// Open and close to fully instantiate
// access 1=READONLY, 2=WRITEONLY, 3=READWRITE, 4=READCAL
extern long UVfullInstantiate (ObitUV* in, long access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  else if (access==4) laccess = OBIT_IO_ReadCal;
  ret = ObitUVOpen (in, laccess, err);
  ret = ObitUVClose (in, err);
  if ((err->error) || (ret!=OBIT_IO_OK)) return 1;
  else return 0;
} // end UVfullInstantiate


extern long UVOpen (ObitUV *in, long access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  else if (access==4) laccess = OBIT_IO_ReadCal;
  ret = ObitUVOpen (in, laccess, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Open

extern long UVRead (ObitUV *in, ObitErr *err) {
  ObitIOCode ret;
  if (in->myIO->access == OBIT_IO_ReadCal) 
    ret = ObitUVReadSelect (in, NULL, err);
  else 
    ret = ObitUVRead (in, NULL, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Read

extern long UVWrite (ObitUV *in, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitUVWrite (in, NULL, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Write

extern long UVRewrite (ObitUV *in, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitUVRewrite (in, NULL, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Write

extern PyObject *UVGetVisBuf (ObitUV *in){
#if PY_MAJOR_VERSION < 3
    return PyBuffer_FromReadWriteMemory(in->buffer,
                                      in->bufferSize * sizeof(ofloat));
#else
    return PyMemoryView_FromMemory((char *)in->buffer,
                                      in->bufferSize * sizeof(ofloat), PyBUF_WRITE);
#endif
}

// force header update 
extern void UVDirty (ObitUV *in) {
  in->myStatus = OBIT_Modified;
} // end Dirty

extern long UVClose (ObitUV *in, ObitErr *err) {
  ObitIOCode ret;
  ret =  ObitUVClose (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Close

extern void UVWeightData (ObitUV *in, ObitErr *err) {
  ObitUVWeightData (in, err);
} // end UVWeightData

extern void UVGetFreq (ObitUV* in, ObitErr *err) {
  ObitUVGetFreq (in, err);
} // end UVGetFreq

extern ObitIOCode UVGetSubA (ObitUV *in, ObitErr *err) {
  ObitIOCode ret;
  ret =  ObitUVGetSubA (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end UVGetSubA

extern ObitUV* UVUnref (ObitUV* in) {
  if (!ObitUVIsA(in)) return NULL;
  return ObitUVUnref(in);
}

extern ObitUV*  UVRef (ObitUV* in) {
  return ObitUVRef(in);
}

extern ObitInfoList* UVGetList (ObitUV* in) {
  return ObitInfoListRef(in->info);
}

extern ObitTableList* UVGetTableList (ObitUV* in) {
  return ObitTableListRef(in->tableList);
}

extern long UVGetHighVer (ObitUV* in, char *tabType) {
  return ObitTableListGetHigh(in->tableList, tabType);
}

extern ObitUVDesc* UVGetDesc (ObitUV* in) {
  return ObitUVDescRef(in->myDesc);
}

extern ObitUVDesc* UVGetIODesc (ObitUV* in) {
  ObitUVDesc *desc=NULL;
  if (in->myIO!=NULL) desc = (ObitUVDesc*)(in->myIO->myDesc);
  return ObitUVDescRef(desc);
}

extern void UVSetDesc (ObitUV* in, ObitUVDesc* desc) {
  in->myDesc = ObitUVDescUnref(in->myDesc);
  in->myDesc = ObitUVDescRef(desc);
}

extern long UVisScratch (ObitUV* in) {
  return (long)in->isScratch;
}

extern long UVIsA (ObitUV* in) {
  return (long)ObitUVIsA(in);
}

extern char* UVGetName (ObitUV* in) {
  if (ObitUVIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

/* UV Utility functions */

// returns float array val, [0]=maximum baseline length (in U,V), [1] = maximum W
extern void UVUtilUVWExtrema(ObitUV* in, ObitErr *err, float *val) {
  ObitUVUtilUVWExtrema (in, &val[0], &val[1], err);
} // end UVUtilUVWExtrema

extern ObitUV* UVUtilCopyZero(ObitUV* in, long scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVUtilCopyZero (in, lscratch, out, err);
} // end UVUtilCopyZero

extern void UVUtilVisDivide(ObitUV* in1, ObitUV *in2, ObitUV *out, ObitErr *err) {
  ObitUVUtilVisDivide (in1, in2, out, err);
} // end UVUtilVisDivide

extern void UVUtilXPolDivide(ObitUV* in, ObitUV *out, ObitErr *err) {
  ObitUVUtilXPolDivide (in, out, err);
} // end UVUtilXPolDivide

extern void UVUtilVisSub(ObitUV* in1, ObitUV *in2, ObitUV *out, ObitErr *err) {
  ObitUVUtilVisSub (in1, in2, out, err);
} // end UVUtilVisSub

extern float UVUtilVisCompare(ObitUV* in1, ObitUV *in2, ObitErr *err) {
  return ObitUVUtilVisCompare (in1, in2, err);
} // end UVUtilVisCompare

extern void UVUtilIndex(ObitUV* inUV, ObitErr *err) {
  return ObitUVUtilIndex (inUV, err);
} // end UVUtilVisIndex

extern void UVUtilQuack(ObitUV* inUV, float begDrop, float endDrop,
	char *Reason, long fgver, ObitErr *err) {
  olong ver;
  ofloat ftemp;
  ObitTableNX *NXTab = NULL;
  ObitTableFG *FGTab = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "UVUtilQuack";

  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  // Tables 
  ver = 1;
  NXTab = newObitTableNXValue ("Index table", (ObitData*)inUV, &ver, 
			       OBIT_IO_ReadOnly, err);
  ver = fgver;
  FGTab = newObitTableFGValue ("Flag table", (ObitData*)inUV, &ver, 
			       OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  // Editing parameters to NXTab 
  dim[0] = 1;
  ftemp = begDrop;
  ObitInfoListAlwaysPut (NXTab->info, "begDrop", OBIT_float, dim, &ftemp);
  ftemp = endDrop;
  ObitInfoListAlwaysPut (NXTab->info, "endDrop", OBIT_float, dim, &ftemp);
  dim[0] = MIN(24, strlen(Reason));
  ObitInfoListAlwaysPut (NXTab->info, "Reason", OBIT_string, dim, Reason);

  // Do editing 
  ObitTableNXUtilQuack (NXTab, FGTab, inUV->mySel, inUV->myDesc->maxAnt, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
} // end UVUtilQuack

extern ObitUV* UVUtilHann(ObitUV* in, long scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVUtilHann (in, lscratch, out, err);
} // end UVUtilHann

extern ObitUV* UVUtilAvgF(ObitUV* in, long scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVUtilAvgF (in, lscratch, out, err);
} // end UVUtilAvgF

extern ObitUV* UVUtilAvgT(ObitUV* in, long scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVUtilAvgT (in, lscratch, out, err);
} // end UVUtilAvgT

extern ObitUV* UVUtilAvgTF(ObitUV* in, long scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVUtilBlAvgTF (in, lscratch, out, err);
} // end UVUtilAvgTF

extern ObitInfoList* UVUtilCount (ObitUV *inUV, float timeInt, ObitErr *err) {
  return ObitUVUtilCount (inUV, (ofloat)timeInt, err);
} // end UVUtilCount

// How to pass array of ObitUV*???
extern void UVUtilSplitCh (ObitUV *inUV, long nOut, ObitUV **outUV, ObitErr *err) {
  return ObitUVUtilSplitCh (inUV, (olong)nOut, outUV, err);
} // end UVUtilSplitCh

extern void UVUtilNoise (ObitUV *inUV, ObitUV *outUV, 
                         float scale, float sigma, ObitErr *err) {
  return ObitUVUtilNoise (inUV, outUV, (ofloat)scale, (ofloat)sigma, err);
} // end UVUtilNoise

extern void UVUtilAppend (ObitUV *inUV, ObitUV *outUV, ObitErr *err) {
  return ObitUVUtilAppend (inUV, outUV, err);
} // end UVUtilAppend

/* UV Edit Utility functions */

extern void UVEditTD(ObitUV* in, ObitUV *out, ObitErr *err) {
  return ObitUVEditTD (in, out, err);
} // end UVEditTD

extern void UVEditFD(ObitUV* in, ObitUV *out, ObitErr *err) {
  return ObitUVEditFD (in, out, err);
} // end UVEditFD

extern void UVEditStokes(ObitUV* in, ObitUV *out, ObitErr *err) {
  return ObitUVEditStokes (in, out, err);
} // end UVEditStokes

extern ObitUV* UVEditClip(ObitUV* in, long scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVEditClip (in, lscratch, out, err);
} // end UVEditClip

extern ObitUV* UVEditClipStokes(ObitUV* in, long scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVEditClipStokes (in, lscratch, out, err);
} // end UVEditClipStokes

extern long UVUtilFlag (ObitUV *inUV, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitUVUtilFlag (inUV, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
}  // end UVUtilFlag

extern ObitUV* UVAvg2One(ObitUV* in, long scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVUtilAvg2One (in, lscratch, out, err);
} // end UVAvg2One

extern ObitUV* UVSmoF(ObitUV* in, long scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVUtilSmoF (in, lscratch, out, err);
} // end UVSmoF

extern void UVVisSub1(ObitUV* in1, ObitUV *in2, ObitUV *out, ObitErr *err) {
  ObitUVUtilVisSub1 (in1, in2, out, err);
} // end UVUtilVisSub1

/* Cal table Utility functions */

extern ObitTable* TableCLGetDummy (ObitUV *inUV,  ObitUV *outUV, long ver, ObitErr *err) {
  return (ObitTable*)ObitTableCLGetDummy (inUV, outUV, (olong)ver, err);
} // end TableCLGetDummy

extern ObitTable* TableSNGetZeroFR (ObitUV *inUV,  ObitUV *outUV, long ver, ObitErr *err) {
  return (ObitTable*)ObitTableSNGetZeroFR (inUV, outUV, (olong)ver, err);
} // end TableCLGetDummy

%}

/* Definitions for Python Shadow class */

