/* $Id: OData.inc,v 1.1 2007/08/22 15:18:32 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitData type                          */
/*                                                                    */
/*;  Copyright (C) 2007,2008                                          */
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
#include "ObitData.h"
%}


%inline %{

extern void ODataSetFITS(ObitData *in, int disk, char *file, ObitErr *err) {
  ObitDataSetFITS(in, disk, file, err);
 } // end ODataSetFITS

extern void ODataSetAIPS(ObitData *in, int disk, int cno, int user, 
                        ObitErr *err) {
  ObitDataSetAIPS(in, disk, cno, user, err);
 }

extern ObitData* ODataScratch (ObitData *in, ObitErr *err) {
  return newObitDataScratch (in, err);
} // end  ODataScratch

extern ObitData* ODataZap  (ObitData *in, ObitErr *err) {
  return ObitDataZap (in, err);
} // end DOataZap

extern void ODataRename  (ObitData *in, ObitErr *err) {
  ObitDataRename (in, err);
} // end ODataRename

extern ObitData* ODataCopy  (ObitData *in, ObitData *out, 
			    ObitErr *err) {
  return ObitDataCopy (in, out, err);
} // end  ODataCopy

extern void ODataClone (ObitData *in, ObitData *out, ObitErr *err) {
   return  ObitDataClone (in, out, err);
} // end  ODataClone

// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
// Table version returned as outValue1
extern ObitTable* newODataTable (ObitData *in, int access, 
			      char *tabType, long *outValue1, ObitErr *err) {
  ObitIOAccess laccess;
  olong  loutValue1 = (olong)*outValue1;
  ObitTable *outTable=NULL;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  outTable =  newObitDataTable (in, laccess, tabType, &loutValue1, err);
  *outValue1 = (long)loutValue1;
  return outTable;
} // end  newDataTable

// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
extern ObitHistory* newODataHistory (ObitData *in, int access, 
			            ObitErr *err) {
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  return newObitDataHistory (in, laccess, err);
} // end  newODataTable

extern int ODataZapTable (ObitData *in, char *tabType, long tabVer, 
			 ObitErr *err) {
  ObitIOCode ret;
  ret = ObitDataZapTable (in, tabType, tabVer, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  ODataZapTable

extern int ODataCopyTables (ObitData *in, ObitData *out, char **exclude,
		  	   char **include, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitDataCopyTables  (in, out, exclude, include, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  ODataCopyTables

extern int ODataUpdateTables (ObitData *in, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitDataUpdateTables (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  ODataUpdateTables

// Open and close to fully instantiate
// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
extern int ODataFullInstantiate (ObitData* in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitDataOpen (in, laccess, err);
  ret = ObitDataClose (in, err);
  if ((err->error) || (ret!=OBIT_IO_OK)) return 1;
  else return 0;
} // end ODataFullInstantiate


extern int ODataOpen (ObitData *in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitDataOpen (in, laccess, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end OOpen

// force header update 
extern void ODataDirty (ObitData *in) {
  in->myStatus = OBIT_Modified;
} // end ODirty

extern int ODataClose (ObitData *in, ObitErr *err) {
  ObitIOCode ret;
  ret =  ObitDataClose (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end OClose

extern ObitData* ODataUnref (ObitData* in) {
  if (!ObitDataIsA(in)) return NULL;
  return ObitDataUnref(in);
}

extern ObitData*  ODataRef (ObitData* in) {
  return ObitDataRef(in);
}

extern ObitInfoList* ODataGetList (ObitData* in) {
  return ObitInfoListRef(in->info);
}

extern ObitTableList* ODataGetTableList (ObitData* in) {
  return ObitTableListRef(in->tableList);
}

extern long ODataGetHighVer (ObitData* in, char *tabType) {
  return ObitTableListGetHigh(in->tableList, tabType);
}

extern int ODataisScratch (ObitData* in) {
  return (int)in->isScratch;
}

extern int ODataIsA (ObitData* in) {
  return ObitDataIsA(in);
}

extern char* ODataGetName (ObitData* in) {
  if (ObitDataIsA(in)) {
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
  ObitData *me;
} OData;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitData *me;
} OData;

%addmethods OData { 
  OData(char *name) {
     OData *out;
     out = (OData *) malloc(sizeof(OData));
     if (strcmp(name, "None")) out->me = newObitData(name);
     else out->me = NULL;
     return out;
   }
  ~OData() { /* Scratch files may be deleted separately*/
   if (self->me->ReferenceCount>0) 
      self->me = ObitDataUnref(self->me);
   free(self);
  }
};

