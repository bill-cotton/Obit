%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTPARDATA.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTPARDATA (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                          
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableGBTPARDATAValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTPARDATAGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTPARDATA *lTab = (ObitTableGBTPARDATA*)inTab;
 PyDict_SetItemString(outDict, "cfgvalid",  PyInt_FromLong((long)lTab->cfgvalid));

  return outDict;
} 

extern void TableGBTPARDATASetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTPARDATA *lTab = (ObitTableGBTPARDATA*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEDATA;

  lTab->cfgvalid = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "cfgvalid"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
