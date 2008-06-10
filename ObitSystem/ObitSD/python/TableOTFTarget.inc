%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableOTFTarget.h"
%}
 
%inline %{
 
extern ObitTable* TableOTFTarget (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableOTFTargetValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableOTFTargetGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableOTFTarget *lTab = (ObitTableOTFTarget*)inTab;
  PyDict_SetItemString(outDict, "FreqID",  PyInt_FromLong((long)lTab->FreqID));

  return outDict;
} 

extern void TableOTFTargetSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableOTFTarget *lTab = (ObitTableOTFTarget*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEOTFTarget;

  lTab->FreqID = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "FreqID"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
