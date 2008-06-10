%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableOTFIndex.h"
%}
 
%inline %{
 
extern ObitTable* TableOTFIndex (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableOTFIndexValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableOTFIndexGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableOTFIndex *lTab = (ObitTableOTFIndex*)inTab;

  return outDict;
} 

extern void TableOTFIndexSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableOTFIndex *lTab = (ObitTableOTFIndex*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEOTFIndex;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
