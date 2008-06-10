%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableOTFFlag.h"
%}
 
%inline %{
 
extern ObitTable* TableOTFFlag (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableOTFFlagValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableOTFFlagGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableOTFFlag *lTab = (ObitTableOTFFlag*)inTab;

  return outDict;
} 

extern void TableOTFFlagSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableOTFFlag *lTab = (ObitTableOTFFlag*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEOTFFlag;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
