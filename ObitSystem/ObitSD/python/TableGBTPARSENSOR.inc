%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTPARSENSOR.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTPARSENSOR (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTPARSENSORValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTPARSENSORGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTPARSENSOR *lTab = (ObitTableGBTPARSENSOR*)inTab;

  return outDict;
} 

extern void TableGBTPARSENSORSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTPARSENSOR *lTab = (ObitTableGBTPARSENSOR*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLESensor;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
