%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTCCBDATA.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTCCBDATA (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTCCBDATAValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTCCBDATAGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTCCBDATA *lTab = (ObitTableGBTCCBDATA*)inTab;

  return outDict;
} 

extern void TableGBTCCBDATASetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTCCBDATA *lTab = (ObitTableGBTCCBDATA*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEDATA;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
