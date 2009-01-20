%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTCCBSTATE.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTCCBSTATE (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTCCBSTATEValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTCCBSTATEGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTCCBSTATE *lTab = (ObitTableGBTCCBSTATE*)inTab;

  return outDict;
} 

extern void TableGBTCCBSTATESetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTCCBSTATE *lTab = (ObitTableGBTCCBSTATE*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLESTATE;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
