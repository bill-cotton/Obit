%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTVEGASSTATE.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTVEGASSTATE (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTVEGASSTATEValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTVEGASSTATEGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTVEGASSTATE *lTab = (ObitTableGBTVEGASSTATE*)inTab;

  return outDict;
} 

extern void TableGBTVEGASSTATESetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTVEGASSTATE *lTab = (ObitTableGBTVEGASSTATE*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLESTATE;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
