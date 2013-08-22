%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTVEGASACT_STATE.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTVEGASACT_STATE (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTVEGASACT_STATEValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTVEGASACT_STATEGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTVEGASACT_STATE *lTab = (ObitTableGBTVEGASACT_STATE*)inTab;

  return outDict;
} 

extern void TableGBTVEGASACT_STATESetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTVEGASACT_STATE *lTab = (ObitTableGBTVEGASACT_STATE*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEACT_STATE;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
