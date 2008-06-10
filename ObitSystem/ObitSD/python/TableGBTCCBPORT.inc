%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTCCBPORT.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTCCBPORT (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTCCBPORTValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTCCBPORTGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTCCBPORT *lTab = (ObitTableGBTCCBPORT*)inTab;

  return outDict;
} 

extern void TableGBTCCBPORTSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTCCBPORT *lTab = (ObitTableGBTCCBPORT*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLERECEIVER;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
