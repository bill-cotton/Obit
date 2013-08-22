%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTVEGASPORT.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTVEGASPORT (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTVEGASPORTValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTVEGASPORTGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTVEGASPORT *lTab = (ObitTableGBTVEGASPORT*)inTab;

  return outDict;
} 

extern void TableGBTVEGASPORTSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTVEGASPORT *lTab = (ObitTableGBTVEGASPORT*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEPORT;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
