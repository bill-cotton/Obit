%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTIF.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTIF (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTIFValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTIFGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTIF *lTab = (ObitTableGBTIF*)inTab;

  return outDict;
} 

extern void TableGBTIFSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTIF *lTab = (ObitTableGBTIF*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEIF;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
