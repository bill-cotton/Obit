%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTANTPOSPF.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTANTPOSPF (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTANTPOSPFValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTANTPOSPFGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTANTPOSPF *lTab = (ObitTableGBTANTPOSPF*)inTab;

  return outDict;
} 

extern void TableGBTANTPOSPFSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTANTPOSPF *lTab = (ObitTableGBTANTPOSPF*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEANTPOSPF;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
