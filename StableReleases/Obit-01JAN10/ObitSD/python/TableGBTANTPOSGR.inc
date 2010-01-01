%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTANTPOSGR.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTANTPOSGR (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTANTPOSGRValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTANTPOSGRGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTANTPOSGR *lTab = (ObitTableGBTANTPOSGR*)inTab;

  return outDict;
} 

extern void TableGBTANTPOSGRSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTANTPOSGR *lTab = (ObitTableGBTANTPOSGR*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEANTPOSGR;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
