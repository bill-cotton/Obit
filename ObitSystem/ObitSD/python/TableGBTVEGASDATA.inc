%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTVEGASDATA.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTVEGASDATA (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTVEGASDATAValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTVEGASDATAGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTVEGASDATA *lTab = (ObitTableGBTVEGASDATA*)inTab;

  return outDict;
} 

extern void TableGBTVEGASDATASetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTVEGASDATA *lTab = (ObitTableGBTVEGASDATA*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEDATA;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
