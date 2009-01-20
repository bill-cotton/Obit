%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTSPRECEIVER.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTSPRECEIVER (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTSPRECEIVERValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTSPRECEIVERGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTSPRECEIVER *lTab = (ObitTableGBTSPRECEIVER*)inTab;

  return outDict;
} 

extern void TableGBTSPRECEIVERSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTSPRECEIVER *lTab = (ObitTableGBTSPRECEIVER*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLERECEIVER;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
