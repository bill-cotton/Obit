%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTDCRRECEIVER.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTDCRRECEIVER (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTDCRRECEIVERValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTDCRRECEIVERGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTDCRRECEIVER *lTab = (ObitTableGBTDCRRECEIVER*)inTab;

  return outDict;
} 

extern void TableGBTDCRRECEIVERSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTDCRRECEIVER *lTab = (ObitTableGBTDCRRECEIVER*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLERECEIVER;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
