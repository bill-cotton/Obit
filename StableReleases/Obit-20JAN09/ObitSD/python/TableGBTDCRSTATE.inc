%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTDCRSTATE.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTDCRSTATE (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTDCRSTATEValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTDCRSTATEGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTDCRSTATE *lTab = (ObitTableGBTDCRSTATE*)inTab;
  PyDict_SetItemString(outDict, "master", PyString_InternFromString(lTab->master));

  return outDict;
} 

extern void TableGBTDCRSTATESetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTDCRSTATE *lTab = (ObitTableGBTDCRSTATE*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLESTATE;

  tstr = PyString_AsString(PyDict_GetItemString(inDict, "master"));
  strncpy (lTab->master, tstr, lstr); lTab->master[lstr-1]=0;

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
