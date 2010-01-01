%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTSPSTATE.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTSPSTATE (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTSPSTATEValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTSPSTATEGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTSPSTATE *lTab = (ObitTableGBTSPSTATE*)inTab;
  PyDict_SetItemString(outDict, "formatid", PyString_InternFromString(lTab->formatid));

  return outDict;
} 

extern void TableGBTSPSTATESetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTSPSTATE *lTab = (ObitTableGBTSPSTATE*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLESTATE;

  tstr = PyString_AsString(PyDict_GetItemString(inDict, "formatid"));
  strncpy (lTab->formatid, tstr, lstr); lTab->formatid[lstr-1]=0;

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
