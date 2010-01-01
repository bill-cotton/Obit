%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTSPDATA.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTSPDATA (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTSPDATAValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTSPDATAGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTSPDATA *lTab = (ObitTableGBTSPDATA*)inTab;
  PyDict_SetItemString(outDict, "ctype1", PyString_InternFromString(lTab->ctype1));
  PyDict_SetItemString(outDict, "ctype2", PyString_InternFromString(lTab->ctype2));

  return outDict;
} 

extern void TableGBTSPDATASetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTSPDATA *lTab = (ObitTableGBTSPDATA*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEDATA;

  tstr = PyString_AsString(PyDict_GetItemString(inDict, "ctype1"));
  strncpy (lTab->ctype1, tstr, lstr); lTab->ctype1[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "ctype2"));
  strncpy (lTab->ctype2, tstr, lstr); lTab->ctype2[lstr-1]=0;

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
