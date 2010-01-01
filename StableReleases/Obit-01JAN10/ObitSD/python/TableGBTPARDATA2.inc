%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTPARDATA2.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTPARDATA2 (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTPARDATA2Value ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTPARDATA2GetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTPARDATA2 *lTab = (ObitTableGBTPARDATA2*)inTab;
 PyDict_SetItemString(outDict, "cfgvalid",  PyInt_FromLong((long)lTab->cfgvalid));

  return outDict;
} 

extern void TableGBTPARDATA2SetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTPARDATA2 *lTab = (ObitTableGBTPARDATA2*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEDATA;

  lTab->cfgvalid = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "cfgvalid"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
