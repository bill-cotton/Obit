%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableSkyModel.h"
%}
 
%inline %{
 
extern ObitTable* TableSkyModel (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableSkyModelValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableSkyModelGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableSkyModel *lTab = (ObitTableSkyModel*)inTab;
  PyDict_SetItemString(outDict, "RA",  PyFloat_FromDouble((double)lTab->RA));
  PyDict_SetItemString(outDict, "Dec",  PyFloat_FromDouble((double)lTab->Dec));
  PyDict_SetItemString(outDict, "Proj", PyString_InternFromString(lTab->Proj));

  return outDict;
} 

extern void TableSkyModelSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableSkyModel *lTab = (ObitTableSkyModel*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLESkyModel;

  lTab->RA = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "RA"));
  lTab->Dec = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "Dec"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "Proj"));
  strncpy (lTab->Proj, tstr, lstr); lTab->Proj[lstr-1]=0;

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
