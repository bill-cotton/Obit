%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableOTFArrayGeom.h"
%}
 
%inline %{
 
extern ObitTable* TableOTFArrayGeom (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableOTFArrayGeomValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableOTFArrayGeomGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableOTFArrayGeom *lTab = (ObitTableOTFArrayGeom*)inTab;
  PyDict_SetItemString(outDict, "TeleX",  PyFloat_FromDouble((double)lTab->TeleX));
  PyDict_SetItemString(outDict, "TeleY",  PyFloat_FromDouble((double)lTab->TeleY));
  PyDict_SetItemString(outDict, "TeleZ",  PyFloat_FromDouble((double)lTab->TeleZ));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));
  PyDict_SetItemString(outDict, "DegDay",  PyFloat_FromDouble((double)lTab->DegDay));
  PyDict_SetItemString(outDict, "PolarX",  PyFloat_FromDouble((double)lTab->PolarX));
  PyDict_SetItemString(outDict, "PolarY",  PyFloat_FromDouble((double)lTab->PolarY));
  PyDict_SetItemString(outDict, "GSTiat0",  PyFloat_FromDouble((double)lTab->GSTiat0));
  PyDict_SetItemString(outDict, "TimeSys", PyString_InternFromString(lTab->TimeSys));

  return outDict;
} 

extern void TableOTFArrayGeomSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableOTFArrayGeom *lTab = (ObitTableOTFArrayGeom*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEOTFArrayGeom;

  lTab->TeleX = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "TeleX"));
  lTab->TeleY = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "TeleY"));
  lTab->TeleZ = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "TeleZ"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;
  lTab->DegDay = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "DegDay"));
  lTab->PolarX = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "PolarX"));
  lTab->PolarY = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "PolarY"));
  lTab->GSTiat0 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "GSTiat0"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "TimeSys"));
  strncpy (lTab->TimeSys, tstr, lstr); lTab->TimeSys[lstr-1]=0;

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
