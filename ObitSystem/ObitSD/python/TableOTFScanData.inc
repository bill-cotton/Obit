%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableOTFScanData.h"
%}
 
%inline %{
 
extern ObitTable* TableOTFScanData (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numDet,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumDet = (oint)numDet;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableOTFScanDataValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumDet,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableOTFScanDataGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableOTFScanData *lTab = (ObitTableOTFScanData*)inTab;
  PyDict_SetItemString(outDict, "numDet",  PyInt_FromLong((long)lTab->numDet));
  PyDict_SetItemString(outDict, "beamSize",  PyFloat_FromDouble((double)lTab->beamSize));
  PyDict_SetItemString(outDict, "diameter",  PyFloat_FromDouble((double)lTab->diameter));
  PyDict_SetItemString(outDict, "OTFType", PyString_InternFromString(lTab->OTFType));

  return outDict;
} 

extern void TableOTFScanDataSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableOTFScanData *lTab = (ObitTableOTFScanData*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEOTFScanData;

  lTab->beamSize = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "beamSize"));
  lTab->diameter = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "diameter"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "OTFType"));
  strncpy (lTab->OTFType, tstr, lstr); lTab->OTFType[lstr-1]=0;

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
