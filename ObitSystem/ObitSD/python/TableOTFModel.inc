%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableOTFModel.h"
%}
 
%inline %{
 
extern ObitTable* TableOTFModel (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numParm,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumParm = (oint)numParm;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableOTFModelValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumParm,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableOTFModelGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableOTFModel *lTab = (ObitTableOTFModel*)inTab;
  PyDict_SetItemString(outDict, "numParm",  PyInt_FromLong((long)lTab->numParm));

  return outDict;
} 

extern void TableOTFModelSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableOTFModel *lTab = (ObitTableOTFModel*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEOTFModel;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
