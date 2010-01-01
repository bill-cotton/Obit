%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTBEAM_OFFSETS.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTBEAM_OFFSETS (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTBEAM_OFFSETSValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTBEAM_OFFSETSGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTBEAM_OFFSETS *lTab = (ObitTableGBTBEAM_OFFSETS*)inTab;

  return outDict;
} 

extern void TableGBTBEAM_OFFSETSSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTBEAM_OFFSETS *lTab = (ObitTableGBTBEAM_OFFSETS*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEBEAM_OFFSETS;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
