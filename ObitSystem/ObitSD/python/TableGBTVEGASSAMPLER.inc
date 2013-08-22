%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGBTVEGASSAMPLER.h"
%}
 
%inline %{
 
extern ObitTable* TableGBTVEGASSAMPLER (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableGBTVEGASSAMPLERValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGBTVEGASSAMPLERGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGBTVEGASSAMPLER *lTab = (ObitTableGBTVEGASSAMPLER*)inTab;

  return outDict;
} 

extern void TableGBTVEGASSAMPLERSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGBTVEGASSAMPLER *lTab = (ObitTableGBTVEGASSAMPLER*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEPORT;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
