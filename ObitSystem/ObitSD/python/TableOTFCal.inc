%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableOTFCal.h"
%}
 
%inline %{
 
extern ObitTable* TableOTFCal (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numDet, int numPoly,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumDet = (oint)numDet;
   oint lnumPoly = (oint)numPoly;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableOTFCalValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumDet, lnumPoly,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableOTFCalGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableOTFCal *lTab = (ObitTableOTFCal*)inTab;
  PyDict_SetItemString(outDict, "numDet",  PyInt_FromLong((long)lTab->numDet));
  PyDict_SetItemString(outDict, "numPoly",  PyInt_FromLong((long)lTab->numPoly));

  return outDict;
} 

extern void TableOTFCalSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableOTFCal *lTab = (ObitTableOTFCal*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEOTFCal;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
