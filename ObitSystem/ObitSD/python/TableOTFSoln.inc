%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableOTFSoln.h"
%}
 
%inline %{
 
extern ObitTable* TableOTFSoln (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableOTFSolnValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumDet, lnumPoly,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableOTFSolnGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableOTFSoln *lTab = (ObitTableOTFSoln*)inTab;
  PyDict_SetItemString(outDict, "numDet",  PyInt_FromLong((long)lTab->numDet));
  PyDict_SetItemString(outDict, "numPoly",  PyInt_FromLong((long)lTab->numPoly));

  return outDict;
} 

extern void TableOTFSolnSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableOTFSoln *lTab = (ObitTableOTFSoln*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEOTFSoln;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
