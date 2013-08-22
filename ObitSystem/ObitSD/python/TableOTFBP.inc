%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableOTFBP.h"
%}
 
%inline %{
 
extern ObitTable* TableOTFBP (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numChan, int numPol, int numFeed,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumChan = (oint)numChan;
   oint lnumPol = (oint)numPol;
   oint lnumFeed = (oint)numFeed;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableOTFBPValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumChan, lnumPol, lnumFeed,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableOTFBPGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableOTFBP *lTab = (ObitTableOTFBP*)inTab;
  PyDict_SetItemString(outDict, "numChan",  PyInt_FromLong((long)lTab->numChan));
  PyDict_SetItemString(outDict, "numPol",  PyInt_FromLong((long)lTab->numPol));
  PyDict_SetItemString(outDict, "numFeed",  PyInt_FromLong((long)lTab->numFeed));

  return outDict;
} 

extern void TableOTFBPSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableOTFBP *lTab = (ObitTableOTFBP*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEOTFBP;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
