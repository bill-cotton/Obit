/* $Id:  $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for parameter file parser utilities        */
/*                                                                    */
/*;  Copyright (C) 2012-2016                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
/*;  This program is free software; you can redistribute it and/or    */
/*;  modify it under the terms of the GNU General Public License as   */
/*;  published by the Free Software Foundation; either version 2 of   */
/*;  the License, or (at your option) any later version.              */
/*;                                                                   */
/*;  This program is distributed in the hope that it will be useful,  */
/*;  but WITHOUT ANY WARRANTY; without even the implied warranty of   */
/*;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    */
/*;  GNU General Public License for more details.                     */
/*;                                                                   */
/*;  You should have received a copy of the GNU General Public        */
/*;  License along with this program; if not, write to the Free       */
/*;  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,     */
/*;  MA 02139, USA.                                                   */
/*;                                                                   */
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{

#include "ObitSurveyUtil.h"
%}


%inline %{
// Print contents of a VL table
extern void OSurveyVLPrint (ObitTable *VLTable, ObitImage *image, 
  char *streamname, ObitErr *err) {
   FILE *outStream=NULL;
  if (!strncmp(streamname, "stdout", 6)) outStream = stdout;
  if (!strncmp(streamname, "stderr", 6)) outStream = stderr;
  ObitSurveyUtilVLPrint ((ObitTableVL*)VLTable, image, outStream, err);
} // end OSurveyVLPrint

// Print Selected contents of an NVSS VL table
extern int OSurveyNVSSPrint (ObitPrinter *printer, ObitData *data, 
  int VLVer, int first, int last, ObitErr *err) {
  gboolean lfirst = first!=0;
  gboolean llast  = last!=0;
  return ObitSurveyNVSSPrint (printer, data, (olong)VLVer, lfirst, llast, err);
} // end OSurveyNVSSPrint

// Print Selected contents of an VLSS VL table
extern int OSurveyVLSSPrint (ObitPrinter *printer, ObitData *data, 
  int VLVer, int first, int last, ObitErr *err) {
  gboolean lfirst = first!=0;
  gboolean llast  = last!=0;
  return ObitSurveyVLSSPrint (printer, data, (olong)VLVer, lfirst, llast, err);
} // end OSurveyVLSSPrint

// Print Selected contents of a generic VL table
extern int OSurveyGenPrint (ObitPrinter *printer, ObitData *data, 
  int VLVer, int first, int last, ObitErr *err) {
  gboolean lfirst = first!=0;
  gboolean llast  = last!=0;
  return ObitSurveyGenPrint (printer, data, (olong)VLVer, lfirst, llast, err);
} // end OSurveyGenPrint

// Get calibrated model parameters and errors
   extern PyObject *OSurveyGenCorErr(PyObject *VLrow, PyObject *calParms) {
     ObitSurveyGenCalParms *calParm = NULL;
     PyObject *outDict = PyDict_New();  
     PyObject *list;
     int i;
     ofloat fluxScale, biasRA, biasDec, calRAEr, calDecEr, ClnBiasAv, ClnBiasEr,
       calAmpEr, calSizeEr, calPolEr, beamMaj, beamMin, beamPA;
     odouble RA, Dec;
     ofloat peak, major, minor, posang, qcent, ucent, pflux, irms, prms,
       flux, eflux, epflux, pang, epang, eRA, eDec, 
       dmajor, dminor, dpa, edmajor, edminor, edpa;
     gboolean rflag[4];

     // Create calParms
     fluxScale = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "fluxScale"));
     biasRA    = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "biasRA"));
     biasDec   = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "biasDec"));
     calRAEr   = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "calRAEr"));
     calDecEr  = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "calDecEr"));
     ClnBiasAv = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "ClnBiasAv"));
     ClnBiasEr = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "ClnBiasEr"));
     calAmpEr  = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "calAmpEr"));
     calSizeEr = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "calSizeEr"));
     calPolEr  = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "calPolEr"));
     beamMaj   = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "beamMaj"));
     beamMin   = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "beamMin"));
     beamPA    = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(calParms, "beamPA"));
     calParm   = ObitSurveyGetCalParms(fluxScale, biasRA, biasDec, calRAEr, calDecEr, 
       ClnBiasAv, ClnBiasEr, calAmpEr, calSizeEr, calPolEr, beamMaj, beamMin, beamPA);

     // info from VLRow
     list   = PyDict_GetItemString(VLrow, "RA(2000)");
     RA     = PyFloat_AsDouble(PyList_GetItem(list, 0));
     list   = PyDict_GetItemString(VLrow, "DEC(2000)");
     Dec    = PyFloat_AsDouble(PyList_GetItem(list, 0));
     list   = PyDict_GetItemString(VLrow, "PEAK INT");
     peak   = (ofloat)PyFloat_AsDouble(PyList_GetItem(list, 0));
     list   = PyDict_GetItemString(VLrow, "MAJOR AX");
     major  = (ofloat)PyFloat_AsDouble(PyList_GetItem(list, 0));
     list   = PyDict_GetItemString(VLrow, "MINOR AX");
     minor  = (ofloat)PyFloat_AsDouble(PyList_GetItem(list, 0));
     list   = PyDict_GetItemString(VLrow, "POSANGLE");
     posang = (ofloat)PyFloat_AsDouble(PyList_GetItem(list, 0));
     list   = PyDict_GetItemString(VLrow, "Q CENTER");
     qcent  = (ofloat)PyFloat_AsDouble(PyList_GetItem(list, 0));
     list   = PyDict_GetItemString(VLrow, "U CENTER");
     ucent  = (ofloat)PyFloat_AsDouble(PyList_GetItem(list, 0));
     list   = PyDict_GetItemString(VLrow, "I RMS");
     irms   = (ofloat)PyFloat_AsDouble(PyList_GetItem(list, 0));
     list   = PyDict_GetItemString(VLrow, "POL RMS");
     prms   = (ofloat)PyFloat_AsDouble(PyList_GetItem(list, 0));

     // calibrate/error analysis
     ObitSurveyGenCorErr(calParm, &RA, &Dec, &peak, &major, &minor, &posang, 
       qcent, ucent, &pflux, irms, prms, &flux, &eflux, &epflux, &pang, &epang, 
       &eRA, &eDec, &dmajor, &dminor, &dpa, &edmajor, &edminor, &edpa, rflag);

     // Copy to output dict
     PyDict_SetItemString(outDict, "RA",      PyFloat_FromDouble((double)RA));
     PyDict_SetItemString(outDict, "Dec",     PyFloat_FromDouble((double)Dec));
     PyDict_SetItemString(outDict, "peak",    PyFloat_FromDouble((double)peak));
     PyDict_SetItemString(outDict, "major",   PyFloat_FromDouble((double)major));
     PyDict_SetItemString(outDict, "minor",   PyFloat_FromDouble((double)minor));
     PyDict_SetItemString(outDict, "posang",  PyFloat_FromDouble((double)posang));
     PyDict_SetItemString(outDict, "qcent",   PyFloat_FromDouble((double)qcent));
     PyDict_SetItemString(outDict, "ucent",   PyFloat_FromDouble((double)ucent));
     PyDict_SetItemString(outDict, "pflux",   PyFloat_FromDouble((double)pflux));
     PyDict_SetItemString(outDict, "irms",    PyFloat_FromDouble((double)irms));
     PyDict_SetItemString(outDict, "prms",    PyFloat_FromDouble((double)prms));
     PyDict_SetItemString(outDict, "flux",    PyFloat_FromDouble((double)flux));
     PyDict_SetItemString(outDict, "eflux",   PyFloat_FromDouble((double)eflux));
     PyDict_SetItemString(outDict, "epflux",  PyFloat_FromDouble((double)epflux));
     PyDict_SetItemString(outDict, "pang",    PyFloat_FromDouble((double)pang));
     PyDict_SetItemString(outDict, "epang",   PyFloat_FromDouble((double)epang));
     PyDict_SetItemString(outDict, "eRA",     PyFloat_FromDouble((double)eRA));
     PyDict_SetItemString(outDict, "eDec",    PyFloat_FromDouble((double)eDec));
     PyDict_SetItemString(outDict, "dmajor",  PyFloat_FromDouble((double)dmajor));
     PyDict_SetItemString(outDict, "dminor",  PyFloat_FromDouble((double)dminor));
     PyDict_SetItemString(outDict, "dpa",     PyFloat_FromDouble((double)dpa));
     PyDict_SetItemString(outDict, "edmajor", PyFloat_FromDouble((double)edmajor));
     PyDict_SetItemString(outDict, "edminor", PyFloat_FromDouble((double)edminor));
     PyDict_SetItemString(outDict, "edpa",    PyFloat_FromDouble((double)edpa));
     list = PyList_New(4);
     for (i=0; i<4; i++) PyList_SetItem(list, i, PyInt_FromLong((long)rflag[i]));
     PyDict_SetItemString(outDict, "rflag",    list);
     Py_DECREF(list);

     if (calParm) g_free(calParm);
     return outDict;
}
// end OSurveyGenCorErr


%}


