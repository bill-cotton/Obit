/* $Id: OASDMData.inc  $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for partial ObitSDMData type               */
/*                                                                    */
/*;  Copyright (C) 2012                                               */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{
#include "ObitSDMData.h"
#include "ObitImage.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanWindow.h"

%}


%inline %{
extern ObitSDMData* OASDMCreate(char* name, char* DataRoot, ObitErr *err) {
  ObitSDMData* asdm=NULL;
  asdm = ObitSDMIntentCreate ((gchar*)name, (gchar*) DataRoot, err);
  // basic check
  if ((asdm==NULL) || (asdm->MainTab==NULL) || (asdm->ConfigDescriptionTab==NULL) ) {
    Obit_log_error(err, OBIT_Error, "ASDM not fully populated", name);
  }

return asdm;
} // end OASDMCreate

// get reference JD
extern float OASDMGetRefJD(ObitSDMData *asdm)
{
  return (float)asdm->refJD;
}

// get array
extern char* OASDMGetArray(ObitSDMData *asdm)
{
  char *EVLA="EVLA", *ALMA="ALMA", *unk="Unknown";
  if (asdm->isEVLA) return EVLA;
  else if (asdm->isALMA) return ALMA;
  else return unk;
}

// Convert Main data to a python list of dicts
static PyObject* OASDMGetMain(ObitSDMData *asdm)
{
  ASDMMainTable* Tab =  ((ObitSDMData*)asdm)->MainTab;
  PyObject *outList, *dict;
  olong irow;

  if (!Tab) return PyList_New(0);  // bad pointer?

   outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "time",                PyFloat_FromDouble((double)Tab->rows[irow]->time));
    PyDict_SetItemString(dict, "interval",            PyFloat_FromDouble((double)Tab->rows[irow]->interval));
    PyDict_SetItemString(dict, "timeSampling",        PyInt_FromLong((long)Tab->rows[irow]->timeSampling));
    PyDict_SetItemString(dict, "scanNumber",          PyInt_FromLong((long)Tab->rows[irow]->scanNumber));
    PyDict_SetItemString(dict, "subscanNumber",       PyInt_FromLong((long)Tab->rows[irow]->subscanNumber));
    PyDict_SetItemString(dict, "numIntegration",      PyInt_FromLong((long)Tab->rows[irow]->numIntegration));
    PyDict_SetItemString(dict, "numAntenna",          PyInt_FromLong((long)Tab->rows[irow]->numAntenna));
    PyDict_SetItemString(dict, "execBlockId",         PyInt_FromLong((long)Tab->rows[irow]->execBlockId));
    PyDict_SetItemString(dict, "configDescriptionId", PyInt_FromLong((long)Tab->rows[irow]->configDescriptionId));
    PyDict_SetItemString(dict, "fieldId",             PyInt_FromLong((long)Tab->rows[irow]->fieldId));
    PyDict_SetItemString(dict, "stateId",             PyInt_FromLong((long)Tab->rows[irow]->stateId));
    if (Tab->rows[irow]->entityId)
      PyDict_SetItemString(dict, "entityId",           PyString_FromString(Tab->rows[irow]->entityId));
    else
      PyDict_SetItemString(dict, "entityId",           PyString_FromString("Unknown"));
    PyDict_SetItemString(dict, "flagRow",             PyBool_FromLong((long)Tab->rows[irow]->flagRow));	
    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetMain

// Convert Scan data to a python list of dicts
static PyObject* OASDMGetScan(ObitSDMData *asdm)
{
  ASDMScanTable* Tab =  ((ObitSDMData*)asdm)->ScanTab;
  PyObject *outList, *dict, *tlist;
  char *Unspecified="Unspecified";
  olong irow, i, cnt;

  if (!Tab) return PyList_New(0);  // bad pointer?

   outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "scanNumber",   PyInt_FromLong((long)Tab->rows[irow]->scanNumber));
    PyDict_SetItemString(dict, "numIntent",    PyInt_FromLong((long)Tab->rows[irow]->numIntent));
    PyDict_SetItemString(dict, "numSubscan",   PyInt_FromLong((long)Tab->rows[irow]->numSubscan));
    PyDict_SetItemString(dict, "execBlockId",  PyInt_FromLong((long)Tab->rows[irow]->execBlockId));
    PyDict_SetItemString(dict, "startTime",    PyFloat_FromDouble((double)Tab->rows[irow]->startTime));
    PyDict_SetItemString(dict, "endTime",      PyFloat_FromDouble((double)Tab->rows[irow]->endTime));
    PyDict_SetItemString(dict, "sourceName",   PyString_InternFromString(Tab->rows[irow]->sourceName));
    PyDict_SetItemString(dict, "flagRow",      PyBool_FromLong((long)Tab->rows[irow]->flagRow));	

    tlist = PyList_New(Tab->rows[irow]->numIntent);
    for (i=0; i<Tab->rows[irow]->numIntent; i++) {
      if (Tab->rows[irow]->scanIntent[i]!=NULL)
        PyList_SetItem(tlist, i, PyString_InternFromString(Tab->rows[irow]->scanIntent[i]));
      else
        PyList_SetItem(tlist, i, PyString_FromString(Unspecified));
    }
    PyDict_SetItemString(dict, "scanIntent",  tlist);

    cnt = 0;   /* Don't know how many, look for NULL */
    for (i=0; i<100; i++) { 
      if (Tab->rows[irow]->calDataType[i]!=NULL) cnt++;
      else break;
    }
    tlist = PyList_New(cnt);
    for (i=0; i<cnt; i++) {
      if (Tab->rows[irow]->calDataType[i]!=NULL)
        PyList_SetItem(tlist, i, PyString_InternFromString(Tab->rows[irow]->calDataType[i]));
      else
        PyList_SetItem(tlist, i, PyString_FromString(Unspecified));
    }
    PyDict_SetItemString(dict, "calDataType",  tlist);

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetScan

// Convert Subscan data to a python list of dicts
static PyObject* OASDMGetSubscan(ObitSDMData *asdm)
{
  ASDMSubscanTable* Tab =  ((ObitSDMData*)asdm)->SubscanTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i;

  if (!Tab) return PyList_New(0);  // bad pointer?

   outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "scanNumber",           PyInt_FromLong((long)Tab->rows[irow]->scanNumber));
    PyDict_SetItemString(dict, "subscanNumber",        PyInt_FromLong((long)Tab->rows[irow]->subscanNumber));
    PyDict_SetItemString(dict, "numberSubintegration", PyInt_FromLong((long)Tab->rows[irow]->numberSubintegration));
    PyDict_SetItemString(dict, "startTime",            PyFloat_FromDouble((double)Tab->rows[irow]->startTime));
    PyDict_SetItemString(dict, "endTime",              PyFloat_FromDouble((double)Tab->rows[irow]->endTime));
    PyDict_SetItemString(dict, "fieldName",            PyString_InternFromString(Tab->rows[irow]->fieldName));
    PyDict_SetItemString(dict, "subscanIntent",        PyString_InternFromString(Tab->rows[irow]->subscanIntent));
    PyDict_SetItemString(dict, "flagRow",              PyBool_FromLong((long)Tab->rows[irow]->flagRow));	

    tlist = PyList_New(Tab->rows[irow]->numberIntegration);
    for (i=0; i<Tab->rows[irow]->numberIntegration; i++) 
       PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->numberSubintegration[i]));
    PyDict_SetItemString(dict, "numberSubintegration",  tlist);

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetSubscan

// Convert Config data to a python list of dicts
static PyObject* OASDMGetConfig(ObitSDMData *asdm)
{
  ASDMConfigDescriptionTable* Tab =  ((ObitSDMData*)asdm)->ConfigDescriptionTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i;

   if (!Tab) return PyList_New(0);  // bad pointer?

   outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "configDescriptionId",   PyInt_FromLong((long)Tab->rows[irow]->configDescriptionId));
    PyDict_SetItemString(dict, "numAntenna",            PyInt_FromLong((long)Tab->rows[irow]->numAntenna));
    PyDict_SetItemString(dict, "numDataDescription",    PyInt_FromLong((long)Tab->rows[irow]->numDataDescription));
    PyDict_SetItemString(dict, "numFeed",               PyInt_FromLong((long)Tab->rows[irow]->numFeed));
    PyDict_SetItemString(dict, "numAtmPhaseCorrection", PyInt_FromLong((long)Tab->rows[irow]->numAtmPhaseCorrection));
    PyDict_SetItemString(dict, "processorId",           PyInt_FromLong((long)Tab->rows[irow]->processorId));
    PyDict_SetItemString(dict, "correlationMode",       PyInt_FromLong((long)Tab->rows[irow]->correlationMode));
    PyDict_SetItemString(dict, "processorType",         PyInt_FromLong((long)Tab->rows[irow]->processorType));
    PyDict_SetItemString(dict, "spectralType",          PyInt_FromLong((long)Tab->rows[irow]->spectralType));

    tlist = PyList_New(Tab->rows[irow]->numAntenna);
    for (i=0; i<Tab->rows[irow]->numAntenna; i++) 
       PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->antennaId[i]));
    PyDict_SetItemString(dict, "antennaId",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numDataDescription);
    for (i=0; i<Tab->rows[irow]->numDataDescription; i++) 
       PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->dataDescriptionId[i]));
    PyDict_SetItemString(dict, "dataDescriptionId",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numFeed);
    for (i=0; i<Tab->rows[irow]->numFeed; i++) 
       PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->feedId[i]));
    PyDict_SetItemString(dict, "feedId",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numAtmPhaseCorrection);
    for (i=0; i<Tab->rows[irow]->numAtmPhaseCorrection; i++) 
      PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->atmPhaseCorrection[i]));
    PyDict_SetItemString(dict, "atmPhaseCorrection",  tlist);

    /* Don't know how many
       tlist = PyList_New(cnt);
       for (i=0; i<cnt; i++) 
       PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->switchCycleId[i]));
       PyDict_SetItemString(dict, "switchCycleId",  tlist); */

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetConfig

// Convert CorrelatorMode data to a python list of dicts
static PyObject* OASDMGetCorrelatorMode(ObitSDMData *asdm)
{
  ASDMCorrelatorModeTable* Tab =  ((ObitSDMData*)asdm)->CorrelatorModeTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i;

   if (!Tab) return PyList_New(0);  // bad pointer?

   outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "correlatorModeId", PyInt_FromLong((long)Tab->rows[irow]->correlatorModeId));
    PyDict_SetItemString(dict, "correlatorName",   PyString_InternFromString(Tab->rows[irow]->correlatorName));
    PyDict_SetItemString(dict, "numBaseband",      PyInt_FromLong((long)Tab->rows[irow]->numBaseband));
    PyDict_SetItemString(dict, "numAxes",          PyInt_FromLong((long)Tab->rows[irow]->numAxes));
    PyDict_SetItemString(dict, "accumMode",        PyInt_FromLong((long)Tab->rows[irow]->accumMode));
    PyDict_SetItemString(dict, "binMode",          PyInt_FromLong((long)Tab->rows[irow]->binMode));

    tlist = PyList_New(Tab->rows[irow]->numBaseband);
    for (i=0; i<Tab->rows[irow]->numBaseband; i++) 
      PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->basebandNames[i]));
    PyDict_SetItemString(dict, "basebandNames",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numBaseband);
    for (i=0; i<Tab->rows[irow]->numBaseband; i++) 
      PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->basebandConfig[i]));
    PyDict_SetItemString(dict, "basebandConfig",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numBaseband);
    for (i=0; i<Tab->rows[irow]->numBaseband; i++) 
      PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->filterMode[i]));
    PyDict_SetItemString(dict, "filterMode",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numAxes);
    for (i=0; i<Tab->rows[irow]->numAxes; i++) 
      PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->axesOrderArray[i]));
    PyDict_SetItemString(dict, "axesOrderArray",  tlist);

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetCorrelatorMode

// Convert DataDescription  data to a python list of dicts
static PyObject* OASDMGetDataDescription(ObitSDMData *asdm)
{
  ASDMDataDescriptionTable* Tab =  ((ObitSDMData*)asdm)->DataDescriptionTab;
  PyObject *outList, *dict;
  olong irow;

  if (!Tab) return PyList_New(0);  // bad pointer?

  outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "dataDescriptionId", PyInt_FromLong((long)Tab->rows[irow]->dataDescriptionId));
    PyDict_SetItemString(dict, "polOrHoloId",       PyInt_FromLong((long)Tab->rows[irow]->polOrHoloId));
    PyDict_SetItemString(dict, "spectralWindowId",  PyInt_FromLong((long)Tab->rows[irow]->spectralWindowId));
    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetDataDescription

// Convert SpectralWindow data to a python list of dicts
static PyObject* OASDMGetSpectralWindow(ObitSDMData *asdm)
{
  ASDMSpectralWindowTable* Tab =  ((ObitSDMData*)asdm)->SpectralWindowTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i;

  if (!Tab) return PyList_New(0);  // bad pointer?

  outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "spectralWindowId",       PyInt_FromLong((long)Tab->rows[irow]->spectralWindowId));
    PyDict_SetItemString(dict, "basebandName",           PyString_InternFromString(Tab->rows[irow]->basebandName));
    PyDict_SetItemString(dict, "netSideband",            PyString_InternFromString(Tab->rows[irow]->netSideband));
    PyDict_SetItemString(dict, "numChan",                PyInt_FromLong((long)Tab->rows[irow]->numChan));
    PyDict_SetItemString(dict, "numAssocValues",         PyInt_FromLong((long)Tab->rows[irow]->numAssocValues));
    PyDict_SetItemString(dict, "sidebandProcessingMode", PyInt_FromLong((long)Tab->rows[irow]->sidebandProcessingMode));
    PyDict_SetItemString(dict, "refFreq",                PyFloat_FromDouble((double)Tab->rows[irow]->refFreq));
    PyDict_SetItemString(dict, "totBandwidth",           PyFloat_FromDouble((double)Tab->rows[irow]->totBandwidth));
    PyDict_SetItemString(dict, "chanFreqStart",          PyFloat_FromDouble((double)Tab->rows[irow]->chanFreqStart));
    PyDict_SetItemString(dict, "chanFreqStep",           PyFloat_FromDouble((double)Tab->rows[irow]->chanFreqStep));
    PyDict_SetItemString(dict, "chanWidth",              PyFloat_FromDouble((double)Tab->rows[irow]->chanWidth));
    PyDict_SetItemString(dict, "effectiveBw",            PyFloat_FromDouble((double)Tab->rows[irow]->effectiveBw));
    PyDict_SetItemString(dict, "resolution",             PyFloat_FromDouble((double)Tab->rows[irow]->resolution));
    PyDict_SetItemString(dict, "windowFunction",         PyInt_FromLong((long)Tab->rows[irow]->windowFunction));
    PyDict_SetItemString(dict, "correlationBit",         PyString_InternFromString(Tab->rows[irow]->correlationBit));
    PyDict_SetItemString(dict, "name",                   PyString_InternFromString(Tab->rows[irow]->name));
    PyDict_SetItemString(dict, "oversampling",           PyBool_FromLong((long)Tab->rows[irow]->oversampling));	
    PyDict_SetItemString(dict, "quantization",           PyBool_FromLong((long)Tab->rows[irow]->quantization));	

    if (Tab->rows[irow]->chanFreqArray!=NULL) {
      tlist = PyList_New(Tab->rows[irow]->numChan);
      for (i=0; i<Tab->rows[irow]->numChan; i++) 
         PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->chanFreqArray[i]));
      PyDict_SetItemString(dict, "chanFreqArray",  tlist);
    }

    if (Tab->rows[irow]->chanWidthArray!=NULL) {
      tlist = PyList_New(Tab->rows[irow]->numChan);
      for (i=0; i<Tab->rows[irow]->numChan; i++) 
         PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->chanWidthArray[i]));
      PyDict_SetItemString(dict, "chanWidthArray",  tlist);
    }

    if (Tab->rows[irow]->effectiveBwArray!=NULL) {
      tlist = PyList_New(Tab->rows[irow]->numChan);
      for (i=0; i<Tab->rows[irow]->numChan; i++) 
         PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->effectiveBwArray[i]));
      PyDict_SetItemString(dict, "effectiveBwArray",  tlist);
    }

    if (Tab->rows[irow]->resolutionArray!=NULL) {
      tlist = PyList_New(Tab->rows[irow]->numChan);
      for (i=0; i<Tab->rows[irow]->numChan; i++) 
         PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->resolutionArray[i]));
      PyDict_SetItemString(dict, "resolutionArray",  tlist);
    }

    if (Tab->rows[irow]->SpecRes!=NULL) {
      tlist = PyList_New(Tab->rows[irow]->numAssocValues);
      for (i=0; i<Tab->rows[irow]->numAssocValues; i++) 
         PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->SpecRes[i]));
      PyDict_SetItemString(dict, "SpecRes",  tlist);
    }

    if (Tab->rows[irow]->assocSpectralWindowId!=NULL) {
      tlist = PyList_New(Tab->rows[irow]->numAssocValues);
      for (i=0; i<Tab->rows[irow]->numAssocValues; i++) 
         PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->assocSpectralWindowId[i]));
      PyDict_SetItemString(dict, "assocSpectralWindowId",  tlist);
    }
    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetSpectralWindow

// Convert Antenna data to a python list of dicts
static PyObject* OASDMGetAntenna(ObitSDMData *asdm)
{
  ASDMAntennaTable* Tab =  ((ObitSDMData*)asdm)->AntennaTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i;

  if (!Tab) return PyList_New(0);  // bad pointer?

   outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "antennaId",    PyInt_FromLong((long)Tab->rows[irow]->antennaId));
    PyDict_SetItemString(dict, "stationId",    PyInt_FromLong((long)Tab->rows[irow]->stationId));
    PyDict_SetItemString(dict, "name",         PyString_InternFromString(Tab->rows[irow]->name));
    PyDict_SetItemString(dict, "antennaMake",  PyInt_FromLong((long)Tab->rows[irow]->antennaMake));
    PyDict_SetItemString(dict, "antennaType",  PyInt_FromLong((long)Tab->rows[irow]->antennaType));
    PyDict_SetItemString(dict, "dishDiameter", PyFloat_FromDouble((double)Tab->rows[irow]->dishDiameter));
    PyDict_SetItemString(dict, "time",         PyFloat_FromDouble((double)Tab->rows[irow]->time));

    tlist = PyList_New(3);
    for (i=0; i<3; i++) 
       PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->position[i]));
    PyDict_SetItemString(dict, "position",  tlist);

    tlist = PyList_New(3);
    for (i=0; i<3; i++) 
       PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->offset[i]));
    PyDict_SetItemString(dict, "offset",  tlist);

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetAntenna

// Convert Station data to a python list of dicts
static PyObject* OASDMGetStation(ObitSDMData *asdm)
{
  ASDMStationTable* Tab =  ((ObitSDMData*)asdm)->StationTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i;

  if (!Tab) return PyList_New(0);  // bad pointer?

  outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "stationId", PyInt_FromLong((long)Tab->rows[irow]->stationId));
    PyDict_SetItemString(dict, "type",      PyInt_FromLong((long)Tab->rows[irow]->type));
    PyDict_SetItemString(dict, "name",      PyString_InternFromString(Tab->rows[irow]->name));

    tlist = PyList_New(3);
    for (i=0; i<3; i++) 
       PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->position[i]));
    PyDict_SetItemString(dict, "position",  tlist);

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetStation

// Convert State data to a python list of dicts
static PyObject* OASDMGetState(ObitSDMData *asdm)
{
  ASDMStateTable* Tab =  ((ObitSDMData*)asdm)->StateTab;
  PyObject *outList, *dict;
  olong irow;

  if (!Tab) return PyList_New(0);  // bad pointer?

  outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "stateId",       PyInt_FromLong((long)Tab->rows[irow]->stateId));
    PyDict_SetItemString(dict, "calDeviceName", PyString_InternFromString(Tab->rows[irow]->calDeviceName));
    PyDict_SetItemString(dict, "sig",           PyBool_FromLong((long)Tab->rows[irow]->sig));	
    PyDict_SetItemString(dict, "ref",           PyBool_FromLong((long)Tab->rows[irow]->ref));	
    PyDict_SetItemString(dict, "onSky",         PyBool_FromLong((long)Tab->rows[irow]->onSky));	

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetState

// Convert ExecBlock data to a python list of dicts
static PyObject* OASDMGetExecBlock(ObitSDMData *asdm)
{
  ASDMExecBlockTable* Tab =  ((ObitSDMData*)asdm)->ExecBlockTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i;

  if (!Tab) return PyList_New(0);  // bad pointer?

  outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "execBlockId",      PyInt_FromLong((long)Tab->rows[irow]->execBlockId));
    PyDict_SetItemString(dict, "execBlockNum",     PyInt_FromLong((long)Tab->rows[irow]->execBlockNum));
    PyDict_SetItemString(dict, "numObservingLog",  PyInt_FromLong((long)Tab->rows[irow]-> numObservingLog));
    PyDict_SetItemString(dict, "numAntenna",       PyInt_FromLong((long)Tab->rows[irow]->numAntenna));
    PyDict_SetItemString(dict, "sbSummaryId;",     PyInt_FromLong((long)Tab->rows[irow]->sbSummaryId));
    PyDict_SetItemString(dict, "ScaleId",          PyInt_FromLong((long)Tab->rows[irow]->ScaleId));
    PyDict_SetItemString(dict, "flagRow",          PyBool_FromLong((long)Tab->rows[irow]->flagRow));	
    PyDict_SetItemString(dict, "configName",       PyString_InternFromString(Tab->rows[irow]->configName));
    PyDict_SetItemString(dict, "telescopeName",    PyString_InternFromString(Tab->rows[irow]->telescopeName));
    PyDict_SetItemString(dict, "observerName",     PyString_InternFromString(Tab->rows[irow]->observerName));
    PyDict_SetItemString(dict, "sessionReference", PyString_InternFromString(Tab->rows[irow]->sessionReference));
    PyDict_SetItemString(dict, "schedulerMode",    PyString_InternFromString(Tab->rows[irow]->schedulerMode));
    PyDict_SetItemString(dict, "baseRangeMin",     PyFloat_FromDouble((double)Tab->rows[irow]->baseRangeMin));
    PyDict_SetItemString(dict, "baseRangeMax",     PyFloat_FromDouble((double)Tab->rows[irow]->baseRangeMax));
    PyDict_SetItemString(dict, "baseRmsMinor",     PyFloat_FromDouble((double)Tab->rows[irow]->baseRmsMinor));
    PyDict_SetItemString(dict, "baseRmsMajor",     PyFloat_FromDouble((double)Tab->rows[irow]->baseRmsMajor));
    PyDict_SetItemString(dict, "basePa",           PyFloat_FromDouble((double)Tab->rows[irow]->basePa));
    PyDict_SetItemString(dict, "siteAltitude",     PyFloat_FromDouble((double)Tab->rows[irow]->siteAltitude));
    PyDict_SetItemString(dict, "siteLongitude",    PyFloat_FromDouble((double)Tab->rows[irow]->siteLongitude));
    PyDict_SetItemString(dict, "siteLatitude",     PyFloat_FromDouble((double)Tab->rows[irow]->siteLatitude));
    PyDict_SetItemString(dict, "aborted",          PyBool_FromLong((long)Tab->rows[irow]->aborted));	

    tlist = PyList_New(Tab->rows[irow]->numAntenna);
    for (i=0; i<Tab->rows[irow]->numAntenna; i++) 
       PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->antennaId[i]));
    PyDict_SetItemString(dict, "antennaId",  tlist);

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetExecBlock

// Convert Source data to a python list of dicts
static PyObject* OASDMGetSource(ObitSDMData *asdm)
{
  ASDMSourceTable* Tab =  ((ObitSDMData*)asdm)->SourceTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i;

  if (!Tab) return PyList_New(0);  // bad pointer?

  outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "sourceId",         PyInt_FromLong((long)Tab->rows[irow]->sourceId));
    PyDict_SetItemString(dict, "sourceNo",         PyInt_FromLong((long)Tab->rows[irow]->sourceNo));
    PyDict_SetItemString(dict, "sourceName",       PyString_InternFromString(Tab->rows[irow]->sourceName));
    PyDict_SetItemString(dict, "numLines",         PyInt_FromLong((long)Tab->rows[irow]->numLines));
    PyDict_SetItemString(dict, "spectralWindowId", PyInt_FromLong((long)Tab->rows[irow]->spectralWindowId));
    PyDict_SetItemString(dict, "code",             PyString_InternFromString(Tab->rows[irow]->code));

    tlist = PyList_New(Tab->rows[irow]->numLines);
    for (i=0; i<Tab->rows[irow]->numLines; i++) 
       PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->restFrequency[i]));
    PyDict_SetItemString(dict, "restFrequency",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numLines);
    for (i=0; i<Tab->rows[irow]->numLines; i++) 
       PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->sysVel[i]));
    PyDict_SetItemString(dict, "sysVel",  tlist);

    tlist = PyList_New(2);
    for (i=0; i<2; i++) 
       PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->direction[i]));
    PyDict_SetItemString(dict, "direction",  tlist);

    tlist = PyList_New(2);
    for (i=0; i<2; i++) 
       PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->properMotion[i]));
    PyDict_SetItemString(dict, "properMotion",  tlist);

    tlist = PyList_New(2);
    PyList_SetItem(tlist, 0, PyFloat_FromDouble((double)Tab->rows[irow]->timeInterval[0]));
    PyList_SetItem(tlist, 1, PyFloat_FromDouble((double)Tab->rows[irow]->timeInterval[1]));
    PyDict_SetItemString(dict, "timeInterval",  tlist);

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetSource

// Convert Field data to a python list of dicts
static PyObject* OASDMGetField(ObitSDMData *asdm)
{
  ASDMFieldTable* Tab =  ((ObitSDMData*)asdm)->FieldTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i;

  if (!Tab) return PyList_New(0);  // bad pointer?

  outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "fieldId",   PyInt_FromLong((long)Tab->rows[irow]->fieldId));
    PyDict_SetItemString(dict, "sourceId",  PyInt_FromLong((long)Tab->rows[irow]->sourceId));
    PyDict_SetItemString(dict, "fieldName", PyString_InternFromString(Tab->rows[irow]->fieldName));
    PyDict_SetItemString(dict, "code",      PyString_InternFromString(Tab->rows[irow]->code));
    PyDict_SetItemString(dict, "time",      PyFloat_FromDouble((double)Tab->rows[irow]->time));
    PyDict_SetItemString(dict, "numPoly",   PyInt_FromLong((long)Tab->rows[irow]->numPoly));

    tlist = PyList_New(2*Tab->rows[irow]->numPoly);
    for (i=0; i<Tab->rows[irow]->numPoly; i+=2) {
       PyList_SetItem(tlist, 2*i,   PyFloat_FromDouble((double)Tab->rows[irow]->delayDir[2*i]));
       PyList_SetItem(tlist, 2*i+1, PyFloat_FromDouble((double)Tab->rows[irow]->delayDir[2*i+1]));
    }
    PyDict_SetItemString(dict, "delayDir",  tlist);

    tlist = PyList_New(2*Tab->rows[irow]->numPoly);
    for (i=0; i<Tab->rows[irow]->numPoly; i+=2) {
       PyList_SetItem(tlist, 2*i,   PyFloat_FromDouble((double)Tab->rows[irow]->phaseDir[2*i]));
       PyList_SetItem(tlist, 2*i+1, PyFloat_FromDouble((double)Tab->rows[irow]->phaseDir[2*i+1]));
    }
    PyDict_SetItemString(dict, "phaseDir",  tlist);

    tlist = PyList_New(2*Tab->rows[irow]->numPoly);
    for (i=0; i<Tab->rows[irow]->numPoly; i+=2) {
       PyList_SetItem(tlist, 2*i,   PyFloat_FromDouble((double)Tab->rows[irow]->referenceDir[2*i]));
       PyList_SetItem(tlist, 2*i+1, PyFloat_FromDouble((double)Tab->rows[irow]->referenceDir[2*i+1]));
    }
    PyDict_SetItemString(dict, "referenceDir",  tlist);

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetField

// Convert Feed data to a python list of dicts
static PyObject* OASDMGetFeed(ObitSDMData *asdm)
{
  ASDMFeedTable* Tab =  ((ObitSDMData*)asdm)->FeedTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i;

  if (!Tab) return PyList_New(0);  // bad pointer?

  outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "receiverId",       PyInt_FromLong((long)Tab->rows[irow]->receiverId));
    PyDict_SetItemString(dict, "antennaId",        PyInt_FromLong((long)Tab->rows[irow]->antennaId));
    PyDict_SetItemString(dict, "spectralWindowId", PyInt_FromLong((long)Tab->rows[irow]->spectralWindowId));
    PyDict_SetItemString(dict, "numReceptor",      PyInt_FromLong((long)Tab->rows[irow]->numReceptor));

    tlist = PyList_New(2);
    PyList_SetItem(tlist, 0, PyFloat_FromDouble((double)Tab->rows[irow]->timeInterval[0]));
    PyList_SetItem(tlist, 1, PyFloat_FromDouble((double)Tab->rows[irow]->timeInterval[1]));
    PyDict_SetItemString(dict, "timeInterval",  tlist);

    tlist = PyList_New(2);
    PyList_SetItem(tlist, 0, PyFloat_FromDouble((double)Tab->rows[irow]->beamOffset[0]));
    PyList_SetItem(tlist, 1, PyFloat_FromDouble((double)Tab->rows[irow]->beamOffset[1]));
    PyDict_SetItemString(dict, "beamOffset",  tlist);

    tlist = PyList_New(2);
    PyList_SetItem(tlist, 0, PyFloat_FromDouble((double)Tab->rows[irow]->focusReference[0]));
    PyList_SetItem(tlist, 1, PyFloat_FromDouble((double)Tab->rows[irow]->focusReference[1]));
    PyDict_SetItemString(dict, "focusReference",  tlist);

    tlist = PyList_New(2);
    PyList_SetItem(tlist, 0, PyInt_FromLong((long)Tab->rows[irow]->polResponse[0]));
    PyList_SetItem(tlist, 1, PyInt_FromLong((long)Tab->rows[irow]->polResponse[1]));
    PyDict_SetItemString(dict, "polResponse",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numReceptor);
    for (i=0; i<Tab->rows[irow]->numReceptor; i++) 
       PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->polarizationTypes[i]));
    PyDict_SetItemString(dict, "polarizationTypes",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numReceptor);
    for (i=0; i<Tab->rows[irow]->numReceptor; i++) 
       PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->receptorAngle[i]));
    PyDict_SetItemString(dict, "receptorAngle",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numReceptor);
    for (i=0; i<Tab->rows[irow]->numReceptor; i++) 
       PyList_SetItem(tlist, i, PyInt_FromLong((long)Tab->rows[irow]->receiverId[i]));
    PyDict_SetItemString(dict, "receiverId",  tlist);

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetFeed

// Convert Receiver data to a python list of dicts
static PyObject* OASDMGetReceiver(ObitSDMData *asdm)
{
  ASDMReceiverTable* Tab =  ((ObitSDMData*)asdm)->ReceiverTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i, n;

  if (!Tab) return PyList_New(0);  // bad pointer?

  outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "receiverId",       PyInt_FromLong((long)Tab->rows[irow]->receiverId));
    PyDict_SetItemString(dict, "spectralWindowId", PyInt_FromLong((long)Tab->rows[irow]->spectralWindowId));
    PyDict_SetItemString(dict, "name",             PyString_InternFromString(Tab->rows[irow]->name));
    PyDict_SetItemString(dict, "receiverSideband", PyString_InternFromString(Tab->rows[irow]->receiverSideband));
    PyDict_SetItemString(dict, "frequencyBand",    PyString_InternFromString(Tab->rows[irow]->frequencyBand));
    PyDict_SetItemString(dict, "numLO",            PyInt_FromLong((long)Tab->rows[irow]->numLO));

    tlist = PyList_New(2);
    PyList_SetItem(tlist, 0, PyFloat_FromDouble((double)Tab->rows[irow]->timeInterval[0]));
    PyList_SetItem(tlist, 1, PyFloat_FromDouble((double)Tab->rows[irow]->timeInterval[1]));
    PyDict_SetItemString(dict, "timeInterval",  tlist);

    n = Tab->rows[irow]->numLO;
    tlist = PyList_New(n);
    for (i=0; i<n; i++) {
      PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->freqLO[i]));
    }
    PyDict_SetItemString(dict, "freqLO",  tlist);

    tlist = PyList_New(n);
    for (i=0; i<n; i++) {
      PyList_SetItem(tlist, i, PyString_InternFromString(Tab->rows[irow]->sidebandLO[i]));
    }
    PyDict_SetItemString(dict, "sidebandLO",  tlist);

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetReceiver

// Convert Polarization data to a python list of dicts
static PyObject* OASDMGetPolarization(ObitSDMData *asdm)
{
  ASDMPolarizationTable* Tab =  ((ObitSDMData*)asdm)->PolarizationTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i;

  if (!Tab) return PyList_New(0);  // bad pointer?

  outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "polarizationId", PyInt_FromLong((long)Tab->rows[irow]->polarizationId));
    PyDict_SetItemString(dict, "numCorr",        PyInt_FromLong((long)Tab->rows[irow]->numCorr));
    PyDict_SetItemString(dict, "flagRow",        PyBool_FromLong((long)Tab->rows[irow]->flagRow));	

    tlist = PyList_New(Tab->rows[irow]->numCorr);
    for (i=0; i<Tab->rows[irow]->numCorr; i++) 
       PyList_SetItem(tlist, i, PyString_InternFromString(Tab->rows[irow]->corrType[i]));
    PyDict_SetItemString(dict, "corrType",  tlist);

    tlist = PyList_New(2*Tab->rows[irow]->numCorr);
    for (i=0; i<Tab->rows[irow]->numCorr; i++) {
       PyList_SetItem(tlist, 2*i,   PyString_InternFromString(Tab->rows[irow]->corrProduct[2*i]));
       PyList_SetItem(tlist, 2*i+1, PyString_InternFromString(Tab->rows[irow]->corrProduct[2*i+1]));
    }
    PyDict_SetItemString(dict, "corrProduct",  tlist);

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetPolarization

// Convert Processor data to a python list of dicts
static PyObject* OASDMGetProcessor(ObitSDMData *asdm)
{
  ASDMProcessorTable* Tab =  ((ObitSDMData*)asdm)->ProcessorTab;
  PyObject *outList, *dict;
  olong irow, cnt;

  if (!Tab) return PyList_New(0);  // bad pointer?

  outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "processorId",      PyInt_FromLong((long)Tab->rows[irow]->processorId));
    PyDict_SetItemString(dict, "modeId",           PyInt_FromLong((long)Tab->rows[irow]->modeId));
    PyDict_SetItemString(dict, "processorType",    PyString_InternFromString(Tab->rows[irow]->processorType));
    PyDict_SetItemString(dict, "processorSubType", PyString_InternFromString(Tab->rows[irow]->processorSubType));

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetProcessor

// Convert SwitchCycle data to a python list of dicts
static PyObject* OASDMGetSwitchCycle(ObitSDMData *asdm)
{
  ASDMSwitchCycleTable* Tab =  ((ObitSDMData*)asdm)->SwitchCycleTab;
  PyObject *outList, *dict, *tlist;
  olong irow, i;

  if (!Tab) return PyList_New(0);  // bad pointer?

  outList = PyList_New(Tab->nrows);

  for (irow=0; irow<Tab->nrows; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "switchCycleId", PyInt_FromLong((long)Tab->rows[irow]->switchCycleId));
    PyDict_SetItemString(dict, "numStep",       PyInt_FromLong((long)Tab->rows[irow]->numStep));

    tlist = PyList_New(Tab->rows[irow]->numStep);
    for (i=0; i<Tab->rows[irow]->numStep; i++) 
       PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->weightArray[i]));
    PyDict_SetItemString(dict, "weightArray",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numStep);
    for (i=0; i<Tab->rows[irow]->numStep; i++) 
       PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->dirOffsetArray[i]));
    PyDict_SetItemString(dict, "dirOffsetArray",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numStep);
    for (i=0; i<Tab->rows[irow]->numStep; i++) 
       PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->freqOffsetArray[i]));
    PyDict_SetItemString(dict, "freqOffsetArray",  tlist);

    tlist = PyList_New(Tab->rows[irow]->numStep);
    for (i=0; i<Tab->rows[irow]->numStep; i++) 
       PyList_SetItem(tlist, i, PyFloat_FromDouble((double)Tab->rows[irow]->stepDurationArray[i]));
    PyDict_SetItemString(dict, "stepDurationArray",  tlist);

    PyList_SetItem(outList, irow, dict);
  }

  return outList;
} // end  GetSwitchCycle

// Get SpectralWindowArray from an ASDM and convert to a python dict
static PyObject* OASDMGetSpectralWindowArray(ObitSDMData *asdm, int mainRow, int SWOrder)
{
  gboolean lSWOrder;
  ASDMSpectralWindowArray* swa=NULL;
  PyObject *outDict, *dict, *tlist, *wlist;
  olong irow, i;

   if (SWOrder) lSWOrder = TRUE;
   else         lSWOrder = FALSE;

  swa = ObitSDMDataGetSWArray(asdm, (olong)mainRow, lSWOrder);

  outDict = PyDict_New();
  PyDict_SetItemString(outDict, "nwinds",  PyInt_FromLong((long)swa->nwinds));
  PyDict_SetItemString(outDict, "refFreq", PyFloat_FromDouble((double)swa->refFreq));
  PyDict_SetItemString(outDict, "refJD",   PyFloat_FromDouble((double)swa->refJD));

  tlist = PyList_New(swa->nwinds);
  for (i=0; i<swa->nwinds; i++) 
    PyList_SetItem(tlist, i, PyInt_FromLong((long)swa->order));
  PyDict_SetItemString(outDict, "order",  tlist);

  wlist = PyList_New(swa->nwinds);
  for (irow=0; irow<swa->nwinds; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "spectralWindowId", PyInt_FromLong((long)swa->winds[irow]->spectralWindowId));
    PyDict_SetItemString(dict, "selected",         PyBool_FromLong((long)swa->winds[irow]->selected));	
    PyDict_SetItemString(dict, "basebandNum",      PyInt_FromLong((long)swa->winds[irow]->basebandNum));
    PyDict_SetItemString(dict, "subbandNum",       PyInt_FromLong((long)swa->winds[irow]->subbandNum));
    PyDict_SetItemString(dict, "numChan",          PyInt_FromLong((long)swa->winds[irow]->numChan));
    PyDict_SetItemString(dict, "nCPoln",           PyInt_FromLong((long)swa->winds[irow]->nCPoln));
    PyDict_SetItemString(dict, "nAPoln" ,          PyInt_FromLong((long)swa->winds[irow]->nAPoln));
    PyDict_SetItemString(dict, "netSideband",      PyString_InternFromString(swa->winds[irow]->netSideband));
    PyDict_SetItemString(dict, "refChan",          PyFloat_FromDouble((double)swa->winds[irow]->refChan));
    PyDict_SetItemString(dict, "refFreq",          PyFloat_FromDouble((double)swa->winds[irow]->refFreq));
    PyDict_SetItemString(dict, "totBandwidth",     PyFloat_FromDouble((double)swa->winds[irow]->totBandwidth));
    PyDict_SetItemString(dict, "chanFreqStart",    PyFloat_FromDouble((double)swa->winds[irow]->chanFreqStart));
    PyDict_SetItemString(dict, "chanFreqStep",     PyFloat_FromDouble((double)swa->winds[irow]->chanFreqStep));
    PyDict_SetItemString(dict, "chanWidth",        PyFloat_FromDouble((double)swa->winds[irow]->chanWidth));
    PyDict_SetItemString(dict, "effectiveBw",      PyFloat_FromDouble((double)swa->winds[irow]->effectiveBw));
    PyDict_SetItemString(dict, "resolution",       PyFloat_FromDouble((double)swa->winds[irow]->resolution));
    PyList_SetItem(wlist, irow, dict);
  }

  PyDict_SetItemString(outDict, "winds",  wlist);

  // cleanup
  swa = ObitSDMDataKillSWArray(swa);

  return outDict;
} // end  GetSpectralWindowArray

// Get Antenna Array from an ASDM and convert to a python dict
static PyObject* OASDMGetAntArray(ObitSDMData *asdm, int mainRow)
{
  ASDMAntennaArray* aa=NULL;
  PyObject *outDict, *dict, *tlist, *alist;
  olong irow, i;
  gchar ctemp[3];

  aa = ObitSDMDataGetAntArray(asdm, (olong)mainRow);

  outDict = PyDict_New();
  PyDict_SetItemString(outDict, "nants",  PyInt_FromLong((long)aa->nants));
  PyDict_SetItemString(outDict, "maxAnt",  PyInt_FromLong((long)aa->maxAnt));
  PyDict_SetItemString(outDict, "arrayName", PyString_InternFromString(aa->arrayName));
  PyDict_SetItemString(outDict, "obsName", PyString_InternFromString(aa->obsName));
  PyDict_SetItemString(outDict, "refFreq", PyFloat_FromDouble((double)aa->refFreq));
  PyDict_SetItemString(outDict, "refJD",   PyFloat_FromDouble((double)aa->refJD));

  alist = PyList_New(aa->nants);
  for (irow=0; irow<aa->nants; irow++) {
    dict = PyDict_New();
    PyDict_SetItemString(dict, "antennaId",    PyInt_FromLong((long)aa->ants[irow]->antennaId));
    PyDict_SetItemString(dict, "antennaNo",    PyInt_FromLong((long)aa->ants[irow]->antennaNo));	
    PyDict_SetItemString(dict, "stationId",    PyInt_FromLong((long)aa->ants[irow]->stationId));
    PyDict_SetItemString(dict, "antName",      PyString_InternFromString(aa->ants[irow]->antName));
    PyDict_SetItemString(dict, "staName",      PyString_InternFromString(aa->ants[irow]->staName));
    PyDict_SetItemString(dict, "numPoln",      PyInt_FromLong((long)aa->ants[irow]->numPoln));
    PyDict_SetItemString(dict, "numPolnCorr",  PyInt_FromLong((long)aa->ants[irow]->numPolnCorr));
    PyDict_SetItemString(dict, "antennaMake",  PyInt_FromLong((long)aa->ants[irow]->antennaMake));
    PyDict_SetItemString(dict, "antennaType",  PyInt_FromLong((long)aa->ants[irow]->antennaType));
    PyDict_SetItemString(dict, "type",         PyInt_FromLong((long)aa->ants[irow]->type));
    PyDict_SetItemString(dict, "dishDiameter", PyFloat_FromDouble((double)aa->ants[irow]->dishDiameter));
    PyDict_SetItemString(dict, "time",         PyFloat_FromDouble((double)aa->ants[irow]->time));
    ctemp[0] = aa->ants[irow]->polnType[0]; ctemp[1] = aa->ants[irow]->polnType[1]; ctemp[2] = 0;
    PyDict_SetItemString(dict, "polnType",     PyString_InternFromString(ctemp));

    tlist = PyList_New(3);
    for (i=0; i<3; i++) 
      PyList_SetItem(tlist, i, PyFloat_FromDouble((double)aa->ants[irow]->antPosition[i]));
    PyDict_SetItemString(dict, "antPosition",  tlist);

    tlist = PyList_New(3);
    for (i=0; i<3; i++) 
      PyList_SetItem(tlist, i, PyFloat_FromDouble((double)aa->ants[irow]->offset[i]));
    PyDict_SetItemString(dict, "offset",  tlist);

     tlist = PyList_New(3);
    for (i=0; i<3; i++) 
      PyList_SetItem(tlist, i, PyFloat_FromDouble((double)aa->ants[irow]->staPosition[i]));
    PyDict_SetItemString(dict, "staPosition",  tlist);

   PyList_SetItem(alist, irow, dict);
  }

  PyDict_SetItemString(outDict, "ants",  alist);

  // cleanup
  aa = ObitSDMDataKillAntArray(aa);

  return outDict;
} // end  GetAntArray

// Get Source Array from an ASDM and convert to a python dict
static PyObject* OASDMGetSourceArray(ObitSDMData *asdm)
{
  ASDMSourceArray* sa=NULL;
  PyObject *outDict, *dict, *tlist, *slist;
  olong irow, i, lastSID=-1, cnt, ocnt;

  sa = ObitSDMDataGetSourceArray(asdm);

  outDict = PyDict_New();
  PyDict_SetItemString(outDict, "nsou",  PyInt_FromLong((long)sa->nsou));

  // count unique SID
  cnt = 0; lastSID=-1;
  for (irow=0; irow<sa->nsou; irow++) {
    if (sa->sou[irow]->sourceId!=lastSID) cnt++;
    lastSID = sa->sou[irow]->sourceId;
  }

  slist = PyList_New(cnt);
  lastSID=-1; ocnt = 0;
  for (irow=0; irow<sa->nsou; irow++) {
	// Already done this sourceID?
    if (sa->sou[irow]->sourceId==lastSID) continue;
    lastSID = sa->sou[irow]->sourceId; ocnt++;
    dict = PyDict_New();
    PyDict_SetItemString(dict, "sourceId",         PyInt_FromLong((long)sa->sou[irow]->sourceId));
    PyDict_SetItemString(dict, "sourceNo",         PyInt_FromLong((long)sa->sou[irow]->sourceNo));	
    PyDict_SetItemString(dict, "sourceName",       PyString_InternFromString(sa->sou[irow]->sourceName));
    PyDict_SetItemString(dict, "spectralWindowId", PyInt_FromLong((long)sa->sou[irow]->spectralWindowId));
    PyDict_SetItemString(dict, "code",             PyString_InternFromString(sa->sou[irow]->code));
    PyDict_SetItemString(dict, "numLines",         PyInt_FromLong((long)sa->sou[irow]->numLines));

    tlist = PyList_New(2);
    PyList_SetItem(tlist, 0, PyFloat_FromDouble((double)sa->sou[irow]->timeInterval[0]));
    PyList_SetItem(tlist, 1, PyFloat_FromDouble((double)sa->sou[irow]->timeInterval[1]));
    PyDict_SetItemString(dict, "timeInterval",  tlist);

    tlist = PyList_New(2);
    PyList_SetItem(tlist, 0, PyFloat_FromDouble((double)sa->sou[irow]->direction[0]));
    PyList_SetItem(tlist, 1, PyFloat_FromDouble((double)sa->sou[irow]->direction[1]));
    PyDict_SetItemString(dict, "direction",  tlist);

    tlist = PyList_New(2);
    PyList_SetItem(tlist, 0, PyFloat_FromDouble((double)sa->sou[irow]->properMotion[0]));
    PyList_SetItem(tlist, 1, PyFloat_FromDouble((double)sa->sou[irow]->properMotion[1]));
    PyDict_SetItemString(dict, "properMotion",  tlist);

    tlist = PyList_New(sa->sou[irow]->numLines);
    for (i=0; i<sa->sou[irow]->numLines; i++) 
      PyList_SetItem(tlist, i, PyFloat_FromDouble((double)sa->sou[irow]->restFrequency[i]));
    PyDict_SetItemString(dict, "restFrequency",  tlist);

    tlist = PyList_New(sa->sou[irow]->numLines);
    for (i=0; i<sa->sou[irow]->numLines; i++) 
      PyList_SetItem(tlist, i, PyFloat_FromDouble((double)sa->sou[irow]->sysVel[i]));
    PyDict_SetItemString(dict, "sysVel",  tlist);

    PyList_SetItem(slist, ocnt-1, dict);
  }

  //  correct number
  PyDict_SetItemString(outDict, "nsou",  PyInt_FromLong((long)cnt));

  PyDict_SetItemString(outDict, "sou",  slist);

  // cleanup
  sa = ObitSDMDataKillSourceArray(sa);

  return outDict;
} // end  GetSourceArray

extern int OASDMIsA (ObitSDMData* in) {
  return ObitSDMDataIsA(in);
} // end  OASDMIsA 

ObitSDMData* OASDMRef (ObitSDMData* in) {
  return ObitSDMDataRef (in);
} // end OASDMRef

ObitSDMData* OASDMUnref (ObitSDMData* in) {
  if (!ObitSDMDataIsA(in)) return NULL;
  return ObitSDMDataUnref (in);
} // end OASDMUnref

extern char* OASDMGetName (ObitSDMData* in) {
  return in->name;
} // end  OASDMGetName



%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitSDMData *me;
} OASDM;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitSDMData *me;
} OASDM;

%addmethods OASDM { 
  OASDM(char* name,  char* DataRoot, ObitErr *err) {
     OASDM *out;
     out = (OASDM *) malloc(sizeof(OASDM));
     if (strcmp(name, "None")) out->me = OASDMCreate(name, DataRoot, err);
     else out->me = NULL;
     return out;
   }
  ~OASDM() {
    self->me = OASDMUnref(self->me);
    free(self);
  }
};

