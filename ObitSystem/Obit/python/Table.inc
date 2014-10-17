/* $Id$   */  
/*--------------------------------------------------------------------*/
/* Swig module description for Table type                             */
/*                                                                    */
/*;  Copyright (C) 2004-2014                                          */
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
#include "ObitTable.h"
#include "ObitTableAN.h"
#include "ObitTableAT.h"
#include "ObitTableBL.h"
#include "ObitTableBP.h"
#include "ObitTableCC.h"
#include "ObitTableCL.h"
#include "ObitTableCQ.h"
#include "ObitTableCP.h"
#include "ObitTableCT.h"
#include "ObitTableFG.h"
#include "ObitTableFQ.h"
#include "ObitTableGC.h"
#include "ObitTableIM.h"
#include "ObitTableMC.h"
#include "ObitTableMF.h"
#include "ObitTableNI.h"
#include "ObitTableNX.h"
#include "ObitTableOB.h"
#include "ObitTableOF.h"
#include "ObitTablePC.h"
#include "ObitTablePD.h"
#include "ObitTablePS.h"
#include "ObitTableSN.h"
#include "ObitTableSU.h"
#include "ObitTableTY.h"
#include "ObitTableVL.h"
#include "ObitTableVZ.h"
#include "ObitTableWX.h"
#include "ObitTableUtil.h"
%}

%inline %{
// Routine to remove trailing blanks from a string 
static void TableDeeBlank (gchar *in) 
{
  olong i;

  for (i=strlen(in)-1; i>=0; i--) {
     if (in[i]==' ') in[i] = 0;
     else if (in[i]!=' ') break;
  }

} // end TableDeeBlank

extern ObitTable* TableCreate (gchar* name) {
  return newObitTable (name);
} // end  TableCreate

extern ObitTable* TableZap  (ObitTable *in, ObitErr *err) {
  return ObitTableZap (in, err);
} // end TableZap

extern ObitTable* TableCopy  (ObitTable *in, ObitTable *out, 
			      ObitErr *err) {
  return ObitTableCopy (in, out, err);
} // end  TableCopy

extern ObitTable* TableClone (ObitTable *in, ObitTable *out) {
   return  ObitTableClone (in, out);
} // end  TableClone

extern void TableConcat  (ObitTable *in, ObitTable *out, ObitErr *err) {
  ObitTableConcat (in, out, err);
} // end  TableCopy

// Open and close to fully instantiate
// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
extern int TablefullInstantiate (ObitTable* in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitTableOpen (in, laccess, err);
  ret = ObitTableClose (in, err);
  if ((err->error) || (ret!=OBIT_IO_OK)) return 1;
  else return 0;
} // end TablefullInstantiate

extern int TableOpen (ObitTable *in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;

  if (!strncmp (in->tabType,"AIPS AN", 7)) {
    ret =  ObitTableANOpen ((ObitTableAN*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS AT", 7)) {
    ret =  ObitTableATOpen ((ObitTableAT*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS BL", 7)) {
    ret =  ObitTableBLOpen ((ObitTableBL*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS BP", 7)) {
    ret =  ObitTableBPOpen ((ObitTableBP*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS CC", 7)) {
    ret =  ObitTableCCOpen ((ObitTableCC*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS CL", 7)) {
    ret =  ObitTableCLOpen ((ObitTableCL*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS CP", 7)) {
    ret =  ObitTableCPOpen ((ObitTableCP*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS CQ", 7)) {
    ret =  ObitTableCQOpen ((ObitTableCQ*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS CT", 7)) {
    ret =  ObitTableCTOpen ((ObitTableCT*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS FG", 7)) {
    ret =  ObitTableFGOpen ((ObitTableFG*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS FQ", 7)) {
    ret =  ObitTableFQOpen ((ObitTableFQ*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS FS", 7)) {
    ret =  ObitTableFSOpen ((ObitTableFS*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS GC", 7)) {
    ret =  ObitTableGCOpen ((ObitTableGC*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS IM", 7)) {
    ret =  ObitTableIMOpen ((ObitTableIM*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS MC", 7)) {
    ret =  ObitTableMCOpen ((ObitTableMC*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS MF", 7)) {
    ret =  ObitTableMFOpen ((ObitTableMF*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS NI", 7)) {
    ret =  ObitTableNIOpen ((ObitTableNI*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS NX", 7)) {
    ret =  ObitTableNXOpen ((ObitTableNX*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS OB", 7)) {
    ret =  ObitTableOBOpen ((ObitTableOB*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS OF", 7)) {
    ret =  ObitTableOFOpen ((ObitTableOF*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS PC", 7)) {
    ret =  ObitTablePCOpen ((ObitTablePC*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS PD", 7)) {
    ret =  ObitTablePDOpen ((ObitTablePD*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS PS", 7)) {
    ret =  ObitTablePSOpen ((ObitTablePS*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS SN", 7)) {
    ret =  ObitTableSNOpen ((ObitTableSN*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS SU", 7)) {
    ret =  ObitTableSUOpen ((ObitTableSU*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS TY", 7)) {
    ret =  ObitTableTYOpen ((ObitTableTY*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS VL", 7)) {
    ret =  ObitTableVLOpen ((ObitTableVL*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS VZ", 7)) {
    ret =  ObitTableVZOpen ((ObitTableVZ*)in, laccess, err);
  }
  else if (!strncmp (in->tabType,"AIPS WX", 7)) {
    ret =  ObitTableWXOpen ((ObitTableWX*)in, laccess, err);
  }
  else
    ret = ObitTableOpen (in, laccess, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Open

// force header update 
extern void TableDirty (ObitTable *in) {
  in->myStatus = OBIT_Modified;
} // end Dirty

extern int TableClose (ObitTable *in, ObitErr *err) {
  ObitIOCode ret;

  if (!strncmp (in->tabType,"AIPS AN", 7)) {
    ret =  ObitTableANClose ((ObitTableAN*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS AT", 7)) {
    ret =  ObitTableATClose ((ObitTableAT*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS BL", 7)) {
    ret =  ObitTableBLClose ((ObitTableBL*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS BP", 7)) {
    ret =  ObitTableBPClose ((ObitTableBP*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS CC", 7)) {
    ret =  ObitTableCCClose ((ObitTableCC*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS CL", 7)) {
    ret =  ObitTableCLClose ((ObitTableCL*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS CP", 7)) {
    ret =  ObitTableCPClose ((ObitTableCP*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS CQ", 7)) {
    ret =  ObitTableCQClose ((ObitTableCQ*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS CT", 7)) {
    ret =  ObitTableCTClose ((ObitTableCT*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS FG", 7)) {
    ret =  ObitTableFGClose ((ObitTableFG*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS FQ", 7)) {
    ret =  ObitTableFQClose ((ObitTableFQ*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS FS", 7)) {
    ret =  ObitTableFSClose ((ObitTableFS*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS GC", 7)) {
    ret =  ObitTableGCClose ((ObitTableGC*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS IM", 7)) {
    ret =  ObitTableIMClose ((ObitTableIM*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS MC", 7)) {
    ret =  ObitTableMCClose ((ObitTableMC*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS MF", 7)) {
    ret =  ObitTableMFClose ((ObitTableMF*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS NI", 7)) {
    ret =  ObitTableNIClose ((ObitTableNI*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS NX", 7)) {
    ret =  ObitTableNXClose ((ObitTableNX*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS OB", 7)) {
    ret =  ObitTableOBClose ((ObitTableOB*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS OF", 7)) {
    ret =  ObitTableOFClose ((ObitTableOF*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS PC", 7)) {
    ret =  ObitTablePCClose ((ObitTablePC*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS PD", 7)) {
    ret =  ObitTablePDClose ((ObitTablePD*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS PS", 7)) {
    ret =  ObitTablePSClose ((ObitTablePS*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS SN", 7)) {
    ret =  ObitTableSNClose ((ObitTableSN*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS SU", 7)) {
    ret =  ObitTableSUClose ((ObitTableSU*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS TY", 7)) {
    ret =  ObitTableTYClose ((ObitTableTY*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS VL", 7)) {
    ret =  ObitTableVLClose ((ObitTableVL*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS VZ", 7)) {
    ret =  ObitTableVZClose ((ObitTableVZ*)in, err);
  }
  else if (!strncmp (in->tabType,"AIPS WX", 7)) {
    ret =  ObitTableWXClose ((ObitTableWX*)in, err);
  }
  else
    ret =  ObitTableClose (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Close

extern  PyObject *TableReadRow (ObitTable *in, int rowno,
                         ObitErr *err) {
  ObitIOCode ret = OBIT_IO_SpecErr;
  ObitTableDesc *desc = in->myDesc;
  olong i, j, j1, j2, k, lrowno = rowno;
  ObitTableRow* row = newObitTableRow (in);
  PyObject *outDict = PyDict_New();
  PyObject *list;
  PyObject *o;
  gshort   *idata;
  olong     *jdata;
  oint     *kdata;
  olong    *ldata;
  gchar    *cdata, *ctemp;
  gboolean *bdata;
  ofloat   *fdata;
  odouble  *ddata;
  gchar *routine = "TableReadRow";

  if (err->error) return outDict;

  ret = ObitTableReadRow (in, lrowno, row, err);
  if ((ret==OBIT_IO_OK) && (!err->error)) {
  // Convert row to python dict
  /* Table name */
  o = PyString_InternFromString(desc->TableName);
  PyDict_SetItemString(outDict, "Table name", o);
  Py_DECREF(o);
  /* number of fields  */
  o = PyInt_FromLong((long)desc->nfield);
  PyDict_SetItemString(outDict, "NumFields", o);
  Py_DECREF(o);
  /* Loop over fields */
  for (i=0; i<in->myDesc->nfield; i++) {
    if (desc->type[i] == OBIT_string)
      list = PyList_New(desc->repeat[i]/MAX (1,desc->dim[i][0]));
    else
      list = PyList_New(desc->repeat[i]);
    /* Fill list by type */
  switch (desc->type[i]) { 
    case OBIT_short:
      idata = ((gshort*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        PyList_SetItem(list, j, PyInt_FromLong((long)idata[j]));
      break;
    case OBIT_int:
      jdata = ((olong*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        PyList_SetItem(list, j, PyInt_FromLong((long)jdata[j]));
      break;
    case OBIT_oint:
      kdata = ((oint*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        PyList_SetItem(list, j, PyInt_FromLong((long)kdata[j]));
      break;
    case OBIT_long:
      ldata = ((olong*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        PyList_SetItem(list, j, PyInt_FromLong((long)ldata[j]));
      break;
    case OBIT_float:
      fdata = ((ofloat*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        PyList_SetItem(list, j, PyFloat_FromDouble((double)fdata[j]));
      break;
    case OBIT_double:
      ddata = ((odouble*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        PyList_SetItem(list, j, PyFloat_FromDouble((double)ddata[j]));
      break;
    case OBIT_string:
      cdata = ((gchar*)row->myRowData)+desc->offset[i];
       /* null terminate string */
      ctemp = g_malloc0(desc->dim[i][0]+1);
      for (j=0; j<desc->repeat[i]/MAX (1,desc->dim[i][0]); j++) {
        for (k=0; k<desc->dim[i][0]; k++) ctemp[k] = cdata[k];  ctemp[k] = 0;
        PyList_SetItem(list, j, PyString_InternFromString(ctemp));
	cdata += desc->dim[i][0];
      }
      g_free (ctemp);
      break;
    case OBIT_bool:
      bdata = ((gboolean*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        PyList_SetItem(list, j, PyInt_FromLong((long)bdata[j]));
      break;
    case OBIT_bits:
      kdata = ((oint*)row->myRowData)+desc->offset[i];
      j1 = 0; j2 = 0;
      for (j=0; j<desc->repeat[i]; j++) {
	if (j2>=sizeof(oint)) {  /* next integer worth? */
          j1++;
          j2 = 0;
        }
	o = PyInt_FromLong((long)((kdata[j1] & (1 << j2)) ? 1 : 0));
        j2++;
        PyList_SetItem(list, j, o);
      }
      break;
    default:
      /* Cannot deal with type */
      Obit_log_error(err, OBIT_Error, 
		   "%s: Cannot deal with Table row data type %d in %s", 
		   routine, desc->type[i], in->name);
    }; /* end switch by type */

    /* Set field in output dict */
    ctemp =  g_strdup(desc->FieldName[i]);
    TableDeeBlank(ctemp);
    PyDict_SetItemString(outDict, ctemp, list);
    Py_DECREF(list);
    g_free (ctemp);
  } /* End loop over fields */

 } /* end of read OK */

  row = ObitTableRowUnref(row);
  return outDict;	
} // end TableReadRow

extern int TableWriteRow (ObitTable *in, int rowno, PyObject *inDict,
                          ObitErr *err) {
  ObitIOCode ret = OBIT_IO_SpecErr;
  ObitTableDesc *desc = in->myDesc;
  olong i, j, j1, j2, k, lrowno = rowno;
  ObitTableRow* row = newObitTableRow (in);
  PyObject *list = NULL, *TabName=NULL;
  gshort   *idata;
  olong     *jdata;
  oint     *kdata;
  olong    *ldata, ltemp;
  gchar    *cdata, *ctemp, *tstr=NULL;
  gboolean *bdata, bad;
  ofloat   *fdata;
  odouble  *ddata;
  gchar *routine = "TableWriteRow";

  if (err->error) return 1;

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return 1;
  }

  /* Check that correct table type */
  TabName = PyDict_GetItemString(inDict, "Table name");
  if (TabName!=NULL)
    tstr = PyString_AsString(TabName);
  if (TabName==NULL || tstr==NULL) {  /* Found it? */
    Obit_log_error(err, OBIT_Error, "%s: Table Name not given", routine);
    return 1;
  }
  if (strncmp (tstr, desc->TableName, strlen(tstr))) {
    Obit_log_error(err, OBIT_Error, "%s: Table type '%s' NOT '%s'",
                   routine, tstr, desc->TableName);
    return 1;
  }
  /* Check number of fields */
  ltemp = PyInt_AsLong(PyDict_GetItemString(inDict, "NumFields"));
  if (ltemp!=desc->nfield) {
    Obit_log_error(err, OBIT_Error, "%s: no. columns %d NOT  %d",
                   routine, ltemp, desc->nfield);
    return 1;
  }

  /* attach row to  output buffer */
  ObitTableSetRow (in, row, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, 1);

  // Convert python dict to row
  /* Loop over fields */
  for (i=0; i<in->myDesc->nfield; i++) {
    bad = FALSE;  /* Check for field size */
    ctemp =  g_strdup(desc->FieldName[i]);
    TableDeeBlank(ctemp);
    list = PyDict_GetItemString(inDict, ctemp);
    if (list==NULL) {
      Obit_log_error(err, OBIT_Error, "%s: %s not given", routine, ctemp);
      return 1;
    }
    /* Is this really a list */
    if (!PyList_Check(list)) {
      Obit_log_error(err, OBIT_Error, "%s: %s member not a PyList",
                     routine, ctemp);
      return 1;
    }
    g_free (ctemp);
    /* Fill list by type */
  switch (desc->type[i]) { 
    case OBIT_short:
      bad = desc->repeat[i] != PyList_Size(list);
      idata = ((gshort*)row->myRowData)+desc->offset[i];
      if (!bad) {
        for (j=0; j<desc->repeat[i]; j++) 
          idata[j] = (gshort)PyInt_AsLong(PyList_GetItem(list, j));
      }
      break;
    case OBIT_int:
      bad = desc->repeat[i] != PyList_Size(list);
      jdata = ((olong*)row->myRowData)+desc->offset[i];
      if (!bad) {
        for (j=0; j<desc->repeat[i]; j++) 
          jdata[j] = (olong)PyInt_AsLong(PyList_GetItem(list, j));
      }
      break;
    case OBIT_oint:
      bad = desc->repeat[i] != PyList_Size(list);
      kdata = ((oint*)row->myRowData)+desc->offset[i];
      if (!bad) {
        for (j=0; j<desc->repeat[i]; j++) 
          kdata[j] = (oint)PyInt_AsLong(PyList_GetItem(list, j));
      }
      break;
    case OBIT_long:
      bad = desc->repeat[i] != PyList_Size(list);
      ldata = ((olong*)row->myRowData)+desc->offset[i];
      if (!bad) {
        for (j=0; j<desc->repeat[i]; j++) 
          ldata[j] = (olong)PyInt_AsLong(PyList_GetItem(list, j));
      }
      break;
    case OBIT_float:
      bad = desc->repeat[i] != PyList_Size(list);
      fdata = ((ofloat*)row->myRowData)+desc->offset[i];
      if (!bad) {
        for (j=0; j<desc->repeat[i]; j++) 
          fdata[j] =  (ofloat)PyFloat_AsDouble(PyList_GetItem(list, j));
      }
      break;
    case OBIT_double:
      bad = desc->repeat[i] != PyList_Size(list);
      ddata = ((odouble*)row->myRowData)+desc->offset[i];
      if (!bad) {
        for (j=0; j<desc->repeat[i]; j++) 
          ddata[j] =  (odouble)PyFloat_AsDouble(PyList_GetItem(list, j));
      }
      break;
    case OBIT_string:
      bad = (desc->repeat[i]/MAX (1,desc->dim[i][0])) != PyList_Size(list);
      if (bad) break;
      cdata = ((gchar*)row->myRowData)+desc->offset[i];
      /* null terminate string */
      for (j=0; j<desc->repeat[i]/MAX (1,desc->dim[i][0]); j++) {
        ctemp = PyString_AsString(PyList_GetItem(list, j));
        for (k=0; k<strlen(ctemp); k++) cdata[k] = ctemp[k];
	cdata += desc->dim[i][0];
      }
      break;
    case OBIT_bool:
      bad = desc->repeat[i] != PyList_Size(list);
      if (bad) break;
      bdata = ((gboolean*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        bdata[j] = (gboolean)PyInt_AsLong(PyList_GetItem(list, j));
      break;
    case OBIT_bits:
      bad = desc->repeat[i] != PyList_Size(list);
      if (bad) break;
      kdata = ((oint*)row->myRowData)+desc->offset[i];
      kdata[0] = 0;
      j1 = 0; j2 = 0;
      for (j=0; j<desc->repeat[i]; j++) {
	if (j2>=sizeof(oint)) {  /* next integer worth? */
          j1++;
          j2 = 0;
        }
        kdata[j1] |= (PyInt_AsLong(PyList_GetItem(list, j2)) ? (1 << j) : 0);
        j2++;
      }
      break;
    default:
      /* Cannot deal with type */
      Obit_log_error(err, OBIT_Error, 
		   "%s: Cannot deal with Table row data type %d in %s", 
		   routine, desc->type[i], in->name);
      return 1;
  }; /* end switch by type */

  /* Check if sizes compatible */
  if (bad) {
      Obit_log_error(err, OBIT_Error, 
		   "%s: wrong size %d %d for %s", 
		   routine, desc->repeat[i], (olong)PyList_Size(list), desc->FieldName[i]);
      return 1;
  }

  } /* End loop over fields */

  ret =  ObitTableWriteRow (in, lrowno, row, err);
  row = ObitTableRowUnref(row);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end TableWriteRow

extern ObitTable* TableUnref (ObitTable* in) {
  if (!ObitTableIsA(in)) return NULL;
  return ObitTableUnref(in);
}

extern ObitTable*  TableRef (ObitTable* in) {
  return ObitTableRef(in);
}

extern ObitInfoList* TableGetList (ObitTable* in) {
  return ObitInfoListRef(in->info);
}

extern ObitInfoList* TableGetIOList (ObitTable* in) {
  ObitInfoList *info=NULL;
  if (in->myIO!=NULL) info = ((ObitTableDesc*)(in->myIO->myDesc))->info;
  return ObitInfoListRef(info);
}

extern ObitTableDesc* TableGetDesc (ObitTable* in) {
  return ObitTableDescRef(in->myDesc);
}

extern ObitTableDesc* TableGetIODesc (ObitTable* in) {
  ObitTableDesc *desc=NULL;
  if (in->myIO!=NULL) desc = (ObitTableDesc*)(in->myIO->myDesc);
  return ObitTableDescRef(desc);
}

extern void TableSetDesc (ObitTable* in, ObitTableDesc* desc) {
  in->myDesc = ObitTableDescUnref(in->myDesc);
  in->myDesc = ObitTableDescRef(desc);
}

extern long TableGetVer (ObitTable* in) {
  return in->myDesc->version;
}

extern int TableIsA (ObitTable* in) {
  return ObitTableIsA(in);
}

extern char* TableGetName (ObitTable* in) {
  if (ObitTableIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

// Table utilities 
extern int  TableUtilSort (ObitTable* in, char *colName, int desc, ObitErr *err) {
  ObitIOCode ret;
  gboolean ldesc;
  ldesc = desc != 0;
  ret =  ObitTableUtilSort (in, colName, ldesc, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitTable *me;
} Table;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitTable *me;
} Table;

%addmethods Table { 
  Table(char *name) {
     Table *out;
     out = (Table *) malloc(sizeof(Table));
     if (strcmp(name, "None")) out->me = TableCreate(name);
     else out->me = NULL;
     return out;
   }
  ~Table() {
    if (!self) return;
    self->me = TableUnref(self->me);
    free(self);
  }
};

