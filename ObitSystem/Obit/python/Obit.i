%module Obit
/* $Id: ObitTypeMaps.swig,v 1.6 2008/05/06 13:22:19 bcotton Exp $ */  
/*--------------------------------------------------------------------*/
/* Swig typemaps for Obit types                                       */
/*                                                                    */
/*   Copyright (C) 2004,2008                                          */
/*   Associated Universities, Inc. Washington DC, USA.                */
/*                                                                    */
/*   This program is free software; you can redistribute it and/or    */
/*   modify it under the terms of the GNU General Public License as   */
/*   published by the Free Software Foundation; either version 2 of   */
/*   the License, or (at your option) any later version.              */
/*                                                                    */
/*   This program is distributed in the hope that it will be useful,  */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of   */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    */
/*   GNU General Public License for more details.                     */
/*                                                                    */
/*   You should have received a copy of the GNU General Public        */
/*   License along with this program; if not, write to the Free       */
/*   Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,     */
/*   MA 02139, USA.                                                   */
/*                                                                    */
/*   Correspondence this software should be addressed as follows:     */
/*          Internet email: bcotton@nrao.edu.                         */
/*          Postal address: William Cotton                            */
/*                          National Radio Astronomy Observatory      */
/*                          520 Edgemont Road                         */
/*                          Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

// Convert list into a long array
%typemap (python,in) long * {
  if (PyList_Check($source)) {
    int size = PyList_Size($source);
    int i = 0;
    $target = (long*) malloc((size+1)*sizeof(long));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyInt_Check(o)) {
         $target[i] = PyInt_AsLong(o);
      } else {
         PyErr_SetString(PyExc_TypeError,"list must contain longs");
         free($target);
         return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// Convert list into an int array
%typemap (python,in) int* {
  if (PyList_Check($source)) {
    int size = PyList_Size($source);
    int i = 0;
    $target = (int*) malloc((size+1)*sizeof(int));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyInt_Check(o)) {
         $target[i] = (int)((PyIntObject*)o)->ob_ival;
      } else {
         PyErr_SetString(PyExc_TypeError,"list must contain ints");
         free($target);
         return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// Convert String into an char array
%typemap (python,in) char* {
  if (PyString_Check($source)) {
    int size = PyString_Size($source);
    char *str;
    int i = 0;
    $target = (char*) malloc((size+1));
    str = PyString_AsString($source);
    for (i = 0; i < size; i++) {
      $target[i] = str[i];
    }
    $target[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError,"not a string");
    return NULL;
  }
}

// Convert list into an gboolean array
%typemap (python,in) gboolean* {
  if (PyList_Check($source)) {
    int size = PyList_Size($source);
    int i = 0;
    $target = (int*) malloc((size+1)*sizeof(gboolean));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyBool_Check(o)) {
         $target[i] = (gboolean)((PyBoolObject*)o)->ob_ival;
      } else if (PyInt_Check(o)) {  /* also allow int for backward compatability */
         $target[i] = (gboolean)((PyIntObject*)o)->ob_ival;
      } else {
         PyErr_SetString(PyExc_TypeError,"list must contain ints");
         free($target);
         return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// Convert list into a float array
%typemap (python,in) float* {
  if (PyList_Check($source)) {
    int size = PyList_Size($source);
    int i = 0;
    $target = (float*) malloc((size+1)*sizeof(float));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyFloat_Check(o))
         $target[i] = (float)((PyFloatObject*)o)->ob_fval;
      else {
         PyErr_SetString(PyExc_TypeError,"list must contain floats");
         free($target);
         return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// Convert list into a double array
%typemap (python,in) double * {
  if (PyList_Check($source)) {
    int size = PyList_Size($source);
    int i = 0;
    $target = (double*) malloc((size+1)*sizeof(double));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyFloat_Check(o))
         $target[i] = (double)((PyFloatObject*)o)->ob_fval;
      else {
         PyErr_SetString(PyExc_TypeError,"list must contain doubles");
         free($target);
         return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// Convert list into a string array
%typemap (python,in) char ** {
  if (PyList_Check($source)) {
    int size2, size = PyList_Size($source);
    int j, i = 0;
    char *tstr;

    $target = (char**) malloc((size+1)*sizeof(char*));
    $target[size] = NULL;  // last string NULL
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyString_Check(o)) {
	 size2 = PyString_Size(o);
         $target[i] = (char*) malloc(size2+1);
	 tstr = PyString_AsString(o);
         for (j=0; j<=size2; j++) $target[i][j] = tstr[j];
      } else {
         PyErr_SetString(PyExc_TypeError,"list must contain Strings");
         free($target);
         return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// Convert list into a ObitImage* array
%typemap (python,in) ObitImage ** {
  if (PyList_Check($source)) {
    int size2, size = PyList_Size($source);
    int j, i = 0;
    char *tstr;

    $target = (ObitImage**) malloc((size+1)*sizeof(ObitImage*));
    $target[size] = NULL;  // last pointer NULL
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyString_Check(o)) {
        if (SWIG_GetPtrObj(o,(void **) &$target[i],"_ObitImage_p")) {
           PyErr_SetString(PyExc_TypeError,"Type error in argument. Expected _ObitImage_p.");
           return NULL;
         }
         if (!ObitImageIsA((ObitImage*)$target[i])) {  // check */
           PyErr_SetString(PyExc_TypeError,"Type error. Expected ObitImage Object.");
           return;
         }
      } else {
         PyErr_SetString(PyExc_TypeError,"list must contain Strings (ObitImage pointers)");
         free($target);
         return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// Convert PyObject into a PyDictObject or PyList
%typemap (python,in) PyObject* {
  if (PyList_Check($source)) {
    $target = PyDict_Copy(PyList_GetItem($source,0));
  } else if (PyDict_Check($source)) {
    $target = PyDict_Copy($source);
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list or dict");
    return NULL;
  }
}

// This cleans up the long * array we malloc'd before the function call
%typemap(python,freearg) long * {
  free((long *) $source);
}
// This cleans up the int * array we malloc'd before the function call
%typemap(python,freearg) int * {
  free((int *) $source);
}
// This cleans up the char * string we malloc'd before the function call
%typemap(python,freearg) char * {
  free((char *) $source);
}
// This cleans up the float * array we malloc'd before the function call
%typemap(python,freearg) gboolean * {
  free((gboolean *) $source);
}
// This cleans up the float * array we malloc'd before the function call
%typemap(python,freearg) float * {
  free((float *) $source);
}
// This cleans up the double * array we malloc'd before the function call
%typemap(python,freearg) double * {
  free((double *) $source);
}
// This cleans up the char ** array we malloc'd before the function call
%typemap(python,freearg) char ** {
  int i = 0;
  
  while ($source[i]!=NULL) { // last string should be NULL
    free ((char *) $source[i++]);
  }
  free((char **) $source);
}
// This cleans up the ObitImage ** array we malloc'd before the function call
%typemap(python,freearg) ObitImage ** {
  free((char **) $source);
}
// This cleans up the PyObject we malloc'd before the function call
%typemap(python,freearg) PyObject * {
  Py_XDECREF ($source);
}

// This allows a C function to return a long* as a Python list with 7 elements
%typemap(python,out) long* {
  int len,i;
  len = 7;
  $target = PyList_New(len);
  for (i = 0; i < len; i++) {
    PyList_SetItem($target,i,PyInt_FromLong($source[i]));
  }
}

// This allows a C function to return a int* as a Python list with 7 elements
%typemap(python,out) int* {
  int len,i;
  len = 7;
  $target = PyList_New(len);
  for (i = 0; i < len; i++) {
    PyList_SetItem($target,i,PyInt_FromLong((long)$source[i]));
  }
}

// This allows a C function to return a float* as a Python list with 7 elements
%typemap(python,out) float* {
  int len,i;
  len = 7;
  $target = PyList_New(len);
  for (i = 0; i < len; i++) {
    PyList_SetItem($target,i,PyFloat_FromDouble((double)$source[i]));
  }
}

// This allows a C function to return a double* as a Python list with 7 elements
%typemap(python,out) double* {
  int len,i;
  len = 7;
  $target = PyList_New(len);
  for (i = 0; i < len; i++) {
    PyList_SetItem($target,i,PyFloat_FromDouble($source[i]));
  }
}

// This allows a C function to return a char** as a Python list with 7 elements
%typemap(python,out) char** {
  int len,i;
  len = 7;
  $target = PyList_New(len);
  for (i = 0; i < len; i++) {
    if (!$source[i]) break;
    PyList_SetItem($target,i,PyString_FromString($source[i]));
  }
}

// This allows a C function to return a PyObject* as a Python dictionary or list or string
%typemap(python,out) PyObject*{
  if (PyList_Check($source) || PyDict_Check($source)
      || PyString_Check($source) || PyBuffer_Check($source)) {
    $target = $source;
  } else {
    PyErr_SetString(PyExc_TypeError,"output PyObject not dict or list");
    return NULL;
  }
}

// This tells SWIG to treat a long * argument with name 'outValue2' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.
// Assume there are always two output values

%typemap(python,argout) long *outValue2 {
  PyObject *o;
  if ((!$target) || ($target == Py_None)) {
    $target = PyList_New(0);
    o = PyInt_FromLong($source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
    o = PyInt_FromLong($source[1]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    o = PyInt_FromLong($source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
    o = PyInt_FromLong($source[1]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   }
}

// This tells SWIG to treat a long * argument with name 'outValue1' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.
// Assume there are always two output values

%typemap(python,argout) long *outValue1 {
  PyObject *o;
  if ((!$target) || ($target == Py_None)) {
    $target = PyList_New(0);
    o = PyInt_FromLong($source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
    o = PyInt_FromLong($source[1]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    o = PyInt_FromLong($source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
    o = PyInt_FromLong($source[1]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   }
}

// This tells SWIG to treat a double * argument with name 'outDbl1' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.

%typemap(python,argout) double *outDbl1 {
  PyObject *o;
  if ((!$target) || ($target == Py_None)) {
    $target = PyList_New(0);
    o = PyFloat_FromDouble($source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    o = PyFloat_FromDouble($source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   }
}

// This tells SWIG to treat a double * argument with name 'outDbl2' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.

%typemap(python,argout) double *outDbl2 {
  PyObject *o;
  if ((!$target) || ($target == Py_None)) {
    $target = PyList_New(0);
    o = PyFloat_FromDouble($source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    o = PyFloat_FromDouble($source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   }
}

// This tells SWIG to treat a double * argument with name 'outFlt1' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.

%typemap(python,argout) float *outFlt1 {
  PyObject *o;
  if ((!$target) || ($target == Py_None)) {
    $target = PyList_New(0);
    o = PyFloat_FromDouble((double)$source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    o = PyFloat_FromDouble((double)$source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   }
}

// This tells SWIG to treat a double * argument with name 'outFlt2' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.

%typemap(python,argout) float *outFlt2 {
  PyObject *o;
  if ((!$target) || ($target == Py_None)) {
    $target = PyList_New(0);
    o = PyFloat_FromDouble((double)$source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    o = PyFloat_FromDouble((double)$source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   }
}

// This tells SWIG to treat an int * argument with name 'outInt1' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.

%typemap(python,argout) int *outInt1 {
  PyObject *o;
  if ((!$target) || ($target == Py_None)) {
    $target = PyList_New(0);
    o = PyInt_FromLong((long)$source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    o = PyInt_FromLong((long)$source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   }
}

// This tells SWIG to treat an int * argument with name 'outInt2' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.

%typemap(python,argout) int *outInt2 {
  PyObject *o;
  if ((!$target) || ($target == Py_None)) {
    $target = PyList_New(0);
    o = PyInt_FromLong((long)$source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    o = PyInt_FromLong((long)$source[0]);
    PyList_Append($target,o);
    Py_XDECREF(o);
   }
}

/* $Id: AIPSDir.inc,v 1.11 2007/09/11 12:38:38 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for AIPS directory utilities               */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitAIPS.h"
#include "ObitAIPSDir.h"
%}

// This cleans up the char * string we malloc'd in AIPSDirInfo
%typemap(python,ret) char *AIPSDirInfo {
  if ($source) g_free($source);
}

%inline %{
extern int AIPSDirFindCNO(int disk, int user,  char *Aname, char *Aclass, 
    		          char *Atype, int seq, ObitErr *err)
{
  gchar LAname[13], LAclass[7], LAtype[3];
  olong i, l;

  /* Init AIPS fixed strings - blank fill until end */
  for (i=0; i<12; i++) LAname[i]  = ' ';LAname[i]  = 0;
  for (i=0; i<6; i++)  LAclass[i] = ' ';LAclass[i] = 0;
  for (i=0; i<2; i++)  LAtype[i]  = ' ';LAtype[i]  = 0;

  /* Copy string input into AIPS fixed strings */
  l = MIN (12, strlen(Aname));
  for (i=0; i<l; i++) LAname[i]  = Aname[i];
  l = MIN (6, strlen(Aclass));
  for (i=0; i<l; i++)  LAclass[i] = Aclass[i];
  l = MIN (2, strlen(Atype));
  for (i=0; i<l; i++)  LAtype[i] = Atype[i];

  return ObitAIPSDirFindCNO(disk, user, LAname, LAclass, LAtype, seq, err);
} /* end AIPSDirFindCNO */

extern int AIPSDirHiSeq(int disk, int user,  char *Aname, char *Aclass, 
    		        char *Atype, ObitErr *err)
{
  gchar LAname[13], LAclass[7], LAtype[3];
  olong i, l;

  /* Init AIPS fixed strings - blank fill until end */
  for (i=0; i<12; i++) LAname[i]  = ' ';LAname[i]  = 0;
  for (i=0; i<6; i++)  LAclass[i] = ' ';LAclass[i] = 0;
  for (i=0; i<2; i++)  LAtype[i]  = ' ';LAtype[i]  = 0;

  /* Copy string input into AIPS fixed strings */
  l = MIN (12, strlen(Aname));
  for (i=0; i<l; i++) LAname[i]  = Aname[i];
  l = MIN (6, strlen(Aclass));
  for (i=0; i<l; i++)  LAclass[i] = Aclass[i];
  l = MIN (2, strlen(Atype));
  for (i=0; i<l; i++)  LAtype[i] = Atype[i];

  return ObitAIPSDirHiSeq(disk, user, LAname, LAclass, LAtype, TRUE, err);
} /* end AIPSDirHiSeq */


extern int AIPSDirAlloc(int disk, int user,  char *Aname, char *Aclass, 
		        char *Atype, int seq, ObitErr *err)
{
  gchar LAname[13], LAclass[7], LAtype[3];
  gboolean exist;
  olong i, l;

  /* Init AIPS fixed strings - blank fill until end */
  for (i=0; i<12; i++) LAname[i]  = ' ';LAname[i]  = 0;
  for (i=0; i<6; i++)  LAclass[i] = ' ';LAclass[i] = 0;
  for (i=0; i<2; i++)  LAtype[i]  = ' ';LAtype[i]  = 0;

  /* Copy string input into AIPS fixed strings */
  l = MIN (12, strlen(Aname));
  for (i=0; i<l; i++) LAname[i]  = Aname[i];
  l = MIN (6, strlen(Aclass));
  for (i=0; i<l; i++)  LAclass[i] = Aclass[i];
  l = MIN (2, strlen(Atype));
  for (i=0; i<l; i++)  LAtype[i] = Atype[i];

  return ObitAIPSDirAlloc(disk, user, LAname, LAclass, LAtype, seq, &exist, err);
} /* end AIPSDirAlloc */

extern void AIPSDirRemoveEntry(int disk, int user, int cno, ObitErr *err)
{
  ObitAIPSDirRemoveEntry(disk, user, cno, err);
} /* end  AIPSDirRemoveEntry */

extern int AIPSDirNumber(int disk, int user, ObitErr *err)
{
  return (int)ObitAIPSDirNumber((olong)disk, (olong)user, err);
} /* end AIPSDirNumber */

/* Returns NULL if no entry */
extern char* AIPSDirInfo(int disk, int user, int cno, ObitErr *err)
{
  ObitAIPSDirCatEntry *entry;
  gchar *out=NULL, stat[5];
  olong nout=60;

  entry = ObitAIPSDirGetEntry ((olong)disk, (olong)user, (olong)cno, err);
  if ((err->error) || (entry==NULL)) return out;
  if (entry->user!=user) {
      g_free (entry);
      return out;
  }

  if (entry->status==0) strncpy (stat, "    ", 4);
  else if (entry->status<0) strncpy (stat, "WRIT", 4);
  else if (entry->status>0) strncpy (stat, "READ", 4);

  out = g_malloc(nout);
  g_snprintf (out, nout, "%-12.12s.%-6.6s. %4d %-2.2s          ",
              entry->name, entry->class, entry->seq, entry->type);
  ObitAIPSDirGetAccess (entry, &out[29]);
  g_snprintf (&out[strlen(out)], nout-strlen(out), " %-4.4s", stat);
  g_free (entry);
  return (char*)out;
} /* end AIPSDirInfo */

extern int AIPSDirStatus(int disk, int user, int cno, int code, ObitErr *err)
{
  return ObitAIPSDirStatus (disk, user, cno, code, err);
} /* end AIPSDirStatus */

extern int AIPSSetDirname(int disk, char *dir, ObitErr *err)
{
   return ObitAIPSSetDirname(disk, dir, err);
} /* end AIPSSetDirname */


%}
/* $Id: BeamShape.inc,v 1.1 2008/05/06 13:20:14 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for BeamShape type                          */
/*                                                                    */
/*;  Copyright (C) 2008                                               */
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
#include "ObitBeamShape.h"
#include "ObitImage.h"
%}


%inline %{
extern ObitBeamShape* newBeamShape (char* name) {
  return newObitBeamShape (name);
} // end  newBeamShape

extern ObitBeamShape* BeamShapeCopy  (ObitBeamShape *in, ObitBeamShape *out, 
				    ObitErr *err) {
  return ObitBeamShapeCopy (in, out, err);
} // end  BeamShapeCopy

extern ObitBeamShape* BeamShapeUnref (ObitBeamShape* in) {
  if (!ObitBeamShapeIsA(in)) return NULL;
  return ObitBeamShapeUnref(in);
}

extern ObitBeamShape*  BeamShapeRef (ObitBeamShape* in) {
  return ObitBeamShapeRef(in);
}

extern ObitBeamShape* BeamShapeCreate (char *name, ObitImage *image, 
				       float pbmin, float antSize, 
				       int doGain) {
  gboolean ldoGain = (doGain!=0);
  return ObitBeamShapeCreate((gchar*)name, image, (ofloat)pbmin, (ofloat)antSize,
			     ldoGain);
} // end BeamShapeCreate

extern float BeamShapeGain (ObitBeamShape* in, double ra, double dec,
		            float parAng ) {
  return ObitBeamShapeGain(in, (odouble)ra, (odouble)dec, (ofloat)parAng);
} // end BeamShapeGain

extern float BeamShapeGainSym (ObitBeamShape* in, double Angle) {
  return ObitBeamShapeGainSym(in, (odouble)Angle);
} // end BeamShapeGainSym

extern char* BeamShapeGetName (ObitBeamShape* in) {
  if (ObitBeamShapeIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int BeamShapeIsA (ObitBeamShape* in) {
  return ObitBeamShapeIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitBeamShape *me;
} BeamShape;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitBeamShape *me;
} BeamShape;

%addmethods BeamShape { 
  BeamShape(char* name, ObitImage *image, float pbmin, float antSize,
            int doGain) {
     BeamShape *out;
     gboolean ldoGain = (doGain!=0);
     out = (BeamShape *) malloc(sizeof(BeamShape));
     if (strcmp(name, "None")) out->me = BeamShapeCreate((gchar*)name, image, 
                                                         (ofloat)pbmin, (ofloat)antSize,
                                			  ldoGain);
     else out->me = NULL;
     return out;
   }
  ~BeamShape() {
   if (!self) return;  // Not defined
   if (self && self->me && self->me->ReferenceCount>0) {
      self->me = BeamShapeUnref(self->me);
      free(self);
   }
  }
};

/* $Id: CArray.inc,v 1.7 2006/10/11 17:24:48 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitCArray type                        */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitCArray.h"
%}


%inline %{
extern ObitCArray* CArrayCreate(char* name, long ndim, long *naxis) {
   olong i, lnaxis[10];
   for (i=0; i<ndim; i++) lnaxis[i] = (olong)naxis[i];
   return  ObitCArrayCreate (name, (olong)ndim, lnaxis);
}

extern ObitCArray* CArrayCopy  (ObitCArray *in, ObitCArray *out, ObitErr *err) {
  return ObitCArrayCopy (in, out, err);
} // end CArrayCopy 

extern int CArrayIsCompatable  (ObitCArray *in1, ObitCArray *in2) {
  return ObitCArrayIsCompatable(in1, in2);
}

extern ObitCArray* CArrayRealloc (ObitCArray* in, long ndim, long *naxis) {
   olong i, lnaxis[10];
   for (i=0; i<ndim; i++) lnaxis[i] = (olong)naxis[i];
   return  ObitCArrayRealloc(in, (olong)ndim, lnaxis);
}

extern void CArrayGetVal(ObitCArray* in, long *pos, float val[2]) {
   float *out;
   olong i, lpos[10];

   for (i=0; i<in->ndim; i++) lpos[i] = (olong)pos[i];
   out = ObitCArrayIndex(in, lpos);
   val[0] = out[0];
   val[1] = out[1];
}

extern void CArraySetVal(ObitCArray* in, long *pos, float val[2]) {
   float *off;
   olong i, lpos[10];

   for (i=0; i<in->ndim; i++) lpos[i] = (olong)pos[i];
   off = ObitCArrayIndex (in, lpos);
   off[0] = val[0];
   off[1] = val[1];
}

extern float CArrayMaxAbs (ObitCArray* in, long *outValue2) {
   olong i,loutValue2[10];
   float val;
   val = (float)ObitCArrayMaxAbs(in, loutValue2);
   for (i=0; i<2; i++) outValue2[i] = (long)loutValue2[i];
	
   return val;
}

extern void CArrayNeg (ObitCArray* in) {
   ObitCArrayNeg(in);
}

extern void CArrayConjg (ObitCArray* in) {
   ObitCArrayConjg(in);
}

extern void CArrayFill (ObitCArray* in, float cmpx[2]) {
   ObitCArrayFill(in, cmpx);
}

extern void CArraySAdd (ObitCArray* in, float scalar) {
   ObitCArraySAdd(in, scalar);
}

extern void CArraySMul (ObitCArray* in, float scalar) {
   ObitCArraySMul(in, scalar);
}

extern void CArrayCSAdd (ObitCArray* in, float scalar[2]) {
   ObitCArrayCSAdd(in, scalar);
}

extern void CArrayCSMul (ObitCArray* in, float scalar[2]) {
   ObitCArrayCSMul(in, scalar);
}

extern void CArrayAdd (ObitCArray* in1, ObitCArray* in2, ObitCArray* out) {
   ObitCArrayAdd (in1, in2, out);
}

extern void CArraySub (ObitCArray* in1, ObitCArray* in2, ObitCArray* out) {
   ObitCArraySub (in1, in2, out);
}

extern void CArrayMul (ObitCArray* in1, ObitCArray* in2, ObitCArray* out) {
   ObitCArrayMul (in1, in2, out);
}

extern void CArrayDiv (ObitCArray* in1, ObitCArray* in2, ObitCArray* out) {
   ObitCArrayDiv (in1, in2, out);
}

extern ObitFArray* CArrayMakeF  (ObitCArray *Cin) {
  return ObitCArrayMakeF  (Cin);
}

extern ObitCArray* CArrayMakeC  (ObitFArray *Fin) {
  return ObitCArrayMakeC  (Fin);
}

extern int CArrayIsFCompatable  (ObitCArray *Cin, ObitFArray *Fin) {
  return ObitCArrayIsFCompatable  (Cin, Fin);
}

extern void CArrayFMul (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out) {
  ObitCArrayFMul (Cin, Fin, out);
}

extern void CArrayFDiv (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out) {
  ObitCArrayFDiv (Cin, Fin, out);
}

extern void CArrayFAdd (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out) {
  ObitCArrayFAdd (Cin, Fin, out);
}

extern void CArrayFRot (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out) {
  ObitCArrayFRot (Cin, Fin, out);
}

extern void CArrayComplex (ObitFArray* Fin1, ObitFArray* Fin2, ObitCArray* out) {
  ObitCArrayComplex (Fin1, Fin2, out);
}

extern void CArrayReal (ObitCArray* in, ObitFArray* out) {
  ObitCArrayReal (in, out);
}

extern void CArrayImag (ObitCArray* in, ObitFArray* out) {
  ObitCArrayImag (in, out);
}

extern void CArrayAmp (ObitCArray* in, ObitFArray* out) {
  ObitCArrayAmp (in, out);
}

extern void CArrayPhase (ObitCArray* in, ObitFArray* out) {
  ObitCArrayPhase (in, out);
}

extern void CArray2DCenter(ObitCArray* in) {
  ObitCArray2DCenter (in);
}

extern ObitCArray* CArrayAddConjg (ObitCArray* in, long numConjRow) {
  return ObitCArrayAddConjg (in, (olong)numConjRow);
}

extern char* CArrayGetName (ObitCArray* in) {
  return in->name;
} // end  CArrayGetName

extern long CArrayGetNdim (ObitCArray* in) {
  return in->ndim;
} // end  CArrayGetNdim

// returns an array of 7 elements no matter
extern PyObject* CArrayGetNaxis (ObitCArray* in) {
  long i;
  PyObject *outList= PyList_New(in->ndim);

  for (i=0; i<in->ndim; i++) {
    PyList_SetItem(outList, i, PyInt_FromLong((long)in->naxis[i]));
  }
  return outList;
  return outList;
} // end  CArrayGetNaxis

extern int CArrayIsA (ObitCArray* in) {
  return ObitCArrayIsA(in);
} // end  CArrayIsA 

ObitCArray* CArrayRef (ObitCArray* in) {
  return ObitCArrayRef (in);
} // end CArrayRef

ObitCArray* CArrayUnref (ObitCArray* in) {
  if (!ObitCArrayIsA(in)) return NULL;
  return ObitCArrayUnref (in);
} // end CArrayUnref



%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitCArray *me;
} CArray;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitCArray *me;
} CArray;

%addmethods CArray { 
  CArray(char* name, long ndim, long *naxis) {
     CArray *out;
     out = (CArray *) malloc(sizeof(CArray));
     if (strcmp(name, "None")) out->me = CArrayCreate(name, ndim, naxis);
     else out->me = NULL;
     return out;
   }
  ~CArray() {
    self->me = CArrayUnref(self->me);
    free(self);
  }
};

/* $Id: Catalog.inc,v 1.4 2007/03/20 14:40:14 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Convolution utilities                  */
/*                                                                    */
/*;  Copyright (C) 2006-2008                                          */
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
#include "ObitTableMFUtil.h"
#include "ObitTableVLUtil.h"
%}


%inline %{
/**  Convert an MF table to a VL Table */
void TableMF2VL (ObitTable *in, ObitTable *out, 
		     ObitImage *image, ObitErr *err)
{
  ObitTableMF *inMF=NULL;
  ObitTableVL *outVL=NULL;

  inMF  = ObitTableMFConvert(in);
  outVL = ObitTableVLConvert(out);
  ObitTableMF2VL (inMF, outVL, image, err);
  inMF  = ObitTableMFUnref(inMF);
  outVL = ObitTableVLUnref(outVL);
}  /* end TableMF2VL */

/**  Write human readable version of an MF table to a FILE */
void TableMFPrint (ObitTable *in, ObitImage *image, char *prtFile, 
		   ObitErr *err)
{
  ObitTableMF *inMF=NULL;
  FILE *file;

  if (strncmp (prtFile, "stdout", 6)) {
    ObitTrimTrail(prtFile);  /* Trim any trailing blanks */
    file = fopen (prtFile, "a");
  } else {
    file = stdout;
  }
  inMF  = ObitTableMFConvert(in);
  ObitTableMFPrint (inMF, image, file, err);
  if (strncmp (prtFile, "stdout", 6)) {
    fclose (file);
  }
  inMF  = ObitTableMFUnref(inMF);
}  /* end TableMFPrint */

/**  Append one VL table to another */
void TableVLAppend (ObitTable *in, ObitTable *out, ObitErr *err)
{
  ObitTableVL *inVL=NULL, *outVL=NULL;

  inVL  = ObitTableVLConvert(in);
  outVL = ObitTableVLConvert(out);
  ObitTableVLAppend (inVL, outVL,  err);
  inVL  = ObitTableVLUnref(inVL);
  outVL = ObitTableVLUnref(outVL);
}  /* end TableVLAppend */

/** Index a VL table */
void TableVLIndex (ObitTable *in, ObitErr *err)
{
  ObitTableVL *inVL=NULL;

  inVL  = ObitTableVLConvert(in);
  ObitTableVLIndex (inVL, err);
  inVL  = ObitTableVLUnref(inVL);
}  /* end TableVLIndex */

/**  Merge overlapping components */
void TableVLMerge (ObitTable *in, ObitErr *err)
{
  ObitTableVL *inVL=NULL;

  inVL  = ObitTableVLConvert(in);
  ObitTableVLMerge (inVL, err);
  inVL  = ObitTableVLUnref(inVL);
}  /* end TableVLMerge */

/**  Select significant components */
void TableVLSelect (ObitTable *in, ObitTable *out, ObitErr *err)
{
  ObitTableVL *inVL=NULL, *outVL=NULL;

  inVL  = ObitTableVLConvert(in);
  outVL = ObitTableVLConvert(out);
  ObitTableVLSelect (inVL, outVL,  err);
  inVL  = ObitTableVLUnref(inVL);
  outVL = ObitTableVLUnref(outVL);
}  /* end TableVLSelect */

/**  Remove entries from a given field */
void TableVLPurge (ObitTable *in, char *field, ObitErr *err)
{
  ObitTableVL *inVL=NULL;

  inVL  = ObitTableVLConvert(in);
  ObitTableVLPurge (inVL, field, err);
  inVL  = ObitTableVLUnref(inVL);
}  /* end TableVLPurge */

/**  Remove redundant entries */
void TableVLRedun (ObitTable *in, ObitTable *out, ObitErr *err)
{
  ObitTableVL *inVL=NULL, *outVL=NULL;

  inVL  = ObitTableVLConvert(in);
  outVL = ObitTableVLConvert(out);
  ObitTableVLRedun (inVL, outVL,  err);
  inVL  = ObitTableVLUnref(inVL);
  outVL = ObitTableVLUnref(outVL);
}  /* end  TableVLRedun*/

/**  Write human readable version of an VL table to a FILE */
void TableVLPrint (ObitTable *in, ObitImage *image, char *prtFile, 
		   ObitErr *err)
{
  FILE *file;
  ObitTableVL *inVL=NULL;

  if (strncmp (prtFile, "stdout", 6)) {
    ObitTrimTrail(prtFile);  /* Trim any trailing blanks */
    file = fopen (prtFile, "a");
  } else {
    file = stdout;
  }
  inVL  = ObitTableVLConvert(in);
  ObitTableVLPrint (inVL, image, file, err);
  if (strncmp (prtFile, "stdout", 6)) {
    fclose (file);
  }
  inVL  = ObitTableVLUnref(inVL);
}  /* end TableVLPrint */

/**  Convert a VL  table to a VZ Table */
ObitTable* TableVL2VZ (ObitTable *in, ObitData *data, 
		     	 ObitErr *err)
{
  ObitTableVL *inVL=NULL;
  ObitTableVZ *outVZ=NULL;

  inVL  = ObitTableVLConvert(in);
  outVZ = ObitTableVL2VZ (inVL, data, err);
  inVL  = ObitTableVLUnref(inVL);
  return (ObitTable*)outVZ;
}  /* end TableVL2VZ */

/**  Select entries in VZ  to a VZ Table */
ObitTable* TableVZSel (ObitTable *in, ObitData *data, 
		       ObitErr *err)
{
  ObitTableVZ *inVZ=NULL;
  ObitTableVZ *outVZ=NULL;

  inVZ  = ObitTableVZConvert(in);
  outVZ = ObitTableVZSel (inVZ, data, err);
  inVZ  = ObitTableVZUnref(inVZ);
  return (ObitTable*)outVZ;
}  /* end TableVL2VZ */


%}
/* $Id: CleanImage.inc,v 1.4 2005/08/04 14:15:26 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for CleanImage type                        */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitDConCleanImage.h"
%}


%inline %{
extern ObitDConCleanImage* newCleanImage (char* name) {
  return newObitDConCleanImage (name);
} // end  newCleanImage

extern ObitDConCleanImage* CleanImageCopy  (ObitDConCleanImage *in, 
	ObitDConCleanImage *out, 
        ObitErr *err) {
  return ObitDConCleanImageCopy (in, out, err);
} // end  CleanImageCopy

extern ObitDConCleanImage* CleanImageUnref (ObitDConCleanImage* in) {
  if (!ObitDConCleanImageIsA(in)) return NULL;
  return ObitDConCleanImageUnref(in);
}

extern ObitDConCleanImage*  CleanImageRef (ObitDConCleanImage* in) {
  return ObitDConCleanImageRef(in);
}

extern ObitInfoList* CleanImageGetList (ObitDConCleanImage* in) {
  return ObitInfoListRef(in->info);
}

extern ObitImageMosaic* CleanImageGetImageMosaic (ObitDConCleanImage* in) {
  return ObitImageMosaicRef(in->mosaic);
}

extern void CleanImageSetImageMosaic (ObitDConCleanImage* in, ObitImageMosaic *mosaic) {
  in->mosaic = ObitImageMosaicUnref(in->mosaic);  /* Out with the old */
  in->mosaic = ObitImageMosaicRef(mosaic);        /* In with the new */
}

extern ObitDConCleanWindow* CleanImageGetWindow (ObitDConCleanImage* in) {
  return ObitDConCleanWindowRef(in->window);
}

extern void CleanImageSetWindow (ObitDConCleanImage* in, ObitDConCleanWindow *window) {
  in->mosaic = ObitDConCleanWindowUnref(in->window);  /* Out with the old */
  in->mosaic = ObitDConCleanWindowRef(window);        /* In with the new */
}

// if (win[0]<0 then this is a round window then win[1]=radius, [2,3] = center
// else rectangular and blc=(win[0],win[1]), trc= blc=(win[2],win[3])
extern void CleanImageAddWindow (ObitDConCleanImage* in, int field, int *win, 
                                 ObitErr *err) {
  olong lfield, window[4];
  ObitDConCleanWindowType type;

  // which type? 
  if (win[0]<0) { // circle
    type = OBIT_DConCleanWindow_round;
    window[0] = win[1];
    window[1] = win[2];
    window[2] = win[3];
  } else { // rectangular
    type = OBIT_DConCleanWindow_rectangle;
    window[0] = win[0];
    window[1] = win[1];
    window[2] = win[2];
    window[3] = win[3];
  }
	
  lfield = field;
  ObitDConCleanWindowAdd (in->window, lfield, type, window, err);
}  // end CleanImageAddWindow 

extern ObitDConCleanImage* CleanImageCreate (char *name, ObitImageMosaic* mosaic, 
                                             ObitErr *err) {
 return ObitDConCleanImageCreate(name, mosaic, err);
}

extern void CleanImageDeconvolve (ObitDConCleanImage* in, ObitErr *err) {
 ObitDConCleanImageDeconvolve((ObitDCon*)in, err);
}

extern void CleanImageDefWindow (ObitDConCleanImage* in, ObitErr *err) {
 ObitDConCleanDefWindow((ObitDConClean*)in, err);
}

extern char* CleanImageGetName (ObitDConCleanImage* in) {
  if (ObitDConCleanImageIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int CleanImageIsA (ObitDConCleanImage* in) {
  return ObitDConCleanImageIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitDConCleanImage *me;
} CleanImage;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDConCleanImage *me;
} CleanImage;

%addmethods CleanImage { 
  CleanImage(char* name) {
     CleanImage *out;
     out = (CleanImage *) malloc(sizeof(CleanImage));
     if (strcmp(name, "None")) out->me = newCleanImage(name);
     else out->me = NULL;
     return out;
   }
  ~CleanImage() {
   if (self->me-> ReferenceCount>0) 
      self->me = CleanImageUnref(self->me);
   free(self);
  }
};

/* $Id: CleanVis.inc,v 1.6 2007/11/06 13:00:14 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for CleanVis type                          */
/*                                                                    */
/*;  Copyright (C) 2005-2008                                          */
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
#include "ObitDConCleanVis.h"
%}


%inline %{
extern ObitDConCleanVis* newCleanVis (char* name) {
  return newObitDConCleanVis (name);
} // end  newCleanVis

extern ObitDConCleanVis* CleanVisCopy  (ObitDConCleanVis *in, 
	ObitDConCleanVis *out, 
        ObitErr *err) {
  return ObitDConCleanVisCopy (in, out, err);
} // end  CleanVisCopy

extern ObitDConCleanVis* CleanVisUnref (ObitDConCleanVis* in) {
  if (!ObitDConCleanVisIsA(in)) return NULL;
  return ObitDConCleanVisUnref(in);
}

extern ObitDConCleanVis*  CleanVisRef (ObitDConCleanVis* in) {
  return ObitDConCleanVisRef(in);
}

extern ObitInfoList* CleanVisGetList (ObitDConCleanVis* in) {
  return ObitInfoListRef(in->info);
}

extern ObitImageMosaic* CleanVisGetImageMosaic (ObitDConCleanVis* in) {
  return ObitImageMosaicRef(in->mosaic);
}

extern void CleanVisSetImageMosaic (ObitDConCleanVis* in, ObitImageMosaic *mosaic) {
  in->mosaic = ObitImageMosaicUnref(in->mosaic);  /* Out with the old */
  in->mosaic = ObitImageMosaicRef(mosaic);        /* In with the new */
}

extern ObitSkyModel* CleanVisGetSkyModel (ObitDConCleanVis* in) {
  return ObitSkyModelRef(in->skyModel);
}

extern void CleanVisSetSkyModel (ObitDConCleanVis* in, ObitSkyModel *skymodel) {
  in->skyModel = ObitSkyModelUnref(in->skyModel);  /* Out with the old */
  in->skyModel = ObitSkyModelRef(skymodel);        /* In with the new */
}

extern ObitDConCleanWindow* CleanVisGetWindow (ObitDConCleanVis* in) {
  return ObitDConCleanWindowRef(in->window);
}

extern void CleanVisSetWindow (ObitDConCleanVis* in, ObitDConCleanWindow *window) {
  in->window = ObitDConCleanWindowUnref(in->window);  /* Out with the old */
  in->window = ObitDConCleanWindowRef(window);        /* In with the new */
}

// if (win[0]<0 then this is a round window then win[1]=radius, [2,3] = center
// else rectangular and blc=(win[0],win[1]), trc= blc=(win[2],win[3])
extern void CleanVisAddWindow (ObitDConCleanVis* in, int field, int *win, 
                                 ObitErr *err) {
  olong lfield, window[4];
  ObitDConCleanWindowType type;

  // which type? 
  if (win[0]<0) { // circle
    type = OBIT_DConCleanWindow_round;
    window[0] = win[1];
    window[1] = win[2];
    window[2] = win[3];
  } else { // rectangular
    type = OBIT_DConCleanWindow_rectangle;
    window[0] = win[0];
    window[1] = win[1];
    window[2] = win[2];
    window[3] = win[3];
  }
	
  lfield = field;
  ObitDConCleanWindowAdd (in->window, lfield, type, window, err);
}  // end CleanVisAddWindow 

extern ObitDConCleanVis* CleanVisCreate (char *name, ObitUV* uvdata, 
                                         ObitErr *err) {
 return ObitDConCleanVisCreate(name, uvdata, err);
}

extern void CleanVisDeconvolve (ObitDConCleanVis* in, ObitErr *err) {
 ObitDConCleanVisDeconvolve((ObitDCon*)in, err);
}

extern void CleanVisDefWindow (ObitDConCleanVis* in, ObitErr *err) {
 ObitDConCleanVisDefWindow((ObitDConClean*)in, err);
}

extern int CleanVisReimage (ObitDConCleanVis* in, ObitUV* uvdata, 
                            ObitErr *err) {
  gboolean bout;
  int out;
  bout = ObitDConCleanVisReimage(in, uvdata, err);
  if (bout) out = 1;
  else out = 0;
  return out;
}

extern char* CleanVisGetName (ObitDConCleanVis* in) {
  if (ObitDConCleanVisIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int CleanVisIsA (ObitDConCleanVis* in) {
  return ObitDConCleanVisIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitDConCleanVis *me;
} CleanVis;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDConCleanVis *me;
} CleanVis;

%addmethods CleanVis { 
  CleanVis(char* name) {
     CleanVis *out;
     out = (CleanVis *) malloc(sizeof(CleanVis));
     if (strcmp(name, "None")) out->me = newCleanVis(name);
     else out->me = NULL;
     return out;
   }
  ~CleanVis() {
   if (self->me->ReferenceCount>0) 
      self->me = CleanVisUnref(self->me);
   free(self);
  }
};

/* $Id: ConvUtil.inc,v 1.1 2006/04/18 20:13:23 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Convolution utilities                  */
/*                                                                    */
/*;  Copyright (C) 2006                                               */
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
#include "ObitConvUtil.h"
%}


%inline %{
void ConvUtilConv (ObitImage *inImage, ObitFArray *convFn, 
		   int doDivide, float rescale,
		   ObitImage *outImage, ObitErr *err) {
  gboolean ldoDivide;
  ofloat lrescale = (ofloat)rescale;
  ldoDivide = doDivide != 0;
  ObitConvUtilConv (inImage, convFn, ldoDivide, lrescale, outImage, err);
} // end ConvUtilConv 

ObitFArray* ConvUtilGaus (ObitImage *inImage, float Maj, float Min, float PA) {
  ofloat Beam[3];
  Beam[0] = Maj;
  Beam[0] = Min;
  Beam[0] = PA;
  return ObitConvUtilGaus (inImage, Beam);
} // end ConvUtilGaus

%}
/* $Id: FArray.inc,v 1.19 2007/03/19 13:45:28 bcotton Exp $           */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitFarray type                        */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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

#include "ObitFArray.h"
%}


%inline %{
extern ObitFArray* FArrayCreate(char* name, long ndim, long *naxis) {
   olong i, lnaxis[10];
   for (i=0; i<ndim; i++) lnaxis[i] = (olong)naxis[i];
   return  ObitFArrayCreate (name,(olong) ndim, (olong*)lnaxis);
}

extern float FArrayGetVal(ObitFArray* in, long *pos) {
   float *off, out[4];
   long *Iout;
   olong i, lpos[10];
   for (i=0; i<in->ndim; i++) lpos[i] = (olong)pos[i];
  // check if in bounds
   off = ObitFArrayIndex (in, lpos);
   if (off==NULL) {
	PyErr_SetString(PyExc_RuntimeError,"Position not in array");
        Iout = (long*)out;
        *Iout = ~0;
        return out[0];
   } else {
     out[0] = *off;
   }
   // return NaN rather than Obit magic value
   if (out[0]==ObitMagicF()) {  // create a word with all bits on
        Iout = (long*)out;
        *Iout = ~0;
   }
   return out[0];
}

extern float FArrayGetBlank(void) {
  return ObitMagicF();
}

extern void FArraySetVal(ObitFArray* in, long *pos, float val) {
   float *off;
   olong i, lpos[10];

   for (i=0; i<in->ndim; i++) lpos[i] = (olong)pos[i];

   off = ObitFArrayIndex (in, lpos);
   // check if in bounds
   if (off==NULL) {
      PyErr_SetString(PyExc_RuntimeError,"Position not in array");
      return;
   }
   *off = val;
}

 extern PyObject *FArrayGetBuf (ObitFArray *in) {
   return PyBuffer_FromReadWriteMemory(in->array,
 				      in->arraySize * sizeof(ofloat));
 }
 
extern ObitFArray* FArrayCopy  (ObitFArray *in, ObitFArray *out, ObitErr *err) {
  return ObitFArrayCopy (in, out, err);
} // end FArrayCopy 

extern ObitFArray* FArrayClone  (ObitFArray *in, ObitFArray *out, ObitErr *err) {
  ObitFArrayClone (in, out, err);
  return out;
} // end FArrayCopy 

extern ObitFArray* FArraySubArr  (ObitFArray *in, long *blc, long *trc, 
				  ObitErr *err) {
  olong i, lblc[10], ltrc[10];

  for (i=0; i<in->ndim; i++) lblc[i] = (olong)blc[i];
  for (i=0; i<in->ndim; i++) ltrc[i] = (olong)trc[i];

  return ObitFArraySubArr (in, lblc, ltrc, err);
} // end FArraySubArr

extern ObitFArray* FArrayTranspose  (ObitFArray *in, long *order, ObitErr *err) {
  olong i, lorder[10];
  for (i=0; i<in->ndim; i++) lorder[i] = (olong)order[i];
  return ObitFArrayTranspose (in, lorder, err);
} // end FArrayTranspose

extern int FArrayIsCompatable  (ObitFArray *in1, ObitFArray *in2) {
  return ObitFArrayIsCompatable(in1, in2);
}

extern ObitFArray* FArrayRealloc (ObitFArray* in, long ndim, long *naxis) {
   olong i, lnaxis[10];
   for (i=0; i<ndim; i++) lnaxis[i] = (olong)naxis[i];
   return  ObitFArrayRealloc(in, (olong)ndim, lnaxis);
}

//extern float* FArrayIndex (ObitFArray* in, long *pos) {
//   return ObitFArrayIndex (in, pos);
//}

extern float FArrayMax (ObitFArray* in, long *outValue2) {
   olong i,loutValue2[10];
   float val;
   val = (float)ObitFArrayMax(in, loutValue2);
   for (i=0; i<2; i++) outValue2[i] = (long)loutValue2[i];
   return val;
}

extern float FArrayMaxAbs (ObitFArray* in, long *outValue2) {
   olong i,loutValue2[10];
   float val;
   val = (float)ObitFArrayMaxAbs(in, loutValue2);
   for (i=0; i<2; i++) outValue2[i] = (long)loutValue2[i];
   return val;
}

extern float FArrayMin (ObitFArray* in, long *outValue2) {
   olong i,loutValue2[10];
   float val;
   val = (float)ObitFArrayMin(in, loutValue2);
   for (i=0; i<2; i++) outValue2[i] = (long)loutValue2[i];
   return val;
}

extern void FArrayDeblank (ObitFArray* in, float scalar) {
   ObitFArrayDeblank (in, (ofloat)scalar);
}

extern float FArrayRMS (ObitFArray* in) {
   float out[4];
   long *Iout;
   out[0] = ObitFArrayRMS(in);
   // return NaN rather than Obit magic value
    if (out[0]==ObitMagicF()) {  // create a word with all bits on
        Iout = (long*)&out;
        *Iout = ~0;
   }
   return out[0];
}

extern float FArrayRMS0 (ObitFArray* in) {
   float out[4];
   long *Iout;
   out[0] = ObitFArrayRMS0(in);
   // return NaN rather than Obit magic value
    if (out[0]==ObitMagicF()) {  // create a word with all bits on
        Iout = (long*)&out;
        *Iout = ~0;
   }
   return out[0];
}

extern float FArrayRawRMS (ObitFArray* in) {
   float out[4];
   long *Iout;
   out[0] = ObitFArrayRawRMS(in);
   // return NaN rather than Obit magic value
    if (out[0]==ObitMagicF()) {  // create a word with all bits on
        Iout = (long*)&out;
        *Iout = ~0;
   }
   return out[0];
}

extern float FArrayMode (ObitFArray* in) {
   float out[4];
   long *Iout;
   out[0] = ObitFArrayMode(in);
   // return NaN rather than Obit magic value
   if (out[0]==ObitMagicF()) {  // create a word with all bits on
        Iout = (long*)&out;
        *Iout = ~0;
   }
   return out[0];
}

extern float FArrayMean (ObitFArray* in) {
   float out[4];
   long *Iout;
   out[0] = ObitFArrayMean(in);
   // return NaN rather than Obit magic value
   if (out[0]==ObitMagicF()) {  // create a word with all bits on
        Iout = (long*)&out;
        *Iout = ~0;
   }
   return out[0];
}

extern void FArrayFill (ObitFArray* in, float scalar) {
   return ObitFArrayFill(in, (ofloat)scalar);
}

extern void FArrayNeg (ObitFArray* in) {
   ObitFArrayNeg(in);
}

extern void FArraySin (ObitFArray* in) {
   ObitFArraySin(in);
}

extern void FArrayCos (ObitFArray* in) {
   ObitFArrayCos(in);
}

extern float FArraySum (ObitFArray* in) {
   return ObitFArraySum(in);
}

extern long FArrayCount (ObitFArray* in) {
   return ObitFArrayCount(in);
}

extern void FArraySAdd (ObitFArray* in, float scalar) {
   ObitFArraySAdd(in, (ofloat)scalar);
}

extern void FArraySMul (ObitFArray* in, float scalar) {
   ObitFArraySMul(in, (ofloat)scalar);
}

extern void FArraySDiv (ObitFArray* in, float scalar) {
   ObitFArraySDiv(in, (ofloat)scalar);
}

extern void FArrayClip (ObitFArray* in, float minVal,float maxVal, float newVal) {
   ObitFArrayClip(in, (ofloat)minVal, (ofloat)maxVal, (ofloat)newVal);
}

extern void FArrayInClip (ObitFArray* in, float minVal,float maxVal, float newVal) {
   ObitFArrayInClip(in, (ofloat)minVal, (ofloat)maxVal, (ofloat)newVal);
}

extern void FArrayDivClip (ObitFArray* in1, ObitFArray* in2, float minVal, 
	ObitFArray* out) {
   ObitFArrayDivClip(in1, in2, (ofloat)minVal, out);
}

extern void FArrayClipBlank (ObitFArray* in, float minVal,float maxVal) {
   ofloat fblank = ObitMagicF();
   ObitFArrayClip(in, (ofloat)minVal, (ofloat)maxVal, (ofloat)fblank);
}

extern void FArrayBlank (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayBlank (in1, in2, out);
}

extern void FArraySumArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArraySumArr (in1, in2, out);
}

extern void FArrayAvgArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayAvgArr (in1, in2, out);
}

extern void FArrayMaxArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayMaxArr (in1, in2, out);
}

extern void FArrayMinArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayMinArr (in1, in2, out);
}

extern void FArrayAdd (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayAdd (in1, in2, out);
}

extern void FArraySub (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArraySub (in1, in2, out);
}

extern void FArrayMul (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayMul (in1, in2, out);
}

extern void FArrayDiv (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayDiv (in1, in2, out);
}

extern float FArrayDot (ObitFArray* in1, ObitFArray* in2) {
	return ObitFArrayDot(in1, in2);
}

extern void FArrayMulColRow (ObitFArray* in, ObitFArray* row, ObitFArray* col,
				  ObitFArray* out) {
   ObitFArrayMulColRow (in, row, col, out);
}

extern void FArray2DCenter (ObitFArray* in) {
   ObitFArray2DCenter (in);
}

extern int FArraySymInv2D (ObitFArray* in) {
   int ierr;
   ObitFArray2DSymInv (in, &ierr);
   return ierr;
}

extern void FArrayCGauss2D (ObitFArray* in, long Cen[2], float FWHM) {
   olong lCen[2] = {Cen[0],Cen[1]};
   ObitFArray2DCGauss (in, lCen, (ofloat)FWHM);
}

extern void FArrayEGauss2D (ObitFArray* in, float amp, float Cen[2], float GauMod[3]) {
   ofloat lCen[2] = {Cen[0],Cen[1]};
   ofloat lGauMod[3] = {GauMod[0], GauMod[1], GauMod[2]};
   ObitFArray2DEGauss (in, (float)amp, lCen, lGauMod);
}

extern void FArrayShiftAdd (ObitFArray* in1, long *pos1,
				 ObitFArray* in2, long *pos2, 
				 float scalar, ObitFArray* out) {
   olong i, lpos1[10], lpos2[10];
   for (i=0; i<in1->ndim; i++) lpos1[i] = (olong)pos1[i];
   for (i=0; i<in2->ndim; i++) lpos2[i] = (olong)pos2[i];
   ObitFArrayShiftAdd (in1, lpos1, in2, lpos2, (float)scalar, out);
} // end FArrayShiftAdd

extern void FArrayPad (ObitFArray* in, ObitFArray* out, float factor) {
   ObitFArrayPad (in, out, (float)factor);
} // end FArrayPad

extern void FArraySelInc  (ObitFArray *in, ObitFArray *out, 
	                   long *blc, long *trc, long *inc, ObitErr *err) {
  olong i, lblc[10], ltrc[10], linc[10];

  for (i=0; i<in->ndim; i++) lblc[i] = (olong)blc[i];
  for (i=0; i<in->ndim; i++) ltrc[i] = (olong)trc[i];
  for (i=0; i<in->ndim; i++) linc[i] = (olong)inc[i];
  ObitFArraySelInc (in, out, lblc, ltrc, linc, err);
} // end FArraySelInc

extern char* FArrayGetName (ObitFArray* in) {
  if (!in) return "Undefined";
  return in->name;
} // end  FArrayGetName

extern long FArrayGetNdim (ObitFArray* in) {
  return in->ndim;
} // end  FArrayGetNdim

// returns an array 
extern  PyObject* FArrayGetNaxis (ObitFArray* in) {
  long i;
  PyObject *outList= PyList_New(in->ndim);

  for (i=0; i<in->ndim; i++) {
    PyList_SetItem(outList, i, PyInt_FromLong((long)in->naxis[i]));
  }
  return outList;
} // end  FArrayGetNaxis

extern int FArrayIsA (ObitFArray* in) {
  return ObitFArrayIsA(in);
} // end  FArrayIsA 

ObitFArray* FArrayRef (ObitFArray* in) {
  return ObitFArrayRef (in);
} // end FArrayRef

ObitFArray* FArrayUnref (ObitFArray* in) {
  if (!ObitFArrayIsA(in)) return NULL;
  if (in && (in->ReferenceCount>0)) in = ObitFArrayUnref (in);
  return in;
} // end FArrayUnref

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitFArray *me;
} FArray;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitFArray *me;
} FArray;

%addmethods FArray { 
  FArray(char* name, long ndim, long *naxis) {
     FArray *out;
     out = (FArray *) malloc(sizeof(FArray));
     if (strcmp(name, "None")) out->me = FArrayCreate(name, ndim, naxis);
     else  out->me = NULL;
     return out;
   }
  ~FArray() {
    self->me = FArrayUnref(self->me);
    free(self);
  }
};

/* $Id: FArrayUtil.inc,v 1.1 2006/04/18 20:28:11 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitFarray type                        */
/*                                                                    */
/*;  Copyright (C) 2005                                               */
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

#include "ObitFArrayUtil.h"
%}


%inline %{
/* Return list with [0]=FWHM, [1] = peak, [2] = cenX, [3]=cenY, 4=RMS residual */
extern PyObject* FArrayUtilFitCGauss (ObitFArray *in, float FWHM, float center[2], 
			      float peak, ObitErr *err)
{
  ofloat lFWHM=FWHM, lcenter[2]={center[0],center[1]}, lpeak=peak;
  ofloat RMS;
  PyObject *o;

  RMS = ObitFArrayUtilFitCGauss (in, &lFWHM, lcenter, &lpeak, err);
  // return list
  o = PyList_New(5);
  PyList_SetItem(o, 0, PyFloat_FromDouble((double)lFWHM));
  PyList_SetItem(o, 1, PyFloat_FromDouble((double)lpeak));
  PyList_SetItem(o, 2, PyFloat_FromDouble((double)lcenter[0]));
  PyList_SetItem(o, 3, PyFloat_FromDouble((double)lcenter[1]));
  PyList_SetItem(o, 4, PyFloat_FromDouble((double)RMS));
  return o;
} // end  FArrayUtilFitCGauss

// Convolution 
ObitFArray* FArrayUtilConvolve (ObitFArray *in1, ObitFArray *in2, 
	 		        ObitErr *err)
{
  return ObitFArrayUtilConvolve (in1, in2, err);
}  // end FArrayUtilConvolve
%}


/* $Id: FFT.inc,v 1.3 2005/02/06 02:00:38 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitFFT type                           */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitFFT.h"
%}


%inline %{
// dir   1 => OBIT_FFT_Forward(R2C), 2 => OBIT_FFT_Reverse(C2R)
// type  1 => OBIT_FFT_FullComplex,  2 => OBIT_FFT_HalfComplex
extern ObitFFT* FFTCreate(char* name, int dir, int type, int rank, int *dim) {
   ObitFFTdir ldir;
   ObitFFTtype ltype;
   if (dir==1) ldir = OBIT_FFT_Forward;
   else ldir = OBIT_FFT_Reverse;
   if (type==1) ltype = OBIT_FFT_FullComplex;
   else ltype = OBIT_FFT_HalfComplex;
   return newObitFFT  (name, ldir, ltype, rank, dim);
}

extern int FFTSuggestSize  (int length) {
  return ObitFFTSuggestSize(length);
}

extern void FFTR2C(ObitFFT* in, ObitFArray *inArray, ObitCArray *outArray) {
   ObitFFTR2C(in, inArray, outArray);
}

extern void FFTC2R(ObitFFT* in, ObitCArray *inArray, ObitFArray *outArray) {
   ObitFFTC2R(in, inArray, outArray);
}

extern void FFTC2C(ObitFFT* in, ObitCArray *inArray, ObitCArray *outArray) {
   ObitFFTC2C(in, inArray, outArray);
}

extern char* FFTGetName (ObitFFT* in) {
  return in->name;
} // end  FFTGetName

extern int FFTGetRank (ObitFFT* in) {
  return in->rank;
} // end  FFTGetRank

// returns an array of 7 elements no matter
extern int* FFTGetDim (ObitFFT* in) {
  return in->dim;
} // end  FFTGetDim

extern int FFTIsA (ObitFFT* in) {
  return ObitFFTIsA(in);
} // end  FFTIsA 

ObitFFT* FFTRef (ObitFFT* in) {
  return ObitFFTRef (in);
} // end FFTRef

ObitFFT* FFTUnref (ObitFFT* in) {
  if (!ObitFFTIsA(in)) return NULL;
  return ObitFFTUnref (in);
} // end FFTUnref



%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitFFT *me;
} FFT;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitFFT *me;
} FFT;

%addmethods FFT { 
  FFT(char* name, int dir, int type, int rank, int *dim) {
     FFT *out;
     /* just create  Python structure here */
     out = (FFT *) malloc(sizeof(FFT));
     if (strcmp(name, "None")) out->me = FFTCreate(name, dir, type, rank, dim);
     else out->me = NULL;
     return out;
   }
  ~FFT() {
    self->me = FFTUnref(self->me);
    free(self);
  }
};

/* $Id: FInterpolate.inc,v 1.3 2005/02/06 02:00:38 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for FInterpolate type                      */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitFInterpolate.h"
%}

%inline %{
extern ObitFInterpolate* 
FInterpolateCreate (char* name, ObitFArray *array, ObitImageDesc *desc, long hwidth) {
  return newObitFInterpolateCreate (name, array, desc, hwidth);
} // end FInterpolateCreate

extern ObitFInterpolate* FInterpolateCopy  (ObitFInterpolate *in, ObitFInterpolate *out, 
			   ObitErr *err) {
  return ObitFInterpolateCopy (in, out, err);
} // end FInterpolateCopy


extern ObitFInterpolate* FInterpolateClone (ObitFInterpolate *in, ObitFInterpolate *out) {
  return ObitFInterpolateClone (in, out);
} // end FInterpolateClone 

extern void FInterpolateReplace (ObitFInterpolate *in, ObitFArray *newArray) {
  ObitFInterpolateReplace (in, newArray);
} // end FInterpolateReplac

extern float FInterpolatePixel (ObitFInterpolate *in, float *pixel, ObitErr *err) {
  return ObitFInterpolatePixel (in, pixel, err);
} // end FInterpolatePixel 

extern float FInterpolate1D (ObitFInterpolate *in, float pixel) {
  return ObitFInterpolate1D (in, pixel);
} // end FInterpolate1D

extern float FInterpolatePosition (ObitFInterpolate *in, double *coord, ObitErr *err) {
  return ObitFInterpolatePosition (in, coord, err);
} // end FInterpolatePosition

extern ObitInfoList* FInterpolateGetList (ObitFInterpolate* in) {
  return ObitInfoListRef(in->info);
}

extern ObitFArray* FInterpolateGetFArray (ObitFInterpolate* in) {
  return ObitFArrayRef(in->myArray);
}

extern ObitImageDesc* FInterpolateGetDesc (ObitFInterpolate* in) {
  return ObitImageDescRef(in->myDesc);
}

extern void FInterpolateSetDesc (ObitFInterpolate* in, ObitImageDesc *desc) {
  in->myDesc = ObitImageDescUnref(in->myDesc);
  in->myDesc = ObitImageDescRef(desc);
}

extern long FInterpolateGetHwidth (ObitFInterpolate* in) {
  return in->hwidth;
}

extern void FInterpolateSetHwidth (ObitFInterpolate* in, long hwidth) {
  in->hwidth = hwidth;
}

ObitFInterpolate* FInterpolateRef (ObitFInterpolate* in) {
  return ObitFInterpolateRef (in);
} // end FInterpolateRef

ObitFInterpolate* FInterpolateUnref (ObitFInterpolate* in) {
  if (!ObitFInterpolateIsA(in)) return NULL;
  return ObitFInterpolateUnref (in);
} // end FInterpolateUnref

extern int FInterpolateIsA (ObitFInterpolate* in) {
  return ObitFInterpolateIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitFInterpolate *me;
} FInterpolate;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitFInterpolate *me;
} FInterpolate;

%addmethods FInterpolate { 
  FInterpolate(char* name, ObitFArray *array, ObitImageDesc *desc, long hwidth) {
     FInterpolate *out;
     /* just create  Python structure here */
     out = (FInterpolate *) malloc(sizeof(FInterpolate));
     if (strcmp(name, "None")) out->me = FInterpolateCreate(name, array, desc, hwidth);
     else out->me = NULL;
     return out;
   }
  ~FInterpolate() {
    self->me = FInterpolateUnref(self->me);
    free(self);
  }
};


/* $Id: FitModel.inc,v 1.3 2007/10/11 13:36:44 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for FitModel type                          */
/*                                                                    */
/*;  Copyright (C) 2007                                               */
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
#include "ObitFitModel.h"
#include "ObitImageFitData.h"
#include "ObitImageDesc.h"
%}


%inline %{
extern ObitFitModel* newFitModel (char* name) {
  return newObitFitModel (name);
} // end  newFitModel

extern ObitFitModel* FitModelCopy  (ObitFitModel *in, ObitFitModel *out, 
				    ObitErr *err) {
  return ObitFitModelCopy (in, out, err);
} // end  FitModelCopy

extern ObitFitModel* FitModelUnref (ObitFitModel* in) {
  if (!ObitFitModelIsA(in)) return NULL;
  return ObitFitModelUnref(in);
}

extern ObitFitModel*  FitModelRef (ObitFitModel* in) {
  return ObitFitModelRef(in);
}

extern ObitFitModel* FitModelCreate (char *name, int type, float Peak, float DeltaX, float DeltaY, 
			             int nparm, float *parms) {
  return ObitFitModelCreate((gchar*)name, (olong)type, (ofloat)Peak, (ofloat)DeltaX, (ofloat)DeltaY, 
			   (olong)nparm, (ofloat*)parms);
}

extern int FitModelGetType(ObitFitModel* in) {
  return (int)in->type;
}  // end FitModelGetType

extern float FitModelGetPeak(ObitFitModel* in) {
  return (float)in->Peak;
}  // end FitModelGetPeak(

extern float FitModelGetDeltaX(ObitFitModel* in) {
  return (float)in->DeltaX;
}  // end FitModelGetDeltaX

extern float FitModelGetDeltaY(ObitFitModel* in) {
  return (float)in->DeltaY;
}  // end FitModelGetDeltaY

extern int FitModelGetNparm(ObitFitModel* in) {
  return (int)in->nparm;
}  // end FitModelGetNparm

extern PyObject*  FitModelGetParms(ObitFitModel* in) {
  PyObject *outList, *o;
  int i;

  outList = PyList_New(in->nparm);
  for (i=0; i<in->nparm; i++) {
    o = PyFloat_FromDouble((double)in->parms[i]);
    PyList_SetItem(outList, i, o);
  }
  return outList;
}  // end FitModelGetParms

extern float FitModelGetePeak(ObitFitModel* in) {
  return (float)in->ePeak;
}  // end FitModelGetePeak

extern float FitModelGeteDeltaX(ObitFitModel* in) {
  return (float)in->eDeltaX;
}  // end FitModelGeteDeltaX

extern float FitModelGeteDeltaY(ObitFitModel* in) {
  return (float)in->eDeltaY;
}  // end FitModelGeteDeltaY

extern PyObject* FitModelGeteParms(ObitFitModel* in) {
  PyObject *outList, *o;
  int i;

  outList = PyList_New(in->nparm);
  for (i=0; i<in->nparm; i++) {
    o = PyFloat_FromDouble((double)in->eparms[i]);
    PyList_SetItem(outList, i, o);
  }
  return outList;
}  // end FitModelGeteParms

extern void FitModelSetType(ObitFitModel* in, int value) {
  in->type = (ObitFitModelCompType)value;
}  // end FitModelSetType

extern void FitModelSetPeak(ObitFitModel* in, float value) {
  in->Peak = (ofloat)value;
}  // end FitModelSetPeak

extern void FitModelSetDeltaX(ObitFitModel* in, float value) {
  in->DeltaX = (ofloat)value;
}  // end FitModelSetDeltaX

extern void FitModelSetDeltaY(ObitFitModel* in, float value) {
  in->DeltaY = (ofloat)value;
}  // end FitModelSetDeltaY

extern void FitModelSetNparm(ObitFitModel* in, int value) {
  in->nparm = (olong)value;
  if (in->parms)  g_free(in->parms);
  if (in->eparms) g_free(in->eparms);
  if (in->nparm>0) {
    in->parms  = g_malloc(in->nparm*sizeof(ofloat));
    in->eparms = g_malloc(in->nparm*sizeof(ofloat));
  } else {
    in->parms  = NULL;
    in->eparms = NULL;
  }
}  // end FitModelSetNparm

extern void FitModelSetParms(ObitFitModel* in, float *value) {
  int i;

  for (i=0; i<in->nparm; i++) {
    in->parms[i] = (ofloat)value[i];
  }
}  // end FitModelSetParms

extern void FitModelSetePeak(ObitFitModel* in, float value) {
  in->ePeak = (ofloat)value;
}  // end FitModelSetePeak

extern void FitModelSeteDeltaX(ObitFitModel* in, float value) {
  in->eDeltaX = (ofloat)value;
}  // end FitModelSeteDeltaX

extern void FitModelSeteDeltaY(ObitFitModel* in, float value) {
  in->eDeltaY = (ofloat)value;
}  // end FitModelSeteDeltaY

extern void FitModelSeteParms(ObitFitModel* in, float *value) {
  int i;

  for (i=0; i<in->nparm; i++) {
    in->eparms[i] = (ofloat)value[i];
  }
}  // end FitModelSeteParms

extern ObitFitModel* DeconGau (ObitFitModel* in, ObitFitModel* out, ObitImageDesc *desc) {
   ofloat dgau[3][3], emaj, emin, epa, xcell;
   olong i, ret;

   /* Check */
    if ((in->type!=OBIT_FitModel_GaussMod) || (in->nparm!=3)) {
 	PyErr_SetString(PyExc_TypeError,"Input model not a Gaussian");
        return NULL;
    }
   xcell = fabs (desc->cdelt[1]);	
   ret = ObitFitModelDeconGau ((ofloat)in->parms[0],   (ofloat)in->parms[1],   (ofloat)(in->parms[2]*RAD2DG), 
			       (ofloat)in->eparms[0],  (ofloat)in->eparms[1],  (ofloat)(in->eparms[2]*RAD2DG),
			       (ofloat)desc->beamMaj/xcell,  (ofloat)desc->beamMin/xcell,  
	                       (ofloat)desc->beamPA, dgau);
  if (ret>1) {
 	//PyErr_SetString(PyExc_RuntimeError,"Deconvolution failed");
	 out->type  = -10;
        return out;
  }
   dgau[0][2]  *= DG2RAD;
   dgau[1][2]  *= DG2RAD;
   dgau[2][2]  *= DG2RAD;
   emaj = fabs(dgau[1][0] - dgau[2][0]) * 0.5;
   emin = fabs(dgau[1][1] - dgau[2][1]) * 0.5;
   epa  = fabs(dgau[1][2] - dgau[2][2]) * 0.5;
   if (epa>G_PI)  epa -= 2.0*G_PI;
   if (epa<-G_PI) epa += 2.0*G_PI;
   out->type    = in->type;
   out->nparm   = in->nparm;
   out->Peak    = in->Peak;
   out->DeltaX  = in->DeltaX;
   out->DeltaY  = in->DeltaY;
   out->ePeak   = in->ePeak;
   out->eDeltaX = in->eDeltaX;
   out->eDeltaY = in->eDeltaY;
   for (i=0; i<3; i++) out->parms[i] = dgau[0][i];
   out->eparms[0] = emaj;
   out->eparms[1] = emin;
   out->eparms[2] = epa;

   return out;
} // end DeconGau

// Determine errors for Gaussian model
extern void FitModelGaussErr (ObitFitModel* in, ObitImageDesc *desc, float irms) {
   ofloat beam[3], xcell;

   /* Check that Gaussian */
    if ((in->type!=OBIT_FitModel_GaussMod) || (in->nparm!=3)) {
 	PyErr_SetString(PyExc_TypeError,"Input model not a Gaussian");
        return;
    }
   xcell = fabs (desc->cdelt[1]);	
   beam[0] = (ofloat)desc->beamMaj/xcell;
   beam[1] = (ofloat)desc->beamMin/xcell;
   beam[2] = (ofloat)desc->beamPA;

   ObitImageFitDataGaussErr ((ofloat)in->Peak,
                             (ofloat)in->parms[0], (ofloat)in->parms[1], (ofloat)(in->parms[2]), 
			     (ofloat)irms,  (ofloat*)&beam[0],
			     (ofloat*)&in->ePeak, (ofloat*)&in->eDeltaX,  (ofloat*)&in->eDeltaY, 
			     (ofloat*)&in->eparms[0], (ofloat*)&in->eparms[1],  (ofloat*)&in->eparms[2]);
   return;
} // end FitModelGaussErr

extern char* FitModelGetName (ObitFitModel* in) {
  if (ObitFitModelIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int FitModelIsA (ObitFitModel* in) {
  return ObitFitModelIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitFitModel *me;
} FitModel;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitFitModel *me;
} FitModel;

%addmethods FitModel { 
  FitModel(char* name, int type, float Peak, float DeltaX, float DeltaY, int nparm, float *parms) {
     FitModel *out;
     out = (FitModel *) malloc(sizeof(FitModel));
     if (strcmp(name, "None")) out->me = FitModelCreate(name, type, Peak, DeltaX, DeltaY, 
			                                nparm, parms);
     else out->me = NULL;
     return out;
   }
  ~FitModel() {
   if (!self) return;  // Not defined
   if (self && self->me && self->me->ReferenceCount>0) {
      self->me = FitModelUnref(self->me);
      free(self);
   }
  }
};

/* $Id: FitRegion.inc,v 1.2 2007/09/20 03:11:58 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for FitRegion type                          */
/*                                                                    */
/*;  Copyright (C) 2007                                               */
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
#include "ObitFitRegion.h"
#include "ObitFitModel.h"
%}


%inline %{
extern ObitFitRegion* newFitRegion (char* name) {
  return newObitFitRegion (name);
} // end  newFitRegion

extern ObitFitRegion* FitRegionCopy  (ObitFitRegion *in, ObitFitRegion *out, 
				    ObitErr *err) {
  return ObitFitRegionCopy (in, out, err);
} // end  FitRegionCopy

extern ObitFitRegion* FitRegionUnref (ObitFitRegion* in) {
  if (!ObitFitRegionIsA(in)) return NULL;
  return ObitFitRegionUnref(in);
}

extern ObitFitRegion*  FitRegionRef (ObitFitRegion* in) {
  return ObitFitRegionRef(in);
}

extern ObitFitRegion* FitRegionCreate (char *name, int corner[2], int dim[2],
		     float peak, float peakResid, float RMSResid,
		     float fluxResid) {
  olong i, lcorner[2], ldim[2];
  ObitFitModel **lmodels=NULL;

  lcorner[0] = (olong)corner[0];
  lcorner[1] = (olong)corner[1];
  ldim[0] = (olong)dim[0];
  ldim[1] = (olong)dim[1];

  return ObitFitRegionCreate((gchar*)name, lcorner, ldim,
		     (ofloat)peak, (ofloat)peakResid, (ofloat)RMSResid,
		     (ofloat)fluxResid, 0, NULL);
}

extern void FitRegionSetCorner(ObitFitRegion* in, int *value){
  in->corner[0] = (olong)value[0];
  in->corner[1] = (olong)value[1];
}  // end FitRegionSetCorner

extern void FitRegionSetDim(ObitFitRegion* in, int *value){
  in->dim[0] = (olong)value[0];
  in->dim[1] = (olong)value[1];
}  // end FitRegionSetDim

extern void FitRegionSetPeak(ObitFitRegion* in, float value){
  in->peak = (ofloat)value;
}  // end FitRegionSetPeak

extern void FitRegionSetPeakResid(ObitFitRegion* in, float value){
  in->peakResid = (ofloat)value;
}  // end FitRegionSetPeakResid

extern void FitRegionSetRMSResid(ObitFitRegion* in, float value){
  in->RMSResid = (ofloat)value;
}  // end FitRegionSetRMSResid

extern void FitRegionSetFluxResid(ObitFitRegion* in, float value){
  in->fluxResid = (ofloat)value;
}  // end FitRegionSetFluxResid

extern void FitRegionSetNmodel(ObitFitRegion* in, int value){
  olong i;
  in->nmodel = (olong)value;
  // resize models
  if (in->models) g_free(in->models);
  in->models = g_malloc0(in->nmodel*sizeof(ObitFitModel*));
  for (i=0; i<in->nmodel; i++) in->models[i] = NULL;
}  // end FitRegionSetNmodel

extern void FitRegionSetModels(ObitFitRegion* in, ObitFitModel *value, int i){
  if (i>=in->nmodel) {
 	PyErr_SetString(PyExc_RuntimeError,"Attempt to set invalid model no.");
        return;
  }
  in->models[i] = ObitFitModelRef(value);
}  // end FitRegionSetModels

extern PyObject*  FitRegionGetCorner(ObitFitRegion* in){
  PyObject *outList, *o;
  int i;

  outList = PyList_New(2);
  for (i=0; i<2; i++) {
    o = PyInt_FromLong((long)in->corner[i]);
    PyList_SetItem(outList, i, o);
  }
  return outList;
}  // end FitRegionGetCorner

extern PyObject*  FitRegionGetDim(ObitFitRegion* in){
  PyObject *outList, *o;
  int i;

  outList = PyList_New(2);
  for (i=0; i<2; i++) {
    o = PyInt_FromLong((long)in->dim[i]);
    PyList_SetItem(outList, i, o);
  }
  return outList;
}  // end FitRegionGetDim

extern float FitRegionGetPeak(ObitFitRegion* in){
  return (float)in->peak;
}  // end FitRegionGetPeak

extern float FitRegionGetPeakResid(ObitFitRegion* in){
  return (float)in->peakResid;
}  // end FitRegionGetPeakResid

extern float FitRegionGetRMSResid(ObitFitRegion* in){
  return (float)in->RMSResid;
}  // end FitRegionGetRMSResid

extern float FitRegionGetFluxResid(ObitFitRegion* in){
  return (float)in->fluxResid;
}  // end FitRegionGetFluxResid

extern int FitRegionGetNmodel(ObitFitRegion* in){
  return (olong)in->nmodel;
}  // end FitRegionGetNmodel

extern ObitFitModel* FitRegionGetModels(ObitFitRegion* in, int i) {
  return ObitFitModelRef(in->models[i]);
}  // end FitRegionGetModels

extern char* FitRegionGetName (ObitFitRegion* in) {
  if (ObitFitRegionIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int FitRegionIsA (ObitFitRegion* in) {
  return ObitFitRegionIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitFitRegion *me;
} FitRegion;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitFitRegion *me;
} FitRegion;

%addmethods FitRegion { 
  FitRegion(char *name, int corner[2], int dim[2], float peak, float peakResid, 
            float RMSResid, float fluxResid) {
     FitRegion *out;
     out = (FitRegion *) malloc(sizeof(FitRegion));
     if (strcmp(name, "None")) out->me = FitRegionCreate(name, corner, dim, peak, 
                                  peakResid, RMSResid, fluxResid);
     else out->me = NULL;
     return out;
   }
  ~FitRegion() {
   if (!self) return;  // Not defined
   if (self->me-> ReferenceCount>0) {
      self->me = FitRegionUnref(self->me);
      free(self);
   }
  }
};

/* $Id: FITSFile.inc,v 1.2 2007/09/11 12:38:38 bcotton Exp $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for FITS file utilities                    */
/*                                                                    */
/*;  Copyright (C) 2007                                               */
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
#include "ObitFileFITS.h"
#include "ObitFITS.h"
%}


%inline %{
extern int FITSFileExist(int disk, char *filename, ObitErr *err)
{
  gchar *fullname=NULL, *fullname2=NULL;
  gboolean exist;
  int lexist=0;

  // Full file name
  fullname = ObitFITSFilename (disk, filename, err);
  if (fullname==NULL) return lexist;

  exist = ObitFileFITSExist (fullname, err);

  // If this didn't work try with ".gz"
  if (!exist) {
    fullname2 = g_strconcat (fullname, ".gz", NULL);
    exist = ObitFileFITSExist (fullname2, err);
  }
  if (exist) lexist = 1;
  else lexist = 0;

  if (fullname)  g_free(fullname);
  if (fullname2) g_free(fullname2);
  return lexist;
} /* end FITSFileExist */

extern int FITSAddDir(char *dir, ObitErr *err)
{
   return ObitFITSAddDir(dir, err);
} /* end FITSSetDirname */

%}
/* $Id: History.inc,v 1.6 2006/01/18 20:02:40 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for History utilities                      */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitHistory.h"
%}

%inline %{

extern ObitHistory* HistoryCreate (gchar* name) {
  return newObitHistory (name);
} // end  HistoryCreate

extern ObitHistory* HistoryZap  (ObitHistory *in, ObitErr *err) {
  return ObitHistoryZap (in, err);
} // end HistoryZap

extern ObitHistory* HistoryCopy  (ObitHistory *in, ObitHistory *out, 
			      ObitErr *err) {

  return ObitHistoryCopy (in, out, err);
} // end  HistoryCopy

extern int HistoryCopyHeader  (ObitHistory *in, ObitHistory *out, 
		 	       ObitErr *err) {

  ObitIOCode ret;
  ret = ObitHistoryCopyHeader (in, out, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  HistoryCopy

extern int HistoryCopy2Header  (ObitHistory *in, ObitHistory *out, 
		 	       ObitErr *err) {

  ObitIOCode ret;
  ret = ObitHistoryCopy2Header (in, out, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  HistoryCopy2Header

extern int HistoryHeader2Header  (ObitHistory *in, ObitHistory *out, 
		 	          ObitErr *err) {

  ObitIOCode ret;
  ret = ObitHistoryHeader2Header (in, out, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  HistoryHeader2Header

extern int HistoryOpen (ObitHistory* in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitHistoryOpen (in, laccess, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
}

extern int HistoryClose (ObitHistory* in, ObitErr *err) {
  ObitIOCode ret;
  ret =  ObitHistoryClose(in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
}

extern PyObject* HistoryReadRec (ObitHistory* in, long recno,  ObitErr *err) {
  ObitIOCode ret;
  gchar hiCard[73];
  ret = ObitHistoryReadRec(in, recno, hiCard, err);
  if (ret==OBIT_IO_OK) return PyString_FromString(hiCard);
  return PyString_FromString("");
}

extern int HistoryWriteRec (ObitHistory* in, long recno, char hiCard[73],
		           ObitErr *err) {
  ObitIOCode ret;
  ret = ObitHistoryWriteRec(in, recno, hiCard, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
}

extern int HistoryTimeStamp (ObitHistory* in, char *label, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitHistoryTimeStamp(in, label, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
}

extern ObitHistory* HistoryUnref (ObitHistory* in) {
  if (!ObitHistoryIsA(in)) return NULL;
  return ObitHistoryUnref(in);
}

extern ObitHistory*  HistoryRef (ObitHistory* in) {
  return ObitHistoryRef(in);
}

extern ObitInfoList* HistoryGetList (ObitHistory* in) {
  return ObitInfoListRef(in->info);
}

extern int HistoryIsA (ObitHistory* in) {
  return ObitHistoryIsA(in);
}

extern char* HistoryGetName (ObitHistory* in) {
  if (ObitHistoryIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitHistory *me;
} History;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitHistory *me;
} History;

%addmethods History { 
  History(char *name, ObitInfoList *info, ObitErr *err) {
    History *out;
    out = (History *) malloc(sizeof(History));
    if (strcmp(name, "None")) out->me = newObitHistoryValue(name, info, err);
    else out->me = NULL;
    return out;
   }
  ~History() {
    self->me = HistoryUnref(self->me);
    free(self);
  }
};

/* $Id: ImageDesc.inc,v 1.9 2007/11/10 22:36:03 bcotton Exp $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageDesc type                         */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitImageDesc.h"
%}

%inline %{
extern ObitImageDesc* ImageDescCreate (char *name) {
  return newObitImageDesc (name);
} // end ImageDescCreate

extern ObitImageDesc* ImageDescCopy (ObitImageDesc* in, 
		              ObitImageDesc* out, ObitErr *err) {
  return ObitImageDescCopy (in, out, err);
} // end ImageDescCopy

extern void ImageDescCopyDesc (ObitImageDesc* in, ObitImageDesc* out,
			ObitErr *err) {
  ObitImageDescCopyDesc  (in, out, err);
} // end ImageDescCopyDesc

extern ObitImageDesc* ImageDescDefault (char *name) {
  return ObitImageDescDefault(name);
} // end ImageDescDefault

extern void ImageDescIndex (ObitImageDesc* in) {
  ObitImageDescIndex (in);
} // end ImageDescIndex

extern int ImageDescOverlap (ObitImageDesc* in1, 
		             ObitImageDesc* in2, ObitErr *err) {
  int tout;
  tout = ObitImageDescOverlap (in1, in2, err);
  return tout;
} // end ImageDescOverlap

// returns an array of 2 floats
extern PyObject* ImageDescCvtPixel (ObitImageDesc* in, ObitImageDesc* out, 
		      float *inPixel, ObitErr *err) {
  ofloat outPixel[7];
  PyObject *outList, *o;
  ObitImageDescCvtPixel (in, out, inPixel, outPixel, err);
  outList = PyList_New(2);
  o = PyFloat_FromDouble((double)outPixel[0]);
  PyList_SetItem(outList, 0, o);
  o = PyFloat_FromDouble((double)outPixel[1]);
  PyList_SetItem(outList, 1, o);	
  return outList;
} // end ImageDescCvtPixel

// returns an array of 2 floats
extern PyObject* ImageDescGetPos (ObitImageDesc* in, 
		      float *inPixel, ObitErr *err) {
  odouble outPos[2];
  PyObject *outList, *o;
  ObitImageDescGetPos (in, (ofloat*)inPixel, outPos, err);
  outList = PyList_New(2);
  o = PyFloat_FromDouble((double)outPos[0]);
  PyList_SetItem(outList, 0, o);
  o = PyFloat_FromDouble((double)outPos[1]);
  PyList_SetItem(outList, 1, o);	
  return outList;
} // end ImageDescGetPos

// returns an array of 2 floats
extern PyObject* ImageDescGetPixel (ObitImageDesc* in, 
		      double *inPos, ObitErr *err) {
  ofloat outPixel[7];
  PyObject *outList, *o;
  ObitImageDescGetPixel (in, (odouble*)inPos, outPixel, err);
  outList = PyList_New(2);
  o = PyFloat_FromDouble((double)outPixel[0]);
  PyList_SetItem(outList, 0, o);
  o = PyFloat_FromDouble((double)outPixel[1]);
  PyList_SetItem(outList, 1, o);	
  return outList;
} // end ImageDescGetPixel

extern ObitInfoList* ImageDescGetList (ObitImageDesc* in) {
  return ObitInfoListRef(in->info);
}
 
extern PyObject *ImageDescGetDict(ObitImageDesc* in) {
  PyObject *outDict = PyDict_New();
  PyObject *list1, *list2, *list3, *list4, *list5, *list6;
  PyObject *value;
  int i, pos = 0;

  PyDict_SetItemString(outDict, "name",    PyString_InternFromString(in->name));
  PyDict_SetItemString(outDict, "object",  PyString_InternFromString(in->object));
  PyDict_SetItemString(outDict, "teles",   PyString_InternFromString(in->teles));
  PyDict_SetItemString(outDict, "instrume",PyString_InternFromString(in->instrument));
  PyDict_SetItemString(outDict, "observer",PyString_InternFromString(in->observer));
  PyDict_SetItemString(outDict, "origin",  PyString_InternFromString(in->origin));
  PyDict_SetItemString(outDict, "date",    PyString_InternFromString(in->date));
  PyDict_SetItemString(outDict, "obsdat",  PyString_InternFromString(in->obsdat));
  PyDict_SetItemString(outDict, "bunit",   PyString_InternFromString(in->bunit));
  PyDict_SetItemString(outDict, "obsra",   PyFloat_FromDouble((double)in->obsra));
  PyDict_SetItemString(outDict, "obsdec",  PyFloat_FromDouble((double)in->obsdec));
  PyDict_SetItemString(outDict, "epoch",   PyFloat_FromDouble((double)in->epoch));
  PyDict_SetItemString(outDict, "equinox", PyFloat_FromDouble((double)in->equinox));
  PyDict_SetItemString(outDict, "beamMaj", PyFloat_FromDouble((double)in->beamMaj));
  PyDict_SetItemString(outDict, "beamMin", PyFloat_FromDouble((double)in->beamMin));
  PyDict_SetItemString(outDict, "beamPA",  PyFloat_FromDouble((double)in->beamPA));
  PyDict_SetItemString(outDict, "minval",  PyFloat_FromDouble((double)in->minval));
  PyDict_SetItemString(outDict, "maxval",  PyFloat_FromDouble((double)in->maxval));
  PyDict_SetItemString(outDict, "xshift",  PyFloat_FromDouble((double)in->xshift));
  PyDict_SetItemString(outDict, "yshift",  PyFloat_FromDouble((double)in->yshift));
  PyDict_SetItemString(outDict, "altCrpix",PyFloat_FromDouble((double)in->altCrpix));
  PyDict_SetItemString(outDict, "altRef",  PyFloat_FromDouble((double)in->altRef));
  PyDict_SetItemString(outDict, "restFreq",PyFloat_FromDouble((double)in->restFreq));
  PyDict_SetItemString(outDict, "areBlanks",PyInt_FromLong((long)in->areBlanks));
  PyDict_SetItemString(outDict, "niter",   PyInt_FromLong((long)in->niter));
  PyDict_SetItemString(outDict, "naxis",   PyInt_FromLong((long)in->naxis));
  PyDict_SetItemString(outDict, "bitpix",  PyInt_FromLong((long)in->bitpix));
  PyDict_SetItemString(outDict, "IOsize",  PyInt_FromLong((long)in->IOsize));
  PyDict_SetItemString(outDict, "VelReference", PyInt_FromLong((long)in->VelReference));
  PyDict_SetItemString(outDict, "VelDef", PyInt_FromLong((long)in->VelDef));
  list1 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list1, i, PyString_InternFromString(in->ctype[i]));
  PyDict_SetItemString(outDict, "ctype", list1);
  list2 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list2, i, PyFloat_FromDouble((double)in->crval[i]));
  PyDict_SetItemString(outDict, "crval", list2);
  list3 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list3, i, PyInt_FromLong((long)in->inaxes[i]));
  PyDict_SetItemString(outDict, "inaxes", list3);
  list4 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list4, i, PyFloat_FromDouble((double)in->cdelt[i]));
  PyDict_SetItemString(outDict, "cdelt", list4);
  list5 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list5, i, PyFloat_FromDouble((double)in->crpix[i]));
  PyDict_SetItemString(outDict, "crpix", list5);
  list6 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list6, i, PyFloat_FromDouble((double)in->crota[i]));
  PyDict_SetItemString(outDict, "crota", list6);
  /* Structural members - mostly read only */
  PyDict_SetItemString(outDict, "jlocr",    PyInt_FromLong((long)in->jlocr));
  PyDict_SetItemString(outDict, "jlocd",    PyInt_FromLong((long)in->jlocd));
  PyDict_SetItemString(outDict, "jlocs",    PyInt_FromLong((long)in->jlocs));
  PyDict_SetItemString(outDict, "jlocf",    PyInt_FromLong((long)in->jlocf));

  /* Discard references to newly created objects. */
  while (PyDict_Next(outDict, &pos, NULL, &value))
    Py_DECREF(value);

  return outDict;
} // end ImageDescGetDict

extern void ImageDescSetDict(ObitImageDesc* in, PyObject *inDict) {
  PyObject *list1, *list2, *list3, *list4, *list5, *list6;
  char *tstr;
  int i, number;

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return;
  }

  tstr = PyString_AsString(PyDict_GetItemString(inDict, "object"));
  strncpy (in->object, tstr, IMLEN_VALUE); in->object[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obsdat"));
  strncpy (in->obsdat, tstr, IMLEN_VALUE); in->obsdat[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "bunit"));
  strncpy (in->bunit, tstr, IMLEN_VALUE); in->bunit[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "teles"));
  strncpy (in->teles, tstr, IMLEN_VALUE); in->teles[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "instrume"));
  strncpy (in->instrument, tstr, IMLEN_VALUE); in->instrument[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "observer"));
  strncpy (in->observer, tstr, IMLEN_VALUE); in->observer[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "origin"));
  strncpy (in->origin, tstr, IMLEN_VALUE); in->origin[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "date"));
  strncpy (in->date, tstr, IMLEN_VALUE); in->date[IMLEN_VALUE-1]=0;
  in->epoch   = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "epoch"));
  in->equinox = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "equinox"));
  in->obsra   = PyFloat_AsDouble(PyDict_GetItemString(inDict, "obsra"));
  in->obsdec  = PyFloat_AsDouble(PyDict_GetItemString(inDict, "obsdec"));
  in->beamMaj = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "beamMaj"));
  in->beamMin = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "beamMin"));
  in->beamPA  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "beamPA"));
  in->restFreq= (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "restFreq"));
  in->minval  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "minval"));
  in->maxval  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "maxval"));
  in->xshift  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "xshift"));
  in->yshift  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "yshift"));
  in->altCrpix= (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "altCrpix"));
  in->altRef  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "altRef"));
  in->niter   = PyInt_AsLong(PyDict_GetItemString(inDict, "niter"));
  in->naxis   = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "naxis"));
  in->bitpix  = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "bitpix"));
  in->IOsize  = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "IOsize"));
  in->VelReference  = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "VelReference"));
  in->VelDef  = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "VelDef"));
  in->areBlanks = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "areBlanks"));
  list1 = PyDict_GetItemString(inDict, "ctype");
  number = MIN (IM_MAXDIM, PyList_Size(list1));
  for (i=0; i<number; i++) {
    tstr = PyString_AsString(PyList_GetItem(list1, i));
    strncpy (in->ctype[i], tstr, IMLEN_KEYWORD);
  }
  list2 = PyDict_GetItemString(inDict, "crval");
  number = MIN (IM_MAXDIM, PyList_Size(list2));
  for (i=0; i<number; i++) in->crval[i] = PyFloat_AsDouble(PyList_GetItem(list2, i));
  list3 = PyDict_GetItemString(inDict, "inaxes");
  number = MIN (IM_MAXDIM, PyList_Size(list3));
  for (i=0; i<number; i++) in->inaxes[i] = (long)PyInt_AsLong(PyList_GetItem(list3, i));
  list4 = PyDict_GetItemString(inDict, "cdelt");
  number = MIN (IM_MAXDIM, PyList_Size(list4));
  for (i=0; i<number; i++) in->cdelt[i] = (float)PyFloat_AsDouble(PyList_GetItem(list4, i));
  list5 = PyDict_GetItemString(inDict, "crpix");
  number = MIN (IM_MAXDIM, PyList_Size(list5));
  for (i=0; i<number; i++) in->crpix[i] = (float)PyFloat_AsDouble(PyList_GetItem(list5, i));
  list6 = PyDict_GetItemString(inDict, "crota");
  number = MIN (IM_MAXDIM, PyList_Size(list6));
  for (i=0; i<number; i++) in->crota[i] = (float)PyFloat_AsDouble(PyList_GetItem(list6, i));
  /* Reindex just to be sure */
  ObitImageDescIndex (in);
} // end ImageDescSetDict

ObitImageDesc* ImageDescRef (ObitImageDesc* in) {
  return ObitImageDescRef (in);
} // end ImageDescRef

ObitImageDesc* ImageDescUnref (ObitImageDesc* in) {
  if (!ObitImageDescIsA(in)) return NULL;
  return ObitImageDescUnref (in);
} // end ImageDescUnref

extern int ImageDescIsA (ObitImageDesc* in) {
  return ObitImageDescIsA(in);
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImageDesc *me;
} ImageDesc;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImageDesc *me;
} ImageDesc;

%addmethods ImageDesc { 
  ImageDesc(char *name) {
     ImageDesc *out;
     out = (ImageDesc *) malloc(sizeof(ImageDesc));
     if (strcmp(name, "None")) out->me = ImageDescCreate (name);
     else out->me = NULL;
     return out;
   }
  ~ImageDesc() {
    self->me = ImageDescUnref(self->me);
    free(self);
  }
};

/* $Id: ImageFit.inc,v 1.1 2007/09/17 15:58:00 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageFit type                          */
/*                                                                    */
/*;  Copyright (C) 2007                                               */
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
#include "ObitImageFit.h"
#include "ObitImage.h"
#include "ObitFitRegion.h"
%}


%inline %{
extern ObitImageFit* newImageFit (char* name) {
  return newObitImageFit (name);
} // end  newImageFit

extern ObitImageFit* ImageFitCopy  (ObitImageFit *in, ObitImageFit *out, 
				    ObitErr *err) {
  return ObitImageFitCopy (in, out, err);
} // end  ImageFitCopy

extern ObitImageFit* ImageFitUnref (ObitImageFit* in) {
  if (!ObitImageFitIsA(in)) return NULL;
  return ObitImageFitUnref(in);
}

extern ObitImageFit*  ImageFitRef (ObitImageFit* in) {
  return ObitImageFitRef(in);
}

extern ObitInfoList* ImageFitGetList (ObitImageFit* in) {
  return ObitInfoListRef(in->info);
}

extern int ImageFitFit (ObitImageFit* in, ObitImage *image, 
                        ObitFitRegion* reg, ObitErr *err) {
 return (int)ObitImageFitFit(in, image, reg, err);
}

extern char* ImageFitGetName (ObitImageFit* in) {
  if (ObitImageFitIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int ImageFitIsA (ObitImageFit* in) {
  return ObitImageFitIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImageFit *me;
} ImageFit;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImageFit *me;
} ImageFit;

%addmethods ImageFit { 
  ImageFit(char* name) {
     ImageFit *out;
     out = (ImageFit *) malloc(sizeof(ImageFit));
     if (strcmp(name, "None")) {
        out->me = newImageFit(name);
     } else out->me = NULL;
     return out;
   }
  ~ImageFit() {
   if (!self) return;  // Not defined
   if (self->me->ReferenceCount>0) 
      self->me = ImageFitUnref(self->me);
   free(self);
  }
};

/* $Id: Image.inc,v 1.14 2008/04/27 20:39:57 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Image type                             */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitImage.h"
#include "ObitData.h"
#include "ObitIOImageAIPS.h"
#include "ObitIOImageFITS.h"
%}


%inline %{

// size=1-> OBIT_IO_byRow, else OBIT_IO_byPlane
extern void ImageSetFITS(ObitImage *in, int size, int disk, char *file, 
                         int blc[7], int trc[7], ObitErr *err) {
  ObitIOSize fsize;

  if (size==1) fsize = OBIT_IO_byRow;
  else fsize = OBIT_IO_byPlane;
  ObitImageSetFITS(in, fsize, disk, file, blc, trc, err);
 }

extern void ImageSetAIPS(ObitImage *in, int size, int disk, int cno, int user, 
                         int blc[7], int trc[7], ObitErr *err) {
  ObitIOSize fsize;

  if (size==1) fsize = OBIT_IO_byRow;
  else fsize = OBIT_IO_byPlane;
  ObitImageSetAIPS(in, fsize, disk, cno, user, blc, trc, err);
 }

extern ObitData* ImageCastData (ObitImage* inImage) {
  return (ObitData*)inImage;
} // end ImageCastData

extern ObitImage* ImageCreate (char* name) {
  return newObitImage (name);
} // end  ImageCreate

extern ObitImage* ImageScratch (ObitImage *in, ObitErr *err) {
  return newObitImageScratch (in, err);
} // end  ImageScratch

extern PyObject* ImageInfo (ObitImage *in, ObitErr *err) {
  ObitIOImageAIPS *AIPSIO=NULL;
  ObitIOImageFITS *FITSIO=NULL;
  PyObject *outDict=PyDict_New();
  PyObject *o=NULL;

  if (err->error) return outDict;

  // Ensure in fully instantiated -assume OK if myIO exists 
  if (!in->myIO) ObitImageFullInstantiate (in, TRUE, err);
  if (err->error) return outDict;

  // Get details and save in dict
  if (ObitIOImageAIPSIsA(in->myIO)) {  // AIPS
    o = PyString_InternFromString("AIPS");
    PyDict_SetItemString(outDict, "type", o);
    AIPSIO = (ObitIOImageAIPS*)in->myIO;
    o = PyInt_FromLong((long)AIPSIO->disk);
    PyDict_SetItemString(outDict, "disk", o);
    o = PyInt_FromLong((long)AIPSIO->CNO);
    PyDict_SetItemString(outDict, "CNO", o);
    o = PyInt_FromLong((long)AIPSIO->UserId);
    PyDict_SetItemString(outDict, "user", o);
  } else if (ObitIOImageFITSIsA(in->myIO)) {  // FITS
    o = PyString_InternFromString("FITS");
    PyDict_SetItemString(outDict, "type", o);
    FITSIO = (ObitIOImageFITS*)in->myIO;
    o = PyInt_FromLong((long)FITSIO->disk);
    PyDict_SetItemString(outDict, "disk", o);
    o = PyString_InternFromString((char*)FITSIO->FileName);
    PyDict_SetItemString(outDict, "filename", o);
  } else {  // Don't know this one
    o = PyString_InternFromString("UNKNOWN");
    PyDict_SetItemString(outDict, "type", o);
  }
  return outDict;
} // end  ImageInfo

extern ObitImage* ImageZap  (ObitImage *in, ObitErr *err) {
  return ObitImageZap (in, err);
} // end ImageZap

extern void ImageRename  (ObitImage *in, ObitErr *err) {
  ObitImageRename (in, err);
} // end ImageRename

extern ObitImage* ImageCopy  (ObitImage *in, ObitImage *out, 
			         ObitErr *err) {
  return ObitImageCopy (in, out, err);
} // end  ImageCopy

extern void ImageClone (ObitImage *in, ObitImage *out, ObitErr *err) {
   return  ObitImageClone (in, out, err);
} // end  ImageClone

extern void ImageClone2 (ObitImage *in1, ObitImage *in2, ObitImage *out, ObitErr *err) {
   return  ObitImageClone2 (in1, in2, out, err);
} // end  ImageClone

extern void ImageCloneMem (ObitImage *in, ObitImage *out, ObitErr *err) {
   return  ObitImageCloneMem (in, out, err);
} // end  ImageCloneMem

// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
// Table verion returned as outValue1
extern ObitTable* ImageGetTable (ObitImage *in, int access, 
			    char *tabType, long *outValue1, ObitErr *err) {
  ObitIOAccess laccess;
  olong loutValue1 = (olong)*outValue1;
  ObitTable *outTable=NULL;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  outTable = newObitImageTable (in, laccess, tabType, &loutValue1, err);
  *outValue1 = (long)loutValue1;
  return outTable;
} // end  ImageGetTable

extern int ImageOpen (ObitImage *in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitImageOpen (in, laccess, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Open

// force header update 
extern void ImageDirty (ObitImage *in) {
  in->myStatus = OBIT_Modified;
} // end Dirty

extern int ImageClose (ObitImage *in, ObitErr *err) {
  ObitIOCode ret;
  ret =  ObitImageClose (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Close

extern int ImageRead (ObitImage *in, ObitErr *err) {
  ObitIOCode ret;
  in->extBuffer = FALSE;
  ret =  ObitImageRead (in, NULL, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Read

extern int ImageWrite (ObitImage *in, ObitErr *err) {
  ObitIOCode ret;
  in->extBuffer = FALSE;
  ret =  ObitImageWrite (in, NULL, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Write

extern int ImageReadFA (ObitImage *in, ObitFArray *array, 
                        ObitErr *err) {
  ObitIOCode ret;
  in->extBuffer = TRUE;
  /* Check comparability with an existing FArray */
  if (in->image) {
     if (!ObitFArrayIsCompatable(in->image, array)) {
       Obit_log_error(err, OBIT_InfoErr, "ImageReadFA: FArray incompatable");
       return 1;
    }
  }
  ret =  ObitImageRead (in, array->array, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end ReadFA

extern int ImageWriteFA (ObitImage *in, ObitFArray *array, 
                         ObitErr *err) {
  ObitIOCode ret;
  in->extBuffer = TRUE;
  /* Check comparability with an existing FArray */
  if (in->image) {
     if (!ObitFArrayIsCompatable(in->image, array)) {
       Obit_log_error(err, OBIT_InfoErr, "ImageWriteFA: FArray incompatable");
       return 1;
    }
  }
  ret =  ObitImageWrite (in, array->array, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Write

extern int ImageGetPlane (ObitImage *in, ObitFArray *outData, int *plane,
			    ObitErr *err) {
  ObitIOCode ret;
  olong i, lplane[5];
  gboolean bad;
  ofloat *data;

  for (i=0; i<5; i++) lplane[i] = plane[i];
  if (outData) data = outData->array;  // buffer given?
  else data = NULL;
  // Check compatibility between outData and in
  bad = FALSE;
  if (outData && (in->image)) bad = !ObitFArrayIsCompatable(outData,in->image);
  else if (outData && (in->myDesc->inaxes[0]>0)) { // check dimensions if possible
    bad = (outData->ndim<2) || (outData->naxis[0]!=in->myDesc->inaxes[0]) ||
          (outData->naxis[1]!=in->myDesc->inaxes[1]);
  }
  if (bad) {
    Obit_log_error(err, OBIT_Error, 
	           "ImageGetPlane: specified FArray incompatable with image");
    return 1;
  }
  ret = ObitImageGetPlane (in, data, plane, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end ImageGetPlane

extern int ImagePutPlane (ObitImage *in, ObitFArray *outData, int *plane,
			  ObitErr *err) {
  ObitIOCode ret;
  olong i, lplane[5];
  gboolean bad;
  ofloat *data;

  for (i=0; i<5; i++) lplane[i] = plane[i];
  if (outData) data = outData->array;  // buffer given?
  else data = NULL;
  // Check compatibility between outData and in
  bad = FALSE;
  if (outData && (in->image)) bad = !ObitFArrayIsCompatable(outData,in->image);
  else if (outData && (in->myDesc->inaxes[0]>0)) { // check dimensions if possible
    bad = (outData->ndim<2) || (outData->naxis[0]!=in->myDesc->inaxes[0]) ||
          (outData->naxis[1]!=in->myDesc->inaxes[1]);
  }
  if (bad) {
    Obit_log_error(err, OBIT_Error, 
	           "ImagePutPlane: specified FArray incompatable with image");
    return 1;
  }
  ret = ObitImagePutPlane (in, data, lplane, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end ImagePutPlane

extern int ImageZapTable (ObitImage *in, char *tabType, long tabVer, 
			    ObitErr *err) {
  ObitIOCode ret;
  ret = ObitImageZapTable (in, tabType, tabVer, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  ImageZapTable

extern int ImageCopyTables (ObitImage *in, ObitImage *out, char **exclude,
		  	        char **include, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitImageCopyTables  (in, out, exclude, include, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  CopyTables

extern int ImageUpdateTables (ObitImage *in, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitImageUpdateTables (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  UpdateTables

// Open and close to fully instantiate
// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
extern int ImagefullInstantiate (ObitImage* in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitImageOpen (in, laccess, err);
  ret = ObitImageClose (in, err);
  if ((err->error) || (ret!=OBIT_IO_OK)) return 1;
  else return 0;
} // end ImagefullInstantiate


extern ObitImage* ImageUnref (ObitImage* in) {
  if (!ObitImageIsA(in)) return NULL;
  if (in && (in->ReferenceCount>0)) in = ObitImageUnref(in);
  return in;
}

extern ObitImage*  ImageRef (ObitImage* in) {
  return ObitImageRef(in);
}

extern ObitInfoList* ImageGetList (ObitImage* in) {
  return ObitInfoListRef(in->info);
}

extern ObitTableList* ImageGetTableList (ObitImage* in) {
  return ObitTableListRef(in->tableList);
}

extern ObitImageDesc* ImageGetDesc (ObitImage* in) {
  return ObitImageDescRef(in->myDesc);
}

extern void ImageSetDesc (ObitImage* in, ObitImageDesc* desc) {
  in->myDesc = ObitImageDescUnref(in->myDesc);
  in->myDesc = ObitImageDescRef(desc);
}

extern ObitFArray* ImageGetFArray (ObitImage* in) {
  return ObitFArrayRef(in->image);
}

extern void ImageSetFArray (ObitImage* in, ObitFArray *image) {
  in->image = ObitFArrayUnref(in->image);
  in->image = ObitFArrayRef(image);
}

extern ObitImage* ImageGetBeam (ObitImage* in) {
  return (ObitImage*)ObitImageRef(in->myBeam);
}

extern void ImageSetBeam (ObitImage* in, ObitImage *beam) {
  in->myBeam = (Obit*)ObitImageUnref(in->myBeam);
  in->myBeam = (Obit*)ObitImageRef(beam);
}

extern long ImageGetHighVer (ObitImage* in, char *tabType) {
  return ObitTableListGetHigh(in->tableList, tabType);
}

extern int ImageisScratch (ObitImage* in) {
  return (int)in->isScratch;
}

extern int ImageIsA (ObitImage* in) {
  return ObitImageIsA(in);
}

extern char* ImageGetName (ObitImage* in) {
  if (ObitImageIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImage *me;
} Image;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImage *me;
} Image;

%addmethods Image { 
  Image(char* name) {
     Image *out;
     out = (Image *) malloc(sizeof(Image));
     if (strcmp(name, "None")) out->me = ImageCreate(name);
     else out->me = NULL;
     return out;
   }
  ~Image() { /* Scratch files may be deleted separately */
  if (self) {
    if (self->me->ReferenceCount>0) 
       self->me = ImageUnref(self->me);
    free(self);
    } 
  }
};

/* $Id: ImageMosaic.inc,v 1.4 2005/02/06 02:00:38 bcotton Exp $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageMosaic type                       */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitImageMosaic.h"
%}


%inline %{
extern ObitImageMosaic* newImageMosaic (char* name, int number) {
  return newObitImageMosaic (name, number);
} // end  newImageMosaic

extern ObitImageMosaic* ImageMosaicCopy  (ObitImageMosaic *in, ObitImageMosaic *out, 
				          ObitErr *err) {
  return ObitImageMosaicCopy (in, out, err);
} // end  ImageMosaicCopy

extern ObitImageMosaic* ImageMosaicUnref (ObitImageMosaic* in) {
  if (!ObitImageMosaicIsA(in)) return NULL;
  return ObitImageMosaicUnref(in);
}

extern ObitImageMosaic*  ImageMosaicRef (ObitImageMosaic* in) {
  return ObitImageMosaicRef(in);
}

extern void  ImageMosaicZapImage(ObitImageMosaic* in, int number, 
                                 ObitErr *err ) {
  return ObitImageMosaicZapImage(in, number, err);
}

extern ObitInfoList* ImageMosaicGetList (ObitImageMosaic* in) {
  return ObitInfoListRef(in->info);
}

extern ObitImage* ImageMosaicGetFullImage (ObitImageMosaic* in, ObitErr *err) {
  return ObitImageMosaicGetFullImage(in, err);
}

extern void ImageMosaicSetFullImage (ObitImageMosaic* in, 
                                 ObitImage *image, ObitErr *err) {
 ObitImageMosaicSetFullImage(in, image, err);
}

extern ObitImage* ImageMosaicGetImage (ObitImageMosaic* in, int number, 
                                       ObitErr *err) {
  return ObitImageMosaicGetImage(in, number, err);
}

extern void ImageMosaicSetImage (ObitImageMosaic* in, int number, 
                                 ObitImage *image, ObitErr *err) {
 ObitImageMosaicSetImage(in, number, image, err);
}

extern ObitImageMosaic* ImageMosaicCreate (char *name, ObitUV *uvData, ObitErr *err) {
 return ObitImageMosaicCreate(name, uvData, err);
}

extern void ImageMosaicDefine (ObitImageMosaic* in, ObitUV *uvData, int doBeam,
                               ObitErr *err) {
 gboolean LdoBeam = doBeam;
 ObitImageMosaicDefine(in, uvData, LdoBeam, err);
}

extern void ImageMosaicFlatten (ObitImageMosaic* in, ObitErr *err) {
 ObitImageMosaicFlatten(in, err);
}

extern char* ImageMosaicGetName (ObitImageMosaic* in) {
  if (ObitImageMosaicIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int ImageMosaicGetNumber (ObitImageMosaic* in) {
  if (ObitImageMosaicIsA(in)) {
    return in->numberImages;
  } else {
    return 0;
  }
}

extern int ImageMosaicIsA (ObitImageMosaic* in) {
  return ObitImageMosaicIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImageMosaic *me;
} ImageMosaic;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImageMosaic *me;
} ImageMosaic;

%addmethods ImageMosaic { 
  ImageMosaic(char* name, int number) {
     ImageMosaic *out;
     out = (ImageMosaic *) malloc(sizeof(ImageMosaic));
     if (strcmp(name, "None")) out->me = newImageMosaic(name, number);
     else out->me = NULL;
     return out;
   }
  ~ImageMosaic() {
   if (self->me-> ReferenceCount>0) 
      self->me = ImageMosaicUnref(self->me);
   free(self);
  }
};

/* $Id: ImageUtil.inc,v 1.4 2006/12/28 15:59:33 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Image utilities                        */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitImageUtil.h"
#include "ObitTable.h"
#include "ObitTableCCUtil.h"
%}



%inline %{
ObitImage* ImageUtilCreateImage (ObitUV *inUV, long fieldNo,
			       int doBeam, ObitErr *err) {
  return ObitImageUtilCreateImage (inUV, fieldNo, doBeam, err);
} // end ImageUtilCreateImag

void ImageUtilMakeImage (ObitUV *inUV, ObitImage *outImage, 
			     long channel, int doBeam, 
			     int doWeight, ObitErr *err) {
  ObitImageUtilMakeImage (inUV, outImage, channel, doBeam, doWeight, err);
} // end ImageUtilMakeImag

void 
ImageUtilInterpolateImage (ObitImage *inImage, ObitImage *outImage, 
			       int *inPlane, int *outPlane,
			       long hwidth, ObitErr *err)
{
  ObitImageUtilInterpolateImage (inImage, outImage, inPlane, outPlane,
			         hwidth, err);
} // end ImageUtilInterpolateImage

void 
ImageUtilPBApply (ObitImage *inImage, ObitImage *pntImage, ObitImage *outImage, 
                      int *inPlane, int *outPlane, float antSize, ObitErr *err)
{
  ObitImageUtilPBApply (inImage, pntImage, outImage, inPlane, outPlane, antSize, err);
} // end ImageUtilPBApply

void 
ImageUtilPBImage (ObitImage *pntImage, ObitImage *outImage, 
                  int *outPlane, float antSize, float minGain, ObitErr *err)
{
  ObitImageUtilPBImage (pntImage, outImage, outPlane, antSize, minGain, err);
} // end ImageUtilPBImage

void 
ImageUtilPBCorr (ObitImage *inImage, ObitImage *pntImage, ObitImage *outImage, 
                      int *inPlane, int *outPlane, float antSize, ObitErr *err)
{
  ObitImageUtilPBCorr (inImage, pntImage, outImage, inPlane, outPlane, antSize, err);
} // end ImageUtilPBCorr

ObitImage * 
ImageUtilQuanFITS(ObitImage *inImage, char *fileName, 
		  int disk, ObitErr *err)
{
  return ObitImageUtilQuanFITS(inImage, (gchar*)fileName, (olong)disk, err);
} // end ImageUtilPBCorr

void 
ImageUtilCCScale (ObitTable *inCCTab, int startComp, int endComp, double scale, 
	          ObitErr *err)
{
  olong lstartComp, lendComp;
  ofloat lscale;
  ObitTableCC *CCTab=NULL;

  if (err->error) return;  // error check 

  lstartComp = startComp;
  lendComp = endComp;
  lscale   = (ofloat)scale;
  CCTab = ObitTableCCConvert(inCCTab);
  ObitTableCCUtilScale (CCTab, lstartComp, lendComp, lscale, err);
  CCTab = ObitTableCCUnref(CCTab);
} // end ImageUtilCCScale

%}
/* $Id: InfoList.inc,v 1.11 2008/02/03 23:06:51 bcotton Exp $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitInfoList type                      */
/*                                                                    */
/*   Copyright (C) 2004-2008                                          */
/*   Associated Universities, Inc. Washington DC, USA.                */
/*                                                                    */
/*   This program is free software; you can redistribute it and/or    */
/*   modify it under the terms of the GNU General Public License as   */
/*   published by the Free Software Foundation; either version 2 of   */
/*   the License, or (at your option) any later version.              */
/*                                                                    */
/*   This program is distributed in the hope that it will be useful,  */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of   */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    */
/*   GNU General Public License for more details.                     */
/*                                                                    */
/*   You should have received a copy of the GNU General Public        */
/*   License along with this program; if not, write to the Free       */
/*   Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,     */
/*   MA 02139, USA.                                                   */
/*                                                                    */
/* Correspondence about this software should be addressed as follows: */
/*          Internet email: bcotton@nrao.edu.                         */
/*          Postal address: William Cotton                            */
/*                          National Radio Astronomy Observatory      */
/*                          520 Edgemont Road                         */
/*                          Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{
#include "ObitInfoList.h"
%}

// This cleans up the InfoListBlob we malloc'd before the function call
%typemap(python,freearg) InfoListBlob * {
   InfoListBlob *me = (InfoListBlob *)$source;
   if (me!=NULL) {
      if (me->name) free(me->name); me->name = NULL;
      if (me->data) free(me->data); me->data = NULL;
      free(me);
   }
}

// This tells SWIG to treat a InfoListBlob * argument with name 'outBlob' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.
// A type of -999 indicates failure

%typemap(python,argout) InfoListBlob *outBlob {
  PyObject *o, *t;
  int i, j, count;
  char tstring[1000];

 // Check for failure
 if ($source->type==-999) {
     o = PyString_FromString("Item not Found");
     PyList_Append($target,o);
     Py_XDECREF(o);
     return NULL;
 }
 o = PyString_FromString($source->name);
  if ((!$target) || ($target == Py_None)) {
      $target = PyList_New(0);
      PyList_Append($target,o);
      Py_XDECREF(o);
  } else {
    if (!PyList_Check($target)) {
      PyObject *o2 = $target;
      $target = PyList_New(0);
      PyList_Append($target,o2);
      Py_XDECREF(o2);
    }
    PyList_Append($target,o);
    Py_XDECREF(o);
  }

  // Type
    o = PyInt_FromLong((long)$source->type);
    PyList_Append($target,o);
    Py_XDECREF(o);

  // Dim (5)
    t = PyList_New(0);
    count = 1;
    for (i=0; i<5; i++) {
      if ($source->dim[i]>0) count *= $source->dim[i];
      o = PyInt_FromLong((long)$source->dim[i]);
      PyList_Append(t, o);
      Py_XDECREF(o);
    }
    PyList_Append($target,t);
    Py_XDECREF(t);

    if ($source->type == OBIT_string)
      count /= $source->dim[0];  /* now number of strings rather than total */

   // Data
    t = PyList_New(0);
    for (i=0; i<count; i++) {
      // a few for now 
      switch ($source->type) {
      case   OBIT_int: 
	o = PyInt_FromLong((long)((int*)$source->data)[i]);
	break;
      case   OBIT_oint:
	o = PyInt_FromLong((long)((oint*)$source->data)[i]);
        break;
      case   OBIT_long:  
	o = PyInt_FromLong(((long*)$source->data)[i]);
        break;
      case   OBIT_float: 
	o = PyFloat_FromDouble((double)((float*)$source->data)[i]);
        break;
      case   OBIT_string:
        for (j=0; j<$source->dim[0]; j++) 
	tstring[j] = ((char*)$source->data)[i*$source->dim[0]+j];
	tstring[j] = 0;
	o = PyString_FromString(tstring);
        break;
      case   OBIT_bool:
	o = PyInt_FromLong((long)((int*)$source->data)[i]);
        break;
      case   OBIT_double:
	o = PyFloat_FromDouble(((double*)$source->data)[i]);
        break;
      default:
        PyErr_SetString(PyExc_TypeError,"Unknown InfoList data type");
        return NULL;
      }; // end switch
      PyList_Append(t, o);
      Py_XDECREF(o);
    }
    PyList_Append($target,t);
    Py_XDECREF(t);
}



// Functions not implemented
//void ObitInfoListPut(ObitInfoList *in, 
//		      char* name, ObitInfoType type, gint32 *dim, 
//		      gconstpointer data, ObitErr *err);
//void ObitInfoListAlwaysPut(ObitInfoList *in, 
//		      char* name, ObitInfoType type, gint32 *dim, 
//		      gconstpointer data);
//int ObitInfoListInfo(ObitInfoList *in, 
//		     char *name, ObitInfoType *type, gint32 *dim, 
//		      ObitErr *err);
//int ObitInfoListGet(ObitInfoList *in, 
//		      char *name, ObitInfoType *type, gint32 *dim, 
//		      gpointer data, ObitErr *err);
//int ObitInfoListGetP(ObitInfoList *in, 
//		     char *name, ObitInfoType *type, gint32 *dim, 
//		      pointer *data);
//int ObitInfoListGetTest(ObitInfoList *in, 
//		      char *name, ObitInfoType *type, gint32 *dim, 
//		      pointer data);
//int ObitInfoListGetNumber(ObitInfoList *in, olong number,
//		      char *name, ObitInfoType *type, gint32 *dim, 
//		      pointer data, ObitErr *err);
//int ObitInfoListGetNumberP(ObitInfoList *in, olong number,
//		       char *name, ObitInfoType *type, gint32 *dim, 
//		       pointer *data);


%inline %{
// Structure for InfoList entry
typedef struct {
  char *name;
  int  type;
  long dim[5];
  void *data;
} InfoListBlob;

extern InfoListBlob* makeInfoListBlob() {
   InfoListBlob *out;
   out = malloc(sizeof(InfoListBlob));
   out->name = NULL;
   out->data = NULL;
   return out;
}

extern void freeInfoListBlob(InfoListBlob* blob) {
   if (!blob) return;
   if (blob->name) free(blob->name);
   if (blob->data) free(blob->data);
   free (blob);
}

ObitInfoList* InfoListCreate (void) {
  return newObitInfoList();
} // end InfoListCreate

ObitInfoList* freeInfoList (ObitInfoList *in) {
   return freeObitInfoList(in);
} // end freeList

ObitInfoList* InfoListCopy (ObitInfoList* in) {
  return ObitInfoListCopy (in);
} // end InfoListCopy

ObitInfoList* InfoListRef (ObitInfoList* in) {
  return ObitInfoListRef (in);
} // end InfoListRef

ObitInfoList* InfoListUnref (ObitInfoList* in) {
  if (!ObitInfoListIsA(in)) return NULL;
  return ObitInfoListUnref (in);
} // end InfoListUnref

ObitInfoList* InfoListCopyData (ObitInfoList* in, ObitInfoList* out) {
  return ObitInfoListCopyData (in, out);
} // end InfoListCopyData

void InfoListRemove (ObitInfoList *in, char *name) {
  return ObitInfoListRemove (in, name);
} // end InfoListRemove

void InfoListItemResize(ObitInfoList *in, 
			    char *name, int type, long *dim) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
  ObitInfoListResize (in, name, type, ddim);
} // end InfoListItemResize

int InfoListIsA (ObitInfoList* in) {
  return  ObitInfoListIsA(in);
} // end InfoListListIsA


// Put by type
extern void InfoListPutInt(ObitInfoList *in, char *name, long *dim, 
	             int* data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_long, ddim, (gpointer)data, err);
} // end InfoListPutInt

extern void InfoListAlwaysPutInt(ObitInfoList *in, char *name, 
	             long *dim, int* data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_long, ddim, (gpointer)data);
} // end InfoListAlwaysPutInp

extern void InfoListPutLong(ObitInfoList *in, char *name, long *dim, 
	             long* data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_long, ddim, (gpointer)data, err);
} // end InfoListPutLong

extern void InfoListAlwaysPutLong(ObitInfoList *in, char *name, 
	             long *dim, long* data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_long, ddim, (gpointer)data);
} // end InfoListAlwaysPutLong

extern void InfoListPutFloat(ObitInfoList *in, char *name, long *dim, 
	             float* data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_float, ddim, (gpointer)data, err);
} // end InfoListPutFloat

extern void InfoListAlwaysPutFloat(ObitInfoList *in, char *name, 
	             long *dim, float* data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_float, ddim, (gpointer)data);
} // end InfoListAlwaysPutFloat

extern void InfoListPutDouble(ObitInfoList *in, char *name, long *dim, 
	             double* data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_double, ddim, (gpointer)data, err);
} // end InfoListPutDouble

extern void InfoListAlwaysPutDouble(ObitInfoList *in, char *name, 
	             long *dim, double* data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_double, ddim, (gpointer)data);
} // end InfoListAlwaysPutDouble

extern void InfoListPutBoolean(ObitInfoList *in, char *name, long *dim, 
	             int* data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_bool, ddim, (gpointer)data, err);
} // end InfoListPutBoolean

extern void InfoListAlwaysPutBoolean(ObitInfoList *in, char *name, 
	             long *dim, int* data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_bool, ddim, (gpointer)data);
} // end InfoListAlwaysPutBoolean

extern void InfoListPutString(ObitInfoList *in, char *name, long *dim, 
	             char** data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_string, ddim, (gpointer)data, err);
} // end InfoListPutString

extern void InfoListAlwaysPutString(ObitInfoList *in, char *name, 
	             long *dim, char** data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_string, ddim, (gpointer)data);
} // end InfoListAlwayPutString

// single string version 
extern void InfoListPutSString(ObitInfoList *in, char *name, long *dim, 
	             char** data, ObitErr *err) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListPut(in, name, OBIT_string, ddim, data[0], err);
} // end InfoListPutSString

// single string version 
extern void InfoListAlwaysPutSString(ObitInfoList *in, char *name, 
	             long *dim, char** data) {
   gint32 i, ddim[5];

   for (i=0; i<5; i++) ddim[i] = dim[i];
   ObitInfoListAlwaysPut(in, name, OBIT_string, ddim, data[0]);
} // end InfoListAlwayPutsString

static int InfoListGetHelper(ObitInfoList *in, char *name, ObitInfoType type,
	gint32 dim[5], void *data, InfoListBlob *outBlob)
{
  gint32 i, j, count;
  float  *fdata;
  long   *ldata;
  int    *idata;
  oint   *odata;
  double *ddata;
  gboolean *bdata;
  char   *cdata, *codata;

  outBlob->name = strdup(name);
  outBlob->type = type;

  count = 1;
  for (i=0; i<5; i++) {
    outBlob->dim[i]  = dim[i];
    if (dim[i]>0) count *= dim[i];
  } 
  switch (type) {
    case   OBIT_int: 
      idata = (int*)data;
      outBlob->data = (void*)malloc(count*sizeof(int));
      for (i=0; i<count; i++) ((int*)outBlob->data)[i] = idata[i];
      break;
    case   OBIT_oint:
      odata = (oint*)data;
      outBlob->data = (void*)malloc(count*sizeof(oint));
      for (i=0; i<count; i++) ((oint*)outBlob->data)[i] = odata[i];
      break;
    case   OBIT_long:  
      ldata = (long*)data;
      outBlob->data = (void*)malloc(count*sizeof(long));
      for (i=0; i<count; i++) ((long*)outBlob->data)[i] = ldata[i];
      break;
    case   OBIT_float: 
      fdata = (float*)data;
      outBlob->data = (void*)malloc(count*sizeof(float));
      for (i=0; i<count; i++) ((float*)outBlob->data)[i] = fdata[i];
      break;
    case   OBIT_double:
      ddata = (double*)data;
      outBlob->data = (void*)malloc(count*sizeof(double));
      for (i=0; i<count; i++) ((double*)outBlob->data)[i] = ddata[i];
      break;
    case   OBIT_bool:
      bdata = (gboolean*)data;
      outBlob->data = (void*)malloc(count*sizeof(gboolean));
      for (i=0; i<count; i++) ((gboolean*)outBlob->data)[i] = bdata[i];
      break;
    case   OBIT_string:
      cdata = (char*)data;
      outBlob->data = (void*)malloc(count+1);
      codata = (char*)outBlob->data;
      count /= dim[0];  /* now number of strings rather than total */
      for (i=0; i<count; i++) {  // add nulls at end of strings
	for (j=0; j<dim[0]; j++)  *(codata++) = *(cdata++);
	*(codata++) = 0;
      }
      break;
    default:
      PyErr_SetString(PyExc_TypeError,"Unknown InfoList data type");
      return 0;
  }; // end switch

  return 0;
} // end InfoListGetHelper

// Convert InfoList data to a python list
static PyObject* Info2List(ObitInfoType type, gint32 dim[5], void *data)
{
  PyObject *outList=NULL, *outList2=NULL, *outList3=NULL, *o=NULL;
  gint32 i, j, k, ii, count[2], ndim;
  float  *fdata;
  long   *ldata;
  int    *idata;
  oint   *odata;
  double *ddata;
  gboolean *bdata;
  char   *tstring=NULL, *cdata;

  // how much data? - do three dimensions for now
  count[0] = dim[0]; ndim = 1;
  count[1] = dim[1]; if (dim[1]>1) ndim=2;
  count[2] = dim[2]; if (dim[2]>1) ndim=3;

  outList = PyList_New(count[0]);  // Output list
  ii = 0;

  switch (type) {
    case   OBIT_int: 
      outList = PyList_New(count[0]);  // Output list
      idata = (int*)data;
      for (i=0; i<count[0]; i++) {
  	if (ndim==1) {  // single dimension
          o = PyInt_FromLong((long)idata[ii++]);
          PyList_SetItem(outList, i, o);
        } else { // 2 or more dimensions
          outList2 = PyList_New(count[1]);
          for (j=0; j<count[1]; j++) {  // loop over second dimension
            if (ndim==2) {  // two dimensions
              o = PyInt_FromLong((long)idata[ii++]);
              PyList_SetItem(outList2, j, o);
            } else { // 3 (or more) dimensions
              outList3 = PyList_New(count[2]);
              for (k=0; k<count[2]; k++) {  // loop over third dimension
                o = PyInt_FromLong((long)idata[ii++]);
                PyList_SetItem(outList3, k, o);
              } // end loop over third dimension
              PyList_SetItem(outList2, j, outList3);
            }  
          }  // end loop over second dimension
          PyList_SetItem(outList, i, outList2);
        }
      } // end loop over first dimension
      break;
    case   OBIT_oint:
      outList = PyList_New(count[0]);  // Output list
      odata = (oint*)data;
      for (i=0; i<count[0]; i++) {
  	if (ndim==1) {  // single dimension
          o = PyInt_FromLong((long)odata[ii++]);
          PyList_SetItem(outList, i, o);
        } else { // 2 or more dimensions
          outList2 = PyList_New(count[1]);
          for (j=0; j<count[1]; j++) {  // loop over second dimension
            if (ndim==2) {  // two dimensions
              o = PyInt_FromLong((long)odata[ii++]);
              PyList_SetItem(outList2, j, o);
            } else { // 3 (or more) dimensions
              outList3 = PyList_New(count[2]);
              for (k=0; k<count[2]; k++) {  // loop over third dimension
                o = PyInt_FromLong((long)odata[ii++]);
                PyList_SetItem(outList3, k, o);
              } // end loop over third dimension
              PyList_SetItem(outList2, j, outList3);
            }  
          }  // end loop over second dimension
          PyList_SetItem(outList, i, outList2);
        }
      } // end loop over first dimension
      break;
    case   OBIT_long:  
      outList = PyList_New(count[0]);  // Output list
      ldata = (long*)data;
      for (i=0; i<count[0]; i++) {
  	if (ndim==1) {  // single dimension
          o = PyInt_FromLong((long)ldata[ii++]);
          PyList_SetItem(outList, i, o);
        } else { // 2 or more dimensions
          outList2 = PyList_New(count[1]);
          for (j=0; j<count[1]; j++) {  // loop over second dimension
            if (ndim==2) {  // two dimensions
              o = PyInt_FromLong((long)ldata[ii++]);
              PyList_SetItem(outList2, j, o);
            } else { // 3 (or more) dimensions
              outList3 = PyList_New(count[2]);
              for (k=0; k<count[2]; k++) {  // loop over third dimension
                o = PyInt_FromLong((long)ldata[ii++]);
                PyList_SetItem(outList3, k, o);
              } // end loop over third dimension
              PyList_SetItem(outList2, j, outList3);
            }  
          }  // end loop over second dimension
          PyList_SetItem(outList, i, outList2);
        }
      } // end loop over first dimension
      break;
    case   OBIT_float: 
      outList = PyList_New(count[0]);  // Output list
      fdata = (float*)data;
      for (i=0; i<count[0]; i++) {
  	if (ndim==1) {  // single dimension
          o = PyFloat_FromDouble((double)fdata[ii++]);
          PyList_SetItem(outList, i, o);
        } else { // 2 or more dimensions
          outList2 = PyList_New(count[1]);
          for (j=0; j<count[1]; j++) {  // loop over second dimension
            if (ndim==2) {  // two dimensions
              o = PyFloat_FromDouble((double)fdata[ii++]);
              PyList_SetItem(outList2, j, o);
            } else { // 3 (or more) dimensions
              outList3 = PyList_New(count[2]);
              for (k=0; k<count[2]; k++) {  // loop over third dimension
                o = PyFloat_FromDouble((double)fdata[ii++]);
                PyList_SetItem(outList3, k, o);
              } // end loop over third dimension
              PyList_SetItem(outList2, j, outList3);
            }  
          }  // end loop over second dimension
          PyList_SetItem(outList, i, outList2);
        }
      } // end loop over first dimension
      break;
    case   OBIT_double:
      outList = PyList_New(count[0]);  // Output list
      ddata = (double*)data;
      for (i=0; i<count[0]; i++) {
  	if (ndim==1) {  // single dimension
          o = PyFloat_FromDouble((double)ddata[ii++]);
          PyList_SetItem(outList, i, o);
        } else { // 2 or more dimensions
          outList2 = PyList_New(count[1]);
          for (j=0; j<count[1]; j++) {  // loop over second dimension
            if (ndim==2) {  // two dimensions
              o = PyFloat_FromDouble((double)ddata[ii++]);
              PyList_SetItem(outList2, j, o);
            } else { // 3 (or more) dimensions
              outList3 = PyList_New(count[2]);
              for (k=0; k<count[2]; k++) {  // loop over third dimension
                o = PyFloat_FromDouble((double)ddata[ii++]);
                PyList_SetItem(outList3, k, o);
              } // end loop over third dimension
              PyList_SetItem(outList2, j, outList3);
            }  
          }  // end loop over second dimension
          PyList_SetItem(outList, i, outList2);
        }
      } // end loop over first dimension
      break;
    case   OBIT_bool:
      outList = PyList_New(count[0]);  // Output list
      bdata = (gboolean*)data;
      for (i=0; i<count[0]; i++) {
  	if (ndim==1) {  // single dimension
          o = PyInt_FromLong((long)bdata[ii++]);
          PyList_SetItem(outList, i, o);
        } else { // 2 or more dimensions
          outList2 = PyList_New(count[1]);
          for (j=0; j<count[1]; j++) {  // loop over second dimension
            if (ndim==2) {  // two dimensions
              o = PyInt_FromLong((long)bdata[ii++]);
              PyList_SetItem(outList2, j, o);
            } else { // 3 (or more) dimensions
              outList3 = PyList_New(count[2]);
              for (k=0; k<count[2]; k++) {  // loop over third dimension
                o = PyInt_FromLong((long)bdata[ii++]);
                PyList_SetItem(outList3, k, o);
              } // end loop over third dimension
              PyList_SetItem(outList2, j, outList3);
            }  
          }  // end loop over second dimension
          PyList_SetItem(outList, i, outList2);
        }
      } // end loop over first dimension
      break;
    case   OBIT_string:
      cdata = (char*)data;
      outList = PyList_New(MAX (1, count[1]));  // Output list
      // First dimension is the length of strings
      tstring	= g_malloc0(count[0]+1);
      if (ndim==1) {  // single dimension - one string
        for (i=0; i<count[0]; i++) tstring[i]=cdata[ii++]; tstring[i]=0;
        o = PyString_InternFromString(tstring);
        PyList_SetItem(outList, 0, o);
      } else { // 2 or more dimensions
        for (j=0; j<count[1]; j++) {  // loop over second dimension
          if (ndim==2) {  // two dimensions
            for (i=0; i<count[0]; i++) tstring[i]=cdata[ii++]; tstring[i]=0;
            o = PyString_InternFromString(tstring);
            PyList_SetItem(outList, j, o);
          } else { // 3 (or more) dimensions
            outList2 = PyList_New(count[2]);
            for (k=0; k<count[2]; k++) {  // loop over third dimension
              for (i=0; i<count[0]; i++) tstring[i]=cdata[ii++]; tstring[i]=0;
              o = PyString_InternFromString(tstring);
              PyList_SetItem(outList2, k, o);
            } // end loop over third dimension
            PyList_SetItem(outList, j, outList2);
          }  
        }  // end loop over second dimension
      } // end multiple dimensions
      if (tstring) g_free(tstring);
      break;
    default:
      PyErr_SetString(PyExc_TypeError,"Unknown InfoList data type");
      return outList;
  }; // end switch

  return outList;
} // end Info2list

// Retrieve value as python list
//       0 - return code, 0=OK else failed
//       1 - name
//       2 - type
//       3 - dimension array
//       4 - data array
extern PyObject* InfoListGet(ObitInfoList *in, char *name) {
  gint32 i, dim[5];
  ObitInfoType type;
  void  *data;
  PyObject *dlist=NULL;
  PyObject *outList=NULL, *tl, *o;

  // output list
  outList = PyList_New(5);

  if (ObitInfoListGetP(in, name, &type, dim, &data)) {
    dlist = Info2List(type, dim, data);
    o = PyInt_FromLong((long)0);
    PyList_SetItem(outList, 0, o);
    o = PyString_InternFromString(name);
    PyList_SetItem(outList, 1, o);	
    o = PyInt_FromLong((long)type);
    PyList_SetItem(outList, 2, o);
    tl = PyList_New(5);
    for (i=0; i<5; i++) {
      o = PyInt_FromLong((long)dim[i]);
      PyList_SetItem(tl, i, o);
    }	
    PyList_SetItem(outList, 3, tl);	
    PyList_SetItem(outList, 4, dlist);	
    return outList;
   } /* end worked OK */

  // failed 
  //PyErr_SetString(PyExc_TypeError,"Item not found");
  o = PyInt_FromLong((long)999);
  PyList_SetItem(outList, 0, o);
  o = PyString_InternFromString(name);
  PyList_SetItem(outList, 1, o);	
  PyList_SetItem(outList, 2, NULL);	
  PyList_SetItem(outList, 3, NULL);	
  PyList_SetItem(outList, 4, NULL);	
  return outList;
} // end InfoListGet

// Retrieve numbered value as python list
//       0 - return code, 0=OK else failed
//       1 - name
//       2 - type
//       3 - dimension array
//       4 - data array
extern PyObject* InfoListGetNumber(ObitInfoList *in, int number)
{
  gint32 i, dim[5];
  char *name;
  ObitInfoType type;
  void  *data;
  PyObject *dlist=NULL;
  PyObject *outList=NULL, *tl, *o;

  // output list
  outList = PyList_New(5);

  if (ObitInfoListGetNumberP(in, number, &name, &type, dim, &data)) {
    dlist = Info2List(type, dim, data);
    o = PyInt_FromLong((long)0);
    PyList_SetItem(outList, 0, o);
    o = PyString_InternFromString(name);
    PyList_SetItem(outList, 1, o);	
    o = PyInt_FromLong((long)type);
    PyList_SetItem(outList, 2, o);
    tl = PyList_New(5);
    for (i=0; i<5; i++) {
      o = PyInt_FromLong((long)dim[i]);
      PyList_SetItem(tl, i, o);
    }	
    PyList_SetItem(outList, 3, tl);	
    PyList_SetItem(outList, 4, dlist);	
    return outList;
   } /* end worked OK */

  // failed 
  PyErr_SetString(PyExc_TypeError,"Item not found");
  o = PyInt_FromLong((long)999);
  PyList_SetItem(outList, 0, o);
  o = PyString_InternFromString(name);
  PyList_SetItem(outList, 1, o);	
  PyList_SetItem(outList, 2, NULL);	
  PyList_SetItem(outList, 3, NULL);	
  PyList_SetItem(outList, 4, NULL);	
  return outList;
} // end InfoListGetNumber

// Retrieve whole infolist as a python dict
// each entry as python list
//       0 - type
//       1 - dimension array
//       2 - data array
extern PyObject* InfoListGetDict(ObitInfoList *in)
{
  PyObject *outDict=PyDict_New();
  gint32 i, iel, dim[5];
  char *name;
  ObitInfoType type;
  void  *data;
  PyObject *dlist=NULL;
  PyObject *tList=NULL, *tl, *o;

  // loop over elements
  for (iel = 1; iel<=in->number; iel++) {
    // output list
    tList = PyList_New(3);
  
    if (ObitInfoListGetNumberP(in, iel, &name, &type, dim, &data)) {
      dlist = Info2List(type, dim, data);
      o = PyInt_FromLong((long)type);
      PyList_SetItem(tList, 0, o);
      tl = PyList_New(5);
      for (i=0; i<5; i++) {
        o = PyInt_FromLong((long)dim[i]);
        PyList_SetItem(tl, i, o);
      }	
      PyList_SetItem(tList, 1, tl);	
      PyList_SetItem(tList, 2, dlist);	
      // add to output
      PyDict_SetItemString(outDict, name, tList);
    } /* end worked OK */
  }	  // end element loop
  return outDict;
} // end InfoListGetDict

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitInfoList *me;
} InfoList;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitInfoList *me;
} InfoList;

%addmethods InfoList { 
  InfoList(void) {
     InfoList *out;
     out = (InfoList *) malloc(sizeof(InfoList));
     out->me = InfoListCreate();
    return out;
   }
  ~InfoList() {
    self->me = InfoListUnref(self->me);
    free(self);
  }
};

/* $Id: IonCal.inc,v 1.3 2007/09/07 12:32:42 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Ionospheric calibration utilities      */
/*                                                                    */
/*;  Copyright (C) 2007-2008                                          */
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
#include "ObitIoN2SolNTable.h"
%}

%inline %{
ObitTable* IoN2SolNTableConvert (ObitUV *inUV, long outSNVer, ObitTable *NITable,
                                 double pos[2], ObitErr *err) {
    ObitTableNI *lNITable = NULL;
    ObitTableSN *loutSN   = NULL;
    ofloat lpos[2];
    olong  loutSNVer = (olong)outSNVer;
    oint numPol=1, numIF=1;

    lNITable  = ObitTableNIConvert(NITable);
    if (inUV->myDesc->jlocif>=0) numIF = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
    if (inUV->myDesc->jlocs>=0) numPol = inUV->myDesc->inaxes[inUV->myDesc->jlocs];

    /* Pre existing SN table? */
    if (outSNVer>0) 
        loutSN = newObitTableSNValue ("Ion SN table", (ObitData*)inUV, &loutSNVer,
   			              OBIT_IO_WriteOnly, numPol, numIF, err);

    lpos[0] = (ofloat)pos[0];
    lpos[1] = (ofloat)pos[1];
    loutSN   = ObitIoN2SolNTableConvert (inUV, lNITable, loutSN, lpos, err);
    lNITable = ObitTableNIUnref(lNITable);
    return (ObitTable*)loutSN;
} // end IoN2SolNTableConvert 


%}
/* $Id: ObitErr.inc,v 1.5 2005/12/08 02:01:20 bcotton Exp $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitErr type                           */
/*                                                                    */
/*;  Copyright (C) 2004                                               */
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
#include "ObitErr.h"
%}

extern ObitErr* ObitErrUnref (ObitErr* in);
extern ObitErr*  ObitErrRef (ObitErr* in);
extern void ObitErrLog (ObitErr* in);
extern void ObitErrClear (ObitErr* in);
extern int ObitErrIsA (ObitErr* in);

%inline %{

extern ObitErr* ObitErrCreate (void) {
  return newObitErr();
}

extern int isError(ObitErr *in) {
   return (int)in->error;
}

extern void SetError(ObitErr *in) {
   in->error = 1;
}

extern void LogError (ObitErr *in, int eCode, char *message) {
   ObitErrCode code;
 /* Should be coordinated with OErr class definition.*/

   switch (eCode) {
   case 0:
     code = OBIT_None;
     break;
   case 1:
     code = OBIT_InfoErr;
     break;
   case 2:
     code = OBIT_InfoWarn;
     break;
   case 3:
     code = OBIT_Traceback;
     break;
   case 4:
     code = OBIT_MildError;
     break;
   case 5:
     code = OBIT_Error;
     break;
   case 6:
     code = OBIT_StrongError;
     break;
   case 7:
     code = OBIT_Fatal;
     break;
   case 8:
     code = OBIT_None;
     break;
   default:
     code = OBIT_None;
   };  

   ObitErrPush (in, code, message);
}  //  end LogError


extern char *OErrMsg(ObitErr *in)
{
  ObitErrCode errLevel;
  gchar *errMsg, *errLevelStr;

/*
 * Human readable versions of the ObitErr codes.
 * Should be coordinated with enum definition.
 */
gchar *ObitErrorLevelString[] = {
  "no message   ", /* OBIT_None        */
  "information  ", /* OBIT_InfoErr     */
  "warning      ", /* OBIT_InfoWarn    */
  "traceback    ", /* OBIT_Traceback   */
  "Mild error   ", /* OBIT_MildError   */
  "Error        ", /* OBIT_Error       */
  "Serious error", /* OBIT_StrongError */
  "Fatal error  ", /* OBIT_Fatal       */
};

  /* error checks */
  g_assert (ObitErrIsA(in));

  ObitErrPop (in, &errLevel, &errMsg);
  if(errMsg) {
    gchar *str;
    /* convert error level to something human readable */
    errLevelStr = ObitErrorLevelString[errLevel];
    str = g_strdup_printf("%s: %s", errLevelStr, errMsg);
    g_free(errMsg);
    errMsg = str;
  }
  return errMsg;
}

/* Force an abort as a really heavyhanded interrupt */
extern void Bomb(void) {
  char ct, *PythonSux = NULL;
  ct = PythonSux[-1000000000];
  /* It that doesn't work - maybe this will */
  abort();
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitErr *me;
} OErr;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitErr *me;
} OErr;

%addmethods OErr { 
  OErr(void) {
     OErr *out;
     out = (OErr *) malloc(sizeof(OErr));
     out->me = ObitErrCreate();
    return out;
   }
  ~OErr() {
    self->me = ObitErrUnref(self->me);
    free(self);
  }
};

/* $Id: ObitSystem.inc,v 1.9 2007/07/26 14:28:42 bcotton Exp $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitSystem type                        */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include <stdio.h>
#include "ObitSystem.h"
#include "ObitMem.h"
%}

%inline %{

extern ObitSystem* Startup (char *pgmName, int pgmNumber, int AIPSuser,
		            int numberAIPSdisk, char* AIPSdir[], 
		            int numberFITSdisk, char* FITSdir[], 
		            int F_TRUE, int F_FALSE, ObitErr *err) {
  oint lF_TRUE  = F_TRUE;
  oint lF_FALSE = F_FALSE;

  return ObitSystemStartup (pgmName, pgmNumber, AIPSuser,
        	            numberAIPSdisk, AIPSdir,
                	    numberFITSdisk, FITSdir,
                            lF_TRUE, lF_FALSE, err);

} //  end Startup

extern ObitSystem* Shutdown (ObitSystem* in) {
  return ObitSystemShutdown(in);
} // end Shutdown

extern int SystemIsInit (void) {
  gboolean retval;
  olong out;
  retval = ObitSystemIsInit ();
  if (retval) out = 1;
  else out = 0;
  return out;
} // end SystemGetPgmName

extern char* SystemToday (void) {
  return ObitToday ();
} // end SystemToday

extern char* SystemGetPgmName (void) {
  return ObitSystemGetPgmName ();
} // end SystemGetPgmName

extern void SystemSetPgmName (char* pgmName) {
  ObitSystemSetPgmName (pgmName);
} // end SystemSetPgmName

extern int SystemGetPgmNumber (void) {
  return ObitSystemGetPgmNumber ();
} // end SystemGetPgmNumber

extern void SystemSetPgmNumber (int pgmNumber) {
  ObitSystemSetPgmNumber (pgmNumber);
} // end SystemSetPgmNumber

extern int SystemGetAIPSuser (void) {
  return ObitSystemGetAIPSuser ();
} // end SystemGetAIPSuser

extern void SystemSetAIPSuser (int user) {
  ObitSystemSetAIPSuser (user);
} // end SystemGetAIPSuser

extern void MemPrint (void) {
  ObitMemPrint(stdout);
} // end MemPrint

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitSystem *me;
} OSystem;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitSystem *me;
} OSystem;

%addmethods OSystem { 
  OSystem(char *pgmName, int pgmNumber, int AIPSuser,
	     int numberAIPSdisk, char* AIPSdir[], 
	     int numberFITSdisk, char* FITSdir[], 
	     int F_TRUE, int F_FALSE, ObitErr *err) {
     OSystem *out;
     gchar **Adir, **Fdir;

     if (numberAIPSdisk>0) Adir = AIPSdir;
     else Adir = NULL;
     if (numberFITSdisk>0) Fdir = FITSdir;
     else Fdir = NULL;
     out = (OSystem *) malloc(sizeof(OSystem));
     out->me = ObitSystemStartup(pgmName, pgmNumber, AIPSuser,
                                 numberAIPSdisk, Adir,
                                 numberFITSdisk, Fdir,
                                 F_TRUE, F_FALSE, err);
     return out;
   }
  ~OSystem() {
    /* best not self->me = ObitSystemShutdown(self->me);*/
    free(self);
  }
};

/* $Id: OData.inc,v 1.1 2007/08/22 15:18:32 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitData type                          */
/*                                                                    */
/*;  Copyright (C) 2007,2008                                          */
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
#include "ObitData.h"
%}


%inline %{

extern void ODataSetFITS(ObitData *in, int disk, char *file, ObitErr *err) {
  ObitDataSetFITS(in, disk, file, err);
 } // end ODataSetFITS

extern void ODataSetAIPS(ObitData *in, int disk, int cno, int user, 
                        ObitErr *err) {
  ObitDataSetAIPS(in, disk, cno, user, err);
 }

extern ObitData* ODataScratch (ObitData *in, ObitErr *err) {
  return newObitDataScratch (in, err);
} // end  ODataScratch

extern ObitData* ODataZap  (ObitData *in, ObitErr *err) {
  return ObitDataZap (in, err);
} // end DOataZap

extern void ODataRename  (ObitData *in, ObitErr *err) {
  ObitDataRename (in, err);
} // end ODataRename

extern ObitData* ODataCopy  (ObitData *in, ObitData *out, 
			    ObitErr *err) {
  return ObitDataCopy (in, out, err);
} // end  ODataCopy

extern void ODataClone (ObitData *in, ObitData *out, ObitErr *err) {
   return  ObitDataClone (in, out, err);
} // end  ODataClone

// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
// Table version returned as outValue1
extern ObitTable* newODataTable (ObitData *in, int access, 
			      char *tabType, long *outValue1, ObitErr *err) {
  ObitIOAccess laccess;
  olong  loutValue1 = (olong)*outValue1;
  ObitTable *outTable=NULL;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  outTable =  newObitDataTable (in, laccess, tabType, &loutValue1, err);
  *outValue1 = (long)loutValue1;
  return outTable;
} // end  newDataTable

// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
extern ObitHistory* newODataHistory (ObitData *in, int access, 
			            ObitErr *err) {
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  return newObitDataHistory (in, laccess, err);
} // end  newODataTable

extern int ODataZapTable (ObitData *in, char *tabType, long tabVer, 
			 ObitErr *err) {
  ObitIOCode ret;
  ret = ObitDataZapTable (in, tabType, tabVer, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  ODataZapTable

extern int ODataCopyTables (ObitData *in, ObitData *out, char **exclude,
		  	   char **include, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitDataCopyTables  (in, out, exclude, include, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  ODataCopyTables

extern int ODataUpdateTables (ObitData *in, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitDataUpdateTables (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  ODataUpdateTables

// Open and close to fully instantiate
// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
extern int ODataFullInstantiate (ObitData* in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitDataOpen (in, laccess, err);
  ret = ObitDataClose (in, err);
  if ((err->error) || (ret!=OBIT_IO_OK)) return 1;
  else return 0;
} // end ODataFullInstantiate


extern int ODataOpen (ObitData *in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitDataOpen (in, laccess, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end OOpen

// force header update 
extern void ODataDirty (ObitData *in) {
  in->myStatus = OBIT_Modified;
} // end ODirty

extern int ODataClose (ObitData *in, ObitErr *err) {
  ObitIOCode ret;
  ret =  ObitDataClose (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end OClose

extern ObitData* ODataUnref (ObitData* in) {
  if (!ObitDataIsA(in)) return NULL;
  return ObitDataUnref(in);
}

extern ObitData*  ODataRef (ObitData* in) {
  return ObitDataRef(in);
}

extern ObitInfoList* ODataGetList (ObitData* in) {
  return ObitInfoListRef(in->info);
}

extern ObitTableList* ODataGetTableList (ObitData* in) {
  return ObitTableListRef(in->tableList);
}

extern long ODataGetHighVer (ObitData* in, char *tabType) {
  return ObitTableListGetHigh(in->tableList, tabType);
}

extern int ODataisScratch (ObitData* in) {
  return (int)in->isScratch;
}

extern int ODataIsA (ObitData* in) {
  return ObitDataIsA(in);
}

extern char* ODataGetName (ObitData* in) {
  if (ObitDataIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}


%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitData *me;
} OData;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitData *me;
} OData;

%addmethods OData { 
  OData(char *name) {
     OData *out;
     out = (OData *) malloc(sizeof(OData));
     if (strcmp(name, "None")) out->me = newObitData(name);
     else out->me = NULL;
     return out;
   }
  ~OData() { /* Scratch files may be deleted separately*/
   if (self->me->ReferenceCount>0) 
      self->me = ObitDataUnref(self->me);
   free(self);
  }
};

/* $Id: ODisplay.inc,v 1.1 2005/08/03 19:45:35 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitDisplay type                       */
/*                                                                    */
/*;  Copyright (C) 2005                                               */
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
#include "ObitDisplay.h"
#include "ObitImage.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanWindow.h"
%}


%inline %{
extern ObitDisplay* ODisplayCreate(char* name, char* serverURL, ObitErr *err) {
   return  ObitDisplayCreate ((gchar*)name, (gchar*)serverURL, err );
}

extern int ODisplayImage  (ObitDisplay *display, ObitImage *image, ObitErr *err) {
  return (int) ObitDisplayShow (display, (Obit*)image, NULL, 1, err);
} // end ODisplayImage

extern int ODisplayMosaic  (ObitDisplay *display, ObitImageMosaic *mosaic, int field, 
                           ObitErr *err) {
  return (int) ObitDisplayShow (display, (Obit*)mosaic, NULL, (olong)field, err);
} // end ODisplayMosaic

extern int ODisplayImageEdit  (ObitDisplay *display, ObitImage *image, 
                              ObitDConCleanWindow *window, ObitErr *err) {
  return (int) ObitDisplayShow (display, (Obit*)image, window, 1, err);
} // end ODisplayImageEdit

extern int ODisplayMosaicEdit  (ObitDisplay *display, ObitImageMosaic *mosaic, int field, 
                               ObitDConCleanWindow *window, ObitErr *err) {
  return (int) ObitDisplayShow (display, (Obit*)mosaic, window, (olong)field, err);
} // end ODisplayMosaicEdit

extern int ODisplayIsA (ObitDisplay* in) {
  return ObitDisplayIsA(in);
} // end  ODisplayIsA 

ObitDisplay* ODisplayRef (ObitDisplay* in) {
  return ObitDisplayRef (in);
} // end ODisplayRef

ObitDisplay* ODisplayUnref (ObitDisplay* in) {
  if (!ObitDisplayIsA(in)) return NULL;
  return ObitDisplayUnref (in);
} // end ODisplayUnref

extern char* ODisplayGetName (ObitDisplay* in) {
  return in->name;
} // end  ODisplayGetName



%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitDisplay *me;
} ODisplay;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDisplay *me;
} ODisplay;

%addmethods ODisplay { 
  ODisplay(char* name, char* serverURL, ObitErr *err) {
     ODisplay *out;
     char *server;
     out = (ODisplay *) malloc(sizeof(ODisplay));
     if ((serverURL==NULL) || (!strncmp(serverURL, "ObitView", 8))) 
        server = "http://localhost:8765/RPC2";
     else server = serverURL;
     if (strcmp(name, "None")) out->me = ODisplayCreate(name, server, err);
     else out->me = NULL;
     return out;
   }
  ~ODisplay() {
    self->me = ODisplayUnref(self->me);
    free(self);
  }
};

/* $Id: OPlot.inc,v 1.3 2008/02/20 14:58:16 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitPlot type                          */
/*                                                                    */
/*;  Copyright (C) 2006,2008                                          */
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
#include "ObitPlot.h"
%}


%inline %{

/** Public: Create plot. */
extern ObitPlot* OPlotCreate(char* name) {
   return  newObitPlot ((gchar*)name);
}

/** Public: Initialize plot. */
extern void PlotInitPlot (ObitPlot* in, char *output, int color, 
                          int nx, int ny, ObitErr *err) {
  gchar *loutput=NULL;

  if (!strncmp (output, "None", 4)) {
    loutput = NULL;	
  } else {
    loutput = (gchar*)output;	
  }
  ObitPlotInitPlot (in, loutput, (olong)color, (olong)nx, (olong)ny, err);
} // end PlotInitPlot

/** Public: Finalize plot. */
extern void PlotFinishPlot (ObitPlot* in,  ObitErr *err) {
  ObitPlotFinishPlot (in, err);
} // end PlotFinishPlot

/** Public: Copy Plot */
extern ObitPlot* PlotCopy (ObitPlot* in, ObitPlot* out, ObitErr *err){
return  ObitPlotCopy (in, out, err);
} // end ObitPlotCopy

/** Public: Simple X-Y plot. */
extern void PlotXYPlot (ObitPlot* in, int symbol, int n, float *x, float *y, 
               	        ObitErr *err) {
  ObitPlotXYPlot (in, (olong)symbol, (olong)n, (ofloat*)x, (ofloat*)y, err);
} // end PlotXYPlot

/** Public: Simple X-Y over plot. */
extern void PlotXYOver (ObitPlot* in, int symbol, int n, float *x, float *y, 
	                ObitErr *err) {
  ObitPlotXYOver (in, (olong)symbol, (olong)n, (ofloat*)x, (ofloat*)y, err);
} // end PlotXYOver
		     
/** Public: Simple X-Y plot with error bars. */
extern void PlotXYErr (ObitPlot* in, int symbol, int n, float *x, float *y, 
	               float *e, ObitErr *err) {
  ObitPlotXYErr (in, (olong)symbol, (olong)n, (ofloat*)x, (ofloat*)y, (ofloat*)e, err);
} // end PlotXYErr


/** Public: Contour plot of image. */
extern void PlotContour (ObitPlot* in, char *label, ObitImage *image, float lev,
	                 float cntfac, ObitErr *err) {
  ObitPlotContour (in, (gchar*)label, image, (ofloat)lev, (ofloat)cntfac, err);
} // end PlotContour

/** Public: Mark positions on Contour plot of image. */
extern void PlotMarkCross (ObitPlot* in, ObitImage *image, int n,
		 	  double *ra, double *dec, float size, 
		          ObitErr *err) {
  ObitPlotMarkCross (in, image, (olong)n, (odouble*)ra, (odouble*)dec, (ofloat)size, err);
} // end PlotMarkCross


/** Public:  set window and viewport and draw labeled frame */
extern  void PlotSetPlot (ObitPlot* in, float xmin, float xmax, float ymin, float ymax, 
		         int just, int axis, ObitErr *err) {
   ObitPlotSetPlot (in, (ofloat)xmin, (ofloat)xmax, (ofloat)ymin, (ofloat)ymax, 
		    (olong)just, (olong)axis, err);
} // end PlotSetPlot

/** Public: write labels for x-axis, y-axis, and top of plot*/
extern  void PlotLabel (ObitPlot* in, char *xlabel, char *ylabel, char *title,
		        ObitErr *err)  {
  ObitPlotLabel (in, (gchar*)xlabel, (gchar*)ylabel, (gchar*)title,  err);
} // end ObitPlotLabel

/** Public: draw labeled frame around viewport */
extern void PlotDrawAxes (ObitPlot* in,  char *xopt, float xtick, int nxsub, 
		    char *yopt,  float ytick, int nysub, 
		    ObitErr *err) {
   ObitPlotDrawAxes (in, (gchar*)xopt, (ofloat)xtick, (olong)nxsub, 
		   (gchar*)yopt,  (ofloat)ytick, (olong)nysub, 
		    err);
} // end  PlotDrawAxes

/** Public: Scaling for characters */
extern  void PlotSetCharSize (ObitPlot* in, float cscale, ObitErr *err) {
   ObitPlotSetCharSize (in, (ofloat)cscale, err);
} // end PlotSetCharSize

/** Public: Set line width */
extern  void PlotSetLineWidth (ObitPlot* in, int lwidth, ObitErr *err) {
   ObitPlotSetLineWidth (in, (olong)lwidth, err);
} // end PlotSetLineWidth

/** Public: Set foreground color */
extern  void PlotSetColor (ObitPlot* in, int color, ObitErr *err) {
   ObitPlotSetColor (in, (olong)color, err);
} // end PlotSetColor

/** Public: Set/advance subpage */
extern  void PlotSetPage (ObitPlot* in, int sub, ObitErr *err) {
   ObitPlotSetPage (in, (olong)sub, err);
} // end PlotSetPage

/** Public: Write text */
extern void PlotText (ObitPlot* in, float x, float y,
		      float dx, float dy,  float fjust, char *text,
		      ObitErr *err) {
  ObitPlotText (in, (ofloat)x, (ofloat)y, (ofloat)dx, (ofloat)dy, 
		(ofloat)fjust, (gchar*)text,
	         err);
} // end PlotText

/** Public: Write text  relative to port */
extern void PlotRelText (ObitPlot* in, char *side, float disp, 
		      float coord, float fjust, char *text,
		       ObitErr *err) {
  ObitPlotRelText (in, (gchar*)side, (ofloat)disp, 
		(ofloat)coord, (ofloat)fjust, (gchar*)text,
	         err);
} // end PlotRelText

/**  Public: Draw a line..*/
extern void  PlotDrawLine (ObitPlot* in, float x1, float y1, 
                           float x2, float y2, ObitErr *err) {
   ObitPlotDrawLine (in, (ofloat)x1, (ofloat)y1, (ofloat)x2, (ofloat)y2, err);
} // end PlotDrawLine

/**  Public: Draw a curve.*/
extern void  PlotDrawCurve (ObitPlot* in, int n, float *x, float *y, 
                           ObitErr *err) {
   ObitPlotDrawCurve (in, (olong)n, (ofloat*)x, (ofloat*)y, err);
} // end PlotDrawCurve

/**  Public: Draw a Symbol.*/
extern void  PlotDrawSymbol (ObitPlot* in, float x, float y, int symbol, 
                             ObitErr *err) {
   ObitPlotDrawSymbol (in, (ofloat)x, (ofloat)y, (olong)symbol, err);
} // end PlotDrawSymbol

extern int OPlotIsA (ObitPlot* in) {
  return ObitPlotIsA(in);
} // end  OPlotIsA 

extern ObitPlot* OPlotRef (ObitPlot* in) {
  return ObitPlotRef (in);
} // end OPlotRef

extern ObitPlot* OPlotUnref (ObitPlot* in) {
  if (!ObitPlotIsA(in)) return NULL;
  return ObitPlotUnref (in);
} // end OPlotUnref

extern char* OPlotGetName (ObitPlot* in) {
  return in->name;
} // end  OPlotGetName

extern ObitInfoList* PlotGetList (ObitPlot* in) {
  return ObitInfoListRef(in->info);
}



%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitPlot *me;
} OPlot;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitPlot *me;
} OPlot;

%addmethods OPlot { 
  OPlot(char* name) {
     OPlot *out;
     out = (OPlot *) malloc(sizeof(OPlot));
     if (strcmp(name, "None")) out->me = OPlotCreate(name);
     else out->me = NULL;
     return out;
   }
  ~OPlot() {
    self->me = OPlotUnref(self->me);
    free(self);
  }
};

/* $Id: OWindow.inc,v 1.2 2005/12/19 00:12:10 bcotton Exp $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitDConCleanWindow type               */
/*                                                                    */
/*;  Copyright (C) 2005,2008                                          */
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
#include "ObitErr.h"
#include "ObitDConClean.h"
#include "ObitDConCleanWindow.h"
%}

%inline %{
extern ObitDConCleanWindow* OWindowCreate (char *name,
                                           ObitImageMosaic *mosaic, 
                                           ObitErr *err) {
  return ObitDConCleanWindowCreate (name, mosaic, err);
} // end OWindowCreate

extern ObitDConCleanWindow* OWindowCreate1 (char *name,
                                            long naxis[2], 
                                            ObitErr *err) {
  olong lnaxis[2];
  lnaxis[0] = (olong)naxis[0];
  lnaxis[1] = (olong)naxis[1];
  return ObitDConCleanWindowCreate1 (name, lnaxis, err);
} // end OWindowCreate1

extern ObitDConCleanWindow* OWindowCopy (ObitDConCleanWindow* in, 
		              ObitDConCleanWindow* out, ObitErr *err) {
  return ObitDConCleanWindowCopy (in, out, err);
} // end OWindowCopy

// Convert list of lists into windows for field 
// list for each window is ID, type, followed by parameters
// previous entries removed
extern void OWindowSetList(ObitDConCleanWindow* in, PyObject *inList,
                      long field, ObitErr *err) {
  PyObject *list=NULL;
  ObitDConCleanWindowType type;
  olong lfield=field, i, n, iD, id, window[10];
  gchar *routine = "OWindowSetList";

  if  (err->error) return;

  /* Clear previous entries */
  ObitDConCleanWindowDel (in, lfield, -1, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  n = PyList_Size(inList);
  for (iD=0; iD<n; iD++) {
    list = PyList_GetItem(inList, iD);
    type = (ObitDConCleanWindowType)PyInt_AsLong(PyList_GetItem(list, 1));
    /* Ignore iD <0 */
    id = (ObitDConCleanWindowType)PyInt_AsLong(PyList_GetItem(list, 0));
    if (id<0) continue;
    switch (type) {
    case OBIT_DConCleanWindow_rectangle:
      for (i=0; i<4; i++) 
         window[i] = (olong)PyInt_AsLong(PyList_GetItem(list, i+2));
      break;
    case OBIT_DConCleanWindow_round:
      for (i=0; i<3; i++) 
         window[i] = (olong)PyInt_AsLong(PyList_GetItem(list, i+2));
      break;
    default:
      g_error ("Undefined Clean window type");
      return;
    }; /* end switch by window type */

    ObitDConCleanWindowAdd (in, lfield, type, window, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }  /* end loop */
} // end OWindowSetList

// Convert windows for field into a list of lists
// list for each window is ID, type, followed by parameters
extern PyObject* OWindowGetList(ObitDConCleanWindow* in, long field, 
                                ObitErr *err) {
  PyObject *outList=NULL;
  PyObject *list=NULL;
  ObitDConCleanWindowType type;
  olong i, n, iD, *window;

  if  (err->error) return outList;

  outList = PyList_New(in->maxId[field-1]);
  n = 0;
  for (iD=1; iD<=in->maxId[field-1]; iD++) {
    if (ObitDConCleanWindowInfo(in, field, iD, &type, &window, err)) {
      switch (type) {
      case OBIT_DConCleanWindow_rectangle:
        list = PyList_New(6);
        PyList_SetItem(list, 0, PyInt_FromLong((long)iD));
        PyList_SetItem(list, 1, PyInt_FromLong((long)type));
        for (i=0; i<4; i++) 
          PyList_SetItem(list, i+2, PyInt_FromLong((long)window[i]));
        break;
      case OBIT_DConCleanWindow_round:
        list = PyList_New(5);
        PyList_SetItem(list, 0, PyInt_FromLong((long)iD));
        PyList_SetItem(list, 1, PyInt_FromLong((long)type));
        for (i=0; i<3; i++) 
          PyList_SetItem(list, i+2, PyInt_FromLong((long)window[i]));
        break;
      default:
        g_error ("Undefined Clean window type");
        return outList;
      }; /* end switch by window type */
      PyList_SetItem(outList, n++, list);
    }  /* end if iD exists */
  } /* end loop over iD */

  /* if any more list entries fill with -1 */
  for (iD=n; iD<in->maxId[field-1]; iD++) {
    list = PyList_New(2);
    PyList_SetItem(list, 0, PyInt_FromLong((long)-1));
    PyList_SetItem(list, 1, PyInt_FromLong((long)-1));
    PyList_SetItem(outList, iD, list);
  }

  return outList;
} // end OWindowGetList

// return iD of new
extern long OWindowAdd (ObitDConCleanWindow* in, long field, long type, long *window, 
                        ObitErr *err) {
  olong lwindow[4] = {(olong)window[0],(olong)window[1],(olong)window[2],(olong)window[3]}; 
  return (long)ObitDConCleanWindowAdd (in, (olong)field, 
                                       (ObitDConCleanWindowType)type,
                                       lwindow, err);
} // end OWindowAdd

extern void OWindowUpdate (ObitDConCleanWindow* in, long field, long iD, 
                           long type, long *window, ObitErr *err) {
  olong lwindow[4] = {(olong)window[0],(olong)window[1],(olong)window[2],(olong)window[3]}; 
  ObitDConCleanWindowUpdate (in, (olong)field, (olong)iD,
                             (ObitDConCleanWindowType)type,
                             lwindow, err);
} // end OWindowUpdate

extern void OWindowDel (ObitDConCleanWindow* in, long field, long iD, ObitErr *err) {
  ObitDConCleanWindowDel (in, (olong)field, (olong)iD, err);
} // end OWindowDel

extern long OWindowGetMaxID (ObitDConCleanWindow* in, long field) {
  return (long)in->maxId[field-1];
} // end OWindowMaxID

extern ObitDConCleanWindow* OWindowRef (ObitDConCleanWindow* in) {
  return ObitDConCleanWindowRef (in);
} // end OWindowRef

extern ObitDConCleanWindow* OWindowUnref (ObitDConCleanWindow* in) {
  if (!ObitDConCleanWindowIsA(in)) return NULL;
  return ObitDConCleanWindowUnref (in);
} // end OWindowUnref

extern int OWindowIsA (ObitDConCleanWindow* in) {
  return ObitDConCleanWindowIsA(in);
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitDConCleanWindow *me;
} OWindow;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDConCleanWindow *me;
} OWindow;

%addmethods OWindow { 
  OWindow(void) {
     OWindow *out;
     out = (OWindow *) malloc(sizeof(OWindow));
     /* Defer creating Obit Object */
     out->me = NULL;
     return out;
   }
  ~OWindow() {
    self->me = OWindowUnref(self->me);
    free(self);
  }
};

/* $Id: ParserUtil.inc,v 1.1 2007/11/06 18:57:38 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for parameter file parser utilities        */
/*                                                                    */
/*;  Copyright (C) 2007                                               */
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

#include "ObitParser.h"
#include "ObitReturn.h"
%}


%inline %{
/* Parse text file, return 0 = OK, else failed */
extern int Parse (char *infile, ObitInfoList *list, ObitErr *err) {
  ObitIOCode retCode;

  retCode = ObitParserParse ((gchar*)infile, list, err);
  if (retCode==OBIT_IO_OK) return 0;
  else return 1;
} // end Parser

/* Dump InfoList to text file, return 0 = OK, else failed */
extern int Dump (char *outfile, ObitInfoList *list, ObitErr *err) {
  ObitIOCode retCode;

  retCode = ObitReturnDump ((gchar*)outfile, list, err);
  if (retCode==OBIT_IO_OK) return 0;
  else return 1;
} // end Parser

%}


/* $Id: SkyGeom.inc,v 1.3 2007/10/30 13:04:20 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Image utilities                        */
/*                                                                    */
/*;  Copyright (C) 2007                                               */
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
#include "ObitSkyGeom.h"
#include "ObitPBUtil.h"
%}



%inline %{

/** Public: Determine shift between two positions */
void SkyGeomShiftXY (double ra, double dec, float rotate,
			double shiftRA, double shiftDec,
			float* outFlt1, float* outFlt2) {
  ObitSkyGeomShiftXY ((odouble)ra, (odouble)dec, (ofloat)rotate,
		     (odouble)shiftRA, (odouble)shiftDec,
		      outFlt1, outFlt2);
}

/** Public: Determine result of a shift to a position */
void SkyGeomXYShift (double ra, double dec, 
			float xShift, float yShift, float rotate,
			double* outDbl1, double* outDbl2) {
  ObitSkyGeomXYShift ((odouble)ra, (odouble)dec, 
		     (ofloat)xShift, (ofloat)yShift, (ofloat)rotate,
                      outDbl1, outDbl2);
}

/** Public: Returns astronomical coordinates given direction cosines, projection */
void 
SkyGeomNewPos (char *type, double ra0, double dec0, double l, double m, 
		   double* outDbl1, double* outDbl2, int* outInt1) {
  ObitSkyGeomProj Proj = OBIT_SkyGeom_SIN;
  /* Projection type */
  if (!strncmp (type, "-TAN", 4)) Proj = OBIT_SkyGeom_TAN;
  if (!strncmp (type, "-ARC", 4)) Proj = OBIT_SkyGeom_ARC;
  if (!strncmp (type, "-NCP", 4)) Proj = OBIT_SkyGeom_NCP;
  if (!strncmp (type, "-GLS", 4)) Proj = OBIT_SkyGeom_GLS;
  if (!strncmp (type, "-MER", 4)) Proj = OBIT_SkyGeom_MER;
  if (!strncmp (type, "-AIT", 4)) Proj = OBIT_SkyGeom_AIT;
  if (!strncmp (type, "-STG", 4)) Proj = OBIT_SkyGeom_STG;
  ObitSkyGeomNewPos (Proj, (odouble)ra0, (odouble)dec0, (odouble)l, (odouble)m, 
		     outDbl1, outDbl2, outInt1);
}

/** Public: accurate position for pixel coordinates */
int 
SkyGeomWorldPos(float xpix, float ypix, double xref, double yref, 
		    float xrefpix, float yrefpix, float xinc, float yinc, 
		    float rot, char *type, double *outDbl1, double *outDbl2) {
  return ObitSkyGeomWorldPos((ofloat)xpix, (ofloat)ypix, (odouble)xref, (odouble)yref, 
		    (ofloat)xrefpix, (ofloat)yrefpix, (ofloat)xinc, (ofloat)yinc, 
		    (ofloat)rot, (gchar*)type, outDbl1, outDbl2);
}

/** Public: Position for pixel coordinates from IRAF style CD matrix */
int 
SkyGeomCDpos(float xpix, float ypix, double xref, double yref,
		 float xrefpix, float yrefpix, float xinc, float yinc, float rot,
		 float cd1[2], float cd2[2], char *type, double *outDbl1, double *outDbl2) {
  ofloat lcd1[2], lcd2[2];
  lcd1[0] = cd1[0];
  lcd1[1] = cd1[1];
  lcd2[0] = cd2[0];
  lcd2[1] = cd2[1];
  return ObitSkyGeomCDpos((ofloat)xpix, (ofloat)ypix, (odouble)xref, (odouble)yref,
		 (ofloat)xrefpix, (ofloat)yrefpix, (ofloat)xinc, (ofloat)yinc, (ofloat)rot,
		 lcd1, lcd2, (gchar*)type, outDbl1, outDbl2);
}

/** Public: Pixel coordinates for an RA and Dec*/
int 
SkyGeomXYpix(double xpos, double ypos, double xref, double yref, 
		 float xrefpix, float yrefpix, float xinc, float yinc, 
		 float rot, char *type, float *outFlt1, float *outFlt2){
  return ObitSkyGeomXYpix((odouble)xpos, (odouble)ypos, (odouble)xref, (odouble)yref, 
		 (ofloat)xrefpix, (ofloat)yrefpix, (ofloat)xinc, (ofloat)yinc, 
		 (ofloat)rot, (gchar*)type, outFlt1, outFlt2);
}

/** Public:pixel coordinates for an RA and Dec from IRAF  style CD matrix. */
int 
SkyGeomCDpix(double xpos, double ypos, double xref, double yref, 
		 float xrefpix, float yrefpix, float xinc, float yinc, float rot,
		 float icd1[2], float icd2[2], char *type, 
	 	 float *outFlt1, float *outFlt2) {
  ofloat licd1[2], licd2[2];
  licd1[0] = icd1[0];
  licd1[1] = icd1[1];
  licd2[0] = icd2[0];
  licd2[1] = icd2[1];
  return ObitSkyGeomCDpix((odouble)xpos, (odouble)ypos, (odouble)xref, (odouble)yref, 
		 (ofloat)xrefpix, (ofloat)yrefpix, (ofloat)xinc, (ofloat)yinc, (ofloat)rot,
		 licd1, licd2, (gchar*)type, outFlt1, outFlt2);
}

/** Public: Position for pixel coordinates from  offsets from the reference position.*/
int 
SkyGeomWorldPosLM(double dx, double dy, double xref, double yref, 
		      float xinc, float yinc, float rot, char *type, 
		      double *outDbl1, double *outDbl2) {
  return ObitSkyGeomWorldPosLM((odouble)dx, (odouble)dy, (odouble)xref, (odouble)yref, 
		      (ofloat)xinc, (ofloat)yinc, (ofloat)rot, (gchar*)type, 
		      outDbl1, outDbl2);
}

/** Public: Coordinate offsets for an RA and Dec   */
int 
SkyGeomXYPixLM(double xpos, double ypos, double xref, double yref, 
		   float xinc, float yinc, float rot, char *type, 
		   double *outDbl1, double *outDbl2) {
  return ObitSkyGeomXYPixLM((odouble)xpos, (odouble)ypos, (odouble)xref, (odouble)yref, 
		   (ofloat)xinc, (ofloat)yinc, (ofloat)rot, (gchar*)type, 
		    outDbl1, outDbl2);
}

/** Public: Precess B1950 to J2000 coordinates  */
void 
SkyGeomBtoJ (double *outDbl1, double *outDbl2) {
  ObitSkyGeomBtoJ (outDbl1, outDbl2);
}

/** Public: Precess J2000 to B1950 coordinates */
void 
SkyGeomJtoB (double *outDbl1, double *outDbl2) {
  ObitSkyGeomJtoB (outDbl1, outDbl2);
}

/** Public: Convert Equatorial (B1950) to Galactic coordinates  */
void SkyGeomEq2Gal (double *outDbl1, double *outDbl2) {
  ObitSkyGeomEq2Gal (outDbl1, outDbl2);
}

/** Public: Convert Galactic to Equatorial (B1950)  */
void SkyGeomGal2Eq (double *outDbl1, double *outDbl2) {
  ObitSkyGeomGal2Eq (outDbl1, outDbl2);
}

/** Public: Convert Equatorial to Ecliptic coordinates */
void SkyGeomEq2Ec (double *outDbl1, double *outDbl2, float epoch){
  ObitSkyGeomEq2Ec (outDbl1, outDbl2, (ofloat)epoch);
}

/** Public: Convert Ecliptic to Equatorial */
void SkyGeomEc2Eq (double *outDbl1, double *outDbl2, float epoch) {
  ObitSkyGeomEc2Eq (outDbl1, outDbl2, (ofloat)epoch);
}

/** Public: Projection to Zernike plane */
void SkyGeomRADec2Zern (double ra, double dec, float xshift, float yshift, 
			    float* outFlt1, float* outFlt2, int *outInt1) {
  ObitSkyGeomRADec2Zern ((odouble)ra, (odouble)dec, (ofloat)xshift, (ofloat)yshift, 
			    outFlt1, outFlt2, outInt1);
}

/** Public: Compute VLA beam shape from a fitted polynomial */
float PBUtilPoly (double Angle, double Freq, float pbmin) {
  ofloat out;
  out = ObitPBUtilPoly ((odouble)Angle, (odouble)Freq, (ofloat)pbmin);
  return (float)out;
}

/** Public:  Compute Antenna beam assuming uniform illumination of an antenna */
float PBUtilJinc (double Angle, double Freq, float antSize, float pbmin) {
  ofloat out;
  out = ObitPBUtilJinc ((odouble)Angle, (odouble)Freq, (ofloat)antSize, (ofloat)pbmin);
  return (float)out;
}

/** Public:  Calculates the relative gain at a reference frequency */
float PBUtilRelPB (double Angle, int nfreq, double *Freq, float antSize, 
                   float pbmin, double refFreq) {
  ofloat out;
  out = ObitPBUtilRelPB ((odouble)Angle, (olong)nfreq, (odouble*)Freq, 
        (ofloat)antSize, (ofloat)pbmin, (odouble)refFreq);
  return (float)out;
}

%}
/* $Id: SkyModel.inc,v 1.5 2007/11/06 13:00:31 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for SkyModel type                          */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitSkyModel.h"
%}


%inline %{
extern ObitSkyModel* newSkyModel (char* name) {
  return newObitSkyModel (name);
} // end  newSkyModel

extern ObitSkyModel* SkyModelCopy  (ObitSkyModel *in, ObitSkyModel *out, 
				    ObitErr *err) {
  return ObitSkyModelCopy (in, out, err);
} // end  SkyModelCopy

extern ObitSkyModel* SkyModelUnref (ObitSkyModel* in) {
  if (!ObitSkyModelIsA(in)) return NULL;
  return ObitSkyModelUnref(in);
}

extern ObitSkyModel*  SkyModelRef (ObitSkyModel* in) {
  return ObitSkyModelRef(in);
}

extern ObitInfoList* SkyModelGetList (ObitSkyModel* in) {
  return ObitInfoListRef(in->info);
}

extern ObitImageMosaic* SkyModelGetImageMosaic (ObitSkyModel* in) {
  return ObitImageMosaicRef(in->mosaic);
}

extern void SkyModelSetImageMosaic (ObitSkyModel* in, ObitImageMosaic *mosaic, 
                              ObitErr *err) {
  in->mosaic = ObitImageMosaicUnref(in->mosaic);  /* Out with the old */
  in->mosaic = ObitImageMosaicRef(mosaic);        /* In with the new */
}

extern ObitSkyModel* SkyModelCreate (char *name, ObitImageMosaic* mosaic) {
 return ObitSkyModelCreate(name, mosaic);
}

extern void SkyModelSubUV (ObitSkyModel* in, ObitUV *inData, ObitUV *outData, ObitErr *err) {
 ObitSkyModelSubUV(in, inData, outData, err);
}

extern void SkyModelDivUV (ObitSkyModel* in, ObitUV *inData, ObitUV *outData, ObitErr *err) {
 ObitSkyModelDivUV(in, inData, outData, err);
}

extern void SkyModelCompressCC (ObitSkyModel* in, ObitErr *err) {
   ObitSkyModelCompressCC (in, err);
}

extern char* SkyModelGetName (ObitSkyModel* in) {
  if (ObitSkyModelIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int SkyModelIsA (ObitSkyModel* in) {
  return ObitSkyModelIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitSkyModel *me;
} SkyModel;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitSkyModel *me;
} SkyModel;

%addmethods SkyModel { 
  SkyModel(char* name) {
     SkyModel *out;
     out = (SkyModel *) malloc(sizeof(SkyModel));
     if (strcmp(name, "None")) out->me = newSkyModel(name);
     else out->me = NULL;
     return out;
   }
  ~SkyModel() {
   if (self->me-> ReferenceCount>0) 
      self->me = SkyModelUnref(self->me);
   free(self);
  }
};

/* $Id:   */  
/*--------------------------------------------------------------------*/
/* Swig module description for SkyModelVMIon type                     */
/*                                                                    */
/*;  Copyright (C) 2007                                               */
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
#include "ObitSkyModelVMIon.h"
%}


%inline %{
extern ObitSkyModelVMIon* newSkyModelVMIon (char* name) {
  return newObitSkyModelVMIon (name);
} // end  newSkyModelVMIon

extern ObitSkyModelVMIon* SkyModelVMIonCopy  (ObitSkyModelVMIon *in, ObitSkyModelVMIon *out, 
				    ObitErr *err) {
  return ObitSkyModelVMIonCopy (in, out, err);
} // end  SkyModelVMIonCopy

extern ObitSkyModelVMIon* SkyModelVMIonUnref (ObitSkyModelVMIon* in) {
  if (!ObitSkyModelVMIonIsA(in)) return NULL;
  return ObitSkyModelVMIonUnref(in);
}

extern ObitSkyModelVMIon*  SkyModelVMIonRef (ObitSkyModelVMIon* in) {
  return ObitSkyModelVMIonRef(in);
}

extern ObitInfoList* SkyModelVMIonGetList (ObitSkyModelVMIon* in) {
  return ObitInfoListRef(in->info);
}

extern ObitSkyModelVMIon* SkyModelVMIonCreate (char *name, ObitImageMosaic* mosaic) {
 return ObitSkyModelVMIonCreate(name, mosaic);
}

extern char* SkyModelVMIonGetName (ObitSkyModelVMIon* in) {
  if (ObitSkyModelVMIonIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int SkyModelVMIonIsA (ObitSkyModelVMIon* in) {
  return ObitSkyModelVMIonIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitSkyModelVMIon *me;
} SkyModelVMIon;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitSkyModelVMIon *me;
} SkyModelVMIon;

%addmethods SkyModelVMIon { 
  SkyModelVMIon(char* name) {
     SkyModelVMIon *out;
     out = (SkyModelVMIon *) malloc(sizeof(SkyModelVMIon));
     if (strcmp(name, "None")) out->me = newSkyModelVMIon(name);
     else out->me = NULL;
     return out;
   }
  ~SkyModelVMIon() {
   if (self->me-> ReferenceCount>0) 
      self->me = SkyModelVMIonUnref(self->me);
   free(self);
  }
};

/* $Id: SpectrumFit.inc,v 1.4 2008/05/16 23:29:02 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for SpectrumFit type                       */
/*                                                                    */
/*;  Copyright (C) 2008                                               */
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
#include "ObitSpectrumFit.h"
#include "ObitImage.h"
%}


%inline %{
extern ObitSpectrumFit* newSpectrumFit (char* name) {
  return newObitSpectrumFit (name);
} // end  newSpectrumFit

extern ObitSpectrumFit* SpectrumFitCopy  (ObitSpectrumFit *in, ObitSpectrumFit *out, 
				    ObitErr *err) {
  return ObitSpectrumFitCopy (in, out, err);
} // end  SpectrumFitCopy

extern ObitSpectrumFit* SpectrumFitUnref (ObitSpectrumFit* in) {
  if (!ObitSpectrumFitIsA(in)) return NULL;
  return ObitSpectrumFitUnref(in);
}

extern ObitSpectrumFit*  SpectrumFitRef (ObitSpectrumFit* in) {
  return ObitSpectrumFitRef(in);
}

extern ObitSpectrumFit* SpectrumFitCreate (char *name, int nterm) {
  return ObitSpectrumFitCreate((gchar*)name, (olong)nterm);
}

extern void SpectrumFitCube (ObitSpectrumFit* in, ObitImage *inImage, 
			     ObitImage *outImage, ObitErr *err) {
  ObitSpectrumFitCube(in, inImage, outImage, err);
}

extern void SpectrumFitImArr (ObitSpectrumFit* in, int nimage, ObitImage **imArr, 
			      ObitImage *outImage, ObitErr *err) {
  ObitSpectrumFitImArr(in, (olong)nimage, imArr, outImage, err);
} // end SpectrumFitImArr 

extern void SpectrumFitEval (ObitSpectrumFit* in, ObitImage *inImage,  
                             double outFreq, ObitImage *outImage, 
                             ObitErr *err) {
  ObitSpectrumFitEval(in, inImage, (odouble)outFreq, outImage, err);
}

extern PyObject* SpectrumFitSingle (int nfreq, int nterm, double *freq, float *flux, float *sigma,
                                 ObitErr *err) {
  ofloat *out=NULL;
  olong i, n;
  PyObject *outList=NULL, *o=NULL;

  out = ObitSpectrumFitSingle((olong)nfreq, (olong)nterm, (odouble*)freq, (ofloat*)flux, (ofloat*)sigma, err);
  if (err->error) {
        ObitErrLog(err);
        PyErr_SetString(PyExc_TypeError,"Spectral Fit failed");
	o = PyString_FromString("FAILED");
        return o;
  }
  n = 1 + nterm*2;
  outList = PyList_New(n); 
  for (i=0; i<n; i++) {
    o = PyFloat_FromDouble((double)out[i]);
    PyList_SetItem(outList, i, o);
  }
  if (out) g_free(out);
  return outList;
}  // end SpectrumFitPixel

extern ObitInfoList* SpectrumFitGetList (ObitSpectrumFit* in) {
  return ObitInfoListRef(in->info);
}

extern char* SpectrumFitGetName (ObitSpectrumFit* in) {
  if (ObitSpectrumFitIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int SpectrumFitIsA (ObitSpectrumFit* in) {
  return ObitSpectrumFitIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitSpectrumFit *me;
} SpectrumFit;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitSpectrumFit *me;
} SpectrumFit;

%addmethods SpectrumFit { 
  SpectrumFit(char* name, int nterm) {
     SpectrumFit *out;
     out = (SpectrumFit *) malloc(sizeof(SpectrumFit));
     if (strcmp(name, "None")) out->me = SpectrumFitCreate(name, nterm);
     else out->me = NULL;
     return out;
   }
  ~SpectrumFit() {
   if (!self) return;  // Not defined
   if (self && self->me && self->me->ReferenceCount>0) {
      self->me = SpectrumFitUnref(self->me);
      free(self);
   }
  }
};

/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableAN.h"
%}
 
%inline %{
 
extern ObitTable* TableAN (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numOrb, int numPCal,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumOrb = (oint)numOrb;
   oint lnumPCal = (oint)numPCal;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableANValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumOrb, lnumPCal,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableANGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableAN *lTab = (ObitTableAN*)inTab;
  PyDict_SetItemString(outDict, "numOrb",  PyInt_FromLong((long)lTab->numOrb));
  PyDict_SetItemString(outDict, "numPCal",  PyInt_FromLong((long)lTab->numPCal));
  PyDict_SetItemString(outDict, "ArrayX",  PyFloat_FromDouble((double)lTab->ArrayX));
  PyDict_SetItemString(outDict, "ArrayY",  PyFloat_FromDouble((double)lTab->ArrayY));
  PyDict_SetItemString(outDict, "ArrayZ",  PyFloat_FromDouble((double)lTab->ArrayZ));
  PyDict_SetItemString(outDict, "GSTiat0",  PyFloat_FromDouble((double)lTab->GSTiat0));
  PyDict_SetItemString(outDict, "DegDay",  PyFloat_FromDouble((double)lTab->DegDay));
  PyDict_SetItemString(outDict, "Freq",  PyFloat_FromDouble((double)lTab->Freq));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));
  PyDict_SetItemString(outDict, "PolarX",  PyFloat_FromDouble((double)lTab->PolarX));
  PyDict_SetItemString(outDict, "PolarY",  PyFloat_FromDouble((double)lTab->PolarY));
  PyDict_SetItemString(outDict, "dataUtc",  PyFloat_FromDouble((double)lTab->dataUtc));
  PyDict_SetItemString(outDict, "TimeSys", PyString_InternFromString(lTab->TimeSys));
  PyDict_SetItemString(outDict, "FreqID",  PyInt_FromLong((long)lTab->FreqID));
  PyDict_SetItemString(outDict, "iatUtc",  PyFloat_FromDouble((double)lTab->iatUtc));
  PyDict_SetItemString(outDict, "polType", PyString_InternFromString(lTab->polType));
  PyDict_SetItemString(outDict, "P_Refant",  PyInt_FromLong((long)lTab->P_Refant));
  PyDict_SetItemString(outDict, "P_Diff01",  PyFloat_FromDouble((double)lTab->P_Diff01));
  PyDict_SetItemString(outDict, "P_Diff02",  PyFloat_FromDouble((double)lTab->P_Diff02));
  PyDict_SetItemString(outDict, "P_Diff03",  PyFloat_FromDouble((double)lTab->P_Diff03));
  PyDict_SetItemString(outDict, "P_Diff04",  PyFloat_FromDouble((double)lTab->P_Diff04));
  PyDict_SetItemString(outDict, "P_Diff05",  PyFloat_FromDouble((double)lTab->P_Diff05));
  PyDict_SetItemString(outDict, "P_Diff06",  PyFloat_FromDouble((double)lTab->P_Diff06));
  PyDict_SetItemString(outDict, "P_Diff07",  PyFloat_FromDouble((double)lTab->P_Diff07));
  PyDict_SetItemString(outDict, "P_Diff08",  PyFloat_FromDouble((double)lTab->P_Diff08));

  return outDict;
} 

extern void TableANSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableAN *lTab = (ObitTableAN*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEAN;

  lTab->ArrayX = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ArrayX"));
  lTab->ArrayY = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ArrayY"));
  lTab->ArrayZ = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ArrayZ"));
  lTab->GSTiat0 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "GSTiat0"));
  lTab->DegDay = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "DegDay"));
  lTab->Freq = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "Freq"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;
  lTab->PolarX = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "PolarX"));
  lTab->PolarY = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "PolarY"));
  lTab->dataUtc = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "dataUtc"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "TimeSys"));
  strncpy (lTab->TimeSys, tstr, lstr); lTab->TimeSys[lstr-1]=0;
  lTab->FreqID = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "FreqID"));
  lTab->iatUtc = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "iatUtc"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "polType"));
  strncpy (lTab->polType, tstr, lstr); lTab->polType[lstr-1]=0;
  lTab->P_Refant = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "P_Refant"));
  lTab->P_Diff01 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "P_Diff01"));
  lTab->P_Diff02 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "P_Diff02"));
  lTab->P_Diff03 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "P_Diff03"));
  lTab->P_Diff04 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "P_Diff04"));
  lTab->P_Diff05 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "P_Diff05"));
  lTab->P_Diff06 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "P_Diff06"));
  lTab->P_Diff07 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "P_Diff07"));
  lTab->P_Diff08 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "P_Diff08"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableAT.h"
%}
 
%inline %{
 
extern ObitTable* TableAT (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numBand,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumBand = (oint)numBand;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableATValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumBand,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableATGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableAT *lTab = (ObitTableAT*)inTab;
  PyDict_SetItemString(outDict, "numBand",  PyInt_FromLong((long)lTab->numBand));
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));
  PyDict_SetItemString(outDict, "obscode", PyString_InternFromString(lTab->obscode));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));
  PyDict_SetItemString(outDict, "polType", PyString_InternFromString(lTab->polType));
  PyDict_SetItemString(outDict, "FreqID",  PyInt_FromLong((long)lTab->FreqID));

  return outDict;
} 

extern void TableATSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableAT *lTab = (ObitTableAT*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEAT;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obscode"));
  strncpy (lTab->obscode, tstr, lstr); lTab->obscode[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "polType"));
  strncpy (lTab->polType, tstr, lstr); lTab->polType[lstr-1]=0;
  lTab->FreqID = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "FreqID"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableBL.h"
%}
 
%inline %{
 
extern ObitTable* TableBL (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numPol, int numIF,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumPol = (oint)numPol;
   oint lnumIF = (oint)numIF;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableBLValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumPol, lnumIF,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableBLGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableBL *lTab = (ObitTableBL*)inTab;
  PyDict_SetItemString(outDict, "numPol",  PyInt_FromLong((long)lTab->numPol));
  PyDict_SetItemString(outDict, "numIF",  PyInt_FromLong((long)lTab->numIF));

  return outDict;
} 

extern void TableBLSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableBL *lTab = (ObitTableBL*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEBL;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableBP.h"
%}
 
%inline %{
 
extern ObitTable* TableBP (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numPol, int numIF, int numChan,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumPol = (oint)numPol;
   oint lnumIF = (oint)numIF;
   oint lnumChan = (oint)numChan;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableBPValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumPol, lnumIF, lnumChan,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableBPGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableBP *lTab = (ObitTableBP*)inTab;
  PyDict_SetItemString(outDict, "numPol",  PyInt_FromLong((long)lTab->numPol));
  PyDict_SetItemString(outDict, "numIF",  PyInt_FromLong((long)lTab->numIF));
  PyDict_SetItemString(outDict, "numChan",  PyInt_FromLong((long)lTab->numChan));
  PyDict_SetItemString(outDict, "startChan",  PyInt_FromLong((long)lTab->startChan));
  PyDict_SetItemString(outDict, "numShifts",  PyInt_FromLong((long)lTab->numShifts));
  PyDict_SetItemString(outDict, "lowShift",  PyInt_FromLong((long)lTab->lowShift));
  PyDict_SetItemString(outDict, "shiftInc",  PyInt_FromLong((long)lTab->shiftInc));

  return outDict;
} 

extern void TableBPSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableBP *lTab = (ObitTableBP*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEBP;

  lTab->startChan = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "startChan"));
  lTab->numShifts = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "numShifts"));
  lTab->lowShift = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "lowShift"));
  lTab->shiftInc = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "shiftInc"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableCC.h"
%}
 
%inline %{
 
extern ObitTable* TableCC (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int noParms,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnoParms = (oint)noParms;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableCCValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnoParms,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableCCGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableCC *lTab = (ObitTableCC*)inTab;
  PyDict_SetItemString(outDict, "noParms",  PyInt_FromLong((long)lTab->noParms));

  return outDict;
} 

extern void TableCCSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableCC *lTab = (ObitTableCC*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLECC;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableCL.h"
%}
 
%inline %{
 
extern ObitTable* TableCL (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numPol, int numIF, int numTerm,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumPol = (oint)numPol;
   oint lnumIF = (oint)numIF;
   oint lnumTerm = (oint)numTerm;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableCLValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumPol, lnumIF, lnumTerm,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableCLGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableCL *lTab = (ObitTableCL*)inTab;
  PyDict_SetItemString(outDict, "numPol",  PyInt_FromLong((long)lTab->numPol));
  PyDict_SetItemString(outDict, "numIF",  PyInt_FromLong((long)lTab->numIF));
  PyDict_SetItemString(outDict, "numTerm",  PyInt_FromLong((long)lTab->numTerm));
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));
  PyDict_SetItemString(outDict, "mGMod",  PyFloat_FromDouble((double)lTab->mGMod));

  return outDict;
} 

extern void TableCLSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableCL *lTab = (ObitTableCL*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLECL;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));
  lTab->mGMod = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "mGMod"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableCQ.h"
%}
 
%inline %{
 
extern ObitTable* TableCQ (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numIF,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumIF = (oint)numIF;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableCQValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumIF,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableCQGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableCQ *lTab = (ObitTableCQ*)inTab;
  PyDict_SetItemString(outDict, "numIF",  PyInt_FromLong((long)lTab->numIF));
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));

  return outDict;
} 

extern void TableCQSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableCQ *lTab = (ObitTableCQ*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLECQ;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableCT.h"
%}
 
%inline %{
 
extern ObitTable* TableCT (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numBand,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumBand = (oint)numBand;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableCTValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumBand,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableCTGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableCT *lTab = (ObitTableCT*)inTab;
  PyDict_SetItemString(outDict, "numBand",  PyInt_FromLong((long)lTab->numBand));
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));
  PyDict_SetItemString(outDict, "obscode", PyString_InternFromString(lTab->obscode));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));

  return outDict;
} 

extern void TableCTSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableCT *lTab = (ObitTableCT*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLECT;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obscode"));
  strncpy (lTab->obscode, tstr, lstr); lTab->obscode[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageDesc type                         */
/*                                                                    */
/*;  Copyright (C) 2005,2008                                          */
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
#include "ObitTableDesc.h"
%}

%inline %{
// Routine to remove trailing blanks from a string 
static void TableDescDeeBlank (gchar *in) 
{
  olong i;

  for (i=strlen(in)-1; i>=0; i--) {
     if (in[i]==' ') in[i] = 0;
     else if (in[i]!=' ') break;
  }

} // end TableDescDeeBlank

extern ObitTableDesc* TableDescCreate (char *name) {
  return newObitTableDesc (name);
} // end TableDescCreate

extern ObitTableDesc* TableDescCopy (ObitTableDesc* in, 
		              ObitTableDesc* out, ObitErr *err) {
  return ObitTableDescCopy (in, out, err);
} // end TableDescCopy

extern void TableDescCopyDesc (ObitTableDesc* in, ObitTableDesc* out,
			ObitErr *err) {
  ObitTableDescCopyDesc  (in, out, err);
} // end TableDescCopyDesc

extern void TableDescIndex (ObitTableDesc* in) {
  ObitTableDescIndex (in);
} // end TableDescIndex

extern ObitInfoList* TableDescGetList (ObitTableDesc* in) {
  return ObitInfoListRef(in->info);
}
 
extern PyObject *TableDescGetDict(ObitTableDesc* in) {
  PyObject *outDict = PyDict_New();
  PyObject *list, *value;
  gchar *ctemp;
  int i, pos = 0;

  /* test validity */
  if (in->nfield<=0) {
    PyErr_SetString(PyExc_TypeError,"Input not fully defined");
    return outDict;
  }

  PyDict_SetItemString(outDict, "Table name", PyString_InternFromString(in->TableName));
  PyDict_SetItemString(outDict, "version",PyInt_FromLong((long)in->version));
  PyDict_SetItemString(outDict, "nrow",   PyInt_FromLong((long)in->nrow));
  PyDict_SetItemString(outDict, "lrow",   PyInt_FromLong((long)in->lrow));
  PyDict_SetItemString(outDict, "sortOrder1", PyInt_FromLong((long)in->sort[0]));
  PyDict_SetItemString(outDict, "sortOrder2", PyInt_FromLong((long)in->sort[1]));

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) PyList_SetItem(list, i, PyInt_FromLong((long)in->repeat[i]));
  PyDict_SetItemString(outDict, "repeat", list);

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) PyList_SetItem(list, i, PyInt_FromLong((long)in->dim[i][0]));
  PyDict_SetItemString(outDict, "dim0", list);

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) PyList_SetItem(list, i, PyInt_FromLong((long)in->dim[i][1]));
  PyDict_SetItemString(outDict, "dim1", list);

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) PyList_SetItem(list, i, PyInt_FromLong((long)in->dim[i][2]));
  PyDict_SetItemString(outDict, "dim2", list);

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) PyList_SetItem(list, i, PyInt_FromLong((long)in->type[i]));
  PyDict_SetItemString(outDict, "type", list);

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) {
    ctemp =  g_strdup(in->FieldName[i]);
    TableDescDeeBlank(ctemp);
    PyList_SetItem(list, i, PyString_InternFromString(ctemp));
    g_free (ctemp);
  }
  PyDict_SetItemString(outDict, "FieldName", list);

  list = PyList_New(in->nfield);
  for (i=0; i<in->nfield; i++) PyList_SetItem(list, i, PyString_InternFromString(in->FieldUnit[i]));
  for (i=0; i<in->nfield; i++) {
    ctemp =  g_strdup(in->FieldUnit[i]);
    TableDescDeeBlank(ctemp);
    PyList_SetItem(list, i, PyString_InternFromString(ctemp));
    g_free (ctemp);
  }
  PyDict_SetItemString(outDict, "FieldUnit", list);

  /* Discard references to newly created objects. */
  while (PyDict_Next(outDict, &pos, NULL, &value))
    Py_DECREF(value);

  return outDict;
} // end TableDescGetDict

extern void TableDescSetDict(ObitTableDesc* in, PyObject *inDict) {
  PyObject *list;
  char *tstr;
  int i, number;

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return;
  }

  in->version  = PyInt_AsLong(PyDict_GetItemString(inDict, "version"));
  in->sort[0]  = PyInt_AsLong(PyDict_GetItemString(inDict, "sortOrder1"));
  in->sort[1]  = PyInt_AsLong(PyDict_GetItemString(inDict, "sortOrder2"));

  list = PyDict_GetItemString(inDict, "FieldName");
  number = MIN (in->nfield, PyList_Size(list));
  for (i=0; i<number; i++) {
    tstr = PyString_AsString(PyList_GetItem(list, i));
    if (in->FieldName[i]) g_free(in->FieldName[i]);
    in->FieldName[i] = g_strdup(tstr);
  }

  list = PyDict_GetItemString(inDict, "FieldUnit");
  number = MIN (in->nfield, PyList_Size(list));
  for (i=0; i<number; i++) {
    tstr = PyString_AsString(PyList_GetItem(list, i));
    if (in->FieldUnit[i]) g_free(in->FieldUnit[i]);
    in->FieldUnit[i] = g_strdup(tstr);
  }

} // end TableDescSetDict

//  Define a descriptor from a python dict
extern ObitTableDesc* TableDescDef(PyObject *inDict) {
  ObitTableDesc *out=NULL;
  PyObject *fieldname, *fieldunit, *dim0, *dim1, *dim2, *repeat, *type;
  PyObject *thing;
  gchar *tstr;
  olong i, nfield;

 
  // Define output
  out = newObitTableDesc(NULL);

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return out;
  }

  // Check input
  fieldname = PyDict_GetItemString(inDict, "FieldName");
  if (!fieldname) {
    PyErr_SetString(PyExc_TypeError,"FieldName Array not found");
    return out;
  }

  nfield = PyList_Size(fieldname);

  fieldunit = PyDict_GetItemString(inDict, "FieldUnit");
  if (!fieldunit) {
    PyErr_SetString(PyExc_TypeError,"FieldUnit Array not found");
    return out;
  }
  if (PyList_Size(fieldunit)!=nfield) {
    PyErr_SetString(PyExc_TypeError,"FieldUnit Array wrong dimension");
    return out;
  }

  repeat = PyDict_GetItemString(inDict, "repeat");
  if (!repeat) {
    PyErr_SetString(PyExc_TypeError,"repeat Array not found");
    return;
  }
  if (PyList_Size(repeat)!=nfield) {
    PyErr_SetString(PyExc_TypeError,"repeat Array wrong dimension");
    return out;
  }

  dim0   = PyDict_GetItemString(inDict, "dim0");
  if (!dim0) {
    PyErr_SetString(PyExc_TypeError,"dim0 Array not found");
    return out;
  }
  if (PyList_Size(dim0)!=nfield) {
    PyErr_SetString(PyExc_TypeError,"dim0 Array wrong dimension");
    return out;
  }

  dim1   = PyDict_GetItemString(inDict, "dim1");
  if (!dim1) {
    PyErr_SetString(PyExc_TypeError,"dim1 Array not found");
    return out;
  }
  if (PyList_Size(dim1)!=nfield) {
    PyErr_SetString(PyExc_TypeError,"dim1 Array wrong dimension");
    return out;
  }

  dim2   = PyDict_GetItemString(inDict, "dim2");
  if (!dim2) {
    PyErr_SetString(PyExc_TypeError,"dim2 Array not found");
    return out;
  }
  if (PyList_Size(dim2)!=nfield) {
    PyErr_SetString(PyExc_TypeError,"dim2 Array wrong dimension");
    return out;
  }

  type   = PyDict_GetItemString(inDict, "type");
  if (!type) {
    PyErr_SetString(PyExc_TypeError,"type Array not found");
    return out;
  }
  if (PyList_Size(type)!=nfield) {
    PyErr_SetString(PyExc_TypeError,"type Array wrong dimension");
    return out;
  }

  // Resize output
  ObitTableDescRealloc (out, nfield);

  thing = PyDict_GetItemString(inDict, "version");
  if (!thing) {
    PyErr_SetString(PyExc_TypeError,"version not found");
    return out;
  }
  out->version  = PyInt_AsLong(thing);
  thing = PyDict_GetItemString(inDict, "sortOrder1");
  if (!thing) {
    PyErr_SetString(PyExc_TypeError,"sortOrder1 not found");
    return out;
  }
  out->sort[0]  = PyInt_AsLong(thing);
  thing = PyDict_GetItemString(inDict, "sortOrder2");
  if (!thing) {
    PyErr_SetString(PyExc_TypeError,"sortOrder2 not found");
    return out;
  }
  out->sort[1]  = PyInt_AsLong(thing);
  thing = PyDict_GetItemString(inDict, "Table name");
  if (!thing) {
    PyErr_SetString(PyExc_TypeError,"Table name not found");
    return out;
  }
  tstr = PyString_AsString(thing);
  if (out->TableName) g_free(out->TableName);
  out->TableName = g_strdup(tstr);

  // field names
  for (i=0; i<nfield; i++) {
    tstr = PyString_AsString(PyList_GetItem(fieldname, i));
    if (out->FieldName[i]) g_free(out->FieldName[i]);
    out->FieldName[i] = g_strdup(tstr);
  }

  // field units
  for (i=0; i<nfield; i++) {
    tstr = PyString_AsString(PyList_GetItem(fieldunit, i));
    if (out->FieldUnit[i]) g_free(out->FieldUnit[i]);
    out->FieldUnit[i] = g_strdup(tstr);
  }

  // other field stuff
  for (i=0; i<nfield; i++) {
    out->repeat[i] = (olong)PyInt_AsLong(PyList_GetItem(repeat, i));
    out->type[i]   = (olong)PyInt_AsLong(PyList_GetItem(type, i));
    out->dim[i][0] = (olong)PyInt_AsLong(PyList_GetItem(dim0, i));
    out->dim[i][1] = (olong)PyInt_AsLong(PyList_GetItem(dim1, i));
    out->dim[i][2] = (olong)PyInt_AsLong(PyList_GetItem(dim2, i));
  }

  ObitTableDescIndex (out);

  return out;

} // end TableDescDef

ObitTableDesc* TableDescRef (ObitTableDesc* in) {
  return ObitTableDescRef (in);
} // end TableDescRef

ObitTableDesc* TableDescUnref (ObitTableDesc* in) {
  if (!ObitTableDescIsA(in)) return NULL;
  return ObitTableDescUnref (in);
} // end TableDescUnref

extern int TableDescIsA (ObitTableDesc* in) {
  return ObitTableDescIsA(in);
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitTableDesc *me;
} TableDesc;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitTableDesc *me;
} TableDesc;

%addmethods TableDesc { 
  TableDesc(char *name) {
     TableDesc *out;
     out = (TableDesc *) malloc(sizeof(TableDesc));
     if (strcmp(name, "None")) out->me = TableDescCreate (name);
     else out->me = NULL;
     return out;
   }
  ~TableDesc() {
    self->me = TableDescUnref(self->me);
    free(self);
  }
};

/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableFG.h"
%}
 
%inline %{
 
extern ObitTable* TableFG (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableFGValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableFGGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableFG *lTab = (ObitTableFG*)inTab;

  return outDict;
} 

extern void TableFGSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableFG *lTab = (ObitTableFG*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEFG;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableFQ.h"
%}
 
%inline %{
 
extern ObitTable* TableFQ (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numIF,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumIF = (oint)numIF;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableFQValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumIF,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableFQGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableFQ *lTab = (ObitTableFQ*)inTab;
  PyDict_SetItemString(outDict, "numIF",  PyInt_FromLong((long)lTab->numIF));

  return outDict;
} 

extern void TableFQSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableFQ *lTab = (ObitTableFQ*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEFQ;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableGC.h"
%}
 
%inline %{
 
extern ObitTable* TableGC (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numBand, int numPol, int numTabs,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumBand = (oint)numBand;
   oint lnumPol = (oint)numPol;
   oint lnumTabs = (oint)numTabs;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableGCValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumBand, lnumPol, lnumTabs,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableGCGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableGC *lTab = (ObitTableGC*)inTab;
  PyDict_SetItemString(outDict, "numBand",  PyInt_FromLong((long)lTab->numBand));
  PyDict_SetItemString(outDict, "numPol",  PyInt_FromLong((long)lTab->numPol));
  PyDict_SetItemString(outDict, "numTabs",  PyInt_FromLong((long)lTab->numTabs));
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));

  return outDict;
} 

extern void TableGCSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableGC *lTab = (ObitTableGC*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEGC;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableHistory.h"
%}
 
%inline %{
 
extern ObitTable* TableHistory (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableHistoryValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableHistoryGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableHistory *lTab = (ObitTableHistory*)inTab;

  return outDict;
} 

extern void TableHistorySetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableHistory *lTab = (ObitTableHistory*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEHistory;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableIDI_ANTENNA.h"
%}
 
%inline %{
 
extern ObitTable* TableIDI_ANTENNA (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int no_band, int numPCal,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lno_band = (oint)no_band;
   oint lnumPCal = (oint)numPCal;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableIDI_ANTENNAValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lno_band, lnumPCal,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableIDI_ANTENNAGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableIDI_ANTENNA *lTab = (ObitTableIDI_ANTENNA*)inTab;
  PyDict_SetItemString(outDict, "no_band",  PyInt_FromLong((long)lTab->no_band));
  PyDict_SetItemString(outDict, "numPCal",  PyInt_FromLong((long)lTab->numPCal));
  PyDict_SetItemString(outDict, "tabrev",  PyInt_FromLong((long)lTab->tabrev));
  PyDict_SetItemString(outDict, "no_stkd",  PyInt_FromLong((long)lTab->no_stkd));
  PyDict_SetItemString(outDict, "stk_1",  PyInt_FromLong((long)lTab->stk_1));
  PyDict_SetItemString(outDict, "no_chan",  PyInt_FromLong((long)lTab->no_chan));
  PyDict_SetItemString(outDict, "ref_freq",  PyFloat_FromDouble((double)lTab->ref_freq));
  PyDict_SetItemString(outDict, "chan_bw",  PyFloat_FromDouble((double)lTab->chan_bw));
  PyDict_SetItemString(outDict, "ref_pixl",  PyInt_FromLong((long)lTab->ref_pixl));
  PyDict_SetItemString(outDict, "obscode", PyString_InternFromString(lTab->obscode));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));
  PyDict_SetItemString(outDict, "polType", PyString_InternFromString(lTab->polType));

  return outDict;
} 

extern void TableIDI_ANTENNASetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableIDI_ANTENNA *lTab = (ObitTableIDI_ANTENNA*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEIDI_ANTENNA;

  lTab->tabrev = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tabrev"));
  lTab->no_stkd = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "no_stkd"));
  lTab->stk_1 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "stk_1"));
  lTab->no_chan = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "no_chan"));
  lTab->ref_freq = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ref_freq"));
  lTab->chan_bw = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "chan_bw"));
  lTab->ref_pixl = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "ref_pixl"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obscode"));
  strncpy (lTab->obscode, tstr, lstr); lTab->obscode[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "polType"));
  strncpy (lTab->polType, tstr, lstr); lTab->polType[lstr-1]=0;

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableIDI_ARRAY_GEOMETRY.h"
%}
 
%inline %{
 
extern ObitTable* TableIDI_ARRAY_GEOMETRY (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int no_band, int numOrb,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lno_band = (oint)no_band;
   oint lnumOrb = (oint)numOrb;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableIDI_ARRAY_GEOMETRYValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lno_band, lnumOrb,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableIDI_ARRAY_GEOMETRYGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableIDI_ARRAY_GEOMETRY *lTab = (ObitTableIDI_ARRAY_GEOMETRY*)inTab;
  PyDict_SetItemString(outDict, "no_band",  PyInt_FromLong((long)lTab->no_band));
  PyDict_SetItemString(outDict, "numOrb",  PyInt_FromLong((long)lTab->numOrb));
  PyDict_SetItemString(outDict, "tabrev",  PyInt_FromLong((long)lTab->tabrev));
  PyDict_SetItemString(outDict, "no_stkd",  PyInt_FromLong((long)lTab->no_stkd));
  PyDict_SetItemString(outDict, "stk_1",  PyInt_FromLong((long)lTab->stk_1));
  PyDict_SetItemString(outDict, "no_chan",  PyInt_FromLong((long)lTab->no_chan));
  PyDict_SetItemString(outDict, "ref_freq",  PyFloat_FromDouble((double)lTab->ref_freq));
  PyDict_SetItemString(outDict, "chan_bw",  PyFloat_FromDouble((double)lTab->chan_bw));
  PyDict_SetItemString(outDict, "ref_pixl",  PyInt_FromLong((long)lTab->ref_pixl));
  PyDict_SetItemString(outDict, "obscode", PyString_InternFromString(lTab->obscode));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));
  PyDict_SetItemString(outDict, "frame", PyString_InternFromString(lTab->frame));
  PyDict_SetItemString(outDict, "Freq",  PyFloat_FromDouble((double)lTab->Freq));
  PyDict_SetItemString(outDict, "TimeSys", PyString_InternFromString(lTab->TimeSys));
  PyDict_SetItemString(outDict, "GSTiat0",  PyFloat_FromDouble((double)lTab->GSTiat0));
  PyDict_SetItemString(outDict, "DegDay",  PyFloat_FromDouble((double)lTab->DegDay));
  PyDict_SetItemString(outDict, "iatUtc",  PyFloat_FromDouble((double)lTab->iatUtc));
  PyDict_SetItemString(outDict, "PolarX",  PyFloat_FromDouble((double)lTab->PolarX));
  PyDict_SetItemString(outDict, "PolarY",  PyFloat_FromDouble((double)lTab->PolarY));
  PyDict_SetItemString(outDict, "ArrayX",  PyFloat_FromDouble((double)lTab->ArrayX));
  PyDict_SetItemString(outDict, "ArrayY",  PyFloat_FromDouble((double)lTab->ArrayY));
  PyDict_SetItemString(outDict, "ArrayZ",  PyFloat_FromDouble((double)lTab->ArrayZ));

  return outDict;
} 

extern void TableIDI_ARRAY_GEOMETRYSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableIDI_ARRAY_GEOMETRY *lTab = (ObitTableIDI_ARRAY_GEOMETRY*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEIDI_ARRAY_GEOMETRY;

  lTab->tabrev = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tabrev"));
  lTab->no_stkd = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "no_stkd"));
  lTab->stk_1 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "stk_1"));
  lTab->no_chan = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "no_chan"));
  lTab->ref_freq = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ref_freq"));
  lTab->chan_bw = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "chan_bw"));
  lTab->ref_pixl = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "ref_pixl"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obscode"));
  strncpy (lTab->obscode, tstr, lstr); lTab->obscode[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "frame"));
  strncpy (lTab->frame, tstr, lstr); lTab->frame[lstr-1]=0;
  lTab->Freq = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "Freq"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "TimeSys"));
  strncpy (lTab->TimeSys, tstr, lstr); lTab->TimeSys[lstr-1]=0;
  lTab->GSTiat0 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "GSTiat0"));
  lTab->DegDay = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "DegDay"));
  lTab->iatUtc = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "iatUtc"));
  lTab->PolarX = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "PolarX"));
  lTab->PolarY = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "PolarY"));
  lTab->ArrayX = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ArrayX"));
  lTab->ArrayY = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ArrayY"));
  lTab->ArrayZ = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ArrayZ"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableIDI_FREQUENCY.h"
%}
 
%inline %{
 
extern ObitTable* TableIDI_FREQUENCY (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int no_band,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lno_band = (oint)no_band;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableIDI_FREQUENCYValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lno_band,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableIDI_FREQUENCYGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableIDI_FREQUENCY *lTab = (ObitTableIDI_FREQUENCY*)inTab;
  PyDict_SetItemString(outDict, "no_band",  PyInt_FromLong((long)lTab->no_band));
  PyDict_SetItemString(outDict, "tabrev",  PyInt_FromLong((long)lTab->tabrev));
  PyDict_SetItemString(outDict, "no_stkd",  PyInt_FromLong((long)lTab->no_stkd));
  PyDict_SetItemString(outDict, "stk_1",  PyInt_FromLong((long)lTab->stk_1));
  PyDict_SetItemString(outDict, "no_chan",  PyInt_FromLong((long)lTab->no_chan));
  PyDict_SetItemString(outDict, "ref_freq",  PyFloat_FromDouble((double)lTab->ref_freq));
  PyDict_SetItemString(outDict, "chan_bw",  PyFloat_FromDouble((double)lTab->chan_bw));
  PyDict_SetItemString(outDict, "ref_pixl",  PyInt_FromLong((long)lTab->ref_pixl));
  PyDict_SetItemString(outDict, "obscode", PyString_InternFromString(lTab->obscode));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));

  return outDict;
} 

extern void TableIDI_FREQUENCYSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableIDI_FREQUENCY *lTab = (ObitTableIDI_FREQUENCY*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEIDI_FREQUENCY;

  lTab->tabrev = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tabrev"));
  lTab->no_stkd = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "no_stkd"));
  lTab->stk_1 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "stk_1"));
  lTab->no_chan = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "no_chan"));
  lTab->ref_freq = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ref_freq"));
  lTab->chan_bw = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "chan_bw"));
  lTab->ref_pixl = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "ref_pixl"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obscode"));
  strncpy (lTab->obscode, tstr, lstr); lTab->obscode[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableIDI_SOURCE.h"
%}
 
%inline %{
 
extern ObitTable* TableIDI_SOURCE (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int no_band,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lno_band = (oint)no_band;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableIDI_SOURCEValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lno_band,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableIDI_SOURCEGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableIDI_SOURCE *lTab = (ObitTableIDI_SOURCE*)inTab;
  PyDict_SetItemString(outDict, "no_band",  PyInt_FromLong((long)lTab->no_band));
  PyDict_SetItemString(outDict, "tabrev",  PyInt_FromLong((long)lTab->tabrev));
  PyDict_SetItemString(outDict, "no_stkd",  PyInt_FromLong((long)lTab->no_stkd));
  PyDict_SetItemString(outDict, "stk_1",  PyInt_FromLong((long)lTab->stk_1));
  PyDict_SetItemString(outDict, "no_chan",  PyInt_FromLong((long)lTab->no_chan));
  PyDict_SetItemString(outDict, "ref_freq",  PyFloat_FromDouble((double)lTab->ref_freq));
  PyDict_SetItemString(outDict, "chan_bw",  PyFloat_FromDouble((double)lTab->chan_bw));
  PyDict_SetItemString(outDict, "ref_pixl",  PyInt_FromLong((long)lTab->ref_pixl));
  PyDict_SetItemString(outDict, "obscode", PyString_InternFromString(lTab->obscode));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));

  return outDict;
} 

extern void TableIDI_SOURCESetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableIDI_SOURCE *lTab = (ObitTableIDI_SOURCE*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEIDI_SOURCE;

  lTab->tabrev = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tabrev"));
  lTab->no_stkd = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "no_stkd"));
  lTab->stk_1 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "stk_1"));
  lTab->no_chan = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "no_chan"));
  lTab->ref_freq = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ref_freq"));
  lTab->chan_bw = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "chan_bw"));
  lTab->ref_pixl = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "ref_pixl"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obscode"));
  strncpy (lTab->obscode, tstr, lstr); lTab->obscode[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableIDI_UV_DATA.h"
%}
 
%inline %{
 
extern ObitTable* TableIDI_UV_DATA (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int no_band, int maxis1, int maxis2, int maxis3, int maxis4, int maxis5,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lno_band = (oint)no_band;
   oint lmaxis1 = (oint)maxis1;
   oint lmaxis2 = (oint)maxis2;
   oint lmaxis3 = (oint)maxis3;
   oint lmaxis4 = (oint)maxis4;
   oint lmaxis5 = (oint)maxis5;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableIDI_UV_DATAValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lno_band, lmaxis1, lmaxis2, lmaxis3, lmaxis4, lmaxis5,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableIDI_UV_DATAGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableIDI_UV_DATA *lTab = (ObitTableIDI_UV_DATA*)inTab;
  PyDict_SetItemString(outDict, "no_band",  PyInt_FromLong((long)lTab->no_band));
  PyDict_SetItemString(outDict, "maxis1",  PyInt_FromLong((long)lTab->maxis1));
  PyDict_SetItemString(outDict, "maxis2",  PyInt_FromLong((long)lTab->maxis2));
  PyDict_SetItemString(outDict, "maxis3",  PyInt_FromLong((long)lTab->maxis3));
  PyDict_SetItemString(outDict, "maxis4",  PyInt_FromLong((long)lTab->maxis4));
  PyDict_SetItemString(outDict, "maxis5",  PyInt_FromLong((long)lTab->maxis5));
  PyDict_SetItemString(outDict, "tabrev",  PyInt_FromLong((long)lTab->tabrev));
  PyDict_SetItemString(outDict, "no_stkd",  PyInt_FromLong((long)lTab->no_stkd));
  PyDict_SetItemString(outDict, "stk_1",  PyInt_FromLong((long)lTab->stk_1));
  PyDict_SetItemString(outDict, "no_chan",  PyInt_FromLong((long)lTab->no_chan));
  PyDict_SetItemString(outDict, "ref_freq",  PyFloat_FromDouble((double)lTab->ref_freq));
  PyDict_SetItemString(outDict, "chan_bw",  PyFloat_FromDouble((double)lTab->chan_bw));
  PyDict_SetItemString(outDict, "ref_pixl",  PyInt_FromLong((long)lTab->ref_pixl));
  PyDict_SetItemString(outDict, "obscode", PyString_InternFromString(lTab->obscode));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));
  PyDict_SetItemString(outDict, "nmatrix",  PyInt_FromLong((long)lTab->nmatrix));
  PyDict_SetItemString(outDict, "maxis",  PyInt_FromLong((long)lTab->maxis));
 PyDict_SetItemString(outDict, "tmatx01",  PyInt_FromLong((long)lTab->tmatx01));
 PyDict_SetItemString(outDict, "tmatx02",  PyInt_FromLong((long)lTab->tmatx02));
 PyDict_SetItemString(outDict, "tmatx03",  PyInt_FromLong((long)lTab->tmatx03));
 PyDict_SetItemString(outDict, "tmatx04",  PyInt_FromLong((long)lTab->tmatx04));
 PyDict_SetItemString(outDict, "tmatx05",  PyInt_FromLong((long)lTab->tmatx05));
 PyDict_SetItemString(outDict, "tmatx06",  PyInt_FromLong((long)lTab->tmatx06));
 PyDict_SetItemString(outDict, "tmatx07",  PyInt_FromLong((long)lTab->tmatx07));
 PyDict_SetItemString(outDict, "tmatx08",  PyInt_FromLong((long)lTab->tmatx08));
 PyDict_SetItemString(outDict, "tmatx09",  PyInt_FromLong((long)lTab->tmatx09));
 PyDict_SetItemString(outDict, "tmatx10",  PyInt_FromLong((long)lTab->tmatx10));
 PyDict_SetItemString(outDict, "tmatx11",  PyInt_FromLong((long)lTab->tmatx11));
 PyDict_SetItemString(outDict, "tmatx12",  PyInt_FromLong((long)lTab->tmatx12));
 PyDict_SetItemString(outDict, "tmatx13",  PyInt_FromLong((long)lTab->tmatx13));
 PyDict_SetItemString(outDict, "tmatx14",  PyInt_FromLong((long)lTab->tmatx14));
 PyDict_SetItemString(outDict, "tmatx15",  PyInt_FromLong((long)lTab->tmatx15));
  PyDict_SetItemString(outDict, "cdelt1",  PyFloat_FromDouble((double)lTab->cdelt1));
  PyDict_SetItemString(outDict, "crpix1",  PyFloat_FromDouble((double)lTab->crpix1));
  PyDict_SetItemString(outDict, "crval1",  PyFloat_FromDouble((double)lTab->crval1));
  PyDict_SetItemString(outDict, "cdelt2",  PyFloat_FromDouble((double)lTab->cdelt2));
  PyDict_SetItemString(outDict, "crpix2",  PyFloat_FromDouble((double)lTab->crpix2));
  PyDict_SetItemString(outDict, "crval2",  PyFloat_FromDouble((double)lTab->crval2));
  PyDict_SetItemString(outDict, "cdelt3",  PyFloat_FromDouble((double)lTab->cdelt3));
  PyDict_SetItemString(outDict, "crpix3",  PyFloat_FromDouble((double)lTab->crpix3));
  PyDict_SetItemString(outDict, "crval3",  PyFloat_FromDouble((double)lTab->crval3));
  PyDict_SetItemString(outDict, "cdelt4",  PyFloat_FromDouble((double)lTab->cdelt4));
  PyDict_SetItemString(outDict, "crpix4",  PyFloat_FromDouble((double)lTab->crpix4));
  PyDict_SetItemString(outDict, "crval4",  PyFloat_FromDouble((double)lTab->crval4));
  PyDict_SetItemString(outDict, "cdelt5",  PyFloat_FromDouble((double)lTab->cdelt5));
  PyDict_SetItemString(outDict, "crpix5",  PyFloat_FromDouble((double)lTab->crpix5));
  PyDict_SetItemString(outDict, "crval5",  PyFloat_FromDouble((double)lTab->crval5));
  PyDict_SetItemString(outDict, "maxis6",  PyInt_FromLong((long)lTab->maxis6));
  PyDict_SetItemString(outDict, "cdelt6",  PyFloat_FromDouble((double)lTab->cdelt6));
  PyDict_SetItemString(outDict, "crpix6",  PyFloat_FromDouble((double)lTab->crpix6));
  PyDict_SetItemString(outDict, "crval6",  PyFloat_FromDouble((double)lTab->crval6));
  PyDict_SetItemString(outDict, "maxis7",  PyInt_FromLong((long)lTab->maxis7));
  PyDict_SetItemString(outDict, "cdelt7",  PyFloat_FromDouble((double)lTab->cdelt7));
  PyDict_SetItemString(outDict, "crpix7",  PyFloat_FromDouble((double)lTab->crpix7));
  PyDict_SetItemString(outDict, "crval7",  PyFloat_FromDouble((double)lTab->crval7));

  return outDict;
} 

extern void TableIDI_UV_DATASetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableIDI_UV_DATA *lTab = (ObitTableIDI_UV_DATA*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEIDI_UV_DATA;

  lTab->tabrev = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tabrev"));
  lTab->no_stkd = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "no_stkd"));
  lTab->stk_1 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "stk_1"));
  lTab->no_chan = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "no_chan"));
  lTab->ref_freq = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "ref_freq"));
  lTab->chan_bw = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "chan_bw"));
  lTab->ref_pixl = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "ref_pixl"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obscode"));
  strncpy (lTab->obscode, tstr, lstr); lTab->obscode[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;
  lTab->nmatrix = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "nmatrix"));
  lTab->maxis = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "maxis"));
  lTab->tmatx01 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx01"));
  lTab->tmatx02 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx02"));
  lTab->tmatx03 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx03"));
  lTab->tmatx04 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx04"));
  lTab->tmatx05 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx05"));
  lTab->tmatx06 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx06"));
  lTab->tmatx07 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx07"));
  lTab->tmatx08 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx08"));
  lTab->tmatx09 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx09"));
  lTab->tmatx10 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx10"));
  lTab->tmatx11 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx11"));
  lTab->tmatx12 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx12"));
  lTab->tmatx13 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx13"));
  lTab->tmatx14 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx14"));
  lTab->tmatx15 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tmatx15"));
  lTab->cdelt1 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "cdelt1"));
  lTab->crpix1 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crpix1"));
  lTab->crval1 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crval1"));
  lTab->cdelt2 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "cdelt2"));
  lTab->crpix2 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crpix2"));
  lTab->crval2 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crval2"));
  lTab->cdelt3 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "cdelt3"));
  lTab->crpix3 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crpix3"));
  lTab->crval3 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crval3"));
  lTab->cdelt4 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "cdelt4"));
  lTab->crpix4 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crpix4"));
  lTab->crval4 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crval4"));
  lTab->cdelt5 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "cdelt5"));
  lTab->crpix5 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crpix5"));
  lTab->crval5 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crval5"));
  lTab->maxis6 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "maxis6"));
  lTab->cdelt6 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "cdelt6"));
  lTab->crpix6 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crpix6"));
  lTab->crval6 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crval6"));
  lTab->maxis7 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "maxis7"));
  lTab->cdelt7 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "cdelt7"));
  lTab->crpix7 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crpix7"));
  lTab->crval7 = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "crval7"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableIM.h"
%}
 
%inline %{
 
extern ObitTable* TableIM (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numBand, int numPol, int npoly,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumBand = (oint)numBand;
   oint lnumPol = (oint)numPol;
   oint lnpoly = (oint)npoly;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableIMValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumBand, lnumPol, lnpoly,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableIMGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableIM *lTab = (ObitTableIM*)inTab;
  PyDict_SetItemString(outDict, "numBand",  PyInt_FromLong((long)lTab->numBand));
  PyDict_SetItemString(outDict, "numPol",  PyInt_FromLong((long)lTab->numPol));
  PyDict_SetItemString(outDict, "npoly",  PyInt_FromLong((long)lTab->npoly));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));
  PyDict_SetItemString(outDict, "obscode", PyString_InternFromString(lTab->obscode));
  PyDict_SetItemString(outDict, "tabrev",  PyInt_FromLong((long)lTab->tabrev));
  PyDict_SetItemString(outDict, "revision",  PyFloat_FromDouble((double)lTab->revision));

  return outDict;
} 

extern void TableIMSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableIM *lTab = (ObitTableIM*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEIM;

  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obscode"));
  strncpy (lTab->obscode, tstr, lstr); lTab->obscode[lstr-1]=0;
  lTab->tabrev = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "tabrev"));
  lTab->revision = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "revision"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id: Table.inc,v 1.16 2007/12/16 20:25:12 bcotton Exp $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for Table type                             */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
  gchar    *cdata, *ctemp, *tstr;
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
      for (j=0; j<desc->repeat[i]; j++) 
        idata[j] = (gshort)PyInt_AsLong(PyList_GetItem(list, j));
      break;
    case OBIT_int:
      bad = desc->repeat[i] != PyList_Size(list);
      jdata = ((olong*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        jdata[j] = (olong)PyInt_AsLong(PyList_GetItem(list, j));
      break;
    case OBIT_oint:
      bad = desc->repeat[i] != PyList_Size(list);
      kdata = ((oint*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        kdata[j] = (oint)PyInt_AsLong(PyList_GetItem(list, j));
      break;
    case OBIT_long:
      bad = desc->repeat[i] != PyList_Size(list);
      ldata = ((olong*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        ldata[j] = (olong)PyInt_AsLong(PyList_GetItem(list, j));
      break;
    case OBIT_float:
      bad = desc->repeat[i] != PyList_Size(list);
      fdata = ((ofloat*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        fdata[j] =  (ofloat)PyFloat_AsDouble(PyList_GetItem(list, j));
      break;
    case OBIT_double:
      bad = desc->repeat[i] != PyList_Size(list);
      ddata = ((odouble*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        ddata[j] =  (odouble)PyFloat_AsDouble(PyList_GetItem(list, j));
      break;
    case OBIT_string:
      bad = (desc->repeat[i]/MAX (1,desc->dim[i][0])) != PyList_Size(list);
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
      bdata = ((gboolean*)row->myRowData)+desc->offset[i];
      for (j=0; j<desc->repeat[i]; j++) 
        bdata[j] = (gboolean)PyInt_AsLong(PyList_GetItem(list, j));
      break;
    case OBIT_bits:
      bad = desc->repeat[i] != PyList_Size(list);
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
		   "%s: wrong size %d %d", 
		   routine, desc->repeat[i], PyList_Size(list));
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

/* $Id: TableList.inc,v 1.4 2007/07/26 14:28:24 bcotton Exp $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageDesc type                         */
/*                                                                    */
/*;  Copyright (C) 2005,2008                                          */
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
#include "ObitTableList.h"
%}

%inline %{
extern ObitTableList* TableListCreate (char *name) {
  return newObitTableList (name);
} // end TableListCreate

extern ObitTableList* TableListCopy (ObitTableList* in, 
		              ObitTableList* out, ObitErr *err) {
  return ObitTableListCopy (in, out, err);
} // end TableListCopy

extern PyObject *TableListGetList(ObitTableList* in, ObitErr *err) {
  PyObject *outList=NULL;
  PyObject *list=NULL;
  ObitTable *table=NULL;
  gchar *tabName=NULL;
  olong tabVer, i;

  if  (err->error) return outList;

  //ObitTableListPrint(in, err);
  outList = PyList_New(in->number);
  for (i=0; i<in->number; i++) {
    ObitTableListGetNumber(in, i+1, &tabName, &tabVer, &table, err);
    list = PyList_New(2);
    PyList_SetItem(list, 0, PyInt_FromLong((long)tabVer));
    PyList_SetItem(list, 1, PyString_InternFromString(tabName));
    PyList_SetItem(outList, i, list);
    if (tabName) g_free(tabName);
    ObitTableUnref(table);
  }
  return outList;
} // end TableListGetList

long TableListGetHigh (ObitTableList* in, char *tabType) {
  return ObitTableListGetHigh (in, tabType);
} // end TableListGetHigh

void TableListPutHi (ObitTableList* in, ObitErr *err) {
  olong lversion = 1;
  gchar *tabType = "AIPS HI";
  ObitTableListPut (in, tabType, &lversion, NULL, err);
} // end TableListPutHi

ObitTableList* TableListRef (ObitTableList* in) {
  return ObitTableListRef (in);
} // end TableListRef

ObitTableList* TableListUnref (ObitTableList* in) {
  if (!ObitTableListIsA(in)) return NULL;
  return ObitTableListUnref (in);
} // end TableListUnref

extern int TableListIsA (ObitTableList* in) {
  return ObitTableListIsA(in);
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitTableList *me;
} TableList;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitTableList *me;
} TableList;

%addmethods TableList { 
  TableList(char *name) {
     TableList *out;
     out = (TableList *) malloc(sizeof(TableList));
     if (strcmp(name, "None")) out->me = TableListCreate (name);
     else out->me = NULL;
     return out;
   }
  ~TableList() {
    self->me = TableListUnref(self->me);
    free(self);
  }
};

/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableMC.h"
%}
 
%inline %{
 
extern ObitTable* TableMC (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numBand, int numPol,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumBand = (oint)numBand;
   oint lnumPol = (oint)numPol;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableMCValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumBand, lnumPol,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableMCGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableMC *lTab = (ObitTableMC*)inTab;
  PyDict_SetItemString(outDict, "numBand",  PyInt_FromLong((long)lTab->numBand));
  PyDict_SetItemString(outDict, "numPol",  PyInt_FromLong((long)lTab->numPol));
  PyDict_SetItemString(outDict, "obscode", PyString_InternFromString(lTab->obscode));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));

  return outDict;
} 

extern void TableMCSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableMC *lTab = (ObitTableMC*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEMC;

  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obscode"));
  strncpy (lTab->obscode, tstr, lstr); lTab->obscode[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;
  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableMF.h"
%}
 
%inline %{
 
extern ObitTable* TableMF (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableMFValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableMFGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableMF *lTab = (ObitTableMF*)inTab;
  PyDict_SetItemString(outDict, "depth1",  PyFloat_FromDouble((double)lTab->depth1));
  PyDict_SetItemString(outDict, "depth2",  PyFloat_FromDouble((double)lTab->depth2));
  PyDict_SetItemString(outDict, "depth3",  PyFloat_FromDouble((double)lTab->depth3));
  PyDict_SetItemString(outDict, "depth4",  PyFloat_FromDouble((double)lTab->depth4));
  PyDict_SetItemString(outDict, "depth5",  PyFloat_FromDouble((double)lTab->depth5));

  return outDict;
} 

extern void TableMFSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableMF *lTab = (ObitTableMF*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEMF;

  lTab->depth1 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "depth1"));
  lTab->depth2 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "depth2"));
  lTab->depth3 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "depth3"));
  lTab->depth4 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "depth4"));
  lTab->depth5 = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "depth5"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableNI.h"
%}
 
%inline %{
 
extern ObitTable* TableNI (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numCoef,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumCoef = (oint)numCoef;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableNIValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumCoef,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableNIGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableNI *lTab = (ObitTableNI*)inTab;
  PyDict_SetItemString(outDict, "numCoef",  PyInt_FromLong((long)lTab->numCoef));
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));
  PyDict_SetItemString(outDict, "heightIon",  PyFloat_FromDouble((double)lTab->heightIon));

  return outDict;
} 

extern void TableNISetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableNI *lTab = (ObitTableNI*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLENI;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));
  lTab->heightIon = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, "heightIon"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableNX.h"
%}
 
%inline %{
 
extern ObitTable* TableNX (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableNXValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableNXGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableNX *lTab = (ObitTableNX*)inTab;

  return outDict;
} 

extern void TableNXSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableNX *lTab = (ObitTableNX*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLENX;


  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableOB.h"
%}
 
%inline %{
 
extern ObitTable* TableOB (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableOBValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableOBGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableOB *lTab = (ObitTableOB*)inTab;
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));

  return outDict;
} 

extern void TableOBSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableOB *lTab = (ObitTableOB*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEOB;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableOF.h"
%}
 
%inline %{
 
extern ObitTable* TableOF (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableOFValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableOFGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableOF *lTab = (ObitTableOF*)inTab;
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));
  PyDict_SetItemString(outDict, "obscode", PyString_InternFromString(lTab->obscode));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));

  return outDict;
} 

extern void TableOFSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableOF *lTab = (ObitTableOF*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEOF;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obscode"));
  strncpy (lTab->obscode, tstr, lstr); lTab->obscode[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTablePC.h"
%}
 
%inline %{
 
extern ObitTable* TablePC (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numPol, int numBand, int numTones,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumPol = (oint)numPol;
   oint lnumBand = (oint)numBand;
   oint lnumTones = (oint)numTones;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTablePCValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumPol, lnumBand, lnumTones,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TablePCGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTablePC *lTab = (ObitTablePC*)inTab;
  PyDict_SetItemString(outDict, "numPol",  PyInt_FromLong((long)lTab->numPol));
  PyDict_SetItemString(outDict, "numBand",  PyInt_FromLong((long)lTab->numBand));
  PyDict_SetItemString(outDict, "numTones",  PyInt_FromLong((long)lTab->numTones));
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));

  return outDict;
} 

extern void TablePCSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTablePC *lTab = (ObitTablePC*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEPC;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTablePS.h"
%}
 
%inline %{
 
extern ObitTable* TablePS (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTablePSValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TablePSGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTablePS *lTab = (ObitTablePS*)inTab;
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));

  return outDict;
} 

extern void TablePSSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTablePS *lTab = (ObitTablePS*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEPS;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableSN.h"
%}
 
%inline %{
 
extern ObitTable* TableSN (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numPol, int numIF,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumPol = (oint)numPol;
   oint lnumIF = (oint)numIF;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableSNValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumPol, lnumIF,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableSNGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableSN *lTab = (ObitTableSN*)inTab;
  PyDict_SetItemString(outDict, "numPol",  PyInt_FromLong((long)lTab->numPol));
  PyDict_SetItemString(outDict, "numIF",  PyInt_FromLong((long)lTab->numIF));
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));
  PyDict_SetItemString(outDict, "numNodes",  PyInt_FromLong((long)lTab->numNodes));
  PyDict_SetItemString(outDict, "mGMod",  PyFloat_FromDouble((double)lTab->mGMod));
 PyDict_SetItemString(outDict, "isApplied",  PyInt_FromLong((long)lTab->isApplied));

  return outDict;
} 

extern void TableSNSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableSN *lTab = (ObitTableSN*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLESN;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));
  lTab->numNodes = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "numNodes"));
  lTab->mGMod = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "mGMod"));
  lTab->isApplied = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "isApplied"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableSU.h"
%}
 
%inline %{
 
extern ObitTable* TableSU (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numIF,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumIF = (oint)numIF;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableSUValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumIF,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableSUGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableSU *lTab = (ObitTableSU*)inTab;
  PyDict_SetItemString(outDict, "numIF",  PyInt_FromLong((long)lTab->numIF));
  PyDict_SetItemString(outDict, "velType", PyString_InternFromString(lTab->velType));
  PyDict_SetItemString(outDict, "velDef", PyString_InternFromString(lTab->velDef));
  PyDict_SetItemString(outDict, "FreqID",  PyInt_FromLong((long)lTab->FreqID));

  return outDict;
} 

extern void TableSUSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableSU *lTab = (ObitTableSU*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLESU;

  tstr = PyString_AsString(PyDict_GetItemString(inDict, "velType"));
  strncpy (lTab->velType, tstr, lstr); lTab->velType[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "velDef"));
  strncpy (lTab->velDef, tstr, lstr); lTab->velDef[lstr-1]=0;
  lTab->FreqID = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "FreqID"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableTY.h"
%}
 
%inline %{
 
extern ObitTable* TableTY (ObitData *inData, long *tabVer,
 	                   int access,
 	                   char *tabName,
                           int numPol, int numIF,
                           ObitErr *err)
 {
   ObitIOAccess laccess;
   /* Cast structural keywords to correct type */
   oint lnumPol = (oint)numPol;
   oint lnumIF = (oint)numIF;
   olong ltabVer = (olong)*tabVer;
   ObitTable *outTable=NULL;
   laccess = OBIT_IO_ReadOnly;
   if (access==2) laccess = OBIT_IO_WriteOnly;
   else if (access==3) laccess = OBIT_IO_ReadWrite;
   outTable = (ObitTable*)newObitTableTYValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                           lnumPol, lnumIF,
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableTYGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableTY *lTab = (ObitTableTY*)inTab;
  PyDict_SetItemString(outDict, "numPol",  PyInt_FromLong((long)lTab->numPol));
  PyDict_SetItemString(outDict, "numIF",  PyInt_FromLong((long)lTab->numIF));
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));

  return outDict;
} 

extern void TableTYSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableTY *lTab = (ObitTableTY*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLETY;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id: TableUtil.inc,v 1.2 2007/08/23 20:49:03 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for Obit Table Utilities                   */
/*                                                                    */
/*;  Copyright (C) 2006                                               */
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

#include "ObitTableCCUtil.h"
%}

%inline %{
olong TableCCUtilMerge (ObitTable *in, ObitTable *out, ObitErr *err)
{
  ObitIOCode ret=1;
  gchar *routine="TableCCUtilMerge";

  /* Check that input TableCCs */
  Obit_retval_if_fail((ObitTableCCIsA(in)), err, ret,
		      "%s: input %s NOT an AIPS CC table", routine, in->name);
  Obit_retval_if_fail((ObitTableCCIsA(out)), err, ret,
		      "%s: output %s NOT an AIPS CC table", routine, out->name);

  ret = ObitTableCCUtilMerge ((ObitTableCC*)in, (ObitTableCC*)out, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end TableCCUtilMerge
%}


/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableVL.h"
%}
 
%inline %{
 
extern ObitTable* TableVL (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableVLValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableVLGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableVL *lTab = (ObitTableVL*)inTab;
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));
  PyDict_SetItemString(outDict, "numIndexed",  PyInt_FromLong((long)lTab->numIndexed));
  PyDict_SetItemString(outDict, "index00",  PyInt_FromLong((long)lTab->index00));
  PyDict_SetItemString(outDict, "index01",  PyInt_FromLong((long)lTab->index01));
  PyDict_SetItemString(outDict, "index03",  PyInt_FromLong((long)lTab->index03));
  PyDict_SetItemString(outDict, "index04",  PyInt_FromLong((long)lTab->index04));
  PyDict_SetItemString(outDict, "index05",  PyInt_FromLong((long)lTab->index05));
  PyDict_SetItemString(outDict, "index06",  PyInt_FromLong((long)lTab->index06));
  PyDict_SetItemString(outDict, "index07",  PyInt_FromLong((long)lTab->index07));
  PyDict_SetItemString(outDict, "index08",  PyInt_FromLong((long)lTab->index08));
  PyDict_SetItemString(outDict, "index09",  PyInt_FromLong((long)lTab->index09));
  PyDict_SetItemString(outDict, "index10",  PyInt_FromLong((long)lTab->index10));
  PyDict_SetItemString(outDict, "index11",  PyInt_FromLong((long)lTab->index11));
  PyDict_SetItemString(outDict, "index12",  PyInt_FromLong((long)lTab->index12));
  PyDict_SetItemString(outDict, "index13",  PyInt_FromLong((long)lTab->index13));
  PyDict_SetItemString(outDict, "index14",  PyInt_FromLong((long)lTab->index14));
  PyDict_SetItemString(outDict, "index15",  PyInt_FromLong((long)lTab->index15));
  PyDict_SetItemString(outDict, "index16",  PyInt_FromLong((long)lTab->index16));
  PyDict_SetItemString(outDict, "index17",  PyInt_FromLong((long)lTab->index17));
  PyDict_SetItemString(outDict, "index18",  PyInt_FromLong((long)lTab->index18));
  PyDict_SetItemString(outDict, "index19",  PyInt_FromLong((long)lTab->index19));
  PyDict_SetItemString(outDict, "index20",  PyInt_FromLong((long)lTab->index20));
  PyDict_SetItemString(outDict, "index21",  PyInt_FromLong((long)lTab->index21));
  PyDict_SetItemString(outDict, "index22",  PyInt_FromLong((long)lTab->index22));
  PyDict_SetItemString(outDict, "index23",  PyInt_FromLong((long)lTab->index23));

  return outDict;
} 

extern void TableVLSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableVL *lTab = (ObitTableVL*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEVL;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));
  lTab->numIndexed = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "numIndexed"));
  lTab->index00 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index00"));
  lTab->index01 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index01"));
  lTab->index03 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index03"));
  lTab->index04 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index04"));
  lTab->index05 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index05"));
  lTab->index06 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index06"));
  lTab->index07 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index07"));
  lTab->index08 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index08"));
  lTab->index09 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index09"));
  lTab->index10 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index10"));
  lTab->index11 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index11"));
  lTab->index12 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index12"));
  lTab->index13 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index13"));
  lTab->index14 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index14"));
  lTab->index15 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index15"));
  lTab->index16 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index16"));
  lTab->index17 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index17"));
  lTab->index18 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index18"));
  lTab->index19 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index19"));
  lTab->index20 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index20"));
  lTab->index21 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index21"));
  lTab->index22 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index22"));
  lTab->index23 = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "index23"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableVZ.h"
%}
 
%inline %{
 
extern ObitTable* TableVZ (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableVZValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableVZGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableVZ *lTab = (ObitTableVZ*)inTab;
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));
  PyDict_SetItemString(outDict, "refFreq",  PyFloat_FromDouble((double)lTab->refFreq));

  return outDict;
} 

extern void TableVZSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableVZ *lTab = (ObitTableVZ*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEVZ;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));
  lTab->refFreq = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, "refFreq"));

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2008                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
%{
#include "Obit.h"
#include "ObitData.h"
#include "ObitTableWX.h"
%}
 
%inline %{
 
extern ObitTable* TableWX (ObitData *inData, long *tabVer,
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
   outTable = (ObitTable*)newObitTableWXValue ((gchar*)tabName, inData, (olong*)&ltabVer,
   			   laccess, 
                          
                           err);
   *tabVer = (long)ltabVer;
   return outTable;
   }
 
extern PyObject* TableWXGetHeadKeys (ObitTable *inTab) {
  PyObject *outDict=PyDict_New();
  ObitTableWX *lTab = (ObitTableWX*)inTab;
  PyDict_SetItemString(outDict, "revision",  PyInt_FromLong((long)lTab->revision));
  PyDict_SetItemString(outDict, "obscode", PyString_InternFromString(lTab->obscode));
  PyDict_SetItemString(outDict, "RefDate", PyString_InternFromString(lTab->RefDate));

  return outDict;
} 

extern void TableWXSetHeadKeys (ObitTable *inTab, PyObject *inDict) {
  ObitTableWX *lTab = (ObitTableWX*)inTab;
  char *tstr;
  int lstr=MAXKEYCHARTABLEWX;

  lTab->revision = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, "revision"));
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obscode"));
  strncpy (lTab->obscode, tstr, lstr); lTab->obscode[lstr-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "RefDate"));
  strncpy (lTab->RefDate, tstr, lstr); lTab->RefDate[lstr-1]=0;

  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) 
    lTab->myStatus = OBIT_Modified;
} 

%}
/* $Id: TimeFilter.inc,v 1.1 2008/01/29 02:24:33 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitTimeFilter type                    */
/*                                                                    */
/*;  Copyright (C) 2008                                               */
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
#include "ObitTimeFilter.h"
%}


%inline %{

/** Public: Create time filter. */
extern ObitTimeFilter* TimeFilterCreate(char* name, int nTime, int nSeries) {
   return  newObitTimeFilter ((gchar*)name, (olong)nTime, (olong)nSeries);
}

/** Public: Resize filter. */
extern void TimeFilterResize (ObitTimeFilter* in, int nTime) {
  ObitTimeFilterResize (in, (olong)nTime);
} // end TimeFilterResize

/** Public: Construct regular time series. */
extern void TimeFilterGridTime (ObitTimeFilter* in, int seriesNo, float dTime, int nTime, 
                        float *times, float *data) {
  ObitTimeFilterGridTime (in, (olong)seriesNo, (ofloat)dTime, (olong)nTime, 
                        (ofloat*)times, (ofloat*)data);
} // end TimeFilterGridTime

/** Public: Return time series. */
extern PyObject* TimeFilterUngridTime (ObitTimeFilter* in, int seriesNo,
			         int nTime, float *times) {
  olong i;
  PyObject *out;
  ofloat *tdata=NULL;

  tdata = g_malloc0(nTime*sizeof(ofloat));
  ObitTimeFilterUngridTime (in, (olong)seriesNo, (olong)nTime, 
                            (ofloat*)times, tdata);
  out = PyList_New(nTime);
  for (i=0; i<nTime; i++) 
    PyList_SetItem(out, i, PyFloat_FromDouble((double)tdata[i]));
  if (tdata) g_free(tdata);
  return out;
} // end TimeFilterUngridTime

/** Public: Compute frequency series. */
extern void TimeFilter2Freq (ObitTimeFilter* in) {
  ObitTimeFilter2Freq (in);
} // end TimeFilter2Freq

/** Public: Compute Time series. */
extern void TimeFilter2Time (ObitTimeFilter* in) {
  ObitTimeFilter2Time (in);
} // end TimeFilter2Time

/** Public: Apply Filter to Frequency series. */
extern void TimeFilterFilter (ObitTimeFilter* in, int seriesNo, int type, 
                              float *parms, ObitErr *err) {
  ObitTimeFilterType ltype=OBIT_TimeFilter_LowPass;
  if (type==1) ltype = OBIT_TimeFilter_HighPass;
  if (type==2) ltype = OBIT_TimeFilter_NotchPass;
  if (type==3) ltype = OBIT_TimeFilter_NotchBlock;

  ObitTimeFilterFilter (in, (olong)seriesNo, ltype, (ofloat*)parms, err);
} // end TimeFilterFilter

/** Public: Apply Filter to Frequency series with physical parameters. */
extern void TimeFilterDoFilter (ObitTimeFilter* in, int seriesNo, int type, 
                                float *parms, ObitErr *err) {
  ObitTimeFilterType ltype=OBIT_TimeFilter_LowPass;
  if (type==1) ltype = OBIT_TimeFilter_HighPass;
  if (type==2) ltype = OBIT_TimeFilter_NotchPass;
  if (type==3) ltype = OBIT_TimeFilter_NotchBlock;

  ObitTimeFilterDoFilter (in, (olong)seriesNo, ltype, (ofloat*)parms, err);
} // end TimeFilterDoFilter

/** Public:  Plot power spectrum. */
extern void TimeFilterPlotPower (ObitTimeFilter* in, int seriesNo, char *label, ObitErr *err) {
  ObitTimeFilterPlotPower (in, (olong)seriesNo, (gchar*)label, err);
} // end TimeFilter

/** Public: Plot Time series.. */
extern void TimeFilterPlotTime (ObitTimeFilter* in, int seriesNo, char *label, ObitErr *err) {
  ObitTimeFilterPlotTime (in, (olong)seriesNo, (gchar*)label, err);
} // end TimeFilter

extern int TimeFilterIsA (ObitTimeFilter* in) {
  return ObitTimeFilterIsA(in);
} // end TimeFilterIsA 

extern ObitTimeFilter* TimeFilterRef (ObitTimeFilter* in) {
  return ObitTimeFilterRef (in);
} // end TimeFilterRef

extern ObitTimeFilter* TimeFilterUnref (ObitTimeFilter* in) {
  if (!ObitTimeFilterIsA(in)) return NULL;
  return ObitTimeFilterUnref (in);
} // end TimeFilterUnref

extern char* TimeFilterGetName (ObitTimeFilter* in) {
  return in->name;
} // end  TimeFilterGetName

// return dict with {"dTime", "time", "data"}
extern PyObject* TimeFilterGetTime (ObitTimeFilter* in, int seriesNo) {
  PyObject *outDict = PyDict_New();
  PyObject *time, *data;
  olong i;

  PyDict_SetItemString(outDict, "dTime", PyFloat_FromDouble((double)in->dTime));
  time = PyList_New(in->nTime);
  data = PyList_New(in->nTime);
  for (i=0; i<in->nTime++; i++) {
    PyList_SetItem(time, i, PyFloat_FromDouble((double)in->times[i]));
    PyList_SetItem(data, i, PyFloat_FromDouble((double)in->timeData[seriesNo][i]));
  }
  PyDict_SetItemString(outDict, "time", time);
  PyDict_SetItemString(outDict, "data", data);
   
  return outDict;
} // end  TimeFilterGetTime

// return dict with {"dFreq", "freq", "data"}
extern PyObject* TimeFilterGetFreq (ObitTimeFilter* in, int seriesNo) {
  PyObject *outDict = PyDict_New();
  PyObject *freq, *data;
  olong i;

  PyDict_SetItemString(outDict, "dFreq", PyFloat_FromDouble((double)in->dFreq));
  freq = PyList_New(in->nFreq);
  data = PyList_New(in->nFreq);

  for (i=0; i<in->nFreq++; i++) {
    PyList_SetItem(freq, i, PyFloat_FromDouble((double)in->freqs[i]));
    PyList_SetItem(data, i, PyComplex_FromDoubles((double)in->freqData[seriesNo][i*2], 
                                                   (double)in->freqData[seriesNo][i*2+1]));
  }
  PyDict_SetItemString(outDict, "freq", freq);
  PyDict_SetItemString(outDict, "data", data);
   
  return outDict;
} // end  TimeFilterGetFreq

// expect dict with {"dTime", "time", "data"}
extern void TimeFilterSetTime (ObitTimeFilter* in, int seriesNo,
                               PyObject *inDict) {
  PyObject *time, *data;
  olong len, i;

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return;
  }

  in->dTime = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "dTime"));
  time = PyDict_GetItemString(inDict, "time");
  data = PyDict_GetItemString(inDict, "data");
  len = PyList_Size(time);
  if (len>in->nTime) {
    PyErr_SetString(PyExc_TypeError,"Data length incompatible with filter");
    return;
  }

  for (i=0; i<in->nTime++; i++) {
    in->times[i]              = (ofloat)PyFloat_AsDouble(PyList_GetItem(time, i));
    in->timeData[seriesNo][i] = (ofloat)PyFloat_AsDouble(PyList_GetItem(data, i));
  }
} // end  TimeFilterSetTime

// expect dict with {"dFreq", "freq", "data"}
extern void TimeFilterSetFreq (ObitTimeFilter* in, int seriesNo,
                               PyObject *inDict) {
  PyObject *freq, *data, *cx;
  olong len, i;

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return;
  }

  in->dFreq = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "dFreq"));
  freq = PyDict_GetItemString(inDict, "freq");
  data = PyDict_GetItemString(inDict, "data");
  len = PyList_Size(freq);
  if (len>in->nFreq) {
    PyErr_SetString(PyExc_TypeError,"Data length incompatible with filter");
    return;
  }

  for (i=0; i<in->nTime++; i++) {
    in->freqs[i] = (ofloat)PyFloat_AsDouble(PyList_GetItem(freq, i));
    cx = PyList_GetItem(data, i);
    in->freqData[seriesNo][i*2]   = (ofloat)PyComplex_RealAsDouble(cx);
    in->freqData[seriesNo][i*2+1] = (ofloat)PyComplex_ImagAsDouble(cx);
  }
   
} // end  TimeFilterSetFreq

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitTimeFilter *me;
} TimeFilter;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitTimeFilter *me;
} TimeFilter;

%addmethods TimeFilter { 
  TimeFilter(char* name, int nTime, int nSeries) {
     TimeFilter *out;
     out = (TimeFilter*) malloc(sizeof(TimeFilter));
     if (strcmp(name, "None")) out->me = TimeFilterCreate(name, nTime, nSeries);
     else out->me = NULL;
     return out;
   }
  ~TimeFilter() {
    self->me = TimeFilterUnref(self->me);
    free(self);
  }
};

/* $Id: UVDesc.inc,v 1.8 2007/11/17 22:33:07 bcotton Exp $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageDesc type                         */
/*                                                                    */
/*;  Copyright (C)2005,2008                                           */
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
#include "ObitUVDesc.h"
%}

%inline %{
extern ObitUVDesc* UVDescCreate (char *name) {
  return newObitUVDesc (name);
} // end UVDescCreate

extern ObitUVDesc* UVDescCopy (ObitUVDesc* in, 
		              ObitUVDesc* out, ObitErr *err) {
  return ObitUVDescCopy (in, out, err);
} // end UVDescCopy

extern void UVDescCopyDesc (ObitUVDesc* in, ObitUVDesc* out,
			ObitErr *err) {
  ObitUVDescCopyDesc  (in, out, err);
} // end UVDescCopyDesc

extern void UVDescIndex (ObitUVDesc* in) {
  ObitUVDescIndex (in);
} // end UVDescIndex

extern ObitInfoList* UVDescGetList (ObitUVDesc* in) {
  return ObitInfoListRef(in->info);
}
 
extern PyObject *UVDescGetDict(ObitUVDesc* in) {
  PyObject *outDict = PyDict_New();
  PyObject *list1, *list2, *list3, *list4, *list5, *list6, *list7;
  PyObject *value;
  int i, pos = 0;

  PyDict_SetItemString(outDict, "name",    PyString_InternFromString(in->name));
  PyDict_SetItemString(outDict, "object",  PyString_InternFromString(in->object));
  PyDict_SetItemString(outDict, "teles",   PyString_InternFromString(in->teles));
  PyDict_SetItemString(outDict, "instrume",PyString_InternFromString(in->instrument));
  PyDict_SetItemString(outDict, "observer",PyString_InternFromString(in->observer));
  PyDict_SetItemString(outDict, "origin",  PyString_InternFromString(in->origin));
  PyDict_SetItemString(outDict, "date",    PyString_InternFromString(in->date));
  PyDict_SetItemString(outDict, "obsdat",  PyString_InternFromString(in->obsdat));
  PyDict_SetItemString(outDict, "isort",   PyString_InternFromString(in->isort));
  PyDict_SetItemString(outDict, "bunit",   PyString_InternFromString(in->bunit));
  PyDict_SetItemString(outDict, "obsra",   PyFloat_FromDouble((double)in->obsra));
  PyDict_SetItemString(outDict, "obsdec",  PyFloat_FromDouble((double)in->obsdec));
  PyDict_SetItemString(outDict, "epoch",   PyFloat_FromDouble((double)in->epoch));
  PyDict_SetItemString(outDict, "equinox", PyFloat_FromDouble((double)in->equinox));
  PyDict_SetItemString(outDict, "xshift",  PyFloat_FromDouble((double)in->xshift));
  PyDict_SetItemString(outDict, "yshift",  PyFloat_FromDouble((double)in->yshift));
  PyDict_SetItemString(outDict, "altCrpix",PyFloat_FromDouble((double)in->altCrpix));
  PyDict_SetItemString(outDict, "altRef",  PyFloat_FromDouble((double)in->altRef));
  PyDict_SetItemString(outDict, "restFreq",PyFloat_FromDouble((double)in->restFreq));
  PyDict_SetItemString(outDict, "JDObs",   PyFloat_FromDouble((double)in->JDObs));
  PyDict_SetItemString(outDict, "naxis",   PyInt_FromLong((long)in->naxis));
  PyDict_SetItemString(outDict, "nvis",    PyInt_FromLong((long)in->nvis));
  PyDict_SetItemString(outDict, "nrparm",  PyInt_FromLong((long)in->nrparm));
  PyDict_SetItemString(outDict, "ncorr",   PyInt_FromLong((long)in->ncorr));
  PyDict_SetItemString(outDict, "VelReference", PyInt_FromLong((long)in->VelReference));
  PyDict_SetItemString(outDict, "VelDef", PyInt_FromLong((long)in->VelDef));
  list1 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list1, i, PyString_InternFromString(in->ctype[i]));
  PyDict_SetItemString(outDict, "ctype", list1);
  list2 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list2, i, PyFloat_FromDouble((double)in->crval[i]));
  PyDict_SetItemString(outDict, "crval", list2);
  list3 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list3, i, PyInt_FromLong((long)in->inaxes[i]));
  PyDict_SetItemString(outDict, "inaxes", list3);
  list4 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list4, i, PyFloat_FromDouble((double)in->cdelt[i]));
  PyDict_SetItemString(outDict, "cdelt", list4);
  list5 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list5, i, PyFloat_FromDouble((double)in->crpix[i]));
  PyDict_SetItemString(outDict, "crpix", list5);
  list6 = PyList_New(IM_MAXDIM);
  for (i=0; i<IM_MAXDIM; i++) PyList_SetItem(list6, i, PyFloat_FromDouble((double)in->crota[i]));
  PyDict_SetItemString(outDict, "crota", list6);
  list7 = PyList_New(UV_MAX_RANP);
  for (i=0; i<UV_MAX_RANP; i++) PyList_SetItem(list7, i, PyString_InternFromString(in->ptype[i]));
  PyDict_SetItemString(outDict, "ptype", list7);
  PyDict_SetItemString(outDict, "firstVis", PyInt_FromLong((long)in->firstVis));
  PyDict_SetItemString(outDict, "numVisBuff", PyInt_FromLong((long)in->numVisBuff));
/* Structural members - mostly read only */
  PyDict_SetItemString(outDict, "ilocu",    PyInt_FromLong((long)in->ilocu));
  PyDict_SetItemString(outDict, "ilocv",    PyInt_FromLong((long)in->ilocv));
  PyDict_SetItemString(outDict, "ilocw",    PyInt_FromLong((long)in->ilocw));
  PyDict_SetItemString(outDict, "iloct",    PyInt_FromLong((long)in->iloct));
  PyDict_SetItemString(outDict, "ilocb",    PyInt_FromLong((long)in->ilocb));
  PyDict_SetItemString(outDict, "ilocsu",   PyInt_FromLong((long)in->ilocsu));
  PyDict_SetItemString(outDict, "ilocfq",   PyInt_FromLong((long)in->ilocfq));
  PyDict_SetItemString(outDict, "ilocit",   PyInt_FromLong((long)in->ilocit));
  PyDict_SetItemString(outDict, "ilocid",   PyInt_FromLong((long)in->ilocid));
  PyDict_SetItemString(outDict, "ilocws",   PyInt_FromLong((long)in->ilocws));
  PyDict_SetItemString(outDict, "jlocc",    PyInt_FromLong((long)in->jlocc));
  PyDict_SetItemString(outDict, "jlocs",    PyInt_FromLong((long)in->jlocs));
  PyDict_SetItemString(outDict, "jlocf",    PyInt_FromLong((long)in->jlocf));
  PyDict_SetItemString(outDict, "jlocr",    PyInt_FromLong((long)in->jlocr));
  PyDict_SetItemString(outDict, "jlocd",    PyInt_FromLong((long)in->jlocd));
  PyDict_SetItemString(outDict, "jlocif",   PyInt_FromLong((long)in->jlocif));
  PyDict_SetItemString(outDict, "incs",     PyInt_FromLong((long)in->incs));
  PyDict_SetItemString(outDict, "incf",     PyInt_FromLong((long)in->incf));
  PyDict_SetItemString(outDict, "incif",    PyInt_FromLong((long)in->incif));
 
  /* Discard references to newly created objects. */
  while (PyDict_Next(outDict, &pos, NULL, &value))
    Py_DECREF(value);

  return outDict;
} // end UVDescGetDict

extern void UVDescSetDict(ObitUVDesc* in, PyObject *inDict) {
  PyObject *list1, *list2, *list3, *list4, *list5, *list6, *list7;
  char *tstr;
  int i, number;

  if (!PyDict_Check(inDict)) {
	PyErr_SetString(PyExc_TypeError,"Input not a Dict");
        return;
  }

  tstr = PyString_AsString(PyDict_GetItemString(inDict, "object"));
  strncpy (in->object, tstr, UVLEN_VALUE); in->object[UVLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "obsdat"));
  strncpy (in->obsdat, tstr, UVLEN_VALUE); in->obsdat[UVLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "bunit"));
  strncpy (in->bunit, tstr, UVLEN_VALUE); in->bunit[UVLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "teles"));
  strncpy (in->teles, tstr, UVLEN_VALUE); in->teles[UVLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "instrume"));
  strncpy (in->instrument, tstr, IMLEN_VALUE); in->instrument[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "observer"));
  strncpy (in->observer, tstr, IMLEN_VALUE); in->observer[IMLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "origin"));
  strncpy (in->origin, tstr, UVLEN_VALUE); in->origin[UVLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "date"));
  strncpy (in->date, tstr, UVLEN_VALUE); in->date[UVLEN_VALUE-1]=0;
  tstr = PyString_AsString(PyDict_GetItemString(inDict, "isort"));
  strncpy (in->isort, tstr, 2); in->isort[2]=0;
  in->epoch   = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "epoch"));
  in->equinox = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "equinox"));
  in->obsra   = PyFloat_AsDouble(PyDict_GetItemString(inDict, "obsra"));
  in->obsdec  = PyFloat_AsDouble(PyDict_GetItemString(inDict, "obsdec"));
  in->restFreq= (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "restFreq"));
  in->JDObs   = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "JDObs"));
  in->xshift  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "xshift"));
  in->yshift  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "yshift"));
  in->altCrpix= (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "altCrpix"));
  in->altRef  = (float)PyFloat_AsDouble(PyDict_GetItemString(inDict, "altRef"));
  in->naxis   = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "naxis"));
  in->nrparm  = PyInt_AsLong(PyDict_GetItemString(inDict, "nrparm"));
  in->nvis    = PyInt_AsLong(PyDict_GetItemString(inDict, "nvis"));
  in->VelReference  = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "VelReference"));
  in->VelDef  = (int)PyInt_AsLong(PyDict_GetItemString(inDict, "VelDef"));
  list1 = PyDict_GetItemString(inDict, "ctype");
  number = MIN (IM_MAXDIM, PyList_Size(list1));
  for (i=0; i<number; i++) {
    tstr = PyString_AsString(PyList_GetItem(list1, i));
    strncpy (in->ctype[i], tstr, UVLEN_KEYWORD);
  }
  list2 = PyDict_GetItemString(inDict, "crval");
  number = MIN (IM_MAXDIM, PyList_Size(list2));
  for (i=0; i<number; i++) in->crval[i] = PyFloat_AsDouble(PyList_GetItem(list2, i));
  list3 = PyDict_GetItemString(inDict, "inaxes");
  number = MIN (IM_MAXDIM, PyList_Size(list3));
  for (i=0; i<number; i++) in->inaxes[i] = (long)PyInt_AsLong(PyList_GetItem(list3, i));
  list4 = PyDict_GetItemString(inDict, "cdelt");
  number = MIN (IM_MAXDIM, PyList_Size(list4));
  for (i=0; i<number; i++) in->cdelt[i] = (float)PyFloat_AsDouble(PyList_GetItem(list4, i));
  list5 = PyDict_GetItemString(inDict, "crpix");
  number = MIN (IM_MAXDIM, PyList_Size(list5));
  for (i=0; i<number; i++) in->crpix[i] = (float)PyFloat_AsDouble(PyList_GetItem(list5, i));
  list6 = PyDict_GetItemString(inDict, "crota");
  number = MIN (IM_MAXDIM, PyList_Size(list6));
  for (i=0; i<number; i++) in->crota[i] = (float)PyFloat_AsDouble(PyList_GetItem(list6, i));
  list7 = PyDict_GetItemString(inDict, "ptype");
  number = MIN (UV_MAX_RANP, PyList_Size(list7));
  for (i=0; i<number; i++) {
    tstr = PyString_AsString(PyList_GetItem(list7, i));
    strncpy (in->ptype[i], tstr, UVLEN_KEYWORD);
  }
  /* Reindex just to be sure */
  ObitUVDescIndex (in);
} // end UVDescSetDict

ObitUVDesc* UVDescRef (ObitUVDesc* in) {
  return ObitUVDescRef (in);
} // end UVDescRef

ObitUVDesc* UVDescUnref (ObitUVDesc* in) {
  if (!ObitUVDescIsA(in)) return NULL;
  return ObitUVDescUnref (in);
} // end UVDescUnref

extern int UVDescIsA (ObitUVDesc* in) {
  return ObitUVDescIsA(in);
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUVDesc *me;
} UVDesc;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUVDesc *me;
} UVDesc;

%addmethods UVDesc { 
  UVDesc(char *name) {
     UVDesc *out;
     out = (UVDesc *) malloc(sizeof(UVDesc));
     if (strcmp(name, "None")) out->me = UVDescCreate (name);
     else out->me = NULL;
     return out;
   }
  ~UVDesc() {
    self->me = UVDescUnref(self->me);
    free(self);
  }
};

/* $Id: UVGSolve.inc,v 1.2 2006/03/13 16:55:30 bcotton Exp $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV data self calibration utilities     */
/*                                                                    */
/*;  Copyright (C) 2006                                               */
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
#include "ObitUVGSolve.h"
#include "ObitUVSoln.h"
%}


%inline %{
extern ObitUVGSolve* newUVGSolve (char* name) {
  return newObitUVGSolve (name);
} // end  newUVGSolve

extern ObitUVGSolve* UVGSolveCreate (char *name) {
 return ObitUVGSolveCreate(name);
}

extern ObitUVGSolve* UVGSolveCopy  (ObitUVGSolve *in, ObitUVGSolve *out, 
				    ObitErr *err) {
  return ObitUVGSolveCopy (in, out, err);
} // end  UVGSolveCopy

extern ObitUVGSolve* UVGSolveUnref (ObitUVGSolve* in) {
  if (!ObitUVGSolveIsA(in)) return NULL;
  return ObitUVGSolveUnref(in);
}

extern ObitUVGSolve*  UVGSolveRef (ObitUVGSolve* in) {
  return ObitUVGSolveRef(in);
}

extern ObitInfoList* UVGSolveGetList (ObitUVGSolve* in) {
  return ObitInfoListRef(in->info);
}

extern int UVGSolveIsA (ObitUVGSolve* in) {
  return ObitUVGSolveIsA(in);
}

extern ObitTable* UVGSolveCal (ObitUVGSolve *in, ObitUV *inUV, ObitUV *outUV, 
                               ObitUVSel *sel, ObitErr *err)
{
  return (ObitTable*)ObitUVGSolveCal (in, inUV, outUV, sel, err);
} // end UVGSolveCal


%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUVGSolve *me;
} UVGSolve;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUVGSolve *me;
} UVGSolve;

%addmethods UVGSolve { 
  UVGSolve(char* name) {
     UVGSolve *out;
     out = (UVGSolve *) malloc(sizeof(UVGSolve));
     if (strcmp(name, "None")) out->me = newUVGSolve(name);
     else out->me = NULL;
     return out;
   }
  ~UVGSolve() {
   if (self->me->ReferenceCount>0) 
      self->me = UVGSolveUnref(self->me);
   free(self);
  }
};

/* $Id: UVImager.inc,v 1.3 2005/02/06 02:00:39 bcotton Exp $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageMosaic type                       */
/*                                                                    */
/*;  Copyright (C) 2005                                               */
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
#include "ObitUVImager.h"
%}


%inline %{
extern ObitUVImager* newUVImager (char* name) {
  return newObitUVImager (name);
} // end  newUVImager

extern ObitUVImager* UVImagerCopy  (ObitUVImager *in, ObitUVImager *out, 
				    ObitErr *err) {
  return ObitUVImagerCopy (in, out, err);
} // end  UVImagerCopy

extern ObitUVImager* UVImagerUnref (ObitUVImager* in) {
  if (!ObitUVImagerIsA(in)) return NULL;
  return ObitUVImagerUnref(in);
}

extern ObitUVImager*  UVImagerRef (ObitUVImager* in) {
  return ObitUVImagerRef(in);
}

extern ObitUV* UVImagerGetUV (ObitUVImager* in) {
  return ObitUVRef(in->uvdata);
}

extern ObitImageMosaic* UVImagerGetMosaic (ObitUVImager* in, ObitErr *err) {
  return ObitUVImagerGetMosaic(in, err);
}

extern ObitUVImager* UVImagerCreate (char *name, ObitUV *uvData, ObitErr *err) {
 return ObitUVImagerCreate(name, uvData, err);
}

extern void UVImagerWeight (ObitUVImager* in, ObitErr *err) {
 ObitUVImagerWeight(in, err);
}

extern void UVImagerImage (ObitUVImager* in, int field, int doWeight, int doBeam, 
                            int doFlatten, ObitErr *err) {
 gboolean LdoWeight=doWeight, LdoBeam=doBeam, LdoFlatten=doFlatten;
 olong Lfield = field;
 ObitUVImagerImage(in, Lfield, LdoWeight, LdoBeam, LdoFlatten, err);
}

extern void UVImagerFlatten (ObitUVImager* in, ObitErr *err) {
 ObitUVImagerFlatten(in, err);
}

extern char* UVImagerGetName (ObitUVImager* in) {
  if (ObitUVImagerIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int UVImagerIsA (ObitUVImager* in) {
  return ObitUVImagerIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUVImager *me;
} UVImager;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUVImager *me;
} UVImager;

%addmethods UVImager { 
  UVImager(char* name, ObitUV *uvData, ObitErr *err) {
     UVImager *out;
     out = (UVImager *) malloc(sizeof(UVImager));
     if (strcmp(name, "None")) out->me = UVImagerCreate(name, uvData, err);
     else out->me = NULL;
     return out;
   }
  ~UVImager() {
   if (self->me->ReferenceCount>0) 
      self->me = UVImagerUnref(self->me);
   free(self);
  }
};

/* $Id: UV.inc,v 1.23 2008/04/27 20:41:16 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV  type                               */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitUV.h"
#include "ObitIOUVAIPS.h"
#include "ObitIOUVFITS.h"
#include "ObitData.h"
#include "ObitUVUtil.h"
#include "ObitUVEdit.h"
#include "ObitUVWeight.h"
%}


%inline %{

extern void UVSetFITS(ObitUV *in, long nvis, int disk, char *file, 
                      ObitErr *err) {
  ObitUVSetFITS(in, nvis, disk, file, err);
 } // end UVSetFITS

extern void UVSetAIPS(ObitUV *in, long nvis, int disk, int cno, int user, 
                      ObitErr *err) {
  ObitUVSetAIPS(in, nvis, disk, cno, user, err);
 }

extern ObitData* UVCastData (ObitUV* inUV) {
  return (ObitData*)inUV;
} // end  UVCastData

extern ObitUV* UVCreate (char* name) {
  return newObitUV (name);
} // end  UVCreate

extern ObitUV* UVScratch (ObitUV *in, ObitErr *err) {
  return newObitUVScratch (in, err);
} // end  UVScratch

extern PyObject* UVInfo (ObitUV *in, ObitErr *err) {
  ObitIOUVAIPS *AIPSIO=NULL;
  ObitIOUVFITS *FITSIO=NULL;
  PyObject *outDict=PyDict_New();
  PyObject *o=NULL;

  if (err->error) return outDict;

  // Ensure in fully instantiated -assume OK if myIO exists 
  if (!in->myIO) ObitUVFullInstantiate (in, TRUE, err);
  if (err->error) return outDict;

  // Get details and save in dict
  if (ObitIOUVAIPSIsA(in->myIO)) {  // AIPS
    o = PyString_InternFromString("AIPS");
    PyDict_SetItemString(outDict, "type", o);
    AIPSIO = (ObitIOUVAIPS*)in->myIO;
    o = PyInt_FromLong((long)AIPSIO->disk);
    PyDict_SetItemString(outDict, "disk", o);
    o = PyInt_FromLong((long)AIPSIO->CNO);
    PyDict_SetItemString(outDict, "CNO", o);
    o = PyInt_FromLong((long)AIPSIO->UserId);
    PyDict_SetItemString(outDict, "user", o);
  } else if (ObitIOUVFITSIsA(in->myIO)) {  // FITS
    o = PyString_InternFromString("FITS");
    PyDict_SetItemString(outDict, "type", o);
    FITSIO = (ObitIOUVFITS*)in->myIO;
    o = PyInt_FromLong((long)FITSIO->disk);
    PyDict_SetItemString(outDict, "disk", o);
    o = PyString_InternFromString((char*)FITSIO->FileName);
    PyDict_SetItemString(outDict, "filename", o);
  } else {  // Don't know this one
    o = PyString_InternFromString("UNKNOWN");
    PyDict_SetItemString(outDict, "type", o);
  }
  return outDict;
} // end  UVInfo

extern ObitUV* UVZap  (ObitUV *in, ObitErr *err) {
  return ObitUVZap (in, err);
} // end UVZap

extern void UVRename  (ObitUV *in, ObitErr *err) {
  ObitUVRename (in, err);
} // end UVRename

extern ObitUV* UVCopy  (ObitUV *in, ObitUV *out, 
			         ObitErr *err) {
  return ObitUVCopy (in, out, err);
} // end  UVCopy

extern void UVClone (ObitUV *in, ObitUV *out, ObitErr *err) {
   return  ObitUVClone (in, out, err);
} // end  UVClone

// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
// Table verion returned as outValue1
extern ObitTable* newUVTable (ObitUV *in, int access, 
			      char *tabType, long *outValue1, ObitErr *err) {
  ObitIOAccess laccess;
  olong ltabVer = (olong)*outValue1;
  ObitTable *outTable=NULL;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  outTable = newObitUVTable (in, laccess, tabType, &ltabVer, err);
  *outValue1 = (long)ltabVer;
  return outTable;
} // end  newUVTable

extern int UVZapTable (ObitUV *in, char *tabType, long tabVer, 
			    ObitErr *err) {
  ObitIOCode ret;
  ret = ObitUVZapTable (in, tabType, tabVer, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  UVZapTable

extern int UVCopyTables (ObitUV *in, ObitUV *out, char **exclude,
		  	        char **include, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitUVCopyTables  (in, out, exclude, include, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  UVCopyTables

extern int UVUpdateTables (ObitUV *in, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitUVUpdateTables (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  UVUpdateTables

// Open and close to fully instantiate
// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
extern int UVfullInstantiate (ObitUV* in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitUVOpen (in, laccess, err);
  ret = ObitUVClose (in, err);
  if ((err->error) || (ret!=OBIT_IO_OK)) return 1;
  else return 0;
} // end UVfullInstantiate


extern int UVOpen (ObitUV *in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitUVOpen (in, laccess, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Open

extern int UVRead (ObitUV *in, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitUVRead (in, NULL, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Read

extern int UVWrite (ObitUV *in, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitUVWrite (in, NULL, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Write

extern int UVRewrite (ObitUV *in, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitUVRewrite (in, NULL, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Write

extern PyObject *UVGetVisBuf (ObitUV *in) {
  return PyBuffer_FromReadWriteMemory(in->buffer,
				      in->bufferSize * sizeof(ofloat));
}

// force header update 
extern void UVDirty (ObitUV *in) {
  in->myStatus = OBIT_Modified;
} // end Dirty

extern int UVClose (ObitUV *in, ObitErr *err) {
  ObitIOCode ret;
  ret =  ObitUVClose (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Close

extern void UVWeightData (ObitUV *in, ObitErr *err) {
  ObitUVWeightData (in, err);
} // end UVWeightData

extern void UVGetFreq (ObitUV* in, ObitErr *err) {
  ObitUVGetFreq (in, err);
} // end UVGetFreq

extern ObitIOCode UVGetSubA (ObitUV *in, ObitErr *err) {
  ObitIOCode ret;
  ret =  ObitUVGetSubA (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end UVGetSubA

extern ObitUV* UVUnref (ObitUV* in) {
  if (!ObitUVIsA(in)) return NULL;
  return ObitUVUnref(in);
}

extern ObitUV*  UVRef (ObitUV* in) {
  return ObitUVRef(in);
}

extern ObitInfoList* UVGetList (ObitUV* in) {
  return ObitInfoListRef(in->info);
}

extern ObitTableList* UVGetTableList (ObitUV* in) {
  return ObitTableListRef(in->tableList);
}

extern long UVGetHighVer (ObitUV* in, char *tabType) {
  return ObitTableListGetHigh(in->tableList, tabType);
}

extern ObitUVDesc* UVGetDesc (ObitUV* in) {
  return ObitUVDescRef(in->myDesc);
}

extern void UVSetDesc (ObitUV* in, ObitUVDesc* desc) {
  in->myDesc = ObitUVDescUnref(in->myDesc);
  in->myDesc = ObitUVDescRef(desc);
}

extern int UVisScratch (ObitUV* in) {
  return (int)in->isScratch;
}

extern int UVIsA (ObitUV* in) {
  return ObitUVIsA(in);
}

extern char* UVGetName (ObitUV* in) {
  if (ObitUVIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

/* UV Utility functions */

// returns float array val, [0]=maximum baseline length (in U,V), [1] = maximum W
extern void UVUtilUVWExtrema(ObitUV* in, ObitErr *err, float val[2]) {
  ObitUVUtilUVWExtrema (in, &val[0], &val[1], err);
} // end UVUtilUVWExtrema

extern ObitUV* UVUtilCopyZero(ObitUV* in, int scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVUtilCopyZero (in, lscratch, out, err);
} // end UVUtilCopyZero

extern void UVUtilVisDivide(ObitUV* in1, ObitUV *in2, ObitUV *out, ObitErr *err) {
  ObitUVUtilVisDivide (in1, in2, out, err);
} // end UVUtilVisDivide

extern void UVUtilVisSub(ObitUV* in1, ObitUV *in2, ObitUV *out, ObitErr *err) {
  ObitUVUtilVisSub (in1, in2, out, err);
} // end UVUtilVisSub

extern float UVUtilVisCompare(ObitUV* in1, ObitUV *in2, ObitErr *err) {
  return ObitUVUtilVisCompare (in1, in2, err);
} // end UVUtilVisCompare

extern void UVUtilIndex(ObitUV* inUV, ObitErr *err) {
  return ObitUVUtilIndex (inUV, err);
} // end UVUtilVisIndex

extern ObitUV* UVUtilAvgF(ObitUV* in, int scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVUtilAvgF (in, lscratch, out, err);
} // end UVUtilAvgF

extern ObitUV* UVUtilAvgT(ObitUV* in, int scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVUtilAvgT (in, lscratch, out, err);
} // end UVUtilAvgT

extern ObitInfoList* UVUtilCount (ObitUV *inUV, float timeInt, ObitErr *err) {
  return ObitUVUtilCount (inUV, (ofloat)timeInt, err);
} // end UVUtilCount

/* UV Edit Utility functions */

extern void UVEditTD(ObitUV* in, ObitUV *out, ObitErr *err) {
  return ObitUVEditTD (in, out, err);
} // end UVEditTD

extern void UVEditFD(ObitUV* in, ObitUV *out, ObitErr *err) {
  return ObitUVEditFD (in, out, err);
} // end UVEditFD

extern void UVEditStokes(ObitUV* in, ObitUV *out, ObitErr *err) {
  return ObitUVEditStokes (in, out, err);
} // end UVEditStokes

extern ObitUV* UVEditClip(ObitUV* in, int scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVEditClip (in, lscratch, out, err);
} // end UVEditClip

extern ObitUV* UVEditClipStokes(ObitUV* in, int scratch, ObitUV *out, ObitErr *err) {
  gboolean lscratch;
  lscratch = scratch!=0;
  return ObitUVEditClipStokes (in, lscratch, out, err);
} // end UVUtilClipStokes


%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUV *me;
} UV;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUV *me;
} UV;

%addmethods UV { 
  UV(char *name) {
     UV *out;
     out = (UV *) malloc(sizeof(UV));
     if (strcmp(name, "None")) out->me = UVCreate(name);
     else out->me = NULL;
     return out;
   }
  ~UV() { /* Scratch files may be deleted separately*/
   if (self->me->ReferenceCount>0) 
      self->me = ObitUVUnref(self->me);
   free(self);
  }
};

/* $Id: UVSelfCal.inc,v 1.4 2006/03/09 19:01:31 bcotton Exp $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV data self calibration utilities     */
/*                                                                    */
/*;  Copyright (C) 2005-2008                                          */
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
#include "ObitUVSelfCal.h"
#include "ObitUVSoln.h"
%}


%inline %{
extern ObitUVSelfCal* newUVSelfCal (char* name) {
  return newObitUVSelfCal (name);
} // end  newUVSelfCal

extern ObitUVSelfCal* UVSelfCalCreate (char *name, ObitSkyModel *skyModel) {
 return ObitUVSelfCalCreate(name, skyModel);
}

extern ObitUVSelfCal* UVSelfCalCopy  (ObitUVSelfCal *in, ObitUVSelfCal *out, 
				    ObitErr *err) {
  return ObitUVSelfCalCopy (in, out, err);
} // end  UVSelfCalCopy

extern ObitUVSelfCal* UVSelfCalUnref (ObitUVSelfCal* in) {
  if (!ObitUVSelfCalIsA(in)) return NULL;
  return ObitUVSelfCalUnref(in);
}

extern ObitUVSelfCal*  UVSelfCalRef (ObitUVSelfCal* in) {
  return ObitUVSelfCalRef(in);
}

extern ObitInfoList* UVSelfCalGetList (ObitUVSelfCal* in) {
  return ObitInfoListRef(in->info);
}

extern ObitSkyModel* UVSelfCalGetSkyModel (ObitUVSelfCal* in) {
  return ObitSkyModelRef(in->skyModel);
}

extern void UVSelfCalSetSkyModel (ObitUVSelfCal* in, ObitSkyModel *skyModel, 
                                     ObitErr *err) {
  in->skyModel = ObitSkyModelUnref(in->skyModel);  /* Out with the old */
  in->skyModel = ObitSkyModelRef(skyModel);        /* In with the new */
}

extern int UVSelfCalIsA (ObitUVSelfCal* in) {
  return ObitUVSelfCalIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUVSelfCal *me;
} UVSelfCal;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUVSelfCal *me;
} UVSelfCal;

%addmethods UVSelfCal { 
  UVSelfCal(char* name) {
     UVSelfCal *out;
     out = (UVSelfCal *) malloc(sizeof(UVSelfCal));
     if (strcmp(name, "None")) out->me = newUVSelfCal(name);
     else out->me = NULL;
     return out;
   }
  ~UVSelfCal() {
   if (self->me->ReferenceCount>0) 
      self->me = UVSelfCalUnref(self->me);
   free(self);
  }
};

/* $Id: UVSoln2Cal.inc,v 1.2 2007/02/05 15:05:17 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Image utilities                        */
/*                                                                    */
/*;  Copyright (C) 2006-2008                                          */
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
#include "ObitUVSoln2Cal.h"
%}

%inline %{
ObitTable* UVSoln2Cal (ObitUV *inUV, ObitUV *outUV, ObitErr *err) {
  return (ObitTable*)ObitUVSoln2Cal (inUV, outUV, err);
} // end UVSoln2Cal


%}
/* $Id: UVSoln.inc,v 1.2 2007/01/31 15:27:35 bcotton Exp $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV data self calibration utilities     */
/*                                                                    */
/*;  Copyright (C) 2006-2008                                          */
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
#include "ObitUVSoln.h"
#include "ObitTableSNUtil.h"
%}


%inline %{
// Utility routines from ObitUVSoln class
extern int UVSolnRefAnt (ObitTable *SNTab, int isuba, int refant, ObitErr* err)
{
   olong lisuba, lrefant;
   lisuba  = isuba;
   lrefant = isuba;
   ObitUVSolnRefAnt ((ObitTableSN*)SNTab, lisuba, &lrefant, err);
   return lrefant;
} // end UVSolnRefAnt

extern void UVSolnSNSmo (ObitTable *SNTab, int isuba, ObitErr* err)
{
   olong lisuba;
   lisuba = isuba;
   ObitUVSolnSNSmo ((ObitTableSN*)SNTab, lisuba, err);
} // end UVSelfSNSmo

extern void UVSolnDeselSN (ObitTable *SNTab, int isuba, int fqid, int nantf, 
                           int *ants, int nsou, int *sources,
                           float timerange[2], ObitErr* err)
{
  olong lisuba, lfqid, lnantf, lnsou, i;
  olong* lants, *lsources; 
  ofloat ltimerange[2];

  lisuba = isuba;
  lfqid  = fqid;
  lnantf = nantf;
  if (ants[0]==0) lnantf = 0;  // 0 => all
  ltimerange[0] = timerange[0]; ltimerange[1] = timerange[1];
  lants = g_malloc0(lnantf*sizeof(olong));
  for (i=0; i<nantf; i++) lants[i] = ants[i];
  lnsou = (olong)nsou;
  lsources = g_malloc0(lnsou*sizeof(olong));
  for (i=0; i<nsou; i++) lsources[i] = sources[i];
  ObitUVSolnDeselSN ((ObitTableSN*)SNTab, lisuba, lfqid, lnantf, lants, 
                      lnsou, lsources, ltimerange, err);
  g_free(lants);
  g_free(lsources);
} // end UVSolnDeselSN

extern void UVSolnDeselCL (ObitTable *CLTab, int isuba, int fqid, int nantf, 
                           int *ants, int nsou, int *sources,
                           float timerange[2], ObitErr* err)
{
  olong lisuba, lfqid, lnantf, lnsou, i;
  olong* lants, *lsources; 
  ofloat ltimerange[2];

  lisuba = isuba;
  lfqid  = fqid;
  lnantf = nantf;
  if (ants[0]==0) lnantf = 0;  // 0 => all
  ltimerange[0] = timerange[0]; ltimerange[1] = timerange[1];
  lants = g_malloc0(lnantf*sizeof(olong));
  for (i=0; i<nantf; i++) lants[i] = ants[i];
  lnsou = (olong)nsou;
  lsources = g_malloc0(lnsou*sizeof(olong));
  for (i=0; i<nsou; i++) lsources[i] = sources[i];
  ObitUVSolnDeselCL ((ObitTableCL*)CLTab, lisuba, lfqid, lnantf, lants, 
                     lnsou, lsources, ltimerange, err);
  g_free(lants);
  g_free(lsources);
} // end UVSolnDeselCL

extern ObitTable* SNInvert (ObitTable *inSN, ObitData *outData, long tabVer,
                            ObitErr* err)
{
  olong ltabVer = tabVer;
  ObitTable *outSN=NULL;

  outSN = (ObitTable*)ObitTableSNUtilInvert ((ObitTableSN*)inSN, outData, 
    &ltabVer, err);

  return outSN;
} // end SNInvert


%}
/* $Id: UVVis.inc,v 1.1 2007/11/17 22:04:43 bcotton Exp $   */  
/*--------------------------------------------------------------------*/
/* Swig module description for UVVis type                             */
/*                                                                    */
/*;  Copyright (C)2007                                                */
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
#include "ObitUV.h"
#include "ObitUVDesc.h"
%}

%inline %{
extern PyObject* UVVisGet (ObitUV* inUV, ObitErr *err) {
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean doCalSelect;
  ObitIOCode iretCode;
  ofloat cbase, *visp=NULL;
  olong ant1, ant2, suid, fqid, icorr;
  PyObject *vis;
  PyObject *vList, *cvis, *cx, *tup;

   vis = PyDict_New();

   doCalSelect = FALSE;
   ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);

    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);

    if (iretCode==OBIT_IO_EOF)  {
      PyDict_SetItemString(vis, "EOF",PyInt_FromLong((long)1));
      return vis;
    }
    if (iretCode!=OBIT_IO_OK) return vis;

    PyDict_SetItemString(vis, "visNo",  PyInt_FromLong((long)inUV->myDesc->firstVis));
    cbase = inUV->buffer[inUV->myDesc->ilocb]; /* Baseline */
    ant1 = (cbase / 256.0) + 0.001;
    ant2 = (cbase - ant1 * 256) + 0.001;
    PyDict_SetItemString(vis, "ant1",    PyInt_FromLong((long)ant1));
    PyDict_SetItemString(vis, "ant2",    PyInt_FromLong((long)ant2));
    PyDict_SetItemString(vis, "u", PyFloat_FromDouble((double)inUV->buffer[inUV->myDesc->ilocu]));
    PyDict_SetItemString(vis, "v", PyFloat_FromDouble((double)inUV->buffer[inUV->myDesc->ilocv]));
    PyDict_SetItemString(vis, "w", PyFloat_FromDouble((double)inUV->buffer[inUV->myDesc->ilocw]));
    PyDict_SetItemString(vis, "time", PyFloat_FromDouble((double)inUV->buffer[inUV->myDesc->iloct]));
    if (inUV->myDesc->ilocsu>=0) {
        suid = (olong)(inUV->buffer[inUV->myDesc->ilocsu]);
        PyDict_SetItemString(vis, "suid", PyInt_FromLong((long)suid));
    }
    if (inUV->myDesc->ilocfq>=0) {
        fqid = (olong)(inUV->buffer[inUV->myDesc->ilocfq]);
        PyDict_SetItemString(vis, "fqid", PyInt_FromLong((long)fqid));
    }
    vList = PyList_New(inUV->myDesc->ncorr);
    visp = &inUV->buffer[inUV->myDesc->nrparm];
    for (icorr=0; icorr<inUV->myDesc->ncorr; icorr++) {
      tup = PyTuple_New(2);
      cx = PyComplex_FromDoubles((double)*visp, (double)*(visp+1));
      PyTuple_SetItem (tup, 0, cx);
      PyTuple_SetItem (tup, 1, PyFloat_FromDouble((double)*(visp+2)));
      PyList_SetItem(vList, icorr, tup);
      /*Py_DECREF(tup);*/
      visp += 3;
    }

    PyDict_SetItemString(vis, "vis", vList);
    Py_DECREF(vList);

   return vis;
} // end UVVisGet

extern void UVVisSet (PyObject* vis, ObitUV* outUV, ObitErr *err) {
  ObitIOCode oretCode;
  ofloat *visp=NULL;
  olong ant1, ant2, icorr, len;
  PyObject *vList, *cvis, *cx, *tup;


    ant1 = (olong)PyInt_AsLong(PyDict_GetItemString(vis, "ant1"));
    ant2 = (olong)PyInt_AsLong(PyDict_GetItemString(vis, "ant2"));
    outUV->buffer[outUV->myDesc->ilocb] = (ofloat)(ant1*256 + ant2);
    outUV->buffer[outUV->myDesc->ilocu] = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(vis, "u"));
    outUV->buffer[outUV->myDesc->ilocv] = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(vis, "v"));
    outUV->buffer[outUV->myDesc->ilocw] = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(vis, "w"));
    outUV->buffer[outUV->myDesc->iloct] = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(vis, "time"));
    if (outUV->myDesc->ilocsu>=0) {
        outUV->buffer[outUV->myDesc->ilocsu] = (ofloat)PyInt_AsLong(PyDict_GetItemString(vis, "suid"));
    }
    if (outUV->myDesc->ilocfq>=0) {
        outUV->buffer[outUV->myDesc->ilocfq] = (ofloat)PyInt_AsLong(PyDict_GetItemString(vis, "fqid"));
    }

    vList = PyDict_GetItemString(vis, "vis");
    len = PyList_Size(vList);
    if (len!=outUV->myDesc->ncorr) {
      PyErr_SetString(PyExc_TypeError,"UVVis incompatible with ObitUV");
      return;
    }
    visp = &outUV->buffer[outUV->myDesc->nrparm];
    for (icorr=0; icorr<len; icorr++) {
      tup = PyList_GetItem(vList, icorr);
      cx  = PyTuple_GetItem (tup, 0);
      *(visp++) = (ofloat)PyComplex_RealAsDouble(cx);
      *(visp++) = (ofloat)PyComplex_ImagAsDouble(cx);
      *(visp++) = (ofloat)PyFloat_AsDouble(PyTuple_GetItem (tup, 1));
    }

    outUV->myDesc->numVisBuff = 1;
    oretCode = ObitUVWrite (outUV, outUV->buffer, err);

} // end UVVisSet

%}

/* $Id: ZernikeUtil.inc,v 1.1 2006/10/26 21:36:42 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for Zernike polynimoal utilities           */
/*                                                                    */
/*;  Copyright (C) 2006                                               */
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

#include "ObitZernike.h"
%}


%inline %{
/* Return Zernike term n for X and Y */
extern float Zernike (int n, float x, float y) {
  olong ln = (olong)n;
  return (float)ObitZernike (ln, (ofloat)x, (ofloat)y);
} // end Zernike

/* Return Zernike term n gradient in X for X and Y */
extern float ZernikeGradX (int n, float x, float y) {
  olong ln = (olong)n;
  return (float)ObitZernikeGradX (ln, (ofloat)x, (ofloat)y);
} // end ZernikeGradX

/* Return Zernike term n gradient in Y for X and Y */
extern float ZernikeGradY (int n, float x, float y) {
  olong ln = (olong)n;
  return (float)ObitZernikeGradY (ln, (ofloat)x, (ofloat)y);
} // end ZernikeGradY

/* Return Zernike term n for Polar coordinates rho and phi */
extern float ZernikePolar (int n, float rho, float phi) {
  olong ln = (olong)n;
  return (float)ObitZernikePolar (ln, (ofloat)rho, (ofloat)phi);
} // end ZernikePolar


%}


