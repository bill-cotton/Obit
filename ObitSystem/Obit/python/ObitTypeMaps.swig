/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig typemaps for Obit types                                       */
/*                                                                    */
/*   Copyright (C) 2004-2013                                          */
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
%{
// typedef for  Py_ssize_t if python < 2.5
#if PY_MAJOR_VERSION==2 && PY_MINOR_VERSION<5
typedef int Py_ssize_t
#endif
%}

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
           return NULL;
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

