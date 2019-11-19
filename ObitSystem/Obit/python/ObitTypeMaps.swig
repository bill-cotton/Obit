/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig typemaps for Obit types                                       */
/*                                                                    */
/*   Copyright (C) 2004-2019                                          */
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
//%include <typemaps.i>
//%apply long *INPUT { long *longarr };
%{
// typedef for  Py_ssize_t if python < 2.5
#if PY_MAJOR_VERSION==2 && PY_MINOR_VERSION<5
typedef int Py_ssize_t;
#endif
%}

// Convert list into a long array
%typemap (in) long * {
  if (PyList_Check($source)) {
    int size = PyList_Size($source);
    int i = 0;
    $1 = (long*) malloc((size+1)*sizeof(long));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyLong_Check(o)) {
         $1[i] = PyLong_AsLong(o);
      } else if (PyInt_Check(o)) {
         $1[i] = PyInt_AsLong(o);
      } else {
         PyErr_SetString(PyExc_TypeError,"list must contain longs or ints");
         free($1);
         return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// This cleans up the long * array we malloc'd before the function call
%typemap(freearg) long * {
  free((long *) $source);
}

// Convert list into a float array
%typemap (in) float* {
  if (PyList_Check($source)) {
    int size = PyList_Size($source);
    int i = 0;
    $1 = (float*) malloc((size+1)*sizeof(float));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyFloat_Check(o))
         $1[i] = (float)((PyFloatObject*)o)->ob_fval;
      else {
         PyErr_SetString(PyExc_TypeError,"list must contain floats");
         free($1);
         return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}
// This cleans up the float * array we malloc'd before the function call
%typemap(freearg) float * {
  free((float *) $source);
}

// Convert list into a double array
%typemap (in) double * {
  if (PyList_Check($source)) {
    int size = PyList_Size($source);
    int i = 0;
    $1 = (double*) malloc((size+1)*sizeof(double));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyFloat_Check(o))
         $1[i] = (double)((PyFloatObject*)o)->ob_fval;
      else {
         PyErr_SetString(PyExc_TypeError,"list must contain doubles");
         free($1);
         return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
} 
// This cleans up the double * array we malloc'd before the function call
%typemap(freearg) double * {
  free((double *) $source);
}

// This allows a C function to return a long* as a Python list with 7 elements
//%typemap(out) long* {
//  int len,i;
//  len = 7;
//  $result = PyList_New(len);
//  for (i = 0; i < len; i++) {
//    PyList_SetItem($result,i,PyLong_FromLong($1[i]));
//  }
//}

// Clone other simple datatypes from long DAMN, this doesn't work
//%apply long* {int*};
//%apply long* {float*};
//%apply long* {double*};
//%apply long* {gboolean*};

// Convert list of bytes or unicode into a string array
%typemap (in) char ** {
  if (PyList_Check($source)) {
    int size2, size = PyList_Size($source);
    int j, i = 0;
    char *tstr;

    $1 = (char**) malloc((size+1)*sizeof(char*));
    $1[size] = NULL;  // last string NULL
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyUnicode_Check(o)) {
        o = PyUnicode_AsUTF8String(o);
      }
      if (PyBytes_Check(o)) {
         size2 = PyBytes_Size(o);
         $1[i] = (char*) malloc(size2+1);
         tstr = PyBytes_AsString(o);
         for (j=0; j<=size2; j++) $1[i][j] = tstr[j];
      } else {
         PyErr_SetString(PyExc_TypeError,"list must contain Strings");
         free($1);
         return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// This cleans up the char ** array we malloc'd before the function call
%typemap(freearg) char ** {
  int i = 0;
  
  while ($source[i]!=NULL) { // last string should be NULL
    free ((char *) $source[i++]);
  }
  free((char **) $source);
}

%typemap(argout) long *outValue2 {
  PyObject *o;
  if ((!$result) || ($result == Py_None)) {
    $result = PyList_New(0);
    o = PyLong_FromLong($1[0]);
    PyList_Append($result,o);
    Py_XDECREF(o);
    o = PyLong_FromLong($1[1]);
    PyList_Append($result,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($result)) {
      PyObject *o2 = $result;
      $result = PyList_New(0);
      PyList_Append($result,o2);
      Py_XDECREF(o2);
    }
    o = PyLong_FromLong($1[0]);
    PyList_Append($result,o);
    Py_XDECREF(o);
    o = PyLong_FromLong($1[1]);
    PyList_Append($result,o);
    Py_XDECREF(o);
   }
}
// This tells SWIG to treat an int * argument with name 'outLong1' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.

%typemap(argout) long *outLong1 {
  PyObject *o;
  if ((!$result) || ($result == Py_None)) {
    $result = PyList_New(0);
    o = PyLong_FromLong((long)$1[0]);
    PyList_Append($result,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($result)) {
      PyObject *o2 = $result;
      $result = PyList_New(0);
      PyList_Append($result,o2);
      Py_XDECREF(o2);
    }
    o = PyLong_FromLong((long)$1[0]);
    PyList_Append($result,o);
    Py_XDECREF(o);
   }
}
// This tells SWIG to treat a double * argument with name 'outDbl1' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.

%typemap(argout) double *outDbl1 {
  PyObject *o;
  if ((!$result) || ($result == Py_None)) {
    $result = PyList_New(0);
    o = PyFloat_FromDouble($1[0]);
    PyList_Append($result,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($result)) {
      PyObject *o2 = $result;
      $result = PyList_New(0);
      PyList_Append($result,o2);
      Py_XDECREF(o2);
    }
    o = PyFloat_FromDouble($1[0]);
    PyList_Append($result,o);
    Py_XDECREF(o);
   }
}

// This tells SWIG to treat a double * argument with name 'outDbl2' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.

%typemap(argout) double *outDbl2 {
  PyObject *o;
  if ((!$result) || ($result == Py_None)) {
    $result = PyList_New(0);
    o = PyFloat_FromDouble($1[0]);
    PyList_Append($result,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($result)) {
      PyObject *o2 = $result;
      $result = PyList_New(0);
      PyList_Append($result,o2);
      Py_XDECREF(o2);
    }
    o = PyFloat_FromDouble($1[0]);
    PyList_Append($result,o);
    Py_XDECREF(o);
   }
}

// This tells SWIG to treat a float * argument with name 'outFlt1' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.

%typemap(argout) float *outFlt1 {
  PyObject *o;
  if ((!$result) || ($result == Py_None)) {
    $result = PyList_New(0);
    o = PyFloat_FromDouble((double)$1[0]);
    PyList_Append($result,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($result)) {
      PyObject *o2 = $result;
      $result = PyList_New(0);
      PyList_Append($result,o2);
      Py_XDECREF(o2);
    }
    o = PyFloat_FromDouble((double)$1[0]);
    PyList_Append($result,o);
    Py_XDECREF(o);
   }
}


// This tells SWIG to treat a float * argument with name 'outFlt2' as
// an output value.  We'll append the value to the current result which
// is guaranteed to be a List object by SWIG.

%typemap(argout) float *outFlt2 {
  PyObject *o;
  if ((!$result) || ($result == Py_None)) {
    $result = PyList_New(0);
    o = PyFloat_FromDouble((double)$1[0]);
    PyList_Append($result,o);
    Py_XDECREF(o);
   } else {
    if (!PyList_Check($result)) {
      PyObject *o2 = $result;
      $result = PyList_New(0);
      PyList_Append($result,o2);
      Py_XDECREF(o2);
    }
    o = PyFloat_FromDouble((double)$1[0]);
    PyList_Append($result,o);
    Py_XDECREF(o);
   }
}
// Convert list into a ObitImage* array
%typemap (in) ObitImage ** {
  if (PyList_Check($source)) {
    int size = PyList_Size($source);
    int i = 0;

    $1 = (ObitImage**) malloc((size+1)*sizeof(ObitImage*));
    $1[size] = NULL;  // last pointer NULL
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyString_Check(o)) {
        if (!SWIG_IsOK(SWIG_ConvertPtr(o,(void **) &$1[i],$1_descriptor,0))) {
           PyErr_SetString(PyExc_TypeError,"Type error in argument. Expected _ObitImage_p.");
           return NULL;
         }
         if (!ObitImageIsA((ObitImage*)$1[i])) {  // check */
           PyErr_SetString(PyExc_TypeError,"Type error. Expected ObitImage Object.");
           return NULL;
         }
      } else {
         PyErr_SetString(PyExc_TypeError,"list must contain Strings (ObitImage pointers)");
         free($1);
         return NULL;
      }
    }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

