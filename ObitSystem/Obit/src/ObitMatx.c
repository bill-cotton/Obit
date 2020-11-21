/* $Id$         */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2020                                               */
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

#include "ObitMatx.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitMatx.c
 * ObitMatx class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitMatx";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitMatxClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitMatxClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitMatxInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitMatxClear (gpointer in);

/** Private: Inatilize Class  */
void ObitMatxClassInit (void);

/** Private: Set Class function pointers. */
static void ObitMatxClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitMatx* newObitMatx (gchar* name)
{
  ObitMatx* out;

  /* Class initialization if needed 
  if (!myClassInfo.initialized) ObitMatxClassInit();*/

  /* allocate/init structure */
  out = ObitMemAlloc0Name(sizeof(ObitMatx), "ObitMatx");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitMatxInit((gpointer)out);

 return out;
} /* end newObitMatx */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitMatxGetClass (void)
{
  /* Class initialization if needed 
  if (!myClassInfo.initialized) ObitMatxClassInit();*/

  return (gconstpointer)&myClassInfo;
} /* end ObitMatxGetClass */

/**
 * Make a deep copy of an ObitMatx.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitMatx* ObitMatxCopy  (ObitMatx *in, ObitMatx *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, j, n, m;
  ofloat val[3];

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  //g_assert (ObitIsA(in, &myClassInfo));
  //if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitMatx(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* arrays */
  n = in->naxis[0]; m = in->naxis[1];  /* Only 2D */
  for (j=0; j<n; j++) {
    for (i=0; j<m; j++) {
      ObitMatxGet(in, i, j, val);
      ObitMatxSet(out, val, i, j);
    }
  }
  return out;
} /* end ObitMatxCopy */

/**
 * Determine if the two input objects have a compatable geometry.
 * Must be same type and geometry
 * \param in1  First object to test.
 * \param in2  Second object to test.
 * \return TRUE if compatable, else FALSE.
 */
gboolean ObitMatxIsCompatible (ObitMatx *in1, ObitMatx *in2)
{
  olong i;

  /* test fails if either input is NULL */
  if (in1==NULL) return FALSE;
  if (in2==NULL) return FALSE;
  /* Same type? */
  if (in1->type!=in2->type) return FALSE;

   /* Same number of non degenerate dimensions? */
  if (in1->ndim!=in2->ndim) {
    /* Don't bother if extra dimensions are length 1 */
    if (in1->ndim>in2->ndim) {  /* in1 has more */
      for (i=in2->ndim; i<in1->ndim; i++)
	if (in1->naxis[i]>1) return FALSE;
    } else {  /* in2 has more */
      for (i=in1->ndim; i<in2->ndim; i++)
	if (in2->naxis[i]>1) return FALSE;
    }
  }  /* end dimension check */

  return TRUE; /* must be OK */
} /* end ObitMatxIsCompatible  */

/**
 * Creates an ObitMatx of a specified type and geometry.
 * \param type  type, OBIT_Real, OBIT_Complex, OBIT_Double
 * \param ndim  Number of dimensions desired, ONLY 2 currently
 * \param naxis Dimensionality along each axis. NULL => don't create array.
 * \return the new object.
 */
ObitMatx* ObitMatxCreate (ObitMatxType type, olong ndim, olong *naxis)
{
  ObitMatx* out;
  olong i, size;

  /* Initialize class if needed */
  if (!myClassInfo.initialized) ObitMatxClassInit();

  /* Create basic structure */
  out = newObitMatx ("Matx");

  /* create data array if wanted */
  if ((ndim<0) || (naxis==NULL)) return out;
  g_assert (ndim==2); /* sanity check */

  /* copy geometry */
  out->ndim = ndim;
  out->naxis = ObitMemAlloc0Name(ndim*sizeof(olong), "Matx naxis");
  size = 1;
  for (i=0; i<ndim; i++) {
    out->naxis[i] = MAX (1, naxis[i]);
    size *= out->naxis[i]; /* total size */
  }

  out->type = type;  /* data type */

  switch (type) {
  case OBIT_Real:
    out->array = ObitMemAlloc0Name(size*sizeof(ofloat) + 
					out->naxis[0]*sizeof(ofloat),
					"Matx array");
    out->flt = (ofloat*)out->array;
   break;
  case OBIT_Complex:
    size *= 2;
    out->array = ObitMemAlloc0Name(size*sizeof(ofloat) + 
					out->naxis[0]*sizeof(ofloat),
					"Matx array");
    out->cpx = (ocomplex*)out->array;
    break;
  case OBIT_Double:
    size *= 2;
    out->array = ObitMemAlloc0Name(size*sizeof(ofloat) + 
					out->naxis[0]*sizeof(ofloat),
					"Matx array");
    out->dbl = (odouble*)out->array;
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };
  out->arraySize = size;

  return out;
} /* end ObitMatxCreate */

/**
 *  Matrix inner multiply
 *  out = in1 * in2
 * in1, in2 and out are square matrices of the same dimension
 * \param in1  Input object
 * \param in2  Input object
 * \param out  Output object
 */
void ObitMatxMult (ObitMatx* in1, ObitMatx* in2, ObitMatx* out)
{
  olong ir, ic, ii, n;
  olong indx1, indx2, indxo;
  ofloat sumr;
  odouble sumd;
  ocomplex tcpx;
  
  /* error checks */
  g_assert (ObitMatxIsCompatible(in1, in2));
  g_assert (ObitMatxIsCompatible(in1, out));
  g_assert (in1->naxis[0]==in1->naxis[1]);  /* Square */

  n = in1->naxis[0];
  switch (in1->type) {
  case OBIT_Real:
    for (ir=0; ir<n; ir++) {
      for (ic=0; ic<n; ic++) {
	sumr = 0.0;
	indxo = ic + ir*n;
	for (ii=0; ii<n; ii++) {
	  indx1 = (ii*n + ic);
	  indx2 = (ii + ir*n);
	  sumr += in1->flt[indx1]*in2->flt[indx2];
	}
	out->flt[indxo]   = sumr;
      }
    }
    break;
  case OBIT_Complex:
    for (ir=0; ir<n; ir++) {
      for (ic=0; ic<n; ic++) {
	indxo = (ic + ir*n);
	COMPLEX_SET (out->cpx[indxo], 0.0, 0.0);
	for (ii=0; ii<n; ii++) {
	  indx1 = (ii*n + ic);
	  indx2 = (ii + ir*n);
	  COMPLEX_MUL2(tcpx, in1->cpx[indx1], in2->cpx[indx2]);
	  COMPLEX_ADD2 (out->cpx[indxo], out->cpx[indxo], tcpx);
	}
      }
    }
    break;
  case OBIT_Double:
    for (ir=0; ir<n; ir++) {
      for (ic=0; ic<n; ic++) {
	sumd = 0.0;
	indxo = ic + ir*n;
	for (ii=0; ii<n; ii++) {
	  indx1 = (ii*n + ic);
	  indx2 = (ii + ir*n);
	  sumd += in1->dbl[indx1]*in2->dbl[indx2];
	}
	out->dbl[indxo]   = sumd;
      }
    }
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };

 /* end row loop */
  
} /* end ObitMatxMatrixMult */

/**
 *  Add corresponding elements of the array.
 *  out = in1 + in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitMatxAdd (ObitMatx* in1, ObitMatx* in2, ObitMatx* out)
{
  olong i;

   /* error checks */
  g_assert (ObitMatxIsCompatible(in1, in2));
  g_assert (ObitMatxIsCompatible(in1, out));

  switch (in1->type) {
  case OBIT_Real:
   for (i=0; i<in1->arraySize; i++)
     {out->flt[i] = in1->flt[i] + in2->flt[i];}
   break;
  case OBIT_Complex:
   for (i=0; i<in1->arraySize/2; i++)
     {COMPLEX_ADD2(out->cpx[i], in1->cpx[i], in2->cpx[i]);}
    break;
  case OBIT_Double:
   for (i=0; i<in1->arraySize/2; i++)
     {out->dbl[i] = in1->dbl[i] + in2->dbl[i];}
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };

} /* end ObitMatxAdd */


/**
 *  Subtract corresponding elements of the arrays.
 *  out = in1 - in2
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output object
 */
void ObitMatxSub (ObitMatx* in1, ObitMatx* in2, ObitMatx* out)
{
  olong i;

   /* error checks */
  g_assert (ObitMatxIsCompatible(in1, in2));
  g_assert (ObitMatxIsCompatible(in1, out));

  switch (in1->type) {
  case OBIT_Real:
   for (i=0; i<in1->arraySize; i++)
     {out->flt[i] = in1->flt[i] - in2->flt[i];}
   break;
  case OBIT_Complex:
   for (i=0; i<in1->arraySize/2; i++)
     {COMPLEX_SUB(out->cpx[i], in1->cpx[i], in2->cpx[i]);}
    break;
  case OBIT_Double:
   for (i=0; i<in1->arraySize/2; i++)
     {out->dbl[i] = in1->dbl[i] - in2->dbl[i];}
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };
 } /* end ObitMatxSub */


/**
 *  Transpose and conjugate a matrix.
 * \param in      Input object
 * \param out     Output object
 */
void ObitMatxCTrans (ObitMatx* in, ObitMatx* out)
{
  olong i, j, n, m;
  ofloat val[3];
  gboolean conjg = in->type==OBIT_Complex;
  /* error checks */
  g_assert (ObitMatxIsCompatible(in, out));
  n = in->naxis[0]; m = in->naxis[1];  /* Only 2D */
  for (j=0; j<n; j++) {
    for (i=0; i<m; i++) {
      ObitMatxGet(in, i, j, val);
      if (conjg) val[1] = -val[1];
      ObitMatxSet(out, val, j, i);
    }
  }
} /* end ObitMatxCTrans */


/**
 * Zero the elements of an array 
 * \param in      Input object to zero
 */
void ObitMatxZero (ObitMatx* in)
{
  olong i;
  ofloat val[3] = {0.0,0.0,0.0};

  /* error checks */
  //g_assert (ObitIsA(in, &myClassInfo));
  //g_assert (in->array != NULL);

  switch (in->type) {
  case OBIT_Real:
    for (i=0; i<in->arraySize; i++) in->flt[i] = 0.0;
   break;
  case OBIT_Complex:
    for (i=0; i<in->arraySize/2; i++) 
      COMPLEX_SET(in->cpx[i], val[0], val[1]); 
    break;
  case OBIT_Double:
    for (i=0; i<in->arraySize/2; i++) in->dbl[i] = 0.0;
    break;
  default:
     g_assert_not_reached(); /* unknown, barf */
   };
} /* end ObitMatxZero */

/**
 * Fills values in a 2x2 complex ObitMatx 
 * \param matx  ObitMatx to fill
 * \param R00   Real part of element (0,0)
 * \param I00   Imaginary part of element (0,0)
 * \param R01   Real part of element (0,1)
 * \param I01   Imaginary part of element (0,1)
 * \param R10   Real part of element (1,0)
 * \param I10   Imaginary part of element (1,0)
 * \param R11   Real part of element (1,1)
 * \param I11   Imaginary part of element (1,1)
 */
void ObitMatxSet2C(ObitMatx *matx, ofloat R00, ofloat I00, ofloat R01, ofloat I01, 
		   ofloat R10, ofloat I10, ofloat R11, ofloat I11)
{
  ofloat val[2];
  val[0] = R00; val[1]= I00; ObitMatxSet(matx, val, 0, 0);
  val[0] = R01; val[1]= I01; ObitMatxSet(matx, val, 0, 1);
  val[0] = R10; val[1]= I10; ObitMatxSet(matx, val, 1, 0);
  val[0] = R11; val[1]= I11; ObitMatxSet(matx, val, 1, 1);
} /* end ObitMatxSet2 */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitMatxClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitMatxClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitMatxClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitMatxClassInfoDefFn (gpointer inClass)
{
  ObitMatxClassInfo *theClass = (ObitMatxClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitMatxClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitMatxClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitMatxGetClass;
  theClass->newObit       = (newObitFP)newObitMatx;
  theClass->ObitCopy      = (ObitCopyFP)ObitMatxCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitMatxClear;
  theClass->ObitInit      = (ObitInitFP)ObitMatxInit;
  theClass->ObitMatxCreate = (ObitMatxCreateFP)ObitMatxCreate;
  theClass->ObitMatxIsCompatible = 
    (ObitMatxIsCompatibleFP)ObitMatxIsCompatible;
  theClass->ObitMatxMult   = (ObitMatxMultFP)ObitMatxMult;
  theClass->ObitMatxAdd    = (ObitMatxAddFP)ObitMatxAdd;
  theClass->ObitMatxSub    = (ObitMatxSubFP)ObitMatxSub;
  theClass->ObitMatxCTrans = (ObitMatxCTransFP)ObitMatxCTrans;
  theClass->ObitMatxZero   = (ObitMatxZeroFP)ObitMatxZero;
  theClass->ObitMatxSet2C  = (ObitMatxSet2CFP)ObitMatxSet2C;
} /* end ObitMatxClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitMatxInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitMatx *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->array        = NULL;
  in->arraySize    = 0;
  in->ndim         = 0;
  in->naxis        = NULL;

} /* end ObitMatxInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitMatx* cast to an Obit*.
 */
void ObitMatxClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitMatx *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->array)  in->array = ObitMemFree(in->array);
  if (in->naxis)  in->naxis = ObitMemFree(in->naxis);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitMatxClear */

