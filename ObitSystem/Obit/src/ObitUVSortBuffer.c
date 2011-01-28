/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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

#include "ObitUVSortBuffer.h"
#include <math.h>
#include <string.h>
#include "glib/gqsort.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVSortBuffer.c
 * ObitUVSortBuffer class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVSortBuffer";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitUVSortBufferClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVSortBufferClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVSortBufferInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVSortBufferClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVSortBufferClassInfoDefFn (gpointer inClass);

/** Private: Sort comparison function for floats */
static gint CompareFloat (gconstpointer in1, gconstpointer in2, 
			  gpointer ncomp);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVSortBuffer* newObitUVSortBuffer (gchar* name)
{
  ObitUVSortBuffer* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVSortBufferClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVSortBuffer));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVSortBufferInit((gpointer)out);

 return out;
} /* end newObitUVSortBuffer */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVSortBufferGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVSortBufferClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVSortBufferGetClass */

/**
 * Make a deep copy of an ObitUVSortBuffer.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitUVSortBuffer* ObitUVSortBufferCopy  (ObitUVSortBuffer *in, 
					 ObitUVSortBuffer *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitUVSortBuffer(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->myUVdata   = ObitUVRef(in->myUVdata);
  out->nvis       = in->nvis;
  out->hiVis      = in->hiVis;
  out->myBuffer   = g_malloc0((in->nvis+2)*in->myUVdata->myDesc->lrec*sizeof(ofloat));
  out->SortStruct = g_malloc0((in->nvis+2)*sizeof(ObitUVSortStruct));
  out->info       = ObitInfoListCopyData(in->info, out->info);
  /* Copy contents of buffer */
  for (i=0; i<in->hiVis*in->myUVdata->myDesc->lrec; i++) 
    out->myBuffer[i] = in->myBuffer[i];

  return out;
} /* end ObitUVSortBufferCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an UVSortBuffer similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitUVSortBufferClone  (ObitUVSortBuffer *in, ObitUVSortBuffer *out, 
			     ObitErr *err)
{
  const ObitClassInfo *ParentClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->myUVdata   = ObitUVRef(in->myUVdata);
  out->nvis       = in->nvis;
  out->hiVis      = 0;
  out->myBuffer   = g_malloc0((in->nvis+2)*in->myUVdata->myDesc->lrec*sizeof(ofloat));
  out->SortStruct = g_malloc0((in->nvis+2)*sizeof(ObitUVSortStruct));
  out->info       = ObitInfoListCopyData(in->info, out->info);

} /* end ObitUVSortBufferClone */

/**
 * Creates an ObitUVSortBuffer 
 * \param name  An optional name for the object.
 * \param inUV  ObitUV from which the buffer is to be created.
 * \param nvis  Size in visibilities of the buffer
 * \param err   Obit error stack object.
 * \return the new object.
 */
ObitUVSortBuffer* ObitUVSortBufferCreate (gchar* name, ObitUV *inUV, 
					  olong nvis, ObitErr *err)
{
  ObitUVSortBuffer* out;

  /* Create basic structure */
  out = newObitUVSortBuffer (name);

  /* Save values */
  out->nvis     = nvis;
  out->myUVdata = ObitUVRef(inUV);

  /* Create buffer */
  out->myBuffer   = g_malloc0((nvis+2)*inUV->myDesc->lrec*sizeof(ofloat));
  /* Sorting structure */
  out->SortStruct = g_malloc0((nvis+2)*sizeof(ObitUVSortStruct));
  out->hiVis = 0;

  return out;
} /* end ObitUVSortBufferCreate */

/**
 * Add a visibility record to the buffer.
 * When the buffer is filled, it is sorted and times up to lastTime
 * written to output.
 * \param in       The object to add to
 * \param vis      Vis record to add (including random parameters)
 * \param lastTime Highest time (days) to write.
 * \param err      Obit error stack object.
 */
void ObitUVSortBufferAddVis (ObitUVSortBuffer *in, ofloat *vis,
			     ofloat lastTime, ObitErr *err)
{
  olong i, lrec, delta, nwrite, nmove, NPIO, bindx, ivis, jvis, ncopy;
  ObitUVSortStruct *sortKeys=NULL;
  gchar *routine = "ObitUVSortBufferAddVis";

  /* error checks */
  if (err->error) return;

  /* Sanity check */
  Obit_return_if_fail((in->hiVis<in->nvis), err, 
 		      "%s: Sort buffer full", routine);

  /* Add to buffer */
  lrec  = in->myUVdata->myDesc->lrec;
  ncopy = lrec*sizeof(ofloat);
  bindx = in->hiVis*lrec;
  memmove(&in->myBuffer[bindx], vis, ncopy);
  in->hiVis++;

  /* Is it full? */
  if (in->hiVis<in->nvis) return;   /* Nope */

  /* Sort */
  ObitUVSortBufferSort(in, err);

  /* How big is I/O buffer? */
  NPIO = in->myUVdata->bufferSize/lrec;
  
  /* Write with times up to lastTime */
  delta = sizeof(ObitUVSortStruct)/sizeof(ofloat);
  ivis = 0;
  while (ivis<in->hiVis) {  /* outer loop */
    nwrite = 0;
    bindx  = 0;
    nmove = MIN(NPIO-2, in->hiVis-ivis);
    
    /* Inner loop copying to I/O buffer */
    for (i=0; i<nmove; i++) {
      sortKeys = (ObitUVSortStruct*)&in->SortStruct[ivis*delta];
      if (sortKeys->key[0]>lastTime) break; /* Past lastTime? */
      /* Copy to I/O buffer */
      jvis = sortKeys->index.itg;
      memmove(&in->myUVdata->buffer[bindx], &in->myBuffer[jvis*lrec], ncopy);
      bindx += lrec;
      nwrite++;
      ivis++;
    } /* End loop copying to IO buffer */
  
    /* Write */
    in->myUVdata->myDesc->numVisBuff = nwrite;
    ObitUVWrite (in->myUVdata, in->myUVdata->buffer, err);
    if(err->error) Obit_traceback_msg (err, routine, in->name);
    
    if (sortKeys->key[0]>lastTime) break; /* Past lastTime? */
  } /* end outer loop */

  /* Copy rest to start of Sort Buffer */
  nmove = 0;
  bindx  = 0;
  while (ivis<in->hiVis) {  /* outer loop */
    sortKeys = (ObitUVSortStruct*)&in->SortStruct[ivis*delta];
    jvis = sortKeys->index.itg;
    memmove(&in->myBuffer[bindx], &in->myBuffer[jvis*lrec], ncopy);
    bindx += lrec;
    nmove ++;
    ivis++;
  }
  in->hiVis = nmove;   /* New number in buffer */
} /* end ObitUVSortBufferAddVis */

/**
 * Sort and write any remaining data in buffer written to output.
 * \param in    The object to flush
 * \param err   Obit error stack object.
 */
void ObitUVSortBufferFlush (ObitUVSortBuffer *in, ObitErr *err)
{
  olong i, ivis, jvis, lrec, delta, nwrite, nmove, NPIO, bindx, ncopy;
  ObitUVSortStruct *sortKeys;
  gchar *routine = "ObitUVSortBufferFlush";

  /* error checks */
  if (err->error) return;

  /* Sort */
  ObitUVSortBufferSort(in, err);

  /* How big is I/O buffer? */
  lrec = in->myUVdata->myDesc->lrec;
  NPIO = in->myUVdata->bufferSize/lrec;
  ncopy = lrec*sizeof(ofloat);
  
  /* Loop writing to I/O buffer and writing */
  delta = sizeof(ObitUVSortStruct)/sizeof(ofloat);
  ivis = 0;
  while (ivis<in->hiVis) {  /* outer loop */
    nwrite = 0;
    bindx  = 0;
    nmove = MIN(NPIO-2, in->hiVis-ivis);
    
    /* Inner loop copying to I/O buffer */
    for (i=0; i<nmove; i++) {
      sortKeys = (ObitUVSortStruct*)&in->SortStruct[ivis*delta];
      /* Copy to I/O buffer */
      jvis = sortKeys->index.itg;
      memmove(&in->myUVdata->buffer[bindx], &in->myBuffer[jvis*lrec], ncopy);
      bindx += lrec;
      nwrite ++;
      ivis++;

    } /* End loop copying to IO buffer */
  
    /* Write */
    in->myUVdata->myDesc->numVisBuff = nwrite;
    ObitUVWrite (in->myUVdata, in->myUVdata->buffer, err);
    if(err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end outer loop */

  in->hiVis = 0;

} /* end ObitUVSortBufferFlush */

/**
 * Sort indices of contents of buffer into time order
 * \param in  The object to sort
 * \param err Obit error stack object.
 */
void ObitUVSortBufferSort (ObitUVSortBuffer *in, ObitErr *err)
{
  olong i, delta, lrec, iloct, ilocb, number, size, ncomp; 
  ObitUVSortStruct *sortKeys;

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Set sort keys */
  delta = sizeof(ObitUVSortStruct)/sizeof(ofloat);
  lrec  = in->myUVdata->myDesc->lrec;
  iloct = in->myUVdata->myDesc->iloct;
  ilocb = in->myUVdata->myDesc->ilocb;
  for (i=0; i<in->hiVis; i++) {
    sortKeys = (ObitUVSortStruct*)&in->SortStruct[i*delta];
    sortKeys->index.itg = i;
    sortKeys->key[0]    = in->myBuffer[i*lrec+iloct];
    sortKeys->key[1]    = in->myBuffer[i*lrec+ilocb];
  }

  /* Sort SortStruct */
  number = in->hiVis;
  size   = sizeof(ObitUVSortStruct);
  ncomp  = 2;
  g_qsort_with_data (in->SortStruct, number, size, CompareFloat, &ncomp);

} /* end ObitUVSortBufferSort */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVSortBufferClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVSortBufferClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVSortBufferClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVSortBufferClassInfoDefFn (gpointer inClass)
{
  ObitUVSortBufferClassInfo *theClass = (ObitUVSortBufferClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVSortBufferClassInit;
  theClass->newObit       = (newObitFP)newObitUVSortBuffer;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVSortBufferClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVSortBufferGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitUVSortBufferCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVSortBufferClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVSortBufferInit;
  theClass->ObitUVSortBufferCreate = (ObitUVSortBufferCreateFP)ObitUVSortBufferCreate;
  theClass->ObitUVSortBufferAddVis = (ObitUVSortBufferAddVisFP)ObitUVSortBufferAddVis;
  theClass->ObitUVSortBufferSort   = (ObitUVSortBufferSortFP)ObitUVSortBufferSort;
  theClass->ObitUVSortBufferFlush  = (ObitUVSortBufferFlushFP)ObitUVSortBufferFlush;


} /* end ObitUVSortBufferClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVSortBufferInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVSortBuffer *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread     = newObitThread();
  in->info       = newObitInfoList();
  in->myUVdata   = NULL;
  in->myBuffer   = NULL;
  in->SortStruct = NULL;

} /* end ObitUVSortBufferInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVSortBuffer* cast to an Obit*.
 */
void ObitUVSortBufferClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVSortBuffer *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->info)       ObitInfoListUnref (in->info); in->info = NULL;
  if (in->thread)     in->thread   = ObitThreadUnref (in->thread);
  if (in->myUVdata)   in->myUVdata = ObitUVUnref(in->myUVdata);
  if (in->myBuffer)   g_free(in->myBuffer);    in->myBuffer   = NULL;
  if (in->SortStruct) g_free(in->SortStruct);  in->SortStruct = NULL;
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVSortBufferClear */

/**
 * Compare two lists of floats
 * Conformant to function type GCompareDataFunc
 * \param in1   First list, preceeded by olong index
 * \param in2   Second list, preceeded by olong index
 * \param ncomp Number of values to compare
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CompareFloat (gconstpointer in1, gconstpointer in2, 
			  gpointer ncomp)
{
  gint out = 0;
  olong nc, i;
  ObitUVSortStruct *float1, *float2;

  /* get correctly typed local values */
  float1 = (ObitUVSortStruct*)in1;
  float2 = (ObitUVSortStruct*)in2;
  nc = *(olong*)ncomp;

  /* List or single value? */
  if (nc==1) {
    if (float1->key[0]<float2->key[0])      out = -1;
    else if (float1->key[0]>float2->key[0]) out = 1;
    else                                    out = 0;
  } else { /* list */
    for (i=0; i<nc; i++) {
      if (float1->key[i]<float2->key[i])      out = -1;
      else if (float1->key[i]>float2->key[i]) out = 1;
      else                                    out = 0;
      if (out) break;   /* stop at first not equal */
    }
  }

  return out;
} /* end CompareFloat */

