/* $Id: ObitFitRegionList.c,v 1.1 2006/05/12 14:05:55 bcotton Exp $        */
/*--------------------------------------------------------------------*/
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

#include "ObitFitRegionList.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFitRegionList.c
 * ObitFitRegionList class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitFitRegionList";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitFitRegionListClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitFitRegionListClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitFitRegionListInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitFitRegionListClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitFitRegionListClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitFitRegionList* newObitFitRegionList (gchar* name)
{
  ObitFitRegionList* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFitRegionListClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitFitRegionList));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitFitRegionListInit((gpointer)out);

 return out;
} /* end newObitFitRegionList */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitFitRegionListGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFitRegionListClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitFitRegionListGetClass */

/**
 * Creates an ObitFitRegionList 
 * \param name  An optional name for the object.
 * \param image ObitImage being described
 * \return the new object.
 */
ObitFitRegionList* ObitFitRegionListCreate (gchar* name, ObitImage *image)
{
  ObitFitRegionList* out;

  /* Create basic structure */
  out = newObitFitRegionList (name);
  out->image = ObitImageRef(image);

  return out;
} /* end ObitFitRegionListCreate */

/**
 * Subtract ObitFitRegionList from its image
 * \param in       List of region models in image.
 * \param outImage Image to write results into
 * \param err      Obit Error/message stack
 */
void ObitFitRegionListSubtract (ObitFitRegionList *in, ObitImage *outImage,
				ObitErr *err)
{
  g_error("WRITE ME");
} /* end ObitFitRegionListSubtract */

/**
 * Append reg to list in
 * \param in   Object with table to add reg to
 * \param reg the element to add. MUST have unique name in list
 */
void ObitFitRegionListAppend (ObitFitRegionList *in, ObitFitRegion *reg)
{
  ObitFitRegion *tmp = NULL;

  /* Make sure it's not already in list */
  tmp =  ObitFitRegionListFind (in, reg->name);
  if (tmp!=NULL) { /* trouble - die in a noisy fashion */
    g_error ("ObitFitRegionAppend: trying to add redundant island: %s ",
	     reg->name);
  }

  /* append to list */
  in->list = g_slist_append  (in->list, reg);
  in->number++;

} /* end ObitFitRegionAppend */

/**
 * Remove region from list in
 * \param in   Object with table to remove elem from
 * \param reg  the element to remove. MUST be in list, unreffed
 */
void ObitFitRegionListRemove (ObitFitRegionList *in, ObitFitRegion *reg)
{
  /* remove from table */
  in->list = g_slist_remove (in->list, reg);

  in->number--; /* keep count */

  ObitFitRegionUnref (reg);

} /* end ObitFitRegionRemove  */

/**
 * Find a region in list in
 * \param in    Object with table to search
 * \param name  Name of region
 * \return pointer to Region containing item, NULL if not found.
 */
ObitFitRegion* ObitFitRegionListFind (ObitFitRegionList *in, gchar *name)
{
  GSList *tmp;
  ObitFitRegion *reg;

  /* error checks */
  g_assert (ObitFitRegionListIsA(in));
  g_assert (name != NULL);

  /* loop through list testing elements */
  tmp = in->list;
  while (tmp!=NULL) {
    reg = (ObitFitRegion*)tmp->data;
    /* check if this is a match */
    if (!strcmp(reg->name, name)) return reg;
    tmp = g_slist_next(tmp);
  }
  return NULL; /* didn't find */
} /* end ObitFitRegionFind */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitFitRegionListClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitFitRegionListClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitFitRegionListClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitFitRegionListClassInfoDefFn (gpointer inClass)
{
  ObitFitRegionListClassInfo *theClass = (ObitFitRegionListClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitFitRegionListClassInit;
  theClass->newObit       = (newObitFP)newObitFitRegionList;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitFitRegionListClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitFitRegionListGetClass;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitFitRegionListClear;
  theClass->ObitInit      = (ObitInitFP)ObitFitRegionListInit;
  theClass->ObitFitRegionListCreate = 
    (ObitFitRegionListCreateFP)ObitFitRegionListCreate;
} /* end ObitFitRegionListClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitFitRegionListInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFitRegionList *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->image  = NULL;
  in->list   = NULL;
  in->number = 0;

} /* end ObitFitRegionListInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitFitRegionList* cast to an Obit*.
 */
void ObitFitRegionListClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFitRegionList *in = inn;
  GSList *tmp;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->image  = ObitImageUnref(in->image);
  if (in->list) {
    /* loop through list deleting elements */
    tmp = in->list;
    while (tmp!=NULL) {
      if (tmp->data) ObitFitRegionUnref(tmp->data);
      tmp = g_slist_next(tmp);
    }
    /* delete members  */
    g_slist_free(in->list);
   }
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitFitRegionListClear */

