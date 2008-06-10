/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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

#include "ObitTableList.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitTableList.c
 * ObitTableList class function definitions.
 *
 * This is a list of associated tables.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitTableList";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/*--------------- File Structure definitions  ----------------*/
/** An element of an ObitTableList */
typedef struct { 
  /**  Table name */
  gchar *name;
  /** Table version number */
  olong version;
  /** ObitTable */
  ObitTable *table;
}  ObitTableListElem; 


/**
 * ClassInfo global structure ObitTableListClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitTableListClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitTableListInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitTableListClear (gpointer in);

/** Private: Find item in a list */
static ObitTableListElem* 
ObitTableListFind(ObitTableList *in, gchar *name, olong *version);

/** Private: Create an ObitTableListElem. */
static ObitTableListElem*
newObitTableListElem (gchar *name, olong version, ObitTable *table);

/** Private: Copy an ObitTableListElem. */
static ObitTableListElem* ObitTableListElemCopy (ObitTableListElem* in);

/** Private: Delete an ObitTableListElem. */
static void freeObitTableListElem (ObitTableListElem* in);

/** Private: Set Class function pointers. */
static void ObitTableListClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \return the new object.
 */
ObitTableList* newObitTableList (gchar* name)
{
  ObitTableList* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableListClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitTableList));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitTableListInit((gpointer)out);

  return out;
} /* end newObitTableList */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitTableListGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableListClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetIOClass */

/**
 * Make a copy of a object.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitTableList* ObitTableListCopy  (ObitTableList *in, ObitTableList *out, 
				   ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  GSList  *tmp;
  ObitTableListElem *elem, *telem;

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
    out = newObitTableList(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  /* number of elements */
  out->number = in->number;
  /* loop through list copying elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (ObitTableListElem*)tmp->data;
    /* make copy */
    telem = ObitTableListElemCopy (elem);
    out->list = g_slist_prepend(out->list, telem); /* add to new list */
    tmp = g_slist_next(tmp);
  }

  /* reverse to get into same order */
  out->list = g_slist_reverse(out->list);

  return out;
} /* end ObitTableListCopy */

/**
 * Store information in the TableList.
 * If the requested name and version does not exist then an entry is created.
 * If a previous entry exists, the information is updated.
 * \param in      Pointer to TableList.
 * \param name    The name of the table to enter.
 * \param version Version number of table, >= 0 -> highest value found +1 and this
 *                value is returned.
 * \param table   ObitTable to reference, may be NULL if table but no object exists.
 *                If NULL and an existing entry has a NonNULL value, the old value is kept.
 */
void 
ObitTableListPut(ObitTableList *in, 
		 gchar* name, olong *version, ObitTable *table, ObitErr *err)
{
  ObitTableListElem *elem;
  gboolean newEntry;
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return; /* error condition exists? */
  g_assert (ObitTableListIsA(in));
  if (table!=NULL) g_assert (ObitTableIsA(table));
  g_assert (name != NULL);
  g_assert (version != NULL);
  
  /* look up */
  newEntry = (*version <= 0); /* force new entry? */
  elem = ObitTableListFind(in, name, version);

  /* New entry forced? */
  if (newEntry) {
    *version = MAX (1, (*version)+1);
    elem = NULL;
  }

  if (elem!=NULL) { /* Found it */
    /* Relink table */
    if (table!=NULL) {
      elem->table = ObitTableUnref(elem->table);
      elem->table = ObitTableRef(table);
    }
    return; /* done */
  } /* end of mismatch */
    
  /* add new one */
  elem =  newObitTableListElem(name, *version, table);
  in->list = g_slist_append (in->list, elem); /* add to list */
  in->number++; /* keep count */

} /* end ObitTableListPut */


/**
 * Retrieve information stored in the TableList by name.
 * \param in      Pointer to TableList.
 * \param name    The name of the table to enter.
 * \param version Version number of table, >= 0 -> highest value found and this
 *                value is returned.
 * \param table   ObitTable to reference, should be Unrefed when done.
 *                May be NULL if table but no object exists.
 * \param err     (output) Obit Error stack for information.
 * \return TRUE if found, else FALSE.
 */
gboolean 
ObitTableListGet(ObitTableList *in, 
		gchar* name, olong *version, ObitTable **table, ObitErr *err)
{
  ObitTableListElem *elem;
  gboolean out = FALSE;

  /* error checks */
  g_assert (ObitTableListIsA(in));
  g_assert (ObitErrIsA(err));
  if (err->error) return out; /* error condition exists? */
  g_assert (name != NULL);
  g_assert (table != NULL);
  g_assert (version != NULL);

  /* look up */
  *table = NULL;
  elem = ObitTableListFind(in, name, version);
  if (elem==NULL) return out;

  /* Is it there and valid? */
  out = (elem->table!=NULL) && (elem->table->ReferenceCount>0);

  /* Set reference to pointer */
  if (out) *table = ObitTableRef(elem->table);

  /* If table is NULL, then it's just not been opened */
  if (elem->table==NULL) out = TRUE;

  /* remove from list if bad */
  if ((*table!=NULL) && ((*table)->ReferenceCount<=0))
    ObitTableListRemove (in, elem->name, elem->version);

  return out;
} /* end ObitTableListGet */

/**
 * Retrieve information stored in the TableList by number.
 * If the TableList contains an actual pointer to the table,
 * a reference is returned (needs to be Unreffed).
 * If the table pointer is NULL, this is returned and the return value 
 * is FALSE (although the table described exists and can be 
 * instantiated with newObit?Table.
 * \param in      Pointer to TableList.
 * \param number  1-rel item number
 * \param name    The name of the table , must be gfreed when done.
 * \param version Version number of table, -9=>invalid table
 * \param table   ObitTable reference, should be Unrefed when done.
 *                May be NULL if table but no object exists.
 * \param err     (output) Obit Error stack for information.
 * \return TRUE   if table object exists in list, else FALSE, 
 *                i.e. needs to be instantiated
 */
gboolean 
ObitTableListGetNumber(ObitTableList *in, olong number,
		      gchar **name, olong *version, ObitTable **table, ObitErr *err)
{
  ObitTableListElem *elem;
  GSList *tmp;
  gboolean found = FALSE;
  olong i;
  gchar *routine = "ObitTableListGetNumber";

  /* error checks */
  g_assert (ObitTableListIsA(in));
  g_assert (ObitErrIsA(err));
  if (err->error) return found;
  g_assert (name != NULL);
  g_assert (version != NULL);
  if ((number<0) || (number>in->number)) { /* number valid */
    Obit_log_error(err, OBIT_Error, 
	"%s: Invalid item number %d, max %d", 
	 routine, number, in->number);
      return found;
  }

  *table = NULL;
  tmp = in->list;
  /* look up */
  for (i=1; i<number; i++) {
    if (i>=number) break;
    tmp = g_slist_next(tmp);
    if (tmp==NULL) break;  /* problem? */
  }
  if (tmp==NULL) { /* not found */
      Obit_log_error(err, OBIT_Error, 
        "%s: I appear to have been corrupted", routine);
      return found;
  }

  elem = (ObitTableListElem*)tmp->data;
  if (elem==NULL) { /* not found */
      Obit_log_error(err, OBIT_Error, 
        "%s: I appear to have been corrupted", routine);
      return found;
  }

  /* Set reference to pointer */
  *name    = g_strdup(elem->name);
  *version = elem->version;
  if ((ObitTableIsA(elem->table) && (elem->table->ReferenceCount>0)) ||
      (elem->table==NULL))
    *table = ObitTableRef(elem->table);
  else {
    /* remove from list if bad */
    if ((elem!=NULL) && (*table!=NULL))
      ObitTableListRemove (in, elem->name, elem->version);
    *table = NULL;
    *version = -9;
  }
  found = *table != NULL;
  return found;
} /* end ObitTableListGetNumber */

/**
 * Find highest version number of a given table type
 * \param in      Pointer to TableList.
 * \param name    The name of the table to be searched for;
 * \return highest version number, 0 if none found
 */
olong ObitTableListGetHigh(ObitTableList *in, gchar *name)
{
  GSList *tmp;
  ObitTableListElem *elem;
  olong high=0;

  /* error checks */
  g_assert (ObitTableListIsA(in));
  g_assert (name != NULL);

  /* loop through list */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (ObitTableListElem*)tmp->data;
    /* is this one? */
    if (!strcmp(elem->name, name)) high = MAX (high, elem->version);
    tmp = g_slist_next(tmp);
  }

  return high;
} /* end ObitTableListGetHigh */

/**
 * Remove the item with keyword name from the list.
 * Item is destroyed; does nothing if item not found.
 * \param in      Pointer to TableList.
 * \param name    The name of the table.
 * \param version Version number of table.
 */
void 
ObitTableListRemove (ObitTableList *in, gchar *name,  olong version)
{
  ObitTableListElem *elem;
  olong tmpVer;

  /* error checks */
  g_assert (ObitTableListIsA(in));
  g_assert (name != NULL);

  /* look up */
  tmpVer = version;
  elem = ObitTableListFind(in, name, &tmpVer);
  if (elem==NULL) return; /* nothing to do */

  /* remove from list */
  in->list = g_slist_remove(in->list, elem);
  in->number--; /* keep count */

  /* delete element */
  freeObitTableListElem(elem);

} /* end ObitTableListRemove */

/**
 *  Public: Print list to stderr
 * \param in      Pointer to TableList.
 * \param err     Obit error/message stack
 */
void ObitTableListPrint  (ObitTableList *in, ObitErr *err)
{
  ObitTableListElem *elem;
  GSList *tmp;
  olong i;

  g_assert (ObitTableListIsA(in));
  g_assert (ObitErrIsA(err));
  if (err->error) return;

  /* loop */
  tmp = in->list;
  i = 0;
  while (tmp) {
    elem = (ObitTableListElem*)tmp->data;
    i++;
    fprintf (stderr," %d %d %s\n",i,elem->version,elem->name);
    tmp = g_slist_next(tmp);
  }
} /* end  ObitTableListPrint */

/**
 *  Public: Check list
 * \param in      Pointer to TableList.to check
 * \param err     Obit error/message stack, sets error if invalid
 */
void ObitTableListCheck  (ObitTableList *in, ObitErr *err)
{
  ObitTableListElem *elem;
  GSList *tmp;
  olong i;
  gboolean bad = FALSE;
  gchar *routine = "ObitTableListCheck";

  g_assert (ObitTableListIsA(in));
  g_assert (ObitErrIsA(err));
  if (err->error) return;

  /* loop */
  tmp = in->list;
  i = 1;
  while (tmp) {
    elem = (ObitTableListElem*)tmp->data;
    i++;
    if (elem->table) {
      bad = elem->table->ObitId!=OBIT_ID;  /* table object bad? */
      bad = bad || 
	((elem->table->ReferenceCount<0) || (elem->table->ReferenceCount>10000));
    }
    if (bad) {
      Obit_log_error(err, OBIT_Error,
		     "%s: I appear to have been corrupted, item %d", routine, i);
      return;
    }
    tmp = g_slist_next(tmp);
  }
} /* end  ObitTableListCheck */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitTableListClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitTableListClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitTableListClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitTableListClassInfoDefFn (gpointer inClass)
{
  ObitTableListClassInfo *theClass = (ObitTableListClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitTableListClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTableListClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitTableListGetClass;
  theClass->newObit       = (newObitFP)newObitTableList;
  theClass->ObitCopy      = (ObitCopyFP)ObitTableListCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitTableListClear;
  theClass->ObitInit      = (ObitInitFP)ObitTableListInit;
  theClass->ObitTableListPut = (ObitTableListPutFP)ObitTableListPut;
  theClass->ObitTableListGet = (ObitTableListGetFP)ObitTableListGet;
  theClass->ObitTableListGetNumber = 
    (ObitTableListGetNumberFP)ObitTableListGetNumber;
  theClass->ObitTableListRemove    = 
    (ObitTableListRemoveFP)ObitTableListRemove;

} /* end ObitTableListClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitTableListInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableList *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->list = NULL;
  in->number = 0;
} /* end ObitTableListInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitTableListClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableList *in = inn;
  GSList *tmp;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  /* loop through list deleting elements */
  tmp = in->list;
  while (tmp!=NULL) {
    freeObitTableListElem(tmp->data); tmp->data = NULL;
    tmp = g_slist_next(tmp);
  }

  /* delete members */
  g_slist_free(in->list);
  
 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitTableListClear */

/** 
 * Find the element with a given name in an ObitTableList
 * \param in      Pointer to TableList.
 * \param name    The label (keyword) of the element to change.
 * \param version Version number <=0 -> highest and returned.
 * \return pointer to ObitTableListElem or NULL if not found.
 */
static ObitTableListElem* 
ObitTableListFind(ObitTableList *in, gchar *name, olong *version)
{
  GSList *tmp;
  ObitTableListElem *elem, *outelem = NULL;
  olong outVer;
  gboolean highVer, match = FALSE;

  /* error checks */
  g_assert (ObitTableListIsA(in));
  g_assert (name != NULL);

  highVer = (*version <= 0);
  outVer = MAX (0, *version); /* output version number */
  /* loop through list testing elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (ObitTableListElem*)tmp->data;
    /* check if this is a match */
    if (!strcmp(elem->name, name)) { /* name matches */
      if (outVer==elem->version) return elem; /* definitely the one */
      /* Is this a possible match (version<=0) */
      if (highVer && (elem->version>outVer)) { /* could be the one */
	outelem = elem;
	outVer  = elem->version;
	match   = TRUE;
      }
    }
    tmp = g_slist_next(tmp); /* next element */
  }

  /* return version? */
  if (highVer && match) *version = outVer;

  /* Make sure table object valid */
  if ((outelem!=NULL) && ((ObitTableIsA(outelem->table)) && 
      (outelem->table->ReferenceCount<=0))) {
    /* remove from list if bad */
    if (outelem) ObitTableListRemove (in, outelem->name, outelem->version);
    outelem = NULL;
  }

  return outelem;
} /* end ObitTableListFind */

/*++++++++++++++++ ObitTableListElem functions ++++++++++++++++++*/

/**
 * Constructor
 * \param label   Name of table
 * \param version Version number of table.
 * \param table   Pointer to ObitTable object. (may be NULL).
 * \return the new ObitTableListElem object.
 */
static ObitTableListElem*
newObitTableListElem (gchar *name, olong version, ObitTable *table)
 { 
  ObitTableListElem *me; 

  /* error checks */
  g_assert (name != NULL);

  me =  g_malloc0(sizeof(ObitTableListElem));
  me->name    = g_strdup(name);
  me->version = version;
  if (table==NULL) me->table = NULL;
  else me->table = ObitTableRef(table);

  return me; 
} /* end newObitTableListElem */ 
    
/**
 * Copy Constructor
 * \param label   Name of table
 * \param version Version number of table.
 * \param table   Pointer to ObitTable object (may be NULL).
 * \return the new ObitTableListElem object.
 */
static ObitTableListElem* ObitTableListElemCopy (ObitTableListElem* in)
 { 
  ObitTableListElem *out; 

  /* error checks */
  g_assert (in != NULL);

  out = newObitTableListElem (in->name, in->version, in->table);

  return out; 
} /* end ObitTableListElemCopy */ 
  
/**
 * Destructor 
 * \param in Object to delete
 */
static void freeObitTableListElem (ObitTableListElem* in)
{ 
  g_assert (in != NULL);

  /* deallocate */
  if (in->name) g_free (in->name);
  if (in->table) in->table = ObitTableUnref(in->table);
  g_free (in);
}/* end of freeObitTableListElem  */ 
  
