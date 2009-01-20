/* $Id$      */
/* g_hash_table version */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
#include <string.h>
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitMem.c
 * ObitMem Module definition file.
 *
 * Obit Memory management class
 * Use compiler switch -DFASTOBITMEM to disable checking or
 * -DMEMWATCH to use MEMWATCH debugging
 */

/*--------------Class definitions-------------------------------------*/
typedef struct { 
  /**  Pointer */
  gpointer mem;
  /** size in bytes */
  gulong size;
  /** name (20 char) */
  gchar name[21];
}  memTableElem; 

/**
 * ObitMem Class structure.
 *
 * This class controls Obit memory management
 */  
typedef struct  {
  /** Have I been initialized  */
  gboolean init;
  /** Threading info member object  */
  ObitThread *thread;
  /** glib singly linked list for registered allocated memory */
  GHashTable* memTable;
  /** How many entries */
  gulong number;
} ObitMemClassInfo; 

/**
 * Structure to use in validity testing.
 */  
typedef struct  {
  /** Has a valid entry been found?  */
  gboolean foundIt;
  /** memory pointer to test  */
  gpointer mem;
} ObitMemValidStruc; 

/** Switch to use fast memory allocation and no accountability */
/* #define FASTOBITMEM */


/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitMemClassInfo myClassInfo = {FALSE};

/**
 * Class structure to use for testing memory validity
 */
static ObitMemValidStruc myValidTest = {FALSE, NULL};

/*---------------Private function prototypes----------------*/
/** Private: Create a memTableElem. */
static memTableElem* newmemTableElem (gpointer mem, gulong size, gchar* name);

/** Private: Delete a memTableElem. */
static void freememTableElem (memTableElem *in);

/** Private: Add an memTableElem to the list. */
static void memTableAdd (ObitMemClassInfo *in, memTableElem *elem);

/** Private: Remove an memTableElem from the list. */
static void memTableRemove (ObitMemClassInfo *in, memTableElem *elem);

/** Private: Find item in a list */
static memTableElem*  memTableFind(ObitMemClassInfo *in, gpointer mem);

/** Private: Is in valid memory? */
static void memTableInValid (gpointer key, gpointer inn, gpointer ttest);

/** Private: Print an memTableElem */
static void  memTablePrint(gpointer key, gpointer inn, gpointer file);

/** Private: Accumulate bytes in an memTableElem */
static void  memTableSum(gpointer key, gpointer inn, gpointer dcount);

/*---------------Public functions---------------------------*/
#ifndef MEMWATCH
/**
 * Allocates a block of memory and returns pointer to it
 * Initializes class if needed on first call.
 * Implementation used g_malloc
 * \param size Number of bytes requested
 * \return pointer to allocated memory, NULL on failure
 */
gpointer ObitMemAlloc (gulong size)
{
#ifdef FASTOBITMEM /* Fast allocation */
  return g_malloc(size);
#else              /* Accountability */
  
  gpointer out = NULL;
  memTableElem* elem = NULL;

  /* Class initialization if needed */
  if (!myClassInfo.init) ObitMemClassInit();

  /* nothing asked for - nothing given */
  if (size==0) return NULL;
  
  /* allocate */
  out = g_malloc(size);

  /* Create list object */
  elem = newmemTableElem (out, size, NULL);

   /* Lock object aginst other threads */
  ObitThreadLock(myClassInfo.thread);

  /* add to list */
  memTableAdd (&myClassInfo, elem);

  /* Unlock object */
  ObitThreadUnlock(myClassInfo.thread);

   return out;
#endif /* FASTOBITMEM  */
} /* end ObitMemAlloc */

/**
 * Allocates a block of memory, zero fills and returns pointer to it
 * Initializes class if needed on first call.
 * Implementation uses g_malloc0
 * \param size Number of bytes requested
 * \return pointer to allocated memory, NULL on failure
 */
gpointer ObitMemAlloc0 (gulong size)
{
#ifdef FASTOBITMEM /* Fast allocation */
  return g_malloc0(size);
#else              /* Accountability */
  gpointer out = NULL;
  memTableElem* elem = NULL;

  /* Class initialization if needed */
  if (!myClassInfo.init) ObitMemClassInit();
  
  /* nothing asked for - nothing given */
  if (size==0) return NULL;
  
  /* allocate */
  out = g_malloc0(size);

  /* Create list object */
  elem = newmemTableElem (out, size, NULL);

   /* Lock object aginst other threads */
  ObitThreadLock(myClassInfo.thread);

  /* add to list */
  memTableAdd (&myClassInfo, elem);

  /* Unlock object */
  ObitThreadUnlock(myClassInfo.thread);

  return out;
#endif /* FASTOBITMEM  */
} /* end ObitMemAlloc0 */

/**
 * Allocates a block of memory and returns pointer to it, giving name.
 * Initializes class if needed on first call.
 * Implementation used g_malloc
 * \param size Number of bytes requested
 * \param name Name for entry, up to 20 char.  Useful for debugging.
 * \return pointer to allocated memory, NULL on failure
 */
gpointer ObitMemAllocName (gulong size, gchar *name)
{
#ifdef FASTOBITMEM /* Fast allocation */
  return g_malloc(size);
#else              /* Accountability */
  gpointer out = NULL;
  memTableElem* elem = NULL;

  /* Class initialization if needed */
  if (!myClassInfo.init) ObitMemClassInit();

  /* nothing asked for - nothing given */
  if (size==0) return NULL;
  
  /* allocate */
  out = g_malloc(size);

  /* Create list object */
  elem = newmemTableElem (out, size, name);

   /* Lock object aginst other threads */
  ObitThreadLock(myClassInfo.thread);

  /* add to list */
  memTableAdd (&myClassInfo, elem);

  /* Unlock object */
  ObitThreadUnlock(myClassInfo.thread);

   return out;
#endif /* FASTOBITMEM  */
} /* end ObitMemAllocName */

/**
 * Allocates a block of memory, zero fills and returns pointer to it
 * Initializes class if needed on first call.
 * Implementation uses g_malloc0
 * \param size Number of bytes requested
 * \param name Name for entry, up to 20 char.  Useful for debugging.
 * \return pointer to allocated memory, NULL on failure
 */
gpointer ObitMemAlloc0Name (gulong size, gchar *name)
{
#ifdef FASTOBITMEM /* Fast allocation */
  return g_malloc0(size);
#else              /* Accountability */
  gpointer out = NULL;
  memTableElem* elem = NULL;

  /* Class initialization if needed */
  if (!myClassInfo.init) ObitMemClassInit();

  /* nothing asked for - nothing given */
  if (size==0) return NULL;
  
  /* allocate */
  out = g_malloc0(size);

  /* Create list object */
  elem = newmemTableElem (out, size, name);

   /* Lock object aginst other threads */
  ObitThreadLock(myClassInfo.thread);

  /* add to list */
  memTableAdd (&myClassInfo, elem);

  /* Unlock object */
  ObitThreadUnlock(myClassInfo.thread);
  /* DEBUG
  if (elem->mem==(gpointer)0x82d3c88) {
    fprintf (stderr, "alloc %9p %8ld %s\n", elem->mem, elem->size, elem->name);
    } */

  return out;
#endif /* FASTOBITMEM  */
}
 /* end ObitMemAlloc0Name */

/**
 * Reallocates a block of memory and returns pointer to it
 * Nothing is changed if the old and new memory sizes are the same
 * Implementation uses g_realloc
 * \param mem Pointer to old memory, if null, allocates anonymous slot
 * \param size Number of bytes requested
 * \return pointer to allocated memory, NULL on failure
 */
gpointer ObitMemRealloc (gpointer mem, gulong size)
{
 #ifdef FASTOBITMEM /* Fast allocation */
  return g_realloc (mem, size);
#else              /* Accountability */
  memTableElem* elem = NULL;
   gpointer out = NULL;
   gboolean isNew = FALSE;

  /* error checks */
  g_assert (myClassInfo.init);
  
  /* Lock object aginst other threads */
  ObitThreadLock(myClassInfo.thread);

  /* is it in the list? */
  elem =  memTableFind (&myClassInfo, mem);

  if (elem==NULL) isNew = TRUE;
  else {  /* previously OK */
    out = elem->mem;
    
    /* Is it the same size? If so just return */
    if (elem->size==size) return out;
  }

  /* Realloc */
  out = g_realloc (out, size);

  /* Is this a new entry */
  if (isNew) {
    /* Create list object */
    elem = newmemTableElem (out, size, NULL);
    /* add to list */
    memTableAdd (&myClassInfo, elem);
  } else { /* just reset pointer, size */
    elem->mem = out;
    elem->size = size;
  }

  /* Unlock object */
  ObitThreadUnlock(myClassInfo.thread);

  return out;
#endif /* FASTOBITMEM  */
} /* end ObitMemRealloc */

/**
 * Deallocates a block of memory and returns NULL.
 * This is a NOP if mem is not to a valid block of memory.
 * Note: ONLY FREE MEMORY ALLOCATED BY ObitMem!!!
 * Implementation used g_free
 * \param mem Pointer to memory to be freed.  
 * \return pointer to allocated memory, NULL on failure
 */
gpointer ObitMemFree (gpointer mem)
{
 #ifdef FASTOBITMEM /* Fast allocation */
  g_free(mem);
  return NULL;
#else              /* Accountability */
 memTableElem* elem = NULL;

  /* error checks */
  g_assert (myClassInfo.init);
  if (mem==NULL) return NULL;  /* needs to point to something */

  /* Lock object aginst other threads */
  ObitThreadLock(myClassInfo.thread);

 /* is it in the list? */
  elem =  memTableFind (&myClassInfo, mem);
  if (elem==NULL) return NULL; /* bag it if it's not there */

  /* DEBUG
  if (elem->mem==(gpointer)0x82d3c88) {
    fprintf (stderr, " free %9p %8ld %s\n", elem->mem, elem->size, elem->name);
  } */

  /* Can it */
  memTableRemove (&myClassInfo, elem);

  /* Unlock object */
  ObitThreadUnlock(myClassInfo.thread);

  /* Deallocate memory */
  g_free(elem->mem); elem->mem = NULL;

  /* destroy list element */
  freememTableElem (elem);

  return NULL;
#endif /* FASTOBITMEM  */
}  /* end ObitMemFree */

/**
 * Determines if a pointer value is inside a valid allocated memory block.
 * \param mem Pointer to memory to be tested.  
 * \return TRUE if in allocated block, else FALSE
 */
gboolean ObitMemValid (gpointer mem)
{
#ifdef FASTOBITMEM /* Fast allocation */
  return TRUE;
#else              /* Accountability */

  /* error check */
  g_assert (myClassInfo.init);

  /* Lock object aginst other threads */
  ObitThreadLock(myClassInfo.thread);

  /* loop through list testing elements */
  myValidTest.foundIt = FALSE;
  myValidTest.mem     = mem;
  g_hash_table_foreach (myClassInfo.memTable, memTableInValid, &myValidTest);
  
  /* Unlock object */
  ObitThreadUnlock(myClassInfo.thread);

  return myValidTest.foundIt;
#endif /* FASTOBITMEM  */
} /* end ObitMemValid */

#endif /* MEMWATCH */

/**
 * Print contents of ObitMem to file (e.g. stdout)
 * \param file to print to.
 */
void ObitMemPrint (FILE *file)
{
#ifdef MEMWATCH /* Using memwatch instead */
  return;  /* No can do */
#endif
#ifdef FASTOBITMEM /* Fast allocation */
  return;  /* No can do */
#else              /* Accountability */
  if (!myClassInfo.init) {
    fprintf (file,"Obit memory allocation system uninitialized\n");
    return;
  }
  fprintf (file,"Obit memory allocation has %d entries\n", 
	   myClassInfo.number);
  fprintf (file,"Address  size(byte) name\n");

  /* loop through list printing elements */
  g_hash_table_foreach (myClassInfo.memTable, memTablePrint, file);
#endif /* FASTOBITMEM  */
} /* end ObitMemPrint */


/**
 * Give Summary of contents of myClassInfo
 * \param number  Number of entries
 * \param total   Total memory allocated in MByte
 */
void ObitMemSummary (olong *number, olong *total)
{
#ifdef MEMWATCH /* Using memwatch instead */
  *number = 0;
  *total  = 0;
  return;  /* No can do */
#endif
#ifdef FASTOBITMEM /* Fast allocation */
  return;  /* No can do */
#else              /* Accountability */
  odouble count = 0.0;

  if (!myClassInfo.init) {
    *number = 0;
    *total  = 0;
    return;
  }
  /* Number */
  *number = myClassInfo.number;

  /* loop through list printing elements */
  g_hash_table_foreach (myClassInfo.memTable, memTableSum, &count);

  /* Total in MByte */
  *total = (olong)(0.5 + (count / (1024.0*1024.0)));
#endif /* FASTOBITMEM  */
} /* end ObitMemSummary */


/**
 * Initialize global ClassInfo Structure.
 */
 void ObitMemClassInit (void)
{

  if (myClassInfo.init) return;  /* only once */
  myClassInfo.init = TRUE;

  /* Threading object */
  myClassInfo.thread = newObitThread();

  /* initialize list */
  myClassInfo.memTable = g_hash_table_new (NULL, NULL);
  myClassInfo.number  = 0;

} /* end ObitMemClassInit */

/*---------------Private functions--------------------------*/

/**
 * memTableElem Constructor
 * \param item  Data
 * \param name Name for entry, up to 20 char.  NULL = none.
 * \return the new  object.
 */
static memTableElem* newmemTableElem (gpointer mem, gulong size, gchar *name)
{
  memTableElem *out=NULL;

  out = g_malloc0(sizeof(memTableElem));
  out->mem  = mem;
  out->size = size;
  if (name) strncpy (out->name, name, 20); out->name[20] = 0;

  return out;
} /* end newmemTableElem */

/**
 * memTableElem Destructor 
 * \param in Object to delete
 */
static void freememTableElem (memTableElem *in)
{
  if (in) {
    in->mem = NULL; /* so it doesn't come back to haunt us */
    g_free(in);
  }
} /* end freememTableElem */

/**
 * Attach elem to table in
 * \param in   Object with table to add elem to
 * \param elem the element to add. MUST NOT be in list
 */
static void memTableAdd (ObitMemClassInfo *in, memTableElem *elem)
{
  memTableElem *tmp = NULL;

  /* Make sure it's not already in list */
  tmp =  memTableFind (&myClassInfo, elem->mem);
  if (tmp!=NULL) { /* trouble - die in a noisy fashion */
    g_error ("ObitMem: trying to allocate memory already allocated: %s new %s",
	     tmp->name, elem->name);
  }

  /* add to table */
  g_hash_table_insert (in->memTable, elem->mem, elem);
  in->number++;

} /* end memTableAdd */

/**
 * Remove elem from table in in
 * \param in   Object with table to remove elem from
 * \param elem the element to remove. MUST be in list
 */
static void memTableRemove (ObitMemClassInfo *in, memTableElem *elem)
{
  /* remove from table */
  if (!g_hash_table_steal (in->memTable, elem->mem)) {
    /* element not removed */
    g_error ("ObitMem: failed to deallocate memory: %s",
	     elem->name);
  }
  in->number--; /* keep count */

} /* end memTableRemove  */

/**
 * Find a pointer is in list in
 * \param in   Object with table to search
 * \param mem item to search for
 * \return pointer to element containing item, NULL if not found.
 */
static memTableElem* memTableFind (ObitMemClassInfo *in, gpointer mem)
{
  memTableElem *out = NULL;
  /* Lookup in hash table */
  out = g_hash_table_lookup (in->memTable, mem);

  /* Guard aginst false detections */
  if (out) 
    if (out->mem!=mem) out = NULL;

  return out;
} /* end memTableFind */

/**
 * Check if a pointer value is inside an allocated block (GHFunc)
 * This function is called from g_hash_table_foreach.
 * \param key   pointer to hash key
 * \param mem   memTableElem* object to test
 * \param test  ObitMemValidStruc* target/result structure 
 */
static void memTableInValid (gpointer key, gpointer inn, gpointer ttest)
{
  memTableElem* in = (memTableElem*)inn;
  ObitMemValidStruc* test = (ObitMemValidStruc*)ttest;

  /* if it's already found just return */
  if (test->foundIt) return;

  /* test if in is in this block */
  test->foundIt = (test->mem>=in->mem) && (test->mem<(in->mem+in->size));
     
} /* end memTableInValid */

/**
 * Print an memTableElem (GHFunc) 
 * This function is called from g_hash_table_foreach.
 * \param key   pointer to hash key
 * \param in    memTableElem* to print
 * \param file  FILE* to write to
 */
static void  memTablePrint(gpointer key, gpointer inn, gpointer file)
{
  memTableElem* in = (memTableElem*)inn;
  fprintf ((FILE*)file, "%9p %8ld %s\n", in->mem, in->size, in->name);
} /* end memTablePrint */

/**
 * AccumulatePrint an memTableElem (GHFunc) 
 * This function is called from g_hash_table_foreach.
 * \param key     pointer to hash key
 * \param in      memTableElem* to print
 * \param dcount  odouble* to accumulate to
 */
static void  memTableSum(gpointer key, gpointer inn, gpointer dcount)
{
  memTableElem* in = (memTableElem*)inn;
  odouble *count = (odouble*)dcount;

  /* Accumulate */
  *count += in->size;
} /* end memTableSum */

