/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2010                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include <string.h>
#include "ObitThread.h"

/** name of the class defined in this file */
static gchar *myClassName = "ObitThread";

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitThreadClassInfo myClassInfo = {FALSE};

/**
 * \file ObitThread.c
 * ObitThread (for multi threading) class function definitions.
 */
/*---------------Public functions---------------------------*/

/**
 * Construct Object.
 * \return pointer to object created.
 */
ObitThread* newObitThread (void)
{
  ObitThread* me;
#ifdef OBIT_THREADS_ENABLED
  GStaticRWLock *outRWLock=NULL;
#endif  /* OBIT_THREADS_ENABLED */

  /* Has g_threads been initialized? */
  if (myClassInfo.initialized==FALSE) {
#ifdef OBIT_THREADS_ENABLED
    g_thread_init(NULL);
#endif  /* OBIT_THREADS_ENABLED */
    myClassInfo.initialized = TRUE;  /* Now initialized */
    myClassInfo.haveThreads = FALSE; /* Don't know about threading yet */
    myClassInfo.nProcessor = 0;      /* Don't know about threading yet */
  }

  /* allocate structure */
  me = g_malloc0(sizeof(ObitThread));

  /* initialize */
  strncpy (me->className, myClassName, 11);
  me->ReferenceCount = 1;
  me->myClassInfo = &myClassInfo;
  me->pool     = NULL;
  me->queue    = NULL;
  me->myMutex  = NULL;
  me->myRWLock = NULL;
#ifdef OBIT_THREADS_ENABLED
  me->myMutex  = g_mutex_new(); 
  outRWLock = g_malloc0(sizeof(GStaticRWLock));
  g_static_rw_lock_init(outRWLock);
  me->myRWLock = outRWLock;
#endif  /* OBIT_THREADS_ENABLED */

  return me;
} /* end newObitThread */

/**
 * Destroy object, deallocation resources.
 * \param in Pointer to object to be destroyed.
 * \return Null pointer.
 */
ObitThread* freeObitThread (ObitThread *in)
{
  /* error checks */
  g_assert (ObitThreadIsA(in));

  /* Delete members */
#ifdef OBIT_THREADS_ENABLED
  if (in->pool)  g_thread_pool_free(in->pool, TRUE, FALSE);
  if (in->queue) g_async_queue_unref (in->queue);
  /*g_mutex_unlock(in->myMutex); *//* Make sure unlocked */
  g_mutex_free(in->myMutex); in->myMutex = NULL;
  /*g_static_rw_lock_reader_unlock(in->myRWLock); *//* Make sure unlocked */
  /*g_static_rw_lock_writer_unlock(in->myRWLock); *//* Make sure unlocked */
  g_free(in->myRWLock); in->myRWLock = NULL;
#endif  /* OBIT_THREADS_ENABLED */

  /* deallocate */
  g_free(in);
  return NULL;
} /* end freeObitThread */


/**
 * Copy constructor.
 * \param in Pointer to object to be copied.
 * \return Pointer to new object.
 */
ObitThread* ObitThreadCopy (ObitThread* in)
{
  ObitThread* out;

  /* error checks */
  g_assert (ObitThreadIsA(in));

  /* allocate output */
  out = newObitThread();

  return out;
} /* end ObitThreadCopy */

/**
 * Reference pointer, Update reference count, return pointer.
 * \param in Pointer to object to be linked.
 * \return Pointer to object.
 */
ObitThread* ObitThreadRef (ObitThread* in)
{
  /* error checks */
  g_assert (ObitThreadIsA(in));

  /* increment reference count */
  in->ReferenceCount++;
  return in;
} /* end ObitThreadRef */

/**
 * Unreference pointer, update reference count and 
 * destroy if zero.
 * This function should be used to dismiss an object.
 * \param in Pointer to object to be unlinked.
 * \return Null Pointer.
 */
ObitThread* ObitThreadUnref (ObitThread* in)
{
  /* error checks */
  if (in==NULL) return NULL;
  g_assert (ObitThreadIsA(in));

  /* decrement reference count, delete if non positive */
  in->ReferenceCount--;

  if (((in->ReferenceCount)) > 1) return NULL;

  if (in->ReferenceCount<=0) freeObitThread(in);

  return NULL; /* new value for pointer */
} /* end ObitThreadUnref */

/**
 * Lock the object so no other thread can access it.
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in Pointer to object to be locked.
 */
void ObitThreadLock (ObitThread *in)
{
  /* error checks 
  g_assert (ObitThreadIsA(in));*/

  /* Lock */
#ifdef OBIT_THREADS_ENABLED
  g_mutex_lock(in->myMutex);
#endif  /* OBIT_THREADS_ENABLED */

} /* end ObitThreadLock */

/**
 * Attempts to lock the object so no other thread can access it.
 * If successful returns TRUE, if the object is already locked, 
 * returns FALSE immediatly.
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in Pointer to object to be locked.
 * \return TRUE if successful
 */
gboolean ObitThreadTryLock (ObitThread *in)
{
  gboolean out = TRUE;
  /* error checks
  g_assert (ObitThreadIsA(in)); */

  /* Lock */
#ifdef OBIT_THREADS_ENABLED
  out = g_mutex_trylock(in->myMutex);
#endif  /* OBIT_THREADS_ENABLED */
  return out;
} /* end ObitThreadTryLock */

/**
 * Unlock the object so other threads can access it.
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in Pointer to object to be unlocked.
 */
void ObitThreadUnlock (ObitThread *in)
{
  /* error checks
  g_assert (ObitThreadIsA(in)); */

  /* Lock */
#ifdef OBIT_THREADS_ENABLED
  g_mutex_unlock(in->myMutex);
#endif  /* OBIT_THREADS_ENABLED */

} /* end ObitThreadUnlock */

/**
 * Lock the RWLock object for read so no other thread can modify it.
 * Other threads are allowed read access
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in Pointer to object to be locked.
 */
void ObitThreadRWReadLock (ObitThread *in)
{
  /* error checks
  g_assert (ObitThreadIsA(in)); */

  /* Lock */
#ifdef OBIT_THREADS_ENABLED
  g_static_rw_lock_reader_lock(in->myRWLock);
#endif  /* OBIT_THREADS_ENABLED */

} /* end ObitThreadRWReadLock */

/**
 * Attempts to lock the RWLock object for read so no other thread
 * can modify it.
 * If successful returns TRUE, if the object is already locked, 
 * returns FALSE immediatly.
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in Pointer to object to be locked.
 * \return TRUE if successful
 */
gboolean ObitThreadRWReadTryLock (ObitThread *in)
{
  gboolean out = TRUE;
  /* error checks
  g_assert (ObitThreadIsA(in)); */

  /* Lock */
#ifdef OBIT_THREADS_ENABLED
  out = g_static_rw_lock_reader_trylock(in->myRWLock);
#endif  /* OBIT_THREADS_ENABLED */
  return out;
} /* end ObitThreadRWReadTryLock */

/**
 * Unlock the RWLock object for read 
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in Pointer to object to be unlocked.
 */
void ObitThreadRWReadUnlock (ObitThread *in)
{
  /* error checks
  g_assert (ObitThreadIsA(in)); */

  /* Lock */
#ifdef OBIT_THREADS_ENABLED
  g_static_rw_lock_reader_unlock(in->myRWLock);
#endif  /* OBIT_THREADS_ENABLED */

} /* end ObitThreadRWReadUnlock */

/**
 * Lock the RWLock object for write so no other 
 * thread can access it.
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in Pointer to object to be locked.
 */
void ObitThreadRWWriteLock (ObitThread *in)
{
  /* error checks 
  g_assert (ObitThreadIsA(in));*/

  /* Lock */
#ifdef OBIT_THREADS_ENABLED
  g_static_rw_lock_writer_lock(in->myRWLock);
#endif  /* OBIT_THREADS_ENABLED */

} /* end ObitThreadRWWriteLock */

/**
 * Attempts to lock the RWLock object for write so no other thread 
 * can modify it.
 * If successful returns TRUE, if the object is already locked, 
 * returns FALSE immediatly.
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in Pointer to object to be locked.
 * \return TRUE if successful
 */
gboolean ObitThreadRWWriteTryLock (ObitThread *in)
{
  gboolean out = TRUE;
  /* error checks
  g_assert (ObitThreadIsA(in)); */

  /* Lock */
#ifdef OBIT_THREADS_ENABLED
  out = g_static_rw_lock_writer_trylock(in->myRWLock);
#endif  /* OBIT_THREADS_ENABLED */
  return out;
} /* end ObitThreadRWWriteTryLock */

/**
 * Unlock the RWLock object for write so other threads can modify it.
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in Pointer to object to be unlocked.
 */
void ObitThreadRWWriteUnlock (ObitThread *in)
{
  /* error checks 
  g_assert (ObitThreadIsA(in));*/

  /* Lock */
#ifdef OBIT_THREADS_ENABLED
  g_static_rw_lock_writer_unlock(in->myRWLock);
#endif  /* OBIT_THREADS_ENABLED */

} /* end ObitThreadRWWriteUnlock */

/**
 * Wait for a specified thread to finish
 * \param in Pointer to object specifying thread to wait for..
 */
void ObitThreadJoin  (ObitThread *in)
{
  /* error checks */
  g_assert (ObitThreadIsA(in));

  /* stubbed */
  g_error ("ObitThreadJoin is not defined - write me");
  g_assert_not_reached(); /* Function not yet defined */
} /* end ObitThreadJoin */

/**
 * Determines if the input object is a member of this class
 * \param in Pointer to object to test.
 * \return TRUE if member else FALSE.
 */
gboolean ObitThreadIsA (ObitThread* in)
{
  gboolean out;

  /* error checks */
  g_assert (in != NULL);

  /* compare class name member */
  out = !strcmp(in->className, myClassName);

  return out;
} /* end ObitThreadIsA */

/**
 * Determines if the Threading is enabled
 * \param in Pointer to a thread object.
 * \return TRUE if threads are enabled
 */
gboolean ObitThreadHaveThreads (ObitThread* in)
{
  /* error checks */
  g_assert (in != NULL);

  return in->myClassInfo->haveThreads;
} /* end ObitThreadHaveThreads */

/**
 * Initializes threading based on the parameter "nThreads" in myInput.
 * No change is made if compilations did not have OBIT_THREADS_ENABLED
 * \param myInput  an ObitInfoList possible containing
 * \li "nThreads"   OBIT_long (1,1,1) Number of threads to attempt per pool.
 */
void ObitThreadInit (ObitInfoList *myInput)
{
  ObitThread *thread=NULL;
  olong nThreads;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  /* Create temporary thread */
  thread = newObitThread ();

  /* Number of threads from infoList */
  nThreads = 1;
  ObitInfoListGetTest(myInput, "nThreads", &type, dim, &nThreads);

  /* Init */
  ObitThreadAllowThreads (thread, nThreads);
  freeObitThread (thread);  /* Cleanup */

} /* end ObitThreadInit*/

/**
 * Allows threading and sets number of threads to use for processing.
 * No change is made if compilations did not have OBIT_THREADS_ENABLED
 * \param in Pointer to a ObitThread object.
 * \param nProc Number of threads
 */
void ObitThreadAllowThreads (ObitThread* in, olong nThreads)
{

  /* error checks */
  g_assert (in != NULL);

  /* Copy inputs if compiled with threads */
#ifdef OBIT_THREADS_ENABLED
  in->myClassInfo->haveThreads = TRUE;
  in->myClassInfo->nProcessor  = nThreads;
#endif  /* OBIT_THREADS_ENABLED */

} /* end ObitThreadAllowThreads */

/**
 * Tells number of processors available for multithreading
 * \param in Pointer to a Thread object
 * \return number of processors which can be threaaded, <=0 => no threading.
 */
olong ObitThreadNumProc (ObitThread* in)
{
  olong out;

  /* error checks */
  g_assert (in != NULL);

  /* set if threading is allowed */
  if (in->myClassInfo->haveThreads)
    out =  in->myClassInfo->nProcessor;
  else
    out = 0;

  return out;
} /* end ObitThreadNumProc */

/** 
 * Initializes Thread Pool and asynchronous queue
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in        Pointer to Thread object
 * \param nthreads  Number of threads to create/run
 * \param func      Function to call to start thread
 * \param args      Array of argument function pointers
*/
void ObitThreadPoolInit (ObitThread* in, olong nthreads,
			 ObitThreadFunc func, gpointer **args)
{
#ifdef OBIT_THREADS_ENABLED
  in->pool  = g_thread_pool_new ((GFunc)func, args, nthreads, FALSE, NULL);
  in->queue = g_async_queue_new ();
#endif  /* OBIT_THREADS_ENABLED */
} /* end ObitThreadPoolInit */

/**
 * Loops over a set of threads, all with same function call
 * If nthreads=1 or threading not allowed, routine called directly.
 * Waits for operations to finish before returning, 10 min timeout.
 * Initializes Thread pool and asynchronous queue (ObitThreadPoolInit)
 * if not already done.
 * When threaded operations are finished, call ObitThreadPoolFree to release
 * Thread pool.
 * \param in Pointer to object
 * \param nthreads  Number of threads to create/run
 * \param func      Function to call to start thread
 *  func should call ObitThreadPoolDone to indicate completion.
 * \param args      Array of argument function pointers
 * \return TRUE if OK else FALSE.
 */
gboolean ObitThreadIterator (ObitThread* in, olong nthreads,
			     ObitThreadFunc func, gpointer **args)
{
  gboolean out = TRUE;
  GTimeVal end_time;
  glong add_time;
  olong i;

  /* error checks */
  g_assert (in != NULL);

  /* If only 1 don't use thread */
  if (nthreads==1) {
    if ((func)(args[0])) {
      /* error check?*/
      out = FALSE;
    }
    return out;
  }
   
  /* If no threading don't use thread - call sequentially */
  if (!in->myClassInfo->haveThreads) {
    for (i=0; i<nthreads; i++) {
      if ((func)(args[i])) {
	/* error check?*/
	out = FALSE;
      }
    }
    return out;
  }

  /* Make sure pool is using the correct function */
  if ((in->pool) && (((GFunc)in->pool->func)!=((GFunc)func))) {
     g_thread_pool_free(in->pool, FALSE, TRUE);  /* Be patient */
     in->pool = NULL;
     if (in->queue) g_async_queue_unref (in->queue);
     in->queue = NULL;
  }
   
  /* Threading allowed - has the Thread Pool been started? */
  if (in->pool==NULL) ObitThreadPoolInit (in, nthreads, func, args);

  /* Submit jobs to the pool */
  for (i=0; i<nthreads; i++) {
    g_thread_pool_push (in->pool, args[i], NULL);
  }

  /* Wait for them to finish, expects each to send a message to the asynchronous 
   queue iff they finish */
  /* 1 min. timeout */
  add_time = 600 * 1000000;
  g_get_current_time (&end_time);
  g_time_val_add (&end_time, add_time); /* add timeout in microseconds */
  for (i=0; i<nthreads; i++) {
    g_async_queue_timed_pop (in->queue, &end_time);
 }

  return out;
} /* end ObitThreadIterator */

/** 
 * Indicates that a thread function is done by sending a message to the 
 * asynchronous queue on in
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in        Pointer to Thread object
 * \param arg       Pointer to message (CANNOT be NULL)
 */
void ObitThreadPoolDone (ObitThread* in, gpointer arg)
{
#ifdef OBIT_THREADS_ENABLED
  if (!in->queue) return;
  g_async_queue_push (in->queue, arg);
#endif  /* OBIT_THREADS_ENABLED */
} /* end ObitThreadPoolDone */

/** 
 * Shuts down Thread Pool and  asynchronous queue on in
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in        Pointer to Thread object
 */
void ObitThreadPoolFree (ObitThread* in)
{
#ifdef OBIT_THREADS_ENABLED
  if (in->pool)  g_thread_pool_free(in->pool, TRUE, TRUE); 
  in->pool = NULL;
  if (in->queue) g_async_queue_unref (in->queue);
  in->queue = NULL;
#endif  /* OBIT_THREADS_ENABLED */
} /* end ObitThreadPoolFree */

/** 
 * Waits for single thread 
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in        Pointer to Thread object
 * \param func      Function to call to start thread
 * \param arg       Function argument pointer
 */
void ObitThreadStart1 (ObitThread* in, ObitThreadFunc func, gpointer args)
{
#ifdef OBIT_THREADS_ENABLED
  in->singleThread = (GThread*)g_thread_create ((GThreadFunc)func, args, TRUE, NULL);
#endif  /* OBIT_THREADS_ENABLED */
} /* end ObitThreadStart1 */

/** 
 * Waits for single thread 
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in        Pointer to Thread object
 * \return pointer to object returned by function
 */
gpointer ObitThreadJoin1 (ObitThread* in)
{
#ifdef OBIT_THREADS_ENABLED
  return g_thread_join(in->singleThread);
#endif  /* OBIT_THREADS_ENABLED */
} /* end ObitThreadJoin1 */

/** 
 * Initializes Thread asynchronous queue
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in        Pointer to Thread object
*/
void ObitThreadQueueInit (ObitThread* in)
{
#ifdef OBIT_THREADS_ENABLED
  in->queue = g_async_queue_new ();
#endif  /* OBIT_THREADS_ENABLED */
} /* end ObitThreadQueueInit */

/** 
 * Check for messages in queue
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in        Pointer to Thread object
 * \param add_time  timeout in  microseconds, <= -> forever
 */
gpointer ObitThreadQueueCheck (ObitThread* in, olong add_time)
{
#ifdef OBIT_THREADS_ENABLED
  GTimeVal end_time;
  
  g_get_current_time (&end_time);
  g_time_val_add (&end_time, add_time);
  if (add_time>0) return g_async_queue_timed_pop (in->queue, &end_time); 
  else return g_async_queue_pop (in->queue);
#endif  /* OBIT_THREADS_ENABLED */
} /* end ObitThreadQueueInit */

/** 
 * Shuts down asynchronous queue on in
 * Noop unless compiled with OBIT_THREADS_ENABLED
 * \param in        Pointer to Thread object
 */
void ObitThreadQueueFree (ObitThread* in)
{
#ifdef OBIT_THREADS_ENABLED
  if (in->queue) g_async_queue_unref (in->queue);
  in->queue = NULL;
#endif  /* OBIT_THREADS_ENABLED */
} /* end ObitThreadQueueFree */

