/* $Id$       */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2008                                          */
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
#ifndef OBITTHREAD_H 
#define OBITTHREAD_H 
#include "ObitErr.h"
#include "ObitInfoList.h"
#include <glib.h>
#include <glib/gthread.h>

/**
 * \file ObitThread.h
 * ObitThread (for multi threading) class definition.
 *
 * This class supports multithreading and provides mutexes for absolute 
 * locking of associated resources as well as RWLocks to allow multiple
 * read accesses but only a single write (and no concurrent reads).
 * Needs OBIT_THREADS_ENABLED defined at compile time and the output of 
 * pkg-config --libs gthread-2.0 added to the libraries.
 */
/*--------------Class definitions-------------------------------------*/
/*-------------------Class Info--------------------------*/

/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any base class
 * (NULL if none) and function pointers.
 */
typedef struct  {
  /** Have I been initialized? */
  gboolean initialized;
  /** Are multiple threads allowed? */
  gboolean haveThreads;
  /** Number of processors for determining number of threads */
  olong nProcessor;
} ObitThreadClassInfo; 


/**
 * ObitThread Class Structure.
 *
 * This class is to facilitate multithreaded operation.
 */  
typedef struct {
  /** class name for verification */
  gchar className[12];
  /** Thread Id */
  gint32 id;
  /** Object reference count. */
  gint32 ReferenceCount; 
  /** class info */
  ObitThreadClassInfo *myClassInfo;
#ifdef OBIT_THREADS_ENABLED
  /** Pool of threads */
  GThreadPool *pool;
  /** Asynchronous message queue */
  GAsyncQueue *queue;
  /** Mutex to lock object */
  GMutex *myMutex;
  /** Read/write lock to lock object */
  GStaticRWLock *myRWLock;
#else  /* No threads */
  gpointer pool;
  gpointer queue;
  gpointer myMutex;
  gpointer myRWLock;
#endif  /* OBIT_THREADS_ENABLED */
} ObitThread;

/** Thread Function template pointer */
typedef GThreadFunc ObitThreadFunc;

/*---------------Public functions---------------------------*/
/** Public: Constructor. */
ObitThread* newObitThread (void);

/** Public: Destructor. */
ObitThread* freeObitThread (ObitThread *in);

/** Public: Copy Thread */
ObitThread* ObitThreadCopy (ObitThread* in);

/** Public: Reference Thread pointer */
ObitThread* ObitThreadRef (ObitThread* in);

/** Public: Unreference Thread pointer */
ObitThread* ObitThreadUnref (ObitThread* in);

/** Public: Lock (mutex) object */
void ObitThreadLock (ObitThread *in);

/** Public: Test Lock (mutex) object */
gboolean ObitThreadTryLock (ObitThread *in);

/** Public: Unlock (mutex) object  */
void ObitThreadUnlock (ObitThread *in);

/** Public: Lock (RWLock) object for reader */
void ObitThreadRWReadLock (ObitThread *in);

/** Public: Test Lock (RWLock) object  for reader */
gboolean ObitThreadRWReadTryLock (ObitThread *in);

/** Public: Unlock (RWLock) object for reader */
void ObitThreadRWReadUnlock (ObitThread *in);

/** Public: Lock (RWLock) object for writer */
void ObitThreadRWWriteLock (ObitThread *in);

/** Public: Test Lock (RWLock) object for writer */
gboolean ObitThreadRWWriteTryLock (ObitThread *in);

/** Public: Unlock (RWLock) object for writer */
void ObitThreadRWWriteUnlock (ObitThread *in);

/** Public: Join thread  */
void ObitThreadJoin  (ObitThread *in);

/** Public: Returns true if input is a ObitThread */
gboolean ObitThreadIsA (ObitThread* in);

/** Public: Returns true if Threads are enabled */
gboolean ObitThreadHaveThreads (ObitThread* in);

/** Public: Initialize Threading */
void ObitThreadInit (ObitInfoList *myInput);

/** Public: Allows Threads and sets number of processors which 
    can be multithreaded */
void ObitThreadAllowThreads (ObitThread* in, olong nThreads);

/** Public: Returns number of processors which can be multithreaded */
olong ObitThreadNumProc (ObitThread* in);

/** Public: Initializes Thread Pool */
void ObitThreadPoolInit (ObitThread* in, olong nthreads,
			 ObitThreadFunc func, gpointer **args);

/** Public: Runs multiple copies of a function in different threads */
gboolean ObitThreadIterator (ObitThread* in, olong nthreads,
			     ObitThreadFunc func, gpointer **args);

/** Public: Indicates that a thread function is done */
void ObitThreadPoolDone (ObitThread* in, gpointer arg);

/** Public: Shuts down Thread Pool */
void ObitThreadPoolFree (ObitThread* in);

#endif /* OBITTHREAD_H */ 

