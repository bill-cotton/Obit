/* $Id$        */
/*--------------------------------------------------------------------*/
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
#ifndef OBITMULTIPROC_H 
#define OBITMULTIPROC_H 

#include "ObitRPC.h"
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitFile.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitMultiProc.h
 *
 * ObitMultiProc template for classes derived from #Obit
 *
 * The ObitMultiProc class enables parallel processing using multiple 
 * asynchronous processes.  These processes are started by the master
 * process and execute FuncContainer, an XMLRPC based compute server.  
 * If multiple processing is not enabled (or only a single job is to be 
 * performed) the execution is sequential and in the main task address space.
 * When the master tasks starts the asynchronous FuncContainer processes, the 
 * critical aspects of the environment (e.g. AIPS directories) are passed.
 * Pieces of computing to be be done in a potentially parallel process are 
 * referred to as "jobs".
 * 
 *   Function arguments are passed to FuncContainer processes and return values, 
 * status and logging is returned by means of ObitInfoLists.  Functions for which 
 * parallel executions are enabled have a simple (ObitThreadFunc) interface which
 * is passed a single pointer to a structure containing the detailed arguments.
 * These functions should be callable either locally or by the appropriate RPC
 * function in a FuncContainer process.  The ObitInfoLists contain the parameters 
 * of the work to be performed.
 *
 *    Implementation is based on the availibity of GThreads.  The asynchronous 
 * processes are each started in a thread.  Threads are then used to manage the 
 * interface between the main task and the auxillary processes.
 * 
 * \section ObitMultiProcaccess Creators and Destructors
 * An ObitMultiProc will usually be created using ObitMultiProcCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitMultiProc should always be made using the
 * #ObitMultiProcRef function which updates the reference count in the object.
 * Then whenever freeing an ObitMultiProc or changing a pointer, the function
 * #ObitMultiProcUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/

/** Typedef for MultiProc enabled functions */
typedef ObitInfoList* (*ObitMultiProcFunc) (ObitInfoList *args, ObitErr *err);

/**
 * MultiProc function argument
 */  
typedef struct {
  /** RPC Client  */
  ObitRPC *client;
  /** URL for FuncContainer  */
  gchar *URL;
  /** Function name to call in FuncContainer for remote RPC call. */
  gchar *MPfunc;
  /** Function pointer to call in local address space,
      the same arguments and functionality as MPfunc. */
  ObitMultiProcFunc localfunc;
  /** Arguments as InfoList  */
  ObitInfoList *args;
  /** Return values as InfoList  */
  ObitInfoList *retVal;
  /** thread object  */
  ObitThread *thread;
  /** thread/stream number  */
  olong ithread;
  /** job number  */
  olong ijob;
  /** Obit error/logging object */
  ObitErr *err;
  /** Arbitrary user data */
  gpointer user_data;
} ObitMultiProcFuncArg;


/** ObitMultiProc Class structure. */
typedef struct {
#include "ObitMultiProcDef.h"   /* this class definition */
} ObitMultiProc;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitMultiProc
 * returns a ObitMultiProc*.
 * in = object to unreference
 */
#define ObitMultiProcUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitMultiProc.
 * returns a ObitMultiProc*.
 * in = object to reference
 */
#define ObitMultiProcRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitMultiProcIsA(in) ObitIsA (in, ObitMultiProcGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitMultiProcClassInit (void);

/** Public: Default Constructor. */
ObitMultiProc* newObitMultiProc (gchar* name);

/** Public: Create/initialize ObitMultiProc structures */
ObitMultiProc* ObitMultiProcCreate (gchar* name, olong njobs,
				    gchar *MPfunc, ObitMultiProcFunc localFunc,
				    ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitMultiProc* (*ObitMultiProcCreateFP) (gchar* name, olong njobs,
						 gchar *MPfunc, ObitMultiProcFunc localFunc,
						 ObitErr *err);

/** Start auxillary FuncContainer processes based on values in Inputs list */
void ObitMultiProcStart(ObitInfoList *myInput, ObitErr *err);

/** Shutdown auxillary FuncContainer processes  */
void ObitMultiProcShutdown(ObitErr *err);

/* Set argument for a given job */
void ObitMultiProcSetFuncArg(ObitMultiProc* in, olong jobNo, ObitInfoList *arg);

/* Set execute flag for a given job */
void ObitMultiProcSetExecFlag(ObitMultiProc* in, olong jobNo, gboolean flag);

/* Execute selected jobs */
void ObitMultiProcExecute (ObitMultiProc* in, ofloat timeout, ObitErr *err);

/* Get pointer to return values for a given job */
ObitInfoList* ObitMultiProcGetFuncRet(ObitMultiProc* in, olong jobNo);

/** Public: ClassInfo pointer */
gconstpointer ObitMultiProcGetClass (void);

/** Public: Copy (deep) constructor. */
ObitMultiProc* ObitMultiProcCopy  (ObitMultiProc *in, ObitMultiProc *out, ObitErr *err);

/** Public: Copy structure. */
void ObitMultiProcClone (ObitMultiProc *in, ObitMultiProc *out, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitMultiProcClassDef.h"
} ObitMultiProcClassInfo; 

#endif /* OBITMULTIPROC_H */ 
