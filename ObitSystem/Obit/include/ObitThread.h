/* $Id$                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2003                                          */
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
#ifndef OBITTHREAD_H 
#define OBITTHREAD_H 
#include <glib.h>

/**
 * \file ObitThread.h
 * ObitThread (for multi threading) class definition.
 */
/*--------------Class definitions-------------------------------------*/
/* basically dummied for now.*/
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
} ObitThread;

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

/** Public: Unlock (mutex) object  */
void ObitThreadUnlock (ObitThread *in);

/** Public: Join thread  */
void ObitThreadJoin  (ObitThread *in);

/** Public: Returns true if input is a  ObitThread* */
gboolean ObitThreadIsA (ObitThread* in);

#endif /* OBITTHREAD_H */ 

