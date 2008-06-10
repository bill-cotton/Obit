/* $Id: ObitThread.c,v 1.1.1.1 2004/07/19 16:42:39 bcotton Exp $                            */
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
#include <string.h>
#include "ObitThread.h"

/** name of the class defined in this file */
static gchar *myClassName = "ObitThread";

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

  /* allocate structure */
  me = g_malloc0(sizeof(ObitThread));

  /* initialize */
  strncpy (me->className, myClassName, 11);
  me->ReferenceCount = 1;

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

  /* DEBUG 
  return NULL;*/

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
 * Lock the object so no other thread can modify it.
 * \param in Pointer to object to be locked.
 */
void ObitThreadLock (ObitThread *in)
{
  /* error checks */
  g_assert (ObitThreadIsA(in));

  /* stubbed
  g_error ("ObitThreadLock is not defined - write me"); */

} /* end ObitThreadLock */

/**
 * Unlock the object so other threads can modify it.
 * \param in Pointer to object to be unlocked.
 */
void ObitThreadUnlock (ObitThread *in)
{
  /* error checks */
  g_assert (ObitThreadIsA(in));

  /* stubbed
  g_error ("ObitThreadUnlock is not defined - write me"); */

} /* end ObitThreadUnlock */

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


