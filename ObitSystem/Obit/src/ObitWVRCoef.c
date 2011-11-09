/* $Id$         */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2011                                               */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include <sys/types.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ObitWVRCoef.h"

/** name of the class defined in this file */
static gchar *myClassName = "ObitWVRCoef";

/**
 * \file ObitWVRCoef.c
 * ObitWVRCoef  WVR calibration coefficient class 
 * Stores calculates sensitivities and uses them to
 * determine the excess path from WVR temperature measurements.
 * Based on software by Bojan Nikolic.
 *
 * \section ObitWVRCoefUsage Usage of member pointers.
 * The Ref and Unref member functions should always be used to make a 
 * copy of an object pointer or to release it.
 * The ref function increments its reference count and returns a pointer.
 * The unref function decrements the reference count, deleted the object
 * if the value is below 1 and returns NULL.
 * Unreferenced pointers should always be NULLed or set to another 
 * valid value.
 */

/*---------------Private function prototypes----------------*/
/** Private: Create ObitWVRCoefElem. */
static ObitWVRCoefElem* 
newObitWVRCoefElem (ofloat timeRange[2], olong sourId, olong ant, 
		    ofloat dTdL[4], ofloat wt);

/** Private: Destroy ObitWVRCoefElem. */
static ObitWVRCoefElem* freeObitWVRCoefElem (ObitWVRCoefElem *in);

/** Private: Destructor. */
static ObitWVRCoef* freeObitWVRCoef (ObitWVRCoef *in);

/** Private: Find item in a list */
static ObitWVRCoefElem* 
ObitWVRCoefFind (ObitWVRCoef *in, ofloat time, olong sourId, olong ant);

/** Time to String */
static void T2String (ofloat time, gchar *msgBuf);
/*---------------Public functions---------------------------*/
/**
 * ObitWVRCoef Constructor.
 * \return the new object.
 */
ObitWVRCoef* newObitWVRCoef (void)
{
  ObitWVRCoef* me;

  /* allocate structure */
  me = ObitMemAlloc0Name(sizeof(ObitWVRCoef),"ObitWVRCoef");

  /* initialize values */
  me->className = g_strdup(myClassName);
  me->number    = 0;
  me->list      = NULL;
  me->ReferenceCount = 1;

  return me;
} /* end newObitWVRCoef */

/**
 * Unconditional ObitWVRCoef destructor.
 * \param in Object to delete
 * \return NULL pointer.
 */
ObitWVRCoef* freeObitWVRCoef (ObitWVRCoef *in)
{
  GSList *tmp;

  /* error checks */
  g_assert (ObitWVRCoefIsA(in));

  /* clear list */
  tmp = in->list;
  while (tmp!=NULL) {
    if (tmp->data) freeObitWVRCoefElem(tmp->data);
    tmp = g_slist_next(tmp);
  }

  /* delete members  */
  g_slist_free(in->list);

  /* delete object */
  ObitMemFree (in);

  return NULL;
} /* end freeObitWVRCoef */


/**
 * To reference pointer, incrementing ReferenceCount and returning 
 * the pointer.
 * This function should always be used to copy pointers as this 
 * will ensure a proper reference count.
 * \param in Pointer to object to link.
 * \return the pointer to in.
 */
ObitWVRCoef* ObitWVRCoefRef (ObitWVRCoef* in)
{
  /* error checks */
  g_assert (ObitWVRCoefIsA(in));

  /* increment reference count */
  in->ReferenceCount++;

  return in;
} /* end ObitWVRCoefRef */

/**
 * Always use this function to dismiss an object as it will
 * ensure that the structure is only deleted when there are no 
 * more  pointers to it.
 * \param  in Pointer to structure to unlink.
 * \return NULL pointer.
 */
ObitWVRCoef* ObitWVRCoefUnref (ObitWVRCoef* in)
{
  if (in==NULL) return NULL;
  /* error checks */
  g_assert (ObitWVRCoefIsA(in));

  /* decrement reference count, delete if non positive */
  in->ReferenceCount--;
  if (in->ReferenceCount<=0) freeObitWVRCoef(in);

  return NULL;
} /* end ObitWVRCoefUnref */

/**
 * Add an item to the list
 * \param in        Pointer to object to add to.
 * \param timeRange Beginning and end times in days
 * \param sourId    Source ID
 * \param ant       Antenna number (1-rel)
 * \param dTdL      Sensitivity coefficients
 * \param wt        Weight for entry
 */
void ObitWVRCoefAdd (ObitWVRCoef *in, ofloat timeRange[2],
		     olong sourId, olong ant, ofloat dTdL[4], 
		     ofloat wt)
{
  ObitWVRCoefElem* elem;

  /* error checks */
  g_assert (ObitWVRCoefIsA(in));

  /* create ObitWVRCoefElem with message */
  elem = newObitWVRCoefElem(timeRange, sourId, ant, dTdL, wt);

  /* add to stack */
  in->list = g_slist_append (in->list, elem); /* add to list */
  in->number++;   /* add one to the count */
} /* end ObitWVRCoefAdd */

/**
 * Determine excess path length from a measurement
 * Uses sensitivities previously added by ObitWVRCoefAdd
 * \param in        Pointer to structure with sensitivities
 * \param time      time in days
 * \param sourId    Source ID, -1=>any
 * \param ant       Antenna number (1-rel)
 * \param T         WVR temps
 * \return excess path length(mm), -1 on failure
 */
ofloat ObitWVRCoefCal (ObitWVRCoef *in, ofloat time,
		       olong sourId, olong ant, ofloat T[4])
{
  ofloat excess = -1.0;
  olong i;
  ObitWVRCoefElem* elem;

  /* error checks */
  g_assert (ObitWVRCoefIsA(in));

  /* Find corresponding sensitivity */
  elem = ObitWVRCoefFind (in, time, sourId, ant);
  if (elem==NULL) return excess;

  /* Average path excess */
  excess = 0.0;
  for (i=0; i<4; i++) {
    excess += T[i] * elem->c[i];
  }

  return excess;
} /* end ObitWVRCoefCal */

/**
 * Print the contents of an InfoList to file
 * \param file  Where to write output
 */
void ObitWVRCoefPrint (ObitWVRCoef *in, FILE *file)
{
  GSList *tmp;
  ObitWVRCoefElem* me;
  gchar line[257], t1[16], t2[16];

  /* error checks */
  g_assert (ObitWVRCoefIsA(in));
  g_assert (file != NULL);

  fprintf (file, "Listing Of ObitWVRCoef\n\n");

  /* loop through list */
  tmp = in->list;
  while (tmp!=NULL) {
    if (tmp->data) {
      me = (ObitWVRCoefElem*)tmp->data;
      T2String (me->timeRange[0], t1);
      T2String (me->timeRange[1], t2);
      g_snprintf (line, 256,
		  "%s-%s sou=%3d ant=%3d dTdL=%6.2f %6.2f %6.2f %6.2f c=%7.4f %7.4f %7.4f %7.4f wt %6.2f",
		  t1, t2, me->sourId, me->ant,
		  me->dTdL[0], me->dTdL[1], me->dTdL[2], me->dTdL[3],
		  me->c[0], me->c[1], me->c[2], me->c[3], me->wt);
      fprintf (file, "%s\n", line);
    }
    tmp = g_slist_next(tmp);
  }
  fprintf (file, "\n");
} /* end ObitWVRCoefPrint  */

/**
 * Determines if the input object is a member of this class
 * \param in Pointer to object to test.
 * \return TRUE if member else FALSE.
 */
gboolean ObitWVRCoefIsA (ObitWVRCoef* in)
{
  gboolean out;

  /* error checks */
  if (in == NULL) return FALSE;
  if (in->className == NULL) return FALSE;

  /* compare class name member */
  out = !strcmp(in->className, myClassName);

  return out;
} /* end ObitWVRCoefIsA */

/*---------------Private functions-------------------------*/
/**
 * ObitWVRCoefElem Constructor.
 * Calculates noise weighted sensitivities.
 * \param timeRange Beginning and end times in days
 * \param sourId    Source ID
 * \param ant       Antenna number (1-rel)
 * \param dTdL      Sensitivity coefficients
 * \return the new object.
 */
static ObitWVRCoefElem*  
newObitWVRCoefElem (ofloat timeRange[2], olong sourId, olong ant, 
		    ofloat dTdL[4], ofloat wt)
{
  olong i;
  odouble noise[4] = { 0.1, 0.08, 0.08, 0.09};
  odouble dLdT, w, sum;
  ObitWVRCoefElem* me;

  /* allocate structure */
  me = ObitMemAlloc0(sizeof(ObitWVRCoefElem));

  /* set values */
  me->timeRange[0] = timeRange[0];
  me->timeRange[1] = timeRange[1];
  me->sourId       = sourId;
  me->ant          = ant;
  me->dTdL[0]      = dTdL[0];
  me->dTdL[1]      = dTdL[1];
  me->dTdL[2]      = dTdL[2];
  me->dTdL[3]      = dTdL[3];
  me->wt           = wt;

  /* Estimate thermal noise */
  sum = 0.0;
  for (i=0; i<4; i++) {
    if (dTdL[i]!=0.0) {
      dLdT      = 1.0 / dTdL[i];
      w         = pow (noise[i] * dLdT, -2);
      me->c[i]  = w*dLdT;
      sum      += w;
    } else { /* bad sensitivity */
      me->c[i]  = 0.0;
    }
  }
  /* Normalize weighted sensitivities */
  for (i=0; i<4; i++) {
    me->c[i]   /= sum;
  }

  return me;
} /* end newObitWVRCoefElem */

/**
 * ObitWVRCoefElem destructor.
 * \param in Object to delete
 * \return NULL.
 */
ObitWVRCoefElem* freeObitWVRCoefElem (ObitWVRCoefElem *in)
{
  g_assert (in != NULL);

  /* deallocate structure */
  ObitMemFree (in);
  
  return NULL;
} /* end freeObitWVRCoefElem */

/** 
 * Find the element with a given time, source and antenna in an ObitWVRCoef
 * \param in     Pointer to WVRCoef.
 * \param time   Time in days, the last
 *               or including valid entry at a time preceeding time
 *               and matching the other criteria is returned.
 * \param sourId source Id, -1 => any
 * \param ant    Antenna number
 * \return pointer to ObitWVRCoefElem or NULL if not found.
 */
static ObitWVRCoefElem* 
ObitWVRCoefFind (ObitWVRCoef *in, ofloat time, olong sourId, olong ant)
{
  GSList *tmp;
  ObitWVRCoefElem *elem;

  /* loop through list testing elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (ObitWVRCoefElem*)tmp->data;
    /* check if this is a match */
    if ((elem->ant==ant) &&
	((elem->sourId==sourId) || (sourId<0)) &&
	(((time>=elem->timeRange[0])&&(time<=elem->timeRange[1])) || 
	 (time>elem->timeRange[1]))) return elem;
    tmp = g_slist_next(tmp);
  }

  /* relax time constraint and try again */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (ObitWVRCoefElem*)tmp->data;
    /* check if this is a match */
    if ((elem->ant==ant) &&
	((elem->sourId==sourId) || (sourId<0))) 
      return elem;
    tmp = g_slist_next(tmp);
  }

  return NULL; /* didn't find any */
} /* end ObitWVRCoefFind */
  
/**
 * Convert a time as time in days to a printable string
 * \param    time  Beginning time, end time in days
 * \msgBuff  Human readable string as "dd/hh:mm:ss.s"
 *           must be allocated at least 13 characters
 */
static void T2String (ofloat time, gchar *msgBuf)
{
  ofloat rtemp, rt1;
  olong   id1, it1, it2;

  id1 = time;
  rtemp = 24.0 * (time - id1);
  it1 = rtemp;
  it1 = MIN (23, it1);
  rtemp = (rtemp - it1)*60.0;
  it2 = rtemp;
  it2 = MIN (59, it2);
  rt1 = (rtemp - it2)*60.0;
  g_snprintf (msgBuf, 30, "%2.2d/%2.2d:%2.2d:%5.2f",
	      id1, it1, it2, rt1);
} /* end of routine T2String   */ 
