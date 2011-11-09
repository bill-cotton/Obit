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
#ifndef OBITWVRCOEF_H 
#define OBITWVRCOEF_H 
#include <glib.h>
#include "ObitTypes.h"
#include "ObitMem.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitWVRCoef.h
 * ObitWVRCoef WVR calibration coefficient class
 * 
 * \section ObitWVRCoefUsage Usage
 * Instances can be obtained using the #newObitWVRCoef constructor or a 
 * pointer duplicated using the #ObitWVRCoefRef function.
 * When an instance is no longer needed, use the #ObitWVRCoefUnref function
 * to release it.
 */

/*---------------Class Structures---------------------------*/
/**  ObitWVRCoef Class structure. */
typedef struct {
  /** class name for verification */
  gchar *className;
  /** Number of entries. */
  olong number;
  /** glib singly linked list */
  GSList* list;
  /** Reference count of pointers to this object. */
  olong ReferenceCount;
} ObitWVRCoef;

/**  ObitWVRCoef Stack Element structure. */
typedef struct {
  /** Timerange (days) */
  ofloat timeRange[2];
  /** Source Id */
  olong sourId;
  /** Antenna number */
  olong ant;
  /** Array of dTdL coefficients */
  ofloat dTdL[4];
  /** Noise weighted dTdL coefficients */
  ofloat c[4];
  /** Weight of solution (c/c_err) */
  ofloat wt;
} ObitWVRCoefElem;


/* Private functions defined only in ObitWVRCoef.c */
/*---------------Public functions---------------------------*/

/** Public: Constructor. */
ObitWVRCoef* newObitWVRCoef (void);

/** Public: Reference to object, update reference count. */
ObitWVRCoef* ObitWVRCoefRef (ObitWVRCoef* in);

/** Public: Unreference object, destroy if no more references. */
ObitWVRCoef* ObitWVRCoefUnref (ObitWVRCoef* in);

/** Public: Add entry. */
void ObitWVRCoefAdd (ObitWVRCoef* in, ofloat timeRange[2],
		     olong sourId, olong ant, ofloat dTdL[4], 
		     ofloat wt);

/** Public: Calibrate a measurement, return excess path */
ofloat ObitWVRCoefCal (ObitWVRCoef* in, ofloat time,
		       olong sourId, olong ant, ofloat T[4]);

/** Public: Print contents to file (e.g. stdout) */
void ObitWVRCoefPrint (ObitWVRCoef* in, FILE *file);

/** Public: Returns TRUE if input is a  ObitWVRCoef* */
gboolean ObitWVRCoefIsA (ObitWVRCoef* in);

#endif /* OBITWVRCOEF_H */ 
