/* $Id$   */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITUTIL_H 
#define OBITUTIL_H 

#include "Obit.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUtil.h
 * ObitUtil utility module definition.
 *
 */

/*---------------Public functions---------------------------*/
/** Public: Get mean value of an array with magic value blanking */
ofloat meanValue (ofloat *array, olong incs, olong n);

/** Public: Get median value of an array with magic value blanking */
ofloat medianValue (ofloat *array, olong incs, olong n);

/** Public: Get average around median value of an array 
    with magic value blanking */
ofloat medianAvg (ofloat *array, olong incs, olong navg, gboolean doWt, olong n);


/** Public: Determine running median and sigma of a float array */ 
void RunningMedian (glong n, olong wind, ofloat *array, ofloat alpha, 
		    ofloat *RMS, ofloat *out, ofloat *work);

/** Public: Median value of an array */
ofloat MedianLevel (olong n, ofloat *value, ofloat alpha);

/** Public: Median sigma of an array */
ofloat MedianSigma (gint n, ofloat *value, ofloat mean);

/** Public: Fit polynomial with magic value blanking */
void  FitPoly (ofloat *poly, olong order, ofloat *x, ofloat *y, ofloat *wt, 
	       olong n);

/** Public: Evaluate polynomial */
ofloat  EvalPoly (olong order, ofloat *poly, ofloat x);
#endif /* OBITUTIL_H */ 
