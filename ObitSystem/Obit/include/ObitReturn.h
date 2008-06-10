/* $Id: ObitReturn.h,v 1.2 2007/02/19 15:07:28 bcotton Exp $  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
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
#ifndef OBITRETURN_H
#define OBITRETURN_H
#include "Obit.h"
#include "ObitFile.h"
#include "ObitErr.h"
#include "ObitInfoList.h"

/**
 * \file ObitReturn.h
 * Output utility dumper
 *
 * This file contains utility functions for writing program output 
 * parameters to a text file.  Parameters are passed in an ObitInfoList.
 * Supported types are:
 * \li integers
 * \li boolean
 * \li floats
 * \li character strings
 *
 * To use do to following:
 * \li Create an ObitInfoList and fill with values
 * \li Use ObitReturnDump to write text file
 *
 * \section output_file output file
 * The text input file is in keyword=value form.
 * An entry consists of a header line beginning with "&Key = ",
 * followed by the name of the entry, a data type code and the dimensionality
 * as a Fortran dimensionality, (n,m) = n x m array with elements on the n axis varying 
 * most rapidly. Values start on the second line of an entry and take as many lines as 
 * necessary with a list of blank separated values.  Character strings are one per line.
 * The type codes are:
 * \li Str String, first element of dimension is number of characters in each string,
 *         Note: these are NOT NULL terminated.
 * \li Int Integer as oint.
 * \li Boo Boolean, 'T' = true, 'F' = false
 * \li Dbl Float
 * \li Flt Double
 * # Example Input file
 * $Key = stringD Str (48,1,1)
 * somefile.fits
 * $Key = integerD Int (1,1,1)
 * 13
 * $Key = floatD Flt (2,1,1)
 * 5.6789 
 * 3.478e5
 * $Key = doubleD Dbl (2,1,1)
 * 5.6789123456
 * $Key = BoolD Boo (1,1,1)
 * T
 */

/*---------------Public functions---------------------------*/
/** Public: Dump return code and ObitInfoList to a text file. */
ObitIOCode ObitReturnDumpRetCode(gint retCode, gchar *outfile, 
				 ObitInfoList *list, ObitErr *err);

/** Public: Dump ObitInfoList to a text file. */
ObitIOCode ObitReturnDump(gchar *outfile, ObitInfoList *list, ObitErr *err);

#endif  /* OBITRETURN_H */
