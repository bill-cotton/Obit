/* $Id: ObitParser.h,v 1.2 2007/08/31 17:24:48 bcotton Exp $                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003                                               */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITPARSER_H
#define OBITPARSER_H
#include "Obit.h"
#include "ObitFile.h"
#include "ObitErr.h"
#include "ObitInfoList.h"

/**
 * \file ObitParser.h
 * Input utility parser
 *
 * This file contains utility functions for parsing program input 
 * parameters from an input text file and storing in an ObitInfoList.
 * Order of the entries in the text file is arbitrary.
 * A set of defined value names with default values are given to the parser.
 * Supported types are:
 * \li integers (as oint)
 * \li boolean
 * \li floats  (as gdouble)
 * \li character strings
 *
 * To use do to following:
 * \li Create an ObitInfoList and fill with default values
 * \li Use ObitParserParse to read text file and store in ObitInfoList
 * \li Obtain values from ObitInfoList
 *
 * \section input_file input file
 * The text input file is a free form with keyword=value form.
 * Comments are preceeded by # and an entire line may be a comment.
 * An entry consists of a header line beginning with "$Key = ",
 * followed by the name of the entry, a data type code and the dimensionality
 * as a Fortran dimensionality, (n,m) = n x m array with elementson the n axis varying 
 * most rapidly. Values start on the second line of an entry and take as many lines as 
 * necessary with a list of blank separated values.  Character strings are one per line.
 * The type codes are:
 * \li Str String, first element of dimension is number of characters in each string,
 *         Note: these are NOT NULL terminated.
 * \li Int Integer as oint.
 * \li Boo Boolean, 'T' = true, 'F' = false
 * \li Flt Floating as gdouble,
 * Example Input file
 * \verbatim
  $Key = stringD Str (48)
  somefile.fits  # string input
  $Key = integerD Int (1)
  13
  $Key = floatD Flt (2)
  5.6789 3.478e5
  $Key = BoolD Boo (1)
  T
 \endverbatim
 */

/*---------------Public functions---------------------------*/
/** Public: Parse text file to ObitInfoList. */
ObitIOCode ObitParserParse(gchar *infile, ObitInfoList *list, ObitErr *err);

#endif  /* OBITPARSER_H */
