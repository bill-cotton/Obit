/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for parameter file parser utilities        */
/*                                                                    */
/*;  Copyright (C) 2007                                               */
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

%{

#include "ObitParser.h"
#include "ObitReturn.h"
%}


%inline %{
/* Parse text file, return 0 = OK, else failed */
extern int Parse (char *infile, ObitInfoList *list, ObitErr *err) {
  ObitIOCode retCode;

  retCode = ObitParserParse ((gchar*)infile, list, err);
  if (retCode==OBIT_IO_OK) return 0;
  else return 1;
} // end Parser

/* Dump InfoList to text file, return 0 = OK, else failed */
extern int Dump (char *outfile, ObitInfoList *list, ObitErr *err) {
  ObitIOCode retCode;

  retCode = ObitReturnDump ((gchar*)outfile, list, err);
  if (retCode==OBIT_IO_OK) return 0;
  else return 1;
} // end Parser

%}


