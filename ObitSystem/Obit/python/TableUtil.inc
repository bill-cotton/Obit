/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for Obit Table Utilities                   */
/*                                                                    */
/*;  Copyright (C) 2006                                               */
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

%{

#include "ObitTableCCUtil.h"
%}

%inline %{
olong TableCCUtilMerge (ObitTable *in, ObitTable *out, ObitErr *err)
{
  ObitIOCode ret=1;
  gchar *routine="TableCCUtilMerge";

  /* Check that input TableCCs */
  Obit_retval_if_fail((ObitTableCCIsA(in)), err, ret,
		      "%s: input %s NOT an AIPS CC table", routine, in->name);
  Obit_retval_if_fail((ObitTableCCIsA(out)), err, ret,
		      "%s: output %s NOT an AIPS CC table", routine, out->name);

  ret = ObitTableCCUtilMerge ((ObitTableCC*)in, (ObitTableCC*)out, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end TableCCUtilMerge
%}


