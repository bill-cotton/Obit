/* $Id:  $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for parameter file parser utilities        */
/*                                                                    */
/*;  Copyright (C) 2012                                               */
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

#include "ObitSurveyUtil.h"
%}


%inline %{
// Print contents of a VL table
extern void OSurveyVLPrint (ObitTable *VLTable, ObitImage *image, 
  char *streamname, ObitErr *err) {
   FILE *outStream=NULL;
  if (!strncmp(streamname, "stdout", 6)) outStream = stdout;
  if (!strncmp(streamname, "stderr", 6)) outStream = stderr;
  ObitSurveyUtilVLPrint ((ObitTableVL*)VLTable, image, outStream, err);
} // end OSurveyVLPrint

// Print Selected contents of an NVSS VL table
extern int OSurveyNVSSPrint (ObitPrinter *printer, ObitData *data, 
  int VLVer, int first, int last, ObitErr *err) {
  gboolean lfirst = first!=0;
  gboolean llast  = last!=0;
  return ObitSurveyNVSSPrint (printer, data, (olong)VLVer, lfirst, llast, err);
} // end OSurveyNVSSPrint

// Print Selected contents of an VLSS VL table
extern int OSurveyVLSSPrint (ObitPrinter *printer, ObitData *data, 
  int VLVer, int first, int last, ObitErr *err) {
  gboolean lfirst = first!=0;
  gboolean llast  = last!=0;
  return ObitSurveyVLSSPrint (printer, data, (olong)VLVer, lfirst, llast, err);
} // end OSurveyVLSSPrint


%}


