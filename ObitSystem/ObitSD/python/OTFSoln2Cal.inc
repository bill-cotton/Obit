/* $Id: OTFSoln2Cal.inc,v 1.2 2005/01/02 22:10:16 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for OTF calibration Utilities              */
/*                                                                    */
/*;  Copyright (C) 2004                                               */
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
#include "ObitOTFSoln2Cal.h"
%}



%inline %{

extern ObitTable* OTFSoln2Cal (ObitOTF *in, ObitOTF *out, ObitErr *err) {
  ObitTableOTFCal *CalTable;

  CalTable =  ObitOTFSoln2Cal (in, out, err);
  return (ObitTable*)CalTable;
} // end OTFSoln2Cal

%}
