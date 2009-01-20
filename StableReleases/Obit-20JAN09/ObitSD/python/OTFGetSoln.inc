/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig module description for OTF calibration Utilities              */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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

%{
#include "ObitOTFGetSoln.h"
%}



%inline %{

extern ObitTable* OTFGetSolnCal (ObitOTF *in, ObitOTF *out, ObitErr *err) {
  ObitTableOTFSoln *SolnTable;

  SolnTable = ObitOTFGetSolnCal (in, out, err);
  return (ObitTable*)SolnTable;
} // end OTFGetSolnCal

extern ObitTable* OTFGetSolnGain (ObitOTF *in, ObitOTF *out, ObitErr *err){
  ObitTableOTFSoln *SolnTable;

  SolnTable = ObitOTFGetSolnGain (in, out, err);
  return (ObitTable*)SolnTable;
} // end  OTFGetSolnGain

extern ObitTable* OTFGetSolnFilter (ObitOTF *in, ObitOTF *out, ObitErr *err) {
  ObitTableOTFSoln *SolnTable;

  SolnTable = ObitOTFGetSolnFilter (in, out, err);
  return (ObitTable*)SolnTable;
} // end OTFGetSolnFilter

extern ObitTable* OTFGetSolnPolyBL (ObitOTF *in, ObitOTF *out, ObitErr *err) {
  ObitTableOTFSoln *SolnTable;

  SolnTable = ObitOTFGetSolnPolyBL (in, out, err);
  return (ObitTable*)SolnTable;
} // end OTFGetSolnPolyBL

extern ObitTable* OTFGetSolnMBBase (ObitOTF *in, ObitOTF *out, ObitErr *err) {
  ObitTableOTFSoln *SolnTable;

  SolnTable = ObitOTFGetSolnMBBase (in, out, err);
  return (ObitTable*)SolnTable;
} // end OTFGetSolnMBBase

extern ObitTable* OTFGetInstCal (ObitOTF *in, ObitOTF *out, ObitErr *err) {
  ObitTableOTFSoln *SolnTable;

  SolnTable = ObitOTFGetInstCal (in, out, err);
  return (ObitTable*)SolnTable;
} // end OTFGetInstCal

extern ObitTable* OTFGetSolnPARGain (ObitOTF *in, ObitOTF *out, ObitErr *err) {
  ObitTableOTFSoln *SolnTable;

  SolnTable = ObitOTFGetSolnPARGain (in, out, err);
  return (ObitTable*)SolnTable;
} // end OTFGetSolnPARGain

extern ObitTable* OTFGetSolnPointTab (ObitOTF *in, ObitOTF *out, ObitErr *err) {
  ObitTableOTFSoln *SolnTable;

  SolnTable = ObitOTFGetSolnPointTab (in, out, err);
  return (ObitTable*)SolnTable;
} // end OTFGetSolnPointTab

extern ObitTable* OTFGetDummyCal (ObitOTF *in, ObitOTF *out, long ver, 
				 int ncoef, ObitErr *err) {
  oint lncoef = ncoef;
  ObitTableOTFCal *CalTable;

  CalTable = ObitOTFGetDummyCal (in, out, ver, lncoef, err);
  return (ObitTable*)CalTable;
} // end  OTFGetDummyCal

extern void OTFGetSolnFlag (ObitOTF *in, ObitImage *model, ObitOTF *out, 
                            long FGver, ObitErr *err){
  ObitOTFGetSolnFlag (in, model, out, (olong)FGver, err);
} // end  OTFGetSolnFlag

extern void OTFGetSolnFlagNoModel (ObitOTF *in, ObitOTF *out, long FGver, ObitErr *err){
  ObitOTFGetSolnFlag (in, NULL, out, (olong)FGver, err);
} // end  OTFGetSolnFlagNoModel

%}
