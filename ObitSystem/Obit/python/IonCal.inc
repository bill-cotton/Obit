/* $Id: IonCal.inc,v 1.3 2007/09/07 12:32:42 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Ionospheric calibration utilities      */
/*                                                                    */
/*;  Copyright (C) 2007-2008                                          */
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
#include "ObitIoN2SolNTable.h"
%}

%inline %{
ObitTable* IoN2SolNTableConvert (ObitUV *inUV, long outSNVer, ObitTable *NITable,
                                 double pos[2], ObitErr *err) {
    ObitTableNI *lNITable = NULL;
    ObitTableSN *loutSN   = NULL;
    ofloat lpos[2];
    olong  loutSNVer = (olong)outSNVer;
    oint numPol=1, numIF=1;

    lNITable  = ObitTableNIConvert(NITable);
    if (inUV->myDesc->jlocif>=0) numIF = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
    if (inUV->myDesc->jlocs>=0) numPol = inUV->myDesc->inaxes[inUV->myDesc->jlocs];

    /* Pre existing SN table? */
    if (outSNVer>0) 
        loutSN = newObitTableSNValue ("Ion SN table", (ObitData*)inUV, &loutSNVer,
   			              OBIT_IO_WriteOnly, numPol, numIF, err);

    lpos[0] = (ofloat)pos[0];
    lpos[1] = (ofloat)pos[1];
    loutSN   = ObitIoN2SolNTableConvert (inUV, lNITable, loutSN, lpos, err);
    lNITable = ObitTableNIUnref(lNITable);
    return (ObitTable*)loutSN;
} // end IoN2SolNTableConvert 


%}
