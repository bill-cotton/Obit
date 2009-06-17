/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitFarray type                        */
/*                                                                    */
/*;  Copyright (C) 2005,2009                                          */
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

#include "ObitFArrayUtil.h"
%}


%inline %{
/* Return list with [0]=FWHM, [1] = peak, [2] = cenX, [3]=cenY, 4=RMS residual */
extern PyObject* FArrayUtilFitCGauss (ObitFArray *in, float FWHM, float center[2], 
			      float peak, ObitErr *err)
{
  ofloat lFWHM=FWHM, lcenter[2]={center[0],center[1]}, lpeak=peak;
  ofloat RMS;
  PyObject *o;

  RMS = ObitFArrayUtilFitCGauss (in, &lFWHM, lcenter, &lpeak, err);
  // return list
  o = PyList_New(5);
  PyList_SetItem(o, 0, PyFloat_FromDouble((double)lFWHM));
  PyList_SetItem(o, 1, PyFloat_FromDouble((double)lpeak));
  PyList_SetItem(o, 2, PyFloat_FromDouble((double)lcenter[0]));
  PyList_SetItem(o, 3, PyFloat_FromDouble((double)lcenter[1]));
  PyList_SetItem(o, 4, PyFloat_FromDouble((double)RMS));
  return o;
} // end  FArrayUtilFitCGauss

// Convolution 
ObitFArray* FArrayUtilConvolve (ObitFArray *in1, ObitFArray *in2, 
	 		        ObitErr *err)
{
  return ObitFArrayUtilConvolve (in1, in2, err);
}  // end FArrayUtilConvolve
 
// Gaussian UV Taper 
ObitFArray* FArrayUtilUVGaus (int *naxis, float *cells, float maprot,
			      float Gaumaj, float Gaumin, float GauPA)
{
  olong lnaxis[2];
  ofloat lcells[2];

  lnaxis[0] = (olong)naxis[0];  lnaxis[1] = (olong)naxis[1]; 
  lcells[0] = (ofloat)cells[0]; lcells[1] = (ofloat)cells[1]; 
  return ObitFArrayUtilUVGaus (&lnaxis[0], &lcells[0], (ofloat)maprot, 
                               (ofloat)Gaumaj, (ofloat)Gaumin, (ofloat)GauPA);
}  // end FArrayUtilConvolve
%}


