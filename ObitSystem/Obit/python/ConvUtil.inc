/* $Id$  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Convolution utilities                  */
/*                                                                    */
/*;  Copyright (C) 2006,2009                                          */
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
#include "ObitConvUtil.h"
%}


%inline %{
extern void ConvUtilConv (ObitImage *inImage, ObitFArray *convFn, 
		   int doDivide, float rescale,
		   ObitImage *outImage, ObitErr *err) {
  gboolean ldoDivide;
  ofloat lrescale = (ofloat)rescale;
  ldoDivide = doDivide != 0;
  ObitConvUtilConv (inImage, convFn, ldoDivide, lrescale, outImage, err);
} // end ConvUtilConv 

extern void ConvUtilConvGauss (ObitImage *inImage, float maj, float min, float pa, float rescale,
		               ObitImage *outImage, ObitErr *err) {

  ObitConvUtilConvGauss (inImage, (ofloat)maj, (ofloat)min, (ofloat)pa, (ofloat)rescale, 
                         outImage, err);
} // end ConvUtilConvGauss

extern ObitFArray* ConvUtilGaus (ObitImage *inImage, float Maj, float Min, float PA) {
  ofloat Beam[3];
  Beam[0] = Maj;
  Beam[1] = Min;
  Beam[2] = PA;
  return ObitConvUtilGaus (inImage, Beam);
} // end ConvUtilGaus

extern PyObject* ConvUtilDeconv (float fMaj,  float fMin,  float fPA, 
                            float cMaj,  float cMin,  float cPA) {
  ofloat rMaj, rMin, rPA;
  PyObject *outList=NULL, *o=NULL;

  ObitConvUtilDeconv (fMaj, fMin, fPA, cMaj, cMin, cPA, &rMaj, &rMin, &rPA);

  // Package results into output list
  outList = PyList_New(3); 
  o = PyFloat_FromDouble((double)rMaj);
  PyList_SetItem(outList, 0, o);
  o = PyFloat_FromDouble((double)rMin);
  PyList_SetItem(outList, 1, o);
  o = PyFloat_FromDouble((double)rPA);
  PyList_SetItem(outList, 2, o);

  return outList;
} // end ConvUtilGaus

%}
