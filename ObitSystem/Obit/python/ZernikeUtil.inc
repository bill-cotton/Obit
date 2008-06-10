/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for Zernike polynimoal utilities           */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{

#include "ObitZernike.h"
%}


%inline %{
/* Return Zernike term n for X and Y */
extern float Zernike (int n, float x, float y) {
  olong ln = (olong)n;
  return (float)ObitZernike (ln, (ofloat)x, (ofloat)y);
} // end Zernike

/* Return Zernike term n gradient in X for X and Y */
extern float ZernikeGradX (int n, float x, float y) {
  olong ln = (olong)n;
  return (float)ObitZernikeGradX (ln, (ofloat)x, (ofloat)y);
} // end ZernikeGradX

/* Return Zernike term n gradient in Y for X and Y */
extern float ZernikeGradY (int n, float x, float y) {
  olong ln = (olong)n;
  return (float)ObitZernikeGradY (ln, (ofloat)x, (ofloat)y);
} // end ZernikeGradY

/* Return Zernike term n for Polar coordinates rho and phi */
extern float ZernikePolar (int n, float rho, float phi) {
  olong ln = (olong)n;
  return (float)ObitZernikePolar (ln, (ofloat)rho, (ofloat)phi);
} // end ZernikePolar


%}


