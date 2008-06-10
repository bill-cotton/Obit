/* $Id: ObitZernike.c,v 1.3 2006/10/29 17:28:46 bcotton Exp $                            */
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

#include "ObitZernike.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitZernike.c
 * ObitZernike function definitions.
 */

/*----------------------Public functions---------------------------*/

/**
 * Zernike term N for X and Y 
 * \param N  1-rel term number, between 1 and 18 supported
 * \param X  "X" rectangular coordinate on unit circle
 * \param Y  "Y" rectangular coordinate on unit circle
 * \return Zernike term
 */
ofloat ObitZernike (gint n, ofloat x, ofloat y)
{
  ofloat r, out = 0.0;
  
  if (n == 1) {
    out = 1.0;
  } else if (n == 2) {
    out = x;
  } else if (n == 3) {
    out = y;
  } else if (n == 4) {
    r = y*y + x*x;
    out = 2.0 * r - 1.0;
  } else if (n == 5) {
    out = x*x - y*y;
  } else if (n == 6) {
    out = 2.0*y*x;
  } else if (n == 7) {
    r = y*y + x*x;
    out = x * (3.0 * r - 2.0);
  } else if (n == 8) {
    r = y*y + x*x;
    out = y * (3.0 * r - 2.0);
  } else if (n == 9) {
    r = y*y + x*x;
    out = 6.0*r*r - 6.0*r + 1.0;
  } else if (n == 10) {
    r = y*y + x*x;
    out = 4.0*x*x*x - 3.0*x*r;
  } else if (n == 11) {
    r = y*y + x*x;
    out = 3.0*y*r - 4.0*y*y*y;
  } else if (n == 12) {
    r = y*y + x*x;
    out = (-3.0+4.0*r)*(x*x-y*y);
  } else if (n == 13) {
    r = y*y + x*x;
    out = 2.0*y*x*(-3.0+4.0*r);
  } else if (n == 14) {
    r = y*y + x*x;
    out = x*(3.0-12.0*r+10.0*r*r);
  } else if (n == 15) {
    r = y*y + x*x;
    out = y*(3.0-12.0*r+10.0*r*r);
  } else if (n == 16) {
    r = y*y + x*x;
    out = -1.0+12.0*r-30.0*r*r+20.0*r*r*r;
  } else if (n == 17) {
    r = y*y + x*x;
    out = 8.0*x*x*x*x - 7.0*r*x*x + r*y*y;
  } else if (n == 18) {
    r = y*y + x*x;
    out = 4.0*r*x*y - 8.0*x*y*y*y;
  } 
  
  return out;
} /* end ObitZernike */

/**
 *  Zernike term N gradient in X for X and Y 
 * \param N  1-rel term number, between 1 and 18 supported
 * \param X  "X" rectangular coordinate on unit circle
 * \param Y  "Y" rectangular coordinate on unit circle
 * \return Zernike term
 */
ofloat ObitZernikeGradX (gint n, ofloat x, ofloat y)
{
  ofloat r=0.0, out = 0.0;

  if (n == 1) {
    out = 0.0;
  } else if (n == 2) {
    out = 1.0;
  } else if (n == 3) {
    out = 0.0;
  } else if (n == 4) {
    out = 4.0 * x;
  } else if (n == 5) {
    out = 2.0 * x;
  } else if (n == 6) {
    out = 2.0 * y;
  } else if (n == 7) {
    r = (x*x + y*y);
    out = 6.0*x*x + 3.0*r - 2.0;
  } else if (n == 8) {
    out = 6.0*y*x;
  } else if (n == 9) {
    r = (x*x + y*y);
    out = 24.0*r*x - 12.0*x;
  } else if (n == 10) {
    r = (x*x + y*y);
    out = 6.0*x*x - 3.0*r;
  } else if (n == 11) {
    r = (x*x + y*y);
    out = 6.0*x*y;
  } else if (n == 12) {
    r = (x*x + y*y);
    out = -8.0*x*y*y - 6.0*x + 8.0*x*x*x + 8.0*r*x;
  } else if (n == 13) {
    r = (x*x + y*y);
    out = -6.0*y + 16.0*x*x*y + 8.0*r*y;
  } else if (n == 14) {
    r = (x*x + y*y);
    out = 3.0-12.0*r - 24.0*x*x + 10.0*r*r + 40.0*r*x*x;
  } else if (n == 15) {
    r = (x*x + y*y);
    out = -24.0*x*y + 40.0*r*x*y;
  } else if (n == 16) {
    r = (x*x + y*y);
    out = 24.0*x - 120.0*r*x + 120.0*r*r*x;
  } else if (n == 17) {
    r = (x*x + y*y);
    out = 18.0*x*x*x - 14.0*r*x + 2.0*x*y*y;
  } else if (n == 18) {
    r = (x*x + y*y);
    out = 6.0*r*y + 8.0*x*x*y - 2.0*r*y - 8.0*y*y*y;
  } 

  return out;
} /* end ObitZernikeGradX */

/**
 * Zernike term N gradient in Y for X and Y 
 * \param N  1-rel term number, between 1 and 18 supported
 * \param X  "X" rectangular coordinate on unit circle
 * \param Y  "Y" rectangular coordinate on unit circle
 * \return Zernike term
 */
ofloat ObitZernikeGradY (gint n, ofloat x, ofloat y)
{
  ofloat r, out = 0.0;

  if (n == 1) {
    out = 0.0;
  } else if (n == 2) {
    out = 0.0;
  } else if (n == 3) {
    out = 1.0;
  } else if (n == 4) {
    out = 4.0*y;
  } else if (n == 5) {
    out = -2.0*y;
  } else if (n == 6) {
    out = 2.0*x;
  } else if (n == 7) {
    out = 6.0*x*y;
  } else if (n == 8) {
    r = (x*x + y*y);
    out = 3.0*r + 6.0*y*y - 2.0;
  } else if (n == 9) {
    r = (x*x + y*y);
    out = 24.0*r*y - 12.0*y;
  } else if (n == 10) {
    out = -6.0*x*y;
  } else if (n == 11) {
    r = (x*x + y*y);
    out = 3.0*r - 6.0*y*y;
  } else if (n == 12) {
    r = (x*x + y*y);
    out = 8.0*x*x*y + 6.0*y - 8.0*r*y - 8.0*y*y*y;
  } else if (n == 13) {
    r = (x*x + y*y);
    out = -6.0*x + 8.0*r*x + 16.0*x*y*y;
  } else if (n == 14) {
    r = (x*x + y*y);
    out = -24.0*x*y + 40.0*r*x*y;
  } else if (n == 15) {
    r = (x*x + y*y);
    out = 3.0 - 12.0*r - 24.0*y*y + 10.0*r*r + 40.0*r*y*y;
  } else if (n == 16) {
    r = (x*x + y*y);
    out = 24.0*y - 120.0*r*y + 120.0*r*r*y;
  } else if (n == 17) {
    r = (x*x + y*y);
    out = -14.0*x*x*y + 2.0*y*y*y + 2.0*r*y;
  } else if (n == 18) {
    r = (x*x + y*y);
    out = 4.0*r*x + 8.0*x*y*y - 24.0*x*y*y;
  } 
  
  return out;
  } /* end ObitZernikeGradY */
  
  /**
 * Zernike term N for Polar coordinates rho and phi
 * taken from http://wyant.opt-sci.arizona.edu/zernikes/zernikes.htm 
 * \param N   1-rel term number, between 1 and 36 supported 
 * \param rho radial coordinate on unit circle
 * \param phi azimuthal coordinate on unit circle (radian)
 * \return Zernike term
 */
ofloat ObitZernikePolar (gint n, ofloat rho, ofloat phi)
{
  ofloat out = 0.0;

  if (n == 1) {
    out = 1.0;
  } else if (n == 2) {
    out = rho * cos(phi);
  } else if (n == 3) {
    out = rho * sin(phi);
  } else if (n == 4) {
    out = 2.0*rho*rho - 1.0;
  } else if (n == 5) {
    out = rho*rho*cos (2.0*phi);
  } else if (n == 6) {
    out = rho*rho*sin (2.0*phi);
  } else if (n == 7) {
    out = rho*(3.0*rho*rho - 2.0) * cos (phi);
  } else if (n == 8) {
    out = rho*(3.0*rho*rho - 2.0) * sin (phi);
  } else if (n == 9) {
    out = 1.0 - 6.0 * rho*rho + 6.0*rho*rho*rho*rho;
  } else if (n == 10) {
    out = rho*rho*rho*cos (3.0*phi);
  } else if (n == 11) {
    out = rho*rho*rho*sin (3.0*phi);
  } else if (n == 12) {
    out = rho*rho*(-3.0+4.0*rho*rho)*cos(2.0*phi);
  } else if (n == 13) {
    out = rho*rho*(-3.0+4.0*rho*rho)*sin(2.0*phi);
  } else if (n == 14) {
    out = rho*(3.0-12.0*rho*rho+10.0*(rho*rho*rho*rho))*cos(phi);
  } else if (n == 15) {
    out = rho*(3.0-12.0*rho*rho+10.0*(rho*rho*rho*rho))*sin(phi);
  } else if (n == 16) {
    out = -1.0+12.0*rho*rho-30.0*(rho*rho*rho*rho)+20.0*(rho*rho*rho*rho*rho*rho);
  } else if (n == 17) {
    out = (rho*rho*rho*rho)*cos(4.0*phi);
  } else if (n == 18) {
    out = (rho*rho*rho*rho)*sin(4.0*phi);
  } else if (n == 19) {
    out = (rho*rho*rho)*(-4.0+5.0*rho*rho)*cos(3.0*phi);
  } else if (n == 20) {
    out = (rho*rho*rho)*(-4.0+5.0*rho*rho)*sin(3.0*phi);
  } else if (n == 21) {
    out = rho*rho*(6.0-20.0*rho*rho+15.0*(rho*rho*rho*rho))*cos(2.0*phi);
  } else if (n == 22) {
    out = rho*rho*(6.0-20.0*rho*rho+15.0*(rho*rho*rho*rho))*sin(2.0*phi);
  } else if (n == 23) {
    out = rho*(-4.0+30.0*rho*rho-60.0*(rho*rho*rho*rho)+35.0*(rho*rho*rho*rho*rho*rho))*cos(phi) ;
  } else if (n == 24) {
    out = rho*(-4.0+30.0*rho*rho-60.0*(rho*rho*rho*rho)+35.0*(rho*rho*rho*rho*rho*rho))*sin(phi) ;
  } else if (n == 25) {
    out = 1.0-20.0*rho*rho+90.0*(rho*rho*rho*rho)-140.0*(rho*rho*rho*rho*rho*rho) + 
      70.0*(rho*rho*rho*rho*rho*rho*rho*rho);
  } else if (n == 26) {
    out = (rho*rho*rho*rho*rho)*cos(5.0*phi);
  } else if (n == 27) {
    out = (rho*rho*rho*rho*rho)*sin(5.0*phi);
  } else if (n == 28) {
    out = (rho*rho*rho*rho)*(-5.0+6.0*rho*rho)*cos(4.0*phi);
  } else if (n == 29) {
    out = (rho*rho*rho*rho)*(-5.0+6.0*rho*rho)*sin(4.0*phi);
  } else if (n == 30) {
    out = (rho*rho*rho)*(10.0-30.0*rho*rho+21.0*(rho*rho*rho*rho)) * cos(3.0*phi) ;
  } else if (n == 31) {
    out = (rho*rho*rho)*(10.0-30.0*rho*rho+21.0*(rho*rho*rho*rho)) * sin(3.0*phi) ;
  } else if (n == 32) {
    out = rho*rho*(-10.0+60.0*rho*rho-105.0*(rho*rho*rho*rho) + 56.0*(rho*rho*rho*rho*rho*rho))*cos(2.0*phi);
  } else if (n == 33) {
    out = rho*rho*(-10.0+60.0*rho*rho-105.0*(rho*rho*rho*rho) + 56.0*(rho*rho*rho*rho*rho*rho))*sin(2.0*phi);
  } else if (n == 34) {
    out = rho*(5.0-60.0*rho*rho+210.0*(rho*rho*rho*rho)-280.0*(rho*rho*rho*rho*rho*rho) + 
	       126.0*(rho*rho*rho*rho*rho*rho*rho*rho))*cos(phi);
  } else if (n == 35) {
    out = rho*(5.0-60.0*rho*rho+210.0*(rho*rho*rho*rho)-280.0*(rho*rho*rho*rho*rho*rho) + 
	       126.0*(rho*rho*rho*rho*rho*rho*rho*rho))*sin(phi);
  } else if (n == 36) {
    out = -1.0 + 30*rho*rho - 210.0*(rho*rho*rho*rho) + 560.0*(rho*rho*rho*rho*rho*rho) - 
      630.0*(rho*rho*rho*rho*rho*rho*rho*rho) + 252.0*(rho*rho*rho*rho*rho*rho*rho*rho*rho*rho);
  } 

  return out;
} /* end  ObitZernikePolar  */

