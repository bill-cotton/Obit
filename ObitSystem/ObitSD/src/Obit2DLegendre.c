/* $Id$                            */
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

#include <math.h>
#include "Obit2DLegendre.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file Obit2DLegendre.c
 * Obit2DLegendre class function definitions.
 */

/*----------------------Public functions---------------------------*/
/**
 * Return requested Legendre Polynomial P(n,m)
 * Up to 5th order is supported.
 * Logical groups are:
 * \li 1; zeroth order
 * \li 3: first order
 * \li 6: second order
 * \li 10: third order
 * \li 15: fourth order
 * \li 21: fifth order
 * 
 * \param n  Coefficient order in "X" dimension
 *           Values 0->5 are allowed.
 * \param x  Location on "X" dimension
 * \param y  Location on "Y" dimension
 */
float Obit2DLegendre (gint n, ofloat x, ofloat y)
{
  ofloat out = 0.0;
  
  /* branch on n */
  switch (n) { 

    /* Zeroth order */
  case 0:         /* P(0,0) */
    out = 1.0;
    break;
    
    /* First order */
  case 1:         /* P(1,0) */
    out = x;
    break;
  case 2:         /* P(1,1) */
    out = y;
    break;
    
    /* second order */
  case 3:         /* P(2,0) */
    out = x*x - 0.33333333;
    break;
    
  case 4:         /* P(2,1) */
    out = x*y;
    break;
    
  case 5:         /* P(2,2) */
    out = y*y - 0.33333333;
    break;
    
    /* Third order */
  case 6:         /* P(3,0) */
    out = x*x*x - x * (3.0 / 5.0);
    break;
    
  case 7:         /* P(3,1) */
    out = x*x*y - y * 0.33333333;
    break;
    
  case 8:         /* P(3,2) */
    out = x*y*y - x * 0.33333333;
    break;
    
  case 9:         /* P(3,3) */
    out = y*y*y - y * (3.0 / 5.0);
    break;
    
    /* Fourth order */
  case 10:         /* P(4,0) */
    out = x*x*x*x - x*x * (6.0 / 7.0) + 3.0 / 35.0;
    break;
    
  case 11:         /* P(4,1) */
    out = x*x*x*y - x*y * (3.0 / 5.0);
    break;
    
  case 12:         /* P(4,2) */
    out = x*x*y*y - (x*x + y*y) * 0.33333333 + 1.0 / 9.0;
    break;
    
  case 13:         /* P(4,3) */
    out = x*y*y*y - x*y * (3.0 / 5.0);
    break;
    
  case 14:         /* P(4,4) */
    out = y*y*y*y - y*y * (6.0 / 7.0) + 3.0 / 35.0;
    break;
    
    /* Fifth order */
  case 15:         /* P(5,0) */
    out = x*x*x*x*x - x*x*x * (10.0 / 9.0) + x * (5.0 / 21.0);
    break;
    
  case 16:         /* P(5,1) */
    out = x*x*x*x*y - x*x*y * (6.0 / 7.0) + y * (3.0 / 35.0);
    break;
    
  case 17:         /* P(5,2) */
    out = x*x*x*y*y - x*x*x * 0.33333333 - x*y*y * (3.0/5.0) + x * (1.0/5.0);
    break;
    
  case 18:         /* P(5,3) */
    out = x*x*y*y*y - x*x*y * (3.0/5.0) - y*y*y * 0.33333333 + y * (1.0/5.0);
    break;
    
  case 19:         /* P(5,4) */
    out = x*y*y*y*y - x*y*y * (6.0 / 7.0) + x * (3.0 / 35.0);
    break;
    
  case 20:         /* P(5,5) */
    out = y*y*y*y*y - y*y*y * (10.0 / 9.0) + y * (5.0 / 21.0);
    break;
    
  default:
    g_assert_not_reached(); /* not allowed, barf */
  };

  return out;
} /* end Obit2DLegendre */
