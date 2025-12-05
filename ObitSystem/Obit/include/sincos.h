/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2025                                               */
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
/**
 * Work around the sincos* glibc extension on other platforms.
 *
 * It might be sufficient to replace sincos with sequential calls
 * to sin and cos, since both gcc and clang with -O2 optimization
 * replace them with a single call to sincos anyway.
 *
 * See:
 * - https://stackoverflow.com/a/61451065
 * - https://discuss.python.org/t/sincos-x-from-math-h-missing/22614/4
 * Hope the compilers can get this right 
*/

#ifndef SINCOS_H
#define SINCOS_H

#include <math.h>

//DAMNinline void sincos2(double x, double* p_sin, double* p_cos) {
static inline void sincos2(double x, double* p_sin, double* p_cos) {
  *p_sin = sin(x);
  *p_cos = cos(x);
}
//damninline void sincos2f(float x, float* p_sinf, float* p_cosf) {
static inline void sincos2f(float x, float* p_sinf, float* p_cosf) {
  *p_sinf = sinf(x);
  *p_cosf = cosf(x);
}

#endif /* SINCOS_H */
