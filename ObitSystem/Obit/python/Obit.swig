/* $Id: Obit.swig,v 1.1.1.1 2004/07/19 16:42:44 bcotton Exp $                            */  
/*--------------------------------------------------------------------*/
/* Swig module description for Obit typeS                             */
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
/* 
Notes:
   python has effectively a different name space for each module dynamically
loaded.  Since Obit uses global class structures extensively it must
be load into a single python module.  This is done by putting the different
Obit class information and the swig typemaps into separate files and appending 
them onto the end of this file to feed to swig.
   The ObitTypeMaps.swig file contains the generic typemaps needed to properly
interface python to Obit.  Obit class specific type maps should have unique 
signatures and be at the beginning of the Obit class (*.inc) files.  
Swig is pretty stupid and the typemaps have to be defined before they are needed.
*/

%module Obit

// Generic Obit object functions
%{
#include "Obit.h"
%}

%inline %{
extern int IsAnObitObject (Obit* in) {
  gconstpointer ObitClassInfo = ObitGetClass();
  return ObitIsA(in, ObitClassInfo);
}

extern char* ObitGetName (Obit* in) {
  gconstpointer ObitClassInfo = ObitGetClass();
  if (ObitIsA(in, ObitClassInfo)) {
    return in->name;
  } else {
    return NULL;
  }
}
%}
