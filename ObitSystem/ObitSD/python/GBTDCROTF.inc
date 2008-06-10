/* $Id: GBTDCROTF.inc,v 1.3 2005/02/06 02:04:17 bcotton Exp $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for GBTDCROTF type                         */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{
#include "ObitGBTDCROTF.h"
%}

%inline %{
ObitGBTDCROTF* GBTDCROTFCreate (char *name) {
  return newObitGBTDCROTF (name);
} // end GBTDCROTFCreate

ObitGBTDCROTF* GBTDCROTFValue (char *name, ObitOTF* outOTF, ObitErr *err) {
  return newObitGBTDCROTFValue (name, outOTF, err);
} // end GBTDCROTFValue

void GBTDCROTFConvert (ObitGBTDCROTF* inGDO, int inDisk, char *scanName, ObitErr *err) {
  ObitGBTDCROTFConvert  (inGDO, inDisk, scanName, err);
} // end GBTDCROTFConvert

ObitGBTDCROTF* GBTDCROTFRef (ObitGBTDCROTF* in) {
  return ObitGBTDCROTFRef (in);
} // end GBTDCROTFRef

ObitGBTDCROTF* GBTDCROTFUnref (ObitGBTDCROTF* in) {
  if (!ObitGBTDCROTFIsA(in)) return NULL;
  return ObitGBTDCROTFUnref (in);
} // end GBTDCROTFUnref

extern int GBTDCROTFIsA (ObitGBTDCROTF* in) {
  return ObitGBTDCROTFIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitGBTDCROTF *me;
} GBTDCROTF;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitGBTDCROTF *me;
} GBTDCROTF;

%addmethods GBTDCROTF { 
  GBTDCROTF(char *name) {
     GBTDCROTF *out;
     out = (GBTDCROTF *) malloc(sizeof(GBTDCROTF));
     if (strcmp(name, "None")) out->me = GBTDCROTFCreate(name);
     else  out->me = NULL;
     return out;
   }
  ~GBTDCROTF() {
    self->me = GBTDCROTFUnref(self->me);
    free(self);
  }
};

