/* $Id:  $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV Doppler correction utilities        */
/*                                                                    */
/*;  Copyright (C) 2012                                               */
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
#include "ObitDoppler.h"
#include "ObitUVDesc.h"
%}


%inline %{
extern ObitDoppler* newDoppler (char* name) {
  return newObitDoppler (name);
} // end  newDoppler

extern ObitDoppler* DopplerCopy  (ObitDoppler *in, ObitDoppler *out, 
				    ObitErr *err) {
  return ObitDopplerCopy (in, out, err);
} // end  DopplerCopy

extern ObitDoppler* DopplerUnref (ObitDoppler* in) {
  if (!ObitDopplerIsA(in)) return NULL;
  return ObitDopplerUnref(in);
}

extern ObitDoppler*  DopplerRef (ObitDoppler* in) {
  return ObitDopplerRef(in);
}

extern ObitInfoList* DopplerGetList (ObitDoppler* in) {
  return ObitInfoListRef(in->info);
}

extern int DopplerIsA (ObitDoppler* in) {
  return ObitDopplerIsA(in);
}

extern float DopplerFreqLSR (double rest, float vlsr, 
			    double ra, double dec, 
			    char *date, double ut, 
			    double x, double y, double z) {
  odouble out, JD;
  olong year, doy;
  ofloat temp;
  // munge date
  ObitUVDescDate2JD ((gchar*)date, &JD);
  ObitDopplerJD2Date (JD, &year, &doy, &temp);
  out = ObitDopplerFreqLSR (rest, vlsr, ra, dec, year, doy, ut,
                            x, y, z);
  return (float)out;
}

extern ObitUV* DopplerCVel (ObitUV *inUV, int scratch, ObitUV *outUV, 
		   	    ObitErr *err) {
 gboolean lscratch = scratch!=0;
  return ObitDopplerCVel (inUV, lscratch, outUV, err);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the st<ruct for c */
%{
typedef struct {
  ObitDoppler *me;
} Doppler;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDoppler *me;
} Doppler;

%addmethods Doppler { 
  Doppler(char* name) {
     Doppler *out;
     out = (Doppler *) malloc(sizeof(Doppler));
     if (strcmp(name, "None")) out->me = newDoppler(name);
     else out->me = NULL;
     return out;
   }
  ~Doppler() {
   if (self->me->ReferenceCount>0) 
      self->me = DopplerUnref(self->me);
   free(self);
  }
};

