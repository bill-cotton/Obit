/* $Id: ObitSystem.inc,v 1.9 2007/07/26 14:28:42 bcotton Exp $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitSystem type                        */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{
#include <stdio.h>
#include "ObitSystem.h"
#include "ObitMem.h"
%}

%inline %{

extern ObitSystem* Startup (char *pgmName, int pgmNumber, int AIPSuser,
		            int numberAIPSdisk, char* AIPSdir[], 
		            int numberFITSdisk, char* FITSdir[], 
		            int F_TRUE, int F_FALSE, ObitErr *err) {
  oint lF_TRUE  = F_TRUE;
  oint lF_FALSE = F_FALSE;

  return ObitSystemStartup (pgmName, pgmNumber, AIPSuser,
        	            numberAIPSdisk, AIPSdir,
                	    numberFITSdisk, FITSdir,
                            lF_TRUE, lF_FALSE, err);

} //  end Startup

extern ObitSystem* Shutdown (ObitSystem* in) {
  return ObitSystemShutdown(in);
} // end Shutdown

extern int SystemIsInit (void) {
  gboolean retval;
  olong out;
  retval = ObitSystemIsInit ();
  if (retval) out = 1;
  else out = 0;
  return out;
} // end SystemGetPgmName

extern char* SystemToday (void) {
  return ObitToday ();
} // end SystemToday

extern char* SystemGetPgmName (void) {
  return ObitSystemGetPgmName ();
} // end SystemGetPgmName

extern void SystemSetPgmName (char* pgmName) {
  ObitSystemSetPgmName (pgmName);
} // end SystemSetPgmName

extern int SystemGetPgmNumber (void) {
  return ObitSystemGetPgmNumber ();
} // end SystemGetPgmNumber

extern void SystemSetPgmNumber (int pgmNumber) {
  ObitSystemSetPgmNumber (pgmNumber);
} // end SystemSetPgmNumber

extern int SystemGetAIPSuser (void) {
  return ObitSystemGetAIPSuser ();
} // end SystemGetAIPSuser

extern void SystemSetAIPSuser (int user) {
  ObitSystemSetAIPSuser (user);
} // end SystemGetAIPSuser

extern void MemPrint (void) {
  ObitMemPrint(stdout);
} // end MemPrint

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitSystem *me;
} OSystem;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitSystem *me;
} OSystem;

%addmethods OSystem { 
  OSystem(char *pgmName, int pgmNumber, int AIPSuser,
	     int numberAIPSdisk, char* AIPSdir[], 
	     int numberFITSdisk, char* FITSdir[], 
	     int F_TRUE, int F_FALSE, ObitErr *err) {
     OSystem *out;
     gchar **Adir, **Fdir;

     if (numberAIPSdisk>0) Adir = AIPSdir;
     else Adir = NULL;
     if (numberFITSdisk>0) Fdir = FITSdir;
     else Fdir = NULL;
     out = (OSystem *) malloc(sizeof(OSystem));
     out->me = ObitSystemStartup(pgmName, pgmNumber, AIPSuser,
                                 numberAIPSdisk, Adir,
                                 numberFITSdisk, Fdir,
                                 F_TRUE, F_FALSE, err);
     return out;
   }
  ~OSystem() {
    /* best not self->me = ObitSystemShutdown(self->me);*/
    free(self);
  }
};

