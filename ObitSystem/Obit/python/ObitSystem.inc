/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitSystem type                        */
/*                                                                    */
/*;  Copyright (C) 2004-2019                                          */
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
#include "ObitThread.h"
%}

/* A copy of the struct for c */
%{
typedef struct OSystem {
  ObitSystem *me;
} OSystem;
%}
/* And a copy of the struct for swig */
typedef struct OSystem {
  ObitSystem *me;
} OSystem;

%inline %{
void CreateOSystem(PyObject *self, char *pgmName, long pgmNumber, long AIPSuser,
             long numberAIPSdisk, char** AIPSdir, 
	     long numberFITSdisk, char** FITSdir, 
	     long F_TRUE, long F_FALSE, ObitErr *err) {
   gchar **Adir, **Fdir;
   oint lF_TRUE  = F_TRUE;
   oint lF_FALSE = F_FALSE;
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_OSystem, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Create: could not recover c struct");
	return;
   }
   if (numberAIPSdisk>0) Adir = AIPSdir;
   else Adir = NULL;
   if (numberFITSdisk>0) Fdir = FITSdir;
   else Fdir = NULL;
   ((OSystem*)ptr)->me = ObitSystemStartup(pgmName, (olong)pgmNumber, (olong)AIPSuser,
   		        (olong)numberAIPSdisk, Adir, (olong)numberFITSdisk, Fdir,
                        lF_TRUE, lF_FALSE, err);
      return;
} // end CreateOSystem

void DeleteOSystem(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_OSystem, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Delete: could not recover c struct");
	return;
   }
   ((OSystem*)ptr)->me = ObitUnref(((OSystem*)ptr)->me);
} // end DeleteOSystem

ObitSystem* OSystem_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_OSystem, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Get_me: could not recover c struct");
	return NULL;
   }
   return ((OSystem*)ptr)->me;
} // end OSystem_Get_me

extern ObitSystem* Shutdown (ObitSystem* in) {
  return ObitSystemShutdown(in);
} // end Shutdown

extern long SystemIsInit (void) {
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

extern long SystemGetPgmNumber (void) {
  return (long)ObitSystemGetPgmNumber ();
} // end SystemGetPgmNumber

extern void SystemSetPgmNumber (long pgmNumber) {
  ObitSystemSetPgmNumber ((olong)pgmNumber);
} // end SystemSetPgmNumber

extern long SystemGetAIPSuser (void) {
  return (long)ObitSystemGetAIPSuser ();
} // end SystemGetAIPSuser

extern void SystemSetAIPSuser (long user) {
  ObitSystemSetAIPSuser ((olong)user);
} // end SystemGetAIPSuser

extern void MemPrint (void) {
  ObitMemPrint(stdout);
} // end MemPrint

extern void SystemAllowThreads (long nThreads) {
  ObitThread *thread=NULL;
 
  /* Create temporary thread */
  thread = newObitThread ();
  /* Init */
  ObitThreadAllowThreads (thread,(olong)nThreads);
  freeObitThread (thread);  /* Cleanup */
 
} // end SystemAllowThreads

extern long SystemGetNoThreads (void) {
  ObitThread *thread=NULL;
  long nThreads;

  /* Create temporary thread */
  thread = newObitThread ();
  nThreads = (int)ObitThreadNumProc (thread);
  freeObitThread (thread);  /* Cleanup */
  return nThreads;
} // end SystemGetNoThreads
%}

