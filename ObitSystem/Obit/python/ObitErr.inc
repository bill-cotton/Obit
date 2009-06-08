/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitErr type                           */
/*                                                                    */
/*;  Copyright (C) 2004,2009                                          */
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
#include "ObitErr.h"
%}

extern ObitErr* ObitErrUnref (ObitErr* in);
extern ObitErr*  ObitErrRef (ObitErr* in);
extern void ObitErrLog (ObitErr* in);
extern void ObitErrClear (ObitErr* in);
extern int ObitErrIsA (ObitErr* in);

%inline %{

extern ObitErr* ObitErrCreate (void) {
  return newObitErr();
}

extern int isError(ObitErr *in) {
   return (int)in->error;
}

extern void SetError(ObitErr *in) {
   in->error = 1;
}

extern void LogError (ObitErr *in, int eCode, char *message) {
   ObitErrCode code;
 /* Should be coordinated with OErr class definition.*/

   switch (eCode) {
   case 0:
     code = OBIT_None;
     break;
   case 1:
     code = OBIT_InfoErr;
     break;
   case 2:
     code = OBIT_InfoWarn;
     break;
   case 3:
     code = OBIT_Traceback;
     break;
   case 4:
     code = OBIT_MildError;
     break;
   case 5:
     code = OBIT_Error;
     break;
   case 6:
     code = OBIT_StrongError;
     break;
   case 7:
     code = OBIT_Fatal;
     break;
   case 8:
     code = OBIT_None;
     break;
   default:
     code = OBIT_None;
   };  

   ObitErrPush (in, code, message);
}  //  end LogError


extern char *OErrMsg(ObitErr *in)
{
  ObitErrCode errLevel;
  gchar *errMsg, *errLevelStr;
  time_t errTimeTag;

/*
 * Human readable versions of the ObitErr codes.
 * Should be coordinated with enum definition.
 */
gchar *ObitErrorLevelString[] = {
  "no message   ", /* OBIT_None        */
  "information  ", /* OBIT_InfoErr     */
  "warning      ", /* OBIT_InfoWarn    */
  "traceback    ", /* OBIT_Traceback   */
  "Mild error   ", /* OBIT_MildError   */
  "Error        ", /* OBIT_Error       */
  "Serious error", /* OBIT_StrongError */
  "Fatal error  ", /* OBIT_Fatal       */
};

  /* error checks */
  g_assert (ObitErrIsA(in));

  ObitErrPop (in, &errLevel, &errMsg, &errTimeTag);
  if(errMsg) {
    gchar *str;
    /* convert error level to something human readable */
    errLevelStr = ObitErrorLevelString[errLevel];
    str = g_strdup_printf("%s: %s", errLevelStr, errMsg);
    g_free(errMsg);
    errMsg = str;
  }
  return errMsg;
}

/* Force an abort as a really heavyhanded interrupt */
extern void Bomb(void) {
  char ct, *PythonSux = NULL;
  ct = PythonSux[-1000000000];
  /* It that doesn't work - maybe this will */
  abort();
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitErr *me;
} OErr;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitErr *me;
} OErr;

%addmethods OErr { 
  OErr(void) {
     OErr *out;
     out = (OErr *) malloc(sizeof(OErr));
     out->me = ObitErrCreate();
    return out;
   }
  ~OErr() {
    self->me = ObitErrUnref(self->me);
    free(self);
  }
};

