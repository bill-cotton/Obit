/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitDisplay type                       */
/*                                                                    */
/*;  Copyright (C) 2005-2016                                          */
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
#include "ObitDisplay.h"
#include "ObitImage.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanWindow.h"
%}


%inline %{
extern ObitDisplay* ODisplayCreate(char* name, char* serverURL, ObitErr *err) {
   return  ObitDisplayCreate ((gchar*)name, (gchar*)serverURL, err );
}

extern int ODisplayImage  (ObitDisplay *display, ObitImage *image, ObitErr *err) {
  return (int) ObitDisplayShow (display, (Obit*)image, NULL, 1, err);
} // end ODisplayImage

extern int ODisplayMosaic  (ObitDisplay *display, ObitImageMosaic *mosaic, int field, 
                           ObitErr *err) {
  return (int) ObitDisplayShow (display, (Obit*)mosaic, NULL, (olong)field, err);
} // end ODisplayMosaic

extern int ODisplayImageEdit  (ObitDisplay *display, ObitImage *image, 
                              ObitDConCleanWindow *window, ObitErr *err) {
  return (int) ObitDisplayShow (display, (Obit*)image, window, 1, err);
} // end ODisplayImageEdit

extern int ODisplayMosaicEdit  (ObitDisplay *display, ObitImageMosaic *mosaic, int field, 
                               ObitDConCleanWindow *window, ObitErr *err) {
  return (int) ObitDisplayShow (display, (Obit*)mosaic, window, (olong)field, err);
} // end ODisplayMosaicEdit

extern int ODisplayMarkPos  (ObitDisplay *display, char *pos, ObitErr *err) {
  return (int) ObitDisplayMarkPos (display, (gchar*)pos, err);
} // end ODisplayMarkPos

extern int ODisplayIsA (ObitDisplay* in) {
  return ObitDisplayIsA(in);
} // end  ODisplayIsA 

ObitDisplay* ODisplayRef (ObitDisplay* in) {
  return ObitDisplayRef (in);
} // end ODisplayRef

ObitDisplay* ODisplayUnref (ObitDisplay* in) {
  if (!ObitDisplayIsA(in)) return NULL;
  return ObitDisplayUnref (in);
} // end ODisplayUnref

extern char* ODisplayGetName (ObitDisplay* in) {
  return in->name;
} // end  ODisplayGetName



%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitDisplay *me;
} ODisplay;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDisplay *me;
} ODisplay;

%addmethods ODisplay { 
  ODisplay(char* name, char* serverURL, ObitErr *err) {
     ODisplay *out;
     char *server;
     out = (ODisplay *) malloc(sizeof(ODisplay));
     if ((serverURL==NULL) || (!strncmp(serverURL, "ObitView", 8))) 
        server = "http://localhost:8765/RPC2";
     else server = serverURL;
     if (strcmp(name, "None")) out->me = ODisplayCreate(name, server, err);
     else out->me = NULL;
     return out;
   }
  ~ODisplay() {
    self->me = ODisplayUnref(self->me);
    free(self);
  }
};

