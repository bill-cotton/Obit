/* $Id: OPrinter.inc  $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitPrinter type                       */
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
#include "ObitPrinter.h"
#include "ObitImage.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanWindow.h"
%}


%inline %{
extern ObitPrinter* OPrinterCreate(char* name, int isInteractive, 
                                   char *streamname, char* fileName) {
  gboolean lisInteractive = isInteractive!=0;
gchar *fn=NULL;
  FILE *outStream=NULL;
  if (!strncmp(streamname, "stdout", 6)) outStream = stdout;
  if (!strncmp(streamname, "stderr", 6)) outStream = stderr;
  if (strncmp(fileName, "None", 4))     fn = (gchar*)fileName;
  return  ObitPrinterCreate ((gchar*)name, lisInteractive, outStream, fn);
} // end OPrinterCreate

extern int OPrinterOpen  (ObitPrinter *printer, int LinesPerPage, 
                          char *Title1, char *Title2, ObitErr *err) {
  ObitPrinterOpen (printer, (olong)LinesPerPage,
                   (gchar*)Title1, (gchar*)Title2, err);
} // end OPrinterOpen

// Returns quit as int
extern int OPrinterWrite  (ObitPrinter *printer, char *line, ObitErr *err) {
  gboolean quit;
   ObitPrinterWrite (printer, (gchar*)line, &quit, err);
  return (int) quit;
} // end OPrinterWrite

extern void OPrinterClose  (ObitPrinter *printer, ObitErr *err) {
  ObitPrinterClose(printer, err);
} // end OPrinterClose

extern int OPrinterNewPage  (ObitPrinter *printer, ObitErr *err) {
  gboolean quit;
  ObitPrinterNewPage (printer, &quit, err);
  return (int) quit;
} // end OPrinterSetNewPage

extern ObitInfoList* OPrinterGetList (ObitPrinter* in) {
  return ObitInfoListRef(in->info);
}
 
extern int OPrinterIsA (ObitPrinter* in) {
  return ObitPrinterIsA(in);
} // end  OPrinterIsA 

ObitPrinter* OPrinterRef (ObitPrinter* in) {
  return ObitPrinterRef (in);
} // end OPrinterRef

ObitPrinter* OPrinterUnref (ObitPrinter* in) {
  if (!ObitPrinterIsA(in)) return NULL;
  return ObitPrinterUnref (in);
} // end OPrinterUnref

extern char* OPrinterGetName (ObitPrinter* in) {
  return in->name;
} // end  OPrinterGetName



%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitPrinter *me;
} OPrinter;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitPrinter *me;
} OPrinter;

%addmethods OPrinter { 
  OPrinter(char* name, int isInteractive, 
           char *streamname, char* fileName) {
     OPrinter *out;
     gboolean lisInteractive = isInteractive!=0;
     out = (OPrinter *) malloc(sizeof(OPrinter));
     if (strcmp(name, "None")) out->me = OPrinterCreate(name, 
                lisInteractive, streamname, (gchar*)fileName);
     else out->me = NULL;
     return out;
   }
  ~OPrinter() {
    self->me = OPrinterUnref(self->me);
    free(self);
  }
};

