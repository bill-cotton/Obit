/* $Id$  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Convolution utilities                  */
/*                                                                    */
/*;  Copyright (C) 2006-2014                                          */
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
#include "ObitTableMFUtil.h"
#include "ObitTableVLUtil.h"
#include "ObitTableFSUtil.h"
%}


%inline %{
/**  Convert an MF table to a VL Table */
void TableMF2VL (ObitTable *in, ObitTable *out, 
		     ObitImage *image, ObitErr *err)
{
  ObitTableMF *inMF=NULL;
  ObitTableVL *outVL=NULL;

  inMF  = ObitTableMFConvert(in);
  outVL = ObitTableVLConvert(out);
  ObitTableMF2VL (inMF, outVL, image, err);
  inMF  = ObitTableMFUnref(inMF);
  outVL = ObitTableVLUnref(outVL);
}  /* end TableMF2VL */

/**  Write human readable version of an MF table to a FILE */
void TableMFPrint (ObitTable *in, ObitImage *image, char *prtFile, 
		   ObitErr *err)
{
  ObitTableMF *inMF=NULL;
  FILE *file;

  if (strncmp (prtFile, "stdout", 6)) {
    ObitTrimTrail(prtFile);  /* Trim any trailing blanks */
    file = fopen (prtFile, "a");
  } else {
    file = stdout;
  }
  inMF  = ObitTableMFConvert(in);
  ObitTableMFPrint (inMF, image, file, err);
  if (strncmp (prtFile, "stdout", 6)) {
    fclose (file);
  }
  inMF  = ObitTableMFUnref(inMF);
}  /* end TableMFPrint */

/**  Append one VL table to another */
void TableVLAppend (ObitTable *in, ObitTable *out, ObitErr *err)
{
  ObitTableVL *inVL=NULL, *outVL=NULL;

  inVL  = ObitTableVLConvert(in);
  outVL = ObitTableVLConvert(out);
  ObitTableVLAppend (inVL, outVL,  err);
  inVL  = ObitTableVLUnref(inVL);
  outVL = ObitTableVLUnref(outVL);
}  /* end TableVLAppend */

/** Index a VL table */
void TableVLIndex (ObitTable *in, ObitErr *err)
{
  ObitTableVL *inVL=NULL;

  inVL  = ObitTableVLConvert(in);
  ObitTableVLIndex (inVL, err);
  inVL  = ObitTableVLUnref(inVL);
}  /* end TableVLIndex */

/**  Merge overlapping components */
void TableVLMerge (ObitTable *in, ObitErr *err)
{
  ObitTableVL *inVL=NULL;

  inVL  = ObitTableVLConvert(in);
  ObitTableVLMerge (inVL, err);
  inVL  = ObitTableVLUnref(inVL);
}  /* end TableVLMerge */

/**  Select significant components */
void TableVLSelect (ObitTable *in, ObitTable *out, ObitErr *err)
{
  ObitTableVL *inVL=NULL, *outVL=NULL;

  inVL  = ObitTableVLConvert(in);
  outVL = ObitTableVLConvert(out);
  ObitTableVLSelect (inVL, outVL,  err);
  inVL  = ObitTableVLUnref(inVL);
  outVL = ObitTableVLUnref(outVL);
}  /* end TableVLSelect */

/**  Remove entries from a given field */
void TableVLPurge (ObitTable *in, char *field, ObitErr *err)
{
  ObitTableVL *inVL=NULL;

  inVL  = ObitTableVLConvert(in);
  ObitTableVLPurge (inVL, field, err);
  inVL  = ObitTableVLUnref(inVL);
}  /* end TableVLPurge */

/**  Remove redundant entries */
void TableVLRedun (ObitTable *in, ObitTable *out, ObitErr *err)
{
  ObitTableVL *inVL=NULL, *outVL=NULL;

  inVL  = ObitTableVLConvert(in);
  outVL = ObitTableVLConvert(out);
  ObitTableVLRedun (inVL, outVL,  err);
  inVL  = ObitTableVLUnref(inVL);
  outVL = ObitTableVLUnref(outVL);
}  /* end  TableVLRedun*/

/**  Write human readable version of an VL table to a FILE */
void TableVLPrint (ObitTable *in, ObitImage *image, char *prtFile, 
		   ObitErr *err)
{
  FILE *file;
  ObitTableVL *inVL=NULL;

  if (strncmp (prtFile, "stdout", 6)) {
    ObitTrimTrail(prtFile);  /* Trim any trailing blanks */
    file = fopen (prtFile, "a");
  } else {
    file = stdout;
  }
  inVL  = ObitTableVLConvert(in);
  ObitTableVLPrint (inVL, image, file, err);
  if (strncmp (prtFile, "stdout", 6)) {
    fclose (file);
  }
  inVL  = ObitTableVLUnref(inVL);
}  /* end TableVLPrint */

/**  Convert a VL  table to a VZ Table */
ObitTable* TableVL2VZ (ObitTable *in, ObitData *data, 
		     	 ObitErr *err)
{
  ObitTableVL *inVL=NULL;
  ObitTableVZ *outVZ=NULL;

  inVL  = ObitTableVLConvert(in);
  outVZ = ObitTableVL2VZ (inVL, data, err);
  inVL  = ObitTableVLUnref(inVL);
  return (ObitTable*)outVZ;
}  /* end TableVL2VZ */

/**  Select entries in VZ  to a VZ Table */
ObitTable* TableVZSel (ObitTable *in, ObitData *data, 
		       ObitErr *err)
{
  ObitTableVZ *inVZ=NULL;
  ObitTableVZ *outVZ=NULL;

  inVZ  = ObitTableVZConvert(in);
  outVZ = ObitTableVZSel (inVZ, data, err);
  inVZ  = ObitTableVZUnref(inVZ);
  return (ObitTable*)outVZ;
}  /* end TableVL2VZ */

/**  Convert an VL table to a FS Table */
void TableVL2FS (ObitTable *in, ObitImage *out, int FSver, ObitErr *err)
{
  ObitTableVL *inVL=NULL;

  inVL  = ObitTableVLConvert(in);
  ObitTableFSUtilVL2FS (inVL, (ObitData*)out, (olong)FSver, err);
  inVL  = ObitTableVLUnref(inVL);
}  /* end TableVL2FS */

/**  Write human readable version of an FS table to a FILE */
void TableFSPrint (ObitTable *in, ObitImage *image, char *prtFile, 
		   ObitErr *err)
{
  FILE *file;
  ObitTableFS *inFS=NULL;

  if (strncmp (prtFile, "stdout", 6)) {
    ObitTrimTrail(prtFile);  /* Trim any trailing blanks */
    file = fopen (prtFile, "a");
  } else {
    file = stdout;
  }
  inFS  = ObitTableFSConvert(in);
  ObitTableFSPrint (inFS, image, file, err);
  if (strncmp (prtFile, "stdout", 6)) {
    fclose (file);
  }
  inFS  = ObitTableFSUnref(inFS);
}  /* end TableFSPrint */

/**  Append one FS table to another */
void TableFSAppend (ObitTable *in, ObitTable *out, ObitErr *err)
{
  ObitTableFS *inFS=NULL, *outFS=NULL;

  inFS  = ObitTableFSConvert(in);
  outFS = ObitTableFSConvert(out);
  ObitTableFSAppend (inFS, outFS,  err);
  inFS  = ObitTableFSUnref(inFS);
  outFS = ObitTableFSUnref(outFS);
}  /* end TableFSAppend */

/** Index a FS table */
void TableFSIndex (ObitTable *in, ObitErr *err)
{
  ObitTableFS *inFS=NULL;

  inFS  = ObitTableFSConvert(in);
  ObitTableFSIndex (inFS, err);
  inFS  = ObitTableFSUnref(inFS);
}  /* end TableFSIndex */

/**  Merge overlapping components */
void TableFSMerge (ObitTable *in, ObitErr *err)
{
  ObitTableFS *inFS=NULL;

  inFS  = ObitTableFSConvert(in);
  ObitTableFSMerge (inFS, err);
  inFS  = ObitTableFSUnref(inFS);
}  /* end TableFSMerge */

/**  Select significant components */
void TableFSSelect (ObitTable *in, ObitTable *out, ObitErr *err)
{
  ObitTableFS *inFS=NULL, *outFS=NULL;

  inFS  = ObitTableFSConvert(in);
  outFS = ObitTableFSConvert(out);
  ObitTableFSSelect (inFS, outFS,  err);
  inFS  = ObitTableFSUnref(inFS);
  outFS = ObitTableFSUnref(outFS);
}  /* end TableFSSelect */

/**  Remove entries from a given field */
void TableFSPurge (ObitTable *in, char *field, ObitErr *err)
{
  ObitTableFS *inFS=NULL;

  inFS  = ObitTableFSConvert(in);
  ObitTableFSPurge (inFS, field, err);
  inFS  = ObitTableFSUnref(inFS);
}  /* end TableFSPurge */

/**  Remove redundant entries */
void TableFSRedun (ObitTable *in, ObitTable *out, ObitErr *err)
{
  ObitTableFS *inFS=NULL, *outFS=NULL;

  inFS  = ObitTableFSConvert(in);
  outFS = ObitTableFSConvert(out);
  ObitTableFSRedun (inFS, outFS,  err);
  inFS  = ObitTableFSUnref(inFS);
  outFS = ObitTableFSUnref(outFS);
}  /* end  TableFSRedun*/

/**  Extract spectra from image cube */
void TableFSGetSpectrum (ObitTable *in, ObitImage *image, ObitErr *err)
{
  ObitTableFS *inFS=NULL;

  inFS  = ObitTableFSConvert(in);
  ObitTableFSGetSpectrum (inFS, image, err);
  inFS  = ObitTableFSUnref(inFS);
}  /* end  TableFSGetSpectrum*/

/**  Filter and fit spectra */
void TableFSFiltVel (ObitTable *in, ObitImage *image, ObitTable *out, 
                     ObitErr *err)
{
  ObitTableFS *inFS=NULL;
  ObitTableFS *outFS=NULL;

  inFS  = ObitTableFSConvert(in);
  outFS = ObitTableFSConvert(out);
  ObitTableFSFiltVel (inFS, image, outFS, err);
  inFS  = ObitTableFSUnref(inFS);
  outFS = ObitTableFSUnref(outFS);
}  /* end  TableFSFiltVel */

%}
