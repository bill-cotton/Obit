/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
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
#include "ObitReturn.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitReturn.c
 * ObitReturn function definitions.
 *
 * This Files contains a utility to dump an infoList to a text file.
 */

/*---------------Private function prototypes----------------*/

/** Private: parse next entry */
static ObitIOCode ObitReturnEntry(ObitFile *myFile, gchar* name, ObitInfoType type, 
				  gint32 *dim, gpointer data, ObitErr *err);


/*----------------------Public functions---------------------------*/
/**
 * Add return code "retCode" to InfoList and dump to text file
 * \param retCode Program return code, 0=> normal termination
 *        -999 at startup, a positive value for detected runtime errors
 * \param outfile Name of the input text file to write, if NULL return
 * \param list    ObitInfoList with values.
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitReturnDumpRetCode (olong retCode, gchar *outfile, 
				  ObitInfoList *list, ObitErr *err)
{
  ObitIOCode ret = OBIT_IO_SpecErr;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  oint   itemp;
  gchar *routine = "ObitReturnDumpRetCode";

  /* error checks */
  if (err->error) return ret;
  g_assert(ObitInfoListIsA(list));

  /* Add/update retCode on list */
  dim[0] = 1; dim[1] = 1;
  itemp = retCode;
  ObitInfoListPut (list, "retCode", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "retCode", ret);

  /* Dump to file */
  ret = ObitReturnDump (outfile, list, err);
  if (err->error) Obit_traceback_val (err, routine, "retCode", ret);

  return ret;
} /* end ObitReturnDumpRetCode */

/**
 * Dump InfoList to text file
 * \param outfile Name of the input text file to write, if NULL return
 * \param list    ObitInfoList with values.
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitReturnDump(gchar *outfile, ObitInfoList *list, ObitErr *err)
{
  ObitFile *myFile=NULL;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean OK = TRUE;
  gchar *nameP;
  olong i;
  gpointer data;
  gchar *routine = "ObitReturnDump";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert(ObitInfoListIsA(list));

  if (!outfile) return OBIT_IO_OK;  /* nothing to do? */

  /* Delete any old versions */
  myFile =  newObitFile(outfile);
  retCode = ObitFileOpen (myFile, outfile, OBIT_IO_WriteOnly, OBIT_IO_Text, 0, err);
  retCode = ObitFileClose (myFile, err);
  myFile = ObitFileZap (myFile, err);
  if (err->error) Obit_traceback_val (err, routine, outfile, retCode);

  /* Open text file */
  myFile =  newObitFile(outfile);
  retCode = ObitFileOpen (myFile, outfile, OBIT_IO_WriteOnly, OBIT_IO_Text, 0, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, outfile, retCode);

  /* Loop through InfoList dumping */
  i = 0;
  while (OK && (retCode==OBIT_IO_OK)) {
    i++;
    OK = ObitInfoListGetNumberP (list, i, &nameP, &type, dim, &data);
    if ((!OK) || (data==NULL)) break;

    /* Dump to file */
    retCode = ObitReturnEntry (myFile, nameP, type, dim, data, err);
    if (retCode!=OBIT_IO_OK) break;

  } /* end loop over file */

  /* Close */
  retCode = ObitFileClose (myFile, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, outfile, retCode);

  /* Cleanup */
  myFile = ObitFileUnref(myFile);
 
  return retCode;
} /* end ObitReturnDump */

/*----------------------Private functions---------------------------*/
/**
 * Dump entry to file
 * \param myFile Open ObitFile to write to
 * \param name    The label (keyword) of the information.
 * \param type    Data type of data element (enum defined in ObitInfoList class.
 * \param dim     Dimensionality of datum. (only 3 dimensions used ).
 *                Note: for strings, the first element is the length in char.
 * \param data Pointer to the data. 
 */
static ObitIOCode ObitReturnEntry(ObitFile *myFile, gchar* name, ObitInfoType type, 
				  gint32 *dim, gpointer data, ObitErr *err)
{
  gchar line[200], typeStr[20];
  olong size, i, j, lstr, nstr;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar    *cdata;
  odouble  *ddata;
  ofloat   *fdata;
  oint     *idata;
  olong     *jdata;
  olong    *kdata;
  gboolean *bdata;
  gchar *routine = "ObitReturnEntry";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitFileIsA(myFile));
  g_assert(name!=NULL);
  g_assert(dim!=NULL);
  g_assert(data!=NULL);

  /* inite output line */
  for (j=0; j<199; j++) line[j] = ' ';  line[j] = 0;

  /* Write header line by type */
  switch (type) {
  case OBIT_string:
    sprintf (typeStr,"Str");
    break;
  case OBIT_oint:
  case OBIT_int:
  case OBIT_long:
    sprintf (typeStr,"Int");
    break;
  case OBIT_bool:
    sprintf (typeStr,"Boo");
    break;
  case OBIT_double:
    sprintf (typeStr,"Dbl");
    break;
  case OBIT_float:
    sprintf (typeStr,"Flt");
    break;
  default:
    break;
  }; /* end switch by type */
  /* Write it */
  sprintf (line,"$Key = %s %s (%d,%d,%d)\n", 
	     name, typeStr, dim[0], dim[1], dim[2]);
  retCode = ObitFileWriteLine (myFile, line, err);
  if (err->error) Obit_traceback_val (err, routine, "Output Dumper", retCode);

  /* How much data? */
  size = MAX (1, dim[0]) * MAX (1, dim[1]) * MAX (1, dim[2]);
  
  /* Dump data by type */
  switch (type) {
  case OBIT_string:
    cdata = (gchar*)data;
    lstr =  MAX (1, dim[0]);
    nstr = size/lstr;
    for (i=0; i<nstr; i++) {
      for (j=0; j<lstr; j++) line[j] = cdata[j];  line[j] = '\n'; line[j+1] = 0;
      cdata += lstr;
      retCode = ObitFileWriteLine (myFile, line, err);
      if (err->error) Obit_traceback_val (err, routine, "Output Dumper", retCode);
    }
    break;
  case OBIT_oint:
    idata = (oint*)data;
    for (i=0; i<size; i++) {
      sprintf (line, "%d \n", idata[i]);
      retCode = ObitFileWriteLine (myFile, line, err);
      if (err->error) Obit_traceback_val (err, routine, "Output Dumper", retCode);
    }
    break;
  case OBIT_long:
    kdata = (olong*)data;
    for (i=0; i<size; i++) {
      sprintf (line, " %d \n", kdata[i]);
      retCode = ObitFileWriteLine (myFile, line, err);
      if (err->error) Obit_traceback_val (err, routine, "Output Dumper", retCode);
    }
    break;
  case OBIT_int:
    jdata = (olong*)data;
    for (i=0; i<size; i++) {
      sprintf (line, "%d \n", jdata[i]);
      retCode = ObitFileWriteLine (myFile, line, err);
      if (err->error) Obit_traceback_val (err, routine, "Output Dumper", retCode);
    }
    break;
  case OBIT_bool:
    bdata = (gboolean*)data;
    for (i=0; i<size; i++) {
      if (bdata[i]) sprintf (line, "T \n");
      else sprintf (line, "F \n");
      retCode = ObitFileWriteLine (myFile, line, err);
      if (err->error) Obit_traceback_val (err, routine, "Output Dumper", retCode);
    }
    break;
  case OBIT_double:
    ddata = (odouble*)data;
    for (i=0; i<size; i++) {
      sprintf (line, "%lf \n", ddata[i]);
      retCode = ObitFileWriteLine (myFile, line, err);
      if (err->error) Obit_traceback_val (err, routine, "Output Dumper", retCode);
    }
    break;
  case OBIT_float:
    fdata = (ofloat*)data;
    for (i=0; i<size; i++) {
      sprintf (line, "%f \n", fdata[i]);
      retCode = ObitFileWriteLine (myFile, line, err);
      if (err->error) Obit_traceback_val (err, routine, "Output Dumper", retCode);
    }
    break;
  default:
    break;
  }; /* end switch by type */
  return retCode;
} /* end ObitReturnEntry */
