/* $Id$  */
/*Debugging tool            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005                                               */
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

#include "ObitImage.h"
#include "ObitSystem.h"
#include "ObitAIPSDir.h"
#include "ObitTableCCUtil.h"

/* Program globals */
gchar *pgmName = "testCCMerge";       /* Program name */
gchar *infile  = "testCCMerge.inp";   /* File with program inputs */
gchar *outfile = "testCCMerge.out";   /* File to contain program outputs */
olong  pgmNumber=1;     /* Program number (like POPS no.) */
olong  AIPSuser=100;    /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=1;         /* Number of FITS directories */
gchar *FITSdirs[]={"../testIt"}; /* List of FITS data directories */

/* Prototypes */
void DbugR (void);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Test Obit program            */
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  /*ObitInfoType type;*/
  ObitErr      *err= NULL;
  /*gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};;*/
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong        nparm, ver;
  olong         inDisk = 1;
  gchar        *inFile   = "testCCMerge.fits";
  ObitImage    *inImage=NULL;
  ObitTableCC  *inCCTable=NULL, *outCCTable=NULL;

  /* Initialize Obit */
  err = newObitErr();
  ierr = 0;
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
    nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); return ierr;}

 /* Set data */
  inImage = newObitImage("input Image");
  ObitImageSetFITS (inImage, OBIT_IO_byPlane, inDisk, inFile, blc, trc, err);
  /* DbugR (); DEBUG */

  ObitImageFullInstantiate (inImage, TRUE, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); return ierr;}

  /* Input table = CC ver 1 */
  ver   = 1;
  nparm = 0;
  inCCTable = newObitTableCCValue ("Input", (ObitData*)inImage, &ver, 
				    OBIT_IO_ReadOnly, nparm, err);
  ver   = 0;
  nparm = inCCTable-> noParms;
  outCCTable = newObitTableCCValue ("Output", (ObitData*)inImage, &ver, 
				    OBIT_IO_WriteOnly, nparm, err);
  /* Merge to new table */
  ObitTableCCUtilMerge (inCCTable, outCCTable, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); return ierr;}
  ObitErrLog(err);

  /* cleanup */
  inImage     = ObitUnref(inImage);
  inCCTable   = ObitUnref(inCCTable);
  outCCTable  = ObitUnref(outCCTable);
  
  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

/** DEBUGGING */
void DbugR (void) 
{
#include "ObitIOImageFITS.h"
  fitsfile *myFptr;
  gchar *FileName = "../testIt/testCCMerge.fits";
  olong status = 0;
    
   /* DEBUG */
  fits_open_file(&myFptr, FileName, READWRITE, &status);
  fprintf (stderr, "DbugR open status %d\n", status);
  fits_close_file (myFptr, &status);
  fprintf (stderr, "DbugR close status %d\n", status);

} /* end DbugR */
