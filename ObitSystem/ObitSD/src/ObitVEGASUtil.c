/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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

#include <time.h>
#include <gsl/gsl_randist.h>
#include "ObitThread.h"
#include "ObitOTFUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitVEGASUtil.c
 * ObitOTF class utility function definitions.
 */

/*---------------Private structures----------------*/

/*---------------Private function prototypes----------------*/
/*----------------------Public functions---------------------------*/
/**
 * Average the frequencies in a GBT/VEGAS OTF
 * \param inVEGAS  Input VEGAS, any calibration and flagging should be applied.
 * \param outVEGAS Output VEGAS, must already be defined
 * \param chAvg  Number of channels to average, -1 => all
 * \param err    Error stack
 */
void ObitVEGASUtilAverage(ObitOTF *inOTF, ObitOTF *outOTF, olong chAvg,
			  ObitErr *err) 
{
  const ObitClassInfo *ParentClass;
  ObitIOCode retCode;
  gboolean doCalSelect, done;
  olong firstRec, navg, nchanIn, nchanOut;
  olong ichan, ochan, ochn, ifeed, nfeed, istok, nstok, indx, ondx;
  ObitInfoType type;
  ObitIOAccess access;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat sum, sumWt, *iData, *oData;
  olong NPIO, iran, irec;
  ObitOTFDesc *iDesc, *oDesc;
  ObitHistory *inHist=NULL, *outHist=NULL;
  /* Don't copy Cal and Soln or data or flag tables */
  gchar *exclude[]={"OTFSoln", "OTFCal", "OTFScanData", "OTFFlag", NULL};
  gchar *routine = "ObitVEGASUtilAverage";

    /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));
  g_assert (ObitOTFIsA(outOTF));
  /* Input and output must be different */
  Obit_return_if_fail ((!ObitOTFSame (inOTF, outOTF, err)), err,
		       "%s: Output cannot be the same as the input", routine);

  /* deep copy any base class members */
  ParentClass = ((ObitClassInfo*)inOTF->ClassInfo)->ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (inOTF, outOTF, err);
  outOTF->mySel     = newObitOTFSel (outOTF->name);
  outOTF->tableList = newObitTableList(outOTF->name);
  outOTF->geom      = newObitOTFArrayGeom(outOTF->name);

  /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, (gint32*)dim, 
		      &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* Make sure NPIO the same */
  NPIO = 1000;
  dim[0] = 1;
  ObitInfoListGetTest (inOTF->info, "nRecPIO", &type, dim,  &NPIO);
  ObitInfoListAlwaysPut (inOTF->info, "nRecPIO", OBIT_long, dim,  &NPIO);
  ObitInfoListAlwaysPut (outOTF->info, "nRecPIO", OBIT_long, dim,  &NPIO);

  /* Open Input Data */
  retCode = ObitOTFOpen (inOTF, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) goto cleanup;

  /* How many channels to average? Defaults to all. */
  iDesc = inOTF->myDesc;
  nchanIn = iDesc->inaxes[iDesc->jlocf];
  if (chAvg>0) navg = chAvg;
  else         navg = nchanIn;
  navg     = MAX (1, MIN (navg, nchanIn));
  nchanOut = MAX (1, (olong)(((ofloat)(nchanIn)/navg)+0.9999));

  /* Number of feeds */
  nfeed = iDesc->inaxes[iDesc->jlocfeed];
  /* NUmber of Stokes */
  nstok = iDesc->inaxes[iDesc->jlocs];

  /* Copy descriptor */
  outOTF->myDesc = (gpointer)ObitOTFDescCopy(iDesc, outOTF->myDesc, err); 
  outOTF->myDesc->nrecord = 0;     /* may not copy all */
  oDesc = outOTF->myDesc;
 
  /* Averaging in frequency */
  oDesc->inaxes[oDesc->jlocf]      = nchanOut;
  oDesc->colRepeat[oDesc->ncol-1] /= navg;
  oDesc->cdelt[oDesc->jlocf]      *= navg;
  oDesc->crpix[oDesc->jlocf]      /= navg;
  ObitOTFDescIndex (oDesc);

  /* copy/average Array Geometry - this time with full information
     NO, VEGAS doesn't have one "detector" per frequency 
  outOTF->geom = ObitOTFArrayGeomAver(inOTF->geom, iDesc, outOTF->geom, oDesc, err);
  if (err->error) goto cleanup; */

  /* Open Output Data */
  retCode = ObitOTFOpen (outOTF, OBIT_IO_WriteOnly, err) ;
  if ((retCode != OBIT_IO_OK) || (err->error>0)) goto cleanup;

  /* Copy any history */
  inHist  = newObitHistoryValue("in history", inOTF->info, err);
  outHist = newObitHistoryValue("out history", outOTF->info, err);
  outHist = ObitHistoryCopy (inHist, outHist, err);
  if (err->error)  goto cleanup;
  inHist  = ObitHistoryUnref(inHist);
  outHist = ObitHistoryUnref(outHist);

   /* Copy tables before data  */
  retCode = ObitOTFCopyTables (inOTF, outOTF, exclude, NULL, err);
  if (err->error) goto cleanup;

  /* Close and reopen input to init calibration which will have 
     been disturbed by the table copy */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) goto cleanup;

  retCode = ObitOTFOpen (inOTF, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) goto cleanup;

  /* Loop over data */
  done = (retCode != OBIT_IO_OK);
  while (!done) {

    /* read buffer */#
    retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) goto cleanup;
    done = (retCode == OBIT_IO_EOF); /* done? */
    if (done) break;
    firstRec = inOTF->myDesc->firstRec;

    iData = inOTF->buffer;   /* Input data pointer */
    oData = outOTF->buffer;  /* Output data pointer */
    
    /* How many? */
    outOTF->myDesc->numRecBuff = inOTF->myDesc->numRecBuff;
    
    /* Loop over buffer */
    for (irec=0; irec<iDesc->numRecBuff; irec++) {

      /* Copy random (descriptive parameters */
      for (iran=0; iran<iDesc->numDesc; iran++) oData[iran] = iData[iran];

      /* Average input to output */
      /* Loop over Stokes */
      for (istok=0; istok<nstok; istok++) {
	/* Loop over Feed */
	for (ifeed=0; ifeed<nfeed; ifeed++) {
	  /* Outer loop over channel */
	  ochn = 0;  /* output channel 0-rel index */
	  for (ochan=0; ochan<nchanIn; ochan+=navg) {
	    /* Inner frequency loop summing */
	    sum = sumWt = 0.0;
	    indx = iDesc->ilocdata + istok*iDesc->incs + ifeed*iDesc->incfeed;
	    ondx = oDesc->ilocdata + istok*oDesc->incs + ifeed*oDesc->incfeed;
	    for (ichan=ochan; ichan<ochan+navg; ichan++) {
	      sum   += iData[indx+ichan*iDesc->incf] * iData[indx+ichan*iDesc->incf+1];
	      sumWt += iData[indx+ichan*iDesc->incf+1];
	    } /* end inner loop */
	    /* Save average to outout */
	    if (sumWt>0.0) {
	      oData[ondx+ochn*oDesc->incf]   = sum/sumWt;
	      oData[ondx+ochn*oDesc->incf+1] = sumWt;
	    } else { /* bad */
	      oData[ondx+ochn*oDesc->incf]   = 0.0;
	      oData[ondx+ochn*oDesc->incf+1] = 0.0;
	    }
	    ochn++;
	  } /* end outer channel loop */
	} /* end feed loop */
      } /* end Stokes loop */
      iData += iDesc->lrec;  /* Update buffer pointers */
      oData += oDesc->lrec;
    } /* end loop over buffer */
    
    /* Write buffer */
    if (outOTF->myDesc->numRecBuff>0) retCode = ObitOTFWrite (outOTF, NULL, err);
    if (err->error) goto cleanup;
  } /* end loop over file */
  
 cleanup:
  /* unset output buffer  */
  outOTF->buffer = NULL;
  outOTF->bufferSize = 0;
  retCode = ObitOTFClose (outOTF, err); /* Close output */

  /* Close input */
  retCode = ObitOTFClose (inOTF, err);

  /* cleanup */
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

} /* end ObitOTFUtilAverage */
