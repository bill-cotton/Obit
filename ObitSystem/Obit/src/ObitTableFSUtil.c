/* $Id:  $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012-2015                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitUtil.h"
#include "ObitTableFSUtil.h"
#include "ObitPosLabelUtil.h"
#include "ObitTableUtil.h"
#include "ObitSkyGeom.h"
#include "ObitFInterpolate.h"
#include "ObitImageFitData.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableFSUtil.c
 * ObitTableFS class utility function definitions.
 */

/*----------------------Private prototypes---------------------------*/
/** Determine if row deselected */
static gboolean FSdesel (ObitTableFSRow *row);

/** Flag (deselect) a row */
static void FSflag (ObitTableFSRow *row);

/** How far from field center */
static ofloat fieldOff(ObitTableFSRow *row, olong center[2]);

/*----------------------Public functions---------------------------*/

/**
 * Print raw contents of FS table
 * \param in        FS table to print
 *   Control parameters in info:
 * \li "minSNR"     OBIT_float (1,1,1) Minimum acceptable SNR [def 5]
 * \param image     Image with description of catalog
 * \param prtFile   Where to write
 * \param *err      ObitErr error stack.
 */
void ObitTableFSPrint (ObitTableFS *in, ObitImage *image, FILE  *prtFile, 
		       ObitErr *err)
{
  ObitTableFSRow *row = NULL;
  olong i, irow, chmax;
  ofloat maj, min, pa, rms;
  gchar rast[19], decst[19], field[9];
  ofloat epeak, errra, errdec, errmaj, errmin, errpa, minSNR;
  odouble glat, glong;
  ofloat beam[3], beamas[2], xcell, flux, eflux;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitTableFSPrint";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableFSIsA(in), err, 
		       "%s input %s not an FS Table", routine, in->name);
  Obit_return_if_fail (ObitImageIsA(image), err, 
		       "%s image %s not an image", routine, image->name);

  ObitTableFSOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Anything to work on? */
  if (in->myDesc->nrow<=0) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: NO entries in FS table", 
		   routine);
    ObitTableFSClose (in,  err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    return;
  }

  /* Get parameters */
  minSNR = 5.0;
  InfoReal.flt = minSNR; type = OBIT_float;
  ObitInfoListGetTest(in->info, "minSNR",  &type, dim,  &InfoReal);
  if (type==OBIT_float) minSNR = InfoReal.flt;
  else if (type==OBIT_double)  minSNR = InfoReal.dbl;

  /* Image stuff */
  xcell = fabs (image->myDesc->cdelt[0]);
  beam[0] = image->myDesc->beamMaj / xcell;  /* cells */
  beam[1] = image->myDesc->beamMin / xcell;
  beam[2] = image->myDesc->beamPA * DG2RAD;
  beamas[0] = image->myDesc->beamMaj * 3600.0;  /* asec */
  beamas[1] = image->myDesc->beamMin * 3600.0;

  fprintf (prtFile,"\n Listing of fitted FS table values\n");
  fprintf (prtFile,"Fitted sizes in asec, Peak in mJy, Velocity in km/s\n");
  fprintf (prtFile,"Error estimates (asec, mJy, deg) given under value\n");
  fprintf (prtFile,"minSNR = %f\n", minSNR);
  fprintf (prtFile,
	   "             RA           Dec        G long  G lat   Peak    Velocity    Width   Fit Maj Fit min   PA    PixX    PixY    PixZ  Field\n");

  row  = newObitTableFSRow (in);

  /* Loop over table printing */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableFSReadRow (in, irow, row, err);
   if (err->error) Obit_traceback_msg (err, routine, in->name);
   if (FSdesel(row)) continue;  /* Skip deselected record */

    ObitPosLabelUtilRA2HMS (row->Ra2000, "RA---SIN", rast);
    ObitPosLabelUtilDec2DMS (row->Dec2000, "DEC--SIN", decst);
    maj = row->MajorAxis * 3600.0;
    min = row->MinorAxis * 3600.0;
    pa  = fmod(row->PosAngle, 360.0);
    for (i=0; i<8; i++) field[i] = row->Field[i]; field[i] = 0;

    /* Galactic coord */
    glat  = row->Dec2000;
    glong = row->Ra2000;
    ObitSkyGeomJtoB (&glong, &glat);
    ObitSkyGeomEq2Gal (&glong, &glat);

    /* Errors */
    chmax = (olong)(row->CenterZ-0.5);
    if (chmax>0) rms = row->RMSCh[chmax];
    else rms = row->PeakInt/5.0;  /*Kludge */
    ObitImageFitDataGaussErr (row->PeakInt, 
			      row->MajorAxis/xcell, row->MinorAxis/xcell, row->PosAngle*DG2RAD,
			      rms, (ofloat*)&beam[0],
			      &epeak, &errra, &errdec, &errmaj, &errmin, &errpa);

    /* Flux - really peak */
    flux  = row->PeakInt;
    eflux = rms;
    /* Use RMS
       eflux = row->IRMS * ((maj/beamas[0]) * (min/beamas[1])); */

    /* minSNR test */
    if (fabs(flux/eflux)<minSNR) continue;

    /* Values */
    fprintf (prtFile," %5d %14s %14s %6.2f %6.2f %8.2f %8.1f %8.3f %7.3f %7.3f %6.1f %7.1f %7.1f %7.1f %8s\n",
	     irow, rast, decst, glong, glat, flux*1000.0, row->Velocity*1.0e-3,row->VelWidth*1.0e-3, maj, min, pa,
	     row->CenterX, row->CenterY, row->CenterZ, field);
    /* errors */
    fprintf (prtFile,"              %7.4f         %6.3f               %8.2f                   %7.3f %7.3f %6.1f \n",
	     errra*xcell*3600.0, errdec*xcell*3600.0, eflux*1000.0, 
	     errmaj*xcell*3600.0, errmin*xcell*3600.0, errpa*RAD2DG);
  } /* end loop over table */
  
  /* Close up */
  ObitTableFSClose (in,  err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  /* release row object */
  row = ObitTableFSRowUnref(row);
} /* end ObitTableFSPrint */

/**
 * Copy the contents of table in to the end of out.
 * Drops deselected entries (both AIPS and FITS).
 * \param in     Input FS table
 * \param out    Output FS table
 * \param *err   ObitErr error stack.
 */
void ObitTableFSAppend (ObitTableFS *in, ObitTableFS *out, ObitErr *err)
{
  ObitTableFSRow *row = NULL;
  olong irow, orow;
  gchar *routine = "ObitTableFSAppend";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableFSIsA(in), err, 
		       "%s input %s not an FS Table", routine, in->name);
  Obit_return_if_fail (ObitTableFSIsA(out), err, 
		       "%s output %s not an FS Table", routine, out->name);
  
  /* Open output */
  ObitTableFSOpen (out, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
  row  = newObitTableFSRow (out);
  ObitTableSetRow ((ObitTable*)out, (ObitTableRow*)row, err);

  /* Will not be sorted */
  out->myDesc->sort[0] = 0;
  out->myDesc->sort[1] = 0;

  /* Open input */
  ObitTableFSOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Copy restoring beam */
  out->BeamMajor = in->BeamMajor;
  out->BeamMinor = in->BeamMinor;
  out->BeamPA    = in->BeamPA;

  /* Loop over input Table */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableFSReadRow (in, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (FSdesel(row)) continue;  /* Skip deselected record */

    /* Make sure RA in range */
    if (row->Ra2000>360.0) row->Ra2000 -= 360.0;
    if (row->Ra2000<0.0)   row->Ra2000 += 360.0;
    
    /* Write row */
    orow = -1;
    ObitTableFSWriteRow (out, orow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, out->name);
  } /* End loop over table */
    
  /* cleanup/close up */
  row = ObitTableFSRowUnref(row);
  ObitTableFSClose (in,  err);
  ObitTableFSClose (out,  err);
} /* end ObitTableFSAppend */

/**
 * Sort to RA order and index the contents of table in
 * \param in        Input FS table
 * \param *err      ObitErr error stack.
 */
void ObitTableFSIndex (ObitTableFS *in, ObitErr *err)
{
  ObitTableFSRow *row = NULL;
  olong irow, i, indx[25], ira;
  gchar *colName = {"RA(2000)"};
  gchar *routine = "ObitTableFSIndex";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableFSIsA(in), err, 
		       "%s input %s not an FS Table", routine, in->name);
  
  /* Sort to RA order */
  ObitTableUtilSort ((ObitTable*)in, colName, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  ObitTableFSOpen (in, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  row  = newObitTableFSRow (in);

  /* Loop over table printing */
  for (i=0; i<24; i++) indx[i] = 1;
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableFSReadRow (in, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (FSdesel(row)) continue;  /* Skip deselected record */

    /* RA bin */
    ira = row->Ra2000/15.0;
    ira = MAX(MIN(ira,23), 0);
    indx[ira+1] = irow;

  } /* End loop over table */

  for (i=1; i<24; i++) indx[i] = MAX (indx[i], indx[i-1]);

  /* Save info */
  in->myStatus = OBIT_Modified;
  in->SortOrder = 1;
  in->numIndexed = in->myDesc->nrow;
  in->index00 = indx[0];
  in->index01 = indx[1];
  in->index02 = indx[2];
  in->index03 = indx[3];
  in->index04 = indx[4];
  in->index05 = indx[5];
  in->index06 = indx[6];
  in->index07 = indx[7];
  in->index08 = indx[8];
  in->index09 = indx[9];
  in->index10 = indx[10];
  in->index11 = indx[11];
  in->index12 = indx[12];
  in->index13 = indx[13];
  in->index14 = indx[14];
  in->index15 = indx[15];
  in->index16 = indx[16];
  in->index17 = indx[17];
  in->index18 = indx[18];
  in->index19 = indx[19];
  in->index20 = indx[20];
  in->index21 = indx[21];
  in->index22 = indx[22];
  in->index23 = indx[23];

  /* Close up */
  ObitTableFSClose (in,  err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  /* release row object */
  row = ObitTableFSRowUnref(row);
} /* end ObitTableFSIndex */

/**
 * Merge overlapping components from a given field
 * Sums the flux of all components within a specified distance 
 * of each component and then merges components weaker than cutoff
 * of the total with the strongest component.  
 * The resultant position is the weighted average and the 
 * flux is the sum.
 * This should be run on a table with all entries derived from 
 * the same image.
 * Adopted from the AIPSish task FSMRG
 * \param in        Input FS table
 * Control parameters on info object:
 * \li "Radius"  OBIT_float (1,1,1) The radius in pixels
 *     within which to sum components. [def 3.01]
 * \li "Cutoff"  OBIT_float (1,1,1) The minimum acceptable fraction 
 *     of summed flux. [def 0.05]
 * \li "begRow"  OBIT_int (1,1,1) First row to include [def 1]
 * \li "endRow"  OBIT_long (1,1,1) last row to include [def or 0 = all]
 * \param *err      ObitErr error stack.
 */
void ObitTableFSMerge (ObitTableFS *in, ObitErr *err)
{
  ObitTableFSRow *row = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat radius, cutoff, maxdis, vlmax, sum, dis;
  olong brow, erow;
  olong i, j, irow, krow, nrow, nvl, cnt1, cnt2;
  olong *vlnum=NULL, *mainvl=NULL;
  ofloat *xpos=NULL, *ypos=NULL, *flux=NULL, *sflux=NULL;
  odouble ra2000, dec2000;
  ofloat peak, peaky;
  gchar *routine = "ObitTableFSMerge";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableFSIsA(in), err, 
		       "%s input %s not an FS Table", routine, in->name);
  /* Control */
  radius = 3.01;
  ObitInfoListGetTest(in->info, "Radius", &type, dim, &radius);
  maxdis = radius*radius;
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  cutoff = 0.05;
  ObitInfoListGetTest(in->info, "Cutoff", &type, dim, &cutoff);

  /* Open input */
  ObitTableFSOpen (in, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  row  = newObitTableFSRow (in);

  /* Rows wanted */
  brow = 1;
  ObitInfoListGetTest(in->info, "begRow", &type, dim, &brow);
  brow = MAX (1, brow);
  erow = in->myDesc->nrow;
  ObitInfoListGetTest(in->info, "endRow", &type, dim, &erow);
  erow = MIN (erow, in->myDesc->nrow);
  if (erow<=0) erow = in->myDesc->nrow;
  nrow = erow-brow+1;

  if (nrow<=0) return;  /* Anything selected? */

  /* Allocate arrays */
  vlnum  = g_malloc0(nrow*sizeof(olong));
  mainvl = g_malloc0(nrow*sizeof(olong));
  xpos   = g_malloc0(nrow*sizeof(ofloat));
  ypos   = g_malloc0(nrow*sizeof(ofloat));
  flux   = g_malloc0(nrow*sizeof(ofloat));
  sflux  = g_malloc0(nrow*sizeof(ofloat));

  /* Loop over input Table - read into internal array */
  nvl = 0;
  for (irow=brow; irow<=erow; irow++) {
    ObitTableFSReadRow (in, irow, row, err);
    if (err->error) goto cleanup;
    if (FSdesel(row)) continue;  /* Skip deselected record */
    vlnum[nvl] = irow;
    xpos[nvl]  = row->CenterX;
    ypos[nvl]  = row->CenterY;
    flux[nvl]  = row->PeakInt;
    nvl++;
  } /* End first loop over table */

  /* Sum fluxes within radius */
  for (i=0; i<nvl; i++) {
    mainvl[i] = i;
    vlmax     = -1.0;
    sum       = 0.0;
    for (j=0; j<nvl; j++) {
      dis = ((xpos[i]-xpos[j])*(xpos[i]-xpos[j]) + 
	     (ypos[i]-ypos[j])*(ypos[i]-ypos[j]));
      if (dis<=maxdis) {
	sum += flux[j];
	/* Save the number of the brightest within radius */
	if (vlmax<flux[j]) {
	  mainvl[i] = vlnum[j];
	  vlmax     = flux[j];
	}
      }
    } /* end inner loop */
    sflux[i] = MAX (sum, 1.0e-10);
  } /* end outer loop */

  /* Determine fractions */
  for (i=0; i<nvl; i++) flux[i] /= sflux[i];

  cnt1 = cnt2 = 0;
  /* Loop over list merging */
  for (irow=0; irow<nvl; irow++) {
    if (flux[irow]<cutoff) {  /* Merge? */
      cnt1++;
      /* Read old */
      krow = vlnum[irow];
      ObitTableFSReadRow (in, krow, row, err);
      if (err->error) goto cleanup;

      if (mainvl[irow] != vlnum[irow]) {
	/* Save some values */
	ra2000  = row->Ra2000;
	dec2000 = row->Dec2000;
	peak    = row->PeakInt;
	
	/* Not a keeper - flag */
	FSflag (row);
	krow = vlnum[irow];
	ObitTableFSWriteRow (in, krow, row, err);
	if (err->error) goto cleanup;
	cnt2++;

	/* Read old main version and update it */
	krow = mainvl[irow];
	ObitTableFSReadRow (in, krow, row, err);
	if (err->error) goto cleanup;
	
	/* New position, peak */
	peaky = peak + row->PeakInt;
	if (peaky>1.0e-10) {
	  /* weighted average */
	  row->Ra2000  = (ra2000*peak  + row->Ra2000*row->PeakInt) / peaky;
	  row->Dec2000 = (dec2000*peak + row->Dec2000*row->PeakInt) / peaky;
	  row->PeakInt = peaky;
	}
	
	/* Update */
	krow = mainvl[irow];
	ObitTableFSWriteRow (in, krow, row, err);
	if (err->error) goto cleanup;
      } else {
	/* This one's a keeper */
      }
   }
  } /* end loop merging */

  /* Tell about results */
  Obit_log_error(err, OBIT_InfoErr, "%s: %d rows affected, %d flagged", 
		 routine, cnt1, cnt2);

 /* cleanup/close up */
  ObitTableFSClose (in,  err);
 cleanup:
  /* Deallocate arrays */
  if (vlnum)  g_free(vlnum);
  if (mainvl) g_free(mainvl);
  if (xpos)   g_free(xpos);
  if (ypos)   g_free(ypos);
  if (flux)   g_free(flux);
  if (sflux)  g_free(sflux);
  row = ObitTableFSRowUnref(row);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitTableFSMerge */

/**
 * Select significant components in a table to out
 * Given a table of radii and fractional fluxes, filter out
 * weak sources in the presence of strong.  For any non zero
 * entry in Steps, if the entry has a peak weaker than fraction 
 * of the total within radius then it is not copied to the output.
 * This should be run on a table with all entries derived from 
 * the same image.
 * Adopted from the AIPSish task FSSEL
 * \param in        Input FS table
 * Control parameters on info object:
 * \li "Steps"  OBIT_float (2,?,1) Pairs of values giving
 *      (radius(cells), fraction), empty entries zero filled
 * Example: {{2.0,0.05}, {3.0,0.02}, {4.0,0.01}}
 * \li "BLC"     OBIT_int (2,1,1) Lowest x,y pixel number selected [def 1,1]
 * \li "TRC"     OBIT_long (2,1,1) Highest pixel number selected [def or 0 = all]
 * \li "begRow"  OBIT_int (1,1,1) First row to include [def 1]
 * \li "endRow"  OBIT_long (1,1,1) last row to include [def or 0 = all]
 * \param out       Output FS table
 * \param *err      ObitErr error stack.
 */
void ObitTableFSSelect (ObitTableFS *in, ObitTableFS *out, ObitErr *err)
{
  ObitTableFSRow *row = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, j, k, irow, orow, nrow, nvl, itemp, cnt1, cnt2;
  olong  brow, erow, blc[7], trc[7], nstep, istep;
  ofloat *steps=NULL, *sum=NULL, *mxdis=NULL, dis;
  ofloat *xpos=NULL, *ypos=NULL, *flux=NULL;
  olong *vlnum=NULL;
  gboolean *want=NULL;
  gchar *routine = "ObitTableFSSelect";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableFSIsA(in), err, 
		       "%s input %s not an FS Table", routine, in->name);
  Obit_return_if_fail (ObitTableFSIsA(out), err, 
		       "%s output %s not an FS Table", routine, out->name);
  
  /* Control */
  blc[0] = blc[1] = 1;
  ObitInfoListGetTest(in->info, "BLC", &type, dim, &blc);
  trc[0] = trc[1] = 10000000;
  ObitInfoListGetTest(in->info, "TRC", &type, dim, &trc);
  if (!ObitInfoListGetP (in->info, "Steps",  &type, dim, (gpointer)&steps)) {
    Obit_log_error(err, OBIT_Error, "%s: No Steps parameter found", routine);
    return;
  }
  nstep = dim[1];
  /* No - really how big? */
  itemp = 0;
  for (istep=0; istep<nstep; istep++) if (steps[istep*2]>0.0) itemp = istep;
  nstep = itemp+1;

 /* Open input */
  ObitTableFSOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  row  = newObitTableFSRow (in);

  /* Rows wanted */
  brow = 1;
  ObitInfoListGetTest(in->info, "begRow", &type, dim, &brow);
  brow = MAX (1, brow);
  erow = in->myDesc->nrow;
  ObitInfoListGetTest(in->info, "endRow", &type, dim, &erow);
  erow = MIN (erow, in->myDesc->nrow);
  if (erow<=0) erow = in->myDesc->nrow;
  nrow = erow-brow+1;

  if (nrow<=0) return;  /* Anything selected? */

   /* Allocate arrays */
  want   = g_malloc0(nrow*sizeof(gboolean));
  vlnum  = g_malloc0(nrow*sizeof(olong));
  xpos   = g_malloc0(nrow*sizeof(ofloat));
  ypos   = g_malloc0(nrow*sizeof(ofloat));
  flux   = g_malloc0(nrow*sizeof(ofloat));
  sum    = g_malloc0(nstep*sizeof(ofloat));
  mxdis  = g_malloc0(nstep*sizeof(ofloat));

  /* Loop over input Table - read into internal array */
  nvl = 0;
  for (irow=brow; irow<=erow; irow++) {
    ObitTableFSReadRow (in, irow, row, err);
    if (err->error) goto cleanup;
    if (FSdesel(row)) continue;  /* Skip deselected records */
    vlnum[nvl] = irow;
    xpos[nvl]  = row->CenterX;
    ypos[nvl]  = row->CenterY;
    flux[nvl]  = row->PeakInt;
    want[nvl]  = ((row->CenterX >= blc[0]) && (row->CenterX <= trc[0])) &&
      ((row->CenterY>=blc[1]) && (row->CenterY<=trc[1]));
    nvl++;
  } /* End first loop over table */

  /* Sum fluxes within radius.  Include components outside of
     guardband (BLC,TRC) in sum */
  for (k=0; k<nstep; k++) mxdis[k] = steps[k*2] * steps[k*2];
  for (i=0; i<nvl; i++) {
    for (k=0; k<nstep; k++) sum[k] = 0.0;
    for (j=0; j<nvl; j++) {
      dis = ((xpos[i]-xpos[j])*(xpos[i]-xpos[j]) + 
	     (ypos[i]-ypos[j])*(ypos[i]-ypos[j]));
      for (k=0; k<nstep; k++) {
	if (dis<=mxdis[k]) sum[k] += flux[j];
      }

    }/* end inner summing loop */
    /* is this one wanted? */
    for (k=0; k<nstep; k++) {
      if ((flux[i]/sum[k]) < steps[k*2+1]) {
	want[i] = FALSE;
	break;
      }
    }
  } /* end outer summing loop */

  /* Open output */
  ObitTableFSOpen (out, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);

  /* copy header */
  out->BeamMajor = in->BeamMajor;
  out->BeamMinor = in->BeamMinor;
  out->BeamPA    = in->BeamPA;
  out->VelRef    = in->VelRef;
  out->VelRPix   = in->VelRPix;
  out->VelDelt   = in->VelDelt;

  /* Output will not be sorted */
  out->myDesc->sort[0] = 0;
  out->myDesc->sort[1] = 0;

  /* Loop over selected entries copying */
  cnt1 = cnt2 = 0;
  for (i=0; i<nvl; i++) {
    if (want[i]) {  /* Want it? */
      cnt1++;
      irow = vlnum[i];
      ObitTableFSReadRow (in, irow, row, err);
      if (err->error) goto cleanup;
      if (FSdesel(row)) continue;

      /* Write row */
      orow = -1;
      ObitTableFSWriteRow (out, orow, row, err);
      if (err->error) goto cleanup;
    } else cnt2++; /* end if wanted */
  } /* End loop copying selected entries */
    
  /* Tell about results */
  Obit_log_error(err, OBIT_InfoErr, "%s: %d rows selected %d dropped", 
		 routine, cnt1, cnt2);

  /* cleanup/close up */
 cleanup:
  /* Deallocate arrays */
  if (vlnum)  g_free(vlnum);
  if (want )  g_free(want);
  if (xpos)   g_free(xpos);
  if (ypos)   g_free(ypos);
  if (flux)   g_free(flux);
  if (sum)    g_free(sum);
  if (mxdis)  g_free(mxdis);
  row = ObitTableFSRowUnref(row);
  ObitTableFSClose (in,  err);
  ObitTableFSClose (out,  err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
} /* end ObitTableFSSelect */

/**
 * Remove FS table entries from a given field
 * \param in        Input FS table
 * \param field     Field name to purge (8 char, blank filled)
 * \param *err      ObitErr error stack.
 */
void ObitTableFSPurge (ObitTableFS *in, gchar *field, ObitErr *err)
{
  ObitTableFSRow *row = NULL;
  olong irow;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong brow, erow, ocount;
  gchar *routine = "ObitTableFSPurge";
  
  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableFSIsA(in), err, 
		       "%s input %s not an FS Table", routine, in->name);
  Obit_return_if_fail (field!=NULL, err, "%s field not defined", routine);
  
  /* Open  */
  ObitTableFSOpen (in, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  row  = newObitTableFSRow (in);

   /* Rows wanted */
  brow = 1;
  ObitInfoListGetTest(in->info, "begRow", &type, dim, &brow);
  brow = MAX (1, brow);
  erow = in->myDesc->nrow;
  ObitInfoListGetTest(in->info, "endRow", &type, dim, &erow);
  erow = MIN (erow, in->myDesc->nrow);
  if (erow<=0) erow = in->myDesc->nrow;

 /* Loop over input Table */
  ocount = 0;
  for (irow=brow; irow<=erow; irow++) {
    ObitTableFSReadRow (in, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (FSdesel(row)) continue;  /* Skip deselected record */

    /* Want to drop this one? */
    if (!strncmp (field, row->Field, 8)) {
      ocount++;
      FSflag (row);

      /* Write row */
      ObitTableFSWriteRow (in, irow, row, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    }
  } /* End loop over table */
    
  /* Tell about results */
  Obit_log_error(err, OBIT_InfoErr, "%s: %d rows dropped", 
		 routine, ocount);

  /* cleanup/close up */
  ObitTableFSClose (in,  err);
  row = ObitTableFSRowUnref(row);
} /* end ObitTableFSPurge */

/**
 * Remove redundant entries from in and write out
 * Search forward from each entry until past time of possible match.
 * If a positional match (maxDist) is found then the one closest
 * to the center of its image (centerPix) is chosen.  The other
 * entry is flagged in the input table.  When the search over the
 * RA range from the original source is finished the final accepted
 * entry is written to the output table.
 * Adopted from the AIPSish task VREDN
 * \param in        Input FS table, will be sorted if needed
 * Control parameters on in->info object:
 * \li "centerPix" OBIT_int (2,1,1) Center pixel in images [def 512,512]
 * \li "maxDist" OBIT_float (1,1,1) How far (") to search for matches [def 15"]
 * \li "begRow"  OBIT_int (1,1,1) First row to include [def 1]
 * \li "endRow"  OBIT_long (1,1,1) last row to include [def or 0 = all]
 * \param out       Output FS table
 * \param *err      ObitErr error stack.
 */
void ObitTableFSRedun (ObitTableFS *in, ObitTableFS *out, ObitErr *err)
{
  ObitTableFSRow *row = NULL, *row2 = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat l, m, maxDist, rad1, rad2;
  olong centerPix[7], brow, erow;
  olong mxbad, ibad, jbad, maxbad, tnobad, *badrow=NULL;
  odouble dist2, ramax, dismax, tramax;
  gboolean isbad, toss1, want1;
  olong irow, jrow, orow, nrow, iold, inWant, ocount;
  gchar Field[9];
  gchar *routine = "ObitTableFSRedun";
  
  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableFSIsA(in), err, 
		       "%s input %s not an FS Table", routine, in->name);
  Obit_return_if_fail (ObitTableFSIsA(out), err,  
		       "%s output %s not an FS Table", routine, out->name);
  
  /* Control */
  centerPix[0] = centerPix[1] = 512;
  ObitInfoListGetTest(in->info, "centerPix", &type, dim, centerPix);
  maxDist = 15.0;
  ObitInfoListGetTest(in->info, "maxDist", &type, dim, &maxDist);
  maxDist /= 3600.0; /* to degrees */
  ramax = maxDist; /* RA range  */
  dismax = maxDist * maxDist;  /* Square of max distance*/

  /* Get rid of any existing output rows */
  ObitTableClearRows ((ObitTable*)out, err);

  /* Sort input to RA order if needed */
  ObitTableUtilSort ((ObitTable*)in, "RA(2000)", FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Open input */
  ObitTableFSOpen (in, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  row  = newObitTableFSRow (in);
  row2 = newObitTableFSRow (in);

  /* Rows wanted */
  brow = 1;
  ObitInfoListGetTest(in->info, "begRow", &type, dim, &brow);
  brow = MAX (1, brow);
  erow = in->myDesc->nrow;
  ObitInfoListGetTest(in->info, "endRow", &type, dim, &erow);
  erow = MIN (erow, in->myDesc->nrow);
  if (erow<=0) erow = in->myDesc->nrow;
  nrow = erow-brow+1;

  if (nrow<=0) return;  /* Anything selected? */

  /* Guess about size of bad row buffer (current list to be dropped) */
  mxbad = nrow;

  /* Allocate arrays */
  badrow  = g_malloc0(mxbad*sizeof(olong));

  /* Open output */
  ObitTableFSOpen (out, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;

  /* copy header */
  out->BeamMajor = in->BeamMajor;
  out->BeamMinor = in->BeamMinor;
  out->BeamPA    = in->BeamPA;
  out->VelRef    = in->VelRef;
  out->VelRPix   = in->VelRPix;
  out->VelDelt   = in->VelDelt;

  /* Output will not be sorted */
  out->myDesc->sort[0] = 0;
  out->myDesc->sort[1] = 0;

  /* Loop over input Table */
  ocount = 0;
  maxbad = 0;
  for (irow=brow; irow<=erow; irow++) {

    /* Flush badrow */
    tnobad = maxbad;
    for (ibad = 0; ibad<maxbad; ibad++) {
      if (ibad>=tnobad) break;
      if (badrow[ibad]<irow) {
	for (jbad=ibad; jbad<maxbad; jbad++) badrow[jbad] = badrow[jbad+1];
	badrow[--tnobad] = 999999999;  /* Pad with big number */
      }
    }
    maxbad = tnobad;

    /* Check if this row in the badrow list */
    isbad = FALSE;
    for (ibad=0; ibad<maxbad; ibad++) {
      if (irow==badrow[ibad]) {isbad=TRUE; break;}
    }
    if (isbad) continue;  /* bad ignore this one */
    
    /* Read row */
    ObitTableFSReadRow (in, irow, row, err);
    if (err->error) goto cleanup;
    if (FSdesel(row)) continue;  /* Skip deselected record */
    /* Have to save Field - pointer in buffer that changes */
    strncpy (Field, row->Field, 8);

    /* Make sure RA in range */
    if (row->Ra2000>360.0) row->Ra2000 -= 360.0;
    if (row->Ra2000<0.0)   row->Ra2000 += 360.0;

    /* How far from center */
    rad1 = fieldOff(row, centerPix);
    iold = irow;

    inWant = irow;  /* keep track or desired output row */

    /* Search following table entries within RA window */
    for (jrow=irow+1; jrow<=erow; jrow++) {
      /* know to be bad? */
      isbad = FALSE;
      for (ibad=0; ibad<maxbad; ibad++) {
	if (jrow==badrow[ibad]) {isbad=TRUE; break;}
      }
      if (isbad) continue;

      /* Read row */
      ObitTableFSReadRow (in, jrow, row2, err);
      if (err->error) goto cleanup;
      if (FSdesel(row)) continue;  /* Skip deselected record */

      /* Make sure RA in range */
      if (row2->Ra2000>360.0) row2->Ra2000 -= 360.0;
      if (row2->Ra2000<0.0)   row2->Ra2000 += 360.0;
      
      /* Is this far enough? */
      tramax = ramax / cos(row->Dec2000*DG2RAD);
      if ((row2->Ra2000-row->Ra2000) > tramax) break;

      /* Determine separation */
      ObitSkyGeomShiftXY (row->Ra2000, row->Dec2000, 0.0, 
			  row2->Ra2000, row2->Dec2000, &l, &m);
      dist2 = l*l + m*m;
      toss1 = dist2<dismax;  /* Does one need to go? */

      /* Only if from separate fields or processings */
      toss1 = toss1 && ((strncmp(Field, row2->Field, 8)));
      /* All may be processed on the same day */
      /* || (row->JDProcess!=row2->JDProcess)); */

      /* Check for complete, exact duplicates */
      toss1 = toss1 || ((row->Ra2000==row2->Ra2000) &&
			(row->Dec2000==row2->Dec2000) &&
			(row->PeakInt==row2->PeakInt));
      
      if (toss1) {
	/* One must die */
	rad2 = fieldOff(row2, centerPix);

	/* Decide if this is the correct field for source */
	want1 = rad2<rad1;
	/* If all else equal take later processing */
	want1 = want1 || ((rad2==rad1));
	/* All may be processed on the same day */
 	/* && (row2->JDProcess>=row->JDProcess));*/

	/* Use this one? */
	if (want1) {
	  inWant = jrow;  /* keep track or desired output row */

	  /* flag old best guess */
	  if (maxbad<mxbad-1) {
	    badrow[maxbad] = iold;
	    maxbad++;
	  } else { /* blew array */
	    Obit_log_error(err, OBIT_Error, "%s: Blew internal array", 
			   routine);
	    goto cleanup;
	  }
	  iold = jrow;
	} /* end if want 1 */

	/* In any case don't reuse this one */
	if (maxbad<mxbad-1) {
	  badrow[maxbad] = jrow;
	  maxbad++;
	} else { /* blew array */
	  Obit_log_error(err, OBIT_Error, "%s: Blew internal array", 
			 routine);
	  goto cleanup;
	}
	
      } /* end of toss1 */
    } /* End forward search for matches */

    /* Reread entry to copy */
    ObitTableFSReadRow (in, inWant, row, err);
    if (err->error) goto cleanup;
    if (FSdesel(row)) continue; 

    /* Make sure RA in range */
    if (row->Ra2000>360.0) row->Ra2000 -= 360.0;
    if (row->Ra2000<0.0)   row->Ra2000 += 360.0;

    /* Write output row */
    orow = -1;
    ocount++;
    ObitTableFSWriteRow (out, orow, row, err);
    if (err->error) goto cleanup;
  } /* End loop over table */

  /* Tell about results */
  Obit_log_error(err, OBIT_InfoErr, "%s:  %d/ %d rows selected ", 
		 routine, ocount, nrow);

  /* Cleanup, shutdown */
 cleanup:
  if (badrow) g_free(badrow);
  ObitTableFSClose (in,  err);
  ObitTableFSClose (out,  err);
  row  = ObitTableFSRowUnref(row);
  row2 = ObitTableFSRowUnref(row2);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
} /* end ObitTableFSRedun */

/**
 * Apply final calibration and do error analysis
 * LOADS OF WORK HERE
 * \param in        Input FS table
 * \param *err      ObitErr error stack.
 */
void ObitTableFSCal (ObitTableFS *in, ObitErr *err)
{
  g_error("stubbed");
} /* end ObitTableFSCal */



/**
 * Convert a VL table to a FS (spectral) table
 * \param in        Input VL table
 * Control parameters on info object:
 * \li "minFlux"    OBIT_float (1,1,1) Minimum acceptable flux density (Jy)
 * \param data      Data object onto which the output FS table it to be attached
 *                  Frequency in data header assumed to be that of the VL table
 * \param FSVer     FS table version number desired
 * \param err       ObitErr error stack.
 */
ObitTableFS* ObitTableFSUtilVL2FS (ObitTableVL *in, ObitData *data, olong FSver,
				   ObitErr *err)
{
  ObitTableFS *outFS=NULL;
  ObitTableFSRow *FSrow=NULL;
  ObitTableVLRow *VLrow=NULL;
  ObitImageDesc *imDesc=NULL;
  ObitUVDesc *uvDesc=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong orow, i,irow, numCh, count = 0;
  odouble RefFreq=1.4e9;
  ofloat  Flux, fblank=ObitMagicF();
  gchar *tname;
  gboolean toEq=FALSE;
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitTableFSUtilVL2FS";

  /* error checks */
  if (err->error) return outFS;
  g_assert (ObitTableVLIsA(in));
  g_assert (ObitDataIsA(data));

  /* Get frequency */
  if (ObitImageIsA(data)) {
    imDesc = ((ObitImage*)data)->myDesc;
    RefFreq = imDesc->crval[imDesc->jlocf];
    /* Need to convert positions to Equatorial? */
    toEq = (!strncmp("GLON", imDesc->ctype[0], 4)) && 
      (!strncmp("GLAT", imDesc->ctype[1], 4));
  } else if (ObitUVIsA(data)) {
    uvDesc = ((ObitUV*)data)->myDesc;
    RefFreq = uvDesc->crval[uvDesc->jlocf];
    /* Need to convert positions to Equatorial? */
    toEq = (!strncmp("GLON", uvDesc->ctype[0], 4)) && 
      (!strncmp("GLAT", uvDesc->ctype[1], 4));
  }

  /* Get parameters */
  Flux = 0.0;
  InfoReal.flt = Flux; type = OBIT_float;
  ObitInfoListGetTest(in->info, "minFlux",  &type, dim,  &InfoReal);
  if (type==OBIT_float) Flux = InfoReal.flt;
  else if (type==OBIT_double)  Flux = InfoReal.dbl;

  /* Set output Table */
  imDesc = ((ObitImage*)data)->myDesc;
  tname = g_strconcat ("FS table for: ", data->name, NULL);
  numCh = imDesc->inaxes[imDesc->jlocf];
  outFS = newObitTableFSValue(tname, data, &FSver, OBIT_IO_ReadWrite, 
			      numCh, err);
  /* Get rid of any existing rows */
  ObitTableClearRows ((ObitTable*)outFS, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, in->name, outFS);

  /* Open tables */
  ObitTableVLOpen (in, OBIT_IO_ReadOnly, err);
  ObitTableFSOpen (outFS, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  VLrow  = newObitTableVLRow (in);
  FSrow  = newObitTableFSRow (outFS);
  ObitTableSetRow ((ObitTable*)outFS, (ObitTableRow*)FSrow, err);
  if (err->error) goto cleanup;

  /* Zero new data */
  FSrow->Velocity = 0.0;
  FSrow->VelWidth = 0.0;
  FSrow->CenterZ  = fblank;
  for (i=0; i<numCh; i++) {
    FSrow->Spectrum[i] = fblank;
    FSrow->RMSCh[i]    = 0.0;
    FSrow->PEAKCh[i]   = 0.0;
  }

  /* copy header */
  outFS->BeamMajor = in->BeamMajor;
  outFS->BeamMinor = in->BeamMinor;
  outFS->BeamPA    = in->BeamPA;

  /* Set Ref Freq 
  outFS->refFreq = RefFreq;*/

  /* Loop over VL Table */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableVLReadRow (in, irow, VLrow, err);
    if (err->error) goto cleanup;

    /* Want this one? */
    if (VLrow->status<0) continue;      /* Skip deselected record */
    if (VLrow->PeakInt<Flux) continue;  /* Too faint? */

    /* Convert to (J2000) Equatorial coordinates? */
    if (toEq) {
      ObitSkyGeomGal2Eq (&VLrow->Ra2000, &VLrow->Dec2000);
      ObitSkyGeomBtoJ   (&VLrow->Ra2000, &VLrow->Dec2000);
    }

    /* Copy row data */
    FSrow->Ra2000    = VLrow->Ra2000;
    FSrow->Dec2000   = VLrow->Dec2000;
    FSrow->PeakInt   = VLrow->PeakInt;
    FSrow->MajorAxis = VLrow->MajorAxis;
    FSrow->MinorAxis = VLrow->MinorAxis;
    FSrow->PosAngle  = VLrow->PosAngle;
    FSrow->QCenter   = VLrow->QCenter;
    FSrow->UCenter   = VLrow->UCenter;
    FSrow->PFlux     = VLrow->PFlux;
    FSrow->IRMS      = VLrow->IRMS;
    FSrow->PolRMS    = VLrow->PolRMS;
    FSrow->ResRMS    = VLrow->ResRMS;
    FSrow->ResPeak   = VLrow->ResPeak;
    FSrow->ResFlux   = VLrow->ResFlux;
    FSrow->CenterX   = VLrow->CenterX;
    FSrow->CenterY   = VLrow->CenterY;
    FSrow->JDProcess = VLrow->JDProcess;
    strncpy(FSrow->Field, VLrow->Field, 8); 

    /* Write row */
    orow = -1;
    count++;
    ObitTableFSWriteRow (outFS, orow, FSrow, err);
    if (err->error) goto cleanup;
  } /* end loop over  VL table */

 cleanup:
  /* Close */
  ObitTableVLClose (in,  err);
  ObitTableFSClose (outFS,  err);
  /* release objects */
  VLrow = ObitTableVLRowUnref(VLrow);
  FSrow = ObitTableFSRowUnref(FSrow);
  if (err->error) Obit_traceback_val (err, routine, in->name, outFS);

  /* Tell how many */
  Obit_log_error(err, OBIT_InfoErr, "%s: %d rows written", 
		 routine, count);

  return outFS;
} /* end  ObitTableFSUtilVL2FS */

/**
 * Extract spectral info from cube, update FS table
 * \param inFS      Input FS table
 * Control parameters on info object:
 * \param im        Data object onto which the output FS table is attached
 *                  Frequency in data header assumed to be that of the VL table
 * \li "RMSsize"    OBIT_float (1,1,1) halfwidth of region to determine RMS [def 50]
 *                  <=0 ->use full image
 * \param err       ObitErr error stack.
 */
void ObitTableFSGetSpectrum(ObitTableFS *inFS, ObitImage *im, ObitErr *err)
{
  ObitTableFSRow *FSrow=NULL;
  ObitFArray *local=NULL;
  ObitIOSize IOBy;
  ObitFInterpolate *interp=NULL;
  ObitInfoType type;
  ofloat RMS, maxAbs, lRMS, lmaxAbs, val, pixel[2];
  odouble coord[2];
  olong irow, RMSsize, iPlane, blc[2], trc[2], Plane[5]={1,1,1,1,1};
  olong pos[2], nx, ny, nch, hwidth;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean toGal;
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitTableFSGetSpectrum";

  /* error checks */
  if (err->error) return;
  g_assert (ObitTableFSIsA(inFS));
  g_assert (ObitImageIsA(im));

  /* Get parameters */
  RMSsize = 50;
  InfoReal.otg = RMSsize; type = OBIT_long;
  ObitInfoListGetTest(im->info, "RMSsize",  &type, dim,  &InfoReal);
  if (type==OBIT_float) RMSsize = InfoReal.flt + 0.5;
  else if (type==OBIT_long)  RMSsize = InfoReal.otg;
  else if (type==OBIT_oint)  RMSsize = InfoReal.itg;

  /* Open table */
  ObitTableFSOpen (inFS, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  FSrow  = newObitTableFSRow (inFS);
  ObitTableSetRow ((ObitTable*)inFS, (ObitTableRow*)FSrow, err);
  if (err->error) goto cleanup;

  /* Anything to work on? */
  if (inFS->myDesc->nrow<=0) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: NO entries in FS table", 
		   routine);
    goto cleanup;
  }

  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (im->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);

  /* Open image */
  if ((ObitImageOpen (im, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, im->name);
    goto cleanup;
  }

  /* Numbers of things */
  nx  = im->myDesc->inaxes[im->myDesc->jlocr];
  ny  = im->myDesc->inaxes[im->myDesc->jlocd];
  nch = im->myDesc->inaxes[im->myDesc->jlocf];

  /* Need to convert position to Galactic? */
  toGal = im->myDesc->ctype[0][0]=='G' && im->myDesc->ctype[0][1]=='L';

  /* Check compatability in no. freq - how many channels actually? */
  if (inFS->numCh<inFS->myDesc->repeat[inFS->SpectrumCol])
    inFS->numCh = inFS->myDesc->repeat[inFS->SpectrumCol];
  Obit_return_if_fail((nch==inFS->numCh), err,
		      "%s: Incompatible no, channels, %d != %d ", 
		      routine, nch, inFS->numCh);

  /* Make interpolator */
  hwidth = 2;
  interp = newObitFInterpolateCreate ("Interpolator", im->image, 
				      im->myDesc, hwidth);
 /* Loop over image planes */
  for (iPlane=1; iPlane<=nch; iPlane++) {
    /* Read input plane */
    Plane[0] = iPlane;
    if ((ObitImageGetPlane (im, NULL, Plane, err)
  	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		     routine, im->name);
      goto cleanup;
    }

    /* Plane statistics */
    RMS    = ObitFArrayRMS (im->image);
    maxAbs = ObitFArrayMaxAbs(im->image, pos);
    
    /* Reset interpolator */
    ObitFInterpolateReplace (interp, im->image);
    
    /* Loop over FS Table */
    for (irow=1; irow<=inFS->myDesc->nrow; irow++) {
      ObitTableFSReadRow (inFS, irow, FSrow, err);
      if (err->error) goto cleanup;
      
      /* Want this one? */
      if (FSrow->status<0) continue;  /* Skip deselected record */
      
      /* Is this in im? */
      coord[0] = FSrow->Ra2000;
      coord[1] = FSrow->Dec2000;

      /* Convert to Galactic? */
      if (toGal) {
	ObitSkyGeomJtoB  (&coord[0], &coord[1]);
	ObitSkyGeomEq2Gal(&coord[0], &coord[1]);
      }

      ObitImageDescGetPixel(im->myDesc, coord, pixel, err);
      if (err->error) goto cleanup;
      if ((pixel[0]<0) || (pixel[1]<0) || (pixel[0]>=nx) || (pixel[1]>=ny)) continue;
      
      val = ObitFInterpolatePosition(interp, coord, err);
      if (err->error) goto cleanup;

      /* Local RMS? */
      if (RMSsize>0) {
	blc[0] = MAX (1, (olong)(pixel[0] - RMSsize + 0.5));
	blc[1] = MAX (1, (olong)(pixel[1] - RMSsize + 0.5));
	trc[0] = MIN (nx-2, (olong)(pixel[0] + RMSsize + 0.5));
	trc[1] = MIN (ny-2, (olong)(pixel[1] + RMSsize + 0.5));
	local = ObitFArraySubArr (im->image, blc, trc, err);
	if (err->error) goto cleanup;
	lRMS    = ObitFArrayRMS (local);
	lmaxAbs = ObitFArrayMaxAbs(local, pos);
	local   = ObitFArrayUnref(local);
      } else {  /* Use plane */
	lRMS    = RMS;
	lmaxAbs = maxAbs;
      }
      
      /* Update */
      FSrow->Spectrum[iPlane-1] = val;
      FSrow->RMSCh[iPlane-1]    = lRMS;
      FSrow->PEAKCh[iPlane-1]   = lmaxAbs;
      
      /* reWrite row */
      ObitTableFSWriteRow (inFS, irow, FSrow, err);
      if (err->error) goto cleanup;
    } /* end loop over  FS table */
  } /* End loop over planes */


cleanup:
  /* Close */
  ObitImageClose (im, err);
  ObitTableFSClose (inFS,  err);
  /* Free image buffer if not memory resident */
  if (im->mySel->FileType!=OBIT_IO_MEM) 
    im->image = ObitFArrayUnref(im->image);
  /* release objects */
  FSrow = ObitTableFSRowUnref(FSrow);
  if (err->error) Obit_traceback_msg (err, routine, im->name);
} /* end ObitTableFSGetSpectrum */

/**
 * Select significant entries and determine velocity.
 * Selected entries written to outFS.
 * \param inFS      Input FS table
 * \param im        Image object onto which the output FS table is attached
 *                  Frequency in data header assumed to be that of the VL table
 *    control parameters on info member:
 * \li "minSNR"     OBIT_float (1,1,1) Minimum acceptable SNR [def 5]
 * \li "minFlux"    OBIT_float (1,1,1) Minimum acceptable flux density (Jy) [def 0]
 * \param outFS     Output FS table
 * \param err       ObitErr error stack.
 */
void ObitTableFSFiltVel(ObitTableFS *inFS, ObitImage *im, ObitTableFS *outFS, 
			ObitErr *err)
{
  ObitTableFSRow *FSrow=NULL;
  olong nx, ny, nch, irow, orow, i, chmax;
  ofloat minSNR, minFlux, median, rms, maxV, tCrpix, *spec=NULL;
  ofloat delnu, refnu, frline, vsign, limit, fblank=ObitMagicF();
  odouble tCrval, dvzero, velite, vel, sum, sum1, sum2;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitImageDesc *imDesc;
  gchar ctype[12];
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitTableFSFiltVel";

  /* error checks */
  if (err->error) return;
  g_assert (ObitTableFSIsA(inFS));
  g_assert (ObitImageIsA(im));

  /* Get parameters */
  minSNR = 5.0;
  InfoReal.flt = minSNR; type = OBIT_float;
  ObitInfoListGetTest(im->info, "minSNR",  &type, dim,  &InfoReal);
  if (type==OBIT_float) minSNR = InfoReal.flt;
  else if (type==OBIT_double)  minSNR = InfoReal.dbl;
  minFlux = 0.0;
  InfoReal.flt = minFlux; type = OBIT_float;
  ObitInfoListGetTest(im->info, "minFlux",  &type, dim,  &InfoReal);
  if (type==OBIT_float) minFlux = InfoReal.flt;
  else if (type==OBIT_double)  minFlux = InfoReal.dbl;

  /* Get rid of any existing output rows */
  ObitTableClearRows ((ObitTable*)outFS, err);

  /* Open tables */
  ObitTableFSOpen (inFS, OBIT_IO_ReadOnly, err);
  ObitTableFSOpen (outFS, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  FSrow  = newObitTableFSRow (inFS);
  ObitTableSetRow ((ObitTable*)outFS, (ObitTableRow*)FSrow, err);
  if (err->error) goto cleanup;

  /* copy header */
  outFS->BeamMajor = inFS->BeamMajor;
  outFS->BeamMinor = inFS->BeamMinor;
  outFS->BeamPA    = inFS->BeamPA;
 
  /* Make sure image descriptor in velocity */
  imDesc = im->myDesc;
  if (strncmp (imDesc->ctype[imDesc->jlocf], "VELO", 4) && 
      strncmp (imDesc->ctype[imDesc->jlocf], "FELO", 4)) {

    /* Set new axis type */
    velite = 2.997924562e8;   /* Speed of light */
    if (imDesc->VelDef==1) {strcpy (ctype, "VELO"); vsign=-1.0;}
    else                   {strcpy (ctype, "FELO"); vsign=+1.0;}
    if (imDesc->VelReference==1)      strcpy (&ctype[4], "-LSR");
    else if (imDesc->VelReference==2) strcpy (&ctype[4], "-HEL");
    else if (imDesc->VelReference==3) strcpy (&ctype[4], "-OBS");
    else                              strcpy (&ctype[4], "-OBS");
    for (i=0; i<8; i++) imDesc->ctype[imDesc->jlocf][i] = ctype[i];
    tCrpix = imDesc->crpix[imDesc->jlocf];
    tCrval = imDesc->crval[imDesc->jlocf];
    imDesc->crpix[imDesc->jlocf] = imDesc->altCrpix;
    imDesc->crval[imDesc->jlocf] = imDesc->altRef;
    delnu   = imDesc->cdelt[imDesc->jlocf];
    refnu   = tCrval;
    frline  = imDesc->altCrpix;
    dvzero  = imDesc->altRef;
    imDesc->cdelt[imDesc->jlocf] = -delnu * (velite + vsign * dvzero) /
      (refnu + delnu * (frline - tCrpix));
    imDesc->altCrpix = tCrpix;
    imDesc->altRef   = tCrval;
  }

  /* Numbers of things */
  nx  = imDesc->inaxes[imDesc->jlocr];
  ny  = imDesc->inaxes[imDesc->jlocd];
  nch = imDesc->inaxes[imDesc->jlocf];

  /* Work array */
  spec = g_malloc0(nch*sizeof(ofloat));

  /* Check compatability in no. freq */
  Obit_return_if_fail((nch==inFS->numCh), err,
		      "%s: Incompatible no, channels, %d != %d ", 
		      routine, nch, inFS->numCh);

  /* Set velocity labeling in FS table */
  outFS->VelRef  = imDesc->crval[imDesc->jlocf];
  outFS->VelRPix = imDesc->crpix[imDesc->jlocf];
  outFS->VelDelt = imDesc->cdelt[imDesc->jlocf];
  
  /* Loop over FS Table */
  for (irow=1; irow<=inFS->myDesc->nrow; irow++) {
    ObitTableFSReadRow (inFS, irow, FSrow, err);
    if (err->error) goto cleanup;
    
    /* Want this one? */
    if (FSrow->status<0) continue;  /* Skip deselected record */

    /* subtract median from spectrum */
    for (i=0; i<nch; i++) spec[i] = FSrow->Spectrum[i];
    median = medianValue(spec, 1, nch);
    if (median==fblank) continue;  /* Any data? */
    /* Get robust RMS */
    rms =  MedianSigma (nch, spec, median);
    for (i=0; i<nch; i++) 
      if (fabs(spec[i]-fblank)>1.0) spec[i] = FSrow->Spectrum[i] - median;
      else spec[i] = fblank;
    /* Find maximum */
    chmax = -1;
    maxV  = -1.0e-20;
    for (i=0; i<nch; i++) {
      if ((fabs(spec[i]-fblank)>10.0) && (spec[i]>maxV)) {
	maxV  = spec[i];
	chmax = i;
      }
    }
    /* Want it? */
    if (chmax<0) continue;
    if (fabs(maxV)<minFlux) continue;
    if (fabs(maxV)<minSNR*rms) continue;
    if (fabs(maxV)<minSNR*FSrow->RMSCh[chmax]) continue;

    /* Find moments above 20% of peak */
    limit = MAX (spec[i]>fabs(maxV)*0.2, 4*MAX(rms,FSrow->RMSCh[chmax]));
    sum1 = sum2 = sum = 0.0;
    for (i=0; i<nch; i++) {
      if ((fabs(spec[i]-fblank)>10.0) && (spec[i]>limit)) {
	sum  += spec[i];
	sum1 += (i - chmax)*spec[i];
	sum2 += (i - chmax)*(i - chmax)*spec[i];
      }
    }

    /* Determine velocity use 1st moment or that of max channel */
    if (sum!=0.0) 
      vel = imDesc->crval[imDesc->jlocf] + imDesc->cdelt[imDesc->jlocf] * 
	(chmax+1+(sum1/sum)-imDesc->crpix[imDesc->jlocf]);
    else
      vel = imDesc->crval[imDesc->jlocf] + imDesc->cdelt[imDesc->jlocf] * 
	(chmax+1-imDesc->crpix[imDesc->jlocf]);
    
    /* Update */
    FSrow->IRMS     = FSrow->RMSCh[chmax];
    FSrow->Velocity = vel;
    /* Second moment */
    if (sum!=0.0) FSrow->VelWidth = 2 * fabs(imDesc->cdelt[imDesc->jlocf]) * sqrt(sum2/sum);  
    FSrow->PeakInt  = maxV;
    FSrow->CenterZ  = (ofloat)(chmax+1.);
    for (i=0; i<nch; i++) FSrow->Spectrum[i] -= median;
    
    /* DEBUG
    vel = imDesc->crval[imDesc->jlocf] + imDesc->cdelt[imDesc->jlocf] * 
      (chmax+1-imDesc->crpix[imDesc->jlocf]);
    fprintf (stderr, "%d %f %lf %lf %lf %lf\n", irow, limit, vel, sum, sum1, sum2); */

    /* reWrite row */
    orow = -1;
    ObitTableFSWriteRow (outFS, orow, FSrow, err);
    if (err->error) goto cleanup;
  } /* end loop over  FS table */


  /* Tell how many */
  Obit_log_error(err, OBIT_InfoErr, "%s: %d rows written", 
		 routine, outFS->myDesc->nrow);
cleanup:
  /* Close */
  ObitTableFSClose (inFS,  err);
  ObitTableFSClose (outFS, err);
  /* release objects */
  FSrow = ObitTableFSRowUnref(FSrow);
  if (spec) g_free(spec);
  if (err->error) Obit_traceback_msg (err, routine, im->name);
} /* end ObitTableFSFiltVel */

/*----------------------Private functions ---------------------------*/
/**
 * Determine if row deselected
 * status <0 or Dec2000<-90.0
 * \param row  Row to test
 */
static gboolean FSdesel (ObitTableFSRow *row)
{
  gboolean out = FALSE;
  out = (row->status<0) || (row->Dec2000<-90.0);

  return out;
} /* end FSdesel */

/**
 * Mark a FS table row as deselected:
 * status =-1,  Dec2000 = -100.0, Ra2000=1000, Peak=0
 * \param row  Row to test
 */
static void FSflag (ObitTableFSRow *row)
{
  row->status  = -1;
  row->Dec2000 = -100.0;
  row->Ra2000  = 1000.0;
  row->PeakInt = 0.0;
} /* end FSflag */

/**
 * Function to determine distance squared from center of the Image
 * (center).  If the position "belongs" to the field the
 * returned value is zero.
 * Adopted from the AIPSish VREDN/FLDOFF which was NVSS specific;
 * concept of "belonging" not implemented.
 * \param row     Row to test
 * \param center center pixel
 * \return distance squared
 */
static ofloat fieldOff (ObitTableFSRow *row, olong center[2])
{
  ofloat off2;

  off2 = (row->CenterX-center[0])*(row->CenterX-center[0]) + 
    (row->CenterY-center[1])*(row->CenterY-center[1]);

  return off2;
} /* end fieldOff */

