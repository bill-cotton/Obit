/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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

#include "ObitTableVLUtil.h"
#include "ObitPosLabelUtil.h"
#include "ObitTableUtil.h"
#include "ObitSkyGeom.h"
#include "ObitImageFitData.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableVLUtil.c
 * ObitTableVL class utility function definitions.
 */

/*----------------------Private prototypes---------------------------*/
/** Determine if row deselected */
static gboolean VLdesel (ObitTableVLRow *row);

/** Flag (deselect) a row */
static void VLflag (ObitTableVLRow *row);

/** How far from field center */
static ofloat fieldOff(ObitTableVLRow *row, olong center[2]);

/** Swallow a VZ table */
static void ObitTableVZSelSwallow (ObitTableVZ *in, olong *NumTab, odouble **TabRA, odouble **TabDec,  
				   ofloat **TabFlux,  ObitErr *err);

/**  Determine if there is another VZ entry within a given distance */
static gboolean
ObitTableVZSelCheck (olong irow, ofloat nearest, ofloat ignore,
		     olong NumTab, odouble *TabRA, odouble *TabDec, ofloat *TabFlux);

/** Merge VZ entries within a given distance. */
static gboolean
ObitTableVZSelAver (olong irow, ofloat nearest, ofloat distAvg, ofloat ignore, 
		    odouble *Ra2000, odouble *Dec2000, ofloat *PeakInt,
		    olong NumTab, odouble *TabRA, odouble *TabDec, ofloat *TabFlux);

/** Quantify how crowded this field is */
static gint
ObitTableVZSelCrowd (olong irow, ofloat crowd, ofloat distAvg, 
		     odouble Ra2000, odouble Dec2000, ofloat PeakInt,
		     olong NumTab, odouble *TabRA, odouble *TabDec, ofloat *TabFlux);
/*----------------------Public functions---------------------------*/

/**
 * Print raw contents of VL table
 * \param in        VL table to print
 * \param prtFile   Where to write
 * \param *err      ObitErr error stack.
 */
void ObitTableVLPrint (ObitTableVL *in, ObitImage *image, FILE  *prtFile, 
		       ObitErr *err)
{
  ObitTableVLRow *row = NULL;
  olong i, irow;
  ofloat maj, min, pa;
  gchar rast[19], decst[19], field[9];
  ofloat epeak, errra, errdec, errmaj, errmin, errpa;
  ofloat beam[3], beamas[2], xcell, flux, eflux;
  gchar *routine = "ObitTableVLPrint";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableVLIsA(in), err, 
		       "%s input %s not an VL Table", routine, in->name);
  Obit_return_if_fail (ObitImageIsA(image), err, 
		       "%s image %s not an image", routine, image->name);

  ObitTableVLOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Image stuff */
  xcell = fabs (image->myDesc->cdelt[0]);
  beam[0] = image->myDesc->beamMaj / xcell;  /* cells */
  beam[1] = image->myDesc->beamMin / xcell;
  beam[2] = image->myDesc->beamPA * DG2RAD;
  beamas[0] = image->myDesc->beamMaj * 3600.0;  /* asec */
  beamas[1] = image->myDesc->beamMin * 3600.0;

  fprintf (prtFile,"\n Listing of fitted VL table values\n");
  fprintf (prtFile,"Fitted sizes in asec, Peak, Flux, IRMS in mJy, residual values relative to Peak\n");
  fprintf (prtFile,"Error estimates (asec, mJy, deg) given under value\n");
  fprintf (prtFile,
	   "             RA           Dec          Peak    Flux    IRMS  Fit Maj Fit min   PA    res. RMS res Peak  PixX    PixY   Field\n");

  row  = newObitTableVLRow (in);

  /* Loop over table printing */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableVLReadRow (in, irow, row, err);
   if (err->error) Obit_traceback_msg (err, routine, in->name);
   if (VLdesel(row)) continue;  /* Skip deselected record */

    ObitPosLabelUtilRA2HMS (row->Ra2000, image->myDesc->ctype[0], rast);
    ObitPosLabelUtilDec2DMS (row->Dec2000, image->myDesc->ctype[1], decst);
    maj = row->MajorAxis * 3600.0;
    min = row->MinorAxis * 3600.0;
    pa  = row->PosAngle;
    for (i=0; i<8; i++) field[i] = row->Field[i]; field[i] = 0;

    /* Errors */
    ObitImageFitDataGaussErr (row->PeakInt, 
			      row->MajorAxis/xcell, row->MinorAxis/xcell, row->PosAngle*DG2RAD,
			      row->IRMS, (ofloat*)&beam[0],
			      &epeak, &errra, &errdec, &errmaj, &errmin, &errpa);

    /* Flux */
    flux = row->PeakInt*((maj/beamas[0]) * (min/beamas[1]));
    eflux = epeak * ((maj/beamas[0]) * (min/beamas[1]));

    /* Values */
    fprintf (prtFile," %5d %14s %14s %8.2f %8.2f %7.3f %7.3f %7.3f %6.1f %8.3f %8.3f %7.1f %7.1f %8s\n",
	     irow, rast, decst, row->PeakInt*1000.0, flux*1000.0, row->IRMS*1000.0, maj, min, pa,
	     row->ResRMS/row->PeakInt, row->ResPeak/row->PeakInt,
	     row->CenterX, row->CenterY, field);
    /* errors */
    fprintf (prtFile,"              %6.2f         %6.2f %8.2f %8.2f         %7.3f %7.3f %6.1f \n",
	     errra*xcell*3600.0, errdec*xcell*3600.0, epeak*1000.0, eflux*1000.0,
	     errmaj*xcell*3600.0, errmin*xcell*3600.0, errpa*RAD2DG);
  } /* end loop over table */
  
  /* Close up */
  ObitTableVLClose (in,  err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  /* release row object */
  row = ObitTableVLRowUnref(row);
} /* end ObitTableVLPrint */

/**
 * Copy the contents of table in to the end of out.
 * Drops deselected entries (both AIPS and FITS).
 * \param in     Input VL table
 * \param out    Output VL table
 * \param *err   ObitErr error stack.
 */
void ObitTableVLAppend (ObitTableVL *in, ObitTableVL *out, ObitErr *err)
{
  ObitTableVLRow *row = NULL;
  olong irow, orow;
  gchar *routine = "ObitTableVLAppend";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableVLIsA(in), err, 
		       "%s input %s not an VL Table", routine, in->name);
  Obit_return_if_fail (ObitTableVLIsA(out), err, 
		       "%s output %s not an VL Table", routine, out->name);
  
  /* Open output */
  ObitTableVLOpen (out, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
  row  = newObitTableVLRow (out);
  ObitTableSetRow ((ObitTable*)out, (ObitTableRow*)row, err);

  /* Will not be sorted */
  out->myDesc->sort[0] = 0;
  out->myDesc->sort[1] = 0;

  /* Open input */
  ObitTableVLOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Copy restoring beam */
  out->BeamMajor = in->BeamMajor;
  out->BeamMinor = in->BeamMinor;
  out->BeamPA    = in->BeamPA;

  /* Loop over input Table */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableVLReadRow (in, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (VLdesel(row)) continue;  /* Skip deselected record */
    /* Write row */
    orow = -1;
    ObitTableVLWriteRow (out, orow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, out->name);
  } /* End loop over table */
    
  /* cleanup/close up */
  row = ObitTableVLRowUnref(row);
  ObitTableVLClose (in,  err);
  ObitTableVLClose (out,  err);
} /* end ObitTableVLAppend */

/**
 * Sort to RA order and index the contents of table in
 * \param in        Input VL table
 * \param *err      ObitErr error stack.
 */
void ObitTableVLIndex (ObitTableVL *in, ObitErr *err)
{
  ObitTableVLRow *row = NULL;
  olong irow, i, indx[25], ira;
  gchar *colName = {"RA(2000)"};
  gchar *routine = "ObitTableVLIndex";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableVLIsA(in), err, 
		       "%s input %s not an VL Table", routine, in->name);
  
  /* Sort to RA order */
  ObitTableUtilSort ((ObitTable*)in, colName, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  ObitTableVLOpen (in, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  row  = newObitTableVLRow (in);

  /* Loop over table printing */
  for (i=0; i<24; i++) indx[i] = 1;
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableVLReadRow (in, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (VLdesel(row)) continue;  /* Skip deselected record */

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
  ObitTableVLClose (in,  err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  /* release row object */
  row = ObitTableVLRowUnref(row);
} /* end ObitTableVLIndex */

/**
 * Convert a VL table to a VZ (short) table
 * \param in        Input VL table
 * Control parameters on info object:
 * \li "minFlux"    OBIT_float (1,1,1) Minimum acceptable flux density (Jy)
 * \param data      Data object onto which the output VZ table it to be attached
 *                  Frequency in data header assumed to be that of the VL table
 * \param err       ObitErr error stack.
 */
ObitTableVZ* ObitTableVL2VZ (ObitTableVL *in, ObitData *data, ObitErr *err)
{
  ObitTableVZ *outVZ=NULL;
  ObitTableVZRow *VZrow=NULL;
  ObitTableVLRow *VLrow=NULL;
  ObitImageDesc *imDesc=NULL;
  ObitUVDesc *uvDesc=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong VZver, orow, irow, count = 0;
  odouble RefFreq=1.4e9;
  ofloat  Flux, peak, flux, beamRat, BiasAv, beamMaj, beamMin;
  gchar *tname;
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitTableVL2VZ";

  /* error checks */
  if (err->error) return outVZ;
  g_assert (ObitTableVLIsA(in));
  g_assert (ObitDataIsA(data));

  /* Get frequency */
  if (ObitImageIsA(data)) {
    imDesc = ((ObitImage*)data)->myDesc;
    RefFreq = imDesc->crval[imDesc->jlocf];
  } else if (ObitUVIsA(data)) {
    uvDesc = ((ObitUV*)data)->myDesc;
    RefFreq = uvDesc->crval[uvDesc->jlocf];
  }

  /* Get parameters */
  Flux = 0.0;
  InfoReal.flt = Flux; type = OBIT_float;
  ObitInfoListGetTest(in->info, "minFlux",  &type, dim,  &InfoReal);
  if (type==OBIT_float) Flux = InfoReal.flt;
  else if (type==OBIT_double)  Flux = InfoReal.dbl;

  /* Set bias */
  BiasAv = 0.0;
  if (fabs(RefFreq-1.4e9)<1.0e6) BiasAv = 0.0004;  /* NVSS */

  /* Set output Table */
  VZver = ObitTableListGetHigh (data->tableList, "AIPS VZ") + 1;
  tname = g_strconcat ("VZ table for: ", data->name, NULL);
  outVZ = newObitTableVZValue(tname, data, &VZver, OBIT_IO_WriteOnly, 
			      err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, in->name, outVZ);

 /* Open tables */
  ObitTableVLOpen (in, OBIT_IO_ReadOnly, err);
  ObitTableVZOpen (outVZ, OBIT_IO_WriteOnly, err);
  if (err->error) goto cleanup;
  VLrow  = newObitTableVLRow (in);
  VZrow  = newObitTableVZRow (outVZ);
  ObitTableSetRow ((ObitTable*)outVZ, (ObitTableRow*)VZrow, err);
  if (err->error) goto cleanup;

  /* Get beam */
  beamMaj = in->BeamMajor;
  beamMin = in->BeamMinor;

  /* Set Ref Freq */
  outVZ->refFreq = RefFreq;

  /* Loop over VL Table */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableVLReadRow (in, irow, VLrow, err);
    if (err->error) goto cleanup;

    /* Get total flux density */
    peak = VLrow->PeakInt;
    peak += BiasAv;  /* Bias correct peak flux */
    beamRat = MAX(1.0, (VLrow->MajorAxis/beamMaj) * (VLrow->MinorAxis/beamMin));
    flux = peak*beamRat;

    /* Want this one? */
    if (VLrow->status<0) continue;  /* Skip deselected record */
    if (flux<Flux) continue;        /* Too faint */

    /* Copy row data */
    VZrow->Ra2000  = VLrow->Ra2000;
    VZrow->Dec2000 = VLrow->Dec2000;
    VZrow->PeakInt = flux;
    VZrow->Quality = 0;

    /* Write row */
    orow = -1;
    count++;
    ObitTableVZWriteRow (outVZ, orow, VZrow, err);
    if (err->error) goto cleanup;
  } /* end loop over  VL table */

 cleanup:
  /* Close */
  ObitTableVLClose (in,  err);
  ObitTableVZClose (outVZ,  err);
  /* release objects */
  VLrow = ObitTableVLRowUnref(VLrow);
  VZrow = ObitTableVZRowUnref(VZrow);
  if (err->error) Obit_traceback_val (err, routine, in->name, outVZ);

  /* Tell how many */
  Obit_log_error(err, OBIT_InfoErr, "%s: %d rows written", 
		 routine, count);

  return outVZ;
} /* end  ObitTableVL2VZ */

/**
 * Select entries in a VZ (short) table
 * Transliterated from the AIPSish SELVZ.FOR:VZTAB
 * \param in        Input VZ table
 * Control parameters on info object:
 * \li "clip"    OBIT_float scalar Minimum acceptable for clipping (Jy)
 *               (used if doClip, no default)
 * \li "nearest" OBIT_float scalar Minimum nearest neighbor to be 
 *               considered isolated (deg)
 *               (used if doIsol, doAver, no default)
 * \li "distAvg" OBIT_float scalar Distance within to average sources (deg)
 *               (used if doAver, no default)
 * \li "ignore"  OBIT_float scalar Minimum acceptable flux density (Jy)
 *               (used if doAver, doIsol, no default)
 * \li "crowd"   OBIT_float scalar  Distance (deg) for crowding test
 *               (used if doAver, default=nearest)
 * \li "doClip"  OBIT_boolean scalar Clipping by minPeak [def FALSE]
 * \li "doIsol"  OBIT_boolean scalar Select isolated entries within nearest 
 *               ignore entries fainter than ignore [def FALSE]
 * \li "doAver"  OBIT_boolean scalar Average entries within distAvg
 *               reject outliers [def FALSE]
 *             If the brightest source within nearest is fainter than 
 *               ignore then the source is passed, else:
 *             If a component is the brightest within nearest 
 *               then all components within distAvg are merged.
 *             If the source has no companions within nearest 
 *               it is passed.
 *             If the component has brighter companions within nearest,
 *                it is dropped.
 * \param data      Data object onto which the output VZ table is to be attached
 * \param *err      ObitErr error stack.
 */
ObitTableVZ* ObitTableVZSel (ObitTableVZ *in,  ObitData *data, ObitErr *err)
{
  ObitTableVZ *outVZ=NULL;
  ObitTableVZRow *inVZrow=NULL;
  ObitTableVZRow *outVZrow=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong VZver, orow, irow, count=0;
  olong qual;
  gboolean doClip, doIsol, doAver;
  ofloat clip=0.0, nearest=0.0, distAvg=0.0, ignore=0.0, crowd=0.0;
  /* Swallowed input table */
  olong NumTab=0;
  odouble *TabRA=NULL, *TabDec=NULL;
  ofloat *TabFlux=NULL;

  gchar *tname;
  gchar *routine = "ObitTableVZSel";

  /* error checks */
  if (err->error) return outVZ;
  g_assert (ObitTableVZIsA(in));
  g_assert (ObitDataIsA(data));

  /* Get parameters */
  doClip = FALSE;
  ObitInfoListGetTest(in->info, "doClip",  &type, dim,  &doClip);
  doIsol = FALSE;
  ObitInfoListGetTest(in->info, "doIsol",  &type, dim,  &doIsol);
  doAver = FALSE;
  ObitInfoListGetTest(in->info, "doAver",  &type, dim,  &doAver);
  /* Anything requested? */
  Obit_retval_if_fail((doClip || doIsol || doAver), err, outVZ,
		      "%s: No operations selected", routine);

  /* Other parameters without defaults */
  if (doClip) {
    ObitInfoListGet(in->info, "clip",     &type, dim,  &clip,    err);
  }
  if (doIsol) {
    ObitInfoListGet(in->info, "nearest",  &type, dim,  &nearest, err);
    ObitInfoListGet(in->info, "ignore",   &type, dim,  &ignore,  err);
  }
  if (doAver) {
    ObitInfoListGet(in->info, "nearest",  &type, dim,  &nearest, err);
    ObitInfoListGet(in->info, "distAvg",  &type, dim,  &distAvg, err);
    ObitInfoListGet(in->info, "ignore",   &type, dim,  &ignore,  err);
    crowd = nearest;
    ObitInfoListGetTest(in->info, "crowd",    &type, dim,  &crowd);
    if (crowd<=0.0) crowd = nearest;
  }
  if (err->error) Obit_traceback_val (err, routine, in->name, outVZ);

  /* Swallow input table if needed */
  if (doIsol || doAver) {
    ObitTableVZSelSwallow (in, &NumTab,  &TabRA,  &TabDec,  &TabFlux, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, outVZ);
  }

  /* Set output Table */
  VZver = ObitTableListGetHigh (data->tableList, "AIPS VZ") + 1;
  tname = g_strconcat ("VZ table for: ", data->name, NULL);
  outVZ = newObitTableVZValue(tname, data, &VZver, OBIT_IO_WriteOnly, 
			      err);
  g_free (tname);
  if (err->error) goto cleanup;

  /* Open files */
  ObitTableVZOpen (in, OBIT_IO_ReadOnly, err);
  inVZrow  = newObitTableVZRow (in);
  ObitTableVZOpen (outVZ, OBIT_IO_WriteOnly, err);
  outVZrow  = newObitTableVZRow (outVZ);
  ObitTableSetRow ((ObitTable*)outVZ, (ObitTableRow*)outVZrow, err);
  if (err->error) goto cleanup;

  /* Set Ref Freq */
  outVZ->refFreq = in->refFreq;

  /* Loop over VZ Table */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableVZReadRow (in, irow, inVZrow, err);
    if (err->error) goto cleanup;

    /* Want this one? */
    if (inVZrow->status<0) continue;  /* Skip deselected record */
    if (doClip && (inVZrow->PeakInt<clip)) continue;  /* Too faint */

    /* Isolated? */
    if (doIsol && !ObitTableVZSelCheck (irow, nearest, ignore,
					NumTab, TabRA, TabDec, TabFlux)) continue;
    
    if (doAver && !ObitTableVZSelAver (irow, nearest, distAvg, ignore, 
				       &inVZrow->Ra2000, &inVZrow->Dec2000, &inVZrow->PeakInt,
				       NumTab, TabRA, TabDec, TabFlux)) continue;
	
    /* OK if it got here - need quality/crowding code? */
    if (doAver) qual = 
		  ObitTableVZSelCrowd (irow, crowd, distAvg, 
				       inVZrow->Ra2000, inVZrow->Dec2000, inVZrow->PeakInt,
				       NumTab, TabRA, TabDec, TabFlux);
    else qual = inVZrow->Quality;

    /* Copy row data */
    outVZrow->Ra2000  = inVZrow->Ra2000;
    outVZrow->Dec2000 = inVZrow->Dec2000;
    outVZrow->PeakInt = inVZrow->PeakInt;
    outVZrow->Quality = qual;

    /* Write row */
    orow = -1;
    count ++;
    ObitTableVZWriteRow (outVZ, orow, outVZrow, err);
    if (err->error) goto cleanup;


  } /* end loop over  VL table */

 cleanup:
  /* Close */
  ObitTableVZClose (in,  err);
  ObitTableVZClose (outVZ, err);
  /* release objects */
  inVZrow  = ObitTableVZRowUnref(inVZrow);
  outVZrow = ObitTableVZRowUnref(outVZrow);
  /* release memory */
  if (TabRA)   g_free(TabRA);
  if (TabDec)  g_free(TabDec);
  if (TabFlux) g_free(TabFlux);
  if (err->error) Obit_traceback_val (err, routine, in->name, outVZ);

  /* Tell how many */
  Obit_log_error(err, OBIT_InfoErr, "%s: %d rows written", 
		 routine, count);

  return outVZ;
} /* end  ObitTableVZSel */

/**
 * Merge overlapping components from a given field
 * Sums the flux of all components within a specified distance 
 * of each component and then merges components weaker than cutoff
 * of the total with the strongest component.  
 * The resultant position is the weighted average and the 
 * flux is the sum.
 * This should be run on a table with all entries derived from 
 * the same image.
 * Adopted from the AIPSish task VLMRG
 * \param in        Input VL table
 * Control parameters on info object:
 * \li "Radius"  OBIT_float (1,1,1) The radius in pixels
 *     within which to sum components. [def 3.01]
 * \li "Cutoff"  OBIT_float (1,1,1) The minimum acceptable fraction 
 *     of summed flux. [def 0.05]
 * \li "begRow"  OBIT_int (1,1,1) First row to include [def 1]
 * \li "endRow"  OBIT_long (1,1,1) last row to include [def or 0 = all]
 * \param *err      ObitErr error stack.
 */
void ObitTableVLMerge (ObitTableVL *in, ObitErr *err)
{
  ObitTableVLRow *row = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat radius, cutoff, maxdis, vlmax, sum, dis;
  olong brow, erow;
  olong i, j, irow, krow, nrow, nvl, cnt1, cnt2;
  olong *vlnum=NULL, *mainvl=NULL;
  ofloat *xpos=NULL, *ypos=NULL, *flux=NULL, *sflux=NULL;
  odouble ra2000, dec2000;
  ofloat peak, peaky;
  gchar *routine = "ObitTableVLMerge";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableVLIsA(in), err, 
		       "%s input %s not an VL Table", routine, in->name);
  /* Control */
  radius = 3.01;
  ObitInfoListGetTest(in->info, "Radius", &type, dim, &radius);
  maxdis = radius*radius;
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  cutoff = 0.05;
  ObitInfoListGetTest(in->info, "Cutoff", &type, dim, &cutoff);

  /* Open input */
  ObitTableVLOpen (in, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  row  = newObitTableVLRow (in);

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
    ObitTableVLReadRow (in, irow, row, err);
    if (err->error) goto cleanup;
    if (VLdesel(row)) continue;  /* Skip deselected record */
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
      ObitTableVLReadRow (in, krow, row, err);
      if (err->error) goto cleanup;

      if (mainvl[irow] != vlnum[irow]) {
	/* Save some values */
	ra2000  = row->Ra2000;
	dec2000 = row->Dec2000;
	peak    = row->PeakInt;
	
	/* Not a keeper - flag */
	VLflag (row);
	krow = vlnum[irow];
	ObitTableVLWriteRow (in, krow, row, err);
	if (err->error) goto cleanup;
	cnt2++;

	/* Read old main version and update it */
	krow = mainvl[irow];
	ObitTableVLReadRow (in, krow, row, err);
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
	ObitTableVLWriteRow (in, krow, row, err);
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
  ObitTableVLClose (in,  err);
 cleanup:
  /* Deallocate arrays */
  if (vlnum)  g_free(vlnum);
  if (mainvl) g_free(mainvl);
  if (xpos)   g_free(xpos);
  if (ypos)   g_free(ypos);
  if (flux)   g_free(flux);
  if (sflux)  g_free(sflux);
  row = ObitTableVLRowUnref(row);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitTableVLMerge */

/**
 * Select significant components in a table to out
 * Given a table of radii and fractional fluxes, filter out
 * weak sources in the presence of strong.  For any non zero
 * entry in Steps, if the entry has a peak weaker than fraction 
 * of the total within radius then it is not copied to the output.
 * This should be run on a table with all entries derived from 
 * the same image.
 * Adopted from the AIPSish task VLSEL
 * \param in        Input VL table
 * Control parameters on info object:
 * \li "Steps"  OBIT_float (2,?,1) Pairs of values giving
 *      (radius(cells), fraction), empty entries zero filled
 * Example: {{2.0,0.05}, {3.0,0.02}, {4.0,0.01}}
 * \li "BLC"     OBIT_int (2,1,1) Lowest x,y pixel number selected [def 1,1]
 * \li "TRC"     OBIT_long (2,1,1) Highest pixel number selected [def or 0 = all]
 * \li "begRow"  OBIT_int (1,1,1) First row to include [def 1]
 * \li "endRow"  OBIT_long (1,1,1) last row to include [def or 0 = all]
 * \param out       Output VL table
 * \param *err      ObitErr error stack.
 */
void ObitTableVLSelect (ObitTableVL *in, ObitTableVL *out, ObitErr *err)
{
  ObitTableVLRow *row = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, j, k, irow, orow, nrow, nvl, itemp, cnt1, cnt2;
  olong  brow, erow, blc[7], trc[7], nstep, istep;
  ofloat *steps=NULL, *sum=NULL, *mxdis=NULL, dis;
  ofloat *xpos=NULL, *ypos=NULL, *flux=NULL;
  olong *vlnum=NULL;
  gboolean *want=NULL;
  gchar *routine = "ObitTableVLSelect";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableVLIsA(in), err, 
		       "%s input %s not an VL Table", routine, in->name);
  Obit_return_if_fail (ObitTableVLIsA(out), err, 
		       "%s output %s not an VL Table", routine, out->name);
  
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
  ObitTableVLOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  row  = newObitTableVLRow (in);

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
    ObitTableVLReadRow (in, irow, row, err);
    if (err->error) goto cleanup;
    if (VLdesel(row)) continue;  /* Skip deselected records */
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
  ObitTableVLOpen (out, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);

  /* Output will not be sorted */
  out->myDesc->sort[0] = 0;
  out->myDesc->sort[1] = 0;

  /* Loop over selected entries copying */
  cnt1 = cnt2 = 0;
  for (i=0; i<nvl; i++) {
    if (want[i]) {  /* Want it? */
      cnt1++;
      irow = vlnum[i];
      ObitTableVLReadRow (in, irow, row, err);
      if (err->error) goto cleanup;
      if (VLdesel(row)) continue;

      /* Write row */
      orow = -1;
      ObitTableVLWriteRow (out, orow, row, err);
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
  row = ObitTableVLRowUnref(row);
  ObitTableVLClose (in,  err);
  ObitTableVLClose (out,  err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
} /* end ObitTableVLSelect */

/**
 * Remove VL table entries from a given field
 * \param in        Input VL table
 * \param field     Field name to purge (8 char, blank filled)
 * \param *err      ObitErr error stack.
 */
void ObitTableVLPurge (ObitTableVL *in, gchar *field, ObitErr *err)
{
  ObitTableVLRow *row = NULL;
  olong irow;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong brow, erow, ocount;
  gchar *routine = "ObitTableVLPurge";
  
  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableVLIsA(in), err, 
		       "%s input %s not an VL Table", routine, in->name);
  Obit_return_if_fail (field!=NULL, err, "%s field not defined", routine);
  
  /* Open  */
  ObitTableVLOpen (in, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  row  = newObitTableVLRow (in);

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
    ObitTableVLReadRow (in, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (VLdesel(row)) continue;  /* Skip deselected record */

    /* Want to drop this one? */
    if (!strncmp (field, row->Field, 8)) {
      ocount++;
      VLflag (row);

      /* Write row */
      ObitTableVLWriteRow (in, irow, row, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    }
  } /* End loop over table */
    
  /* Tell about results */
  Obit_log_error(err, OBIT_InfoErr, "%s: %d rows dropped", 
		 routine, ocount);

  /* cleanup/close up */
  ObitTableVLClose (in,  err);
  row = ObitTableVLRowUnref(row);
} /* end ObitTableVLPurge */

/**
 * Remove redundant entries from in and write out
 * Search forward from each entry until past time of possible match.
 * If a positional match (maxDist) is found then the one closest
 * to the center of its image (centerPix) is chosen.  The other
 * entry is flagged in the input table.  When the search over the
 * RA range from the original source is finished the final accepted
 * entry is written to the output table.
 * Adopted from the AIPSish task VREDN
 * \param in        Input VL table, will be sorted if needed
 * Control parameters on in->info object:
 * \li "centerPix" OBIT_int (2,1,1) Center pixel in images [def 512,512]
 * \li "maxDist" OBIT_float (1,1,1) How far (") to search for matches [def 15"]
 * \li "begRow"  OBIT_int (1,1,1) First row to include [def 1]
 * \li "endRow"  OBIT_long (1,1,1) last row to include [def or 0 = all]
 * \param out       Output VL table
 * \param *err      ObitErr error stack.
 */
void ObitTableVLRedun (ObitTableVL *in, ObitTableVL *out, ObitErr *err)
{
  ObitTableVLRow *row = NULL, *row2 = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat l, m, maxDist, rad1, rad2;
  olong centerPix[7], brow, erow;
  olong mxbad, ibad, jbad, maxbad, tnobad, *badrow=NULL;
  odouble dist2, ramax, dismax;
  gboolean isbad, toss1, want1;
  olong irow, jrow, orow, nrow, iold, inWant, ocount;
  gchar *routine = "ObitTableVLRedun";
  
  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableVLIsA(in), err, 
		       "%s input %s not an VL Table", routine, in->name);
  Obit_return_if_fail (ObitTableVLIsA(out), err,  
		       "%s output %s not an VL Table", routine, out->name);
  
  /* Control */
  centerPix[0] = centerPix[1] = 512;
  ObitInfoListGetTest(in->info, "centerPix", &type, dim, centerPix);
  maxDist = 15.0;
  ObitInfoListGetTest(in->info, "maxDist", &type, dim, &maxDist);
  maxDist /= 3600.0; /* to degrees */
  ramax = maxDist; /* RA range  */
  dismax = maxDist * maxDist;  /* Square of max distance*/

  /* Sort input to RA order if needed */
  ObitTableUtilSort ((ObitTable*)in, "RA(2000)", FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Open input */
  ObitTableVLOpen (in, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  row  = newObitTableVLRow (in);
  row2 = newObitTableVLRow (in);

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
  mxbad = nrow/10;

  /* Allocate arrays */
  badrow  = g_malloc0(mxbad*sizeof(olong));

  /* Open output */
  ObitTableVLOpen (out, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;

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
    ObitTableVLReadRow (in, irow, row, err);
    if (err->error) goto cleanup;
    if (VLdesel(row)) continue;  /* Skip deselected record */

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
      ObitTableVLReadRow (in, jrow, row2, err);
      if (err->error) goto cleanup;
      if (VLdesel(row)) continue;  /* Skip deselected record */

      /* Is this far enough? */
      if ((row2->Ra2000-row->Ra2000) > ramax) break;

      /* Determine separation */
      ObitSkyGeomShiftXY (row->Ra2000, row->Dec2000, 0.0, 
			  row2->Ra2000, row2->Dec2000, &l, &m);
      dist2 = l*l + m*m;
      toss1 = dist2<dismax;  /* Does one need to go? */

      /* Only if from separate fields or processings */
      toss1 = toss1 && ((strncmp(row->Field, row2->Field, 8)) ||
	(row->JDProcess!=row2->JDProcess));

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
	want1 = want1 || ((rad2==rad1) && (row2->JDProcess>=row->JDProcess));

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

    /* Reread entry to copy if needed */
    if (inWant!=irow) {
      ObitTableVLReadRow (in, inWant, row, err);
      if (err->error) goto cleanup;
      if (VLdesel(row)) continue; 
    }

    /* Write output row */
    orow = -1;
    ocount++;
    ObitTableVLWriteRow (out, orow, row, err);
    if (err->error) goto cleanup;
  } /* End loop over table */

  /* Tell about results */
  Obit_log_error(err, OBIT_InfoErr, "%s:  %d/ %d rows selected ", 
		 routine, ocount, nrow);

  /* Cleanup, shutdown */
 cleanup:
  if (badrow) g_free(badrow);
  ObitTableVLClose (in,  err);
  ObitTableVLClose (out,  err);
  row  = ObitTableVLRowUnref(row);
  row2 = ObitTableVLRowUnref(row2);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
} /* end ObitTableVLRedun */

/**
 * Apply final calibration and do error analysis
 * LOADS OF WORK HERE
 * \param in        Input VL table
 * \param *err      ObitErr error stack.
 */
void ObitTableVLCal (ObitTableVL *in, ObitErr *err)
{
  g_error("stubbed");
} /* end ObitTableVLCal */


/*----------------------Private functions ---------------------------*/
/**
 * Determine if row deselected
 * status <0 or Dec2000<-90.0
 * \param row  Row to test
 */
static gboolean VLdesel (ObitTableVLRow *row)
{
  gboolean out = FALSE;
  out = (row->status<0) || (row->Dec2000<-90.0);

  return out;
} /* end VLdesel */

/**
 * Mark a VL table row as deselected:
 * status =-1,  Dec2000 = -100.0, Ra2000=1000, Peak=0
 * \param row  Row to test
 */
static void VLflag (ObitTableVLRow *row)
{
  row->status  = -1;
  row->Dec2000 = -100.0;
  row->Ra2000  = 1000.0;
  row->PeakInt = 0.0;
} /* end VLflag */

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
static ofloat fieldOff (ObitTableVLRow *row, olong center[2])
{
  ofloat off2;

  off2 = (row->CenterX-center[0])*(row->CenterX-center[0]) + 
    (row->CenterY-center[1])*(row->CenterY-center[1]);

  return off2;
} /* end fieldOff */

/**
 * Swallow a VZ table
 * Arrays are allocates, should be g_freeed when done
 * Transliterated from the AIPSish SELVZ.FOR:GETTBX
 * \param in        Input VZ table
 * \param NumTab   [out] Number of entries
 * \param TabRA    [out] RAs of entries
 * \param TabDec   [out] Declinations of entries
 * \param TabFlux  [out] Flux densities of entries
 * \param err       ObitErr error stack.
 */
static void ObitTableVZSelSwallow (ObitTableVZ *in, olong *NumTab, odouble **TabRA, odouble **TabDec,  
				   ofloat **TabFlux,  ObitErr *err)
{
  olong numTab=0, irow;
  odouble *tabRA=NULL, *tabDec=NULL;
  ofloat  *tabFlux=NULL;
  ObitTableVZRow *VZrow=NULL;
  gchar *routine = "ObitTableVZSelSwallow";
  
  /* initialize output */
  *NumTab  = numTab;
  *TabRA   = tabRA;
  *TabDec  = tabDec;
  *TabFlux = tabFlux;
  
  /* Open input */
  ObitTableVZOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  VZrow  = newObitTableVZRow (in);
  
  /* Allocate arrays */
  numTab = in->myDesc->nrow;
  tabRA   = g_malloc(numTab*sizeof(odouble));
  tabDec  = g_malloc(numTab*sizeof(odouble));
  tabFlux = g_malloc(numTab*sizeof(ofloat));
		     
  /* set output */
  *NumTab  = numTab;
  *TabRA   = tabRA;
  *TabDec  = tabDec;
  *TabFlux = tabFlux;
  
  /* Loop over VZ Table */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableVZReadRow (in, irow, VZrow, err);
    if (err->error) goto cleanup;
    if (VZrow->status<0) continue;  /* Skip deselected record */
    tabRA[irow-1]   = VZrow->Ra2000;
    tabDec[irow-1]  = VZrow->Dec2000;
    tabFlux[irow-1] = VZrow->PeakInt;
  } /* end loop over VZ table */
  
 cleanup:
  ObitTableVZClose (in,  err);
  VZrow = ObitTableVLRowUnref(VZrow);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitTableVZSelSwallow */

/**
 * Determine if there is another entry within a given distance
 * Returns true if component is isolated
 * Transliterated from the AIPSish SELVZ.FOR:CHKTBX
 * \param irow     Table row number
 * \param ignore   Min. brightest  object to consider (MINPEK)
 * \param nearest  Minimum distance (deg) (MXDIST)
 * \param NumTab   Number of entries
 * \param TabRA    RAs of entries
 * \param TabDec   Declinations of entries
 * \param TabFlux  Flux densities of entries
 */
static gboolean
ObitTableVZSelCheck (olong irow, ofloat nearest, ofloat ignore,
		     olong NumTab, odouble *TabRA, odouble *TabDec, ofloat *TabFlux)
{
  gboolean want=TRUE;
  olong i;
  odouble cosdec, radiff, dra, ddec, dif, mxd2;

  mxd2   = nearest * nearest;
  cosdec = cos (TabDec[irow-1] / 57.2957795);
  radiff = nearest / cosdec;

  /* Loop over table */
  for (i=0; i<NumTab; i++) {
    /* Pay attention to this one? */
    if (TabFlux[i] <= ignore) continue;

    dra = (TabRA[i]-TabRA[irow-1]);

    /* Gone far enough? */
    if (dra > radiff) return want;

    /* RA in range? */
    if ((i != (irow-1)) && (fabs(dra) < radiff)) {
      /* check full separation */
      dra *= cosdec;
      ddec = (TabDec[i]-TabDec[irow-1]);
      dif = dra*dra + ddec*ddec;
      if (dif <= mxd2) return FALSE; /* Found one */
    }
  } /* end loop over table */

  return want;
} /* end ObitTableVZSelCheck */

/**
 * Merge entries within a given distance.
 * If this is the strongest source within nearest, all components within 
 * distAvg are summed and the positions weighted averaged.
 * If this is not the strongest source, it is dropped.
 * If nothing exceeds ignore then TRUE is returned.
 * Returns TRUE is this entry wanted.
 * Transliterated from the AIPSish SELVZ.FOR:MRGTBX 
 * \param irow     Table row number
 * \param nearest  Minimum distance (deg) (MXDIST)
 * \param distAvg  Distance to which to sum components. (SUMDIS)
 * \param ignore   Min. brightest  object to consider (MINPEK)
 * \param Ra2000   RA (2000) deg, may be modified output
 * \param Dec2000  Dec (2000) deg, may be modified output
 * \param PeakInt  IPol flux density Jy, may be modified output
 * \param NumTab   Number of entries
 * \param TabRA    RAs of entries
 * \param TabDec   Declinations of entries
 * \param TabFlux  Flux densities of entries
 */
static gboolean
ObitTableVZSelAver (olong irow, ofloat nearest, ofloat distAvg, ofloat ignore, 
		    odouble *Ra2000, odouble *Dec2000, ofloat *PeakInt,
		    olong NumTab, odouble *TabRA, odouble *TabDec, ofloat *TabFlux)
{
  gboolean want=TRUE;
  olong i, beg=0, last=0;
  gboolean isolated=TRUE, isBrightest=TRUE;
  odouble cosdec, radiff, dra, ddec, dif, mxd2, smds2, sumra, sumdec, sumflux;
  ofloat maxFlux;

  mxd2 = nearest * nearest;
  cosdec = cos (TabDec[irow-1] / 57.2957795);
  radiff = nearest / cosdec;
  beg = -1;
  maxFlux = *PeakInt;

  for (i=0; i<NumTab; i++) {  /* Loop over table */

    dra = (TabRA[i]-TabRA[irow-1]);
    /* Gone far enough? */
    if (dra > radiff) break;

    /* RA in range? */
    if ((i != (irow-1)) && ((fabs(dra) < radiff))) {
      /* check full separation */
      dra *= cosdec;
      ddec = (TabDec[i]-TabDec[irow-1]);
      dif = dra*dra + ddec*ddec;
      if (dif <= mxd2) {
	/* Found one */
	maxFlux = MAX (maxFlux, TabFlux[i]);
	isolated = FALSE;
	isBrightest = isBrightest && (*PeakInt > TabFlux[i]);
	if (beg < 0) beg = i;
	last = i;
	/*  debug
	    IF (TabFlux[i] > ignore) {
	    WRITE (MSGTXT,1666) I, TabFlux[i], TabRA[i], TabDec[i]
	    1666             FORMAT (I7,F10.3, 2F12.6)
	    CALL MSGWRT (5)
	    } */
      }
    }
  } /* end loop over table */
  
  /* If it's isolated, you're done */
  if (isolated) return TRUE;

  /* If nothing exceeds the limit, you're done */
  if (maxFlux < ignore) return TRUE;

  /* If it's not the strongest, drop it */
  if (!isBrightest) return FALSE;

  /* Sum everything within dist. */
  sumra   = (*Ra2000)*(*PeakInt);	
  sumdec  = (*Dec2000)*(*PeakInt);
  sumflux = (*PeakInt);
  smds2   = distAvg*distAvg;

  for (i=beg; i<=last; i++) {

    dra = (TabRA[i]-TabRA[irow-1]);
    /* Gone far enough? */;
    if (dra > radiff) break;

    /* RA in range? */
    if ((i != (irow-1)) && (fabs(dra) < radiff)) {
      /* check full separation */
      dra *= cosdec;
      ddec = (TabDec[i]-TabDec[irow-1]);
      dif = dra*dra + ddec*ddec;
      if (dif <= smds2) {
	/* Found one */
	sumra   += TabRA[i] * TabFlux[i];
	sumdec  += TabDec[i] * TabFlux[i];
	sumflux += TabFlux[i];
      }
    }
  } /* end loop */

  /* Weighted average position */
  *Ra2000   = sumra  / sumflux;
  *Dec2000  = sumdec / sumflux;
  *PeakInt  = sumflux;

  return want;
} /* end ObitTableVZSelAver */

/**
 * Quantify how crowded this field is
 *   0 = nothing else within crowd
 *   1 = >10% of peak other flux within crowd
 *   2 = >30% of peak other flux within crowd
 *   3 = >70% of peak other flux within crowd
 *   4 = >100% of peak other flux within crowd
 *   5 = >200% of peak other flux within crowd
 * Arrays are allocates, should be g_freeed when done
 * Transliterated from the AIPSish SELVZ.FOR:HOCROW
 * \param irow     Table row number
 * \param crowd    Minimum distance to consider for crowding (deg) (MXDIST)
 *                 Sources outside this distance will be ignored.
 * \param distAvg  Distance to which to sum components. (SUMDIS)
 * \param Ra2000   RA (2000) deg
 * \param Dec2000  Dec (2000) deg
 * \param PeakInt  IPol flux density Jy
 * \param NumTab   Number of entries
 * \param TabRA    RAs of entries
 * \param TabDec   Declinations of entries
 * \param TabFlux  Flux densities of entries
 */
static gint
ObitTableVZSelCrowd (olong irow, ofloat crowd, ofloat distAvg, 
		     odouble Ra2000, odouble Dec2000, ofloat PeakInt,
		     olong NumTab, odouble *TabRA, odouble *TabDec, ofloat *TabFlux)
{
  olong i, qual=0;
  gboolean isolated=TRUE, isumed;
  ofloat sumFlux, ratio, sum2;
  odouble cosdec, radiff, dra, ddec, dif, mxd2;

  mxd2 = crowd * crowd;
  cosdec = cos (TabDec[irow-1] / 57.2957795);
  radiff = crowd / cosdec;
  sumFlux = 0.0;
  sum2 = TabFlux[irow-1] + 1.0e-20;

  /* Loop over table */
  for (i=0; i<NumTab; i++) {
    dra = (TabRA[i]-Ra2000);

    /* Gone far enough? */
    if (dra > radiff) break;

    /* RA in range? */
    if ((i != (irow-1)) && (fabs(dra) < radiff)) {
      /* check full separation */
      dra *= cosdec;
      ddec = (TabDec[i]-Dec2000);
      dif = dra*dra + ddec*ddec;

      if (dif<=mxd2) {
	/* Is this going to be in the average? */
	dra  = (TabRA[i]-TabRA[irow-1]) * cosdec;
	ddec = (TabDec[i]-TabDec[irow-1]);
	dif = sqrt (dra*dra + ddec*ddec);
	isumed = (dif <= distAvg);
	if (isumed) {
	  sum2 += TabFlux[i];    /* Sum within average */
	} else {
	  sumFlux += TabFlux[i]; /* Sum not within average */
	}
	isolated = FALSE; /* Other things nearby */
      }
    }
  } /* end loop over table */

  /* If it's isolated, you're done */
  if (isolated) return qual;

  /*  Set return value by ratio of fluxes */
  ratio = sumFlux / sum2;
  if (ratio < 0.1) {
    qual = 0;
  }  else if (ratio < 0.3) {
    qual = 1;
  }  else if (ratio < 0.7) {
    qual = 2;
  }  else if (ratio < 1.0) {
    qual = 3;
  } else if (ratio < 2.0) {
    qual = 4;
  } else qual = 5;

  return qual;
} /* end ObitTableVZSelCrowd */
