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

#include "ObitTableMFUtil.h"
#include "ObitPosLabelUtil.h"
#include "ObitPBUtil.h"
#include "ObitUV.h"
#include "ObitUVDesc.h"
#include "ObitFInterpolate.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableMFUtil.c
 * ObitTableMF class utility function definitions.
 */

/*----------------------Private prototypes---------------------------*/
/* Average Polarization over a Gaussian */
static ofloat AverPFlux (ObitFArray *QData, ObitFArray *UData, ofloat pixel[2],
			 ofloat Gaus[3], ofloat PRMS, ofloat beamarea);
/* Polarization bias correction */
static ofloat PBias (ofloat p, ofloat RMS);

/* RMS in box about pixel */
static ofloat RMSbox (ObitFArray *Data, ofloat pixel[2], olong RMSsize, 
		      ObitErr *err);
/*----------------------Public functions---------------------------*/

/**
 * Convert an ObitFitRegionList to entries in an MF table
 * \param MFTable   MF table to write
 * \param regList   RegionList to translate
 * \param image     Image being described
 * \param *err      ObitErr error stack.
 */

void ObitTableMFRegions2MF (ObitTableMF *MFTable, ObitFitRegionList *regList, 
			    ObitImage *image, ObitErr *err)
{
  ObitTableMFRow *row = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong imBLC[7];
  olong i, im, irow;
  ofloat offx, offy;
  ObitFitRegion *reg=NULL;
  ofloat fblank = ObitMagicF();
  ofloat beamarea, cells;
  ofloat bmaj, bmin, bpa, ebmaj, ebmin, ebpa, cbmaj, cbmin, cbpa, dgau[3][3];
  gboolean bad;
  gchar *regname=NULL;
  gchar *routine = "ObitTableMFRegions2MF";

  /* error checks */
  if (err->error) return;

  /* Image information */
  cells = sqrt (fabs(image->myDesc->cdelt[0]*image->myDesc->cdelt[1]));
  /* Subimage in image? */
  for (i=0; i<7; i++) imBLC[i] = 1;
  ObitInfoListGetTest(image->info, "BLC", &type, dim, imBLC);
  for (i=0; i<7; i++) imBLC[i] = MAX (1,imBLC[i]) - 1;

  /* Clean beam area in pixels */
  if (image->myDesc->beamMaj>0.0) {
    beamarea = 1.1331 * (image->myDesc->beamMaj/cells) * 
      (image->myDesc->beamMin/cells);
  } else beamarea = 1.0;

  ObitTableMFOpen (MFTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, MFTable->name);

  /* Mark as unsorted */
  MFTable->myDesc->sort[0] = 0;
  MFTable->myDesc->sort[1] = 0;

  row  = newObitTableMFRow (MFTable);

  /* Loop over regions */
  for (i=1; i<=regList->number; i++) {
    /* Get region */
    regname   = ObitFitRegionName(i);
    reg = ObitFitRegionListFind (regList, regname);
    if (regname) g_free(regname);
    if (!reg) continue;

    /* region information */
    offx = reg->corner[0] + 1 - image->myDesc->crpix[0] + imBLC[0];
    offy = reg->corner[1] + 1 - image->myDesc->crpix[1] + imBLC[1];

    /* Convert to MF table - Loop over models */
    for (im=0; im<reg->nmodel; im++) {

      /* Ignore if zero peak */
      if (reg->models[im]->Peak==0.0) continue;

      /* Deconvolve */
      bmaj = reg->models[im]->parms[0] * cells;
      bmin = reg->models[im]->parms[1] * cells;
      bpa  = reg->models[im]->parms[2] * RAD2DG;
      ebmaj = reg->models[im]->eparms[0] * cells;
      ebmin = reg->models[im]->eparms[1] * cells;
      ebpa  = reg->models[im]->eparms[2] * RAD2DG;
      cbmaj = image->myDesc->beamMaj;
      cbmin = image->myDesc->beamMin;
      cbpa  = image->myDesc->beamPA;
      ObitFitModelDeconGau (bmaj, bmin, bpa, ebmaj, ebmin, ebpa, 
			    cbmaj, cbmin, cbpa, dgau);

      row->plane      = 1;
      row->Peak       = reg->models[im]->Peak;
      row->IFlux      = reg->models[im]->Peak * bmaj*bmin/(cbmaj*cbmin);
      row->DeltaX     = (offx + reg->models[im]->DeltaX) * image->myDesc->cdelt[0];
      row->DeltaY     = (offy + reg->models[im]->DeltaY) * image->myDesc->cdelt[1];
      row->MajorAx    = bmaj;
      row->MinorAx    = bmin;
      row->PosAngle   = bpa;
      row->QFlux      = fblank;
      row->UFlux      = fblank;
      row->VFlux      = fblank;
      row->errPeak    = reg->models[im]->ePeak;
      row->errIFlux   = reg->models[im]->ePeak * bmaj*bmin/(cbmaj*cbmin);
      row->errDeltaX  = reg->models[im]->eDeltaX * cells;
      row->errDeltaY  = reg->models[im]->eDeltaY * cells;
      row->errMajorAx = ebmaj;
      row->errMinorAx = ebmin;
      row->errPosAngle= ebpa;
      row->errQFlux   = fblank;
      row->errUFlux   = fblank;
      row->errVFlux   = fblank;
      row->TypeMod    = 1;
      row->D0Major    = dgau[0][0];
      row->D0Minor    = dgau[0][1];
      row->D0PosAngle = dgau[0][2];
      row->DmMajor    = dgau[2][0];
      row->DmMinor    = dgau[2][1];
      row->DmPosAngle = dgau[2][2];
      row->DpMajor    = dgau[1][0];
      row->DpMinor    = dgau[1][1];
      row->DpPosAngle = dgau[1][2];
      row->ResRMS     = reg->RMSResid;
      row->ResPeak    = reg->peakResid;
      row->ResFlux    = reg->fluxResid / beamarea;
      row->PixelCenterX = reg->corner[0] + reg->models[im]->DeltaX + imBLC[0];
      row->PixelCenterY = reg->corner[1] + reg->models[im]->DeltaY + imBLC[1];
      row->PixelMajorAxis = reg->models[im]->parms[0];
      row->PixelMinorAxis = reg->models[im]->parms[1];
      row->PixelPosAngle  = reg->models[im]->parms[2] * RAD2DG;

      /* Reject screwy fits */
      bad = FALSE;
      bad = fabs(reg->peakResid) > 2.0*reg->peak;
      bad = bad || fabs(reg->RMSResid) > 2.0*reg->peak;

      /* Write row */
      irow = -1;
      if (!bad) ObitTableMFWriteRow (MFTable, irow, row, err);
      if (err->error) Obit_traceback_msg (err, routine, MFTable->name);
    } /* end loop over models */
  } /* end loop over regions */

  /* Close up */
  ObitTableMFClose (MFTable,  err);
  
  /* release row object */
  row = ObitTableMFRowUnref(row);
} /* end ObitTableMFRegions2MF */

/**
 * Convert contents in an ObitTableMF to entries in an VL table
 * \param MFTable   MF table to copy
 * \param VLTable   VL table to write
 * \param regList   RegionList to translate
 * \param image     Image being described, following from info:
 * \li "doPBCorr" OBIT_bool  (1,1,1) If true make PB correction [FALSE]
 * \li "asize"    OBIT_float (1,1,1) Antenna diameter in meters [25]
 *                Used in Primary beam correction.
 * \li "RMSsize"  OBIT_int  (1,1,1) Halfwidth of RMS box [wole image]
 * \param *err      ObitErr error stack.
 */

void ObitTableMF2VL (ObitTableMF *MFTable, ObitTableVL *VLTable, 
		     ObitImage *image, ObitErr *err)
{
  ObitTableMFRow *MFrow = NULL;
  ObitTableVLRow *VLrow = NULL;
  ObitFInterpolate *QInterp = NULL, *UInterp = NULL;
  ObitFArray *IData=NULL, *QData=NULL, *UData=NULL;
  ObitIOSize IOBy;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong  RMSsize, plane[5] = {1,1,1,1,1};
  olong i, irow, orow;
  ofloat pbf, asize, Gaus[3], fblank = ObitMagicF();
  ofloat beamarea, cells, IRMS, QRMS, URMS, QCen, UCen, PFlux;
  gboolean doRMSbox, doPBCorr, GetQU, doJinc, badGain;
  ofloat pixel[2];
  odouble jd, ra, dec, pos[2], zz, dist, raPnt, decPnt, Freq;
  gchar *fieldName, *today;
  gchar *routine = "ObitTableMF2VL";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableMFIsA(MFTable), err, 
		       "%s MFTable %s not an MF Table", routine, MFTable->name);
  Obit_return_if_fail (ObitTableVLIsA(VLTable), err, 
		       "%s VLTable %s not an VL Table", routine, VLTable->name);
  Obit_return_if_fail (ObitImageIsA(image), err, 
		       "%s image %s not an image", routine, image->name);

  /* Image information */
  /* Read IPol image -  Set to read whole plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (image->info, "IOBy", OBIT_long, dim, &IOBy);
  dim[0] = 7;
  for (i=0; i<IM_MAXDIM; i++) {blc[i] = 1; trc[i] = 0;}
  ObitInfoListAlwaysPut (image->info, "BLC", OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (image->info, "TRC", OBIT_long, dim, trc);
  /* Read and extract image data*/
  ObitImageGetPlane (image, NULL, plane, err);
  IData = ObitFArrayCopy(image->image, IData, err);
  if (err->error) goto cleanup;

  /* Other information */
  cells = sqrt (fabs(image->myDesc->cdelt[0]*image->myDesc->cdelt[1]));
  fieldName = image->myDesc->object;   /* field name = "object"*/
  
  /* Processing date - today */
  today = ObitToday();
  ObitUVDescDate2JD (today, &jd);
  g_free(today);
  
  /* Primary beam corr? */
  doPBCorr = FALSE;
  ObitInfoListGetTest(image->info, "doPBCorr",  &type, dim,  &doPBCorr);

  /* antenna size */
  asize = 25.0;
  ObitInfoListGetTest(image->info, "asize",  &type, dim,  &asize);

  /* RMS box */
  RMSsize = 0;
  ObitInfoListGetTest(image->info, "RMSsize",  &type, dim,  &RMSsize);
  doRMSbox = RMSsize>1;
  
  /* Pointing position */
  ObitImageDescGetPoint (image->myDesc, &raPnt, &decPnt);
  raPnt *= DG2RAD; decPnt *= DG2RAD;  /* To radians */
  
  /* Clean beam area in pixels */
  if (image->myDesc->beamMaj>0.0) {
    beamarea = 1.1331 * (image->myDesc->beamMaj/cells) * 
      (image->myDesc->beamMin/cells);
  } else beamarea = 1.0;

  /* which beam model to use */
  Freq = image->myDesc->crval[image->myDesc->jlocf];
  doJinc = (Freq >= 1.0e9);
  
  /* Is there Q,U data? If the Stokes axis yas at least 3 pixels and
     Stokes = I,Q,U get Q, U */
  GetQU = (image->myDesc->jlocs>0) && 
    (image->myDesc->inaxes[image->myDesc->jlocs]>=3) &&
    (image->myDesc->crval[image->myDesc->jlocs]>0.0);
  if (GetQU) {
    /* Q */
    plane[0] = 2;
    ObitImageGetPlane (image, NULL, plane, err);
    QData = ObitFArrayCopy(image->image, QData, err);
    if (err->error) goto cleanup;
    QInterp = newObitFInterpolateCreate ("Q", QData, image->myDesc, 2);
    /* U */
    plane[0] = 3;
    ObitImageGetPlane (image, NULL, plane, err);
    UData = ObitFArrayCopy(image->image, UData, err);
    UInterp = newObitFInterpolateCreate ("U", UData, image->myDesc, 2);
    if (err->error) goto cleanup;
  }
  /* No longer need image buffer */
  image->image = ObitFArrayUnref(image->image);
  
   /* Image statistics if not using a local box */
  if (doRMSbox) {
    IRMS = -1.0;
    QRMS = URMS = fblank;
  } else { /* Whole image */
    IRMS = ObitFArrayRMSQuant(IData);
    if (GetQU) {
      QRMS = ObitFArrayRMSQuant(QData);
      URMS = ObitFArrayRMSQuant(UData);
      QRMS = 0.5 * (QRMS + URMS);
    } else {
      QRMS = fblank;
      URMS = fblank;
    }
  } /* End of whole image RMSes */
  
 /* Open tables */
  ObitTableVLOpen (VLTable, OBIT_IO_ReadWrite, err);
  ObitTableMFOpen (MFTable, OBIT_IO_ReadOnly, err);
  if (err->error) goto cleanup;
  VLrow  = newObitTableVLRow (VLTable);
  ObitTableSetRow ((ObitTable*)VLTable, (ObitTableRow*)VLrow, err);
  if (err->error) goto cleanup;
  MFrow  = newObitTableMFRow (MFTable);
  
  /* Mark as unsorted */
  VLTable->myDesc->sort[0] = 0;
  VLTable->myDesc->sort[1] = 0;

  /* Save restoring beam */
  VLTable->BeamMajor = image->myDesc->beamMaj;
  VLTable->BeamMinor = image->myDesc->beamMin;
  VLTable->BeamPA    = image->myDesc->beamPA;
 
  /* Loop over MF Table */
  for (irow=1; irow<=MFTable->myDesc->nrow; irow++) {
    ObitTableMFReadRow (MFTable, irow, MFrow, err);
    if (err->error) goto cleanup;
    
    /* Get position */
    pixel[0] = MFrow->PixelCenterX;
    pixel[1] = MFrow->PixelCenterY;
    ObitImageDescGetPos(image->myDesc, pixel, pos, err);
    ra  = pos[0] * DG2RAD;
    dec = pos[1] * DG2RAD;
    
    /* Need primary beam corrections? */
    if (doPBCorr) {
      /* separation from pointing center */
      zz = sin (dec) * sin (decPnt) + cos (dec) * cos (decPnt) * cos (ra-raPnt);
      zz = MIN (zz, 1.000);
      dist = acos (zz) * RAD2DG;
      
      /* primary beam correction to flux density */
      if (doJinc) {
	pbf = ObitPBUtilJinc (dist, Freq, asize, 0.0);
      } else {
	pbf = ObitPBUtilPoly (dist, Freq, 0.0);
      } 
      if (pbf!=0.0) pbf = 1.0 / pbf;
      else pbf = 1.0e-20;
    } else pbf = 1.0;  /* No correction */
    
    /* PB Gain out of bounds? (> 5% */
    badGain = pbf > 20.0;
    if (badGain) pbf = 1.0;

    /* Got Stokes? */
    if (GetQU) { 
      pixel[0] = MFrow->PixelCenterX; 
      pixel[1] = MFrow->PixelCenterY;
      QCen  = ObitFInterpolatePixel (QInterp, pixel, err);
      UCen  = ObitFInterpolatePixel (UInterp, pixel, err);
      if (err->error) goto cleanup;
      if (QCen !=fblank) QCen *= pbf;
      if (UCen !=fblank) UCen *= pbf;
      Gaus[0] = MFrow->PixelMajorAxis;
      Gaus[1] = MFrow->PixelMinorAxis;
      Gaus[2] = MFrow->PixelPosAngle;
      PFlux = AverPFlux (QData, UData, pixel, Gaus, QRMS, beamarea);
      if (PFlux!=fblank) PFlux *= pbf;
    } else {
      QCen  = fblank;
      UCen  = fblank;
      PFlux = fblank;
    }
    /* Local RMSes? */
    if (doRMSbox) {
      IRMS = RMSbox (IData, pixel, RMSsize, err);
      if (GetQU) {
	QRMS = RMSbox (QData, pixel, RMSsize, err);
	URMS = RMSbox (UData, pixel, RMSsize, err);
	if (err->error) goto cleanup;
	QRMS = 0.5 * (QRMS + URMS);
      }
    } /* end local RMS boxes */


    /* Set output values */
    VLrow->Ra2000    = pos[0];
    VLrow->Dec2000   = pos[1];
    VLrow->PeakInt   = pbf * MFrow->Peak;
    VLrow->MajorAxis = MFrow->MajorAx;
    VLrow->MinorAxis = MFrow->MinorAx;
    VLrow->PosAngle  = MFrow->PosAngle;
    VLrow->QCenter   = QCen;
    VLrow->UCenter   = UCen;
    VLrow->PFlux     = PFlux;
    VLrow->IRMS      = pbf * IRMS;
    if (QRMS==fblank)  VLrow->PolRMS = fblank;
    else VLrow->PolRMS = pbf * QRMS;
    VLrow->ResRMS    = pbf * MFrow->ResRMS;
    VLrow->ResPeak   = pbf * MFrow->ResPeak;
    VLrow->ResFlux   = pbf * MFrow->ResFlux;
    VLrow->CenterX   = MFrow->PixelCenterX;
    VLrow->CenterY   = MFrow->PixelCenterY;
    VLrow->JDProcess = jd;
    strncpy (VLrow->Field, fieldName, 8);
    
    /* Write row - don't write if bad gain */
    orow = -1;
    if (!badGain) ObitTableVLWriteRow (VLTable, orow, VLrow, err);
    if (err->error) goto cleanup;
  } /* end loop over  MF table */

  /* Close up */
  ObitTableMFClose (MFTable,  err);
  ObitTableVLClose (VLTable,  err);
  
 cleanup:
  
  /* release objects */
  MFrow = ObitTableMFRowUnref(MFrow);
  VLrow = ObitTableVLRowUnref(VLrow);
  IData = ObitFArrayUnref(IData);
  QData = ObitFArrayUnref(QData);
  UData = ObitFArrayUnref(UData);
  QInterp = ObitFInterpolateUnref(QInterp);
  UInterp = ObitFInterpolateUnref(UInterp);
  if (err->error) Obit_traceback_msg (err, routine, MFTable->name);
} /* end ObitTableMF2VL */

/**
 * Print contents of MF table
 * \param in        MF table to print
 * \param prtFile   Where to write
 * \param *err      ObitErr error stack.
 */
void ObitTableMFPrint (ObitTableMF *in, ObitImage *image, FILE  *prtFile, 
		       ObitErr *err)
{
  ObitTableMFRow *row = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, imBLC[7];
  olong irow;
  odouble pos[2], ra, dec;
  ofloat cells, maj, min, pa, pixel[2];
  gchar rast[19], decst[19];
  gchar *routine = "ObitTableMFPrint";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableMFIsA(in), err, 
		       "%s input %s not an MF Table", routine, in->name);
  Obit_return_if_fail (ObitImageIsA(image), err,  
		       "%s image %s not an image", routine, image->name);

  /* Image information */
  cells = sqrt (fabs(image->myDesc->cdelt[0]*image->myDesc->cdelt[1]));
  /* Subimage in image? */
  for (i=0; i<7; i++) imBLC[i] = 1;
  ObitInfoListGetTest(image->info, "BLC", &type, dim, imBLC);
  for (i=0; i<7; i++) imBLC[i] = MAX (1,imBLC[i]) - 1;

  ObitTableMFOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  fprintf (prtFile,"\n Listing of MF table\n");
  fprintf (prtFile,"Sizes in asec, Peak, Flux in mJy, residual values relative to Peak\n");
  fprintf (prtFile,
	   "             RA           Dec          Peak  Fit Maj Fit min   PA    Flux  decon maj  min      PA   res. RMS res Peak    PixX    PixY\n");

  row  = newObitTableMFRow (in);

  /* Loop over table printing */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableMFReadRow (in, irow, row, err);

    pixel[0] = row->PixelCenterX - imBLC[0];
    pixel[1] = row->PixelCenterY - imBLC[1];
    ObitImageDescGetPos(image->myDesc, pixel, pos, err);
    ra  = pos[0];
    dec = pos[1];
    ObitPosLabelUtilRA2HMS (ra, image->myDesc->ctype[0], rast);
    /*dec = image->myDesc->crval[1] + row->DeltaY;*/
    ObitPosLabelUtilDec2DMS (dec, image->myDesc->ctype[1], decst);
    maj = row->D0Major * 3600.0;
    min = row->D0Minor * 3600.0;
    pa  = row->D0PosAngle;

    fprintf (prtFile," %5d %14s %14s %8.2f %7.3f %7.3f %6.1f %8.2f %7.3f %7.3f %6.1f %8.3f %8.3f %7.1f %7.1f\n",
	     irow, rast, decst, row->Peak*1000.0,
	     row->MajorAx*3600.0,row->MinorAx*3600.0, row->PosAngle,
	     row->IFlux*1000.0,  maj, min, pa,
	     row->ResRMS/row->Peak, row->ResPeak/row->Peak,
	     row->PixelCenterX, row->PixelCenterY);
  } /* end loop over table */

  /* Close up */
  ObitTableMFClose (in,  err);
   if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  /* release row object */
  row = ObitTableMFRowUnref(row);
} /* end ObitTableMFPrint */

/*----------------------Private functions ---------------------------*/
/**
 * Routine translated from the AIPSish VDPAVG.FOR/VDPAVG  
 * Determines integrated polarized flux density for a Gaussian  
 * component by doing a weighted sum and then normalizing by the beam  
 * area.  Data bias corected before averaging.
 * Routine translated from the AIPSish VDPAVG.FOR/VDPAVG  
 * \param QData   Q pol array
 * \param UData   U pol array MUST have sane geometry as QData
 * \param pixel   1-rel pixel coordinates
 * \param Gaus    Gaussian (maj (FWHM)< minor, PA) (pix, pix, deg)
 * \param rms     Polarization (Q,U) RMS. 
 * \param barea   Beam area of synthesized beam in pixels 
 * return Integrated polarization in map units*beam, fblank if all blanked
 */
/* Average Polarization over a Gaussian */
static ofloat AverPFlux (ObitFArray *QData, ObitFArray *UData, ofloat pixel[2],
			 ofloat Gaus[3], ofloat PRMS, ofloat beamarea)
{
  ofloat sumpol;
  olong   ix, iy, ixbeg, ixend, iybeg, iyend, nx, ny;
  olong pos[2];
  ofloat *q, *u, gmaj, gmin, pa, fblank = ObitMagicF();
  ofloat xpix, ypix, sum, sumwt, wt, spa2, cpa2, s2pa, xmaj2, xmin2, 
    a, b, c, x, y, arg, box, ppol;
  ofloat con = 2.772589;

  /* Checks */
  g_assert (ObitFArrayIsCompatable(QData, UData));

  /* Input to local */
  xpix = pixel[0] - 1;  /* 0-rel */
  ypix = pixel[1] - 1;
  gmaj = Gaus[0];
  gmin = Gaus[1];
  pa   = Gaus[2];
  nx   = QData->naxis[0];
  ny   = QData->naxis[1];
  sumpol = fblank;

  /* Outside window? */
  if ((xpix < 0.5)  ||  (ypix < 0.5)  ||  
      (xpix > nx-1.5)  ||  (ypix > ny-1.5)) return sumpol;

  /* Consider area major axis  from center */
  box = MAX (gmaj, gmin);
  ixbeg = xpix - box + 0.5;
  ixend = xpix + box + 0.5;
  iybeg = ypix - box + 0.5;
  iyend = ypix + box + 0.5;
  ixbeg = MAX (0, ixbeg);
  iybeg = MAX (0, iybeg);
  ixend = MIN (nx-1, ixend);
  iyend = MIN (ny-1, iyend);
  sum   = 0.0;
  sumwt = 0.0;
  
  /* Set Gaussian parameters */
  spa2 = sin (pa * 1.745329e-2)*sin (pa * DG2RAD);
  cpa2 = cos (pa * 1.745329e-2)*cos (pa *  DG2RAD);
  s2pa = - sin (2.0 * pa *  DG2RAD);
  xmaj2 = gmaj * gmaj / con;
  xmin2 = gmin * gmin / con;
  a = (cpa2 / xmaj2) + (spa2 / xmin2);
  b = (spa2 / xmaj2) + (cpa2 / xmin2);
  c = s2pa * ((1.0 / xmin2) - (1.0/xmaj2));

  /* Loop over region */
  for (iy= iybeg; iy<=iyend; iy++) { /* loop 200 */
    y = (iy - ypix);
    pos[0] = 0; pos[1] = iy;
    q = ObitFArrayIndex (QData, pos);
    u = ObitFArrayIndex (UData, pos);
    for (ix= ixbeg; ix<=ixend; ix++) { /* loop 100 */
      if ((q[ix-1] != fblank)  &&  (u[ix-1] != fblank)) {
	x = (ix - xpix);
	arg = a*x*x + b*y*y + c*x*y;

	/* Don't bother with < 10^-2 */
	if (arg < 4.6) {
	  wt = exp (-arg);
	  /* Get pol. amp. removing bias. */
	  ppol = sqrt (q[ix-1]*q[ix-1] + u[ix-1]*u[ix-1]);

	  /* Correct bias using method from  POLCO. */
	  PBias (ppol, PRMS);

	  /* Weight by square of Gaussian */
	  sum   += ppol * wt;
	  sumwt += wt * wt;
	} 
      } 
    } /* end loop  L100: */
  } /* end loop  L200: */

  /* Enough data? */
  if (sumwt > 0.5) {
    sumpol = (gmaj * gmin * 1.1331 * sum / sumwt) / beamarea;
  } else {
    sumpol = fblank;
  }
  return sumpol;
} /* end AverPFlux */

/**
 * Estimates the polarization bias in a polarization amplitude, P,  
 * measured in the presence of Q and U RMS noise, RMS.  
 * Returns the corrected value.  
 * The bias correction is such that the average bias is removed;  
 * thus the average in the absence of a signal is zero.  Does table  
 * look of values calculated by J. Condon in the range of P/RMS of 1.253  
 * (the no signal limit) and 4 (the high SNR regime).  Does second  
 * order Lagrange interpolation.  At lower values of P/RMS the bias is  
 * a constant 1.253*RMS. Above a signal-to-noise ratio of 4, use the  
 * formula:     
 * normalized bias = 1 / (2 * s) + 1 / (8 * s**3),   
 * where s is the true normalized flux density, iterating once to  
 * estimate s from the normalized map flux density.  "Normalized" means  
 * divided by the rms noise in the q and u maps.  
 * Routine translated from the AIPSish PBIAS.FOR/PDBIAS  
 * \param p     P is the estimated intrinsic total polarized intensity.  
 * \param rms   The standard deviation of the (assumed equal) 
 *              Gaussian distributions of the Stokes Q or U maps. 
 * \return bias corrected value
 */
/* Polarization bias correction */
static ofloat PBias (ofloat p, ofloat rms)
{
  olong   i, index, i1, i2, i3;
  ofloat pnorm, bias, d1, d2, d3, wt1, wt2, wt3, sum, sumwt;
  ofloat table[40][2] = {
    {1.253,1.253},  {1.256,1.156},  {1.266,1.066},  {1.281,0.9814},
    {1.303,0.9030}, {1.330,0.8304}, {1.364,0.7636}, {1.402,0.7023},
    {1.446,0.6462}, {1.495,0.5951}, {1.549,0.5486}, {1.606,0.5064},
    {1.668,0.4683}, {1.734,0.4339}, {1.803,0.4028}, {1.875,0.3749},
    {1.950,0.3498}, {2.027,0.3273}, {2.107,0.3070}, {2.189,0.2888},
    {2.272,0.2724}, {2.358,0.2576}, {2.444,0.2442}, {2.532,0.2321},
    {2.621,0.2212}, {2.711,0.2112}, {2.802,0.2021}, {2.894,0.1938},
    {2.986,0.1861}, {3.079,0.1791}, {3.173,0.1726}, {3.267,0.1666},
    {3.361,0.1610}, {3.456,0.1557}, {3.551,0.1509}, {3.646,0.1463},
    {3.742,0.1420}, {3.838,0.1380}, {3.934,0.1342}, {4.031,0.1306}
  };

  /* Check RMS */
  if (rms <= 0.0) return p;
  pnorm = p / rms;
  
  /* Which regime? */
  if (pnorm <= table[0][0]) {
    /* Low (no) SNR case */
    p -= table[0][1] * rms;
  } else if (pnorm >= table[39][0]) {
    /* High SNR */
    bias = 1.0 / (2.0 * pnorm) + 1.0 / (8.0 * pnorm*pnorm*pnorm);
    pnorm -= bias;
    bias = 1.0 / (2.0 * pnorm) + 1.0 / (8.0 * pnorm*pnorm*pnorm);
    
    /* Correct for bias */
    p -= bias * rms;
  } else {
    /* Middle, interpolate in table */
    index = 1;
    for (i= 3; i<=39; i++) { /* loop 20 */
      if (pnorm < table[i-1][0]) break;
      index = i;
    } /* end loop  L20:  */;
    
    /* Lagrange interpolation */
    i1 = index - 1;
    i2 = index;
    i3 = index + 1;
    d1 = (table[i1-1][0] - table[i2-1][0]) * (table[i1-1][0] - table[i3-1][0]);
    d2 = (table[i2-1][0] - table[i1-1][0]) * (table[i2-1][0] - table[i3-1][0]);
    d3 = (table[i3-1][0] - table[i1-1][0]) * (table[i3-1][0] - table[i2-1][0]);
    wt1 = (pnorm - table[i2-1][0]) * (pnorm - table[i3-1][0]) / d1;
    wt2 = (pnorm - table[i1-1][0]) * (pnorm - table[i3-1][0]) / d2;
    wt3 = (pnorm - table[i1-1][0]) * (pnorm - table[i2-1][0]) / d3;
    sum = table[i1-1][1] * wt1 + table[i2-1][1] * wt2 + table[i3-1][1] * wt3;
    sumwt = wt1 + wt2 + wt3;
    if (sumwt > 0.0) {
      bias = sum / sumwt;
    } else {
      /* Shouldn't ever get here but do  something reasonable. */
      bias = table[i2-1][1];
    } 
    /* Correct for bias */
    p -= bias * rms;
  } /* end interpolate */
  
  return p;
} /* end PBias */
  
/**
 * Determines the RMS in Data in a window of half width RMSsize
 * ObitFArray and copy the values.
 * \param Data    Pixel array
 * \param pixel   Desired center pixel (1-rel)
 * \param RMSsize half width of box around pixel in Data
 * \param err     Obit error stack object.
 * \return the RMS from a histogram analysis of a possibly quantized image, 
 * -1 on failure.
 */
/* RMS in box about pixel */
static ofloat RMSbox (ObitFArray *Data, ofloat pixel[2], olong RMSsize,
		      ObitErr *err)
{
  ofloat out = -1.0;
  olong blc[MAXFARRAYDIM] = {1,1,1,1,1,1,1};
  olong trc[MAXFARRAYDIM] = {0,0,0,0,0,0,0};
  ObitFArray *Box=NULL;
  gchar *routine = "RMSbox";

   /* error checks */
  if (err->error) return out;

  /* Set window (convert to 0-rel)*/
  blc[0] = (olong)(pixel[0]+0.5) - 1 - RMSsize;
  blc[1] = (olong)(pixel[1]+0.5) - 1 - RMSsize;
  trc[0] = (olong)(pixel[0]+0.5) - 1 + RMSsize;
  trc[1] = (olong)(pixel[1]+0.5) - 1 + RMSsize;

  /* Trim at edges */
  blc[0] = MAX (0, blc[0]);
  blc[1] = MAX (0, blc[1]);
  trc[0] = MIN (Data->naxis[0]-1, trc[0]);
  trc[1] = MIN (Data->naxis[1]-1, trc[1]);

  /* Get box */
  Box = ObitFArraySubArr(Data, blc, trc, err);
  if (err->error) {
    Box = ObitFArrayUnref(Box);  /* Cleanup */
    Obit_traceback_val (err, routine, Data->name, out);
  }
  out = ObitFArrayRMSQuant (Box);
  Box = ObitFArrayUnref(Box);  /* Cleanup */
  return out;
} /* end RMSbox */

