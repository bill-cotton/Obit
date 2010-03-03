/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2010                                          */
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

#include "glib/gqsort.h"
#include "ObitTableCCUtil.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableCCUtil.c
 * ObitTableCC class utility function definitions.
 */

/*----------------------Private function prototypes----------------------*/
/** Private: Form sort structure for a table */
static ofloat* 
MakeCCSortStruct (ObitTableCC *in, olong *number, olong *size, olong *ncomp,
		  ofloat *parms, ObitErr *err);

/** Private: Form sort structure for a table with selection by row  */
static ofloat* 
MakeCCSortStructSel (ObitTableCC *in, olong startComp, olong endComp, 
		     olong *size, olong *number, olong *ncomp, ofloat *parms, 
		     ObitErr *err);

/** Private: Sort comparison function for positions */
static gint CCComparePos (gconstpointer in1, gconstpointer in2, 
			  gpointer ncomp);

/** Private: Sort comparison function for Flux density */
static gint CCCompareFlux (gconstpointer in1, gconstpointer in2, 
		     gpointer ncomp);

/** Private: Merge entries in Sort structure */
static void CCMerge (ofloat *base, olong size, olong number); 

/** Private: Merge spectral entries in Sort structure */
static void CCMergeSpec (ofloat *base, olong size, olong number, 
			 gboolean doSpec, gboolean doTSpec); 

/** Private: reorder table based on Sort structure */
static ObitIOCode 
ReWriteTable(ObitTableCC *out, ofloat *base, olong size, olong number, 
	     ofloat *parms, ObitErr *err);
/*----------------------Public functions---------------------------*/

/**
 * Grid components as points onto grid.  The image from which the components is
 * derived is described in desc.  The output grid is padded by a factor OverSample.
 * If the components are Gaussians, their parameters are returned in gaus.
 * \param in         Table to grid
 * \param OverSample Expansion factor for output image
 * \param first      First component (1-rel) to include, 0=>1, filled in if changed
 * \param last       Last component (1-rel) to include, 0=>all, filled in if changed
 * \param noNeg      Ignore first negative flux component and after
 * \param factor     factor to multiply timec fluxes
 * \param minFlux    Minimum abs. value flux density to include (before factor)
 * \param maxFlux    Maximum abs. value flux density to include (before factor)
 * \param desc       Descriptor for image from which components derived
 * \param grid       [out] filled in array, created, resized if necessary
 * \param gparm      [out] Gaussian parameters (major, minor, PA all in deg) if
 *                   the components in in are Gaussians, else, -1.
 * \param ncomp      [out] number of components gridded.
 * \param err        ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableCCUtilGrid (ObitTableCC *in, olong OverSample, 
				olong *first, olong *last, gboolean noNeg,
				ofloat factor, ofloat minFlux, ofloat maxFlux,
				ObitImageDesc *desc, ObitFArray **grid, 
				ofloat gparm[3], olong *ncomp, 
				ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableCCRow *CCRow = NULL;
  olong itestX, itestY;
  ofloat ftestX, ftestY, maxX, minX, maxY, minY;
  ofloat *array, xpoff, ypoff;
  olong j, irow, xPix, yPix, iAddr;
  ofloat iCellX, iCellY, fNx, fNy;
  olong ndim, naxis[2], nx, ny, count = 0, badCnt = 0;
  gchar *routine = "ObitTableCCGrid";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableCCIsA(in));

  gparm[0] = gparm[1] = gparm[2] = -1.0; /* init Gaussian */
  
  /* Create/resize output if necessary */
  ndim = 2;
  naxis[0] = OverSample*desc->inaxes[desc->jlocr];
  naxis[1] = OverSample*desc->inaxes[desc->jlocd];
  /* (re)allocate memory for plane */
  if (*grid!=NULL) *grid = ObitFArrayRealloc(*grid, ndim, naxis);
  else *grid = ObitFArrayCreate("ModelImage", ndim, naxis);

  /* Zero fill */
  ObitFArrayFill (*grid, 0.0);

  /* Get pointer to in->plane data array */
  naxis[0] = 0; naxis[1] = 0; 
  array = ObitFArrayIndex(*grid, naxis);
  
  /* Image size as float */
  nx  = OverSample*desc->inaxes[desc->jlocr];
  ny  = OverSample*desc->inaxes[desc->jlocd];
  fNx = (ofloat)nx;
  fNy = (ofloat)ny;
  /* allowed range of X */
  minX = (-fNx/2.0) * fabs(desc->cdelt[desc->jlocr]);
  maxX = ((fNx/2.0) - 1.0) * fabs(desc->cdelt[desc->jlocr]);
  /* allowed range of Y */
  minY = (-fNy/2.0) * fabs(desc->cdelt[desc->jlocd]);
  maxY = ((fNy/2.0) - 1.0) * fabs(desc->cdelt[desc->jlocd]);

  /* Open CC table */
  retCode = ObitTableCCOpen (in, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);
  
  /* Create table row */
  if (!CCRow) CCRow = newObitTableCCRow (in);

  /* Field specific stuff */
  /* Inverse Cell spacings */
  if (desc->cdelt[desc->jlocr]!=0.0) iCellX = 1.0/desc->cdelt[desc->jlocr];
  else iCellX = 1.0;
  if (desc->cdelt[desc->jlocd]!=0.0) iCellY = 1.0/desc->cdelt[desc->jlocd];
  else iCellY = 1.0;

  /*    Get reference pixel offsets from (nx/2+1, ny/2+1) */
  xpoff = (desc->crpix[desc->jlocr] - (desc->inaxes[desc->jlocr]/2) - 1) *
    desc->cdelt[desc->jlocr];
  ypoff = (desc->crpix[desc->jlocd] - (desc->inaxes[desc->jlocd]/2) - 1) *
    desc->cdelt[desc->jlocd];

  /* loop over CCs */
  count = 0;  /* Count of components */
  if (*first<=0) *first = 1;
  if (*last<=0) *last = in->myDesc->nrow;
  *last = MIN (*last, in->myDesc->nrow);
  for (j=*first; j<=*last; j++) {
    irow = j;
    retCode = ObitTableCCReadRow (in, irow, CCRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, retCode);
    
    /* Get any Gaussian parameters on first, else check */
    if (j==*first) {
      /* Is this a Gaussian component? */
      if ((in->parmsCol>=0) &&
	  (in->myDesc->dim[in->parmsCol][0]>4) && 
	  (CCRow->parms[3]==1.0)) {
	gparm[0] = CCRow->parms[0];
	gparm[1] = CCRow->parms[1];
	gparm[2] = CCRow->parms[2];
      }
    } else if (gparm[0]>0.0) {
      /* All Gaussians MUST be the same */
      if ((CCRow->parms[0]!=gparm[0]) || (CCRow->parms[1]!=gparm[1]) || 
	  (CCRow->parms[2]!=gparm[2])) {
	Obit_log_error(err, OBIT_Error,"%s: All Gaussians MUST have same size",
		       routine);
	return retCode ;
      }
    } /* end of Gaussian Check */

    /* Only to first negative? */
    if (noNeg && (CCRow->Flux<0.0)) break;
  
    /* Process component */
    CCRow->Flux *= factor;     /* Apply factor */
    CCRow->DeltaX += xpoff;    /* Reference pixel offset from  (nx/2,ny/2)*/
    CCRow->DeltaY += ypoff;

    /* Component wanted? larger than in->minFlux and not zero */
    if ((fabs(CCRow->Flux)<minFlux)  || (CCRow->Flux==0.0)) continue;
    /* Nothing too big */
    if (fabs(CCRow->Flux)>maxFlux) continue;
    
    /* Check that comps are on cells */
    ftestX = CCRow->DeltaX * iCellX; /* Convert to cells */
    ftestY = CCRow->DeltaY * iCellY;
    if (ftestX>0.0) itestX = (olong)(ftestX + 0.5);
    else itestX = (olong)(ftestX - 0.5);
    if (ftestY>0.0) itestY = (olong)(ftestY + 0.5);
    else itestY = (olong)(ftestY - 0.5);
  
    /* Count bad cells */
    if ((fabs((ftestX-itestX)>0.1)) || (fabs((ftestY-itestY)>0.1))) {
      badCnt++;
      /* Warn but keep going */
      if (badCnt<50) {
	Obit_log_error(err, OBIT_InfoWarn, "%s Warning: Bad cell %f %f", 
		       routine, CCRow->DeltaX, CCRow->DeltaY);
      }
    }

    /* Clip range of X,Y */
    CCRow->DeltaX = MIN (maxX, MAX (CCRow->DeltaX, minX));
    CCRow->DeltaY = MIN (maxY, MAX (CCRow->DeltaY, minY));

    /* X,Y to cells */
    CCRow->DeltaX *= iCellX;
    CCRow->DeltaY *= iCellY;
    /* 0-rel pixel numbers */
    xPix = (olong)(CCRow->DeltaX + nx/2 + 0.5);
    yPix = (olong)(CCRow->DeltaY + ny/2 + 0.5);

    /* Sum into image */
    iAddr = xPix + nx * yPix;
    array[iAddr] += CCRow->Flux;

    count++;        /* how many */
  } /* end loop over components */

  /* How many? */
  *ncomp = count;
  
  /* Close Table */
  retCode = ObitTableCCClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);
  
  /* Release table/row */
  CCRow   = ObitTableCCRowUnref (CCRow);
  
  return retCode;
} /* end ObitTableCCUtilGrid */

/**
 * Grid spectral components as points onto grid.  The image from which the components is
 * derived is described in desc.  The output grid is padded by a factor OverSample.
 * If the components are Gaussians, their parameters are returned in gaus.
 * Output image is spectral term iterm times the flux density.
 * Works for both parameterized spectra (Param[3} 10-19) 
 * or tabulated spectra (Param[3} 20-29).
 * \param in         Table to grid
 * \param OverSample Expansion factor for output image
 * \param iterm      Spectral term to grid, 0=flux, 1=flux*si, 2=flux*curve...
 *                   For tabulated spectra this is the spectrum
 * \param first      First component (1-rel) to include, 0=>1, filled in if changed
 * \param last       Last component (1-rel) to include, 0=>all, filled in if changed
 * \param noNeg      Ignore first negative flux component and after
 * \param factor     factor to multiply timec fluxes
 * \param minFlux    Minimum abs. value flux density to include (before factor)
 * \param maxFlux    Maximum abs. value flux density to include (before factor)
 * \param desc       Descriptor for image from which components derived
 * \param grid       [out] filled in array, created, resized if necessary
 * \param gparm      [out] Gaussian parameters (major, minor, PA all in deg) if
 *                   the components in in are Gaussians, else, -1.
 * \param ncomp      [out] number of components gridded.
 * \param err        ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableCCUtilGridSpect (ObitTableCC *in, olong OverSample, olong iterm,
				     olong *first, olong *last, gboolean noNeg,
				     ofloat factor, ofloat minFlux, ofloat maxFlux,
				     ObitImageDesc *desc, ObitFArray **grid, 
				     ofloat gparm[3], olong *ncomp, 
				     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableCCRow *CCRow = NULL;
  olong itestX, itestY;
  ofloat ftestX, ftestY, maxX, minX, maxY, minY;
  ofloat *array, xpoff, ypoff;
  olong j, irow, xPix, yPix, iAddr;
  ofloat iCellX, iCellY, fNx, fNy, spectTerm;
  olong ndim, naxis[2], nx, ny, parmoff, count = 0, badCnt = 0;
  gboolean doSpec=TRUE, doTSpec=FALSE;
  gchar *routine = "ObitTableCCGridSpect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableCCIsA(in));

  gparm[0] = gparm[1] = gparm[2] = -1.0; /* init Gaussian */
  
  /* Create/resize output if necessary */
  ndim = 2;
  naxis[0] = OverSample*desc->inaxes[desc->jlocr];
  naxis[1] = OverSample*desc->inaxes[desc->jlocd];
  /* (re)allocate memory for plane */
  if (*grid!=NULL) *grid = ObitFArrayRealloc(*grid, ndim, naxis);
  else *grid = ObitFArrayCreate("ModelImage", ndim, naxis);

  /* Zero fill */
  ObitFArrayFill (*grid, 0.0);

  /* Get pointer to in->plane data array */
  naxis[0] = 0; naxis[1] = 0; 
  array = ObitFArrayIndex(*grid, naxis);
  
  /* Image size as float */
  nx  = OverSample*desc->inaxes[desc->jlocr];
  ny  = OverSample*desc->inaxes[desc->jlocd];
  fNx = (ofloat)nx;
  fNy = (ofloat)ny;
  /* allowed range of X */
  minX = (-fNx/2.0) * fabs(desc->cdelt[desc->jlocr]);
  maxX = ((fNx/2.0) - 1.0) * fabs(desc->cdelt[desc->jlocr]);
  /* allowed range of Y */
  minY = (-fNy/2.0) * fabs(desc->cdelt[desc->jlocd]);
  maxY = ((fNy/2.0) - 1.0) * fabs(desc->cdelt[desc->jlocd]);

  /* Open CC table */
  retCode = ObitTableCCOpen (in, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);
  
  /* Create table row */
  if (!CCRow) CCRow = newObitTableCCRow (in);

  /* Field specific stuff */
  /* Inverse Cell spacings */
  if (desc->cdelt[desc->jlocr]!=0.0) iCellX = 1.0/desc->cdelt[desc->jlocr];
  else iCellX = 1.0;
  if (desc->cdelt[desc->jlocd]!=0.0) iCellY = 1.0/desc->cdelt[desc->jlocd];
  else iCellY = 1.0;

  /*    Get reference pixel offsets from (nx/2+1, ny/2+1) */
  xpoff = (desc->crpix[desc->jlocr] - (desc->inaxes[desc->jlocr]/2) - 1) *
    desc->cdelt[desc->jlocr];
  ypoff = (desc->crpix[desc->jlocd] - (desc->inaxes[desc->jlocd]/2) - 1) *
    desc->cdelt[desc->jlocd];

  /* Where is spectral term? */
  parmoff = 3;

  /* loop over CCs */
  count = 0;  /* Count of components */
  if (*first<=0) *first = 1;
  if (*last<=0) *last = in->myDesc->nrow;
  *last = MIN (*last, in->myDesc->nrow);
  for (j=*first; j<=*last; j++) {
    irow = j;
    retCode = ObitTableCCReadRow (in, irow, CCRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, retCode);
    
    /* Make sure this has a spectrum */
    Obit_retval_if_fail((CCRow->parms && (CCRow->parms[3]>=10.)), err, retCode, 
			"%s: CCs do not contain spectra", routine);
    

    /* Get any Gaussian parameters on first, else check */
    if (j==*first) {
      /* Is this a Gaussian component? */
      if ((in->parmsCol>=0) &&
	  (in->myDesc->dim[in->parmsCol][0]>4) && 
	  ((CCRow->parms[3]==1.0) || (CCRow->parms[3]==11.0) || (CCRow->parms[3]==21.0))) {
	gparm[0] = CCRow->parms[0];
	gparm[1] = CCRow->parms[1];
	gparm[2] = CCRow->parms[2];
      }
    } else if (gparm[0]>0.0) {
      /* All Gaussians MUST be the same */
      if ((CCRow->parms[0]!=gparm[0]) || (CCRow->parms[1]!=gparm[1]) || 
	  (CCRow->parms[2]!=gparm[2])) {
	Obit_log_error(err, OBIT_Error,"%s: All Gaussians MUST have same size",
		       routine);
	return retCode ;
      }
    } /* end of Gaussian Check */

      /* Get spectrum type */
    doSpec  = (CCRow->parms[3]>=9.9)  && (CCRow->parms[3]<=19.0);
    doTSpec = (CCRow->parms[3]>=19.9) && (CCRow->parms[3]<=29.0);
    
    /* Only to first negative? */
    if (noNeg && (CCRow->Flux<0.0)) break;
  
    /* Process component */
    CCRow->Flux *= factor;     /* Apply factor */
    CCRow->DeltaX += xpoff;    /* Reference pixel offset from  (nx/2,ny/2)*/
    CCRow->DeltaY += ypoff;

    /* Component wanted? larger than in->minFlux and not zero */
    if ((fabs(CCRow->Flux)<minFlux)  || (CCRow->Flux==0.0)) continue;
    /* Nothing too big */
    if (fabs(CCRow->Flux)>maxFlux) continue;
    
    /* Check that comps are on cells */
    ftestX = CCRow->DeltaX * iCellX; /* Convert to cells */
    ftestY = CCRow->DeltaY * iCellY;
    if (ftestX>0.0) itestX = (olong)(ftestX + 0.5);
    else itestX = (olong)(ftestX - 0.5);
    if (ftestY>0.0) itestY = (olong)(ftestY + 0.5);
    else itestY = (olong)(ftestY - 0.5);
  
    /* Count bad cells */
    if ((fabs((ftestX-itestX)>0.1)) || (fabs((ftestY-itestY)>0.1))) {
      badCnt++;
      /* Warn but keep going */
      if (badCnt<50) {
	Obit_log_error(err, OBIT_InfoWarn, "%s Warning: Bad cell %f %f", 
		       routine, CCRow->DeltaX, CCRow->DeltaY);
      }
    }

    /* Clip range of X,Y */
    CCRow->DeltaX = MIN (maxX, MAX (CCRow->DeltaX, minX));
    CCRow->DeltaY = MIN (maxY, MAX (CCRow->DeltaY, minY));

    /* X,Y to cells */
    CCRow->DeltaX *= iCellX;
    CCRow->DeltaY *= iCellY;
    /* 0-rel pixel numbers */
    xPix = (olong)(CCRow->DeltaX + nx/2 + 0.5);
    yPix = (olong)(CCRow->DeltaY + ny/2 + 0.5);

    /* Spectral term parameter or tabulated */
    if (doSpec) {
      spectTerm = 1.0;
      if (iterm==1) spectTerm = CCRow->parms[1+parmoff];
      if (iterm==2) spectTerm = CCRow->parms[2+parmoff];
      if (iterm==3) spectTerm = CCRow->parms[3+parmoff];
    } else if (doTSpec) {
      spectTerm = CCRow->parms[iterm+parmoff];
    } else  spectTerm = 0.0;  /* Something went wrong */

    /* Sum into image */
    iAddr = xPix + nx * yPix;
    if (doSpec) {
      array[iAddr] += CCRow->Flux*spectTerm;
    } else if (doTSpec) {
     array[iAddr] += spectTerm;
    }

    count++;        /* how many */
  } /* end loop over components */

  /* How many? */
  *ncomp = count;
  
  /* Close Table */
  retCode = ObitTableCCClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);
  
  /* Release table/row */
  CCRow   = ObitTableCCRowUnref (CCRow);
  
  return retCode;
} /* end ObitTableCCUtilGridSpect */

/**
 * Return an ObitFArray containing the list of components in the CC table
 * from one image which appear in another.
 * Returned array has component values on the first axis and one row per
 * overlapping component.  (X cell (0-rel), Y cell (0-rel), flux).  
 * CCs on cells within 0.5 pixels of outDesc are included.
 * Components in the same cell are combined.
 * If the components are Gaussians, their parameters are returned in gparm.
 * \param in         Table of CCs
 * \param inDesc     Descriptor for image from which components derived
 * \param outDesc    Descriptor for output image 
 * \param grid       [out] filled in array, created, resized if necessary
 * \param gparm      [out] Gaussian parameters (major, minor, PA (all deg)) 
 *                   if the components in in are Gaussians, else, -1.
 *                   These are the values from the first CC.
 * \param ncomp      [out] number of components in output list (generally less 
 *                   than size of FArray).
 * \param err        ObitErr error stack.
 * \return pointer to list of components, may be NULL on failure, 
 *  MUST be Unreffed.
 */
ObitFArray* 
ObitTableCCUtilCrossList (ObitTableCC *inCC, ObitImageDesc *inDesc,  
			  ObitImageDesc *outDesc, ofloat gparm[3], 
			  olong *ncomps, ObitErr *err)
{
  ObitFArray *outArray=NULL;
  ObitIOCode retCode;
  ObitTableCCRow *CCRow = NULL;
  olong i, count, number, nout, irow, lrec;
  olong ncomp;
  olong nrow, ndim, naxis[2], size, fsize, tsize;
  ofloat *table, inPixel[2], outPixel[2], *SortStruct = NULL;
  ofloat *entry;
  gboolean wanted, doSpec=TRUE, doTSpec=FALSE;
  gchar *outName;
  ObitCCCompType modType;
  gchar *routine = "ObitTableCCUtilCrossList";
  
  /* error checks */
  if (err->error) return outArray;

  /* Open */
  retCode = ObitTableCCOpen (inCC, OBIT_IO_ReadOnly, err);
  /* If this fails try ReadWrite */
  if (err->error) { 
    ObitErrClearErr(err);  /* delete failure messages */
    retCode = ObitTableCCOpen (inCC, OBIT_IO_ReadWrite, err);
  }
  if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inCC->name, outArray);

  /* Create sortStruct 
     element size */
  nrow = inCC->myDesc->nrow;
  /* Normal CCs or with spectra? */
  doSpec = (inCC->noParms>4);
  if (doSpec) fsize = 4;
  else fsize = 3;
  size = fsize * sizeof(ofloat);
  /*   Total size of structure in case all rows valid */
  tsize = size * nrow;
  /* create output structure */
  SortStruct = ObitMemAlloc0Name (tsize, "CCSortStructure");

  /* Create table row */
  CCRow = newObitTableCCRow (inCC);

  /* Initialize */
  gparm[0] = gparm[1] = gparm[2] = -1.0;  /* No Gaussian yet */
  count = 0;         /* How many CCs accepted */
  /* If only 3 col, or parmsCol 0 size then this is a point model */
  if ((inCC->myDesc->nfield==3) || 
      (inCC->parmsCol<0) ||
      (inCC->myDesc->dim[inCC->parmsCol]<=0)) 
    modType = OBIT_CC_PointMod;
  else  
    modType = OBIT_CC_Unknown; /* Model type not yet known */

  /* Get spectrum type */
  doSpec  = (CCRow->parms[3]>=9.9)  && (CCRow->parms[3]<=19.0);
  doTSpec = (CCRow->parms[3]>=19.9) && (CCRow->parms[3]<=29.0);
  
  /* Loop over table reading CCs */
  for (i=1; i<=nrow; i++) {

    irow = i;
    retCode = ObitTableCCReadRow (inCC, irow, CCRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

    /* Get model type  */
    if (modType == OBIT_CC_Unknown) {
      modType = (olong)(CCRow->parms[3] + 0.5);
      /* If Gaussian take model */
      if ((modType==OBIT_CC_GaussMod)     || (modType==OBIT_CC_CGaussMod) ||
	  (modType==OBIT_CC_GaussModSpec) || (modType==OBIT_CC_CGaussModSpec)) {
	gparm[0] = CCRow->parms[0];
	gparm[1] = CCRow->parms[1];
	gparm[2] = CCRow->parms[2];
      }
      /* If neither a point nor Gaussian - barf */
      if ((modType!=OBIT_CC_GaussMod) && (modType!=OBIT_CC_CGaussMod) && 
	  (modType!=OBIT_CC_PointMod) && (modType!=OBIT_CC_GaussModSpec) && 
	  (modType!=OBIT_CC_CGaussModSpec) && (modType!=OBIT_CC_PointModSpec) &&
	  (modType!=OBIT_CC_GaussModTSpec) && (modType!=OBIT_CC_CGaussModTSpec) && 
	  (modType!=OBIT_CC_PointModTSpec)) {
	Obit_log_error(err, OBIT_Error,
		       "%s: Model type %d neither point nor Gaussian in %s",
		       routine, modType, inCC->name);
	goto cleanup;
      }
    } /* end model type checking */

    /* Is this one within 3 pixels of outDesc? */
    inPixel[0] = CCRow->DeltaX / inDesc->cdelt[0] + inDesc->crpix[0];
    inPixel[1] = CCRow->DeltaY / inDesc->cdelt[1] + inDesc->crpix[1];
    wanted = ObitImageDescCvtPixel (inDesc, outDesc, inPixel, outPixel, err);
    if (err->error) goto cleanup;

    if (wanted) { /* yes */
      /* add to structure */
      entry = (ofloat*)(SortStruct + count * fsize);  /* set pointer to entry */
      entry[0] = outPixel[0] - 1.0;  /* Make zero rel. pixels */
      entry[1] = outPixel[1] - 1.0;
      entry[2] = CCRow->Flux;
      count++; /* How many? */
    }
  } /* end loop over TableCC */

  /* Close */
  retCode = ObitTableCCClose (inCC, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Release table row */
  CCRow = ObitTableCCRowUnref (CCRow);

  /* Catch anything? */
  *ncomps = count;
  if (count<=0) goto cleanup;
    
  /* Sort */
  number = count; /* Total number of entries */
  ncomp  = 2;     /* number of values to compare */
  g_qsort_with_data (SortStruct, number, size, CCComparePos, &ncomp);

  /* Merge entries - Normal or with spectra? */
  doSpec = (inCC->noParms>4);
  if (doSpec || doTSpec) 
    CCMergeSpec (SortStruct, fsize, number, doSpec, doTSpec);
  else
    CCMerge (SortStruct, fsize, number);
  
  /* Sort to descending merged flux densities */
  ncomp = 1;
  g_qsort_with_data (SortStruct, number, size, CCCompareFlux, &ncomp);

  /* Count number of valid entries left */
  entry = SortStruct;
  count = 0;
  for (i=0; i<number; i++) {
    if (entry[0]>-1.0e19) count++;
    entry += fsize;  /* pointer in table */
  }

  /* Create FArray list large enough for merged CCs */
  ndim = 2; naxis[0] = 3; naxis[1] = count;
  nout = naxis[1];
  lrec = naxis[0];   /* size of table row */
  nout = naxis[1];   /* Size of output array */
  outName =  g_strconcat ("CC List:", inCC->name, NULL);
  outArray = ObitFArrayCreate (outName, ndim, naxis);
  g_free (outName);  /* deallocate name */

  /* Get pointer to array */
  naxis[0] = naxis[1] = 0;
  table = ObitFArrayIndex (outArray, naxis);

  /* Copy structure to output array */
  entry = SortStruct;
  count = 0;
  for (i=0; i<number; i++) {

    /* Deleted? */
    if (entry[0]>-1.0e19) {
      /* Check that array not blown */
      if (count>nout) {
	Obit_log_error(err, OBIT_Error,"%s: Internal array overrun",
		       routine);
	goto cleanup;
      }
      /* copy to out */
      table[0] = entry[0];
      table[1] = entry[1];
      table[2] = entry[2];
      table += lrec;
      count++;
    } /* end of contains value */
    entry += fsize;  /* pointer in table */
  } /* end loop over array */
  
  /* How many? */
  *ncomps = count;

  /* Cleanup */
 cleanup:
  if (SortStruct) ObitMemFree(SortStruct);
  if (err->error) Obit_traceback_val (err, routine, inCC->name, outArray);

  return outArray;
} /*  end ObitTableCCUtilCrossList */

/**
 * Return an ObitFArray containing the list of spectral components in the 
 * CC table from one image which appear in another.
 * Returned array has component values on the first axis and one row per
 * overlapping component.  (X cell (0-rel), Y cell (0-rel), flux, spectral terms).  
 * CCs on cells within 0.5 pixels of outDesc are included.
 * Components in the same cell are combined.
 * If the components are Gaussians, their parameters are returned in gparm.
 * Works for both parameterized spectra (Param[3} 10-19) 
 * or tabulated spectra (Param[3} 20-29).
 * For tabulated spectra, selected channel components are summed and saved as 
 * element 4 in the output array.
 * For paramertized spectra, element 4 is flux weighted spectral term, except for 
 * iterm=0 in which it's the sum of the flux density .
 * \param in         Table of CCs
 * \param inDesc     Descriptor for image from which components derived
 * \param outDesc    Descriptor for output image 
 * \param grid       [out] filled in array, created, resized if necessary
 * \param gparm      [out] Gaussian parameters (major, minor, PA (all deg)) 
 *                   if the components in in are Gaussians, else, -1.
 *                   These are the values from the first CC.
 * \param ncomp      [out] number of components in output list (generally less 
 *                   than size of FArray).
 * \param iterm      Select spectral term, 0=flux, 1=si, 2=curvature...
 *                   For tabulated spectra this is the spectrum
 * \param err        ObitErr error stack.
 * \return pointer to list of components, may be NULL on failure, 
 *  MUST be Unreffed.
 */
ObitFArray* 
ObitTableCCUtilCrossListSpec (ObitTableCC *inCC, ObitImageDesc *inDesc,  
			      ObitImageDesc *outDesc, ofloat gparm[3], 
			      olong *ncomps, olong iterm, ObitErr *err)
{
  ObitFArray *outArray=NULL;
  ObitIOCode retCode;
  ObitTableCCRow *CCRow = NULL;
  olong i, count, number, nout, irow, lrec;
  olong ncomp, parmoff=3;
  olong nrow, ndim, naxis[2], size, fsize, tsize;
  ofloat *table, inPixel[2], outPixel[2], *SortStruct = NULL;
  ofloat *entry, spectTerm;
  gboolean wanted;
  gboolean doSpec=TRUE, doTSpec=FALSE;
  gchar *outName;
  ObitCCCompType modType;
  gchar *routine = "ObitTableCCUtilCrossListSpec";
  
  /* error checks */
  if (err->error) return outArray;

  /* Open */
  retCode = ObitTableCCOpen (inCC, OBIT_IO_ReadOnly, err);
  /* If this fails try ReadWrite */
  if (err->error) { 
    ObitErrClearErr(err);  /* delete failure messages */
    retCode = ObitTableCCOpen (inCC, OBIT_IO_ReadWrite, err);
  }
  if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inCC->name, outArray);

  /* Create sortStruct 
     element size */
  nrow = inCC->myDesc->nrow;
  fsize = 4;
  size = fsize * sizeof(ofloat);
  /*   Total size of structure in case all rows valid */
  tsize = size * nrow;
  /* create output structure */
  SortStruct = ObitMemAlloc0Name (tsize, "CCSortStructure");

  /* Create table row */
  CCRow = newObitTableCCRow (inCC);

  /* Initialize */
  gparm[0] = gparm[1] = gparm[2] = -1.0;  /* No Gaussian yet */
  count = 0;         /* How many CCs accepted */
  /* If only 3 col, or parmsCol 0 size then this is a point model */
  if ((inCC->myDesc->nfield==3) || 
      (inCC->parmsCol<0) ||
      (inCC->myDesc->dim[inCC->parmsCol]<=0)) 
    modType = OBIT_CC_PointMod;
  else  
    modType = OBIT_CC_Unknown; /* Model type not yet known */
  
  /* Loop over table reading CCs */
  for (i=1; i<=nrow; i++) {

    irow = i;
    retCode = ObitTableCCReadRow (inCC, irow, CCRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

    /* Get model type  */
    if (modType == OBIT_CC_Unknown) {
      modType = (olong)(CCRow->parms[3] + 0.5);
      /* Get spectrum type */
      doSpec  = (CCRow->parms[3]>=9.9)  && (CCRow->parms[3]<=19.0);
      doTSpec = (CCRow->parms[3]>=19.9) && (CCRow->parms[3]<=29.0);
      /* If Gaussian take model */
      if ((modType==OBIT_CC_GaussMod)      || (modType==OBIT_CC_CGaussMod) ||
	  (modType==OBIT_CC_GaussModSpec)  || (modType==OBIT_CC_CGaussModSpec) ||
	  (modType==OBIT_CC_GaussModTSpec) || (modType==OBIT_CC_CGaussModTSpec)) {
	gparm[0] = CCRow->parms[0];
	gparm[1] = CCRow->parms[1];
	gparm[2] = CCRow->parms[2];
      }
      /* If neither a point nor Gaussian - barf */
      if ((modType!=OBIT_CC_PointMod)       && (modType!=OBIT_CC_PointModSpec)  &&
	  (modType!=OBIT_CC_PointModTSpec)  &&
	  (modType!=OBIT_CC_CGaussMod)      && (modType!=OBIT_CC_GaussMod)      && 
	  (modType!=OBIT_CC_CGaussModSpec)  && (modType!=OBIT_CC_GaussModSpec)  && 
	  (modType!=OBIT_CC_CGaussModTSpec) && (modType!=OBIT_CC_GaussModTSpec)) {
	Obit_log_error(err, OBIT_Error,
		       "%s: Model type %d neither point nor Gaussian in %s",
		       routine, modType, inCC->name);
	goto cleanup;
      }
    } /* end model type checking */

    /* Is this one within 3 pixels of outDesc? */
    inPixel[0] = CCRow->DeltaX / inDesc->cdelt[0] + inDesc->crpix[0];
    inPixel[1] = CCRow->DeltaY / inDesc->cdelt[1] + inDesc->crpix[1];
    wanted = ObitImageDescCvtPixel (inDesc, outDesc, inPixel, outPixel, err);
    if (err->error) goto cleanup;

    if (wanted) { /* yes */
    /* Spectral term parameter or tabulated */
    if (doSpec) {
      spectTerm = 1.0;
      if (iterm==1) spectTerm = CCRow->parms[1+parmoff];
      if (iterm==2) spectTerm = CCRow->parms[2+parmoff];
      if (iterm==3) spectTerm = CCRow->parms[3+parmoff];
    } else if (doTSpec) {
      spectTerm = CCRow->parms[iterm+parmoff];
    } else  spectTerm = 0.0;  /* Something went wrong */

      /* add to structure */
      entry = (ofloat*)(SortStruct + count * fsize);  /* set pointer to entry */
      entry[0] = outPixel[0] - 1.0;  /* Make zero rel. pixels */
      entry[1] = outPixel[1] - 1.0;
      entry[2] = CCRow->Flux;
      entry[3] = spectTerm;
      count++; /* How many? */
    }
  } /* end loop over TableCC */

  /* Close */
  retCode = ObitTableCCClose (inCC, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Release table row */
  CCRow = ObitTableCCRowUnref (CCRow);

  /* Catch anything? */
  *ncomps = count;
  if (count<=0) goto cleanup;
    
  /* Sort */
  number = count; /* Total number of entries */
  ncomp  = 2;     /* number of values to compare */
  g_qsort_with_data (SortStruct, number, size, CCComparePos, &ncomp);

  /* Merge entries */
  CCMergeSpec (SortStruct, fsize, number, doSpec, doTSpec);
  
  /* Sort to descending merged flux densities */
  ncomp = 1;
  g_qsort_with_data (SortStruct, number, size, CCCompareFlux, &ncomp);

  /* Count number of valid entries left */
  entry = SortStruct;
  count = 0;
  for (i=0; i<number; i++) {
    if (entry[0]>-1.0e19) count++;
    entry += fsize;  /* pointer in table */
  }

  /* Create FArray list large enough for merged CCs */
  ndim = 2; naxis[0] = fsize; naxis[1] = count;
  nout = naxis[1];
  lrec = naxis[0];   /* size of table row */
  nout = naxis[1];   /* Size of output array */
  outName =  g_strconcat ("CC List:", inCC->name, NULL);
  outArray = ObitFArrayCreate (outName, ndim, naxis);
  g_free (outName);  /* deallocate name */

  /* Get pointer to array */
  naxis[0] = naxis[1] = 0;
  table = ObitFArrayIndex (outArray, naxis);

  /* Copy structure to output array */
  entry = SortStruct;
  count = 0;
  for (i=0; i<number; i++) {

    /* Deleted? */
    if (entry[0]>-1.0e19) {
      /* Check that array not blown */
      if (count>nout) {
	Obit_log_error(err, OBIT_Error,"%s: Internal array overrun",
		       routine);
	goto cleanup;
      }
      /* copy to out */
      table[0] = entry[0];
      table[1] = entry[1];
      table[2] = entry[2];
      table[3] = entry[3];
      table += lrec;
      count++;
    } /* end of contains value */
    entry += fsize;  /* pointer in table */
  } /* end loop over array */

  /* Flux weighted flux is sum of flux */
  if (iterm==0) {
    table = ObitFArrayIndex (outArray, naxis);
    for (i=0; i<count; i++) {
      table[3] = table[2];
      table += lrec;
    }
  } /* end of flux only */
  
  /* How many? */
  *ncomps = count;

  /* Cleanup */
 cleanup:
  if (SortStruct) ObitMemFree(SortStruct);
  if (err->error) Obit_traceback_val (err, routine, inCC->name, outArray);

  return outArray;
} /*  end ObitTableCCUtilCrossListSpec */

/**
 * Merge elements of an ObitTableCC on the same position.
 * First sorts table, collapses, sorts to desc. flux
 * \param in      Table to sort
 * \param out     Table to write output to
 * \param err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableCCUtilMerge (ObitTableCC *in, ObitTableCC *out, 
				 ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong size, fsize, number=0, ncomp;
  ofloat parms[20];
  ofloat *SortStruct = NULL;
  gboolean doSpec=TRUE, doTSpec=FALSE;
  gchar *routine = "ObitTableCCUtilMerge";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableCCIsA(in));

  /* Open table */
  retCode = ObitTableCCOpen (in, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);
  
  /* Must be something in the table else just return */
  if (in->myDesc->nrow<=0) {
    retCode = ObitTableCCClose (in, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);
    return OBIT_IO_OK;
  }

  /* build sort structure from table */
  SortStruct = MakeCCSortStruct (in, &size, &number, &ncomp, parms, err);
  if (err->error) goto cleanup;

  /* Close table */
  retCode = ObitTableCCClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Sort */
  g_qsort_with_data (SortStruct, number, size, CCComparePos, &ncomp);

  /* Get spectrum type */
  doSpec  = (parms[3]>=9.9)  && (parms[3]<=19.0);
  doTSpec = (parms[3]>=19.9) && (parms[3]<=29.0);

  /* Merge entries - Normal or with spectra? */
  fsize = size/sizeof(ofloat);
  if (doSpec || doTSpec) 
    CCMergeSpec (SortStruct, fsize, number, doSpec, doTSpec);
  else
    CCMerge (SortStruct, fsize, number);
  
  /* Sort to descending merged flux densities */
  ncomp = 1;
  g_qsort_with_data (SortStruct, number, size, CCCompareFlux, &ncomp);

  /* Clone output table from input */
  out = (ObitTableCC*)ObitTableClone ((ObitTable*)in, (ObitTable*)out);

  /* Write output table */
  retCode = ReWriteTable (out, SortStruct, fsize, number, parms, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Cleanup */
 cleanup:
  if (SortStruct) ObitMemFree(SortStruct);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  return OBIT_IO_OK;
} /* end ObitTableCCUtilMerge */

/**
 * Merge elements of an ObitTableCC on the same position.
 * with selection by row number.
 * First sorts table, collapses, sorts to desc. flux
 * \param in        Table to sort
 * \param startComp First component to select 
 * \param endComp   Last component to select, 0=> all
 * \param parms     [out] Non-point parameters (MUST all be the same ) dimen. at least 4.
 *                  parms[3] = type
 *                  0 => point, no parameters
 *                  1 = Sky Gaussian, [0:3]=maj, min, PA
 *                  2 = Convolved Gaussian, [0:3]=maj axis, min axis, PA (all deg)
 *                  3 = Uniform Sphere [0] = radius (deg)
 * \param err       ObitErr error stack.
 * \return FArray containing merged CC table contents; MUST be Unreffed.
 *                Will contain flux, X, Y, + any spectral terms
 * \li Flux
 * \li Delta X
 * \li Delta Y
 */
ObitFArray* ObitTableCCUtilMergeSel (ObitTableCC *in, olong startComp, 
				     olong endComp, ofloat *parms, 
				     ObitErr *err)
{
  ObitFArray *out = NULL;
  olong i, j, count, lout, nout, ndim, naxis[2];
  ObitIOCode retCode;
  olong size, fsize, number=0, ncomp, nterms;
  ofloat lparms[20];
  ofloat *entry, *outArray, *SortStruct = NULL;
  gboolean doSpec=TRUE, doTSpec=FALSE;
  gchar *routine = "ObitTableCCUtilMergeSel";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitTableCCIsA(in));

  /* Open table */
  retCode = ObitTableCCOpen (in, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Must be something in the table else just return */
  if (in->myDesc->nrow<=0) {
    retCode = ObitTableCCClose (in, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    return OBIT_IO_OK;
  }

  /* Check range of components */
  startComp = MAX (1, startComp);
  endComp   = MIN (in->myDesc->nrow, endComp);

  /* build sort structure from table */
  SortStruct = MakeCCSortStructSel (in, startComp, endComp, 
				    &size, &number, &ncomp, lparms, 
				    err);
  if (err->error) goto cleanup;

  /* Close table */
  retCode = ObitTableCCClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Sort */
  g_qsort_with_data (SortStruct, number, size, CCComparePos, &ncomp);

  /* Get spectrum type */
  doSpec  = (lparms[3]>=9.9)  && (lparms[3]<=19.0);
  doTSpec = (lparms[3]>=19.9) && (lparms[3]<=29.0);

  /* Merge entries */
  fsize = size/sizeof(ofloat);
  if (doSpec || doTSpec) 
    CCMergeSpec (SortStruct, fsize, number, doSpec, doTSpec);
  else
    CCMerge (SortStruct, fsize, number);
  
  /* Sort to descending merged flux densities */
  ncomp = 1;
  g_qsort_with_data (SortStruct, number, size, CCCompareFlux, &ncomp);

  /* Count number of valid entries left */
  entry = SortStruct;
  count = 0;
  for (i=0; i<number; i++) {
    if (entry[0]>-1.0e19) count++;
    entry += fsize;  /* pointer in table */
  }


  /* Create output Array */
  lout = 3;
  /* Need room for spectral terms? */
  if (in->noParms>4)  nterms = in->noParms-4;
  else nterms = 0;
  if (in->noParms>4) lout += nterms;

  ndim = 2; naxis[0] = lout; naxis[1] = count;
  nout = count;
  out      = ObitFArrayCreate ("MergedCC", ndim, naxis);
  naxis[0] = naxis[1] = 0;
  outArray = ObitFArrayIndex(out, naxis);

  /* Any model parameters */
  for (i=0; i<4; i++) parms[i] = 0.0;
  for (i=0; i<MIN (4,in->noParms); i++) parms[i] = lparms[i];

  /* Copy structure to output array */
  entry = SortStruct;
  count = 0;
  for (i=0; i<number; i++) {

    /* Deleted? */
    if (entry[0]>-1.0e19) {
      /* Check that array not blown */
      if (count>nout) {
	Obit_log_error(err, OBIT_Error,"%s: Internal array overrun",
		       routine);
	goto cleanup;
      }
      /* copy to out */
      outArray[0] = entry[2];
      outArray[1] = entry[0];
      outArray[2] = entry[1];
      for (j=0; j<nterms; j++) outArray[3+j] = entry[3+j];
      outArray += lout;
      count++;
    } /* end of contains value */
    entry += fsize;  /* pointer in table */
  } /* end loop over array */
  
  /* Cleanup */
 cleanup:
  if (SortStruct) ObitMemFree(SortStruct);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  return out;
} /* end ObitTableCCUtilMergeSel */


/**
 * Merge spectral elements of an ObitTableCC on the same position.
 * with selection by row number.
 * Spectral components are flux weighted average
 * First sorts table, collapses, sorts to desc. flux
 * \param in        Table to sort
 * \param startComp First component to select 
 * \param endComp   Last component to select, 0=> all
 * \param parms     [out] Non-point parameters (MUST all be the same ) dimen. at least 4.
 *                  parms[3] = type
 *                  0 => point, no parameters
 *                  1 = Sky Gaussian, [0:3]=maj, min, PA
 *                  2 = Convolved Gaussian, [0:3]=maj axis, min axis, PA (all deg)
 *                  3 = Uniform Sphere [0] = radius (deg)
 * \param err       ObitErr error stack.
 * \return FArray containing merged CC table contents; MUST be Unreffed.
 *                Will contain flux, X, Y, + any spectral terms
 * \li Flux
 * \li Delta X
 * \li Delta Y
 * \li Parms if [3]>=10 [4-?] are spectral components
 */
ObitFArray* ObitTableCCUtilMergeSelSpec (ObitTableCC *in, olong startComp, 
					  olong endComp, ofloat *parms, 
					  ObitErr *err)
{
  ObitFArray *out = NULL;
  olong i, j, count, lout, nout, ndim, naxis[2];
  ObitIOCode retCode;
  olong size, fsize, number=0, ncomp, nterms;
  ofloat lparms[20];
  gboolean doSpec=TRUE, doTSpec=FALSE;
  ofloat *entry, *outArray, *SortStruct = NULL;
  gchar *routine = "ObitTableCCUtilMergeSelSpec";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitTableCCIsA(in));

  /* Open table */
  retCode = ObitTableCCOpen (in, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Must be something in the table else just return */
  if (in->myDesc->nrow<=0) {
    retCode = ObitTableCCClose (in, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    return OBIT_IO_OK;
  }

  /* Check range of components */
  startComp = MAX (1, startComp);
  endComp   = MIN (in->myDesc->nrow, endComp);

  /* build sort structure from table */
  SortStruct = MakeCCSortStructSel (in, startComp, endComp, 
				    &size, &number, &ncomp, lparms, 
				    err);
  if (err->error) goto cleanup;

  /* Close table */
  retCode = ObitTableCCClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Sort */
  g_qsort_with_data (SortStruct, number, size, CCComparePos, &ncomp);

  /* Get spectrum type */
  doSpec  = (lparms[3]>=9.9)  && (lparms[3]<=19.0);
  doTSpec = (lparms[3]>=19.9) && (lparms[3]<=29.0);

  /* Merge entries */
  fsize = size/sizeof(ofloat);
  CCMergeSpec (SortStruct, fsize, number, doSpec, doTSpec);
  
  /* Sort to descending merged flux densities */
  ncomp = 1;
  g_qsort_with_data (SortStruct, number, size, CCCompareFlux, &ncomp);

  /* Count number of valid entries left */
  entry = SortStruct;
  count = 0;
  for (i=0; i<number; i++) {
    if (entry[0]>-1.0e19) count++;
    entry += fsize;  /* pointer in table */
  }

  /* Create output Array */
  lout = 3;
  /* Need room for spectral terms? */
  if (in->noParms>4)  nterms = in->noParms-4;
  else nterms = 0;
  if (in->noParms>4) lout += nterms;

  ndim = 2; naxis[0] = lout; naxis[1] = count;
  nout = count;
  out      = ObitFArrayCreate ("MergedCC", ndim, naxis);
  naxis[0] = naxis[1] = 0;
  outArray = ObitFArrayIndex(out, naxis);

  /* Any model parameters */
  for (i=0; i<4; i++) parms[i] = 0.0;
  for (i=0; i<MIN (4,in->noParms); i++) parms[i] = lparms[i];

  /* Copy structure to output array */
  entry = SortStruct;
  count = 0;
  for (i=0; i<number; i++) {

    /* Deleted? */
    if (entry[0]>-1.0e19) {
      /* Check that array not blown */
      if (count>nout) {
	Obit_log_error(err, OBIT_Error,"%s: Internal array overrun",
		       routine);
	goto cleanup;
      }
      /* copy to out */
      outArray[0] = entry[2];
      outArray[1] = entry[0];
      outArray[2] = entry[1];
      for (j=0; j<nterms; j++) outArray[3+j] = entry[3+j];
      outArray += lout;
      count++;
    } /* end of contains value */
    entry += fsize;  /* pointer in table */
  } /* end loop over array */
  
  /* Cleanup */
 cleanup:
  if (SortStruct) ObitMemFree(SortStruct);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  return out;
} /* end ObitTableCCUtilMergeSelSpec */


/**
 * Merge selected entries in a CC Table then write values with 
 * absolute values of the flux densities into an output table
 * \param image     input image with input CC table
 * \param inCCver   input CC table
 * \param outCCver  Desired output CC table on image, if 0 then new
 *                  value used returned.
 * \param startComp First component to select 
 * \param endComp   Last component to select, 0=> all
 * \param range     Max and min abs value of flux densities to accept
 * \param err       ObitErr error stack.
 */
ObitTableCC* 
ObitTableCCUtilMergeSel2Tab (ObitImage *image, olong inCCver, olong *outCCver,
			     olong startComp, olong endComp, 
			     ofloat range[2], ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitFArray *Merged=NULL;
  ObitTable *tempTable=NULL;
  ObitTableCC *inCCTable = NULL, *outCCTable = NULL;
  ObitTableCCRow *CCRow = NULL;
  gchar *tabType = "AIPS CC";
  olong naxis[2], warray, larray, i, j, ver, orow, offset, nSpec=0;
  ofloat *array, parms[20];
  gboolean doSpec=FALSE, doTSpec=FALSE;
  gchar *routine = "bitTableCCUtilMergeSel2Tab";

   /* error checks */
  if (err->error) return outCCTable;
  g_assert (ObitImageIsA(image));

  /* Get input CC table */
  ver = inCCver;
  tempTable = newObitImageTable (image,OBIT_IO_ReadOnly, tabType, &ver, err);
  if ((tempTable==NULL) || (err->error)) 
     Obit_traceback_val (err, routine, image->name, outCCTable);
  inCCTable = ObitTableCCConvert(tempTable);
  tempTable = ObitTableUnref(tempTable);
  if (err->error) Obit_traceback_val (err, routine, image->name, outCCTable);
  
  /* Open input */
  retCode = ObitTableCCOpen (inCCTable, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, image->name, outCCTable);

  /* Select and merge input table */
  Merged = ObitTableCCUtilMergeSel (inCCTable, startComp, endComp, parms, err);
  if (err->error) goto cleanup;
  naxis[0] = 0; naxis[1]=0; 
  array = ObitFArrayIndex(Merged, naxis);
  warray = Merged->naxis[0];
  larray = Merged->naxis[1];
  
  /* Get spectrum type */
  if (inCCTable->noParms>4) {
    doSpec  = (parms[3]>=9.9)  && (parms[3]<=19.0);
    doTSpec = (parms[3]>=19.9) && (parms[3]<=29.0);
    nSpec = inCCTable->noParms - 4;
  }
  
  /* Create output CC table */
  ver = *outCCver;
  outCCTable = newObitTableCCValue ("Merged/selected", (ObitData*)image,
				    &ver, OBIT_IO_WriteOnly, inCCTable->noParms, 
				    err);
  if (err->error) goto cleanup;
  *outCCver = ver;  /* save if defaulted (0) */
      
  /* Open output */
  retCode = ObitTableCCOpen (outCCTable, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Create table row */
  CCRow = newObitTableCCRow (outCCTable);

  /* Copy Parms to row */
  for (j=0; j<MIN(4,inCCTable->noParms); j++) CCRow->parms[j] = parms[j];
  
  offset = 4;
  /* loop over table */
  for (j=1; j<=larray; j++) {
    /* Want this one? */
    if ((fabs(array[0])>=range[0]) && (fabs(array[0])<=range[1])) {
      /* Copy row data */
      CCRow->Flux   = array[0];
      CCRow->DeltaX = array[1];
      CCRow->DeltaY = array[2];
      /* Copy any spectra */
      if (doSpec) {
	for (i=0; i<nSpec; i++) CCRow->parms[i+offset] = array[3+i];
      } else if (doTSpec) {
	for (i=0; i<nSpec; i++) CCRow->parms[i+offset] = array[3+i];
      }
      
      /* Write output */
      orow = -1;
      retCode = ObitTableCCWriteRow (outCCTable, orow, CCRow, err);
    }
    if  (err->error) goto cleanup;
    
    /* Update */
    array += warray;
  } /* end loop over table */

  /* Close/cleanup */
 cleanup:
  retCode   = ObitTableCCClose (outCCTable, err);
  inCCTable = ObitTableUnref(inCCTable);
  CCRow     = ObitTableRowUnref(CCRow);
  Merged    = ObitFArrayUnref(Merged);
  if  (err->error) Obit_traceback_val (err, routine, image->name, outCCTable);

 return outCCTable;
} /* end ObitTableCCUtilMergeSel2Tab */

/**
 * Scale flux densities of  elements of an ObitTableCC.
 * \param in        Table to scale
 * \param startComp First component to select 
 * \param endComp   Last component to select, 0=> all
 * \param scale     Factor to multiply times flux densities
 * \param err       ObitErr error stack.
 */
void ObitTableCCUtilScale (ObitTableCC *in, olong startComp, 
			   olong endComp, ofloat scale, ObitErr *err)
{
  ObitTableCCRow *row=NULL;
  olong irow;
  ObitIOCode retCode;
  gchar *routine = "ObitTableCCUtilScale";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitTableCCIsA(in));

  /* Open table */
  retCode = ObitTableCCOpen (in, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Must be something in the table else just return */
  if (in->myDesc->nrow<=0) {
    retCode = ObitTableCCClose (in, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    return;
  }

  /* Check range of components */
  startComp = MAX (1, startComp);
  endComp   = MIN (in->myDesc->nrow, endComp);
  if (endComp<=0) endComp = in->myDesc->nrow;

  if (err->error) goto cleanup;

  /* Create table row */
  row = newObitTableCCRow (in);

  /* Loop over table */
  for (irow=startComp; irow<=endComp; irow++) {
    retCode = ObitTableCCReadRow (in, irow, row, err);
    if (row->status<0) continue;  /* Skip deselected record */
    if (err->error) goto cleanup;
    row->Flux *= scale;
    retCode = ObitTableCCWriteRow (in, irow, row, err);
    if (err->error) goto cleanup;
  } /* End loop over table */

  /* Close table */
  retCode = ObitTableCCClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Cleanup */
 cleanup:
  row = ObitTableCCRowUnref (row); /* Release table row */
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitTableCCUtilScale */

/**
 * Append the selected components of one table onto another
 * \param inCC      Table to copy from
 * \param outCC     Table to copy to
 * \param startComp First component to select , 0=>1
 * \param endComp   Last component to select, 0=> all
 * \param err       ObitErr error stack.
 */
void ObitTableCCUtilAppend  (ObitTableCC *inCC, ObitTableCC *outCC, 
			     olong startComp, olong endComp, ObitErr *err)
{
  ObitTableCCRow *row=NULL;
  olong irow, orow;
  ObitIOCode retCode;
  gchar *routine = "ObitTableCCUtilAppend";

  /* error checks */
  if (err->error) return;
  g_assert (ObitTableCCIsA(inCC));
  g_assert (ObitTableCCIsA(outCC));

  /* Open input table */
  retCode = ObitTableCCOpen (inCC, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Must be something in the table else just return */
  if (inCC->myDesc->nrow<=0) {
    retCode = ObitTableCCClose (inCC, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    return;
  }

  /* Open output table */
  retCode = ObitTableCCOpen (outCC, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Check range of components */
  startComp = MAX (1, startComp);
  endComp   = MIN (inCC->myDesc->nrow, endComp);
  if (endComp<=0) endComp = inCC->myDesc->nrow;

  /* Create table row */
  row = newObitTableCCRow (outCC);
  ObitTableCCSetRow (outCC, row, err);
  if (err->error) goto cleanup;

  /* Loop over table */
  for (irow=startComp; irow<=endComp; irow++) {
    retCode = ObitTableCCReadRow (inCC, irow, row, err);
    if (row->status<0) continue;  /* Skip deselected record */
    if (err->error) goto cleanup;
    orow = -1;
    retCode = ObitTableCCWriteRow (outCC, orow, row, err);
    if (err->error) goto cleanup;
  } /* End loop over table */

  /* Close tables */
  retCode = ObitTableCCClose (inCC, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  retCode = ObitTableCCClose (outCC, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Cleanup */
 cleanup:
  row = ObitTableCCRowUnref (row); /* Release table row */
  if (err->error) Obit_traceback_msg (err, routine, inCC->name);
} /* end ObitTableCCUtilAppend */

/**
 * Filter out weak, isolated components in a table 
 * Zeroes any component for which the sum of the flux densities of all
 * components within radius is less than minFlux.
 * \param CCTab     CC table object to filter. 
 * \param radius    Radius within which to consider components. (deg) 
 * \param minFlux   Minimum acceptable summed flux
 * \param err       Error/message stack
 * \return True if any components zeroed
 */
gboolean ObitTableCCUtilFiltCC (ObitTableCC *CCTab, ofloat radius, ofloat minFlux, 
				ObitErr* err)  
{
  gboolean out = FALSE;
  ObitTableCCRow *CCRow = NULL;
  olong   irow, nrow, ncc, i, j;
  ofloat mxdis, sum, dis;
  ofloat *xpos=NULL, *ypos=NULL, *flux=NULL;
  olong  *row=NULL;
  gboolean checkMore;
  gchar *routine = "ObitTableCCUtilFiltCC";

  /* Error checks */
  if (err->error) return out;  /* previous error? */

  /* Open table */
  ObitTableCCOpen (CCTab, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;

  /* allocate storage */
  nrow = CCTab->myDesc->nrow;
  if (nrow<1) goto closeit;
  xpos = g_malloc0(nrow*sizeof(ofloat));
  ypos = g_malloc0(nrow*sizeof(ofloat));
  flux = g_malloc0(nrow*sizeof(ofloat));
  row  = g_malloc0(nrow*sizeof(olong));

  /* Copy table to arrays */
  CCRow = newObitTableCCRow (CCTab);
  ncc = 0;
  nrow = CCTab->myDesc->nrow;
  for (irow=1; irow<=nrow; irow++) {
    ObitTableCCReadRow (CCTab, irow, CCRow, err);
    if (err->error) goto cleanup;
    if (CCRow->status==-1) continue;  /* Deselected? */
    row[ncc]    = irow;
    xpos[ncc]   = CCRow->DeltaX;
    ypos[ncc]   = CCRow->DeltaY;
    flux[ncc++] = CCRow->Flux;
  } /* end loop reading table */

  /* Loop to convergence */
  mxdis = radius * radius;
  checkMore = TRUE;
  while (checkMore) {
    /* Sum fluxes within radius */
    checkMore = FALSE;
    for (i=0; i<ncc; i++) { 
      if (flux[i]==0.0) continue;    /* need to check? */
      sum = 0.0;
      for (j=0; j<ncc; j++) {
	if (flux[j]==0.0) continue;  /* need to check? */
	dis = (xpos[i]-xpos[j])*(xpos[i]-xpos[j]) + (ypos[i]-ypos[j])*(ypos[i]-ypos[j]);
	if (dis <= mxdis) sum += flux[j];
      } /* end loop  summing */
      if (sum < minFlux)  {  /* If too small, replace with zero flux */
	out       = TRUE;
	checkMore = TRUE;
	flux[i]   = 0.0;
	ObitTableCCReadRow (CCTab, row[i], CCRow, err);
	CCRow->Flux = 0.0;
	ObitTableCCWriteRow (CCTab, row[i], CCRow, err);
      } 
    } /* end loop over table */
  } /* End iteration */

  /* If out, mark as unsorted */
  if (out) {
    CCTab->myDesc->sort[0] = 0;
    CCTab->myDesc->sort[1] = 0;
  }

  /* Close table */
 closeit:
  ObitTableCCClose (CCTab, err);
  if (err->error) goto cleanup;
      
 cleanup:
   /* deallocate storage */
  CCRow = ObitTableCCRowUnref(CCRow);  
  if (xpos) g_free(xpos); 
  if (ypos) g_free(ypos); 
  if (flux) g_free(flux); 
  if (row)  g_free(row); 
  if (err->error) Obit_traceback_val (err, routine, CCTab->name, out);
  return out;
} /* end of routine ObitTableCCUtilFiltCC */ 

/*----------------------Private functions---------------------------*/
/**
 * Create/fill sort structure for a CC table
 * The sort structure has one "entry" per row which contains 
 * \li Delta X
 * \li Delta Y
 * \li Flux
 * \li [optional] spectral terms +
 *
 * Each valid row in the table has an entry.
 * \param in     Table to sort, assumed already open;
 * \param size   [out] Number of bytes in entry
 * \param number [out] Number of entries
 * \param ncomp  [out] Number of values to compare
 * \param parms  [out] Parms of first element if they exist
 * \param err     ObitErr error stack.
 * \return sort structure, should be ObitMemFreeed when done.
 */
static ofloat* 
MakeCCSortStruct (ObitTableCC *in, olong *size, olong *number, olong *ncomp,
		  ofloat *parms, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ofloat *out = NULL;
  ObitTableCCRow *row = NULL;
  ofloat *entry;
  olong irow, nrow, tsize, count, i, j;
  olong nterms, fsize;
  gchar *routine = "MakeCCSortStruct";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitTableCCIsA(in));

  /* Get table info */
  nrow = in->myDesc->nrow;

  /* element size */
  fsize = 3;
  /* Need room for spectral terms? */
  if (in->noParms>4)  nterms = in->noParms-4;
  else nterms = 0;
  fsize += nterms;
  *size = fsize * sizeof(ofloat);

  /* Total size of structure in case all rows valid */
  tsize = (*size) * (nrow);
  /* create output structure */
  out = ObitMemAlloc0Name (tsize, "CCSortStructure");
  
  /* Compare 2  (X, Y pos) */
  *ncomp = 2;

  /* Create table row */
  row = newObitTableCCRow (in);

 /* loop over table */
  irow = 0;
  count = 0;
  retCode = OBIT_IO_OK;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableCCReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, out);
    if (row->status<0) continue;  /* Skip deselected record */

    /* add to structure */
    entry = (ofloat*)(out + count * fsize);  /* set pointer to entry */
    entry[0] = row->DeltaX;
    entry[1] = row->DeltaY;
    entry[2] = row->Flux;
    /* First 4 parms are model, following are spectral parameters */
    for (j=0; j<nterms; j++) entry[3+j] = row->parms[4+j];

    /* Save parms if any for first record */
    if ((count<=0) && (in->noParms>0)) {
      for (i=0; i<MIN (4,in->noParms); i++) parms[i] = row->parms[i];
    }
    
    count++;  /* How many valid */
  } /* end loop over file */
  
  /* check for errors */
  if ((retCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);
  
  /* Release table row */
  row = ObitTableCCRowUnref (row);
  
  /* Actual number */
  *number = count;

  return out;
} /* end MakeCCSortStruc */ 

/**
 * Create/fill sort structure for a CC table selecting by row
 * The sort structure has one "entry" per row which contains 
 * \li Delta X
 * \li Delta Y
 * \li Delta Flux
 *
 * Each valid row in the table has an entry.
 * \param in        Table to sort, assumed already open;
 * \param startComp First component to select 
 * \param endComp   Last component to select, 0=> all
 * \param size      [out] Number of bytes in entry
 * \param number    [out] Number of entries
 * \param ncomp     [out] Number of values to compare
 * \param parms     [out] Parms of first element if they exist
 * \param err        ObitErr error stack.
 * \return sort structure, should be ObitMemFreeed when done.
 */
static ofloat* 
MakeCCSortStructSel (ObitTableCC *in, olong startComp, olong endComp, 
		     olong *size, olong *number, olong *ncomp, ofloat *parms, 
		     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ofloat *out = NULL;
  ObitTableCCRow *row = NULL;
  ofloat *entry;
  olong irow, nrow, tsize, count, i, j;
  olong nterms, fsize;
  gchar *routine = "MakeCCSortStruct";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitTableCCIsA(in));

  /* Get table info */
  nrow = in->myDesc->nrow;

  /* element size */
  fsize = 3;
  /* Need room for spectral terms? */
  if (in->noParms>4)  nterms = in->noParms-4;
  else nterms = 0;
  fsize += nterms;
  *size = fsize * sizeof(ofloat);

  /* Total size of structure in case all rows valid */
  tsize = (*size) * (nrow);
  /* create output structure */
  out = ObitMemAlloc0Name (tsize, "CCSortStructure");
  
  /* Compare 2  (X, Y pos) */
  *ncomp = 2;

  /* Create table row */
  row = newObitTableCCRow (in);

 /* loop over table */
  irow = startComp-1;
  count = 0;
  retCode = OBIT_IO_OK;
  while ((irow<endComp) && (retCode==OBIT_IO_OK)) {
    irow++;
    retCode = ObitTableCCReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, out);
    if (row->status<0) continue;  /* Skip deselected record */

    /* add to structure */
    entry = (ofloat*)(out + count * fsize);  /* set pointer to entry */
    entry[0] = row->DeltaX;
    entry[1] = row->DeltaY;
    entry[2] = row->Flux;
    /* First 4 parms are model, following are spectral parameters */
    for (j=0; j<nterms; j++) entry[3+j] = row->parms[4+j];

    /* Save 1st 4 parms if any for first record */
    if ((count<=0) && (in->noParms>0)) {
      for (i=0; i<MIN (4,in->noParms); i++) parms[i] = row->parms[i];
    }
    
    count++;  /* How many valid */
  } /* end loop over file */
  
  /* check for errors */
  if ((retCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);
  
  /* Release table row */
  row = ObitTableCCRowUnref (row);
  
  /* Actual number */
  *number = count;

  return out;
} /* end MakeSortStrucSel */ 

/**
 * Compare two lists of floats
 * Conformant to function type GCompareDataFunc
 * \param in1   First list
 * \param in2   Second list
 * \param ncomp Number of values to compare (2)
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CCComparePos (gconstpointer in1, gconstpointer in2, 
			  gpointer ncomp)
{
  gint out = 0;
  ofloat *float1, *float2;

  /* get correctly typed local values */
  float1 = (float*)(in1);
  float2 = (float*)(in2);

  if (float1[0]<float2[0])      out = -1;
  else if (float1[0]>float2[0]) out = 1;
  else                          out = 0;
  if (!out) { /* compare second needed? */
    if (float1[1]<float2[1])      out = -1;
    else if (float1[1]>float2[1]) out = 1;
    else                          out = 0;
  }

  return out;
} /* end CCComparePos */


/**
 * Compare fluxes, to give descending abs order.
 * Conformant to function type GCompareDataFunc
 * \param in1   First list
 * \param in2   Second list
 * \param ncomp Number of values to compare (1)
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CCCompareFlux (gconstpointer in1, gconstpointer in2, 
			   gpointer ncomp)
{
  gint out = 0;
  ofloat *float1, *float2;

  /* get correctly typed local values */
  float1 = (float*)(in1 + 2*sizeof(ofloat));
  float2 = (float*)(in2 + 2*sizeof(ofloat));
  if (fabs(*float1)<fabs(*float2))      out =  1;
  else if (fabs(*float1)>fabs(*float2)) out = -1;
  else                          out =  0;
  return out;
} /* end CCCompareFlux */

/**
 * Merge entries in sort structure
 * leaves "X" entry in defunct rows -1.0e20
 * table and then copies over the input table.
 * \param base    Base address of sort structure
 * \param size    Size in gfloats of a sort element
 * \param number  Number of sort elements
 */
static void CCMerge (ofloat *base, olong size, olong number)
{
  olong i, j;
  ofloat *array = base;
  
  i = 0;
  while (i<number) {
    j=i+1;
    while (j<number) {
      if ((array[j*size]!=array[i*size]) || (array[j*size+1]!=array[i*size+1]))
	break;
      array[i*size+2] += array[j*size+2];  /* sum fluxes */
      array[j*size] = -1.0e20;             /* Don't need any more */
      j++;
    } /* end finding matches */
    i = j;   /* move on */
  } /* end loop over table */

} /* end CCMerge */

/**
 * Merge Spectral entries in sort structure
 * leaves "X" posn entry in defunct rows -1.0e20
 * table and then copies over the input table.
 * For parameterized spectra:
 * Takes flux weighted average of spectral components,
 * assumed to be entries 3+
 * For tabulated spectra:
 * Takes sums spectral components assumed to be entries 3+
 * \param base    Base address of sort structure
 * \param size    Size in gfloats of an element
 * \param number  Number of sort elements
 * \param doSpec  TRUE if parameterized spectra
 * \param doTSpec TRUE if tabulated spectra
 */
static void CCMergeSpec (ofloat *base, olong size, olong number, 
			 gboolean doSpec, gboolean doTSpec)
{
  olong i, j, k;
  ofloat *array = base;
  
  /* Multiply parameterized spectral terms by flux */
  if (doSpec) {
    j = 0;
    while (j<number) {
      for (k=3; k<size; k++) 
	array[j*size+k] *=  array[j*size+2];
      j++;
    }
  }

  i = 0;
  while (i<number) {
    j=i+1;
    while (j<number) {
      if ((array[j*size]!=array[i*size]) || (array[j*size+1]!=array[i*size+1]))
	break;
      /* Sum spectral components  flux */
      for (k=3; k<size; k++) 
	array[i*size+k] += array[j*size+k]; 
      array[i*size+2] += array[j*size+2];  /* sum fluxes */
      array[j*size] = -1.0e20;             /* Don't need any more */
      j++;
    } /* end finding matches */
    i = j;   /* move on */
  } /* end loop over table */

  /* Normalize parameterized spectra by sum of flux */
  if (doSpec) {
    i = 0;
    while (i<number) {
      if ((array[i*size]>-1.0e-19) && (fabs(array[i*size+3])>0.0)) {
	for (k=3; k<size; k++) 
	  array[i*size+k] /= array[i*size+2]; 
      }
      i++;
    }
  }

} /* end CCMergeSpec */

/**
 * Write valid entries in sort structure
 * \param out     Table write
 * \param base    Base address of sort structure
 * \param size    Size in floats of a sort element
 * \param number  Number of sort elements
 * \param parms   Parms of components
 * \param err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
static ObitIOCode 
ReWriteTable(ObitTableCC *out, ofloat *base, olong size, olong number, 
	     ofloat *parms, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableCCRow *row = NULL;
  ofloat *entry;
  olong irow, i, nterms, count;
  gchar *routine = "ReWriteTable";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableCCIsA(out));

  /* Open table */
  retCode = ObitTableCCOpen (out, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, out->name, retCode);

  /* Mark as sorted by descending flux */
  out->myDesc->sort[0] = -(out->FluxCol+257);
  out->myDesc->sort[1] = 0;

  /* Fooey!  This one counts */
 ((ObitTableDesc*)out->myIO->myDesc)->sort[0] = -(out->FluxCol+257);
  ((ObitTableDesc*)out->myIO->myDesc)->sort[1] = 0;
  
  /* Need to copy any  spectral terms? */
  if (out->noParms>4) nterms = out->noParms-4;
  else nterms = 0;

  /* Create row structure */
  row = newObitTableCCRow (out);

  /* loop over table */
  retCode = OBIT_IO_OK;
  entry = (ofloat*)base;
  irow = 0;
  count = 0;
  while (count<number) {

    /* Deleted? */
    if (entry[0]>-1.0e19) {

      /* copy to row */
      row->DeltaX = entry[0];
      row->DeltaY = entry[1];
      row->Flux   = entry[2];
            /* copy any model parms - only one of these */
      if (out->noParms>0) {
	for (i=0; i<MAX(4, out->noParms); i++) row->parms[i] = parms[i];
      }
      /* any spectral terms - one per entry */
      for (i=0; i<nterms; i++) row->parms[4+i] = entry[3+i];

      /* Write */
      irow++;
      retCode = ObitTableCCWriteRow (out, irow, row, err);
      if ((retCode != OBIT_IO_OK) || (err->error)) 
	Obit_traceback_val (err, routine, out->name, retCode);
    }
    count++;
    entry += size;  /* pointer in table */
  } /* end loop over file */
  
  /* Close table */
  retCode = ObitTableCCClose (out, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, out->name, retCode);

  /* Release table row */
  row = ObitTableCCRowUnref (row);

  /* Tell what you've done */
  Obit_log_error(err, OBIT_InfoErr,
		 "Merged %d CC components into %d for %s",
		 number, irow, out->name);

  return retCode;
} /* end ReWriteTable */

