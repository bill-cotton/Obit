/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2007-2011                                          */
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

#include "ObitUVPeelUtil.h"
#include "ObitUVUtil.h"
#include "ObitUVImager.h"
#include "ObitUVImagerIon.h"
#include "ObitUVImagerMF.h"
#include "ObitImageMosaic.h"
#include "ObitSkyGeom.h"
#include "ObitSkyModel.h"
#include "ObitSkyModelVM.h"
#include "ObitSkyModelVMSquint.h"
#include "ObitSkyModelVMIon.h"
#include "ObitUVSelfCal.h"
#include "ObitTableCCUtil.h"
#include "ObitTableSNUtil.h"
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVPeelUtil.c
 * ObitUVPeelUtil Utility function definitions.
 */

/*---------------Private function prototypes----------------*/
void Convert2Dto3D (ObitImage *image, ofloat ccShift[2], ObitErr* err);
void MoveBack (ObitTableCC *inCC, ObitTableCC *outCC, ofloat ccShift[2], 
	       ObitErr* err);
 
/*----------------------Public functions---------------------------*/
/**
 * Loops over sources to be peeled.
 * This routine should only be run after all self-calibration and 
 * autocentering has been done, at least one previous CLEAN is required.
 * Each loop picks the strongest field with peak above PeelFlux and 
 * subtracts all others using the SkyModel on myClean and then 
 * self-calibrates the residual data to obtain a best model and 
 * calibration for that field.
 * A model uv data set derived from the self-calibrated model and 
 * inverse of the self cal calibration is used to corrupt the model 
 * data which is then permanently subtracted from the  input data.
 * The imaging is then redone using the myClean setup.
 * The components peeled are written to CC Table 2 on the output peeled image.
 * After all peeling is done, the components on the CC tables 2 are 
 * appended to CC tables 1 so that these tables contain all components 
 * subtracted.
 * On return, ncomp contains the number of components in each CC table 1
 * which were NOT peeled, i.e. should be subtracted from inUV to produce a 
 * residual data set.
 * No component restoration or flattening is done.
 * \param myInput InfoList with control parameters (most have defaults):
 *  \li PeelFlux      f Minimum level for Peeling (peal in CCs)
 *  \li PeelLoop      i max. number of self cal loops
 *  \li PeelSolInt    f Peel SC Solution interval (min)
 *  \li PeelType      s Peel SC Solution type '  ', 'L1'
 *  \li PeelMode      s Peel SC Solution mode:'A&P', 'P', 'P!A',
 *  \li PeelNiter     i Niter for peel.CLEAN
 *  \li PeelMinFlux   f Minimum Peel Clean component (Jy)
 *  \li PeelRefAnt    i Peel SC Reference antenna
 *  \li PeelSNRMin    f Min. allowed SNR in peel selfcal
 *  \li PeelAvgPol    b Avg. poln in peel self cal?
 *  \li PeelAvgIF     b Avg. IFs in peel self cal?
 *  \li PBCor         b Apply Frequency PB Corr?
 *  \li antSize       f Diameter of ant. for PBCor (m)
 *  \li Robust        f Robustness power
 *  \li nuGrid        i Size in u of weighting grid 
 *  \li nvGrid        i Size in v of weighting grid 
 *  \li WtBox         i Additional rows and columns in weighting
 *  \li WtFunc        i Box function type when WtBox
 *  \li UVTaper       f [2] (U,V) Gaussian taper klambda
 *  \li WtPower       f Power to raise weights to
 *  \li MaxBaseline   f maximum baseline length in wavelengths.
 *  \li MinBaseline   f minimum baseline length in wavelengths.
 *  \li rotate        f rotation of images
 *  \li xCells        f Image cell spacing in X in asec.
 *  \li yCells        f Image cell spacing in Y in asec.
 *  \li Gain          f CLEAN loop gain
 *  \li minPatch      i Min. BEAM half-width.
 *  \li autoWindow    b If true, automatically set windows
 *  \li WtUV          f Weighting to use outside of basic uv range
 *  \li minNo         i Min. allowed no. antennas in selfcal
 *  \li doSmoo        b If true interpolate failed solutions
 *  \li prtLv         i Print level in selfcal, 0=>none
 * \param inUV    Data to be peeled, on return all SN tables will be removed.
 *                and peeled source will have been subtracted
 * \param myClean Clean object which has previously been CLEANed and which 
 *                has a field with a CC peak in excess of PeelFlux.
 *                Need for peeling based on values of myClean->peakFlux and
 *                PeelFlux.  The test for individual fields is the maximum as 
 *                determined from the CLEAN components by 
 *                #ObitImageMosaicMaxField
 *                On output, this will include info on last CLEAN
 * \param  nfield [out] Number of entries in ncomp
 * \param  ncomp  [out] Array of number of components in data (i.e. after any peel)
 *                -1 => none
 *                per field,  this array should be g_freeed after use
 * \param err     Error/message stack
 */
void ObitUVPeelUtilLoop (ObitInfoList* myInput, ObitUV* inUV, 
			 ObitDConCleanVis *myClean, 
			 olong *nfield, olong **ncomp,  ObitErr* err)
{
  ObitTableCC  *peelCCTable=NULL, *outCCTable=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  oint         noParms;
  olong        ver, peelfield, *peeled, i;
  ofloat       PeelFlux;
  gchar *routine = "ObitUVPeelUtilLoop";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inUV));

  /* Initialize output - In case, free old ncomp */
  if (*ncomp!=NULL) g_free(*ncomp); *ncomp = NULL;
  *nfield = myClean->mosaic->numberImages;  /* Keep track of components not peeled */
  *ncomp = g_malloc0((*nfield)*sizeof(olong));

  /* Set prtLv on err */
  ObitInfoListGetTest(myInput, "prtLv", &type, dim, &err->prtLv); 

  /* Peeling trip level */
  PeelFlux = 1.0e20;
  ObitInfoListGetTest(myInput, "PeelFlux", &type, dim, &PeelFlux); 

  /* Anything to do ? */
  if (myClean->peakFlux<=PeelFlux) return;

  /* To keep track of peeled fields */
  peeled = g_malloc0(myClean->mosaic->numberImages*sizeof(olong));
  for (i=0; i< myClean->mosaic->numberImages; i++) peeled[i] = 0;

  /* Loop while sources still brighter than PeelFlux */
  while (myClean->peakFlux>PeelFlux) {
    peelfield = ObitUVPeelUtilPeel (myInput, inUV, myClean, peeled, err);
    if (err->error) Obit_traceback_msg (err, routine, myClean->name);
    if (peelfield>0) {
      peeled[peelfield-1] = 1;  /* keep track of peeled fields */
    } else break;  
  }

  /* Tell why quit */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Brightest facet after peeling %g Jy/bm", myClean->peakFlux);
  ObitErrLog(err); 

  /* Loop appending components in CC Table 2 of peeled fields to Table 1 */
  for (i=0; i<myClean->mosaic->numberImages; i++) {
    /* Keep track of unpeeled components */
    (*ncomp)[i] = myClean->skyModel->endComp[i];
    /* If zero use -1 instead */
    if ((*ncomp)[i]==0) (*ncomp)[i] = -1;

    if (peeled[i]<=0) continue;  /* ingnore fields not peeled */
    /* Copy CCs from peeled table (CC 2)  to CC 1 on output image */
    ver = 2;
    noParms = 0;
    peelCCTable = newObitTableCCValue ("Peeled CC", (ObitData*)myClean->mosaic->images[i],
				       &ver, OBIT_IO_ReadOnly, noParms, 
				       err);
    /* Make sure  created */
    if (peelCCTable==NULL) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: No CC table 2 found on %s", routine, myClean->mosaic->images[i]->name);
      return;
    }
    ver = 1;
    outCCTable = newObitTableCCValue ("outCC", (ObitData*)myClean->mosaic->images[i],
				      &ver, OBIT_IO_ReadWrite, peelCCTable->noParms, 
				      err);
    /* Make sure  object created */
    if (outCCTable==NULL) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: No CC table  found on %s", routine, 
		     myClean->mosaic->images[i]->name);
      return;
    }

    /* Append CCs */
    ObitTableCCUtilAppend (peelCCTable, outCCTable, 1, 0, err);
    if (err->error) Obit_traceback_msg (err, routine, myClean->name);
    peelCCTable = ObitTableCCUnref(peelCCTable);
    outCCTable  = ObitTableCCUnref(outCCTable);
    
  } /* end loop copying peeled CCs */
  if (peeled) g_free(peeled); peeled = NULL;  /* Done with array */

} /* end ObitUVPeelUtilLoop */

/**
 * Peel a strong source from a data set based on previous CLEAN.
 * Picks the strongest field with peak above PeelFlux and subtracts all
 * others using the SkyModel on myClean and then self calibrates the 
 * residual data to obtain a best model and calibration for that field.
 * A model uv data set derived from the selfcalibrated model and 
 * inverse of the self cal calibration is used to corrupt the model 
 * data which is then subtracted from the  input data.
 * The imaging is then redone using the myClean setup .
 * The components peeled are written to CC Table 2 on the output Peeled image.
 * Note: If multiple peels are done, the components from previously peeled 
 * fields will not be in Table 1 and the components from any CC table 2 need 
 * to be copied to table 1 when the peeling is finished.
 * No restoration or flattening is done.
 * \param myInput InfoList with control parameters (most have defaults):
 *  \li PeelFlux      f Minimum level for Peeling (peal in CCs)
 *  \li PeelLoop      i max. number of self cal loops
 *  \li PeelSolInt    f Peel SC Solution interval (min)
 *  \li PeelType      s Peel SC Solution type '  ', 'L1'
 *  \li PeelMode      s Peel SC Solution mode:'A&P', 'P', 'P!A',
 *  \li PeelNiter     i Niter for peel.CLEAN
 *  \li PeelMinFlux   f Minimum Peel Clean component (Jy)
 *  \li PeelRefAnt    i Peel SC Reference antenna
 *  \li PeelSNRMin    f Min. allowed SNR in peel selfcal
 *  \li PeelAvgPol    b Avg. poln in peel self cal?
 *  \li PeelAvgIF     b Avg. IFs in peel self cal?
 *  \li PBCor         b Apply Frequency PB Corr?
 *  \li antSize       f Diameter of ant. for PBCor (m)
 *  \li Robust        f Robustness power
 *  \li nuGrid        i Size in u of weighting grid 
 *  \li nvGrid        i Size in v of weighting grid 
 *  \li WtBox         i Additional rows and columns in weighting
 *  \li WtFunc        i Box function type when WtBox
 *  \li UVTaper       f [2] (U,V) Gaussian taper klambda
 *  \li WtPower       f Power to raise weights to
 *  \li MaxBaseline   f maximum baseline length in wavelengths.
 *  \li MinBaseline   f minimum baseline length in wavelengths.
 *  \li rotate        f rotation of images
 *  \li xCells        f Image cell spacing in X in asec.
 *  \li yCells        f Image cell spacing in Y in asec.
 *  \li Gain          f CLEAN loop gain
 *  \li minPatch      i Min. BEAM half-width.
 *  \li autoWindow    b If true, automatically set windows
 *  \li WtUV          f Weighting to use outside of basic uv range
 *  \li minNo         i Min. allowed no. antennas in selfcal
 *  \li doSmoo        b If true interpolate failed solutions
 *  \li prtLv         i Print level in selfcal, 0=>none
 * \param inUV    Data to be peeled, on return all SN tables will be removed.
 *                and peeled source will have been subtracted
 * \param myClean Clean object which has previously been CLEANed and which 
 *                has a field with a CC peak in excess of PeelFlux.
 *                Need for peeling based on values of myClean->peakFlux and
 *                PeelFlux.  The test for individual fields is the maximum as 
 *                determined from the CLEAN components by 
 *                #ObitImageMosaicMaxField
 *                On output, this will include info on last CLEAN
 * \param donePeel List of flags for fields already peeled, 1 in posn -> peeled
 * \param err     Error/message stack
 * \return the 1-rel field number of the peeled field or -1 if no peel
 */
olong ObitUVPeelUtilPeel (ObitInfoList* myInput, ObitUV* inUV, 
			  ObitDConCleanVis *myClean, olong *donePeel, ObitErr* err)
{
  olong             peeled=-1;
  ObitInfoType      type;
  const ObitSkyModelClassInfo *skyModelClass;
  const ObitUVImagerClassInfo *imagerClass;
  ObitUVImager      *tmpImager=NULL; 
  ObitImageMosaic   *tmpMosaic=NULL; 
  ObitUV            *scrUV=NULL, *tmpUV=NULL;
  ObitSkyModel      *tmpSkyModel=NULL; 
  ObitDConCleanVis  *tmpClean=NULL; 
  ObitUVSelfCal     *selfCal = NULL;
  ObitTableCC       *peelCCTable=NULL, *outCCTable=NULL, *CCTab=NULL;
  ObitTableSN       *inSNTable=NULL, *outSNTable=NULL, *SNInver=NULL;
  olong        i,  j, k, jtemp, peelField, *bcomp=NULL, *ecomp=NULL, *ignore=NULL;
  olong        dft, iter, MaxSCLoop=1, nx, ny, win[4];
  oint         noParms, numPol, numIF;
  olong        ver;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat       PeelFlux, ftemp, xCells, yCells, ccShift[2];
  gboolean     converged, didSC=FALSE, Fl=FALSE, Tr=TRUE, init, noSCNeed, btemp; 
  gchar        stemp[5], *pmode = "P   ";
  gchar        *imgParms[] = {  /* Imaging, weighting parameters */
    "PBCor", "antSize", 
    "Robust", "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "WtPower",
    "MaxBaseline", "MinBaseline", "rotate", "Beam",  "xCells", "yCells", 
    "nx", "ny",  "nxBeam", "nyBeam", "do3D", /*"dispURL",*/
    NULL
  };
  gchar        *CLEANParms[] = {  /* Clean parameters */
    "Gain", "minFlux", "Niter", "minPatch", "Beam", 
    "Mode", "CCFilter", "maxPixel", "autoWindow",
    NULL
  };
  gchar        *SCParms[] = {  /* Self cal parameters */
    "RefAnt", "WtUV", "avgPol", "avgIF", "doMGM", "SNRmin", 
    "minNo", "doSmoo", "modelFlux", "modelPos", "modelParm",
    /* "dispURL", "prtLv", */
    NULL
  };
  gchar        *peelParms[] = {  /* peeling parameters */
    "PBCor", "antSize", 
    NULL
  };
  gchar solmod[5], soltyp[5];
  gchar *routine = "ObitUVPeelUtilPeel";

  /* error checks */
  if (err->error) return peeled;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inUV));

  /* Peeling trip level */
  PeelFlux = 1.0e20;
  ObitInfoListGetTest(myInput, "PeelFlux", &type, dim, &PeelFlux);
 
  /* Need to peel bright sources? */
  if (myClean->peakFlux>PeelFlux) {

    /* Compress CC files */
    ObitSkyModelCompressCC (myClean->skyModel, err);
    if (err->error) Obit_traceback_val (err, routine, myClean->name, peeled);

    /* Make list of fields already peeled to ignore */
    ignore = g_malloc0((myClean->mosaic->numberImages+1)*sizeof(olong));
    k = 0;
    for (j=0; j<=myClean->mosaic->numberImages; j++) ignore[j] = 0;
    for (j=0; j<myClean->mosaic->numberImages; j++) {
      if (donePeel[j]>0) ignore[k++] = j+1;  /* 1-rel */
    }

    /* Get ImageMosaic for brightest field if any exceeds PeelFlux */
    tmpMosaic = ObitImageMosaicMaxField (myClean->mosaic, PeelFlux, 
					 ignore, &peelField, err);
    if (ignore) g_free(ignore); ignore = NULL;
    if (err->error) Obit_traceback_val (err, routine, myClean->name, peeled);
    if (tmpMosaic==NULL) goto DonePeel;

    /* Convert 2D model to 3D if necessary */
    Convert2Dto3D (tmpMosaic->images[0], ccShift, err);
    if (err->error) goto cleanup;

    Obit_log_error(err, OBIT_InfoErr, 
		   " ******  Peeling strong source from field %d", peelField);
    ObitErrLog(err); 

    /* Subtract all model except for peelField */
    bcomp = g_malloc0(myClean->nfield*sizeof(olong));
    ecomp = g_malloc0(myClean->nfield*sizeof(olong));
    for (i=0; i<myClean->nfield; i++) {
      if (i==(peelField-1)) { /* Don't subtract the peel field */
	bcomp[i] = 2;
	ecomp[i] = 1;
      } else if (myClean->mosaic->isShift[i]==peelField) {	
	/* or the shifted counterpart of 2D autoCentered field */
	bcomp[i] = 2;
	ecomp[i] = 1;
      } else {  /* All of other fields */
	bcomp[i] = 1;
	ecomp[i] = 0;
      }
    }
    dim[0] = myClean->nfield; dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myClean->skyModel->info, "BComp", OBIT_long, dim, bcomp);
    ObitInfoListAlwaysPut (myClean->skyModel->info, "EComp", OBIT_long, dim, ecomp);
    if (bcomp) g_free(bcomp); bcomp = NULL;
    if (ecomp) g_free(ecomp); ecomp = NULL;

    /* Need all stokes */
    dim[0] = 4;
    stemp[0] = stemp[1] = stemp[2] = stemp[3] = ' ';  stemp[4] = 0;
    ObitInfoListAlwaysPut (inUV->info, "Stokes", OBIT_string, dim, stemp);
    ObitUVOpen (inUV, OBIT_IO_ReadCal, err);
    ObitUVClose (inUV, err);
    if (err->error) goto cleanup;

    /* Scratch file */
    scrUV = newObitUVScratch (inUV, err);
    /* Give more sensible name */
    if (scrUV->name) g_free(scrUV->name);
    scrUV->name = g_strdup("Peel data");

    /* Subtract */
    Obit_log_error(err, OBIT_InfoErr, " ******  Subtract non-peel sources from uv data");
    ObitErrLog(err); 
    ObitSkyModelSubUV (myClean->skyModel, inUV, scrUV, err);
    if (err->error) goto cleanup;
    /* If no data in output file simply copy - this may occur if there are
       no components in the non peel fields */
    if (scrUV->myDesc->nvis<=0) {
      scrUV  = ObitUVCopy (inUV, scrUV, err);
      if (err->error) goto cleanup;
    }

    /* Copy imaging control info */
    ObitInfoListCopyList (myInput, scrUV->info, imgParms);
    dim[0] = dim[1] = dim[2] = 1;
    xCells = fabs (tmpMosaic->xCells)*3600.0;
    yCells = fabs (tmpMosaic->yCells)*3600.0;
    ObitInfoListAlwaysPut (scrUV->info, "xCells",     OBIT_long, dim, &xCells);
    ObitInfoListAlwaysPut (scrUV->info, "yCells",     OBIT_long, dim, &yCells);
    ObitInfoListAlwaysPut (scrUV->info, "nx",         OBIT_long, dim, tmpMosaic->nx);
    ObitInfoListAlwaysPut (scrUV->info, "ny",         OBIT_long, dim, tmpMosaic->ny);

    /* Temporary Imager*/
    /* Trap Ion variation */
    if (ObitUVImagerIonIsA(myClean->imager)) {
      imagerClass = (const ObitUVImagerClassInfo*)ObitUVImagerGetClass();
    } else {
      imagerClass = (const ObitUVImagerClassInfo*)myClean->imager->ClassInfo;
    }
    if(ObitUVImagerMFIsA(myClean->imager)) {
      /* Needs extra parameter */
      tmpImager   = (ObitUVImager*)ObitUVImagerMFCreate2("Peel imager", 
							 ((ObitImageMosaicMF*)myClean->mosaic)->maxOrder, 
							 scrUV, (ObitImageMosaicMF*)tmpMosaic, err);
    } else {
      tmpImager   = imagerClass->ObitUVImagerCreate2("Peel imager", scrUV, tmpMosaic, err);
    }
    if (err->error) goto cleanup;

    /* Trap Ion variation  SkyModel */
    if (ObitSkyModelVMIonIsA(myClean->skyModel)) {
      skyModelClass = (const ObitSkyModelClassInfo*)ObitSkyModelGetClass();
    /* Trap VMSquint variations 
       I don't know why this seemed needed */
    } else if (ObitSkyModelVMSquintIsA(myClean->skyModel)) {
      skyModelClass = (const ObitSkyModelClassInfo*)ObitSkyModelGetClass();
    } else {
      skyModelClass = (const ObitSkyModelClassInfo*)myClean->skyModel->ClassInfo;
    }
    tmpSkyModel   = skyModelClass->ObitSkyModelCreate("Peel SkyModel", tmpMosaic);
    if (err->error) goto cleanup;

    /* Use DFT model */
    dim[0] = dim[1] = 1;
    dft = (olong)OBIT_SkyModel_DFT;
    ObitInfoListAlwaysPut (tmpSkyModel->info, "Mode", OBIT_long, dim, &dft);

    /* Make temp CleanVis */
    tmpClean = ObitDConCleanVisCreate2("Peel Clean Object", scrUV, 
				      (ObitUVImager*)tmpImager, 
				      (ObitSkyModel*)tmpSkyModel, err);
    if (err->error) goto cleanup;
    /* Share display with myClean */
    tmpClean->display = ObitDisplayRef(myClean->display);

    /* Get input parameters from myInput, copy to tmpClean */
    ObitInfoListCopyList (myInput, tmpClean->info, CLEANParms);
    /* Maximum number of iterations */
    jtemp = 200; dim[0] = dim[1]= dim[2] = 1;
    ObitInfoListGetTest (myInput, "PeelNiter", &type, dim, &jtemp);
    ObitInfoListAlwaysPut(tmpClean->info, "Niter", OBIT_long, dim, &jtemp);
    /* Min Peel image flux density */
    ftemp = 1.0; dim[0] = dim[1]= dim[2] = 1;
    ObitInfoListGetTest (myInput, "PeelMinFlux", &type, dim, &ftemp);
    ObitInfoListAlwaysPut(tmpClean->info, "minFlux", OBIT_float, dim, &ftemp);

    /* No restore or flatten on peeled image */
    dim[0] = dim[1] = 1;
    ObitInfoListAlwaysPut(tmpClean->info, "doRestore", OBIT_bool, dim, &Fl);
    ObitInfoListAlwaysPut(tmpClean->info, "doFlatten", OBIT_bool, dim, &Fl);
    /* Explicitly do weighting  */
    ObitInfoListAlwaysPut(tmpClean->info, "doWeight", OBIT_bool, dim, &Tr);

    if (tmpClean->CCver<=0) tmpClean->CCver = 2;  /* Specify CC ver */
 
    /* Use DFT model */
    dim[0] = dim[1] = 1;
    dft = (olong)OBIT_SkyModel_DFT;
    ObitInfoListAlwaysPut (tmpClean->info, "Mode", OBIT_long, dim, &dft);
    
    /* Create self cal object */
    selfCal = ObitUVSelfCalCreate ("SelfCal", tmpSkyModel);

    /* Set Parameters  */   
    strncpy (solmod, "P   ",4); dim[0] = 4; dim[1] = dim[2] = 1;
    ObitInfoListGetTest (myInput, "PeelMode", &type, dim, solmod);
    ObitInfoListAlwaysPut(selfCal->info, "solMode", OBIT_string, dim, solmod);
    strncpy (soltyp, "    ",4); dim[0] = 4; dim[1] = dim[2] = 1;
    ObitInfoListGetTest (myInput, "PeelType", &type, dim, soltyp);
    ObitInfoListAlwaysPut(selfCal->info, "solType", OBIT_string, dim, soltyp);
    dim[0] = 1;
    /* Solution interval */
    ftemp = 1.0; dim[0] = dim[1]= dim[2] = 1;
    ObitInfoListGetTest (myInput, "PeelSolInt", &type, dim, &ftemp);
    ObitInfoListAlwaysPut(selfCal->info, "solInt", OBIT_float, dim, &ftemp);
    /* Min acceptable peak for Self cal */
    ftemp = PeelFlux*0.5;dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListPut (selfCal->info, "minFluxPSC", OBIT_float, dim, &ftemp, err);
    ObitInfoListPut (selfCal->info, "minFluxASC", OBIT_float, dim, &ftemp, err);
    /* Minimum SNR */
    ftemp = 3.0; dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListGetTest (myInput, "PeelSNRMin", &type, dim, &ftemp);
    ObitInfoListAlwaysPut(selfCal->info, "SNRMin", OBIT_float, dim, &ftemp);
    /* Average IF? */
    btemp = FALSE; dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListGetTest (myInput, "PeelAvgIF", &type, dim, &btemp);
    ObitInfoListAlwaysPut(selfCal->info, "avgIF", OBIT_bool, dim, &btemp);
    /* Average Poln? */
    btemp = FALSE; dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListGetTest (myInput, "PeelAvgPol", &type, dim, &btemp);
    ObitInfoListAlwaysPut(selfCal->info, "avgPol", OBIT_bool, dim, &btemp);
    /* reference antenna? */
    jtemp = i; dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListGetTest (myInput, "PeelRefAnt", &type, dim, &jtemp);
    ObitInfoListAlwaysPut(selfCal->info, "refAnt", OBIT_long, dim, &jtemp);
    /* Always smooth */
    ObitInfoListAlwaysPut(selfCal->info, "doSmoo", OBIT_bool, dim, &Tr);
    /* Copy control info */
    ObitInfoListCopyList (myInput, selfCal->info, SCParms);

    /* Add outer window */
    nx = tmpMosaic->images[0]->myDesc->inaxes[0]; 
    ny = tmpMosaic->images[0]->myDesc->inaxes[1]; 
    win[0] = (nx/2)-10; win[1] = nx/2; win[2] = ny/2; 
    ObitDConCleanWindowOuter (tmpClean->window, 1, OBIT_DConCleanWindow_round,
			      win, err);
    if (err->error) goto cleanup;

    /* Link myClean display if given */
    if (myClean->display) selfCal->display = ObitDisplayRef(myClean->display);
    
    /* Set vis vs baseline histogram */
    ObitUVSelfCalFluxHist(selfCal, scrUV, err);
    if (err->error) goto cleanup;
    
    /* How many loops? */
    MaxSCLoop=1;
    ObitInfoListGetTest (myInput, "PeelLoop", &type, dim, &MaxSCLoop);

    /* Self cal loop */
    converged = FALSE;
    init = TRUE;
    for (iter=0; iter<MaxSCLoop; iter++) {

      Obit_log_error(err, OBIT_InfoErr, " ******  Peel Self Calibration number %d", iter+1);
      ObitErrLog(err); 

      /* Use all components on peel field */
      jtemp = 1; dim[0] = dim[1] = dim[2] = 1;
      ObitInfoListAlwaysPut (myClean->skyModel->info, "BComp", OBIT_long, dim, &jtemp);
      jtemp = myClean->Pixels->iterField[peelField-1];
      ObitInfoListAlwaysPut (myClean->skyModel->info, "EComp", OBIT_long, dim, &jtemp);
      
      /* Phase only except for last iteration */
      dim[0] = 4; dim[1] = dim[2] = dim[3] = 1;
      if (iter<(MaxSCLoop-1)) {
	ObitInfoListAlwaysPut(selfCal->info, "solMode", OBIT_string, dim, pmode);
	Obit_log_error(err, OBIT_InfoErr, "Peel self in %s mode",pmode);
     } else {
	ObitInfoListAlwaysPut(selfCal->info, "solMode", OBIT_string, dim, solmod);
  	Obit_log_error(err, OBIT_InfoErr, "Peel self in %s mode",solmod);
    }

      /* Self calibrate using only that field */
      converged = ObitUVSelfCalSelfCal (selfCal, scrUV, init, &noSCNeed, 
					tmpClean->window, err);
      /* If something went wrong - cleanup and bail */
      if (err->error) {
	ObitErrClearErr(err);   /* Clear error messages */
  	Obit_log_error(err, OBIT_InfoWarn, "Peel Self Cal failed - skip peel");
	ObitErrLog(err);
	peeled = -1;  /* Like it never even happened */
	goto cleanup;
      }  /* end failure recovery */
      init = FALSE;
      if (converged || noSCNeed) break;  /* Done? */
      didSC = TRUE;  /* Have done a Selfcal */
      /*  Image Stokes I */ 
      dim[0] = 4;
      stemp[0] = 'I'; stemp[1] = stemp[2] = stemp[3] = ' ';  stemp[4] = 0;
      ObitInfoListAlwaysPut (scrUV->info, "Stokes", OBIT_string, dim, stemp);
      
      /* Image only peel field - using residual data */
      ObitDConCleanVisDeconvolve ((ObitDCon*)tmpClean, err);
      if (err->error) goto cleanup;
    } /* end self cal loop */

    /* Did a self cal occur?  If not just cleanup and go home */
    if (!didSC) goto DonePeel;
  
    /* Copy SN table to inUV */
    numPol = numIF = 0;
    ver = 0;
    inSNTable = newObitTableSNValue ("Peeled SN", (ObitData*)scrUV,
				     &ver, OBIT_IO_ReadOnly, numPol, numIF, 
				     err);
    /* Make sure  created */
    if (inSNTable==NULL) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: No SN table found on %s", routine, scrUV->name);
      goto cleanup;
    }
    ver = 0;
    outSNTable = newObitTableSNValue ("outSN", (ObitData*)inUV,
				      &ver, OBIT_IO_WriteOnly, 
				      inSNTable->numPol, inSNTable->numIF, err);
    /* Make sure  created */
    if (outSNTable==NULL) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: No SN table generated for %s", routine, inUV->name);
      goto cleanup;
    }

    outSNTable = ObitTableSNCopy (inSNTable, outSNTable, err);
    if (err->error) goto cleanup;

    scrUV = ObitUVUnref(scrUV); /* Get rid of scrUV */

    /* Create model dataset with the model corrupted by the fitted SN soln */
    /* Work file with zeroed data */
    /* Pass all input data */
    dim[0] = 1; 
    ObitInfoListAlwaysPut (inUV->info, "passAll", OBIT_bool, dim, &Tr);
    tmpUV = ObitUVUtilCopyZero(inUV, TRUE, NULL, err);
    ObitInfoListAlwaysPut (inUV->info, "passAll", OBIT_bool, dim, &Fl);
    if (err->error) goto cleanup;

    /* Use all components on peel field */
    dim[0] = 1; jtemp = 1;
    ObitInfoListAlwaysPut (tmpSkyModel->info, "BComp", OBIT_long, dim, &jtemp);
    dim[0] = 1; jtemp = 0;
    ObitInfoListAlwaysPut (tmpSkyModel->info, "EComp", OBIT_long, dim, &jtemp);
    /* Model Mode DFT */
    dft = (olong)OBIT_SkyModel_DFT;
    ObitInfoListAlwaysPut (tmpSkyModel->info, "Mode", OBIT_long, dim, &dft);
    /* CCVer 2 */
    dim[0] = 1; jtemp = 2;
    ObitInfoListAlwaysPut (tmpSkyModel->info, "CCVer", OBIT_long, dim, &jtemp);
    /* Replace data */
    ObitInfoListAlwaysPut (tmpSkyModel->info, "REPLACE", OBIT_bool, dim, &Tr);
    /* Factor */
    dim[0] = 1; ftemp = 1.0;
    ObitInfoListAlwaysPut (tmpSkyModel->info, "Factor", OBIT_float, dim, &ftemp);

    /* No translation in Stokes */ 
    dim[0] = 4;
    stemp[0] = stemp[1] = stemp[2] = stemp[3] = ' ';  stemp[4] = 0;
    ObitInfoListAlwaysPut (tmpUV->info, "Stokes", OBIT_string, dim, stemp);
    /* Don't apply calibration */
    dim[0] = 1; jtemp = -1;
    ObitInfoListAlwaysPut (tmpUV->info, "doCalib", OBIT_long, dim, &jtemp);
    jtemp = 0;
    ObitInfoListAlwaysPut (tmpUV->info, "gainUse",  OBIT_long, dim, &jtemp);
    /* Parameters from myInput  */
    ObitInfoListCopyList (myInput, tmpUV->info, peelParms);

    /* Calculate Model data */
    ObitSkyModelSubUV (tmpSkyModel, tmpUV, tmpUV, err);
    if (err->error) goto cleanup;
    
    /* Invert SN table from self cal */
    ver = 0;
    SNInver = ObitTableSNUtilInvert (outSNTable, (ObitData*)tmpUV, &ver, TRUE, err);
    if (err->error) goto cleanup;

    /* Subtract tmpUV from inUV applying inverse of selfcalibration to corrupt */
    dim[0] = 1; jtemp = 1;
    ObitInfoListAlwaysPut (tmpUV->info, "doCalib", OBIT_long, dim, &jtemp);
    jtemp = ver;
    ObitInfoListAlwaysPut (tmpUV->info, "gainUse",  OBIT_long, dim, &jtemp);
    /* Select */
    ObitInfoListAlwaysPut (tmpUV->info, "doCalSelect", OBIT_bool, dim, &Tr);
    /* Pass all data */
    ObitInfoListAlwaysPut (tmpUV->info, "passAll", OBIT_bool, dim, &Tr);
 
    /* Subtract */
    Obit_log_error(err, OBIT_InfoErr, " ******  Subtract Peeled model from UV Data");
    ObitErrLog(err); 
    ver = 0;
    ObitUVUtilVisSub(inUV, tmpUV, inUV, err);
    if (err->error) goto cleanup;
    ObitUVUnref(tmpUV);  /* Done with tmpUV */

    peeled = peelField;  /* Peel occured on field peelField */
        
    /*  Reimage Stokes I */ 
    dim[0] = 4;
    stemp[0] = 'I'; stemp[1] = stemp[2] = stemp[3] = ' ';  stemp[4] = 0;
    ObitInfoListAlwaysPut (inUV->info, "Stokes", OBIT_string, dim, stemp);
    /* Don't Apply calibration */
    dim[0] = 1; jtemp = -1;
    ObitInfoListAlwaysPut (inUV->info, "doCalib", OBIT_long, dim, &jtemp);
    /* Don't need to remake beams  */
    dim[0] = 1;dim[1] = 1;
    ObitInfoListAlwaysPut(myClean->info, "doBeam", OBIT_bool, dim, &Fl);

    /* Start CLEAN fresh */
    noParms = 0;
    ver = 1;
    for (j=0; j<myClean->mosaic->numberImages; j++) {
      CCTab = newObitTableCCValue ("Temp CC", (ObitData*)myClean->mosaic->images[j],
				   &ver, OBIT_IO_ReadWrite, noParms, err);
      ObitTableClearRows ((ObitTable*)CCTab, err); /* Remove the entries and redo */
      CCTab = ObitTableCCUnref(CCTab);
      if  (err->error) Obit_traceback_val (err, routine, myClean->mosaic->images[j]->name, peeled);
    }

    /* reImage/Clean All*/
    Obit_log_error(err, OBIT_InfoErr, " ******  Reimage");
    ObitErrLog(err); 
    ObitDConCleanVisDeconvolve ((ObitDCon*)myClean, err);
    if (err->error) goto cleanup;

    /* Make sure CLEAN indicates some cleaning */
    myClean->Pixels->currentIter = MAX (1, myClean->Pixels->currentIter);

    /* Copy CCs from peeled table (CC # 2) */
    ver = 2;
    noParms = 0;
    peelCCTable = newObitTableCCValue ("Peeled CC", (ObitData*)tmpMosaic->images[0],
				       &ver, OBIT_IO_ReadOnly, noParms, 
				       err);
    /* Make sure  created */
    if (peelCCTable==NULL) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: No CC table 1 found on %s", routine, tmpMosaic->images[0]->name);
      goto cleanup;
    }
    /* Copy Peeled CC table to CC 2 on output image */
    ver = 2;
    outCCTable = newObitTableCCValue ("outCC", (ObitData*)myClean->mosaic->images[peelField-1],
				      &ver, OBIT_IO_WriteOnly, peelCCTable->noParms, 
				      err);
    /* Make sure object created */
    if (outCCTable==NULL) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: No CC table 2 created on %s", routine, 
		     myClean->mosaic->images[peelField-1]->name);
      goto cleanup;
    }

    /* Clear any existing entries in output table */
    ObitTableClearRows((ObitTable*)outCCTable, err);
    if (err->error) goto cleanup;

    /*ObitTableCCUtilAppend (peelCCTable, outCCTable, 1, 0, err);*/
    /* Copy Peel CC to CC2 reverting any 2D to 3D shift */
    MoveBack (peelCCTable, outCCTable, ccShift, err);
    peelCCTable = ObitTableCCUnref(peelCCTable);
    outCCTable  = ObitTableCCUnref(outCCTable);

    /* Cleanup - delete temp Imager and CLEAN */
  DonePeel:   
    ver = -1;
    /* ObitUVZapTable (inUV, "AIPS SN", ver, err);          delete Peel solutions */
 
    PeelFlux = 1.0e20;  /* only once */
  } /* End peeling */

  /* Cleanup - delete temp Imager and CLEAN */
 cleanup:   
  if (tmpClean) ObitImageMosaicZapImage (tmpClean->mosaic, -1, err); /* Delete temp images */
  tmpImager   = ObitUVImagerUnref(tmpImager);
  tmpSkyModel = ObitSkyModelUnref(tmpSkyModel);
  tmpClean    = ObitDConCleanVisUnref(tmpClean);
  selfCal     = ObitUVSelfCalUnref(selfCal);
  inSNTable   = ObitTableSNUnref(inSNTable);
  outSNTable  = ObitTableSNUnref(outSNTable);
  SNInver     = ObitTableSNUnref(SNInver);
  if (err->error) Obit_traceback_val (err, routine, myClean->name, peeled);
  
  return peeled;
} /* end ObitUVPeelUtilPeel */

/*----------------------Private functions---------------------------*/
/**
 * Convert a 2d IMAGE TO 3d
 * Fiddle descriptor and the CC table 1
 * \param image   Image to convert
 * \parms ccShift [out] x,y pixel shifts added to CC components
 * \param err     Error/message stack
 */
void Convert2Dto3D (ObitImage *image,ofloat ccShift[2],  ObitErr* err)
{
  ofloat xPxOff, yPxOff, xCen, yCen, xShift, yShift, xPix, yPix;
  odouble ra, dec;
  olong nx, ny, CCVer, noParms, irow;
  ObitTableCC *CCTable=NULL;
  ObitTableCCRow *CCRow=NULL;
  gchar *routine = "Convert2Dto3D";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageIsA(image));

  /* Init output */
  ccShift[0] = 0.0;
  ccShift[1] = 0.0;

  /* Open image */
  ObitImageOpen (image, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);

  /* If already 3D just return */
  xPxOff = image->myDesc->xPxOff;
  yPxOff = image->myDesc->yPxOff;
  if ((image->myDesc->do3D) && (xPxOff==0.0) && (yPxOff==0.0)) return;

  /* Old image stuff */
  nx   = image->myDesc->inaxes[0];
  ny   = image->myDesc->inaxes[1];
  xCen = image->myDesc->crpix[0];
  yCen = image->myDesc->crpix[1];
  xPix = 1.0 + nx/2.0;
  yPix = 1.0 + ny/2.0;
  ObitSkyGeomWorldPos(xPix, yPix, 
		      image->myDesc->crval[0], image->myDesc->crval[1], 
		      image->myDesc->crpix[0], image->myDesc->crpix[1],
		      image->myDesc->cdelt[0], image->myDesc->cdelt[1],
		      image->myDesc->crota[image->myDesc->jlocd], 
		      &image->myDesc->ctype[image->myDesc->jlocr][4],
		      &ra, &dec);

  /* new values -
     New reference position = center */
  image->myDesc->crval[0] = ra;
  image->myDesc->crval[1] = dec;
  image->myDesc->crpix[0] = 1 + nx/2;
  image->myDesc->crpix[1] = 1 + ny/2;
  image->myDesc->xPxOff   = 0.0;
  image->myDesc->yPxOff   = 0.0;
  image->myDesc->do3D     = TRUE;

  /* How much to shift (add to) CCs */
  xShift = xPxOff * image->myDesc->cdelt[0];
  yShift = yPxOff * image->myDesc->cdelt[1];

  /* Save shift  */
  ccShift[0] = xShift;
  ccShift[1] = yShift;

  /* Close image */
  image->myStatus = OBIT_Modified;  /* Force update */
  ObitImageClose (image, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);

  /* Fix Clean components */
  CCVer = 1;
  noParms = 0;
  CCTable = newObitTableCCValue ("Peeled CC", (ObitData*)image,
				 &CCVer, OBIT_IO_ReadWrite, noParms, 
				 err);
  /* Open CC table */
  ObitTableCCOpen (CCTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
  /* Create table row */
  CCRow = newObitTableCCRow (CCTable);

  /* Loop over table */
  for (irow=1; irow<=CCTable->myDesc->nrow; irow++) {
    ObitTableCCReadRow (CCTable, irow, CCRow, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);

    /* Update row */
    CCRow->DeltaX += xShift;
    CCRow->DeltaY += yShift;

    ObitTableCCWriteRow (CCTable, irow, CCRow, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
  } /* end loop over table */

  /* Close Table */
  ObitTableCCClose (CCTable, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
  /* Cleanup */
  CCTable = ObitTableCCUnref(CCTable);
  CCRow   = ObitTableCCRowUnref (CCRow);
  
} /* end Convert2Dto3D */

/**
 * Copy contents on inCC to the end of outCC adjusting positions.
 * \param inCC    input CC table
 * \param outCC   output CC table
 * \parms ccShift x,y pixel shifts to subtract from CC components
 * \param err     Error/message stack
 */
void MoveBack (ObitTableCC *inCC, ObitTableCC *outCC, ofloat ccShift[2],  
	       ObitErr* err)
{
  olong irow, orow;
  ObitTableCCRow *CCRow=NULL;
  gchar *routine = "MoveBack";

  /* error checks */
  if (err->error) return;
  g_assert (ObitTableCCIsA(inCC));
  g_assert (ObitTableCCIsA(outCC));

  /* Open CC tables */
  ObitTableCCOpen (inCC,  OBIT_IO_ReadOnly, err);
  ObitTableCCOpen (outCC, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inCC->name);
  
  /* Create table row */
  CCRow = newObitTableCCRow (outCC);

  /* Loop over table */
  for (irow=1; irow<=inCC->myDesc->nrow; irow++) {
    ObitTableCCReadRow (inCC, irow, CCRow, err);
    if (err->error) Obit_traceback_msg (err, routine, inCC->name);

    /* Update row */
    CCRow->DeltaX -= ccShift[0];
    CCRow->DeltaY -= ccShift[1];

    orow = -1;
    ObitTableCCWriteRow (outCC, orow, CCRow, err);
    if (err->error) Obit_traceback_msg (err, routine, inCC->name);
  } /* end loop over table */

  /* Close Tables */
  ObitTableCCClose (inCC,  err);
  ObitTableCCClose (outCC, err);
  if (err->error) Obit_traceback_msg (err, routine, inCC->name);
  
  /* Cleanup */
  CCRow   = ObitTableCCRowUnref (CCRow);
  
} /* end MoveBack */

