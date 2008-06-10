/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2007                                               */
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
#include "ObitImageMosaic.h"
#include "ObitSkyModel.h"
#include "ObitUVSelfCal.h"
#include "ObitTableCCUtil.h"
#include "ObitTableSNUtil.h"
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVPeelUtil.c
 * ObitUVPeelUtil Utility function definitions.
 */

/*---------------Private function prototypes----------------*/
/*----------------------Public functions---------------------------*/
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
 *  \li PeelNiter     i Niter for peel.CLEAN
 *  \li PeelMinFlux   f Minimum Peel Clean component (Jy)
 *  \li refAnt        i Peel SC Reference antenna
 *  \li SNRMin        f Min. allowed SNR in peel selfcal
 *  \li AvgPol        b Avg. poln in peel self cal?
 *  \li AvgIF         b Avg. IFs in peel self cal?
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
 *  \li dispURL       s The URL of the display server to use
 * \param inUV    Data to be peeled, on return all SN tables will be removed.
 *                and peeled source will have been subtracted
 * \param myClean Clean object which has previously been CLEANed and which 
 *                has a field with a CC peak in excess of PeelFlux.
 *                Need for peeling based on values of myClean->peakFlux and
 *                PeelFlux.  The test for individual fields is the maximum as 
 *                determined from the CLEAN components by 
 *                #ObitImageMosaicMaxField
 *                On output, this will include info on last CLEAN
 * \param err     Error/message stack
 * \return the 1-rel field number of the peeled field or -1 if no peel
 */
olong ObitUVPeelUtilPeel (ObitInfoList* myInput, ObitUV* inUV, 
			 ObitDConCleanVis *myClean, ObitErr* err)
{
  olong              peeled=-1;
  ObitInfoType      type;
  ObitUVImager      *tmpImager=NULL; 
  ObitImageMosaic   *tmpMosaic=NULL; 
  ObitUV            *scrUV=NULL, *tmpUV=NULL;
  ObitSkyModel      *tmpSkyModel=NULL; 
  ObitDConCleanVis  *tmpClean=NULL; 
  ObitUVSelfCal     *selfCal = NULL;
  ObitTableCC       *peelCCTable=NULL, *outCCTable=NULL;
  ObitTableSN       *inSNTable=NULL, *outSNTable=NULL, *SNInver=NULL;
  olong         i,  jtemp, peelField, *bcomp=NULL, *ecomp=NULL, ignore[1]={0};
  olong         dft, iter, MaxSCLoop=1;
  oint         noParms, numPol, numIF;
  olong        ver;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat       PeelFlux, ftemp, xCells, yCells;
  gboolean     converged, didSC=FALSE, Fl=FALSE, Tr=TRUE, init, noSCNeed; 
  gchar        stemp[5];
  gchar        *imgParms[] = {  /* Imaging, weighting parameters */
    "PBCor", "antSize", 
    "Robust", "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "WtPower",
    "MaxBaseline", "MinBaseline", "rotate", "Beam",  "xCells", "yCells", 
    "nx", "ny",  "nxBeam", "nyBeam", "dispURL",
    NULL
  };
  gchar        *CLEANParms[] = {  /* Clean parameters */
    "Gain", "minFlux", "Niter", "minPatch", "Beam", 
    "Mode", "CCFilter", "maxPixel", "dispURL", "autoWindow",
    NULL
  };
  gchar        *SCParms[] = {  /* Self cal parameters */
    "refAnt", "WtUV", "avgPol", "avgIF", "doMGM", "SNRmin", 
    "minNo", "doSmoo", "modelFlux", "modelPos", "modelParm",
    "dispURL", /* "prtLv", */
    NULL
  };
  gchar        *peelParms[] = {  /* peeling parameters */
    "PBCor", "antSize", 
  };
  gchar *solmod = "P   ", *soltyp = "L1  ";
  gchar *routine = "ObitPeelUtilPeel";

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

    /* Get ImageMosaic for brightest field if any exceeds PeelFlux */
    tmpMosaic = ObitImageMosaicMaxField (myClean->mosaic, PeelFlux, 
					 ignore, &peelField, err);
    if (err->error) Obit_traceback_val (err, routine, myClean->name, peeled);
    if (tmpMosaic==NULL) goto donePeel;

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

    /* Scratch file */
    scrUV = newObitUVScratch (inUV, err);
    /* Give more sensible name */
    if (scrUV->name) g_free(scrUV->name);
    scrUV->name = g_strdup("Peel data");
    dim[0] = 4;
    stemp[0] = stemp[1] = stemp[2] = stemp[3] = ' ';  stemp[4] = 0;
    ObitInfoListAlwaysPut (inUV->info, "Stokes", OBIT_string, dim, stemp);

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

    /* Temporary Imager */
    tmpImager = ObitUVImagerCreate2("Peel imager", scrUV, tmpMosaic, err);
    if (err->error) goto cleanup;

    /* Create temp SkyModel */
    tmpSkyModel = ObitSkyModelCreate("Peel SkyModel", tmpMosaic);
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
 
    /* Use DFT model */
    dim[0] = dim[1] = 1;
    dft = (olong)OBIT_SkyModel_DFT;
    ObitInfoListAlwaysPut (tmpClean->info, "Mode", OBIT_long, dim, &dft);
    
    /* Create self cal object */
    selfCal = ObitUVSelfCalCreate ("SelfCal", tmpSkyModel);
    /* Set Parameters  */
    dim[0] = 4;
    ObitInfoListAlwaysPut(selfCal->info, "solMode", OBIT_string, dim, solmod);
    ObitInfoListAlwaysPut(selfCal->info, "solType", OBIT_string, dim, soltyp);
    dim[0] = 1;
    /* Solution interval */
    ftemp = 1.0; dim[0] = dim[1]= dim[2] = 1;
    ObitInfoListGetTest (myInput, "PeelSolInt", &type, dim, &ftemp);
    ObitInfoListAlwaysPut(selfCal->info, "solInt", OBIT_float, dim, &ftemp);
    /* Min acceptable peak for Self cal */
    ftemp = PeelFlux*0.5;
    ObitInfoListPut (selfCal->info, "minFluxPSC", OBIT_float, dim, &ftemp, err);
    /* Minimum SNR */
    ftemp = 3.0; dim[0] = dim[1]= dim[2] = 1;
    ObitInfoListGetTest (myInput, "PeelSolInt", &type, dim, &ftemp);
    ObitInfoListAlwaysPut(selfCal->info, "SNRMin", OBIT_float, dim, &ftemp);
    ObitInfoListAlwaysPut(selfCal->info, "avgPol", OBIT_bool, dim, &Tr);
    ObitInfoListAlwaysPut(selfCal->info, "doSmoo", OBIT_bool, dim, &Tr);
    /* Copy control info */
    ObitInfoListCopyList (myInput, selfCal->info, SCParms);
    
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
      
      /* Self calibrate using only that field */
      converged = ObitUVSelfCalSelfCal (selfCal, scrUV, init, &noSCNeed, 
					tmpClean->window, err);
      if (err->error) goto cleanup;
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
    if (!didSC) goto donePeel;
  
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
    tmpUV = ObitUVUtilCopyZero(inUV, TRUE, NULL, err);
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
    SNInver = ObitTableSNUtilInvert (outSNTable, (ObitData*)tmpUV, &ver, err);
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

    /* reImage/Clean All*/
    Obit_log_error(err, OBIT_InfoErr, " ******  Reimage");
    ObitErrLog(err); 
    ObitDConCleanVisDeconvolve ((ObitDCon*)myClean, err);
    if (err->error) goto cleanup;

    /* Copy CCs from peeled table to CC 1 on output image */
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
    ObitTableCCUtilAppend (peelCCTable, outCCTable, 1, 0, err);
    peelCCTable = ObitTableCCUnref(peelCCTable);
    outCCTable  = ObitTableCCUnref(outCCTable);

    /* Cleanup - delete temp Imager and CLEAN */
  donePeel:   
    ver = -1;
    ObitUVZapTable (inUV, "AIPS SN", ver, err);          /* delete Peel solutions */
 
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
