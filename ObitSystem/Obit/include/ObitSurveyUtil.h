/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012,2016                                          */
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
#ifndef OBITSURVEYUTIL_H 
#define OBITSURVEYUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitTableVL.h"
#include "ObitTableVZ.h"
#include "ObitFitRegionList.h"
#include "ObitPrinter.h"
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSurveyUtil.h
 * ObitSurvey class utility routine definition.
 */

/** Public: Define structure for Generic survey calibration constants  */
typedef struct {
  /** Flux density scale as factor */
  ofloat fluxScale;
  /** RA, dec bias (deg) */
  ofloat biasRA, biasDec;
  /** Calibration component of  position error (squared) */
  ofloat  calRAEr, calDecEr;
  /** Mean and uncertainty of CLEAN  bias  */
  ofloat   ClnBiasAv, ClnBiasEr;
  /** Amplitude calibration uncertainty as fraction */
  ofloat calAmpEr;
  /** Axis size calibration uncertainty as fraction */
  ofloat calSizeEr;
  /** Polarization calibration error  */
  ofloat calPolEr;
  /** CLEAN convolving beam, major, minor, pa in deg  */
  ofloat beamMaj, beamMin, beamPA;
} ObitSurveyGenCalParms;

/*---------------Public functions---------------------------*/
/** Public: Print all contents of a VL table */
void ObitSurveyUtilVLPrint (ObitTableVL *in, ObitImage *image, FILE  *prtFile, 
			    ObitErr *err);
/** Public: Printed selected Generic survey catalog entries */
gboolean ObitSurveyGenPrint (ObitPrinter *printer, ObitData *data, olong VLVer, 
			      gboolean first, gboolean last, ObitErr* err);

/** Public: Printed selected NVSS survey catalog entries */
gboolean ObitSurveyNVSSPrint (ObitPrinter *printer, ObitData *data, olong VLVer, 
			      gboolean first, gboolean last, ObitErr* err);

/** Public: Printed selected VLSS survey catalog entries */
gboolean ObitSurveyVLSSPrint (ObitPrinter *printer, ObitData *data, olong VLVer, 
			      gboolean first, gboolean last, ObitErr* err);

/** Create default Generic survey calibration constants */
ObitSurveyGenCalParms* 
ObitSurveyGetCalParms(ofloat fluxScale, ofloat biasRA, ofloat biasDec, 
		      ofloat calRAEr, ofloat calDecEr, 
		      ofloat ClnBiasAv, ofloat ClnBiasEr, 
		      ofloat calAmpEr, ofloat calSizeEr, ofloat calPolEr,
		      ofloat beamMaj, ofloat beamMin, ofloat beamPA);

/** Public: Generic Survey (VL table format) error calculation */
void ObitSurveyGenCorErr (ObitSurveyGenCalParms *calParms, odouble *ra, odouble *dec, 
			  ofloat *peak, ofloat *major, ofloat *minor, ofloat *posang,
			  ofloat qcent, ofloat ucent,  ofloat *pflux, 
			  ofloat irms,  ofloat prms,   
			  ofloat *flux, ofloat *eflux, ofloat *epflux, 
			  ofloat *pang, ofloat *epang,  ofloat *errra, ofloat *errdec, 
			  ofloat *cmajor, ofloat *cminor, ofloat *cpa, 
			  ofloat *emajor, ofloat *eminor, ofloat *epa, gboolean rflag[4]);

#endif /* OBITSURVEYUTIL_H */ 
