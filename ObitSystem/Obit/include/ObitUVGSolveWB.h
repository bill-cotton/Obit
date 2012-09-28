/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010-2012                                          */
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
#ifndef OBITUVGSOLVEWB_H 
#define OBITUVGSOLVEWB_H 

#include "ObitUV.h"
#include "ObitTableSN.h"
#include "ObitSkyModel.h"
#include "ObitUVGSolve.h"
#if HAVE_GSL==1  /* GSL stuff */
#include "gsl/gsl_blas.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_multifit_nlin.h"
#endif /* GSL stuff */

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVGSolveWB.h
 * ObitUVGSolveWB Class for wide band uv data self calibration
 *
 * This class is derived from the #ObitUVGSolve class.
 *
 * Routines determine calibration for an ObitUV and write an SN table.
 * Control parameters are on the info member.
 * \li "subA"    OBIT_int   (1,1,1) Selected subarray (default 1)
 * \li "solInt"  OBIT_float (1,1,1) Solution interval (min). (default scan)
 * \li "refAnt"  OBIT_int   (1,1,1) Ref ant to use. (default 1)
 * \li "refAnts" OBIT_int   (?,1,1) list of Ref ants to use. (default (1))
 * \li "avgPol"  OBIT_bool  (1,1,1) True if RR and LL to be averaged (false)
 * \li "avgIF"   OBIT_bool  (1,1,1) True if all IFs to be averaged (false)
 *                                  otherwise individual IF fits.
 * \li "doTwo"   OBIT_bool  (1,1,1) Use 2 BL combinations as well as 1 BL. (true)
 * \li "minSNR"  OBIT_float (1,1,1) Minimum acceptable SNR (5)
 * \li "doMGM"   OBIT_bool  (1,1,1) True then find the mean gain modulus (true)
 * \li "elevMGM" OBIT_float (1,1,1) Min. elevation to include in mean gain modulus
 * \li "solType" OBIT_string (4,1,1 Solution type '  ', 'L1',  (' ')
 * \li "solMode" OBIT_string (4,1,1 Solution mode: 'A&P', 'P', 'P!A', 'GCON' ('P')
 * \li "minNo"   OBIT_int   (1,1,1) Min. no. antennas. (default 4)
 * \li "antWt"   OBIT_float (*,1,1) Antenna weights. (default 1.0)
 * \li "UVR_Full"OBIT_float (2,1,1) Range of baseline lengths with full weight
 *                                  (kilolamda). If none is given then 
 *                                  derive one if possible.
 * \li "WtUV"    OBIT_float (1,1,1) Weight outside of UVRANG. (default 1.0)
 * \li "prtLv"   OBIT_int   (1,1,1) Print level (default no print)
 * 
 * \section ObitUVGSolveWBaccess Creators and Destructors
 * An ObitUVGSolveWB will usually be created using ObitUVGSolveWBCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitUVGSolveWB should always be made using the
 * #ObitUVGSolveWBRef function which updates the reference count in the object.
 * Then whenever freeing an ObitUVGSolveWB or changing a pointer, the function
 * #ObitUVGSolveWBUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*---------------Private structures----------------*/
/* Baseline data structure */
typedef struct {
  /* Baseline antenna numbers */
  olong ant1, ant2;
  /* Number of IFs */
  olong numIF;
  /* Number of poln */
  olong numPoln;
  /* Number of frequency averages per IF */
  olong numFreq;
  /* 1-rel IF number for data arrays */
  olong *IFindex;
  /* 1-rel Stokes number for data arrays */
  olong *Stokeindex;
  /* Array of weights per IF/poln (IF varies fastest )*/
  ofloat *WtPolnIF;
  /* Array of IF frequency offset (GHz) arrays from reference Frequency,
   fblank => no data*/
  ofloat **dFreq;
  /* Array of accumulated weights per IF/poln, dimensioned frequency, 
     IF varies faster than poln */
  ObitFArray **wtArray;
  /* Array of accumulated visibilities per IF/poln, dimensioned frequency, 
     IF varies faster than poln */
  ObitCArray **visArray;
  /* Array of accumulated visibility phases per IF/poln, dimensioned frequency  
     IF varies faster than poln */
  ObitFArray **phArray;
} BaselineData;

/* Scan data structure */
typedef struct {
  /* Maximum antenna number (1-rel) */
  olong maxAnt;
  /* Number of baselines */
  olong numBase;
  /* Number of spectral channels */
  olong numChan;
  /* Reference channel */
  ofloat refChan;
  /* Time averaging (days) - bin size */
  ofloat timeAvg;
  /* Frequency averaging (Hz) - signed bin size, per IF */
  ofloat *freqAvg;
  /* Average polarization?  If TRUE, only one poln. in data */
  gboolean avgPoln;
  /* Average IF?  otherwise single IF fits */
  gboolean avgIF;
  /* Use two baseline combinations as well as 1 baseline */
  gboolean doTwo;
  /* Center time of observations (days) */
  odouble timec;
  /* Actual time interval (days) */
  ofloat timei;
  /* RMS Phase residual */
  ofloat RMSRes;
  /* Subarray number */
  olong suba;
  /* Source Id if present else -1 */
  olong sid;
  /* FQ id if present else -1 */
  olong fqid;
  /* Array of baseline structures, antenna 1 varies fastest */
  BaselineData **BLData;
  /* Current (1-rel) antenna number being solved for */
  olong curAnt;
  /* Current (0-rel) polarization being solved for */
  olong curPoln;
  /* Current (0-rel) IF being solved for */
  olong curIF;
  /* Array with Frequency offsets for antenna stacked data arrays (GHz) */
  ObitFArray *antFreqOff;
  /* Array of stacked antenna phases per antenna/poln (ant varies fastest) */
  ObitFArray **antStackPh;
  /* Array of stacked antenna Weights per antenna/poln (ant varies fastest) */
  ObitFArray **antStackWt;
} ScanData;

/* Fringe fit threaded function argument */
typedef struct {
  /* scan data object */
  ScanData     *data;
  /* First (1-rel) baseline to process this thread */
  olong        firstBL;
  /* Highest (1-rel) baseline to process this thread  */
  olong        lastBL;
  /* zero rel polarization being solved for  */
  olong        iPoln;
  /* Total number of data points  */
  olong        ndata;
  /*  Total number of coefficients being solved for */
  olong        ncoef;
  /*  number of coefficients per antenna */
  olong        np;
  /*  reference antenna */
  olong        refAnt;
  /* Work coefficient array */
  ofloat       *lcoef;
  /* thread number , >0 -> no threading  */
  olong        ithread;
  /* Thread */
  ObitThread *thread;
#if HAVE_GSL==1  /* GSL stuff */
  /* Coefficient array */
  gsl_vector   *coef;
  /* Residual array */
  gsl_vector   *f;
  /* Jacobean matrix */
  gsl_matrix   *J;
#endif /* GSL stuff */
} FringeFitFuncArg;
/*--------------Class definitions-------------------------------------*/
/** ObitUVGSolveWB Class structure. */
typedef struct {
#include "ObitUVGSolveWBDef.h"   /* this class definition */
} ObitUVGSolveWB;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUVGSolveWB
 * returns a ObitUVGSolveWB*.
 * in = object to unreference
 */
#define ObitUVGSolveWBUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUVGSolveWB.
 * returns a ObitUVGSolveWB*.
 * in = object to reference
 */
#define ObitUVGSolveWBRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVGSolveWBIsA(in) ObitIsA (in, ObitUVGSolveWBGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVGSolveWBClassInit (void);

/** Public: Default Constructor. */
ObitUVGSolveWB* newObitUVGSolveWB (gchar* name);

/** Public: Create/initialize ObitUVGSolveWB structures */
ObitUVGSolveWB* ObitUVGSolveWBCreate (gchar* name);
/** Typedef for definition of class pointer structure */
typedef ObitUVGSolveWB* (*ObitUVGSolveWBCreateFP) (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitUVGSolveWBGetClass (void);

/** Public: Copy (deep) constructor. */
ObitUVGSolveWB* ObitUVGSolveWBCopy  (ObitUVGSolveWB *in, ObitUVGSolveWB *out, 
				   ObitErr *err);

/** Public: Copy structure. */
void ObitUVGSolveWBClone (ObitUVGSolveWB *in, ObitUVGSolveWB *out, ObitErr *err);

/** Public: Determine calibration from UV data divided by model. */
ObitTableSN* ObitUVGSolveWBCal (ObitUVGSolveWB *in, ObitUV *inUV, ObitUV *outUV, 
				ObitUVSel *sel, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVGSolveWBClassDef.h"
} ObitUVGSolveWBClassInfo; 

#endif /* OBITFUVGSOLVEWB_H */ 
