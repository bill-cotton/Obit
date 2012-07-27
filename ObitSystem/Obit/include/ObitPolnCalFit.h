/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012                                               */
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
#ifndef OBITPOLNCALFIT_H 
#define OBITPOLNCALFIT_H

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitSource.h"
#include "ObitSourceList.h"
#include "ObitAntennaList.h"
#include "ObitComplex.h"
#include "ObitTableCP.h"
#include "ObitTablePD.h"
#include "ObitTableBP.h"
#ifdef HAVE_GSL
#include <gsl/gsl_multifit_nlin.h>
#endif /* HAVE_GSL */ 

/*-------- Obit:  Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPolnCalFit.h
 *
 * ObitPolnCalFit Class for fitting polarization calibration to UV data.
 *
 * This class does least squares fitting of a full nonlinear model
 * To instrumental polarization and calibration terms and optionally source 
 * polarization parameters.
 * Results are stored in a combination of AIPS PD (instrumental poln),
 * AIPS CP (source poln) and AIPS BP (gain corrections including R-L phase).
 * 
 * \section ObitPolnCalFitaccess Creators and Destructors
 * An ObitPolnCalFit will usually be created using ObitPolnCalFitCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitPolnCalFit should always be made using the
 * #ObitPolnCalFitRef function which updates the reference count in the object.
 * Then whenever freeing an ObitPolnCalFit or changing a pointer, the function
 * #ObitPolnCalFitUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */
/*---------------Private structures----------------*/
/**
 * \enum polnParmType
 * enum for I/O file type.
 * This specifies the type of underlying data file.
 */
enum polnParmType {
  /** Unspecified */
  polnParmUnspec=0, 
  /** Source parameter */
  polnParmSou,
  /** Antenna parameter */
  polnParmAnt,
  /** Phase difference  */
  polnParmPD
}; /* end enum polnParmType */

/** typedef for enum for polnParmType . */
typedef enum polnParmType PolnParmType;

/* Threaded function argument */
typedef struct {
  /* Obit error stack object */
  ObitErr      *err;
  /** thread    */
  ObitThread    *thread;
  /** thread number, >0 -> no threading   */
  olong        ithread;
  /** Do error analysis? */
  gboolean doError;
  /** First 0-rel datum in Data/Wt arrays */
  olong lo;
  /** Last 0-rel datum  in Data/Wt arrays */
  olong hi;
  /** selected antenna, 0-> all */
  olong selAnt;
  /** selected source, 0-> all */
  olong selSou;
 /** Number of visibilities being fitted */
  olong nvis;
 /** Number of data points being fitted (nvis*4) */
  olong ndata;
  /** Array of input data arrays, [10 x nvis]
      each row: 2*parallactic angle (rad), blank, 
      RR_r, RR_i, LL_r, LL_i, RL_r, RL_i,  LR_r, LR_i   */
  ofloat *inData;
  /** Array of input data weights (~ 1/var), [4 x nvis]
      each row: RR, LL, RL,  LR   */
  ofloat *inWt;
  /** Source no of visibility, in range [0,nsou-1] */
  olong *souNo;
  /** Antenna numbers  of visibility, in range [0,nant-1] */
  olong *antNo;
  /** Number of antennas */
  olong nant;
  /** Reference antenna */
  olong refAnt;
  /** Fit only reference antenna R-L phase */
  gboolean doFitRL;
  /** R-L (or X-Y) phase difference */
  odouble PD;
  /** Error estimate R-L (or X-Y) phase difference */
  odouble PDerr;
  /** Fit Stokes I poln, per cal */
  gboolean *doFitI;
  /** Fit fit source linear poln, per cal */
  gboolean *doFitPol;
  /** Fit fit Stokes V poln, per cal */
  gboolean *doFitV;
  /** Antenna parameters 4 x nant, 
      each row: D1_r, D1_i, D2_r, D2_i */
  odouble *antParm;
  /** Antenna error estimates 4 x nant */ 
  odouble *antErr;
  /** Antenna parameters fit flags, 4 x nant */
  gboolean **antFit;
  /** Antenna parameters number, 4 x nant */
  olong **antPNumb;
  /** Number of calibrator sources */
  olong nsou;
  /** Source parameters 4 x nsou, 
     each row: I, Poln_amp, Poln_RL (rad), V */
  odouble *souParm;
  /** Source error estimates 4 x nsou, */
  odouble *souErr;
  /** Source parameters fit flags, 4 x nsou */
  gboolean **souFit;
  /** Source parameters number, 4 x nsou */
  olong **souPNumb;
  /** Source Ids in SU table, Data */
  olong *souIDs;
  /** Inverse Source Ids lookup, index of data source ID-1 
      gives calibrator number */
  olong *isouIDs;
  /** Number of parameters being fitted */
  olong nparam;
  /** Number of valid observations */
  olong nobs;
  /** Chi squared of fit */
  ofloat ChiSq;
  /** Sum of parallel residuals */
  odouble sumParResid;
  /** Sum of cross residuals */
  odouble sumXResid;
  /** Parameter type */
  PolnParmType paramType;
  /** Parameter number */
  olong paramNumber;
  /** Sum of first derivatives */
  odouble sumDeriv;
  /** Sum of second derivatives */
  odouble sumDeriv2;
  /** Sum of weights */
  odouble sumWt;
  /** Frequency of data */
  odouble freq;
  /** R-L Phase of calibrators (rad) at Freq */
  ofloat *RLPhase;
  /** Central channel (1-rel) of fit */
  olong Chan;
  /** IF number (1-rel) */
  olong IFno;
  /** complex work arrays for antenna based parameters, max 100 ant */
  dcomplex RS[100], RD[100], LS[100], LD[100],
    RSc[100], RDc[100], LSc[100], LDc[100],
    PR[100], PRc[100], PL[100], PLc[100];
  odouble SR[100], DR[100], SL[100], DL[100];
  /** Max antenna number */
  olong maxAnt;
  /** Solution method */
  gchar solnType[5];
} PolnFitArg;

/*--------------Class definitions-------------------------------------*/
/** ObitPolnCalFit Class structure. */
typedef struct {
#include "ObitPolnCalFitDef.h"   /* this class definition */
} ObitPolnCalFit;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitPolnCalFit
 * returns a ObitPolnCalFit*.
 * in = object to unreference
 */
#define ObitPolnCalFitUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitPolnCalFit.
 * returns a ObitPolnCalFit*.
 * in = object to reference
 */
#define ObitPolnCalFitRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitPolnCalFitIsA(in) ObitIsA (in, ObitPolnCalFitGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitPolnCalFitClassInit (void);

/** Public: Default Constructor. */
ObitPolnCalFit* newObitPolnCalFit (gchar* name);

/** Public: Create/initialize ObitPolnCalFit structures */
ObitPolnCalFit* ObitPolnCalFitCreate (gchar* name);
/** Typedef for definition of class pointer structure */
typedef ObitPolnCalFit* (*ObitPolnCalFitCreateFP) (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitPolnCalFitGetClass (void);

/** Public: Copy (deep) constructor. */
ObitPolnCalFit* ObitPolnCalFitCopy  (ObitPolnCalFit *in, 
				     ObitPolnCalFit *out, ObitErr *err);

/** Public: Copy structure. */
void ObitPolnCalFitClone (ObitPolnCalFit *in, ObitPolnCalFit *out, 
			  ObitErr *err);

/** Public: Fit a UV data */
void ObitPolnCalFitFit (ObitPolnCalFit* in, ObitUV *inUV, 
			ObitUV *outUV, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void(*ObitPolnCalFitFitFP) (ObitPolnCalFit* in, ObitUV *inUV, 
				    ObitUV *outUV, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitPolnCalFitClassDef.h"
} ObitPolnCalFitClassInfo; 

#endif /* OBITPOLNCALFIT_H */ 
