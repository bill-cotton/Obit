/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012,2015                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
/*  Define the basic components of the ObitPolnCalFit structure      */
/**
 * \file ObitPolnCalFitDef.h
 * ObitPolnCalFit structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** Do Error analysis: */
gboolean doError;
/** Fit global R-L phase: */
gboolean doFitRL;
/** blank failed solns? else default parameters: */
gboolean doBlank;
/** Fit global X & Y feed gains?: */
gboolean doFitGain;
/** Fit Feed orientations? */
gboolean doFitOri;
/** Are the feeds circularly polarized? */
gboolean isCircFeed;
/** R-L (or X-Y) phase difference */
odouble PD;
/** Error estimate R-L (or X-Y) phase difference */
odouble PDerr;
/** PD parameter number */
olong PDPNumb;
/** Fit source IPol per source */
gboolean *doFitI;
/** Fit source VPol per source: */
gboolean *doFitV;
/** Fit source polarization per source */
gboolean *doFitPol;
/** Input UV descriptor */
ObitUVDesc *inDesc;
/** Output UV descriptor */
ObitUVDesc *outDesc;
/** Number of channels in poln soln  */
olong ChWid;
/** Spacing of  poln soln  */
olong ChInc;
/** Source (CP) table to write */
olong CPSoln;
/** Instrumental (PD) table to write */
olong PDSoln;
/** Bandpass (BP) table to write */
olong BPSoln;
/** If >0 then BP table BPVer was applied to input */
olong doBand;
/** input BP table  */
olong BPVer;
/** First IF selected in output  */
olong BIF;
/** First Channel selected in output  */
olong BChan;
/** Current channel selected in output  */
olong Chan;
/** Current IF selected in output  */
olong IFno;
/** Reference Antenna (1-rel) number */
olong refAnt;
/** Number of antennas */
olong nant;
/** Antenna parameters 4 x nant, 
    each row: OriR/X, ElipR/X, OriL/Y, ElipL/Y
*/
odouble *antParm;
/** Antenna error estimates 4 x nant */ 
odouble *antErr;
/* Antenna parameters fit flags, 4 x nant */
gboolean **antFit;
/* Have data for antenna?*/
gboolean *gotAnt;
/** Antenna parameters number, 4 x nant */
olong **antPNumb;
/** Antenna gains 2 x nant, each row: Gain_X, gain_Y */
odouble *antGain;
/** Antenna gains error estimates 2 x nant */ 
odouble *antGainErr;
/* Antenna parameters fit flags, 2 x nant */
gboolean **antGainFit;
/** Antenna gain parameters number, 2 x nant */
olong **antGainPNumb;
/* Number of calibrator sources */
olong nsou;
/** Source parameters 4 x nsou, 
    each row: I, Poln_amp, Poln_RL (rad), V */
odouble *souParm;
/** Source error estimates 4 x nsou, */
odouble *souErr;
/** Previous valid source parameters  */
odouble *lastSouParm;
/* Source parameters fit flags, 4 x nsou */
gboolean **souFit;
/** Source parameters number, 8 x nant */
olong **souPNumb;
/** Source Ids in SU table, data */
olong *souIDs;
/** Inverse Source Ids lookup, index of data source ID-1 
    gives calibrator number */
olong *isouIDs;
/** reference frequency of data */
odouble refFreq;
/**  frequency of current data */
odouble freq;
/**  Initial Chi**2 */
odouble initChiSq;
/**  Chi**2 of model fit */
odouble ChiSq;
/** Parallel, cross hand RMS residual */
ofloat ParRMS, XRMS;
/** Input R-L Phase of calibrators (deg) at refFreq */
ofloat *RLPhaseIn;
/** R-L Phase of calibrators (deg) at freq */
ofloat *RLPhase;
/** Rotation measures of calibrators (rad/m^2) */
ofloat *RM;
/** Fractional linear polarization per calibrator */
ofloat *PPol;
/** Flux density of calibrators at reference freq */
ofloat *IFlux0;
/** Spectral index of calibrators at reference freq */
ofloat *IFlux1;
/** Curvature (in ln(nu)) of calibrators */
ofloat *IFlux2;
/** Array of input data arrays, [10 x nvis]
    each row: ant 1, 2 2*parallactic angle (rad),
    RR_r, RR_i, LL_r, LL_i, RL_r, RL_i,  LR_r, LR_i   */
ofloat *inData;
/** Array of input data weights (~ 1/var), [4 x nvis]
    each row: RR, LL, RL,  LR   */
ofloat *inWt;
/** Source no per visibility, in range [0,nsou-1] */
olong *souNo;
/** Antenna numbers (pair) per visibility, in range [0,nant-1] */
olong *antNo;
/** Number of parameters being fitted */
olong nparam;
/** Number of data poiints being fitted */
olong ndata;
/** Number of visibilities being fitted */
olong nvis;
/** List of source structures if multiple sources */
ObitSourceList *SouList;
/** Single source */
ObitSource *oneSource;
/** Antenna Lists */
ObitAntennaList **AntLists;
/** Number of subarrays */
olong numSubA;
/** selected antenna, 0-> all */
olong selAnt;
/** selected source, 0-> all */
olong selSou;
/** Source poln table */
ObitTableCP *CPTable;
/** Antenna poln table */
ObitTablePD *PDTable;
/** Bandpass table for gain corrections */
ObitTableBP *BPTable;
/** Print Level */
olong prtLv;
/** Max. number of Threads */
olong nThread;
/** Number of Threads used */
olong nTh;
/** Thread argument array */
PolnFitArg **thArgs;
/** complex work arrays for antenna based parameters, max 100 ant */
dcomplex RS[100], RD[100], LS[100], LD[100], 
  RSc[100], RDc[100], LSc[100], LDc[100],
  PR[100], PRc[100], PL[100], PLc[100];
odouble SR[100], DR[100], SL[100], DL[100];
/** Max antenna number */
olong maxAnt;
/** Solution method */
gchar solnType[5];
#ifdef HAVE_GSL
  /** Fitting solver */
  gsl_multifit_fdfsolver *solver;
  /** Fitting solver function structure */
  gsl_multifit_function_fdf *funcStruc;
  /** Fitting work vector */
  gsl_vector *work;
  /** Covariance matrix  */
  gsl_matrix *covar;
#endif /* HAVE_GSL */ 
