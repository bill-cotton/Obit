/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012-2013                                          */
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

#include "ObitPolnCalFit.h"
#include "ObitSpectrumFit.h"
#include "ObitThread.h"
#include "ObitTableAN.h"
#include "ObitTableANUtil.h"
#include "ObitTableSU.h"
#include "ObitTableSUUtil.h"
#include "ObitPrecess.h"
#ifdef HAVE_GSL
#include <gsl/gsl_blas.h>
#endif /* HAVE_GSL */ 
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPolnCalFit.c
 * ObitPolnCalFit class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitPolnCalFit";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitPolnCalFitClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitPolnCalFitClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitPolnCalFitInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitPolnCalFitClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitPolnCalFitClassInfoDefFn (gpointer inClass);

/** Private: Write output image */
static void WriteOutput (ObitPolnCalFit* in, ObitUV *outUV, 
			 gboolean isOK, ObitErr *err);

/** Private: Read data to thread arguments, init soln */
static void  ReadData(ObitPolnCalFit *in, ObitUV *inUV, olong *iChan, 
		      olong EChan, olong *iIF, olong EIF, 
		      gboolean first,  ObitErr *err);

/* Private: Initialize new Source poln.(CP) table */
static void InitSourceTab(ObitPolnCalFit *in, ObitErr *err);

/* Private: Initialize new Instrumental poln.(PD) table */
static void InitInstrumentalTab(ObitPolnCalFit *in, ObitErr *err);

/* Private: Initialize new Bandpass.(BP) table */
static void InitBandpassTab(ObitPolnCalFit *in, ObitErr *err);

/* Private: Update Source poln.(CP) table */
static void UpdateSourceTab(ObitPolnCalFit *in, ObitErr *err);

/* Private:  Update Instrumental poln.(PD) table */
static void UpdateInstrumentalTab(ObitPolnCalFit *in, gboolean isOK, ObitErr *err);

/* Private: Update Bandpass.(BP) table */
static void UpdateBandpassTab(ObitPolnCalFit *in, gboolean isOK, ObitErr *err);

/** Private: Fit spectra in SU list */
static void FitSpectra (ObitPolnCalFit *in, ObitErr *err);

/** Private: Fit data */
static gboolean doFitFast(ObitPolnCalFit *in, ObitErr *err);
static void doFitGSL (ObitPolnCalFit *in, ObitErr *err);

/** Private: Get Chi Sq & derivatives of current model */
static odouble GetChi2 (olong nThreads, ObitPolnCalFit *in, 
			PolnParmType paramType, olong paramNumber,
			ofloat *ParRMS, ofloat *XRMS, 
			odouble *dChi2, odouble *d2Chi2,
			ObitErr *err);

/** Private: Fit R/L Phase difference */
static ofloat FitRLPhase (ObitPolnCalFit *in, ObitErr *err);

/** Private: Make Threaded args */
static void MakePolnFitFuncArgs (ObitPolnCalFit *in, ObitErr *err);

/** Private: Delete Threaded args */
static void KillPolnFitFuncArgs (ObitPolnCalFit *in);

/** Private: Threaded Chi**2 R/L evaluator */
static gpointer ThreadPolnFitRLChi2 (gpointer arg);

/** Private: Threaded Chi**2 X/Y evaluator */
static gpointer ThreadPolnFitXYChi2 (gpointer arg);

/** Private: Check for crazy antenna solutions */
static gboolean CheckCrazy(ObitPolnCalFit *in, ObitErr *err);

#ifdef HAVE_GSL
/** Private: Circular feed Solver function evaluation */
static int PolnFitFuncOERL (const gsl_vector *x, void *params, 
			    gsl_vector *f);

/** Private: Circular feed Solver Jacobian evaluation */
static int PolnFitJacOERL (const gsl_vector *x, void *params, 
			   gsl_matrix *J);

/** Private: Circular feed Solver function + Jacobian evaluation */
static int PolnFitFuncJacOERL (const gsl_vector *x, void *params, 
			       gsl_vector *f, gsl_matrix *J);

/** Private: Linear feed Solver function evaluation */
static int PolnFitFuncOEXY (const gsl_vector *x, void *params, 
			    gsl_vector *f);

/** Private: Linear feed Solver Jacobian evaluation */
static int PolnFitJacOEXY (const gsl_vector *x, void *params, 
			   gsl_matrix *J);

/** Private: Linear feed Solver function + Jacobian evaluation */
static int PolnFitFuncJacOEXY (const gsl_vector *x, void *params, 
			       gsl_vector *f, gsl_matrix *J);

#endif /* HAVE_GSL */ 

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitPolnCalFit* newObitPolnCalFit (gchar* name)
{
  ObitPolnCalFit* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPolnCalFitClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitPolnCalFit));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitPolnCalFitInit((gpointer)out);

 return out;
} /* end newObitPolnCalFit */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitPolnCalFitGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPolnCalFitClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitPolnCalFitGetClass */

/**
 * Make a deep copy of an ObitPolnCalFit.
 * NOT YET IMPLEMENTED
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitPolnCalFit* ObitPolnCalFitCopy  (ObitPolnCalFit *in, ObitPolnCalFit *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitPolnCalFit(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */

  /* Arrays */
    
  /* reference this class members */
  return out;
} /* end ObitPolnCalFitCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an PolnCalFit similar to the input one.
 * NOT YET IMPLEMENTED
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitPolnCalFitClone  (ObitPolnCalFit *in, ObitPolnCalFit *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  
  /* Arrays */

  /* reference this class members */
} /* end ObitPolnCalFitClone */

/**
 * Creates an ObitPolnCalFit 
 * \param name   An optional name for the object.
 * \return the new object.
 */
ObitPolnCalFit* ObitPolnCalFitCreate (gchar* name)
{
  ObitPolnCalFit* out;

  /* Create basic structure */
  out = newObitPolnCalFit (name);

  return out;
} /* end ObitPolnCalFitCreate */

/**
 * Calculate circular feed model using a thread argument
 * \param args   argument
 * \param Rarray output array
 * \param idata  0-rel data number
 */
static void calcmodelRL (ObitPolnCalFit *args, ofloat Rarray[8], olong idata)
{
  odouble    *antParm   = args->antParm;
  odouble    *souParm   = args->souParm;
  ofloat     *data      = args->inData;
  dcomplex   *RS        = args->RS;
  dcomplex   *RD        = args->RD;
  dcomplex   *LS        = args->LS;
  dcomplex   *LD        = args->LD;
  dcomplex   *RSc       = args->RSc;
  dcomplex   *RDc       = args->RDc;
  dcomplex   *LSc       = args->LSc;
  dcomplex   *LDc       = args->LDc;

  odouble ipol=0.0, qpol=0.0, upol=0.0, vpol=0.0;
  ofloat  PD;
  olong i, ia1, ia2, isou;

  dcomplex PRref, PLref, PPRL, PPLR, PA1, PA2, PA1c, PA2c, ct1, ct2;
  dcomplex S[4], VRR, VRL, VLR, VLL, MC1, MC2, MC3, MC4;
  ofloat root2, chi1, chi2;

  /* R-L phase difference  at reference antenna */
  if (args->doFitRL) {
    PD = args->PD;
  } else PD = 0.0;

  COMPLEX_SET (MC1, 0.0, 0.0);  /* Other stuff */
  COMPLEX_SET (MC2, 0.0, 0.0);
  COMPLEX_SET (MC3, 0.0, 0.0);
  COMPLEX_SET (MC4, 0.0, 0.0);
  chi1  = data[idata*10+0];   /* parallactic angle ant 1 */
  chi2  = data[idata*10+1];   /* parallactic angle ant 2 */
  COMPLEX_EXP (PA1, 2*chi1);
  COMPLEX_EXP (PA2, 2*chi2);
  COMPLEX_CONJUGATE (PA1c, PA1);
  COMPLEX_CONJUGATE (PA2c, PA2);
  
  isou  = MAX (0, args->souNo[idata]);    /* Source number */
  
  /* Source parameters */
  ipol = souParm[isou*4+0];
  qpol = souParm[isou*4+1];
  upol = souParm[isou*4+2];
  vpol = souParm[isou*4+3];
  /* Complex Stokes array */
  COMPLEX_SET (S[0], ipol+vpol, 0.0);
  COMPLEX_SET (S[1], qpol,  upol);
  COMPLEX_SET (S[2], qpol, -upol);
  COMPLEX_SET (S[3], ipol-vpol, 0.0);
  
  /* Injest model factorize into antenna components - 
     data in order Orientation R/X, Elipticity R/X, Orientation L/Y, Elipticity L/Y*/
  root2 = 1.0 / sqrt(2.0);
  /* Elipticity, Orientation terms */
  for (i=0; i<args->nant; i++) {
    COMPLEX_SET(RS[i], root2*(cos(antParm[i*4+1]) + sin(antParm[i*4+1])), 0);
    COMPLEX_SET(ct1,   root2*(cos(antParm[i*4+1]) - sin(antParm[i*4+1])), 0);
    COMPLEX_EXP(ct2, 2*antParm[i*4+0]);
    COMPLEX_MUL2 (RD[i], ct1, ct2);
    COMPLEX_SET (ct1, root2*(cos(antParm[i*4+3]) + sin(antParm[i*4+3])), 0);
    COMPLEX_EXP (ct2, -2*antParm[i*4+2]);
    COMPLEX_MUL2 (LS[i], ct1, ct2);
    COMPLEX_SET (LD[i], root2*(cos(antParm[i*4+3]) - sin(antParm[i*4+3])), 0);
    COMPLEX_CONJUGATE (RSc[i], RS[i]);
    COMPLEX_CONJUGATE (RDc[i], RD[i]);
    COMPLEX_CONJUGATE (LSc[i], LS[i]);
    COMPLEX_CONJUGATE (LDc[i], LD[i]);
  }

  /* Reference antenna phase terms */
  COMPLEX_EXP (PRref,  antParm[(args->refAnt-1)*4+0]);
  COMPLEX_EXP (PLref, -antParm[(args->refAnt-1)*4+2]+PD);
  COMPLEX_CONJUGATE (ct1, PLref);
  COMPLEX_MUL2 (PPRL, PRref, ct1);
  COMPLEX_CONJUGATE (ct1, PRref);
  COMPLEX_MUL2 (PPLR, PLref, ct1);
  
  ia1    = args->antNo[idata*2+0];
  ia2    = args->antNo[idata*2+1]; 
  
  /* VRR = S[0] * RS[ia1] * RSc[ia2] +        
           S[1] * RS[ia1] * RDc[ia2] * PA2c + 
	   S[2] * RD[ia1] * RSc[ia2] * PA1  + 
	   S[3] * RD[ia1] * RDc[ia2] * PA1  * PA2c; */
  COMPLEX_MUL2 (MC1, RS[ia1], RSc[ia2]);
  COMPLEX_MUL2 (VRR, S[0], MC1);
  COMPLEX_MUL3 (MC2, RS[ia1], RDc[ia2],  PA2c);
  COMPLEX_MUL2 (ct1, S[1], MC2);
  COMPLEX_ADD2 (VRR, VRR,  ct1);
  COMPLEX_MUL3 (MC3, RD[ia1], RSc[ia2], PA1);
  COMPLEX_MUL2 (ct1, S[2], MC3);
  COMPLEX_ADD2 (VRR, VRR,  ct1);
  COMPLEX_MUL4 (MC4, RD[ia1], RDc[ia2], PA1, PA2c);
  COMPLEX_MUL2 (ct1, S[3], MC4);
  COMPLEX_ADD2 (VRR, VRR,  ct1);
  Rarray[0] = VRR.real;
  Rarray[1] = VRR.imag;

  /* VLL = S[0] * LS[ia1] * LSc[ia2] * PA1c * PA2 +	
           S[1] * LS[ia1] * LDc[ia2] * PA1c +
	   S[2] * LD[ia1] * LSc[ia2] * PA2  +
	   S[3] * LD[ia1] * LDc[ia2]; */
  COMPLEX_MUL4 (MC1, LS[ia1], LSc[ia2], PA1c, PA2);
  COMPLEX_MUL2 (VLL, S[0], MC1);
  COMPLEX_MUL3 (MC2, LS[ia1], LDc[ia2], PA1c);
  COMPLEX_MUL2 (ct1, S[1], MC2);
  COMPLEX_ADD2 (VLL, VLL,  ct1);
  COMPLEX_MUL3 (MC3, LD[ia1], LSc[ia2], PA2);
  COMPLEX_MUL2 (ct1, S[2], MC3);
  COMPLEX_ADD2 (VLL, VLL,  ct1);
  COMPLEX_MUL2 (MC4, LD[ia1], LDc[ia2]);
  COMPLEX_MUL2 (ct1, S[3], MC4);
  COMPLEX_ADD2 (VLL, VLL,  ct1);
  Rarray[6] = VLL.real;
  Rarray[7] = VLL.imag;

  /* 	    RL */
  /* VRL = PPRL * S[0] * RS[ia1] * LSc[ia2] * PA2 +
           PPRL * S[1] * RS[ia1] * LDc[ia2] + 
	   PPRL * S[2] * RD[ia1] * LSc[ia2] * PA1 * PA2 +
	   PPRL * S[3] * RD[ia1] * LDc[ia2] * PA1; */
  COMPLEX_MUL4 (MC1, PPRL, RS[ia1], LSc[ia2], PA2);
  COMPLEX_MUL2 (VRL, S[0], MC1);
  COMPLEX_MUL3 (MC2, PPRL, RS[ia1], LDc[ia2]);
  COMPLEX_MUL2 (ct1, S[1], MC2);
  COMPLEX_ADD2 (VRL, VRL,  ct1);
  COMPLEX_MUL5 (MC3, PPRL, RD[ia1], LSc[ia2],  PA1,  PA2);
  COMPLEX_MUL2 (ct1, S[2], MC3);
  COMPLEX_ADD2 (VRL, VRL,  ct1);
  COMPLEX_MUL4 (MC4, PPRL, RD[ia1], LDc[ia2],  PA1);
  COMPLEX_MUL2 (ct1, S[3], MC4);
  COMPLEX_ADD2 (VRL, VRL,  ct1);
  Rarray[2] = VRL.real;
  Rarray[3] = VRL.imag;

  /*        LR */
  /* VLR = PPLR * S[0] * LS[ia1] * RSc[ia2] * PA1c +
           PPLR * S[1] * LS[ia1] * RDc[ia2] * PA1c * PA2c +
	   PPLR * S[2] * LD[ia1] * RSc[ia2] +
	   PPLR * S[3] * LD[ia1] * RDc[ia2] * PA2c */
  COMPLEX_MUL4 (MC1, PPLR, LS[ia1], RSc[ia2], PA1c);
  COMPLEX_MUL2 (VLR, S[0], MC1);
  COMPLEX_MUL5 (MC2, PPLR, LS[ia1], RDc[ia2], PA1c,  PA2c);
  COMPLEX_MUL2 (ct1, S[1], MC2);
  COMPLEX_ADD2 (VLR, VLR,  ct1);
  COMPLEX_MUL3 (MC3, PPLR, LD[ia1], RSc[ia2]);
  COMPLEX_MUL2 (ct1, S[2], MC3);
  COMPLEX_ADD2 (VLR, VLR,  ct1);
  COMPLEX_MUL4 (MC4, PPLR, LD[ia1], RDc[ia2],  PA2c);
  COMPLEX_MUL2 (ct1, S[3], MC4);
  COMPLEX_ADD2 (VLR, VLR, ct1);
  Rarray[4] = VLR.real;
  Rarray[5] = VLR.imag;
} /* end calcmodelRL */

/**
 * Calculate linear feed model using a thread argument
 * R Perley version
 * \param args   argument
 * \param Rarray output array
 * \param idata  0-rel data number
 */
static void calcmodelXY (ObitPolnCalFit *args, ofloat Rarray[8], olong idata)
{
  odouble    *antParm   = args->antParm;
  odouble    *antGain   = args->antGain;
  odouble    *souParm   = args->souParm;
  ofloat     *data      = args->inData;
  dcomplex   *CX        = args->RS;
  dcomplex   *SX        = args->RD;
  dcomplex   *CY        = args->LS;
  dcomplex   *SY        = args->LD;
  dcomplex   *CXc       = args->RSc;
  dcomplex   *SXc       = args->RDc;
  dcomplex   *CYc       = args->LSc;
  dcomplex   *SYc       = args->LDc;

  odouble ipol=0.0, qpol=0.0, upol=0.0, vpol=0.0;
  olong i, ia1, ia2, isou;

  dcomplex Jp, Jm, SPA, DPA, SPAc, DPAc, ct1, ct2;
  dcomplex S[4], VXX, VXY, VYX, VYY, MC1, MC2, MC3, MC4;
  dcomplex SM1, SM2, SM3, SM4, ggPD;
  ofloat root2, chi1, chi2, PD;

  /* Init working variables */
  COMPLEX_SET (MC1, 0.0, 0.0);  /* Muller matrix */
  COMPLEX_SET (MC2, 0.0, 0.0);
  COMPLEX_SET (MC3, 0.0, 0.0);
  COMPLEX_SET (MC4, 0.0, 0.0);
  COMPLEX_SET (Jp,  0.0, 1.0);
  COMPLEX_SET (Jm,  0.0,-1.0);
  chi1  = data[idata*10+0];   /* parallactic angle ant 1 */
  chi2  = data[idata*10+1];   /* parallactic angle ant 2 */
  COMPLEX_EXP (SPA,chi1+chi2);
  /*COMPLEX_SET (SPA, 1.0, 0.0); DEBUG */
  COMPLEX_EXP (DPA,chi1-chi2);
  /*COMPLEX_SET (DPA, 1.0, 0.0);  DEBUG */
  COMPLEX_CONJUGATE (SPAc, SPA);
  COMPLEX_CONJUGATE (DPAc, DPA);
  
  isou  = MAX (0, args->souNo[idata]);    /* Source number */
  
   /* X-Y phase difference  at reference antenna */
  if (args->doFitRL) {
    PD = args->PD;
  } else PD = 0.0;

 /* Source parameters */
  ipol = souParm[isou*4+0];
  qpol = souParm[isou*4+1];
  upol = souParm[isou*4+2];
  vpol = souParm[isou*4+3];
  /* Complex Stokes array as R/L correlations */
  COMPLEX_SET (S[0], ipol+vpol, 0.0);
  COMPLEX_SET (S[1], qpol,  upol);
  COMPLEX_SET (S[2], qpol, -upol);
  COMPLEX_SET (S[3], ipol-vpol, 0.0);
  
  /* Injest model factorize into antenna components - 
     data in order 0: Orientation X, 1: Elipticity X, 2: Orientation Y, 3: Elipticity Y */
  root2 = 1.0 / sqrt(2.0);
  /* Elipticity, Orientation terms */
  for (i=0; i<args->nant; i++) {
    COMPLEX_EXP (ct1, -antParm[i*4+0]);
    /*COMPLEX_EXP (ct1, -antParm[i*4+0]-chi1); DEBUG */
    COMPLEX_SET (ct2, cos(G_PI*0.25+antParm[i*4+1]), 0.0);
    COMPLEX_MUL3 (CX[i], Jp, ct1, ct2);
    COMPLEX_EXP (ct1, antParm[i*4+0]);
    /* COMPLEX_EXP (ct1, antParm[i*4+0]+chi1);  DEBUG*/
    COMPLEX_SET (ct2, sin(G_PI*0.25+antParm[i*4+1]), 0.0);
    COMPLEX_MUL2 (SX[i], ct1, ct2);
    COMPLEX_EXP (ct1, antParm[i*4+2]);
    /*COMPLEX_EXP (ct1, antParm[i*4+2]+chi1); DEBUG */
    COMPLEX_SET (ct2, cos(G_PI*0.25-antParm[i*4+3]), 0.0);
    COMPLEX_MUL2 (CY[i], ct1, ct2);
    COMPLEX_EXP (ct1, -antParm[i*4+2]);
    /*COMPLEX_EXP (ct1, -antParm[i*4+2]-chi1); DEBUG */
    COMPLEX_SET (ct2, sin(G_PI*0.25-antParm[i*4+3]), 0.0);
    COMPLEX_MUL3 (SY[i], Jp, ct1, ct2);
    COMPLEX_CONJUGATE (CXc[i], CX[i]);
    COMPLEX_CONJUGATE (SXc[i], SX[i]);
    COMPLEX_CONJUGATE (CYc[i], CY[i]);
    COMPLEX_CONJUGATE (SYc[i], SY[i]);
  }

  ia1    = args->antNo[idata*2+0];
  ia2    = args->antNo[idata*2+1]; 
  
  /* VXX = {S[0] * CX[ia1] * CXc[ia2] * DPAc  +
            S[1] * CX[ia1] * SXc[ia2] * SPAc  +
	    S[2] * SX[ia1] * CXc[ia2] * SPA   + 
	    S[3] * SX[ia1] * SXc[ia2] * DPA} * g1X * g2X ;
  */
  COMPLEX_MUL3 (MC1, CX[ia1], CXc[ia2], DPAc);
  COMPLEX_MUL3 (MC2, CX[ia1], SXc[ia2], SPAc);
  COMPLEX_MUL3 (MC3, SX[ia1], CXc[ia2], SPA);
  COMPLEX_MUL3 (MC4, SX[ia1], SXc[ia2], DPA);
  COMPLEX_MUL2 (VXX, S[0], MC1);
  COMPLEX_MUL2 (ct1, S[1], MC2);
  COMPLEX_ADD2 (VXX, VXX,  ct1);
  COMPLEX_MUL2 (ct1, S[2], MC3);
  COMPLEX_ADD2 (VXX, VXX,  ct1);
  COMPLEX_MUL2 (ct1, S[3], MC4);
  COMPLEX_ADD2 (VXX, VXX,  ct1);
  COMPLEX_SET (ct1,  antGain[ia1*2+0]*antGain[ia2*2+0], 0);
  COMPLEX_MUL2 (VXX, VXX, ct1);
  Rarray[0] = VXX.real;
  Rarray[1] = VXX.imag;

  /* VYY = {S[0] * SY[ia1] * SYc[ia2] * DPAc +       
            S[1] * SY[ia1] * CYc[ia2] * SPAc +
	    S[2] * CY[ia1] * SYc[ia2] * SPA  + 
	    S[3] * CY[ia1] * CYc[ia2] * DPA} * g1Y * g2Y ;
  */
  COMPLEX_MUL3 (MC1, SY[ia1], SYc[ia2], DPAc);
  COMPLEX_MUL3 (MC2, SY[ia1], CYc[ia2], SPAc);
  COMPLEX_MUL3 (MC3, CY[ia1], SYc[ia2], SPA);
  COMPLEX_MUL3 (MC4, CY[ia1], CYc[ia2], DPA);
  COMPLEX_MUL2 (VYY, S[0], MC1);
  COMPLEX_MUL2 (ct1, S[1], MC2);
  COMPLEX_ADD2 (VYY, VYY,  ct1);
  COMPLEX_MUL2 (ct1, S[2], MC3);
  COMPLEX_ADD2 (VYY, VYY,  ct1);
  COMPLEX_MUL2 (ct1, S[3], MC4);
  COMPLEX_ADD2 (VYY, VYY,  ct1);
  COMPLEX_SET (ct1,  antGain[ia1*2+1]*antGain[ia2*2+1], 0);
  COMPLEX_MUL2 (VYY, VYY, ct1);
  Rarray[6] = VYY.real;
  Rarray[7] = VYY.imag;

  /* VXY = {S[0] * CX[ia1] * SYc[ia2] * DPAc +       
            S[1] * CX[ia1] * CYc[ia2] * SPAc +
	    S[2] * SX[ia1] * SYc[ia2] * SPA  + 
	    S[3] * SX[ia1] * CYc[ia2] * DPA}} * g1X * g2Y * exp(i PD);
  */
  COMPLEX_MUL3 (MC1, CX[ia1], SYc[ia2], DPAc);
  COMPLEX_MUL3 (MC2, CX[ia1], CYc[ia2], SPAc);
  COMPLEX_MUL3 (MC3, SX[ia1], SYc[ia2], SPA);
  COMPLEX_MUL3 (MC4, SX[ia1], CYc[ia2], DPA);
  COMPLEX_MUL2 (SM1, S[0], MC1);
  COMPLEX_MUL2 (SM2, S[1], MC2);
  COMPLEX_MUL2 (SM3, S[2], MC3);
  COMPLEX_MUL2 (SM4, S[3], MC4);
  COMPLEX_ADD4 (VXY, SM1, SM2, SM3, SM4);
  COMPLEX_SET (ct1,  antGain[ia1*2+0]*antGain[ia2*2+1], 0);
  COMPLEX_EXP (ct2, PD);
  COMPLEX_MUL2 (ggPD, ct1, ct2);
  COMPLEX_MUL2 (VXY, VXY, ggPD);
  Rarray[2] = VXY.real;
  Rarray[3] = VXY.imag;

  /* VYX = {S[0] * SY[ia1] * CXc[ia2] * DPAc +       
            S[1] * SY[ia1] * SXc[ia2] * SPAc +
	    S[2] * CY[ia1] * CXc[ia2] * SPA  + 
	    S[3] * CY[ia1] * SXc[ia2] * DPA} * g1Y * g2X * exp(-i PD);
  */
  COMPLEX_MUL3 (MC1, SY[ia1], CXc[ia2], DPAc);
  COMPLEX_MUL3 (MC2, SY[ia1], SXc[ia2], SPAc);
  COMPLEX_MUL3 (MC3, CY[ia1], CXc[ia2], SPA);
  COMPLEX_MUL3 (MC4, CY[ia1], SXc[ia2], DPA);
  COMPLEX_MUL2 (SM1, S[0], MC1);
  COMPLEX_MUL2 (SM2, S[1], MC2);
  COMPLEX_MUL2 (SM3, S[2], MC3);
  COMPLEX_MUL2 (SM4, S[3], MC4);
  COMPLEX_ADD4 (VYX, SM1, SM2, SM3, SM4);
  COMPLEX_SET (ct1,  antGain[ia1*2+1]*antGain[ia2*2+0], 0);
  COMPLEX_EXP (ct2, -PD);
  COMPLEX_MUL2 (ggPD, ct1, ct2);
  COMPLEX_MUL2 (VYX, VYX, ggPD);
  Rarray[4] = VYX.real;
  Rarray[5] = VYX.imag;
  /* DEBUG 
      fprintf (stdout, "model %4d %d %d XY M %8.3f %8.3f SM %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ggPD %8.3f %8.3f \n",
	       idata,ia1, ia2, VYX.real, VYX.imag, 
               SM1.real, SM1.imag, SM2.real, SM2.imag, SM3.real, SM3.imag, SM4.real, SM4.imag, ggPD.real, ggPD.imag);
    end DEBUG */
} /* end calcmodelXY */

/**
 * Fit selected parameters to a UV data.
 * \param in       Spectral fitting object
 *                 Potential parameters on in->info:
 * \li nCal       OBIT_int scalar Number of calibrator sources [1]
 * \li solnType   OBIT_string [4] Solution type:'FULD'
 * \li Sources    OBIT_string [16,*] Calibrator sources 
 * \li Qual       OBIT_long scalar Source qualifier -1=>any
 * \li souCode    OBIT_string [4] Calibrator code '    '=>any
 * \li doFitI     OBIT_boolean [*] If True fit Stokes I poln, per cal [F]
 * \li doFitPol   OBIT_boolean [*] If True fit source linear poln, per cal  [T]
 * \li doFitV     OBIT_boolean [*] If True fit Stokes V poln, per cal  [F]
 * \li doFitRL    OBIT_boolean  If True fit R-L phase difference  [F]
 * \li RLPhase    OBIT_float scalar R-L phase difference @ ref Freq per calibrator
 *                <= -999 => fit value [-1000.] in deg.
 * \li RM         OBIT_float scalar Rotation measure (rad/m^2) to be used with RLPhase
 * \li PPol       OBIT_float scalar Fractional linear polarization per calibrator
 * \li ChWid      OBIT_long scalar Number of channels in poln soln  [1]
 * \li ChInc      OBIT_long scalar  Spacing of poln soln  [1]
 * \li CPSoln     OBIT_long scalar Source (CP) table to write, 0=> create new  [0]
 * \li PDSoln     OBIT_long scalar Instrumental (PD) table to write, 0=> create new  [0]
 * \li BPSoln     OBIT_long scalar Bandpass (BP) table to write, 0=> create new  [0]
 * \li doBlank    OBIT_boolean scalar If True blanked failed solns, else default values [T]
 * \li doBand     OBIT_long scalar If >0 then BP table BPVer was applied to input  [0]
 * \li BPVer      OBIT_long scalar Input BP table  [0]
 * \li refAnt     OBIT_long scalar Reference antenna  [1]
 * \li BIF        OBIT_long scalar First IF (1-rel) selected in outUV  [1]
 * \li BChan      OBIT_long scalar First channel (1-rel) selected in outUV  [1]
 * \param inUV  Averaged/calibrated/divided UV data to be fitted
 *              Should be averaged to the solInt time.
 *              MUST be a multisource file
 * \param outUV UV data with fitted results in tables
 * \param err      Obit error stack object.
 */
void ObitPolnCalFitFit (ObitPolnCalFit* in, ObitUV *inUV, 
			ObitUV *outUV, ObitErr *err)
{
  olong nCal=1;
  ObitIOCode retCode;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1}, sdim[MAXINFOELEMDIM];
  ObitInfoType type;
  ObitTableAN    *ANTable = NULL;
  ObitTableSU    *SUTable = NULL;
  olong i, j, k, iant, isou, numAnt, BChan, EChan, iChan, BIF, EIF, iIF;
  olong nchleft, first, highANver, iANver, ver, numIF, Qual, suba=1, nvis;
  olong highSUver, size, oldNparam, oldNdata, iarr[5]={0,0,0,0,0};
  ofloat *farr;
  olong number, maxSUId=-1;
  gboolean done, xselect=TRUE, *barr, isOK;
  gchar solnType[5], *Sources, souCode[5];
  ofloat  endChi2, ParRMS, XRMS;
#ifdef HAVE_GSL
  const gsl_multifit_fdfsolver_type *T=NULL;
#endif /* HAVE_GSL */ 
  /* DEBUG */
  ofloat RArray[8];
  odouble RRr, RRi, LLr, LLi, RLr, RLi, LRr, LRi;
  olong ia1, ia2, refAnt;
  gchar *routine = "ObitPolnCalFitFit";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert(ObitUVIsA(inUV));
  g_assert(ObitUVIsA(outUV));
  /* Multisource file? */
  Obit_return_if_fail (( inUV->myDesc->ilocsu>=0), err, 
		       "%s Input %s NOT a multisource file",  
		       routine, inUV->name);  

  /* How many calibrators? */
  ObitInfoListGetTest(in->info, "nCal", &type, dim, &nCal);

  /* Allocate calibrator arrays */
  in->RLPhaseIn= g_malloc0(nCal*sizeof(ofloat));
  in->RLPhase  = g_malloc0(nCal*sizeof(ofloat));
  in->RM       = g_malloc0(nCal*sizeof(ofloat));
  in->PPol     = g_malloc0(nCal*sizeof(ofloat));
  in->doFitI   = g_malloc0(nCal*sizeof(olong));
  in->doFitV   = g_malloc0(nCal*sizeof(olong));
  in->doFitPol = g_malloc0(nCal*sizeof(olong));
  in->IFlux0   = g_malloc0(nCal*sizeof(ofloat));
  in->IFlux1   = g_malloc0(nCal*sizeof(ofloat));
  in->IFlux2   = g_malloc0(nCal*sizeof(ofloat));

  /* Control parameters - arrays may be defined larger than nCal */
  if (ObitInfoListGetP(in->info, "RLPhase", &type, dim, (gpointer)&farr))
    for (i=0; i<nCal; i++) in->RLPhaseIn[i] = farr[i] * DG2RAD;
  else
    for (i=0; i<nCal; i++) in->RLPhaseIn[i] =- 1000.0;
  if (ObitInfoListGetP(in->info, "RM", &type, dim, (gpointer)&farr))
    for (i=0; i<nCal; i++) in->RM[i] = farr[i];
  else
    for (i=0; i<nCal; i++)  in->RM[i] = 0.0;
  if (ObitInfoListGetP(in->info, "PPol", &type, dim, (gpointer)&farr))
    for (i=0; i<nCal; i++) in->PPol[i] = farr[i];
  else
    for (i=0; i<nCal; i++)  in->PPol[i] = 0.0;
  if (ObitInfoListGetP(in->info, "doFitI", &type, dim, (gpointer)&barr))
    for (i=0; i<nCal; i++) in->doFitI[i] = barr[i];
  else
    for (i=0; i<nCal; i++) in->doFitI[i] = FALSE;
  if (ObitInfoListGetP(in->info, "doFitV", &type, dim, (gpointer)&barr))
    for (i=0; i<nCal; i++) in->doFitV[i] = barr[i];
  else
    for (i=0; i<nCal; i++) in->doFitV[i] = FALSE;
  if (ObitInfoListGetP(in->info, "doFitPol", &type, dim, (gpointer)&barr))
    for (i=0; i<nCal; i++) in->doFitPol[i] = barr[i];
  else
    for (i=0; i<nCal; i++) in->doFitPol[i] = TRUE;
  ObitInfoListGetTest(in->info, "doFitRL",  &type, dim, &in->doFitRL);
  ObitInfoListGetTest(in->info, "doBlank",  &type, dim, &in->doBlank);
  ObitInfoListGetTest(in->info, "doFitGain",&type, dim, &in->doFitGain);
  ObitInfoListGetTest(in->info, "ChWid",    &type, dim, &in->ChWid);
  ObitInfoListGetTest(in->info, "ChInc",    &type, dim, &in->ChInc);
  ObitInfoListGetTest(in->info, "CPSoln",   &type, dim, &in->CPSoln);
  ObitInfoListGetTest(in->info, "PDSoln",   &type, dim, &in->PDSoln);
  ObitInfoListGetTest(in->info, "BPSoln",   &type, dim, &in->BPSoln);
  ObitInfoListGetTest(in->info, "doBand",   &type, dim, &in->doBand);
  ObitInfoListGetTest(in->info, "BPVer",    &type, dim, &in->BPVer);
  ObitInfoListGetTest(in->info, "BIF",      &type, dim, &in->BIF);
  ObitInfoListGetTest(in->info, "BChan",    &type, dim, &in->BChan);
  ObitInfoListGetTest(in->info, "refAnt",   &type, dim, &in->refAnt);
  ObitInfoListGetTest(in->info, "prtLv",    &type, dim, &in->prtLv);
  ObitInfoListGetTest(in->info, "Qual",     &type, dim, &Qual);
  ObitInfoListGetTest(in->info, "souCode",  &type, dim, souCode);
  strcpy (solnType, "LM  ");
  ObitInfoListGetTest(in->info, "solnType", &type, dim, solnType);
  strncpy (in->solnType, solnType, 4);
  ObitInfoListGetP(in->info, "Sources", &type, sdim, (gpointer)&Sources);

  /* Open input data to get info */
  retCode = ObitUVOpen (inUV, OBIT_IO_ReadOnly, err);
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inUV->name);

  /* Save in/output descriptor */
  in->inDesc  = ObitUVDescRef(inUV->myDesc);
  in->outDesc = ObitUVDescRef(outUV->myDesc);

  /* Is this circular (linear) feed data? */
  in->isCircFeed = 
    (in->outDesc->crval[in->outDesc->jlocs]<0.0) && 
    (in->outDesc->crval[in->outDesc->jlocs]>-1.5);

  /* Gain fitting only for linear feeeds */
  if (in->isCircFeed) in->doFitGain = FALSE;

  /* Number of antennas */
  numAnt  = inUV->myDesc->numAnt[suba-1];/* actually highest antenna number */
  in->nant = numAnt;
  in->nsou = nCal;

  /* How many AN tables (subarrays) */
  highANver = ObitTableListGetHigh (inUV->tableList, "AIPS AN");
  
  /* Antenna lists */
  in->AntLists = g_malloc0(highANver*sizeof(ObitAntennaList*));
  in->numSubA  = highANver;
 /* Read Info from AN tables  */
  for (i=0; i<highANver; i++) {
    iANver = i+1;
    /* Get table */
    ANTable = newObitTableANValue (inUV->name, (ObitData*)inUV, &iANver, 
				   OBIT_IO_ReadOnly, 0, 0, 0, err);
    if (err->error) Obit_traceback_msg (err, routine, inUV->name);
    
    in->AntLists[i] = ObitTableANGetList (ANTable, err);
    if (err->error) Obit_traceback_msg (err, routine, inUV->name);

    /* release table object */
    ANTable = ObitTableANUnref(ANTable);
  } /* End loop over subarrays */

  /* Allocate source lookup table */
  in->souIDs = g_malloc0(in->nsou*sizeof(olong));

  /* Convert SU table into Source List if there is an SU table */
  highSUver = ObitTableListGetHigh (outUV->tableList, "AIPS AN");
  if (highSUver>0) {
    ver   = 1;
    numIF = 0;
    SUTable = newObitTableSUValue (outUV->name, (ObitData*)outUV, &ver, numIF, 
				   OBIT_IO_ReadOnly, err);
    if (SUTable) in->SouList = ObitTableSUGetList (SUTable, err);
    /* Look up source IDs - one at a time to get correct order */
    number = 1; sdim[1] = 1;
    for (i=0; i<in->nsou; i++) {
      j = sdim[0]*i;
      ObitTableSULookup (SUTable, sdim, &Sources[j], Qual, souCode, iarr, 
			 &xselect, &number, err);
      if (err->error) Obit_traceback_msg (err, routine, inUV->name);
      in->souIDs[i] = MAX (1, iarr[0]);
    } /* End loop over calibrators */
    /* Get maximum source id */
    for (i=0; i<in->nsou; i++) maxSUId = MAX(maxSUId, in->souIDs[i]);
    maxSUId = MAX (1, maxSUId);
    /* Inverse lookup table */
    in->isouIDs = g_malloc0(maxSUId*sizeof(olong));
    for (j=0; j<maxSUId; j++)  in->isouIDs[j] = -1000;
    for (i=0; i<in->nsou; i++) in->isouIDs[in->souIDs[i]-1] = i;

    /* Fit spectra to IPol in table */
    FitSpectra (in, err);
    if (err->error) Obit_traceback_msg (err, routine, inUV->name);
    SUTable = ObitTableSUUnref(SUTable);
    
  } else {  /* Create a single source object */
    in->oneSource = newObitSource("Single");
    strncpy (in->oneSource->SourceName, inUV->myDesc->object, MIN(20,UVLEN_VALUE));
    in->oneSource->equinox = inUV->myDesc->equinox;
    in->oneSource->RAMean  = inUV->myDesc->crval[inUV->myDesc->jlocr];
    SUTable = ObitTableSUUnref(SUTable);
    in->oneSource->DecMean = inUV->myDesc->crval[inUV->myDesc->jlocd];
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (inUV->myDesc, in->oneSource);
    in->souIDs[0] = 1;  /* Source ID */
    /* Inverse lookup table */
    in->isouIDs = g_malloc0(1*sizeof(olong));
    in->isouIDs[0] = 1;
  }

  /* Allocate Fitting flag/order arrays */
  in->gotAnt = g_malloc0(in->nant*sizeof(gboolean*));
  in->antFit = g_malloc0(in->nant*sizeof(gboolean*));
  for (i=0; i<in->nant; i++) in->antFit[i] = g_malloc0(4*sizeof(gboolean));
  in->souFit = g_malloc0(in->nsou*sizeof(gboolean*));
  for (i=0; i<in->nsou; i++) in->souFit[i] = g_malloc0(4*sizeof(gboolean));

  in->antPNumb = g_malloc0(in->nant*sizeof(gboolean*));
  for (i=0; i<in->nant; i++) in->antPNumb[i] = g_malloc0(4*sizeof(gboolean));
  in->souPNumb = g_malloc0(in->nsou*sizeof(gboolean*));
  for (i=0; i<in->nsou; i++) in->souPNumb[i] = g_malloc0(4*sizeof(gboolean));
  oldNparam = 0;
  oldNdata  = 0;

  /* Parameter/error arrays */
  in->antParm     = g_malloc0(in->nant*4*sizeof(odouble));
  in->antErr      = g_malloc0(in->nant*4*sizeof(odouble));
  in->souParm     = g_malloc0(in->nsou*4*sizeof(odouble));
  in->souErr      = g_malloc0(in->nsou*4*sizeof(odouble));
  in->lastSouParm = g_malloc0(in->nsou*4*sizeof(odouble));

  /* Antenna gains for linear feeds? */
  if (!in->isCircFeed ) {
    in->antGain     = g_malloc0((in->nant+2)*2*sizeof(odouble));
    for (i=0; i<in->nant*2; i++) in->antGain[i] = 1.0;
    in->antGainErr  = g_malloc0(in->nant*2*sizeof(odouble));
    in->antGainFit  = g_malloc0(in->nant*sizeof(gboolean*));
    for (i=0; i<in->nant; i++) {
      in->antGainFit[i] = g_malloc0(2*sizeof(gboolean));
      in->antGainFit[i][0] = in->antGainFit[i][1] = in->doFitGain;
    }
    in->antGainPNumb = g_malloc0(in->nant*sizeof(gboolean*));
    for (i=0; i<in->nant; i++) 
      in->antGainPNumb[i] = g_malloc0(2*sizeof(gboolean));

  } /* End create feed gain arrays */

  /* Fill Fitting flag arrays */
  /* Order of antenna parameters:
     0,1 ori, elip (r/x),
     2,3 ori, elip  (l/y),
   */
  for (i=0; i<in->nant; i++) {
    /* Initialize FALSE */
    in->gotAnt[i] = FALSE;

    /* Antenna terms */
    for (j=0; j<4; j++) in->antFit[i][j] = TRUE;
    /* Don't fit reference antenna terms for R */
    if  (((i+1)==in->refAnt) && in->isCircFeed) {
      in->antFit[i][0] = FALSE;
    }
  } /* end filling antenna fitting flags */

  /* Order of source parameters:
     0, Stokes I,
     1, Stokes Q,
     2, Stokes U
     3, Stokes V
   */
  for (i=0; i<in->nsou; i++) {
    if (in->doFitI[i])    in->souFit[i][0] = TRUE;
    else                  in->souFit[i][0] = FALSE;
    if (in->doFitPol[i]) {in->souFit[i][1] = TRUE;  in->souFit[i][2] = TRUE;}
    else                 {in->souFit[i][1] = FALSE; in->souFit[i][2] = FALSE;}
    if (in->doFitV[i])    in->souFit[i][3] = TRUE;
    else                  in->souFit[i][3] = FALSE;
    in->souParm[i*4+0] = in->souParm[i*4+1] = in->souParm[i*4+2] = in->souParm[i*4+3] = 0.0;
    /* Last valid source parameters, I<-1000 =>none */
    in->lastSouParm[i*4+0] = -1000.0;
    in->lastSouParm[i*4+1] = in->lastSouParm[i*4+2] = in->lastSouParm[i*4+3] = 0.0;
  } /* end filling source fitting flags */

  /* Get numbers of channels/IFs */
  BChan   = 1;
  EChan   = in->inDesc->inaxes[in->inDesc->jlocf];
  iChan   = 1 + in->ChInc/2;
  nchleft = EChan;
  BIF     = 1;
  EIF     = in->inDesc->inaxes[in->inDesc->jlocif];
  iIF     = 1;
  done    = FALSE;

  /* Data tables */
  size = in->inDesc->nvis;
  in->inData = g_malloc0(10*size*sizeof(ofloat));
  in->inWt   = g_malloc0(4*size*sizeof(ofloat));
  in->antNo  = g_malloc0(2*size*sizeof(olong));
  in->souNo  = g_malloc0(size*sizeof(olong));

  /* Initialize threads*/
  MakePolnFitFuncArgs (in, err); 

  /* Loop over data in blocks of Channels */
  while (!done) {

    /* Read data, init solutions */
    first = (iIF!=in->IFno) || (iChan==BChan);   /* First of an IF */
    ReadData (in, inUV, &iChan, EChan, &iIF, EIF, first, err);
    if (err->error) Obit_traceback_msg (err, routine, inUV->name);

    if (in->prtLv>=2) {
      Obit_log_error(err, OBIT_InfoErr, "Process IF %d Channel %d no chan %d",
		     in->IFno, in->Chan+in->BChan-1, in->ChWid);
      ObitErrLog(err); 
   }

    /* Set order and count number of parameters being fitted */
    in->nparam = 0;
    /* get model parameters - first antenna */
    for (iant=0; iant<in->nant; iant++) {
      if (in->gotAnt[iant]) {
	/* Loop over parameters */
	for (k=0; k<4; k++) {
	  /* Fitting? */
	  if (in->antFit[iant][k]) {
	    in->antPNumb[iant][k] = in->nparam++;
	  } else in->antPNumb[iant][k] = -999;
	} /* end loop over parameters */
      } /* end if gotAnt */
    } /* end loop over antennas */

    /* Antenna gains if needed */
    if (!in->isCircFeed && in->doFitGain) {
      for (iant=0; iant<in->nant; iant++) {
	if (in->gotAnt[iant]) {
	  /* Loop over parameters */
	  for (k=0; k<2; k++) {
	    /* Fitting? */
	    if (in->antGainFit[iant][k]) {
	      in->antGainPNumb[iant][k] = in->nparam++;
	    } else in->antGainPNumb[iant][k] = -999;
	  } /* end loop over parameters */
	} /* end if gotAnt */
      } /* end loop over antennas */
    } /* End if antenna gains */
   
    /* now source */
    for (isou=0; isou<in->nsou; isou++) {
      /* Loop over parameters */
      for (k=0; k<4; k++) {
	/* Fitting? */
	if (in->souFit[isou][k]) {
	  in->souPNumb[isou][k] = in->nparam++;
	} else in->souPNumb[isou][k] = -999;
      } /* end loop over parameters */
    } /* end loop over sources */

    /* Phase difference parameter number if being fitted */
    if (in->doFitRL) in->PDPNumb = in->nparam++;
    
    /* Init solver stuff first pass */
    nvis       = inUV->myDesc->nvis;
    in->ndata  = nvis*4*2;   /* 4 complex visibililties */
#ifdef HAVE_GSL
    /* Need to rebuild? */
    if ((in->solver!=NULL) || (oldNparam!=in->nparam) || (oldNdata!=in->ndata))  {
      if( in->solver)    gsl_multifit_fdfsolver_free (in->solver); in->solver = NULL;
      if (in->funcStruc) g_free(in->funcStruc);                    in->funcStruc = NULL;
      if (in->work)      gsl_vector_free(in->work);                in->work = NULL;
      if (in->covar)     gsl_matrix_free(in->covar);               in->covar = NULL;
    }
    if (in->solver==NULL) {
      oldNparam = in->nparam;
      oldNdata  = in->ndata;
      T = gsl_multifit_fdfsolver_lmder;
      in->solver    = gsl_multifit_fdfsolver_alloc(T, in->ndata, in->nparam);
      in->funcStruc = g_malloc0(sizeof(gsl_multifit_function_fdf));
      in->work      = gsl_vector_alloc(in->nparam);
      in->covar     = gsl_matrix_alloc(in->nparam, in->nparam);
      /* Function by feed type  */
      if (in->isCircFeed) {
	/* Circular feeds */
	in->funcStruc->f   = &PolnFitFuncOERL;
	in->funcStruc->df  = &PolnFitJacOERL;
	in->funcStruc->fdf = &PolnFitFuncJacOERL;
      } else {
	/* Linear feeds  */
	in->funcStruc->f   = &PolnFitFuncOEXY;
	in->funcStruc->df  = &PolnFitJacOEXY;
	in->funcStruc->fdf = &PolnFitFuncJacOEXY;
      }
    }
  
#endif /* HAVE_GSL */ 

    /* Do actual fitting  */
    /* Do fit - fast method to get better soln, then possibly GSL */
    isOK = doFitFast (in, err);
    if (isOK) {
      if (!strncmp(in->solnType, "LM  ", 4)) doFitGSL (in, err);
      
      endChi2 = GetChi2 (in->nThread, in, polnParmUnspec, 0, 
			 &ParRMS, &XRMS, NULL, NULL, err);  /* Final value */
      if (err->prtLv>=2) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "Final Chi2=%g Par RMS %g X RMS %g ", endChi2, ParRMS, XRMS);
	if (err->error) Obit_traceback_msg (err, routine, inUV->name);
      }
      
    } /* End fit OK */
    
    /* Done? */
    done = iIF > EIF;
    /* Diagnostics */
    if (isOK && (err->prtLv>=4)) {
      /* BL 11-18 sou 1 */
      refAnt = in->refAnt-1;
      fprintf (stderr, "Source 0 IF %d Chan %d\n", in->IFno, in->Chan);
      /* Header by feed type  */
      if (in->isCircFeed) {
	/* Circular feeds */
	fprintf (stderr, " bl     PA          RRr (o,c)         RRi                LLr                LLi                RLr                RLi                LRr                LRi \n");
      } else  {
	/* Linear feeds  */
	fprintf (stderr, " bl     PA          XXr (o,c)         XXi                YYr                YYi                XYr                XYi                YXr                YXi \n");
      }
      for (i=0; i<in->nvis; i++) {
	isou   = MAX(0, in->souNo[i]);    /* Source number */
	ia1    = in->antNo[i*2+0];
	ia2    = in->antNo[i*2+1]; 
	/* DEBUG if ((isou==0) && (ia1==10) && (ia2==17)) {  11-18 */
	/* DEBUG if ((isou==0)&& (ia1==2) && (ia2==3)) {  first source 3-4 */
	if (isou==0) {   /* first source */
	  /* Function by feed type  */
	  if (in->isCircFeed ) {
	    /* Circular feeds */
	    calcmodelRL (in, RArray, i);
	  } else {
	    /* Linear feeds  */
	    calcmodelXY (in, RArray, i);
	  }
	  RRr = RArray[0]; RRi = RArray[1]; 
	  LLr = RArray[6]; LLi = RArray[7]; 
	  RLr = RArray[2]; RLi = RArray[3]; 
	  LRr = RArray[4]; LRi = RArray[5]; 
	  fprintf (stderr, "%2d-%2d %8.5f %8.5f %8.5f  %8.5f %8.5f  %8.5f %8.5f  %8.5f %8.5f  %8.5f %8.5f  %8.5f %8.5f  %8.5f %8.5f  %8.5f %8.5f \n",
		   ia1+1,ia2+1,in->inData[i*10+0],
		   in->inData[i*10+2],  RRr, in->inData[i*10+3], RRi, in->inData[i*10+4], LLr, in->inData[i*10+5], LLi,
		   in->inData[i*10+6],  RLr, in->inData[i*10+7], RLi, in->inData[i*10+8], LRr, in->inData[i*10+9], LRi);
	}
      }
    } /* End debug diagnostics */
    
    /* Write output to Tables  if some valid XPol data */
    isOK = isOK && (XRMS>0.0);
    WriteOutput(in, outUV, isOK, err);
    if (err->error) Obit_traceback_msg (err, routine, inUV->name);

    } /* end loop over blocks of channels */
  
} /* end ObitPolnCalFitFit */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitPolnCalFitClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitPolnCalFitClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitPolnCalFitClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitPolnCalFitClassInfoDefFn (gpointer inClass)
{
  ObitPolnCalFitClassInfo *theClass = (ObitPolnCalFitClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitPolnCalFitClassInit;
  theClass->newObit       = (newObitFP)newObitPolnCalFit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitPolnCalFitClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitPolnCalFitGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitPolnCalFitCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitPolnCalFitClear;
  theClass->ObitInit      = (ObitInitFP)ObitPolnCalFitInit;
  theClass->ObitPolnCalFitCreate = (ObitPolnCalFitCreateFP)ObitPolnCalFitCreate;
  theClass->ObitPolnCalFitFit    = (ObitPolnCalFitFitFP)ObitPolnCalFitFit;
} /* end ObitPolnCalFitClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitPolnCalFitInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
ObitPolnCalFit *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread     = newObitThread();
  in->info       = newObitInfoList(); 
  in->inData     = NULL;
  in->inWt       = NULL;
  in->nant       = 0;
  in->antNo      = NULL;
  in->antParm    = NULL;
  in->antErr     = NULL;
  in->gotAnt     = NULL;
  in->antFit     = NULL;
  in->antGain    = NULL;
  in->antGainErr = NULL;
  in->antGain    = NULL;
  in->antGainFit = NULL;
  in->AntLists   = NULL;
  in->nsou       = 0;
  in->souNo      = NULL;
  in->souParm    = NULL;
  in->souErr     = NULL;
  in->lastSouParm= NULL;
  in->souFit     = NULL;
  in->antPNumb   = NULL;
  in->SouList    = NULL;
  in->souIDs     = NULL;
  in->isouIDs    = NULL;
  in->inDesc     = NULL;
  in->outDesc    = NULL;
  in->RLPhaseIn  = NULL;
  in->RLPhase    = NULL;
  in->RM         = NULL;
  in->PPol       = NULL;
  in->doFitI     = NULL;
  in->doFitV     = NULL;
  in->doFitPol   = NULL;
  in->IFlux0     = NULL;
  in->IFlux1     = NULL;
  in->IFlux2     = NULL;
  in->CPTable    = NULL;
  in->PDTable    = NULL;
  in->BPTable    = NULL;
  in->doFitRL    = FALSE;
  in->doBlank    = TRUE;
  in->doFitGain  = TRUE;
  in->BIF        = 1;
  in->BChan      = 1;
  in->ChWid      = 1;
  in->ChInc      = 1;
  in->CPSoln     = 0;
  in->PDSoln     = 0;
  in->BPSoln     = 0;
  in->doBand     = 0;
  in->BPVer      = 0;
  in->doError    = TRUE;
  in->isCircFeed = TRUE;
  in->maxAnt   = 100;
} /* end ObitPolnCalFitInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitPolnCalFit* cast to an Obit*.
 */
void ObitPolnCalFitClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPolnCalFit *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->antFit) {
    for (i=0; i<in->nant; i++) if (in->antFit[i]) g_free(in->antFit[i]);
    g_free(in->antFit);
  }
  if (in->antGainFit) {
    for (i=0; i<in->nant; i++) if (in->antGainFit[i]) g_free(in->antGainFit[i]);
    g_free(in->antGainFit);
  }
  if (in->souFit) {
    for (i=0; i<in->nsou; i++) if (in->souFit[i]) g_free(in->souFit[i]);
    g_free(in->souFit);
  }
  if (in->antPNumb) {
    for (i=0; i<in->nant; i++) if (in->antPNumb[i]) g_free(in->antPNumb[i]);
    g_free(in->antPNumb);
  }
  if (in->antGainPNumb) {
    for (i=0; i<in->nant; i++) if (in->antGainPNumb[i]) g_free(in->antGainPNumb[i]);
    g_free(in->antGainPNumb);
  }
  if (in->souPNumb) {
    for (i=0; i<in->nsou; i++) if (in->souPNumb[i]) g_free(in->souPNumb[i]);
    g_free(in->souPNumb);
  }
  if (in->AntLists) {
    for (i=0; i<in->numSubA; i++) if (in->AntLists[i]) ObitAntennaListUnref(in->AntLists[i]);
    g_free(in->AntLists);
  }
  if (in->inWt)     g_free(in->inWt);
  if (in->inData)   g_free(in->inData);
  if (in->antNo)    g_free(in->antNo);
  if (in->antParm)  g_free(in->antParm);
  if (in->gotAnt)   g_free(in->gotAnt);
  if (in->antErr)   g_free(in->antErr);
  if (in->antGain)  g_free(in->antGain);
  if (in->antGainErr) g_free(in->antGainErr);
  if (in->souNo)    g_free(in->souNo);
  if (in->souParm)  g_free(in->souParm);
  if (in->souErr)   g_free(in->souErr);
  if (in->lastSouParm) g_free(in->lastSouParm);
  if (in->souIDs)   g_free(in->souIDs);
  if (in->isouIDs)  g_free(in->isouIDs);
  if (in->RLPhaseIn)g_free(in->RLPhaseIn);
  if (in->RLPhase)  g_free(in->RLPhase);
  if (in->RM)       g_free(in->RM);
  if (in->PPol)     g_free(in->PPol);
  if (in->doFitI)   g_free(in->doFitI);
  if (in->doFitV)   g_free(in->doFitV);
  if (in->doFitPol) g_free(in->doFitPol);
  if (in->IFlux0)   g_free(in->IFlux0);
  if (in->IFlux1)   g_free(in->IFlux1);
  if (in->IFlux2)   g_free(in->IFlux2);
  if (in->thArgs) KillPolnFitFuncArgs (in);
  in->inDesc    = ObitUVDescUnref(in->inDesc);
  in->outDesc   = ObitUVDescUnref(in->outDesc);
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  in->SouList   = ObitSourceListUnref(in->SouList);
  in->oneSource = ObitSourceUnref(in->oneSource);
  in->CPTable   = ObitTableUnref(in->CPTable);
  in->PDTable   = ObitTableUnref(in->PDTable);
  in->BPTable   = ObitTableUnref(in->BPTable);
#ifdef HAVE_GSL
    if( in->solver)    gsl_multifit_fdfsolver_free (in->solver);
    if (in->funcStruc) g_free(in->funcStruc);
    if (in->work)      gsl_vector_free(in->work);
    if (in->covar)     gsl_matrix_free(in->covar);
#endif /* HAVE_GSL */ 

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitPolnCalFitClear */

/**
 * Write contents on in to tables on outUV
 * Tables initialized on call with first thread for chan=IF=1
 * (including BChan, BIF)
 * \param in       Fitting object
 * \param outUV    UV date with results in tables
 * \param isOK     Something worth writing
 * \param err      Obit error stack object.
 */
static void WriteOutput (ObitPolnCalFit* in, ObitUV *outUV, 
			 gboolean isOK, ObitErr *err)
{
  oint numPol, numIF, numChan;
  ObitIOCode retCode;
  ObitTableBP *oldBPTab=NULL;
  gboolean sameBP=FALSE, crazy;
  ObitIOAccess access;
  ObitUVDesc *IODesc;
  gchar *routine = "ObitPolnCalFit:WriteOutput";

  /* error checks */
  if (err->error) return;


  /* Need to initialize? First channel and IF */
  if (((in->Chan+in->BChan-1-in->ChInc/2)==1) && 
      ((in->IFno+in->BIF-1)==1)) {
    /* Open output */
    retCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) 
      Obit_traceback_msg (err, routine, outUV->name);

    /* How many of things? Use underlying sizes */
    IODesc = (ObitUVDesc*)outUV->myIO->myDesc;
    numPol  = MIN (2,IODesc->inaxes[IODesc->jlocs]);
    numIF   = IODesc->inaxes[IODesc->jlocif];
    numChan = IODesc->inaxes[IODesc->jlocf];

   /* Source table */
    in->CPTable = newObitTableCPValue ("Source", (ObitData*)outUV, &in->CPSoln, 
				       OBIT_IO_WriteOnly, numIF, numChan, err);
    InitSourceTab(in, err);

    /* Instrumental poln table */
    in->PDTable = newObitTablePDValue ("Instrum", (ObitData*)outUV, &in->PDSoln, 
				       OBIT_IO_WriteOnly, numPol, numIF, numChan, 
				       err);
    InitInstrumentalTab(in, err);

    /* BP table used if fitting R/L X/Y phase difference or X/Y gain*/
    if (in->doFitRL || in->doFitGain) {
      /* Bandpass table - copy and update if doBand>0 */
      if (in->doBand>0) {
	oldBPTab = newObitTableBPValue ("Temp BP", (ObitData*)outUV, &in->BPVer, 
					OBIT_IO_ReadOnly, numPol, numIF, numChan, 
					err);
	/* Check that you're not about to clobber the input */
	sameBP = in->BPVer==in->BPSoln;
	if (sameBP) access = OBIT_IO_ReadWrite;
	access = OBIT_IO_WriteOnly;
	/* Warn if overwriting */
	if (sameBP) {
	  Obit_log_error(err, OBIT_InfoWarn, "Modifying input BP Table");
	  ObitErrLog(err); 
	}
      } else { /* No prior bandpass cal */
	if (in->BPSoln<=0) access = OBIT_IO_WriteOnly;
	else access = OBIT_IO_ReadWrite;
      }
      
      /* Create or open output */
      in->BPTable = newObitTableBPValue ("Bandpass", (ObitData*)outUV, &in->BPSoln, 
					 access , numPol, numIF, numChan, 
					 err);
      /* Copy old to new */
      if ((in->doBand>0) && !sameBP) {
	in->BPTable = ObitTableBPCopy (oldBPTab, in->BPTable, err);
	oldBPTab = ObitTableBPUnref(oldBPTab);
      } else {
	InitBandpassTab(in, err);
      }
    } /* end setup BP table */

    /* Close output */
    retCode = ObitUVClose (outUV, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) 
      Obit_traceback_msg (err, routine, outUV->name);
    /* end init */
  } else { /* Make sure output tables defined */
    /* How many of things? */
    numPol  = MIN (2,in->outDesc->inaxes[in->outDesc->jlocs]);
    numIF   = in->outDesc->inaxes[in->outDesc->jlocif];
    numChan = in->outDesc->inaxes[in->outDesc->jlocf];

    if (!in->CPTable) 
      in->CPTable = newObitTableCPValue ("Source", (ObitData*)outUV, &in->CPSoln, 
					 OBIT_IO_ReadWrite, numIF, numChan, err);
    if (!in->PDTable) 
      in->PDTable = newObitTablePDValue ("Instrum", (ObitData*)outUV, &in->PDSoln, 
					 OBIT_IO_ReadWrite, numPol, numIF, numChan, 
					 err);
    if ((!in->BPTable) && (in->doFitRL))
      in->BPTable = newObitTableBPValue ("Bandpass", (ObitData*)outUV, &in->BPSoln, 
					 OBIT_IO_ReadWrite, numPol, numIF, numChan, 
					 err);
  }

   /* Check for crazy antenna solutions and reset defaults */
  crazy = CheckCrazy(in, err);
  if (crazy) isOK = FALSE;

  /* Update Tables */
  if (isOK) UpdateSourceTab(in, err);
  UpdateInstrumentalTab(in, isOK, err);
  if (in->doFitRL || in->doFitGain) UpdateBandpassTab(in, isOK, err);
  if (err->error)  Obit_traceback_msg (err, routine, outUV->name);

} /* end WriteOutput */

/**
 * Read data and load arrays.
 * Averages ChWid channels centered on each channel.
 * Require all 4 correlations to be valid in each channel
 * \param in    Fitting object
 * \param inUV  Averaged/calibrated/divided UV data to be fitted
 * \param iChan Current first channel (1-rel) in thread group, 
 *              updated to value for next call on exit;
 * \param EChan Highest channel
 * \param iIF   Current IF number (1-rel) , updated to value for next call
 * \param EIF   Highest IF
 * \param init  If True, solutions also initialized
 * \param err   Obit error stack object.
 * \return number of threads with data loaded.
 */
static void ReadData (ObitPolnCalFit *in, ObitUV *inUV, olong *iChan, 
		      olong EChan, olong *iIF, olong EIF, 
		      gboolean init,  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  ObitSource *curSource=NULL;
  olong ivis, ant1, ant2, isou, jsou, lastSou=-999, jvis, istok, suba;
  olong nChan, bch, ech, cnt, i, j, indx, indx1, indx2, indx3, indx4;
  olong ChInc=in->ChInc, jChan, jIF;
  ofloat *buffer, cbase, curPA1=0.0, curPA2=0.0, curTime=0.0, lastTime=-1.0e20;
  ofloat sumRe, sumIm, sumWt; 
  odouble lambdaRef, lambda, ll, sum;
  gboolean OK, allBad;
  gchar *routine = "ObitPolnCalFit:ReadData";

  /* Previous error? */
  if (err->error) return;

  lambdaRef = VELIGHT/in->inDesc->freq;  /* Wavelength at reference freq */

  /* initial center channel/IF */
  jChan = *iChan;
  jIF   = *iIF;
  nChan = in->inDesc->inaxes[in->inDesc->jlocf];

  /* Set channel/IF in thread */
  ChInc = MAX (1,ChInc);
  in->Chan  = (*iChan);
  in->IFno = *iIF;
  (*iChan) += ChInc;   /* Next channel */
  if (((*iChan)>EChan) && ((*iIF)<EIF)) {  /* Done with IF? */
    (*iChan) = 1 + ChInc/2;
    (*iIF)++;
  }
  
  /* Keep track of antennas with data */
  for (i=0; i<in->nant; i++) in->gotAnt[i] = FALSE;

  /* Frequency */
  indx = (jIF-1)*nChan + (jChan-1);
  in->freq = in->inDesc->freqArr[indx];
  /* Set RL phases of calibrators */
  lambda = VELIGHT/in->freq;
  /* Loop over calibrators */
  for (i=0; i<in->nsou; i++) {
    /* Set R-L phase if given */
    if (in->RLPhaseIn[i]>-900.0) {
      in->RLPhase[i] = 	in->RLPhaseIn[i] + 
	2.0 * (lambda*lambda - lambdaRef*lambdaRef) * in->RM[i];
    } else in->RLPhase[i] = 0.0;
  }

  /* Beginnning and end channel */
  bch = in->Chan - in->ChWid/2;
  ech = in->Chan + in->ChWid/2;
  /* Force to range */
  bch = MAX (1, bch);
  ech = MIN (nChan, ech);
  /* Convert to 0-rel */
  bch--;
  ech--;

  /* Open input */
  retCode = ObitUVOpen (inUV, OBIT_IO_ReadOnly, err);
  if (err->error) goto cleanup;

  /* loop loading data */
  jvis = 0;
  while (retCode == OBIT_IO_OK) {
    
    /* read buffer */
    retCode = ObitUVRead (inUV, NULL, err);
    if (retCode == OBIT_IO_EOF) break; /* done? */
    if (err->error) goto cleanup;
    
    /* Loop over buffer */
    buffer = inUV->buffer;
    for (ivis=0; ivis<inUV->myDesc->numVisBuff; ivis++) {
      
      /* Check array bounds */
      if (jvis>=in->inDesc->nvis) {
	Obit_log_error(err, OBIT_Error, "%s: Exceeded vis buffer size", routine);
	goto cleanup;
      }

      /* Check that all Stokes correlations in each channel are present,
	 if only partial, flag the rest */
      allBad = TRUE;
      for (i=bch; i<=ech; i++) {
	indx = inUV->myDesc->nrparm + i*inUV->myDesc->incf
	  + (jIF-1)*inUV->myDesc->incif + 2;
	indx1 = indx + 0*inUV->myDesc->incs;
	indx2 = indx + 1*inUV->myDesc->incs;
	indx3 = indx + 2*inUV->myDesc->incs;
	indx4 = indx + 3*inUV->myDesc->incs;
	OK = (buffer[indx1]>0.0) && (buffer[indx2]>0.0) && 
	  (buffer[indx3]>0.0) && (buffer[indx4]>0.0);
	if (!OK) {  /* Kill 'em all */
	  buffer[indx1] = buffer[indx2] = buffer[indx3] = buffer[indx4] = 0.0;
	} else allBad = FALSE;  /* Some OK */
      } /* end checking data loop */
      
      /* if all data flagged skip to end */
      if (allBad) goto endloop;

      /* Get info - source ID */
      if (inUV->myDesc->ilocsu>=0) isou = buffer[inUV->myDesc->ilocsu]+0.5;
      else isou = in->souIDs[0];
      jsou = in->isouIDs[isou-1];  /* Calibrator number */
      /* Antennas */
      cbase = buffer[inUV->myDesc->ilocb]; 
      ant1 = (cbase / 256.0) + 0.001;
      ant2 = (cbase - ant1 * 256) + 0.001;
      suba = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 1.5);
      
      /* Source change? */
      if (isou!=lastSou) {
	lastSou = isou;
	if (in->SouList) curSource = in->SouList->SUlist[isou-1];
	else curSource = in->oneSource;
      }
      
      /* Get parallactic angle when time changes */
      curTime = buffer[inUV->myDesc->iloct];
      if (curTime>lastTime) {
	curPA1 = ObitAntennaListParAng (in->AntLists[suba-1], ant1, curTime, curSource);
	curPA2 = ObitAntennaListParAng (in->AntLists[suba-1], ant2, curTime, curSource);
	lastTime = curTime;
      }
      
      /* Save info */
      in->souNo[jvis]     = jsou;  /* relative to list, not data ID */
      in->antNo[jvis*2+0] = ant1-1;
      in->antNo[jvis*2+1] = ant2-1;
      /* Parallactic angle */
      in->inData[jvis*10]    = curPA1;
      in->inData[jvis*10+1]  = curPA2;
      /* getAnt flags */
      in->gotAnt[ant1-1] = TRUE;
      in->gotAnt[ant2-1] = TRUE;
      
      /* Loop over Stokes */
      for (istok=0; istok<4; istok++) {
	
	/* Average vis */
	sumRe = sumIm = sumWt = 0; cnt = 0;
	for (i=bch; i<=ech; i++) {
	  /* Where in buffer? */
	  indx = inUV->myDesc->nrparm + i*inUV->myDesc->incf
	    + (jIF-1)*inUV->myDesc->incif + istok*inUV->myDesc->incs;
	  if (buffer[indx+2]>0) {
	    cnt++;
	    sumRe += buffer[indx];
	    sumIm += buffer[indx+1];
	    sumWt += buffer[indx+2];
	  }
	} /* end average loop */
	  /* Average */
	if ((cnt>0) && (sumWt>0.0)) {
	  sumRe /= cnt;
	  sumIm /= cnt;
	} else {  /* invalid */
	  sumRe = 0.0;
	  sumIm = 0.0;
	  sumWt = 0.0;
	}
	/* Weight */
	in->inWt[jvis*4+istok] = sumWt;
	  /* Data */
	in->inData[jvis*10+(istok+1)*2]   = sumRe;
	in->inData[jvis*10+(istok+1)*2+1] = sumIm;
      } /* end Stokes loop */
      
      jvis++;  /* Data array visibility index */
    endloop:  /* to here if data bad */
      buffer += inUV->myDesc->lrec;
    } /* end loop over vis */
  } /* end loop reading data */

  in->nvis = jvis;  /* How many good data? */
  /* Initialize solutions if init */
  if (init) {
    /* first zero everything 
       for (i=0; i<in->nsou; i++) {
       for (j=0; j<4; j++) in->souParm[i*4+j] = 0.0;
       }*/
    for (i=0; i<in->nant; i++) {
      for (j=0; j<4; j++) in->antParm[i*4+j] = 0.0;
      if (in->isCircFeed) {
	/* Circular feeds ori_r, elip_r, ori_l, elip_l */
	in->antParm[i*4+0] = 0.0;
	in->antParm[i*4+1] = G_PI/4.0;
	in->antParm[i*4+2] =  0.0;
	in->antParm[i*4+3] =  -G_PI/4.0;
      } else  {
 	/* Linear feeds, if the antenna angles on sky are given in the AntennaList[0] use then, 
           otherwise assume Feeds are X, Y
	   ori_x, elip_x, ori_y, elip_y,  */
	in->antParm[i*4+1] = 0.0;  /* ellipticity */
	in->antParm[i*4+3] = 0.0;
	if (fabs (in->AntLists[0]->ANlist[i]->FeedAPA-in->AntLists[0]->ANlist[i]->FeedBPA)<1.0) {
	  /* Assume X, Y */
	  in->antParm[i*4+0] = 0.0;
	  in->antParm[i*4+2] = +G_PI/2.0;
	  /* DEBUG KAT 
	  in->antParm[i*4+2] = 0.0;
	  in->antParm[i*4+0] = +G_PI/2.0;*/
	} else {  /* Use what's in the AN table as initial convert to radians */
	  in->antParm[i*4+0] = in->AntLists[0]->ANlist[i]->FeedAPA*DG2RAD;
	  in->antParm[i*4+2] = in->AntLists[0]->ANlist[i]->FeedBPA*DG2RAD;
	} /* end if antenna angle given */
      }
    }
  } /* end init */

  /* Some parameters always need updating */
  in->PD = 0.0;  /* R-L (X/Y) phase difference */
  for (i=0; i<in->nsou; i++) {
    
    /* Stokes I - evaluate spectrum in ln(nu/nu_ref) */
    ll = log(in->freq/in->inDesc->freq);
    sum = in->IFlux1[i]*ll + in->IFlux2[i]*ll*ll;
    in->souParm[i*4] = exp(sum) * in->IFlux0[i];
    /* Fixed linear poln? */
    if (in->PPol[i]>0.0) {
      in->souParm[i*4+1] = in->PPol[i] * in->souParm[i*4] * cos(in->RLPhase[i]);
      in->souParm[i*4+2] = in->PPol[i] * in->souParm[i*4] * sin(in->RLPhase[i]);
    }
    /* Add 1% linear poln if fitting on init 
       if (init && in->souFit[i][1]) {
       in->souParm[i*4+1] = 0.01*in->souParm[i*4];
       }*/
   
  } /* end loop over source */
  
  /* Cleanup */
 cleanup:
  /* Close data */
  ObitUVClose (inUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Set iChan, *iIF for next */
  if ((*iChan)>nChan) {  /* Done with IF? */
    (*iChan) = 1;
    (*iIF)++;
  }

  return; 
} /* end ReadData */

/**
 * Initialize new Source poln.(CP) table
 * All zero entries
 * \param in    Fitting object
 * \param err   Obit error stack object.
 */
static void InitSourceTab(ObitPolnCalFit* in, ObitErr *err) 
{
  olong i, irow, nif, nchan, isou;
  ObitTableCPRow *row=NULL;
  gchar *routine = "ObitPolnCalFit:InitSourceTab";
  
  /* Clear any existing rows */
  ObitTableClearRows ((ObitTable*)in->CPTable, err);

  /* Open  */
  ObitTableCPOpen (in->CPTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* define row */
  row = newObitTableCPRow(in->CPTable);
  ObitTableCPSetRow (in->CPTable, row, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  /* Set header values */
  in->CPTable->FreqID = 0;

  nif   = in->outDesc->inaxes[in->outDesc->jlocif];
  nchan = in->outDesc->inaxes[in->outDesc->jlocf];

 /* Initialize Row */
  row->SourID = 0;
  for (i=0; i<nif*nchan; i++) row->IPol[i] = 0.0;
  for (i=0; i<nif*nchan; i++) row->QPol[i] = 0.0;
  for (i=0; i<nif*nchan; i++) row->UPol[i] = 0.0;
  for (i=0; i<nif*nchan; i++) row->VPol[i] = 0.0;

  /* Loop over sources writing */
  for (isou=0; isou<in->nsou; isou++) {
    row->SourID = in->souIDs[isou];
    irow = -1;
    ObitTableCPWriteRow (in->CPTable, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* End loop onver source writing */

  /* Close */
  ObitTableCPClose (in->CPTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  row = ObitTableCPRowUnref(row);  /* Cleanup */

} /* end InitSourceTab */

/**
 * Initialize new Instrumental poln.(PD) table
 * All zero entries
 * \param in    Fitting object
 * \param err   Obit error stack object.
 */
static void InitInstrumentalTab(ObitPolnCalFit* in, ObitErr *err) 
{
  olong i, irow, iant, npol, nif, nchan, nant;
  ObitTablePDRow *row=NULL;
  gchar *routine = "ObitPolnCalFit:InitInstrumentalTab";
  
  /* Clear any existing rows */
  ObitTableClearRows ((ObitTable*)in->PDTable, err);

  /* Open Write only */
  ObitTablePDOpen (in->PDTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* define row */
  row = newObitTablePDRow(in->PDTable);
  ObitTablePDSetRow (in->PDTable, row, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  /* Set header values */
  in->PDTable->numAnt = in->nant;
  /* Mark soln type */
  strncpy (in->PDTable->polType, "ORI-ELP", MAXKEYCHARTABLEPD);

  nant  = in->nant;
  nif   = in->outDesc->inaxes[in->outDesc->jlocif];
  nchan = in->outDesc->inaxes[in->outDesc->jlocf];
  npol  = MIN (2,in->outDesc->inaxes[in->outDesc->jlocs]);

 /* Initialize Row */
  /*  row->SourID  = 0;*/
  row->antNo   = 0;
  row->SubA    = 0;
  row->FreqID  = 0;
  row->RefAnt  = in->refAnt;
  for (i=0; i<nif*nchan; i++) row->RLPhase[i] = 0.0;
  /* Circular feed default */
  if (in->isCircFeed) {
    for (i=0; i<nif*nchan; i++) row->Real1[i] = G_PI/4.0;
    for (i=0; i<nif*nchan; i++) row->Imag1[i] = 0.0;
    if (npol>1) {
      for (i=0; i<nif*nchan; i++) row->Real2[i] = -G_PI/4.0;
      for (i=0; i<nif*nchan; i++) row->Imag2[i] = 0.0;
    }
  } else {
    /* Linear feeds  ori_x, elip_x, ori_y, elip_y,  */
    for (i=0; i<nif*nchan; i++) row->Real1[i] = G_PI/4.0;
    for (i=0; i<nif*nchan; i++) row->Imag1[i] = G_PI/2.0;
    if (npol>1) {
      for (i=0; i<nif*nchan; i++) row->Real2[i] = 0.0;
      for (i=0; i<nif*nchan; i++) row->Imag2[i] = -G_PI/2.0;
    }
  } 

  /* Loop over antennas writing */
  for (iant=0; iant<in->nant; iant++) {
    row->antNo = iant+1;
    irow = -1;
    ObitTablePDWriteRow (in->PDTable, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* End loop over source writing */

  /* Close */
  ObitTablePDClose (in->PDTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  row = ObitTablePDRowUnref(row);  /* Cleanup */

} /* end InitInstrumentalTab */

/**
 * Initialize new Bandpass.(BP) table 
 * All zero entries
 * \param in    Fitting object
 * \param err   Obit error stack object.
 */
static void InitBandpassTab(ObitPolnCalFit* in, ObitErr *err) 
{
  olong i, irow, iant, npol, nif, nchan, nant;
  ObitUVDesc *desc;
  ObitTableBPRow *row=NULL;
  gchar *routine = "ObitPolnCalFit:InitBandpassTab";
  
  /* Initialize */
  nant  = in->nant;
  nif   = in->outDesc->inaxes[in->outDesc->jlocif];
  nchan = in->outDesc->inaxes[in->outDesc->jlocf];
  npol  = MIN (2,in->outDesc->inaxes[in->outDesc->jlocs]);

  /* Clear any existing rows */
  ObitTableClearRows ((ObitTable*)in->BPTable, err);

  /* Open  */
  ObitTableBPOpen (in->BPTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* define row */
  row = newObitTableBPRow(in->BPTable);
  ObitTableBPSetRow (in->BPTable, row, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  /* Set header values */
  in->BPTable->numAnt = nant;

  /* Initialize Row */
  desc = (ObitUVDesc*)in->outDesc;
  row->BW           = desc->cdelt[desc->jlocf];
  row->ChanShift[0] = 0.0;
  row->ChanShift[1] = 0.0;
  row->RefAnt1      = in->refAnt;
  row->RefAnt2      = in->refAnt;
  row->SubA         = 0;
  row->FreqID       = 0;
  row->TimeI        = 24.0;
  row->SourID       = 0;
  for (i=0; i<nif; i++) row->ChanShift[i] = 0.0;
  for (i=0; i<nif; i++) row->Weight1[i]   = 1.0;
  if (npol>1) for (i=0; i<nif; i++) row->Weight2[i] = 1.0;
  for (i=0; i<nchan*nif; i++) { 
    row->Real1[i]   = 1.0;
    row->Imag1[i]   = 0.0;
    if (npol>1) {
      row->Real2[i]   = 1.0;
      row->Imag2[i]   = 0.0;
    }
  }
  /* Loop over antennas writing */
  for (iant=0; iant<in->nant; iant++) {
    row->antNo = iant+1;
    irow = -1;
    ObitTableBPWriteRow (in->BPTable, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* End loop onver source writing */

  /* Close */
  ObitTableBPClose (in->BPTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  row = ObitTableBPRowUnref(row);  /* Cleanup */


} /* end InitBandpassTab */

/**
 * Update Source poln.(CP) table
 * Results obtained from in
 * \param in    Fitting object
 * \param err   Obit error stack object.
 */
static void UpdateSourceTab(ObitPolnCalFit* in, ObitErr *err) 
{
  olong irow, nif, nchan, indx, ich, iif, isou, jsou;
  olong chanOff=in->BChan-1, ifOff=in->BIF-1, chans[2];
  ObitTableCPRow *row=NULL;
  gchar *routine = "ObitPolnCalFit:UpdateSourceTab";
  
  /* Open */
  ObitTableCPOpen (in->CPTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* define row */
  row = newObitTableCPRow(in->CPTable);
  ObitTableCPSetRow (in->CPTable, row, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  nif   = in->outDesc->inaxes[in->outDesc->jlocif];
  nchan = in->outDesc->inaxes[in->outDesc->jlocf];

  /* 0-rel Channels covered in ChInc */
  chans[0] = in->Chan-1+chanOff - in->ChInc/2;
  chans[0] = MAX (0, chans[0]);
  chans[1] = in->Chan-1+chanOff + in->ChInc/2;
  chans[1] = MIN (nchan-1, chans[1]);
  
  /* Loop over row updating */
  for (irow=1; irow<=in->CPTable->myDesc->nrow; irow++) {
    /* Read previous */
    ObitTableCPReadRow (in->CPTable, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    isou = row->SourID - 1;
    jsou = in->isouIDs[isou];  /* Calibrator number */

    /* Loop over chans */
    for (ich=chans[0]; ich<=chans[1]; ich++) {
      iif = in->IFno-1+ifOff;    /* 0 rel IF */
      indx = iif*nchan + ich;
      row->IPol[indx] = in->souParm[jsou*4+0];
      row->QPol[indx] = in->souParm[jsou*4+1];
      row->UPol[indx] = in->souParm[jsou*4+2];
      row->VPol[indx] = in->souParm[jsou*4+3];
    } /* end loop over channels */

    /* Rewrite */
    ObitTableCPWriteRow (in->CPTable, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* End loop onver source writing */

  /* Close */
  ObitTableCPClose (in->CPTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  row = ObitTableCPRowUnref(row);  /* Cleanup */

} /* end UpdateSourceTab */

/**
 * Update Instrumental poln.(PD) table
 * if in->doBlank, blank failed solutions, else noop
 * \param in    Fitting object
 * \param isOK  was fitting successful?
 * \param err   Obit error stack object.
 */
static void UpdateInstrumentalTab(ObitPolnCalFit* in, gboolean isOK, 
				  ObitErr *err) 
{
  olong irow, npol, nif, nchan, indx, ich, iif, iant;
  ofloat fblank = ObitMagicF();
  olong chanOff=in->BChan-1, ifOff=in->BIF-1, chans[2];
  ObitTablePDRow *row=NULL;
  gchar *routine = "ObitPolnCalFit:UpdateInstrumentalTab";
  
  /* If doBlank FALSE and this one failed, simple return */
  if (!isOK && (in->doBlank==FALSE)) return;
  
  /* Open */
  ObitTablePDOpen (in->PDTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* define row */
  row = newObitTablePDRow(in->PDTable);
  ObitTablePDSetRow (in->PDTable, row, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  nif   = in->outDesc->inaxes[in->outDesc->jlocif];
  nchan = in->outDesc->inaxes[in->outDesc->jlocf];
  npol  = MIN (2,in->outDesc->inaxes[in->outDesc->jlocs]);

  /* 0-rel Channels covered in ChInc */
  chans[0] = in->Chan-1+chanOff - in->ChInc/2;
  chans[0] = MAX (0, chans[0]);
  chans[1] = in->Chan-1+chanOff + in->ChInc/2;
  chans[1] = MIN (nchan-1, chans[1]);
  
  /* Loop over row updating */
  for (irow=1; irow<=in->PDTable->myDesc->nrow; irow++) {
    /* Read previous */
    ObitTablePDReadRow (in->PDTable, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    iant = row->antNo - 1;

    /* Update */
    /* Loop over chans */
    for (ich=chans[0]; ich<=chans[1]; ich++) {
      iif = in->IFno-1+ifOff;  /* 0 rel IF */
      indx = iif*nchan + ich;
      
         /* OK solution? */
      if (isOK && in->gotAnt[iant]) {
	row->RLPhase[indx] = in->PD*RAD2DG;  /* R-L phase difference */
	row->Real1[indx]   = in->antParm[iant*4+1];
	row->Imag1[indx]   = in->antParm[iant*4+0];
	if (npol>1) {
	  row->Real2[indx] = in->antParm[iant*4+3];
	  row->Imag2[indx] = in->antParm[iant*4+2];
	}
      } else { /* failed solution */
	row->RLPhase[indx] = fblank;
	row->Real1[indx]   = fblank;
	row->Imag1[indx]   = fblank;
	if (npol>1) {
	  row->Real2[indx] = fblank;
	  row->Imag2[indx] = fblank;
	}
      }
    } /* end loop over channels */
    /* Rewrite */
    ObitTablePDWriteRow (in->PDTable, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* End loop onver source writing */

  /* Close */
  ObitTablePDClose (in->PDTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  row = ObitTablePDRowUnref(row);  /* Cleanup */

} /* end UpdateInstrumentalTab */

/**
 * Update Bandpass (BP) table
 * if in->doBlank, blank failed solutions, else noop
 * \param in    Fitting object
 * \param isOK  was fitting successful?
 * \param err   Obit error stack object.
 */
static void UpdateBandpassTab(ObitPolnCalFit* in, gboolean isOK, ObitErr *err) 
{
  olong irow, npol, nif, nchan, indx, ich, iif, iant;
  ofloat amp1, amp2, phase, fblank = ObitMagicF();
  ObitTableBPRow *row=NULL;
  olong chanOff=in->BChan-1, ifOff=in->BIF-1, chans[2];
  gchar *routine = "ObitPolnCalFit:UpdateBandpassTab";
  
  /* If doBlank FALSE and this one failed, simple return */
  if (!isOK && (in->doBlank==FALSE)) return;
  
  /* Open */
  ObitTableBPOpen (in->BPTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* define row */
  row = newObitTableBPRow(in->BPTable);
  ObitTableBPSetRow (in->BPTable, row, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  nif   = in->outDesc->inaxes[in->outDesc->jlocif];
  nchan = in->outDesc->inaxes[in->outDesc->jlocf];
  npol  = MIN (2,in->outDesc->inaxes[in->outDesc->jlocs]);

  /* 0-rel Channels covered in ChInc */
  chans[0] = in->Chan-1+chanOff - in->ChInc/2;
  chans[0] = MAX (0, chans[0]);
  chans[1] = in->Chan-1+chanOff + in->ChInc/2;
  chans[1] = MIN (nchan-1, chans[1]);

  /* Fitted values - as correction */
  if (in->doFitRL) phase = -in->PD;
  else             phase = 0.0;
  
  /* Loop over row updating */
  for (irow=1; irow<=in->BPTable->myDesc->nrow; irow++) {
    /* Read previous */
    ObitTableBPReadRow (in->BPTable, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    iant = row->antNo - 1;
    /* Amplitudes of gains */
    if (in->doFitGain) 
      {amp1 = in->antGain[iant*2];amp2 = in->antGain[iant*2+1];}
    else {amp1 = 1.0; amp2 = 1.0;}

    /* Update */
    /* Loop over chans */
    for (ich=chans[0]; ich<=chans[1]; ich++) {
      iif = in->IFno-1+ifOff;  /* 0 rel IF */
      indx = iif*nchan + ich;
      
      if (isOK && in->gotAnt[iant]) {
	row->Real1[indx] = amp1;
	row->Imag1[indx] = 0.0;
	if (npol>1) {
	  row->Real2[indx] = amp2 * cos(phase);
	  row->Imag2[indx] = amp2 * sin(phase);
	}
      } else { /* failed solution */
	row->Real1[indx] = fblank;
	row->Imag1[indx] = fblank;
	if (npol>1) {
	  row->Real2[indx] = fblank;
	  row->Imag2[indx] = fblank;
	}
      }
      } /* end loop over channels */
    /* Rewrite */
    ObitTableBPWriteRow (in->BPTable, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* End loop onver source writing */

  /* Close */
  ObitTableBPClose (in->BPTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  row = ObitTableBPRowUnref(row);  /* Cleanup */

} /* end UpdateBandpassTab */

/**
 * Fit spectra to calibrator IPol flux densities in Source list
 * \param in       Fitting object
 * \param err      Obit error stack object.
 */
static void FitSpectra (ObitPolnCalFit *in, ObitErr *err)
  {
    ofloat *parms=NULL, *sigma=NULL;
    olong i, isou, jsou, nterm=3, nfreq;
    gchar *routine="ObitPolnCalFit:FitSpectra";

    if (err->error) return;  /* Error exists? */

    /* Dummy error array */
    nfreq = in->inDesc->inaxes[in->inDesc->jlocif];
    sigma = g_malloc(nfreq*sizeof(ofloat));
    for (i=0; i<nfreq; i++) sigma[i] = 0.01;
    nterm = MIN (3, nfreq);

    /* Loop over calibrators */
    for (isou=0; isou<in->nsou; isou++) {
      jsou = MAX (1, in->souIDs[isou]);  /* Source ID in data */
      parms = ObitSpectrumFitSingle (nfreq, nterm, in->inDesc->freq, 
				     in->inDesc->freqIF, 
				     in->SouList->SUlist[jsou-1]->IFlux, sigma, 
				     FALSE, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      in->IFlux0[isou] = parms[0];
      in->IFlux1[isou] = parms[1];
      in->IFlux2[isou] = parms[2];
      g_free(parms);
    }

    if (sigma) g_free(sigma);
 } /* end FitSpectra */

/**
 * Determine Fit Poln parameters
 * Uses GSL nonlinear solver
 * \param in    Fitting object
 * \param err   Obit error stack object.
 */
static void doFitGSL (ObitPolnCalFit *in, ObitErr *err)
{
  /*odouble difParam, difChi2=0.0, endChi2=0.0;*/
  ofloat fpol, fpa;
  olong isou=0, iant, i, k=0, iter;
  olong nparam, ndata, nvalid, nvis, j;
  odouble sumwt, chi2Test=0.0; 
  double epsabs=1.0e-5, epsrel=1.0e-4;  /* Stopping criteria */
  olong maxIter=10;                     /* Stopping criteria */
  int status;
  ofloat chi2, initChi2, ParRMS, XRMS;
  gchar *routine="ObitPolnCalFit:doFitGSL";

#ifdef HAVE_GSL
  gsl_multifit_fdfsolver *solver = in->solver;
  gsl_matrix *covar              = in->covar;
  gsl_vector *work               = in->work;

  if (err->error) return;  /* Error exists? */
  
  /* Set initial parameters */
  /* first antenna */
  for (iant=0; iant<in->nant; iant++) {
    /* Loop over parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if ((in->antFit[iant][k]) && (in->gotAnt[iant])) {
	j = in->antPNumb[iant][k];
	gsl_vector_set(work, j, in->antParm[iant*4+k]);
      }
    } /* end loop over parameters */
  } /* end loop over antennas */

    /* Antenna gains if needed */
  if (!in->isCircFeed && in->doFitGain) {
    for (iant=0; iant<in->nant; iant++) {
      if (in->gotAnt[iant]) {
	/* Loop over parameters */
	for (k=0; k<2; k++) {
	  /* Fitting? */
	  if (in->antGainFit[iant][k]) {
	    j = in->antGainPNumb[iant][k];
	    gsl_vector_set(work, j, in->antGain[iant*2+k]);
	  } 
	} /* end loop over parameters */
      } /* end if gotAnt */
    } /* end loop over antennas */
  } /* End if antenna gains */
  
  /* now source */
  for (isou=0; isou<in->nsou; isou++) {
    /* Loop over parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if (in->souFit[isou][k]) {
	j = in->souPNumb[isou][k];
	gsl_vector_set(work, j, in->souParm[isou*4+k]);
      }
    } /* end loop over parameters */
  } /* end loop over sources */

  /* Phase difference  if being fitted */
  if (in->doFitRL) {
    j = in->PDPNumb;
    gsl_vector_set(work, j, in->PD);
  }

  /* Set up fitting */
  nparam                = in->nparam;
  ndata                 = in->ndata;
  nvis                  = in->nvis;
  in->funcStruc->n      = ndata;
  in->funcStruc->p      = nparam;
  in->funcStruc->params = in;
  iter = 0;

  /* init solver */
  gsl_multifit_fdfsolver_set (solver, in->funcStruc, work);

  /* initial Chi2 */
  if (err->prtLv>=3) {
    initChi2 = GetChi2 (in->nThread, in, polnParmUnspec, 0, 
			&ParRMS, &XRMS, NULL, NULL, err);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Initial LM Chi2=%g Par RMS %g X RMS %g ", 
		   initChi2, ParRMS, XRMS);
    if (err->error) Obit_traceback_msg (err, routine, "diagnostic");
  } /* end initial Chi2 */

  /* iteration loop */
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    /* current  Chi2 */
    if (err->prtLv>=3) {
      chi2 = GetChi2 (in->nThread, in, polnParmUnspec, 0, 
		      &ParRMS, &XRMS, NULL, NULL, err);
    Obit_log_error(err, OBIT_InfoErr, 
		   "current LM Chi2=%g Par RMS %g X RMS %g", 
		   chi2, ParRMS, XRMS);
    if (err->error) Obit_traceback_msg (err, routine, "diagnostic");
  } /* end Chi2 */

    /* Diagnostics */
    if (err->prtLv>=3) {
      sumwt = (ofloat)gsl_blas_dnrm2(solver->f);
      /*Obit_log_error(err, OBIT_InfoErr,"iter %d status %d ant 1 %g %g %g %g sumwt %f", 
	iter, status, in->antParm[0], in->antParm[1], in->antParm[2], in->antParm[3], sumwt); */
      Obit_log_error(err, OBIT_InfoErr,"iter %d status %d grad norm %f", 
		     iter, status, sumwt);
      ObitErrLog(err); 
    }     /* end  Diagnostics */

    /* Minimum of two iterations */
    if (iter<2) status = GSL_CONTINUE;

    /* Convergence test */    
    if (iter>1)
      status = gsl_multifit_test_delta (solver->dx, solver->x, 
					epsabs, epsrel);
  } while ((status==GSL_CONTINUE) && (iter<maxIter));
  
 /* If it didn't work - bail */
  if ((status!=GSL_SUCCESS) && (status!=GSL_CONTINUE)) {
    fprintf (stderr, "Failed, status = %s\n", gsl_strerror(status));
    goto done;
  }

  /* Get fitting results - fixed parameters will still be in the arrays */
  /* first antennas */
  for (iant=0; iant<in->nant; iant++) {
    /* Loop over parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if ((in->antFit[iant][k])  && (in->gotAnt[iant])) {
	j = in->antPNumb[iant][k];
	in->antParm[iant*4+k] = gsl_vector_get(solver->x, j);
	/* if ((k==0) || (k==2))  undo wraps */
	in->antParm[iant*4+k] = fmod(in->antParm[iant*4+k], 2*G_PI);
	if (in->antParm[iant*4+k]>G_PI) in->antParm[iant*4+k] -= 2*G_PI;
      }
    } /* end loop over parameters */
  } /* end loop over antennas */

  /* Antenna gains if needed */
  if (!in->isCircFeed && in->doFitGain) {
    for (iant=0; iant<in->nant; iant++) {
      if (in->gotAnt[iant]) {
	/* Loop over parameters */
	for (k=0; k<2; k++) {
	  /* Fitting? */
	  if (in->antGainFit[iant][k]) {
	    j = in->antGainPNumb[iant][k];
	    in->antGain[iant*2+k] = gsl_vector_get(solver->x, j);
	  }
	} /* end loop over parameters */
      } /* end if gotAnt */
    } /* end loop over antennas */
  } /* End if antenna gains */
  
  /* now source */
  for (isou=0; isou<in->nsou; isou++) {
    /* Loop over parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if (in->souFit[isou][k]) {
	j = in->souPNumb[isou][k];
	in->souParm[isou*4+k] = gsl_vector_get(solver->x, j);
      }
    } /* end loop over parameters */
  } /* end loop over sources */

  if (in->doFitRL) {
    j = in->PDPNumb;
    in->PD = gsl_vector_get(solver->x, j);
  }
  /* Errors */
  if (in->doError) {
    /* Get covariance matrix - extract diagonal terms */
    gsl_multifit_covar (solver->J, 0.0, covar);
    for (iant=0; iant<in->nant; iant++) {
      /* Loop over antenna parameters */
      for (k=0; k<4; k++) {
	/* Fitting? */
	if ((in->antFit[iant][k])  && (in->gotAnt[iant])) {
	  j = in->antPNumb[iant][k];
	  in->antErr[iant*4+k] =  sqrt(gsl_matrix_get(covar, j, j));
	} else in->antErr[iant*4+k] =  -1.0;
      } /* end loop over parameters */
    } /* end loop over antennas */
    
    /* Antenna gains if needed */
    if (!in->isCircFeed && in->doFitGain) {
      for (iant=0; iant<in->nant; iant++) {
	if (in->gotAnt[iant]) {
	  /* Loop over parameters */
	  for (k=0; k<2; k++) {
	    /* Fitting? */
	    if (in->antGainFit[iant][k]) {
	      j = in->antGainPNumb[iant][k];
	      in->antGainErr[iant*2+k] =  sqrt(gsl_matrix_get(covar, j, j));
	    } else in->antGainErr[iant*2+k] =  -1.0;
	  } /* end loop over parameters */
	} /* end if gotAnt */
      } /* end loop over antennas */
    } /* End if antenna gains */

    /* now source */
    for (isou=0; isou<in->nsou; isou++) {
      /* Loop over parameters */
      for (k=0; k<4; k++) {
	/* Fitting? */
	if (in->souFit[isou][k]) {
	  j = in->souPNumb[isou][k];
	  in->souErr[isou*4+k] =  sqrt(gsl_matrix_get(covar, j, j));
	} else in->souErr[isou*4+k] = -1.0;
      } /* end loop over parameters */
    } /* end loop over sources */

    /* Now phase difference */
    if (in->doFitRL) {
      j = in->PDPNumb;
      in->PDerr = sqrt(gsl_matrix_get(covar, j, j));
    }

    /* Count valid data */
    nvalid = 0;
    for (i=0; i<nvis*4; i++) if (in->inWt[i]>0.0) nvalid++;
    
    /* normalized Chi squares */
    if (nvalid>nparam) {
      sumwt = (ofloat)gsl_blas_dnrm2(solver->f);
      chi2Test = (sumwt*sumwt)/(nvalid-nparam);
    } else chi2Test = -1.0;
    
    in->ChiSq = chi2Test;
  } /* end do error */

    /* Diagnostics */
 done:
  if (err->prtLv>=2) {
    Obit_log_error(err, OBIT_InfoErr, "%d iter LM fit", iter);
    for (isou=0; isou<in->nsou; isou++) {
      fpol = sqrt (in->souParm[isou*4+1]*in->souParm[isou*4+1] + in->souParm[isou*4+2]*in->souParm[isou*4+2]) /
	in->souParm[isou*4+0];
      fpa = 57.296*atan2(in->souParm[isou*4+2], in->souParm[isou*4+1]);
      Obit_log_error(err, OBIT_InfoErr, 
		     "sou %3d %8.4f (%8.4f) %8.4f (%8.4f) %8.4f (%8.4f) %8.4f (%8.4f)", 
		     isou+1, in->souParm[isou*4+0], in->souErr[isou*4+0], 
		     in->souParm[isou*4+1], in->souErr[isou*4+1],
		     in->souParm[isou*4+2], in->souErr[isou*4+2], 
		     in->souParm[isou*4+3], in->souErr[isou*4+3]);
      Obit_log_error(err, OBIT_InfoErr, 
		     "        (fpol %6.4f %s %6.2f", fpol, "@", fpa);
    }
    
    for (iant=0; iant<in->nant; iant++) {
      if (in->gotAnt[iant]) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "ant %3d %8.2f (%8.2f) %8.2f (%8.2f) %8.2f (%8.2f) %8.2f (%8.2f)", 
		       iant+1, in->antParm[iant*4+0]*57.296, MAX(-1.0, in->antErr[iant*4+0]*57.296), 
		       in->antParm[iant*4+1]*57.296, MAX(-1.0, in->antErr[iant*4+1]*57.296), 
		       in->antParm[iant*4+2]*57.296, MAX(-1.0, in->antErr[iant*4+2]*57.296), 
		       in->antParm[iant*4+3]*57.296, MAX(-1.0, in->antErr[iant*4+3]*57.296));
      }
    }
    /* Phase differrence */
    if (in->doFitRL) {
      Obit_log_error(err, OBIT_InfoErr, 
		     "Phase difference %8.2f (%8.2f)",in->PD*57.296,in->PDerr*57.296);
    }

    /* Antenna gain if fitted */
    if (!in->isCircFeed && in->doFitGain) {
      for (iant=0; iant<in->nant; iant++) {
	if (in->gotAnt[iant]) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "ant %3d gain X %6.3f (%6.3f) Y %6.3f (%6.3f)", 
			 iant+1, in->antGain[iant*2+0], in->antGainErr[iant*2+0],
			 in->antGain[iant*2+1], in->antGainErr[iant*2+1]);
	}
      }
    } /* end antenna gain */
  } /* end diagnostics */
  ObitErrLog(err); 

#endif /* HAVE_GSL */ 
 return;
 } /* end doFitGSL */

/**
 * Determine Fit Poln parameters
 * Solution uses a relaxation method from Fred Schwab:  
 * Pn+1 = Pn + atan2 (dChi2/dP), (d2Chi2/dP2))  
 * for each parameter P where n or n+1 indicates a given iteration,   
 * dChi2/dP is the first partial derivative of Chi squared wrt P,  
 * d2Chi2/d2P is the second partial derivative of Chi squared wrt P,  
 * Chi2 = Sum (w (model-obs)**2)
 * \param in    Fitting object
 * \param err   Obit error stack object.
 * \return TRUE if worked, else FALSE.
 */
static gboolean doFitFast (ObitPolnCalFit *in, ObitErr *err)
{
  odouble begChi2, endChi2=0.0, iChi2, tChi2, deriv=0.0, deriv2=1.0;
  odouble difParam, difChi2=0.0, hiChi2=0.0, dChi2, d2Chi2;
  ofloat tParam, sParam, delta, sdelta = 0.01;
  ofloat ParRMS, XRMS, fpol, fpa;
  olong isou=0, iant, i, k=0, iter;
  olong pas=0, ptype=0, pnumb=-1;
  gchar *routine="ObitPolnCalFit:doFitFast";

  if (err->error) return TRUE;  /* Error exists? */
  
  iter = 0;
  while (iter<300) { 
    in->selSou = -1;
    in->selAnt = -1;
    begChi2 = GetChi2 (in->nThread, in,  polnParmUnspec, 0,
		       &ParRMS, &XRMS, NULL, NULL, err);  /* Initial value */
    /* Test for valid data */
    if (begChi2<=0.0) {
      Obit_log_error(err, OBIT_InfoErr, "No valid data");
      ObitErrLog(err); 
      return FALSE;
   }
    hiChi2 = MAX (hiChi2, begChi2);
    difParam = 0.0; ptype = -1; pnumb = -1; 
    if ((iter==0) &&(err->prtLv>=2)) {
      /*hiChi2 = begChi2;  WHAT???*/
      Obit_log_error(err, OBIT_InfoErr, "Initial Chi2=%g Parallel RMS %g Cross RMS %g", 
		     begChi2, ParRMS, XRMS);
   }

    /* Fit phase difference? */
    if (in->doFitRL) {
      iChi2 = GetChi2 (in->nThread, in, polnParmPD, k,
		       &ParRMS, &XRMS, &dChi2, &d2Chi2, err);     /* Initial value */
      sParam = in->PD;
      deriv  = dChi2;
      deriv2 = d2Chi2;
      if (fabs(deriv2)>(fabs(deriv)*1.0e-4)) {      /* Bother with this one?*/
	delta = -0.5*atan2(deriv, deriv2);
	/* Don't go wild */
	if (delta>0.) delta = MIN (delta,  20*sdelta);
	if (delta<0.) delta = MAX (delta, -20*sdelta);
	/* Loop decreasing delta until fit improves */
	for (i=0; i<10; i++) {
	  in->PD = sParam + delta;
	  tChi2 = GetChi2 (in->nThread, in, polnParmUnspec, k, 
			   &ParRMS, &XRMS, NULL, NULL, err);
	  if (tChi2<iChi2) {
	    /* Got new value? */
	    sParam += delta;
	    if (fabs(difParam)<fabs(delta)) {
	      ptype = 0; pnumb = 0; pas = 0;
	      difParam = delta;
	    }
	    break;
	  }
	  delta *= -0.7;  /* Test signs */
	  if (fabs(delta)<1.0e-6) break;
	} /* end loop decreasing */
      } 
      in->PD = sParam;   /* Set parameter to new or old value */
    } /* end  Fit phase difference */
    else in->PD = 0.0;

    /* Loop over sources */
    for (isou=0; isou<in->nsou; isou++) {
      in->selSou = isou;
      /* Loop over parameters ipol, qpol, upol, vpol*/
      for (k=0; k<4; k++) {
	/* Fitting? */
	if (in->souFit[isou][k]) {
	  iChi2 = GetChi2 (in->nThread, in, polnParmSou, k,
	     &ParRMS, &XRMS, &dChi2, &d2Chi2, err); /* Initial value */
	  sParam = in->souParm[isou*4+k];
	  deriv  = dChi2;
	  deriv2 = d2Chi2;
	  if (fabs(deriv2)>(fabs(deriv)*1.0e-4))      /* Bother with this one?*/
	    delta = -0.5*atan2(deriv, deriv2);
	  else {in->souParm[isou*4+k] = sParam; continue;}
	  if (delta>0.) delta = MIN (delta, 20*sdelta);
	  if (delta<0.) delta = MAX (delta, -20*sdelta);
	  /* Loop decreasing delta until fit improves */
	  for (i=0; i<10; i++) {
	    in->souParm[isou*4+k] = sParam + delta;
	    tChi2 = GetChi2 (in->nThread, in, polnParmUnspec, k, 
			     &ParRMS, &XRMS, NULL, NULL, err);
	    if (tChi2<iChi2) {
	      /* Got new value */
	      sParam += delta;
	      if (fabs(difParam)<fabs(delta)) {
		ptype = 1; pnumb = k; pas = isou;
		difParam = delta;
	      }
	      break;
	    }
	    delta *= -0.7;   /* Test signs */
	    if (fabs(delta)<1.0e-6) break;
	  } /* end loop decreasing */
	  in->souParm[isou*4+k] = sParam;   /* Set parameter to new or old value */
	} /* end if fitting parameter */ 
	else if (k==1) { /* Fixed Q poln? */
	  if (in->PPol[isou]>0.0)
	    in->souParm[isou*4+k] = 
	      in->PPol[isou] * in->souParm[isou*4] * cos(in->RLPhase[isou]);
	} /* end Fixed Q poln */
	else if (k==2) { /* Fixed U poln? */
	  if (in->PPol[isou]>0.0)
	    in->souParm[isou*4+k] = 
	      in->PPol[isou] * in->souParm[isou*4] * sin(in->RLPhase[isou]);
	} /* end Fixed U poln */
	if (err->error) Obit_traceback_val (err, routine, in->name, FALSE);
      } /* end loop over parameters */
    } /* end loop over sources */
    in->selSou = -1;
    
    /* Antenna gain if fitted */
    if (!in->isCircFeed && in->doFitGain) {
      for (iant=0; iant<in->nant; iant++) {
	if (!in->gotAnt[iant]) continue;  /* Have data? */
	in->selAnt = iant;
	/* Loop over parameters, gain_X, Gain_Y */
	for (k=0; k<2; k++) {
	  /* Fitting? */
	  if (in->antGainFit[iant][k]) {
	    iChi2 = GetChi2 (in->nThread, in, polnParmGain, k,
			     &ParRMS, &XRMS, &dChi2, &d2Chi2, err);  /* Initial value */
	    if (iChi2<=0.0) continue;
	    sParam = in->antGain[iant*2+k];
	    deriv  = dChi2;
	    deriv2 = d2Chi2;
	    if ((fabs(deriv2)>(fabs(deriv)*1.0e-4)) && (fabs(deriv)>(fabs(deriv2)*1.0e-3)))  /* Bother with this one?*/
	      delta = -0.5*atan2(deriv, deriv2);
	    else {in->antGain[iant*2+k] = sParam; continue;}
	    /* Loop decreasing delta until fit improves */
	    for (i=0; i<10; i++) {
	      tParam = sParam + delta;
	      in->antGain[iant*2+k] = tParam;
	      tChi2 = GetChi2 (in->nThread, in, polnParmUnspec, k, 
			       &ParRMS, &XRMS, NULL, NULL, err);
	      if (tChi2<iChi2) {
		/* Got new value */
		sParam += delta;
		if (fabs(difParam)<fabs(delta)) {
		  ptype = 2; pnumb = k; pas = iant;
		  difParam = delta;
		}
		break;
	      }
	      delta *= -0.7;  /* Test signs */
	      if (fabs(delta)<1.0e-6) break;
	    } /* end loop decreasing */
	    in->antGain[iant*2+k] = sParam;   /* Set parameter to new or old value */
	  } /* end if fitting parameter */ 
	  if (err->error) Obit_traceback_val (err, routine, in->name, FALSE);
	} /* end loop over parameters */
      } /* end loop over antennas */
    } /* end antenna gain */

    /* Don't change antenna non-gain parameters the first 3 loops */
    if (iter<3) {iter++; continue;}

    /* Loop over antennas */
    for (iant=0; iant<in->nant; iant++) {
      if (!in->gotAnt[iant]) continue;  /* Have data? */
      in->selAnt = iant;
      /* Loop over parameters 
	 0,1 ori, elip (r/x),
	 2,3 ori, elip  (l/y),
      */
      for (k=0; k<4; k++) {
	/* Fitting? */
	if (in->antFit[iant][k]) {
	  iChi2 = GetChi2 (in->nThread, in, polnParmAnt, k,
			   &ParRMS, &XRMS, &dChi2, &d2Chi2, err);  /* Initial value */
	  if (iChi2<=0.0) continue;
	  sParam = in->antParm[iant*4+k];
	  deriv  = dChi2;
	  deriv2 = d2Chi2;
	  if ((fabs(deriv2)>(fabs(deriv)*1.0e-4)) && (fabs(deriv)>(fabs(deriv2)*1.0e-3)))  /* Bother with this one?*/
	    delta = -0.5*atan2(deriv, deriv2);
	  else {in->antParm[iant*4+k] = sParam; continue;}
	  /* Loop decreasing delta until fit improves */
	  for (i=0; i<10; i++) {
	    tParam = sParam + delta;
	    in->antParm[iant*4+k] = tParam;
	    tChi2 = GetChi2 (in->nThread, in, polnParmUnspec, k, 
			     &ParRMS, &XRMS, NULL, NULL, err);
	    if (tChi2<iChi2) {
	      /* Got new value */
	      sParam += delta;
	      sParam = fmod (sParam, 2*G_PI);
	      if (sParam> G_PI) sParam -= 2*G_PI;
	      if (sParam<-G_PI) sParam += 2*G_PI;
	      if (fabs(difParam)<fabs(delta)) {
		ptype = 2; pnumb = k; pas = iant;
		difParam = delta;
	      }
	      break;
	    }
	    delta *= -0.7;  /* Test signs */
	    if (fabs(delta)<1.0e-6) break;
	  } /* end loop decreasing */
	  in->antParm[iant*4+k] = sParam;   /* Set parameter to new or old value */
	} /* end if fitting parameter */ 
	if (err->error) Obit_traceback_val (err, routine, in->name, FALSE);
      } /* end loop over parameters */
    } /* end loop over antennas */
    in->selAnt = -1;
    
    in->selSou = -1;
    in->selAnt = -1;
    endChi2 = GetChi2 (in->nThread, in, polnParmUnspec, 0, 
		       &ParRMS, &XRMS, NULL, NULL, err);  /* Final value */

    /* Convergence test */
    difChi2 = fabs(begChi2 - endChi2);

    if ((fabs(difParam)<1.0e-6) || 
      ((fabs(difParam)<5.0e-3)&&(difChi2<=1.0e-5*hiChi2))) break; 

    /* Diagnostics */
    if (err->prtLv>=4) {
      Obit_log_error(err, OBIT_InfoErr, 
		     "%d iter Chi2=%g Param %g type %d numb %d sou/ant %d dChi2 %g", 
		     iter+1, endChi2, difParam, ptype, pnumb, pas, difChi2);
      ObitErrLog(err); 
   }
    
    iter++;
  } /* end iteration loop */
  /* Diagnostics */
  if (err->prtLv>=2) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "%d iter relaxation Chi2=%g Par RMS %g X RMS %g ", 
		   iter+1, endChi2, ParRMS, XRMS);
    /* Don't give fit is using GSL */
    
    if (strncmp(in->solnType, "LM  ", 4)) {
      Obit_log_error(err, OBIT_InfoErr, "Phase difference %8.2f",  in->PD*57.296);
      for (isou=0; isou<in->nsou; isou++) {
	fpol = sqrt (in->souParm[isou*4+1]*in->souParm[isou*4+1] + in->souParm[isou*4+2]*in->souParm[isou*4+2]) /
	  in->souParm[isou*4+0];
	fpa = 57.296*atan2(in->souParm[isou*4+2], in->souParm[isou*4+1]);
	Obit_log_error(err, OBIT_InfoErr, 
		       "sou %3d %8.4f %8.4f %8.4f %8.4f (fpol %6.4f %s %6.2f)", 
		       isou+1, in->souParm[isou*4+0], in->souParm[isou*4+1], 
		       in->souParm[isou*4+2], in->souParm[isou*4+3], fpol, "@", fpa);
      }
		       
      for (iant=0; iant<in->nant; iant++)
	if (in->gotAnt[iant]) 
	  Obit_log_error(err, OBIT_InfoErr, 
			 "ant %3d %8.2f %8.2f %8.2f %8.2f", 
			 iant+1, in->antParm[iant*4+0]*57.296, in->antParm[iant*4+1]*57.296, 
			 in->antParm[iant*4+2]*57.296, in->antParm[iant*4+3]*57.296);
      /* Antenna gain if fitted */
      if (!in->isCircFeed && in->doFitGain) {
	for (iant=0; iant<in->nant; iant++) {
	  if (in->gotAnt[iant]) {
	    Obit_log_error(err, OBIT_InfoErr, 
			   "ant %3d gain X %6.3f Y %6.3f ", 
			   iant+1, in->antGain[iant*2+0], in->antGain[iant*2+1]);
	  }
	}
      } /* end antenna gain */
    } /* end if not LM */
  } /* end diagnostics */
  ObitErrLog(err); 
  return TRUE;
 } /* end doFitFast */

/**
 * Determine Chi Sq and RMS residuals of current model on selected data
 * Also optionally computes the first and second derivative wrt a given 
 * parameter which is specified by paramType, paramNumber
 * \param nThreads Number of threads to use
 * \param in           Fitting object
 * \param paramType    Parameter type
 * \li polnParmUnspec  Unspecified = don't compute derivatives
 * \li polnParmSou     Source parameter
 * \li polnParmAnt     Antenna parameter
 * \li polnParmGain    Antenna gains
 * \li polnParmPD      Phase difference
 * \param paramNumber  Parameter number,
 *                     Sou: 0=Ipol, 1=Qpol, 2=Upol, 3=VPol
 *                     Gain 0=X, 1=Y
 *                     Ant: 0= Ori R/X, 1=Elp R/X, 2= Ori L/Y, 1=Elp L/Y,       
 * \param ParRMS       [out] Parallel hand RMS
 * \param XRMS         [out] Cross hand RMS
 * \param dChi2        [out] First derivative of Chi2 wrt parameter
 *                     May be NULL if not needed
 * \param d2Chi2       [out] Second derivative of Chi2 wrt parameter
 *                     May be NULL if not needed
 * \param err          Obit error stack object.
 * \return             Chi^2
 */
static odouble GetChi2 (olong nThreads, ObitPolnCalFit *in, 
			PolnParmType paramType, olong paramNumber,
			ofloat *ParRMS, ofloat *XRMS, 
			odouble *dChi2, odouble *d2Chi2,
			ObitErr *err)
{
  odouble Chi2 = -1.0;
  ObitThreadFunc func=NULL;
  odouble sumParResid, sumXResid, sumWt, ldChi2, ld2Chi2;
  olong iTh, nPobs, nXobs, nVisPerThread;
  gboolean OK;
  gchar *routine="ObitPolnCalFit:GetChi2";

  if (err->error) return Chi2;  /* Error exists? */
  if (dChi2!=NULL)  *dChi2  = 0.0;
  if (d2Chi2!=NULL) *d2Chi2 = 1.0;

  /* Feed polarization type? */
  if (in->isCircFeed) {
    /* Circular feeds */
    func = (ObitThreadFunc)ThreadPolnFitRLChi2;
  } else {
    /* Linear feeds  */
    func = (ObitThreadFunc)ThreadPolnFitXYChi2;
  }

  /* Init threads */
  nVisPerThread = in->nvis / in->nThread;
  for (iTh=0; iTh<nThreads; iTh++) {
    in->thArgs[iTh]->selSou = in->selSou;
    in->thArgs[iTh]->selAnt = in->selAnt;
    in->thArgs[iTh]->inData = in->inData;
    in->thArgs[iTh]->inWt   = in->inWt;
    in->thArgs[iTh]->antNo  = in->antNo;
    in->thArgs[iTh]->souNo  = in->souNo;
    in->thArgs[iTh]->PD     = in->PD;
    in->thArgs[iTh]->nvis   = in->nvis;
    in->thArgs[iTh]->paramType   = paramType;
    in->thArgs[iTh]->paramNumber = paramNumber;
    in->thArgs[iTh]->lo          = iTh*nVisPerThread;     /* Zero rel first */
    in->thArgs[iTh]->hi          = (iTh+1)*nVisPerThread; /* Zero rel last */
  }
  /* Make sure do all data */
  iTh = in->nThread-1;
  in->thArgs[iTh]->hi = MAX (in->thArgs[iTh]->hi, in->nvis);

  OK = ObitThreadIterator (in->thread, nThreads, func, (gpointer)in->thArgs);
  if (!OK) {
    Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
    return Chi2;
  }

  /* Sum over threads */
  Chi2 = sumParResid = sumXResid = sumWt = ldChi2 = ld2Chi2 = 0.0;
  nPobs = nXobs = 0;
  for (iTh=0; iTh<nThreads; iTh++) {
    Chi2        += in->thArgs[iTh]->ChiSq;
    sumParResid += in->thArgs[iTh]->sumParResid;
    sumXResid   += in->thArgs[iTh]->sumXResid;
    sumWt       += in->thArgs[iTh]->sumWt;
    nPobs       += in->thArgs[iTh]->nPobs;
    nXobs       += in->thArgs[iTh]->nXobs;
    if (paramType!=polnParmUnspec) {
      ldChi2  += in->thArgs[iTh]->sumDeriv;
      ld2Chi2 += in->thArgs[iTh]->sumDeriv2;
    }
  }
  if (sumWt>0.0) Chi2 /= sumWt;

  /* output RMSes */
  *ParRMS = sqrt (sumParResid / nPobs);
  *XRMS   = sqrt (sumXResid / nXobs);
  if ((paramType!=polnParmUnspec) && (dChi2!=NULL) && (sumWt>0.0))  
    *dChi2  = ldChi2 / sumWt;
  if ((paramType!=polnParmUnspec) && (d2Chi2!=NULL) && (sumWt>0.0))  
    *d2Chi2 = ld2Chi2 / sumWt;

  return Chi2;
 } /* end GetChi2 */

/**
 * Determine average R-L phase difference
 * \param in           Fitting object
 * \param err          Obit error stack object.
 * \return             Chi^2
 */
static ofloat FitRLPhase (ObitPolnCalFit *in, ObitErr *err)
{
  ofloat PD = 0.0;
  ofloat *RLr=NULL, *RLi=NULL, *LRr=NULL, *LRi=NULL, *RLwt=NULL, *LRwt=NULL;
  ofloat *data, *wt, iipol, Re1, Re2, Im1, Im2, PPol1, PPol2, Pang1, Pang2;
  odouble sumr, sumi, sumw;
  olong i, isou, idata;
  gboolean OK;
  /*gchar *routine="ObitPolnCalFit:FitRLPhase";*/

  if (err->error) return PD;  /* Error exists? */

  /* Alias work arrays */
  RLr = (ofloat*)in->SR;
  RLi = (ofloat*)in->DR;
  LRr = (ofloat*)in->SL;
  LRi = (ofloat*)in->DL;
  RLwt = (ofloat*)in->RS;
  LRwt = (ofloat*)in->RS;

  /* Init sums */
  OK = FALSE;
  for (i=0; i<in->nsou; i++) {
    RLr[i] = RLi[i] = LRr[i] = LRi[i] = RLwt[i] = LRwt[i] = 0.0;
    /* Check for suitable calibrators */
    OK = OK || (!in->souFit[i][1] && !in->souFit[i][2] && 
		(in->PPol[i]>0.0001));
  }

  /* Any suitable calibrators? */
  if (!OK) return PD;
  wt   = in->inWt;
  data = in->inData;

  /* Loop over data */
  for (idata=0; idata<in->nvis; idata++) {
    isou  = MAX (0, in->souNo[idata]);    /* Source number */
    /* This one usable? */
    if (in->souFit[isou][1] || in->souFit[isou][2] || 
	(in->PPol[isou]<=0.0001)) continue;
    /* Normalize data py 1/ipol */
    iipol = 1.0 / in->souParm[isou*4];
    if (wt[idata*4+2]>0.0) {
      RLr[isou]  += data[idata*10+6] * iipol * wt[idata*4+2];
      RLi[isou]  += data[idata*10+7] * iipol * wt[idata*4+2];
      RLwt[isou] += wt[idata*4+2];
    }
    if (wt[idata*4+3]>0.0) {
      LRr[isou] += data[idata*10+8] * iipol * wt[idata*4+3];
      LRi[isou] += data[idata*10+9] * iipol * wt[idata*4+3];
      LRwt[isou] += wt[idata*4+3];
    }
  } /* end loop over data */

  /* Loop over sources differences phase with model - 
     weight average by polarized flux,
     average as real and imaginary parts */
  sumr = sumi = sumw = 0.0;
  /* Signs here need to be checked */
  for (isou=0; isou<in->nsou; isou++) {
    if (in->souFit[isou][1] || in->souFit[isou][2] || 
	(in->PPol[isou]<=0.0001)) continue;
    /* Get weighted average RL and LR per source */
    if (RLwt[isou]>0.0) {
      Re1   = RLr[isou] / RLwt[isou];
      Im1   = RLi[isou] / RLwt[isou];
      PPol1 = in->souParm[isou*4] * sqrt (Re1*Re1 + Im1*Im1);
    } else {Re1 = 1.0; Im1=0.0; PPol1 = 0.0;}
    Pang1 = atan2 (Im1, Re1);
    sumr += PPol1 * cos (Pang1 - in->RLPhase[isou]);
    sumi += PPol1 * sin (Pang1 - in->RLPhase[isou]);
    sumw += PPol1;
    if (LRwt[isou]>0.0) {
      Re2   = LRr[isou] / LRwt[isou];
      Im2   = LRi[isou] / LRwt[isou];
      PPol2 = in->souParm[isou*4] * sqrt (Re2*Re2 + Im2*Im2);
    } else {Re2 = 1.0; Im2=0.0; PPol2 = 0.0;}
    Pang2 = atan2 (Im2, Re2);
    sumr += PPol2 * cos (Pang2 + in->RLPhase[isou]);
    sumi -= PPol2 * sin (Pang2 + in->RLPhase[isou]);
    sumw += PPol2;
    }

  if ((sumr!=0.0) || (sumi!=0.0))
    PD = atan2 (sumi, sumr); /* Weighted phase difference */

  return PD;
} /* end FitRLPhase */

/**
 * Circular polarization version.
 * Threaded Chi**2 evaluator for polarization fitting
 * Evaluates sum [(model-observed) / sigma] and derivatives
 * Parallel hands count 0.3 of cross in sums
 * If selAnt or selSou are set, then only data involving that 
 * source or antenna (or ref ant) is included.
 * \param arg   PolnFitArg pointer with elements:
 * \li lo       First 0-rel datum in Data/Wt arrays
 * \li hi       Last 0-rel datum  in Data/Wt arrays
 * \li selAnt   selected antenna, -1-> all
 * \li selSou   selected source,  -1> all
 * \li paramType       Parameter type
 * \li polnParmUnspec  Unspecified = don't compute derivatives
 * \li polnParmAnt     Antenna parameter
 * \li polnParmSou     Source parameter
 * \li polnParmPD      Phase difference
 * \li paramNumber  Parameter number,
 *                     Sou: 0=Ipol, 1=Qpol, 2=Upol, 3=VPol
                       Ant: 0= Ori R/X, 1=Elp R/X, 2= Ori L/Y, 1=Elp L/Y,       
 * \li ChiSq        [out] computed Chi2
 * \li nPobs        [out] Number of valid parallel measurements
 * \li nXobs        [out] Number of valid cross pol measurements
 * \li ParRMS       [out] Parallel hand RMS
 * \li sumParResid  [out] Cross hand RMS
 * \li sumXResid    [out] First derivative of Chi2 wrt parameter
 * \li d2Chi2       [out] Second derivative of Chi2 wrt parameter
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadPolnFitRLChi2 (gpointer arg)
{
  PolnFitArg *args = (PolnFitArg*)arg;
  odouble    *antParm   = args->antParm;
  gboolean   **antFit   = args->antFit;
  odouble    *souParm   = args->souParm;
  gboolean   **souFit   = args->souFit;
  ofloat     *data      = args->inData;
  ofloat     *wt        = args->inWt;
  odouble    *SR        = args->SR;
  odouble    *DR        = args->DR;
  odouble    *SL        = args->SL;
  odouble    *DL        = args->DL;
  dcomplex   *PR        = args->PR;
  dcomplex   *PRc       = args->PRc;
  dcomplex   *PL        = args->PL;
  dcomplex   *PLc       = args->PLc;
  dcomplex   *RS        = args->RS;
  dcomplex   *RD        = args->RD;
  dcomplex   *LS        = args->LS;
  dcomplex   *LD        = args->LD;
  dcomplex   *RSc       = args->RSc;
  dcomplex   *RDc       = args->RDc;
  dcomplex   *LSc       = args->LSc;
  dcomplex   *LDc       = args->LDc;
  PolnParmType paramType = args->paramType;
  olong paramNumber      = args->paramNumber;
  olong selSou           = args->selSou;
  olong selAnt           = args->selAnt;

  odouble ipol=0.0, qpol=0.0, upol=0.0, vpol=0.0;
  odouble residR=0.0, residI=0.0, isigma=0.0;
  odouble sumParResid, sumXResid;
  ofloat PD, chi1, chi2;
  olong nPobs, nXobs, ia1, ia2, isou, idata, isouLast=-999;
  gboolean isAnt1, isAnt2;
  size_t i;
  odouble sum=0.0, sumwt=0.0, sumd, sumd2;

  dcomplex PRref, PLref, PPRL, PPLR, PA1, PA2, PA1c, PA2c;
  dcomplex ct1, ct2, ct3, ct4, ct5, ct6, dt1, dt2;
  dcomplex S[4], VRR, VRL, VLR, VLL, MC1, MC2, MC3, MC4, DFDP, DFDP2;
  ofloat root2;

  COMPLEX_SET (S[0], 0.0, 0.0);  /* Initialize poln vector */
  COMPLEX_SET (S[1], 0.0, 0.0);
  COMPLEX_SET (S[2], 0.0, 0.0);
  COMPLEX_SET (S[3], 0.0, 0.0);
  COMPLEX_SET (MC1, 0.0, 0.0);  /* Other stuff */
  COMPLEX_SET (MC2, 0.0, 0.0);
  COMPLEX_SET (MC3, 0.0, 0.0);
  COMPLEX_SET (MC4, 0.0, 0.0);
  COMPLEX_SET (VRR, 0.0, 0.0);
  COMPLEX_SET (VLL, 0.0, 0.0);
  COMPLEX_SET (VLR, 0.0, 0.0);
  COMPLEX_SET (VRL, 0.0, 0.0);
  COMPLEX_SET (DFDP,  0.0, 0.0);
  COMPLEX_SET (DFDP2, 0.0, 0.0);
  COMPLEX_SET (dt1, 0.0, 0.0);
  COMPLEX_SET (dt2, 0.0, 0.0);
 
  /* RMS sums and counts */
  sumParResid = sumXResid = 0.0;
  sumd = sumd2 = 0.0;
  nPobs = nXobs = 0;
  /* R-L phase difference  at reference antenna */
  if (args->doFitRL) {
    PD = args->PD;
  } else PD = 0.0;

  /* Injest model, factorize into antenna components - 
     data in order Orientation R/X, Elipticity R/X, Orientation L/Y, Elipticity L/Y */
  root2 = 1.0 / sqrt(2.0);
  /* Elipticity, Orientation terms */
  for (i=0; i<args->nant; i++) {
    SR[i] = cos(antParm[i*4+1]) + sin(antParm[i*4+1]);
    DR[i] = cos(antParm[i*4+1]) - sin(antParm[i*4+1]);
    SL[i] = cos(antParm[i*4+3]) + sin(antParm[i*4+3]);
    DL[i] = cos(antParm[i*4+3]) - sin(antParm[i*4+3]);
    COMPLEX_SET  (RS[i], root2*SR[i], 0.);
    COMPLEX_SET  (ct1,   root2*DR[i], 0.);
    COMPLEX_EXP  (PR[i], 2*antParm[i*4+0]);
    COMPLEX_CONJUGATE (PRc[i], PR[i]);
    COMPLEX_MUL2 (RD[i], ct1, PR[i]);
    COMPLEX_SET  (ct1, root2*SL[i], 0.);
    COMPLEX_EXP  (PL[i], -2*antParm[i*4+2]);
    COMPLEX_CONJUGATE (PLc[i], PL[i]);
    COMPLEX_MUL2 (LS[i], ct1, PL[i]);
    COMPLEX_SET  (LD[i], root2*DL[i], 0.);
    COMPLEX_CONJUGATE (RSc[i], RS[i]);
    COMPLEX_CONJUGATE (RDc[i], RD[i]);
    COMPLEX_CONJUGATE (LSc[i], LS[i]);
    COMPLEX_CONJUGATE (LDc[i], LD[i]);
  }

  /* Reference antenna phase terms */
  COMPLEX_EXP (PRref,  antParm[(args->refAnt-1)*4+0]);
  COMPLEX_EXP (PLref, -antParm[(args->refAnt-1)*4+2]+PD);
  COMPLEX_CONJUGATE (ct1, PLref);
  COMPLEX_MUL2(PPRL, PRref, ct1);
  COMPLEX_CONJUGATE (ct1, PRref);
  COMPLEX_MUL2(PPLR, PLref, ct1);

  /* Loop over data */
  i = 0;
  for (idata=args->lo; idata<args->hi; idata++) {
    /* Parallactic angle terms */
    chi1  = data[idata*10+0];   /* parallactic angle ant 1 */
    chi2  = data[idata*10+1];   /* parallactic angle ant 2 */
    COMPLEX_EXP (PA1,2*chi1);
    COMPLEX_EXP (PA2,2*chi2);
    COMPLEX_CONJUGATE (PA1c, PA1);
    COMPLEX_CONJUGATE (PA2c, PA2);

    isou  = MAX (0, args->souNo[idata]);    /* Source number */
    /* Selected source? */
    if ((selSou>=0) && (selSou!=isou)) continue;

    /* New source? get parameters */
    if (isou!=isouLast) {
      isouLast = isou;
      /* Source parameters */
      ipol = souParm[isou*4+0];      /* Fitting or fixed? */
      if (args->souFit[isou][1]) 
	qpol = souParm[isou*4+1];
      else
	qpol = args->PPol[isou]*ipol*cos(args->RLPhase[isou]);
      if (args->souFit[isou][2]) 
	upol = souParm[isou*4+2];
      else
	upol = args->PPol[isou]*ipol*sin(args->RLPhase[isou]);
      vpol = souParm[isou*4+3];
      /* Complex Stokes array */
      COMPLEX_SET (S[0], ipol+vpol, 0.0);
      COMPLEX_SET (S[1], qpol,  upol);
      COMPLEX_SET (S[2], qpol, -upol);
      COMPLEX_SET (S[3], ipol-vpol, 0.0);
    }

    /* Antenna parameters (0 ref) */
    ia1    = args->antNo[idata*2+0];
    ia2    = args->antNo[idata*2+1]; 
    /* Selected source? */
    if ((selAnt>=0) && 
	((selAnt!=ia1) && (ia1!=args->refAnt)) && 
	 (selAnt!=ia2) && (ia2!=args->refAnt)) continue;
  /* Which antenna is the selected one in the baseline? */
  if (selAnt==ia1) {isAnt1 = TRUE; isAnt2 = FALSE;}
  else if (selAnt==ia2) {isAnt2 = TRUE; isAnt1 = FALSE;}
  else {isAnt1 = FALSE; isAnt2 = FALSE;}  /* Only refant */

  /* Calculate residals - Note different order for LL */
  /*      RR */
  if (wt[idata*4]>0.0) {
    isigma = 0.3*wt[idata*4]; /* downweight parallel */
    
    /* VRR = S[0] * RS[ia1] * RSc[ia2] +        
             S[1] * RS[ia1] * RDc[ia2] * PA2c + 
	     S[2] * RD[ia1] * RSc[ia2] * PA1  + 
	     S[3] * RD[ia1] * RDc[ia2] * PA1  * PA2c; */
    COMPLEX_MUL2 (MC1, RS[ia1], RSc[ia2]);
    COMPLEX_MUL2 (VRR, S[0], MC1);
    COMPLEX_MUL3 (MC2, RS[ia1], RDc[ia2],  PA2c);
    COMPLEX_MUL2 (ct1, S[1], MC2);
    COMPLEX_ADD2 (VRR, VRR,  ct1);
    COMPLEX_MUL3 (MC3, RD[ia1], RSc[ia2], PA1);
    COMPLEX_MUL2 (ct1, S[2], MC3);
    COMPLEX_ADD2 (VRR, VRR,  ct1);
    COMPLEX_MUL4 (MC4, RD[ia1], RDc[ia2], PA1, PA2c);
    COMPLEX_MUL2 (ct1, S[3], MC4);
    COMPLEX_ADD2 (VRR, VRR,  ct1);
    residR = VRR.real - data[idata*10+2];
    sum += isigma * residR * residR; sumwt += isigma;
    residI = VRR.imag - data[idata*10+3];
    sum += isigma * residI * residI; sumwt += isigma;
    nPobs++;
    sumParResid += residR * residR + residI * residI;
    /* DEBUG 
       if ((paramType==polnParmUnspec) && ((fabs(residR)>10.0) || ((fabs(residI)>10.0))))
       fprintf (stderr, "vis %d RR %d-%d %f %f\n", idata, ia1, ia2, residR, residI);
       end DEBUG */
    /* Derivatives */
    if (paramType==polnParmAnt) {         /* Antenna parameters */
      /* Default partials */
      COMPLEX_SET (DFDP,  0.0, 0.0);
      COMPLEX_SET (DFDP2, 0.0, 0.0);
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* RR wrt Or */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* (0, 2 i) * (S[2]*MC3 + S[3]*MC4) */
	    COMPLEX_MUL2(ct1, S[2], MC3);
	    COMPLEX_MUL2(ct2, S[3], MC4);
	    COMPLEX_ADD2(ct3, ct1, ct2);
	    COMPLEX_SET (ct1, 0.0, 2.0);
	    COMPLEX_MUL2(DFDP,  ct1, ct3);
	    COMPLEX_MUL2(DFDP2, ct1, DFDP);
	  }
	} else if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* (0,-2 i) * (S[1]*MC2 + S[3]*MC4) */
	    COMPLEX_MUL2(ct1, S[1], MC2);
	    COMPLEX_MUL2(ct2, S[3], MC4);
	    COMPLEX_ADD2(ct3, ct1, ct2);
	    COMPLEX_SET (ct1, 0.0, -2.0);
	    COMPLEX_MUL2(DFDP,  ct1, ct3);
	    COMPLEX_MUL2(DFDP2, ct1, DFDP);
	  }
	} 
	break;
      case 1:     /* RR wrt Er */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* part = r2 * DR[ia1] * (S[0] * RSc[ia2] + S[1] * RDc[ia2] * PA2c) -
	              r2 * SR[ia1] * PR[ia1] * (S[2] * RSc[ia2] * PA1 + 
		                                S[3] * RDc[ia2] * PA1  * PA2c) */
	    COMPLEX_MUL2(ct1, S[0], RSc[ia2]);
	    COMPLEX_MUL3(ct2, S[1], RDc[ia2], PA2c);
	    COMPLEX_ADD2(dt1, ct1, ct2);
	    COMPLEX_SET (ct4, root2*DR[ia1], 0.0);
	    COMPLEX_MUL2(ct5, ct4, dt1);
	    COMPLEX_MUL3(ct1, S[2], RSc[ia2], PA1);
	    COMPLEX_MUL4(ct2, S[3], RDc[ia2], PA1, PA2c);
	    COMPLEX_ADD2(dt2, ct1, ct2);
	    COMPLEX_SET (ct4, root2*SR[ia1], 0.0);
	    COMPLEX_MUL3(ct6, ct4, PR[ia1], dt2);
	    COMPLEX_SUB (DFDP, ct5, ct6);
	    /* part2 = -r2 * SR[ia1] * (S[0] * RSc[ia2] + S[1] * RDc[ia2] * PA2c) -
	                r2 * DR[ia1] * PR[ia1] * (S[2] * RSc[ia2] * PA1 + 
                                                  S[3] * RDc[ia2] * PA1  * PA2c) */
	    COMPLEX_SET (ct4, -root2*SR[ia1], 0.0);
	    COMPLEX_MUL2(ct5, ct4, dt1);
	    COMPLEX_SET (ct4, root2*DR[ia1], 0.0);
	    COMPLEX_MUL3(ct6, ct4, PR[ia1], dt2);
	    COMPLEX_SUB (DFDP2, ct5, ct6);
	  }
	} else if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* part = r2 * DR[ia2] * (S[0] * RS[ia1] + S[2] * RD[ia1] * PA1) -
	              r2 * SR[ia2] * PRc[ia2] * PA2c  * (S[1] * RS[ia1] * PA2c w+ 
                                                 S[3] * RD[ia1] * PA1) */
	    COMPLEX_MUL2(ct1, S[0], RS[ia1]);
	    COMPLEX_MUL3(ct2, S[2], RD[ia1], PA1);
	    COMPLEX_ADD2(dt2, ct1, ct2);
	    COMPLEX_SET (ct4, root2*DR[ia2], 0.0);
	    COMPLEX_MUL2(ct5, ct4, dt1);
	    COMPLEX_MUL2(ct1, S[1], RS[ia1]);
	    COMPLEX_MUL3(ct2, S[3], RD[ia1], PA1);
	    COMPLEX_ADD2(dt2, ct1, ct2);
	    COMPLEX_SET (ct4, root2*SR[ia2], 0.0);
	    COMPLEX_MUL4(ct6, ct4, PRc[ia2], PA2c, dt2);
	    COMPLEX_SUB (DFDP, ct5, ct6);
	    /* part2 = -r2 * SR[ia2] * (S[0] * RS[ia1] + S[2] * RD[ia1] * PA1) -
	                r2 * DR[ia2] * PRc[ia2] * PA2c * (S[1] * RS[ia1] + 
                                                          S[3] * RD[ia1] * PA1) */
	    COMPLEX_SET (ct4, -root2*SR[ia2], 0.0);
	    COMPLEX_MUL2(ct5, ct4, dt1);
	    COMPLEX_SET (ct4, root2*DR[ia2], 0.0);
	    COMPLEX_MUL4(ct6, ct4, PRc[ia2], PA2c, dt2);
	    COMPLEX_SUB (DFDP2, ct5, ct6);
	  }
	}
	break;
      case 2:     /* RR wrt Ol - nope */
	break;
      case 3:     /* RR wrt El - nope */
	break;
      default:
	break;
      }; /* end antenna parameter switch */
      /* end antenna param */
    } else if (paramType==polnParmSou) {   /* Source parameters */
      /* Default partials  */
      COMPLEX_SET (DFDP,  0.0, 0.0);
      COMPLEX_SET (DFDP2, 0.0, 0.0);
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* RR wrt I */
	if (souFit[isou][paramNumber]) {
	  /* part = MC1 + MC4 */
	  COMPLEX_ADD2(DFDP, MC1, MC4);
	}
	break;
      case 1:     /* RR wrt QPol */
	if (souFit[isou][paramNumber]) {
	  /* part = MC2 + MC3 */
	  COMPLEX_ADD2(DFDP, MC2, MC3);
	}
	break;
      case 2:     /* RR wrt UPol */
	if (souFit[isou][paramNumber]) {
	  /* part = i (MC2 - MC3) */
	  COMPLEX_SUB (ct3, MC2, MC3);
	  COMPLEX_SET(ct4, 0.0, 1.0);
	  COMPLEX_MUL2(DFDP, ct4, ct3);
	}
	break;
      case 3:     /* RR wrt Vpol */
	if (souFit[isou][paramNumber]) {
	  /* part = MC1 - MC4 */
	  COMPLEX_SUB (DFDP, MC1, MC4);
	}
	break;  
      default:
	break;
      }; /* end source parameter switch */
      /* end source param */
    } else if (paramType==polnParmPD) {   /* R-L Phase difference */
      /* No dependence on RR */
      COMPLEX_SET (DFDP,  0.0, 0.0);
      COMPLEX_SET (DFDP2, 0.0, 0.0);      
    } /* end R-L phase difference */
    /* Accumulate partials */
    if (paramType!=polnParmUnspec) {
      sumd  += 2.0 * isigma * (residR*DFDP.real + residI*DFDP.imag);
      sumd2 += 2.0 * isigma * (DFDP.real*DFDP.real + DFDP.imag*DFDP.imag +
			       residR*DFDP2.real + residI*DFDP2.imag);
    } /* end set partials */
  } /* end valid data */
  
    /*      LL */
  if (wt[idata*4+1]>0.0) {
    isigma = 0.3*wt[idata*4+1];  /* downweight parallel */
    /* VLL = S[0] * LS[ia1] * LSc[ia2] * PA1c * PA2 +	
             S[1] * LS[ia1] * LDc[ia2] * PA1c +
	     S[2] * LD[ia1] * LSc[ia2] * PA2  +
	     S[3] * LD[ia1] * LDc[ia2]; */
    COMPLEX_MUL4 (MC1, LS[ia1], LSc[ia2], PA1c, PA2);
    COMPLEX_MUL2 (VLL, S[0], MC1);
    COMPLEX_MUL3 (MC2, LS[ia1], LDc[ia2], PA1c);
    COMPLEX_MUL2 (ct1, S[1], MC2);
    COMPLEX_ADD2 (VLL, VLL,  ct1);
    COMPLEX_MUL3 (MC3, LD[ia1], LSc[ia2], PA2);
    COMPLEX_MUL2 (ct1, S[2], MC3);
    COMPLEX_ADD2 (VLL, VLL,  ct1);
    COMPLEX_MUL2 (MC4, LD[ia1], LDc[ia2]);
    COMPLEX_MUL2 (ct1, S[3], MC4);
    COMPLEX_ADD2 (VLL, VLL,  ct1);
    residR = VLL.real - data[idata*10+4];
    sum += isigma * residR * residR; sumwt += isigma; 
    residI = VLL.imag - data[idata*10+5];
    sum += isigma * residI * residI; sumwt += isigma; 
    nPobs++;
    sumParResid += residR * residR + residI * residI;
    /* DEBUG 
       if ((paramType==polnParmUnspec) && ((fabs(residR)>5.0) || ((fabs(residI)>5.0))))
       fprintf (stderr, "vis %d LL %d-%d %f %f\n", idata, ia1, ia2, residR, residI);
       end DEBUG */
    /* Derivatives */
    if (paramType==polnParmAnt) {         /* Antenna parameters */
      /* Default partials */
      COMPLEX_SET (DFDP,  0.0, 0.0);
      COMPLEX_SET (DFDP2, 0.0, 0.0);
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* LL wrt Or - nope */
	break;
      case 1:     /* LL wrt Er - nope*/
	break;
      case 2:     /* LL wrt Ol */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* (0,-2 i) * (S[0]*MC1 + S[1]*MC2) */
	    COMPLEX_MUL2(ct1, S[0], MC1);
	    COMPLEX_MUL2(ct2, S[1], MC2);
	    COMPLEX_ADD2(ct3, ct1, ct2);
	    COMPLEX_SET (ct1, 0.0, -2.0);
	    COMPLEX_MUL2(DFDP,  ct1, ct3);
	    COMPLEX_MUL2(DFDP2, ct1, DFDP);
	  }
	} else if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* (0, 2 i) * (S[0]*MC1 + S[2]*MC3) */
	    COMPLEX_MUL2(ct1, S[0], MC1);
	    COMPLEX_MUL2(ct2, S[2], MC3);
	    COMPLEX_ADD2(ct3, ct1, ct2);
	    COMPLEX_SET (ct1, 0.0, 2.0);
	    COMPLEX_MUL2(DFDP,  ct1, ct3);
	    COMPLEX_MUL2(DFDP2, ct1, DFDP);
	  }
	}
	break;
      case 3:     /* LL wrt El */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* part = r2 * DL[ia1] * PL[ia1] * PA1c * (S[0] * LSc[ia2] * PA2 + 
	                                               S[1] * LDc[ia2) -
	              r2 * SL[ia1] * (S[2] * LSc[ia2] * PA2  + S[3] * LDc[ia2]) */
	    COMPLEX_MUL3(ct1, S[0], LSc[ia2], PA2);
	    COMPLEX_MUL2(ct2, S[1], LDc[ia2]);
	    COMPLEX_ADD2(dt1, ct1, ct2);
	    COMPLEX_SET (ct4, root2*DL[ia1], 0.0);
	    COMPLEX_MUL4(ct5, ct4, PL[ia1], PA1c, dt1);
	    COMPLEX_MUL3(ct1, S[2], LSc[ia2], PA2);
	    COMPLEX_MUL2(ct2, S[3], LDc[ia2]);
	    COMPLEX_ADD2(dt2, ct1, ct2);
	    COMPLEX_SET (ct4, root2*SL[ia1], 0.0);
	    COMPLEX_MUL2(ct6, ct4, dt2);
	    COMPLEX_SUB (DFDP, ct5, ct6);
	    /* part2 = -r2 * SL[ia1] * PL[ia1] * PA1c * (S[0] * LSc[ia2] * PA2 + 
                                                         S[1] * LDc[ia2]) -
	                r2 * DL[ia1] * (S[2] * LSc[ia2] * PA2  + S[3] * LDc[ia2]) */
	    COMPLEX_SET (ct4, -root2*SL[ia1], 0.0);
	    COMPLEX_MUL4(ct5, ct4, PL[ia1], PA1c, dt1);
	    COMPLEX_SET (ct4, root2*DL[ia1], 0.0);
	    COMPLEX_MUL2(ct6, ct4, dt2);
	    COMPLEX_SUB (DFDP2, ct5, ct6);
	  }
	} else if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* part = r2 * DL[ia2] * PLc[ia2] * PA2 * (S[0] * LS[ia1] * PA1c  + 
	                                               S[2] * LD[ia1]) -
	              r2 * SL[ia2] * (S[1] * LS[ia1] * PA1c + S[3] * LD[ia1]) */
	    COMPLEX_MUL3(ct1, S[0], LS[ia1], PA1c);
	    COMPLEX_MUL2(ct2, S[2], LD[ia1]);
	    COMPLEX_ADD2(dt1, ct1, ct2);
	    COMPLEX_SET (ct4, root2*DL[ia2], 0.0);
	    COMPLEX_MUL4(ct5, ct4, PLc[ia2], PA2, dt1);
	    COMPLEX_MUL3(ct1, S[1], LS[ia1], PA1c);
	    COMPLEX_MUL2(ct2, S[3], LD[ia1]);
	    COMPLEX_ADD2(dt2, ct1, ct2);
	    COMPLEX_SET (ct4, root2*SL[ia2], 0.0);
	    COMPLEX_MUL2(ct6, ct4, dt2);
	    COMPLEX_SUB (DFDP, ct5, ct6);
	    /* part2 = -r2 * SL[ia2] * PLc[ia2] * PA2 * (S[0] * LS[ia1] * PA1c  + 
	                                                 S[2] * LD[ia1) -
	                r2 * DL[ia2] * (S[1] * LS[ia1] * PA1c + S[3] * LD[ia1]) */
	    COMPLEX_SET (ct4, -root2*SL[ia2], 0.0);
	    COMPLEX_MUL4(ct5, ct4, PLc[ia2], PA2, dt1);
	    COMPLEX_SET (ct4, root2*DL[ia2], 0.0);
	    COMPLEX_MUL2(ct6, ct4, dt2);
	    COMPLEX_SUB (DFDP2, ct5, ct6);
	  }
	}
	break;
      default:
	break;
      }; /* end antenna parameter switch */
      /* end antenna param */
    } else if (paramType==polnParmSou) {   /* Source parameters */
      /* Default partials  */
      COMPLEX_SET (DFDP,  0.0, 0.0);
      COMPLEX_SET (DFDP2, 0.0, 0.0);
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* LL wrt I */
	if (souFit[isou][paramNumber]) {
	  /* part = MC1 + MC4 */
	  COMPLEX_ADD2(DFDP, MC1, MC4);
	}
	break;
      case 1:     /* LL wrt QPol */
	if (souFit[isou][paramNumber]) {
	  /* part = MC2 + MC3 */
	  COMPLEX_ADD2(DFDP, MC2, MC3);
	}
	break;
      case 2:     /* LL wrt UPol */
	if (souFit[isou][paramNumber]) {
	  /* part = i (MC2 - MC3) */
	  COMPLEX_SUB (ct3, MC2, MC3);
	  COMPLEX_SET(ct4, 0.0, 1.0);
	  COMPLEX_MUL2(DFDP, ct4, ct3);
	}
	break;
      case 3:     /* LL wrt Vpol */
	if (souFit[isou][paramNumber]) {
	  /* part = MC1 - MC4 */
	  COMPLEX_SUB (DFDP, MC1, MC4);
	}
	break;  
      default:
	break;
      }; /* end source parameter switch */
      /* end source param */
    } else if (paramType==polnParmPD) {   /* R-L Phase difference */
      /* No dependence on LL */
      COMPLEX_SET (DFDP,  0.0, 0.0);
      COMPLEX_SET (DFDP2, 0.0, 0.0);      
    } /* end R-L phase difference */
    /* Accumulate partials */
    if (paramType!=polnParmUnspec) {
      sumd  += 2.0 * isigma * (residR*DFDP.real + residI*DFDP.imag);
      sumd2 += 2.0 * isigma * (DFDP.real*DFDP.real + DFDP.imag*DFDP.imag +
			       residR*DFDP2.real + residI*DFDP2.imag);
    } /* end set partials */
  } /* end valid data */
  
    /* 	    RL */
  if (wt[idata*4+2]>0.0) {
    isigma = wt[idata*4+2];
    /* VRL = PPRL * S[0] * RS[ia1] * LSc[ia2] * PA2 +
             PPRL * S[1] * RS[ia1] * LDc[ia2] + 
	     PPRL * S[2] * RD[ia1] * LSc[ia2] * PA1 * PA2 +
	     PPRL * S[3] * RD[ia1] * LDc[ia2] * PA1; */
    COMPLEX_MUL4 (MC1, PPRL, RS[ia1], LSc[ia2], PA2);
    COMPLEX_MUL2 (VRL, S[0], MC1);
    COMPLEX_MUL3 (MC2, PPRL, RS[ia1], LDc[ia2]);
    COMPLEX_MUL2 (ct1, S[1], MC2);
    COMPLEX_ADD2 (VRL, VRL,  ct1);
    COMPLEX_MUL5 (MC3, PPRL, RD[ia1], LSc[ia2],  PA1,  PA2);
    COMPLEX_MUL2 (ct1, S[2], MC3);
    COMPLEX_ADD2 (VRL, VRL,  ct1);
    COMPLEX_MUL4 (MC4, PPRL, RD[ia1], LDc[ia2],  PA1);
    COMPLEX_MUL2 (ct1, S[3], MC4);
    COMPLEX_ADD2 (VRL, VRL,  ct1);
    residR = VRL.real - data[idata*10+6];
    sum += isigma * residR * residR; sumwt += isigma;
    residI = VRL.imag - data[idata*10+7];
    sum += isigma * residI * residI; sumwt += isigma;
    nXobs++;
    sumXResid += residR * residR + residI * residI;
    /* Derivatives */
    if (paramType==polnParmAnt) {         /* Antenna parameters */
      /* Default partials  */
      COMPLEX_SET (DFDP,  0.0, 0.0);
      COMPLEX_SET (DFDP2, 0.0, 0.0);
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* RL wrt Or */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* (0, 2 i) * (S[2]*MC3 + S[3]*MC4) */
	    COMPLEX_MUL2(ct1, S[2], MC3);
	    COMPLEX_MUL2(ct2, S[3], MC4);
	    COMPLEX_ADD2(ct3, ct1, ct2);
	    COMPLEX_SET (ct1, 0.0,  2.0);
	    COMPLEX_MUL2(DFDP, ct1, ct3);
	    COMPLEX_MUL2(DFDP2, ct1, DFDP);
	  }
	} 
	/* If ia1==refant) */
	if (ia1==args->refAnt) {
	  COMPLEX_SET (ct1, 0.0,  1.0);
	  COMPLEX_MUL2(ct2, ct1, VRL);
	  COMPLEX_ADD2(DFDP, DFDP, ct2);
	  COMPLEX_MUL2(ct2, ct1, DFDP);
	  COMPLEX_ADD2(DFDP2, DFDP2, ct2);
	}
	break;
      case 1:     /* RL wrt Er */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* part = r2 * DR[ia1] * PPRL * (S[0] * LSc[ia2] * PA2 +S[1] * LDc[ia2]) -
	              r2 * SR[ia1] * PR[ia1] * PPRL * PA1 * (S[2] * LSc[ia2] * PA2 +
	                                                     S[3] * LDc[ia2]) */
	    COMPLEX_MUL3(ct1, S[0], LSc[ia2], PA2);
	    COMPLEX_MUL2(ct2, S[1], LDc[ia2]);
	    COMPLEX_ADD2(dt1, ct1, ct2);
	    COMPLEX_SET (ct4, root2*DR[ia1], 0.0);
	    COMPLEX_MUL3(ct5, ct4, PPRL, dt1);
	    COMPLEX_MUL3(ct1, S[2], LSc[ia2], PA2);
	    COMPLEX_MUL2(ct2, S[3], LDc[ia2]);
	    COMPLEX_ADD2(dt2, ct1, ct2);
	    COMPLEX_SET (ct4, root2*SR[ia1], 0.0);
	    COMPLEX_MUL5(ct6, ct4, PR[ia1], PPRL, PA1, dt2);
	    COMPLEX_SUB (DFDP, ct5, ct6);
	    /* part2 = -r2 * SR[ia1] * PPRL * (S[0] * LSc[ia2] * PA2 + S[1] * LDc[ia2]) -
	                r2 * DR[ia1] * PR[ia1] * PPRL * PA1 * (S[2] * LSc[ia2] * PA2 + 
	                                                       S[3] * LDc[ia2] ) */
	    COMPLEX_SET (ct4, -root2*SR[ia1], 0.0);
	    COMPLEX_MUL3(ct5, ct4, PPRL, dt1);
	    COMPLEX_SET (ct4, root2*DR[ia1], 0.0);
	    COMPLEX_MUL5(ct6, ct4, PR[ia1], PPRL, PA1, dt2);
	    COMPLEX_SUB (DFDP2, ct5, ct6);
	  }
	}
	break;
      case 2:     /* RL wrt Ol */
	if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* (0, 2 i) * (S[0]*MC1 + S[2]*MC3) */
	    COMPLEX_MUL2(ct1, S[0], MC1);
	    COMPLEX_MUL2(ct2, S[2], MC3);
	    COMPLEX_ADD2(ct3, ct1, ct2);
	    COMPLEX_SET (ct1, 0.0,  2.0);
	    COMPLEX_MUL2(DFDP, ct1, ct3);
	    COMPLEX_MUL2(DFDP2, ct1, DFDP);
	  }
	} 
	/* If ia2==refant) */
	if (ia2==args->refAnt) {
	  COMPLEX_SET (ct1, 0.0, 1.0);
	  COMPLEX_MUL2(ct2, ct1, VRL);
	  COMPLEX_ADD2(DFDP, DFDP, ct2);
	  COMPLEX_MUL2(ct2, ct1, DFDP);
	  COMPLEX_ADD2(DFDP2, DFDP2, ct2);
	}
	break;
      case 3:     /* RL wrt El */
	if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* part = r2 * DL[ia2] * PLc[ia2] * PPRL * PA2 * (S[0] * RS[ia1] + 
	                                                      S[2] * RD[ia1] * PA1) -
	              r2 * SL[ia2] * PPRL * (S[1] * RS[ia1] + S[3] * RD[ia1] * PA1) */
	    COMPLEX_MUL2(ct1, S[0], RS[ia1]);
	    COMPLEX_MUL3(ct2, S[2], RD[ia1], PA1);
	    COMPLEX_ADD2(dt1, ct1, ct2);
	    COMPLEX_SET (ct4, root2*DL[ia2], 0.0);
	    COMPLEX_MUL5(ct5, ct4, PLc[ia2], PPRL, PA2, dt1);
	    COMPLEX_MUL2(ct1, S[1], RS[ia1]);
	    COMPLEX_MUL3(ct2, S[3], RD[ia1], PA1);
	    COMPLEX_ADD2(dt2, ct1, ct2);
	    COMPLEX_SET (ct4, root2*SL[ia2], 0.0);
	    COMPLEX_MUL3(ct6, ct4, PPRL, dt2);
	    COMPLEX_SUB (DFDP, ct5, ct6);
	    /* part2 = -r2 * SL[ia2] * PLc[ia2] * PPRL * PA2 * (S[0] * RS[ia1] + 
	                                                        S[2] * RD[ia1] * PA1) -
	                r2 * DL[ia2] * PPRL * S[1] * RS[ia1] + S[3] * RD[ia1] * PA1) */
	    COMPLEX_SET (ct4, -root2*SL[ia2], 0.0);
	    COMPLEX_MUL5(ct5, ct4, PLc[ia2], PPRL, PA2, dt1);
	    COMPLEX_SET (ct4, root2*DL[ia2], 0.0);
	    COMPLEX_MUL3(ct6, ct4, PPRL, dt2);
	    COMPLEX_SUB (DFDP2, ct5, ct6);
	  }
	}
	break;
      default:
	break;
      }; /* end antenna parameter switch */
      /* end antenna param */
    } else if (paramType==polnParmSou) {   /* Source parameters */
      /* Default partials  */
      COMPLEX_SET (DFDP,  0.0, 0.0);
      COMPLEX_SET (DFDP2, 0.0, 0.0);
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* RL wrt I */
	if (souFit[isou][paramNumber]) {
	  /* part = MC1 + MC4 */
	  COMPLEX_ADD2(DFDP, MC1, MC4);
	}
	break;
      case 1:     /* RL wrt QPol */
	if (souFit[isou][paramNumber]) {
	  /* part = MC2 + MC3 */
	  COMPLEX_ADD2(DFDP, MC2, MC3);
	}
	break;
      case 2:     /* RL wrt UPol */
	if (souFit[isou][paramNumber]) {
	  /* part = i (MC2 - MC3) */
	  COMPLEX_SUB (ct3, MC2, MC3);
	  COMPLEX_SET(ct4, 0.0, 1.0);
	  COMPLEX_MUL2(DFDP, ct4, ct3);
	}
	break;
      case 3:     /* RL wrt Vpol */
	if (souFit[isou][paramNumber]) {
	  /* part = MC1 - MC4 */
	  COMPLEX_SUB (DFDP, MC1, MC4);
	}
	break;  
      default:
	break;
      }; /* end source parameter switch */
      /* end source param */
    } else if (paramType==polnParmPD) {   /* R-L Phase difference */
      COMPLEX_SET (ct1, 0.0,  -1.0);
      COMPLEX_MUL2(DFDP, ct1, VRL);
      COMPLEX_MUL2(DFDP2, ct1, DFDP);
    } /* end R-L phase difference */
    /* Accumulate partials */
    if (paramType!=polnParmUnspec) {
      sumd  += 2.0 * isigma * (residR*DFDP.real + residI*DFDP.imag);
      sumd2 += 2.0 * isigma * (DFDP.real*DFDP.real + DFDP.imag*DFDP.imag +
			       residR*DFDP2.real + residI*DFDP2.imag);
   } /* end set partials */
  } /* end valid data */
  
    /*        LR */
  if (wt[idata*4+3]>0.0) {
    isigma = wt[idata*4+3];
    /* VLR = PPLR * S[0] * LS[ia1] * RSc[ia2] * PA1c +
             PPLR * S[1] * LS[ia1] * RDc[ia2] * PA1c * PA2c +
	     PPLR * S[2] * LD[ia1] * RSc[ia2] +
	     PPLR * S[3] * LD[ia1] * RDc[ia2] * PA2c */
    COMPLEX_MUL4 (MC1, PPLR, LS[ia1], RSc[ia2], PA1c);
    COMPLEX_MUL2 (VLR, S[0], MC1);
    COMPLEX_MUL5 (MC2, PPLR, LS[ia1], RDc[ia2], PA1c,  PA2c);
    COMPLEX_MUL2 (ct1, S[1], MC2);
    COMPLEX_ADD2 (VLR, VLR,  ct1);
    COMPLEX_MUL3 (MC3, PPLR, LD[ia1], RSc[ia2]);
    COMPLEX_MUL2 (ct1, S[2], MC3);
    COMPLEX_ADD2 (VLR, VLR,  ct1);
    COMPLEX_MUL4 (MC4, PPLR, LD[ia1], RDc[ia2],  PA2c);
    COMPLEX_MUL2 (ct1, S[3], MC4);
    COMPLEX_ADD2 (VLR, VLR, ct1);
    residR = VLR.real - data[idata*10+8];
    sum += isigma * residR * residR; sumwt += isigma;
    residI = VLR.imag - data[idata*10+9];
    sum += isigma * residI * residI; sumwt += isigma;
    nXobs++;
    sumXResid += residR * residR + residI * residI;
    /* Derivatives */
    if (paramType==polnParmAnt) {         /* Antenna parameters */
      /* Default partials if only refant on baseline */
      COMPLEX_SET (DFDP,  0.0, 0.0);
      COMPLEX_SET (DFDP2, 0.0, 0.0);
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* LR wrt Or */
	if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* (0,-2 i) * (S[1]*MC2 + S[3]*MC4) */
	    COMPLEX_MUL2(ct1, S[1], MC2);
	    COMPLEX_MUL2(ct2, S[3], MC4);
	    COMPLEX_ADD2(ct3, ct1, ct2);
	    COMPLEX_SET (ct1, 0.0, -2.0);
	    COMPLEX_MUL2(DFDP, ct1, ct3);
	    COMPLEX_MUL2(DFDP2, ct1, DFDP);
	  }
	} 
	/* If ia2==refant */
	if (ia2==args->refAnt) {
	  COMPLEX_SET (ct1, 0.0, -1.0);
	  COMPLEX_MUL2(ct2, ct1, VLR);
	  COMPLEX_ADD2(DFDP, DFDP, ct2);
	  COMPLEX_MUL2(ct2, ct1, DFDP);
	  COMPLEX_ADD2(DFDP2, DFDP2, ct2);
	}
	break;
      case 1:     /* LR wrt Er */
	if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* part = r2 * DR[ia2] * (PPLR * S[0] * LS[ia1] * PA1c + PPLR * S[2] * LD[ia1]) -
	              r2 * SR[ia2] * PRc[ia2] * (PPLR * S[1] * LS[ia1] * PA1c * PA2c + 
	                                         PPLR * S[3] * LD[ia1] * PA2c) */
	    COMPLEX_MUL3(ct1, S[0], LS[ia1], PA1c);
	    COMPLEX_MUL2(ct2, S[2], LD[ia1]);
	    COMPLEX_ADD2(dt1, ct1, ct2);
	    COMPLEX_SET (ct4, root2*DR[ia2], 0.0);
	    COMPLEX_MUL3(ct5, ct4, PPLR, dt1);
	    COMPLEX_MUL3(ct1, S[1], LS[ia1], PA1c);
	    COMPLEX_MUL2(ct2, S[3], LD[ia1]);
	    COMPLEX_ADD2(dt2, ct1, ct2);
	    COMPLEX_SET (ct4, root2*SR[ia2], 0.0);
	    COMPLEX_MUL5(ct6, ct4, PRc[ia2], PPLR, PA2c, dt2);
	    COMPLEX_SUB (DFDP, ct5, ct6);
	    /* part2 = -r2 * SR[ia2] * (PPLR * S[0] * LS[ia1] * PA1c + PPLR * S[2] * LD[ia1]) -
	                r2 * DR[ia2] * PRc[ia2] * (PPLR * S[1] * LS[ia1] * PA1c * PA2c + 
	                                           PPLR * S[3] * LD[ia1] * PA2c) */
	    COMPLEX_SET (ct4, -root2*SR[ia2], 0.0);
	    COMPLEX_MUL3(ct5, ct4, PPLR, dt1);
	    COMPLEX_SET (ct4, root2*DR[ia2], 0.0);
	    COMPLEX_MUL5(ct6, ct4, PRc[ia2], PPLR, PA2c, dt2);
	    COMPLEX_SUB (DFDP2, ct5, ct6);
	  }
	}
	break;
      case 2:     /* LR wrt Ol */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* (0,-2 i) * (S[0]*MC1 + S[1]*MC2) */
	    COMPLEX_MUL2(ct1, S[0], MC1);
	    COMPLEX_MUL2(ct2, S[1], MC2);
	    COMPLEX_ADD2(ct3, ct1, ct2);
	    COMPLEX_SET (ct1, 0.0, -2.0);
	    COMPLEX_MUL2(DFDP, ct1, ct3);
	    COMPLEX_MUL2(DFDP2, ct1, DFDP);
	  }
	  /* If ia1==refant) */
	} else  if (ia1==args->refAnt) {
	  COMPLEX_SET (ct1, 0.0, -1.0);
	  COMPLEX_MUL2(ct2, ct1, VLR);
	  COMPLEX_ADD2(DFDP, DFDP, ct2);
	  COMPLEX_MUL2(ct2, ct1, DFDP);
	  COMPLEX_ADD2(DFDP2, DFDP2, ct2);
	}
	break;
      case 3:     /* LR wrt El */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* part = r2 * DL[ia1] * PL[ia1] * PPLR * PA1c * (S[0] * RSc[ia2 + 
	                                                      S[1] * RDc[ia2] * PA2c) -
	              r2 * SL[ia1] * PPLR * (S[2] * RSc[ia2] + S[3] * RDc[ia2] * PA2c) */
	    COMPLEX_MUL2(ct1, S[0], RSc[ia2]);
	    COMPLEX_MUL3(ct2, S[1], RDc[ia2], PA1c);
	    COMPLEX_ADD2(dt1, ct1, ct2);
	    COMPLEX_SET (ct4, root2*DL[ia1], 0.0);
	    COMPLEX_MUL5(ct5, ct4, PL[ia1], PPLR, PA1c, dt1);
	    COMPLEX_MUL2(ct1, S[2], RSc[ia2]);
	    COMPLEX_MUL3(ct2, S[3], RDc[ia2], PA2c);
	    COMPLEX_ADD2(dt2, ct1, ct2);
	    COMPLEX_SET (ct4, root2*SL[ia1], 0.0);
	    COMPLEX_MUL3(ct6, ct4, PPLR, dt2);
	    COMPLEX_SUB (DFDP, ct5, ct6);
	    /* part2 = -r2 * SL[ia1] * PL[ia1] * PPLR * PA1c * (S[0] * RSc[ia2]  + 
	                                                        S[1] * RDc[ia2] * PA2c) -
	                r2 * DL[ia1] * PPLR * (S[2] * RSc[ia2] + S[3] * RDc[ia2] * PA2c) */
	    COMPLEX_SET (ct4, -root2*SL[ia1], 0.0);
	    COMPLEX_MUL5(ct5, ct4, PL[ia1], PPLR, PA1c, dt1);
	    COMPLEX_SET (ct4, root2*DL[ia1], 0.0);
	    COMPLEX_MUL3(ct6, ct4, PPLR, dt2);
	    COMPLEX_SUB (DFDP2, ct5, ct6);
	  }
	}
	break;
      default:
	break;
      }; /* end antenna parameter switch */
      /* end antenna param */
    } else if (paramType==polnParmSou) {   /* Source parameters */
      /* Default partials  */
      COMPLEX_SET (DFDP,  0.0, 0.0);
      COMPLEX_SET (DFDP2, 0.0, 0.0);
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* LR wrt I */
	if (souFit[isou][paramNumber]) {
	  /* part = MC1 + MC4 */
	  COMPLEX_ADD2(DFDP, MC1, MC4);
	}
	break;
      case 1:     /* LR wrt QPol */
	if (souFit[isou][paramNumber]) {
	  /* part = MC2 + MC3 */
	  COMPLEX_ADD2(DFDP, MC2, MC3);
	}
	break;
      case 2:     /* LR wrt UPol */
	if (souFit[isou][paramNumber]) {
	  /* part = i (MC2 - MC3) */
	  COMPLEX_SUB (ct3, MC2, MC3);
	  COMPLEX_SET(ct4, 0.0, 1.0);
	  COMPLEX_MUL2(DFDP, ct4, ct3);
	}
	break;
      case 3:     /* LR wrt Vpol */
	if (souFit[isou][paramNumber]) {
	  /* part = MC1 - MC4 */
	  COMPLEX_SUB (DFDP, MC1, MC4);
	}
	break;  
      default:
	break;
      }; /* end source parameter switch */
      /* end source param */
    } else if (paramType==polnParmPD) {   /* R-L Phase difference */
      /* No dependence on MORE HERE */
      COMPLEX_SET (ct1, 0.0, 1.0);
      COMPLEX_MUL2(DFDP, ct1, VRL);
      COMPLEX_MUL2(DFDP2, ct1, DFDP);
    } /* end R-L phase difference */
    /* Accumulate partials */
    if (paramType!=polnParmUnspec) {
      sumd  += 2.0*isigma * (residR*DFDP.real + residI*DFDP.imag);
      sumd2 += 2.0*isigma * (DFDP.real*DFDP.real + DFDP.imag*DFDP.imag +
			      residR*DFDP2.real + residI*DFDP2.imag);
      } /* end set partials */
    }  /* end valid data */
  } /* End loop over visibilities */
  
  if (sumwt<=0.0) sumwt = 1.0;  /* Trap no data */
  args->ChiSq       = sum;   /* Save results  */
  args->sumParResid = sumParResid;
  args->sumXResid   = sumXResid;
  args->sumWt       = sumwt;
  args->nPobs       = nPobs;
  args->nXobs       = nXobs;
  if (paramType!=polnParmUnspec) args->sumDeriv  = sumd;
  if (paramType!=polnParmUnspec) args->sumDeriv2 = sumd2;

 /* Indicate completion if threaded */
  if (args->ithread>=0)
    ObitThreadPoolDone (args->thread, (gpointer)&args->ithread);

  return NULL;

} /*  end ThreadPolnFitRLChi2 */

/**
 * Linear polarization version.
 * Threaded Chi**2 evaluator for polarization fitting
 * Evaluates sum [(model-observed) / sigma] and derivatives
 * Parallel hands count 0.3 of cross in sums
 * If selAnt or selSou are set, then only data involving that 
 * source or antenna (or ref ant) is included.
 * \param arg   PolnFitArg pointer with elements:
 * \li lo       First 0-rel datum in Data/Wt arrays
 * \li hi       Last 0-rel datum  in Data/Wt arrays
 * \li selAnt   selected antenna, -1-> all
 * \li selSou   selected source,  -1-> all
 * \li paramType       Parameter type
 * \li polnParmUnspec  Unspecified = don't compute derivatives
 * \li polnParmAnt     Antenna parameter
 * \li polnParmGain    Antenna gains
 * \li polnParmSou     Source parameter
 * \li polnParmPD      Phase difference
 * \li paramNumber  Parameter number,
 *                     Sou: 0=Ipol, 1=Qpol, 2=Upol, 3=VPol
                       Ant: 0= Ori R/X, 1=Elp R/X, 2= Ori L/Y, 1=Elp L/Y,       
 * \li ChiSq        [out] computed Chi2
 * \li nPobs        [out] Number of valid parallel measurements
 * \li nXobs        [out] Number of valid cross pol measurements
 * \li ParRMS       [out] Parallel hand RMS
 * \li sumParResid  [out] Cross hand RMS
 * \li sumXResid    [out] First derivative of Chi2 wrt parameter
 * \li d2Chi2       [out] Second derivative of Chi2 wrt parameter
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadPolnFitXYChi2 (gpointer arg)
{
  PolnFitArg *args = (PolnFitArg*)arg;
  odouble    *antParm   = args->antParm;
  gboolean   **antFit   = args->antFit;
  odouble    *antGain   = args->antGain;
  gboolean   **antGainFit= args->antGainFit;
  odouble    *souParm   = args->souParm;
  gboolean   **souFit   = args->souFit;
  ofloat     *data      = args->inData;
  ofloat     *wt        = args->inWt;
  dcomplex   *CX        = args->RS;
  dcomplex   *SX        = args->RD;
  dcomplex   *CY        = args->LS;
  dcomplex   *SY        = args->LD;
  dcomplex   *CXc       = args->RSc;
  dcomplex   *SXc       = args->RDc;
  dcomplex   *CYc       = args->LSc;
  dcomplex   *SYc       = args->LDc;
  PolnParmType paramType = args->paramType;
  olong paramNumber      = args->paramNumber;
  olong selSou           = args->selSou;
  olong selAnt           = args->selAnt;

  odouble ipol=0.0, qpol=0.0, upol=0.0, vpol=0.0;
  odouble residR=0.0, residI=0.0, isigma=0.0;
  odouble sumParResid, sumXResid;
  ofloat PD, chi1, chi2;
  olong nPobs, nXobs, ia1, ia2, isou, idata, isouLast=-999;
  gboolean isAnt1=FALSE, isAnt2=FALSE;
  size_t i;
  odouble sum=0.0, sumwt=0.0, sumd, sumd2;

  dcomplex  SPA, DPA, SPAc, DPAc, ggPD;
  dcomplex ct1, ct2, ct3, ct4, ct5, dt1, dt2, Jm, Jp;
  dcomplex S[4], VXX, VXY, VYX, VYY, MC1, MC2, MC3, MC4, DFDP, DFDP2;
  dcomplex SM1, SM2, SM3, SM4;

  COMPLEX_SET (S[0], 0.0, 0.0);  /* Initialize poln vector */
  COMPLEX_SET (S[1], 0.0, 0.0);
  COMPLEX_SET (S[2], 0.0, 0.0);
  COMPLEX_SET (S[3], 0.0, 0.0);
  COMPLEX_SET (MC1, 0.0, 0.0);  /* Other stuff */
  COMPLEX_SET (MC2, 0.0, 0.0);
  COMPLEX_SET (MC3, 0.0, 0.0);
  COMPLEX_SET (MC4, 0.0, 0.0);
  COMPLEX_SET (VXX, 0.0, 0.0);
  COMPLEX_SET (VYY, 0.0, 0.0);
  COMPLEX_SET (VYX, 0.0, 0.0);
  COMPLEX_SET (VXY, 0.0, 0.0);
  COMPLEX_SET (DFDP,  0.0, 0.0);
  COMPLEX_SET (DFDP2, 0.0, 0.0);
  COMPLEX_SET (dt1, 0.0, 0.0);
  COMPLEX_SET (dt2, 0.0, 0.0);
  COMPLEX_SET (Jm,  0.0,-1.0);
  COMPLEX_SET (Jp,  0.0, 1.0);
 
  /* RMS sums and counts */
  sumParResid = sumXResid = 0.0;
  sumd = sumd2 = 0.0;
  nPobs = nXobs = 0;
  /* X-Y phase difference  at reference antenna */
  if (args->doFitRL) {
    PD = args->PD;
  } else PD = 0.0;

  /* Injest model, factorize into antenna components - 
     data in order Orientation R/X, Elipticity R/X, Orientation L/Y, Elipticity L/Y */
  /* Elipticity, Orientation terms */
  for (i=0; i<args->nant; i++) {
    COMPLEX_EXP (ct1, -antParm[i*4+0]);
    COMPLEX_SET (ct2, cos(G_PI*0.25+antParm[i*4+1]), 0.0);
    COMPLEX_MUL3 (CX[i], Jp, ct1, ct2);
    COMPLEX_EXP (ct1, antParm[i*4+0]);
    COMPLEX_SET (ct2, sin(G_PI*0.25+antParm[i*4+1]), 0.0);
    COMPLEX_MUL2 (SX[i], ct1, ct2);
    COMPLEX_EXP (ct1, antParm[i*4+2]);
    COMPLEX_SET (ct2, cos(G_PI*0.25-antParm[i*4+3]), 0.0);
    COMPLEX_MUL2 (CY[i], ct1, ct2);
    COMPLEX_EXP (ct1, -antParm[i*4+2]);
    COMPLEX_SET (ct2, sin(G_PI*0.25-antParm[i*4+3]), 0.0);
    COMPLEX_MUL3 (SY[i], Jp, ct1, ct2);
    COMPLEX_CONJUGATE (CXc[i], CX[i]);
    COMPLEX_CONJUGATE (SXc[i], SX[i]);
    COMPLEX_CONJUGATE (CYc[i], CY[i]);
    COMPLEX_CONJUGATE (SYc[i], SY[i]);
  }

  /* Loop over data */
  i = 0;
  for (idata=args->lo; idata<args->hi; idata++) {
    /* Parallactic angle terms */
    chi1  = data[idata*10+0];   /* parallactic angle ant 1 */
    chi2  = data[idata*10+1];   /* parallactic angle ant 2 */
    COMPLEX_EXP (SPA,chi1+chi2);
    COMPLEX_EXP (DPA,chi1-chi2);
    COMPLEX_CONJUGATE (SPAc, SPA);
    COMPLEX_CONJUGATE (DPAc, DPA);

    isou  = MAX (0, args->souNo[idata]);    /* Source number */
    /* Selected source? */
    if ((selSou>=0) && (selSou!=isou)) continue;

    /* New source? get parameters */
    if (isou!=isouLast) {
      isouLast = isou;
      /* Source parameters */
      ipol = souParm[isou*4+0];
      /* Fitting or fixed? */
      if (args->souFit[isou][1]) 
	qpol = souParm[isou*4+1];
      else
	qpol = args->PPol[isou]*ipol*cos(args->RLPhase[isou]);
      if (args->souFit[isou][2]) 
	upol = souParm[isou*4+2];
      else
	upol = args->PPol[isou]*ipol*sin(args->RLPhase[isou]);
      vpol = souParm[isou*4+3];
      /* Complex Stokes array */
      COMPLEX_SET (S[0], ipol+vpol, 0.0);
      COMPLEX_SET (S[1], qpol,  upol);
      COMPLEX_SET (S[2], qpol, -upol);
      COMPLEX_SET (S[3], ipol-vpol, 0.0);
    }

    /* Antenna parameters (0 ref) */
    ia1    = args->antNo[idata*2+0];
    ia2    = args->antNo[idata*2+1]; 
    /* Selected source? */
    if ((selAnt>=0) && (selAnt!=ia1) && (selAnt!=ia2)) continue;
    /* Which antenna is the selected one in the baseline? */
    if (selAnt==ia1) {isAnt1 = TRUE; isAnt2 = FALSE;}
    else if (selAnt==ia2) {isAnt2 = TRUE; isAnt1 = FALSE;}
    
    /* Calculate residals -  XX */
  if (wt[idata*4]>0.0) {
    isigma = wt[idata*4];
    
    /* VXX = {S[0] * CX[ia1] * CXc[ia2] * DPAc  +
              S[1] * CX[ia1] * SXc[ia2] * SPAc  +
	      S[2] * SX[ia1] * CXc[ia2] * SPA   + 
	      S[3] * SX[ia1] * SXc[ia2] * DPA} * g1X * g2X ;
    */
    COMPLEX_MUL3 (MC1, CX[ia1], CXc[ia2], DPAc);
    COMPLEX_MUL3 (MC2, CX[ia1], SXc[ia2], SPAc);
    COMPLEX_MUL3 (MC3, SX[ia1], CXc[ia2], SPA);
    COMPLEX_MUL3 (MC4, SX[ia1], SXc[ia2], DPA);
    COMPLEX_MUL2 (SM1, S[0], MC1);
    COMPLEX_MUL2 (SM2, S[1], MC2);
    COMPLEX_MUL2 (SM3, S[2], MC3);
    COMPLEX_MUL2 (SM4, S[3], MC4);
    COMPLEX_ADD4 (VXX, SM1, SM2, SM3, SM4);
    COMPLEX_SET (ggPD,  antGain[ia1*2+0]*antGain[ia2*2+0], 0);
    COMPLEX_MUL2 (VXX, VXX, ggPD);
    residR = VXX.real - data[idata*10+2];
    sum += isigma * residR * residR; sumwt += isigma;
    residI = VXX.imag - data[idata*10+3];
    sum += isigma * residI * residI; sumwt += isigma;
    /* DEBUG  
       if ((paramType==polnParmAnt) && (paramNumber==1) && (selAnt==2) && (ia1==2) && (ia2==3)) {
       fprintf (stdout, "vis %4d %d %d XX chi2 %8.3f M %8.3f %8.3f O %8.3f %8.3f R %8.3f %8.3f \n",
       idata,ia1, ia2, isigma*(residR*residR+residI*residI),
       VXX.real, VXX.imag, data[idata*10+2], data[idata*10+3], residR, residI);
       if (idata==0) {
       fprintf (stdout, "ant %d parm %f %f %f %f PD %f\n", 
       ia1, antParm[ia1*4+0],  antParm[ia1*4+1], antParm[ia1*4+2], antParm[ia1*4+3], PD);
       fprintf (stdout, "ant %d parm %f %f %f %f PA %f %f\n", 
       ia2, antParm[ia2*4+0],  antParm[ia2*4+1], antParm[ia2*4+2], antParm[ia2*4+3], chi1, chi2);
       }
       }   End Debug  */
    nPobs++;
    sumParResid += residR * residR + residI * residI;
    /* Derivatives */
    /* Default partials */
    COMPLEX_SET (DFDP,  0.0, 0.0);
    COMPLEX_SET (DFDP2, 0.0, 0.0);
    if (paramType==polnParmAnt) {         /* Antenna parameters */
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* XX wrt Ox */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* {(0, -1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
	        (0,  1) * S[2] * MC3 + (0,  1) * S[3] * MC4}  * gX[ia1] * gX[ia2]  */
	    COMPLEX_ADD2 (ct1, SM1, SM2);
	    COMPLEX_MUL2 (ct2, Jm, ct1);
	    COMPLEX_ADD2 (ct1, SM3, SM4);
	    COMPLEX_MUL2 (ct3, Jp, ct1);
	    COMPLEX_ADD2 (ct2, ct2, ct3);
	    COMPLEX_MUL2 (DFDP, ct2, ggPD);
	    /* -VXX */
	    COMPLEX_NEGATE(DFDP2, VXX);
	  }
	} else if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
	        (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4}  * gX[ia1] * gX[ia2]   */
	    COMPLEX_ADD2 (ct1, SM1, SM3);
	    COMPLEX_MUL2 (ct2, Jp, ct1);
	    COMPLEX_ADD2 (ct1, SM2, SM4);
	    COMPLEX_MUL2 (ct3, Jm, ct1);
	    COMPLEX_ADD2 (ct2, ct2, ct3);
	    COMPLEX_MUL2 (DFDP, ct2, ggPD);
	    /*  -VXX  */
	    COMPLEX_NEGATE(DFDP2, VXX);
	  }
	} 
	break;
      case 1:     /* XX wrt Ex */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* part = (0, -1) * {S[0] * CXc[ia2] * DPAc * SXc[ia1] +
	                         S[1] * SXc[ia2] * SPAc * SXc[ia1] +
				 S[2] * CXc[ia2] * SPA  * CXc[ia1] +	       
				 S[3] * SXc[ia2] * DPA  * CXc[ia1]}  * gX[ia1] * gX[ia2]  */
	    COMPLEX_MUL4(ct1, S[0], CXc[ia2], DPAc, SXc[ia1]);
	    COMPLEX_MUL4(ct2, S[1], SXc[ia2], SPAc, SXc[ia1]);
	    COMPLEX_MUL4(ct3, S[2], CXc[ia2], SPA,  SXc[ia1]);
	    COMPLEX_MUL4(ct4, S[3], SXc[ia2], DPA,  CXc[ia1]);
	    COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
	    COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
	    /* part2 = -VXX */
	    COMPLEX_NEGATE(DFDP2, VXX);
	  }
	} else if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* part = (0, -1) * {S[0] * CX[ia1] * DPAc * SX[ia2] +
                                 S[1] * CX[ia1] * SPAc * CX[ia2] +
				 S[2] * SX[ia1] * SPA  * SX[ia2] +
				 S[3] * SX[ia1] * DPA  * CX[ia2]} * gX[ia1] * gX[ia2]  */
	    COMPLEX_MUL4(ct1, S[0], CX[ia1], DPAc, SX[ia2]);
	    COMPLEX_MUL4(ct2, S[1], CX[ia1], SPAc, CX[ia2]);
	    COMPLEX_MUL4(ct3, S[2], SX[ia1], SPA,  SX[ia2]);
	    COMPLEX_MUL4(ct4, S[3], SX[ia1], DPA,  CX[ia2]);
	    COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
	    COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
	    /* part2 = -VXX */
	    COMPLEX_NEGATE(DFDP2, VXX);
	  }
	}
	break;
      case 2:     /* XX wrt Oy - nope */
	break;
      case 3:     /* XX wrt Ey - nope */
	break;
      default:
	break;
	}; /* end antenna parameter switch */
	/* end antenna param */
      } else if (paramType==polnParmSou) {   /* Source parameters */
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* XX wrt I */
	if (souFit[isou][paramNumber]) {
	  /* part =   (MC1 + MC4) * gX[ia1] * gX[ia2] */
	  COMPLEX_ADD2(ct1, MC1, MC4);
	  COMPLEX_MUL2(DFDP, ct1, ggPD);
	}
	break;
      case 1:     /* XX wrt QPol */
	if (souFit[isou][paramNumber]) {
	  /* part = (MC2 + MC3) * gX[ia1] * gX[ia2] */
	  COMPLEX_ADD2(ct1, MC2, MC3);
	  COMPLEX_MUL2(DFDP, ct1, ggPD);
	}
	break;
      case 2:     /* XX wrt UPol */
	if (souFit[isou][paramNumber]) {
	  /* i (MC2 - MC3) * gX[ia1] * gX[ia2] */
	  COMPLEX_SUB(ct1, MC2, MC3);
	  COMPLEX_MUL3(DFDP, Jp, ct1, ggPD);
	}
	break;
      case 3:     /* XX wrt Vpol */
	if (souFit[isou][paramNumber]) {
	  /* (MC1 - MC4) * gX[ia1] * gX[ia2] */
	  COMPLEX_SUB(ct1, MC1, MC4);
	  COMPLEX_MUL2(DFDP, ct1, ggPD);
	}
	break;  
      default:
	break;
      }; /* end source parameter switch */
      /* end source param */
    } else if (paramType==polnParmGain) {   /* Antenna Gains */
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* XX wrt gX */
	if (isAnt1 && (antGainFit[ia1][paramNumber])) {
	  /* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gX[ia2] */
	  COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct2, antGain[ia2*2+0], 0);
	  COMPLEX_MUL2 (DFDP, ct1, ct2);
	  /* part2 = 0 */
	} else if (isAnt2 && (antGainFit[ia2][paramNumber])) {
	  /* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gX[ia1] */
	  COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct2, antGain[ia1*2+0], 0);
	  COMPLEX_MUL2 (DFDP, ct1, ct2);
	  /* part2 = 0 */
	}
	break;
      case 1:     /* XX wrt gY - nope */
	break;
       default:
	break;
      }  /* end gain switch */
    } else if (paramType==polnParmPD) {   /* X-Y phase difference */
    } /* end parameter types */
    /* Accumulate partials */
    if (paramType!=polnParmUnspec) {
      sumd  += 2.0 * isigma * (residR*DFDP.real + residI*DFDP.imag); 
      sumd2 += 2.0 * isigma * (DFDP.real*DFDP.real + DFDP.imag*DFDP.imag +
			       residR*DFDP2.real + residI*DFDP2.imag);
    } /* end set partials */
  } /* end valid data */
  
  /*      YY */
  if (wt[idata*4+1]>0.0) {
    isigma = wt[idata*4+1];
    /* VYY = {S[0] * SY[ia1] * SYc[ia2] * DPAc +       
              S[1] * SY[ia1] * CYc[ia2] * SPAc +
	      S[2] * CY[ia1] * SYc[ia2] * SPA  + 
	      S[3] * CY[ia1] * CYc[ia2] * DPA} * g1Y * g2Y ;
    */
    COMPLEX_MUL3 (MC1, SY[ia1], SYc[ia2], DPAc);
    COMPLEX_MUL3 (MC2, SY[ia1], CYc[ia2], SPAc);
    COMPLEX_MUL3 (MC3, CY[ia1], SYc[ia2], SPA);
    COMPLEX_MUL3 (MC4, CY[ia1], CYc[ia2], DPA);
    COMPLEX_MUL2 (SM1, S[0], MC1);
    COMPLEX_MUL2 (SM2, S[1], MC2);
    COMPLEX_MUL2 (SM3, S[2], MC3);
    COMPLEX_MUL2 (SM4, S[3], MC4);
    COMPLEX_ADD4 (VYY, SM1, SM2, SM3, SM4);
    COMPLEX_SET (ggPD,  antGain[ia1*2+1]*antGain[ia2*2+1], 0);
    COMPLEX_MUL2 (VYY, VYY, ggPD);
    residR = VYY.real - data[idata*10+4];
    sum += isigma * residR * residR; sumwt += isigma; 
    residI = VYY.imag - data[idata*10+5];
    sum += isigma * residI * residI; sumwt += isigma; 
    /* DEBUG 
       if ((paramType==polnParmAnt) && (paramNumber==1) && (selAnt==2) && (ia1==2) && (ia2==3)) {
       fprintf (stdout, "vis %4d %d %d YY chi2 %8.3f M %8.3f %8.3f O %8.3f %8.3f R %8.3f %8.3f \n",
       idata,ia1, ia2, isigma*(residR*residR+residI*residI),
       VYY.real, VYY.imag, data[idata*10+4], data[idata*10+5], residR, residI);
       }  End Debug */
    nPobs++;
    sumParResid += residR * residR + residI * residI;
    /* Derivatives */
    /* Default partials  */
    COMPLEX_SET (DFDP,  0.0, 0.0);
    COMPLEX_SET (DFDP2, 0.0, 0.0);
    if (paramType==polnParmAnt) {         /* Antenna parameters */
      /* Default partials */
      COMPLEX_SET (DFDP,  0.0, 0.0);
      COMPLEX_SET (DFDP2, 0.0, 0.0);
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* YY wrt Ox - nope */
	break;
      case 1:     /* YY wrt Ex - nope*/
	break;
      case 2:     /* YY wrt Oy */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* part = {(0, -1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
	               (0,  1) * S[2] * MC3 + (0,  1) * S[3] * MC4} * gY[ia1] * gY[ia2] */
	    COMPLEX_ADD2 (ct1, SM1, SM2);
	    COMPLEX_MUL2 (ct2, Jm, ct1);
	    COMPLEX_ADD2 (ct1, SM3, SM4);
	    COMPLEX_MUL2 (ct3, Jp, ct1);
	    COMPLEX_ADD2 (ct2, ct2, ct3);
	    COMPLEX_MUL2 (DFDP, ct2, ggPD);
	    /* part2 = -VYY */
	    COMPLEX_NEGATE(DFDP2, VYY);
	  }
	} else if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* part = {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
                       (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4} * gY[ia1] * gY[ia2]  */
	    COMPLEX_ADD2 (ct1, SM1, SM3);
	    COMPLEX_MUL2 (ct2, Jp, ct1);
	    COMPLEX_ADD2 (ct1, SM2, SM4);
	    COMPLEX_MUL2 (ct3, Jm, ct1);
	    COMPLEX_ADD2 (ct2, ct2, ct3);
	    COMPLEX_MUL2 (DFDP,  ct2, ggPD);
	    /* part2 = -VYY */
	    COMPLEX_NEGATE(DFDP2, VYY);
	  }
	}
	break;
      case 3:     /* YY wrt Ey */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* part = (0, -1) * {S[0] * SYc[ia2] * DPAc * CYc[ia1] +
                                 S[1] * CYc[ia2] * SPAc * CYc[ia1] +
				 S[2] * SYc[ia2] * SPA  * SYc[ia1] +	       
				 S[3] * CYc[ia2] * DPA  * SYc[ia1]}  * gX[ia1] * gX[ia2]  */
	    COMPLEX_MUL4(ct1, S[0], SYc[ia2], DPAc, CYc[ia1]);
	    COMPLEX_MUL4(ct2, S[1], CYc[ia2], SPAc, CYc[ia1]);
	    COMPLEX_MUL4(ct3, S[2], SYc[ia2], SPA,  SYc[ia1]);
	    COMPLEX_MUL4(ct4, S[3], CYc[ia2], DPA,  SYc[ia1]);
	    COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
	    COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
	    /* part2 = -VYY */
	    COMPLEX_NEGATE(DFDP2, VYY);
	  }
	} else if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* part = (0, -1) * {S[0] * SY[ia1] * DPAc * CY[ia2] + 
                                 S[1] * SY[ia1] * SPAc * SY[ia2] +
				 S[2] * CY[ia1] * SPA  * CY[ia2]  + 
				 S[3] * CY[ia1] * DPA  * SY[ia2]} * gY[ia1] * gY[ia2] */
	    COMPLEX_MUL4(ct1, S[0], SY[ia1], DPAc, CY[ia2]);
	    COMPLEX_MUL4(ct2, S[1], SY[ia1], SPAc, SY[ia2]);
	    COMPLEX_MUL4(ct3, S[2], CY[ia1], SPA,  CY[ia2]);
	    COMPLEX_MUL4(ct4, S[3], CY[ia1], DPA,  SY[ia2]);
	    COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
	    COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
	    /* part2 = -VYY */
	    COMPLEX_NEGATE(DFDP2, VYY);
	  }
	}
	break;
      default:
	break;
      }; /* end antenna parameter switch */
      /* end antenna param */
    } else if (paramType==polnParmSou) {   /* Source parameters */
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* YY wrt IPol */
	if (souFit[isou][paramNumber]) {
	  /* part =   (MC1 + MC4) * gY[ia1] * gY[ia2] */
	  COMPLEX_ADD2(ct1, MC1, MC4);
	  COMPLEX_MUL2(DFDP, ct1, ggPD);
	}
	break;
      case 1:     /* YY wrt QPol */
	if (souFit[isou][paramNumber]) {
	  /* part = (MC2 + MC3) * gY[ia1] * gY[ia2] */
	  COMPLEX_ADD2(ct1, MC2, MC3);
	  COMPLEX_MUL2(DFDP, ct1, ggPD);
	}
	break;
      case 2:     /* YY wrt UPol */
	if (souFit[isou][paramNumber]) {
	  /* i (MC2 - MC3) * gY[ia1] * gY[ia2] */
	  COMPLEX_SUB(ct1, MC2, MC3);
	  COMPLEX_MUL3(DFDP, Jp, ct1, ggPD);
	}
	break;
      case 3:     /* YY wrt Vpol */
	if (souFit[isou][paramNumber]) {
	  /* (MC1 - MC4) * gY[ia1] * gY[ia2] */
	  COMPLEX_SUB (ct1, MC1, MC4);
	  COMPLEX_MUL2(DFDP, ct1, ggPD);
	}
	break;  
      default:
	break;
      }; /* end source parameter switch */
      /* end source param */
    } else if (paramType==polnParmGain) {   /* Antenna Gains */
     switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* YY wrt gX - nope */
	break;
      case 1:     /* YY wrt gY */
	if (isAnt1 && (antGainFit[ia1][paramNumber])) {
	  /* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gY[ia2] */
	  COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct2, antGain[ia2*2+1], 0);
	  COMPLEX_MUL2 (DFDP, ct1, ct2);
	  /* part2 = 0 */
	} else 	if (isAnt2 && (antGainFit[ia2][paramNumber])) {
	  /* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gY[ia1] */
	  COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct2, antGain[ia1*2+1], 0);
	  COMPLEX_MUL2 (DFDP, ct1, ct2);
	  /* part2 = 0 */
	}
	break;
       default:
	break;
      } /* end gain switch */
    } else if (paramType==polnParmPD) {   /* X-Y phase difference */
      /* Nope */
    } /* end parameter types */

    /* Accumulate partials */
    if (paramType!=polnParmUnspec) {
      sumd  += 2.0 * isigma * (residR*DFDP.real + residI*DFDP.imag); 
      sumd2 += 2.0 * isigma * (DFDP.real*DFDP.real + DFDP.imag*DFDP.imag +
			       residR*DFDP2.real + residI*DFDP2.imag);
    } /* end set partials */
  } /* end valid data */
  
    /* 	    XY */
  if (wt[idata*4+2]>0.0) {
    isigma = wt[idata*4+2];
    /* VXY = {S[0] * CX[ia1] * SYc[ia2] * DPAc +       
              S[1] * CX[ia1] * CYc[ia2] * SPAc +
	      S[2] * SX[ia1] * SYc[ia2] * SPA  + 
	      S[3] * SX[ia1] * CYc[ia2] * DPA}} * g1X * g2Y * exp(i PD);
    */
    COMPLEX_MUL3 (MC1, CX[ia1], SYc[ia2], DPAc);
    COMPLEX_MUL3 (MC2, CX[ia1], CYc[ia2], SPAc);
    COMPLEX_MUL3 (MC3, SX[ia1], SYc[ia2], SPA);
    COMPLEX_MUL3 (MC4, SX[ia1], CYc[ia2], DPA);
    COMPLEX_MUL2 (SM1, S[0], MC1);
    COMPLEX_MUL2 (SM2, S[1], MC2);
    COMPLEX_MUL2 (SM3, S[2], MC3);
    COMPLEX_MUL2 (SM4, S[3], MC4);
    COMPLEX_ADD4 (VXY, SM1, SM2, SM3, SM4);
    COMPLEX_SET (ct1,  antGain[ia1*2+0]*antGain[ia2*2+1], 0);
    COMPLEX_EXP (ct2, PD);
    COMPLEX_MUL2 (ggPD, ct1, ct2);
    COMPLEX_MUL2 (VXY, VXY, ggPD);
    residR = VXY.real - data[idata*10+6];
    sum += isigma * residR * residR; sumwt += isigma;
    residI = VXY.imag - data[idata*10+7];
    sum += isigma * residI * residI; sumwt += isigma;
    /* DEBUG 
       if ((paramType==polnParmAnt) && (paramNumber==1) && (selAnt==2) && (ia1==2) && (ia2==3)) {
       fprintf (stdout, "vis %4d %d %d XY chi2 %8.3f M %8.3f %8.3f O %8.3f %8.3f R %8.3f %8.3f SM %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ggPD %8.3f %8.3f \n",
       idata,ia1, ia2, isigma*(residR*residR+residI*residI),
       VXY.real, VXY.imag, data[idata*10+6], data[idata*10+7], residR, residI,
       SM1.real, SM1.imag, SM2.real, SM2.imag, SM3.real, SM3.imag, SM4.real, SM4.imag, ggPD.real, ggPD.imag);
       }  End Debug */
    nXobs++;
    sumXResid += residR * residR + residI * residI;
    /* Derivatives */
    /* Default partials  */
    COMPLEX_SET (DFDP,  0.0, 0.0);
    COMPLEX_SET (DFDP2, 0.0, 0.0);
     if (paramType==polnParmAnt) {         /* Antenna parameters */
     switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* XY wrt Ox */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* part = {(0, -1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
	               (0,  1) * S[2] * MC3 + (0,  1) * S[3] * MC4} *
		                        gX[ia1] * gY[ia2] * exp(i PD) */
	    COMPLEX_ADD2 (ct1, SM1, SM2);
	    COMPLEX_MUL2 (ct2, Jm, ct1);
	    COMPLEX_ADD2 (ct1, SM3, SM4);
	    COMPLEX_MUL2 (ct3, Jp, ct1);
	    COMPLEX_ADD2 (ct2, ct2, ct3);
	    COMPLEX_MUL2 (DFDP, ct2, ggPD);
	    /* part2 = -VXY */
	    COMPLEX_NEGATE(DFDP2, VXY);
	  }
	} 
	break;
      case 1:     /* XY wrt Ex */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* part =  (0, -1) * {S[0] * SYc[ia2] * DPAc * SXc[ia1]  + 
                                  S[1] * CYc[ia2] * SPAc * SXc[ia1]  +
				  S[2] * SYc[ia2] * SPA  * CXc[ia1]  + 
				  S[3] * CYc[ia2] * DPA  * CXc[ia1]  } * 
				         gX[ia1] * gY[ia2] * exp(i PD) */
	    COMPLEX_MUL4(ct1, S[0], SYc[ia2], DPAc, SXc[ia1]);
	    COMPLEX_MUL4(ct2, S[1], CYc[ia2], SPAc, SXc[ia1]);
	    COMPLEX_MUL4(ct3, S[2], SYc[ia2], SPA,  CXc[ia1]);
	    COMPLEX_MUL4(ct4, S[3], CYc[ia2], DPA,  CXc[ia1]);
	    COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
	    COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
	    /*   part2 = -VXY  */
	    COMPLEX_NEGATE(DFDP2, VXY);
	  }
	}
	break;
      case 2:     /* XY wrt Oy */
	if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* part = {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
	               (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4} * 
                             gX[ia1] * gY[ia2] * exp(i PD) */
	    COMPLEX_ADD2 (ct1, SM1, SM3);
	    COMPLEX_MUL2 (ct2, Jp, ct1);
	    COMPLEX_ADD2 (ct1, SM2, SM4);
	    COMPLEX_MUL2 (ct3, Jm, ct1);
	    COMPLEX_ADD2 (ct2, ct2, ct3);
	    COMPLEX_MUL2 (DFDP,  ct2, ggPD);
	    /* part2 = -VXY */
	    COMPLEX_NEGATE(DFDP2, VXY);
	  }
	} 
	break;
      case 3:     /* XY wrt Ey */
	if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* part = (0, -1) *{S[0] * CX[ia1] * DPAc * CY[ia2] +       
	                        S[1] * CX[ia1] * SPAc * SY[ia2] +
				S[2] * SX[ia1] * SPA  * CY[ia2] + 
				S[3] * SX[ia1] * DPA  * SY[ia2]} * 
			   	       gX[ia1] * gY[ia2] * exp(i PD) */
	    COMPLEX_MUL4(ct1, S[0], CX[ia1], DPAc, CY[ia2]);
	    COMPLEX_MUL4(ct2, S[1], CX[ia1], SPAc, SY[ia2]);
	    COMPLEX_MUL4(ct3, S[2], SX[ia1], SPA,  CY[ia2]);
	    COMPLEX_MUL4(ct4, S[3], SX[ia1], DPA,  SY[ia2]);
	    COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
	    COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
	    /* part2 = -VXY */
	    COMPLEX_NEGATE(DFDP2, VXY);
	  }
	}
	break;
      default:
	break;
      }; /* end antenna parameter switch */
      /* end antenna param */
    } else if (paramType==polnParmSou) {   /* Source parameters */
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* XY wrt I */
	if (souFit[isou][paramNumber]) {
	  /* part =   (MC1 + MC4) * gX[ia1] * gY[ia2] * exp(i PD) */
	  COMPLEX_ADD2(ct1, MC1, MC4);
	  COMPLEX_MUL2(DFDP, ct1, ggPD);
	}
	break;
      case 1:     /* XY wrt QPol */
	if (souFit[isou][paramNumber]) {
	  /* part = (MC2 + MC3) * gX[ia1] * gY[ia2] * exp(i PD) */
	  COMPLEX_ADD2(ct1, MC2, MC3);
	  COMPLEX_MUL2(DFDP, ct1, ggPD);
	}
	break;
      case 2:     /* XY wrt UPol */
	if (souFit[isou][paramNumber]) {
	  /* i (MC2 - MC3) * gX[ia1] * gY[ia2] * exp(i PD) */
	  COMPLEX_SUB(ct1, MC2, MC3);
	  COMPLEX_MUL3(DFDP, Jp, ct1, ggPD);
	}
	break;
      case 3:     /* XY wrt Vpol */
	if (souFit[isou][paramNumber]) {
	  /* (MC1 - MC4) * gX[ia1] * gY[ia2] * exp(i PD) */
	  COMPLEX_SUB (ct1, MC1, MC4);
	  COMPLEX_MUL2(DFDP, ct1, ggPD);
	}
	break;  
      default:
	break;
      }; /* end source parameter switch */
      /* end source param */
    } else if (paramType==polnParmGain) {   /* Antenna Gains */
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* XY wrt gX */
	if (isAnt1 && (antGainFit[ia1][paramNumber])) {
	  /* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) 
                           * gY[ia2]  * exp(i PD) */
	  COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct2, antGain[ia2*2+1], 0);
	  COMPLEX_EXP (ct3, PD);
	  COMPLEX_MUL3 (DFDP, ct1, ct2, ct3);
	  /* part2 = 0 */
	}
	break;
      case 1:     /* XY wrt gY */
	if (isAnt2 && (antGainFit[ia2][paramNumber])) {
	  /* part = {S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) 
	                    * gX[ia1] * exp(i PD) */
	  COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct2, antGain[ia1*2+0], 0);
	  COMPLEX_EXP (ct3, PD);
	  COMPLEX_MUL3 (DFDP, ct1, ct2, ct3);
	  /* part2 = 0 */
	}
	break;
       default:
	break;
      } /* end gain switch */
  } else if (paramType==polnParmPD) {   /* X-Y phase difference */
      /* part = (0,  1) * VXY */   
      COMPLEX_MUL2 (DFDP, Jp, VXY);
      /* part2 = -VXY */
      COMPLEX_NEGATE(DFDP2, VXY);
  } /* end parameter types */

   /* Accumulate partials */
    if (paramType!=polnParmUnspec) {
      sumd  += 2.0 * isigma * (residR*DFDP.real + residI*DFDP.imag);
      sumd2 += 2.0 * isigma * (DFDP.real*DFDP.real + DFDP.imag*DFDP.imag +
			       residR*DFDP2.real + residI*DFDP2.imag);
   } /* end set partials */
  } /* end valid data */
  
    /*        YX */
  if (wt[idata*4+3]>0.0) {
    isigma = wt[idata*4+3];
    /* VYX = {S[0] * SY[ia1] * CXc[ia2] * DPAc +       
              S[1] * SY[ia1] * SXc[ia2] * SPAc +
	      S[2] * CY[ia1] * CXc[ia2] * SPA  + 
	      S[3] * CY[ia1] * SXc[ia2] * DPA} * g1Y * g2X * exp(-i PD)
    */
    COMPLEX_MUL3 (MC1, SY[ia1], CXc[ia2], DPAc);
    COMPLEX_MUL3 (MC2, SY[ia1], SXc[ia2], SPAc);
    COMPLEX_MUL3 (MC3, CY[ia1], CXc[ia2], SPA);
    COMPLEX_MUL3 (MC4, CY[ia1], SXc[ia2], DPA);
    COMPLEX_MUL2 (SM1, S[0], MC1);
    COMPLEX_MUL2 (SM2, S[1], MC2);
    COMPLEX_MUL2 (SM3, S[2], MC3);
    COMPLEX_MUL2 (SM4, S[3], MC4);
    COMPLEX_ADD4 (VYX, SM1, SM2, SM3, SM4);
    COMPLEX_SET (ct1,  antGain[ia1*2+1]*antGain[ia2*2+0], 0);
    COMPLEX_EXP (ct2, -PD);
    COMPLEX_MUL2 (ggPD, ct1, ct2);
    COMPLEX_MUL2 (VYX, VYX, ggPD);
    residR = VYX.real - data[idata*10+8];
    sum += isigma * residR * residR; sumwt += isigma;
    residI = VYX.imag - data[idata*10+9];
    sum += isigma * residI * residI; sumwt += isigma;
    /* DEBUG
       if ((paramType==polnParmAnt) && (paramNumber==1) && (selAnt==2) && (ia1==2) && (ia2==3)) {
       fprintf (stdout, "vis %4d %d %d YX chi2 %8.3f M %8.3f %8.3f O %8.3f %8.3f R %8.3f %8.3f SM %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ggPD %8.3f %8.3f \n",
       idata,ia1, ia2, isigma*(residR*residR+residI*residI),
       VYX.real, VYX.imag, data[idata*10+8], data[idata*10+9], residR, residI,
       SM1.real, SM1.imag, SM2.real, SM2.imag, SM3.real, SM3.imag, SM4.real, SM4.imag, ggPD.real, ggPD.imag);
       }   End Debug */
    nXobs++;
    sumXResid += residR * residR + residI * residI;
    /* Derivatives */
    /* Default partials  */
    COMPLEX_SET (DFDP,  0.0, 0.0);
    COMPLEX_SET (DFDP2, 0.0, 0.0);      
    if (paramType==polnParmAnt) {         /* Antenna parameters */
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* YX wrt Ox */
	if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* part = {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
                       (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4} * 
		       gY[ia1] * gX[ia2]  * exp(-i PD) */
	    COMPLEX_ADD2 (ct1, SM1, SM3);
	    COMPLEX_MUL2 (ct2, Jp, ct1);
	    COMPLEX_ADD2 (ct1, SM2, SM4);
	    COMPLEX_MUL2 (ct3, Jm, ct1);
	    COMPLEX_ADD2 (ct2, ct2, ct3);
	    COMPLEX_MUL2 (DFDP,  ct2, ggPD);
	    /*   part2 = -VYX  */
	    COMPLEX_NEGATE(DFDP2, VYX);
	  }
	} 
	break;
      case 1:     /* YX wrt Ex */
	if (isAnt2) {
	  if (antFit[ia2][paramNumber]) {
	    /* part = (0, -1) * {S[0] * SY[ia1] * DPAc * SX[ia2] + 
	                         S[1] * SY[ia1] * SPAc * CX[ia2] +
				 S[2] * CY[ia1] * SPA  * SX[ia2] + 
				 S[3] * CY[ia1] * DPA  * CX[ia2]} * 
				        gY[ia1] * gX[ia2] * exp(-i PD)*/
	    COMPLEX_MUL4(ct1, S[0], SY[ia1], DPAc, SX[ia2]);
	    COMPLEX_MUL4(ct2, S[1], SY[ia1], SPAc, CX[ia2]);
	    COMPLEX_MUL4(ct3, S[2], CY[ia1], SPA,  SX[ia2]);
	    COMPLEX_MUL4(ct4, S[3], CY[ia1], DPA,  CX[ia2]);
	    COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
	    COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
	    /*   part2 = -VYX  */
	    COMPLEX_NEGATE(DFDP2, VYX);
	  }
	}
	break;
      case 2:     /* YX wrt Oy */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* part = {(0, -1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
	               (0,  1) * S[2] * MC3 + (0,  1) * S[3] * MC4} * 
		       gX[ia1] * gY[ia2] * exp(-i PD) */
	    COMPLEX_ADD2 (ct1, SM1, SM2);
	    COMPLEX_MUL2 (ct2, Jm, ct1);
	    COMPLEX_ADD2 (ct1, SM3, SM4);
	    COMPLEX_MUL2 (ct3, Jp, ct1);
	    COMPLEX_ADD2 (ct2, ct2, ct3);
	    COMPLEX_MUL2 (DFDP,  ct2, ggPD);
	    /* part2 = -VYX */
	    COMPLEX_NEGATE(DFDP2, VYX);
	  }
	}
	break;
      case 3:     /* YX wrt Ey */
	if (isAnt1) {
	  if (antFit[ia1][paramNumber]) {
	    /* part = (0, -1) * {S[0] * CXc[ia2] * DPAc * CYc[ia1] +       
	                         S[1] * SXc[ia2] * SPAc * CYc[ia1] +
				 S[2] * CXc[ia2] * SPA  * SYc[ia1]  + 
				 S[3] * SXc[ia2] * DPA  * SYc[ia1]} * 
				 gY[ia1] * gX[ia2] * exp(-i PD) */
	    COMPLEX_MUL4(ct1, S[0], CXc[ia2], DPAc, CYc[ia1]);
	    COMPLEX_MUL4(ct2, S[1], SXc[ia2], SPAc, CYc[ia1]);
	    COMPLEX_MUL4(ct3, S[2], CXc[ia2], SPA,  SYc[ia1]);
	    COMPLEX_MUL4(ct4, S[3], SXc[ia2], DPA,  SYc[ia1]);
	    COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
	    COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
	    /* part2 = -VYX */
	    COMPLEX_NEGATE(DFDP2, VYX);
	  }
	}
	break;
      default:
	break;
      }; /* end antenna parameter switch */
      /* end antenna param */
    } else if (paramType==polnParmSou) {   /* Source parameters */
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* YX wrt IPol */
	if (souFit[isou][paramNumber]) {
	  /* part =   (MC1 + MC4) * gY[ia1] * gX[ia2] * exp(-i PD) */
	  COMPLEX_ADD2(ct1, MC1, MC4);
	  COMPLEX_MUL2(DFDP, ct1, ggPD);
	}
	break;
      case 1:     /* YX wrt QPol */
	if (souFit[isou][paramNumber]) {
	  /* part = (MC2 + MC3) * gY[ia1] * gX[ia2] * exp(-i PD)  */
	  COMPLEX_ADD2(ct1, MC2, MC3);
	  COMPLEX_MUL2(DFDP, ct1, ggPD);
	}
	break;
      case 2:     /* YX wrt UPol */
	if (souFit[isou][paramNumber]) {
	  /* i (MC2 - MC3) * gY[ia1] * gX[ia2] * exp(-i PD)  */
	  COMPLEX_SUB(ct1, MC2, MC3);
	  COMPLEX_MUL3(DFDP, Jp, ct1, ggPD);
	}
	break;
      case 3:     /* YX wrt Vpol */
	if (souFit[isou][paramNumber]) {
	  /* (MC1 - MC4) * gY[ia1] * gX[ia2] * exp(-i PD)  */
	  COMPLEX_SUB (ct1, MC1, MC4);
	  COMPLEX_MUL2(DFDP, ct1, ggPD);
	}
	break;  
      default:
	break;
      }; /* end source parameter switch */
      /* end source param */
    } else if (paramType==polnParmGain) {   /* Antenna Gains */
      switch (paramNumber) {   /* Switch over parameter */
      case 0:     /* YX wrt gX */
	if (isAnt2 && (antGainFit[ia2][paramNumber])) {
	  /* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) 
                           * gY[ia1]  * exp(-i PD) */
	  COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct2, antGain[ia1*2+1], 0);
	  COMPLEX_EXP (ct3, -PD);
	  COMPLEX_MUL3 (DFDP, ct1, ct2, ct3);
	  /* part2 = 0 */
	}
	break;
      case 1:     /* YX wrt gY */
	if  (isAnt1 && (antGainFit[ia1][paramNumber])) {
	  /* part =  (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * 
	                 gX[ia2] * exp(-i PD) */
	  COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct2, antGain[ia2*2+0], 0);
	  COMPLEX_EXP (ct3, -PD);
	  COMPLEX_MUL3 (DFDP, ct1, ct2, ct3);
	}
	break;
       default:
	break;
      } /* end gain switch */
  } else if (paramType==polnParmPD) {   /* X-Y phase difference */
      /* part = (0, -1) * VYX */   
      COMPLEX_MUL2 (DFDP, Jm, VYX);
      /* part2 = -VYX */
      COMPLEX_NEGATE(DFDP2, VYX);
  } /* end parameter types */

    /* Accumulate partials */
    if (paramType!=polnParmUnspec) {
      sumd  += 2.0 * isigma * (residR*DFDP.real + residI*DFDP.imag);
      sumd2 += 2.0*isigma * (DFDP.real*DFDP.real + DFDP.imag*DFDP.imag +
			      residR*DFDP2.real + residI*DFDP2.imag);
      } /* end set partials */
  }  /* end valid data */
  } /* End loop over visibilities */
  
  if (sumwt<=0.0) sumwt = 1.0;  /* Trap no data */
  args->ChiSq       = sum;   /* Save results  */
  args->sumParResid = sumParResid;
  args->sumXResid   = sumXResid;
  args->sumWt       = sumwt;
  args->nPobs       = nPobs;
  args->nXobs       = nXobs;
  if (paramType!=polnParmUnspec) args->sumDeriv  = sumd;
  if (paramType!=polnParmUnspec) args->sumDeriv2 = sumd2;

 /* Indicate completion if threaded */
  if (args->ithread>=0)
    ObitThreadPoolDone (args->thread, (gpointer)&args->ithread);

  return NULL;

} /*  end ThreadPolnFitXYChi2 */

/**
 * Create threading arguments
 * \param in       Fitting object
 * \param err      Obit error stack object.
 */
static void MakePolnFitFuncArgs (ObitPolnCalFit *in, ObitErr *err)
{
  olong i, nVisPerThread;

  /* If they already exist, first delete */
  if (in->thArgs) KillPolnFitFuncArgs(in);

  /* How many threads? */
  in->nThread   = MAX (1, ObitThreadNumProc(in->thread));
  nVisPerThread = in->nvis / in->nThread;

  /* Initialize threadArg array */
  in->thArgs = g_malloc0(in->nThread*sizeof(PolnFitArg*));
  for (i=0; i<in->nThread; i++) 
    in->thArgs[i] = g_malloc0(sizeof(PolnFitArg)); 
  for (i=0; i<in->nThread; i++) {
    in->thArgs[i]->thread   = ObitThreadRef(in->thread);
    in->thArgs[i]->err      = ObitErrRef(err);
    in->thArgs[i]->ithread  = i;
    if (in->nThread<=1) in->thArgs[i]->ithread = -1;
    in->thArgs[i]->doError  = in->doError;
    in->thArgs[i]->lo       = i*nVisPerThread;     /* Zero rel first */
    in->thArgs[i]->hi       = (i+1)*nVisPerThread; /* Zero rel last */
    in->thArgs[i]->ndata    = in->nvis*8;
    in->thArgs[i]->inData   = in->inData;
    in->thArgs[i]->inWt     = in->inWt;
    in->thArgs[i]->souNo    = in->souNo;
    in->thArgs[i]->antNo    = in->antNo;
    in->thArgs[i]->nant     = in->nant;
    in->thArgs[i]->refAnt   = in->refAnt;
    in->thArgs[i]->doFitRL  = in->doFitRL;
    in->thArgs[i]->doFitI   = in->doFitI;
    in->thArgs[i]->doFitPol = in->doFitPol;
    in->thArgs[i]->doFitV   = in->doFitV;
    in->thArgs[i]->doFitGain= in->doFitGain;
    in->thArgs[i]->isCircFeed= in->isCircFeed;
    in->thArgs[i]->antParm  = in->antParm;
    in->thArgs[i]->antErr   = in->antErr;
    in->thArgs[i]->antFit   = in->antFit;
    in->thArgs[i]->antPNumb = in->antPNumb;
    in->thArgs[i]->antGain      = in->antGain;
    in->thArgs[i]->antGainErr   = in->antGainErr;
    in->thArgs[i]->antGainFit   = in->antGainFit;
    in->thArgs[i]->antGainPNumb = in->antGainPNumb;
    in->thArgs[i]->nsou     = in->nsou;
    in->thArgs[i]->souParm  = in->souParm;
    in->thArgs[i]->souErr   = in->souErr;
    in->thArgs[i]->souFit   = in->souFit;
    in->thArgs[i]->souPNumb = in->souPNumb;
    in->thArgs[i]->souIDs   = in->souIDs;
    in->thArgs[i]->isouIDs  = in->isouIDs;
    in->thArgs[i]->nparam   = in->nparam;
    in->thArgs[i]->freq     = in->inDesc->freq;
    in->thArgs[i]->RLPhase  = in->RLPhase;
    in->thArgs[i]->PPol     = in->PPol;
    in->thArgs[i]->Chan     = 0;
    in->thArgs[i]->IFno     = 0;
    in->thArgs[i]->ChiSq    = 0.0;
    in->thArgs[i]->maxAnt   = 100;
  } /* end loop over thread args */

  /* Make sure do all data */
  i = in->nThread-1;
  in->thArgs[i]->hi = MAX (in->thArgs[i]->hi, in->nvis);

} /* end MakePolnFitFuncArgs */

/**
 * Delete threading arguments
 * \param in       Fitting object
 */
static void KillPolnFitFuncArgs (ObitPolnCalFit *in)
{
  olong i;

  /* If they already exist? */
  if (in->thArgs==NULL) return;

  /* Loop over threadArg arrays */
  for (i=0; i<in->nThread; i++) {
    ObitThreadUnref(in->thArgs[i]->thread);
    ObitErrUnref(in->thArgs[i]->err);
    g_free(in->thArgs[i]);  in->thArgs[i] = NULL;
  } /* end loop over thread args */
  g_free(in->thArgs );  in->thArgs = NULL;
} /* end  KillPolnFitFuncArgs*/

/**
 * Check for crazy antenna solutions, reset defaults if so.
 * If the RMS deviation from the default elipticity exceeds 10 deg
 * the fit is deemed crazy.
 * Uses last valid set of source parameters in reset
 * \param in           Fitting object
 * \param err          Obit error stack object.
 * \return             TRUE if reset defaults, else FALSE
 */
static gboolean CheckCrazy(ObitPolnCalFit *in, ObitErr *err)
{
  gboolean crazy = FALSE;
  odouble sum, RMS;
  olong i, j;

  if (err->error) return crazy;  /* Error detected? */

  /* Get mean square difference from default elipticity */
  sum = 0.0;
  for (i=0; i<in->nant; i++) {
    if (in->isCircFeed) {
      /* Circular feeds ori_r, elip_r, ori_l, elip_l */
     sum += (in->antParm[i*4+1] - G_PI/4.0)*(in->antParm[i*4+1] - G_PI/4.0) + 
            (in->antParm[i*4+3] + G_PI/4.0)*(in->antParm[i*4+3] + G_PI/4.0);
    } else {
      /* Linear feeds  ori_x, elip_x, ori_y, elip_y,  */
      sum += (in->antParm[i*4+1])*(in->antParm[i*4+1]) + 
             (in->antParm[i*4+3])*(in->antParm[i*4+3]);
    }
  } /* end loop over antennas */

  /* mean */
  sum /= in->nant*2.0;
  RMS =  sqrt(sum)*RAD2DG;
  crazy = RMS > 10.0;

  /* Crazy? -> reset */
  if (crazy) {
    Obit_log_error(err, OBIT_InfoWarn, "Antenna solution crazy, RMS %lf, reseting defaults", RMS);
    for (i=0; i<in->nant; i++) {
      for (j=0; j<4; j++) in->antParm[i*4+j] = 0.0;
      if (in->isCircFeed) {
	/* Circular feeds ori_r, elip_r, ori_l, elip_l */
	in->antParm[i*4+0] = 0.0;
	in->antParm[i*4+1] = G_PI/4.0;
	in->antParm[i*4+2] =  0.0;
	in->antParm[i*4+3] =  -G_PI/4.0;
       } else {
	/* Linear feeds, if the antenna angles on sky are given in the AntennaList[0] use then, 
           otherwise assume Feeds are X, Y
	   ori_x, elip_x, ori_y, elip_y,  */
	if (fabs (in->AntLists[0]->ANlist[i]->FeedAPA-in->AntLists[0]->ANlist[i]->FeedBPA)<1.0) {
	  /* Assume X, Y */
	  in->antParm[i*4+0] = 0.0;
	  in->antParm[i*4+2] = +G_PI/2.0;
	  /* DEBUG KAT 
	  in->antParm[i*4+2] = 0.0;
	  in->antParm[i*4+0] = +G_PI/2.0; */
	} else {  /* Use what's in the AN table as initial convert to radians */
	  in->antParm[i*4+0] = in->AntLists[0]->ANlist[i]->FeedAPA*DG2RAD;
	  in->antParm[i*4+2] = in->AntLists[0]->ANlist[i]->FeedBPA*DG2RAD;
	} /* end if antenna angle given */
      }
    } /* end loop over antennas */

    /* Reset source parameters */
    for (i=0; i<in->nsou; i++) {
      for (j=0; j<4; j++) {
	if (in->lastSouParm[i*4+j]>0.0) 
	  in->souParm[i*4+j] = in->lastSouParm[i*4+j];
	else 	in->souParm[i*4+j] = 0.0;
      } 
    }
  } /* end reset */
  else {  /* OK - save source parameters */
    for (i=0; i<in->nsou; i++) {
      for (j=0; j<4; j++) in->lastSouParm[i*4+j] = in->souParm[i*4+j];
    }
  } /* end save parameters */
  
  return crazy;
} /* end CheckCrazy */

#ifdef HAVE_GSL
/**
 * Circular feed function evaluator for polarization fitting solver
 * Orientation/Ellipticity version
 * Evaluates (model-observed) / sigma
 * Function from 
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of polarization_terms
 * \param param   Function parameter structure (ObitPolnCalFit)
 * \param f       Vector of (model-obs)/sigma for data points
 * \return completion code GSL_SUCCESS=OK
 */
static int PolnFitFuncOERL (const gsl_vector *x, void *params, 
			    gsl_vector *f)
{
  ObitPolnCalFit *args = (ObitPolnCalFit*)params;
  ofloat *data, *wt;
  gboolean **antFit  = args->antFit;
  olong   **antPNumb = args->antPNumb;
  odouble  *antParm  = args->antParm;
  gboolean **souFit  = args->souFit;
  odouble  *souParm  = args->souParm;
  olong   **souPNumb = args->souPNumb;
  dcomplex   *RS     = args->RS;
  dcomplex   *RD     = args->RD;
  dcomplex   *LS     = args->LS;
  dcomplex   *LD     = args->LD;
  dcomplex   *RSc    = args->RSc;
  dcomplex   *RDc    = args->RDc;
  dcomplex   *LSc    = args->LSc;
  dcomplex   *LDc    = args->LDc;
  ofloat PD, chi1, chi2;
  double val;
  odouble ipol=0.0, qpol=0.0, upol=0.0, vpol=0.0;
  odouble residR=0.0, residI=0.0, isigma=0.0;
  olong k, iant, ia1, ia2, isou, idata, refAnt;
  olong isouLast=-999;
  dcomplex PRref, PLref, PPRL, PPLR, PA1, PA2, PA1c, PA2c, ct1, ct2;
  dcomplex S[4], VRR, VRL, VLR, VLL;
  ofloat root2;
  size_t i, j;

  /* Initialize output */
  val = 0.0;
  for (i=0; i<args->ndata; i++) gsl_vector_set(f, i, val);

  COMPLEX_SET (S[0], 0.0, 0.0);  /* Initialize poln vector */
  COMPLEX_SET (S[1], 0.0, 0.0);
  COMPLEX_SET (S[2], 0.0, 0.0);
  COMPLEX_SET (S[3], 0.0, 0.0);
  
  /* R-L phase difference  at reference antenna */
  if (args->doFitRL) {
    j = args->PDPNumb;
    PD = gsl_vector_get(x, j);
  } else PD = args->PD;

  /* get model parameters - first antenna */
  for (iant=0; iant<args->nant; iant++) {
    /* Loop over parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if ((antFit[iant][k]) && (args->gotAnt[iant])) {
	j = antPNumb[iant][k];
	antParm[iant*4+k] = gsl_vector_get(x, j);
      }
    } /* end loop over parameters */
  } /* end loop over antennas */

  /* Ref antenna - 0 rel */
  refAnt = MAX(0, args->refAnt-1);

  /* now source */
  for (isou=0; isou<args->nsou; isou++) {
    /* Loop over parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if (souFit[isou][k]) {
	j = souPNumb[isou][k];
	souParm[isou*4+k] = gsl_vector_get(x, j);
      }
    } /* end loop over parameters */
  } /* end loop over sources */

  /* data & wt pointers */
  data = args->inData;
  wt   = args->inWt;

  /* Injest model factorize into antenna components - 
     data in order Orientation R/X, Elipticity R/X, Orientation L/Y, Elipticity L/Y*/
  root2 = 1.0 / sqrt(2.0);
  /* Elipticity, Orientation terms */
  for (i=0; i<args->nant; i++) {
    COMPLEX_SET(RS[i], root2*(cos(antParm[i*4+1]) + sin(antParm[i*4+1])), 0.);
    COMPLEX_SET(ct1,   root2*(cos(antParm[i*4+1]) - sin(antParm[i*4+1])), 0.);
    COMPLEX_EXP(ct2, 2*antParm[i*4+0]);
    COMPLEX_MUL2 (RD[i], ct1, ct2);
    COMPLEX_SET (ct1, root2*(cos(antParm[i*4+3]) + sin(antParm[i*4+3])), 0.);
    COMPLEX_EXP (ct2, -2*antParm[i*4+2]);
    COMPLEX_MUL2 (LS[i], ct1, ct2);
    COMPLEX_SET (LD[i], root2*(cos(antParm[i*4+3]) - sin(antParm[i*4+3])), 0.);
    COMPLEX_CONJUGATE (RSc[i], RS[i]);
    COMPLEX_CONJUGATE (RDc[i], RD[i]);
    COMPLEX_CONJUGATE (LSc[i], LS[i]);
    COMPLEX_CONJUGATE (LDc[i], LD[i]);
  }

  /* Reference antenna phase terms */
  COMPLEX_EXP (PRref,  antParm[(args->refAnt-1)*4+0]);
  COMPLEX_EXP (PLref, -antParm[(args->refAnt-1)*4+2]+PD);
  COMPLEX_CONJUGATE (ct1, PLref);
  COMPLEX_MUL2 (PPRL, PRref, ct1);
  COMPLEX_CONJUGATE (ct1, PRref);
  COMPLEX_MUL2 (PPLR, PLref, ct1);

  /* Loop over data */
  i = 0;
  for (idata=0; idata<args->nvis; idata++) {
    /* Parallactic angle terms */
    chi1  = data[idata*10+0];   /* parallactic angle ant 1 */
    chi2  = data[idata*10+1];   /* parallactic angle ant 2 */
    COMPLEX_EXP (PA1, 2*chi1);
    COMPLEX_EXP (PA2, 2*chi2);
    COMPLEX_CONJUGATE (PA1c, PA1);
    COMPLEX_CONJUGATE (PA2c, PA2);

    isou  = MAX (0, args->souNo[idata]);    /* Source number */
    /* New source? get parameters */
    if (isou!=isouLast) {
      isouLast = isou;
      /* Source parameters */
      ipol = souParm[isou*4+0];
      /* Fitting or fixed? */
      if (args->souFit[isou][1]) 
	qpol = souParm[isou*4+1];
      else
	qpol = args->PPol[isou]*ipol*cos(args->RLPhase[isou]);
      if (args->souFit[isou][2]) 
	upol = souParm[isou*4+2];
      else
	upol = args->PPol[isou]*ipol*sin(args->RLPhase[isou]);
      vpol = souParm[isou*4+3];
      /* Complex Stokes array */
      COMPLEX_SET (S[0], ipol+vpol, 0.0);
      COMPLEX_SET (S[1], qpol,  upol);
      COMPLEX_SET (S[2], qpol, -upol);
      COMPLEX_SET (S[3], ipol-vpol, 0.0);
    }

    /* Antenna parameters (0 ref) */
    ia1    = args->antNo[idata*2+0];
    ia2    = args->antNo[idata*2+1]; 
    
    /* Loop over correlations calcularing residuals */
    for (k=0; k<4; k++) {
      isigma = wt[idata*4+k];
      if (k<2) isigma *= 0.3;  /* Downweight parallel hand */
      switch (k) { 
      case 0:     /* RR */
	if (wt[idata*4+k]>0.0) {
	  /* VRR = S[0] * RS[ia1] * RSc[ia2] +        
	           S[1] * RS[ia1] * RDc[ia2] * PA2c + 
		   S[2] * RD[ia1] * RSc[ia2] * PA1  + 
		   S[3] * RD[ia1] * RDc[ia2] * PA1  * PA2c; */
	  COMPLEX_MUL3 (VRR, S[0], RS[ia1],  RSc[ia2]);
	  COMPLEX_MUL4 (ct1, S[1], RS[ia1],  RDc[ia2],  PA2c);
	  COMPLEX_ADD2 (VRR, VRR,  ct1);
	  COMPLEX_MUL4 (ct1, S[2], RD[ia1], RSc[ia2], PA1);
	  COMPLEX_ADD2 (VRR, VRR,  ct1);
	  COMPLEX_MUL5 (ct1, S[3], RD[ia1], RDc[ia2], PA1, PA2c);
	  COMPLEX_ADD2 (VRR, VRR,  ct1);
	  residR = VRR.real - data[idata*10+(k+1)*2];
	  residI = VRR.imag - data[idata*10+(k+1)*2+1];
	} else  residR = residI = 0.0; /* Invalid data */
	break;
      case 1:     /* LL */
	if (wt[idata*4+k]>0.0) {
	  /* VLL = S[0] * LS[ia1] * LSc[ia2] * PA1c * PA2 +	
	           S[1] * LS[ia1] * LDc[ia2] * PA1c +
		   S[2] * LD[ia1] * LSc[ia2] * PA2  +
		   S[3] * LD[ia1] * LDc[ia2]; */
	  COMPLEX_MUL5 (VLL, S[0], LS[ia1], LSc[ia2], PA1c, PA2);
	  COMPLEX_MUL4 (ct1, S[1], LS[ia1], LDc[ia2], PA1c);
	  COMPLEX_ADD2 (VLL, VLL,  ct1);
	  COMPLEX_MUL4 (ct1, S[2], LD[ia1], LSc[ia2], PA2);
	  COMPLEX_ADD2 (VLL, VLL,  ct1);
	  COMPLEX_MUL3 (ct1, S[3], LD[ia1], LDc[ia2]);
	  COMPLEX_ADD2 (VLL, VLL,  ct1);
	  residR = VLL.real - data[idata*10+(k+1)*2];
	  residI = VLL.imag - data[idata*10+(k+1)*2+1];
	} else  residR = residI = 0.0; /* Invalid data */
	break;
      case 2:     /* RL */
	if (wt[idata*4+k]>0.0) {
	  /* VRL = PPRL * S[0] * RS[ia1] * LSc[ia2] * PA2 +
	           PPRL * S[1] * RS[ia1] * LDc[ia2] + 
		   PPRL * S[2] * RD[ia1] * LSc[ia2] * PA1 * PA2 +
		   PPRL * S[3] * RD[ia1] * LDc[ia2] * PA1; */
	  COMPLEX_MUL4 (VRL, S[0], RS[ia1], LSc[ia2], PA2);
	  COMPLEX_MUL3 (ct1, S[1], RS[ia1], LDc[ia2]);
	  COMPLEX_ADD2 (VRL, VRL,  ct1);
	  COMPLEX_MUL5 (ct1, S[2], RD[ia1], LSc[ia2],  PA1,  PA2);
	  COMPLEX_ADD2 (VRL, VRL,  ct1);
	  COMPLEX_MUL4 (ct1, S[3], RD[ia1], LDc[ia2],  PA1);
	  COMPLEX_ADD2 (ct2, VRL,  ct1);
	  COMPLEX_MUL2 (VRL, PPRL, ct2);
	  residR = VRL.real - data[idata*10+(k+1)*2];
	  residI = VRL.imag - data[idata*10+(k+1)*2+1];
	} else  residR = residI = 0.0; /* Invalid data */
	break;
      case 3:     /* LR */
	if (wt[idata*4+k]>0.0) {
	  /* VLR = PPLR * S[0] * LS[ia1] * RSc[ia2] * PA1c +
	           PPLR * S[1] * LS[ia1] * RDc[ia2] * PA1c * PA2c +
		   PPLR * S[2] * LD[ia1] * RSc[ia2] +
		   PPLR * S[3] * LD[ia1] * RDc[ia2] * PA2c */
	  COMPLEX_MUL4 (VLR, S[0], LS[ia1], RSc[ia2], PA1c);
	  COMPLEX_MUL5 (ct1, S[1], LS[ia1], RDc[ia2], PA1c,  PA2c);
	  COMPLEX_ADD2 (VLR, VLR,  ct1);
	  COMPLEX_MUL3 (ct1, S[2], LD[ia1], RSc[ia2]);
	  COMPLEX_ADD2 (VLR, VLR,  ct1);
	  COMPLEX_MUL4 (ct1, S[3], LD[ia1], RDc[ia2], PA2c);
	  COMPLEX_ADD2 (ct2, VLR,  ct1);
	  COMPLEX_MUL2 (VLR, PPLR, ct2);
	  residR = VLR.real - data[idata*10+(k+1)*2];
	  residI = VLR.imag - data[idata*10+(k+1)*2+1];
	} else  residR = residI = 0.0; /* Invalid data */
	break;
      default:
	break;
      }; /* end switch */
      gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
      gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
      i++;  /* Update datum number */
    } /* end loop over correlations */
  } /* End loop over visibilities */
  
  return GSL_SUCCESS;
} /*  end PolnFitFuncOERL */

/**
 * Circular feed Jacobian evaluator for polarization fitting solver
 * Orientation/Ellipticity version
 * Evaluates partial derivatives of model wrt each parameter
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of polarization_terms
 * \param param   Function parameter structure (ObitPolnCalFit)
 * \param J       Jacobian matrix J[data_point, parameter]
 * \return completion code GSL_SUCCESS=OK
 */
static int PolnFitJacOERL (const gsl_vector *x, void *params, 
			 gsl_matrix *J)
{
  ObitPolnCalFit *args = (ObitPolnCalFit*)params;
  ofloat *data, *wt;
  gboolean **antFit  = args->antFit;
  olong   **antPNumb = args->antPNumb;
  odouble  *antParm  = args->antParm;
  gboolean **souFit  = args->souFit;
  odouble  *souParm  = args->souParm;
  olong   **souPNumb = args->souPNumb;
  odouble    *SR     = args->SR;
  odouble    *DR     = args->DR;
  odouble    *SL     = args->SL;
  odouble    *DL     = args->DL;
  dcomplex   *PR     = args->PR;
  dcomplex   *PRc    = args->PRc;
  dcomplex   *PL     = args->PL;
  dcomplex   *PLc    = args->PLc;
  dcomplex   *RS     = args->RS;
  dcomplex   *RD     = args->RD;
  dcomplex   *LS     = args->LS;
  dcomplex   *LD     = args->LD;
  dcomplex   *RSc    = args->RSc;
  dcomplex   *RDc    = args->RDc;
  dcomplex   *LSc    = args->LSc;
  dcomplex   *LDc    = args->LDc;
  ofloat PD, chi1, chi2;
  double val;
  odouble ipol=0.0, qpol=0.0, upol=0.0, vpol=0.0;
  odouble modelR, modelI, residR, residI, gradR, gradI, isigma;
  olong k, kk, iant, ia1, ia2, isou, idata, refAnt;
  olong isouLast=-999;
  dcomplex PRref, PLref, PPRL, PPLR, mPPRL, mPPLR, PA1, PA2, PA1c, PA2c;
  dcomplex ct1, ct2, ct3, ct4, ct5, ct6;
  dcomplex S[4], VRR, VRL, VLR, VLL, DFDP, MC1, MC2, MC3, MC4;
  ofloat root2;
  size_t i, j;

   /* Initialize output */
  val = 0.0;
  for (i=0; i<args->ndata; i++) {
    for (j=0; j<args->nparam; j++) gsl_matrix_set(J, i, j, val);
  }

  COMPLEX_SET (S[0], 0.0, 0.0);  /* Initialize poln vector */
  COMPLEX_SET (S[1], 0.0, 0.0);
  COMPLEX_SET (S[2], 0.0, 0.0);
  COMPLEX_SET (S[3], 0.0, 0.0);
  COMPLEX_SET (MC1, 0.0, 0.0);  /* Other stuff */
  COMPLEX_SET (MC2, 0.0, 0.0);
  COMPLEX_SET (MC3, 0.0, 0.0);
  COMPLEX_SET (MC4, 0.0, 0.0);
  COMPLEX_SET (VRR, 0.0, 0.0);
  COMPLEX_SET (VLL, 0.0, 0.0);
  COMPLEX_SET (VLR, 0.0, 0.0);
  COMPLEX_SET (VRL, 0.0, 0.0);
  
  /* R-L phase difference  at reference antenna */
  if (args->doFitRL) {
    j = args->PDPNumb;
    PD = gsl_vector_get(x, j);
  } else PD = args->PD;
  
  /* get model parameters - first antenna */
  for (iant=0; iant<args->nant; iant++) {
    /* Loop over antenna parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if ((antFit[iant][k]) && (args->gotAnt[iant])) {
	j = antPNumb[iant][k];
	antParm[iant*4+k] = gsl_vector_get(x, j);
      }
    } /* end loop over antenna parameters */
  } /* end loop over antennas */
  
  /* Ref antenna - 0 rel */
  refAnt = MAX(0, args->refAnt-1);
  
  /* now source */
  for (isou=0; isou<args->nsou; isou++) {
    /* Loop over source parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if (souFit[isou][k]) {
	j = souPNumb[isou][k];
	souParm[isou*4+k] = gsl_vector_get(x, j);
      }
    } /* end loop over source parameters */
  } /* end loop over sources */
  
  /* data & wt pointers */
  data = args->inData;
  wt   = args->inWt;

  /* Injest model factorize into antenna components - 
     data in order Orientation R/X, Elipticity R/X, Orientation L/Y, Elipticity L/Y*/
  root2 = 1.0 / sqrt(2.0);
  /* Elipticity, Orientation terms */
  for (i=0; i<args->nant; i++) {
    SR[i] = cos(antParm[i*4+1]) + sin(antParm[i*4+1]);
    DR[i] = cos(antParm[i*4+1]) - sin(antParm[i*4+1]);
    SL[i] = cos(antParm[i*4+3]) + sin(antParm[i*4+3]);
    DL[i] = cos(antParm[i*4+3]) - sin(antParm[i*4+3]);
    COMPLEX_SET  (RS[i], root2*SR[i], 0.);
    COMPLEX_SET  (ct1,   root2*DR[i], 0.);
    COMPLEX_EXP  (PR[i], 2*antParm[i*4+0]);
    COMPLEX_CONJUGATE (PRc[i], PR[i]);
    COMPLEX_MUL2 (RD[i], ct1, PR[i]);
    COMPLEX_SET  (ct1, root2*SL[i], 0.);
    COMPLEX_EXP  (PL[i], -2*antParm[i*4+2]);
    COMPLEX_CONJUGATE (PLc[i], PL[i]);
    COMPLEX_MUL2 (LS[i], ct1, PL[i]);
    COMPLEX_SET  (LD[i], root2*DL[i], 0.);
    COMPLEX_CONJUGATE (RSc[i], RS[i]);
    COMPLEX_CONJUGATE (RDc[i], RD[i]);
    COMPLEX_CONJUGATE (LSc[i], LS[i]);
    COMPLEX_CONJUGATE (LDc[i], LD[i]);
  }

  /* Reference antenna phase terms */
  COMPLEX_EXP (PRref,  antParm[(args->refAnt-1)*4+0]);
  COMPLEX_EXP (PLref, -antParm[(args->refAnt-1)*4+2]+PD);
  COMPLEX_CONJUGATE (ct1, PLref);
  COMPLEX_MUL2 (PPRL, PRref, ct1);
  COMPLEX_CONJUGATE (ct1, PRref);
  COMPLEX_MUL2 (PPLR, PLref, ct1);
  COMPLEX_NEGATE(mPPRL, PPRL);
  COMPLEX_NEGATE(mPPLR, PPLR);

  /* Loop over data */
  i = 0;
  for (idata=0; idata<args->nvis; idata++) {
    /* Parallactic angle terms */
    chi1  = data[idata*10+0];   /* parallactic angle ant 1 */
    chi2  = data[idata*10+1];   /* parallactic angle ant 2 */
    COMPLEX_EXP (PA1, 2*chi1);
    COMPLEX_EXP (PA2, 2*chi2);
    COMPLEX_CONJUGATE (PA1c, PA1);
    COMPLEX_CONJUGATE (PA2c, PA2);

    isou  = MAX (0, args->souNo[idata]);    /* Source number */
    /* New source? get parameters */
    if (isou!=isouLast) {
      isouLast = isou;
      /* Source parameters */
      ipol = souParm[isou*4+0];
      /* Fitting or fixed? */
      if (args->souFit[isou][1]) 
	qpol = souParm[isou*4+1];
      else
	qpol = args->PPol[isou]*ipol*cos(args->RLPhase[isou]);
      if (args->souFit[isou][2]) 
	upol = souParm[isou*4+2];
      else
	upol = args->PPol[isou]*ipol*sin(args->RLPhase[isou]);
      vpol = souParm[isou*4+3];
      /* Complex Stokes array */
      COMPLEX_SET (S[0], ipol+vpol, 0.0);
      COMPLEX_SET (S[1], qpol,  upol);
      COMPLEX_SET (S[2], qpol, -upol);
      COMPLEX_SET (S[3], ipol-vpol, 0.0);
    }

    /* Antenna parameters */
    ia1    = args->antNo[idata*2+0];
    ia2    = args->antNo[idata*2+1]; 
  
    /* i = datum number */
    /* Loop over correlations calculating derivatives */
    for (kk=0; kk<4; kk++) {
      isigma = wt[idata*4+kk];
      if (kk<2) isigma *= 0.3;  /* Downweight parallel hand */
      switch (kk) { 
      case 0:     /* RR */
	if (wt[idata*4+kk]>0.0) {
	  /* VRR = S[0] * RS[ia1] * RSc[ia2] +        
	           S[1] * RS[ia1] * RDc[ia2] * PA2c + 
		   S[2] * RD[ia1] * RSc[ia2] * PA1  + 
		   S[3] * RD[ia1] * RDc[ia2] * PA1  * PA2c; */
	  COMPLEX_MUL2 (MC1, RS[ia1], RSc[ia2]);
	  COMPLEX_MUL2 (VRR, S[0], MC1);
	  COMPLEX_MUL3 (MC2, RS[ia1], RDc[ia2],  PA2c);
	  COMPLEX_MUL2 (ct1, S[1], MC2);
	  COMPLEX_ADD2 (VRR, VRR,  ct1);
	  COMPLEX_MUL3 (MC3, RD[ia1], RSc[ia2], PA1);
	  COMPLEX_MUL2 (ct1, S[2], MC3);
	  COMPLEX_ADD2 (VRR, VRR,  ct1);
	  COMPLEX_MUL4 (MC4, RD[ia1], RDc[ia2], PA1, PA2c);
	  COMPLEX_MUL2 (ct1, S[3], MC4);
	  COMPLEX_ADD2 (VRR, VRR,  ct1);
	  modelR = VRR.real; modelI = VRR.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	} else  residR = residI = 0.0; /* Invalid data */
	
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (0, 2 i) * (S[2]*MC3 + S[3]*MC4) */
		COMPLEX_MUL2(ct1, S[2], MC3);
		COMPLEX_MUL2(ct2, S[3], MC4);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0, 2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		gradR = DFDP.real;    /* RrrR wrt Or1 */
		gradI = DFDP.imag;    /* RrrI wrt Or1 */
	      } else gradR = gradI =0.0;     /* Invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Er1 */
	    /* Fitting? */
 	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt Er1 */
		/* part = r2 * DR[ia1] * (S[0] * RSc[ia2] + S[1] * RDc[ia2] * PA2c) -
		          r2 * SR[ia1] * PR[ia1] * (S[2] * RSc[ia2] * PA1 + S[3] * RDc[ia2] * PA1  * PA2c) */
		COMPLEX_MUL2(ct1, S[0], RSc[ia2]);
		COMPLEX_MUL3(ct2, S[1], RDc[ia2], PA2c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DR[ia1], 0.0);
		COMPLEX_MUL2(ct5, ct4, ct3);
		COMPLEX_MUL3(ct1, S[2], RSc[ia2], PA1);
		COMPLEX_MUL4(ct2, S[3], RDc[ia2], PA1, PA2c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SR[ia1], 0.0);
		COMPLEX_MUL3(ct6, ct4, PR[ia1], ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;    /* RrrR wrt Er1 */
		gradI = DFDP.imag;    /* RrrI wrt Er1 */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt OL1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt EL1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt Or1 */
		/* (0,-2 i) * (S[1]*MC2 + S[3]*MC4) */
		COMPLEX_MUL2(ct1, S[1], MC2);
		COMPLEX_MUL2(ct2, S[3], MC4);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0, -2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		gradR = DFDP.real;    /* RrrR wrt Ol2 */
		gradI = DFDP.imag;    /* RrrI wrt Ol2 */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Er2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt Er2 */
		/* part = r2 * DR[ia2] * (S[0] * RS[ia1] + S[2] * RD[ia1] * PA1) -
		          r2 * SR[ia2] * PRc[ia2] * (S[1] * RS[ia1] * PA2c + S[3] * RD[ia1] * PA1  * PA2c) */
		COMPLEX_MUL2(ct1, S[0], RS[ia1]);
		COMPLEX_MUL3(ct2, S[2], RD[ia1], PA1);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DR[ia2], 0.0);
		COMPLEX_MUL2(ct5, ct4, ct3);
		COMPLEX_MUL3(ct1, S[1], RS[ia1], PA2c );
		COMPLEX_MUL4(ct2, S[3], RD[ia1], PA1, PA2c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SR[ia2], 0.0);
		COMPLEX_MUL3(ct6, ct4, PRc[ia2], ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;    /* RrrR wrt Ol2 */
		gradI = DFDP.imag;    /* RrrI wrt Ol2 */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	}  /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt I */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = RS[ia1] * RSc[ia2] + RD[ia1] * RDc[ia2] * PA1 * PA2c */
		COMPLEX_MUL2(ct1, RS[ia1], RSc[ia2]);
		COMPLEX_MUL4(ct2, RD[ia1], RDc[ia2], PA1, PA2c);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RrrR wrt ipol */
		gradI = DFDP.imag;    /* RrrI wrt ipol */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt qpol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = RS[ia1] * RDc[ia2] * PA2c + RD[ia1] * RSc[ia2] * PA1 */
		COMPLEX_MUL3(ct1, RS[ia1], RDc[ia2], PA2c);
		COMPLEX_MUL3(ct2, RD[ia1], RSc[ia2], PA1);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RrrR wrt qpol */
		gradI = DFDP.imag;    /* RrrI wrt qpol */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt upol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = i (RS[ia1] * RDc[ia2] * PA2c - RD[ia1] * RSc[ia2] * PA1) */
		COMPLEX_MUL3(ct1, RS[ia1], RDc[ia2], PA2c);
		COMPLEX_MUL3(ct2, RD[ia1], RSc[ia2], PA1);
		COMPLEX_SUB (ct3, ct1, ct2);
 		COMPLEX_SET(ct4, 0.0, 1.0);
		COMPLEX_MUL2(DFDP, ct4, ct3);
		gradR = DFDP.real;    /* RrrR wrt upol */
		gradI = DFDP.imag;    /* RrrI wrt upol */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt V */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = RS[ia1] * RSc[ia2] - RD[ia1] * RDc[ia2] * PA1 * PA2c */
		COMPLEX_MUL2(ct1, RS[ia1], RSc[ia2]);
		COMPLEX_MUL4(ct2, RD[ia1], RDc[ia2], PA1, PA2c);
		COMPLEX_SUB (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RrrR wrt vpol */
		gradI = DFDP.imag;    /* RrrI wrt vpol */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */
      
	/* gradient wrt PD = 0 */
	if (args->doFitRL) {
	  gradR = gradI = 0.0;
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	}
	
	break;  /* End RR */
	    
      case 1:     /* LL */
	if (wt[idata*4+kk]>0.0) {
	  /* VLL = S[0] * LS[ia1] * LSc[ia2] * PA1c * PA2 +	
	           S[1] * LS[ia1] * LDc[ia2] * PA1c +
		   S[2] * LD[ia1] * LSc[ia2] * PA2  +
		   S[3] * LD[ia1] * LDc[ia2]; */
	  COMPLEX_MUL4 (MC1, LS[ia1], LSc[ia2], PA1c, PA2);
	  COMPLEX_MUL2 (VLL, S[0], MC1);
	  COMPLEX_MUL3 (MC2, LS[ia1], LDc[ia2], PA1c);
	  COMPLEX_MUL2 (ct1, S[1], MC2);
	  COMPLEX_ADD2 (VLL, VLL,  ct1);
	  COMPLEX_MUL3 (MC3, LD[ia1], LSc[ia2], PA2);
	  COMPLEX_MUL2 (ct1, S[2], MC3);
	  COMPLEX_ADD2 (VLL, VLL,  ct1);
	  COMPLEX_MUL2 (MC4, LD[ia1], LDc[ia2]);
	  COMPLEX_MUL2 (ct1, S[3], MC4);
	  COMPLEX_ADD2 (VLL, VLL,  ct1);
	  modelR = VLL.real; modelI = VLL.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	} else  residR = residI = 0.0; /* Invalid data */
	  
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt OR1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Er1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt Or1 */
		/* (0,-2 i) * (S[0]*MC1 + S[1]*MC2) */
		COMPLEX_MUL2(ct1, S[0], MC1);
		COMPLEX_MUL2(ct2, S[1], MC2);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0, -2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		gradR = DFDP.real;    /* RllR wrt Ol1 */
		gradI = DFDP.imag;    /* RllI wrt Ol1 */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt El1 */
		/* part = r2 * DL[ia1] * PL[ia1] * (S[0] * LSc[ia2] * PA1c * PA2 + 
                                                    S[1] * LDc[ia2] * PA1c) -
                          r2 * SL[ia1] * (S[2] * LSc[ia2] * PA2  + S[3] * LDc[ia2]) */
		COMPLEX_MUL4(ct1, S[0], LSc[ia2], PA1c, PA2);
		COMPLEX_MUL3(ct2, S[1], LDc[ia2], PA1c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DL[ia1], 0.0);
		COMPLEX_MUL3(ct5, ct4, PL[ia1], ct3);
		COMPLEX_MUL3(ct1, S[2], LSc[ia2], PA2);
		COMPLEX_MUL2(ct2, S[3], LDc[ia2]);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SL[ia1], 0.0);
		COMPLEX_MUL2(ct6, ct4, ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;    /* RllR wrt El1 */
		gradI = DFDP.imag;    /* RllI wrt El1 */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* Rll wrt E2r = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* Rll wrt Ol2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt Ol2 */
		/* (0, 2 i) * (S[0]*MC1 + S[2]*MC3) */
		COMPLEX_MUL2(ct1, S[0], MC1);
		COMPLEX_MUL2(ct2, S[2], MC3);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0, 2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		gradR = DFDP.real;    /* RllR wrt Ol2 */
		gradI = DFDP.imag;    /* RllI wrt Ol2 */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El2  */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt El2 */
		/* part = r2 * DL[ia2] * PLc[ia2] * (S[0] * LS[ia1] * PA1c * PA2 + 
                                                     S[2] * LD[ia1] * PA2) -
                          r2 * SL[ia2] * (S[1] * LS[ia1] * PA1c + S[3] * LD[ia1]) */
		COMPLEX_MUL4(ct1, S[0], LS[ia1], PA1c, PA2);
		COMPLEX_MUL3(ct2, S[2], LD[ia1], PA2);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DL[ia2], 0.0);
		COMPLEX_MUL3(ct5, ct4, PLc[ia2], ct3);
		COMPLEX_MUL3(ct1, S[1], LS[ia1], PA1c);
		COMPLEX_MUL2(ct2, S[3], LD[ia1]);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SL[ia2], 0.0);
		COMPLEX_MUL2(ct6, ct4, ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;    /* RllR wrt El2 */
		gradI = DFDP.imag;    /* RllI wrt El2 */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	} /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt I */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = LS[ia1] * LSc[ia2] * PA1c * PA2 + LD[ia1] * LDc[ia2] */
		COMPLEX_MUL4(ct1, LS[ia1], LSc[ia2], PA1c, PA2);
		COMPLEX_MUL2(ct2, LD[ia1], LDc[ia2]);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RLLR wrt ipol */
		gradI = DFDP.imag;    /* RLLI wrt ipol */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Qpol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (LS[ia1] * LDc[ia2] * PA1c + LD[ia1] * LSc[ia2] * PA2) */
		COMPLEX_MUL3(ct1, LS[ia1], LDc[ia2], PA1c);
		COMPLEX_MUL3(ct2, LD[ia1], LSc[ia2], PA2);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RllR wrt qpol */
		gradI = DFDP.imag;    /* RllI wrt qpol */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =  i( LS[ia1] * LDc[ia2] * PA1c - LD[ia1] * LSc[ia2] * PA2) */
		COMPLEX_MUL3(ct1, LS[ia1], LDc[ia2], PA1c);
		COMPLEX_MUL3(ct2, LD[ia1], LSc[ia2], PA2);
		COMPLEX_SUB (ct3, ct1, ct2);
		COMPLEX_SET(ct4, 0.0, 1.0);
		COMPLEX_MUL2(DFDP, ct4, ct3);
		gradR = DFDP.real;    /* RllR wrt upol */
		gradI = DFDP.imag;    /* RllI wrt upol */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt V */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = LS[ia1] * LSc[ia2] * PA1c * PA2 - LD[ia1] * LDc[ia2] */
		COMPLEX_MUL4(ct1, LS[ia1], LSc[ia2], PA1c, PA2);
		COMPLEX_MUL2(ct2, LD[ia1], LDc[ia2]);
		COMPLEX_SUB (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RllR wrt vpol */
		gradI = DFDP.imag;    /* RllI wrt vpol */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */

	/* gradient wrt PD - no effect */
	if (args->doFitRL) {
	  gradR = gradI = 0.0;
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	}
	
	break;  /* End LL */
	break;  /* End LL */

        case 2:     /* RL */
	if (wt[idata*4+kk]>0.0) {
	  /* VRL = PPRL * S[0] * RS[ia1] * LSc[ia2] * PA2 +
	           PPRL * S[1] * RS[ia1] * LDc[ia2] + 
		   PPRL * S[2] * RD[ia1] * LSc[ia2] * PA1 * PA2 +
		   PPRL * S[3] * RD[ia1] * LDc[ia2] * PA1; */
	  COMPLEX_MUL4 (MC1, PPRL, RS[ia1], LSc[ia2], PA2);
	  COMPLEX_MUL2 (VRL, S[0], MC1);
	  COMPLEX_MUL3 (MC2, PPRL, RS[ia1], LDc[ia2]);
	  COMPLEX_MUL2 (ct1, S[1], MC2);
	  COMPLEX_ADD2 (VRL, VRL,  ct1);
	  COMPLEX_MUL5 (MC3, PPRL, RD[ia1], LSc[ia2],  PA1,  PA2);
	  COMPLEX_MUL2 (ct1, S[2], MC3);
	  COMPLEX_ADD2 (VRL, VRL,  ct1);
	  COMPLEX_MUL4 (MC4, PPRL, RD[ia1], LDc[ia2],  PA1);
	  COMPLEX_MUL2 (ct1, S[3], MC4);
	  COMPLEX_ADD2 (VRL, VRL,  ct1);
	  modelR = VRL.real; modelI = VRL.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	} else  residR = residI = 0.0; /* Invalid data */

	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (0, 2 i) * (S[2]*MC3 + S[3]*MC4) */
		COMPLEX_MUL2(ct1, S[2], MC3);
		COMPLEX_MUL2(ct2, S[3], MC4);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0,  2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		/* If ia1==refant) */
		if (ia1==refAnt) {
		  COMPLEX_MUL2(ct2, ct1, VRL);
		  COMPLEX_ADD2(DFDP, DFDP, ct2);
		}
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Er1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = r2 * DR[ia1] * (PPRL * S[0] * LSc[ia2] * PA2 + PPRL * S[1] * LDc[ia2]) -
                          r2 * SR[ia1] * PR[ia1] * (PPRL * S[2] * LSc[ia2] * PA1 * PA2 +
			                           PPRL * S[3] * LDc[ia2] * PA1) */
		COMPLEX_MUL4(ct1, PPRL, S[0], LSc[ia2], PA2);
		COMPLEX_MUL3(ct2, PPRL, S[1], LDc[ia2]);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DR[ia1], 0.0);
		COMPLEX_MUL2(ct5, ct4, ct3);
		COMPLEX_MUL5(ct1, PPRL, S[2], LSc[ia2], PA1, PA2);
		COMPLEX_MUL4(ct2, PPRL, S[3], LDc[ia2], PA1);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SR[ia1], 0.0);
		COMPLEX_MUL3(ct6, ct4, PR[ia1], ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol1 =0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
      
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Er2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (0, 2 i) * (S[0]*MC1 + S[2]*MC3) */
		COMPLEX_MUL2(ct1, S[0], MC1);
		COMPLEX_MUL2(ct2, S[2], MC3);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0,  2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		/* If ia2==refant) */
		if (ia2==refAnt) {
		  COMPLEX_SET (ct1, 0.0, 2.0);
		  COMPLEX_MUL2(ct2, ct1, VRL);
		  COMPLEX_ADD2(DFDP, DFDP, ct2);
		}
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = r2 * DL[ia2] * PLc[ia2] * (PPRL * S[0] * RS[ia1] * PA2 + 
                                                     PPRL * S[2] * RD[ia1] * PA1 * PA2) -
			  r2 * SL[ia2] * (PPRL * S[1] * RS[ia1] + PPRL * S[3] * RD[ia1] * PA1) */
		COMPLEX_MUL4(ct1, PPRL, S[0], RS[ia1], PA2);
		COMPLEX_MUL5(ct2, PPRL, S[2], RD[ia1], PA1, PA2);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DL[ia2], 0.0);
		COMPLEX_MUL3(ct5, ct4, PLc[ia2], ct3);
		COMPLEX_MUL3(ct1, PPRL, S[1], RS[ia1]);
		COMPLEX_MUL4(ct2, PPRL, S[3], RD[ia1], PA1);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SL[ia2], 0.0);
		COMPLEX_MUL2(ct6, ct4, ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	} /* end loop over second antenna parameters */
  
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt I */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =  PPRL * RS[ia1] * LSc[ia2] * PA2 + PPRL * RD[ia1] * LDc[ia2] * PA1*/
		COMPLEX_MUL4(ct1, PPRL, RS[ia1], LSc[ia2], PA2 );
		COMPLEX_MUL4(ct2, PPRL, RD[ia1], LDc[ia2], PA1);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt QPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (PPRL * RS[ia1] * LDc[ia2] + PPRL * RD[ia1] * LSc[ia2] * PA1 * PA2) */
		COMPLEX_MUL3(ct1, PPRL, RS[ia1], LDc[ia2]);
		COMPLEX_MUL5(ct2, PPRL, RD[ia1], LSc[ia2], PA1, PA2);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = i(PPRL * RS[ia1] * LDc[ia2] - PPRL * RD[ia1] * LSc[ia2] * PA1 * PA2) */
		COMPLEX_MUL2(ct1, RS[ia1], LDc[ia2]);
		COMPLEX_MUL4(ct2, RD[ia1], LSc[ia2], PA1, PA2);
		COMPLEX_SUB (ct3, ct1, ct2);
		COMPLEX_SET(ct4, 0.0, 1.0);
		COMPLEX_MUL3(DFDP, ct4, PPRL, ct3);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt V */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = PPRL * RS[ia1] * LSc[ia2] * PA2 - PPRL * RD[ia1] * LDc[ia2] * PA1 */
		COMPLEX_MUL4(ct1, PPRL, RS[ia1], LSc[ia2], PA2);
		COMPLEX_MUL4(ct2, PPRL, RD[ia1], LDc[ia2], PA1);
		COMPLEX_SUB (DFDP, ct1, ct2);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */
      
	/* gradient wrt PD */
	if (args->doFitRL) {
	  if (wt[idata*4+kk]>0.0) {
	    COMPLEX_SET(ct1, 0.0, 1.0);
	    COMPLEX_MUL2(DFDP, ct1, VRL);
	    gradR = DFDP.real;
	    gradI = DFDP.imag;
	  } else gradR = gradI = 0.0;    /* invalid data */
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	}
	
	break;  /* End RL */
      case 3:     /* LR */
	if (wt[idata*4+kk]>0.0) {
	  /* VLR = PPLR * S[0] * LS[ia1] * RSc[ia2] * PA1c +
	           PPLR * S[1] * LS[ia1] * RDc[ia2] * PA1c * PA2c +
		   PPLR * S[2] * LD[ia1] * RSc[ia2] +
		   PPLR * S[3] * LD[ia1] * RDc[ia2] * PA2c */
	  COMPLEX_MUL4 (MC1, PPLR, LS[ia1], RSc[ia2], PA1c);
	  COMPLEX_MUL2 (VLR, S[0], MC1);
	  COMPLEX_MUL5 (MC2, PPLR, LS[ia1], RDc[ia2], PA1c,  PA2c);
	  COMPLEX_MUL2 (ct1, S[1], MC2);
	  COMPLEX_ADD2 (VLR, VLR,  ct1);
	  COMPLEX_MUL3 (MC3, PPLR, LD[ia1], RSc[ia2]);
	  COMPLEX_MUL2 (ct1, S[2], MC3);
	  COMPLEX_ADD2 (VLR, VLR,  ct1);
	  COMPLEX_MUL4 (MC4, PPLR, LD[ia1], RDc[ia2],  PA2c);
	  COMPLEX_MUL2 (ct1, S[3], MC4);
	  COMPLEX_ADD2 (VLR, VLR, ct1);
	  modelR = VLR.real; modelI = VLR.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	} else  residR = residI = 0.0; /* Invalid data */
	
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Er1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (0,-2 i) * (S[0]*MC1 + S[1]*MC2) */
		COMPLEX_MUL2(ct1, S[0], MC1);
		COMPLEX_MUL2(ct2, S[1], MC2);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0, -2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		/* If ia1==refant) */
		if (ia1==refAnt) {
		  COMPLEX_SET (ct1, 0.0, -2.0);
		  COMPLEX_MUL2(ct2, ct1, VLR);
		  COMPLEX_ADD2(DFDP, DFDP, ct2);
		}
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = r2 * DL[ia1] * PL[ia1] * (PPLR * S[0] * RSc[ia2] * PA1c + 
                                                    PPLR * S[1] * RDc[ia2] * PA1c * PA2c) -
                          r2 * SL[ia1] * (PPLR * S[2] * RSc[ia2] + PPLR * S[3] * RDc[ia2] * PA2c) */
		COMPLEX_MUL4(ct1, PPLR, S[0], RSc[ia2], PA1c);
		COMPLEX_MUL5(ct2, PPLR, S[1], RDc[ia2], PA1c, PA2c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DL[ia1], 0.0);
		COMPLEX_MUL3(ct5, ct4, PL[ia1], ct3);
		COMPLEX_MUL3(ct1, PPLR, S[2], RSc[ia2]);
		COMPLEX_MUL4(ct2, PPLR, S[3], RDc[ia2], PA2c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SL[ia1], 0.0);
		COMPLEX_MUL2(ct6, ct4, ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (0,-2 i) * (S[1]*MC2 + S[3]*MC4) */
		COMPLEX_MUL2(ct1, S[1], MC2);
		COMPLEX_MUL2(ct2, S[3], MC4);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0, -2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		/* If ia2==refant */
		if (ia2==refAnt) {
		  COMPLEX_MUL2(ct2, ct1, VLR);
		  COMPLEX_ADD2(DFDP, DFDP, ct2);
		}
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Er2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = r2 * DR[ia2] * (PPLR * S[0] * LS[ia1] * PA1c + PPLR * S[2] * LD[ia1]) -
                          r2 * SR[ia2] * PRc[ia2] * (PPLR * S[1] * LS[ia1] * PA1c * PA2c + 
                                                     PPLR * S[3] * LD[ia1] * PA2c) */
		COMPLEX_MUL4(ct1, PPLR, S[0], LS[ia1], PA1c);
		COMPLEX_MUL3(ct2, PPLR, S[2], LD[ia1]);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DR[ia2], 0.0);
		COMPLEX_MUL2(ct5, ct4, ct3);
		COMPLEX_MUL5(ct1, PPLR, S[1], LS[ia1], PA1c, PA2c);
		COMPLEX_MUL4(ct2, PPLR, S[3], LD[ia1], PA2c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SR[ia2], 0.0);
		COMPLEX_MUL3(ct6, ct4, PRc[ia2], ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	} /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt I */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = PPLR * LS[ia1] * RSc[ia2] * PA1c + PPLR * LD[ia1] * RDc[ia2] * PA2c */
		COMPLEX_MUL4(ct1, PPLR, LS[ia1], RSc[ia2], PA1c);
		COMPLEX_MUL4(ct2, PPLR, LD[ia1], RDc[ia2], PA2c);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt QPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (PPLR * LS[ia1] * RDc[ia2] * PA1c * PA2c + PPLR * LD[ia1] * RSc[ia2]) */
		COMPLEX_MUL5(ct1, PPLR, LS[ia1], RDc[ia2], PA1c, PA2c);
		COMPLEX_MUL3(ct2, PPLR, LD[ia1], RSc[ia2]);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = i (PPLR * LS[ia1] * RDc[ia2] * PA1c * PA2c - PPLR * LD[ia1] * RSc[ia2]) */
		COMPLEX_MUL4(ct1, LS[ia1], RDc[ia2], PA1c, PA2c);
		COMPLEX_MUL2(ct2, LD[ia1], RSc[ia2]);
		COMPLEX_SUB (ct3, ct1, ct2);
		COMPLEX_SET(ct4, 0.0, 1.0);
		COMPLEX_MUL3(DFDP, ct4, PPLR, ct3);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt Vpol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = PPLR * LS[ia1] * RSc[ia2] * PA1c - PPLR * LD[ia1] * RDc[ia2] * PA2c */
		COMPLEX_MUL4(ct1, PPLR, LS[ia1], RSc[ia2], PA1c);
		COMPLEX_MUL4(ct2, PPLR, LD[ia1], RDc[ia2], PA2c);
		COMPLEX_SUB (DFDP, ct1, ct2);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	break;
      default:
	break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */

	/* gradient wrt PD */
	if (args->doFitRL) {
	  if (wt[idata*4+kk]>0.0) {
	    COMPLEX_SET(ct1, 0.0, -1.0);
	    COMPLEX_MUL2(DFDP, ct1, VLR);
	    gradR = DFDP.real;
	    gradI = DFDP.imag;
	  } else gradR = gradI = 0.0;    /* invalid data */
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	}
	/* End LR */
      }; /* end switch over data correlation */
      i++;  /* Update complex datum number */
    } /* end loop over correlations */
  } /* End loop over visibilities */

  return GSL_SUCCESS;
} /*  end PolnFitJacOERL */

/**
 * Circular feed Function & Jacobian evaluator for polarization fitting solver
 * Orientation/Ellipticity for circular feeds version
 * Evaluates partial derivatives of model wrt each parameter
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of polarization_terms
 * \param param   Function parameter structure (ObitPolnCalFit)
 * \param f       Vector of (model-obs)/sigma for data points
 * \param J       Jacobian matrix J[data_point, parameter]
 * \return completion code GSL_SUCCESS=OK
 */
static int PolnFitFuncJacOERL (const gsl_vector *x, void *params, 
			       gsl_vector *f, gsl_matrix *J)
{
  ObitPolnCalFit *args = (ObitPolnCalFit*)params;
  ofloat *data, *wt;
  gboolean **antFit  = args->antFit;
  olong   **antPNumb = args->antPNumb;
  odouble  *antParm  = args->antParm;
  gboolean **souFit  = args->souFit;
  odouble  *souParm  = args->souParm;
  olong   **souPNumb = args->souPNumb;
  odouble    *SR     = args->SR;
  odouble    *DR     = args->DR;
  odouble    *SL     = args->SL;
  odouble    *DL     = args->DL;
  dcomplex   *PR     = args->PR;
  dcomplex   *PRc    = args->PRc;
  dcomplex   *PL     = args->PL;
  dcomplex   *PLc    = args->PLc;
  dcomplex   *RS     = args->RS;
  dcomplex   *RD     = args->RD;
  dcomplex   *LS     = args->LS;
  dcomplex   *LD     = args->LD;
  dcomplex   *RSc    = args->RSc;
  dcomplex   *RDc    = args->RDc;
  dcomplex   *LSc    = args->LSc;
  dcomplex   *LDc    = args->LDc;
  ofloat PD, chi1, chi2;
  double val;
  odouble ipol=0.0, qpol=0.0, upol=0.0, vpol=0.0;
  odouble residR, residI, gradR, gradI, modelR, modelI, isigma;
  olong k, kk, iant, ia1, ia2, isou, idata, refAnt;
  olong isouLast=-999;
  dcomplex PRref, PLref, PPRL, PPLR, mPPRL, mPPLR, PA1, PA2, PA1c, PA2c;
  dcomplex ct1, ct2, ct3, ct4, ct5, ct6;
  dcomplex S[4], VRR, VRL, VLR, VLL, DFDP, MC1, MC2, MC3, MC4;
  ofloat root2;
  size_t i, j;

   /* Initialize output */
  val = 0.0;
  for (i=0; i<args->ndata; i++) {
    gsl_vector_set(f, i, val);
    for (j=0; j<args->nparam; j++) gsl_matrix_set(J, i, j, val);
  }

  COMPLEX_SET (S[0], 0.0, 0.0);  /* Initialize poln vector */
  COMPLEX_SET (S[1], 0.0, 0.0);
  COMPLEX_SET (S[2], 0.0, 0.0);
  COMPLEX_SET (S[3], 0.0, 0.0);
  COMPLEX_SET (MC1, 0.0, 0.0);  /* Other stuff */
  COMPLEX_SET (MC2, 0.0, 0.0);
  COMPLEX_SET (MC3, 0.0, 0.0);
  COMPLEX_SET (MC4, 0.0, 0.0);
  COMPLEX_SET (VRR, 0.0, 0.0);
  COMPLEX_SET (VLL, 0.0, 0.0);
  COMPLEX_SET (VLR, 0.0, 0.0);
  COMPLEX_SET (VRL, 0.0, 0.0);
  
  /* R-L phase difference  at reference antenna */
  if (args->doFitRL) {
    j = args->PDPNumb;
    PD = gsl_vector_get(x, j);
  } else PD = args->PD;
  
  /* get model parameters - first antenna */
  for (iant=0; iant<args->nant; iant++) {
    /* Loop over antenna parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if ((antFit[iant][k]) && (args->gotAnt[iant])) {
	j = antPNumb[iant][k];
	antParm[iant*4+k] = gsl_vector_get(x, j);
      }
    } /* end loop over antenna parameters */
  } /* end loop over antennas */
  
  /* Ref antenna - 0 rel */
  refAnt = MAX(0, args->refAnt-1);
  
  /* now source */
  for (isou=0; isou<args->nsou; isou++) {
    /* Loop over source parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if (souFit[isou][k]) {
	j = souPNumb[isou][k];
	souParm[isou*4+k] = gsl_vector_get(x, j);
      }
    } /* end loop over source parameters */
  } /* end loop over sources */
  
  /* data & wt pointers */
  data = args->inData;
  wt   = args->inWt;

  /* Injest model factorize into antenna components - 
     data in order Orientation R/X, Elipticity R/X, Orientation L/Y, Elipticity L/Y*/
  root2 = 1.0 / sqrt(2.0);
  /* Elipticity, Orientation terms */
  for (i=0; i<args->nant; i++) {
    SR[i] = cos(antParm[i*4+1]) + sin(antParm[i*4+1]);
    DR[i] = cos(antParm[i*4+1]) - sin(antParm[i*4+1]);
    SL[i] = cos(antParm[i*4+3]) + sin(antParm[i*4+3]);
    DL[i] = cos(antParm[i*4+3]) - sin(antParm[i*4+3]);
    COMPLEX_SET  (RS[i], root2*SR[i], 0.);
    COMPLEX_SET  (ct1,   root2*DR[i], 0.);
    COMPLEX_EXP  (PR[i], 2*antParm[i*4+0]);
    COMPLEX_CONJUGATE (PRc[i], PR[i]);
    COMPLEX_MUL2 (RD[i], ct1, PR[i]);
    COMPLEX_SET  (ct1, root2*SL[i], 0.);
    COMPLEX_EXP  (PL[i], -2*antParm[i*4+2]);
    COMPLEX_CONJUGATE (PLc[i], PL[i]);
    COMPLEX_MUL2 (LS[i], ct1, PL[i]);
    COMPLEX_SET  (LD[i], root2*DL[i], 0.);
    COMPLEX_CONJUGATE (RSc[i], RS[i]);
    COMPLEX_CONJUGATE (RDc[i], RD[i]);
    COMPLEX_CONJUGATE (LSc[i], LS[i]);
    COMPLEX_CONJUGATE (LDc[i], LD[i]);
  }

  /* Reference antenna phase terms */
  COMPLEX_EXP (PRref,  antParm[(args->refAnt-1)*4+0]);
  COMPLEX_EXP (PLref, -antParm[(args->refAnt-1)*4+2]+PD);
  COMPLEX_CONJUGATE (ct1, PLref);
  COMPLEX_MUL2 (PPRL, PRref, ct1);
  COMPLEX_CONJUGATE (ct1, PRref);
  COMPLEX_MUL2 (PPLR, PLref, ct1);
  COMPLEX_NEGATE(mPPRL, PPRL);
  COMPLEX_NEGATE(mPPLR, PPLR);

  /* Loop over data */
  i = 0;
  for (idata=0; idata<args->nvis; idata++) {
    /* Parallactic angle terms */
    chi1  = data[idata*10+0];   /* parallactic angle ant 1 */
    chi2  = data[idata*10+1];   /* parallactic angle ant 2 */
    COMPLEX_EXP (PA1, 2*chi1);
    COMPLEX_EXP (PA2, 2*chi2);
    COMPLEX_CONJUGATE (PA1c, PA1);
    COMPLEX_CONJUGATE (PA2c, PA2);

    isou  = MAX (0, args->souNo[idata]);    /* Source number */
    /* New source? get parameters */
    if (isou!=isouLast) {
      isouLast = isou;
      /* Source parameters */
      ipol = souParm[isou*4+0];
      /* Fitting or fixed? */
      if (args->souFit[isou][1]) 
	qpol = souParm[isou*4+1];
      else
	qpol = args->PPol[isou]*ipol*cos(args->RLPhase[isou]);
      if (args->souFit[isou][2]) 
	upol = souParm[isou*4+2];
      else
	upol = args->PPol[isou]*ipol*sin(args->RLPhase[isou]);
      vpol = souParm[isou*4+3];
      /* Complex Stokes array */
      COMPLEX_SET (S[0], ipol+vpol, 0.0);
      COMPLEX_SET (S[1], qpol,  upol);
      COMPLEX_SET (S[2], qpol, -upol);
      COMPLEX_SET (S[3], ipol-vpol, 0.0);
    }

    /* Antenna parameters */
    ia1    = args->antNo[idata*2+0];
    ia2    = args->antNo[idata*2+1]; 
    
    /* i = datum number */
    /* Loop over correlations calculating derivatives */
    for (kk=0; kk<4; kk++) {
      isigma = wt[idata*4+kk];
      if (kk<2) isigma *= 0.3;  /* Downweight parallel hand */
      switch (kk) { 
	
      case 0:     /* RR */
	if (wt[idata*4+kk]>0.0) {
	  /* VRR = S[0] * RS[ia1] * RSc[ia2] +        
	           S[1] * RS[ia1] * RDc[ia2] * PA2c + 
		   S[2] * RD[ia1] * RSc[ia2] * PA1  + 
		   S[3] * RD[ia1] * RDc[ia2] * PA1  * PA2c; */
	  COMPLEX_MUL2 (MC1, RS[ia1], RSc[ia2]);
	  COMPLEX_MUL2 (VRR, S[0], MC1);
	  COMPLEX_MUL3 (MC2, RS[ia1], RDc[ia2],  PA2c);
	  COMPLEX_MUL2 (ct1, S[1], MC2);
	  COMPLEX_ADD2 (VRR, VRR,  ct1);
	  COMPLEX_MUL3 (MC3, RD[ia1], RSc[ia2], PA1);
	  COMPLEX_MUL2 (ct1, S[2], MC3);
	  COMPLEX_ADD2 (VRR, VRR,  ct1);
	  COMPLEX_MUL4 (MC4, RD[ia1], RDc[ia2], PA1, PA2c);
	  COMPLEX_MUL2 (ct1, S[3], MC4);
	  COMPLEX_ADD2 (VRR, VRR,  ct1);
	  modelR = VRR.real; modelI = VRR.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
	  gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
	} else  residR = residI = 0.0; /* Invalid data */
	
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (0, 2 i) * (S[2]*MC3 + S[3]*MC4) */
		COMPLEX_MUL2(ct1, S[2], MC3);
		COMPLEX_MUL2(ct2, S[3], MC4);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0, 2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		gradR = DFDP.real;    /* RrrR wrt Or1 */
		gradI = DFDP.imag;    /* RrrI wrt Or1 */
	      } else gradR = gradI =0.0;     /* Invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Er1 */
	    /* Fitting? */
 	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt Er1 */
		/* part = r2 * DR[ia1] * (S[0] * RSc[ia2] + S[1] * RDc[ia2] * PA2c) -
		          r2 * SR[ia1] * PR[ia1] * (S[2] * RSc[ia2] * PA1 + S[3] * RDc[ia2] * PA1  * PA2c) */
		COMPLEX_MUL2(ct1, S[0], RSc[ia2]);
		COMPLEX_MUL3(ct2, S[1], RDc[ia2], PA2c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DR[ia1], 0.0);
		COMPLEX_MUL2(ct5, ct4, ct3);
		COMPLEX_MUL3(ct1, S[2], RSc[ia2], PA1);
		COMPLEX_MUL4(ct2, S[3], RDc[ia2], PA1, PA2c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SR[ia1], 0.0);
		COMPLEX_MUL3(ct6, ct4, PR[ia1], ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;    /* RrrR wrt Er1 */
		gradI = DFDP.imag;    /* RrrI wrt Er1 */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      /* gradR = -gradR; gradI = -gradI; DEBUG */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt OL1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt EL1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt Or1 */
		/* (0,-2 i) * (S[1]*MC2 + S[3]*MC4) */
		COMPLEX_MUL2(ct1, S[1], MC2);
		COMPLEX_MUL2(ct2, S[3], MC4);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0, -2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		gradR = DFDP.real;    /* RrrR wrt Ol2 */
		gradI = DFDP.imag;    /* RrrI wrt Ol2 */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Er2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt Er1 */
		/* part = r2 * DR[ia2] * (S[0] * RS[ia1] + S[2] * RD[ia1] * PA1) -
		          r2 * SR[ia2] * PRc[ia2] * (S[1] * RS[ia1] * PA2c + S[3] * RD[ia1] * PA1  * PA2c) */
		COMPLEX_MUL2(ct1, S[0], RS[ia1]);
		COMPLEX_MUL3(ct2, S[2], RD[ia1], PA1);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DR[ia2], 0.0);
		COMPLEX_MUL2(ct5, ct4, ct3);
		COMPLEX_MUL3(ct1, S[1], RS[ia1], PA2c );
		COMPLEX_MUL4(ct2, S[3], RD[ia1], PA1, PA2c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SR[ia2], 0.0);
		COMPLEX_MUL3(ct6, ct4, PRc[ia2], ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;    /* RrrR wrt Ol2 */
		gradI = DFDP.imag;    /* RrrI wrt Ol2 */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	}  /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt IPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = RS[ia1] * RSc[ia2] + RD[ia1] * RDc[ia2] * PA1 * PA2c */
		COMPLEX_MUL2(ct1, RS[ia1], RSc[ia2]);
		COMPLEX_MUL4(ct2, RD[ia1], RDc[ia2], PA1, PA2c);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RrrR wrt ipol */
		gradI = DFDP.imag;    /* RrrI wrt ipol */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt QPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (RS[ia1] * RDc[ia2] * PA2c + RD[ia1] * RSc[ia2] * PA1) */
		COMPLEX_MUL3(ct1, RS[ia1], RDc[ia2], PA2c);
		COMPLEX_MUL3(ct2, RD[ia1], RSc[ia2], PA1);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RrrR wrt qpol */
		gradI = DFDP.imag;    /* RrrI wrt qpol */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = i (RS[ia1] * RDc[ia2] * PA2c - RD[ia1] * RSc[ia2] * PA1) */
		COMPLEX_MUL3(ct1, RS[ia1], RDc[ia2], PA2c);
		COMPLEX_MUL3(ct2, RD[ia1], RSc[ia2], PA1);
		COMPLEX_SUB (ct3, ct1, ct2);
 		COMPLEX_SET(ct4, 0.0, 1.0);
		COMPLEX_MUL2(DFDP, ct4, ct3);
 		gradR = DFDP.real;    /* RrrR wrt upol */
		gradI = DFDP.imag;    /* RrrI wrt upol */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt VPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = RS[ia1] * RSc[ia2] - RD[ia1] * RDc[ia2] * PA1 * PA2c */
		COMPLEX_MUL2(ct1, RS[ia1], RSc[ia2]);
		COMPLEX_MUL4(ct2, RD[ia1], RDc[ia2], PA1, PA2c);
		COMPLEX_SUB (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RrrR wrt vpol */
		gradI = DFDP.imag;    /* RrrI wrt vpol */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */

	/* gradient wrt PD = 0 */
	if (args->doFitRL) {
	  gradR = gradI = 0.0;
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	}
	
	break;  /* End RR */
	
      case 1:     /* LL */
	if (wt[idata*4+kk]>0.0) {
	  /* VLL = S[0] * LS[ia1] * LSc[ia2] * PA1c * PA2 +	
	           S[1] * LS[ia1] * LDc[ia2] * PA1c +
		   S[2] * LD[ia1] * LSc[ia2] * PA2  +
		   S[3] * LD[ia1] * LDc[ia2]; */
	  COMPLEX_MUL4 (MC1, LS[ia1], LSc[ia2], PA1c, PA2);
	  COMPLEX_MUL2 (VLL, S[0], MC1);
	  COMPLEX_MUL3 (MC2, LS[ia1], LDc[ia2], PA1c);
	  COMPLEX_MUL2 (ct1, S[1], MC2);
	  COMPLEX_ADD2 (VLL, VLL,  ct1);
	  COMPLEX_MUL3 (MC3, LD[ia1], LSc[ia2], PA2);
	  COMPLEX_MUL2 (ct1, S[2], MC3);
	  COMPLEX_ADD2 (VLL, VLL,  ct1);
	  COMPLEX_MUL2 (MC4, LD[ia1], LDc[ia2]);
	  COMPLEX_MUL2 (ct1, S[3], MC4);
	  COMPLEX_ADD2 (VLL, VLL,  ct1);
	  modelR = VLL.real; modelI = VLL.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
	  gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
	} else  residR = residI = 0.0; /* Invalid data */
	  
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt OR1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Er1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt Or1 */
		/* (0,-2 i) * (S[0]*MC1 + S[1]*MC2) */
		COMPLEX_MUL2(ct1, S[0], MC1);
		COMPLEX_MUL2(ct2, S[1], MC2);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0, -2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		gradR = DFDP.real;    /* RllR wrt Ol1 */
		gradI = DFDP.imag;    /* RllI wrt Ol1 */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt El1 */
		/* part = r2 * DL[ia1] * PL[ia1] * (S[0] * LSc[ia2] * PA1c * PA2 + 
                                                    S[1] * LDc[ia2] * PA1c) -
                          r2 * SL[ia1] * (S[2] * LSc[ia2] * PA2  + S[3] * LDc[ia2]) */
		COMPLEX_MUL4(ct1, S[0], LSc[ia2], PA1c, PA2);
		COMPLEX_MUL3(ct2, S[1], LDc[ia2], PA1c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DL[ia1], 0.0);
		COMPLEX_MUL3(ct5, ct4, PL[ia1], ct3);
		COMPLEX_MUL3(ct1, S[2], LSc[ia2], PA2);
		COMPLEX_MUL2(ct2, S[3], LDc[ia2]);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SL[ia1], 0.0);
		COMPLEX_MUL2(ct6, ct4, ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;    /* RllR wrt El1 */
		gradI = DFDP.imag;    /* RllI wrt El1 */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* Rll wrt E2r = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* Rll wrt Ol2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt Ol2 */
		/* (0, 2 i) * (S[0]*MC1 + S[2]*MC3) */
		COMPLEX_MUL2(ct1, S[0], MC1);
		COMPLEX_MUL2(ct2, S[2], MC3);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0, 2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		gradR = DFDP.real;    /* RllR wrt Ol2 */
		gradI = DFDP.imag;    /* RllI wrt Ol2 */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El2  */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* Derivative of model wrt El2 */
		/* part = r2 * DL[ia2] * PLc[ia2] * (S[0] * LS[ia1] * PA1c * PA2 + 
                                                     S[2] * LD[ia1] * PA2) -
                          r2 * SL[ia2] * (S[1] * LS[ia1] * PA1c + S[3] * LD[ia1]) */
		COMPLEX_MUL4(ct1, S[0], LS[ia1], PA1c, PA2);
		COMPLEX_MUL3(ct2, S[2], LD[ia1], PA2);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DL[ia2], 0.0);
		COMPLEX_MUL3(ct5, ct4, PLc[ia2], ct3);
		COMPLEX_MUL3(ct1, S[1], LS[ia1], PA1c);
		COMPLEX_MUL2(ct2, S[3], LD[ia1]);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SL[ia2], 0.0);
		COMPLEX_MUL2(ct6, ct4, ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;    /* RllR wrt El2 */
		gradI = DFDP.imag;    /* RllI wrt El2 */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	} /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt IPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = LS[ia1] * LSc[ia2] * PA1c * PA2 + LD[ia1] * LDc[ia2] */
		COMPLEX_MUL4(ct1, LS[ia1], LSc[ia2], PA1c, PA2);
		COMPLEX_MUL2(ct2, LD[ia1], LDc[ia2]);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RLLR wrt ipol */
		gradI = DFDP.imag;    /* RLLI wrt ipol */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Qpol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (LS[ia1] * LDc[ia2] * PA1c + LD[ia1] * LSc[ia2] * PA2) */
		COMPLEX_MUL3(ct1, LS[ia1], LDc[ia2], PA1c);
		COMPLEX_MUL3(ct2, LD[ia1], LSc[ia2], PA2);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RllR wrt qpol */
		gradI = DFDP.imag;    /* RllI wrt qpol */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =  i (LS[ia1] * LDc[ia2] * PA1c - LD[ia1] * LSc[ia2] * PA2) */
		COMPLEX_MUL3(ct1, LS[ia1], LDc[ia2], PA1c);
		COMPLEX_MUL3(ct2, LD[ia1], LSc[ia2], PA2);
		COMPLEX_SUB (ct3, ct1, ct2);
 		COMPLEX_SET(ct4, 0.0, 1.0);
		COMPLEX_MUL2(DFDP, ct4, ct3);
		gradR = DFDP.real;    /* RllR wrt upol */
		gradI = DFDP.imag;    /* RllI wrt upol */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt VPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = LS[ia1] * LSc[ia2] * PA1c * PA2 - LD[ia1] * LDc[ia2] */
		COMPLEX_MUL4(ct1, LS[ia1], LSc[ia2], PA1c, PA2);
		COMPLEX_MUL2(ct2, LD[ia1], LDc[ia2]);
		COMPLEX_SUB(DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RllR wrt vpol */
		gradI = DFDP.imag;    /* RllI wrt vpol */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */

	/* gradient wrt PD - no effect */
	if (args->doFitRL) {
	  gradR = gradI = 0.0;
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	}
	
	break;  /* End LL */
	
      case 2:     /* RL */
	if (wt[idata*4+kk]>0.0) {
	  /* VRL = PPRL * S[0] * RS[ia1] * LSc[ia2] * PA2 +
	           PPRL * S[1] * RS[ia1] * LDc[ia2] + 
		   PPRL * S[2] * RD[ia1] * LSc[ia2] * PA1 * PA2 +
		   PPRL * S[3] * RD[ia1] * LDc[ia2] * PA1; */
	  COMPLEX_MUL4 (MC1, PPRL, RS[ia1], LSc[ia2], PA2);
	  COMPLEX_MUL2 (VRL, S[0], MC1);
	  COMPLEX_MUL3 (MC2, PPRL, RS[ia1], LDc[ia2]);
	  COMPLEX_MUL2 (ct1, S[1], MC2);
	  COMPLEX_ADD2 (VRL, VRL,  ct1);
	  COMPLEX_MUL5 (MC3, PPRL, RD[ia1], LSc[ia2],  PA1,  PA2);
	  COMPLEX_MUL2 (ct1, S[2], MC3);
	  COMPLEX_ADD2 (VRL, VRL,  ct1);
	  COMPLEX_MUL4 (MC4, PPRL, RD[ia1], LDc[ia2],  PA1);
	  COMPLEX_MUL2 (ct1, S[3], MC4);
	  COMPLEX_ADD2 (VRL, VRL,  ct1);
	  modelR = VRL.real; modelI = VRL.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
	  gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
	} else  residR = residI = 0.0; /* Invalid data */
	
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (0, 2 i) * (S[2]*MC3 + S[3]*MC4) */
		COMPLEX_MUL2(ct1, S[2], MC3);
		COMPLEX_MUL2(ct2, S[3], MC4);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0,  2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		/* If ia1==refant) */
		if (ia1==refAnt) {
		  COMPLEX_MUL2(ct2, ct1, VRL);
		  COMPLEX_ADD2(DFDP, DFDP, ct2);
		}
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      /* DEBUG 
	      if (ia1==22) fprintf (stdout, "ant 23 grad RL wrt Or1 %g %g\n",gradR, gradI);*/
	    }
	    break;
	  case 1:     /* wrt Er1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = r2 * DR[ia1] * (PPRL * S[0] * LSc[ia2] * PA2 + PPRL * S[1] * LDc[ia2]) -
                          r2 * SR[ia1] * PR[ia1] * (PPRL * S[2] * LSc[ia2] * PA1 * PA2 +
			                           PPRL * S[3] * LDc[ia2] * PA1) */
		COMPLEX_MUL4(ct1, PPRL, S[0], LSc[ia2], PA2);
		COMPLEX_MUL3(ct2, PPRL, S[1], LDc[ia2]);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DR[ia1], 0.0);
		COMPLEX_MUL2(ct5, ct4, ct3);
		COMPLEX_MUL5(ct1, PPRL, S[2], LSc[ia2], PA1, PA2);
		COMPLEX_MUL4(ct2, PPRL, S[3], LDc[ia2], PA1);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SR[ia1], 0.0);
		COMPLEX_MUL3(ct6, ct4, PR[ia1], ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol1 =0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Er2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (0, 2 i) * (S[0]*MC1 + S[2]*MC3) */
		COMPLEX_MUL2(ct1, S[0], MC1);
		COMPLEX_MUL2(ct2, S[2], MC3);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0,  2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		/* If ia2==refant) */
		if (ia2==refAnt) {
		  COMPLEX_SET (ct1, 0.0, 2.0);
		  COMPLEX_MUL2(ct2, ct1, VRL);
		  COMPLEX_ADD2(DFDP, DFDP, ct2);
		}
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = r2 * DL[ia2] * PLc[ia2] * (PPRL * S[0] * RS[ia1] * PA2 + 
                                                     PPRL * S[2] * RD[ia1] * PA1 * PA2) -
			  r2 * SL[ia2] * (PPRL * S[1] * RS[ia1] + PPRL * S[3] * RD[ia1] * PA1) */
		COMPLEX_MUL4(ct1, PPRL, S[0], RS[ia1], PA2);
		COMPLEX_MUL5(ct2, PPRL, S[2], RD[ia1], PA1, PA2);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DL[ia2], 0.0);
		COMPLEX_MUL3(ct5, ct4, PLc[ia2], ct3);
		COMPLEX_MUL3(ct1, PPRL, S[1], RS[ia1]);
		COMPLEX_MUL4(ct2, PPRL, S[3], RD[ia1], PA1);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SL[ia2], 0.0);
		COMPLEX_MUL2(ct6, ct4, ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	  } /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt IPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =  PPRL * RS[ia1] * LSc[ia2] * PA2 + PPRL * RD[ia1] * LDc[ia2] * PA1*/
		COMPLEX_MUL4(ct1, PPRL, RS[ia1], LSc[ia2], PA2);
		COMPLEX_MUL4(ct2, PPRL, RD[ia1], LDc[ia2], PA1);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt QPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (PPRL * RS[ia1] * LDc[ia2] - PPRL * RD[ia1] * LSc[ia2] * PA1 * PA2) */
		COMPLEX_MUL3(ct1, PPRL, RS[ia1], LDc[ia2]);
		COMPLEX_MUL5(ct2, PPRL, RD[ia1], LSc[ia2], PA1, PA2);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = i (PPRL * RS[ia1] * LDc[ia2] - PPRL * RD[ia1] * LSc[ia2] * PA1 * PA2) */
		COMPLEX_MUL2(ct1, RS[ia1], LDc[ia2]);
		COMPLEX_MUL4(ct2, RD[ia1], LSc[ia2], PA1, PA2);
		COMPLEX_SUB(ct3, ct1, ct2);
		COMPLEX_SET(ct4, 0.0, 1.0);
		COMPLEX_MUL3(DFDP, ct4, PPRL, ct3);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt VPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = PPRL * RS[ia1] * LSc[ia2] * PA2 - PPRL * RD[ia1] * LDc[ia2] * PA1 */
		COMPLEX_MUL4(ct1, PPRL, RS[ia1], LSc[ia2], PA2);
		COMPLEX_MUL4(ct2, PPRL, RD[ia1], LDc[ia2], PA1);
		COMPLEX_SUB (DFDP, ct1, ct2);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	   } 
	    break;
	  default:
	    break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */
	
	/* gradient wrt PD */
	if (args->doFitRL) {
	  if (wt[idata*4+kk]>0.0) {
	    COMPLEX_SET(ct1, 0.0, 1.0);
	    COMPLEX_MUL2(DFDP, ct1, VRL);
	    gradR = DFDP.real;
	    gradI = DFDP.imag;
	  } else gradR = gradI = 0.0;    /* invalid data */
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	}
	
	break;  /* End RL */
	
      case 3:     /* LR */
	if (wt[idata*4+kk]>0.0) {
	  /* VLR = PPLR * S[0] * LS[ia1] * RSc[ia2] * PA1c +
	           PPLR * S[1] * LS[ia1] * RDc[ia2] * PA1c * PA2c +
		   PPLR * S[2] * LD[ia1] * RSc[ia2] +
		   PPLR * S[3] * LD[ia1] * RDc[ia2] * PA2c */
	  COMPLEX_MUL4 (MC1, PPLR, LS[ia1], RSc[ia2], PA1c);
	  COMPLEX_MUL2 (VLR, S[0], MC1);
	  COMPLEX_MUL5 (MC2, PPLR, LS[ia1], RDc[ia2], PA1c,  PA2c);
	  COMPLEX_MUL2 (ct1, S[1], MC2);
	  COMPLEX_ADD2 (VLR, VLR,  ct1);
	  COMPLEX_MUL3 (MC3, PPLR, LD[ia1], RSc[ia2]);
	  COMPLEX_MUL2 (ct1, S[2], MC3);
	  COMPLEX_ADD2 (VLR, VLR,  ct1);
	  COMPLEX_MUL4 (MC4, PPLR, LD[ia1], RDc[ia2], PA2c);
	  COMPLEX_MUL2 (ct1, S[3], MC4);
	  COMPLEX_ADD2 (VLR, VLR, ct1);
	  modelR = VLR.real; modelI = VLR.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
	  gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
	} else  residR = residI = 0.0; /* Invalid data */
	
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Er1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (0,-2 i) * (S[0]*MC1 + S[1]*MC2) */
		COMPLEX_MUL2(ct1, S[0], MC1);
		COMPLEX_MUL2(ct2, S[1], MC2);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0, -2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		/* If ia1==refant) */
		if (ia1==refAnt) {
		  COMPLEX_MUL2(ct2, ct1, VLR);
		  COMPLEX_ADD2(DFDP, DFDP, ct2);
		}
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = r2 * DL[ia1] * PL[ia1] * (PPLR * S[0] * RSc[ia2] * PA1c + 
                                                    PPLR * S[1] * RDc[ia2] * PA1c * PA2c) -
                          r2 * SL[ia1] * (PPLR * S[2] * RSc[ia2] + PPLR * S[3] * RDc[ia2] * PA2c) */
		COMPLEX_MUL4(ct1, PPLR, S[0], RSc[ia2], PA1c);
		COMPLEX_MUL5(ct2, PPLR, S[1], RDc[ia2], PA1c, PA2c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DL[ia1], 0.0);
		COMPLEX_MUL3(ct5, ct4, PL[ia1], ct3);
		COMPLEX_MUL3(ct1, PPLR, S[2], RSc[ia2]);
		COMPLEX_MUL4(ct2, PPLR, S[3], RDc[ia2], PA2c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SL[ia1], 0.0);
		COMPLEX_MUL2(ct6, ct4, ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (0,-2 i) * (S[1]*MC2 + S[3]*MC4) */
		COMPLEX_MUL2(ct1, S[1], MC2);
		COMPLEX_MUL2(ct2, S[3], MC4);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct1, 0.0, -2.0);
		COMPLEX_MUL2(DFDP, ct1, ct3);
		/* If ia2==refant */
		if (ia2==refAnt) {
		  COMPLEX_MUL2(ct2, ct1, VLR);
		  COMPLEX_ADD2(DFDP, DFDP, ct2);
		}
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      /* DEBUG 
	      if (ia2==22) fprintf (stdout, "ant 23 grad LR wrt Or2 %g %g\n",gradR, gradI);*/
	    }
	    break;
	  case 1:     /* wrt Er2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = r2 * DR[ia2] * (PPLR * S[0] * LS[ia1] * PA1c + PPLR * S[2] * LD[ia1]) -
                          r2 * SR[ia2] * PRc[ia2] * (PPLR * S[1] * LS[ia1] * PA1c * PA2c + 
                                                     PPLR * S[3] * LD[ia1] * PA2c) */
		COMPLEX_MUL4(ct1, PPLR, S[0], LS[ia1], PA1c);
		COMPLEX_MUL3(ct2, PPLR, S[2], LD[ia1]);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*DR[ia2], 0.0);
		COMPLEX_MUL2(ct5, ct4, ct3);
		COMPLEX_MUL5(ct1, PPLR, S[1], LS[ia1], PA1c, PA2c);
		COMPLEX_MUL4(ct2, PPLR, S[3], LD[ia1], PA2c);
		COMPLEX_ADD2(ct3, ct1, ct2);
		COMPLEX_SET (ct4, root2*SR[ia2], 0.0);
		COMPLEX_MUL3(ct6, ct4, PRc[ia2], ct3);
		COMPLEX_SUB (DFDP, ct5, ct6);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	} /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt IPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = PPLR * LS[ia1] * RSc[ia2] * PA1c + PPLR * LD[ia1] * RDc[ia2] * PA2c */
		COMPLEX_MUL4(ct1, PPLR, LS[ia1], RSc[ia2], PA1c);
		COMPLEX_MUL4(ct2, PPLR, LD[ia1], RDc[ia2], PA2c);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt QPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (PPLR * LS[ia1] * RDc[ia2] * PA1c * PA2c + PPLR * LD[ia1] * RSc[ia2]) */
		COMPLEX_MUL5(ct1, PPLR, LS[ia1], RDc[ia2], PA1c, PA2c);
		COMPLEX_MUL3(ct2, PPLR, LD[ia1], RSc[ia2]);
		COMPLEX_ADD2(DFDP, ct1, ct2);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = i (PPLR * LS[ia1] * RDc[ia2] * PA1c * PA2c - PPLR * LD[ia1] * RSc[ia2]) */
		COMPLEX_MUL4(ct1, LS[ia1], RDc[ia2], PA1c, PA2c);
		COMPLEX_MUL2(ct2, LD[ia1], RSc[ia2]);
		COMPLEX_SUB (ct3, ct1, ct2);
		COMPLEX_SET(ct4, 0.0, 1.0);
		COMPLEX_MUL3(DFDP, ct4, PPLR, ct3);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt VPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = PPLR * LS[ia1] * RSc[ia2] * PA1c - PPLR * LD[ia1] * RDc[ia2] * PA2c */
		COMPLEX_MUL4(ct1, PPLR, LS[ia1], RSc[ia2], PA1c);
		COMPLEX_MUL4(ct2, PPLR, LD[ia1], RDc[ia2], PA2c);
		COMPLEX_SUB (DFDP, ct1, ct2);
		gradR = DFDP.real;
		gradI = DFDP.imag;
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	break;
      default:
	break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */
	
	/* gradient wrt PD */
	if (args->doFitRL) {
	  if (wt[idata*4+kk]>0.0) {
	    COMPLEX_SET(ct1, 0.0, -1.0);
	    COMPLEX_MUL2(DFDP, ct1, VLR);
	    gradR = DFDP.real;
	    gradI = DFDP.imag;
	  } else gradR = gradI = 0.0;    /* invalid data */
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	}
	/* End LR */	
      }; /* end switch over data correlation */
      
      i++;  /* Update complex datum number */
    } /* end loop over correlations */
  } /* End loop over visibilities */

  return GSL_SUCCESS;
} /*  end PolnFitFuncJacOERL */

/**
 * Linear feed function evaluator for polarization fitting solver
 * Orientation/Ellipticity version
 * Evaluates (model-observed) / sigma
 * Function from 
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of polarization_terms
 * \param param   Function parameter structure (ObitPolnCalFit)
 * \param f       Vector of (model-obs)/sigma for data points
 * \return completion code GSL_SUCCESS=OK
 */
static int PolnFitFuncOEXY (const gsl_vector *x, void *params, 
			    gsl_vector *f)
{
  ObitPolnCalFit *args = (ObitPolnCalFit*)params;
  ofloat *data, *wt;
  gboolean  **antFit     = args->antFit;
  olong     **antPNumb   = args->antPNumb;
  odouble    *antGain    = args->antGain;
  gboolean  **antGainFit = args->antGainFit;
  olong     **antGainPNumb  = args->antGainPNumb;
  odouble    *antParm    = args->antParm;
  gboolean  **souFit     = args->souFit;
  odouble    *souParm    = args->souParm;
  olong     **souPNumb   = args->souPNumb;
  dcomplex   *CX     = args->RS;
  dcomplex   *SX     = args->RD;
  dcomplex   *CY     = args->LS;
  dcomplex   *SY     = args->LD;
  dcomplex   *CXc    = args->RSc;
  dcomplex   *SXc    = args->RDc;
  dcomplex   *CYc    = args->LSc;
  dcomplex   *SYc    = args->LDc;
  ofloat PD, chi1, chi2;
  double val;
  odouble ipol=0.0, qpol=0.0, upol=0.0, vpol=0.0;
  odouble residR, residI, modelR, modelI, isigma;
  olong k, kk, iant, ia1, ia2, isou, idata, refAnt;
  olong isouLast=-999;
  dcomplex  SPA, DPA, SPAc, DPAc, ggPD;
  dcomplex ct1, ct2, Jm, Jp;
  dcomplex S[4], VXX, VXY, VYX, VYY, MC1, MC2, MC3, MC4;
  dcomplex SM1, SM2, SM3, SM4;
  size_t i, j;

   /* Initialize output */
  val = 0.0;
  for (i=0; i<args->ndata; i++) {
    gsl_vector_set(f, i, val);
  }

  COMPLEX_SET (S[0], 0.0, 0.0);  /* Initialize poln vector */
  COMPLEX_SET (S[1], 0.0, 0.0);
  COMPLEX_SET (S[2], 0.0, 0.0);
  COMPLEX_SET (S[3], 0.0, 0.0);
  COMPLEX_SET (MC1, 0.0, 0.0);  /* Other stuff */
  COMPLEX_SET (MC2, 0.0, 0.0);
  COMPLEX_SET (MC3, 0.0, 0.0);
  COMPLEX_SET (MC4, 0.0, 0.0);
  COMPLEX_SET (VXX, 0.0, 0.0);
  COMPLEX_SET (VYY, 0.0, 0.0);
  COMPLEX_SET (VYX, 0.0, 0.0);
  COMPLEX_SET (VXY, 0.0, 0.0);
  COMPLEX_SET (Jm,  0.0,-1.0);
  COMPLEX_SET (Jp,  0.0, 1.0);
  
  /* R-L phase difference  at reference antenna */
  if (args->doFitRL) {
    j = args->PDPNumb;
    PD = gsl_vector_get(x, j);
  } else PD = args->PD;
  
  /* get model parameters - first antenna */
  for (iant=0; iant<args->nant; iant++) {
    /* Loop over antenna parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if ((antFit[iant][k]) && (args->gotAnt[iant])) {
	j = antPNumb[iant][k];
	antParm[iant*4+k] = gsl_vector_get(x, j);
      }
    } /* end loop over antenna parameters */
  } /* end loop over antennas */
  
  /* antenna gain */
  for (iant=0; iant<args->nant; iant++) {
    /* Loop over antenna gains */
    for (k=0; k<2; k++) {
      /* Fitting? */
      if ((antGainFit[iant][k]) && (args->gotAnt[iant])) {
	j = antGainPNumb[iant][k];
	antGain[iant*2+k] = gsl_vector_get(x, j);
      }
    } /* end loop over antenna parameters */
  } /* end loop over antennas */
  
  /* Ref antenna - 0 rel */
  refAnt = MAX(0, args->refAnt-1);
  
  /* now source */
  for (isou=0; isou<args->nsou; isou++) {
    /* Loop over source parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if (souFit[isou][k]) {
	j = souPNumb[isou][k];
	souParm[isou*4+k] = gsl_vector_get(x, j);
      }
    } /* end loop over source parameters */
  } /* end loop over sources */
  

  /* data & wt pointers */
  data = args->inData;
  wt   = args->inWt;

   /* Injest model, factorize into antenna components - 
     data in order Orientation R/X, Elipticity R/X, Orientation L/Y, Elipticity L/Y */
  /* Elipticity, Orientation terms */
  for (i=0; i<args->nant; i++) {
    COMPLEX_EXP (ct1, -antParm[i*4+0]);
    COMPLEX_SET (ct2, cos(G_PI*0.25+antParm[i*4+1]), 0.0);
    COMPLEX_MUL3 (CX[i], Jp, ct1, ct2);
    COMPLEX_EXP (ct1, antParm[i*4+0]);
    COMPLEX_SET (ct2, sin(G_PI*0.25+antParm[i*4+1]), 0.0);
    COMPLEX_MUL2 (SX[i], ct1, ct2);
    COMPLEX_EXP (ct1, antParm[i*4+2]);
    COMPLEX_SET (ct2, cos(G_PI*0.25-antParm[i*4+3]), 0.0);
    COMPLEX_MUL2 (CY[i], ct1, ct2);
    COMPLEX_EXP (ct1, -antParm[i*4+2]);
    COMPLEX_SET (ct2, sin(G_PI*0.25-antParm[i*4+3]), 0.0);
    COMPLEX_MUL3 (SY[i], Jp, ct1, ct2);
    COMPLEX_CONJUGATE (CXc[i], CX[i]);
    COMPLEX_CONJUGATE (SXc[i], SX[i]);
    COMPLEX_CONJUGATE (CYc[i], CY[i]);
    COMPLEX_CONJUGATE (SYc[i], SY[i]);
  }

  /* Loop over data */
  i = 0;
  for (idata=0; idata<args->nvis; idata++) {
    /* Parallactic angle terms */
    chi1  = data[idata*10+0];   /* parallactic angle ant 1 */
    chi2  = data[idata*10+1];   /* parallactic angle ant 2 */
    COMPLEX_EXP (SPA,chi1+chi2);
    COMPLEX_EXP (DPA,chi1-chi2);
    COMPLEX_CONJUGATE (SPAc, SPA);
    COMPLEX_CONJUGATE (DPAc, DPA);

    isou  = MAX (0, args->souNo[idata]);    /* Source number */
    /* New source? get parameters */
    if (isou!=isouLast) {
      isouLast = isou;
      /* Source parameters */
      ipol = souParm[isou*4+0];
      /* Fitting or fixed? */
      if (args->souFit[isou][1]) 
	qpol = souParm[isou*4+1];
      else
	qpol = args->PPol[isou]*ipol*cos(args->RLPhase[isou]);
      if (args->souFit[isou][2]) 
	upol = souParm[isou*4+2];
      else
	upol = args->PPol[isou]*ipol*sin(args->RLPhase[isou]);
      vpol = souParm[isou*4+3];
      /* Complex Stokes array */
      COMPLEX_SET (S[0], ipol+vpol, 0.0);
      COMPLEX_SET (S[1], qpol,  upol);
      COMPLEX_SET (S[2], qpol, -upol);
      COMPLEX_SET (S[3], ipol-vpol, 0.0);
    }

    /* Antenna parameters */
    ia1    = args->antNo[idata*2+0];
    ia2    = args->antNo[idata*2+1]; 
    
    /* i = datum number */
    /* Loop over correlations calculating derivatives */
    for (kk=0; kk<4; kk++) {
      switch (kk) { 
	
      case 0:     /* XX */
	if (wt[idata*4+kk]>0.0) {
	  isigma = wt[idata*4+kk];
	  /* VXX = {S[0] * CX[ia1] * CXc[ia2] * DPAc  +
	            S[1] * CX[ia1] * SXc[ia2] * SPAc  +
		    S[2] * SX[ia1] * CXc[ia2] * SPA   + 
		    S[3] * SX[ia1] * SXc[ia2] * DPA} * g1X * g2X ;
	  */
	  COMPLEX_MUL3 (MC1, CX[ia1], CXc[ia2], DPAc);
	  COMPLEX_MUL3 (MC2, CX[ia1], SXc[ia2], SPAc);
	  COMPLEX_MUL3 (MC3, SX[ia1], CXc[ia2], SPA);
	  COMPLEX_MUL3 (MC4, SX[ia1], SXc[ia2], DPA);
	  COMPLEX_MUL2 (SM1, S[0], MC1);
	  COMPLEX_MUL2 (SM2, S[1], MC2);
	  COMPLEX_MUL2 (SM3, S[2], MC3);
	  COMPLEX_MUL2 (SM4, S[3], MC4);
	  COMPLEX_ADD4 (VXX, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ggPD,  antGain[ia1*2+0]*antGain[ia2*2+0], 0);
	  COMPLEX_MUL2 (VXX, VXX, ggPD);
	  modelR = VXX.real; modelI = VXX.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  /* residR = -residR; residI = - residI;   DEBUG */
	  gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
	  gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
	} 
	
	break;  /* End XX */
	
      case 1:     /* YY */
	if (wt[idata*4+kk]>0.0) {
	  isigma = wt[idata*4+kk];
	  /* VYY = {S[0] * SY[ia1] * SYc[ia2] * DPAc +       
	            S[1] * SY[ia1] * CYc[ia2] * SPAc +
		    S[2] * CY[ia1] * SYc[ia2] * SPA  + 
		    S[3] * CY[ia1] * CYc[ia2] * DPA} * g1Y * g2Y ;
	  */
	  COMPLEX_MUL3 (MC1, SY[ia1], SYc[ia2], DPAc);
	  COMPLEX_MUL3 (MC2, SY[ia1], CYc[ia2], SPAc);
	  COMPLEX_MUL3 (MC3, CY[ia1], SYc[ia2], SPA);
	  COMPLEX_MUL3 (MC4, CY[ia1], CYc[ia2], DPA);
	  COMPLEX_MUL2 (SM1, S[0], MC1);
	  COMPLEX_MUL2 (SM2, S[1], MC2);
	  COMPLEX_MUL2 (SM3, S[2], MC3);
	  COMPLEX_MUL2 (SM4, S[3], MC4);
	  COMPLEX_ADD4 (VYY, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ggPD,  antGain[ia1*2+1]*antGain[ia2*2+1], 0);
	  COMPLEX_MUL2 (VYY, VYY, ggPD);
	  modelR = VYY.real; modelI = VYY.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  /* residR = -residR; residI = - residI;   DEBUG */
	  gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
	  gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
	} 
	  
	break;  /* End YY */
	
      case 2:     /* XY */
	if (wt[idata*4+kk]>0.0) {
	  isigma = wt[idata*4+kk];
	  /* VXY = {S[0] * CX[ia1] * SYc[ia2] * DPAc +       
	            S[1] * CX[ia1] * CYc[ia2] * SPAc +
		    S[2] * SX[ia1] * SYc[ia2] * SPA  + 
		    S[3] * SX[ia1] * CYc[ia2] * DPA}} * g1X * g2Y * exp(i PD);
	  */
	  COMPLEX_MUL3 (MC1, CX[ia1], SYc[ia2], DPAc);
	  COMPLEX_MUL3 (MC2, CX[ia1], CYc[ia2], SPAc);
	  COMPLEX_MUL3 (MC3, SX[ia1], SYc[ia2], SPA);
	  COMPLEX_MUL3 (MC4, SX[ia1], CYc[ia2], DPA);
	  COMPLEX_MUL2 (SM1, S[0], MC1);
	  COMPLEX_MUL2 (SM2, S[1], MC2);
	  COMPLEX_MUL2 (SM3, S[2], MC3);
	  COMPLEX_MUL2 (SM4, S[3], MC4);
	  COMPLEX_ADD4 (VXY, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct1,  antGain[ia1*2+0]*antGain[ia2*2+1], 0);
	  COMPLEX_EXP (ct2, PD);
	  COMPLEX_MUL2 (ggPD, ct1, ct2);
	  COMPLEX_MUL2 (VXY, VXY, ggPD);
	  modelR = VXY.real; modelI = VXY.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  /* residR = -residR; residI = - residI;   DEBUG */
	  gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
	  gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
	} 
	
	break;  /* End XY */
      	
      case 3:     /* YX */
	if (wt[idata*4+kk]>0.0) {
	  isigma = wt[idata*4+kk];
	  /* VYX = {S[0] * SY[ia1] * CXc[ia2] * DPAc +       
	            S[1] * SY[ia1] * SXc[ia2] * SPAc +
		    S[2] * CY[ia1] * CXc[ia2] * SPA  + 
		    S[3] * CY[ia1] * SXc[ia2] * DPA} * g1Y * g2X * exp(-i PD)
	  */
	  COMPLEX_MUL3 (MC1, SY[ia1], CXc[ia2], DPAc);
	  COMPLEX_MUL3 (MC2, SY[ia1], SXc[ia2], SPAc);
	  COMPLEX_MUL3 (MC3, CY[ia1], CXc[ia2], SPA);
	  COMPLEX_MUL3 (MC4, CY[ia1], SXc[ia2], DPA);
	  COMPLEX_MUL2 (SM1, S[0], MC1);
	  COMPLEX_MUL2 (SM2, S[1], MC2);
	  COMPLEX_MUL2 (SM3, S[2], MC3);
	  COMPLEX_MUL2 (SM4, S[3], MC4);
	  COMPLEX_ADD4 (VYX, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct1,  antGain[ia1*2+1]*antGain[ia2*2+0], 0);
	  COMPLEX_EXP (ct2, -PD);
	  COMPLEX_MUL2 (ggPD, ct1, ct2);
	  COMPLEX_MUL2 (VYX, VYX, ggPD);
	  modelR = VYX.real; modelI = VYY.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  /* residR = -residR; residI = - residI;   DEBUG */
	  gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
	  gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
	} 
		
	break;  /* End YX */
	
      default:
	  break;
	}; /* end switch over data correlation */
      
      i++;  /* Update complex datum number */
    } /* end loop over correlations */
  } /* End loop over visibilities */

  return GSL_SUCCESS;
} /*  end PolnFitFuncOEXY */

/**
 * Linear feed Jacobian evaluator for polarization fitting solver
 * Orientation/Ellipticity version
 * Evaluates partial derivatives of model wrt each parameter
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of polarization_terms
 * \param param   Function parameter structure (ObitPolnCalFit)
 * \param J       Jacobian matrix J[data_point, parameter]
 * \return completion code GSL_SUCCESS=OK
 */
static int PolnFitJacOEXY (const gsl_vector *x, void *params, 
			   gsl_matrix *J)
{
  ObitPolnCalFit *args = (ObitPolnCalFit*)params;
  ofloat *data, *wt;
  gboolean  **antFit     = args->antFit;
  olong     **antPNumb   = args->antPNumb;
  odouble    *antParm    = args->antParm;
  odouble    *antGain    = args->antGain;
  gboolean  **antGainFit = args->antGainFit;
  olong     **antGainPNumb  = args->antGainPNumb;
  gboolean  **souFit     = args->souFit;
  odouble    *souParm    = args->souParm;
  olong     **souPNumb   = args->souPNumb;
  dcomplex   *CX     = args->RS;
  dcomplex   *SX     = args->RD;
  dcomplex   *CY     = args->LS;
  dcomplex   *SY     = args->LD;
  dcomplex   *CXc    = args->RSc;
  dcomplex   *SXc    = args->RDc;
  dcomplex   *CYc    = args->LSc;
  dcomplex   *SYc    = args->LDc;
  ofloat PD, chi1, chi2;
  double val;
  odouble ipol=0.0, qpol=0.0, upol=0.0, vpol=0.0;
  odouble residR=0.0, residI=0.0, gradR=0.0, gradI=0.0, modelR=0.0, modelI=0.0, isigma=0.0;
  olong k, kk, iant, ia1, ia2, isou, idata, refAnt;
  olong isouLast=-999;
  dcomplex  SPA, DPA, SPAc, DPAc, ggPD;
  dcomplex ct1, ct2, ct3, ct4, ct5, Jm, Jp;
  dcomplex S[4], VXX, VXY, VYX, VYY, MC1, MC2, MC3, MC4, DFDP;
  dcomplex SM1, SM2, SM3, SM4;
  size_t i, j;

   /* Initialize output */
  val = 0.0;
  for (i=0; i<args->ndata; i++) {
    for (j=0; j<args->nparam; j++) gsl_matrix_set(J, i, j, val);
  }

  COMPLEX_SET (S[0], 0.0, 0.0);  /* Initialize poln vector */
  COMPLEX_SET (S[1], 0.0, 0.0);
  COMPLEX_SET (S[2], 0.0, 0.0);
  COMPLEX_SET (S[3], 0.0, 0.0);
  COMPLEX_SET (MC1, 0.0, 0.0);  /* Other stuff */
  COMPLEX_SET (MC2, 0.0, 0.0);
  COMPLEX_SET (MC3, 0.0, 0.0);
  COMPLEX_SET (MC4, 0.0, 0.0);
  COMPLEX_SET (VXX, 0.0, 0.0);
  COMPLEX_SET (VYY, 0.0, 0.0);
  COMPLEX_SET (VYX, 0.0, 0.0);
  COMPLEX_SET (VXY, 0.0, 0.0);
  COMPLEX_SET (Jm,  0.0,-1.0);
  COMPLEX_SET (Jp,  0.0, 1.0);
  COMPLEX_SET (SM1, 0.0, 0.0);
  COMPLEX_SET (SM2, 0.0, 0.0);
  COMPLEX_SET (SM3, 0.0, 0.0);
  COMPLEX_SET (SM4, 0.0, 0.0);
  COMPLEX_SET (ggPD, 0.0, 0.0);
  
  /* R-L phase difference  at reference antenna */
  if (args->doFitRL) {
    j = args->PDPNumb;
    PD = gsl_vector_get(x, j);
  } else PD = args->PD;
  
  /* get model parameters - first antenna */
  for (iant=0; iant<args->nant; iant++) {
    /* Loop over antenna parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if ((antFit[iant][k]) && (args->gotAnt[iant])) {
	j = antPNumb[iant][k];
	antParm[iant*4+k] = gsl_vector_get(x, j);
      }
    } /* end loop over antenna parameters */
  } /* end loop over antennas */
  
  /* antenna gain */
  for (iant=0; iant<args->nant; iant++) {
    /* Loop over antenna gains */
    for (k=0; k<2; k++) {
      /* Fitting? */
      if ((antGainFit[iant][k]) && (args->gotAnt[iant])) {
	j = antGainPNumb[iant][k];
	antGain[iant*2+k] = gsl_vector_get(x, j);
      }
    } /* end loop over antenna parameters */
  } /* end loop over antennas */
  
  /* Ref antenna - 0 rel */
  refAnt = MAX(0, args->refAnt-1);
  
  /* now source */
  for (isou=0; isou<args->nsou; isou++) {
    /* Loop over source parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if (souFit[isou][k]) {
	j = souPNumb[isou][k];
	souParm[isou*4+k] = gsl_vector_get(x, j);
      }
    } /* end loop over source parameters */
  } /* end loop over sources */
  
  /* data & wt pointers */
  data = args->inData;
  wt   = args->inWt;

   /* Injest model, factorize into antenna components - 
     data in order Orientation R/X, Elipticity R/X, Orientation L/Y, Elipticity L/Y */
  /* Elipticity, Orientation terms */
  for (i=0; i<args->nant; i++) {
    COMPLEX_EXP (ct1, -antParm[i*4+0]);
    COMPLEX_SET (ct2, cos(G_PI*0.25+antParm[i*4+1]), 0.0);
    COMPLEX_MUL3 (CX[i], Jp, ct1, ct2);
    COMPLEX_EXP (ct1, antParm[i*4+0]);
    COMPLEX_SET (ct2, sin(G_PI*0.25+antParm[i*4+1]), 0.0);
    COMPLEX_MUL2 (SX[i], ct1, ct2);
    COMPLEX_EXP (ct1, antParm[i*4+2]);
    COMPLEX_SET (ct2, cos(G_PI*0.25-antParm[i*4+3]), 0.0);
    COMPLEX_MUL2 (CY[i], ct1, ct2);
    COMPLEX_EXP (ct1, -antParm[i*4+2]);
    COMPLEX_SET (ct2, sin(G_PI*0.25-antParm[i*4+3]), 0.0);
    COMPLEX_MUL3 (SY[i], Jp, ct1, ct2);
    COMPLEX_CONJUGATE (CXc[i], CX[i]);
    COMPLEX_CONJUGATE (SXc[i], SX[i]);
    COMPLEX_CONJUGATE (CYc[i], CY[i]);
    COMPLEX_CONJUGATE (SYc[i], SY[i]);
  }

  /* Loop over data */
  i = 0;
  for (idata=0; idata<args->nvis; idata++) {
    /* Parallactic angle terms */
    chi1  = data[idata*10+0];   /* parallactic angle ant 1 */
    chi2  = data[idata*10+1];   /* parallactic angle ant 2 */
    COMPLEX_EXP (SPA,chi1+chi2);
    COMPLEX_EXP (DPA,chi1-chi2);
    COMPLEX_CONJUGATE (SPAc, SPA);
    COMPLEX_CONJUGATE (DPAc, DPA);

    isou  = MAX (0, args->souNo[idata]);    /* Source number */
    /* New source? get parameters */
    if (isou!=isouLast) {
      isouLast = isou;
      /* Source parameters */
      ipol = souParm[isou*4+0];
      /* Fitting or fixed? */
      if (args->souFit[isou][1]) 
	qpol = souParm[isou*4+1];
      else
	qpol = args->PPol[isou]*ipol*cos(args->RLPhase[isou]);
      if (args->souFit[isou][2]) 
	upol = souParm[isou*4+2];
      else
	upol = args->PPol[isou]*ipol*sin(args->RLPhase[isou]);
      vpol = souParm[isou*4+3];
      /* Complex Stokes array */
      COMPLEX_SET (S[0], ipol+vpol, 0.0);
      COMPLEX_SET (S[1], qpol,  upol);
      COMPLEX_SET (S[2], qpol, -upol);
      COMPLEX_SET (S[3], ipol-vpol, 0.0);
    }

    /* Antenna parameters */
    ia1    = args->antNo[idata*2+0];
    ia2    = args->antNo[idata*2+1]; 
    
    /* i = datum number */
    /* Loop over correlations calculating derivatives */
    for (kk=0; kk<4; kk++) {
      switch (kk) { 
	
      case 0:     /* XX */
	if (wt[idata*4+kk]>0.0) {
	  isigma = wt[idata*4+kk];
	  /* VXX = {S[0] * CX[ia1] * CXc[ia2] * DPAc  +
	            S[1] * CX[ia1] * SXc[ia2] * SPAc  +
		    S[2] * SX[ia1] * CXc[ia2] * SPA   + 
		    S[3] * SX[ia1] * SXc[ia2] * DPA} * g1X * g2X ;
	  */
	  COMPLEX_MUL3 (MC1, CX[ia1], CXc[ia2], DPAc);
	  COMPLEX_MUL3 (MC2, CX[ia1], SXc[ia2], SPAc);
	  COMPLEX_MUL3 (MC3, SX[ia1], CXc[ia2], SPA);
	  COMPLEX_MUL3 (MC4, SX[ia1], SXc[ia2], DPA);
	  COMPLEX_MUL2 (SM1, S[0], MC1);
	  COMPLEX_MUL2 (SM2, S[1], MC2);
	  COMPLEX_MUL2 (SM3, S[2], MC3);
	  COMPLEX_MUL2 (SM4, S[3], MC4);
	  COMPLEX_ADD4 (VXX, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ggPD,  antGain[ia1*2+0]*antGain[ia2*2+0], 0);
	  COMPLEX_MUL2 (VXX, VXX, ggPD);
	  modelR = VXX.real; modelI = VXX.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  /* residR = -residR; residI = - residI;   DEBUG */
	} else  residR = residI = 0.0; /* Invalid data */
	
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* {(0, -1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		    (0,  1) * S[2] * MC3 + (0,  1) * S[3] * MC4}  * gX[ia1] * gX[ia2]  */
		COMPLEX_ADD2 (ct1, SM1, SM2);
		COMPLEX_MUL2 (ct2, Jm, ct1);
		COMPLEX_ADD2 (ct1, SM3, SM4);
		COMPLEX_MUL2 (ct3, Jp, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP, ct2, ggPD);
		gradR = DFDP.real;    /* RxxR wrt Ox1 */
		gradI = DFDP.imag;    /* RxxI wrt Ox1 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI =0.0;     /* Invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex1 */
	    /* Fitting? */
 	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) * {S[0] * CXc[ia2] * DPAc * SXc[ia1] +
		                     S[1] * SXc[ia2] * SPAc * SXc[ia1] +
				     S[2] * CXc[ia2] * SPA  * CXc[ia1] +	       
				     S[3] * SXc[ia2] * DPA  * CXc[ia1]}  * gX[ia1] * gX[ia2]  */
		COMPLEX_MUL4(ct1, S[0], CXc[ia2], DPAc, SXc[ia1]);
		COMPLEX_MUL4(ct2, S[1], SXc[ia2], SPAc, SXc[ia1]);
		COMPLEX_MUL4(ct3, S[2], CXc[ia2], SPA,  SXc[ia1]);
		COMPLEX_MUL4(ct4, S[3], SXc[ia2], DPA,  CXc[ia1]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;    /* RxxR wrt Ex1 */
		gradI = DFDP.imag;    /* RxxI wrt Ex1 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt OY1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt EY1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		    (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4}  * gX[ia1] * gX[ia2]   */
		COMPLEX_ADD2 (ct1, SM1, SM3);
		COMPLEX_MUL2 (ct2, Jp, ct1);
		COMPLEX_ADD2 (ct1, SM2, SM4);
		COMPLEX_MUL2 (ct3, Jm, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP, ct2, ggPD);
		gradR = DFDP.real;    /* RxxR wrt Ox2 */
		gradI = DFDP.imag;    /* RxxI wrt Ox2 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) * {S[0] * CX[ia1] * DPAc * SX[ia2] +
		                     S[1] * CX[ia1] * SPAc * CX[ia2] +
				     S[2] * SX[ia1] * SPA  * SX[ia2] +
				     S[3] * SX[ia1] * DPA  * CX[ia2]} * gX[ia1] * gX[ia2]  */
		COMPLEX_MUL4(ct1, S[0], CX[ia1], DPAc, SX[ia2]);
		COMPLEX_MUL4(ct2, S[1], CX[ia1], SPAc, CX[ia2]);
		COMPLEX_MUL4(ct3, S[2], SX[ia1], SPA,  SX[ia2]);
		COMPLEX_MUL4(ct4, S[3], SX[ia1], DPA,  CX[ia2]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;    /* RxxR wrt Ox2 */
		gradI = DFDP.imag;    /* RxxI wrt Ox2 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Oy2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt Ey2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	}  /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt IPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =   (MC1 + MC4) * gX[ia1] * gX[ia2] */
		COMPLEX_ADD2(ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;    /* RxxR wrt ipol */
		gradI = DFDP.imag;    /* RxxI wrt ipol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt QPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (MC2 + MC3) * gX[ia1] * gX[ia2] */
		COMPLEX_ADD2(ct1, MC2, MC3);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;    /* RxxR wrt qpol */
		gradI = DFDP.imag;    /* RxxI wrt qpol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* i (MC2 - MC3) * gX[ia1] * gX[ia2] */
		COMPLEX_SUB(ct1, MC2, MC3);
		COMPLEX_MUL3(DFDP, Jp, ct1, ggPD);
 		gradR = DFDP.real;    /* RxxR wrt upol */
		gradI = DFDP.imag;    /* RxxI wrt upol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt VPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (MC1 - MC4) * gX[ia1] * gX[ia2] */
		COMPLEX_SUB(ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;    /* RxxR wrt vpol */
		gradI = DFDP.imag;    /* RxxI wrt vpol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */

	/* gradient wrt PD = 0 */
	if (args->doFitRL) {
	  gradR = gradI = 0.0;
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	}

	/* Loop over antenna gains */
	for (k=0; k<2; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt gX1 */
	    /* Fitting? */
	    if (antGainFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gX[ia2] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia2*2+0], 0);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RxxR wrt gX1 */
		gradI = DFDP.imag;    /* RxxI wrt gX1 */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia1][0];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  case 1:  /* wrt gX2 */
	    if (antGainFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gX[ia1] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia1*2+0], 0);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RxxR wrt gX2 */
		gradI = DFDP.imag;    /* RxxI wrt gX2 */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia2][0];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  default:
	    break;
	  }; /* end gain switch */
	} /* end loop over gains */
	break;  /* End XX */
	
      case 1:     /* YY */
	if (wt[idata*4+kk]>0.0) {
	  isigma = wt[idata*4+kk];
	  /* VYY = {S[0] * SY[ia1] * SYc[ia2] * DPAc +       
	            S[1] * SY[ia1] * CYc[ia2] * SPAc +
		    S[2] * CY[ia1] * SYc[ia2] * SPA  + 
		    S[3] * CY[ia1] * CYc[ia2] * DPA} * g1Y * g2Y ;
	  */
	  COMPLEX_MUL3 (MC1, SY[ia1], SYc[ia2], DPAc);
	  COMPLEX_MUL3 (MC2, SY[ia1], CYc[ia2], SPAc);
	  COMPLEX_MUL3 (MC3, CY[ia1], SYc[ia2], SPA);
	  COMPLEX_MUL3 (MC4, CY[ia1], CYc[ia2], DPA);
	  COMPLEX_MUL2 (SM1, S[0], MC1);
	  COMPLEX_MUL2 (SM2, S[1], MC2);
	  COMPLEX_MUL2 (SM3, S[2], MC3);
	  COMPLEX_MUL2 (SM4, S[3], MC4);
	  COMPLEX_ADD4 (VYY, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ggPD,  antGain[ia1*2+1]*antGain[ia2*2+1], 0);
	  COMPLEX_MUL2 (VYY, VYY, ggPD);
	  modelR = VYY.real; modelI = VYY.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  /* residR = -residR; residI = - residI;   DEBUG */
	} else  residR = residI = 0.0; /* Invalid data */
	  
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = {(0, -1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		           (0,  1) * S[2] * MC2 + (0,  1) * S[3] * MC4} * gY[ia1] * gY[ia2] */
		COMPLEX_ADD2 (ct1, SM1, SM2);
		COMPLEX_MUL2 (ct2, Jm, ct1);
		COMPLEX_ADD2 (ct1, SM3, SM4);
		COMPLEX_MUL2 (ct3, Jp, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP, ct2, ggPD);
		gradR = DFDP.real;    /* RyyR wrt Oy1 */
		gradI = DFDP.imag;    /* RyyI wrt Oy1 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) * {S[0] * SYc[ia2] * DPAc * CYc[ia1] +
		                     S[1] * CYc[ia2] * SPAc * CYc[ia1] +
				     S[2] * SYc[ia2] * SPA  * SYc[ia1] +	       
				     S[3] * CYc[ia2] * DPA  * SYc[ia1]}  * gX[ia1] * gX[ia2]  */
		COMPLEX_MUL4(ct1, S[0], SYc[ia2], DPAc, CYc[ia1]);
		COMPLEX_MUL4(ct2, S[1], CYc[ia2], SPAc, CYc[ia1]);
		COMPLEX_MUL4(ct3, S[2], SYc[ia2], SPA,  SYc[ia1]);
		COMPLEX_MUL4(ct4, S[3], CYc[ia2], DPA,  SYc[ia1]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;    /* RyyR wrt Ey1 */
		gradI = DFDP.imag;    /* RyyI wrt Ey1 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* Ryy wrt Ex2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* Rll wrt Oy2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		           (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4} * gY[ia1] * gY[ia2]  */
		COMPLEX_ADD2 (ct1, SM1, SM3);
		COMPLEX_MUL2 (ct2, Jp, ct1);
		COMPLEX_ADD2 (ct1, SM2, SM4);
		COMPLEX_MUL2 (ct3, Jm, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP,  ct2, ggPD);
		gradR = DFDP.real;    /* RyyR wrt Oy2 */
		gradI = DFDP.imag;    /* RyyI wrt Oy2 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El2  */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) * {S[0] * SY[ia1] * DPAc * CY[ia2] + 
		                     S[1] * SY[ia1] * SPAc * SY[ia2] +
				     S[2] * CY[ia1] * SPA  * CY[ia2]  + 
				     S[3] * CY[ia1] * DPA  * SY[ia2]} * gY[ia1] * gY[ia2] */
		COMPLEX_MUL4(ct1, S[0], SY[ia1], DPAc, CY[ia2]);
		COMPLEX_MUL4(ct2, S[1], SY[ia1], SPAc, SY[ia2]);
		COMPLEX_MUL4(ct3, S[2], CY[ia1], SPA,  CY[ia2]);
		COMPLEX_MUL4(ct4, S[3], CY[ia1], DPA,  SY[ia2]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;    /* RyyR wrt Ey2 */
		gradI = DFDP.imag;    /* RyyI wrt Ey2 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	} /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt IPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =   (MC1 + MC4) * gY[ia1] * gY[ia2] */
		COMPLEX_ADD2(ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;    /* RyyR wrt ipol */
		gradI = DFDP.imag;    /* RyyI wrt ipol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Qpol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (MC2 + MC3) * gY[ia1] * gY[ia2] */
		COMPLEX_ADD2(ct1, MC2, MC3);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;    /* RyyR wrt qpol */
		gradI = DFDP.imag;    /* RyyI wrt qpol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* i (MC2 - MC3) * gY[ia1] * gY[ia2] */
		COMPLEX_SUB(ct1, MC2, MC3);
		COMPLEX_MUL3(DFDP, Jp, ct1, ggPD);
		gradR = DFDP.real;    /* RyyR wrt upol */
		gradI = DFDP.imag;    /* RyyI wrt upol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt VPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (MC1 - MC4) * gY[ia1] * gY[ia2] */
		COMPLEX_SUB (ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;    /* RyyR wrt vpol */
		gradI = DFDP.imag;    /* RyyI wrt vpol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */

	/* gradient wrt PD - no effect */
	if (args->doFitRL) {
	  gradR = gradI = 0.0;
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	}

	/* Loop over antenna gains */
	for (k=0; k<2; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt gY1 */
	    /* Fitting? */
	    if (antGainFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gY[ia2] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia2*2+1], 0);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RyyR wrt gY1 */
		gradI = DFDP.imag;    /* RyyI wrt gY1 */
	 	/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia1][1];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  case 1:  /* wrt gY2 */
	    if (antGainFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gY[ia1] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia1*2+0], 1);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RyyR wrt gY2 */
		gradI = DFDP.imag;    /* RyyI wrt gY2 */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia2][1];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  default:
	    break;
	  }; /* end gain switch */
	} /* end loop over gains */
	
	break;  /* End YY */
	
      case 2:     /* XY */
	if (wt[idata*4+kk]>0.0) {
	  isigma = wt[idata*4+kk];
	  /* VXY = {S[0] * CX[ia1] * SYc[ia2] * DPAc +       
	            S[1] * CX[ia1] * CYc[ia2] * SPAc +
		    S[2] * SX[ia1] * SYc[ia2] * SPA  + 
		    S[3] * SX[ia1] * CYc[ia2] * DPA}} * g1X * g2Y * exp(i PD);
	  */
	  COMPLEX_MUL3 (MC1, CX[ia1], SYc[ia2], DPAc);
	  COMPLEX_MUL3 (MC2, CX[ia1], CYc[ia2], SPAc);
	  COMPLEX_MUL3 (MC3, SX[ia1], SYc[ia2], SPA);
	  COMPLEX_MUL3 (MC4, SX[ia1], CYc[ia2], DPA);
	  COMPLEX_MUL2 (SM1, S[0], MC1);
	  COMPLEX_MUL2 (SM2, S[1], MC2);
	  COMPLEX_MUL2 (SM3, S[2], MC3);
	  COMPLEX_MUL2 (SM4, S[3], MC4);
	  COMPLEX_ADD4 (VXY, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct1,  antGain[ia1*2+0]*antGain[ia2*2+1], 0);
	  COMPLEX_EXP (ct2, PD);
	  COMPLEX_MUL2 (ggPD, ct1, ct2);
	  COMPLEX_MUL2 (VXY, VXY, ggPD);
	  modelR = VXY.real; modelI = VXY.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  /* residR = -residR; residI = - residI;   DEBUG */
	} else  residR = residI = 0.0; /* Invalid data */
	
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Or1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = {(0, -1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		           (0,  1) * S[2] * MC3 + (0,  1) * S[3] * MC4} *
			   gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_ADD2 (ct1, SM1, SM2);
		COMPLEX_MUL2 (ct2, Jm, ct1);
		COMPLEX_ADD2 (ct1, SM3, SM4);
		COMPLEX_MUL2 (ct3, Jp, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP, ct2, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =  (0, -1) * {S[0] * SYc[ia2] * DPAc * SXc[ia1]  + 
		                      S[1] * CYc[ia2] * SPAc * SXc[ia1]  +
				      S[2] * SYc[ia2] * SPA  * CXc[ia1]  + 
				      S[3] * CYc[ia2] * DPA  * CXc[ia1]  } * 
				              gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_MUL4(ct1, S[0], SYc[ia2], DPAc, SXc[ia1]);
		COMPLEX_MUL4(ct2, S[1], CYc[ia2], SPAc, SXc[ia1]);
		COMPLEX_MUL4(ct3, S[2], SYc[ia2], SPA,  CXc[ia1]);
		COMPLEX_MUL4(ct4, S[3], CYc[ia2], DPA,  CXc[ia1]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Ol1 =0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt El1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Oy2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		           (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4} * 
                                 gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_ADD2 (ct1, SM1, SM3);
		COMPLEX_MUL2 (ct2, Jp, ct1);
		COMPLEX_ADD2 (ct1, SM2, SM4);
		COMPLEX_MUL2 (ct3, Jm, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP,  ct2, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt Ey2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) *{S[0] * CX[ia1] * DPAc * CY[ia2] +       
		                    S[1] * CX[ia1] * SPAc * SY[ia2] +
				    S[2] * SX[ia1] * SPA  * CY[ia2] + 
				    S[3] * SX[ia1] * DPA  * SY[ia2]} * 
			   	          gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_MUL4(ct1, S[0], CX[ia1], DPAc, CY[ia2]);
		COMPLEX_MUL4(ct2, S[1], CX[ia1], SPAc, SY[ia2]);
		COMPLEX_MUL4(ct3, S[2], SX[ia1], SPA,  CY[ia2]);
		COMPLEX_MUL4(ct4, S[3], SX[ia1], DPA,  SY[ia2]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	  } /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt IPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =   (MC1 + MC4) * gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_ADD2(ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt QPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (MC2 + MC3) * gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_ADD2(ct1, MC2, MC3);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* i (MC2 - MC3) * gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_SUB(ct1, MC2, MC3);
		COMPLEX_MUL3(DFDP, Jp, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt VPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (MC1 - MC4) * gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_SUB (ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	   } 
	    break;
	  default:
	    break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */
	
	/* gradient wrt PD */
	if (args->doFitRL) {
	  if (wt[idata*4+kk]>0.0) {
	    /* part = (0,  1) * VXY */   
	    COMPLEX_MUL2 (DFDP, Jp, VXY);
	    gradR = DFDP.real;
	    gradI = DFDP.imag;
	    /* gradR = -gradR; gradI = -gradI;   DEBUG */
	  } else gradR = gradI = 0.0;    /* invalid data */
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	}
	
	/* Loop over antenna gains */
	for (k=0; k<2; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt gX1 */
	    /* Fitting? */
	    if (antGainFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gY[ia2] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia2*2+1], 0);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RxyR wrt gX1 */
		gradI = DFDP.imag;    /* RxyI wrt gX1 */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia1][0];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  case 1:  /* wrt gY2 */
	    if (antGainFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gX[ia1] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia1*2+0], 1);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RxyR wrt gY2 */
		gradI = DFDP.imag;    /* RxyI wrt gY2 */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia2][1];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  default:
	    break;
	  }; /* end gain switch */
	} /* end loop over gains */
	
	break;  /* End XY */
      	
      case 3:     /* YX */
	if (wt[idata*4+kk]>0.0) {
	  isigma = wt[idata*4+kk];
	  /* VYX = {S[0] * SY[ia1] * CXc[ia2] * DPAc +       
	            S[1] * SY[ia1] * SXc[ia2] * SPAc +
		    S[2] * CY[ia1] * CXc[ia2] * SPA  + 
		    S[3] * CY[ia1] * SXc[ia2] * DPA} * g1Y * g2X * exp(-i PD)
	  */
	  COMPLEX_MUL3 (MC1, SY[ia1], CXc[ia2], DPAc);
	  COMPLEX_MUL3 (MC2, SY[ia1], SXc[ia2], SPAc);
	  COMPLEX_MUL3 (MC3, CY[ia1], CXc[ia2], SPA);
	  COMPLEX_MUL3 (MC4, CY[ia1], SXc[ia2], DPA);
	  COMPLEX_MUL2 (SM1, S[0], MC1);
	  COMPLEX_MUL2 (SM2, S[1], MC2);
	  COMPLEX_MUL2 (SM3, S[2], MC3);
	  COMPLEX_MUL2 (SM4, S[3], MC4);
	  COMPLEX_ADD4 (VYX, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct1,  antGain[ia1*2+1]*antGain[ia2*2+0], 0);
	  COMPLEX_EXP (ct2, -PD);
	  COMPLEX_MUL2 (ggPD, ct1, ct2);
	  COMPLEX_MUL2 (VYX, VYX, ggPD);
	  modelR = VYX.real; modelI = VYX.imag;
	} else  residR = residI = 0.0; /* Invalid data */
	
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Oy1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		           (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4} * 
			        gY[ia1] * gX[ia2]  * exp(-i PD) */
		COMPLEX_ADD2 (ct1, SM1, SM3);
		COMPLEX_MUL2 (ct2, Jp, ct1);
		COMPLEX_ADD2 (ct1, SM2, SM4);
		COMPLEX_MUL2 (ct3, Jm, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP,  ct2, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt Ey1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) * {S[0] * CXc[ia2] * DPAc * CYc[ia1] +       
		                     S[1] * SXc[ia2] * SPAc * CYc[ia1] +
				     S[2] * CXc[ia2] * SPA  * SYc[ia1]  + 
				     S[3] * SXc[ia2] * DPA  * SYc[ia1]} * 
				           gY[ia1] * gX[ia2] * exp(-i PD) */
		COMPLEX_MUL4(ct1, S[0], CXc[ia2], DPAc, CYc[ia1]);
		COMPLEX_MUL4(ct2, S[1], SXc[ia2], SPAc, CYc[ia1]);
		COMPLEX_MUL4(ct3, S[2], CXc[ia2], SPA,  SYc[ia1]);
		COMPLEX_MUL4(ct4, S[3], SXc[ia2], DPA,  SYc[ia1]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		           (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4} * 
			         gY[ia1] * gX[ia2]  * exp(-i PD) */
		COMPLEX_ADD2 (ct1, SM1, SM3);
		COMPLEX_MUL2 (ct2, Jp, ct1);
		COMPLEX_ADD2 (ct1, SM2, SM4);
		COMPLEX_MUL2 (ct3, Jm, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP,  ct2, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) * {S[0] * CXc[ia2] * DPAc * CYc[ia1] +       
		                     S[1] * SXc[ia2] * SPAc * CYc[ia1] +
				     S[2] * CXc[ia2] * SPA  * SYc[ia1]  + 
				     S[3] * SXc[ia2] * DPA  * SYc[ia1]} * 
				          gY[ia1] * gX[ia2] * exp(-i PD) */
		COMPLEX_MUL4(ct1, S[0], CXc[ia2], DPAc, CYc[ia1]);
		COMPLEX_MUL4(ct2, S[1], SXc[ia2], SPAc, CYc[ia1]);
		COMPLEX_MUL4(ct3, S[2], CXc[ia2], SPA,  SYc[ia1]);
		COMPLEX_MUL4(ct4, S[3], SXc[ia2], DPA,  SYc[ia1]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Oy2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt Ey2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	} /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt IPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =   (MC1 + MC4) * gY[ia1] * gX[ia2] */
		COMPLEX_ADD2(ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt QPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (MC2 + MC3) * gY[ia1] * gX[ia2] */
		COMPLEX_ADD2(ct1, MC2, MC3);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* i (MC2 - MC3) * gY[ia1] * gX[ia2] */
		COMPLEX_SUB(ct1, MC2, MC3);
		COMPLEX_MUL3(DFDP, Jp, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt VPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (MC1 - MC4) * gY[ia1] * gX[ia2] */
		COMPLEX_SUB (ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	break;
      default:
	break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */
	
	/* gradient wrt PD */
	if (args->doFitRL) {
	  if (wt[idata*4+kk]>0.0) {
	    /* part = (0, -1) * VYX */   
	    COMPLEX_MUL2 (DFDP, Jm, VYX);
	    gradR = DFDP.real;
	    gradI = DFDP.imag;
	    /* gradR = -gradR; gradI = -gradI;   DEBUG */
	  } else gradR = gradI = 0.0;    /* invalid data */
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	}
	/* Loop over antenna gains */
	for (k=0; k<2; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt gY1 */
	    /* Fitting? */
	    if (antGainFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gX[ia2] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia2*2+0], 0);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RyxR wrt gY1 */
		gradI = DFDP.imag;    /* RyxI wrt gY1 */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia1][0];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  case 1:  /* wrt gX2 */
	    if (antGainFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gY[ia1] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia1*2+1], 1);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RyxR wrt gX2 */
		gradI = DFDP.imag;    /* RyxI wrt gX2 */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia2][1];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  default:
	    break;
	  }; /* end gain switch */
	} /* end loop over gains */
	/* End YX */	
	
      default:
	  break;
	}; /* end switch over data correlation */
      
      i++;  /* Update complex datum number */
    } /* end loop over correlations */
  } /* End loop over visibilities */

  return GSL_SUCCESS;
} /*  end PolnFitJacOEXY */

/**
 * Linear feed Function & Jacobian evaluator for polarization fitting solver
 * Orientation/Ellipticity for linear feeds version
 * Evaluates partial derivatives of model wrt each parameter
 * \param x       Vector of parameters to be fitted
 *                Flux,array_of polarization_terms
 * \param param   Function parameter structure (ObitPolnCalFit)
 * \param f       Vector of (model-obs)/sigma for data points
 * \param J       Jacobian matrix J[data_point, parameter]
 * \return completion code GSL_SUCCESS=OK
 */
static int PolnFitFuncJacOEXY (const gsl_vector *x, void *params, 
			       gsl_vector *f, gsl_matrix *J)
{
  ObitPolnCalFit *args = (ObitPolnCalFit*)params;
  ofloat *data, *wt;
  gboolean   **antFit    = args->antFit;
  olong      **antPNumb  = args->antPNumb;
  odouble     *antParm   = args->antParm;
  odouble    *antGain    = args->antGain;
  gboolean   **antGainFit= args->antGainFit;
  olong     **antGainPNumb  = args->antGainPNumb;
  gboolean   **souFit    = args->souFit;
  odouble    *souParm    = args->souParm;
  olong      **souPNumb  = args->souPNumb;
  dcomplex   *CX     = args->RS;
  dcomplex   *SX     = args->RD;
  dcomplex   *CY     = args->LS;
  dcomplex   *SY     = args->LD;
  dcomplex   *CXc    = args->RSc;
  dcomplex   *SXc    = args->RDc;
  dcomplex   *CYc    = args->LSc;
  dcomplex   *SYc    = args->LDc;
  ofloat PD, chi1, chi2;
  double val;
  odouble ipol=0.0, qpol=0.0, upol=0.0, vpol=0.0;
  odouble residR=0.0, residI=0.0, gradR=0.0, gradI=0.0, modelR=0.0, modelI=0.0, isigma=0.0;
  olong k, kk, iant, ia1, ia2, isou, idata, refAnt;
  olong isouLast=-999;
  dcomplex  SPA, DPA, SPAc, DPAc, ggPD;
  dcomplex ct1, ct2, ct3, ct4, ct5, Jm, Jp;
  dcomplex S[4], VXX, VXY, VYX, VYY, MC1, MC2, MC3, MC4, DFDP;
  dcomplex SM1, SM2, SM3, SM4;
  size_t i, j;

   /* Initialize output */
  val = 0.0;
  for (i=0; i<args->ndata; i++) {
    gsl_vector_set(f, i, val);
    for (j=0; j<args->nparam; j++) 
      gsl_matrix_set(J, i, j, val);
  }

  COMPLEX_SET (S[0], 0.0, 0.0);  /* Initialize poln vector */
  COMPLEX_SET (S[1], 0.0, 0.0);
  COMPLEX_SET (S[2], 0.0, 0.0);
  COMPLEX_SET (S[3], 0.0, 0.0);
  COMPLEX_SET (MC1, 0.0, 0.0);  /* Other stuff */
  COMPLEX_SET (MC2, 0.0, 0.0);
  COMPLEX_SET (MC3, 0.0, 0.0);
  COMPLEX_SET (MC4, 0.0, 0.0);
  COMPLEX_SET (VXX, 0.0, 0.0);
  COMPLEX_SET (VYY, 0.0, 0.0);
  COMPLEX_SET (VYX, 0.0, 0.0);
  COMPLEX_SET (VXY, 0.0, 0.0);
  COMPLEX_SET (Jm,  0.0,-1.0);
  COMPLEX_SET (Jp,  0.0, 1.0);
  COMPLEX_SET (SM1, 0.0, 0.0);
  COMPLEX_SET (SM2, 0.0, 0.0);
  COMPLEX_SET (SM3, 0.0, 0.0);
  COMPLEX_SET (SM4, 0.0, 0.0);
  COMPLEX_SET (ggPD, 0.0, 0.0);
  
  /* R-L phase difference */
  if (args->doFitRL) {
    j = args->PDPNumb;
    PD = gsl_vector_get(x, j);
  } else PD = args->PD;
  
  /* get model parameters - first antenna */
  for (iant=0; iant<args->nant; iant++) {
    /* Loop over antenna parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if ((antFit[iant][k]) && (args->gotAnt[iant])) {
	j = antPNumb[iant][k];
	antParm[iant*4+k] = gsl_vector_get(x, j);
      }
    } /* end loop over antenna parameters */
  } /* end loop over antennas */
  
  /* antenna gain */
  for (iant=0; iant<args->nant; iant++) {
    /* Loop over antenna gains */
    for (k=0; k<2; k++) {
      /* Fitting? */
      if ((antGainFit[iant][k]) && (args->gotAnt[iant])) {
	j = antGainPNumb[iant][k];
	antGain[iant*2+k] = gsl_vector_get(x, j);
      }
    } /* end loop over antenna parameters */
  } /* end loop over antennas */
  
  /* Ref antenna - 0 rel */
  refAnt = MAX(0, args->refAnt-1);
  
  /* now source */
  for (isou=0; isou<args->nsou; isou++) {
    /* Loop over source parameters */
    for (k=0; k<4; k++) {
      /* Fitting? */
      if (souFit[isou][k]) {
	j = souPNumb[isou][k];
	souParm[isou*4+k] = gsl_vector_get(x, j);
      }
    } /* end loop over source parameters */
  } /* end loop over sources */
  
  /* data & wt pointers */
  data = args->inData;
  wt   = args->inWt;

   /* Injest model, factorize into antenna components - 
     data in order Orientation R/X, Elipticity R/X, Orientation L/Y, Elipticity L/Y */
  /* Elipticity, Orientation terms */
  for (i=0; i<args->nant; i++) {
    COMPLEX_EXP (ct1, -antParm[i*4+0]);
    COMPLEX_SET (ct2, cos(G_PI*0.25+antParm[i*4+1]), 0.0);
    COMPLEX_MUL3 (CX[i], Jp, ct1, ct2);
    COMPLEX_EXP (ct1, antParm[i*4+0]);
    COMPLEX_SET (ct2, sin(G_PI*0.25+antParm[i*4+1]), 0.0);
    COMPLEX_MUL2 (SX[i], ct1, ct2);
    COMPLEX_EXP (ct1, antParm[i*4+2]);
    COMPLEX_SET (ct2, cos(G_PI*0.25-antParm[i*4+3]), 0.0);
    COMPLEX_MUL2 (CY[i], ct1, ct2);
    COMPLEX_EXP (ct1, -antParm[i*4+2]);
    COMPLEX_SET (ct2, sin(G_PI*0.25-antParm[i*4+3]), 0.0);
    COMPLEX_MUL3 (SY[i], Jp, ct1, ct2);
    COMPLEX_CONJUGATE (CXc[i], CX[i]);
    COMPLEX_CONJUGATE (SXc[i], SX[i]);
    COMPLEX_CONJUGATE (CYc[i], CY[i]);
    COMPLEX_CONJUGATE (SYc[i], SY[i]);
  }

  /* Loop over data */
  i = 0;
  for (idata=0; idata<args->nvis; idata++) {
    /* Parallactic angle terms */
    chi1  = data[idata*10+0];   /* parallactic angle ant 1 */
    chi2  = data[idata*10+1];   /* parallactic angle ant 2 */
    COMPLEX_EXP (SPA,chi1+chi2);
    COMPLEX_EXP (DPA,chi1-chi2);
    COMPLEX_CONJUGATE (SPAc, SPA);
    COMPLEX_CONJUGATE (DPAc, DPA);

    isou  = MAX (0, args->souNo[idata]);    /* Source number */
    /* New source? get parameters */
    if (isou!=isouLast) {
      isouLast = isou;
      /* Source parameters */
      ipol = souParm[isou*4+0];
      /* Fitting or fixed? */
      if (args->souFit[isou][1]) 
	qpol = souParm[isou*4+1];
      else
	qpol = args->PPol[isou]*ipol*cos(args->RLPhase[isou]);
      if (args->souFit[isou][2]) 
	upol = souParm[isou*4+2];
      else
	upol = args->PPol[isou]*ipol*sin(args->RLPhase[isou]);
      vpol = souParm[isou*4+3];
      /* Complex Stokes array */
      COMPLEX_SET (S[0], ipol+vpol, 0.0);
      COMPLEX_SET (S[1], qpol,  upol);
      COMPLEX_SET (S[2], qpol, -upol);
      COMPLEX_SET (S[3], ipol-vpol, 0.0);
    }

    /* Antenna parameters */
    ia1    = args->antNo[idata*2+0];
    ia2    = args->antNo[idata*2+1]; 
    
    /* i = datum number */
    /* Loop over correlations calculating derivatives */
    for (kk=0; kk<4; kk++) {
      switch (kk) { 
	
      case 0:     /* XX */
	if (wt[idata*4+kk]>0.0) {
	  isigma = wt[idata*4+kk];
	  /* VXX = {S[0] * CX[ia1] * CXc[ia2] * DPAc  +
	            S[1] * CX[ia1] * SXc[ia2] * SPAc  +
		    S[2] * SX[ia1] * CXc[ia2] * SPA   + 
		    S[3] * SX[ia1] * SXc[ia2] * DPA} * g1X * g2X ;
	  */
	  COMPLEX_MUL3 (MC1, CX[ia1], CXc[ia2], DPAc);
	  COMPLEX_MUL3 (MC2, CX[ia1], SXc[ia2], SPAc);
	  COMPLEX_MUL3 (MC3, SX[ia1], CXc[ia2], SPA);
	  COMPLEX_MUL3 (MC4, SX[ia1], SXc[ia2], DPA);
	  COMPLEX_MUL2 (SM1, S[0], MC1);
	  COMPLEX_MUL2 (SM2, S[1], MC2);
	  COMPLEX_MUL2 (SM3, S[2], MC3);
	  COMPLEX_MUL2 (SM4, S[3], MC4);
	  COMPLEX_ADD4 (VXX, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ggPD,  antGain[ia1*2+0]*antGain[ia2*2+0], 0);
	  COMPLEX_MUL2 (VXX, VXX, ggPD);
	  modelR = VXX.real; modelI = VXX.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  /* residR = -residR; residI = - residI;   DEBUG */
	  gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
	  gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
	} else  residR = residI = 0.0; /* Invalid data */
	
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* {(0, -1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		    (0,  1) * S[2] * MC3 + (0,  1) * S[3] * MC4}  * gX[ia1] * gX[ia2]  */
		COMPLEX_ADD2 (ct1, SM1, SM2);
		COMPLEX_MUL2 (ct2, Jm, ct1);
		COMPLEX_ADD2 (ct1, SM3, SM4);
		COMPLEX_MUL2 (ct3, Jp, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP, ct2, ggPD);
		gradR = DFDP.real;    /* RxxR wrt Ox1 */
		gradI = DFDP.imag;    /* RxxI wrt Ox1 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI =0.0;     /* Invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex1 */
	    /* Fitting? */
 	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) * {S[0] * CXc[ia2] * DPAc * SXc[ia1] +
		                     S[1] * SXc[ia2] * SPAc * SXc[ia1] +
				     S[2] * CXc[ia2] * SPA  * CXc[ia1] +	       
				     S[3] * SXc[ia2] * DPA  * CXc[ia1]}  * gX[ia1] * gX[ia2]  */
		COMPLEX_MUL4(ct1, S[0], CXc[ia2], DPAc, SXc[ia1]);
		COMPLEX_MUL4(ct2, S[1], SXc[ia2], SPAc, SXc[ia1]);
		COMPLEX_MUL4(ct3, S[2], CXc[ia2], SPA,  SXc[ia1]);
		COMPLEX_MUL4(ct4, S[3], SXc[ia2], DPA,  CXc[ia1]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;    /* RxxR wrt Ex1 */
		gradI = DFDP.imag;    /* RxxI wrt Ex1 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt OY1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt EY1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		    (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4}  * gX[ia1] * gX[ia2]   */
		COMPLEX_ADD2 (ct1, SM1, SM3);
		COMPLEX_MUL2 (ct2, Jp, ct1);
		COMPLEX_ADD2 (ct1, SM2, SM4);
		COMPLEX_MUL2 (ct3, Jm, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP, ct2, ggPD);
		gradR = DFDP.real;    /* RxxR wrt Ox2 */
		gradI = DFDP.imag;    /* RxxI wrt Ox2 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) * {S[0] * CX[ia1] * DPAc * SX[ia2] +
		                     S[1] * CX[ia1] * SPAc * CX[ia2] +
				     S[2] * SX[ia1] * SPA  * SX[ia2] +
				     S[3] * SX[ia1] * DPA  * CX[ia2]} * gX[ia1] * gX[ia2]  */
		COMPLEX_MUL4(ct1, S[0], CX[ia1], DPAc, SX[ia2]);
		COMPLEX_MUL4(ct2, S[1], CX[ia1], SPAc, CX[ia2]);
		COMPLEX_MUL4(ct3, S[2], SX[ia1], SPA,  SX[ia2]);
		COMPLEX_MUL4(ct4, S[3], SX[ia1], DPA,  CX[ia2]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;    /* RxxR wrt Ox2 */
		gradI = DFDP.imag;    /* RxxI wrt Ox2 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Oy2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt Ey2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	}  /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt IPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =   (MC1 + MC4) * gX[ia1] * gX[ia2] */
		COMPLEX_ADD2(ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;    /* RxxR wrt ipol */
		gradI = DFDP.imag;    /* RxxI wrt ipol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt QPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (MC2 + MC3) * gX[ia1] * gX[ia2] */
		COMPLEX_ADD2(ct1, MC2, MC3);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;    /* RxxR wrt qpol */
		gradI = DFDP.imag;    /* RxxI wrt qpol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* i (MC2 - MC3) * gX[ia1] * gX[ia2] */
		COMPLEX_SUB(ct1, MC2, MC3);
		COMPLEX_MUL3(DFDP, Jp, ct1, ggPD);
 		gradR = DFDP.real;    /* RxxR wrt upol */
		gradI = DFDP.imag;    /* RxxI wrt upol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt VPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (MC1 - MC4) * gX[ia1] * gX[ia2] */
		COMPLEX_SUB(ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;    /* RxxR wrt vpol */
		gradI = DFDP.imag;    /* RxxI wrt vpol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */

	/* gradient wrt PD = 0 */
	if (args->doFitRL) {
	  gradR = gradI = 0.0;
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	}

	/* Loop over antenna gains */
	for (k=0; k<2; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt gX1 */
	    /* Fitting? */
	    if (antGainFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gX[ia2] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia2*2+0], 0);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RxxR wrt gX1 */
		gradI = DFDP.imag;    /* RxxI wrt gX1 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
		 /* gradR = 0.0; gradI = 0.0; DEBUG GAIN */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia1][0];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  case 1:  /* wrt gX2 */
	    if (antGainFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gX[ia1] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia1*2+0], 0);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RxxR wrt gX2 */
		gradI = DFDP.imag;    /* RxxI wrt gX2 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia2][0];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  default:
	    break;
	  }; /* end gain switch */
	} /* end loop over gains */
	break;  /* End XX */
	
      case 1:     /* YY */
	if (wt[idata*4+kk]>0.0) {
	  isigma = wt[idata*4+kk];
	  /* VYY = {S[0] * SY[ia1] * SYc[ia2] * DPAc +       
	            S[1] * SY[ia1] * CYc[ia2] * SPAc +
		    S[2] * CY[ia1] * SYc[ia2] * SPA  + 
		    S[3] * CY[ia1] * CYc[ia2] * DPA} * g1Y * g2Y ;
	  */
	  COMPLEX_MUL3 (MC1, SY[ia1], SYc[ia2], DPAc);
	  COMPLEX_MUL3 (MC2, SY[ia1], CYc[ia2], SPAc);
	  COMPLEX_MUL3 (MC3, CY[ia1], SYc[ia2], SPA);
	  COMPLEX_MUL3 (MC4, CY[ia1], CYc[ia2], DPA);
	  COMPLEX_MUL2 (SM1, S[0], MC1);
	  COMPLEX_MUL2 (SM2, S[1], MC2);
	  COMPLEX_MUL2 (SM3, S[2], MC3);
	  COMPLEX_MUL2 (SM4, S[3], MC4);
	  COMPLEX_ADD4 (VYY, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ggPD,  antGain[ia1*2+1]*antGain[ia2*2+1], 0);
	  COMPLEX_MUL2 (VYY, VYY, ggPD);
	  modelR = VYY.real; modelI = VYY.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  /* residR = -residR; residI = - residI;   DEBUG */
	  gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
	  gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
	} else  residR = residI = 0.0; /* Invalid data */
	  
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Oy1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = {(0, -1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		           (0,  1) * S[2] * MC2 + (0,  1) * S[3] * MC4} * gY[ia1] * gY[ia2] */
		COMPLEX_ADD2 (ct1, SM1, SM2);
		COMPLEX_MUL2 (ct2, Jm, ct1);
		COMPLEX_ADD2 (ct1, SM3, SM4);
		COMPLEX_MUL2 (ct3, Jp, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP, ct2, ggPD);
		gradR = DFDP.real;    /* RyyR wrt Oy1 */
		gradI = DFDP.imag;    /* RyyI wrt Oy1 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt Ey1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) * {S[0] * SYc[ia2] * DPAc * CYc[ia1] +
		                     S[1] * CYc[ia2] * SPAc * CYc[ia1] +
				     S[2] * SYc[ia2] * SPA  * SYc[ia1] +	       
				     S[3] * CYc[ia2] * DPA  * SYc[ia1]}  * gX[ia1] * gX[ia2]  */
		COMPLEX_MUL4(ct1, S[0], SYc[ia2], DPAc, CYc[ia1]);
		COMPLEX_MUL4(ct2, S[1], CYc[ia2], SPAc, CYc[ia1]);
		COMPLEX_MUL4(ct3, S[2], SYc[ia2], SPA,  SYc[ia1]);
		COMPLEX_MUL4(ct4, S[3], CYc[ia2], DPA,  SYc[ia1]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;    /* RyyR wrt Ey1 */
		gradI = DFDP.imag;    /* RyyI wrt Ey1 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* Ryy wrt Ex2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* Rll wrt Oy2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		           (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4} * gY[ia1] * gY[ia2]  */
		COMPLEX_ADD2 (ct1, SM1, SM3);
		COMPLEX_MUL2 (ct2, Jp, ct1);
		COMPLEX_ADD2 (ct1, SM2, SM4);
		COMPLEX_MUL2 (ct3, Jm, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP,  ct2, ggPD);
		gradR = DFDP.real;    /* RyyR wrt Oy2 */
		gradI = DFDP.imag;    /* RyyI wrt Oy2 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt Ey2  */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) * {S[0] * SY[ia1] * DPAc * CY[ia2] + 
		                     S[1] * SY[ia1] * SPAc * SY[ia2] +
				     S[2] * CY[ia1] * SPA  * CY[ia2]  + 
				     S[3] * CY[ia1] * DPA  * SY[ia2]} * gY[ia1] * gY[ia2] */
		COMPLEX_MUL4(ct1, S[0], SY[ia1], DPAc, CY[ia2]);
		COMPLEX_MUL4(ct2, S[1], SY[ia1], SPAc, SY[ia2]);
		COMPLEX_MUL4(ct3, S[2], CY[ia1], SPA,  CY[ia2]);
		COMPLEX_MUL4(ct4, S[3], CY[ia1], DPA,  SY[ia2]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;    /* RyyR wrt Ey2 */
		gradI = DFDP.imag;    /* RyyI wrt Ey2 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	} /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt IPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =   (MC1 + MC4) * gY[ia1] * gY[ia2] */
		COMPLEX_ADD2(ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;    /* RyyR wrt ipol */
		gradI = DFDP.imag;    /* RyyI wrt ipol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Qpol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (MC2 + MC3) * gY[ia1] * gY[ia2] */
		COMPLEX_ADD2(ct1, MC2, MC3);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;    /* RyyR wrt qpol */
		gradI = DFDP.imag;    /* RyyI wrt qpol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* i (MC2 - MC3) * gY[ia1] * gY[ia2] */
		COMPLEX_SUB(ct1, MC2, MC3);
		COMPLEX_MUL3(DFDP, Jp, ct1, ggPD);
		gradR = DFDP.real;    /* RyyR wrt upol */
		gradI = DFDP.imag;    /* RyyI wrt upol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt VPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (MC1 - MC4) * gY[ia1] * gY[ia2] */
		COMPLEX_SUB (ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;    /* RyyR wrt vpol */
		gradI = DFDP.imag;    /* RyyI wrt vpol */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */

	/* gradient wrt PD - no effect */
	if (args->doFitRL) {
	  gradR = gradI = 0.0;
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	}

	/* Loop over antenna gains */
	for (k=0; k<2; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt gY1 */
	    /* Fitting? */
	    if (antGainFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gY[ia2] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia2*2+1], 0);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RyyR wrt gY1 */
		gradI = DFDP.imag;    /* RyyI wrt gY1 */
		/* gradR = 0.0; gradI = 0.0; DEBUG GAIN */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia1][1];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  case 1:  /* wrt gY2 */
	    if (antGainFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gY[ia1] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia1*2+0], 1);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RyyR wrt gY2 */
		gradI = DFDP.imag;    /* RyyI wrt gY2 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia2][1];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  default:
	    break;
	  }; /* end gain switch */
	} /* end loop over gains */
	
	break;  /* End YY */
	
      case 2:     /* XY */
	if (wt[idata*4+kk]>0.0) {
	  isigma = wt[idata*4+kk];
	  /* VXY = {S[0] * CX[ia1] * SYc[ia2] * DPAc +       
	            S[1] * CX[ia1] * CYc[ia2] * SPAc +
		    S[2] * SX[ia1] * SYc[ia2] * SPA  + 
		    S[3] * SX[ia1] * CYc[ia2] * DPA}} * g1X * g2Y * exp(i PD);
	  */
	  COMPLEX_MUL3 (MC1, CX[ia1], SYc[ia2], DPAc);
	  COMPLEX_MUL3 (MC2, CX[ia1], CYc[ia2], SPAc);
	  COMPLEX_MUL3 (MC3, SX[ia1], SYc[ia2], SPA);
	  COMPLEX_MUL3 (MC4, SX[ia1], CYc[ia2], DPA);
	  COMPLEX_MUL2 (SM1, S[0], MC1);
	  COMPLEX_MUL2 (SM2, S[1], MC2);
	  COMPLEX_MUL2 (SM3, S[2], MC3);
	  COMPLEX_MUL2 (SM4, S[3], MC4);
	  COMPLEX_ADD4 (VXY, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct1,  antGain[ia1*2+0]*antGain[ia2*2+1], 0);
	  COMPLEX_EXP (ct2, PD);
	  COMPLEX_MUL2 (ggPD, ct1, ct2);
	  COMPLEX_MUL2 (VXY, VXY, ggPD);
	  modelR = VXY.real; modelI = VXY.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  /* residR = -residR; residI = - residI;   DEBUG */
	  gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
	  gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
	} else  residR = residI = 0.0; /* Invalid data */
	
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = {(0, -1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		           (0,  1) * S[2] * MC3 + (0,  1) * S[3] * MC4} *
			   gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_ADD2 (ct1, SM1, SM2);
		COMPLEX_MUL2 (ct2, Jm, ct1);
		COMPLEX_ADD2 (ct1, SM3, SM4);
		COMPLEX_MUL2 (ct3, Jp, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP, ct2, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =  (0, -1) * {S[0] * SYc[ia2] * DPAc * SXc[ia1]  + 
		                      S[1] * CYc[ia2] * SPAc * SXc[ia1]  +
				      S[2] * SYc[ia2] * SPA  * CXc[ia1]  + 
				      S[3] * CYc[ia2] * DPA  * CXc[ia1]  } * 
				              gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_MUL4(ct1, S[0], SYc[ia2], DPAc, SXc[ia1]);
		COMPLEX_MUL4(ct2, S[1], CYc[ia2], SPAc, SXc[ia1]);
		COMPLEX_MUL4(ct3, S[2], SYc[ia2], SPA,  CXc[ia1]);
		COMPLEX_MUL4(ct4, S[3], CYc[ia2], DPA,  CXc[ia1]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Oy1 =0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt Ey1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Oy2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		           (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4} * 
                                 gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_ADD2 (ct1, SM1, SM3);
		COMPLEX_MUL2 (ct2, Jp, ct1);
		COMPLEX_ADD2 (ct1, SM2, SM4);
		COMPLEX_MUL2 (ct3, Jm, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP,  ct2, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt Ey2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) *{S[0] * CX[ia1] * DPAc * CY[ia2] +       
		                    S[1] * CX[ia1] * SPAc * SY[ia2] +
				    S[2] * SX[ia1] * SPA  * CY[ia2] + 
				    S[3] * SX[ia1] * DPA  * SY[ia2]} * 
			   	          gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_MUL4(ct1, S[0], CX[ia1], DPAc, CY[ia2]);
		COMPLEX_MUL4(ct2, S[1], CX[ia1], SPAc, SY[ia2]);
		COMPLEX_MUL4(ct3, S[2], SX[ia1], SPA,  CY[ia2]);
		COMPLEX_MUL4(ct4, S[3], SX[ia1], DPA,  SY[ia2]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	  } /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt IPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =   (MC1 + MC4) * gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_ADD2(ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt QPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (MC2 + MC3) * gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_ADD2(ct1, MC2, MC3);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* i (MC2 - MC3) * gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_SUB(ct1, MC2, MC3);
		COMPLEX_MUL3(DFDP, Jp, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt VPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (MC1 - MC4) * gX[ia1] * gY[ia2] * exp(i PD) */
		COMPLEX_SUB (ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	   } 
	    break;
	  default:
	    break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */
	
	/* gradient wrt PD */
	if (args->doFitRL) {
	  if (wt[idata*4+kk]>0.0) {
	    /* part = (0,  1) * VXY */   
	    COMPLEX_MUL2 (DFDP, Jp, VXY);
	    gradR = DFDP.real;
	    gradI = DFDP.imag;
	    /* gradR = -gradR; gradI = -gradI;   DEBUG */
	  } else gradR = gradI = 0.0;    /* invalid data */
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	}
	
	/* Loop over antenna gains */
	for (k=0; k<2; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt gX1 */
	    /* Fitting? */
	    if (antGainFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gY[ia2] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia2*2+1], 0);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RxyR wrt gX1 */
		gradI = DFDP.imag;    /* RxyI wrt gX1 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia1][0];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  case 1:  /* wrt gY2 */
	    if (antGainFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gX[ia1] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia1*2+0], 1);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RxyR wrt gY2 */
		gradI = DFDP.imag;    /* RxyI wrt gY2 */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia2][1];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  default:
	    break;
	  }; /* end gain switch */
	} /* end loop over gains */
	
	break;  /* End XY */
      	
      case 3:     /* YX */
	if (wt[idata*4+kk]>0.0) {
	  isigma = wt[idata*4+kk];
	  /* VYX = {S[0] * SY[ia1] * CXc[ia2] * DPAc +       
	            S[1] * SY[ia1] * SXc[ia2] * SPAc +
		    S[2] * CY[ia1] * CXc[ia2] * SPA  + 
		    S[3] * CY[ia1] * SXc[ia2] * DPA} * g1Y * g2X * exp(-i PD)
	  */
	  COMPLEX_MUL3 (MC1, SY[ia1], CXc[ia2], DPAc);
	  COMPLEX_MUL3 (MC2, SY[ia1], SXc[ia2], SPAc);
	  COMPLEX_MUL3 (MC3, CY[ia1], CXc[ia2], SPA);
	  COMPLEX_MUL3 (MC4, CY[ia1], SXc[ia2], DPA);
	  COMPLEX_MUL2 (SM1, S[0], MC1);
	  COMPLEX_MUL2 (SM2, S[1], MC2);
	  COMPLEX_MUL2 (SM3, S[2], MC3);
	  COMPLEX_MUL2 (SM4, S[3], MC4);
	  COMPLEX_ADD4 (VYX, SM1, SM2, SM3, SM4);
	  COMPLEX_SET (ct1,  antGain[ia1*2+1]*antGain[ia2*2+0], 0);
	  COMPLEX_EXP (ct2, -PD);
	  COMPLEX_MUL2 (ggPD, ct1, ct2);
	  COMPLEX_MUL2 (VYX, VYX, ggPD);
	  modelR = VYX.real; modelI = VYX.imag;
	  residR = modelR - data[idata*10+(kk+1)*2];
	  residI = modelI - data[idata*10+(kk+1)*2+1];
	  /* residR = -residR; residI = - residI;   DEBUG */
	  gsl_vector_set(f, i*2,   residR*isigma); /* Save function resids */
	  gsl_vector_set(f, i*2+1, residI*isigma); /* Save function resids */
	} else  residR = residI = 0.0; /* Invalid data */
	
	/* Loop over first antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex1 = 0 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Oy1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		           (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4} * 
			        gY[ia1] * gX[ia2]  * exp(-i PD) */
		COMPLEX_ADD2 (ct1, SM1, SM3);
		COMPLEX_MUL2 (ct2, Jp, ct1);
		COMPLEX_ADD2 (ct1, SM2, SM4);
		COMPLEX_MUL2 (ct3, Jm, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP,  ct2, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma); /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt Ey1 */
	    /* Fitting? */
	    if (antFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) * {S[0] * CXc[ia2] * DPAc * CYc[ia1] +       
		                     S[1] * SXc[ia2] * SPAc * CYc[ia1] +
				     S[2] * CXc[ia2] * SPA  * SYc[ia1]  + 
				     S[3] * SXc[ia2] * DPA  * SYc[ia1]} * 
				           gY[ia1] * gX[ia2] * exp(-i PD) */
		COMPLEX_MUL4(ct1, S[0], CXc[ia2], DPAc, CYc[ia1]);
		COMPLEX_MUL4(ct2, S[1], SXc[ia2], SPAc, CYc[ia1]);
		COMPLEX_MUL4(ct3, S[2], CXc[ia2], SPA,  SYc[ia1]);
		COMPLEX_MUL4(ct4, S[3], SXc[ia2], DPA,  SYc[ia1]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia1][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end first antenna parameter switch */
	} /* end loop over first antenna parameters */
	
	/* Loop over second antenna parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt Ox2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = {(0,  1) * S[0] * MC1 + (0, -1) * S[1] * MC2 + 
		           (0,  1) * S[2] * MC3 + (0, -1) * S[3] * MC4} * 
			         gY[ia1] * gX[ia2]  * exp(-i PD) */
		COMPLEX_ADD2 (ct1, SM1, SM3);
		COMPLEX_MUL2 (ct2, Jp, ct1);
		COMPLEX_ADD2 (ct1, SM2, SM4);
		COMPLEX_MUL2 (ct3, Jm, ct1);
		COMPLEX_ADD2 (ct2, ct2, ct3);
		COMPLEX_MUL2 (DFDP,  ct2, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt Ex2 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (0, -1) * {S[0] * CXc[ia2] * DPAc * CYc[ia1] +       
		                     S[1] * SXc[ia2] * SPAc * CYc[ia1] +
				     S[2] * CXc[ia2] * SPA  * SYc[ia1]  + 
				     S[3] * SXc[ia2] * DPA  * SYc[ia1]} * 
				          gY[ia1] * gX[ia2] * exp(-i PD) */
		COMPLEX_MUL4(ct1, S[0], CXc[ia2], DPAc, CYc[ia1]);
		COMPLEX_MUL4(ct2, S[1], SXc[ia2], SPAc, CYc[ia1]);
		COMPLEX_MUL4(ct3, S[2], CXc[ia2], SPA,  SYc[ia1]);
		COMPLEX_MUL4(ct4, S[3], SXc[ia2], DPA,  SYc[ia1]);
		COMPLEX_ADD4(ct5, ct1, ct2, ct3, ct4);
		COMPLEX_MUL3(DFDP, Jm, ct5, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt Oy2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt Ey2 = 0 */
	    /* Fitting? */
	    if (antFit[ia2][k]) {
	      gradR = gradI = 0.0;
	      j = antPNumb[ia2][k];
	      gsl_matrix_set(J, i*2,   j, gradR);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI);  /* Save Jacobian */
	    }
	    break;
	  default:
	    break;
	  }; /* end second antenna parameter switch */
	} /* end loop over second antenna parameters */
	
	/* Loop over source parameters */
	for (k=0; k<4; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt IPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part =   (MC1 + MC4) * gY[ia1] * gX[ia2] */
		COMPLEX_ADD2(ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 1:     /* wrt QPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (MC2 + MC3) * gY[ia1] * gX[ia2] */
		COMPLEX_ADD2(ct1, MC2, MC3);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 2:     /* wrt UPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* i (MC2 - MC3) * gY[ia1] * gX[ia2] */
		COMPLEX_SUB(ct1, MC2, MC3);
		COMPLEX_MUL3(DFDP, Jp, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	    break;
	  case 3:     /* wrt VPol */
	    /* Fitting? */
	    if (souFit[isou][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* (MC1 - MC4) * gY[ia1] * gX[ia2] */
		COMPLEX_SUB (ct1, MC1, MC4);
		COMPLEX_MUL2(DFDP, ct1, ggPD);
		gradR = DFDP.real;
		gradI = DFDP.imag;
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
	      } else gradR = gradI = 0.0;   /* invalid data */
	      j = souPNumb[isou][k];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	    }
	break;
      default:
	break;
	  }; /* end source parameter switch */
	} /* end loop over source parameters */
	
	/* gradient wrt PD */
	if (args->doFitRL) {
	  if (wt[idata*4+kk]>0.0) {
	    /* part = (0, -1) * VYX */   
	    COMPLEX_MUL2 (DFDP, Jm, VYX);
	    gradR = DFDP.real;
	    gradI = DFDP.imag;
	    /* gradR = -gradR; gradI = -gradI;   DEBUG */
	  } else gradR = gradI = 0.0;    /* invalid data */
	  j = args->PDPNumb;
	  gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	  gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	}
	/* Loop over antenna gains */
	for (k=0; k<2; k++) {
	  switch (k) {   /* Switch over parameter */
	  case 0:     /* wrt gY1 */
	    /* Fitting? */
	    if (antGainFit[ia1][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gX[ia2] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia2*2+0], 0);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RyxR wrt gY1 */
		gradI = DFDP.imag;    /* RyxI wrt gY1 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia1][0];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  case 1:  /* wrt gX2 */
	    if (antGainFit[ia2][k]) {
	      if (wt[idata*4+kk]>0.0) {
		/* part = (S[0]*MC1 + S[1]*MC2 + S[2]*MC3 + S[3]*MC4) * gY[ia1] */
		COMPLEX_ADD4(ct1, SM1, SM2, SM3, SM4);
		COMPLEX_SET (ct2, antGain[ia1*2+1], 1);
		COMPLEX_MUL2 (DFDP, ct1, ct2);
		gradR = DFDP.real;    /* RyxR wrt gX2 */
		gradI = DFDP.imag;    /* RyxI wrt gX2 */
		/* gradR = -gradR; gradI = -gradI;   DEBUG */
		/* gradR = 0.0; gradI = 0.0;  DEBUG GAIN */
	      } else gradR = gradI = 0.0;    /* invalid data */
	      j = antGainPNumb[ia2][1];
	      gsl_matrix_set(J, i*2,   j, gradR*isigma);  /* Save Jacobian */
	      gsl_matrix_set(J, i*2+1, j, gradI*isigma);  /* Save Jacobian */
	      }
	    break;
	  default:
	    break;
	  }; /* end gain switch */
	} /* end loop over gains */
	/* End YX */	
	
      default:
	  break;
	}; /* end switch over data correlation */
      
      i++;  /* Update complex datum number */
    } /* end loop over correlations */
  } /* End loop over visibilities */

  return GSL_SUCCESS;
} /*  end PolnFitFuncJacOEXY */

#endif /* HAVE_GSL */ 
