/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2013                                          */
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

#include "ObitUVSel.h"
#include "ObitUVCalSelect.h"
#include "ObitUVDesc.h"
#include "ObitUVCal.h"
#include "ObitUVCalSelect.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalSelect.c
 * ObitUVCal utilities for selecting data and translating Stokes.
 */


/*---------------Private function prototypes----------------*/
/** Private: Init data selection for Stokes I */
static void 
ObitUVCalSelectIInit (ObitUVCal *in, olong kstoke0, olong nstoke, 
		      olong ivpnt, gboolean formalI, ObitErr *err);

/** Private: Init data selection for Stokes Q */
static void 
ObitUVCalSelectQInit (ObitUVCal *in, olong kstoke0, olong nstoke, 
		      olong ivpnt, ObitErr *err);

/** Private: Init data selection for Stokes U */
static void 
ObitUVCalSelectUInit (ObitUVCal *in, olong kstoke0, olong nstoke, 
		      olong ivpnt, ObitErr *err);

/** Private: Init data selection for Stokes V */
static void 
ObitUVCalSelectVInit (ObitUVCal *in, olong kstoke0, olong nstoke, 
		      olong ivpnt, ObitErr *err);

/** Private: Init data selection for Stokes RR */
static void 
ObitUVCalSelectRRInit (ObitUVCal *in, olong kstoke0, olong nstoke, 
		      olong ivpnt, ObitErr *err);

/** Private: Init data selection for Stokes LL */
static void 
ObitUVCalSelectLLInit (ObitUVCal *in, olong kstoke0, olong nstoke, 
		      olong ivpnt, ObitErr *err);

/** Private: Init data selection for Stokes RL */
static void 
ObitUVCalSelectRLInit (ObitUVCal *in, olong kstoke0, olong nstoke, 
		      olong ivpnt, ObitErr *err);

/** Private: Init data selection for Stokes LR */
static void 
ObitUVCalSelectLRInit (ObitUVCal *in, olong kstoke0, olong nstoke, 
		      olong ivpnt, ObitErr *err);

/** Private: Init data selection for Stokes XX */
static void 
ObitUVCalSelectXXInit (ObitUVCal *in, olong kstoke0, olong nstoke, 
		      olong ivpnt, ObitErr *err);

/** Private: Init data selection for Stokes YY */
static void 
ObitUVCalSelectYYInit (ObitUVCal *in, olong kstoke0, olong nstoke, 
		      olong ivpnt, ObitErr *err);

/** Private: Init data selection for Stokes XY */
static void 
ObitUVCalSelectXYInit (ObitUVCal *in, olong kstoke0, olong nstoke, 
		      olong ivpnt, ObitErr *err);

/** Private: Init data selection for Stokes YX */
static void 
ObitUVCalSelectYXInit (ObitUVCal *in, olong kstoke0, olong nstoke, 
		      olong ivpnt, ObitErr *err);

/*----------------------Public functions---------------------------*/

/**
 * Initialize structures for data selection.
 * Output descriptor modified to reflect data selection.
 * Stokes 'F' = formal I (both RR+LL or XX+YY needed)
 * On output data are ordered Stokes, Freq, IF if Stokes translation
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in   Object to initialize.
 * \li mySel->transPol  if TRUE translate stokes
 * \li mySel->bothCorr  if TRUE all data needed in combination
 * \li jadr  (0-rel) indices  to the first and second visibility
 *     input records to be used in the output record. 
 *     If jadr(1,n) is negative (Stokes U from RL,LR) use
 *     abs(jadr(1,n)) and multiply the visibility by i (=sqrt(-1)).  
 *     If jadr(2,n) < 0, (RL,LR from Q,U) use abs and  construct (1) + i*(2)
 * \li selFact   Factors to be multiplied by the first and
 *      second input vis's to make the output vis.
 * \param sel     Data selector.
 * \param inDesc  Input  data descriptor.
 * \param outDesc Output data descriptor (after transformations/selection).
 * \param err     ObitError stack.
 */
void ObitUVCalSelectInit (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *inDesc, 
			  ObitUVDesc *outDesc, ObitErr *err)
{
  olong   i, ifb, ife, icb, ice, nmode=21;
  olong   ierr, pmode, nchan, nif, nstok, ivpnt, kstoke0;
  olong   knpoln, knchan, knif;
  ofloat crpixpoln=0.0, crpixchan=0.0, crpixif=0.0, crotapoln=0.0, crotachan=0.0, crotaif=0.0;
  ofloat cdeltpoln=0.0, cdeltchan=0.0, cdeltif=0.0;
  gchar  ctypepoln[9], ctypechan[9], ctypeif[9];
  odouble crvalpoln=0.0, crvalchan=0.0, crvalif=0.0;
  gchar chmode[21][5] = 
    {"I   ","V   ","Q   ", "U   ", "IQU ", "IQUV", 
     "IV  ", "RR  ", "LL  ", "RL  ", "LR  ", "HALF", "FULL",
     "F   ", "FQU ", "FQUV", "IV  ",
     "XX  ", "YY  ", "XY  ", "YX  "};
  gboolean formalI;
  odouble crval=0.0;
  ofloat cdelt=0.0;

  /* data selection parameters on in 
  in->jadr;
  in->selFact;*/

  /* Update output descriptor */
  /* Channel selection */
  outDesc->inaxes[outDesc->jlocf] = in->mySel->numberChann;
  outDesc->crpix[outDesc->jlocf] -= in->mySel->startChann - 1;
  outDesc->altCrpix -= in->mySel->startChann - 1;

  /* IF selection */
  if (outDesc->jlocif>=0) {
    outDesc->inaxes[outDesc->jlocif] = in->mySel->numberIF;
    outDesc->crpix[outDesc->jlocif] -= in->mySel->startIF - 1;
    /* Update frequency */
    outDesc->crval[outDesc->jlocf] = inDesc->freqIF[in->mySel->startIF-1];
    outDesc->cdelt[outDesc->jlocf] = inDesc->chIncIF[in->mySel->startIF-1];
  }
  
  /* Reorder axes to Poln, Freq, IF - save */
  if (outDesc->jlocs>=0) {
    knpoln   = outDesc->inaxes[outDesc->jlocs];
    for (i=0; i<8; i++) 
      ctypepoln[i] = outDesc->ctype[outDesc->jlocs][i];
    crpixpoln= outDesc->crpix[outDesc->jlocs];
    crotapoln= outDesc->crota[outDesc->jlocs];
    cdeltpoln= outDesc->cdelt[outDesc->jlocs];
    crvalpoln= outDesc->crval[outDesc->jlocs];
  } else knpoln = 0;

  if (outDesc->jlocf>=0) {
    knchan   = outDesc->inaxes[outDesc->jlocf];
    for (i=0; i<8; i++) 
      ctypechan[i] = outDesc->ctype[outDesc->jlocf][i];
    crpixchan= outDesc->crpix[outDesc->jlocf];
    crotachan= outDesc->crota[outDesc->jlocf];
    cdeltchan= outDesc->cdelt[outDesc->jlocf];
    crvalchan= outDesc->crval[outDesc->jlocf];
  } else knchan = 0;

  if (outDesc->jlocif>=0) {
    knif   = outDesc->inaxes[outDesc->jlocif];
    for (i=0; i<8; i++) 
      ctypeif[i] = outDesc->ctype[outDesc->jlocif][i];
    crpixif= outDesc->crpix[outDesc->jlocif];
    crotaif= outDesc->crota[outDesc->jlocif];
    cdeltif= outDesc->cdelt[outDesc->jlocif];
    crvalif= outDesc->crval[outDesc->jlocif];
  } else knif = 0;

  /* Put back in order */
  if (knpoln>0) {
    outDesc->jlocs = 1;
    for (i=0; i<8; i++) 
      outDesc->ctype[outDesc->jlocs][i] = ctypepoln[i];
    outDesc->inaxes[outDesc->jlocs] = knpoln;
    outDesc->crpix[outDesc->jlocs] = crpixpoln;
    outDesc->crota[outDesc->jlocs] = crotapoln;
    outDesc->cdelt[outDesc->jlocs] = cdeltpoln;
    outDesc->crval[outDesc->jlocs] = crvalpoln;
  } 

  if (knchan>0) {
    outDesc->jlocf = outDesc->jlocs+1;
    for (i=0; i<8; i++) 
      outDesc->ctype[outDesc->jlocf][i] = ctypechan[i];
    outDesc->inaxes[outDesc->jlocf] = knchan;
    outDesc->crpix[outDesc->jlocf] = crpixchan;
    outDesc->crota[outDesc->jlocf] = crotachan;
    outDesc->cdelt[outDesc->jlocf] = cdeltchan;
    outDesc->crval[outDesc->jlocf] = crvalchan;
  } 

  if (knif>0) {
    outDesc->jlocif = outDesc->jlocf+1;
    for (i=0; i<8; i++) 
      outDesc->ctype[outDesc->jlocif][i] = ctypeif[i];
    outDesc->inaxes[outDesc->jlocif] = knif;
    outDesc->crpix[outDesc->jlocif] = crpixif;
    outDesc->crota[outDesc->jlocif] = crotaif;
    outDesc->cdelt[outDesc->jlocif] = cdeltif;
    outDesc->crval[outDesc->jlocif] = crvalif;
  } 

  /* Reindex output descriptor */
  ObitUVDescIndex (outDesc);
  sel->lrecUC   = outDesc->lrec;
  sel->nrparmUC = outDesc->nrparm;

  /* do we need to translate data? No if Stokes ="    " */
  in->mySel->transPol = strncmp(sel->Stokes, "    ", 4);
  if (!in->mySel->transPol) return;

 /* Compute translation parameters, find poln. mode */
  pmode = -1;
  for (i= 0; i<nmode; i++) 
    if (!strncmp(sel->Stokes, chmode[i], 4)) pmode = i;
  formalI = pmode>12;
	
  /* Unrecognized Stokes' */
  if (pmode < 0) {
    Obit_log_error(err, OBIT_Error, 
		   "Unknown Stokes request %s for %s", sel->Stokes, in->name);
    return;
  } 
  
  /* Get first stokes index */
  if (inDesc->crval[inDesc->jlocs]> 0.0) {
    kstoke0 = inDesc->crval[inDesc->jlocs] + 
      (1.0-inDesc->crpix[inDesc->jlocs]) * inDesc->cdelt[inDesc->jlocs] + 0.5;
  } else {
    kstoke0 = inDesc->crval[inDesc->jlocs] + 
      (1.0-inDesc->crpix[inDesc->jlocs]) * inDesc->cdelt[inDesc->jlocs] - 0.5;
  } 

  /* linear polarized data (x-y), assume that it is being 
     calibrated and will be changed  to rr,ll,rl,lr data.
  if (kstoke0 <= -5) kstoke0 = kstoke0 + 4;  Maybe not */

  /* get start, end channels, IFs from Selector */
  ifb = sel->startIF;
  ife = sel->startIF + sel->numberIF - 1;
  icb = sel->startChann;
  ice = sel->startChann + sel->numberChann - 1;

  /* number channels, IFs */
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocf>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;

  /* Number of Stokes parameters in input */
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;

  ierr = 1;
  /* check channel; default = all */
  if (icb <= 0)  ice = nchan;
  icb = MAX (icb, 1);
  ice = MAX (ice, 1);
  /* Force to range */
  icb = MIN (icb, nchan);
  ice = MIN (ice, nchan);
  icb = MIN (icb, ice);
 
  /* IF.  default is all */
  if (ifb <= 0) ife = nif;
  ifb = MAX (ifb, 1);
  ife = MAX (ife, 1);
  /* Force to range */
  ifb = MIN (ifb, nif);
  ife = MIN (ife, nif);
  ifb = MIN (ifb, ife);

  /* Are several correlations needed?  */
  in->mySel->bothCorr = FALSE;  /* Unless determined otherwise */

  /* initialize visibility index */
  ivpnt = 0;

  /* branch by pmode */
  switch (pmode) {
  case 0:  /* I    */
  case 13: /* F    */
    ObitUVCalSelectIInit(in, kstoke0, nstok, ivpnt, formalI, err);
    crval = 1.0; /* Stokes index */
    cdelt = 1.0; /* stokes increment */
    ivpnt++;
    break;

  case 1: /* V    */
    ObitUVCalSelectVInit (in, kstoke0, nstok, ivpnt, err);
    crval = 4.0; /* Stokes index */
    cdelt = 1.0; /* stokes increment */
    ivpnt++;
    break;
    
  case 2: /* Q    */
    ObitUVCalSelectQInit (in, kstoke0, nstok, ivpnt, err);
    crval = 2.0; /* Stokes index */
    cdelt = 1.0; /* stokes increment */
    ivpnt++;
    break;
    
  case 3: /* U    */
    ObitUVCalSelectUInit (in, kstoke0, nstok, ivpnt, err);
    crval = 3.0; /* Stokes index */
    cdelt = 1.0; /* stokes increment */
    ivpnt++;
    break;
    
  case 4:  /* IQU  */
  case 14: /* FQU   */
    crval = 1.0; /* Stokes index */
    cdelt = 1.0; /* stokes increment */
    ObitUVCalSelectIInit (in, kstoke0, nstok, ivpnt, formalI, err);
    ivpnt++;
    ObitUVCalSelectQInit (in, kstoke0, nstok, ivpnt, err);
    ivpnt++;
    ObitUVCalSelectUInit (in, kstoke0, nstok, ivpnt, err);
    ivpnt++;
    break;
    
  case 5:  /* IQUV */
  case 15: /* FQUV   */
   crval = 1.0; /* Stokes index */
    cdelt = 1.0; /* stokes increment */
    ObitUVCalSelectIInit (in, kstoke0, nstok, ivpnt, formalI, err);
    ivpnt++;
    ObitUVCalSelectQInit (in, kstoke0, nstok, ivpnt, err);
    ivpnt++;
    ObitUVCalSelectUInit (in, kstoke0, nstok, ivpnt, err);
    ivpnt++;
    ObitUVCalSelectVInit (in, kstoke0, nstok, ivpnt, err);
    ivpnt++;
    break;
    
  case 6:  /* IV   */
  case 16: /* FV   */
    ObitUVCalSelectIInit (in, kstoke0, nstok, ivpnt, formalI, err);
    crval = 1.0; /* Stokes index */
    cdelt = 4.0; /* stokes increment */
    ivpnt++;
    ObitUVCalSelectVInit (in, kstoke0, nstok, ivpnt, err);
    ivpnt++;
    break;
    
  case 7: /* RR   */
    ObitUVCalSelectRRInit (in, kstoke0, nstok, ivpnt, err);
    crval = -1.0; /* Stokes index */
    cdelt =  1.0; /* stokes increment */
    ivpnt++;
    break;

  case 8: /* LL   */
    ObitUVCalSelectLLInit (in, kstoke0, nstok, ivpnt, err);
    crval = -2.0; /* Stokes index */
    cdelt =  1.0; /* stokes increment */
    ivpnt++;
    break;

  case 9: /* RL   */
    ObitUVCalSelectRLInit (in, kstoke0, nstok, ivpnt, err);
    crval = -3.0; /* Stokes index */
    cdelt =  1.0; /* stokes increment */
    ivpnt++;
    break;

  case 10: /* LR  */
    ObitUVCalSelectLRInit (in, kstoke0, nstok, ivpnt, err);
    crval = -4.0; /* Stokes index */
    cdelt =  1.0; /* stokes increment */
    ivpnt++;
    break;

  case 11: /* HALF = RR,LL  */
    ObitUVCalSelectRRInit (in, kstoke0, nstok, ivpnt, err);
    crval = -1.0; /* Stokes index */
    cdelt =  1.0; /* stokes increment */
    ivpnt++;
    ObitUVCalSelectLLInit (in, kstoke0, nstok, ivpnt, err);
    ivpnt++;
    break;

  case 12: /* FULL */
    ObitUVCalSelectRRInit (in, kstoke0, nstok, ivpnt, err);
    crval = -1.0; /* Stokes index */
    cdelt =  1.0; /* stokes increment */
    ivpnt++;
    ObitUVCalSelectLLInit (in, kstoke0, nstok, ivpnt, err);
    ivpnt++;
    ObitUVCalSelectRLInit (in, kstoke0, nstok, ivpnt, err);
    ivpnt++;
    ObitUVCalSelectLRInit (in, kstoke0, nstok, ivpnt, err);
    ivpnt++;
    break;

  case 17: /* XX */
    ObitUVCalSelectXXInit (in, kstoke0, nstok, ivpnt, err);
    crval = -5.0; /* Stokes index */
    cdelt =  1.0; /* stokes increment */
    ivpnt++;
    break;

  case 18: /* YY */
    ObitUVCalSelectYYInit (in, kstoke0, nstok, ivpnt, err);
    crval = -6.0; /* Stokes index */
    cdelt =  1.0; /* stokes increment */
    ivpnt++;
    break;

  case 19: /* XY */
    ObitUVCalSelectXYInit (in, kstoke0, nstok, ivpnt, err);
    crval = -7.0; /* Stokes index */
    cdelt =  1.0; /* stokes increment */
    ivpnt++;
    break;

  case 20: /* YX */
    ObitUVCalSelectYXInit (in, kstoke0, nstok, ivpnt, err);
    crval = -8.0; /* Stokes index */
    cdelt =  1.0; /* stokes increment */
    ivpnt++;
    break;

  default: /* should never get here */
    g_assert_not_reached();
    break;
  }; /* end switch on poln mode */

  /* How many polarizations out */
  in->mySel->numberPoln  = ivpnt;

  /* How many total visibilities? */
  in->mySel->numberVis  = 
    ivpnt * sel->numberIF * sel->numberChann; 

  /* Update output descriptor */
  /* Stokes selection */
  outDesc->inaxes[outDesc->jlocs] = ivpnt;
  outDesc->crpix[outDesc->jlocs] = 1;
  outDesc->crval[outDesc->jlocs] = crval;
  outDesc->cdelt[outDesc->jlocs] = cdelt;

  /* Reindex output descriptor */
  ObitUVDescIndex (outDesc);

  /* Change Selector to reflect true output data size */
  sel->lrecUC   = outDesc->lrec;
  sel->nrparmUC = outDesc->nrparm;

} /* end ObitUVCalSelectInit */

/**
 * Select data and translate to desired Stokes parameter
 * Adapted from the AIPSish DGGET.FOR.
 * U, V, W scales if selecting in IF. (done in ObitUVCal)
 * \param in     Calibration Object.
 * \li mySel->transPol  if TRUE translate stokes
 * \li mySel->bothCorr  if TRUE all data needed in combination
 * \li jadr  (0-rel) indices  to the first and second visibility
 *     input records to be used in the output record. 
 *     If jadr(1,n) is negative (Stokes U from RL,LR) use
 *     abs(jadr(1,n)) and multiply the visibility by i (=sqrt(-1)).  
 *     If jadr(2,n) < 0, (RL,LR from Q,U) use abs and  construct (1) + i*(2)
 * \li selFact   Factors to be multiplied by the first and
 *     second input vis's to make the output vis.
 * \param RP     Random parameters array.
 * \param visIn  input visibility as an array of floats
 * \param visOut output visibility as an array of floats
 * \param err    ObitError stack.
 * \returns TRUE if at least some data is valid
 */
gboolean ObitUVCalSelect (ObitUVCal *in, ofloat *RP, ofloat *visIn, ofloat *visOut, 
			  ObitErr *err)
{
  olong   i, ip1, ip2, op, lf, lc, lp, lfoff, ioff, incf, incif, incs;
  olong channInc, IFInc, maxChan, maxIF;
  gboolean   good;
  ofloat wt1, wt2, temp, xfact1, xfact2;
  ObitUVSel *sel = in->mySel;
  ObitUVDesc *desc = in->myDesc;

  /* existing error? */
  if (err->error) return FALSE;
  
  good = FALSE;
  /* see if weight in data */

  /* check input. */
  if (sel->numberVis <= 0) return good;

  /* Data will always be uncompressed when it gets here, 
     make sure data increments are correct */
  incs  = 3 * desc->incs  / desc->inaxes[0];
  incf  = 3 * desc->incf  / desc->inaxes[0];
  incif = 3 * desc->incif / desc->inaxes[0];
  channInc = MAX (1, sel->channInc);
  maxChan  = MIN (sel->numberChann*channInc, desc->inaxes[desc->jlocf]);
  IFInc    = MAX (1, sel->IFInc);
  maxIF    =  desc->inaxes[desc->jlocif];

  /* loop checking  and summing weights and getting visibilities. */
  i = -1; /* output visibility number */

  /* initial visibility offset in input */
  lfoff = (sel->startChann-1-channInc) * incf - incif;

  /* loop over IF */
  for (lf=0; lf<maxIF; lf++) { /* loop 70 */
    lfoff += incif;
    ioff   = lfoff;
    /* This one wanted? */
    if (!sel->IFSel[lf]) continue;

    /* Loop over channel */
    for (lc=0; lc<maxChan; lc += channInc) { /* loop 60 */
      ioff += incf * channInc;

      /* translating stokes? */
      if (sel->transPol) {
	for (lp= 0; lp<sel->numberPoln; lp++) { /* loop 40 */
	  i++; /* output correlation number */

	  /* initialize output */
	  visOut[3*i]   = 0.0;
	  visOut[3*i+1] = 0.0;
	  visOut[3*i+2] = 0.0;

	  /* set input visibility indices. */
	  ip1 = abs (in->jadr[lp][0]) + ioff;
	  ip2 = abs (in->jadr[lp][1]) + ioff;

	  /* get weights. */
	  wt1 = MAX (0.0, visIn[ip1+2]) * 0.5;
	  wt2 = MAX (0.0, visIn[ip2+2]) * 0.5;

	  /* Do we need all data and only one present? */
	  if (!((sel->bothCorr)  &&  ((wt1 <= 0.0)  ||  (wt2 <= 0.0)))) {

	    /* set visibility. */
	    visOut[3*i+2] = wt1 + wt2;
	    if (visOut[3*i+2] > 0.0) {
	      xfact1 = 0.0;
	      xfact2 = 0.0;
	      if (wt1 > 0.0) xfact1 = in->selFact[lp][0];
	      if (wt2 > 0.0) xfact2 = in->selFact[lp][1];
	      good = TRUE;
	      
	      /* RL, LR from q, u */
	      if (in->jadr[lp][1] < 0) {
		visOut[3*i]   = xfact1 * visIn[ip1] - xfact2 * visIn[ip2+1];
		visOut[3*i+1] = xfact1 * visIn[ip1+1] + xfact2 * visIn[ip2];
	      } else {
		visOut[3*i]   = xfact1 * visIn[ip1] + xfact2 * visIn[ip2];
		visOut[3*i+1] = xfact1 * visIn[ip1+1] + xfact2 * visIn[ip2+1];
		/* upol, mult. by i */
		if (in->jadr[lp][0] < 0) {
		  temp = visOut[3*i];
		  visOut[3*i]   = - visOut[3*i+1];
		  visOut[3*i+1] = temp;
		} 
	      } 
	      /* if one wt <= 0  double vis. */
	      if ((wt1 <= 0.0)  ||  (wt2 <= 0.0)) {
		visOut[3*i]   = 2.0 * visOut[3*i];
		visOut[3*i+1] = 2.0 * visOut[3*i+1];
	      } 
	    } 
	  } /* end loop  L40:  */;
	} /* end need all weights */

	/* no translation - just copy and reorder data */
      } else {
	for (lp=0; lp<in->numStok; lp++) { /* loop 50 */
	  i++;
	  /* ip1 = abs (in->jadr[lp][0]) + ioff;*/
	  /* set input visibility index . */
	  ip1 = ioff +  lp * incs;
	  op  = i*3; 

	  /* get weight. */
	  wt1 = MAX (0.0, visIn[ip1+2]);

	  /* set visibility. */
	  /* Does Stokes come first in output? */
	  if (incs==3) {
	    visOut[3*i+2] = wt1;
	    if (wt1 > 0.0) good = TRUE;
	    visOut[3*i]   = visIn[ip1];
	    visOut[3*i+1] = visIn[ip1+1];
	  } else { /* no reording */
	    /* This may not always be correct */
	    visOut[op+2] = wt1;
	    if (wt1 > 0.0) good = TRUE;
	    visOut[op]   = visIn[ip1];
	    visOut[op+1] = visIn[ip1+1];
	  }
	} /* end loop  L50:  */;
      } 
    } /* end loop  L60:  */;
  } /* end loop  L70:  */;

  return good;
} /* end ObitUVCalSelect */

/*---------------Private functions--------------------------*/

/**
 * Initialize structures for data selection for Ipol.
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in      Object to initialize.
 * \param kstoke0 Stokes index of first correlation
 *                (1-4) => I,Q,U,V; 
 *                (-1 -> -4)=> RR, LL, RL, LR
 *                (-5 -> -8)=> XX, YY, XY, YX
 * \param nstok  Number of Stokes correlations
 * \param ivpnt   Index in data selection structures
 *                Updated on output.
 * \param formalI If true formal I (both RR+LL or XX+YY) needed
 * \param err     ObitError stack.
 */
static void 
ObitUVCalSelectIInit (ObitUVCal *in, olong kstoke0, olong nstok, 
		      olong ivpnt, gboolean formalI, ObitErr *err)
{
  gboolean missing;
  olong incs;

  /* error check */
  if (err->error) return;

  /* Make sure data there */
  missing = ((kstoke0 > 1) || (kstoke0 < -2));
  /* Linear poloarization */
  if (kstoke0<=-5) missing = (kstoke0 < -6);
  if (formalI && (kstoke0<0))  
    missing = missing || ((nstok<2) || ((kstoke0!=-1) && (kstoke0!=-5)));
  if (missing) {
    if (formalI)
      Obit_log_error(err, OBIT_Error, 
		     "Stokes Formal I not available for %s", in->name);
    else
      Obit_log_error(err, OBIT_Error, 
		     "Stokes I not available for %s", in->name);
    return;
  }
    
  /* Data will always be uncompressed by the time it gets here */
  incs = 3 * in->myDesc->incs / in->myDesc->inaxes[0];
  if (kstoke0 > 0) {
    /* true Stokes parameters */
    in->jadr[ivpnt][0] = (1-kstoke0) * incs;
    in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
    in->selFact[ivpnt][0] = 1.0;
    in->selFact[ivpnt][1] = 0.0;
  } else {
    /* have RR, LL correlations */
    in->mySel->bothCorr = formalI;  /* Need both? */
    in->jadr[ivpnt][0] = 0;
    in->jadr[ivpnt][1] = in->jadr[ivpnt][0] + incs;
    in->selFact[ivpnt][0] = 0.5;
    in->selFact[ivpnt][1] = 0.5;
    /* check if only RR or LL and if so use it. */
    if ((nstok < 2)  ||  (kstoke0 != -1)) 
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
  }
} /* end ObitUVCalSelectIInit */

/**
 * Initialize structures for data selection for Qpol.
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in      Object to initialize.
 * \param kstoke0 Stokes index of first correlation
 *                (1-4) => I,Q,U,V; 
 *                (-1 -> -4)=> RR, LL, RL, LR
 * \param nstok   Number of Stokes correlations
 * \param ivpnt   Index in data selection structures
 *                Updated on output.
 * \param err     ObitError stack.
 */
static void 
ObitUVCalSelectQInit (ObitUVCal *in, olong kstoke0, olong nstok, 
		      olong ivpnt, ObitErr *err)
{
  gboolean missing;
  olong incs;

  /* error check */
  if (err->error) return;
    
    /* Add Stokes Q -  check if in data */
    missing = ((kstoke0>0) && ((kstoke0>2) || (kstoke0+nstok-1<2)));
    missing = missing ||
      ((kstoke0<0) && ((nstok<2) || (kstoke0<-3) || (kstoke0-nstok+1>-4)));
    if (missing) {
      Obit_log_error(err, OBIT_Error, 
		     "Stokes Q not available for %s", in->name);
      return;
    }

    /* Data will always be uncompressed by the time it gets here */
    incs = 3 * in->myDesc->incs / in->myDesc->inaxes[0];
    if (kstoke0 >= 0) {
      /* true Stokes. */
      in->jadr[ivpnt][0] = (2-kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 0.0;
      
    } else {
      /* RL, LR  */
      in->mySel->bothCorr = TRUE;  /* Need both */
      in->jadr[ivpnt][0] = (3+kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0] + incs;
      in->selFact[ivpnt][0] = 0.5;
      in->selFact[ivpnt][1] = 0.5;
    } 

} /* end ObitUVCalSelectQInit */

/**
 * Initialize structures for data selection for Upol.
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in      Object to initialize.
 * \param kstoke0 Stokes index of first correlation
 *                (1-4) => I,Q,U,V; 
 *                (-1 -> -4)=> RR, LL, RL, LR
 * \param nstok   Number of Stokes correlations
 * \param ivpnt   Index in data selection structures
 *                Updated on output.
 * \param err     ObitError stack.
 */
static void 
ObitUVCalSelectUInit (ObitUVCal *in, olong kstoke0, olong nstok, 
		      olong ivpnt, ObitErr *err)
{
  gboolean missing;
  olong incs;

  /* error check */
  if (err->error) return;

    /* Add Stokes U -  check if in data */
    missing = ((kstoke0>0) && ((kstoke0>3) || (kstoke0+nstok-1<3)));
    missing = missing ||
      ((kstoke0<0) && ((nstok<2) || (kstoke0<-3) || (kstoke0-nstok+1>-4)));
    if (missing) {
      Obit_log_error(err, OBIT_Error, 
		     "Stokes U not available for %s", in->name);
      return;
    }

    /* Data will always be uncompressed by the time it gets here */
    incs = 3 * in->myDesc->incs / in->myDesc->inaxes[0];
    
    if (kstoke0 >= 0) {
      /* true stokes. check. */
      in->jadr[ivpnt][0] = (3-kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 0.0;
      
    } else {
      /* RL, LR. */
      in->mySel->bothCorr = TRUE;  /* Need both */
      in->jadr[ivpnt][0] = (3+kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0] + incs;
      /* make jadr(1,ivpnt) negative to indicate to multiply 
	 by i (= sqrt (-1)) */
      in->jadr[ivpnt][0] = -in->jadr[ivpnt][0];
      in->selFact[ivpnt][0] = -0.5;
      in->selFact[ivpnt][1] = 0.5;
    } 
} /* end ObitUVCalSelectUInit */

/**
 * Initialize structures for data selection for Vpol.
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in      Object to initialize.
 * \param kstoke0 Stokes index of first correlation
 *                (1-4) => I,Q,U,V; 
 *                (-1 -> -4)=> RR, LL, RL, LR
 * \param nstok   Number of Stokes correlations
 * \param ivpnt   Index in data selection structures
 *                Updated on output.
 * \param err     ObitError stack.
 */
static void 
ObitUVCalSelectVInit (ObitUVCal *in, olong kstoke0, olong nstok, 
		      olong ivpnt, ObitErr *err)
{
  gboolean missing;
  olong incs;

  /* error check */
  if (err->error) return;

    
    /* Add Stokes V -  check if in data */
    missing = ((kstoke0>0) && ((kstoke0>4) || (kstoke0+nstok-1<4)));
    missing = missing ||
      ((kstoke0<0) && ((nstok<2) || (kstoke0<-1) || (kstoke0-nstok+1>-2)));
    if (missing) {
      Obit_log_error(err, OBIT_Error, 
		     "Stokes V not available for %s", in->name);
      return;
    }
    
    /* Data will always be uncompressed by the time it gets here */
    incs = 3 * in->myDesc->incs / in->myDesc->inaxes[0];
    
    if (kstoke0 >= 0) {
      /* true Vpol. */
      in->jadr[ivpnt][0] = (4-kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 0.0;
      
    } else {
      /* RR, LL */
      in->mySel->bothCorr = TRUE;  /* Need both */
      in->jadr[ivpnt][0] = (1+kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0] + incs;
      in->selFact[ivpnt][0] = 0.5;
      in->selFact[ivpnt][1] = -0.5;
    } 
} /* end ObitUVCalSelectVInit */

/**
 * Initialize structures for data selection for RR.
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in      Object to initialize.
 * \param kstoke0 Stokes index of first correlation
 *                (1-4) => I,Q,U,V; 
 *                (-1 -> -4)=> RR, LL, RL, LR
 * \param nstok   Number of Stokes correlations
 * \param ivpnt   Index in data selection structures
 *                Updated on output.
 * \param err     ObitError stack.
 */
static void 
ObitUVCalSelectRRInit (ObitUVCal *in, olong kstoke0, olong nstok, 
		      olong ivpnt, ObitErr *err)
{
  gboolean missing;
  olong incs;

  /* error check */
  if (err->error) return;

     /* Add Stokes RR -  check if in data */
    missing = ((kstoke0>0) && ((nstok<4) || (kstoke0>1) || (kstoke0+nstok-1<4)));
    missing = missing ||
      ((kstoke0<0) && ((kstoke0<-1) || (kstoke0-nstok+1>-1)));
    if (missing) {
      Obit_log_error(err, OBIT_Error, 
		     "Stokes RR not available for %s", in->name);
      return;
    }

    /* Data will always be uncompressed by the time it gets here */
    incs = 3 * in->myDesc->incs / in->myDesc->inaxes[0];
    if (kstoke0 > 0) {
      /* true Stokes (need IQUV for this to work) */
      in->mySel->bothCorr = TRUE;  /* Need both */
      in->jadr[ivpnt][0] = (1-kstoke0) * incs;
      in->jadr[ivpnt][1] = (4-kstoke0) * incs;
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 1.0;

    } else {
      /* RR */
      in->jadr[ivpnt][0] = 0;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 0.0;
    } 
} /* end ObitUVCalSelectRRInit */

/**
 * Initialize structures for data selection for LL.
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in      Object to initialize.
 * \param kstoke0 Stokes index of first correlation
 *                (1-4) => I,Q,U,V; 
 *                (-1 -> -4)=> RR, LL, RL, LR
 * \param nstok   Number of Stokes correlations
 * \param ivpnt   Index in data selection structures
 *                Updated on output.
 * \param err     ObitError stack.
 */
static void 
ObitUVCalSelectLLInit (ObitUVCal *in, olong kstoke0, olong nstok, 
		      olong ivpnt, ObitErr *err)
{
  gboolean missing;
  olong incs;

  /* error check */
  if (err->error) return;

    /* Add Stokes LL -  check if in data */
    missing = ((kstoke0>0) && ((nstok<4) || (kstoke0>1) || (kstoke0+nstok-1<4)));
    missing = missing ||
      ((kstoke0<0) && ((kstoke0<-2) || (kstoke0-nstok+1>-2)));
    if (missing) {
      Obit_log_error(err, OBIT_Error, 
		     "Stokes LL not available for %s", in->name);
      return;
    }
    /* true Stokes.(need IQUV for this to work) */

    /* Data will always be uncompressed by the time it gets here */
    incs = 3 * in->myDesc->incs / in->myDesc->inaxes[0];
    
    if (kstoke0 > 0) {
      in->mySel->bothCorr = TRUE;  /* Need both */
      in->jadr[ivpnt][0] = (1-kstoke0) * incs;
      in->jadr[ivpnt][1] = (4-kstoke0) * incs;
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = -1.0;
    } else {
      /* LL. */
      in->jadr[ivpnt][0] = (2+kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 0.0;
    } 
} /* end ObitUVCalSelectLLInit */

/**
 * Initialize structures for data selection for RL.
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in      Object to initialize.
 * \param kstoke0 Stokes index of first correlation
 *                (1-4) => I,Q,U,V; 
 *                (-1 -> -4)=> RR, LL, RL, LR
 * \param nstok   Number of Stokes correlations
 * \param ivpnt   Index in data selection structures
 *                Updated on output.
 * \param err     ObitError stack.
 */
static void 
ObitUVCalSelectRLInit (ObitUVCal *in, olong kstoke0, olong nstok, 
		      olong ivpnt, ObitErr *err)
{
  gboolean missing;
  olong incs;

  /* error check */
  if (err->error) return;

    /* Add Stokes LL -  check if in data */
    missing = ((kstoke0>0) && ((nstok<4) || (kstoke0>1) || (kstoke0+nstok-1<4)));
    missing = missing ||
      ((kstoke0<0) && ((kstoke0<-3) || (kstoke0-nstok+1>-2)));
    if (missing) {
      Obit_log_error(err, OBIT_Error, 
		     "Stokes LL not available for %s", in->name);
      return;
    }

    /* Data will always be uncompressed by the time it gets here */
    incs = 3 * in->myDesc->incs / in->myDesc->inaxes[0];
    
    if (kstoke0 > 0) {
      /* true Stokes (need QU for this to work)*/
      in->mySel->bothCorr = TRUE;  /* Need both */
      in->jadr[ivpnt][0] = (2-kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0] + incs;
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 1.0;
      /* q + iu code */
      in->jadr[ivpnt][1] = -in->jadr[ivpnt][1];
    } else {
      /* RL */
      in->jadr[ivpnt][0] = (3+kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 0.0;
    } 
} /* end ObitUVCalSelectRLInit */

/**
 * Initialize structures for data selection for LR.
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in      Object to initialize.
 * \param kstoke0 Stokes index of first correlation
 *                (1-4) => I,Q,U,V; 
 *                (-1 -> -4)=> RR, LL, RL, LR
 * \param nstok   Number of Stokes correlations
 * \param ivpnt   Index in data selection structures
 *                Updated on output.
 * \param err     ObitError stack.
 */
static void 
ObitUVCalSelectLRInit (ObitUVCal *in, olong kstoke0, olong nstok, 
		      olong ivpnt, ObitErr *err)
{
  gboolean missing;
  olong incs;

  /* error check */
  if (err->error) return;

    /* Add Stokes LR -  check if in data */
    missing = ((kstoke0>0) && ((nstok<4) || (kstoke0>1) || (kstoke0+nstok-1<4)));
    missing = missing ||
      ((kstoke0<0) && ((kstoke0<-4) || (kstoke0-nstok+1>-4)));
    if (missing) {
      Obit_log_error(err, OBIT_Error, 
		     "Stokes LR not available for %s", in->name);
      return;
    }

    /* Data will always be uncompressed by the time it gets here */
    incs = 3 * in->myDesc->incs / in->myDesc->inaxes[0];
    
    if (kstoke0 > 0) {
      /* true Stokes (need QU for this to work) */
      in->mySel->bothCorr = TRUE;  /* Need both */
      in->jadr[ivpnt][0] = (2-kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0] + incs;
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = -1.0;
      /* q + iu code */
      in->jadr[ivpnt][1] = -in->jadr[ivpnt][1];
    } else {
      /* LR. */
      in->jadr[ivpnt][0] = (4+kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 0.0;
    } 
} /* end ObitUVCalSelectLRInit */

/**
 * Initialize structures for data selection for RR.
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in      Object to initialize.
 * \param kstoke0 Stokes index of first correlation
 *                (1-4) => I,Q,U,V; 
 *                (-5 -> -8)=> XX, YY, XY, YX
 * \param nstok   Number of Stokes correlations
 * \param ivpnt   Index in data selection structures
 *                Updated on output.
 * \param err     ObitError stack.
 */
static void 
ObitUVCalSelectXXInit (ObitUVCal *in, olong kstoke0, olong nstok, 
		      olong ivpnt, ObitErr *err)
{
  gboolean missing;
  olong incs;

  /* error check */
  if (err->error) return;

     /* Add Stokes XX -  check if in data */
    missing = ((kstoke0>0) && ((nstok<4) || (kstoke0>1) || (kstoke0+nstok-1<4)));
    missing = missing ||
      ((kstoke0<0) && ((kstoke0<-5) || (kstoke0-nstok+1>-5)));
    if (missing) {
      Obit_log_error(err, OBIT_Error, 
		     "Stokes XX not available for %s", in->name);
      return;
    }

    /* Data will always be uncompressed by the time it gets here */
    incs = 3 * in->myDesc->incs / in->myDesc->inaxes[0];
    if (kstoke0 > 0) {
      /* true Stokes (need IQUV for this to work) */
      g_error ("FIX ME");
      in->mySel->bothCorr = TRUE;  /* Need both */
      in->jadr[ivpnt][0] = (1-kstoke0) * incs;
      in->jadr[ivpnt][1] = (4-kstoke0) * incs;
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 1.0;

    } else {
      /* XX */
      in->jadr[ivpnt][0] = 0;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 0.0;
    } 
} /* end ObitUVCalSelectXXInit */

/**
 * Initialize structures for data selection for YY.
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in      Object to initialize.
 * \param kstoke0 Stokes index of first correlation
 *                (1-4) => I,Q,U,V; 
 *                (-5 -> -8)=> XX, YY, XY, YX
 * \param nstok   Number of Stokes correlations
 * \param ivpnt   Index in data selection structures
 *                Updated on output.
 * \param err     ObitError stack.
 */
static void 
ObitUVCalSelectYYInit (ObitUVCal *in, olong kstoke0, olong nstok, 
		      olong ivpnt, ObitErr *err)
{
  gboolean missing;
  olong incs;

  /* error check */
  if (err->error) return;

    /* Add Stokes YY -  check if in data */
    missing = ((kstoke0>0) && ((nstok<4) || (kstoke0>1) || (kstoke0+nstok-1<4)));
    missing = missing ||
      ((kstoke0<0) && ((kstoke0<-6) || (kstoke0-nstok+1>-6)));
    if (missing) {
      Obit_log_error(err, OBIT_Error, 
		     "Stokes YY not available for %s", in->name);
      return;
    }
    /* true Stokes.(need IQUV for this to work) */

    /* Data will always be uncompressed by the time it gets here */
    incs = 3 * in->myDesc->incs / in->myDesc->inaxes[0];
    
    if (kstoke0 > 0) {
      g_error ("FIX ME"); /* Convert */
      in->mySel->bothCorr = TRUE;  /* Need both */
      in->jadr[ivpnt][0] = (1-kstoke0) * incs;
      in->jadr[ivpnt][1] = (4-kstoke0) * incs;
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = -1.0;
    } else {
      /* YY. */
      in->jadr[ivpnt][0] = (6+kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 0.0;
    } 
} /* end ObitUVCalSelectYYInit */

/**
 * Initialize structures for data selection for XY.
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in      Object to initialize.
 * \param kstoke0 Stokes index of first correlation
 *                (1-4) => I,Q,U,V; 
 *                (-5 -> -8)=> XX, YY, XY, YX
 * \param nstok   Number of Stokes correlations
 * \param ivpnt   Index in data selection structures
 *                Updated on output.
 * \param err     ObitError stack.
 */
static void 
ObitUVCalSelectXYInit (ObitUVCal *in, olong kstoke0, olong nstok, 
		      olong ivpnt, ObitErr *err)
{
  gboolean missing;
  olong incs;

  /* error check */
  if (err->error) return;

    /* Add Stokes LL -  check if in data */
    missing = ((kstoke0>0) && ((nstok<4) || (kstoke0>1) || (kstoke0+nstok-1<4)));
    missing = missing ||
      ((kstoke0<0) && ((kstoke0<-7) || (kstoke0-nstok+1>-7)));
    if (missing) {
      Obit_log_error(err, OBIT_Error, 
		     "Stokes XY not available for %s", in->name);
      return;
    }

    /* Data will always be uncompressed by the time it gets here */
    incs = 3 * in->myDesc->incs / in->myDesc->inaxes[0];
    
    if (kstoke0 > 0) {
      /* true Stokes (need QU for this to work)*/
      g_error ("FIX ME"); /* Convert */
      in->mySel->bothCorr = TRUE;  /* Need both */
      in->jadr[ivpnt][0] = (2-kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0] + incs;
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 1.0;
      /* q + iu code */
      in->jadr[ivpnt][1] = -in->jadr[ivpnt][1];
    } else {
      /* XY */
      in->jadr[ivpnt][0] = (7+kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 0.0;
    } 
} /* end ObitUVCalSelectXYInit */

/**
 * Initialize structures for data selection for YX.
 * Adapted from the AIPSish DGINIT.FOR.
 * \param in      Object to initialize.
 * \param kstoke0 Stokes index of first correlation
 *                (1-4) => I,Q,U,V; 
 *                (-5 -> -8)=> XX, YY, XY, YX
 * \param nstok   Number of Stokes correlations
 * \param ivpnt   Index in data selection structures
 *                Updated on output.
 * \param err     ObitError stack.
 */
static void 
ObitUVCalSelectYXInit (ObitUVCal *in, olong kstoke0, olong nstok, 
		      olong ivpnt, ObitErr *err)
{
  gboolean missing;
  olong incs;

  /* error check */
  if (err->error) return;

    /* Add Stokes YX -  check if in data */
    missing = ((kstoke0>0) && ((nstok<4) || (kstoke0>1) || (kstoke0+nstok-1<4)));
    missing = missing ||
      ((kstoke0<0) && ((kstoke0<-8) || (kstoke0-nstok+1>-8)));
    if (missing) {
      Obit_log_error(err, OBIT_Error, 
		     "Stokes YX not available for %s", in->name);
      return;
    }

    /* Data will always be uncompressed by the time it gets here */
    incs = 3 * in->myDesc->incs / in->myDesc->inaxes[0];
    
    if (kstoke0 > 0) {
      /* true Stokes (need QU for this to work) */
      g_error ("FIX ME"); /* COnvert */
      in->mySel->bothCorr = TRUE;  /* Need both */
      in->jadr[ivpnt][0] = (2-kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0] + incs;
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = -1.0;
      /* q + iu code */
      in->jadr[ivpnt][1] = -in->jadr[ivpnt][1];
    } else {
      /* YX. */
      in->jadr[ivpnt][0] = (8+kstoke0) * incs;
      in->jadr[ivpnt][1] = in->jadr[ivpnt][0];
      in->selFact[ivpnt][0] = 1.0;
      in->selFact[ivpnt][1] = 0.0;
    } 
} /* end ObitUVCalSelectYXInit */

