/* $Id: ObitGBTCCBUtil.c,v 1.1 2006/03/30 18:38:00 bcotton Exp $                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitOTF.h"
#include "ObitGBTCCBUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGBTCCBUtil.c
 * GBT CalTech Continuum Backend (CCB) utility routine definitions.
 */


/*----------------------Public functions---------------------------*/
/**
 * Given a set of cal measurements, average the valid ones and replace all
 * with these values.
 * If there is only one state in the data then it is beamswitched data and only
 * the abs values if taken of each value.
 * The CCB system injects cals into either RCP/X or LCP/Y; in additions a given 
 * port switches both feed and polarization.  
 * For each port average the valid cals and replace all with the average.
 * \param desc    Descriptor for data
 * \param calFlag Cal indicator, >0 => RCP/X and <0 => LCP/Y
 * \param ndetect Number of detectors (data values)
 * \param cal     Cal measurements, will be replaces by the port averages
 */
void ObitGBTCCBAvgCal (ObitOTFDesc *desc, ofloat calFlag, olong ndetect, ofloat *cal)
{
  olong i, ifreq, nfreq, ipoln, npoln, istate, nstate, ifeed, nfeed;
  olong jpoln, jfeed, incfreq, incpoln, incstate, incfeed;
  olong indx;
  ofloat sum, count, avg, fblank = ObitMagicF();

  if (ndetect<=0) return;
  
  /* How many of things */
  if (desc->jlocf>=0) {
    nfreq   = desc->inaxes[desc->jlocf];
    incfreq = desc->incf;
  } else {nfreq = 1;  incfreq = 1;}
  if (desc->jlocs>=0) {
    npoln   = desc->inaxes[desc->jlocs];
    incpoln = desc->incs;
  } else {npoln= 1;  incpoln = 1;}
  if (desc->jlocfeed>=0) {
    nfeed   = desc->inaxes[desc->jlocfeed];
    incfeed = desc->incfeed;
  } else {nfeed= 1;  incfeed = 1;}
  if (desc->jlocstate>=0) {
    nstate   = desc->inaxes[desc->jlocstate];
    incstate = desc->incstate;
  } else {nstate= 1; incstate = 1;}

  /* If only one state then data has been beamswitched - only take abs */
  if (nstate<=1) {
    for (i=0; i<ndetect; i++) cal[i] = fabs(cal[i]);
    return;
  }
  
  /* loop over frequency */
  for (ifreq=0; ifreq<nfreq; ifreq++) {
    
    /* (L1,R2) (L2,R1) are the switching pairs */
    for (jfeed=0; jfeed<nfeed; jfeed++) {
      for (jpoln=0; jpoln<npoln; jpoln++) {
	sum = count = 0.0;
	if (((calFlag>0.0) && (jpoln==0)) || ((calFlag<0.0) && (jpoln==1))) {
	  /* Cal in sig */      
	  istate = 0;
	  ifeed  = jfeed;
	  ipoln  = 0;
	  indx = npoln * (nfeed * (istate*nfreq + ifreq) + ifeed) + ipoln;
	  if (cal[indx]!=fblank) {sum += cal[indx]; count++;}
	  istate = 3;
	  indx = npoln * (nfeed * (istate*nfreq + ifreq) + ifeed) + ipoln;
	  if (cal[indx]!=fblank) {sum += cal[indx]; count++;}
	} else {
	  /* Cal in ref */      
	  istate = 1;
	  ifeed  = 1-jfeed;
	  ipoln  = 1;
	  indx = npoln * (nfeed * (istate*nfreq + ifreq) + ifeed) + ipoln;
	  if (cal[indx]!=fblank) {sum += cal[indx]; count++;}
	  istate = 2;
	  indx = npoln * (nfeed * (istate*nfreq + ifreq) + ifeed) + ipoln;
	  if (cal[indx]!=fblank) {sum += cal[indx]; count++;}
	}
	
	/* Average */
	if (count>0.0) avg = sum / count;
	else avg = fblank;
	
	/* Replace all for this port */
	istate = 0;  /* sig */
	ifeed  = jfeed;
	ipoln  = jpoln;
	indx = npoln * (nfeed * (istate*nfreq + ifreq) + ifeed) + ipoln;
	cal[indx] = avg;
	istate = 3;
	indx = npoln * (nfeed * (istate*nfreq + ifreq) + ifeed) + ipoln;
	cal[indx] = avg;
	istate = 1;    /* ref */
	ifeed  = 1-jfeed;
	ipoln  = 1-jpoln;
	indx = npoln * (nfeed * (istate*nfreq + ifreq) + ifeed) + ipoln;
	cal[indx] = avg;
	istate = 2;
	indx = npoln * (nfeed * (istate*nfreq + ifreq) + ifeed) + ipoln;
	cal[indx] = avg;
      } /* end loop over poln */
    } /* end loop over feed */
  } /* end loop over freq */
} /* end ObitGBTCCBAvgCal */

