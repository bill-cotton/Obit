/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV data self calibration utilities     */
/*                                                                    */
/*;  Copyright (C) 2006-2008                                          */
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

%{
#include "ObitUVSoln.h"
#include "ObitTableSNUtil.h"
%}


%inline %{
// Utility routines from ObitUVSoln class
extern int UVSolnRefAnt (ObitTable *SNTab, int isuba, int refant, ObitErr* err)
{
   olong lisuba, lrefant;
   lisuba  = isuba;
   lrefant = isuba;
   ObitUVSolnRefAnt ((ObitTableSN*)SNTab, lisuba, &lrefant, err);
   return lrefant;
} // end UVSolnRefAnt

extern void UVSolnSNSmo (ObitTable *SNTab, int isuba, ObitErr* err)
{
   olong lisuba;
   lisuba = isuba;
   ObitUVSolnSNSmo ((ObitTableSN*)SNTab, lisuba, err);
} // end UVSelfSNSmo

extern void UVSolnDeselSN (ObitTable *SNTab, int isuba, int fqid, int nantf, 
                           int *ants, int nsou, int *sources,
                           float timerange[2], ObitErr* err)
{
  olong lisuba, lfqid, lnantf, lnsou, i;
  olong* lants, *lsources; 
  ofloat ltimerange[2];

  lisuba = isuba;
  lfqid  = fqid;
  lnantf = nantf;
  if (ants[0]==0) lnantf = 0;  // 0 => all
  ltimerange[0] = timerange[0]; ltimerange[1] = timerange[1];
  lants = g_malloc0(lnantf*sizeof(olong));
  for (i=0; i<nantf; i++) lants[i] = ants[i];
  lnsou = (olong)nsou;
  lsources = g_malloc0(lnsou*sizeof(olong));
  for (i=0; i<nsou; i++) lsources[i] = sources[i];
  ObitUVSolnDeselSN ((ObitTableSN*)SNTab, lisuba, lfqid, lnantf, lants, 
                      lnsou, lsources, ltimerange, err);
  g_free(lants);
  g_free(lsources);
} // end UVSolnDeselSN

extern void UVSolnDeselCL (ObitTable *CLTab, int isuba, int fqid, int nantf, 
                           int *ants, int nsou, int *sources,
                           float timerange[2], ObitErr* err)
{
  olong lisuba, lfqid, lnantf, lnsou, i;
  olong* lants, *lsources; 
  ofloat ltimerange[2];

  lisuba = isuba;
  lfqid  = fqid;
  lnantf = nantf;
  if (ants[0]==0) lnantf = 0;  // 0 => all
  ltimerange[0] = timerange[0]; ltimerange[1] = timerange[1];
  lants = g_malloc0(lnantf*sizeof(olong));
  for (i=0; i<nantf; i++) lants[i] = ants[i];
  lnsou = (olong)nsou;
  lsources = g_malloc0(lnsou*sizeof(olong));
  for (i=0; i<nsou; i++) lsources[i] = sources[i];
  ObitUVSolnDeselCL ((ObitTableCL*)CLTab, lisuba, lfqid, lnantf, lants, 
                     lnsou, lsources, ltimerange, err);
  g_free(lants);
  g_free(lsources);
} // end UVSolnDeselCL

extern ObitTable* SNInvert (ObitTable *inSN, ObitData *outData, long tabVer,
                            int doRepl, ObitErr* err)
{
  olong ltabVer = tabVer;
  gboolean ldoRepl = (gboolean)doRepl;
  ObitTable *outSN=NULL;

  outSN = (ObitTable*)ObitTableSNUtilInvert ((ObitTableSN*)inSN, outData, 
    &ltabVer, ldoRepl, err);

  return outSN;
} // end SNInvert


%}
