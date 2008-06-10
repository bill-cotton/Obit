/* $Id$  
/*--------------------------------------------------------------------*/
/* Swig module description for OTF Utilities                          */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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
#include "ObitFInterpolate.h"
#include "ObitOTFUtil.h"
#include "ObitOTFCalUtil.h"
%}



%inline %{

//extern void OTFUtilSubSkyModelBuff (ObitOTF *in, ObitOTFSkyModel *sky, float factor) {
//  ObitOTFUtilSubSkyModelBuff (in, sky, factor);
//} // end OTFUtilSubSkyModelBuff

extern void OTFUtilSubImage(ObitOTF *inOTF, ObitOTF *outOTF, ObitFArray *image, 
			 ObitImageDesc *desc, ObitErr *err) {
  ObitOTFUtilSubImage(inOTF, outOTF, image, desc, err);
} // end OTFUtilSubImage

extern void OTFUtilModelImage(ObitOTF *inOTF, ObitOTF *outOTF, ObitFArray *image, 
			 ObitImageDesc *desc, ObitErr *err) {
  ObitOTFUtilModelImage(inOTF, outOTF, image, desc, err);
} // end OTFUtilModelImage

extern void OTFUtilScale(ObitOTF *inOTF, ObitOTF *outOTF, float scale, float offset,
			 ObitErr *err) {
  ObitOTFUtilScale(inOTF, outOTF, scale, offset, err);
} // end OTFUtilScale

extern void OTFUtilNoise(ObitOTF *inOTF, ObitOTF *outOTF, float scale, float offset,
			 float sigma, ObitErr *err) {
  ObitOTFUtilNoise(inOTF, outOTF, scale, offset, sigma, err);
} // end OTFUtilNoise

extern void OTFUtilSubImageBuff (ObitOTF *in, ObitFInterpolate *image, float factor, 
                              ObitErr *err) {
  ObitOTFUtilSubImageBuff (in, image, factor, err);
} // end OTFUtilSubImageBuff

extern ObitImage* OTFUtilCreateImage (ObitOTF *inOTF, ObitErr *err) {
  return ObitOTFUtilCreateImage (inOTF, err);
} // end OTFUtilCreateImage

// doBeam = true if "Beam" desired
// doWt = true if "Weight image" desired
extern void OTFUtilMakeImage (ObitOTF *inOTF,  ObitImage *outImage, int doBeam, 
	ObitImage *beam, ObitImage *Wt, ObitErr *err) {
  gboolean ldoBeam = doBeam!=0;
  // May not be passed a beam, name of beam will be "NoBeam"
  ObitImage *lbeam=NULL, *lWt=NULL;
  if (beam && strncmp(beam->name,"NoBeam",6)) lbeam = beam;
  // May not be passed a Wt, name of Wt will be "NoWt"
  if (Wt && strncmp(Wt->name,"NoWt",4)) lWt = Wt;
  ObitOTFUtilMakeImage( inOTF, outImage, ldoBeam, lbeam, lWt, err);
} // end OTFUtilMakeImage

extern void OTFUtilIndex (ObitOTF *inOTF, ObitErr *err) {
  ObitOTFUtilIndex (inOTF, err);
} // end OTFUtilIndex

extern void OTFUtilFitCal (ObitOTF *inOTF, int detect, ObitErr *err) {
  ObitOTFCalUtilFitCal (inOTF, detect, err);
}  // end OTFUtilFitCal

extern void OTFUtilFitOnOff (ObitOTF *inOTF, int detect, ObitErr *err) {
  ObitOTFCalUtilFitOnOff (inOTF, detect, err);
}  // end OTFUtilFitOnOff

extern void OTFCalUtilFitTip (ObitOTF *inOTF, ObitErr *err) {
  ObitOTFCalUtilFitTip (inOTF, err);
}  // end OTFCalUtilFitTip

extern void OTFCalUtilFitNod (ObitOTF *inOTF, int detect, ObitErr *err) {
  ObitOTFCalUtilFitNod (inOTF, detect, err);
}  // end OTFCalUtilFitNod

extern void OTFUtilDiffNod (ObitOTF *inOTF, int scan, ObitErr *err) {
  ObitOTFUtilDiffNod (inOTF, scan, err);
}  // end OTFUtilDiffNod

extern int OTFCalUtilFlag (ObitOTF *inOTF, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitOTFCalUtilFlag (inOTF, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
}  // end OTFCalUtilFlag

extern ObitFArray* OTFUtilConvBeam (ObitTable *inTab, ObitImage *Beam, 
  ObitFArray *Template, ObitErr *err) {
  ObitFArray *out = NULL;
  ObitTableCC *CCTab=NULL;
  CCTab = ObitTableCCConvert(inTab);
  out = ObitOTFUtilConvBeam (CCTab, Beam, Template, err);
  CCTab = ObitTableCCUnref(CCTab);
  return out;
} // end OTFUtilConvBeam
%}
