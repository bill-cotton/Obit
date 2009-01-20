/* $Id$  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Image utilities                        */
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
#include "ObitImageUtil.h"
#include "ObitTable.h"
#include "ObitTableCCUtil.h"
%}



%inline %{
ObitImage* ImageUtilCreateImage (ObitUV *inUV, long fieldNo,
			       int doBeam, ObitErr *err) {
  return ObitImageUtilCreateImage (inUV, fieldNo, doBeam, err);
} // end ImageUtilCreateImag

void ImageUtilMakeImage (ObitUV *inUV, ObitImage *outImage, 
			     long channel, int doBeam, 
			     int doWeight, ObitErr *err) {
  ObitImageUtilMakeImage (inUV, outImage, channel, doBeam, doWeight, err);
} // end ImageUtilMakeImag

void 
ImageUtilInterpolateImage (ObitImage *inImage, ObitImage *outImage, 
			       int *inPlane, int *outPlane,
			       long hwidth, ObitErr *err)
{
  ObitImageUtilInterpolateImage (inImage, outImage, inPlane, outPlane,
			         hwidth, err);
} // end ImageUtilInterpolateImage

void 
ImageUtilPBApply (ObitImage *inImage, ObitImage *pntImage, ObitImage *outImage, 
                      int *inPlane, int *outPlane, float antSize, ObitErr *err)
{
  ObitImageUtilPBApply (inImage, pntImage, outImage, inPlane, outPlane, antSize, err);
} // end ImageUtilPBApply

void 
ImageUtilPBImage (ObitImage *pntImage, ObitImage *outImage, 
                  int *outPlane, float antSize, float minGain, ObitErr *err)
{
  ObitImageUtilPBImage (pntImage, outImage, outPlane, antSize, minGain, err);
} // end ImageUtilPBImage

void 
ImageUtilPBCorr (ObitImage *inImage, ObitImage *pntImage, ObitImage *outImage, 
                      int *inPlane, int *outPlane, float antSize, ObitErr *err)
{
  ObitImageUtilPBCorr (inImage, pntImage, outImage, inPlane, outPlane, antSize, err);
} // end ImageUtilPBCorr

ObitImage * 
ImageUtilQuanFITS(ObitImage *inImage, char *fileName, 
		  int disk, ObitErr *err)
{
  return ObitImageUtilQuanFITS(inImage, (gchar*)fileName, (olong)disk, err);
} // end ImageUtilPBCorr

void 
ImageUtilCCScale (ObitTable *inCCTab, int startComp, int endComp, double scale, 
	          ObitErr *err)
{
  olong lstartComp, lendComp;
  ofloat lscale;
  ObitTableCC *CCTab=NULL;

  if (err->error) return;  // error check 

  lstartComp = startComp;
  lendComp = endComp;
  lscale   = (ofloat)scale;
  CCTab = ObitTableCCConvert(inCCTab);
  ObitTableCCUtilScale (CCTab, lstartComp, lendComp, lscale, err);
  CCTab = ObitTableCCUnref(CCTab);
} // end ImageUtilCCScale

void 
ImageUtilUVFilter (ObitImage *inImage, ObitImage *outImage, double radius, 
	          ObitErr *err)
{
  ofloat lradius;

  if (err->error) return;  // error check 

  lradius   = (ofloat)radius;
  ObitImageUtilUVFilter (inImage, outImage, lradius, err);
} // end ImageUtilUVFilter

%}
