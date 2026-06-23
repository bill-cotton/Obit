# Obit numpy access to data, plotting images
# $Id$
""" 
Utilities for numpy access to Obit image data, plotting using matplotlib, astropy

Includes wcs image axis labeling
Needs numpy for data access, also matplotlib, astropy for plotting
"""
#-----------------------------------------------------------------------
#  Copyright (C) 2026
#  Associated Universities, Inc. Washington DC, USA.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#  MA 02139, USA.
#
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------

# Python shadow class to ObitFArray class
from __future__ import absolute_import
from __future__ import print_function
import Obit, _Obit, InfoList, Image, FArray, OErr

def PGetSubImage (inImage, err, blc=[1,1,1,1], trc=[0,0,0,0]):
    """
    Return a memory resident image with a subsection of one plane of inImage
        
    * inImage   = ObitImage object
    * err       = Obit error/message object
    * blc       = Bottom left corner pixel (1-rel)
    * trc       = Top right corner pixel (1-rel), 0=> all
    * returns Memory resident ObitImage
    """
    ################################################################
    if ('myClass' in inImage.__dict__) and (inImage.myClass=='AIPSImage'):
        raise TypeError("Function unavailable for "+inImage.myClass)
    inCast = inImage.cast('ObitImage')
    # The simpler ways of doing this don't work
    inCast.List.set("BLC",blc)
    inCast.List.set('TRC',trc)
    inCast.FreeBuffer(err)
    z=inCast.FullInstantiate(Image.READONLY, err)
    tmp = Image.Image('subimage')
    Image.PCloneMem(inCast,tmp,err)
    z=tmp.FullInstantiate(Image.READWRITE, err)
    tmp.FreeBuffer(err)
    tmp.FArray = inCast.ReadPlane(err, blc=blc, trc=trc) 
    return tmp
# end PGetSubImage

# Only if numpy is available
try:
    import numpy as np
    
    def PGetImageNPArray (inImage, err, blc=[1,1,1,1], trc=[0,0,0,0]):
        """
        Return a numpy array for the data in the image buffer in inImage
        
        * inImage   = Python ObitImage (or ObitImageMF) object,
        * err       = Obit error/message object
        * blc       = Bottom left corner pixel (1-rel)
        * trc       = Top right corner pixel (1-rel), 0=> all
        * returns numpy array
        """
        ################################################################
        if ('myClass' in inImage.__dict__) and (inImage.myClass=='AIPSImage'):
            raise TypeError("Function unavailable for "+inImage.myClass)
        tmp = PGetSubImage (inImage, err, blc, trc)
        nx,ny=tmp.Desc.Dict['inaxes'][0:2]
        return np.frombuffer(tmp.PixBuf,dtype=np.float32).reshape(ny,nx,order='F')
    # end PGetImageNPArray
    
    def PSetImageNPArray (NPArr, inImage):
        """
        Copy a numpy array to the image buffer in inImage
        
        * NPArray   = input Numpy array
        * inImage   = Python ObitImage (or ObitImageMF) object
        """
        ################################################################
        if ('myClass' in inImage.__dict__) and (inImage.myClass=='AIPSImage'):
            raise TypeError("Function unavailable for "+inImage.myClass)
        if not isinstance(NPArr, np.ndarray):
            raise TypeError("First argument not an numpy ndarray ")
        # Check compatibility - NB: data in Fortran order
        inCast = inImage.cast('ObitImage')
        nx,ny = inCast.FArray.Naxis[0:2]
        if (NPArr.shape[0]!=ny) or (NPArr.shape[1]!=nx):
            raise RuntimeError("Incompatible sizes"+str((ny,nx))+" != "+str(NPArr.shape))
        # Copy
        np.copyto(NPArr, np.frombuffer(inCast.PixBuf,dtype=np.float32).reshape(ny,nx,order='F'))
        return 
    # end PSetImageNPArray

    # astropy stuff
    try:
        from astropy.wcs import WCS
        from astropy.io import fits
        from astropy.visualization import ImageNormalize, AsinhStretch

        def PGetWCS (inImage):
            """
            Create astropy wcs object for an image
            
            Useful for labeling plots, may not get image rotation correct
            * inImage   = Python Image object
            * returns  astropy wcs object 
            """
            ################################################################
            if ('myClass' in inImage.__dict__) and (inImage.myClass=='AIPSImage'):
                raise TypeError("Function unavailable for "+inImage.myClass)
            inCast = inImage.cast('ObitImage')
            d = inCast.Desc.Dict
            cards = []
            cards.append(fits.Card("SIMPLE","T","file does conform to FITS standard"))
            cards.append(fits.Card("BITPIX","-32","IEEE float"))
            cards.append(fits.Card("NAXIS",2,"Number of axes"))
            cards.append(fits.Card("NAXIS1",d['inaxes'][0],"Length of axis 1"))
            cards.append(fits.Card("NAXIS2",d['inaxes'][1],"Length of axis 2"))
            cards.append(fits.Card("CTYPE1",d['ctype'][0],"Type of axis 1"))
            cards.append(fits.Card("CTYPE2",d['ctype'][1],"Type of axis 2"))
            cards.append(fits.Card("CRPIX1",d['crpix'][0],"Reference pixel of axis 1"))
            cards.append(fits.Card("CRPIX2",d['crpix'][1],"Reference pixel of axis 2"))
            cards.append(fits.Card("CRVAL1",d['crval'][0],"Coordinate of axis 1"))
            cards.append(fits.Card("CRVAL2",d['crval'][1],"Coordinate of axis 2"))
            cards.append(fits.Card("CDELT1",d['cdelt'][0],"Coordinate increment of axis 1"))
            cards.append(fits.Card("CDELT2",d['cdelt'][1],"Coordinate increment of axis 2"))
            cards.append(fits.Card("CROTA1",d['crota'][0],"Coordinate rotation of axis 1"))
            cards.append(fits.Card("CROTA2",d['crota'][1],"Coordinate rotation of axis 2"))
            hh = fits.Header(cards)
            wcs = WCS(hh);
            return wcs
        # end PGetWCS

        # matplotlib stuff
        try:
            import matplotlib
            import matplotlib.pyplot as plt

            def PPlotImage (inImage, plotfile, err, \
                            color='gray', title=None, scale=1.0, vmin=None, vmax=None, \
                            blc=[1,1,1,1], trc=[0,0,0,0], doColorBar=False, doAsinh=False, \
                            barLocation='right', barLabel='Jy/beam', dpi=100, fontsize=10):
                """
                Plot an image in a pdf file
                
                * inImage   = Python ObitImage (or ObitImageMF) object
                * plotfile  = root of plot file, ".pdf" added
                * err       = Obit error/message object
                * color     = scheme "gray", "plasma", "inferno"
                              import matplotlib.pyplot as plt
                              see help(plt.colormaps)
                * title     = plot title, defaults to image object
                * scale     = Scale factor for image
                * vmin      = min pixel value, defaults to image min
                              after applying scale
                * vmax      = max pixel value, defaults to image max
                * blc       = Bottom left corner pixel (1-rel)
                * trc       = Top right corner pixel (1-rel), 0=> all
                * doColorBar= Show colorbar?
                * doAsinh   = Use Asinh (nonlinear) stretch?
                * barLocation = location of color bar "top", "right","lerft","bottom"
                * barLabel  = label for colorbar
                * dpi       = Output resolution in dots per inch
                * fontsize  = font size in points for labels
                """
                ################################################################
                if ('myClass' in inImage.__dict__) and (inImage.myClass=='AIPSImage'):
                    raise TypeError("Function unavailable for "+inImage.myClass)
                # Get subimage
                tmp = PGetSubImage (inImage, err, blc, trc)
                d = tmp.Desc.Dict
                nx,ny=d['inaxes'][0:2]; object = d['object']
                ttitle = title
                if not title:
                    ttitle = object
                ff=tmp.FArray;  FArray.PDeblank(ff, 0.0);  # Get pixel FArray, remove any blanks
                FArray.PSMul(ff,scale)  # Apply scale
                # Get max/min if needed
                vvmin = vmin; vvmax = vmax
                if not vmin:
                    pos = [0,0]
                    vvmin = FArray.PMin(ff,pos)
                if not vmax:
                    pos = [0,0]
                    vvmax = FArray.PMax(ff,pos)
                # Clip data
                FArray.PInClip(ff, -1.0e10, vvmin, vvmin)
                FArray.PInClip(ff, vvmax, 1.0e10, vvmax)
                # Get numpy array
                s=np.frombuffer(FArray.PGetBuf(ff),dtype=np.float32).reshape(nx,ny,order='F')
                ss=s.transpose()  #Get it right way around
                wcs =  PGetWCS(inImage)  # Get WCS info
                # Generate plot
                plt.rcParams.update({'font.size': fontsize})
                xsize=nx/dpi; ysize=ny/dpi
                fig = plt.figure(figsize=[xsize,ysize])
                ax = fig.add_subplot(111, projection=wcs)
                cblabel = barLabel
                if doAsinh:
                    # Create asinh normalization
                    norm = ImageNormalize(ss, stretch=AsinhStretch(a=0.1))
                    im=ax.imshow(ss, norm=norm,cmap=color,origin='lower')
                else:
                    im=ax.imshow(ss, cmap=color,origin='lower',vmin=vvmin,vmax=vvmax)
                labx = d['ctype'][0][0:4].replace('-','')+" (J2000)";
                laby = d['ctype'][1][0:4].replace('-','')+" (J2000)";
                z=plt.xlabel(labx); z=plt.ylabel(laby); z=plt.title(ttitle);
                if doColorBar:
                    z=plt.colorbar(im,label=cblabel,shrink=0.8,location=barLocation)
                matplotlib.pyplot.savefig(plotfile+".pdf",bbox_inches="tight",dpi=dpi)
                plt.close()  # Free resources
            # end  PPlotImage 
 
        except Exception as exception:
            print(exception)
            print ("Sorry, matplotlib unavailable")
    # end if astropy available
    except Exception as exception:
        print(exception)
        print ("Sorry, astropy unavailable")
    # end if astropy available
# end if numpy available
except Exception as exception:
    print(exception)
    print ("Sorry, numpy unavailable")


        
