# Obit numpy access to data, plotting images
# $Id$
#exec(open('ONumpyPlot.py').read())
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
PGetSubImage=None; PGetImageNPArray=None;  PSetImageNPArray=None; PGetWCS=None;
PPlotImage=None; PPlotHueInt=None

del PGetSubImage
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
    
    del PGetImageNPArray
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
    
    del PSetImageNPArray
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
        from astropy.visualization.wcsaxes import SphericalCircle
        from astropy.coordinates import SkyCoord
        import astropy.units as u
 
        del PGetWCS
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
            import matplotlib.colors as mcolors

            del PPlotImage 
            def PPlotImage (inImage, plotfile, err, \
                            color='gray', title=None, scale=1.0, vmin=None, vmax=None, \
                            blc=[1,1,1,1], trc=[0,0,0,0], doColorBar=False, doAsinh=False, \
                            knee=None, barLocation='right', barLabel='Jy/beam', \
                            dpi=100, fontsize=10):
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
                * knee      = Asinh transition from linear to logarithmic
                              default 10% of the way from the minimum to maximum.
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
                # default knee for asinh
                if knee==None:
                    knee = vvmin + 0.10*(vvmax-vvmin)
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
                labx = d['ctype'][0][0:4].replace('-','')+" (J2000)";
                laby = d['ctype'][1][0:4].replace('-','')+" (J2000)";
                z=plt.xlabel(labx); z=plt.ylabel(laby); z=plt.title(ttitle);
                if doAsinh:
                    # Create asinh version
                    linear_width = knee/vvmax
                    stretched_data = np.arcsinh(ss / linear_width)
                    im=ax.imshow(stretched_data, cmap=color,origin='lower')
                    if doColorBar:
                        # 4. Set up the colorbar (google suggests)
                        cbar = plt.colorbar(im)
                        # 5. FIX THE LABELS: Define ideal tick marks based on your ORIGINAL data values
                        # Pick values that represent the scale of your real data
                        original_ticks = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300] 
                        
                        # Map those original values into the stretched space to find where they belong
                        stretched_ticks = np.arcsinh(np.array(original_ticks) / linear_width)
                        
                        # Apply the positions and the original strings to the colorbar
                        cbar.set_ticks(stretched_ticks)
                        cbar.set_ticklabels([str(x) for x in original_ticks])
                        cbar.set_label(cblabel+' (asinh stretch)')
                else:
                    im=ax.imshow(ss, cmap=color,origin='lower',vmin=vvmin,vmax=vvmax)
                    if doColorBar:
                        z=plt.colorbar(im,label=cblabel,shrink=0.8,location=barLocation)
                matplotlib.pyplot.savefig(plotfile+".pdf",bbox_inches="tight",dpi=dpi)
                plt.close()  # Free resources
            # end  PPlotImage
            
            del PPlotHueInt
            def PPlotHueInt (inInt, inHue, plotfile, err, \
                             color='gray', hcolor='rainbow', title=None, scale=1.0, \
                             vminI=None, vmaxI=None,  vminH=None, vmaxH=None, \
                             blc=[1,1,1,1], trc=[0,0,0,0], doAsinh=False, \
                             doColorBar=True, barLabelI='Jy/beam', barLabelH='Spectral Index ($\\alpha$)', \
                             xtick_min=4, ytick_deg=20,  dpi=100, fontsize=10):
                """
                Plot a non-cyclic hue/intensity image in a pdf file
                
                Much help from google on this.
                * inInt     = Intensity image as python ObitImage (or ObitImageMF) object
                * inHue     = Hue image python ObitImage (or ObitImageMF) object
                * plotfile  = root of plot file, ".pdf" added
                * err       = Obit error/message object
                * color     = intensity scheme "gray", "plasma", "inferno"; "gray" best
                * hcolor    = hue color scheme "twilight_shifted" for cyclic, e,g, EVPA
                              import matplotlib.pyplot as plt
                              see help(plt.colormaps)
                * title     = plot title, defaults to image object
                * scale     = Scale factor for inInt
                * vminI     = inInt min pixel value, defaults to image min
                              after applying scale
                * vmaxI     = inInt max pixel value, defaults to image max
                * vminH     = inHue min value, defaults to image min
                * vmaxH     = inHue max value, defaults to image max
                * blc       = Bottom left corner pixel (1-rel)
                * trc       = Top right corner pixel (1-rel), 0=> all
                * doAsinh   = Use Asinh (nonlinear) stretch for inInt ?
                * doColorBar= Show colorbar? for hue
                * barLocation = location of color bar "top", "right","left","bottom"
                * barLabelI = label for intensity colorbar
                * barLabelH = label for hue colorbar
                * xtick_min = number of x ticks per minute of time
                * ytick_deg = number of x ticks per degree.
                * dpi       = Output resolution in dots per inch
                * fontsize  = font size in points for labels
                """
                ################################################################
                if ('myClass' in inInt.__dict__) and (inInt.myClass=='AIPSImage'):
                    raise TypeError("Function unavailable for "+inInt.myClass)
                if ('myClass' in inHue.__dict__) and (inHue.myClass=='AIPSImage'):
                    raise TypeError("Function unavailable for "+inHue.myClass)
                # Get subimages
                tmpInt = PGetSubImage (inInt, err, blc, trc)
                dInt = tmpInt.Desc.Dict
                nx,ny=dInt['inaxes'][0:2]; object = dInt['object']
                tmpHue = PGetSubImage (inHue, err, blc, trc)
                dHue = tmpHue.Desc.Dict
                ttitle = title
                if not title:
                    ttitle = object
                # Process inInt
                ff=tmpInt.FArray;  FArray.PDeblank(ff, 0.0);  # Get pixel FArray, remove any blanks
                FArray.PSMul(ff,scale)  # Apply scale
                # Deblank inHue
                fff=tmpHue.FArray;  FArray.PDeblank(fff, 0.0); 
                # Get inInt max/min if needed
                Ivmin = vminI; Ivmax = vmaxI
                if vminI==None:
                    pos = [0,0]
                    Ivmin = FArray.PMin(ff,pos)
                if vmaxI==None:
                    pos = [0,0]
                    Ivmax = FArray.PMax(ff,pos)
                # Clip inInt data
                FArray.PInClip(ff, -1.0e10, Ivmin, Ivmin)
                FArray.PInClip(ff, Ivmax, 1.0e10, Ivmax)
                # Get inInt numpy array
                sInt=np.frombuffer(FArray.PGetBuf(ff),dtype=np.float32).reshape(nx,ny,order='F')
                IntArr=sInt.transpose()  # Get it right way around
                wcs =  PGetWCS(inInt)    # Get WCS info
                wcs.wcs.cunit = ['deg', 'deg']  # Force sumbitch
                wcs.wcs.pc = [[1.0, 0.0], [0.0, 1.0]] # DAMN
                # Get inHue max/min if needed
                vmin_hue = vminH; vmax_hue = vmaxH
                if vminH==None:
                    pos = [0,0]
                    vmin_hue = FArray.PMin(fff,pos)
                if vmaxH==None:
                    pos = [0,0]
                    vmax_hue = FArray.PMax(fff,pos)
                # Clip inHue data
                FArray.PInClip(fff, -1.0e10,  vmin_hue,  vmin_hue)
                FArray.PInClip(fff, vmax_hue, 1.0e10, vmax_hue)
                # Get inHue numpy array
                sHue=np.frombuffer(FArray.PGetBuf(fff),dtype=np.float32).reshape(nx,ny,order='F')
                HueArr=sHue.transpose()  # Get it right way around
                # Hue color mapping
                cmap_hue = plt.get_cmap(hcolor).reversed()
                norm_hue = mcolors.Normalize(vmin=vmin_hue, vmax=vmax_hue)
                rgb_colors = cmap_hue(norm_hue(HueArr))[:, :, :3]
                defective = False
                # Using asinh stretch?
                if doAsinh:
                    # Create asinh normalization -  avoid matplotlib version issues using numpy
                    # 'linear_width' defines the scale where the transition from linear to log happens.
                    # Set it close to your background noise standard deviation (e.g., 0.1).
                    # --- 3. Intensity Mapping (Asinh via NumPy Math) ---
                    linear_width = 0.1 
                    
                    # Perform the raw mathematical arcsinh scaling
                    asinh_flux = np.arcsinh(IntArr / linear_width)
                    vmin_asinh = np.arcsinh(Ivmin / linear_width)
                    vmax_asinh = np.arcsinh(Ivmax / linear_width)
                    
                    # Scale precisely between 0.0 and 1.0 for image composite blending
                    int_scaled = (asinh_flux - vmin_asinh) / (vmax_asinh - vmin_asinh)
                    int_scaled = np.clip(int_scaled, 0, 1)
                else:
                    # linear
                    norm_intensity = mcolors.Normalize(vmin=Ivmin, vmax=Ivmax)
                    int_scaled = np.clip(norm_intensity(IntArr), 0, 1)
                # Apply intensity weight to your RGB colors
                final_int = rgb_colors * int_scaled[:, :, np.newaxis]
                
                # Generate plot
                plt.rcParams.update({'font.size': fontsize})
                xsize=nx/dpi; ysize=ny/dpi
                plt.clf()  # Donno why
                fig = plt.figure(figsize=[xsize,ysize])
                ax = fig.add_subplot(111)
                ax.set_aspect('equal', adjustable='box')
                # Plot image first a simple one (otherwise blows axis labeliing)
                img = ax.imshow(final_int, origin='lower')
                # Force axis labeling, astropy can't handle the complex image
                ny, nx = IntArr.shape
                
                # 2. Get absolute sky coordinates at the boundaries of your image footprint
                corner_sky_low = wcs.pixel_to_world(0, 0)
                corner_sky_high = wcs.pixel_to_world(nx - 1, ny - 1)
                
                # Extract raw degrees
                ra_min_deg, ra_max_deg = min(corner_sky_low.ra.deg, corner_sky_high.ra.deg), max(corner_sky_low.ra.deg, corner_sky_high.ra.deg)
                dec_min_deg, dec_max_deg = min(corner_sky_low.dec.deg, corner_sky_high.dec.deg), max(corner_sky_low.dec.deg, corner_sky_high.dec.deg)
                
                # =========================================================================
                # 3. GENERATE INTEGRAL MINUTE/ARCMINUTE BOUNDARIES
                # =========================================================================
                # Convert RA degree footprints to clean minutes/4 of time integers
                #xtick_min = 4
                ra_min_min = np.floor(ra_min_deg * xtick_min*4)
                ra_max_min = np.ceil(ra_max_deg * xtick_min*4)
                # Generate clean integers (e.g., step by 15 seconds)
                clean_ra_mins = np.arange(ra_min_min, ra_max_min + 1, 1) 
                clean_ra_degs = clean_ra_mins / (xtick_min*4)
                # Dec ticks every 
                #ytick_deg = 20
                dec_min_arcmin = np.floor(dec_min_deg * ytick_deg)
                dec_max_arcmin = np.ceil(dec_max_deg * ytick_deg)
                clean_dec_arcmins = np.arange(dec_min_arcmin, dec_max_arcmin + 1, 1) 
                clean_dec_degs = clean_dec_arcmins / (ytick_deg)
                # 4. REVERSE MATH: Convert clean degrees to precise pixel positions on your plot
                # Keep the crossing coordinate centered (using CRVAL) to guarantee valid mapping lines
                ra_sky_ticks = SkyCoord(ra=clean_ra_degs * u.deg, dec=np.ones_like(clean_ra_degs) * wcs.wcs.crval[1] * u.deg, frame='icrs')
                dec_sky_ticks = SkyCoord(ra=np.ones_like(clean_dec_degs) * wcs.wcs.crval[0] * u.deg, dec=clean_dec_degs * u.deg, frame='icrs')
                
                x_pixel_positions, _ = wcs.world_to_pixel(ra_sky_ticks)
                _, y_pixel_positions = wcs.world_to_pixel(dec_sky_ticks)
                
                # 5. FILTER: Only keep the ticks that actually fall inside your image frame
                valid_x_mask = (x_pixel_positions >= 0) & (x_pixel_positions < nx)
                x_ticks = x_pixel_positions[valid_x_mask]
                ra_labels = ra_sky_ticks[valid_x_mask].ra.to_string(unit='hour', sep=('h', 'm', 's'), fields=3, precision=0, pad=True)
                
                valid_y_mask = (y_pixel_positions >= 0) & (y_pixel_positions < ny)
                y_ticks = y_pixel_positions[valid_y_mask]
                dec_labels = dec_sky_ticks[valid_y_mask].dec.to_string(unit='deg', sep=('°', "'"), fields=2, precision=0, pad=True)
                
                # 6. INJECT PERFECTLY CLEAN COORDINATES
                ax.set_xticks(x_ticks)
                #ax.set_xticklabels(ra_labels, rotation=15, ha='right')
                ax.set_xticklabels(ra_labels, rotation=0, ha='center') # Labels can be too big
                
                ax.set_yticks(y_ticks)
                ax.set_yticklabels(dec_labels)

                ax.set_aspect('equal')

                # Labeling
                labx = dInt['ctype'][0][0:4].replace('-','')+" (J2000)";
                laby = dInt['ctype'][1][0:4].replace('-','')+" (J2000)";
                z=plt.xlabel(labx); z=plt.ylabel(laby); z=plt.title(ttitle);
                
                if doColorBar:
                    # --- 5. Hue Colorbar (Explicit Bound Enforcement) ---
                    sm_hue = plt.cm.ScalarMappable(cmap=cmap_hue, norm=norm_hue)
                    sm_hue.set_array([]) 

                    #cbar_hue = fig.colorbar(sm_hue, ax=ax, orientation='vertical', pad=0.03, shrink=0.8)
                    cbar_hue = fig.colorbar(sm_hue, ax=ax, orientation='vertical', pad=0.03, shrink=0.6)
                    cbar_hue.set_label(barLabelH, fontsize=11)
                    
                    # --- 6. Intensity Colorbar (Manual Non-Linear Tick Placement) ---
                    # Create a simple, neutral monochrome scale for the intensity reference
                    cmap_int = plt.get_cmap('gray')
                    norm_int = mcolors.Normalize(vmin=0, vmax=1) # The image array maps 0 to 1 under the hood
                    
                    sm_int = plt.cm.ScalarMappable(cmap=cmap_int, norm=norm_int)
                    sm_int.set_array([])
                    
                    cbar_int = fig.colorbar(sm_int, ax=ax, orientation='horizontal', pad=0.10, shrink=0.75)
                    #cbar_int = fig.colorbar(sm_int, ax=ax, orientation='horizontal', pad=0.15, shrink=0.8)
                    if doAsinh:
                        cbar_int.set_label(barLabelI+" (asinh scale)", fontsize=fontsize)
                    else:
                        cbar_int.set_label(barLabelI+" fract. of max.", fontsize=fontsize)
                   
                    # DEFINE CHOSEN TICK LOCATIONS (In your original physical data units)
                    if doAsinh:
                        physical_ticks = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, \
                                          50.0, 100.0, 200.0, 300.0, 500.0, 1000.0]
                    else:
                        # Linear
                        physical_ticks = [-0.1, 0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
                    
                    # Convert those physical thresholds into the corresponding 0-1 colorbar positions
                    colorbar_tick_positions = []
                    for val in physical_ticks:
                        if doAsinh:
                            # Asinh
                            asinh_val = np.arcsinh(val / linear_width)
                            scaled_val = (asinh_val - vmin_asinh) / (vmax_asinh - vmin_asinh)
                            colorbar_tick_positions.append(scaled_val)
                        else:
                             colorbar_tick_positions.append(val)
                    # Update the intensity colorbar with mapped ticks and raw text labels
                    cbar_int.set_ticks(colorbar_tick_positions)
                    cbar_int.set_ticklabels([str(t) for t in physical_ticks])
                    # end colorbars
                plt.tight_layout()  # tighter layout
                matplotlib.pyplot.savefig(plotfile+".pdf",bbox_inches="tight",dpi=dpi)
                print ("Plotted ",plotfile+".pdf")
                plt.close()  # Free resources
            # end  PPlotHueInt 
 
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


        
