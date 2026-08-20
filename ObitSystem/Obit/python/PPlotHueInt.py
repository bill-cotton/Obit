# Plot hue intensity image
#exec(open('PPlotHueInt.py').read())
PPlotHueInt=None
# Only if numpy, matplotlib are available

try:
    import numpy as np
    # astropy stuff
    try:
        from astropy.wcs import WCS
        from astropy.io import fits
        from astropy.visualization import ImageNormalize, AsinhStretch
        from astropy.visualization.wcsaxes import SphericalCircle
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        # matplotlib stuff
        try:
            import matplotlib
            import matplotlib.pyplot as plt
            import matplotlib.colors as mcolors

            # Stuff from ONumpyPlot
            from ONumpyPlot import PGetWCS, PGetImageNPArray, PGetSubImage

            del PPlotHueInt
            def PPlotHueInt (inInt, inHue, plotfile, err, \
                             color='gray', hueColor='rainbow', title=None, scale=1.0, \
                             vminI=None, vmaxI=None,  vminH=None, vmaxH=None, \
                             blc=[1,1,1,1], trc=[0,0,0,0], doAsinh=False, \
                             doColorBar=True, barLabelI='Jy/beam', barLabelH='Spectral Index ($\\alpha$)', \
                             barScale=1.0, xtick_min=4, ytick_deg=20,  dpi=100, fontsize=10):
                """
                Plot a hue/intensity image in a pdf file
                
                Much help from google on this.
                * inInt     = Intensity image as python ObitImage (or ObitImageMF) object
                * inHue     = Hue image python ObitImage (or ObitImageMF) object
                * plotfile  = root of plot file, ".pdf" added
                * err       = Obit error/message object
                * color     = intensity scheme "gray", "plasma", "inferno"
                              import matplotlib.pyplot as plt
                              see help(plt.colormaps), best with 'gray'
                * hueColor  = hue colormap, 'coolwarm' or 'rainbow' work well
                              "twilight" for cyclic hues
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
                * doColorBar= Show colorbars? 
                * barLabelI = label for intensity colorbar, accepts LaTex notation, $^{-2}$
                * barLabelH = label for hue colorbar
                * barScale  = scaling factor for colorbars
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
                cmap_hue = plt.get_cmap(hueColor).reversed()
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

                    cbar_hue = fig.colorbar(sm_hue, ax=ax, orientation='vertical', pad=0.03*barScale, shrink=0.6*barScale)
                    cbar_hue.set_label(barLabelH, fontsize=fontsize)
                    
                    # --- 6. Intensity Colorbar (Manual Non-Linear Tick Placement) ---
                    # Create a simple, neutral monochrome scale for the intensity reference
                    cmap_int = plt.get_cmap('gray')
                    norm_int = mcolors.Normalize(vmin=0, vmax=1) # The image array maps 0 to 1 under the hood
                    
                    sm_int = plt.cm.ScalarMappable(cmap=cmap_int, norm=norm_int)
                    sm_int.set_array([])
                    
                    cbar_int = fig.colorbar(sm_int, ax=ax, orientation='horizontal', pad=0.10*barScale, shrink=0.75*barScale)
                    if doAsinh:
                        cbar_int.set_label(barLabelI+" (asinh scale)", fontsize=fontsize)
                    else:
                        cbar_int.set_label(barLabelI, fontsize=fontsize)
                   
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
