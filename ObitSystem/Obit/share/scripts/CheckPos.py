# Check if a pointing is in a Mosaic
import ImageDesc, SkyGeom, OErr
def CheckPos(p, g, cat, galactic, prtLv, err, maxd=0.60):
    """ 
    Check if a pointing p overlaps mosaic g
    
    No correction for corners.
    returns True if overlap, else False
    
    * p        = Pointing image name (e.g. "T1R01C05")
    * g        = Mosaic image
    * cat      = catalog of pointings as dict (Txxxxxxxx:(RA,Dec))
    * galactic = if True convert to Galactic
    * prtLv    = diagnostic print level >=2 -> some
    * err      = Obit error stack
    * maxd     = max distance from edge (deg), 
    """
    ################################################################
    (pgra,pgdec) = cat[p][0:2]                # Pointing position, J2000
    if galactic:
        (ra1950, dec1950) = SkyGeom.PJtoB(pgra, pgdec)
        (lon, lat)        = SkyGeom.PEq2Gal (ra1950, dec1950);
    else:
        lon = pgra; lat - pgdec
    (nx,ny)        = g.Desc.Dict['inaxes'][0:2] # mosaic size
    # Is it inside mosaic?  Grumble, Grumble, throws exception when too far
    try:
        (xp,yp) = ImageDesc.PGetPixel(g.Desc, [lon,lat], err)
    except Exception as exception:
        print (exception)
        print ("Bombed in CheckPos")
        err.Clear()
        return False
    if OErr.PIsErr(err):
        print( 'Check error',p,overlapx,overlapy, xp,yp)
        err.Clear()
        return False
    if (xp>0.0) and (xp<=nx) and (yp>0.0) and (yp<=ny):
        if prtLv>=2:
            print ('Check:',p, "inside accumulators")
        return True
    overlapx = (xp>0.0) and (xp<=nx)
    overlapy = (yp>0.0) and (yp<=ny)
    (gglong,gglat) = g.Desc.Dict['crval'][0:2]   # mosaic position
    (dx,dy)        = g.Desc.Dict['cdelt'][0:2]   # cell spacing
    # check distance from top 
    if abs(((yp-ny)*dy))<maxd:
        overlapy = True
    # bottom
    if abs(((-yp)*dy))<maxd:
        overlapy = True
    # left
    if abs(((-xp)*dx))<maxd:
        overlapx = True
    # right
    if abs(((xp-nx)*dx))<maxd:
        overlapx = True
    if prtLv>=2:
        print ('Check overlap',p,'X',overlapx,'Y',overlapy, xp,yp, dx,dy,maxd)
    return overlapx and overlapy
# end CheckMKGP
