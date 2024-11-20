# Parameters defining output mosaics
stktrans = 'I'          # Stokes or line transition name
galactic = False        # Output in Galactic coordinates?
size    = [1200,1200]   # size of master accumulations/mosaic
cells   = 1.5           # cellsize (asec)
# triplets of mosaic name root, center position (deg, deg) and 
#             beam (asec,asec,deg) or none
# Mosaic FITS files have names <mosaic name root>.<stktrans>_Mosaic.fits
# Example
targets =  [ \
             ('myMosaic_',[91.5, 19.0], [8.0,8.0,0.0]), \
]
#exec(open("Mosaic.py").read()) # to execute top level script 
