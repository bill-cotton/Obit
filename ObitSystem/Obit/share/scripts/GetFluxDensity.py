# Script for getting the flux density of an extended source by 
# blanking everything else in a FITS image
file_root = 'J1234+5678'  # Root of file name, file_root+".fits"

pfile = file_root+'_Window.pickle' # Pickle file for saved window
plane = 1                          # Plane to blank
import Image, FArray, OErr, Blank

img = Image.newPFImage('input', file_root+'.fits',0,True,err)
blank = Blank.Blank("blank",img,err)
OErr.printErr(err)
# Set mask
blank.set_mask(disp,err)
# Save 
blank.save_window(pfile,err)
OErr.printErr(err)
print ('Wrote pickle file',pfile)
# If need to reedit 
#blank.fetch_window(pfile,err)
#blank.set_mask(disp,err)
# Make mask
blank.make_mask(err)
# Apply mask
blank.blank(img, err, plane=plane)
OErr.printErr(err)
# Read image into FArray
img.GetPlane(None, [plane,1,1,1,1], err)
# Save masked image
mimg = Image.PFArray2FITS(img.FArray, file_root+'_masked.fits',err, \
                          outDisk=0, oDesc=img.Desc)
OErr.printErr(err)
print ('Wrote image',file_root+'_masked.fits')
# Do integral
imstat(mimg)

