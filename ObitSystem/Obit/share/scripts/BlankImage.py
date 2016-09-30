# interactively blank portions of an image
#x=imlod('Arcade_A_wide_3secHG.fits',0,'Arcade_A','3secHG',1,1,err)
#x=imlod('ArcadeAB3comb.fits',0,'Arcade_BC','Blank',1,1,err)
import OErr, Image, OWindow, ODisplay, FArray
BlankImage = None

# Define image
Aname = 'Arcade_A'; Aclass='3secHG'; disk = 1; seq = 1; offs = [0,0]
Aname = 'Arcade_BC'; Aclass='Blank'; disk = 1; seq = 1; offs = [-139,-169]
plane = [1,1,1,1,1]
image = Image.newPAImage("Im", Aname, Aclass, disk, seq, True, err)
OErr.printErr(err)

# Read image
image.GetPlane(None, plane, err)
OErr.printErr(err)

# Define blanking interactively defining a OWindow
naxis = image.Desc.Dict['inaxes'][0:2]
win = OWindow.PCreate1('wind', naxis, err)
ODisplay.PImage(disp, image, err,window=win)
OErr.printErr(err)
wlist = OWindow.PGetList(win, 1, err)
OErr.printErr(err)

# Apply
BlankImage(image,OWindow.PGetList(win,1,err),off=offs)
OErr.printErr(err)

# Update image
image.PutPlane( None, plane, err)
OErr.printErr(err)

# Blank Image function
del BlankImage
def BlankImage(image, win, off = [0,0]):
    """
    Blank portion of image described by win
    
    Blank rectangular or round regions specified in win
    image  Image to blank, FArray should be filled in
    win    List containing windows (OWindow.PGetList)
    off    offset in image
    """
    fblank = FArray.fblank
    nx = image.Desc.Dict['inaxes'][0]-1; ny = image.Desc.Dict['inaxes'][1]-1
    for w in win:
        if w[1]==0:  # Rectangular
            ix1 = w[2]-1; ix2 = w[4]; iy1 = w[3]-1; iy2 = w[5];
            ix1 = max(0,ix1); iy1 = max(0,iy1); ix2 = min(nx,ix2); iy2 = min(ny,iy2); 
            for ix in range(ix1,ix2):
                for iy in range(iy1,iy2):
                    image.FArray.set(fblank, ix+off[0], iy+off[1])
        elif w[1]==1:  # Round
            r = float(w[2]); r2 = r*r; xc = float(w[3]); yc = float(w[4])
            ix1 = w[3]-w[2]-1; ix2 = w[3]+w[2]; iy1 = w[4]-w[2]-1; iy2 = w[4]+w[2];
            ix1 = max(0,ix1); iy1 = max(0,iy1); ix2 = min(nx,ix2); iy2 = min(ny,iy2); 
            for ix in range(ix1,ix2):
                for iy in range(iy1,iy2):
                    d = (ix-xc)**2 + (iy-yc)**2
                    if d<=r2:
                        image.FArray.set(fblank, ix+off[0], iy+off[1])
    # end loop over win
# end BlankImage
