# Use Obit image "CLEAN" window to edit source lists
# exec(open('CleanWindowEdit.py').read())
import Image, ImageDesc, Table, OWindow, ODisplay, OTObit, OErr
import pickle
GetWindowList=None; IsInWindow=None; SaveWindowList=None; FetchWindowList=None;
WindowList2Mask=None; reEditWindowList=None
err  = OErr.OErr()
disp = ODisplay.ODisplay("ObitView", "ObitView", err)

del GetWindowList
def GetWindowList(inIm, err, disp=disp):
    """
    Create window on an image via editing on display disp
    
    Adds image Descriptor.Dict to win.Desc
    returns the window list:
    list of [id, type(0=rect, 1=round), parameters[blc,trc] or [rad,cenx,ceny])
    * inIm  = image for window
    * err   = Obit message/error object
    * disp  = image display
    """
    win = OTObit.window(inIm)
    ODisplay.PImage(disp, inIm, err, win)
    winList = OWindow.PGetList(win, 1, err)
    return winList
# end GetWindowList

del IsInWindow
def IsInWindow(x, y, inIm, winList,err):
    """
    Determine if a position is in the window list

    Returns True or False
    * x       = RA or long (deg)
    * y       = Dec or lat (deg)
    * inIm    = image for window
    * winList = list from GetWindowList
    * err     = Obit message/error object
    """
    if y<-90.0:
        return True  # Already flagged
    try:
        pixel = ImageDesc.PGetPixel(inIm.Desc, [x,y],err)
    except Exception as exception:
        print (exception)
        print ("failed for position",x,y)
        return True  # Flag to be sure
    #print (pixel)
    # loop over winList
    for w in winList:
        #print (win)
        if w[1]==0: # Rectangle [(2,3),(4,5)]
            if (pixel[0]>=w[2]) and (pixel[0]<=w[4]) and (pixel[1]>=w[3]) and (pixel[1]<=w[5]):
                #print('rec',w)
                return True
        else:       # Circle (rad=2, x=3, y=4]
            del2 = (pixel[0]-w[3])**2 + (pixel[1]-w[4])**2 
            if del2<w[2]*w[2]:
                #print('cir',w)
                return True
    return False  # No match
# end IsInWindow

del SaveWindowList
def SaveWindowList(winList, pfile):
    """
    Save WindowList to pickle jar

    * winList = OWindow to save
    * pfile   = root name of pickle file- ".pickle" added
    """
    tfile = pfile+'.pickle'
    fd = open(tfile,"wb")
    pickle.dump(winList, fd)
    fd.close()
# end SaveWindowList

del FetchWindowList
def FetchWindowList(pfile):
    """
    Get WindowList from pickle jar

    returns WindowList
    * pfile   = root name of pickle file- ".pickle" added
    """
    tfile = pfile+'.pickle'
    fd = open(tfile,"rb")
    winList = pickle.load(fd)
    fd.close()
    return winList
# end FetchWindowList

del WindowList2Mask
def WindowList2Mask(inIm, winList, mfile, err):
    """
    Write a CLEAN mask file from a WindowList

    * inIm    = Image being used
    * winList = WindowList to write as mask file
    * mfile   = name of CLEAN mask file
    """
    fd = open(mfile,"a")
    for w in winList:
        if w[1]==0: # Rectangle [(2,3),(4,5)]
            pass # ignore for now
        else:       # Circle (rad=2, x=3, y=4]
            (ra,dec) = ImageDesc.PGetPos(inIm.Desc, [float(w[3]), float(w[4])], err)
            ras = ImageDesc.PRA2HMS(ra); decs = ImageDesc.PDec2DMS(dec)
            line = "%s %s %d\n"%(ras,decs,w[2])
            fd.write(line)
         # end loop over list
    fd.close()
# end WindowList2Mask

del reEditWindowList
def reEditWindowList(inIm, winList, err):
    """
    Convert a winList to a OWindow and edit

    returns new WindowList
    * inIm    = Image being used
    * winList = WindowList to edit
    * err     = Obit message/error object
    """
    win = OTObit.window(inIm)
    OWindow.PSetList(win, winList, 1, err)
    ODisplay.PImage(disp, inIm, err, win)
    winList = OWindow.PGetList(win, 1, err)
    return winList
# end reEditWindowList

# Sample usage for making a mask file
# inIm=getFITS('../Abell3395_100_IPol.fits',0)
# winList = GetWindowList(inIm, err)
# mfile='Abell3395_S.mask'
# WindowList2Mask(inIm, winList, mfile, err)
# pfile='Abell3395_S.WindowList'
# SaveWindowList(winList, pfile)
# winList2=reEditWindowList(inIm, winList, err)

# win = OTObit.window(inIm)
# ODisplay.PImage(disp, inIm, err, win)
# nwin=OWindow.PGetMaxID(win, 1)
# lst = OWindow.PGetList(win, 1, err)
# list of [id, type(0=rect, 1=round), parameters[blc,trc], [rad,cenx,ceny])
# Deleted???
