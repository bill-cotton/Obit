""" Obit Plotting class

Create a plot object using newOPlot which allows specifying the output
and background color.  If no output is specified this information
will be prompted.
   Next, the plotting region must be specified using either PSetPlot,
one of the XY plotting routines (PXYPlot, PXYOver, or PXYErr) or
PContour.  Then additional lines, curves, text or symbols may be added.
When all has been added to the plot, use PShow to finalize it.

   Notes: on text strings in PLPlot installations

   If the Obit installation uses PLPlot for plotting the following
can be used in text strings:
 - Greek letters, A #g immediately prior to a Latin character will cause
   the Greek equivalent to be used, e.g. #ga will be a lower case alpha.
 - Subscripts: Characters between a #d and #u will be written as subscripts
 - Superscripts: Characters between a #u and #d will be written as
   superscripts
"""
# $Id: OPlot.py,v 1.5 2008/05/15 21:01:59 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C) 2006,2008
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

# Python shadow class to ObitPlot class
import Obit, InfoList, Image
import math

class OPlotPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.OPlotUnref(Obit.OPlot_me_get(self.this))
            # In with the new
            Obit.OPlot_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != OPlot:
            return
        if name == "me" : 
            return Obit.OPlot_me_get(self.this)
        # Functions to return members
        if name=="List":
            return PGetList(self)
        raise AttributeError,name
    def __repr__(self):
        if self.__class__ != OPlot:
            return
        return "<C OPlot instance> " + Obit.OPlotGetName(self.me)

class OPlot(OPlotPtr):
    """ Python Obit interface to display server
    
    This class is for creating and using the interface to a plot
    Image Members with python interfaces:
    InfoList  - used to pass instructions to processing
                Member List 
    """
    def __init__(self, name):
        self.this = Obit.new_OPlot(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_OPlot(self.this)

# Foreground Colors
BLACK       = 0
RED         = 1
YELLOW      = 2
GREEN       = 3
AQUAMARINE  = 4
PINK        = 5
WHEAT       = 6
GRAY        = 7
BROWN       = 8
BLUE        = 9
BLUEVIOLET  = 10
CYAN        = 11
TURQUOISE   = 12
MAGENTA     = 13
SALMON      = 14
WHITE       = 15

def newOPlot(name, err, output="None", bgcolor=BLACK, nx=1, ny=1 ):
    """ Create and initialize an ObitPlot

    name     = name desired for object (labeling purposes)
    err      = Python Obit Error/message stack
    output   = name and type of output device:
               "None"  interactive prompt
               "xwin"  X-Window (Xlib)
               "gcw"   Gnome Canvas Widget (interacts with ObitTalk)
               "ps"    PostScript File (monochrome)
               "psc"   PostScript File (color)
               "xfig"  Fig file
               "png"   PNG file
               "jpeg"  JPEG file
               "gif"   GIF file
               "null"  Null device
    bgcolor   = background color index (1-15), symbolic names:
                BLACK, RED(default), YELLOW, GREEN, 
                AQUAMARINE, PINK, WHEAT, GRAY, BROWN,
                BLUE, BLUEVIOLET, CYAN, TURQUOISE,
                MAGENTA, SALMON, WHITE
    nx        = Number of horizontal subpages
    ny        = Number of vertical subpages
    """
    ################################################################
    out = OPlot(name)
    Obit.PlotInitPlot(out.me, output, bgcolor, nx, ny, err.me)
    return out 
    # end newOPlot

def PXYPlot (plot, symbol, x, y, err):
    """ Simple XY Plot

    Plot X vs Y using symbol.
    Plot should be finalized and displayed with PShow
    This routine draws the frame and adds labels, to only overplot data
    on the same frame, use ObitPlotXYOver

    plot    = plot
    symbol  = Symbol index to use for plotting
              values in the range [1,12] are usable 
              if negative, use abs value and connect points
        0 = line only
        1 = dot
        2 = plus
        3 = *
        4 = open circle
        5 = x
        6 = open square
        7 = open triangle
        8 = open star
        9 = filled triangle
        10 = filled square
        11 = filled circle
        12 = filled star
    x   =    Independent variable, if None use index
    y   =    Dependent variable
    err =    ObitErr error stack
    
    Optional parameters on plot InfoList
    XMAX (float) maximum X value (defaults to actual value)
    XMIN (float) minimum X value (defaults to actual value)
    YMAX (float) maximum Y value (defaults to actual value)
    YMIN (float) minimum Y value (defaults to actual value)
    TITLE (string)  Label for the plot (defaults to none), max 120
    XLABEL (string) Label for horizontal axis (defaults to none)
    XOPT   (string) Options for horizontal axis (default "BCNTS")
          See PDrawAxes for details.
    YLABEL (string) Label for vertical axis (defaults to none)
    YOPT   (string) Options for  vertical axis (default "BCNTS")
          See PDrawAxes for details.
    XTICK (float) world coordinate interval between major tick marks
          on X axis. If xtick=0.0 [def], the interval is chosen.
    NXSUB (int) the number of subintervals to divide the major
          coordinate interval into. If xtick=0.0 or nxsub=0,
          the number is chosen. [def 0]
    YTICK  (float)  like xtick for the Y axis.
    NYSUB  (int)    like nxsub for the Y axis
    CSIZE  (int)    Scaling factor for characters(default = 1)
    SSIZE  (int)    Scaling factor for symbols(default = 1)
    LWIDTH (int)    Line width (default = 1)
    JUST   (int)    If !=0 then force X and Y axis scaling to be the same
    """
    ################################################################
    # Checks
    if not PIsA(plot):
        print "Actually ",plot.__class__
        raise TypeError,"plot MUST be a Python Obit Plot"
    n = len(y)  # How many points?
    Obit.PlotXYPlot (plot.me, symbol, n, x, y, err.me)
    # end PXYPlot


def PXYOver (plot, symbol, x, y, err):
    """ Overplot X vs Y 

    Overplot X vs Y using symbol.
    Plot should be finalized and displayed with PShow

    plot    = plot
    symbol  = Symbol index to use for plotting
              values in the range [1,12] are usable 
              if negative, use abs value and connect points
        0 = line only
        1 = dot
        2 = plus
        3 = *
        4 = open circle
        5 = x
        6 = open square
        7 = open triangle
        8 = open star
        9 = filled triangle
        10 = filled square
        11 = filled circle
        12 = filled star
    x   =    Independent variable, if None use index
    y   =    Dependent variable
    err =    ObitErr error stack
    """
    ################################################################
    # Checks
    if not PIsA(plot):
        print "Actually ",plot.__class__
        raise TypeError,"plot MUST be a Python Obit Plot"
    n = len(y)  # How many points?
    Obit.PlotXYOver (plot.me, symbol, n, x, y, err.me)
    # end PXYOver

def PXYErr (plot, symbol, x, y, e, err):
    """ Simple XY Plot with error bars

    Plot X vs Y using symbol and error bars.
    Plot should be finalized and displayed with PShow
    This routine draws the frame and adds labels, to only overplot data
    on the same frame, use ObitPlotXYOver

    plot     = plot
    symbol  = Symbol index to use for plotting
              values in the range [1,12] are usable 
              if negative, use abs value and connect points
        1 = dot
        2 = plus
        3 = *
        4 = open circle
        5 = x
        6 = open square
        7 = open triangle
        8 = open star
        9 = filled triangle
        10 = filled square
        11 = filled circle
        12 = filled star
    x   =    Independent variable, if None use index
    y   =    Dependent variable
    e   =    if nonNone, error in y
    err =    ObitErr error stack
    
    Optional parameters on plot InfoList
    XMAX (float) maximum X value (defaults to actual value)
    XMIN (float) minimum X value (defaults to actual value)
    YMAX (float) maximum Y value (defaults to actual value)
    YMIN (float) minimum Y value (defaults to actual value)
    TITLE (string)  Label for the plot (defaults to none), max 120
    XLABEL (string) Label for horizontal axis (defaults to none)
    XOPT   (string) Options for horizontal axis (default "BCNTS")
          See PDrawAxes for details.
    YLABEL (string) Label for vertical axis (defaults to none)
    YOPT   (string) Options for  vertical axis (default "BCNTS")
         See PDrawAxes for details.
    XTICK (float) world coordinate interval between major tick marks
          on X axis. If xtick=0.0 [def], the interval is chosen.
    NXSUB (int) the number of subintervals to divide the major
          coordinate interval into. If xtick=0.0 or nxsub=0,
          the number is chosen. [def 0]
    YTICK  (float)  like xtick for the Y axis.
    NYSUB  (int)    like nxsub for the Y axis
    CSIZE  (int)    Scaling factor for characters(default = 1)
    SSIZE  (int)    Scaling factor for symbols(default = 1)
    LWIDTH (int)    Line width (default = 1)
    JUST   (int)    If !=0 then force X and Y axis scaling to be the same
    """
    ################################################################
    # Checks
    if not PIsA(plot):
        print "Actually ",plot.__class__
        raise TypeError,"plot MUST be a Python Obit Plot"
    n = len(y)  # How many points?
    Obit.PlotXYErr (plot.me, symbol, n, x, y, e, err.me)
    # end PXYErr


def PContour (plot, label, image, lev, cntfac, err):
    """ Contour plot of image

    Contours at lev times powers of cntfac
    Plot should be finalized and displayed with PShow
    plot     = plot
    label    = Label for plot
    image    = ObitImage to plot
    lev      = basic contour level (def 0.1 peak)
    cntfac   = factor for spacing between contours (def sqrt(2)
    err      =    ObitErr error stack
    
    Optional parameters on plot InfoList
    XTICK (float) world coordinate interval between major tick marks
          on X axis. If xtick=0.0 [def], the interval is chosen.
    NXSUB (int) the number of subintervals to divide the major
          coordinate interval into. If xtick=0.0 or nxsub=0,
          the number is chosen. [def 0]
    YTICK  (float)  like xtick for the Y axis.
    NYSUB  (int)    like nxsub for the Y axis
    CSIZE  (int)    Scaling factor for characters(default = 1)
    LWIDTH (int)    Line width (default = 1)
    """
    ################################################################
    # Checks
    if not PIsA(plot):
        print "Actually ",plot.__class__
        raise TypeError,"plot MUST be a Python Obit Plot"
    if not Image.PIsA(image):
        print "Actually ",image.__class__
        raise TypeError,"image MUST be a Python Obit Image"
    Obit.PlotContour (plot.me, label, image.me, lev, cntfac, err.me)
    # end PContour


def PMarkCross (plot, image, ra, dec, err, size=5.0):
    """ Mark positions on Contour plot of image

    Place cross at positions.
    Plot should be finalized and displayed with PShow
    plot     = plot
    image    = ObitImage to plot
    ra       = list of RAs (deg)
    dec      = list of Declinations (deg)
    err      = ObitErr error stack
    size     = size of cross in pixels
    
    Optional parameters on plot InfoList
    CSIZE  (int)    Scaling factor for characters(default = 1)
    LWIDTH (int)    Line width (default = 1)
    """
    ################################################################
    # Checks
    if not PIsA(plot):
        print "Actually ",plot.__class__
        raise TypeError,"plot MUST be a Python Obit Plot"
    if not Image.PIsA(image):
        print "Actually ",image.__class__
        raise TypeError,"image MUST be a Python Obit Image"
    n = len(ra)
    Obit.PlotMarkCross (plot.me, image.me, n, ra, dec, size, err.me)
    # end PMarkCross

def PShow (plot, err):
    """ Display plot

    plot   = Python Plot object
    err =    ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    Obit.PlotFinishPlot(plot.me, err.me)
    # end  PShow
    
def PSetPlot (plot, xmin, xmax, ymin, ymax, just, axis, err):
    """ Define plotting area

    plot   = Python Plot object
    xmin   = the world x-coordinate at the bottom left corner of the viewport.
    xmax   = the world x-coordinate at the top right corner of the viewport 
             (note XMAX may be less than XMIN).
    ymin   = the world y-coordinate at the bottom left corner
             of the viewport.
    ymax   = the world y-coordinate at the top right corner
             of the viewport (note YMAX may be less than YMIN)
    just   = if JUST=1, the scales of the x and y axes (in
             world coordinates per inch) will be equal,
             otherwise they will be scaled independently.
    axis   = controls the plotting of axes, tick marks, etc:
         axis = -2 : draw no box, axes or labels;
         axis = -1 : draw box only;
         axis =  0 : draw box and label it with coordinates;
         axis =  1 : same as axis=0, but also draw the
                     coordinate axes (X=0, Y=0);
         axis =  2 : same as axis=1, but also draw grid lines
                     at major increments of the coordinates;
         axis = 10 : draw box and label X-axis logarithmically;
         axis = 20 : draw box and label Y-axis logarithmically;
         axis = 30 : draw box and label both axes logarithmically.
   err =    ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    Obit.PlotSetPlot(plot.me, xmin, xmax, ymin, ymax, just, axis, err.me)
    # end  PSetPlot

def PLabel (plot, xlabel, ylabel, title, err):
    """ Display plot

    plot   = Python Plot object
    xlabel  =  a label for the x-axis (centered below the  viewport).
    ylabel  =  a label for the y-axis (centered to the left
               of the viewport, drawn vertically)
    title   =  a label for the entire plot (centered above the viewport)
    err     =  ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    Obit.PlotLabel(plot.me, xlabel, ylabel, title, err.me)
    # end  PLabel

def PDrawAxes(plot, xopt, xtick, nxsub, yopt, ytick, nysub, err):
    """  Draw axes for a plot, label

    plot   = Python Plot object
    xopt   = string of options for X (horizontal) axis of plot.
             Options are single letters, and may be in any order 
             Axis options:
               A : draw Axis (X axis is horizontal line Y=0, Y axis is vertical line X=0).
               B : draw bottom (X) or left (Y) edge of frame.
               C : draw top (X) or right (Y) edge of frame.
               G : draw Grid of vertical (X) or horizontal (Y) lines
               I : Invert the tick marks; ie draw them outside the viewport instead of inside.
               L : label axis Logarithmically
               N : write Numeric labels in the conventional location below the
                   viewport (X) or to the left of the viewport (Y).
               M : write numeric labels in the unconventional location above the
                   viewport (X) or to the right of the viewport (Y).
               P : extend ("Project") major tick marks outside the box (ignored if
                   option I is specified)
               T : draw major Tick marks at the major coordinate interval.
               S : draw minor tick marks (Subticks).
    xtick =  World coordinate interval between major tick marks
             on X axis. If xtick=0.0, the interval is chosen.
    nxsub =  The number of subintervals to divide the major coordinate interval
             into. If xtick=0.0 or nxsub=0, the number is chosen.
    yopt  =  string of options for Y (vertical) axis of plot.
             Coding is the same as for xopt.
    ytick =  like xtick for the Y axis.
    nysub =  like nxsub for the Y axis
    err =    ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    Obit.PlotDrawAxes(plot.me, xopt, xtick, nxsub, yopt, ytick, nysub, err.me)
    # end  DrawAxes

def PSetCharSize (plot,cscale, err):
    """ Set scaling for characters

    plot   = Python Plot object
    cscale =  new character size (integer multiple of the default size).
    err    =  ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    Obit.PlotSetCharSize (plot.me, cscale, err.me)
    # end  PSetCharSize 

def PSetLineWidth (plot, lwidth, err):
    """ Set line width

    plot   = Python Plot object
    lwidth = Width of line (integer multiple of the default size).
    err    =    ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    Obit.PlotSetLineWidth(plot.me, lwidth, err.me)
    # end  PetLineWidth

def PSetColor (plot, color, err):
    """ Set foreground color

    plot   = Python Plot object
    color  =  color index (1-15), symbolic names
              BLACK, RED(default), YELLOW, GREEN, 
              AQUAMARINE, PINK, WHEAT, GRAY, BROWN,
              BLUE, BLUEVIOLET, CYAN, TURQUOISE,
              MAGENTA, SALMON, WHITE
   err     =  ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    Obit.PlotSetColor(plot.me, color, err.me)
    # end  PSetColor

def PSetPage (plot, sub, err):
    """ Set or advance sub page

    Note: some functions such as PContour advance the page
    plot   = Python Plot object
    sub    = if <=0 advance page, if >0 set current subpage to sub
             numbering starts at the top left at 1 and increases along
             rows and columns
    err    = ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    Obit.PlotSetPage(plot.me, sub, err.me)
    # end  PSetPage

def PText (plot, x, y, angle, just, text, err):
    """ Write text on plot

    plot    = Python Plot object
    x       = Plot x in world coordinates
    y       = Plot y in world coordinates
    angle   = Orientation of the text in deg, 0=horizontal
    just    = Controls justification of the string parallel to
              the specified edge of the viewport. If
              FJUST = 0.0, the left-hand end of the string will
              be placed at (x,y); if JUST = 0.5, the center of
              the string will be placed at (x,y); if JUST = 1.0,
              the right-hand end of the string will be placed at
              at (x,y). Other values between 0 and 1 give intermediate
              placing, but they are not very useful.
    text    = The text string to be plotted. Trailing spaces are
              ignored when justifying the string, but leading spaces are
              significant.
    err =    ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    dx = math.cos(angle/57.296)
    dy = math.sin(angle/57.296)
    Obit.PlotText(plot.me, x, y, dx, dy, just, text, err.me)
    # end  PText

def PRelText (plot, side, disp, coord, fjust, text, err):
    """ Write text on plot relative to port 

    plot   =  Python Plot object
    side   =  Must include one of the characters 'B', 'L', 'T',
              or 'R' signifying the Bottom, Left, Top, or Right
              margin of the viewport. If it includes 'LV' or
              'RV', the string is written perpendicular to the
              frame rather than parallel to it.
    disp   =  The displacement of the character string from the
              specified edge of the viewport, measured outwards
              from the viewport in units of the character
              height. Use a negative value to write inside the
              viewport, a positive value to write outside.
    coord  = The location of the character string along the
             specified edge of the viewport, as a fraction of
             the length of the edge.
    just   = Controls justification of the string parallel to
             the specified edge of the viewport. If
             just = 0.0, the left-hand end of the string will
             be placed at COORD; if JUST = 0.5, the center of
             the string will be placed at COORD; if JUST = 1.0,
             the right-hand end of the string will be placed at
             at COORD. Other values between 0 and 1 give
             intermediate placing, but they are not very useful.
    text   = The text string to be plotted. Trailing spaces are
             ignored when justifying the string, but leading
             spaces are significant.
    err =    ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    Obit.PlotRelText(plot.me, side, disp, coord, fjust, text, err.me)
    # end  PRelText

def PDrawLine (plot, x1, y1, x2, y2, err):
    """ Draw a line.

    plot   = Python Plot object
    x1      = world x-coordinate of the new pen position.
    y1      = world y-coordinate of the new pen position.
    x2      = world x-coordinate of the new pen position.
    y2      = world y-coordinate of the new pen position.
    err     = ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    Obit.PlotDrawLine(plot.me, x1, y1, x2, y2, err.me)
    # end  PDrawLine

def PDrawCurve (plot, x, y, err):
    """ Draw a curve.

    plot   = Python Plot object
    x      =  Array of world x-coordinates of points
    y      =  Array of world y-coordinates of points
    err    =  ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    n = len(x)
    Obit.PlotDrawCurve (plot.me, n, x, y, err.me)
    # end  PDrawCurve 

def PDrawSymbol (plot, x, y, symbol, err):
    """ Draw a Symbol

    plot   = Python Plot object
    x      =  world x-coordinate of the center of the symbol
    y       = world y-coordinate of the center of the symbol
    symbol  = Symbol index to use for plotting
              values in the range [1,12] are usable 
        1 = dot
        2 = plus
        3 = *
        4 = open circle
        5 = x
        6 = open square
        7 = open triangle
        8 = open star
        9 = filled triangle
        10 = filled square
        11 = filled circle
        12 = filled star
    err =    ObitErr error stack
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    Obit.PlotDrawSymbol(plot.me, x, y, symbol, err.me)
    # end  PDrawSymbol

def PGetList (plot):
    """ Return the member InfoList

    returns InfoList
    plot   = Python Obit Plot object
    """
    ################################################################
     # Checks
    if not PIsA(plot):
        raise TypeError,"plot MUST be a Python Obit plot"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.PlotGetList(plot.me)
    return out
    # end PGetList

def PIsA (plot):
    """ Tells if the input is a Python ObitPlot

    returns true or false (1,0)
    plot = Python Obit Plot to test
    """
    ################################################################
      # Checks
    if plot.__class__ != OPlot:
        return 0
    return Obit.OPlotIsA(plot.me)
    # end PIsA


