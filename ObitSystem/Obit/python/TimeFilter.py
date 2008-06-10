""" Time domain filtering class
"""
# $Id: TimeFilter.py,v 1.1 2008/01/29 02:24:33 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C) 2008
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

# Python shadow class to ObitTimeFilter class
import Obit, FFT

class TimeFilterPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.TimeFilterUnref(Obit.TimeFilter_me_get(self.this))
            # In with the new
            Obit.TimeFilter_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != TimeFilter:
            return
        if name == "me" : 
            return Obit.TimeFilter_me_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        if self.__class__ != TimeFilter:
            return
        return "<C TimeFilter instance> " + Obit.TimeFilterGetName(self.me)

class TimeFilter(TimeFilterPtr):
    """ Python Obit interface to display server
    
    This class is for creating and using the interface to a plot
    TimeFilter Members with python interfaces:

    """
    def __init__(self, name, nTime, nSeries):
        self.this = Obit.new_TimeFilter(name, nTime, nSeries)
    def __del__(self):
        if Obit!=None:
            Obit.delete_TimeFilter(self.this)

# Define filter types
LowPass    = 0
HighPass   = 1
NotchPass  = 2
NotchBlock = 3

def newTimeFilter(name, nTime, nSeries):
    """ Create and initialize an ObitTimeFilter

    name     = name desired for object (labeling purposes)
    nTime    = Number of elements in time series
    nSeries  = Number of time serieses
    """
    ################################################################
    out = TimeFilter(name, nTime, nSeries)
    return out 
    # end newTimeFilter

def PTimeFilterResize (filter, nTime):
    """ Resize Filter

    Change number of time samples in filter
    Note: the actual size will be the result of FFT.PSuggestSize(nTime)
    filter   = time filter
    nTime    = minimum number of time elements
    """
    ################################################################
    # Checks
    if not PIsA(filter):
        print "Actually ",filter.__class__
        raise TypeError,"filter MUST be a Python Obit TimeFilter"
    n = FFT.PSuggestSize(nTime);  # How many actual points?
    Obit.TimeFilterResize (filter.me, n)
    # end PXYTimeFilter


def PGridTime (filter, seriesNo, dTime, nTime, times, data):
    """ Convert time-ordered data stream into regular time grid.

    Will resize in if needed.
    Data will be averaging into time bins.

    filter     = time filter
    seriesNo   = Which time/frequency series to apply to (0-rel)
    dTime      = Increment of desired time grid (days)
    nTime      = Number of times in times, data
    times      = list of times (days)
    data       = list of data elements corresponding to times.
    """
    ################################################################
    # Checks
    if not PIsA(filter):
        print "Actually ",filter.__class__
        raise TypeError,"filter MUST be a Python Obit TimeFilter"
    Obit.TimeFilterGridTime (filter.me, seriesNo, dTime, nTime, times, data)
    # end PGridTime

def PUngridTime (filter, seriesNo, nTime, times):
    """ Copy time series to external form

    Time series data are copied by nearest time stamp
    Returns list of data elements corresponding to times.
    filter     = time filter
    seriesNo   = Which time/frequency series to apply to (0-rel)
    nTime      = Number of times in times, data
    times      = list of times (days)
    """
    ################################################################
    # Checks
    if not PIsA(filter):
        print "Actually ",filter.__class__
        raise TypeError,"filter MUST be a Python Obit TimeFilter"
    return Obit.TimeFilterUngridTime (filter.me, seriesNo, nTime, times)
    # end PUngridTime

def P2Freq (filter):
    """ Fourier transform all time series to frequency series

    Any blanked values in the time series are interpolated and then the time series 
    is FFTed to the frequency domain.  A linear interpolation between the 
    last valid point and the first valid point is made to reduce the wraparound edge
    effects.
    filter     = time filter to FFT to frequency
    """
    ################################################################
    # Checks
    if not PIsA(filter):
        print "Actually ",filter.__class__
        raise TypeError,"filter MUST be a Python Obit TimeFilter"
    Obit.TimeFilter2Freq (filter.me)
    # end P2Freq


def P2Time (filter):
    """ Fourier transform filtered frequency series to time series

    filter     = time filter to FFT to frequency
    """
    ################################################################
    # Checks
    if not PIsA(filter):
        print "Actually ",filter.__class__
        raise TypeError,"filter MUST be a Python Obit TimeFilter"
    Obit.TimeFilter2Time (filter.me)
    # end P2Time


def PFilter (filter, seriesNo, type, freq, err):
    """ Apply specified filter to specified time series

    filter   = filter
    seriesNo = Which time/frequency series to apply to (0-rel), <0 => all
    type     = Filter type to apply
               LowPass: Zeroes frequencies above freq[0] (Hz)
               HighPass: Zeroes frequencies below freq[0] (Hz)
               NotchPass: Zeroes frequencies not in frequency range freq[0]->freq[1]
               NotchBlock:Zeroes frequencies in frequency range freq[0]->freq[1]
    freq     = array of frequencies (Hz) for filter, meaning depends on type.
    err      =    ObitErr error stack
    """
    ################################################################
    # Checks
    if not PIsA(filter):
        print "Actually ",filter.__class__
        raise TypeError,"filter MUST be a Python Obit TimeFilter"
    if freq.__class__==list:
        lfreq = freq
    else:
        lfreq=[freq]
    Obit.TimeFilterDoFilter (filter.me, seriesNo, type, lfreq, err.me)
    # end PFilter


def PPlotPower (filter, seriesNo, label, err):
    """ Plot the power spectrum of a selected frequency series

    Output of plot will be prompted for
    filter   = Time filter to plot
    seriesNo = Which time/frequency series to plot (0-rel)
    label    = Label for plot
    err      = ObitErr error stack
    """
    ################################################################
    # Checks
    if not PIsA(filter):
        print "Actually ",filter.__class__
        raise TypeError,"filter MUST be a Python Obit TimeFilter"
    Obit.TimeFilterPlotPower (filter.me, seriesNo, label, err.me)
    # end PPlotPower

def PPlotTime (filter, seriesNo, label, err):
    """ Plot the selected time series

    Output of plot will be prompted for
    filter   = Time filter to plot
    seriesNo = Which time/frequency series to plot (0-rel)
    label    = Label for plot
    err      = ObitErr error stack
    """
    ################################################################
    # Checks
    if not PIsA(filter):
        print "Actually ",filter.__class__
        raise TypeError,"filter MUST be a Python Obit TimeFilter"
    Obit.TimeFilterPlotTime (filter.me, seriesNo, label, err.me)
    # end PPlotTime

def PGetTime (filter, seriesNo):
    """ Returns time data for a selected series

    Returns dict with members:
        "dTime" = time increment in days
        "time"  = array of times in days
        "data"  = array of data corresponding to "time"
    filter   = Python TimeFilter object
    seriesNo = Which time/frequency series to return (0-rel)
    """
    ################################################################
     # Checks
    if not PIsA(filter):
        raise TypeError,"filter MUST be a Python Obit filter"
    #
    Obit.TimeFilterGetTime(filter.me, seriesNo)
    # end  PGetTime

def PGetFreq (filter, seriesNo):
    """ Returns frequency spectra for a selected series

    Returns dict with members:
        "dFreq" = frequency increment in Hz
        "freq"  = array of frequencies in Hz
        "data"  = array of complex data corresponding to "freq"
    filter   = Python TimeFilter object
    seriesNo = Which time/frequency series to return (0-rel)
    """
    ################################################################
     # Checks
    if not PIsA(filter):
        raise TypeError,"filter MUST be a Python Obit filter"
    #
    Obit.TimeFilterGetFreq(filter.me, seriesNo)
    # end  PGetTime

def PSetTime (filter, seriesNo, inDict):
    """ Sets time data for a selected series

    filter   = Python TimeFilter object
    seriesNo = Which time/frequency series to return (0-rel)
    inDict   = dict with members: (as returned by PGetTime)
        "dTime" = time increment in days
        "time"  = array of times in days
        "data"  = array of data corresponding to "time"
   """
    ################################################################
     # Checks
    if not PIsA(filter):
        raise TypeError,"filter MUST be a Python Obit filter"
    #
    Obit.TimeFilterGetTime(filter.me, seriesNo, inDict)
    # end  PSetTime

def PSetFreq (filter, seriesNo, inDict):
    """ Returns frequency spectra for a selected series

    filter   = Python TimeFilter object
    seriesNo = Which time/frequency series to return (0-rel)
    inDict   = dict with members: (as returned by PGetTime)
    Returns dict with members:
        "dFreq" = frequency increment in Hz
        "freq"  = array of frequencies in Hz
        "data"  = array of complex data corresponding to "freq"
  """
    ################################################################
     # Checks
    if not PIsA(filter):
        raise TypeError,"filter MUST be a Python Obit filter"
    #
    Obit.TimeFilterGetFreq(filter.me, seriesNo)
    # end  PSetTime

def PIsA (filter):
    """ Tells if the input is a Python ObitTimeFilter

    returns true or false (1,0)
    filter = Python Obit TimeFilter to test
    """
    ################################################################
      # Checks
    if filter.__class__ != TimeFilter:
        return 0
    return Obit.TimeFilterIsA(filter.me)
    # end PIsA


