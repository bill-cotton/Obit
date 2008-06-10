# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2006
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

# Python interface to GBT CCB utilities
import Obit, Image, Table, ImageUtil, ImageDesc, OTF, FArray, History, OErr
import OTFUtil, OTFDesc, InfoList
import math

def PNod (inOTF, scan, err):
    """ Process nodding scan

    Differences on source and offsource in a nodding scan and returns
    a dictionary object containing the channel averages,
    error estimates and frequencies.
    Target table must contain position
    
    inOTF     = Input Python OTF, any calibration should be soecified
    scan      = Scan number
    err       = Python Obit Error/message stack
    returns distionary with entries:
    "avg"  = array of channel average differences (0.5 the on-off)
    "rms"  = array of rms about avg
    "freq" = array of frequencies of channels
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    OTFUtil.PDiffNod (inOTF, scan, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error analysing noddong data")
    # Get results
    inInfo = OTF.PGetList(inOTF)
    onOff  = InfoList.PGet (inInfo, "OnOff")[4]
    #print "DEBUG onOff", onOff

    # Header info
    inDesc = OTF.PGetDesc(inOTF)
    inDescDict = OTFDesc.PGetDict(inDesc)
    Freq0 = inDescDict["crval"][inDescDict["jlocf"]] - inDescDict["cdelt"][inDescDict["jlocf"]] * \
            (inDescDict["crpix"][inDescDict["jlocf"]] - 1.0)
    dFreq = inDescDict["cdelt"][inDescDict["jlocf"]]
    
    # Number of channels
    nchan = len(onOff)/4
    # Average, get rms and frequency
    avg = []
    rms = []
    freq = []
    for i in range(0,nchan):
        sum = onOff[4*i+0] + onOff[4*i+1] + onOff[4*i+2] + onOff[4*i+3]
        ave =sum / 4
        avg.append(ave*0.5)
        # RMS
        t = onOff[4*i]
        sum = (t - ave) * (t - ave)
        t =  onOff[4*i+1]
        sum = sum + (t - ave) * (t - ave)
        t =  onOff[4*i+2]
        sum = sum + (t - ave) * (t - ave)
        t =  onOff[4*i+3]
        sum = sum + (t - ave) * (t - ave)
        r = math.sqrt( 0.5 * sum/3)
        rms.append(r)
        # Frequency
        freq.append(Freq0+dFreq*i)

    # Set output distionary
    out = {"avg":avg, "rms":rms, "freq":freq, "onOff":onOff}
    return out
    # end PCCScale

