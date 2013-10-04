""" Utility to add a Source table to a UV data

"""
# $Id: 
#-----------------------------------------------------------------------
#  Copyright (C) 2013
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
import UV, Table, OErr
import numpy      
from numpy import numarray

def UVAddSource (inUV, outUV, err):
    """ 
    Add a Source table to a UV data
    
    Source information is derived from the UV descriptor and a Source
    random parameter added if needed.
    * inUV        = input Obit UV object
    * outUV       = output Obit UV object, defined but not instantiated
                    Onlyt really works for AIPS
    * err         = Obit error/message stack
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not UV.PIsA(outUV):
        raise TypeError,"outUV MUST be a defined Python Obit UV"
    # Patch UV Descriptor
    DescAddSource (inUV, outUV, err)
    OErr.printErrMsg(err,"Error updating Descriptor")
    CreateSU (inUV, outUV, err)
    OErr.printErrMsg(err,"Error creating SU Table")
    # Copy data
    CopyData (inUV, outUV, err)
    # Update
    # outUV.UpdateDesc(err)
    outUV.Header(err)   # show results
    OErr.printErrMsg(err,"Error copying")

    # end UVAddSource 
    
def DescAddSource (inUV, outUV, err):
    """ 
    Update outUV descriptor for combined file
       
    * inUV        = input Obit UV object
    * outUV       = output Obit UV object, defined but not instantiated
    * err         = Obit error/message stack
    """
    ################################################################
    inUV.Clone(outUV, err)   # Clone input
    # Add Source id random parameter
    d = outUV.Desc.Dict
    if d["ilocsu"]<0:
        d["ptype"][d["nrparm"]] = "SOURCE"
        d["nrparm"] += 1
        outUV.Desc.Dict = d
        UV.PGetIODesc(outUV).Dict = d  # And your little dog too
    # Update
    outUV.UpdateDesc(err)
    outUV.Open(UV.WRITEONLY,err)
    outUV.Close(err)
    #outUV.Header(err)
    OErr.printErrMsg(err,"Error converting Descriptor")
    # end DescAddSource
    
def CreateSU(inUV, outUV, err):
    """
    Write data in meta to SU table
    
    * inUV       = list of input Obit UV objects (freq order)
    * outUV       = output Obit UV object, defined but not instantiated
    * err         = Obit error/message stack
    """
    ###############################################################
    # Create SU table
    idd = inUV.Desc.Dict   # Info from input header
    nIF = idd["inaxes"][idd["jlocif"]]   # how many IFs
    sutab = outUV.NewTable(Table.READWRITE, "AIPS SU",1, err ,numIF=nIF)
    if err.isErr:
        OErr.printErrMsg(err, "Error with SU table")
    sutab.Open(Table.READWRITE, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening SU table")
    # Update header
    sutab.keys['RefDate'] = idd["obsdat"]
    sutab.keys['Freq']    = idd["crval"][idd["jlocf"]]
    Table.PDirty(sutab)  # Force update
    row = sutab.ReadRow(1,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error reading SU table")
    irow = 1
    IFZero = [] #  Array of NIF zeroes
    for i in range(0,nIF):
        IFZero.append(0.0)

    row['ID. NO.']   = [1]
    row['SOURCE']    = [idd["object"]]
    row['RAEPO']     = [idd["crval"][idd["jlocr"]]]
    row['DECEPO']    = [idd["crval"][idd["jlocd"]]]
    row['RAOBS']     = [idd["obsra"]]
    row['DECOBS']    = [idd["obsdec"]]
    row['EPOCH']     = [idd["equinox"]]
    row['RAAPP']     = [0.0]
    row['DECAPP']    = [0.0]
    row['BANDWIDTH'] = [idd["cdelt"][idd["jlocf"]]]
    row['Table name']= 'AIPS SU'
    row['IFLUX']    = IFZero
    row['QFLUX']    = IFZero
    row['UFLUX']    = IFZero
    row['VFLUX']    = IFZero
    row['FREQOFF']  = IFZero
    row['RESTFREQ'] = IFZero
    row['LSRVEL']   = IFZero
    row['CALCODE']  = ['    ']
    row['NumFields'] = 22
    row['QUAL']     = [0]
    row['PMRA']     = [0.0]
    row['PMDEC']    = [0.0]
    row['_status']  = [1]
    sutab.WriteRow(irow, row,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error writing SU table")
    sutab.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing SU table")
    # end CreateSU

def CopyData (inUV, outUV, err):
    """ 
    Copy the raw visibility records from inUV to outUV
    
    * inUV        = input Obit UV object
    * outUV       = output Obit UV object, defined but not instantiated
    * err         = Obit error/message stack
    """
    ################################################################
    od = outUV.Desc.Dict
    odd = UV.PGetIODesc(outUV)
    idd = inUV.Desc.Dict
    nvis   = idd['nvis']                                     # Number of records
    ilrec  = idd['nrparm']+idd['ncorr']*idd["inaxes"][0]     # Size of input record in floats
    olrec  = od['nrparm']+od['ncorr']*od["inaxes"][0]        # Size of output record in floats
    inrparm = idd['nrparm']                                  # number of input random parameters
    onrparm = od['nrparm']                                   # number of output random parameters
    lcorr  = idd['ncorr']*idd["inaxes"][0]                   # length of correlations
    ilocsu = idd['ilocsu']                                   # Input source ID offset
    olocsu = od['ilocsu']                                    # output source ID offset

    # Set data to read one vis per IO
    inUV.List.set("nVisPIO", 1)
    outUV.List.set("nVisPIO", 1)
    
    # Open files
    zz = outUV.Open(UV.READWRITE, err)
    zz = inUV.Open(UV.READONLY, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening UVs")
    # Get IO buffers as numpy arrays
    shape = len(inUV.VisBuf) / 4
    ibuffer = numarray.array(sequence=inUV.VisBuf,
                                      type=numarray.Float32, shape=shape)
    shape = len(outUV.VisBuf) / 4
    obuffer =  numarray.array(sequence=outUV.VisBuf,
                             type=numarray.Float32, shape=shape)
    # Loop over data, copying 1 vis at a time
    for ivis in range(1,nvis+1):
        inUV.Read(err, firstVis=ivis)   # Read input
        if err.isErr:
            OErr.printErrMsg(err, "Error reading data")
        # Random parameters from first
        obuffer[0:inrparm] = ibuffer[0:inrparm]
        # Source id if making it up
        if ilocsu<0:
            obuffer[olocsu] = 1.0
        # Copy record parts from input
        obuffer[onrparm:olrec] = ibuffer[inrparm:ilrec]
        od = outUV.Desc.Dict
        od['numVisBuff'] = 1;
        outUV.Desc.Dict   = od   # Gotta tell it to write 1 vis
        outUV.Write(err, firstVis=ivis)   # Write output
        if err.isErr:
            OErr.printErrMsg(err, "Error writing data")
        # end loop over data
    #   close up
    outUV.Close(err)
    inUV.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing data")
 
