""" Utility to generate multiple IFs from single IF Data
"""
# $Id: 
#-----------------------------------------------------------------------
#  Copyright (C) 2012
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

def UVAddIF (inUV, outUV, nIF, err):
    """ 
    Create outUV like inUV but divided into nIF IFs
    
    
    * inUV        = input Obit UV object
    * outUV       = output Obit UV object, defined but not instantiated
                    Onlyt really works for AIPS
    * nIF         = number of desired output IFs
                    MUST be the same number of channels per IF
    * err         = Obit error/message stack
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not UV.PIsA(outUV):
        raise TypeError,"outUV MUST be a defined Python Obit UV"
    # Input can have 1 or no IFs defined
    jlocif = inUV.Desc.Dict["jlocif"]
    if jlocif>=0 and inUV.Desc.Dict["inaxes"][jlocif]>1:
        raise RuntimeError,"Input UV has excessive IFs already:" \
              +str(inUV.Desc.Dict["inaxes"][jlocif])
    # Check number of requested IFs
    if nIF<2:
        raise RuntimeError,"Too few output IFs requested: " + str(nIF)
    # Must be an even number of channels per output IF
    jlocf = inUV.Desc.Dict["jlocf"]
    nchan = inUV.Desc.Dict["inaxes"][jlocf]
    if (nchan%nIF) != 0:
        raise RuntimeError,"Unequal numbers of channels per "+str(nIF)+" IFs"

    # Patch UV Descriptor
    DescAddIF (inUV, outUV, nIF, err)
    # Convert FG Table
    UpdateFQ (inUV, outUV, nIF, err)
    # Convert AN Table
    UpdateAN (inUV, outUV, nIF, err)
    # Convert SU Table
    UpdateSU (inUV, outUV, nIF, err)
    # Copy data
    CopyData (inUV, outUV, err)
    # Update
    outUV.UpdateDesc(err)
    OErr.printErrMsg(err,"Error updating output")
     # end UVAddIF 
    
def DescAddIF (inUV, outUV, nIF, err):
    """ 
    Update outUV descriptor like inUV but divided into nIF IFs
    
    
    * inUV        = input Obit UV object
    * outUV       = output Obit UV object, defined but not instantiated
    * nIF         = number of desired output IFs
                    MUST be the same number of channels per IF
    * err         = Obit error/message stack
    """
    ################################################################
    inUV.Clone(outUV, err)
    d = inUV.Desc.Dict
    # Have IF axis?
    if d["jlocif"]>=0:
        jlocif = d["jlocif"]
    else:  # create one
        jlocif = d['naxis']
        d['naxis'] += 1
        d['ctype'][jlocif] = "IF"
        d['crval'][jlocif] = 1.0
        d['crpix'][jlocif] = 1.0
        d['cdelt'][jlocif] = 1.0
        d['crota'][jlocif] = 0.0
        
    jlocf = d["jlocf"]
    nchan = d["inaxes"][jlocf]/nIF
    d["inaxes"][jlocif] = nIF
    d["inaxes"][jlocf]  = nchan
    outUV.Desc.Dict = d
    UV. PGetIODesc(outUV).Dict = d  # And your little dog too
    # Update
    outUV.UpdateDesc(err)
    outUV.Open(UV.WRITEONLY,err)
    outUV.Close(err)
    #outUV.Header(err)
    OErr.printErrMsg(err,"Error converting Descriptor")
    # end DescAddIF
    
def UpdateFQ (inUV, outUV, nIF, err):
    """ 
    Convert FQ table in inUV (1 fqid)
    
    * inUV        = input Obit UV object
    * outUV       = output Obit UV object, defined but not instantiated
    * nIF         = number of desired output IFs
                    MUST be the same number of channels per IF
    * err         = Obit error/message stack
    """
    ################################################################
    iFQTab = inUV.NewTable(Table.READONLY, "AIPS FQ",1,err)
    # zap FQ old
    outUV.ZapTable("AIPS FQ",1,err)
    oFQTab = outUV.NewTable(Table.WRITEONLY, "AIPS FQ",1,err,numIF=nIF)
    # Input info
    d = inUV.Desc.Dict
    jlocf = d["jlocf"]
    reffreq = d["crval"][jlocf]
    delfreq = d["cdelt"][jlocf]
    freqpix = d["crpix"][jlocf]
    nchan   = d["inaxes"][jlocf]/nIF
    # Can use row from input table
    iFQTab.Open(Table.READONLY, err)
    oFQTab.Open(Table.WRITEONLY, err)
    row = iFQTab.ReadRow(1,err)
    iFQTab.Close(err)
    
    # Update row for nif IFs
    freqarr = []
    chw = row['CH WIDTH'][0]
    chwarr = []
    tbw = row['TOTAL BANDWIDTH'][0] / nIF
    tbwarr = []
    sideband = row['SIDEBAND'][0]
    sbarr = []
    rxc = row['RXCODE'][0]
    rxcarr = ""
    for iIF in range(1,nIF+1):
        freqarr.append((iIF-freqpix)*delfreq*nchan )
        chwarr.append(chw)
        tbwarr.append(tbw)
        sbarr.append(sideband)
        rxcarr += rxc
    row['IF FREQ']         = freqarr
    row['CH WIDTH']        = chwarr
    row['TOTAL BANDWIDTH'] = tbwarr
    row['SIDEBAND']        = sbarr
    row['RXCODE']          = [rxcarr]
        
    # Write output
    oFQTab.WriteRow(1,row,err)
    oFQTab.Close(err)
    # Update
    outUV.UpdateDesc(err)
    OErr.printErrMsg(err,"Error converting FQ Table")
    # end UpdateFQ 
    
def UpdateAN (inUV, outUV, nIF, err):
    """ 
    Convert AN table in inUV to outUV with nIF IFs
    
    * inUV        = input Obit UV object
    * outUV       = output Obit UV object, defined but not instantiated
    * nIF         = number of desired output IFs
    * err         = Obit error/message stack
    """
    ################################################################
    iANTab = inUV.NewTable(Table.READONLY, "AIPS AN",1,err)
    # zap AN old
    outUV.ZapTable("AIPS AN",1,err)
    oANTab = outUV.NewTable(Table.WRITEONLY, "AIPS AN",1,err,numIF=nIF)
    
    iANTab.Open(Table.READONLY, err)
    oANTab.Open(Table.WRITEONLY, err)

    nrow = iANTab.Desc.Dict['nrow']  # How many rows?
    for irow in range(1,nrow+1):
        row = iANTab.ReadRow(irow,err)  # Read input row
        
        # Update  row
        pca0 = row['POLCALA'][0]
        pca1 = row['POLCALA'][1]
        pcb0 = row['POLCALB'][0]
        pcb1 = row['POLCALB'][1]
        bm   = row['BEAMFWHM'][0]
        Beama    = []
        PolcalAa = []
        PolcalBa = []
        for i in range (0,nIF):
            Beama.append(bm)
            PolcalAa.append(pca0)
            PolcalAa.append(pca1)
            PolcalBa.append(pcb0)
            PolcalBa.append(pcb1)
        
        row['BEAMFWHM'] = Beama
        row['POLCALA']  = PolcalAa
        row['POLCALB']  = PolcalBa
        
        # Write output
        oANTab.WriteRow(irow,row,err)
        # End loop over rows
    iANTab.Close(err)
    oANTab.Close(err)
    # Update
    outUV.UpdateDesc(err)
    OErr.printErrMsg(err,"Error converting AN Table")
    # end UpdateAN 
    
def UpdateSU (inUV, outUV, nIF, err):
    """ 
    Convert SU table in inUV to outUV with nIF IFs
    
    * inUV        = input Obit UV object
    * outUV       = output Obit UV object, defined but not instantiated
    * nIF         = number of desired output IFs
    * err         = Obit error/message stack
    """
    ################################################################
    iSUTab = inUV.NewTable(Table.READONLY, "AIPS SU",1,err)
    # zap SU old
    outUV.ZapTable("AIPS SU",1,err)
    oSUTab = outUV.NewTable(Table.WRITEONLY, "AIPS SU",1,err,numIF=nIF)

    iSUTab.Open(Table.READONLY, err)
    oSUTab.Open(Table.WRITEONLY, err)
    nrow = iSUTab.Desc.Dict['nrow']  # How many rows?

    for irow in range(1,nrow+1):
        row = iSUTab.ReadRow(irow,err)  # Read input row
        
        # Update row
        fo   = row['FREQOFF'][0]
        bw   = row['BANDWIDTH'][0]
        iflx = row['IFLUX'][0]
        qflx = row['QFLUX'][0]
        uflx = row['UFLUX'][0]
        vflx = row['VFLUX'][0]
        lsr  = row['LSRVEL'][0]
        rest = row['RESTFREQ'][0]
        iflxa = []; qflxa = []; uflxa = []; vflxa = [];
        foa = []; lsra = []; rfa = []
        for i in range (0,nIF):
            iflxa.append(iflx)
            qflxa.append(qflx)
            uflxa.append(uflx)
            vflxa.append(vflx)
            foa.append(fo)
            lsra.append(lsr)
            rfa.append(rest)
         
        row['IFLUX']     = iflxa
        row['QFLUX']     = qflxa
        row['UFLUX']     = uflxa
        row['VFLUX']     = vflxa
        row['FREQOFF']   = foa
        row['LSRVEL']    = lsra
        row['RESTFREQ']  = rfa
        
        # Write output
        oSUTab.WriteRow(irow,row,err)
        # end loop over rows
    iSUTab.Close(err)
    oSUTab.Close(err)
    # Update
    outUV.UpdateDesc(err)
    OErr.printErrMsg(err,"Error converting SU Table")
    # end UpdateSU 
    
def CopyData (inUV, outUV, err):
    """ 
    Copy the raw visibility records from inUV to outUV
    
    * inUV        = input Obit UV object
    * outUV       = output Obit UV object, defined but not instantiated
    * err         = Obit error/message stack
    """
    ################################################################
    # Checks, sizes should be the same
    id = inUV.Desc.Dict
    od = outUV.Desc.Dict
    odd = UV.PGetIODesc(outUV)
    if id['nrparm'] != od['nrparm']:
        raise RuntimeError,"Input and output have different numbers of random parameters: " \
              +str(id['nrparm'])+" != "+str(od['nrparm'])
    if id['ncorr'] != od['ncorr']:
        raise RuntimeError,"Input and output have different numbers of correlations: " \
              +str(id['ncorr'])+" != "+str(od['ncorr'])

    nvis = inUV.Desc.Dict['nvis']                     # Number of records
    lrec = id['nrparm'] + id['ncorr']*id["inaxes"][0]  # Size of record in floats

    # Set data to read one vis per IO
    inUV.List.set("nVisPIO", 1)
    outUV.List.set("nVisPIO", 1)
    
    # Open files
    zz = inUV.Open(UV.READONLY, err)
    zz = outUV.Open(UV.READWRITE, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening UV")
    # Get IO buffers as numpy arrays
    shape = len(inUV.VisBuf) / 4
    ibuffer =  numarray.array(sequence=inUV.VisBuf,
                             type=numarray.Float32, shape=shape)
    obuffer =  numarray.array(sequence=outUV.VisBuf,
                             type=numarray.Float32, shape=shape)
    od['numVisBuff'] = 1; # Gotta tell it to write 1 vis
    outUV.Desc.Dict   = od
    outUV.IODesc.Dict = od

    # Loop over data, copying 1 vis at a time
    for ivis in range(1,nvis+1):
        inUV.Read(err, firstVis=ivis)   # Read input
        if err.isErr:
            OErr.printErrMsg(err, "Error reading data")
        # Copy record
        obuffer[0:lrec] = ibuffer[0:lrec]
        outUV.Desc.Dict = od
        outUV.Write(err, firstVis=ivis)   # Write output
        if err.isErr:
            OErr.printErrMsg(err, "Error writing data")
        # end loop over data
    inUV.Close(err)
    outUV.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing data")
 
