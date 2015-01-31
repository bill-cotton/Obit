""" Utility to convert KAT-7 HDF format data to AIPS (or FITS)

This module requires katdal and katpoint and their dependencies
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2014
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
try:
    import katdal as katfile
    import katpoint
except Exception, exception:
    print exception
    print "KAT software not available"
    raise  RuntimeError, "KAT software unavailable"
else:
    pass
import time,math,os
import UV, UVVis, OErr, UVDesc, Table, History
from OTObit import day2dhms
import numpy
from numpy import numarray

def KAT2AIPS (katdata, outUV, disk, fitsdisk, err, \
              calInt=1.0, **kwargs):
    """Convert KAT-7 HDF 5 data set to an Obit UV.

    This module requires katdal and katpoint and their dependencies
    contact Ludwig Schwardt <schwardt@ska.ac.za> for details.

    Parameters
    ----------
    katdata : string
        input katfile object
    outUV : ??
        Obit UV object, should be a KAT template for the
        appropriate number of IFs and poln.
    disk  : int
        AIPS Disk number
    fitsdisk: int
        FITS Disk number
    err : ??
        Obit error/message stack
    calInt : 
        Calibration interval in min.
    targets : list, optinal
        List of targetnames to extract from the file
    """
    ################################################################
    OErr.PLog(err, OErr.Info, "Converting h5 data to AIPS UV format.")
    OErr.printErr(err)
    print "Converting h5 data to AIPS UV format.\n"

    # Extract metadata
    meta = GetKATMeta(katdata, err)

    # Extract AIPS parameters of the uv data to the metadata
    meta["Aproject"] = outUV.Aname
    meta["Aclass"] = outUV.Aclass
    meta["Aseq"] = outUV.Aseq
    meta["Adisk"] = disk
    meta["calInt"] = calInt
    meta["fitsdisk"] = fitsdisk
    # Update descriptor
    UpdateDescriptor (outUV, meta, err)
    # Write AN table
    WriteANTable (outUV, meta, err)
    # Write FQ table
    WriteFQTable (outUV, meta, err)
    # Write SU table
    WriteSUTable (outUV, meta, err)

    # Convert data
    ConvertKATData(outUV, katdata, meta, err)

    # Index data
    OErr.PLog(err, OErr.Info, "Indexing data")
    OErr.printErr(err)
    UV.PUtilIndex (outUV, err)

    # Open/close UV to update header
    outUV.Open(UV.READONLY,err)
    outUV.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, message="Update UV header failed")

    # initial CL table
    OErr.PLog(err, OErr.Info, "Create Initial CL table")
    OErr.printErr(err)
    print "Create Initial CL table\n"
    UV.PTableCLGetDummy(outUV, outUV, 1, err, solInt=calInt)
    outUV.Open(UV.READONLY,err)
    outUV.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, message="Update UV header failed")

    # History
    outHistory = History.History("outhistory", outUV.List, err)
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp("Convert KAT7 HDF 5 data to Obit", err)
    outHistory.WriteRec(-1,"datafile = "+katdata.name, err)
    outHistory.WriteRec(-1,"calInt   = "+str(calInt), err)
    outHistory.Close(err)
    outUV.Open(UV.READONLY,err)
    outUV.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, message="Update UV header failed")
    # Return the metadata for the pipeline
    return meta
   # end KAT2AIPS

def GetKATMeta(katdata, err):
    """
    Get KAT metadata and return as a dictionary.

    Parameters
    ----------
     * katdata  = input KAT dataset
     * err      = Python Obit Error/message stack to init

     Returns : dictionary
     "spw"     Spectral window array as tuple (nchan, freq0, chinc)
               nchan=no. channels, freq0 = freq of channel 0,
               chinc=signed channel increment, one tuple per SPW
     "targets" Array of target tuples:
               (index, name, ra2000, dec2000, raapp, decapp)
     "bpcal"   List of source indices of Bandpass calibrators
     "gaincal" List of source indices of Gain calibrators
     "source" List of source indices of imaging targets
     "targLookup" dict indexed by source number with source index
     "tinteg"   Integration time in seconds
     "obsdate"  First day as YYYY-MM-DD
     "observer" name of observer
     "ants"     Array of antenna tuples (index, name, X, Y, Z, diameter)
     "nstokes"  Number of stokes parameters
     "products" Tuple per data product (ant1, ant2, offset)
                where offset is the index on the Stokes axis (XX=0...)
    """
    ################################################################
    # Checks
    out = {}
    # Spectral windows
    sw = []
    for s in katdata.spectral_windows:
        sw.append((s.num_chans, s.channel_freqs[0], s.channel_freqs[1]-s.channel_freqs[0]))
    out["spw"] = [sw[0]]
    # targets
    tl = []
    tb = []
    tg = []
    tt = []
    td = {}
    i = 0

    for ti in katdata.target_indices:
        t = katdata.catalogue.targets[ti]
        #Aips doesn't like spaces in names!!
        name = (t.name+"                ")[0:16]
        ras, decs = t.radec()
        dec = UVDesc.PDMS2Dec(str(decs).replace(':',' '))
        ra  = UVDesc.PHMS2RA(str(ras).replace(':',' '))
        # Apparent posn
        ras, decs = t.apparent_radec()
        deca = UVDesc.PDMS2Dec(str(decs).replace(':',' '))
        raa  = UVDesc.PHMS2RA(str(ras).replace(':',' '))
        i += 1
        tl.append((i, name, ra, dec, raa, deca))
        if 'bpcal' in t.tags:
            tb.append(t)
        elif 'gaincal' in t.tags:
            tg.append(t)
        else:                   # Assume all other targets are for imaging
            tt.append(t)
        td[name.rstrip()] = i
    out["targets"] = tl
    out["targLookup"] = td
    out["bpcal"] = tb
    out["gaincal"] = tg
    out["source"] = tt
    # Antennas
    al = []
    alook = {}
    i = 0
    for a in katdata.ants:
        name  = a.name
        x,y,z = a.position_ecef
        diam  = a.diameter
        i += 1
        al.append((i, name, x, y, z, diam))
        alook[name] = i
    out["ants"] = al
    out["antLookup"] = alook
    # Data products
    dl = []
    nstokes = 1
    for d in katdata.corr_products:
        a1 = out["antLookup"][d[0][:4]]
        a2 = out["antLookup"][d[1][:4]]
        if d[0][4:]=='h' and d[1][4:]=='h':
            dp = 0
        elif d[0][4:]=='v' and d[1][4:]=='v':
            dp = 1
        elif d[0][4:]=='h' and d[1][4:]=='v':
            dp = 2
        else:
            dp = 3
        dl.append((a1, a2, dp))
        nstokes = max (nstokes,dp+1)
    out["products"] = dl
    out["nstokes"]  = nstokes
    # integration time
    out["tinteg"] = katdata.dump_period
    # observing date
    start=time.gmtime(katdata.timestamps[0])
    out["obsdate"] = time.strftime('%Y-%m-%d', start)
    # Observer's name
    out["observer"] = katdata.observer
    # Number of channels
    numchan = len(katdata.channels)
    out["numchan"] = numchan
    # Correlator mode (assuming 1 spectral window KAT-7)
    out["corrmode"] = katdata.spectral_windows[0].product
    # Central frequency (in Hz)
    out["centerfreq"] = katdata.channel_freqs[numchan //2]
    # Expose all KAT-METADATA to calling script
    out["katdata"] = katdata
    return out
    # end GetKATMeta

def UpdateDescriptor (outUV, meta, err):
    """
    Update information in data descriptor

    NB: Cannot change geometry of visibilities
    * outUV    = Obit UV object
    * meta     = dict with data meta data
    * err      = Python Obit Error/message stack to init
    """
    ################################################################
    chinc   =  meta["spw"][0][2]   # Frequency increment
    reffreq =  meta["spw"][0][1]   # reference frequency
    nchan   =  meta["spw"][0][0]   # number of channels
    nif     = len(meta["spw"])     # Number of IFs
    nstok   = meta["nstokes"]      # Number of Stokes products
    desc = outUV.Desc.Dict
    outUV.Desc.Dict = desc
    desc['obsdat']   = meta["obsdate"]
    desc['observer'] = meta["observer"]
    desc['JDObs']    = UVDesc.PDate2JD(meta["obsdate"])
    desc['naxis']    = 6
    desc['inaxes']   = [3,nstok,nchan,nif,1,1,0]
    desc['cdelt']    = [1.0,-1.0,chinc, 1.0, 0.0, 0.0, 0.0]
    desc['crval']    = [1.0, -5.0,reffreq, 1.0, 0.0, 0.0, 0.0]
    desc['crota']    = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    outUV.Desc.Dict = desc
    outUV.UpdateDesc(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error updating UV descriptor")
    # end UpdateDescriptor

def WriteANTable(outUV, meta, err):
    """
    Write data in meta to AN table

     * outUV    = Obit UV object
     * meta     = dict with data meta data
     * err      = Python Obit Error/message stack to init
    """
    ################################################################
    antab = outUV.NewTable(Table.READWRITE, "AIPS AN",1,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error with AN table")
    antab.Open(Table.READWRITE, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening AN table")
    # Update header
    antab.keys['RefDate'] = meta["obsdate"]
    antab.keys['Freq']    = meta["spw"][0][1]
    JD                    = UVDesc.PDate2JD(meta["obsdate"])
    antab.keys['GSTiat0'] = UVDesc.GST0(JD)*15.0
    antab.keys['DEGPDY']  = UVDesc.ERate(JD)*360.0
    Table.PDirty(antab)  # Force update

    row = antab.ReadRow(1,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error reading AN table")
    OErr.printErr(err)
    irow = 0
    for ant in meta["ants"]:
        irow += 1
        row['NOSTA']    = [ant[0]]
        row['ANNAME']   = [ant[1]+"    "]
        row['STABXYZ']  = [ant[2],ant[3],ant[4]]
        row['DIAMETER'] = [ant[5]]
        row['POLAA']    = [90.0]
        antab.WriteRow(irow, row,err)
        if err.isErr:
            OErr.printErrMsg(err, "Error writing AN table")
    antab.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing AN table")
    # end WriteANTable

def WriteFQTable(outUV, meta, err):
    """
    Write data in meta to FQ table
    An old FQ table is deleted

     * outUV    = Obit UV object
     * meta     = dict with data meta data
     * err      = Python Obit Error/message stack to init
    """
    ################################################################
    # If an old table exists, delete it
    if outUV.GetHighVer("AIPS FQ")>0:
        zz = outUV.ZapTable("AIPS FQ", 1, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error zapping old FQ table")
    reffreq =  meta["spw"][0][1]   # reference frequency
    noif    = 1     # Number of IFs (1 always for KAT7)
    fqtab = outUV.NewTable(Table.READWRITE, "AIPS FQ",1,err,numIF=noif)
    if err.isErr:
        OErr.printErrMsg(err, "Error with FQ table")
    fqtab.Open(Table.READWRITE, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening FQ table")
    # Update header
    fqtab.keys['NO_IF'] = 1  # Structural so no effect
    Table.PDirty(fqtab)  # Force update
    # Create row
    row = {'FRQSEL': [1], 'CH WIDTH': [0.0], 'TOTAL BANDWIDTH': [0.0], \
           'RXCODE': ['L'], 'SIDEBAND': [-1], 'NumFields': 7, 'Table name': 'AIPS FQ', \
           '_status': [0], 'IF FREQ': [0.0]}
    if err.isErr:
        OErr.printErrMsg(err, "Error reading FQ table")
    OErr.printErr(err)
    irow = 0
    for sw in meta["spw"]:
        irow += 1
        row['FRQSEL']    = [irow]
        row['IF FREQ']   = [sw[1] - reffreq]
        row['CH WIDTH']  = [sw[2]]
        row['TOTAL BANDWIDTH']  = [abs(sw[2])*sw[0]]
        row['RXCODE']  = ['L']
        if sw[2]>0.0:
            row['SIDEBAND']  = [1]
        else:
            row['SIDEBAND']  = [-1]
        fqtab.WriteRow(irow, row,err)
        if err.isErr:
            OErr.printErrMsg(err, "Error writing FQ table")
    fqtab.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing FQ table")
    # end WriteFQTable

def WriteFGTable(outUV, katdata, meta, err):
    """
    Get the flags from the h5 file and convert to FG table format.
    UNUSED- This implimentation is too slow!
    outUV    = Obit UV object
    meta     =  dict with data meta data
    err      = Python Obit Error/message stack to init
    """
    ###############################################################

    # work out Start time in unix sec
    tm = katdata.timestamps[1:2]
    tx = time.gmtime(tm[0])
    time0   = tm[0] - tx[3]*3600.0 - tx[4]*60.0 - tx[5]
    int_time = katdata.dump_period

    #Loop through scans in h5 file
    row = 0
    for scan, state, target in katdata.scans():
        name=target.name.replace(' ','_')
        if state != 'track':
            continue
        tm = katdata.timestamps[:]
        nint=len(tm)
        el=target.azel(tm[int(nint/2)])[1]*180./math.pi
        if el<15.0:
            continue
        row+=1
        flags = katdata.flags()[:]
        numflag = 0
        for t, chan_corr in enumerate(flags):
            for c, chan in enumerate(chan_corr):
                cpflagged=[]
                for p, cp in enumerate(chan):
                #for cpaverage in meta['pairLookup']:
                    flag=cp #numpy.any(chan[meta['pairLookup'][cpaverage]])
                    product=meta['products'][p]
                    if product[0] == product[1]:
                        continue
                    if flag and (not product[0:2] in cpflagged):
                        cpflagged.append(product[0:2])
                        numflag+=1
                        starttime=float((tm[t]-time0 - (int_time/2))/86400.0)
                        endtime=float((tm[t]-time0 + (int_time/2))/86400.0)
                        UV.PFlag(outUV,err,timeRange=[starttime,endtime], flagVer=1, \
                                     Ants=[product[0],product[1]], \
                                     Chans=[c+1,c+1], IFs=[1,1], Stokes='1111', Reason='Online flag')
        numvis=t*c*(p/meta["nstokes"])
        msg = "Scan %4d %16s   Online flags: %7s of %8s vis"%(row,name,numflag,numvis)
        OErr.PLog(err, OErr.Info, msg);
        OErr.printErr(err)

def WriteSUTable(outUV, meta, err):
    """
    Write data in meta to SU table

     * outUV    = Obit UV object
     * meta     = dict with data meta data
     * err      = Python Obit Error/message stack to init
    """
    ################################################################
    sutab = outUV.NewTable(Table.READWRITE, "AIPS SU",1,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error with SU table")
    sutab.Open(Table.READWRITE, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening SU table")
    # Update header
    sutab.keys['RefDate'] = meta["obsdate"]
    sutab.keys['Freq']    = meta["spw"][0][1]
    Table.PDirty(sutab)  # Force update
    row = sutab.ReadRow(1,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error reading SU table")
    OErr.printErr(err)
    irow = 0
    for tar in meta["targets"]:
        irow += 1
        row['ID. NO.']   = [tar[0]]
        row['SOURCE']    = [tar[1]]
        row['RAEPO']     = [tar[2]]
        row['DECEPO']    = [tar[3]]
        row['RAOBS']     = [tar[2]]
        row['DECOBS']    = [tar[3]]
        row['EPOCH']     = [2000.0]
        row['RAAPP']     = [tar[4]]
        row['DECAPP']    = [tar[5]]
        row['BANDWIDTH'] = [meta["spw"][0][2]]
        sutab.WriteRow(irow, row,err)
        if err.isErr:
            OErr.printErrMsg(err, "Error writing SU table")
    sutab.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing SU table")
    # end WriteSUtable

def StopFringes(visData,freqData,wData,polProd):
    """
    Stop the fringes for a KAT-7 antenna/polarisation pair.

    * visData  = input array of visibilities
    * freqData = the frequencies corresponding to each channel in visData
    * wData    = the w term for each channel in VisData
    * polProd  = the polarisation product which is being stopped

    Outputs an array of the same size as visData with fringes stopped
    """

    # KAT-7 Antenna Delays (From h5toms.py)
    # NOTE: This should be checked before running (only for w stopping) to see how up to date the cable delays are !!!
    delays = {}
    delays['h'] = {'ant1': 2.32205060e-05, 'ant2': 2.32842541e-05, 'ant3': 2.34093761e-05, 'ant4': 2.35162232e-05, 'ant5': 2.36786287e-05, 'ant6': 2.37855760e-05, 'ant7': 2.40479534e-05}
    delays['v'] = {'ant1': 2.32319854e-05, 'ant2': 2.32902574e-05, 'ant3': 2.34050180e-05, 'ant4': 2.35194585e-05, 'ant5': 2.36741915e-05, 'ant6': 2.37882216e-05, 'ant7': 2.40424086e-05}
    #updated by mattieu 21 Oct 2011 (for all antennas - previously antennas 2 and 4 not updated)
    cable_delay = delays[polProd[1][-1]][polProd[1][:-1]] - delays[polProd[0][-1]][polProd[0][:-1]]
    # Result of delay model in turns of phase. This is now frequency dependent so has shape (tstamps, channels)
    turns = numpy.outer((wData / katpoint.lightspeed) - cable_delay, freqData)
    outVisData = visData*numpy.exp(-2j * numpy.pi * turns)

    return outVisData


def ConvertKATData(outUV, katdata, meta, err):
    """
    Read KAT HDF data and write Obit UV

     * outUV    = Obit UV object
     * katdata  = input KAT dataset
     * meta     = dict with data meta data
     * err      = Python Obit Error/message stack to init
    """
    ################################################################
    reffreq =  meta["spw"][0][1]    # reference frequency
    lamb    = 2.997924562e8/reffreq # wavelength of reference freq
    nchan   =  meta["spw"][0][0]    # number of channels
    nif     = len(meta["spw"])      # Number of IFs
    nstok   = meta["nstokes"]       # Number of Stokes products
    p       = meta["products"]      # baseline stokes indices
    nprod   = len(p)                # number of correlations/baselines
    # work out Start time in unix sec
    tm = katdata.timestamps[1:2]
    tx = time.gmtime(tm[0])
    time0   = tm[0] - tx[3]*3600.0 - tx[4]*60.0 - tx[5]

    # Set data to read one vis per IO
    outUV.List.set("nVisPIO", 1)

    # Open data
    zz = outUV.Open(UV.READWRITE, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening output UV")
    # visibility record offsets
    d = outUV.Desc.Dict
    ilocu   = d['ilocu']
    ilocv   = d['ilocv']
    ilocw   = d['ilocw']
    iloct   = d['iloct']
    ilocb   = d['ilocb']
    ilocsu  = d['ilocsu']
    nrparm  = d['nrparm']
    jlocc   = d['jlocc']
    jlocs   = d['jlocs']
    jlocf   = d['jlocf']
    jlocif  = d['jlocif']
    naxes   = d['inaxes']
    count = 0.0
    visno = 0
    # Get IO buffers as numpy arrays
    shape = len(outUV.VisBuf) / 4
    buffer =  numarray.array(sequence=outUV.VisBuf,
                             type=numarray.Float32, shape=shape)

    # Template vis
    vis = outUV.ReadVis(err, firstVis=1)
    first = True
    firstVis = 1
    numflags=0
    numvis=0
    # Do we need to stop Fringes
    try:
        autodelay=[int(ad) for ad in katdata.sensor['DBE/auto-delay']]
        autodelay=all(autodelay)
    except:
        autodelay=False
    if not autodelay:
        msg = "W term in UVW coordinates will be used to stop the fringes."
        OErr.PLog(err, OErr.Info, msg)
        OErr.printErr(err)
        print msg
    for scan, state, target in katdata.scans():
        # Fetch data
        tm = katdata.timestamps[:]
        nint = len(tm)
        vs = katdata.vis[:]
        wt = katdata.weights()[:]
        fg = katdata.flags()[:]
        #Get target suid
        # Only on targets in the input list
        try:
            suid = meta["targLookup"][target.name[0:16]]
        except:
            continue
        # Negate the weights that are online flagged (ie. apply the online flags here)
        wt = numpy.where(fg,-wt,wt)
        numflags += numpy.sum(fg)
        numvis += fg.size
        uu = katdata.u
        vv = katdata.v
        ww = katdata.w
        # Number of integrations
        msg = "Scan:%4d Int: %4d %16s Start %s"%(scan,nint,target.name,day2dhms((tm[0]-time0)/86400.0)[0:12])
        OErr.PLog(err, OErr.Info, msg);
        OErr.printErr(err)
        print msg
        # Loop over integrations
        for iint in range(0,nint):
            # loop over data products/baselines
            for iprod in range(0,nprod):
                thisvis=vs[iint:iint+1,:,iprod:iprod+1]
                thisw=ww[iint:iint+1,iprod]
                # Fringe stop the data if necessary
                if not autodelay:
                    thisvis=StopFringes(thisvis[:,:,0],katdata.channel_freqs,thisw,katdata.corr_products[iprod])
                # Copy slices
                indx = nrparm+(p[iprod][2])*3
                buffer[indx:indx+(nchan+1)*nstok*3:nstok*3] = thisvis.real.flatten()
                indx += 1
                buffer[indx:indx+(nchan+1)*nstok*3:nstok*3] = thisvis.imag.flatten()
                indx += 1
                buffer[indx:indx+(nchan+1)*nstok*3:nstok*3] = wt[iint:iint+1,:,iprod:iprod+1].flatten()
                # Write if Stokes index >= next or the last
                if (iprod==nprod-1) or (p[iprod][2]>=p[iprod+1][2]):
                    # Random parameters
                    buffer[ilocu]  = uu[iint][iprod]/lamb
                    buffer[ilocv]  = vv[iint][iprod]/lamb
                    buffer[ilocw]  = ww[iint][iprod]/lamb
                    buffer[iloct]  = (tm[iint]-time0)/86400.0 # Time in days
                    buffer[ilocb]  = p[iprod][0]*256.0 + p[iprod][1]
                    buffer[ilocsu] = suid
                    outUV.Write(err, firstVis=visno)
                    visno += 1
                    buffer[3]= -3.14159
                    #print visno,buffer[0:5]
                    firstVis = None  # Only once
                    # initialize visibility
                    first = True
        # end loop over integrations
        if err.isErr:
            OErr.printErrMsg(err, "Error writing data")
    # end loop over scan
    if numvis>0:
        msg= "Applied %s online flags to %s visibilities (%.3f%%)"%(numflags,numvis,float(numflags)/float(numvis))
        OErr.PLog(err, OErr.Info, msg)
        OErr.printErr(err)
    outUV.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing data")
    # end ConvertKATData

