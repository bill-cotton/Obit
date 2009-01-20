# routines to run Ringer and enter values in ST tables
# Should be run from ObitTalk

import Obit, ObitTask, Table, TableSTar

def run (in1, in2, logFile, err, \
         Cutoff=0.03, Center=[0.0,0.0], Spacing=2.0, NSearch=51, prtLv=0):
    """ Run Ringer to fit the diameters and widths of SiO maser rings

    Run Ringer and if successful, write ST table
    also log fitted values.
    This is fed pairs of images at the same time and geometry and with
    the same geometries, use Ringer to fit rings.
    in1       = First image of pair
    in2       = Second image of pair
    logFile   = Text file to write results into, if None no log
    err       = Python Obit Error/message stack
    Cutoff    = minimum pixel value to use
    Center    = Initial center pixel for search
    Spacing   = Spacing (pixels) of annular sums
    NSearch   = Number of Spacings to search
    prtLv     = Diagnostic print level
    """
    ################################################################
    ring = ObitTask.ObitTask("Ringer")
    ring.DataType= in1.FileType
    ring.inDisk  = in1.Disk
    ring.in2Disk = in2.Disk
    if (ring.DataType=="AIPS"):
        ring.inName   = in1.Aname
        ring.inClass  = in1.Aclass
        ring.inSeq    = in1.Aseq
        ring.in2Name  = in2.Aname
        ring.in2Class = in2.Aclass
        ring.in2Seq   = in2.Aseq
    else:
        ring.inFile   = in1.FileName
        ring.in2File  = in2.FileName
    ring.Cutoff  = Cutoff
    ring.Center  = Center
    ring.Spacing = Spacing
    ring.NSearch = NSearch
    ring.prtLv   = prtLv
    #ring.debug = True  # DEBUG
    ring.g    # Run
    #
    # Save to logging file
    if logFile:
        log = open(logFile, "a")
        line = ring.inName+" "+ring.in2Name+" "+str(ring.Center)+" "+str(ring.Radius)+" "+\
               str(ring.Width)+" "+str(ring.Frac)+" "+in1.Desc.Dict["obsdat"]+"\n"
        print line
        log.write(line)
        log.close()
    #
    # Write ST tables
    st=TableSTar.PCreate(in1, err)
    st.Open(Table.WRITEONLY,err)
    TableSTar.PWriteCirc(st, in1, ring.Center, 2*ring.Radius[0], err)
    st.Close(err)
    # Second image
    st=TableSTar.PCreate(in2, err)
    st.Open(Table.WRITEONLY,err)
    TableSTar.PWriteCirc(st, in2, ring.Center, 2*ring.Radius[1], err)
    st.Close(err)
    # end run


