"""
Fix VLBA data loaded from IDI files produced with the old correlator. IDI files
produced by the DiffX correlator do not require this fix. There are three
potential problems with VLBA IDI data, all of which are fixed here.

1. If the IDI files are not loaded in the correct order, the data must first be
sorted in time-baseline (TB) order.  This is done with AIPS task UVSRT. 
2. Table contents may also be in the incorrect order; this is fixed using
function MergeCal.MergeCal, which employs AIPS task TAMRG. 
3. Finally, an NX file must be generated using AIPS task INDXR.
"""

import OTObit, AIPSTask, MergeCal, UV

def IDIFix( uv, err ):
    """
Fix VLBA data loaded from IDI files.

* uv = UV data set to be sorted
* err = OErr object

Returns the corrected UV data object.
    """
    # Sort data
    uvsrt = AIPSTask.AIPSTask('uvsrt')
    OTObit.setname( uv, uvsrt )
    uvsrt.outname = uvsrt.inname
    uvsrt.outclass = uvsrt.inclass
    uvsrt.outseq = uvsrt.inseq + 1
    uvsrt.outdisk = uvsrt.indisk
    uvsrt.sort = 'TB'
    uvsrt.go()

    # Get output UV data object
    uvFixed = UV.newPAUV( "uvFixed", uvsrt.outname, uvsrt.outclass, 
        int(uvsrt.outdisk), int(uvsrt.outseq), True, err)

    # Fix tables
    MergeCal.MergeCal( uvFixed, err )

    # Create NX table
    indxr = AIPSTask.AIPSTask('indxr')
    OTObit.setname( uvFixed, indxr )
    indxr.go()

    return uvFixed
