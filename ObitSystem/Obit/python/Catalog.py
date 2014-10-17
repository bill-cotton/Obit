"""  
Catalog manipulation utilities

This module contains the python interfaces to  utility routines to manupulate
ObitTableMF, ObitTableVL and ObitTableFS tables containing source catalogs.
These tables can be generated using Obit Task FndSou or the AIPS tasks
VSAD or SAD.  FS Tables also contain spectra.
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2006-2014
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

# Python utility package for source catalog manipulation
import Obit, Image, ImageDesc, Table, InfoList, OErr, History, SkyGeom

def PMF2VL(MFTable, VLTable, image, err):
    """
    Convert contents in an ObitTableMF to entries in an VL table
    
    Converts an MF table produced by FndSou (or AIPS tasks SAD, VSAD)
    to a VL table.

    * MFTable   = Input Obit MF Table
    * VLTable   = Output Obit VL Table
    * image     = Image being cataloged
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(MFTable):
        print "Actually ",MFTable.__class__
        raise TypeError,"MFTable MUST be a Python Obit Table"
    if not Table.PIsA(VLTable):
        print "Actually ",VLTable.__class__
        raise TypeError,"VLTable MUST be a Python Obit Table"
    if not Image.PIsA(image):
        print "Actually ",image.__class__
        raise TypeError,"image MUST be a Python Obit Image"
    #
    Obit.TableMF2VL (MFTable.me, VLTable.me, image.me, err.me);
    # end PMF2VL

def PMFPrint (MFTable, image, err, file="stdout"):
    """
    Write human readable version of an MF table to a file
    
    Print contents of MF table

    * MFTable   = Input Obit MF Table
    * image     = Image being cataloged
    * err       = Python Obit Error/message stack
    * file      = Name of a file or "stdout" for terminal display
    """
    ################################################################
    # Checks
    if not Table.PIsA(MFTable):
        print "Actually ",MFTable.__class__
        raise TypeError,"MFTable MUST be a Python Obit Table"
    if not Image.PIsA(image):
        print "Actually ",image.__class__
        raise TypeError,"image MUST be a Python Obit Image"
    #
    Obit.TableMFPrint (MFTable.me, image.me, file, err.me);
    # end PMFPrint

def PVLAppend (inVL, outVL, err):
    """
    Append contents of one VL table to another

    * inVL      = Input Obit VL Table
    * outVL     = Output Obit VL Table
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(inVL):
        print "Actually ",inVL.__class__
        raise TypeError,"inVL MUST be a Python Obit Table"
    if not Table.PIsA(outVL):
        print "Actually ",outVL.__class__
        raise TypeError,"VLTable MUST be a Python Obit Table"
    #
    Obit.TableVLAppend (inVL.me, outVL.me, err.me);
    # end PVLAppend 
    
def PVLIndex (inVL,err):
    """
    Index a VL table
    
    Create an index for faster access to contents

    * inVL     = Output Obit VL Table
    * err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(inVL):
        print "Actually ",inVL.__class__
        raise TypeError,"inVL MUST be a Python Obit Table"
    #
    Obit.TableVLIndex (inVL.me, err.me);
    # end PVLIndex 
    
def PVLMerge (inVL, err, Radius = 3, Cutoff=0.05):
    """
    Merge contents of VL table
    
    Merge overlapping components from a given field
    Sums the flux of all components within a specified distance 
    of each component and then merges components weaker than cutoff
    of the total with the strongest component.  
    The resultant position is the weighted average and the 
    flux is the sum.
    This should be run on a table with all entries derived from 
    the same image.

    * inVL      = Input Obit VL Table
    * err       = Python Obit Error/message stack
    * Radius    = The radius in pixels within which to sum components.
    * Cutoff    = The minimum acceptable fraction  of summed flux.
    """
    ################################################################
    # Checks
    if not Table.PIsA(inVL):
        print "Actually ",inVL.__class__
        raise TypeError,"inVL MUST be a Python Obit Table"
    #
    # Set control values
    inInfo = inVL.List
    dim = [1,1,1,1,1]
    InfoList.PPutFloat   (inInfo, "Radius", dim, [Radius], err)
    InfoList.PPutFloat   (inInfo, "Cutoff", dim, [Cutoff], err)
    #
    Obit.TableVLMerge (inVL.me, err.me);
    # end PVLMerge 

def PVLSelect (inVL, outVL, err, BLC=[10,10], TRC=[100000,100000], \
               Steps=[[2.0,0.05], [3.0,0.02], [4.0,0.01]]):
    """
    Select significant components in a VL table
    
    Select significant components in a table to out
    Given a table of radii and fractional fluxes, filter out
    weak sources in the presence of strong.  For any non zero
    entry in Steps, if the entry has a peak weaker than fraction 
    of the total within radius then it is not copied to the output.
    This should be run on a table with all entries derived from 
    the same image.

    * inVL      = Input Obit VL Table
    * outVL     = Output Obit VL Table
    * err       = Python Obit Error/message stack
    * BLC       = Lowest x,y pixel number selected
    * TRC       = Highest pixel number selected
    * Steps     = Pairs of values giving (radius(cells), fraction),
      empty entries zero filled
    """
    ################################################################
    # Checks
    if not Table.PIsA(inVL):
        print "Actually ",inVL.__class__
        raise TypeError,"inVL MUST be a Python Obit Table"
    if not Table.PIsA(outVL):
        print "Actually ",outVL.__class__
        raise TypeError,"VLTable MUST be a Python Obit Table"
    #
    # Set control values
    inInfo = inVL.List
    dim = [2,1,1,1,1]
    InfoList.PPutInt   (inInfo, "BLC", dim, BLC, err)
    InfoList.PPutInt   (inInfo, "TRC", dim, TRC, err)
    # Convert Steps to 1-D
    sss = []
    for s in Steps:
        sss.append(s[0])
        sss.append(s[1])
    dim = [2, len(Steps),1,1,1]
    InfoList.PPutFloat   (inInfo, "Steps", dim, sss, err)
    #
    Obit.TableVLSelect (inVL.me, outVL.me, err.me);
    # end PVLSelect 

def PVLPurge (inVL, field, err):
    """
    Remove VL table entries from a given field
    
    Entries from "field" are deselected

    * inVL      = Input Obit VL Table
    * field     = Name of field to purge
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(inVL):
        print "Actually ",inVL.__class__
        raise TypeError,"inVL MUST be a Python Obit Table"
    #
    Obit.TableVLPurge (inVL.me, outVL.me, err.me);
    # end PVLPurge 

def PVLRedun (inVL, outVL, err, centerPix=[512,512], maxDist=15.0):
    """
    Remove redundant entries from a VL table
    
    Remove redundant entries from in and write out
    Search forward from each entry until past time of possible match.
    If a positional match (maxDist) is found then the one closest
    to the center of its image (centerPix) is chosen.  The other
    entry is flagged in the input table.  When the search over the
    RA range from the original source is finished the final accepted
    entry is written to the output table.

    * inVL      = Input Obit VL Table
    * outVL     = Output Obit VL Table
    * err       = Python Obit Error/message stack
    * centerPix = Center pixel in images
    * maxDist   = How far (asec) to search for matches
    """
    ################################################################
    # Checks
    if not Table.PIsA(inVL):
        print "Actually ",inVL.__class__
        raise TypeError,"inVL MUST be a Python Obit Table"
    if not Table.PIsA(outVL):
        print "Actually ",outVL.__class__
        raise TypeError,"VLTable MUST be a Python Obit Table"
    # Set control values
    inInfo = inVL.List
    dim = [1,1,1,1,1]
    InfoList.PPutFloat (inInfo, "maxDist", dim, [maxDist], err)
    dim = [2,1,1,1,1]
    InfoList.PPutInt   (inInfo, "centerPix", dim, centerPix, err)
    #
    Obit.TableVLRedun (inVL.me, outVL.me, err.me);
    # end PVLRedun 

def PVLFilter (inVL, outVL, err, \
                   col1="PEAK INT", col2="I RMS", minRatio=5.0):
    """
    Select entries in a VL table based on the ratio of the values
    
    Copy table rows whose ratio of the values in two columns exceed minRatio
    value(col1)/value(col2) > minRatio are copied

    * inVL      = Input Obit VL Table
    * outVL     = Output Obit VL Table, new data appended to end of old
    * err       = Python Obit Error/message stack
    * col1      = label of first column
    * col2      = label of second column
    * minRatio  = minimum acceptable ratio
    """
    ################################################################
    # Checks
    if not Table.PIsA(inVL):
        print "Actually ",inVL.__class__
        raise TypeError,"inVL MUST be a Python Obit Table"
    if not Table.PIsA(outVL):
        print "Actually ",outVL.__class__
        raise TypeError,"VLTable MUST be a Python Obit Table"
    inVL.Open(Table.READONLY,err)
    outVL.Open(Table.READWRITE,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error Opening tables")
    nrow = inVL.Desc.Dict['nrow']
    orow = -1;
    for irow in range(1,nrow+1):
        row = inVL.ReadRow(irow,err)
        ratio = row[col1][0]/row[col2][0]
        if ratio>minRatio:
            outVL.WriteRow(orow, row, err)
            if err.isErr:
                break
    inVL.Close(err)
    outVL.Close(err)

    # end PVLFilter 

def PVLScale (inImage, inVL, outVL, scale, err,
              exclude=None):
    """
    Scale the flux density values in a VL table
    
    Copies beam parameters and indexes when done.
    * inImage   = Input image with VL tables
    * inVL      = Input Obit VL Table
    * outVL     = Output Obit VL Table, new data appended to end of old
    * scale     = scale factor for flux densities
    * err       = Python Obit Error/message stack
    * exclude   = list of tuples with (ra, deg, size) all in deg
                  describing regions not to rescale
    """
    ################################################################
    # Checks
    if not Table.PIsA(inVL):
        print "Actually ",inVL.__class__
        raise TypeError,"inVL MUST be a Python Obit Table"
    if not Table.PIsA(outVL):
        print "Actually ",outVL.__class__
        raise TypeError,"VLTable MUST be a Python Obit Table"
    fblank = Obit.FArrayGetBlank ()  # Blanked value
    inVL.Open(Table.READONLY,err)
    outVL.Open(Table.READWRITE,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error Opening tables")
    # Copy beam info
    outVL.keys["BeamMajor"] = inVL.keys["BeamMajor"]
    outVL.keys["BeamMinor"] = inVL.keys["BeamMinor"]
    outVL.keys["BeamPA"]    = inVL.keys["BeamPA"]
    #print "debug",outVL.keys
    nrow = inVL.Desc.Dict['nrow']
    orow = -1;
    for irow in range(1,nrow+1):
        row = inVL.ReadRow(irow,err)
        if err.isErr:
            break
        doScale = True    # Need to scale?
        if exclude:
            for pos in exclude:
                (dra,ddec) = SkyGeom.PShiftXY(row["RA(2000)"][0], row["DEC(2000)"][0], 0.0, pos[0], pos[1])
                if (abs(dra)<=pos[2]) and (abs(ddec)<=pos[2]):
                    doScale = False
                    print "exclude", irow, ImageDesc.PRA2HMS(row["RA(2000)"][0]), ImageDesc.PDec2DMS(row["DEC(2000)"][0]), \
                        row["PEAK INT"][0]
                    break
        if doScale:
            row["PEAK INT"][0] *= scale
            row["I RMS"][0]    *= scale
            row["RES RMS"][0]  *= scale
            row["RES PEAK"][0] *= scale
            row["RES FLUX"][0] *= scale
            # Have poln data?
            if (row["Q CENTER"][0]!=fblank) and (row["U CENTER"][0]!=fblank):
                row["Q CENTER"][0] *= scale
                row["U CENTER"][0] *= scale
                row["P Flux"][0]   *= scale
                row["POL RMS"][0]  *= scale
        # write
        outVL.WriteRow(orow, row, err)
        if err.isErr:
            break
    inVL.Close(err)
    outVL.Close(err)

    # History
    #outVL.Open(Table.READWRITE,err)
    inHistory  = History.History("history", inImage.List, err)
    inHistory.Open(History.READWRITE, err)
    inHistory.TimeStamp("Start PVLScale", err)
    inHistory.WriteRec(-1,"PVLScale / scale = "+str(scale),err)
    if exclude:
        i = 0
        for pos in exclude:
            i += 1
            hi = "PVLScale exclude[%d] = (%9.5f,%9.5f,%8.6f) \ no Scale"%(i,pos[0], pos[1],pos[2])
            inHistory.WriteRec(-1,hi,err)
    inHistory.Close(err)
    #outVL.Close(err)

    # Index
    print "Index output"
    PVLIndex (outVL,err)

    # end PVLScale

def PVLPrint (VLTable, image, err, file="stdout"):
    """
    Write human readable version of an VL table to a file
    
    Print contents of VL table

    * VLTable   = Input Obit VL Table
    * image     = Image being cataloged
    * err       = Python Obit Error/message stack
    * file      = Name of a file or "stdout" for terminal display
    """
    ################################################################
    # Checks
    if not Table.PIsA(VLTable):
        print "Actually ",VLTable.__class__
        raise TypeError,"VLTable MUST be a Python Obit Table"
    if not Image.PIsA(image):
        print "Actually ",image.__class__
        raise TypeError,"image MUST be a Python Obit Image"
    #
    Obit.TableVLPrint (VLTable.me, image.me, file, err.me);
    # end PVLPrint

def PVL2VZ(VLTable, image, err):
    """
    Convert contents in an ObitTableVL to entries in an VZ table
    
    Converts VL table into forme used to specify outliers, calibrators
    Returns VZ table

    * VLTable   = Input Obit VL Table. Control parameters on info object:

      ======  =================  ==================================== 
      "Flux"  OBIT_float scalar  Minimum acceptable flux density (Jy)
      ======  =================  ==================================== 

    * image     = Image to which VZ table to be attached
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(VLTable):
        print "Actually ",VLTable.__class__
        raise TypeError,"VLTable MUST be a Python Obit Table"
    if not Image.PIsA(image):
        print "Actually ",image.__class__
        raise TypeError,"image MUST be a Python Obit Image"
    #
    # Get data
    id =  Obit.ImageCastData(image.me)
    outVZ  = Table.Table("None")
    outVZ.me = Obit.TableVL2VZ (VLTable.me, id, err.me);
    return outVZ
    # end PVL2VZ

def PVZSel(VZTable, image, err):
    """
    Select entries in VZ table
    
    Select and possibly average entries in a VZ table and determine
    crowding quality code.
    Returns new VZ table

    * VZTable   = Input Obit VZ Table.

        Control parameters on info object:

        ========= =================== =========================================
        "clip"    OBIT_float scalar   Minimum acceptable for clipping (Jy)
                                      (used if doClip, no default)
        "nearest" OBIT_float scalar   Minimum nearest neighbor to be 
                                      considered isolated (deg) (used if doIsol,
                                      doAver, no default)
        "distAvg" OBIT_float scalar   Distance within to average sources (deg)
                                      (used if doAver, no default)
        "ignore"  OBIT_float scalar   Minimum acceptable flux density (Jy)
                                      (used if doAver, no default)
        "crowd"   OBIT_float scalar   Distance (deg) for crowding test
                                      (used if doAver doIsol, no default).  The
                                      flux of sources inside of crowd are summed
                                      and the crowing quality code is determined
                                      from the ratio of the sum of the (averaged)
                                      flux density of the source to the sum
                                      inside the field.

                                      * 0 = nothing else within crowd
                                      * 1 = >10% of peak other flux within crowd
                                      * 2 = >30% of peak other flux within crowd
                                      * 3 = >70% of peak other flux within crowd
                                      * 4 = >100% of peak other flux within crowd
                                      * 5 = >200% of peak other flux within crowd
                                    
        "doClip"  OBIT_boolean scalar Clipping by minPeak [def FALSE]
        "doIsol"  OBIT_boolean scalar Select isolated entries within nearest 
                                      ignore entries fainter than ignore [def FALSE]
        "doAver"  OBIT_boolean scalar Average entries within distAvg
                                      reject outliers [def FALSE]

                                      * If the brightest source within *nearest* is fainter than 
                                        *ignore* then the source is passed, else:
                                      * If a component is the brightest within *nearest* 
                                        then all components within *distAvg* are merged.
                                      * If the source has no companions within nearest 
                                        it is passed.
                                      * If the component has brighter companions within *nearest*,
                                        it is dropped.
        ========= =================== =========================================

    * image     = Image to which VZ table to be attached
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(VZTable):
        print "Actually ",VZTable.__class__
        raise TypeError,"VZTable MUST be a Python Obit Table"
    if not Image.PIsA(image):
        print "Actually ",image.__class__
        raise TypeError,"image MUST be a Python Obit Image"
    #
    # Get data
    id =  Obit.ImageCastData(image.me)
    outVZ  = Table.Table("None")
    outVZ.me = Obit.TableVZSel (VZTable.me, id, err.me);
    return outVZ
    # end PVZSel

def PVL2FS(VLTable, image, FSver, err):
    """
    Convert contents in an ObitTableVL to entries in an FS table
    
    Converts an VL table produced by FndSou (or AIPS tasks SAD, VSAD)
    to a FS table (allows spectra).

    * VLTable   = Input Obit VL Table, control parameters on List:
                  minFlux = Minimum acceptable flux density (Jy) [def 0]
    * image     = Image cube being cataloged
    * FSver     = FS version
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(VLTable):
        print "Actually ",VLTable.__class__
        raise TypeError,"VLTable MUST be a Python Obit Table"
    if not Image.PIsA(image):
        print "Actually ",image.__class__
        raise TypeError,"image MUST be a Python Obit Image"
    #
    Obit.TableVL2FS (VLTable.me, image.me, FSver, err.me);
    # end PVL2FS

def PFSAppend (inFS, outFS, err):
    """
    Append contents of one FS table to another

    * inFS      = Input Obit FS Table
    * outFS     = Output Obit FS Table
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(inFS):
        print "Actually ",inFS.__class__
        raise TypeError,"inFS MUST be a Python Obit Table"
    if not Table.PIsA(outFS):
        print "Actually ",outFS.__class__
        raise TypeError,"FSTable MUST be a Python Obit Table"
    #
    Obit.TableFSAppend (inFS.me, outFS.me, err.me);
    # end PFSAppend 
    
def PFSIndex (inFS,err):
    """
    Index a FS table
    
    Create an index for faster access to contents

    * inFS     = Output Obit FS Table
    * err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(inFS):
        print "Actually ",inFS.__class__
        raise TypeError,"inFS MUST be a Python Obit Table"
    #
    Obit.TableFSIndex (inFS.me, err.me);
    # end PFSIndex 
    
def PFSMerge (inFS, err, Radius = 3, Cutoff=0.05):
    """
    Merge contents of FS table
    
    Merge overlapping components from a given field
    Sums the flux of all components within a specified distance 
    of each component and then merges components weaker than cutoff
    of the total with the strongest component.  
    The resultant position is the weighted average and the 
    flux is the sum.
    This should be run on a table with all entries derived from 
    the same image.

    * inFS      = Input Obit FS Table, additional control parameters on List:
                  begRow = 1st (1-rel) row [def 1]
                  endRow = last (1-rel) row [def all]
    * err       = Python Obit Error/message stack
    * Radius    = The radius in pixels within which to sum components.
    * Cutoff    = The minimum acceptable fraction  of summed flux.
    """
    ################################################################
    # Checks
    if not Table.PIsA(inFS):
        print "Actually ",inFS.__class__
        raise TypeError,"inFS MUST be a Python Obit Table"
    #
    # Set control values
    inInfo = inFS.List
    dim = [1,1,1,1,1]
    InfoList.PPutFloat   (inInfo, "Radius", dim, [Radius], err)
    InfoList.PPutFloat   (inInfo, "Cutoff", dim, [Cutoff], err)
    #
    Obit.TableFSMerge (inFS.me, err.me);
    # end PFSMerge 

def PFSSelect (inFS, outFS, err, BLC=[10,10], TRC=[100000,100000], \
               Steps=[[2.0,0.05], [3.0,0.02], [4.0,0.01]]):
    """
    Select significant components in a FS table
    
    Select significant components in a table to out
    Given a table of radii and fractional fluxes, filter out
    weak sources in the presence of strong.  For any non zero
    entry in Steps, if the entry has a peak weaker than fraction 
    of the total within radius then it is not copied to the output.
    This should be run on a table with all entries derived from 
    the same image.

    * inFS      = Input Obit FS Table, additional control parameters on List:
                  begRow = 1st (1-rel) row [def 1]
                  endRow = last (1-rel) row [def all]
    * outFS     = Output Obit FS Table
    * err       = Python Obit Error/message stack
    * BLC       = Lowest x,y pixel number selected
    * TRC       = Highest pixel number selected
    * Steps     = Pairs of values giving (radius(cells), fraction),
      empty entries zero filled
    """
    ################################################################
    # Checks
    if not Table.PIsA(inFS):
        print "Actually ",inFS.__class__
        raise TypeError,"inFS MUST be a Python Obit Table"
    if not Table.PIsA(outFS):
        print "Actually ",outFS.__class__
        raise TypeError,"FSTable MUST be a Python Obit Table"
    #
    # Set control values
    inInfo = inFS.List
    dim = [2,1,1,1,1]
    InfoList.PPutInt   (inInfo, "BLC", dim, BLC, err)
    InfoList.PPutInt   (inInfo, "TRC", dim, TRC, err)
    # Convert Steps to 1-D
    sss = []
    for s in Steps:
        sss.append(s[0])
        sss.append(s[1])
    dim = [2, len(Steps),1,1,1]
    InfoList.PPutFloat   (inInfo, "Steps", dim, sss, err)
    #
    Obit.TableFSSelect (inFS.me, outFS.me, err.me);
    # end PFSSelect 

def PFSPurge (inFS, field, err):
    """
    Remove FS table entries from a given field
    
    Entries from "field" are deselected

    * inFS      = Input Obit FS Table
    * field     = Name of field to purge
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(inFS):
        print "Actually ",inFS.__class__
        raise TypeError,"inFS MUST be a Python Obit Table"
    #
    Obit.TableFSPurge (inFS.me, outFS.me, err.me);
    # end PFSPurge 

def PFSRedun (inFS, outFS, err, centerPix=[512,512], maxDist=15.0):
    """
    Remove redundant entries from a FS table
    
    Search forward from each entry until past time of possible match.
    If a positional match (maxDist) is found then the one closest
    to the center of its image (centerPix) is chosen.  The other
    entry is flagged in the input table.  When the search over the
    RA range from the original source is finished the final accepted
    entry is written to the output table.

    * inFS      = Input Obit FS Table, additional control parameters on List:
                  begRow = 1st (1-rel) row [def 1]
                  endRow = last (1-rel) row [def all]
    * outFS     = Output Obit FS Table
    * err       = Python Obit Error/message stack
    * centerPix = Center pixel in images
    * maxDist   = How far (asec) to search for matches
    """
    ################################################################
    # Checks
    if not Table.PIsA(inFS):
        print "Actually ",inFS.__class__
        raise TypeError,"inFS MUST be a Python Obit Table"
    if not Table.PIsA(outFS):
        print "Actually ",outFS.__class__
        raise TypeError,"FSTable MUST be a Python Obit Table"
    # Set control values
    inInfo = inFS.List
    dim = [1,1,1,1,1]
    InfoList.PPutFloat (inInfo, "maxDist", dim, [maxDist], err)
    dim = [2,1,1,1,1]
    InfoList.PPutInt   (inInfo, "centerPix", dim, centerPix, err)
    #
    Obit.TableFSRedun (inFS.me, outFS.me, err.me);
    # end PFSRedun 

def PFSPrint (FSTable, image, err, file="stdout"):
    """
    Write human readable version of an FS table to a file
    
    Print contents of FS table

    * FSTable   = Input Obit FS Table, List controls:
    *             minSNR = min SNR to select
    * image     = Image being cataloged
    * err       = Python Obit Error/message stack
    * file      = Name of a file or "stdout" for terminal display
    """
    ################################################################
    # Checks
    if not Table.PIsA(FSTable):
        print "Actually ",FSTable.__class__
        raise TypeError,"FSTable MUST be a Python Obit Table"
    if not Image.PIsA(image):
        print "Actually ",image.__class__
        raise TypeError,"image MUST be a Python Obit Image"
    #
    Obit.TableFSPrint (FSTable.me, image.me, file, err.me);
    # end PFSPrint

def PFSGetSpectrum (inFS, image, err):
    """
    Extract spectra for FS table entries from image
    
    * inFS      = Input Obit FS Table
    * image     = Image being cataloged, control parameters on List:
                  RMSsize =  halfwidth of region to determine RMS [def 50]
                             <0 => use full image
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(inFS):
        print "Actually ",inFS.__class__
        raise TypeError,"inFS MUST be a Python Obit Table"
    if not Image.PIsA(image):
        print "Actually ",image.__class__
        raise TypeError,"image MUST be a Python Obit Image"
    #
    Obit.TableFSGetSpectrum (FSTable.me, image.me, err.me);
    # end PFSGetSpectrum

def PFSFiltVel (inFS, image, outFS, err):
    """
    Filter and fit spectra
    
    * inFS      = Input Obit FS Table
    * image     = Image cube being cataloged, control parameters on List:
                  minSNR = Minimum acceptable SNR [def 5]
                  minFlux = Minimum acceptable flux density (Jy) [def 0]
    * outFS     = output Obit FS Table
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(inFS):
        print "Actually ",inFS.__class__
        raise TypeError,"inFS MUST be a Python Obit Table"
    if not Table.PIsA(outFS):
        print "Actually ",outFS.__class__
        raise TypeError,"outFS MUST be a Python Obit Table"
    if not Image.PIsA(image):
        print "Actually ",image.__class__
        raise TypeError,"image MUST be a Python Obit Image"
    #
    Obit.TableFSFiltVel(inFS.me, image.me, outFS.me, err.me);
    # end PFSFiltVel


