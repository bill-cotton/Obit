""" Python Obit AIPS directory utilities
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004,2007
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

def ehex(n, width=0, padding=None):
    """Convert a number into "extended hex".

    Returns the extended hex presentation for N, optionally padding it
    up to WIDTH with PADDING."""

    ehex_digits = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    result = ''

    while n > 0:
        result = ehex_digits[n % len(ehex_digits)] + result
        n /= len(ehex_digits)
        width -= 1
        continue

    # Pad if requested to do so.
    if padding != None:
        while width > 0:
            result = str(padding) + result
            width -= 1
            continue
        pass

    return result
# end ehex

# Python interface to AIPS directory utilities
import Obit, OSystem, OErr, pydoc, string, re, os, Image, UV
global AIPSdisks, nAIPS
nAIPS = 0
AIPSdisks = []
# Some routines need ObitTalk AIPS info
try:
    import AIPS
    # Get AIPS disks Info from ObitTalk
    nAIPS = len(AIPS.disks)-1
    for i in range(1,nAIPS+1):
        AIPSdisks.append(AIPS.disks[i].dirname)
except:
    # Lookup disk info
    for disk in xrange(1, 35):
        area = 'DA' + ehex(disk, 2, '0')
        dir = os.getenv(area)
        if dir:
            AIPSdisks.append(dir)
            nAIPS = len(AIPSdisks)
else:
    pass

# Get list of AIPS disks
# Catalog types
AllType = 0
MAType  = 1
UVType  = 2

def Amatch(tmpl, strn):
    """ Test for AIPS match of string

    A "?" matches one of any character, "*" matches any string
    all blanks matches anything
    Returns True or False
    tmpl     = Template string
    strng    = String to search for occurance of tmpl
    """
    # All blank?
    if re.match("^ *$",tmpl):
        return True
    # Full test
    t   = re.escape(tmpl.strip())
    tt  = t.replace("\\?",".")
    ttt = tt.replace("\\*",".*")
    pat = "^"+ttt+"$"
    m = re.match(pat,strn.strip())
    return m!=None
    # end Amatch

def WantDir(line, type='  ', Aname=None, Aclass=None, Aseq=0):
    """ Test if Catalog entry desired

    Compare PInfo catalog entry to see if it's selected
    Strings use AIPS wild cards:
        blank => any
        '?'   => one of any character
        "*"   => arbitrary string
    Returns True or False
    line     = catalog entry description from PInfo
    type     = AIPS data type '  ' (all) 'MA' Image, 'UV' UV data
    Aname     = desired name, using AIPS wildcards, None -> don't check
    Aclass    = desired class, using AIPS wildcards, None -> don't check
    Aseq      = desired sequence, 0=> any
    """
    # Valid entry?
    if line==None:
        return False
    # Type match
    if (type!="  ") and (type!=line[26:28]) :
        return False
    # Sequence match
    if (Aseq>0) and (Aseq!=int(line[20:25])):
        return False
    # Name and class defaults
    want = True
    if Aname:
        want = Amatch(Aname, line[0:12])
    if Aclass and want:
        want = Amatch(Aclass, line[13:19])
    return want
    # end WantDir

def PFindCNO(disk, user, Aname, Aclass, Atype, seq, err):
    """ Lookup AIPS catalog slot number

    returns AIPS cno
    disk     = AIPS disk number
    user     = AIPS user number
    Aname    = AIPS file name
    Aclass   = AIPS class name
    Atype    = 'MA' or 'UV' for image or uv data
    seq      = AIPS sequence number
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.AIPSDirFindCNO(disk, user, Aname, Aclass, Atype, seq, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "AIPS Catalog error")
    return ret
    # end PFindCNO

def PTestCNO(disk, user, Aname, Aclass, Atype, seq, err):
    """ Test if AIPS file exists

    returns AIPS cno, -1 => not found
    disk     = AIPS disk number
    user     = AIPS user number
    Aname    = AIPS file name
    Aclass   = AIPS class name
    Atype    = 'MA' or 'UV' for image or uv data
    seq      = AIPS sequence number
    err      = Python Obit Error/message stack, 
    """
    ################################################################
    # Checks
    if err.isErr:
        return -1
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Print message stack to clear
    OErr.printErr(err)
    ret = Obit.AIPSDirFindCNO(disk, user, Aname, Aclass, Atype, seq, err.me)
    if err.isErr:
        return ret
    # Clear any couldn't find message
    OErr.PClear(err)
    return ret
    # end PTestCNO

def PHiSeq(disk, user, Aname, Aclass, Atype, err):
    """ find highest sequence number for AIPS name, class...

    returns highest sequence number, -1 => not found
    disk     = AIPS disk number
    user     = AIPS user number
    Aname    = AIPS file name
    Aclass   = AIPS class name
    Atype    = 'MA' or 'UV' for image or uv data
    err      = Python Obit Error/message stack, 
    """
    ################################################################
    # Checks
    if err.isErr:
        return -1
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.AIPSDirHiSeq(disk, user, Aname, Aclass, Atype, err.me)
    return ret
    # end PHiSeq


def PAlloc(disk, user, Aname, Aclass, Atype, seq, err):
    """ Allocate AIPS catalog slot number

    returns AIPS cno
    disk     = AIPS disk number
    user     = AIPS user number
    Aname    = AIPS file name
    Aclass   = AIPS class name
    Atype    = 'MA' or 'UV' for image or uv data
    seq      = AIPS sequence number
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.AIPSDirAlloc(disk, user, Aname, Aclass, Atype, seq, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "AIPS catalog allocation error")
    return ret
    # end PAlloc

def PRemove(disk, user, cno, err):
    """ Deallocate AIPS catalog slot number

    disk     = AIPS disk number
    user     = AIPS user number
    cno      = AIPS catalog slot to deassign
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.AIPSDirRemoveEntry(disk, user, cno, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error removing AIPS catalog entry")
    return ret
    # end PRemove
    
def PNumber(disk, user, err):
    """ Return highest current allowed AIPS catalog Slot number

    disk     = AIPS disk number
    user     = AIPS user number
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.AIPSDirNumber(disk, user, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error obtaining highest AIPS catalog number")
    return ret
    # end PNumber

def PInfo(disk, user, cno, err):
    """ Get information for a  AIPS catalog slot number

    Returned string s:
    Aname  = s[0:12]
    Aclass = s[13:19]
    Aseq   = int(s[20:25])
    Atype  = s[26:28]

    returns string describing entry, None=>empty
    disk     = AIPS disk number
    user     = AIPS user number
    cno      = AIPS catalog slot to deassign
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.AIPSDirInfo(disk, user, cno, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error reading AIPS catalog header")
    return ret
    # end PInfo

def PSetDir(disk, newName, err, URL=None):
    """ Set the directory name for a given AIPS directory

    returns disk number actually used
    disk     = AIPS disk number, <=0 => assign
    newName  = new directory path
    err      = Python Obit Error/message stack
    URL      = URL if on a remote host (Only if using OTObit/ParselTongue)
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    retDisk = Obit.AIPSSetDirname(disk, newName, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error setting AIPS directory name")
    AIPSdisks.append(newName)
    nAIPS = len(AIPSdisks)
    if (disk<=0):
        try:
            AIPS.AIPS.disks.append(AIPS.AIPSDisk(URL, retDisk, newName))
        except:
            pass
        else:
            pass
    return retDisk;
    # end PSetDir


def PListDir(disk, err, type = "  ", first=1, last=1000,
             Aname=None, Aclass=None, Aseq=0, giveList=False):
    """ List AIPS directory

    Entries can be selected using Aname, Aclass, Aseq, using AIPS wildcards
    A "?" matches one of any character, "*" matches any string
    all blanks matches anything
    If giveList then return list of CNOs
    disk     = AIPS disk number
    err      = Python Obit Error/message stack
    type     = optional file type
    Aname    = desired name, using AIPS wildcards, None -> don't check
    Aclass   = desired class, using AIPS wildcards, None -> don't check
    Aseq     = desired sequence, 0=> any
    first    = optional first slot number (1-rel)
    last     = optional last slot number
    giveList = If true, return list of CNOs matching
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    user = OSystem.PGetAIPSuser()
    ncno = PNumber (disk, user, err)
    OErr.printErrMsg(err, "Error getting number of cnos")
    # Init output
    if giveList:
        olist = []
    else:
        olist = None

    mincno = first;
    maxcno = min (ncno, last)
    dirlist = "AIPS Directory listing for disk "+str(disk)+"\n"
    for cno in range(mincno, maxcno):
        line=PInfo(disk, user, cno, err);
        if WantDir(line, type=type, Aname=Aname, Aclass=Aclass, Aseq=Aseq):
            OErr.printErrMsg(err, "Error reading entry")
            dirlist = dirlist+string.rjust(str(cno),3)+" "+line+"\n"
            if giveList:
                olist.append(cno)
    if err.isErr:
        OErr.printErrMsg(err, "Error getting AIPS directory listing")
    # User pager 
    pydoc.ttypager(dirlist)
    return olist
    # end PListDir

def PListCat(cat, disk, type = "  ", first=1, last=1000,
             Aname=None, Aclass=None, Aseq=0, giveList=False):
    """ List AIPS directory given as entries in cat

    Entries can be selected using Aname, Aclass, Aseq, using AIPS wildcards
    A "?" matches one of any character, "*" matches any string
    all blanks matches anything
    If giveList then return list of CNOs
    cat      = list of catalog entries as (cno,s)
               s is a string consisting of:
               Aname  = s[0:12]
               Aclass = s[13:19]
               Aseq   = int(s[20:25])
               Atype  = s[26:28]
    disk     = AIPS disk number
    type     = optional file type
    Aname    = desired name, using AIPS wildcards, None -> don't check
    Aclass   = desired class, using AIPS wildcards, None -> don't check
    Aseq     = desired sequence, 0=> any
    first    = optional first slot number (1-rel)
    last     = optional last slot number
    giveList = If true, return list of CNOs matching (no terminal output)
    """
    ################################################################
    # Init output
    if giveList:
        olist = []
    else:
        olist = None
    ncno = len(cat)
    mincno = first;
    maxcno = min (ncno, last)
    dirlist = "AIPS Directory listing for disk "+str(disk)+"\n"
    for (cno,line) in cat:
        if WantDir(line, type=type, Aname=Aname, Aclass=Aclass, Aseq=Aseq):
            dirlist = dirlist+string.rjust(str(cno),3)+" "+line+"\n"
            if giveList:
                olist.append(cno)
    # Use pager if giveList is False
    if not giveList:
        pydoc.ttypager(dirlist)
    return olist
    # end PListCat

def PAllDest(disk, err, Atype = "  ",  Aname=None, Aclass=None, Aseq=0):
    """ Delete selected AIPS catalog entries

    Entries can be selected using Atype, Aname, Aclass, Aseq,
    using AIPS wildcards
    A "?" matches one of any character, "*" matches any string
    all blanks matches anything
    disk     = AIPS disk number, 0=>all
    err      = Python Obit Error/message stack
    type     = optional file type
    Aname    = desired name, using AIPS wildcards, None -> don't check
    Aclass   = desired class, using AIPS wildcards, None -> don't check
    Aseq     = desired sequence, 0=> any
    first    = optional first slot number (1-rel)
    last     = optional last slot number
    giveList = If true, return list of CNOs matching
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    if (disk<0) or (disk>(len(AIPS.AIPS.disks)+1)):
        raise RuntimeError,"Disk "+str(disk)+" out of range"
    # disks
    if disk>0:
        disks=[disk]
    else:
        disks= range(1,len(AIPS.AIPS.disks))
    user = OSystem.PGetAIPSuser()
    # loop over disks
    for dsk in disks:
        ncno = PNumber (dsk, user, err)
        OErr.printErrMsg(err, "Error getting number of cnos")
        
        for cno in range(1, ncno):
            line=PInfo(dsk, user, cno, err);
            OErr.printErrMsg(err, "Error reading entry")
            if WantDir(line, type=Atype, Aname=Aname, Aclass=Aclass,
                       Aseq=Aseq):
                # parse directory string
                Tname = line[0:12]
                Tclass = line[13:19]
                Tseq = int(line[20:25])
                Ttype = line[26:28]
                z = None
                if Ttype == 'MA':
                    z = Image.newPAImage("Zap image", Tname, Tclass, dsk, Tseq, True, err, \
                               verbose=False)
                    print "Zap AIPS Image",Tname, Tclass, dsk, Tseq
                elif Ttype == 'UV':
                    z = UV.newPAUV("Zap UV data", Tname, Tclass, dsk, Tseq, True, err, \
                                   verbose=False)
                    print "Zap AIPS UV",Tname, Tclass, dsk, Tseq
                # Delete entry
                if z:
                    z.Zap(err)
                    del z
        # end loop over cno
    # end loop over disk
# end PAllDest

