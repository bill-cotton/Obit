""" Python Obit Interferometric gain solution class

This class solves for interferometer complex gains

UVGSolve Members with python interfaces:
InfoList    - used to pass instructions to processing
Member List  (readonly)
"""
# $Id: UVGSolve.py,v 1.2 2008/01/03 15:31:36 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C)2006,2008
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

# Interferometric selfcal class
import Obit, OErr, InfoList, UV, types

# Python shadow class to ObitUVGSolve class
 
class UVGSolvePtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.UVGSolveUnref(Obit.UVGSolve_me_get(self.this))
            # In with the new
            Obit.UVGSolve_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != UVGSolve:
            return
        if name == "me" : 
            return Obit.UVGSolve_me_get(self.this)
        # Virtual members
        if name=="List":
            return PGetList(self)
        raise AttributeError,str(name)
    def __repr__(self):
        if self.__class__ != UVGSolve:
            return
        return "<C UVGSolve instance> " + Obit.UVGSolveGetName(self.me)
#
class UVGSolve(UVGSolvePtr):
    """ Python Obit Interferometric gain solution class

    This class solves for interferometer complex gains

    UVGSolve Members with python interfaces:
    InfoList  - used to pass instructions to processing
                Member List
                (readonly)
    """
    def __init__(self, name) :
        self.this = Obit.new_UVGSolve(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_UVGSolve(self.this)


def input(inputDict):
    """ Print the contents of an input Dictionary

    inputDict = Python Dictionary containing the parameters for a routine
    There should be a member of the dictionary ('structure') with a value
    being a list containing:
    1) The name for which the input is intended (string)
    2) a list of tuples consisting of (parameter name, doc string)
       with an entry for each parameter in the dictionary.
       The display of the the inputs dictionary will be in the order of
       the tuples and display the doc string after the value.
       An example:
       Soln2CalInput={'structure':['Soln2Cal',[('InData','Input OTF'),
                                               ('soln','input soln table version'),
                                               ('oldCal','input cal table version, -1=none'),
                                               ('newCal','output cal table')]],
                      'InData':None, 'soln':0, 'oldCal':-1, 'newCal':0}
    """
    ################################################################
    structure = inputDict['structure']  # Structure information
    print 'Inputs for ',structure[0]
    for k,v in structure[1]:
        print '  ',k,' = ',inputDict[k],' : ',v
        
    # end input

def newObit(name, err):
    """ Create and initialize an UVGSolve structure

    Create sky model object
    Returns the Python UVGSolve object
    name     = name desired for object (labeling purposes)
    err      = Python Obit Error/message stack
    """
    ################################################################
    out = UVGSolve (name)
    return out      # seems OK
    # end newObit

def PCopy (inUVGSolve, outUVGSolve, err):
    """ Make a shallow copy of input object.

    Makes structure the same as inUVGSolve, copies pointers
    inUVGSolve  = Python UVGSolve object to copy
    outUVGSolve = Output Python UVGSolve object, must be defined
    err         = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inUVGSolve):
        raise TypeError,"inUVGSolve MUST be a Python Obit UVGSolve"
    if not PIsA(outUVGSolve):
        raise TypeError,"outUVGSolve MUST be a Python Obit UVGSolve"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.UVGSolveCopy (inUVGSolve.me, outUVGSolve.me, err.me)
    if err.isErr:
        printErrMsg(err, "Error copying UVGSolve")
    # end PCopy

def PGetList (inUVGSolve):
    """ Return the member InfoList

    returns InfoList
    inUVGSolve  = Python UVGSolve object
    """
    ################################################################
     # Checks
    if not PIsA(inUVGSolve):
        raise TypeError,"inUVGSolve MUST be a Python Obit UVGSolve"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.UVGSolveGetList(inUVGSolve.me)
    return out
    # end PGetList

# Define UVGSolve - most parameters can be defined here
GSolveInput={'structure':['SelfCal',[('subA ','Selected subarray (default 1)'),
                                     ('solInt','Solution interval (min). (default 1 sec)'),
                                     ('refAnt','Ref ant to use. (default 1)'),
                                     ('avgPol','True if RR and LL to be averaged (false)'),
                                     ('avgIF','True if all IFs to be averaged (false)'),
                                     ('minSNR','Minimum acceptable SNR (5)'),
                                     ('doMGM','True then find the mean gain modulus (false)'),
                                     ('solType',"Solution type '  ', 'L1',  (' ')"),
                                     ('solMode',"Solution mode: 'A&P', 'P', 'P!A', (default 'P')"),
                                     ('minNo','Min. no. antennas. (default 4)'),
                                     ('WtUV','Weight outside of UVRANG. (default 0.0)'),
                                     ('antWt','Weight per antenna, def all 1'),
                                     ('prtLv','Print level (default no print)')]],
             # defaults
             'subA':1,
             'solInt':0.01666667,
             'refAnt':1,
             'avgPol':False,
             'avgIF':False,
             'minSNR':5,
             'doMGM':False,
             'solType':'    ',
             'solMode':'    ',
             'minNo':4,
             'WtUV':0.0,
             'prtLv':0}

def PCreate (name, err, input=GSolveInput):
    """ Create the parameters and underlying structures of a UVGSolve.

    Returns UVGSolve created.
    name        = Name to be given to object
                  Most control parameters are in InfoList member
    err         = Python Obit Error/message stack
    input       = control parameters:
    subA        = Selected subarray (default 1)
    solInt      = Solution interval (min). (default 1 sec)
    refAnt      = Ref ant to use. (default 1)
    avgPol      = True if RR and LL to be averaged (false)
    avgIF       = True if all IFs to be averaged (false)
    minSNR      = Minimum acceptable SNR (5)
    doMGM       = True then find the mean gain modulus (true)
    solType     = Solution type '  ', 'L1',  (' ')
    solMode     = Solution mode: 'A&P', 'P', 'P!A', 'GCON' ('P')
    minNo       = Min. no. antennas. (default 4)
    WtUV        = Weight outside of UV_Full. (default 1.0)
    antWt       = Weight per antenna, def all 1
    prtLv       = Print level (default no print)
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Create
    out = UVGSolve(name);
    out.me = Obit.UVGSolveCreate(name, skyModel.me)
    # Set SelfCal control values on out
    inInfo = PGetList(out)    # 
    dim = [1,1,1,1,1]
    # Set control values on SelfCal
    dim[0] = 1;
    inInfo = UV.PGetList(skyModel)    #
    InfoList.PPutInt     (inInfo, "subA",   dim, [input["subA"]],  err)
    InfoList.PPutInt     (inInfo, "refAnt", dim, [input["refAnt"]],err)
    InfoList.PPutInt     (inInfo, "minNo",  dim, [input["minNo"]], err)
    InfoList.PPutInt     (inInfo, "prtLv",  dim, [input["prtLv"]], err)
    InfoList.PPutBoolean (inInfo, "avgPol", dim, [input["avgPol"]],err)
    InfoList.PPutBoolean (inInfo, "avgIF",  dim, [input["avgIF"]], err) 
    InfoList.PPutBoolean (inInfo, "doMGM",  dim, [input["doMGM"]], err)
    InfoList.PPutFloat   (inInfo, "solInt", dim, [input["solInt"]],err)
    InfoList.PPutFloat   (inInfo, "minSNR", dim, [input["minSNR"]],err)
    InfoList.PPutFloat   (inInfo, "WtUV",   dim, [input["WtUV"]],  err)
    dim[0] = len(input["solType"])
    InfoList.PAlwaysPutString  (inInfo, "solType",dim, [input["solType"]])
    dim[0] = len(input["solMode"])
    InfoList.PAlwaysPutString  (inInfo, "solMode",dim, [input["solMode"]])
    # antWt if given 
    if input["antWt"].__class__!=None.__class__:
        dim[0] = len(input["antWt"])
        InfoList.PAlwaysPutFloat (inInfo, "antWt",   dim, input["antWt"])
    # show any errors 
    #OErr.printErrMsg(err, "UVGSolveCreate: Error setting parameters")
    #
    return out;
    # end PCreate

# Self calibrate
UVGSolveCalInput={
    'structure':['UVGSolveCal',[('InData',  'Input UV data'),
                                ('OutData', 'Output UV data'),
                                ('subA ','Selected subarray (default 1)'),
                                ('solInt','Solution interval (min). (default 1 sec)'),
                                ('refAnt','Ref ant to use. (default 1)'),
                                ('avgPol','True if RR and LL to be averaged (false)'),
                                ('avgIF','True if all IFs to be averaged (false)'),
                                ('minSNR','Minimum acceptable SNR (5)'),
                                ('doMGM','True then find the mean gain modulus (false)'),
                                ('solType',"Solution type '  ', 'L1',  (' ')"),
                                ('solMode',"Solution mode: 'A&P', 'P', 'P!A', (default 'P')"),
                                ('minNo','Min. no. antennas. (default 4)'),
                                ('WtUV','Weight outside of UVRANG. (default 0.0)'),
                                ('antWt','Weight per antenna, def all 1'),
                                ('prtLv','Print level (default no print)')]],
    # defaults
    'InData':None,
    'OutData':None,
    'subA':1,
    'solInt':0.0,
    'refAnt':1,
    'avgPol':False,
    'avgIF':False,
    'minSNR':5,
    'doMGM':False,
    'solType':'    ',
    'solMode':'    ',
    'minNo':4,
    'WtUV':0.0,
    'prtLv':0}
def PCal (inUVGSolve, err, input=UVGSolveCalInput):
    """ Solve for gains

    Determine a gain table from a set of data assumed divided by model.
    Solution table attached to  OutData
    inUVGSolve  = Gain solver object
    err         = Python Obit Error/message stack
    input       = input parameter dictionary
    
    Input dictionary entries:
    InData   = Input Python UV data from which to determine calibration
               Should have been divided by model if not a point source.
    OutData  = UV data to which output SN table will be attached.
    subA     = Selected subarray (default 1)
    solInt   = Solution interval (min). (default 1 sec)
    refAnt   = Ref ant to use. (default 1)
    avgPol   = True if RR and LL to be averaged (false)
    avgIF    = True if all IFs to be averaged (false)
    minSNR   = Minimum acceptable SNR (5)
    doMGM    = True then find the mean gain modulus (true)
    solType  = Solution type '  ', 'L1',  (' ')
    solMode  = Solution mode: 'A&P', 'P', 'P!A', 'GCON' ('P')
    minNo    = Min. no. antennas. (default 4)
    WtUV     = Weight outside of UV_Full. (default 1.0)
    antWt    = Weight per antenna, def all 1
    prtLv    = Print level (default no print)
    returns  = Output SN Table with solution
    """
    ################################################################
    # Get input parameters
    InData    = input["InData"]
    OutData   = input["OutData"]
    #
    # Checks
    if not PIsA(inUVGSolve):
        raise TypeError, 'Bad input selfcalibrator'
    if not UV.PIsA(InData):
        raise TypeError, 'Bad input UV data'
    if not UV.PIsA(OutData):
        raise TypeError, 'Bad output UV data'
    dim = [1,1,1,1,1]
    # Set control values on SelfCal
    dim[0] = 1;
    inInfo = PGetList(inUVGSolve)    # 
    InfoList.PPutInt     (inInfo, "subA",   dim, [input["subA"]],  err)
    InfoList.PPutInt     (inInfo, "refAnt", dim, [input["refAnt"]],err)
    InfoList.PPutInt     (inInfo, "minNo",  dim, [input["minNo"]], err)
    InfoList.PPutInt     (inInfo, "prtLv",  dim, [input["prtLv"]], err)
    InfoList.PPutBoolean (inInfo, "avgPol", dim, [input["avgPol"]],err)
    InfoList.PPutBoolean (inInfo, "avgIF",  dim, [input["avgIF"]], err) 
    InfoList.PPutBoolean (inInfo, "doMGM",  dim, [input["doMGM"]], err)
    InfoList.PPutFloat   (inInfo, "solInt", dim, [input["solInt"]],  err)
    InfoList.PPutFloat   (inInfo, "minSNR", dim, [input["minSNR"]],  err)
    InfoList.PPutFloat   (inInfo, "WtUV",   dim, [input["WtUV"]],    err)
    dim[0] = len(input["solType"])
    InfoList.PAlwaysPutString  (inInfo, "solType",dim, [input["solType"]])
    dim[0] = len(input["solMode"])
    InfoList.PAlwaysPutString  (inInfo, "solMode",dim, [input["solMode"]])
    # antWt if given 
    if input["antWt"].__class__!=None.__class__:
        dim[0] = len(input["antWt"])
        InfoList.PAlwaysPutFloat (inInfo, "antWt",   dim, input["antWt"])
    # show any errors 
    OErr.printErrMsg(err, "UVGSolveCreate: Error setting parameters")
    #
    # Create output
    outTable    = Table.Table("None")
    outTable.me = UVGSolveCal(inUVGSolve.me, InData.me, OutData.me, err)
    # show any errors 
    OErr.printErrMsg(err, "UVGSolveCal: Error Self-calibrating")
    #
    return outTable
    # end  PCal

def PRefAnt (SNTab, isuba, refant, err):
    """ rereference phases of subarray isuba in a SN table to refant

    SNTab     = Python AIPS SN Table
    isuba     = subarray [def. 1]
    refant    = reference antenna, if 0 -> pick one
    err       = Python Obit Error/message stack
    Returns actual reference antenna used
    """
    ################################################################
    # Checks
    if not Table.PIsA(SNTab):
        raise TypeError,"SNTab MUST be a Python Obit Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.UVSolnRefAnt (SNTab.me, isuba, refant, err.me)
    if err.isErr:
        printErrMsg(err, "Error rereferencing phases")
    return ret
    # end PRefAnt


# SN table smoothing
UVGSolveSNSmoInput={
    'structure':['UVSelfSNSmo',[('InData',  'Input UV data'),
                                ('InTable', 'Input SN Table'),
                                ('isuba',   'Desired subarray, 0=> 1 '),
                                ('smoType', 'Smoothing type MWF, GAUS or BOX, def BOX'),
                                ('amoAmp',  'Amplitude smoothing time in min'),
                                ('smpPhase','Phase smoothing time in min. (0 => fix failed only')]],
    # defaults
    'InData':None,
    'OutTable':None,
    'isuba':0,
    'smoType':'BOX',
    'smoAmp':0.0,
    'smpPhase':0.0}
def PSNSmo (err, input=UVGSolveSNSmoInput):
    """ Smooth SN table possibly replacing blanked soln.

    err        = Python Obit Error/message stack
    input      = input parameter dictionary
    
    Input dictionary entries:
    InData   = Input Python UV data 
    InTable  = Input SN table
    isuba    = Desired subarray, 0=> 1 
    smoType  = Smoothing type MWF, GAUS or BOX, def BOX
    smoAmp   = Amplitude smoothing time in min
    smpPhase = Phase smoothing time in min. (0 => fix failed only')
    """
    ################################################################
    # Get input parameters
    InData    = input["InData"]
    InTable   = input["InTable"]
    isuba     = input["isuba"]
    #
    # Checks
    if not UV.PIsA(InData):
        raise TypeError, 'PCal: Bad input UV data'
    if not Table.PIsA(InTable):
        raise TypeError, 'PCal: Bad input table'
    # Set control values on UV 
    dim[0] = 1;
    inInfo = UV.PGetList(InData)  # Add control to UV data
    dim[0] = len(input["smoType"]);
    InfoList.PAlwaysPutString  (inInfo, "smoType", dim, [input["smoType"]])
    dim[0] = 1;
    InfoList.PPutFloat   (inInfo, "smoAmp",  dim, [input["smoAmp"]],  err)
    InfoList.PPutFloat   (inInfo, "smoPhase",dim, [input["smoPhase"]],err)
    # Smooth
    Obit.UVSolnSNSmo(InData.me, InTable.me, isuba, err.me)
    if err.isErr:
        printErrMsg(err, "Error smoothing SN table")
    # end  PSNSmo


def PDeselSN (SNTab, isuba, fgid, ants, timerange, err):
    """ Deselect entries in an SN table

    Routine to deselect records in an SN table if they match a given
    subarray, have a selected FQ id, appear on a list of antennas and
    are in a given timerange.
    SNTab     = Python AIPS SN Table
    isuba     = subarray, <=0 -> any
    fqid      = Selected FQ id, <=0 -> any
    ants      = array of integer antenna numbers, 0->all
    timerange = timerange (days)
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inUVGSolve):
        raise TypeError, 'Bad input gain solver'
    if not Table.PIsA(SNTab):
        raise TypeError,"SNTab MUST be a Python Obit Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    nantf = len(ants)
    Obit.UVSolnDeselSN (SNTab.me, isuba, fqid, nantf, ants, timerange,
                        err.me)
    if err.isErr:
        printErrMsg(err, "Error deselecting solutions")
    # end PDeselSN

def PDeselCL (CLTab, isuba, fgid, ants, timerange, err):
    """ Deselect entries in an CL table

    Routine to deselect records in an CL table if they match a given
    subarray, have a selected FQ id, appear on a list of antennas and
    are in a given timerange.
    CLTab     = Python AIPS CL Table
    isuba     = subarray, <=0 -> any
    fqid      = Selected FQ id, <=0 -> any
    ants      = array of integer antenna numbers, 0->all
    timerange = timerange (days)
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inUVGSolve):
        raise TypeError, 'Bad input gain solver'
    if not Table.PIsA(CLTab):
        raise TypeError,"CLTab MUST be a Python Obit Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    nantf = len(ants)
    Obit.UVSolnDeselCL (CLTab.me, isuba, fqid, nantf, ants, timerange,
                        err.me)
    if err.isErr:
        printErrMsg(err, "Error deselecting solutions")
    # end PDeselCL

def PGetName (inUVGSolve):
    """ Tells Image object name (label)

    returns name as character string
    inUVGSolve  = Python UVGSolve object
    """
    ################################################################
     # Checks
    if not PIsA(inUVGSolve):
        raise TypeError,"inUVGSolve MUST be a Python Obit UVGSolve"
    #
    return Obit.UVGSolveGetName(inUVGSolve.me)
    # end PGetName

def PIsA (inUVGSolve):
    """ Tells if input really a Python Obit UVGSolve

    return true, false (1,0)
    inUVGSolve   = Python UVGSolve object
    """
    ################################################################
    # Checks
    if inUVGSolve.__class__ != UVGSolve:
        return 0
    return Obit.UVGSolveIsA(inUVGSolve.me)
    # end PIsA
