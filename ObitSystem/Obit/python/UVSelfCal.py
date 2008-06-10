# $Id: UVSelfCal.py,v 1.11 2007/01/31 15:27:35 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C)2005,2006,2007
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
import Obit, OErr, InfoList, UV, Table, types

# Python shadow class to ObitUVSelfCal class
 
class UVSelfCalPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.UVSelfCalUnref(Obit.UVSelfCal_me_get(self.this))
            # In with the new
            Obit.UVSelfCal_me_set(self.this,value)
            return
        if name=="SkyModel":
            PSetSkyModel(self, value)
            return 
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != UVSelfCal:
            return
        if name == "me" : 
            return Obit.UVSelfCal_me_get(self.this)
        # Virtual members
        if name=="List":
            return PGetList(self)
        if name=="SkyModel":
            return PGetSkyModel(self)
        raise AttributeError,str(name)
    def __repr__(self):
        if self.__class__ != UVSelfCal:
            return
        return "<C UVSelfCal instance> " + Obit.UVSelfCalGetName(self.me)
#
class UVSelfCal(UVSelfCalPtr):
    """ Python Obit Image class

    This class does visibility-based (Cotton-Schwab) CLEAN

    UVSelfCal Members with python interfaces:
    InfoList  - used to pass instructions to processing
                Member List
                (readonly)
    skyModel    - SkyModel, use PGetSkyModel
    """
    def __init__(self, name) :
        self.this = Obit.new_UVSelfCal(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_UVSelfCal(self.this)


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
       Soln2CalInput={'structure':['Soln2Cal',[('InData','Input UV'),
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
    """ Create and initialize an UVSelfCal structure

    Create sky model object
    Returns the Python UVSelfCal object
    name     = name desired for object (labeling purposes)
    err      = Python Obit Error/message stack
    """
    ################################################################
    out = UVSelfCal (name)
    return out      # seems OK
    # end newObit

def PCopy (inUVSelfCal, outUVSelfCal, err):
    """ Make a shallow copy of input object.

    Makes structure the same as inUVSelfCal, copies pointers
    inUVSelfCal  = Python UVSelfCal object to copy
    outUVSelfCal = Output Python UVSelfCal object, must be defined
    err         = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inUVSelfCal):
        raise TypeError,"inUVSelfCal MUST be a Python Obit UVSelfCal"
    if not PIsA(outUVSelfCal):
        raise TypeError,"outUVSelfCal MUST be a Python Obit UVSelfCal"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.UVSelfCalCopy (inUVSelfCal.me, outUVSelfCal.me, err.me)
    if err.isErr:
        printErrMsg(err, "Error copying UVSelfCal")
    # end PCopy

def PGetList (inUVSelfCal):
    """ Return the member InfoList

    returns InfoList
    inUVSelfCal  = Python UVSelfCal object
    """
    ################################################################
     # Checks
    if not PIsA(inUVSelfCal):
        raise TypeError,"inUVSelfCal MUST be a Python Obit UVSelfCal"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.UVSelfCalGetList(inUVSelfCal.me)
    return out
    # end PGetList

def PGetSkyModel (inUVSelfCal):
    """ Return the member sky model

    returns ImageMosaic
    inUVSelfCal  = Python UVSelfCal object
    """
    ################################################################
     # Checks
    if not PIsA(inUVSelfCal):
        raise TypeError,"inUVSelfCal MUST be a Python Obit UVSelfCal"
    #
    out    = SkyModel.SkyModel("None")
    out.me = Obit.UVSelfCalGetSkyModel(inUVSelfCal.me)
    return out
    # end PGetSkyModel

def PSetSkyModel (inUVSelfCal, skyModel):
    """ Replace the SkyModel in the UVSelfCal

    inUVSelfCal = Python UVSelfCal object
    skyModel    = Python SkyModel to attach
    """
    ################################################################
    # Checks
    if not PIsA(inUVSelfCal):
        raise TypeError,"inUVSelfCal MUST be a Python ObitUVSelfCal"
    if not SkyModel.PIsA(skyModel):
        raise TypeError,"array MUST be a Python Obit SkyModel"
    #
    Obit.UVSelfCalSetSkyModel(inUVSelfCal.me, mosaic.me)
    # end PSetSkyModel

# Define UVSelfCal - most parameters can be defined here
SelfCalInput={'structure':['SelfCal',[('subA ','Selected subarray (default 1)'),
                                      ('solInt','Solution interval (min). (default 1 sec)'),
                                      ('refAnt','Ref ant to use. (default 1)'),
                                      ('avgPol','True if RR and LL to be averaged (false)'),
                                      ('avgIF','True if all IFs to be averaged (false)'),
                                      ('minSNR','Minimum acceptable SNR (5)'),
                                      ('doMGM','True then find the mean gain modulus (false)'),
                                      ('solType','Solution type " ", "L1",  (" ")'),
                                      ('solMode','Solution mode: "A&P", "P", "P!A", ("P")'),
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

def PCreate (name, skyModel, err, input=SelfCalInput):
    """ Create the parameters and underlying structures of a UVSelfCal.

    Returns UVSelfCal created.
    name      = Name to be given to object
                Most control parameters are in InfoList member
    skyModel    = Python uv data from which image is to be made
    err       = Python Obit Error/message stack
    input     = control parameters:
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
    antWt    = Weight per antenna, def all 1
    prtLv       = Print level (default no print)
    """
    ################################################################
    # Checks
    if not SkyModel.PIsA(skyModel):
        raise TypeError,"skyModel MUST be a Python Obit SkyModel"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Create
    out = UVSelfCal(name);
    out.me = Obit.UVSelfCalCreate(name, skyModel.me)
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
    #OErr.printErrMsg(err, "UVSelfCalCreate: Error setting parameters")
    #
    return out;
    # end PCreate

# Self calibrate
UVSelfCalCalInput={
    'structure':['UVSelfCalCal',[('InData',  'Input UV data'),
                                 ('OutData', 'Output UV data'),
                                 ('subA ','Selected subarray (default 1)'),
                                 ('solInt','Solution interval (min). (default 1 sec)'),
                                 ('refAnt','Ref ant to use. (default 1)'),
                                 ('avgPol','True if RR and LL to be averaged (false)'),
                                 ('avgIF','True if all IFs to be averaged (false)'),
                                 ('minSNR','Minimum acceptable SNR (5)'),
                                 ('doMGM','True then find the mean gain modulus (false)'),
                                 ('solType','Solution type "  ", "L1",  (" ")'),
                                 ('solMode','Solution mode: "A&P", "P", "P!A", ("P")'),
                                 ('minNo','Min. no. antennas. (default 4)'),
                                 ('WtUV','Weight outside of UVRANG. (default 0.0)'),
                                 ('antWt','Weight per antenna, def all 1'),
                                 ('prtLv','Print level (default no print)')]],
    # defaults
    'InData':None,
    'OutData':None,
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
def PCal (inSC, err, input=UVSelfCalCalInput):
    """ Self calibrate data

    Determine a gain table from a set of data assumed divided by model.
    Solution table attached to  OutData
    inSC    = Selfcal object
    err     = Python Obit Error/message stack
    input   = input parameter dictionary
    
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
    if not PIsA(inSC):
        raise TypeError, 'Bad input selfcalibrator'
    if not UV.PIsA(InData):
        raise TypeError, 'Bad input UV data'
    if not UV.PIsA(OutData):
        raise TypeError, 'Bad output UV data'
    dim = [1,1,1,1,1]
    # Create solver
    solver = UVGSolve.newObit("Selfcal Solver", err)
    # Set control values on solver
    dim[0] = 1;
    inInfo = PGetList(solver)    # 
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
    OErr.printErrMsg(err, "UVSelfCalCreate: Error setting parameters")
    #
    # Create output
    outTable    = Table.Table("None")
    # Solve
    outTable.me = UVGSolveCal(solver.me, InData.me, OutData.me, err)
    #outTable.me = UVSelfCalCal(inSC.me, InData.me, OutData.me, err)
    # show any errors 
    OErr.printErrMsg(err, "UVSelfCalCal: Error Self-calibrating")
    #
    return outTable
    # end  PCal

def PRef (inSC, SNTab, isuba, refant, err):
    """ rereference phases of subarray isuba in a SN table to refant

    inSC      = Selfcal object
    SNTab     = Python AIPS SN Table
    isuba     = subarray [def. 1]
    refant    = reference antenna, if 0 -> pick one
    err       = Python Obit Error/message stack
    Returns actual reference antenna used
    """
    ################################################################
    # Checks
    if not PIsA(inSC):
        raise TypeError, 'Bad input selfcalibrator'
    if not Table.PIsA(SNTab):
        raise TypeError,"SNTab MUST be a Python Obit Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.UVSolnRefAnt (SNTab.me, isuba, refant, err.me)
    if err.isErr:
        printErrMsg(err, "Error rereferencing phases")
    return ret
    # end PRef


# Self calibrate
UVSelfSNSmoInput={
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
def PSNSmo (inSC, err, input=UVSelfSNSmoInput):
    """ Smooth SN table possibly replacing blanked soln.

    inSC    = Selfcal object
    err     = Python Obit Error/message stack
    input   = input parameter dictionary
    
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
    Obit.UVSolnSNSmo(InTable.me, isuba, err.me)
    if err.isErr:
        printErrMsg(err, "Error smoothing SN table")
    # show any errors 
    #OErr.printErrMsg(err, "UVSelfCalCal: Error Smoothing SN table")
    # end  PSNSmo


def PDeselSN (inSC, SNTab, isuba, fgid, ants, timerange, err):
    """ Deselect entries in an SN table

    Routine to deselect records in an SN table if they match a given
    subarray, have a selected FQ id, appear on a list of antennas and
    are in a given timerange.
    inSC      = Selfcal object
    SNTab     = Python AIPS SN Table
    isuba     = subarray, <=0 -> any
    fqid      = Selected FQ id, <=0 -> any
    ants      = array of integer antenna numbers, 0->all
    timerange = timerange (days)
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inSC):
        raise TypeError, 'Bad input selfcalibrator'
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

def PInvertSN (SNTab, outUV, outVer, err):
    """ Invert the calibration in an SN table

    Routine to reverse the effects of the calibration in the
    input SN table and create a new SN table on outUV
    Returns new SN table
    SNTab     = Input Python AIPS SN Table to invert
    outUV     = output UV data to which to attach the new SN table
    outVer    = Output SN table version number, 0=> create new
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(SNTab):
        raise TypeError,"SNTab MUST be a Python Obit Table"
    if not UV.PIsA(outUV):
        raise TypeError, 'outUV Must be UV data '
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Create output
    outSN  = Table.Table("None")
    # Cast uv to Data
    outData = Obit.UVCastData(outUV.me)
    # invert table
    outSN.me = Obit.SNInvert (SNTab.me, outData, outVer, err.me)
    if err.isErr:
        printErrMsg(err, "Error inverting solutions")
    return outSN
    # end PInvertSN

def PGetName (inUVSelfCal):
    """ Tells Image object name (label)

    returns name as character string
    inUVSelfCal  = Python UVSelfCal object
    """
    ################################################################
     # Checks
    if not PIsA(inUVSelfCal):
        raise TypeError,"inUVSelfCal MUST be a Python Obit UVSelfCal"
    #
    return Obit.UVSelfCalGetName(inUVSelfCal.me)
    # end PGetName

def PIsA (inUVSelfCal):
    """ Tells if input really a Python Obit UVSelfCal

    return true, false (1,0)
    inUVSelfCal   = Python UVSelfCal object
    """
    ################################################################
    # Checks
    if inUVSelfCal.__class__ != UVSelfCal:
        return 0
    return Obit.UVSelfCalIsA(inUVSelfCal.me)
    # end PIsA
