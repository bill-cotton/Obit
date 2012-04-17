""" Python Obit interface to ASDM descriptive data

This class is for creating and using the interface to an ASDM
The information meant to derive contents and intent of a BDF data set.
"""
# $Id$
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
#  GNU General Public License f more details.
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

# Python shadow class to partial ObitSDMData class
import Obit

class OASDMPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.OASDMUnref(Obit.OASDM_me_get(self.this))
            # In with the new
            Obit.OASDM_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != OASDM:
            return
        if name == "me" : 
            return Obit.OASDM_me_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        if self.__class__ != OASDM:
            return
        return "<C OASDM instance> " + Obit.OASDMGetName(self.me)

class OASDM(OASDMPtr):
    """
    Python Obit interface to ASDM intent related data
    
    This class is for creating and using the interface to a partial ASDM

    OASDM Members with python interfaces:

    * refJD - reference JD (propably not good)
    * Array - name of array ('EVLA' 'ALMA', recognized)
      ASDM XML tables, see ASDM documentation:
    * Main            - Main table as list of dict
    * Scan            - Scan table as list of dict
    * Config          - configuration table as list of dict
    * CorrelatorMode  - correlator table as list of dict
    * DataDescription - Data description table as list of dict
    * SpectralWindow  - Spectral window table as list of dict
    * Antenna         - Antenna table as list of dict
    * Station         - Station table as list of dict
    * State           - State table as list of dict
    * ExecBlock       - ExecBlock table as list of dict
    * Source          - Source table as list of dict
    * Field           - Field table as list of dict
    * Feed            - Feed table as list of dict
    * Polarization    - Polarization table as list of dict
    * Processor       - Processor table as list of dict
    * SwitchCycle     - SwitchCycle  table as list of dict

    """
    def __init__(self, err, name="ASDM", DataRoot="None"):
        """
        Create ASDM object
        
        DataRoot      = directory path to root of ASDM/BDF
        err           = Obit error/message object
        """
        self.this = Obit.new_OASDM(name, DataRoot, err.me)
    def __del__(self):
        if Obit!=None:
            Obit.delete_OASDM(self.this)
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.OASDM_me_get(self.this)
        # Functions to return members
        if name=="Scan":
            out = Obit.OASDMGetScan(self.me)
            return out
        if name=="Main":
            out = Obit.OASDMGetMain(self.me)
            return out
        if name=="Config":
            out = Obit.OASDMGetConfig(self.me)
            return out
        if name=="CorrelatorMode":
            out = Obit.OASDMGetCorrelatorMode(self.me)
            return out
        if name=="DataDescription":
            out = Obit.OASDMGetDataDescription(self.me)
            return out
        if name=="SpectralWindow":
            out = Obit.OASDMGetSpectralWindow(self.me)
            return out
        if name=="Antenna":
            out = Obit.OASDMGetAntenna(self.me)
            return out
        if name=="Station":
            out = Obit.OASDMGetStation(self.me)
            return out
        if name=="State":
            out = Obit.OASDMGetState(self.me)
            return out
        if name=="ExecBlock":
            out = Obit.OASDMGetExecBlock(self.me)
            return out
        if name=="Source":
            out = Obit.OASDMGetSource(self.me)
            return out
        if name=="Field":
            out = Obit.OASDMGetField(self.me)
            return out
        if name=="Feed":
            out = Obit.OASDMGetFeed(self.me)
            return out
        if name=="Polarization":
            out = Obit.OASDMGetPolarization(self.me)
            return out
        if name=="Processor":
            out = Obit.OASDMGetProcessor(self.me)
            return out
        if name=="SwitchCycle":
            out = Obit.OASDMGetSwitchCycle(self.me)
            return out
        if name=="refJD":
            return Obit.OASDMGetRefJD(self.me)
        if name=="Array":
            return Obit.OASDMGetArray(self.me)
        raise AttributeError,str(name)
    
    def GetSpectralWindowArray(self, mainrow=0, SWOrder=True):
        """ Partially digested Spectral window array """
        out = Obit.OASDMGetSpectralWindowArray(self.me, mainrow, SWOrder)
        return out
    def GetAntArray(self, mainrow=0):
        """ Partially digested Antenna array """
        out = Obit.OASDMGetAntArray(self.me, mainrow)
        return out
    def GetSourceArray(self):
        """ Partially digested Source array """
        out = Obit.OASDMGetSourceArray(self.me)
        return out
    
    # Enumerations
    # Bands
    ASDMBand_Any = 0
    ASDMBand_4 = ASDMBand_Any+1
    ASDMBand_P = ASDMBand_4+1+1
    ASDMBand_L = ASDMBand_P+1
    ASDMBand_S = ASDMBand_L+1
    ASDMBand_C = ASDMBand_S+1
    ASDMBand_X = ASDMBand_C+1
    ASDMBand_Ku = ASDMBand_X+1
    ASDMBand_K =  ASDMBand_Ku+1
    ASDMBand_Ka = ASDMBand_K+1
    ASDMBand_Q =  ASDMBand_Ka+1
    ASDMBand_A3 = ASDMBand_Q+1
    ASDMBand_A4 = ASDMBand_A3+1
    ASDMBand_A5 = ASDMBand_A4+1
    ASDMBand_A6 = ASDMBand_A5+1
    ASDMBand_A7 = ASDMBand_A6+1
    ASDMBand_A8 = ASDMBand_A7+1
    ASDMBand_A9 = ASDMBand_A8+1
    ASDMBand_A10 = ASDMBand_A9+1
    ASDMBand_A11 = ASDMBand_A10+1

    def GetPhaseCal (self, config):
        """
        Return list of phase calibrator names
    
        * self   = ASDM object
        """
        out = []
        scan = self.Scan
        main = self.Main
        for m in main:
            if m["configDescriptionId"]==config:
                for s in scan:
                    want = False
                    # Check intents
                    for int in s['scanIntent']:
                        if int=='CALIBRATE_PHASE' and s['sourceName'] not in out:
                            out.append(s['sourceName'])
        del scan, main
        return out
    # end GetPhaseCal 
    
    def GetAmpCal (self, config):
        """
        Return list of amplitude calibrator names
    
        * self   = ASDM object
         * config = configuration ID
       """
        out = []
        scan = self.Scan
        main = self.Main
        for m in main:
            if m["configDescriptionId"]==config:
                for s in scan:
                    # Check intents
                    for int in s['scanIntent']:
                        if int=='CALIBRATE_AMPLI' and s['sourceName'] not in out:
                            out.append(s['sourceName'])
        del scan, main
        return out
    # end GetPhaseCal 
                    
    def GetBandpassCal (self, config):
        """
        Return list of bandpass calibrator names
    
        * self   = ASDM object
        * config = configuration ID
        """
        out = []
        scan = self.Scan
        main = self.Main
        for m in main:
            if m["configDescriptionId"]==config:
                for s in scan:
                    # Check intents
                    for int in s['scanIntent']:
                        if int=='CALIBRATE_BANDPASS' and s['sourceName'] not in out:
                            out.append(s['sourceName'])
        del scan, main
        return out
    # end GetBandpassCal 
                    
    def GetTargets (self, config):
        """
        Return list of bandpass calibrator names
    
        * self   = ASDM object
        * config = configuration ID
        """
        out = []
        scan = self.Scan
        main = self.Main
        for m in main:
            if m["configDescriptionId"]==config:
                for s in scan:
                    # Check intents
                    for int in s['scanIntent']:
                        if int=='OBSERVE_TARGET' and s['sourceName'] not in out:
                            out.append(s['sourceName'])
        del scan, main
        return out
    # end GetTargets 
                    
    def GetArrayConfig (self):
        """
        Return Array configuration name ("A", "B", "CnD"...)
    
        * self   = ASDM object
        """
        eb = self.ExecBlock
        out = eb[0]["configName"]
        del eb
        return out
    # end GetTargets 
                    
    def GetConfigs (self):
        """
        Return list of configurarion info dicts

        `Each entry contains configDescriptionId, numAntenna, correlationMode,
        spectralType, avgRefFreq, nspWinds a list of the number of different numbers of channels
        and a list (spWinds)of spectral windows with: 
            spectralWindowId, refFreq, chanFreqStep, numChan
        correlationMode values 0 = cross only, 1 = auto only, 2=cross & auto
        spectralType values 0 = channel average, 1= baseband wide, 2=full resolution
        * self   = ASDM object
        """
        out = []
        sumFreq = 0.0; cntFreq = 0
        config = self.Config
        sw     = self.SpectralWindow
        dd     = self.DataDescription
        for c in config:
            dict = {"configDescriptionId":c["configDescriptionId"], \
                    "numAntenna":c["numAntenna"],                   \
                    "correlationMode":c["correlationMode"],         \
                    "spectralType":c["spectralType"]}
            # lookup spectral Windows, screwy design
            swlist = []
            chlist = []
            for dId in c["dataDescriptionId"]:
                for ddata in dd:
                    if ddata["dataDescriptionId"]==dId:
                        for swin in sw:
                            if swin["spectralWindowId"]==ddata["spectralWindowId"]:
                                swdict = {"spectralWidowId":swin["spectralWindowId"], \
                                          "refFreq":swin["refFreq"],                  \
                                          "chanFreqStep":swin["chanFreqStep"],        \
                                          "numChan":swin["numChan"]}
                                # Keep track of channels used
                                if swin["numChan"] not in chlist:
                                    chlist.append(swin["numChan"])
                                sumFreq += swin["refFreq"]; cntFreq += 1
                swlist.append(swdict)
            dict["spWinds"]  = swlist
            dict["nspWinds"] = len(swlist)
            dict["nchands"]  = chlist
            # Average reference frequency
            dict["avgRefFreq"] = sumFreq/cntFreq
            out.append(dict)
        del sw, config, dd
        return out
    # end GetConfigs 
                    
    def GetSourceInfo (self, source):
        """
        Return dict with source info
    
        Entries:
          codes    = list of codes used
          IDs      = list out Source IDs used
          intents  = list of stated intents
          position = [RA, Dec]
        * self   = ASDM object
        * source = source name
        """
        codes = []
        pos   = []
        IDs   = []
        field = self.Field
        for f in field:
            # Check intents
            if f["fieldName"]==source:
                codes.append(f["code"])
                IDs.append(f["sourceId"])
                pos = f["referenceDir"]
        del field
        intents = []
        scan = self.Scan
        for s in scan:
            # Check intents
            for int in s['scanIntent']:
                if int not in intents:
                    intents.append(int)
        del scan
        return {"codes":codes, "IDs":IDs, "intents":intents, "position":pos}
    # end GetSourceInfo 
                    
    def Get1stBandpassScan (self, config):
        """
        Return info on first bandpass calibrator scan
    
        Return dict with entries:
          source    = source name
          timerange = timerange in days wrt first scan
        * self   = ASDM object
        * config = configuration ID
        """
        scan = self.Scan
        main = self.Main
        refJD = int(scan[0]['startTime'])
        source = None; timeRange = [0.0, 1000.]
        for m in main:
            if m["configDescriptionId"]==config:
                s = scan[m["scanNumber"]-1]
                # Check intents
                for intn in s['scanIntent']:
                    if intn=='CALIBRATE_BANDPASS':
                        source = s['sourceName']
                        timeRange = [s['startTime']-refJD, s['endTime']-refJD]
                        break;
                if source!=None:
                    break;
        del scan,main
        return {"source":source,"timeRange":timeRange}
    # end Get1stBandpassScan

    # end class OASDM
    
def Freq2Band (freq):
    """ Determine Band corresponding to a given frequency (only approximate)

    Returns one of 
      ASDMBand_Any,  ASDMBand_4,  ASDMBand_P,  ASDMBand_L,  ASDMBand_S,
      ASDMBand_C,  ASDMBand_X,  ASDMBand_Ku,  ASDMBand_K,  ASDMBand_Ka,
      ASDMBand_Q,  ASDMBand_A3,  ASDMBand_A4,  ASDMBand_A5,  ASDMBand_A6,
      ASDMBand_A7,  ASDMBand_A8,  ASDMBand_A9,  ASDMBand_A10,  ASDMBand_A11
    * freq    = frequency (Hz)
    """
    if (freq<100.0e6):
        return ASDMBand_4;
    if (freq<900.0e6):
        return ASDMBand_P;
    if (freq<2.0e9):
        return ASDMBand_L;
    if (freq<3.7e9):
        return ASDMBand_S;
    if (freq<7.5e9):
        return ASDMBand_C;
    if (freq<12.0e9):
        return ASDMBand_X;
    if (freq<18.0e9):
        return ASDMBand_Ku;
    if (freq<26.5e9):
        return ASDMBand_K;
    if (freq<40.0e9):
        return ASDMBand_Ka;
    if (freq<50.0e9):
        return ASDMBand_Q;
    if (freq<117.0e9):
        return ASDMBand_A3;
    if (freq<163.0e9):
        return ASDMBand_A4;
    if (freq<211.0e9):
        return ASDMBand_A5;
    if (freq<275.0e9):
        return ASDMBand_A6;
    if (freq<375.0e9):
        return ASDMBand_A7;
    if (freq<510.0e9):
        return ASDMBand_A8;
    if (freq<730.0e9):
        return ASDMBand_A9;
    if (freq<960.0e9):
        return ASDMBand_A10;
    if (freq<2000.0e9):
        return ASDMBand_A11;
    return out;
    # end Freq2Band
    
def PSWChanSel(spWin, selChan, selIF, band=OASDM.ASDMBand_Any):
    """ Select Spectral windows by number of channels/band

    Returns True if any windows are selected, else False
    * spWin    = Spectral window array returned by GetSpectralWindowArray
    * selChan  = Selected no. channels
    * selIF    = Selected no. IFs
    * band     = band code, one of
      ASDMBand_Any,  ASDMBand_4,  ASDMBand_P,  ASDMBand_L,  ASDMBand_S,
      ASDMBand_C,  ASDMBand_X,  ASDMBand_Ku,  ASDMBand_K,  ASDMBand_Ka,
      ASDMBand_Q,  ASDMBand_A3,  ASDMBand_A4,  ASDMBand_A5,  ASDMBand_A6,
      ASDMBand_A7,  ASDMBand_A8,  ASDMBand_A9,  ASDMBand_A10,  ASDMBand_A11
    """
    out = False  # Until shown otherwise

    nwinds = spWin["nwinds"]
    for iSW in range(0,nwinds):
        wind = spWin["winds"][iSW]
        OK = (wind["numChan"]==selChan) and \
             ((Freq2Band (wind["refFreq"])== band) or (band==ASDMBand_Any)) and \
             (nwinds==selIF)
        wind["selected"] = wind["selected"] or OK
        if wind["selected"]:
            out = TRUE
    return out
# end  PSWChanSel

def PIsA (asdm):
    """
    Tells if the input is a Python OASDM
    
    returns True or False (1,0)
    * asdm = Python Obit ASDM to test
    """
    ################################################################
      # Checks
    if printer.__class__ != OASDM:
        return False
    return Obit.OASDMIsA(asdm.me)
    # end PIsA


