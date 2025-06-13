# Template MeerKAT project parameter file 
# Generate parameter file using MKPrepare
#
# Substitutions surrounded by 'at' characters
# PROJECT     Project name (up to 12 char)
# SESSION     Session code
# BAND        Band code
# DATAFILE    Full or relative path to Raw data uvtab file
# DCALFILE    Full or relative path to DelayCal Raw data uvtab file
# DOPOL       True if polarization cal. and imaging wanted (else False)
# BPCAL       Bandpass calibrator list
# GAINCAL     Phase calibrator list
# AMPCAL      Amplitude calibrator list
# DLYCAL      Delay calibrator list
# UNPOLCAL    Instrumental polarization calibrator list
# POLCAL      Known polarized calibrator list
# REFANT      Reference antenna
# PLOTSRC     Diagnostic plot source name or None
# PLOTTIIME   Diagnostic plot timerange
# TARGETS     List of target sources to image
# STOKES      Stokes parameters to image

# DESTDIR     Output directory 
#--------------------------------------------------------------------------------------------
# Project specific parameter values for EVLAPipeline
project       = "@PROJECT@"           # Project name (12 char or less, used as AIPS Name)
session       = "@SESSION@"           # Project session code
band          = "@BAND@"              # Observing band
dataClass     = band+"Band"           # AIPS class of raw uv data
logFile       = project+"_"+session+"_"+band+".log"  # Processing log file

# Archive data parameters
parms["doLoadArchive"] = True       # Load from archive uvtab?
parms["DataFile"]      = '@DATAFILE@' # Name of Raw Data uvtab file
parms["DCalFile"]      = '@DCALFILE@' # Name of DelayCal Raw Data uvtab file
# Get metadata from data
inUV = UV.newPFUV('Raw',parms["DataFile"], 0, True, err)
meta = MKGetMeta(inUV, {}, "", err)
parms['MKFreq'] = meta["MKFreq"]  

# Calibration sources/models
parms["BPCal"]       = @BPCAL@      # Bandpass calibrator


from MeerKATCal import MKCalModel,MKStdModel
# Amp/phase calibration
calist = @GAINCAL@
# List of Gain calibrators
ACals = []
GCalList = []
for cal in calist:
    ACals.append(MKCalModel(cal))
    GCalList.append(cal)
parms['GCalList'] = GCalList
# Sources for phase calibration
PCals = []
for cal in calist:
    PCals.append(MKCalModel(cal))
# Also Amp cals on phase cal list
calist = @AMPCAL@
for cal in calist:
    PCals.append(MKCalModel(cal))

# Amplitude calibration
calist = @AMPCAL@
for cal in calist:
    ACals.append(MKCalModel(cal))
# Check for standard model
ACals = MKStdModel(ACals, parms['MKFreq']*1.0e-6)
parms["ACals"] = ACals  # Amplitude calibrators

# Bandpass Calibration
calist = @BPCAL@
BPCals = []
for cal in calist:
    BPCals.append(MKCalModel(cal))
    PCals.append(MKCalModel(cal))  # Add to phase calibration
# Check for standard model
BPCals = MKStdModel(BPCals, parms['MKFreq']*1.0e-6)
parms["BPCals"]  = BPCals      # Bandpass calibrator(s)

# Delay calibration
calist = @DLYCAL@
DCals = []
for cal in calist:
    DCals.append(MKCalModel(cal))
    PCals.append(MKCalModel(cal))  # Add to phase calibration
# Check for standard model
DCals = MKStdModel(DCals, parms['MKFreq']*1.0e-6)
parms["DCals"] = DCals      # delay calibrators

parms["refAnt"]        = @REFANT@   # Reference antenna

# Sample spectra
parms["plotSource"]    = @PLOTSRC@        # Source name or None
parms["plotTime"]      = @PLOTTIME@       # timerange
parms["doRawSpecPlot"] = @PLOTSRC@!=None  # Plot Raw spectrum
parms["doSpecPlot"]    = @PLOTSRC@!=None  # Plot spectrum at various stages of processing

# Poln  Cal
parms["doPol"]         = @DOPOL@      # Do polarization calibration?
parms["UnPolCal"]      = @UNPOLCAL@   # List of instrumental poln calibrators
parms["doPolCal"]      = len(parms["UnPolCal"])>0   # Do polarization (PCal) calibration?
parms["doPolSpecPlot"] = parms["doPol"] and @POLCAL@!=None  # Plot polarized spectrum?


# X-Y phase/delay calibration
parms["XYDCal"]    = @POLCAL@     # X-Y delay calibrator list
parms["xyrefAnt"]  = @REFANT@     # Reference antenna for X-Y cal, defaults to refAnt
parms["doXYDelay"] = parms["doPol"] and (len(parms["XYDCal"])>0)  # Determine X-Y delay? If calibrator given
parms["doPhsCal"]  = parms["doPol"] and ((len(parms["XYDCal"])>0)  or (len(parms["UnPolCal"])>0))
XYCals = []
if parms["doXYDelay"]:
    for cal in parms["XYDCal"]:
        XYCals.append(MKPolModel(cal))
    parms["XYCals"]   = XYCals         # XY delay calibrator(s)
else:
     parms["XYCals"]   = None

# Check for standard model for phase calibration
PCals = MKStdModel(PCals, parms["MKFreq"]*1.0e-6)
parms["PCals"] = PCals   # Phase calibrator(s)

# Imaging
parms["targets"] = @TARGETS@     # targets, empty = all
parms["Stokes"]  = "@STOKES@"    # Stokes to image
parms["doMB"] = True # MeerKAT always wideband

################## The following might need fiddling #######################
parms["doFD1"]       = True         # Do initial frequency domain flagging
parms["FD1widMW"]    = 55           # Width of the initial FD median window
parms["FD1maxRes"]   = 10.0         # Clipping level in sigma
parms["FD1TimeAvg"]  = 2.0          # time averaging in min. for initial FD flagging
parms["FD1baseSel"]   = [0,0,0,0]   # Channels for baseline fit (start, end, increment, IF)

parms["doMedn"]      = True         # Median editing?
parms["mednSigma"]   = 10.0         # Median sigma clipping level
parms["timeWind"]    = 2.0          # Median window width in min for median flagging
parms["avgTime"]     = 10.0/60.     # Averaging time in min
parms["avgFreq"]     = 1            # 1=>avg chAvg chans, 2=>avg all chan, 3=> avg chan and IFs
parms["chAvg"]       = 2            # number of channels to average

parms["doRMSAvg"]    = True         # Edit calibrators by RMSAvg?
parms["RMSAvg"]      = 5.0          # AutoFlag Max RMS/Avg for time domain RMS filtering
parms["RMSTimeAvg"]  = 1.0          # AutoFlag time averaging in min.

# Special editing list
doEditList  = False        # Edit using editList?
editFG      = 2            # Table to apply edit list to
# Channel numbers after Hanning if any
editList = [
    #{"timer":("0/00:00:0.0","5/00:00:0.0"),"Ant":[ 1,0],"IFs":[1,0],"Chans":[1,0],  "Stokes":'1111',"Reason":"No Rcvr"},
    ]

################## The following flags control the script executation #######################
# Control, mark items as F to disable
T   = True
F   = False
parms["nThreads"]      = 16       # number of threads to allow, overridden AIPSSetup
check                  = F        # Only check script, don't execute tasks
debug                  = F        # run tasks debug
parms["doHann"]        = parms["doHann"] # Load data/Hanning smooth?
parms["doClearTab"]    = T        # Clear cal/edit tables
parms["doEditList"]    = T        # Special editing
parms["doCopyFG"]      = T        # Copy FG 1 to FG 2
parms["doShadow"]      = parms["doShadow"] # Flag shadowed data?
parms["doMedn"]        = T        # Median editing?
parms["doFD1"]         = T        # Do initial frequency domain flagging
parms["doRMSAvg"]      = T        # Do RMS/Mean editing for calibrators MAY NEED THIS
parms["doRawSpecPlot"] = parms["doRawSpecPlot"]  # Plot sample raw spectra?
parms["doNDCal"]       = T        # Noise Diode calibration?  Only for doPol
parms["doDelayCal"]    = T        # Group Delay calibration?
parms["doBPCal"]       = T        # Determine Bandpass calibration
parms["doAmpPhaseCal"] = T        # Amplitude/phase calibration
parms["doAutoFlag"]    = T        # Autoflag editing after final calibration?
parms["doClipCals"]    = T        # Autoflag Clipping on Calibrators
parms["doRecal"]       = T        # Redo calibration after editing
parms["doNDCal2"]      = T        # 2nd  Noise Diode calibration?  Only for doPol
parms["doDelayCal2"]   = T        # Group Delay calibration of averaged data?, 2nd pass
parms["doBPCal2"]      = T        # Determine Bandpass calibration, 2nd pass
parms["doAmpPhaseCal2"]= T        # Amplitude/phase calibration, 2nd pass
parms["doAutoFlag2"]   = T        # Autoflag editing after final calibration?
parms["doCalAvg"]      = "BL"     # Calibrate and baseline dependent average data
                                  # "BL"=> bl dependent, "Splat"=> no time averaging.
parms["doPhsCal"]      = parms["doPhsCal"]  # Phase calibrate poln calibrators?
parms["doPolCal"]      = parms["doPolCal"]  # Do instr. polarization calibration?
parms["doXYDelay"]     = parms["doXYDelay"] # Determine X-Y delay?
parms["doSaveTab"]     = T        # Save UV tables to FITS
parms["doSaveUV"]      = T        # Save calibrated UV data to FITS
parms["doImage"]       = T        # Image targets
parms["doSaveImg"]     = T        # Save results to FITS
parms["doCleanup"]     = T        # Destroy AIPS files

# diagnostics
parms["doSNPlot"]      = T        # Plot SN tables etc
parms["doPolSpecPlot"] = parms["doPolSpecPlot"]  # Plot sample Polarization spectra?
parms["doSpecPlot"]    = parms["doSpecPlot"]     # Plot sample calibrated/edited spectra?
parms["doBPPlot"]      = parms["doBPPlot"]       # Plot XY phase cal BP table?
parms["doReport"]      = T        # Individual source report

