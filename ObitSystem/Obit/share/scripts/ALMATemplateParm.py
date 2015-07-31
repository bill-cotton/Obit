# Template project parameter file for ALMA projects
# Generate parameter file using ALMACal.ALMAMakeParmFile
#
# Substitutions surrounded by 'at' characters
# PROJECT     Project name (up to 12 char)
# SESSION     Session code
# BAND        Band code
# ALMAFREQ    Frequency in Hz
# SPANBW      Spanned Frequency in Hz
# MAXBL       ALMA maximum baseline (km)
# DATAROOT    Root of archive data directory
# CONFIG      SDM Configuration
# SELCHAN     Number of channels for selected configuration
# BPCAL       Bandpass calibrator list
# PHSCAL      Phase calibrator list
# AMPCAL      Amplitude calibrator list
# DLYCAL      Delay calibrator list
# PINCAL      Instrumental polarization calibrator list
# REFANT      Reference antenna
# PLOTSRC     Diagnostic plot source name or None
# PLOTTIIME   Diagnostic plot timerange
# XYGSRC      Set X/Y Gain source name or None
# XYGTIIME    Set X/Y Gain timetange
# XYDSRC      Set X/Y Delay source name or None
# XYDTIIME    Set X/Y Delay timetange
# TARGET      List of target sources
#--------------------------------------------------------------------------------------------
# Project specific parameter values for ALMAPipe
parms["seq"]           = 1                     # Sequence number for AIPS files
parms["project"]       = "@PROJECT@"           # Project name (12 char or less, used as AIPS Name)
parms["session"]       = "@SESSION@"           # Project session code
parms["band"]          = "@BAND@"              # Observing band
parms["dataClass"]     = "@BAND@Band"          # AIPS class of raw uv data
parms["ALMAFreq"]      = @ALMAFREQ@            # Representive frequency
parms["ALMAMaxBl"]     = @MAXBL@               # ALMA maximum baseline (km)

# Archive parameters
parms["doLoadArchive"] = True         # Load from archive?
parms["archRoots"]     = @DATAROOTS@  # List of Roots of ASDM/BDF data
parms["selBand"]       = "@BAND@"     # Selected band, def = first  
parms["selConfig"]     = @CONFIG@     # Selected frequency config, def = first  
parms["selNIF"]        = 0            # Selected number of IFs, def = first  
parms["selChan"]       = @SELCHAN@    # Selected number of channels, def = first  

parms["doFD1"]       = True         # Do initial frequency domain flagging
parms["FD1widMW"]    = 31           # Width of the initial FD median window
parms["FD1maxRes"]   = 5.0          # Clipping level in sigma
parms["FD1TimeAvg"]  = 1.0          # time averaging in min. for initial FD flagging

parms["doMedn"]      = True         # Median editing?
parms["mednSigma"]   = 5.0          # Median sigma clipping level
parms["timeWind"]    = 1.0          # Median window width in min for median flagging
parms["avgTime"]     = 10.0/60.     # Averaging time in min
parms["avgFreq"]     = 0            # 1=>avg chAvg chans, 2=>avg all chan, 3=> avg chan and IFs
parms["chAvg"]       = 1            # number of channels to average

parms["doRMSAvg"]    = True         # Edit calibrators by RMSAvg?
parms["RMSAvg"]      = 3.0          # AutoFlag Max RMS/Avg for time domain RMS filtering
parms["RMSTimeAvg"]  = 1.0          # AutoFlag time averaging in min.

# Special editing list
parms["doEditList"]  = False        # Edit using editList?
parms["editFG"]      = 2            # Table to apply edit list to
# Channel numbers after Hanning if any
parms["editList"] = [
    #{"timer":("0/00:00:0.0","5/00:00:0.0"),"Ant":[ 1,0],"IFs":[1,0],"Chans":[1,0],  "Stokes":'1111',"Reason":"No Rcvr"},
    ]

from ALMACal import ALMACalModel,ALMAStdModel
freq = parms["ALMAFreq"]
# Bandpass Calibration
calist = @BPCAL@
BPCals = []
for cal in calist:
    BPCals.append(ALMACalModel(cal))
# Check for standard model
ALMAStdModel(BPCals, freq)
parms["BPCals"]     = BPCals  # Bandpass calibrator(s)
parms["bpChWid2"]   = 5       # Number of channels to smooth in BP calibration                                                                              

# Phase only and secondary calibrators
calist = @PHSCAL@
PCals = []
tcals = []
for cal in calist:
    if not cal in tcals:
        PCals.append(ALMACalModel(cal))
        tcals.append(cal)
# Check for standard model
ALMAStdModel(PCals, freq)
parms["PCals"]          = PCals   # Phase calibrator(s)

#  Amp/phase calibration = secondary amplitude calibrators
calist = @PHSCAL@
APCals = []
tcals = []
for cal in calist:
    if not cal in tcals:
        APCals.append(ALMACalModel(cal))
        tcals.append(cal)
# Check for standard model
ALMAStdModel(APCals, freq)
parms["APCals"]          = APCals   # calibrator(s)

calist = @AMPCAL@
ACals = []
tcals = []
for cal in calist:
    if not cal in tcals:
        ACals.append(ALMACalModel(cal))
        tcals.append(cal)
# Check for standard model
ALMAStdModel(ACals, freq)
# Set model calibrator flux density
ACals[0]['CalModelFlux'] = 1.0   # ******************* Set this **********************
parms["ACals"]          = ACals  # Amplitude calibrators

calist = @DLYCAL@
DCals = []
tcals = []
for cal in calist:
    if not cal in tcals:
        DCals.append(ALMACalModel(cal))
        tcals.append(cal)
# Check for standard model
ALMAStdModel(DCals, freq)
parms["DCals"]          = DCals      # delay calibrators

parms["refAnt"]        = @REFANT@   # Reference antenna

# Sample spectra
parms["plotSource"]    = @PLOTSRC@          # Source name or None
parms["plotTime"]      = @PLOTTIME@         # timerange

# X/Y Gain calibration
parms["XYGainSource"]    = @XYGSRC@          # Source name or None
parms["XYGainTime"]      = @XYGTIME@         # timerange

# X/Y Delay calibration
parms["XYDelaySource"]    = @XYDSRC@          # Source name or None
parms["XYDelayTime"]      = @XYDTIME@         # timerange

# Instrumental Poln  Cal
PClist                 = @PINCAL@  # List of instrumental poln calibrators
parms["PCInsCals"]     = []
# Remove redundancies 
tcals = []
for cal in PClist:
    if not cal in tcals:
        parms["PCInsCals"].append(cal)
        tcals.append(cal)
parms["doPolCal"]      = len(parms["PCInsCals"])>0  # Do polarization calibration?
parms["doPol"]         = parms["doPolCal"]

# Imaging
parms["targets"] = @TARGET@     # targets, empty = all
parms["Stokes"]  = "I"          # Stokes to image
# Multi frequency or narrow band?
SpanBW = @SPANBW@
if SpanBW<=parms['ALMAFreq']*parms["MBmaxFBW"]*2.0:
    parms["doMB"] = False
parms["MBmaxFBW"] = -parms["MBmaxFBW"]  # IF bins
# May need setting:
parms["IClip"]      = None   # AutoFlag Stokes I clipping, None=>default
parms["CAchAvg"]    = 1      # Channels to average
parms["CalAvgTime"] = None   # Time (min) for averaging calibrated uv data

# Control, mark items as F to disable
T   = True
F   = False
check                  = parms["check"]      # Only check script, don't execute tasks
debug                  = parms["debug"]      # run tasks debug
parms["doLoadArchive"] = T                   # Load from archive?
parms["doHann"]        = parms["doHann"]     # Apply Hanning?
parms["doClearTab"]    = T                   # Clear cal/edit tables
parms["doCopyFG"]      = T                   # Copy FG 1 to FG 2
parms["doEditList"]    = parms["doEditList"] # Edit using editList?
parms["doQuack"]       = T                   # Quack data?
parms["doShad"]        = parms["doShad"]     # Flag shadowed data?
parms["doOnlineCal"]   = T                   # Apply online calibration
parms["doMedn"]        = T                   # Median editing?
parms["doFD1"]         = False               # Do initial frequency domain flagging, Need/dangerous for ALMA?
parms["doRMSAvg"]      = T                   # Do RMS/Mean editing for calibrators
parms["doDelayCal"]    = T                   # Group Delay calibration?
parms["doBPCal"]       = T                   # Determine Bandpass calibration
parms["doXYFixGain"]   = T                   # set X/Y gains and initial calibration
parms["doImgCal"]      = T                   # Self calibrate calibrators
parms["doPhaseCal"]    = T                   # Phase calibration
parms["doAmpPhaseCal"] = T                   # Amplitude/phase calibration
parms["doAutoFlag"]    = T                   # Autoflag editing after final calibration?
parms["doCalAvg"]      = T                   # calibrate and average data
parms["doXYDelay"]     = T                   # X/Y Delay calibration
parms["doPolCal"]      = parms["doPolCal"]   # Instrumental polarization calibration
parms["doImage"]       = T                   # Image targets
parms["doSaveImg"]     = T                   # Save results to FITS
parms["doSaveUV"]      = T                   # Save calibrated UV data to FITS
parms["doSaveTab"]     = T                   # Save UV tables to FITS
parms["doKntrPlots"]   = T                   # Contour plots
parms["doDiagPlots"]   = T                   # Source diagnostic plots
parms["doMetadata"]    = T                   # Generate metadata dictionaries
parms["doHTML"]        = T                   # Generate HTML report
parms["doVOTable"]     = T                   # VOTable report
parms["doCleanup"]     = T                   # Destroy AIPS files

# diagnostics
parms["doSNPlot"]      = T                   # Plot SN tables etc
parms["doBPPlot"]      = T                   # Plot BP tables etc
parms["doReport"]      = T                   # Individual source report
parms["doRawSpecPlot"] = @PLOTSRC@!='None'   # Plot Raw spectrum
parms["doSpecPlot"]    = @PLOTSRC@!='None'   # Plot spectrum at various stages of processing
                                 
