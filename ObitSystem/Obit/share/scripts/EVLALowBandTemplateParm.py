# Template project parameter file for EVLA LowBand (P band) projects
# Generate parameter file using EVLACal.EVLALowBandMakeParmFile
#
# Substitutions surrounded by 'at' characters
# PROJECT     Project name (up to 12 char)
# SESSION     Session code
# BAND        Band code
# VLAFREQ     Frequency in Hz
# SPANBW      Spanned Frequency in Hz
# VLACFG      VLA configuraton ("A", "B", "CnD"...)
# DATAROOT    Root of archive data directory
# CONFIG      Configuration
# SELCHAN     Number of channels for selected configuration
# BPCAL       Bandpass calibrator list
# PHSCAL      Phase calibrator list
# AMPCAL      Amplitude calibrator list
# DLYCAL      Delay calibrator list
# PINCAL      Instrumental polarization calibrator list
# REFANT      Reference antenna
# PLOTSRC     Diagnostic plot source name or None
# PLOTTIIME   Diagnostic plot timerange
# TARGET      List of target sources
#--------------------------------------------------------------------------------------------
# Project specific parameter values for EVLALowBandPipe
parms["seq"]           = 1                     # Sequence number for AIPS files
parms["project"]       = "@PROJECT@"           # Project name (12 char or less, used as AIPS Name)
parms["session"]       = "@SESSION@"           # Project session code
parms["band"]          = "@BAND@"              # Observing band
parms["dataClass"]     = "@BAND@Band"          # AIPS class of raw uv data
parms["VLAFreq"]       = @VLAFREQ@             # Representive frequency
parms["VLACfg"]        = "@VLACFG@ "           # VLA configuraton ("A", "B", "CnD"...)

# Archive parameters
parms["doLoadArchive"] = True         # Load from archive?
parms["archRoot"]      = "@DATAROOT@" # Root of ASDM/BDF data
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

from EVLACal import EVLACalModel,EVLAStdModel
freq = parms["VLAFreq"]
# Bandpass Calibration
calist = @BPCAL@
BPCals = []
for cal in calist:
    BPCals.append(EVLACalModel(cal))
# Check for standard model
EVLAStdModel(BPCals, freq)
parms["BPCals"]     = BPCals  # Bandpass calibrator(s)
parms["bpChWid2"]   = 5       # Number of channels to smooth in BP calibration                                                                              

# Amp/phase calibration
calist = @PHSCAL@
PCals = []
tcals = []
for cal in calist:
    if not cal in tcals:
        PCals.append(EVLACalModel(cal))
        tcals.append(cal)
# Check for standard model
EVLAStdModel(PCals, freq)
parms["PCals"]          = PCals   # Phase calibrator(s)

calist = @AMPCAL@
ACals = []
tcals = []
for cal in calist:
    if not cal in tcals:
        ACals.append(EVLACalModel(cal))
        tcals.append(cal)
# Check for standard model
EVLAStdModel(ACals, freq)
parms["ACals"]          = ACals  # Amplitude calibrators

calist = @DLYCAL@
DCals = []
tcals = []
for cal in calist:
    if not cal in tcals:
        DCals.append(EVLACalModel(cal))
        tcals.append(cal)
# Check for standard model
EVLAStdModel(DCals, freq)
parms["DCals"]          = DCals      # delay calibrators

parms["refAnt"]        = @REFANT@   # Reference antenna

# Sample spectra
parms["plotSource"]    = @PLOTSRC@          # Source name or None
parms["plotTime"]      = @PLOTTIME@         # timerange

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
if SpanBW<=@VLAFREQ@*parms["MBmaxFBW"]*1.5:
    parms["doMB"] = False

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
parms["doMedn"]        = T                   # Median editing?
parms["doFD1"]         = T                   # Do initial frequency domain flagging
parms["doRMSAvg"]      = T                   # Do RMS/Mean editing for calibrators
parms["doSYGain"]      = T                   # Determine Gains from SY table
parms["doDelayCal"]    = T                   # Group Delay calibration?
parms["doBPCal"]       = T                   # Determine Bandpass calibration
parms["doAmpPhaseCal"] = T                   # Amplitude/phase calibration
parms["doAutoFlag"]    = T                   # Autoflag editing after final calibration?
parms["doSrvrEdt"]     = T                   # Survivor editing
parms["doRecal"]       = T                   # Redo calibration after editing
parms["doMedn2"]       = T                   # redo Median editing?, 2nd pass
parms["doFD2"]         = T                   # Redo initial frequency domain flagging, 2nd pass
parms["doDelayCal2"]   = T                   # Group Delay calibration of averaged data?, 2nd pass
parms["doBPCal2"]      = T                   # Determine Bandpass calibration, 2nd pass
parms["doAmpPhaseCal2"]= T                   # Amplitude/phase calibration, 2nd pass
parms["doAutoFlag2"]   = T                   # Autoflag editing after final calibration?
parms["doSrvrEdt2"]    = T                   # Survivor editing, 2nd pass
parms["doCalAvg"]      = T                   # calibrate and average data
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

parms["XClip"]         = F                   # No X pol clipping                                                                                                       
