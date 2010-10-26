# Project parameter file for EVLAPipeline
project       = "BL149"                    # Project name (12 char or less, used as AIPS Name)
session       = "AJ"                       # Project session code
band          = "Ku"                       # Observing band
logFile       = project+"_"+session+"_"+band+".log"  # Processing log file

dataInUVF     = "BL149AJ.uvfits"       # Input uvfits data file name
# NOTE: this REALLY HAS TO BE IN $FITS!!!!!!
calInt        = 10.0/60.                # Calibration table interval in min.
Compress      = True                    # Use compressed UV data?

# Special editing list
doEditList  = False        # Edit using editList?
editFG      = 2           # Table to apply edit list to
editList = [
    # Ft Davis, LPol all, dead
    #{"timer":("0/0:0:0.0","10/0:0:0.0"),"Ant":[2,0],"IFs":[1,2],"Chans":[1,0],"Stokes":'0111',"Reason":"dead rcvr"},
    ]

# Quantization correction
doQuantCor  = True      # Do quantization correction
QuantSmo    = 0.5       # Smoothing time (hr) for quantization corrections
QuantFlag   = 0.8       # If >0, flag solutions < QuantFlag

# Find good calibration data
doFindCal    = True         # Search for good calibration/reference antenna
contCals     = ["1253-055","1510-089","2005+403","2200+420","2251+158"]       # List of cals
contCalModel = None                                                     # List of models
targets     = []               # targets
goodCalPicklefile = "./"+project+"_"+session+"_"+band+"GoodCal.pickle"  # Where results saved
refAnt        = 4                       # Reference antenna (used?) 
refAnts       = [4,5,8,9]               # List of acceptable reference antennas for fringe fitting 

# 'Manual" phase cal - even to tweak up PCals
doManPCal      = True      # Determine and apply manual phase cals?

# Bandpass Calibration
doBPCal       = True                    # Bandpass calibration
bpsolint1     = 20.0/60.0               # BPass phase correction solution in min
bpsolint2     = 10.0                    # BPass bandpass solution in min
bpBChan1      = 1                       # Low freq. channel,  initial cal
bpEChan1      = 0                       # Highest freq channel, initial cal, 0=>all
bpDoCenter1   = 0.5                     # Fraction of  channels in 1st, overrides bpBChan1, bpEChan1
bpBChan2      = 1                       # Low freq. channel for BP cal
bpEChan2      = 0                       # Highest freq channel for BP cal,  0=>all 
bpChWid2      = 31                      # Number of channels in running mean BP soln
bpdoAuto      = True                    # Use autocorrelations rather than cross?

# Editing
doClearTab    = True                    # Clear cal/edit tables
doGain        = True                    # Clear SN and CL tables >1
doFlag        = True                    # Clear FG tables > 1
doBP          = True                    # Clear BP tables

doCopyFG      = True                    # Copy FG 1 to FG 2

doQuack       = False                   # Quack data?
doQuack       = True                    # Quack data?
quackBegDrop  = 5.0/60.                 # Time to drop from start of each scan in min
quackEndDrop  = 0.0                     # Time to drop from end of each scan in min

# Imaging sources
doImgCal    = True         # Image calibrators
Robust      = 0.0          # Weighting robust parameter
FOV         = 20.0e-3/3600 # Field of view radius in deg.
Niter       = 500          # Max number of clean iterations
minSNR      = 3.0          # Minimum Allowed SNR

# Delay calibration
doDelayCal  = True         # Determine/apply delays from contCals

# Amplitude calibration
doAmpCal    = True         # Determine/smooth/apply amplitudes from contCals

# Check for calibrator models
contCalModel = VLBAImageModel(contCals, outIclass, disk, seq, err)

# Phase calibration of all targets in averaged calibrated data 
doPhaseCal    = True       # Phase calibrate all data with self-cal?

# Instrumental polarization cal?
doInstPol     = True      # determination instrumental polarization from instPolCal
instPolCal    = None      # Defaults to contCals

# Right-Left phase (EVPA) calibration 
doRLCal      = True          # Set RL phases from RLCal - also needs RLCal
# if given, a triplet, (name, R-L phase(deg@1GHz), RM)
# interpolated from VLA calibration
RLCal = [  \
    ("1253-055", -8.,0.0), \
        ("2200+420", 16.,0.0), \
        ("2251+158",120.,0.0), \
        ]

# Final Image/Clean
doImgFullTarget = True    # Final Image/Clean/selfcal
Stokes          = "IQU"   # Stokes to image

# Control
T   = True
F   = False
check           = F       # Only check script, don't execute tasks
debug           = F       # run tasks debug
doLoadIDI       = F       # Load data from IDI FITS?, else already in AIPS?
doLoadUVF       = F       # Load the "AIPS Friendly" UV FITS  version?
doClearTab      = F       # Clear cal/edit tables
doCopyFG        = F       # Copy FG 1 to FG 2
doEditList      = F       # Edit using editList?
doQuack         = F       # Quack data?
doQuantCor      = F       # Quantization correction/flagging?
doPACor         = F       # Make parallactic angle correction?
doOpacCor       = F       # Make Opacity/Tsys/gain correction?
doFindCal       = F       # Search for best calibration/reference antenna
doPCcor         = F       # Apply PC table?
doManPCal       = F       # Determine and apply manual phase cals?
doBPCal         = F       # Determine Bandpass calibration
doImgCal        = F       # Image calibrators
doDelayCal      = F       # Determine/apply delays from contCals
doAmpCal        = F       # Determine/smooth/apply amplitudes from contCals
doCalAvg        = T       # calibrate and average
doImgTarget     = T       # Image targets?
doPhaseCal      = T       # Phase calibrate all data with self-cal?
doInstPol       = T       # determination instrumental polarization from instPolCal 
doRLCal         = T       # Set RL phases from RLCal - also needs RLCal
doImgFullTarget = F       # Final Image/Clean/selfcal 
doSaveUV        = F       # Save UV (calibrated and averaged) results
doSaveImg       = F       # Save image results
doSaveTab       = F       # Save calibration and editing tablesT
doCleanup       = F       # Cleanup AIPS direstories?

# diagnostics
doSNPlot        = T       # Plot SN tables etc
doPCPlot        = T       # Plot PC results?

