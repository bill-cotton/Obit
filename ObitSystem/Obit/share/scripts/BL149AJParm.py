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

# Control
T   = True
F   = False
check         = F       # Only check script, don't execute tasks
debug         = F       # run tasks debug
doLoadIDI     = F       # Load data from IDI FITS?, else already in AIPS?
doLoadUVF     = T       # Load the "AIPS Friendly" UV FITS  version?
doClearTab    = T       # Clear cal/edit tables
doCopyFG      = T       # Copy FG 1 to FG 2
doEditList    = T       # Edit using editList?
doQuack       = T       # Quack data?
doQuantCor    = T       # Quantization correction/flagging?
doPACor       = T       # Make parallactic angle correction?
doOpacCor     = T       # Make Opacity/Tsys/gain correction?
doFindCal     = T       # Search for best calibration/reference antenna
doPCcor       = T       # Apply PC table?
doManPCal     = T       # Determine and apply manual phase cals?
doBPCal       = T       # Determine Bandpass calibration
doImgCal      = T       # Image calibrators
doDelayCal    = T       # Determine/apply delays from contCals
doAmpCal      = T       # Determine/smooth/apply amplitudes from contCals
doImgTarget   = T       # Image targets?
doCalAvg      = T       # calibrate and average
doSaveUV      = T       # Save UV (calibrated and averaged) results
doSaveImg     = T       # Save image results
doSaveTab     = T       # Save calibration and editing tables

# diagnostics
doCleanup     = T       # Cleanup AIPS direstories?
doSNPlot      = T       # Plot SN tables etc
doPCPlot      = T       # Plot PC results?
