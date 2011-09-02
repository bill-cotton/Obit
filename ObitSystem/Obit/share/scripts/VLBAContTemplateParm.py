# Template project parameter file VLBAContPipe
# Generate parameter file using VLBACal.VLBAMakeParmFile
#
# Substitutions surrounded by 'at' characters
# PROJECT     Project name (up to 12 char)
# SESSION     Session code
# BAND        Band code
# UVFITS      Name of uvfits file in $FITS
# IDIFITS     List of IDIFITS files in $FITS 
# CALINT      CL table interval in min
# DESTDIR     Output directory 
# ARCHFILEID  Archive file ID
# DOLOADIDI   Load IDI FITS  files
# DOLOADFITS  Load UV FITS files
#--------------------------------------------------------------------------------------------
project       = "@PROJECT@"                # Project name (12 char or less, used as AIPS Name)
session       = "@SESSION@"                # Project session code
band          = "@BAND@"                   # Observing band
archFileID    = @ARCHFILEID@               # Archive File ID
logFile       = project+"_"+session+"_"+band+".log"  # Processing log file
# Destination directory for copying output files.  empty string -> do not copy
parms["copyDestDir"]   = "@DESTDIR@/" + project + session + band 

# NOTE: these files must be in $FITS directory!
dataInUVF     = "@UVFITS@"    # UVFITS data file name
dataInIDI     = @IDIFITS@     # List of IDIFITS data files

calInt        = @CALINT@            # Calibration table interval in min.
Compress      = True                # Use compressed UV data?

# Quantization correction
parms["QuantFlag"]   = 0.8          # If >0, flag solutions < QuantFlag (use 0.9 for 1 bit, 0.8 for 2 bit)

# Specify calibration/target sources. (When initialized, both parameters
# default to all sources, so leave commented here until ready to use.)
#parms["contCals"]    = []           # list of calibrators
#parms["targets"]     = []           # targets, empty = all
parms["contCalModel"] = None        # List of models
parms["refAnts"]      = [2,4,5,8,9] # List of acceptable reference antennas for fringe fitting 

# Final Image/Clean
parms["Stokes"]          = "I"      # Stokes to image

# Control
T   = True
F   = False
check                    = F       # Only check script, don't execute tasks
debug                    = F       # run tasks debug
doLoadIDI                = @DOLOADIDI@       # Load data from IDI FITS?, else already in AIPS?
doLoadUVF                = @DOLOADUVF@       # Load the "AIPS Friendly" UV FITS  version?
parms["doClearTab"]      = T       # Clear cal/edit tables
parms["doCopyFG"]        = T       # Copy FG 1 to FG 2
parms["doEditList"]      = T       # Edit using editList?
parms["doQuack"]         = T       # Quack data?
parms["doQuantCor"]      = T       # Quantization correction/flagging?
parms["doPACor"]         = T       # Make parallactic angle correction?
parms["doOpacCor"]       = T       # Make Opacity/Tsys/gain correction?
parms["doFindOK"]        = T       # Search for OK cals if contCals not given
parms["doFindCal"]       = T       # Search for best calibration/reference antenna
parms["doPCcor"]         = T       # Apply PC table?
parms["doManPCal"]       = T       # Determine and apply manual phase cals?
parms["doBPCal"]         = T       # Determine Bandpass calibration
parms["doImgCal"]        = T       # Image calibrators
parms["doDelayCal"]      = T       # Determine/apply delays from contCals
parms["doAmpCal"]        = T       # Determine/smooth/apply amplitudes from contCals
parms["doCalAvg"]        = T       # calibrate and average
parms["doImgTarget"]     = T       # Image targets?
parms["doPhaseCal"]      = T       # Phase calibrate all data with self-cal?
parms["doImgFullTarget"] = T       # Final Image/Clean/selfcal ?
parms["doSaveUV"]        = T       # Save UV (calibrated and averaged) results?
parms["doSaveImg"]       = T       # Save image results?
parms["doSaveTab"]       = T       # Save calibration and editing tables?
parms["doCleanup"]       = T       # Cleanup AIPS directories?

# diagnostics/reports
parms["doSNPlot"]        = T       # Plot SN tables etc
parms["doPCPlot"]        = T       # Plot PC results?
parms["doSpecPlot"]      = T       # Plot the amp. and phase across the spectrum
parms["doDiagPlots"]     = T       # Source plots
parms["doKntrPlots"]     = T       # Contour plots
parms["doMetadata"]      = T       # Generate metadata dictionaries
parms["doVOTable"]       = T       # VOTable report
parms["doHTML"]          = T       # Generate HTML report?
