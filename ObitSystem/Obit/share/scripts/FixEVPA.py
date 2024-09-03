# Generate a BP table to correct errors in EVPA derived from a scan on 3C48
# Works on data in a circular basis
# Needs:
# - Q and U CLEAN FITS images of 3C48 with enhanced spectral resolution
#   derived applying a previously applied, but inadequate, polarization
#   calibration. 
# - Partially poln calibrated UV data set in AIPS FITAB format with a
#   BP table to which corrections can be applied. 
# Execute as
#exec(open('FixEVPA.py').read())
import UV, Image, OErr, FArray, ImageDesc, UVPolnUtil, History
import math
err = OErr.OErr()
FBlank = FArray.fblank  # Blanked pixel value

# UV data file name in cwd, AIPS UVTAB format
uv_file='My_Data.uvtab'
BPinVer  = 2  # Input BP table version
BPoutVer = 3  # Output BP table version

# 3C48 Q&U FITS images in cwd
Q_file='3C48_QPol.fits'
U_file='3C48_UPol.fits'

# 3C48 Polarization model from Perley & Butler 2013, ApJS, 206,16
cal = '3C48'
EVPA_cal = math.radians(122.)  # EVPA (deg) @ 0 wavelength
RM_cal   = -68.                # Rotation measure (rad/m^2)
ra  = ImageDesc.PHMS2RA("01:37:41.299")    # Calibrator RA
dec = ImageDesc.PDMS2Dec("33:09:35.133")   # Calibrator Dec

# Shouldn't need to edit below here:
#-----------------------------------------------------------------------
# files to Obit objects
uv = UV.newPFUV("in UV", uv_file, 0, True, err)
# Open and close to check
uv.Open(UV.READONLY, err); uv.Close(err);
uv_desc = uv.Desc.Dict   # UV data header as python dict
bpin  = uv.NewTable(Table.READONLY, "AIPS BP", BPinVer, err)
numPol = min(2, uv_desc['inaxes'][uv_desc['jlocs']])
numIF = uv_desc['inaxes'][uv_desc['jlocif']]
numChan=uv_desc['inaxes'][uv_desc['jlocf']]
bpout = uv.NewTable(Table.WRITEONLY, "AIPS BP", BPoutVer, err,
                    numPol=numPol, numIF=numIF, numChan=numChan)
bpout = Table.PClone(bpin, bpout)  # Clone input table
OErr.printErrMsg(err,message='Error with uv data')  # Error check

# Images
im_Q = Image.newPFImage("Q", Q_file, 0, True, err)
im_U = Image.newPFImage("U", U_file, 0, True, err)
# Open and close to check
im_Q.Open(Image.READONLY, err); im_Q.Close(err);
im_U.Open(Image.READONLY, err); im_U.Close(err);
OErr.printErrMsg(err,message='Error with image')  # Error check
Q_desc = im_Q.Desc.Dict       # Q image header
Q_info = im_Q.Desc.List.Dict  # List of frequency bins
nterm = Q_info['NTERM'][2][0]  # no. spectral terms in images
nspec = Q_info['NSPEC'][2][0]  # Number of spectral planes
# Pixel number of position of calibrator
pix = ImageDesc.PGetPixel(im_Q.Desc, [ra,dec], err)
pixel = [int(pix[0]+0.5), int(pix[1]+0.5)]  # Round to nearest pixel

# Get list of subband bins as lambda^2, 
# low, center, high (freq), Q, U, corr
clight = 2.997924562e8   # Speed of light m/s
bins = []
for i in range(0,nspec):
    key = "FREL%4.4d"%(i+1);  lof = Q_info[key][2][0]
    key = "FREQ%4.4d"%(i+1); cenf = Q_info[key][2][0]
    key = "FREH%4.4d"%(i+1);  hif = Q_info[key][2][0]
    # center to wavelength^2
    cenl2 = (clight/cenf)**2
    bins.append([lof,cenl2,hif, FBlank, FBlank, FBlank])

# end loop

# Get subband flux densities at Q,U position of peak
for i in range(0,nspec):
    ip = i+nterm+1  # 1-rel Plane number in cube
    im_Q.GetPlane(None, [ip,1,1,1,1], err)
    Q_flux = im_Q.FArray.get(pixel[0], pixel[1])
    im_U.GetPlane(None, [ip,1,1,1,1], err)
    U_flux = im_U.FArray.get(pixel[0], pixel[1])
    # Valid data?
    if (im_Q.FArray.RMS!=0.0) and (im_U.FArray.RMS!=0.0):
        bins[i][3] = Q_flux; bins[i][4] = U_flux
        # Get correction (rad)
        evpa_obs = 0.5*math.atan2(U_flux,Q_flux)
        evpa_mod = EVPA_cal + bins[i][1]*RM_cal
        bins[i][5] = -2*(evpa_obs-evpa_mod)  # Correction 
        #print (i, "correction=",math.degrees(bins[i][5]), bins[i][1]) # debug
# End loop
OErr.printErrMsg(err,message='Error reading images')  # Error check

# Fill in missing values, first valid
last_corr = FBlank
for i in range(0,nspec):
    if bins[i][5]!=FBlank:
        last_corr = bins[i][5]
        break

# Fill in with last_corr
for i in range(0,nspec):
    if bins[i][5]==FBlank:
        bins[i][5] = last_corr
    else:
        last_corr = bins[i][5]

# List of Frequencies per IF/channel
freqs = UVPolnUtil.GetFreqArr(uv, err)
OErr.printErrMsg(err,message='Error reading frequencies')  # Error check

# Get list of corrections (rlcorr) per channel
nch = len(freqs); last_bin = 0
rlcorr_c = nch*[0.0]; rlcorr_s = nch*[0.0]
for i in range(0,nch):
    f = freqs[i]
    if f<bins[last_bin][0]:
        # in an earlier bin, search from start
        last_bin = 0
    # Search if needed
    while not ((f>=bins[last_bin][0]) and (f<=bins[last_bin][2])):
        last_bin += 1
    rlcorr_c[i] = math.cos(bins[last_bin][5]);
    rlcorr_s[i] = math.sin(bins[last_bin][5]);

# Adjust BP table BPinVer to BPoutVer
bpin.Open(Table.READONLY, err); bpout.Open(Table.WRITEONLY,err)
nrow = bpin.Desc.Dict['nrow']
for irow in range(1, nrow+1):
    # Update row REAL 2 and IMAG 2 (R-L phase)
    row=bpin.ReadRow(irow,err)
    for i in range(0,nch):
        if (row['REAL 2'][i]!=FBlank):
            g_in  = complex(row['REAL 2'][i], row['IMAG 2'][i]);
            g_add = complex(rlcorr_c[i],rlcorr_s[i])
            g_out = g_in * g_add   # Update
            row['REAL 2'][i] = g_out.real; row['IMAG 2'][i] = g_out.imag;
        # End if input valid
    # End loop over channels
    bpout.WriteRow(irow,row,err)
# End loop over rows

bpin.Close(err); bpout.Close(err)
OErr.printErrMsg(err,message='Error updating BP table')  # Error check

# Add history
pgm = "FixEVPA"
outHistory = History.History("history", uv.List, err)
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit "+pgm,err)
outHistory.WriteRec(-1,pgm+" BPIn   = "+str(BPinVer),err)
outHistory.WriteRec(-1,pgm+" BPOut  = "+str(BPoutVer),err)
outHistory.WriteRec(-1,pgm+" Q_file = "+Q_file,err)
outHistory.WriteRec(-1,pgm+" U_file = "+U_file,err)
outHistory.WriteRec(-1,pgm+" cal    = "+cal,err)
outHistory.WriteRec(-1,pgm+" EVPA_cal = "+str(math.degrees(EVPA_cal)),err)
outHistory.WriteRec(-1,pgm+" RM_cal   = "+str(RM_cal),err)
outHistory.Close(err)

