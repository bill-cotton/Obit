# average Stokes I channel weights in a set of images
#exec(open('AvgWeights.py').read())
# Saves information in dict saved in AvgWeights.pickle
# weights, one per subband, RMSes:[mosaic_name,RMS], 
# Freqs:[center freq, MHz], EffFreq: Effective frequency,
# dir:data directory, post: ending of FITS file name
print("\nDetermining weighting, effective frequency")

# Assumes MFImage output
dir = "../data/"      # Data directory
# file name roots
posn = ["Field_1" ,"Field_2", "Field_3", "Field_4"]
post = '_I_Mosaic.fits' # file=dir+<name root>+post
doLin = True  # If True weight by 1/sigma, else 1/sigma^2
nthreads = 16  # Number of threads to use

import Image, OErr, OSystem, FArray
OSystem.PAllowThreads(nthreads)  # with threads
from PipeUtil import SaveObject, FetchObject
import math, pickle

# Function to get RMSHist for plane planeno
def imrms(inImage, planeno=1):
    """ Get plane RMSHist
    Returns dictionary with statistics of selected region with entries:
    * inImage   = Python Image object, created with getname, getFITS
    """
    inImage.GetPlane(None,[planeno,1,1,1,1],err)
    rms = inImage.FArray.RMS
    return rms
# end imrms

# Get information from first file
x=getFITS(dir+posn[0]+post,0)
# test that MFImage output
if x.Desc.Dict['ctype'][2]!='SPECLNMF':
    raise RuntimeError("not MFImage output")

nterm  = x.Desc.List.Dict['NTERM'][2][0]
nspec  = x.Desc.List.Dict['NSPEC'][2][0]
nplane = x.Desc.Dict['inaxes'][2]

# Get frequencies
freqs = []
for i in range(1,nspec+1):
    key = "FREQ%4.4d"%i
    freqs.append(x.Desc.List.Dict[key][2][0]*1.0e-6)

# end channel loop

# Channel RMSes per image stored in result[p]
result = {}
for i in range (0,len(posn)):
    p=posn[i]
    x=getFITS(dir+p+post,0)
    print ("File",dir+p+post)
    rms = []
    for j in range(nterm+1,nplane+1):
        s = imrms(x,j)
        rms.append(s)
    
    # end channel loop
    result[p] = rms
# end image loop

# broadband rmses
ff = []
for i in range (0,len(posn)):
    p=posn[i];  x=getFITS(dir+p+post,0)
    s = imrms(x, 1)
    ff.append((p,s))

# Average weights
for i in range (0,len(posn)):
    p=posn[i];  wt = []; rms = result[p]
    for j in range(0,len(rms)):
        if rms[j]>0.0:
            if doLin:
                wt.append(1.0e-6*(rms[j]**-1))  # weight by 1/sigma
            else:
                wt.append(1.0e-6*(rms[j]**-2)) # weight by 1/sigma^2
        else:
            wt.append(0.0)

# end channel loop
# Normalize by sum of weights
sumwt = sum(wt)
for i in range(0,len(wt)):
    wt[i] = 100*wt[i]/sumwt

# Print subband weights
print ("Subband  Freq   Weight")
for i in range(0,len(wt)):
    print ("%3d %9.2f %8.2f"%(i+1,freqs[i],wt[i]))

# Effective frequency
sum1 = 0.0; sumwt=0.0
for i in range(0,len(wt)):
    sum1 += freqs[i]*wt[i]; sumwt += wt[i];

EffFreq = sum1/sumwt
print ("\nEffective Frequency = %9.2f"%EffFreq,"MHz")

# Save to pickle jar
# weights, one per subband, RMSes:[mosaic_name,RMS], Freqs:[center freq, MHz],
# dir:data directory, post: ending of FITS file name
stuff = {"weights":wt, "RMSes":ff, "Freqs":freqs, "EffFreq":EffFreq, \
         "dir":dir, "post":post}
SaveObject(stuff, "AvgWeights.pickle", True)
print ("Wrote AvgWeights.pickle" )

