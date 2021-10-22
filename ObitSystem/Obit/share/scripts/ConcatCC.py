# Script to append with shifts from a set of facets of a mosaic to the output
import sys, Obit, Image, UV, OErr, TableUtil, ObitTask

name = 'NGC1532'; 
uvFile='NGC1532_Ref59.uvtab'; 
outFile='NGC1532.IPol.testCC.fits'
nfield=268; dir = 'Model/'

uvData = UV.newPFUV('data',uvFile, 0,True,err)
outImage = Image.newPFImage('out', outFile, 0,True,err)
# Copy first
taco=ObitTask.ObitTask('TabCopy'); taco.inTab='AIPS CC'; taco.DataType='FITS'
taco.inFile = dir+name+'.IM%4.4d.fits'%1
taco.outFile=outFile
taco.g

inCCVer=1; outCCVer=1
for i in range(2,nfield+1):
    f = dir+name+'.IM%4.4d.fits'%i
    x = Image.newPFImage('in',f,0,True,err)
    TableUtil.PTableCCAppendShift(x, inCCVer, outImage, outCCVer, uvData, err)

