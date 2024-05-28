# Descimate an image
xinc=10; yinc=10
inIm=getn(12); outIm=getn(15) # G53.6-2.2
inIm=getn(16); outIm=getn(19) # G36.6+2.6
inIm=getn(17); outIm=getn(20) # G55.7+3.4
inIm=getn(18); outIm=getn(21) # G57.2+0.8
inIm.GetPlane(None,[1,1,1,1,1], err); outIm.GetPlane(None,[1,1,1,1,1], err); 
(nxi,nyi) =  inIm.Desc.Dict['inaxes'][0:2]
(nxo,nyo) = outIm.Desc.Dict['inaxes'][0:2]
import FArray
inArr = inIm.FArray; outArr = outIm.FArray
iiy = -1
for iy in range(0,nyi,yinc):
    iiy += 1; iix = -1
    if iiy>=nyo:
        break
    for ix in range(0,nxi,xinc):
        iix += 1;
        if iix>=nxo:
            break
        val = inArr.get(ix,iy)
        outArr.set(val,iix,iiy)

# end loop
outIm.PutPlane(None,[1,1,1,1,1], err); 
OErr.printErr(err)
