# Polarization processing and plotting of MeerKAT Poln data
# Cut and paste relevant sections
import RMFit, Image, ImageUtil, OErr, OSystem, History

# J1709-2226
src = 'J1709-2226';
blc=[421,173,1];trc=[695,649,1]; what="" # In Stokes I
blc2=[421,173,1];trc2=[695,649,1]        # In Stokes I
dir = './'; pdir = './plot/'; sigma_qu = 15.0e-6 # guess
rotate=90.      # Vector rotation angle to "B" vec
rotate=0.       # Vector rotation anglr to "B" vec
PxRange=[40.0e-6, 100.0e-6]    # ppolc range
RMRange=[-40.0,0.0]            # RM range
IPxRange=[-50.0e-6,2000.0e-6]; # Stokes I range
SIRange=[-1.5,0.0]             # Spectral index range
FPRange=[0.0,0.75]             # FPol range 
PxRange=[5.0e-6,40.0e-6]       # PPol range 
clip=50.0e-6; pclip=30.0e-6; iclip=200.0e-6 # pol and IPol clipping levels
pcut = 0.0003; icut=200.0e-6; xinc = 7.; yinc = 7.; factor =20. 

# Load to AIPS from MFImage & RM cubes
si=ObitTask('SubImage');si.DataType='FITS';si.outDType='AIPS'
xfi =  Image.newPFImage('I',dir+iname,0,True,err)
ceni=xfi.Desc.Dict['crpix']
si.inFile=xfi.Fname; si.outName=src;si.outSeq=13;si.outDisk=1
si.BLC=blc;si.TRC=trc; si.outClass='I'; 
si.g
# Total intensity
xi   = Image.newPAImage('I', si.outName, si.outClass, si.outDisk, si.outSeq,True,err)
# Spectral index
si.BLC[2]=2;si.TRC[2]=2; si.outClass='SI'; 
si.g
# Total intensity
xsi   = Image.newPAImage('SI', si.outName, si.outClass, si.outDisk, si.outSeq,True,err)

# Extract polarization planes frm MFImage & RM cubes
xfrm = Image.newPFImage('RM',dir+rmname,0,True,err)
cenrm=xfrm.Desc.Dict['crpix']
si.inFile=xfrm.Fname; si.outClass='PPol';
si.BLC=[blc[0]-int(ceni[0]-cenrm[0]), blc[1]-int(ceni[1]-cenrm[1]), 3]
si.TRC=[trc[0]-int(ceni[0]-cenrm[0]), trc[1]-int(ceni[1]-cenrm[1]), 3]
si.g  # PPol
ppol = Image.newPAImage('PPol', si.outName, si.outClass, si.outDisk, si.outSeq,True,err)
si.inFile=xfrm.Fname; si.outClass='EVPA';
si.BLC=[blc[0]-int(ceni[0]-cenrm[0]), blc[1]-int(ceni[1]-cenrm[1]), 2]
si.TRC=[trc[0]-int(ceni[0]-cenrm[0]), trc[1]-int(ceni[1]-cenrm[1]), 2]
si.g  # EVPA
evpa = Image.newPAImage('EVPA', si.outName, si.outClass, si.outDisk, si.outSeq,True,err)
ImageUtil.PScaleImage(evpa,math.degrees(1),err)  # to degrees
si.inFile=xfrm.Fname; si.outClass='RM';
si.BLC=[blc[0]-int(ceni[0]-cenrm[0]), blc[1]-int(ceni[1]-cenrm[1]), 1]
si.TRC=[trc[0]-int(ceni[0]-cenrm[0]), trc[1]-int(ceni[1]-cenrm[1]), 1]
si.g # RM
xrm = Image.newPAImage('RM', si.outName, si.outClass, si.outDisk, si.outSeq,True,err)

# From individual files
xi  = imlod(src+'_I.fits',0,src,"I",1,1,err)
ppol= imlod(src+'_PPol.fits',0,src,"PPol",1,1,err)
xrm = imlod(src+'_RM.fits',0,src,"RM",1,1,err)
evpa= imlod(src+'_EVPA.fits',0,src,"EVPA",1,1,err)

# if already loaded - MORE WORK
xi   = Image.newPAImage('I',     src+what, "I",    1,1,True,err)
ppol = Image.newPAImage('PPol',  src+what, "PPol", 1,1,True,err)
evpa = Image.newPAImage('EVPA',  src+what, "EVPA", 1,1,True,err)
xrm  = Image.newPAImage('RM',    src+what, "RM",   1,1,True,err)
ppolc= Image.newPAImage('PPolc', src+what, "POLCO",1,1,True,err)
fpol = Image.newPAImage('FPol',  src+what, "FPol", 1,1,True,err)

# Clip RM,EVPA at PPol=pclip, IPol=clip
ppol.GetPlane(None, [1,1,1,1,1], err)
xi.GetPlane(None, [1,1,1,1,1], err)
# Make mask in IPolPix
FArray.PInClip(xi.FArray, -1.0e6, clip, FArray.fblank)
# Make mask in PPolPix
FArray.PInClip(ppol.FArray, -1.0e6, pclip, FArray.fblank)
FArray.PBlank(ppol.FArray, xi.FArray, ppol.FArray)
ppol.PutPlane(None, [1,1,1,1,1], err)
xrm.GetPlane(None, [1,1,1,1,1], err)
# Blank
FArray.PBlank(xrm.FArray, xi.FArray, xrm.FArray)
FArray.PBlank(xrm.FArray, ppol.FArray, xrm.FArray)
xrm.PutPlane(None, [1,1,1,1,1], err)
# Add history
outHistory = History.History("history", xrm.List, err)
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit PolClip",err)
outHistory.WriteRec(-1,"IClip"+" Clip = "+str(clip),err)
outHistory.WriteRec(-1,"PolClip"+" Clip = "+str(pclip),err)
outHistory.Close(err)
# EVPA
evpa.GetPlane(None, [1,1,1,1,1], err);FArray.PBlank(evpa.FArray, ppol.FArray, evpa.FArray)
FArray.PBlank(evpa.FArray, xi.FArray, evpa.FArray); evpa.PutPlane(None, [1,1,1,1,1], err)
# Add history
outHistory = History.History("history", evpa.List, err)
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit PolClip",err); outHistory.WriteRec(-1,"PolClip"+" Clip = "+str(pclip),err)
outHistory.WriteRec(-1,"IClip"+" Clip = "+str(clip),err)
outHistory.Close(err)
# Clip SI below iclip - already clipped
#xi.GetPlane(None, [1,1,1,1,1], err)
# Make mask in IPol at iclip
#FArray.PInClip(xi.FArray, -1.0e6, iclip, FArray.fblank)
#xsi.GetPlane(None, [1,1,1,1,1], err)
# Blank
#FArray.PBlank(xsi.FArray, xi.FArray, xsi.FArray)
#xsi.PutPlane(None, [1,1,1,1,1], err)
# Add history
#outHistory = History.History("history", xsi.List, err)
#outHistory.Open(History.READWRITE, err)
#outHistory.TimeStamp(" Start Obit IClip",err)
#outHistory.WriteRec(-1,"IClip"+" IClip = "+str(iclip),err)
#outHistory.Close(err)

# Polco Q/U for polarization bias correction - lotsa SCREAMING from mama AIPS
# may need si=ObitTask("SubImage");si.outName=src; si.outDisk=1;si.outSeq=1
polco=AIPSTask('polco')
setname(ppol,polco); polco.pixstd=sigma_qu; polco.pcut=1.0; polco.g
ppolc = Image.newPAImage('PPol', si.outName, 'POLCO', si.outDisk, si.outSeq,True,err)
# Fractional poln
comb=AIPSTask('comb'); comb.opcode='DIV';comb.indisk=si.outDisk; comb.outseq=si.outSeq
setname(ppol,comb);set2name(xi,comb); comb.outclass='FPol'
comb.doalign=0. 
comb.aparm[10]=clip # clip on IPol
comb.g
fpol = Image.newPAImage('FPol', si.outName, comb.outclass, si.outDisk, si.outSeq,True,err)
# Blank FPol where xrm blanked
fpol.GetPlane(None, [1,1,1,1,1], err)
xrm.GetPlane(None, [1,1,1,1,1], err)
# Blank
FArray.PBlank(fpol.FArray, xrm.FArray, fpol.FArray)
fpol.PutPlane(None, [1,1,1,1,1], err)

# Save new files (EVPA, RM)
xf=imtab(xi,src+what+".I.fits",0,err)
xf=imtab(xsi,src+what+".SI.fits",0,err)
xf=imtab(ppol,src+what+".PPol.fits",0,err)
xf=imtab(evpa,src+what+".EVPA.fits",0,err)
xf=imtab(xrm,src+what+".RM.fits",0,err)
xf=imtab(fpol,src+what+".FPol.fits",0,err)

# Plots - poln
kntr=AIPSTask('kntr');lwpla=AIPSTask('lwpla')
kntr.dogrey=-1.;kntr.docont=1.;kntr.dovect=1.
setname(xi,kntr); setname(xi,lwpla)
set3name(fpol, kntr); set4name(evpa,kntr)
kntr.ltype=3; kntr.clev=50.0e-6; kntr.rotate=rotate
kntr.levs[1:3]=[-1.,1.]
for i in range(3,15):
    kntr.levs[i] = kntr.levs[i-1]*2.0

kntr.blc[1:4]=blc2; kntr.trc[1:4]=trc2
kntr.pcut=pcut; kntr.icut=icut; kntr.xinc=xinc;kntr.yinc=yinc;kntr.factor=factor
lwpla.outfile=pdir+src+what+'_PolVec.ps'
kntr.g;lwpla.g

# vectors on phlame
kntr.pcut=pcut; kntr.icut=icut; kntr.xinc=xinc;kntr.yinc=yinc;kntr.factor=factor
kntr.dogrey=1.;kntr.docont=-1.;kntr.dovect=1.
kntr.ofmfile ='./RYPHLAME.000'
kntr.pixrange[1:3]=IPxRange
kntr.g;lwpla.g

# RM
kntr=AIPSTask('kntr');lwpla=AIPSTask('lwpla')
kntr.dogrey=2.;kntr.docont=1.;kntr.dovect=-1.
setname(xi,kntr); setname(xi,lwpla)
set2name(xrm, kntr); 
kntr.ltype=3; kntr.clev=100.0e-6; 
kntr.levs[1:3]=[-1.,1.]
for i in range(3,10):
    kntr.levs[i] = kntr.levs[i-1]*4.

kntr.blc[1:4]=blc2; kntr.trc[1:4]=trc2
kntr.pixrange[1:3]=RMRange; kntr.ofmfile='./RAINBOW.001'
lwpla.outfile=pdir+src+what+'_RM.ps'
kntr.clev=50.0e-6; 
kntr.g;lwpla.g

# Poln intensity w/ Stokes I contours, FPOL
kntr=AIPSTask('kntr');lwpla=AIPSTask('lwpla')
kntr.dogrey=2.;kntr.docont=1.;kntr.dovect=-1.
setname(xi,kntr); setname(xi,lwpla)
set2name(fpol, kntr); 
kntr.ltype=3; kntr.clev=100.0e-6; kntr.ofmfile ='./RYPHLAME.000'
kntr.ofmfile='./RAINBOW.001'
kntr.levs[1:3]=[-1.,1.]
for i in range(3,10):
    kntr.levs[i] = kntr.levs[i-1]*4.0

kntr.blc[1:4]=blc2; kntr.trc[1:4]=trc2
kntr.pixrange[1:3]=FPRange; 
lwpla.outfile=pdir+src+what+'_FPol.ps'
kntr.clev=50.0e-6;
kntr.g;lwpla.g

# Stokes I contours w/ poln vec
kntr=AIPSTask('kntr');lwpla=AIPSTask('lwpla')
kntr.dogrey=-1.;kntr.docont=1.;kntr.dovect=1.
setname(xi,kntr); setname(xi,lwpla)
set3name(fpol, kntr); set4name(evpa,kntr)
kntr.pcut=pcut; kntr.icut=icut; kntr.xinc=xinc;kntr.yinc=yinc;kntr.factor=factor
kntr.ltype=3; kntr.clev=25.0e-6; kntr.rotate=rotate
kntr.levs[1:3]=[-1.,1.]
for i in range(3,20):
    kntr.levs[i] = kntr.levs[i-1]*(2.0)**0.5

kntr.blc[1:4]=blc2; kntr.trc[1:4]=trc2
kntr.pixrange[1:3]=PxRange; 
lwpla.outfile=pdir+src+what+'_ContVec.ps'
kntr.g;lwpla.g

# SI
kntr.dogrey=2.;kntr.docont=1.;kntr.dovect=-1.
setname(xi,kntr); setname(xi,lwpla)
set2name(xsi, kntr); 
kntr.ltype=3; kntr.clev=100.0e-6; 
kntr.levs[1:3]=[-1.,1.]
for i in range(3,10):
    kntr.levs[i] = kntr.levs[i-1]*4.

kntr.blc[1:4]=blc2; kntr.trc[1:4]=trc2
kntr.pixrange[1:3]=SIRange; kntr.ofmfile='./RAINBOW.001'
lwpla.outfile=pdir+src+what+'_SI.ps'
kntr.clev=50.0e-6; 
kntr.g;lwpla.g

# color image - all
kntr=AIPSTask('kntr');lwpla=AIPSTask('lwpla')
kntr.dogrey=1.;kntr.docont=-1.;kntr.dovect=-1.
setname(xi,kntr); setname(xi,lwpla)
#kntr.ltype=3; kntr.functype='NG'
kntr.ltype=3; kntr.functype='  '
kntr.blc[1:4]=blc2; kntr.trc[1:4]=trc2
kntr.pixrange[1:3]=IPxRange; kntr.ofmfile=''
kntr.ofmfile='./RYPHLAME.000'
lwpla.outfile=pdir+src+what+'_Color.ps'
#lwpla.outfile=pdir+src+what+'_BW.ps'
kntr.g;lwpla.g


# Stars
star=AIPSTask('stars'); 
setname(xi,star);
star.intext='./software/Slice_stars.text'
star.stvers=1; kntr.stvers = star.stvers=1; kntr.stfactor=-1.
