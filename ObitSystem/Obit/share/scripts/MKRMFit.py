# exec(open('MKRMFit.py').read())
# RM analysis for Q and U cubes produced by MFImage
import OErr, OSystem, UV, AIPS, FITS
# Do RM fitting
src = ['Abell_85_']; post=''

import Image,RMFit,OSystem,OErr
err = OErr.OErr()
OSystem.PAllowThreads(24)  # 24 threads
for s in src:
    fit = RMFit.PCreate('Fitter')
    fit.List.set('doError',True);   fit.List.set('doRMSyn',True)
    #fit.List.set('minRMSyn',-1500.); fit.List.set('maxRMSyn',+1500.)
    #fit.List.set('minRMSyn',-400.); fit.List.set('maxRMSyn',+400.)
    fit.List.set('minRMSyn',-150.); fit.List.set('maxRMSyn',+150.)
    fit.List.set('delRMSyn',0.5);   fit.List.set('maxChi2',10000.0)
    fit.List.set('minQUSNR',1.0);   fit.List.set('minFrac',0.25);   
    fit.List.set('refLamb2',1.0e-6);  
    print ('do',s)
    inQ = Image.newPFImage('Q', s+'Q'+post+'.fits', 0, True,err)
    inU = Image.newPFImage('U', s+'U'+post+'.fits', 0, True,err)
    outRM = Image.newPFImage('RM', s+'RM_RMSyn'+post+'.fits',0,False,err)
    #fit.List.set('doRMSyn',False)
    #outRM = Image.newPFImage('RM', s+'RM_LSQ'+post+'.fits',0,False,err)
    fit.Cube(inQ, inU, outRM, err)
    del fit, inQ, inU, outRM

