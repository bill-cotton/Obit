from __future__ import absolute_import
from __future__ import print_function
import OErr
import OSystem
import AIPS
import Image
import OTObit

import numpy as np


err = OErr.OErr()
user = 100

ObitSys = OSystem.OSystem('Pipeline', 1, user, 1, ['./aipsdisk'], \
                          1, ['./'], True, False, err)
#OSystem.PAllowThreads(8)
#ObitTalkUtil.SetEnviron(aipsdirs=[(None, './aipsdisk')], fitsdirs=[(None, './')])
user = OSystem.PGetAIPSuser()
AIPS.userno = user

x = Image.newPAImage("TEST", "TEST", "FIT", 1, 1, True, err)
x.GetPlane(None, [1,1,1,1,1], err)
buf=x.FArray.Buf
nparr = np.frombuffer(buf, dtype=np.float32)
nparr[np.where(nparr==3140.8928)]=np.nan
print("numpy Mean", np.nanmean(nparr))
z=OTObit.imstat(x) # should give same mean
print("imstat Mean", z['Mean'])
OErr.printErr(err)
OSystem.Shutdown(ObitSys)
