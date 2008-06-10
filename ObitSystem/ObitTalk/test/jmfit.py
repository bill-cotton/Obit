from AIPS import AIPS
from AIPSTask import AIPSTask
from AIPSData import AIPSImage

AIPS.userno = 1999

image = AIPSImage('MANDELBROT', 'MANDL', 1, 1)
if image.exists():
    image.zap()

mandl = AIPSTask('mandl')
mandl.outdata = image
mandl.imsize[1:] = [ 64, 64 ]
mandl.go()

try:
    jmfit = AIPSTask('jmfit')
    jmfit.indata = image
    jmfit.ngauss = 4
    jmfit.domax[1:] = [1, 0, 0, 0]
    jmfit.go()
    print 'Peak values:', jmfit.fmax[1:]
    for fmax in jmfit.fmax[1:]:
        assert(fmax)
finally:
    image.zap()
