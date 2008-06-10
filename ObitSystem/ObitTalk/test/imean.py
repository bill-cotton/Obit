from AIPS import AIPS
from AIPSTask import AIPSTask
from AIPSData import AIPSImage

AIPS.userno = 1999

image = AIPSImage('MANDELBROT', 'MANDL', 1, 1)
if image.exists():
    image.zap()

mandl = AIPSTask('mandl')
mandl.outdata = image
mandl.imsize[1:] = [ 512, 512 ]
mandl.go()

try:
    imean = AIPSTask('imean')
    imean.indata = image
    imean.go()
    print 'Average: %f, RMS noise: %f' % (imean.pixavg, imean.pixstd)
    assert(imean.pixavg)
    assert(imean.pixstd)
finally:
    image.zap()
