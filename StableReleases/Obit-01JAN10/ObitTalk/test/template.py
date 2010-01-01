from AIPS import AIPS
from AIPSTask import AIPSTask
from AIPSData import AIPSImage
from ObitTask import ObitTask

AIPS.userno = 103

image = AIPSImage('MANDELBROT', 'MANDL', 1, 1)

mandl = AIPSTask('mandl')
mandl.outdata = image
mandl.imsize[1:] = [ 512, 512 ]
mandl.go()

try:
    template = ObitTask('Template')
    template.DataType = 'AIPS'
    template.inName = image.name
    template.inClass = image.klass
    template.inDisk = image.disk
    template.inSeq = image.seq
    template.go()
finally:
    image.zap()
