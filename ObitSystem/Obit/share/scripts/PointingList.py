# Convert pointing list to catalog of pointings
# generates dict cat name:(ra, dec, factor, offx, offy)
# FITS file names <file root>.<stktrans>.fits (stktrans from Targ.py)
# FITS input can be gzipped
#  Pointing center positions are as 'hh:mm:ss.s", "dd:mm:ss"
#  factor is an addition weight scaling factor
#  offx, offy offsets in ra,dec (asec) (observed-unbiased)
# Example:
#           file root       RA(J2000)         Dec(J2000) factor off_ra off_dec
data = [ \
         ('06000+17576', '06:00:00.00000', '17:57:36.0000 ', 1.0, 0.0, 0.0), \
         ('06000+18243', '06:00:00.00000', '18:24:18.0000 ', 1.0, 0.0, 0.0), \
         ('06000+19178', '06:00:00.00000', '19:17:48.0000 ', 1.0, 0.0, 0.0), \
         ('06000+19446', '06:00:00.00000', '19:44:36.0000 ', 1.0, 0.0, 0.0), \
         ('06000+18510', '06:00:00.00000', '18:51:00.0000 ', 1.0, 0.0, 0.0), \
         ('06000+20383', '06:00:00.00000', '20:38:18.0000',  1.0, 0.0, 0.0), \
         ('06015+19580', '06:01:30.00000', '19:58:00.0000 ', 1.0, 0.0, 0.0), \
         ('06015+19312', '06:01:30.00000', '19:31:12.0000 ', 1.0, 0.0, 0.0), \
         ('06015+19044', '06:01:30.00000', '19:04:24.0000 ', 1.0, 0.0, 0.0), \
         ('06015+18377', '06:01:30.00000', '18:37:42.0000 ', 1.0, 0.0, 0.0), \
         ('06015+18109', '06:01:30.00000', '18:10:54.0000 ', 1.0, 0.0, 0.0), \
         ('06030+17576', '06:02:60.00000', '17:57:36.0000 ', 1.0, 0.0, 0.0), \
         ('06030+18243', '06:02:60.00000', '18:24:18.0000 ', 1.0, 0.0, 0.0), \
         ('06030+18510', '06:02:60.00000', '18:51:00.0000 ', 1.0, 0.0, 0.0), \
         ('06030+19178', '06:02:60.00000', '19:17:48.0000 ', 1.0, 0.0, 0.0), \
         ('06030+19446', '06:02:60.00000', '19:44:36.0000 ', 1.0, 0.0, 0.0), \
]
cat = {}
for d in data:
    cat[d[0]] = (ImageDesc.PHMS2RA(d[1]), ImageDesc.PDMS2Dec(d[2]), d[3], d[4],d[5])

