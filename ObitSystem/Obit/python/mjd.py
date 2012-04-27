'''
Created on Jan 29, 2010

@author: jonathan.hirschauer

Class for converting dates between Gregorian and MJD formats. [Module
obtained from Google code repository for project python-plotter-qwt.]
'''

import time
from math import floor

class mjd:
    '''
    Functions for converting to/from MJD
    '''
    def __init__(self):
        '''
        Constructor
        '''
        pass
        
    def mjd_now(self, m=None, d=None, y=None, h=None, mi=None, s=None):
        '''
        Returns floating point MJD to the nearest second
        if m==None then returns the current UTC time in MJD
        if m=month, d=day, y=year, h=hour, mi=minute, s=seconds
        '''
        if m==None:
            t_time = time.gmtime()
            y = t_time[0]   # Year
            m = t_time[1]   # Month
            d = t_time[2]   # Day
            h = t_time[3]   # Hour
            mi = t_time[4]   # Minute
            s = t_time[5]   # Second

        mjd_frac = float(0.0)
        mjd_day = 0.0
        mjd_frac += (h*3600+mi*60+s)/float(86400)
        #print mjd_frac
        mjd_day += 367*y - 7*(y+(m+9)/12)/4 \
                         - 3*((y+(m-9)/7)/100+1)/4 \
                         + 275*m/9 \
                         + d + 1721028 - 2400000 
        return mjd_day + mjd_frac

    def mjd_to_datetime(self, mjd):
        '''
        returns a tuple (y,mo,d,h,mi,s)
        '''
        jd = floor(mjd)+2400000.5
        jdi = floor(jd)
        jdf = jd-jdi+0.5
        
        if jdf >= 1.0:
            jdf = jdf-1.0
            jdi += 1
        
        l = jdi + 68569
        n = floor(4*l/146097)
        
        l = floor(l) - floor((146097*n+3)/4)
        year = floor(4000*(l+1)/1461001)
       
        l = l-floor(1461*year/4) + 31
        month = floor(80*l/2447)
        day = l-floor(2447*month/80)
        l = month/11
        month = floor(month+2-(12*l))
        year = floor(100*(n-49)+year+l)

        mjdseconds = (mjd-floor(mjd))*86400
        hour = floor(mjdseconds/3600)
        minute = floor((mjdseconds-hour*3600)/60)
        second = floor(mjdseconds-hour*3600-minute*60)
        return (year, month, day, hour, minute, second)


if __name__ == "__main__":
    t = mjd()
    print t.mjd_now(9,16,2010,20,03,30)
    print t.mjd_now()
    print t.mjd_to_datetime(t.mjd_now())
    print time.asctime(time.gmtime())
    
    
