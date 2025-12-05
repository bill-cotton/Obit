# Cosine square beam shape function
from math import radians, cos, pi
Cos2Beam = None
del Cos2Beam
def Cos2Beam(rho, nu, D, fudge=1.052237):
    """
    Calculate cosine sq. beam shape (Condon & Ransom, Essential Radio Astronomy eq 3.94)
    
    Return power gain of circularly symmetric beam
    * rho   = offset from center (degrees)
    * nu    = Frequency (Hz)
    * D     = Antenna diameter (m)
    * fudge  = Fudge factor for rho*D/lambda, default=MeerKAT
    """
    ################################################################
    c = 2.997924562e8;  # Velocity of light m/s
    arg = fudge*radians(rho)*D*nu/c
    div = (1.-4.*(arg**2))
    if abs(div)<1.0e-5:
        div = 1.0e-5
    gain = (cos(pi*arg)/div)**2
    return gain
# end MKCosBeam
