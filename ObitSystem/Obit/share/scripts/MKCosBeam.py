# MeerKAT Cosine beam shape function
from math import radians, cos, pi
MKCosBeam = None
del MKCosBeam
def MKCosBeam(rho, nu):
    """
    Calculate cosine beam shape (Condon & Ransom, Essential Radio Astronomy eq 3.95)
    
    Return power gain of circularly symmetric beam
    * rho   = offset from center (degrees)
    * nu    = Frequency (Hz)
    * D     = Antenna diameter (m)
    * beam  = beam FWHM (amin)
    """
    ################################################################
    #theta_b = radians(57.5/60) * (1.5e9/nu)
    theta_b = 0.0167261 * (1.5e9/nu)
    rhor = 1.18896*radians(rho)/theta_b
    div = (1.-4.*(rhor**2))
    if abs(div)<1.0e-5:
        div = 1.0e-5
    gain = (cos(pi*rhor)/div)**2
    return gain
# end MKCosBeam
