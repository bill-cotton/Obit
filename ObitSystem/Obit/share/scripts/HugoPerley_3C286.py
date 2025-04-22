# Hugo - Perley polarization model for 3C286 2024
#exec(open('HugoPerley_3C286.py').read())
HPAlias = ["3C286","J1331+3030"] # Aliases
import math
def HugoPerley_3C286 (nu):
    """
    Calculate polarization model 
    
    from SARAO memo SSA-0004E-001 B
    returns (P [fract lin. pol.], EVPA[deg])
    nu = frequency (Hz)
    """
    clight = 2.997924562e8   # Speed of light m/s
    lamb = (clight/nu); lamb2 = lamb**2   # Wavelength and squared
    nug = nu*1.0e-9          # GHz
    P = 10.0; EVPA = 33.0  # Default
    # EVPA By frequency range
    if (nug>=1.7) and (nug<=12.0):
        EVPA = 32.64 + 85.37*lamb2
    elif nug<1.7:
        c1 = math.log10(nug)**3
        EVPA = 29.53 + lamb2 *  (4005.88 *c1 - 39.38)
    # P By frequency range
    c2 =  math.log10(lamb2)
    if (nug>=1.1) and (nug<=12.0):
        P = 0.080 - 0.053*lamb2 - 0.015*c2
    elif nug<1.1:
        P = 0.029- 0.172*lamb2 - 0.067*c2
    return (P,EVPA)
# end HugoPerley_3C286 

def HugoPerley_3C286_l2 (lamb2):
    """
    Calculate polarization model 
    
    from SARAO memo SSA-0004E-001 B
    returns (P [fract. lin. pol.], EVPA[deg])
    lamb2 = wavelength^2 in m^2
    """
    lamb = lamb2**0.5
    clight = 2.997924562e8   # Speed of light m/s
    nu = clight/lamb
    nug = nu*1.0e-9          # GHz
    P = 10.0; EVPA = 33.0  # Default
    # EVPA By frequency range
    if (nug>=1.7) and (nug<=12.0):
        EVPA = 32.64 + 85.37*lamb2
    elif nug<1.7:
        c1 = math.log10(nug)**3
        EVPA = 29.53 + lamb2 *  (4005.88 *c1 - 39.38)
    # P By frequency range
    c2 =  math.log10(lamb2)
    if (nug>=1.1) and (nug<=12.0):
        P = 0.080 - 0.053*lamb2 - 0.015*c2
    elif nug<1.1:
        P = 0.029- 0.172*lamb2 - 0.067*c2
    return (P,EVPA)
# end HugoPerley_3C286_l2
