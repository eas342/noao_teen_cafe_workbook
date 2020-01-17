import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table
import batman

t = np.linspace(-1,1,100) ## time

def initial_lightcurve():
    params = batman.TransitParams()
    params.t0 = 0.                       #time of inferior conjunction
    params.per = 48.                      #orbital period
    params.rp = 0.1                      #planet radius (in units of stellar radii)
    params.a = 15.                       #semi-major axis (in units of stellar radii)
    params.inc = 87.                     #orbital inclination (in degrees)
    params.ecc = 0.                      #eccentricity
    params.w = 90.                       #longitude of periastron (in degrees)
    params.u = [0.3, 0.0]                #limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"       #limb darkening model
    
    return params
    
