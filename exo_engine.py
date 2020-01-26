import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table
import batman
from ipywidgets import interactive, fixed
from IPython.display import display

plt.rcParams.update({'font.size': 18})

t = np.linspace(-1,1,100) ## time

def initial_lightcurve(radius=0.1):
    params = batman.TransitParams()
    params.t0 = 0.                       #time of inferior conjunction
    params.per = 48.                     #orbital period
    params.rp = radius                   #planet radius (in units of stellar radii)
    params.a = 15.                       #semi-major axis (in units of stellar radii)
    params.inc = 87.                     #orbital inclination (in degrees)
    params.ecc = 0.                      #eccentricity
    params.w = 90.                       #longitude of periastron (in degrees)
    params.u = [0.2, 0.0]                #limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"       #limb darkening model
    
    return params


def plot_lc(radius=0.1,color='blue'):
    params = initial_lightcurve(radius=radius)
    t = np.linspace(-1, 1, 100)
    m = batman.TransitModel(params, t)    #initializes model
    flux = m.light_curve(params)          #calculates light curve
    plt.plot(t, flux)
    plt.xlabel("Time from central transit")
    plt.ylabel("Relative brightness")
    plt.ylim(0.985,1.002)
    return radius

def plot_initial_lc():
    plot_lc(radius=0.1)

def plot_interactive_rad():
    lc_rad = interactive(plot_lc, radius=(0.0,0.12,0.005),color=fixed('blue'))
    #rad = interact(plot_lc, radius=(0.0,0.12,0.005))
    display(lc_rad)
    return lc_rad


slopeStart = -0.3
slopeEnd = 0.0

class spectral_lc:
    def __init__(self):
        """
        Class for spectroscopic lightcuves
        """
        self.wavelengths =  np.array([  6.0   ,   5.0   ,  4.0   ,   3.0   , 2.0     , 1.0 ])
        self.calc_radii()
        #self.radius_array = np.array([  0.08  ,  0.085  , 0.090  , 0.095   ,0.100   ,  0.105])
        self.colors_array = np.array([  'red' ,'orange','yellow' ,'green',  'blue',  'violet'])
        #self.limb_dark_arr = [   0.1 , 0.15   ,   0.2,   0.25     , 0.3     , 0.35]

    
    def calc_radii(self,slope=-0.3):
        self.radius_array = 0.1 + slope * (self.wavelengths - 2.5) / 50.
    
    def plot_lc_multicolor_loop(self,slope=slopeStart):
        self.calc_radii(slope=slope)
        for one_radius,one_color in zip(self.radius_array,self.colors_array):
            plot_lc(one_radius,color=one_color)
    
    def plot_lc_multicolor(self,slope=-0.3):
        lc_multi = interactive(self.plot_lc_multicolor_loop, slope=(slopeStart,slopeEnd,0.02))
        display(lc_multi)
    
    def spectrum_plot(self,slope=slopeStart):
        self.calc_radii(slope=slope)
        plt.plot(self.wavelengths,self.radius_array * 10.)
        plt.xlabel('Wavelength (microns)')
        plt.ylabel('Size (Jupiter Radii)')
        plt.ylim(0.8,1.1)
        
    def spectrum_plot_i(self):
        spec_p = interactive(self.spectrum_plot,slope=(slopeStart,slopeEnd,0.02))
        display(spec_p)
        
    def visualize_colors(self,slope=slopeStart):
        self.calc_radii(slope=slope)
        
        fig, ax = plt.subplots(figsize=(8,8))
        ax.set_aspect('equal')
        ax.set_xlim(-1.3,1.3)
        ax.set_ylim(-1.3,1.3)
        for one_radius,one_color in zip(self.radius_array,self.colors_array):
            circlePatch = plt.Circle((0., 0.), one_radius * 10.,linewidth=2,
                                     edgecolor=one_color,facecolor='none')
            ax.add_artist(circlePatch)

        for ind,angle in zip([0,-1],[0,1,0.7]):
            ax.plot([0,self.radius_array[ind] * np.cos(angle) * 10.],
                    [0,self.radius_array[ind] * np.sin(angle) * 10.],color=self.colors_array[ind])
            ax.text(0.2,np.sin(angle)* 0.5,"{} microns".format(self.wavelengths[ind]))
        ax.set_xlabel('X Size in Jupiters')
        ax.set_ylabel('Y Size in Jupiters')
    
    def visualize_colors_i(self):
        size_p = interactive(self.visualize_colors,slope=(slopeStart,slopeEnd,0.02))
        display(size_p)
        
    