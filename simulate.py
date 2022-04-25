from add_ons.Make_mocks import *

from astropy.wcs import WCS
from reproject import reproject_interp

def mock_it(t, RA, Dec, Mass, z):
    """
    Input: 
        - t (float): on source time
        - RA (float): RA in degree [J200]
        - Dec (float): Dec in degree [J200]
        - Mass (float): mass in terms in x10^(14) M_\odot
        - z (float): redshift
    """
    
    c = SkyCoord(RA, Dec, frame='icrs', unit='deg')
    direc = 'J2000 {}'.format(c.to_string('hmsdms'))
    
    ttaca    = ['{}h'.format(t),'{}h'.format(int(t))]
    ttalma   = ['60.16min','020min']
    doaca    = True
    doalma   = False
    
    synthatic_obs = Mockobs(direc, doaca, doalma, ttaca, ttalma, '1')
    
    # change this dictionairy according to MD, A10-UPP, or CC.
    popt = {'A10Pressure_0': {'RA':          RA, 
                              'Dec':         Dec, 
                              'log10':       np.log10(Mass*1e14), 
                              'c500':        1.1770E+00, 
                              'e':           0.0, 
                              'Angle':       0.0, 
                              'Offset':      0.0, 
                              'Temperature': 0.0, 
                              'Alpha':       1.0510E+00, 
                              'Beta':        5.4905E+00, 
                              'Gamma':       3.0810E-01, 
                              'P0':          8.4030E+00, 
                              'Alpha_p':     1.2000E-01, 
                              'z':           z, 
                              'bias':        0.0}
           }
    
    SZmap_noisy, SZmap_noisless, header = synthatic_obs.run(popt)
    std  = np.nanstd(SZmap_noisy - SZmap_noisless)    
    return SZmap_noisy, SZmap_noisless, std, header



RA   = ...
Dec  = ...
M500 = ...
z    = ...


SZmap_noisy, SZmap_noisless, std, header = mock_it(2.1, RA, Dec, M500, z)
SNR_want                                 = abs(np.nanmin(SZmap_noisless)/93.71266889746519e-3)
t_integration                            = (8/SNR_want)**2 * 2.1    

print("Found Integration time: t = {:.1f}h, SNR = {:.2f}".format(t_integration, 8))
print()
