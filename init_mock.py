from add_ons.Make_mocks import *

###########################
# Set observation details #
###########################

direc = 'J2000 20h14m51.60s -24d30m23.678s'

doaca    = False
doalma   = True
ttaca    = ['5.2h','05h']
ttalma   = ['60.16min','020min']
almaconf = '1'
noisy    = True

synthatic_obs = Mockobs(direc, doaca, doalma, ttaca, ttalma, noisy, almaconf)

############################
# Initialize the SZ-params #
############################

popt = {'A10Pressure_0': {'RA':         303.7155549310805, 
                          'Dec':        -24.506378097715704, 
                          'log10':       14.751315458971625, 
                          'c500':        1.128, 
                          'e':           0.0, 
                          'Angle':       0.0, 
                          'Offset':      0.0, 
                          'Temperature': 0.0, 
                          'Alpha':       1.2223, 
                          'Beta':        5.4905, 
                          'Gamma':       0.7736, 
                          'P0':          3.412942589563412, 
                          'Alpha_p':    -0.1, 
                          'z':           0.1555, 
                          'bias':        0.0},
        
        'A10Pressure_1': {'RA':         303.7155549310805 + 20/60/60, 
                          'Dec':        -24.506378097715704 + 20/60/60, 
                          'log10':       14.751315458971625 - 0.5, 
                          'c500':        1.128, 
                          'e':           0.0, 
                          'Angle':       0.0, 
                          'Offset':      0.0, 
                          'Temperature': 0.0, 
                          'Alpha':       1.2223, 
                          'Beta':        5.4905, 
                          'Gamma':       0.7736, 
                          'P0':          3.412942589563412, 
                          'Alpha_p':    -0.1, 
                          'z':           0.1555, 
                          'bias':        0.0}
       }


SZmap = synthatic_obs.run(popt)