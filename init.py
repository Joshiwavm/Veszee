from plots.scripts.plot_Im import *
from plots.scripts.plot_UV import *

##add on


telescp      = ['com07m','com12m']
filename_samples   = './data/RXC_J2014.8_ACA_and_ALMA_cc_pointandgaus_pickle'

filename_ms  = [
    './data/RXC_J2014.8_ACA_timebin30s_SZ_concat.ms',
    './data/RXC_J2014.8_ALMA_timebin30s_SZ_concat.ms'
]

plot_it_uv(filename_samples, filename_ms, telescp)
plot_it_im(filename_samples, filename_ms, telescp)