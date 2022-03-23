from add_ons.Jack_Knife import *

vis                 = './data/RXC_J2014.8_ACA_timebin30s_SZ_concat.ms'
typ                 = 'com07m'
model_name          = './output/com07m/model/Model_RXC_J2014.8_ACA_and_ALMA_cc_pointandgaus.im.fits'

imager              = Test_SZonly(vis, typ, model_name)
Noise, SZ, Combined = imager.run()
