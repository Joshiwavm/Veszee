import scipy
from scipy.optimize import curve_fit
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse  
plt.style.use('./plots/scripts/thesis.mplstyle')
matplotlib.rc('image', origin='lower')

from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
from astropy import units as u
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

import casatools
from casatasks import *

####### Own Tools #######
from src.MsManager import *
from src.Settings import *
from src.UnitTransform import *

Tcmb    = 2.7255
mec2    = ((const.m_e*const.c*const.c).to(u.keV)).value
clight  = const.c.value
kboltz  = const.k_B.value
hplanck = const.h.value

def circle_mask(im, xc, yc, rcirc):
        ny, nx = im.shape
        y,x = np.mgrid[0:nx,0:ny]
        r = np.sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc))
        return ( (r < rcirc))

def r_theta(im, xc, yc):
    # returns the radius rr and the angle phi for point (xc,yc)
    ny, nx = im.shape
    yp, xp = np.mgrid[0:ny,0:nx]
    yp = yp - yc
    xp = xp - xc
    rr = np.sqrt(np.power(yp,2.) + np.power(xp,2.))
    phi = np.arctan2(yp, xp)
    return(rr, phi)



class Mockobs:
    def __init__(self, direc, doaca, doalma, ttaca, ttalma, noisy, almaconf = '1'):
        self.doaca    = doaca
        self.doalma   = doalma
        self.ttaca    = ttaca
        self.ttalma   = ttalma
        self.noisy    = noisy
        self.almaconf = almaconf
        
        if self.doaca:
            self.imsize = 256
            self.imcell = '1.50arcsec'
        elif self.doalma:
            self.imsize = 1024
            self.imcell = '0.15arcsec'

        self.direc = direc
                    
    def _make_grid(self, popt):

        typ, RA, Dec = self.direc.split(' ')
        c = SkyCoord(RA, Dec, frame='icrs') #don't make it direc but from popt

        x, y  = np.meshgrid(np.arange(-self.imsize/2, self.imsize/2, 1.), np.arange(-self.imsize/2, self.imsize/2, 1.))
        x += 0.5
        y -= 0.5
                
        RA  = x * float(self.imcell.split('arcsec')[0])/60/60 + c.ra.value
        Dec = y * float(self.imcell.split('arcsec')[0])/60/60 + c.dec.value

        r  = (Dec - popt['Dec'])**2 
        r += ((RA - popt['RA'])*np.cos(np.deg2rad(c.dec.value)))**2
        r  = r**0.5 
        
        return r

    def _make_modelimage(self, popt):

        image = np.zeros((self.imsize, self.imsize))

        for t  in popt.keys(): 
            if COMPONENTS[t.split('_')[0]]['make_image'] == True:     
                
                grid  = self._make_grid(popt[t])

                transform = TransformInput(popt[t], t)
                input_par = transform.generate() 

                info = getinfo()
                
                coord = grid * info.cosmo.angular_diameter_distance(popt[t]['z'])
                coord = coord.value/input_par['major']

                if t.split('_')[0] == 'gnfwPressure':
                    rs = info.linmesh[1:]
                else:
                    rs = info.linmesh

                integrated_P = COMPONENTS[t.split('_')[0]]['function'](rs, **input_par)
                image += np.interp(coord, rs, integrated_P) * info.ysznorm.value
        return image 

    def _make_model(self, popt):
        Ymap      = self._make_modelimage(popt)
        SZmap     = Ymap * comptonToJyPix(95e9,                                        #freq
                                          float(self.imcell.split('arcsec')[0])/60/60, #pixelsize
                                          float(self.imcell.split('arcsec')[0])/60/60) #pixelsize

        if self.doalma:
            fits.writeto("./output/add_ons/synthatic_SZmodel_com12m.fits", SZmap, overwrite = True)
        elif self.doaca:
            fits.writeto("./output/add_ons/synthatic_SZmodel_com07m.fits", SZmap, overwrite = True)
            
    def _make_mocks(self):
        file = open('ptgfile.txt','w')
        file.write(self.direc)
        file.close()

        names = []
        if self.doaca:
            simobserve(project = 'test_{0}_{1}'.format(self.ttaca[1] if self.doaca else 'noaca', self.ttalma[1] if self.doalma else 'noalma'),
                      skymodel = './output/add_ons/synthatic_SZmodel_com07m.fits',
                  setpointings = False,
                       ptgfile = 'ptgfile.txt',
                     overwrite = True,
                     totaltime = self.ttaca[0],
                      inbright = '',
                        incell = '1.50arcsec',
                   indirection = self.direc,
                      incenter = '95GHz',
                       inwidth =  '8GHz',
                   antennalist = './add_ons/confc8/aca.cycle7.cfg',
                      graphics = 'none')

        if self.doalma:
            simobserve(project = 'test_{0}_{1}'.format(self.ttaca[1] if self.doaca else 'noaca',self.ttalma[1] if self.doalma else 'noalma'),
                      skymodel = './output/add_ons/synthatic_SZmodel_com12m.fits',
                  setpointings = False,
                       ptgfile = 'ptgfile.txt',
                     overwrite = True,
                     totaltime = self.ttalma[0],
                      inbright = '',
                        incell = '0.15arcsec',
                   indirection = self.direc,
                      incenter = '95GHz',
                       inwidth =  '8GHz',
                   antennalist = './add_ons/confc8/alma.cycle7.{0}.cfg'.format(self.almaconf),
                      graphics = 'none')
        
        if os.path.exists('./output/add_ons/'+'test_{0}_{1}'.format(self.ttaca[1] if self.doaca else 'noaca',self.ttalma[1] if self.doalma else 'noalma')):
            shutil.rmtree('./output/add_ons/'+'test_{0}_{1}'.format(self.ttaca[1] if self.doaca else 'noaca',self.ttalma[1] if self.doalma else 'noalma'))
            
        shutil.move('test_{0}_{1}'.format(self.ttaca[1] if self.doaca else 'noaca',self.ttalma[1] if self.doalma else 'noalma'), './output/add_ons/')
        os.system('rm -vf ptgfile.txt')

    def _image_mocks(self):
        name = 'test_{0}_{1}'.format(self.ttaca[1] if self.doaca else 'noaca',self.ttalma[1] if self.doalma else 'noalma')
        
        if self.noisy:
            if self.doalma:
                name = './output/add_ons/' + name + '/' + name + '.alma.cycle7.'+self.almaconf+'.noisy.ms'
            elif self.doaca:
                name = './output/add_ons/' + name + '/' + name + '.aca.cycle7.noisy.ms'
        else:
            if self.doalma:
                name = './output/add_ons/' + name + '/' + name + '.alma.cycle7.'+self.almaconf+'.ms'
            elif self.doaca:
                name = './output/add_ons/' + name + '/' + name + '.aca.cycle7.ms'

        tclean( vis         =                                name,
                imagename   =          name.replace('.ms', '.im'), 
                niter       =                                   0,
                imsize      =                         self.imsize,
                cell        =                         self.imcell,
                gridder     =                          'standard',
                weighting   =                           'natural',
                specmode    =                               'mfs',
                parallel    =                               False)

        exportfits(imagename=name.replace('.ms', '.im.image'),
                   fitsimage=name.replace('.ms', '.im.fits'),
                   overwrite=True)

        image_names_synthatic = name.replace('.ms', '.im.fits')
        return image_names_synthatic

    def _show_images(self, image_names_synthatic):
        Synthatic_SZimage, header_SZ        = fits.getdata(image_names_synthatic, header=True)
        Synthatic_SZimage = Synthatic_SZimage[0,0] *1e3

        FOV     = [int(self.imsize/4 *1.3),int(self.imsize - self.imsize/4*1.3)]
        mask_rms = circle_mask(Synthatic_SZimage, self.imsize//2, self.imsize//2, 0.2*self.imsize)

        fig, ax = plt.subplots(1,1)
        fig.suptitle("Synthatic SZ Observations")

        im0 = ax.imshow(Synthatic_SZimage[FOV[0]:FOV[1],FOV[0]:FOV[1]])
        cbar = plt.colorbar(im0, ax = ax)
        cbar.set_label(r'$S_{\nu}$ [mJy/Beam]')
        ax.text(0.05*(FOV[1]-FOV[0]),0.90*(FOV[1]-FOV[0]), 
                '{} = {:.4f} [mJy/Beam]'.format(r"$\sigma_{S_\nu}$", np.nanstd(Synthatic_SZimage[~mask_rms])), 
                color='white')
            
        if self.doalma:
            x = np.arange(0.05*(FOV[1]-FOV[0]), 0.05*(FOV[1]-FOV[0])+60, 1)
        elif  self.doaca:
            x = np.arange(0.05*(FOV[1]-FOV[0]), 0.05*(FOV[1]-FOV[0])+6, 1)
            
        y = np.copy(x)
        y[:] = 0.05*(FOV[1]-FOV[0])

        ax.plot(x, y, '-', linewidth=3, color='white')
        ax.text(0.07*(FOV[1]-FOV[0]), 0.10*(FOV[1]-FOV[0]), r'9"', 
                verticalalignment='top', color='white', size=11)

        ax.set_xticks([])
        ax.set_yticks([])
            
        self.modelcont = np.nanstd(Synthatic_SZimage[~mask_rms])*np.array([-7,-6, -5, -4, -3, -2, 2, 3, 4])
        ax.contour(Synthatic_SZimage[FOV[0]:FOV[1],FOV[0]:FOV[1]], self.modelcont)

        plt.savefig('./plots/add_ons/Synthatic_SZobs_{0}_{1}_noisy{2}.pdf'.format(self.ttaca[1] if self.doaca else 'noaca', self.ttalma[1] if self.doalma else 'noalma', str(self.noisy)))
        plt.close()
        
        return Synthatic_SZimage
        
    def run(self, popt):
        
        print('.. Initialize Model')
        self._make_model(popt)
        
        print()
        print('.. Make Synthatic Observations')
        self._make_mocks()

        print()
        print('.. Image mocks')
        names = self._image_mocks()
        
        print()
        print('.. Visualize')
        SZmap = self._show_images(names)
        
        os.system('rm -vf *.last')
        os.system('rm -vf *.log')
        
        return SZmap