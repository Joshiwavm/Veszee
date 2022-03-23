from src.Modeler import *
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar
import numpy as np

plt.style.use('./plots/scripts/thesis.mplstyle')
matplotlib.rc('image', origin='lower')


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

def plot_it(image_names, telescp, filename_samples, imsize, imcell): 
    
    _, he = fits.getdata(image_names[0], header = True)
    FOV = [int(imsize/4 *1.6),int(imsize - imsize/4*1.6)]
    
    arcsec = 7.423904974016333
    pixel_width = int(arcsec/float(imcell.split('arcsec')[0]))
    x = np.arange(int(imsize*0.01), int(imsize*0.01)+pixel_width, 1)
    y = np.copy(x)
    y[:] = int(imsize*0.01)

    #observations
    dat0 = fits.getdata(image_names[0])[0,0]
    dat0_large = np.copy(dat0)
    dat0       = dat0[FOV[0]:FOV[1],FOV[0]:FOV[1]] 

    dat1 = fits.getdata(image_names[1])[0,0]
    dat1 = dat1[FOV[0]:FOV[1],FOV[0]:FOV[1]]

    dat2 = fits.getdata(image_names[2])[0,0]
    dat2_large = np.copy(dat2)
    dat2 = dat2[FOV[0]:FOV[1],FOV[0]:FOV[1]]

    vmin = np.nanstd(dat0) * -4 * 1e3 #THINK ABOUT COLOR SCALE
    vmax = np.nanmax(dat0)* 1e3

    #contours
    mask = circle_mask(dat0_large, imsize//2 ,imsize//2, 0.12*imsize)
    
    rms_2 = np.nanstd(dat2_large[~mask])
    cont_2 = np.array([-8, -7, -6,-5,-4,-3,-2, 2, 3, 4])*rms_2

    #Figure
    wpanel = 20
    wspace1 = -9.5
    wspace2 = -4.0
    wspace3 = -1.5
    wbar = 0.8

    fig = plt.figure(figsize = (15,4))

    gs = GridSpec(1,9, wspace=0 , width_ratios=[wpanel,wspace1,wpanel,wspace2,wbar,wspace3,wpanel,wspace2,wbar], height_ratios = [1])

    if telescp == 'com07m':
        fig.suptitle('ACA Observations')    
    elif telescp == 'com12m':
        fig.suptitle('ALMA Observations')
        
    #first panel == Observation
    ax0 = fig.add_subplot(gs[0])
    im0 = ax0.imshow(dat0*1e3, vmin = vmin, vmax = vmax)  

    ax0.text(int(len(dat0)*0.05), int(len(dat0)*0.9), 'Observation', c = 'white')    
    ax0.set_xticklabels([])
    ax0.set_xticks([]) 
    ax0.set_yticklabels([])
    ax0.set_yticks([]) 

    ax0.plot(x, y, '-', linewidth=3, color='white')
    ax0.text(int(imsize*0.01),int(imsize*0.02), r'20 kpc', verticalalignment='top', color='white', size=11)

    #second panel == Model
    ax1 = fig.add_subplot(gs[2])
    ax1.imshow(dat1*1e3, vmin = vmin, vmax = vmax)    

    ax1.text(int(len(dat1)*0.05), int(len(dat1)*0.9), 'Model', c = 'white')    
    ax1.set_xticklabels([])
    ax1.set_xticks([]) 
    ax1.set_yticklabels([])
    ax1.set_yticks([]) 
    
    ax1.plot(x, y, '-', linewidth=3, color='white')
    ax1.text(int(imsize*0.01),int(imsize*0.02), r'20 kpc', verticalalignment='top', color='white', size=11)
    
    #first bar = colorbar
    ax2 = fig.add_subplot(gs[4])
    cbar1 = fig.colorbar(im0, cax=ax2, aspect = 20)
    cbar1.set_label(r'S$_\nu$ [mJy/beam]')

    #third panel == Residue
    ax3 = fig.add_subplot(gs[6]) 
    im3 = ax3.imshow(dat2*1e3)    

    ax3.text(int(len(dat2)*0.05), int(len(dat2)*0.9), 'Residuals', c = 'white' )    
    ax3.set_xticklabels([])
    ax3.set_xticks([]) 
    ax3.set_yticklabels([])
    ax3.set_yticks([]) 
 
    ax3.plot(x, y, '-', linewidth=3, color='white')
    ax3.text(int(imsize*0.01),int(imsize*0.02), r'20 kpc', verticalalignment='top', color='white', size=11)

    ax3.contour(dat2*1e3, cont_2*1e3, colors='k', linestyles = 'solid', linewidths = 0.5)

    #second bar = colorbar
    ax4 = fig.add_subplot(gs[8])
    cbar2 = fig.colorbar(im3, cax=ax4, aspect = 20)
    cbar2.set_label(r'S$_\nu$ [mJy/beam]')

    plt.tight_layout()
    plt.savefig('./plots/outputs/SZimage_plots_'+(filename_samples.split('/')[-1]).split('_pickle')[0]+ '_'+telescp+'.pdf', bbox_inches='tight')
    plt.close()

def clean_it(ms_files, imsize, imcell):
    
    image_names = []
    for j in range(len(ms_files)):

        tclean( vis         =                         ms_files[j],
                imagename   =    ms_files[j].replace('.ms','.im'),
                niter       =                                   0,
                imsize      =                              imsize,
                cell        =                              imcell,
                gridder     =                          'standard',
                weighting   =                           'natural',
                specmode    =                               'mfs',
                parallel    =                               False)

        exportfits(imagename=ms_files[j].replace('.ms','.im.image'), 
                   fitsimage=ms_files[j].replace('.ms','.im.fits'), 
                   overwrite=True)
        
        image_names.append(ms_files[j].replace('.ms','.im.fits'))
        
        os.system('rm -rf {0}.pb'.format(ms_files[j].replace('.ms','.im')))
        os.system('rm -rf {0}.psf'.format(ms_files[j].replace('.ms','.im')))

        os.system('rm -rf {0}.image'.format(ms_files[j].replace('.ms','.im')))
        os.system('rm -rf {0}.model'.format(ms_files[j].replace('.ms','.im')))
        os.system('rm -rf {0}.sumwt'.format(ms_files[j].replace('.ms','.im')))
        os.system('rm -rf {0}.weight'.format(ms_files[j].replace('.ms','.im')))
        os.system('rm -rf {0}.residual'.format(ms_files[j].replace('.ms','.im')))

    os.system('rm -vf *.last')
    os.system('rm -vf *.log')
    
    return image_names 
    
def plot_it_im(filename_samples, filename_ms, telescp):

    print()
    print("Creating image-plane figures")

    for i in range(len(telescp)):
        filedir_save = './output/{1}/model/'.format('3',telescp[i])
        ms_files = [
            filename_ms[i], 
            filedir_save + 'Model_'+ (filename_samples.split('/')[-1]).split('_pickle')[-2] + '.ms',
            filedir_save + 'Model-residue_'+ (filename_samples.split('/')[-1]).split('_pickle')[-2] + '.ms'
        ]

        if telescp[i] =='com07m':
            imsize = 256
            imcell = '1.50arcsec'
        elif telescp[i] =='com12m':
            imsize = 1024
            imcell = '0.15arcsec'

        image_names = clean_it(ms_files, imsize, imcell)
        plot_it(image_names, telescp[i], filename_samples, imsize, imcell)