from src.Modeler import *
import matplotlib.pyplot as plt
import numpy as np

import matplotlib

plt.style.use('./plots/scripts/thesis.mplstyle')
matplotlib.rc('image', origin='lower')

def uv_load(name, telescp, i):

    MS = MsReader(name, telescp[i])
    uvdist, uv_real, uv_imag, uv_wghts = MS.uvdata_loader()
        
    return uv_real.flatten(), uv_imag.flatten(), uv_wghts.flatten(), uvdist.flatten()

def bin_it(UVreal, UVimag, uvwghts, uvdist, nbins):
    bins = np.logspace(np.log10(np.nanmin(uvdist)), np.log10(np.nanmax(uvdist)), nbins+1)

    UVrealbinned = np.empty(nbins)
    UVrealerrors = np.empty(nbins)
    UVimagbinned = np.empty(nbins)
    UVimagerrors = np.empty(nbins)
    UVabsbinned  = np.empty(nbins)
    UVabserrors  = np.empty(nbins)
    
    for i in range(nbins):
        mask = (uvdist>bins[i]) & (uvdist <= bins[i+1])
        UVrealbinned[i], UVrealerrors[i] = np.average(UVreal[mask], weights = uvwghts[mask], returned = True) 
        UVimagbinned[i], UVimagerrors[i] = np.average(UVimag[mask], weights = uvwghts[mask], returned = True)
        UVabsbinned[i],  UVabserrors[i]  = np.average((UVreal[mask]**2+UVimag[mask]**2)**0.5, weights = uvwghts[mask], returned = True) 

        UVrealerrors[i] = UVrealerrors[i]**(-0.5)
        UVimagerrors[i] = UVimagerrors[i]**(-0.5)
        UVabserrors[i]  = UVabserrors[i]**(-0.5)
        
    return UVrealbinned, UVrealerrors, UVimagbinned, UVimagerrors, UVabsbinned, UVabserrors, bins
    
def plot_it_uv(filename_samples, filename_ms, telescp): 

    for j, res in enumerate([False, True]):
        print()
        print("Creating UV-plot, Residue Image? " + str(res))
        
        uvdist = np.empty(0)
        uvreal = np.empty(0)
        uvimag = np.empty(0)
        uvwght = np.empty(0)
        uvtype = np.empty(0, dtype = int)
        us     = np.empty(0)
        vs     = np.empty(0)
        
        for i in range(len(telescp)):
        
            filedir_save = './output/{1}/model/'.format('3',telescp[i])
            
            ms_files = [
            filedir_save + 'Model_'+ (filename_samples.split('/')[-1]).split('_pickle')[-2] + '.ms',
            filedir_save + 'Model-residue_'+ (filename_samples.split('/')[-1]).split('_pickle')[-2] + '.ms'
            ]

            if j == 0:
                obj = Modeler(filename_samples, filename_ms[i], telescp[i], obs_amount = len(telescp), save = True, load_pb = False)
                _, uvdata = obj.run()
                u, v, freqs = uvdata    
            else:
                pass
                
            r_uv_real, r_uv_imag, uv_wghts, uv_dist = uv_load(ms_files[j], telescp, i)
            r_uv_real *= 1e3
            r_uv_imag *= 1e3

            uvtype = np.append(uvtype, np.ones_like(u, dtype = int)*i)
            uvdist = np.append(uvdist, uv_dist)
            uvreal = np.append(uvreal, r_uv_real)
            uvimag = np.append(uvimag, r_uv_imag)
            uvwght = np.append(uvwght, uv_wghts)
            us     = np.append(us, u)
            vs     = np.append(vs, v)
    
        UVrealbinned, UVrealerrors, UVimagbinned, UVimagerrors, UVabsbinned, UVabserrors, bins = bin_it(uvreal, uvimag, uvwght, uvdist, nbins = 12)        
        plot_it(UVrealbinned, UVrealerrors, UVimagbinned, UVimagerrors, UVabsbinned, UVabserrors, bins, res, telescp, filename_samples,                uvdist, uvreal, uvimag, uvwght, uvtype, us, vs)
    
    
def plot_it(UVrealbinned, UVrealerrors, UVimagbinned, UVimagerrors, UVabsbinned, UVabserrors, bins, residue, telescp, filename_samples, uvdist, uvreal, uvimag, uvwght, uvtype, us, vs):

    select = np.random.uniform(low = 0, high = len(uvdist), size = 200000).astype(int)
    # select = np.arange(0, len(uvdist), 1).astype(int)
        
    fig, ax = plt.subplots(1)
    fig.set_figwidth(6)
    fig.set_figheight(4)
    
    # uvdist - Re(V)
    ax.errorbar((bins[1:]+bins[:-1])/2, UVrealbinned, xerr = (bins[1:] - bins[0:-1])/2, yerr = UVrealerrors, c = 'C0', ls = '', marker = 'o')
    ax.axhline(0, ls = 'dashed', c = 'gray')
    
    if residue == False:   
        np.save('/scigarfs/home/jvmarrewijk/eszee/UVspace_operations/Binned_SZall.npy', np.array([UVrealbinned, UVimagbinned, bins]), allow_pickle=True)
        UVrealbinnedSZ, UVimagbinnedSZ,binsSZ = np.load('/scigarfs/home/jvmarrewijk/eszee/UVspace_operations/Binned_SZall.npy',  allow_pickle=True)
        ax.plot((binsSZ[1:] + binsSZ[0:-1])/2,UVrealbinnedSZ, c ='C1', label = 'SZ-Profile')
    
    ax.axis(ymin = np.nanmin(UVrealbinned)-0.5, ymax = np.nanmax(UVrealbinned)+0.5)        
    ax.set_ylabel('Fourier Amplitude [mJy]')
    ax.set_xlabel('uv distance [k$\lambda$]')
    ax.set_xscale('log')
    plt.legend()
    
    plt.tight_layout()
    if residue == False:
        plt.savefig('./plots/outputs/SZuvflux_model_'+(filename_samples.split('/')[-1]).split('_pickle')[0]+ '.png', dpi = 300)
    else:
        plt.savefig('./plots/outputs/SZuvflux_residue_'+(filename_samples.split('/')[-1]).split('_pickle')[0]+ '.png', dpi = 300)
    plt.close()
    
    
###################################################################################################    
    
    fig, ax = plt.subplots(2,2)
    fig.set_figwidth(10)
    fig.set_figheight(8)
    
    # uvdist - Re(V)
    if len(np.unique(uvtype))>1:
        ax[0,0].scatter(uvdist[select][uvtype[select] == 0], uvreal[select][uvtype[select] == 0], s = .1, label = telescp[0])
        ax[0,0].scatter(uvdist[select][uvtype[select] == 1], uvreal[select][uvtype[select] == 1], s = .1, label = telescp[1])
    else:
        ax[0,0].scatter(uvdist[select][uvtype[select] == 0], uvreal[select][uvtype[select] == 0], s = .1, label = telescp[0])

    ax[0,0].errorbar((bins[1:]+bins[:-1])/2, UVrealbinned, xerr = (bins[1:] - bins[0:-1])/2, yerr = UVrealerrors, c = 'C2', ls = '', marker = 'o')
    ax[0,0].axhline(0, ls = 'dashed', c = 'gray')    
    ax[0,0].axis(ymin = np.nanmin(UVrealbinned)-0.5, ymax = np.nanmax(UVrealbinned)+0.5)        

    if residue == False:
        ax[0,0].set_ylabel('Re(v)$_{model}$ [mJy]')
    else:
        ax[0,0].set_ylabel('Re(v)$_{obs}$ - Re(v)$_{model}$ [mJy]')
        
    ax[0,0].set_xlabel('uv distance [k$\lambda$]')
    ax[0,0].set_xscale('log')
    ax[0,0].legend(loc=3)

    # uvdist - phase(V)
    if len(np.unique(uvtype))>1:
        ax[0,1].scatter(uvdist[select][uvtype[select] == 0], 
                        np.angle(uvreal[select][uvtype[select] == 0] + uvimag[select][uvtype[select] == 0]*1j, deg=True), 
                        s = .1, label = telescp[0])
        ax[0,1].scatter(uvdist[select][uvtype[select] == 1], 
                        np.angle(uvreal[select][uvtype[select] == 1] + uvimag[select][uvtype[select] == 1]*1j, deg=True), 
                        s = .1, label = telescp[1])
    else:
        ax[0,1].scatter(uvdist[select][uvtype[select] == 0], 
                        np.angle(uvreal[select][uvtype[select] == 0] + uvimag[select][uvtype[select] == 0]*1j, deg=True), 
                        s = .1, label = telescp[0])

    ax[0,1].axhline(0, ls = 'dashed', c = 'gray')
    ax[0,1].axis(xmin = 1e0)
    
    if residue == False:
        ax[0,1].set_ylabel(r'$\phi$(V)$_{model}$ [deg]')
    else:
        ax[0,1].set_ylabel(r'$\phi$(V)$_{obs}$ - $\phi$(V)$_{model}$ [deg]')
        
    ax[0,1].set_xlabel('uv distance [k$\lambda$]') 
    ax[0,1].set_xscale('log')
    ax[0,1].legend(loc=3)
    
    # uvdist - |V|
    if len(np.unique(uvtype))>1:
        ax[1,0].scatter(uvdist[select][uvtype[select] == 0], 
                        (uvreal[select][uvtype[select] == 0]**2 + uvimag[select][uvtype[select] == 0]**2)**0.5, 
                        s = .1, label = telescp[0])
        ax[1,0].scatter(uvdist[select][uvtype[select] == 1], 
                        (uvreal[select][uvtype[select] == 1]**2 + uvimag[select][uvtype[select] == 1]**2)**0.5, 
                        s = .1, label = telescp[1])
    else:
        ax[1,0].scatter(uvdist[select][uvtype[select] == 0], 
                        (uvreal[select][uvtype[select] == 0]**2 + uvimag[select][uvtype[select] == 0]**2)**0.5, 
                        s = .1, label = telescp[0])

    ax[1,0].errorbar((bins[1:]+bins[:-1])/2, UVabsbinned, xerr = (bins[1:] - bins[0:-1])/2, yerr = UVabserrors, c = 'C2', ls = '', marker = 'o')
    ax[1,0].axis(ymax = np.amax(UVabsbinned) + 2*np.amax(UVabserrors), ymin = np.amin(UVabsbinned) -2* np.amax(UVabserrors))   

    if residue == False:
        ax[1,0].set_ylabel('|V|$_{model}$ [mJy]')   
    else:
        ax[1,0].set_ylabel('|V|$_{obs}$ - |V|$_{model}$ [mJy]')
        
    ax[1,0].set_xlabel('uv distance [k$\lambda$]') 
    ax[1,0].set_xscale('log')
    ax[1,0].legend(loc=3)
    
    # u-v coverage with amplitude as colorcoding --> Goes wrong with colormap
    if len(np.unique(uvtype))>1:
        pos = ax[1,1].scatter(us[select][uvtype[select] == 0]/1e3, 
                              vs[select][uvtype[select] == 0]/1e3, 
                              c = uvreal[select][uvtype[select] == 0], s = 0.1, marker = '.', label = telescp[0])
        ax[1,1].scatter(us[select][uvtype[select] == 1]/1e3, 
                        vs[select][uvtype[select] == 1]/1e3, 
                        c = uvreal[select][uvtype[select] == 1], s = 0.1, marker = 's', label = telescp[1])
    else:
        pos = ax[1,1].scatter(us[select][uvtype[select] == 0]/1e3, 
                              vs[select][uvtype[select] == 0]/1e3, 
                              c = uvreal[select][uvtype[select] == 0], s = 0.1, marker = '.', label = telescp[0])
    ax[1,1].axhline(0, c = 'gray', ls = 'dashed')
    ax[1,1].axvline(0, c = 'gray', ls = 'dashed')
    ax[1,1].set_ylabel('v [k$\lambda$]')
    ax[1,1].set_xlabel('u [k$\lambda$]') 
    ax[1,1].legend(loc = 1)
    cbar = fig.colorbar(pos, ax=ax[1,1])
    if residue == False:
        cbar.set_label(r'Re(V)$_{obs}$ [mJy]')
    else:
        cbar.set_label(r'Re(V)$_{obs}$ - Re(V)$_{model}$ [mJy]')

    plt.tight_layout()
    if residue == False:
        plt.savefig('./plots/outputs/SZuvplots_model_'+(filename_samples.split('/')[-1]).split('_pickle')[0]+ '.png', dpi = 300)
    else:
        plt.savefig('./plots/outputs/SZuvplots_residue_'+(filename_samples.split('/')[-1]).split('_pickle')[0]+ '.png', dpi = 300)
    plt.close()
    
    
    ###################################

    
    
    