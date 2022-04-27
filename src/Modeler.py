from src.Wrapper import *
from src.MsManager import *
from src.UnitTransform import *
from src.FourierTransform import *

from astropy.io import fits
from astropy.wcs import WCS

class Modeler:
    def __init__(self, 
                 filename_samples, 
                 filename_ms, 
                 obs_type, 
                 obs_amount,
                 outputdir = './output/',
                 load_pb = False,
                 save = True):
        
        self.file = filename_samples
        self.msfile = filename_ms
        self.outputdir = outputdir
        self.save = save
        
        reader = FileReader()
        self.popt = reader.execute(self.file)
        self.info = getinfo()
                
        self.obs_type = obs_type
        self.obs_amount = obs_amount

        self.msreader = MsReader(self.msfile, self.obs_type)        
        if load_pb == False: self.msreader.prep_pb()
    
    def _load_pb(self, spw, field):
        outvis = self.msreader.binvis.replace('-fid','-{0}'.format(field))
        outvis = outvis.replace('-sid','-{0}'.format(spw))
                
        pb, he_pb = fits.getdata('{0}.pbeam.fits'.format(outvis.replace('.ms','.im')), header = True)
        pb_im = pb[0,0] 
        pb_im[np.isnan(pb_im)] = 0.0
        
        self.pb_im = pb_im
        self.pb_hb = he_pb

    ################################################################################

    def _partouv(self):
        spwpoints = np.zeros(np.shape(self.uvdata[0]),dtype=np.complex128)
        
        for i in np.unique(self.popt['Link']):
            model_type, spectrum_type = np.array(list(self.popt)[:len(self.popt['Link'])])[(np.array(self.popt['Link']) == i)]
            
            if COMPONENTS[model_type.split('_')[0]]['make_image'] == True:
                continue
                        
            spwpoints = self.FT.partouv(spwpoints, self.popt[model_type], self.popt[spectrum_type], model_type, spectrum_type)
        return spwpoints
    
    def _get_analytical_model(self): 
        vis_model = np.empty(0)
        datas    = np.empty((3,0))
        
        print('Starting to build Parameter based model for {0}'.format(self.obs_type))
        for field in self.msreader.fields:
            print('- Processing field {0}'.format(field))
            for spw in self.msreader.spws:
                print('-- Spectral window {0}'.format(spw))

                self._load_pb(spw, field)
                self.uvdata = self.msreader.uvloader(spw, field) 
                self.FT     = FT(self.uvdata, self.pb_hb, self.pb_im, self.info)
                
                vis_model = np.append(vis_model, self._partouv())
                datas     = np.hstack((datas, self.uvdata))
                
        print()
        return vis_model, datas
    
    ################################################################################
        
    def _make_grid(self, model):
        wcs_pb = WCS(self.pb_hb, naxis = 2)
        
        x, y  = np.meshgrid(np.arange(0, len(self.pb_im), 1.), np.arange(0, len(self.pb_im), 1.))
        x += 0.5
        y -= 0.5
        
        sky = wcs_pb.pixel_to_world(x, y)
        RA  = sky.ra*u.deg
        Dec = sky.dec*u.deg
        
        r  = (Dec.value -  -1*model['Dec'])**2 
        r += ((RA.value -     model['RA'])*np.cos(np.deg2rad(-1*model['Dec'])))**2
        r  = r**0.5 
        
        return r
    
    def _make_modelimage(self, vis_model, index):

        for i in np.unique(self.popt['Link']):
            model_type, spectrum_type = np.array(list(self.popt)[:len(self.popt['Link'])])[(np.array(self.popt['Link']) == i)]                
            if COMPONENTS[model_type.split('_')[0]]['make_image'] == True:   
                
                grid  = self._make_grid(self.popt[model_type])
                
                transform = TransformInput(self.popt[model_type], model_type)
                input_par = transform.generate() 
                
                coord = np.deg2rad(grid) * info.cosmo.angular_diameter_distance(popt[t]['z'])
                coord = coord.value/input_par['major']
                                
                if model_type.split('_')[0] == 'gnfwPressure':
                    rs = self.info.linmesh[1:]
                else:
                    rs = self.info.linmesh
                    
                integrated_P = COMPONENTS[model_type.split('_')[0]]['function'](rs, **input_par)
                image        = np.interp(coord, rs, integrated_P) * self.info.ysznorm.value
                
                vis_model[index:index+len(self.uvdata[0])] = vis_model[index:index+len(self.uvdata[0])] + self.FT.imtouv(image, 
                                                                                                                         self.popt[model_type], 
                                                                                                                         self.popt[spectrum_type], 
                                                                                                                         model_type, 
                                                                                                                         spectrum_type)

        return vis_model
    
    def _get_imgplane_model(self, vis_model):
        index = 0
        print('Starting to build Image plane model for {0}'.format(self.obs_type))
        for field in self.msreader.fields:
            print('- Processing field {0}'.format(field))
            for spw in self.msreader.spws:
                print('-- Spectral window {0}'.format(spw))
                self._load_pb(spw, field)
                self.uvdata = self.msreader.uvloader(spw, field) 
                self.FT     = FT(self.uvdata, self.pb_hb, self.pb_im, self.info)
                vis_model = self._make_modelimage(vis_model, index)
                index += len(self.uvdata[0])

        return vis_model
    
    ################################################################################
    
    def _save_module(self, vis_model, savetype):
        if savetype == 'par': ## Quickfix + hardcoded for only two visibillities
            if self.obs_amount > 1:
                if self.obs_type == 'com07m':
                    scale = self.popt['Scaling']['scale0']
                elif self.obs_type == 'com12m':
                    scale =  self.popt['Scaling']['scale1']
            else:
                scale = self.popt['Scaling']['scale0']

            print("Saving model to ms-file")
            self.msreader.model_to_ms(vis_model, scale, todo = 'replace', savename = (self.file.split('/')[-1]).split('_pickle')[-2])
            print("Saving residue to ms-file")
            self.msreader.model_to_ms(vis_model, scale, todo = 'subtract', savename = (self.file.split('/')[-1]).split('_pickle')[-2])
            
        elif savetype == 'img':
            pass
        
    
    def run(self, onlySZ = True):        
        vis_model, datas = self._get_analytical_model()
        
        if onlySZ:
            vis_model        = self._get_imgplane_model(np.zeros_like(vis_model, dtype = np.complex128))
        else:
            vis_model        = self._get_imgplane_model(vis_model)

        if self.save: self._save_module(vis_model, savetype = 'par')
            
        return vis_model, datas
