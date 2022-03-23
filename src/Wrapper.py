from src.Settings import *
import scipy
 
class FileReader:

    def _mapper(self, params):

        mapped_types = {}
        count = 0
        for idx, t in enumerate(self.types):
            mapped_params = {}
            for key, value in zip(self.labels[idx], params[count: count + len(self.labels[idx]) + 1]):
                mapped_params[key] = value
 
            count += len(self.labels[idx])
            if t == 'Scaling':
                mapped_types[t] = mapped_params
            else:
                mapped_types[t+'_'+str(self.link[idx])] = mapped_params

        return mapped_types
        
    def _read_combinedata(self, results, popt, fixd):
    
        params = []
        count = 0
        been_there = np.zeros(len(np.unique(self.types)))
    
        for idx, lbl in enumerate(self.labels):
            for jdx, l in enumerate(lbl):
                if fixd[idx][jdx]: 
                    params.append(popt[count])
                    count +=1
                else:
                    if COMPONENTS[self.types[idx]]['spectrum'] == False:
                        comp = 'model'
                    elif  (COMPONENTS[self.types[idx]]['spectrum'] == True) & (self.types[idx] != 'Scaling'):
                        comp = 'spectrum'
                    elif  (COMPONENTS[self.types[idx]]['spectrum'] == True) & (self.types[idx] == 'Scaling'):
                        comp = 'scale'
                    
                    if comp != 'scale':
                        index = int(idx/2)
                        params.append(results['pars'][index][comp]['guess'][jdx])    
        
                    elif comp == 'scale':
                        params.append(results['scales'][jdx]) 
        
        return params
        
    def _read_fixedvalues(self, results, popt):
        
        fixed = []
        for var in results['vary'][:-1]:
            for v in var['values']:
                fixed.append(list(var['values'][v]['vary']))
        
        scaling = []
        for scle in results['vary'][-1]['values']['vary']:
            scaling.append(scle)
        
        fixed.append(list(scaling))
        return fixed
    
    def _read_popt(self, results):
        samples = np.copy(results['samples']['samples'])
        
        weights = results['samples']['logwt']-scipy.special.logsumexp(results['samples']['logwt']-results['samples']['logz'][-1])
        weights = np.exp(weights-results['samples']['logz'][-1])
        popt    = np.average(samples, weights = weights, axis = 0)
        return popt
    
    
    def _read_labels(self,results):
        self.labels = []
        for t in self.types[:-1]:
            self.labels.append(COMPONENTS[t]["variables"])
    
        scaling = []
        for i in range(len(results['vary'][-1]['values']['vary'])):
            scaling.append('scale'+str(i))
            
        self.labels.append(scaling)
        
    def _read_types(self, results):
        self.link  = [] 
        self.types = []
        for idx, res in enumerate(results['pars']):
            for r in res:
                self.types.append(res[r]['type'])
                self.link.append(idx)

        self.types.append('Scaling')
        
    def _read(self, filename):
        pickled_file = np.load(filename, allow_pickle=True)   
        self._read_types(pickled_file) 
        self._read_labels(pickled_file) 

        if pickled_file['type'] == 'dynesty':
            popt = self._read_popt(pickled_file)
            fixd = self._read_fixedvalues(pickled_file, popt)
            params = self._read_combinedata(pickled_file, popt, fixd)
            
            return params
        else:
            printError('Provide correct sampler type')

    def execute(self, filename):
        params = self._read(filename)
        mapped_data = self._mapper(params)
        mapped_data['Link'] = self.link
        return mapped_data