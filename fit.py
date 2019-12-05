import numpy as np
import scipy as sp
from spectraClass import NmrSpectra
from scipy.optimize import minimize


class Library:

    def __init__(self, assignmemts, freq):

        self.library = []
        self.freq = freq

        for a in assignmemts:
            self.library.append(NmrSpectra(a, name=a['name'], frequency=freq))



    def update_library(self, shifts):

        for i in range(shifts):
            for j in range(len(shifts[i])):
                self.library[i]['peak{}'.format(j)] += shifts[i][j]


    def generate_sum_spectra(self, ALPHA):

        spectra = np.zeros(int(12*self.freq/10**6))

        
        for i in range(len(self.library)):
            spectra += ALPHA[i] * np.array(spectra)

        
        return spectra
        


class FittingModel:

    def __init__(self, spectra, library, meta, frequency, ShiftsMatrix=None, ConcMatrix=None):

        self.spectra = spectra
        self.library = Library(library, freq=frequency)
        self.metainfo = meta

        

        if ShiftsMatrix is not None:
            self.PHI = ShiftsMatrix

        if ConcMatrix is not None:
            self.ALPHA = ConcMatrix
    

    def _fit_function(self, INPUT):

        PHI, ALPHA = INPUT[0], INPUT[1]

        self.library.update_library(PHI)

        return self.library.generate_sum_spectra(ALPHA)

    def optimize(self):

        ALPHA = np.zeros(len(self.library)) + 1
        PHI = [[0 for j in range(len(self.library.library[i]))] for i in range(len(self.library))]

        res = minimize(self._fit_function, (PHI, ALPHA), method='nelder-mead')

        return res




