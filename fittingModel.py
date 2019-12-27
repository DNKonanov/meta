import numpy as np
import scipy as sp
from spectraClass_mod import NmrSpectra
from scipy.optimize import minimize
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

class Library:

    def __init__(self, assignments, freq):

        self.library = []
        self.freq = freq

        for a in assignments:
            self.library.append(NmrSpectra(assignments[a], name=a, frequency=freq))



    def update_library(self, shifts):

        for i in range(len(shifts)):
            for j in range(len(shifts[i])):

                self.library[i].assignment['peak{}'.format(j+1)]['center'] += shifts[i][j]


    def generate_sum_spectra(self, ALPHA):

        spectra = np.zeros(int(12*self.freq/10**6))

        for i in range(len(self.library)):
            spectra += ALPHA[i] * np.array(self.library[i].generate_spectra())

        
        return spectra[::-1]
        


class FittingModel:

    def __init__(self, spectra, library, meta, frequency, ShiftsMatrix=None, ConcMatrix=None):

        self.spectra = spectra
        self.library = Library(library, freq=frequency)
        self.metainfo = meta

        

        if ShiftsMatrix is not None:
            self.PHI = ShiftsMatrix

        if ConcMatrix is not None:
            self.ALPHA = ConcMatrix


    def _identify_compounds(self):

        P = max(self.spectra)/100
        peaks, _ = find_peaks(self.spectra, height=P, prominence=P)

        print(peaks)
        return peaks
    

    def _fit_function(self, INPUT, alpha_len, phi_structure):

        ALPHA, linePHI = INPUT[:alpha_len], INPUT[alpha_len:]

        sl = 0
        PHI = []
        for i in phi_structure:
            #PHI.append([0 for i in range(i)])
            PHI.append(linePHI[sl:sl+i])
            sl += i
        print(ALPHA, PHI)

        self.library.update_library(PHI)
        
        result = np.sum((self.spectra - self.library.generate_sum_spectra(ALPHA))**2)
        result += 1000*np.sum(linePHI)
        print(result)
        return result

    def optimize(self, method='nelder-mead'):

        initial_ALPHA = list(np.zeros(len(self.library.library)) + 1.5)

        initial_PHI = [[0 for j in range(len(self.library.library[i]))] for i in range(len(self.library.library))]

        phi_structure = [len(i) for i in initial_PHI]

        INPUT = initial_ALPHA + [0 for i in range(sum(phi_structure))]

        res = minimize(
            self._fit_function, INPUT, 
            args=(len(initial_ALPHA), phi_structure), 
            method=method,
            bounds=[
                (0, np.inf) for i in range(len(initial_ALPHA))
            ] + [
                (-0.1, 0.1) for i in range(sum(phi_structure))
            ])

        plt.plot(self.library.generate_sum_spectra(res.x))
        plt.show()

        return res





#from libdict import library


#reflib = Library(library, 500_000_000)
#spectra = reflib.generate_sum_spectra([1 for i in range(len(reflib.library))])

#fitModel = FittingModel(spectra, library, '', 500_000_000)

#fitModel.optimize(method='trust-constr')

#SLSQP is FASTEST