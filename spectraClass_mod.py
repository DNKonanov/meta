import numpy as np
import matplotlib.pyplot as plt

#test data
test_data = {
    'name': 'test_name',
    'peak1': {
        'sum_intensity': 3,
        'center': 2,
        'sigma': 2,
        'J': [20,20, 15]
    },


    'peak2': {
        'sum_intensity': 1,
        'center': 4,
        'sigma': 3,
        'J': [15,20]
    }
}


class NmrSpectra:

    def __init__(self, assignment, name='DefaultName', frequency=None):

        self.name = name
        self.frequency = frequency
        self.assignment = assignment

    
    @staticmethod
    def compute_peaks(center, J):

        pl, pr = center - J, center + J

        Il, Ir = 1, 1

        return pl, pr, Il, Ir


    def generate_spectra(self):

        spectra = np.zeros(int(12*self.frequency/10**6))

        for peak_key in self.assignment:
            if peak_key == 'name':
                continue
            
            peak = self.assignment[peak_key]
            sum_peak = []
            ppm = self.frequency - peak['center']*self.frequency/10**6
            
            sum_peak.append(
                (ppm, peak['sum_intensity'])
            )

            for J in peak['J']:

                new_sum_peak = []

                for p in sum_peak:

                    pl, pr, Il, Ir = self.compute_peaks(p[0], J)
                    

                    new_sum_peak.extend( 
                        (
                            (pl, Il),
                            (pr, Ir),
                        ) 
                    )
                
                
                sum_peak = new_sum_peak

            for i in sum_peak:

                intensity = i[1]*peak['sum_intensity']/len(sum_peak)
                
                coord = int(self.frequency - i[0])
                
                for j in range(-10*peak['sigma'], 10*peak['sigma']):

                    try:
                        spectra[j+coord] += self.lorenz(coord, peak['sigma'], j+coord)*intensity
                    except IndexError:
                        pass

        return spectra

    def __repr__(self):
        return str(self.assignment) + '\n'

    def __len__(self):
        return len(self.assignment)
        
    @staticmethod
    def normal(mu, sigma, coord):

        return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (coord - mu)**2 / (2 * sigma**2))

    @staticmethod
    def lorenz(mu, sigma, coord):

        return 1/(np.pi*(1 + (coord-mu)**2/sigma))

    def plot(self):

        spectra = self.generate_spectra()

        plt.style.use('seaborn')
        plt.plot(

            [i*10**6/self.frequency for i in range(len(spectra))],

            spectra)

        plt.xlabel('Chemical shift, ppm')
        plt.ylabel('Intensity')
        plt.show()




    
    