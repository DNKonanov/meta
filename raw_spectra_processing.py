import nmrglue as ng
import numpy as np
def open_spectra(bruker_folder):

    print('Opening spectra from {}'.format(bruker_folder))
    dic, data = ng.bruker.read(bruker_folder)

    data = ng.bruker.remove_digital_filter(dic, data)

    #EXTRACT SCALING
    udic = ng.bruker.guess_udic(dic, data)
    uc = ng.fileiobase.uc_from_udic(udic)
    ppm_scale = uc.ppm_scale()

    data = ng.proc_base.fft(data)
    print('Phase correction...')
    data = phase_correction(data)
    #data = ng.proc_base.di(data)
    data = ng.proc_base.rev(data)
    print('Baseline correction...')
    data = ng.process.proc_bl.baseline_corrector(data)

    return ppm_scale, data.real/np.max(data.real)

def resize_data(ppm_scale, data, frequency=500_000_000):

    for i in range(1, len(ppm_scale)):

        if ppm_scale[i-1] > 12 >= ppm_scale[i]:
            left = i
        elif ppm_scale[i-1] > 0 >= ppm_scale[i]:
            right = i

    ppm_scale = ppm_scale[left:right+1]
    data = data[left:right + 1]
    N = int(frequency*12/1000000)
    M = len(data)
    nscale = []
    ndata = []
    for i in range(N):
        nscale.append(ppm_scale[int(i*M/N)])
        ndata.append(data[int(i*M/N)])

    return np.array(nscale), np.array(ndata) 



def phase_correction(data):
    # autocorrection should be implemented!!!
    data = ng.proc_base.ps_exp(data, p0=-55)
    return data

def delete_water(data):

    ndata = data
    ndata[15000:20000] = 0
    return data

def write():
    return