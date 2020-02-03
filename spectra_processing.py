import argparse
import nmrglue as ng
import matplotlib.pyplot as plt
import os
import subprocess
from scipy.signal import find_peaks
import numpy as np
from phase_correction import auto_phase_correction

import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()

parser.add_argument('-d', '--dir', type=str, required=True, help='path to folder')
parser.add_argument('-o', '--out', type=str, required=True, help='output csv file')
parser.add_argument('-r', '--results', type=str, required=True, help='final results file')
parser.add_argument('-ncores', type=int, default=1, help='num of cores for quantification')
parser.add_argument('-library', type=str, default='default', help='path to preset file')

args = parser.parse_args()


files = os.listdir(args.dir)

files.sort()

spectras = {file : None for file in  files}

lines = ['']


_N = len(files)
n = 1
for file in files:
    
    print('{} processing... ({} from {})'.format(file, n, _N))
    n += 1

    dic, data = ng.bruker.read('{}/{}'.format(args.dir, file))
    data = ng.bruker.remove_digital_filter(dic, data)

    udic = ng.bruker.guess_udic(dic, data)
    uc = ng.fileiobase.uc_from_udic(udic)
    ppm_scale = uc.ppm_scale()

    data = ng.proc_base.fft(data)
    
    data = auto_phase_correction(data)
    data = ng.process.proc_bl.baseline_corrector(data)
    
    lines[0] += file + ' '
    N = len(ppm_scale) - 1
    
    for i in range(len(data)):
        try:
            lines[i+1] += str(np.real(data[i])) + ' '
        except IndexError:
            lines.append('"{}" '.format(ppm_scale[N - i]))
            lines[i+1] += str(np.real(data[i])) + ' '
    print()
            

outfile = open(args.out, 'w')

print()
print('CSV file writing...')

for l in lines:
    try:
        ppm = float(l.split(' ')[0][1:-1])
    except:
        ppm = 1
        pass
    
    outfile.write(l[:-1] + '\n')
    

outfile.close()
print('Complete!')
print()
print('Call quantify_metabolites.r script...')

splitDir = args.out.split('/')
outFolder = '/'.join(splitDir[:-1])
outCSV = splitDir[-1]

subprocess.call('Rscript quantify_metabolites.r -dir {} -csv {} -out {} -ncores {} -library {}'.format(
    outFolder, 
    outCSV, 
    args.results, 
    args.ncores,
    args.library), 
    shell=True)