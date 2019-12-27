import argparse
from fittingModel import FittingModel
from raw_spectra_processing import open_spectra, delete_water, resize_data
from libdict import library
import sys

parser = argparse.ArgumentParser()

parser.add_argument(
    '-ifolder', 
    type=str, 
    help='path to folder with Bruker fid file', 
    required=True)
parser.add_argument(
    '-library', 
    type=str, 
    default='library/library', 
    help='path library file. Format is described in README. Default is AQuA library')
parser.add_argument(
    '-frequency',
    type=int,
    default=500_000_000    
)


args = parser.parse_args()

ppm_scale, spectra = open_spectra(args.ifolder)
spectra = delete_water(spectra)

ppm_scale, spectra = resize_data(ppm_scale, spectra, frequency=args.frequency)
print(spectra)



model = FittingModel(spectra, library, [], args.frequency)

c = model._identify_compounds()
 
cp = [ppm_scale[i] for i in c]

print(cp)

import matplotlib.pyplot as plt
plt.plot(ppm_scale, spectra)
plt.scatter(cp, [-0.001 for i in range(len(cp))])



plt.show()



sys.exit()

res = model.optimize()


print(res)

