import pandas as pd
import numpy as np
from spectraClass_mod import NmrSpectra
import matplotlib.pyplot as plt

f = open('library/library')


types = {
    's': (1,),
    'd': (2,),
    't': (3,),
    'dd': (2,2),
    'dt': (2,3),
    'q': (4,),
    'dq': (2,4,),
    'sept': (7,),
    'qvd': (5,2,),
    'o': (8,),
    'm': (),
    'six': (6,),
    'non': (9,),
    'septd': (7,2,),
}



skip = 0
library = {}
for line in f:
    string = line[:-1].split(' ')

    ref_id, ref_name, ref_peaks = string[0][:-1], string[1], ''.join(string[2:])

    peaks = ref_peaks.split(',')

    library[ref_name] = {
    }

    peak_num = 1

    for peak in peaks:
        ppm, properties = peak.split('(')[0], ''.join(peak.split('(')[1][:-1])

        _string = properties.split(';')

        peaktype, Js, intensity = _string[0], _string[1:-1], _string[-1][:-1]

        if peaktype == 'm':
            continue

        elif peaktype == 's':
            J_list = []

        else:
            J_list = []
            print (peak)
            peaktype = types[peaktype]

            for i in range(len(peaktype)):
                J = float(Js[i].split('=')[1])
                for _ in range(peaktype[i] - 1):
                        J_list.append(J)



        library[ref_name]['peak{}'.format(peak_num)] = {
            'center': float(ppm),
            'sum_intensity': int(intensity),
            'sigma': 2, #CHANGGGGGEEEEE!!!
            'J': J_list
        } 
        peak_num += 1


    
         

    print(ref_name)
    print(library[ref_name])

print(library)