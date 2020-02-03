from scipy.optimize import minimize
import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt

def auto_phase_correction(data):
    print('Phase correction...')

    f = minimize(_diff_score, x0=0, args=data)
    
    data = ng.proc_base.ps_exp(data, p0=f.x[0])

    print('Optimal phase correction angle = {}'.format(f.x[0]))
        
    return data
    
def _diff_score(p0, data):    
    data = ng.proc_base.ps_exp(data, p0=p0)
    
    fragment = np.real(data[6000:7000]) #extract DSS solvent
    fragment = fragment/np.max(fragment) #normilize data
    
    fragment = np.log(fragment - np.min(fragment)+10)
    index = np.where(fragment == np.amax(fragment))[0][0]
    
    left_side = fragment[index-50:index+1]
    right_side = fragment[index:index + 51][::-1]
   
    score = np.sum((left_side[:20] - right_side[:20])**2) # compute symmetry score
    return score
    