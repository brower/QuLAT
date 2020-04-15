'''
Created on Mar 27, 2020

@author: kwibu
'''
import numpy as np

def expectation(counts, shots):
    """
    Calculate the empirical expectation value from the measurement counts
    
    Arguments:
    counts: the frequency of the measured values of a quantum circuit obtained by execute(***).result().get_counts()
    shots: int, number of shots
    """
    
    keys = list(counts.keys())
    coefs = np.ones((len(keys)))
    for i, key in enumerate(keys):
        coefs[i] = 1 if key.count('1')%2 == 0 else -1
        
    return sum([counts.get(key, 0)/shots*coefs[i] for i, key in enumerate(keys)])