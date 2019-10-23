'''
Created on Sep 27, 2019

@author: kwibu
'''
import numpy as np
from scipy import sparse

Sigma = {"I": sparse.csr_matrix(np.eye(2)),
         "X": sparse.csr_matrix(np.array([[0, 1], [1, 0]])),
         "Y": sparse.csr_matrix(np.array([[0, -1j], [1j, 0]])),
         "Z": sparse.csr_matrix(np.array([[1, 0], [0, -1]]))}


