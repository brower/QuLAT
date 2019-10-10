'''
Created on Sep 27, 2019

@author: kwibu
'''
import numpy as np

Sigma = {"I": np.eye(2),
         "X": np.array([[0, 1], [1, 0]]),
         "Y": np.array([[0, -1j], [1j, 0]]),
         "Z": np.array([[1, 0], [0, -1]])}


