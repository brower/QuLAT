'''
Created on Oct 10, 2019

@author: Hiroki
'''
import sys
import matplotlib.pyplot as plt
import numpy as np
from operators.triangle_plaquette_hamiltonian import TrianglePlaquetteHamiltonian
np.set_printoptions(threshold=sys.maxsize)

plaqH = TrianglePlaquetteHamiltonian(1, 1, 2)
print(plaqH.matrix_form())

plt.matshow(np.abs(plaqH.matrix_form()))
plt.show()