'''
Created on Oct 10, 2019

@author: Hiroki
'''
import matplotlib.pyplot as plt
import numpy as np
from operators.triangle_plaquette_hamiltonian import TrianglePlaquetteHamiltonian


plaqH = TrianglePlaquetteHamiltonian(1, 0.5, 2)

plt.matshow(np.abs(plaqH.matrix_form()))
plt.show()