'''
Created on Oct 10, 2019

@author: Hiroki
'''
import sys
sys.path.append("..")
import matplotlib.pyplot as plt
import numpy as np
from operators.pauli_hamiltonian import PauliHamiltonian
np.set_printoptions(threshold=sys.maxsize)

Jz = 1.0
Jx = 1.0
H = PauliHamiltonian([-Jz,-Jx,-Jx],[{0 : "Z",1 : "Z"},{0 : "X"},{1 : "X"}],2)
print(H.matrix_form())

plt.matshow(np.abs(H.matrix_form()))
plt.show()