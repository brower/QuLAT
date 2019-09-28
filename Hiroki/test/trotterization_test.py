'''
Created on Sep 27, 2019

@author: kwibu
'''
import numpy as np
from scipy.linalg import expm, eig
from qiskit import QuantumCircuit, QuantumRegister, BasicAer, execute
from quantum_circuit.trotterization import trotter_circuit
from operators.pauli_hamiltonian import PauliOperator, PauliHamiltonian

# 2-site Heisenberg
H = PauliHamiltonian([1., -0.5, -0.5], [{0: "Z", 1: "Z"}, {0: "X"}, {1: "X"}], n_sites=2)
#H = PauliHamiltonian([1.], [{0: "Z", 1: "X"}], n_sites=2)
backend = BasicAer.get_backend('unitary_simulator')
qr = QuantumRegister(2, 'qr')
circ = QuantumCircuit(qr)
circ = trotter_circuit(circ, qr, H, 0.1, 10)
circ.draw(filename='unitary.jpg', output='mpl')
job = execute(circ, backend)
print(expm(-1j*H.matrix_form()*0.1))
print(job.result().get_unitary(circ, decimals=3))