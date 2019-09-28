'''
Created on Sep 27, 2019

@author: kwibu
'''
import numpy as np
from qiskit import QuantumRegister, QuantumCircuit, BasicAer, execute

backend = BasicAer.get_backend('unitary_simulator')

qr = QuantumRegister(1)
circ = QuantumCircuit(qr)

circ.u1(np.pi/3, qr[0])
circ.u3(np.pi, 0, np.pi, qr[0])

circ.u1(-np.pi/3, qr[0])
circ.u3(np.pi, 0, np.pi, qr[0])
circ.draw(filename='unitary.jpg', output='mpl')
job = execute(circ, backend)
print(job.result().get_unitary(circ, decimals=3))

