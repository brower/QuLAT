'''
Created on Mar 27, 2020

@author: kwibu
'''
import numpy as np
def change_basis_measurement(circ, i, trg):
    """
    Change the basis of measurement by adding H gate for the X basis or RX(Pi/2) gate for the Y basis
    Arguments:
    circ: qiskit.QuantumCircuit, the target circuit
    i: int, the basis, 0: identity, 1: X, 2: Y, 3: Z
    trg: the index of the target qubit
    """
    assert(i < 4), "Invalid basis index is given: %d"%i
    if i == 1:
        circ.u2(0, np.pi, trg)
    elif i == 2:
        circ.u3(np.pi/2, -np.pi/2, np.pi/2, trg) # RX(Pi/2)
        
def change_basis_input(circ, i, trg):
    """
    Change the basis of the input state: RY(+-Pi/2):|0>->|+->, RX(+-Pi/2):|0>->|R L>, X:|0>->|1>
    Arguments:
    circ: qiskit.QuantumCircuit, the target circuit
    i: int, the basis, 0: identity, 1: X, 2: Y, 3: Z
    trg: the index of the target qubit
    """
    assert(i < 6), "Invalid basis index is given: %d"%i
    if i == 0:
        circ.u2(0, np.pi, trg) # H
    elif i == 1:
        circ.u3(-np.pi/2, 0, 0, trg) # RY(-Pi/2)
    elif i == 2:
        circ.u3(np.pi/2, -np.pi/2, np.pi/2, trg) # RX(Pi/2)
    elif i == 3:
        circ.u3(-np.pi/2, -np.pi/2, np.pi/2, trg) # RX(-Pi/2)
    elif i == 5:
        circ.u3(np.pi, 0, np.pi, trg) # X