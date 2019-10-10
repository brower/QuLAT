'''
Created on Sep 26, 2019

@author: kwibu
'''
import numpy as np
import qiskit
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister

def trotter_pauli(q_circuit, qr, pauli_op, deltaT):
    """
    Trotter step for one pauli term
    Arguments:
    q_circuit: qiskit.QuantumCircuit, the target circuit
    qr: qiskit.QunatumRegister, the input qubit state
    pauli_op: operators.pauli_hamiltonian.PauliOperator, the target pauli operator
    deltaT: small time step
    Return:
    qiskit.QuantumCircuit, the circuit after added the trotterization step.
    """ 
    target_site_list = list(sorted(pauli_op.pauli.keys()))
    # Add a H or Rx(pi/2) gate to qubits corresponding to X or Y operators, respectively. 
    for j in target_site_list:
        if pauli_op.pauli[j] == "X":
            q_circuit.h(qr[j])
        elif pauli_op.pauli[j] == "Y":
            q_circuit.rx(np.pi/2, qr[j])
    
    # Add cnot gates to target qubits
    for i in range(len(target_site_list)-1):
        q_circuit.cx(qr[target_site_list[i]], qr[target_site_list[i+1]])
    # Add e^(-i*coef*dt*Z) to the last qubit
    q_circuit.u1(pauli_op.coef*deltaT, qr[target_site_list[-1]])
    q_circuit.u3(np.pi, 0, np.pi, qr[target_site_list[-1]])
    q_circuit.u1(-pauli_op.coef*deltaT, qr[target_site_list[-1]])
    q_circuit.u3(np.pi, 0, np.pi, qr[target_site_list[-1]])
    
    # Add cnot gates to target qubits
    for i in range(len(target_site_list)-1)[::-1]:
        q_circuit.cx(qr[target_site_list[i]], qr[target_site_list[i+1]])
    
    # Uncompute
    for j in target_site_list:
        if pauli_op.pauli[j] == "X":
            q_circuit.h(qr[j])
        elif pauli_op.pauli[j] == "Y":
            q_circuit.rx(-np.pi/2, qr[j])
    return q_circuit
        
        
    
def trotter_circuit(q_circuit, qr, pauli_h, T, n_steps):
    """
    Add a quantum circuit for Trotterization for pauli_h Hamiltonian with T time evolution to q_circuit. (e^(-iHT))
    Arguments:
    q_circuit: qiskit.QuantumCircuit, the target circuit
    qr: qiskit.QunatumRegister, the input qubit state
    pauli_h: operators.pauli_hamiltonian.PauliHamiltonian, the target hamiltonian
    T: float, evolution time
    n_steps: int, number of Trotterization steps
    
    Return:
    qiskit.QuantumCircuit, the circuit after added the trotterization operation. 
    """
    for d in range(n_steps):
        for p in pauli_h.pauli_op_list:
            q_circuit = trotter_pauli(q_circuit, qr, p, T/n_steps)
    return q_circuit
    
        
    