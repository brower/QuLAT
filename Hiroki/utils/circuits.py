import numpy as np

def add_crx(q_circuit, theta, ctr, tgt):
    q_circuit.u1(np.pi/2, tgt)
    q_circuit.cx(ctr, tgt)
    q_circuit.u3(-theta, 0, 0, tgt)
    q_circuit.cx(ctr, tgt)
    q_circuit.u3(theta, 0, 0, tgt)
    q_circuit.u1(-np.pi/2, tgt)
    return q_circuit
    
def add_ccrx(q_circuit, theta, ctr1, ctr2, tgt):
    q_circuit = add_crx(q_circuit, theta/2, ctr2, tgt)
    q_circuit.cx(ctr1, ctr2)
    q_circuit = add_crx(q_circuit, -theta/2, ctr2, tgt)
    q_circuit.cx(ctr1, ctr2)
    q_circuit = add_crx(q_circuit, theta/2, ctr1, tgt)
    return q_circuit