'''
Created on Jan 14, 2020

@author: kwibu
'''
import tensornetwork as tn
import numpy as np

"""
Compute the ground state of 1-d 4-site Ising model with periodic boundary condition using MPS. 

Hamiltonian:
H = -Jz sum_i Z_i Z_{i+1} -J sum_i X_i
"""

# create nodes
nodes = []
for i in range(4):
    nodes.append(tn.Node(np.zeros((2, 2, 2))))
    

# construct the graph
edges = []
for i in range(4):
    edges.append(nodes[i][2]^nodes[(i+1)%4][0])
    