'''
Created on Sep 29, 2019

@author: Hiroki
'''
from operators.pauli_hamiltonian import PauliHamiltonian

class TrianglePlaquetteHamiltonian(PauliHamiltonian):
    def __init__(self, g, alpha, n_layers):
        """
        Triangle plaquette Hamiltonian for simulation of single triangle quantum link model. 
        
        Arguments:
        g: float, coefficient for the electric and plaquette terms
        alpha: float, coefficient for the coupling term
        n_layers: int, number of layers
        """
        n_sites = n_layers*3
        coef_list = []
        pauli_list = []
        for s in range(n_layers):
            # Electric terms
            for j in range(3):
                coef_list.append(g**2/2)
                pauli_list.append({j+(s%n_layers)*3: "Z", j+((s+1)%n_layers)*3: "Z"})
                
            # Coupling terms
            for j in range(3):
                coef_list += [alpha, alpha]
                pauli_list.append({j+(s%n_layers)*3: "X", j+((s+1)%n_layers)*3: "X"})
                pauli_list.append({j+(s%n_layers)*3: "Y", j+((s+1)%n_layers)*3: "Y"})
                
            # Plaqutte terms
            coef_list += [-1/(2*g**2), 1/(2*g**2), 1/(2*g**2), 1/(2*g**2)]
            pauli_list += [{0+(s%n_layers)*3: "X", 1+(s%n_layers)*3: "X", 2+(s%n_layers)*3: "X"}, 
                           {0+(s%n_layers)*3: "X", 1+(s%n_layers)*3: "Y", 2+(s%n_layers)*3: "Y"}, 
                           {0+(s%n_layers)*3: "Y", 1+(s%n_layers)*3: "X", 2+(s%n_layers)*3: "Y"},
                           {0+(s%n_layers)*3: "Y", 1+(s%n_layers)*3: "Y", 2+(s%n_layers)*3: "X"}]
        super(TrianglePlaquetteHamiltonian, self).__init__(coef_list, pauli_list, n_sites=n_sites)
        
