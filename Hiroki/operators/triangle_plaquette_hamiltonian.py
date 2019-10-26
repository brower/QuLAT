'''
Created on Sep 29, 2019

@author: Hiroki
'''
import numpy as np
import scipy.sparse
import itertools
from operators.pauli_hamiltonian import PauliHamiltonian

class TrianglePlaquetteHamiltonian(PauliHamiltonian):
    def __init__(self, g, alpha, n_layers):
        """
        Triangle plaquette Hamiltonian for simulation of single triangle quantum link model. 
        
        Arguments:
        g: float, coefficient for the electric and plaquette terms
        alpha: float, coefficient for the coupling term
        n_layers: int, number of layers
        
        members:
        permuted_matrix: permuted Hamiltonian. Notice that it is not computed until gauge_rotation_basis function is run. 
        """
        self.n_layers = n_layers
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
            coef_list += [-1/(8*g**2), 1/(8*g**2), 1/(8*g**2), 1/(8*g**2)]
            pauli_list += [{0+s*3: "X", 1+s*3: "X", 2+s*3: "X"}, 
                           {0+s*3: "X", 1+s*3: "Y", 2+s*3: "Y"}, 
                           {0+s*3: "Y", 1+s*3: "X", 2+s*3: "Y"},
                           {0+s*3: "Y", 1+s*3: "Y", 2+s*3: "X"}]
        super(TrianglePlaquetteHamiltonian, self).__init__(coef_list, pauli_list, n_sites=n_sites)
        self.permuted_matrix = None
        
    def gauge_rotation_basis(self, sparse=False):
        """
        Compute Hamiltonian permutated by the values of eigenvalues of the gauge transformation operators: J12^2+J23^2+J31^2.
        
        Arguments:
        sparse: bool, True if you want a sparse matrix as return. 
        
        Return:
        scipy.sparse.csr_matrix(if sparse=False), or numpy array(if sparse=True), 
            the permutated Hamiltonian sparse = False
        """
        if self.permuted_matrix is None: 
            # Define gauge rotation operators
            J12 = PauliHamiltonian([1, -1]*self.n_layers, 
                                   list(itertools.chain.from_iterable([[{0+s*3: "Z"}, {1+s*3: "Z"}] for s in range(self.n_layers)])), n_sites=self.n_sites)
            J23 = PauliHamiltonian([1, -1]*self.n_layers, 
                                   list(itertools.chain.from_iterable([[{1+s*3: "Z"}, {2+s*3: "Z"}] for s in range(self.n_layers)])), n_sites=self.n_sites)
            J31 = PauliHamiltonian([1, -1]*self.n_layers, 
                                   list(itertools.chain.from_iterable([[{2+s*3: "Z"}, {0+s*3: "Z"}] for s in range(self.n_layers)])), n_sites=self.n_sites)
            
            # Get diagonal elements (eigenvalues) of rotation operators
            D12 = J12.matrix_form(True).diagonal()
            D23 = J23.matrix_form(True).diagonal()
            D31 = J31.matrix_form(True).diagonal()
            
            dtype=[('val', 'float'), ('ind', 'int')]
            squared_sum = np.array([(val, ind) for val, ind in enumerate(np.square(D12) + np.square(D23) + np.square(D31))], dtype=dtype)
            perm = np.argsort(squared_sum, order=['ind', 'val'])
            self.gauge_eigs = (np.square(D12) + np.square(D23) + np.square(D31))[perm]
            
            if self.matrix is None:
                self._compute_matrix()
            # Permute the Hamiltonian matrix
            coo_matrix = self.matrix.tocoo()
            coo_matrix.row = np.array([np.where(perm == ind)[0][0] for ind in coo_matrix.row])
            coo_matrix.col = np.array([np.where(perm == ind)[0][0] for ind in coo_matrix.col])
            self.permuted_matrix = coo_matrix.tocsr()
        return self.permuted_matrix if sparse else self.permuted_matrix.toarray()
    
    def block_sectors(self, sparse = False):
        """
        Return block sectors of Hamiltonian permutated taking account of eigenvalues of the gauge transformation operators: J12^2+J23^2+J31^2.
        
        Arguments:
        sparse: bool, True if you want a sparse matrix as return. 
        
        Return:
        scipy.sparse.csr_matrix(if sparse=False), or numpy array(if sparse=True), 
            the permutated Hamiltonian sparse = False
        """
        
        if self.permuted_matrix is None:
            self.gauge_rotation_basis(sparse=True)
        
        sectors = []
        for p in np.unique(self.gauge_eigs):
            ind = np.where(self.gauge_eigs == p)[0]
            sec = self.permuted_matrix[ind[:, None], ind]
            sectors.append(sec if sparse else sec.toarray())
            
        return sectors
    