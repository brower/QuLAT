'''
Created on Sep 27, 2019

@author: kwibu
'''
import itertools
from operators.sigma_operators import Sigma
from utils.matrix_operations import tensor_prod
from scipy.linalg import eig
class PauliOperator:
    def __init__(self, coef, pauli, n_sites=None):
        """
        PauliOperator class. Defines a multi-site pauli operator.
        
        members:
        coef: float, coefficient corresponding to the pauli operator.
        paulis: dictionary of {int: string}, the key integer represents the site the sigma is applied, 
                and string ("X", "Y", or "Z") represents the type of the sigma. 
        n_sites (optional): int, number of lattice sites in the system. If None, the biggest site number in the pauli_list.
        matrix: matrix form of the operator. Notice that it is not computed until _compute_matrix function is run. 
        """
        self.coef = coef
        self.pauli = pauli
        self.n_sites = n_sites if n_sites is not None else max(pauli.keys())+1
        self.matrix = None
        
    def _compute_matrix(self):
        """
        Private
        Compute the matrix form of the operator and set it to the class member self.matrix.
        """
        mat_list = []
        # Operators on lower sites go to right hand side
        # E.g. X_2 otimes X_1, not the other way around.
        for j in range(self.n_sites)[::-1]:
            if j in self.pauli.keys():
                mat_list.append(Sigma[self.pauli[j]])
            else:
                mat_list.append(Sigma["I"])
        self.matrix = self.coef*tensor_prod(mat_list, sparse=True)
        self.matrix.eliminate_zeros()
        
    def matrix_form(self, sparse=False):
        """
        Return the matrix form of the operator.
        
        Return:
        numpy array, the matrix form of the pauli
        """
        if self.matrix is None:
            self._compute_matrix()
        return self.matrix if sparse else self.matrix.toarray()
    
    def eigensystem(self):
        """
        Return the eigenvalues and eigenvectors of the operator.
        
        Return:
        eigenvalues, eigenvectors: numpy array of eigenvalues and eigenvectors
        """
        if self.matrix is None:
            self._compute_matrix()
        return eig(self.matrix.toarray())
        
        
class PauliHamiltonian(object):
    def __init__(self, coef_list, pauli_list, n_sites = None):
        """
        PauliHamiltonian class. Defines a Hamiltonian in pauli basis (sum_k(h_k*P_k)). 
        
        Arguments:
        coef_list: list of float, list of coefficients for each pauli terms. 
        pauli_list: list of dictionary of {int: string}, the key integer of the dictionary represents 
                the site the sigma is applied, and string ("X", "Y", or "Z") represents the type of the sigma.
        n_sites (optional): int, number of lattice sites in the system. If None, the biggest site number in the pauli_list.
        matrix: matrix form of the hamiltonian. Notice that it is not computed until _compute_matrix function is run.  
        """
        assert(len(coef_list) == len(pauli_list)), "The length of the list of coefficients (%d) and the list of pauli operators (%d) do not match" %(len(coef_list), len(pauli_list))
        self.n_sites = n_sites if n_sites is not None else max(itertools.chain.from_iterable([list(p.keys()) for p in pauli_list]))+1
        self.pauli_op_list = []
        for coef, pauli in zip(coef_list, pauli_list):
            self.pauli_op_list.append(PauliOperator(coef, pauli, self.n_sites))
        self.matrix = None
        
    def _compute_matrix(self):
        """
        Private
        Compute the matrix form of the Hamiltonian and set it to the class member self.matrix.
        """
        self.matrix = sum([p.matrix_form(sparse=True) for p in self.pauli_op_list])
        
    def matrix_form(self, sparse=False):
        """
        Return the matrix form of the operator.
        
        Return:
        numpy array, the matrix form of the Hamiltonian
        """
        if self.matrix is None:
            self._compute_matrix()
        return self.matrix if sparse else self.matrix.toarray()
    
    def eigensystem(self):
        """
        Return the eigenvalues and eigenvectors of the Hamiltonian.
        
        Return:
        eigenvalues, eigenvectors: numpy array of eigenvalues and eigenvectors
        """
        if self.matrix is None:
            self._compute_matrix()
        return eig(self.matrix.toarray())
    
