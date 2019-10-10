import numpy as np

def commutator(A, B, op=np.matmul):
    """
    Compute [A, B] = A*B - B*A where * is defined as op
    
    Arguments:
    A: any python object
    B: any python object
    op: function (default=np.matmul), the binary operator *. Must be defined to take two arguments with the same class as A & B.
    
    Return:
    [A, B] 
    """
    return op(A, B) - op(B, A)


def anticommutator(A, B, op=np.matmul):
    """
    Compute {A, B} = A*B + B*A where * is defined as op
    
    Arguments:
    A: any python object
    B: any python object
    op: function (default=np.matmul), the binary operator *. Must be defined to take two arguments with the same class as A & B.
    
    Return:
    {A, B} 
    """
    return op(A, B) + op(B, A)

