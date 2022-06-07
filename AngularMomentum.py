# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 14:16:15 2016
@author: Ariel Sommer

Finds angular momentum matrices for a single spin
"""
import numpy as np
pi = np.pi

def _getN(J):
    #Find the number of states and check J
    if not(2*J == int(2*J)):
        raise ValueError("J must be half of an integer")
    return int(2*J+1)
    
def AngIndex(J,m):
    """
    Finds the index of m in the standard basis for spin J
    """
    ind_f = m+J #formula to go from m to index, but potentially gives floats
    ind_i = int(ind_f) #make it an int
    assert ind_f == ind_i #make sure we were given valid inputs
    return ind_i

def AngState(J,m):
    """
    Returns the state vector for m in the standard basis for spin J
    """
    N = _getN(J)
    state = np.zeros(N, dtype=np.complex)
    state[AngIndex(J,m)] = 1
    return state

#Angular momentum operators (divided by hbar)
def AngLadder(J,sign):
    """
    Finds the ladder operator for spin J
    sign = 1: raising
    sign =-1: lowering
    
    assumes the basis |J m_J> order starting with m_J = -J
    """    
    if sign != 1 and sign != -1:
        raise ValueError("sign must be +- 1")
    
    N = _getN(J)
    Mat = np.zeros((N,N),dtype=np.complex)
    mRange = np.arange(-J,J) if sign == 1 else np.arange(-J+1,J+1)
    for m in mRange:
        Mat[AngIndex(J,m+sign),AngIndex(J,m)] = np.sqrt(J*(J+1)-m*(m+sign))
    return Mat

def AngZ(J):
    """
    Operator matrix for z component of angular momentum in standard basis
    """
    N = _getN(J)
    Mat = np.zeros((N,N),dtype=np.complex)
    for m in np.arange(-J,J+1):
        Mat[AngIndex(J,m),AngIndex(J,m)] = m
    return Mat
    
def Id(J):
    return np.eye(_getN(J))

def AngX(J):
    return 0.5*(AngLadder(J,1) + AngLadder(J,-1))

def AngY(J):
    return -1j*0.5*(AngLadder(J,1) - AngLadder(J,-1))    

def Dot(I,J):
    """
    dot product operator for spins with magnitudes I and J
    """
    return np.kron(AngX(I), AngX(J)) + np.kron(AngY(I), AngY(J))+np.kron(AngZ(I), AngZ(J))

if __name__ == "__main__":
    J = .5
    state = AngState(J,-.5)
    print( AngLadder(J,1).dot(state))
    print( AngX(J))
    print( AngY(J))
    print( AngZ(J))
#    import matplotlib.pyplot as plt
    
