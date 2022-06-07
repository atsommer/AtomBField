# -*- coding: utf-8 -*-
"""
Created on Mon May 22 10:32:50 2017

@author: Ariel Sommer

Hyperfine structure of non-rotating molecule

"""
#import matplotlib.pyplot as plt    
import numpy as np
import AngularMomentum as Ang
try:
    import numericalunits as units
except ImportError:
    print "numericalunits module not found, using simple unit system" 
    class Units:
        pass
    units = Units()
    units.MHz = 1.
    units.G = 1.

"""
H = S dot (a1 I1 + a2 I2)
basis: (I1, I2), S
"""
muB = 1.3996245*units.MHz/units.G #Bohr magneton in frequency units

def dag(A):
    return np.conjugate(np.transpose(A))

def sorted_eigs(H):
    vals, vecs = np.linalg.eig(H)
    valsr = [v.real for v in vals]
    vecs = np.transpose(vecs)
    pairs = sorted(zip(valsr,vecs), key=lambda x: x[0])
    vals2 = [p[0] for p in pairs]
    N=len(vals)
    vecs2 = [p[1].reshape(N,1) for p in pairs]
    return vals2, vecs2

def kron3(A,B,C):
    return np.kron(np.kron(A,B),C)    

class Molecule(object):
    def __init__(self,S,I1,I2,gS,gI1,gI2,a1,a2):
        self.S=S
        self.I1=I1
        self.I2=I2
        self.gS=gS
        self.gI1=gI1
        self.gI2=gI2
        self.a1=a1
        self.a2=a2

        eyeI1 = np.eye(Ang._getN(I1))
        eyeI2 = np.eye(Ang._getN(I2))
        eyeS =np.eye( Ang._getN(S))
        eyeI = np.kron(eyeI1, eyeI2)
        
        I1x = kron3(Ang.AngX(I1), eyeI2, eyeS)
        I1y = kron3(Ang.AngY(I1), eyeI2, eyeS)
        I1z = kron3(Ang.AngZ(I1), eyeI2, eyeS)
        
        I2x = kron3(eyeI1, Ang.AngX(I2), eyeS)
        I2y = kron3(eyeI1, Ang.AngY(I2), eyeS)
        I2z = kron3(eyeI1, Ang.AngZ(I2), eyeS)
        
        self.Sx = Sx = np.kron(eyeI, Ang.AngX(S))
        self.Sy = Sy = np.kron(eyeI, Ang.AngY(S))
        self.Sz = Sz = np.kron(eyeI, Ang.AngZ(S))
        
        
        self.muB=muB
        # S-I
        self.SdotI1 = I1x.dot(Sx) + I1y.dot(Sy) + I1z.dot(Sz)
        self.SdotI2 = I2x.dot(Sx) + I2y.dot(Sy) + I2z.dot(Sz)
#        self.SdotI1 = np.kron(I1x, Ang.AngX(S)) +np.kron(I1y, Ang.AngY(S)) +np.kron(I1z, Ang.AngZ(S))
#        self.SdotI2 = np.kron(I2x, Ang.AngX(S)) +np.kron(I2y, Ang.AngY(S)) +np.kron(I2z, Ang.AngZ(S))
        
        self.F1x = F1x = I1x + Sx
        self.F1y = F1y = I1y + Sy
        self.F1z = F1z = I1z + Sz
        
        self.F2x = F2x = I2x + Sx
        self.F2y = F2y = I2y + Sy
        self.F2z = F2z = I2z + Sz
        
        self.Fx = Fx = I1x + I2x + Sx
        self.Fy = Fy = I1y + I2y + Sy
        self.Fz = Fz = I1z + I2z + Sz
        
        self.F1squared = F1x.dot(F1x) + F1y.dot(F1y) + F1z.dot(F1z)
        self.F2squared = F2x.dot(F2x) + F2y.dot(F2y) + F2z.dot(F2z)
        self.Fsquared = Fx.dot(Fx) + Fy.dot(Fy) + Fz.dot(Fz)
        
        self.N = Ang._getN(S)*Ang._getN(I1)*Ang._getN(I2)
        self.inited=0
        self.init()
    
    def init(self):
        self.H_HF = self.a1*self.SdotI1 + self.a2*self.SdotI2
        self.mu = self.muB*self.gS*self.Sz
        self.inited=1
        
    def eig(self, B):
        # Find eigenvalues and eigenvectors at magnetic field B
        self.init()
        H = self.H_HF + B*self.mu
        self.vals, self.vecs = sorted_eigs(H)
        return self.vals, self.vecs
    
    def S_matels(self, vecs, ind1, ind2):
        """
        Finds the matrix elements of the electron spin operator
        """
        Xmatel = dag(vecs[ind1]).dot(self.Sx).dot(vecs[ind2])[0,0]    
        Ymatel = dag(vecs[ind1]).dot(self.Sy).dot(vecs[ind2])[0,0]
        Zmatel = dag(vecs[ind1]).dot(self.Sz).dot(vecs[ind2])[0,0]
        return Xmatel, Ymatel, Zmatel
    
    def Fsquared_matel(self, vecs, ind1, ind2):
        """
        Finds the matrix element of the F^2 operator
        """
        return dag(vecs[ind1]).dot(self.Fsquared).dot(vecs[ind2])[0,0]
        
    def F1squared_matel(self, vecs, ind1, ind2):
        """
        Finds the matrix element of the F1^2 operator
        """
        return dag(vecs[ind1]).dot(self.F1squared).dot(vecs[ind2])[0,0]
        
    def F2squared_matel(self, vecs, ind1, ind2):
        """
        Finds the matrix element of the F2^2 operator
        """
        return dag(vecs[ind1]).dot(self.F2squared).dot(vecs[ind2])[0,0]