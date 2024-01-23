# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 14:16:15 2016
@author: Ariel Sommer

Compute the energy eigenvalues of an alkali atom in a magnetic field
compare values to find imaging frequencies
eventually also find matrix elements

Model: fixed n,L,J
- essentially an electron spin J and nuclear spin I in a B field

"""
import matplotlib.pyplot as plt
import numpy as np
import AngularMomentum as Ang

try:
    import numericalunits as units
except ImportError:
    print( "numericalunits module not found, using simple unit system" )
    class Units:
        pass
    units = Units()
    units.MHz = 1.
    units.G = 1.

def landeG(L,L1,L2,g1,g2):
    def term(g,a,b):
        return g*(L*(L+1.)-b*(b+1.)+a*(a+1.))/(2.*L*(L+1.))
    return term(g1,L1,L2) + term(g2,L2,L1)

muB = 1.399624604*units.MHz/units.G #Bohr magneton / h
class AtomIJ(object):
    """
    Models an atom as a nuclear spin I and electron spin J
    Assumes we have chosen n,L
    """
    def __init__(self,name,L_label,I,J,gI, gJ, A, B=0, muBref=muB):
        self.name = name
        self.L_label = L_label
        self.I = I
        self.J = J
        self.gJ = gJ #electron Lande g factor
        self.gI = gI #nuclear g factor
        self.A = A #hyperfine constant A
        self.B = B #hyperfine constant B (electric quadrupole)
        self.muB = muBref #Bohr magneton
        self._init()

    def _init(self):
        #initialize the Hamiltonian
        I = self.I; J=self.J
        NI = Ang._getN(I)
        NJ = Ang._getN(J)
        eyeI = np.eye(NI); eyeJ = np.eye(NJ)
        eye = np.kron(eyeI, eyeJ)

        IdotJ = np.kron(Ang.AngX(I), Ang.AngX(J)) + np.kron(Ang.AngY(I), Ang.AngY(J))+np.kron(Ang.AngZ(I), Ang.AngZ(J))
        if J != .5 and I != .5 and self.B != 0:
            Bterm = self.B*(3*np.dot(IdotJ,IdotJ)+ 1.5*IdotJ  - eye*I*(I+1)*J*(J+1))/(2*I*(2*I-1)*J*(2*J-1))
        else:
            Bterm = 0*IdotJ

        self.HF = self.A*IdotJ + Bterm #hyperfine Hamiltonian / h

        mJ = np.kron(eyeI, Ang.AngZ(J))
        mI = np.kron(Ang.AngZ(I), eyeJ)
        self.mF_operator = mI+mJ
        self.mu = self.muB*(self.gJ*mJ + self.gI*mI) #total magnetic moment operator / h

        if J == int(J):
            Jlabel = "$_{%d}$" % int(J)
        else:
            Jlabel = "$_{%d/%d}$" % (int(2*J), 2)
        self.term = self.L_label + Jlabel
        
        #Calculate allowed values of F
        self.F_vals = [F for F in np.arange(abs(I-J),I+J+1)]
    
        #Calculate hyperfine g factors gF
        self.gF_vals = {}
        for F in self.F_vals:
            self.gF_vals[F] = landeG(F,I,J,self.gI,self.gJ)

    def eig(self,Bz):
        """
        find the energies at field Bz
        returns vals, vecs
        """
        H = self.HF + self.mu*Bz
        vals, vecs = np.linalg.eig(H)
        #sort:
        vecs =  np.array([x for (y,x) in sorted(zip(vals,np.transpose(vecs)), key=lambda pair: pair[0].real)])
        vals = np.array(sorted(vals, key=lambda x: x.real))

        return vals.real, vecs

    def S_MatrixElements(self, ind1, ind2, Bz, return_eig = False):
        #Calculates matrix elements:
        #<ind1 |S_x| ind2>, <ind1 |S_y| ind2>, <ind1 |S_z| ind2>
        #Omitting a factor of hbar

        vals, vecs = self.eig(Bz)
        Sx = np.kron(Ang.Id(self.I),Ang.AngX(self.J)) #electron spin in x direction
        Sy = np.kron(Ang.Id(self.I),Ang.AngY(self.J)) #electron spin in y direction
        Sz = np.kron(Ang.Id(self.I),Ang.AngZ(self.J)) #electron spin in z direction            
        
        Sx_el = vecs[ind1].dot(Sx).dot(vecs[ind2])
        Sy_el = vecs[ind1].dot(Sy).dot(vecs[ind2])        
        Sz_el = vecs[ind1].dot(Sz).dot(vecs[ind2])
        if return_eig:
            return Sx_el, Sy_el, Sz_el, vals, vecs
        else:
            return Sx_el, Sy_el, Sz_el

    #def 
    
    def stateString(self, vec):
        """
        Construct a string representing the given vector in Dirac notation
        """
        s=""
        I=self.I
        J=self.J
        for ind_mI in range(int(2*I+1)):
            for ind_mJ in range(int(2*J+1)):
                mI=ind_mI-I
                mJ=ind_mJ-J
                basisState=np.kron(Ang.AngState(I,mI),Ang.AngState(J,mJ))
                basisStateDag = np.conj(basisState)
                coef=basisStateDag.dot(vec)
                if coef:
                    if s != "":
                        s+=" + "
                    s+=str(coef)+"|mI="+str(mI)+", mJ="+str(mJ)+">"
        return s


# End of AtomIJ class

#### Utility functions ####

def plotEnergies(atom, Bmax,Bmin=0,show=True):
    Bfields = np.linspace(Bmin,Bmax,1001)
    energies = []
    for Bz in Bfields:
        vals, vecs = atom.eig(Bz)
        energies.append(sorted(vals))
    energies_plot = np.transpose(np.array(energies))
    #Plot
    plt.figure()
    for ind,line in reversed(list(enumerate(energies_plot))):
        plt.plot(Bfields/units.G, line/units.MHz,label=ind)
    plt.xlabel("B (G)")
    plt.ylabel("E/h (MHz)")
    plt.title(atom.name+" "+atom.term+" Energies")
    plt.legend()
    plt.tight_layout()
    if show:
        plt.show()

def findDiffFreq(atom_gnd, atom_exc, ind_g, ind_e, Bz):
    """
    Find the transition frequency between ground and excited states
    Gives the answer relative to the difference between the lines' centers of mass
    i.e. don't make any correction for overall energy offset (relevant for transitions between manifolds)
    """
    vals_g, vecs_g = atom_gnd.eig(Bz)
    vals_g = sorted(vals_g)

    vals_e, vecs_e = atom_exc.eig(Bz)
    vals_e = sorted(vals_e)
    return vals_e[ind_e] - vals_g[ind_g]


def plotImgFrequencies(atom_gnd, atom_exc, ind_g, ind_e, ref, Bmax, NB=1001, filepath=None, divide=1):
    """
    Plot transition frequency relative to a reference transition frequency (ref)
    ref is the "bare" frequency, i.e. that comes out of findDiffFreq()
    """
    Bfields = np.linspace(0,Bmax,NB)
    energies = []
    for Bz in Bfields:
        dE = findDiffFreq(atom_gnd, atom_exc, ind_g, ind_e, Bz)
        energies.append(dE)

    energies_plot = (np.array(energies) - ref)/divide

    if filepath:
        save_array = np.transpose(np.array([Bfields/units.G, energies_plot/units.MHz]))
        np.savetxt(filepath, save_array, delimiter=",",header="B(G), ImgFreq (MHz)/%d" % (divide), comments="")

    #Plot
    plt.figure()
    plt.plot(Bfields/units.G, energies_plot/units.MHz)
    plt.xlabel("B (G)")
    plt.ylabel("$\\Delta E$ (MHz)")
    plt.title("Image Freq/%d" % divide)
    plt.tight_layout()
    plt.show()

def imgProgram(gnd, exc, ref_g, ref_e, ref_offset_freq, ind_g, ind_e, Bz, filepath=None, Bmax = 1000*units.G, Bref=0,doplot=True,divide=1,Bmin=0,NB=1001):

    nucyc = findDiffFreq(gnd, exc, ind_g=ref_g, ind_e=ref_e, Bz=Bref)
    ref = nucyc + ref_offset_freq

    dE = findDiffFreq(gnd, exc, ind_g, ind_e, Bz) - ref
    msg = "Image Frequency at %g G:" % (Bz/units.G)
    if divide !=1:
        msg += str(divide) + "x"
    msg += "%g MHz from reference" % (dE/units.MHz/divide)
    print(msg)

    if doplot:
        plotEnergies(gnd, Bmax, Bmin)
        plotEnergies(exc, Bmax, Bmin)
        plotImgFrequencies(gnd, exc, ind_g, ind_e, ref=ref, Bmax=Bmax, filepath=filepath, divide=divide, NB=NB)


if __name__ == "__main__":
    # Sodium ground manifold S_1/2
    L=0;    S= 0.5; J = 0.5; I = 1.5
    gL=0.99997613;     gS=2.0023193043622;     gI = -0.00080461080 # g factors
    gJ = landeG(J, L, S, gL, gS)
    A = 885.81306440*units.MHz
    B = 0
    Na_gnd = AtomIJ("$^{23}$Na","S",I,J,gI, gJ, A, B)

    plotEnergies(Na_gnd, 3000*units.G)
