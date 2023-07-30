# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 2023

@author: ariel
"""
from __future__ import division
import AtomBField
from Li6Bfield import gnd
# import AngularMomentum as Ang
# import numericalunits as u
import matplotlib.pyplot as plt
import numpy as np

"""
Calculates RF frequencies exactly and compares to linear prediction
for low field F=1/2-->F=3/2
"""
print(plt.rcParams["font.size"])
plt.rcParams.update({'font.size': 8})

print(plt.rcParams["lines.linewidth"])
plt.rcParams.update({"lines.linewidth":1})

for F in gnd.F_vals:
    print("g_F(F={}): {}".format(F,gnd.gF_vals[F]))

def RF_transition(atom, ind1, ind2, Bz, verbose=False):
    """
    Finds the properties of a given RF transition 
    within the given electronic state
    """
    #Matrix elements and eigenvals/vecs
    Sx, Sy, Sz, vals, vecs = atom.S_MatrixElements(ind1, ind2, Bz,return_eig=True)
    
    #Transition frequency
    nu_MHz = vals[ind2]-vals[ind1]
    
    if verbose:
        print("Sx={}".format(Sx))
        print("Sy={}".format(Sy))
        print("Sz={}".format(Sz))
        print("nu_MHz={}".format(nu_MHz))
    
    return Sx, Sy, Sz, nu_MHz, vals, vecs

"""
Loop over the possible transitions, find their frequency and "strength"
then plot those vs. magnetic field
"""
exc_range={} #allowed F=3/2 states from a given F=1/2 state
#These numbers are for Li-6
gnd_range = [0,1] #indices of the F=1/2 states
exc_range[0] = [3,4,5] #allowed final state indices for changing F
exc_range[1] = [2,3,4]
Bmin=0.01 #Must be non-zero to break degeneracy
Bmax=10.0#Gauss
B_vals=np.linspace(Bmin,Bmax,20)#Gauss
freqs = {}
strengths={}
fmts={0:"-",1:".--"} #plot formats for frequency plot
plt.figure(figsize=(12,5),dpi=120)
plt.subplot(1,3,2)
for ind1 in gnd_range:
    for ind2 in exc_range[ind1]:
        freqs[ind1,ind2]=[]
        strengths[ind1,ind2]=[]
        for Bz in B_vals:
            Sx, Sy, Sz, nu_MHz, vals, vecs = RF_transition(gnd,ind1,ind2,Bz)
            # S_dict = {}
            # S_dict["Sx"]=Sx
            # S_dict["Sy"]=Sy
            # S_dict["Sz"]=Sz
            # if Bz==Bmin:
            #     for key in S_dict.keys():
            #         if S_dict[key] != 0.0:
            #             print(ind1,ind2,": {}={}".format(key,S_dict[key]))
            #             print(ind1,":",gnd.stateString(vecs[ind1]))
            #             print(ind2,":",gnd.stateString(vecs[ind2]))
            strength = np.abs(Sx)**2 + np.abs(Sy)**2 + np.abs(Sz)**2 
            if strength == 0.0:
                nu_MHz = None
            freqs[ind1,ind2].append(nu_MHz)
            strengths[ind1,ind2].append(strength)
        #plt.plot(B_vals, freqs[ind1,ind2],fmts[ind1],label="{}->{}".format(ind1,ind2))
        # print(strengths[ind1,ind2])
        plt.scatter(B_vals, freqs[ind1,ind2],s=50*np.array(strengths[ind1,ind2]),
                    label="{}->{}".format(ind1,ind2))
#first-order prediction
freqs_linear = {}
for mF1 in [-.5, .5]:
    for mF2 in np.arange(mF1-1,mF1+1+1):
        freqs_linear[mF1,mF2]=gnd.A*1.5 + 2./3*gnd.muB*(mF1+mF2)*B_vals
        plt.plot(B_vals, freqs_linear[mF1,mF2],'--',color="gray")
plt.title("Transition Frequencies")
plt.legend()
plt.ylabel("Freq (MHz)")
plt.xlabel("B (G)")

###############################
### Plot transition strengths 
plt.subplot(1,3,3)
for ind1 in gnd_range:
    for ind2 in exc_range[ind1]:
        plt.plot(B_vals, strengths[ind1,ind2],label="{}->{}".format(ind1,ind2))
plt.legend()

plt.title("Transition Strengths")
plt.xlabel("B (G)")
plt.ylabel("strength")

#########################
## Plot the spectrum

plt.subplot(1,3,1)
AtomBField.plotEnergies(gnd, Bmax,Bmin,show=False)
plt.tight_layout()
plt.savefig("Li6RF_LowField_output.png",dpi=120)
plt.show()

