# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 17:00:09 2016

@author: Ariel Sommer

Defines data for Sodium-23 in B-field
Plots energies, image frequency, saves data
"""
import numpy as np
from AtomBField import AtomIJ, landeG, units,imgProgram,findDiffFreq
#plotImgFrequencies,findDiffFreq, plotEnergies
#import numericalunits as units


#All Sodium:
I = 1.5; S = 0.5
gL=0.99997613;     gS=2.0023193043622;     gI = -0.00080461080 # g factors

#ground state
L_g=0; J_g = 0.5
gJ_g = landeG(J_g, L_g, S, gL, gS) 
A_gnd = 885.81306440*units.MHz
B_gnd = 0 
Na_gnd = AtomIJ("$^{23}$Na","S",I, J_g, gI, gJ_g, A_gnd, B_gnd)

#excited state P_1/2
L_e1=1; J_e1 = 0.5
gJ_e1 = landeG(J_e1, L_e1, S, gL, gS) 
A_D1 = 94.44*units.MHz
Na_exc1 = AtomIJ("$^{23}$Na","P",I, J_e1, gI, gJ_e1, A_D1, B=0)

#excited state P_3/2
L_e2=1; J_e2 = 1.5
gJ_e2 = landeG(J_e2, L_e2, S, gL, gS) 
A_D2 = 18.534*units.MHz
B_D2 = 2.724*units.MHz
Na_exc2 = AtomIJ("$^{23}$Na","P",I, J_e2, gI, gJ_e2, A_D2, B_D2)


if __name__=="__main__":

    #reference transition
    ref_g = -1 #highest energy hyperfine state (3)
    ref_e = -1 #highest energy hyperfine state (2)
    ref_offset = 765*units.MHz #detuning from reference transition before double-pass AOM
    Bref = 0*units.G #reference is defined at zero field

    Bz= 14*units.G #field at which to print the image frequency to the console
    ind_g = 0 #lowest energy ground state at high field
    ind_e = 0#lowest energy excited state at high field
  
    Bmax=1000*units.G    
    
    #RF spectroscopy
    nu_rf = findDiffFreq(Na_gnd, Na_gnd, 0, 1, Bz)
    div=2
    print("RF freq at %g G: %gx%g MHz" % (Bz/units.G,div, nu_rf/units.MHz/div))
    
    #low field RF transitions
    vals, vegs = Na_gnd.eig(13*units.G)
    vals = np.array(sorted(vals))
    for ind in range(4,8):
        print((vals[ind]-vals[ind-1])*units.G)

    import os
    folder = os.getcwd()
    filepath = os.path.join(folder,"Na_23_imgfreq.csv")

    imgProgram(Na_gnd, Na_exc2, ref_g, ref_e, ref_offset, ind_g, ind_e, Bz, filepath, Bmax, Bref,doplot=1, divide=2)
#    
    
    
