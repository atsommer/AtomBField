# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 17:00:09 2016

@author: Ariel Sommer

Defines data for Li-6 in B-field
Plots energies, image frequency, saves data
"""

from AtomBField import AtomIJ, landeG, units, imgProgram,findDiffFreq
#import numericalunits as units


#All Lithium:
I = 1; S = 0.5
# g factors
gL=0.99999587; #Gehm, uncited
gI=-0.0004476540 #Arimondo 1977
#gS=2.0023193043737;#from Michael Gehm citing CODATA (bare electron?)
gS=2.0023010 #from Arimondo et al 1977

name = "$^{6}$Li"

#ground state
L_g=0; J_g = 0.5
gJ_g = gS#landeG(J_g, L_g, S, gL, gS) 
A_gnd = 152.1368407*units.MHz
B_gnd = 0 
gnd = AtomIJ(name,"S",I, J_g, gI, gJ_g, A_gnd, B_gnd)

#excited state P_1/2
L_e1=1; J_e1 = 0.5
gJ_e1 = landeG(J_e1, L_e1, S, gL, gS) 
A_D1 = 17.386*units.MHz
exc1 = AtomIJ(name,"P",I, J_e1, gI, gJ_e1, A_D1, B=0)

#excited state P_3/2
L_e2=1; J_e2 = 1.5
gJ_e2 = landeG(J_e2, L_e2, S, gL, gS) 
A_D2 = -1.155*units.MHz
B_D2 = -0.10*units.MHz
exc2 = AtomIJ(name,"P",I, J_e2, gI, gJ_e2, A_D2, B_D2)

#
if __name__=="__main__":
    #reference transition
    ref_g = -1 #highest energy hyperfine state (3/2)
    ref_e = 0 #lowest energy hyperfine state (5/2)
    Bref = 0*units.G #reference is defined at zero field
    
    #beat lock setup
    probe_offset = -200.81*units.MHz #Imaging AOM
    #probe_offset = 0.0*units.MHz #Imaging AOM
    beat_ref_offset = -28*units.MHz #MOT detuning from reference transition
    ref_offset = probe_offset+beat_ref_offset #Total detuning from reference transition
    
    
#    Bz=13.9*units.G#653.5*units.G #field at which to print the image frequency to the console
    Bz=834.2*units.G
    ind_g = 0 #lowest energy ground state at high field
    ind_e = 0#lowest energy excited state at high field
    
    nu_rf = findDiffFreq(gnd, gnd, 0, 1, Bz)
    print("RF freq at %g G: %g MHz" % (Bz/units.G, nu_rf/units.MHz))
    Bmax=1000*units.G    
    
    import os
    folder = os.getcwd()
    filepath = os.path.join(folder,"Li_6_imgfreq.csv")

    imgProgram(gnd, exc2, ref_g, ref_e, ref_offset, ind_g, ind_e, Bz, filepath, Bmax, Bref, doplot=False)
    
