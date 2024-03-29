# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 17:00:09 2016

@author: Ariel Sommer

Quick calculation of transitions frequencies
=== Use this file in the lab ====
Plots energies, image frequency, saves data
"""

from AtomBField import units, imgProgram, findDiffFreq
from Li6Bfield import gnd, exc1, exc2
import os


if __name__=="__main__":
    
    Bz= 286.05*units.G #field at which to print the image frequency to the console
    
    # vals_g, vecs_g = gnd.eig(Bz)
    # vals_g=sorted(vals_g)
    # print(vals_g[0]/units.MHz)

    # RF Spectroscopy
    #nu_rf = findDiffFreq(gnd, gnd, 0, 1, Bz)
    #print("RF freq at %g G: %g MHz" % (Bz/units.G, nu_rf/units.MHz))
    
    #High field imaging
    #reference transition (low-field cycling)
    ref_g = -1 #highest energy hyperfine state (3/2)
    ref_e = 0 #lowest energy hyperfine state (5/2)
    Bref = 0*units.G #reference is defined at zero field
    #beat lock setup
    probe_offset = 0#-200.81*units.MHz #Imaging AOM
    #probe_offset = 0.0*units.MHz #Imaging AOM
    beat_ref_offset =0# -28*units.MHz #MOT detuning from reference transition
    ref_offset = probe_offset+beat_ref_offset #Total detuning from reference transition
    
    Bmax=1000*units.G    
    ind_g = 1 #target ground state at high field
    ind_e = 2-ind_g#target excited state at high field for sigma- transition
    
    
    folder = os.getcwd()
    filepath = os.path.join(folder,"Li_6_imgfreq.csv")

    imgProgram(gnd, exc2, ref_g, ref_e, ref_offset, ind_g, ind_e, Bz, filepath, Bmax, Bref, doplot=1,NB=5001)
    
