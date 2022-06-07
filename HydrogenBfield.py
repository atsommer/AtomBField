# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 17:00:09 2016

@author: Ariel Sommer

Defines data for Li-6 in B-field
Plots energies, image frequency, saves data
"""

from AtomBField import AtomIJ, plotEnergies,landeG
import numericalunits as units

#Hydrogen
L=0;    S= 0.5; J = 0.5; I = 0.5
gL=1;     gS=2.0;     gI = 0.0 # g factors
gJ = landeG(J, L, S, gL, gS)
A = 1420*units.MHz
B = 0
gnd = AtomIJ("$^{1}$H","S",I,J,gI, gJ, A, B)

#plotEnergies(gnd, 3000*units.G)


#Two spin 1/2 electrons
J = 0.5; I = 0.5
gJ=2.0;     gI = 2.0 # g factors

A = 1420*units.MHz
B = 0
gnd = AtomIJ("$^{1}$H","S",I,J,gI, gJ, A, B)

plotEnergies(gnd, 3000*units.G)

print(landeG(2,1.5,.5,0.,2./3))