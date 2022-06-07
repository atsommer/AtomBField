# -*- coding: utf-8 -*-
"""
Demo using Li6
Author: Ariel Sommer
"""
from __future__ import division
from Li6Bfield import gnd
from AtomBField import units
from AtomBField import findDiffFreq, plotEnergies
import numpy as np

"""
Calculate RF transitions from the ground state
"""
Bz=832*units.G #DC magnetic field in z direction

#Set the indices of the two states in the transition.
#States are numbered starting from 0 (lowest energy)
ind1 = 0
ind2 = 1

# Calculate the transition frequency
nu_rf = findDiffFreq(gnd, gnd, ind1, ind2, Bz)
print("RF freq at %g G: %g MHz" % (Bz/units.G, nu_rf/units.MHz))

# Calculate the matrix elements for RF acting on the electron spin
S_elements = gnd.S_MatrixElements(ind1,ind2,Bz)
print("<%g|Sx|%g>/hbar=%g" % (ind1, ind2, np.abs(S_elements[0])))

#plot the energy levels
plotEnergies(gnd, 1000*units.G)
