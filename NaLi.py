# -*- coding: utf-8 -*-
"""
Created on Thu May 25 17:39:17 2017

@author: BEC3

Parameters for NaLi ground state
"""
from Molecule_ground import units, Molecule
    
S = 1. #electron spin
I1 = 1.5 #1.5 Sodium nucleus
I2 = 1. #1 Lithium nucleus

gS = 2. #electron spin g factor
gI1 = -0.008 # sodium nuclear g factor
gI2 = -0.004 # lithium nuclear g factor

#Hyperfine constants from PRL 119, 143001 (2017)
a1 = 433.20*units.MHz #a_Na
a2 = 74.61*units.MHz #a_Li

def NewNaLi():
    return Molecule(S,I1,I2,gS,gI1,gI2,a1,a2)

NaLi = NewNaLi()
