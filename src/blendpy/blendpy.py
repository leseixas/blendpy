#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# file: blendpy.py

# This code is part of quasigraph.
# MIT License
#
# Copyright (c) 2025 Leandro Seixas Rocha <leandro.fisica@gmail.com> 
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

'''
Module blendpy
'''

import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.optimize import BFGS, MDMin, FIRE
from ase.filters import UnitCellFilter

class Alloy(Atoms):
    def __init__(self, alloy_basis=[], supercell=[1,1,1], suballoy = False):
        super().__init__()
        self.alloy_basis: list = alloy_basis
        self.supercell = supercell
        self.suballoy = suballoy

    def get_alloy_basis(self):
        atoms_alloy_basis = []
        for alloy in self.alloy_basis: 
            atoms = read(alloy)
            atoms_alloy_basis.append(atoms)
        return atoms_alloy_basis

    def get_supercells(self):
        atoms_supercells = []
        for alloy in self.atoms_alloy_basis:
            alloy_sc = alloy.repeat(self.supercell) 
            atoms_supercells.append(alloy_sc)
        return atoms_supercells

class Blendy(Alloy):
    def __init__(self,  method='dsi', calculator = None, optimizer = None):
        super().__init__()
        self.method = method
        self.calculator = calculator
        self.optimizer = optimizer


    # TODO
    def get_enthalpy(self):
        pass

    # TODO
    def get_entropy(self):
        pass



if __name__ == '__main__':
    import warnings
    warnings.filterwarnings("ignore")
#    from gpaw import GPAW, PW, FermiDirac, Davidson, Mixer
#    calc_gpaw = GPAW(mode=PW(400),
#                     xc='PBE',
#                     kpts=(7,7,7),
#                     occupations=FermiDirac(5),
#                     eigensolver=Davidson(5),
#                     spinpol=False,
#                     mixer=Mixer(0.05, 5, 100))

    from mace.calculators import mace_mp
    calc_mace = mace_mp(model="small",
                        dispersion=False,
                        default_dtype="float32",
                        device='cpu')

    alloy = Alloy(alloy_basis=['../../test/Au.vasp', '../../Pt.vasp'],
                    calculator = calc_mace,
                    supercell=[3,3,3])
