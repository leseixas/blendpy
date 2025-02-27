#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# file: blendpy.py

# This code is part of blendpy.
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
    def __init__(self, alloy_basis=[], supercell=[1,1,1], sublattice_alloy = None):
        """
        Initializes the Alloy object.
        
        Parameters:
            alloy_basis (list): A list of filenames (e.g., POSCAR, extxyz, or CIF).
            supercell (list): A list representing the supercell dimensions, e.g., [3, 3, 3].
        """
        super().__init__(symbols=[], positions=[]) # super().__init__()
        self.alloy_basis = alloy_basis
        self.supercell = supercell
        self._supercells = []         # To store the supercell Atoms objects
        self._chemical_elements = []  # To store the unique chemical elements for each file
        self._create_supercells()
        self._store_chemical_elements()
        self.sublattice_alloy = sublattice_alloy

    def _create_supercells(self):
        """
        Reads each file in alloy_basis as an ASE Atoms object,
        applies the repeat (supercell) transformation,
        and stores the resulting supercell.
        """
        for filename in self.alloy_basis:
            # Read the structure from file (ASE infers file type automatically)
            atoms = read(filename)
            # Create the supercell using the repeat method
            supercell_atoms = atoms.repeat(self.supercell)
            self._supercells.append(supercell_atoms)

    def _store_chemical_elements(self):
        """
        For each supercell, retrieve the chemical symbols using the 
        inherited get_chemical_symbols method, convert them to a set 
        to list unique elements, and store them in _chemical_elements.
        """
        for atoms in self._supercells:
            elements = atoms.get_chemical_symbols()
            self._chemical_elements.append(elements)

    def get_supercells(self):
        """
        Returns the list of supercell ASE Atoms objects.
        """
        return self._supercells

    def get_chemical_elements(self):
        """
        Returns the list of unique chemical elements (as sets) for each file.
        """
        return self._chemical_elements


class Blendpy(Alloy):
    def __init__(self, alloy_basis, supercell):
        """
        Initializes the Blendpy object.
        In addition to the Alloy initialization, it computes the dilute alloys.
        
        Parameters:
            alloy_basis (list): List of filenames (e.g., POSCAR, extxyz, or CIF).
            supercell (list): Supercell dimensions, e.g., [3, 3, 3].
        """
        super().__init__(alloy_basis, supercell)
        self.dilute_alloys = self._create_dilute_alloys()

    # def __init__(self, method='dsi', calculator = None, optimizer = None):
    #     super().__init__()
    #     self.method = method
    #     self.calculator = calculator
    #     self.optimizer = optimizer

    def _create_dilute_alloys(self):
        """
        Creates and returns a list of diluted alloy supercells.
        
        For each supercell (base), a copy is made and its last atom's chemical symbol is replaced 
        by the last atom's symbol from every other supercell.
        
        For example, with alloys_basis=['Au.cif', 'Ag.cif', 'Pt.cif'] and supercell=[2,2,2]:
          - From Au8, two alloys are made: Au7Ag1 and Au7Pt1.
          - From Ag8, two alloys: Ag7Au1 and Ag7Pt1.
          - From Pt8, two alloys: Pt7Au1 and Pt7Ag1.
        """
        dilute_supercells = []
        n = len(self._supercells)
        if n < 2:
            raise ValueError("Need at least two supercells to create dilute alloys.")
        
        # Iterate over all pairs (i, j) with i != j
        for i in range(n):
            for j in range(n):
                if i != j:
                    # Copy the base supercell from index i
                    new_atoms = self._supercells[i].copy()
                    # Get the replacement symbol from the last atom of supercell j
                    replacement_symbol = self._supercells[j].get_chemical_symbols()[-1]
                    # Replace the last atom's symbol in the copy
                    new_atoms[-1].symbol = replacement_symbol
                    dilute_supercells.append(new_atoms)
        
        return dilute_supercells

    # TODO
    def get_enthalpy(self):
        pass

    # TODO
    def get_entropy(self):
        pass


# Example usage:
if __name__ == '__main__':
    import warnings
    warnings.filterwarnings("ignore")

    # alloy_files = ['../../test/Au.vasp', '../../test/Pt.vasp']
    # supercell = [2,2,2]

    # Example for binary alloys:
    print("Binary Alloy Example:")
    alloy_files_binary = ['../../test/Au.vasp', '../../test/Pt.vasp']
    supercell = [2, 2, 2]  # This will result in supercells with 8 atoms each.
    
    blendpy_binary = Blendpy(alloy_files_binary, supercell)
    
    print("Original supercells:")
    for i, sc in enumerate(blendpy_binary.get_supercells()):
        print(f"File: {alloy_files_binary[i]} -> {sc}")
    print("\nDilute alloys:")
    for da in blendpy_binary.dilute_alloys:
        print(da)
    
    # Example for ternary alloys:
    print("\nTernary Alloy Example:")
    alloy_files_ternary = ['../../test/Au.vasp', '../../test/Ag.vasp', '../../test/Pt.vasp']
    blendpy_ternary = Blendpy(alloy_files_ternary, supercell)
    
    print("Original supercells:")
    for i, sc in enumerate(blendpy_ternary.get_supercells()):
        print(f"File: {alloy_files_ternary[i]} -> {sc}")
    print("\nDilute alloys:")
    for da in blendpy_ternary.dilute_alloys:
        print(da)


# GPAW calculator
#    from gpaw import GPAW, PW, FermiDirac, Davidson, Mixer
#    calc_gpaw = GPAW(mode=PW(400),
#                     xc='PBE',
#                     kpts=(7,7,7),
#                     occupations=FermiDirac(5),
#                     eigensolver=Davidson(5),
#                     spinpol=False,
#                     mixer=Mixer(0.05, 5, 100))

# MACE calculator
#    from mace.calculators import mace_mp
#    calc_mace = mace_mp(model="small",
#                        dispersion=False,
#                        default_dtype="float32",
#                        device='cpu')


