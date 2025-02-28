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
from ase.optimize import BFGS, BFGSLineSearch, CellAwareBFGS, MDMin, FIRE, FIRE2, GPMin, LBFGS, LBFGSLineSearch, ODE12r, GoodOldQuasiNewton
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
    def __init__(self, alloy_basis, supercell, calculator=None):
        """
        Initializes the Blendpy object.
        
        Parameters:
            alloy_basis (list): List of filenames (e.g., POSCAR, extxyz, or CIF).
            supercell (list): Supercell dimensions, e.g., [3, 3, 3].
            calculator (optional): A calculator instance to attach to all Atoms objects.
        """
        super().__init__(alloy_basis, supercell)
        self.supercells = self.get_supercells()
        self.dilute_alloys = self._create_dilute_alloys()

        # If a calculator is provided, attach it to each Atoms object.
        if calculator is not None:
            for atoms in self.supercells:
                atoms.calc = calculator
                energy = atoms.get_potential_energy()
                atoms.info['energy'] = energy
            for atoms in self.dilute_alloys:
                atoms.calc = calculator
                energy = atoms.get_potential_energy()
                atoms.info['energy'] = energy


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
        
        dopant = [atoms.get_chemical_symbols()[-1] for atoms in self._supercells]

        # Iterate over all pairs (i, j) with i != j
        for i in range(n):
            for j in range(n):
                if i != j:
                    # Copy the base supercell from index i.
                    new_atoms = self._supercells[i].copy()
                    # Replace the last atom's symbol with the last symbol from supercell j.
                    new_atoms[-1].symbol = dopant[j]
                    dilute_supercells.append(new_atoms)
        
        return dilute_supercells
    

    def get_energies(self):
        """
        Calculates the total energies for each Atoms object in both the supercells
        and dilute_alloys lists. For each structure, the potential energy is stored
        in the Atoms object's info dictionary with key 'energy'.

        Returns:
            list: A list containing two lists: [energy_supercells, energy_dilute_alloys].
        """
        energy_supercells = []
        for atoms in self.supercells:
            energy = atoms.get_potential_energy()
            atoms.info['energy'] = energy
            energy_supercells.append(energy)

        energy_dilute_alloys = []
        for atoms in self.dilute_alloys:
            energy = atoms.get_potential_energy()
            atoms.info['energy'] = energy
            energy_dilute_alloys.append(energy)

        return [energy_supercells, energy_dilute_alloys]


    def optimize(self, method=BFGSLineSearch, fmax=0.01, steps=500, mask = [1,1,1,1,1,1]):
        """
        Optimizes all Atoms objects in both supercells and dilute_alloys lists.
        
        For each Atoms object, a UnitCellFilter is applied and a BFGS optimizer is created.
        The optimizer is then run with the provided fmax and steps.
        
        Parameters:
            fmax (float): The maximum force criteria.
            steps (int): The maximum number of optimization steps.
        """
        # Optimize original supercells.
        for atoms in self.supercells:
            ucf = UnitCellFilter(atoms, mask=mask)
            optimizer = method(ucf)
            optimizer.run(fmax=fmax, steps=steps)
            energy = atoms.get_potential_energy()
            atoms.info['energy'] = energy

        # Optimize dilute alloy structures.
        for atoms in self.dilute_alloys:
            ucf = UnitCellFilter(atoms, mask=mask)
            optimizer = method(ucf)
            optimizer.run(fmax=fmax, steps=steps)
            energy = atoms.get_potential_energy()
            atoms.info['energy'] = energy


    # TODO
    def get_diluting_parameters(self):
        number_atoms = [len(atoms) for atoms in self.supercells]
        if len(set(number_atoms)) != 1:
            raise ValueError(f"Not all supercells have the same number of atoms: {number_atoms}.")
        
        dilution = 1/number_atoms[0]
        # m12 = 
        # diluting_matrix = 
        return dilution


    # TODO
    def get_entropy(self, method='dsi'):
        pass


# Example usage:
if __name__ == '__main__':
    import warnings
    warnings.filterwarnings("ignore")

    # MACE calculator
    from mace.calculators import mace_mp
    calc_mace = mace_mp(model="small",
                        dispersion=False,
                        default_dtype="float32",
                        device='cpu')

    # Example for binary alloys:
    print("Binary Alloy Example:")
    alloy_files = ['../../test/Au.vasp', '../../test/Pt.vasp']
    supercell = [2, 2, 2]  # This will result in supercells with 8 atoms each.

    blendpy = Blendpy(alloy_files, supercell, calculator=calc_mace)
    
    print("Before optimization:")
    for i, sc in enumerate(blendpy.supercells):
        print(f"Supercell from {alloy_files[i]}: {sc.info['energy']}")
    for da in blendpy.dilute_alloys:
        print("Dilute alloy:", da.info['energy'])
    
    # Optimize all structures.
    blendpy.optimize(method=BFGSLineSearch, fmax=0.01, steps=500)
    
    print("\nAfter optimization:")
    for i, sc in enumerate(blendpy.supercells):
        print(f"Optimized supercell from {alloy_files[i]}: {sc.info['energy']}")
    for da in blendpy.dilute_alloys:
        print("Optimized dilute alloy:", da.info['energy'])


# GPAW calculator
#    from gpaw import GPAW, PW, FermiDirac, Davidson, Mixer
#    calc_gpaw = GPAW(mode=PW(400),
#                     xc='PBE',
#                     kpts=(7,7,7),
#                     occupations=FermiDirac(0.1),
#                     eigensolver=Davidson(5),
#                     spinpol=False,
#                     mixer=Mixer(0.05, 5, 100))




