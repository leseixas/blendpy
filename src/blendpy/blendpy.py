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

version = '25.2.4'

import numpy as np
# import pandas as pd
from ase.io import read
from ase import Atoms
from ase.optimize import BFGS, BFGSLineSearch, CellAwareBFGS, MDMin, FIRE, FIRE2, GPMin, LBFGS, LBFGSLineSearch, ODE12r, GoodOldQuasiNewton
from ase.filters import UnitCellFilter

class Alloy(Atoms):
    def __init__(self, alloy_basis=[], sublattice_alloy = None):
        """
        Initializes the Alloy object.
        
        Parameters:
            alloy_basis (list): A list of filenames (e.g., POSCAR, extxyz, or CIF).
        """
        super().__init__(symbols=[], positions=[])
        self.alloy_basis = alloy_basis
        self._chemical_elements = []  # To store the unique chemical elements for each file
        self._store_chemical_elements()
        self.sublattice_alloy = sublattice_alloy

    def _store_chemical_elements(self):
        """
        For each supercell, retrieve the chemical symbols using the 
        inherited get_chemical_symbols method, convert them to a set 
        to list unique elements, and store them in _chemical_elements.
        """
        for filename in self.alloy_basis:
            atoms = read(filename)
            elements = atoms.get_chemical_symbols()
            self._chemical_elements.append(elements)

    def get_chemical_elements(self):
        """
        Returns the list of unique chemical elements (as sets) for each file.
        """
        return set(self._chemical_elements)


class DSIModel(Alloy):
    def __init__(self, alloy_basis, supercell=[1,1,1], calculator=None):
        """
        Initializes the Dilute Solution Interpolation (DSI) Model object.
        
        Parameters:
            alloy_basis (list): List of filenames (e.g., POSCAR, extxyz, or CIF).
            supercell (list): Supercell dimensions, e.g., [3, 3, 3].
            calculator (optional): A calculator instance to attach to all Atoms objects.
        """

        super().__init__(alloy_basis)
        self.n_components = len(alloy_basis)
        self.supercell = supercell
        self._supercells = []         # To store the supercell Atoms objects
        self._create_supercells()
        self.dilute_alloys = self._create_dilute_alloys()

        # Show blendpy initial banner
        self.banner()

        # If a calculator is provided, attach it to each Atoms object.
        if calculator is not None:
            for row in self.dilute_alloys:
                for atoms in row:
                    atoms.calc = calculator
                    energy = atoms.get_potential_energy()
                    atoms.info['energy'] = energy


    def banner(self):
        print("                                                ")
        print("   _      _                   _                 ")
        print("  | |__  | |  ___  _ __    __| | _ __   _   _   ")
        print("  | '_ \\ | | / _ \\| '_ \\  / _` || '_ \\ | | | |  ")
        print("  | |_) || ||  __/| | | || (_| || |_) || |_| |  ")
        print("  |_.__/ |_| \\___||_| |_| \\__,_|| .__/  \\__, |  ")
        print("                                |_|     |___/   ")
        print("                                                ")
        print(f"                 version: {version}                 ")
        print("                                                ")


    def _create_supercells(self):
        """
        Reads each file in alloy_basis as an ASE Atoms object, applies the repeat (supercell) transformation, and stores the resulting supercell.
        """
        for filename in self.alloy_basis:
            # Read the structure from file (ASE infers file type automatically)
            atoms = read(filename)
            # Create the supercell using the repeat method
            supercell_atoms = atoms.repeat(self.supercell)
            self._supercells.append(supercell_atoms)


    def get_supercells(self):
        """
        Returns the list of supercell ASE Atoms objects.
        """
        return self._supercells
    

    def _create_dilute_alloys(self):
        """
        Creates and returns a list of diluted alloy supercells.
        """
        n = len(self._supercells)
        if n < 2:
            raise ValueError("Need at least two elements to create an alloy.")
        
        dopant = [atoms.get_chemical_symbols()[0] for atoms in self._supercells]

        # Iterate over all pairs (i, j)
        dilute_supercells_matrix = []
        for i in range(n):
            dilute_matrix_row = []
            for j in range(n):
                # Copy the base supercell from index i.
                new_atoms = self._supercells[i].copy()
                # Replace the first atom's symbol with the first symbol from supercell j.
                new_atoms[0].symbol = dopant[j]
                dilute_matrix_row.append(new_atoms)
            dilute_supercells_matrix.append(dilute_matrix_row)
        
        return dilute_supercells_matrix


    def optimize(self, method=BFGSLineSearch, fmax=0.01, steps=500, logfile='optimization.log', mask = [1,1,1,1,1,1]):
        """
        Atoms objects are optimized according to the specified optimization method and parameters.
        
        Parameters:
            method (class): The method to optimize the Atoms object. (Default: BFGSLineSearch)
            fmax (float): The maximum force criteria. (Default: 0.01 eV/ang)
            steps (int): The maximum number of optimization steps. (Default: 500)
            logfile (string): Specifies the file name where the computed optimization forces will be recorded. (Default: 'optimization.log')
            mask (list): A list of directions and angles in Voigt notation that can be optimized.
                         A value of 1 enables optimization, while a value of 0 fixes it. (Default: [1,1,1,1,1,1])
        """

        for row in self.dilute_alloys:
            for atoms in row:
                ucf = UnitCellFilter(atoms, mask=mask)
                optimizer = method(ucf, logfile=logfile)
                optimizer.run(fmax=fmax, steps=steps)
                energy = atoms.get_potential_energy()
                atoms.info['energy'] = energy


    def get_energy_matrix(self):
        n  = self.n_components
        energy_matrix = np.zeros((n,n), dtype=float)
        for i, row in enumerate(self.dilute_alloys):
            for j, atoms in enumerate(row):
                energy_matrix[i,j] = atoms.info['energy']
        return energy_matrix


    def get_diluting_parameters(self):
        number_atoms_list = [ len(atoms) for row in self.dilute_alloys for atoms in row ]
        if len(set(number_atoms_list)) != 1:
            raise ValueError(f"Not all supercells have the same number of atoms: {number_atoms_list}.")
        n  = self.n_components
        x = 1/number_atoms_list[0] # dilution parameter

        m_dsi = np.zeros((n,n), dtype=float)
        energy = self.get_energy_matrix()
        for i, row in enumerate(self.dilute_alloys):
            for j, atoms in enumerate(row):
                m_dsi[i,j] = energy[i,j] - ((1-x)*energy[i,i] + x * energy[j,j])
        return m_dsi * (96.4853321233100184) # converting value to kJ/mol


    def get_enthalpy_of_mixing(self, A=0, B=1, npoints=21):
        x = np.linspace(0,1,npoints) # molar fraction
        m_dsi = self.get_diluting_parameters()
        enthalpy = m_dsi[A,B] * x * (1-x)**2 + m_dsi[B,A] * x**2 * (1-x)
        return enthalpy


    # TODO
    def get_structural_energy_transition(self):
        pass


    # TODO
    def get_entropy_of_mixing(self):
        pass


    # TODO
    def get_spinodal_decomposition(self):
        pass


    # TODO
    def get_phase_diagram(self):
        pass


# Example usage:
if __name__ == '__main__':
    import warnings
    warnings.filterwarnings("ignore")

    ### MACE calculator
    from mace.calculators import mace_mp
    calc_mace = mace_mp(model="small",
                        dispersion=False,
                        default_dtype="float32",
                        device='cpu')

    ### GPAW calculator
    # from gpaw import GPAW, PW, FermiDirac, Davidson, Mixer
    # calc_gpaw = GPAW(mode=PW(500),
    #                  xc='PBE',
    #                  kpts=(7,7,7),
    #                  occupations=FermiDirac(0.1),
    #                  eigensolver=Davidson(5),
    #                  spinpol=False,
    #                  mixer=Mixer(0.05, 5, 100))

    # Example:
    alloy_files = ['../../test/Au.vasp', '../../test/Pt.vasp']
    supercell = [2,2,2]  # This will result in supercells with 8 atoms each.

    blendpy = DSIModel(alloy_files, supercell, calculator=calc_mace)

    # Optimize all structures.
    blendpy.optimize(method=BFGSLineSearch, fmax=0.01, steps=500)
    
    enthalpy = blendpy.get_enthalpy_of_mixing(A=0, B=1, npoints=21)
    print(enthalpy)


