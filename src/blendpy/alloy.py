# -*- coding: utf-8 -*-
# file: alloy.py

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
Module alloy
'''

import numpy as np
from ase.io import read
from ase.atoms import Atoms
from .constants import R, convert_eVatom_to_kJmol
from ase.optimize import BFGS, BFGSLineSearch, LBFGS, LBFGSLineSearch, MDMin, GPMin, FIRE, FIRE2, ODE12r, GoodOldQuasiNewton
from ase.filters import UnitCellFilter


class Alloy(Atoms):
    """
    A class to represent an alloy composed of multiple components.
    Methods:
        __init__(alloy_components: list):
        _store_chemical_elements():
        get_chemical_elements():
        get_configurational_entropy(eps: float = 1.e-4, npoints: int = 101):
    """
    def __init__(self, alloy_components: list, calculator = None):
        """
        Initialize a new instance of the Alloy class.

        Parameters:
            alloy_components (list): A list of alloy components.

        Attributes:
            alloy_components (list): Stores the alloy components.
            _chemical_elements (list): Stores the unique chemical elements for each file.
        """
        super().__init__(symbols=[], positions=[])
        self.alloy_components = alloy_components
        self.n_components = len(self.alloy_components)
        self._chemical_elements = []  # To store the chemical elements for each file
        self._alloy_atoms = []
        self._store_from_atoms()

        # If a calculator is provided, attach it to each Atoms object.
        if calculator is not None:
            for atoms in self._alloy_atoms:
                atoms.calc = calculator


    def _store_from_atoms(self):
        """
        Reads atomic structures from `alloy_components` and stores them.
        """
        list_atoms = []
        list_elements = []
        for filename in self.alloy_components:
            atoms = read(filename)
            list_atoms.append(atoms)
            list_elements.append(atoms.get_chemical_symbols())
        self._alloy_atoms = list_atoms
        self._chemical_elements = list_elements


    def get_chemical_elements(self):
        """
        Retrieve the list of chemical elements in the alloy.
        Returns:
            list: A list containing the chemical elements present in the alloy.
        """
        return self._chemical_elements
    
    
    def get_energies(self):
        """
        Calculate and return the potential energies of all alloy atoms.

        This method iterates over all atoms in the alloy, calculates their potential energy,
        stores the energy in the atoms' info dictionary, prints the energy along with the 
        chemical formula of the atoms, and appends the energy to a list.

        Returns:
            list: A list of potential energies for each atom in the alloy.
        """
        energies = []
        for atoms in self._alloy_atoms:
            energy = atoms.get_potential_energy()
            atoms.info['energy'] = energy
            # print(f"    Total energy ({atoms.get_chemical_formula()}) [Non-relaxed]: {energy} eV")
            energies.append(energy)
        return energies


    def optimize(self, method=BFGSLineSearch, fmax: float = 0.01, steps: int = 500, logfile = None, mask: list = [1,1,1,1,1,1]):
        """
        Optimize the atomic structure of the alloy using the specified optimization method.

        Parameters:
        method (class, optional): The optimization method to use. Default is BFGSLineSearch.
        fmax (float, optional): The maximum force convergence criterion. Default is 0.01.
        steps (int, optional): The maximum number of optimization steps. Default is 500.
        logfile (str, optional): The name of the file to log the optimization process. Default is None.
        mask (list, optional): A list indicating which degrees of freedom are allowed to relax. Default is [1, 1, 1, 1, 1, 1].

        Returns:
        None

        Prints:
        The total energy of the relaxed atomic structure for each alloy component.
        """
        if method not in [BFGS, BFGSLineSearch, LBFGS, LBFGSLineSearch, MDMin, GPMin, FIRE, FIRE2, ODE12r, GoodOldQuasiNewton]:
            raise ValueError("Invalid optimization method.")
        if isinstance(fmax, float) == False:
            raise ValueError("fmax must be a float.")
        if isinstance(steps, int) == False:
            raise ValueError("steps must be an integer.")
        if isinstance(mask, list) == False:
            raise ValueError("mask must be a list.")
        if len(mask) != 6:
            raise ValueError("mask must have 6 elements.")
        if all(isinstance(i, int) for i in mask) == False:
            raise ValueError("All elements in mask must be integers.")
        if all(i in [0,1] for i in mask) == False:
            raise ValueError("All elements in mask must be either 0 or 1.")
        if logfile is not None and isinstance(logfile, str) == False:
            raise ValueError("logfile must be a string.")

        for atoms in self._alloy_atoms:
            ucf = UnitCellFilter(atoms, mask=mask)
            optimizer = method(ucf, logfile=logfile)
            optimizer.run(fmax=fmax, steps=steps)
            # print(f"    Total energy ({atoms.get_chemical_formula()}) [Relaxed]: {atoms.get_potential_energy()} eV")


    def get_structural_energy_transition(self, method=BFGSLineSearch, fmax: float = 0.01, steps: int = 500, logfile = None, mask: list = [1,1,1,1,1,1]):
        """
        This method calculates the energy difference per atom between two
        structure after optimizing their structures. The result is converted from
        eV/atom to kJ/mol.
        Returns:
            float: The structural energy transition in kJ/mol.
        Raises:
            ValueError: If the alloy does not have exactly two components.
        """
        if len(self._alloy_atoms) != 2:
            raise ValueError("The alloy must have exactly two components to calculate the structural energy transition.")
        
        self.optimize(method=method, fmax=fmax, steps=steps, logfile=logfile, mask=mask)

        [energy_alpha, energy_beta] = self.get_energies()

        num_atoms_alpha = len(self._alloy_atoms[0])
        num_atoms_beta = len(self._alloy_atoms[1])
        delta_energy = energy_beta/num_atoms_beta - energy_alpha/num_atoms_alpha
        return delta_energy * convert_eVatom_to_kJmol # converting value to kJ/mol


    def get_configurational_entropy(self, eps: float = 1.e-4, npoints: int = 101):
        """
        Calculate the configurational entropy of a binary mixture.

        Parameters:
            eps (float): A small value to avoid division by zero in logarithm calculations (Default: 1.e-4).
            npoints (int): Number of points in the molar fraction range to calculate the entropy (Default: 101).

        Returns:
            numpy.ndarray: Array of configurational entropy values for the given molar fraction range.
        """
        x = np.linspace(0,1,npoints)
        entropy = - R * ( (1-x-eps)*np.log(1-(x-eps)) + (x+eps)*np.log(x+eps) )
        return entropy