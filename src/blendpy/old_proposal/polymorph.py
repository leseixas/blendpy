# -*- coding: utf-8 -*-
# file: polymorph.py

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

from ase.io import read
from ase.atoms import Atoms
from ase.optimize import BFGSLineSearch
from ase.constraints import UnitCellFilter


class Polymorph(Atoms):
    def __init__(self, alpha: str, beta: str, calculator = None):
        """
        Initialize a Polymorph object.
        Parameters:
        alpha (str): The file path to the alpha polymorph structure.
        beta (str): The file path to the beta polymorph structure.
        calculator: An optional calculator object to be used for calculations. 
                    If not provided, a ValueError will be raised.
        Raises:
        ValueError: If the calculator is not provided.
        Attributes:
        alpha: The alpha polymorph structure read from the provided file path.
        beta: The beta polymorph structure read from the provided file path.
        calculator: The calculator object used for calculations.
        polymorphs: A list containing the alpha and beta polymorph structures.
        """

        super().__init__()
        self.alpha = read(alpha)
        self.beta = read(beta)
        if calculator is None:
            raise ValueError("Polymorph object need to have a calculator.")
        self.calculator = calculator
        self.polymorphs = [self.alpha, self.beta]
        for atoms in self.polymorphs:
            atoms.calc = self.calculator

    
    def get_energies(self):
        """
        Calculate and return the potential energies of polymorphs.
        This method iterates over the polymorphs, calculates the potential energy
        for each set of atoms, stores the energy in the atoms' info dictionary,
        and appends the energy to a list.
        Returns:
            list: A list of potential energies for each set of atoms in polymorphs.
        """

        energies = []
        for atoms in self.polymorphs:
            energy = atoms.get_potential_energy()
            atoms.info['energy'] = energy
            energies.append(energy)
        return energies


    def optimize(self,
                 method=BFGSLineSearch,
                 fmax: float = 0.01,
                 steps: int = 500,
                 logfile: str = 'optimization.log',
                 mask: list = [1,1,1,1,1,1]):
        """
        Atoms objects are optimized according to the specified optimization method and parameters.
        
        Parameters:
            method (class): The method to optimize the Atoms object (Default: BFGSLineSearch).
            fmax (float): The maximum force criteria (Default: 0.01 eV/ang).
            steps (int): The maximum number of optimization steps (Default: 500).
            logfile (string): Specifies the file name where the computed optimization forces will be recorded (Default: 'optimize.log').
            mask (list): A list of directions and angles in Voigt notation that can be optimized.
                         A value of 1 enables optimization, while a value of 0 fixes it. (Default: [1,1,1,1,1,1])
        """
        print("    Optimization method:", method.__name__)
        print("    Maximum force criteria:", fmax, "eV/ang")
        print("    Maximum number of steps:", steps)
        print("    Logfile:", logfile)
        print("    Mask:", mask)

        for atoms in self.polymorphs:
            ucf = UnitCellFilter(atoms, mask=mask)
            optimizer = method(ucf, logfile=logfile)
            optimizer.run(fmax=fmax, steps=steps)
            atoms.info['energy'] = atoms.get_potential_energy()


    def get_structural_energy_transition(self):
        """
        Calculates and returns the difference in structural energy between the alpha and beta phases.
        This method computes the energy difference (in kJ/mol) between the beta and alpha phases of a structure.
        The energy difference is calculated as:
            delta_energy = (energy(beta) / num_atoms(beta)) - (energy(alpha) / num_atoms(alpha))
        The result is then converted to kJ/mol.
        Returns:
            float: The energy difference between the beta and alpha phases in kJ/mol.
        """
        energy_alpha = self.alpha.info['energy']
        energy_beta = self.beta.info['energy']
        num_atoms_alpha = len(self.alpha)
        num_atoms_beta = len(self.beta)
        delta_energy = energy_beta/num_atoms_beta - energy_alpha/num_atoms_alpha
        return delta_energy * (96.4853321233100184) # converting value to kJ/mol
    
