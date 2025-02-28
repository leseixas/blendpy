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

version = '25.2.8'

import numpy as np
import pandas as pd
from ase.io import read
from ase import Atoms
from ase.optimize import BFGS, BFGSLineSearch, CellAwareBFGS, MDMin, FIRE, FIRE2, GPMin, LBFGS, LBFGSLineSearch, ODE12r, GoodOldQuasiNewton
from ase.filters import UnitCellFilter
from .alloy import Alloy
from .polymorph import Polymorph



class DSIModel(Alloy):
    def __init__(self, alloy_components: list, supercell: list = [1,1,1], calculator = None):
        """
        Initialize the DSIModel class with alloy components, supercell dimensions, and an optional calculator.

        Parameters:
        alloy_components (list): List of alloy components.
        supercell (list, optional): Dimensions of the supercell. (Default: [1, 1, 1])
        calculator (optional): Calculator to attach to each Atoms object. Defaults to None.

        Attributes:
        n_components (int): Number of alloy components.
        supercell (list): Dimensions of the supercell.
        _supercells (list): List to store the supercell Atoms objects.
        dilute_alloys (list): List of dilute alloy configurations.

        Methods:
        banner(): Display the initial banner for blendpy.
        _create_supercells(): Create supercell configurations.
        _create_dilute_alloys(): Create dilute alloy configurations.
        """
        super().__init__(alloy_components)
        self.n_components = len(alloy_components)
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
        """
        Prints a banner with the application name and version.

        The banner includes a stylized text representation of the application name
        and the current version number.

        Returns:
            None
        """
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
        Creates supercells for each alloy component and appends them to the _supercells list.

        This method reads the atomic structure from each file in the alloy_components list,
        creates a supercell by repeating the atomic structure according to the supercell attribute,
        and appends the resulting supercell to the _supercells list.

        Returns:
            None
        """
        for filename in self.alloy_components:
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
        Create a matrix of dilute alloys from the provided supercells.
        This method generates a matrix where each element is a supercell with the 
        first atom's symbol replaced by the first atom's symbol of another supercell.
        The resulting matrix has dimensions n x n, where n is the number of supercells.
        Returns:
            list: A 2D list (matrix) of supercells with diluted alloys.
        Raises:
            ValueError: If there are fewer than two supercells provided.
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


    def optimize(self, method=BFGSLineSearch, fmax: float = 0.01, steps: int = 500, logfile: str = 'optimization.log', mask: list = [1,1,1,1,1,1]):
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
        """
        Computes and returns the energy matrix for the dilute alloys.

        The energy matrix is a square matrix of size `n_components` x `n_components`,
        where each element (i, j) represents the energy of the alloy at position (i, j)
        in the `dilute_alloys` array.

        Returns:
            np.ndarray: A 2D numpy array of shape (n_components, n_components) containing
                        the energy values of the dilute alloys.
        """
        n  = self.n_components
        energy_matrix = np.zeros((n,n), dtype=float)
        for i, row in enumerate(self.dilute_alloys):
            for j, atoms in enumerate(row):
                energy_matrix[i,j] = atoms.info['energy']
        return energy_matrix


    def get_diluting_parameters(self):
        """
        Calculate the diluting parameters for the given dilute alloys.

        This method computes the diluting parameters matrix (m_dsi) for the dilute alloys
        based on the energy differences between the alloys and their components.

        Returns:
            np.ndarray: A 2D numpy array containing the diluting parameters in kJ/mol.

        Raises:
            ValueError: If not all supercells have the same number of atoms.
        """
        number_atoms_list = [ len(atoms) for row in self.dilute_alloys for atoms in row ]
        if len(set(number_atoms_list)) != 1:
            raise ValueError(f"Not all supercells have the same number of atoms: {number_atoms_list}.")
        n  = self.n_components
        x = 1/number_atoms_list[0] # dilution parameter

        m_dsi = np.zeros((n,n), dtype=float)
        energy = self.get_energy_matrix()
        for i, row in enumerate(self.dilute_alloys):
            for j in range(len(row)):
                m_dsi[i,j] = energy[i,j] - ((1-x)*energy[i,i] + x * energy[j,j])
        return m_dsi * (96.4853321233100184) # converting value to kJ/mol


    def get_enthalpy_of_mixing(self, A: int = 0, B: int = 1, npoints: int = 21, slope: list = [0,0]):
        """
        Calculate the enthalpy of mixing for a binary mixture.

        Parameters:
        A (int): Index of the first component in the mixture. (Default: 0)
        B (int): Index of the second component in the mixture. (Default: 1)
        npoints (int): Number of points to calculate along the molar fraction range. (Default: 21)
        slope (list): List containing the slope values for the linear term in the enthalpy calculation. (Default: [0, 0])

        Returns:
        numpy.ndarray: Array of enthalpy values corresponding to the molar fraction range.
        """
        x = np.linspace(0,1,npoints) # molar fraction
        m_dsi = self.get_diluting_parameters()
        enthalpy = m_dsi[A,B] * x * (1-x)**2 + m_dsi[B,A] * x**2 * (1-x) + (1-x) * slope[0] + x * slope[1]
        return enthalpy


    def get_configurational_entropy(self, eps: float = 1.e-4, npoints: int = 21):
        """
        Calculate the configurational entropy of a binary mixture.

        Parameters:
        eps (float): A small value to avoid division by zero in logarithm calculations. (Default: 1.e-4)
        npoints (int): Number of points in the molar fraction range to calculate the entropy. (Default: 21)

        Returns:
        numpy.ndarray: Array of configurational entropy values for the given molar fraction range.

        Notes:
        - The gas constant R is used in kJ/(mol*K).
        - The molar fraction x ranges from 0 to 1.
        """
        R = 8.314/1000
        x = np.linspace(0,1,npoints) # molar fraction
        eps = eps
        entropy = - R * ( (1-x-eps)*np.log(1-(x-eps)) + (x+eps)*np.log(x+eps) )
        return entropy


    def get_spinodal_decomposition(self, A: int = 0, B: int = 1, eps: float = 1.e-4, temperatures = np.arange(600, 2501, 5), npoints: int = 201):
        """
        Calculate the spinodal decomposition curve for a binary mixture.
        Parameters:
        -----------
        A : int, optional
            Index of the first component in the mixture (Default: 0).
        B : int, optional
            Index of the second component in the mixture (Default: 1).
        eps : float, optional
            Small value to avoid division by zero in entropy calculation (Default: 1.e-4).
        temperatures : array-like, optional
            Array of temperatures at which to calculate the spinodal decomposition (Default: np.arange(600, 2501, 5)).
        npoints : int, optional
            Number of points to use in the calculation (Default: 201).
        Returns:
        --------
        df_spinodal : pandas.DataFrame
            DataFrame containing the spinodal decomposition curve with columns "x" (molar fraction) and "t" (temperature).
        Notes:
        ------
        The function calculates the Gibbs free energy as a function of temperature and molar fraction, and then determines
        the spinodal points where the second derivative of the Gibbs free energy with respect to molar fraction changes sign.
        """
        A = A
        B = B
        eps = eps
        temperatures = temperatures
        npoints = npoints
        x = np.linspace(0,1,npoints) # molar fraction

        enthalpy = self.get_enthalpy_of_mixing(A, B, npoints)
        entropy = self.get_configurational_entropy(eps, npoints)

        spinodal = []
        for t in temperatures:
            gibbs = enthalpy - t * entropy
            dx = 1/npoints
            diff_gibbs = np.gradient(gibbs, dx)
            diff2_gibbs = np.gradient(diff_gibbs, dx)
            idx = np.argwhere(np.diff(np.sign(diff2_gibbs - np.zeros(npoints)))).flatten()
            data = [t, x[idx]]
            flattened_array = np.concatenate([np.atleast_1d(item) for item in data])
            spinodal.append(flattened_array)
        
        df0 = pd.DataFrame(spinodal)
        df1 = df0[[1,0]]
        df1.columns = ["x","t"]
        df2 = df0[[2,0]]
        df2.columns = ["x","t"]
        reversed_df2 = df2.iloc[::-1].reset_index(drop=True)
        df_result = pd.concat([df1, reversed_df2], axis=0, ignore_index=True)
        df_spinodal = df_result.dropna()
        return df_spinodal


    # TODO
    def get_binodal_curve(self):
        pass


    # TODO
    def get_phase_diagram(self):
        pass



class BlendpyIO(Atoms):
    def get_energies(self):
        pass

    def get_dsimodel(self):
        pass

    def plot(self):
        pass



class SQS(Atoms):
    def __init__(self):
        pass



class ACE(Atoms):
    def __init__(self):
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
    blendpy.optimize(method=BFGSLineSearch, fmax=0.01, steps=500, logfile=None)

    enthalpy = blendpy.get_enthalpy_of_mixing(A=0, B=1, npoints=101)
    print(enthalpy)

    df_spinodal = blendpy.get_spinodal_decomposition(temperatures = np.arange(300, 3000, 1), npoints = 501)
    df_spinodal.to_csv("../../test/spinodal.csv", index=False, header=True, sep=',')

    # blendpy = Polymorph(alpha='../../test/Pt_fcc.vasp', beta='../../test/Pt_bcc.vasp', calculator = calc_mace)
    # blendpy.optimize()
    # print("Difference between alpha and beta phases:")
    # print(blendpy.get_structural_energy_transition())



