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

# import os
import numpy as np
# import pandas as pd
from ase.io import read
# from ase.atoms import Atoms
from ase.optimize import BFGS, BFGSLineSearch, LBFGS, LBFGSLineSearch, MDMin, GPMin, FIRE, FIRE2, ODE12r, GoodOldQuasiNewton
from ase.filters import UnitCellFilter
from ase.parallel import parprint as print
from .alloy import Alloy
from .constants import *


class DSIModel(Alloy):
    def __init__(self, alloy_components: list = [], supercell: list = [1,1,1],
                 calculator = None, doping_site: int = 0,
                 x0 = None, energy_matrix = None):
        """
        Initialize the DSI Model.

        Parameters:
            alloy_components (list): List of alloy components.
            supercell (list): Dimensions of the supercell. Default is [1, 1, 1].
            calculator: Calculator object for computing energies. Default is None.
            doping_site (int): Index of the doping site. Default is 0.
            x0: Minimum dilution factor. Default is None.
            energy_matrix: Energy matrix. Default is None.
        Attributes:
            n_components (int): Number of alloy components.
            supercell (list): Dimensions of the supercell.
            _supercells (list): List to store the supercell Atoms objects.
            doping_site (int): Index of the doping site.
            _dilute_alloys (list): List of dilute alloy Atoms objects.
            x0 (float): Minimum dilution factor.
            _energy_matrix: Energy matrix.
            _diluting_parameters: Matrix of diluting parameters.
        """

        print("-----------------------------------------------")
        print("\033[36mDSI Model initialized\033[0m")
        print("-----------------------------------------------")
        super().__init__(alloy_components)
        self.n_components = len(alloy_components)
        print("    Number of components:", self.n_components)
        self.supercell = supercell                          # List: [2,2,2]
        print("    Supercell dimensions:", self.supercell)
        self._supercells = []                               # To store the supercell Atoms objects, like [Atoms("Au32"), Atoms("Pt32")].
        self.doping_site = doping_site
        print("    Doping site:", self.doping_site)
        self._create_supercells()                           # Create supercells from the alloy components
        self._dilute_alloys = self._create_dilute_alloys()  # List [Atoms("Au32"), Atoms("Au31Pt1"), Atoms("Au1Pt31"), Atoms("Pt32")]
        
        if len(self._supercells) > 0:
            n_atoms = len(self._supercells[0])
        else:
            n_atoms = 0
        print("    Number of atoms in the supercell:", n_atoms)

        if x0 is not None:
                self.x0 = x0
        else:
            self.x0 = 1 / n_atoms if n_atoms > 0 else None

        print("    Minimum dilution factor:", self.x0)

        # To store energy_matrix
        self._energy_matrix = energy_matrix

        # To store the diluting parameters matrix
        self._diluting_parameters = None                    # M_DSI matrix

        # If a calculator is provided, attach it to each Atoms object.
        if calculator is not None and len(self._dilute_alloys) > 0:
            for row in self._dilute_alloys:
                for atoms in row:
                    atoms.calc = calculator
                    energy = atoms.get_potential_energy()
                    atoms.info['energy'] = energy
                    print(f"    Total energy ({atoms.get_chemical_formula()}) [Non-relaxed]: {energy} eV")
        
        
    def _create_supercells(self):
        """
        Creates supercells for each alloy component and appends them to the _supercells list.

        This method reads the atomic structure from each file in the alloy_components list,
        creates a supercell by repeating the atomic structure according to the supercell attribute,
        and appends the resulting supercell to the _supercells list.

        Returns:
        --------
            None
        """
        if len(self.alloy_components) > 0:
            for filename in self.alloy_components:
                atoms = read(filename)                          # Read the atomic structure from the file
                supercell_atoms = atoms.repeat(self.supercell)  # Create the supercell using the repeat method
                self._supercells.append(supercell_atoms)

    def get_supercells(self):
        """
        Retrieve the list of supercells.

        Returns:
        --------
            list: A list containing the supercells.
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
        """
        n = len(self._supercells)
        if n < 2:
            return None
        else:
            dopant = [atoms.get_chemical_symbols()[self.doping_site] for atoms in self._supercells]
            print("    Dopant atoms:", dopant)

            list_alloys = []
            # Iterate over all pairs (i, j)
            dilute_supercells_matrix = []
            for i in range(n):
                dilute_matrix_row = []
                for j in range(n):
                    # Copy the base supercell from index i.
                    new_atoms = self._supercells[i].copy()
                    new_atoms[self.doping_site].symbol = dopant[j]
                    list_alloys.append(new_atoms.get_chemical_formula())
                    dilute_matrix_row.append(new_atoms)
                dilute_supercells_matrix.append(dilute_matrix_row)

            print("    Listing dilute alloys:", list_alloys)   # Example: ['Au32', 'Au31Pt1', 'Au1Pt31', 'Pt32']
            return dilute_supercells_matrix


    def optimize(self, method=BFGSLineSearch, fmax: float = 0.01, steps: int = 500, logfile: str = 'optimize.log', mask: list = [1,1,1,1,1,1], verbose: bool = True):
        """
        Optimize the structure of dilute alloys using the specified optimization method.

        Parameters:
            method (Optimizer): The optimization method to use (default is BFGSLineSearch).
            fmax (float): The maximum force criteria for convergence in eV/angstrom (default is 0.01).
            steps (int): The maximum number of optimization steps (default is 500).
            logfile (str): The name of the logfile to store optimization details (default is 'optimize.log').
            mask (list): A list of integers specifying which degrees of freedom to optimize (default is [1, 1, 1, 1, 1, 1]).
            verbose (bool): If True, prints detailed information during optimization (default is True).

        Returns:
        None
        """
        if verbose:
            print("-----------------------------------------------")
            print("\033[36mDilute alloys optimization\033[0m")
            print("-----------------------------------------------")
            print("    Optimization method:", method.__name__)
            print("    Maximum force criteria:", fmax, "eV/ang")
            print("    Maximum number of steps:", steps)
            print("    Logfile:", logfile)
            print("    Mask:", mask)

        for row in self._dilute_alloys:
            for atoms in row:
                ucf = UnitCellFilter(atoms, mask=mask)
                optimizer = method(ucf, logfile=logfile)
                optimizer.run(fmax=fmax, steps=steps)
                if verbose:
                    print(f"    Total energy ({atoms.get_chemical_formula()}) [Relaxed]: {atoms.get_potential_energy()} eV")


    def optimize_nostress(self, fmax: float = 0.01, steps: int = 500, logfile: str = 'optimize.log', verbose: bool = True):
        """
        Optimize the structure of dilute alloys without stress.

        Parameters:
        -----------
        fmax : float, optional
            Maximum force criteria for the optimization in eV/angstrom. Default is 0.01.
        steps : int, optional
            Maximum number of optimization steps. Default is 500.
        logfile : str, optional
            Name of the file where the optimization log will be saved. Default is 'optimize.log'.
        verbose : bool, optional
            If True, prints detailed information about the optimization process. Default is True.

        Returns:
        --------
        None
        """
        if verbose:
            print("-----------------------------------------------")
            print("\033[36mDilute alloys optimization\033[0m")
            print("-----------------------------------------------")
            print("    Maximum force criteria:", fmax, "eV/ang")
            print("    Maximum number of steps:", steps)
            print("    Logfile:", logfile)
        for row in self._dilute_alloys:
            for atoms in row:
                optimizer = LBFGS(atoms, logfile=logfile)
                optimizer.run(fmax=fmax, steps=steps)
                if verbose:
                    print(f"    Total energy ({atoms.get_chemical_formula()}) [Relaxed]: {atoms.get_potential_energy()} eV")


    def set_x0(self, x0: float):
        """
        Set the minimum dilution factor for the model.

        Parameters:
        -----------
        x0 : float
            Minimum dilution factor for the model.

        Raises:
        -------
        ValueError
            If the minimum dilution factor is not a float.
            If the minimum dilution factor is not between 0 and 1.
        """
        if not isinstance(x0, float):
            raise ValueError("The minimum dilution factor must be a float.")
        if x0 < 0 or x0 > 1:
            raise ValueError("The minimum dilution factor must be between 0 and 1.")
        self.x0 = x0


    def get_x0(self):
        """
        Retrieve the minimum dilution factor for the model.

        Returns:
        --------
        float: The minimum dilution factor.
        """
        return self.x0


    def set_energy_matrix(self, energy_matrix: np.ndarray):
        """
        Set the energy matrix for the model.

        Parameters
        ----------
        energy_matrix : np.ndarray
            A 2D numpy array representing the energy matrix. It must be of shape 
            (n_components, n_components) and contain floating point numbers.

        Raises
        ------
        ValueError
            If the energy matrix is not a 2D numpy array.
            If the energy matrix does not contain floating point numbers.
            If the shape of the energy matrix is not (n_components, n_components).
        """
        energy_matrix = np.array(energy_matrix)
        if energy_matrix.ndim != 2:
            raise ValueError("The energy matrix must be a 2D numpy array.")
        if not np.issubdtype(energy_matrix.dtype, np.floating):
            raise ValueError("The energy matrix must be a nd.array of floats.")
        if energy_matrix.shape != (self.n_components, self.n_components):
            raise ValueError("The energy matrix must be a square matrix.")
        self._energy_matrix = energy_matrix

    
    def get_energy_matrix(self, verbose: bool = True) -> np.ndarray:
        """
        Computes and returns the energy matrix for the dilute alloys.

        The energy matrix is a square matrix of size `n_components` x `n_components`,
        where each element (i, j) represents the energy of the alloy at position (i, j)
        in the `dilute_alloys` array.

        Returns:
            np.ndarray: A 2D numpy array of shape (n_components, n_components) containing
                        the energy values of the dilute alloys.
        """
        if self._energy_matrix is not None:
            if verbose:
                print("    Loading energy_matrix...")
                print("    Energy matrix:")
                print(self._energy_matrix)
            return self._energy_matrix
        else:
            if verbose:
                print("    Calculating energy_matrix...")
            energy_matrix = np.zeros((self.n_components, self.n_components), dtype=float)
            for i, row in enumerate(self._dilute_alloys):
                for j, atoms in enumerate(row):
                    if 'energy' not in atoms.info:
                        # print("WARNING: 'energy' is not in atoms.info. Calculating this now in get_energy_matrix method.")
                        atoms.info['energy'] = atoms.get_potential_energy()
                    energy_matrix[i,j] = atoms.info['energy']
            
            # Store energy_matrix as DSIModel attribute
            self._energy_matrix = energy_matrix
            return self._energy_matrix


    def set_diluting_parameters(self, m_dsi: np.ndarray):
        """
        Set the diluting parameters matrix for the model.

        Parameters
        ----------
        m_dsi : np.ndarray
            A 2D numpy array of floats representing the diluting parameters matrix. 
            The shape of the matrix must be (n_components, n_components).

        Raises
        ------
        ValueError
            If the input matrix is not a 2D numpy array.
            If the input matrix does not contain floats.
            If the shape of the input matrix is not (n_components, n_components).
        """
        m_dsi = np.array(m_dsi)
        if m_dsi.ndim != 2:
            raise ValueError("The diluting parameters matrix must be a 2D numpy array.")
        if not np.issubdtype(m_dsi.dtype, np.floating):
            raise ValueError("The diluting parameters matrix must be a nd.array of floats.")
        if m_dsi.shape != (self.n_components, self.n_components):
            raise ValueError("The diluting parameters matrix must be a square matrix.")
        self._diluting_parameters = m_dsi


    def get_diluting_parameters(self, verbose: bool = True) -> np.ndarray:
        """
        Calculate and return the diluting parameters matrix in kJ/mol.
        This method computes the diluting parameters matrix based on the energy matrix
        and stores it as an attribute of the DSIModel instance. If the diluting parameters
        matrix is already computed and stored, it returns the stored matrix.
        Args:
            verbose (bool): If True, prints detailed information about the computation process.
                            Default is True.
        Returns:
            np.ndarray: The diluting parameters matrix in kJ/mol.
        Raises:
            NotImplementedError: If not all supercells have the same number of atoms.
        """

        if verbose:
            print("-----------------------------------------------")
            print("\033[36mDiluting parameters matrix (in kJ/mol)\033[0m")
            print("-----------------------------------------------")

        if self._diluting_parameters is not None:
            return self._diluting_parameters
        else:
            dilute_alloys_flatten = [ atoms for row in self._dilute_alloys for atoms in row]
            number_atoms_list = [ len(atoms) for atoms in dilute_alloys_flatten ]

            if len(set(number_atoms_list)) != 1:
                raise NotImplementedError(f"Not all supercells have the same number of atoms.")
            
            n = self.n_components

            m_dsi = np.zeros((n,n), dtype=float)
            energy = self.get_energy_matrix()
            for i in range(n):
                for j in range(n):
                    m_dsi[i,j] = energy[i,j] - ( (1-self.x0) * energy[i,i] + self.x0 * energy[j,j] )

            m_dsi_kjmol = m_dsi  * convert_eVatom_to_kJmol    # converting value to kJ/mol

            if verbose:
                print("    Energy matrix:")
                print(energy)
                print("    Diluting parameters matrix (in eV):")
                print(m_dsi)
                print("    Diluting parameters matrix (in kJ/mol):")
                print(m_dsi_kjmol)
            
            self._diluting_parameters = m_dsi_kjmol         # Store diluting_parameters as DSIModel attribute
            return m_dsi_kjmol
            

    def get_enthalpy_of_mixing(self, A: int = 0, B: int = 1, slope: list = [0,0], npoints: int = 101) -> np.ndarray:
        """
        Calculate the enthalpy of mixing for a binary mixture.

        Parameters:
        A (int): Index of the first component in the mixture (Default: 0).
        B (int): Index of the second component in the mixture (Default: 1).
        slope (list): List containing the slope values for the linear term in the enthalpy calculation (Default: [0, 0]).
        npoints (int): Number of points to calculate along the molar fraction range (Default: 101).
        
        Returns:
        numpy.ndarray: Array of enthalpy values corresponding to the molar fraction range.
        """
        if not (isinstance(A, int) and isinstance(B, int)):
            raise ValueError("The component indices must be integers.")
        if A >= self.n_components or B >= self.n_components:
            raise ValueError("The component indices must be less than the number of components.")
        if len(slope) != 2: 
            raise ValueError("The slope parameter must have two values.")
        if npoints < 2:
            raise ValueError("The number of points must be greater than 1.")
        
        x = np.linspace(0, 1, npoints)

        m_dsi = self.get_diluting_parameters(verbose=False)

        enthalpy = m_dsi[A,B] * x * (1-x)**2 + m_dsi[B,A] * x**2 * (1-x) + (1-x) * slope[0] + x * slope[1]
        return np.array(enthalpy)


