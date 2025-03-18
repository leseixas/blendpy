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

# from warnings import deprecated
import numpy as np
from ase.io import read
from ase.atoms import Atoms
from .constants import R, convert_eVatom_to_kJmol
from ase.optimize import BFGS, BFGSLineSearch, LBFGS, LBFGSLineSearch, MDMin, GPMin, FIRE, FIRE2, ODE12r, GoodOldQuasiNewton
from ase.filters import UnitCellFilter
from ase.parallel import parprint as print


class Alloy(Atoms):
    """
    The Alloy class represents a collection of atomic structures (alloy components) and provides methods to manipulate and analyze them.
        n_components (int): The number of alloy components.
        _alloy_atoms (list): Stores the atomic structures read from alloy components.
    Methods:
        __init__(alloy_components: list, calculator=None):
        _store_from_atoms():
        get_chemical_elements():
        get_energies():
        optimize(method=BFGSLineSearch, fmax: float = 0.01, steps: int = 500, logfile=None, mask: list = [1,1,1,1,1,1]):
        get_structural_energy_transition(method=BFGSLineSearch, fmax: float = 0.01, steps: int = 500, logfile=None, mask: list = [1,1,1,1,1,1]):
            Calculate the energy difference per atom between two structures after optimizing their structures.
        get_configurational_entropy(eps: float = 1.e-4, npoints: int = 101):
    """
    def __init__(self, alloy_components: list = [], calculator = None):
        """
        Initialize an Alloy object with given alloy components and an optional calculator.

        Parameters:
        alloy_components (list): A list of file names representing the alloy components.
                                 Example: ['Au4.cif', 'Pd4.cif']
        calculator (optional): A calculator object to be attached to each Atoms object.
                               Default is None.

        Attributes:
        alloy_components (list): Stores the provided alloy components.
        n_components (int): Number of alloy components.
        _chemical_elements (list): List of lists containing the chemical elements of each component.
                                   Example: [['Au', 'Au', 'Au', 'Au'], ['Pd', 'Pd', 'Pd', 'Pd']]
        _alloy_atoms (list): List of Atoms objects created from the alloy components.
                             Example: [Atoms('Au4'), Atoms('Pd4')]

        Methods:
        _store_from_atoms: Stores information from the Atoms objects.
        """
        super().__init__(symbols=[], positions=[])
        self.alloy_components = alloy_components           # Example: ['Au4.cif', 'Pd4.cif']
        self.n_components = len(self.alloy_components)     # Example: 2
        self._chemical_elements = []                       # Example: [['Au', 'Au', 'Au', 'Au'], ['Pd', 'Pd', 'Pd', 'Pd']]
        self._alloy_atoms = []                             # Example: [Atoms('Au4'), Atoms('Pd4')]
        self._store_from_atoms()
        self._simplex_space = None

        # If a calculator is provided, attach it to each Atoms object.
        if calculator is not None:
            for atoms in self._alloy_atoms:
                atoms.calc = calculator


    def _store_from_atoms(self):
        """
        Reads atomic data from files specified in `self.alloy_components` and stores 
        the atomic structures and their chemical elements.

        This method reads atomic structures from files listed in `self.alloy_components`, 
        appends the atomic structures to `list_atoms`, and appends the chemical symbols 
        of the atoms to `list_elements`. The results are then stored in the instance 
        variables `_alloy_atoms` and `_chemical_elements`.

        Returns:
            None
        """
        list_atoms = []
        list_elements = []
        if len(self.alloy_components) > 0:
            for filename in self.alloy_components:
                if isinstance(filename, str) == False:
                    raise ValueError("The alloy components must be strings.")
                atoms = read(filename)
                list_atoms.append(atoms)
                list_elements.append(atoms.get_chemical_symbols())
        self._alloy_atoms = list_atoms
        self._chemical_elements = list_elements


    def get_chemical_elements(self) -> list:
        """
        Retrieve the list of chemical elements in the alloy.
        Returns:
            list: A list containing the chemical elements present in the alloy.
        """
        return self._chemical_elements                     # Example: [['Au', 'Au', 'Au', 'Au'], ['Pd', 'Pd', 'Pd', 'Pd']]
    
    
    def get_energies(self, verbose: bool = False) -> list:
        """
        Calculate and return the potential energies of all alloy atoms.

        This method iterates over all atoms in the alloy, calculates their potential energy,
        stores the energy in the atoms' info dictionary, prints the energy along with the 
        chemical formula of the atoms, and appends the energy to a list.

        Returns:
            list: A list of potential energies for each atom in the alloy.
        """
        energies = []
        if len(self._alloy_atoms) > 0:
            for atoms in self._alloy_atoms:
                energy = atoms.get_potential_energy()
                atoms.info['energy'] = energy
                if verbose:
                    print(f"    Total energy ({atoms.get_chemical_formula()}) [Non-relaxed]: {energy} eV")
                energies.append(energy)
        return energies                                     # Example: [-12.4, -10.2]


    def optimize(self, method=BFGSLineSearch, fmax: float = 0.01,
                 steps: int = 500, logfile = None,
                 mask: list = [1,1,1,1,1,1], verbose: bool = False):
        """
        Optimize the atomic structure using the specified optimization method.

        Parameters:
        method (class, optional): The optimization method to use. Default is BFGSLineSearch.
                                  Must be one of [BFGS, BFGSLineSearch, LBFGS, LBFGSLineSearch, MDMin, GPMin, FIRE, FIRE2, ODE12r, GoodOldQuasiNewton].
        fmax (float, optional): The maximum force convergence criterion. Default is 0.01.
        steps (int, optional): The maximum number of optimization steps. Default is 500.
        logfile (str, optional): The file to log the optimization process. Default is None.
        mask (list, optional): A list of six integers (0 or 1) specifying which degrees of freedom to optimize. Default is [1, 1, 1, 1, 1, 1].

        Raises:
        ValueError: If an invalid optimization method is provided.
        ValueError: If fmax is not a float.
        ValueError: If steps is not an integer.
        ValueError: If mask is not a list.
        ValueError: If mask does not have 6 elements.
        ValueError: If any element in mask is not an integer.
        ValueError: If any element in mask is not 0 or 1.
        ValueError: If logfile is not a string.

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
            if verbose:
                print(f"    Total energy ({atoms.get_chemical_formula()}) [Relaxed]: {atoms.get_potential_energy()} eV")


    def get_structural_energy_transition(self, method=BFGSLineSearch,
                                         fmax: float = 0.01,steps: int = 500,
                                         logfile = None, mask: list = [1,1,1,1,1,1]) -> float:
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


    def create_simplex_space(self, resolution: int = 50):
        """
        Create a simplex space for the compositional grid of an alloy.
        This method generates a grid of points in an N-dimensional simplex space,
        where N is the number of components in the alloy. The points in the simplex
        space are normalized such that they sum to 1.
        Parameters:
        resolution (int): The resolution of the grid. Must be a positive integer.
        Raises:
        ValueError: If resolution is not an integer.
        ValueError: If resolution is less than or equal to zero.
        ValueError: If the number of components (self.n_components) < 2.
        Returns:
        None: The generated simplex space is stored in the instance variable
        self._simplex_space.
        """

        if isinstance(resolution, int) == False:
            raise ValueError("resolution must be an integer.")
        if resolution <= 0:
            raise ValueError("resolution must be greater than zero.")
        if self.n_components < 2:
            raise ValueError("The compositional space grid cannot be created for n_components less than 2.")

        # Create an N-dimensional grid of indices from 0 to resolution
        grid_shape = (resolution + 1,) * self.n_components
        indices = np.indices(grid_shape).reshape(self.n_components, -1).T  # shape: (total_points, n_components)
        
        # Filter the indices that sum to resolution
        mask = np.sum(indices, axis=1) == resolution
        simplex_indices = indices[mask]
        
        # Normalize the indices to get points on the simplex (they sum to 1)
        simplex_points = simplex_indices / resolution

        # Store the simplex space in the instance variable
        self._simplex_space = simplex_points


    def get_simplex_space(self, resolution: int = 50):
        """
        Get the compositional space grid for the alloy.

        Parameters:
        resolution (int, optional): The resolution of the grid. Default is 50.

        Returns:
        np.ndarray: A grid of compositional space points.
        """
        if self._simplex_space is None:
            self.create_simplex_space(resolution)
        return self._simplex_space


    def get_ideal_configurational_entropy(self, resolution: int = 50,
                                          eps: float = 1.e-3) -> np.ndarray:
        """
        Calculate the ideal configurational entropy of the system.

        Parameters:
        resolution (int): The resolution of the simplex space. Default is 50.
        eps (float): A small value to avoid log(0). Default is 1.e-3.

        Returns:
        np.ndarray: The calculated ideal configurational entropy.
        """
        if eps <= 0:
            raise ValueError("eps must be greater than zero.")

        x = self.get_simplex_space(resolution)
        entropy = -R * np.sum( (x+eps) * np.log(x + eps), axis=-1 )
        return entropy

