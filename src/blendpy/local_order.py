# -*- coding: utf-8 -*-
# file: local_order.py

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

from ase import Atoms
from ase.neighborlist import neighbor_list
from collections import Counter
import numpy as np
from .constants import R


class LocalOrder(Atoms):
    def __init__(self, atoms: Atoms, mask = None):
        """
        Initialize the LocalOrder class with a given set of atoms.
        """
        self.atoms = atoms
        self.mask = mask                                                    # To select only a sub-lattice in the atomic structure
        self.symbols = self.atoms.get_chemical_symbols()                    # Example: ['Co', 'Co', 'Co', ... ]
        self.N = len(self.symbols)                                          # Example: 32
        unique, counts = np.unique(self.symbols, return_counts=True)        # Example: ['Co', 'Cr', 'Ni'], [8, 4, 20]
        
        if self.N > 0:
            self._concentrations = dict(zip(unique, counts / float(self.N)))    # {'Ni': 0.625, 'Co': 0.25, 'Cr': 0.125}


    # def calculate_concentrations(self):
    #     """
    #     Calculate the concentrations of each element in the atomic structure.
    #     """
    #     elements = Counter(self.atoms.get_chemical_symbols())                              # Example: {'Ni': 20, 'Co': 8, 'Cr': 4}                                                      # Example: 32
    #     self._concentrations = np.array([count / self.N for count in elements.values()])   # Example: [0.625, 0.25, 0.125]


    # def get_concentrations(self):
    #     """
    #     Retrieve the concentrations.
    #     """
    #     if self._concentrations is not None:
    #         return self._concentrations
    #     else:
    #         self.calculate_concentrations()
    #         return self._concentrations


    def configurational_entropy(self):
        """
        Calculate the ideal configurational entropy.
        
        S_config = - R * sum_i (X_i * ln(X_i))
        
        Returns:
          S_config_total : float
              The configurational entropy.
        """
        s_config = - R * sum(X * np.log(X) for X in self._concentrations.values() if X > 0)
        return s_config


    def nearest_neighbor_correlation(self, cutoff: float = 5.0):
        """
        Calculate the nearest-neighbor correlation.
        
        Using ASE's neighbor list with the provided cutoff, we first count each unique
        bond (avoiding double counting by considering only pairs with i < j). Then we
        compute the correction as:
        
          S_corr = R * (total_pairs) * sum_{pairs} p_ij * ln(p_ij / (X_i * X_j))
          
        Returns:
          S_corr_total : float
              The nearest-neighbor correlation entropy correction.
        """
        i_list, j_list, _ = neighbor_list('ijd', self.atoms, cutoff)
        
        # Count unique nearest-neighbor pairs
        pair_counts = {}
        for i, j in zip(i_list, j_list):
            if i < j:  # Avoid double counting
                pair = tuple(sorted([self.symbols[i], self.symbols[j]]))
                pair_counts[pair] = pair_counts.get(pair, 0) + 1         # Example: {('Co', 'Cr'): 2, ('Co', 'Ni'): 4, ('Cr', 'Ni'): 4}
        
        total_pairs = sum(pair_counts.values())
        s_corr = 0.0
        for pair, count in pair_counts.items():
            p_ij = count / total_pairs
            prod = self.concentrations[pair[0]] * self.concentrations[pair[1]]      # prod = x_i * x_j
            if prod > 0 and p_ij > 0:
                s_corr += p_ij * np.log(p_ij / prod)
        
        S_corr_total = R * total_pairs * s_corr
        return S_corr_total


    def kikuchi_entropy(self):
        """
        Calculate the (total) Kikuchi entropy as the difference between the ideal configurational
        entropy and the nearest-neighbor correlation entropy.
        
        S_total = S_config - S_corr
        
        Returns:
          S_total : float
              The total Kikuchi entropy.
        """
        S_config = self.configurational_entropy()
        S_corr = self.nearest_neighbor_correlation()
        return S_config - S_corr


    def get_chemical_short_range_order(self):
        pass


    def get_chemical_long_range_order(self):
        pass


    def get_chemical_ordering_energy(self):
        pass