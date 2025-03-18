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
from collections import Counter
import numpy as np

class LocalOrder(Atoms):
    def __init__(self, atoms: Atoms):
        """
        Initialize the LocalOrder class with a given set of atoms.

        Parameters
        ----------
        atoms : Atoms
            An instance of the Atoms class representing the atomic structure.
        """
        self._atoms = atoms
        self._concentrations = None
        self._pair_correlation_function = None


    def calculate_concentrations(self):
        """
        Calculate the concentrations of each element in the atomic structure.

        This method counts the occurrences of each chemical symbol in the atomic 
        structure and calculates their concentrations as the ratio of the count 
        of each element to the total number of atoms.

        Attributes:
            _atoms (ase.Atoms): An ASE Atoms object containing the atomic structure.
            _concentrations (np.ndarray): A numpy array storing the calculated 
                                          concentrations of each element.

        Returns:
            None
        """
        elements = Counter(self._atoms.get_chemical_symbols())
        num_atoms = len(self._atoms)
        self._concentrations = np.array([count / num_atoms for count in elements.values()])


    def get_concentrations(self):
        """
        Retrieve the concentrations.

        If the concentrations have already been calculated and stored, return them.
        Otherwise, calculate the concentrations, store them, and then return them.

        Returns:
            dict: The concentrations.
        """
        if self._concentrations is not None:
            return self._concentrations
        else:
            self.calculate_concentrations()
            return self._concentrations
        
    def pair_correlation_function(self):
        pass

    def get_kikuchi_entropy(self):
        pass

    def get_chemical_short_range_order(self):
        pass

    def get_chemical_long_range_order(self):
        pass

    def get_chemical_ordering_energy(self):
        pass