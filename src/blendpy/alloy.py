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

from ase.io import read
from ase.atoms import Atoms

class Alloy(Atoms):
    '''
    A class representing an alloy, inheriting from the Atoms class.
    Methods:
        __init__(alloy_components: list, sublattice_alloy=None):
        _store_chemical_elements():
        get_chemical_elements():
    '''
    def __init__(self, alloy_components: list, sublattice_alloy = None):
        """
        Initialize a new instance of the Alloy class.

        Parameters:
        alloy_components (list): A list of alloy components.
        sublattice_alloy (optional): An optional parameter for sublattice alloy. Default is None.

        Attributes:
        alloy_components (list): Stores the alloy components.
        _chemical_elements (list): Stores the unique chemical elements for each file.
        sublattice_alloy: Stores the sublattice alloy if provided.
        """
        super().__init__(symbols=[], positions=[])
        self.alloy_components = alloy_components
        self._chemical_elements = []  # To store the unique chemical elements for each file
        self._store_chemical_elements()
        self.sublattice_alloy = sublattice_alloy


    def _store_chemical_elements(self):
        """
        Reads chemical elements from files and stores them in the instance variable.
        This method iterates over the filenames stored in `self.alloy_components`, reads the atomic
        structure from each file, extracts the chemical symbols of the atoms, and appends these symbols
        to the instance variable `self._chemical_elements`.
        Note:
            This method assumes that the filenames in `self.alloy_components` can be read using the `read`
            function, which returns an object with a `get_chemical_symbols` method.
        Raises:
            Any exceptions raised by the `read` function or the `get_chemical_symbols` method will propagate
            up to the caller.
        """
        for filename in self.alloy_components:
            atoms = read(filename)
            elements = atoms.get_chemical_symbols()
            self._chemical_elements.append(elements)


    def get_chemical_elements(self):
        """
        Retrieve the set of chemical elements present in the alloy.
        Returns:
            set: A set containing the chemical elements.
        """
        return set(self._chemical_elements)