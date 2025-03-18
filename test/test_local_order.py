import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import pytest
import numpy as np
from ase import Atoms
from blendpy.local_order import LocalOrder


def test_local_order_initialization():
    atoms = Atoms('H2O')
    local_order = LocalOrder(atoms)
    assert local_order._atoms == atoms
    assert local_order._concentrations is None
    assert local_order._pair_correlation_function is None


def test_calculate_concentrations_single_element():
    atoms = Atoms('H2')
    local_order = LocalOrder(atoms)
    local_order.calculate_concentrations()
    expected_concentrations = np.array([1.0])
    np.testing.assert_array_almost_equal(local_order.get_concentrations(), expected_concentrations)


def test_calculate_concentrations_multiple_elements():
    atoms = Atoms('H2O')
    local_order = LocalOrder(atoms)
    local_order.calculate_concentrations()
    expected_concentrations = np.array([2/3, 1/3])
    np.testing.assert_array_almost_equal(local_order.get_concentrations(), expected_concentrations)


def test_calculate_concentrations_complex_structure():
    atoms = Atoms('H2O2')
    local_order = LocalOrder(atoms)
    local_order.calculate_concentrations()
    expected_concentrations = np.array([0.5, 0.5])
    np.testing.assert_array_almost_equal(local_order.get_concentrations(), expected_concentrations)


def test_calculate_concentrations_no_atoms():
    atoms = Atoms('')
    local_order = LocalOrder(atoms)
    local_order.calculate_concentrations()
    expected_concentrations = np.array([])
    np.testing.assert_array_almost_equal(local_order.get_concentrations(), expected_concentrations)

