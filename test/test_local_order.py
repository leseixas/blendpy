import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import pytest
import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.neighborlist import neighbor_list
from blendpy.constants import R
from blendpy.local_order import LocalOrder


@pytest.fixture
def setup_atoms():
    # Create a fcc bulk structure of Ni
    atoms = bulk('Ni', a=3.52, cubic=True)

    # Change some atoms to Co and Cr
    atoms[0].symbol = 'Co'
    atoms[1].symbol = 'Co'
    atoms[2].symbol = 'Cr'
    return atoms

def test_local_order_initialization(setup_atoms):
    atoms = setup_atoms

    # Initialize LocalOrder with the atoms
    local_order = LocalOrder(atoms)
    
    # Check if the atoms are correctly assigned
    assert local_order.atoms == atoms
    
    # Check if the mask is None
    assert local_order.mask is None
    
    # Check if the number of atoms is correct
    assert local_order.N == len(atoms)
    
    # Check if the concentrations are correctly calculated
    unique, counts = np.unique(local_order.symbols, return_counts=True)
    expected_concentrations = dict(zip(unique, counts / float(local_order.N)))
    assert local_order._concentrations == expected_concentrations


def test_configurational_entropy(setup_atoms):
    atoms = setup_atoms
    
    # Initialize LocalOrder with the atoms
    local_order = LocalOrder(atoms)
    
    # Calculate the configurational entropy
    s_config = local_order.configurational_entropy()
    
    # Calculate the expected configurational entropy manually
    unique, counts = np.unique(local_order.symbols, return_counts=True)
    concentrations = dict(zip(unique, counts / float(local_order.N)))
    expected_s_config = - R * sum(X * np.log(X) for X in concentrations.values() if X > 0)
    
    # Check if the calculated configurational entropy matches the expected value
    assert np.isclose(s_config, expected_s_config, atol=1.e-3), f"Expected {expected_s_config}, but got {s_config}"


def test_nearest_neighbor_correlation(setup_atoms):
    atoms = setup_atoms
    
    # Initialize LocalOrder with the atoms
    local_order = LocalOrder(atoms)
    
    # Calculate the nearest-neighbor correlation entropy correction
    s_corr = local_order.nearest_neighbor_correlation(cutoff=2.5)
    
    # Manually calculate the expected nearest-neighbor correlation entropy correction
    i_list, j_list, _ = neighbor_list('ijd', atoms, 2.5)
    pair_counts = {}
    for i, j in zip(i_list, j_list):
        if i < j:
            pair = tuple(sorted([local_order.symbols[i], local_order.symbols[j]]))
            pair_counts[pair] = pair_counts.get(pair, 0) + 1
    
    total_pairs = sum(pair_counts.values())
    s_corr_manual = 0.0
    for pair, count in pair_counts.items():
        p_ij = count / total_pairs
        prod = local_order._concentrations[pair[0]] * local_order._concentrations[pair[1]]
        if prod > 0 and p_ij > 0:
            s_corr_manual += p_ij * np.log(p_ij / prod)
    
    expected_s_corr = R * total_pairs * s_corr_manual
    
    # Check if the calculated nearest-neighbor correlation entropy correction matches the expected value
    assert np.isclose(s_corr, expected_s_corr, atol=1.e-3), f"Expected {expected_s_corr}, but got {s_corr}"


def test_kikuchi_entropy(setup_atoms):
    atoms = setup_atoms
    
    # Initialize LocalOrder with the atoms
    local_order = LocalOrder(atoms)
    
    # Calculate the Kikuchi entropy
    s_total = local_order.kikuchi_entropy()
    
    # Calculate the expected Kikuchi entropy manually
    s_config = local_order.configurational_entropy()
    s_corr = local_order.nearest_neighbor_correlation()
    expected_s_total = s_config - s_corr
    
    # Check if the calculated Kikuchi entropy matches the expected value
    assert np.isclose(s_total, expected_s_total, atol=1.e-3), f"Expected {expected_s_total}, but got {s_total}"


def test_warren_cowley_parameter(setup_atoms):
    atoms = setup_atoms
    
    # Initialize LocalOrder with the atoms
    local_order = LocalOrder(atoms)
    local_order.cutoff = 3.0  # Set the cutoff distance for neighbor list
    
    # Calculate the Warren–Cowley parameters
    alpha = local_order.warren_cowley_parameter()
    
    # Manually calculate the expected Warren–Cowley parameters
    i_list, j_list, _ = neighbor_list('ijd', atoms, local_order.cutoff)
    ordered_counts = {}
    total_neighbors = {spec: 0 for spec in local_order._concentrations.keys()}
    
    for i, j in zip(i_list, j_list):
        spec_i = local_order.symbols[i]
        spec_j = local_order.symbols[j]
        ordered_counts[(spec_i, spec_j)] = ordered_counts.get((spec_i, spec_j), 0) + 1
        total_neighbors[spec_i] += 1
        ordered_counts[(spec_j, spec_i)] = ordered_counts.get((spec_j, spec_i), 0) + 1
        total_neighbors[spec_j] += 1
    
    expected_alpha = {}
    for (spec_i, spec_j), count in ordered_counts.items():
        if total_neighbors[spec_i] > 0:
            P_ij = count / total_neighbors[spec_i]
            c_j = local_order._concentrations[spec_j]
            expected_alpha[(spec_i, spec_j)] = 1 - (P_ij / c_j)
    
    # Check if the calculated Warren–Cowley parameters match the expected values
    for pair in expected_alpha:
        assert np.isclose(alpha[pair], expected_alpha[pair], atol=1.e-3), f"Expected {expected_alpha[pair]} for pair {pair}, but got {alpha[pair]}"


