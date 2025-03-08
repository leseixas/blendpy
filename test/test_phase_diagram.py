import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import pytest
import numpy as np
import pandas as pd
from blendpy.phase_diagram import PhaseDiagram

def test_phase_diagram_init_valid_input():
    enthalpy = np.array([1, 2, 3])
    entropy = np.array([0.1, 0.2, 0.3])
    temperatures = np.array([300, 400, 500])
    
    phase_diagram = PhaseDiagram(enthalpy, entropy, temperatures)
    
    assert np.array_equal(phase_diagram.enthalpy, enthalpy)
    assert np.array_equal(phase_diagram.entropy, entropy)
    assert np.array_equal(phase_diagram.temperatures, temperatures)
    assert phase_diagram.npoints == len(enthalpy)
    assert phase_diagram.enthalpy.ndim == 1
    assert np.array_equal(phase_diagram.x, np.linspace(0, 1, len(enthalpy)))


def test_phase_diagram_init_default_temperatures():
    enthalpy = np.array([1, 2, 3])
    entropy = np.array([0.1, 0.2, 0.3])
    
    phase_diagram = PhaseDiagram(enthalpy, entropy)
    
    assert np.array_equal(phase_diagram.enthalpy, enthalpy)
    assert np.array_equal(phase_diagram.entropy, entropy)
    assert np.array_equal(phase_diagram.temperatures, np.arange(300, 3001, 50))
    assert phase_diagram.npoints == len(enthalpy)
    assert np.array_equal(phase_diagram.x, np.linspace(0, 1, len(enthalpy)))


def test_phase_diagram_init_invalid_enthalpy_type():
    enthalpy = np.array(['a', 'b', 'c'])
    entropy = np.array([0.1, 0.2, 0.3])
    temperatures = np.array([300, 400, 500])
    
    with pytest.raises(TypeError, match="Enthalpy should be a numeric array."):
        PhaseDiagram(enthalpy, entropy, temperatures)


def test_phase_diagram_init_invalid_entropy_type():
    enthalpy = np.array([1, 2, 3])
    entropy = np.array(['a', 'b', 'c'])
    temperatures = np.array([300, 400, 500])
    
    with pytest.raises(TypeError, match="Entropy should be a numeric array."):
        PhaseDiagram(enthalpy, entropy, temperatures)


def test_phase_diagram_init_invalid_temperatures_type():
    enthalpy = np.array([1, 2, 3])
    entropy = np.array([0.1, 0.2, 0.3])
    temperatures = np.array(['a', 'b', 'c'])
    
    with pytest.raises(TypeError, match="Temperatures should be a numeric array."):
        PhaseDiagram(enthalpy, entropy, temperatures)


# GIBBS FREE ENERGY
def test_get_gibbs_free_energy():
    enthalpy = np.array([1, 2, 3])
    entropy = np.array([0.1, 0.2, 0.3])
    temperatures = np.array([300, 400, 500])
    
    phase_diagram = PhaseDiagram(enthalpy, entropy, temperatures)
    gibbs_free_energy = phase_diagram.get_gibbs_free_energy()
    
    expected_gibbs_free_energy = enthalpy - temperatures[:, None] * entropy
    assert np.array_equal(gibbs_free_energy, expected_gibbs_free_energy)
    assert np.array_equal(phase_diagram.gibbs_free_energy, expected_gibbs_free_energy)
    assert phase_diagram.gibbs_free_energy.ndim == 2
    assert phase_diagram.gibbs_free_energy.shape == (len(temperatures), len(enthalpy))


def test_get_gibbs_free_energy_default_temperatures():
    enthalpy = np.array([1, 2, 3])
    entropy = np.array([0.1, 0.2, 0.3])
    
    phase_diagram = PhaseDiagram(enthalpy, entropy)
    gibbs_free_energy = phase_diagram.get_gibbs_free_energy()
    
    temperatures = np.arange(300, 3001, 50)
    expected_gibbs_free_energy = enthalpy - temperatures[:, None] * entropy
    assert np.array_equal(gibbs_free_energy, expected_gibbs_free_energy)
    assert np.array_equal(phase_diagram.gibbs_free_energy, expected_gibbs_free_energy)
    assert phase_diagram.gibbs_free_energy.ndim == 2
    assert phase_diagram.gibbs_free_energy.shape == (len(temperatures), len(enthalpy))




def test_phase_diagram_init_valid_input():
    enthalpy = np.array([1, 2, 3])
    entropy = np.array([0.1, 0.2, 0.3])
    temperatures = np.array([300, 400, 500])
    
    phase_diagram = PhaseDiagram(enthalpy, entropy, temperatures)
    
    assert np.array_equal(phase_diagram.enthalpy, enthalpy)
    assert np.array_equal(phase_diagram.entropy, entropy)
    assert np.array_equal(phase_diagram.temperatures, temperatures)
    assert phase_diagram.npoints == len(enthalpy)
    assert phase_diagram.enthalpy.ndim == 1
    assert np.array_equal(phase_diagram.x, np.linspace(0, 1, len(enthalpy)))


def test_phase_diagram_init_default_temperatures():
    enthalpy = np.array([1, 2, 3])
    entropy = np.array([0.1, 0.2, 0.3])
    
    phase_diagram = PhaseDiagram(enthalpy, entropy)
    
    assert np.array_equal(phase_diagram.enthalpy, enthalpy)
    assert np.array_equal(phase_diagram.entropy, entropy)
    assert np.array_equal(phase_diagram.temperatures, np.arange(300, 3001, 50))
    assert phase_diagram.npoints == len(enthalpy)
    assert np.array_equal(phase_diagram.x, np.linspace(0, 1, len(enthalpy)))


def test_phase_diagram_init_invalid_enthalpy_type():
    enthalpy = np.array(['a', 'b', 'c'])
    entropy = np.array([0.1, 0.2, 0.3])
    temperatures = np.array([300, 400, 500])
    
    with pytest.raises(TypeError, match="Enthalpy should be a numeric array."):
        PhaseDiagram(enthalpy, entropy, temperatures)


def test_phase_diagram_init_invalid_entropy_type():
    enthalpy = np.array([1, 2, 3])
    entropy = np.array(['a', 'b', 'c'])
    temperatures = np.array([300, 400, 500])
    
    with pytest.raises(TypeError, match="Entropy should be a numeric array."):
        PhaseDiagram(enthalpy, entropy, temperatures)


def test_phase_diagram_init_invalid_temperatures_type():
    enthalpy = np.array([1, 2, 3])
    entropy = np.array([0.1, 0.2, 0.3])
    temperatures = np.array(['a', 'b', 'c'])
    
    with pytest.raises(TypeError, match="Temperatures should be a numeric array."):
        PhaseDiagram(enthalpy, entropy, temperatures)


# GIBBS FREE ENERGY
def test_get_gibbs_free_energy():
    enthalpy = np.array([1, 2, 3])
    entropy = np.array([0.1, 0.2, 0.3])
    temperatures = np.array([300, 400, 500])
    
    phase_diagram = PhaseDiagram(enthalpy, entropy, temperatures)
    gibbs_free_energy = phase_diagram.get_gibbs_free_energy()
    
    expected_gibbs_free_energy = enthalpy - temperatures[:, None] * entropy
    assert np.array_equal(gibbs_free_energy, expected_gibbs_free_energy)
    assert np.array_equal(phase_diagram.gibbs_free_energy, expected_gibbs_free_energy)
    assert phase_diagram.gibbs_free_energy.ndim == 2
    assert phase_diagram.gibbs_free_energy.shape == (len(temperatures), len(enthalpy))


def test_get_gibbs_free_energy_default_temperatures():
    enthalpy = np.array([1, 2, 3])
    entropy = np.array([0.1, 0.2, 0.3])
    
    phase_diagram = PhaseDiagram(enthalpy, entropy)
    gibbs_free_energy = phase_diagram.get_gibbs_free_energy()
    
    temperatures = np.arange(300, 3001, 50)
    expected_gibbs_free_energy = enthalpy - temperatures[:, None] * entropy
    assert np.array_equal(gibbs_free_energy, expected_gibbs_free_energy)
    assert np.array_equal(phase_diagram.gibbs_free_energy, expected_gibbs_free_energy)
    assert phase_diagram.gibbs_free_energy.ndim == 2
    assert phase_diagram.gibbs_free_energy.shape == (len(temperatures), len(enthalpy))


# SPINODAL DECOMPOSITION


