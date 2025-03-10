import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import pytest
import numpy as np
import pandas as pd
from blendpy.phase_diagram import PhaseDiagram


# INITIALIZATION
def test_phase_diagram_init_valid_input():
    """
    Test the initialization of the PhaseDiagram class with valid input arrays.
    This test verifies that the PhaseDiagram object is correctly initialized
    with given enthalpy, entropy, and temperature arrays. It checks that the
    attributes of the PhaseDiagram object match the input arrays and that
    additional properties such as the number of points and the dimensionality
    of the enthalpy array are correctly set.
    Assertions:
        - The enthalpy attribute of the PhaseDiagram object matches the input enthalpy array.
        - The entropy attribute of the PhaseDiagram object matches the input entropy array.
        - The temperatures attribute of the PhaseDiagram object matches the input temperatures array.
        - The npoints attribute of the PhaseDiagram object equals the length of the enthalpy array.
        - The enthalpy attribute is a 1-dimensional array.
        - The x attribute of the PhaseDiagram object is a linearly spaced array from 0 to 1 with the same length as the enthalpy array.
    """
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
    """
    Test the initialization of the PhaseDiagram class with default temperatures.
    This test verifies that the PhaseDiagram object is correctly initialized with
    the provided enthalpy and entropy arrays, and that the default temperature
    range is set correctly. It also checks that the number of points and the x
    values are correctly calculated.
    Assertions:
        - The enthalpy array in the PhaseDiagram object matches the input enthalpy array.
        - The entropy array in the PhaseDiagram object matches the input entropy array.
        - The temperatures array in the PhaseDiagram object matches the default range
          from 300 to 3000 with a step of 50.
        - The number of points (npoints) in the PhaseDiagram object matches the length
          of the input enthalpy array.
        - The x values in the PhaseDiagram object are correctly calculated as a
          linearly spaced array from 0 to 1 with the same length as the input enthalpy array.
    """
    enthalpy = np.array([1, 2, 3])
    entropy = np.array([0.1, 0.2, 0.3])
    
    phase_diagram = PhaseDiagram(enthalpy, entropy)
    
    assert np.array_equal(phase_diagram.enthalpy, enthalpy)
    assert np.array_equal(phase_diagram.entropy, entropy)
    assert np.array_equal(phase_diagram.temperatures, np.arange(300, 3001, 50))
    assert phase_diagram.npoints == len(enthalpy)
    assert np.array_equal(phase_diagram.x, np.linspace(0, 1, len(enthalpy)))


def test_phase_diagram_init_invalid_enthalpy_type():
    """
    Test the initialization of PhaseDiagram with invalid enthalpy type.
    This test checks that a TypeError is raised when the enthalpy array contains
    non-numeric values. The enthalpy array is initialized with string values,
    which should trigger the error. The error message should match the expected
    message indicating that the enthalpy should be a numeric array.
    Raises:
        TypeError: If the enthalpy array contains non-numeric values.
    """
    enthalpy = np.array(['a', 'b', 'c'])
    entropy = np.array([0.1, 0.2, 0.3])
    temperatures = np.array([300, 400, 500])
    
    with pytest.raises(TypeError, match="Enthalpy should be a numeric array."):
        PhaseDiagram(enthalpy, entropy, temperatures)


def test_phase_diagram_init_invalid_entropy_type():
    """
    Test the initialization of the PhaseDiagram class with an invalid entropy type.
    This test checks that a TypeError is raised when the entropy array contains
    non-numeric values. The error message should indicate that the entropy should
    be a numeric array.
    Raises:
        TypeError: If the entropy array contains non-numeric values.
    """
    enthalpy = np.array([1, 2, 3])
    entropy = np.array(['a', 'b', 'c'])
    temperatures = np.array([300, 400, 500])
    
    with pytest.raises(TypeError, match="Entropy should be a numeric array."):
        PhaseDiagram(enthalpy, entropy, temperatures)


def test_phase_diagram_init_invalid_temperatures_type():
    """
    Test the initialization of the PhaseDiagram class with invalid temperature types.
    This test checks that a TypeError is raised when the temperatures array contains
    non-numeric values. The temperatures array in this test contains string values,
    which should trigger the error.
    Raises:
        TypeError: If the temperatures array contains non-numeric values.
    """
    enthalpy = np.array([1, 2, 3])
    entropy = np.array([0.1, 0.2, 0.3])
    temperatures = np.array(['a', 'b', 'c'])
    
    with pytest.raises(TypeError, match="Temperatures should be a numeric array."):
        PhaseDiagram(enthalpy, entropy, temperatures)


# GIBBS FREE ENERGY
def test_get_gibbs_free_energy():
    """
    Test the get_gibbs_free_energy method of the PhaseDiagram class.
    This test verifies that the Gibbs free energy is correctly calculated
    using the provided enthalpy, entropy, and temperature arrays. It checks
    the following:
    - The calculated Gibbs free energy matches the expected values.
    - The gibbs_free_energy attribute of the PhaseDiagram instance is correctly set.
    - The gibbs_free_energy attribute has the correct number of dimensions (2).
    - The shape of the gibbs_free_energy attribute matches the expected shape.
    The Gibbs free energy is calculated using the formula:
    G = H - T * S
    where:
    - G is the Gibbs free energy
    - H is the enthalpy
    - T is the temperature
    - S is the entropy
    """
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
    """
    Test the `get_gibbs_free_energy` method of the `PhaseDiagram` class with default temperatures.
    This test checks the following:
    - The Gibbs free energy is correctly calculated using the provided enthalpy and entropy arrays.
    - The Gibbs free energy is calculated for a range of temperatures from 300K to 3000K in steps of 50K.
    - The calculated Gibbs free energy matches the expected values.
    - The Gibbs free energy attribute of the `PhaseDiagram` instance is correctly set.
    - The Gibbs free energy array has the correct dimensions and shape.
    The Gibbs free energy is calculated using the formula:
    G = H - T * S
    where:
    - G is the Gibbs free energy
    - H is the enthalpy
    - T is the temperature
    - S is the entropy
    """
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


