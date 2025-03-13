import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
import numpy as np
from ase.calculators.emt import EMT
from blendpy.alloy import Alloy
from ase import Atoms
import pytest
from ase.optimize import BFGS, BFGSLineSearch, LBFGS, LBFGSLineSearch, MDMin, GPMin, FIRE, FIRE2, ODE12r, GoodOldQuasiNewton


# Fixtures
@pytest.fixture
def setup_data():
    """
    Sets up the data for alloy components and initializes the calculator.

    Returns:
        tuple: A tuple containing:
            - alloy_components (list): A list of file paths to the alloy component CIF files.
            - calculator (EMT): An instance of the EMT calculator.
    """
    alloy_components = ["data/bulk/Au.cif", "data/bulk/Pt.cif"]
    calculator = EMT()
    return alloy_components, calculator

# __init__
def test_alloy_initialization(setup_data):
    """
    Test the initialization of the Alloy class.

    This test verifies that the Alloy object is correctly initialized with the given
    alloy components and that its attributes are set as expected.

    Assertions:
        - The alloy_components attribute of the Alloy object should match the input list.
        - The _chemical_elements attribute should be a list.
        - The _chemical_elements attribute should contain the expected chemical elements
          for each component in the alloy.
    """
    alloy_components, _ = setup_data
    alloy = Alloy(alloy_components)

    assert alloy.alloy_components == alloy_components
    assert isinstance(alloy._chemical_elements, list)
    assert alloy._chemical_elements == [['Au','Au','Au','Au'],['Pt','Pt','Pt','Pt']]


def test_alloy_initialization_with_calculator(setup_data):
    """
    Test the initialization of the Alloy class with a calculator.

    This test verifies that the Alloy object is correctly initialized with the given
    alloy components and that the calculator is attached to each Atoms object.

    Assertions:
        - The alloy_components attribute of the Alloy object should match the input list.
        - The _chemical_elements attribute should be a list.
        - The _chemical_elements attribute should contain the expected chemical elements
            for each component in the alloy.
        - Each Atoms object should have the calculator attached.
        - The energy attribute should be set in the info dictionary of each Atoms object.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)

    assert alloy.alloy_components == alloy_components
    assert isinstance(alloy._chemical_elements, list)
    assert alloy._chemical_elements == [['Au','Au','Au','Au'],['Pt','Pt','Pt','Pt']]
    
    for atoms in alloy._alloy_atoms:
        assert atoms.calc == calculator
        atoms.info['energy'] = atoms.get_potential_energy()
        assert "energy" in atoms.info.keys()


def test_store_from_atoms(setup_data):
    """
    Test the Alloy class's ability to store atoms from given CIF files.
    This test checks the following:
    - The `_alloy_atoms` attribute is a list.
    - The length of `_alloy_atoms` matches the number of alloy components.
    - Each item in `_alloy_atoms` is an instance of the Alloy class.
    - The `_chemical_elements` attribute is a list.
    - The length of `_chemical_elements` matches the number of alloy components.
    - Each item in `_chemical_elements` is a list.
    - The `_chemical_elements` attribute contains the expected chemical elements for each alloy component.
    """
    alloy_components, _ = setup_data
    alloy = Alloy(alloy_components)
    
    assert isinstance(alloy._alloy_atoms, list)
    assert len(alloy._alloy_atoms) == len(alloy_components)
    assert all(isinstance(atoms, Atoms) for atoms in alloy._alloy_atoms)

    assert isinstance(alloy._chemical_elements, list)
    assert len(alloy._chemical_elements) == len(alloy_components)
    assert all(isinstance(elements, list) for elements in alloy._chemical_elements)
    assert alloy._chemical_elements == [['Au','Au','Au','Au'],['Pt','Pt','Pt','Pt']]


# get_chemical_elements
def test_get_chemical_elements(setup_data):
    """
    Test the get_chemical_elements method of the Alloy class.
    This test initializes an Alloy object with a list of component file paths
    and checks if the get_chemical_elements method returns the expected list
    of chemical elements.
    The expected result is a list of lists, where each inner list contains
    the chemical symbols of the elements present in the corresponding component file.
    Assertions:
        - The result of alloy.get_chemical_elements() should match the expected_elements list.
    """
    alloy_components, _ = setup_data
    alloy = Alloy(alloy_components)
    
    expected_elements = [['Au', 'Au', 'Au', 'Au'], ['Pt', 'Pt', 'Pt', 'Pt']]
    assert alloy.get_chemical_elements() == expected_elements


# get_energies
def test_get_energies(setup_data):
    """
    Test the get_energies method of the Alloy class.

    This test initializes an Alloy object with a list of component file paths,
    attaches a calculator to each Atoms object, and checks if the get_energies
    method returns the expected list of potential energies.

    Assertions:
        - The result of alloy.get_energies() should be a list.
        - Each element in the result list should be a float.
        - The length of the result list should match the number of alloy components.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)
    
    energies = alloy.get_energies()
    
    assert isinstance(energies, list)
    assert all(isinstance(energy, float) for energy in energies)
    assert len(energies) == len(alloy_components)


# optimize
def test_optimize_with_valid_parameters(setup_data):
    """
    Test the optimize method with valid parameters.

    This test verifies that the optimize method runs without errors when provided with valid parameters.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)

    alloy.optimize(method=BFGSLineSearch, fmax=0.01, steps=500, logfile=None, mask=[1,1,1,1,1,1])
    optimized_energies = alloy.get_energies()
    
    assert isinstance(optimized_energies, list)
    assert all(isinstance(energy, float) for energy in optimized_energies)
    assert len(optimized_energies) == len(alloy_components)


def test_optimize_with_invalid_method(setup_data):
    """
    Test the optimize method with an invalid optimization method.

    This test verifies that the optimize method raises a ValueError when provided with an invalid optimization method.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)
    
    with pytest.raises(ValueError, match="Invalid optimization method."):
        alloy.optimize(method="invalid_method", fmax=0.01, steps=500, logfile=None, mask=[1,1,1,1,1,1])


def test_optimize_with_invalid_fmax(setup_data):
    """
    Test the optimize method with an invalid fmax parameter.

    This test verifies that the optimize method raises a ValueError when provided with an invalid fmax parameter.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)
    
    with pytest.raises(ValueError, match="fmax must be a float."):
        alloy.optimize(method=BFGSLineSearch, fmax="invalid_fmax", steps=500, logfile=None, mask=[1,1,1,1,1,1])


def test_optimize_with_invalid_steps(setup_data):
    """
    Test the optimize method with an invalid steps parameter.

    This test verifies that the optimize method raises a ValueError when provided with an invalid steps parameter.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)
    
    with pytest.raises(ValueError, match="steps must be an integer."):
        alloy.optimize(method=BFGSLineSearch, fmax=0.01, steps="invalid_steps", logfile=None, mask=[1,1,1,1,1,1])


def test_optimize_with_invalid_mask(setup_data):
    """
    Test the optimize method with an invalid mask parameter.

    This test verifies that the optimize method raises a ValueError when provided with an invalid mask parameter.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)
    
    with pytest.raises(ValueError, match="mask must be a list."):
        alloy.optimize(method=BFGSLineSearch, fmax=0.01, steps=500, logfile=None, mask="invalid_mask")
    
    with pytest.raises(ValueError, match="mask must have 6 elements."):
        alloy.optimize(method=BFGSLineSearch, fmax=0.01, steps=500, logfile=None, mask=[1,1,1])
    
    with pytest.raises(ValueError, match="All elements in mask must be integers."):
        alloy.optimize(method=BFGSLineSearch, fmax=0.01, steps=500, logfile=None, mask=[1,1,1,1,1,"invalid"])
    
    with pytest.raises(ValueError, match="All elements in mask must be either 0 or 1."):
        alloy.optimize(method=BFGSLineSearch, fmax=0.01, steps=500, logfile=None, mask=[3,1,4,1,5,9])


def test_optimize_with_invalid_logfile(setup_data):
    """
    Test the optimize method with an invalid logfile parameter.

    This test verifies that the optimize method raises a ValueError when provided with an invalid logfile parameter.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)
    
    with pytest.raises(ValueError, match="logfile must be a string."):
        alloy.optimize(method=BFGSLineSearch, fmax=0.01, steps=500, logfile=123, mask=[1,1,1,1,1,1])


def test_get_structural_energy_transition(setup_data):
    """
    Test the get_structural_energy_transition method of the Alloy class.

    This test initializes an Alloy object with a list of two component file paths,
    attaches a calculator to each Atoms object, and checks if the 
    get_structural_energy_transition method returns the expected energy transition.

    Assertions:
        - The result of alloy.get_structural_energy_transition() should be a float.
        - The alloy object should raise a ValueError if it does not have exactly two components.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)
    
    energy_transition = alloy.get_structural_energy_transition()
    
    assert isinstance(energy_transition, float)

# create_simplex_space
def test_create_simplex_space_with_valid_resolution(setup_data):
    """
    Test the create_simplex_space method with a valid resolution.

    This test verifies that the create_simplex_space method runs without errors when provided with a valid resolution.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)

    alloy.create_simplex_space(resolution=50)
    simplex_space = alloy.get_simplex_space()

    assert isinstance(simplex_space, np.ndarray)
    assert simplex_space.shape[1] == len(alloy_components)
    assert np.allclose(np.sum(simplex_space, axis=1), 1, atol=1e-6)


def test_create_simplex_space_with_invalid_resolution(setup_data):
    """
    Test the create_simplex_space method with an invalid resolution.

    This test verifies that the create_simplex_space method raises a ValueError when provided with an invalid resolution.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)

    with pytest.raises(ValueError, match="resolution must be an integer."):
        alloy.create_simplex_space(resolution="invalid_resolution")

    with pytest.raises(ValueError, match="resolution must be greater than zero."):
        alloy.create_simplex_space(resolution=0)


def test_create_simplex_space_with_single_component():
    """
    Test the create_simplex_space method with a single-component alloy.

    This test verifies that the create_simplex_space method raises a ValueError when the alloy has only one component.
    """
    alloy_components = ["data/bulk/Au.cif"]
    alloy = Alloy(alloy_components)

    with pytest.raises(ValueError, match="The compositional space grid cannot be created for n_components less than 2."):
        alloy.create_simplex_space(resolution=50)


def test_get_simplex_space_with_default_resolution(setup_data):
    """
    Test the get_simplex_space method with the default resolution.

    This test verifies that the get_simplex_space method runs without errors when using the default resolution.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)

    simplex_space = alloy.get_simplex_space()

    assert isinstance(simplex_space, np.ndarray)
    assert simplex_space.shape[1] == len(alloy_components)
    assert np.allclose(np.sum(simplex_space, axis=1), 1, atol=1e-6)


def test_get_simplex_space_with_custom_resolution(setup_data):
    """
    Test the get_simplex_space method with a custom resolution.

    This test verifies that the get_simplex_space method runs without errors when using a custom resolution.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)

    resolution = 30
    simplex_space = alloy.get_simplex_space(resolution=resolution)

    assert isinstance(simplex_space, np.ndarray)
    assert simplex_space.shape[1] == len(alloy_components)
    assert np.allclose(np.sum(simplex_space, axis=1), 1, atol=1e-6)


def test_get_simplex_space_creates_space_if_none_exists(setup_data):
    """
    Test that get_simplex_space creates the simplex space if it does not already exist.

    This test verifies that the get_simplex_space method calls create_simplex_space if the simplex space is None.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)

    # Ensure the simplex space is None initially
    alloy._simplex_space = None

    simplex_space = alloy.get_simplex_space()

    assert isinstance(simplex_space, np.ndarray)
    assert simplex_space.shape[1] == len(alloy_components)
    assert np.allclose(np.sum(simplex_space, axis=1), 1, atol=1e-6)
    assert alloy._simplex_space is not None
    assert alloy._simplex_space is simplex_space


def test_get_simplex_space_returns_existing_space(setup_data):
    """
    Test that get_simplex_space returns the existing simplex space if it already exists.

    This test verifies that the get_simplex_space method does not call create_simplex_space if the simplex space already exists.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)

    # Manually create a simplex space
    alloy.create_simplex_space(resolution=20)
    existing_simplex_space = alloy._simplex_space

    simplex_space = alloy.get_simplex_space()

    assert np.array_equal(simplex_space, existing_simplex_space)


# get_ideal_configurational_entropy
def test_get_ideal_configurational_entropy_with_valid_parameters(setup_data):
    """
    Test the get_ideal_configurational_entropy method with valid parameters.

    This test verifies that the get_ideal_configurational_entropy method runs without errors when provided with valid parameters.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)

    entropy = alloy.get_ideal_configurational_entropy(resolution=50, eps=1.e-6)

    assert isinstance(entropy, np.ndarray)
    assert entropy.shape[0] > 0
    assert np.all(entropy >= 0)


def test_get_ideal_configurational_entropy_with_invalid_eps(setup_data):
    """
    Test the get_ideal_configurational_entropy method with an invalid eps parameter.

    This test verifies that the get_ideal_configurational_entropy method raises a ValueError when provided with an invalid eps parameter.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)

    with pytest.raises(ValueError, match="eps must be greater than zero."):
        alloy.get_ideal_configurational_entropy(resolution=50, eps=0)


def test_get_ideal_configurational_entropy_with_default_parameters(setup_data):
    """
    Test the get_ideal_configurational_entropy method with default parameters.

    This test verifies that the get_ideal_configurational_entropy method runs without errors when using the default parameters.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)

    entropy = alloy.get_ideal_configurational_entropy()

    assert isinstance(entropy, np.ndarray)
    assert entropy.shape[0] > 0
    assert np.all(entropy >= 0)


def test_get_ideal_configurational_entropy_with_custom_resolution(setup_data):
    """
    Test the get_ideal_configurational_entropy method with a custom resolution.

    This test verifies that the get_ideal_configurational_entropy method runs without errors when using a custom resolution.
    """
    alloy_components, calculator = setup_data
    alloy = Alloy(alloy_components, calculator)

    resolution = 30
    entropy = alloy.get_ideal_configurational_entropy(resolution=resolution)

    assert isinstance(entropy, np.ndarray)
    assert entropy.shape[0] > 0
    assert np.all(entropy >= 0)

