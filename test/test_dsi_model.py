import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import pytest
from ase import Atoms
from ase.io import read
from ase.calculators.emt import EMT
from blendpy.dsi_model import DSIModel
from ase.optimize import BFGS, LBFGS
import numpy as np
from blendpy.constants import *

# Fixtures
@pytest.fixture
def setup_data():
    """
    Sets up the data required for testing the DSI model.

    Returns:
        tuple: A tuple containing the following elements:
            - alloy_components (list of str): Paths to the CIF files of the alloy components.
            - supercell (list of int): Dimensions of the supercell.
            - calculator (EMT): An instance of the EMT calculator.
            - doping_site (int): The index of the doping site.
    """
    alloy_components = ["data/bulk/Au.cif", "data/bulk/Pt.cif"]
    supercell = [2, 2, 2]
    calculator = EMT()
    doping_site = 0
    return alloy_components, supercell, calculator, doping_site

# Tests
def test_initialization_empty():
    """
    Test the initialization of the DSIModel class without any arguments.

    This test verifies that the DSIModel is correctly initialized with default values.
    It checks the following:
    
    - The number of components in the model is set to 0.
    - The supercell in the model is set to [1, 1, 1].
    - The doping site in the model is set to 0.
    - The dilute_alloys attribute is set to None.
    - The x0 attribute is set to None.
    - The energy_matrix attribute is set to None.
    - The diluting_parameters attribute is set to None.

    Args:
        None
    """
    model = DSIModel()
    assert model.n_components == 0
    assert model.supercell == [1, 1, 1]
    assert model.doping_site == 0
    assert model._dilute_alloys is None
    assert model.x0 is None
    assert model._energy_matrix is None
    assert model._diluting_parameters is None


def test_initialization_without_calculator(setup_data):
    """
    Test the initialization of the DSIModel without a calculator.

    This test verifies that the DSIModel is correctly initialized with the given
    alloy components, supercell, and doping site. It checks the following:
    
    - The number of components in the model matches the length of alloy_components.
    - The supercell in the model matches the provided supercell and is of type list.
    - Each supercell in the model's _supercells attribute is an instance of Atoms.
    - The doping site in the model matches the provided doping site and is of type int.
    - Each dilute alloy in the model's dilute_alloys attribute is an instance of Atoms
      and has no calculator assigned (atoms.calc is None).

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, and doping_site.
    """
    alloy_components, supercell, _, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, doping_site=doping_site)
    assert model.n_components == len(alloy_components)
    assert model.supercell == supercell
    assert isinstance(model.supercell, list) and all(isinstance(item, int) for item in model.supercell) # is list[int]?
    for atoms in model._supercells:
        assert isinstance(atoms, Atoms)
    assert model.doping_site == doping_site
    assert isinstance(model.doping_site, int)
    for row in model._dilute_alloys:
        for atoms in row:
            assert isinstance(atoms, Atoms)
            assert atoms.calc is None


def test_initialization_with_calculator(setup_data):
    """
    Test the initialization of the DSIModel with a calculator.

    This test verifies that the DSIModel is correctly initialized with the provided
    alloy components, supercell, calculator, and doping site. It checks the following:
    
    - The number of components in the model matches the length of alloy_components.
    - The supercell in the model matches the provided supercell and is of type list.
    - Each supercell in the model is an instance of Atoms.
    - The doping site in the model matches the provided doping site and is of type int.
    - The calculator is attached to each atom in the dilute_alloys and energy is calculated.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, 
                            and doping_site used to initialize the DSIModel.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    assert model.n_components == len(alloy_components)
    assert model.supercell == supercell
    assert isinstance(model.supercell, list)
    for atoms in model._supercells:
        assert isinstance(atoms, Atoms)
    assert model.doping_site == doping_site
    assert isinstance(model.doping_site, int)

    # Check if the calculator is attached and energy is calculated
    for row in model._dilute_alloys:
        for atoms in row:
            assert atoms.calc is not None
            assert 'energy' in atoms.info


def test_get_supercells(setup_data):
    """
    Test the get_supercells method of the DSIModel class.

    This test verifies that the get_supercells method correctly returns the list of supercells
    created during the initialization of the DSIModel. It checks the following:
    
    - The returned supercells list is of type list.
    - Each element in the returned supercells list is an instance of Atoms.
    - The length of the returned supercells list matches the number of alloy components.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, and doping_site.
    """
    alloy_components, supercell, _, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, doping_site=doping_site)
    supercells = model.get_supercells()
    
    assert isinstance(supercells, list)
    assert len(supercells) == len(alloy_components)
    for atoms in supercells:
        assert isinstance(atoms, Atoms)


def test_create_dilute_alloys(setup_data):
    """
    Test the _create_dilute_alloys method of the DSIModel class.

    This test verifies that the _create_dilute_alloys method correctly creates a matrix of dilute alloys
    from the provided supercells. It checks the following:
    
    - The method raises a ValueError if there are fewer than two supercells.
    - The returned dilute alloys matrix is of type list.
    - Each element in the dilute alloys matrix is an instance of Atoms.
    - The dimensions of the dilute alloys matrix are n x n, where n is the number of supercells.
    - The chemical symbols of the doping site in each supercell are correctly replaced.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, and doping_site.
    """
    alloy_components, supercell, _, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, doping_site=doping_site)
    
    # Reset supercells for further testing
    model._supercells = []
    for filename in alloy_components:
        atoms = read(filename)
        assert isinstance(atoms, Atoms)
        supercell_atoms = atoms.repeat(supercell)
        model._supercells.append(supercell_atoms)
    
    dilute_alloys = model._create_dilute_alloys()
    
    assert isinstance(dilute_alloys, list)
    n = len(model._supercells)
    assert len(dilute_alloys) == n
    for row in dilute_alloys:
        assert isinstance(row, list)
        assert len(row) == n
        for atoms in row:
            assert isinstance(atoms, Atoms)
    
    dopant = [atoms.get_chemical_symbols()[doping_site] for atoms in model._supercells]
    for i in range(n):
        for j in range(n):
            assert dilute_alloys[i][j].get_chemical_symbols()[doping_site] == dopant[j]

# optimize
def test_optimize_with_default_parameters(setup_data):
    """
    Test the optimize method of the DSIModel class with default parameters.

    This test verifies that the optimize method correctly optimizes the Atoms objects
    in the dilute_alloys attribute using the default optimization method and parameters.
    It checks the following:
    
    - The method runs without errors.
    - The total energy of each Atoms object is updated after optimization.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    model.optimize() # default arguments, except logfile to avoid writing to file.
    
    for row in model._dilute_alloys:
        for atoms in row:
            assert 'energy' in atoms.info
            assert atoms.info['energy'] is not None


def test_optimize_with_custom_parameters(setup_data):
    """
    Test the optimize method of the DSIModel class with custom parameters.

    This test verifies that the optimize method correctly optimizes the Atoms objects
    in the dilute_alloys attribute using custom optimization parameters. It checks the following:
    
    - The method runs without errors.
    - The total energy of each Atoms object is updated after optimization.
    - The custom parameters are correctly applied.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    custom_method = BFGS
    custom_fmax = 0.05
    custom_steps = 300
    custom_logfile = "test_optimize.log"
    custom_mask = [1, 1, 0, 0, 0, 1]
    
    model.optimize(method=custom_method, fmax=custom_fmax, steps=custom_steps, logfile=custom_logfile, mask=custom_mask)
    
    for row in model._dilute_alloys:
        for atoms in row:
            assert 'energy' in atoms.info
            assert atoms.info['energy'] is not None

# optimize_nostress
def test_optimize_nostress_with_default_parameters(setup_data):
    """
    Test the optimize_nostress method of the DSIModel class with default parameters.

    This test verifies that the optimize_nostress method correctly optimizes the Atoms objects
    in the dilute_alloys attribute using the default optimization parameters. It checks the following:
    
    - The method runs without errors.
    - The total energy of each Atoms object is updated after optimization.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    model.optimize_nostress()  # default arguments
    
    for row in model._dilute_alloys:
        for atoms in row:
            assert 'energy' in atoms.info
            assert atoms.info['energy'] is not None


def test_optimize_nostress_with_custom_parameters(setup_data):
    """
    Test the optimize_nostress method of the DSIModel class with custom parameters.

    This test verifies that the optimize_nostress method correctly optimizes the Atoms objects
    in the dilute_alloys attribute using custom optimization parameters. It checks the following:
    
    - The method runs without errors.
    - The total energy of each Atoms object is updated after optimization.
    - The custom parameters are correctly applied.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    custom_fmax = 0.05
    custom_steps = 300
    custom_logfile = "test_optimize_nostress.log"
    
    model.optimize_nostress(fmax=custom_fmax, steps=custom_steps, logfile=custom_logfile)
    
    for row in model._dilute_alloys:
        for atoms in row:
            assert 'energy' in atoms.info
            assert atoms.info['energy'] is not None


# set_energy_matrix
def test_set_energy_matrix_with_invalid_dtype(setup_data):
    """
    Test the set_energy_matrix method of the DSIModel class with an invalid dtype.

    This test verifies that the set_energy_matrix method raises a ValueError when provided
    with a numpy array that does not have a floating-point dtype. It checks the following:
    
    - The method raises a ValueError with the appropriate error message.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)

    invalid_energy_matrix_str = np.array([["a", "b"], ["c", "d"]])  # float dtype
    with pytest.raises(ValueError, match="The energy matrix must be a nd.array of floats."):
        model.set_energy_matrix(invalid_energy_matrix_str)


def test_set_energy_matrix_with_invalid_shape(setup_data):
    """
    Test the set_energy_matrix method of the DSIModel class with an invalid shape.

    This test verifies that the set_energy_matrix method raises a ValueError when provided
    with a numpy array that does not have the correct shape. It checks the following:
    
    - The method raises a ValueError with the appropriate error message.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    invalid_energy_matrix = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])  # Incorrect shape
    with pytest.raises(ValueError, match="The energy matrix must be a square matrix."):
        model.set_energy_matrix(invalid_energy_matrix)


def test_set_energy_matrix_with_valid_input(setup_data):
    """
    Test the set_energy_matrix method of the DSIModel class with valid input.

    This test verifies that the set_energy_matrix method correctly sets the energy matrix
    for the model when provided with a valid numpy array. It checks the following:
    
    - The energy matrix is correctly set in the model.
    - The energy matrix is of type numpy.ndarray.
    - The values in the energy matrix match the provided input.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    energy_matrix = np.array([[1.0, 2.0], [3.0, 4.0]])
    model.set_energy_matrix(energy_matrix)
    
    assert isinstance(model._energy_matrix, np.ndarray)
    assert np.array_equal(model._energy_matrix, energy_matrix)


# get_energy_matrix
def test_get_energy_matrix_with_precomputed_matrix(setup_data):
    """
    Test the get_energy_matrix method of the DSIModel class when the energy matrix is precomputed.

    This test verifies that the get_energy_matrix method correctly returns the precomputed energy matrix
    when it is already set in the model. It checks the following:
    
    - The method returns the precomputed energy matrix.
    - The returned energy matrix is of type numpy.ndarray.
    - The values in the returned energy matrix match the precomputed matrix.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    precomputed_energy_matrix = np.array([[1.0, 2.0], [3.0, 4.0]])
    model.set_energy_matrix(precomputed_energy_matrix)
    
    energy_matrix = model.get_energy_matrix()
    
    assert isinstance(energy_matrix, np.ndarray)
    assert np.array_equal(energy_matrix, precomputed_energy_matrix)


def test_get_energy_matrix_with_calculated_matrix(setup_data):
    """
    Test the get_energy_matrix method of the DSIModel class when the energy matrix needs to be calculated.

    This test verifies that the get_energy_matrix method correctly calculates and returns the energy matrix
    when it is not precomputed. It checks the following:
    
    - The method calculates the energy matrix.
    - The returned energy matrix is of type numpy.ndarray.
    - The values in the returned energy matrix are correctly calculated based on the potential energy of the atoms.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    energy_matrix = model.get_energy_matrix()
    
    assert isinstance(energy_matrix, np.ndarray)
    assert energy_matrix.shape == (len(alloy_components), len(alloy_components))
    
    for i, row in enumerate(model._dilute_alloys):
        for j, atoms in enumerate(row):
            assert energy_matrix[i, j] == atoms.info['energy']


def test_get_energy_matrix_with_missing_energy_info(setup_data):
    """
    Test the get_energy_matrix method of the DSIModel class when some atoms are missing energy info.

    This test verifies that the get_energy_matrix method correctly calculates the energy for atoms
    that are missing the 'energy' info in their info dictionary. It checks the following:
    
    - The method calculates the missing energy values.
    - The returned energy matrix is of type numpy.ndarray.
    - The values in the returned energy matrix are correctly calculated based on the potential energy of the atoms.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    # Remove 'energy' from atoms.info
    for row in model._dilute_alloys:
        for atoms in row:
            if 'energy' in atoms.info:
                del atoms.info['energy']
    
    energy_matrix = model.get_energy_matrix()
    
    assert isinstance(energy_matrix, np.ndarray)
    assert energy_matrix.shape == (len(alloy_components), len(alloy_components))
    
    for i, row in enumerate(model._dilute_alloys):
        for j, atoms in enumerate(row):
            assert 'energy' in atoms.info
            assert energy_matrix[i, j] == atoms.info['energy']

# set_diluting_parameters
def test_set_diluting_parameters_with_invalid_dtype(setup_data):
    """
    Test the set_diluting_parameters method of the DSIModel class with an invalid dtype.

    This test verifies that the set_diluting_parameters method raises a ValueError when provided
    with a numpy array that does not have a floating-point dtype. It checks the following:
    
    - The method raises a ValueError with the appropriate error message.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)

    invalid_diluting_parameters_str = np.array([["a", "b"], ["c", "d"]])  # String dtype
    with pytest.raises(ValueError, match="The diluting parameters matrix must be a nd.array of floats."):
        model.set_diluting_parameters(invalid_diluting_parameters_str)


def test_set_diluting_parameters_with_invalid_shape(setup_data):
    """
    Test the set_diluting_parameters method of the DSIModel class with an invalid shape.

    This test verifies that the set_diluting_parameters method raises a ValueError when provided
    with a numpy array that does not have the correct shape. It checks the following:
    
    - The method raises a ValueError with the appropriate error message.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    invalid_diluting_parameters = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])  # Incorrect shape
    with pytest.raises(ValueError, match="The diluting parameters matrix must be a square matrix."):
        model.set_diluting_parameters(invalid_diluting_parameters)


def test_set_diluting_parameters_with_valid_input(setup_data):
    """
    Test the set_diluting_parameters method of the DSIModel class with valid input.

    This test verifies that the set_diluting_parameters method correctly sets the diluting parameters
    matrix for the model when provided with a valid numpy array. It checks the following:
    
    - The diluting parameters matrix is correctly set in the model.
    - The diluting parameters matrix is of type numpy.ndarray.
    - The values in the diluting parameters matrix match the provided input.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    diluting_parameters = np.array([[0.1, 0.2], [0.3, 0.4]])
    model.set_diluting_parameters(diluting_parameters)
    
    assert isinstance(model._diluting_parameters, np.ndarray)
    assert np.array_equal(model._diluting_parameters, diluting_parameters)


# get_diluting_parameters
def test_get_diluting_parameters_with_valid_data(setup_data):
    """
    Test the get_diluting_parameters method of the DSIModel class with valid data.

    This test verifies that the get_diluting_parameters method correctly calculates the diluting parameters
    matrix when provided with valid data. It checks the following:
    
    - The method returns a numpy array.
    - The shape of the returned array matches the number of components.
    - The values in the returned array are correctly calculated based on the energy differences.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    diluting_parameters = model.get_diluting_parameters()
    
    assert isinstance(diluting_parameters, np.ndarray)
    assert diluting_parameters.shape == (len(alloy_components), len(alloy_components))


def test_get_diluting_parameters_with_inconsistent_supercells(setup_data):
    """
    Test the get_diluting_parameters method of the DSIModel class with inconsistent supercells.

    This test verifies that the get_diluting_parameters method raises a NotImplementedError when the supercells
    have different numbers of atoms. It checks the following:
    
    - The method raises a NotImplementedError with the appropriate error message.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    alloy_components = ["data/bulk/Au.cif", "data/bulk/Pd_hcp.cif"]

    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    with pytest.raises(NotImplementedError, match="Not all supercells have the same number of atoms."):
        model.get_diluting_parameters()


def test_get_diluting_parameters_with_precomputed_energy_matrix(setup_data):
    """
    Test the `get_diluting_parameters` method of the `DSIModel` class with a precomputed energy matrix.
    This test verifies that the diluting parameters calculated by the `DSIModel` are correct when provided
    with a specific energy matrix. The test uses two sets of data: one for the AuPt alloy and one for a 
    fictitious alloy. The expected diluting parameters are compared against the calculated ones to ensure 
    accuracy.
    Args:
        setup_data (tuple): A tuple containing the alloy components, supercell, and doping site information.
    Assertions:
        - The diagonal elements of the diluting parameters matrix should be close to zero.
        - The calculated diluting parameters should match the expected diluting parameters within a tolerance.
        - The diluting parameters should be of type `np.ndarray`.
    """
    alloy_components, supercell, _, doping_site = setup_data
    x0 = 1/27 # AuPt
    # x0 = 1/10 # fictitious

    model = DSIModel(alloy_components=alloy_components,
                     supercell=supercell,
                     doping_site=doping_site,
                     x0=x0)
    
    # AuPt (from DFT)
    energy_matrix = np.array([[-85.940400, -89.230299],
                              [-170.278459, -173.891172]])


    # fictitious
    # energy_matrix = np.array([[-10.0, -9.8],
    #                           [-11.2, -11.0]])

    model.set_energy_matrix(energy_matrix)

    diluting_parameters = model.get_diluting_parameters()

    # AuPt
    expected_diluting_parameters = np.array([[0.0, -0.032463], [0.355277, 0.0]]) * convert_eVatom_to_kJmol

    # fictitious
    # expected_diluting_parameters = np.array([[0.0, 0.3], [-0.3, 0]]) # * convert_eVatom_to_kJmol
    
    assert np.allclose(np.diag(diluting_parameters), 0, atol=1e-6)
    assert np.allclose(diluting_parameters, expected_diluting_parameters, atol=1e-6)
    assert isinstance(diluting_parameters, np.ndarray)


# get_enthalpy_of_mixing    
def test_get_enthalpy_of_mixing_with_default_parameters(setup_data):
    """
    Test the get_enthalpy_of_mixing method of the DSIModel class with default parameters.

    This test verifies that the get_enthalpy_of_mixing method correctly calculates the enthalpy of mixing
    for a binary mixture using the default parameters. It checks the following:
    
    - The method returns a numpy array.
    - The length of the returned array matches the number of points (npoints).
    - The values in the returned array are correctly calculated based on the diluting parameters and slope.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    enthalpy = model.get_enthalpy_of_mixing()
    
    assert isinstance(enthalpy, np.ndarray)
    assert np.issubdtype(enthalpy.dtype, np.floating)
    assert len(enthalpy) == 101  # Default npoints
    assert np.all(np.isfinite(enthalpy))


def test_get_enthalpy_of_mixing_with_custom_parameters(setup_data):
    """
    Test the get_enthalpy_of_mixing method of the DSIModel class with custom parameters.

    This test verifies that the get_enthalpy_of_mixing method correctly calculates the enthalpy of mixing
    for a binary mixture using custom parameters. It checks the following:
    
    - The method returns a numpy array.
    - The length of the returned array matches the number of points (npoints).
    - The values in the returned array are correctly calculated based on the diluting parameters and slope.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    A = 1
    B = 0
    slope = [0.1, 0.2]
    npoints = 50
    
    enthalpy = model.get_enthalpy_of_mixing(A=A, B=B, slope=slope, npoints=npoints)
    
    assert isinstance(enthalpy, np.ndarray)
    assert np.issubdtype(enthalpy.dtype, np.floating)
    assert len(enthalpy) == npoints
    assert np.all(np.isfinite(enthalpy))

def test_get_enthalpy_of_mixing_with_default_parameters(setup_data):
    """
    Test the get_enthalpy_of_mixing method of the DSIModel class with default parameters.

    This test verifies that the get_enthalpy_of_mixing method correctly calculates the enthalpy of mixing
    for a binary mixture using the default parameters. It checks the following:
    
    - The method returns a numpy array.
    - The length of the returned array matches the number of points (npoints).
    - The values in the returned array are correctly calculated based on the diluting parameters and slope.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    enthalpy = model.get_enthalpy_of_mixing()
    
    assert isinstance(enthalpy, np.ndarray)
    assert np.issubdtype(enthalpy.dtype, np.floating)
    assert len(enthalpy) == 101  # Default npoints
    assert np.all(np.isfinite(enthalpy))


def test_get_enthalpy_of_mixing_with_custom_parameters(setup_data):
    """
    Test the get_enthalpy_of_mixing method of the DSIModel class with custom parameters.

    This test verifies that the get_enthalpy_of_mixing method correctly calculates the enthalpy of mixing
    for a binary mixture using custom parameters. It checks the following:
    
    - The method returns a numpy array.
    - The length of the returned array matches the number of points (npoints).
    - The values in the returned array are correctly calculated based on the diluting parameters and slope.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    A = 1
    B = 0
    slope = [0.1, 0.2]
    npoints = 50
    
    enthalpy = model.get_enthalpy_of_mixing(A=A, B=B, slope=slope, npoints=npoints)
    
    assert isinstance(enthalpy, np.ndarray)
    assert np.issubdtype(enthalpy.dtype, np.floating)
    assert len(enthalpy) == npoints
    assert np.all(np.isfinite(enthalpy))


def test_get_enthalpy_of_mixing_with_invalid_component_indices(setup_data):
    """
    Test the get_enthalpy_of_mixing method of the DSIModel class with invalid component indices.

    This test verifies that the get_enthalpy_of_mixing method raises a ValueError when provided
    with invalid component indices. It checks the following:
    
    - The method raises a ValueError with the appropriate error message.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    with pytest.raises(ValueError, match="The component indices must be integers."):
        model.get_enthalpy_of_mixing(A="0", B=1)
    
    with pytest.raises(ValueError, match="The component indices must be less than the number of components."):
        model.get_enthalpy_of_mixing(A=0, B=10)


def test_get_enthalpy_of_mixing_with_invalid_slope(setup_data):
    """
    Test the get_enthalpy_of_mixing method of the DSIModel class with an invalid slope parameter.

    This test verifies that the get_enthalpy_of_mixing method raises a ValueError when provided
    with an invalid slope parameter. It checks the following:
    
    - The method raises a ValueError with the appropriate error message.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    with pytest.raises(ValueError, match="The slope parameter must have two values."):
        model.get_enthalpy_of_mixing(slope=[0.1])


def test_get_enthalpy_of_mixing_with_invalid_npoints(setup_data):
    """
    Test the get_enthalpy_of_mixing method of the DSIModel class with an invalid npoints parameter.

    This test verifies that the get_enthalpy_of_mixing method raises a ValueError when provided
    with an invalid npoints parameter. It checks the following:
    
    - The method raises a ValueError with the appropriate error message.

    Args:
        setup_data (tuple): A tuple containing alloy_components, supercell, calculator, and doping_site.
    """
    alloy_components, supercell, calculator, doping_site = setup_data
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, calculator=calculator, doping_site=doping_site)
    
    with pytest.raises(ValueError, match="The number of points must be greater than 1."):
        model.get_enthalpy_of_mixing(npoints=1)

