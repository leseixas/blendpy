import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import pytest
from ase import Atoms
from ase.calculators.emt import EMT
from blendpy.dsi_model import DSIModel

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
    alloy_components = ["../data/bulk/Au.cif", "../data/bulk/Pt.cif"]
    supercell = [2, 2, 2]
    calculator = EMT()
    doping_site = 0
    return alloy_components, supercell, calculator, doping_site

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
    for row in model.dilute_alloys:
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
    for row in model.dilute_alloys:
        for atoms in row:
            assert atoms.calc is not None
            assert 'energy' in atoms.info

def test_initialization_with_diluting_parameters(setup_data):
    """
    Test the initialization of the DSIModel with diluting parameters.

    This test ensures that the DSIModel is correctly initialized with the provided
    diluting parameters and that the diluting parameters are stored as a list. It 
    also verifies that each element in the `dilute_alloys` attribute is an instance 
    of the `Atoms` class and that the `calc` attribute of each `Atoms` instance is 
    None.

    Args:
        setup_data (tuple): A tuple containing the alloy components, supercell, 
                            and doping site used for initializing the DSIModel.

    Asserts:
        - The `diluting_parameters` attribute of the model is equal to the provided 
          diluting parameters.
        - The `diluting_parameters` attribute is a list.
        - Each element in the `dilute_alloys` attribute is an instance of the `Atoms` 
          class.
        - The `calc` attribute of each `Atoms` instance in `dilute_alloys` is None.
    """
    alloy_components, supercell, _, doping_site = setup_data
    diluting_parameters = [[0.1, 0.2], [0.3, 0.4]]
    model = DSIModel(alloy_components=alloy_components, supercell=supercell, diluting_parameters=diluting_parameters, doping_site=doping_site)
    assert model.diluting_parameters == diluting_parameters
    assert isinstance(model.diluting_parameters, list)
    for row in model.dilute_alloys:
        for atoms in row:
            assert isinstance(atoms, Atoms)
            assert atoms.calc is None


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



