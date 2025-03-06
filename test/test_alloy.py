import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from blendpy.alloy import Alloy

def test_alloy_initialization():
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
    alloy_components = ["data/bulk/Au.cif", "data/bulk/Pt.cif"]
    alloy = Alloy(alloy_components)

    assert alloy.alloy_components == alloy_components
    assert isinstance(alloy._chemical_elements, list)
    assert alloy._chemical_elements == [['Au','Au','Au','Au'],['Pt','Pt','Pt','Pt']]

def test_get_chemical_elements():
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
    alloy_components = ["data/bulk/Au.cif", "data/bulk/Pt.cif"]
    alloy = Alloy(alloy_components)
    
    expected_elements = [['Au', 'Au', 'Au', 'Au'], ['Pt', 'Pt', 'Pt', 'Pt']]
    assert alloy.get_chemical_elements() == expected_elements



