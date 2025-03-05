import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
import pytest

from blendpy.alloy import Alloy

def test_alloy_initialization():
    alloy_components = ["../data/bulk/Au.cif", "../data/bulk/Pt.cif"]
    alloy = Alloy(alloy_components)

    assert alloy.alloy_components == alloy_components
    assert isinstance(alloy._chemical_elements, list)
    assert alloy._chemical_elements == [['Au','Au','Au','Au'],['Pt','Pt','Pt','Pt']]

def test_get_chemical_elements():
    alloy_components = ["../data/bulk/Au.cif", "../data/bulk/Pt.cif"]
    alloy = Alloy(alloy_components)
    
    expected_elements = [['Au', 'Au', 'Au', 'Au'], ['Pt', 'Pt', 'Pt', 'Pt']]
    assert alloy.get_chemical_elements() == expected_elements



