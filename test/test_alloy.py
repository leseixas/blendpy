import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
import pytest

from blendpy.alloy import Alloy

def test_alloy_initialization():
    alloy_components = ["../data/bulk/Au.cif", "../data/bulk/Pt.cif"]
    sublattice_alloy = "sublattice"
    alloy = Alloy(alloy_components, sublattice_alloy)

    assert alloy.alloy_components == alloy_components
    assert alloy.sublattice_alloy == sublattice_alloy
    assert isinstance(alloy._chemical_elements, list)
    assert alloy._chemical_elements != []
