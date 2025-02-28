<p align="center">
<img src="https://raw.githubusercontent.com/leseixas/blendpy/refs/heads/main/logo.png" style="height: 150px"></p>

[![License: MIT](https://img.shields.io/github/license/leseixas/blendpy?color=green&style=for-the-badge)](LICENSE)    [![PyPI](https://img.shields.io/pypi/v/blendpy?color=red&label=version&style=for-the-badge)](https://pypi.org/project/blendpy/)

# blendpy
**Blendpy** uses atomistic simulations with ASE calculators to compute alloy properties like enthalpy of mixing. It supports binary and multicomponent systems, including alloys and pseudoalloys.

## Installation

Install blendpy easily using pip, Python’s package manager:
```bash
$ pip install blendpy
```

## Getting started

```python
from blendpy import Blendpy

# Calculator
from mace.calculators import mace_mp
calc_mace = mace_mp(model="small",
                    dispersion=False,
                    default_dtype="float32",
                    device='cpu')

# The alloy is created by combining two key components.                
alloy_files = ['Au.cif', 'Pt.cif']

# Supercell to create the dilution.
supercell = [2,2,2]

# Create a Blendpy object.
blendpy = Blendpy(alloy_files, supercell, calculator=calc_mace)

# Optimize the structures.
blendpy.optimize(method=BFGSLineSearch, fmax=0.01, steps=500)

# Calculate the enthalpy of mixing for the AuPt alloy.
enthalpy_of_mixing = blendpy.get_enthalpy_of_mixing(A=0, B=1, npoints=21)
print(enthalpy_of_mixing)
```

## License

This is an open source code under [MIT License](LICENSE.txt).