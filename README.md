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

First, import the necessary modules from ASE and MACE:
```python
from ase.io import write
from ase.build import bulk
from ase.optimize import BFGSLineSearch
from ase.filters import UnitCellFilter
```

Next, create `Atoms` objects for gold (Au) and platinum (Pt) using the `bulk` function:
```python
# Create Au and Pt Atoms object
gold = bulk("Au", cubic=True)
platinum = bulk("Pt", cubic=True)
```

Create a MACE calculator object to optimize the structures:
```python
# Initialize the MACE calculator
from mace.calculators import mace_mp
calc_mace = mace_mp(model="small",
                    dispersion=False,
                    default_dtype="float32",
                    device='cpu')
```

Assign the calculator to the `Atoms` objects:
```python
# Assign the calculator to the Atoms objects
gold.calc = calc_mace
platinum.calc = calc_mace
```

Optimize the unit cells of Au and Pt using the `BFGSLineSearch` optimizer:
```python
# Optimize Au and Pt unit cells
optimizer_gold = BFGSLineSearch(UnitCellFilter(gold))
optimizer_platinum = BFGSLineSearch(UnitCellFilter(platinum))
optimizer_gold.run(fmax=0.01)
optimizer_platinum.run(fmax=0.01)
```

Save the optimized unit cells to CIF files:
```python
# Save the optimized unit cells for Au and Pt
write("Au.cif", gold)
write("Pt.cif", platinum)

```

Now, import the `DSIModel` from blendpy and create a `DSIModel` object using the optimized structures:
```python
from blendpy import DSIModel

# Create a DSIModel object
blendpy = DSIModel(alloy_components = ['Au.cif', 'Pt.cif'],
                   supercell = [2,2,2],
                   calculator=calc_mace)
```

Optimize the structures within the `DSIModel` object:
```python
# Optimize the structures
blendpy.optimize(method=BFGSLineSearch, fmax=0.01)
```

Calculate the enthalpy of mixing for the AuPt alloy:
```python
# Calculate the enthalpy of mixing
enthalpy_of_mixing = blendpy.get_enthalpy_of_mixing(npoints=101)
```

Plottingthe enthalpy of mixing
```python
import matplotlib.pyplot as plt

# TODO: Add code to plot the enthalpy of mixing

```

## Spinodal decomposition curve

```python

# TODO: Add code to calculate and plot the spinodal decomposition curve

```

## Phase diagram

```python

# TODO: Add code to calculate and plot the phase diagram

```


## License

This is an open source code under [MIT License](LICENSE).

## Acknowledgements

We thank financial support from FAPESP [(Grant No. 2022/14549-3)](https://bvs.fapesp.br/pt/auxilios/111791/materiais-de-alta-entropia-inteligiveis-desenvolvendo-modelos-dados-e-aplicacoes/), INCT Materials Informatics (Grant No. 406447/2022-5), and CNPq (Grant No. 311324/2020-7).