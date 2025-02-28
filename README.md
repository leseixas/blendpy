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
from blendpy import DSIModel

# Create a calculator object to optimize structures.
from mace.calculators import mace_mp
calc_mace = mace_mp(model="small",
                    dispersion=False,
                    default_dtype="float32",
                    device='cpu')

# Create a DSIModel object.
blendpy = DSIModel(alloy_components = ['Au.cif', 'Pt.cif'],
                   supercell = [2,2,2],
                   calculator=calc_mace)

# Optimize the structures.
blendpy.optimize(method=BFGSLineSearch, fmax=0.01)

# Calculate the enthalpy of mixing for the AuPt alloy.
enthalpy_of_mixing = blendpy.get_enthalpy_of_mixing(npoints=101)
```

### Plotting the enthalpy of mixing

```python
import matplotlib.pyplot as plt

TODO

```

## Spinodal decomposition curve

```python

TODO

```

## Phase diagram

```python

TODO

```


## License

This is an open source code under [MIT License](LICENSE).

## Acknowledgements

We thank financial support from FAPESP [(Grant No. 2022/14549-3)](https://bvs.fapesp.br/pt/auxilios/111791/materiais-de-alta-entropia-inteligiveis-desenvolvendo-modelos-dados-e-aplicacoes/), INCT Materials Informatics (Grant No. 406447/2022-5), and CNPq (Grant No. 311324/2020-7).