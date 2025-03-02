<h1 align="center" style="margin-top:20px; margin-bottom:50px;">
<img src="https://raw.githubusercontent.com/leseixas/blendpy/refs/heads/main/logo.png" style="height: 70px"></h1>

[![License: MIT](https://img.shields.io/github/license/leseixas/blendpy?color=green&style=for-the-badge)](LICENSE)    [![PyPI](https://img.shields.io/pypi/v/blendpy?color=red&label=version&style=for-the-badge)](https://pypi.org/project/blendpy/)
![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/leseixas/blendpy/total?style=for-the-badge&logo=github&label=github&color=blue)


**Blendpy** uses atomistic simulations with ASE calculators to compute alloy properties like enthalpy of mixing. It supports binary and multicomponent systems, including alloys and pseudoalloys.

# Installation

Install blendpy easily using pip, Python’s package manager:
```bash
$ pip install blendpy
```

# Getting started

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
optimizer_gold.run(fmax=0.01)

optimizer_platinum = BFGSLineSearch(UnitCellFilter(platinum))
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
                   calculator = calc_mace)
```

Optimize the structures within the `DSIModel` object:
```python
# Optimize the structures
blendpy.optimize(method=BFGSLineSearch, fmax=0.01, logfile=None)
```

Calculate the enthalpy of mixing for the AuPt alloy:
```python
# Calculate the enthalpy of mixing
enthalpy_of_mixing = blendpy.get_enthalpy_of_mixing(npoints=101)
```

Plotting the enthalpy of mixing
```python
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1,1, figsize=(5,5))

x = np.linspace(0, 1, 101)
ax.set_xlabel("$x$", fontsize=20)
ax.set_ylabel("$\Delta H_{mix}$ (kJ/mol)", fontsize=20)
ax.set_xlim(0,1)
ax.set_ylim(-7,7)
ax.set_xticks(np.linspace(0,1,6))
ax.set_yticks(np.arange(-6,7,2))

# Plot the data
ax.plot(x, enthalpy_of_mixing, color='#d53e4f', linewidth=3, zorder=2)
ax.scatter(x[::10], enthalpy_of_mixing[::10], color='#d53e4f', s=80, zorder=2, label="Au$_{1-x}$Pt$_{x}$")

# Reference: H. Okamoto and T.B. Massalski, Bull. Alloy Phase Diagrams 1 (1985) 46.
df_exp = pd.read_csv("data/experimental/exp_AuPt.csv")
ax.plot(df_exp['x'][::2], df_exp['enthalpy'][::2], 's', color='#000000', markersize=8, label="Exp. Data", zorder=1)
ax.legend(loc="best", fontsize=16)

ax.tick_params(direction='in', axis='both', which='major', labelsize=20, width=3, length=8)
ax.set_box_aspect(1)
for spine in ax.spines.values():
    spine.set_linewidth(3)

plt.tight_layout()
# plt.savefig("enthalpy_of_mixing.png", dpi=600, format='png', bbox_inches='tight') # uncomment this if you want to save the figure
plt.show()
```

<p align="center">
<img src="https://raw.githubusercontent.com/leseixas/blendpy/refs/heads/main/figs/enthalpy_of_mixing.png" style="height: 400px"></p>

<p align="center"><a name="fig1">Figure 1</a> - Enthalpy of mixing of the Au-Pt alloy computed using the DSI model and MACE interatomic potentials.</p>


# Phase diagram

By analyzing the mixing enthalpies and entropies, we can calculate the Gibbs free energy of the Au–Pt alloy mixture and determine both the spinodal and binodal (solvus) decomposition curves. These curves, which form key features of the alloy's phase diagram, delineate regions of differing stability: below the binodal curve, the solid solution (Au, Pt) is metastable, whereas it becomes unstable beneath the spinodal curve.

We begin by defining a temperature range over which to calculate the spinodal and binodal curves. Optionally, the results can be saved in CSV files.
```python

temperatures = np.arange(300, 3001, 5)

# spinodal curve
df_spinodal = blendpy.get_spinodal_decomposition(temperatures = temperatures, npoints = 501)
df_spinodal.to_csv("data/phase_diagram/spinodal_AuPt.csv", index=False, header=True, sep=',')

# binodal curve
df_binodal = blendpy.get_binodal_curve(temperatures = temperatures, npoints=501)
df_binodal.to_csv("data/phase_diagram/binodal_AuPt.csv", index=False, header=True, sep=',')
```

To plot the phase diagram featuring the spinodal and binodal decomposition curves, we proceed as follows:
```python
import pandas as pd

# Create figure and axis
fig, ax = plt.subplots(1,1, figsize=(8,8))

x = np.linspace(0, 1, 101)

# Configure axis labels and limits
ax.set_xlabel("$x$", fontsize=20)
ax.set_ylabel("$T$ (K)", fontsize=20)
ax.set_xlim(0,1)
ax.set_ylim(300, 2500)
ax.set_xticks(np.linspace(0,1,6))

# Plot the data
ax.plot(df_spinodal['x'], df_spinodal['t'], color='#d53e4f', linestyle='--', linewidth=3, label="Spinodal curve")
ax.plot(df_binodal['x'], df_binodal['t'], color='#d53e4f', linewidth=3, label="Binodal curve")

# Fill below the curves with transparency (alpha=0.3 means 30% opacity)
ax.fill_between(df_spinodal['x'], df_spinodal['t'], 300, color='#d53e4f', alpha=0.3)
ax.fill_between(df_binodal['x'], df_binodal['t'], 300, color='#d53e4f', alpha=0.3)
ax.legend(loc="best", fontsize=20)

# Add text annotations
ax.text(0.2, 1500, "Stable", fontsize=20, ha='center', va='center')
ax.text(0.4, 950, "Metastable", fontsize=20, ha='center', va='center', rotation=60)
ax.text(0.7, 700, "Unstable", fontsize=20, ha='center', va='center')

# Customize tick parameters
ax.tick_params(direction='in', axis='both', which='major', labelsize=20, width=3, length=8)
ax.set_box_aspect(1)
for spine in ax.spines.values():
    spine.set_linewidth(3)

plt.tight_layout()
# plt.savefig("phase_diagram.png", dpi=600, format='png', bbox_inches='tight') # uncomment this if you want to save the figure
plt.show()
```

<p align="center">
<img src="https://raw.githubusercontent.com/leseixas/blendpy/refs/heads/main/figs/phase_diagram.png" style="height: 500px"></p>

<p align="center"><a name="fig1">Figure 2</a> - Phase diagram of the Au–Pt alloy computed using the DSI model and MACE interatomic potentials.</p>


# Enthalpy of mixing with polymorphism

```python

# TODO: 

```

# License

This is an open source code under [MIT License](LICENSE).

# Acknowledgements

We thank financial support from FAPESP [(Grant No. 2022/14549-3)](https://bvs.fapesp.br/pt/auxilios/111791/materiais-de-alta-entropia-inteligiveis-desenvolvendo-modelos-dados-e-aplicacoes/), INCT Materials Informatics (Grant No. 406447/2022-5), and CNPq (Grant No. 311324/2020-7).