#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# file: blendpy.py

import numpy as np
from ase.io import read
from ase import Atoms
from ase.optimize import BFGS, MDMin, FIRE
from ase.filters import UnitCellFilter

class Blendpy(Atoms):
  def __init__(self, method='dsi', alloy_basis=[], calculator = None, optimizer = None, supercell=[1,1,1], suballoy = False):
    super().__init__()
    self.method = method
    self.alloy_basis: list = alloy_basis
    self.calculator = calculator
    self.optimizer = optimizer
    self.supercell = supercell
    self.suballoy = suballoy
    self.atoms_alloy_basis = None # self.get_alloy_basis()
    self.atoms_supercells = None # self.get_supercells()

  def get_alloy_basis(self):
    atoms_alloy_basis = []
    for alloy in self.alloy_basis: 
      atoms = read(alloy)
      atoms_alloy_basis.append(atoms)
    return atoms_alloy_basis

  def get_supercells(self):
    atoms_supercells = []
    for alloy in self.atoms_alloy_basis:
      alloy_sc = alloy.repeat(self.supercell) 
      atoms_supercells.append(alloy_sc)
    return atoms_supercells


if __name__ == '__main__':
  import warnings
  warnings.filterwarnings("ignore")
#  from gpaw import GPAW, PW, FermiDirac, Davidson, Mixer
#  calc_gpaw = GPAW(mode=PW(400),
#                   xc='PBE',
#                   kpts=(7,7,7),
#                   occupations=FermiDirac(5),
#                   eigensolver=Davidson(5),
#                   spinpol=False,
#                   mixer=Mixer(0.05, 5, 100))

  from mace.calculators import mace_mp
  calc_mace = mace_mp(model="small",
                      dispersion=False,
                      default_dtype="float32",
                      device='cpu')

  alloy = Blendpy(alloy_basis=['../../test/Au.vasp', '../../Pt.vasp'],
                  calculator = calc_mace,
                  supercell=[3,3,3])
