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


if __name__ == '__main__':
  from gpaw import GPAW, PW, FermiDirac, Davidson, Mixer
  calc = GPAW(mode=PW(400),
              xc='PBE',
              kpts=(7,7,7),
              occupations=FermiDirac(5),
              eigensolver=Davidson(5),
              spinpol=False,
              mixer=Mixer(0.05, 5, 100))

  alloy = Blendpy(alloy_basis=['Au_fcc.vasp', 'Pt_fcc.vasp'],
                  calculator = calc,
                  supercell=[3,3,3])
