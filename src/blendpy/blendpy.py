#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# file: blendpy.py

import numpy as np
from ase.io import read
from ase import Atoms

class Blendpy(Atoms):
  def __init__(self, method='dsi', alloy_basis=[], calculator = None, supercell=[1,1,1], suballoy = False):
    super().__init__()
    self.method = method
    self.alloy_basis: list = alloy_basis
    self.calculator = calculator
    self.supercell = supercell
    self.suballoy = suballoy


if __name__ == '__main__':
  alloy = Blendpy(calculator = calc)
