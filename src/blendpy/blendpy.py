#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# file: blendpy.py

# This code is part of blendpy.
# MIT License
#
# Copyright (c) 2025 Leandro Seixas Rocha <leandro.fisica@gmail.com> 
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

'''
Module blendpy
'''

import numpy as np
import pandas as pd
from ase.io import read
from ase import Atoms
from ase.optimize import BFGS, BFGSLineSearch, LBFGS, LBFGSLineSearch, MDMin, GPMin, FIRE, FIRE2, ODE12r, GoodOldQuasiNewton
from ase.filters import UnitCellFilter

# Example usage:
if __name__ == '__main__':
    from alloy import Alloy
    from dsi_model import DSIModel
    from polymorph import Polymorph

    import warnings
    warnings.filterwarnings("ignore")

    ### MACE calculator
    from mace.calculators import mace_mp
    calc_mace = mace_mp(model="small",
                        dispersion=False,
                        default_dtype="float32",
                        device='cpu')

    ### GPAW calculator
    # from gpaw import GPAW, PW, FermiDirac, Davidson, Mixer
    # calc_gpaw = GPAW(mode=PW(500),
    #                  xc='PBE',
    #                  kpts=(7,7,7),
    #                  occupations=FermiDirac(0.1),
    #                  eigensolver=Davidson(5),
    #                  spinpol=False,
    #                  mixer=Mixer(0.05, 5, 100))

    # Example:

    blendpy = DSIModel(alloy_components = ['../../test/Au.vasp', '../../test/Pt.vasp'],
                       supercell = [2,2,2],
                       calculator = calc_mace)

    # Optimize all structures.
    blendpy.optimize(method=BFGSLineSearch, fmax=0.01, steps=500, logfile=None)

    temperatures = np.arange(300, 3001, 5)
    df_spinodal = blendpy.get_spinodal_decomposition(temperatures = temperatures, A=0, B=1, npoints = 501)
    df_spinodal.to_csv("../../test/spinodal_AuPt.csv", index=False, header=True, sep=',')

    df_binodal = blendpy.get_binodal_curve(temperatures = temperatures, A=0, B=1, npoints=501)
    df_binodal.to_csv("../../test/binodal_AuPt.csv", index=False, header=True, sep=',')

    # blendpy = Polymorph(alpha='../../test/Pt_fcc.vasp', beta='../../test/Pt_bcc.vasp', calculator = calc_mace)
    # blendpy.optimize()
    # print("Difference between alpha and beta phases:")
    # print(blendpy.get_structural_energy_transition())
