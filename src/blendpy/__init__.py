# -*- coding: utf-8 -*-
# file: __init__.py

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

import warnings
warnings.filterwarnings("ignore")

import os
from datetime import datetime
from socket import gethostname
from sys import version as __python_version__
from sys import executable as __python_executable__
from ase import __version__ as __ase_version__
from numpy import __version__ as __numpy_version__
from pandas import __version__ as __pandas_version__

from ._version import __version__
from .constants import R, convert_eVatom_to_kJmol
from .alloy import Alloy
from .dsi_model import DSIModel
from .polymorph import Polymorph
from .phase_diagram import PhaseDiagram
# from .ace import ACE
# from .sqs import SQS
# from .intermetallics import Intermetallics

__all__ = ['Alloy', 'DSIModel', 'Polymorph', 'PhaseDiagram']

def banner_large():
    print("                                                  ")
    print("  \033[36m  _     _                _             \033[0m")
    print("  \033[36m | |   | |              | |            \033[0m")
    print("  \033[36m | |__ | | ___ _ __   __| |_ __  _   _ \033[0m")
    print("  \033[36m | '_ \\| |/ _ \\ '_ \\ / _` | '_ \\| | | |\033[0m")
    print("  \033[36m | |_) | |  __/ | | | (_| | |_) | |_| |\033[0m")
    print("  \033[36m |_.__/|_|\\___|_| |_|\\__,_| .__/ \\__, |\033[0m")
    print("  \033[36m                          | |     __/ |\033[0m")
    print("  \033[36m                          |_|    |___/ \033[0m")
    print("                                                  ")


def banner_small():
    print("                                   ")
    print("\033[36m     _   _           _             \033[0m")
    print("\033[36m    | |_| |___ ___ _| |___ _ _     \033[0m")
    print("\033[36m    | . | | -_|   | . | . | | |    \033[0m")
    print("\033[36m    |___|_|___|_|_|___|  _|_  |    \033[0m")
    print("\033[36m                      |_| |___|    \033[0m")
    print("                                   ")


banner_small()

print(f"\033[36m    version:\033[0m {__version__}                             ")
print("                                                  ")
print("\033[36m    developed by:\033[0m Leandro Seixas, PhD             ")
print("                                                  ")
print("-----------------------------------------------")
print("                                               ")
print("System:")
print(f"├── \033[36muser:\033[0m {os.environ['USER']}")
print(f"├── \033[36mhostname:\033[0m {gethostname()}")
print(f"├── \033[36mcwd:\033[0m {os.getcwd()}")
print(f"└── \033[36mdate:\033[0m {datetime.today().strftime("%Y-%m-%d, %H:%M:%S")}")
print("                                               ")
print("Python:")
print(f"├── \033[36mversion:\033[0m {__python_version__}      ")
print(f"└── \033[36mexecutable:\033[0m {__python_executable__}      ")
print("                                               ")
print("Dependencies:")
print(f"├── \033[36mase version:\033[0m {__ase_version__}")
print(f"├── \033[36mnumpy version:\033[0m {__numpy_version__}")
print(f"└── \033[36mpandas version:\033[0m {__pandas_version__}")
print("                                               ")
