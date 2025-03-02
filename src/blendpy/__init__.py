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

from ._version import __version__
from .alloy import Alloy
from .dsi_model import DSIModel
from .polymorph import Polymorph

from sys import version as __python_version__
from ase import __version__ as __ase_version__
from numpy import __version__ as __numpy_version__
from pandas import __version__ as __pandas_version__

__all__ = ['Alloy', 'DSIModel', 'Polymorph']


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
print(f"\033[36m    version:\033[0m {__version__}                             ")
print("                                                  ")
print("-----------------------------------------------")
print("                                               ")
print(f"\033[36mpython version:\033[0m {__python_version__}      ")
print("                                               ")
print("Dependencies:")
print(f"    \033[36mase version:\033[0m {__ase_version__}")
print(f"    \033[36mnumpy version:\033[0m {__numpy_version__}")
print(f"    \033[36mpandas version:\033[0m {__pandas_version__}")
print("                                               ")