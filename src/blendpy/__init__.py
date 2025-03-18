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
import platform
from datetime import datetime
from socket import gethostname
from sys import version as __python_version__
from sys import executable as __python_executable__
from ase import __version__ as __ase_version__
from ase import __file__ as __ase_file__
from numpy import __version__ as __numpy_version__
from numpy import __file__ as __numpy_file__
from pandas import __version__ as __pandas_version__
from pandas import __file__ as __pandas_file__
from pytest import __version__ as __pytest_version__
from pytest import __file__ as __pytest_file__
from yaml import __version__ as __yaml_version__
from yaml import __file__ as __yaml_file__
from tqdm import __version__ as __tqdm_version__
from tqdm import __file__ as __tqdm_file__
from mpltern import __version__ as __mpltern_version__
from mpltern import __file__ as __mpltern_file__

from ase.parallel import parprint as print
from ._version import __version__
from .constants import R, convert_eVatom_to_kJmol
from .alloy import Alloy
from .dsi_model import DSIModel
from .phase_diagram import PhaseDiagram
from .local_order import LocalOrder
# from .ace import ACE
# from .sqs import SQS
# from .intermetallics import Intermetallics

__all__ = ['Alloy', 'DSIModel', 'PhaseDiagram', 'LocalOrder']

ansi_colors = {'black': ["\033[30m", "\033[0m"],
               'red': ["\033[31m", "\033[0m"],
               'green': ["\033[32m", "\033[0m"],
               'yellow': ["\033[33m", "\033[0m"],
               'blue': ["\033[34m", "\033[0m"],
               'magenta': ["\033[35m", "\033[0m"],
               'cyan': ["\033[36m", "\033[0m"],
               'white': ["\033[37m", "\033[0m"],
               'empty': ["", ""],
               'bold': ["\033[1m", "\033[0m"],
               'underline': ["\033[4m", "\033[0m"],
               'blink': ["\033[5m", "\033[0m"],
               'reverse': ["\033[7m", "\033[0m"],
               'invisible': ["\033[8m", "\033[0m"]}


def starter(color='cyan'):
    color_start, color_end = ansi_colors[color]
    print("                                                  ")
    print(f"{color_start}     _   _           _             {color_end}")
    print(f"{color_start}    | |_| |___ ___ _| |___ _ _     {color_end}")
    print(f"{color_start}    | . | | -_|   | . | . | | |    {color_end}")
    print(f"{color_start}    |___|_|___|_|_|___|  _|_  |    {color_end}")
    print(f"{color_start}                      |_| |___|    {color_end}")
    print("                                                  ")
    print(f"{color_start}    version:{color_end} {__version__}                              ")
    print(f"{color_start}    developed by:{color_end} Leandro Seixas, PhD             ")
    print(f"{color_start}    homepage:{color_end} https://github.com/leseixas/blendpy")
    print("                                                  ")
    print("--------------------------------------------------")
    print("                                                  ")
    print("System:")
    print(f"├── {color_start}architecture:{color_end} {platform.machine()}")
    print(f"├── {color_start}platform:{color_end} {platform.system()}")
    print(f"├── {color_start}user:{color_end} {os.environ['USER']}")
    print(f"├── {color_start}hostname:{color_end} {gethostname()}")
    print(f"├── {color_start}cwd:{color_end} {os.getcwd()}")
    print(f"├── {color_start}date:{color_end} {datetime.today().strftime("%Y-%m-%d, %H:%M:%S")}")
    print(f"└── {color_start}PID:{color_end} {os.getpid()}")
    print("                                               ")
    print("Python:")
    print(f"├── {color_start}version:{color_end} {__python_version__}      ")
    print(f"└── {color_start}executable:{color_end} {__python_executable__}      ")
    print("                                               ")
    print("Dependencies:")
    print(f"├── {color_start}ase version:{color_end} {__ase_version__}    [{__ase_file__[:-11]}]")
    print(f"├── {color_start}numpy version:{color_end} {__numpy_version__}    [{__numpy_file__[:-11]}]")
    print(f"├── {color_start}pandas version:{color_end} {__pandas_version__}    [{__pandas_file__[:-11]}]")
    print(f"├── {color_start}pytest version:{color_end} {__pytest_version__}    [{__pytest_file__[:-11]}]")
    print(f"├── {color_start}yaml version:{color_end} {__yaml_version__}    [{__yaml_file__[:-11]}]")
    print(f"├── {color_start}tqdm version:{color_end} {__tqdm_version__}    [{__tqdm_file__[:-11]}]")
    print(f"└── {color_start}mpltern version:{color_end} {__mpltern_version__}    [{__mpltern_file__[:-11]}]")
    print("                                               ")


starter()