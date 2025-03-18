import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import pytest
import numpy as np
from ase import Atoms
from blendpy.local_order import LocalOrder


