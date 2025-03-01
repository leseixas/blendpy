# -*- coding: utf-8 -*-
# file: __init__.py

from ._version import __version__
from .alloy import Alloy
from .dsi_model import DSIModel
from .polymorph import Polymorph

__all__ = ['Alloy', 'DSIModel', 'Polymorph']

print("                                                ")
print("\033[36m   _      _                   _                 \033[0m")
print("\033[36m  | |__  | |  ___  _ __    __| | _ __   _   _   \033[0m")
print("\033[36m  | '_ \\ | | / _ \\| '_ \\  / _` || '_ \\ | | | |  \033[0m")
print("\033[36m  | |_) || ||  __/| | | || (_| || |_) || |_| |  \033[0m")
print("\033[36m  |_.__/ |_| \\___||_| |_| \\__,_|| .__/  \\__, |  \033[0m")
print("\033[36m                                |_|     |___/   \033[0m")
print("\033[36m                                                \033[0m")
print(f"                 version: {__version__}                 ")
print("                                                ")