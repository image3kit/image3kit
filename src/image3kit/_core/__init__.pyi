"""
The _core (PyBind11 wrapper) of image3 module of containing sirun and VxlImg submodules.
"""
from __future__ import annotations
from . import sirun
from . import voxlib
__all__: list[str] = ['sirun', 'voxlib']
__version__: str = '0.0.1'
