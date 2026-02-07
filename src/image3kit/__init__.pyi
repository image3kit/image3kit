"""
The _core (PyBind11 wrapper) of image3 module of containing sirun and VxlImg submodules.
"""

from __future__ import annotations

from image3kit._core.sirun import Input
from image3kit._core.voxlib import VxlImgF32, VxlImgI32, VxlImgU8, VxlImgU16, cube, cylinder, sphere

__all__: list = [
    "Input",
    "VxlImgF32",
    "VxlImgI32",
    "VxlImgU8",
    "VxlImgU16",
    "__doc__",
    "__version__",
    "cube",
    "cylinder",
    "sphere",
]
__version__: str = "0.0.1"
