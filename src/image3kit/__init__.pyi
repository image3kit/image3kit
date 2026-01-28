"""
The _core (PyBind11 wrapper) of image3 module of containing sirun and VxlImg submodules.
"""
from __future__ import annotations
from image3kit._core.sirun import Input
from image3kit._core.voxlib import VxlImgU8
from image3kit._core.voxlib import cube
from image3kit._core.voxlib import cylinder
from image3kit._core.voxlib import readImage
from image3kit._core.voxlib import readImageU8
from image3kit._core.voxlib import sphere
from . import _core
__all__: list = ['__doc__', '__version__', 'VxlImgU8', 'Input', 'sphere', 'cylinder', 'cube', 'readImage', 'readImageU8']
__version__: str = '0.0.1'
