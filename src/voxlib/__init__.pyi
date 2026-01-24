"""
 -----------------------

        .. currentmodule:: VoxelImage

        .. autosummary::
           :toctree: _generate

           add
           subtract
    
"""
from __future__ import annotations
from voxlib._core import Input
from voxlib._core import VxlImgU8
from voxlib._core import cube
from voxlib._core import cylinder
from voxlib._core import sphere
from . import _core
__all__: list = ['__doc__', '__version__', 'VxlImgU8', 'Input', 'sphere', 'cylinder', 'cube']
__version__: str = '0.0.1'
