"""
Pybind11 example plugin
-----------------------

.. currentmodule:: voxlib

.. autosummary::
    :toctree: _generate

    VxlImgU8
    Input
    shape
    cylinder
    sphere
    ...
"""

class Input: pass

class shape:
    pass

class cylinder(shape):
    def __init__(self, p1: tuple[float], p2: tuple[float], inside: int) -> None:
        pass
class sphere(shape):
    def __init__(self, center: tuple[float], radius: float, inside: int) -> None:
        pass
class cube(shape):
    def __init__(self, p1: tuple[float], size: tuple[float], inside: int) -> None:
        pass

class VxlImgU8:
    """
    Voxel (3D) image of type uint8.
    """
    def __init__(self, n1: tuple[int], value: int) -> None:
        pass

    def Paint(self, sh: shape):
        pass

    def nx() -> int:
        pass
    def ny() -> int:
        pass
    def nz() -> int:
        pass
    def dx() -> float:
        pass

    def printInfo():
        pass
    def write():
        pass
    def writeNoHdr():
        pass
    def readFromHeader():
        pass
    def readAscii():
        pass
    def cropD():
        pass
    def growBox():
        pass
    def shrinkBox():
        pass
    def fillHoles():
        pass
    def threshold101():
        pass
    def writeAConnectedPoreVoxel():
        pass
    def AND():
        pass
    def NOT():
        pass
    def OR():
        pass
    def XOR():
        pass
    def resampleMode():
        pass
    def resampleMax():
        pass
    def resampleMean():
        pass
    def resliceZ():
        pass
    def Paint():
        pass
    def PaintAdd():
        pass
    def PaintBefore():
        pass
    def PaintAfter():
        pass
    def PaintAddBefore():
        pass
    def PaintAddAfter():
        pass
    def writeContour():
        pass
    def sliceToPng(normalAxis: str, fnam: str, iSlice: int, bgnv: int, endv: int, color: str):
        pass