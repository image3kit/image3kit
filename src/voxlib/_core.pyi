"""
 -----------------------

        .. currentmodule:: VoxelImage

        .. autosummary::
           :toctree: _generate

           add
           subtract
    
"""
from __future__ import annotations
import collections.abc
import numpy
import numpy.typing
import typing
__all__: list[str] = ['Input', 'VxlImgU8', 'cube', 'cylinder', 'dbl3', 'int3', 'readImage', 'shape', 'sphere', 'voxelImageTBase']
class Input:
    def __init__(self, arg0: dict) -> None:
        ...
    def add(self, arg0: str, arg1: str) -> None:
        ...
    def echoKeywords(self) -> None:
        ...
    def get(self, arg0: str, arg1: str) -> str:
        ...
    def renameKeys(self, arg0: str, arg1: str) -> None:
        ...
    def set(self, arg0: str, arg1: str) -> None:
        ...
    def setDefault(self, arg0: str, arg1: str) -> None:
        ...
class VxlImgU8:
    def AND(self, other: VxlImgU8) -> None:
        """
        Pixel-wise AND operation.
        """
    def FaceMedian06(self, arg0: typing.SupportsInt, arg1: typing.SupportsInt) -> int:
        ...
    def NOT(self, arg0: VxlImgU8) -> None:
        """
        Pixel-wise NOT operation.
        """
    def OR(self, other: VxlImgU8) -> None:
        """
        Pixel-wise OR operation.
        """
    def Paint(self, shape: shape) -> None:
        """
        Paint a shape into the image.
        """
    def PaintAdd(self, shape: shape) -> None:
        """
        Add a shape's values to the image.
        """
    def PaintAddAfter(self, shape: shape) -> None:
        """
        Add paint after...?
        """
    def PaintAddBefore(self, shape: shape) -> None:
        """
        Add paint before...?
        """
    def PaintAfter(self, shape: shape) -> None:
        """
        Paint after...?
        """
    def PaintBefore(self, shape: shape) -> None:
        """
        Paint before...?
        """
    def PointMedian032(self, arg0: typing.SupportsInt, arg1: typing.SupportsInt, arg2: typing.SupportsInt, arg3: typing.SupportsInt) -> None:
        ...
    def XOR(self, other: VxlImgU8) -> None:
        """
        Pixel-wise XOR operation.
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    def __init__(self, tpl: tuple, value: typing.SupportsInt = 0) -> None:
        """
        Initialize from a size tuple (nz, ny, nx) with an optional fill value.
        """
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    def addSurfNoise(self, mask1: typing.SupportsInt, mask2: typing.SupportsInt, threshold: typing.SupportsInt, seed: typing.SupportsInt) -> None:
        """
        Add surface noise.
        """
    def averageWith(self, arg0: str) -> bool:
        ...
    def averageWith_mBE(self, arg0: str) -> bool:
        ...
    def circleOut(self, x: typing.SupportsInt, y: typing.SupportsInt, r: typing.SupportsInt, d: str, val: typing.SupportsInt) -> None:
        """
        Circle out operation.
        """
    def cropD(self, begin: int3, end: int3, emptyLayers: typing.SupportsInt = 0, emptyLayersValue: typing.SupportsInt = 1, verbose: bool = False) -> None:
        """
        Crop the image by a specified depth.
        """
    def data(self) -> numpy.typing.NDArray[numpy.uint8]:
        """
        Get the raw data buffer as a numpy array.
        """
    def delense032(self, x: typing.SupportsInt, y: typing.SupportsInt, r: typing.SupportsInt, d: str, val: typing.SupportsInt) -> None:
        """
        Delense operation.
        """
    def direction(self, arg0: str) -> None:
        """
        Get direction?
        """
    def faceMedNgrowToFrom(self, labelTo: typing.SupportsInt, labelFrom: typing.SupportsInt, nDiff: typing.SupportsInt) -> None:
        """
        Face median grow to/from labels.
        """
    def fillHoles(self, maxHoleRadius: typing.SupportsInt) -> None:
        """
        Fill closed holes in the image.
        """
    def grow0(self) -> None:
        """
        Grow pore phase (0).
        """
    def growBox(self, layers: typing.SupportsInt) -> None:
        """
        Expand the image boundaries.
        """
    def growLabel(self, arg0: typing.SupportsInt) -> None:
        ...
    def keepLargest0(self) -> None:
        ...
    def mapFrom(self, sourceImage: VxlImgU8, vmin: typing.SupportsInt, vmax: typing.SupportsInt, scale: typing.SupportsFloat, shift: typing.SupportsFloat) -> None:
        """
        Map values from another image.
        """
    def medianFilter(self) -> None:
        """
        Apply median filter.
        """
    def medianX(self) -> None:
        """
        Apply median X filter.
        """
    def modeNSames(self, nSameNeighbors: typing.SupportsInt) -> int:
        """
        Apply mode filter based on neighbor count.
        """
    def nx(self) -> int:
        ...
    def ny(self) -> int:
        ...
    def nz(self) -> int:
        ...
    def plotAll(self, arg0: str) -> bool:
        ...
    def printInfo(self) -> None:
        ...
    def rangeTo(self, min: typing.SupportsInt, max: typing.SupportsInt, val: typing.SupportsInt) -> None:
        """
        Set values in range [min, max] to val.
        """
    def read(self, filename: str) -> None:
        """
        Read image from file.
        """
    def readAscii(self, filename: str) -> bool:
        """
        Read image data from an ASCII file.
        """
    def readFromHeader(self, header_file: str, processKeys: typing.SupportsInt = 1) -> None:
        """
        Read image dimensions/metadata from a header file.
        """
    def redirect(self, arg0: str) -> None:
        """
        Rotate/Redirect image.
        """
    def replaceRange(self, min: typing.SupportsInt, max: typing.SupportsInt, val: typing.SupportsInt) -> None:
        """
        Replace values in range [min, max] with val.
        """
    def resample(self, factor: typing.SupportsFloat) -> None:
        """
        Resample image by factor f.
        """
    def resampleMax(self, rate: typing.SupportsFloat) -> VxlImgU8:
        """
        Resample the image using max interpolation.
        """
    def resampleMean(self, rate: typing.SupportsFloat) -> VxlImgU8:
        """
        Resample the image using mean interpolation.
        """
    def resampleMode(self, rate: typing.SupportsFloat) -> VxlImgU8:
        """
        Resample the image using mode interpolation.
        """
    def rescaleValues(self, min: typing.SupportsInt, max: typing.SupportsInt) -> None:
        """
        Rescale image values to [min, max].
        """
    def resliceZ(self, rate: typing.SupportsFloat) -> VxlImgU8:
        """
        Reslice along Z axis.
        """
    def segment2(self, nSegs: typing.SupportsInt = 2, th: collections.abc.Sequence[typing.SupportsInt] = [], minSizs: collections.abc.Sequence[typing.SupportsInt] = [], noisev: typing.SupportsFloat = 2.0, localF: typing.SupportsFloat = 800.0, flatnes: typing.SupportsFloat = 0.1, resolution: typing.SupportsFloat = 2.0, gradFactor: typing.SupportsFloat = 0.0, krnl: typing.SupportsInt = 2, nItrs: typing.SupportsInt = 13, writedumps: typing.SupportsInt = 0) -> bool:
        ...
    def setOffset(self, offset: dbl3) -> None:
        """
        Set the spatial offset (origin).
        """
    def shape(self) -> tuple:
        ...
    def shrink0(self) -> None:
        """
        Shrink pore phase (0).
        """
    @typing.overload
    def shrinkBox(self, layers: typing.SupportsInt) -> None:
        """
        Shrink the image boundaries.
        """
    @typing.overload
    def shrinkBox(self, layers: typing.SupportsInt) -> None:
        """
        Shrink the image boundaries.
        """
    def sliceToPng(self, normalAxis: str, filename: str, sliceIndex: typing.SupportsInt, val_min: typing.SupportsInt, val_max: typing.SupportsInt, color_map: str = 'gray') -> None:
        """
        Save a 2D slice as a PNG image.
        """
    def threshold101(self, min: typing.SupportsInt, max: typing.SupportsInt) -> None:
        """
        Apply a threshold to convert to 0/1.
        """
    def write(self, filename: str) -> None:
        """
        Write the image to a file.
        """
    def write8bit(self, filename: str, min: typing.SupportsFloat, max: typing.SupportsFloat) -> None:
        """
        Write as 8-bit image scaled between min and max.
        """
    def writeAConnectedPoreVoxel(self, filename: str) -> None:
        """
        Write a specific connected pore voxel to file.
        """
    def writeContour(self, outSurf: str) -> None:
        """
        Write contour extraction to a surface file.
        """
    def writeNoHdr(self, filename: str) -> None:
        """
        Write the raw image data without a header.
        """
class cube(shape):
    def __init__(self, arg0: tuple, arg1: tuple, arg2: typing.SupportsInt) -> None:
        ...
class cylinder(shape):
    def __init__(self, arg0: tuple, arg1: tuple, arg2: typing.SupportsFloat, arg3: typing.SupportsInt) -> None:
        ...
class dbl3:
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: tuple) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def x(self) -> float:
        ...
    @x.setter
    def x(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def y(self) -> float:
        ...
    @y.setter
    def y(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def z(self) -> float:
        ...
    @z.setter
    def z(self, arg0: typing.SupportsFloat) -> None:
        ...
class int3:
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: typing.SupportsInt, arg1: typing.SupportsInt, arg2: typing.SupportsInt) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: tuple) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def x(self) -> int:
        ...
    @x.setter
    def x(self, arg0: typing.SupportsInt) -> None:
        ...
    @property
    def y(self) -> int:
        ...
    @y.setter
    def y(self, arg0: typing.SupportsInt) -> None:
        ...
    @property
    def z(self) -> int:
        ...
    @z.setter
    def z(self, arg0: typing.SupportsInt) -> None:
        ...
class shape:
    pass
class sphere(shape):
    def __init__(self, arg0: tuple, arg1: typing.SupportsFloat, arg2: typing.SupportsInt) -> None:
        ...
class voxelImageTBase:
    pass
def readImage(filename: str, processKeys: typing.SupportsInt = 1) -> voxelImageTBase:
    """
    Global helper to read an image from a file.
    """
__version__: str = '0.0.1'
