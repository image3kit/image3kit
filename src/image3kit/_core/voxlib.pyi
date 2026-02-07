"""
The VxlImg template classes
"""

from __future__ import annotations

import collections.abc
import typing

import numpy
import numpy.typing

import image3kit._core.sirun

__all__: list[str] = [
    "VxlImgF32",
    "VxlImgI32",
    "VxlImgU8",
    "VxlImgU16",
    "cube",
    "cylinder",
    "labelImage",
    "readImageBase",
    "shape",
    "sphere",
    "voxelImageTBase",
]

class VxlImgF32:
    def AND(self, image2: VxlImgF32) -> None:
        """
        Voxel-by-voxel inplace AND operation.
        """
    def FaceMedian06(self, nAdj0: typing.SupportsInt, nAdj1: typing.SupportsInt) -> int:
        """
        Set voxel value to 0/1 if it has more than nAdj0/1 neighbours with value 0/1, in its 6 nearest voxels
        """
    def NOT(self, image2: VxlImgF32) -> None:
        """
        Voxel-by-voxel inplace NOT operation, alias for img = img & !img2.
        """
    def OR(self, image2: VxlImgF32) -> None:
        """
        Voxel-by-voxel inplace OR operation.
        """
    def XOR(self, image2: VxlImgF32) -> None:
        """
        Voxel-by-voxel inplace XOR operation.
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def __init__(self, shape: tuple = (0, 0, 0), value: typing.SupportsFloat = 0) -> None:
        """
        Initialize a new image of size tuple (nx, ny, nz) with the fill value.
        """
    @typing.overload
    def __init__(self, filepath: typing.Any, processKeys: bool = True) -> None:
        """
        Read image dimensions/metadata from a (header) file. SUpported file types are .am, .raw
        """
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    def __repr__(self) -> str: ...
    def addSurfNoise(
        self,
        mask1: typing.SupportsInt,
        mask2: typing.SupportsInt,
        threshold: typing.SupportsInt,
        seed: typing.SupportsInt,
    ) -> None:
        """
        Add surface noise.
        """
    def adjustBrightnessWith(self, image_file: str) -> bool: ...
    def adjustSliceBrightness(
        self,
        mask_a: VxlImgU8,
        mask_b: VxlImgU8,
        ref_image: VxlImgF32,
        smooth_iter: typing.SupportsInt = 3,
        smooth_kernel: typing.SupportsInt = 20,
    ) -> bool: ...
    def averageWith(self, arg0: str) -> bool: ...
    def averageWith_mBE(self, arg0: str) -> bool: ...
    def bilateralGauss(
        self,
        iterations: typing.SupportsInt = 1,
        kernel_radius: typing.SupportsInt = 1,
        sigma_val: typing.SupportsFloat = 16.0,
        sharpness: typing.SupportsFloat = 0.1,
        sigma_spatial: typing.SupportsFloat = 2.0,
    ) -> bool: ...
    def bilateralX(
        self,
        iterations: typing.SupportsInt = 1,
        kernel_radius: typing.SupportsInt = 1,
        x_step: typing.SupportsInt = 2,
        sigma_val: typing.SupportsFloat = 16.0,
        sharpness: typing.SupportsFloat = 0.1,
        sigma_spatial: typing.SupportsFloat = 2.0,
    ) -> bool:
        """
        Bilateral filter with Xtra large kernel radius, actual kernel size is: kernel_radius * x_step cubed.
        """
    def circleOut(
        self, x: typing.SupportsInt, y: typing.SupportsInt, r: typing.SupportsInt, d: str, val: typing.SupportsFloat
    ) -> None:
        """
        Circle out operation.
        """
    def copy(self) -> VxlImgF32:
        """
        duplicate the image data
        """
    def cropD(
        self,
        begin: tuple,
        end: tuple,
        emptyLayers: typing.SupportsInt = 0,
        emptyLayersValue: typing.SupportsFloat = 1,
        verbose: bool = False,
    ) -> None:
        """
        Crop the image (inplace) from begin index tupe ix,iy,iz (inclusive) to and and end index (not inclusive) tuple.
        """
    def cutOutside(
        self,
        axis: str = "z",
        extra_out: typing.SupportsInt = 0,
        threshold: typing.SupportsInt = -1,
        cut_highs: typing.SupportsInt = 0,
        shift_x: typing.SupportsInt = 0,
        shift_y: typing.SupportsInt = 0,
        fill_val: typing.SupportsInt = 0,
    ) -> None: ...
    def data(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Get the raw data buffer as a numpy array.
        """
    def delense032(
        self,
        iterations: typing.SupportsInt,
        nAdj0: typing.SupportsInt,
        nAdj1: typing.SupportsInt,
        lbl0: typing.SupportsFloat,
        lbl1: typing.SupportsFloat,
    ) -> None:
        """
        Delense operation.
        """
    def dering(
        self,
        x0: typing.SupportsInt,
        y0: typing.SupportsInt,
        x1: typing.SupportsInt,
        y1: typing.SupportsInt,
        min_val: typing.SupportsInt = 0,
        max_val: typing.SupportsInt = 255,
        nr: typing.SupportsInt = 0,
        ntheta: typing.SupportsInt = 18,
        nz: typing.SupportsInt = 0,
        scale_dif_val: typing.SupportsFloat = 1,
        bilateral_sharpen: typing.SupportsFloat = 0.05,
        nGrowBox: typing.SupportsInt = 10,
        write_dumps: bool = True,
    ) -> None: ...
    def direction(self, arg0: str) -> None:
        """
        Get direction?
        """
    def faceMedNgrowToFrom(
        self, labelTo: typing.SupportsFloat, labelFrom: typing.SupportsFloat, nDiff: typing.SupportsInt
    ) -> None:
        """
        Face median grow to/from labels.
        """
    def fillHoles(self, maxHoleRadius: typing.SupportsInt) -> None:
        """
        Fill closed holes in the image.
        """
    def flipEndian(self) -> None: ...
    def grow0(self) -> None:
        """
        Grow pore phase (voxel values of 0).
        """
    def growBox(self, num_layers: typing.SupportsInt) -> None:
        """
        Expand the image boundaries, increasing its size by `num_layers` in all directions
        """
    def growLabel(self, arg0: typing.SupportsFloat) -> None: ...
    def growingThreshold(
        self,
        startMin: typing.SupportsFloat,
        startMax: typing.SupportsFloat,
        finalMin: typing.SupportsFloat,
        finalMax: typing.SupportsFloat,
        iterations: typing.SupportsInt = 4,
    ) -> None: ...
    def keepLargest(self, min: typing.SupportsFloat, max: typing.SupportsFloat) -> None:
        """
        Keep largest singly-connected region with values in [min, max].
        """
    def mapFrom(
        self,
        sourceImage: VxlImgF32,
        vmin: typing.SupportsFloat,
        vmax: typing.SupportsFloat,
        scale: typing.SupportsFloat,
        shift: typing.SupportsFloat,
    ) -> None:
        """
        Map values from another image.
        """
    def meanWide(
        self,
        width: typing.SupportsInt = 0,
        noise_val: typing.SupportsInt = 4,
        average: typing.SupportsInt = 0,
        delta: typing.SupportsInt = 20,
        iterations: typing.SupportsInt = 15,
        smooth_image: str = "",
    ) -> bool:
        """
        computes a background image, used to correct for lens artifacts
        """
    def medianFilter(self) -> None:
        """
        Apply a 1+6-neighbour median filter.
        """
    def medianX(self) -> None:
        """
        Apply median filter with kernel size of 1 voxels in x-direction
        """
    def medianY(self) -> None:
        """
        Apply median filter with kernel size of 1 voxels in y-direction
        """
    def medianZ(self) -> None:
        """
        Apply median filter with kernel size of 1 voxels in z-direction
        """
    def mode26(self, arg0: typing.SupportsInt) -> None: ...
    def modeNSames(self, nSameNeighbors: typing.SupportsInt) -> int:
        """
        Apply mode filter based on nearest 6 neighbor voxels.
        """
    def nx(self) -> int: ...
    def ny(self) -> int: ...
    def nz(self) -> int: ...
    def origin(self) -> tuple:
        """
        Get the origin value (x0, y0, z0).
        """
    def otsu_threshold(
        self, min_val: typing.SupportsInt = 0, max_val: typing.SupportsInt = 256
    ) -> typing.Annotated[list[float], "FixedSize(5)"]: ...
    def paint(self, shape: shape) -> None:
        """
        Paint (set values of) a shape into the image.
        """
    def paintAdd(self, shape: shape) -> None:
        """
        Add (+) a shape's value to the image.
        """
    def paintAddAfter(self, shape: shape) -> None:
        """
        Add (+) a shape's value after the shape (plane...)
        """
    def paintAddBefore(self, shape: shape) -> None:
        """
        Add (+) a shape's value before the shape (plane...)
        """
    def paintAfter(self, shape: shape) -> None:
        """
        Paint after the shape (plane...)
        """
    def paintBefore(self, shape: shape) -> None:
        """
        Paint before the shape (plane...)
        """
    def plotAll(
        self,
        filename: str = "pltAll",
        min_val: typing.SupportsInt = 0,
        max_val: typing.SupportsInt = -1000001,
        slice_index: typing.SupportsInt = -1000000,
        histogram_bins: typing.SupportsInt = 128,
        normal_axis: str = "xyz",
        grey: bool = True,
        color: bool = True,
        histogram: bool = True,
        z_profile: bool = True,
        alpha_image: VxlImgF32 = None,
        alpha_min: typing.SupportsInt = 0,
        alpha_max: typing.SupportsInt = -1000001,
    ) -> bool:
        """
        Plot all visualizations (Histogram, ZProfile, Slices) with various options, hackish for debugging
        """
    def plotHistogramSvg(
        self,
        filename: str = "aa.svg",
        bins: typing.SupportsInt = 128,
        min_val: typing.SupportsFloat = 3e38,
        max_val: typing.SupportsFloat = -3e38,
    ) -> int: ...
    def plotSlice(
        self,
        filename: str,
        normal_axis: str = "xyz",
        slice_index: typing.SupportsInt = -1000000,
        min_val: typing.SupportsInt = 0,
        max_val: typing.SupportsInt = -1000001,
        color_map: str = "gray",
    ) -> None:
        """
        Save a 2D slice as a PNG image.
        normalAxis can be 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'.
        color_map can be 'gray' or 'RGB'
        """
    def plotZProfileSvg(
        self, filename: str = "aa.svg", min_val: typing.SupportsFloat = 0, max_val: typing.SupportsFloat = 255
    ) -> int: ...
    def pointMedian032(
        self,
        nAdj0: typing.SupportsInt,
        nAdj1: typing.SupportsInt,
        lbl0: typing.SupportsFloat,
        lbl1: typing.SupportsFloat,
    ) -> None:
        """
        Set voxel value to lbl0/1 if it has more than nAdj0/1 neighbours with value lbl0/1, in its 6+26 nearest voxels
        """
    def printInfo(self) -> None: ...
    def rangeTo(self, min: typing.SupportsFloat, max: typing.SupportsFloat, val: typing.SupportsFloat) -> None:
        """
        Set values in range [min, max] to val.
        """
    def readAscii(self, filename: str) -> bool:
        """
        Read image data from an ASCII file.
        """
    def readBin(self, filename: str, nSkipBytes: typing.SupportsInt = 0) -> None:
        """
        Read image data from  a .raw, .raw/.raw.gz, or reset and read from a .am or .tif file.
        """
    def readFromFloat(
        self, header: str, scale: typing.SupportsFloat = 1.0, shift: typing.SupportsFloat = 0.0
    ) -> bool: ...
    def readFromHeader(self, filename: str) -> None:
        """
        Reset and read image from file.
        """
    def redirect(self, arg0: str) -> None:
        """
        Swap X axis with the given axis (y or z).
        """
    def replaceByImageRange(
        self, min_val: typing.SupportsFloat, max_val: typing.SupportsFloat, image_file: str
    ) -> None: ...
    def replaceOutSideValue(
        self, val_old: typing.SupportsInt = 0, val_new: typing.SupportsInt = 2, hole_size: typing.SupportsInt = 5
    ) -> bool: ...
    def replaceRange(self, min: typing.SupportsFloat, max: typing.SupportsFloat, val: typing.SupportsFloat) -> None:
        """
        Replace values in range [min, max] with val.
        """
    def replaceRangeByImage(
        self, min_val: typing.SupportsFloat, max_val: typing.SupportsFloat, image_file: str
    ) -> None: ...
    def resample(self, factor: typing.SupportsFloat) -> None:
        """
        Resample image by factor f, using averaging (downsampling, f>1) or nearest when upsampling (f<1)
        """
    def resampleMax(self, factor: typing.SupportsFloat) -> VxlImgF32:
        """
        Downsample the image, setting voxel values to maximum of original encompassing voxel values.
        """
    def resampleMean(self, factor: typing.SupportsFloat) -> VxlImgF32:
        """
        Resample the image using mean interpolation.
        """
    def resampleMode(self, factor: typing.SupportsFloat) -> VxlImgF32:
        """
        Downsample the image, setting voxel values to mode of original encompassing voxel values.
        """
    def rescaleValues(self, min: typing.SupportsFloat, max: typing.SupportsFloat) -> None:
        """
        Rescale image values to [min, max].
        """
    def resliceZ(self, factor: typing.SupportsFloat) -> VxlImgF32:
        """
        Reslice along the Z axis.
        """
    def scaleDx(self, scale: typing.SupportsFloat) -> None:
        """
        Scale the voxel size (dx, dy, dz) and origin by a factor.
        """
    def setOffset(self, offset: image3kit._core.sirun.dbl3) -> None:
        """
        Set the spatial offset (origin).
        """
    def setOrigin(self, origin: tuple) -> None:
        """
        Set the origin value (x0, y0, z0).
        """
    def setVoxelSize(self, voxelSize: tuple) -> None:
        """
        Set the voxel size (dx, dy, dz).
        """
    def shape(self) -> tuple: ...
    def shrink0(self) -> None:
        """
        Shrink pore phase (voxel values of 0).
        """
    def shrinkBox(self, num_layers: typing.SupportsInt) -> None:
        """
        Shrink the image boundaries, decreasing its size by the given num_layers in all directions
        """
    def smooth(
        self,
        iterations: typing.SupportsInt = 1,
        kernel_radius: typing.SupportsInt = 1,
        sigma_val: typing.SupportsFloat = 16.0,
        sharpness: typing.SupportsFloat = 0.1,
    ) -> bool:
        """
        bilateral smoothing filter
        """
    def threshold101(self, min: typing.SupportsFloat, max: typing.SupportsFloat) -> None:
        """
        Apply a threshold to binarize the image, set voxel-values to convert to 0 in between the min and max thresholds and 1 outside of it
        """
    def variance(self, min_val: typing.SupportsInt = 0, max_val: typing.SupportsInt = 255) -> float:
        """
        Set outer tubing of a circular core-holder image to fill_val
        """
    def voxelSize(self) -> tuple:
        """
        Get the voxel size (dx, dy, dz).
        """
    def write(self, filename: str) -> None:
        """
        Write the image to a file (.mhd, .raw, .ra.gz formats).
        """
    def write8bit(self, filename: str, min: typing.SupportsFloat = 0.0, max: typing.SupportsFloat = -0.5) -> None:
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
    def writeNoHeader(self, filename: str) -> None:
        """
        Write the raw image data without a header.
        """

class VxlImgI32:
    def AND(self, image2: VxlImgI32) -> None:
        """
        Voxel-by-voxel inplace AND operation.
        """
    def FaceMedian06(self, nAdj0: typing.SupportsInt, nAdj1: typing.SupportsInt) -> int:
        """
        Set voxel value to 0/1 if it has more than nAdj0/1 neighbours with value 0/1, in its 6 nearest voxels
        """
    def NOT(self, image2: VxlImgI32) -> None:
        """
        Voxel-by-voxel inplace NOT operation, alias for img = img & !img2.
        """
    def OR(self, image2: VxlImgI32) -> None:
        """
        Voxel-by-voxel inplace OR operation.
        """
    def XOR(self, image2: VxlImgI32) -> None:
        """
        Voxel-by-voxel inplace XOR operation.
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def __init__(self, shape: tuple = (0, 0, 0), value: typing.SupportsInt = 0) -> None:
        """
        Initialize a new image of size tuple (nx, ny, nz) with the fill value.
        """
    @typing.overload
    def __init__(self, filepath: typing.Any, processKeys: bool = True) -> None:
        """
        Read image dimensions/metadata from a (header) file. SUpported file types are .am, .raw
        """
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    def __repr__(self) -> str: ...
    def addSurfNoise(
        self,
        mask1: typing.SupportsInt,
        mask2: typing.SupportsInt,
        threshold: typing.SupportsInt,
        seed: typing.SupportsInt,
    ) -> None:
        """
        Add surface noise.
        """
    def adjustBrightnessWith(self, image_file: str) -> bool: ...
    def adjustSliceBrightness(
        self,
        mask_a: VxlImgU8,
        mask_b: VxlImgU8,
        ref_image: VxlImgI32,
        smooth_iter: typing.SupportsInt = 3,
        smooth_kernel: typing.SupportsInt = 20,
    ) -> bool: ...
    def averageWith(self, arg0: str) -> bool: ...
    def averageWith_mBE(self, arg0: str) -> bool: ...
    def bilateralGauss(
        self,
        iterations: typing.SupportsInt = 1,
        kernel_radius: typing.SupportsInt = 1,
        sigma_val: typing.SupportsFloat = 16.0,
        sharpness: typing.SupportsFloat = 0.1,
        sigma_spatial: typing.SupportsFloat = 2.0,
    ) -> bool: ...
    def bilateralX(
        self,
        iterations: typing.SupportsInt = 1,
        kernel_radius: typing.SupportsInt = 1,
        x_step: typing.SupportsInt = 2,
        sigma_val: typing.SupportsFloat = 16.0,
        sharpness: typing.SupportsFloat = 0.1,
        sigma_spatial: typing.SupportsFloat = 2.0,
    ) -> bool:
        """
        Bilateral filter with Xtra large kernel radius, actual kernel size is: kernel_radius * x_step cubed.
        """
    def circleOut(
        self, x: typing.SupportsInt, y: typing.SupportsInt, r: typing.SupportsInt, d: str, val: typing.SupportsInt
    ) -> None:
        """
        Circle out operation.
        """
    def copy(self) -> VxlImgI32:
        """
        duplicate the image data
        """
    def cropD(
        self,
        begin: tuple,
        end: tuple,
        emptyLayers: typing.SupportsInt = 0,
        emptyLayersValue: typing.SupportsInt = 1,
        verbose: bool = False,
    ) -> None:
        """
        Crop the image (inplace) from begin index tupe ix,iy,iz (inclusive) to and and end index (not inclusive) tuple.
        """
    def cutOutside(
        self,
        axis: str = "z",
        extra_out: typing.SupportsInt = 0,
        threshold: typing.SupportsInt = -1,
        cut_highs: typing.SupportsInt = 0,
        shift_x: typing.SupportsInt = 0,
        shift_y: typing.SupportsInt = 0,
        fill_val: typing.SupportsInt = 0,
    ) -> None: ...
    def data(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        Get the raw data buffer as a numpy array.
        """
    def delense032(
        self,
        iterations: typing.SupportsInt,
        nAdj0: typing.SupportsInt,
        nAdj1: typing.SupportsInt,
        lbl0: typing.SupportsInt,
        lbl1: typing.SupportsInt,
    ) -> None:
        """
        Delense operation.
        """
    def dering(
        self,
        x0: typing.SupportsInt,
        y0: typing.SupportsInt,
        x1: typing.SupportsInt,
        y1: typing.SupportsInt,
        min_val: typing.SupportsInt = 0,
        max_val: typing.SupportsInt = 255,
        nr: typing.SupportsInt = 0,
        ntheta: typing.SupportsInt = 18,
        nz: typing.SupportsInt = 0,
        scale_dif_val: typing.SupportsFloat = 1,
        bilateral_sharpen: typing.SupportsFloat = 0.05,
        nGrowBox: typing.SupportsInt = 10,
        write_dumps: bool = True,
    ) -> None: ...
    def direction(self, arg0: str) -> None:
        """
        Get direction?
        """
    def faceMedNgrowToFrom(
        self, labelTo: typing.SupportsInt, labelFrom: typing.SupportsInt, nDiff: typing.SupportsInt
    ) -> None:
        """
        Face median grow to/from labels.
        """
    def fillHoles(self, maxHoleRadius: typing.SupportsInt) -> None:
        """
        Fill closed holes in the image.
        """
    def flipEndian(self) -> None: ...
    def grow0(self) -> None:
        """
        Grow pore phase (voxel values of 0).
        """
    def growBox(self, num_layers: typing.SupportsInt) -> None:
        """
        Expand the image boundaries, increasing its size by `num_layers` in all directions
        """
    def growLabel(self, arg0: typing.SupportsInt) -> None: ...
    def growingThreshold(
        self,
        startMin: typing.SupportsInt,
        startMax: typing.SupportsInt,
        finalMin: typing.SupportsInt,
        finalMax: typing.SupportsInt,
        iterations: typing.SupportsInt = 4,
    ) -> None: ...
    def keepLargest(self, min: typing.SupportsInt, max: typing.SupportsInt) -> None:
        """
        Keep largest singly-connected region with values in [min, max].
        """
    def mapFrom(
        self,
        sourceImage: VxlImgI32,
        vmin: typing.SupportsInt,
        vmax: typing.SupportsInt,
        scale: typing.SupportsFloat,
        shift: typing.SupportsFloat,
    ) -> None:
        """
        Map values from another image.
        """
    def meanWide(
        self,
        width: typing.SupportsInt = 0,
        noise_val: typing.SupportsInt = 4,
        average: typing.SupportsInt = 0,
        delta: typing.SupportsInt = 20,
        iterations: typing.SupportsInt = 15,
        smooth_image: str = "",
    ) -> bool:
        """
        computes a background image, used to correct for lens artifacts
        """
    def medianFilter(self) -> None:
        """
        Apply a 1+6-neighbour median filter.
        """
    def medianX(self) -> None:
        """
        Apply median filter with kernel size of 1 voxels in x-direction
        """
    def medianY(self) -> None:
        """
        Apply median filter with kernel size of 1 voxels in y-direction
        """
    def medianZ(self) -> None:
        """
        Apply median filter with kernel size of 1 voxels in z-direction
        """
    def mode26(self, arg0: typing.SupportsInt) -> None: ...
    def modeNSames(self, nSameNeighbors: typing.SupportsInt) -> int:
        """
        Apply mode filter based on nearest 6 neighbor voxels.
        """
    def nx(self) -> int: ...
    def ny(self) -> int: ...
    def nz(self) -> int: ...
    def origin(self) -> tuple:
        """
        Get the origin value (x0, y0, z0).
        """
    def otsu_threshold(
        self, min_val: typing.SupportsInt = 0, max_val: typing.SupportsInt = 256
    ) -> typing.Annotated[list[float], "FixedSize(5)"]: ...
    def paint(self, shape: shape) -> None:
        """
        Paint (set values of) a shape into the image.
        """
    def paintAdd(self, shape: shape) -> None:
        """
        Add (+) a shape's value to the image.
        """
    def paintAddAfter(self, shape: shape) -> None:
        """
        Add (+) a shape's value after the shape (plane...)
        """
    def paintAddBefore(self, shape: shape) -> None:
        """
        Add (+) a shape's value before the shape (plane...)
        """
    def paintAfter(self, shape: shape) -> None:
        """
        Paint after the shape (plane...)
        """
    def paintBefore(self, shape: shape) -> None:
        """
        Paint before the shape (plane...)
        """
    def plotAll(
        self,
        filename: str = "pltAll",
        min_val: typing.SupportsInt = 0,
        max_val: typing.SupportsInt = -1000001,
        slice_index: typing.SupportsInt = -1000000,
        histogram_bins: typing.SupportsInt = 128,
        normal_axis: str = "xyz",
        grey: bool = True,
        color: bool = True,
        histogram: bool = True,
        z_profile: bool = True,
        alpha_image: VxlImgI32 = None,
        alpha_min: typing.SupportsInt = 0,
        alpha_max: typing.SupportsInt = -1000001,
    ) -> bool:
        """
        Plot all visualizations (Histogram, ZProfile, Slices) with various options, hackish for debugging
        """
    def plotHistogramSvg(
        self,
        filename: str = "aa.svg",
        bins: typing.SupportsInt = 128,
        min_val: typing.SupportsFloat = 3e38,
        max_val: typing.SupportsFloat = -3e38,
    ) -> int: ...
    def plotSlice(
        self,
        filename: str,
        normal_axis: str = "xyz",
        slice_index: typing.SupportsInt = -1000000,
        min_val: typing.SupportsInt = 0,
        max_val: typing.SupportsInt = -1000001,
        color_map: str = "gray",
    ) -> None:
        """
        Save a 2D slice as a PNG image.
        normalAxis can be 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'.
        color_map can be 'gray' or 'RGB'
        """
    def plotZProfileSvg(
        self, filename: str = "aa.svg", min_val: typing.SupportsFloat = 0, max_val: typing.SupportsFloat = 255
    ) -> int: ...
    def pointMedian032(
        self, nAdj0: typing.SupportsInt, nAdj1: typing.SupportsInt, lbl0: typing.SupportsInt, lbl1: typing.SupportsInt
    ) -> None:
        """
        Set voxel value to lbl0/1 if it has more than nAdj0/1 neighbours with value lbl0/1, in its 6+26 nearest voxels
        """
    def printInfo(self) -> None: ...
    def rangeTo(self, min: typing.SupportsInt, max: typing.SupportsInt, val: typing.SupportsInt) -> None:
        """
        Set values in range [min, max] to val.
        """
    def readAscii(self, filename: str) -> bool:
        """
        Read image data from an ASCII file.
        """
    def readBin(self, filename: str, nSkipBytes: typing.SupportsInt = 0) -> None:
        """
        Read image data from  a .raw, .raw/.raw.gz, or reset and read from a .am or .tif file.
        """
    def readFromFloat(
        self, header: str, scale: typing.SupportsFloat = 1.0, shift: typing.SupportsFloat = 0.0
    ) -> bool: ...
    def readFromHeader(self, filename: str) -> None:
        """
        Reset and read image from file.
        """
    def redirect(self, arg0: str) -> None:
        """
        Swap X axis with the given axis (y or z).
        """
    def replaceByImageRange(
        self, min_val: typing.SupportsFloat, max_val: typing.SupportsFloat, image_file: str
    ) -> None: ...
    def replaceOutSideValue(
        self, val_old: typing.SupportsInt = 0, val_new: typing.SupportsInt = 2, hole_size: typing.SupportsInt = 5
    ) -> bool: ...
    def replaceRange(self, min: typing.SupportsInt, max: typing.SupportsInt, val: typing.SupportsInt) -> None:
        """
        Replace values in range [min, max] with val.
        """
    def replaceRangeByImage(
        self, min_val: typing.SupportsFloat, max_val: typing.SupportsFloat, image_file: str
    ) -> None: ...
    def resample(self, factor: typing.SupportsFloat) -> None:
        """
        Resample image by factor f, using averaging (downsampling, f>1) or nearest when upsampling (f<1)
        """
    def resampleMax(self, factor: typing.SupportsFloat) -> VxlImgI32:
        """
        Downsample the image, setting voxel values to maximum of original encompassing voxel values.
        """
    def resampleMean(self, factor: typing.SupportsFloat) -> VxlImgI32:
        """
        Resample the image using mean interpolation.
        """
    def resampleMode(self, factor: typing.SupportsFloat) -> VxlImgI32:
        """
        Downsample the image, setting voxel values to mode of original encompassing voxel values.
        """
    def rescaleValues(self, min: typing.SupportsInt, max: typing.SupportsInt) -> None:
        """
        Rescale image values to [min, max].
        """
    def resliceZ(self, factor: typing.SupportsFloat) -> VxlImgI32:
        """
        Reslice along the Z axis.
        """
    def scaleDx(self, scale: typing.SupportsFloat) -> None:
        """
        Scale the voxel size (dx, dy, dz) and origin by a factor.
        """
    def setOffset(self, offset: image3kit._core.sirun.dbl3) -> None:
        """
        Set the spatial offset (origin).
        """
    def setOrigin(self, origin: tuple) -> None:
        """
        Set the origin value (x0, y0, z0).
        """
    def setVoxelSize(self, voxelSize: tuple) -> None:
        """
        Set the voxel size (dx, dy, dz).
        """
    def shape(self) -> tuple: ...
    def shrink0(self) -> None:
        """
        Shrink pore phase (voxel values of 0).
        """
    def shrinkBox(self, num_layers: typing.SupportsInt) -> None:
        """
        Shrink the image boundaries, decreasing its size by the given num_layers in all directions
        """
    def smooth(
        self,
        iterations: typing.SupportsInt = 1,
        kernel_radius: typing.SupportsInt = 1,
        sigma_val: typing.SupportsFloat = 16.0,
        sharpness: typing.SupportsFloat = 0.1,
    ) -> bool:
        """
        bilateral smoothing filter
        """
    def threshold101(self, min: typing.SupportsInt, max: typing.SupportsInt) -> None:
        """
        Apply a threshold to binarize the image, set voxel-values to convert to 0 in between the min and max thresholds and 1 outside of it
        """
    def variance(self, min_val: typing.SupportsInt = 0, max_val: typing.SupportsInt = 255) -> float:
        """
        Set outer tubing of a circular core-holder image to fill_val
        """
    def voxelSize(self) -> tuple:
        """
        Get the voxel size (dx, dy, dz).
        """
    def write(self, filename: str) -> None:
        """
        Write the image to a file (.mhd, .raw, .ra.gz formats).
        """
    def write8bit(self, filename: str, min: typing.SupportsFloat = 0.0, max: typing.SupportsFloat = -0.5) -> None:
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
    def writeNoHeader(self, filename: str) -> None:
        """
        Write the raw image data without a header.
        """

class VxlImgU16:
    def AND(self, image2: VxlImgU16) -> None:
        """
        Voxel-by-voxel inplace AND operation.
        """
    def FaceMedian06(self, nAdj0: typing.SupportsInt, nAdj1: typing.SupportsInt) -> int:
        """
        Set voxel value to 0/1 if it has more than nAdj0/1 neighbours with value 0/1, in its 6 nearest voxels
        """
    def NOT(self, image2: VxlImgU16) -> None:
        """
        Voxel-by-voxel inplace NOT operation, alias for img = img & !img2.
        """
    def OR(self, image2: VxlImgU16) -> None:
        """
        Voxel-by-voxel inplace OR operation.
        """
    def XOR(self, image2: VxlImgU16) -> None:
        """
        Voxel-by-voxel inplace XOR operation.
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def __init__(self, shape: tuple = (0, 0, 0), value: typing.SupportsInt = 0) -> None:
        """
        Initialize a new image of size tuple (nx, ny, nz) with the fill value.
        """
    @typing.overload
    def __init__(self, filepath: typing.Any, processKeys: bool = True) -> None:
        """
        Read image dimensions/metadata from a (header) file. SUpported file types are .am, .raw
        """
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    def __repr__(self) -> str: ...
    def addSurfNoise(
        self,
        mask1: typing.SupportsInt,
        mask2: typing.SupportsInt,
        threshold: typing.SupportsInt,
        seed: typing.SupportsInt,
    ) -> None:
        """
        Add surface noise.
        """
    def adjustBrightnessWith(self, image_file: str) -> bool: ...
    def adjustSliceBrightness(
        self,
        mask_a: VxlImgU8,
        mask_b: VxlImgU8,
        ref_image: VxlImgU16,
        smooth_iter: typing.SupportsInt = 3,
        smooth_kernel: typing.SupportsInt = 20,
    ) -> bool: ...
    def averageWith(self, arg0: str) -> bool: ...
    def averageWith_mBE(self, arg0: str) -> bool: ...
    def bilateralGauss(
        self,
        iterations: typing.SupportsInt = 1,
        kernel_radius: typing.SupportsInt = 1,
        sigma_val: typing.SupportsFloat = 16.0,
        sharpness: typing.SupportsFloat = 0.1,
        sigma_spatial: typing.SupportsFloat = 2.0,
    ) -> bool: ...
    def bilateralX(
        self,
        iterations: typing.SupportsInt = 1,
        kernel_radius: typing.SupportsInt = 1,
        x_step: typing.SupportsInt = 2,
        sigma_val: typing.SupportsFloat = 16.0,
        sharpness: typing.SupportsFloat = 0.1,
        sigma_spatial: typing.SupportsFloat = 2.0,
    ) -> bool:
        """
        Bilateral filter with Xtra large kernel radius, actual kernel size is: kernel_radius * x_step cubed.
        """
    def circleOut(
        self, x: typing.SupportsInt, y: typing.SupportsInt, r: typing.SupportsInt, d: str, val: typing.SupportsInt
    ) -> None:
        """
        Circle out operation.
        """
    def copy(self) -> VxlImgU16:
        """
        duplicate the image data
        """
    def cropD(
        self,
        begin: tuple,
        end: tuple,
        emptyLayers: typing.SupportsInt = 0,
        emptyLayersValue: typing.SupportsInt = 1,
        verbose: bool = False,
    ) -> None:
        """
        Crop the image (inplace) from begin index tupe ix,iy,iz (inclusive) to and and end index (not inclusive) tuple.
        """
    def cutOutside(
        self,
        axis: str = "z",
        extra_out: typing.SupportsInt = 0,
        threshold: typing.SupportsInt = -1,
        cut_highs: typing.SupportsInt = 0,
        shift_x: typing.SupportsInt = 0,
        shift_y: typing.SupportsInt = 0,
        fill_val: typing.SupportsInt = 0,
    ) -> None: ...
    def data(self) -> numpy.typing.NDArray[numpy.uint16]:
        """
        Get the raw data buffer as a numpy array.
        """
    def delense032(
        self,
        iterations: typing.SupportsInt,
        nAdj0: typing.SupportsInt,
        nAdj1: typing.SupportsInt,
        lbl0: typing.SupportsInt,
        lbl1: typing.SupportsInt,
    ) -> None:
        """
        Delense operation.
        """
    def dering(
        self,
        x0: typing.SupportsInt,
        y0: typing.SupportsInt,
        x1: typing.SupportsInt,
        y1: typing.SupportsInt,
        min_val: typing.SupportsInt = 0,
        max_val: typing.SupportsInt = 255,
        nr: typing.SupportsInt = 0,
        ntheta: typing.SupportsInt = 18,
        nz: typing.SupportsInt = 0,
        scale_dif_val: typing.SupportsFloat = 1,
        bilateral_sharpen: typing.SupportsFloat = 0.05,
        nGrowBox: typing.SupportsInt = 10,
        write_dumps: bool = True,
    ) -> None: ...
    def direction(self, arg0: str) -> None:
        """
        Get direction?
        """
    def faceMedNgrowToFrom(
        self, labelTo: typing.SupportsInt, labelFrom: typing.SupportsInt, nDiff: typing.SupportsInt
    ) -> None:
        """
        Face median grow to/from labels.
        """
    def fillHoles(self, maxHoleRadius: typing.SupportsInt) -> None:
        """
        Fill closed holes in the image.
        """
    def flipEndian(self) -> None: ...
    def grow0(self) -> None:
        """
        Grow pore phase (voxel values of 0).
        """
    def growBox(self, num_layers: typing.SupportsInt) -> None:
        """
        Expand the image boundaries, increasing its size by `num_layers` in all directions
        """
    def growLabel(self, arg0: typing.SupportsInt) -> None: ...
    def growingThreshold(
        self,
        startMin: typing.SupportsInt,
        startMax: typing.SupportsInt,
        finalMin: typing.SupportsInt,
        finalMax: typing.SupportsInt,
        iterations: typing.SupportsInt = 4,
    ) -> None: ...
    def keepLargest(self, min: typing.SupportsInt, max: typing.SupportsInt) -> None:
        """
        Keep largest singly-connected region with values in [min, max].
        """
    def mapFrom(
        self,
        sourceImage: VxlImgU16,
        vmin: typing.SupportsInt,
        vmax: typing.SupportsInt,
        scale: typing.SupportsFloat,
        shift: typing.SupportsFloat,
    ) -> None:
        """
        Map values from another image.
        """
    def meanWide(
        self,
        width: typing.SupportsInt = 0,
        noise_val: typing.SupportsInt = 4,
        average: typing.SupportsInt = 0,
        delta: typing.SupportsInt = 20,
        iterations: typing.SupportsInt = 15,
        smooth_image: str = "",
    ) -> bool:
        """
        computes a background image, used to correct for lens artifacts
        """
    def medianFilter(self) -> None:
        """
        Apply a 1+6-neighbour median filter.
        """
    def medianX(self) -> None:
        """
        Apply median filter with kernel size of 1 voxels in x-direction
        """
    def medianY(self) -> None:
        """
        Apply median filter with kernel size of 1 voxels in y-direction
        """
    def medianZ(self) -> None:
        """
        Apply median filter with kernel size of 1 voxels in z-direction
        """
    def mode26(self, arg0: typing.SupportsInt) -> None: ...
    def modeNSames(self, nSameNeighbors: typing.SupportsInt) -> int:
        """
        Apply mode filter based on nearest 6 neighbor voxels.
        """
    def nx(self) -> int: ...
    def ny(self) -> int: ...
    def nz(self) -> int: ...
    def origin(self) -> tuple:
        """
        Get the origin value (x0, y0, z0).
        """
    def otsu_threshold(
        self, min_val: typing.SupportsInt = 0, max_val: typing.SupportsInt = 256
    ) -> typing.Annotated[list[float], "FixedSize(5)"]: ...
    def paint(self, shape: shape) -> None:
        """
        Paint (set values of) a shape into the image.
        """
    def paintAdd(self, shape: shape) -> None:
        """
        Add (+) a shape's value to the image.
        """
    def paintAddAfter(self, shape: shape) -> None:
        """
        Add (+) a shape's value after the shape (plane...)
        """
    def paintAddBefore(self, shape: shape) -> None:
        """
        Add (+) a shape's value before the shape (plane...)
        """
    def paintAfter(self, shape: shape) -> None:
        """
        Paint after the shape (plane...)
        """
    def paintBefore(self, shape: shape) -> None:
        """
        Paint before the shape (plane...)
        """
    def plotAll(
        self,
        filename: str = "pltAll",
        min_val: typing.SupportsInt = 0,
        max_val: typing.SupportsInt = -1000001,
        slice_index: typing.SupportsInt = -1000000,
        histogram_bins: typing.SupportsInt = 128,
        normal_axis: str = "xyz",
        grey: bool = True,
        color: bool = True,
        histogram: bool = True,
        z_profile: bool = True,
        alpha_image: VxlImgU16 = None,
        alpha_min: typing.SupportsInt = 0,
        alpha_max: typing.SupportsInt = -1000001,
    ) -> bool:
        """
        Plot all visualizations (Histogram, ZProfile, Slices) with various options, hackish for debugging
        """
    def plotHistogramSvg(
        self,
        filename: str = "aa.svg",
        bins: typing.SupportsInt = 128,
        min_val: typing.SupportsFloat = 3e38,
        max_val: typing.SupportsFloat = -3e38,
    ) -> int: ...
    def plotSlice(
        self,
        filename: str,
        normal_axis: str = "xyz",
        slice_index: typing.SupportsInt = -1000000,
        min_val: typing.SupportsInt = 0,
        max_val: typing.SupportsInt = -1000001,
        color_map: str = "gray",
    ) -> None:
        """
        Save a 2D slice as a PNG image.
        normalAxis can be 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'.
        color_map can be 'gray' or 'RGB'
        """
    def plotZProfileSvg(
        self, filename: str = "aa.svg", min_val: typing.SupportsFloat = 0, max_val: typing.SupportsFloat = 255
    ) -> int: ...
    def pointMedian032(
        self, nAdj0: typing.SupportsInt, nAdj1: typing.SupportsInt, lbl0: typing.SupportsInt, lbl1: typing.SupportsInt
    ) -> None:
        """
        Set voxel value to lbl0/1 if it has more than nAdj0/1 neighbours with value lbl0/1, in its 6+26 nearest voxels
        """
    def printInfo(self) -> None: ...
    def rangeTo(self, min: typing.SupportsInt, max: typing.SupportsInt, val: typing.SupportsInt) -> None:
        """
        Set values in range [min, max] to val.
        """
    def readAscii(self, filename: str) -> bool:
        """
        Read image data from an ASCII file.
        """
    def readBin(self, filename: str, nSkipBytes: typing.SupportsInt = 0) -> None:
        """
        Read image data from  a .raw, .raw/.raw.gz, or reset and read from a .am or .tif file.
        """
    def readFromFloat(
        self, header: str, scale: typing.SupportsFloat = 1.0, shift: typing.SupportsFloat = 0.0
    ) -> bool: ...
    def readFromHeader(self, filename: str) -> None:
        """
        Reset and read image from file.
        """
    def redirect(self, arg0: str) -> None:
        """
        Swap X axis with the given axis (y or z).
        """
    def replaceByImageRange(
        self, min_val: typing.SupportsFloat, max_val: typing.SupportsFloat, image_file: str
    ) -> None: ...
    def replaceOutSideValue(
        self, val_old: typing.SupportsInt = 0, val_new: typing.SupportsInt = 2, hole_size: typing.SupportsInt = 5
    ) -> bool: ...
    def replaceRange(self, min: typing.SupportsInt, max: typing.SupportsInt, val: typing.SupportsInt) -> None:
        """
        Replace values in range [min, max] with val.
        """
    def replaceRangeByImage(
        self, min_val: typing.SupportsFloat, max_val: typing.SupportsFloat, image_file: str
    ) -> None: ...
    def resample(self, factor: typing.SupportsFloat) -> None:
        """
        Resample image by factor f, using averaging (downsampling, f>1) or nearest when upsampling (f<1)
        """
    def resampleMax(self, factor: typing.SupportsFloat) -> VxlImgU16:
        """
        Downsample the image, setting voxel values to maximum of original encompassing voxel values.
        """
    def resampleMean(self, factor: typing.SupportsFloat) -> VxlImgU16:
        """
        Resample the image using mean interpolation.
        """
    def resampleMode(self, factor: typing.SupportsFloat) -> VxlImgU16:
        """
        Downsample the image, setting voxel values to mode of original encompassing voxel values.
        """
    def rescaleValues(self, min: typing.SupportsInt, max: typing.SupportsInt) -> None:
        """
        Rescale image values to [min, max].
        """
    def resliceZ(self, factor: typing.SupportsFloat) -> VxlImgU16:
        """
        Reslice along the Z axis.
        """
    def scaleDx(self, scale: typing.SupportsFloat) -> None:
        """
        Scale the voxel size (dx, dy, dz) and origin by a factor.
        """
    def segment(
        self,
        n_segments: typing.SupportsInt = 2,
        thresholds: collections.abc.Sequence[typing.SupportsInt] = [],
        min_sizes: collections.abc.Sequence[typing.SupportsInt] = [],
        smooth_image: str = "",
        noise_val: typing.SupportsFloat = 16.0,
        resolution_sq: typing.SupportsFloat = 2.0,
        write_dumps: typing.SupportsInt = 0,
    ) -> bool: ...
    def segment2(
        self,
        thresholds: collections.abc.Sequence[typing.SupportsInt] = [],
        min_sizes: collections.abc.Sequence[typing.SupportsInt] = [],
        noise_val: typing.SupportsFloat = 2.0,
        local_factor: typing.SupportsFloat = 0.05,
        flatnes: typing.SupportsFloat = 0.1,
        effective_resolution: typing.SupportsFloat = 2.0,
        gradient_factor: typing.SupportsFloat = 0.0,
        kernel_radius: typing.SupportsInt = 2,
        n_iterations: typing.SupportsInt = 13,
        write_dumps: typing.SupportsInt = 0,
    ) -> bool: ...
    def setOffset(self, offset: image3kit._core.sirun.dbl3) -> None:
        """
        Set the spatial offset (origin).
        """
    def setOrigin(self, origin: tuple) -> None:
        """
        Set the origin value (x0, y0, z0).
        """
    def setVoxelSize(self, voxelSize: tuple) -> None:
        """
        Set the voxel size (dx, dy, dz).
        """
    def shape(self) -> tuple: ...
    def shrink0(self) -> None:
        """
        Shrink pore phase (voxel values of 0).
        """
    def shrinkBox(self, num_layers: typing.SupportsInt) -> None:
        """
        Shrink the image boundaries, decreasing its size by the given num_layers in all directions
        """
    def smooth(
        self,
        iterations: typing.SupportsInt = 1,
        kernel_radius: typing.SupportsInt = 1,
        sigma_val: typing.SupportsFloat = 16.0,
        sharpness: typing.SupportsFloat = 0.1,
    ) -> bool:
        """
        bilateral smoothing filter
        """
    def threshold101(self, min: typing.SupportsInt, max: typing.SupportsInt) -> None:
        """
        Apply a threshold to binarize the image, set voxel-values to convert to 0 in between the min and max thresholds and 1 outside of it
        """
    def variance(self, min_val: typing.SupportsInt = 0, max_val: typing.SupportsInt = 255) -> float:
        """
        Set outer tubing of a circular core-holder image to fill_val
        """
    def voxelSize(self) -> tuple:
        """
        Get the voxel size (dx, dy, dz).
        """
    def write(self, filename: str) -> None:
        """
        Write the image to a file (.mhd, .raw, .ra.gz formats).
        """
    def write8bit(self, filename: str, min: typing.SupportsFloat = 0.0, max: typing.SupportsFloat = -0.5) -> None:
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
    def writeNoHeader(self, filename: str) -> None:
        """
        Write the raw image data without a header.
        """

class VxlImgU8:
    def AND(self, image2: VxlImgU8) -> None:
        """
        Voxel-by-voxel inplace AND operation.
        """
    def FaceMedian06(self, nAdj0: typing.SupportsInt, nAdj1: typing.SupportsInt) -> int:
        """
        Set voxel value to 0/1 if it has more than nAdj0/1 neighbours with value 0/1, in its 6 nearest voxels
        """
    def NOT(self, image2: VxlImgU8) -> None:
        """
        Voxel-by-voxel inplace NOT operation, alias for img = img & !img2.
        """
    def OR(self, image2: VxlImgU8) -> None:
        """
        Voxel-by-voxel inplace OR operation.
        """
    def XOR(self, image2: VxlImgU8) -> None:
        """
        Voxel-by-voxel inplace XOR operation.
        """
    def __buffer__(self, flags):
        """
        Return a buffer object that exposes the underlying memory of the object.
        """
    @typing.overload
    def __init__(self, shape: tuple = (0, 0, 0), value: typing.SupportsInt = 0) -> None:
        """
        Initialize a new image of size tuple (nx, ny, nz) with the fill value.
        """
    @typing.overload
    def __init__(self, filepath: typing.Any, processKeys: bool = True) -> None:
        """
        Read image dimensions/metadata from a (header) file. SUpported file types are .am, .raw
        """
    def __release_buffer__(self, buffer):
        """
        Release the buffer object that exposes the underlying memory of the object.
        """
    def __repr__(self) -> str: ...
    def addSurfNoise(
        self,
        mask1: typing.SupportsInt,
        mask2: typing.SupportsInt,
        threshold: typing.SupportsInt,
        seed: typing.SupportsInt,
    ) -> None:
        """
        Add surface noise.
        """
    def adjustBrightnessWith(self, image_file: str) -> bool: ...
    def adjustSliceBrightness(
        self,
        mask_a: VxlImgU8,
        mask_b: VxlImgU8,
        ref_image: VxlImgU8,
        smooth_iter: typing.SupportsInt = 3,
        smooth_kernel: typing.SupportsInt = 20,
    ) -> bool: ...
    def averageWith(self, arg0: str) -> bool: ...
    def averageWith_mBE(self, arg0: str) -> bool: ...
    def bilateralGauss(
        self,
        iterations: typing.SupportsInt = 1,
        kernel_radius: typing.SupportsInt = 1,
        sigma_val: typing.SupportsFloat = 16.0,
        sharpness: typing.SupportsFloat = 0.1,
        sigma_spatial: typing.SupportsFloat = 2.0,
    ) -> bool: ...
    def bilateralX(
        self,
        iterations: typing.SupportsInt = 1,
        kernel_radius: typing.SupportsInt = 1,
        x_step: typing.SupportsInt = 2,
        sigma_val: typing.SupportsFloat = 16.0,
        sharpness: typing.SupportsFloat = 0.1,
        sigma_spatial: typing.SupportsFloat = 2.0,
    ) -> bool:
        """
        Bilateral filter with Xtra large kernel radius, actual kernel size is: kernel_radius * x_step cubed.
        """
    def circleOut(
        self, x: typing.SupportsInt, y: typing.SupportsInt, r: typing.SupportsInt, d: str, val: typing.SupportsInt
    ) -> None:
        """
        Circle out operation.
        """
    def copy(self) -> VxlImgU8:
        """
        duplicate the image data
        """
    def cropD(
        self,
        begin: tuple,
        end: tuple,
        emptyLayers: typing.SupportsInt = 0,
        emptyLayersValue: typing.SupportsInt = 1,
        verbose: bool = False,
    ) -> None:
        """
        Crop the image (inplace) from begin index tupe ix,iy,iz (inclusive) to and and end index (not inclusive) tuple.
        """
    def cutOutside(
        self,
        axis: str = "z",
        extra_out: typing.SupportsInt = 0,
        threshold: typing.SupportsInt = -1,
        cut_highs: typing.SupportsInt = 0,
        shift_x: typing.SupportsInt = 0,
        shift_y: typing.SupportsInt = 0,
        fill_val: typing.SupportsInt = 0,
    ) -> None: ...
    def data(self) -> numpy.typing.NDArray[numpy.uint8]:
        """
        Get the raw data buffer as a numpy array.
        """
    def delense032(
        self,
        iterations: typing.SupportsInt,
        nAdj0: typing.SupportsInt,
        nAdj1: typing.SupportsInt,
        lbl0: typing.SupportsInt,
        lbl1: typing.SupportsInt,
    ) -> None:
        """
        Delense operation.
        """
    def dering(
        self,
        x0: typing.SupportsInt,
        y0: typing.SupportsInt,
        x1: typing.SupportsInt,
        y1: typing.SupportsInt,
        min_val: typing.SupportsInt = 0,
        max_val: typing.SupportsInt = 255,
        nr: typing.SupportsInt = 0,
        ntheta: typing.SupportsInt = 18,
        nz: typing.SupportsInt = 0,
        scale_dif_val: typing.SupportsFloat = 1,
        bilateral_sharpen: typing.SupportsFloat = 0.05,
        nGrowBox: typing.SupportsInt = 10,
        write_dumps: bool = True,
    ) -> None: ...
    def direction(self, arg0: str) -> None:
        """
        Get direction?
        """
    def distMapExtrude(
        self,
        distMapDict: dict = {},
        offset: typing.SupportsFloat = 0.5,
        scale: typing.SupportsFloat = 1.0,
        power: typing.SupportsFloat = 1.0,
    ) -> None:
        """
        Extrude proportional to distance map
        """
    def faceMedNgrowToFrom(
        self, labelTo: typing.SupportsInt, labelFrom: typing.SupportsInt, nDiff: typing.SupportsInt
    ) -> None:
        """
        Face median grow to/from labels.
        """
    def fillHoles(self, maxHoleRadius: typing.SupportsInt) -> None:
        """
        Fill closed holes in the image.
        """
    def flipEndian(self) -> None: ...
    def grow0(self) -> None:
        """
        Grow pore phase (voxel values of 0).
        """
    def growBox(self, num_layers: typing.SupportsInt) -> None:
        """
        Expand the image boundaries, increasing its size by `num_layers` in all directions
        """
    def growLabel(self, arg0: typing.SupportsInt) -> None: ...
    def growingThreshold(
        self,
        startMin: typing.SupportsInt,
        startMax: typing.SupportsInt,
        finalMin: typing.SupportsInt,
        finalMax: typing.SupportsInt,
        iterations: typing.SupportsInt = 4,
    ) -> None: ...
    def keepLargest(self, min: typing.SupportsInt, max: typing.SupportsInt) -> None:
        """
        Keep largest singly-connected region with values in [min, max].
        """
    def mapFrom(
        self,
        sourceImage: VxlImgU8,
        vmin: typing.SupportsInt,
        vmax: typing.SupportsInt,
        scale: typing.SupportsFloat,
        shift: typing.SupportsFloat,
    ) -> None:
        """
        Map values from another image.
        """
    def meanWide(
        self,
        width: typing.SupportsInt = 0,
        noise_val: typing.SupportsInt = 4,
        average: typing.SupportsInt = 0,
        delta: typing.SupportsInt = 20,
        iterations: typing.SupportsInt = 15,
        smooth_image: str = "",
    ) -> bool:
        """
        computes a background image, used to correct for lens artifacts
        """
    def medianFilter(self) -> None:
        """
        Apply a 1+6-neighbour median filter.
        """
    def medianX(self) -> None:
        """
        Apply median filter with kernel size of 1 voxels in x-direction
        """
    def medianY(self) -> None:
        """
        Apply median filter with kernel size of 1 voxels in y-direction
        """
    def medianZ(self) -> None:
        """
        Apply median filter with kernel size of 1 voxels in z-direction
        """
    def mode26(self, arg0: typing.SupportsInt) -> None: ...
    def modeNSames(self, nSameNeighbors: typing.SupportsInt) -> int:
        """
        Apply mode filter based on nearest 6 neighbor voxels.
        """
    def nx(self) -> int: ...
    def ny(self) -> int: ...
    def nz(self) -> int: ...
    def origin(self) -> tuple:
        """
        Get the origin value (x0, y0, z0).
        """
    def otsu_threshold(
        self, min_val: typing.SupportsInt = 0, max_val: typing.SupportsInt = 256
    ) -> typing.Annotated[list[float], "FixedSize(5)"]: ...
    def paint(self, shape: shape) -> None:
        """
        Paint (set values of) a shape into the image.
        """
    def paintAdd(self, shape: shape) -> None:
        """
        Add (+) a shape's value to the image.
        """
    def paintAddAfter(self, shape: shape) -> None:
        """
        Add (+) a shape's value after the shape (plane...)
        """
    def paintAddBefore(self, shape: shape) -> None:
        """
        Add (+) a shape's value before the shape (plane...)
        """
    def paintAfter(self, shape: shape) -> None:
        """
        Paint after the shape (plane...)
        """
    def paintBefore(self, shape: shape) -> None:
        """
        Paint before the shape (plane...)
        """
    def plotAll(
        self,
        filename: str = "pltAll",
        min_val: typing.SupportsInt = 0,
        max_val: typing.SupportsInt = -1000001,
        slice_index: typing.SupportsInt = -1000000,
        histogram_bins: typing.SupportsInt = 128,
        normal_axis: str = "xyz",
        grey: bool = True,
        color: bool = True,
        histogram: bool = True,
        z_profile: bool = True,
        alpha_image: VxlImgU8 = None,
        alpha_min: typing.SupportsInt = 0,
        alpha_max: typing.SupportsInt = -1000001,
    ) -> bool:
        """
        Plot all visualizations (Histogram, ZProfile, Slices) with various options, hackish for debugging
        """
    def plotHistogramSvg(
        self,
        filename: str = "aa.svg",
        bins: typing.SupportsInt = 128,
        min_val: typing.SupportsFloat = 3e38,
        max_val: typing.SupportsFloat = -3e38,
    ) -> int: ...
    def plotSlice(
        self,
        filename: str,
        normal_axis: str = "xyz",
        slice_index: typing.SupportsInt = -1000000,
        min_val: typing.SupportsInt = 0,
        max_val: typing.SupportsInt = -1000001,
        color_map: str = "gray",
    ) -> None:
        """
        Save a 2D slice as a PNG image.
        normalAxis can be 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'.
        color_map can be 'gray' or 'RGB'
        """
    def plotZProfileSvg(
        self, filename: str = "aa.svg", min_val: typing.SupportsFloat = 0, max_val: typing.SupportsFloat = 255
    ) -> int: ...
    def pointMedian032(
        self, nAdj0: typing.SupportsInt, nAdj1: typing.SupportsInt, lbl0: typing.SupportsInt, lbl1: typing.SupportsInt
    ) -> None:
        """
        Set voxel value to lbl0/1 if it has more than nAdj0/1 neighbours with value lbl0/1, in its 6+26 nearest voxels
        """
    def printInfo(self) -> None: ...
    def rangeTo(self, min: typing.SupportsInt, max: typing.SupportsInt, val: typing.SupportsInt) -> None:
        """
        Set values in range [min, max] to val.
        """
    def readAscii(self, filename: str) -> bool:
        """
        Read image data from an ASCII file.
        """
    def readBin(self, filename: str, nSkipBytes: typing.SupportsInt = 0) -> None:
        """
        Read image data from  a .raw, .raw/.raw.gz, or reset and read from a .am or .tif file.
        """
    def readFromFloat(
        self, header: str, scale: typing.SupportsFloat = 1.0, shift: typing.SupportsFloat = 0.0
    ) -> bool: ...
    def readFromHeader(self, filename: str) -> None:
        """
        Reset and read image from file.
        """
    def redirect(self, arg0: str) -> None:
        """
        Swap X axis with the given axis (y or z).
        """
    def replaceByImageRange(
        self, min_val: typing.SupportsFloat, max_val: typing.SupportsFloat, image_file: str
    ) -> None: ...
    def replaceOutSideValue(
        self, val_old: typing.SupportsInt = 0, val_new: typing.SupportsInt = 2, hole_size: typing.SupportsInt = 5
    ) -> bool: ...
    def replaceRange(self, min: typing.SupportsInt, max: typing.SupportsInt, val: typing.SupportsInt) -> None:
        """
        Replace values in range [min, max] with val.
        """
    def replaceRangeByImage(
        self, min_val: typing.SupportsFloat, max_val: typing.SupportsFloat, image_file: str
    ) -> None: ...
    def resample(self, factor: typing.SupportsFloat) -> None:
        """
        Resample image by factor f, using averaging (downsampling, f>1) or nearest when upsampling (f<1)
        """
    def resampleMax(self, factor: typing.SupportsFloat) -> VxlImgU8:
        """
        Downsample the image, setting voxel values to maximum of original encompassing voxel values.
        """
    def resampleMean(self, factor: typing.SupportsFloat) -> VxlImgU8:
        """
        Resample the image using mean interpolation.
        """
    def resampleMode(self, factor: typing.SupportsFloat) -> VxlImgU8:
        """
        Downsample the image, setting voxel values to mode of original encompassing voxel values.
        """
    def rescaleValues(self, min: typing.SupportsInt, max: typing.SupportsInt) -> None:
        """
        Rescale image values to [min, max].
        """
    def resliceZ(self, factor: typing.SupportsFloat) -> VxlImgU8:
        """
        Reslice along the Z axis.
        """
    def scaleDx(self, scale: typing.SupportsFloat) -> None:
        """
        Scale the voxel size (dx, dy, dz) and origin by a factor.
        """
    def segment(
        self,
        n_segments: typing.SupportsInt = 2,
        thresholds: collections.abc.Sequence[typing.SupportsInt] = [],
        min_sizes: collections.abc.Sequence[typing.SupportsInt] = [],
        smooth_image: str = "",
        noise_val: typing.SupportsFloat = 16.0,
        resolution_sq: typing.SupportsFloat = 2.0,
        write_dumps: typing.SupportsInt = 0,
    ) -> bool: ...
    def segment2(
        self,
        thresholds: collections.abc.Sequence[typing.SupportsInt] = [],
        min_sizes: collections.abc.Sequence[typing.SupportsInt] = [],
        noise_val: typing.SupportsFloat = 2.0,
        local_factor: typing.SupportsFloat = 0.05,
        flatnes: typing.SupportsFloat = 0.1,
        effective_resolution: typing.SupportsFloat = 2.0,
        gradient_factor: typing.SupportsFloat = 0.0,
        kernel_radius: typing.SupportsInt = 2,
        n_iterations: typing.SupportsInt = 13,
        write_dumps: typing.SupportsInt = 0,
    ) -> bool: ...
    def setOffset(self, offset: image3kit._core.sirun.dbl3) -> None:
        """
        Set the spatial offset (origin).
        """
    def setOrigin(self, origin: tuple) -> None:
        """
        Set the origin value (x0, y0, z0).
        """
    def setVoxelSize(self, voxelSize: tuple) -> None:
        """
        Set the voxel size (dx, dy, dz).
        """
    def shape(self) -> tuple: ...
    def shrink0(self) -> None:
        """
        Shrink pore phase (voxel values of 0).
        """
    def shrinkBox(self, num_layers: typing.SupportsInt) -> None:
        """
        Shrink the image boundaries, decreasing its size by the given num_layers in all directions
        """
    def smooth(
        self,
        iterations: typing.SupportsInt = 1,
        kernel_radius: typing.SupportsInt = 1,
        sigma_val: typing.SupportsFloat = 16.0,
        sharpness: typing.SupportsFloat = 0.1,
    ) -> bool:
        """
        bilateral smoothing filter
        """
    def threshold101(self, min: typing.SupportsInt, max: typing.SupportsInt) -> None:
        """
        Apply a threshold to binarize the image, set voxel-values to convert to 0 in between the min and max thresholds and 1 outside of it
        """
    def variance(self, min_val: typing.SupportsInt = 0, max_val: typing.SupportsInt = 255) -> float:
        """
        Set outer tubing of a circular core-holder image to fill_val
        """
    def voxelSize(self) -> tuple:
        """
        Get the voxel size (dx, dy, dz).
        """
    def write(self, filename: str) -> None:
        """
        Write the image to a file (.mhd, .raw, .ra.gz formats).
        """
    def write8bit(self, filename: str, min: typing.SupportsFloat = 0.0, max: typing.SupportsFloat = -0.5) -> None:
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
    def writeNoHeader(self, filename: str) -> None:
        """
        Write the raw image data without a header.
        """

class cube(shape):
    def __init__(self, p1: tuple, size: tuple, val: typing.SupportsInt) -> None:
        """
        p1: first point, size: size of cuboid sides, val: paint value
        """

class cylinder(shape):
    def __init__(self, p1: tuple, p2: tuple, r: typing.SupportsFloat, val: typing.SupportsInt) -> None:
        """
        p1: first point on axis, p2: second point on axis, r: radius, val: paint value
        """

class shape:
    pass

class sphere(shape):
    def __init__(self, arg0: tuple, arg1: typing.SupportsFloat, arg2: typing.SupportsInt) -> None: ...

class voxelImageTBase:
    pass

def labelImage(arg0: VxlImgU8, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat) -> VxlImgI32: ...
def readImageBase(filename: typing.Any, processKeys: typing.SupportsInt = 1) -> voxelImageTBase:
    """
    Global helper to read an image from a file.
    """
