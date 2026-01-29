
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybind11/pytypes.h"
#include "voxelImage.h"
#include "voxelImageI.h"
#include "shapeToVoxel.h"
#include "InputFile.h"
#include "VxlStrips.h"
#include "voxelImageProcess.h"
#include "voxelNoise.h"

 

namespace py = pybind11;

using namespace VoxLib;


inline InputFile pyCastInput(py::dict dic) {
    InputFile inp;
    for (const auto &kv : dic) {
        inp.add(kv.first.cast<std::string>(), py::str(kv.second).cast<std::string>());
    }
    return inp;
}


template<typename VxTyp>
void addDodgyFuncsInt(py::class_<voxelImageT<VxTyp>> &m) requires(sizeof(VxTyp)>=3) {}

template<typename VxTyp> 
void addDodgyFuncsInt(py::class_<voxelImageT<VxTyp>> &m) requires(sizeof(VxTyp)<=2) {
    m.def("segment2", [](voxelImageT<VxTyp> &m, int nSegs, std::vector<intOr<VxTyp>> th, std::vector<int> minSizs,
                        double noisev, double localF, double flatnes, double resolution, double gradFactor,
                        int krnl, int nItrs, int writedumps) {
            return segment2(m, nSegs, th, minSizs, noisev, localF, flatnes, resolution, gradFactor, krnl, nItrs, writedumps);
        },
        py::arg("nSegs")=2, py::arg("th")=std::vector<int>(), py::arg("minSizs")=std::vector<int>(),
        py::arg("noisev")=2., py::arg("localF")=800., py::arg("flatnes")=0.1, py::arg("resolution")=2.,
        py::arg("gradFactor")=0., py::arg("krnl")=2, py::arg("nItrs")=13, py::arg("writedumps")=0
    )
    .def("segment", [](voxelImageT<VxTyp> &m, int nSegs, std::vector<int> trshlds, std::vector<int> minSizs, std::string smoot, double noisev, double resolutionSqr, int writedumps) {
         return MCTProcessing::segment(m, nSegs, trshlds, minSizs, smoot, noisev, resolutionSqr, writedumps);
        }, py::arg("n_segments")=2, py::arg("thresholds"), py::arg("min_sizes"), py::arg("smooth_image")="", 
        py::arg("noise_val")=16.0, py::arg("resolution_sq")=2.0, py::arg("write_dumps")=0
    )
    ;
}

template<typename VxTyp>
void addDodgyFuncsU8(py::class_<voxelImageT<VxTyp>> &m) requires(sizeof(VxTyp)>=2) {
}
template<typename VxTyp> 
void addDodgyFuncsU8(py::class_<voxelImageT<VxTyp>> &m) requires(sizeof(VxTyp)<=1) {
    m
    .def("distMapExtrude", [](voxelImageT<VxTyp> &m, py::dict dic, double offsetFactor, double scaleR, double powR) { // returns copy
        m = distMapExtrude(m, pyCastInput(dic), offsetFactor, scaleR, powR, true); },
        py::arg("distMapDict")=py::dict(), py::arg("offset")=0.5, py::arg("scale")=1.0, py::arg("power")=1.0,
        "Extrude proportional to distance map"
    )
    ;
}

template<typename VxTyp> 
void bind_VxlImg(py::module &mod, const char* VxTypS) {

auto clas = py::class_<voxelImageT<VxTyp>>(mod, VxTypS, py::buffer_protocol())
.def_buffer([](voxelImageT<VxTyp> &m) -> py::buffer_info {
        return py::buffer_info(
            m.data(),                             // Pointer to buffer
            sizeof(VxTyp),                          // Size of one scalar
            py::format_descriptor<VxTyp>::format(), // Python struct-style format descriptor
            3,                                    // Number of dimensions
            { size_t(m.nx()), size_t(m.ny()), size_t(m.nz()) },             // Dimensions
            { long(sizeof(VxTyp)), long(sizeof(VxTyp) *m.nx()), long(sizeof(VxTyp) * m.nxy()) } // Strides (in bytes)
        );
    })
    .def(py::init([](py::tuple tpl, VxTyp value) { return voxelImageT<VxTyp>(tpl[0].cast<int>(), tpl[1].cast<int>(), tpl[2].cast<int>(), value); }),
         py::arg("shape")=py::make_tuple(0,0,0), py::arg("value") = 0, "Initialize from a size tuple (nz, ny, nx) with an optional fill value.")
    .def("__repr__", [VxTypS](const voxelImageT<VxTyp> &m) {
        return "<" + std::string(VxTypS) + " shape=(" + 
               _s(m.size3()) + ")>";
    })
//.def(py::init<int, int, int, VxTyp>())
    .def(py::init(
        [](py::object filepath, bool processKeys)  { return voxelImageT<VxTyp>(py::str(filepath).cast<std::string>(), processKeys? readOpt::procAndConvert : readOpt::justRead); }), 
        py::arg("filepath"), py::arg("processKeys")=true, "Read image dimensions/metadata from a (header) file.") // readConvertFromHeader
.def("data", [](voxelImageT<VxTyp> &m) {
    return py::array_t<VxTyp>(py::buffer_info{
        m.data(),
        sizeof(VxTyp),
        py::format_descriptor<VxTyp>::format(),
        3,
        { size_t(m.nx()), size_t(m.ny()), size_t(m.nz()) },             // Dimensions
        { long(sizeof(VxTyp)), long(sizeof(VxTyp) *m.nx()), long(sizeof(VxTyp) * m.nxy()) } // Strides (in bytes)
    }); }, "Get the raw data buffer as a numpy array.")
.def("nx", &voxelImageT<VxTyp>::nx)
.def("ny", &voxelImageT<VxTyp>::ny)
.def("nz", &voxelImageT<VxTyp>::nz)
.def("printInfo", &voxelImageT<VxTyp>::printInfo)
.def("shape", [&](voxelImageT<VxTyp> &m) { return py::make_tuple(m.nx(), m.ny(), m.nz()); })
.def("write", &voxelImageT<VxTyp>::write, py::arg("filename"), "Write the image to a file (.mhd, .raw, .ra.gz formats).")
.def("writeNoHdr", &voxelImageT<VxTyp>::writeNoHdr, py::arg("filename"), "Write the raw image data without a header.")
// .def("writeHeader", &voxelImageT<VxTyp>::writeHeader)
.def("readAscii", &voxelImageT<VxTyp>::readAscii, py::arg("filename"), "Read image data from an ASCII file.")
// .def("readRLE", &voxelImageT<VxTyp>::readRLE)
.def("cropD", &voxelImageT<VxTyp>::cropD, py::arg("begin"), py::arg("end"), py::arg("emptyLayers")=0, py::arg("emptyLayersValue")=1, py::arg("verbose")=false, "Crop the image by a specified depth.")
.def("growBox", &voxelImageT<VxTyp>::growBox, py::arg("layers"), "Expand the image boundaries.")
.def("shrinkBox", &voxelImageT<VxTyp>::shrinkBox, py::arg("layers"), "Shrink the image boundaries.")
.def("shrinkBox", &voxelImageT<VxTyp>::shrinkBox, py::arg("layers"), "Shrink the image boundaries.")
.def("fillHoles", &voxelImageT<VxTyp>::fillHoles, py::arg("maxHoleRadius"), "Fill closed holes in the image.")
.def("threshold101", &voxelImageT<VxTyp>::threshold101, py::arg("min"), py::arg("max"), "Apply a threshold to convert to 0/1.")
.def("writeAConnectedPoreVoxel", &voxelImageT<VxTyp>::writeAConnectedPoreVoxel, py::arg("filename"), "Write a specific connected pore voxel to file.")
.def("AND", &voxelImageT<VxTyp>::AND, py::arg("other"), "Pixel-wise AND operation.")
.def("NOT", &voxelImageT<VxTyp>::NOT, "Pixel-wise NOT operation.")
.def("OR", &voxelImageT<VxTyp>::OR, py::arg("other"), "Pixel-wise OR operation.")
.def("XOR", &voxelImageT<VxTyp>::XOR, py::arg("other"), "Pixel-wise XOR operation.")
.def("resampleMode", [&](voxelImageT<VxTyp> &m, double nReSampleNotSafe) { return resampleMode(m,nReSampleNotSafe); }, py::arg("rate"), "Resample the image using mode interpolation.")
.def("resampleMax", [&](voxelImageT<VxTyp> &m, double nReSampleNotSafe) { return resampleMax(m,nReSampleNotSafe); }, py::arg("rate"), "Resample the image using max interpolation.")
.def("resampleMean", [&](voxelImageT<VxTyp> &m, double nReSampleNotSafe) { return resampleMean(m,nReSampleNotSafe); }, py::arg("rate"), "Resample the image using mean interpolation.")
.def("resliceZ", [&](voxelImageT<VxTyp> &m, double nReSampleNotSafe) { return resliceZ(m,nReSampleNotSafe); }, py::arg("rate"), "Reslice along Z axis.")
.def("Paint"         , [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, setIn); }, py::arg("shape"), "Paint a shape into the image.")
.def("PaintAdd"      , [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, addTo); }, py::arg("shape"), "Add a shape's values to the image.")
.def("PaintBefore"   , [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, setBefor); }, py::arg("shape"), "Paint before...?")
.def("PaintAfter"    , [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, setAfter); }, py::arg("shape"), "Paint after...?")
.def("PaintAddBefore", [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, addBefor); }, py::arg("shape"), "Add paint before...?")
.def("PaintAddAfter" , [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, addAfter); }, py::arg("shape"), "Add paint after...?")
.def("writeContour"  , [&](voxelImageT<VxTyp> &m, const string& outSurf) { // TODO
    InputFile inp;
    if (outSurf.size())  inp.set("outputSurface", outSurf);
    // vxlToSurfWrite(m, inp);
    }, py::arg("outSurf"), "Write contour extraction to a surface file.")
.def("sliceToPng"    , [&](voxelImageT<VxTyp> &m, const string& normalAxis, const string& fnam, int iSlice, int bgnv, int endv, const string& color) {
    ::sliceToPng(m, normalAxis, fnam, iSlice, bgnv, endv, color);
    }, py::arg("normalAxis"), py::arg("filename"), py::arg("sliceIndex"), py::arg("val_min"), py::arg("val_max"), py::arg("color_map")="gray", "Save a 2D slice as a PNG image.")
.def("sliceFromPng"    , [&](voxelImageT<VxTyp> &m, const string& normalAxis, const string& fnam, int iSlice, int bgnv, int endv) {
    ::sliceFromPng(m, normalAxis, fnam, iSlice, bgnv, endv);
    }, py::arg("normalAxis"), py::arg("filename"), py::arg("sliceIndex"), py::arg("val_min"), py::arg("val_max"), "Read a slice from a Png image")
// Member functions and wrappers
.def("rescaleValues", [](voxelImageT<VxTyp> &m, VxTyp min, VxTyp max) { rescaleValues(m, min, max); }, py::arg("min"), py::arg("max"), "Rescale image values to [min, max].")
.def("setOffset", [](voxelImageT<VxTyp> &m, dbl3 off) { m.X0Ch() = off; }, py::arg("offset"), "Set the spatial offset (origin).")
.def("redirect", &voxelImageT<VxTyp>::rotate, "Rotate/Redirect image.")
.def("direction", &voxelImageT<VxTyp>::rotate, "Get direction?")
.def("grow0", &voxelImageT<VxTyp>::growPore, "Grow pore phase (0).")
.def("shrink0", &voxelImageT<VxTyp>::shrinkPore, "Shrink pore phase (0).")
.def("resample", [](voxelImageT<VxTyp> &m, double f) { m = resampleMean(m, f); }, py::arg("factor"), "Resample image by factor f.")
.def("rangeTo", [](voxelImageT<VxTyp> &m, VxTyp min, VxTyp max, VxTyp v) { replaceRange(m, min, max, v); }, py::arg("min"), py::arg("max"), py::arg("val"), "Set values in range [min, max] to val.")
.def("replaceRange", [](voxelImageT<VxTyp> &m, VxTyp min, VxTyp max, VxTyp v) { replaceRange(m, min, max, v); }, py::arg("min"), py::arg("max"), py::arg("val"), "Replace values in range [min, max] with val.")
.def("write8bit", [](voxelImageT<VxTyp> &m, std::string outName, double minv, double maxv) { _write8bit(m, outName, minv, maxv); }, py::arg("filename"), py::arg("min"), py::arg("max"), "Write as 8-bit image scaled between min and max.")
.def("read", [](voxelImageT<VxTyp> &m, std::string s) { m.reset(int3(0), 0); m.readFromHeader(s); }, py::arg("filename"), "Read image from file.")
// //.def("readAtZ", &voxelImageT<VxTyp>::readAtZ)
.def("modeNSames", [](voxelImageT<VxTyp> &m, const short nSameNei) { return modeNSames(m, nSameNei); }, py::arg("nSameNeighbors"), "Apply mode filter based on neighbor count.")
.def("medianFilter", [](voxelImageT<VxTyp> &m) { m = median(m); }, "Apply median filter.")
.def("medianX", [](voxelImageT<VxTyp> &m) { m = medianx(m); }, "Apply median X filter.")
.def("FaceMedian06", &voxelImageT<VxTyp>::FaceMedian06)
.def("PointMedian032", &voxelImageT<VxTyp>::PointMedian032)
.def("faceMedNgrowToFrom", [](voxelImageT<VxTyp> &m, VxTyp lblTo, VxTyp lblFrm, int ndif) { FaceMedGrowToFrom(m, lblTo, lblFrm, ndif); }, py::arg("labelTo"), py::arg("labelFrom"), py::arg("nDiff"), "Face median grow to/from labels.")
.def("delense032", [](voxelImageT<VxTyp> &m, int x, int y, int r, char d, VxTyp v) { _delense032(m, x, y, r, d, v); }, py::arg("x"), py::arg("y"), py::arg("r"), py::arg("d"), py::arg("val"), "Delense operation.")
.def("circleOut", [](voxelImageT<VxTyp> &m, int x, int y, int r, char d, VxTyp v) { circleOut(m, x, y, r, d, v); }, py::arg("x"), py::arg("y"), py::arg("r"), py::arg("d"), py::arg("val"), "Circle out operation.")
.def("growLabel", &voxelImageT<VxTyp>::growLabel)
.def("keepLargest0", [](voxelImageT<VxTyp> &m) { keepLargest0(m); })
//.def("maskWriteFraction", &voxelImageT<VxTyp>::maskWriteFraction)
.def("mapFrom", [](voxelImageT<VxTyp>& m, const voxelImageT<VxTyp>& vimg2, VxTyp vmin, VxTyp vmax, double scale, double shift) { mapToFrom(m, vimg2, vmin, vmax, scale, shift); }, py::arg("sourceImage"), py::arg("vmin"), py::arg("vmax"), py::arg("scale"), py::arg("shift"), "Map values from another image.") 
.def("addSurfNoise", [](voxelImageT<VxTyp> &m, const int randMask1, const int randMask2, int trsh, int randseed) { addSurfNoise(m, randMask1, randMask2, trsh, randseed); }, py::arg("mask1"), py::arg("mask2"), py::arg("threshold"), py::arg("seed"), "Add surface noise.")

// // Individual bindings for convenience
.def("averageWith", [](voxelImageT<VxTyp> &m, std::string args) { std::stringstream ss(args); return MCTProcessing::averageWith(ss, m); })
.def("averageWith_mBE", [](voxelImageT<VxTyp> &m, std::string args) { std::stringstream ss(args); return MCTProcessing::averageWith_mBE(ss, m); })
    .def("plotAll", [](voxelImageT<VxTyp> &m, std::string fnam_, int minv, int maxv, int iSlice_, int nBins,
                        std::string normalAxis, bool grey, bool color, bool histogram, bool zProfile,
                        voxelImageT<VxTyp>* img2Ptr, int mina, int maxa) {
        int colrGreyHistZprofile = grey*1 + color*2 + histogram*4 + zProfile*8;
        return MCTProcessing::plotAll(m, minv, maxv, iSlice_, nBins, normalAxis, fnam_, colrGreyHistZprofile, img2Ptr, mina, maxa);
    },  py::arg("name")="pltAll",
        py::arg("minv")=0, py::arg("maxv")=-1000001, py::arg("sliceIndex")=-1000000, py::arg("nBins")=128,
        py::arg("normalAxis")="xyz",
        py::arg("grey")=true, py::arg("color")=true, py::arg("histogram")=true, py::arg("zProfile")=true,
        py::arg("alphaImage")=nullptr, py::arg("alphaMin")=0, py::arg("alphaMax")=-1000001,
        "Plot all visualizations (Histogram, ZProfile, Slices) with various options, hackish for debugging")
    .def("mode26", [](voxelImageT<VxTyp> &m, int nMinD) { return mode26(m,nMinD); })
    .def("growingThreshold", [](voxelImageT<VxTyp> &m, VxTyp t1, VxTyp t2, VxTyp t3, VxTyp t4, int nIter) {
        ::growingThreshold(m, t1, t2, t3, t4, nIter);
    }, py::arg("startMin"), py::arg("startMax"), py::arg("finalMin"), py::arg("finalMax"), py::arg("iterations")=4)
     .def("replaceOutSideValue", [](voxelImageT<VxTyp> &m, int vo, int vnew, int nHoleSize) {
        return MCTProcessing::replaceOutSideValue(m, vo, vnew, nHoleSize);
    }, py::arg("val_old")=0, py::arg("val_new")=2, py::arg("hole_size")=5)
    .def("smooth", [](voxelImageT<VxTyp> &m, int nItrs, int kernRad, double sigmavv, double sharpFact) {
        return MCTProcessing::smooth(m, nItrs, kernRad, sigmavv, sharpFact);
    }, py::arg("iterations")=1, py::arg("kernel_radius")=1, py::arg("sigma_val")=16.0, py::arg("sharpness")=0.1)
    .def("svgHistogram", [](voxelImageT<VxTyp> &m, std::string fnam, int nBins, double minV, double maxV) {
        return MCTProcessing::svgHistogram(m, fnam, nBins, minV, maxV);
    }, py::arg("filename")="aa.svg", py::arg("bins")=128, py::arg("min_val")=3e38, py::arg("max_val")=-3e38)
    .def("svgZProfile", [](voxelImageT<VxTyp> &m, std::string fnam, double minV, double maxV) {
        return MCTProcessing::svgZProfile(m, fnam, VxTyp(minV), VxTyp(maxV));
    }, py::arg("filename")="aa.svg", py::arg("min_val")=0, py::arg("max_val")=255) // Assuming 255 default for maxV based on typical usage, though maxT(T) is ideal
    .def("flipEndian", [](voxelImageT<VxTyp> &m) {
        ::flipEndian(m);
    })
    .def("replaceRangeByImage", [](voxelImageT<VxTyp> &m, double minv, double maxv, std::string fnam) {
        if(fnam.empty()) throw std::runtime_error("no image name provided");
        ::replaceRangeByImage(m, VxTyp(minv), VxTyp(maxv), voxelImageT<VxTyp>(fnam, readOpt::procAndConvert));
    }, py::arg("min_val"), py::arg("max_val"), py::arg("image_file"))
    .def("replaceByImageRange", [](voxelImageT<VxTyp> &m, double minv, double maxv, std::string fnam) {
        if(fnam.empty()) throw std::runtime_error("no image name provided");
        ::replaceByImageRange(m, VxTyp(minv), VxTyp(maxv), voxelImageT<VxTyp>(fnam, readOpt::procAndConvert));
    }, py::arg("min_val"), py::arg("max_val"), py::arg("image_file"))
    .def("readFromFloat", [](voxelImageT<VxTyp> &m, std::string header, float scale, float shift) {
        return MCTProcessing::readFromFloat(m, header, scale, shift);
    }, py::arg("header"), py::arg("scale")=1.0f, py::arg("shift")=0.0f)
    .def("bilateralX", [](voxelImageT<VxTyp> &m, int nItrs, int kernRad, int Xstp, double sigmavv, double sharpFact, double sigmadd) {
         return ::bilateralX(m, nItrs, kernRad, Xstp, sigmavv, sharpFact, sigmadd);
    }, py::arg("iterations")=1, py::arg("kernel_radius")=1, py::arg("x_step")=2, py::arg("sigma_val")=16.0, py::arg("sharpness")=0.1, py::arg("sigma_spatial")=2.0)
    .def("bilateralGauss", [](voxelImageT<VxTyp> &m, int nItrs, int kernRad, double sigmavv, double sharpFact, double sigmadd) {
         return ::bilateralGauss(m, nItrs, kernRad, sigmavv, sharpFact, sigmadd);
    }, py::arg("iterations")=1, py::arg("kernel_radius")=1, py::arg("sigma_val")=16.0, py::arg("sharpness")=0.1, py::arg("sigma_spatial")=2.0)
    .def("meanWide", [](voxelImageT<VxTyp> &m, int nW, int noisev, int avg, int delta, int nItrs, std::string smoothImg) {
         return MCTProcessing::meanWide(m, nW, noisev, avg, delta, nItrs, smoothImg);
    }, py::arg("width")=0, py::arg("noise_val")=4, py::arg("average")=0, py::arg("delta")=20, py::arg("iterations")=15, py::arg("smooth_image")="")
    .def("otsu_th", [](voxelImageT<VxTyp> &m, int minv, int maxv) {
         return ::otsu_th(m, minv, maxv);
    }, py::arg("min_val")=0, py::arg("max_val")=256)
    .def("dering", [](voxelImageT<VxTyp> &m, int X0, int Y0, int X1, int Y1, int minV, int maxV, int nr, int ntheta, int nz) {
         return MCTProcessing::dering(m, X0, Y0, X1, Y1, minV, maxV, nr, ntheta, nz);
    }, py::arg("x0"), py::arg("y0"), py::arg("x1"), py::arg("y1"), py::arg("min_val")=0, py::arg("max_val")=255, py::arg("nr")=0, py::arg("ntheta")=18, py::arg("nz")=0)
    .def("adjustBrightnessWith", [](voxelImageT<VxTyp> &m, std::string imgName) {
         return MCTProcessing::adjustBrightnessWith(m, imgName);
    }, py::arg("image_file"))
    .def("adjustSliceBrightness", [](voxelImageT<VxTyp> &m, voxelImageT<unsigned char>& mskA, voxelImageT<unsigned char>& mskB, voxelImageT<VxTyp>& img2, int nSmoothItr, int nSmoothKrnl) {
         return MCTProcessing::adjustSliceBrightness(m, mskA, mskB, img2, nSmoothItr, nSmoothKrnl);
    }, py::arg("mask_a"), py::arg("mask_b"), py::arg("ref_image"), py::arg("smooth_iter")=3, py::arg("smooth_kernel")=20)
    .def("cutOutside", [](voxelImageT<VxTyp> &m, char dir, int nExtraOut, int threshold, int cuthighs, int nShiftX, int nShiftY, int outVal) {
         ::cutOutside(m, dir, nExtraOut, threshold, cuthighs, nShiftX, nShiftY, VxTyp(outVal));
         return true;
    }, py::arg("axis")='z', py::arg("extra_out")=0, py::arg("threshold")=-1, py::arg("cut_highs")=0, py::arg("shift_x")=0, py::arg("shift_y")=0, py::arg("fill_val")=0)
    .def("variance", [](voxelImageT<VxTyp> &m, int minV, int maxV) {
         return ::varianceDbl(m, minV, maxV);
    }, py::arg("min_val")=0, py::arg("max_val")=255)
    ;
    addDodgyFuncsInt(clas);
    addDodgyFuncsU8(clas);
}
