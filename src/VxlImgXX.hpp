
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybind11/pytypes.h"
#include "voxelEndian.h"
#include "voxelImage.h"
#include "voxelImageI.h"
#include "shapeToVoxel.h"
#include "InputFile.h"
#include "VxlStrips.h"
#include "voxelImageProcess.h"
#include "voxelNoise.h"

 
namespace VxlPy {
namespace py = pybind11;
using py::arg;
using namespace VoxLib;


inline InputFile pyCastInput(py::dict dic) {
    InputFile inp;
    for (const auto &kv : dic) {
        inp.add(kv.first.cast<std::string>(), py::str(kv.second).cast<std::string>());
    }
    return inp;
}


template<typename VxT>
void addDodgyFuncsInt(py::class_<voxelImageT<VxT>> &m) requires(sizeof(VxT)>=3) {}

template<typename VxT> 
void addDodgyFuncsInt(py::class_<voxelImageT<VxT>> &m) requires(sizeof(VxT)<=2) {

    using SelfT = voxelImageT<VxT>;

    m.def("segment2", [](SelfT &m, int nSegs, std::vector<intOr<VxT>> th, std::vector<int> minSizs,
                        double noisev, double localF, double flatnes, double resolution, double gradFactor,
                        int krnl, int nItrs, int writedumps) {
            return segment2(m, nSegs, th, minSizs, noisev, localF, flatnes, resolution, gradFactor, krnl, nItrs, writedumps);
        },
        arg("nSegs")=2, arg("th")=std::vector<int>(), arg("minSizs")=std::vector<int>(),
        arg("noisev")=2., arg("localF")=800., arg("flatnes")=0.1, arg("resolution")=2.,
        arg("gradFactor")=0., arg("krnl")=2, arg("nItrs")=13, arg("writedumps")=0
    )
    .def("segment", [](SelfT &m, int nSegs, std::vector<int> trshlds, std::vector<int> minSizs, std::string smoot, double noisev, double resolutionSqr, int writedumps) {
         return MCTProcessing::segment(m, nSegs, trshlds, minSizs, smoot, noisev, resolutionSqr, writedumps);
        }, arg("n_segments")=2, arg("thresholds"), arg("min_sizes"), arg("smooth_image")="", 
        arg("noise_val")=16.0, arg("resolution_sq")=2.0, arg("write_dumps")=0
    )
    ;
}

template<typename VxT>
void addDodgyFuncsU8(py::class_<voxelImageT<VxT>> &m) requires(sizeof(VxT)>=2) {
}
template<typename VxT> 
void addDodgyFuncsU8(py::class_<voxelImageT<VxT>> &m) requires(sizeof(VxT)<=1) {
    m
    .def("distMapExtrude", [](voxelImageT<VxT> &m, py::dict dic, double offsetFactor, double scaleR, double powR) { // returns copy
        m = distMapExtrude(m, pyCastInput(dic), offsetFactor, scaleR, powR, true); },
        arg("distMapDict")=py::dict(), arg("offset")=0.5, arg("scale")=1.0, arg("power")=1.0,
        "Extrude proportional to distance map"
    )
    ;
}

template<typename VxT> 
void bind_VxlImg(py::module &mod, const char* VxTypS) {

using SelfT = voxelImageT<VxT>;

auto clas = py::class_<SelfT>(mod, VxTypS, py::buffer_protocol())
.def_buffer([](SelfT &m) -> py::buffer_info {
        return py::buffer_info(
            m.data(), sizeof(VxT),  py::format_descriptor<VxT>::format(),
            3, { size_t(m.nx()), size_t(m.ny()), size_t(m.nz()) },             // Dimensions
            { long(sizeof(VxT)), long(sizeof(VxT)*m.nx()), long(sizeof(VxT) * m.nxy()) } // Strides (in bytes)
        );
    })
    .def(py::init([](py::tuple tpl, VxT value) { return SelfT(tpl[0].cast<int>(), tpl[1].cast<int>(), tpl[2].cast<int>(), value); }),
         arg("shape")=py::make_tuple(0,0,0), arg("value") = 0, "Initialize a new image of size tuple (nx, ny, nz) with the fill value.")
    .def("__repr__", [VxTypS](const SelfT &m) {
        return "<" + std::string(VxTypS) + " shape=("+_s(m.size3())+")>";
    })
    .def(py::init( // readConvertFromHeader
        [](py::object filepath, bool processKeys)  { return SelfT(py::str(filepath).cast<std::string>(), processKeys? readOpt::procAndConvert : readOpt::justRead); }), 
        arg("filepath"), arg("processKeys")=true, 
        "Read image dimensions/metadata from a (header) file. SUpported file types are .am, .raw")
    .def("data", [](SelfT &m) {
        return py::array_t<VxT>(py::buffer_info{
            m.data(), sizeof(VxT), py::format_descriptor<VxT>::format(), 3,
            { size_t(m.nx()), size_t(m.ny()), size_t(m.nz()) },             // Dimensions
            { long(sizeof(VxT)), long(sizeof(VxT) *m.nx()), long(sizeof(VxT) * m.nxy()) } // Strides (in bytes)
        }); }, "Get the raw data buffer as a numpy array.")
    .def("nx", &SelfT::nx)
    .def("ny", &SelfT::ny)
    .def("nz", &SelfT::nz)
    .def("printInfo", &SelfT::printInfo)
    .def("shape", [&](SelfT &m) { return py::make_tuple(m.nx(), m.ny(), m.nz()); })
    .def("write", &SelfT::write, arg("filename"), "Write the image to a file (.mhd, .raw, .ra.gz formats).")
    .def("writeNoHeader", &SelfT::writeNoHdr, arg("filename"), "Write the raw image data without a header.")
    // .def("writeHeader", &SelfT::writeHeader)
    .def("readAscii", &SelfT::readAscii, arg("filename"), "Read image data from an ASCII file.")
    // .def("readRLE", &SelfT::readRLE)
    .def("cropD", &SelfT::cropD, arg("begin"), arg("end"), arg("emptyLayers")=0, arg("emptyLayersValue")=1, arg("verbose")=false, 
         "Crop the image (inplace) from begin index tupe ix,iy,iz (inclusive) to and end index tuple.")
    .def("growBox", &SelfT::growBox, arg("num_layers"), 
         "Expand the image boundaries, increasing its size by `num_layers` in all directions")
    .def("shrinkBox", &SelfT::shrinkBox, arg("num_layers"), 
         "Shrink the image boundaries, decreasing its size by the given num_layers in all directions")
    .def("fillHoles", &SelfT::fillHoles, arg("maxHoleRadius"), "Fill closed holes in the image.")
    .def("threshold101", &SelfT::threshold101, arg("min"), arg("max"),
         "Apply a threshold to binarize the image, set voxel-values to convert to 0 in between the min and max thresholds and 1 outside of it")
    .def("writeAConnectedPoreVoxel", &SelfT::writeAConnectedPoreVoxel, arg("filename"), "Write a specific connected pore voxel to file.")
    .def("AND", &SelfT::AND, arg("other"), "Voxel-by-voxel AND operation.")
    .def("NOT", &SelfT::NOT, "Voxel-by-voxel NOT operation.")
    .def("OR", &SelfT::OR, arg("other"), "Voxel-by-voxel OR operation.")
    .def("XOR", &SelfT::XOR, arg("other"), "Voxel-by-voxel XOR operation.")
    .def("resampleMode", [&](SelfT &m, double nReSampleNotSafe) { return resampleMode(m,nReSampleNotSafe); }, arg("factor"), 
        "Downsample the image, setting voxel values to mode of original encompassing voxel values.")
    .def("resampleMax", [&](SelfT &m, double nReSampleNotSafe) { return resampleMax(m,nReSampleNotSafe); }, arg("factor"), 
        "Downsample the image, setting voxel values to maximum of original encompassing voxel values.")
    .def("resampleMean", [&](SelfT &m, double nReSampleNotSafe) { return resampleMean(m,nReSampleNotSafe); }, arg("factor"), "Resample the image using mean interpolation.")
    .def("resliceZ", [&](SelfT &m, double nReSampleNotSafe) { return resliceZ(m,nReSampleNotSafe); }, arg("factor"), "Reslice along the Z axis.")
    .def("Paint"         , [&](SelfT &m, const shape& sh) { _SHAPERATEPy(m, sh, setIn);    }, arg("shape"), "Paint (set values of) a shape into the image.")
    .def("PaintAdd"      , [&](SelfT &m, const shape& sh) { _SHAPERATEPy(m, sh, addTo);    }, arg("shape"), "Add (+) a shape's value to the image.")
    .def("PaintBefore"   , [&](SelfT &m, const shape& sh) { _SHAPERATEPy(m, sh, setBefor); }, arg("shape"), "Paint before the shape (plane...)")
    .def("PaintAfter"    , [&](SelfT &m, const shape& sh) { _SHAPERATEPy(m, sh, setAfter); }, arg("shape"), "Paint after the shape (plane...)")
    .def("PaintAddBefore", [&](SelfT &m, const shape& sh) { _SHAPERATEPy(m, sh, addBefor); }, arg("shape"), "Add (+) a shape's value before the shape (plane...)")
    .def("PaintAddAfter" , [&](SelfT &m, const shape& sh) { _SHAPERATEPy(m, sh, addAfter); }, arg("shape"), "Add (+) a shape's value after the shape (plane...)")
    .def("writeContour"  , [&](SelfT &m, const string& outSurf) { // TODO
        InputFile inp;
        if (outSurf.size())  inp.set("outputSurface", outSurf);
        // vxlToSurfWrite(m, inp);
        }, arg("outSurf"), "Write contour extraction to a surface file.")
    .def("sliceToPng"    , [&](SelfT &m, const string& normalAxis, const string& fnam, int iSlice, int bgnv, int endv, const string& color) {
        ::sliceToPng(m, normalAxis, fnam, iSlice, bgnv, endv, color);
        }, arg("normalAxis"), arg("filename"), arg("sliceIndex")=-1000000, arg("val_min")=0, arg("val_max")=-1000001, arg("color_map")="gray",
        "Save a 2D slice as a PNG image, color_map can be 'gray' or 'RGB'")
    // .def("sliceFromPng"    , [&](SelfT &m, const string& normalAxis, const string& fnam, int iSlice, int bgnv, int endv) {
    //     ::sliceFromPng(m, normalAxis, fnam, iSlice, bgnv, endv);
    //     }, arg("normalAxis"), arg("filename"), arg("sliceIndex")=0, arg("val_min")=0, arg("val_max")=255, "Read a slice from a Png image")
    // Member functions and wrappers
    .def("rescaleValues", [](SelfT &m, VxT min, VxT max) { rescaleValues(m, min, max); }, arg("min"), arg("max"), "Rescale image values to [min, max].")
    .def("setOffset", [](SelfT &m, dbl3 off) { m.X0Ch() = off; }, arg("offset"), "Set the spatial offset (origin).")
    .def("redirect", &SelfT::rotate, "Swap X axis with the given axis (y or z).")
    .def("direction", &SelfT::rotate, "Get direction?")
    .def("grow0", &SelfT::growPore, "Grow pore phase (voxel values of 0).")
    .def("shrink0", &SelfT::shrinkPore, "Shrink pore phase (voxel values of 0).")
    .def("resample", [](SelfT &m, double f) { m = resampleMean(m, f); }, arg("factor"), 
         "Resample image by factor f, using averaging (downsampling, f>1) or nearest when upsampling (f<1)")
    .def("rangeTo", [](SelfT &m, VxT min, VxT max, VxT val) { replaceRange(m, min, max, val); }, arg("min"), arg("max"), arg("val"), 
         "Set values in range [min, max] to val.")
    .def("replaceRange", [](SelfT &m, VxT min, VxT max, VxT val) { replaceRange(m, min, max, val); }, arg("min"), arg("max"), arg("val"),
         "Replace values in range [min, max] with val.")
    .def("write8bit", [](SelfT &m, std::string outName, double minv, double maxv) { _write8bit(m, outName, minv, maxv); }, arg("filename"), arg("min"), arg("max"),
         "Write as 8-bit image scaled between min and max.")
    .def("read", [](SelfT &m, std::string s) { m.reset(int3(0), 0); m.readFromHeader(s); }, arg("filename"), "Read image from file.")
    // //.def("readAtZ", &SelfT::readAtZ)
    .def("modeNSames", [](SelfT &m, const short nSameNei) { return modeNSames(m, nSameNei); }, arg("nSameNeighbors"), 
         "Apply mode filter based on nearest 6 neighbor voxels.")
    .def("medianFilter", [](SelfT &m) { m = median(m); }, "Apply a 1+6-neighbour median filter.")
    .def("medianX", [](SelfT &m) { m = medianx(m); }, "Apply median filter with kernel size of 1 in x-direction to reduce compressed file size")
    .def("FaceMedian06", &SelfT::FaceMedian06, arg("nAdj0"), arg("nAdj1"), 
         "Set voxel value to 0/1 if it has more than nAdj0/1 neighbours with value 0/1, in its 6 nearest voxels")
    .def("pointMedian032", &SelfT::PointMedian032, arg("nAdj0"), arg("nAdj1"), arg("lbl0"), arg("lbl1"), 
         "Set voxel value to lbl0/1 if it has more than nAdj0/1 neighbours with value lbl0/1, in its 6+26 nearest voxels")
    .def("faceMedNgrowToFrom", [](SelfT &m, VxT lblTo, VxT lblFrm, int ndif) { FaceMedGrowToFrom(m, lblTo, lblFrm, ndif); }, arg("labelTo"), arg("labelFrom"), arg("nDiff"), 
         "Face median grow to/from labels.")
    .def("delense032", [](SelfT &m, int x, int y, int r, char d, VxT v) { _delense032(m, x, y, r, d, v); }, arg("x"), arg("y"), arg("r"), arg("d"), arg("val"), "Delense operation.")
    .def("circleOut", [](SelfT &m, int x, int y, int r, char d, VxT v) { circleOut(m, x, y, r, d, v); }, arg("x"), arg("y"), arg("r"), arg("d"), arg("val"), "Circle out operation.")
    .def("growLabel", &SelfT::growLabel)
    .def("keepLargest", [](SelfT &m, VxT minvv, VxT maxvv) { keepLargestvv(m, minvv, maxvv); }, arg("min"), arg("max"), "Keep largest singly-connected region with values in [min, max].")
    //.def("maskWriteFraction", &SelfT::maskWriteFraction)
    .def("mapFrom", [](SelfT& m, const SelfT& vimg2, VxT vmin, VxT vmax, double scale, double shift) { mapToFrom(m, vimg2, vmin, vmax, scale, shift); }, arg("sourceImage"), arg("vmin"), arg("vmax"), arg("scale"), arg("shift"), "Map values from another image.") 
    .def("addSurfNoise", [](SelfT &m, const int randMask1, const int randMask2, int trsh, int randseed) { addSurfNoise(m, randMask1, randMask2, trsh, randseed); }, arg("mask1"), arg("mask2"), arg("threshold"), arg("seed"), "Add surface noise.")

    // // Individual bindings for convenience
    .def("averageWith", [](SelfT &m, std::string args) { std::stringstream ss(args); return MCTProcessing::averageWith(ss, m); })
    .def("averageWith_mBE", [](SelfT &m, std::string args) { std::stringstream ss(args); return MCTProcessing::averageWith_mBE(ss, m); })
    .def("plotAll", [](SelfT &m, std::string fnam_, int minv, int maxv, int iSlice_, int nBins,
                        std::string normalAxis, bool grey, bool color, bool histogram, bool zProfile,
                        SelfT* img2Ptr, int mina, int maxa) {
                        int colrGreyHistZprofile = grey*1 + color*2 + histogram*4 + zProfile*8;
                        return MCTProcessing::plotAll(m, minv, maxv, iSlice_, nBins, normalAxis, fnam_, colrGreyHistZprofile, img2Ptr, mina, maxa);
                    },  
        arg("name")="pltAll", arg("minv")=0, arg("maxv")=-1000001, arg("sliceIndex")=-1000000, arg("nBins")=128,
        arg("normalAxis")="xyz", arg("grey")=true, arg("color")=true, arg("histogram")=true, arg("zProfile")=true,
        arg("alphaImage")=nullptr, arg("alphaMin")=0, arg("alphaMax")=-1000001,
        "Plot all visualizations (Histogram, ZProfile, Slices) with various options, hackish for debugging, all args except are optional")
    .def("mode26", [](SelfT &m, int nMinD) { return mode26(m,nMinD); })
    .def("growingThreshold", [](SelfT &m, VxT t1, VxT t2, VxT t3, VxT t4, int nIter) {
        ::growingThreshold(m, t1, t2, t3, t4, nIter);
    }, arg("startMin"), arg("startMax"), arg("finalMin"), arg("finalMax"), arg("iterations")=4)
    .def("replaceOutSideValue", [](SelfT &m, int vo, int vnew, int nHoleSize) {
        return MCTProcessing::replaceOutSideValue(m, vo, vnew, nHoleSize);
    }, arg("val_old")=0, arg("val_new")=2, arg("hole_size")=5)
    .def("smooth", [](SelfT &m, int nItrs, int kernRad, double sigmavv, double sharpFact) {
        return MCTProcessing::smooth(m, nItrs, kernRad, sigmavv, sharpFact);
    }, arg("iterations")=1, arg("kernel_radius")=1, arg("sigma_val")=16.0, arg("sharpness")=0.1)
    .def("svgHistogram", [](SelfT &m, std::string fnam, int nBins, double minV, double maxV) {
        return MCTProcessing::svgHistogram(m, fnam, nBins, minV, maxV);
    }, arg("filename")="aa.svg", arg("bins")=128, arg("min_val")=3e38, arg("max_val")=-3e38)
    .def("svgZProfile", [](SelfT &m, std::string fnam, double minV, double maxV) {
            return MCTProcessing::svgZProfile(m, fnam, VxT(minV), VxT(maxV));
        }, arg("filename")="aa.svg", arg("min_val")=0, arg("max_val")=255) // Assuming 255 default for maxV based on typical usage, though maxT(T) is ideal
    .def("flipEndian", [](SelfT &m) { ::flipEndian(m); })
    .def("replaceRangeByImage", [](SelfT &m, double minv, double maxv, std::string fnam) {
        if(fnam.empty()) throw std::runtime_error("no image name provided");
        ::replaceRangeByImage(m, VxT(minv), VxT(maxv), SelfT(fnam, readOpt::procAndConvert));
    }, arg("min_val"), arg("max_val"), arg("image_file"))
    .def("replaceByImageRange", [](SelfT &m, double minv, double maxv, std::string fnam) {
        if(fnam.empty()) throw std::runtime_error("no image name provided");
        ::replaceByImageRange(m, VxT(minv), VxT(maxv), SelfT(fnam, readOpt::procAndConvert));
    }, arg("min_val"), arg("max_val"), arg("image_file"))
    .def("readFromFloat", [](SelfT &m, std::string header, float scale, float shift) {
        return MCTProcessing::readFromFloat(m, header, scale, shift);
    }, arg("header"), arg("scale")=1.0f, arg("shift")=0.0f)
    .def("bilateralX", [](SelfT &m, int nItrs, int kernRad, int Xstp, double sigmavv, double sharpFact, double sigmadd) {
         return ::bilateralX(m, nItrs, kernRad, Xstp, sigmavv, sharpFact, sigmadd);
    }, arg("iterations")=1, arg("kernel_radius")=1, arg("x_step")=2, arg("sigma_val")=16.0, arg("sharpness")=0.1, arg("sigma_spatial")=2.0)
    .def("bilateralGauss", [](SelfT &m, int nItrs, int kernRad, double sigmavv, double sharpFact, double sigmadd) {
         return ::bilateralGauss(m, nItrs, kernRad, sigmavv, sharpFact, sigmadd);
    }, arg("iterations")=1, arg("kernel_radius")=1, arg("sigma_val")=16.0, arg("sharpness")=0.1, arg("sigma_spatial")=2.0)
    .def("meanWide", [](SelfT &m, int nW, int noisev, int avg, int delta, int nItrs, std::string smoothImg) {
         return MCTProcessing::meanWide(m, nW, noisev, avg, delta, nItrs, smoothImg);
    }, arg("width")=0, arg("noise_val")=4, arg("average")=0, arg("delta")=20, arg("iterations")=15, arg("smooth_image")="")
    .def("otsu_threshold", [](SelfT &m, int minv, int maxv) { return ::otsu_th(m, minv, maxv); }, arg("min_val")=0, arg("max_val")=256)
    .def("dering", [](SelfT &m, int X0, int Y0, int X1, int Y1, int minV, int maxV, int nr, int ntheta, int nz) {
         return MCTProcessing::dering(m, X0, Y0, X1, Y1, minV, maxV, nr, ntheta, nz);
    }, arg("x0"), arg("y0"), arg("x1"), arg("y1"), arg("min_val")=0, arg("max_val")=255, arg("nr")=0, arg("ntheta")=18, arg("nz")=0)
    .def("adjustBrightnessWith", [](SelfT &m, std::string imgName) { return MCTProcessing::adjustBrightnessWith(m, imgName); }, arg("image_file"))
    .def("adjustSliceBrightness", [](SelfT &m, voxelImage& mskA, voxelImage& mskB, SelfT& img2, int nSmoothItr, int nSmoothKrnl) {
         return MCTProcessing::adjustSliceBrightness(m, mskA, mskB, img2, nSmoothItr, nSmoothKrnl);
    }, arg("mask_a"), arg("mask_b"), arg("ref_image"), arg("smooth_iter")=3, arg("smooth_kernel")=20)
    .def("cutOutside", [](SelfT &m, char dir, int nExtraOut, int threshold, int cuthighs, int nShiftX, int nShiftY, int outVal) {
            VoxLib::cutOutside(m, dir, nExtraOut, threshold, cuthighs, nShiftX, nShiftY, VxT(outVal)); }, 
         arg("axis")='z', arg("extra_out")=0, arg("threshold")=-1, arg("cut_highs")=0, arg("shift_x")=0, arg("shift_y")=0, arg("fill_val")=0)
    .def("variance", [](SelfT &m, int minV, int maxV) { return ::varianceDbl(m, minV, maxV); }, arg("min_val")=0, arg("max_val")=255)
    ;
    addDodgyFuncsInt(clas);
    addDodgyFuncsU8(clas);
}

} // namespace VxlPy 