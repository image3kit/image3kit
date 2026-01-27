#ifndef VxTypS
#include "pybind11/pytypes.h"
static_assert(false, "please define VxTyp and VxTypS");
#define VxTyp unsigned char
#define VxTypS "U8"
namespace py = pybind11;
#endif

py::class_<voxelImageT<VxTyp>>(mod, VxTypS, py::buffer_protocol())
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
    .def("__repr__", [](const voxelImageT<VxTyp> &m) {
        return "<" + std::string(VxTypS) + " shape=(" + 
               _s(m.size3()) + ")>";
    })
//.def(py::init<int, int, int, VxTyp>())
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
.def("readFromHeader", &voxelImageT<VxTyp>::readFromHeader, py::arg("header_file"), py::arg("processKeys")=1, "Read image dimensions/metadata from a header file.")
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
.def("distMapExtrude", [](voxelImageT<VxTyp> &m, py::dict dic) { m = distMapExtrude(m, pyCastInput(dic), true); }, py::arg("distMapDict"), "Extrude proportional to distance map.")
// // Process wrappers (string interface)

// // Individual bindings for convenience
.def("averageWith", [](voxelImageT<VxTyp> &m, std::string args) { std::stringstream ss(args); return MCTProcessing::averageWith(ss, m); })
.def("averageWith_mBE", [](voxelImageT<VxTyp> &m, std::string args) { std::stringstream ss(args); return MCTProcessing::averageWith_mBE(ss, m); })
.def("segment2", [](voxelImageT<VxTyp> &m, int nSegs, std::vector<int> th, std::vector<int> minSizs,
                    double noisev, double localF, double flatnes, double resolution, double gradFactor,
                    int krnl, int nItrs, int writedumps) {
        return segment2(m, nSegs, th, minSizs, noisev, localF, flatnes, resolution, gradFactor, krnl, nItrs, writedumps);
    },
    py::arg("nSegs")=2, py::arg("th")=std::vector<int>(), py::arg("minSizs")=std::vector<int>(),
    py::arg("noisev")=2., py::arg("localF")=800., py::arg("flatnes")=0.1, py::arg("resolution")=2.,
    py::arg("gradFactor")=0., py::arg("krnl")=2, py::arg("nItrs")=13, py::arg("writedumps")=0
)
.def("plotAll", [](voxelImageT<VxTyp> &m, std::string args) { std::stringstream ss(args); return MCTProcessing::plotAll(ss, m); })
;

// py::class_<instream>(mod, "instream", py::buffer_protocol())
// .def(py::init([](py::object kwargs) { return new instream(kwargs); }))
// ;
