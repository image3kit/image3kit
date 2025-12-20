#ifndef VxTypS
error please define VxTyp and VxTypS
#endif

py::class_<voxelImageT<VxTyp>>(mod, VxTypS, py::buffer_protocol())
.def_buffer([](voxelImageT<VxTyp> &m) -> py::buffer_info {
        return py::buffer_info(
            m.data(),                             // Pointer to buffer
            sizeof(VxTyp),                          // Size of one scalar
            py::format_descriptor<VxTyp>::format(), // Python struct-style format descriptor
            3,                                    // Number of dimensions
            { size_t(m.nz()), size_t(m.ny()), size_t(m.nx()) },             // Dimensions
            { long(sizeof(VxTyp) * m.nxy()), long(sizeof(VxTyp) *m.nx()), long(sizeof(VxTyp)) } // Strides (in bytes)
        );
    })
    .def(py::init([](py::tuple tpl, VxTyp value) { return voxelImageT<VxTyp>(tpl[0].cast<int>(), tpl[1].cast<int>(), tpl[2].cast<int>(), value); }))
//.def(py::init<int, int, int, VxTyp>())
.def("data", [](voxelImageT<VxTyp> &m) {
    return py::array_t<VxTyp>(py::buffer_info{
        m.data(),
        sizeof(VxTyp),
        py::format_descriptor<VxTyp>::format(),
        3,
        { size_t(m.nz()), size_t(m.ny()), size_t(m.nx()) },             // Dimensions
        { long(sizeof(VxTyp) * m.nxy()), long(sizeof(VxTyp) *m.nx()), long(sizeof(VxTyp)) } // Strides (in bytes)
    }); })
.def("nx", &voxelImageT<VxTyp>::nx)
.def("ny", &voxelImageT<VxTyp>::ny)
.def("nz", &voxelImageT<VxTyp>::nz)
.def("printInfo", &voxelImageT<VxTyp>::printInfo)
.def("shape", [&](voxelImageT<VxTyp> &m) { return py::make_tuple(m.nz(), m.nz(), m.nz()); })
.def("write", &voxelImageT<VxTyp>::write)
.def("writeNoHdr", &voxelImageT<VxTyp>::writeNoHdr)
// .def("writeHeader", &voxelImageT<VxTyp>::writeHeader)
.def("readFromHeader", &voxelImageT<VxTyp>::readFromHeader)
.def("readAscii", &voxelImageT<VxTyp>::readAscii)
// .def("readRLE", &voxelImageT<VxTyp>::readRLE)
.def("cropD", &voxelImageT<VxTyp>::cropD)
.def("growBox", &voxelImageT<VxTyp>::growBox)
.def("shrinkBox", &voxelImageT<VxTyp>::shrinkBox)
.def("shrinkBox", &voxelImageT<VxTyp>::shrinkBox)
.def("fillHoles", &voxelImageT<VxTyp>::fillHoles)
.def("threshold101", &voxelImageT<VxTyp>::threshold101)
.def("writeAConnectedPoreVoxel", &voxelImageT<VxTyp>::writeAConnectedPoreVoxel)
.def("AND", &voxelImageT<VxTyp>::AND)
.def("NOT", &voxelImageT<VxTyp>::NOT)
.def("OR", &voxelImageT<VxTyp>::OR)
.def("XOR", &voxelImageT<VxTyp>::XOR)
.def("resampleMode", [&](voxelImageT<VxTyp> &m, double nReSampleNotSafe) { return resampleMode(m,nReSampleNotSafe); })
.def("resampleMax", [&](voxelImageT<VxTyp> &m, double nReSampleNotSafe) { return resampleMax(m,nReSampleNotSafe); })
.def("resampleMean", [&](voxelImageT<VxTyp> &m, double nReSampleNotSafe) { return resampleMean(m,nReSampleNotSafe); })
.def("resliceZ", [&](voxelImageT<VxTyp> &m, double nReSampleNotSafe) { return resliceZ(m,nReSampleNotSafe); })
.def("Paint", [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, setIn); })
.def("Paint"         , [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, setIn); })
.def("PaintAdd"      , [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, addTo); })
.def("PaintBefore"   , [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, setBefor); })
.def("PaintAfter"    , [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, setAfter); })
.def("PaintAddBefore", [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, addBefor); })
.def("PaintAddAfter" , [&](voxelImageT<VxTyp> &m, const shape& sh) { _SHAPERATEPy(m, sh, addAfter); })
.def("writeContour"  , [&](voxelImageT<VxTyp> &m, const string& outSurf) {
    InputFile inp;
    if (outSurf.size())  inp.set("outputSurface", outSurf);
    // vxlToSurfWrite(m, inp);
    })
.def("sliceToPng"    , [&](voxelImageT<VxTyp> &m, const string& normalAxis, const string& fnam, int iSlice, int bgnv, int endv, const string& color) {
    ::sliceToPng(m, normalAxis, fnam, iSlice, bgnv, endv, color);
    })
;

mod.def("readImage", &readImage);

// py::class_<instream>(mod, "instream", py::buffer_protocol())
// .def(py::init([](py::object kwargs) { return new instream(kwargs); }))
// ;
