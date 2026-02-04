
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "shapeToVoxel.h"
#include "VxlImgXX.hpp"

namespace py = pybind11;

// Explicit instantiation to ensure symbol is generated, WTF: crap macOS clang++
template void voxelField<unsigned char>::reset(var3<int>, unsigned char);
template void voxelField<unsigned short>::reset(var3<int>, unsigned short);
template void voxelField<int>::reset(var3<int>, int);
template void voxelField<float>::reset(var3<int>, float);
template void voxelField<var3<float>>::reset(var3<int>, var3<float>);

template voxelImageT<unsigned char>::voxelImageT(const std::string&, readOpt);
template voxelImageT<unsigned short>::voxelImageT(const std::string&, readOpt);
template voxelImageT<int>::voxelImageT(const std::string&, readOpt);
template voxelImageT<float>::voxelImageT(const std::string&, readOpt);

PYBIND11_MODULE(_core, mod, py::mod_gil_not_used()) {
    using namespace VxlPy;

    // **************** sirun submodule ***************

    auto sirun = mod.def_submodule("sirun", "The sirun submodule");

    // using namespace MCTProcessing;

    py::class_<var3<int>>(sirun, "int3")
    .def(py::init<>())
    .def(py::init<int, int, int>())
    .def(py::init([](py::tuple t) { return tov3<int>(t); }))
    .def_readwrite("x", &var3<int>::x)
    .def_readwrite("y", &var3<int>::y)
    .def_readwrite("z", &var3<int>::z)
    .def("__getitem__", [](const var3<int> &v, int iVal) {
        int i = iVal;
        if (i < 0) i += 3;
        if (i < 0 || i >= 3) throw py::index_error();
        return v[i];
    })
    .def("__setitem__", [](var3<int> &v, int i, int val) { v[i] = val; })
    .def("__len__", [](const var3<int> &) { return 3; })
    .def("__repr__", [](const var3<int> &v) { return "int3(" + _s(v) + ")"; })
    ;

    py::class_<var3<double>>(sirun, "dbl3")
    .def(py::init<>())
    .def(py::init<double, double, double>())
    .def(py::init([](py::tuple t) { return tov3<double>(t); }))
    .def_readwrite("x", &var3<double>::x)
    .def_readwrite("y", &var3<double>::y)
    .def_readwrite("z", &var3<double>::z)
    .def("__getitem__", [](const var3<double> &v, int iVal) {
        int i = iVal;
        if (i < 0) i += 3;
        if (i < 0 || i >= 3) throw py::index_error();
        return v[i];
    })
    .def("__setitem__", [](var3<double> &v, int i, double val) { v[i] = val; })
    .def("__len__", [](const var3<double> &) { return 3; })
    .def("__repr__", [](const var3<double> &v) { return "dbl3(" + _s(v) + ")"; })
    ;


    py::class_<InputFile>(sirun, "Input")
    .def(py::init([](py::dict dic) { return pyCastInput(dic); }))
    .def("add", [](InputFile &inp, std::string key, std::string val) { inp.add(key, val); })
    .def("set", [](InputFile &inp, std::string key, std::string val) { inp.set(key, val); })
    .def("get", [](InputFile &inp, std::string key) { return inp.kwrd(key); })
    .def("setDefault", &InputFile::setDefault)
    .def("echoKeywords", [](InputFile &inp) { inp.echoKeywords(); })
    .def("renameKeys", &InputFile::renameKeys)
    ;


    // **************** sirun submodule ***************

    auto voxlib = mod.def_submodule("voxlib", "The VxlImg template classes");

    auto _voxelImageTBase = py::class_<voxelImageTBase>(voxlib, "voxelImageTBase")
    ;

    auto _shape = py::class_<shape>(voxlib, "shape");

    py::class_<sphere,shape>(voxlib, "sphere")
    .def(py::init([](py::tuple tpl, double r, int val) {
        return sphere(tov3<double>(tpl), r, val); }))
    ;

    py::class_<cylinder,shape>(voxlib, "cylinder")
    .def(py::init([](py::tuple p1, py::tuple p2, double r, int val) {
        return cylinder(tov3<double>(p1), tov3<double>(p2), r, val); }),
        py::arg("p1"), py::arg("p2"), py::arg("r"), py::arg("val"),
        "p1: first point on axis, p2: second point on axis, r: radius, val: paint value")
    ;

    py::class_<kube,shape>(voxlib, "cube")
    .def(py::init([](py::tuple p1, py::tuple size, int val) {
        std::cout << " p1: " << tov3<double>(p1)  << " size: " << tov3<double>(size)  << " val: " << val << " ";
        return kube(tov3<double>(p1), tov3<double>(size), val); }),
        py::arg("p1"), py::arg("size"), py::arg("val"),
        "p1: first point, size: size of cuboid sides, val: paint value")
    ;

    // TODO switch to int32_t...
    bind_VxlImg<unsigned char>(voxlib, "VxlImgU8");
    bind_VxlImg<unsigned short>(voxlib, "VxlImgU16");
    bind_VxlImg<int>(voxlib, "VxlImgI32");
    bind_VxlImg<float>(voxlib, "VxlImgF32");

    voxlib.def("labelImage", [](voxelImage &m, double minvv, double maxvv) {
        auto lbls = labelImage(m, (unsigned char)(minvv), (unsigned char)(maxvv));
	    compressLabelImage(lbls);
        return lbls;
    });

    voxlib.def("readImageBase", [](py::object filename, int processKeys) {
        return readImage(py::str(filename).cast<std::string>(), processKeys);
    }, py::arg("filename"), py::arg("processKeys")=1, "Global helper to read an image from a file.");



    // Bind docstrings or versions to the main module or submodules as needed
    mod.doc() = "The _core (PyBind11 wrapper) of image3 module of containing sirun and VxlImg submodules.";

#ifdef VERSION_INFO
    mod.attr("__version__") = TOSTRING(VERSION_INFO);
#else
    mod.attr("__version__") = "dev";
#endif
}
