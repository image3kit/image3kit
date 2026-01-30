
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "shapeToVoxel.h"
#include "InputFile.h"
#include "vxlPyUX.hpp"

namespace py = pybind11;

inline dbl3 tpl2d3(py::tuple v) { return dbl3(v[0].cast<double>(), v[1].cast<double>(), v[2].cast<double>()); }

PYBIND11_MODULE(_core, mod, py::mod_gil_not_used()) {
    using namespace VxlPy;

    // **************** sirun submodule ***************

    auto sirun = mod.def_submodule("sirun", "The sirun submodule");

    // using namespace MCTProcessing;

    py::class_<var3<int>>(sirun, "int3")
    .def(py::init<>())
    .def(py::init<int, int, int>())
    .def(py::init([](py::tuple t) { return var3<int>(t[0].cast<int>(), t[1].cast<int>(), t[2].cast<int>()); }))
    .def_readwrite("x", &var3<int>::x)
    .def_readwrite("y", &var3<int>::y)
    .def_readwrite("z", &var3<int>::z)
    .def("__repr__", [](const var3<int> &v) { return "int3(" + _s(v) + ")"; })
    ;

    py::class_<var3<double>>(sirun, "dbl3")
    .def(py::init<>())
    .def(py::init<double, double, double>())
    .def(py::init([](py::tuple t) { return var3<double>(t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>()); }))
    .def_readwrite("x", &var3<double>::x)
    .def_readwrite("y", &var3<double>::y)
    .def_readwrite("z", &var3<double>::z)
    .def("__repr__", [](const var3<double> &v) { return "dbl3(" + _s(v) + ")"; })
    ;


    py::class_<InputFile>(sirun, "Input")
    .def(py::init([](py::dict dic) { return pyCastInput(dic); }))
    .def("add", [](InputFile &inp, std::string key, std::string val) { inp.add(key, val); })
    .def("set", [](InputFile &inp, std::string key, std::string val) { inp.set(key, val); })
    .def("get", [](InputFile &inp, std::string key, std::string val) { return inp.kwrd(key); })
    .def("setDefault", &InputFile::setDefault)
    .def("echoKeywords", [](InputFile &inp) { inp.echoKeywords(); })
    .def("renameKeys", &InputFile::renameKeys)
    ;


    // **************** sirun submodule ***************

    auto voxlib = mod.def_submodule("voxlib", "The VxlImg template classes");

    py::class_<voxelImageTBase>(voxlib, "voxelImageTBase")
    ;

    py::class_<shape>(voxlib, "shape");

    py::class_<sphere,shape>(voxlib, "sphere")
    .def(py::init([](py::tuple tpl, double r, int val) {
        return sphere(tpl2d3(tpl), r, val); }))
    ;

    py::class_<cylinder,shape>(voxlib, "cylinder")
    .def(py::init([](py::tuple p1, py::tuple p2, double r, int val) {
        return cylinder(tpl2d3(p1), tpl2d3(p2), r, val); }))
    ;

    py::class_<kube,shape>(voxlib, "cube")
    .def(py::init([](py::tuple p1, py::tuple size, int val) {
        return kube(tpl2d3(p1), tpl2d3(size), val); }))
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
