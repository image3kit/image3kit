
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "voxelImageI.h"
#include "shapeToVoxel.h"
#include "InputFile.h"
#include "VxlStrips.h"

#include "voxelImageProcess.h"
#include "voxelNoise.h"

namespace py = pybind11;

inline dbl3 tpl2d3(py::tuple v) { return dbl3(v[0].cast<double>(), v[1].cast<double>(), v[2].cast<double>()); }

InputFile pyCastInput(py::dict dic) {
    InputFile inp;
    for (const auto &kv : dic) {
        inp.add(kv.first.cast<std::string>(), py::str(kv.second).cast<std::string>());
    }
    return inp;
}

PYBIND11_MODULE(_core, mod, py::mod_gil_not_used()) {
    using namespace MCTProcessing;

    py::class_<var3<int>>(mod, "int3")
    .def(py::init<>())
    .def(py::init<int, int, int>())
    .def(py::init([](py::tuple t) { return var3<int>(t[0].cast<int>(), t[1].cast<int>(), t[2].cast<int>()); }))
    .def_readwrite("x", &var3<int>::x)
    .def_readwrite("y", &var3<int>::y)
    .def_readwrite("z", &var3<int>::z)
    .def("__repr__", [](const var3<int> &v) {
        return "int3(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ")";
    });

    py::class_<var3<double>>(mod, "dbl3")
    .def(py::init<>())
    .def(py::init<double, double, double>())
    .def(py::init([](py::tuple t) { return var3<double>(t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>()); }))
    .def_readwrite("x", &var3<double>::x)
    .def_readwrite("y", &var3<double>::y)
    .def_readwrite("z", &var3<double>::z)
    .def("__repr__", [](const var3<double> &v) {
        return "dbl3(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ")";
    });

    py::class_<voxelImageTBase>(mod, "voxelImageTBase");

    //py::class_<dbl3>(mod, "dbl3")
    //.def(py::init<double, double, double>())
    //.def(py::init<>())
    //.def(py::init([](py::tuple tuple) {
        //return dbl3(tuple[0].cast<double>(), tuple[1].cast<double>(), tuple[2].cast<double>());
    //}));

    py::class_<shape>(mod, "shape");

    py::class_<sphere,shape>(mod, "sphere")
    .def(py::init([](py::tuple tpl, double r, int val) {
        return sphere(tpl2d3(tpl), r, val); }))
    ;

    py::class_<cylinder,shape>(mod, "cylinder")
    .def(py::init([](py::tuple p1, py::tuple p2, double r, int val) {
        return cylinder(tpl2d3(p1), tpl2d3(p2), r, val); }))
    ;

    py::class_<kube,shape>(mod, "cube")
    .def(py::init([](py::tuple p1, py::tuple size, int val) {
        return kube(tpl2d3(p1), tpl2d3(size), val); }))
    ;


    py::class_<InputFile>(mod, "Input")
    .def(py::init([](py::dict dic) { return pyCastInput(dic); }))
    .def("add", [](InputFile &inp, std::string key, std::string val) { inp.add(key, val); })
    .def("set", [](InputFile &inp, std::string key, std::string val) { inp.set(key, val); })
    .def("get", [](InputFile &inp, std::string key, std::string val) { return inp.kwrd(key); })
    .def("setDefault", &InputFile::setDefault)
    .def("echoKeywords", [](InputFile &inp) { inp.echoKeywords(); })
    .def("renameKeys", &InputFile::renameKeys)
    ;

#define VxTyp unsigned char
#define VxTypS "VxlImgU8"
#include "./vxlPyUX.hpp"

    mod.def("readImage", [](py::object filename, int processKeys) {
        return readImage(py::str(filename).cast<std::string>(), processKeys);
    }, py::arg("filename"), py::arg("processKeys")=1, "Global helper to read an image from a file.");

    mod.def("readImageU8", [](py::dict dic) { return readImageU8(pyCastInput(dic)); }, py::arg("dict"), "Global helper to read an uint8 (3D) image from a file.");

    mod.doc() = //"voxelImage of type " VxTypS
       R"pbdoc( -----------------------

        .. currentmodule:: VoxelImage

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

#ifdef VERSION_INFO
    mod.attr("__version__") = TOSTRING(VERSION_INFO);
#else
    mod.attr("__version__") = "dev";
#endif
}
