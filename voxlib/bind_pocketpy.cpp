/*-------------------------------------------------------------------------*\

Embeds pocketpy and exposes image3kit's real VxlImg bindings (VxlImgXX.hpp) to it, so
`.mhd` headers can embed a `VxlPy { <script> }` block that runs against the image being
loaded. Replaces the old voxelplugins.cpp/vxlPro*.cpp keyword-command DSL.

\*-------------------------------------------------------------------------*/

#include "VxlImgXX.hpp"
#include "shapeToVoxel.h"
#include <sstream>
#include <vector>
#include <algorithm>

using namespace VxlPy;
using namespace VoxLib;
namespace py = pybind11;

PYBIND11_EMBEDDED_MODULE(VxlPy, mod) {
    py::class_<var3<int>>(mod, "int3")
    .def(py::init<>())
    .def(py::init<int, int, int>())
    .def(py::init([](py::tuple s) { return tov3<int>(s); }))
    .def_readwrite("x", &var3<int>::x)
    .def_readwrite("y", &var3<int>::y)
    .def_readwrite("z", &var3<int>::z)
    ;

    py::class_<var3<double>>(mod, "dbl3")
    .def(py::init<>())
    .def(py::init<double, double, double>())
    .def(py::init([](py::tuple s) { return tov3<double>(s); }))
    .def_readwrite("x", &var3<double>::x)
    .def_readwrite("y", &var3<double>::y)
    .def_readwrite("z", &var3<double>::z)
    ;

    py::class_<VoxelImagesBase>(mod, "VoxelImagesBase")
    .def("write", &VoxelImagesBase::write, py::arg("filename"))
    .def("print_info", &VoxelImagesBase::printInfo)
    ;

    py::class_<shape>(mod, "shape");

    py::class_<sphere,shape>(mod, "sphere")
    .def(py::init([](py::tuple tpl, double r, int val) {
        return sphere(tov3<double>(tpl), r, val); }),
        py::arg("center"), py::arg("r"), py::arg("val"))
    ;

    py::class_<cylinder,shape>(mod, "cylinder")
    .def(py::init([](py::tuple p1, py::tuple p2, double r, int val) {
        return cylinder(tov3<double>(p1), tov3<double>(p2), r, val); }),
        py::arg("p1"), py::arg("p2"), py::arg("r"), py::arg("val"),
        "p1: first point on axis, p2: second point on axis, r: radius, val: paint value")
    ;

    py::class_<kube,shape>(mod, "cube")
    .def(py::init([](py::tuple p1, py::tuple size, int val) {
        return kube(tov3<double>(p1), tov3<double>(size), val); }),
        py::arg("p1"), py::arg("size"), py::arg("val"),
        "p1: first point, size: size of cuboid sides, val: paint value")
    ;

    py::class_<triangular,shape>(mod, "triangular")
    .def(py::init([](py::tuple po, double L1, double L2, double h, double Lt, double ch, int val) {
        return triangular(tov3<double>(po), L1, L2, h, Lt, ch, val); }),
        py::arg("po"), py::arg("L1"), py::arg("L2"), py::arg("h"), py::arg("Lt"), py::arg("ch"), py::arg("val"),
        "po: apex point, L1/L2: half-widths, h: height, Lt: throat length, ch: contraction ratio, val: paint value")
    ;

    VxlPy::bind_VxlImg<unsigned char>(mod, "VxlImgU8");    VxlPy::bind_funcs<unsigned char>(mod);
    VxlPy::bind_VxlImg<unsigned short>(mod, "VxlImgU16");  VxlPy::bind_funcs<unsigned short>(mod);
    VxlPy::bind_VxlImg<int>(mod, "VxlImgI32");             VxlPy::bind_funcs<int>(mod);
    VxlPy::bind_VxlImg<float>(mod, "VxlImgF32");           VxlPy::bind_funcs<float>(mod);
}

namespace VxlPy {

//! lazily initialize a persistent pocketpy interpreter (avoid repeated init/finalize
//! cost across many image loads within one process). Not thread-safe: header-script
//! processing happens per-process in observed usage (MPI-parallel, not OpenMP-parallel).
static void ensureInterpreter() {
    static py::scoped_interpreter guard{};
}

//! InputFile's brace-capture keeps a "VxlPy { ... }" block's text verbatim, including
//! whatever leading whitespace the header author used for readability - Python (and
//! pocketpy) reject a top-level script whose first line is indented with no enclosing
//! block, so strip the common leading-whitespace prefix shared by all non-blank lines
//! (same idea as Python's textwrap.dedent) before executing.
std::string dedent(const std::string& script) {
    std::istringstream in(script);
    std::vector<std::string> lines;
    std::string line;
    size_t commonIndent = std::string::npos;
    while (std::getline(in, line)) {
        lines.push_back(line);
        size_t firstNonSpace = line.find_first_not_of(" \t");
        if (firstNonSpace != std::string::npos)
            commonIndent = std::min(commonIndent, firstNonSpace);
    }
    if (commonIndent == std::string::npos || commonIndent == 0) return script;
    std::string result;
    for (const auto& l : lines) {
        result += (l.size() > commonIndent ? l.substr(commonIndent) : "");
        result += "\n";
    }
    return result;
}

template<typename T>
static bool execHeaderScriptT(const std::string& script, VoxelImagesBase* imgPtr) {
    if (auto img = dynamic_cast<VoxelImageT<T>*>(imgPtr)) {
        ensureInterpreter();
        py::dict locals;
        locals["img"] = py::cast(img, py::return_value_policy::reference);
        // bring VxlImgU8/sphere/cylinder/cube/int3/dbl3/etc. into scope unqualified,
        // so header scripts can write `img.paint(sphere(...))` without an import line.
        py::exec("from VxlPy import *", py::globals(), locals);
        py::exec(script, py::globals(), locals);
        return true;
    }
    return false;
}

//! Replaces the old vxlProcess<InpT, SupportedVoxTyps>(...) dispatch: runs `script`
//! (already-extracted text, e.g. from a "VxlPy { ... }" block) as a pocketpy script
//! against the given image, bound as `img`.
void execScript(const std::string& rawScript, VoxelImagesBase* imgPtr, const std::string&) {
    if (rawScript.empty()) return;
    std::string script = dedent(rawScript);
    bool ran = execHeaderScriptT<unsigned char>(script, imgPtr)
            || execHeaderScriptT<unsigned short>(script, imgPtr)
            || execHeaderScriptT<int>(script, imgPtr)
            || execHeaderScriptT<float>(script, imgPtr);
    if (!ran) std::cout << "Unknown image type." << std::endl;
}

//! Finds the "VxlPy" keyword's raw (brace-captured, verbatim multi-line) text in a
//! header's InputFile and executes it via execScript() - the direct replacement for
//! readFromHeader()'s old vxlProcess<InputFile, SupportedVoxTyps>(...) call.
void execHeaderScript(const InputFile& inp, VoxelImagesBase* imgPtr, const std::string& nam) {
    execScript(inp.kwrd("VxlPy"), imgPtr, nam);
}

} // namespace VxlPy
