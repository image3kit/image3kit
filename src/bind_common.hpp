#pragma once
#include <pybind11/pybind11.h>
#ifndef _POCKETPY
#include <pybind11/numpy.h>
#include "pybind11/pytypes.h"
#endif
#include <pybind11/stl.h>
#ifdef min
#undef min
#undef max
#endif
#include "typses.h"
#include "InputFile.h"

namespace VxlPy {
namespace py = pybind11;


#ifdef _POCKETPY
// pkbind has no py::sequence / py::make_tuple; VxlImgXX.hpp's only caller already passes
// a py::tuple, which is-a py::sequence for the real-pybind11 build below.
template<typename T> inline var3<T> tov3(py::tuple v) { return var3<T>(v[0].cast<T>(), v[1].cast<T>(), v[2].cast<T>()); }
template<typename T> inline py::tuple tot3(var3<T> v) { py::tuple t(3); t[0]=v.x; t[1]=v.y; t[2]=v.z; return t; }
inline py::tuple make_tuple3(int x, int y, int z) { py::tuple t(3); t[0]=x; t[1]=y; t[2]=z; return t; }
#else
template<typename T> inline var3<T> tov3(py::sequence v) { return var3<T>(v[0].cast<T>(), v[1].cast<T>(), v[2].cast<T>()); }
template<typename T> inline py::tuple tot3(var3<T> v) { return py::make_tuple(v.x, v.y, v.z); }
inline py::tuple make_tuple3(int x, int y, int z) { return py::make_tuple(x, y, z); }
#endif


inline InputFile pyCastInput(py::dict dic) {
  InputFile inp;
  for (const auto &kv : dic) {
        inp.add(kv.first.cast<std::string>(), py::str(kv.second).cast<std::string>());
  }
  return inp;
}

} // namespace VxlPy
