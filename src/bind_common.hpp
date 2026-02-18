#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybind11/pytypes.h"
#ifdef min
#undef min
#undef max
#endif
#include "typses.h"
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


template<typename T> inline var3<T> tov3(py::tuple v) { return var3<T>(v[0].cast<T>(), v[1].cast<T>(), v[2].cast<T>()); }
template<typename T> inline py::tuple to3(var3<T> v) { return py::make_tuple(v.x, v.y, v.z); }


inline InputFile pyCastInput(py::dict dic) {
  InputFile inp;
  for (const auto &kv : dic) {
        inp.add(kv.first.cast<std::string>(), py::str(kv.second).cast<std::string>());
  }
  return inp;
}

} // namespace VxlPy
