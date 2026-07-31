#include "VxlImgXX.hpp"

template VoxelImageT<float>::VoxelImageT(const std::string&, readOpt, int);

void bind_VxlImgF32(pybind11::module &m, const char* name) {
    VxlPy::bind_VxlImg<float>(m, name);
    VxlPy::bind_funcs<float>(m);
}
