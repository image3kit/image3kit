#include "VxlImgXX.hpp"

template VoxelImageT<unsigned short>::VoxelImageT(const std::string&, readOpt, int);

void bind_VxlImgU16(pybind11::module &m, const char* name) {
    VxlPy::bind_VxlImg<unsigned short>(m, name);
    VxlPy::bind_funcs<unsigned short>(m);
}
