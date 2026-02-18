#include "VxlImgXX.hpp"

template void voxelField<int>::reset(var3<int>, int);
template voxelImageT<int>::voxelImageT(const std::string&, readOpt);

void bind_VxlImgI32(pybind11::module &m, const char* name) {
    VxlPy::bind_VxlImg<int>(m, name);
}
