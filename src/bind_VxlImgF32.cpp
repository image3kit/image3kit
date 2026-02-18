#include "VxlImgXX.hpp"

template void voxelField<float>::reset(var3<int>, float);
template voxelImageT<float>::voxelImageT(const std::string&, readOpt);
template void voxelField<var3<float>>::reset(var3<int>, var3<float>);

void bind_VxlImgF32(pybind11::module &m, const char* name) {
    VxlPy::bind_VxlImg<float>(m, name);
}
