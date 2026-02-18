#include "VxlImgXX.hpp"

template void voxelField<unsigned char>::reset(var3<int>, unsigned char);
template voxelImageT<unsigned char>::voxelImageT(const std::string&, readOpt);

void bind_VxlImgU8(pybind11::module &m, const char* name) {
    VxlPy::bind_VxlImg<unsigned char>(m, name);
}
