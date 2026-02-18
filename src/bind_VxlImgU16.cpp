#include "VxlImgXX.hpp"

template void voxelField<unsigned short>::reset(var3<int>, unsigned short);
template voxelImageT<unsigned short>::voxelImageT(const std::string&, readOpt);

void bind_VxlImgU16(pybind11::module &m, const char* name) {
    VxlPy::bind_VxlImg<unsigned short>(m, name);
}
