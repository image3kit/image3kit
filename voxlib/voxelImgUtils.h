#pragma once

#include "typses.h"

#include <array>
#include <string>

namespace VoxLib {

const std::string& imgExt(const std::string& defSuffix=""); //!< get and set default image format
void fileNameToImgInfo(const std::string& fname, int3& nnn, dbl3& dx_);
std::string fileNameToElemType(const std::string& fname);
std::string getAmiraDataType(const std::string& fnam);
void getAmiraHeaderSize(const std::string& fnam, int3& nnn, dbl3& dx_, dbl3& X0_, int& nSkipBytes, int& RLE);
void getMhdHeaderData(const std::string& hdrNam, std::string& fnam, int3& nnn, dbl3& dx_, dbl3& X0_,
                      int& nSkipBytes, std::string& BinaryData, std::string& flipSigByt,
                      double& unit_, bool& X0read, bool& dxread, bool& autoUnit,
                      int maxNz, size_t elemSize = 0);

std::array<double,5> otsu_minAvThresholdAvMx(const piece<double>& hist, int  ibgn, int  iend, double shift, double scale);

}

using VoxLib::imgExt;
using VoxLib::otsu_minAvThresholdAvMx;
