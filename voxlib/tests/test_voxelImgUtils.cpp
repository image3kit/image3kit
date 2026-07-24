#include "utest.h"
#include "voxelImgUtils.h"
#include <fstream>
#include <cstdio>
#include <vector>
#include <string>

UTEST(VoxelImgUtils, FileNameToElemTypeStaticNames) {
    EXPECT_STREQ("MET_UCHAR", VoxLib::fileNameToElemType("test_8bit.raw").c_str());
    EXPECT_STREQ("MET_USHORT", VoxLib::fileNameToElemType("test_16bit.raw").c_str());
    EXPECT_STREQ("MET_INT", VoxLib::fileNameToElemType("test_32bit.raw").c_str());
    EXPECT_STREQ("MET_INT", VoxLib::fileNameToElemType("test_32bit_3p4um.raw").c_str());
    EXPECT_STREQ("MET_FLOAT", VoxLib::fileNameToElemType("test_float.raw").c_str());
    EXPECT_STREQ("MET_INT", VoxLib::fileNameToElemType("sample_VElems.raw").c_str());
}

UTEST(VoxelImgUtils, MhdHeaderParsing) {
    const char* mhdPath = "test_synth_header.mhd";
    {
        std::ofstream mhd(mhdPath);
        mhd << "ObjectType = Image\n";
        mhd << "NDims = 3\n";
        mhd << "ElementType = MET_SHORT\n";
    }

    std::string elemType = VoxLib::fileNameToElemType(mhdPath);
    std::remove(mhdPath);

    EXPECT_STREQ("MET_SHORT", elemType.c_str());
}

UTEST(VoxelImgUtils, RawSizeDetection) {
    const char* rawPath = "test_01_synth_10x10x10.raw";
    {
        // 10*10*10 = 1000 voxels; 2 bytes per voxel => 2000 bytes total size
        std::ofstream raw(rawPath, std::ios::binary);
        std::vector<char> buffer(2000, 0);
        raw.write(buffer.data(), buffer.size());
    }

    std::string elemType = VoxLib::fileNameToElemType(rawPath);
    std::remove(rawPath);

    EXPECT_STREQ("MET_USHORT", elemType.c_str());
}

UTEST(VoxelImgUtils, FileNameToImgInfo) {
    int3 nnn{}; dbl3 dx{};
    VoxLib::fileNameToImgInfo("01_tst_2D_Extruded_240x200x28_5p0um.raw", nnn, dx);
    EXPECT_EQ(nnn, int3(240,200,28));
    EXPECT_NEAR(dx[0], 5.0, 1e-7);
}
