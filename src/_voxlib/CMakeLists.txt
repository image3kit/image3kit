
cmake_minimum_required(VERSION 3.28 FATAL_ERROR)
project(voxlib LANGUAGES CXX)

set(USE_ZLIB 1)
set(USE_TIFF 1)
set(USE_SVG 1)
set(USE_OMP 1)
set(USE_PYB 1)
set(USE_PNG 1)
set(CSD ${CMAKE_CURRENT_LIST_DIR})
include(${CSD}/../../src/script/config.cmake)

set(SRCS vxlPro1.cpp  vxlPro2.cpp  vxlPro.cpp  voxelImage.cpp)

add_exe(voxelImageProcess  voxelImageProcess.cpp  ${SRCS})
set_inc(voxelImageProcess  ${msSrc}/include ${CSD}/../primitives  ${CSD}/../optim)

add_exe(voxelImageProcPy  voxelImageProcess.cpp  ${SRCS})
set_inc(voxelImageProcPy  ${msSrc}/include ${CSD}/../primitives  ${CSD}/../optim)

add_exe(voxelGrad  voxelGrad.cpp ${SRCS})
set_inc(voxelGrad  ${msSrc}/include )
set_def(voxelGrad  _VoxBasic )

add_exe(performance  performance.cpp)
set_inc(performance  ${msSrc}/include )

get_test_property(${PROJECT_NAME} WORKING_DIRECTORY test_dir)
message("My test's working directory: ${test_dir}")
