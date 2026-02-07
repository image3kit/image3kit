cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(voxlib LANGUAGES CXX)

set(USE_ZLIB 1)
set(USE_TIFF 1)
set(USE_SVG 1)
set(USE_OMP 1)
set(USE_PYB 1)
set(USE_PNG 1)
set(CSD ${CMAKE_CURRENT_LIST_DIR})

set(SRCS vxlPro1.cpp vxlPro2.cpp vxlPro.cpp voxelImage.cpp)

get_test_property(${PROJECT_NAME} WORKING_DIRECTORY test_dir)
message("My test's working directory: ${test_dir}")
