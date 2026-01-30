# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 4.0
cmake_minimum_required(VERSION 3.15...4.0)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

# Find pybind11 and Python
find_package(pybind11 CONFIG REQUIRED)

find_package(ZLIB REQUIRED)

# Add a library using FindPython's tooling (pybind11 also provides a helper like this)
file(GLOB PyBind_SRCS CONFIGURE_DEPENDS "src/*.cpp")
set(VOXLIB_SRC 
   src/_voxlib/voxelImage.cpp
   src/_voxlib/vxlPro2.cpp
   src/_voxlib/voxelGrad.cpp
   src/_voxlib/vxlPro1.cpp
   src/_voxlib/vxlPro.cpp
   src/_voxlib/VxlStrips.cpp
)
python_add_library(_core MODULE ${VOXLIB_SRC} ${PyBind_SRCS} WITH_SOABI)
target_compile_definitions(_core PRIVATE
    SVG
    LPNG
    ZLIB
    DUMMY=1
)
# find_package(pybind11 CONFIG REQUIRED PATHS ./build)
target_link_libraries(_core PRIVATE pybind11::headers)

target_include_directories(_core PRIVATE src/_include)
target_include_directories(_core PRIVATE src/_voxlib)
target_include_directories(_core PRIVATE pkgs/optim)
target_include_directories(_core PRIVATE pkgs/optim/primitives)
target_include_directories(_core PRIVATE pkgs/stb)
target_include_directories(_core PRIVATE pkgs/zlib)
target_include_directories(_core PRIVATE pkgs/svplot)

# This is passing in the version as a define just as an example
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

set_property(TARGET _core PROPERTY CXX_STANDARD 23)
# Toolset: gcc-linux:
set(CMAKE_CXX_FLAGS_INIT "-Werror -Wall -Wextra -Wpedantic")
# target_link_options(_core PRIVATE "-Wl,-z,defs")
target_link_libraries(_core PRIVATE Python::Module ZLIB::ZLIB)

# The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION image3kit)


# Add targets for local includes to ensure they appear in compile_commands.json
# This fixes tooling errors in headers like typses.h
add_subdirectory(src/_include)
