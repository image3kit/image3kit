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


if(MSVC)
  add_compile_definitions(NOMINMAX)
endif()

include(FetchContent)

find_package(ZLIB)
if(NOT ZLIB_FOUND)
  message(STATUS "ZLIB not found, fetching from source...")
  FetchContent_Declare(
    ZLIB
    URL https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz
    EXCLUDE_FROM_ALL
  )
  FetchContent_MakeAvailable(ZLIB)
  if(NOT TARGET ZLIB::ZLIB)
    if(TARGET zlibstatic)
      add_library(ZLIB::ZLIB ALIAS zlibstatic)
    elseif(TARGET zlib)
      add_library(ZLIB::ZLIB ALIAS zlib)
    endif()
  endif()
endif()

if(APPLE)
  if(NOT DEFINED OpenMP_ROOT AND DEFINED ENV{OpenMP_ROOT})
    set(OpenMP_ROOT $ENV{OpenMP_ROOT})
  endif()
  
  find_package(OpenMP)
  
  if(NOT OpenMP_FOUND)
    # Hint for Homebrew libomp
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    if (EXISTS "${OpenMP_ROOT}/lib/libomp.dylib")
      set(OpenMP_omp_LIBRARY "${OpenMP_ROOT}/lib/libomp.dylib")
    elseif(EXISTS "/opt/homebrew/opt/libomp/lib/libomp.dylib")
       set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    elseif(EXISTS "/usr/local/opt/libomp/lib/libomp.dylib")
       set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib")
    endif()
    
    if (DEFINED OpenMP_omp_LIBRARY)
       message(STATUS "Manual config for OpenMP: ${OpenMP_omp_LIBRARY}")
       set(OpenMP_FOUND TRUE)
       add_library(OpenMP::OpenMP_CXX UNKNOWN IMPORTED)
       set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
         IMPORTED_LOCATION "${OpenMP_omp_LIBRARY}"
         INTERFACE_COMPILE_OPTIONS "${OpenMP_CXX_FLAGS}"
         INTERFACE_LINK_LIBRARIES "${OpenMP_omp_LIBRARY}"
       )
    endif()
  endif()
endif()

if(NOT OpenMP_FOUND)
  find_package(OpenMP REQUIRED)
endif()
FetchContent_Declare(
  TIFF
  URL https://gitlab.com/libtiff/libtiff/-/archive/v4.6.0/libtiff-v4.6.0.tar.gz
  EXCLUDE_FROM_ALL
)
set(tiff-tools OFF CACHE BOOL "Build libtiff tools")
set(tiff-tests OFF CACHE BOOL "Build libtiff tests")
set(tiff-contrib OFF CACHE BOOL "Build libtiff contrib")
set(tiff-docs OFF CACHE BOOL "Build libtiff docs")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
FetchContent_MakeAvailable(TIFF)

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
    TIFLIB
    ZLIB
    OpenMP
    VMMLIB
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

if(MSVC)
  target_compile_options(_core PRIVATE /openmp:llvm /diagnostics:caret)  # /W4
  target_compile_definitions(_core PRIVATE NOMINMAX)
else()
  target_compile_options(_core PRIVATE -Wall -Wextra -Wpedantic)
  # target_compile_options(_core PRIVATE -Werror)
endif()

target_link_libraries(_core PRIVATE Python::Module ZLIB::ZLIB TIFF::tiff OpenMP::OpenMP_CXX)

# The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION image3kit)


# Add targets for local includes to ensure they appear in compile_commands.json
# This fixes tooling errors in headers like typses.h
add_subdirectory(src/_include)
