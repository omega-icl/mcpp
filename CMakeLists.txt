cmake_minimum_required(VERSION 3.18)
project(mcpp LANGUAGES CXX)

# Prefer BoostConfig.cmake with newer CMake/Boost versions.
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 1. Compilation Modes
message(STATUS "--------------------------------------------------------")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build" FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Setting build type to 'Release'")
else()
  message(STATUS "Setting build type to 'Debug'")
endif()
message(STATUS "--------------------------------------------------------")

# Set specific flags if needed, although standard CMake flags are usually
# sufficient. Debug: -g is standard. -O0 is often default for Debug but we can
# ensure it. Release: -O2 is standard. User asked for high optimization. Match
# the classical Make setup: Debug uses -O0 -g; Release uses -O2 without defining
# NDEBUG. CMake's default Release flags often include -O3 -DNDEBUG, so set the
# cache entries explicitly.
set(CMAKE_CXX_FLAGS_DEBUG
    "-O0 -g"
    CACHE STRING "Debug flags" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE
    "-O2"
    CACHE STRING "Release flags" FORCE)

# Options
option(ENABLE_HSL "Enable HSL libraries" OFF)
option(ENABLE_EXAMPLES "Build examples" OFF)
option(ENABLE_TORCH "Enable Torch interface" OFF)
option(MC__USE_FADBAD "Use the FADBAD++ forward-AD implementation" OFF)

set(MC_INTERVAL_LIBRARY
    "BOOST"
    CACHE STRING
          "Interval library backend: BOOST, PROFIL, FILIB, or NONVERIFIED")
set_property(CACHE MC_INTERVAL_LIBRARY PROPERTY STRINGS BOOST PROFIL FILIB
                                                NONVERIFIED)

# Dependencies

# Python 3 Configuration CUSTOM_PYTHON_PATH may be either: * an executable, e.g.
# /path/to/.venv/bin/python * a virtual-environment/root directory, e.g.
# /path/to/.venv Leave empty to let CMake find Python.
set(CUSTOM_PYTHON_PATH
    ""
    CACHE STRING "Custom Python root directory or Python executable")

if(CUSTOM_PYTHON_PATH)
  if(EXISTS "${CUSTOM_PYTHON_PATH}" AND NOT IS_DIRECTORY
                                        "${CUSTOM_PYTHON_PATH}")
    set(Python3_EXECUTABLE
        "${CUSTOM_PYTHON_PATH}"
        CACHE FILEPATH "Python executable" FORCE)
    message(
      STATUS
        "Custom Python path recognized as executable: ${Python3_EXECUTABLE}")
  elseif(IS_DIRECTORY "${CUSTOM_PYTHON_PATH}")
    set(Python3_ROOT_DIR
        "${CUSTOM_PYTHON_PATH}"
        CACHE PATH "Python root directory" FORCE)
    find_program(
      Python3_EXECUTABLE_HINT
      NAMES python3 python
      PATHS "${CUSTOM_PYTHON_PATH}/bin" "${CUSTOM_PYTHON_PATH}"
      NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
    if(Python3_EXECUTABLE_HINT)
      set(Python3_EXECUTABLE
          "${Python3_EXECUTABLE_HINT}"
          CACHE FILEPATH "Python executable" FORCE)
      message(
        STATUS "Found Python executable in custom root: ${Python3_EXECUTABLE}")
    else()
      message(
        FATAL_ERROR
          "CUSTOM_PYTHON_PATH='${CUSTOM_PYTHON_PATH}' is a directory, but no python/python3 executable was found in it"
      )
    endif()
  else()
    message(
      FATAL_ERROR "CUSTOM_PYTHON_PATH='${CUSTOM_PYTHON_PATH}' does not exist")
  endif()
endif()

# Find Python with the selected executable locked in. Development.Module is the
# appropriate component for building extension modules such as pymc; it avoids
# unnecessarily requesting the embedding library.
set(Python3_FIND_VIRTUALENV FIRST)
find_package(
  Python3
  COMPONENTS Interpreter Development.Module
  REQUIRED)

execute_process(
  COMMAND "${Python3_EXECUTABLE}" -c
          "import sysconfig; print(sysconfig.get_path('platlib') or '')"
  OUTPUT_VARIABLE Python3_SITEARCH
  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)

message(STATUS "--------------------------------------------------------")
message(STATUS "Python3 Interpreter: ${Python3_EXECUTABLE}")
message(STATUS "Python3 Module libs: ${Python3_LIBRARIES}")
message(STATUS "Python3 Includes:    ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3 site-packages: ${Python3_SITEARCH}")
message(STATUS "--------------------------------------------------------")

# pybind11 Check if extern/pybind11 is populated
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/CMakeLists.txt")
  message(
    FATAL_ERROR
      "pybind11 submodule not found. Please run 'git submodule update --init --recursive'"
  )
endif()
set(PYBIND11_FINDPYTHON ON)
add_subdirectory(extern/pybind11)

# Boost Always find Boost because it might be used for other things or by
# default
find_package(Boost CONFIG REQUIRED)

# Armadillo
find_package(Armadillo REQUIRED)

# LAPACK / BLAS
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Torch (Optional)
#
# Important: do not call find_package(Torch) here.  PyTorch's TorchConfig.cmake
# can enable CUDA and force CMake to identify nvcc even when MC++ only needs
# CPU-side libtorch symbols.  For MC++ we only link against torch_cpu and c10.
#
# TORCH_PYTHON_PREFIX may be: * the torch Python package directory, e.g.
# /path/to/.venv/lib/python3.14/site-packages/torch * the CMake prefix inside
# that package, e.g. /path/to/.../site-packages/torch/share/cmake * the
# TorchConfig directory, e.g. /path/to/.../site-packages/torch/share/cmake/Torch
# Leave empty to import torch from Python3_EXECUTABLE and discover the package.
set(TORCH_PYTHON_PREFIX
    ""
    CACHE
      PATH
      "Optional Torch Python package directory, CMake prefix, or TorchConfig directory"
)

if(ENABLE_TORCH)
  message(STATUS "Enabling Torch interface")

  set(MC_TORCH_PACKAGE_DIR "")
  set(MC_TORCH_LIB_DIR "")
  set(MC_TORCH_INCLUDE_DIRS "")
  set(MC_TORCH_COMPILE_DEFINITIONS "")

  if(TORCH_PYTHON_PREFIX)
    if(EXISTS "${TORCH_PYTHON_PREFIX}/include" AND EXISTS
                                                   "${TORCH_PYTHON_PREFIX}/lib")
      # User supplied the torch Python package directory.
      set(MC_TORCH_PACKAGE_DIR "${TORCH_PYTHON_PREFIX}")
    elseif(EXISTS "${TORCH_PYTHON_PREFIX}/TorchConfig.cmake")
      # User supplied .../share/cmake/Torch.
      get_filename_component(_MC_TORCH_CMAKE_DIR "${TORCH_PYTHON_PREFIX}"
                             DIRECTORY)
      get_filename_component(_MC_TORCH_SHARE_DIR "${_MC_TORCH_CMAKE_DIR}"
                             DIRECTORY)
      get_filename_component(MC_TORCH_PACKAGE_DIR "${_MC_TORCH_SHARE_DIR}"
                             DIRECTORY)
    elseif(EXISTS "${TORCH_PYTHON_PREFIX}/Torch/TorchConfig.cmake")
      # User supplied .../share/cmake.
      get_filename_component(_MC_TORCH_SHARE_DIR "${TORCH_PYTHON_PREFIX}"
                             DIRECTORY)
      get_filename_component(MC_TORCH_PACKAGE_DIR "${_MC_TORCH_SHARE_DIR}"
                             DIRECTORY)
    else()
      message(
        FATAL_ERROR
          "TORCH_PYTHON_PREFIX='${TORCH_PYTHON_PREFIX}' was provided, but it does not look like a Torch Python package directory, "
          "a Torch CMake prefix, or a TorchConfig directory.")
    endif()
  else()
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" -c
              "import os, torch; print(os.path.dirname(torch.__file__))"
      RESULT_VARIABLE MC_TORCH_DETECT_RESULT
      OUTPUT_VARIABLE MC_TORCH_PACKAGE_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_VARIABLE MC_TORCH_DETECT_ERROR)
    if(NOT MC_TORCH_DETECT_RESULT EQUAL 0 OR NOT MC_TORCH_PACKAGE_DIR)
      message(
        FATAL_ERROR
          "ENABLE_TORCH=ON, but torch could not be imported from '${Python3_EXECUTABLE}'.\n"
          "Install torch in that uv environment or pass -DTORCH_PYTHON_PREFIX=/path/to/site-packages/torch."
      )
    endif()
  endif()

  if(NOT EXISTS "${MC_TORCH_PACKAGE_DIR}/include"
     OR NOT EXISTS "${MC_TORCH_PACKAGE_DIR}/lib")
    message(
      FATAL_ERROR
        "Resolved Torch package directory '${MC_TORCH_PACKAGE_DIR}' does not contain both include/ and lib/ directories."
    )
  endif()

  set(MC_TORCH_LIB_DIR "${MC_TORCH_PACKAGE_DIR}/lib")

  list(APPEND MC_TORCH_INCLUDE_DIRS "${MC_TORCH_PACKAGE_DIR}/include"
       "${MC_TORCH_PACKAGE_DIR}/include/torch/csrc/api/include")

  foreach(_MC_TORCH_OPTIONAL_INCLUDE_DIR "${MC_TORCH_PACKAGE_DIR}/include/TH"
                                         "${MC_TORCH_PACKAGE_DIR}/include/THC")
    if(EXISTS "${_MC_TORCH_OPTIONAL_INCLUDE_DIR}")
      list(APPEND MC_TORCH_INCLUDE_DIRS "${_MC_TORCH_OPTIONAL_INCLUDE_DIR}")
    endif()
  endforeach()

  find_library(
    TORCH_CPU_LIBRARY
    NAMES torch_cpu
    PATHS "${MC_TORCH_LIB_DIR}"
    NO_DEFAULT_PATH)
  if(NOT TORCH_CPU_LIBRARY)
    message(FATAL_ERROR "Could not find libtorch_cpu in '${MC_TORCH_LIB_DIR}'")
  endif()

  find_library(
    C10_LIBRARY
    NAMES c10
    PATHS "${MC_TORCH_LIB_DIR}"
    NO_DEFAULT_PATH)
  if(NOT C10_LIBRARY)
    message(FATAL_ERROR "Could not find libc10 in '${MC_TORCH_LIB_DIR}'")
  endif()

  # Mirror PyTorch's libstdc++ ABI setting without using TorchConfig.cmake.
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c
            "import torch; print(int(torch._C._GLIBCXX_USE_CXX11_ABI))"
    RESULT_VARIABLE MC_TORCH_ABI_RESULT
    OUTPUT_VARIABLE MC_TORCH_CXX11_ABI
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(MC_TORCH_ABI_RESULT EQUAL 0 AND MC_TORCH_CXX11_ABI MATCHES "^[01]$")
    list(APPEND MC_TORCH_COMPILE_DEFINITIONS
         "_GLIBCXX_USE_CXX11_ABI=${MC_TORCH_CXX11_ABI}")
  else()
    message(
      WARNING
        "Could not determine PyTorch _GLIBCXX_USE_CXX11_ABI setting; using compiler default"
    )
  endif()

  message(STATUS "Torch package dir: ${MC_TORCH_PACKAGE_DIR}")
  message(STATUS "Torch include dirs: ${MC_TORCH_INCLUDE_DIRS}")
  message(STATUS "Torch lib dir:      ${MC_TORCH_LIB_DIR}")
  message(STATUS "Torch CPU library:  ${TORCH_CPU_LIBRARY}")
  message(STATUS "C10 library:        ${C10_LIBRARY}")
  message(STATUS "Torch ABI defs:     ${MC_TORCH_COMPILE_DEFINITIONS}")
endif()

# HSL (Optional)
if(ENABLE_HSL)
  message(STATUS "Enabling HSL support")
  find_library(LIB_MC13 NAMES mc13)
  find_library(LIB_MC21 NAMES mc21)
  find_library(LIB_MC33 NAMES mc33)

  set(HSL_LIBRARIES "")
  if(LIB_MC13)
    message(STATUS "  Found MC13: ${LIB_MC13}")
    list(APPEND HSL_LIBRARIES ${LIB_MC13})
  else()
    message(
      STATUS
        "  MC13 not found via find_library, assuming 'mc13' is in linker path")
    list(APPEND HSL_LIBRARIES mc13)
  endif()

  if(LIB_MC21)
    message(STATUS "  Found MC21: ${LIB_MC21}")
    list(APPEND HSL_LIBRARIES ${LIB_MC21})
  else()
    message(
      STATUS
        "  MC21 not found via find_library, assuming 'mc21' is in linker path")
    list(APPEND HSL_LIBRARIES mc21)
  endif()

  if(LIB_MC33)
    message(STATUS "  Found MC33: ${LIB_MC33}")
    list(APPEND HSL_LIBRARIES ${LIB_MC33})
  else()
    message(
      STATUS
        "  MC33 not found via find_library, assuming 'mc33' is in linker path")
    list(APPEND HSL_LIBRARIES mc33)
  endif()

  list(APPEND HSL_LIBRARIES gfortran)
  add_compile_definitions(MC__USE_HSL)
endif()

# Interval Library Configuration
message(STATUS "--------------------------------------------------------")
message(STATUS "Active Interval Library Backend: ${MC_INTERVAL_LIBRARY}")
message(STATUS "--------------------------------------------------------")

set(INTERVAL_LIBRARIES "")
set(INTERVAL_INCLUDE_DIRS "")
set(INTERVAL_DEFINITIONS "")
set(INTERVAL_COMPILE_OPTIONS "")

if(MC_INTERVAL_LIBRARY STREQUAL "BOOST")
  list(APPEND INTERVAL_DEFINITIONS MC__USE_BOOST)
  # Boost headers are already in Boost_INCLUDE_DIRS, which we add globally or to
  # target
elseif(MC_INTERVAL_LIBRARY STREQUAL "PROFIL")
  if(NOT DEFINED PROFIL_HOME)
    if(DEFINED ENV{PROFIL_HOME})
      set(PROFIL_HOME $ENV{PROFIL_HOME})
    else()
      message(
        FATAL_ERROR
          "PROFIL_HOME not set. Please set PROFIL_HOME environment variable or cmake variable."
      )
    endif()
  endif()
  list(APPEND INTERVAL_INCLUDE_DIRS ${PROFIL_HOME}/include)
  list(
    APPEND
    INTERVAL_LIBRARIES
    -L${PROFIL_HOME}/lib
    ProfilPackages
    Profil
    Bias
    lr)
  list(APPEND INTERVAL_DEFINITIONS MC__USE_PROFIL)
elseif(MC_INTERVAL_LIBRARY STREQUAL "FILIB")
  if(NOT DEFINED FILIB_HOME)
    if(DEFINED ENV{FILIB_HOME})
      set(FILIB_HOME $ENV{FILIB_HOME})
    else()
      message(
        FATAL_ERROR
          "FILIB_HOME not set. Please set FILIB_HOME environment variable or cmake variable."
      )
    endif()
  endif()
  list(APPEND INTERVAL_INCLUDE_DIRS ${FILIB_HOME}/include
       ${FILIB_HOME}/include/interval)
  list(APPEND INTERVAL_LIBRARIES -L${FILIB_HOME}/lib prim)
  list(APPEND INTERVAL_DEFINITIONS MC__USE_FILIB)
  list(APPEND INTERVAL_COMPILE_OPTIONS -frounding-math)
elseif(MC_INTERVAL_LIBRARY STREQUAL "NONVERIFIED")
  # No interval library
else()
  message(
    FATAL_ERROR
      "Invalid MC_INTERVAL_LIBRARY: ${MC_INTERVAL_LIBRARY}. Must be one of BOOST, PROFIL, FILIB, NONVERIFIED."
  )
endif()

# MC Library (Interface)
add_library(mc INTERFACE)

# Include directories
target_include_directories(
  mc
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/mc>
    $<INSTALL_INTERFACE:include>
    # ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/boost
    #${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/fadbad++
    # ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/cpplapack-2015.05.11-1/include
    ${Boost_INCLUDE_DIRS}
    ${ARMADILLO_INCLUDE_DIRS}
    ${INTERVAL_INCLUDE_DIRS})

# Compile definitions
target_compile_definitions(
  mc INTERFACE MC__USE_THREAD MC__USE_TADIFF MC__USE_ARMADILLO
               ${INTERVAL_DEFINITIONS})

if(MC__USE_FADBAD)
  target_compile_definitions(mc INTERFACE MC__USE_FADBAD)
endif()

# Match the classical Make warning flags and add backend-specific compile
# options such as -frounding-math only when required, e.g. with FILIB.
target_compile_options(
  mc INTERFACE -Wall -Wno-misleading-indentation -Wno-unknown-pragmas
               -Wno-parentheses -Wno-return-type ${INTERVAL_COMPILE_OPTIONS})

# Link libraries
target_link_libraries(mc INTERFACE ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES}
                                   ${LAPACK_LIBRARIES} ${INTERVAL_LIBRARIES})

if(ENABLE_HSL)
  target_link_libraries(mc INTERFACE ${HSL_LIBRARIES})
endif()

if(ENABLE_TORCH)
  target_include_directories(mc INTERFACE ${MC_TORCH_INCLUDE_DIRS})
  target_compile_definitions(mc INTERFACE MC__USE_TORCH
                                          ${MC_TORCH_COMPILE_DEFINITIONS})
  target_link_libraries(mc INTERFACE ${TORCH_CPU_LIBRARY} ${C10_LIBRARY})
endif()

# PyMC Library

# Explicitly list sources as per pymc.mk
set(PYMC_SOURCES
    src/pymc/mcfunc.cpp
    src/pymc/interval.cpp
    src/pymc/mccormick.cpp
    src/pymc/specbnd.cpp
    src/pymc/tmodel.cpp
    src/pymc/cmodel.cpp
    src/pymc/scmodel.cpp
    src/pymc/sicmodel.cpp
    src/pymc/supmodel.cpp
    src/pymc/polimage.cpp
    src/pymc/ellimage.cpp
    src/pymc/smon.cpp
    src/pymc/spoly.cpp
    src/pymc/ffdep.cpp
    src/pymc/ffinv.cpp
    src/pymc/ffunc.cpp
    src/pymc/ffmon.cpp
    src/pymc/ffpoly.cpp
    src/pymc/slift.cpp
    src/pymc/fflin.cpp
    src/pymc/ffmlp.cpp
    src/pymc/ffdagext.cpp
    src/pymc/ffvect.cpp
    src/pymc/ffcustom.cpp
    src/pymc/pymc.cpp)

pybind11_add_module(pymc ${PYMC_SOURCES})

# Keep the vendored pybind11 include directory ahead of Torch's bundled copy.
target_include_directories(
  pymc BEFORE PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/include")

target_link_libraries(pymc PRIVATE mc)

if(ENABLE_TORCH
   AND MC_TORCH_LIB_DIR
   AND EXISTS "${MC_TORCH_LIB_DIR}")
  set_property(
    TARGET pymc
    APPEND
    PROPERTY BUILD_RPATH "${MC_TORCH_LIB_DIR}")
  set_property(
    TARGET pymc
    APPEND
    PROPERTY INSTALL_RPATH "${MC_TORCH_LIB_DIR}")
endif()

# Installation

# Direct CMake installs default to: ${CMAKE_INSTALL_PREFIX}/include   for MC++
# headers ${CMAKE_INSTALL_PREFIX}/lib       for the compiled pymc extension
# module ${CMAKE_INSTALL_PREFIX}/notebook  for notebooks and Python scripts
#
# For wheel builds with scikit-build-core, override PYMC_INSTALL_DIR with
# -DPYMC_INSTALL_DIR=. and select only the python_modules component.
set(PYMC_INSTALL_DIR
    "lib"
    CACHE PATH "Installation directory for pymc module")
set(MCPP_NOTEBOOK_INSTALL_DIR
    "notebook"
    CACHE PATH "Installation directory for notebooks and Python scripts")

# Display the computed installation path for user confirmation
message(STATUS "--------------------------------------------------------")
message(
  STATUS "Installation Prefix (CMAKE_INSTALL_PREFIX): ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Headers will be installed to: ${CMAKE_INSTALL_PREFIX}/include")
message(
  STATUS
    "pymc module will be installed to: ${CMAKE_INSTALL_PREFIX}/${PYMC_INSTALL_DIR}"
)
message(
  STATUS
    "Notebooks/scripts will be installed to: ${CMAKE_INSTALL_PREFIX}/${MCPP_NOTEBOOK_INSTALL_DIR}"
)
message(STATUS "--------------------------------------------------------")

# Install headers Explicitly list headers as per mc.mk
set(MC_HEADERS
    mcfunc.hpp
    mctime.hpp
    # mclapack.hpp
    mcop.hpp
    mcboost.hpp
    mcprofil.hpp
    mcfilib.hpp
    fdiff.hpp
    bdiff.hpp
    tdiff.hpp
    mcfadbad.hpp
    ocbase.hpp
    interval.hpp
    mccormick.hpp
    specbnd.hpp
    supmodel.hpp
    pwcu.hpp
    pwlu.hpp
    ellipsoid.hpp
    ellimage.hpp
    polimage.hpp
    polymodel.hpp
    tmodel.hpp
    cmodel.hpp
    smon.hpp
    scmodel.hpp
    sicmodel.hpp
    ffdep.hpp
    ffinv.hpp
    ffunc.hpp
    ffexpr.hpp
    mchsl.hpp
    slift.hpp
    spoly.hpp
    ocbase.hpp
    fflin.hpp
    ffspol.hpp
    ffvect.hpp
    ffextern.hpp
    ffdagext.hpp
    ffmlp.hpp
    ffmlpreg.hpp
    ffcustom.hpp)

# Prefix headers with source path
set(MC_HEADERS_FULL "")
foreach(HEADER ${MC_HEADERS})
  list(APPEND MC_HEADERS_FULL "src/mc/${HEADER}")
endforeach()

install(
  FILES ${MC_HEADERS_FULL}
  DESTINATION include
  COMPONENT headers)

# Install pymc into ${CMAKE_INSTALL_PREFIX}/${PYMC_INSTALL_DIR}.
install(
  TARGETS pymc
  DESTINATION ${PYMC_INSTALL_DIR}
  COMPONENT python_modules)

# Install notebooks and Python scripts, if the notebook/ directory exists.
install(
  DIRECTORY notebook/
  DESTINATION ${MCPP_NOTEBOOK_INSTALL_DIR}
  COMPONENT notebooks
  OPTIONAL FILES_MATCHING
  PATTERN "*.ipynb"
  PATTERN "*.py"
  PATTERN "*.md"
  PATTERN "*.pt"
  PATTERN "*.toml"
  PATTERN "*.lock"
  PATTERN "*.python-version")

# Examples
if(ENABLE_EXAMPLES)
  file(
    GLOB_RECURSE EXAMPLE_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/test/*.cpp"
  )

  foreach(SOURCE_FILE IN LISTS EXAMPLE_SOURCES)

    # FADBAD-specific examples require MC__USE_FADBAD.
    if(SOURCE_FILE MATCHES "/test/FADBAD/" AND NOT MC__USE_FADBAD)
      continue()
    endif()

    get_filename_component(EXEC_NAME "${SOURCE_FILE}" NAME_WE)
    add_executable("${EXEC_NAME}" "${SOURCE_FILE}")
    target_link_libraries("${EXEC_NAME}" PRIVATE mc)
  endforeach()
endif()
