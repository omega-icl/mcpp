cmake_minimum_required(VERSION 3.18)
project(mcpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 1. Compilation Modes
message(STATUS "--------------------------------------------------------")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Setting build type to 'Release'")
else()
  message(STATUS "Setting build type to 'Debug'")
endif()
message(STATUS "--------------------------------------------------------")

# Set specific flags if needed, although standard CMake flags are usually sufficient.
# Debug: -g is standard. -O0 is often default for Debug but we can ensure it.
# Release: -O2 is standard. User asked for high optimization.
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG")

# Options
option(ENABLE_HSL "Enable HSL libraries" OFF)
option(ENABLE_EXAMPLES "Build examples" OFF)
option(ENABLE_TORCH "Enable Torch interface" OFF)

set(MC_INTERVAL_LIBRARY "BOOST" CACHE STRING "Interval library backend: BOOST, PROFIL, FILIB, or NONVERIFIED")
set_property(CACHE MC_INTERVAL_LIBRARY PROPERTY STRINGS BOOST PROFIL FILIB NONVERIFIED)

# Dependencies

# Python 3 Configuration
# CUSTOM_PYTHON_PATH may be either:
#   * an executable, e.g. /path/to/.venv/bin/python
#   * a virtual-environment/root directory, e.g. /path/to/.venv
# Leave empty to let CMake find Python.
set(CUSTOM_PYTHON_PATH "" CACHE STRING "Custom Python root directory or Python executable")

if(CUSTOM_PYTHON_PATH)
    if(EXISTS "${CUSTOM_PYTHON_PATH}" AND NOT IS_DIRECTORY "${CUSTOM_PYTHON_PATH}")
        set(Python3_EXECUTABLE "${CUSTOM_PYTHON_PATH}" CACHE FILEPATH "Python executable" FORCE)
        message(STATUS "Custom Python path recognized as executable: ${Python3_EXECUTABLE}")
    elseif(IS_DIRECTORY "${CUSTOM_PYTHON_PATH}")
        set(Python3_ROOT_DIR "${CUSTOM_PYTHON_PATH}" CACHE PATH "Python root directory" FORCE)
        find_program(Python3_EXECUTABLE_HINT
            NAMES python3 python
            PATHS "${CUSTOM_PYTHON_PATH}/bin" "${CUSTOM_PYTHON_PATH}"
            NO_DEFAULT_PATH
            NO_CMAKE_ENVIRONMENT_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
        )
        if(Python3_EXECUTABLE_HINT)
            set(Python3_EXECUTABLE "${Python3_EXECUTABLE_HINT}" CACHE FILEPATH "Python executable" FORCE)
            message(STATUS "Found Python executable in custom root: ${Python3_EXECUTABLE}")
        else()
            message(FATAL_ERROR "CUSTOM_PYTHON_PATH='${CUSTOM_PYTHON_PATH}' is a directory, but no python/python3 executable was found in it")
        endif()
    else()
        message(FATAL_ERROR "CUSTOM_PYTHON_PATH='${CUSTOM_PYTHON_PATH}' does not exist")
    endif()
endif()

# Find Python with the selected executable locked in.
# Development.Module is the appropriate component for building extension modules
# such as pymc; it avoids unnecessarily requesting the embedding library.
set(Python3_FIND_VIRTUALENV FIRST)
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)

execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('platlib') or '')"
    OUTPUT_VARIABLE Python3_SITEARCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

message(STATUS "--------------------------------------------------------")
message(STATUS "Python3 Interpreter: ${Python3_EXECUTABLE}")
message(STATUS "Python3 Module libs: ${Python3_LIBRARIES}")
message(STATUS "Python3 Includes:    ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3 site-packages: ${Python3_SITEARCH}")
message(STATUS "--------------------------------------------------------")

# pybind11
# Check if extern/pybind11 is populated
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/CMakeLists.txt")
    message(FATAL_ERROR "pybind11 submodule not found. Please run 'git submodule update --init --recursive'")
endif()
set(PYBIND11_FINDPYTHON ON)
add_subdirectory(extern/pybind11)

# Boost
# Always find Boost because it might be used for other things or by default
find_package(Boost REQUIRED)

# Armadillo
find_package(Armadillo REQUIRED)

# LAPACK / BLAS
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Torch (Optional)
# TORCH_PYTHON_PREFIX may be any of the following:
#   * the Python package directory, e.g. /path/to/.venv/lib/python3.14/site-packages/torch
#   * the CMake prefix returned by torch.utils.cmake_prefix_path
#   * the TorchConfig directory, e.g. .../share/cmake/Torch
set(TORCH_PYTHON_PREFIX "" CACHE PATH "Optional Torch Python package directory, CMake prefix, or TorchConfig directory")

if(ENABLE_TORCH)
    message(STATUS "Enabling Torch interface")

    set(MC_TORCH_CMAKE_PREFIX "")
    set(MC_TORCH_DIR "")
    set(MC_TORCH_PYTHON_LIB_DIR "")

    if(TORCH_PYTHON_PREFIX)
        if(EXISTS "${TORCH_PYTHON_PREFIX}/TorchConfig.cmake")
            # User supplied .../share/cmake/Torch directly.
            set(MC_TORCH_DIR "${TORCH_PYTHON_PREFIX}")
        elseif(EXISTS "${TORCH_PYTHON_PREFIX}/Torch/TorchConfig.cmake")
            # User supplied .../share/cmake.
            set(MC_TORCH_CMAKE_PREFIX "${TORCH_PYTHON_PREFIX}")
        elseif(EXISTS "${TORCH_PYTHON_PREFIX}/share/cmake/Torch/TorchConfig.cmake")
            # User supplied the torch Python package directory.
            set(MC_TORCH_CMAKE_PREFIX "${TORCH_PYTHON_PREFIX}/share/cmake")
            if(EXISTS "${TORCH_PYTHON_PREFIX}/lib")
                set(MC_TORCH_PYTHON_LIB_DIR "${TORCH_PYTHON_PREFIX}/lib")
            endif()
        else()
            message(FATAL_ERROR
                "TORCH_PYTHON_PREFIX='${TORCH_PYTHON_PREFIX}' was provided, but no TorchConfig.cmake was found. "
                "Use the torch package directory, torch.utils.cmake_prefix_path, or .../share/cmake/Torch.")
        endif()
    else()
        execute_process(
            COMMAND "${Python3_EXECUTABLE}" -c "import torch; print(torch.utils.cmake_prefix_path)"
            RESULT_VARIABLE MC_TORCH_DETECT_RESULT
            OUTPUT_VARIABLE MC_TORCH_CMAKE_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_VARIABLE MC_TORCH_DETECT_ERROR
        )
        if(NOT MC_TORCH_DETECT_RESULT EQUAL 0 OR NOT MC_TORCH_CMAKE_PREFIX)
            message(FATAL_ERROR
                "ENABLE_TORCH=ON, but Torch could not be imported from '${Python3_EXECUTABLE}'.\n"
                "Either install torch in that uv environment, or pass -DTORCH_PYTHON_PREFIX=/path/to/site-packages/torch.")
        endif()

        execute_process(
            COMMAND "${Python3_EXECUTABLE}" -c "import os, torch; print(os.path.join(os.path.dirname(torch.__file__), 'lib'))"
            OUTPUT_VARIABLE MC_TORCH_PYTHON_LIB_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
    endif()

    if(MC_TORCH_DIR)
        set(Torch_DIR "${MC_TORCH_DIR}" CACHE PATH "Directory containing TorchConfig.cmake" FORCE)
        message(STATUS "Using explicit Torch_DIR: ${Torch_DIR}")
    endif()

    if(MC_TORCH_CMAKE_PREFIX)
        list(PREPEND CMAKE_PREFIX_PATH "${MC_TORCH_CMAKE_PREFIX}")
        message(STATUS "Using Torch CMake prefix: ${MC_TORCH_CMAKE_PREFIX}")
    endif()

    find_package(Torch REQUIRED)

    set(TORCH_CXX_FLAGS_LIST "")
    if(TORCH_CXX_FLAGS)
        separate_arguments(TORCH_CXX_FLAGS_LIST NATIVE_COMMAND "${TORCH_CXX_FLAGS}")
    endif()

    message(STATUS "Torch include dirs: ${TORCH_INCLUDE_DIRS}")
    message(STATUS "Torch libraries:    ${TORCH_LIBRARIES}")
    message(STATUS "Torch CXX flags:    ${TORCH_CXX_FLAGS}")
    message(STATUS "Torch Python lib:   ${MC_TORCH_PYTHON_LIB_DIR}")
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
        message(STATUS "  MC13 not found via find_library, assuming 'mc13' is in linker path")
        list(APPEND HSL_LIBRARIES mc13)
    endif()

    if(LIB_MC21)
        message(STATUS "  Found MC21: ${LIB_MC21}")
        list(APPEND HSL_LIBRARIES ${LIB_MC21})
    else()
        message(STATUS "  MC21 not found via find_library, assuming 'mc21' is in linker path")
        list(APPEND HSL_LIBRARIES mc21)
    endif()

    if(LIB_MC33)
        message(STATUS "  Found MC33: ${LIB_MC33}")
        list(APPEND HSL_LIBRARIES ${LIB_MC33})
    else()
        message(STATUS "  MC33 not found via find_library, assuming 'mc33' is in linker path")
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

if(MC_INTERVAL_LIBRARY STREQUAL "BOOST")
    list(APPEND INTERVAL_DEFINITIONS MC__USE_BOOST)
    # Boost headers are already in Boost_INCLUDE_DIRS, which we add globally or to target
elseif(MC_INTERVAL_LIBRARY STREQUAL "PROFIL")
    if(NOT DEFINED PROFIL_HOME)
        if(DEFINED ENV{PROFIL_HOME})
            set(PROFIL_HOME $ENV{PROFIL_HOME})
        else()
             message(FATAL_ERROR "PROFIL_HOME not set. Please set PROFIL_HOME environment variable or cmake variable.")
        endif()
    endif()
    list(APPEND INTERVAL_INCLUDE_DIRS ${PROFIL_HOME}/include)
    list(APPEND INTERVAL_LIBRARIES -L${PROFIL_HOME}/lib ProfilPackages Profil Bias lr)
    list(APPEND INTERVAL_DEFINITIONS MC__USE_PROFIL)
elseif(MC_INTERVAL_LIBRARY STREQUAL "FILIB")
    if(NOT DEFINED FILIB_HOME)
         if(DEFINED ENV{FILIB_HOME})
            set(FILIB_HOME $ENV{FILIB_HOME})
        else()
             message(FATAL_ERROR "FILIB_HOME not set. Please set FILIB_HOME environment variable or cmake variable.")
        endif()
    endif()
    list(APPEND INTERVAL_INCLUDE_DIRS ${FILIB_HOME}/include ${FILIB_HOME}/include/interval)
    list(APPEND INTERVAL_LIBRARIES -L${FILIB_HOME}/lib prim)
    list(APPEND INTERVAL_DEFINITIONS MC__USE_FILIB)
    add_compile_options(-frounding-math)
elseif(MC_INTERVAL_LIBRARY STREQUAL "NONVERIFIED")
    # No interval library
else()
    message(FATAL_ERROR "Invalid MC_INTERVAL_LIBRARY: ${MC_INTERVAL_LIBRARY}. Must be one of BOOST, PROFIL, FILIB, NONVERIFIED.")
endif()


# MC Library (Interface)
add_library(mc INTERFACE)

# Include directories
target_include_directories(mc INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/mc>
    $<INSTALL_INTERFACE:include>
    #${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/boost
    ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/fadbad++
    #${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/cpplapack-2015.05.11-1/include
    ${Boost_INCLUDE_DIRS}
    ${ARMADILLO_INCLUDE_DIRS}
    ${INTERVAL_INCLUDE_DIRS}
)

# Compile definitions
target_compile_definitions(mc INTERFACE
    MC__USE_THREAD
    MC__USE_TADIFF
    MC__USE_ARMADILLO
    ${INTERVAL_DEFINITIONS}
)

# Link libraries
target_link_libraries(mc INTERFACE
    ${ARMADILLO_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${INTERVAL_LIBRARIES}
)

if(ENABLE_HSL)
    target_link_libraries(mc INTERFACE ${HSL_LIBRARIES})
endif()

if(ENABLE_TORCH)
    target_include_directories(mc INTERFACE ${TORCH_INCLUDE_DIRS})
    target_compile_definitions(mc INTERFACE MC__USE_TORCH)
    target_compile_options(mc INTERFACE -frounding-math ${TORCH_CXX_FLAGS_LIST})
    target_link_libraries(mc INTERFACE ${TORCH_LIBRARIES})
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
    src/pymc/pymc.cpp
)

pybind11_add_module(pymc ${PYMC_SOURCES})

# Keep the vendored pybind11 include directory ahead of Torch's bundled copy.
target_include_directories(pymc BEFORE PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/include")
target_link_libraries(pymc PRIVATE mc)

if(ENABLE_TORCH AND MC_TORCH_PYTHON_LIB_DIR AND EXISTS "${MC_TORCH_PYTHON_LIB_DIR}")
    set_property(TARGET pymc APPEND PROPERTY BUILD_RPATH "${MC_TORCH_PYTHON_LIB_DIR}")
    set_property(TARGET pymc APPEND PROPERTY INSTALL_RPATH "${MC_TORCH_PYTHON_LIB_DIR}")
endif()

# Installation

# Display the computed installation path for user confirmation
message(STATUS "--------------------------------------------------------")
message(STATUS "Installation Prefix (CMAKE_INSTALL_PREFIX): ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Headers will be installed to: ${CMAKE_INSTALL_PREFIX}/include")
message(STATUS "Libraries will be installed to: ${CMAKE_INSTALL_PREFIX}/lib")
message(STATUS "--------------------------------------------------------")

# Install headers
# Explicitly list headers as per mc.mk
set(MC_HEADERS
    mcfunc.hpp
    mctime.hpp
    #mclapack.hpp
    mcop.hpp
    mcboost.hpp
    mcprofil.hpp
    mcfilib.hpp
    mcfadbad.hpp
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
    #ocbase.hpp
    #ocenv.hpp
    #ffipdae.hpp
    fflin.hpp
    ffspol.hpp
    ffvect.hpp
    ffextern.hpp
    ffdagext.hpp
    ffmlp.hpp
    ffmlpreg.hpp
    ffcustom.hpp
)

# Prefix headers with source path
set(MC_HEADERS_FULL "")
foreach(HEADER ${MC_HEADERS})
    list(APPEND MC_HEADERS_FULL "src/mc/${HEADER}")
endforeach()

install(FILES ${MC_HEADERS_FULL} DESTINATION include COMPONENT headers)

# Install pymc.  For uv builds, default to the selected environment's site-packages
# directory when it is available; otherwise fall back to ${CMAKE_INSTALL_PREFIX}/lib.
if(Python3_SITEARCH)
    set(PYMC_INSTALL_DIR "${Python3_SITEARCH}" CACHE PATH "Installation directory for pymc module")
else()
    set(PYMC_INSTALL_DIR "lib" CACHE PATH "Installation directory for pymc module")
endif()
install(TARGETS pymc DESTINATION ${PYMC_INSTALL_DIR} COMPONENT python_modules)

# Examples
if(ENABLE_EXAMPLES)
    file(GLOB_RECURSE EXAMPLE_SOURCES "test/*.cpp")
    foreach(SOURCE_FILE ${EXAMPLE_SOURCES})
        get_filename_component(EXEC_NAME ${SOURCE_FILE} NAME_WE)
        add_executable(${EXEC_NAME} ${SOURCE_FILE})
        target_link_libraries(${EXEC_NAME} PRIVATE mc)
    endforeach()
endif()
