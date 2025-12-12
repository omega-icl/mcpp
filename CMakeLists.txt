cmake_minimum_required(VERSION 3.14)
project(mcpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 1. Compilation Modes
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# Set specific flags if needed, although standard CMake flags are usually sufficient.
# Debug: -g is standard. -O0 is often default for Debug but we can ensure it.
# Release: -O3 is standard. User asked for high optimization.
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")


# Options
option(ENABLE_HSL "Enable HSL libraries" OFF)
option(ENABLE_EXAMPLES "Build examples" OFF)

set(MC_INTERVAL_LIBRARY "BOOST" CACHE STRING "Interval library backend: BOOST, PROFIL, FILIB, or NONVERIFIED")
set_property(CACHE MC_INTERVAL_LIBRARY PROPERTY STRINGS BOOST PROFIL FILIB NONVERIFIED)

# Dependencies

# Python 3
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# pybind11
# Check if extern/pybind11 is populated
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/CMakeLists.txt")
    message(FATAL_ERROR "pybind11 submodule not found. Please run 'git submodule update --init --recursive'")
endif()
add_subdirectory(extern/pybind11)

# Boost
# Always find Boost because it might be used for other things or by default
find_package(Boost REQUIRED)

# Armadillo
find_package(Armadillo REQUIRED)

# LAPACK / BLAS
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

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
set(INTERVAL_LIBRARIES "")
set(INTERVAL_INCLUDE_DIRS "")
set(INTERVAL_DEFINITIONS "")

# Always include BOOST_UBLAS_NO_STD_CERR
list(APPEND INTERVAL_DEFINITIONS BOOST_UBLAS_NO_STD_CERR)

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
    ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/fadbad++
    ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/cpplapack-2015.05.11-1/include
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

# PyMC Library

# Explicitly list sources as per pymc.mk
set(PYMC_SOURCES
    src/pymc/mcfunc.cpp
    src/pymc/interval.cpp
    src/pymc/mccormick.cpp
    src/pymc/supmodel.cpp
    src/pymc/ffunc.cpp
    src/pymc/fflin.cpp
    src/pymc/ffvect.cpp
    src/pymc/ffcustom.cpp
    src/pymc/pymc.cpp
)

pybind11_add_module(pymc ${PYMC_SOURCES})

target_link_libraries(pymc PRIVATE mc)

# Installation

# Install headers
# Explicitly list headers as per mc.mk
set(MC_HEADERS
    mcfunc.hpp mctime.hpp mclapack.hpp
    mcop.hpp mcboost.hpp mcprofil.hpp mcfilib.hpp mcfadbad.hpp
    interval.hpp mccormick.hpp specbnd.hpp
    supmodel.hpp pwcu.hpp pwlu.hpp
    ellipsoid.hpp ellimage.hpp polimage.hpp
    polymodel.hpp tmodel.hpp cmodel.hpp smon.hpp scmodel.hpp sicmodel.hpp
    spoly.hpp squad.hpp
    ffdep.hpp ffinv.hpp ffunc.hpp ffexpr.hpp mchsl.hpp sred.hpp slift.hpp selim.hpp
    fflin.hpp ffspol.hpp ffvect.hpp ffextern.hpp ffdagext.hpp ffmlp.hpp ffmlpreg.hpp ffcustom.hpp
)

# Prefix headers with source path
set(MC_HEADERS_FULL "")
foreach(HEADER ${MC_HEADERS})
    list(APPEND MC_HEADERS_FULL "src/mc/${HEADER}")
endforeach()

install(FILES ${MC_HEADERS_FULL} DESTINATION include COMPONENT headers)

# Install pymc
set(PYMC_INSTALL_DIR "lib" CACHE STRING "Installation directory for pymc module")
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
