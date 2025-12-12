cmake_minimum_required(VERSION 3.14)
project(mcpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Options
option(ENABLE_HSL "Enable HSL libraries" OFF)

set(MC_INTERVAL_LIBRARY "BOOST" CACHE STRING "Interval library backend: BOOST, PROFIL, FILIB, or NONE")
set_property(CACHE MC_INTERVAL_LIBRARY PROPERTY STRINGS BOOST PROFIL FILIB NONE)

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
    # Assuming system libraries for HSL
    set(HSL_LIBRARIES mc13 mc21 mc33 gfortran)
    add_compile_definitions(MC__USE_HSL)
endif()

# Interval Library Configuration
set(INTERVAL_LIBRARIES "")
set(INTERVAL_INCLUDE_DIRS "")
set(INTERVAL_DEFINITIONS "")

if(MC_INTERVAL_LIBRARY STREQUAL "BOOST")
    list(APPEND INTERVAL_DEFINITIONS MC__USE_BOOST BOOST_UBLAS_NO_STD_CERR)
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
elseif(MC_INTERVAL_LIBRARY STREQUAL "NONE")
    # No interval library
else()
    message(FATAL_ERROR "Invalid MC_INTERVAL_LIBRARY: ${MC_INTERVAL_LIBRARY}. Must be one of BOOST, PROFIL, FILIB, NONE.")
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
file(GLOB PYMC_SOURCES src/pymc/*.cpp)
pybind11_add_module(pymc ${PYMC_SOURCES})

target_link_libraries(pymc PRIVATE mc)

# Installation

# Install headers
# The legacy makefile installs specific headers from src/mc to include/
# We will install all headers in src/mc to include/
install(DIRECTORY src/mc/ DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)

# Install pymc
install(TARGETS pymc DESTINATION lib)
