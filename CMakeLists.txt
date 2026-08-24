cmake_minimum_required(VERSION 3.16)
project(CRAPPYSOLVER LANGUAGES CXX)

# ---- Compiler / flags (match your g++ line) ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to a build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# -O0 -g -w
add_compile_options(-w -g -O0)

# ---- Executable name ----
set(EXECUTABLE_NAME CrappySolver)   # <-- change this to whatever you want the binary called

# ---- Paths (edit if needed) ----
set(PROJ_ROOT ${CMAKE_CURRENT_SOURCE_DIR})  # e.g., /Users/username/Desktop/CrappySolver

# ---- Sources (test.cpp + BB/*.cpp) ----
file(GLOB BB_SOURCES CONFIGURE_DEPENDS "${PROJ_ROOT}/BB/*.cpp")

add_executable(${EXECUTABLE_NAME}
    "${PROJ_ROOT}/test.cpp"
    ${BB_SOURCES}
)

# ---- Include directories (match your -I list) ----
target_include_directories(${EXECUTABLE_NAME} PRIVATE
    "${PROJ_ROOT}"
    "${PROJ_ROOT}/src/3rdparty/fadbad++"
    "${PROJ_ROOT}/include"
    "${PROJ_ROOT}/src/3rdparty/cpplapack-2015.05.11-1/include"

)

# ---- Library search paths / linking (edit if needed) ----
# target_link_directories(${EXECUTABLE_NAME} PRIVATE "${PROJ_ROOT}/some/lib/path")
# target_link_libraries(${EXECUTABLE_NAME} PRIVATE some_lib_name)

# Optional: if your cluster needs it for static_pic objects
# target_link_options(${EXECUTABLE_NAME} PRIVATE -Wl,--no-as-needed)