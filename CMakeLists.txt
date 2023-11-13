# =========================================================
# Numerical Environment Basics

cmake_minimum_required(VERSION 3.20)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "cmake")

set(ProjectName "nuenv-basics")
set(ProjectNamespace "nuenv")

project("${ProjectName}")

set(CMAKE_CXX_STANDARD 23)

add_compile_options(-Wall -Werror -Wno-error=uninitialized
        $<$<CONFIG:RELEASE>:-O3>
        $<$<CONFIG:RELEASE>:-DNDEBUG>
        $<$<CONFIG:RELWITHDEBINFO>:-O2>
        $<$<CONFIG:RELWITHDEBINFO>:-DNDEBUG>
        $<$<CONFIG:RELWITHDEBINFO>:-g>
        $<$<CONFIG:MINSIZEREL>:-Oz>
        $<$<CONFIG:MINSIZEREL>:-DNDEBUG>
        $<$<CONFIG:DEBUG>:-Og>
        $<$<CONFIG:DEBUG>:-g>
)

# =========================================================
# gtest

include(FetchContent)
FetchContent_Declare(googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.13.0.zip)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set_property(TARGET gmock PROPERTY FOLDER "gtest")
set_property(TARGET gmock_main PROPERTY FOLDER "gtest")
set_property(TARGET gtest PROPERTY FOLDER "gtest")
set_property(TARGET gtest_main PROPERTY FOLDER "gtest")

# =========================================================
# Log4cxx

find_package(log4cxx 0.11)

# =========================================================
# Eigen

include(FetchContent)
FetchContent_Declare(Eigen
        URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip)

set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(EIGEN_FORCE_SHARED_CRT ON CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
FetchContent_MakeAvailable(Eigen)

# =========================================================
# Targets

# --- project

add_library(${ProjectName} INTERFACE)
add_library(${ProjectName}::${ProjectNamespace} ALIAS ${ProjectName})

target_include_directories(${ProjectName} INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>)

#target_include_directories("${ProjectName}" PRIVATE $<TARGET_PROPERTY:log4cxx,INTERFACE_INCLUDE_DIRECTORIES>)
#target_link_libraries( "${ProjectName}" PRIVATE log4cxx)
target_link_libraries("${ProjectName}" INTERFACE Eigen3::Eigen)

set_target_properties("${ProjectName}" PROPERTIES LINKER_LANGUAGE CXX)

# --- main

add_executable("${ProjectName}-main" examples/main.cpp)
target_link_libraries("${ProjectName}-main" PUBLIC "${ProjectName}")

# --- gtest

enable_testing()

add_executable("${ProjectName}-test"
        test/Algorithm/Search.cpp
        test/Algorithm/Space.cpp
        test/Integrate/Quadrature.cpp
        test/Integrate/RK4.cpp
        test/Interpolate/Interp1d.cpp
        test/Optimize/DiffEvolution.cpp
)

target_link_libraries("${ProjectName}-test" gtest_main "${ProjectName}")

include(GoogleTest)
gtest_discover_tests("${ProjectName}-test")

# =========================================================
