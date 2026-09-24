# SPDX-License-Identifier: MIT
cmake_minimum_required(VERSION 3.22...4.2)
project(Maxent VERSION 2.0.0 LANGUAGES CXX)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are not supported. Use a separate build directory.")
endif()

include(GNUInstallDirs)
# Installed programs keep the directories of their shared dependencies
# (ALPSCore, HDF5, Boost) in their RPATH.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

option(MAXENT_BUILD_TESTS "Build unit tests and the regression suite" ON)
option(MAXENT_BUILD_UTILITIES "Build the kk and legendre_convert utilities" ON)
option(MAXENT_USE_LAPACK "Use LAPACK for the SVD of the kernel instead of Eigen" OFF)
option(MAXENT_WERROR "Treat compiler warnings as errors" OFF)
option(MAXENT_REGRESSION_FULL "Also register the full-size regression cases (about 1 min)" OFF)
option(MAXENT_USE_SYSTEM_GTEST "Use an installed GoogleTest instead of downloading it" OFF)

# ---------------------------------------------------------------- dependencies
# ALPSCore first: it reports the Boost version it was built with, and Maxent
# must use the same one.
find_package(ALPSCore 2.0 REQUIRED)
if(DEFINED ALPSCore_BOOST_VERSION)
  find_package(Boost ${ALPSCore_BOOST_VERSION} EXACT CONFIG)
else()
  find_package(Boost 1.70 CONFIG)
endif()
if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost ${ALPSCore_BOOST_VERSION} (the version ALPSCore was built with) was not "
    "found. Pass -DBoost_DIR=<dir containing BoostConfig.cmake> or add its prefix to CMAKE_PREFIX_PATH.")
endif()
find_package(Eigen3 3.3 CONFIG REQUIRED)

set(HAVE_BLAS 0)
set(HAVE_LAPACK 0)
if(MAXENT_USE_LAPACK)
  find_package(LAPACK REQUIRED)  # also finds BLAS
  set(HAVE_BLAS 1)
  set(HAVE_LAPACK 1)
  message(STATUS "Using LAPACK for the SVD of the kernel")
endif()
configure_file(cmake/maxent_config.hpp.in ${PROJECT_BINARY_DIR}/config/maxent_config.hpp)

# ---------------------------------------------------------------- warnings
add_library(maxent_warnings INTERFACE)
target_compile_options(maxent_warnings INTERFACE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
  $<$<AND:$<BOOL:${MAXENT_WERROR}>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-Werror>)

# ---------------------------------------------------------------- library and program
add_library(maxent_core STATIC
  src/default_model.cpp
  src/maxent_backcont.cpp
  src/maxent_grid.cpp
  src/maxent_helper.cpp
  src/maxent_kernel.cpp
  src/maxent_params.cpp
  src/maxent_simulation.cpp)
add_library(maxent::core ALIAS maxent_core)
target_include_directories(maxent_core PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/config>)
target_link_libraries(maxent_core
  PUBLIC ${ALPSCore_LIBRARIES} Boost::headers Eigen3::Eigen
  PRIVATE maxent_warnings)
if(MAXENT_USE_LAPACK)
  target_link_libraries(maxent_core PUBLIC LAPACK::LAPACK)
endif()

add_executable(maxent src/maxent.cpp)
target_link_libraries(maxent PRIVATE maxent::core maxent_warnings)

# ---------------------------------------------------------------- tests and utilities
if(MAXENT_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

if(MAXENT_BUILD_UTILITIES)
  add_subdirectory(kk)
  add_subdirectory(legendre_convert)
endif()

# ---------------------------------------------------------------- install
install(TARGETS maxent RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
