# SPDX-License-Identifier: MIT
#
# Stand-alone build of dump_components against an existing Maxent build tree.
# (It becomes a regular target of the main build in step 2.1.)
#
#   cmake -S test/regression/components -B <dir> \
#     -DMAXENT_SOURCE_DIR=<maxent source> -DMAXENT_BUILD_DIR=<maxent build> \
#     -DALPSCore_DIR=... -DBoost_DIR=...
cmake_minimum_required(VERSION 3.18)
project(maxent_dump_components CXX)

find_package(Boost CONFIG REQUIRED)
find_package(ALPSCore 2.0 REQUIRED)
find_package(Eigen3 REQUIRED CONFIG)
find_library(GSL_LIB gsl REQUIRED)
find_library(GSLCBLAS_LIB gslcblas REQUIRED)
find_package(LAPACK)  # needed when the Maxent build used USE_LAPACK
find_library(MAXENT_LIB NAMES libmaxent.a libmaxent.dylib libmaxent.so
             PATHS ${MAXENT_BUILD_DIR} NO_DEFAULT_PATH REQUIRED)

add_executable(dump_components dump_components.cpp)
target_include_directories(dump_components PRIVATE ${MAXENT_SOURCE_DIR}/src ${MAXENT_BUILD_DIR}/config)
target_link_libraries(dump_components PRIVATE ${MAXENT_LIB} ${ALPSCore_LIBRARIES} Eigen3::Eigen
                      ${GSL_LIB} ${GSLCBLAS_LIB} ${LAPACK_LIBRARIES})
