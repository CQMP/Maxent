cmake_minimum_required(VERSION 2.8.12)
project (MAXENT)

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
find_package(ALPSCore REQUIRED)

#note: LAPACK also finds BLAS
find_package(LAPACK)
if(LAPACK_FOUND AND USE_LAPACK)
    set(HAVE_BLAS 1)
    set(HAVE_LAPACK 1)
    message(STATUS "Using LAPACK and BLAS routines for SVD")
else()
    set(HAVE_BLAS 0)
    set(HAVE_LAPACK 0)
endif()
configure_file(
	"${PROJECT_SOURCE_DIR}/cmake/maxent_config.hpp.in"
	"${PROJECT_BINARY_DIR}/config/maxent_config.hpp"
    )

include_directories("${PROJECT_BINARY_DIR}/config")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
message(STATUS "Finding packages, please stand by...")
include(UseGtest)
find_package(GSL REQUIRED)
find_package (Eigen3 3.1 REQUIRED)

set (GTEST_INCLUDE_DIR ${GTEST_ROOT}/include )
set (GTEST_SRC_DIR ${GTEST_ROOT}/src)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG -Wno-return-type-c-linkage")

#let gcc take advantage of Eigen3 vectorization
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
endif()

include_directories(${CMAKE_SOURCE_DIR}/gtest/include)
include_directories(${GTEST_INCLUDE_DIR})
include_directories(${GTEST_ROOT})
include_directories(${GSL_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(maxent 
maxent.cpp
maxent_helper.cpp
maxent_parms.cpp
maxent_parmsTest.cpp
maxent_grid.cpp
maxent_gridTest.cpp
maxent_kernel.cpp
maxent_legendre_util.cpp
maxent_simulation.cpp
maxent_simulationTest.cpp
maxent_backcont.cpp
maxent_backcontTest.cpp
default_modelTest.cpp
default_model.cpp
${GTEST_SRC_DIR}/gtest-all.cc
)
target_link_libraries(maxent
	${ALPSCore_LIBRARIES} ${GSL_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )

#add companion utilities
add_subdirectory(legendre_convert)
add_subdirectory(kk) 
add_subdirectory(pade/pade_arbitrary_degree)

