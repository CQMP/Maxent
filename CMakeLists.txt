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

set(LIB_FILES 
  ./src/maxent_helper.cpp
  ./src/maxent_parms.cpp
  ./src/maxent_grid.cpp
  ./src/maxent_kernel.cpp
  ./src/maxent_legendre_util.cpp
  ./src/maxent_simulation.cpp
  ./src/maxent_backcont.cpp
  ./src/default_model.cpp
  )

ADD_LIBRARY(libmaxent ${LIB_FILES})
target_link_libraries(libmaxent
	${ALPSCore_LIBRARIES} ${GSL_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )

#executable
add_executable(maxent 
  ./src/maxent.cpp
  )
target_link_libraries(maxent libmaxent ${ALPSCore_LIBRARIES} ${GSL_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

#testing setup
option(Testing "Enable testing" ON)
include(EnableGtests) #defined in ./cmake
include_directories("test")

list(APPEND LINK_ALL libmaxent ${ALPSCore_LIBRARIES} ${GSL_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
enable_testing(test)
set(test_src
  default_modelTest
  backcontTest
  gridTest
  parmsTest
  simulationTest
  paramFailureTest
)
foreach(test ${test_src})
  add_gtest(${test} test)
endforeach(test)

#add companion utilities
add_subdirectory(legendre_convert)
add_subdirectory(kk) 
if(PADE)
  add_subdirectory(pade/pade_arbitrary_degree)
endif(PADE)
