cmake_minimum_required(VERSION 2.8.12)
project (MAXENT)

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
find_package(ALPSCore REQUIRED)
#include(${ALPS_USE_FILE})
#BLAS and LAPACK find doesn't work well in c++ on apple
#hard coding -framework Accelerate until better solution found
#find_package(BLAS QUIET)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(UseGtest)

#set (GTEST_ROOT gtest )
set (GTEST_INCLUDE_DIR ${GTEST_ROOT}/include )
set (GTEST_SRC_DIR ${GTEST_ROOT}/src)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2") 

include_directories(${CMAKE_SOURCE_DIR}/gtest/include)
include_directories(${GTEST_INCLUDE_DIR})
include_directories(${GTEST_ROOT})
include_directories("/usr/local/Cellar/gsl/1.16/include")
link_directories("/usr/local/Cellar/gsl/1.16/lib")

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
default_modelTest.cpp
default_model.cpp
${GTEST_SRC_DIR}/gtest-all.cc
)
target_link_libraries(maxent
	${ALPSCore_LIBRARIES} gsl ${LAPACK_LIBRARIES} "-framework Accelerate")

