cmake_minimum_required(VERSION 2.8)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

set (GTEST_ROOT gtest )
set (GTEST_INCLUDE_DIR ${GTEST_ROOT}/include )
set (GTEST_SRC_DIR ${GTEST_ROOT}/src)

include_directories(${CMAKE_SOURCE_DIR}/gtest/include)
include_directories(${GTEST_INCLUDE_DIR})
include_directories(${GTEST_ROOT})

project (MAXENT)

add_executable(maxent 
maxent.cpp
maxent_helper.cpp
maxent_parms.cpp
maxent_simulation.cpp
default_modelTest.cpp
default_model.cpp
${GTEST_SRC_DIR}/gtest-all.cc
)
target_link_libraries(maxent
${ALPS_LIBRARIES} )

