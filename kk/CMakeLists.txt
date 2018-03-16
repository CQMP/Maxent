cmake_minimum_required(VERSION 2.8)
#include(${ALPS_USE_FILE})

project (kk)
find_package(OpenMP)
find_package(GSL)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
find_package(ALPSCore REQUIRED)
include_directories ("/opt/local/include")
link_directories("/opt/local/lib")
#link_directories("/opt/local/lib/gcc48/")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3  -DNDEBUG")


add_executable(kk kk.cpp )
target_link_libraries(kk ${ALPSCore_LIBRARIES} ${GSL_LIBRARIES} ${Boost_LIBRARIES})
install(TARGETS kk DESTINATION bin)
