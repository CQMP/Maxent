cmake_minimum_required(VERSION 2.8)

project (kk)
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
find_package(ALPSCore REQUIRED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3  -DNDEBUG")


add_executable(kk kk.cpp )
target_link_libraries(kk ${ALPSCore_LIBRARIES} ${GSL_LIBRARIES})
install(TARGETS kk DESTINATION bin)
