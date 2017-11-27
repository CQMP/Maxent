cmake_minimum_required(VERSION 2.8)

project (kk)
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
find_package(Boost 1.56 REQUIRED COMPONENTS program_options)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3  -DNDEBUG")


add_executable(kk kk.cpp )
target_link_libraries(kk ${Boost_LIBRARIES} ${GSL_LIBRARIES})
target_include_directories(kk PUBLIC ${Boost_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})

install(TARGETS kk DESTINATION bin)
