cmake_minimum_required(VERSION 2.8)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

project (kk)
include_directories ("/opt/local/include")
link_directories("/opt/local/lib")
link_directories("/opt/local/lib/gcc48/")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3  -DNDEBUG")


add_executable(kk kk.cpp )
target_link_libraries(kk alps boost gsl)
