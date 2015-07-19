cmake_minimum_required(VERSION 2.8)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

SET(CMAKE_CXX_FLAGS "-O3 -DNDEBUG")
#SET(CMAKE_C_FLAGS "-O3 -DNDEBUG")
#set(CMAKE_CXX_COMPILER g++)
#set(CMAKE_C_COMPILER gcc)
include_directories("/opt/local/include")
link_directories("/opt/local/lib")
project (pade)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost -fopenmp -mkl -stdlib=libstdc++")

add_executable(pade pade.cpp pade_grid.cpp pade_real.cpp pade_imag.cpp pade_interpolator.cpp pade_solver.cpp)
target_link_libraries(pade ${ALPS_LIBRARIES} mpfr gmpxx gmp)
