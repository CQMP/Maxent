cmake_minimum_required(VERSION 2.8)
project (pade)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(ALPSCore REQUIRED COMPONENTS hdf5 accumulators mc params)
find_package (Eigen3 3.1 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
link_directories("/opt/local/lib")

add_executable(pade pade.cpp pade_grid.cpp pade_real.cpp pade_imag.cpp pade_interpolator.cpp pade_solver.cpp)
target_link_libraries(pade ${ALPSCore_LIBRARIES} gmpxx gmp)
install(TARGETS pade DESTINATION bin)
