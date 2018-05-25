cmake_minimum_required(VERSION 2.8)
project (pade)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(ALPSCore REQUIRED COMPONENTS hdf5 accumulators mc params)
find_package (Eigen3 3.1 REQUIRED)
set(boost_req program_options)
find_package (Boost 1.54.0 COMPONENTS ${boost_req} REQUIRED)
message(STATUS "Found Boost includes: ${Boost_INCLUDE_DIRS}" )
message(STATUS "Found Boost libs: ${Boost_LIBRARIES}" )

include_directories(${EIGEN3_INCLUDE_DIR})
link_directories("/opt/local/lib")

add_executable(pade pade.cpp pade_grid.cpp pade_real.cpp pade_imag.cpp pade_interpolator.cpp pade_solver.cpp)
target_link_libraries(pade ${ALPSCore_LIBRARIES} gmpxx gmp ${Boost_LIBRARIES})
install(TARGETS pade DESTINATION bin)
