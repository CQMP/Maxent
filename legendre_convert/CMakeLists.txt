cmake_minimum_required(VERSION 2.8)
project(legendre_convert)

add_executable(legendre_convert legendre_convert.cpp)

set(CMAKE_CXX_FLAGS  "-W -Wno-sign-compare -O2 -g")
find_package(Boost 1.56 REQUIRED COMPONENTS program_options)
target_link_libraries(legendre_convert ${Boost_LIBRARIES})
target_include_directories(legendre_convert PUBLIC ${Boost_INCLUDE_DIRS})
install(TARGETS legendre_convert DESTINATION bin)
