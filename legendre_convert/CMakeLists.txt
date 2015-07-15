cmake_minimum_required(VERSION 2.8)

project(legendre_convert)

include_directories("/opt/local/include/")
link_directories("/opt/local/lib/")
add_executable(legendre_convert legendre_convert.cpp)

set(CMAKE_CXX_FLAGS  "-W -Wno-sign-compare -O2 -g")
#target_link_libraries("/opt/local/lib/")
target_link_libraries(legendre_convert boost_program_options-mt boost_random-mt)
