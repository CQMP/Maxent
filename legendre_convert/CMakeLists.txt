cmake_minimum_required(VERSION 2.8)
project(legendre_convert)

add_executable(legendre_convert legendre_convert.cpp)

set(CMAKE_CXX_FLAGS  "-W -Wno-sign-compare -O2 -g")
# use ALPSCore for transitive dependencies on Boost:
target_link_libraries(legendre_convert ${ALPSCore_LIBRARIES})
install(TARGETS legendre_convert DESTINATION bin)
