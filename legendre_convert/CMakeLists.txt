cmake_minimum_required(VERSION 2.8)
project(legendre_convert)

add_executable(legendre_convert legendre_convert.cpp)

set(CMAKE_CXX_FLAGS  "-W -Wno-sign-compare -O2 -g")
set(boost_req program_options)
find_package (Boost 1.54.0 COMPONENTS ${boost_req} REQUIRED)
message(STATUS "Found Boost includes: ${Boost_INCLUDE_DIRS}" )
message(STATUS "Found Boost libs: ${Boost_LIBRARIES}" )
include_directories(${Boost_INCLUDE_DIRS})

target_link_libraries(legendre_convert ${ALPSCore_LIBRARIES} ${Boost_LIBRARIES})
install(TARGETS legendre_convert DESTINATION bin)
