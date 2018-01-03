cmake_minimum_required(VERSION 2.8)

project (kk)
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# find_package(ALPSCore REQUIRED) # does it use ALPSCore?
set(boost_req program_options)
find_package (Boost COMPONENTS ${boost_req} REQUIRED)
message(STATUS "Found Boost includes: ${Boost_INCLUDE_DIRS}" )
message(STATUS "Found Boost libs: ${Boost_LIBRARIES}" )


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3  -DNDEBUG")


add_executable(kk kk.cpp )
target_include_directories(kk PUBLIC ${Boost_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})
target_link_libraries(kk ${GSL_LIBRARIES} ${Boost_LIBRARIES})
install(TARGETS kk DESTINATION bin)
