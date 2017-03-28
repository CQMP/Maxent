cmake_minimum_required(VERSION 2.8.12)
project (MAXENT)

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
find_package(ALPSCore REQUIRED)

#note: LAPACK also finds BLAS
find_package(LAPACK)
if(LAPACK_FOUND AND USE_LAPACK)
    set(HAVE_BLAS 1)
    set(HAVE_LAPACK 1)
    message(STATUS "Using LAPACK and BLAS routines for SVD")
else()
    set(HAVE_BLAS 0)
    set(HAVE_LAPACK 0)
endif()
configure_file(
	"${PROJECT_SOURCE_DIR}/cmake/maxent_config.hpp.in"
	"${PROJECT_BINARY_DIR}/config/maxent_config.hpp"
    )

include_directories("${PROJECT_BINARY_DIR}/config")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
message(STATUS "Finding packages, please stand by...")
find_package(GSL REQUIRED)
find_package (Eigen3 3.1 REQUIRED)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-Wno-return-type-c-linkage SUPPORTS_FLAG)
if(SUPPORTS_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-return-type-c-linkage")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fno-omit-frame-pointer -g")

#let gcc take advantage of Eigen3 vectorization
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
endif()

include_directories(${GSL_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

set(LIB_FILES 
  ./src/maxent_helper.cpp
  ./src/maxent_parms.cpp
  ./src/maxent_svd.cpp
  ./src/maxent_grid.cpp
  ./src/maxent_kernel.cpp
  ./src/maxent_simulation.cpp
  ./src/maxent_backcont.cpp
  ./src/default_model.cpp
  ./src/eigen_hdf5.cpp
  ./src/SpM_ADMM.cpp
  )

ADD_LIBRARY(libmaxent ${LIB_FILES})
target_link_libraries(libmaxent
	${ALPSCore_LIBRARIES} ${GSL_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
#remove default "lib" prefix
set_target_properties(libmaxent PROPERTIES PREFIX "")
#executable
add_executable(maxent 
  ./src/maxent.cpp
)
target_link_libraries(maxent libmaxent ${ALPSCore_LIBRARIES} ${GSL_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(SpM
  ./src/SpM_method.cpp ./src/SpM_simulation.cpp 
)
target_link_libraries(SpM libmaxent ${ALPSCore_LIBRARIES} ${GSL_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

#testing setup
option(Testing "Enable testing" ON)


list(APPEND LINK_ALL libmaxent ${ALPSCore_LIBRARIES} ${GSL_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
enable_testing(test)
add_subdirectory(test)

#add companion utilities
add_subdirectory(legendre_convert)
add_subdirectory(kk) 
if(PADE)
  add_subdirectory(pade/pade_arbitrary_degree)
endif(PADE)

#install
install(TARGETS maxent DESTINATION bin)
