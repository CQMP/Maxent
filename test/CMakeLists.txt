# SPDX-License-Identifier: MIT
if(MAXENT_USE_SYSTEM_GTEST)
  find_package(GTest 1.18 CONFIG REQUIRED)
else()
  include(FetchContent)
  FetchContent_Declare(googletest
    URL https://github.com/google/googletest/releases/download/v1.18.0/googletest-1.18.0.tar.gz
    URL_HASH SHA256=6e3191c1455468b3fc35a417fb565c1c5071aee1b7e7f85e30cf48a98d37d8b5
    DOWNLOAD_EXTRACT_TIMESTAMP ON)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()
include(GoogleTest)

set(unit_tests
  backcontTest
  default_modelTest
  gridTest
  paramFailureTest
  paramsTest
  simulationTest)

foreach(test IN LISTS unit_tests)
  add_executable(${test} ${test}.cpp)
  # GoogleTest first: the ALPSCore install ships its own (older) gtest headers in
  # its include directory, which must not shadow the ones matching this library.
  target_link_libraries(${test} PRIVATE GTest::gtest_main maxent::core maxent_warnings)
  gtest_discover_tests(${test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    PROPERTIES LABELS unit)
endforeach()

add_subdirectory(regression)
