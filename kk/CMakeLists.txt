# SPDX-License-Identifier: MIT
find_package(Boost 1.70 CONFIG REQUIRED COMPONENTS program_options)
find_package(OpenMP COMPONENTS CXX)

add_executable(kk kk.cpp)
target_link_libraries(kk PRIVATE Boost::program_options GSL::gsl maxent_warnings)
if(OpenMP_CXX_FOUND)
  target_link_libraries(kk PRIVATE OpenMP::OpenMP_CXX)
endif()
install(TARGETS kk RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
