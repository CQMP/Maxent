cmake_minimum_required(VERSION 2.8)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

project (MAXENT)

add_executable(maxent 
maxent.cpp
maxent_helper.cpp
maxent_parms.cpp
maxent_simulation.cpp)
target_link_libraries(maxent
${ALPS_LIBRARIES} )

