cmake_minimum_required(VERSION 3.20)

project("AirSpacerModel")

# Check dependencies
find_package(ViennaLS REQUIRED)
find_package(LAPACK REQUIRED)

# Set default build type
set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${DEFAULT_BUILD_TYPE}"
      CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel"
                                               "RelWithDebInfo")
endif()

# Enable more warning levels
if (MSVC)
    add_compile_options(/W3)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

# Add include directories
list(APPEND PROJECT_INCLUDE_DIRS ${VIENNALS_INCLUDE_DIRS} include/)
# Add libraries
list(APPEND PROJECT_LIBRARIES ${VIENNALS_LIBRARIES} ${LAPACK_LIBRARIES})

# The main executable
set(EMULATION AirSpacerEmulation)
add_executable(${EMULATION} src/${EMULATION}.cpp)
target_include_directories(${EMULATION} PUBLIC ${PROJECT_INCLUDE_DIRS})
target_link_libraries(${EMULATION} PRIVATE ${PROJECT_LIBRARIES})

set(L2V lvst2vtp)
add_executable(${L2V} src/${L2V}.cpp)
target_include_directories(${L2V} PUBLIC ${PROJECT_INCLUDE_DIRS})
target_link_libraries(${L2V} PRIVATE ${PROJECT_LIBRARIES})

find_package(fmt QUIET)
if(fmt_FOUND)
  list(APPEND PROJECT_LIBRARIES fmt::fmt)

  set(FR FillingRatio)
  add_executable(${FR} src/${FR}.cpp)
  target_include_directories(${FR} PUBLIC ${PROJECT_INCLUDE_DIRS})
  target_link_libraries(${FR} PRIVATE ${PROJECT_LIBRARIES})

  set(SPLINE TestSpline)
  add_executable(${SPLINE} src/${SPLINE}.cpp)
  target_include_directories(${SPLINE} PUBLIC ${PROJECT_INCLUDE_DIRS})
  target_link_libraries(${SPLINE} PRIVATE ${PROJECT_LIBRARIES})
endif()

# Copy data and config files to the build directory
configure_file(data/data.csv ${CMAKE_CURRENT_BINARY_DIR}/data.csv COPYONLY)
configure_file(data/config.txt ${CMAKE_CURRENT_BINARY_DIR}/config.txt COPYONLY)

# Optional executables for generating the dataset as well as for testing
# the extraction and reconstruction functionality
set(VIENNAPS_VERBOSE OFF)
find_package(ViennaPS QUIET)
if(ViennaPS_FOUND)
  list(APPEND PROJECT_INCLUDE_DIRS ${VIENNAPS_INCLUDE_DIRS})
  list(APPEND PROJECT_LIBRARIES ${VIENNAPS_LIBRARIES})

  set(CREATE_DATA CreateData)
  add_executable(${CREATE_DATA} src/${CREATE_DATA}.cpp)
  target_include_directories(${CREATE_DATA} PUBLIC ${PROJECT_INCLUDE_DIRS})
  target_link_libraries(${CREATE_DATA} PRIVATE ${PROJECT_LIBRARIES})

  set(RECONSTRUCTION GeometryReconstruction)
  add_executable(${RECONSTRUCTION} src/${RECONSTRUCTION}.cpp)
  target_include_directories(${RECONSTRUCTION} PUBLIC ${PROJECT_INCLUDE_DIRS})
  target_link_libraries(${RECONSTRUCTION} PRIVATE ${PROJECT_LIBRARIES})

  set(SIMULATION AirSpacerSimulation)
  add_executable(${SIMULATION} src/${SIMULATION}.cpp)
  target_include_directories(${SIMULATION} PUBLIC ${PROJECT_INCLUDE_DIRS})
  target_link_libraries(${SIMULATION} PRIVATE ${PROJECT_LIBRARIES})
endif()
