cmake_minimum_required(VERSION 3.14)

project("AirSpacerModel")

# Check dependencies
find_package(ViennaLS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(fmt QUIET)
set(VIENNAPS_VERBOSE OFF)
find_package(ViennaPS QUIET)

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
list(APPEND PROJECT_LIBRARIES stdc++fs ${VIENNALS_LIBRARIES} ${LAPACK_LIBRARIES})

set(EMULATION AirSpacerEmulation)
add_executable(${EMULATION} src/${EMULATION}.cpp)
target_include_directories(${EMULATION} PUBLIC ${PROJECT_INCLUDE_DIRS})
target_link_libraries(${EMULATION} PRIVATE ${PROJECT_LIBRARIES})
configure_file(data/airspacer_dataset.csv ${CMAKE_CURRENT_BINARY_DIR}/airspacer_dataset.csv COPYONLY)
configure_file(data/airspacer.txt ${CMAKE_CURRENT_BINARY_DIR}/airspacer.txt COPYONLY)

set(TRENCHFILL TrenchFill)
add_executable(${TRENCHFILL} src/${TRENCHFILL}.cpp)
target_include_directories(${TRENCHFILL} PUBLIC ${PROJECT_INCLUDE_DIRS})
target_link_libraries(${TRENCHFILL} PRIVATE ${PROJECT_LIBRARIES})
configure_file(data/trenchfill_dataset.csv ${CMAKE_CURRENT_BINARY_DIR}/trenchfill_dataset.csv COPYONLY)
configure_file(data/trenchfill.txt ${CMAKE_CURRENT_BINARY_DIR}/trenchfill.txt COPYONLY)
if(ViennaPS_FOUND)
  target_compile_definitions(${TRENCHFILL} PRIVATE WITH_VIENNAPS)
  target_include_directories(${TRENCHFILL} PUBLIC ${VIENNAPS_INCLUDE_DIRS})
  target_link_libraries(${TRENCHFILL} PRIVATE ${VIENNAPS_LIBRARIES})
endif()

# Optional executables for generating the dataset as well as for testing
# the extraction and reconstruction functionality
if(ViennaPS_FOUND)
  # The reference simulation to check the results of the AirSpacerEmulation
  set(SIMULATION AirSpacerSimulation)
  add_executable(${SIMULATION} src/${SIMULATION}.cpp)
  target_include_directories(${SIMULATION} PUBLIC ${PROJECT_INCLUDE_DIRS} ${VIENNAPS_INCLUDE_DIRS})
  target_link_libraries(${SIMULATION} PRIVATE ${PROJECT_LIBRARIES} ${VIENNAPS_LIBRARIES})

  # The dataset generation application
  set(CREATE_DATA CreateData)
  add_executable(${CREATE_DATA} src/${CREATE_DATA}.cpp)
  target_include_directories(${CREATE_DATA} PUBLIC ${PROJECT_INCLUDE_DIRS} ${VIENNAPS_INCLUDE_DIRS})
  target_link_libraries(${CREATE_DATA} PRIVATE ${PROJECT_LIBRARIES} ${VIENNAPS_LIBRARIES})
  configure_file(data/create_airspacer_dataset.txt ${CMAKE_CURRENT_BINARY_DIR}/create_airspacer_dataset.txt COPYONLY)
  configure_file(data/create_trenchfill_dataset.txt ${CMAKE_CURRENT_BINARY_DIR}/create_trenchfill_dataset.txt COPYONLY)

  # Program that test the extraction and reconstruction mechanism
  set(RECONSTRUCTION GeometryReconstruction)
  add_executable(${RECONSTRUCTION} src/${RECONSTRUCTION}.cpp)
  target_include_directories(${RECONSTRUCTION} PUBLIC ${PROJECT_INCLUDE_DIRS} ${VIENNAPS_INCLUDE_DIRS})
  target_link_libraries(${RECONSTRUCTION} PRIVATE ${PROJECT_LIBRARIES} ${VIENNAPS_LIBRARIES})
endif()

# Utility for converting lvst to vtp files
set(L2V lvst2vtp)
add_executable(${L2V} src/${L2V}.cpp)
target_include_directories(${L2V} PUBLIC ${PROJECT_INCLUDE_DIRS})
target_link_libraries(${L2V} PRIVATE ${PROJECT_LIBRARIES})

if(fmt_FOUND)
  list(APPEND PROJECT_LIBRARIES fmt::fmt)

  # Application for demoing the spline interpolation functionality
  set(SPLINE TestSpline)
  add_executable(${SPLINE} src/${SPLINE}.cpp)
  target_include_directories(${SPLINE} PUBLIC ${PROJECT_INCLUDE_DIRS})
  target_link_libraries(${SPLINE} PRIVATE ${PROJECT_LIBRARIES})
endif()
