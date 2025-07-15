# AMG solver library and test program
cmake_minimum_required(VERSION 3.10)

# Add this directory to the module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../../cmake")

# Define locations of required libraries
set(KINDS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../../src/Utilities/Utility/kinds")
set(GEOMETRY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../../src/MOTOR-Repository/StructuredGrid/geometry")
set(GZM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../../src/MOTOR-Repository/StructuredGrid/gzm")
set(STATE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../../src/MOTOR-Repository/StructuredGrid/state")

# Add dependencies as subdirectories if they're part of the project
add_subdirectory(${KINDS_DIR} kinds_build)
add_subdirectory(${GEOMETRY_DIR} geometry_build)
add_subdirectory(${GZM_DIR} gzm_build)
add_subdirectory(${STATE_DIR} state_build)

# Create the AMG library
add_library(amg2d_lib
    amg2d_m.F90
)

# Set the library properties
set_target_properties(amg2d_lib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
    POSITION_INDEPENDENT_CODE ON
)

# Create the test executable
add_executable(test_amg2d
    test_amg2d.F90
)

# Link the test program with the library
target_link_libraries(test_amg2d
    PRIVATE
    amg2d_lib
    kinds_lib
    geometry_lib
    gzm_lib
    state_lib
)

# Include directories
target_include_directories(amg2d_lib
    PUBLIC
    ${CMAKE_BINARY_DIR}/modules
    ${KINDS_DIR}
    ${GEOMETRY_DIR}
    ${GZM_DIR}
    ${STATE_DIR}
)

target_include_directories(test_amg2d
    PUBLIC
    ${CMAKE_BINARY_DIR}/modules
)

# Installation rules
install(TARGETS amg2d_lib test_amg2d
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
)

install(FILES
    ${CMAKE_BINARY_DIR}/modules/amg2d_m.mod
    DESTINATION include
)
