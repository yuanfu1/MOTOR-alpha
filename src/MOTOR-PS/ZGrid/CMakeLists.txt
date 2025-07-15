cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(ZGrid)

#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
include_directories($ENV{COMMON_MODS})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
link_directories($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include(srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
add_library(zgrid STATIC ${zgrid_SRC})
target_link_libraries(zgrid Field State MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys State2NC)

# Create the GZM library
add_library(gzm_lib
    gzm_Xie.F90
)

# Set properties
set_target_properties(gzm_lib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
    POSITION_INDEPENDENT_CODE ON
)

# Link dependencies
target_link_libraries(gzm_lib
    PUBLIC
    Geometry
    State
)

# Include directories
target_include_directories(gzm_lib
    PUBLIC
    ${CMAKE_BINARY_DIR}/modules
)

# Installation
install(TARGETS gzm_lib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES
    ${CMAKE_BINARY_DIR}/modules/gzm_m.mod
    DESTINATION include
)

# 2. Build test executables:
# set(EXE1 "Test_gzm.exe")

# 3. Add codes to the executables:
# add_executable(${EXE1} Test_gzm.F90)
# target_link_libraries(${EXE1} zgrid PoissonSolver MPDD Geometry Utility FLog ${MPI_Fortran_LIBRARIES} PreSys State2NC)

# Add required dependencies - only keep MPI which is external
find_package(MPI REQUIRED)

# Add Laplacian adjoint test
# add_executable(test_laplacian_adjoint test_laplacian_adjoint.F90)
# target_link_libraries(test_laplacian_adjoint
#     PRIVATE
#     gzm_lib
#     State
#     Geometry
#     Field
#     MPDD
#     Utility
#     YAMLRead
#     PreSys
#     State2NC
#     ${MPI_Fortran_LIBRARIES}
#     ${NETCDF_LIBRARIES}
# )

# # Make sure directories are included first
# target_include_directories(test_laplacian_adjoint
#     PRIVATE
#     ${CMAKE_BINARY_DIR}/modules
#     ${MPI_Fortran_INCLUDE_PATH}
#     $ENV{COMMON_MODS}
# )

# Test configuration migrated to UnitTests/MOTOR-PS/CMakeLists.txt
# add_test(NAME LaplacianAdjoint 
#          COMMAND test_laplacian_adjoint
#          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# 4. Add tests:
# ADD_TEST(NAME "Test_gzm" 
#         COMMAND mpirun -np 4 --oversubscribe Test_gzm.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})

# 5. Set test criteria:
# SET_TESTS_PROPERTIES(Test_gzm PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
# SET_TESTS_PROPERTIES(Test_gzm PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
