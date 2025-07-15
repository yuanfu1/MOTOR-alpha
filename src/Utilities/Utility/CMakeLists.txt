cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
PROJECT(Utility)

INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

# Add netcdf libraries.
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})
LINK_DIRECTORIES(${NETCDF_LIBRARY})

# Add netcdf libraries.
# INCLUDE_DIRECTORIES(${SLINT_Fortran_INCLUDE_DIR})
# LINK_DIRECTORIES(${SLINT_LIBRARY})
#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Source codes to this build:
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(Utility STATIC ${Utility_SRC})
TARGET_LINK_LIBRARIES(Utility PUBLIC lapack blas ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_Shared} ${NETCDF_LIBRARY_C_Shared})

# Add the AMG solver
#add_subdirectory(AMG)

# 2. Build test executables:
#SET (EXE2 "Test_GeoTools.exe")
# SET (EXE3 "interp1d_test.exe")
#SET (EXE4 "Test_AdvTime.exe")
#SET (EXE5 "Test_RegularGrid.exe")
#SET (EXE6 "Test_Netcdf.exe")
#SET (EXE7 "Test_Derivative.exe")
#SET (EXE8 "Test_GMRES.exe")
#SET (EXE9 "Test_ECMIO.exe")

# 3. Add code to the executables:
#ADD_EXECUTABLE(${EXE2} Test_GeoTools.F90)
# ADD_EXECUTABLE(${EXE3} interp1d_test.F90)
#ADD_EXECUTABLE(${EXE4} Test_AdvTime.F90)
#ADD_EXECUTABLE(${EXE5} Test_RegularGrid.F90)
#ADD_EXECUTABLE(${EXE6} Test_Netcdf.F90)
#ADD_EXECUTABLE(${EXE7} Test_Derivative.F90)
#ADD_EXECUTABLE(${EXE8} test_gmres.F90)
#ADD_EXECUTABLE(${EXE9} Test_ECMIO.F90)
#  LINK_DIRECTORIES(dir1 dir2 ...)
#-------------------------------------------------------------------------------
#  To link another lib like -l in cmake:
#-------------------------------------------------------------------------------
#  TARGET_LINK_LIBRARIES(target <list of libnames>)
# TARGET_LINK_LIBRARIES(${EXE2} Utility)
# TARGET_LINK_LIBRARIES(${EXE3} Utility)
# TARGET_LINK_LIBRARIES(${EXE4} Utility)
# TARGET_LINK_LIBRARIES(${EXE5} Utility)
# TARGET_LINK_LIBRARIES(${EXE6} Utility)
# TARGET_LINK_LIBRARIES(${EXE7} Utility) 
# TARGET_LINK_LIBRARIES(${EXE7} Utility)
# TARGET_LINK_LIBRARIES(${EXE8} Utility)
# TARGET_LINK_LIBRARIES(${EXE9} Utility)

# # 4. Add tests: (MIGRATED TO UnitTests/Utilities/CMakeLists.txt)
# ADD_TEST(NAME Utility_Interp1D
#          COMMAND interp1d_test.exe 
#          WORKING_DIRECTORY $ENV{BIN_DIR})
# ADD_TEST(NAME Utility_RegularGrid
#          COMMAND Test_RegularGrid.exe 
#          WORKING_DIRECTORY $ENV{BIN_DIR})
#ADD_TEST(NAME DataStructures
#         COMMAND TestDataStructures.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})

# # 5. Set test criteria: (MIGRATED TO UnitTests/Utilities/CMakeLists.txt)
# SET_TESTS_PROPERTIES(Utility_Interp1D PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed")
# SET_TESTS_PROPERTIES(Utility_Interp1D PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
# SET_TESTS_PROPERTIES(Utility_RegularGrid PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed")
# SET_TESTS_PROPERTIES(Utility_RegularGrid PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
# SET_TESTS_PROPERTIES(DataStructures PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed")
# SET_TESTS_PROPERTIES(DataStructures PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
