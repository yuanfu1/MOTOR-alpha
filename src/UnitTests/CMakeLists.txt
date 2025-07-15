cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(UnitTests)

# Enable CTest
ENABLE_TESTING()

#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

# Add Netcdf includes
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})

# Add rttov libraries.
INCLUDE_DIRECTORIES(${RTTOV_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${RTTOV_MOD_DIR})

# Add HDF5 libraries.
IF (CMAKE_Fortran_COMPILER MATCHES "ifort")
	INCLUDE_DIRECTORIES(${HDF5_Fortran_INCLUDE_DIR})
ENDIF ()

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Source codes to this build:
# include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:

# 1. Build test executables:
SET (EXE1  "TestShallowWater.exe")
SET (EXE2  "TestMixedSolver.exe")
SET (EXE3  "TestPoisson.exe")

# 2. Add codes to the executables:
# # EXE1
ADD_EXECUTABLE(${EXE1} TestShallowWater.F90)
ADD_EXECUTABLE(${EXE2} TestMixedSolver.F90)
ADD_EXECUTABLE(${EXE3} Test_Poisson.F90)
TARGET_LINK_LIBRARIES(${EXE1} rhsBase TimeIntegration C2M BMatrix UV2W BkgdMG IOGrapes IOWRF lapack blas)
TARGET_LINK_LIBRARIES(${EXE2} PoissonSolver)
TARGET_LINK_LIBRARIES(${EXE3} PoissonSolver)

# 4. Add subdirectories for different test modules:
ADD_SUBDIRECTORY(MOTOR-DA)
ADD_SUBDIRECTORY(MOTOR-PS)
ADD_SUBDIRECTORY(MOTOR-QC)
ADD_SUBDIRECTORY(MOTOR-DP)
ADD_SUBDIRECTORY(MOTOR-Repository)
ADD_SUBDIRECTORY(Utilities)
ADD_SUBDIRECTORY(Application)
ADD_SUBDIRECTORY(Template)

# # 5. Add tests:
# ADD_TEST(NAME TestShallowWater
# COMMAND TestShallowWater.exe
# WORKING_DIRECTORY $ENV{BIN_DIR})

# ADD_TEST(NAME TestMixedSolver
# COMMAND TestMixedSolver.exe
# WORKING_DIRECTORY $ENV{BIN_DIR})

# ADD_TEST(NAME TestPoisson
# COMMAND TestPoisson.exe
# WORKING_DIRECTORY $ENV{BIN_DIR})

# # 6. Set test criteria:
# SET_TESTS_PROPERTIES(TestShallowWater PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed")
# SET_TESTS_PROPERTIES(TestShallowWater PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")

# SET_TESTS_PROPERTIES(TestMixedSolver PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed")
# SET_TESTS_PROPERTIES(TestMixedSolver PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")

# SET_TESTS_PROPERTIES(TestPoisson PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed")
# SET_TESTS_PROPERTIES(TestPoisson PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
