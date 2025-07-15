cmake_minimum_required(VERSION 3.12)

# Build Application Unit Tests
PROJECT(Application_UnitTests)

#-------------------------------------------------------------------------------
#  Include paths and library directories
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

#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Build test executables for Application
#-------------------------------------------------------------------------------


# Application Test Executables

# Application Tests
# ADD_EXECUTABLE(Test_ImmersedBoundary.exe Test_ImmersedBoundary.F90)
# TARGET_LINK_LIBRARIES(Test_ImmersedBoundary.exe Application Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# Conventional Analysis Tests
ADD_EXECUTABLE(UnitTest_sfc_randoms.exe UnitTest_sfc_randoms.F90)
TARGET_LINK_LIBRARIES(UnitTest_sfc_randoms.exe Application Geometry MPDD YAMLRead Utility State Field ObsSet ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(UnitTest_Convention.exe UnitTest_Convention.F90)
# TARGET_LINK_LIBRARIES(UnitTest_Convention.exe Application Geometry MPDD YAMLRead Utility State Field ObsSet ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(UnitTest_Convention_MultiObsVcycle.exe UnitTest_Convention_MultiObsVcycle.F90)
# TARGET_LINK_LIBRARIES(UnitTest_Convention_MultiObsVcycle.exe Application Geometry MPDD YAMLRead Utility State Field ObsSet ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(UnitTest_Conventions.exe UnitTest_Conventions.F90)
# TARGET_LINK_LIBRARIES(UnitTest_Conventions.exe Application Geometry MPDD YAMLRead Utility State Field ObsSet ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(UnitTest_Convention_MultiObs.exe UnitTest_Convention_MultiObs.F90)
# TARGET_LINK_LIBRARIES(UnitTest_Convention_MultiObs.exe Application Geometry MPDD YAMLRead Utility State Field ObsSet ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(UnitTest_Convention_Adv.exe UnitTest_Convention_Adv.F90)
# TARGET_LINK_LIBRARIES(UnitTest_Convention_Adv.exe Application Geometry MPDD YAMLRead Utility State Field ObsSet ${MPI_Fortran_LIBRARIES})

#-------------------------------------------------------------------------------
# CTest configurations for Application
#-------------------------------------------------------------------------------

# Add test from original ConventionalAnalysis
ADD_TEST(NAME "Test_sfc_random" COMMAND mpirun -np 4 --oversubscribe UnitTest_sfc_randoms.exe
	WORKING_DIRECTORY $ENV{BIN_DIR})

# Set test criteria
SET_TESTS_PROPERTIES(Test_sfc_random PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
SET_TESTS_PROPERTIES(Test_sfc_random PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")

