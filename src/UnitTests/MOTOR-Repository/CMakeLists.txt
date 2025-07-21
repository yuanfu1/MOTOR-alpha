cmake_minimum_required(VERSION 3.12)

# Build MOTOR-Repository Unit Tests
PROJECT(MOTOR-Repository_UnitTests)

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
# Build test executables for MOTOR-Repository
#-------------------------------------------------------------------------------


# MOTOR-Repository Test Executables

# Geometry Tests
ADD_EXECUTABLE(Test_Halo.exe Test_Halo.F90)
TARGET_LINK_LIBRARIES(Test_Halo.exe IOGrapes State2NC Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_Geometry.exe Test_Geometry.F90)
TARGET_LINK_LIBRARIES(Test_Geometry.exe Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# Multi-Grid Tests
ADD_EXECUTABLE(test_mgGen.exe test_mgGen.F90)
TARGET_LINK_LIBRARIES(test_mgGen.exe MGGen Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(test_icos.exe test_icos.F90)
TARGET_LINK_LIBRARIES(test_icos.exe MGGen Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(TestEdgeSetting.exe TestEdgeSetting.F90)
TARGET_LINK_LIBRARIES(TestEdgeSetting.exe MGGen Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

#-------------------------------------------------------------------------------
# CTest configurations for MOTOR-Repository
#-------------------------------------------------------------------------------

# Add test from original Geometry
ADD_TEST(NAME Test_Geometry
        COMMAND mpirun -np 4 --oversubscribe Test_Geometry.exe
        WORKING_DIRECTORY $ENV{BIN_DIR})

# Set test criteria
SET_TESTS_PROPERTIES(Test_Geometry PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")

