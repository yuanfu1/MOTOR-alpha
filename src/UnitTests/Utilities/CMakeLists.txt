cmake_minimum_required(VERSION 3.12)

# Build Utilities Unit Tests
PROJECT(Utilities_UnitTests)

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
# Build test executables for Utilities
#-------------------------------------------------------------------------------


# Utilities Test Executables

# Utility Function Tests
ADD_EXECUTABLE(Test_FLog.exe Test_FLog.F90)
TARGET_LINK_LIBRARIES(Test_FLog.exe FLog Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_updatebdy.exe Test_updatebdy.F90)
TARGET_LINK_LIBRARIES(Test_updatebdy.exe UpdateBdy Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_YAMLRead.exe Test_YAMLRead.F90)
TARGET_LINK_LIBRARIES(Test_YAMLRead.exe YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_GrapesIO.exe Test_GrapesIO.F90)
TARGET_LINK_LIBRARIES(Test_GrapesIO.exe GrapesIO Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_qcqrIO.exe Test_qcqrIO.F90)
TARGET_LINK_LIBRARIES(Test_qcqrIO.exe GrapesIO Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_GrapesModelvarIO_500m.exe Test_GrapesModelvarIO_500m.F90)
TARGET_LINK_LIBRARIES(Test_GrapesModelvarIO_500m.exe GrapesIO Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_GrapesPostvarIO.exe Test_GrapesPostvarIO.F90)
# TARGET_LINK_LIBRARIES(Test_GrapesPostvarIO.exe GrapesIO Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_GrapesModelvarIO.exe Test_GrapesModelvarIO.F90)
TARGET_LINK_LIBRARIES(Test_GrapesModelvarIO.exe GrapesIO Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_NMLRead.exe Test_NMLRead.F90)
TARGET_LINK_LIBRARIES(Test_NMLRead.exe NMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_IntepHP.exe Test_IntepHP.F90)
# TARGET_LINK_LIBRARIES(Test_IntepHP.exe InterpHP Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_VerticalInterp.exe Test_VerticalInterp.F90)
TARGET_LINK_LIBRARIES(Test_VerticalInterp.exe Utility MPDD ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_Derivative.exe Test_Derivative.F90)
TARGET_LINK_LIBRARIES(Test_Derivative.exe Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_Specials.exe Test_Specials.F90)
TARGET_LINK_LIBRARIES(Test_Specials.exe Geometry ZGrid Dyn specialFunctions Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_RegularGrid.exe Test_RegularGrid.F90)
TARGET_LINK_LIBRARIES(Test_RegularGrid.exe Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_Netcdf.exe Test_Netcdf.F90)
TARGET_LINK_LIBRARIES(Test_Netcdf.exe Utility ${MPI_Fortran_LIBRARIES} ${NETCDF_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_AdvTime.exe Test_AdvTime.F90)
TARGET_LINK_LIBRARIES(Test_AdvTime.exe Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_ECMIO.exe Test_ECMIO.F90)
TARGET_LINK_LIBRARIES(Test_ECMIO.exe Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_GeoTools.exe Test_GeoTools.F90)
TARGET_LINK_LIBRARIES(Test_GeoTools.exe Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(test_amg2d.exe test_amg2d.F90)
# TARGET_LINK_LIBRARIES(test_amg2d.exe AMG Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(test_gmres.exe test_gmres.F90)
TARGET_LINK_LIBRARIES(test_gmres.exe Utility ${MPI_Fortran_LIBRARIES})

# Unit Test Dynamic Functions

ADD_EXECUTABLE(interp1d_test.exe interp1d_test.F90)
TARGET_LINK_LIBRARIES(interp1d_test.exe Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(TestDataStructures.exe TestDataStructures.F90)
# TARGET_LINK_LIBRARIES(TestDataStructures.exe Utility ${MPI_Fortran_LIBRARIES})

#-------------------------------------------------------------------------------
# CTest configurations for Utilities
#-------------------------------------------------------------------------------

# Add tests from original NMLRead
ADD_TEST(NAME Test_NMLRead 
        COMMAND Test_NMLRead.exe
        WORKING_DIRECTORY ${BIN_DIR})

# Add tests from original Utility
ADD_TEST(NAME Utility_Interp1D
         COMMAND interp1d_test.exe 
         WORKING_DIRECTORY $ENV{BIN_DIR})
ADD_TEST(NAME Utility_RegularGrid
         COMMAND Test_RegularGrid.exe 
         WORKING_DIRECTORY $ENV{BIN_DIR})
# Add a test: 
ADD_TEST(NAME verticalInterp
         COMMAND Test_VerticalInterp.exe
         WORKING_DIRECTORY $ENV{BIN_DIR})

# Set test criteria
SET_TESTS_PROPERTIES(Test_NMLRead PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed")
SET_TESTS_PROPERTIES(Test_NMLRead PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")

SET_TESTS_PROPERTIES(Utility_Interp1D PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed")
SET_TESTS_PROPERTIES(Utility_Interp1D PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
SET_TESTS_PROPERTIES(Utility_RegularGrid PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed")
SET_TESTS_PROPERTIES(Utility_RegularGrid PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")

SET_TESTS_PROPERTIES(verticalInterp PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed")
SET_TESTS_PROPERTIES(verticalInterp PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
