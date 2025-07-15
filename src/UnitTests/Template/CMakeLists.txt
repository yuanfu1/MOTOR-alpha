cmake_minimum_required(VERSION 3.12)

# Build Template Unit Tests
PROJECT(Template_UnitTests)

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
# Build test executables for Template
#-------------------------------------------------------------------------------


# Template Test Executables

# Template Tests
# ADD_EXECUTABLE(Test_YM.exe Test_YM.F90)
# TARGET_LINK_LIBRARIES(Test_YM.exe Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_Temp.exe Test_Temp.F90)
# TARGET_LINK_LIBRARIES(Test_Temp.exe Utility ${MPI_Fortran_LIBRARIES})

