cmake_minimum_required(VERSION 3.12)

# Build MOTOR-DP Unit Tests
PROJECT(MOTOR-DP_UnitTests)

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
# Build test executables for MOTOR-DP
#-------------------------------------------------------------------------------


# MOTOR-DP Test Executables

# Satellite Processing Tests
# ADD_EXECUTABLE(Test_wrfinput_App.exe Test_wrfinput_App.F90)
# TARGET_LINK_LIBRARIES(Test_wrfinput_App.exe DP_Satellite Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES} ${NETCDF_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_Raw_fy3_4_mwts2.exe Test_Raw_fy3_4_mwts2.F90)
# TARGET_LINK_LIBRARIES(Test_Raw_fy3_4_mwts2.exe DP_Satellite Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_Raw_fy4_giirs.exe Test_Raw_fy4_giirs.F90)
# TARGET_LINK_LIBRARIES(Test_Raw_fy4_giirs.exe DP_Satellite Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_Raw_fy3_4_mwhs2.exe Test_Raw_fy3_4_mwhs2.F90)
# TARGET_LINK_LIBRARIES(Test_Raw_fy3_4_mwhs2.exe DP_Satellite Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_Raw_fy4_agri.exe Test_Raw_fy4_agri.F90)
# TARGET_LINK_LIBRARIES(Test_Raw_fy4_agri.exe DP_Satellite Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_Raw_fy3_4_hiras.exe Test_Raw_fy3_4_hiras.F90)
# TARGET_LINK_LIBRARIES(Test_Raw_fy3_4_hiras.exe DP_Satellite Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

