cmake_minimum_required(VERSION 3.12)

# Build MOTOR-QC Unit Tests
PROJECT(MOTOR-QC_UnitTests)

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
# Build test executables for MOTOR-QC
#-------------------------------------------------------------------------------


# MOTOR-QC Test Executables

# Observation Tests
ADD_EXECUTABLE(Test_ObsSuper.exe Test_ObsSuper.F90)
TARGET_LINK_LIBRARIES(Test_ObsSuper.exe Application Minimization Geometry BMatrix C2O
ObsConvention Geometry MPDD YAMLRead Utility State Field ObsSet ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_obsControl.exe Test_obsControl.F90)
TARGET_LINK_LIBRARIES(Test_obsControl.exe ObsConvention Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_ObsRead.exe Test_ObsRead.F90)
TARGET_LINK_LIBRARIES(Test_ObsRead.exe ObsConvention Geometry MPDD YAMLRead Utility ObsUtilities ${MPI_Fortran_LIBRARIES})

# Radar Tests
# ADD_EXECUTABLE(Test_RefVar.exe Test_RefVar.F90)
# TARGET_LINK_LIBRARIES(Test_RefVar.exe ObsRadar Geometry MPDD YAMLRead Utility State Field ObsSet ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_RefThinningOutput.exe Test_RefThinningOutput.F90)
TARGET_LINK_LIBRARIES(Test_RefThinningOutput.exe ObsRadar Geometry MPDD YAMLRead Utility State Field ObsSet ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_VelUV2RWndVar.exe Test_VelUV2RWndVar.F90)
TARGET_LINK_LIBRARIES(Test_VelUV2RWndVar.exe Application Minimization BMatrix C2O
ObsRadar Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_RadarRAW.exe Test_RadarRAW.F90)
TARGET_LINK_LIBRARIES(Test_RadarRAW.exe ObsRadar Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_RADPreProcessing.exe Test_RADPreProcessing.F90)
TARGET_LINK_LIBRARIES(Test_RADPreProcessing.exe ObsRadar Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_VelThinningOutput.exe Test_VelThinningOutput.F90)
TARGET_LINK_LIBRARIES(Test_VelThinningOutput.exe Application Minimization BMatrix C2O
ObsRadar Geometry MPDD YAMLRead Utility State Field ObsSet ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_VelUVW2RwndTLAD.exe Test_VelUVW2RwndTLAD.F90)
TARGET_LINK_LIBRARIES(Test_VelUVW2RwndTLAD.exe Application Minimization BMatrix C2O
ObsRadar Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_VarRadarRefNormal.exe Test_VarRadarRefNormal.F90)
# TARGET_LINK_LIBRARIES(Test_VarRadarRefNormal.exe Application Minimization BMatrix C2O
# ObsRadar Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_VelUVW2RWndVar.exe Test_VelUVW2RWndVar.F90)
TARGET_LINK_LIBRARIES(Test_VelUVW2RWndVar.exe Application Minimization BMatrix C2O
ObsRadar Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# GNSS Tests
ADD_EXECUTABLE(Test_ObsGNSS.exe Test_ObsGNSS.F90)
TARGET_LINK_LIBRARIES(Test_ObsGNSS.exe ObsGNSS Geometry MPDD YAMLRead Utility State Field ObsSet ${MPI_Fortran_LIBRARIES})

# Utilities Tests
ADD_EXECUTABLE(Test_netCDFReader.exe Test_netCDFReader.F90)
TARGET_LINK_LIBRARIES(Test_netCDFReader.exe ObsUtilities Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES} ${NETCDF_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_ObsConcat.exe Test_ObsConcat.F90)
TARGET_LINK_LIBRARIES(Test_ObsConcat.exe Application Minimization BMatrix C2O
ObsUtilities Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

