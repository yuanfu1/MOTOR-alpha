cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(Application)

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
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(Application STATIC ${Application_SRC})

TARGET_LINK_LIBRARIES(Application PUBLIC ${NETCDF_LIBRARY} ${MPI_Fortran_LIBRARIES} IOGrapes IOWRF ObsConvention ObsUtilities C2O C2M RMatrix GeosBal BMatrix Mock BkgdMG ObsMG
	MGOpts Obs2State Minimization JFunc ObsRadar ObsSatellite ObsSatob PostProc lapack blas ObsGNSS IOERA5)


add_subdirectory(ConventionalAnalysis)

# 1. Build test executables:
SET (EXE7  "App_3DVAR_MultiGrid_WRF_MultiObs.exe")
SET (EXE8  "App_CMA_GD.exe")
SET (EXE9  "App_CMA_GD_EnLoc.exe")
SET (EXE10 "App_CalculateBEC_EnLoc.exe")
SET (EXE11 "App_SurfaceAnalysis.exe")
SET (EXE12 "App_CMA_GD_Opr.exe")
SET (EXE13 "App_3DAnalysis.exe")
SET (EXE14 "App_Get_BcgAtObs.exe")
SET (EXE15 "App_CMA_GD_EnLoc_BEC.exe")
SET (EXE16 "App_CMA_GD_EnLoc_Modelvar.exe")
SET (EXE17 "Driver_App_4DVAR.exe")
# SET (EXE18 "App_CMA_GD_HybInc.exe")

# 2. Add codes to the executables:


# # EXE8
ADD_EXECUTABLE(${EXE8} App_CMA_GD.F90)
TARGET_LINK_LIBRARIES(${EXE8} Application Minimization Field Geometry State ObsSet BMatrix C2O
	Mock State2NC JFunc IOWRF ObsBase RMatrix ObsConvention ObsUtilities MGOpts Utility Obs2State
	ObsRadar RTTOV ${RTTOV_LIBRARY} YAMLRead ${MPI_Fortran_LIBRARIES} ObsSatellite PostProc IOECM)

# EXE9
ADD_EXECUTABLE(${EXE9} App_CMA_GD_EnLoc.F90)
TARGET_LINK_LIBRARIES(${EXE9} Application Minimization Field Geometry State ObsSet BMatrix C2O
	Mock State2NC JFunc IOWRF ObsBase RMatrix ObsConvention ObsUtilities MGOpts Utility Obs2State
	ObsRadar RTTOV ${RTTOV_LIBRARY} YAMLRead ${MPI_Fortran_LIBRARIES} ObsSatellite PostProc lapack blas IOECM)

 # App_SurfaceAnalysis
ADD_EXECUTABLE(${EXE11} App_SurfaceAnalysis.F90)
TARGET_LINK_LIBRARIES(${EXE11} Application Minimization Field Geometry State ObsSet BMatrix C2O
	Mock State2NC JFunc IOWRF ObsBase RMatrix ObsConvention ObsUtilities MGOpts Utility Obs2State
	ObsRadar RTTOV ${RTTOV_LIBRARY} YAMLRead ${MPI_Fortran_LIBRARIES} ObsSatellite PostProc lapack blas IOECM)

ADD_EXECUTABLE(${EXE12} App_CMA_GD_Opr.F90)
TARGET_LINK_LIBRARIES(${EXE12} Application Minimization Field Geometry State ObsSet BMatrix C2O
	Mock State2NC JFunc IOWRF ObsBase RMatrix ObsConvention ObsUtilities MGOpts Utility Obs2State
	ObsRadar RTTOV ${RTTOV_LIBRARY} YAMLRead ${MPI_Fortran_LIBRARIES} ObsSatellite PostProc lapack blas IOECM)

ADD_EXECUTABLE(${EXE13} App_3DAnalysis.F90)
TARGET_LINK_LIBRARIES(${EXE13} Application Minimization Field Geometry State ObsSet BMatrix C2O
	Mock State2NC JFunc IOWRF ObsBase RMatrix ObsConvention ObsUtilities MGOpts Utility Obs2State
	ObsRadar RTTOV ${RTTOV_LIBRARY} YAMLRead ${MPI_Fortran_LIBRARIES} ObsSatellite PostProc lapack blas IOECM)

ADD_EXECUTABLE(${EXE14} App_Get_BcgAtObs.F90)
TARGET_LINK_LIBRARIES(${EXE14} Application Minimization Field Geometry State ObsSet BMatrix C2O
	Mock State2NC JFunc IOWRF ObsBase RMatrix ObsConvention ObsUtilities MGOpts Utility Obs2State
	ObsRadar RTTOV ${RTTOV_LIBRARY} YAMLRead ${MPI_Fortran_LIBRARIES} ObsSatellite PostProc lapack blas IOECM)

# 3. Add tests: Temporarily turned it off waiting for how to fit background fields
#ADD_TEST(NAME "App_CMA_GD" COMMAND mpirun -np 4 App_CMA_GD.exe $ENV{STATIC_DIR}/Application/App_CMA_GD.yaml
#	WORKING_DIRECTORY $ENV{BIN_DIR})

# 4. Set test criteria:
#SET_TESTS_PROPERTIES(App_CMA_GD PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
#SET_TESTS_PROPERTIES(App_CMA_GD PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")

# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE8} App_3DVAR_MultiGrid_WRF_EnLoc.F90)
# TARGET_LINK_LIBRARIES(${EXE8} EnLoc Minimization Field Geometry State ObsSet BMatrix C2O Mock State2NC JFunc IOWRF RMatrix ObsBase ObsConvention MGOpts Utility Obs2State)

# EXE15
ADD_EXECUTABLE(${EXE15} App_CMA_GD_EnLoc_BEC.F90)
TARGET_LINK_LIBRARIES(${EXE15} Application Minimization Field Geometry State ObsSet BMatrix C2O
	Mock State2NC JFunc IOWRF ObsBase RMatrix ObsConvention ObsUtilities MGOpts Utility Obs2State
	ObsRadar RTTOV ${RTTOV_LIBRARY} YAMLRead ${MPI_Fortran_LIBRARIES} ObsSatellite PostProc lapack blas IOECM)

# EXE16
ADD_EXECUTABLE(${EXE16} App_CMA_GD_EnLoc_Modelvar.F90)
TARGET_LINK_LIBRARIES(${EXE16} Application Minimization Field Geometry State ObsSet BMatrix C2O
	Mock State2NC JFunc IOWRF ObsBase RMatrix ObsConvention ObsUtilities MGOpts Utility Obs2State
	ObsRadar RTTOV ${RTTOV_LIBRARY} YAMLRead ${MPI_Fortran_LIBRARIES} ObsSatellite PostProc lapack blas IOECM)

# EXE17
ADD_EXECUTABLE(${EXE17} Driver_App_4DVAR.F90)
TARGET_LINK_LIBRARIES(${EXE17} Application Minimization Field Geometry State ObsSet BkgdMG ObsMG BMatrix C2M
	Mock State2NC JFunc IOWRF ObsBase RMatrix ObsConvention ObsUtilities MGOpts Utility Obs2State
	ObsRadar RTTOV ${RTTOV_LIBRARY} YAMLRead ${MPI_Fortran_LIBRARIES} ObsSatellite PostProc lapack blas IOECM)

# EXE18
# ADD_EXECUTABLE(${EXE18} Test_ImmersedBoundary.F90)
# TARGET_LINK_LIBRARIES(${EXE18} Utility Geometry BkgdMG IOGrapes IOWRF State2NC YAMLRead ${MPI_Fortran_LIBRARIES})

