cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(ObsRadar)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

# Add Slint libraries.
INCLUDE_DIRECTORIES(${SLINT_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
# Source codes to this build:
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(ObsRadar STATIC ${RadarRAW_SRC})
TARGET_LINK_LIBRARIES(ObsRadar ${NETCDF_LIBRARY} Utility ObsBase YAMLRead MPDD M2CV OprRadar Obs2State)

# 2. Build test executables:
#SET (EXE1 "Test_RadarRAW.exe")
#SET (EXE2 "Test_RefThinningOutput.exe")
#SET (EXE3 "Test_VelThinningOutput.exe")
#SET (EXE4 "Test_RefVar.exe")
#SET (EXE5 "Test_VarRadarRefNormal.exe")
#SET (EXE6 "Test_VelUVW2RWndVar.exe")
#SET (EXE7 "Test_VelUVW2RwndTLAD.exe")
#SET (EXE8 "Test_VelUV2RWndVar.exe")

# # 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_RadarRAW.F90)
# TARGET_LINK_LIBRARIES(${EXE1} ObsRadar Utility YAMLRead)

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE2} Test_RefThinningOutput.F90)
# TARGET_LINK_LIBRARIES(${EXE2} Geometry MPDD Utility C2O MPObs State ObsSet State2NC RMatrix ObsRadar ${MPI_Fortran_LIBRARIES} YAMLRead)

# # 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE3} Test_VelThinningOutput.F90)
# TARGET_LINK_LIBRARIES(${EXE3} Geometry MPDD Utility C2O MPObs State ObsSet State2NC RMatrix ObsRadar ${MPI_Fortran_LIBRARIES}
# RTTOV ${RTTOV_LIBRARY} YAMLRead IOGrapes IOERA5)

# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE4} Test_RefVar.F90)
# TARGET_LINK_LIBRARIES(${EXE4} MGOpts Geometry MPDD Utility Minimization
# JFunc BMatrix RMatrix C2O MPObs State ObsSet State2NC ObsRadar RTTOV ${MPI_Fortran_LIBRARIES}
# ${RTTOV_LIBRARY} YAMLRead)

# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE5} Test_VarRadarRefNormal.F90)
# TARGET_LINK_LIBRARIES(${EXE5} MGOpts Geometry MPDD IOGrapes Utility Minimization
# JFunc RMatrix BMatrix C2O MPObs State ObsSet State2NC ObsRadar RTTOV ${MPI_Fortran_LIBRARIES}
# ${RTTOV_LIBRARY} YAMLRead)

# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE6} Test_VelUVW2RWndVar.F90)
# TARGET_LINK_LIBRARIES(${EXE6} MGOpts Geometry MPDD IOGrapes Utility Minimization
# JFunc RMatrix BMatrix C2O MPObs State ObsSet State2NC ObsRadar RTTOV ${MPI_Fortran_LIBRARIES}
# ${RTTOV_LIBRARY} YAMLRead)

# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE7} Test_VelUVW2RwndTLAD.F90)
# TARGET_LINK_LIBRARIES(${EXE7} Geometry MPDD Utility C2O MPObs State ObsSet State2NC RMatrix ObsRadar RTTOV ${MPI_Fortran_LIBRARIES}
# ${RTTOV_LIBRARY} YAMLRead)

# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE8} Test_VelUV2RWndVar.F90)
# TARGET_LINK_LIBRARIES(${EXE8} MGOpts Geometry MPDD IOGrapes Utility Minimization RTTOV
# JFunc RMatrix BMatrix C2O MPObs State ObsSet State2NC ObsRadar ${MPI_Fortran_LIBRARIES}
# ${RTTOV_LIBRARY} YAMLRead)
#ADD_EXECUTABLE("Test_RADPreProcessing.exe" Test_RADPreProcessing.F90)
#TARGET_LINK_LIBRARIES("Test_RADPreProcessing.exe" Geometry MPDD Utility C2O MPObs State ObsSet State2NC RMatrix ObsRadar ${MPI_Fortran_LIBRARIES}
    # RTTOV ${RTTOV_LIBRARY} YAMLRead IOGrapes)

# # 4. Add tests:
# ADD_TEST(NAME "Test_VelUVW2RwndTLAD" COMMAND mpirun -np 4 Test_VelUVW2RwndTLAD.exe
# 	WORKING_DIRECTORY $ENV{BIN_DIR})

# # 5. Set test criteria:
# SET_TESTS_PROPERTIES(Test_VelUVW2RwndTLAD PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
# SET_TESTS_PROPERTIES(Test_VelUVW2RwndTLAD PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
