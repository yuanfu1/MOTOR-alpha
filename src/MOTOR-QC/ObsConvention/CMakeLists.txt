cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(ObsConvention)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
# INCLUDE_DIRECTORIES(/Users/yjzx-xieyuanfu/developments/mods)
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})
#LINK_DIRECTORIES(${NETCDF_LIBRARY})

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
ADD_LIBRARY(ObsConvention STATIC ${ObsConvention_SRC})
TARGET_LINK_LIBRARIES(ObsConvention ${NETCDF_LIBRARY} FLog Utility M2ODirect State ObsSet ${SLINT_LIBRARY} ObsBase linpack)

# 2. Build test executables:
#SET (EXE1 "Test_ObsRead.exe")
#SET (EXE2 "Test_ObsSuper.exe")
# SET (EXE3 "Test_ObsRead_VWPW.exe")
# SET (EXE4 "Test_ObsSuper_VWPW.exe")
#SET (EXE5 "Test_ObsRead_LBUOY.exe")
#SET (EXE6 "Test_ObsRead_SING.exe")
#SET (EXE7 "Test_ObsRead_GWST.exe")

#SET(EXE8 "Test_ObsControl.exe")

# 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} Test_ObsRead.F90)
# TARGET_LINK_LIBRARIES(${EXE1} Geometry MPDD YAMLRead Field State ObsSet ObsBase ObsConvention C2O Mock State2NC Obs2State MGOpts lapack blas)

# ADD_EXECUTABLE(${EXE2} Test_ObsSuper.F90)
# TARGET_LINK_LIBRARIES(${EXE2} MGOpts MPObs MPDD Geometry Utility Minimization 
# Field State ObsSet ObsBase ObsConvention BMatrix C2O Mock State2NC JFunc MGOpts Obs2State
# RTTOV ${RTTOV_LIBRARY} lapack blas)

# ADD_EXECUTABLE(${EXE3} Test_ObsRead_VWPW.F90)
# TARGET_LINK_LIBRARIES(${EXE3} Geometry MPDD YAMLRead Field State ObsSet ObsBase ObsConvention C2O Mock State2NC Obs2State MGOpts)

# ADD_EXECUTABLE(${EXE4} Test_ObsSuper_VWPW.F90)
# TARGET_LINK_LIBRARIES(${EXE4} MGOpts MPObs MPDD Geometry Utility Minimization 
# Field State ObsSet ObsBase ObsConvention BMatrix C2O Mock State2NC JFunc MGOpts Obs2State
# RTTOV ${RTTOV_LIBRARY} )

# ADD_EXECUTABLE(${EXE5} Test_ObsRead_LBUOY.F90)
# TARGET_LINK_LIBRARIES(${EXE5} Geometry MPDD YAMLRead Field State ObsSet ObsBase ObsConvention C2O Mock State2NC Obs2State MGOpts lapack blas)

# ADD_EXECUTABLE(${EXE6} Test_ObsRead_SING.F90)
# TARGET_LINK_LIBRARIES(${EXE6} Geometry MPDD YAMLRead Field State ObsSet ObsBase ObsConvention C2O Mock State2NC Obs2State MGOpts lapack blas)

# ADD_EXECUTABLE(${EXE7} Test_ObsRead_GWST.F90)
# TARGET_LINK_LIBRARIES(${EXE7} Geometry MPDD YAMLRead Field State ObsSet ObsBase ObsConvention C2O Mock State2NC Obs2State MGOpts lapack blas)

# ADD_EXECUTABLE(${EXE8} Test_obsControl.F90)
# TARGET_LINK_LIBRARIES(${EXE8} Geometry MPDD YAMLRead Field State ObsSet ObsBase ObsConvention C2O Mock State2NC Obs2State MGOpts lapack blas)
