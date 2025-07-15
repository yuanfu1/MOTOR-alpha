cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(ObsUtilities)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
# INCLUDE_DIRECTORIES(/Users/yjzx-xieyuanfu/developments/mods)
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})
# LINK_DIRECTORIES(${NETCDF_LIBRARY})

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
ADD_LIBRARY(ObsUtilities STATIC ${ObsUtilities_SRC})
TARGET_LINK_LIBRARIES(ObsUtilities ${NETCDF_LIBRARY} FLog Utility M2ODirect State ObsSet ObsBase fortran_stdlib ${SLINT_LIBRARY})

# # 2. Build test executables:
# SET (EXE1 "Test_ObsConcat.exe")
# SET (EXE2 "Test_netCDFReader.exe")

# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} Test_ObsConcat.F90)
# ADD_EXECUTABLE(${EXE2} Test_netCDFReader.F90)
# TARGET_LINK_LIBRARIES(${EXE1} MGOpts MPObs MPDD Geometry Utility Minimization Field State ObsSet ObsBase ObsConvention ObsUtilities BMatrix C2O Mock State2NC JFunc MGOpts Obs2State RTTOV ${RTTOV_LIBRARY} lapack blas)
# TARGET_LINK_LIBRARIES(${EXE2} Utility ObsUtilities)
