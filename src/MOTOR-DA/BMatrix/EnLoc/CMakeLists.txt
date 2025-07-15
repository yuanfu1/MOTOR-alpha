cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(EnLoc)

# Include the common template:
IF (IS_ABSOLUTE $ENV{CTEMPLATES})
  MESSAGE(STATUS "CTEMPLATES path is set!")
ELSE ($ENV{CTEMPLATES} LESS 1)
  MESSAGE(FATAL_ERROR "Env CTEMPLATES is not set, please set and rerun!")
ENDIF (IS_ABSOLUTE $ENV{CTEMPLATES})
include($ENV{CTEMPLATES}/cmakeCommon.cmake)

#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
# INCLUDE_DIRECTORIES(/Users/yjzx-xieyuanfu/developments/mods)
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})
LINK_DIRECTORIES(${NETCDF_LIBRARY})

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
ADD_LIBRARY(EnLoc STATIC ${EnLoc_SRC})
TARGET_LINK_LIBRARIES(EnLoc PUBLIC ${NETCDF_LIBRARY} Utility State YAMLRead blas eispack linpack lapack)
# TARGET_LINK_LIBRARIES(EnLoc PUBLIC ${NETCDF_LIBRARY} Utility State YAMLRead eispack linpack lapack)

# # 2. Build test executables:
# SET (EXE1 "Test_EnLoc_model.exe")
# SET (EXE2 "Test_EnLoc_state.exe")
# SET (EXE3 "Test_EnLoc_Multigrid.exe")
# SET (EXE4 "Test_QRrecover.exe")
# SET (EXE5 "Test_QR_Multigrid.exe")
# SET (EXE6 "Test_QRdecomp_invR.exe")
# SET (EXE7 "Test_QRdecomp_invRv2.exe")
# SET (EXE8 "Test_BKErr_model_FWD_ADJ_Binv.exe")

# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} Test_EnLoc_model.F90)
# TARGET_LINK_LIBRARIES(${EXE1} EnLoc eispack linpack lapack blas)

# ADD_EXECUTABLE(${EXE2} Test_EnLoc_state.F90)
# TARGET_LINK_LIBRARIES(${EXE2} EnLoc MPDD Geometry IOWRF YAMLRead State2NC ${MPI_Fortran_LIBRARIES} eispack linpack lapack blas ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_C}
# ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

# ADD_EXECUTABLE(${EXE3} Test_EnLoc_Multigrid.F90)
# TARGET_LINK_LIBRARIES(${EXE3} MGOpts Geometry IOWRF State2NC EnLoc eispack linpack lapack blas)

# ADD_EXECUTABLE(${EXE4} Test_QRrecover.F90)
# TARGET_LINK_LIBRARIES(${EXE4} EnLoc eispack linpack lapack blas)

# ADD_EXECUTABLE(${EXE5} Test_QR_Multigrid.F90)
# TARGET_LINK_LIBRARIES(${EXE5} EnLoc MPDD Geometry IOWRF YAMLRead State2NC MGOpts Field ${MPI_Fortran_LIBRARIES} eispack linpack lapack blas ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_C}
# ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

# ADD_EXECUTABLE(${EXE6} Test_QRdecomp_invR.F90)
# TARGET_LINK_LIBRARIES(${EXE6} EnLoc eispack linpack lapack blas)

# ADD_EXECUTABLE(${EXE7} Test_QRdecomp_invRv2.F90)
# TARGET_LINK_LIBRARIES(${EXE7} EnLoc eispack linpack lapack blas)


# ADD_EXECUTABLE(${EXE8} Test_BKErr_model_FWD_ADJ_Binv.F90)
# TARGET_LINK_LIBRARIES(${EXE8} EnLoc eispack linpack lapack blas)
