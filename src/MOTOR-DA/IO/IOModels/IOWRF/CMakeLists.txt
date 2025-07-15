cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(IOWRF)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
# INCLUDE_DIRECTORIES("/opt/hpc/software/mathlib/netcdf/4.4.1/gcc/include/") 
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})
# LINK_DIRECTORIES("/opt/hpc/software/mathlib/netcdf/4.4.1/gcc/lib")

SET (CMAKE_Fortran_FLAGS "-cpp -g")
#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

# #-------------------------------------------------------------------------------
# # Generate lib file:
ADD_LIBRARY(IOWRF STATIC ${IOWRF_SRC})
TARGET_LINK_LIBRARIES(IOWRF Utility MPDD State IOModel YAMLRead InterpHP ModelCoupler ${MPI_Fortran_LIBRARIES})

# 2. Build test executables:
#SET (EXE1 "Test_IOWRF.exe")
#SET (EXE2 "Test_IOWRF_Ens.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_IOWRF.F90)
# TARGET_LINK_LIBRARIES(${EXE1} Utility MPDD Geometry IOWRF State YAMLRead State2NC ${MPI_Fortran_LIBRARIES} )
#ADD_EXECUTABLE(${EXE2} Test_IOWRF_Ens.F90)
# TARGET_LINK_LIBRARIES(${EXE2} Utility MPDD Geometry IOWRF State YAMLRead State2NC ${MPI_Fortran_LIBRARIES} )
