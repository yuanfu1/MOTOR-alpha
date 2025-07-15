cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(IOModel)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
# INCLUDE_DIRECTORIES("/opt/hpc/software/mathlib/netcdf/4.4.1/gcc/include/")

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})
# LINK_DIRECTORIES("/opt/hpc/software/mathlib/netcdf/4.4.1/gcc/lib")

SET (CMAKE_Fortran_FLAGS "-cpp")
#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

# #-------------------------------------------------------------------------------
# # Generate lib file:
 ADD_LIBRARY(IOModel STATIC ${IOModel_SRC})
 TARGET_LINK_LIBRARIES(IOModel Utility MPDD State ${MPI_Fortran_LIBRARIES})

