cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(MPDD)


#-------------------------------------------------------------------------------
# Add subdirectories to this build:

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

# SET (utility $ENV{COMMON_LIBS}/libutility.a)
#-------------------------------------------------------------------------------
# Generate lib file:
# ADD_LIBRARY(my_project_name STATIC ${my_project_name_SRC})

#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
ADD_LIBRARY(MPDD STATIC ${MPDD_SRC})
TARGET_LINK_LIBRARIES(MPDD Utility ${MPI_Fortran_LIBRARIES})


