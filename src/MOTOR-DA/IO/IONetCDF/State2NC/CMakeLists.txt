cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(State2NC)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(State2NC STATIC ${State2NC_SRC})
TARGET_LINK_LIBRARIES(State2NC Utility State)

# 2. Build test executables:
#SET (EXE2 "Test_State2NC.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE2} Test_State2NC.F90)
# TARGET_LINK_LIBRARIES(${EXE2} Utility MPDD Geometry State State2NC ${MPI_Fortran_LIBRARIES} YAMLRead )