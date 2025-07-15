cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(IOERA5)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
# INCLUDE_DIRECTORIES("/opt/hpc/software/mathlib/netcdf/4.4.1/gcc/include/") 
# INCLUDE_DIRECTORIES(${SLINT_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})
# LINK_DIRECTORIES("/opt/hpc/software/mathlib/netcdf/4.4.1/gcc/lib")

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

# #-------------------------------------------------------------------------------
# # Generate lib file:
ADD_LIBRARY(IOERA5 STATIC ${IOERA5_SRC})
TARGET_LINK_LIBRARIES(IOERA5 PUBLIC 
GrapesIO Utility MPDD State IOModel YAMLRead InterpHP ModelCoupler M2CV ${MPI_Fortran_LIBRARIES} YAMLRead MGOpts)

# # 2. Build test executables:
#SET (EXE2 "Test_IOERA5.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE2} Test_IOERA5.F90)
# TARGET_LINK_LIBRARIES(${EXE2} Utility GrapesIO MPDD Geometry IOERA5 State State2NC ${MPI_Fortran_LIBRARIES} MGOpts IOERA5 IOGrapes)

