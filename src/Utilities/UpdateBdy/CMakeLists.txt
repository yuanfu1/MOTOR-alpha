cmake_minimum_required(VERSION 3.5)

PROJECT(UpdateBdy)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

SET (CMAKE_Fortran_FLAGS "-DUSE_MPI -fpic -O3 -cpp -ffree-line-length-none -D_REAL8_ -fopenmp -DDEBUG")
#-------------------------------------------------------------------------------
# Source codes to this build:

# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

# #-------------------------------------------------------------------------------
# # Generate lib file:
ADD_LIBRARY(UpdateBdy STATIC ${UpdateBdy_SRC})
TARGET_LINK_LIBRARIES(UpdateBdy GrapesIO Utility NMLRead YAMLRead ${MPI_Fortran_LIBRARIES})
 
# 2. Build test executables:
#SET (EXE1 "Test_updatebdy.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_updatebdy.F90)
TARGET_LINK_LIBRARIES(${EXE1} UpdateBdy GrapesIO Utility NMLRead YAMLRead ${MPI_Fortran_LIBRARIES})
