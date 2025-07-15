cmake_minimum_required(VERSION 3.12)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(BMatrix)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(BMatrix STATIC ${BMatrix_SRC})
TARGET_LINK_LIBRARIES(BMatrix PUBLIC Utility AuxType BField M2CV linpack blas)

# 2. Build test executables:
#SET (EXE1 "Test_BMatrix.exe")

# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} Test_BMatrix.F90)
# TARGET_LINK_LIBRARIES(${EXE1} BMatrix Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES} linpack lapack blas)

#ADD_EXECUTABLE(Test_BMatTLAD.exe Test_BMatTLAD.F90)
#TARGET_LINK_LIBRARIES(Test_BMatTLAD.exe IOGrapes BMatrix Geometry MPDD YAMLRead Utility State Field EnLoc ${MPI_Fortran_LIBRARIES})
