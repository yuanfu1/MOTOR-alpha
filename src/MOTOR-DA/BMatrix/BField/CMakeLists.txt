cmake_minimum_required(VERSION 3.12)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(BField)


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
ADD_LIBRARY(BField STATIC ${BField_SRC})
TARGET_LINK_LIBRARIES(BField Utility EnLoc blas)

# 2. Build test executables:
#SET (EXE1 "Test_BField.exe")
#SET (EXE2 "Test_EnLoc_ADcheck.exe")

# 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} Test_BField.F90)
# TARGET_LINK_LIBRARIES(${EXE1} BField Utility)

#ADD_EXECUTABLE(${EXE2} Test_EnLoc_ADcheck.F90)
# TARGET_LINK_LIBRARIES(${EXE2} Utility EnLoc linpack lapack blas M2CV)
