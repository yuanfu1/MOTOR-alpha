cmake_minimum_required(VERSION 3.12)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(AuxType)


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
ADD_LIBRARY(AuxType STATIC ${AuxType_SRC})
TARGET_LINK_LIBRARIES(AuxType Utility MPDD MPObs Geometry)

# 2. Build test executables:
# SET (EXE1 "Test_AuxType.exe")

# 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} Test_AuxType.F90)
# TARGET_LINK_LIBRARIES(${EXE1} AuxType Utility)
