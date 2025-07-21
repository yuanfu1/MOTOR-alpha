cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(OprRadar)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Source codes to this build:
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(OprRadar STATIC ${OprRadar_SRC})
TARGET_LINK_LIBRARIES(OprRadar ObsSet State Field Geometry Utility TransBase)

# 2. Build test executables:
#SET (EXE1 "Test_RadarVel.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_RadarVel.F90)
# TARGET_LINK_LIBRARIES(${EXE1} OprRadar)