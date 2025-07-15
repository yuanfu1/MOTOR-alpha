cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(C2O)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
INCLUDE_DIRECTORIES(${RTTOV_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${RTTOV_MOD_DIR})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Source codes to this build:
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(C2O STATIC ${C2O_SRC})
TARGET_LINK_LIBRARIES(C2O Utility M2ODirect State ObsSet M2CV OprRadar UV2W GeosBal Cumulus RTTOV UV2DirSpd OprGNSS )

ADD_DEPENDENCIES(C2O GeosBal)

# 2. Build test executables:
#SET (EXE1 "Test_C2O.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_C2O.F90)
# TARGET_LINK_LIBRARIES(${EXE1} GeosBal Geometry MPDD YAMLRead Field State ObsSet C2O Mock State2NC Cumulus RTTOV UV2DirSpd)
