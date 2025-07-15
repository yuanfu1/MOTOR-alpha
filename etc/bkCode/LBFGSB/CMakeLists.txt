cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
PROJECT(LBFGSB)

# Include the common template:
IF (IS_ABSOLUTE $ENV{CTEMPLATES})
  MESSAGE(STATUS "CTEMPLATES path is set!")
ELSE ($ENV{CTEMPLATES} LESS 1)
  MESSAGE(FATAL_ERROR "Env CTEMPLATES is not set, please set and rerun!")
ENDIF (IS_ABSOLUTE $ENV{CTEMPLATES})
include($ENV{CTEMPLATES}/cmakeCommon.cmake)

# Add netcdf libraries.
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})
LINK_DIRECTORIES(${NETCDF_LIBRARY})

#-------------------------------------------------------------------------------
# Source codes to this build:
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(LBFGSB STATIC ${LBFGSB_SRC})

# 2. Build test executables:
SET (EXE1 "drver1.exe")
# SET (EXE2 "test_rotate.exe")

# 3. Add code to the executables:
ADD_EXECUTABLE(${EXE1} driver1.f90)
TARGET_LINK_LIBRARIES(${EXE1} LBFGSB)
