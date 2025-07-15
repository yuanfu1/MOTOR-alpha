cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
PROJECT(laps_lib)

#-------------------------------------------------------------------------------
# Source codes to this build:
include (srcFiles.txt)

#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES(../../include)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(laps_lib STATIC ${laps_lib_SRC})
