cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(CFD)


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
SET (CFD_SRC	
	F3_MtmW.F90
	)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(CFD STATIC ${CFD_SRC})
TARGET_LINK_LIBRARIES(CFD PUBLIC State Field Geometry Utility TransBase lapack blas YAMLRead)

