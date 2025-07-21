cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(PostProc)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

# -------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
SET(PostProc_SRC
	Export2SelDomain.F90
	PostProcTools.F90
	Export2HASCoordInSelDomain.F90
)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(PostProc STATIC ${PostProc_SRC})
TARGET_LINK_LIBRARIES(PostProc Geometry State State2NC)

# 2. Build test executables:
#SET(EXE1 "Test_Export2SelDomain.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_Export2SelDomain.F90)
# TARGET_LINK_LIBRARIES(${EXE1} PostProc)
