cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(yourModule)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# 1. Generate lib file:
ADD_LIBRARY(yourModule STATIC ${yourModule_SRC})
TARGET_LINK_LIBRARIES(yourModule Utility)

# 2. Build test executables:
SET (EXE1 "Test_YM.exe")

# 3. Add codes to the executables:
ADD_EXECUTABLE(${EXE1} Test_YM.F90)
TARGET_LINK_LIBRARIES(${EXE1} yourModule YAMLRead Utility)
