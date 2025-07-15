cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(mgGen)


#-------------------------------------------------------------------------------
# Add subdirectories to this build:
#ADD_SUBDIRECTORY(grids)
#ADD_SUBDIRECTORY(utils)
#ADD_SUBDIRECTORY(initialConditions)

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (FilesUsed.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
# ADD_LIBRARY(my_project_name STATIC ${my_project_name_SRC})
ADD_LIBRARY(MGGen STATIC ${MGGen_SRC})
TARGET_LINK_LIBRARIES(MGGen YAMLRead)
# ADD_LIBRARY(utility STATIC IMPORTED)
# SET_TARGET_PROPERTIESADD_LIBRARY(icos STATIC IMPORTED)
# ADD_LIBRARY(itmesh STATIC IMPORTED)

#-------------------------------------------------------------------------------
# Add dependencies:
# ADD_DEPENDENCIES(gzm gzGrids)

#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
IF (IS_ABSOLUTE $ENV{COMMON_MODS})
  MESSAGE(STATUS "COMMON_MODS path is set!")
ELSE ($ENV{COMMON_MODS} LESS 1)
  MESSAGE(FATAL_ERROR "Env COMMON_MODS is not set, please set and rerun!")
ENDIF (IS_ABSOLUTE $ENV{COMMON_MODS})
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

# Setting FLAG to turn on and off MACROs in F90 codes:
SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG}")
#SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -DMID_EDGE")
MESSAGE(STATUS "CMAKE_Fortran_FLAGS: " ${CMAKE_Fortran_FLAGS})

#-------------------------------------------------------------------------------
IF (IS_ABSOLUTE $ENV{COMMON_LIBS})
  MESSAGE(STATUS "COMMON_LIBS path is set!")
ELSE ($ENV{CTEMPLATES} LESS 1)
  MESSAGE(FATAL_ERROR "Env COMMON_LIBS is not set, please set and rerun!")
ENDIF (IS_ABSOLUTE $ENV{COMMON_LIBS})
#-------------------------------------------------------------------------------
# For libraries saved in a place, import them:
# ADD_LIBRARY(libname STATIC IMPORTED)
# SET_TARGET_PROPERTIES(libname PROPERTIES IMPORTED_LOCATION <fullpath/libname.a>)
# Set up libraries
# SET_TARGET_PROPERTIES(utility PROPERTIES IMPORTED_LOCATION $ENV{COMMON_LIBS}/libutility.a)
#SET_TARGET_PROPERTIES(icos PROPERTIES IMPORTED_LOCATION $ENV{COMMON_LIBS}/libicos.a)
#SET_TARGET_PROPERTIES(itmesh PROPERTIES IMPORTED_LOCATION $ENV{COMMON_LIBS}/libitmesh.a)

#-------------------------------------------------------------------------------
#  The following is to set up ctest:
#-------------------------------------------------------------------------------
# 1. Copy test scripts: for multiple testings, simply add more of the follow
# lines
# FILE(COPY <test_script> DESTINATION ${PROJECT_BINARY_DIR})

# 2. Build test executables:
#SET (EXE1 "test_mgGen.exe")
#SET (EXE2 "test_icos.exe")
#SET (EXE3 "TestEdgeSetting.exe")
# SET (utility $ENV{COMMON_LIBS}/libutility.a)

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} unitTestLatLon.F90 test_mgGen.F90)
# ADD_EXECUTABLE(${EXE2} test_icos.F90)
# ADD_EXECUTABLE(${EXE3} TestEdgeSetting.F90)

#-------------------------------------------------------------------------------
#  For libraries saved in the same directory of COMMON_LIBS:
#  To point to another lib like -L in cmake:
#-------------------------------------------------------------------------------
#  LINK_DIRECTORIES(dir1 dir2 ...)
#-------------------------------------------------------------------------------
#  To link another lib like -l in cmake:
#-------------------------------------------------------------------------------
#  TARGET_LINK_LIBRARIES(target <list of libnames>)
# TARGET_LINK_LIBRARIES(${EXE1} MGGen Utility)
# TARGET_LINK_LIBRARIES(${EXE2} MGGen Utility)
# TARGET_LINK_LIBRARIES(${EXE3} MGGen Utility)

# # 4. Add tests:
# ADD_TEST("unit_test1" test_mgGen.exe)
# ADD_TEST("unit_test2" test_icos.exe)

# # 5. Set test criteria:
# SET_TESTS_PROPERTIES(unit_test1 PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
# SET_TESTS_PROPERTIES(unit_test1 PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
# SET_TESTS_PROPERTIES(unit_test2 PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
# SET_TESTS_PROPERTIES(unit_test2 PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
