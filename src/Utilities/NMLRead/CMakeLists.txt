cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(NMLRead)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

# SET (CMAKE_Fortran_FLAGS "-cpp -g")
#-------------------------------------------------------------------------------
# Source codes to this build:

# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

# #-------------------------------------------------------------------------------
# # Generate lib file:
ADD_LIBRARY(NMLRead STATIC ${NMLRead_SRC})
TARGET_LINK_LIBRARIES(NMLRead Utility)
 
# 2. Build test executables:
#SET (EXE1 "Test_NMLRead.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_NMLRead.F90)
# TARGET_LINK_LIBRARIES(${EXE1} Utility NMLRead)

# 4. Add tests: (MIGRATED TO UnitTests/Utilities/CMakeLists.txt)
# ADD_TEST(NAME Test_NMLRead 
#         COMMAND Test_NMLRead.exe
#         WORKING_DIRECTORY ${BIN_DIR})

# 5. Set test criteria: (MIGRATED TO UnitTests/Utilities/CMakeLists.txt)
# SET_TESTS_PROPERTIES(Test_NMLRead PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed")
# SET_TESTS_PROPERTIES(Test_NMLRead PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
