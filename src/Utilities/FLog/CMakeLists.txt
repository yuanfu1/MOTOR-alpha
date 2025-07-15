cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(FLog)


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
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

# #-------------------------------------------------------------------------------
# # Generate lib file:
ADD_LIBRARY(FLog STATIC ${FLog_SRC})
TARGET_LINK_LIBRARIES(FLog YAMLRead)

# 2. Build test executables:
#SET (EXE1 "Test_FLog.exe")
# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_FLog.F90)
TARGET_LINK_LIBRARIES(${EXE1} FLog YAMLRead ${MPI_Fortran_LIBRARIES})

# 4. Add tests:
# ADD_TEST(NAME Test_FLog 
#         COMMAND mpirun -np 4 --oversubscribe Test_FLog.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})

# # 5. Set test criteria:
# SET_TESTS_PROPERTIES(Test_FLog PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
