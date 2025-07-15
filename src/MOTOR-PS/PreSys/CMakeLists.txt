cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(PreSys)


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

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(PreSys STATIC ${PreSys_SRC})
TARGET_LINK_LIBRARIES(PreSys MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} zgrid )

# 2. Build test executables:
#SET (EXE1 "Test_interpvalue.exe")
#SET (EXE2 "Test_bdyidx.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_InterpValue.F90)
TARGET_LINK_LIBRARIES(${EXE1} PreSys MPDD Geometry Utility FLog ${MPI_Fortran_LIBRARIES} zgrid)

#ADD_EXECUTABLE(${EXE2} Test_bdyidx.F90)
TARGET_LINK_LIBRARIES(${EXE2} PreSys MPDD Geometry Utility FLog ${MPI_Fortran_LIBRARIES} YAMLRead)

# # 4. Add tests:
# ADD_TEST(NAME "Test_PossionSolver" 
#         COMMAND mpirun -np 4 --oversubscribe Test_PS.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})

# # 5. Set test criteria:
# SET_TESTS_PROPERTIES(Test_PossionSolver PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
# SET_TESTS_PROPERTIES(Test_PossionSolver PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
