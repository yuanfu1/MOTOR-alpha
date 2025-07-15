cmake_minimum_required(VERSION 3.12)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(State_Xm)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(State STATIC ${State_SRC})
TARGET_LINK_LIBRARIES(State Field MPDD Geometry YAMLRead Utility AuxType ${MPI_Fortran_LIBRARIES})

# 2. Build test executables:
#SET (EXE1 "Test_State.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_State.F90)
TARGET_LINK_LIBRARIES(${EXE1} State MPDD Geometry Utility ${MPI_Fortran_LIBRARIES} YAMLRead MGGen fortran-yaml-cpp)

# 4. Add tests: (MIGRATED TO UnitTests/MOTOR-DA/CMakeLists.txt)
# ADD_TEST(NAME "Test_States" 
#         COMMAND mpirun -np 4 --oversubscribe Test_State.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})

# 5. Set test criteria: (MIGRATED TO UnitTests/MOTOR-DA/CMakeLists.txt)
# SET_TESTS_PROPERTIES(Test_States PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
# SET_TESTS_PROPERTIES(Test_States PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
