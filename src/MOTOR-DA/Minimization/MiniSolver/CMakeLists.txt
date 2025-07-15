cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(Minimization)


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
ADD_LIBRARY(Minimization STATIC ${Minimization_SRC})
TARGET_LINK_LIBRARIES(Minimization Utility LBFGSB Field State ObsSet BMatrix C2O Mock State2NC JFunc)

# 2. Build test executables:
#SET (EXE1 "Test_MiniSolver.exe")
#SET (EXE2 "Test_SolverFRCG.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_MiniSolver.F90)
# TARGET_LINK_LIBRARIES(${EXE1} Minimization Field State ObsSet BMatrix C2O Mock State2NC JFunc RMatrix ${MPI_Fortran_LIBRARIES}
# RTTOV ${RTTOV_LIBRARY} lapack blas)

#ADD_EXECUTABLE(${EXE2} Test_SolverFRCG.F90)
# TARGET_LINK_LIBRARIES(${EXE2} Minimization Field State ObsSet BMatrix C2O Mock State2NC JFunc RMatrix ${MPI_Fortran_LIBRARIES}
# Minimization Field State ObsSet BMatrix C2O Mock State2NC JFunc RMatrix
# RTTOV ${RTTOV_LIBRARY} lapack blas)


# 4. Add tests: (MIGRATED TO UnitTests/MOTOR-DA/CMakeLists.txt)
# ADD_TEST(NAME Test_Solver_FRCG 
#         COMMAND mpirun -np 4 --oversubscribe Test_SolverFRCG.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})

# 5. Set test criteria: (MIGRATED TO UnitTests/MOTOR-DA/CMakeLists.txt)
# SET_TESTS_PROPERTIES(Test_Solver_FRCG PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
# SET_TESTS_PROPERTIES(Test_Solver_FRCG PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")

# 4. Add tests: (MIGRATED TO UnitTests/MOTOR-DA/CMakeLists.txt)
# ADD_TEST(NAME Test_Solver_LBFGS 
#         COMMAND mpirun -np 4 --oversubscribe Test_MiniSolver.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})

# 5. Set test criteria: (MIGRATED TO UnitTests/MOTOR-DA/CMakeLists.txt)
# SET_TESTS_PROPERTIES(Test_Solver_LBFGS PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
# SET_TESTS_PROPERTIES(Test_Solver_LBFGS PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
