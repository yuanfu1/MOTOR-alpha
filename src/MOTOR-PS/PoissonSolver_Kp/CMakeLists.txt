cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(PoissonSolver)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
#  Check the STATIC_DIR:
IF (IS_ABSOLUTE $ENV{STATIC_DIR}) 
        MESSAGE(STATUS "env STATIC_DIR is set!") 
ELSE ($ENV{COMMON_LIBS} LESS 1)
        MESSAGE(FATAL_ERROR "Env STATIC is not set, please set and rerun!")
ENDIF (IS_ABSOLUTE $ENV{STATIC_DIR})

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(PoissonSolver STATIC ${PS_SRC})
TARGET_LINK_LIBRARIES(PoissonSolver MGOpts MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

# 2. Build test executables:
#SET(EXE1 "Test_PSForward.exe")
SET (EXE2 "jacobi.exe")
#SET (EXE3 "TestMixedSolver.exe")
SET (EXE4 "PlotAState.exe")
#SET (EXE5 "TestAnalytics.exe")
#SET (EXE6 "TestPsiChi2UV.exe")
#SET (EXE7 "Test_PSAdjoint.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_PSForward.F90)
ADD_EXECUTABLE(${EXE2} jacobi_test.f90)
#ADD_EXECUTABLE(${EXE3} Analytics.F90 TestMixedSolver.F90)
ADD_EXECUTABLE(${EXE4} PlotAState.F90)
# ADD_EXECUTABLE(${EXE5} TestAnalytics.F90)
# ADD_EXECUTABLE(${EXE6} TestPsiChi2UV.F90)
#ADD_EXECUTABLE(${EXE7} Test_PSAdjoint.F90)

# TARGET_LINK_LIBRARIES(${EXE1} PoissonSolver MPDD Geometry Utility FLog ${MPI_Fortran_LIBRARIES} PreSys zgrid)
TARGET_LINK_LIBRARIES(${EXE2} Utility PoissonSolver)
#TARGET_LINK_LIBRARIES(${EXE3} ${NETCDF_LIBRARY} State2NC Utility MPDD Geometry PoissonSolver)
TARGET_LINK_LIBRARIES(${EXE4} ${NETCDF_LIBRARY} State2NC Utility MPDD Geometry)
# TARGET_LINK_LIBRARIES(${EXE5} ${NETCDF_LIBRARY} State2NC Utility MPDD Geometry PoissonSolver)
# TARGET_LINK_LIBRARIES(${EXE6} ${NETCDF_LIBRARY} State2NC Utility MPDD Geometry PoissonSolver)
# TARGET_LINK_LIBRARIES(${EXE7} PoissonSolver MPDD Geometry Utility FLog ${MPI_Fortran_LIBRARIES} PreSys zgrid)

# 4. Add tests: (MIGRATED TO UnitTests/MOTOR-PS/CMakeLists.txt)
# ADD_TEST(NAME "Test_PossionSolver_Forward"
#         COMMAND mpirun -np 4 --oversubscribe Test_PSForward.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})
# ADD_TEST(NAME "Test_PossionSolver_Adjoint"
#         COMMAND mpirun -np 4 --oversubscribe Test_PSAdjoint.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})
# ADD_TEST(NAME "Test_LaplaceOpr"
#         COMMAND mpirun -np 4 --oversubscribe TestAnalytics.exe $ENV{STATIC_DIR}/Application/App_4DVAR.yaml
#         WORKING_DIRECTORY $ENV{BIN_DIR})
# ADD_TEST(NAME "TestPsiChi2UV"
#         COMMAND mpirun -np 4 --oversubscribe TestPsiChi2UV.exe $ENV{STATIC_DIR}/Application/App_4DVAR.yaml
#         WORKING_DIRECTORY $ENV{BIN_DIR})
        
# 5. Set test criteria: (MIGRATED TO UnitTests/MOTOR-PS/CMakeLists.txt)
# SET_TESTS_PROPERTIES(Test_PossionSolver_Forward PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
# SET_TESTS_PROPERTIES(Test_PossionSolver_Forward PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
# SET_TESTS_PROPERTIES(Test_PossionSolver_Adjoint PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
# SET_TESTS_PROPERTIES(Test_PossionSolver_Adjoint PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
# SET_TESTS_PROPERTIES(Test_LaplaceOpr PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
# SET_TESTS_PROPERTIES(Test_LaplaceOpr PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
# SET_TESTS_PROPERTIES(TestPsiChi2UV PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
# SET_TESTS_PROPERTIES(TestPsiChi2UV PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")

