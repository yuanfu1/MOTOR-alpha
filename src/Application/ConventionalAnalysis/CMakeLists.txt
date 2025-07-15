cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(ConventionalAnalysis)


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
#include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
#ADD_LIBRARY(ConventionalAnalysis STATIC ${ConventionalAnalysis_SRC})
#TARGET_LINK_LIBRARIES(ConventionalAnalysis ${NETCDF_LIBRARY} MGOpts Utility M2ODirect State ObsSet)

# 1. Build test executables:
SET (EXE1 "plotAnalytic.exe")
#SET (EXE2 "UnitTest_sfc_randoms.exe")
#SET (EXE3 "UnitTest_Convention.exe")
#SET (EXE4 "UnitTest_Convention_Adv.exe")
#SET (EXE5 "UnitTest_Convention_MultiObs.exe")
#SET (EXE6 "UnitTest_Conventions.exe")
#SET (EXE7 "UnitTest_Convention_MultiObsVcycle.exe")

# 2. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} plotAnalytic.F90)
#ADD_EXECUTABLE(${EXE2} UnitTest_sfc_randoms.F90)
# ADD_EXECUTABLE(${EXE3} UnitTest_Convention.F90)
# ADD_EXECUTABLE(${EXE4} UnitTest_Convention_Adv.F90)
# ADD_EXECUTABLE(${EXE5} UnitTest_Convention_MultiObs.F90)
# ADD_EXECUTABLE(${EXE6} UnitTest_Conventions.F90)
# ADD_EXECUTABLE(${EXE7} UnitTest_Convention_MultiObsVcycle.F90)


# TARGET_LINK_LIBRARIES(${EXE1} Application Utility Field State State2NC  Obs2State)
# TARGET_LINK_LIBRARIES(${EXE2} MGOpts MPObs MPDD Geometry Application Utility Minimization Field State ObsSet ObsBase ObsConvention BMatrix RMatrix C2O Mock State2NC JFunc MGOpts Obs2State RTTOV ${RTTOV_LIBRARY} ${MPI_Fortran_LIBRARIES})
# TARGET_LINK_LIBRARIES(${EXE3} MGOpts MPObs MPDD Geometry Application Utility Minimization Field State ObsSet ObsBase ObsConvention BMatrix RMatrix C2O Mock State2NC JFunc MGOpts Obs2State RTTOV ${RTTOV_LIBRARY} ${MPI_Fortran_LIBRARIES})
# TARGET_LINK_LIBRARIES(${EXE4} MGOpts MPObs MPDD Geometry Application Utility Minimization Field State ObsSet ObsBase ObsConvention BMatrix RMatrix C2O Mock State2NC JFunc MGOpts Obs2State RTTOV ${RTTOV_LIBRARY} ${MPI_Fortran_LIBRARIES})
# TARGET_LINK_LIBRARIES(${EXE5} MGOpts MPObs MPDD Geometry Application Utility Minimization Field State ObsSet ObsBase ObsConvention ObsUtilities BMatrix RMatrix C2O Mock State2NC JFunc MGOpts Obs2State RTTOV ${RTTOV_LIBRARY} ${MPI_Fortran_LIBRARIES})
# TARGET_LINK_LIBRARIES(${EXE6} MGOpts MPObs MPDD Geometry Application Utility Minimization Field State ObsSet ObsBase ObsConvention ObsUtilities BMatrix ObsRadar ObsSatellite RMatrix C2O Mock State2NC JFunc MGOpts Obs2State Application RTTOV ${RTTOV_LIBRARY} ${MPI_Fortran_LIBRARIES} PostProc)
# TARGET_LINK_LIBRARIES(${EXE7} MGOpts MPObs MPDD Geometry Application Utility Minimization Field State ObsSet ObsBase ObsConvention ObsUtilities BMatrix ObsRadar RMatrix C2O Mock State2NC JFunc MGOpts Obs2State Application RTTOV ${RTTOV_LIBRARY} ${MPI_Fortran_LIBRARIES})

# 3. Add tests: (MIGRATED TO UnitTests/Application/CMakeLists.txt)
# ADD_TEST(NAME "Test_sfc_random" COMMAND mpirun -np 4 --oversubscribe UnitTest_sfc_randoms.exe
#	WORKING_DIRECTORY $ENV{BIN_DIR})
#ADD_TEST(NAME "Test_Convention_MultiObs" COMMAND mpirun -np 4 --oversubscribe UnitTest_Convention_MultiObs.exe
#	WORKING_DIRECTORY $ENV{BIN_DIR})

# 4. Set test criteria: (MIGRATED TO UnitTests/Application/CMakeLists.txt)
# SET_TESTS_PROPERTIES(Test_sfc_random PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
# SET_TESTS_PROPERTIES(Test_sfc_random PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
#SET_TESTS_PROPERTIES(Test_Convention_MultiObs PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
#SET_TESTS_PROPERTIES(Test_Convention_MultiObs PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
