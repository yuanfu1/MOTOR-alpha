cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(TimeIntegration)

# Include the common template:
IF (IS_ABSOLUTE $ENV{CTEMPLATES})
  MESSAGE(STATUS "CTEMPLATES path is set!")
ELSE ($ENV{CTEMPLATES} LESS 1)
  MESSAGE(FATAL_ERROR "Env CTEMPLATES is not set, please set and rerun!")
ENDIF (IS_ABSOLUTE $ENV{CTEMPLATES})
include($ENV{CTEMPLATES}/cmakeCommon.cmake)

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
ADD_LIBRARY(TimeIntegration STATIC ${TimeIntegration_SRC})
TARGET_LINK_LIBRARIES(TimeIntegration MPDD Geometry Utility M2CV ${MPI_Fortran_LIBRARIES} PoissonSolver)

# 2. Build test executables:
# SET (EXE1 "Test_rk4.exe")
# SET (EXE2 "Test_WarmBubble.exe")



# 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} Test_rk4.F90)
# TARGET_LINK_LIBRARIES(${EXE1} RK4 zgrid State State2NC MPDD Geometry Utility ${MPI_Fortran_LIBRARIES} Dyn PoissonSolver PreSys YAMLRead)

# ADD_EXECUTABLE(${EXE2} Test_warm_bubble.F90)
# TARGET_LINK_LIBRARIES(${EXE2} RK4 zgrid State State2NC MPDD Geometry Utility ${MPI_Fortran_LIBRARIES} Dyn PoissonSolver PreSys YAMLRead)

#  ADD_EXECUTABLE(${EXE2} Cumulus/Test_Cum.F90)
#  TARGET_LINK_LIBRARIES(${EXE2} Cumulus MPDD Geometry Utility State State2NC ${MPI_Fortran_LIBRARIES})

#  ADD_EXECUTABLE(${EXE3} Cumulus/Test_Cum_AD.F90)
#  TARGET_LINK_LIBRARIES(${EXE3} Cumulus MPDD Geometry Utility State State2NC ${MPI_Fortran_LIBRARIES} TapenadeAux YAMLRead IOGrapes UV2W MGOpts)

#  ADD_EXECUTABLE(${EXE4} Cumulus/Test_Cum_interface.F90)
#  TARGET_LINK_LIBRARIES(${EXE4} Cumulus MPDD Geometry Utility State State2NC ${MPI_Fortran_LIBRARIES} TapenadeAux YAMLRead IOGrapes UV2W MGOpts)

# #  ADD_EXECUTABLE(${EXE5} Cumulus/Test_CumJfun.F90)
# #  TARGET_LINK_LIBRARIES(${EXE5} Minimization Cumulus MPDD Geometry Utility State State2NC ${MPI_Fortran_LIBRARIES} 
# #  TapenadeAux BMatrix C2O JFunc MGOpts YAMLRead Field State ObsSet Mock Obs2State OprRadar RTTOV ${RTTOV_LIBRARY} )

#  ADD_EXECUTABLE(${EXE6} Cumulus/Test_cum_GRAPES.F90)
#  TARGET_LINK_LIBRARIES(${EXE6} Cumulus MPDD Geometry Utility State State2NC ${MPI_Fortran_LIBRARIES} TapenadeAux YAMLRead IOGrapes UV2W)

# #  ADD_EXECUTABLE(${EXE7} shared/pushpop_test.F90)
# #  TARGET_LINK_LIBRARIES(${EXE7} Cumulus ${MPI_Fortran_LIBRARIES} TapenadeAux)

# #  ADD_EXECUTABLE(${EXE8} shared/Test_interp1d_AD.F90)
# #  TARGET_LINK_LIBRARIES(${EXE8} Cumulus MPDD Geometry Utility State State2NC ${MPI_Fortran_LIBRARIES} TapenadeAux YAMLRead IOGrapes)
 



#  # # 4. Add tests:
#  ADD_TEST(NAME Test_Cumulus
#          COMMAND Test_Cum_AD.exe Test_Cum_interface.exe
#          WORKING_DIRECTORY $ENV{BIN_DIR})

# # # 5. Set test criteria:
#  SET_TESTS_PROPERTIES(Test_Cumulus PROPERTIES FAIL_REGULAR_EXPRESSION "tested failed")
#  SET_TESTS_PROPERTIES(Test_Cumulus PROPERTIES PASS_REGULAR_EXPRESSION "tested passed")


