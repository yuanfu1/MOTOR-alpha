cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(ObsSatob)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

# Add Slint libraries.
INCLUDE_DIRECTORIES(${SLINT_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
# Source codes to this build:
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(ObsSatob STATIC ${ObsSatob_SRC})
TARGET_LINK_LIBRARIES(ObsSatob ${NETCDF_LIBRARY} FLog Utility ObsBase NMLRead MPDD M2CV Obs2State)

# 2. Build test executables:
#SET (EXE1 "Test_ObsSatobPrep.exe")
#SET (EXE2 "Test_SatobSuper.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_ObsSatobPrep.F90)
# TARGET_LINK_LIBRARIES(${EXE1} IOGrapes C2O Geometry MPDD NMLRead Field State ObsSet
# Mock State2NC Utility ObsRadar ObsSatob ObsBase ObsUtilities RTTOV Obs2State ${MPI_Fortran_LIBRARIES})
#TARGET_LINK_LIBRARIES(${EXE1} ObsSatob Field State ObsSet Utility)

## 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE2} Test_ObsSatobSuper.F90)
#TARGET_LINK_LIBRARIES(${EXE2} Geometry MPDD Utility C2O MPObs Field State ObsSet State2NC RMatrix ObsSatob ${MPI_Fortran_LIBRARIES})



# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE2} Test_ObsSatobSuper.F90)
#TARGET_LINK_LIBRARIES(${EXE2} Geometry MPDD IOGrapes Utility Minimization NMLRead
#JFunc RMatrix BMatrix C2O MPObs State ObsSet State2NC meteoro_constants extreme_value Satob gtsqc interhp RTTOV ${MPI_Fortran_LIBRARIES} 
#${RTTOV_LIBRARY} )

# # 4. Add tests:
# ADD_TEST(NAME "Test_VelUVW2RwndTLAD" COMMAND mpirun -np 4 Test_VelUVW2RwndTLAD.exe
# 	WORKING_DIRECTORY $ENV{BIN_DIR})meteoro_constants extreme_value 

# # 5. Set test criteria:
# SET_TESTS_PROPERTIES(Test_VelUVW2RwndTLAD PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
# SET_TESTS_PROPERTIES(Test_VelUVW2RwndTLAD PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
