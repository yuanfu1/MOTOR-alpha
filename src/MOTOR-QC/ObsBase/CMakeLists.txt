cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(ObsBase)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
# INCLUDE_DIRECTORIES(/Users/yjzx-xieyuanfu/developments/mods)
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})
# LINK_DIRECTORIES(${NETCDF_LIBRARY})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

# Add Slint libraries.
INCLUDE_DIRECTORIES(${SLINT_Fortran_INCLUDE_DIR})
# LINK_DIRECTORIES(${SLINT_LIBRARY})

#-------------------------------------------------------------------------------
# Source codes to this build:
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(ObsBase STATIC ${ObsBase_SRC})
TARGET_LINK_LIBRARIES(ObsBase ${NETCDF_LIBRARY} FLog Utility M2ODirect State ObsSet ${SLINT_LIBRARY})
# TARGET_LINK_LIBRARIES(ObsBase ${NETCDF_LIBRARY} FLog Utility M2ODirect State ObsSet ${SLINT_LIBRARY} RTTOV )

# # 2. Build test executables:
# SET (EXE1 "Test_griddedObsMap.exe")

# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} Test_griddedObsMap.F90)
# TARGET_LINK_LIBRARIES(${EXE1} Geometry MPDD YAMLRead Field State ObsSet ObsBase ObsConvention C2O Mock State2NC Obs2State MGOpts)

# # 4. Add tests:
# ADD_TEST(NAME "Test_griddedObsMap" COMMAND mpirun -np 4 --oversubscribe Test_griddedObsMap.exe
# 	WORKING_DIRECTORY $ENV{BIN_DIR})

# 5. Set test criteria:
# SET_TESTS_PROPERTIES(Test_griddedObsMap PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
# SET_TESTS_PROPERTIES(Test_griddedObsMap PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
