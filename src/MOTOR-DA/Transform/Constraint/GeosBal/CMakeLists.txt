cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(GeosBal)

# Include the common template:
IF (IS_ABSOLUTE $ENV{CTEMPLATES})
  MESSAGE(STATUS "CTEMPLATES path is set!")
ELSE ($ENV{CTEMPLATES} LESS 1)
  MESSAGE(FATAL_ERROR "Env CTEMPLATES is not set, please set and rerun!")
ENDIF (IS_ABSOLUTE $ENV{CTEMPLATES})
include($ENV{CTEMPLATES}/cmakeCommon.cmake)

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
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(GeosBal STATIC ${GeosBal_SRC})
TARGET_LINK_LIBRARIES(GeosBal PUBLIC State State2NC Field Geometry Utility TransBase lapack blas)

# # 2. Build test executables:
# #wangjia SET (EXE1 "Test_QRSolver.exe")
# SET (EXE1 "Test_NewGeosBal.exe")
# SET (EXE2 "Test_GeosBalTLAD.exe")

# # 3. Add codes to the executables:
# #wangjia ADD_EXECUTABLE(${EXE1} Test_QRSolver.F90)
# #wangjia TARGET_LINK_LIBRARIES(${EXE2} GeosBal lapack blas Minimization Geometry Utility JFunc State MPDD)

# ADD_EXECUTABLE(${EXE1} Test_NewGeosBal.F90)
# ADD_EXECUTABLE(${EXE2} Test_GeosBalTLAD.F90)
# TARGET_LINK_LIBRARIES(${EXE1} GeosBal Geometry Application Minimization Utility Field State State2NC BMatrix RMatrix ObsConvention ObsRadar RTTOV ${RTTOV_LIBRARY} ObsSatellite ObsBase ObsSet ObsUtilities Mock MPDD State2NC C2O MGOpts ${MPI_Fortran_LIBRARIES} IOGrapes PostProc)
# TARGET_LINK_LIBRARIES(${EXE2} GeosBal Geometry Utility State MPDD State2NC ${MPI_Fortran_LIBRARIES} IOGrapes)

# # 4. Add tests:
# ADD_TEST(NAME "Test_GeosBal"
#         COMMAND mpirun -np 4 --oversubscribe Test_NewGeosBal.exe $ENV{STATIC_DIR}/UnitTest/Test_NewGeosBal.yaml
#         WORKING_DIRECTORY $ENV{BIN_DIR})

# # 5. Set test criteria:
# SET_TESTS_PROPERTIES(Test_GeosBal PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
# SET_TESTS_PROPERTIES(Test_GeosBal PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
