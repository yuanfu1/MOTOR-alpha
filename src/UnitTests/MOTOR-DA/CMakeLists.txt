cmake_minimum_required(VERSION 3.12)

# Build MOTOR-DA Unit Tests
PROJECT(MOTOR-DA_UnitTests)

#-------------------------------------------------------------------------------
#  Include paths and library directories
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

# Add Netcdf includes
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})

# Add rttov libraries.
INCLUDE_DIRECTORIES(${RTTOV_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${RTTOV_MOD_DIR})

# Add HDF5 libraries.
IF (CMAKE_Fortran_COMPILER MATCHES "ifort")
	INCLUDE_DIRECTORIES(${HDF5_Fortran_INCLUDE_DIR})
ENDIF ()

#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Build test executables for MOTOR-DA
#-------------------------------------------------------------------------------


# MOTOR-DA Test Executables

# BMatrix Tests
ADD_EXECUTABLE(Test_BMatrix.exe Test_BMatrix.F90)
TARGET_LINK_LIBRARIES(Test_BMatrix.exe BMatrix Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES} linpack lapack blas)

ADD_EXECUTABLE(Test_BMatTLAD.exe Test_BMatTLAD.F90)
TARGET_LINK_LIBRARIES(Test_BMatTLAD.exe IOGrapes BMatrix Geometry MPDD YAMLRead Utility State Field EnLoc ${MPI_Fortran_LIBRARIES})

# BField Tests
# ADD_EXECUTABLE(Test_BField_Laplace_Terrain.exe Test_BField_Laplace_Terrain.F90)
# TARGET_LINK_LIBRARIES(Test_BField_Laplace_Terrain.exe BField Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_EnLoc_ADcheck.exe Test_EnLoc_ADcheck.F90)
TARGET_LINK_LIBRARIES(Test_EnLoc_ADcheck.exe EnLoc BField Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# EnLoc Tests
ADD_EXECUTABLE(Test_QRrecover.exe Test_QRrecover.F90)
TARGET_LINK_LIBRARIES(Test_QRrecover.exe EnLoc Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_EnLoc_state.exe Test_EnLoc_state.F90)
# TARGET_LINK_LIBRARIES(Test_EnLoc_state.exe EnLoc Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_QR_Multigrid.exe Test_QR_Multigrid.F90)
# TARGET_LINK_LIBRARIES(Test_QR_Multigrid.exe EnLoc Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_QRdecomp_invR.exe Test_QRdecomp_invR.F90)
TARGET_LINK_LIBRARIES(Test_QRdecomp_invR.exe EnLoc Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_BKErr_model_FWD_ADJ_Binv.exe Test_BKErr_model_FWD_ADJ_Binv.F90)
# TARGET_LINK_LIBRARIES(Test_BKErr_model_FWD_ADJ_Binv.exe EnLoc Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_QRdecomp_invRv2.exe Test_QRdecomp_invRv2.F90)
TARGET_LINK_LIBRARIES(Test_QRdecomp_invRv2.exe EnLoc Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_EnLoc_Multigrid.exe Test_EnLoc_Multigrid.F90)
# TARGET_LINK_LIBRARIES(Test_EnLoc_Multigrid.exe EnLoc Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_EnLoc_model.exe Test_EnLoc_model.F90)
# TARGET_LINK_LIBRARIES(Test_EnLoc_model.exe EnLoc Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# Cumulus Tests
# ADD_EXECUTABLE(Test_cum_GRAPES.exe Test_cum_GRAPES.F90)
# TARGET_LINK_LIBRARIES(Test_cum_GRAPES.exe Cumulus IOGrapes Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_CumJfun.exe Test_CumJfun.F90)
# TARGET_LINK_LIBRARIES(Test_CumJfun.exe Cumulus Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_Cum_AD.exe Test_Cum_AD.F90)
TARGET_LINK_LIBRARIES(Test_Cum_AD.exe Cumulus Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_Cum_interface.exe Test_Cum_interface.F90)
# TARGET_LINK_LIBRARIES(Test_Cum_interface.exe Cumulus Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_Cum.exe Test_Cum.F90)
# TARGET_LINK_LIBRARIES(Test_Cum.exe Cumulus Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# IO Tests
ADD_EXECUTABLE(Test_State2NC.exe Test_State2NC.F90)
TARGET_LINK_LIBRARIES(Test_State2NC.exe State2NC IOGrapes Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES} ${NETCDF_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_IOWRF_Ens.exe Test_IOWRF_Ens.F90)
TARGET_LINK_LIBRARIES(Test_IOWRF_Ens.exe IOWRF Geometry MPDD YAMLRead Utility State State2NC Field ${MPI_Fortran_LIBRARIES} ${NETCDF_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_IOWRF.exe Test_IOWRF.F90)
TARGET_LINK_LIBRARIES(Test_IOWRF.exe IOWRF Geometry MPDD YAMLRead Utility State State2NC Field ${MPI_Fortran_LIBRARIES} ${NETCDF_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_IOGrapesPostvar.exe Test_IOGrapesPostvar.F90)
TARGET_LINK_LIBRARIES(Test_IOGrapesPostvar.exe IOGrapes Geometry MPDD YAMLRead Utility State State2NC Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_IOGrapesModelvar.exe Test_IOGrapesModelvar.F90)
TARGET_LINK_LIBRARIES(Test_IOGrapesModelvar.exe IOGrapes Geometry MPDD YAMLRead Utility State State2NC Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_IOGrapesModelvar_Ens.exe Test_IOGrapesModelvar_Ens.F90)
TARGET_LINK_LIBRARIES(Test_IOGrapesModelvar_Ens.exe IOGrapes Geometry MPDD YAMLRead Utility State State2NC Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_IOGrapesModelvar_500m_to_nc.exe Test_IOGrapesModelvar_500m_to_nc.F90)
TARGET_LINK_LIBRARIES(Test_IOGrapesModelvar_500m_to_nc.exe PostProc IOGrapes Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES} ${NETCDF_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_IOGrapes.exe Test_IOGrapes.F90)
TARGET_LINK_LIBRARIES(Test_IOGrapes.exe IOGrapes Geometry MPDD YAMLRead Utility State State2NC Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_IOERA5.exe Test_IOERA5.F90)
TARGET_LINK_LIBRARIES(Test_IOERA5.exe IOERA5 IOGrapes Geometry MPDD YAMLRead Utility State State2NC Field ${MPI_Fortran_LIBRARIES} ${NETCDF_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_IOECM.exe Test_IOECM.F90)
TARGET_LINK_LIBRARIES(Test_IOECM.exe IOECM Geometry MPDD YAMLRead Utility State State2NC Field ${MPI_Fortran_LIBRARIES} ${NETCDF_Fortran_LIBRARIES})

# Other DA Tests
ADD_EXECUTABLE(Test_Export2SelDomain.exe Test_Export2SelDomain.F90)
TARGET_LINK_LIBRARIES(Test_Export2SelDomain.exe PostProc Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_ObsSet.exe Test_ObsSet.F90)
TARGET_LINK_LIBRARIES(Test_ObsSet.exe ObsSet Geometry MPDD YAMLRead Utility State Field Mock ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_C2O.exe Test_C2O.F90)
TARGET_LINK_LIBRARIES(Test_C2O.exe C2O Geometry MPDD YAMLRead Utility State Field ObsSet Mock ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_State.exe Test_State.F90)
TARGET_LINK_LIBRARIES(Test_State.exe State Geometry MPDD YAMLRead Utility Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_SolverFRCG.exe Test_SolverFRCG.F90)
TARGET_LINK_LIBRARIES(Test_SolverFRCG.exe Minimization JFunc Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_MiniSolver.exe Test_MiniSolver.F90)
TARGET_LINK_LIBRARIES(Test_MiniSolver.exe Minimization JFunc Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# Transform Tests
ADD_EXECUTABLE(Test_OprGNSSRefrac.exe Test_OprGNSSRefrac.F90)
TARGET_LINK_LIBRARIES(Test_OprGNSSRefrac.exe OprGNSS ObsBase ObsGNSS Geometry MPDD YAMLRead Utility State Field ObsSet C2O IOGrapes ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_UV2DirSpd.exe Test_UV2DirSpd.F90)
TARGET_LINK_LIBRARIES(Test_UV2DirSpd.exe UV2DirSpd Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_RadarVel.exe Test_RadarVel.F90)
TARGET_LINK_LIBRARIES(Test_RadarVel.exe OprRadar Geometry MPDD YAMLRead Utility State Field ObsSet ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_UV2WForward_SurfInte.exe Test_UV2WForward_SurfInte.F90)
TARGET_LINK_LIBRARIES(Test_UV2WForward_SurfInte.exe IOGrapes State2NC UV2W Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_UVW2DivgTLAD.exe Test_UVW2DivgTLAD.F90)
TARGET_LINK_LIBRARIES(Test_UVW2DivgTLAD.exe UV2W IOGrapes State2NC Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_QRSolver.exe Test_QRSolver.F90)
TARGET_LINK_LIBRARIES(Test_QRSolver.exe UV2W Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_UV2WForward_Poisson.exe Test_UV2WForward_Poisson.F90)
TARGET_LINK_LIBRARIES(Test_UV2WForward_Poisson.exe IOGrapes UV2W PoissonSolver Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_UV2WTLAD.exe Test_UV2WTLAD.F90)
TARGET_LINK_LIBRARIES(Test_UV2WTLAD.exe UV2W IOGrapes State2NC Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_UV2WParaTest.exe Test_UV2WParaTest.F90)
# TARGET_LINK_LIBRARIES(Test_UV2WParaTest.exe UV2W Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_NewGeosBal.exe Test_NewGeosBal.F90)
TARGET_LINK_LIBRARIES(Test_NewGeosBal.exe Application IOECM GeosBal Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_GeosBalTLAD.exe Test_GeosBalTLAD.F90)
TARGET_LINK_LIBRARIES(Test_GeosBalTLAD.exe GeosBal Geometry MPDD YAMLRead Utility State Field ${MPI_Fortran_LIBRARIES})

#-------------------------------------------------------------------------------
# CTest configurations for MOTOR-DA
#-------------------------------------------------------------------------------

# Add tests from original ObsSet
ADD_TEST(NAME Test_ObsSet 
        COMMAND mpirun -np 4 --oversubscribe  Test_ObsSet.exe
        WORKING_DIRECTORY $ENV{BIN_DIR})

# Add tests from original States
ADD_TEST(NAME "Test_States" 
        COMMAND mpirun -np 4 --oversubscribe Test_State.exe
        WORKING_DIRECTORY $ENV{BIN_DIR})

# Add tests from original MiniSolver
ADD_TEST(NAME Test_Solver_FRCG 
        COMMAND mpirun -np 4 --oversubscribe Test_SolverFRCG.exe
        WORKING_DIRECTORY $ENV{BIN_DIR})

ADD_TEST(NAME Test_Solver_LBFGS 
        COMMAND mpirun -np 1 --oversubscribe Test_MiniSolver.exe
        WORKING_DIRECTORY $ENV{BIN_DIR})

# Set test criteria
SET_TESTS_PROPERTIES(Test_ObsSet PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
SET_TESTS_PROPERTIES(Test_ObsSet PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")

SET_TESTS_PROPERTIES(Test_States PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
SET_TESTS_PROPERTIES(Test_States PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")

SET_TESTS_PROPERTIES(Test_Solver_FRCG PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
SET_TESTS_PROPERTIES(Test_Solver_FRCG PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")

SET_TESTS_PROPERTIES(Test_Solver_LBFGS PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
SET_TESTS_PROPERTIES(Test_Solver_LBFGS PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")

