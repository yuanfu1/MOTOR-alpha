cmake_minimum_required(VERSION 3.12)

# Build MOTOR-PS Unit Tests
PROJECT(MOTOR-PS_UnitTests)

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
# Build test executables for MOTOR-PS
#-------------------------------------------------------------------------------


# MOTOR-PS Test Executables

# Poisson Solver Tests
ADD_EXECUTABLE(Test_PSAdjoint.exe Test_PSAdjoint.F90)
TARGET_LINK_LIBRARIES(Test_PSAdjoint.exe PoissonSolver Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_PS_chen.exe Test_PS_chen.F90)
TARGET_LINK_LIBRARIES(Test_PS_chen.exe PoissonSolver Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_PSForward.exe Test_PSForward.F90)
TARGET_LINK_LIBRARIES(Test_PSForward.exe PoissonSolver Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(TestPsiChi2UV.exe TestPsiChi2UV.F90)
TARGET_LINK_LIBRARIES(TestPsiChi2UV.exe PoissonSolver Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(TestAnalytics.exe TestAnalytics.F90)
TARGET_LINK_LIBRARIES(TestAnalytics.exe PoissonSolver Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# PreSys Tests
# ADD_EXECUTABLE(Test_CalVerDer.exe Test_CalVerDer.F90)
# TARGET_LINK_LIBRARIES(Test_CalVerDer.exe PreSys Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_bdyidx.exe Test_bdyidx.F90)
TARGET_LINK_LIBRARIES(Test_bdyidx.exe PreSys Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_InterpValue.exe Test_InterpValue.F90)
TARGET_LINK_LIBRARIES(Test_InterpValue.exe PreSys Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# Dynamics Tests
ADD_EXECUTABLE(Test_tenW.exe Test_tenW.F90)
TARGET_LINK_LIBRARIES(Test_tenW.exe Dyn Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_tendiv.exe Test_tendiv.F90)
TARGET_LINK_LIBRARIES(Test_tendiv.exe Dyn Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_tenRho.exe Test_tenRho.F90)
TARGET_LINK_LIBRARIES(Test_tenRho.exe Dyn Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(Test_tenvor.exe Test_tenvor.F90)
TARGET_LINK_LIBRARIES(Test_tenvor.exe Dyn Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ZGrid Tests
ADD_EXECUTABLE(Test_gzm.exe Test_gzm.F90)
TARGET_LINK_LIBRARIES(Test_gzm.exe ZGrid Dyn Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

ADD_EXECUTABLE(test_laplacian_adjoint test_laplacian_adjoint.F90)
TARGET_LINK_LIBRARIES(test_laplacian_adjoint
	PRIVATE
	gzm_lib
	State
	Geometry
	Field
	MPDD
	Utility
	YAMLRead
	PreSys
	State2NC
	${MPI_Fortran_LIBRARIES}
	${NETCDF_LIBRARIES}
)

# ADD_EXECUTABLE(Test_gzm_use_Dyn_Function.exe Test_gzm_use_Dyn_Function.F90)
# TARGET_LINK_LIBRARIES(Test_gzm_use_Dyn_Function.exe ZGrid Dyn Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_gzm_Xie.exe Test_gzm_Xie.F90)
# TARGET_LINK_LIBRARIES(Test_gzm_Xie.exe ZGrid Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_gzm_Chen.exe Test_gzm_Chen.F90)
# TARGET_LINK_LIBRARIES(Test_gzm_Chen.exe ZGrid Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# GZM Update Tests
# ADD_EXECUTABLE(Test_GZM_TL_AD.exe Test_GZM_TL_AD.F90)
# TARGET_LINK_LIBRARIES(Test_GZM_TL_AD.exe gzm_upd Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_PreCal_TL_AD.exe Test_PreCal_TL_AD.F90)
# TARGET_LINK_LIBRARIES(Test_PreCal_TL_AD.exe gzm_upd Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_CalVerDer_tlm.exe Test_CalVerDer_tlm.F90)
# TARGET_LINK_LIBRARIES(Test_CalVerDer_tlm.exe gzm_upd Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_CalVerDer_TL_AD.exe Test_CalVerDer_TL_AD.F90)
# TARGET_LINK_LIBRARIES(Test_CalVerDer_TL_AD.exe gzm_upd Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_CalVerDer_AD_without_TL.exe Test_CalVerDer_AD_without_TL.F90)
# TARGET_LINK_LIBRARIES(Test_CalVerDer_AD_without_TL.exe gzm_upd Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# Time Integration Tests
# ADD_EXECUTABLE(Test_warm_bubble.exe Test_warm_bubble.F90)
# TARGET_LINK_LIBRARIES(Test_warm_bubble.exe TimeIntegration Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

# ADD_EXECUTABLE(Test_rk4.exe Test_rk4.F90)
# TARGET_LINK_LIBRARIES(Test_rk4.exe TimeIntegration Geometry MPDD YAMLRead Utility ${MPI_Fortran_LIBRARIES})

#-------------------------------------------------------------------------------
# CTest configurations for MOTOR-PS
#-------------------------------------------------------------------------------

# Add tests from original PoissonSolver_Kp
ADD_TEST(NAME "Test_PossionSolver_Forward"
        COMMAND mpirun -np 4 --oversubscribe Test_PSForward.exe
        WORKING_DIRECTORY $ENV{BIN_DIR})
ADD_TEST(NAME "Test_PossionSolver_Adjoint"
        COMMAND mpirun -np 4 --oversubscribe Test_PSAdjoint.exe
        WORKING_DIRECTORY $ENV{BIN_DIR})
ADD_TEST(NAME "Test_LaplaceOpr"
        COMMAND mpirun -np 4 --oversubscribe TestAnalytics.exe $ENV{STATIC_DIR}/Application/App_4DVAR.yaml
        WORKING_DIRECTORY $ENV{BIN_DIR})
ADD_TEST(NAME "TestPsiChi2UV"
        COMMAND mpirun -np 4 --oversubscribe TestPsiChi2UV.exe $ENV{STATIC_DIR}/Application/App_4DVAR.yaml
        WORKING_DIRECTORY $ENV{BIN_DIR})

# Add test from original ZGrid
ADD_TEST(NAME LaplacianAdjoint 
         COMMAND test_laplacian_adjoint
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Set test criteria for PoissonSolver tests
SET_TESTS_PROPERTIES(Test_PossionSolver_Forward PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
SET_TESTS_PROPERTIES(Test_PossionSolver_Forward PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
SET_TESTS_PROPERTIES(Test_PossionSolver_Adjoint PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
SET_TESTS_PROPERTIES(Test_PossionSolver_Adjoint PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
SET_TESTS_PROPERTIES(Test_LaplaceOpr PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
SET_TESTS_PROPERTIES(Test_LaplaceOpr PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
SET_TESTS_PROPERTIES(TestPsiChi2UV PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
SET_TESTS_PROPERTIES(TestPsiChi2UV PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")

