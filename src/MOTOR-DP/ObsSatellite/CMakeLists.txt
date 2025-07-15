cmake_minimum_required(VERSION 3.5)

PROJECT(DpSatellite)


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

# Add rttov libraries.
INCLUDE_DIRECTORIES(${RTTOV_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${RTTOV_MOD_DIR})

# Add HDF5 libraries.
INCLUDE_DIRECTORIES(${HDF5_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(DP_Satellite STATIC ${DpSatellite_SRC})
TARGET_LINK_LIBRARIES(DP_Satellite ${NETCDF_LIBRARY} Geometry MPDD YAMLRead Field State ObsSet 
Mock State2NC Utility RTTOV ObsBase
${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

# 2. Build test executables:
# Test for Read Raw FY4 AGRI OBS
SET (EXE1 "Process_agri.exe")
SET (EXE2 "Process_giirs.exe")
SET (EXE3 "Process_hiras.exe")
SET (EXE4 "Process_mwhs.exe")
SET (EXE5 "Process_mwts.exe")

#IF ( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )

  # 3. Add codes to the executables:

  ADD_EXECUTABLE(${EXE1} Process_agri.F90)
  TARGET_LINK_LIBRARIES(${EXE1} IOGrapes C2O MGOpts Geometry MPDD YAMLRead
  Field State ObsSet DP_Satellite State2NC Utility ObsBase lapack RTTOV Obs2State ${MPI_Fortran_LIBRARIES}
  ${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
  ${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
  ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

  ADD_EXECUTABLE(${EXE2} Process_giirs.F90)
  TARGET_LINK_LIBRARIES(${EXE2} IOGrapes C2O MGOpts Geometry MPDD YAMLRead
  Field State ObsSet DP_Satellite State2NC Utility ObsBase lapack RTTOV Obs2State ${MPI_Fortran_LIBRARIES}
  ${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
  ${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
  ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

  ADD_EXECUTABLE(${EXE3} Process_hiras.F90)
  TARGET_LINK_LIBRARIES(${EXE3} IOGrapes C2O MGOpts Geometry MPDD YAMLRead
  Field State ObsSet DP_Satellite State2NC Utility ObsBase lapack RTTOV Obs2State ${MPI_Fortran_LIBRARIES}
  ${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
  ${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
  ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

  ADD_EXECUTABLE(${EXE4} Process_mwhs.F90)
  TARGET_LINK_LIBRARIES(${EXE4} IOGrapes C2O MGOpts Geometry MPDD YAMLRead
  Field State ObsSet DP_Satellite State2NC Utility ObsBase lapack RTTOV Obs2State ${MPI_Fortran_LIBRARIES}
  ${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
  ${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
  ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

  ADD_EXECUTABLE(${EXE5} Process_mwts.F90)
  TARGET_LINK_LIBRARIES(${EXE5} IOGrapes C2O MGOpts Geometry MPDD YAMLRead
  Field State ObsSet DP_Satellite State2NC Utility ObsBase lapack RTTOV Obs2State ${MPI_Fortran_LIBRARIES}
  ${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
  ${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
  ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})
  
  # # 4. Add tests:
  # ADD_TEST(NAME "Test_Raw_fy4_agri" 
  #         COMMAND mpirun -np 4 --oversubscribe Test_Raw_fy4_agri.exe
  #         WORKING_DIRECTORY $ENV{BIN_DIR})

  # ADD_TEST(NAME "FY4_AGRI_App" 
  #         COMMAND mpirun -np 4 --oversubscribe FY4_AGRI_App.exe
  #         WORKING_DIRECTORY $ENV{BIN_DIR})

  # ADD_TEST(NAME "Test_AGRI_thinning" 
  #         COMMAND mpirun -np 4 --oversubscribe Test_AGRI_thinning.exe
  #         WORKING_DIRECTORY $ENV{BIN_DIR})

  # ADD_TEST(NAME "Test_wrfinput_App" 
  #         COMMAND mpirun -np 4 --oversubscribe Test_wrfinput_App.exe
  #         WORKING_DIRECTORY $ENV{BIN_DIR})

#   # 5. Set test criteria:
#   SET_TESTS_PROPERTIES(Test_Raw_fy4_agri PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
#   SET_TESTS_PROPERTIES(Test_Raw_fy4_agri PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")

#   SET_TESTS_PROPERTIES(FY4_AGRI_App PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
#   SET_TESTS_PROPERTIES(FY4_AGRI_App PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")

#   SET_TESTS_PROPERTIES(Test_AGRI_thinning PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
#   SET_TESTS_PROPERTIES(Test_AGRI_thinning PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")

#   SET_TESTS_PROPERTIES(Test_wrfinput_App PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
#   SET_TESTS_PROPERTIES(Test_wrfinput_App PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")

# #ENDIF()
