cmake_minimum_required(VERSION 3.5)

PROJECT(ObsGNSS)


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
#INCLUDE_DIRECTORIES(${RTTOV_INCLUDE_DIR})
#INCLUDE_DIRECTORIES(${RTTOV_MOD_DIR})

# Add HDF5 libraries.
INCLUDE_DIRECTORIES(${HDF5_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(ObsGNSS STATIC ${ObsGNSS_SRC})
TARGET_LINK_LIBRARIES(ObsGNSS ${NETCDF_LIBYRAR} Geometry MPDD YAMLRead Field State ObsSet 
Mock State2NC Utility ObsRadar ObsBase
${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

# 2. Build test executables:
#SET (EXE2 "Test_ObsGNSS.exe")
#ADD_EXECUTABLE(${EXE2} Test_ObsGNSS.F90)
# TARGET_LINK_LIBRARIES(${EXE2} C2O MGOpts Geometry MPDD YAMLRead
# Field State ObsSet ObsSatellite Mock State2NC Utility ObsRadar ObsBase ObsGNSS ObsUtilities lapack Obs2State ${MPI_Fortran_LIBRARIES}
# ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

#IF ( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )

  # 3. Add codes to the executables:
  # ADD_EXECUTABLE(${EXE1} Test_Raw_fy4_agri.F90)
  # TARGET_LINK_LIBRARIES(${EXE1} C2O MGOpts Geometry MPDD YAMLRead
  # Field State ObsSet ObsSatellite Mock State2NC Utility ObsRadar ObsBase ObsUtilities lapack RTTOV Obs2State ${MPI_Fortran_LIBRARIES}
  # ${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
  # ${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
  # ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

  # ADD_EXECUTABLE(${EXE2} FY4_AGRI_App.F90)
  # TARGET_LINK_LIBRARIES(${EXE2} Minimization BMatrix IOGrapes C2O JFunc IOWRF MGOpts Geometry MPDD YAMLRead
  # Field State ObsSet ObsSatellite Mock State2NC Utility ObsRadar ObsBase ObsUtilities lapack RTTOV Obs2State ${MPI_Fortran_LIBRARIES}
  # ${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
  # ${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
  # ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

  # ADD_EXECUTABLE(${EXE3} Test_AGRI_thinning.F90)
  # TARGET_LINK_LIBRARIES(${EXE3} Minimization BMatrix IOGrapes C2O JFunc IOWRF MGOpts Geometry MPDD YAMLRead
  # Field State ObsSet ObsSatellite Mock State2NC Utility ObsRadar ObsBase ObsUtilities lapack RTTOV Obs2State ${MPI_Fortran_LIBRARIES}
  # ${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
  # ${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
  # ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

  # ADD_EXECUTABLE(${EXE4} Test_Raw_fy4_giirs.F90)
  # TARGET_LINK_LIBRARIES(${EXE4} C2O MGOpts Geometry MPDD YAMLRead
  # Field State ObsSet ObsSatellite Mock State2NC Utility ObsRadar ObsBase ObsUtilities lapack RTTOV Obs2State ${MPI_Fortran_LIBRARIES}
  # ${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
  # ${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
  # ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

  # ADD_EXECUTABLE(${EXE5} FY4_GIIRS_App.F90)
  # TARGET_LINK_LIBRARIES(${EXE5} Minimization BMatrix IOGrapes C2O JFunc IOWRF MGOpts Geometry MPDD YAMLRead
  # Field State ObsSet ObsSatellite Mock State2NC Utility ObsRadar ObsBase ObsUtilities lapack RTTOV Obs2State ${MPI_Fortran_LIBRARIES}
  # ${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
  # ${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
  # ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

  # ADD_EXECUTABLE(${EXE6} Test_GIIRS_thinning.F90)
  # TARGET_LINK_LIBRARIES(${EXE6} Minimization BMatrix IOGrapes C2O JFunc IOWRF MGOpts Geometry MPDD YAMLRead
  # Field State ObsSet ObsSatellite Mock State2NC Utility ObsRadar ObsBase ObsUtilities lapack RTTOV Obs2State ${MPI_Fortran_LIBRARIES}
  # ${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
  # ${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
  # ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

#   ADD_EXECUTABLE(${EXE4} Test_wrfinput_App.F90)
#   TARGET_LINK_LIBRARIES(${EXE4} Minimization BMatrix IOGrapes C2O JFunc IOWRF MGOpts Geometry MPDD YAMLRead
#   Field State ObsSet ObsSatellite Mock State2NC Utility ObsRadar ObsBase ObsUtilities RTTOV Obs2State ${MPI_Fortran_LIBRARIES}
#   ${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
#   ${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
#   ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})
  
#   # 4. Add tests:
#   ADD_TEST(NAME "Test_Raw_fy4_agri" 
#           COMMAND mpirun -np 4 --oversubscribe Test_Raw_fy4_agri.exe
#           WORKING_DIRECTORY $ENV{BIN_DIR})

#   ADD_TEST(NAME "FY4_AGRI_App" 
#           COMMAND mpirun -np 4 --oversubscribe FY4_AGRI_App.exe
#           WORKING_DIRECTORY $ENV{BIN_DIR})

#   ADD_TEST(NAME "Test_AGRI_thinning" 
#           COMMAND mpirun -np 4 --oversubscribe Test_AGRI_thinning.exe
#           WORKING_DIRECTORY $ENV{BIN_DIR})

#   ADD_TEST(NAME "Test_wrfinput_App" 
#           COMMAND mpirun -np 4 --oversubscribe Test_wrfinput_App.exe
#           WORKING_DIRECTORY $ENV{BIN_DIR})

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
