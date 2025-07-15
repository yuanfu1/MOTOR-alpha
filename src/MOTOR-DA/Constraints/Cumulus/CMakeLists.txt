cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(Cumulus)

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
ADD_LIBRARY(Cumulus STATIC ${Cumulus_SRC})
TARGET_LINK_LIBRARIES(Cumulus MPDD Geometry Utility M2CV ${MPI_Fortran_LIBRARIES} TapenadeAux)

# # 2. Build test executables:
# SET (EXE1 "Tendency_test.exe")
# SET (EXE2 "Test_Cum.exe")
# SET (EXE3 "Test_Cum_AD.exe")
# SET (EXE4 "Test_Cum_interface.exe")
# SET (EXE5 "Test_CumJfun.exe")
# SET (EXE6 "Test_cum_GRAPES.exe")


# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} Tendency_test.F90)
# TARGET_LINK_LIBRARIES(${EXE1} Cumulus MPDD Geometry Utility State State2NC ${MPI_Fortran_LIBRARIES} ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_C}
# ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

#  ADD_EXECUTABLE(${EXE2} Test_Cum.F90)
#  TARGET_LINK_LIBRARIES(${EXE2} Cumulus MPDD Geometry Utility State State2NC ${MPI_Fortran_LIBRARIES} ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_C}
# ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

#  ADD_EXECUTABLE(${EXE3} Test_Cum_AD.F90)
#  TARGET_LINK_LIBRARIES(${EXE3} Cumulus MPDD Geometry Utility State State2NC ${MPI_Fortran_LIBRARIES} TapenadeAux ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_C}
# ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

#  ADD_EXECUTABLE(${EXE4} Test_Cum_interface.F90)
#  TARGET_LINK_LIBRARIES(${EXE4} Cumulus MPDD Geometry Utility State State2NC ${MPI_Fortran_LIBRARIES} TapenadeAux ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_C}
# ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

#  ADD_EXECUTABLE(${EXE5} Test_CumJfun.F90)
#  TARGET_LINK_LIBRARIES(${EXE5} Minimization Cumulus MPDD Geometry Utility State State2NC ${MPI_Fortran_LIBRARIES} 
#  TapenadeAux BMatrix C2O JFunc MGOpts YAMLRead Field State ObsSet Mock Obs2State OprRadar RTTOV ${RTTOV_LIBRARY} ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_C}
# ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

#  ADD_EXECUTABLE(${EXE6} Test_cum_GRAPES.F90)
#  TARGET_LINK_LIBRARIES(${EXE6} Cumulus MPDD Geometry Utility State State2NC ${MPI_Fortran_LIBRARIES} TapenadeAux YAMLRead IOGrapes UV2W ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_C}
# ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

#  # # 4. Add tests:
#  ADD_TEST(NAME Test_Cumulus
#          COMMAND Tendency_test.exe Test_Cum.exe Test_Cum_AD.exe Test_Cum_interface.exe Test_CumJfun.exe
#          WORKING_DIRECTORY $ENV{BIN_DIR})

# # # 5. Set test criteria:
#  SET_TESTS_PROPERTIES(Test_Cumulus PROPERTIES FAIL_REGULAR_EXPRESSION "tested failed")
#  SET_TESTS_PROPERTIES(Test_Cumulus PROPERTIES PASS_REGULAR_EXPRESSION "tested passed")


