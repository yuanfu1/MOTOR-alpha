cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(RTTOV)


# Add rttov libraries.
INCLUDE_DIRECTORIES(${RTTOV_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${RTTOV_MOD_DIR})

#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

IF (CMAKE_Fortran_COMPILER MATCHES "ifort")
        # # Add hdf5 libraries.
        INCLUDE_DIRECTORIES(${HDF5_Fortran_INCLUDE_DIR})
        LINK_DIRECTORIES(/public/software/mathlib/hdf5/intel/1.8.12/lib)
ENDIF()

# SET (CMAKE_Fortran_FLAGS "-DUSE_MPI -fpic -O3 -cpp -ffree-line-length-none -D_REAL8_ -fopenmp")
# Add Slint libraries.
# INCLUDE_DIRECTORIES(${SLINT_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# 1. Generate lib file:
ADD_LIBRARY(RTTOV STATIC ${RTTOV_SRC})
TARGET_LINK_LIBRARIES(RTTOV PUBLIC MPDD Geometry Utility ObsSet State YAMLRead M2CV ${MPI_Fortran_LIBRARIES} ${ZLIB_LIBRARY}
${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

# # 2. Build test executables:
# SET (EXE1 "NLfeatures_RTTOV_T.exe")
# SET (EXE2 "NLfeatures_RTTOV_Q.exe")
# SET (EXE3 "NLfeatures_RTTOV_TQ.exe")
# SET (EXE4 "NLfeatures_RTTOV_QC.exe")
# SET (EXE5 "NLfeatures_RTTOV_QI.exe")
# SET (EXE6 "NLfeatures_RTTOV_QCQI.exe")
# SET (EXE7 "UnitTest_RTTOV.exe")
# SET (EXE8 "TLAD_check_RTTOV.exe")
# SET (EXE9 "Diag_Jac_WeightFunc.exe")

# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} NLfeatures_RTTOV_T.F90)
# TARGET_LINK_LIBRARIES(${EXE1} Minimization BMatrix C2O JFunc IOWRF MGOpts Geometry M2CV ${MPI_Fortran_LIBRARIES}
# MPDD YAMLRead Field State ObsSet ObsSatellite Mock State2NC Utility RTTOV Obs2State OprRadar lapack
# ${RTTOV_LIBRARY} )

# ADD_EXECUTABLE(${EXE2} NLfeatures_RTTOV_Q.F90)
# TARGET_LINK_LIBRARIES(${EXE2} Minimization BMatrix C2O JFunc IOWRF MGOpts Geometry M2CV ${MPI_Fortran_LIBRARIES}
# MPDD YAMLRead Field State ObsSet ObsSatellite Mock State2NC Utility RTTOV Obs2State OprRadar lapack
# ${RTTOV_LIBRARY} )

# ADD_EXECUTABLE(${EXE3} NLfeatures_RTTOV_TQ.F90)
# TARGET_LINK_LIBRARIES(${EXE3} Minimization BMatrix C2O JFunc IOWRF MGOpts Geometry M2CV ${MPI_Fortran_LIBRARIES}
# MPDD YAMLRead Field State ObsSet ObsSatellite Mock State2NC Utility RTTOV Obs2State OprRadar lapack
# ${RTTOV_LIBRARY} )

# ADD_EXECUTABLE(${EXE4} NLfeatures_RTTOV_QC.F90)
# TARGET_LINK_LIBRARIES(${EXE4} Minimization BMatrix C2O JFunc IOWRF MGOpts Geometry M2CV ${MPI_Fortran_LIBRARIES}
# MPDD YAMLRead Field State ObsSet ObsSatellite Mock State2NC Utility RTTOV Obs2State OprRadar lapack
# ${RTTOV_LIBRARY} )

# ADD_EXECUTABLE(${EXE5} NLfeatures_RTTOV_QI.F90)
# TARGET_LINK_LIBRARIES(${EXE5} Minimization BMatrix C2O JFunc IOWRF MGOpts Geometry M2CV  ${MPI_Fortran_LIBRARIES}
# MPDD YAMLRead Field State ObsSet ObsSatellite Mock State2NC Utility RTTOV Obs2State OprRadar lapack
# ${RTTOV_LIBRARY} )

# ADD_EXECUTABLE(${EXE6} NLfeatures_RTTOV_QCQI.F90)
# TARGET_LINK_LIBRARIES(${EXE6} Minimization BMatrix C2O JFunc IOWRF MGOpts Geometry M2CV ${MPI_Fortran_LIBRARIES}
# MPDD YAMLRead Field State ObsSet ObsSatellite Mock State2NC Utility RTTOV Obs2State OprRadar lapack
# ${RTTOV_LIBRARY} )

# ADD_EXECUTABLE(${EXE7} UnitTest_RTTOV.F90)
# TARGET_LINK_LIBRARIES(${EXE7} Minimization BMatrix C2O JFunc IOWRF MGOpts Geometry M2CV ${MPI_Fortran_LIBRARIES}
# MPDD YAMLRead Field State ObsSet ObsSatellite Mock State2NC Utility RTTOV Obs2State OprRadar lapack
# ${RTTOV_LIBRARY} )

# ADD_EXECUTABLE(${EXE8} TLAD_check_RTTOV.F90)
# TARGET_LINK_LIBRARIES(${EXE8} Minimization BMatrix C2O JFunc IOWRF MGOpts Geometry ${MPI_Fortran_LIBRARIES} 
# MPDD YAMLRead Field State ObsSet ObsSatellite Mock State2NC Utility RTTOV Obs2State OprRadar  lapack
# ${RTTOV_LIBRARY} ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

# ADD_EXECUTABLE(${EXE9} Diag_Jac_WeightFunc.F90)
# TARGET_LINK_LIBRARIES(${EXE9} Minimization BMatrix C2O JFunc IOWRF MGOpts Geometry ${MPI_Fortran_LIBRARIES} 
# MPDD YAMLRead Field State ObsSet ObsSatellite Mock State2NC Utility RTTOV Obs2State OprRadar  lapack
# ${RTTOV_LIBRARY} ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

# # 4. Add tests:
# ADD_TEST(NAME "UnitTest_RTTOV" 
#         COMMAND mpirun -np 1 --oversubscribe UnitTest_RTTOV.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})
# ADD_TEST(NAME "TLAD_check_RTTOV" 
#         COMMAND mpirun -np 1 --oversubscribe TLAD_check_RTTOV.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})
# # 5. Set test criteria:
# SET_TESTS_PROPERTIES(UnitTest_RTTOV PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
# SET_TESTS_PROPERTIES(UnitTest_RTTOV PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
# SET_TESTS_PROPERTIES(TLAD_check_RTTOV PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed!")
# SET_TESTS_PROPERTIES(TLAD_check_RTTOV PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
