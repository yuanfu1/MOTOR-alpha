cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(UV2VortDivg)


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
ADD_LIBRARY(UV2VortDivg STATIC ${UV2VortDivg_SRC})
SET (UV2VortDivg_SRC	
	UV2VortDivg.F90 
	)
TARGET_LINK_LIBRARIES(UV2VortDivg PUBLIC State Field Geometry Utility TransBase lapack blas YAMLRead)

# # 2. Build test executables:
# SET (EXE1 "Test_QRSolver.exe")
# SET (EXE2 "Test_UV2WTLAD.exe")
# SET (EXE3 "Test_UV2WForward.exe")
# SET (EXE4 "Test_UV2WParaTest.exe")

# # 3. Add codes to the executables:
# # ADD_EXECUTABLE(${EXE1} Test_QRSolver.F90)
# # TARGET_LINK_LIBRARIES(${EXE1} UV2W lapack blas Minimization Geometry Utility JFunc State MPDD)

# # ADD_EXECUTABLE(${EXE2} Test_UV2WTLAD.F90)
# # TARGET_LINK_LIBRARIES(${EXE2} UV2W Geometry Utility State MPDD State2NC IOGrapes ${MPI_Fortran_LIBRARIES} ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_C} /public/software/mathlib/netcdf/intel/4.4.0/lib/libnetcdff.so /public/software/mathlib/netcdf/intel/4.4.0/lib/libnetcdf.so 
# # ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY} MGGen MPObs)

# ADD_EXECUTABLE(${EXE3} Test_UV2WForward.F90)
# TARGET_LINK_LIBRARIES(${EXE3} UV2W Geometry Utility State MPDD State2NC IOGrapes ${MPI_Fortran_LIBRARIES} MGGen MPObs)

# # ADD_EXECUTABLE(${EXE4} Test_UV2WParaTest.F90)
# # TARGET_LINK_LIBRARIES(${EXE4} UV2W Geometry Utility State MPDD State2NC IOGrapes ${MPI_Fortran_LIBRARIES} ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_C} /public/software/mathlib/netcdf/intel/4.4.0/lib/libnetcdff.so /public/software/mathlib/netcdf/intel/4.4.0/lib/libnetcdf.so 
# # ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY} MGGen MPObs)
