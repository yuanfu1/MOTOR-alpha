cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(InterpHP)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES(${SLINT_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

# Add Slint libraries.
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# 1. Generate lib file:
ADD_LIBRARY(InterpHP STATIC ${InterpHP_SRC})
TARGET_LINK_LIBRARIES(InterpHP MPDD Geometry Utility ${SLINT_LIBRARY} ${MPI_Fortran_LIBRARIES})

# 2. Build test executables:
#SET (EXE1 "Test_InterpHP.exe")

# 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} Test_IntepHP.F90)
# TARGET_LINK_LIBRARIES(${EXE1} Utility GrapesIO MPDD Geometry IOGrapes State State2NC ${MPI_Fortran_LIBRARIES} ${SLINT_LIBRARY} ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_C}
# ${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

# # 4. Add tests:
# ADD_TEST(NAME Test_InterpHP
#         COMMAND mpirun -np 4 Test_InterpHP.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})

# # 5. Set test criteria:
# SET_TESTS_PROPERTIES(Test_InterpHP PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")